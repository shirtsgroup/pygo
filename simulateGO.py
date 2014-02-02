import datetime
import numpy
from optparse import OptionParser
from sys import stdout, exit
import profile
import scipy.linalg
from scipy.misc import comb
import os
import pp
import pdb
import cPickle

numpy.random.seed(10)

import writetopdb
import moveset
import energyfunc
import replicaexchange
import simulationobject
import surfacesimulation
import umbrellasimulation

from simulationobject import Simulation
from surfacesimulation import SurfaceSimulation
from umbrellasimulation import UmbrellaSimulation

t1=datetime.datetime.now()
kb = 0.0019872041 #kcal/mol/K

def parse_args():
    parser=OptionParser()
    
    parser.add_option("-f", "--files", dest="filename", default="GO_1PGB.pdb", help="protein .pdb file")
    parser.add_option("-p", "--paramfile", dest="paramfile", default='GO_1PGB.param', help="protein .param file")
    parser.add_option("-t", "--temprange", nargs=2, default=[200.,400.], type="float", dest="temprange", help="temperature range of replicas")
    parser.add_option("--tfile", dest="tfile", default="", help="file of temperatures")
    parser.add_option("-r", "--nreplicas", default=8, type="int",dest="nreplicas", help="number of replicas")
    parser.add_option("-n", "--moves", dest="totmoves", type="int", default='10000', help="total number of moves")
    parser.add_option("-s", "--save", dest="save", type="int", default='1000', help="number of moves between save operations")
    parser.add_option("-k", "--swapstep", dest="swap", type="int", default='1000', help="number of moves between swap moves")
    parser.add_option("-w", "--writetraj", dest="writetraj", action="store_true", default=False, help="flag to write out trajectory (default: False)")
    parser.add_option("--swapsper", nargs=1, dest="swapsper", type="int",default=500, help="number of swaps to perform (default: 500)")
    parser.add_option("--id", nargs=1, dest="id", type="int", default=0, help="the simlog id number or umbrella id number")
    parser.add_option("--freq", nargs=7, dest="freq", type="float", default=[0,0,1,3,3,3,10], help="ratio of move frequencies (tr:ro:an:di:gc:pr:md)")
    parser.add_option("--md", nargs=2, default=[45,50],type="float", dest="md", help="step size (fs) and nsteps for MD move")
 
    parser.add_option("--surf", action="store_true", default=False, help="surface simulation flag")
    parser.add_option("--surfparamfile", dest="surfparamfile", default='avgsurfparam.npy', help="surface param file")
    parser.add_option("--scale", type="float", default=1, help="scaling of surface attraction strength")

    parser.add_option("--umbrella", type="float", default=0., help="umbrella simulation flag and distance of pinning")
    parser.add_option("--umbrellak", type="float", default=1, help="umbrella spring constant (default: 1)")

    parser.add_option("--Qumbrella", default='', help="Q umbrella simulation flag and file of Q_pins")
    parser.add_option("--Qumbrellak", type="float", default=10, help="Q umbrella spring constant (default: 10)")

    parser.add_option("--cluster", action="store_true", default=False, help="flag for running on cluster")
    parser.add_option("--restart", action="store_true", default=False, help="restart from a checkpoint")
    parser.add_option("--extend", nargs=1, dest="extend", type="int", default=0, help="id number of existing simulation to start simulation from")
    (args,_)=parser.parse_args()
    return args

def get_temperature(args):
    '''Returns array of temperatures'''
    if args.tfile:
        T = numpy.loadtxt(args.tfile)
        assert(len(T)==args.nreplicas)
    elif args.nreplicas == 1:
        if args.temprange[0] != args.temprange[1]:
            print 'WARNING: using lower temperature bound for one replica simulation'
        T = [args.temprange[0]]
    elif args.nreplicas == 2:
        T = args.temprange
    else:
        # Calculate temperature distribution
        T = numpy.empty(args.nreplicas)
        alpha = (args.temprange[1] / args.temprange[0])**(1 / float(args.nreplicas - 1))
        T[0] = args.temprange[0]
        for i in range(1, args.nreplicas):
            T[i] = T[i-1] * alpha
    return T

def get_movefreq(args):
    '''
    Converts move frequency from ratio to cumulative probability
    Ex: [0 0 1 4 3 1 1] => [0 0 .1 .5 .8 .9]
    '''
    if args.surf:
        nmoves = 7 # 7 possible moves for a surface simulation
        start = 0
    else:
        nmoves = 5
        start = 2
        if args.freq[0] or args.freq[1]:
            print 'WARNING: not performing any translation or rotation moves'
            args.freq[0] = 0
            args.freq[1] = 0
    percentmove = numpy.zeros(nmoves-1) 
    tot = float(numpy.sum(args.freq))
    percentmove[0] = args.freq[start]/tot
    for i in range(1,nmoves-1):
        percentmove[i] = percentmove[i-1]+args.freq[start+i]/tot
    return percentmove

def set_up_dir(args):
    # to do: change this before putting on github
    parentdirc = './replicaexchange'
    if not os.path.exists(parentdirc):
        os.mkdir(parentdirc)
    if args.umbrella:
        direc = './replicaexchange/umbrella%i' % args.id
        if not os.path.exists(direc):
            os.mkdir(direc)
        direc = './replicaexchange/umbrella%i/%i' %(args.id,int(args.umbrella))
        if not os.path.exists(direc):
            os.mkdir(direc)
    else:
        direc = './replicaexchange/simlog%i' % args.id
        if not os.path.exists(direc):
            os.mkdir(direc)
    return direc

def pprun(job_server, replicas, moves, dict):
    '''parallel python run'''
    if len(replicas) == 1:
        newreplicas = [simulationobject.run(replicas[0], moves, dict)]
    else:
        jobs = [job_server.submit(simulationobject.run, (replica, moves, dict), 
                (simulationobject.update_energy, simulationobject.save), 
                ("numpy","energyfunc", "moveset", "writetopdb","cPickle")) for replica in replicas]
        newreplicas = [job() for job in jobs]
    return newreplicas

def pprun_surf(job_server, replicas, moves, dict):
    '''parallel python run with surface'''
    if len(replicas) == 1:
        newreplicas = [surfacesimulation.run_surf(replicas[0], moves, dict)]
    else:
        jobs = [job_server.submit(surfacesimulation.run_surf, (replica, moves, dict), 
                (surfacesimulation.update_energy, surfacesimulation.save), 
                ("numpy","energyfunc", "moveset", "writetopdb","cPickle")) for replica in replicas]
        newreplicas = [job() for job in jobs]
    return newreplicas

def pprun_Q(job_server, replicas, moves, dict):
    jobs = [job_server.submit(surfacesimulation.run_surf, (replica, moves, dict), 
                (update_energy_Q, surfacesimulation.save), 
                ("numpy","energyfunc", "moveset", "writetopdb","cPickle")) for replica in replicas]
    newreplicas = [job() for job in jobs]
    return newreplicas

def pprun_Q_surf(job_server, replicas, moves, dict):
    jobs = [job_server.submit(surfacesimulation.run_surf, (replica, moves, dict), 
                (surfacesimulation.update_energy, surfacesimulation.save), 
                ("numpy","energyfunc", "moveset", "writetopdb","cPickle")) for replica in replicas]
    newreplicas = [job() for job in jobs]
    return newreplicas

def savestate(args, replicas, protein_location):
    output = open('%s/cptstate.pkl' % args.direc, 'wb')
    cPickle.dump(replicas[0].move, output)
    cPickle.dump(protein_location, output)
    output.close()
    for i in range(args.nreplicas):
        replicas[i].saveenergy(False)
        replicas[i].savenc(False)
        replicas[i].savecoord()
        if args.surf:
            replicas[i].savesurfenergy(False)
        if umbrella:
            replicas[i].save_z()

def loadstate(args, replicas, protein_location):
    input = open('%s/cptstate.pkl' % args.direc, 'rb')
    move = cPickle.load(input)
    protein_location = cPickle.load(input)
    input.close()
    for i in range(args.nreplicas):
       	replicas[i].loadstate()
        replicas[i].move = move
        replicas[protein_location[i][-1]].whoami = i 
    return move, replicas, protein_location
 
def main():
    print 'Start time', t1
    print ''

    # --- process inputs --- #
    args = parse_args()
    coord,_ = writetopdb.get_coord(args.filename)
    numbeads=len(coord)
    mass = energyfunc.getmass('%stop' % (args.paramfile[0:-5]), numbeads)
    T = get_temperature(args)
    beta = 1/(kb*T)
    percentmove = get_movefreq(args)
    tsize, tsteps = args.md
    tsize /= 100. # input is in fs
    tsteps = int(tsteps)
    direc = set_up_dir(args)

    # --- report inputs --- #
    print ''
    print '-----Inputs-----'
    print 'Output directory:', os.path.abspath(direc)
    print 'System:', args.filename
    if args.umbrella:
        print '    Running surface umbrella sampling simulation with pinning at %f A' % args.umbrella
        assert(args.surf == True)
    if args.surf:
        print '    Running surface simulation'
    print 'Number of temperature replicas:', args.nreplicas
    print 'Temperature(s):'
    print T
    
    print 'Total number of moves:', args.totmoves
    print 'Save interval:', args.save
    print 'Replica exchange interval:', args.swap
    print 'There are %i swaps at each exchange point' %(args.swapsper)
    print 'Ratio of moves frequencies (tr:rot:ang:dih:crank:parrot:MD) is', args.freq
    print 'MD steps per move is', args.md[1]
    print 'MD time step is', args.md[0],' fs'
    print ''

    # --- get parameters from .param file --- #
    angleparam = energyfunc.getangleparam(args.paramfile, numbeads)
    torsparam = energyfunc.gettorsionparam(args.paramfile, numbeads)
    
    # --- pregenerate list of interactions --- #
    numint = numpy.around(comb(numbeads, 2)) # number of interactions
    numint = numint - 2 * (numbeads - 2) - 1 # don't count 12 and 13 neighbors
    
    # --- get native LJ parameter --- #
    nativeparam = energyfunc.getnativefix_n(args.paramfile, numint, numbeads) # [ones and zeros, native epsilon, native sigma]
    totnc = numpy.sum(nativeparam[:,0]) #total native contacts
    nativecutoff2 = 1.2**2
    nsigma2 = nativecutoff2 * nativeparam[:,2] * nativeparam[:,2]
    
    # --- get nonnative LJ parameter --- #
    (nonnativesig, nnepsil) = energyfunc.getLJparam_n(args.paramfile, numbeads, numint) #[nonnative sigmas for every interaction, epsilon (one value)]
    nonnatindex = -1 * (nativeparam[:,0] - 1) # array of ones and zeros
    nonnativeparam = numpy.column_stack((nonnatindex, nonnativesig)) #[ones and zeros, nonnative sigma]
    
    # --- set Simulation class variables --- #
    Simulation.angleparam = angleparam
    Simulation.torsparam = torsparam
    Simulation.totnc = totnc
    Simulation.nativeparam = nativeparam
    Simulation.nsigma2 = nsigma2
    Simulation.nonnativeparam = nonnativeparam
    Simulation.nnepsil = nnepsil
    Simulation.totmoves = args.totmoves
    Simulation.save = args.save
    Simulation.numbeads = numbeads
    Simulation.numint = numint
    Simulation.mass = mass
    Simulation.percentmove = percentmove
    Simulation.tsize = tsize
    Simulation.tsteps = tsteps
    
    # --- put class variables in a dictionary for pprun --- #
    dict = {'tsize':tsize,
            'tsteps':tsteps,
            'percentmove':percentmove,
            'numbeads':numbeads,
            'save':args.save, 
            'totmoves':args.totmoves, 
            'numint':numint, 
            'angleparam':angleparam, 
            'torsparam':torsparam, 
            'nativeparam':nativeparam, 
            'nonnativeparam':nonnativeparam, 
            'nnepsil':nnepsil, 
            'nsigma2':nsigma2, 
            'writetraj':args.writetraj, 
            'mass':mass}

    # --- set up surface --- #
    if args.surf:
        # to do: not have surface be hardcoded
        xlength = 135
        ylength = 135
        spacing = 7
        yspacing = spacing*3.**.5
        surface = surfacesimulation.getsurf(xlength+15,ylength+15,spacing)
        surfacesimulation.writesurf('surface.pdb',surface)
        nsurf = len(surface)
        nspint = nsurf*numbeads # surface-protein interactions
        sfile = args.paramfile[0:args.paramfile.find('GO_')]+args.paramfile[args.paramfile.find('GO_')+3:-5] + 'pdb'
        surfparam = energyfunc.getsurfparam(sfile, args.surfparamfile, numbeads, nsurf, nspint, args.scale)
        SurfaceSimulation.scale = args.scale
        SurfaceSimulation.surface = surface
        SurfaceSimulation.nsurf = nsurf
        SurfaceSimulation.nspint = nspint
        SurfaceSimulation.surfparam = surfparam
        dict.update({'nspint':nspint, 
                    'nsurf':nsurf, 
                    'surfparam':surfparam, 
                    'surface':surface, 
                    'xlength':xlength, 
                    'ylength':ylength, 
                    'spacing':spacing, 
                    'yspacing':yspacing, 
                    'scale':args.scale})
        print '-----Surface Details-----'
        print 'Surface is %i by %i array of leucine residues with spacing %i A' % (xlength, ylength, spacing)
        print 'Surface energy parameters scaled by %f' % args.scale
        print ''


    # --- instantiate replicas --- #
    replicas = []
    for i in range(args.nreplicas):
        name = 'replica%i' % i
        if args.umbrella:
            replicas.append(UmbrellaSimulation(name, os.path.abspath(direc), coord, T[i], surface, args.umbrella, mass, args.umbrellak))
        elif args.surf:
            replicas.append(SurfaceSimulation(name, os.path.abspath(direc), coord, T[i], surface))
        else:
            replicas.append(Simulation(name, os.path.abspath(direc), coord, T[i]))
        replicas[i].whoami = i
        if args.writetraj:
            f = open('%s/trajectory%i' %(replicas[i].out, int(replicas[i].T)), 'wb')
            numpy.save(f,replicas[i].coord)
            f.close()

    # --- begin simulation --- #
    if args.umbrella:
    	print 'Starting umbrella simulation... at %f A' % args.umbrella
        PPRUN = pprun_surf
    elif args.surf:
    	print 'Starting surface simulation...'
        PPRUN = pprun_surf
    else:
    	print 'Starting simulation...'
        PPRUN = pprun
    move = 0
    swapaccepted = numpy.zeros(args.nreplicas-1)
    swaprejected = numpy.zeros(args.nreplicas-1)
    protein_location = [[i] for i in range(args.nreplicas)]

    # --- setup ppserver --- #
    if args.cluster:
    	print '    Running on a cluster...'
    	if args.umbrella:
    		try:
    			f = open('nodefile%i-%s.txt'% (id,umbrella),'r')
    		except:
    			f = open('nodefile%i-%i.txt'% (id,umbrella),'r')
    	else:
    		f = open('nodefile%i.txt'% id,'r')
    	ppservers = f.read().split("\n")
    	f.close()
    	ppservers = filter(None,ppservers)
    	ppservers = [x+':43334' for x in ppservers]
    	ppservers = tuple(ppservers)
    	job_server = pp.Server(0,ppservers=ppservers)
    	print 'Running pp on: '
    	print ppservers
    	print "Starting pp with", job_server.get_ncpus(), "workers"
    else:
    	# running on one machine
    	job_server = pp.Server(ppservers=())
    	print "Starting pp with", job_server.get_ncpus(), "workers"
 
    # --- set up simulation extension --- #
    if args.extend:
    	print '    Extending simulation %i...' % args.extend
    	if args.umbrella:
    		extenddirec = os.getcwd()+'/replicaexchange/umbrella%i/%i' % (args.extend,int(args.umbrella))
    	else:
    		extenddirec = os.getcwd()+'/replicaexchange/simlog%i' % args.extend
    	if not os.path.exists(extenddirec):
    		sys.exit('Simulation %i does not exist at %s' % (args.extend, extenddirec))
    	input = open('%s/protein_location.pkl' %extenddirec, 'rb')
    	protein_location = cPickle.load(input)
    	protein_location = [[protein_location[i][-1]] for i in range(args.nreplicas)]	
    	input.close()
    	for i in range(args.nreplicas):
    		replicas[i].loadextend(extenddirec)
    		replicas[protein_location[i][-1]].whoami = i
    
    # --- load existing simulation --- #
    if args.restart:
    	print '    Restarting from last checkpoint...'
    	move = loadstate()

    ti = datetime.datetime.now()
    tcheck = ti
    for i in xrange(move/args.swap,args.totmoves/args.swap):
        replicas = PPRUN(job_server,replicas, args.swap, dict)
        job_server.wait()
        if args.swap!=args.totmoves:
        	swapaccepted, swaprejected, protein_location = replicaexchange.tryrepeatedswaps(args, replicas, swapaccepted, swaprejected, protein_location, beta)
        tnow = datetime.datetime.now()
        t_remain = (tnow - ti)/(i+1)*(args.totmoves/args.swap - i - 1)
        if not args.cluster:
    	    stdout.write(str(t_remain) + '\r')
    	    stdout.flush()
        # checkpoint
        if tnow-tcheck > datetime.timedelta(seconds=900): #every 15 minutes
            savestate()
            f = open('%s/status.txt' % direc, 'w')
            f.write('Completed %i moves out of %i moves\n' %(replicas[0].move, args.totmoves))
            f.write('%i swaps performed in %s\n' %(i, str(tnow-ti))) #useful for the MD vs. MC comparison
            f.close()
            tcheck=tnow

    # --- output --- #
    job_server.print_stats()
    output = open('%s/protein_location.pkl' % direc, 'wb')
    cPickle.dump(protein_location, output)
    output.close()
    for i in range(args.nreplicas):
        replicas[i].output()
        replicas[i].saveenergy()
        replicas[i].savenc()
        print 'The average Q is %f' %(numpy.average(replicas[i].nc)/Simulation.totnc)
        replicas[i].savecoord()
        if args.surf:
                replicas[i].savesurfenergy()
        if args.umbrella:
    	    replicas[i].save_z()
    
    if args.swap!=args.totmoves:
        Q_trajec_singleprot = numpy.zeros((args.nreplicas, args.totmoves/args.save+1))
        k=args.swap/args.save
        for i in xrange(len(protein_location[0])):
                for j in range(args.nreplicas):
                    rep = protein_location[j][i]
                    Q_trajec_singleprot[j,k*i+1:k*(i+1)+1] = replicas[rep].nc[k*i+1:k*(i+1)+1]
        Q_trajec_singleprot[:,0] = totnc
        Q_trajec_singleprot = Q_trajec_singleprot/totnc
        numpy.save('%s/Qtraj_singleprot.npy' % direc, Q_trajec_singleprot)
    for i in range(args.nreplicas-1):
        print 'swaps accepted between replica%i and replica%i: %3.2f percent' % (i, i+1, (swapaccepted[i] / float(swapaccepted[i] + swaprejected[i]) * 100))
    print 'Total swaps accepted: %i' % numpy.sum(swapaccepted)
    print 'Total swaps rejected: %i' % numpy.sum(swaprejected)
    
if __name__ == '__main__':
    main()
    t2 = datetime.datetime.now()
    print 'Simulation time: '+str(t2-t1)
    
    
