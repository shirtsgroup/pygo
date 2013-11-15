#=======================================================================================================
# IMPORTS
#=======================================================================================================

import datetime
import numpy
from optparse import OptionParser
try:
	# does not work on the cluster
	import matplotlib.pyplot as plt
except:
	plot=False
	import os
	os.environ['TMPDIR']
from sys import stdout, exit
import profile
import scipy.linalg
from scipy.misc import comb
import os
import pp
import pdb
import cPickle

from writetopdb import *
from moveset import *
from energyfunc import *
#from simulationobject import *

t1=datetime.datetime.now()

parser=OptionParser()
parser.add_option("-f", "--files", dest="datafile", default="GO_1PGB.pdb", help="protein .pdb file")
parser.add_option("-p", "--parameterfile", dest="paramfile", default='GO_1PGB.param', help="protein .param file")
parser.add_option("-t", "--temprange", nargs=2, default=[300.,450.], type="float", dest="temprange", help="temperature range of replicas")
parser.add_option("-r", "--replicas", default=8, type="int",dest="replicas", help="number of replicas")
parser.add_option("-v", "--verbose", action="store_false", default=True, help="more verbosity")
parser.add_option("-c", "--addconnect", action="store_true", default=False, help="add bonds to .pdb file")
parser.add_option("-n", "--moves", dest="totmoves", type="int", default='10000', help="total number of moves")
parser.add_option("-s", "--stepsize", dest="step", type="int", default='1000', help="number of moves between save operations")
parser.add_option("-k", "--swapstep", dest="swap", type="int", default='1000', help="number of moves between swap moves")
parser.add_option("-g", "--histogram", dest="histname", default='', help="name histogram of conformational energies, if desired")
parser.add_option("-a", "--plot", action="store_true", default=False, help="plots energy, rmsd, fractional nativeness")
parser.add_option("-b", "--writepdb", action="store_true", default=False, help="the output pdb file")
parser.add_option("-e", "--swaps", nargs=1, dest="e", type="int",default=500, help="number of swaps to perform")
parser.add_option("--id", nargs=1, dest="id", type="int", default=0, help="the simlog id number or umbrella id number")
parser.add_option("--freq", nargs=6, dest="freq", type="float", default=[.2,.4,.7,1.,1.,1.], help="move frequencies")
parser.add_option("--surf", action="store_true", default=False, help="surface simulation flag")
parser.add_option("--umbrella", type="float", default=0., help="umbrella simulation flag and distance of pinning")
parser.add_option("--scale", type="float", default=1, help="scaling of surface attraction strength")
parser.add_option("--cluster", action="store_true", default=False, help="flag for running on cluster")
parser.add_option("--restart", action="store_true", default=False, help="restart from a checkpoint")
parser.add_option("--md", nargs=2, default=[.05,360],type="float", dest="md", help="step size and nsteps for MD move")
parser.add_option("--tfile", dest="tfile", default="", help="file of temperatures")
parser.add_option("--extend", nargs=1, dest="extend", type="int", default=0, help="id number of existing simulation to start simulation from")
#parser.add_option("--random", type="int", dest="random", default=10, help="random seed for reproducability")

(options,args)=parser.parse_args()

#======================================================================================================
# CONSTANTS
#======================================================================================================
verbose = options.verbose
scale = options.scale
surf = options.surf
umbrella = options.umbrella
if surf:
    from surfacesimulation import *
elif umbrella:
    from umbrellasimulation import *
    surf = True
else:
    from simulationobject import *
trange = options.temprange # Kelvin
tfile = options.tfile
numreplicas = options.replicas
totmoves = options.totmoves
step = options.step
swap = options.swap
filename = options.datafile
paramfile = options.paramfile
plot = options.plot
histname = options.histname
pdbfile = options.writepdb
e=options.e
addbonds = options.addconnect
id = options.id # simlog id number, used to make directory for results
extend = options.extend #simlog id number to extend simulation of
percentmove = options.freq
cluster = options.cluster
restart = options.restart
[tsize,tsteps] = options.md
tsteps = int(tsteps)
if restart:
	print 'Restarting from last checkpoint'
if extend:
	print 'Extending simulation %i' % extend
#if options.random:
#	random.seed(options.random)
#random.seed(10)
kb = 0.0019872041 #kcal/mol/K
nativecutoff2 = 1.2**2

# read .pdb to get number of beads (numbeads)
file = open(filename, 'r')
numbeads = 0
while 1:
	line = file.readline()
	if not line:
		break
	splitline = line.split('  ')
	if splitline[0] == 'ATOM':
		numbeads += 1

if (addbonds):
	addconnect(pdbfile, numbeads)
	
# gets bead coordinates from .pdb file
file.seek(0)
coord = numpy.empty((numbeads, 3))
i = 0
k = 0
wordtemplate = [] #template of all words
positiontemplate = [] #template of position lines
ATOMlinenum = [] # array of location of ATOM lines
while 1:
    line = file.readline()
    if not line:
        break
    wordtemplate.append(line)
    splitline = line.split('  ')
    if splitline[0] == 'ATOM':
        positiontemplate.append(line)
        ATOMlinenum.append(k)
        coord[i,0] = float(line[31:38])
        coord[i,1] = float(line[39:46])
        coord[i,2] = float(line[47:54])
        i += 1
    k += 1
file.close()

untransform = getmovietransform(coord)
transform = transpose(untransform)
#coord_nat = moviecoord(coord,transform)

#Get parameters from .param file
angleparam = getangleparam(paramfile, numbeads)
torsparam = gettorsionparam(paramfile, numbeads)
mass = getmass('%stop' % (paramfile[0:-5]), numbeads)

#pregenerate list of interactions
numint = around(comb(numbeads, 2)) # number of interactions
numint = numint - 2 * (numbeads - 2) - 1 # don't count 12 and 13 neighbors

#native LJ parameter getting
nativeparam_n = getnativefix_n(paramfile, numint, numbeads) # [ones and zeros, native epsilon, native sigma]

totnc = sum(nativeparam_n[:,0]) #total native contacts
nsigma2 = nativecutoff2 * nativeparam_n[:,2] * nativeparam_n[:,2]

#nonnative LJ parameter getting
[nonnativesig, nnepsil] = getLJparam_n(paramfile, numbeads, numint) #[nonnative sigmas for every interaction, epsilon (one value)]
nonnatindex = -1 * (nativeparam_n[:,0] - 1) # array of ones and zeros
nonnativeparam = column_stack((nonnatindex, nonnativesig)) #[ones and zeros, nonnative sigma]

# set Simulation class variables
Simulation.angleparam = angleparam
Simulation.torsparam = torsparam
Simulation.totnc = totnc
Simulation.nativeparam_n = nativeparam_n
Simulation.nsigma2 = nsigma2
Simulation.nonnativeparam = nonnativeparam
Simulation.nnepsil = nnepsil
Simulation.totmoves = totmoves
Simulation.step = step
Simulation.numbeads = numbeads
Simulation.numint = numint
Simulation.wordtemplate = wordtemplate
Simulation.positiontemplate = positiontemplate
Simulation.ATOMlinenum = ATOMlinenum
Simulation.transform = transform
#Simulation.coord_nat = coord_nat
Simulation.mass = mass
Simulation.percentmove = percentmove
Simulation.tsize = tsize
Simulation.tsteps = tsteps

# put class variables in a dictionary for pprun
dict = {'tsize':tsize,'tsteps':tsteps,'percentmove':percentmove,'numbeads':numbeads, 'step':step, 'totmoves':totmoves, 'numint':numint, 'angleparam':angleparam, 'torsparam':torsparam, 'nativeparam_n':nativeparam_n, 'nonnativeparam':nonnativeparam, 'nnepsil':nnepsil, 'nsigma2':nsigma2, 'transform':transform, 'positiontemplate':positiontemplate, 'pdbfile':pdbfile, 'mass':mass}

# Calculate temperature distribution
if numreplicas == 1:
    if trange[0] != trange[1]:
        print 'WARNING: using lower temperature bound for one replica simulation'
    T = [trange[0]]
elif numreplicas == 2:
    T = trange
elif not tfile:
    T = numpy.empty(numreplicas)
    alpha = (trange[1] / trange[0])**(1 / float(numreplicas - 1))
    T[0] = trange[0]
    for i in range(1, numreplicas):
        T[i] = T[i-1] * alpha
else:
    T = numpy.loadtxt(tfile)
    assert(len(T)==numreplicas)

if verbose:
    if cluster:
	print 'Running on a cluster'
    output = ['Verbosity is %s' %(str(verbose)),
    'The total number of moves is %d' %(totmoves),
    'The save interval is %d moves' %(step),
    'The replica exchange interval is %d moves' %(swap),
    'There are %i swaps at each exchange point' %(e),
    'There are %d replicas at temperatures:' %(numreplicas),
    str(T),
    'There are %d residues in %s' %(numbeads,filename),
    'Percent move is',
    str(Simulation.percentmove),'']
    print "\r\n".join(output)
    print 'MD steps per move is '+str(Simulation.tsteps)
    print 'MD time step is '+str(Simulation.tsize*100)+' fs'

if surf:
    xlength = 135
    ylength = 135
    spacing = 7
    print 'Surface is %i by %i with spacing %i (Angstroms)' %( xlength, ylength, spacing)
    yspacing = spacing*3.**.5
    surface = getsurf(xlength+15,ylength+15,spacing)
    writesurf('surface.pdb',surface)
    nsurf = len(surface)
    nspint = nsurf*numbeads # surface-protein interactions
    sfile = paramfile[0:paramfile.find('GO_')]+paramfile[paramfile.find('GO_')+3:-5] + 'pdb'
    surfparam = getsurfparam(sfile, numbeads, nsurf, nspint, scale)
    #/surfparam[:,0] = surfparam[:,0]*scale
    print 'Surface energy parameters scaled by %f' % scale
    SurfaceSimulation.scale = scale
    SurfaceSimulation.surface = surface
    SurfaceSimulation.nsurf = nsurf
    SurfaceSimulation.nspint = nspint
    SurfaceSimulation.surfparam = surfparam
    dict.update({'nspint':nspint, 'nsurf':nsurf, 'surfparam':surfparam, 'surface':surface, 'xlength':xlength, 'ylength':ylength, 'spacing':spacing, 'yspacing':yspacing, 'scale':scale})
#======================================================================================================
# SUBROUTINES
#======================================================================================================

#def tryswap(Replicas, Swapaccepted, Swaprejected, reploc, start):
def tryswap(Replicas, Swapaccepted, Swaprejected, start):
    if numreplicas == 1:
        pass
    else:
        for i in xrange(start, numreplicas-1, 2):
            P=numpy.exp((Replicas[i].u0-Replicas[i+1].u0)*(1/(kb*T[i])-1/(kb*T[i+1])))
            if(random.random()<P):
                Swapaccepted[i]+=1
                Replicas[i].coord, Replicas[i+1].coord = Replicas[i+1].coord, Replicas[i].coord
                Replicas[i].u0, Replicas[i+1].u0 = Replicas[i+1].u0, Replicas[i].u0
                Replicas[i].r2, Replicas[i+1].r2 = Replicas[i+1].r2, Replicas[i].r2
                Replicas[i].torsE, Replicas[i+1].torsE = Replicas[i+1].torsE, Replicas[i].torsE
                Replicas[i].angE, Replicas[i+1].angE = Replicas[i+1].angE, Replicas[i].angE
                Replicas[i].whoami, Replicas[i+1].whoami = Replicas[i+1].whoami, Replicas[i].whoami # whoami keeps track of the individual protein in each Replica 
 #               reploc[i], reploc[i+1] = reploc[i+1], reploc[i]
                if surf:
                    Replicas[i].surfE, Replicas[i+1].surfE = Replicas[i+1].surfE, Replicas[i].surfE
		if umbrella:
		    Replicas[i].z_array, Replicas[i+1].z_array = Replicas[i+1].z_array, Replicas[i].z_array

		#Whoiswhere[Replicas[i].whoami].append(i)
                #Whoiswhere[Replicas[i+1].whoami].append(i+1)
            else:
                Swaprejected[i]+=1
                #Whoiswhere[Replicas[i].whoami].append(Whoiswhere[Replicas[i].whoami][-1])
                #Whoiswhere[Replicas[i+1].whoami].append(Whoiswhere[Replicas[i+1].whoami][-1])
    #if numreplicas%2 == 1:
        #Whoiswhere[Replicas[numreplicas-1].whoami].append(Whoiswhere[Replicas[numreplicas-1].whoami][-1])
#    return Swapaccepted, Swaprejected, reploc
    return Swapaccepted, Swaprejected

#def tryrepeatedswaps(Replicas, Swapaccepted, Swaprejected, protein_location, transmat):
def tryrepeatedswaps(Replicas, Swapaccepted, Swaprejected, protein_location):
	if numreplicas == 1:
#		return Swapaccepted, Swaprejected, protein_location, transmat
		return Swapaccepted, Swaprejected, protein_location
#	replica_location = numpy.arange(numreplicas) # keeps track of where the replicas go after 100 (default) swaps in order to properly increment the transition matrix
	fairweather = 0
	for k in range(e): # default: try swapping 100 times
		if fairweather:
#			Swapaccepted, Swaprejected, replica_location = tryswap(Replicas, Swapaccepted, Swaprejected, replica_location, 0) #type 1 switch
			Swapaccepted, Swaprejected = tryswap(Replicas, Swapaccepted, Swaprejected, 0) #type 1 switch
		else:	
#			Swapaccepted, Swaprejected, replica_location = tryswap(Replicas, Swapaccepted, Swaprejected, replica_location, 1) #type 2 switch
			Swapaccepted, Swaprejected = tryswap(Replicas, Swapaccepted, Swaprejected, 1) #type 2 switch
		fairweather = -(fairweather-1)
#	replica_location_indexed = numpy.zeros(numreplicas) # keeps track of which replica went where
	for i in range(numreplicas):
#		replica_location_indexed[replica_location[i]] = i #extracts which replica went where
		protein_location[Replicas[i].whoami].append(i) #extracts which protein went where
#	for i in range(numreplicas):
#		transmat[i,replica_location_indexed[i]] += 1 #increments the transition matrix
#	return Swapaccepted, Swaprejected, protein_location, transmat
	return Swapaccepted, Swaprejected, protein_location

#def update_energy(self, torschange, angchange):
    #self.newtorsE = energyfunc.ctorsionenergy(self.newcoord, self.torsE, Simulation.torsparam, torschange)
    #self.newangE = energyfunc.cangleenergy(self.newcoord, self.angE, Simulation.angleparam, angchange)
    #self.u1 = sum(self.newtorsE)+sum(self.newangE) + energyfunc.cgetLJenergy(self.newcoord, Simulation.numint, Simulation.numbeads, Simulation.nativeparam_n, Simulation.nonnativeparam, Simulation.nnepsil)
    #return self

if not surf:
    def pprun(Replicas, Moves, Dict):
        if len(Replicas) == 1:
            newReplicas = [run(Replicas[0], Moves, Dict)]
        else:
            jobs = [job_server.submit(run, (replica, Moves, Dict), (update_energy, save), ("random", "numpy",
                "energyfunc", "moveset", "writetopdb","cPickle")) for replica in Replicas]
            newReplicas = [job() for job in jobs]
        return newReplicas
else:
    def pprun(Replicas, Moves, Dict):
        if len(Replicas) == 1:
            newReplicas = [run_surf(Replicas[0], Moves, Dict)]
        else:
            jobs = [job_server.submit(run_surf, (replica, Moves, Dict), (update_energy, save), ("random", "numpy",
                "energyfunc", "moveset", "writetopdb","cPickle")) for replica in Replicas]
            newReplicas = [job() for job in jobs]
        return newReplicas
#def energy(mpos, rsquare, torsE, angE):
    #"""Returns the total energy of a configurat"""
    #energy = sum(angE) + sum(torsE) + energyfunc.LJenergy_n(rsquare, Simulation.nativeparam_n,
             #Simulation.nonnativeparam, Simulation.nnepsil)
    #return energy

def savestate():
    output = open('%s/cptstate.pkl' % direc, 'wb')
    cPickle.dump(replicas[0].move, output)
    cPickle.dump(protein_location, output)
#    cPickle.dump(transmat, output)
    output.close()
    for i in range(len(replicas)):
        replicas[i].saveenergy(False)
        replicas[i].savenc(False)
	replicas[i].savecoord()
#       replicas[i].savermsd(plot)
#       replicas[i].savehist(plot)
        if surf:
            replicas[i].savesurfenergy(False)
        if umbrella:
	    replicas[i].save_z()

def loadstate():
    input = open('%s/cptstate.pkl' % direc, 'rb')
    move = cPickle.load(input)
    protein_location = cPickle.load(input)
#    try:
#	transmat = cPickle.load(input)
#    except:
#	pass
    input.close()
    for i in range(numreplicas):
       	replicas[i].loadstate()
	replicas[i].move = move
	replicas[protein_location[i][-1]].whoami = i 
    return move

#=======================================================================================================
# SIMULATE
#=======================================================================================================
move = 0
swapaccepted = numpy.zeros(numreplicas-1)
swaprejected = numpy.zeros(numreplicas-1)
protein_location = [[i] for i in range(numreplicas)]
#transmat = [zeros((numreplicas,numreplicas)) for i in range(10)]
#if swap!=totmoves:
#    assert(totmoves / swap % 10 == 0)
#    transmat_index = totmoves/swap/10

# instantiate replicas
replicas = []
parentdirc = './replicaexchange'
if not os.path.exists(parentdirc):
    os.mkdir(parentdirc)

if umbrella:
    direc = './replicaexchange/umbrella%i' %id
    if not os.path.exists(direc):
        os.mkdir(direc)
    direc = './replicaexchange/umbrella%i/%i' %(id,int(umbrella))
    if not os.path.exists(direc):
        os.mkdir(direc)

else:
    direc = './replicaexchange/simlog%i' % id
    if not os.path.exists(direc):
        os.mkdir(direc)


for i in range(len(T)):
        name = 'replica%i' % i
	if umbrella:
		replicas.append(UmbrellaSimulation(name, os.path.abspath(direc), coord, T[i], surface, umbrella, mass))
	elif surf:
		replicas.append(SurfaceSimulation(name, os.path.abspath(direc), coord, T[i], surface))
	else:
        	replicas.append(Simulation(name, os.path.abspath(direc), coord, T[i]))
        replicas[i].whoami = i
        if pdbfile:
           # mcoord = moviecoord(coord, transform)
           # writepdb(replicas[i].coord, wordtemplate, ATOMlinenum, 0, '%s/trajectory%i.pdb' % (replicas[i].out, int(replicas[i].T)))
	     f = open('%s/trajectory%i' %(replicas[i].out, int(replicas[i].T)), 'wb')
	     numpy.save(f,replicas[i].coord)
	     f.close()

if extend:
	if umbrella:
		extenddirec = os.getcwd()+'/replicaexchange/umbrella%i/%i' % (extend,int(umbrella))
	else:
		extenddirec = os.getcwd()+'/replicaexchange/simlog%i' % extend
	if not os.path.exists(extenddirec):
		sys.exit('Simulation %i does not exist at %s' % (extend, extenddirec))
	input = open('%s/protein_location.pkl' %extenddirec, 'rb')
	protein_location = cPickle.load(input)
	protein_location = [[protein_location[i][-1]] for i in range(numreplicas)]	
	input.close()
	for i in range(numreplicas):
		replicas[i].loadextend(extenddirec)
		replicas[protein_location[i][-1]].whoami = i

if cluster:
	# running on the cluster
	if umbrella:
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
	#import time
	#time.sleep(35)
	print 'Running pp on: '
	print ppservers
	print "Starting pp with", job_server.get_ncpus(), "workers"
else:
	# running on one machine
	job_server = pp.Server(ppservers=())
	print "Starting pp with", job_server.get_ncpus(), "workers"
if restart:
	move = loadstate()

ti = datetime.datetime.now()
tcheck = ti
if umbrella:
	print 'Starting umbrella simulation... at %f' % umbrella
elif surf:
	print 'Starting surface simulation...'
else:
	print 'Starting simulation...'
for i in xrange(move/swap,totmoves/swap):
    replicas = pprun(replicas, swap, dict)
    job_server.wait()
    if swap!=totmoves:
    	#swapaccepted, swaprejected, protein_location, transmat[i/transmat_index] = tryrepeatedswaps(replicas, swapaccepted, swaprejected, protein_location, transmat[i/transmat_index])
    	swapaccepted, swaprejected, protein_location = tryrepeatedswaps(replicas, swapaccepted, swaprejected, protein_location)
    tnow = datetime.datetime.now()
    t_remain = (tnow - ti)/(i+1)*(totmoves/swap - i - 1)
    if not cluster:
	stdout.write(str(t_remain) + '\r')
	stdout.flush()
    # checkpoint
    if tnow-tcheck > datetime.timedelta(seconds=900): #every 15 minutes
	savestate()
	f = open('%s/status.txt' % direc, 'w')
	f.write('Completed %i moves out of %i moves\n' %(replicas[0].move, totmoves))
	f.write('%i swaps performed in %s\n' %(i, str(tnow-ti))) #useful for the MD vs. MC comparison
	f.close()
	tcheck=tnow
#=======================================================================================================
# POST ANALYSIS
#=======================================================================================================
#if swap!=totmoves and numreplicas > 1:
    #calculate mixing time
#    mixtime = zeros(10)
#    fulltransmat = zeros((numreplicas,numreplicas))
#    for i in range(10):
#            fulltransmat += transmat[i]
#            transmat[i] = transmat[i]/transmat_index #normalize to probabilities
#            transmatrix = numpy.matrix(transmat[i])
#            eig = scipy.linalg.eigvals(transmatrix)
#            eig = numpy.sort(eig)
#            print '-----Transition Matrix '+str(i)+' Eigenvalues-----'
#            print eig
#            mixtime[i] = 1/(1-eig[-2])
#            print mixtime[i]
#    fulltransmat = fulltransmat/totmoves*swap
#    try:
#            fulltransmat = numpy.matrix(fulltransmat)
#            eig = scipy.linalg.eigvals(fulltransmat)
#            eig = numpy.sort(eig)
#            print '-----Full Transition Matrix Eigenvalues-----'
#            print eig
#            print 'Full Tmat miximg time: ' + str(1/(1-eig[-2]))
#            print 'Average mixing time from submatrices: ' + str(average(mixtime))
#            print 'Standard deviation in mixing time: ' + str(std(mixtime))
#    except:
#            pdb.set_trace()
#=======================================================================================================
# OUTPUT
#=======================================================================================================
job_server.print_stats()
output = open('%s/protein_location.pkl' % direc, 'wb')
cPickle.dump(protein_location, output)
output.close()
for i in range(len(replicas)):
    replicas[i].output(verbose)
    replicas[i].saveenergy(plot)
    replicas[i].savenc(plot)
    print 'The average Q is %f' %(numpy.average(replicas[i].nc)/Simulation.totnc)
    replicas[i].savecoord()
#    replicas[i].savermsd(plot)
#    replicas[i].savehist(plot)
    if surf:
            replicas[i].savesurfenergy(plot)
    if umbrella:
	    replicas[i].save_z()

if swap!=totmoves:
    Q_trajec_singleprot = numpy.zeros((numreplicas, totmoves/step+1))
    # assumes swap interval is larger than or equal to step (save interval)
    # assumes swap/step is an integer
#    for i in xrange(0, totmoves, swap):
 #           for j in range(numreplicas):
  #                  rep = protein_location[j][i/swap]
   #                 Q_trajec_singleprot[j,i:i+swap] = replicas[rep].nc[i:i+swap]
    k=swap/step
    for i in xrange(len(protein_location[0])):
            for j in range(numreplicas):
		    rep = protein_location[j][i]
                    Q_trajec_singleprot[j,k*i+1:k*(i+1)+1] = replicas[rep].nc[k*i+1:k*(i+1)+1]
    Q_trajec_singleprot[:,0] = totnc
    Q_trajec_singleprot = Q_trajec_singleprot/totnc
    numpy.save('%s/Qtraj_singleprot.npy' % direc, Q_trajec_singleprot)
if verbose:
    for i in range(numreplicas-1):
        print 'Swaps accepted between replica%i and replica%i: %3.2f percent' % (i, i+1, (swapaccepted[i] / float(swapaccepted[i] + swaprejected[i]) * 100))
    print 'Total swaps accepted: %i' % sum(swapaccepted)
    print 'Total swaps rejected: %i' % sum(swaprejected)


t2 = datetime.datetime.now()
print 'Simulation time: '+str(t2-t1)


