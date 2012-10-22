#=======================================================================================================
# IMPORTS
#=======================================================================================================

from datetime import datetime
import numpy
from optparse import OptionParser
try:
	import matplotlib.pyplot as plt
except:
	plot=False
from sys import stdout
import random
import profile
import scipy.linalg
from scipy.misc import comb
import os
import pp
import pdb

from writetopdb import *
from moveset import *
from energyfunc import *
#from simulationobject import *

t1=datetime.now()

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
parser.add_option("-a", "--plot", action="store_false", default=True, help="plots energy, rmsd, fractional nativeness")
parser.add_option("-b", "--writepdb", action="store_true", default=False, help="the output pdb file")
parser.add_option("-e", "--percentmove", nargs=1, dest="e", type="int",default=100, help="the output pdb file")
parser.add_option("--id", nargs=1, dest="id", type="int", default=0, help="the simlog id number")
parser.add_option("--freq", nargs=4, dest="freq", type="float", default=[.2,.4,.7,1.], help="move frequencies")
parser.add_option("--surf", action="store_true", default=False, help="surface simulation flag")
parser.add_option("--umbrella", type="float", default=0., help="umbrella simulation flag")

(options,args)=parser.parse_args()

#======================================================================================================
# CONSTANTS
#======================================================================================================
verbose = options.verbose
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
percentmove = options.freq

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
coord_nat = moviecoord(coord,transform)


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
Simulation.coord_nat = coord_nat
Simulation.mass = mass
Simulation.percentmove = percentmove

# put class variables in a dictionary for pprun
dict = {'percentmove':percentmove,'numbeads':numbeads, 'step':step, 'totmoves':totmoves, 'numint':numint, 'angleparam':angleparam, 'torsparam':torsparam, 'nativeparam_n':nativeparam_n, 'nonnativeparam':nonnativeparam, 'nnepsil':nnepsil, 'nsigma2':nsigma2, 'transform':transform, 'coord_nat':coord_nat, 'positiontemplate':positiontemplate, 'pdbfile':pdbfile, 'mass':mass}

# Calculate temperature distribution
if numreplicas == 1:
    if trange[0] != trange[1]:
        print 'WARNING: using lower temperature bound for one replica simulation'
    T = [trange[0]]
elif numreplicas == 2:
    T = trange
else:
    T = numpy.empty(numreplicas)
    alpha = (trange[1] / trange[0])**(1 / float(numreplicas - 1))
    T[0] = trange[0]
    for i in range(1, numreplicas):
        T[i] = T[i-1] * alpha

if (verbose):
    output = ['verbosity is %s' %(str(verbose)),
    'total number of moves is %d' %(totmoves),
    'save interval is %d moves' %(step),
    'there are %d replicas' %(numreplicas),
    'Temperature:',
    str(T),
    'the replica exchange interval is %d steps' %(swap),
    'there are %i swaps at each exchange point' %(e),
    'There are %d residues in %s' %(numbeads,filename),
    'percent move is',
    str(Simulation.percentmove),'']
    print "\r\n".join(output)


if surf:
    xlength = 135
    ylength = 135
    spacing = 10
    yspacing = spacing*3.**.5
    surface = getsurf(xlength+15,ylength+15,spacing)
    writesurf('surface.pdb',surface)
    nsurf = len(surface)
    nspint = nsurf*numbeads # surface-protein interactions
    surfparam = getsurfparam('%spdb' % (paramfile[3:-5]), numbeads, nsurf, nspint)
    scale = 1
    surfparam[:,0] = surfparam[:,0]*scale
    print 'Surface energy parameters scaled by %f' % scale
    SurfaceSimulation.surface = surface
    SurfaceSimulation.nsurf = nsurf
    SurfaceSimulation.nspint = nspint
    SurfaceSimulation.surfparam = surfparam
    dict.update({'nspint':nspint, 'nsurf':nsurf, 'surfparam':surfparam, 'surface':surface, 'xlength':xlength, 'ylength':ylength, 'spacing':spacing, 'yspacing':yspacing})
#======================================================================================================
# SUBROUTINES
#======================================================================================================

#type 1 switches
def tryswap1(Replicas, Swapaccepted, Swaprejected, Whoiswhere):
    if numreplicas == 1:
        pass
    else:
        for i in xrange(0, numreplicas-1, 2):
            P=numpy.exp((Replicas[i].u0-Replicas[i+1].u0)*(1/(kb*T[i])-1/(kb*T[i+1])))
            if(random.random()<P):
                Swapaccepted[i]+=1
                Replicas[i].coord, Replicas[i+1].coord = Replicas[i+1].coord, Replicas[i].coord
                Replicas[i].u0, Replicas[i+1].u0 = Replicas[i+1].u0, Replicas[i].u0
                Replicas[i].r2, Replicas[i+1].r2 = Replicas[i+1].r2, Replicas[i].r2
                Replicas[i].torsE, Replicas[i+1].torsE = Replicas[i+1].torsE, Replicas[i].torsE
                Replicas[i].angE, Replicas[i+1].angE = Replicas[i+1].angE, Replicas[i].angE
                Replicas[i].whoami, Replicas[i+1].whoami = Replicas[i+1].whoami, Replicas[i].whoami
                Whoiswhere[i], Whoiswhere[i+1] = Whoiswhere[i+1], Whoiswhere[i]
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
    return Swapaccepted, Swaprejected, Whoiswhere

#type 2 switches
def tryswap2(Replicas, Swapaccepted, Swaprejected, Whoiswhere):
    if numreplicas == 1:
        pass
    else:
        for i in xrange(1,numreplicas-1, 2):
            P=numpy.exp((Replicas[i].u0-Replicas[i+1].u0)*(1/(kb*T[i])-1/(kb*T[i+1])))
            if(random.random()<P):
                Swapaccepted[i]+=1
                Replicas[i].coord, Replicas[i+1].coord = Replicas[i+1].coord, Replicas[i].coord
                Replicas[i].u0, Replicas[i+1].u0 = Replicas[i+1].u0, Replicas[i].u0
                Replicas[i].r2, Replicas[i+1].r2 = Replicas[i+1].r2, Replicas[i].r2
                Replicas[i].torsE, Replicas[i+1].torsE = Replicas[i+1].torsE, Replicas[i].torsE
                Replicas[i].angE, Replicas[i+1].angE = Replicas[i+1].angE, Replicas[i].angE
                Replicas[i].whoami, Replicas[i+1].whoami = Replicas[i+1].whoami, Replicas[i].whoami
                if surf:
                    Replicas[i].surfE, Replicas[i+1].surfE = Replicas[i+1].surfE, Replicas[i].surfE
		if umbrella:
		    Replicas[i].z_array, Replicas[i+1].z_array = Replicas[i+1].z_array, Replicas[i].z_array

                Whoiswhere[i], Whoiswhere[i+1] = Whoiswhere[i+1], Whoiswhere[i]
	        #Whoiswhere[Replicas[i].whoami].append(i)
                #Whoiswhere[Replicas[i+1].whoami].append(i+1)
            else:
                Swaprejected[i] += 1
                #Whoiswhere[Replicas[i].whoami].append(Whoiswhere[Replicas[i].whoami][-1])
                #Whoiswhere[Replicas[i+1].whoami].append(Whoiswhere[Replicas[i+1].whoami][-1])
    #Whoiswhere[Replicas[0].whoami].append(Whoiswhere[Replicas[0].whoami][-1])
    #if numreplicas%2 == 0:
        #Whoiswhere[Replicas[numreplicas-1].whoami].append(Whoiswhere[Replicas[numreplicas-1].whoami][-1])
    return Swapaccepted, Swaprejected, Whoiswhere

def tryrepeatedswaps(Replicas, Swapaccepted, Swaprejected, protein_location, transmat):
	if numreplicas == 1:
		return Swapaccepted, Swaprejected, protein_location, transmat
	replica_location = numpy.arange(numreplicas)
	for k in range(e): # try swapping 100 times
		if random.random() < .5:
			Swapaccepted, Swaprejected, replica_location = tryswap1(Replicas, Swapaccepted, Swaprejected, replica_location)
		else:	
			Swapaccepted, Swaprejected, replica_location = tryswap2(Replicas, Swapaccepted, Swaprejected, replica_location)
	replica_location_indexed = numpy.zeros(numreplicas)
	for i in range(numreplicas):
		replica_location_indexed[replica_location[i]] = i
		protein_location[Replicas[i].whoami].append(i)
	for i in range(numreplicas):
		transmat[i,replica_location_indexed[i]] += 1
	return Swapaccepted, Swaprejected, protein_location, transmat

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
                "energyfunc", "moveset", "writetopdb")) for replica in Replicas]
            newReplicas = [job() for job in jobs]
        return newReplicas
else:
    def pprun(Replicas, Moves, Dict):
        if len(Replicas) == 1:
            newReplicas = [run_surf(Replicas[0], Moves, Dict)]
        else:
            jobs = [job_server.submit(run_surf, (replica, Moves, Dict), (update_energy, save), ("random", "numpy",
                "energyfunc", "moveset", "writetopdb")) for replica in Replicas]
            newReplicas = [job() for job in jobs]
        return newReplicas
#def energy(mpos, rsquare, torsE, angE):
    #"""Returns the total energy of a configurat"""
    #energy = sum(angE) + sum(torsE) + energyfunc.LJenergy_n(rsquare, Simulation.nativeparam_n,
             #Simulation.nonnativeparam, Simulation.nnepsil)
    #return energy


#=======================================================================================================
# SIMULATE
#=======================================================================================================
move = 0
swapaccepted = numpy.zeros(numreplicas-1)
swaprejected = numpy.zeros(numreplicas-1)
protein_location = [[i] for i in range(numreplicas)]
transmat = [zeros((numreplicas,numreplicas)) for i in range(10)]
if swap!=totmoves:
    assert(totmoves / swap % 10 == 0)
    transmat_index = totmoves/swap/10

# instantiate replicas
replicas = []
direc = './replicaexchange/simlog%i' % id
if not os.path.exists(direc):
    os.mkdir(direc)


for i in range(len(T)):
        name = 'replica%i' % i
	if umbrella:
		replicas.append(UmbrellaSimulation(name, direc, coord, T[i], surface, umbrella, mass))
	elif surf:
		replicas.append(SurfaceSimulation(name, direc, coord, T[i], surface))
	else:
        	replicas.append(Simulation(name, direc, coord, T[i]))
        replicas[i].whoami = i
        if pdbfile:
            mcoord = moviecoord(coord, transform)
            writepdb(mcoord, wordtemplate, ATOMlinenum, 0, '%s/trajectory%i.pdb' % (replicas[i].out, int(replicas[i].T)))



try:
	# running on the cluster
	f = open('nodefile.txt','r')
	ppservers = f.readline()
	ppservers = ppservers.split(' ')
	ppservers[-1] = ppservers[-1][0:-1]
	ppservers = [server+':23335' for server in ppservers]
	ppservers = tuple(ppservers)
	print ppservers
	job_server = pp.Server(ppservers=ppservers)
except:
	# running on one machine
	job_server = pp.Server(ppservers=())

ti = datetime.now()
if umbrella:
	print 'Starting umbrella simulation...'
elif surf:
	print 'Starting surface simulation...'
else:
	print 'Starting simulation...'
for i in xrange(totmoves/swap):
    replicas = pprun(replicas, swap, dict)
    job_server.wait()
    if swap!=totmoves:
    	swapaccepted, swaprejected, protein_location, transmat[i/transmat_index] = tryrepeatedswaps(replicas, swapaccepted, swaprejected, protein_location, transmat[i/transmat_index])
    t_remain = (datetime.now() - ti)/(i+1)*(totmoves/swap - i - 1)
    stdout.write(str(t_remain) + '\r')
    stdout.flush()


#=======================================================================================================
# POST ANALYSIS
#=======================================================================================================
if swap!=totmoves:
    #calculate mixing time
    mixtime = zeros(10)
    fulltransmat = zeros((numreplicas,numreplicas))
    for i in range(10):
            fulltransmat += transmat[i]
            transmat[i] = transmat[i]/transmat_index #normalize to probabilities
            transmatrix = numpy.matrix(transmat[i])
            eig = scipy.linalg.eigvals(transmatrix)
            eig = numpy.sort(eig)
            print '-----Transition Matrix '+str(i)+' Eigenvalues-----'
            print eig
            mixtime[i] = 1/(1-eig[-2])
            print mixtime[i]
    fulltransmat = fulltransmat/totmoves*swap
    try:
            fulltransmat = numpy.matrix(fulltransmat)
            eig = scipy.linalg.eigvals(fulltransmat)
            eig = numpy.sort(eig)
            print '-----Full Transition Matrix Eigenvalues-----'
            print eig
            print 'Full Tmat miximg time: ' + str(1/(1-eig[-2]))
            print 'Average mixing time from submatrices: ' + str(average(mixtime))
            print 'Standard deviation in mixing time: ' + str(std(mixtime))
    except:
            pdb.set_trace()
#=======================================================================================================
# OUTPUT
#=======================================================================================================
job_server.print_stats()
for i in range(len(replicas)):
    replicas[i].output(verbose)
    replicas[i].saveenergy(plot)
    replicas[i].savermsd(plot)
    replicas[i].savenc(plot)
    replicas[i].savehist(plot)
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
    
    for i in xrange(len(protein_location[0])):
            for j in range(numreplicas):
		    rep = protein_location[j][i]
                    Q_trajec_singleprot[j,swap*i+1:swap*(i+1)+1] = replicas[rep].nc[swap*i+1:swap*(i+1)+1]
    Q_trajec_singleprot[:,0] = totnc
    Q_trajec_singleprot = Q_trajec_singleprot/totnc
    numpy.savetxt('%s/Qtraj_singleprot.txt' % direc, Q_trajec_singleprot)
    numpy.savetxt('%s/protein_location.txt' % direc, numpy.array(protein_location))
         
if verbose:
    for i in range(numreplicas-1):
        print 'Swaps accepted between replica%i and replica%i: %3.2f percent' % (i, i+1, (swapaccepted[i] / float(swapaccepted[i] + swaprejected[i]) * 100))
    print 'Total swaps accepted: %i' % sum(swapaccepted)
    print 'Total swaps rejected: %i' % sum(swaprejected)


t2 = datetime.now()
print 'Simulation time: '+str(t2-t1)


