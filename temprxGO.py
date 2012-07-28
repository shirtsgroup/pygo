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
from scipy.misc import comb
import os
import pp
import pdb

from writetopdb import *
from moveset import *
from energyfunc import *
from simulationobject import *

t1=datetime.now()

parser=OptionParser()
parser.add_option("-f", "--files", dest="datafile", default="GO_1PGB.pdb", help="protein .pdb file")
parser.add_option("-p", "--parameterfile", dest="paramfile", default='GO_1PGB.param', help="protein .param file")
parser.add_option("-t", "--temprange", nargs=2, default=[300.0,300.0], type="float", dest="temprange", help="temperature range of replicas")
parser.add_option("-r", "--replicas", default=1, type="int",dest="replicas", help="number of replicas")
parser.add_option("-v", "--verbose", action="store_false", default=True, help="more verbosity")
parser.add_option("-c", "--addconnect", action="store_true", default=False, help="add bonds to .pdb file")
parser.add_option("-n", "--moves", dest="totmoves", type="int", default='100', help="total number of moves")
parser.add_option("-s", "--stepsize", dest="step", type="int", default='100', help="number of moves between save operations")
parser.add_option("-k", "--swapstep", dest="swap", type="int", default='1000', help="number of moves between swap moves")
parser.add_option("-g", "--histogram", dest="histname", default='', help="name histogram of conformational energies, if desired")
parser.add_option("-a", "--plot", action="store_false", default=True, help="plots energy, rmsd, fractional nativeness")
parser.add_option("-b", "--writepdb", action="store_true", default=False, help="the output pdb file")
#parser.add_option("-e", "--percentmove", nargs=2, dest="percentmove", type="float",default=[.33,.66], help="the output pdb file")

(options,args)=parser.parse_args()

#======================================================================================================
# CONSTANTS
#======================================================================================================
verbose = options.verbose
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
#percentmove=options.percentmove
addbonds = options.addconnect

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

# put class variables in a dictionary for threading
dict = {'numbeads':numbeads, 'step':step, 'totmoves':totmoves, 'numint':numint, 'angleparam':angleparam, 'torsparam':torsparam, 'nativeparam_n':nativeparam_n, 'nonnativeparam':nonnativeparam, 'nnepsil':nnepsil, 'nsigma2':nsigma2, 'transform':transform, 'coord_nat':coord_nat, 'positiontemplate':positiontemplate, 'pdbfile':pdbfile, 'mass':mass}

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
    print 'verbosity is %s' %(str(verbose))
    print 'total number of moves is %d' %(totmoves)
    print 'autocorrelation step size is %d moves' %(step)
    print 'there are %d replicas' %(numreplicas)
    print 'Temperature: '
    print T
    print 'the replica exchange interval is %d steps' %(swap)
    print 'There are %d residues in %s' %(numbeads,filename)
    print 'percent move is'
    print Simulation.percentmove

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
                Whoiswhere[Replicas[i].whoami].append(i)
                Whoiswhere[Replicas[i+1].whoami].append(i+1)
            else:
                Swaprejected[i]+=1
                Whoiswhere[Replicas[i].whoami[-1]].append(Whoiswhere[Replicas[i].whoami[-1]][-1])
                Whoiswhere[Replicas[i+1].whoami[-1]].append(Whoiswhere[Replicas[i+1].whoami[-1]][-1])
    if numreplicas%2 == 1:
        Whoiswhere[Replicas[numreplicas-1].whoami].append(Whoiswhere[Replicas[numreplicas-1].whoami][-1])
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
                Whoiswhere[Replicas[i].whoami].append(i)
                Whoiswhere[Replicas[i+1].whoami].append(i+1)
            else:
                Swaprejected[i] += 1
                Whoiswhere[Replicas[i].whoami].append(Whoiswhere[Replicas[i].whoami][-1])
                Whoiswhere[Replicas[i+1].whoami].append(Whoiswhere[Replicas[i+1].whoami][-1])
    Whoiswhere[Replicas[0].whoami].append(Whoiswhere[Replicas[0].whoami][-1])
    if numreplicas%2 == 0:
        Whoiswhere[Replicas[numreplicas-1].whoami].append(Whoiswhere[Replicas[numreplicas-1].whoami][-1])
    return Swapaccepted, Swaprejected, Whoiswhere

def tryrepeatedswaps(Replicas, Swapaccepted, Swaprejected, Whoiswhere):
	if numreplicas == 1:
		return Swapaccepted, Swaprejected, Whoiswhere
	for k in range(50): # try swapping 100 times
		Swapaccepted, Swaprejected, Whoiswhere = tryswap1(Replicas, Swapaccepted, Swaprejected, Whoiswhere)
		Swapaccepted, Swaprejected, Whoiswhere = tryswap2(Replicas, Swapaccepted, Swaprejected, Whoiswhere)
	return Swapaccepted, Swaprejected, Whoiswhere

#def update_energy(self, torschange, angchange):
    #self.newtorsE = energyfunc.ctorsionenergy(self.newcoord, self.torsE, Simulation.torsparam, torschange)
    #self.newangE = energyfunc.cangleenergy(self.newcoord, self.angE, Simulation.angleparam, angchange)
    #self.u1 = sum(self.newtorsE)+sum(self.newangE) + energyfunc.cgetLJenergy(self.newcoord, Simulation.numint, Simulation.numbeads, Simulation.nativeparam_n, Simulation.nonnativeparam, Simulation.nnepsil)
    #return self

def pprun(Replicas, Moves, Dict):
    if len(Replicas) == 1:
        newReplicas = [run(Replicas[0], Moves, Dict)]
    else:
        jobs = [job_server.submit(run, (replica, Moves, Dict), (update_energy, save), ("random", "numpy",
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

# instantiate replicas
replicas = []
for i in range(len(T)):
    direc = './replicaexchange/replica%i' % i
    if not os.path.exists(direc):
        os.mkdir(direc)
    name = 'replica%i' % i
    replicas.append(Simulation(name, direc, coord, T[i]))
    replicas[i].whoami = i
    if pdbfile:
        mcoord = moviecoord(coord, transform)
        writepdb(mcoord, wordtemplate, ATOMlinenum, 0, '%s/trajectory%i.pdb' % (replicas[i].out, int(replicas[i].T)))
	

move = 0
swapaccepted = numpy.zeros(numreplicas-1)
swaprejected = numpy.zeros(numreplicas-1)

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


whoiswhere = [[i] for i in range(numreplicas)]

########SIMULATE##########
#for i in whoiswhere:
    #whoiswhere[i] = [i]

#for i in range(totmoves/swap):
	#for replica in replicas: run(replica,swap,dict)
	#tryswap1(replicas,swapaccepted,swaprejected,whoiswhere)
	#print i
	

ti = datetime.now()
print 'Starting simulation...'
for i in xrange(totmoves/swap):
    replicas = pprun(replicas, swap, dict)
    job_server.wait()
    swapaccepted, swaprejected, whoiswhere = tryrepeatedswaps(replicas, swapaccepted, swaprejected, whoiswhere)
    #if i%2 == 0:
        ##type 1 switch
        #[swapaccepted, swaprejected] = tryswap1(replicas, swapaccepted, swaprejected)
    #else:
        ##type 2 switch
        #[swapaccepted, swaprejected] = tryswap2(replicas, swapaccepted, swaprejected)
    t_remain = (datetime.now() - ti)/(i+1)*(totmoves/swap - i - 1)
    stdout.write(str(t_remain) + '\r')
    stdout.flush()


#########OUTPUT##########
job_server.print_stats()
numpy.savetxt('./replicaexchange/whoiswhere.txt', whoiswhere)
for i in range(len(replicas)):
    replicas[i].output(verbose)
    replicas[i].saveenergy(plot)
    replicas[i].savermsd(plot)
    replicas[i].savenc(plot)
    replicas[i].savehist(plot)
    #print replica.whoami
    plt.figure(5)
    plt.plot(range(len(whoiswhere[i])), whoiswhere[i])
plt.xlabel('swap')
plt.ylabel('replica')
plt.savefig('./replicaexchange/whoiswhere.png')

if verbose:
    for i in range(numreplicas-1):
        print 'Swaps accepted between replica%i and replica%i: %3.2f percent' % (i, i+1, (swapaccepted[i] / float(swapaccepted[i] + swaprejected[i]) * 100))
    print 'Total swaps accepted: %i' % sum(swapaccepted)
    print 'Total swaps rejected: %i' % sum(swaprejected)


t2 = datetime.now()
print 'Simulation time: '+str(t2-t1)


