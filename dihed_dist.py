
#=======================================================================================================
# IMPORTS
#=======================================================================================================

from datetime import datetime
import numpy
from optparse import OptionParser
import matplotlib.pyplot as plt
from sys import stdout
import random
import profile
import scipy.misc
import os
import pp
import pdb

from writetopdb import *
from moveset import *
from energyfunc import *
from simulationobject import *

t1=datetime.now()

parser=OptionParser()
parser.add_option("-f", "--files", dest="datafile", default="GO_protein.pdb", help="protein .pdb file")
parser.add_option("-p", "--parameterfile", dest="paramfile", default='GO_protein.param', help="protein .param file")
parser.add_option("-t", "--temprange", nargs=2, default=[300.0,300.0], type="float", dest="temprange", help="temperature range of replicas")
parser.add_option("-r", "--replicas", default=3, type="int",dest="replicas", help="number of replicas")
parser.add_option("-v", "--verbose", action="store_false", default=True, help="more verbosity")
parser.add_option("-c", "--addconnect", action="store_true", default=False, help="add bonds to .pdb file")
parser.add_option("-n", "--moves", dest="totmoves", type="int", default='100', help="total number of moves")
parser.add_option("-s", "--stepsize", dest="step", type="int", default='100', help="number of moves between save operations")
parser.add_option("-k", "--swapstep", dest="swap", type="int", default='1000', help="number of moves between swap moves")
parser.add_option("-g", "--histogram", dest="histname", default='', help="name histogram of conformational energies, if desired")
parser.add_option("-a", "--plot", action="store_true", default=True, help="plots energy, rmsd, fractional nativeness")
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



# set Simulation class variables
Simulation.totmoves = totmoves
Simulation.step = step
Simulation.numbeads = numbeads


# put class variables in a dictionary for threading
dict = {'numbeads':numbeads, 'step':step, 'totmoves':totmoves}

# Calculate temperature distribution
#if numreplicas == 1:
    #if trange[0] != trange[1]:
        #print 'WARNING: using lower temperature bound for one replica simulation'
    #T = [trange[0]]
#elif numreplicas == 2:
    #T = trange
#else:
    #T = numpy.empty(numreplicas)
    #alpha = (trange[1] / trange[0])**(1 / float(numreplicas - 1))
    #T[0] = trange[0]
    #for i in range(1, numreplicas):
        #T[i] = T[i-1] * alpha
T=range(2)

if (verbose):
    print 'verbosity is %s' %(str(verbose))
    print 'total number of moves is %d' %(totmoves)
    print 'autocorrelation step size is %d moves' %(step)
    print 'there are %d replicas' %(numreplicas)
    print 'Temperature: '
    print T
    print 'the replica exchange interval is %d steps' %(swap)
    print 'There are %d residues in %s' %(numbeads,filename)


def pprun(Replicas, Moves, Dict):
    jobs = [job_server.submit(run_ff, (replica, Moves, Dict), (), ("random", "numpy",
            "energyfunc", "moveset", "writetopdb", "hmcGO")) for replica in Replicas]
    newReplicas = [job() for job in jobs]
    return newReplicas



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
replicas[0].percentmove=[0,0,0,0] # anglebend, axistorsion, crankshaft, localmove
replicas[1].percentmove=[0,0,0,0]
replicas[0].a=0 # no jacobian
replicas[0].b=1
replicas[1].a=1 # yes jacobian
replicas[1].b=0
replicas[0].dihedarr=array(dihedral(coord))
replicas[1].dihedarr=array(dihedral(coord))
replicas[0].angarr=array(angle(coord))
replicas[1].angarr=array(angle(coord))

    #if pdbfile:
        #mcoord = moviecoord(coord, transform)
        #writepdb(mcoord, wordtemplate, ATOMlinenum, 0, '%s/trajectory%i.pdb' % (replicas[i].out, int(replicas[i].T)))
	
#result=replicas[0].run(1000,replicas[0].coord)
#job_server=pp.Server()
#print job_server
#jobs=job_server.submit(replicas[0].run,(1000,replicas[0].coord,numbeads,step,totmoves,numint,angleparam,torsparam,nativeparam_n,nonnativeparam,nnepsil,nsigma2,transform,coord_nat),(energy,energyprint),modules=('random','numpy','energyfunc','moveset','writetopdb'))
#result=jobs()
#print result

move = 0
job_server = pp.Server(ppservers=())

########SIMULATE##########


	
for i in xrange(totmoves/swap):
    stdout.write(str(i) + '\r')
    stdout.flush()
    replicas = pprun(replicas, swap, dict)
    job_server.wait()

#########OUTPUT##########
job_server.print_stats()

for i in range(len(replicas)):
    plt.figure(i+1)
    plt.hist(replicas[i].dihedarr,40,histtype='barstacked')
    plt.xlabel('dihedrals')
    plt.savefig('%s/dihedralhist%i.png' % (replicas[i].out, i))
    plt.figure(i+5)
    plt.hist(replicas[i].angarr,40,histtype='barstacked')
    plt.xlabel('angle')
    plt.savefig('%s/anglehist%i.png' % (replicas[i].out, i))
#plt.show()

t2 = datetime.now()
print 'Simulation time: '+str(t2-t1)


