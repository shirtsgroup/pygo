
#========================================================================================================
# IMPORTS
#========================================================================================================

from datetime import datetime
t1=datetime.now()
from numpy import *
from writetopdb import *
from moveset import *
from energyfunc import *
from optparse import OptionParser
import matplotlib.pyplot as plt
from sys import stdout
import random
import profile
import scipy.misc
from simulationobject import *
import os
import pp
import pdb

parser=OptionParser()
parser.add_option("-f", "--files", dest="datafile", default="GO_protein.pdb", help="protein .pdb file")
parser.add_option("-p", "--parameterfile", dest="paramfile", default='GO_protein.param', help="protein .param file")
parser.add_option("-t", "--temprange", nargs=2, default=[300.0,300.0], type="float", dest="temprange", help="temperature range of replicas")
parser.add_option("-r", "--replicas", default=1, type="int",dest="replicas", help="number of replicas")
parser.add_option("-v", "--verbose", action="store_false", default=True, help="more verbosity")
parser.add_option("-c", "--addconnect", action="store_true", default=False, help="add bonds to .pdb file")
parser.add_option("-n", "--moves", dest="totmoves", type="int", default='100', help="total number of moves")
parser.add_option("-s", "--stepsize", dest="step", type="int", default='100', help="number of moves between save operations")
parser.add_option("-k", "--swapstep", dest="swap", type="int", default='1000', help="number of moves between swap moves")
parser.add_option("-g", "--histogram", dest="histname", default='', help="name histogram of conformational energies, if desired")
parser.add_option("-a", "--plot", action="store_true", default=True, help="plots energy, rmsd, fractional nativeness")
parser.add_option("-b", "--writepdb", dest="pdbfile", default='', help="the output pdb file")
#parser.add_option("-e", "--percentmove", nargs=2, dest="percentmove", type="float",default=[.33,.66], help="the output pdb file")

(options,args)=parser.parse_args()
#========================================================================================================
# CONSTANTS
#========================================================================================================


verbose=options.verbose
trange=options.temprange#Kelvin
numreplicas=options.replicas
totmoves=options.totmoves
step=options.step
swap=options.swap
filename=options.datafile
paramfile=options.paramfile
plot=options.plot
histname=options.histname
pdbfile=options.pdbfile
#percentmove=options.percentmove
addbonds=options.addconnect

kb=0.0019872041 #kcal/mol/K
nativecutoff2=1.2**2

# read .pdb to get number of beads (numbeads)
file=open(filename,'r')
numbeads=0
while 1:
	line=file.readline()
	if not line:
		break
	splitline=line.split('  ')
	if splitline[0]=='ATOM':
		numbeads += 1


if (addbonds):
	addconnect(pdbfile,numbeads)
	
# gets bead coordinates from .pdb file
file.seek(0)
coord=empty((numbeads,3)) #3 dimensional
i=0 # index for coordinate matrix
k=0 # index for line number
wordtemplate=[] #template of all words
positiontemplate=[] #template of position lines
ATOMlinenum=[] # array of location of ATOM lines
while 1:
    line = file.readline()
    if not line:
        break
    wordtemplate.append(line)
    splitline=line.split('  ')
    if splitline[0]=='ATOM':
        positiontemplate.append(line)
        ATOMlinenum.append(k)
        coord[i][0]=float(line[31:38])
        coord[i][1]=float(line[39:46])
        coord[i][2]=float(line[47:54])
        i=i+1
    k=k+1
file.close()



#Get native conformation
#if (rmsdfig != ""):

untransform=getmovietransform(coord)
transform=transpose(untransform)
Simulation.transform=transform
coord_nat=moviecoord(coord,transform)
Simulation.coord_nat=coord_nat

#Get parameters from .param file
angleparam=getangleparam(paramfile,numbeads)
torsparam=gettorsionparam(paramfile,numbeads)
Simulation.angleparam=angleparam
Simulation.torsparam=torsparam

#pregenerate list of interactions
numint=int(scipy.misc.comb(numbeads,2)+1) # number of interactions
numint= numint - 2*(numbeads-2)-1 # don't count 12 and 13 neighbors

#native LJ parameter getting
nativeparam_n=getnativefix_n(paramfile,numint,numbeads) # [ones and zeros, native epsilon, native sigma]
Simulation.nativeparam_n=nativeparam_n
Simulation.totnc=sum(nativeparam_n[:,0]) #total native contacts
nsigma2=nativecutoff2*nativeparam_n[:,2]*nativeparam_n[:,2]
Simulation.nsigma2=nsigma2

#nonnative LJ parameter getting
nonnativesig=getLJparam_n(paramfile,numbeads,numint) #[nonnative sigmas for every interaction, epsilon (one value)]
nnepsil=-nonnativesig[-1] # last element
nonnativesig=delete(nonnativesig,-1) #remove epsilon from sigma array
nonnatindex=-1*(nativeparam_n[:,0]-1) # array of ones and zeros
nonnativeparam=column_stack((nonnatindex,nonnativesig)) #[ones and zeros, nonnative sigma]
Simulation.nonnativeparam=nonnativeparam
Simulation.nnepsil=nnepsil
# set additional Simulation class variables
Simulation.totmoves=totmoves
Simulation.step=step
Simulation.numbeads=numbeads

Simulation.numint=numint
Simulation.wordtemplate=wordtemplate
Simulation.positiontemplate=positiontemplate
Simulation.ATOMlinenum=ATOMlinenum

# Calculate temperature distribution
if(numreplicas==1):
	if(trange[0]!=trange[1]):
		print 'WARNING: using lower temperature bound for one replica simulation'
	T=[trange[0]]
elif(numreplicas==2):
	T=trange
elif(numreplicas==3):
	T=empty(numreplicas)
	alpha=(trange[1]/trange[0])**(1/float(numreplicas-1))
	T[0]=trange[0]
	for i in range(1,numreplicas):
		T[i]=T[i-1]*alpha

if (verbose):
	print 'verbosity is %s' %(str(verbose))
	print 'total number of moves is %d' %(totmoves)
	print 'autocorrelation step size is %d moves' %(step)
	print 'there are %d replicas' %(numreplicas)
	print 'Temperature: '
	print T
	print 'the replica exchange interval is %d steps' %(swap)
	print 'There are %d residues in %s' %(numbeads,filename)


def tryswap(Replicas):
	if(numreplicas==1):
		pass
	else:
		i=randint(0,numreplicas-2) #tests swap between replica i and i+1
		P=exp((Replicas[i].u0-Replicas[i+1].u0)*(1/(kb*T[i])-1/(kb*T[i+1])))
		print P
		if(random.random()<P):
			Replicas[i].coord, Replicas[i+1].coord = Replicas[i+1].coord, Replicas[i].coord
			Replicas[i].u0, Replicas[i+1].u0 = Replicas[i+1].u0, Replicas[i].u0
			Replicas[i].r2, Replicas[i+1].r2 = Replicas[i+1].r2, Replicas[i].r2
			Replicas[i].torsE, Replicas[i+1].torsE = Replicas[i+1].torsE, Replicas[i].torsE
			Replicas[i].angE, Replicas[i+1].angE = Replicas[i+1].angE, Replicas[i].angE
			print 'Swapping replica'+str(i)+' and replica'+str(i+1)
		else:
			print 'not swapping replica'+str(i)+' and replica'+str(i+1)

def pprun(Replicas,moves):
	job_server=pp.Server()
	jobs=[job_server.submit(run,(Replicas[i],moves,),(),("random","numpy","energyfunc","moveset","writetopdb")) for i in range(numreplicas)]
	for job in jobs: print job()
	#jobs=[job_server.submit(run,(moves,Simulations[i].coord,numbeads,step,totmoves,numint,angleparam,torsparam,nativeparam_n,nonnativeparam,nnepsil,nsigma2,transform,coord_nat),(energy,energyprint),('random','numpy','energyfunc','moveset','writetopdb')) for i in range(numreplicas)]
	#lastenergy=[job() for job in jobs]
	#return lastenergy

#def energyprint(mpos,rsquare,torsE,angE):
	#LJ=energyfunc.LJenergy_n(rsquare,Simulation.nativeparam_n,Simulation.nonnativeparam,Simulation.nnepsil)
	#energy=sum(angE)+sum(torsE)+LJ
	#print 'Angle energy: ' + str(sum(angE))
	#print 'Torsion energy: '+str(sum(torsE))
	#print 'LJ energy: '+str(LJ)
	#print 'Total energy: '+str(energy)
	#return energy

#def energy(mpos,rsquare,torsE,angE):
	#energy=sum(angE)+sum(torsE)+energyfunc.LJenergy_n(rsquare,Simulation.nativeparam_n,Simulation.nonnativeparam,Simulation.nnepsil)
	#return energy


#========================================================================================================
# SIMULATE
#========================================================================================================

# instantiate replicas
replicas=[]
for i in range(len(T)):
	direc='./replicaexchange/replica'+str(i)
	if not os.path.exists(direc):
		os.mkdir(direc)
	name='replica'+str(i)
	replicas.append(Simulation(name,direc,coord,T[i]))

#result=replicas[0].run(1000,replicas[0].coord)
#job_server=pp.Server()
#print job_server
#jobs=job_server.submit(replicas[0].run,(1000,replicas[0].coord,numbeads,step,totmoves,numint,angleparam,torsparam,nativeparam_n,nonnativeparam,nnepsil,nsigma2,transform,coord_nat),(energy,energyprint),modules=('random','numpy','energyfunc','moveset','writetopdb'))
#result=jobs()
#print result

move=0
########SIMULATE##########
for i in range(totmoves/swap):
	for replica in replicas: run(replica,swap)
	tryswap(replicas)

#job_server=pp.Server()
#for i in range(totmoves/swap):
	##pdb.set_trace()
	#pprun(replicas,swap)
	#tryswap(replicas)

#########OUTPUT##########
for replica in replicas:
	#pdb.set_trace()
	replica.output(verbose)
	replica.saveenergy(plot)
	replica.savermsd(plot)
	replica.savenc(plot)
	replica.savehist(plot)

#job_server=pp.Server()
#move=0
#print move

#for i in range(1):#range(totmoves/swap):
	#currentenergy=pprun(replicas,swap)
	#tryswap(replicas,currentenergy)
	#move += totmoves/swap
	#print move

#e=pprun(replicas,totmoves%swap)
#print e
#move += totmoves%swap
#print move

#output=''
#for replica in replicas:
	#output+=replica.output(verbose)
	#replica.saveenergy(plot)
	#replica.savermsd(plot)
	#replica.savenc(plot)
#print output
#f=open('./replicaexchange/output.txt','w')
#f.write(output)
#f.close	
	
	
t2=datetime.now()

#if (pdbfile != ''):
	#f=open(pdbfile,'a')
	#f.write('END\r\n')
	#f.close
	#print 'wrote trajectory to %s' %(pdbfile)


##========================================================================================================
## OUTPUT
##========================================================================================================

#fraction=nc/totnc
#print 'excluding first '+str(len(fraction)/5)+' native contact values from average'
#fraction=fraction[len(fraction)/20:-1]
#fractionfile='fractionnative'+str(int(T))+'.txt'
#savetxt(fractionfile,fraction)
#print 'wrote ever %d conformation fraction native contacts to %s' %(step,fractionfile)
#print sum(fraction)/len(fraction)


print('Simulation time: '+str(t2-t1))

#x=range(len(rmsd_array))
#plt.plot(x,rmsd_array)
#plt.ylabel('RMSD')
#plt.xlabel('move/100')
#plt.title('RMSD from native conformation at %f Kelvin taken every %d moves' %(T,step))
#plt.savefig(rmsdname)
