
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
parser.add_option("-b", "--writepdb", action="store_true", default=False, help="the output pdb file")
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
pdbfile=options.writepdb
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

coord_nat=moviecoord(coord,transform)


#Get parameters from .param file
angleparam=getangleparam(paramfile,numbeads)
torsparam=gettorsionparam(paramfile,numbeads)
mass=getmass(paramfile[0:-5]+'top',numbeads)

#pregenerate list of interactions
numint=around(scipy.misc.comb(numbeads,2)) # number of interactions
numint= numint - 2*(numbeads-2)-1 # don't count 12 and 13 neighbors

#native LJ parameter getting
nativeparam_n=getnativefix_n(paramfile,numint,numbeads) # [ones and zeros, native epsilon, native sigma]

totnc=sum(nativeparam_n[:,0]) #total native contacts
nsigma2=nativecutoff2*nativeparam_n[:,2]*nativeparam_n[:,2]


#nonnative LJ parameter getting
[nonnativesig,nnepsil]=getLJparam_n(paramfile,numbeads,numint) #[nonnative sigmas for every interaction, epsilon (one value)]
nonnatindex=-1*(nativeparam_n[:,0]-1) # array of ones and zeros
nonnativeparam=column_stack((nonnatindex,nonnativesig)) #[ones and zeros, nonnative sigma]

# set Simulation class variables
Simulation.angleparam=angleparam
Simulation.torsparam=torsparam
Simulation.totnc=totnc
Simulation.nativeparam_n=nativeparam_n
Simulation.nsigma2=nsigma2
Simulation.nonnativeparam=nonnativeparam
Simulation.nnepsil=nnepsil
Simulation.totmoves=totmoves
Simulation.step=step
Simulation.numbeads=numbeads
Simulation.numint=numint
Simulation.wordtemplate=wordtemplate
Simulation.positiontemplate=positiontemplate
Simulation.ATOMlinenum=ATOMlinenum
Simulation.transform=transform
Simulation.coord_nat=coord_nat
Simulation.mass=mass

# put class variables in a dictionary for threading
dict={'numbeads':numbeads,'step':step,'totmoves':totmoves,'numint':numint,'angleparam':angleparam,'torsparam':torsparam,'nativeparam_n':nativeparam_n,'nonnativeparam':nonnativeparam,'nnepsil':nnepsil,'nsigma2':nsigma2,'transform':transform,'coord_nat':coord_nat,'positiontemplate':positiontemplate,'pdbfile':pdbfile,'mass':mass}

# Calculate temperature distribution
if(numreplicas==1):
	if(trange[0]!=trange[1]):
		print 'WARNING: using lower temperature bound for one replica simulation'
	T=[trange[0]]
elif(numreplicas==2):
	T=trange
else:
	T=empty(numreplicas)
	alpha=(trange[1]/trange[0])**(1/float(numreplicas-1))
	T[0]=trange[0]
	for i in range(1,numreplicas):
		T[i]=T[i-1]*alpha

##########temporary#############
#T[1]=295.0
#T[2:7]=arange(310,335,5)
#trange=[330.0,trange[1]]
#alpha=(trange[1]/trange[0])**(1/9.0)
#for i in range(7,numreplicas):
	#T[i]=T[i-1]*alpha
###################################3

if (verbose):
	print 'verbosity is %s' %(str(verbose))
	print 'total number of moves is %d' %(totmoves)
	print 'autocorrelation step size is %d moves' %(step)
	print 'there are %d replicas' %(numreplicas)
	print 'Temperature: '
	print T
	print 'the replica exchange interval is %d steps' %(swap)
	print 'There are %d residues in %s' %(numbeads,filename)
	



#type 1 switches
def tryswap1(Replicas,Swapaccepted,Swaprejected,Whoiswhere):
	if(numreplicas==1):
		pass
	else:
		#i=randint(0,numreplicas-2) #tests swap between replica i and i+1
		for i in xrange(0,numreplicas-1,2):
			P=exp((Replicas[i].u0-Replicas[i+1].u0)*(1/(kb*T[i])-1/(kb*T[i+1])))
			if(random.random()<P):
				Replicas[i].coord, Replicas[i+1].coord = Replicas[i+1].coord, Replicas[i].coord
				Replicas[i].u0, Replicas[i+1].u0 = Replicas[i+1].u0, Replicas[i].u0
				Replicas[i].r2, Replicas[i+1].r2 = Replicas[i+1].r2, Replicas[i].r2
				Replicas[i].torsE, Replicas[i+1].torsE = Replicas[i+1].torsE, Replicas[i].torsE
				Replicas[i].angE, Replicas[i+1].angE = Replicas[i+1].angE, Replicas[i].angE
				#print 'Swapping replica'+str(i)+' and replica'+str(i+1)
				Swapaccepted[i]+=1
				Replicas[i].whoami.append(Replicas[i+1].whoami[-1])
				Replicas[i+1].whoami.append(Replicas[i].whoami[-2])
				Whoiswhere[Replicas[i].whoami[-1]].append(i)
				Whoiswhere[Replicas[i+1].whoami[-1]].append(i+1)
			else:
				#print 'not swapping replica'+str(i)+' and replica'+str(i+1)
				Swaprejected[i]+=1
				Replicas[i].whoami.append(Replicas[i].whoami[-1])
				Replicas[i+1].whoami.append(Replicas[i+1].whoami[-1])
				Whoiswhere[Replicas[i].whoami[-1]].append(Whoiswhere[Replicas[i].whoami[-1]][-1])
				Whoiswhere[Replicas[i+1].whoami[-1]].append(Whoiswhere[Replicas[i+1].whoami[-1]][-1])
	if numreplicas%2==1:
		Replicas[numreplicas-1].whoami.append(Replicas[numreplicas-1].whoami[-1])
		Whoiswhere[Replicas[numreplicas-1].whoami[-1]].append(Whoiswhere[Replicas[numreplicas-1].whoami[-1]][-1])
	return [Swapaccepted,Swaprejected,Whoiswhere]

#type 2 switches
def tryswap2(Replicas,Swapaccepted,Swaprejected,Whoiswhere):
	if(numreplicas==1):
		pass
	else:
		#i=randint(0,numreplicas-2) #tests swap between replica i and i+1
		for i in xrange(1,numreplicas-1,2):
			P=exp((Replicas[i].u0-Replicas[i+1].u0)*(1/(kb*T[i])-1/(kb*T[i+1])))
			if(random.random()<P):
				Replicas[i].coord, Replicas[i+1].coord = Replicas[i+1].coord, Replicas[i].coord
				Replicas[i].u0, Replicas[i+1].u0 = Replicas[i+1].u0, Replicas[i].u0
				Replicas[i].r2, Replicas[i+1].r2 = Replicas[i+1].r2, Replicas[i].r2
				Replicas[i].torsE, Replicas[i+1].torsE = Replicas[i+1].torsE, Replicas[i].torsE
				Replicas[i].angE, Replicas[i+1].angE = Replicas[i+1].angE, Replicas[i].angE
				#print 'Swapping replica'+str(i)+' and replica'+str(i+1)
				Swapaccepted[i]+=1
				Replicas[i].whoami.append(Replicas[i+1].whoami[-1])
				Replicas[i+1].whoami.append(Replicas[i].whoami[-2])
				Whoiswhere[Replicas[i].whoami[-1]].append(i)
				Whoiswhere[Replicas[i+1].whoami[-1]].append(i+1)
			else:
				#print 'not swapping replica'+str(i)+' and replica'+str(i+1)
				Swaprejected[i]+=1
				Replicas[i].whoami.append(Replicas[i].whoami[-1])
				Replicas[i+1].whoami.append(Replicas[i+1].whoami[-1])
				Whoiswhere[Replicas[i].whoami[-1]].append(Whoiswhere[Replicas[i].whoami[-1]][-1])
				Whoiswhere[Replicas[i+1].whoami[-1]].append(Whoiswhere[Replicas[i+1].whoami[-1]][-1])
	Replicas[0].whoami.append(Replicas[0].whoami[-1])
	Whoiswhere[Replicas[0].whoami[-1]].append(Whoiswhere[Replicas[0].whoami[-1]][-1])
	if len(Replicas[numreplicas-1].whoami) != len(Replicas[0].whoami):
		Replicas[numreplicas-1].whoami.append(Replicas[numreplicas-1].whoami[-1])
		Whoiswhere[Replicas[numreplicas-1].whoami[-1]].append(Whoiswhere[Replicas[numreplicas-1].whoami[-1]][-1])
	return [Swapaccepted,Swaprejected,Whoiswhere]

def pprun(Replicas,moves,Dict):
	jobs=[job_server.submit(run,(Replicas[i],moves,Dict),(energy,),("random","numpy","energyfunc","moveset","writetopdb","hmcGO")) for i in range(numreplicas)]
	#jobs=[job_server.submit(run,(moves,Simulations[i].coord,numbeads,step,totmoves,numint,angleparam,torsparam,nativeparam_n,nonnativeparam,nnepsil,nsigma2,transform,coord_nat),(energy,energyprint),('random','numpy','energyfunc','moveset','writetopdb')) for i in range(numreplicas)]
	newReplicas=[job() for job in jobs]
	return newReplicas

def energy(mpos,rsquare,torsE,angE):
	energy=sum(angE)+sum(torsE)+energyfunc.LJenergy_n(rsquare,Simulation.nativeparam_n,Simulation.nonnativeparam,Simulation.nnepsil)
	return energy


#========================================================================================================
# SIMULATE
#========================================================================================================

# instantiate replicas
replicas=[]
for i in range(len(T)):
	direc='./replicaexchange/replica%s' % (str(i))
	if not os.path.exists(direc):
		os.mkdir(direc)
	name='replica%s' % (str(i))
	replicas.append(Simulation(name,direc,coord,T[i]))
	replicas[i].whoami.append(i)
	if (pdbfile):
		mcoord=moviecoord(coord,transform)
		writepdb(mcoord,wordtemplate,ATOMlinenum,0,'%s/trajectory%s.pdb' %(replicas[i].out,str(int(replicas[i].T))))
	
#result=replicas[0].run(1000,replicas[0].coord)
#job_server=pp.Server()
#print job_server
#jobs=job_server.submit(replicas[0].run,(1000,replicas[0].coord,numbeads,step,totmoves,numint,angleparam,torsparam,nativeparam_n,nonnativeparam,nnepsil,nsigma2,transform,coord_nat),(energy,energyprint),modules=('random','numpy','energyfunc','moveset','writetopdb'))
#result=jobs()
#print result

move=0
swapaccepted=zeros(numreplicas-1)
swaprejected=zeros(numreplicas-1)
job_server=pp.Server(ppservers=())
whoiswhere=range(numreplicas)

########SIMULATE##########
for i in whoiswhere:
	whoiswhere[i]=[i]

#for i in range(totmoves/swap):
	#for replica in replicas: run(replica,swap,dict)
	#tryswap1(replicas,swapaccepted,swaprejected,whoiswhere)
	#print i
	

	
for i in xrange(totmoves/swap):
	stdout.write(str(i)+'\r')
	stdout.flush()
	replicas=pprun(replicas,swap,dict)
	job_server.wait()
	if i%2==0:
		#type 1 switch
		[swapaccepted,swaprejected,whoiswhere]=tryswap1(replicas,swapaccepted,swaprejected,whoiswhere)
	else:
		#type 2 switch
		[swapaccepted,swaprejected,whoiswhere]=tryswap2(replicas,swapaccepted,swaprejected,whoiswhere)

#########OUTPUT##########
job_server.print_stats()
savetxt('./replicaexchange/whoiswhere.txt',whoiswhere)
for i in range(len(replicas)):
	#pdb.set_trace()
	replicas[i].output(verbose)
	replicas[i].saveenergy(plot)
	replicas[i].savermsd(plot)
	replicas[i].savenc(plot)
	replicas[i].savehist(plot)
	#print replica.whoami
	plt.figure(5)
	plt.plot(range(len(whoiswhere[i])),whoiswhere[i])
plt.xlabel('swap')
plt.ylabel('replica')
plt.savefig('./replicaexchange/whoiswhere.png')

if(verbose):
	for i in range(numreplicas-1):
		print 'Swaps accepted betweeen replica'+str(i)+' and replica'+str(i+1)+': '+str('%3.2f') %(swapaccepted[i]/float(swapaccepted[i]+swaprejected[i])*100)+'%'
	print 'Total swaps accepted: '+str(sum(swapaccepted))
	print 'Total swaps rejected: '+str(sum(swaprejected))


t2=datetime.now()
#plt.show()
#if (pdbfile != ''):
	#f=open(pdbfile,'a')
	#f.write('END\r\n')
	#f.close
	#print 'wrote trajectory to %s' %(pdbfile)


print('Simulation time: '+str(t2-t1))


