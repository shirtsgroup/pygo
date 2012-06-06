#========================================================================================================
# IMPORTS
#========================================================================================================

from datetime import datetime
t1=datetime.now()
import numpy as np
from writetopdb import *
from moveset import *
from energyfunc import *
from optparse import OptionParser
import matplotlib.pyplot as plt
from sys import stdout
from random import *
import profile
import scipy.misc
import pdb
from hmcGO import *

parser=OptionParser()
parser.add_option("-f", "--files", dest="datafiles", default='GO_protein.pdb', help="protein .pdb file")
parser.add_option("-p", "--parameterfile", dest="paramfiles", default='GO_protein.param', help="protein .param file")
parser.add_option("-d", "--directory", dest="datafile_directory", default='./', help="the directory the data files are in")
parser.add_option("-t", "--temperature", default='300', dest="T", type="float", help="temperature")
parser.add_option("-v", "--verbose", action="store_false", default=True, help="more verbosity")
parser.add_option("-c", "--addconnect", action="store_true", default=False, help="add bonds to .pdb file")
parser.add_option("-n", "--moves", dest="totmoves", type="int", default='100', help="total number of moves")
parser.add_option("-s", "--stepsize", dest="step", type="int", default='100', help="number of moves between save operations")
parser.add_option("-o", "--outputfiles", dest="outputfiles", default='conformation_energies.txt', help="the output file for every [step] energy")
parser.add_option("-g", "--histogram", dest="histname", default='', help="name histogram of conformational energies, if desired")
parser.add_option("-a", "--energyplot", dest="plotname", default='', help="name of plot of accepted conformation energies, if desired")
parser.add_option("-b", "--writepdb", dest="pdbfile", default='', help="the output pdb file")
#parser.add_option("-e", "--percentmove", nargs=2, dest="percentmove", type="float",default=[.33,.66], help="the output pdb file")
parser.add_option("-r", "--rmsd", dest="rmsdfig", default='', help="name of rmsd figure, if desired")

(options,args)=parser.parse_args()
#========================================================================================================
# CONSTANTS
#========================================================================================================


verbose=options.verbose
T=options.T #Kelvin
totmoves=options.totmoves
step=options.step
filename=options.datafile_directory + '/' + options.datafiles
print filename
paramfile=options.datafile_directory + '/' + options.paramfiles
outputfiles=options.outputfiles
histname=options.histname
plotname=options.plotname
pdbfile=options.pdbfile
#percentmove=options.percentmove
rmsdfig=options.rmsdfig
addbonds=options.addconnect


energyarray=zeros(totmoves/step+1)
nc=zeros(totmoves/step+1)
if (rmsdfig != ""):
	rmsd_array=zeros(totmoves/step+1)
	

kb=0.0019872041 #kcal/mol/K
percentmove=[.2,.3,.5] #% bend,% axis torsion,% crankshaft, %local move
#percentmove=[0,0,0]
maxtheta=[10*T/300.,10*T/200.,20*T/250.,5.] # bend, axistorsion, crankshaft
nativecutoff=1.2


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
	addconnect(filename,numbeads)
	
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
if (rmsdfig != ""):
	untransform=getmovietransform(coord)
	transform=transpose(untransform)
	coord_nat=moviecoord(coord,transform)

#Get parameters from .param file
angleparam=getangleparam(paramfile,numbeads)
torsparam=gettorsionparam(paramfile,numbeads)


#speed up terms
numint=around(scipy.misc.comb(numbeads,2)) # number of interactions

numint= numint - 2*(numbeads-2)-1 # don't count 12 and 13 neighbors

#native LJ parameter getting
nativeparam_n=getnativefix_n(paramfile,numint,numbeads) # [ones and zeros, native epsilon, native sigma]
totnc=sum(nativeparam_n[:,0]) #total native contacts
nsigma2=nativecutoff*nativecutoff*nativeparam_n[:,2]*nativeparam_n[:,2]

#nonnative LJ parameter getting
nonnativesig=getLJparam_n(paramfile,numbeads,numint) #[nonnative sigmas for every interaction, epsilon (one value)]
nnepsil=-nonnativesig[-1] # last element
nonnativesig=delete(nonnativesig,-1) #remove epsilon from sigma array
nonnatindex=-1*(nativeparam_n[:,0]-1) # array of ones and zeros
nonnativeparam=column_stack((nonnatindex,nonnativesig)) #[ones and zeros, nonnative sigma]

if (verbose):
	print 'verbosity is %s' %(str(verbose))
	print 'temperature is %f' %(T)
	print 'total number of moves is %d' %(totmoves)
	print 'autocorrelation step size is %d moves' %(step)
	print 'There are %d residues in %s' %(numbeads,filename)

def energyprint(mpos,rsquare,torsE,angE):
	LJ=LJenergy_n(rsquare,nativeparam_n,nonnativeparam,nnepsil)
	energy=sum(angE)+sum(torsE)+LJ
	print 'Angle energy: ' + str(sum(angE))
	print 'Torsion energy: '+str(sum(torsE))
	print 'LJ energy: '+str(LJ)
	return energy

def energy(mpos,rsquare,torsE,angE):
	energy=sum(angE)+sum(torsE)+LJenergy_n(rsquare,nativeparam_n,nonnativeparam,nnepsil)
	return energy

		

r2=getLJr2(coord,numint,numbeads)
torsE=torsionenergy_nn(coord,zeros(numbeads-3),torsparam,arange(numbeads-3))
angE=angleenergy_n(coord,zeros(numbeads-2),angleparam,arange(numbeads-2))
u0=energyprint(coord,r2,torsE,angE)

writepdb(coord,wordtemplate,ATOMlinenum,0,pdbfile)

h=.04
nsteps=100
m=120.368 # average mass of all residues
tol=1e-10



mpos=coord.copy()
vel=np.random.normal(0,sqrt(.0083144621*T/120.3684),(numbeads,3)) #in nm/ps, uses average residue mass
bonds=bond(mpos)
d2=sum(bonds**2,axis=1)
force=bondedforces(mpos,torsparam,angleparam,bonds)+nonbondedforces(mpos,numint,numbeads,nativeparam_n,nonnativeparam,nnepsil)

print 'Potential: '+str(u0)
print 'Kinetic: '+str(.5*sum(vel**2))



for i in range(nsteps):
	
	# finding r(t+dt)
	conv=ones(len(mpos)-1)
	q=vel+h/(2*m)*force
	while sum(conv)!=0:
		for i in range(len(mpos)-1): # loops through all bonds/constraints
			s=bonds[i]+h*(q[i,:]-q[i+1,:])
			diff=sum(s**2)-d2[i]
		#	print diff
			if (abs(diff)<tol):
				conv[i]=0
			else:
				g=diff/(4*h*dot(s,bonds[i]))*m
				q[i,:]-=g/m*bonds[i]
				q[i+1,:]+=g/m*bonds[i]
		print sum(conv)
	mpos+=h*q #r(t+dt)
	vel=q # new v(t)
	bonds=bond(mpos) #rij(t+dt)
	force=bondedforces(mpos,torsparam,angleparam,bonds)+nonbondedforces(mpos,numint,numbeads,nativeparam_n,nonnativeparam,nnepsil)
	
	#finding v(t+dt)
	conv=ones(len(mpos)-1)
	vel+=h/(2*m)*force
	while sum(conv)!=0:
		for i in range(len(mpos)-1):
			if (abs(dot(bonds[i,:],vel[i,:]-vel[i+1,:]))<tol):
				conv[i]=0
			else:
				k=dot(bonds[i,:],(vel[i,:]-vel[i+1,:])/d2[i]*m/2)
				vel[i] -= k/m*bonds[i,:]
				vel[i+1] += k/m*bonds[i,:]
				
	addtopdb(mpos,positiontemplate,i+1,pdbfile)
	stdout.write(str(i)+'\r')
	stdout.flush()
print sum(bonds**2,axis=1)-d2

#for i in range(nsteps):
	#mpos_new=mpos+v*tstep+.5*a*tstep**2
	#bonds=bond(mpos_new)
        #fnew=bondedforces(mpos_new,torsparam,angleparam,bonds)+nonbondedforces(mpos_new,numint,numbeads,nativeparam_n,nonnativeparam,nnepsil)
	##fnew=nonbondedforces(mpos_new,numint,numbeads,nativeparam_n,nonnativeparam,nnepsil)
	#fnew=fnew*4.184 # convert from kcal to kJ
        #anew=fnew/120.368 # use average residue mass
	#vnew=v+.5*tstep*(anew+a)
	#v=vnew
	#a=anew
	#mpos=mpos_new
	#addtopdb(mpos,positiontemplate,i+1,pdbfile)
	##print i

coord=mpos
r2=getLJr2(coord,numint,numbeads)
torsE=torsionenergy_nn(coord,zeros(numbeads-3),torsparam,arange(numbeads-3))
angE=angleenergy_n(coord,zeros(numbeads-2),angleparam,arange(numbeads-2))
u0=energyprint(coord,r2,torsE,angE)
print 'Potential: '+str(u0)
print 'Kinetic: '+str(.5*sum(vel**2))

f=open(pdbfile,'a')
f.write('END\r\n')
f.close

print 'Simulation time: '+str(datetime.now()-t1)
