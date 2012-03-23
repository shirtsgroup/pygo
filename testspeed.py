
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
from random import *
import profile
import scipy.misc
import pdb

parser=OptionParser()
parser.add_option("-f", "--files", dest="datafiles", default='GO_protein.pdb', help="protein .pdb file")
parser.add_option("-p", "--parameterfile", dest="paramfiles", default='GO_protein.param', help="protein .param file")
parser.add_option("-d", "--directory", dest="datafile_directory", default='./', help="the directory the data files are in")
parser.add_option("-t", "--temperature", default='300', dest="T", type="float", help="temperature")
parser.add_option("-v", "--verbose", action="store_false", default=True, help="more verbosity")
parser.add_option("-c", "--addconnect", action="store_true", default=False, help="add bonds to .pdb file")
parser.add_option("-n", "--moves", dest="totmoves", type="int", default='100', help="total number of moves")
parser.add_option("-s", "--stepsize", dest="step", type="int", default='50', help="number of moves between each sample")
parser.add_option("-o", "--outputfiles", dest="outputfiles", default='conformation_energies.txt', help="the output file for every [step] energy")
parser.add_option("-g", "--histogram", dest="histname", default='', help="name histogram of conformational energies, if desired")
parser.add_option("-a", "--energyplot", dest="plotname", default='', help="name of plot of accepted conformation energies, if desired")
parser.add_option("-b", "--writepdb", action="store_true", default=False, help="write pdbs")

(options,args)=parser.parse_args()
#========================================================================================================
# CONSTANTS
#========================================================================================================

verbose=options.verbose
T=options.T #Kelvin
totmoves=options.totmoves
step=options.step
energyarray=zeros(totmoves/step+1)
rmsd_array=zeros(totmoves/step+1)
filename=options.datafile_directory + '/' + options.datafiles
paramfile=options.datafile_directory + '/' + options.paramfiles
outputfiles=options.outputfiles
histname=options.histname
plotname=options.plotname
writepdb=options.writepdb

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
coord_nat=coord.copy()


#Get parameters from .param file
angleparam=getangleparam(paramfile,numbeads)
torsparam=gettorsionparam(paramfile,numbeads)
LJparam=getLJparam(paramfile,numbeads)
nativeparam=getnativefix(paramfile)

if (verbose):
	print 'verbosity is %s' %(str(verbose))
	print 'temperature is %f' %(T)
	print 'total number of moves is %d' %(totmoves)
	print 'autocorrelation step size is %d moves' %(step)
	print 'There are %d residues in %s' %(numbeads,filename)


###############################
######## SPEED UP TERMS########
###############################
#numint=scipy.misc.comb(numbeads,2) # number of interactions
#numint= numint - 110 # don't count 12 and 13 neighbors
##native LJ parameter getting
#nativeparam_n=getnativefix_n(paramfile,numint,numbeads) # [ones and zeros, native epsilon, native sigma]
##nonnative LJ parameter getting
#nonnativesig=getLJparam_n(paramfile,numbeads,numint) #[nonnative sigmas for every interaction, epsilon (one value)]

#nnepsil=-nonnativesig[-1] # last element
#nonnativesig=delete(nonnativesig,-1) #remove epsilon from sigma array
#nonnatindex=-1*(nativeparam_n[:,0]-1) # array of ones and zeros
#nonnativeparam=column_stack((nonnatindex,nonnativesig)) #[ones and zeros, nonnative sigma]

##temporary
def energy(mpos,torsE,change):
	oldtorsE=torsionenergy_n(mpos,torsE,torsparam,change)
	newtorsE=torsionenergy_nn(mpos,torsE,torsparam,change)
	energyold=sum(oldtorsE)
	energynew=sum(newtorsE)
	return [newtorsE,energyold,energynew]


## calculate square distances for each interaction
#r2=getLJr2(coord,numint,numbeads)
#r2_new=getLJr2_n(coord,numint,numbeads)


#print energy(coord,r2)-energy(coord,r2_new)


move=0
torsE=zeros(54)
change=arange(54)
[torsE,old,new]=energy(coord,torsE,change)
print old-new

#while move<totmoves:
	#randdir=random()
	#m=randint(4,52)
	#coord=torsion(coord,m,randdir)
	#if randdir<.5:
		#change=arange(m-1,m+2)
	#else:
		#change=arange(m-3,m)
	#[newangE,o,n]=energy(coord,angE,change)
	#print o-n
	#move +=1
	#angE=newangE
	
#while move<totmoves:
	#m=randint(4,52)
	#coord=crankshaft(coord,m)
	#change=[m-2,m]
	#[newangE,o,n]=energy(coord,angE,change)
	#print o-n
	#move +=1
	#angE=newangE

#while move<totmoves:
	#randdir=random()
	#m=randint(4,52)
	#coord=axistorsion(coord,m,randdir)
	#if randdir<.5:
		#change=[]
	#else:
		#change=[]
	#[newangE,o,n]=energy(coord,angE,change)
	#print o-n
	#move +=1
	#angE=newangE
	
#dihedarr=dihedral(coord,numbeads)
#newcoord=axistorsion(coord,5,.7)
#dihedarr_n=dihedral(newcoord,numbeads)
#print dihedarr_n-dihedarr

#angles=angle(coord)
#newcoord=torsion(coord,5,.7)
#newangle=angle(newcoord)
#print len(angles)
#print newangle-angles

torsE=zeros(54)
change=arange(54)
[torsE,old,new]=energy(coord,torsE,change)
t1=datetime.now()
for i in range(10000):
	oldtorsE=torsionenergy_n(coord,torsE,torsparam,change)
	energyold=sum(oldtorsE)
t2=datetime.now()
for i in range(10000):
	newtorsE=torsionenergy_nn(coord,torsE,torsparam,change)
	energynew=sum(newtorsE)
print datetime.now()-t2
print t2-t1
