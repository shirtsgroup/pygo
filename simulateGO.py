
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
energyarray=zeros(totmoves/step+1)
filename=options.datafile_directory + '/' + options.datafiles
paramfile=options.datafile_directory + '/' + options.paramfiles
outputfiles=options.outputfiles
histname=options.histname
plotname=options.plotname
pdbfile=options.pdbfile
#percentmove=options.percentmove
rmsdfig=options.rmsdfig
addbonds=options.addconnect



if (rmsdfig != ""):
	rmsd_array=zeros(totmoves/step+1)
	

kb=0.0019872041 #kcal/mol/K
percentmove=[.33,.66] #% bend,% axis torsion,% crankshaft
maxtheta=[10*T/300.,10*T/300.,20*T/300.] # bend, axistorsion, crankshaft

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
numint=int(scipy.misc.comb(numbeads,2)+1) # number of interactions
numint= numint - 2*(numbeads-2)-1 # don't count 12 and 13 neighbors

#native LJ parameter getting
nativeparam_n=getnativefix_n(paramfile,numint,numbeads) # [ones and zeros, native epsilon, native sigma]

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

#========================================================================================================
# SIMULATE
#========================================================================================================
r2=getLJr2(coord,numint,numbeads)
torsE=torsionenergy_nn(coord,zeros(numbeads-3),torsparam,arange(numbeads-3))
angE=angleenergy_n(coord,zeros(numbeads-2),angleparam,arange(numbeads-2))
u0=energyprint(coord,r2,torsE,angE)
print u0
energyarray[0]=u0
if (rmsdfig != ""):
	rmsd_array[0]=rmsd(coord_nat,coord_nat)

if(pdbfile != ''):
	untransform=getmovietransform(coord)
	transform=transpose(untransform)
	mcoord=moviecoord(coord,transform)
	#writeseqpdb(mcoord,wordtemplate,ATOMlinenum,0)
	writepdb(mcoord,wordtemplate,ATOMlinenum,0,pdbfile)
	print 'writing trajectory to %s...' %(pdbfile)

# constants for move stats
angmoves=0
crankmoves=0
atormoves=0
acceptedat=0
accepteda=0
acceptedc=0
accepted=0
rejected=0
movetype=''
move=0

while move<totmoves:
        randmove=random()
	randdir=random()
	m=randint(1,numbeads-2) #random bead, not end ones
	
	#bend
	if randmove < percentmove[0]:
	    theta=maxtheta[0]/180.*pi-random()*maxtheta[0]*pi/180.*2
            newcoord=anglebend(coord,m,randdir,theta)
	    angmoves += 1
	    movetype='a'
	    change=[]
	    angchange=[m-1]
	    #if randdir<.5:
		    #change=arange(m-2,m+2)
		    #angchange=arange(m-1,m+2)
		    #if m==1:
			    #change=arange(0,3)
		    #elif m==numbeads-4: #53
			    #change=arange(m-2,m+1)
		    #elif m>numbeads-4: #53
			    #change=change[change<(numbeads-3)]
			    #angchange=angchange[angchange<(numbeads-2)]
	    #else:
		    #change=arange(m-4,m)
		    #angchange=arange(m-3,m)
		    #if m==3:
			    #change=arange(0,m)
		    #elif m<3:
			    #change=change[change>-1]
			    #angchange=angchange[angchange>-1]
		    #elif m==numbeads-2:
			    #change=arange(m-4,m-1)

	#axis torsion
	elif randmove < percentmove[1]:
	    theta=maxtheta[1]/180.*pi-random()*pi*maxtheta[1]/180.*2
	    newcoord=axistorsion(coord,m,randdir,theta)
	    movetype='at'
	    atormoves += 1
	    angchange=[]
	    if randdir<.5:
		    change=[m-2]
		    if m<2:
			    change=[]
	    elif m==numbeads-2:
		    change=[]
	    else:
		    change=[m-1]
	
	#crankshaft
	else:
	    theta=maxtheta[2]/180.*pi-random()*maxtheta[2]*pi/180.*2
	    newcoord=crankshaft(coord,m,theta)
	    crankmoves += 1
	    movetype='c'
	    change=arange(m-3,m+1)
	    angchange=[m-2,m]
	    if m==2:
		    change=arange(m-2,m+1)
	    elif m<2:
		    change=change[change>-1]
		    angchange=[m]
 	    elif m==numbeads-4 or m==numbeads-3:
		    change=change[change<(numbeads-3)]
	    elif m==numbeads-2:
		    change=arange(m-3,m-1)
		    angchange=[m-2]

	r2new=getLJr2(newcoord,numint,numbeads)
	newangE=angleenergy_n(newcoord,angE,angleparam,angchange)
	newtorsE=torsionenergy_nn(newcoord,torsE,torsparam,change)
	u1=energy(newcoord,r2new,newtorsE,newangE)
        move += 1
	stdout.write(str(move)+'\r')
	stdout.flush()
        boltz=exp(-(u1-u0)/(kb*T))
	if u1< u0 or random() < boltz:
        	accepted += 1
		if movetype=='a':
			accepteda += 1
    		elif movetype=='r':
			acceptedr += 1
    		elif movetype=='at':
			acceptedat +=1
			#accepted_angle.append(theta)
    		else: 
        		acceptedc += 1
    		if (pdbfile != ''):
			mcoord=moviecoord(newcoord,transform)
			#writeseqpdb(mcoord,wordtemplate,ATOMlinenum,accepted)
			addtopdb(mcoord,positiontemplate,move,pdbfile)
		r2=r2new
		coord=newcoord
		torsE=newtorsE
		angE=newangE
    		u0=u1
	else:
	    rejected += 1
	if move%step==0:
	    energyarray[move/step]=u0
	    if (rmsdfig != ''):
	    	mcoord=moviecoord(coord,transform)
	    	rmsd_array[move/step]=rmsd(coord_nat,mcoord)
t2=datetime.now()

if (pdbfile != ''):
	f=open(pdbfile,'a')
	f.write('END\r\n')
	f.close
	print 'wrote trajectory to %s' %(pdbfile)


#========================================================================================================
# OUTPUT
#========================================================================================================

savetxt(outputfiles,energyarray)
print 'wrote every %d conformation energies to %s' %(step,outputfiles)

#savetxt('rmsd.txt',rmsd_array)
#print 'wrote rmsd array to rmsd.txt'

if plotname != '':
	print 'generating conformational energy plot...'
	plt.figure(1)
	plt.plot(range(len(energyarray)),energyarray)
	plt.xlabel('move/%d' %(step))
	plt.ylabel('energy (kcal/mol)')
	plt.title('Go-like model monte carlo simulation at '+str(T)+' K')
	plt.savefig(plotname)
	print 'conformational energy plot saved to %s' %(plotname)

if rmsdfig != '':
	print 'generating rmsd plot...'
	plt.figure(1)
	plt.plot(range(len(rmsd_array)),rmsd_array)
	plt.xlabel('move/%d' %(step))
	plt.ylabel('rmsd')
	plt.title('Go-like model monte carlo simulation at '+str(T)+' K')
	plt.savefig(rmsdfig)
	print 'rmsd plot saved to %s' %(rmsdfig)

#if rmsdname != '':
	#print 'generating rmsd plot...'
	#plt.figure(3)
	#plt.plot(range(len(rmsdname)),rmsd)

if histname != '':
	print 'generating conformation energy histogram'
	plt.figure(2)
	plt.hist(energyarray,40)
	plt.title('conformation energies at %f Kelvin taken every %d moves' %(T,step))
	plt.xlabel('conformation energies (kcal/mol)')
	plt.savefig(histname)
	print 'conformation energy histogram saved to %s' %(histname)

if(verbose):
	print 'total accepted moves: %d' %(accepted)
	print 'total rejected moves: %d' %(rejected)
	print 'angle bend: %d moves accepted out of %d tries' %(accepteda,angmoves)
	print 'crankshaft: %d moves accepted out of %d tries' %(acceptedc,crankmoves)
	print 'torsion: %d moves accepted out of %d tries' %(acceptedat,atormoves)
	print('Simulation time: '+str(t2-t1))

#x=range(len(rmsd_array))
#plt.plot(x,rmsd_array)
#plt.ylabel('RMSD')
#plt.xlabel('move/100')
#plt.title('RMSD from native conformation at %f Kelvin taken every %d moves' %(T,step))
#plt.savefig(rmsdname)
