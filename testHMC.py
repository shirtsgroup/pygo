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
from HMCforce import *

parser=OptionParser()
parser.add_option("-f", "--files", dest="datafiles", default='GO_1PGB.pdb', help="protein .pdb file")
parser.add_option("-p", "--parameterfile", dest="paramfiles", default='GO_1PGB.param', help="protein .param file")
parser.add_option("-d", "--directory", dest="datafile_directory", default='./', help="the directory the data files are in")
parser.add_option("-t", "--temperature", default='300', dest="T", type="float", help="temperature")
parser.add_option("-v", "--verbose", action="store_false", default=True, help="more verbosity")
parser.add_option("-c", "--addconnect", action="store_true", default=False, help="add bonds to .pdb file")
parser.add_option("-n", "--moves", dest="totmoves", type="int", default='100', help="total number of moves")
parser.add_option("-s", "--stepsize", dest="step", type="float", default='100', help="number of moves between save operations")
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
        coord[i][2]=float(line[46:54])
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
[nonnativesig,nnepsil]=getLJparam_n(paramfile,numbeads,numint) #[nonnative sigmas for every interaction, epsilon (one value)]
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

		

r2=cgetLJr2(coord,numint,numbeads)
torsE=torsionenergy_nn(coord,zeros(numbeads-3),torsparam,arange(numbeads-3))
angE=angleenergy_n(coord,zeros(numbeads-2),angleparam,arange(numbeads-2))
u0=energyprint(coord,r2,torsE,angE)

if pdbfile:
    writepdb(coord,wordtemplate,ATOMlinenum,0,pdbfile)

h=step
nsteps=totmoves
#m=120.368 # average mass of all residues
m=getmass('GO_1PGB.top',numbeads)
tol=1e-12
maxloop=1000

mpos=coord.copy()
vel=empty((numbeads,3))
for i in range(numbeads):
	vel[i,:]=np.random.normal(0,sqrt(4.184*kb*T/m[i]),3) #in nm/ps, uses average residue mass
bonds=bond(mpos)
d2=sum(bonds**2,axis=1)
d=d2**.5
#force=bondedforces(mpos,torsparam,angleparam,bonds,d2,d,numbeads)+nonbondedforces(mpos,numint,numbeads,nativeparam_n,nonnativeparam,nnepsil)
force=cangleforces(coord,angleparam,bonds,d,numbeads)+cdihedforces(torsparam,bonds,d2,d,numbeads)+nonbondedforces(mpos,numint,numbeads,nativeparam_n,nonnativeparam,nnepsil)

a=transpose(force)/m
vel, conv = HMCforce.crattle(bonds, vel, m, d2, maxloop, numbeads, tol)
print vel
print 'Potential: '+str(u0)
ke=.5*m/4.184*sum(vel**2,axis=1)
print 'Kinetic: '+str(sum(ke))
h0=u0+sum(ke)
print 'Hamiltonian: '+str(h0)
pot=zeros(nsteps+1)
pot[0]=u0
K=zeros(nsteps+1)
K[0]=sum(ke)

for e in range(nsteps):
    #pdb.set_trace()
    v_half = vel + h / 2 * numpy.transpose(a) # unconstrained v(t+dt/2)
    #v_half, conv = HMCforce.shake(bonds, v_half, h, m, d2, maxloop, numbeads, tol) # constrained v(t+dt/2)
    v_half, conv = HMCforce.cshake(bonds, v_half, h, m, d2, maxloop, numbeads, tol)
    if not conv:
            print 'not converging'
            break
    mpos += h * v_half #constrained r(t+dt)
    bonds = mpos[0:numbeads-1,:]-mpos[1:numbeads,:] #rij(t+dt)
    force=cangleforces(mpos,angleparam,bonds,d,numbeads)+cdihedforces(torsparam,bonds,d2,d,numbeads)+nonbondedforces(mpos,numint,numbeads,nativeparam_n,nonnativeparam,nnepsil)
    #force = HMCforce.bondedforces(mpos, torsparam, angleparam, bonds, d2, d, numbeads) +	HMCforce.nonbondedforces(mpos, numint, numbeads, nativeparam_n, nonnativeparam, nnepsil)
    a = transpose(force)/m
    vel = v_half + h/2*transpose(a)
    vel, conv = HMCforce.crattle(bonds, vel, m, d2, maxloop, numbeads, tol)
    if not conv:
            print 'not converging'
            break
    if pdbfile:
        addtopdb(mpos,positiontemplate,i+1,pdbfile)
    stdout.write(str(e)+'\r')
    stdout.flush()
    
    coord=mpos
    r2=cgetLJr2(coord,numint,numbeads)
    torsE=ctorsionenergy(coord,zeros(numbeads-3),torsparam,arange(numbeads-3))
    angE=cangleenergy(coord,zeros(numbeads-2),angleparam,arange(numbeads-2))
    u1=energy(coord,r2,torsE,angE)
    ke=.5*m/4.184*sum(vel**2,axis=1)
    
    
    pot[e+1]=u1
    K[e+1]=sum(ke)
	
print sum(bonds**2,axis=1)-d2

coord=mpos
r2=cgetLJr2(coord,numint,numbeads)
torsE=torsionenergy_nn(coord,zeros(numbeads-3),torsparam,arange(numbeads-3))
angE=angleenergy_n(coord,zeros(numbeads-2),angleparam,arange(numbeads-2))
u1=energyprint(coord,r2,torsE,angE)
print '----------------------------'
print 'Potential: '+str(u1)
ke=.5*m/4.184*sum(vel**2,axis=1)
print 'Kinetic: '+str(sum(ke))
h1=u1+sum(ke)
print '----------------------------'
print 'Hamiltonian: '+str(h1)

print 'Difference in H ' + str(h1-h0)
print 'H Boltzmann Factor ' + str(exp(-(h1-h0)/kb/T))

print '----------------------------'
print 'Difference in U ' + str(u1-u0)
print 'U Boltzmann Factor ' + str(exp(-(u1-u0)/kb/T))
print '----------------------------'
print 'Variance in H ' + str(std(pot+K))
coeff = polyfit(arange(nsteps+1), pot+K, 1)
print 'Slope of line fit to H ' + str(coeff[0])
print 'Variance / Slope ' + str(std(pot+K)/coeff[0])
#f=open(pdbfile,'a')
#f.write('END\r\n')
#f.close

print 'Simulation time: '+str(datetime.now()-t1)

plt.figure(1)
plt.plot(range(nsteps+1),K,label='KE')
plt.plot(range(nsteps+1),pot,label='PE')
plt.plot(range(nsteps+1),K+pot,label='Total Energy')
plt.legend()
plt.xlabel('step (h=%f)'%(h))
plt.ylabel('Energy (kcal/mol)')
plt.title('MD simulation at %i K for %f ps' %(T, h*nsteps))

plt.figure(2)
plt.plot(range(nsteps+1),pot+K)
plt.show()
