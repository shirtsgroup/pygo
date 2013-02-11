#!/usr/bin/python

# todo -- simplify the total energy read in, the kinetic energy read-in, temperature read-in
#=========================================================
#IMPORTS
#=========================================================
import sys
import numpy
import pymbar # for MBAR analysis
import timeseries # for timeseries analysis
import os
import os.path
import pdb  # for debugging
from optparse import OptionParser

parser=OptionParser()
parser.add_option("-t", "--temprange", nargs=2, default=[300.0,450.0], type="float", dest="temprange", help="temperature range of replicas")
parser.add_option("-r", "--replicas", default=8, type="int",dest="replicas", help="number of replicas")
parser.add_option("--direc", dest="datafile", default="replicaexchange/simlog0/", help="Qtraj_singleprot.txt file location")

(options,args) = parser.parse_args()
trange = options.temprange
numreplicas = options.replicas
direc = options.datafile

#========================================================
# CONSTANTS
#=========================================================
#simulation = sys.argv[1]   #This program takes the argument of the simulation top directory
kB = 0.00831447/4.184  #Boltzmann constant (Gas constant) in kJ/(mol*K)
#TE_COL_NUM = 11  #The column number of the total energy in ener_box#.output

#NumTemps = 24          # Last TEMP # + 1 (start counting at 1)
NumIterations = 500  # The number of energies to be taken and analyzed, starting from the last
                  # Extra data will be ignored
dT = 1.25              # Temperature increment for calculating Cv(T)

###########################################################
# For Cv vs T    _____
#               /     \
#  ____________/       \____________
#
############################################################

#=========================================================
# SUBROUTINES
#=========================================================

def read_total_energies(pathname):
       """Reads in the TEMP#/ener_box#.output file and parses it, returning an array of energies 
       ARGUMENTS
             filename (string) - the path to the folder of the simulation
       """

       print "--Reading total energies from %s/..." % pathname


       # Initialize Return variables
       E_kn = numpy.zeros([NumTemps,NumIterations], numpy.float64)
       
       #Read files
       for k in range(NumTemps):
             #Construct each TEMP#/ener_box#.output name and read in the file
             filename = os.path.join(pathname,'TEMP'+str(k), 'ener_box'+str(k)+'.output')
             infile = open(filename, 'r')
             lines = infile.readlines()
             infile.close()
             numLines = len(lines)

             #Initialize arrays for E
             E_from_file = numpy.zeros(NumIterations, numpy.float64)      
       
             #Parse lines in each file
             for n in range(NumIterations):
                  m = numLines - 2 - n #Count down (the 2 is for index purposes(1) and to not use the double-counted last line (1))
                  elements = lines[m].split()
                  E_from_file[n] = float(elements[TE_COL_NUM])
             
             #Add in the E's for each timestep (n) at this temperature (k)
             E_kn[k] = E_from_file;
       return E_kn   

def read_simulation_temps(pathname):
       """Reads in the various temperatures from each TEMP#/simul.output file by knowing
             beforehand the total number of temperatures (parameter at top)
       """
       
       print "--Reading temperatures from %s/..." % pathname       

       # Initialize return variable
       temps_from_file = numpy.zeros(NumTemps, numpy.float64)

       for k in range(NumTemps):
             infile = open(os.path.join(pathname,'TEMP'+ str(k), 'simul'+str(k)+'.output'), 'r')
             lines = infile.readlines()
             infile.close()

             l = len(lines)-1
             test = 'test'
             while (test != 'Temperature_av:'):
                  line = lines[l].split()
                  if (len(line) > 1) : test = line[1]
                  l = l-1 
             temps_from_file[k] = float(line[2])
       
       return temps_from_file

#========================================================================
# MAIN
#========================================================================

#------------------------------------------------------------------------
# Read Data From File
#------------------------------------------------------------------------
print("")
print("Preparing data:")
#T_from_file = read_simulation_temps(simulation)
#E_from_file = read_total_energies(simulation)
#K = len(T_from_file)

#trange = [300.0, 450.0]
#numreplicas = 24
T = numpy.empty(numreplicas)
alpha = (trange[1]/trange[0])**(1/float(numreplicas-1))
T[0]=trange[0]
for i in range(1,numreplicas):
	T[i]=T[i-1]*alpha
print T
files=[]
surf=[]
for i in range(len(T)):
	files.append(direc+'/energy'+str(int(T[i]))+'.txt')

	#files.append('replicaexchange/replica'+str(i)+'/energy'+str(int(T[i]))+'.txt')
#	files.append('replicaexchange/simlog99/energy'+str(int(T[i]))+'.txt')
#	files.append('surface_replica_exchange/replica'+str(i)+'/energy'+str(int(T[i]))+'.txt')

nc=numpy.loadtxt(files[0])

for i,file in enumerate(files):
	nctemp=numpy.loadtxt(file)
#	ncsurf=numpy.loadtxt(surf[i])
	#nctemp -= ncsurf # uncomment this line to get protein-only energies
	nc=numpy.vstack((nc,nctemp))
nc=nc[1:numreplicas+1,:]
temporary = True
if temporary:
	nc = nc[:,::10]
nc = nc[:,-NumIterations:-1]
T_from_file = T
E_from_file = nc.copy()
K = numreplicas
files=[]
for i in range(len(T)):
	files.append(direc+'/fractionnative'+str(int(T[i]))+'.txt')
#	files.append('replicaexchange/simlog99/fractionnative'+str(int(T[i]))+'.txt')

	#files.append('replicaexchange/replica'+str(i)+'/fractionnative'+str(int(T[i]))+'.txt')
#	files.append('surface_replica_exchange/replica'+str(i)+'/fractionnative'+str(int(T[i]))+'.txt')


nc=numpy.loadtxt(files[0])

for file in files:
	nctemp=numpy.loadtxt(file)
	nc=numpy.vstack((nc,nctemp))
nc=nc[1:numreplicas+1,:]
nc = nc[:,-NumIterations:-1]

N_k = numpy.zeros(K,numpy.int32)
g = numpy.zeros(K,numpy.float64)
for k in range(K):  # subsample the energies
   g[k] = timeseries.statisticalInefficiency(E_from_file[k])#,suppress_warning=True)
   indices = numpy.array(timeseries.subsampleCorrelatedData(E_from_file[k],g=g[k])) # indices of uncorrelated samples
   N_k[k] = len(indices) # number of uncorrelated samples
   E_from_file[k,0:N_k[k]] = E_from_file[k,indices]
   nc[k,0:N_k[k]] = nc[k,indices]

#------------------------------------------------------------------------
# Insert Intermediate T's and corresponding blank U's and E's
#------------------------------------------------------------------------
# Set up variables
Temp_k = T_from_file
currentT = T_from_file[0] + dT
maxT = T_from_file[len(T_from_file) - 1]
i = 1

print("--Inserting intermediate temperatures...")

# Loop, inserting T's at which we are interested in the properties
while (currentT < maxT) :
       if (currentT < Temp_k[i]):
             Temp_k = numpy.insert(Temp_k, i, currentT)
             currentT = currentT + dT
       else:
             currentT = Temp_k[i] + dT
             i = i + 1
             
# Update number of states
K = len(Temp_k)

print("--Inserting blank energies to match up with inserted temperatures...")

# Loop, inserting E's into blank matrix (leaving blanks only where new Ts are inserted)

Nall_k = numpy.zeros([K], numpy.int32) # Number of samples (n) for each state (k) = number of iterations/energies
E_kn = numpy.zeros([K, NumIterations], numpy.float64)
Q_kn = numpy.zeros([K, NumIterations], numpy.float64)
i = 0

for k in range(K):
       if (Temp_k[k] == T_from_file[i]):
             E_kn[k,0:N_k[i]] = E_from_file[i,0:N_k[i]]
	     Q_kn[k,0:N_k[i]] = nc[i,0:N_k[i]]
             Nall_k[k] = N_k[i]
             i = i + 1
#------------------------------------------------------------------------
# Compute inverse temperatures
#------------------------------------------------------------------------
beta_k = 1 / (kB * Temp_k)

#------------------------------------------------------------------------
# Compute reduced potential energies
#------------------------------------------------------------------------

print "--Computing reduced energies..."

u_kln = numpy.zeros([K,K,NumIterations], numpy.float64) # u_kln is reduced pot. ener. of segment n of temp k evaluated at temp l

for k in range(K):
       for l in range(K):
             u_kln[k,l,0:Nall_k[k]] = beta_k[l] * E_kn[k,0:Nall_k[k]]

#------------------------------------------------------------------------
# Initialize MBAR
#------------------------------------------------------------------------

# Initialize MBAR with Newton-Raphson
print ""
print "Initializing MBAR:"
print "--K = number of Temperatures"
print "--L = number of Temperatures"
print "--N = number of Energies per Temperature"

# Use Adaptive Method (Both Newton-Raphson and Self-Consistent, testing which is better)
mbar = pymbar.MBAR(u_kln, Nall_k, method = 'adaptive', verbose=True, relative_tolerance=1e-12)

#------------------------------------------------------------------------
# Compute Expectations for E_kt and E2_kt as E_expect and E2_expect
#------------------------------------------------------------------------
print ""
print "Computing Expectations for E..."
(E_expect, dE_expect) = mbar.computeExpectations(u_kln)*(beta_k)**(-1)
#print "Computing Expectations for E^2..."
#(E2_expect,dE2_expect) = mbar.computeExpectations(u_kln*u_kln)*(beta_k)**(-2)

print "Computing Expectations for Q..."
(Q,dQ) = mbar.computeExpectations(Q_kn)

#------------------------------------------------------------------------
# Compute Cv for NVT simulations as <E^2> - <E>^2 / (RT^2)
#------------------------------------------------------------------------
#print ""
#print "Computing Heat Capacity as ( <E^2> - <E>^2 ) / ( R*T^2 )..."

#Cv_expect = numpy.zeros([K], numpy.float64)
#dCv_expect = numpy.zeros([K], numpy.float64)

#for i in range(K):
#       Cv_expect[i] = (E2_expect[i] - (E_expect[i]*E_expect[i])) / ( kB * Temp_k[i] * Temp_k[i])
#       dCv_expect[i] = 2*dE_expect[i]**2 / (kB *Temp_k[i]*Temp_k[i])   # from propagation of error

#print "Temperature  dA         <E> +/- d<E>       <E^2> +/- d<E^2>       Cv +/- dCv"     
#print "-------------------------------------------------------------------------------"
#for k in range(K):
#       print "%8.3f %8.3f %9.3f +/- %5.3f  %9.1f +/- %5.1f   %7.4f +/- %6.4f" % (Temp_k[k],mbar.f_k[k],E_expect[k],dE_expect[k],E2_expect[k],dE2_expect[k],Cv_expect[k], dCv_expect[k])
#numpy.savetxt('/home/edz3fz/Qsurf_int.txt',Q)
#numpy.savetxt('/home/edz3fz/dQsurf_int.txt',dQ)
#numpy.savetxt('/home/edz3fz/dQsol.txt',dQ)

#numpy.savetxt('/home/edz3fz/Qtemp.tt',Temp_k)
import matplotlib.pyplot as plt
ncavg = numpy.average(nc, axis=1)

plt.figure(1)
plt.plot(T, ncavg, 'ko')
plt.plot(Temp_k,Q,'k')
plt.errorbar(Temp_k, Q, yerr=dQ)
plt.xlabel('Temperature (K)')
plt.ylabel('Q fraction native')
#plt.title('Heat Capacity from Go like model MC simulation of 1BSQ')
plt.savefig(direc+'/foldingcurve.png')
plt.show()
