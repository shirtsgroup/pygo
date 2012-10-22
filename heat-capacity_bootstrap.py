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

#========================================================
# CONSTANTS
#=========================================================
#simulation = sys.argv[1]   #This program takes the argument of the simulation top directory
kB = 0.00831447/4.184  #Boltzmann constant (Gas constant) in kJ/(mol*K)
#TE_COL_NUM = 11  #The column number of the total energy in ener_box#.output

NumTemps = 8          # Last TEMP # + 1 (start counting at 1)
NumIterations = 1000  # The number of energies to be taken and analyzed, starting from the last
                  # Extra data will be ignored
dT = 2.5              # Temperature increment for calculating Cv(T)
nBoots = 50

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

trange = [300.0, 450.0]
numreplicas = 8
T = numpy.empty(numreplicas)
alpha = (trange[1]/trange[0])**(1/float(numreplicas-1))
T[0]=trange[0]
for i in range(1,numreplicas):
	T[i]=T[i-1]*alpha
print T
files=[]
for i in range(len(T)):
	files.append('results/1PGB/surface/simlog4/energy'+str(int(T[i]))+'.txt')
#	files.append('surface_replica_exchange/replica'+str(i)+'/energy'+str(int(T[i]))+'.txt')


nc=numpy.loadtxt(files[0])

for file in files:
	nctemp=numpy.loadtxt(file)
	nc=numpy.vstack((nc,nctemp))
nc=nc[1:numreplicas+1,:]
nc = nc[:,-1000:-1]
T_from_file = T
E_from_file = nc.copy()
K = numreplicas


N_k = numpy.zeros(K,numpy.int32)
g = numpy.zeros(K,numpy.float64)
for k in range(K):  # subsample the energies
   g[k] = timeseries.statisticalInefficiency(E_from_file[k])
   indices = numpy.array(timeseries.subsampleCorrelatedData(E_from_file[k],g=g[k])) # indices of uncorrelated samples
   N_k[k] = len(indices) # number of uncorrelated samples
   E_from_file[k,0:N_k[k]] = E_from_file[k,indices]

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
i = 0

for k in range(K):
       if (Temp_k[k] == T_from_file[i]):
             E_kn[k,0:N_k[i]] = E_from_file[i,0:N_k[i]]
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
E_kn_samp = numpy.zeros([K,NumIterations], numpy.float64) # u_kln is reduced pot. ener. of segment n of temp k evaluated at temp l

allCv_expect = numpy.zeros([K,2,nBoots], numpy.float64)
allE_expect = numpy.zeros([K,nBoots], numpy.float64)
allE2_expect = numpy.zeros([K,nBoots], numpy.float64)

for n in range(nBoots):
       print "Bootstrap: %d/%d" % (n,nBoots)
       for k in range(K):
       # resample the results:
             if Nall_k[k] > 0:
                  booti=numpy.random.randint(Nall_k[k],size=Nall_k[k])
                  E_kn_samp[k,0:Nall_k[k]] = E_kn[k,booti] 

       for k in range(K):
             for l in range(K):
                  u_kln[k,l,0:Nall_k[k]] = beta_k[l] * E_kn_samp[k,0:Nall_k[k]]

#------------------------------------------------------------------------
# Initialize MBAR
#------------------------------------------------------------------------

# Initialize MBAR with Newton-Raphson
#print ""
#print "Initializing MBAR:"
#print "--K = number of Temperatures"
#print "--L = number of Temperatures"
#print "--N = number of Energies per Temperature"

# Use Adaptive Method (Both Newton-Raphson and Self-Consistent, testing which is better)
       mbar = pymbar.MBAR(u_kln, Nall_k, method = 'adaptive', verbose=True, relative_tolerance=1e-12)

        #------------------------------------------------------------------------
        # Compute Expectations for E_kt and E2_kt as E_expect and E2_expect
        #------------------------------------------------------------------------
       print ""
       print "Computing Expectations for E..."
       (E_expect, dE_expect) = mbar.computeExpectations(u_kln)*(beta_k)**(-1)
       allE_expect[:,n] = E_expect[:]
        # expectations for the differences, which we need for numerical derivatives
       (DeltaE_expect, dDeltaE_expect) = mbar.computeExpectations(u_kln,output='differences')*(beta_k)**(-1)
       print "Computing Expectations for E^2..."
       (E2_expect,dE2_expect) = mbar.computeExpectations(u_kln*u_kln)*(beta_k)**(-2)
       allE2_expect[:,n] = E2_expect[:]

        #------------------------------------------------------------------------
       # Compute Cv for NVT simulations as <E^2> - <E>^2 / (RT^2)
        #------------------------------------------------------------------------
       print ""
       print "Computing Heat Capacity as ( <E^2> - <E>^2 ) / ( R*T^2 )..."

       for i in range(K):
             # heat capacity by fluctutations
             allCv_expect[i,0,n] = (E2_expect[i] - (E_expect[i]*E_expect[i])) / ( kB * Temp_k[i] * Temp_k[i])

               # heat capacity by T-differences
             im = i-1
             ip = i+1
             if (i==0):
                  im = 0
             if (i==K-1):
                  ip = K-1
             allCv_expect[i,1,n] = (DeltaE_expect[im,ip])/(Temp_k[ip]-Temp_k[im])

Cv_boot = numpy.zeros([K,2])
dCv_boot = numpy.zeros([K,2])

E_boot = numpy.zeros([K,nBoots], numpy.float64)
dE_boot = numpy.zeros([K,nBoots], numpy.float64)

E2_boot = numpy.zeros([K,nBoots], numpy.float64)
dE2_boot = numpy.zeros([K,nBoots], numpy.float64)

for k in range(K):
       for i in range(2):
             Cv_boot[k,i] = numpy.mean(allCv_expect[k,i,:])
             dCv_boot[k,i] = numpy.std(allCv_expect[k,i,:])
       E_boot[k] = numpy.mean(allE_expect[k,:])
       dE_boot[k] = numpy.std(allE_expect[k,:])
       E2_boot[k] = numpy.mean(allE2_expect[k,:])
       dE2_boot[k] = numpy.std(allE2_expect[k,:])

print "Temperature  dA         <E> +/- d<E>      <E^2> +/- d<E^2>    Cv +/- dCv (var)   Cv +/- dCv (dT)"    
print "------------------------------------------------------------------------------------------------------"
for k in range(K):
       print "%8.3f %8.3f %9.3f +/- %5.3f  %9.1f +/- %5.1f  %7.4f +/- %6.4f %7.4f +/- %6.4f" % (Temp_k[k],mbar.f_k[k],E_expect[k],dE_expect[k],E2_expect[k],dE2_expect[k],Cv_boot[k,0],dCv_boot[k,0],Cv_boot[k,1],dCv_boot[k,1])

numpy.savetxt('/home/edz3fz/Csurf_weak.txt', Cv_boot)
numpy.savetxt('/home/edz3fz/dCsurf_weak.txt', dCv_boot)

import matplotlib.pyplot as plt
plt.figure(1)
plt.plot(Temp_k,Cv_boot,'k', lw=1.5)
plt.plot(Temp_k, Cv_boot+dCv_boot, '0.5')
plt.plot(Temp_k, Cv_boot-dCv_boot, '0.5')
#plt.errorbar(Temp_k, Cv_expect, yerr=dCv_expect)
plt.xlabel('Temperature (K)')
plt.ylabel('C (kcal/mol/K)')
plt.xlim(300,450)
#plt.title('Heat Capacity from Go like model MC simulation of 1PBG.pdb')
plt.savefig('heatcap.png')
plt.show()


