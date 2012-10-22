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
dT = 1.25              # Temperature increment for calculating Cv(T)

###########################################################
# For Cv vs T    _____
#               /     \
#  ____________/       \____________
#
############################################################
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
	files.append('results/1PGB/surface/simlog6/energy'+str(int(T[i]))+'.txt')


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
   g[k] = timeseries.statisticalInefficiency(E_from_file[k])#,suppress_warning=True)
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
print "Computing Expectations for E^2..."
(E2_expect,dE2_expect) = mbar.computeExpectations(u_kln*u_kln)*(beta_k)**(-2)

#------------------------------------------------------------------------
# Compute Cv for NVT simulations as <E^2> - <E>^2 / (RT^2)
#------------------------------------------------------------------------
print ""
print "Computing Heat Capacity as ( <E^2> - <E>^2 ) / ( R*T^2 )..."

Cv_expect = numpy.zeros([K], numpy.float64)
dCv_expect = numpy.zeros([K], numpy.float64)

for i in range(K):
       Cv_expect[i] = (E2_expect[i] - (E_expect[i]*E_expect[i])) / ( kB * Temp_k[i] * Temp_k[i])
       dCv_expect[i] = 2*dE_expect[i]**2 / (kB *Temp_k[i]*Temp_k[i])   # from propagation of error

print "Temperature  dA         <E> +/- d<E>       <E^2> +/- d<E^2>       Cv +/- dCv"     
print "-------------------------------------------------------------------------------"
for k in range(K):
       print "%8.3f %8.3f %9.3f +/- %5.3f  %9.1f +/- %5.1f   %7.4f +/- %6.4f" % (Temp_k[k],mbar.f_k[k],E_expect[k],dE_expect[k],E2_expect[k],dE2_expect[k],Cv_expect[k], dCv_expect[k])


import matplotlib.pyplot as plt
plt.figure(1)
#plt.plot(Temp_k,E_expect)
#plt.errorbar(Temp_k, E_expect, yerr=dE_expect)
plt.plot(Temp_k,E2_expect)
plt.errorbar(Temp_k, E2_expect, yerr=dE2_expect)
#plt.plot(Temp_k,Cv_expect)
#plt.errorbar(Temp_k, Cv_expect, yerr=dCc_expect)
plt.xlabel('Temperature (K)')
#plt.ylabel('Ave E (kcal/mol/K)')
#plt.title('Ave E from Go like model MC simulation of 1PBG.pdb')
plt.ylabel('Ave E^2 (kcal/mol/K)')
plt.title('Ave E^2 from Go like model MC simulation of 1PBG.pdb')
#plt.ylabel('Heat Capacity (kcal/mol/K)')
#plt.title('Heat Capacity from Go like model MC simulation of 1PBG.pdb')

#plt.savefig('heatcap.png')
#plt.show()

# single point heat capacities

K = len(T)
vfile = numpy.zeros(K)
dvfile = numpy.zeros(K)
for k in range(K):
       #vfile[k] = numpy.average(E_from_file[k,0:N_k[k]])
       #dvfile[k] = numpy.std(E_from_file[k,0:N_k[k]])/numpy.sqrt(N_k[k])
       vfile[k] = numpy.average((E_from_file[k,0:N_k[k]])**2)
       dvfile[k] = numpy.std((E_from_file[k,0:N_k[k]])**2)/numpy.sqrt(N_k[k])
       #vfile[k] = (numpy.average(E_from_file[k,0:N_k[k]]**2)-numpy.average(E_from_file[k,0:N_k[k])**2)/(kB*T**2)
       #defile[k] = ?
pdb.set_trace()
plt.plot(T,vfile, 'bo')
plt.errorbar(T, vfile, yerr=dvfile)
#plt.savefig('aveE_singlepoint.png')
plt.savefig('aveE2_singlepoint.png')
#plt.savefig('heatcap_singlepoint.png')
plt.show()
