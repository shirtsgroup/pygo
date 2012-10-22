# Example illustrating the application of MBAR to compute a 1D PMF from an umbrella sampling simulation.
#
# The data represents an umbrella sampling simulation for the chi torsion of a valine sidechain in lysozyme L99A with benzene bound in the cavity.
# 
# REFERENCE
# 
# D. L. Mobley, A. P. Graves, J. D. Chodera, A. C. McReynolds, B. K. Shoichet and K. A. Dill, "Predicting absolute ligand binding free energies to a simple model site," Journal of Molecular Biology 371(4):1118-1134 (2007).
# http://dx.doi.org/10.1016/j.jmb.2007.06.002

import numpy # numerical array library
#from math import *
import pymbar # multistate Bennett acceptance ratio
import timeseries # timeseries analysis
import pdb
from optparse import OptionParser


parser=OptionParser()
parser.add_option("-t", "--temp", nargs=1, default=300, type="int", dest="temp", help="temperature")

(options,args) = parser.parse_args()
temperature = options.temp
# Constants.
#kB = 1.381e-23 * 6.022e23 / 1000.0 # Boltzmann constant in kJ/mol/K
kB = 0.00831447/4.184 # Boltzmann constant in kcal/mol/K

#temperature = 300 # assume a single temperature -- can be overridden with data from center.dat 
# Parameters
K = 9 # number of umbrellas
N_max = 1001 # maximum number of snapshots/simulation
T_k = numpy.ones(K,float)*temperature # inital temperatures are all equal 
beta = 1.0 / (kB * temperature) # inverse temperature of simulations (in 1/(kJ/mol))
z_min = 17. # min for PMF
z_max = 24. # max for PMF
nbins = 36 # number of bins for 1D PMF

# Allocate storage for simulation data
N_k = numpy.ones([K], numpy.int32)*N_max # N_k[k] is the number of snapshots from umbrella simulation k
K_k = numpy.zeros([K], numpy.float64) # K_k[k] is the spring constant (in kJ/mol/deg**2) for umbrella simulation k
z0_k = numpy.zeros([K], numpy.float64) # chi0_k[k] is the spring center location (in deg) for umbrella simulation k
z_kn = numpy.zeros([K,N_max], numpy.float64) # chi_kn[k,n] is the torsion angle (in deg) for snapshot n from umbrella simulation k
u_kn = numpy.zeros([K,N_max], numpy.float64) # u_kn[k,n] is the reduced potential energy without umbrella restraints of snapshot n of umbrella simulation k
g_k = numpy.zeros([K],numpy.float32);

z0_k = numpy.array([18,19,20,21,22,23.5,25,26.5,28])
K_k = numpy.array([2.5,2.5,5.,2.5,2.5,2.5,2.5,2.5,2.5])

z18 = numpy.loadtxt('replicaexchange/simlog45/z_traj18_%i.txt' % (temperature))
z19 = numpy.loadtxt('replicaexchange/simlog40/z_traj19_%i.txt' % (temperature))
z20 = numpy.loadtxt('replicaexchange/simlog38/z_traj20_%i.txt' % (temperature))
z21 = numpy.loadtxt('replicaexchange/simlog39/z_traj21_%i.txt' % (temperature))
z22 = numpy.loadtxt('replicaexchange/simlog41/z_traj22_%i.txt' % (temperature))
z23 = numpy.loadtxt('replicaexchange/simlog46/z_traj23_%i.txt' % (temperature))
z25 = numpy.loadtxt('replicaexchange/simlog47/z_traj25_%i.txt' % (temperature))
z26 = numpy.loadtxt('replicaexchange/simlog48/z_traj26_%i.txt' % (temperature))
z28 = numpy.loadtxt('replicaexchange/simlog49/z_traj28_%i.txt' % (temperature))
u18 = numpy.loadtxt('replicaexchange/simlog45/energy%i.txt' % (temperature))
u19 = numpy.loadtxt('replicaexchange/simlog40/energy%i.txt' % (temperature))
u20 = numpy.loadtxt('replicaexchange/simlog38/energy%i.txt' % (temperature))
u21 = numpy.loadtxt('replicaexchange/simlog39/energy%i.txt' % (temperature))
u22 = numpy.loadtxt('replicaexchange/simlog41/energy%i.txt' % (temperature))
u23 = numpy.loadtxt('replicaexchange/simlog46/energy%i.txt' % (temperature))
u25 = numpy.loadtxt('replicaexchange/simlog47/energy%i.txt' % (temperature))
u26 = numpy.loadtxt('replicaexchange/simlog48/energy%i.txt' % (temperature))
u28 = numpy.loadtxt('replicaexchange/simlog49/energy%i.txt' % (temperature))

z_kn[0,:] = z18[-N_max::]
z_kn[1,:] = z19[-N_max::]
z_kn[2,:] = z20[-N_max::]
z_kn[3,:] = z21[-N_max::]
z_kn[4,:] = z22[-N_max::]
z_kn[5,:] = z23[-N_max::]
z_kn[6,:] = z25[-N_max::]
z_kn[7,:] = z26[-N_max::]
z_kn[8,:] = z28[-N_max::]

u_kn[0,:] = u18[-N_max::]
u_kn[1,:] = u19[-N_max::]
u_kn[2,:] = u20[-N_max::]
u_kn[3,:] = u21[-N_max::]
u_kn[4,:] = u22[-N_max::]
u_kn[5,:] = u23[-N_max::]
u_kn[6,:] = u25[-N_max::]
u_kn[7,:] = u26[-N_max::]
u_kn[8,:] = u28[-N_max::]
z_min = numpy.min(z_kn) # min for PMF
z_max = numpy.max(z_kn) # max for PMF
#check if energy.txt includes surface energy, i don't think it does

for i in range(K):
	u_kn[i,:] -= K_k[i]*(z_kn[i,:]-z0_k[i])**2



## Read in umbrella spring constants and centers.
#infile = open('data/centers.dat', 'r')
#lines = infile.readlines()
#infile.close()
#for k in range(K):
    ## Parse line k.
    #line = lines[k]
    #tokens = line.split()
    #chi0_k[k] = float(tokens[0]) # spring center locatiomn (in deg)
    #K_k[k] = float(tokens[1]) * (numpy.pi/180)**2 # spring constant (read in kJ/mol/rad**2, converted to kJ/mol/deg**2)    
    #if len(tokens) > 2:
        #T_k[k] = float(tokens[2])  # temperature the kth simulation was run at.

#beta_k = 1.0/(kB*T_k)   # beta factor for the different temperatures
#DifferentTemperatures = True
#if (min(T_k) == max(T_k)):
    #DifferentTemperatures = False            # if all the temperatures are the same, then we don't have to read in energies.
## Read the simulation data
#for k in range(K):
    ## Read torsion angle data.
    #filename = 'data/prod%d_dihed.xvg' % k
    #print "Reading %s..." % filename
    #infile = open(filename, 'r')
    #lines = infile.readlines()
    #infile.close()
    ## Parse data.
    #n = 0
    #for line in lines:
        #if line[0] != '#' and line[0] != '@':
            #tokens = line.split()
            #chi = float(tokens[1]) # torsion angle
            ## wrap chi_kn to be within [-180,+180)
            #while(chi < -180.0):
                #chi += 360.0
            #while(chi >= +180.0):
                #chi -= 360.0
            #chi_kn[k,n] = chi
            
            #n += 1
    #N_k[k] = n

    #if (DifferentTemperatures):  # if different temperatures are specified the metadata file, 
                                 ## then we need the energies to compute the PMF
        ## Read energies
        #filename = 'data/prod%d_energies.xvg' % k
        #print "Reading %s..." % filename
        #infile = open(filename, 'r')
        #lines = infile.readlines()
        #infile.close()
        ## Parse data.
        #n = 0
        #for line in lines:
            #if line[0] != '#' and line[0] != '@':
                #tokens = line.split()            
                #u_kn[k,n] = beta_k[k] * (float(tokens[2]) - float(tokens[1])) # reduced potential energy without umbrella restraint
                #n += 1

    ## Compute correlation times for potential energy and chi
    ## timeseries.  If the temperatures differ, use energies to determine samples; otherwise, use the cosine of chi
            
    #if (DifferentTemperatures):        
        #g_k[k] = timeseries.statisticalInefficiency(u_kn[k,:], u_kn[k,0:N_k[k]])
        #print "Correlation time for set %5d is %10.3f" % (k,g_k[k])
        #indices = timeseries.subsampleCorrelatedData(u_kn[k,0:N_k[k]])
    #else:
        #chi_radians = chi_kn[k,0:N_k[k]]/(180.0/numpy.pi)
        #g_cos = timeseries.statisticalInefficiency(numpy.cos(chi_radians))
        #g_sin = timeseries.statisticalInefficiency(numpy.sin(chi_radians))
        #print "g_cos = %.1f | g_sin = %.1f" % (g_cos, g_sin)
        #g_k[k] = max(g_cos, g_sin)
        #print "Correlation time for set %5d is %10.3f" % (k,g_k[k])
        #indices = timeseries.subsampleCorrelatedData(chi_radians, g=g_k[k]) 
    ## Subsample data.
    #N_k[k] = len(indices)
    #u_kn[k,0:N_k[k]] = u_kn[k,indices]
    #chi_kn[k,0:N_k[k]] = chi_kn[k,indices]

#N_max = numpy.max(N_k) # shorten the array size
u_kln = numpy.zeros([K,K,N_max], numpy.float64) # u_kln[k,l,n] is the reduced potential energy of snapshot n from umbrella simulation k evaluated at umbrella l

# Set zero of u_kn -- this is arbitrary.
u_kn -= u_kn.min()

# Construct torsion bins
print "Binning data..."
delta = (z_max - z_min) / float(nbins)
# compute bin centers
bin_center_i = numpy.zeros([nbins], numpy.float64)
for i in range(nbins):
    bin_center_i[i] = z_min + delta/2 + delta * i
# Bin data
bin_kn = numpy.zeros([K,N_max], numpy.int32)
for k in range(K):
    for n in range(N_k[k]):
        # Compute bin assignment.
        bin_kn[k,n] = int((z_kn[k,n] - z_min) / delta)
# Evaluate reduced energies in all umbrellas
print "Evaluating reduced potential energies..."
for k in range(K):
    for n in range(N_k[k]):
        # Compute minimum-image torsion deviation from umbrella center l
        dchi = z_kn[k,n] - z0_k
        # Compute energy of snapshot n from simulation k in umbrella potential l
        u_kln[k,:,n] = u_kn[k,n] + (K_k) * dchi**2
# Initialize MBAR.
print "Running MBAR..."
mbar = pymbar.MBAR(u_kln, N_k, verbose = True, method = 'adaptive')

# Compute PMF in unbiased potential (in units of kT).
(f_i, df_i) = mbar.computePMF(u_kn, bin_kn, nbins)

# Write out PMF
print "PMF (in units of kT)"
print "%8s %8s %8s" % ('bin', 'f', 'df')
for i in range(nbins):
    print "%8.1f %8.3f %8.3f" % (bin_center_i[i], f_i[i], df_i[i])

import matplotlib.pyplot as plt
plt.errorbar(bin_center_i,f_i,df_i)
plt.xlabel('z distance from surface (A)')
plt.ylabel('PMF (kcal/mol)')
plt.show()
