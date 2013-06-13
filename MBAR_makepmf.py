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
import matplotlib.pyplot as plt


parser=OptionParser()
parser.add_option("-t", "--temp", nargs=1, default=200, type="int", dest="temp", help="temperature")
parser.add_option("--direc", dest="datafile",  help="directory of files")
parser.add_option("-n", default=8, type="int",dest="n", help="number of umbrellas")
parser.add_option("-z", "--zrange", nargs=2, default=[9,31.5], type="float", dest="zrange", help="range of z pinnings for umbrellas, assuming spacing of 2")
parser.add_option("-k", dest="k", type="float",default=1, help="harmonic constant")

(options,args) = parser.parse_args()
T = options.temp
# Constants.
#kB = 1.381e-23 * 6.022e23 / 1000.0 # Boltzmann constant in kJ/mol/K
kB = 0.00831447/4.184 # Boltzmann constant in kcal/mol/K

#temperature = 300 # assume a single temperature -- can be overridden with data from center.dat 
# Parameters
K = options.n # number of umbrellas
k = options.k
zrange = options.zrange
direc = options.datafile
N_max = 5001 # maximum number of snapshots/simulation
T_k = numpy.ones(K,float)*T # inital temperatures are all equal 
beta = 1.0 / (kB * T) # inverse temperature of simulations (in 1/(kJ/mol))
nbins = 25 # number of bins for 1D PMF
# Allocate storage for simulation data
N_k = numpy.ones([K], numpy.int32)*N_max # N_k[k] is the number of snapshots from umbrella simulation k
K_k = k*numpy.ones(K, numpy.float64) # K_k[k] is the spring constant (in kJ/mol/deg**2) for umbrella simulation k
z0_k = numpy.arange(zrange[0],zrange[1],1.5, numpy.float64) # chi0_k[k] is the spring center location (in deg) for umbrella simulation k
z_kn = numpy.zeros([K,N_max], numpy.float64) # chi_kn[k,n] is the torsion angle (in deg) for snapshot n from umbrella simulation k
u_kn = numpy.zeros([K,N_max], numpy.float64) # u_kn[k,n] is the reduced potential energy without umbrella restraints of snapshot n of umbrella simulation k

for index, i in enumerate(z0_k):
	zfile = '%s/%i/z_traj%i_%i.npy' %(direc, int(i), int(i), T)
	z_kn[index,:]=numpy.load(zfile)[-N_max::]
	plt.hist(z_kn[index,:])
	ufile = '%s/%i/energy%i.npy' %(direc, int(i), T)
	u_kn[index,:]=numpy.load(ufile)[-N_max::]

plt.show()	

z_min = numpy.min(z_kn) # min for PMF
z_max = numpy.max(z_kn) # max for PMF
#check if energy.txt includes surface energy, i don't think it does

for i in range(K):
	u_kn[i,:] -= K_k[i]*(z_kn[i,:]-z0_k[i])**2


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
pdb.set_trace()
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
f_i -= f_i[-1]
# Write out PMF
print "PMF (in units of kT)"
print "%8s %8s %8s" % ('bin', 'f', 'df')
for i in range(nbins):
    print "%8.1f %8.3f %8.3f" % (bin_center_i[i], f_i[i], df_i[i])

plt.errorbar(bin_center_i,f_i,df_i)
plt.xlabel('z distance from surface (A)')
plt.ylabel('PMF (kcal/mol)')
save = numpy.zeros((nbins,3))
save[:,0] = bin_center_i
save[:,1] = f_i
save[:,2] = df_i
numpy.save(direc+'/PMF%i' %T, save)
df = numpy.sqrt(df_i[numpy.argmin(f_i)]**2 + df_i[numpy.argmax(f_i)]**2)
print 'The well depth is %f +- %f kcal/mol' % (numpy.min(f_i), df)
plt.savefig(direc+'/PMF%i.png' %T)
plt.show()
