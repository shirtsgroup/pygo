import numpy
from math import *
import pymbar # for MBAR analysis
import timeseries # for timeseries analysis
import commands
import os
import os.path
import pdb
from optparse import OptionParser
import matplotlib.pyplot as plt


#===================================================================================================
# CONSTANTS
#===================================================================================================

kB = 1.3806503 * 6.0221415 / 4184.0 # Boltzmann constant in kcal/mol/K
N_max = 10001 # maximum number of snapshots/state
#===================================================================================================
# PARAMETERS
#===================================================================================================

parser=OptionParser()
parser.add_option("-t", "--temp", nargs=1, default=300, type="int", dest="temp", help="target temperature")
parser.add_option("--tfile", dest="tfile", default="T.txt", help="temperature file T.txt")
parser.add_option("--direc", dest="datafile",  help="directory of files")
parser.add_option("-k", dest="k", type="float",default=1, help="harmonic constant")
parser.add_option("-z", dest="z", type="float",default=30, help="umbrella")

(options,args) = parser.parse_args()
target_temperature = options.temp
nbins = 40
tfile = options.tfile
k = options.k
direc = options.datafile
z0 = options.z

#data_directory = 'data/' # directory containing the parallel tempering data
#temperature_list_filename = os.path.join(data_directory, 'temperatures') # file containing temperatures in K
#potential_energies_filename = os.path.join(data_directory, 'energies', 'potential-energies') # file containing total energies (in kcal/mol) for each temperature and snapshot
#trajectory_segment_length = 20 # number of snapshots in each contiguous trajectory segment
#niterations = 500 # number of iterations to use
#target_temperature = 302 # target temperature for 2D PMF (in K)
#nbins_per_torsion = 10 # number of bins per torsion dimension

#===================================================================================================
# SUBROUTINES
#===================================================================================================
def logSum(log_terms):
   max_log_term = log_terms.max()
   terms = numpy.exp(log_terms - max_log_term)
   log_sum = log( terms.sum() ) + max_log_term
   return log_sum

def histogram_wham(beta_k, U_kn, N_k, nbins = 100, bin_width = None, maximum_iterations = 5000, relative_tolerance = 1.0e-8, initial_f_k = None):
   K = N_k.size
   N_max = N_k.max()
   mask_kn = numpy.zeros([K,N_max], dtype=numpy.bool)
   for k in range(0,K):
      mask_kn[k,0:N_k[k]] = True
   sample_indices = numpy.where(mask_kn)
   M = nbins # number of energy bins
   SMALL = 1.0e-6
   U_min = U_kn[sample_indices].min()
   U_max = U_kn[sample_indices].max()
   U_max += (U_max - U_min) * SMALL # increment by a bit
   delta_U = (U_max - U_min) / float(M)
   if (bin_width != None):
      delta_U = bin_width
      M = int(numpy.ceil((U_max - U_min) / bin_width))
      print "Using %d bins to achieve bin width of %f" % (M, delta_U)
   else:
      print "Bin width is %f energy units to achieve nbins = %d" % (delta_U, M)
   U_m = U_min + delta_U * (0.5 + numpy.arange(0,M,dtype = numpy.float64))
   H_m = numpy.zeros([M], numpy.float64)
   bin_kn = numpy.zeros([K,N_max], numpy.int32) # bin_kn[k,n] is bin index of U_kn[k,n]
   bin_kn[sample_indices] = numpy.array(((U_kn[sample_indices] - U_min) / delta_U), numpy.int32)
   H_m = numpy.bincount(bin_kn[sample_indices])
   LOG_ZERO = -1000.0 # replacement for log(0)
   log_H_m = numpy.ones([M], numpy.float64) * LOG_ZERO
   for m in range(M):
      if (H_m[m] > 0):
         log_H_m[m] = log(H_m[m])
   log_N_k = numpy.ones([K], numpy.float64) * LOG_ZERO
   for k in range(K):
      if (N_k[k] > 0):
         log_N_k[k] = log(N_k[k])
   f_k = numpy.zeros([K], numpy.float64)
   if (initial_f_k != None):
      f_k = initial_f_k.copy()
   f_k_new = numpy.zeros([K], numpy.float64)
   max_delta = 0.0
   for iteration in range(maximum_iterations):
      log_w_mi = numpy.zeros([M,K], numpy.float64)
      for m in range(M):
         exp_arg_k = log_N_k[:] + f_k[:] - beta_k[:]*U_m[m]
         log_w_mi[m,:] = exp_arg_k[:] - logSum(exp_arg_k[:])
      w_mi = numpy.zeros([M,K], numpy.float64)
      for m in range(M):
         exp_arg_k = f_k[:] - beta_k[:]*U_m[m]
         max_arg = exp_arg_k.max()
         numerators_i = N_k[:] * numpy.exp(exp_arg_k[:] - max_arg)
         w_mi[m,:] = numerators_i[:] / numerators_i.sum()
      for i in range(K):
         f_k_new[i] = f_k[i] + log(N_k[i]) - logSum(log_H_m[:] + log_w_mi[:,i])
      f_k_new -= f_k_new[0]
      Delta_f_k = f_k_new - f_k
      f_k = f_k_new.copy()
      if numpy.all(f_k == 0.0):
        break
      max_delta = (numpy.abs(Delta_f_k) / (numpy.abs(f_k_new)).max()).max()
      print "iteration %8d relative max_delta = %8e" % (iteration, max_delta)
      if numpy.isnan(max_delta) or (max_delta < relative_tolerance):
         break
   if (iteration == maximum_iterations):
      print "Did not converge in %d iterations (final relative tolerance %e)" % (maximum_iterations, max_delta)
      print "f_k = "
      print f_k
      return f_k
   print "Converged to relative tolerance of %e (convergence tolerance %e) in %d iterations" % (max_delta, relative_tolerance, iteration)
   print "f_k = "
   print f_k
   return f_k

def write_free_energies(filename, f_k):
   K = f_k.size
   outfile = open(filename, 'w')
   for k in range(K):
      outfile.write("%f " % f_k[k])
   outfile.write("\n")
   outfile.close()
   return
#===================================================================================================
# MAIN
#===================================================================================================

#===================================================================================================
# Read temperatures
#===================================================================================================

temperature_k = numpy.loadtxt(tfile) # temperature_k[k] is temperature of temperature index k in K
K = len(temperature_k)
beta_k = (kB * temperature_k)**(-1) # compute inverse temperatures

# Read list of temperatures.
#lines = read_file(temperature_list_filename)
# Construct list of temperatures
#temperatures = lines[0].split()
# Create numpy array of temperatures
#K = len(temperatures)
#temperature_k = numpy.zeros([K], numpy.float32) # temperature_k[k] is temperature of temperature index k in K
#for k in range(K):
#   temperature_k[k] = float(temperatures[k])
# Compute inverse temperatures
#beta_k = (kB * temperature_k)**(-1) 

# Define other constants
#T = trajectory_segment_length * niterations # total number of snapshots per temperature

#===================================================================================================
# Read potential eneriges
#===================================================================================================

z_kn = numpy.empty([K,N_max], numpy.float64) # z_kn[k,n] is the z dist from the surface (in A) for snapshot n from temperature simulation k
U_kn = numpy.empty([K,N_max], numpy.float64)
Q_kn = numpy.empty([K,N_max], numpy.float64)

print "Reading data..."
for i, t in enumerate(temperature_k):
	zfile = '%s/%i/z_traj%i_%i.npy' %(direc, z0, z0, t)
	z_kn[i,:] = numpy.load(zfile)[-N_max::]
	ufile = '%s/%i/energy%i.npy' %(direc, z0, t)
	U_kn[i,:] = numpy.load(ufile)[-N_max::]
	Qfile = '%s/%i/fractionnative%i.npy' %(direc, z0, t)
	Q_kn[i,:] = numpy.load(Qfile)[-N_max::]
for i in range(K):
	U_kn[i,:] -= k*(z_kn[i,:] - z0)**2


#===================================================================================================
# Compute reduced potential energy of all snapshots at all temperatures
#===================================================================================================

print "Computing reduced potential energies..."
u_kln = numpy.zeros([K,K,N_max], numpy.float32) # u_kln[k,l,n] is reduced potential energy of trajectory segment n of temperature k evaluated at temperature l
for k in range(K):
   for l in range(K):
      u_kln[k,l,:] = beta_k[l] * U_kn[k,:]
N_k = numpy.zeros([K], numpy.int32)
N_k[:] = N_max
#===================================================================================================
# Bin torsions into histogram bins for PMF calculation
#===================================================================================================

# Here, we bin the (phi,psi) samples into bins in a 2D histogram.
# We assign indices 0...(nbins-1) to the bins, even though the histograms are in two dimensions.
# All bins must have at least one sample in them.
# This strategy scales to an arbitrary number of dimensions.
print "Binning Q..."
# Determine torsion bin size (in degrees)
Q_min = 0.
Q_max = 1.
dx = (Q_max - Q_min) / float(nbins)
# Assign torsion bins
bin_kn = numpy.zeros([K,N_max], numpy.int16) # bin_kn[k,n] is the index of which histogram bin sample n from temperature index k belongs to
nbins2 = nbins
nbins = 0
bin_counts = list()
bin_centers = list() # bin_centers[i] is a (phi,psi) tuple that gives the center of bin i
for i in range(nbins2):
      Q = Q_min + dx * (i + 0.5)
      in_bin = (Q-dx/2 <= Q_kn) & (Q_kn < Q+dx/2)
      bin_count = in_bin.sum()
      if (bin_count > 0):
         bin_centers.append( Q )
         bin_counts.append( bin_count )
         bin_kn[in_bin] = nbins
         nbins += 1
print "Using WHAM to generate histogram-based initial guess of dimensionless free energies f_k..."
f_k = histogram_wham(beta_k, U_kn, N_k)
print "Initializing MBAR (will estimate free energy differences first time)..."
mbar = pymbar.MBAR(u_kln, N_k, verbose = True, initial_f_k = f_k)
#write_free_energies(free_energies_filename, mbar.f_k)

#===================================================================================================
# Compute PMF at the desired temperature.
#===================================================================================================

print "Computing potential of mean force..."

# Compute reduced potential energies at the temperaure of interest
target_beta = 1.0 / (kB * target_temperature)
u_kn = target_beta * U_kn
# Compute PMF at this temperature, returning dimensionless free energies and uncertainties.
# f_i[i] is the dimensionless free energy of bin i (in kT) at the temperature of interest
# d2f_ij[i,j] is an estimate of the covariance in the estimate of (f_i[i] - f_j[j])
(f_i, d2f_i) = mbar.computePMF_states(u_kn, bin_kn, nbins)

# Find index of bin with lowest free energy.
imin = f_i.argmin()

# Show free energy and uncertainty of each occupied bin relative to lowest free energy
print "PMF"
print ""
print "%8s %6s %8s %10s %10s" % ('bin', 'Q', 'N', 'f', 'df')
error = numpy.zeros(nbins)
for i in range(nbins):
   print '%8d %6.1f %8d %10.3f %10.3f' % (i, bin_centers[i], bin_counts[i], f_i[i], sqrt(d2f_i[i,imin]))
   error[i] = sqrt(d2f_i[i,imin])

plt.errorbar(bin_centers,f_i,error)
plt.xlabel('z distance from surface (A)')
plt.ylabel('PMF (kcal/mol)')
save = numpy.zeros((nbins,3))
save[:,0] = bin_centers
save[:,1] = f_i
save[:,2] = error
numpy.save(direc+'/PMF%i_Q%i' % ( z0, target_temperature), save)
#df = numpy.sqrt(df_i[numpy.argmin(f_i)]**2 + df_i[numpy.argmax(f_i)]**2)
#print 'The well depth is %f +- %f kcal/mol' % (numpy.min(f_i), df)
plt.savefig(direc+'/PMF%i_Q%i.png' % (z0,target_temperature))
plt.show()
