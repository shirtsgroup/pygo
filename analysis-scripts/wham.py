#!/usr/bin/python2.4

import numpy
from math import *
import pymbar # for MBAR analysis
import timeseries # for timeseries analysis
import commands
import os
import os.path

#===================================================================================================
# CONSTANTS
#===================================================================================================

kB = 1.3806503 * 6.0221415 / 4184.0 # Boltzmann constant in kcal/mol/K

#===================================================================================================
# SUBROUTINES
#===================================================================================================

def logSum(log_terms):
   """Compute the log of a sum of terms whose logarithms are provided.

   REQUIRED ARGUMENTS  
      log_terms is the array (possibly multidimensional) containing the logs of the terms to be summed.

   RETURN VALUES
      log_sum is the log of the sum of the terms.

   """

   # compute the maximum argument
   max_log_term = log_terms.max()

   # compute the reduced terms
   terms = numpy.exp(log_terms - max_log_term)

   # compute the log sum
   log_sum = log( terms.sum() ) + max_log_term

   # return the log sum
   return log_sum

def histogram_wham(beta_k, U_kn, N_k, nbins = 100, bin_width = None, maximum_iterations = 5000, relative_tolerance = 1.0e-8, initial_f_k = None):
   """Construct an initial guess of the f_k by histogram reweighting (specifically, WHAM [2]).

   ARGUMENTS
     beta_k (numpy K array) - inverse temperatures (in units of 1/energy)
     U_kn (numpy K x N_max array) - potential energies (in units of energy)
     N_k (numpy K array of numpy.int32) - number of samples per states, N_max = N_k.max()

   OPTIONAL ARGUMENTS
     nbins (int) - if specified, sets the number of bins to use (default: 100)
     bin_width (float) - if specified, sets the bin width (overrides nbins) (defulat: None)
     maximum_iterations (int) - maximum number of iterations to use
     relative_tolerance (floeat) - relative convergence tolerance (default: 1.0e-8)

   RETURNS
     f_k (numpy K array) - guess at initial state dimensionless free energies

   REFERENCE
     [2] Kumar S, Bouzida D, Swensen RH, Kollman PA, and Rosenberg JM. The weighted histogram analysis method
     for free-energy calculations on biomolecules. I. The Method. J. Comput Chem. 13:1011, 1992.
   
   """

   # Get sizes
   K = N_k.size
   N_max = N_k.max()

   # Create a list of indices of all configurations in kn-indexing.
   mask_kn = numpy.zeros([K,N_max], dtype=numpy.bool)
   for k in range(0,K):
      mask_kn[k,0:N_k[k]] = True
   # Create a list from this mask.
   sample_indices = numpy.where(mask_kn)

   # Construct histogram bins
   M = nbins # number of energy bins
   SMALL = 1.0e-6
   U_min = U_kn[sample_indices].min()
   U_max = U_kn[sample_indices].max()
   U_max += (U_max - U_min) * SMALL # increment by a bit
   delta_U = (U_max - U_min) / float(M)
   if (bin_width != None):
      delta_U = bin_width
      M = int(numpy.ceil((U_max - U_min) / bin_width))
      print "Using %d bins to achieve bin width of %f" % (M, delta_U)
   else:
      print "Bin width is %f energy units to achieve nbins = %d" % (delta_U, M)
   U_m = U_min + delta_U * (0.5 + numpy.arange(0,M,dtype = numpy.float64))
   H_m = numpy.zeros([M], numpy.float64)
   # assign snapshots to energy bins
   bin_kn = numpy.zeros([K,N_max], numpy.int32) # bin_kn[k,n] is bin index of U_kn[k,n]
   bin_kn[sample_indices] = numpy.array(((U_kn[sample_indices] - U_min) / delta_U), numpy.int32)
   H_m = numpy.bincount(bin_kn[sample_indices])
   
   # compute logs of various quantities
   LOG_ZERO = -1000.0 # replacement for log(0)
   log_H_m = numpy.ones([M], numpy.float64) * LOG_ZERO
   for m in range(M):
      if (H_m[m] > 0):
         log_H_m[m] = log(H_m[m])
   log_N_k = numpy.ones([K], numpy.float64) * LOG_ZERO
   for k in range(K):
      if (N_k[k] > 0):
         log_N_k[k] = log(N_k[k])

   # initialize free energies
   f_k = numpy.zeros([K], numpy.float64)
   if (initial_f_k != None):
      f_k = initial_f_k.copy()

   # iterate
   f_k_new = numpy.zeros([K], numpy.float64)
   max_delta = 0.0
   for iteration in range(maximum_iterations):
      # print "iteration %d" % iteration
      
      # Form auxiliary matrices, used in both self-consistent iteration and Newtom-Raphson.
      # log space
      log_w_mi = numpy.zeros([M,K], numpy.float64)
      for m in range(M):
         # denominator = \sum_k N_k exp[f_k - \beta_k U_m]
         exp_arg_k = log_N_k[:] + f_k[:] - beta_k[:]*U_m[m]
         log_w_mi[m,:] = exp_arg_k[:] - logSum(exp_arg_k[:])
      # real space
      w_mi = numpy.zeros([M,K], numpy.float64)
      for m in range(M):
         exp_arg_k = f_k[:] - beta_k[:]*U_m[m]
         max_arg = exp_arg_k.max()
         numerators_i = N_k[:] * numpy.exp(exp_arg_k[:] - max_arg)
         w_mi[m,:] = numerators_i[:] / numerators_i.sum()
         
      # Compute new estimates of log weights using self-consistent iteration.
      for i in range(K):
         # Compute new estimate of log weights.
         f_k_new[i] = f_k[i] + log(N_k[i]) - logSum(log_H_m[:] + log_w_mi[:,i])
      
      # shift to ensure f_k_new[0] = 0.0
      f_k_new -= f_k_new[0]
      # store difference
      Delta_f_k = f_k_new - f_k
      # update f_k
      f_k = f_k_new.copy()

      # If all f_k are zero, terminate
      if numpy.all(f_k == 0.0):
        break
      
      # Terminate when max((f - fold) / f) < relative_tolerance for all nonzero f.
      max_delta = (numpy.abs(Delta_f_k) / (numpy.abs(f_k_new)).max()).max()
      print "iteration %8d relative max_delta = %8e" % (iteration, max_delta)
      if numpy.isnan(max_delta) or (max_delta < relative_tolerance):
         break

   # if unconverged
   if (iteration == maximum_iterations):
      print "Did not converge in %d iterations (final relative tolerance %e)" % (maximum_iterations, max_delta)
      print "f_k = "
      print f_k
      return f_k
   
   # summary
   print "Converged to relative tolerance of %e (convergence tolerance %e) in %d iterations" % (max_delta, relative_tolerance, iteration)
   print "f_k = "
   print f_k

   # return the estimate of the dimensionless free energies
   return f_k


