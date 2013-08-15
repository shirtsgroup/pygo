#!/usr/bin/python2.4

import numpy
from math import *
import pymbar
import timeseries
import commands
import os
import pdb
import matplotlib.pyplot as plt
import optparse
import wham

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


