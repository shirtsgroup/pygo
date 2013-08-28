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

def parse_args():
	parser = optparse.OptionParser(description='Calculates the PMF(Q) at various input temperatures')
	parser.add_option('-t','--temp', dest= 'temp', nargs = 2, type = 'float', help = 'desired temperatures')
	parser.add_option('--tfile', dest='tfile', default='T.txt', help = 'simulation temperature file')
	parser.add_option('--direc', dest='direc', help='directory of simulation data')
	parser.add_option("-n", "--N_max", default=100000, type="int",dest="N_max", help="number of data points to read in (default: 100k)")
	parser.add_option("-s", "--skip", default=1, type="int",dest="skip", help="skip every n data points")
	(options,args) = parser.parse_args()
	return options

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

if __name__ == "__main__":
	options = parse_args()
	target_temperatures = list(options.temp)
	
	direc = options.direc
	tfile = options.tfile
	N_max = options.N_max
	skip = options.skip

	kB = 0.00831447/4.184
	nbins = 30

	T = numpy.loadtxt(tfile)
	
	print T
	
	K = len(T)
	U_kn = numpy.empty([K,N_max/skip], numpy.float64)
	Q_kn = numpy.empty([K,N_max/skip], numpy.float64)
	
	print "Reading data..."
	for i, t in enumerate(T):
		ufile = '%s/energy%i.npy' %(direc, t)
		data = numpy.load(ufile)[-N_max::]
		U_kn[i,:] = data[::skip]
		Qfile = '%s/fractionnative%i.npy' %(direc, t)
		data = numpy.load(Qfile)[-N_max::]
		Q_kn[i,:] = data[::skip]
	
	N_max /= skip
	
	N_k = numpy.zeros(K,numpy.int32)
	g = numpy.zeros(K,numpy.float64)
	for k in range(K):  # subsample the energies
	   g[k] = timeseries.statisticalInefficiency(Q_kn[k])#,suppress_warning=True)
	   indices = numpy.array(timeseries.subsampleCorrelatedData(Q_kn[k],g=g[k])) # indices of uncorrelated samples
	   N_k[k] = len(indices) # number of uncorrelated samplesadsf
	   U_kn[k,0:N_k[k]] = U_kn[k,indices]
	   Q_kn[k,0:N_k[k]] = Q_kn[k,indices]

	print 'Binning Q'
	Q_min = 0
	Q_max = 1
	dx = (Q_max - Q_min) / float(nbins)
	bin_kn = numpy.zeros([K,N_max],numpy.int16)
	bin_counts = list()
	bin_centers = list()
	bin = 0
	for i in range(nbins):
		Q = Q_min + dx * (i + 0.5)
		in_bin = (Q-dx/2 <= Q_kn) & (Q_kn < Q+dx/2)
		bin_count = in_bin.sum()
		if bin_count > 0:
			bin_centers.append(Q)
			bin_counts.append(bin_count)
			bin_kn[in_bin] = bin
			bin += 1
	
	beta_k = 1 / (kB * T)
	print "--Computing reduced energies..."	
	u_kln = numpy.zeros([K,K,N_max], numpy.float64) # u_kln is reduced pot. ener. of segment n of temp k evaluated at temp l
	for k in range(K):
	       for l in range(K):
	             u_kln[k,l,0:N_k[k]] = beta_k[l] * U_kn[k,0:N_k[k]]
	f_k = histogram_wham(beta_k, U_kn, N_k)
	mbar = pymbar.MBAR(u_kln, N_k, initial_f_k = f_k)
	f_i = numpy.zeros((nbins,len(target_temperatures)))
	df_i = numpy.zeros((nbins,len(target_temperatures)))
	for i,temp in enumerate(target_temperatures):
		target_beta = 1.0 / (kB * temp)
		u_kn = target_beta * U_kn
		f_i[:,i], d2f_i = mbar.computePMF_states(u_kn, bin_kn, nbins)
		imin = f_i.argmin()
		for j in range(nbins):
			df_i[j,i] = sqrt(d2f_i[j,imin]) # uncertainty relative to lowest free energy
	
	#plotting
	for i in range(len(target_temperatures)):
		plt.errorbar(bin_centers,f_i[:,i],df_i[:,i])
	plt.xlabel('Q')
	plt.ylabel('PMF (kcal/mol)')
	plt.show()
