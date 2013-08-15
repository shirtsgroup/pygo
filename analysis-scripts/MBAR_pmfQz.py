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

def parse_args():
	parser = optparse.OptionParser(description='Calculates the PMF(Q) at various input temperatures')
#	parser.add_option('-t','--temp', dest= 'temp', nargs = 2, type = 'float', help = 'desired temperatures')
	parser.add_option('--tfile', dest='tfile', default='T.txt', help = 'simulation temperature file')
	parser.add_option('--direc', dest='direc', help='directory of simulation data')
	parser.add_option("-n", "--N_max", default=100000, type="int",dest="N_max", help="number of data points to read in (default: 100k)")
	parser.add_option("-s", "--skip", default=1, type="int",dest="skip", help="skip every n data points")
	parser.add_option('--f_file', dest='f_file', default='', help='free energy filename, if it exists')
	(options,args) = parser.parse_args()
	return options

def main():
	options = parse_args()
#	target_temperatures = list(options.temp)
	
	direc = options.direc
	tfile = options.tfile
	N_max = options.N_max
	skip = options.skip
	f_file = options.f_file

	kB = 0.00831447/4.184
	nbins_per = 30
	spring_constant = 1.5 #check this 

	T = numpy.loadtxt(tfile)
	beta_k = 1 / (kB * T)
	print 'temperature states are', T
	Z = numpy.arange(9,31.5,1.5)
	print 'distance states are', z	

	K = len(T)*len(Z)

	U_kn = numpy.empty([K,N_max/skip], numpy.float64)
	Q_kn = numpy.empty([K,N_max/skip], numpy.float64)
	z_kn = numpy.empty([K,N_max/skip], numpy.float64)
	
	print "Reading data..."
	i = 0
	for z in Z: 
		for t in T:
			ufile = '%s/%i/energy%i.npy' %(direc, z, t)
			data = numpy.load(ufile)[-N_max::]
			U_kn[i,:] = data[::skip]
			Qfile = '%s/%i/fractionnative%i.npy' %(direc, z, t)
			data = numpy.load(Qfile)[-N_max::]
			Q_kn[i,:] = data[::skip]
			zfile = '%s/%i/z_traj%i_%i.npy' %(direc, z, z, t)
			data = numpy.load(zfile)[-N_max::]
			z_kn[i,:] = data[::skip]
			U_kn[i,:] -= spring_constant*(z_kn[i,:] - z)**2
			i += 1
	N_max /= skip
	
	# test for statistical inefficiencies
	assume_uncorrelated = True
	if assume_uncorrelated:
		'Assuming data is uncorrelated...'
		N_k = numpy.zeros(K, numpy.int32)
		N_k[:] = N_max
	else:	
		'Subsampling the data...'
		N_k = numpy.zeros(K,numpy.int32)
		g = numpy.zeros(K,numpy.float64)
		for k in range(K):  # subsample the energies
		   g[k] = timeseries.statisticalInefficiency(Q_kn[k])#,suppress_warning=True)
		   indices = numpy.array(timeseries.subsampleCorrelatedData(Q_kn[k],g=g[k])) # indices of uncorrelated samples
		   N_k[k] = len(indices) # number of uncorrelated samplesadsf
		   U_kn[k,0:N_k[k]] = U_kn[k,indices]
		   Q_kn[k,0:N_k[k]] = Q_kn[k,indices]
		   z_kn[k,0:N_k[k]] = z_kn[k,indices]

	# generate a list of indices of all configurations in kn-indicing
	mask_kn = numpy.zeros([K,N_max], dtype=numpy.bool)
	for k in range(0,K):
		mask_kn[k,0:N_k[k]] = True
	indices = nump.where(mask_kn)
	
	# compute reduced potential energy of all snapshots at all temperatures and distances
	print 'Computing reduced potential energies...'
	u_kln = numpy.zeros([K,K,N_max], numpy.float32)
	for k in range(K):
		for l in range(K):
			z_index = l/len(T) # z is outer dimension
			T_index = l%len(T) # T is inner dimension
			dz = z_kn[k,0:N_k[k]] - Z[z_index]
			u_kln[k,l,0:N_k[k]] = beta_k[T_index] * U_kn[k,0:N_k[k]] + spring_constant*(dz)**2
	temp_file = 'u_kln.npy'
	print 'Saving u_kln matrix to %s...' % temp_file
	numpy.save(temp_file, u_kln)

	# bin data for PMF calculation
	print 'Binning data...'
	Q_min = 0
	Q_max = 1
	dQ = (Q_max - Q_min) / float(nbins_per)
	z_min = numpy.min(z_kn[indices])
	z_max = numpy.max(z_kn[indices])
	dz = (z_max - z_in) / float(nbins_per)
	bin_kn = numpy.zeros([K,N_max],numpy.int16)
	bin_counts = []
	bin_centers = []
	nbins = 0
	for i in range(nbins_per):
		for j in range(nbins_per):
			z = z_min + dz * (i + 0.5)
			Q = Q_min + dQ * (j + 0.5)
			in_bin = (Q-dQ/2 <= Q_kn[indices]) and (Q_kn[indices] < Q+dQ/2) and (z-dz/2 <= z_kn[indices]) and (z_kn[indices] < z+dz/2)
			bin_count = in_bin.sum()
			if bin_count > 0:
				bin_centers.append((z,Q))
				bin_counts.append(bin_count)
				bin_kn[in_bin] = bin
				nbins += 1
	print '%i bins were populated:' %nbins
	for i in range(nbins):
		print 'bin %5i (%6.1f, %6.1f) %12i conformations' % (i, bin_centers[i][0], bin_centers[i][1], bin_counts[i])

	# read free energies if they exist
	# otherwise use WHAM to quickly compute an initial guess of dimensionless free energies f_k
	# then initialize MBAR
	if f_file:
		f_k = numpy.load(f_file)
		print 'Reading free energies from %s...' % f_file
		mbar = pymbar.MBAR(u_kln, N_k, maximum_iterations = 0, verbose = True, initial_f_k = f_k)
	else:
		print 'Using WHAM to generate historgram-based initial guess of dimensionless free energies f_k...'
		f_k = histogram_wham(beta_k, U_kn, N_k)
		print 'Initializing MBAR...'
		mbar = pymbar.MBAR(u_kln, N_k, initial_f_k = f_k)
		f_file = 'f_k.pkl'
		print 'Writing free energies to %s' % f_file
		numpy.save(f_file, mbar.f_k)
	
	mbar_file = 'mbar.pkl'
	with open(mbar_file,'wb'):
		print 'Saving mbar object to %s' mbar_file
		pickle.dump(mbar_file, mbar)
	

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

