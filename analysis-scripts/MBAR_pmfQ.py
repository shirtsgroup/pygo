#!/usr/bin/python2.4

import numpy
from math import *
import pymbar
import timeseries
import commands
import os
#import pdb
#import matplotlib.pyplot as plt
import optparse
import wham
import cPickle

def parse_args():
	parser = optparse.OptionParser(description='Calculates the PMF(Q)')
#	parser.add_option('-t','--temp', dest= 'temp', nargs = 2, type = 'float', help = 'desired temperatures')
	parser.add_option('--tfile', dest='tfile', default='T.txt', help = 'simulation temperature file')
	parser.add_option('--direc', dest='direc', help='directory of simulation data')
	parser.add_option("-n", "--N_max", default=100000, type="int",dest="N_max", help="number of data points to read in (default: 100k)")
	parser.add_option("-s", "--skip", default=1, type="int",dest="skip", help="skip every n data points")
	(options,args) = parser.parse_args()
	return options

def read_data(args,K,T):
    print "Reading data..."
    U_kn = numpy.empty([K,args.N_max/args.skip], numpy.float64)
    Q_kn = numpy.empty([K,args.N_max/args.skip], numpy.float64)
    for i, t in enumerate(T):
        ufile = '%s/energy%i.npy' %(args.direc, t)
        data = numpy.load(ufile)[-args.N_max::]
        U_kn[i,:] = data[::args.skip]
        Qfile = '%s/fractionnative%i.npy' %(args.direc, t)
        data = numpy.load(Qfile)[-args.N_max::]
        Q_kn[i,:] = data[::args.skip]
    N_max = args.N_max/args.skip
    return U_kn, Q_kn, N_max

def subsample(U_kn,Q_kn,K,N_max):
    assume_uncorrelated = False
    if assume_uncorrelated:
        print 'Assuming data is uncorrelated'
        N_k = numpy.zeros(K, numpy.int32)
        N_k[:] = N_max
    else:	
        print 'Subsampling the data...'
        N_k = numpy.zeros(K,numpy.int32)
        g = numpy.zeros(K,numpy.float64)
        for k in range(K):  # subsample the energies
            g[k] = timeseries.statisticalInefficiency(Q_kn[k])#,suppress_warning=True)
            indices = numpy.array(timeseries.subsampleCorrelatedData(Q_kn[k],g=g[k])) # indices of uncorrelated samples
            N_k[k] = len(indices) # number of uncorrelated samplesadsf
            U_kn[k,0:N_k[k]] = U_kn[k,indices]
            Q_kn[k,0:N_k[k]] = Q_kn[k,indices]
    return U_kn, Q_kn, N_k

def get_bins(nbins,K,N_max,Q_kn):
    print 'Binning Q'
    Q_min = 0
    Q_max = 1
    dQ = (Q_max - Q_min) / float(nbins)
    bin_kn = numpy.zeros([K,N_max],numpy.int16)
    bin_counts = list()
    bin_centers = list()
    bin = 0
    for i in range(nbins):
        Q = Q_min + dQ * (i + 0.5)
        in_bin = (Q-dQ/2 <= Q_kn) & (Q_kn < Q+dQ/2)
        bin_count = in_bin.sum()
        if bin_count > 0:
            bin_centers.append(Q)
            bin_counts.append(bin_count)
            bin_kn[in_bin] = bin
            bin += 1
    return bin, bin_centers, bin_counts, bin_kn

def get_ukln(N_max,K,U_kn,N_k,beta_k):
    print 'Computing reduced potential energies...'
    u_kln = numpy.zeros([K,K,N_max], numpy.float64) # u_kln is reduced pot. ener. of segment n of temp k evaluated at temp l
    for k in range(K):
	       for l in range(K):
	             u_kln[k,l,0:N_k[k]] = beta_k[l] * U_kn[k,0:N_k[k]]
    return u_kln

def get_mbar(beta_k,U_kn,N_k,u_kln):
    print 'Initializing mbar...'
    f_k = wham.histogram_wham(beta_k, U_kn, N_k)
    mbar = pymbar.MBAR(u_kln, N_k, initial_f_k = f_k, verbose=True)
    return mbar 

def main():
    options = parse_args()
	
    # set constants
    kB = 0.00831447/4.184
    nbins = 25

    # get temperature states
    T = numpy.loadtxt(options.tfile)
    beta_k = 1 / (kB * T)
    print 'temperature states are\n', T
    K = len(T)
	
    # read in data
    U_kn, Q_kn, N_max = read_data(options,K,T)

    # test for statistical inefficiencies
    U_kn, Q_kn, N_k = subsample(U_kn, Q_kn, K, N_max)

	# binning data
    nbins, bin_centers, bin_counts, bin_kn = get_bins(nbins,K,N_max,Q_kn)
    print '%i bins were populated:' %nbins
    for i in range(nbins):
        print 'bin %5i (%6.1f) %12i conformations' % (i, bin_centers[i], bin_counts[i])

    # compute reduced potential energy of all snapshots at all temperatures
    u_kln = get_ukln(N_max,K,U_kn,N_k,beta_k)

    # initialize mbar
    mbar = get_mbar(beta_k,U_kn,N_k,u_kln)

    target_temperatures = numpy.arange(295.,360.,5.)
#    target_temperatures = [300, 325, 350]
    print 'Calculating the PMF at', target_temperatures

#	f_i = numpy.zeros((nbins,len(target_temperatures)))
#	df_i = numpy.zeros((nbins,len(target_temperatures)))
    for i,temp in enumerate(target_temperatures):
        target_beta = 1.0 / (kB * temp)
        u_kn = target_beta * U_kn
        f_i, d2f_i = mbar.computePMF_states(u_kn, bin_kn, nbins)
#       imin = f_i.argmin()
#       for j in range(nbins):
#			df_i[j,i] = sqrt(d2f_i[j,imin]) # uncertainty relative to lowest free energy
        pmf_file = '%s/pmf_%i.pkl' % (options.direc, temp)
        f = file(pmf_file, 'wb')
        print 'Saving target temperature, bin centers, f_i, df_i to %s' % pmf_file
        cPickle.dump(temp,f)
        cPickle.dump(bin_centers,f)
        cPickle.dump(f_i,f)
        cPickle.dump(d2f_i,f)
        f.close()
	
if __name__ == '__main__':
    main()

