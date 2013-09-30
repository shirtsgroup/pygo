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
import MBAR_pmfQ
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

def get_2_state_bins(bin_centers,K,N_max,Q_kn):
    print 'Binning Q'
    dQ = .15
    bin_kn = numpy.zeros([K,N_max],numpy.int16)
    bin_counts = list()
    for i,Q in enumerate(bin_centers):
        in_bin = (Q-dQ/2 <= Q_kn) & (Q_kn < Q+dQ/2)
        bin_count = in_bin.sum()
        bin_counts.append(bin_count)
        bin_kn[in_bin] = i
    return bin_counts, bin_kn

def get_2_state_bins_all(Q_cutoff,K,N_max,Q_kn):
    print 'Binning Q'
    bin_kn = numpy.zeros([K,N_max],numpy.int16)
    bin_counts = list()

    # unfolded
    in_bin = (Q_kn <= Q_cutoff)
    bin_count = in_bin.sum()
    bin_counts.append(bin_count)
    bin_kn[in_bin] = 0

    # folded
    in_bin = (Q_kn > Q_cutoff)
    bin_count = in_bin.sum()
    bin_counts.append(bin_count)
    bin_kn[in_bin] = 1
    
    return bin_counts, bin_kn

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
    U_kn, Q_kn, N_max = MBAR_pmfQ.read_data(options,K,T)

    # test for statistical inefficiencies
    U_kn, Q_kn, N_k = MBAR_pmfQ.subsample(U_kn, Q_kn, K, N_max)

	# binning data
    nbins = 2
    bin_centers = [.225,.925] # dummy
    bin_counts, bin_kn = get_2_state_bins_all(0.6,K,N_max,Q_kn)
    print '%i bins were populated:' %nbins
    for i in range(nbins):
        print 'bin %5i (%6.1f) %12i conformations' % (i, bin_centers[i], bin_counts[i])

    # compute reduced potential energy of all snapshots at all temperatures
    u_kln = MBAR_pmfQ.get_ukln(N_max,K,U_kn,N_k,beta_k)

    # initialize mbar
    mbar = MBAR_pmfQ.get_mbar(beta_k,U_kn,N_k,u_kln)

    target_temperatures = numpy.arange(295.,360.,5.)
    f_i = numpy.zeros((len(target_temperatures),nbins))
    df_i = []
    for i,temp in enumerate(target_temperatures):
        print 'Calculating the PMF at', temp
        target_beta = 1.0 / (kB * temp)
        u_kn = target_beta * U_kn
        f_i[i,:], d2f_i = mbar.computePMF_states(u_kn, bin_kn, nbins)
        df_i.append(d2f_i)

    results_file = '%s/dG_raw.pkl' % options.direc
    f = file(results_file,'wb')
    print 'Saving target temperatures, bin centers, f_i, df_i to %s' % results_file
    cPickle.dump(target_temperatures,f)
    cPickle.dump(bin_centers,f)
    cPickle.dump(f_i,f)
    cPickle.dump(df_i,f)
    f.close()
	
if __name__ == '__main__':
    main()

