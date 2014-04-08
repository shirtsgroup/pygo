#!/usr/bin/python2.4

import numpy
from math import *
import pymbar
import timeseries
import commands
import os
import pdb
#import matplotlib.pyplot as plt
import optparse
import wham
import cPickle
import MBAR_pmfQz
import MBAR_pmfQ
import MBAR_4_state_pmf

def parse_args():
    parser = optparse.OptionParser(description='Calculates the PMF(Q,z)')
#	parser.add_option('-t','--temp', dest= 'temp', nargs = 2, type = 'float', help = 'desired temperatures')
    parser.add_option('--tfile', dest='tfile', default='T.txt', help = 'simulation temperature file')
    parser.add_option('--direc', dest='direc', help='directory of simulation data')
    parser.add_option("-n", "--N_max", default=100000, type="int",dest="N_max", help="number of data points to read in (default: 100k)")
    parser.add_option("-s", "--skip", default=1, type="int",dest="skip", help="skip every n data points")
#	parser.add_option('--f_file', dest='f_file', default='', help='free energy filename, if it exists')
    parser.add_option('--cpt', action="store_true", default=False, help="use checkpoint files, if they exist")
    (options,args) = parser.parse_args()
    return options
def main():
    # read in parameters
    options = parse_args()

    # set constants
    kB = 0.00831447/4.184
    nbins_per = 25
    spring_constant = 1  

    # get temperature and distance states
    T = numpy.loadtxt(options.tfile)
#    T = numpy.array([305.,320.,330.,340.]) # smaller subset for testing purposes
    beta_k = 1 / (kB * T)
    print 'temperature states are\n', T
    Z = numpy.arange(9,31.5,1.5)
    Z = numpy.concatenate((Z,numpy.array([33,36,39,42,45,48])))
#    Z = numpy.array([15,16.5,18]) # smaller subset for testing purposes
    print 'distance states are\n', Z	
    K = len(T)*len(Z)

    # read in data
    U_kn, Q_kn, z_kn, N_max = MBAR_pmfQz.read_data(options, K, Z, T, spring_constant)    

	# test for statistical inefficiencies
    U_kn, Q_kn, z_kn, N_k = MBAR_pmfQz.subsample(U_kn, Q_kn, z_kn, K, N_max)

	# generate a list of indices of all configurations in kn-indicing
    mask_kn = numpy.zeros([K,N_max], dtype=numpy.bool)
    for k in range(0,K):
        mask_kn[k,0:N_k[k]] = True
    indices = numpy.where(mask_kn)
	
	# compute reduced potential energy of all snapshots at all temperatures and distances
    u_kln = MBAR_pmfQz.get_ukln(options, N_max, K, Z, T, spring_constant, U_kn, z_kn, N_k, beta_k)

	# bin data for PMF calculation
    nbins, bin_centers, bin_counts, bin_kn = MBAR_pmfQz.get_bins(options, nbins_per, K, N_max, indices, Q_kn, z_kn)
    print '%i bins were populated:' %nbins
    for i in range(nbins):
        print 'bin %5i (%6.1f, %6.1f) %12i conformations' % (i, bin_centers[i][0], bin_centers[i][1], bin_counts[i])

	# use WHAM to quickly compute an initial guess of dimensionless free energies f_k
	# then initialize MBAR
    mbar = MBAR_pmfQz.get_mbar(options, beta_k, Z, U_kn, N_k, u_kln)

    print 'Computing PMF(Q) at 325 K'
    nbins = 25
    target_temperature = 325
    target_beta = 1.0/(kB*target_temperature)
    nbins, bin_centers, bin_counts, bin_kn = MBAR_pmfQ.get_bins(nbins,K,N_max,Q_kn)
    u_kn = target_beta*U_kn
    f_i, d2f_i = mbar.computePMF_states(u_kn, bin_kn, nbins)
    pmf_file = '%s/pmfQ_umbrella_%i.pkl' % (options.direc, target_temperature)
    f = file(pmf_file, 'wb')
    print 'Saving target temperature, bin centers, f_i, df_i to %s' % pmf_file
    cPickle.dump(target_temperature,f)
    cPickle.dump(bin_centers,f)
    cPickle.dump(f_i,f)
    cPickle.dump(d2f_i,f)
    f.close()
#
if __name__ == '__main__':
    main()



