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
import MBAR_4_state_pmf

def parse_args():
    parser = optparse.OptionParser(description='Calculates the PMF(Q,z)')
    parser.add_option('--tfile', dest='tfile', default='T.txt', help = 'simulation temperature file')
    parser.add_option('--direc', dest='direc', help='directory of simulation data')
    parser.add_option("-n", "--N_max", default=100000, type="int",dest="N_max", help="number of data points to read in (default: 100k)")
    parser.add_option("-s", "--skip", default=1, type="int",dest="skip", help="skip every n data points")
    parser.add_option('--cpt', action="store_true", default=False, help="use checkpoint files, if they exist")
    (options,args) = parser.parse_args()
    return options

def get_mbar(args, beta_k, Z, U_kn, N_k, u_kln):
    print 'Reading in free energies from %s/f_k.npy' % args.direc
    f_k = numpy.load('%s/f_k.npy' % args.direc)
    mbar = pymbar.MBAR(u_kln,N_k,initial_f_k = f_k, verbose=True,use_optimized=1)
    return mbar

def main():
    # read in parameters
    options = parse_args()

    # set constants
    kB = 0.00831447/4.184
    nbins_per = 25
    spring_constant = 1  

    # get temperature and distance states
    T = numpy.loadtxt(options.tfile)
    beta_k = 1 / (kB * T)
    print 'temperature states are\n', T
    Z = numpy.arange(9,31.5,1.5)
    Z = numpy.concatenate((Z,numpy.array([33,36,39,42,45,48])))
    print 'distance states are\n', Z	
    K = len(T)*len(Z)

    # read in data
    U_kn, Q_kn, z_kn, N_max = MBAR_pmfQz.read_data(options, K, Z, T, spring_constant)    

    # produce a histogram of Q and z
#    plt.hist(numpy.reshape(z_kn,N_max*K),400)
#    plt.savefig('%s/z_hist.png' % options.direc)

	# test for statistical inefficiencies
    U_kn, Q_kn, z_kn, N_k = MBAR_pmfQz.subsample(U_kn, Q_kn, z_kn, K, N_max)

	# generate a list of indices of all configurations in kn-indicing (for binning)
    mask_kn = numpy.zeros([K,N_max], dtype=numpy.bool)
    for k in range(0,K):
        mask_kn[k,0:N_k[k]] = True
    indices = numpy.where(mask_kn)
	
    # bootstrap
    boots = range(3,13) # change this as more bootstraps are generated
    for boot in boots:
        print 'Bootstrap %i' % boot
        for k in range(K):
            if N_k[k] > 0:
                booti = numpy.random.randint(N_k[k],size=N_k[k])
                U_kn[k,0:N_k[k]] = U_kn[k,booti]
                Q_kn[k,0:N_k[k]] = Q_kn[k,booti]
                z_kn[k,0:N_k[k]] = z_kn[k,booti]
    
    	# compute reduced potential energy of all snapshots at all temperatures and distances
        u_kln = MBAR_pmfQz.get_ukln(options, N_max, K, Z, T, spring_constant, U_kn, z_kn, N_k, beta_k)
    
    	# bin data for PMF calculation
        nbins = 4
        #bin_centers = [(10.5,.225),(13.5,.925),(28.5,.225),(28.5,.925)]
        bin_centers = [(13.5,.225),(13.5,.925),(40.5,.225),(40.5,.925)]
        Q_cutoff = 0.6
        #bin_counts, bin_kn = get_4_state_bins_varz(Q_cutoff, K, N_max, indices, Q_kn, z_kn)
        bin_counts, bin_kn = MBAR_4_state_pmf.get_4_state_bins_alldata(Q_cutoff, 30, K, N_max, indices, Q_kn, z_kn)
        print '%i bins were populated:' %nbins
        for i in range(nbins):
            print 'bin %5i (%6.1f, %6.1f) %12i conformations' % (i, bin_centers[i][0], bin_centers[i][1], bin_counts[i])
    
        # use free energies from original analysis to initialize resampled MBAR
        mbar = get_mbar(options, beta_k, Z, U_kn, N_k, u_kln)
    
    
        # calculate PMF at the target temperatures
        target_temperatures = numpy.arange(295.,360.,5)
        
        f_i = numpy.zeros((len(target_temperatures),nbins))
        df_i = [] 
        for i,temp in enumerate(target_temperatures):
            print 'Calculating the PMF at', temp
            target_beta = 1.0 / (kB * temp)
            u_kn = target_beta * U_kn
            f_i[i,:], d2f_i = mbar.computePMF_states(u_kn, bin_kn, nbins)
    #        imin = f_i.argmin()
    #        for j in range(nbins):
    #           df_i[j,i] = sqrt(d2f_i[j,imin]) # uncertainty relative to lowest free energy
            df_i.append(d2f_i)
    
    
        results_file = '%s/bootstrap/dG_raw_boot%i.pkl' % (options.direc,boot)
        f = file(results_file,'wb')
        print 'Saving target temperatures, bin centers, f_i, df_i to %s' % results_file
        cPickle.dump(target_temperatures,f)
        cPickle.dump(bin_centers,f)
        cPickle.dump(f_i,f)
        cPickle.dump(df_i,f)
        f.close()

if __name__ == '__main__':
    main()
    #plot in separate script by reading in dG_raw.pkl



