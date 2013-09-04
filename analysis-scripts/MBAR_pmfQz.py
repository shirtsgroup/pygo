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

def get_ukln(args,N_max,K,Z,T,spring_constant,U_kn,z_kn,N_k,beta_k):
    if args.cpt:
        if os.path.exists('%s/u_kln.npy' % args.direc):
            print 'Reading in u_kln matrix from %s/u_kln.npy' % args.direc
            u_kln = numpy.load('%s/u_kln.npy' % args.direc)
            return u_kln
    else:
    	print 'Computing reduced potential energies...'
    	u_kln = numpy.zeros([K,K,N_max], numpy.float32)
    	for k in range(K):
    		for l in range(K):
    			z_index = l/len(T) # z is outer dimension
    			T_index = l%len(T) # T is inner dimension
    			dz = z_kn[k,0:N_k[k]] - Z[z_index]
    			u_kln[k,l,0:N_k[k]] = beta_k[T_index] * U_kn[k,0:N_k[k]] + spring_constant*(dz)**2
    	temp_file = '%s/u_kln.npy' % args.direc
    	print 'Saving u_kln matrix to %s...' % temp_file
    	numpy.save(temp_file, u_kln)
        return u_kln

def get_mbar(args, beta_k, Z, U_kn, N_k, u_kln):
    if args.cpt:
        if os.path.exists('%s/mbar.pkl' % args.direc):
            print 'Reading in mbar object from %s/mbar.pkl' % args.direc
            f = file('%s/mbar.pkl' % args.direc,'rb')
            mbar = cPickle.load(f)
            f.close()
            return mbar
    else:
    	print 'Using WHAM to generate historgram-based initial guess of dimensionless free energies f_k...'
    	beta_k = numpy.array(beta_k.tolist()*len(Z)) 
    	f_k = wham.histogram_wham(beta_k, U_kn, N_k)
    	print 'Initializing MBAR...'
    	mbar = pymbar.MBAR(u_kln, N_k, initial_f_k = f_k, use_optimized=0, verbose=True)
    	mbar_file = '%s/mbar.pkl' % args.direc
    	f = file(mbar_file,'wb')
    	print 'Saving mbar object to %s' % mbar_file
        saving = False
        if saving:
    	    cPickle.dump(mbar,f)
        print 'saving mbar did not work'
        f.close()	
        return mbar

def get_bins(args,nbins_per,K,N_max,indices,Q_kn,z_kn):
    if args.cpt:
        if os.path.exists('%s/bin_data.pkl' % args.direc):
            print 'Reading in bins from %s/bin_data.pkl' % args.direc
            f = file('%s/bin_data.pkl' % args.direc,'rb')
            nbins = cPickle.load(f)
            bin_centers = cPickle.load(f)
            bin_counts = cPickle.load(f)
            bin_kn = cPickle.load(f)
            f.close()
            return nbins, bin_centers, bin_counts, bin_kn
    else:
    	print 'Binning data...'
    	Q_min = 0
    	Q_max = 1
    	dQ = (Q_max - Q_min) / float(nbins_per)
    	z_min = numpy.min(z_kn[indices])
    	z_max = numpy.max(z_kn[indices])
    	dz = (z_max - z_min) / float(nbins_per)
    	bin_kn = numpy.zeros([K,N_max],numpy.int16)
    	bin_counts = []
    	bin_centers = []
    	nbins = 0
    	for i in range(nbins_per):
    		for j in range(nbins_per):
    			z = z_min + dz * (i + 0.5)
    			Q = Q_min + dQ * (j + 0.5)
    			in_bin = (Q-dQ/2 <= Q_kn[indices]) & (Q_kn[indices] < Q+dQ/2) & (z-dz/2 <= z_kn[indices]) & (z_kn[indices] < z+dz/2)
    			bin_count = in_bin.sum()
    			indices_in_bin = (indices[0][in_bin], indices[1][in_bin])
    			if bin_count > 0:
    				bin_centers.append((z,Q))
    				bin_counts.append(bin_count)
    				bin_kn[indices_in_bin] = nbins
    				nbins += 1
        bin_file = '%s/bin_data.pkl' % args.direc
        f = open(bin_file,'wb')
        cPickle.dump(nbins,f)
        cPickle.dump(bin_centers,f)
        cPickle.dump(bin_counts,f)
        cPickle.dump(bin_kn,f)
        f.close()
        return nbins, bin_centers, bin_counts, bin_kn

def read_data(args,K,Z,T,spring_constant):
    U_kn = numpy.empty([K,args.N_max/args.skip], numpy.float64)
    Q_kn = numpy.empty([K,args.N_max/args.skip], numpy.float64)
    z_kn = numpy.empty([K,args.N_max/args.skip], numpy.float64)	
    print "Reading data..."
    i = 0
    for z in Z: 
        for t in T:
            ufile = '%s/%i/energy%i.npy' %(args.direc, z, t)
            data = numpy.load(ufile)[-args.N_max::]
            U_kn[i,:] = data[::args.skip]
            Qfile = '%s/%i/fractionnative%i.npy' %(args.direc, z, t)
            data = numpy.load(Qfile)[-args.N_max::]
            Q_kn[i,:] = data[::args.skip]
            zfile = '%s/%i/z_traj%i_%i.npy' %(args.direc, z, z, t)
            data = numpy.load(zfile)[-args.N_max::]
            z_kn[i,:] = data[::args.skip]
            sfile = '%s/%i/surfenergy%i.npy' %(args.direc,z,t)
            data = numpy.load(sfile)[-args.N_max::]
            if numpy.shape(data)==(args.N_max,2):
                if numpy.all(data[:,0]==data[:,1]):
                    data = data[:,0]
                else:
                    data = numpy.sum(data,axis=1)
            U_kn[i,:] -= data[::args.skip]
            U_kn[i,:] -= spring_constant*(z_kn[i,:] - z)**2
            i += 1
    N_max = args.N_max/args.skip
    return U_kn, Q_kn, z_kn, N_max

def subsample(U_kn,Q_kn,z_kn,K,N_max):
    assume_uncorrelated = False
    if assume_uncorrelated:
        print 'Assuming data is uncorrelated...'
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
            z_kn[k,0:N_k[k]] = z_kn[k,indices]
    return U_kn, Q_kn, z_kn, N_k


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
#    Z = numpy.array([15,16.5,18]) # smaller subset for testing purposes
    print 'distance states are\n', Z	
    K = len(T)*len(Z)

    # read in data
    U_kn, Q_kn, z_kn, N_max = read_data(options, K, Z, T, spring_constant)    

    # produce a histogram of Q and z
#    plt.hist(numpy.reshape(z_kn,N_max*K),400)
#    plt.savefig('%s/z_hist.png' % options.direc)

	# test for statistical inefficiencies
    U_kn, Q_kn, z_kn, N_k = subsample(U_kn, Q_kn, z_kn, K, N_max)

	# generate a list of indices of all configurations in kn-indicing
    mask_kn = numpy.zeros([K,N_max], dtype=numpy.bool)
    for k in range(0,K):
        mask_kn[k,0:N_k[k]] = True
    indices = numpy.where(mask_kn)
	
	# compute reduced potential energy of all snapshots at all temperatures and distances
    u_kln = get_ukln(options, N_max, K, Z, T, spring_constant, U_kn, z_kn, N_k, beta_k)

	# bin data for PMF calculation
    nbins, bin_centers, bin_counts, bin_kn = get_bins(options, nbins_per, K, N_max, indices, Q_kn, z_kn)
    print '%i bins were populated:' %nbins
    for i in range(nbins):
        print 'bin %5i (%6.1f, %6.1f) %12i conformations' % (i, bin_centers[i][0], bin_centers[i][1], bin_counts[i])

	# use WHAM to quickly compute an initial guess of dimensionless free energies f_k
	# then initialize MBAR
    mbar = get_mbar(options, beta_k, Z, U_kn, N_k, u_kln)


    # calculate PMF at the target temperatures
    target_temperatures = [300,325,350]
    print 'Calculating the PMF at', target_temperatures
    
#    f_i = numpy.zeros((nbins,len(target_temperatures)))
#    df_i = numpy.zeros((nbins,len(target_temperatures)))
#    df_i = [] 
    for i,temp in enumerate(target_temperatures):
        target_beta = 1.0 / (kB * temp)
        u_kn = target_beta * U_kn
        f_i, d2f_i = mbar.computePMF_states(u_kn, bin_kn, nbins)
#        imin = f_i.argmin()
#        for j in range(nbins):
#           df_i[j,i] = sqrt(d2f_i[j,imin]) # uncertainty relative to lowest free energy

        pmf_file = '%s/pmf_%i.pkl' % (options.direc, temp)
        f = file(pmf_file,'wb')
        print 'Saving target temperatures, bin centers, f_i, df_i to %s' % pmf_file
        cPickle.dump(temp,f)
        cPickle.dump(bin_centers,f)
        cPickle.dump(f_i,f)
        cPickle.dump(d2f_i,f)
        f.close()

if __name__ == '__main__':
    main()
    #plot in separate script by reading in pmf.pkl



