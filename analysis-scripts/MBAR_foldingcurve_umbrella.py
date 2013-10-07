#!/usr/bin/python2.4

import sys
import numpy
import pymbar # for MBAR analysis
import timeseries # for timeseries analysis
import os
import os.path
import pdb  # for debugging
from optparse import OptionParser
import MBAR_pmfQz
import wham
import MBAR_pmfQ
import cPickle

def parse_args():
    parser=OptionParser()
    #parser.add_option("-t", "--temprange", nargs=2, default=[300.0,450.0], type="float", dest="temprange", help="temperature range of replicas")
    parser.add_option("-r", "--replicas", default=24, type="int",dest="replicas", help="number of replicas (default: 24)")
    parser.add_option("-n", "--N_max", default=100000, type="int",dest="N_max", help="number of data points to read in (default: 100k)")
    parser.add_option("-s", "--skip", default=1, type="int",dest="skip", help="skip every n data points")
    parser.add_option("--direc", dest="direc", help="Qtraj_singleprot.txt file location")
    parser.add_option("--tfile", dest="tfile", default="/home/edz3fz/proteinmontecarlo/T.txt", help="file of temperatures (default: T.txt)")
    parser.add_option('--cpt', action="store_true", default=False, help="use checkpoint files, if they exist")
    (options,args) = parser.parse_args()
    return options

def get_ukln(args,N_max,K,Z,beta_k,spring_constant,U_kn,z_kn,N_k):
    print 'Computing reduced potential energies...'
    u_kln = numpy.zeros([K,K,N_max], numpy.float32)
    for k in range(K):
   		for l in range(K):
   			#z_index = l/(len(T)) # z is outer dimension
   			#T_index = l%(len(T)) # T is inner dimension
   			dz = z_kn[k,0:N_k[k]] - Z[l]
   			u_kln[k,l,0:N_k[k]] = beta_k[l] * (U_kn[k,0:N_k[k]] + spring_constant[l]*(dz)**2)
    return u_kln

def get_mbar(args, beta_k, Z, U_kn, N_k, u_kln):
    if args.cpt:
        if os.path.exists('%s/f_k_foldingcurve.npy' % args.direc):
            print 'Reading in free energies from %s/f_k.npy' % args.direc
            f_k = numpy.load('%s/f_k.npy' % args.direc)
            mbar = pymbar.MBAR(u_kln,N_k,initial_f_k = f_k, maximum_iterations=0,verbose=True,use_optimized=1)
            return mbar
    print 'Using WHAM to generate historgram-based initial guess of dimensionless free energies f_k...'
    #beta_k = numpy.array(beta_k.tolist()*len(Z)) 
    #f_k = wham.histogram_wham(beta_k, U_kn, N_k)
    print 'Initializing MBAR...'
    mbar = pymbar.MBAR(u_kln, N_k, #initial_f_k = f_k, 
                        use_optimized='', verbose=True)
    mbar_file = '%s/f_k_foldingcurve.npy' % args.direc
    print 'Saving free energies to %s' % mbar_file
    saving = True
    if saving:
    	numpy.save(mbar_file, mbar.f_k)
    return mbar

def main():
    options = parse_args()
    kB = 0.00831447/4.184  #Boltzmann constant
    
    T = numpy.loadtxt(options.tfile)
    Z = numpy.arange(9,31.5,1.5)
    print 'Initial temperature states are', T
    print 'Distance states are', Z
    K = len(T)*len(Z)
    spring_constant = numpy.ones(K)
  
    # read in data
    U_kn, Q_kn, z_kn, N_max = MBAR_pmfQz.read_data(options, K, Z, T, spring_constant[0])

    # subsample the data
    U_kn, Q_kn, z_kn, N_k = MBAR_pmfQz.subsample(U_kn,Q_kn,z_kn,K,N_max)
 
    # insert unweighted states
    T_new = numpy.arange(200,410,10)
    T_new = numpy.array([200,225,250,275,300,305,310,315,320,325,330,335,340,345,350,375,400])
    Z_new = numpy.zeros(len(T_new))
    K_new = len(T_new)
    print 'inserting unweighted temperature states', T_new
    
    # update states
    print 'Inserting blank states'
    Z = Z.tolist()
    Z = [x for x in Z for _ in range(len(T))]
    Z = numpy.concatenate((numpy.array(Z),Z_new))
    T = numpy.array(T.tolist()*(K/len(T)))
    T = numpy.concatenate((T,T_new))
    K += K_new
    spring_constant = numpy.concatenate((spring_constant,numpy.zeros(K_new)))
    print 'all temperature states are ', T
    print 'all surface states are ', Z
    print 'there are a total of %i states' % K
    N_k = numpy.concatenate((N_k,numpy.zeros(K_new)))
    U_kn = numpy.concatenate((U_kn,numpy.zeros([K_new,N_max])))
    Q_kn = numpy.concatenate((Q_kn,numpy.zeros([K_new,N_max])))
    z_kn = numpy.concatenate((z_kn,numpy.zeros([K_new,N_max])))

    beta_k = 1/(kB*T)

    u_kln = get_ukln(options, N_max, K, Z, beta_k, spring_constant, U_kn, z_kn, N_k)
    
    print "Initializing MBAR..."
    # Use Adaptive Method (Both Newton-Raphson and Self-Consistent, testing which is better)
    mbar = get_mbar(options,beta_k,Z,U_kn,N_k,u_kln)
    
    print "Computing Expectations for E..."
    (E_expect, dE_expect) = mbar.computeExpectations(u_kln)*(beta_k)**(-1)
    
    print "Computing Expectations for E^2..."
    (E2_expect,dE2_expect) = mbar.computeExpectations(u_kln*u_kln)*(beta_k)**(-2)
    
    print "Computing Expectations for Q..."
    (Q,dQ) = mbar.computeExpectations(Q_kn)
    
    print "Computing Heat Capacity as ( <E^2> - <E>^2 ) / ( R*T^2 )..."
    Cv = numpy.zeros([K], numpy.float64)
    dCv = numpy.zeros([K], numpy.float64)
    for i in range(K):
           Cv[i] = (E2_expect[i] - (E_expect[i]*E_expect[i])) / ( kB * T[i] * T[i])
           dCv[i] = 2*dE_expect[i]**2 / (kB *T[i]*T[i])   # from propagation of error

    numpy.save(options.direc+'/foldingcurve_umbrella',numpy.array([T, Q, dQ]))
    numpy.save(options.direc+'/heatcap_umbrella',numpy.array([T, Cv, dCv]))
    
#    pdb.set_trace()
# 
#    print 'Computing PMF(Q) at 325 K'
#    nbins = 25
#    target_temperature = 325
#    target_beta = 1.0/(kB*target_temperature)
#    nbins, bin_centers, bin_counts, bin_kn = get_bins(nbins,K,N_max,Q_kn)
#    u_kn = target_beta*U_kn
#    f_i, d2f_i = mbar.computePMF_states(u_kn, bin_kn, nbins)
#    pmf_file = '%s/pmfQ_umbrella_%i.pkl' % (options.direc, target_temperature)
#    f = file(pmf_file, 'wb')
#    print 'Saving target temperature, bin centers, f_i, df_i to %s' % pmf_file
#    cPickle.dump(target_temperature,f)
#    cPickle.dump(bin_centers,f)
#    cPickle.dump(f_i,f)
#    cPickle.dump(d2f_i,f)
#    f.close()
#
#    try:
#        import matplotlib.pyplot as plt
#        plt.figure(1)
#        plt.plot(T,Q,'k')
#        plt.errorbar(T, Q, yerr=dQ)
#        plt.xlabel('Temperature (K)')
#        plt.ylabel('Q fraction native contacts')
#        plt.savefig(options.direc+'/foldingcurve_umbrella.png')
#        plt.show()
#    except:
#        pass    
#
if __name__ == '__main__':
    main()
