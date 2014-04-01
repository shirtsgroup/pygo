# Ellen Zhong
# ellen.zhong@virginia.edu
# 03/08/2014

import sys
import numpy
import pymbar # for MBAR analysis
import timeseries # for timeseries analysis
import os
import os.path
import pdb 
import wham
from optparse import OptionParser

def parse_args():
    parser=OptionParser()
    parser.add_option("-r", "--replicas", default=24, type="int",dest="replicas", help="number of replicas (default: 24)")
    parser.add_option("-n", "--N_max", default=100000, type="int",dest="N_max", help="number of data points to read in (default: 100k)")
    parser.add_option("-s", "--skip", default=1, type="int",dest="skip", help="skip every n data points")
    parser.add_option("--direc", dest="direc", help="Qtraj_singleprot.txt file location")
    parser.add_option('-t', "--tfile", dest="tfile", default="/home/edz3fz/proteinmontecarlo/T32.txt", help="file of temperatures (default: T32.txt)")
    parser.add_option('-Q', "--Qfile", dest="Qfile", default="/home/edz3fz/proteinmontecarlo/Q32.txt", help="file of Qpins (default: Q32.txt)")
    parser.add_option("--k_Qpin", type="float", default=10, help="Q umbrella spring constant (default: 10)")
    parser.add_option('--show', action="store_true", default=False, help="show plot at end")
    (args,_) = parser.parse_args()
    return args

def read_data(args, T, Q, K):
    U_kn = numpy.empty([K,args.N_max/args.skip], numpy.float64)
    Q_kn = numpy.empty([K,args.N_max/args.skip], numpy.float64)
    print "Reading data..."
    for i in range(len(T)):
        suffix = '%i_%2.2f' % (int(T[i]), Q[i])
    	ufile = '%s/energy%s.npy' % (args.direc, suffix)
    	data = numpy.load(ufile)[-args.N_max::]
    	U_kn[i,:] = data[::args.skip]
    	Qfile = '%s/fractionnative%s.npy' %(args.direc, suffix)
    	data = numpy.load(Qfile)[-args.N_max::]
    	Q_kn[i,:] = data[::args.skip]
#        if args.surf:
#            sfile = '%s/surfenergy%i.npy' %(args.direc, t)
#            data = numpy.load(sfile)[-args.N_max::]
#            if numpy.shape(data) == (N_max,2):
#                if data[:,0]==data[:,1]:
#                    data = data[:,0]
#                else:
#                    data = numpy.sum(data,axis=1)
#            U_kn[i,:] -= data[::args.skip]
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

def get_ukln(args, N_max, K, Qpin, beta_k, k_Qpin, U_kn, Q_kn, N_k):
    print 'Computing reduced potential energies...'
    u_kln = numpy.zeros([K,K,N_max], numpy.float32)
    for k in range(K):
   		for l in range(K):
   			u_kln[k,l,0:N_k[k]] = beta_k[l] * (U_kn[k,0:N_k[k]] - k_Qpin[k]*(Q_kn[k,0:N_k[k]]-Qpin[k])**2 + k_Qpin[l]*(Q_kn[k,0:N_k[k]]-Qpin[l])**2)
    return u_kln

def get_mbar(beta_k, U_kn, N_k, u_kln):
    print 'Initializing mbar...'
    #f_k = wham.histogram_wham(beta_k, U_kn, N_k)
    try:
        f_k = numpy.loadtxt('f.k.out')
        assert(len(f_k)==len(beta_k))
        mbar = pymbar.MBAR(u_kln, N_k, initial_f_k = f_k, verbose=True)
    except:
        mbar = pymbar.MBAR(u_kln, N_k, verbose=True)
    #mbar = pymbar.MBAR(u_kln, N_k, initial_f_k = f_k, verbose=True)
    return mbar 

def main():
    args = parse_args()
    kB = 0.00831447/4.184  #Boltzmann constant (Gas constant) in kJ/(mol*K)
    dT = 2.5              # Temperature increment for calculating Cv(T)
    
    T = numpy.loadtxt(args.tfile)
    K = len(T)
    Qpin = numpy.loadtxt(args.Qfile)
    k_Qpin = args.k_Qpin*numpy.ones(K)
    print 'Initial temperature states are', T
  
    U_kn, Q_kn, N_max = read_data(args, T, Qpin, K)
    U_kn, Q_kn, N_k = subsample(U_kn, Q_kn, K, N_max)

    # Define new states without Q biasing
    T_new = numpy.arange(250,325,5)
    K_new = len(T_new)
    # Update states
    T = numpy.concatenate((T, T_new))
    Qpin = numpy.concatenate((Qpin, numpy.zeros(K_new)))
    k_Qpin = numpy.concatenate((k_Qpin, numpy.zeros(K_new)))
    K += K_new
    N_k = numpy.concatenate((N_k,numpy.zeros(K_new)))
    U_kn = numpy.concatenate((U_kn,numpy.zeros([K_new,N_max])))
    Q_kn = numpy.concatenate((Q_kn,numpy.zeros([K_new,N_max])))

    beta_k = 1/(kB*T)
    pdb.set_trace()

    u_kln = get_ukln(args, N_max, K, Qpin, beta_k, k_Qpin, U_kn, Q_kn, N_k)
    
    print "Initializing MBAR..."
    # Use Adaptive Method (Both Newton-Raphson and Self-Consistent, testing which is better)
    mbar = get_mbar(beta_k, U_kn, N_k, u_kln)
    
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

    numpy.save(args.direc+'/foldingcurve_umbrella',numpy.array([T, Q, dQ]))
    numpy.save(args.direc+'/heatcap_umbrella',numpy.array([T, Cv, dCv]))
 
    import matplotlib.pyplot as plt
    #ncavg = numpy.average(Q_fromfile, axis=1)


    
    plt.figure(1)
    #plt.plot(T, ncavg, 'ko')
    plt.plot(T[-K_new::],Q[-K_new::],'k')
    plt.errorbar(T[-K_new::], Q[-K_new::], yerr=dQ[-K_new::])
    plt.xlabel('Temperature (K)')
    plt.ylabel('Q fraction native contacts')
    #plt.title('Heat Capacity from Go like model MC simulation of 1BSQ')
    plt.savefig(args.direc+'/foldingcurve.png')
    numpy.save(args.direc+'/foldingcurve',numpy.array([T, Q, dQ]))
    numpy.save(args.direc+'/heatcap',numpy.array([T, Cv, dCv]))
    if args.show:
        plt.show()

if __name__ == '__main__':
    main()
