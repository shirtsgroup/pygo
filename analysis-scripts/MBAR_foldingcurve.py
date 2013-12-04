#!/usr/bin/python2.4

import sys
import numpy
import pymbar # for MBAR analysis
import timeseries # for timeseries analysis
import os
import os.path
import pdb  # for debugging
import wham
from optparse import OptionParser

def parse_args():
    parser=OptionParser()
    #parser.add_option("-t", "--temprange", nargs=2, default=[300.0,450.0], type="float", dest="temprange", help="temperature range of replicas")
    parser.add_option("-r", "--replicas", default=24, type="int",dest="replicas", help="number of replicas (default: 24)")
    parser.add_option("-n", "--N_max", default=100000, type="int",dest="N_max", help="number of data points to read in (default: 100k)")
    parser.add_option("-s", "--skip", default=1, type="int",dest="skip", help="skip every n data points")
    parser.add_option("--direc", dest="direc", help="Qtraj_singleprot.txt file location")
    parser.add_option("--tfile", dest="tfile", default="/home/edz3fz/proteinmontecarlo/T.txt", help="file of temperatures (default: T.txt)")
    parser.add_option('--show', action="store_true", default=False, help="show plot at end")
    (options,args) = parser.parse_args()
    return options

def read_data(args,T,K):
    U_kn = numpy.empty([K,args.N_max/args.skip], numpy.float64)
    Q_kn = numpy.empty([K,args.N_max/args.skip], numpy.float64)
    print "Reading data..."
    for i, t in enumerate(T):
    	ufile = '%s/energy%i.npy' %(args.direc, t)
    	data = numpy.load(ufile)[-args.N_max::]
    	U_kn[i,:] = data[::args.skip]
    	Qfile = '%s/fractionnative%i.npy' %(args.direc, t)
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
    return U_kn,Q_kn,N_max

def main():
    options = parse_args()
    kB = 0.00831447/4.184  #Boltzmann constant (Gas constant) in kJ/(mol*K)
    dT = 2.5              # Temperature increment for calculating Cv(T)
    
    T = numpy.loadtxt(options.tfile)
    print 'Initial temperature states are', T
    K = len(T)
  
    U_kn, Q_kn, N_max = read_data(options,T,K)

    print 'Subsampling Q...' 
    N_k = numpy.zeros(K,numpy.int32)
    g = numpy.zeros(K,numpy.float64)
    for k in range(K):  # subsample the energies
       g[k] = timeseries.statisticalInefficiency(Q_kn[k])#,suppress_warning=True)
       indices = numpy.array(timeseries.subsampleCorrelatedData(Q_kn[k],g=g[k])) # indices of uncorrelated samples
       N_k[k] = len(indices) # number of uncorrelated samplesadsf
       print '%i uncorrelated samples out of %i total samples' %(len(indices),options.N_max/options.skip)
       U_kn[k,0:N_k[k]] = U_kn[k,indices]
       Q_kn[k,0:N_k[k]] = Q_kn[k,indices]

    insert = True
    if insert: 
        #------------------------------------------------------------------------
        # Insert Intermediate T's and corresponding blank U's and E's
        #------------------------------------------------------------------------
        # Set up variables
        Temp_k = T
        currentT = T[0] + dT
        maxT = T[-1]
        i = 1
        
        print("--Inserting intermediate temperatures...")
        
        # Loop, inserting T's at which we are interested in the properties
        while (currentT < maxT) :
               if (currentT < Temp_k[i]):
                     Temp_k = numpy.insert(Temp_k, i, currentT)
                     currentT = currentT + dT
               else:
                     currentT = Temp_k[i] + dT
                     i = i + 1
                     
        # Update number of states
        K = len(Temp_k)
        
        print("--Inserting blank energies to match up with inserted temperatures...")
        
        # Loop, inserting E's into blank matrix (leaving blanks only where new Ts are inserted)
        Q_fromfile = Q_kn
        Nall_k = numpy.zeros([K], numpy.int32) # Number of samples (n) for each state (k) = number of iterations/energies
        E_kn = numpy.zeros([K, N_max], numpy.float64)
        Q_kn = numpy.zeros([K, N_max], numpy.float64)
        i = 0
        
        for k in range(K):
            if (Temp_k[k] == T[i]):
                E_kn[k,0:N_k[i]] = U_kn[i,0:N_k[i]]
                Q_kn[k,0:N_k[i]] = Q_fromfile[i,0:N_k[i]]
                Nall_k[k] = N_k[i]
                i = i + 1
    else:
        print 'Not inserting intermediate temperatures'
        Temp_k = T
        E_kn = U_kn
        Nall_k = N_k

    #------------------------------------------------------------------------
    # Compute inverse temperatures
    #------------------------------------------------------------------------
    beta_k = 1 / (kB * Temp_k)
    
    #------------------------------------------------------------------------
    # Compute reduced potential energies
    #------------------------------------------------------------------------
    
    print "--Computing reduced energies..."
    
    u_kln = numpy.zeros([K,K,N_max], numpy.float64) # u_kln is reduced pot. ener. of segment n of temp k evaluated at temp l
    
    for k in range(K):
           for l in range(K):
                 u_kln[k,l,0:Nall_k[k]] = beta_k[l] * E_kn[k,0:Nall_k[k]]
    
    #------------------------------------------------------------------------
    # Initialize MBAR
    #------------------------------------------------------------------------
    
    # Initialize MBAR with Newton-Raphson
    print ""
    print "Initializing MBAR:"
    print "--K = number of Temperatures"
    print "--L = number of Temperatures"
    print "--N = number of Energies per Temperature"
    
    # Use Adaptive Method (Both Newton-Raphson and Self-Consistent, testing which is better)
    if insert:
        mbar = pymbar.MBAR(u_kln, Nall_k, method = 'adaptive', verbose=True, relative_tolerance=1e-12)
    else:
        f_k = wham.histogram_wham(beta_k, U_kn, Nall_k, relative_tolerance = 1.0e-4)
        mbar = pymbar.MBAR(u_kln, Nall_k, initial_f_k = f_k, verbose=True)
    #------------------------------------------------------------------------
    # Compute Expectations for E_kt and E2_kt as E_expect and E2_expect
    #------------------------------------------------------------------------
    print ""
    print "Computing Expectations for E..."
    (E_expect, dE_expect) = mbar.computeExpectations(u_kln)*(beta_k)**(-1)
    print "Computing Expectations for E^2..."
    (E2_expect,dE2_expect) = mbar.computeExpectations(u_kln*u_kln)*(beta_k)**(-2)
    
    print "Computing Expectations for Q..."
    (Q,dQ) = mbar.computeExpectations(Q_kn)
    
    #------------------------------------------------------------------------
    # Compute Cv for NVT simulations as <E^2> - <E>^2 / (RT^2)
    #------------------------------------------------------------------------
    #print ""
    #print "Computing Heat Capacity as ( <E^2> - <E>^2 ) / ( R*T^2 )..."
    
    Cv_expect = numpy.zeros([K], numpy.float64)
    dCv_expect = numpy.zeros([K], numpy.float64)
    
    for i in range(K):
           Cv_expect[i] = (E2_expect[i] - (E_expect[i]*E_expect[i])) / ( kB * Temp_k[i] * Temp_k[i])
           dCv_expect[i] = 2*dE_expect[i]**2 / (kB *Temp_k[i]*Temp_k[i])   # from propagation of error
    
    #print "Temperature  dA         <E> +/- d<E>       <E^2> +/- d<E^2>       Cv +/- dCv"     
    #print "-------------------------------------------------------------------------------"
    #for k in range(K):
    #       print "%8.3f %8.3f %9.3f +/- %5.3f  %9.1f +/- %5.1f   %7.4f +/- %6.4f" % (Temp_k[k],mbar.f_k[k],E_expect[k],dE_expect[k],E2_expect[k],dE2_expect[k],Cv_expect[k], dCv_expect[k])
    #numpy.savetxt('/home/edz3fz/Qsurf_int.txt',Q)
    #numpy.savetxt('/home/edz3fz/dQsurf_int.txt',dQ)
    #numpy.savetxt('/home/edz3fz/dQsol.txt',dQ)
    
    #numpy.savetxt('/home/edz3fz/Qtemp.tt',Temp_k)
    import matplotlib.pyplot as plt
    #ncavg = numpy.average(Q_fromfile, axis=1)
    
    plt.figure(1)
    #plt.plot(T, ncavg, 'ko')
    plt.plot(Temp_k,Q,'k')
    plt.errorbar(Temp_k, Q, yerr=dQ)
    plt.xlabel('Temperature (K)')
    plt.ylabel('Q fraction native contacts')
    #plt.title('Heat Capacity from Go like model MC simulation of 1BSQ')
    plt.savefig(options.direc+'/foldingcurve.png')
    numpy.save(options.direc+'/foldingcurve',numpy.array([Temp_k, Q, dQ]))
    numpy.save(options.direc+'/heatcap',numpy.array([Temp_k, Cv_expect, dCv_expect]))
    if options.show:
        plt.show()

if __name__ == '__main__':
    main()
