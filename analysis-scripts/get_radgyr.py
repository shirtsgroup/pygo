#!/usr/bin/python2.4

# Ellen Zhong
# ellen.zhong@virginia.edu
# 10/10/2013

import numpy
#from optparse import OptionParser
import pdb
import matplotlib.pyplot as plt
import matplotlib
import pdb

def read_data(direc,K,N,T):
    Q_kn = numpy.empty([K,N])
    print "Reading Q data..."
    for i,t in enumerate(T):
        Qfile = '%s/fractionnative%i.npy' %(direc, t)
        Q_kn[i,:] = numpy.load(Qfile)[-N::]
    return Q_kn

def get_Rg(direc,K,N,N_max,T):
    Rg_kn = numpy.empty([K,N])
    print "Reading traj data..."
    for i,t in enumerate(T):
        print 'Reading traj at %i' % t
        file = open('%s/trajectory%i' %(direc,t),'r')
        for j in range(N_max - N): # skip through beginning
            _ = numpy.load(file)
        nbeads = len(_) # get numbeads
        for j in range(N):
            coord = numpy.load(file)
            Rg_kn[i,j] = Rg(coord,nbeads)
        file.close()
    return Rg_kn

def Rg(coord,nbeads):
    com = numpy.sum(coord,axis=0) / nbeads
    Rg = coord - com
    Rg2 = numpy.sum(Rg**2) / nbeads
    return Rg2

def bin_Rg(nbins,K,N,Q_kn,Rg_kn):
    print 'Binning Q'
    Q_min = 0
    Q_max = 1
    dQ = (Q_max - Q_min) / float(nbins)
    Rg_Q = []
    bin_centers = []
    bin = 0 # number of bins that are populated
    colors = matplotlib.cm.cool(numpy.linspace(0,1,nbins))
    for i in range(nbins):
        Q = Q_min + dQ * (i + 0.5)
        in_bin = (Q-dQ/2 <= Q_kn) & (Q_kn < Q+dQ/2)
        bin_count = in_bin.sum()
        if bin_count > 0:
            data_in_bin = Rg_kn[in_bin]
            plt.hist(data_in_bin,30,color = colors[i])
            print '%5.5f    %5.5f' % (Q, numpy.average(data_in_bin))
            bin_centers.append(Q)
            Rg_Q.append(numpy.average(data_in_bin))
            bin += 1
    return bin, bin_centers, Rg_Q


def main():
    N_max = 3001 # total number of snapshots saved
    N = 2500 # use last 2000 snapshots
    direc = '/home/edz3fz/proteinmontecarlo/results/1PGB/simlog229'
    tfile = '/home/edz3fz/proteinmontecarlo/T.txt'
    T = numpy.loadtxt(tfile)
    K = len(T)
    Q_kn = read_data(direc,K,N,T)
    Rg_kn = get_Rg(direc,K,N,N_max,T)
    Rg_kn = Rg_kn**.5
    nbins = 50
    nbins, bin_centers, Rg_Q = bin_Rg(nbins,K,N,Q_kn,Rg_kn)
    plt.figure(2)
    plt.plot(bin_centers,Rg_Q)
    plt.show()
    numpy.savetxt('test.out', (numpy.array(bin_centers),Rg_Q))

if __name__ == '__main__':
    main()

