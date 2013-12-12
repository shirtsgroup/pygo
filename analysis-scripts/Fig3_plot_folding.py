#!/usr/bin/python2.4

import numpy
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib
import cPickle
import optparse
import plot_dG
import plot_dG_solution
import pdb

def main():
    f,ax1 = plt.subplots()
    plt.rc('text',usetex=True)
    matplotlib.rc('font', family = 'serif')
    #plt.rc('text.latex', preamble='\usepackage{sfmath}')
    font = {'family' : 'serif',
            'size'   : 'larger'}

    Q = '/home/edz3fz/proteinmontecarlo/results/1PGB/solution/foldingcurve.npy'
    data = numpy.load(Q)
    ln1 = ax1.errorbar(data[0,:],data[1,:],data[2,:], color='k',label='Q')
    ax1.set_ylabel('Q', fontdict = font)
    ax1.set_xlabel('temperature (K)', fontdict = font)

    ax2 = ax1.twinx()
    Cv = '/home/edz3fz/proteinmontecarlo/results/1PGB/solution/heatcap.npy'
    data = numpy.load(Cv)
    ln2 = ax2.errorbar(data[0,:],data[1,:],data[2,:], color='k',ls='--',label=r'C_{v}')
    ax2.set_ylabel(r'C_{v}  \textnormal{ (kcal/mol$\cdot$K)}', fontdict = font)

    ln = ln1 + ln2
    plt.legend(ln,['Q',r'C_{v}'])#,prop={'size':11})
    plt.savefig('/home/edz3fz/proteinmontecarlo/manuscripts/figures/Fig3_solution.pdf')

    plt.show()
 
if __name__ == '__main__':
    main()
