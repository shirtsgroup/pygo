#!/usr/bin/python2.4

import numpy
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib
import cPickle
import optparse
import plot_dG
import plot_dG_solution

def main():
    lam = [.1, .15, .2, .25, .3, .35, .4, .45, .5, .55, .6]#, .65, .7]
    PMF_files = ['/home/edz3fz/proteinmontecarlo/results/1PGB/surface/umbrella_lambda%s/pmfQ_umbrella_325.pkl' % str(x)[1::] for x in lam]
    colors = cm.spring(numpy.linspace(0,1,len(lam)))
    
    f = plt.figure(1)
    plt.rc('text',usetex=True)
    matplotlib.rc('font', family = 'serif')
    font = {'family' : 'serif',
            'size'   : 'larger'}
    for i in range(len(lam)):
        data = numpy.load(PMF_files[i])
        plt.errorbar(data[0,-x::],data[1,-x::],data[2,-x::],label=r'$\lambda$ = %s' % lam[i], color=colors[i])
    soln = '/home/edz3fz/proteinmontecarlo/results/1PGB/solution/pmf_325.pkl'
    data = numpy.load(soln)
    plt.errorbar(data[0,:],data[1,:],data[2,:],label='solution', color='k')
    plt.ylabel('PMF')
    plt.xlabel('Q')
    plt.legend(prop={'size':8})
    #plt.savefig('/home/edz3fz/proteinmontecarlo/manuscripts/figures/Fig7_foldingcurve.eps')
    plt.show()
 
if __name__ == '__main__':
    main()
