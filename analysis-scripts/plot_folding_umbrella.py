#!/usr/bin/python2.4

import numpy
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import cPickle
import optparse
import plot_dG
import plot_dG_solution

def main():
    lam = [.1, .2, .3, .35, .4, .45, .5, .6, .7]
    lam = [.1, .2, .3, .4, .5, .675, .7]
    lam = [.1, .15, .2, .25, .3, .35, .4, .45, .5]#, .55, .6, .65, .7]
    files = ['/home/edz3fz/proteinmontecarlo/results/1PGB/surface/umbrella_lambda%s/foldingcurve_umbrella.npy' % str(x)[1::] for x in lam]
    colors = cm.cool(numpy.linspace(0,1,len(lam)))
    
    plt.figure(1)
    plt.rc('text',usetex=True)
    
    for i in range(len(lam)):
        data = numpy.load(files[i])
        x=17
        plt.errorbar(data[0,-x::],data[1,-x::],data[2,-x::],label=r'$\lambda$ = %s' % lam[i], color=colors[i])
    soln = '/home/edz3fz/proteinmontecarlo/results/1PGB/solution/foldingcurve.npy'
    data = numpy.load(soln)
    plt.errorbar(data[0,:],data[1,:],data[2,:],label='solution', color='k')
    plt.xlabel('temperature (K)')
    plt.ylabel('Q')
    plt.legend(prop={'size':8})
    plt.show()
 
if __name__ == '__main__':
    main()
