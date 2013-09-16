#!/usr/bin/python2.4

import numpy
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import cPickle
import optparse
import plot_dG
import plot_dG_solution
def main():
    lam = [.1, .2, .3, .4, .5, .525, .55, .6, .625, .675, .7]
    lam = [.1, .2, .3, .35, .4, .45, .5, .6, .7]
    lam = [.2, .3, .35, .4, .45, .5, .55, .6, .7]
    files = ['/home/edz3fz/proteinmontecarlo/results/1PGB/surface/umbrella_lambda%s/dG_raw.pkl' % str(x)[1::] for x in lam]
    colors = cm.summer(numpy.linspace(0,1,len(lam)))
    for i in range(len(lam)):
        plot_dG.main(files[i],lam[i],colors[i])
    for i in range(1,7):
        plt.figure(i)
        plt.legend()
    soln = '/home/edz3fz/proteinmontecarlo/results/1PGB/solution/dG_raw.pkl'
    plot_dG_solution.solution(soln)
    plt.show()
    
if __name__ == '__main__':
    main()
