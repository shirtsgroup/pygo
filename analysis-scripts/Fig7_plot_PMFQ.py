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
    lam = [.1, .15, .2, .25, .3, .35, .45, .5, .55, .6]#, .65, .7]
    PMF_files = ['/home/edz3fz/proteinmontecarlo/results/1PGB/surface/umbrella_lambda%s/pmfQ_umbrella_325.pkl' % str(x)[1::] for x in lam]
    PMF_files[-2] = '/home/edz3fz/proteinmontecarlo/results/1PGB/surface/lambda.55/pmf_325.pkl'
    #PMF_files[-1] = '/home/edz3fz/proteinmontecarlo/results/1PGB/surface/lambda.6/pmf_325.pkl'
    colors = cm.cool(numpy.linspace(0,1,len(lam)))
    
    f = plt.figure(1)
    plt.rc('text',usetex=True)
    matplotlib.rc('font', family = 'serif')
    font = {'family' : 'serif',
            'size'   : 'larger'}
    for i in range(len(lam)):
        f = open(PMF_files[i],'rb')
        temp = cPickle.load(f)
        bin_centers = cPickle.load(f)
        f_i = cPickle.load(f)
        df_i = cPickle.load(f)
        f.close()
        idx = numpy.argmin(f_i[int(.7*len(f_i))::]) + int(.7*len(f_i))
        f_i -= f_i[idx]
        plt.errorbar(bin_centers[1::],f_i[1::],numpy.array(df_i[0,:])[0,1::],label=r'$\lambda$ = %s' % lam[i], color=colors[i])
 
    soln = '/home/edz3fz/proteinmontecarlo/results/1PGB/solution/pmf_325.pkl'
    f = open(soln,'rb')
    temp = cPickle.load(f)
    bin_centers = cPickle.load(f)
    f_i = cPickle.load(f)
    df_i = cPickle.load(f)
    f.close()
    idx = numpy.argmin(f_i[int(.7*len(f_i))::]) + int(.7*len(f_i))
    f_i -= f_i[idx]
    plt.errorbar(bin_centers[1::],f_i[1::],numpy.array(df_i[0,:])[0,1::],label='solution', color='k')
 
    plt.ylabel('PMF')
    plt.xlabel('Q')
    plt.legend(prop={'size':10},loc=4)
    plt.savefig('/home/edz3fz/proteinmontecarlo/manuscripts/figures/Fig7_pmfQ.pdf')
    plt.show()
 
if __name__ == '__main__':
    main()
