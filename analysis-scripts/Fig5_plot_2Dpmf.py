#!/usr/bin/python2.4

import numpy
import matplotlib.pyplot as plt
import matplotlib
import pickle
import optparse
import pdb

def plot_2D_PMF(subplot,temp,lam,bin_centers,f):
    plt.subplot(subplot)
    X = [data[0] for data in bin_centers]
    Y = [data[1] for data in bin_centers]
    Z = f 
	
    plt.hexbin(X,Y,C=Z, gridsize=20)
#    plt.xlabel('COM distance from surface')
#    plt.ylabel('Q fraction native contacts')
#    plt.title('temperature = %i, lambda = %3.1f' % (temp, lam),fontsize=6)

def main():
    fig = plt.figure(1,(10*.8,9*.8))
    matplotlib.rc('text', usetex = True)
    matplotlib.rc('font', family = 'serif')
    font = {'family' : 'serif',
            'size'   : 'larger'}
    files = []
    lam = numpy.array([.1,.35,.6])
    for l in lam:
        for t in numpy.array([300, 325, 350]):
            files.append('/home/edz3fz/proteinmontecarlo/results/1PGB/surface/umbrella_lambda%s/pmf_%i.pkl' % (str(l)[1::],  t))
    for i,file in enumerate(files): 
        f = open(file,'rb')
        temp = pickle.load(f)
        bin_centers = pickle.load(f)
        f_i = pickle.load(f)
        df_i = pickle.load(f)
        f.close()
        plot_2D_PMF('33%i' %(i+1),temp,lam[i/3],bin_centers,f_i)
    fig.text(.5,.04,'z', ha='center',va='center',fontdict=font)
    fig.text(.06,.5,'Q',ha='center',va='center',rotation='vertical',fontdict=font)
    font = {'family' : 'serif'}
           # 'size'   : 'normal'}
    fig.text(.24,.93, '300 K', ha='center',va='center',fontdict=font)
    fig.text(.52,.93, '325 K', ha='center',va='center',fontdict=font)
    fig.text(.79,.93, '350 K', ha='center',va='center',fontdict=font)
    fig.text(.98,.2, r'$\lambda$ = 0.6', ha='right',va='center',fontdict=font)
    fig.text(.98,.5, r'$\lambda$ = 0.35', ha='right',va='center',fontdict=font)
    fig.text(.98,.8, r'$\lambda$ = 0.1', ha='right',va='center',fontdict=font)
    #file = options.file[0:options.file.rfind('/')]
    #plt.savefig('%s/%i.png' % (file,temp))	
    plt.savefig('/home/edz3fz/proteinmontecarlo/manuscripts/figures/Fig5_plot_2DPMF.eps')

    plt.show()

if __name__ == '__main__':
    main()
