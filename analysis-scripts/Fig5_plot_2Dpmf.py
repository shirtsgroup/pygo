#!/usr/bin/python2.4

import numpy
import matplotlib.pyplot as plt
import matplotlib
import pickle
import optparse
import pdb

def plot_2D_PMF(subplot,temp,lam,bin_centers,f):
    ax = plt.subplot(subplot)

    # make free energies relative to desorbed,folded state
    f_max = 0
    for i,bin in enumerate(bin_centers):
        if bin[0] > 40 and bin[1] > .6 and f[i] > f_max:
            f_max = f[i]
    f -= f_max

    X = [data[0] for data in bin_centers]
    Y = [data[1] for data in bin_centers]
    Z = f 
    plt.hexbin(X,Y,C=Z, gridsize=20)
    plt.xticks(numpy.arange(5,55,10))

def main():
    fig = plt.figure(1,(10*.9,9*.8))
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
    fig.text(.23,.93, '300 K', ha='center',va='center',fontdict=font)
    fig.text(.47,.93, '325 K', ha='center',va='center',fontdict=font)
    fig.text(.70,.93, '350 K', ha='center',va='center',fontdict=font)
    fig.text(.87,.2, r'$\lambda$ = 0.6', ha='right',va='center',fontdict=font)
    fig.text(.88,.5, r'$\lambda$ = 0.35', ha='right',va='center',fontdict=font)
    fig.text(.87,.8, r'$\lambda$ = 0.1', ha='right',va='center',fontdict=font)
    #file = options.file[0:options.file.rfind('/')]
    #plt.savefig('%s/%i.png' % (file,temp))	
    
    fig.subplots_adjust(right=0.8)
    cbar_ax = fig.add_axes([.89,.15,.035,.7])
    #matplotlib.colorbar.ColorbarBase(cbar_ax)
    cb = plt.colorbar(cax=cbar_ax)
    cb.set_label('kcal/mol')
    
#    cax,_ = matplotlib.colorbar.make_axes(ax,fraction=.05,pad=.025,shrink=1.8)
#    cbar = matplotlib.colorbar.ColorbarBase(cax)#, ticks = [0,.5,1])
#    cbar.ax.set_yticklabels(['300 K','325 K','350 K'])
    
    plt.savefig('/home/edz3fz/proteinmontecarlo/manuscripts/figures/Fig5_plot_2DPMF.eps')

    plt.show()

if __name__ == '__main__':
    main()
