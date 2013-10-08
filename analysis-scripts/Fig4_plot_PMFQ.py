#!/usr/bin/python2.4

import numpy
import cPickle
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import matplotlib

def plot_solution_pmf(direc):
    # differs from plot_dG_solution method only by the color 
    matplotlib.rc('text', usetex = True)
    matplotlib.rc('font', family = 'serif')
    font = {'family' : 'serif',
            'size'   : 'larger'}
    T = numpy.arange(300,355,5)
    colors = cm.cool(numpy.linspace(0,1,len(T)))
    fig = plt.figure(1)
    ax = fig.add_subplot(1,1,1)
    for i,t in enumerate(T):
        f = open('%s/pmf_%i.pkl' %(direc,t),'rb')
        temp = cPickle.load(f)
        bin_centers = cPickle.load(f)
        f_i = cPickle.load(f)
        df_i = cPickle.load(f)
        f.close()
        idx = numpy.argmin(f_i[int(.7*len(f_i))::]) + int(.7*len(f_i))
        f_i -= f_i[idx]
        plt.errorbar(bin_centers[1::],f_i[1::],numpy.array(df_i[0,:])[0,1::],label='T = %i K' % t, color=colors[i])
    #plt.legend(loc=0, prop={'size':8})
    plt.xlabel('Q')
    plt.ylabel('PMF (kcal/mol)')
    
    #ax = fig.add_axes([.9,.25,.03,.5])
    cax,_ = matplotlib.colorbar.make_axes(ax,fraction=.05,pad=.025,shrink=.8)
    cbar = matplotlib.colorbar.ColorbarBase(cax,cmap=cm.cool, ticks = [0,.5,1])
    cbar.ax.set_yticklabels(['300 K','325 K','350 K'])

def main():
    file = '/home/edz3fz/proteinmontecarlo/results/1PGB/solution'
    plot_solution_pmf(file)
    plt.savefig('/home/edz3fz/proteinmontecarlo/manuscripts/figures/Fig4_PMFQ.eps')
    plt.show()

if __name__=='__main__':
    main()
