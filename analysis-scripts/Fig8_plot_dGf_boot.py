#!/usr/bin/python2.4

import numpy
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib
import cPickle
import optparse
import plot_dG_solution
import pdb

def read_boot(file):
    f = open(file,'rb')
    temp = cPickle.load(f)
    X = cPickle.load(f)
    dX = cPickle.load(f)
    f.close()
    return temp,X,dX

def main():
    lam = [.1, .15, .2, .25, .3, .35, .4, .45, .5, .55, .6]
    colors = cm.cool(numpy.linspace(0,1,len(lam)))
    direc = ['/home/edz3fz/proteinmontecarlo/results/1PGB/surface/umbrella_lambda%s/bootstrap' % str(x)[1::] for x in lam]
    
    # data w/o bootstrap first #
    lam = [.1, .15]
#    files = ['/home/edz3fz/proteinmontecarlo/results/1PGB/surface/umbrella_lambda%s/dG_raw_varz.pkl' % str(x)[1::] for x in lam]
    files = ['/home/edz3fz/proteinmontecarlo/results/1PGB/surface/umbrella_lambda%s/dG_raw_noint_2.pkl' % str(x)[1::] for x in lam]
    plt.rc('text',usetex=True)
    matplotlib.rc('font', family = 'serif', size=20)
    fig=plt.figure(1,(7,10))
    soln = '/home/edz3fz/proteinmontecarlo/results/1PGB/solution/dG_raw.pkl'
    plot_dG_solution.solution(soln)
    for i in range(len(lam)):
        print 'Reading %s' % files[i]
        f = open(files[i],'rb')
        target_temperatures = cPickle.load(f)
        bin_centers = cPickle.load(f)
        print 'The 4 states are centered at', bin_centers
        dG = cPickle.load(f)
        ddG = cPickle.load(f)
        f.close()
        # Free energies of folding
        dGf = dG[:,1] - dG[:,0] # free energy of folding while adsorbed
        dGf_des = dG[:,3] - dG[:,2] # free eneryg of folding while desorbed
    
        # Entropies of folding
        # dS(t) = -dG/dT = -(G[t+.5h] - G[t-.5h])/h
        dSf = (dGf[2::] - dGf[0:-2])/-10 # finite difference
        dSf_des = (dGf_des[2::] - dGf_des[0:-2])/-10
        temp_sub = target_temperatures[1:-1]
    
        # Enthalpies of folding
        # dG = dH - TdS
        # dH = dG + TdS
        dHf = dGf[1:-1] + temp_sub*dSf
        dHf_des = dGf_des[1:-1] + temp_sub*dSf_des
    
        # Uncertainties - just propagation of error for now
        ddGf = [numpy.sqrt(x[1,0]) for x in ddG]
        ddGf = numpy.array(ddGf)
        ddSf = (ddGf[0:-2]**2 + ddGf[2::]**2)**.5/10
        ddHf = (ddGf[1:-1]**2 + (temp_sub*ddSf)**2)**.5
    
        ddGf_des = [numpy.sqrt(x[3,2]) for x in ddG]
        ddGf_des = numpy.array(ddGf_des)
        ddSf_des = (ddGf_des[0:-2]**2 + ddGf_des[2::]**2)**.5/10
        ddHf_des = (ddGf_des[1:-1]**2 + (temp_sub*ddSf_des)**2)**.5
     
        fig.subplots_adjust(hspace=0,left=.2)
        
        ax1 = plt.subplot(311)
        ax1.errorbar(target_temperatures,dGf,ddGf,label=r'$\lambda$ = %s' % lam[i], color=colors[i])
        plt.xlim((300,350))
        plt.ylabel(r'$\Delta$G_{folding}$')
        box = ax1.get_position()
        ax1.set_position([box.x0,box.y0,box.width*.82,box.height])
        plt.setp(ax1.get_xticklabels(), visible=False)
        
        ax3 = plt.subplot(313)
        plt.xlim((300,350))
        ax3.errorbar(temp_sub,dHf,ddHf,label=r'$\lambda$ = %s' % lam[i], color=colors[i])
        plt.xlabel(r'temperature (K)')
        plt.ylabel(r'$\Delta$H_{folding}$')
        plt.yticks(numpy.arange(-100,-10,10))
        box = ax3.get_position()
        ax3.set_position([box.x0,box.y0,box.width*.82,box.height])
       
        ax2 = plt.subplot(312)
        plt.xlim((300,350))
        ax2.errorbar(temp_sub,dSf,ddSf,label=r'$\lambda$ = %s' % lam[i], color=colors[i])
        plt.setp(ax2.get_xticklabels(), visible=False)
        plt.ylabel(r'$\Delta$S_{folding}$')
        box = ax2.get_position()
        ax2.set_position([box.x0,box.y0,box.width*.82,box.height])
#        lgd = ax2.legend(bbox_to_anchor=(1.4, 1.08), prop={'size':14})
 
    lam = [.1, .15, .2, .25, .3, .35, .4, .45, .5, .55,.6]
    for i in range(2,len(lam)):
        fig.subplots_adjust(hspace=0)
       
        target_temperatures,dGf,ddGf = read_boot('%s/dGf_boot.pkl' % direc[i]) 
        ax1 = plt.subplot(311)
        ax1.errorbar(target_temperatures,dGf,ddGf,label=r'$\lambda$ = %s' % lam[i], color=colors[i])
        plt.xlim((300,350))
        plt.ylabel(r'$\Delta$G_{folding}$')
        box = ax1.get_position()
        ax1.set_position([box.x0,box.y0,box.width*.82,box.height])
        plt.setp(ax1.get_xticklabels(), visible=False)
        
        temp_sub,dHf,ddHf = read_boot('%s/dHf_boot.pkl' % direc[i]) 
        ax3 = plt.subplot(313)
        plt.xlim((300,350))
        ax3.errorbar(temp_sub,dHf,ddHf,label=r'$\lambda$ = %s' % lam[i], color=colors[i])
        plt.xlabel(r'temperature (K)')
        plt.ylabel(r'$\Delta$H_{folding}$')
        plt.ylim((-100,-20))
        plt.yticks(numpy.arange(-100,-20,20))
        box = ax3.get_position()
        ax3.set_position([box.x0,box.y0,box.width*.82,box.height])
       
        temp_sub,dSf,ddSf = read_boot('%s/dSf_boot.pkl' % direc[i]) 
        ax2 = plt.subplot(312)
        plt.xlim((300,350))
        ax2.errorbar(temp_sub,dSf,ddSf,label=r'$\lambda$ = %s' % lam[i], color=colors[i])
        plt.setp(ax2.get_xticklabels(), visible=False)
        plt.ylim((-.3,-.1))
        plt.yticks(numpy.arange(-.3,-.1,.04))
        plt.ylabel(r'$\Delta$S_{folding}$')
        box = ax2.get_position()
        ax2.set_position([box.x0,box.y0,box.width*.82,box.height])
        lgd = ax2.legend(bbox_to_anchor=(1.35, 1.15), prop={'size':12})
    

    fig.savefig('/home/edz3fz/proteinmontecarlo/manuscripts/figures/Fig8_dGf_boot.pdf')
    fig.savefig('/home/edz3fz/proteinmontecarlo/manuscripts/figures/Fig8_dGf_boot.png')
    plt.show()
    
if __name__ == '__main__':
    main()
