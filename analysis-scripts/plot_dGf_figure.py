#!/usr/bin/python2.4

import numpy
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib as mpl
import cPickle
import optparse
import plot_dG_solution

def main():
    lam = [.1, .15, .2, .25, .3, .35, .4, .45, .5, .55]
    files = ['/home/edz3fz/proteinmontecarlo/results/1PGB/surface/umbrella_lambda%s/dG_raw.pkl' % str(x)[1::] for x in lam]
    colors = cm.jet(numpy.linspace(0,1,len(lam)))
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
     
        f=plt.figure(1)
        plt.rc('text',usetex=True)
        
        ax1 = plt.subplot(311)
        ax1.errorbar(target_temperatures,dGf,ddGf,label=r'$\lambda$ = %s' % lam[i], color=colors[i])
        plt.xlim((300,350))
        plt.ylabel(r'$\Delta$G_{folding}$')
        plt.setp(ax1.get_xticklabels(), visible=False)
        
        ax2 = plt.subplot(312)
        plt.xlim((300,350))
        ax2.errorbar(temp_sub,dSf,ddSf,label=r'$\lambda$ = %s' % lam[i], color=colors[i])
        plt.setp(ax2.get_xticklabels(), visible=False)
        plt.ylabel(r'$\Delta$S_{folding}$')
        
        ax3 = plt.subplot(313)
        plt.xlim((300,350))
        ax3.errorbar(temp_sub,dHf,ddHf,label=r'$\lambda$ = %s' % lam[i], color=colors[i])
        plt.xlabel(r'temperature (K)')
        plt.ylabel(r'$\Delta$H_{folding}$')
    
        f.subplots_adjust(hspace=0)



    soln = '/home/edz3fz/proteinmontecarlo/results/1PGB/solution/dG_raw.pkl'
    plot_dG_solution.solution(soln)
    plt.show()
    
if __name__ == '__main__':
    main()
