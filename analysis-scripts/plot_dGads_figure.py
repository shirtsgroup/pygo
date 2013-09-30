#!/usr/bin/python2.4

import numpy
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib as mpl
import cPickle
import optparse
import plot_dG_solution

def main():
    lam = [.1, .15, .2, .3, .35, .4, .45, .5, .55, .6]
    files = ['/home/edz3fz/proteinmontecarlo/results/1PGB/surface/umbrella_lambda%s/dG_raw.pkl' % str(x)[1::] for x in lam]
    colors = cm.spring(numpy.linspace(0,1,len(lam)))

    dGf = numpy.zeros(len(lam))
    dHf = numpy.zeros(len(lam))
    dSf = numpy.zeros(len(lam))
    dGu = numpy.zeros(len(lam))
    dHu = numpy.zeros(len(lam))
    dSu = numpy.zeros(len(lam))
    ddGf = numpy.zeros(len(lam))
    ddHf = numpy.zeros(len(lam))
    ddSf = numpy.zeros(len(lam))
    ddGu = numpy.zeros(len(lam))
    ddHu = numpy.zeros(len(lam))
    ddSu = numpy.zeros(len(lam))

    for i in range(len(lam)):
        print 'Reading %s' % files[i]
        f = open(files[i],'rb')
        target_temperatures = cPickle.load(f)
        bin_centers = cPickle.load(f)
        print 'The 4 states are centered at', bin_centers
        dG = cPickle.load(f)
        ddG = cPickle.load(f)
        f.close()
    
        # Free energies of adsorption
        dGads_folded = dG[:,1] - dG[:,3] # free energy of adsorption while adsorbed
        dGads_unfolded = dG[:,0] - dG[:,2] # free eneryg of adsorption while desorbed
    
        # Entropies of adsorption
        # dS(t) = -dG/dT = -(G[t+.5h] - G[t-.5h])/h
        dSads_folded = (dGads_folded[2::] - dGads_folded[0:-2])/-10 # finite difference
        dSads_unfolded = (dGads_unfolded[2::] - dGads_unfolded[0:-2])/-10
        temp_sub = target_temperatures[1:-1]
    
        # Enthalpies of adsorption
        # dG = dH - TdS
        # dH = dG + TdS
        dHads_folded = dGads_folded[1:-1] + temp_sub*dSads_folded
        dHads_unfolded = dGads_unfolded[1:-1] + temp_sub*dSads_unfolded
    
        # Uncertainties - just propagation of error for now
        ddGads_folded = [numpy.sqrt(x[1,0]) for x in ddG]
        ddGads_folded = numpy.array(ddGads_folded)
        ddSads_folded = (ddGads_folded[0:-2]**2 + ddGads_folded[2::]**2)**.5/10
        ddHads_folded = (ddGads_folded[1:-1]**2 + (temp_sub*ddSads_folded)**2)**.5
    
        ddGads_unfolded = [numpy.sqrt(x[3,2]) for x in ddG]
        ddGads_unfolded = numpy.array(ddGads_unfolded)
        ddSads_unfolded = (ddGads_unfolded[0:-2]**2 + ddGads_unfolded[2::]**2)**.5/10
        ddHads_unfolded = (ddGads_unfolded[1:-1]**2 + (temp_sub*ddSads_unfolded)**2)**.5
     
        dGf[i] = dGads_folded[6]        
        dHf[i] = dHads_folded[5]        
        dSf[i] = dSads_folded[5]        
        dGu[i] = dGads_unfolded[6]        
        dHu[i] = dHads_unfolded[5]        
        dSu[i] = dSads_unfolded[5]        
        ddGf[i] = ddGads_folded[6]        
        ddHf[i] = ddHads_folded[5]        
        ddSf[i] = ddSads_folded[5]        
        ddGu[i] = ddGads_unfolded[6]        
        ddHu[i] = ddHads_unfolded[5]        
        ddSu[i] = ddSads_unfolded[5]        

        f=plt.figure(1)
        plt.rc('text',usetex=True)
        
        ax1 = plt.subplot(311)
        ax1.errorbar(target_temperatures,dGads_folded,ddGads_folded,label=r'$\lambda$ = %s' % lam[i], color=colors[i])
        plt.xlim((300,350))
        plt.ylabel(r'$\Delta$G_{adsorption}$')
        plt.setp(ax1.get_xticklabels(), visible=False)
        plt.legend(prop={'size':6})
        plt.title('Thermodynamics of adsorption, folded')        

        ax2 = plt.subplot(312)
        plt.xlim((300,350))
        ax2.errorbar(temp_sub,dSads_folded,ddSads_folded,label=r'$\lambda$ = %s' % lam[i], color=colors[i])
        plt.setp(ax2.get_xticklabels(), visible=False)
        plt.ylabel(r'$\Delta$S_{adsorption}$')
        plt.legend(prop={'size':6})
        
        ax3 = plt.subplot(313)
        plt.xlim((300,350))
        ax3.errorbar(temp_sub,dHads_folded,ddHads_folded,label=r'$\lambda$ = %s' % lam[i], color=colors[i])
        plt.xlabel(r'temperature (K)')
        plt.ylabel(r'$\Delta$H_{adsorption}$')
        plt.legend(prop={'size':6})
    
        f.subplots_adjust(hspace=0)

        f=plt.figure(3)
        plt.rc('text',usetex=True)
        
        ax1 = plt.subplot(311)
        ax1.errorbar(target_temperatures,dGads_folded-dGads_unfolded,ddGads_unfolded,label=r'$\lambda$ = %s' % lam[i], color=colors[i])
        plt.xlim((300,350))
        plt.ylabel(r'$\Delta$G_{adsorption}$')
        plt.setp(ax1.get_xticklabels(), visible=False)
        plt.legend(prop={'size':6})
        plt.title('Relative thermodynamics of adsorption (folded-unfolded)')        
        
        ax2 = plt.subplot(312)
        plt.xlim((300,350))
        ax2.errorbar(temp_sub,dSads_folded-dSads_unfolded,ddSads_unfolded,label=r'$\lambda$ = %s' % lam[i], color=colors[i])
        plt.setp(ax2.get_xticklabels(), visible=False)
        plt.ylabel(r'$\Delta$S_{adsorption}$')
        plt.legend(prop={'size':6})
        
        ax3 = plt.subplot(313)
        plt.xlim((300,350))
        ax3.errorbar(temp_sub,dHads_folded-dHads_unfolded,ddHads_unfolded,label=r'$\lambda$ = %s' % lam[i], color=colors[i])
        plt.xlabel(r'temperature (K)')
        plt.ylabel(r'$\Delta$H_{adsorption}$')
        plt.legend(prop={'size':6})
    
        f.subplots_adjust(hspace=0)

        f=plt.figure(2)
        plt.rc('text',usetex=True)
        
        ax1 = plt.subplot(311)
        ax1.errorbar(target_temperatures,dGads_unfolded,ddGads_unfolded,label=r'$\lambda$ = %s' % lam[i], color=colors[i])
        plt.xlim((300,350))
        plt.ylabel(r'$\Delta$G_{adsorption}$')
        plt.setp(ax1.get_xticklabels(), visible=False)
        plt.legend(prop={'size':6})
        plt.title('Thermodynamics of adsorption, unfolded')        
        
        ax2 = plt.subplot(312)
        plt.xlim((300,350))
        ax2.errorbar(temp_sub,dSads_unfolded,ddSads_unfolded,label=r'$\lambda$ = %s' % lam[i], color=colors[i])
        plt.setp(ax2.get_xticklabels(), visible=False)
        plt.ylabel(r'$\Delta$S_{adsorption}$')
        plt.legend(prop={'size':6})
        
        ax3 = plt.subplot(313)
        plt.xlim((300,350))
        ax3.errorbar(temp_sub,dHads_unfolded,ddHads_unfolded,label=r'$\lambda$ = %s' % lam[i], color=colors[i])
        plt.xlabel(r'temperature (K)')
        plt.ylabel(r'$\Delta$H_{adsorption}$')
        plt.legend(prop={'size':6})
    
        f.subplots_adjust(hspace=0)


    #soln = '/home/edz3fz/proteinmontecarlo/results/1PGB/solution/dG_raw.pkl'
    #plot_dG_solution.solution(asoln)
    
    f=plt.figure(4)
    plt.rc('text',usetex=True)
    ax1 = plt.subplot(311)
    ax1.errorbar(lam,dGf,ddGf,label='folded')    
    ax1.errorbar(lam,dGu,ddGu,label='unfolded')    
    plt.ylabel(r'$\Delta$G_{adsorption}$')
    plt.setp(ax1.get_xticklabels(), visible=False)
    plt.title('Thermodynamics of adsorption')        
    plt.legend(prop={'size':6})

    ax2 = plt.subplot(312)
    ax2.errorbar(lam,dSf,ddSf,label='folded')    
    ax2.errorbar(lam,dSu,ddSu,label='unfolded')    
    plt.setp(ax2.get_xticklabels(), visible=False)
    plt.ylabel(r'$\Delta$S_{adsorption}$')
    plt.legend(prop={'size':6})
    
    ax3 = plt.subplot(313)
    ax3.errorbar(lam,dHf,ddHf,label='folded')    
    ax3.errorbar(lam,dHu,ddHu,label='unfolded')    
    plt.xlabel(r'$\lambda$')
    plt.ylabel(r'$\Delta$H_{adsorption}$')
    plt.legend(prop={'size':6})

    f.subplots_adjust(hspace=0)

    plt.show()
    
if __name__ == '__main__':
    main()
