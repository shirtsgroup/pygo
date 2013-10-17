#!/usr/bin/python2.4

import numpy
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib
import cPickle
import optparse
import plot_dG_solution

def main():
    lam = [.1, .15, .2, .25, .3, .35, .4, .45, .5, .55, .6]
    #lam = [.1, .15, .2, .25, .3, .35, .4]
    files = ['/home/edz3fz/proteinmontecarlo/results/1PGB/surface/umbrella_lambda%s/dG_raw_varz.pkl' % str(x)[1::] for x in lam]
    files = ['/home/edz3fz/proteinmontecarlo/results/1PGB/surface/umbrella_lambda%s/dG_raw_noint_2.pkl' % str(x)[1::] for x in lam]
    colors = cm.spring(numpy.linspace(0,1,len(lam)))
    matplotlib.rc('font', family = 'serif')
    font = {'family' : 'serif',
            'size'   : 'larger'}

    dGf = numpy.zeros((len(lam),13))
    dHf = numpy.zeros((len(lam),11))
    dSf = numpy.zeros((len(lam),11))
    dGu = numpy.zeros((len(lam),13))
    dHu = numpy.zeros((len(lam),11))
    dSu = numpy.zeros((len(lam),11))
    ddGf = numpy.zeros((len(lam),13))
    ddHf = numpy.zeros((len(lam),11))
    ddSf = numpy.zeros((len(lam),11))
    ddGu = numpy.zeros((len(lam),13))
    ddHu = numpy.zeros((len(lam),11))
    ddSu = numpy.zeros((len(lam),11))

    fig = plt.figure(10,(7,10))
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
     
        dGf[i,:] = dGads_folded        
        dHf[i,:] = dHads_folded        
        dSf[i,:] = dSads_folded        
        dGu[i,:] = dGads_unfolded        
        dHu[i,:] = dHads_unfolded        
        dSu[i,:] = dSads_unfolded        
        ddGf[i,:] = ddGads_folded        
        ddHf[i,:] = ddHads_folded        
        ddSf[i,:] = ddSads_folded        
        ddGu[i,:] = ddGads_unfolded        
        ddHu[i,:] = ddHads_unfolded        
        ddSu[i,:] = ddSads_unfolded        

        plt.rc('text',usetex=True)
        
        fig.subplots_adjust(hspace=0)
        
        ax1 = plt.subplot(311)
        ax1.errorbar(target_temperatures,dGads_folded-dGads_unfolded,ddGads_unfolded,label=r'$\lambda$ = %s' % lam[i], color=colors[i])
        plt.xlim((300,350))
        plt.ylabel(r'$\Delta$G_{adsorption}$')
        box = ax1.get_position()
        ax1.set_position([box.x0,box.y0,box.width*.82,box.height])
        plt.setp(ax1.get_xticklabels(), visible=False)
        
        ax3 = plt.subplot(313)
        plt.xlim((300,350))
        ax3.errorbar(temp_sub,dHads_folded-dHads_unfolded,ddHads_unfolded,label=r'$\lambda$ = %s' % lam[i], color=colors[i])
        plt.xlabel(r'temperature (K)')
        plt.ylabel(r'$\Delta$H_{adsorption}$')
        plt.yticks(numpy.arange(-10,50,10))
        box = ax3.get_position()
        ax3.set_position([box.x0,box.y0,box.width*.82,box.height])
        
        ax2 = plt.subplot(312)
        plt.xlim((300,350))
        ax2.errorbar(temp_sub,dSads_folded-dSads_unfolded,ddSads_unfolded,label=r'$\lambda$ = %s' % lam[i], color=colors[i])
        plt.setp(ax2.get_xticklabels(), visible=False)
        plt.ylabel(r'$\Delta$S_{adsorption}$')
        plt.yticks(numpy.arange(-.05,.15,.05))
        box = ax2.get_position()
        ax2.set_position([box.x0,box.y0,box.width*.82,box.height])
        lgd = ax2.legend(bbox_to_anchor=(1.29,.98),prop={'size':10})
    

    plt.savefig('/home/edz3fz/proteinmontecarlo/manuscripts/figures/Fig10_ddGads.pdf')

    f=plt.figure(9,(.85*7,10))
    f.subplots_adjust(hspace=0)
    plt.rc('text',usetex=True)
    
    ax1 = plt.subplot(311)
    i=6
    ax1.errorbar(lam,dGf[:,i],ddGf[:,i],color='k',label='folded')    
    ax1.errorbar(lam,dGu[:,i],ddGu[:,i],color='k',ls='--',label='unfolded')    
    plt.ylabel(r'$\Delta$G_{adsorption}$')
    plt.setp(ax1.get_xticklabels(), visible=False)
    #plt.title('Thermodynamics of adsorption')        
    plt.legend(prop={'size':9})
#    box = ax1.get_position()
#    ax1.set_position([box.x0,box.y0,box.width*.82,box.height])
    plt.xlim((.1,.6))

    ax2 = plt.subplot(312)
    i -= 1
    ax2.errorbar(lam,dSf[:,i],ddSf[:,i],color='k',label='folded')    
    ax2.errorbar(lam,dSu[:,i],ddSu[:,i],color='k',ls='--',label='unfolded')    
    plt.setp(ax2.get_xticklabels(), visible=False)
    plt.ylabel(r'$\Delta$S_{adsorption}$')
#    plt.legend(prop={'size':9})
    plt.yticks(numpy.arange(-.2,.05,.05))
#    box = ax2.get_position()
#    ax2.set_position([box.x0,box.y0,box.width*.82,box.height])
    plt.xlim((.1,.6))

    ax3 = plt.subplot(313)
    ax3.errorbar(lam,dHf[:,i],ddHf[:,i],color='k',label='folded')    
    ax3.errorbar(lam,dHu[:,i],ddHu[:,i],color='k',ls='--',label='unfolded')    
    plt.xlabel(r'$\lambda$')
    plt.ylabel(r'$\Delta$H_{adsorption}$')
#    lgd = plt.legend(prop={'size':9})
    plt.yticks(numpy.arange(-100,20,20))
#    lgd = ax2.legend(bbox_to_anchor=(.5,-.1),prop={'size':10},ncol=2)
#    box = ax3.get_position()
#    ax3.set_position([box.x0,box.y0,box.width*.82,box.height])

    plt.xlim((.1,.6))
    
    plt.savefig('/home/edz3fz/proteinmontecarlo/manuscripts/figures/Fig9_dGads.pdf')
    
#    f=plt.figure(4)
#    plt.rc('text',usetex=True)
#    ax1 = plt.subplot(311)
#    for i in range(len(dGf[:,0])):
#        ax1.errorbar(lam,dGf[:,i]-dGu[:,i],ddGf[:,i]-ddGu[:,i],label='T = %i' % target_temperatures[i])    
#    plt.ylabel(r'$\Delta$G_{adsorption}$')
#    plt.setp(ax1.get_xticklabels(), visible=False)
#    plt.title('Thermodynamics of adsorption')        
#    plt.legend(prop={'size':6})
#
#    ax2 = plt.subplot(312)
#    for i in range(len(dSf[:,0])):
#        ax2.errorbar(lam,dSf[:,i]-dSu[:,i],ddSf[:,i]-ddSu[:,i],label='T = %i' % temp_sub[i])    
#    plt.setp(ax2.get_xticklabels(), visible=False)
#    plt.ylabel(r'$\Delta$S_{adsorption}$')
#    plt.legend(prop={'size':6})
#    
#    ax3 = plt.subplot(313)
#    for i in range(len(dSf[:,0])):
#        ax3.errorbar(lam,dHf[:,i]-dHu[:,i],ddHf[:,i]-ddHu[:,i],label='T = %i' % temp_sub[i])    
#    plt.xlabel(r'$\lambda$')
#    plt.ylabel(r'$\Delta$H_{adsorption}$')
#    plt.legend(prop={'size':6})
#
#    f.subplots_adjust(hspace=0)
#
#    f=plt.figure(2)
#    plt.rc('text',usetex=True)
#    ax1 = plt.subplot(311)
#    for i in range(len(dGf[:,0])):
#        ax1.errorbar(lam,dGf[:,i],ddGf[:,i],label='T = %i' % target_temperatures[i])    
#    plt.ylabel(r'$\Delta$G_{adsorption}$')
#    plt.setp(ax1.get_xticklabels(), visible=False)
#    plt.title('Thermodynamics of adsorption')        
#    plt.legend(prop={'size':6})
#
#    ax2 = plt.subplot(312)
#    for i in range(len(dSf[:,0])):
#        ax2.errorbar(lam,dSf[:,i],ddSf[:,i],label='T = %i' % temp_sub[i])    
#    plt.setp(ax2.get_xticklabels(), visible=False)
#    plt.ylabel(r'$\Delta$S_{adsorption}$')
#    plt.legend(prop={'size':6})
#    
#    ax3 = plt.subplot(313)
#    for i in range(len(dSf[:,0])):
#        ax3.errorbar(lam,dHf[:,i],ddHf[:,i],label='T = %i' % temp_sub[i])    
#    plt.xlabel(r'$\lambda$')
#    plt.ylabel(r'$\Delta$H_{adsorption}$')
#    plt.legend(prop={'size':6})
#
#    f.subplots_adjust(hspace=0)
#
#
#    f=plt.figure(3)
#    plt.rc('text',usetex=True)
#    ax1 = plt.subplot(311)
#    for i in range(len(dGf[:,0])):
#        ax1.errorbar(lam,dGu[:,i],ddGu[:,i],label='T = %i' % target_temperatures[i])    
#    plt.ylabel(r'$\Delta$G_{adsorption}$')
#    plt.setp(ax1.get_xticklabels(), visible=False)
#    plt.title('Thermodynamics of adsorption')        
#    plt.legend(prop={'size':6})
#
#    ax2 = plt.subplot(312)
#    for i in range(len(dSf[:,0])):
#        ax2.errorbar(lam,dSu[:,i],ddSu[:,i],label='T = %i' % temp_sub[i])    
#    plt.setp(ax2.get_xticklabels(), visible=False)
#    plt.ylabel(r'$\Delta$S_{adsorption}$')
#    plt.legend(prop={'size':6})
#    
#    ax3 = plt.subplot(313)
#    for i in range(len(dSf[:,0])):
#        ax3.errorbar(lam,dHu[:,i],ddHu[:,i],label='T = %i' % temp_sub[i])    
#    plt.xlabel(r'$\lambda$')
#    plt.ylabel(r'$\Delta$H_{adsorption}$')
#    plt.legend(prop={'size':6})
#
#    f.subplots_adjust(hspace=0)
#


    
if __name__ == '__main__':
    main()
    plt.show()
