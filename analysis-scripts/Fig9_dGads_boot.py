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
    #lam = [.1, .15, .2, .25, .3, .35, .4]
    files = ['/home/edz3fz/proteinmontecarlo/results/1PGB/surface/umbrella_lambda%s/dG_raw_varz.pkl' % str(x)[1::] for x in lam]
    files = ['/home/edz3fz/proteinmontecarlo/results/1PGB/surface/umbrella_lambda%s/dG_raw_noint_2.pkl' % str(x)[1::] for x in lam]
    direc = ['/home/edz3fz/proteinmontecarlo/results/1PGB/surface/umbrella_lambda%s/bootstrap' % str(x)[1::] for x in lam]
    colors = cm.cool(numpy.linspace(0,1,len(lam)))
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
    for i in range(2):
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
        ddGads_folded = [numpy.sqrt(x[1,3]) for x in ddG]
        ddGads_folded = numpy.array(ddGads_folded)
        ddSads_folded = (ddGads_folded[0:-2]**2 + ddGads_folded[2::]**2)**.5/10
        ddHads_folded = (ddGads_folded[1:-1]**2 + (temp_sub*ddSads_folded)**2)**.5
    
        ddGads_unfolded = [numpy.sqrt(x[0,2]) for x in ddG]
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
 
    for i in range(2,len(lam)):
   
        target_temperature, dGf[i,:], ddGf[i,:] = read_boot('%s/dGadsf_boot.pkl' % direc[i])
        temp_sub, dSf[i,:], ddSf[i,:] = read_boot('%s/dSadsf_boot.pkl' % direc[i])
        temp_sub, dHf[i,:], ddHf[i,:] = read_boot('%s/dHadsf_boot.pkl' % direc[i])
        target_temperature, dGu[i,:], ddGu[i,:] = read_boot('%s/dGadsu_boot.pkl' % direc[i])
        temp_sub, dSu[i,:], ddSu[i,:] = read_boot('%s/dSadsu_boot.pkl' % direc[i])
        temp_sub, dHu[i,:], ddHu[i,:] = read_boot('%s/dHadsu_boot.pkl' % direc[i])

        target_temperature, dGrel, ddGrel = read_boot('%s/ddGads_boot.pkl' % direc[i])
        temp_sub, dSrel, ddSrel = read_boot('%s/ddSads_boot.pkl' % direc[i])
        temp_sub, dHrel, ddHrel = read_boot('%s/ddHads_boot.pkl' % direc[i])
        plt.rc('text',usetex=True)
        
        fig.subplots_adjust(hspace=0)
        
        ax1 = plt.subplot(311)
        ax1.errorbar(target_temperatures,dGrel,ddGrel,label=r'$\lambda$ = %s' % lam[i], color=colors[i])
        plt.xlim((300,350))
        plt.ylabel(r'$\Delta$G_{adsorption}$')
        box = ax1.get_position()
        ax1.set_position([box.x0,box.y0,box.width*.82,box.height])
        plt.setp(ax1.get_xticklabels(), visible=False)
        
        ax3 = plt.subplot(313)
        plt.xlim((300,350))
        ax3.errorbar(temp_sub,dHrel,ddHrel,label=r'$\lambda$ = %s' % lam[i], color=colors[i])
        plt.xlabel(r'temperature (K)')
        plt.ylabel(r'$\Delta$H_{adsorption}$')
        plt.yticks(numpy.arange(-10,50,10))
        box = ax3.get_position()
        ax3.set_position([box.x0,box.y0,box.width*.82,box.height])
        
        ax2 = plt.subplot(312)
        plt.xlim((300,350))
        ax2.errorbar(temp_sub,dSrel,ddSrel,label=r'$\lambda$ = %s' % lam[i], color=colors[i])
        plt.setp(ax2.get_xticklabels(), visible=False)
        plt.ylabel(r'$\Delta$S_{adsorption}$')
        plt.yticks(numpy.arange(-.05,.15,.05))
        box = ax2.get_position()
        ax2.set_position([box.x0,box.y0,box.width*.82,box.height])
        lgd = ax2.legend(bbox_to_anchor=(1.39,1.08),prop={'size':14})
    

   # plt.savefig('/home/edz3fz/proteinmontecarlo/manuscripts/figures/Fig10_ddGads_boot.pdf')
   # plt.savefig('/home/edz3fz/proteinmontecarlo/manuscripts/figures/Fig10_ddGads_boot.png')

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
    plt.legend(prop={'size':14})
#    box = ax1.get_position()
#    ax1.set_position([box.x0,box.y0,box.width*.82,box.height])
    plt.xlim((.1,.6))

    ax2 = plt.subplot(312)
    i -= 1
    ddSf[0:4,i] = ddSf[5,i]
    ddSu[0:4,i] = ddSu[5,i]
    ax2.errorbar(lam,dSf[:,i],ddSf[:,i],color='k',label='folded')    
    ax2.errorbar(lam,dSu[:,i],ddSu[:,i],color='k',ls='--',label='unfolded')    
    plt.setp(ax2.get_xticklabels(), visible=False)
    plt.ylabel(r'$\Delta$S_{adsorption}$')
#    plt.legend(prop={'size':9})
    plt.ylim((-.2,0))
    plt.yticks(numpy.arange(-.2,0,.04))
#    box = ax2.get_position()
#    ax2.set_position([box.x0,box.y0,box.width*.82,box.height])
    plt.xlim((.1,.6))

    ax3 = plt.subplot(313)
    ddHf[0:4,i] = ddHf[5,i]
    ddHu[0:4,i] = ddHu[5,i]
    ax3.errorbar(lam,dHf[:,i],ddHf[:,i],color='k',label='folded')    
    ax3.errorbar(lam,dHu[:,i],ddHu[:,i],color='k',ls='--',label='unfolded')    
    plt.xlabel(r'$\lambda$')
    plt.ylabel(r'$\Delta$H_{adsorption}$')
#    lgd = plt.legend(prop={'size':9})
    plt.ylim((-100,0))
    plt.yticks(numpy.arange(-100,0,20))
#    lgd = ax2.legend(bbox_to_anchor=(.5,-.1),prop={'size':10},ncol=2)
#    box = ax3.get_position()
#    ax3.set_position([box.x0,box.y0,box.width*.82,box.height])

    plt.xlim((.1,.6))
    
    plt.savefig('/home/edz3fz/proteinmontecarlo/manuscripts/figures/Fig9_dGads_boot.pdf')
    plt.savefig('/home/edz3fz/proteinmontecarlo/manuscripts/figures/Fig9_dGads_boot.png')

    
if __name__ == '__main__':
    main()
    plt.show()
