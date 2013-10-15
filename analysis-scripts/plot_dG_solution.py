#!/usr/bin/python2.4

import numpy
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib as mpl
import cPickle
import optparse
import pdb

def parse_args():
	parser = optparse.OptionParser(description='Plots dG, dS, dH of folding and adsorption')
	parser.add_option('-f','--file', dest='file', help = 'filename of pickled results')
	(options,args) = parser.parse_args()
	return options

def plot(fig,x,y,dy,xlabel,ylabel,label):
    plt.figure(fig)
    plt.errorbar(x,y,dy,label=label,color="k")
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)

def solution(file):
    results_file = file
    print 'Reading %s' % file
    f = open(results_file,'rb')
    target_temperatures = cPickle.load(f)
    bin_centers = cPickle.load(f)
    print 'The 4 states are centered at', bin_centers
    dG = cPickle.load(f)
    ddG = cPickle.load(f)
    f.close()

    # Free energies of folding
    dGf = dG[:,1] - dG[:,0]

    # Entropies of folding
    # dS(t) = -dG/dT = -(G[t+.5h] - G[t-.5h])/h
    dSf = (dGf[2::] - dGf[0:-2])/-10 # finite difference
    temp_sub = target_temperatures[1:-1]

    # Enthalpies of folding
    # dG = dH - TdS
    # dH = dG + TdS
    dHf = dGf[1:-1] + temp_sub*dSf

    # Uncertainties - just propagation of error for now
    ddGf = [numpy.sqrt(x[1,0]) for x in ddG]
    ddGf = numpy.array(ddGf)
    ddSf = (ddGf[0:-2]**2 + ddGf[2::]**2)**.5/10
    ddHf = (ddGf[1:-1]**2 + (temp_sub*ddSf)**2)**.5

    # plotting
#    plot(1,target_temperatures,dGf,ddGf,'temperature (K)','dG_folding (kcal/mol)','')
#    plot(3,temp_sub,dSf,ddSf,'temperature (K)','dS_folding (kcal/molK)','')
#    plot(5,temp_sub,dHf,ddHf,'temperature (K)','dH_folding (kcal/mol)','')

    # making publication figure
    f=plt.figure(1)
    plt.rc('text',usetex=True)
     
    ax1 = plt.subplot(311)
    ax1.errorbar(target_temperatures,dGf,ddGf,color='k',label='solution')
    plt.xlim((300,350))
    plt.ylabel(r'$\Delta$G_{folding}$')
    plt.setp(ax1.get_xticklabels(), visible=False)
 #   plt.legend(prop={'size':6})
 
    ax2 = plt.subplot(312)
    plt.xlim((300,350))
    ax2.errorbar(temp_sub,dSf,ddSf,color='k',label='solution')
    plt.setp(ax2.get_xticklabels(), visible=False)
    plt.ylabel(r'$\Delta$S_{folding}$')
#    plt.legend(prop={'size':6})
 
    ax3 = plt.subplot(313)
    plt.xlim((300,350))
    ax3.errorbar(temp_sub,dHf,ddHf,color='k',label='solution')
    plt.xlabel(r'temperature (K)')
    plt.ylabel(r'$\Delta$H_{folding}$')
#    plt.legend(prop={'size':6})

    f.subplots_adjust(hspace=0)

def solution_pmf(direc):
    T = numpy.arange(300,355,5)
    colors = cm.cool(numpy.linspace(0,1,len(T)))
    plt.figure(7)
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
    plt.xlabel('Q')
    plt.ylabel('PMF (kcal/mol)')
    plt.legend(loc=0, prop={'size':8})


if __name__ == '__main__':
#    options = parse_args()
    direc = '/home/edz3fz/proteinmontecarlo/results/1PGB/solution'
    solution_pmf(direc)
    solution('%s/dG_raw.pkl' % direc)
    plt.show()
