#!/usr/bin/python2.4

import numpy
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import cPickle
import optparse

def parse_args():
	parser = optparse.OptionParser(description='Plots dG, dS, dH of folding and adsorption')
	parser.add_option('-f','--file', dest='file', help = 'filename of pickled results')
	parser.add_option('-l','--lam', type='float', dest='lam', help = 'lambda value of surface attractiveness')
	(options,args) = parser.parse_args()
	return options

def plot_2D_PMF(temp,lam,bin_centers,f):
	X = [data[0] for data in bin_centers]
	Y = [data[1] for data in bin_centers]
	Z = f 
	
	plt.hexbin(X,Y,C=Z)
	plt.xlabel('COM distance from surface')
	plt.ylabel('Q fraction native contacts')
	plt.title('temperature = %i, lambda = %3.1f' % (temp, lam))

def plot(fig,x,y,dy,xlabel,ylabel,label,color):
    plt.figure(fig)
    plt.errorbar(x,y,dy,label=label, color=color)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)

def solution(file, lam):
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
    plot(1,target_temperatures,dGf,ddGf,'temperature (K)','dG_folding (kcal/mol)','adsorbed, lam = %s' % lam)
    plot(2,temp_sub,dSf,ddSf,'temperature (K)','dS_folding (kcal/molK)','adsorbed, lam = %s' % lam)
    plot(3,temp_sub,dHf,ddHf,'temperature (K)','dH_folding (kcal/mol)','adsorbed, lam = %s' % lam)

def main(file, lam, color):
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
    dGf_ads = dG[:,1] - dG[:,0] # free energy of folding while adsorbed
    dGf_des = dG[:,3] - dG[:,2] # free eneryg of folding while desorbed

    # Entropies of folding
    # dS(t) = -dG/dT = -(G[t+.5h] - G[t-.5h])/h
    dSf_ads = (dGf_ads[2::] - dGf_ads[0:-2])/-10 # finite difference
    dSf_des = (dGf_des[2::] - dGf_des[0:-2])/-10
    temp_sub = target_temperatures[1:-1]

    # Enthalpies of folding
    # dG = dH - TdS
    # dH = dG + TdS
    dHf_ads = dGf_ads[1:-1] + temp_sub*dSf_ads
    dHf_des = dGf_des[1:-1] + temp_sub*dSf_des

    # Uncertainties - just propagation of error for now
    ddGf_ads = [numpy.sqrt(x[1,0]) for x in ddG]
    ddGf_ads = numpy.array(ddGf_ads)
    ddSf_ads = (ddGf_ads[0:-2]**2 + ddGf_ads[2::]**2)**.5/10
    ddHf_ads = (ddGf_ads[1:-1]**2 + (temp_sub*ddSf_ads)**2)**.5

    ddGf_des = [numpy.sqrt(x[3,2]) for x in ddG]
    ddGf_des = numpy.array(ddGf_des)
    ddSf_des = (ddGf_des[0:-2]**2 + ddGf_des[2::]**2)**.5/10
    ddHf_des = (ddGf_des[1:-1]**2 + (temp_sub*ddSf_des)**2)**.5
    
    # plotting
    plot(1,target_temperatures,dGf_ads,ddGf_ads,'temperature (K)','dG_folding (kcal/mol)','adsorbed, lam = %s' % lam, color)
    plot(2,target_temperatures,dGf_des,ddGf_des,'temperature (K)','dG_folding (kcal/mol)','desorbed, lam = %s' % lam, color)
    plot(3,temp_sub,dSf_ads,ddSf_ads,'temperature (K)','dS_folding (kcal/molK)','adsorbed, lam = %s' % lam, color)
    plot(4,temp_sub,dSf_des,ddSf_des,'temperature (K)','dS_folding (kcal/molK)','desorbed, lam = %s' % lam, color)
    plot(5,temp_sub,dHf_ads,ddHf_ads,'temperature (K)','dH_folding (kcal/mol)','adsorbed, lam = %s' % lam, color)
    plot(6,temp_sub,dHf_des,ddHf_des,'temperature (K)','dH_folding (kcal/mol)','desorbed, lam = %s' % lam, color)


if __name__ == '__main__':
    options = parse_args()
    main(options.file, options.lam,"k")
    plt.show()
