#!/usr/bin/python2.4

import numpy
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import cPickle
import optparse

def parse_args():
	parser = optparse.OptionParser(description='Plots dG, dS, dH of adsorption and adsorption')
	parser.add_option('-f','--file', dest='file', help = 'filename of pickled results')
	parser.add_option('-l','--lam', type='float', dest='lam', help = 'lambda value of surface attractiveness')
	(options,args) = parser.parse_args()
	return options

def plot(fig,x,y,dy,xlabel,ylabel,label,color):
    plt.figure(fig)
    plt.errorbar(x,y,dy,label=label, color=color)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)

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
    
    # plotting
    plot(1,target_temperatures,dGads_folded,ddGads_folded,'temperature (K)','dG_adsorption (kcal/mol)','lam = %s' % lam, color)
    plot(2,target_temperatures,dGads_unfolded,ddGads_unfolded,'temperature (K)','dG_adsorption (kcal/mol)','lam = %s' % lam, color)
    plot(3,temp_sub,dSads_folded,ddSads_folded,'temperature (K)','dS_adsorption (kcal/molK)','lam = %s' % lam, color)
    plot(4,temp_sub,dSads_unfolded,ddSads_unfolded,'temperature (K)','dS_adsorption (kcal/molK)','lam = %s' % lam, color)
    plot(5,temp_sub,dHads_folded,ddHads_folded,'temperature (K)','dH_adsorption (kcal/mol)','lam = %s' % lam, color)
    plot(6,temp_sub,dHads_unfolded,ddHads_unfolded,'temperature (K)','dH_adsorption (kcal/mol)','lam = %s' % lam, color)

    plt.figure(1)
    plt.title('Free energy of adsorption, folded')
    plt.figure(2)
    plt.title('Free energy of adsorption, unfolded')
    plt.figure(3)
    plt.title('Entropy of adsorption, folded')
    plt.figure(4)
    plt.title('Entropy of adsorption, unfolded')
    plt.figure(5)
    plt.title('Enthalpy of adsorption, folded')
    plt.figure(6)
    plt.title('Enthalpy of adsorption, unfolded')

if __name__ == '__main__':
    options = parse_args()
    main(options.file, options.lam,"k")
    plt.show()
