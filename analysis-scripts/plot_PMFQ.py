#!/usr/bin/python2.4

import numpy
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib as mpl
import cPickle
import optparse
import pdb
import plot_dG_solution

def parse_args():
	parser = optparse.OptionParser(description='Plots dG, dS, dH of folding and adsorption')
	parser.add_option('-f','--file', dest='file', default='', help = 'filename of pickled results')
	parser.add_option('-t','--temp', dest='temp', default='', help = 'temperature')
	(options,args) = parser.parse_args()
	return options

def pmf_t(temp):
    lam = [.1,.15,.2,.25,.3,.35,.4,.45,.5,.55]
    colors = cm.cool(numpy.linspace(0,1,len(lam)))
    plt.figure(8)
    plt.rc('text',usetex=True)
    direc = '/home/edz3fz/proteinmontecarlo/results/1PGB/surface'
    for i,l in enumerate(lam):
        f = open('%s/lambda%s/pmf_%i.pkl' %(direc,str(l)[1::],int(temp)),'rb')
        temp = cPickle.load(f)
        bin_centers = cPickle.load(f)
        f_i = cPickle.load(f)
        df_i = cPickle.load(f)
        f.close()
        idx = numpy.argmin(f_i[int(.7*len(f_i))::]) + int(.7*len(f_i))
        f_i -= f_i[idx]
        plt.errorbar(bin_centers[1::],f_i[1::],numpy.array(df_i[0,:])[0,1::],label=r'$\lambda$ = %s' % l, color=colors[i])
    
    f = open('/home/edz3fz/proteinmontecarlo/results/1PGB/solution/pmf_%i.pkl' %(int(temp)),'rb')
    temp = cPickle.load(f)
    bin_centers = cPickle.load(f)
    f_i = cPickle.load(f)
    df_i = cPickle.load(f)
    f.close()
    idx = numpy.argmin(f_i[int(.7*len(f_i))::]) + int(.7*len(f_i))
    f_i -= f_i[idx]
    plt.errorbar(bin_centers[1::],f_i[1::],numpy.array(df_i[0,:])[0,1::],label='solution', color='k')
    
    plt.xlabel('Q')
    plt.ylabel('PMF (kcal/mol)')
    plt.title('T = %s K' % temp)
    plt.legend(loc=0, prop={'size':8})


if __name__=='__main__':
    options = parse_args()
    if options.file:
        plot_dG_solution.solution_pmf(options.file)
        plt.show()
    if options.temp:
        pmf_t(options.temp)
        plt.show() 
