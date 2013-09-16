#!/usr/bin/python2.4

import numpy
import matplotlib.pyplot as plt
import pickle
import optparse

def parse_args():
	parser = optparse.OptionParser(description='Plots the PMF(Q,z)')
	parser.add_option('-f','--file', dest='file', help = 'filename of pickled results (pmf_temp.pkl)')
	parser.add_option('-l','--lam', type='float', dest='lam', help = 'lambda value of surface attractiveness')
	(options,args) = parser.parse_args()
	return options

def plot_2D_PMF(temp,lam,bin_centers,f):
	X = [data[0] for data in bin_centers]
	Y = [data[1] for data in bin_centers]
	Z = f 
	
	plt.hexbin(X,Y,C=Z, gridsize=20)
	plt.xlabel('COM distance from surface')
	plt.ylabel('Q fraction native contacts')
	plt.title('temperature = %i, lambda = %3.1f' % (temp, lam))

def main():
    options = parse_args()
    results_file = options.file 	
    f = open(results_file,'rb')
    temp = pickle.load(f)
    bin_centers = pickle.load(f)
    f_i = pickle.load(f)
    df_i = pickle.load(f)
    f.close()

    plot_2D_PMF(temp,options.lam,bin_centers,f_i)
    file = options.file[0:options.file.rfind('/')]
    plt.savefig('%s/%i.png' % (file,temp))	
    plt.show()

if __name__ == '__main__':
    main()
