#!/usr/bin/python2.4

import numpy
import matplotlib.pyplot as plt
import pickle
import optparse
import pdb

def plot_2D_PMF(subplot,temp,lam,bin_centers,f):
    plt.subplot(subplot)
    X = [data[0] for data in bin_centers]
    Y = [data[1] for data in bin_centers]
    Z = f 
	
    plt.hexbin(X,Y,C=Z, gridsize=20)
#    plt.xlabel('COM distance from surface')
#    plt.ylabel('Q fraction native contacts')
    plt.title('temperature = %i, lambda = %3.1f' % (temp, lam),fontsize=6)

def main():
    files = []
    lam = numpy.array([.1,.35,.6])
    for l in lam:
        for t in numpy.array([300, 325, 350]):
            files.append('/home/edz3fz/proteinmontecarlo/results/1PGB/surface/umbrella_lambda%s/pmf_%i.pkl' % (str(l)[1::],  t))
    for i,file in enumerate(files): 
        f = open(file,'rb')
        temp = pickle.load(f)
        bin_centers = pickle.load(f)
        f_i = pickle.load(f)
        df_i = pickle.load(f)
        f.close()
        plot_2D_PMF('33%i' %(i+1),temp,lam[i/3],bin_centers,f_i)
    #file = options.file[0:options.file.rfind('/')]
    #plt.savefig('%s/%i.png' % (file,temp))	
    plt.show()

if __name__ == '__main__':
    main()
