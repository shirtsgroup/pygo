#!/usr/bin/python2.4

# Ellen Zhong
# ellen.zhong@virginia.edu
# 10/06/2013

# This script is used to generate Figure 1 in the HMC paper: 
# 4 most average Q trajectories for each HMC, HMD, MC, and MD
# and give the folding rates for Table 1

import numpy
from optparse import OptionParser
import pdb
import matplotlib.pyplot as plt
import matplotlib
import foldingrate_Fig1

def plot(Q):
    matplotlib.rc('text', usetex = True)
    matplotlib.rc('font', family = 'serif')
    font = {'family' : 'serif',
            'size'   : 'larger'}
    f = plt.figure(1,(5,6))
    for i in range(len(Q)):
        plt.subplot(12,2,i+1) # 2 by 2 subplot
        plt.plot(Q[i,:])
        plt.ylim(0,1)
    f.text(.5,.04,'simulation length', ha='center',va='center',fontdict=font)
    f.text(.06,.5,'Q, fraction native',ha='center',va='center',rotation='vertical',fontdict=font)
    #f.text(.5,.96, title, ha='center',va='center',fontdict=font)
    #f.subplots_adjust(hspace=0)
    #f.subplots_adjust(vspace=0)
    plt.setp([a.get_xticklabels() for a in f.axes], visible=False)
    plt.setp([a.get_yticklabels() for a in f.axes], visible=False)

def main():
    print '---------HMC--------'
    HMC_file = '/home/edz3fz/proteinmontecarlo/results/1PGB/comparison/HMC/Qtraj_singleprot.npy'
    HMC_Q, HMC_count, HMC_count_smooth = foldingrate_Fig1.get_rate(HMC_file)
    plot(HMC_Q)
    plt.savefig('/home/edz3fz/proteinmontecarlo/manuscripts/figures/Fig2_HMC.eps')
    plt.show()

if __name__ == '__main__':
    main()

