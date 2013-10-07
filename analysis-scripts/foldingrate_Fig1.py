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

def get_rate(file):
    smoothing = 50
    N_max = 120000 # make sure amount of data is the same
    Q = numpy.load(file)[0:N_max]
    Q -= 0.5 # used to count transitions
    # smooth data
    Q_smooth = numpy.zeros((len(Q),len(Q[0,:])/smoothing))
    for rep in range(len(Q_smooth)):
        for i in xrange(len(Q_smooth[0,:])):
		    Q_smooth[rep,i] = numpy.average(Q[rep,smoothing*i:(smoothing*i+smoothing)])
    # count transitions 
    count = numpy.zeros(len(Q))
    count_smooth = numpy.zeros(len(Q))
    for rep in range(len(Q)):
        # count number of transitions for raw data
        for i in xrange(len(Q[0,:])-1):
            a = Q[rep,i]
            b = Q[rep,i+1]
            if a*b < 0:
                count[rep] += 1
        # count number of transitions for smoothed data
        for i in xrange(len(Q_smooth[0,:])-1):
            a = Q_smooth[rep,i]
            b = Q_smooth[rep,i+1]
            if a*b < 0:
                count_smooth[rep] += 1
    return Q_smooth + 0.5, count, count_smooth

def get_4_index(count):
    avg = numpy.average(count)
    idx = numpy.argsort(numpy.abs(count-avg))
    return idx[0:4]

def plot(fig,Q,idx,title):
    matplotlib.rc('text', usetex = True)
    matplotlib.rc('font', family = 'serif')
    font = {'family' : 'serif',
            'size'   : 'larger'}
    f = plt.figure(fig,(6,4))
    for i,index in enumerate(idx):
        plt.subplot(2,2,i+1) # 2 by 2 subplot
        plt.plot(Q[index,:])
        plt.ylim(0,1)
    f.text(.5,.04,'simulation length', ha='center',va='center',fontdict=font)
    f.text(.06,.5,'Q, fraction native',ha='center',va='center',rotation='vertical',fontdict=font)
    f.text(.5,.96, title, ha='center',va='center',fontdict=font)
    #f.subplots_adjust(hspace=0)
    #f.subplots_adjust(vspace=0)
    plt.setp([a.get_xticklabels() for a in f.axes], visible=False)
#    plt.setp([a.get_yticklabels() for a in f.axes[::2]], visible=False)


def main():
    print '---------HMD--------'
    HMD_file = '/home/edz3fz/proteinmontecarlo/results/1PGB/comparison/HMD/Qtraj_singleprot.npy'
    HMD_Q, HMD_count, HMD_count_smooth = get_rate(HMD_file)
    print 'Number of transitions for raw data:'
    print HMD_count
    print numpy.average(HMD_count)
    print 'Number of transitions for smoothed data:'
    print HMD_count_smooth
    print numpy.average(HMD_count_smooth)
    HMD_index = get_4_index(HMD_count_smooth)
    plot(1,HMD_Q,HMD_index,'HMD')
    plt.savefig('/home/edz3fz/proteinmontecarlo/manuscripts/figures/Fig1a_HMD.eps')

    print '---------HMC--------'
    HMC_file = '/home/edz3fz/proteinmontecarlo/results/1PGB/comparison/HMC/Qtraj_singleprot.npy'
    HMC_Q, HMC_count, HMC_count_smooth = get_rate(HMC_file)
    print 'Number of transitions for raw data:'
    print HMC_count
    print numpy.average(HMC_count)
    print 'Number of transitions for smoothed data:'
    print HMC_count_smooth
    print numpy.average(HMC_count_smooth)
    HMC_index = get_4_index(HMC_count_smooth)
    plot(2,HMC_Q,HMC_index,'HMC')
    plt.savefig('/home/edz3fz/proteinmontecarlo/manuscripts/figures/Fig1b_HMC.eps')
    
    print '---------MD---------'
    MD_file = '/home/edz3fz/proteinmontecarlo/results/1PGB/comparison/MD/Qtraj_singleprot.npy'
    MD_Q, MD_count, MD_count_smooth = get_rate(MD_file)
    print 'Number of transitions for raw data:'
    print MD_count
    print numpy.average(MD_count)
    print 'Number of transitions for smoothed data:'
    print MD_count_smooth
    print numpy.average(MD_count_smooth)
    MD_index = get_4_index(MD_count_smooth)
    plot(3,MD_Q,MD_index,'MD')
    plt.savefig('/home/edz3fz/proteinmontecarlo/manuscripts/figures/Fig1c_MD.eps')
    
    print '---------MC---------'
    MC_file = '/home/edz3fz/proteinmontecarlo/results/1PGB/comparison/MC/Qtraj_singleprot.npy'
    MC_Q, MC_count, MC_count_smooth = get_rate(MC_file)
    print 'Number of transitions for raw data:'
    print MC_count
    print numpy.average(MC_count)
    print 'Number of transitions for smoothed data:'
    print MC_count_smooth
    print numpy.average(MC_count_smooth)
    MC_index = get_4_index(MC_count_smooth)
    plot(4,MC_Q,MC_index,'MC')
    plt.savefig('/home/edz3fz/proteinmontecarlo/manuscripts/figures/Fig1d_MC.eps')

    plt.show()

if __name__ == '__main__':
    main()

