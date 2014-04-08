#!/usr/bin/python2.4

# Ellen Zhong
# ellen.zhong@virginia.edu
# 10/06/2013

import numpy
from optparse import OptionParser
import pdb
import matplotlib.pyplot as plt
import matplotlib
import Fig2_foldingrate

def get_rate(file):
    smoothing = 50
#    N_max = 120000 # make sure amount of data is the same
    Q = numpy.load(file)#[0:N_max]
    Q -= 0.6 # used to count transitions
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
    return Q_smooth + 0.6, count, count_smooth

def plot(i,Q):
    matplotlib.rc('text', usetex = True)
    matplotlib.rc('font', family = 'serif')
    font = {'family' : 'serif',
            'size'   : 'larger'}
    f = plt.figure(i,(5,6))
    colors = matplotlib.cm.jet(numpy.linspace(0,1,len(Q)))
    for j in range(len(Q[0])):
        plt.subplot(12,2,j+1) # 2 by 2 subplot
        index = 0
        for i in range(len(Q)):
            length = len(Q[i][j,:])
            plt.plot(range(index,index+length),Q[i][j,:],color=colors[i])
            index += length
        plt.ylim(0,1)
    f.text(.5,.04,'simulation length', ha='center',va='center',fontdict=font)
    f.text(.06,.5,'Q, fraction native',ha='center',va='center',rotation='vertical',fontdict=font)
    #f.text(.5,.96, title, ha='center',va='center',fontdict=font)
    #f.subplots_adjust(hspace=0)
    #f.subplots_adjust(vspace=0)
    plt.setp([a.get_xticklabels() for a in f.axes], visible=False)
    plt.setp([a.get_yticklabels() for a in f.axes], visible=False)


def main():
    lam = [.1, .15, .2, .25, .3, .35, .4, .45, .5, .55, .6]
    files = ['/home/edz3fz/proteinmontecarlo/results/1PGB/surface/umbrella_lambda%s' % str(x)[1::] for x in lam]
    
    Z = numpy.arange(9,31.5,1.5)
    Z = numpy.concatenate((Z,numpy.array([33,36,39,42,45,48])))
    Z = numpy.array([9,10.5,12,13.5,15,16.5])

    
    colors = matplotlib.cm.jet(numpy.linspace(0,1,len(Z)))   
    counts = numpy.zeros(len(files))
    uncert = numpy.zeros(len(files))
    Q_all = [[] for x in range(len(files))]
    for i,direc in enumerate(files):
        print 'Analyzing', direc
        for z in Z:
    #        print 'Analyzing', z
            f = '%s/%i/Qtraj_singleprot.npy' % (direc,z)
            q, count, count_smooth = get_rate(f)
            Q_all[i].append(q)
            counts[i] += numpy.average(count_smooth)
            uncert[i] += numpy.std(count_smooth)**2
            #Q_all[i].append(numpy.zeros((24,24)))
            #print 'Number of transitions for raw data:', numpy.average(count)
            #print 'Number of transitions for smoothed data:', numpy.average(count_smooth)
        #plot(i,Q_all[i])
    print counts
    data = numpy.zeros((len(files),2))
    data[:,0] = counts
    data[:,1] = uncert**.5
    numpy.save('adsorbed_rate.out',data)

    Z = numpy.array([33,36,39,42,45,48])
    counts = numpy.zeros(len(files))
    uncert = numpy.zeros(len(files))
    Q_all = [[] for x in range(len(files))]
    for i,direc in enumerate(files):
        print 'Analyzing', direc
        for z in Z:
    #        print 'Analyzing', z
            f = '%s/%i/Qtraj_singleprot.npy' % (direc,z)
            q, count, count_smooth = get_rate(f)
            Q_all[i].append(q)
            counts[i] += numpy.average(count_smooth)
            uncert[i] += numpy.std(count_smooth)**2
            #Q_all[i].append(numpy.zeros((24,24)))
            #Q_all[i].append(numpy.zeros((24,24)))
            #print 'Number of transitions for raw data:', numpy.average(count)
            #print 'Number of transitions for smoothed data:', numpy.average(count_smooth)
        #plot(i,Q_all[i])
    print counts
    data = numpy.zeros((len(files),2))
    data[:,0] = counts
    data[:,1] = uncert**.5
    numpy.save('desorbed_rate.out',data)

   

    #plt.show()
if __name__ == '__main__':
    main()

