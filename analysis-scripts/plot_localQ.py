#!/usr/bin/python2.4

# Ellen Zhong
# ellen.zhong@virginia.edu
# 10/10/2013

import numpy
from optparse import OptionParser
import pdb
import matplotlib.pyplot as plt
import matplotlib
import localQ
import timeseries

def get_peptides(contacts,Q):
    peptides = [[0,10],[0,18],[0,23],[11,17],[0,27],[28,53],[41,53],[72,80],[82,92],[82,94],[87,94],[95,103],[95,104],[102,113],[122,132],[122,147],[130,148],[131,148],[133,148],[133,155]]
    assert(len(peptides)==20)
    Q_pep = numpy.zeros(len(peptides))
    for i,pep in enumerate(peptides):
        ncontact = 0
        for j in range(pep[0],pep[1]+1):
            Q_pep[i] += len(contacts[j])*Q[j]
            ncontact += len(contacts[j])
        Q_pep[i] /= float(ncontact)
    return Q_pep

def subsample(Q_n,localQ):
    print 'Subsampling the data'
    g = timeseries.statisticalInefficiency(Q_n)
    indices = numpy.array(timeseries.subsampleCorrelatedData(Q_n,g))
    print '%i uncorrelated samples found of %i original samples' %(len(indices),len(Q_n))
    localQ = localQ[:,indices]
    return localQ
        
def main():
    ID = [225,227,230,232,235]
    ID = [232,238]
#    ID = [232,235]
    T = [274]
    files = ['/home/edz3fz/proteinmontecarlo/results/1BSQ/simlog%s/localQ%s.npy' %(id,t) for id in ID for t in T]
    files = ['/home/edz3fz/proteinmontecarlo/results/1BSQ/simlog%s/localQ%s.npy' %(id,t) for t in T for id in ID]
    
    print files

    colors = matplotlib.cm.cool(numpy.linspace(0,1,len(ID)))
    
    numbeads = 162
    paramfile = '/home/edz3fz/proteinmontecarlo/GO_1BSQ.param'
    contacts, sig_cutoff = localQ.get_contacts(paramfile,numbeads)

    #ID = ['solution','surface']
    i=0
    for j,t in enumerate(T):
        plt.figure(j+1)
        plt.title('T=%s' % t)
        for k,id in enumerate(ID):
            data = numpy.load(files[i])
            #data = data[:,-20000::]
            #data = data[:,::5]

            correlated = True
            if correlated:
                Qfile = '%s/fractionnative%i.npy' %(files[i][0:files[i].rfind('/')],t)
                Q_n = numpy.load(Qfile)
                Q_n = Q_n[0:len(data[0,:])]
                data = subsample(Q_n,data)
            data = numpy.average(data,axis=1)
            pep = get_peptides(contacts,data)
            #plt.plot(data,'o',color=colors[k])
            plt.plot(1-pep,color=colors[k],label=str(id))
            plt.plot(1-pep,'o',color=colors[k],label=str(id))
            numpy.savetxt('%s.txt' %id, pep)
            i += 1
        plt.legend(loc=3)
    plt.show() 

if __name__ == '__main__':
    main()
