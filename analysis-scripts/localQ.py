#!/usr/bin/python2.4

# Ellen Zhong
# ellen.zhong@virginia.edu
# 10/10/2013

import numpy
from optparse import OptionParser
import pdb
import matplotlib.pyplot as plt
import matplotlib

def parse_args():
    parser=OptionParser()
    parser.add_option('-n', '--numbeads', default=162, type='int',dest='numbeads', help='number of residues (default: 162)')
    parser.add_option('-f','--file', dest='file', help='trajectory file')
    parser.add_option('-p','--param', dest='paramfile', help='GO model parameter file')
    (options,args) = parser.parse_args()
    return options

def get_contacts(paramfile,numbeads):
    contacts = [[] for x in range(numbeads)]
    sig_cutoff = [[] for x in range(numbeads)]
    f = open(paramfile,'r')
    while 1:
        line = f.readline()
        if 'NBFIX' in line:
            break
    while 1:
        line = f.readline()
        if not line:
            break
        if 'G' in line:
            i = int(line[1:4])
            j = int(line[9:12])
            s = float(line[32:-1])
            s = 1.2*1.2*s**2 # cutoff is 1.2sig
            contacts[i-1].append(j-1) # 0 indexing
            contacts[j-1].append(i-1) 
            sig_cutoff[i-1].append(s)
            sig_cutoff[j-1].append(s)
    return contacts,sig_cutoff

def get_localQ(coord,contacts,sig_cutoff):
    Q = numpy.zeros(len(coord))
    for i,contact in enumerate(contacts):
        if contact: # if there are native contacts
            for j_ind,j in enumerate(contact):
                if numpy.sum((coord[i,:] - coord[j,:])**2) < sig_cutoff[i][j_ind]:
                    Q[i] += 1
            Q[i] /= float(len(contact))
    return Q
                

def get_traj_localQ(file,numbeads,N,contacts,sig_cutoff):
    localQ_kn = numpy.empty([numbeads,N])
    print 'Reading traj data for %s...' % file
    file = open(file,'r')
#    for i in range(N_max - N): # skip through beginning
#        _ = numpy.load(file)
#    assert(numbeads==len(_))
    for i in range(N):
        coord = numpy.load(file)
        localQ_kn[:,i] = get_localQ(coord,contacts,sig_cutoff)
    file.close()
    return localQ_kn


def main():
    args = parse_args()
    direc = args.file[0:args.file.rfind('/')]
    print direc
    contacts,sig_cutoff = get_contacts(args.paramfile,args.numbeads)
    N = 1000
    localQ_kn = get_traj_localQ(args.file,args.numbeads,N,contacts,sig_cutoff)
    numpy.save(direc+'/localQ.npy',localQ_kn)

def main_noargs():
    ID = [225,227,230,232,233,234,235]
    ID = [235,232]
    T = [271,274,277,280,283]
    T = [286,289,292]
    T = [260,265,268]
    files = ['/home/edz3fz/proteinmontecarlo/results/1BSQ/simlog%s/trajectory%s' %(id,t) for id in ID for t in T]
    T = 6*T
    numbeads = 162
    N = 60000
    paramfile = '/home/edz3fz/proteinmontecarlo/GO_1BSQ.param'

    contacts, sig_cutoff = get_contacts(paramfile,numbeads)
    for i,file in enumerate(files):
        localQ_kn = get_traj_localQ(file,numbeads,N,contacts,sig_cutoff)
        direc = file[0:file.rfind('/')]
        numpy.save('%s/localQ%s.npy' %(direc,T[i]),localQ_kn)

if __name__ == '__main__':
    main_noargs()
