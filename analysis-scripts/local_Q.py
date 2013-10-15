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
    parser.add_option("-r", "--replicas", default=32, type="int",dest="replicas", help="number of replicas (default: 32)")
    parser.add_option("-n", "--N_max", default=100000, type="int",dest="N_max", help="number of data points to read in (default: 100k)")
    parser.add_option("-s", "--skip", default=1, type="int",dest="skip", help="skip every n data points")
    parser.add_option("--direc", dest="direc", help="Qtraj_singleprot.txt file location")
    parser.add_option("--tfile", dest="tfile", default="/home/edz3fz/proteinmontecarlo/T.txt", help="file of temperatures (default: T.txt)")
    (options,args) = parser.parse_args()
    return options

def get_contacts(paramfile,numbeads):
    contacts = [[] for x in range(numbeads)]
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
            contacts[i-1].append(j-1) # 0 indexing
            contacts[j-1].append(i-1) 
    return contacts

def get_Rg(direc,K,N,N_max,T):
    Rg_kn = numpy.empty([K,N])
    print "Reading traj data..."
    for i,t in enumerate(T):
        print 'Reading traj at %i' % t
        file = open('%s/trajectory%i' %(direc,t),'r')
        for j in range(N_max - N): # skip through beginning
            _ = numpy.load(file)
        nbeads = len(_) # get numbeads
        for j in range(N):
            coord = numpy.load(file)
            Rg_kn[i,j] = Rg(coord,nbeads)
        file.close()
    return Rg_kn


def main():
    args = parse_args()
    T = numpy.loadtxt(args.tfile)
    
