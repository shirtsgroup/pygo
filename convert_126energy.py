#!/usr/bin/python2.4

import numpy
import writetopdb_docking as pdb
import optparse
import energyfunc
from pdb import *

def parse_args():
    parser = optparse.OptionParser(description='Calculates prot-surf energy if energy function were 12-6')
    parser.add_option('--surf', dest='surffile', help='surface pdb file')
    parser.add_option('--pdb', dest='pdbfile', help='original protein pdb file')
    parser.add_option('--traj', dest='trajfile', default='', help='trajectory file of protein')
    parser.add_option('-n', dest='n', type='int',default=100, help='number of frames to re-analyze')
    (options,args) = parser.parse_args()
    return options

def convert_pdb(args):
    _,pdb_text = pdb.get_coord(args.pdbfile)
    f = open(args.trajfile,'r')
    print 'reading trajectory from %s' % args.trajfile
    print 'saving trajectory to %s.pdb' % args.trajfile
    i = 0
    while 1:
        try:
            coord = numpy.load(f)
            pdb.write_coord(args.trajfile+'pdb',i,coord,pdb_text)
            i += 1
        except:
            break
    print 'trajectory converted'
    print 'done'

def main():
    args = parse_args()
    surf,_ = pdb.get_coord(args.surffile)
    numbeads = 162
    scale = 1
    param = energyfunc.getsurfparam_old(args.pdbfile, numbeads, len(surf), numbeads*len(surf), scale)
    energy = numpy.zeros((args.n,2))
    f = open(args.trajfile,'r')
    for i in range(args.n):
        coord = numpy.load(f)
        energy[i,:] = energyfunc.csurfenergy_old(coord, surf, numbeads, numbeads*len(surf), param, scale)
    f.close()
    print energy
    file = '%s/surfenergy%i_126.npy' % (args.trajfile[0:args.trajfile.rfind('/')], int(args.trajfile[-7:-4]))
    numpy.save(file,energy)


if __name__ == '__main__':
    main()    
