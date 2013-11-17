#!/usr/bin/python2.4

import numpy
import writetopdb_docking as pdb
import optparse
from pdb import *
import energyfunc

def parse_args():
    parser = optparse.OptionParser(description='converts a npy trajectory into a pdb trajectory for VMD')
    parser.add_option('--pdb', dest='pdbfile', help='original pdb file')
    parser.add_option('--traj', dest='trajfile', default='', help='trajectory file to convert')
    (options,args) = parser.parse_args()
    return options

def get_seqres(file):
    f = open(file, 'r')
    missingres = []
    res = []
    while 1:
		line = f.readline()
		if not line:
			break
		if line.startswith('REMARK 465'):
			try:
				missingres.append(int(line[21:-1])-1)
			except: pass
		if line.startswith('SEQRES'):
			words = line.split(' ')
			words = [x for x in words if x]
			res.extend(words[4:-1])
    res = [res[i] for i in range(len(res)) if i not in missingres]
    f.close()
    return res	

def get_beta(seqres):
    index = [energyfunc.getindex(residue) for residue in seqres]
    scale = []
    hydropathy = numpy.loadtxt('/home/edz3fz/proteinmontecarlo/hydropathy.txt')
    for i in index:
        scale.append(hydropathy[i])
    return scale

def convert_pdb(args):
    _,pdb_text = pdb.get_coord(args.pdbfile)
    
    # add residue identity to pdb file
    seqres = get_seqres(args.pdbfile.replace('GO_',''))
    assert(len(seqres) == len(_))
    scale = get_beta(seqres)
    scale = numpy.array(scale)*100
    count = 0
    for i in range(len(pdb_text)):
        line = pdb_text[i]
        if 'ATOM' in line:
            n = len(line)
            extra_spaces = ' '*(60 - n)
            pdb_text[i] = '%s%s%s%5f%s' %(line[0:17],seqres[count],line[20:61],scale[count],line[66::])
            count += 1

    # convert npy trajectory to pdb trajectory
    f = open(args.trajfile,'r')
    print 'reading trajectory from %s' % args.trajfile
    print 'saving trajectory to %s.pdb' % args.trajfile
    i = 0
    while 1:
        try:
            coord = numpy.load(f)
            pdb.write_coord(args.trajfile+'.pdb',i,coord,pdb_text)
            i += 1
        except:
            break
    print 'trajectory converted'
    print 'done'

def main():
    args = parse_args()
    convert_pdb(args)

if __name__ == '__main__':
    main()    
