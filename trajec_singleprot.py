from datetime import datetime
import numpy
#from optparse import OptionParser
import pdb
from sys import stdout
#parser=OptionParser()
##parser.add_option("-f", "--files", dest="datafile", default="GO_1PGB.pdb", help="Qtraj_singleprot.txt file")
#parser.add_option("-n", "--moves", dest="average", type="int", default='10', help="number of data points to average")
#(options,args)=parser.parse_args()
t1 = datetime.now()

def getmodel(fi, num, ofile):
	line = ''
	while(line!='MODEL        '+str(num)+'\r\n'):
		line = fi.readline()
	w = line
	while(line!='ENDMDL\r\n'):
		line = fi.readline()
		w += line
	fo = file(ofile,'a')
	fo.write(w)
	fo.close()
	

swap = 40
numrep = 8
trange = [300., 450.]

T = numpy.empty(numrep)
alpha = (trange[1] / trange[0])**(1 / float(numrep - 1))
T[0] = trange[0]
for i in range(1, numrep):
        T[i] = T[i-1] * alpha

id = '92'
direc = 'replicaexchange/simlog'+id+'/'

prot_loc = numpy.loadtxt(direc+'protein_location.txt')
ifiles = [direc+'trajectory'+str(int(T[i]))+'.pdb' for i in range(numrep)]
ofiles = [direc+'singleprot_trajectory'+str(i)+'.pdb' for i in range(numrep)]
ifiles_obj = [file(ifiles[i],'r') for i in range(numrep)]
# write Model 0
for i in range(numrep):
	getmodel(ifiles_obj[i],0,ofiles[i])

for i in xrange(len(prot_loc[0])-1):
	for j in range(numrep):
		rep = int(prot_loc[j][i])
		for k in range(swap*i+1,swap*(i+1)+1):
			getmodel(ifiles_obj[rep],k,ofiles[j])
			stdout.write(str(k)+'\r')
			stdout.flush()
	
t2=datetime.now()
print t2-t1

