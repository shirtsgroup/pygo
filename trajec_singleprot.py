from datetime import datetime
import numpy
#from optparse import OptionParser
import pdb

#parser=OptionParser()
##parser.add_option("-f", "--files", dest="datafile", default="GO_1PGB.pdb", help="Qtraj_singleprot.txt file")
#parser.add_option("-n", "--moves", dest="average", type="int", default='10', help="number of data points to average")
#(options,args)=parser.parse_args()
t1 = datetime.now()

def getmodel(ifile, num, ofile):
	fi = file(ifile,'r')
	line = ''
	while(line!='MODEL        '+str(num)+'\r\n'):
		line = fi.readline()
	w = line
	while(line!='ENDMDL\r\n'):
		line = fi.readline()
		w += line
	fi.close()
	fo = file(ofile,'a')
	fo.write(w)
	fo.close()
	

swap = 1000
numrep = 8
trange = [300., 450.]

T = numpy.empty(numrep)
alpha = (trange[1] / trange[0])**(1 / float(numrep - 1))
T[0] = trange[0]
for i in range(1, numrep):
        T[i] = T[i-1] * alpha

id = '91'
direc = 'replicaexchange/simlog'+id+'/'

prot_loc = numpy.loadtxt(direc+'protein_location.txt')
ifiles = [direc+'trajectory'+str(int(T[i]))+'.pdb' for i in range(numrep)]
ofiles = [direc+'singleprot_trajectory'+str(i)+'.pdb' for i in range(numrep)]

# write Model 0
for i in range(numrep):
	getmodel(ifiles[i],0,ofiles[i])

for i in xrange(len(prot_loc[0])-1):
	for j in range(numrep):
		rep = int(prot_loc[j][i])
		for k in range(swap*i+1,swap*(i+1)+1):
			getmodel(ifiles[rep],k,ofiles[j])
			print k
	
t2=datetime.now()
print t2-t1

