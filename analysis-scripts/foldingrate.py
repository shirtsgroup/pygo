#!/usr/bin/python2.4


import numpy
from optparse import OptionParser
import pdb
import matplotlib.pyplot as plt
import matplotlib

font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : 22}

parser=OptionParser()
parser.add_option("--direc", dest="datafile", help="Qtraj_singleprot.txt file location")
parser.add_option("-n", "--smooth-traj", dest="average", type="int", default='10', help="number of data points to average")
parser.add_option("-s", "--save-interval", dest="save", type="int", default='1000', help="number of moves between save operations")
#parser.add_option("-k", "--swapsize", dest="swap", type="int", default='1000', help="number of moves between swap operations")
parser.add_option("-r", "--replicas", default=24, type="int",dest="replicas", help="number of replicas (default=24)")
parser.add_option("--show", default=False, action="store_true", dest="show", help="show plot (default: False)")
parser.add_option("-o", "--output", default='', dest="ofile", help="Plot filename to save to")
(options,args)=parser.parse_args()


file = options.datafile+'/Qtraj_singleprot.npy'
n = options.average

Q = numpy.load(file)
Q = Q - .5
numrep = options.replicas
save = options.save
#swap = options.swap
show = options.show
ofile = options.ofile

fig =plt.figure(1)
#plt.xlabel('move/%i'% (n*save))
#plt.ylabel('Q fraction native')
Q_avg = numpy.zeros((numrep,len(Q[0,:])/n))
for rep in range(numrep):
	for i in range(len(Q[0,:])/n): # Average every 5 points
		Q_avg[rep,i] = numpy.average(Q[rep,n*i:(n*i+n)])
	plt.subplot(numrep/2,2,rep+1)
	plt.plot(numpy.arange(len(Q_avg[0,:])),Q_avg[rep,:]+.5)
	plt.ylim(0,1)
fig.text(.5,.04, 'move/%i' % (n*save), ha='center',va='center', fontdict=font)
fig.text(.5,.96, 'Q trajectories for single proteins', ha='center',va='center',fontdict=font)
fig.text(.06,.5,'Q fraction native',ha='center',va='center',rotation='vertical',fontdict=font)


count = numpy.zeros(numrep)
for rep in range(numrep):
	for i in xrange(len(Q[0,:])-1):
		a = Q[rep,i]
		b = Q[rep,i+1]
		if a*b < 0: # sign change
			count[rep] += 1

count_new = numpy.zeros(numrep)
for rep in range(numrep):
	for i in xrange(len(Q_avg[0,:])-1):
		a = Q_avg[rep,i]
		b = Q_avg[rep,i+1]
		if a*b < 0: # sign change
			count_new[rep] += 1
print 'Original data'
print count
print numpy.average(count)

print 'Average every '+str(n)+' data points'
print count_new
print numpy.average(count_new)

if ofile:
	plt.savefig(ofile+'.png')
	plt.savefig(ofile+'.pdf')
	plt.savefig(ofile+'.eps')

orgdataplt = False
if orgdataplt:
	plt.figure(2)
	import cPickle
	pt = cPickle.load(open(options.datafile+'protein_location.pkl','rb'))
	pt = numpy.array(pt)
#	pt = numpy.loadtxt(options.datafile+'protein_location.txt')
	pt /=numrep
	for i in range(numrep):
		plt.subplot(numrep/2,2,i+1)
		plt.plot(numpy.arange(len(Q[0,:])),Q[i,:]+.5)
#		plt.plot(numpy.arange(len(Q[0,:])),Q[i,:]+.5,'o')
		plt.plot(numpy.arange(0,len(Q[0,:]),swap/save),pt[i])
	plt.xlabel('moves/%i' % save)
	plt.ylabel('Q fraction native')
	plt.title('Q trajectories for single proteins')
	plt.savefig('%s/Qtraj_singleprot.png' % options.datafile)
if show:
	plt.show()

