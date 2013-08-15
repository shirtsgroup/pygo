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
(options,args)=parser.parse_args()


file = options.datafile+'/Qtraj_singleprot.npy'
n = options.average

Q = numpy.load(file)
numrep = options.replicas
save = options.save
#swap = options.swap
show = options.show


fig =plt.figure(1)
#plt.xlabel('move/%i'% (n*save))
#plt.ylabel('Q fraction native')
Q_avg = numpy.zeros((numrep,len(Q[0,:])/n))
for rep in range(numrep):
	for i in range(len(Q[0,:])/n): # Average every 5 points
		Q_avg[rep,i] = numpy.average(Q[rep,n*i:(n*i+n)])
	plt.subplot(numrep/2,2,rep+1)
	plt.plot(numpy.arange(len(Q_avg[0,:])),Q_avg[rep,:])
	plt.ylim(0,1)
fig.text(.5,.04, 'move/%i' % (n*save), ha='center',va='center', fontdict=font)
fig.text(.5,.96, 'Q trajectories for single proteins', ha='center',va='center',fontdict=font)
fig.text(.06,.5,'Q fraction native',ha='center',va='center',rotation='vertical',fontdict=font)


mixing_time = [[] for x in range(numrep)]
for rep in range(numrep):
	highQ_index = 0
	lowQ_index = 0
	for i in xrange(len(Q[0,:])-1):
		if Q[rep,i] > 0.95:
			if lowQ_index > highQ_index:
				mixing_time[rep].append(i-highQ_index)
			highQ_index = i
		if Q[rep,i] < 0.2:
			lowQ_index = i


mixing_time_avg = [[] for x in range(numrep)]
for rep in range(numrep):
	highQ_index = 0
	lowQ_index = 0
	for i in xrange(len(Q_avg[0,:])-1):
		if Q_avg[rep,i] > 0.95:
			if lowQ_index > highQ_index:
				mixing_time_avg[rep].append(n*(i-highQ_index))
			highQ_index = i
		if Q_avg[rep,i] < 0.2:
			lowQ_index = i

print 'Original data'
print mixing_time
print [numpy.average(x) for x in mixing_time]
print numpy.average([numpy.average(x) for x in mixing_time if not numpy.isnan(numpy.average(x))])

print 'Average every '+str(n)+' data points'
print mixing_time_avg
print [numpy.average(x) for x in mixing_time_avg]
print numpy.average([numpy.average(x) for x in mixing_time_avg if not numpy.isnan(numpy.average(x))])

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

