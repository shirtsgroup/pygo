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
parser.add_option("--direc", dest="datafile", default="replicaexchange/simlog0/", help="Qtraj_singleprot.txt file location")
parser.add_option("-n", "--moves", dest="average", type="int", default='10', help="number of data points to average")
#parser.add_option("-t", "--temprange", nargs=2, default=[300.0,450.0], type="float", dest="temprange", help="temperature range of replicas")
parser.add_option("-s", "--stepsize", dest="step", type="int", default='1000', help="number of moves between save operations")
parser.add_option("-k", "--swapsize", dest="swap", type="int", default='1000', help="number of moves between swap operations")
parser.add_option("-r", "--replicas", default=8, type="int",dest="replicas", help="number of replicas")
(options,args)=parser.parse_args()


file = options.datafile+'/Qtraj_singleprot.npy'
n = options.average

#file = '/home/edz3fz/proteinmontecarlo/replicaexchange/simlog24/Qtraj_singleprot.txt'
Q_full = numpy.load(file)

HMD = [1,2,16,20]

Q = numpy.zeros((4,len(Q_full[0,:])))
for i in range(4):
	Q[i,:] = Q_full[HMD[i],:]

Q = Q - .5
numrep = options.replicas
numrep = 4
save = options.step
swap = options.swap
#Q = Q[:,::10]
# Averages every n points to reduce noise
#n = 5
fig =plt.figure(1)
#plt.axes(frameon=False)
#plt.xticks([])
#plt.yticks([])
plt.xlabel('move/%i'% (n*save))
plt.ylabel('Q fraction native')
Q_new = numpy.zeros((numrep,len(Q[0,:])/n))
for rep in range(numrep):
	for i in range(len(Q[0,:])/n): # Average every 5 points
		Q_new[rep,i] = numpy.average(Q[rep,n*i:(n*i+n)])
	plt.subplot(numrep/2,2,rep+1)
	plt.plot(numpy.arange(len(Q_new[0,:])),Q_new[rep,:]+.5)
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
	for i in xrange(len(Q_new[0,:])-1):
		a = Q_new[rep,i]
		b = Q_new[rep,i+1]
		if a*b < 0: # sign change
			count_new[rep] += 1
print 'Original data'
print count
print numpy.average(count)

print 'Average every '+str(n)+' data points'
print count_new
print numpy.average(count_new)
pdb.set_trace()
plt.savefig(options.datafile+'/Qtraj_singleprot_avg.png')

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
plt.show()

