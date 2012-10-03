import numpy
from optparse import OptionParser
import pdb
import matplotlib.pyplot as plt

parser=OptionParser()
parser.add_option("-f", "--files", dest="datafile", default="GO_1PGB.pdb", help="Qtraj_singleprot.txt file")
parser.add_option("-n", "--moves", dest="average", type="int", default='10', help="number of data points to average")
(options,args)=parser.parse_args()


file = options.datafile+'/Qtraj_singleprot.txt'
n = options.average

#file = '/home/edz3fz/proteinmontecarlo/replicaexchange/simlog24/Qtraj_singleprot.txt'
Q = numpy.loadtxt(file)
Q = Q - .5
numrep = 8

# Averages every n points to reduce noise
#n = 5
plt.figure(1)
Q_new = numpy.zeros((numrep,len(Q[0,:])/n))
for rep in range(numrep):
	for i in range(len(Q[0,:])/n): # Average every 5 points
		Q_new[rep,i] = numpy.average(Q[rep,n*i:(n*i+n)])
	plt.subplot(numrep/2,2,rep+1)
	plt.plot(numpy.arange(len(Q_new[0,:])),Q_new[rep,:])
plt.xlabel('moves/step')
plt.ylabel('Q fraction native')
plt.title('Q trajectories for single proteins')



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

plt.savefig(options.datafile+'/Qtraj_singleprot_avg.png')

orgdataplt = True
if orgdataplt:
	plt.figure(2)
	for i in range(numrep):
		plt.subplot(numrep/2,2,i+1)
		plt.plot(numpy.arange(len(Q[0,:])),Q[i,:])
	plt.xlabel('moves/step')
	plt.ylabel('Q fraction native')
	plt.title('Q trajectories for single proteins')
	plt.savefig('%s/Qtraj_singleprot.png' % options.datafile)
plt.show()

