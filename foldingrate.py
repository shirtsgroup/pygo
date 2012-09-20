import numpy

file = '/home/edz3fz/proteinmontecarlo/replicaexchange/simlog24/Qtraj_singleprot.txt'
Q = numpy.loadtxt(file)
Q = Q - .5
numrep = 8

count = numpy.zeros(numrep)
for rep in range(numrep):
	for i in xrange(len(Q[0,:])-1):
		a = Q[rep,i]
		b = Q[rep,i+1]
		if a*b < 0: # sign change
			count[rep] += 1

print count
print numpy.average(count)	
