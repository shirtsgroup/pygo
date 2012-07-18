from numpy import *
import matplotlib.pyplot as plt

trange=[300.0,600.0]
numreplicas=16
T=empty(numreplicas)
alpha=(trange[1]/trange[0])**(1/float(numreplicas-1))
T[0]=trange[0]
for i in range(1,numreplicas):
	T[i]=T[i-1]*alpha
print T
files=[]
for i in range(len(T)):
	files.append('replicaexchange/replica'+str(i)+'/fractionnative'+str(int(T[i]))+'.txt')


nc=loadtxt(files[0])

for file in files:
	nctemp=loadtxt(file)
	nc=vstack((nc,nctemp))
nc=nc[1:numreplicas+1,:]
print shape(nc)
nc=nc[:,8000:-1]
ncavg=average(nc,axis=1)
print ncavg
#ncavg=ncavg/151.
plt.figure(1)
plt.plot(T,ncavg)
plt.plot(T,ncavg,'bo')
plt.xlabel('Temperature (Kelvin)')
plt.ylabel('Q, fractional nativeness')
plt.title('Go-like model MC simulation of 2QUG.pdb')
plt.savefig('ncplot1R69.png')
plt.show()

##histogram!

#Efiles=[]
#for i in range(len(T)):
	#Efiles.append('replicaexchange/replica'+str(i)+'/energy'+str(T[i])+'.txt')


#hist=loadtxt(Efiles[0])

#for file in Efiles:
	#histtemp=loadtxt(file)
	#hist=vstack((hist,histtemp))
#hist=hist[1:9,:]
##hist=hist[:,25000:35001]
#plt.figure(2)
#for i in range(numreplicas):
	##plt.figure(i)
	#plt.hist(hist[i,:],40,label=str(T[i])+' K')
#plt.xlabel('Energy (kcal/mol)')
#plt.legend(loc=2)
#plt.title('Go-like model MC simulation of 1R69.pdb')
#plt.savefig('histplot1R69.png')
#plt.show()
