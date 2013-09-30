from numpy import *
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from optparse import OptionParser

parser=OptionParser()
parser.add_option("-t", "--temprange", nargs=2, default=[300.0,450.0], type="float", dest="temprange", help="temperature range of replicas")
parser.add_option("-s", "--stepsize", dest="step", type="int", default='1000', help="number of moves between save operations")
parser.add_option("-r", "--replicas", default=8, type="int",dest="replicas", help="number of replicas")
parser.add_option("-a", "--energy", action="store_true", default=False, help="plots energy trajectory")
parser.add_option("-b", "--histogram", action="store_true", default=False, help="plots histogram")
parser.add_option("-c", "--native", action="store_true", default=False, help="plots fractional nativeness trajectory")
parser.add_option("-d", "--rmsd", action="store_true", default=False, help="plots rmsd trajectory")
parser.add_option("-e", "--surfE", action="store_true", default=False, help="plots protein-surface trajectory")
parser.add_option("-f", "--nosurfhist", action="store_true", default=False, help="plots protein energy histogram")
parser.add_option("--id", nargs=1, dest="id", type="int", default=0, help="the simlog id number")
parser.add_option("--direc", dest="datafile", default="replicaexchange/simlog0/", help="Qtraj_singleprot.txt file location")
parser.add_option("--show", action="store_false", default=True, help="shows plots")
parser.add_option("--tfile", dest="tfile", default="", help="file of temperatures")
parser.add_option("--subplot", action="store_false", default=True, help="shows plots")


(options,args) = parser.parse_args()
trange = options.temprange #Kelvin
tfile = options.tfile
step = options.step
numreplicas = options.replicas
subplot = options.subplot
energy = options.energy
histogram = options.histogram
native = options.native
rmsd = options.rmsd
surfE = options.surfE
nosurfEhist = options.nosurfhist
id = options.id
direc = options.datafile
show = options.show
T=empty(numreplicas)
alpha=(trange[1]/trange[0])**(1/float(numreplicas-1))
T[0]=trange[0]
for i in range(1,numreplicas):
	T[i]=T[i-1]*alpha
if tfile:
	T = loadtxt(tfile)
print T
colors = cm.gist_rainbow(linspace(0,1,len(T)))[::-1]

if energy:
	files = []
	for i in range(len(T)):
		#files.append('./surface_replica_exchange/replica'+str(i)+'/energy'+str(int(T[i]))+'.txt')
		files.append(direc+'/energy'+str(int(T[i]))+'.npy')
		#files.append('./replicaexchange/replica'+str(i)+'/energy'+str(int(T[i]))+'.txt')
	nc=load(files[0])
	for file in files:
		nctemp=load(file)
		nc=vstack((nc,nctemp))
	nc=nc[1:numreplicas+1,:]
	plt.figure(1)
	for i,e in enumerate(nc):
		plt.plot(range(len(e)), e, label = str('%3.2f') % T[i] + 'K',color=colors[i])
	plt.ylabel('energy (kcal/mol)')
	plt.xlabel('move/%d' %(step))
	plt.title('Go-like model MC simulation')
	plt.legend(loc=2, prop={'size':8})
	plt.savefig(direc+'/energy.png')
	if subplot:
		fig = plt.figure(7)
		for i,e in enumerate(nc):
			plt.subplot(numreplicas/2,2,i+1)
			plt.plot(range(len(e)),e,label=str('3.2f') % T[i] + 'K')
		fig.text(.5, .04, 'move/%i' % step, ha='center', va='center')
		fig.text(.06,.5,'energy (kcal/mol)', ha='center', va='center', rotation='vertical')
		plt.savefig(direc+'/energy_sub.png')
totE=nc.copy()
if histogram:
	if energy:
		plt.figure(2)
		for i,e in enumerate(nc):
			plt.hist(e, 40, label=str('%3.2f') % T[i] + 'K',color=colors[i])
		plt.xlabel('energy (kcal/mol)')
		plt.title('Go-like model MC simulation')
		plt.legend(loc=2, prop={'size':8})
		plt.savefig(direc+'/energyhist.png')
	else:
		print 'still needs to be coded!'
if native:
	files = []
	for i in range(len(T)):
		#files.append('./surface_replica_exchange/replica'+str(i)+'/fractionnative'+str(int(T[i]))+'.txt')
		#files.append('./replicaexchange/replica'+str(i)+'/fractionnative'+str(int(T[i]))+'.txt')
		files.append(direc+'/fractionnative'+str(int(T[i]))+'.npy')
	nc=load(files[0])
	for file in files:
		nctemp=load(file)
		nc=vstack((nc,nctemp))
	nc=nc[1:numreplicas+1,:]
	plt.figure(3)
	for i,e in enumerate(nc):
		plt.plot(range(len(e)), e, label = str('%3.2f') % T[i] + 'K',color=colors[i])
	plt.ylabel('Q, fractional nativeness')
	plt.xlabel('move/%d' %(step))
	plt.title('Go-like model MC simulation')
	plt.legend(loc=2, prop={'size':8})
	plt.savefig(direc+'/fractionnative.png')
	
	if subplot:
		fig = plt.figure(8)
		for i,e in enumerate(nc):
			plt.subplot(numreplicas/2,2,i+1)
			plt.plot(range(len(e)),e,label=str('3.2f') % T[i] + 'K')
		fig.text(.5, .04, 'move/%i' % step, ha='center', va='center')
		fig.text(.06,.5,'energy (kcal/mol)', ha='center', va='center', rotation='vertical')
		plt.savefig(direc+'/fractionnative_sub.png')
if rmsd:
	files = []
	for i in range(len(T)):
		#files.append('./surface_replica_exchange/replica'+str(i)+'/rmsd'+str(int(T[i]))+'.txt')
		files.append(direc+'/rmsd'+str(int(T[i]))+'.npy')
		#files.append('./replicaexchange/replica'+str(i)+'/rmsd'+str(int(T[i]))+'.txt')
	nc=load(files[0])
	for file in files:
		nctemp=load(file)
		nc=vstack((nc,nctemp))
	nc=nc[1:numreplicas+1,:]
	plt.figure(4)
	for i,e in enumerate(nc):
		plt.plot(range(len(e)), e, label = str('%3.2f') % T[i] + 'K')
	plt.ylabel('RMSD')
	plt.xlabel('move/%d' %(step))
	plt.title('Go-like model MC simulation')
	plt.legend(loc=2, prop={'size':8})
	plt.savefig('rmsd.png')

if surfE:
	files = []
	for i in range(len(T)):
		files.append(direc+'/surfenergy'+str(int(T[i]))+'.npy')
	nc=load(files[0])
	for file in files:
		nctemp=load(file)
		nc=vstack((nc,nctemp))
	nc=nc[1:numreplicas+1,:]
	plt.figure(5)
	for i,e in enumerate(nc):
		plt.plot(range(len(e)), e, label = str('%3.2f') % T[i] + 'K')
	plt.ylabel('protein-surface energy (kcal/mol)')
	plt.xlabel('move/%d' %(step))
	plt.title('Go-like model MC simulation')
	plt.legend(loc=2, prop={'size':8})
	plt.savefig('surfE.png')

if nosurfEhist:
	nc = totE-nc
	plt.figure(6)
	for i,e in enumerate(nc):
		plt.hist(e, 40, label=str('%3.2f') % T[i] + 'K')
	plt.xlabel('energy (kcal/mol)')
	plt.title('Go-like model MC simulation')
	plt.legend(loc=2, prop={'size':8})
	plt.savefig('nosurf_energyhist.png')

if show:
	plt.show()
