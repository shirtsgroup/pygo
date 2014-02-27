#!/usr/bin/python2.4

from numpy import *
from checkensemble import *
import matplotlib.pyplot as plt
from optparse import OptionParser
import timeseries

parser=OptionParser()
#parser.add_option("-t", "--temprange", nargs=2, default=[300.0,450.0], type="float", dest="temprange", help="temperature range of replicas")
parser.add_option("--tfile", dest="tfile", default="/home/edz3fz/proteinmontecarlo/T.txt", help="file of temperatures (default: T.txt)")
#parser.add_option("-s", "--stepsize", dest="step", type="int", default='1000', help="number of moves between save operations")
#parser.add_option("-r", "--replicas", default=8, type="int",dest="replicas", help="number of replicas")
parser.add_option("--direc", dest="datafile", default="GO_1PGB.pdb", help="Qtraj_singleprot.txt file")


(options,args) = parser.parse_args()
#trange = options.temprange #Kelvin
#step = options.step
#numreplicas = options.replicas
direc = options.datafile
tfile = options.tfile

T = loadtxt(tfile)

T_k = T[-2::]
files = ['%s/energy%i.npy' % (direc, T[-2]), '%s/energy%i.npy' % (direc, T[-1])]
#file=['/home/edz3fz/checkensemble_high/CE_high.txt','/home/edz3fz/checkensemble_low/CE_low.txt']
#file=[direc+'/energy426.txt',direc+'/energy442.txt']
#file = ['/home/edz3fz/surface_replica_exchange/replica0/energy300.txt', '/home/edz3fz/surface_replica_exchange/replica3/energy356.txt']
down=load(files[0])
up=load(files[1])
length = len(down)
down = down[length/2::]
up = up[length/2::]
#up=up[-50000::]
#down=down[-50000::]
#up=up[::100]
#down=down[::100]

g_up = timeseries.statisticalInefficiency(up)
indices_up = numpy.array(timeseries.subsampleCorrelatedData(up,g=g_up))
print len(indices_up), 'samples'

g_down = timeseries.statisticalInefficiency(down)
indices_down = numpy.array(timeseries.subsampleCorrelatedData(up,g=g_down))
print len(indices_down), 'samples'



type='total'
U_kn=zeros([2,len(up)])
U_kn[0,0:len(indices_down)] = down[indices_down]
U_kn[1,0:len(indices_up)] = up[indices_up]
#T_k=array([300.,336.8472786])
#T_k=array([426.81933819,442.13650313])
#T_k=array([424.67492585,450])
#T_k=array([437.99897735,450])
N_k=[len(indices_up),len(indices_down)]
ProbabilityAnalysis(N_k=N_k, T_k=T_k, U_kn=U_kn, kB=0.0019872041, eunits='kcal/mol', figname='figure', title='name')
#ProbabilityAnalysis(N_k=N_k, T_k=T_k, U_kn=U_kn, eunits='kcal/mol', figname='figure',title='name')
plt.show()
