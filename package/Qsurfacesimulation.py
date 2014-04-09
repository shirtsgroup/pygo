from Qsimulationobject import QSimulation
from surfacesimulation import SurfaceSimulation
from simulationobject import Simulation
import cPickle
import numpy
import pdb
import energyfunc
import moveset
import writetopdb

numpy.random.seed(10)

class QSurfaceSimulation(SurfaceSimulation, QSimulation): # left to right, depth first
    kb = 0.0019872041 #kcal/mol/K

    def __init__(self, name, outputdirectory, coord, temp, surf_coord, Qpin):
        self.Qpin = Qpin
        SurfaceSimulation.__init__(self, name, outputdirectory, coord, temp, surf_coord)
        self.suffix = '%i_%2.2f' % (int(temp), Qpin)

    def setenergy(self):
        SurfaceSimulation.setenergy(self)
        self.Q = energyfunc.nativecontact(self.r2, Simulation.nativeparam, Simulation.nsigma2) / Simulation.totnc
        self.u0 += QSimulation.k_Qpin*(self.Q - self.Qpin)**2

    def update_energy(self, torschange, angchange, dict):
        QSimulation.k_Qpin = dict['k_Qpin']
        Simulation.nsigma2 = dict['nsigma2']
        Simulation.totnc = dict['totnc']
        SurfaceSimulation.update_energy(self, torschange, angchange, dict)
        self.newQ = energyfunc.nativecontact(self.r2new, Simulation.nativeparam, Simulation.nsigma2) / Simulation.totnc
        self.u1 += QSimulation.k_Qpin*(self.newQ - self.Qpin)**2

    def accept_state(self):
        SurfaceSimulation.accept_state(self)
        self.Q = self.newQ

    def output(self):
        QSimulation.output(self)
        if self.trmoves:
            print 'translation:       %d percent acceptance (%i/%i)' %(float(self.acceptedtr)/float(self.trmoves)*100, self.acceptedtr, self.trmoves)
        if self.rmoves:
            print 'rotation:          %d percent acceptance (%i/%i)' %(float(self.acceptedr)/float(self.rmoves)*100, self.acceptedr, self.rmoves)
