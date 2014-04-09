from simulationobject import Simulation
import scipy.misc
import numpy
import writetopdb
import moveset
import energyfunc
import pdb
from sys import stdout
import cPickle
numpy.random.seed(10)
        
class QSimulation(Simulation):
    kb = 0.0019872041 #kcal/mol/K

    def __init__(self, name, outputdirectory, coord, temp, Qpin):
        self.Qpin = Qpin
        Simulation.__init__(self, name, outputdirectory, coord, temp)
        self.suffix = '%i_%2.2f' % (int(temp), Qpin)

    def setenergy(self):
	# sets the u0, r2, torsE, angE from the current configuration
	# called when restarting from a checkpoint
        Simulation.setenergy(self)
        self.Q = energyfunc.nativecontact(self.r2, Simulation.nativeparam, Simulation.nsigma2) / Simulation.totnc
        self.u0 += QSimulation.k_Qpin*(self.Q - self.Qpin)**2
	
    def update_energy(self, torschange, angchange, dict):
        QSimulation.k_Qpin = dict['k_Qpin']
        Simulation.nsigma2 = dict['nsigma2']
        Simulation.totnc = dict['totnc']
        Simulation.update_energy(self, torschange, angchange, dict)
        self.newQ = energyfunc.nativecontact(self.r2new, Simulation.nativeparam, Simulation.nsigma2) / Simulation.totnc
        self.u1 += QSimulation.k_Qpin*(self.newQ - self.Qpin)**2
    
    def accept_state(self):
        Simulation.accept_state(self)
        self.Q = self.newQ

    def output(self):
        print '-------- %s Simulation Results --------' % (self.name)
        print 'temperature:', self.T
        print 'Qpin:', self.Qpin
        print 'k_Qpin:', QSimulation.k_Qpin
        print 'total accepted moves: %d' % (self.accepted)
        print 'total rejected moves: %d' % (self.rejected)
        if self.amoves:
            print 'angle bend:        %d percent acceptance (%i/%i)' %(float(self.accepteda)/float(self.amoves)*100, self.accepteda, self.amoves)
        if self.atmoves:
            print 'torsion:           %d percent acceptance (%i/%i)' %(float(self.acceptedat)/float(self.atmoves)*100, self.acceptedat, self.atmoves)
        if self.gcmoves:
            print 'crankshaft:        %d percent acceptance (%i/%i)' %(float(self.acceptedgc)/float(self.gcmoves)*100, self.acceptedgc, self.gcmoves)
        if self.pmoves:
            print 'ParRot:            %d percent acceptance (%i/%i)' %(float(self.acceptedp)/float(self.pmoves)*100, self.acceptedp, self.pmoves)
            print '                   %i ParRot moves rejected due to chain closure' %(self.pclosure)
        if self.mdmoves:
            print 'MD:                %d percent acceptance (%i/%i)' %(float(self.acceptedmd)/float(self.mdmoves)*100, self.acceptedmd, self.mdmoves)
