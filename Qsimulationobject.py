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
   
    def setenergy(self):
	# sets the u0, r2, torsE, angE from the current configuration
	# called when restarting from a checkpoint
        Simulation.setenergy(self)
        self.Q = energyfunc.nativecontact(self.r2, Simulation.nativeparam, Simulation.nsigma2) / Simulation.totnc
        self.u0 += QSimulation.k_Qpin*(self.Q - self.Qpin)**2
	
    def loadstate(self):
        self.coord = numpy.load('%s/coord%i.npy' %(self.out, int(self.T)))
        self.setenergy()
        self.energyarray = numpy.load('%s/energy%i.npy' %(self.out, int(self.T)))
        self.nc = numpy.load('%s/fractionnative%i.npy' %(self.out, int(self.T)))

    def loadextend(self,extenddirec):
        self.coord = numpy.load('%s/coord%i.npy' %(extenddirec, int(self.T)))
        self.setenergy()
        self.energyarray[0]=self.u0
        self.nc[0] = self.Q

    def savenc(self):
        filename = '%s/fractionnative%i' % (self.out, int(self.T))
        numpy.save(filename, self.nc)

    def update_energy(self, torschange, angchange, dict):
        QSimulation.k_Qpin = dict['k_Qpin']
        Simulation.nsigma2 = dict['nsigma2']
        Simulation.totnc = dict['totnc']
        Simulation.update_energy(self, torschange, angchange, dict)
        self.newQ = energyfunc.nativecontact(self.r2new, Simulation.nativeparam, Simulation.nsigma2) / Simulation.totnc
        self.u1 += QSimulation.k_Qpin*(self.newQ - self.Qpin)**2
    
    def save_state(self, dict):
        Simulation.save = dict['save']
        Simulation.nativeparam = dict['nativeparam']
        Simulation.nsigma2 = dict['nsigma2']
        Simulation.writetraj = dict['writetraj']
        self.energyarray[self.move/Simulation.save] = self.u0
        self.nc[self.move/Simulation.save] = self.Q 
        if (Simulation.writetraj):
            f = open('%s/trajectory%i' %(self.out, int(self.T)), 'ab')
            numpy.save(f,self.coord)
            f.close()

    def accept_state(self):
        Simulation.accept_state(self)
        self.Q = self.newQ

#    def run(self, nummoves, dict):
#        #differs from Simulation.run only by self.Q updating
#        Simulation.percentmove = dict['percentmove']   
#        Simulation.numbeads = dict['numbeads']
#        Simulation.save = dict['save']
#        Simulation.totmoves = dict['totmoves']
#        Simulation.numint = dict['numint']
#        mass = dict['mass']
#        Simulation.tsteps = dict['tsteps']
#        Simulation.tsize = dict['tsize']
#    
#        for i in xrange(nummoves):
#            randmove = numpy.random.random()
#            randdir = numpy.random.random()
#            m = numpy.random.randint(1, Simulation.numbeads-1) #random bead, not end ones
#            uncloseable = False
#            
#            # angle bend
#            if randmove < Simulation.percentmove[0]:    
#                theta = self.maxtheta[0] - 2 * self.maxtheta[0] * numpy.random.random()
#                self.newcoord, jac = moveset.canglebend(self.coord, m, randdir, theta)
#                self.amoves += 1
#                movetype = 'a'
#                torschange = numpy.array([])
#                angchange = numpy.array([m-1])
#    
#            # axis torsion
#            elif randmove < Simulation.percentmove[1]:
#                jac = 1
#                theta = self.maxtheta[1] - 2 * self.maxtheta[1] * numpy.random.random()
#                self.newcoord = moveset.caxistorsion(self.coord, m, randdir, theta)
#                movetype = 'at'
#                self.atmoves += 1
#                angchange = numpy.array([])
#                if randdir < .5:
#                        torschange = numpy.array([m-2])
#                        if m < 2:
#                                torschange = numpy.array([])
#                elif m == Simulation.numbeads - 2:
#                        torschange = numpy.array([])
#                else:
#                        torschange = numpy.array([m-1])
#    
#            # global crankshaft
#            elif randmove < Simulation.percentmove[2]:
#                jac = 1
#                theta = numpy.random.normal(0,self.maxtheta[2],Simulation.numbeads-2)
#                self.newcoord = moveset.cglobalcrank(self.coord,theta)
#                self.gcmoves += 1
#                movetype = 'gc'
#                torschange = numpy.arange(Simulation.numbeads-3)
#                angchange = numpy.arange(Simulation.numbeads-2)
#    	
#    	    # parrot move
#            elif randmove < Simulation.percentmove[3]:
#                theta = self.maxtheta[3] - 2 * self.maxtheta[3] * numpy.random.random()
#                movetype = 'p'
#                self.pmoves += 1
#                angchange = numpy.array([])
#                if m == 1:
#                    self.newcoord = moveset.caxistorsion(self.coord, m, 1, theta)
#                    torschange = numpy.array([0])
#                    jac = 1
#                elif m == Simulation.numbeads - 2:
#                    self.newcoord = moveset.caxistorsion(self.coord, m, 0, theta)
#                    torschange = numpy.array([m-2])
#                    jac = 1
#                else:
#                    self.newcoord, jac = moveset.parrot(self.coord, m, randdir, theta)
#                    if numpy.any(numpy.isnan(self.newcoord)):
#                        uncloseable = True
#                        self.rejected += 1
#                        self.pclosure += 1
#                        self.move += 1
#                        del self.newcoord
#                        del jac
#                    elif randdir < .5:
#                        torschange = numpy.arange(m-2,m+2)
#                        torschange = torschange[torschange<(Simulation.numbeads-3)]
#                    else:
#                        torschange = numpy.arange(Simulation.numbeads-m-5,Simulation.numbeads-m-1)
#                        torschange = torschange[torschange>-1]
#                        torschange = torschange[torschange<(Simulation.numbeads-3)]
#    
#            # run molecular dynamics
#            else:
#                self=moveset.runMD(self, Simulation.tsteps, Simulation.tsize, dict)
#                movetype = 'md'
#                self.mdmoves += 1
#                torschange = numpy.arange(Simulation.numbeads - 3)
#                angchange = numpy.arange(Simulation.numbeads - 2)
#    
#            # accept or reject
#            if not uncloseable and movetype=='md':
#                self.update_energy(torschange, angchange, dict)
#                self.newH=self.u1+.5/4.184*numpy.sum(mass*numpy.sum(self.vel**2,axis=1)) # in kcal/mol
#                self.move += 1
#                boltz = numpy.exp(-(self.newH-self.oldH)/(Simulation.kb*self.T))
#                if numpy.random.random() < boltz:
#                    self.accepted += 1
#                    self.acceptedmd += 1
#                    self.r2 = self.r2new
#                    self.coord = self.newcoord
#                    self.torsE = self.newtorsE
#                    self.angE = self.newangE
#                    self.u0 = self.u1
#                    self.Q = self.newQ
#                else:
#                    self.rejected += 1
#            elif not uncloseable:
#                self.update_energy(torschange, angchange, dict)
#                self.move += 1
#                boltz = jac*numpy.exp(-(self.u1-self.u0)/(Simulation.kb*self.T))
#                if numpy.random.random() < boltz:
#                    self.accepted += 1
#                    if movetype == 'a':
#                        self.accepteda += 1
#                    elif movetype == 'at':
#                        self.acceptedat += 1
#                    elif movetype == 'gc':
#                        self.acceptedgc += 1
#                    else:
#                        self.acceptedp += 1
#                    self.r2 = self.r2new
#                    self.coord = self.newcoord
#                    self.torsE = self.newtorsE
#                    self.angE = self.newangE
#                    self.u0 = self.u1
#                    self.Q = self.newQ
#                else:
#                    self.rejected += 1
#            if self.move % Simulation.save == 0:
#                self.save_state(dict)
#        return self
