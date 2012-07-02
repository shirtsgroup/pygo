#from numpy import *
#from writetopdb import *
#from moveset import *
#from energyfunc import *
import matplotlib.pyplot as plt
#from sys import stdout
#from random import *
import random
import profile
import scipy.misc
import numpy
import writetopdb
import moveset
import energyfunc
import matplotlib.pyplot as plt
import pdb
from sys import stdout


def energy(mpos, rsquare, torsE, angE):
    energy = sum(angE) + sum(torsE) + energyfunc.cLJenergy(rsquare, Simulation.nativeparam_n, Simulation.nonnativeparam, Simulation.nnepsil)
    return energy


        
        
class Simulation:
    kb = 0.0019872041 #kcal/mol/K
    percentmove = [.3, .6, .4, .4, .9] # % bend,% axis torsion,% crankshaft, %local move, %ParRot move, %MD

    def __init__(self, name, outputdirectory, coord, temp):
        self.name = name
        self.out = outputdirectory
        self.coord = coord
        self.energyarray = numpy.empty(Simulation.totmoves/Simulation.step + 1)
        self.nc = numpy.empty(Simulation.totmoves/Simulation.step + 1) #native contacts
        self.rmsd_array = numpy.empty(Simulation.totmoves/Simulation.step + 1)
        self.move = 0
        self.T = temp
        self.maxtheta = numpy.array([10. * self.T / 300, # bend
                        10. * self.T / 200, # torsion
                        20. * self.T / 250, # crankshaft
                        5. * self.T / 300, # local move
                        5. * self.T/300]) # ParRot move
	self.maxtheta = self.maxtheta * numpy.pi / 180
        self.r2 = energyfunc.cgetLJr2(self.coord, Simulation.numint, Simulation.numbeads)
        self.torsE = energyfunc.ctorsionenergy(self.coord, numpy.zeros(Simulation.numbeads - 3), Simulation.torsparam, numpy.arange(Simulation.numbeads - 3))
        self.angE = energyfunc.angleenergy_n(self.coord, numpy.zeros(Simulation.numbeads - 2), Simulation.angleparam, numpy.arange(Simulation.numbeads - 2))
        self.u0 = energy(self.coord, self.r2, self.torsE, self.angE)
        self.energyarray[0] = self.u0
        self.rmsd_array[0] = 0.0
        self.nc[0] = energyfunc.nativecontact(self.r2, Simulation.nativeparam_n, Simulation.nsigma2)

        # Instantiate constants for move stats
        self.angmoves = 0
        self.crankmoves = 0
        self.atormoves = 0
        self.lmoves = 0
        self.mdmove = 0
        self.pmoves = 0
        self.acceptedp = 0
        self.acceptedmd = 0
        self.acceptedlm = 0
        self.acceptedat = 0
        self.accepteda = 0
        self.acceptedc = 0
        self.accepted = 0
        self.rejected = 0
        self.closure = 0
	self.pclosure = 0
        self.movetype = ''
        #self.whoami = []

    def output(self, verbose):
        write = ['-------- %s Simulation Results --------\r\n' % (self.name),
                'total accepted moves: %d \r\n' % (self.accepted),
                'total rejected moves: %d \r\n' % (self.rejected),
		'ParRot move: %d moves accepted out of %d tries \r\n' % (self.acceptedp, self.pmoves),
		'%d ParRot moves rejected due to chain closure \r\n' % (self.pclosure),
                'local move: %d moves accepted out of %d tries \r\n' % (self.acceptedlm,self.lmoves),
                '%d local moves rejected due to chain closure \r\n' % (self.closure),
                'angle bend: %d moves accepted out of %d tries \r\n' % (self.accepteda,self.angmoves),
                'crankshaft: %d moves accepted out of %d tries \r\n' % (self.acceptedc,self.crankmoves),
                'torsion: %d moves accepted out of %d tries \r\n' % (self.acceptedat,self.atormoves),
                'MD moves: %d moves accepted out of %d tries \r\n' % (self.acceptedmd,self.mdmove)]
        if verbose:
            print "".join(write)
        return write

    def saveenergy(self, plot):
        filename='%s/energy%i.txt' % (self.out, int(self.T))
        numpy.savetxt(filename, self.energyarray)
        #print 'wrote every %d conformation energies to %s' %(Simulation.step,filename)
        if plot:
            plotname = '%s/energy%i.png' % (self.out, int(self.T))
            #print 'generating conformational energy plot...'
            fig = plt.figure(1)
            plt.plot(range(len(self.energyarray)), self.energyarray, label=str('%3.2f') % self.T +' K')
            plt.xlabel('move/%d' % (Simulation.step))
            plt.ylabel('energy (kcal/mol)')
            plt.title('Go-like model monte carlo simulation at %d K' % (self.T))
            plt.grid(True)
            lgd = plt.legend(loc=2, prop={'size':8})
            fig.savefig(plotname)
            #print 'conformational energy plot saved to %s' %(plotname)

    def savermsd(self, plot):
        filename='%s/rmsd%i.txt' % (self.out, int(self.T))
        numpy.savetxt(filename, self.rmsd_array)
        #print 'wrote every %d rmsd values to %s' %(Simulation.step,filename)
        if plot:
            plotname = '%s/rmsd%i.png' % (self.out, int(self.T))
            #print 'generating rmsd plot...'
            fig = plt.figure(2)
            plt.plot(range(len(self.rmsd_array)), self.rmsd_array, label=str('%3.2f') % self.T +' K')
            plt.xlabel('move/%d' % (Simulation.step))
            plt.ylabel('rmsd')
            plt.title('Go-like model monte carlo simulation at %d K' % (self.T))
            plt.grid(True)
            lgd = plt.legend(loc=2, prop={'size':8})
            fig.savefig(plotname)
            #print 'rmsd plot saved to %s' %(plotname)

    def savenc(self, plot):
        fraction = self.nc / Simulation.totnc
        fractionfile = '%s/fractionnative%i.txt' % (self.out, int(self.T))
        numpy.savetxt(fractionfile, fraction)
        #print 'wrote every %d fractional nativeness values to %s' %(Simulation.step,fractionfile)
        print 'excluding first %i out of %i fractional nativeness values from average' % (len(fraction)/5, len(fraction))
        fractioncut = fraction[len(fraction)/5:-1]
        print 'the average fractional nativeness is %f' % (sum(fractioncut) / len(fractioncut))
        if plot:
            plotname = '%s/fractionnative%i.png' % (self.out, int(self.T))
            #print 'generating fractional nativeness plot...'
            fig = plt.figure(3)
            plt.plot(range(len(fraction)), fraction, label=str('%3.2f') % self.T +' K')
            plt.xlabel('move/%d' % (Simulation.step))
            plt.ylabel('Q, fractional nativeness')
            plt.title('Go-like model monte carlo simulation at %d K' % (self.T))
            plt.grid(True)
            lgd = plt.legend(loc=3, prop={'size':8})
            fig.savefig(plotname)
            #print 'fractional nativeness plot saved to %s' %(plotname)

    def savehist(self, plot):
        if plot:
            plotname = '%s/energyhist%i.png' % (self.out, int(self.T))
            #print 'generating energy histogram...'
            fig = plt.figure(4)
            plt.hist(self.energyarray, 40, label=str('%3.2f') % self.T +' K')
            plt.xlabel('energy (kcal/mol)')
            plt.title('Go-like model MC simulation at %d K' % (self.T))
            plt.legend(loc=2, prop={'size':8})
            fig.savefig(plotname)
            #print 'energy histogram saved to %s' %(plotname)
	
class SurfaceSimulation(Simulation):
    def addsurface(self, surf_coord, dist_to_surf):
        self.surface = surf_coord
        ## randomly orient protein
        rotx = numpy.arccos(2*random.random()-1)
        rotz = 2*numpy.pi*random.random()
        com = sum(self.coord, axis=0)/float(numbeads)
        self.coord -= com
        self.coord = self.coord.transpose()
        self.coord = numpy.array([self.coord[0,:]*numpy.cos(rotz)+self.coord[1,:]*numpy.sin(rotz), -self.coord[0,:]*numpy.sin(rotz)+self.coord[1,:]*numpy.cos(rotz), self.coord[2,:]])
        self.coord = numpy.array([self.coord[0,:], numpy.cos(rotx)*self.coord[1,:]+numpy.sin(rotx)*self.coord[2,:], -numpy.sin(rotx)*self.coord[1,:] + numpy.cos(rotx)*self.coord[2,:]])
        #self.coord = self.coord.transpose() + com
        self.coord[:,2] += (dist_to_surf - min(self.coord[:,2]))
        
def run(self, nummoves, dict):
    Simulation.numbeads = dict['numbeads']
    Simulation.step = dict['step']
    Simulation.totmoves = dict['totmoves']
    Simulation.numint = dict['numint']
    Simulation.angleparam = dict['angleparam']
    Simulation.torsparam = dict['torsparam']
    Simulation.nativeparam_n = dict['nativeparam_n']
    Simulation.nonnativeparam = dict['nonnativeparam']
    Simulation.nnepsil = dict['nnepsil']
    Simulation.nsigma2 = dict['nsigma2']
    Simulation.transform = dict['transform']
    Simulation.coord_nat = dict['coord_nat']
    Simulation.positiontemplate = dict['positiontemplate']
    Simulation.pdbfile = dict['pdbfile']
    mass = dict['mass']

    for i in xrange(nummoves):
        self.randmove = random.random()
        self.randdir = random.random()
        self.m = random.randint(1, Simulation.numbeads-2) #random bead, not end ones
        self.uncloseable = False
        
        # angle bend
        if self.randmove < Simulation.percentmove[0]:
	    self.theta = self.maxtheta[0] - 2 * self.maxtheta[0] * random.random()
            self.newcoord, self.jac = moveset.canglebend(self.coord, self.m, self.randdir, self.theta)
            self.angmoves += 1
            self.movetype = 'a'
            self.change = []
            self.angchange = [self.m-1]

        # axis torsion
        elif self.randmove < Simulation.percentmove[1]:
            self.jac = 1
	    self.theta = self.maxtheta[1] - 2 * self.maxtheta[1] * random.random()
            self.newcoord = moveset.caxistorsion(self.coord, self.m, self.randdir, self.theta)
            self.movetype = 'at'
            self.atormoves += 1
            self.angchange = []
            if self.randdir < .5:
                    self.change = [self.m-2]
                    if self.m < 2:
                            self.change = []
            elif self.m == Simulation.numbeads - 2:
                    self.change = []
            else:
                    self.change = [self.m-1]

        # crankshaft
        elif self.randmove < Simulation.percentmove[2]:
            self.jac = 1
	    self.theta = self.maxtheta[2] - 2 * self.maxtheta[2] * random.random()
            self.newcoord = moveset.crankshaft(self.coord, self.m, self.theta)
            self.crankmoves += 1
            self.movetype = 'c'
            self.change = numpy.arange(self.m-3, self.m+1)
            self.angchange = [self.m-2, self.m]
            if self.m == 2:
                    self.change = numpy.arange(self.m-2, self.m+1)
            elif self.m < 2:
                    self.change = self.change[self.change > -1]
                    self.angchange = [self.m]
            elif self.m == Simulation.numbeads-4 or self.m == Simulation.numbeads-3:
                    self.change = self.change[self.change < (Simulation.numbeads-3)]
            elif self.m == Simulation.numbeads-2:
                    self.change = numpy.arange(self.m-3, self.m-1)
                    self.angchange = [self.m-2]

        # local move
        elif self.randmove < Simulation.percentmove[3]:
	    self.jac = 1
            self.theta = self.maxtheta[3] - 2 * self.maxtheta[3] * random.random()
            self.movetype = 'lm'
            self.lmoves += 1
            if self.m < 5 and self.randdir > .5:
                self.newcoord = moveset.caxistorsion(self.coord, self.m, self.randdir, self.theta)
                self.angchange = []
                self.change = [self.m-1]
            elif self.m > Simulation.numbeads-6 and self.randdir < .5:
                self.newcoord = moveset.caxistorsion(self.coord, self.m, self.randdir, self.theta)
                self.angchange = []
                self.change = [self.m-2]
            else:
                self.newcoord = moveset.localmove(self.coord, self.m, self.randdir, self.theta)
                if numpy.any(numpy.isnan(self.newcoord)):
                    self.uncloseable = True
                    self.rejected += 1
                    self.closure += 1
                    self.move += 1
                    del self.newcoord
                elif self.randdir < .5: #toward end of chain
                    self.angchange = numpy.arange(self.m+2, self.m+5)
                    self.change = numpy.array([self.m-2, self.m+1, self.m+2, self.m+3, self.m+4])
                    self.angchange = self.angchange[self.angchange < (Simulation.numbeads-2)]
                    self.change = self.change[self.change < (Simulation.numbeads-3)]
                    self.angchange = self.angchange[self.angchange > -1]
                    self.change = self.change[self.change > -1]
                else: # toward beginning of chain
                    self.angchange = numpy.arange(self.m-6, self.m-3)
                    self.change = numpy.array([self.m-7, self.m-6, self.m-5, self.m-4, self.m-1])
                    self.angchange = self.angchange[self.angchange > -1]
                    self.change = self.change[self.change > -1]
                    self.angchange = self.angchange[self.angchange < (Simulation.numbeads-2)]
                    self.change = self.change[self.change < (Simulation.numbeads-3)]
	
	# parrot move
	elif self.randmove < Simulation.percentmove[4]:
	    self.theta = self.maxtheta[4] - 2 * self.maxtheta[4] * random.random()
	    self.movetype = 'p'
	    self.pmoves += 1
	    self.angchange = []
	    if self.m == 1:
		self.newcoord = moveset.caxistorsion(self.coord, self.m, 1, self.theta)
		self.change = [0]
		self.jac = 1
	    elif self.m == Simulation.numbeads - 2:
		self.newcoord = moveset.caxistorsion(self.coord, self.m, 0, self.theta)
		self.change = [self.m-2]
		self.jac = 1
	    else:
		self.newcoord, self.jac = moveset.parrot(self.coord, self.m, self.randdir, self.theta)
		if numpy.any(numpy.isnan(self.newcoord)):
		    self.uncloseable = True
		    self.rejected += 1
		    self.pclosure += 1
		    self.move += 1
		    del self.newcoord
		    del self.jac
	        elif self.randdir < .5:
		    self.change = numpy.arange(self.m-2,self.m+2)
		    self.change = self.change[self.change<(Simulation.numbeads-3)]
		else:
		    #self.change = numpy.arange(self.m-4,self.m)
                    self.change = numpy.arange(Simulation.numbeads-1-self.m-4,Simulation.numbeads-1-self.m)
		    self.change = self.change[self.change>-1]
		    self.change = self.change[self.change<(Simulation.numbeads-3)]

        # run molecular dynamics
        else:
            self=moveset.runMD(self, 5, .01, dict)
            self.movetype = 'md'
            self.move += 1
            self.mdmove += 1

        # accept or reject
        if self.uncloseable == False and self.movetype == 'md':
            self.r2new = energyfunc.cgetLJr2(self.newcoord, Simulation.numint, Simulation.numbeads)
            self.newtorsE = energyfunc.ctorsionenergy(self.newcoord, numpy.zeros(Simulation.numbeads-3),Simulation.torsparam,numpy.arange(Simulation.numbeads-3))
            self.newangE = energyfunc.angleenergy_n(self.newcoord, numpy.zeros(Simulation.numbeads-2),Simulation.angleparam,numpy.arange(Simulation.numbeads-2))
            self.u1 = energy(self.newcoord, self.r2new, self.newtorsE, self.newangE)
            #self.newH  =self.u1 + .5 / 4.184 * numpy.sum(mass * numpy.sum(self.vel**2, axis=1))
            self.boltz = numpy.exp(-(self.u1-self.u0)/(Simulation.kb*self.T))
            if random.random() < self.boltz:
                self.accepted += 1
                self.acceptedmd += 1
                self.r2 = self.r2new
                self.coord = self.newcoord
                self.torsE = self.newtorsE
                self.angE = self.newangE
                self.u0 = self.u1
            else:
                self.rejected += 1
        elif(self.uncloseable == False and self.movetype != 'md'):
            self.r2new = energyfunc.cgetLJr2(self.newcoord, Simulation.numint, Simulation.numbeads)
            self.newangE = energyfunc.angleenergy_n(self.newcoord, self.angE, Simulation.angleparam,self.angchange)
            self.newtorsE = energyfunc.ctorsionenergy(self.newcoord, self.torsE, Simulation.torsparam,self.change)
            self.u1 = energy(self.newcoord, self.r2new, self.newtorsE, self.newangE)
            self.move += 1
            self.boltz = self.jac*numpy.exp(-(self.u1-self.u0)/(Simulation.kb*self.T))
            if random.random() < self.boltz:
                self.accepted += 1
                if self.movetype == 'a':
                    self.accepteda += 1
                elif self.movetype == 'r':
                    self.acceptedr += 1
                elif self.movetype == 'at':
                    self.acceptedat += 1
                    #accepted_angle.append(theta)
                elif self.movetype == 'c': 
                    self.acceptedc += 1
                elif self.movetype == 'lm':
                    self.acceptedlm += 1
		else:
		    self.acceptedp += 1
                #if (pdbfile != ''):
                    #mcoord=moviecoord(newcoord,transform)
                    ##writeseqpdb(mcoord,wordtemplate,ATOMlinenum,accepted)
                    #addtopdb(mcoord,positiontemplate,move,pdbfile)
                self.r2 = self.r2new
                self.coord = self.newcoord
                self.torsE = self.newtorsE
                self.angE = self.newangE
                self.u0 = self.u1
            else:
                self.rejected += 1
        if self.move%Simulation.step == 0:
            # energy array
            self.energyarray[self.move/Simulation.step] = self.u0
            # native contact array
            self.nc[self.move/Simulation.step] = energyfunc.nativecontact(self.r2, Simulation.nativeparam_n, Simulation.nsigma2)
            # rmsd array
            self.mcoord = writetopdb.moviecoord(self.coord, Simulation.transform)
            self.rmsd_array[self.move/Simulation.step] = energyfunc.rmsd(Simulation.coord_nat, self.mcoord)
            if (Simulation.pdbfile):
                pass
                #writetopdb.addtopdb(self.coord,Simulation.positiontemplate,self.move/Simulation.step,'%s/trajectory%i.pdb' % (self.out,int(self.T)))
    return self

def run_ff(self, nummoves, dict):
    # no force field
    Simulation.numbeads = dict['numbeads']
    Simulation.step = dict['step']
    Simulation.totmoves = dict['totmoves']

    for i in xrange(nummoves):
        self.randmove = random.random()
        self.randdir = random.random()
        self.m = random.randint(1, Simulation.numbeads-2) #random bead, not end ones
        self.uncloseable = False
        #bend
        if self.randmove < self.percentmove[0]:
            self.theta = self.maxtheta[0] - 2 * self.maxtheta[0] * random.random()
            self.newcoord, jac = moveset.canglebend(self.coord, self.m, self.randdir, self.theta)
            self.angmoves += 1
            self.movetype = 'a'

        #axis torsion
        elif self.randmove < self.percentmove[1]:
            self.theta = self.maxtheta[1] - 2 * self.maxtheta[1] * random.random()
            self.newcoord = moveset.caxistorsion(self.coord, self.m, self.randdir, self.theta)
            self.movetype = 'at'
            self.atormoves += 1

        #crankshaft
        elif self.randmove < self.percentmove[2]:
            self.theta = self.maxtheta[2] - 2 * self.maxtheta[2] * random.random()
            self.newcoord = moveset.crankshaft(self.coord, self.m, self.theta)
            self.crankmoves += 1
            self.movetype = 'c'

        #local move
        elif self.randmove < self.percentmove[3]:
            self.theta = self.maxtheta[3] - 2 * self.maxtheta[3] * random.random()
            self.movetype = 'lm'
            self.lmoves += 1
            if self.m < 5 and self.randdir > .5:
                self.newcoord = moveset.caxistorsion(self.coord, self.m, self.randdir, self.theta)
            elif self.m > Simulation.numbeads-6 and self.randdir < .5:
                self.newcoord = moveset.caxistorsion(self.coord, self.m, self.randdir, self.theta)
            else:
                self.newcoord = moveset.localmove(self.coord, self.m, self.randdir, self.theta)
                if numpy.any(numpy.isnan(self.newcoord)):
                    self.uncloseable = True
                    self.rejected += 1
                    self.closure += 1
                    del self.newcoord

        # parrot
        elif self.randmove < self.percentmove[4]:
            self.theta = self.maxtheta[3] - 2 * self.maxtheta[3] * random.random()
            self.movetype = 'p'
            self.pmoves += 1
            if self.m == 1:
                self.newcoord = moveset.caxistorsion(self.coord, self.m, 1, self.theta)
                jac = 1
            elif self.m == Simulation.numbeads-2:
                self.newcoord = moveset.caxistorsion(self.coord, self.m, 0, self.theta)
                jac = 1
            else:
                self.newcoord, jac = moveset.parrot(self.coord, self.m, self.randdir, self.theta)
                if numpy.any(numpy.isnan(self.newcoord)):
                    self.uncloseable = True
                    self.rejected += 1
                    self.closure += 1
                    del self.newcoord
                    del jac
        # run MD
        else:
            self=moveset.runMD(self, 20, .1, dict)
            self.movetype = 'md'
            self.move += 1
            self.mdmove += 1
            jac = 1
        # accept or reject
        if self.uncloseable == False:
            self.move += 1
            stdout.write(str(self.move)+'\r')
            stdout.flush()
            if random.random() < jac:
                self.accepted += 1
                if self.movetype == 'a':
                    self.accepteda += 1
                elif self.movetype == 'r':
                    self.acceptedr += 1
                elif self.movetype == 'at':
                    self.acceptedat += 1
                    #accepted_angle.append(theta)
                elif self.movetype == 'c': 
                    self.acceptedc += 1
                else:
                    self.acceptedlm += 1
                #if (pdbfile != ''):
                    #mcoord=moviecoord(newcoord,transform)
                    ##writeseqpdb(mcoord,wordtemplate,ATOMlinenum,accepted)
                    #addtopdb(mcoord,positiontemplate,move,pdbfile)
                self.coord = self.newcoord
        if self.move%Simulation.step == 0:
            dihed=numpy.array(energyfunc.cdihedral(self.coord))
            self.dihedarr=numpy.vstack((self.dihedarr,dihed))
	    ang=numpy.array(energyfunc.angle(self.coord))
	    self.angarr=numpy.vstack((self.angarr, ang))
    return self




























	
