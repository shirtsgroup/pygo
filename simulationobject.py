import random
random.seed(10)
import profile
import scipy.misc
import numpy
import writetopdb
import moveset
import energyfunc
try:
	import matplotlib.pyplot as plt
except: pass
import pdb
from sys import stdout


def energy(coord):
	r2, u0 = energyfunc.cgetLJenergy(coord, Simulation.numint, Simulation.numbeads, Simulation.nativeparam_n, Simulation.nonnativeparam, Simulation.nnepsil)
        torsE = energyfunc.ctorsionenergy(coord, numpy.zeros(Simulation.numbeads - 3), Simulation.torsparam, numpy.arange(Simulation.numbeads - 3))
        angE = energyfunc.cangleenergy(coord, numpy.zeros(Simulation.numbeads - 2), Simulation.angleparam, numpy.arange(Simulation.numbeads - 2))
	energy = u0 + sum(torsE) + sum(angE)
	return energy
#def energy(mpos, rsquare, torsE, angE):
#    energy = numpy.sum(angE) + numpy.sum(torsE) + energyfunc.cLJenergy(rsquare, Simulation.nativeparam_n, Simulation.nonnativeparam, Simulation.nnepsil)
#    return energy

def update_energy(self, torschange, angchange):
    self.newtorsE = energyfunc.ctorsionenergy(self.newcoord, self.torsE, Simulation.torsparam, torschange)
    self.newangE = energyfunc.cangleenergy(self.newcoord, self.angE, Simulation.angleparam, angchange)
    self.r2new, self.u1 = energyfunc.cgetLJenergy(self.newcoord, Simulation.numint, Simulation.numbeads, Simulation.nativeparam_n, Simulation.nonnativeparam, Simulation.nnepsil)
    self.u1 += sum(self.newtorsE)+sum(self.newangE)
    return self

def save(self):
    try:
    	self.energyarray[self.move/Simulation.step] = self.u0
    except:
	print self.move
	print Simulation.step
    self.nc[self.move/Simulation.step] = energyfunc.nativecontact(self.r2, Simulation.nativeparam_n, Simulation.nsigma2)
    self.mcoord = writetopdb.moviecoord(self.coord, Simulation.transform)
    self.rmsd_array[self.move/Simulation.step] = energyfunc.rmsd(Simulation.coord_nat, self.mcoord)
    if (Simulation.pdbfile):
        writetopdb.addtopdb(self.coord,Simulation.positiontemplate,self.move/Simulation.step,'%s/trajectory%i.pdb' % (self.out,int(self.T)))
    return self

        
class Simulation:
    kb = 0.0019872041 #kcal/mol/K
    tsteps=360
    tsize=.05
#    percentmove = [0.25, 0.5, 0.75, 1.] # % bend,% axis torsion, % global crankshaft, %ParRot move, %MD

    def __init__(self, name, outputdirectory, coord, temp):
        self.name = name
        self.out = outputdirectory
        self.coord = coord
        self.energyarray = numpy.empty(Simulation.totmoves/Simulation.step + 1)
        self.nc = numpy.empty(Simulation.totmoves/Simulation.step + 1) #native contacts
        self.rmsd_array = numpy.empty(Simulation.totmoves/Simulation.step + 1)
        self.move = 0
        self.T = temp
        self.maxtheta = numpy.array([5., # bend
                        10., # torsion
                        1./self.T * 500, # global crankshaft
                        5.]) # ParRot move
	self.maxtheta = self.maxtheta * numpy.pi / 180 * self.T**1.5 / 5250. * 50 / Simulation.numbeads
	self.setenergy()
        self.energyarray[0] = self.u0
        self.nc[0] = Simulation.totnc
#       self.rmsd_array[0] = 0.0

        # Instantiate constants for move stats
        self.amoves = 0
        self.atmoves = 0
	self.gcmoves = 0
        self.pmoves = 0
        self.mdmoves = 0
	self.acceptedgc =0
        self.acceptedp = 0
        self.acceptedmd = 0
        self.acceptedat = 0
        self.accepteda = 0
	self.pclosure = 0
        self.accepted = 0
        self.rejected = 0
        self.movetype = ''
        #self.whoami = ''
    
    def setenergy(self):
	# sets the u0, r2, torsE, angE from the current configuration
	# called when restarting from a checkpoint
        self.r2, self.u0 = energyfunc.cgetLJenergy(self.coord, Simulation.numint, Simulation.numbeads, Simulation.nativeparam_n, Simulation.nonnativeparam, Simulation.nnepsil)
        self.torsE = energyfunc.ctorsionenergy(self.coord, numpy.zeros(Simulation.numbeads - 3), Simulation.torsparam, numpy.arange(Simulation.numbeads - 3))
        self.angE = energyfunc.cangleenergy(self.coord, numpy.zeros(Simulation.numbeads - 2), Simulation.angleparam, numpy.arange(Simulation.numbeads - 2))
        self.u0 += numpy.sum(self.angE) + numpy.sum(self.torsE)
	

    def output(self, verbose):
	
        print '-------- %s Simulation Results --------' % (self.name)
	print 'total accepted moves: %d' % (self.accepted)
	print 'total rejected moves: %d' % (self.rejected)
	if verbose:
		if self.amoves:
			print 'angle bend:        %d percent acceptance (%i/%i)' %(float(self.accepteda)/float(self.amoves)*100, self.accepteda, self.amoves)
		if self.atmoves:
			print 'torsion:           %d percent acceptance (%i/%i)' %(float(self.acceptedat)/float(self.atmoves)*100, self.acceptedat, self.atmoves)
		if self.gcmoves:
			print 'global crankshaft: %d percent acceptance (%i/%i)' %(float(self.acceptedgc)/float(self.gcmoves)*100, self.acceptedgc, self.gcmoves)
		if self.pmoves:
			print 'ParRot:            %d percent acceptance (%i/%i)' %(float(self.acceptedp)/float(self.pmoves)*100, self.acceptedp, self.pmoves)
			print '                   %i ParRot moves rejected due to chain closure' %(self.pclosure)
		if self.mdmoves:
			print 'MD:                %d percent acceptance (%i/%i)' %(float(self.acceptedmd)/float(self.mdmoves)*100, self.acceptedmd, self.mdmoves)

    def loadstate(self):
	self.coord = numpy.loadtxt('%s/coord%i.cpt' %(self.out, int(self.T)))
	self.setenergy()
	self.energyarray = numpy.loadtxt('%s/energy%i.txt' %(self.out, int(self.T)))
	self.nc = Simulation.totnc*numpy.loadtxt('%s/fractionnative%i.txt' %(self.out, int(self.T)))


    def savecoord(self):
	filename = '%s/coord%i.cpt' % (self.out, int(self.T))
	numpy.savetxt(filename, self.coord)

    def saveenergy(self, plot):
        filename = '%s/energy%i.txt' % (self.out, int(self.T))
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
        filename = '%s/fractionnative%i.txt' % (self.out, int(self.T))
        numpy.savetxt(filename, fraction)
        #print 'wrote every %d fractional nativeness values to %s' %(Simulation.step,fractionfile)
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

def run(self, nummoves, dict):
    Simulation.percentmove = dict['percentmove']   

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
        randmove = random.random()
        randdir = random.random()
        m = random.randint(1, Simulation.numbeads-2) #random bead, not end ones
        uncloseable = False
        
        # angle bend
        if randmove < Simulation.percentmove[0]:
	    theta = self.maxtheta[0] - 2 * self.maxtheta[0] * random.random()
            self.newcoord, jac = moveset.canglebend(self.coord, m, randdir, theta)
            self.amoves += 1
            movetype = 'a'
            torschange = numpy.array([])
            angchange = numpy.array([m-1])

        # axis torsion
        elif randmove < Simulation.percentmove[1]:
            jac = 1
	    theta = self.maxtheta[1] - 2 * self.maxtheta[1] * random.random()
            self.newcoord = moveset.caxistorsion(self.coord, m, randdir, theta)
            movetype = 'at'
            self.atmoves += 1
            angchange = numpy.array([])
            if randdir < .5:
                    torschange = numpy.array([m-2])
                    if m < 2:
                            torschange = numpy.array([])
            elif m == Simulation.numbeads - 2:
                    torschange = numpy.array([])
            else:
                    torschange = numpy.array([m-1])

        # global crankshaft
	elif randmove < Simulation.percentmove[2]:
		jac = 1
		theta = numpy.random.normal(0,self.maxtheta[2],Simulation.numbeads-2)
		self.newcoord = moveset.cglobalcrank(self.coord,theta)
		self.gcmoves += 1
		movetype = 'gc'
		torschange = numpy.arange(Simulation.numbeads-3)
		angchange = numpy.arange(Simulation.numbeads-2)
	
	# parrot move
	elif randmove < Simulation.percentmove[3]:
	    theta = self.maxtheta[3] - 2 * self.maxtheta[3] * random.random()
	    movetype = 'p'
	    self.pmoves += 1
	    angchange = numpy.array([])
	    if m == 1:
		self.newcoord = moveset.caxistorsion(self.coord, m, 1, theta)
		torschange = numpy.array([0])
		jac = 1
	    elif m == Simulation.numbeads - 2:
		self.newcoord = moveset.caxistorsion(self.coord, m, 0, theta)
		torschange = numpy.array([m-2])
		jac = 1
	    else:
		self.newcoord, jac = moveset.parrot(self.coord, m, randdir, theta)
		if numpy.any(numpy.isnan(self.newcoord)):
		    uncloseable = True
		    self.rejected += 1
		    self.pclosure += 1
		    self.move += 1
		    del self.newcoord
		    del jac
	        elif randdir < .5:
		    torschange = numpy.arange(m-2,m+2)
		    torschange = torschange[torschange<(Simulation.numbeads-3)]
		else:
                    torschange = numpy.arange(Simulation.numbeads-m-5,Simulation.numbeads-m-1)
		    torschange = torschange[torschange>-1]
		    torschange = torschange[torschange<(Simulation.numbeads-3)]

        # run molecular dynamics
        else:
            self=moveset.runMD(self, Simulation.tsteps, Simulation.tsize, dict)
            movetype = 'md'
            self.mdmoves += 1
            torschange = numpy.arange(Simulation.numbeads - 3)
            angchange = numpy.arange(Simulation.numbeads - 2)

        # accept or reject
        if not uncloseable and movetype=='md':
            self = update_energy(self, torschange, angchange)
            self.newH=self.u1+.5/4.184*numpy.sum(mass*numpy.sum(self.vel**2,axis=1)) # in kcal/mol
            self.move += 1
            boltz = numpy.exp(-(self.newH-self.oldH)/(Simulation.kb*self.T))
	    alwaysMD = False 
	    if random.random() < boltz or alwaysMD:
                self.accepted += 1
                self.acceptedmd += 1
                self.r2 = self.r2new
                self.coord = self.newcoord
                self.torsE = self.newtorsE
                self.angE = self.newangE
                self.u0 = self.u1
            else:
                self.rejected += 1
        elif not uncloseable:
            self = update_energy(self, torschange, angchange)
            self.move += 1
            boltz = jac*numpy.exp(-(self.u1-self.u0)/(Simulation.kb*self.T))
            if random.random() < boltz:
                self.accepted += 1
                if movetype == 'a':
                    self.accepteda += 1
                elif movetype == 'at':
                    self.acceptedat += 1
                elif movetype == 'gc':
                    self.acceptedgc += 1
		else:
		    self.acceptedp += 1
                self.r2 = self.r2new
                self.coord = self.newcoord
                self.torsE = self.newtorsE
                self.angE = self.newangE
                self.u0 = self.u1
            else:
                self.rejected += 1
        if self.move % Simulation.step == 0:
            self = save(self)
    return self

def run_ff(self, nummoves, dict):
    # no force field
    Simulation.numbeads = dict['numbeads']
    Simulation.step = dict['step']
    Simulation.totmoves = dict['totmoves']

    for i in xrange(nummoves):
        randmove = random.random()
        randdir = random.random()
        m = random.randint(1, Simulation.numbeads-2) #random bead, not end ones
        uncloseable = False
        #bend
        if randmove < self.percentmove[0]:
            theta = self.maxtheta[0] - 2 * self.maxtheta[0] * random.random()
            self.newcoord, jac = moveset.canglebend(self.coord, m, randdir, theta)
            self.amoves += 1
            movetype = 'a'

        #axis torsion
        elif randmove < self.percentmove[1]:
            theta = self.maxtheta[1] - 2 * self.maxtheta[1] * random.random()
            self.newcoord = moveset.caxistorsion(self.coord, m, randdir, theta)
            movetype = 'at'
            self.atmoves += 1
        
        # global crankshaft
        elif randmove < self.percentmove[3]:
            theta = numpy.random.normal(0,self.maxtheta[3],Simulation.numbeads-2)
            self.newcoord = moveset.globalcrank(self.coord, theta)
            self.gcmoves += 1
            movetype = 'gc'
            jac = 1
        
        # parrot
        elif randmove < self.percentmove[5]:
            theta = self.maxtheta[5] - 2 * self.maxtheta[5] * random.random()
            movetype = 'p'
            self.pmoves += 1
            if self.m == 1:
                self.newcoord = moveset.caxistorsion(self.coord, self.m, 1, self.theta)
                jac = 1
            elif self.m == Simulation.numbeads-2:
                self.newcoord = moveset.caxistorsion(self.coord, self.m, 0, self.theta)
                jac = 1
            else:
                self.newcoord, jac = moveset.parrot(self.coord, m, randdir, self.theta)
                if numpy.any(numpy.isnan(self.newcoord)):
                    self.uncloseable = True
                    self.rejected += 1
                    self.closure += 1
                    del self.newcoord
                    del jac
        # run MD
        else:
            self=moveset.runMD(self, 1000, .001, dict)
            self.movetype = 'md'
            self.move += 1
            self.mdmoves += 1
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

























	
