from simulationobject import Simulation
import numpy
import pdb
import random
import energyfunc
import moveset
import writetopdb
try:
	import matplotlib.pyplot as plt
except: pass

def getsurf(xsize, ysize, spacing):
    """Generates coordinates for a surface given surface parameters; uses hexangonal packing"""
    xlen = xsize / spacing
    ylen = ysize / (spacing * 3.**.5 / 2.) # more dense from hexagonal packing
    surfcoord = numpy.zeros((int(xlen)*int(ylen),3))
    for j in range(1,xlen):
        surfcoord[j,:] = surfcoord[j-1,:] + numpy.array([spacing,0,0])
    k = 0
    for i in range(xlen,len(surfcoord),xlen):
        surfcoord[i,:] = surfcoord[i-xlen,:] + numpy.array([(-1)**k*spacing/2., spacing*3.**.5/2., 0]) 
        k += 1
        for j in range(1,xlen):
            surfcoord[i+j,:] = surfcoord[i+j-1,:] + numpy.array([spacing,0,0])
    com = numpy.sum(surfcoord, axis=0)/len(surfcoord)
    surfcoord -= com # centers surface at (0,0,0)
    return surfcoord

def writesurf(filename, surf_coord):
    fout = open(filename, 'w')
    fout.write('HEADER   This pdb file contains the coordinates for a Go model of surface\r\n')
    whitespace = "".join([" " for i in range(25)])
    for row in surf_coord:
        line = ['ATOM', whitespace, '%8.3f' % row[0], '%8.3f' % row[1], '%8.3f' % row[2], '\r\n']
        fout.write("".join(line))
    fout.close()

def update_energy(self, torschange, angchange):
    self.newtorsE = energyfunc.ctorsionenergy(self.newcoord, self.torsE, Simulation.torsparam, torschange)
    self.newangE = energyfunc.cangleenergy(self.newcoord, self.angE, Simulation.angleparam, angchange)
    self.newsurfE = energyfunc.csurfenergy(self.newcoord, SurfaceSimulation.surface, Simulation.numbeads, SurfaceSimulation.nspint, SurfaceSimulation.surfparam)
    self.r2new, self.u1 = energyfunc.cgetLJenergy(self.newcoord, Simulation.numint, Simulation.numbeads, Simulation.nativeparam_n, Simulation.nonnativeparam, Simulation.nnepsil)
    self.u1 += numpy.sum(self.newtorsE)+numpy.sum(self.newangE)+self.newsurfE
    return self

def save(self):
    index = self.move/Simulation.step
    self.energyarray[index] = self.u0
    self.surfE_array[index] = self.surfE
    self.nc[index] = energyfunc.nativecontact(self.r2, Simulation.nativeparam_n, Simulation.nsigma2)
    #self.radgyr[self.move/Simulation.step] = energyfunc.radgyr(self.coord)
    self.mcoord = writetopdb.moviecoord(self.coord, Simulation.transform)
    self.rmsd_array[index] = energyfunc.rmsd(self.coord_nat, self.mcoord)
    if (Simulation.pdbfile):
        writetopdb.addtopdb(self.coord,Simulation.positiontemplate,index,'%s/trajectory%i.pdb' % (self.out,int(self.T)))
    return self
                

class SurfaceSimulation(Simulation):
    kb = 0.0019872041 #kcal/mol/K
    percentmove = [.1, .2, .3, .4, .75, 1] # % translation, % rotation % bend,% axis torsion,% global crankshaft, %ParRot move, %MD

    def __init__(self, name, outputdirectory, coord, temp, surf_coord):
        Simulation.__init__(self, name, outputdirectory, coord, temp)
        self.addsurface(surf_coord)
        self.surfE = energyfunc.csurfenergy(self.coord, surf_coord, Simulation.numbeads, SurfaceSimulation.nspint, SurfaceSimulation.surfparam)
        self.surfE_array = numpy.empty(Simulation.totmoves/Simulation.step + 1)
        self.surfE_array[0] = self.surfE
        #self.radgyr = numpy.empty(Simulation.totmoves/Simulation.step + 1)
        #self.radgyr[0] = energyfunc.radgyr(self.coord)
        self.moveparam = numpy.array([2., # translation
                        2., # rotation
                        10., # bend
                        10., # torsion
                        1., # global crankshaft
                        5.]) # ParRot move
        self.moveparam = self.moveparam * numpy.pi / 180 * self.T / 300 * 50 / Simulation.numbeads
	self.u0 += self.surfE
        self.energyarray[0] = self.u0
        self.rmsd_array[0] = 0.0
        self.coord_nat = writetopdb.moviecoord(self.coord, Simulation.transform)
        self.nc[0] = Simulation.totnc
        self.trmoves = 0
        self.acceptedtr = 0
        self.rotmoves = 0
        self.acceptedrot = 0
    
    def addsurface(self, surf_coord):
        self.surface = surf_coord
        ## randomly orient protein
        rotx = numpy.arccos(2*random.random()-1)
        rotz = 2*numpy.pi*random.random()
        com = numpy.sum(self.coord, axis=0)/float(Simulation.numbeads)
        self.coord -= com
        self.coord = self.coord.transpose()
        self.coord = numpy.array([self.coord[0,:]*numpy.cos(rotz)+self.coord[1,:]*numpy.sin(rotz), -self.coord[0,:]*numpy.sin(rotz)+self.coord[1,:]*numpy.cos(rotz), self.coord[2,:]])
        self.coord = numpy.array([self.coord[0,:], numpy.cos(rotx)*self.coord[1,:]+numpy.sin(rotx)*self.coord[2,:], -numpy.sin(rotx)*self.coord[1,:] + numpy.cos(rotx)*self.coord[2,:]])
        self.coord = self.coord.transpose()
        self.coord[:,2] += 10 - numpy.min(self.coord[:,2]) # protein is 1 nm from surface
        
    def output(self, verbose):
        try: percenta = float(self.accepteda)/float(self.angmoves)*100
	except:	percenta = 0
        try: percentat = float(self.acceptedat)/float(self.atormoves)*100
	except:	percentat = 0
        try: percentgc = float(self.acceptedgc)/float(self.gcmoves)*100
	except:	percentgc = 0
        try: percentp = float(self.acceptedp)/float(self.pmoves)*100
	except:	percentp = 0
        try: percentmd = float(self.acceptedmd)/float(self.mdmove)*100
	except:	percentmd = 0
        try: percenttr = float(self.acceptedtr)/float(self.trmoves)*100
	except:	percenttr = 0
        try: percentrot = float(self.acceptedrot)/float(self.rotmoves)*100
	except:	percentrot = 0
        write = ['-------- %s Simulation Results --------\r\n' % (self.name),
                'total accepted moves: %d \r\n' % (self.accepted),
                'total rejected moves: %d \r\n' % (self.rejected),
		'Translation: %d moves accepted out of %d tries... that is %d percent \r\n' % (self.acceptedtr, self.trmoves, percenttr),
		'Rotation: %d moves accepted out of %d tries... that is %d percent \r\n' % (self.acceptedrot, self.rotmoves, percentrot),
                'angle bend: %d moves accepted out of %d tries... that is %d percent \r\n' % (self.accepteda, self.angmoves, percenta),
                'torsion: %d moves accepted out of %d tries... that is %d percent \r\n' % (self.acceptedat, self.atormoves, percentat),
		'global crankshaft: %d moves accepted out of %d tries... that is %d percent \r\n' % (self.acceptedgc, self.gcmoves, percentgc),
		'ParRot move: %d moves accepted out of %d tries... that is %d percent \r\n' % (self.acceptedp, self.pmoves, percentp),
		'%d ParRot moves rejected due to chain closure \r\n' % (self.pclosure),
                'MD moves: %d moves accepted out of %d tries... that is %d percent \r\n' % (self.acceptedmd, self.mdmove, percentmd)]
        if verbose:
            print "".join(write)
        return write
        
    def savesurfenergy(self, plot):
        filename='%s/surfenergy%i.txt' % (self.out, int(self.T))
        numpy.savetxt(filename, self.surfE_array)
        #print 'wrote every %d conformation energies to %s' %(Simulation.step,filename)
        if plot:
            plotname = '%s/surfenergy%i.png' % (self.out, int(self.T))
            #print 'generating conformational energy plot...'
            fig = plt.figure(5)
            plt.plot(range(len(self.surfE_array)), self.surfE_array, label=str('%3.2f') % self.T +' K')
            plt.xlabel('move/%d' % (Simulation.step))
            plt.ylabel('protein-surface energy (kcal/mol)')
            plt.title('Go-like model monte carlo simulation at %d K' % (self.T))
            plt.grid(True)
            lgd = plt.legend(loc=2, prop={'size':8})
            fig.savefig(plotname)
            #print 'conformational energy plot saved to %s' %(plotname)

    def saveradgyr(self, plot):
        filename='%s/radgyr%i.txt' % (self.out, int(self.T))
        numpy.savetxt(filename, self.radgyr)
        #print 'wrote every %d rmsd values to %s' %(Simulation.step,filename)
        if plot:
            plotname = '%s/radgyr%i.png' % (self.out, int(self.T))
            #print 'generating rmsd plot...'
            fig = plt.figure(6)
            plt.plot(range(len(self.radgyr)), self.radgyr, label=str('%3.2f') % self.T +' K')
            plt.xlabel('move/%d' % (Simulation.step))
            plt.ylabel('radius of gyration (Angstroms)')
            plt.title('Go-like model monte carlo simulation at %d K' % (self.T))
            plt.grid(True)
            lgd = plt.legend(loc=2, prop={'size':8})
            fig.savefig(plotname)
            #print 'rmsd plot saved to %s' %(plotname)

def run_surf(self, nummoves, dict):
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
    Simulation.positiontemplate = dict['positiontemplate']
    Simulation.pdbfile = dict['pdbfile']
    SurfaceSimulation.surface = dict['surface']
    SurfaceSimulation.nspint = dict['nspint']
    SurfaceSimulation.nsurf = dict['nsurf']
    SurfaceSimulation.surfparam = dict['surfparam']
    xlength = dict['xlength']
    ylength = dict['ylength']
    spacing = dict['spacing']
    yspacing = dict['yspacing']
    mass = dict['mass']
    
    for i in xrange(nummoves):
        randmove = random.random()
        randdir = random.random()
        m = random.randint(1, Simulation.numbeads-2) #random bead, not end ones
        uncloseable = False

        # translation
        if randmove < SurfaceSimulation.percentmove[0]:
            jac = 1
            [x,y,z] = numpy.random.normal(0, self.moveparam[0], 3)
            self.newcoord = moveset.translation(self.coord, x, y, z)
            self.trmoves += 1
            movetype = 'tr'
            torschange = numpy.array([])
            angchange = numpy.array([])
        
        # rotation
        elif randmove < SurfaceSimulation.percentmove[1]:
            jac = 1
            [a,b,g] = numpy.random.normal(0, self.moveparam[1], 3)
            self.newcoord = moveset.rotation(self.coord, a, b, g)
            self.rotmoves += 1
            movetype = 'rot'
            torschange = numpy.array([])
            angchange = numpy.array([])

            
        # angle bend
        elif randmove < SurfaceSimulation.percentmove[2]:
	    #theta = self.moveparam[2] - 2 * self.moveparam[2] * random.random()
            theta = numpy.random.normal(0, self.moveparam[2])
            self.newcoord, jac = moveset.canglebend(self.coord, m, randdir, theta)
            self.angmoves += 1
            movetype = 'a'
            torschange = numpy.array([])
            angchange = numpy.array([m-1])

        # axis torsion
        elif randmove < SurfaceSimulation.percentmove[3]:
            jac = 1
	    #theta = self.moveparam[3] - 2 * self.moveparam[3] * random.random()
            theta = numpy.random.normal(0, self.moveparam[3])
            self.newcoord = moveset.caxistorsion(self.coord, m, randdir, theta)
            movetype = 'at'
            self.atormoves += 1
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
	elif randmove < SurfaceSimulation.percentmove[4]:
		jac = 1
		theta = numpy.random.normal(0,self.moveparam[4],Simulation.numbeads-2)
		self.newcoord = moveset.cglobalcrank(self.coord,theta)
		self.gcmoves += 1
		movetype = 'gc'
		torschange = numpy.arange(Simulation.numbeads-3)
		angchange = numpy.arange(Simulation.numbeads-2)

	# parrot move
	elif randmove < SurfaceSimulation.percentmove[5]:
	    #theta = self.moveparam[5] - 2 * self.moveparam[5] * random.random()
            theta = numpy.random.normal(0, self.moveparam[5])
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
            self=moveset.runMD_surf(self, 18, .1, dict)
            movetype = 'md'
            self.mdmove += 1
            torschange = numpy.arange(Simulation.numbeads-3)
            angchange = numpy.arange(Simulation.numbeads-2)
            jac = 1
        
        # check boundary conditions
        if not uncloseable:
		# x direction
		max = numpy.max(self.newcoord[:,0])
		min = numpy.min(self.newcoord[:,0])
		if max  > .5*xlength:
			self.newcoord[:,0] -= numpy.ceil((max-xlength*.5)/spacing)*spacing
		elif min < -.5*xlength:
			self.newcoord[:,0] += numpy.ceil((-.5*xlength-min)/spacing)*spacing
		# y direction
		max = numpy.max(self.newcoord[:,1])
		min = numpy.min(self.newcoord[:,1]) 
		if max  > .5*ylength:
			self.newcoord[:,1] -= numpy.ceil((max-ylength*.5)/yspacing)*yspacing

		elif min < -.5*ylength:
			self.newcoord[:,1] += numpy.ceil((-.5*ylength-min)/yspacing)*yspacing

        # accept or reject
        if not uncloseable and movetype=='md':
            self = update_energy(self, torschange, angchange)
            self.newH=self.u1+.5/4.184*numpy.sum(mass*numpy.sum(self.vel**2,axis=1)) # in kcal/mol
            self.move += 1
            boltz = numpy.exp(-(self.newH-self.oldH)/(Simulation.kb*self.T))
            if random.random() < boltz:
                self.accepted += 1
                self.acceptedmd += 1
                self.r2 = self.r2new
                self.surfE = self.newsurfE
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
                if movetype == 'tr':
                    self.acceptedtr += 1
                elif movetype == 'rot':
                    self.acceptedrot += 1
                elif movetype == 'a':
                    self.accepteda += 1
                elif movetype == 'at': 
                    self.acceptedat += 1
                elif movetype == 'gc':
                    self.acceptedgc += 1
		else:
		    self.acceptedp += 1
                self.r2 = self.r2new
                self.surfE = self.newsurfE
                self.coord = self.newcoord
                self.torsE = self.newtorsE
                self.angE = self.newangE
                self.u0 = self.u1
            else:
                self.rejected += 1
        if self.move % Simulation.step == 0:
            self = save(self)
    return self





