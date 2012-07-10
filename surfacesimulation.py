from simulationobject import Simulation
import numpy
import pdb
import random
import energyfunc
import moveset
import writetopdb

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

def energy(rsquare, torsE, angE, surfE):
    """Returns the total energy of a configuration"""
    energy = sum(angE) + sum(torsE) + energyfunc.LJenergy_n(rsquare, Simulation.nativeparam_n,
             Simulation.nonnativeparam, Simulation.nnepsil) + surfE
    return energy

    
class SurfaceSimulation(Simulation):
    kb = 0.0019872041 #kcal/mol/K
    percentmove = [.1, .2, .3, .4, .75, 1] # % translation, % rotation % bend,% axis torsion,% global crankshaft, %ParRot move, %MD

    def __init__(self, name, outputdirectory, coord, temp, surf_coord):
        Simulation.__init__(self, name, outputdirectory, coord, temp)
        self.addsurface(surf_coord)
        self.r2surf = energyfunc.getr2surf(self.coord, surf_coord, Simulation.numbeads, SurfaceSimulation.nspint)
        self.surfE = energyfunc.surfenergy(self.r2surf, SurfaceSimulation.surfparam)
        self.surfE_array = numpy.empty(Simulation.totmoves/Simulation.step + 1)
        self.surfE_array[0] = self.surfE
        self.comx = numpy.empty(Simulation.totmoves/Simulation.step + 1)
        self.comy = numpy.empty(Simulation.totmoves/Simulation.step + 1)
        self.comz = numpy.empty(Simulation.totmoves/Simulation.step + 1)
        self.comx[0] = 0
        self.comy[0] = 0
        self.comz[0] = numpy.sum(self.coord[:,2])/Simulation.numbeads
        
        self.moveparam = numpy.array([5. * self.T / 350 * 50 / Simulation.numbeads, # translation
                        .5 * self.T / 400 * 50 / Simulation.numbeads, # rotation
                        10. * self.T / 250 * 60 / Simulation.numbeads * numpy.pi / 180, # bend
                        10. * self.T / 400 * 60 / Simulation.numbeads * numpy.pi / 180, # torsion
                        3. * self.T / 400 * 60 / Simulation.numbeads * numpy.pi / 180, # global crankshaft
                        5. * self.T / 300 * 60 / Simulation.numbeads * numpy.pi / 180]) # ParRot move
        self.u0 = energy(self.r2, self.torsE, self.angE, self.surfE)
        self.energyarray[0] = self.u0
        self.rmsd_array[0] = 0.0
        self.coord_nat = writetopdb.moviecoord(self.coord, Simulation.transform)
        self.nc[0] = energyfunc.nativecontact(self.r2, Simulation.nativeparam_n, Simulation.nsigma2)

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
        self.coord[:,2] += energyfunc.radgyr(self.coord)*1.5
        
    def output(self, verbose):
        pdb.set_trace()
        write = ['-------- %s Simulation Results --------\r\n' % (self.name),
                'total accepted moves: %d \r\n' % (self.accepted),
                'total rejected moves: %d \r\n' % (self.rejected),
		'ParRot move: %d moves accepted out of %d tries... that is %d percent \r\n' % (self.acceptedp, self.pmoves, float(self.acceptedp)/float(self.pmoves)*100),
		'%d ParRot moves rejected due to chain closure \r\n' % (self.pclosure),
                'angle bend: %d moves accepted out of %d tries... that is %d percent \r\n' % (self.accepteda, self.angmoves, float(self.accepteda)/float(self.angmoves)*100),
		'global crankshaft: %d moves accepted out of %d tries... that is %d percent \r\n' % (self.acceptedgc, self.gcmoves, float(self.acceptedgc)/float(self.gcmoves)*100),
                'torsion: %d moves accepted out of %d tries... that is %d percent \r\n' % (self.acceptedat, self.atormoves, float(self.acceptedat)/float(self.atormoves)*100),
                'MD moves: %d moves accepted out of %d tries... that is 0 percent \r\n' % (self.acceptedmd, self.mdmove)] #, float(self.acceptedmd)/self.mdmove*100)]
        if verbose:
            print "".join(write)
        return write
        
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
    mass = dict['mass']
    surface = dict['surface']
    nspint = dict['nspint']
    surfparam = dict['surfparam']

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
	    theta = self.moveparam[2] - 2 * self.moveparam[2] * random.random()
            self.newcoord, jac = moveset.canglebend(self.coord, m, randdir, theta)
            self.angmoves += 1
            movetype = 'a'
            torschange = numpy.array([])
            angchange = numpy.array([m-1])

        # axis torsion
        elif randmove < SurfaceSimulation.percentmove[3]:
            jac = 1
	    theta = self.moveparam[3] - 2 * self.moveparam[3] * random.random()
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
	    theta = self.moveparam[5] - 2 * self.moveparam[5] * random.random()
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
            self=moveset.runMD_surf(self, 5, .01, dict)
            movetype = 'md'
            self.move += 1
            self.mdmove += 1

        # accept or reject
        if uncloseable == False and movetype == 'md':
            self.r2new = energyfunc.cgetLJr2(self.newcoord, Simulation.numint, Simulation.numbeads)
            self.r2surfnew = energyfunc.getr2surf(self.newcoord, surface, Simulation.numbeads, nspint)
            
            self.newtorsE = energyfunc.ctorsionenergy(self.newcoord, numpy.zeros(Simulation.numbeads-3), Simulation.torsparam, numpy.arange(Simulation.numbeads-3))
            self.newangE = energyfunc.cangleenergy(self.newcoord, numpy.zeros(Simulation.numbeads-2), Simulation.angleparam, numpy.arange(Simulation.numbeads-2))
            self.newsurfE = energyfunc.surfenergy(self.r2surfnew, surfparam)
            
            self.u1 = energy(self.r2new, self.newtorsE, self.newangE, self.newsurfE)
            #self.newH  =self.u1 + .5 / 4.184 * numpy.sum(mass * numpy.sum(self.vel**2, axis=1))
            boltz = numpy.exp(-(self.u1-self.u0)/(Simulation.kb*self.T))
            if random.random() < boltz:
                self.accepted += 1
                self.acceptedmd += 1
                self.r2 = self.r2new
                self.r2surf = self.r2surfnew
                self.surfE = self.newsurfE
                self.coord = self.newcoord
                self.torsE = self.newtorsE
                self.angE = self.newangE
                self.u0 = self.u1
            else:
                self.rejected += 1
        elif(uncloseable == False and movetype != 'md'):
            self.r2new = energyfunc.cgetLJr2(self.newcoord, Simulation.numint, Simulation.numbeads)
            self.r2surfnew = energyfunc.getr2surf(self.newcoord, surface, Simulation.numbeads, nspint)
            
            self.newangE = energyfunc.cangleenergy(self.newcoord, self.angE, Simulation.angleparam, angchange)
            self.newtorsE = energyfunc.ctorsionenergy(self.newcoord, self.torsE, Simulation.torsparam, torschange)
            self.newsurfE = energyfunc.surfenergy(self.r2surfnew, surfparam)
            
            self.u1 = energy(self.r2new, self.newtorsE, self.newangE, self.newsurfE)
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
                self.r2surf = self.r2surfnew
                self.surfE = self.newsurfE
                self.coord = self.newcoord
                self.torsE = self.newtorsE
                self.angE = self.newangE
                self.u0 = self.u1
            else:
                self.rejected += 1
        if self.move % Simulation.step == 0:
            # energy array
            self.energyarray[self.move/Simulation.step] = self.u0
            self.comx[self.move/Simulation.step] = numpy.sum(self.coord[:,0])/Simulation.numbeads
            self.comy[self.move/Simulation.step] = numpy.sum(self.coord[:,1])/Simulation.numbeads 
            self.comz[self.move/Simulation.step] = numpy.sum(self.coord[:,2])/Simulation.numbeads
            self.surfE_array[self.move/Simulation.step] = self.surfE
            
            # native contact array
            self.nc[self.move/Simulation.step] = energyfunc.nativecontact(self.r2, Simulation.nativeparam_n, Simulation.nsigma2)
            # rmsd array
            self.mcoord = writetopdb.moviecoord(self.coord, Simulation.transform)
            self.rmsd_array[self.move/Simulation.step] = energyfunc.rmsd(self.coord_nat, self.mcoord)
            if (Simulation.pdbfile):
                writetopdb.addtopdb(self.coord,Simulation.positiontemplate,self.move/Simulation.step,'%s/trajectory%i.pdb' % (self.out,int(self.T)))
    return self





