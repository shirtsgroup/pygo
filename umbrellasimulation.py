from surfacesimulation import SurfaceSimulation, getsurf, writesurf, run_surf
from simulationobject import Simulation
import numpy
import pdb
#numpy.random.seed(10)
import energyfunc
import moveset
import writetopdb
import cPickle

def update_energy(self, torschange, angchange):
    self.newtorsE = energyfunc.ctorsionenergy(self.newcoord, self.torsE, Simulation.torsparam, torschange)
    self.newangE = energyfunc.cangleenergy(self.newcoord, self.angE, Simulation.angleparam, angchange)
    self.newsurfE = energyfunc.csurfenergy(self.newcoord, SurfaceSimulation.surface, Simulation.numbeads, SurfaceSimulation.nspint, SurfaceSimulation.surfparam,SurfaceSimulation.scale)
    self.r2new, self.u1 = energyfunc.cgetLJenergy(self.newcoord, Simulation.numint, Simulation.numbeads, Simulation.nativeparam, Simulation.nonnativeparam, Simulation.nnepsil)
    self.u1 += numpy.sum(self.newtorsE)+numpy.sum(self.newangE)+numpy.sum(self.newsurfE)
    self.u1 += energyfunc.umbrellaenergy(self.newcoord, self.z_pin, self.mass, self.totmass, self.k)
    return self

def save(self):
    index = self.move/Simulation.save
    self.energyarray[index] = self.u0
    self.surfE_array[index,:] = self.surfE
    self.nc[index] = energyfunc.nativecontact(self.r2, Simulation.nativeparam, Simulation.nsigma2)
    self.z_array[index] = numpy.sum(self.mass*self.coord[:,2])/self.totmass
    if (Simulation.writetraj):
        f = open('%s/trajectory%i' %(self.out, int(self.T)), 'ab')
        numpy.save(f,self.coord)
        f.close()
    return self

class UmbrellaSimulation(SurfaceSimulation):
    def __init__(self, name, outputdirectory, coord, temp, surf_coord, z_pin, mass, k):
        self.z_pin = z_pin
        self.k = k
        SurfaceSimulation.__init__(self, name, outputdirectory, coord, temp, surf_coord)
        self.mass = mass
        self.totmass = numpy.sum(mass)
        self.z_array = numpy.empty(Simulation.totmoves/Simulation.save + 1)
        self.z_array[0] = numpy.sum(mass*self.coord[:,2])/self.totmass
        self.u0 += energyfunc.umbrellaenergy(self.coord, self.z_pin, mass, self.totmass, self.k)
        self.energyarray[0] = self.u0

    def addsurface(self, surf_coord):
        self.surface = surf_coord
        self.coord = moveset.rotation(self.coord,numpy.random.random(3)) # randomly orient protein
        self.coord[:,2] += self.z_pin - numpy.min(self.coord[:,2]) # protein is minimum z_pin distance from surface

    def save_z(self):
        filename = '%s/z_traj%i_%i' % (self.out, int(self.z_pin), int(self.T))
        numpy.save(filename, self.z_array)

    def loadstate(self):
        SurfaceSimulation.loadstate(self)
        self.z_array = numpy.load('%s/z_traj%i_%i.npy' %(self.out, int(self.z_pin), int(self.T)))
        self.u0 += energyfunc.umbrellaenergy(self.coord, self.z_pin, self.mass, self.totmass, self.k)

    def loadextend(self,extenddirec):
        SurfaceSimulation.loadextend(self,extenddirec)
        self.z_array[0] = numpy.sum(self.mass*self.coord[:,2])/self.totmass
        self.u0 += energyfunc.umbrellaenergy(self.coord, self.z_pin, self.mass, self.totmass, self.k)
        self.energyarray[0] = self.u0
	
