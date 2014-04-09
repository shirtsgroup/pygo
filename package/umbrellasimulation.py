from surfacesimulation import SurfaceSimulation, getsurf, writesurf
from simulationobject import Simulation
import numpy
import pdb
#numpy.random.seed(10)
import energyfunc
import moveset
import writetopdb
import cPickle

class UmbrellaSimulation(SurfaceSimulation):
    def __init__(self, name, outputdirectory, coord, temp, surf_coord, z_pin, mass):
        self.z_pin = z_pin
        SurfaceSimulation.__init__(self, name, outputdirectory, coord, temp, surf_coord)
        self.mass = mass
        self.totmass = numpy.sum(mass)
        self.z_array = numpy.empty(Simulation.totmoves/Simulation.save + 1)
        self.z_array[0] = numpy.sum(mass*self.coord[:,2])/self.totmass
        self.u0 += energyfunc.umbrellaenergy(self.coord, self.z_pin, mass, self.totmass)
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
        self.u0 += energyfunc.umbrellaenergy(self.coord, self.z_pin, self.mass, self.totmass)

    def loadextend(self,extenddirec):
        SurfaceSimulation.loadextend(self,extenddirec)
        self.z_array[0] = numpy.sum(self.mass*self.coord[:,2])/self.totmass
        self.u0 += energyfunc.umbrellaenergy(self.coord, self.z_pin, self.mass, self.totmass)
        self.energyarray[0] = self.u0
    
    def update_energy(self, torschange, angchange, dict):
        SurfaceSimulation.update_energy(self, torschange, angchange, dict)
        self.u1 += energyfunc.umbrellaenergy(self.newcoord, self.z_pin, self.mass, self.totmass)
    
    def save_state(self, dict):
        Simulation.save = dict['save']
        SurfaceSimulation.save_state(self, dict)
        self.z_array[self.move/Simulation.save] = numpy.sum(self.mass*self.coord[:,2])/self.totmass
