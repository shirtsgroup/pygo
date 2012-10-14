from surfacesimulation import SurfaceSimulation, getsurf, writesurf, run_surf
#from simulationobject import Simulation
import numpy
import pdb
import random
import energyfunc
import moveset
import writetopdb
try:
	import matplotlib.pyplot as plt
except: pass

def update_energy(self, torschange, angchange):
    self.newtorsE = energyfunc.ctorsionenergy(self.newcoord, self.torsE, Simulation.torsparam, torschange)
    self.newangE = energyfunc.cangleenergy(self.newcoord, self.angE, Simulation.angleparam, angchange)
    self.newsurfE = energyfunc.csurfenergy(self.newcoord, SurfaceSimulation.surface, Simulation.numbeads, SurfaceSimulation.nspint, SurfaceSimulation.surfparam)
    self.r2new, self.u1 = energyfunc.cgetLJenergy(self.newcoord, Simulation.numint, Simulation.numbeads, Simulation.nativeparam_n, Simulation.nonnativeparam, Simulation.nnepsil)
    self.u1 += numpy.sum(self.newtorsE)+numpy.sum(self.newangE)+self.newsurfE
    self.u1 += energyfunc.umbrellaenergy(self.newcoord, self.z_pin, self.mass, self.totmass)
    return self

def save(self):
    index = self.move/Simulation.step
    self.energyarray[index] = self.u0
    self.surfE_array[index] = self.surfE
    self.nc[index] = energyfunc.nativecontact(self.r2, Simulation.nativeparam_n, Simulation.nsigma2)
    #self.radgyr[self.move/Simulation.step] = energyfunc.radgyr(self.coord)
    self.mcoord = writetopdb.moviecoord(self.coord, Simulation.transform)
    self.rmsd_array[index] = energyfunc.rmsd(self.coord_nat, self.mcoord)
    self.z_array[index] = numpy.sum(self.mass*self.coord[:,2])/Simulation.numbeads
    if (Simulation.pdbfile):
        writetopdb.addtopdb(self.coord,Simulation.positiontemplate,index,'%s/trajectory%i.pdb' % (self.out,int(self.T)))
    return self

class UmbrellaSimulation(SurfaceSimulation):
    def __init__(self, name, outputdirectory, coord, temp, surf_coord, z_pin, mass):
	SurfaceSimulation.__init__(self, name, outputdirectory, coord, temp, surf_coord)
	self.z_pin = z_pin
	self.mass = mass
	self.totmass = numpy.sum(mass)
	self.z_array = numpy.empty(Simulation.totmovces/Simulation.step + 1)
	self.z_array[0] = numpy.sum(mass*self.coord[:,2])/Simulation.numbeads
	self.u0 += energyfunc.umbrellaenergy(self.coord, self.z_pin, mass, self.totmass)
	self.energyarray[0] = self.u0
    def save_z(self):
	filename = '%s/z_traj%i_%i.txt' % (self.out, int(self.z_pin), int(self.T))
	numpy.savetxt(filename, self.z_array)


