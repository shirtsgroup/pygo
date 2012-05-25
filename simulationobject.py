from numpy import *
from writetopdb import *
from moveset import *
from energyfunc import *
import matplotlib.pyplot as plt
from sys import stdout
from random import *
import random
import profile
import scipy.misc
import numpy
import writetopdb
import moveset
import energyfunc
import matplotlib.pyplot as plt
import pdb

#def energyprint(mpos,rsquare,torsE,angE):
	#LJ=energyfunc.LJenergy_n(rsquare,Simulation.nativeparam_n,Simulation.nonnativeparam,Simulation.nnepsil)
	#energy=sum(angE)+sum(torsE)+LJ
	#print 'Angle energy: ' + str(sum(angE))
	#print 'Torsion energy: '+str(sum(torsE))
	#print 'LJ energy: '+str(LJ)
	#print 'Total energy: '+str(energy)
	#return energy

def energy(mpos,rsquare,torsE,angE):
	energy=sum(angE)+sum(torsE)+energyfunc.LJenergy_n(rsquare,Simulation.nativeparam_n,Simulation.nonnativeparam,Simulation.nnepsil)
	return energy


class Simulation:
	kb=0.0019872041 #kcal/mol/K
	percentmove=[.33,.66] #% bend,% axis torsion,% crankshaft
	
	def __init__(self,name,outputdirectory,coord,temp):
		self.name=name
		self.out=outputdirectory
		self.coord=coord
		self.energyarray=numpy.empty(Simulation.totmoves/Simulation.step+1)
		self.nc=numpy.empty(Simulation.totmoves/Simulation.step+1) #native contacts
		self.rmsd_array=numpy.empty(Simulation.totmoves/Simulation.step+1)
		self.move=0
		self.T=temp
		self.r2=energyfunc.getLJr2(self.coord,Simulation.numint,Simulation.numbeads)
		self.torsE=energyfunc.torsionenergy_nn(self.coord,numpy.zeros(Simulation.numbeads-3),Simulation.torsparam,numpy.arange(Simulation.numbeads-3))
		self.angE=energyfunc.angleenergy_n(self.coord,numpy.zeros(Simulation.numbeads-2),Simulation.angleparam,numpy.arange(Simulation.numbeads-2))
		self.u0=energy(self.coord,self.r2,self.torsE,self.angE)
		self.energyarray[0]=self.u0
		self.rmsd_array[0]=energyfunc.rmsd(Simulation.coord_nat,Simulation.coord_nat)
		self.nc[0]=energyfunc.nativecontact(self.r2,Simulation.nativeparam_n,Simulation.nsigma2)
		self.maxtheta=[10*self.T/300.,10*self.T/200.,20*self.T/250.] # bend, axistorsion, crankshaft
		# constants for move stats
		self.angmoves=0
		self.crankmoves=0
		self.atormoves=0
		self.acceptedat=0
		self.accepteda=0
		self.acceptedc=0
		self.accepted=0
		self.rejected=0
		self.movetype=''
		self.whoami=[]
	
	
	
	def output(self,verbose):
		write='-------- %s Simulation Results --------\r\n' %(self.name) 
		write += 'total accepted moves: %d \r\n' %(self.accepted)
		write += 'total rejected moves: %d \r\n' %(self.rejected)
		#write += 'angle bend: %d moves accepted out of %d tries \r\n' %(self.accepteda,self.angmoves)
		#write += 'crankshaft: %d moves accepted out of %d tries \r\n' %(self.acceptedc,self.crankmoves)
		#write += 'torsion: %d moves accepted out of %d tries \r\n' %(self.acceptedat,self.atormoves)
		if(verbose):
			print write
		return write
	
	def saveenergy(self,plot):
		filename=self.out+'/energy'+str(int(self.T))+'.txt'
		savetxt(filename,self.energyarray)
		#print 'wrote every %d conformation energies to %s' %(Simulation.step,filename)
		if(plot):
			plotname=self.out+'/energy'+str(int(self.T))+'.png'
			#print 'generating conformational energy plot...'
			fig=plt.figure(1)
			plt.plot(range(len(self.energyarray)),self.energyarray,label=str('%3.2f') % self.T +' K')
			plt.xlabel('move/%d' %(Simulation.step))
			plt.ylabel('energy (kcal/mol)')
			plt.title('Go-like model monte carlo simulation at %d K' %(self.T))
			plt.grid(True)
			lgd=plt.legend(loc=2,prop={'size':8})
			fig.savefig(plotname)
			#print 'conformational energy plot saved to %s' %(plotname)
	
	def savermsd(self,plot):
		filename=self.out+'/rmsd'+str(int(self.T))+'.txt'
		savetxt(filename,self.rmsd_array)
		#print 'wrote every %d rmsd values to %s' %(Simulation.step,filename)
		if(plot):
			plotname=self.out+'/rmsd'+str(int(self.T))+'.png'
			#print 'generating rmsd plot...'
			fig=plt.figure(2)
			plt.plot(range(len(self.rmsd_array)),self.rmsd_array,label=str('%3.2f') % self.T +' K')
			plt.xlabel('move/%d' %(Simulation.step))
			plt.ylabel('rmsd')
			plt.title('Go-like model monte carlo simulation at %d K' %(self.T))
			plt.grid(True)
			lgd=plt.legend(loc=2,prop={'size':8})
			fig.savefig(plotname)
			#print 'rmsd plot saved to %s' %(plotname)
		
	
	def savenc(self,plot):
		fraction=self.nc/Simulation.totnc
		fractionfile=self.out+'/fractionnative'+str(int(self.T))+'.txt'
		savetxt(fractionfile,fraction)
		#print 'wrote every %d fractional nativeness values to %s' %(Simulation.step,fractionfile)
		print 'excluding first '+str(len(fraction)/5)+' out of ' + str(len(fraction))+' fractional nativeness values from average'
		fractioncut=fraction[len(fraction)/5:-1]
		print 'the average fractional nativeness is '+str(sum(fractioncut)/len(fractioncut))
		if(plot):
			plotname=self.out+'/fractionnative'+str(int(self.T))+'.png'
			#print 'generating fractional nativeness plot...'
			fig=plt.figure(3)
			plt.plot(range(len(fraction)),fraction,label=str('%3.2f') % self.T +' K')
			plt.xlabel('move/%d' %(Simulation.step))
			plt.ylabel('Q, fractional nativeness')
			plt.title('Go-like model monte carlo simulation at %d K' %(self.T))
			plt.grid(True)
			lgd=plt.legend(loc=3,prop={'size':8})
			fig.savefig(plotname)
			#print 'fractional nativeness plot saved to %s' %(plotname)
	def savehist(self,plot):
		if(plot):
			plotname=self.out+'/energyhist'+str(int(self.T))+'.png'
			#print 'generating energy histogram...'
			fig=plt.figure(4)
			plt.hist(self.energyarray,40,label=str('%3.2f') % self.T +' K')
			plt.xlabel('energy (kcal/mol)')
			plt.title('Go-like model MC simulation at %d K' %(self.T))
			plt.legend(loc=2,prop={'size':8})
			fig.savefig(plotname)
			#print 'energy histogram saved to %s' %(plotname)

		
def run(self,nummoves,dict):
	#,numbeads,step,totmoves,numint,angleparam,torsparam,nativeparam_n,nonnativeparam,nnepsil,nsigma2,transform,coord_nat):
	Simulation.numbeads=dict['numbeads']
	Simulation.step=dict['step']
	Simulation.totmoves=dict['totmoves']
	Simulation.numint=dict['numint']
	Simulation.angleparam=dict['angleparam']
	Simulation.torsparam=dict['torsparam']
	Simulation.nativeparam_n=dict['nativeparam_n']
	Simulation.nonnativeparam=dict['nonnativeparam']
	Simulation.nnepsil=dict['nnepsil']
	Simulation.nsigma2=dict['nsigma2']
	Simulation.transform=dict['transform']
	Simulation.coord_nat=dict['coord_nat']
	for i in range(nummoves):
		self.randmove=random.random()
		self.randdir=random.random()
		self.m=random.randint(1,Simulation.numbeads-2) #random bead, not end ones
		
		#bend
		if self.randmove < Simulation.percentmove[0]:
			self.theta=self.maxtheta[0]/180.*numpy.pi-random.random()*self.maxtheta[0]*numpy.pi/180.*2
			self.newcoord=moveset.anglebend(self.coord,self.m,self.randdir,self.theta)
			self.angmoves += 1
			self.movetype='a'
			self.change=[]
			self.angchange=[self.m-1]
		
		#axis torsion
		elif self.randmove < Simulation.percentmove[1]:
			self.theta=self.maxtheta[1]/180.*numpy.pi-random.random()*numpy.pi*self.maxtheta[1]/180.*2
			self.newcoord=moveset.axistorsion(self.coord,self.m,self.randdir,self.theta)
			self.movetype='at'
			self.atormoves += 1
			self.angchange=[]
			if self.randdir<.5:
				self.change=[self.m-2]
				if self.m<2:
					self.change=[]
			elif self.m==Simulation.numbeads-2:
				self.change=[]
			else:
				self.change=[self.m-1]
			
		#crankshaft
		else:
			self.theta=self.maxtheta[2]/180.*numpy.pi-random.random()*self.maxtheta[2]*numpy.pi/180.*2
			self.newcoord=moveset.crankshaft(self.coord,self.m,self.theta)
			self.crankmoves += 1
			self.movetype='c'
			self.change=numpy.arange(self.m-3,self.m+1)
			self.angchange=[self.m-2,self.m]
			if self.m==2:
				self.change=numpy.arange(self.m-2,self.m+1)
			elif self.m<2:
				self.change=self.change[self.change>-1]
				self.angchange=[self.m]
			elif self.m==Simulation.numbeads-4 or self.m==Simulation.numbeads-3:
				self.change=self.change[self.change<(Simulation.numbeads-3)]
			elif self.m==Simulation.numbeads-2:
				self.change=numpy.arange(self.m-3,self.m-1)
				self.angchange=[self.m-2]
		self.r2new=energyfunc.getLJr2(self.newcoord,Simulation.numint,Simulation.numbeads)
		self.newangE=energyfunc.angleenergy_n(self.newcoord,self.angE,Simulation.angleparam,self.angchange)
		self.newtorsE=energyfunc.torsionenergy_nn(self.newcoord,self.torsE,Simulation.torsparam,self.change)
		self.u1=energy(self.newcoord,self.r2new,self.newtorsE,self.newangE)
		self.move += 1
		#stdout.write(str(self.move)+'\r')
		#stdout.flush()
		self.boltz=numpy.exp(-(self.u1-self.u0)/(Simulation.kb*self.T))
		#stdout.write(str(self.boltz)+'\r')
		#stdout.flush()
		if self.u1< self.u0 or random.random() < self.boltz:
			self.accepted += 1
			if self.movetype=='a':
				self.accepteda += 1
			elif self.movetype=='r':
				self.acceptedr += 1
			elif self.movetype=='at':
				self.acceptedat +=1
				#accepted_angle.append(theta)
			else: 
				self.acceptedc += 1
			#if (pdbfile != ''):
				#mcoord=moviecoord(newcoord,transform)
				##writeseqpdb(mcoord,wordtemplate,ATOMlinenum,accepted)
				#addtopdb(mcoord,positiontemplate,move,pdbfile)
			self.r2=self.r2new
			self.coord=self.newcoord
			self.torsE=self.newtorsE
			self.angE=self.newangE
			self.u0=self.u1
		else:
			self.rejected += 1
		if self.move%Simulation.step==0:
			# energy array
			self.energyarray[self.move/Simulation.step]=self.u0
			
			# native contact array
			self.nc[self.move/Simulation.step]=energyfunc.nativecontact(self.r2,Simulation.nativeparam_n,Simulation.nsigma2)
				
			# rmsd array
			self.mcoord=writetopdb.moviecoord(self.coord,Simulation.transform)
			self.rmsd_array[self.move/Simulation.step]=energyfunc.rmsd(Simulation.coord_nat,self.mcoord)
	#print str(self.name)+' '+str(self.u0)
	return self
	#return [self.energyarray,self.nc,self.rmsd_array,self.coord,self.move,self.accepted,self.rejected,self.acceptedc,self.acceptedat,self.accepteda,self.u0,self.r2,self.angE,self.torsE,self.crankmoves,self.atormoves,self.angmoves]
	
		
	
	