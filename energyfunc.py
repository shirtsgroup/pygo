from numpy import *
import numpy
from random import *


def getangleparam(paramfile,numbeads):
    f=open(paramfile,'r')
    param=empty((numbeads-2,2))
    while 1:
	line=f.readline()
        if "ANGLE" in line:
            break
    for i in range(numbeads-2): #no angle formed for end beads
	line=f.readline()
	param[i][0]=float(line[26:36]) #somewhat hardcoded, may apply to all/most go model param files?
	param[i][1]=float(line[37:-1])
    f.close()
    return param

def gettorsionparam(paramfile, numbeads):
    f=open(paramfile,'r')
    param=empty((4*(numbeads-3),3)) #4 torsional potentials per 4 molecules
    while 1:
	line=f.readline()
        if "DIHEDRAL" in line:
            break
    for i in range(4*(numbeads-3)):
	line=f.readline()
	param[i][0]=float(line[22:30]) #somewhat hardcoded
	param[i][1]=float(line[32:33])
	param[i][2]=float(line[35:45])
    f.close()
    return param

#def gettorsionparam_n(paramfile, numbeads):
    #f=open(paramfile,'r')
    #param=empty((numbeads-3,12)) #4 torsional potentials per 4 molecules
    #while 1:
	#line=f.readline()
        #if "DIHEDRAL" in line:
            #break
    #for i in range(numbeads-3):
	#line1=f.readline()
	#line2=f.readline()
	#line3=f.readline()
	#line4=f.readline()
	#param[i][0]=float(line1[22:30]) #somewhat hardcoded
	#param[i][1]=float(line1[32:33])
	#param[i][2]=float(line1[35:45])
	#param[i][3]=float(line2[22:30]) #somewhat hardcoded
	#param[i][4]=float(line2[32:33])
	#param[i][5]=float(line2[35:45])
	#param[i][6]=float(line3[22:30]) #somewhat hardcoded
	#param[i][7]=float(line3[32:33])
	#param[i][8]=float(line3[35:45])
	#param[i][9]=float(line4[22:30]) #somewhat hardcoded
	#param[i][10]=float(line4[32:33])
	#param[i][11]=float(line4[35:45])    
    #f.close()
    
    #return param


#def getLJparam(paramfile,numbeads):
    #f=open(paramfile,'r')
    #param=empty((numbeads,2)) #two nonzero parameters, ignored first column (all 0)
    #while 1:
	#line=f.readline()
        #if "NONBONDED" in line:
            #f.readline() # two lines between header and parameters
	    #f.readline()
	    #break
    #for i in range(numbeads): # param for each bead, use combining rules for interaction
	#line=f.readline()
	#param[i,0]=float(line[14:23]) #somewhat hardcoded
	#param[i,1]=float(line[25:33])
    #f.close()
    #return param

#def getnativefix(paramfile):
    #f=open(paramfile,'r')
    #param=zeros((0,4)) #will be appending
    #while 1:
	#line=f.readline()
        #if "NBFIX" in line:
            #break
    #while 1:
	#line=f.readline()
	#if not line:
		#break
	#if "G" in line:
		#currentparam=[int(line[1:3]),int(line[9:11]),float(line[19:28]),float(line[32:-1])]
		#param=vstack((param,currentparam))
    #f.close()
    #return param
	
# speed up version
def getLJparam_n(paramfile,numbeads,numint):
	f=open(paramfile,'r')
    	param=empty(numbeads) #two nonzero parameters, ignored first column (all 0)
    	while 1:
		line=f.readline()
        	if "NONBONDED" in line:
            		f.readline() # two lines between header and parameters
	    		f.readline()
	    		break
    	for i in range(numbeads): # param for each bead, use combining rules for interaction
		line=f.readline()
		epsil=float(line[14:23]) #somewhat hardcoded
		param[i]=float(line[25:33])
    	f.close()
	sigarray=zeros(numint)
	index=arange(numbeads)
	k=0
	for i in index:
		vdw=index[index>i+2]
		for j in vdw:
			sigarray[k]=param[i]+param[j]
			k += 1
	sigarray=append(sigarray,epsil)
    	return sigarray

#speed up version
def getnativefix_n(paramfile,numint,numbeads):
	param=zeros((numint,3))
	f=open(paramfile,'r')
    	while 1:
		line=f.readline()
        	if "NBFIX" in line:
            		break
    	while 1:
		line=f.readline()
		if not line:
			break
		if "G" in line:
			[i,j,ep,sig]=[int(line[1:3]),int(line[9:11]),float(line[19:28]),float(line[32:-1])]
			intindex=sum(numbeads-arange(1,i))+(j-i)-1-2*i
			param[intindex,:]=[1,-ep,sig]
	return param

#speed up version
def getLJr2(mpos,numint,numbeads):
	r2array=empty(numint)
	k=0
	for i in range(numbeads):
		BC=mpos[i,:]-mpos[i+3:numbeads,:]
		knew=k+numbeads-(i+3)
		r2array[k:knew]=sum(BC**2,axis=1)
		k=knew
	return r2array #r^2 values for every interaction



#speed up version
def LJenergy_n(r2,natparam,nonnatparam,nnepsil):
	#native calculation
	nE=natparam[:,0]*natparam[:,2]*natparam[:,2]/r2 #sigma2/r2
	nE6=nE*nE*nE
	nE=natparam[:,1]*(13*nE6*nE6-18*nE6*nE*nE+4*nE6)
	#nonnative calculation
	nnE=nonnatparam[:,0]*nonnatparam[:,1]*nonnatparam[:,1]/r2 #simga2/r2
	nnE=nnE*nnE #sigma4/r4
	nnE=nnepsil*nnE*nnE*nnE
	energy=sum(nE)+sum(nnE)
	return energy

#def LJenergy(mpos, param, paramnative):
    #energy=0.0
    #index=arange(len(mpos))
    #for i in index: #from 0 to numbeads-1
    	#vdw=index[index>i+2] #index of beads excluding 12 and 13 neighbors in one direction (don't want to overcount)
        #for j in vdw:
            #k=isNative(i+1,j+1,paramnative) #returns row index of native parameters or -1 if not native
	    #r=((mpos[i,0]-mpos[j,0])**2+(mpos[i,1]-mpos[j,1])**2+(mpos[i,2]-mpos[j,2])**2)**.5
	    #if k==-1: #not native
		    #sigma=(param[i,1]+param[j,1]) #average of the two, but parameter is sigma/2 already
		    #epsil=(param[i,0]*param[j,0])**.5
		    #energy += epsil*(sigma/r)**12
		    ##print(term)
	    #else:
		    #sigma=paramnative[k,3]
		    #epsil=-1*paramnative[k,2]
	    	    #energy += epsil*(13*(sigma/r)**12-18*(sigma/r)**10+4*(sigma/r)**6)

    #return energy

#def isNative(i,j,nativeloc):
	## checks if i and j beads are native
	## nativeloc: first two columns of paramnative file returned by getnativefix
	## if native, returns index of paramnative
	## if not native, returns -1
	#index1=where(nativeloc[:,0]==i)
	#index2=where(nativeloc[:,1]==j)
	#index1=index1[0][:]
	#index2=index2[0][:]
	#loc= set(index1) & set(index2)
	#index=-1
	#for i in loc:
		#index=i
	#return index


#def angleenergy(mpos, param):
    #energy=0.0
    #for i in range(1,len(mpos)-1):
        #ktheta=param[i-1,0] # param file goes from 0 to len(mpos)-2
	#optangle=pi/180*param[i-1,1]
	#BA=mpos[i-1,:]-mpos[i,:]
        #BC=mpos[i+1,:]-mpos[i,:]
        #angle=arccos(dot(BA,BC)/(dot(BA,BA)**.5*dot(BC,BC)**.5)) #in radians
        #energy+=ktheta*(angle-optangle)**2
    ##print('angle energy: '+str(energy))
    #return energy
    
def angleenergy_n(mpos, oldE,param,change):
    newE=oldE.copy()
    for i in change:
        ktheta=param[i,0] # param file goes from 0 to len(mpos)-2
	optangle=pi/180*param[i,1]
	BA=mpos[i,:]-mpos[i+1,:]
        BC=mpos[i+2,:]-mpos[i+1,:]
	dotBABC=BA[0]*BC[0]+BA[1]*BC[1]+BA[2]*BC[2]
	dotBA=BA[0]*BA[0]+BA[1]*BA[1]+BA[2]*BA[2]
	dotBC=BC[0]*BC[0]+BC[1]*BC[1]+BC[2]*BC[2]
        angle=arccos(dotBABC/(dotBA**.5*dotBC**.5)) #in radians
        newE[i]=ktheta*(angle-optangle)*(angle-optangle)
    #print('angle energy: '+str(energy))
    return newE


def angle(mpos):
	angle=zeros(len(mpos)-2)
	for i in range(1,len(mpos)-1):
		BA=mpos[i-1,:]-mpos[i,:]
		BC=mpos[i+1,:]-mpos[i,:]
		angle[i-1]=arccos(dot(BA,BC)/(dot(BA,BA)**.5*dot(BC,BC)**.5))
	return angle

#def torsionenergy(mpos, param):
    #energy=0.0
    #for i in range(0,len(mpos)-3):
        #AB=mpos[:][i+1]-mpos[:][i]
        #BC=mpos[:][i+2]-mpos[:][i+1]
        #CD=mpos[:][i+3]-mpos[:][i+2]
        #plane1=cross(BC,AB)
        #plane2=cross(CD,BC)
        #dihedral=arccos((plane1[0]*plane2[0]+plane1[1]*plane2[1]+plane1[2]*plane2[2])/((plane1[0]**2+plane1[1]**2+plane1[2]**2)**.5*(plane2[0]**2+plane2[1]**2+plane2[2]**2)**.5))
        #energy1=param[4*i,0]*(1+cos(param[4*i,1]*dihedral-pi/180*param[4*i,2]))
	#energy2=param[4*i+1,0]*(1+cos(param[4*i+1,1]*dihedral-pi/180*param[4*i+1,2]))
	#energy3=param[4*i+2,0]*(1+cos(param[4*i+2,1]*dihedral-pi/180*param[4*i+2,2]))
	#energy4=param[4*i+3,0]*(1+cos(param[4*i+3,1]*dihedral-pi/180*param[4*i+3,2]))
	#energy += energy1+energy2+energy3+energy4
    ##print('torsion energy: '+str(energy))
    #return energy

def dihedral(mpos,param):
	#newdihed=olddihed.copy()
	#newE=oldE.copy()
	newdihed=zeros(len(mpos)-3)
	for i in range(len(mpos)-3):
		AB=mpos[:][i+1]-mpos[:][i]
        	BC=mpos[:][i+2]-mpos[:][i+1]
        	CD=mpos[:][i+3]-mpos[:][i+2]
        	plane1=cross(BC,AB)
        	plane2=cross(CD,BC)
        	newdihed[i]=arccos((plane1[0]*plane2[0]+plane1[1]*plane2[1]+plane1[2]*plane2[2])/((plane1[0]**2+plane1[1]**2+plane1[2]**2)**.5*(plane2[0]**2+plane2[1]**2+plane2[2]**2)**.5))
		
	return newdihed

#def torsionenergy_n(mpos,oldE,param,change):
    #newE=oldE.copy()
    #for i in change:
        #AB=mpos[:][i+1]-mpos[:][i]
        #BC=mpos[:][i+2]-mpos[:][i+1]
        #CD=mpos[:][i+3]-mpos[:][i+2]
        #plane1=cross(BC,AB)
        #plane2=cross(CD,BC)
        #dihedral=arccos((plane1[0]*plane2[0]+plane1[1]*plane2[1]+plane1[2]*plane2[2])/((plane1[0]**2+plane1[1]**2+plane1[2]**2)**.5*(plane2[0]**2+plane2[1]**2+plane2[2]**2)**.5))
        #energy1=param[4*i,0]*(1+cos(param[4*i,1]*dihedral-pi/180*param[4*i,2]))
	#energy2=param[4*i+1,0]*(1+cos(param[4*i+1,1]*dihedral-pi/180*param[4*i+1,2]))
	#energy3=param[4*i+2,0]*(1+cos(param[4*i+2,1]*dihedral-pi/180*param[4*i+2,2]))
	#energy4=param[4*i+3,0]*(1+cos(param[4*i+3,1]*dihedral-pi/180*param[4*i+3,2]))
	#newE[i]=energy1+energy2+energy3+energy4
    #return newE
    
def torsionenergy_nn(mpos,oldE,param,change):
    newE=oldE.copy()
    for i in change:
        AB=mpos[:][i+1]-mpos[:][i]
        BC=mpos[:][i+2]-mpos[:][i+1]
        CD=mpos[:][i+3]-mpos[:][i+2]
        plane1=[BC[1]*AB[2]-BC[2]*AB[1],BC[2]*AB[0]-BC[0]*AB[2],BC[0]*AB[1]-BC[1]*AB[0]] #cross(BC,AB)
        plane2=[CD[1]*BC[2]-CD[2]*BC[1],CD[2]*BC[0]-CD[0]*BC[2],CD[0]*BC[1]-CD[1]*BC[0]]#cross(CD,BC)
        dotplane1=plane1[0]*plane1[0]+plane1[1]*plane1[1]+plane1[2]*plane1[2]
	dotplane2=plane2[0]*plane2[0]+plane2[1]*plane2[1]+plane2[2]*plane2[2]
	dihedral=arccos((plane1[0]*plane2[0]+plane1[1]*plane2[1]+plane1[2]*plane2[2])/(dotplane1*dotplane2)**.5)
        energy=param[4*i:4*i+4,0]*(1+cos(param[4*i:4*i+4,1]*dihedral-pi/180*param[4*i:4*i+4,2]))
	newE[i]=sum(energy)
    return newE

#used in simulatepolymer (no amino acid interactions)
def polymerenergy(mpos):
    energy=0.0 #potential energy
    sig=4.6 #angstroms for polyethylene
    e= .42 #kcal/mol for polyethylene
    ktheta= .82 #kcal/mol
    A=5.22 #torsional parameter
    B=2.88 #torsional parameter
    C=1.95 #torsional parameter
    index=arange(len(mpos))
    # 6-12 LJ potential 
    for i in index:
        low=index[index<i-2]
        high=index[index>i+2]
        vdw=append(low,high) #index of beads excluding 12 and 13 neighbors
        for j in vdw:
            r=((mpos[i][0]-mpos[j][0])**2+(mpos[i][1]-mpos[j][1])**2+(mpos[i][2]-mpos[j][2])**2)**.5
            energy=energy+2*e*((sig/r)**12-(sig/r)**6) #divided by two since count each interaction twice
    # angle potential
    for i in range(1,len(mpos)-1):
        BA=mpos[:][i-1]-mpos[:][i]
        BC=mpos[:][i+1]-mpos[:][i]
        angle=arccos(dot(BA,BC)/(dot(BA,BA)**.5*dot(BC,BC)**.5)) #in radians
        energy=energy+ktheta/2*(angle-pi)**2
    # torsional potential
    for i in range(0,len(mpos)-3):
        AB=mpos[:][i+1]-mpos[:][i]
        BC=mpos[:][i+2]-mpos[:][i+1]
        CD=mpos[:][i+3]-mpos[:][i+2]
        plane1=cross(BC,AB)
        plane2=cross(CD,BC)
        dihedral=arccos((plane1[0]*plane2[0]+plane1[1]*plane2[1]+plane1[2]*plane2[2])/((plane1[0]**2+plane1[1]**2+plane1[2]**2)**.5*(plane2[0]**2+plane2[1]**2+plane2[2]**2)**.5))
        energy=energy+A+B*cos(dihedral)+C*cos(3*dihedral)
    return energy
    
def rmsd(crds1, crds2):
  	"""Returns RMSD between 2 sets of [nx3] numpy array"""
 	assert(crds1.shape[1] == 3)
 	assert(crds1.shape == crds2.shape)
 	n_vec = numpy.shape(crds1)[0]
 	correlation_matrix = numpy.dot(numpy.transpose(crds1), crds2)
 	v, s, w_tr = numpy.linalg.svd(correlation_matrix)
 	is_reflection = (numpy.linalg.det(v) * numpy.linalg.det(w_tr)) < 0.0
 	if is_reflection:
  		s[-1] = - s[-1]
  	E0 = sum(sum(crds1 * crds1)) + sum(sum(crds2 * crds2))
  	rmsd_sq = (E0 - 2.0*sum(s)) / float(n_vec)
  	rmsd_sq = max([rmsd_sq, 0.0])
 	return numpy.sqrt(rmsd_sq)


