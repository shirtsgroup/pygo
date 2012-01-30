from numpy import *
from random import *

def getangleparam(paramfile,numbeads):
    f=open(paramfile,'r')
    param=zeros((numbeads-2,2))
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
    param=zeros((4*(numbeads-3),3)) #4 torsional potentials per 4 molecules
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

def getLJparam(paramfile,numbeads):
    f=open(paramfile,'r')
    param=zeros((numbeads,2)) #two nonzero parameters, ignored first column (all 0)
    while 1:
	line=f.readline()
        if "NONBONDED" in line:
            f.readline() # two lines between header and parameters
	    f.readline()
	    break
    for i in range(numbeads): # param for each bead, use combining rules for interaction
	line=f.readline()
	param[i,0]=float(line[14:23]) #somewhat hardcoded
	param[i,1]=float(line[25:33])
    f.close()
    return param

def getnativefix(paramfile):
    f=open(paramfile,'r')
    param=zeros((0,4)) #will be appending
    while 1:
	line=f.readline()
        if "NBFIX" in line:
            break
    while 1:
	line=f.readline()
	if not line:
		break
	if "G" in line:
		currentparam=[int(line[1:3]),int(line[9:11]),float(line[19:28]),float(line[32:-1])]
		param=vstack((param,currentparam))
    f.close()
    return param

def LJenergy(mpos, param, paramnative):
    energy=0.0
    index=arange(len(mpos))
    for i in index: #from 0 to numbeads-1
    	vdw=index[index>i+2] #index of beads excluding 12 and 13 neighbors in one direction (don't want to overcount)
        for j in vdw:
            k=isNative(i+1,j+1,paramnative) #returns row index of native parameters or -1 if not native
	    r=((mpos[i,0]-mpos[j,0])**2+(mpos[i,1]-mpos[j,1])**2+(mpos[i,2]-mpos[j,2])**2)**.5
	    if k==-1: #not native
		    sigma=(param[i,1]+param[j,1]) #average of the two, but parameter is sigma/2 already
		    epsil=(param[i,0]*param[j,0])**.5
		    energy += epsil*(sigma/r)**12
		    #print(term)
	    else:
		    sigma=paramnative[k,3]
		    epsil=-1*paramnative[k,2]
	    	    energy += epsil*(13*(sigma/r)**12-18*(sigma/r)**10+4*(sigma/r)**6)
    ##calculate repulsive radius term for 12 and 13 neighbor in one direction
    #for i in range(len(mpos)-2):
	#for m in range(2):
		#j=i+m+1
		#k=isNative(i,j,paramnative)
		#if k==-1:
			#epsil=1*(param[i,0]*param[j,0])**.5
			#r=((mpos[i,0]-mpos[j,0])**2+(mpos[i,1]-mpos[j,1])**2+(mpos[i,2]-mpos[j,2])**2)**.5
		    	#term=epsil*((radii[i]+radii[j])*.5/r)**12
			#energy += term
			##print(term)
    ## calculate 12 neighbor for last two beads
    #i=len(mpos)-2
    #j=len(mpos)-1
    #k=isNative(i,j,paramnative)
    #if k==-1:
	#epsil=1*(param[i,0]*param[j,0])**.5
	#r=((mpos[i,0]-mpos[j,0])**2+(mpos[i,1]-mpos[j,1])**2+(mpos[i,2]-mpos[j,2])**2)**.5
	#term=epsil*((radii[i]+radii[j])*.5/r)**12
	#energy += term
    #print('LJ: '+str(energy))
    return energy

def isNative(i,j,nativeloc):
	# checks if i and j beads are native
	# nativeloc: first two columns of paramnative file returned by getnativefix
	# if native, returns index of paramnative
	# if not native, returns -1
	index1=where(nativeloc[:,0]==i)
	index2=where(nativeloc[:,1]==j)
	index1=index1[0][:]
	index2=index2[0][:]
	loc= set(index1) & set(index2)
	index=-1
	for i in loc:
		index=i
	return index

#def calcrepulsiver(mpos,paramnative): 
	## calculates the repulsive radius as defined as the shortest distance to a nonnative residue
	#n=len(mpos)
	#radii=zeros(n)
	#r=10000*ones((n,n));
	#for i in range(n):
		#rmin=10000 #meh, is 10,000 good enough?
		#for j in range(i+1,n):
			#if isNative(i,j,paramnative) == -1:
				#r[i,j]=((mpos[i,0]-mpos[j,0])**2+(mpos[i,1]-mpos[j,1])**2+(mpos[i,2]-mpos[j,2])**2)**.5
				#if r[i,j]<rmin:
					#rmin=r[i,j]
			#else:
				#pass
		#for k in range(i-1,-1,-1):
			#if r[k,i]<rmin:
				#rmin=r[k,i]
		#radii[i]=rmin
	#return radii

# this one is slower
#def calcrepulsiver2(mpos,paramnative): 
	## calculates the repulsive radius as defined as the shortest distance to a nonnative residue
	#n=len(mpos)
	#radii=zeros(n)
	#r=10000*ones((n,n));
	#for i in range(n):
		#rmin=10000 #meh, is 10,000 good enough?
		#for j in range(n):
			#if isNative(i,j,paramnative) != -1 or i==j:
				#pass
			#elif j > i:
				#r[i,j]=((mpos[i,0]-mpos[j,0])**2+(mpos[i,1]-mpos[j,1])**2+(mpos[i,2]-mpos[j,2])**2)**.5
				#if r[i,j] < rmin:
					#rmin = r[i,j]
			#else: #i < j
				#if r[j,i]< rmin: # has already been calculated, note r[j,i]=r[i,j]
					#rmin = r[j,i]
		#radii[i]=rmin
	#return radii

def angleenergy(mpos, param):
    energy=0.0
    ktheta=0.0
    optangle=0.0
    for i in range(1,len(mpos)-1):
        ktheta=param[i-1][0] # param file goes from 0 to len(mpos)-2
	optangle=pi/180*param[i-1,1]
	BA=mpos[i-1,:]-mpos[i,:]
        BC=mpos[i+1,:]-mpos[i,:]
        angle=arccos(dot(BA,BC)/(dot(BA,BA)**.5*dot(BC,BC)**.5)) #in radians
        energy=energy+ktheta*(angle-optangle)**2
    #print('angle energy: '+str(energy))
    return energy

def torsionenergy(mpos, param):
    energy=0.0
    for i in range(0,len(mpos)-3):
        AB=mpos[:][i+1]-mpos[:][i]
        BC=mpos[:][i+2]-mpos[:][i+1]
        CD=mpos[:][i+3]-mpos[:][i+2]
        plane1=cross(BC,AB)
        plane2=cross(CD,BC)
        dihedral=arccos((plane1[0]*plane2[0]+plane1[1]*plane2[1]+plane1[2]*plane2[2])/((plane1[0]**2+plane1[1]**2+plane1[2]**2)**.5*(plane2[0]**2+plane2[1]**2+plane2[2]**2)**.5))
        energy1=param[4*i,0]*(1+cos(param[4*i,1]*dihedral-pi/180*param[4*i,2]))
	energy2=param[4*i+1,0]*(1+cos(param[4*i+1,1]*dihedral-pi/180*param[4*i+1,2]))
	energy3=param[4*i+2,0]*(1+cos(param[4*i+2,1]*dihedral-pi/180*param[4*i+2,2]))
	energy4=param[4*i+3,0]*(1+cos(param[4*i+3,1]*dihedral-pi/180*param[4*i+3,2]))
	energy=energy+energy1+energy2+energy3+energy4
    #print('torsion energy: '+str(energy))
    return energy

#used in simulatepolymer (no amino acid interactions)
#def energy(mpos):
    #energy=0.0 #potential energy
    #sig=4.6 #angstroms for polyethylene
    #e= .42 #kcal/mol for polyethylene
    #ktheta= .82 #kcal/mol
    #A=5.22 #torsional parameter
    #B=2.88 #torsional parameter
    #C=1.95 #torsional parameter
    #index=arange(len(mpos))
    ## 6-12 LJ potential 
    #for i in index:
        #low=index[index<i-2]
        #high=index[index>i+2]
        #vdw=append(low,high) #index of beads excluding 12 and 13 neighbors
        #for j in vdw:
            #r=((mpos[i][0]-mpos[j][0])**2+(mpos[i][1]-mpos[j][1])**2+(mpos[i][2]-mpos[j][2])**2)**.5
            #energy=energy+2*e*((sig/r)**12-(sig/r)**6) #divided by two since count each interaction twice
    ## angle potential
    #for i in range(1,len(mpos)-1):
        #BA=mpos[:][i-1]-mpos[:][i]
        #BC=mpos[:][i+1]-mpos[:][i]
        #angle=arccos(dot(BA,BC)/(dot(BA,BA)**.5*dot(BC,BC)**.5)) #in radians
        #energy=energy+ktheta/2*(angle-pi)**2
    ## torsional potential
    #for i in range(0,len(mpos)-3):
        #AB=mpos[:][i+1]-mpos[:][i]
        #BC=mpos[:][i+2]-mpos[:][i+1]
        #CD=mpos[:][i+3]-mpos[:][i+2]
        #plane1=cross(BC,AB)
        #plane2=cross(CD,BC)
        #dihedral=arccos((plane1[0]*plane2[0]+plane1[1]*plane2[1]+plane1[2]*plane2[2])/((plane1[0]**2+plane1[1]**2+plane1[2]**2)**.5*(plane2[0]**2+plane2[1]**2+plane2[2]**2)**.5))
        #energy=energy+A+B*cos(dihedral)+C*cos(3*dihedral)
    #return energy

