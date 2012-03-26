from numpy import *
from random import *
import pdb

def crankshaft(mpos123,m,theta): #moleculeposition, will need bonds later
    mpos=mpos123.copy()

    AB=mpos[m,:]-mpos[m-1,:]
    AC=mpos[m+1,:]-mpos[m-1,:]

    x1=AC/dot(AC,AC)**.5
    y1=AB-dot(AB,AC)/dot(AC,AC)*AC
    y1=y1/dot(y1,y1)**.5
    z1=cross(x1,y1)
    transform=array([x1,y1,z1])
    untransform=transpose(transform)

    rotate=array([[1,0,0],[0,cos(theta),sin(theta)],[0,-sin(theta),cos(theta)]])

    AB1=dot(transform,AB.transpose())
    AB2=dot(rotate,AB1)
    AB=dot(untransform,AB2)
    mpos[m,:]=mpos[m-1,:]+AB.transpose()
    return mpos

#def reptation(mpos123):
    #mpos=mpos123.copy()
    #theta=arccos(1-2*random())
    #phi=2*pi*random()
    #rho=3.86470703286 #unhardcode this
    #vec=[rho*sin(theta)*cos(phi),rho*sin(theta)*sin(phi),rho*cos(theta)]    
    #n=len(mpos)
    #if random() < .5:    
	#for i in range(n-1):
             #mpos[i,:]=mpos[i+1,:]
    	#mpos[n-1,:]=mpos[n-2][:]+vec
    #else:
	#for i in range(n-1,0,-1):
	     #mpos[i,:]=mpos[i-1,:]
	#mpos[0,:]=mpos[1,:]+vec
    #return mpos

#new
def bend_n(mpos123,m,rand,theta,phi):
	mpos=mpos123.copy()
	polar=[cos(theta)*sin(phi),sin(theta)*sin(phi),cos(phi)]
	if rand<.5:
		x=1
		end=len(mpos)-1
	else:
		x=-1
		end=0
	for i in range(m,end,x):
		BC=mpos123[i+x,:]-mpos123[i,:]
		z1=BC/dot(BC,BC)**.5
		k=(BC[0]**2+BC[1]**2)**.5 # + or - ?
		y1=[-BC[1]/k,BC[0]/k,0] #or y1=[BC[1]/k,-BC[0]/k,0]
		x1=[-z1[2]*y1[1],z1[2]*y1[0],z1[0]*y1[1]-z1[1]*y1[0]]
		transform=array([x1,y1,z1])
		untransform=transpose(transform)
		BC_polar=dot(transform,BC)
		rho=BC_polar[2]
		BC_polar=rho*polar
		BC_new=dot(untransform,BC_polar)
		mpos[i+x,:]=mpos[i,:]+BC_new
	return mpos 

#old
def bend(mpos123,m,rand,theta,phi): 
	#m=randint(1,len(mpos123)-2)
	mpos=mpos123.copy()
	#theta=2*pi*random()
	#phi=3./180*pi*random()
	x=[1,0,0]
	y=[0,1,0]
	z=[0,0,1]
	AB=[1,1,1]
	if rand<.5:
		for i in range(m,len(mpos)-1):
			BC=mpos123[i+1,:]-mpos123[i,:]
			z1=BC/dot(BC,BC)**.5
			y1=AB-dot(AB,BC)/dot(BC,BC)*BC
			y1=y1/dot(y1,y1)**.5
			x1=cross(z1,y1)
			untransform=[[dot(x,x1),dot(x,y1),dot(x,z1)],[dot(y,x1),dot(y,y1),dot(y,z1)],[dot(z,x1),dot(z,y1),dot(z,z1)]]
			transform=transpose(untransform)
			BC_polar=dot(transform,BC)
			rho=BC_polar[2]
			BC_polar[0]=rho*cos(theta)*sin(phi)
			BC_polar[1]=rho*sin(theta)*sin(phi)
			BC_polar[2]=rho*cos(phi)
			BC_new=dot(untransform,BC_polar)
			mpos[i+1,:]=mpos[i,:]+BC_new
	else:
		for i in range(m,0,-1):
			BC=mpos123[i-1,:]-mpos123[i,:]
			z1=BC/dot(BC,BC)**.5
			y1=AB-dot(AB,BC)/dot(BC,BC)*BC
			y1=y1/dot(y1,y1)**.5
			x1=cross(z1,y1)
			untransform=[[dot(x,x1),dot(x,y1),dot(x,z1)],[dot(y,x1),dot(y,y1),dot(y,z1)],[dot(z,x1),dot(z,y1),dot(z,z1)]]
			transform=transpose(untransform)
			BC_polar=dot(transform,BC)
			rho=BC_polar[2]
			BC_polar[0]=rho*cos(theta)*sin(phi)
			BC_polar[1]=rho*sin(theta)*sin(phi)
			BC_polar[2]=rho*cos(phi)
			BC_new=dot(untransform,BC_polar)
			mpos[i-1,:]=mpos[i,:]+BC_new
	return mpos
	    
def axistorsion(mpos123,m,rand,theta):
	mpos=mpos123.copy()
	rotate=array([[1,0,0],[0,cos(theta),sin(theta)],[0,-sin(theta),cos(theta)]])
	if rand <.5:
		x=1
	        end=len(mpos)-1
	else:
		x=-1
		end=0
	AB=mpos[m,:]-mpos[m-x,:]
	BC=mpos[m+x,:]-mpos[m,:]
	x1=AB/dot(AB,AB)**.5
        y1=BC-dot(BC,AB)/dot(AB,AB)*AB
        y1=y1/dot(y1,y1)**.5
        z1=cross(x1,y1)
        transform=array([x1,y1,z1])
	untransform=transpose(transform)
	for i in range(m,end,x):
		   bond=mpos123[i+x,:]-mpos123[i,:]
		   bond=dot(transform,bond.transpose())
		   bond=dot(rotate,bond)
		   bond=dot(untransform,bond)
		   mpos[i+x,:]=mpos[i,:]+bond
	return mpos 

#def enddist(mpos):
    #distvec=mpos[:][0]-mpos[:][7] #this is hardcoded
    #return dot(distvec,distvec)**.5

