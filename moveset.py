from numpy import *
from random import *

def crankshaft(mpos123): #moleculeposition, will need bonds later
    mpos=mpos123.copy()
    m=randint(1,len(mpos)-2) #random molecule, not end ones

    posb=mpos[m,:] # middle molecule, will get crankshafted
    posa=mpos[m-1,:] #one 'before' it
    posc=mpos[m+1,:] # one 'after' it

    AB=posb-posa
    AC=posc-posa
    
    x=[1,0,0]
    y=[0,1,0]
    z=[0,0,1]
    x1=AC/dot(AC,AC)**.5
    y1=AB-dot(AB,AC)/dot(AC,AC)*AC
    y1=y1/dot(y1,y1)**.5
    z1=cross(x1,y1)
    untransform=[[dot(x,x1),dot(x,y1),dot(x,z1)],[dot(y,x1),dot(y,y1),dot(y,z1)],[dot(z,x1),dot(z,y1),dot(z,z1)]]
    transform=transpose(untransform)
    
    
    theta=2*pi*random()
    rotate=array([[1,0,0],[0,cos(theta),sin(theta)],[0,-sin(theta),cos(theta)]])

    AB1=dot(transform,AB.transpose())
    AB2=dot(rotate,AB1)
    AB=dot(untransform,AB2)
    mpos[m,:]=posa+AB.transpose()
    return mpos

def reptation(mpos123):
    mpos=mpos123.copy()
    theta=arccos(1-2*random())
    phi=2*pi*random()
    rho=3.86470703286 #unhardcode this
    vec=[rho*sin(theta)*cos(phi),rho*sin(theta)*sin(phi),rho*cos(theta)]    
    n=len(mpos)
    if random() < .5:    
	for i in range(n-1):
             mpos[i,:]=mpos[i+1,:]
    	mpos[n-1,:]=mpos[n-2][:]+vec
    else:
	for i in range(n-1,0,-1):
	     mpos[i,:]=mpos[i-1,:]
	mpos[0,:]=mpos[1,:]+vec
    return mpos

def torsion(mpos123):
	mpos=mpos123.copy()
	dtheta=arccos(1-2*random())
	dphi=2*pi*random()
	m=randint(1,len(mpos)-2)
	if random() <.5:
	    for i in range(m,len(mpos)-1):
		BC=mpos123[i+1,:]-mpos123[i,:]
        	r=dot(BC,BC)**.5
		theta=arccos(BC[2]/r)+dtheta #new theta
		phi=arctan(BC[1]/BC[0])+dphi #new phi
		BC[0]=r*sin(theta)*cos(phi)
		BC[1]=r*sin(theta)*sin(phi)
		BC[2]=r*cos(theta)
		mpos[i+1,:]=mpos[i,:]+BC
	else:
	    for i in range(m,0,-1):
                BC=mpos123[i-1,:]-mpos123[i,:]
                r=dot(BC,BC)**.5
                theta=arccos(BC[2]/r)+dtheta #new theta
                phi=arctan(BC[1]/BC[0])+dphi #new phi
                BC[0]=r*sin(theta)*cos(phi)
                BC[1]=r*sin(theta)*sin(phi)
                BC[2]=r*cos(theta)
                mpos[i-1,:]=mpos[i,:]+BC
	return mpos

def axistorsion(mpos123):
	mpos=mpos123.copy()
	theta=2*pi*random() #not spherical coord theta, just arbitrary angle
	rotate=array([[1,0,0],[0,cos(theta),sin(theta)],[0,-sin(theta),cos(theta)]])
	m=randint(1,len(mpos)-2)
	posb=mpos[m,:] # position of randomly chosen bead
   	posa=mpos[m-1,:] #one 'before' it
    	posc=mpos[m+1,:] # one 'after' it
	AB=posb-posa
        BC=posc-posb
	x=[1,0,0]
        y=[0,1,0]
        z=[0,0,1]
	if random() <.5:
	    x1=AB/dot(AB,AB)**.5
            y1=BC-dot(BC,AB)/dot(AB,AB)*AB
            y1=y1/dot(y1,y1)**.5
            z1=cross(x1,y1)
            untransform=[[dot(x,x1),dot(x,y1),dot(x,z1)],[dot(y,x1),dot(y,y1),dot(y,z1)],[dot(z,x1),dot(z,y1),dot(z,z1)]]
            transform=transpose(untransform)
	    for i in range(m,len(mpos)-1):
		   bond=mpos123[i+1,:]-mpos123[i,:]
		   bond=dot(transform,bond.transpose())
		   bond=dot(rotate,bond)
		   bond=dot(untransform,bond)
		   mpos[i+1,:]=mpos[i,:]+bond
	else:
	    x1=BC/dot(BC,BC)**.5
            y1=AB-dot(AB,BC)/dot(BC,BC)*BC
            y1=y1/dot(y1,y1)**.5
            z1=cross(x1,y1)
            untransform=[[dot(x,x1),dot(x,y1),dot(x,z1)],[dot(y,x1),dot(y,y1),dot(y,z1)],[dot(z,x1),dot(z,y1),dot(z,z1)]]
            transform=transpose(untransform)
	    for i in range(m,0,-1):
		   bond=mpos123[i-1,:]-mpos123[i,:]
		   bond=dot(transform,bond.transpose())
		   bond=dot(rotate,bond)
		   bond=dot(untransform,bond)
		   mpos[i-1,:]=mpos[i,:]+bond
	return mpos 

def enddist(mpos):
    distvec=mpos[:][0]-mpos[:][7] #this is hardcoded
    return dot(distvec,distvec)**.5

