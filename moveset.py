from numpy import *
from random import *

def crankshaft(mpos123,m,theta): 
    mpos=mpos123.copy()

    AB=mpos[m,:]-mpos[m-1,:]
    AC=mpos[m+1,:]-mpos[m-1,:]
    
    dotAC=AC[0]*AC[0]+AC[1]*AC[1]+AC[2]*AC[2]
    
    x=AC/dotAC**.5
    y=AB-(AB[0]*AC[0]+AB[1]*AC[1]+AB[2]*AC[2])/dotAC*AC
    y=y/(y[0]*y[0]+y[1]*y[1]+y[2]*y[2])**.5
    z=[x[1]*y[2]-x[2]*y[1],x[2]*y[0]-x[0]*y[2],x[0]*y[1]-x[1]*y[0]]
    #transform=array([x,y,z])
    #untransform=transpose(transform)
    
    c=cos(theta)
    s=sin(theta)
    #rotate=array([[1,0,0],[0,c,s],[0,-s,c]])

    AB1=[x[0]*AB[0]+x[1]*AB[1]+x[2]*AB[2],y[0]*AB[0]+y[1]*AB[1]+y[2]*AB[2],z[0]*AB[0]+z[1]*AB[1]+z[2]*AB[2]] #dot(transform,AB.transpose())
    AB2=[AB1[0],c*AB1[1]+s*AB1[2],-s*AB1[1]+c*AB1[2]] #dot(rotate,AB1)
    AB=[x[0]*AB2[0]+y[0]*AB2[1]+z[0]*AB2[2],x[1]*AB2[0]+y[1]*AB2[1]+z[1]*AB2[2],x[2]*AB2[0]+y[2]*AB2[1]+z[2]*AB2[2]] #dot(untransform,AB2)
    mpos[m,:]=mpos[m-1,:]+AB
    return mpos

def axistorsion(mpos123,m,rand,theta):
	mpos=mpos123.copy()
	c=cos(theta)
	s=sin(theta)
	rotate=array([[1,0,0],[0,c,s],[0,-s,c]])
	if rand <.5:
		n=1
	        end=len(mpos)-1
	else:
		n=-1
		end=0
	AB=mpos[m,:]-mpos[m-n,:]
	BC=mpos[m+n,:]-mpos[m,:]
	
	dotAB=AB[0]*AB[0]+AB[1]*AB[1]+AB[2]*AB[2]
	x=AB/dotAB**.5
        y=BC-(BC[0]*AB[0]+BC[1]*AB[1]+BC[2]*AB[2])/dotAB*AB
        y=y/(y[0]*y[0]+y[1]*y[1]+y[2]*y[2])**.5
        z=[x[1]*y[2]-x[2]*y[1],x[2]*y[0]-x[0]*y[2],x[0]*y[1]-x[1]*y[0]]
        transform=array([x,y,z])
	untransform=transpose(transform)
	for i in range(m,end,n):
		   bond=mpos123[i+n,:]-mpos123[i,:]
		   bond=dot(transform,bond.transpose())
		   bond=dot(rotate,bond)
		   bond=dot(untransform,bond)
		   mpos[i+n,:]=mpos[i,:]+bond
	return mpos 

#def axistorsion_n(mpos123,m,rand,theta):
	#mpos=mpos123.copy()
	#c=cos(theta)
	#s=sin(theta)
	#rotate=array([[1,0,0],[0,c,s],[0,-s,c]])
	#if rand <.5:
		#n=1
	        #end=len(mpos)-1
	#else:
		#n=-1
		#end=0
	#AB=mpos[m,:]-mpos[m-n,:]
	#BC=mpos[m+n,:]-mpos[m,:]
	
	#dotAB=AB[0]*AB[0]+AB[1]*AB[1]+AB[2]*AB[2]
	#x=AB/dotAB**.5
        #y=BC-(BC[0]*AB[0]+BC[1]*AB[1]+BC[2]*AB[2])/dotAB*AB
        #y=y/(y[0]*y[0]+y[1]*y[1]+y[2]*y[2])**.5
        #z=[x[1]*y[2]-x[2]*y[1],x[2]*y[0]-x[0]*y[2],x[0]*y[1]-x[1]*y[0]]
        #transform=array([x,y,z])
	#untransform=transpose(transform)
	#for i in range(m,end,n):
		   #bond=mpos123[i+n,:]-mpos123[i,:]
    		   #bond1=[x[0]*bond[0]+x[1]*bond[1]+x[2]*bond[2],y[0]*bond[0]+y[1]*bond[1]+y[2]*bond[2],z[0]*bond[0]+z[1]*bond[1]+z[2]*bond[2]] #dot(transform,AB.transpose())
   		   #bond2=[bond1[0],c*bond1[1]+s*bond1[2],-s*bond1[1]+c*bond1[2]] #dot(rotate,AB1)
  		   #bond=[x[0]*bond2[0]+y[0]*bond2[1]+z[0]*bond2[2],x[1]*bond2[0]+y[1]*bond2[1]+z[1]*bond2[2],x[2]*bond2[0]+y[2]*bond2[1]+z[2]*bond2[2]] #dot(untransform,AB2)
		   #mpos[i+n,:]=mpos[i,:]+bond
	#return mpos 

def anglebend(mpos123,m,rand,theta):
	mpos=mpos123.copy()
	c=cos(theta)
	s=sin(theta)
	rotate=array([[1,0,0],[0,c,s],[0,-s,c]])
	if rand <.5:
		n=1
	        end=len(mpos)-1
	else:
		n=-1
		end=0	
	AB=mpos[m,:]-mpos[m-n,:]
	BC=mpos[m+n,:]-mpos[m,:]
	dotAB=AB[0]*AB[0]+AB[1]*AB[1]+AB[2]*AB[2]
	z=AB/dotAB**.5
        y=BC-(BC[0]*AB[0]+BC[1]*AB[1]+BC[2]*AB[2])/dotAB*AB
        y=y/(y[0]*y[0]+y[1]*y[1]+y[2]*y[2])**.5
        x=[z[1]*y[2]-z[2]*y[1],z[2]*y[0]-z[0]*y[2],z[0]*y[1]-z[1]*y[0]]
        transform=array([x,y,z])
	untransform=transpose(transform)
	for i in range(m,end,n):
		   bond=mpos123[i+n,:]-mpos123[i,:]
		   bond=dot(transform,bond.transpose())
		   bond=dot(rotate,bond)
		   bond=dot(untransform,bond)
		   mpos[i+n,:]=mpos[i,:]+bond
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

#new
#def bend_n(mpos123,m,rand,theta,phi):
	#mpos=mpos123.copy()
	#polar=array([cos(theta)*sin(phi),sin(theta)*sin(phi),cos(phi)])
	#if rand<.5:
		#x=1
		#end=len(mpos)-1
	#else:
		#x=-1
		#end=0
	#for i in range(m,end,x):
		#BC=mpos123[i+x,:]-mpos123[i,:]
		#z1=BC/dot(BC,BC)**.5
		#k=(BC[0]**2+BC[1]**2)**.5 # + or - ?
		#y1=[-BC[1]/k,BC[0]/k,0] #or y1=[BC[1]/k,-BC[0]/k,0]
		#x1=[-z1[2]*y1[1],z1[2]*y1[0],z1[0]*y1[1]-z1[1]*y1[0]]
		#transform=array([x1,y1,z1])
		#untransform=transpose(transform)
		#BC_polar=dot(transform,BC)
		#rho=BC_polar[2]
		#BC_polar=rho*polar
		#BC_new=dot(untransform,BC_polar.transpose())
		#mpos[i+x,:]=mpos[i,:]+BC_new
	#return mpos 

##old
#def bend(mpos123,m,rand,theta,phi): 
	##m=randint(1,len(mpos123)-2)
	#mpos=mpos123.copy()
	##theta=2*pi*random()
	##phi=3./180*pi*random()
	#x=[1,0,0]
	#y=[0,1,0]
	#z=[0,0,1]
	#AB=[1,1,1]
	#if rand<.5:
		#for i in range(m,len(mpos)-1):
			#BC=mpos123[i+1,:]-mpos123[i,:]
			#z1=BC/dot(BC,BC)**.5
			#y1=AB-dot(AB,BC)/dot(BC,BC)*BC
			#y1=y1/dot(y1,y1)**.5
			#x1=cross(z1,y1)
			#untransform=[[dot(x,x1),dot(x,y1),dot(x,z1)],[dot(y,x1),dot(y,y1),dot(y,z1)],[dot(z,x1),dot(z,y1),dot(z,z1)]]
			#transform=transpose(untransform)
			#BC_polar=dot(transform,BC)
			#rho=BC_polar[2]
			#BC_polar[0]=rho*cos(theta)*sin(phi)
			#BC_polar[1]=rho*sin(theta)*sin(phi)
			#BC_polar[2]=rho*cos(phi)
			#BC_new=dot(untransform,BC_polar)
			#mpos[i+1,:]=mpos[i,:]+BC_new
	#else:
		#for i in range(m,0,-1):
			#BC=mpos123[i-1,:]-mpos123[i,:]
			#z1=BC/dot(BC,BC)**.5
			#y1=AB-dot(AB,BC)/dot(BC,BC)*BC
			#y1=y1/dot(y1,y1)**.5
			#x1=cross(z1,y1)
			#untransform=[[dot(x,x1),dot(x,y1),dot(x,z1)],[dot(y,x1),dot(y,y1),dot(y,z1)],[dot(z,x1),dot(z,y1),dot(z,z1)]]
			#transform=transpose(untransform)
			#BC_polar=dot(transform,BC)
			#rho=BC_polar[2]
			#BC_polar[0]=rho*cos(theta)*sin(phi)
			#BC_polar[1]=rho*sin(theta)*sin(phi)
			#BC_polar[2]=rho*cos(phi)
			#BC_new=dot(untransform,BC_polar)
			#mpos[i-1,:]=mpos[i,:]+BC_new
	#return mpos
	    

def enddist(mpos):
    distvec=mpos[:][0]-mpos[:][7] #this is hardcoded
    return dot(distvec,distvec)**.5

