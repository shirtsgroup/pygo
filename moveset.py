from numpy import *
from random import *
import scipy.optimize
import numpy
import energyfunc
import HMCforce
import pdb
import writetopdb

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
	    
def conrot(mpos123,m,rand,delw):
	mpos=mpos123.copy()
	solutions=[]
	jacobians=[]
	#get bond lengths and angles
	l=empty(6)
	u=empty((6,3))
	theta=empty(6)
	index=0
	for i in range(m,m+7):
		BA=mpos[i-1,:]-mpos[i,:]
		BC=mpos[i+1,:]-mpos[i,:]
		theta[index]=arccos(dot(BA,BC)/dot(BA,BA)**.5*dot(BC,BC)**.5)
		l[index]=dot(BA,BA)**.5
		u[index,:]=-BA/l[index]
		index += 1
	#get driver torsion w0
	AB=mpos[m-1,:]-mpos[m-2,:]
	BC=mpos[m,:]-mpos[m-1,:]
	CD=mpos[m+1,:]-mpos[m,:]
	plane1=cross(AB,BC)
        plane2=cross(BC,CD)
	w0=arccos((plane1[0]*plane2[0]+plane1[1]*plane2[1]+plane1[2]*plane2[2])/((plane1[0]**2+plane1[1]**2+plane1[2]**2)**.5*(plane2[0]**2+plane2[1]**2+plane2[2]**2)**.5))
	if ((plane1[0]*CD[0]+plane1[1]*CD[1]+plane1[2]*CD[2])<0):
		w0=-abs(w0)+2*pi
	else:
		w0=abs(w0)
	q1=array([[l[3]+l[2]*cos(theta[2])],[l[2]*sin(theta[2])],[0]])
	q2=array([[l[5]+l[4]*cos(theta[4])],[l[4]*sin(theta[4])],[0]])
	q11=array([[l[2]+l[3]*cos(theta[2])],[l[3]*sin(theta[2])],[0]])
	q22=array([[l[4]+l[5]*cos(theta[4])],[l[5]*sin(theta[4])],[0]])
	T1lab=zeros((3,3))
	T1lab[:,0]=u[1]
	T1lab[:,2]=cross(u[1],u[0])/sin(theta[1])
	T1lab[:,1]=cross(T1lab[:,2],T1lab[:,0])
	r51=linalg.solve(T1lab,mpos[m+5]-mpos[m])
	u61=linalg.solve(T1lab,u[6])
	def func1(w1):
		#branch 1
		T1=array([[cos(theta[1]),sin(theta[1]),0],[sin(theta[1])*cos(w1),-cos(theta[1])*cos(w1),sin(w1)],[sin(theta[1])*sin(w1),-cos(theta[1])*sin(w1),-cos(w1)]])
		t=linalg.solve(T1,r51-array([[l[1]],[0],[0]]))
		w=(sum(q1**2)-sum(q2**2)+sum(t**2)-2*t[0]*q11[0])/(2*q11[1])
		cosw2=(t[1]*w+t[2]*(sum(t**2)-t[0]**2-w**2)**.5)/(sum(t**2)-t[0]**2)
		sinw2=(t[2]*w-t[1]*(sum(t**2)-t[0]**2-w**2)**.5)/(sum(t**2)-t[0]**2)
		T2=array([[cos(theta[2]),sin(theta[2]),0],[sin(theta[2])*cosw2,-cos(theta[2])*cosw2,sinw2],[sin(theta[2])*sinw2,-cos(theta[2])*sinw2,-cosw2]])
		q3=linalg.solve(T2,t)
		q3=q3-q1
		cosw4=(q3[0]-cos(theta[3])*q22[0])/(sin(theta[3])*q22[1])
		sinw4=sin(arccos(cosw4))
		A=q22[2]*sinw4
		B=sin(theta[3])*q22[0]-cos(theta[3])*q22[1]*cosw4
		sinw3=(q33[1]+B/A*q33[2])/(B**2/A+A)
		cosw3=(B*sinw3-q33[2])/A
		T3=array([[cos(theta[3]),sin(theta[3]),0],[sin(theta[3])*cosw3,-cos(theta[3])*cosw3,sinw3],[sin(theta[3])*sinw3,-cos(theta[3])*sinw3,-cosw3]])
		T4=array([[cos(theta[4]),sin(theta[4]),0],[sin(theta[4])*cosw4,-cos(theta[4])*cosw4,sinw4],[sin(theta[4])*sinw4,-cos(theta[4])*sinw4,-cosw4]])
		ans=dot(T4,array([[1],[0],[0]]))
		ans=dot(T3,ans)
		ans=dot(T2,ans)
		ans=dot(T1,ans)
		ans=dot(u61.transpose(),ans)
		ans=ans-cos(theta[5])
		return ans
	def func2(w1):
		#branch 2
		T1=array([[cos(theta[1]),sin(theta[1]),0],[sin(theta[1])*cos(w1),-cos(theta[1])*cos(w1),sin(w1)],[sin(theta[1])*sin(w1),-cos(theta[1])*sin(w1),-cos(w1)]])
		t=linalg.solve(T1,r51-array([[l[1]],[0],[0]]))
		w=(sum(q1**2)-sum(q2**2)+sum(t**2)-2*t[0]*q11[0])/(2*q11[1])
		cosw2=(t[1]*w+t[2]*(sum(t**2)-t[0]**2-w**2)**.5)/(sum(t**2)-t[0]**2)
		sinw2=(t[2]*w-t[1]*(sum(t**2)-t[0]**2-w**2)**.5)/(sum(t**2)-t[0]**2)
		T2=array([[cos(theta[2]),sin(theta[2]),0],[sin(theta[2])*cosw2,-cos(theta[2])*cosw2,sinw2],[sin(theta[2])*sinw2,-cos(theta[2])*sinw2,-cosw2]])
		q3=linalg.solve(T2,t)
		q3=q3-q1
		cosw4=(q3[0]-cos(theta[3])*q22[0])/(sin(theta[3])*q22[1])
		sinw4=sin(-arccos(cosw4))
		A=q22[2]*sinw4
		B=sin(theta[3])*q22[0]-cos(theta[3])*q22[1]*cosw4
		sinw3=(q33[1]+B/A*q33[2])/(B**2/A+A)
		cosw3=(B*sinw3-q33[2])/A
		T3=array([[cos(theta[3]),sin(theta[3]),0],[sin(theta[3])*cosw3,-cos(theta[3])*cosw3,sinw3],[sin(theta[3])*sinw3,-cos(theta[3])*sinw3,-cosw3]])
		T4=array([[cos(theta[4]),sin(theta[4]),0],[sin(theta[4])*cosw4,-cos(theta[4])*cosw4,sinw4],[sin(theta[4])*sinw4,-cos(theta[4])*sinw4,-cosw4]])
		ans=dot(T4,array([[1],[0],[0]]))
		ans=dot(T3,ans)
		ans=dot(T2,ans)
		ans=dot(T1,ans)
		ans=dot(u61.transpose(),ans)
		ans=ans-cos(theta[5])
		return ans
	def func3(w1):
		#branch 1
		T1=array([[cos(theta[1]),sin(theta[1]),0],[sin(theta[1])*cos(w1),-cos(theta[1])*cos(w1),sin(w1)],[sin(theta[1])*sin(w1),-cos(theta[1])*sin(w1),-cos(w1)]])
		t=linalg.solve(T1,r51-array([[l[1]],[0],[0]]))
		w=(sum(q1**2)-sum(q2**2)+sum(t**2)-2*t[0]*q11[0])/(2*q11[1])
		cosw2=(t[1]*w-t[2]*(sum(t**2)-t[0]**2-w**2)**.5)/(sum(t**2)-t[0]**2)
		sinw2=(t[2]*w+t[1]*(sum(t**2)-t[0]**2-w**2)**.5)/(sum(t**2)-t[0]**2)
		T2=array([[cos(theta[2]),sin(theta[2]),0],[sin(theta[2])*cosw2,-cos(theta[2])*cosw2,sinw2],[sin(theta[2])*sinw2,-cos(theta[2])*sinw2,-cosw2]])
		q3=linalg.solve(T2,t)
		q3=q3-q1
		cosw4=(q3[0]-cos(theta[3])*q22[0])/(sin(theta[3])*q22[1])
		sinw4=sin(arccos(cosw4))
		A=q22[2]*sinw4
		B=sin(theta[3])*q22[0]-cos(theta[3])*q22[1]*cosw4
		sinw3=(q33[1]+B/A*q33[2])/(B**2/A+A)
		cosw3=(B*sinw3-q33[2])/A
		T3=array([[cos(theta[3]),sin(theta[3]),0],[sin(theta[3])*cosw3,-cos(theta[3])*cosw3,sinw3],[sin(theta[3])*sinw3,-cos(theta[3])*sinw3,-cosw3]])
		T4=array([[cos(theta[4]),sin(theta[4]),0],[sin(theta[4])*cosw4,-cos(theta[4])*cosw4,sinw4],[sin(theta[4])*sinw4,-cos(theta[4])*sinw4,-cosw4]])
		ans=dot(T4,array([[1],[0],[0]]))
		ans=dot(T3,ans)
		ans=dot(T2,ans)
		ans=dot(T1,ans)
		ans=dot(u61.transpose(),ans)
		ans=ans-cos(theta[5])
		return ans
	def func4(w1):
		#branch 1
		T1=array([[cos(theta[1]),sin(theta[1]),0],[sin(theta[1])*cos(w1),-cos(theta[1])*cos(w1),sin(w1)],[sin(theta[1])*sin(w1),-cos(theta[1])*sin(w1),-cos(w1)]])
		t=linalg.solve(T1,r51-array([[l[1]],[0],[0]]))
		w=(sum(q1**2)-sum(q2**2)+sum(t**2)-2*t[0]*q11[0])/(2*q11[1])
		cosw2=(t[1]*w-t[2]*(sum(t**2)-t[0]**2-w**2)**.5)/(sum(t**2)-t[0]**2)
		sinw2=(t[2]*w+t[1]*(sum(t**2)-t[0]**2-w**2)**.5)/(sum(t**2)-t[0]**2)
		T2=array([[cos(theta[2]),sin(theta[2]),0],[sin(theta[2])*cosw2,-cos(theta[2])*cosw2,sinw2],[sin(theta[2])*sinw2,-cos(theta[2])*sinw2,-cosw2]])
		q3=linalg.solve(T2,t)
		q3=q3-q1
		cosw4=(q3[0]-cos(theta[3])*q22[0])/(sin(theta[3])*q22[1])
		sinw4=sin(-arccos(cosw4))
		A=q22[2]*sinw4
		B=sin(theta[3])*q22[0]-cos(theta[3])*q22[1]*cosw4
		sinw3=(q33[1]+B/A*q33[2])/(B**2/A+A)
		cosw3=(B*sinw3-q33[2])/A
		T3=array([[cos(theta[3]),sin(theta[3]),0],[sin(theta[3])*cosw3,-cos(theta[3])*cosw3,sinw3],[sin(theta[3])*sinw3,-cos(theta[3])*sinw3,-cosw3]])
		T4=array([[cos(theta[4]),sin(theta[4]),0],[sin(theta[4])*cosw4,-cos(theta[4])*cosw4,sinw4],[sin(theta[4])*sinw4,-cos(theta[4])*sinw4,-cosw4]])
		ans=dot(T4,array([[1],[0],[0]]))
		ans=dot(T3,ans)
		ans=dot(T2,ans)
		ans=dot(T1,ans)
		ans=dot(u61.transpose(),ans)
		ans=ans-cos(theta[5])
		return ans
	# prerotation solutions
	pass
	# postrotation solutions
	

def localmove(mpos123,m,rand,theta):
	mpos=mpos123.copy()
	c=cos(theta)
	s=sin(theta)
	rotate=array([[1,0,0],[0,c,s],[0,-s,c]])
	if rand <.5:
		n=1
	        #end=len(mpos)-1
	else:
		n=-1
		#end=0
	AB=mpos[m,:]-mpos[m-n,:]
	BC=mpos[m+n,:]-mpos[m,:]
	
	dotAB=AB[0]*AB[0]+AB[1]*AB[1]+AB[2]*AB[2]
	x=AB/dotAB**.5
        y=BC-(BC[0]*AB[0]+BC[1]*AB[1]+BC[2]*AB[2])/dotAB*AB
        y=y/(y[0]*y[0]+y[1]*y[1]+y[2]*y[2])**.5
        z=[x[1]*y[2]-x[2]*y[1],x[2]*y[0]-x[0]*y[2],x[0]*y[1]-x[1]*y[0]]
        transform=array([x,y,z])
	untransform=transpose(transform)
	for i in range(m,m+n*3,n):
		bond=mpos123[i+n,:]-mpos123[i,:]
		bond=dot(transform,bond.transpose())
		bond=dot(rotate,bond)
		bond=dot(untransform,bond)
		mpos[i+n,:]=mpos[i,:]+bond
	AB=mpos[m+n*5,:]-mpos[m+n*3,:]
	BC=mpos[m+n*4,:]-mpos[m+n*3,:]
	dotAB=AB[0]*AB[0]+AB[1]*AB[1]+AB[2]*AB[2]
	x=AB/dotAB**.5
        y=BC-(BC[0]*AB[0]+BC[1]*AB[1]+BC[2]*AB[2])/dotAB*AB
        y=y/(y[0]*y[0]+y[1]*y[1]+y[2]*y[2])**.5
        z=[x[1]*y[2]-x[2]*y[1],x[2]*y[0]-x[0]*y[2],x[0]*y[1]-x[1]*y[0]]
        transform=array([x,y,z])
	untransform=transpose(transform)
	l42=dot(mpos123[m+n*4,:]-mpos123[m+n*3,:],mpos123[m+n*4,:]-mpos123[m+n*3,:])
	l4=l42**.5
	l52=dot(mpos123[m+n*5,:]-mpos123[m+n*4,:],mpos123[m+n*5,:]-mpos123[m+n*4,:])	
	q2=dot(mpos[m+n*5,:]-mpos[m+n*3,:],mpos[m+n*5,:]-mpos[m+n*3,:])
	cosine=(l42+q2-l52)/(2*l4*q2**.5)
	if abs(cosine)>1:
		return float('nan')
	else:
		sine=sin(arccos(cosine))
		bond=array([l4*cosine,l4*sine,0])
		mpos[m+n*4,:]=dot(untransform,bond)+mpos[m+n*3,:]
		return mpos 
	
def runMD(self,nsteps,h,dict):
	numbeads=dict['numbeads']
	numint=dict['numint']
	angleparam=dict['angleparam']
	torsparam=dict['torsparam']
	nativeparam_n=dict['nativeparam_n']
	nonnativeparam=dict['nonnativeparam']
	nnepsil=dict['nnepsil']
	m=dict['mass']
	positiontemplate=dict['positiontemplate']
	
	#m=120.368 # average mass of all residues
	tol=1e-8
	maxloop=100
	
	self.newcoord=self.coord.copy()
	self.vel=empty((numbeads,3))
	for i in range(numbeads): self.vel[i,:]=numpy.random.normal(0,(4.184*self.kb*self.T/m[i])**.5,3) #in nm/ps, uses average residue mass
	bonds=self.coord[0:numbeads-1,:]-self.coord[1:numbeads,:]
	d2=numpy.sum(bonds**2,axis=1)
	d=d2**.5
	force=HMCforce.bondedforces(self.coord,torsparam,angleparam,bonds,d2,d)+HMCforce.nonbondedforces(self.coord,numint,numbeads,nativeparam_n,nonnativeparam,nnepsil)
	self.oldH=self.u0+.5/4.184*numpy.sum(m*numpy.sum(self.vel**2,axis=1)) # in kcal/mol
	
	for e in range(nsteps):
		#finding r(t+dt)
		loops=0
		conv=numpy.ones(numbeads-1)
		a=numpy.transpose(force)/m
		q=self.vel+h/2*transpose(a)
		while numpy.sum(conv)!=0 and loops<maxloop:
			for i in range(numbeads-1): # loops through all bonds/constraints
				s=bonds[i]+h*(q[i,:]-q[i+1,:])
				diff=s[0]*s[0]+s[1]*s[1]+s[2]*s[2]-d2[i]
				#print diff
				if (numpy.abs(diff)<tol):
					conv[i]=0
				else:
					g=diff/(2*h*(s[0]*bonds[i,0]+s[1]*bonds[i,1]+s[2]*bonds[i,2])*(1/m[i]+1/m[i+1]))
					q[i,:]-=g/m[i]*bonds[i]
					q[i+1,:]+=g/m[i+1]*bonds[i]
			loops += 1
		if (loops==maxloop):
			abc=q
			pdb.set_trace()
			break
		self.newcoord+=h*q #r(t+dt)
		self.vel=q # new v(t)
		bonds=self.newcoord[0:numbeads-1,:]-self.newcoord[1:numbeads,:] #rij(t+dt)
		force=HMCforce.bondedforces(self.newcoord,torsparam,angleparam,bonds,d2,d)+HMCforce.nonbondedforces(self.newcoord,numint,numbeads,nativeparam_n,nonnativeparam,nnepsil)
		
		#finding v(t+dt)
		loops=0
		conv=numpy.ones(numbeads-1)
		a=transpose(force)/m
		self.vel+=h/2*transpose(a)
		while numpy.sum(conv)!=0 and loops<maxloop:
			for i in range(numbeads-1):
				vij=self.vel[i,:]-self.vel[i+1,:]
				diff=bonds[i,0]*vij[0]+bonds[i,1]*vij[1]+bonds[i,2]*vij[2]
				if (numpy.abs(diff)<tol):
					conv[i]=0
				else:
					k=diff/(d2[i]*(1/m[i]+1/m[i+1]))
					self.vel[i] -= k/m[i]*bonds[i,:]
					self.vel[i+1] += k/m[i+1]*bonds[i,:]
			loops += 1
		if (loops==maxloop):
			break
		writetopdb.addtopdb(self.newcoord,positiontemplate,self.move*nsteps+e,'%s/trajectory%i.pdb' % (self.out,int(self.T)))
	if(loops==maxloop):
		self.uncloseable=True
		self.rejected += 1
	return self

def enddist(mpos):
    distvec=mpos[0,:]-mpos[-1,:]
    return dot(distvec,distvec)**.5

