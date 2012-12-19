from numpy import *
from random import *
import numpy
import energyfunc
import HMCforce
import pdb
import writetopdb
from scipy import weave
#seed(10)
#numpy.random.seed(10)

def crankshaft(mpos123, m, theta): 
    mpos = mpos123.copy()

    AB = mpos[m,:] - mpos[m-1,:]
    AC = mpos[m+1,:] - mpos[m-1,:]
    
    dotAC = AC[0]*AC[0] + AC[1]*AC[1] + AC[2]*AC[2]
    
    x = AC / dotAC**.5
    y = AB - (AB[0]*AC[0] + AB[1]*AC[1] + AB[2]*AC[2]) / dotAC * AC
    y = y / (y[0]*y[0] + y[1]*y[1] + y[2]*y[2])**.5
    z = [x[1]*y[2]-x[2]*y[1], x[2]*y[0]-x[0]*y[2], x[0]*y[1]-x[1]*y[0]]
    #transform=array([x,y,z])
    #untransform=transpose(transform)
    
    c = cos(theta)
    s = sin(theta)
    #rotate=array([[1,0,0],[0,c,s],[0,-s,c]])

    AB1 = [x[0]*AB[0]+x[1]*AB[1]+x[2]*AB[2],y[0]*AB[0]+y[1]*AB[1]+y[2]*AB[2],z[0]*AB[0]+z[1]*AB[1]+z[2]*AB[2]] #dot(transform,AB.transpose())
    AB2 = [AB1[0],c*AB1[1]+s*AB1[2],-s*AB1[1]+c*AB1[2]] #dot(rotate,AB1)
    AB = [x[0]*AB2[0]+y[0]*AB2[1]+z[0]*AB2[2],x[1]*AB2[0]+y[1]*AB2[1]+z[1]*AB2[2],x[2]*AB2[0]+y[2]*AB2[1]+z[2]*AB2[2]] #dot(untransform,AB2)
    mpos[m,:] = mpos[m-1,:] + AB
    return mpos

def globalcrank(mpos123,theta):
    mpos = mpos123.copy()
    for m in range(1, len(mpos)-1):
        AB = mpos[m,:] - mpos[m-1,:]
        AC = mpos[m+1,:] - mpos[m-1,:]
        dotAC = AC[0]*AC[0] + AC[1]*AC[1] + AC[2]*AC[2]
        x = AC / dotAC**.5
        y = AB - (AB[0]*AC[0] + AB[1]*AC[1] + AB[2]*AC[2]) / dotAC * AC
        y = y / (y[0]*y[0] + y[1]*y[1] + y[2]*y[2])**.5
        z = [x[1]*y[2]-x[2]*y[1], x[2]*y[0]-x[0]*y[2], x[0]*y[1]-x[1]*y[0]]
        c = cos(theta[m-1])
        s = sin(theta[m-1])
        AB1 = [x[0]*AB[0]+x[1]*AB[1]+x[2]*AB[2],y[0]*AB[0]+y[1]*AB[1]+y[2]*AB[2],z[0]*AB[0]+z[1]*AB[1]+z[2]*AB[2]] #dot(transform,AB.transpose())
        AB2 = [AB1[0],c*AB1[1]+s*AB1[2],-s*AB1[1]+c*AB1[2]] #dot(rotate,AB1)
        AB = [x[0]*AB2[0]+y[0]*AB2[1]+z[0]*AB2[2],x[1]*AB2[0]+y[1]*AB2[1]+z[1]*AB2[2],x[2]*AB2[0]+y[2]*AB2[1]+z[2]*AB2[2]] #dot(untransform,AB2)
        mpos[m,:] = mpos[m-1,:] + AB
    return mpos
    
def cglobalcrank(mpos123,theta):
    mpos = mpos123.copy()
    numbeads = len(mpos)
    code = """
    double c, s, dotAC, dotABBC, mag;
    double xAB, yAB, zAB, xAC, yAC, zAC;
    double xx, xy, xz, yx, yy, yz, zx, zy, zz;
    double btx, bty, btz;
    double brx, bry, brz, bxx, byy, bzz;
    for ( int m = 1; m < numbeads-1; m++){
        xAB = MPOS2(m,0) - MPOS2(m-1,0); yAB = MPOS2(m,1) - MPOS2(m-1,1); zAB = MPOS2(m,2) - MPOS2(m-1,2);
        xAC = MPOS2(m+1,0) - MPOS2(m-1,0); yAC = MPOS2(m+1,1) - MPOS2(m-1,1); zAC = MPOS2(m+1,2) - MPOS2(m-1,2);
        dotAC = xAC*xAC + yAC*yAC + zAC*zAC;
        dotABBC = xAB*xAC + yAB*yAC + zAB*zAC;
        mag = sqrt(dotAC);
        xx = xAC/mag; xy = yAC/mag; xz = zAC/mag;
        mag = dotABBC/dotAC;
        yx = xAB - mag*xAC; yy = yAB - mag*yAC; yz = zAB - mag*zAC;
        mag = sqrt(yx*yx+yy*yy+yz*yz);
        yx = yx/mag; yy = yy/mag; yz = yz/mag;
        zx = xy*yz - xz*yy; zy = xz*yx - xx*yz; zz = xx*yy - xy*yx;
	c = cos(THETA1(m-1));
	s = sin(THETA1(m-1));
	btx = xx*xAB + xy*yAB + xz*zAB; bty = yx*xAB + yy*yAB + yz*zAB; btz = zx*xAB + zy*yAB + zz*zAB;
	brx = btx; bry = c*bty + s*btz; brz = -s*bty + c*btz;
	bxx = xx*brx + yx*bry + zx*brz; byy = xy*brx + yy*bry + zy*brz; bzz = xz*brx + yz*bry + zz*brz;
	MPOS2(m,0) = MPOS2(m-1,0) + bxx; MPOS2(m,1) = MPOS2(m-1,1) + byy; MPOS2(m,2) = MPOS2(m-1,2) + bzz;
    }
    
    """
    info = weave.inline(code, ['mpos','theta','numbeads'], headers = ['<math.h>', '<stdlib.h>'])
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

def caxistorsion(mposs,m,rand,theta):
    mpos = mposs.copy()
    numbeads = len(mpos)
    code = """
    int n, end;
    double c, s, dotAB, dotABBC, mag, angle;
    double xAB, yAB, zAB, xBC, yBC, zBC;
    double xx, xy, xz, yx, yy, yz, zx, zy, zz;
    double bx, by, bz, btx, bty, btz;
    double brx, bry, brz, bxx, byy, bzz;
    c = cos(theta);
    s = sin(theta);
    if (rand < .5){
        n = 1;
        end = numbeads - 1;
    }
    else{
        n = -1;
        end = 0;
    }
    xAB = MPOS2(m,0) - MPOS2(m-n,0); yAB = MPOS2(m,1) - MPOS2(m-n,1); zAB = MPOS2(m,2) - MPOS2(m-n,2);
    xBC = MPOS2(m+n,0) - MPOS2(m,0); yBC = MPOS2(m+n,1) - MPOS2(m,1); zBC = MPOS2(m+n,2) - MPOS2(m,2);
    dotAB = xAB*xAB + yAB*yAB + zAB*zAB;
    dotABBC = xAB*xBC + yAB*yBC + zAB*zBC;
    mag = sqrt(dotAB);
    xx = xAB/mag; xy = yAB/mag; xz = zAB/mag;
    mag = dotABBC/dotAB;
    yx = xBC - mag*xAB; yy = yBC - mag*yAB; yz = zBC - mag*zAB;
    mag = sqrt(yx*yx+yy*yy+yz*yz);
    yx = yx/mag; yy = yy/mag; yz = yz/mag;
    zx = xy*yz - xz*yy; zy = xz*yx - xx*yz; zz = xx*yy - xy*yx;
    for( int i = m; i != end; i += n){
            bx = MPOSS2(i+n,0) - MPOSS2(i,0); by = MPOSS2(i+n,1) - MPOSS2(i,1); bz = MPOSS2(i+n,2) - MPOSS2(i,2);
            btx = xx*bx + xy*by + xz*bz; bty = yx*bx + yy*by + yz*bz; btz = zx*bx + zy*by + zz*bz;
            brx = btx; bry = c*bty + s*btz; brz = -s*bty + c*btz;
            bxx = xx*brx + yx*bry + zx*brz; byy = xy*brx + yy*bry + zy*brz; bzz = xz*brx + yz*bry + zz*brz;
            MPOS2(i+n,0) = MPOS2(i,0) + bxx; MPOS2(i+n,1) = MPOS2(i,1) + byy; MPOS2(i+n,2) = MPOS2(i,2) + bzz;
    }
    """
    info = weave.inline(code,['mposs','mpos','m','rand','theta','numbeads'], headers=['<math.h>', '<stdlib.h>'])
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
	angle = pi - arccos((AB[0]*BC[0]+AB[1]*BC[1]+AB[2]*BC[2])/(dotAB*(BC[0]*BC[0]+BC[1]*BC[1]+BC[2]*BC[2]))**.5)
        jac = sin(angle-theta)/sin(angle)
        return mpos, jac

def canglebend(mposs, m, rand, theta):
        mpos = mposs.copy()
	jac = numpy.array([0.])
	numbeads = len(mposs)
        code = """
        int n, end;
        double c, s, dotAB, dotABBC, mag, angle;
        double xAB, yAB, zAB, xBC, yBC, zBC;
	double xx, xy, xz, yx, yy, yz, zx, zy, zz;
	double bx, by, bz, btx, bty, btz;
        double brx, bry, brz, bxx, byy, bzz;
        c = cos(theta);
	s = sin(theta);
	if (rand < .5){
		n = 1;
		end = numbeads - 1;
	}
	else{
		n = -1;
		end = 0;
	}	

	xAB = MPOS2(m,0) - MPOS2(m-n,0); yAB = MPOS2(m,1) - MPOS2(m-n,1); zAB = MPOS2(m,2) - MPOS2(m-n,2);
        xBC = MPOS2(m+n,0) - MPOS2(m,0); yBC = MPOS2(m+n,1) - MPOS2(m,1); zBC = MPOS2(m+n,2) - MPOS2(m,2);
	dotAB = xAB*xAB + yAB*yAB + zAB*zAB;
	dotABBC = xAB*xBC + yAB*yBC + zAB*zBC;
	mag = sqrt(dotAB);
	zx = xAB/mag; zy = yAB/mag; zz = zAB/mag;
	mag = dotABBC/dotAB;
	yx = xBC - mag*xAB; yy = yBC - mag*yAB; yz = zBC - mag*zAB;
	mag = sqrt(yx*yx+yy*yy+yz*yz);
	yx = yx/mag; yy = yy/mag; yz = yz/mag;
	xx = zy*yz - zz*yy; xy = zz*yx - zx*yz; xz = zx*yy - zy*yx;
	for( int i = m; i != end; i += n){
		bx = MPOSS2(i+n,0) - MPOSS2(i,0); by = MPOSS2(i+n,1) - MPOSS2(i,1); bz = MPOSS2(i+n,2) - MPOSS2(i,2);
		btx = xx*bx + xy*by + xz*bz; bty = yx*bx + yy*by + yz*bz; btz = zx*bx + zy*by + zz*bz;
		brx = btx; bry = c*bty + s*btz; brz = -s*bty + c*btz;
		bxx = xx*brx + yx*bry + zx*brz; byy = xy*brx + yy*bry + zy*brz; bzz = xz*brx + yz*bry + zz*brz;
		MPOS2(i+n,0) = MPOS2(i,0) + bxx; MPOS2(i+n,1) = MPOS2(i,1) + byy; MPOS2(i+n,2) = MPOS2(i,2) + bzz;
	}
	angle = M_PI - acos(dotABBC/sqrt(dotAB*(xBC*xBC+yBC*yBC+zBC*zBC)));
	JAC1(0) = sin(angle-(double)theta)/sin(angle);
	"""
	info = weave.inline(code,['mposs','mpos','m','rand','theta','jac','numbeads'], headers=['<math.h>', '<stdlib.h>'])
	return mpos, jac[0]

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
	#mpos=mpos12u3.copy()
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
	
def localmove_find(mpos123, m, rand, theta):
	"""Performs localmove_find with the reverse proximity criterion"""
	newcoord = localmove_find(mpos123, m, rand, theta)
	oldcoord = localmove_find(newcoord, m, rand, -theta)
	print oldcoord - mpos123

def localmove(mpos123, m, rand, theta):
	mpos = mpos123.copy()
	c = cos(theta)
	s = sin(theta)
	rotate = array([[1,0,0],[0,c,s],[0,-s,c]])
	if rand < .5:
		n = 1
	else:
		n=-1
	AB = mpos[m,:] - mpos[m-n,:]
	BC = mpos[m+n,:] - mpos[m,:]
	
	dotAB = AB[0]*AB[0] + AB[1]*AB[1] + AB[2]*AB[2]
	x = AB / dotAB**.5
        y = BC - (BC[0]*AB[0] + BC[1]*AB[1] + BC[2]*AB[2]) / dotAB * AB
        y = y / (y[0]*y[0] + y[1]*y[1] + y[2]*y[2])**.5
        z = [x[1]*y[2] - x[2]*y[1], x[2]*y[0] - x[0]*y[2], x[0]*y[1] - x[1]*y[0]]
        transform = array([x,y,z])
	untransform = transpose(transform)
	for i in range(m,m+n*3,n):
		bond = mpos123[i+n,:] - mpos123[i,:]
		bond = dot(transform,bond.transpose())
		bond = dot(rotate,bond)
		bond = dot(untransform,bond)
		mpos[i+n,:] = mpos[i,:] + bond
	AB = mpos[m+n*5,:] - mpos[m+n*3,:]
	BC = mpos[m+n*4,:] - mpos[m+n*3,:]
	dotAB = AB[0]*AB[0] + AB[1]*AB[1] + AB[2]*AB[2]
        x = AB / dotAB**.5
        y = BC - (BC[0]*AB[0]+BC[1]*AB[1]+BC[2]*AB[2]) / dotAB * AB
        y = y / (y[0]*y[0] + y[1]*y[1] + y[2]*y[2])**.5
        z = [x[1]*y[2] - x[2]*y[1], x[2]*y[0] - x[0]*y[2], x[0]*y[1] - x[1]*y[0]]
        transform=array([x,y,z])
	untransform=transpose(transform)
	l42=dot(mpos123[m+n*4,:]-mpos123[m+n*3,:],mpos123[m+n*4,:]-mpos123[m+n*3,:])
	l4=l42**.5
	l52=dot(mpos123[m+n*5,:]-mpos123[m+n*4,:],mpos123[m+n*5,:]-mpos123[m+n*4,:])	
	q2=dot(mpos[m+n*5,:]-mpos[m+n*3,:],mpos[m+n*5,:]-mpos[m+n*3,:])
	cosine=(l42+q2-l52)/(2*l4*q2**.5)
	if abs(cosine)>1:
		return nan
	else:
		theta = arccos(cosine)
                sine = sin(theta)
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
	maxloop=1000


	
	self.newcoord=self.coord.copy()
	self.vel=empty((numbeads,3))
	for i in range(numbeads): self.vel[i,:]=numpy.random.normal(0,(4.184*self.kb*self.T/m[i])**.5,3) #in nm/ps, uses average residue mass
	
        bonds=self.coord[0:numbeads-1,:]-self.coord[1:numbeads,:]
	d2=numpy.sum(bonds**2,axis=1)
	d=d2**.5
	force = HMCforce.cangleforces(self.coord, angleparam,bonds,d,numbeads) + HMCforce.cdihedforces(torsparam, bonds, d2, d, numbeads) + HMCforce.cnonbondedforces(self.coord,numint,numbeads,nativeparam_n,nonnativeparam,nnepsil)
	a = numpy.transpose(force) / m
	self.vel, conv = HMCforce.crattle(bonds, self.vel, m, d2, maxloop, numbeads, tol)
	self.vold = self.vel
        self.oldH=self.u0+.5/4.184*numpy.sum(m*numpy.sum(self.vel**2,axis=1)) # in kcal/mol
	
	temporary = False
	if temporary:
		Hplot = numpy.zeros(nsteps+1)
		Hplot[0] = self.oldH

	for e in range(nsteps):
		#finding r(t+dt)
		#loops = 0
		#conv = numpy.ones(numbeads-1)
		v_half = self.vel + h / 2 * numpy.transpose(a) # unconstrained v(t+dt/2)
		v_half, conv = HMCforce.cshake(bonds, v_half, h, m, d2, maxloop, numbeads, tol) # 
		if not conv:
			print 'MD not converging, reject'
			self.uncloseable=True
			self.rejected += 1
			break
		self.newcoord += h * v_half #constrained r(t+dt)
		bonds = self.newcoord[0:numbeads-1,:]-self.newcoord[1:numbeads,:] #rij(t+dt)
		force = HMCforce.cangleforces(self.newcoord, angleparam,bonds,d,numbeads) + HMCforce.cdihedforces(torsparam, bonds, d2, d, numbeads) + HMCforce.cnonbondedforces(self.newcoord,numint,numbeads,nativeparam_n,nonnativeparam,nnepsil)
		a=transpose(force)/m
		self.vel = v_half + h/2*transpose(a) # unconstrained v(t+dt)
		self.vel, conv = HMCforce.crattle(bonds, self.vel, m, d2, maxloop, numbeads, tol)
		if not conv:
			print 'MD not converging, reject'
			self.uncloseable=True
			self.rejected += 1
			break
		if temporary:
			coord=self.newcoord
			from energyfunc import *
    			r2=cgetLJr2(coord,numint,numbeads)
    			torsE=torsionenergy_nn(coord,zeros(numbeads-3),torsparam,arange(numbeads-3))
   			angE=angleenergy_n(coord,zeros(numbeads-2),angleparam,arange(numbeads-2))
   			u1=sum(angE)+sum(torsE)+LJenergy_n(r2,nativeparam_n,nonnativeparam,nnepsil)
			ke=.5*m/4.184*sum(self.vel**2,axis=1)
			Hplot[e+1] = u1 + sum(ke)
	
		#writetopdb.addtopdb(self.newcoord,positiontemplate,self.move*nsteps+e,'%s/trajectory%i.pdb' % (self.out,int(self.T)))
	#if(loops==maxloop):
		#self.uncloseable=True
		#self.rejected += 1
	#assert all(sum(bonds**2,axis=1)-d2< 1e-7)
	if temporary:
		import matplotlib.pyplot as plt
		plt.plot(range(nsteps+1),Hplot)
		plt.show()
	return self

def runMD_surf(self,nsteps,h,dict):
	numbeads=dict['numbeads']
	numint=dict['numint']
	angleparam=dict['angleparam']
	torsparam=dict['torsparam']
	nativeparam_n=dict['nativeparam_n']
	nonnativeparam=dict['nonnativeparam']
	nnepsil=dict['nnepsil']
	m=dict['mass']
        surface = dict['surface']
        surfparam = dict['surfparam']
        nspint = dict['nspint']
	positiontemplate=dict['positiontemplate']
	#m=120.368 # average mass of all residues
	tol=1e-8
	maxloop=1000
	
	self.newcoord=self.coord.copy()
	self.vel=empty((numbeads,3))
	for i in range(numbeads): self.vel[i,:]=numpy.random.normal(0,(4.184*self.kb*self.T/m[i])**.5,3) #in nm/ps, uses average residue mass
	
        bonds=self.coord[0:numbeads-1,:]-self.coord[1:numbeads,:]
	d2=numpy.sum(bonds**2,axis=1)
	d=d2**.5
        force = HMCforce.cangleforces(self.coord, angleparam,bonds,d,numbeads) + HMCforce.cdihedforces(torsparam, bonds, d2, d, numbeads) + HMCforce.cnonbondedforces(self.coord,numint,numbeads,nativeparam_n,nonnativeparam,nnepsil) + HMCforce.cgetsurfforce(self.coord, surface, nspint, numbeads, surfparam)
	a = numpy.transpose(force) / m
	self.vel, conv = HMCforce.crattle(bonds, self.vel, m, d2, maxloop, numbeads, tol)
        self.oldH=self.u0+.5/4.184*numpy.sum(m*numpy.sum(self.vel**2,axis=1)) # in kcal/mol
	
	for e in range(nsteps):
		#finding r(t+dt)
		#loops = 0
		#conv = numpy.ones(numbeads-1)
		v_half = self.vel + h / 2 * numpy.transpose(a) # unconstrained v(t+dt/2)
		v_half, conv = HMCforce.cshake(bonds, v_half, h, m, d2, maxloop, numbeads, tol) # 
		if not conv:
			self.uncloseable=True
			self.rejected += 1
			break
		self.newcoord += h * v_half #constrained r(t+dt)
		bonds = self.newcoord[0:numbeads-1,:]-self.newcoord[1:numbeads,:] #rij(t+dt)
		force = HMCforce.cangleforces(self.newcoord, angleparam,bonds,d,numbeads) + HMCforce.cdihedforces(torsparam, bonds, d2, d, numbeads) + HMCforce.cnonbondedforces(self.newcoord,numint,numbeads,nativeparam_n,nonnativeparam,nnepsil) + HMCforce.cgetsurfforce(self.coord, surface, nspint, numbeads, surfparam)
		a=transpose(force)/m
		self.vel = v_half + h/2*transpose(a) # unconstrained v(t+dt)
		self.vel, conv = HMCforce.crattle(bonds, self.vel, m, d2, maxloop, numbeads, tol)
		if not conv:
			self.uncloseable=True
			self.rejected += 1
			break
	return self



def runMD_noreplica(nsteps, h, coord, numbeads, numint, angleparam, torsparam, nativeparam_n, nonnativeparam, nnepsil, mass, kb, T, tol, maxloop):
	newcoord = coord.copy()
	vel = empty((numbeads,3))
	for i in range(numbeads): vel[i,:]=numpy.random.normal(0,(4.184*kb*T/mass[i])**.5,3) #in nm/ps, uses average residue mass
        bonds=coord[0:numbeads-1,:]-coord[1:numbeads,:]
	d2=numpy.sum(bonds**2,axis=1)
	d=d2**.5
	force = HMCforce.cangleforces(coord, angleparam, bonds, d, numbeads) + HMCforce.cdihedforces(torsparam, bonds, d2, d, numbeads) + HMCforce.cnonbondedforces(coord, numint, numbeads, nativeparam_n, nonnativeparam, nnepsil)
	a = numpy.transpose(force) / mass
	vel, conv = HMCforce.crattle(bonds, vel, mass, d2, maxloop, numbeads, tol)
        #self.oldH=self.u0+.5/4.184*numpy.sum(m*numpy.sum(self.vel**2,axis=1)) # in kcal/mol
	
	for e in range(nsteps):
		#finding r(t+dt)
		v_half = vel + h / 2 * numpy.transpose(a) # unconstrained v(t+dt/2)
		v_half, conv = HMCforce.cshake(bonds, v_half, h, mass, d2, maxloop, numbeads, tol) # 
		if not conv:
			return nan
		newcoord += h * v_half #constrained r(t+dt)
		bonds = newcoord[0:numbeads-1,:]-newcoord[1:numbeads,:] #rij(t+dt)
		force = HMCforce.cangleforces(newcoord, angleparam, bonds, d, numbeads) + HMCforce.cdihedforces(torsparam, bonds, d2, d, numbeads) + HMCforce.cnonbondedforces(newcoord, numint, numbeads, nativeparam_n, nonnativeparam, nnepsil)
		a = transpose(force)/mass
		vel = v_half + h/2*transpose(a) # unconstrained v(t+dt)
		vel, conv = HMCforce.crattle(bonds, vel, mass, d2, maxloop, numbeads, tol)
		if not conv:
			return nan
		#writetopdb.addtopdb(self.newcoord,positiontemplate,self.move*nsteps+e,'%s/trajectory%i.pdb' % (self.out,int(self.T)))
	return newcoord

def enddist(mpos):
    distvec=mpos[0,:]-mpos[-1,:]
    return dot(distvec,distvec)**.5

def makeT(c,s,cp,sp):
    """Generates transformation matrix from coordinate frame of bond i+1 to coordinate frame of bond i"""
    T = array([[-c,s,0],[-cp*s,-cp*c,-sp],[-sp*s,-sp*c,cp]])
    #T = array([[-cos(angle),sin(angle),0], [-cos(phi)*sin(angle),-cos(phi)*cos(angle),-sin(phi)], [-sin(phi)*sin(angle),-sin(phi)*cos(angle),cos(phi)]])
    return T

#def parrot_older(mpos123, m, rand, theta):
    #"""Performs a parallel-rotation move"""
    #mpos = mpos123.copy()
    #soln = ones((2,4)) # two possible sets of solutions
    #if rand > .5:
        #mpos = mpos[::-1]
    #angle = energyfunc.angle(mpos)
    #dihed = energyfunc.dihedral(mpos)
    #bonds = mpos[1:len(mpos),:]-mpos[0:-1,:]
    #phi_old=dihed[m-2:m+2]
    #soln[:,0]=dihed[m-2]+theta # new phi0
    
    ## find phi2_new
    #T0 = makeT(angle[m-1],soln[0,0]) # T(phi0_new)
    #try:
        #ulab = bonds[m+2,:]/sum(bonds[m+2,:]**2)**.5
    #except IndexError: # need fictitious bonds/angles
        #assert(m == len(mpos) - 3)
        #ulab = [0,0,1]
        #angle = append(angle,[68*pi/180, 68*pi/180])
        #phi_old = append(phi_old, [0,0])
    #jac_old = parrot_jac(ulab, bonds[m], bonds[m+1])

    #x=bonds[m-1,:]/sum(bonds[m-1,:]**2)**.5
    #z1=bonds[m-2,:]/sum(bonds[m-2,:]**2)**.5
    #z=cross(x,z1)/sin(pi-angle[m-2])
    #y=cross(z,x)
    #Tlab=vstack((x,y,z))
    #Tlab=transpose(Tlab) # transformation matrix from local coordinate frame of bond 0 to the laboratory frame
    #u = linalg.solve(Tlab,ulab) # find u in reference frame of bond 0
    #v = linalg.solve(T0,u) 
    #cosphi2 = (cos(angle[m])*cos(angle[m+1])-v[0]) / (sin(angle[m+1])*sin(angle[m]))
    #if abs(cosphi2) > 1:
        #return nan, nan # no solutions
    #soln[0,2] = arccos(cosphi2) 
    #soln[1,2] = 2*pi-arccos(cosphi2)
    
    ## find unique solution of phi1new
    #soln[0,1] = findphi1(soln[0,2], angle, m, v)
    #soln[1,1] = findphi1(soln[1,2], angle, m, v)

    ## test solutions 
    #try:
        #oldT = [makeT(angle[m-1],phi_old[0]), makeT(angle[m],phi_old[1]), makeT(angle[m+1],phi_old[2]), makeT(angle[m+2],phi_old[3])]
    #except IndexError:
	#assert(m == len(mpos) - 4)
	#angle = append(angle, 68*pi/180)
	#phi_old = append(phi_old, 0)
        #oldT = [makeT(angle[m-1],phi_old[0]), makeT(angle[m],phi_old[1]), makeT(angle[m+1],phi_old[2]), makeT(angle[m+2],phi_old[3])]
    #newT = [[makeT(angle[m-1],soln[0,0]), makeT(angle[m],soln[0,1]), makeT(angle[m+1],soln[0,2]), makeT(angle[m+2],soln[0,3])]
           #,[makeT(angle[m-1],soln[1,0]), makeT(angle[m],soln[1,1]), makeT(angle[m+1],soln[1,2]), makeT(angle[m+2],soln[1,3])]]
    
    #if not checku(oldT, newT[0]):
	#assert(m==len(mpos)-3)
	##pdb.set_trace()
	##checku(oldT, newT[0])
	##findphi1(soln[0,2], angle, m, v)
	##soln = delete(soln,s_[0],axis=0)
    #if not checku(oldT, newT[1]):
        #assert(m==len(mpos)-3)
	##soln = delete(soln,s_[1],axis=0)
    
    #pos = [findpos(soln[i,:], mpos, m, newT[i], Tlab, bonds) for i in range(len(soln))]

    #dihednew = [energyfunc.dihedral(mpos) for mpos in pos]
    #for i in range(len(dihednew)):
	#try:
	    #soln[i,3] = dihednew[i][m+1] # set new phi3
    	#except IndexError:
	    #assert(m==len(mpos)-3 or m==len(mpos)-4)
    #angnew = [energyfunc.angle(mpos) for mpos in pos]
    #for i,T in enumerate(newT):
	#try:
	    #T[3]=makeT(angle[m+2],soln[i,3])
 	#except IndexError:
	    #assert(m==len(mpos)-3 or m==len(mpos)-4)
	#if not checkuT(oldT,T):
	    #assert(m==len(mpos)-3 or m==len(mpos)-4)
	    ##soln = delete(soln, s_[i], axis=0)
	#checkdihed(dihed,dihednew[i])
	#checkang(angle,angnew[i])
    #jac_new = parrot_jac(ulab, mpos[m+1]-mpos[m], mpos[m+2]-mpos[m])
    #jac = jac_new/jac_old
    #if rand > .5:
	#pos[0]=pos[0][::-1]
    #return pos[0], jac

def parrot(mpos123, m, rand, theta):
    """Performs a parallel-rotation move"""
    mpos = mpos123.copy()
    soln = ones(3) # two possible sets of solutions
    if rand > .5:
        mpos = mpos[::-1]
    angle = energyfunc.angle(mpos, [m-2, m-1, m, m+1])
    c = cos(angle)
    s = sin(angle)
    dihed = energyfunc.dihedral(mpos, [m-2])
    bonds = mpos[1:len(mpos),:]-mpos[0:-1,:]
    soln[0]=dihed[0]+theta # new phi0
    
    # find phi2_new
    T0 = makeT(c[1],s[1],cos(soln[0]),sin(soln[0])) # T(phi0_new)
    try:
        ulab = bonds[m+2,:]/(bonds[m+2,0]**2+bonds[m+2,1]**2+bonds[m+2,2]**2)**.5 # sum(bonds[m+2,:]**2)**.5
    except IndexError: # need fictitious bonds/angles
        assert(m == len(mpos) - 3)
        ulab = [0,0,1]
        #angle = append(angle,[68*pi/180, 68*pi/180])
        c = append(c, [cos(68*pi/180)]*2)
        s = append(s, [sin(68*pi/180)]*2)
    jac_old = parrot_jac(ulab, bonds[m], bonds[m+1])
    try:
        dihed = energyfunc.dihedral(mpos, [m])
        if cos(dihed) > .999:
            Nmn = 1.
        else:
            Nmn = 2.
    except IndexError:
        Nmn = 2.
    x = bonds[m-1,:]/(bonds[m-1,0]**2+bonds[m-1,1]**2+bonds[m-1,2]**2)**.5
    z1 = bonds[m-2,:]/(bonds[m-2,0]**2+bonds[m-2,1]**2+bonds[m-2,2]**2)**.5
    z = array([x[1]*z1[2]-x[2]*z1[1],x[2]*z1[0]-x[0]*z1[2],x[0]*z1[1]-x[1]*z1[0]])/s[0]
    y = array([z[1]*x[2]-z[2]*x[1],z[2]*x[0]-z[0]*x[2],z[0]*x[1]-z[1]*x[0]])
    Tlab = array([[x[0],y[0],z[0]],[x[1],y[1],z[1]],[x[2],y[2],z[2]]])  # transformation matrix from local coordinate frame of bond 0 to the laboratory frame
    u = linalg.solve(Tlab,ulab) # find u in reference frame of bond 0
    v = linalg.solve(T0,u) 
    cosphi2 = (c[2]*c[3]-v[0]) / (s[3]*s[2])
    if abs(cosphi2) > 1:
        return nan, nan # no solutions
    elif abs(cosphi2) > .999:
        Nnm = 1.
    else:
        Nnm = 2.
    
    if random() < .5:
        soln[2] = arccos(cosphi2)
    else:
        soln[2] = 2*pi - arccos(cosphi2)
    
    # find unique solution of phi1new
    soln[1] = findphi1(soln[2], c, s, v)
    
    # test solutions 
    try:
        newT = [T0, makeT(c[2],s[2],cos(soln[1]),sin(soln[1])), makeT(c[3],s[3],cosphi2,sin(soln[2]))]
    except IndexError:
        assert(m == len(mpos) - 4)
	#angle = append(angle, 68*pi/180)
        c = append(c, [cos(68*pi/180)])
        s = append(s, [sin(68*pi/180)])
	newT = [T0, makeT(c[2],s[2],cos(soln[1]),sin(soln[1])), makeT(c[3],s[3],cosphi2,sin(soln[2]))]
    
    mpos = findpos(soln, mpos, m, newT, Tlab, bonds)
    jac_new = parrot_jac(ulab, mpos[m+1]-mpos[m], mpos[m+2]-mpos[m+1])
    jac = jac_old/jac_new*Nmn/Nnm
    if rand > .5:
	mpos=mpos[::-1]
    return mpos, jac


def parrot_jac(u, u1, u2):
    u1 = u1 / (u1[0]**2+u1[1]**2+u2[2]**2)**.5 # sum(u1**2)**.5
    u2 = u2 / (u2[0]**2+u2[1]**2+u2[2]**2)**.5 # sum(u2**2)**.5
    jac = array([u1[1]*u2[2]-u1[2]*u2[1],u1[2]*u2[0]-u1[0]*u2[2],u1[0]*u2[1]-u1[1]*u2[0]])
    jac = abs(u[0]*jac[0]+u[1]*jac[1]+u[2]*jac[2])
    return jac

def findphi1(phi2, c, s, v):
    a = s[2]*c[3] + cos(phi2)*c[2]*s[3]
    b = s[3]*sin(phi2)
    cosphi1 = (a*v[1] - b*v[2])/(1-v[0]**2)
    sinphi1 = (b*v[1] + a*v[2])/(1-v[0]**2)
    a = array([arccos(cosphi1), 2*pi - arccos(cosphi1)])
    b = array([arcsin(sinphi1), pi - arcsin(sinphi1)])
    b[b<0] += 2*pi
    decimal=14
    phi1=intersect1d(around(a,decimal),around(b,decimal))
    while not phi1.any() and decimal >  5:
	decimal -= 1
        phi1=intersect1d(around(a,decimal),around(b,decimal))
    try:
   	phi1=phi1[0]
    except IndexError:
	print a
	print b
	print phi1
	return nan
    return phi1

def checku(Told, Tnew):
    """Checks whether u = T0*T1*T2*i == T0n*T1n*T2n*i in ParRot move"""
    u_old = dot(dot(dot(Told[0],Told[1]),Told[2]),array([1,0,0]))
    u_new = dot(dot(dot(Tnew[0],Tnew[1]),Tnew[2]),array([1,0,0]))
    diff = u_old-u_new
    if sum(abs(diff)) < .000001:
        return True
    else:
        return False

def checkuT(Told, Tnew):
    """Checks whether uT = T0*T1*T2*T3*j == T0n*T1n*T2n*T3n*j in ParRot move"""
    uT_old = dot(dot(dot(dot(Told[0],Told[1]),Told[2]),Told[3]),array([0,1,0]))
    uT_new = dot(dot(dot(dot(Tnew[0],Tnew[1]),Tnew[2]),Tnew[3]),array([0,1,0]))
    diff = uT_old - uT_new
    if sum(abs(diff)) < .000001:
        return True
    else:
        return False

def findpos(phi, mpos123, m, T, Tlab, bonds):
    """Finds new bead positions from new phis for ParRot move"""
    mpos = mpos123.copy()
    bond=dot(Tlab,T[0])
    mpos[m+1,:] = mpos[m,:] + bond[:,0]*(bonds[m,0]**2+bonds[m,1]**2+bonds[m,2]**2)**.5 #dot(bond,array([sqrt(sum(bonds[m,:]**2)),0,0]))
    bond=dot(bond,T[1])
    mpos[m+2,:]=mpos[m+1,:] + bond[:,0]*(bonds[m+1,0]**2+bonds[m+1,1]**2+bonds[m+1,2]**2)**.5 # dot(bond,array([sqrt(sum(bonds[m+1,:]**2)),0,0]))
    bond=dot(bond,T[2])
    try:
	mpos[m+3,:]=mpos[m+2,:] + bond[:,0]*(bonds[m+2,0]**2+bonds[m+2,1]**2+bonds[m+2,2]**2)**.5 # dot(bond,array([sqrt(sum(bonds[m+2,:]**2)),0,0]))
    except IndexError:
	assert(m==len(mpos)-3)
    #T3=makeT(angle[m+2],dihed[m+1])
    #bond=dot(bond,T3)
    #bond=dot(bond,array([sqrt(sum(bonds[m+3,:]**2)),0,0]))
    #phi3=arccos(dot(bond,bonds[m+3])/(dot(bond,bond)*dot(bonds[m+3],bonds[m+3]))**.5)
    for i in range(m+3,len(mpos)-1):
        mpos[i+1,:] = mpos[i,:] + bonds[i,:]
    return mpos

def checkdihed(dihedold, dihednew):
    diff = dihedold - dihednew
    if sum(abs(diff) > .00000001) > 4:
	pdb.set_trace()
    return

def checkang(angold, angnew):
    try:
        diff = angold - angnew
    except ValueError:
	return
    if sum(abs(diff) > .00000001) > 0:
	pdb.set_trace()
    return

def translation(coord_old, x, y, z):
    """Used in SurfaceSimulation to translate protein"""
    coord = coord_old.copy()
    coord[:,0] += x
    coord[:,1] += y
    coord[:,2] += z
    return coord

def rotation(coord_old, alpha, beta, gamma):
    """Used in SurfaceSimulation to rotate protein"""
    coord = coord_old.copy()
    com = sum(coord, axis = 0)/len(coord) # center of mass
    coord -= com # rotation performed around CoM 
    ca = cos(alpha) 
    sa = sin(alpha)
    cb = cos(beta)
    sb = sin(beta)
    cg = cos(gamma)
    sg = sin(gamma)
    rotation_matrix = numpy.array([[cb*cg, cb*sg, -sb], [sa*sb*cg-ca*sg, sa*sb*sg+ca*cg, sa*cb], [ca*sb*cg+sa*sg, ca*sb*sg-sa*cg, ca*cb]])
    coord = dot(coord,rotation_matrix)
    coord = coord + com
    return coord
    
