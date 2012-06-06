from numpy import *
from energyfunc import *
import pdb

def bondedforces(mpos,torsparam,angleparam,bonds):
	forces=zeros((len(mpos),3))
	#torsion forces
	dihed = dihedral(mpos)
	for i in range(len(mpos)-3):
		rij=bonds[i]
		rkj=-bonds[i+1]
		rkl=bonds[i+2]
		m=cross(rij,rkj)
		n=cross(rkj,rkl)
		dV=-torsparam[4*i:4*i+4,0]*torsparam[4*i:4*i+4,1]*sin(torsparam[4*i:4*i+4,1]*dihed[i]-pi/180*torsparam[4*i:4*i+4,2])
		dV=sum(dV)
		Fi=-dV*sqrt(sum(rkj**2))
		Fl=Fi*n/sum(n**2)
		Fi=Fi*m/sum(m**2)
		Fj=dot(rij,rkj)/sum(rkj**2)*Fi-dot(rkl,rkj)/sum(rkj**2)*Fl
		Fk=-Fj-Fl
		Fj=-Fi+Fj
		forces[i,:]+=Fi
		forces[i+1,:]+=Fj
		forces[i+2,:]+=Fk	
		forces[i+3,:]+=Fl
	#angle forces
	angles = angle(mpos)
	for i in range(len(mpos)-2):
		rki=mpos[i,:]-mpos[i+2,:]
		rji=matrix(bonds[i])
		rjk=-matrix(bonds[i+1])
		dV=2*angleparam[i,0]*(angles[i]-pi/180*angleparam[i,1])
		fprime=-dV/sin(angles[i])/float(dot(rjk,rjk.T))*rki
		Fi=eye(3)-dot(rji.T,rji)/(rji[0,0]**2+rji[0,1]**2+rji[0,2]**2)
		Fi=dot(Fi,fprime)
		Fk=eye(3)-dot(rjk.T,rjk)/(rjk[0,0]**2+rjk[0,1]**2+rjk[0,2]**2)
		Fk=-dot(Fk,fprime)
		Fj=-Fi-Fk
		forces[i,:]+=array(Fi)[0,:]
		forces[i+1,:]+=array(Fj)[0,:]
		forces[i+2,:]+=array(Fk)[0,:]
	return forces*4.184

def nonbondedforces(mpos,numint,numbeads,natparam,nonnatparam,nnepsil):
	forces=zeros((len(mpos),3))
	#nonbonded forces
	rvec=getforcer(mpos,numint,numbeads)
	r2=sum(rvec**2, axis=1)
	rr=r2**.5
	ndV=natparam[:,0]*natparam[:,2]*natparam[:,2]/r2
	ndV6=ndV*ndV*ndV
	ndV=natparam[:,1]*(-156*ndV6*ndV6/rr+180*ndV6*ndV*ndV/rr-24*ndV6/rr)
	nndV=nonnatparam[:,0]*nonnatparam[:,1]*nonnatparam[:,1]/r2
	nndV=nndV*nndV*nndV
	nndV=-12*nnepsil*nndV*nndV/rr	
	k=0
	F=-(ndV+nndV)/rr
	for i in range(len(F)): rvec[i,:]=F[i]*rvec[i,:]
	F=rvec
	for i in range(len(mpos)):
		for j in range(i+3,len(mpos)):
			forces[i,:]+=F[k]
			forces[j,:]+=-F[k]
			k+=1
	return forces*4.184


def runMD(mpos123,tstep,nsteps,T,numbeads,numint,natparam,nonnatparam,nnepsil,angleparam,torsparam):
	mpos=mpos123.copy()
	v=random.normal(0,sqrt(.0083144621*T/120.3684),(numbeads,3)) #in nm/ps, uses average residue mass
	f=bondedforces(mpos,torsparam,angleparam)+nonbondedforces(mpos,numint,numbeads,natparam,nonnatparam,nnepsil)
	f=f*4.184 # convert from kcal to kJ
	a=f/120.368 # use average residue mass
	for i in range(nsteps):
		mpos_new=mpos+v*tstep+.5*a*tstep**2
                fnew=bondedforces(mpos_new,torsparam,angleparam)+nonbondedforces(mpos_new,numint,numbeads,natparam,nonnatparam,nnepsil)
		fnew=fnew*4.184 # convert from kcal to kJ
                anew=fnew/120.368 # use average residue mass
		vnew=v+.5*tstep*(anew+a)
		v=vnew
		a=anew
		mpos=mpos_new
	return mpos
	
def rattle(mpos,force,h,vel,d2,bonds):
	m=120.368 # average mass of all residues
	tol=1e-10
	
	# finding r(t+dt)
	conv=ones(len(mpos)-1)
	q=vel+h/(2*m)*force
	while sum(conv)!=0:
		for i in range(len(mpos)-1): # loops through all bonds/constraints
			s=bonds[i]+h*(q[i,:]-q[i+1,:])
			diff=sum(s**2)-d2[i]
			if (diff<tol):
				conv[i]=0
			else:
				g=diff/(h*dot(s,bonds[i]))*m
				q[i,:]-=g/m*bonds[i]
				q[i+1,:]+=g/m*bonds[i]
	mpos+=h*q #r(t+dt)
	vel=q # new v(t)
	bonds=bond(mpos) #rij(t+dt)
	#finding v(t+dt)
	conv=ones(len(mpos)-1)
	vel+=h/(2*m)*force
	while sum(conv)!=0:
		for i in range(len(mpos)-1):
			if (dot(bonds[i,:],vel[i,:]-vel[i+1,:])<tol):
				conv[i]=0
			else:
				k=dot(bonds[i,:],(vel[i,:]-vel[i+1,:])/d2[i]*2*m)
				vel[i] -= k/m*bonds[i,:]
				vel[i+1] += k/m*bonds[i,:]
	return [mpos,vel]	
