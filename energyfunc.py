from numpy import *
from random import *

def getangleparam(paramfile):
    findline=0
    linestart=0
    for line in open(paramfile):
        findline += 1
        if "ANGLE" in line:
            linestart=findline
    
    
def getLJparam(paramfile)
    f=open(paramfile,'r')
    pass

def gettorsionparam(paramfile)
    f=open(paramfile,'r')
    pass

def LJenergy(mpos, param)

    pass

def angleenergy(mpos, param)
    pass

def torsionenergy(mpos, param)
    pass

def energy(mpos):
    energy=0.0; #potential energy
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

