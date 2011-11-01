from numpy import *
from random import *

def crankshaft(mpos): #moleculeposition, will need bonds later
    m=randint(1,6) #random molecule, not end ones

    posb=mpos[m][:] # middle molecule, will get crankshafted
    posa=mpos[m-1][:] #one 'before' it
    posc=mpos[m+1][:] # one 'after' it

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
    AB1=dot(transform,AB.transpose())
    #print(AB1)
    if randint(1,3) ==1:    
	#180 deg crank
        AB1[1]=-AB1[1]
    elif randint(1,3)==2:
	#90 deg crank
	AB1[1],AB1[2]=AB1[2],AB1[1]
    else:
        #90 deg other direction
        AB1[1],AB1[2]=AB1[2],AB1[1]*-1
    AB=dot(untransform,AB1)
    mpos[m][:]=posa+AB
    return mpos

def reptation(mpos):
    theta=pi*random()
    phi=2*pi*random()
    rho=1.53710311951 #unhardcode this later
    vec=[rho*sin(theta)*cos(phi),rho*sin(theta)*sin(phi),rho*cos(theta)]    
    n=len(mpos)
    if random() < .5:    
	for i in range(n-1):
             mpos[i][:]=mpos[i+1][:]
    	mpos[n-1][:]=mpos[n-2][:]+vec
    else:
	for i in range(n-1,0,-1):
	     mpos[i][:]=mpos[i-1][:]
	mpos[0][:]=mpos[1][:]+vec
    return mpos
    
def writepdb(mpos,text,posline,move):
    j=0
    for i in posline:
        words=text[i].split('  ')
        words[len(words)-3]=str('%6.3f') % mpos[j][0]
        words[len(words)-2]=str('%6.3f') % mpos[j][1]
        words[len(words)-1]=str('%6.3f') % mpos[j][2]+'\r\n'
        j=j+1
        text[i]='  '.join(words)
    f=file(str(move)+'.pdb','w')
    write=''
    for k in range(len(text)):
        write=write+text[k]
    f.write(write)
    f.close

def energy(mpos): #need bond for more complex molecules
    energy=0.0; #potential energy
    sig=4.6 #angstroms for polyethylene
    e= .42 #kcal/mol for polyethylene
    index=arange(len(mpos))
    for i in index:
        low=index[index<i-2]
        high=index[index>i+2]
        vdw=append(low,high) #index of beads excluding 12 and 13 neighbors
        for j in vdw:
            r=((mpos[i][0]-mpos[j][0])**2+(mpos[i][1]-mpos[j][1])**2+(mpos[i][2]-mpos[j][2])**2)**.5
            energy=energy+4*e*((sig/r)**12-(sig/r)**6)
    return energy/2

