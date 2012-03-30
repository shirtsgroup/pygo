from numpy import *
from random import *

def writeseqpdb(mpos,text,posline,move):
    # multi file pdb output
    j=0
    write=''
    for i in posline:
        words=text[i][0:30]
        coordstr=''
        coordstr=coordstr+str('%8.3f') % mpos[j][0]
        coordstr=coordstr+str('%8.3f') % mpos[j][1]
        coordstr=coordstr+str('%8.3f') % mpos[j][2]
        coordstr=coordstr+'\r\n'
        j=j+1
        text[i]=words+coordstr
    f=file(str(move)+'.pdb','w')
    for k in range(len(text)): #don't want 'END'
        write=write+text[k]
    f.write(write)
    f.close


def writepdb(mpos,text,posline,move,filename):
    # 1 file pdb output
    j=0
    for i in posline:
        words=text[i][0:30]
	coordstr=''
	coordstr=coordstr+str('%8.3f') % mpos[j][0]
	coordstr=coordstr+str('%8.3f') % mpos[j][1]
	coordstr=coordstr+str('%8.3f') % mpos[j][2]
	coordstr=coordstr+'\r\n'
        j=j+1
        text[i]=words+coordstr
    f=file(filename,'w')
    write='MODEL        '+str(move)+'\r\n' #check moves here
    for k in range(len(text)-1): #don't want 'END'
        write=write+text[k]
    f.write(write)
    f.write('ENDMDL\r\n')
    f.close

def addtopdb(mpos,coordtext,move,filename):
    # 1 file pdb output
    write='MODEL        '+str(move)+'\r\n'
    for i in range(len(coordtext)):
        words=coordtext[i][0:30]
        coordstr=''
        coordstr=coordstr+str('%8.3f') % mpos[i,0]
	coordstr=coordstr+str('%8.3f') % mpos[i,1]
	coordstr=coordstr+str('%8.3f') % mpos[i,2]
	coordstr=coordstr+'\r\n'
        coordtext[i]=words+coordstr
    f=file(filename,'a')
    for k in range(len(coordtext)):
        write=write+coordtext[k]
    f.write(write)
    f.write('ENDMDL\r\n')
    f.close

def addconnect(filename,k):
	#filename = .pdb file for linear chain polymer/protein without bonds
	#k = number of beads in polymer/protein
	f=open(filename,'a')
	text=''
	text=text+'CONECT    1    2\r\n'
	for i in range(2,k):
		text=text+'CONECT  '
		text=text+str('%3.0f') % i
    		text=text+'  '
    		text=text+str('%3.0f') % (i-1)
    		text=text+'  '
    		text=text+str('%3.0f') % (i+1)
    		text=text+'\r\n'
	text=text+'CONECT  '+str(k)+'  '+str(k-1)+'\r\nEND\r\n'
	f.write(text)
	f.close()
	
def getmovietransform(nativecoord):
	nc=nativecoord.copy()
	center=len(nc)/2
	translate= nc[center,:]
	translate=translate.copy()
	for i in range(len(nc)):
		nc[i,:] -= translate
	BC=nc[center+1,:]
	x1=BC/dot(BC,BC)**.5
	AB=[.5,.5,.5]; #random, but constant for all simulations
	y1=AB-dot(AB,BC)/dot(BC,BC)*BC
	y1=y1/dot(y1,y1)**.5
	z1=cross(x1,y1)
	return [x1,y1,z1]


def moviecoord(mpos123,transform):
	mpos=mpos123.copy()
	center=len(mpos)/2
	translate=mpos[center,:]
	#translate=translate.copy()
	for i in range(len(mpos)):
		mpos[i,:] -= translate
	for i in range(center,len(mpos)-1):
		BC=mpos123[i+1,:]-mpos123[i,:]
		BCnew=dot(transform,BC.transpose())
		mpos[i+1,:]=mpos[i,:]+BCnew
	for i in range(center,0,-1):
		BC=mpos123[i-1,:]-mpos123[i,:]
		BCnew=dot(transform,BC.transpose())
		mpos[i-1,:]=mpos[i,:]+BCnew
	return mpos
	
#def moviecoord(mpos123):
	#mpos=mpos123.copy()
	#center=len(mpos)/2
	#translate=mpos[center,:]
	#translate=translate.copy()
	#for i in range(len(mpos)):
		#mpos[i,:] -= translate
	#return mpos

