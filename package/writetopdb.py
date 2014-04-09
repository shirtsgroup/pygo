import numpy
import pdb

'''This module contains various I/O helper methods'''

def get_coord(filename):
    f = file(filename,'r')
    # first count number of atoms to instantiate coord array
    n = 0
    for line in f:
        if 'ATOM' in line:
            n += 1
    coord = numpy.empty((n,3))
    pdb_text = [] # will hold the pdb file text
    i = 0
    f.seek(0)
    for line in f:
        if 'ATOM' in line:
            pdb_text.append(line)
            coord[i,0] = float(line[31:38])
            coord[i,1] = float(line[39:46])
            coord[i,2] = float(line[47:54])
            i += 1
    f.close()
    assert(i == n)
    return coord, pdb_text

def writeseqpdb(mpos,text,posline,move):
    '''Deprecated'''
    # multi file pdb output
    j=0
    write=''
    for i in posline:
        words=text[i,][0:30]
        coordstr=''
        coordstr=coordstr+str('%8.3f') % mpos[j,0]
        coordstr=coordstr+str('%8.3f') % mpos[j,1]
        coordstr=coordstr+str('%8.3f') % mpos[j,2]
        coordstr=coordstr+'\r\n'
        j+=1
        text[i]=words+coordstr
    f=file(str(move)+'.pdb','w')
    for k in range(len(text)): #don't want 'END'
        write=write+text[k]
    f.write(write)
    f.close

def writepdb(mpos,text,posline,move,filename):
    '''Deprecated'''
    # 1 file pdb output
    j=0
    for i in posline:
        words=[text[i][0:30],'%8.3f' %(mpos[j,0]),'%8.3f' %(mpos[j,1]),'%8.3f'%(mpos[j,2]),'\r\n']
        j=j+1
        text[i]="".join(words)
    f=file(filename,'w')
    f.write('MODEL        %i\r\n' % (move)) #check moves here
    write="".join(text[0:-1])
    f.write(write)
    f.write('ENDMDL\r\n')
    f.close

def addtopdb(mpos,coordtext,move,filename):
    '''Deprecated'''
    # 1 file pdb output
    for i in range(len(coordtext)):
        words=[coordtext[i][0:30],'%8.3f' %(mpos[i,0]),'%8.3f' %(mpos[i,1]),'%8.3f'%(mpos[i,2]),'\r\n']
	coordtext[i]="".join(words)
    f=file(filename,'a')
    f.write('MODEL        %i\r\n' % (move))
    write="".join(coordtext)
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
    '''Deprecated'''
    nc=nativecoord.copy()
    translate= nc[0,:]
    nc -= translate
    BC = nc[1,:]
    x1 = BC/numpy.dot(BC,BC)**.5
    AB = numpy.array([.5,.5,.5]); #random, but constant for all simulations
    y1 = AB-numpy.dot(AB,BC)/numpy.dot(BC,BC)*BC
    y1 = y1/numpy.sum(y1**2)**.5
    z1 = numpy.cross(x1,y1)
    return numpy.array([x1,y1,z1])

def getmovietransform_old(nativecoord):
    '''Deprecated'''
    nc = nativecoord.copy()
    center = len(nc)/2
    translate = nc[center,:]
    translate = translate.copy()
    for i in range(len(nc)):
		nc[i,:] -= translate
    BC = nc[center+1,:]
    x1 = BC/numpy.dot(BC,BC)**.5
    AB = [.5,.5,.5]; #random, but constant for all simulations
    y1 = AB-numpy.dot(AB,BC)/numpy.dot(BC,BC)*BC
    y1 = y1/numpy.dot(y1,y1)**.5
    z1 = numpy.cross(x1,y1)
    return [x1,y1,z1]

def moviecoord(mpos123,transform):
    '''Deprecated'''
    mpos = mpos123.copy()
    mpos[0,:] = numpy.zeros(3)
    bond = mpos123[1:len(mpos123),:]-mpos123[0:-1,:]
    bond = numpy.dot(bond,transform)
    for i in xrange(len(mpos)-1):
		mpos[i+1,:] = mpos[i,:]+bond[i,:]
    return mpos

def moviecoord_old(mpos123,transform):
    '''Deprecated'''
    mpos=mpos123.copy()
    center=len(mpos)/2
    translate=mpos[center,:]
    mpos-=translate
    for i in range(center,len(mpos)-1):
        BC=mpos123[i+1,:]-mpos123[i,:]
        BCnew=dot(transform,BC.transpose())
        mpos[i+1,:]=mpos[i,:]+BCnew
    for i in range(center,0,-1):
        BC=mpos123[i-1,:]-mpos123[i,:]
        BCnew=dot(transform,BC.transpose())
        mpos[i-1,:]=mpos[i,:]+BCnew
    return mpos
	
