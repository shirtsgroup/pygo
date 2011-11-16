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
    f=file('simulate.pdb','w')
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
        coordstr=coordstr+str('%8.3f') % mpos[i][0]
	coordstr=coordstr+str('%8.3f') % mpos[i][1]
	coordstr=coordstr+str('%8.3f') % mpos[i][2]
	coordstr=coordstr+'\r\n'
        coordtext[i]=words+coordstr
    f=file('simulate.pdb','a')
    for k in range(len(coordtext)):
        write=write+coordtext[k]
    f.write(write)
    f.write('ENDMDL\r\n')
    f.close



