from numpy import *
from crank import *

openfile=open('polyethConly.pdb','r')
i=0 # index for coordinate matrix
k=0 # index for line number
coord=zeros((8,3))#unhardcode this later
wtemp=[]
hetatm=[]
while 1:
    line = openfile.readline()
    if not line:
        break
    wtemp.append(line)
    splitline=line.split('  ')
    if splitline[0]=='HETATM':
        hetatm.append(k)
        n=len(splitline)
        coord[i][0]=float(splitline[n-3])
        coord[i][1]=float(splitline[n-2])
        word=splitline[n-1]
        coord[i][2]=float(word[:len(word)-2])
        i=i+1
    k=k+1

u0=energy(coord)
print(u0)

for move in range(1):
    print(move)
    if random() < .5:
    	coord=crankshaft(coord)
    	print('crank')
    else:
        coord=reptation(coord)
        print('reptation')
    writepdb(coord,wtemp,hetatm,move)
    u1=energy(coord)
    print(u1)
