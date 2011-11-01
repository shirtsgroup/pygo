from numpy import *
from crank import *

T=10000 #Kelvin

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
writepdb(coord,wtemp,hetatm,0)

for move in range(1,10):
    print(move)
    while(1):        
        if random() < .5:
            newcoord=crankshaft(coord)
        else:
            newcoord=reptation(coord)
	u1=energy(newcoord)
        if u1< u0:
            break
        kb=0.0019872041 #in kcal/mol
        boltz=exp(-u1/(kb*T))
        #print(boltz)
        if random()<boltz:
            break
    writepdb(newcoord,wtemp,hetatm,move)
    coord=newcoord
    u0=u1
    print(u0)



    
