from numpy import *
from writetopdb import *
from moveset import *
from energyfunc import *

T=300 #Kelvin
totmoves=100
energyarray=zeros(totmoves)
filename='simulate.pdb'

openfile=open('polyethConly.pdb','r')
i=0 # index for coordinate matrix
k=0 # index for line number
coord=zeros((8,3))#unhardcode this later
wtemp=[] #template of all words
ptemp=[] #template of position lines
hetatm=[]
while 1:
    line = openfile.readline()
    if not line:#checks for end of file
        break
    wtemp.append(line)
    splitline=line.split('  ')
    if splitline[0]=='HETATM':
        ptemp.append(line)
        hetatm.append(k)
        n=len(splitline)
        coord[i][0]=float(splitline[n-3])
        coord[i][1]=float(splitline[n-2])
        word=splitline[n-1]
        coord[i][2]=float(word[:len(word)-2])
        i=i+1
    k=k+1

u0=energy(coord)
writepdb(coord,wtemp,hetatm,1,filename)
energyarray[0]=u0

for move in range(2,totmoves+1):
    print(move)
    while(1):        
        rand=random()
	if rand < 1:#.3333333:
            newcoord=torsion(coord)
        elif rand < .6666667:
            newcoord=reptation(coord)
	else:
	    newcoord=crankshaft(coord)
	u1=energy(newcoord)
        if u1< u0:
            break
        kb=0.0019872041 #in kcal/mol
        boltz=exp(-u1/(kb*T))
        if random()<boltz:
            break
    addtopdb(newcoord,ptemp,move,filename)
    coord=newcoord
    u0=u1
    energyarray[move-1]=u0

f=open(filename,'a')
f.write('END\r\n')
f.close


import matplotlib.pyplot as plt
plt.plot(range(totmoves),energyarray)
plt.xlabel('moves')
plt.ylabel('energy (kcal/mol)')
plt.title('8 monomer polyethylene monte carlo simulation at '+str(T)+' K')
plt.show()
#plt.hist(randcoildist,500,normed=1)
#plt.show()


