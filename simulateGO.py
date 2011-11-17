from numpy import *
from writetopdb import *
from moveset import *
from energyfunc import *

T=300 #Kelvin
totmoves=100
energyarray=zeros(totmoves)

openfile=open('GO_protein.pdb','r')
i=0 # index for coordinate matrix
k=0 # index for line number
numbeads=57
coord=zeros((numbeads,3))# unhardcode this later
wtemp=[] #template of all words
ptemp=[] #template of position lines
hetatm=[] # array of location of ATOM lines, used to be HETATM
while 1:
    line = openfile.readline()
    if not line:#checks for end of file
        break
    wtemp.append(line)
    splitline=line.split('  ')
    if splitline[0]=='ATOM':
        ptemp.append(line)
        hetatm.append(k)
        n=len(splitline)
        coord[i][0]=float(line[31:38])
        coord[i][1]=float(line[39:46])
        coord[i][2]=float(line[47:54])
        i=i+1
    k=k+1

#Get parameters from GO_protein.param
file="GO_protein.param"
angleparam=getangleparam(file,numbeads)
torsparam=gettorsionparam(file,numbeads)
LJparam=getLJparam(file,numbeads)
nativeparam=getnativefix(file)

def energy(mpos):
	energy=angleenergy(mpos,angleparam)+torsionenergy(mpos,torsparam)+LJenergy(mpos,LJparam,nativeparam)
	return energy

u0=energy(coord)
print('0')
print(u0)
writeseqpdb(coord,wtemp,hetatm,0)
energyarray[0]=u0

for move in range(1,totmoves):
    print(move)
    while(1):        
        rand=random()
	if rand < .3333333:
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
    #addtopdb(newcoord,ptemp,move,filename)
    writeseqpdb(newcoord,wtemp,hetatm,move)
    coord=newcoord
    print(u1)
    u0=u1
    energyarray[move]=u0

#f=open(filename,'a')
#f.write('END\r\n')
#f.close


import matplotlib.pyplot as plt
plt.plot(range(totmoves),energyarray)
plt.xlabel('moves')
plt.ylabel('energy (kcal/mol)')
plt.title('Go-like model monte carlo simulation at '+str(T)+' K')
plt.show()
#plt.hist(randcoildist,500,normed=1)
#plt.show()


