from numpy import *
from crank import *

T=300 #Kelvin
totmoves=10
energyarray=zeros(totmoves)

openfile=open('halfGO_testgo.pdb','r')
i=0 # index for coordinate matrix
k=0 # index for line number
coord=zeros((374,3))#unhardcode this later
wtemp=[] #template of all words
ptemp=[] #template of position lines
hetatm=[]
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
            print('torsion')
	elif rand < .6666667:
            newcoord=reptation(coord)
	    print('snake')
	else:
	    newcoord=crankshaft(coord)
	    print('crank')
	u1=energy(newcoord)
	print(u1)
        if u1< u0:
            break
	if u1> u0:
	    break # temporary, curious
        kb=0.0019872041 #in kcal/mol
        boltz=exp(-u1/(kb*T))
        if random()<boltz:
            break
    #addtopdb(newcoord,ptemp,move,filename)
    writeseqpdb(newcoord,wtemp,hetatm,move)
    coord=newcoord
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


