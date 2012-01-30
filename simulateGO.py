from numpy import *
from writetopdb import *
from moveset import *
from energyfunc import *
from datetime import datetime
t1=datetime.now()

T=250 #Kelvin
totmoves=100
energyarray=zeros(totmoves)

openfile=open('GO_protein.pdb','r')
i=0 # index for coordinate matrix
k=0 # index for line number
numbeads=57# unhardcode this later
coord=zeros((numbeads,3)) #3 dimensional
wtemp=[] #template of all words
ptemp=[] #template of position lines
hetatm=[] # array of location of ATOM lines, used to be HETATM

# fills array of position vectors
while 1:
    line = openfile.readline()
    if not line:#checks for end of file
        break
    wtemp.append(line)
    splitline=line.split('  ')
    if splitline[0]=='ATOM': #too hard coded?
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
	radii=calcrepulsiver(mpos,nativeparam)
	energy=angleenergy(mpos,angleparam)+torsionenergy(mpos,torsparam)+LJenergy(mpos,LJparam,nativeparam,radii)
	return energy

u0=energy(coord)
print('0')
print(u0)
writeseqpdb(coord,wtemp,hetatm,0)
energyarray[0]=u0
torsmoves=0 #temporary, to check what kind of moves are accepted
reptmoves=0
crankmoves=0
atormoves=0
acceptedat=0
acceptedt=0
acceptedr=0
acceptedc=0
accepted=0
rejected=0
movetype='t'
for move in range(1,totmoves):
    print(move)
    while(1):        
        rand=random()
	if rand < .01:
            newcoord=torsion(coord)
	    torsmoves += 1
	    movetype='t'
	elif rand < 0:
            newcoord=reptation(coord)
	    reptmoves += 1
	    movetype='r'
	elif rand < .85:
	    newcoord=axistorsion(coord)
	    movetype='at'
	    atormoves += 1
	else:
	    newcoord=crankshaft(coord)
	    crankmoves += 1
	    movetype='c'
	u1=energy(newcoord)
        if u1< u0:
            accepted += 1
	    break
        kb=0.0019872041 #in kcal/mol/K
        boltz=exp(-u1/(kb*T))
	print(boltz)
        if random()<boltz:
            accepted += 1
	    break
	rejected += 1
    if movetype=='t':
	acceptedt += 1
    elif movetype=='r':
	acceptedr += 1
    elif movetype=='at':
	acceptedat +=1
    else: 
        acceptedc += 1
    
    
    #addtopdb(newcoord,ptemp,move,filename)
    writeseqpdb(newcoord,wtemp,hetatm,move)
    coord=newcoord
    print(u1)
    u0=u1
    energyarray[move]=u0

print('bend/torsion: '+str(acceptedt) +' out of '+str(torsmoves))
print('reptation: '+str(acceptedr) +' out of '+str(reptmoves))
print('crankshaft: '+str(acceptedc) +' out of '+str(crankmoves))
print('twist/torsion: '+str(acceptedat)+' out of '+str(atormoves))
print('accepted: '+str(accepted))
print('rejected: '+str(rejected))

#import matplotlib.pyplot as plt
#plt.plot(range(totmoves),energyarray)
#plt.xlabel('moves')
#plt.ylabel('energy (kcal/mol)')
#plt.title('Go-like model monte carlo simulation at '+str(T)+' K')
#plt.show()
#plt.hist(randcoildist,500,normed=1)
#plt.show()

t2=datetime.now()
print(t2-t1)
