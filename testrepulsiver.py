# this script will test the two repulsive radius methods

from energyfunc import *
from datetime import *

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
    if splitline[0]=='ATOM':
        ptemp.append(line)
        hetatm.append(k)
        n=len(splitline)
        coord[i][0]=float(line[31:38])
        coord[i][1]=float(line[39:46])
        coord[i][2]=float(line[47:54])
        i=i+1
    k=k+1

file = "GO_protein.param"
nativeparam=getnativefix(file)

t1=datetime.now()
radius1=calcrepulsiver(coord, nativeparam)
t2=datetime.now()
radius2=calcrepulsiver2(coord,nativeparam)
t3=datetime.now()

print(radius1)
print(radius2)
print(radius1==radius2)
r1t=t2-t1
r2t=t3-t2
print(r1t)
print(r2t)