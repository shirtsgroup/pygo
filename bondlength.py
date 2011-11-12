from numpy import *
openfile=open('GO_protein.pdb','r')
mpos=zeros((57,3))
dist=zeros(56)
i=0
while 1:
	line=openfile.readline()
	if not line:
		break
	if line[0:4]=='ATOM':
		mpos[i][0]=line[31:38]
		mpos[i][1]=line[39:46]
		mpos[i][2]=line[47:54]
		i=i+1
for i in range(len(mpos)-1):
	distvec=mpos[i+1][:]-mpos[i][:] #this is hardcoded
	dist[i]=dot(distvec,distvec)**.5
print(dist)
print(mean(dist))

	
