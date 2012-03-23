from numpy import *
from time import time
from random import *

t=time()
count=0
n=10000
numbeads=57
while count<n:
	dir=random()
	m=randint(1,numbeads-3)
	if dir<.5:
		change=[m-2]
	else:
		change=[m-1]
	change=change[change>-1]
	#print change
	count+=1
print time()-t

t=time()
count=0
n=10000
numbeads=57
while count<n:
	dir=random()
	m=randint(1,numbeads-3)
	if dir<.5:
		change=[m-2]
		if m<2:
			change=[]
	else:
		change=[m-1]
	#print change
	count+=1
print time()-t
	