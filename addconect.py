f=open('GO_protein.pdb','a')
text=''
k=57
text=text+'CONECT    1    2\r\n'
for i in range(2,k):
    text=text+'CONECT  '
    text=text+str('%3.0f') % i
    text=text+'  '
    text=text+str('%3.0f') % (i-1)
    text=text+'  '
    text=text+str('%3.0f') % (i+1)
    text=text+'\r\n'
text=text+'CONECT  57  56\r\nEND\r\n'
f.write(text)
f.close()

