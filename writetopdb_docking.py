import numpy

def get_coord(filename):
    f = file(filename,'r')
    # first count number of atoms to instantiate coord array
    n = 0
    for line in f:
        if 'ATOM' in line:
            n += 1
    coord = numpy.empty((n,3))
    pdb_text = [] # will hold the pdb file text
    i = 0
    f.seek(0)
    for line in f:
        if 'ATOM' in line:
            pdb_text.append(line)
            coord[i,0] = float(line[31:38])
            coord[i,1] = float(line[39:46])
            coord[i,2] = float(line[47:54])
            i += 1
    f.close()
    assert(i == n)
    return coord, pdb_text

def write_coord(filename,frame_num,coord,pdb_text):
    for i in range(len(pdb_text)):
        line = [pdb_text[i][0:30],'%8.3f' % coord[i,0],'%8.3f' % coord[i,1],'%8.3f' % coord[i,2],pdb_text[i][54::]]
        pdb_text[i] = "".join(line)
    if frame_num == 0 :
        f = file(filename,'w')
    else:
        f = file(filename,'a')
    f.write('MODEL %i\r\n' % frame_num)
    f.write("".join(pdb_text))
    f.write('ENDMDL\r\n')
    f.close()

def get_indices(filename,target_names):
    f = file(filename,'r')
    index = []
    for target in target_names: # target is a tuple of the string and its vwd radius
        i = 0
        for line in f:
            if 'ATOM' in line:
                if target[0] in line:
                    index.append((i,target[1])) 
                    break
                i += 1
        f.seek(0)
    assert(len(index) == len(target_names))
    return index

