from scipy import weave
import numpy
import random
import pdb
#==========================================
# PARAMETER READ IN METHODS
#==========================================

def getangleparam(paramfile, numbeads):
    f = open(paramfile, 'r')
    param = numpy.empty((numbeads-2, 2))
    while 1:
	line = f.readline()
        if "ANGLE" in line:
            break
    for i in xrange(numbeads-2): #no angle formed for end beads
	line = f.readline()
	param[i,0] = float(line[26:36]) #somewhat hardcoded, may apply to all/most go model param files?
	param[i,1] = numpy.pi/180*float(line[37:-1])
    f.close()
    return param

def getmass(topfile, numbeads):
    f = open(topfile, 'r')
    mass = numpy.empty(numbeads)
    while 1:
        line = f.readline()
        if "MASS" in line:
            break
    for i in xrange(numbeads):
        mass[i] = float(line[18:-1])
        line = f.readline()
    f.close()
    return mass

def gettorsionparam_old(paramfile, numbeads):
    f=open(paramfile, 'r')
    param = numpy.empty((4*(numbeads-3), 3)) #4 torsional potentials per 4 molecules
    while 1:
	line = f.readline()
        if "DIHEDRAL" in line:
            break
    for i in xrange(4*(numbeads-3)):
	line = f.readline()
	param[i,0] = float(line[22:30]) #somewhat hardcoded
	param[i,1] = float(line[32:33])
	param[i,2] = float(line[35:45])
    f.close()
    return param

def gettorsionparam(paramfile, numbeads):
    f = open(paramfile, 'r')
    param = numpy.empty((4*(numbeads-3),3)) #4 torsional potentials per 4 molecules
    while 1:
        line = f.readline()
        if "DIHEDRAL" in line:
            break
    i = 0
    while 1:
        line = f.readline()
        if line[1:4] == '':
                break
        if i/4+1 == int(line[1:4]):
            param[i,0] = float(line[22:30]) #somewhat hardcoded
            param[i,1] = float(line[32:33])
            param[i,2] = numpy.pi / 180 * float(line[35:45])
        else:
            param[i,:]=[0.0, 0.0, 0.0]
            f.seek(-len(line), 1)
            print 'Dihedral parameter missing'
        i += 1
    f.close()
    return param
	
# speed up version
def getLJparam_n(paramfile, numbeads, numint):
    f = open(paramfile, 'r')
    param = numpy.empty(numbeads) #two nonzero parameters, ignored first column (all 0)
    while 1:
        line = f.readline()
        if "NONBONDED" in line:
            f.readline() # two lines between header and parameters
            f.readline()
            break
    for i in xrange(numbeads): # param for each bead, use combining rules for interaction
        line = f.readline()
        epsil = -float(line[14:23]) #somewhat hardcoded
        param[i] = float(line[25:33])
    f.close()
    sigarray = numpy.zeros(numint)
    index = numpy.arange(numbeads)
    k=0
    for i in index:
        vdw = index[index > i+2]
        for j in vdw:
            sigarray[k] = param[i] + param[j]
            k += 1
    return [sigarray, epsil]

#speed up version
def getnativefix_n(paramfile, numint, numbeads):
    param = numpy.zeros((numint, 3))
    f = open(paramfile, 'r')
    while 1:
        line = f.readline()
        if "NBFIX" in line:
            break
    while 1:
        line = f.readline()
        if not line:
            break
        if "G" in line:
            [i, j, ep, sig] = [int(line[1:4]), int(line[9:12]), float(line[19:28]), float(line[32:-1])]
            intindex = numpy.sum(numbeads-numpy.arange(1,i)) + (j-i) - 1 - 2 * i
            param[intindex,:] = numpy.array([1,-ep,sig])
    return param

#==========================================
# ENERGY CALCULATION METHODS
#==========================================

def getforcer(mpos, numint, numbeads):
    r2array = numpy.empty((numint,3)) # distance vectors
    k = 0
    for i in xrange(numbeads):
        BC = mpos[i,:] - mpos[i+3:numbeads,:]
        knew = k + numbeads - (i + 3)
        r2array[k:knew,:] = BC
        k = knew
    return r2array #r^2 values for every interaction

def getLJr2(mpos, numint, numbeads):
    r2array = numpy.empty(numint)
    k = 0
    for i in xrange(numbeads):
        BC = mpos[i,:] - mpos[i+3:numbeads,:]
        knew = k + numbeads - (i + 3)
        r2array[k:knew] = numpy.sum(BC**2, axis=1)
        k = knew
    return r2array #r^2 values for every interaction

def cgetLJr2(mpos, numint, numbeads):
    r2array = numpy.empty(numint)
    code = """
    int k = 0;
    double x, y, z;
    for ( int i = 0; i < numbeads; i++){
        for ( int j = i+3; j < numbeads; j++){
	    x = MPOS2(i,0) - MPOS2(j,0);
	    y = MPOS2(i,1) - MPOS2(j,1);
	    z = MPOS2(i,2) - MPOS2(j,2);
	    R2ARRAY1(k) = x*x + y*y + z*z;
	    k++;
	}
    }
    """
    info = weave.inline(code, ['mpos','numint','numbeads','r2array'], headers=['<math.h>', '<stdlib.h>'])
    return r2array

def cgetforcer(mpos, numint, numbeads):
    r2array = numpy.empty((numint,3))
    code = """
    int k = 0;
    double x, y, z;
    for ( int i = 0; i < numbeads; i++){
        for ( int j = i+3; j < numbeads; j++){
	    R2ARRAY2(k,0) = MPOS2(i,0) - MPOS2(j,0);
	    R2ARRAY2(k,1) = MPOS2(i,1) - MPOS2(j,1);
	    R2ARRAY2(k,2) = MPOS2(i,2) - MPOS2(j,2);
	    k++;
	}
    }
    """
    info = weave.inline(code, ['mpos','numint','numbeads','r2array'], headers=['<math.h>', '<stdlib.h>'])
    return r2array

#speed up version
def LJenergy(r2, natparam, nonnatparam, nnepsil):
    #native calculation
    nE = natparam[:,0] * natparam[:,2] * natparam[:,2] / r2 #sigma2/r2
    nE6 = nE * nE * nE
    nE = natparam[:,1] * (13 * nE6 * nE6 - 18 * nE6 * nE * nE + 4 * nE6)
    #nonnative calculation
    nnE = nonnatparam[:,0] * nonnatparam[:,1] * nonnatparam[:,1] / r2 #simga2/r2
    nnE6 = nnE * nnE * nnE
    nnE = nnepsil * (13 * nnE6 * nnE6 - 18 * nnE6 * nnE * nnE + 4 * nnE6)
    energy = numpy.sum(nE) + numpy.sum(nnE)
    return energy

def LJenergy_n(r2, natparam, nonnatparam, nnepsil):
    #native calculation
    nE = natparam[:,0] * natparam[:,2] * natparam[:,2] / r2 #sigma2/r2
    nE6 = nE * nE * nE
    nE = natparam[:,1] * (13 * nE6 * nE6 - 18 * nE6 * nE * nE + 4 * nE6)
    #nonnative calculation
    nnE = nonnatparam[:,0] * nonnatparam[:,1] * nonnatparam[:,1] / r2 #simga2/r2
    nnE = nnE * nnE #sigma4/r4
    nnE = nnepsil * nnE * nnE * nnE
    energy = numpy.sum(nE) + numpy.sum(nnE)
    return energy

def LJenergy_CHARMM(r2, natparam, nonnatparam, nnepsil):
    #native calculation
    nE = natparam[:,0] * natparam[:,2] * natparam[:,2] / r2 #sigma2/r2
    nE6 = nE * nE * nE
    nE = natparam[:,1] * (nE6 * nE6 - 2 * nE6)
    #nonnative calculation
    nnE = nonnatparam[:,0] * nonnatparam[:,1] * nonnatparam[:,1] /r2 #simga2/r2
    nnE = nnE * nnE * nnE #sigma6/r6
    nnE = nnepsil * (nnE * nnE - 2 * nnE)
    energy = numpy.sum(nE)+numpy.sum(nnE)
    return energy
    
def angleenergy_n(mpos, oldE, param, change):
    newE = oldE.copy()
    for i in change:
        ktheta = param[i,0] # param file goes from 0 to len(mpos)-2
	optangle = param[i,1]
	BA = mpos[i,:] - mpos[i+1,:]
        BC = mpos[i+2,:] - mpos[i+1,:]
	dotBABC = BA[0]*BC[0]+BA[1]*BC[1]+BA[2]*BC[2]
	dotBA = BA[0]*BA[0]+BA[1]*BA[1]+BA[2]*BA[2]
	dotBC = BC[0]*BC[0]+BC[1]*BC[1]+BC[2]*BC[2]
        angle = numpy.arccos(dotBABC / (dotBA * dotBC)**.5) #in radians
        newE[i] = ktheta * (angle - optangle)**2
    #print('angle energy: '+str(energy))
    return newE
    
def torsionenergy_nn(mpos, oldE, param, change):
    newE=oldE.copy()
    for i in change:
        AB = mpos[i+1,:] - mpos[i,:]
        BC = mpos[i+2,:] - mpos[i+1,:]
	CD = mpos[i+3,:] - mpos[i+2,:]
        plane1 = numpy.array([AB[1]*BC[2]-AB[2]*BC[1], AB[2]*BC[0]-AB[0]*BC[2], AB[0]*BC[1]-AB[1]*BC[0]])  #cross(AB,BC)
        plane2 = numpy.array([BC[1]*CD[2]-BC[2]*CD[1], BC[2]*CD[0]-BC[0]*CD[2], BC[0]*CD[1]-BC[1]*CD[0]]) #cross(CD,BC)
        AB = (BC[0]**2+BC[1]**2+BC[2]**2)**.5 * AB
        a = AB[0]*plane2[0] + AB[1]*plane2[1] + AB[2]*plane2[2]
        b = plane1[0]*plane2[0] + plane1[1]*plane2[1] + plane1[2]*plane2[2]
        dihedral = numpy.arctan2(a,b)
        if dihedral < 0:
           dihedral += 2*numpy.pi
        energy = param[4*i:4*i+4,0] * (1 + numpy.cos(param[4*i:4*i+4,1]*dihedral-param[4*i:4*i+4,2]))
	newE[i] = energy[0]+energy[1]+energy[2]+energy[3]
    return newE

def ctorsionenergy(mpos, oldE, param, change):
    newE = oldE.copy()
    n = len(change)
    change = numpy.array(change)
    #bonds = mpos[0:-1,:]-mpos[1:len(mpos),:]
    code = """
    int i;
    double x1, x2, x3, y1, y2, y3, z1, z2, z3;
    double mx, my, mz, nx, ny, nz;
    double magBC, a, b, dihed;
    for ( int j = 0; j < n; j++){
        i = CHANGE1(j);
        x1 = MPOS2(i,0) - MPOS2(i+1,0); y1 = MPOS2(i,1) - MPOS2(i+1,1); z1 = MPOS2(i,2) - MPOS2(i+1,2);
        x2 = MPOS2(i+2,0) - MPOS2(i+1,0); y2 = MPOS2(i+2,1) - MPOS2(i+1,1); z2 = MPOS2(i+2,2) - MPOS2(i+1,2);
        x3 = MPOS2(i+2,0) - MPOS2(i+3,0); y3 = MPOS2(i+2,1) - MPOS2(i+3,1); z3 = MPOS2(i+2,2) - MPOS2(i+3,2);
        mx = y1*z2 - z1*y2; my = z1*x2 - x1*z2; mz = x1*y2 - y1*x2;
        nx = y2*z3 - z2*y3; ny = z2*x3 - x2*z3; nz = x2*y3 - y2*x3;
        magBC = sqrt(x2*x2 + y2*y2 + z2*z2);
        a = magBC*(x1*nx + y1*ny + z1*nz);
        b = mx*nx + my*ny + mz*nz;
        dihed = atan2(a,b);
        if ( dihed < 0){
            dihed += 2*M_PI;
        }
        NEWE1(i) = PARAM2(4*i,0)*(1+cos(PARAM2(4*i,1)*dihed-PARAM2(4*i,2))) + PARAM2(4*i+1,0)*(1+cos(PARAM2(4*i+1,1)*dihed-PARAM2(4*i+1,2))) + PARAM2(4*i+2,0)*(1+cos(PARAM2(4*i+2,1)*dihed-PARAM2(4*i+2,2))) + PARAM2(4*i+3,0)*(1+cos(PARAM2(4*i+3,1)*dihed-PARAM2(4*i+3,2)));
    }
    """
    info = weave.inline(code, ['param', 'newE', 'mpos', 'change','n'], headers=['<math.h>', '<stdlib.h>'])
    return newE

#used in simulatepolymer
def polymerenergy(mpos):
    energy = 0.0 # potential energy
    sig = 4.6 # angstroms for polyethylene
    e = .42 # kcal/mol for polyethylene
    ktheta = .82 # kcal/mol
    A = 5.22 # torsional parameter
    B = 2.88 # torsional parameter
    C = 1.95 # torsional parameter
    index = numpy.arange(len(mpos))
    # 6-12 LJ potential 
    for i in index:
        low = index[index<i-2]
        high = index[index>i+2]
        vdw = numpy.append(low, high) #index of beads excluding 12 and 13 neighbors
        for j in vdw:
            r = ((mpos[i][0]-mpos[j][0])**2+(mpos[i][1]-mpos[j][1])**2+(mpos[i][2]-mpos[j][2])**2)**.5
            energy = energy + 2*e*((sig/r)**12-(sig/r)**6) #divided by two since count each interaction twice
    # angle potential
    for i in range(1, len(mpos)-1):
        BA = mpos[:][i-1] - mpos[:][i]
        BC = mpos[:][i+1] - mpos[:][i]
        angle = numpy.arccos(dot(BA,BC)/(dot(BA,BA)**.5*dot(BC,BC)**.5)) #in radians
        energy = energy + ktheta / 2 * (angle-numpy.pi)**2
    # torsional potential
    for i in range(len(mpos)-3):
        AB = mpos[:][i+1] - mpos[:][i]
        BC = mpos[:][i+2] - mpos[:][i+1]
        CD = mpos[:][i+3] - mpos[:][i+2]
        plane1 = numpy.cross(BC,AB)
        plane2 = numpy.cross(CD,BC)
        dihedral  = numpy.arccos((plane1[0]*plane2[0]+plane1[1]*plane2[1]+plane1[2]*plane2[2]) / ((plane1[0]**2+plane1[1]**2+plane1[2]**2)**.5*(plane2[0]**2+plane2[1]**2+plane2[2]**2)**.5))
        energy = energy + A + B * numpy.cos(dihedral) + C * numpy.cos(3*dihedral)
    return energy

#==========================================
# CONFIGURATION ANALYSIS METHODS
#==========================================
def rmsd(crds1, crds2):
  	"""Returns RMSD between 2 sets of [nx3] numpy array"""
 	assert(crds1.shape[1] == 3)
 	assert(crds1.shape == crds2.shape)
 	n_vec = numpy.shape(crds1)[0]
 	correlation_matrix = numpy.dot(numpy.transpose(crds1), crds2)
 	v, s, w_tr = numpy.linalg.svd(correlation_matrix)
 	is_reflection = (numpy.linalg.det(v) * numpy.linalg.det(w_tr)) < 0.0
 	if is_reflection:
  		s[-1] = - s[-1]
  	E0 = sum(sum(crds1 * crds1)) + sum(sum(crds2 * crds2))
  	rmsd_sq = (E0 - 2.0*sum(s)) / float(n_vec)
  	rmsd_sq = max([rmsd_sq, 0.0])
 	return rmsd_sq**.5
	
def nativecontact(r2,nativeparam,nsigma2):
	# native contact if ij pair is native and if rij < 1.2 sigmaij
	r2 = r2 * nativeparam[:,0]
	nc = nsigma2 - r2
	nc = nc[nc>0]
	return len(nc)

def bond(mpos):
    bonds = mpos[0:len(mpos)-1,:] - mpos[1:len(mpos),:] #bond=rij=ri-rj
    return bonds		

def angle(mpos, rnge=None):
    if not rnge:
        rnge=xrange(len(mpos)-2)
    angle = numpy.zeros(len(rnge))
    try:    
        for index,i in enumerate(rnge):
            BA = mpos[i,:] - mpos[i+1,:]
            BC = mpos[i+2,:] - mpos[i+1,:]
            dotBA = BA[0]*BA[0] + BA[1]*BA[1] + BA[2]*BA[2]
            dotBC = BC[0]*BC[0] + BC[1]*BC[1] + BC[2]*BC[2]
            angle[index] = numpy.arccos((BA[0]*BC[0]+BA[1]*BC[1]+BA[2]*BC[2])/(dotBA*dotBC)**.5)
    except IndexError:
        angle=angle[0:-1]
        assert(angle[-1]!=0)
    return angle

def anglem(mpos,i):
    BA = mpos[i-1,:] - mpos[i,:]
    BC = mpos[i+1,:] - mpos[i,:]
    angle = arccos(numpy.dot(BA,BC)/(numpy.dot(BA,BA)**.5*numpy.dot(BC,BC)**.5))
    return angle


def dihedral(mpos, rnge=None):
    if not rnge:
        rnge=xrange(len(mpos)-3)
    newdihed = numpy.zeros(len(rnge))
    for index,i in enumerate(rnge):
        AB = mpos[i,:] - mpos[i+1,:] #rij
        BC = mpos[i+2,:] - mpos[i+1,:] #rkj
        CD = mpos[i+2,:] - mpos[i+3,:] #rkl
        plane1 = numpy.array([AB[1]*BC[2]-AB[2]*BC[1], AB[2]*BC[0]-AB[0]*BC[2], AB[0]*BC[1]-AB[1]*BC[0]])  #cross(AB,BC)
        plane2 = numpy.array([BC[1]*CD[2]-BC[2]*CD[1], BC[2]*CD[0]-BC[0]*CD[2], BC[0]*CD[1]-BC[1]*CD[0]]) #cross(CD,BC)
        AB = (BC[0]**2+BC[1]**2+BC[2]**2)**.5 * AB
        a = AB[0]*plane2[0] + AB[1]*plane2[1] + AB[2]*plane2[2]
        b = plane1[0]*plane2[0] + plane1[1]*plane2[1] + plane1[2]*plane2[2]
        newdihed[index] = numpy.arctan2(a,b)
        if newdihed[index] < 0:
            newdihed[index] += 2*numpy.pi
    return newdihed
 
def cdihedral(mpos, rnge=None):
    if not rnge:
        rnge = numpy.arange(len(mpos)-3)
    else:
        rnge = numpy.array(rnge)
    n = len(rnge)
    newdihed = numpy.zeros(n)
    code = """
    int i;
    double x1, x2, x3, y1, y2, y3, z1, z2, z3;
    double mx, my, mz, nx, ny, nz;
    double magBC, a, b;
    for ( int j = 0; j < n; j++){
        i = RNGE1(j);
        x1 = MPOS2(i,0) - MPOS2(i+1,0); y1 = MPOS2(i,1) - MPOS2(i+1,1); z1 = MPOS2(i,2) - MPOS2(i+1,2);
        x2 = MPOS2(i+2,0) - MPOS2(i+1,0); y2 = MPOS2(i+2,1) - MPOS2(i+1,1); z2 = MPOS2(i+2,2) - MPOS2(i+1,2);
        x3 = MPOS2(i+2,0) - MPOS2(i+3,0); y3 = MPOS2(i+2,1) - MPOS2(i+3,1); z3 = MPOS2(i+2,2) - MPOS2(i+3,2);
        mx = y1*z2 - z1*y2; my = z1*x2 - x1*z2; mz = x1*y2 - y1*x2;
        nx = y2*z3 - z2*y3; ny = z2*x3 - x2*z3; nz = x2*y3 - y2*x3;
        magBC = sqrt(x2*x2 + y2*y2 + z2*z2);
        a = magBC*(x1*nx + y1*ny + z1*nz);
        b = mx*nx + my*ny + mz*nz;
        NEWDIHED1(i) = atan2(a,b);
        if ( NEWDIHED1(i) < 0){
            NEWDIHED1(i) += 2*M_PI;
        }
    }
    """
    info = weave.inline(code, ['mpos', 'n', 'rnge'], headers=['<math.h>', '<stdlib.h>'])
    return newdihed

