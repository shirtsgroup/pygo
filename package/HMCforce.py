import numpy
import energyfunc
import pdb
from scipy import weave
#==========================================
# FORCE CALCULATION METHODS
#==========================================

def cdihedforces(torsparam, bonds, dsq, d, numbeads):
    forces = numpy.zeros((numbeads,3))
    code = """
    double x1, x2, x3, y1, y2, y3, z1, z2, z3;
    double mx, my, mz, m2, nx, ny, nz, n2;
    double magBC, a, b, dihed, dV;
    double Fix, Fiy, Fiz, Fjx, Fjy, Fjz, Fkx, Fky, Fkz, Flx, Fly, Flz;
    for (int i = 0; i < numbeads-3; i++){
        x1 = BONDS2(i,0); y1 = BONDS2(i,1); z1 = BONDS2(i,2);
        x2 = -BONDS2(i+1,0); y2 = -BONDS2(i+1,1); z2 = -BONDS2(i+1,2);
        x3 = BONDS2(i+2,0); y3 = BONDS2(i+2,1); z3 = BONDS2(i+2,2);
        mx = y1*z2 - z1*y2; my = z1*x2 - x1*z2; mz = x1*y2 - y1*x2;
        nx = y2*z3 - z2*y3; ny = z2*x3 - x2*z3; nz = x2*y3 - y2*x3;
        m2 = mx*mx + my*my + mz*mz;
        n2 = nx*nx + ny*ny + nz*nz;
        magBC = sqrt(x2*x2 + y2*y2 + z2*z2);
        a = magBC*(x1*nx + y1*ny + z1*nz);
        b = mx*nx + my*ny + mz*nz;
        dihed = atan2(a,b);
        if ( dihed < 0){
            dihed += 2*M_PI;
        }
        dV = -TORSPARAM2(4*i,0)*TORSPARAM2(4*i,1)*sin(TORSPARAM2(4*i,1)*dihed-TORSPARAM2(4*i,2)) - TORSPARAM2(4*i+1,0)*TORSPARAM2(4*i+1,1)*sin(TORSPARAM2(4*i+1,1)*dihed-TORSPARAM2(4*i+1,2)) - TORSPARAM2(4*i+2,0)*TORSPARAM2(4*i+2,1)*sin(TORSPARAM2(4*i+2,1)*dihed-TORSPARAM2(4*i+2,2)) - TORSPARAM2(4*i+3,0)*TORSPARAM2(4*i+3,1)*sin(TORSPARAM2(4*i+3,1)*dihed-TORSPARAM2(4*i+3,2));
        Fix = -dV * D1(i+1) * mx / m2; Fiy = -dV * D1(i+1) * my / m2; Fiz = -dV * D1(i+1) * mz / m2;
        Flx = dV * D1(i+1) * nx / n2; Fly = dV * D1(i+1) * ny / n2; Flz = dV * D1(i+1) * nz / n2;
        a = (x1*x2 + y1*y2 + z1*z2) / DSQ1(i+1);
        b = (x2*x3 + y2*y3 + z2*z3) / DSQ1(i+1);
        Fjx = (a-1)*Fix - b*Flx; Fjy = (a-1)*Fiy - b*Fly; Fjz = (a-1)*Fiz - b*Flz;
        Fkx = -a*Fix + (b-1)*Flx; Fky = -a*Fiy + (b-1)*Fly; Fkz = -a*Fiz + (b-1)*Flz;
        FORCES2(i,0) += Fix; FORCES2(i,1) += Fiy; FORCES2(i,2) += Fiz;
        FORCES2(i+1,0) += Fjx; FORCES2(i+1,1) += Fjy; FORCES2(i+1,2) += Fjz;
        FORCES2(i+2,0) += Fkx; FORCES2(i+2,1) += Fky; FORCES2(i+2,2) += Fkz;
        FORCES2(i+3,0) += Flx; FORCES2(i+3,1) += Fly; FORCES2(i+3,2) += Flz;
    }
    """
    info = weave.inline(code, ['forces', 'torsparam', 'numbeads', 'bonds', 'dsq', 'd'], headers=['<math.h>', '<stdlib.h>'])
    return forces*4.184

def cangleforces(mpos, angleparam, bonds, d, numbeads):
    rki = mpos[0:-2,:] - mpos[2:len(mpos),:]
    forces = numpy.zeros((numbeads,3))
    code = """
    double xba, yba, zba, xbc, ybc, zbc;
    double angle, dV;
    double fpx, fpy, fpz, ba2, xy, xz, yz;
    double Fix, Fiy, Fiz, Fkx, Fky, Fkz;
    for ( int i = 0; i < numbeads - 2; i++){
        xba = BONDS2(i,0); yba = BONDS2(i,1); zba = BONDS2(i,2);
        xbc = -BONDS2(i+1,0); ybc = -BONDS2(i+1,1); zbc = -BONDS2(i+1,2);
        angle = acos((xba*xbc + yba*ybc + zba*zbc) / (D1(i)*D1(i+1)));
        //dV = 2 * ANGLEPARAM2(i,0) * (cos(angle) - cos(ANGLEPARAM2(i,1))) / (sin(ANGLEPARAM2(i,1))*sin(ANGLEPARAM2(i,1)));
        dV = 2*ANGLEPARAM2(i,0)*(angle-ANGLEPARAM2(i,1))/-sin(angle);
        fpx = dV / (D1(i)*D1(i+1)) * RKI2(i,0); fpy = dV / (D1(i)*D1(i+1)) * RKI2(i,1); fpz = dV / (D1(i)*D1(i+1)) * RKI2(i,2);
        ba2 = xba*xba + yba*yba + zba*zba;
        xy = xba*yba/ba2; xz = xba*zba/ba2; yz = yba*zba/ba2;
        Fix = fpx*(1 - xba*xba/ba2) - fpy*xy - fpz*xz;
        Fiy = -fpx*xy + fpy*(1 - yba*yba/ba2) - fpz*yz;
        Fiz = -fpx*xz - fpy*yz + fpz*(1 - zba*zba/ba2);
        ba2 = xbc*xbc + ybc*ybc + zbc*zbc;
        xy = xbc*ybc/ba2; xz = xbc*zbc/ba2; yz = ybc*zbc/ba2;
        Fkx = -fpx*(1 - xbc*xbc/ba2) + fpy*xy + fpz*xz;
        Fky = fpx*xy - fpy*(1 - ybc*ybc/ba2) + fpz*yz;
        Fkz = fpx*xz + fpy*yz - fpz*(1 - zbc*zbc/ba2);
        FORCES2(i,0) += Fix; FORCES2(i,1) += Fiy; FORCES2(i,2) += Fiz;
        FORCES2(i+1,0) += -Fix-Fkx; FORCES2(i+1,1) += -Fiy-Fky; FORCES2(i+1,2) += -Fiz-Fkz;
        FORCES2(i+2,0) += Fkx; FORCES2(i+2,1) += Fky; FORCES2(i+2,2) += Fkz;
    }
    """
    info = weave.inline(code, ['forces', 'angleparam', 'numbeads', 'bonds', 'rki', 'd'], headers=['<math.h>', '<stdlib.h>'])
    return forces*4.184
                
def bondedforces(mpos, torsparam, angleparam, bonds, d2, d, numbeads):
	"""Returns the bonded forces (angle + torsion) from a given configuration"""
	forces = numpy.zeros((numbeads,3))
	arccos = numpy.arccos
        sin = numpy.sin
        cos = numpy.cos
        array = numpy.array
        
        # dihedral force calculation
	for i in range(numbeads-3):
		# find dihedral angle
		rij = bonds[i]
		rkj = -bonds[i+1]
		rkl = bonds[i+2]
		m = array([rij[1]*rkj[2]-rij[2]*rkj[1], rij[2]*rkj[0]-rij[0]*rkj[2], rij[0]*rkj[1]-rij[1]*rkj[0]])
		n = array([rkj[1]*rkl[2]-rkj[2]*rkl[1], rkj[2]*rkl[0]-rkj[0]*rkl[2], rkj[0]*rkl[1]-rkj[1]*rkl[0]])
		m2 = m[0]*m[0] + m[1]*m[1] + m[2]*m[2]
		n2 = n[0]*n[0] + n[1]*n[1] + n[2]*n[2]
		dihed = arccos((m[0]*n[0]+m[1]*n[1]+m[2]*n[2]) / (m2*n2)**.5)
		if ((m[0]*rkl[0]+m[1]*rkl[1]+m[2]*rkl[2]) < 0):
			dihed = -abs(dihed) + 2 * numpy.pi
		else:
			dihed = abs(dihed)
		# calculate gradient of dihedral potential and resulting forces on i j k l
		dV = -torsparam[4*i:4*i+4,0] * torsparam[4*i:4*i+4,1] * sin(torsparam[4*i:4*i+4,1] * dihed - torsparam[4*i:4*i+4,2])
		dV = dV[0]+dV[1]+dV[2]+dV[3]
		Fi = -dV * d[i+1]
		Fl = -Fi * n / n2
		Fi = Fi * m / m2
		Fj = (rij[0]*rkj[0]+rij[1]*rkj[1]+rij[2]*rkj[2]) / d2[i+1] * Fi - (rkl[0]*rkj[0]+rkl[1]*rkj[1]+rkl[2]*rkj[2]) / d2[i+1] * Fl
		Fk = -Fj - Fl
		Fj = -Fi + Fj
		# add forces to total force
		forces[i,:] += Fi
		forces[i+1,:] += Fj
		forces[i+2,:] += Fk	
		forces[i+3,:] += Fl

	#angle force calculation
		# find angle
		rki = mpos[i,:] - mpos[i+2,:]
		rji = bonds[i]
		rjk = -bonds[i+1]
		dotBABC = rji[0]*rjk[0] + rji[1]*rjk[1] + rji[2]*rjk[2]
		dotBA = rji[0]*rji[0] + rji[1]*rji[1] + rji[2]*rji[2]
		dotBC = rjk[0]*rjk[0] + rjk[1]*rjk[1] + rjk[2]*rjk[2]
        	angle = numpy.arccos((rji[0]*rjk[0] + rji[1]*rjk[1] + rji[2]*rjk[2]) / (d[i]*d[i+1])) #in radians
		# calculate gradient of angle potential and resulting forces on i j k
		dV = 2 * angleparam[i,0] * (cos(angle) - cos(angleparam[i,1])) / sin(angleparam[i,1])**2
		fprime = dV / (d[i]*d[i+1]) * rki
                Fi = array(projdot(rji,fprime))
                Fk = -array(projdot(rjk,fprime))
		Fj = -Fi - Fk
		# add forces to total force
		forces[i,:] += Fi
		forces[i+1,:] += Fj
		forces[i+2,:] += Fk
	# angle force calculations for three end beads not included in loop
	i = numbeads - 3
	rki = mpos[i,:] - mpos[i+2,:]
	rji = bonds[i]
	rjk = -bonds[i+1]
	dotBABC = rji[0]*rjk[0] + rji[1]*rjk[1] + rji[2]*rjk[2]
	dotBA = rji[0]*rji[0] + rji[1]*rji[1] + rji[2]*rji[2]
	dotBC = rjk[0]*rjk[0] + rjk[1]*rjk[1] + rjk[2]*rjk[2]
	angle = arccos((rji[0]*rjk[0] + rji[1]*rjk[1] + rji[2]*rjk[2]) / (d[i]*d[i+1])) #in radians
	dV = 2 * angleparam[i,0] * (cos(angle) - cos(angleparam[i,1])) / sin(angleparam[i,1])**2
	fprime = dV / (d[i]*d[i+1]) * rki
        Fi = array(projdot(rji,fprime))
        Fk = -array(projdot(rjk,fprime))
        Fj = -Fi - Fk
        # add forces to total force
        forces[i,:] += Fi
        forces[i+1,:] += Fj
        forces[i+2,:] += Fk
	return forces *4.184 # convert forces to kJ/mol/K

def projdot(p,d):
    p2 = p[0]*p[0] + p[1]*p[1] + p[2]*p[2]
    x2 = 1 - p[0]**2/p2
    y2 = 1 - p[1]**2/p2
    z2 = 1 - p[2]**2/p2
    xy = p[0] * p[1] / p2
    xz = p[0] * p[2] / p2
    yz = p[1] * p[2] / p2
    xcomp = d[0]*x2 - d[1]*xy - d[2]*xz
    ycomp = -d[0]*xy + d[1]*y2 - d[2]*yz
    zcomp = -d[0]*xz - d[1]*yz + d[2]*z2
    return [xcomp,ycomp,zcomp]

def nonbondedforces(mpos, numint, numbeads, natparam, nonnatparam, nnepsil):
	"""Returns the nonbonded forces of a given configuration"""
	forces = numpy.zeros((numbeads,3))
	# get distances, square distances, and magnitude distances for all interactions 
	rvec = energyfunc.cgetforcer(mpos, numint, numbeads) # excludes 12 and 13 neightbors
	r2 = numpy.sum(rvec**2, axis=1)
	rr = r2**.5
	# calculates potential energy gradient for native and nonnative interactions
	ndV = natparam[:,0] * natparam[:,2] * natparam[:,2] / r2
	ndV6 = ndV * ndV * ndV
	ndV = natparam[:,1] * (-156*ndV6*ndV6/rr + 180*ndV6*ndV*ndV/rr - 24*ndV6/rr)
	nndV = nonnatparam[:,0] * nonnatparam[:,1] * nonnatparam[:,1] / r2
	nndV = nndV * nndV * nndV
	nndV = -12 * nnepsil * nndV * nndV / rr
	# add forces to total force
	k = 0
	F = -(ndV+nndV) / rr
	F = numpy.transpose(F) * numpy.transpose(rvec) # now 3xN instead of Nx3
        F = F.transpose()
        #for i in range(numbeads):
		#for j in range(i+3,numbeads):
			#forces[i,:] += F[:,k]
			#forces[j,:] += -F[:,k]
			#k += 1
                        
        code = """
        int k = 0;
        for ( int i = 0; i < numbeads; i++){
            for ( int j = i+3; j < numbeads; j++){
                FORCES2(i,0) += F2(k,0);
                FORCES2(i,1) += F2(k,1);
                FORCES2(i,2) += F2(k,2);
                FORCES2(j,0) += -F2(k,0);
                FORCES2(j,1) += -F2(k,1);
                FORCES2(j,2) += -F2(k,2);
                k++;
            }
        }
        """
        info = weave.inline(code, ['forces', 'F', 'numbeads'], headers=['<math.h>', '<stdlib.h>'])
	return forces * 4.184 # converts force to kJ/mol/K

def cnonbondedforces(mpos, numint, numbeads, natparam, nonnatparam, nnepsil):
	"""Returns the nonbonded forces of a given configuration"""
	forces = numpy.zeros((numbeads,3))
	# get distances, square distances, and magnitude distances for all interactions 
	#rvec = energyfunc.cgetforcer(mpos, numint, numbeads) # excludes 12 and 13 neightbors
	code = """
        int k = 0;
        double r2, r, ndV, ndV6, F,x,y,z;
        for ( int i = 0; i < numbeads; i++){
            for ( int j = i+3; j < numbeads; j++){
                x = MPOS2(i,0) - MPOS2(j,0);
                y = MPOS2(i,1) - MPOS2(j,1);
                z = MPOS2(i,2) - MPOS2(j,2);
                r2 = x*x + y*y + z*z;
                r = sqrt(r2);
                if (NATPARAM2(k,0) == 1){
                    ndV = NATPARAM2(k,2)*NATPARAM2(k,2)/r2;
                    ndV6 = ndV*ndV*ndV;
                    ndV = NATPARAM2(k,1)*(-156*ndV6*ndV6/r + 180*ndV6*ndV*ndV/r - 24*ndV6/r);
                }
                else{
                    ndV = NONNATPARAM2(k,1)*NONNATPARAM2(k,1)/r2;
                    ndV = ndV*ndV*ndV;
                    ndV = -12*nnepsil*ndV*ndV/r;
                }
                F = -ndV/r;
                FORCES2(i,0) += F*x;
                FORCES2(i,1) += F*y;
                FORCES2(i,2) += F*z;
                FORCES2(j,0) += -F*x;
                FORCES2(j,1) += -F*y;
                FORCES2(j,2) += -F*z;
                k++;
            }
        }
        """
        info = weave.inline(code, ['forces', 'mpos', 'numbeads', 'natparam', 'nonnatparam', 'nnepsil'], headers=['<math.h>', '<stdlib.h>'])

	return forces * 4.184 # converts force to kJ/mol/K

def getsurfforce(prot_coord, surf_coord, numint, numbeads, param):
    ep = param[0]
    sig = param[1]
    scale = param[2]
    rvec = numpy.zeros((numint,3))
    for i in range(len(surf_coord)):
        rvec[i*numbeads:i*numbeads+numbeads] = surf_coord[i,:] - prot_coord
    r2 = numpy.sum(rvec**2,axis = 1)
    r = r2**.5
    ndV = sig*sig/r2
    ndV6 = ndV*ndV*ndV
    ndV12 = ndV6*ndV6
    ndV = -12*ep*ndV12/r+scale[:,0]*ep*(-144*ndV12/r+180*ndV6*ndV*ndV/r - 24*ndV6/r)
    F = -ndV/r
    F = numpy.transpose(F)*numpy.transpose(rvec)
    F = F.transpose()
    forces = numpy.zeros((numbeads,3))
    for i in range(len(surf_coord)):
   	forces += -F[i*numbeads:i*numbeads+numbeads,:]
    return forces*4.184

def cgetsurfforce(prot_coord, surf_coord, numint, numbeads, param, _):
    forces = numpy.zeros((numbeads,3))
    ep_in = param[0]
    sig_in = param[1]
    scale = param[2][:,0]
    code = """
    double ep = ep_in;
    double sig = sig_in;
    double x, y, z, r2, r, dV, dV6, F;
    for ( int i = 0; i < numint; i++){
        x = SURF_COORD2(i/numbeads,0) - PROT_COORD2(i % numbeads, 0);
        y = SURF_COORD2(i/numbeads,1) - PROT_COORD2(i % numbeads, 1);
        z = 0 - PROT_COORD2(i % numbeads, 2);
        r2 = x*x + y*y + z*z;
	if(r2<400){
            r = sqrt(r2);
            dV = sig*sig/r2;
            dV6 = dV*dV*dV;
            dV = SCALE1(i)*ep*(-(12/SCALE1(i)+144)*dV6*dV6/r+180*dV6*dV*dV/r-24*dV6/r);
            F = -dV/r;
            x = F*x; y = F*y; z = F*z;
            FORCES2(i % numbeads, 0) -= x;
            FORCES2(i % numbeads, 1) -= y;
            FORCES2(i % numbeads, 2) -= z;
	}
    }
    """
    info = weave.inline(code, ['forces', 'prot_coord', 'surf_coord', 'numbeads', 'numint', 'ep_in', 'sig_in', 'scale'], headers=['<math.h>', '<stdlib.h>'])
    return forces*4.184

def getsurfforce_old(prot_coord, surf_coord, numint, numbeads, param):
    rvec = numpy.zeros((numint,3))
    for i in range(len(surf_coord)):
        rvec[i*numbeads:i*numbeads+numbeads] = surf_coord[i,:] - prot_coord
    r2 = numpy.sum(rvec**2,axis = 1)
    r = r2**.5
    ndV = param[:,1]*param[:,1]/r2
    ndV6 = ndV*ndV*ndV
    #ndV = param[:,0]*(-156*ndV6*ndV6/r+180*ndV6*ndV*ndV/r - 24*ndV6/r)
    ndV = param[:,0]*(-12*ndV6*ndV6/r + 12*ndV6/r)
    F = -ndV/r
    F = numpy.transpose(F)*numpy.transpose(rvec)
    F = F.transpose()
    forces = numpy.zeros((numbeads,3))
    for i in range(len(surf_coord)):
   	forces += -F[i*numbeads:i*numbeads+numbeads,:]
    return forces*4.184

def cgetsurfforce_old(prot_coord, surf_coord, numint, numbeads, param, scale):
    forces = numpy.zeros((numbeads,3))
    code = """
    double x, y, z, r2, r, dV, F;
    for ( int i = 0; i < numint; i++){
        x = SURF_COORD2(i/numbeads,0) - PROT_COORD2(i % numbeads, 0);
        y = SURF_COORD2(i/numbeads,1) - PROT_COORD2(i % numbeads, 1);
        z = 0 - PROT_COORD2(i % numbeads, 2);
        r2 = x*x + y*y + z*z;
	if(r2<400){
            r = sqrt(r2);
            dV = PARAM2(i,1)*PARAM2(i,1)/r2;
            dV = dV*dV*dV;
            dV = PARAM2(i,0)*(-12*dV*dV/r + 12*scale*dV/r);
            F = -dV/r;
            x = F*x; y = F*y; z = F*z;
            FORCES2(i % numbeads, 0) -= x;
            FORCES2(i % numbeads, 1) -= y;
            FORCES2(i % numbeads, 2) -= z;
	}
    }
    """
    info = weave.inline(code, ['forces', 'prot_coord', 'surf_coord', 'numbeads', 'numint', 'param','scale'], headers=['<math.h>', '<stdlib.h>'])
    return forces*4.184

#==========================================
# BOND CONSTRAINT METHODS
#==========================================

def shake(bonds, v_half, h, m, d2, maxloop, numbeads, tol):
    """Performs SHAKE algorithm to constrain positions"""
    #pdb.set_trace()
    loops = 0
    conv = numbeads -1
    while conv != 0 and loops < maxloop:
        conv = numbeads -1
        for i in range(numbeads-1):
            s = bonds[i] + h * (v_half[i,:]-v_half[i+1,:])
            diff = s[0]*s[0] + s[1]*s[1] + s[2]*s[2] - d2[i]
            if numpy.abs(diff) < tol:
                conv -= 1
            else:
                g=diff/(2*h*(s[0]*bonds[i,0]+s[1]*bonds[i,1]+s[2]*bonds[i,2])*(1/m[i]+1/m[i+1]))
                v_half[i,:]-=g/m[i]*bonds[i]
                v_half[i+1,:]+=g/m[i+1]*bonds[i]
        loops += 1
    conv = True
    if loops == maxloop:
        conv = False
    return v_half, conv

def rattle(bonds, vel, m, d2, maxloop, numbeads, tol):
    """Performs RATTLE algorithm to constrain velocities"""
    conv = numbeads - 1 #numpy.ones(numbeads-1)
    loops = 0
    while conv != 0 and loops < maxloop:
        conv = numbeads - 1
        for i in range(numbeads-1):
            vij = vel[i,:] - vel[i+1,:]
            diff = bonds[i,0]*vij[0] + bonds[i,1]*vij[1] + bonds[i,2]*vij[2]
            if numpy.abs(diff) < tol:
                conv -= 1
            else:
                k = diff / (d2[i] * (1/m[i]+1/m[i+1]))
                vel[i] -= k / m[i] * bonds[i,:]
                vel[i+1] += k / m[i+1] * bonds[i,:]
        loops += 1
    conv = True
    if loops == maxloop:
        conv = False
    return vel, conv

def cshake(bonds, v_half, h, m, dsq, maxloop, numbeads, tol):
    """Performs SHAKE algorithm to constrain positions"""
    loops = numpy.array([0])
    code="""
    int conv = numbeads - 1;
    double x, y, z, diff, g;
    while ( conv != 0 && LOOPS1(0) < maxloop) {
        conv = numbeads - 1;
        for ( int i = 0; i < numbeads-1; i++){
            x = BONDS2(i,0) + h * (V_HALF2(i,0)-V_HALF2(i+1,0));
            y = BONDS2(i,1) + h * (V_HALF2(i,1)-V_HALF2(i+1,1));
            z = BONDS2(i,2) + h * (V_HALF2(i,2)-V_HALF2(i+1,2));
            diff = x * x + y * y + z * z - DSQ1(i);
            if (fabs(diff) < tol){
                conv -= 1;
            }
            else{
                g = diff / (2.0*h*(x*BONDS2(i,0)+y*BONDS2(i,1)+z*BONDS2(i,2))*(1.0/M1(i)+1.0/M1(i+1)));
                V_HALF2(i,0) -= g/M1(i)*BONDS2(i,0);
                V_HALF2(i,1) -= g/M1(i)*BONDS2(i,1);
                V_HALF2(i,2) -= g/M1(i)*BONDS2(i,2);
                V_HALF2(i+1,0) += g/M1(i+1)*BONDS2(i,0);
                V_HALF2(i+1,1) += g/M1(i+1)*BONDS2(i,1);
                V_HALF2(i+1,2) += g/M1(i+1)*BONDS2(i,2);
            }
        }
        LOOPS1(0) += 1;
    }
    """
    info = weave.inline(code, ['bonds', 'v_half', 'h', 'm', 'dsq', 'maxloop', 'numbeads', 'tol','loops'], headers=['<math.h>', '<stdlib.h>'])
    #print "cshake iterations " + str(loops[0])
    conv = True
    if loops[0] == maxloop:
        conv = False
    return v_half, conv

def crattle(bonds, vel, m, dsq, maxloop, numbeads, tol):
    """Performs RATTLE algorithm to constrain velocities"""
    loops = numpy.array([0])
    code = """
    int conv = numbeads - 1;
    double diff, k, x, y, z;
    while ( conv != 0 && LOOPS1(0) < maxloop){
        conv = numbeads - 1;
        for ( int i = 0; i < numbeads-1; i++){
            x = VEL2(i,0) - VEL2(i+1,0);
            y = VEL2(i,1) - VEL2(i+1,1);
            z = VEL2(i,2) - VEL2(i+1,2); 
            diff = BONDS2(i,0)*x + BONDS2(i,1)*y + BONDS2(i,2)*z;
            if (fabs(diff) < tol){
                conv -= 1;
            }
            else{
                k = diff / (DSQ1(i) * (1.0/M1(i)+1.0/M1(i+1)));
                VEL2(i,0) -= k/M1(i)*BONDS2(i,0);
                VEL2(i,1) -= k/M1(i)*BONDS2(i,1);
                VEL2(i,2) -= k/M1(i)*BONDS2(i,2);
                VEL2(i+1,0) += k/M1(i+1)*BONDS2(i,0);
                VEL2(i+1,1) += k/M1(i+1)*BONDS2(i,1);
                VEL2(i+1,2) += k/M1(i+1)*BONDS2(i,2);
            }
        }
        LOOPS1(0)++;
    }
    """
    info = weave.inline(code, ['bonds', 'vel', 'm', 'dsq', 'maxloop', 'numbeads', 'tol', 'loops'], headers=['<math.h>', '<stdlib.h>'])
    conv = True
    if loops[0] == maxloop:
        conv = False
    return vel, conv
    
    
    
