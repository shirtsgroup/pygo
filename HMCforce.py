import numpy
import energyfunc
import pdb

#==========================================
# FORCE CALCULATION METHODS
#==========================================

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
	return forces * 4.184 # convert forces to kJ/mol/K

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
	rvec = energyfunc.getforcer(mpos, numint, numbeads) # excludes 12 and 13 neightbors
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
	for i in range(numbeads):
		for j in range(i+3,numbeads):
			forces[i,:] += F[:,k]
			forces[j,:] += -F[:,k]
			k += 1
	return forces * 4.184 # converts force to kJ/mol/K

#==========================================
# BOND CONSTRAINT METHODS
#==========================================

def shake(bonds, v_half, h, m, d2, maxloop, numbeads, tol):
    """Performs SHAKE algorithm to constrain positions"""
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