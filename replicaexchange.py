import pdb
import numpy

def Q_reweight(replica, k_Qpin, Qpin_new):
        return replica.u0 - k_Qpin*(replica.Q - replica.Qpin)**2 + k_Qpin*(replica.Q - Qpin_new)**2

def tryswap(args, replicas, swapaccepted, swaprejected, start, beta, Q, wiw):
    ''' helper function in tryrepeatedswaps() '''
    for i in xrange(start, args.nreplicas-1, 2):
        if args.Qfile:
            # I believe form is correct
            P = numpy.exp((replicas[wiw[i]].u0 - Q_reweight(replicas[wiw[i+1]], args.k_Qpin, Q[i]))*beta[i]
                        - (Q_reweight(replicas[wiw[i]], args.k_Qpin, Q[i+1]) - replicas[wiw[i+1]].u0)*beta[i+1])
            # incorrect:
            #P = numpy.exp((replicas[i].u0 - Q_reweight(replicas[i], args.k_Qpin, Q[i+1]))*beta[i]
            #            - (Q_reweight(replicas[i+1], args.k_Qpin, Q[i]) - replicas[i+1].u0)*beta[i+1])
        else:
            P = numpy.exp((replicas[wiw[i]].u0-replicas[wiw[i+1]].u0)*(beta[i]-beta[i+1]))
        if numpy.random.random() < P:
            swapaccepted[i] += 1
            wiw[i], wiw[i+1] = wiw[i+1], wiw[i]
        else:
            swaprejected[i] += 1
    return swapaccepted, swaprejected, wiw

def swap(args, replicas, wiw):
    ''' helper function in tryrepeatedswaps() '''
    wiw_start = range(args.nreplicas)
    for i, j in enumerate(wiw):
        if wiw_start[i] != j:
            X = wiw_start.index(j)
            wiw_start[i], wiw_start[X] = wiw_start[X], wiw_start[i]
            replicas[i].coord, replicas[X].coord = replicas[X].coord, replicas[i].coord
            replicas[i].u0, replicas[X].u0 = replicas[X].u0, replicas[i].u0
            replicas[i].r2, replicas[X].r2 = replicas[X].r2, replicas[i].r2
            replicas[i].torsE, replicas[X].torsE = replicas[X].torsE, replicas[i].torsE
            replicas[i].angE, replicas[X].angE = replicas[X].angE, replicas[i].angE
            replicas[i].whoami, replicas[X].whoami = replicas[X].whoami, replicas[i].whoami # whoami keeps track of the individual protein in each Replica 
            if args.surf:
                replicas[i].surfE, replicas[X].surfE = replicas[X].surfE, replicas[i].surfE
            if args.umbrella:
                replicas[i].z_array, replicas[X].z_array = replicas[X].z_array, replicas[i].z_array
            if args.Qfile:
                replicas[i].Q, replicas[X].Q = replicas[X].Q, replicas[i].Q
    assert numpy.all(wiw_start == wiw)     

def tryrepeatedswaps(args, replicas, swapaccepted, swaprejected, protein_location, beta, Q=''):
    if args.nreplicas == 1:
		return swapaccepted, swaprejected, protein_location
    wiw = range(args.nreplicas)
    switch = 0
    for k in range(args.swapsper):
		if switch:
			swapaccepted, swaprejected, wiw = tryswap(args, replicas, swapaccepted, swaprejected, 0, beta, Q, wiw) #type 1 swap
		else:	
			swapaccepted, swaprejected, wiw = tryswap(args, replicas, swapaccepted, swaprejected, 1, beta, Q, wiw) #type 2 swap
		switch = -(switch-1)
    swap(args, replicas, wiw)
    for i in range(args.nreplicas):
		protein_location[replicas[i].whoami].append(i) #extracts which protein went where
    return swapaccepted, swaprejected, protein_location
