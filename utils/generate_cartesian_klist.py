import numpy as np
import wien2k
from wien2k.utils.remove_duplicates import remove_duplicates
from wien2k.utils.reduce_ibz import reduce_ibz
from wien2k.utils.expand_ibz import expand_ibz

import pdb

def generate_cartesian_klist(points, sym_mats=None, reduce_to_ibz=True, outfile=None, verbose=False):

    tot_points = points[0]*points[1]*points[2]
    i_mid_zone = points[0]/2.0
    j_mid_zone = points[1]/2.0
    k_mid_zone = points[2]/2.0

    if verbose == True:
        print 'Generating Kmesh ...'

    klist_buffer = []
    klist = None
    for i, j, k in np.ndindex(points[0], points[1], points[2]):
        klist_buffer.append([i-i_mid_zone, j-j_mid_zone, k-k_mid_zone])

    klist = np.array(klist_buffer)
    del(klist_buffer)


    if verbose == True:
        print '%d points in the full zone ...' % len(klist)

    if reduce_to_ibz == True:
        # Reduce to irreducible Brillouin zone
        if verbose == True:
            print 'Reducing to irreducible Brillouin zone ...'
        klist = reduce_ibz(klist, sym_mats=sym_mats)
        if verbose == True:
            print '%d points in reduced zone ...' % len(klist)

    # Generate the ids and add the ids
    if verbose == True:
        print 'Generating the ids ...'
    ids = np.arange(1, len(klist)+1).astype(np.int32) # Have it move to 32 bit to accommodate huge integers
    ids = ids.reshape((-1,1))
    klist = np.concatenate((ids, klist), axis=1)

    # Add the denominator
    if verbose == True:
        print 'Adding the denominators ...'
    denoms = np.ones((len(klist), 1)) * points[0] # FIXME: Find out what the format is for klist when multiple points allowed
    klist = np.concatenate((klist, denoms), axis=1)

    if reduce_to_ibz == True:
        # Expand out so can find the weights
        if verbose == True:
            print 'calculating the k point weights ...'
        full_zone_ids = expand_ibz(ibz_data=klist, \
          sym_mats=sym_mats, bz_dims=points)[:,0]
        weights = np.bincount(full_zone_ids.astype(np.int32))
        len_weights = len(weights)
        del(full_zone_ids)
        klist = np.concatenate((klist, np.ones((len(klist), 1))), axis=1)
        for i in xrange(len(klist)):
            id = klist[i,0]
            if id >= len_weights:
                klist[i,-1] = 0 # FIXME: Sureley this should not happen? (but it does)       
            else:
                klist[i,-1] = weights[id]        
    else:
        # Weights are 1
        klist = np.concatenate((klist, np.ones((len(klist), 1))), axis=1)

    if outfile is not None:
        # Write to the .klist file    
        if verbose == True:
            print 'Building the .klist file ...'
        kw = wien2k.KlistWriter(outfile)
        kw.total_number_k_points = tot_points
        kw.num_rlv_points = [points[0], points[1], points[2]]
        if verbose == True:
            print 'Writing the .klist file ...'
        kw.data = klist
        kw.write()
        if verbose == True:
            print 'File written to %s' % outfile

    # Return the data
    return klist
