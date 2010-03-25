
__all__ = ['expand_ibz']

import numpy as np
from wien2k.utils.remove_duplicates import remove_duplicates

def expand_ibz(klist_rdr=None, \
        struct_rdr=None, sym_mats=None, band=None, sort_by=None, \
        ibz_data=None, bz_dims=None, constrain_to_bz=True, cols=[1,2,3]):
    '''
    Expands a data set containing irreducible k points into a full Brillouin
    zone of k points

    The following parameters need to be passed,

    sym_mats   A list SymMat instances
    ibz_data            A NxM>=4 array of N k points representing the
                        irreducible Brillouin zone with each row of form
                        id,kx,ky,kz,[klist denom],...

    Optional parameters are as follows,

    cols                The list of columns in ibz_data that contain the 
                        i, j, k values respectively. (default: [1,2,3])
    constrain_to_bz     Bool. (default: True) Symmetry matrices
                        may map points outside of Brillouin zone
                        (i.e. extended zone scheme), if True
                        then the function will map the points
                        back inside the Brillouin zone
    bz_dims             The dimensions of the Brillouin Zone used for
                        constraining the zone, if not specified this
                        will be read from the KlistReader or the
                        ibz_data, if not included will raise a ValueError
    sort_by             A list of columns to sort by in order (default: None)

    The above can be covered with the following objects,

    struct_rdr          A StructReader instance will cover the symmetry matrices
    klist_rdr           A KlistReader instance will cover the ibz_data and the bz_dims arguments
    band                A Band instance will cover the ibz_data

    Returns: An OxM array of k points with each row of form id,kx,ky,kz,...

    Example:

    Manually specifying each value and 'expanding' using the identity
    matrix (note that the list is no longer in order)

    >>> import wien2k
    >>> ibz_data = np.array([ 
    ...   [1,0,0,0,1.1],  
    ...   [2,0,0,1,1.2],  
    ...   [3,0,1,0,1.0],  
    ...   [4,0,1,1,1.1],  
    ...   [5,1,0,0,1.2],  
    ...   [6,1,0,1,1.2],  
    ...   [7,1,1,0,1.2], 
    ...   [8,1,1,1,1.0]])
    >>> sms = [wien2k.SymMat(matrix=np.identity(3))]
    >>> bz_dims = [1,1,1]
    >>> expand_ibz(sym_mats=sms, \
        bz_dims=bz_dims, ibz_data=ibz_data, sort_by=[0])
    array([[ 1. ,  0. ,  0. ,  0. ,  1.1],
           [ 2. ,  0. ,  0. ,  1. ,  1.2],
           [ 3. ,  0. ,  1. ,  0. ,  1. ],
           [ 4. ,  0. ,  1. ,  1. ,  1.1],
           [ 5. ,  1. ,  0. ,  0. ,  1.2],
           [ 6. ,  1. ,  0. ,  1. ,  1.2],
           [ 7. ,  1. ,  1. ,  0. ,  1.2],
           [ 8. ,  1. ,  1. ,  1. ,  1. ]])

    '''

    # == CALCULATING THE RECTANGULAR LATTICE ==
    #
    # Let k be a k point vector
    # Let M be the .struct symmetry matrices
    # Let N be the .outputkgen symmetry matrices
    # Let R be the rectangular lattice vectors as read from .outputkgen after normalisation
    #
    # The k point in rectagular co-ordinates can be found from either:
    #     M.k     or
    #     N.R^(-1).k
    #

    # Assign the parameters depending on how the function was called
    if klist_rdr is not None:
        if ibz_data is None:
            ibz_data = klist_rdr.data
        if bz_dims is None:
            klist_denom = np.unique(klist_rdr.denominators)[0]
            bz_dims = 3*[klist_denominator]
    if struct_rdr is not None:
        if sym_mats is not None:
            raise ValueError('Multiple sets of symmetry matrices specified - struct reader and a set of symmetry matrices')
        else:
            sym_mats = sym_mats or struct_rdr.sym_mats
    if band is not None:
        # Don't check for ibz_data as Klist may have set it
        ibz_data = band.data

##     if band != None:
##         ibz_data = ibz_data or band.data
##     if outputkgen_rdr != None:
##         if sym_mats = None:
##             # See note above about having to operate the matrices on the
##             # inverse of the normalised rectangular lattice vectors
##             # FIXME, may not be needed for hexagonal lattices - investigate this
##             rlvs = outputkgen_rdr.rlvs
##             norm_rlvs = rlvs / rlvs.max(axis=0)
##             sym_mats = []
##             for sm in outputkgen_rdr.sym_mats:
##                 sym_mats.append(np.dot(sm,np.linalg.inv(norm_rlvs)))

    # Complain a bit if necessary
    if (ibz_data is None) and (sym_mats is None):
        raise ValueError('One of the parameters was not specified')
    if (constrain_to_bz == True) and \
        ((bz_dims is None) or (len(bz_dims) != 3)):
        raise ValueError('Need to specify the dimensions of the Brillouin Zone if points are to mapped back inside the Brillouin zone')
    if (len(sym_mats) == 0):
        raise ValueError('No symmetry matrices in list')

    # Use the Struct vectors to build the full Brollouin zone
    full_bz = None
    for symmetry_matrix in sym_mats:
        kpoints_buffer = symmetry_matrix.map(ibz_data, cols=cols)
        # Map transforms back into the unit cell if required
        if constrain_to_bz == True:
            i_col, j_col, k_col = kpoints_buffer[:,cols].transpose()
            while i_col.max() > bz_dims[0]:
                i_col[i_col > bz_dims[0]] -= bz_dims[0]
            while j_col.max() > bz_dims[1]:
                j_col[j_col > bz_dims[1]] -= bz_dims[1]
            while k_col.max() > bz_dims[0]:
                k_col[k_col > bz_dims[2]] -= bz_dims[2]
            kpoints_buffer[:,cols[0]] = i_col
            kpoints_buffer[:,cols[1]] = j_col
            kpoints_buffer[:,cols[2]] = k_col
        # Add buffered data to the end of the list of k points
        if full_bz is None:
            full_bz = kpoints_buffer
        else:
            full_bz = np.concatenate((full_bz, kpoints_buffer))
    # Remove duplicate k points
    full_bz = remove_duplicates(full_bz)

##     kx,ky,kz = full_bz[:,1:4].transpose() # A way to sort by kx, then ky,kz,id
##     sorted_indices = np.lexsort(keys=(kz,ky,kx))
##     full_bz = np.take(full_bz, sorted_indices, axis=0)
##     full_bz = full_bz.tolist()
##     full_bz = [kp for i,kp in enumerate(full_bz) \
##             if i == 0 or kp[1:4] != full_bz[i-1][1:4]]

    # Sort the results if required
    if sort_by is not None:
        sorted_inds = np.lexsort([full_bz[:,c] for c in sort_by])
        full_bz = full_bz[sorted_inds,:]

    # Returns a copy
    return full_bz.copy()



if __name__ == '__main__':
    import doctest
    doctest.testmod()
    
