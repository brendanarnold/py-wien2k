
__all__ = ['expand_irreducible_brillouin_zone']

import numpy as np

def expand_irreducible_brillouin_zone(klist_rdr=None, \
        struct_rdr=None, symmetry_matrices=None, \
        ibz_data=None, klist_denominator=None):
    '''
    Expands a data set containing irreducible k points into a full Brillouin
    zone of k points

    The following parameters need to be passed,

    symmetry_matrices           : A list SymmetryMatrix instances
    ibz_data                    : A NxM>=4 array of N k points representing the
                                  irreducible Brillouin zone with each row of form id,kx,ky,kz,...
    klist_denominator           : The single denominator for the klist vectors

    Typically a StructReader instance and a KlistReader instance cover all
    these parameters

    Returns: An OxM array of k points with each row of form id,kx,ky,kz,...

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
    if klist_rdr != None:
        ibz_data = ibz_data or klist_rdr.data
        klist_denominator = klist_denominator or \
                np.unique(klist_rdr.denominators)
    if struct_rdr != None:
        symmetry_matrices = symmetry_matrices or struct_rdr.symmetry_matrices
##     if band != None:
##         ibz_data = ibz_data or band.data
##     if outputkgen_rdr != None:
##         if symmetry_matrices = None:
##             # See note above about having to operate the matrices on the
##             # inverse of the normalised rectangular lattice vectors
##             # FIXME, may not be needed for hexagonal lattices - investigate this
##             rlvs = outputkgen_rdr.reciprocal_lattice_vectors
##             norm_rlvs = rlvs / rlvs.max(axis=0)
##             symmetry_matrices = []
##             for sm in outputkgen_rdr.symmetry_matrices:
##                 symmetry_matrices.append(np.dot(sm,np.linalg.inv(norm_rlvs)))

    if not (len(ibz_data) and klist_denominator and len(symmetry_matrices)):
        raise ValueError('One of the parameters was not specified')

    # Use the Struct vectors to build the full Brollouin zone
    ibz_k_points = ibz_data[:,1:4]
    full_bz_data = np.array([])
    for symmetry_matrix in symmetry_matrices:
        # n.b. Normally use <Operator>.<Vector> to transform but for Numpy to
        # iterate over a list of vectors this has to be reversed and hence the
        # operator transposed
        k_points_buffer = np.dot(ibz_k_points, \
                symmetry_matrix.matrix.transpose()) + symmetry_matrix.tau_offsets
        # Map transforms back into the unit cell
        k_points_buffer = np.mod(k_points_buffer, klist_denominator)
        ## TEST REMOVE
##         k_points_buffer = np.dot(k_points_buffer, np.linalg.inv(rect_latt_vectors).transpose())
        ## END TEST
        # Slot the transformed vectors into a copy of the original data array
        ibz_data_buffer = ibz_data.copy()
        ibz_data_buffer[:,1:4] = k_points_buffer
        # Add buffered data to the end of the list of k points
        ibz_data_buffer_height, ibz_data_buffer_width = ibz_data_buffer.shape
        new_full_bz_data_height = full_bz_data.shape[0] + ibz_data_buffer_height
        new_full_bz_data_width = ibz_data_buffer_width
        full_bz_data = np.resize(full_bz_data, (new_full_bz_data_height, new_full_bz_data_width)) 
        full_bz_data[-ibz_data_buffer_height:,:] = ibz_data_buffer
    # Remove duplicate k points
    kx,ky,kz = full_bz_data[:,1:4].transpose() # A way to sort by kx, then ky,kz,id
    sorted_indices = np.lexsort(keys=(kz,ky,kx))
    full_bz_data = np.take(full_bz_data, sorted_indices, axis=0)
    full_bz_data = full_bz_data.tolist()
    full_bz_data = [kp for i,kp in enumerate(full_bz_data) \
            if i == 0 or kp[1:4] != full_bz_data[i-1][1:4]]
    # Return as floats since rest of data array may be float (i.e. energy from EnergyReader)
    full_bz_data = np.array(full_bz_data, dtype=float)
    return full_bz_data
    
