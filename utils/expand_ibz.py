__all__ = ['expand_ibz']

import sys
import numpy as np
from wien2k.utils.remove_duplicates import remove_duplicates

def expand_ibz(klist_rdr=None, \
        outputkgen_rdr=None, sym_mats=None, band=None, \
        sort_by=None, ibz_data=None, bz_dims=None, bz_centre=None, \
        constrain_to_bz=True, cols=[1,2,3], rlvs=None):
    '''
    Expands a data set containing irreducible k points into a full Brillouin
    zone of k points

    The following parameters need to be passed,

    sym_mats            A list SymMat instances
    rlvs                A 3x3 Numpy array of reciprocal lattice vectors,
                        default is the identity matrix
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
                        will be read from the KlistReader or the if not 
                        included will default to [1,1,1]
    bz_centre           The co-ordinates of the Brillouin zone centre in 
                        terms of the bz_dims (default: Half the bz_dims values)
    sort_by             A list of columns to sort by in order (default: [0])

    The above can be covered with the following objects,

    outputkgen_rdr      An OutputkgenReader instance will cover the
                        symmetry matrices and the rlvs arguments
    klist_rdr           A KlistReader instance will cover the ibz_data
                        and the bz_dims
    band                A Band instance will cover the ibz_data

    Returns: An OxM array of k points with each row of form id,kx,ky,kz,...

    EXAMPLES:

    Expanding a .klist

    >>> import wien2k
    >>> klist_rdr = wien2k.KlistReader(TiC_klist_filename)
    >>> outputkgen_rdr = wien2k.OutputkgenReader(TiC_outputkgen_filename)
    >>> klist_rdr.data.shape     # There is 47 values in the IBZ
    (47, 6)
    >>> klist_rdr.bz_shape       # Expanded we get 10x10x10 = 1000 values
    (10, 10, 10)
    >>> full_zone_data = expand_ibz(outputkgen_rdr=outputkgen_rdr, klist_rdr=klist_rdr)
    >>> full_zone_data.shape     # We have 1000 values
    (1000, 6)

    Expanding a Band read from a .energy file

    >>> energy_rdr = wien2k.EnergyReader(TiC_energy_filename)
    >>> band = energy_rdr.bands[6]  # A band that crosses the Fermi energy
    >>> band.data.shape
    (47, 5)
    >>> full_zone_data = expand_ibz(outputkgen_rdr=outputkgen_rdr, band=band)
    >>> full_zone_data.shape
    (1000, 5)

    '''

    # == CALCULATING THE RECTANGULAR LATTICE ==
    #
    # Let k be a k point vector
    # Let M be the .struct symmetry matrices
    # Let N be the .outputkgen symmetry matrices
    # Let R be the rectangular lattice vectors as read from .outputkgen after normalisation
    #
    # The k point in rectagular co-ordinates can be found from:
    #     N.R^(-1).k
    #
    # Don't know how to do this with the struct reader

    # Assign the parameters depending on how the function was called
    if band is not None:
        ibz_data = band.data
    if klist_rdr is not None:
        if ibz_data is None:
            ibz_data = klist_rdr.data
        if bz_dims is None:
            bz_dims = klist_rdr.bz_shape 
    if outputkgen_rdr is not None:
        if sym_mats is None:
            sym_mats = outputkgen_rdr.sym_mats
        if rlvs is None:
            # Must normalise the reciprocal lattice vectors from 
            # outputkgen - don't know why ...
            rlvs = outputkgen_rdr.rlvs / outputkgen_rdr.rlvs.max()
    # Set the default values if needed
    if bz_dims is None:
        bz_dims = [1.,1.,1.]
    if rlvs is None:
        rlvs = np.identity(3)
    if sort_by is None:
        sort_by = [0]
    if (bz_centre is None) and (bz_dims is not None):
        bz_centre = [bz_dims[0]/2.0, bz_dims[1]/2.0, bz_dims[2]/2.0]

    # Complain a bit if necessary
    if (ibz_data is None) or (sym_mats is None):
        raise ValueError('One of the parameters was not specified')
    if (constrain_to_bz == True) and \
        ((bz_dims is None) or (len(bz_dims) != 3)):
        raise ValueError('Need to specify the dimensions of the Brillouin Zone if points are to mapped back inside the Brillouin zone')
    if (len(sym_mats) == 0):
        raise ValueError('No symmetry matrices in list')

    # Build the full Brollouin zone
    full_bz = None
    for symmetry_matrix in sym_mats:
        kpoints_buffer = ibz_data.copy()
        # FIXME, rlv operation may not be needed for hexagonal lattices - investigate this
        kpoints_buffer[:,cols] = np.dot(kpoints_buffer[:,cols], np.linalg.inv(rlvs).transpose())
        kpoints_buffer = symmetry_matrix.map(kpoints_buffer, cols=cols)
        # Map transforms back into the unit cell if required
        if constrain_to_bz == True:
            i_col, j_col, k_col = kpoints_buffer[:,cols].transpose()
            max_i = bz_dims[0]/2.0 + bz_centre[0] 
            max_j = bz_dims[1]/2.0 + bz_centre[1]
            max_k = bz_dims[2]/2.0 + bz_centre[2]
            min_i = bz_centre[0] - bz_dims[0]/2.0
            min_j = bz_centre[1] - bz_dims[1]/2.0
            min_k = bz_centre[2] - bz_dims[2]/2.0
            # Map coords that are too big inside Brillouin zone
            while i_col.max() > max_i:
                i_col[i_col > max_i] -= bz_dims[0]
            while j_col.max() > max_j:
                j_col[j_col > max_j] -= bz_dims[1]
            while k_col.max() > max_k:
                k_col[k_col > max_k] -= bz_dims[2]
            # Map coords that are too small inside Brillouin zone
            while i_col.min() < min_i:
                i_col[i_col < min_i] += bz_dims[0]
            while j_col.min() < min_j:
                j_col[j_col < min_j] += bz_dims[1]
            while k_col.min() < min_k:
                k_col[k_col < min_k] += bz_dims[2]
            kpoints_buffer[:,cols[0]] = i_col
            kpoints_buffer[:,cols[1]] = j_col
            kpoints_buffer[:,cols[2]] = k_col
        # Add buffered data to the end of the list of k points
        if full_bz is None:
            full_bz = kpoints_buffer
        else:
            full_bz = np.concatenate((full_bz, kpoints_buffer))

    full_bz = remove_duplicates(full_bz, dp_tol=6, cols=cols)

    # Remove duplicate k points, with tolerance that allows up to one
    # million k points
    # Sort the results if required
    if sort_by is not None:
        sorted_inds = np.lexsort([full_bz[:,c] for c in sort_by])
        full_bz = full_bz[sorted_inds,:]

    # Returns a copy
    return full_bz.copy()


if __name__ == '__main__':
    # Set up the testing
    import doctest
    import os
    import sys
    # Try to find the test files as robustly as possible
    TiC_klist_filename = os.path.join(sys.path[0], '..', 'tests', 'TiC', 'TiC.klist')
    TiC_outputkgen_filename = os.path.join(sys.path[0], '..', 'tests', 'TiC', 'TiC.outputkgen')
    TiC_energy_filename = os.path.join(sys.path[0], '..', 'tests', 'TiC', 'TiC.energy')
    globs = {
        'TiC_klist_filename' : TiC_klist_filename,
        'TiC_outputkgen_filename' : TiC_outputkgen_filename,
        'TiC_energy_filename' : TiC_energy_filename,
        'expand_ibz' : expand_ibz,
    }
    doctest.testfile(os.path.join(sys.path[0], '..', 'tests', 'expand_ibz_test.txt'), globs=globs)
    doctest.testmod(globs=globs)
    
