__all__ = ['reduce_to_irreducible_brillouin_zone']

import numpy as np

def reduce_to_irreducible_brillouin_zone(klist=None,
        struct_rdr=None, symmetry_matrices=None,
        kmesh=None, tolerance=None):
    '''
    Reduces a full Brillouin zone into its irreducible counterpart - c.f.
    expand_irreducible_brillouin_zone
    '''
    # Rudimetary check of the parameters
    if klist is kmesh is None:
        raise ValueError('One of klist or kmesh must be passed')
    if (klist is not None) and (kmesh is not None):
        raise ValueError('Only one of klist or kmesh can be passed')
    if struct_rdr is symmetry_matrices is None:
        raise ValueError('One of struct_rdr or symmetry_matrices must be passed')
    if (symmetry_matrices is not None) and (struct_rdr is not None):
        raise ValueError('Only one of symmetry_matrices or struct_rdr can be passed')
    klist = klist or kmesh.kpoints
    sms = symmetry_matrices or struct_rdr.symmetry_matrices
    kcoords = klist[:,1:4]
    tmp_klist = klist.copy()
    reduced_kcoords = np.zero((0, klist.shape[1]))
    for sm in sms:
        # Send all points to reduced zone
        tmp_klist[:,1:4] = np.dot(kcoords, np.linalg.inv(sm).transpose())
        # Copy back into the master klist array
        reduced_kcoords = np.concatenate(reduced_kcoords, tmp_klist)
    # Remove duplicate kcoords
    # TODO: The rest


