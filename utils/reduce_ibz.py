__all__ = ['reduce_ibz']

import numpy as np
from wien2k.utils import remove_duplicates

def reduce_ibz(klist=None,
        struct_rdr=None, sym_mats=None,
        kmesh=None, tolerance=None):
    '''
    Reduces a full Brillouin zone into its irreducible counterpart - c.f.
    expand_ibz
    '''
    # Rudimetary check of the parameters
    if klist is kmesh is None:
        raise ValueError('One of klist or kmesh must be passed')
    if (klist is not None) and (kmesh is not None):
        raise ValueError('Only one of klist or kmesh can be passed')
    if struct_rdr is sym_mats is None:
        raise ValueError('One of struct_rdr or sym_mats must be passed')
    if (sym_mats is not None) and (struct_rdr is not None):
        raise ValueError('Only one of sym_mats or struct_rdr can be passed')
    if klist is None:
        klist = kmesh.kpoints
    if sym_mats is None:
        sym_mats = struct_rdr.sym_mats
    # Check to see if the ids are appended, is not then set a flag
    no_ids = False
    if klist.shape[1] == 3:
        no_ids = True
    if no_ids == True:
        kcoords = klist[:,:]
    else:
        kcoords = klist[:,1:4]
    tmp_klist = klist.copy()
    reduced_kcoords = None
    for sm in sym_mats:
        # Send all points to reduced zone
        if no_ids == True:
            tmp_klist[:,:3] = np.dot(kcoords, np.linalg.inv(sm.matrix).transpose())
        else:
            tmp_klist[:,1:4] = np.dot(kcoords, np.linalg.inv(sm.matrix).transpose())
        # Copy back into the master klist array
        if reduced_kcoords is None:
            reduced_kcoords = tmp_klist
        else:
            reduced_kcoords = np.concatenate((reduced_kcoords, tmp_klist))
    # Remove duplicate kcoords
    if no_ids == True:
        reduced_kcoords = remove_duplicates(reduced_kcoords, cols=[0,1,2])
    else:
        reduced_kcoords = remove_duplicates(reduced_kcoords, cols=[1,2,3])
    return reduced_kcoords
