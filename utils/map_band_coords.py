'''
Generates a Band.data array Mx(id,kx,ky,kz,energy) given a Band object and an array of points of form Mx(id,kx,ky,kz[,...])
'''

import numpy as np

def map_band_coords(band, new_coords):
    '''
    Generates a Band.data array Mx(id,kx,ky,kz,energy) given a Band object and
    an array of points of form Mx(id,kx,ky,kz[,...])
    '''
    sorted_new_coords = new_coords[:,:4].tolist()
    sorted_new_coords.sort()
    sorted_new_coords = np.array(sorted_new_coords)
    indexes = np.searchsorted(sorted_new_coords[:,0], band.k_point_ids)
    new_band = np.column_stack((sorted_new_coords[indexes,:], band.energies))
    return new_band
