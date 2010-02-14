'''
Generates a Band.data array Mx(id,kx,ky,kz,energy) given a Band object and an array of points of form Mx(id,kx,ky,kz[,...])
'''

import numpy as np

def map_band_coords(band, new_coords):
    '''
    Generates a Band.data array Mx(id,kx,ky,kz,energy) given a Band object and
    an array of points of form Mx(id,kx,ky,kz[,...])
    '''
    new_coords = new_coords[:,:4]
    ids = new_coords[:,0]
    # Have to cast to a list to sort whilst preserving rows - stupid
    energy_lookup = band.data[:,(0,4)].tolist()
    energy_lookup.sort()
    energy_lookup = np.array(energy_lookup)
    indexes = np.searchsorted(energy_lookup[:,0], ids)
    new_band = np.column_stack((new_coords, energy_lookup[indexes,1]))
    return new_band
