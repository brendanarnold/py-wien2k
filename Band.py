'''Band.py

A class representing an energy band
'''

__all__ = ['Band']

import numpy as np

class Band(object):
    '''An object representing the bands in the calculation
    '''
    def __init__(self, id=None, character='', data=None):
        self.id = id
        self.character = character
        self.data = data
        
    def k_points(self):
        if self.data is not None:
            return self.data[:,:4]
        else:
            return None

    def k_point_ids(self):
        if self.data is not None:
            return np.array(self.data[:,0], dtype=int)
        else:
            return None

    def i_vals(self):
        if self.data is not None:
            return self.data[:,1]
        else:
            return None
        
    def j_vals(self):
        if self.data is not None:
            return self.data[:,2]
        else:
            return None
        
    def k_vals(self):
        if self.data is not None:
            return self.data[:,3]
        else:
            return None

    def energies(self):
        if self.data is not None:
            return self.data[:,4]
        else:
            return None

    # Overwrite the function labels with properties
    k_points = property(k_points)
    k_point_ids = property(k_point_ids)
    i_vals = property(i_vals)
    j_vals = property(j_vals)
    k_vals = property(k_vals)
    energies = property(energies)
    
    def map_coords(self, new_coords):
        '''
        Given an array of form Nx(id,kx,ky,kz[, ...]) will exchange the band
        objects coordinate system
        '''
        new_coords = new_coords[:,:4]
        ids = new_coords[:,0]
        # Have to cast to a list to sort whilst preserving rows - stupid
        energy_lookup = self.data[:,(0,4)].tolist()
        energy_lookup.sort()
        energy_lookup = np.array(energy_lookup)
        indexes = np.searchsorted(energy_lookup[:,0], ids)
        self.data = np.column_stack((new_coords, energy_lookup[indexes,1]))
