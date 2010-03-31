'''
Module containg the Kmesh class
'''

__all__ = ['Kmesh']

import numpy as np
import wien2k

class Kmesh(object):
    '''
    Pass in a band istance or a band_data array (an Nx5 array of id, i, j, k,
    energy values, as from a Band object) containing k points evenly spaced
    over a Cartesian grid (although points may be missing)

    A full cubic 3D mesh of energies is created (missing points are 'masked')
    correspdonging i, j and k values are also created.
    '''
    def __init__(self, band_data):
        self.i_series_spacing = None
        self.j_series_spacing = None
        self.k_series_spacing = None
        self.i_series_offset = None
        self.j_series_offset = None
        self.k_series_offset = None
        # Allow a band to be passed in
        if isinstance(band_data, wien2k.Band):
            band_data = band_data.data
        self._build_mesh(band_data)

    def i_vals(self):
        return np.arange(self.mesh.shape[0]) * \
          self.i_series_spacing + self.i_series_offset
    i_vals = property(i_vals)

    def j_vals(self):
        return np.arange(self.mesh.shape[1]) * \
          self.j_series_spacing + self.j_series_offset
    j_vals = property(j_vals)

    def k_vals(self):
        return np.arange(self.mesh.shape[2]) * \
          self.k_series_spacing + self.k_series_offset
    k_vals = property(k_vals)

    def indexes(self):
        '''
        Returns an Nx4 array of the i,j,k,energy values as indexes
        '''
        ind_list = np.array([ind + (en,) for ind, en in np.ndenumerate(self.mesh)])
        ids = np.array([id for ind, id in np.ndenumerate(self.mesh_ids)])
        return np.concatenate((ids.reshape((ids.shape[0], 1)), ind_list), axis=1)
    indexes = property(indexes)

    def kpoints(self):
        '''
        Returns an Nx4 array of the i,j,k,energy values as k values
        '''
        ind_list = self.indexes
        ind_list[:,1] = ind_list[:,1] * self.i_series_spacing + self.i_series_offset
        ind_list[:,2] = ind_list[:,2] * self.j_series_spacing + self.j_series_offset
        ind_list[:,3] = ind_list[:,3] * self.k_series_spacing + self.k_series_offset
        return ind_list
    kpoints = property(kpoints)

    def centre_point(self):
        ''' Returns the cente point of the zone '''
        i_centre = (self.i_vals.max() + self.i_vals.min()) / 2.0
        j_centre = (self.j_vals.max() + self.j_vals.min()) / 2.0
        k_centre = (self.k_vals.max() + self.k_vals.min()) / 2.0
        return np.array([i_centre, j_centre, k_centre])
    centre_point = property(centre_point)

    def _build_mesh(self, band_data):
        '''Builds a 3D array of energies'''
        band_data = band_data.copy()
        i_vals = np.unique(band_data[:,1])
        j_vals = np.unique(band_data[:,2])
        k_vals = np.unique(band_data[:,3])
        j_vals.sort()
        j_vals.sort()
        k_vals.sort()
        self.i_series_offset, self.i_series_spacing = \
            self._find_arithmetic_series_formula(i_vals)
        self.j_series_offset, self.j_series_spacing = \
            self._find_arithmetic_series_formula(j_vals)
        self.k_series_offset, self.k_series_spacing = \
            self._find_arithmetic_series_formula(k_vals)
        # Find the dimensions of the 3D array - bearing in mind that there may
        # not be a point at every spacing
        if self.i_series_spacing != 0.0:
            i_dimension = int((i_vals.max() - i_vals.min()) / self.i_series_spacing) + 1
        else:
            i_dimension = 1
        if self.j_series_spacing != 0.0:
            j_dimension = int((j_vals.max() - j_vals.min()) / self.j_series_spacing) + 1
        else:
            j_dimension = 1
        if self.k_series_spacing != 0.0:
            k_dimension = int((k_vals.max() - k_vals.min()) / self.k_series_spacing) + 1
        else:
            k_dimension = 1
        self.mesh = np.ma.zeros((i_dimension, j_dimension, k_dimension))
        self.mesh_ids = np.ma.zeros((i_dimension, j_dimension, k_dimension))
        self.mesh.mask = True
        self.mesh_ids.mask = True
        # Convert the band_data co-ordinates to indexes
##         if self.i_series_spacing == 0:
##             band_data[:,1] = 0
##         else:
##             band_data[:,1] = (band_data[:,1]-self.i_series_offset)/self.i_series_spacing
##         if self.j_series_spacing == 0:
##             band_data[:,2] = 0
##         else:
##             band_data[:,2] = (band_data[:,2]-self.j_series_offset)/self.j_series_spacing
##         if self.k_series_spacing == 0:
##             band_data[:,3] = 0
##         else:
##             band_data[:,3] = (band_data[:,3]-self.k_series_offset)/self.k_series_spacing
        # Now populate the mesh
        for k in band_data:
            id, i, j, k, energy = k
            if self.i_series_spacing == 0:
                i = 0
            else:
                i = (i-self.i_series_offset)/self.i_series_spacing
            if self.j_series_spacing == 0:
                j = 0
            else:
                j = (j-self.j_series_offset)/self.j_series_spacing
            if self.k_series_spacing == 0:
                k = 0
            else:
                k = (k-self.k_series_offset)/self.k_series_spacing
            id = int(id)
            i = int(i)
            j = int(j)
            k = int(k)
            self.mesh[i, j, k] = energy
            self.mesh_ids[i, j, k] = id
            self.mesh.mask[i, j, k] = False
            self.mesh_ids.mask[i, j, k] = False

    def _find_arithmetic_series_formula(self, series):
        '''Returns the formula for an incomplete arithmetic progression
        i.e. returns a and b for x_n = a + b*n'''
        copy_series = np.unique(series[:])
        copy_series.sort()
        if len(copy_series) > 1:
            differences = copy_series[1:] - copy_series[:-1]
        else:
            differences = np.array([0.0])
        a = copy_series[0]
        b = differences.min()
        return (a, b)
