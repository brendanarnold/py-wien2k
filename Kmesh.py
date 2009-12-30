'''
Module containg the Kmesh class
'''
import numpy as np

class Kmesh(object):
    '''
    Pass in a band_data array (an Nx5 array of id, i, j, k, energy values, as
    from a Band object) containing k points evenly spaced over a Cartesian grid
    (although points may be missing)

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
        self._build_mesh(band_data)

    def i_vals(self):
        return np.arange(self.mesh.shape[0]) * self.i_series_spacing + self.i_series_offset

    def j_vals(self):
        return np.arange(self.mesh.shape[1]) * self.j_series_spacing + self.k_series_offset

    def k_vals(self):
        return np.arange(self.mesh.shape[2]) * self.j_series_spacing + self.k_series_offset

    def _build_mesh(self, band_data):
        '''Builds a 3D array of energies'''
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
        # Find the dimensions of the 3D array
        i_dimension = int((i_vals - self.i_series_offset) / self.i_series_spacing)
        j_dimension = int((j_vals - self.j_series_offset) / self.j_series_spacing)
        k_dimension = int((k_vals - self.k_series_offset) / self.k_series_spacing)
        self.mesh = np.ma.zeros(i_dimension, j_dimension, k_dimension)
        self.mesh.mask = True
        k_point_indexes = np.zeros((4, band_data.shape[1]), dtype=int)
        k_point_indexes[:,0] = (band_data[:,1] - self.i_series_offset) / self.i_series_spacing
        k_point_indexes[:,1] = (band_data[:,2] - self.j_series_offset) / self.j_series_spacing
        k_point_indexes[:,2] = (band_data[:,3] - self.k_series_offset) / self.k_series_spacing
        k_point_indexes[:,3] = band_data[:,4]
        for k in k_point_indexes:
            self.mesh[k[0], k[1], k[2]] = k[3]
            self.mesh.mask[k[0], k[1], k[2]] = False

    def _find_arithmetic_series_formula(self, series):
        '''Returns the formula for an incomplete arithemtic progression
        i.e. returns a and b for x_n = a + b*n'''
        differences = series[1:] - series[:-1]
        a = series[0]
        b = differences.min()
        return (a, b)
