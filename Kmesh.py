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
    corresponding i, j and k values are also created.

    Retrieve some band data and expand to full Brillouin zone

    >>> import wien2k
    >>> from wien2k.utils import expand_ibz
    >>> b = wien2k.EnergyReader('TiC.energy').bands[0]
    >>> sym_mats = wien2k.StructReader('TiC.struct').sym_mats
    >>> full_zone = expand_ibz(ibz_data=b.data, sym_mats=sym_mats)
    
    Now use the Kmesh object

    >>> km = Kmesh(full_zone)
    >>> km.shape
        (10, 10, 10)
    >>> km2 = km[:,:,5]
    >>> km2.shape
        (10, 10)
    '''
    def __init__(self, band_data=None):
        self.i_spacing = None
        self.j_spacing = None
        self.k_spacing = None
        self.i_offset = None
        self.j_offset = None
        self.k_offset = None
        # Allow a band to be passed in
        if band_data is not None:
            if isinstance(band_data, wien2k.Band):
                band_data = band_data.data
            self._build_mesh(band_data)
        else:
            self.energies = None
            self.ids = None

    def i_vals(self):
        return np.arange(self.energies.shape[0]) * \
          self.i_spacing + self.i_offset
    i_vals = property(i_vals)

    def j_vals(self):
        return np.arange(self.energies.shape[1]) * \
          self.j_spacing + self.j_offset
    j_vals = property(j_vals)

    def k_vals(self):
        return np.arange(self.energies.shape[2]) * \
          self.k_spacing + self.k_offset
    k_vals = property(k_vals)

    def shift_centre(self, shift_by, rel=True):
        '''
        Makes the assumption that the mesh is a repeating cell and moves
        the values over by a specified amount

        INPUT:

        shift_by    An iterable containing the i, j,k values that define
                    the shift vector
        rel         Use normalised units if True (i.e. 0.5 moves the
                    centre to one edge of the cell), otherwise use k
                    values (default: True)

        BEWARE: This is limited by the fineness of the mesh, it will
        round to the nearest appropriate mesh point, i.e. if the centre is
        shifted by half and there are an odd number of points in the mesh then
        the choice of the centre point will be always rounded up
        '''
        if self.energies.shape != self.ids.shape:
            raise ValueError('Energy mesh and id mesh do not match in size! Cannot shift ...')
        spacings = np.array([self.i_spacing, self.j_spacing, self.k_spacing])
        if rel == True:
            shift_places = np.ceil(np.array(shift_by) * (np.array(self.energies.shape)-1)).astype(int)
        else:
            shift_places = np.ceil(np.array(shift_by) / spacings).astype(int)
        # Roll the matrices around
        self.energies = np.roll(self.energies, shift_places[0], axis=0)
        self.energies = np.roll(self.energies, shift_places[1], axis=1)
        self.energies = np.roll(self.energies, shift_places[2], axis=2)
        self.ids = np.roll(self.ids, shift_places[0], axis=0)
        self.ids = np.roll(self.ids, shift_places[1], axis=1)
        self.ids = np.roll(self.ids, shift_places[2], axis=2)
        # Adjust the offsets
        self.i_offset = self.i_offset + shift_places[0] * self.i_spacing
        self.j_offset = self.j_offset + shift_places[1] * self.j_spacing
        self.k_offset = self.k_offset + shift_places[2] * self.k_spacing

    def indexes(self):
        '''
        Returns an Nx4 array of the id,i,j,k,energy values as indexes
        '''
        ind_list = np.array([ind + (en,) for ind, en in np.ndenumerate(self.energies) if (self.energies.mask[ind] == False)])
        ids = np.array([id for ind, id in np.ndenumerate(self.ids) if (self.ids.mask[ind] == False)])
        return np.concatenate((ids.reshape((ids.shape[0], 1)), ind_list), axis=1)
    indexes = property(indexes)

    def kpoints(self):
        '''
        Returns an Nx4 array of the id,i,j,k,energy values as k values
        '''
        ind_list = self.indexes
        ind_list[:,1] = ind_list[:,1] * self.i_spacing + self.i_offset
        ind_list[:,2] = ind_list[:,2] * self.j_spacing + self.j_offset
        ind_list[:,3] = ind_list[:,3] * self.k_spacing + self.k_offset
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
        self.i_offset, self.i_spacing = \
            self._find_arithmetic_series_formula(i_vals)
        self.j_offset, self.j_spacing = \
            self._find_arithmetic_series_formula(j_vals)
        self.k_offset, self.k_spacing = \
            self._find_arithmetic_series_formula(k_vals)
        # Find the dimensions of the 3D array - bearing in mind that there may
        # not be a point at every spacing
        if self.i_spacing != 0.0:
            i_dimension = int((i_vals.max() - i_vals.min()) / self.i_spacing) + 1
        else:
            i_dimension = 1
        if self.j_spacing != 0.0:
            j_dimension = int((j_vals.max() - j_vals.min()) / self.j_spacing) + 1
        else:
            j_dimension = 1
        if self.k_spacing != 0.0:
            k_dimension = int((k_vals.max() - k_vals.min()) / self.k_spacing) + 1
        else:
            k_dimension = 1
        self.energies = np.ma.zeros((i_dimension, j_dimension, k_dimension))
        self.ids = np.ma.zeros((i_dimension, j_dimension, k_dimension))
        self.energies.mask = True
        self.ids.mask = True
        # Convert the band_data co-ordinates to indexes
##         if self.i_spacing == 0:
##             band_data[:,1] = 0
##         else:
##             band_data[:,1] = (band_data[:,1]-self.i_offset)/self.i_spacing
##         if self.j_spacing == 0:
##             band_data[:,2] = 0
##         else:
##             band_data[:,2] = (band_data[:,2]-self.j_offset)/self.j_spacing
##         if self.k_spacing == 0:
##             band_data[:,3] = 0
##         else:
##             band_data[:,3] = (band_data[:,3]-self.k_offset)/self.k_spacing
        # Now populate the mesh
        for k in band_data:
            id, i, j, k, energy = k
            if self.i_spacing == 0:
                i = 0
            else:
                i = (i-self.i_offset)/self.i_spacing
            if self.j_spacing == 0:
                j = 0
            else:
                j = (j-self.j_offset)/self.j_spacing
            if self.k_spacing == 0:
                k = 0
            else:
                k = (k-self.k_offset)/self.k_spacing
            id = int(id)
            i = int(i)
            j = int(j)
            k = int(k)
            self.energies[i, j, k] = energy
            self.ids[i, j, k] = id
            self.energies.mask[i, j, k] = False
            self.ids.mask[i, j, k] = False

    def _find_arithmetic_series_formula(self, series):
        '''Returns the formula for an incomplete arithmetic progression
        i.e. returns a and b for x_n = a + b*n'''
        # Test if series is iterable, if not then cast as np array
        try:
            tst = iter(series)
        except TypeError:
            series = np.array([series])
        copy_series = np.unique(series[:])
        copy_series.sort()
        if len(copy_series) > 1:
            differences = copy_series[1:] - copy_series[:-1]
        else:
            differences = np.array([0.0])
        # Pick a difference that is at least slightly sensible
        cut_off = copy_series.min() * 1e-10
        i = 0
        a = copy_series[0]
        for diff in differences:
            if diff > cut_off:
                b = diff
                break
        return (a, b)

    def shape(self):
        return self.energies.shape
    shape = property(shape)

    def __len__(self):
        return len(self.energies)

    def __getitem__(self, *args):
        '''
        Allow Kmesh objects to be generated from slices called on the
        object itself (raw data can be got from the attributes directly)
        '''
        slices = args[0]
        new_km = Kmesh()
        # Copy across the meshes
        new_km.energies = self.energies[slices].copy()
        new_km.ids = self.ids[slices].copy()
        # Redefine the spacings and offsets
        if type(slices) != tuple:
            slices = tuple([slices])
        if len(slices) > 0:
            i_vals = self.i_vals
            i_vals = i_vals[slices[0]]
            new_km.i_offset, new_km.i_spacing = self._find_arithmetic_series_formula(i_vals)
        if len(slices) > 1:
            j_vals = self.j_vals
            j_vals = j_vals[slices[1]]
            new_km.j_offset, new_km.j_spacing = self._find_arithmetic_series_formula(j_vals)
        if len(slices) > 2:
            k_vals = self.k_vals
            k_vals = k_vals[slices[0]]
            new_km.k_offset, new_km.k_spacing = self._find_arithmetic_series_formula(k_vals)
        return new_km

    def query(self):
        '''
        Returns a bunch of useful stuff when using interactively
        '''
        out = '''
        Kmesh object:
          i_spacing = %f
          i_offset = %f
          j_spacing = %f
          j_offset = %f
          k_spacing = %f
          k_offset = %f
          energies.shape = %s
          ids.shape = %s
          len(kpoints) = %d
          energy min = %f
          energy max = %f
          id min = %d
          id max = %d
          total points = %d
          num masked = %d
        ''' % (
            self.i_spacing,
            self.i_offset,
            self.j_spacing,
            self.j_offset,
            self.k_spacing,
            self.k_offset,
            str(self.energies.shape),
            str(self.ids.shape),
            len(self.kpoints),
            self.energies.min(),
            self.energies.max(),
            self.ids.min(),
            self.ids.max(),
            self.energies.shape[0] * self.energies.shape[1] * self.energies.shape[2],
            self.energies.mask.sum()
        )
        print out
    q = property(query)
        

if __name__ == '__main__':
    '''Do some testing ...'''
    import doctest
    import os
    import sys
    orig_dir = os.getcwd()
    os.chdir(os.path.join(sys.path[0], 'tests', 'TiC'))
    doctest.testmod()
    doctest.testfile(os.path.join('..', 'Kmesh_test.txt'))
    os.chdir(orig_dir)
