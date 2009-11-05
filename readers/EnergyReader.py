'''
EnergyReader.py

A class for reading Wien2k .energy files
'''

__all__ = ['EnergyReader']

import numpy as np
from wien2k.Band import Band

# Some file format information
kpoint_line_lengths = (85, 88)
band_line_length = 37

class EnergyReader(object):
    '''An object which reads WIEN2k .energy files and
    places the energies into bands

    Paramters,

    filename:               The filename of the .energy(so) file to parse
    spin_orbit_direction:   If is a spin orbit calculation, specify either 'up' or 'down' - default: None
    
    Results in,
    
    bands:                  A list of Band objects for this energy file
    '''
    def __init__(self, filename, spin_orbit_direction=None):
        self.filename = filename
        self.spin_orbit_direction = spin_orbit_direction
        self.bands = []
        tmp_bands = []
        file_handle = open(filename, 'r')
        for line in file_handle:
            if len(line) in kpoint_line_lengths:
                i, j, k, k_id, unknown, num_bands, k_weight = \
                    [float(x.strip()) for x in line.split(' ') if x.strip() != '']
            elif len(line) == band_line_length:
                band_num_id, energy = [x for x in line.split(' ') if x.strip() != '']
                energy = float(energy)
                band_num_id = int(band_num_id)
                while band_num_id > len(tmp_bands):
                    tmp_bands.append([])
                tmp_bands[band_num_id - 1].append([i, j, k, energy])
        file_handle.close()

        # Now cast into Band objects
        for i in range(len(tmp_bands)):
            self.bands.append(Band(str(i)))
            data = np.array(tmp_bands[i])
            self.bands[i].kpoints = data[:,:3] # This slice notation is supremely gay
            self.bands[i].energies = data[:,3]
        del tmp_bands # Frees up memory?
        #del data