'''
EnergyReader.py

A class for reading Wien2k .energy files
'''

__all__ = ['EnergyReader']

import numpy as np
from wien2k.Band import Band
from wien2k.errors import UnexpectedFileFormat

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
        line_num = 0
        for line in file_handle:
            line_num = line_num + 1
            if len(line) in kpoint_line_lengths:
                try:
                    i, j, k, k_id, unknown, num_bands, k_weight = \
                        [float(x.strip()) for x in line.split(' ') if x.strip() != '']
                except:
                    raise UnexpectedFileFormat('A non-float was parsed from a line identified as a k-point line (line: %d)' % line_num)
            elif len(line) == band_line_length:
                try:
                    band_num_id, energy = [x for x in line.split(' ') if x.strip() != '']
                    energy = float(energy)
                    band_num_id = int(band_num_id)
                except:
                    raise UnexpectedFileFormat('A non-number was parsed from a line identified as a band energy line (line: %d)' % line_num)
                while band_num_id > len(tmp_bands):
                    tmp_bands.append([])
                tmp_bands[band_num_id - 1].append([i, j, k, energy])
        file_handle.close()

        # Now cast into Band objects
        for i in range(len(tmp_bands)):
            self.bands.append(Band(str(i)))
            self.bands[i].data = np.array(tmp_bands[i])
        del tmp_bands # Frees up memory?
        #del data