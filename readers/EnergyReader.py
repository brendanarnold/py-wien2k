'''
EnergyReader.py

A class for reading Wien2k .energy files
'''

__all__ = ['EnergyReader']

import re
import numpy as np
from wien2k.Band import Band
from wien2k.errors import UnexpectedFileFormat

fmt = {
    # k Point line contains details on kpoint and is of format,
    # 0.000000000000E+00 0.000000000000E+00 0.000000000000E+00         1   455    50  1.0
    'kpoint_line_lengths' : (85, 88),
    'num_kpoint_vals' : 7,
    'kpoint_line' : re.compile('''
        \s*
        (-?[-+E\d\.]+)      # i
        \s+(-?[-+E\d\.]+)   # j
        \s+(-?[-+E\d\.]+)   # k
        \s+(-?[-+E\d\.]+)   # k point id
        \s+(-?[-+E\d\.]+)   # Unknown
        \s+(-?[-+E\d\.]+)   # Number of bands
        \s+(-?[-+E\d\.]+)   # Weight
        \s*
    ''', re.VERBOSE),
    # Band line contains details on band energy at a particular k point and is of format,
    #           3  -3.30918392086173
    'band_line_length' : 37,
    'num_band_vals' : 2,
    'band_line' : re.compile('''
        \s*
        (-?[-+E\d\.]+)
        \s+(-?[-+E\d\.]+)
        \s*
    ''', re.VERBOSE)
}

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
            if len(line) in fmt['kpoint_line_lengths']:
                try:
                    kpoint_vals = fmt['kpoint_line'].match(line).groups()
                    if len(kpoint_vals) != fmt['num_kpoint_vals']:
                        raise Exception
                    i, j, k, k_id, unknown, num_bands, k_weight = [float(x.strip()) for x in kpoint_vals]
                    num_bands = int(num_bands)
                    k_id = int(k_id)
                except:
                    raise UnexpectedFileFormat('A non-float was parsed from a line identified as a k-point line (line: %d)' % line_num)
            elif len(line) == fmt['band_line_length']:
                try:
                    band_vals = fmt['band_line'].match(line).groups()
                    if len(band_vals) != fmt['num_band_vals']:
                        raise Exception
                    band_num_id, energy = band_vals
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