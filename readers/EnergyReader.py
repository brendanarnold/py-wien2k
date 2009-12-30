'''
EnergyReader.py

A class for reading Wien2k .energy files
'''

__all__ = ['EnergyReader']

import re
import numpy as np
from wien2k.Band import Band
from wien2k.Kpoint import Kpoint
from wien2k.errors import UnexpectedFileFormat

fmt = {
    # Begins with expansions energy (E_J) for atom 'I'
    # WRITE(11,'(100(f9.5))') (E(J,I),J=1,LMAX)
    # Also has expansion energy for local orbitals (ELO_J) fro atom 'I'
    # WRITE(11,'(100(f9.5))') ((ELO(J,k,I),J=0,LOMAX),k=1,nloat)
    
    # k Point line contains details on kpoint and is of format,
    # 0.000000000000E+00 0.000000000000E+00 0.000000000000E+00         1   455    50  1.0
    # Actual format is
    # WRITE(11,'(3e19.12,a10,2i6,f5.1,a3)') SX, SY, SZ, KNAME, NV, NE, WEIGHT, IPGR
    # From line 46 of tapewf.f in SRC_lapw1
    'kpoint_line_lengths' : (85, 88),
    'num_kpoint_vals' : 7,
    'kpoint_line' : re.compile('''
        \s*
        (-?[-+E\d\.]+)      # i
        \s+(-?[-+E\d\.]+)   # j
        \s+(-?[-+E\d\.]+)   # k
        \s+(-?[-+E\d\.]+)   # k point id
        \s+(-?[-+E\d\.]+)   # Number plane waves (tot num considered recip. latt. vecs.) plus number of local orbitals
        \s+(-?[-+E\d\.]+)   # Number of bands
        \s+(-?[-+E\d\.]+)   # Weight
        \s*                 # Point group (not captured)
    ''', re.VERBOSE),
    # Band line contains details on band energy at a particular k point and is of format,
    #           3  -3.30918392086173
    # Actual format is
    # WRITE(11,*) I, E(I)
    # From line 50 of tapewf.f in SRC_lapw1
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
        self.kpoints = []
        tmp_bands = []
        file_handle = open(filename, 'r')
        line_num = 0
        for line in file_handle:
            line_num = line_num + 1
            if len(line) in fmt['kpoint_line_lengths']:
                try:
                    kpoint_vals = fmt['kpoint_line'].match(line).groups()
                    if len(kpoint_vals) != fmt['num_kpoint_vals']:
                        raise UnexpectedFileFormat('The specified number of values was not parsed from the k point line (line: %d)' % line_num)
                    i, j, k, kpoint_id, unknown, num_bands, k_weight = [float(x.strip()) for x in kpoint_vals]
                    num_bands = int(num_bands)
                    kpoint_id = int(kpoint_id)
                except:
                    raise UnexpectedFileFormat('A non-float was parsed from a line identified as a k-point line (line: %d)' % line_num)
                self.kpoints.append(Kpoint(kpoint_id, i, j, k))
            elif len(line) == fmt['band_line_length']:
                try:
                    band_vals = fmt['band_line'].match(line).groups()
                    if len(band_vals) != fmt['num_band_vals']:
                        raise UnexpectedFileFormat('A line describing a band contains less values than expected (line: %d)' % line_num)
                    band_id, energy = band_vals
                    energy = float(energy)
                    band_id = int(band_id)
                except:
                    raise UnexpectedFileFormat('A non-number was parsed from a line identified as a band energy line (line: %d)' % line_num)
                while band_id > len(tmp_bands):
                    tmp_bands.append([])
                tmp_bands[band_id - 1].append([kpoint_id, i, j, k, energy])
        file_handle.close()

        # Now cast into Band objects
        for i in range(len(tmp_bands)):
            self.bands.append(Band(id=i+1, data=np.array(tmp_bands[i])))
        del tmp_bands # Frees up memory?
