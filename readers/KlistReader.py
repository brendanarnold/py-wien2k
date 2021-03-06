'''
KlistReader.py

A class to read in .klist files from WIEN2k
'''

__all__ = ['KlistReader']

import numpy as np
import re
from wien2k.errors import UnexpectedFileFormat


# fmt describes the format of the WIEN2k .klist file
fmt = {
    # Note that file ends with line of FORTRAN format statment from line 286 of main.f
    # format('END',/)
    'file_terminator' : re.compile(r'\s*END\s*'),
    # Meta data line typically looks like,
    #         1         0         0         0        21  1.0 -7.0  1.5     10000 k, div: ( 21 21 21)
    # FORTRAN format statement found in main.f on line 287
    # FORMAT(I10,4I10,3f5.1,4x,i6,' k, div: (',3i3,')')  
    'meta_data_line_num' : 1,    # Specify a lne to make it quicker
    'num_meta_data_vals' : 12,    # Number of values read from the meta line - used to check parsing successful
    'meta_data_line' : re.compile(r'''
        \s*(-?[\d\.]+)      # k point id
        \s+(-?[\d\.]+)      # i' value
        \s+(-?[\d\.]+)      # j' value
        \s+(-?[\d\.]+)      # k' value
        \s+(-?[\d\.]+)      # denominator (i.e. divide i',j',k' by this to get true i,j,k)
        \s+(-?[\d\.]+)      # weight (no. of times this k point should be replicated in full BZ)
        \s+(-?[\d\.]+)      # Unknown number (i.e. '-7.0') - actually hardcoded!
        \s+(-?[\d\.]+)      # Unknown number (i.e. '1.5') - actually hardcoded!
        \s+(-?[\d\.]+)      # 'Number of k points in whole cell' (i.e. '10000') - poss. IBZ points expanded equals this number
        .*                  # Unknown string (i.e. 'k, div:')
        \(\s*([\d\s]{3})    # The number to divide the x recip. latt. vec. by to get the spacings of the mesh points in x dirn.
        \s*([\d\s]{3})      # The number to divide the y recip. latt. vec. by to get the spacings of the mesh points in y dirn
        \s*([\d\s]{3})\s*\) # The number to divide the z recip. latt. vec. by to get the spacings of the mesh points in z dirn
        .*
    ''', re.VERBOSE),
    # K point line typically looks like
    #         2         1         1         0        21  4.0
    # FORTRAN format statement found in main.f on line 285
    # FORMAT(I10,4I10,f5.1)
    'num_k_point_vals' : 6,   # Number of values read from k point line - used to check parsing successful
    'k_point_line' : re.compile(r'''
        \s+(-?[\d\.]+) # k point id
        \s+(-?[\d\.]+) # i' value
        \s+(-?[\d\.]+) # j' value
        \s+(-?[\d\.]+) # k' value
        \s+(-?[\d\.]+) # denominator (i.e. divide i',j',k' by this to get true i,j,k)
        \s+(-?[\d\.]+) # weight (no. of times this k point should be replicated in full BZ)
        \s*
    ''', re.VERBOSE)
}

xcd_band_fmt = {
    # The case.klist_band file (at least as generated by XCrysden 1.4.1) is of
    # a different format to the .klist file generated by WIEN2k. This
    # describes the file
    'file_terminator' : re.compile(r'\s*END\s*'),
    # Meta data line typically looks like,
    #N             0  -10    0   10  2.0-8.00 8.00    k-list generated by XCrySDen
    # Guess a format statement of
    # FORMAT('A10,4I5,3f5.1,A')  
    'meta_data_line_num' : 1,    # Specify a line to make it quicker
    'num_meta_data_vals' : 8,   # Number of values read from the meta line - used to check parsing successful
    'meta_data_line' : re.compile(r'''
        (.{10})               # Point name i.e. 'N' or 'GAMMA'
        ([-+\s\d]{5})         # 'i' value
        ([-+\s\d]{5})         # 'j' value
        ([-+\s\d]{5})         # 'k' value
        ([-+\s\d]{5})         # denominator (i.e. divide i',j',k' by this to get true i,j,k)
        ([-+\s\dDdEe\.]{5})   # weight (no. of times this k point should be replicated in full BZ)
        ([-+\s\dDdEe\.]{5})   # Unknown number (i.e. '8.0') - actually hardcoded!
        ([-+\s\dDdEe\.]{5})   # Unknown number (i.e. '8.0') - actually hardcoded!
        .*                    # Any leftover string i.e. 'k-list generated by XCrySDen' (not captured)
    ''', re.VERBOSE),
    # K point line typically looks like
    #              0   -8    0   10  2.0
    # Guess a format statement of
    # FORMAT(A10,4I5,f5.1)
    'num_k_point_vals' : 5,   # Number of values read from k point line - used to check parsing successful
    'k_point_line' : re.compile(r'''
        .{10}                 # Blanks (not captured)
        ([-+\s\d]{5})         # 'i' value
        ([-+\s\d]{5})         # 'j' value
        ([-+\s\d]{5})         # 'k' value
        ([-+\s\d]{5})         # denominator (i.e. divide i',j',k' by this to get true i,j,k)
        ([-+\s\dDdEe\.]{5})   # weight (no. of times this k point should be replicated in full BZ)
        .*                    # Any leftover string (not captured)
    ''', re.VERBOSE)    
}

class KlistReader(object):
    '''A class to read in .klist files from WIEN2k'''
    def __init__(self, filename):
        self.filename = filename
        self.data = None
        self.bz_shape = []
        self.is_bandlist = False
        self._load_values()
    def _load_values(self):
        # Loads the object with values from the .klist file
        file_handle = open(self.filename, 'r')
        line_num = 0
        tmp_data = []
        for line in file_handle:
            line_num = line_num + 1
            # Quit if reach file terminator string
            ##if line.strip().startswith(file_terminator):
            if fmt['file_terminator'].match(line):
                break
            # Skip if line empty
            if line.strip() == '':
                continue
            # Try to determine if a file is a .klist or a .klist_band
            if line_num == fmt['meta_data_line_num']:
                if fmt['meta_data_line'].match(line) is None:
                    assert line_num == xcd_band_fmt['meta_data_line_num']
                    if xcd_band_fmt['meta_data_line'].match(line) is None:
                        raise UnexpectedFileFormat('Could not determine if file is .klist or .klist_band')
                    else:
                        self.is_bandlist = True
                else:
                    self.is_bandlist = False
            if self.is_bandlist == True:
                # Parse out meta data line
                if line_num == xcd_band_fmt['meta_data_line_num']:
                    # Extract all the matches to the regular expression
                    meta_vals = xcd_band_fmt['meta_data_line'].match(line).groups()
                    if len(meta_vals) != xcd_band_fmt['num_meta_data_vals']:
                        raise UnexpectedFileFormat('Meta line parsed wrong number of elements - check that the file is correct or update format (line: %d))' % line_num)
                    tmp_data.append([0.0] + [float(x) for x in meta_vals[1:6]]) # Line also contains a k point - add this
                # Otherwise parse out the k point line
                else:
                    # Extract all the matches to the regular expression
                    k_point_vals = xcd_band_fmt['k_point_line'].match(line).groups()
                    if len(k_point_vals) != xcd_band_fmt['num_k_point_vals']:
                        raise UnexpectedFileFormat('k point line parsed wrong number of elements - check that the file is correct or update format (line: %d))' % line_num)
                    # i', j', k', denominator, weight
                    tmp_data.append([0.0] + [float(x) for x in k_point_vals])
            else:
                # Parse out meta data line
                if line_num == fmt['meta_data_line_num']:
                    # Extract all the matches to the regular expression
                    meta_vals = fmt['meta_data_line'].match(line).groups()
                    if len(meta_vals) != fmt['num_meta_data_vals']:
                        raise UnexpectedFileFormat('Meta line parsed wrong number of elements - check that the file is correct or update format (line: %d))' % line_num)
                    self.bz_shape = tuple([int(x) for x in meta_vals[9:12]])
                    tmp_data.append([float(x) for x in meta_vals[:6]]) # Line also contains a k point - add this
                # Otherwise parse out the k point line
                else:
                    # Extract all the matches to the regular expression
                    k_point_vals = fmt['k_point_line'].match(line).groups()
                    if len(k_point_vals) != fmt['num_k_point_vals']:
                        raise UnexpectedFileFormat('k point line parsed wrong number of elements - check that the file is correct or update format (line: %d))' % line_num)
                    # i', j', k', denominator, weight
                    tmp_data.append([float(x) for x in k_point_vals[:6]])
        self.data = np.array(tmp_data)
        file_handle.close()

    def ids(self):
        if self.data != None:
            return np.array(self.data[:,0], dtype=int)
        else:
            return None

    def i_vals(self):
        if self.data != None:
            return self.data[:,1]
        else:
            return None
    
    def j_vals(self):
        if self.data != None:
            return self.data[:,2]
        else:
            return None
    
    def k_vals(self):
        if self.data != None:
            return self.data[:,3]
        else:
            return None
    
    def denominators(self):
        if self.data != None:
            return self.data[:,4]
        else:
            return None
    
    def weights(self):
        if self.data != None:
            return self.data[:,5]
        else:
            return None

    def k_points(self):
        return self.data[:,:4]
    
    # Overwrite the function labels with properties
    ids = property(ids)
    i_vals = property(i_vals)
    j_vals = property(j_vals)
    k_vals = property(k_vals)
    denominators = property(denominators)
    weights = property(weights)
    k_points = property(k_points)
    
