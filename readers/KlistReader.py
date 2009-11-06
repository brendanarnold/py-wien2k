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
    'file_terminator' : re.compile(r'\s*END\s*'),
    # Meta data line typically looks like,
    #         1         0         0         0        21  1.0 -7.0  1.5     10000 k, div: ( 21 21 21)
    'meta_data_line' : 1,    # Specify a lne to make it quicker
    'num_meta_data_vals' : 12,    # Number of values read from the meta line - used to check parsing successful
    'meta_data_line' : re.compile(r'''
        \s*(-?[\d\.]+)      # k point id
        \s+(-?[\d\.]+)      # i' value
        \s+(-?[\d\.]+)      # j' value
        \s+(-?[\d\.]+)      # k' value
        \s+(-?[\d\.]+)      # denominator (i.e. divide i',j',k' by this to get true i,j,k)
        \s+(-?[\d\.]+)      # weight (no. of times this k point should be replicated in full BZ)
        \s+(-?[\d\.]+)      # Unknown number (i.e. '-7.0')
        \s+(-?[\d\.]+)      # Unknown number (i.e. '1.5')
        \s+(-?[\d\.]+)      # Unknown number (i.e. '10000') - poss. requested number k points?
        .*                  # Unknown string (i.e. 'k, div:')
        \(\s*(-?[\d\.]+)    # Dimensions(?)
        \s+(-?[\d\.]+)      # Dimensions(?)
        \s+(-?[\d\.]+)\s*\) # Dimensions(?)
        .*
    ''', re.VERBOSE),
    # K point line typically looks like
    #         2         1         1         0        21  4.0
    'num_kpoint_vals' : 6,   # Number of values read from k point line - used to check parsing successful
    'kpoint_line' : re.compile(r'''
        \s+(-?[\d\.]+) # k point id
        \s+(-?[\d\.]+) # i' value
        \s+(-?[\d\.]+) # j' value
        \s+(-?[\d\.]+) # k' value
        \s+(-?[\d\.]+) # denominator (i.e. divide i',j',k' by this to get true i,j,k)
        \s+(-?[\d\.]+) # weight (no. of times this k point should be replicated in full BZ)
        \s*
    ''', re.VERBOSE)
}

class KlistReader(object):
    '''A class to read in .klist files from WIEN2k'''
    def __init__(self, filename):
        self.filename = filename
        self.data = None
        self.dimensions = []
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
            # Parse out meta data line
            if line_num == fmt['meta_data_line']:
                # Extract all the matches to the regular expression
                try:
                    meta_vals = fmt['meta_data_line'].match(line).groups()
                    if len(meta_vals) != fmt['num_meta_data_vals']:
                        raise Exception
                except:
                    raise UnexpectedFileFormat('Meta line parsed wrong number of elements - check that the file is correct or update format (line: %d))' % line_num)
                self.dimensions = [float(x) for x in meta_vals[9:12]]
                tmp_data.append([float(x) for x in meta_vals[1:6]]) # Line also contains a k point - add this
            # Otherwise parse out the k point line
            else:
                # Extract all the matches to the regular expression
                try:
                    kpoint_vals = fmt['kpoint_line'].match(line).groups()
                    if len(kpoint_vals) != fmt['num_kpoint_vals']:
                        raise Exception
                except:
                    raise UnexpectedFileFormat('k point line parsed wrong number of elements - check that the file is correct or update format (line: %d))' % line_num)
                # i', j', k', denominator, weight
                tmp_data.append([float(x) for x in kpoint_vals[1:6]])
        self.data = np.array(tmp_data)
        file_handle.close()

    def i_vals(self):
        return self.data[:,0]/self.data[:,4]
    
    def j_vals(self):
        return self.data[:,1]/self.data[:,4]
    
    def k_vals(self):
        return self.data[:,2]/self.data[:,4]
    
    def denominators(self):
        return self.data[:,3]
    
    def weights(self):
        return self.data[:,4]
    
    # Overwrite the function labels with properties
    i_vals = property(i_vals)
    j_vals = property(j_vals)
    k_vals = property(k_vals)
    denominators = property(denominators)
    weights = property(weights)
    