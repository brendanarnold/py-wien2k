'''
KlistReader.py

A class to read in .klist files from WIEN2k
'''

__all__ = ['KlistReader']

import numpy as np
from wien2k.errors import UnexpectedFileFormat

file_terminator = 'END'
meta_data_line = 1

##class Kpoint(object):
##    def __init__(object, id, i, j, k, weight):
##        self.id = id
##        self.i = i
##        self.j = j
##        self.k = k
##        self.weight = weight

class KlistReader(object):
    '''A class to read in .klist files from WIEN2k'''
    def __init__(self, filename):
        self.filename = filename
        self.data = None
        file_handle = open(self.filename, 'r')
        line_num = 0
        tmp_data = []
        for line in file_handle:
            line_num = line_num + 1
            # Quit if reach file terminator string
            if line.strip().startswith(file_terminator):
                break
            vals = [x.strip() for x in line.split(' ') if x.strip != '']
            try:
                # i, j, k, denominator, weight
                tmp_data.append([float(x) for x in vals[1:6]])
            except:
                raise UnexpectedFileFormat('Cannot parse k-point data from file (line: %d))' % line_num)
            if line_num == meta_data_line:
                # TODO: Read in metadata as well
                pass
        self.data = np.array(tmp_data)
                
        file_handle.close()

    def i_vals(self):
        return self.data[:,0]/self.data[:4]
    
    def j_vals(self):
        return self.data[:,1]/self.data[:4]
    
    def k_vals(self):
        return self.data[:,2]/self.data[:4]
    
    def denominators(self):
        return self.data[:,3]
    
    def weights(self):
        return self.data[:,4]
    
    # Overwrite the function labels with properties
    i_vals = property(i_vals)
    j_vals = property(j_vals)
    k_vals = property(k_vals)