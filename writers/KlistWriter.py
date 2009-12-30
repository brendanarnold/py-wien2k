'''
KlistWriter.py

A class to write .klist file for use in WIEN2k calculations
'''

__all__ = ['KlistWriter']

import numpy as np

# FORTRAN format statement found in main.f on line 287
# FORMAT(I10,4I10,3f5.1,4x,i6,' k, div: (',3i3,')')  
# FORTRAN format statement found in main.f on line 285
# FORMAT(I10,4I10,f5.1)
# Note that file ends with line of FORTRAN format statment from line 286 of main.f
# format('END',/)

class KlistWriter(object):
    '''A class to write .klist files for WIEN2k calculations'''
    def __init__(self, filename=None):
        self.filename = filename
        self.data = None
        self.total_number_k_points = 0
        self.number_points_along_rlvs = []

    def write(self):
        filehandle = open(self.filename, 'w')
        if self.data != None:
            for line_num, k in enumerate(self.data):
                if line_num == 0:
                    if self.number_points_along_rlvs == []:
                        number_points_along_rlvs = [0, 0, 0]
                    else:
                        number_points_along_rlvs = self.number_points_along_rlvs
                    filehandle.write('%10d%10d%10d%10d%10d%5.1f%5.1f%5.1f    %6d k, div: (%3d%3d%3d)\n' % \
                        (k[0], k[1], k[2], k[3], k[4], k[5], -7.0, 1.5, self.total_number_k_points, \
                        number_points_along_rlvs[0], number_points_along_rlvs[1], number_points_along_rlvs[2]))
                else:
                    filehandle.write('%10d%10d%10d%10d%10d%5.1f\n' % (k[0], k[1], k[2], k[3], k[4], k[5]))
        filehandle.write('END\n')
        filehandle.close()

    def ids(self):
        if self.data != None:
            return np.array(self.data[:,0], dtype=int)
        else:
            return None

    def i_vals(self):
        if self.data != None:
            return self.data[:,1]/self.data[:,4]
        else:
            return None
    
    def j_vals(self):
        if self.data != None:
            return self.data[:,2]/self.data[:,4]
        else:
            return None
    
    def k_vals(self):
        if self.data != None:
            return self.data[:,3]/self.data[:,4]
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
