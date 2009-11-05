'''
OutputkgenReader.py

Reads from the .outputkgen file which contains symmetry matrices
'''

__all__ = ['OutputkgenReader']

import numpy as np
##from wien2k.errors import UnexpectedFileFormat

skip_lines = 11
matrix_header_test_string = 'SYMMETRY MATRIX NR.' # A string that denotes the line before the symmetry matrices

class OutputkgenReader(object):
    '''A class which reads from the .outputkgen file which contains symmetry vectors'''
    def __init__(self, filename):
        self.filename = filename
        self.symmetry_matrices = []
        # Read in the symmetry matrices
        file_handle = open(filename, 'r')
        line_num = 0
        while True:
            line = file_handle.readline()
            line_num = line_num + 1
            # Test for EOF
            if line == '':
                break
            # Skip a certain number of lines - first row of matrices contain unknown data
            if line_num < skip_lines:
                continue
            # Test to see if a set of matrices is coming up
            if line.strip().startswith(matrix_header_test_string):
                tmp_matrices = []
                # Read the matrices off the next three lines into a buffer
                for i in [1,2,3]:
                    line = file_handle.readline()
                    vals = [int(x.strip()) for x in line.split(' ') if x.strip() != '']
                    # Count how many matrices required (may not be three)
                    num_matrices = len(vals) / 3
                    for j in range(num_matrices):
                        tmp_matrices.append([])
                    for j in range(num_matrices):
                        tmp_matrices[j].append(vals[3*j:3*(j+1)])
                # Read the buffered matrices into the object as a Numpy array object only
                # if the determinant is not zero
                for i in range(num_matrices):
                    tmp_matrix = np.array(tmp_matrices[i])
                    if np.linalg.det(tmp_matrix) == 0:
                        continue
                    else:
                        self.symmetry_matrices.append(tmp_matrix)
##        if len(self.symmetry_matrices) == 0:
##            raise UnexpectedFileFormat('No symmetry matrices found in file: %s' % self.filename)
