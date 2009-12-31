'''
OutputkgenReader.py

Reads from the .outputkgen file which contains symmetry matrices
'''

__all__ = ['OutputkgenReader']

import numpy as np
from wien2k.errors import UnexpectedFileFormat
from wien2k.SymmetryMatrix import SymmetryMatrix

fmt = {
    # Skip these line until find out what they are
    'skip_lines' : 11,
    # A string that denotes the line before the symmetry matrices
    'matrix_header_test_string' : 'SYMMETRY MATRIX NR.',
    # A string that denotes the line before the rectangular lattice vectors
    'rec_lattice_vectors_header_test_string' : 'G1        G2        G3'
}

 


class OutputkgenReader(object):
    '''A class which reads from the .outputkgen file which contains symmetry vectors'''
    def __init__(self, filename):
        self.filename = filename
        self.symmetry_matrices = []
        self.rectangular_lattice_vectors = []
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
            if line_num < fmt['skip_lines']:
                continue
            # Test to see if a set of matrices is coming up
            if line.strip().startswith(fmt['matrix_header_test_string']):
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
                        self.symmetry_matrices.append(SymmetryMatrix(matrix = tmp_matrix))
            # Test to see if is line preceding the rectangular lattice vectors
            if line.strip().startswith(fmt['rec_lattice_vectors_header_test_string']):
                g1 = []
                g2 = []
                g3 = []
                for i in [1,2,3]:
                    line_num = line_num + 1
                    line = file_handle.readline()
                    try:
                        vals = [float(x.strip()) for x in line.split(' ') if x.strip() != '']
                        g1.append(vals[0])
                        g2.append(vals[1])
                        g3.append(vals[2])
                    except:
                        raise UnexpectedFileFormat('Could not parse rectangular lattice vectors from line (line: %d)' % line_num)
                self.rectangular_lattice_vectors.append(np.array(g1))
                self.rectangular_lattice_vectors.append(np.array(g2))
                self.rectangular_lattice_vectors.append(np.array(g3))
                # Only take the first set of lattice vectors (second set is just 2*pi*a)
                break
            # TODO: Parse out other values
        file_handle.close()
