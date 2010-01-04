'''
OutputkgenReader.py

Reads from the .outputkgen file which contains symmetry matrices
'''

__all__ = ['OutputkgenReader']

import numpy as np
from wien2k.errors import UnexpectedFileFormat
from wien2k.SymmetryMatrix import SymmetryMatrix
import wien2k.CONSTANTS as CNST

fmt = {
    'ortho_line_num' : 1,
    'r1_vector_line_num' : 2, 
    'iarb_line_num' : 5,
    'point_group_symmetries_line_num' : 6,
    'symmetry_matrix_header_line_startswith' : 'SYMMETRY MATRIX NR.',
    'reciprocal_lattice_vectors_header_line_startswith' : 'G1        G2        G3',
    'length_reciprocal_lattice_vectors_line_startswith' : 'length of reciprocal lattice vectors:',
    'submesh_shift_line_startswith' : 'SUBMESH SHIFTED; SHIFT:',
    'number_mesh_points_line_startswith' : 'NO. OF MESH POINTS IN THE BRILLOUIN ZONE =',
    'reciprocal_lattice_vector_intervals_line_startswith' : 'DIVISION OF RECIPROCAL LATTICE VECTORS (INTERVALS)=',
    'num_inequiv_k_points_line_startswith' : 'NO. OF INEQUIVALENT K-POINTS',
    'len_bloch_vector_line' : 74,
    'tetrahedra_to_sort_line_startswith' : 'tetrahedra to sort:',
    'owork_line_startswith' : 'owork:',
    'number_different_tetrahedra_line_startwith' : 'NUMBER OF DIFFERENT TETRAHEDRA :',
    'number_k_points_div_afact_line_startswith' : 'NKP,NDIV,afact',
    'len_tretrahedra_point_line' : 85,
}

class OutputkgenReader(object):
    '''
    A class which reads from the .outputkgen file which contains symmetry vectors

    symmetry_matrices:         A set of symmetry matrices that 
    '''

    def __init__(self, filename):
        self.filename = filename
        self.reciprocal_lattice_vectors = None
        self.reciprocal_lattice_vectors_by_2pi = None
        self.reciprocal_lattice_vectors_in_inv_angs = None
        self.symmetry_matrices = []
        self.point_group_symmetry_matrices = []
        self.bloch_vectors = None
        self.tetrahedra = None
        self._load_values()

    def _load_values(self):
        self._file_handle = open(self.filename, 'r')
        self._line_num = 0
        for line in self._file_handle:
            self._line_num = self._line_num + 1
            line = line.strip()

            if self._line_num == fmt['ortho_line_num']:
                #TODO
                pass

            elif self._line_num == fmt['r1_vector_line_num']:
                #TODO
                line2 = self._file_handle.next()
                line3 = self._file_handle.next()
                self._line_num = self._line_num + 2

            elif self._line_num == fmt['iarb_line_num']:
                #TODO
                pass

            elif self._line_num == fmt['point_group_symmetries_line_num']:
                #TODO read in number of point group symmetry matrices
                self._file_handle.next() # Throw away a buffer line
                self.point_group_symmetry_matrices.extend(self._read_symmetry_matrices())
                pass

            elif line.startswith(fmt['symmetry_matrix_header_line_startswith']):
                self.symmetry_matrices.extend(self._read_symmetry_matrices())

            elif line.startswith(fmt['reciprocal_lattice_vectors_header_line_startswith']):
                tmp_vectors = np.zeros((3,3))
                for i in [0,1,2]:
                    self._line_num = self._line_num + 1
                    try:
                        line = self._file_handle.next()
                    except StopIteration:
                        raise UnexpectedFileFormat('Less than 3 lines were specified for the reciprocal lattice vectors (line: %d)' % self._line_num)
                    try:
                        tmp_vectors[:,i] = [float(x.strip()) for x in line.split(' ') if x.strip() != '']
                    except:
                        raise UnexpectedFileFormat('Could not parse reciprocal lattice vectors from line (line: %d)' % self._line_num)
                # First set of reciprocal lattic vetors is followed by an
                # identical set multiplied by 2*pi
                if self.reciprocal_lattice_vectors == None:
                    self.reciprocal_lattice_vectors = tmp_vectors.copy()
                    self.reciprocal_lattice_vectors_in_inv_angs = self.reciprocal_lattice_vectors * (2. * np.pi / CNST.BOHR_RADIUS_IN_ANGSTROM)
                else:
                    self.reciprocal_lattice_vectors_by_2pi = tmp_vectors.copy()

            elif line.startswith(fmt['length_reciprocal_lattice_vectors_line_startswith']):
                #TODO
                pass

            elif line.startswith(fmt['submesh_shift_line_startswith']):
                #TODO
                pass

            elif line.startswith(fmt['number_mesh_points_line_startswith']):
                #TODO
                pass

            elif line.startswith(fmt['reciprocal_lattice_vector_intervals_line_startswith']):
                #TODO
                pass

            elif line.startswith(fmt['num_inequiv_k_points_line_startswith']):
                #TODO
                # fmt['len_bloch_vector_line']
                pass

            elif line.startswith(fmt['tetrahedra_to_sort_line_startswith']):
                #TODO
                pass

            elif line.startswith(fmt['owork_line_startswith']):
                #TODO
                pass

            elif line.startswith(fmt['number_different_tetrahedra_line_startwith']):
                #TODO
                pass

            elif line.startswith(fmt['number_k_points_div_afact_line_startswith']):
                #TODO
                #'len_tretrahedra_point_line' : 85,
                pass

        self._file_handle.close()

    def _read_symmetry_matrices(self):
        tmp_matrices = []
        # Read the matrices off the next three lines into a buffer
        for i in [1,2,3]:
            try:
                line = self._file_handle.next()
                self._line_num = self._line_num + 1
            except StopIteration:
                raise UnexpectedFileFormat( \
                    'Symmetry matrices do not span 3 lines (line: %d)' % self._line_num)
            vals = [int(x.strip()) for x in line.split(' ') if x.strip() != '']
            # Count how many matrices required (may not be three)
            num_matrices = len(vals) / 3
            for j in range(num_matrices):
                tmp_matrices.append([])
                tmp_matrices[j].append(vals[3*j:3*(j+1)])
        # Read the buffered matrices into the object as a Numpy array object only
        # if the determinant is not zero
        returned_matrices = []
        for i in range(num_matrices):
            tmp_matrix = np.array(tmp_matrices[i])
            if np.linalg.det(tmp_matrix) == 0:
                continue
            else:
                returned_matrices.append(SymmetryMatrix(matrix = tmp_matrix))
        return returned_matrices
