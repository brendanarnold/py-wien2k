'''
SymMat.py

An object to store symmetry matrix data
'''

__all__ = ['SymMat']

import numpy as np

class SymMat(object):
    '''
    An object to store symmetry matrix data as well as perform symmetry
    transforms

    Input:
    matrix          The symmetry matrix itself, optional constructor
                    argument, by default will be a zero filled Nump 2D array
    tau_offsets     The vector offset to apply the points to push them
                    off high symmetry points, optional constructor
                    argument, default a zero filled Numpy array of length 3 

    Example:

    Creating an identity symmetry matrix

    >>> sm = SymMat(np.identity(3))
    >>> sm.matrix
    array([[ 1.,  0.,  0.],
           [ 0.,  1.,  0.],
           [ 0.,  0.,  1.]])

    '''
    def __init__(self, matrix=None, tau_offsets=None, id=None):
        if matrix == None:
            matrix = np.zeros((3,3)) 
        if tau_offsets == None:
            tau_offsets = np.zeros((3)) 
        self.matrix = matrix
        self.tau_offsets = tau_offsets
        self.id = id

    def __mul__(self, other):
        '''Allows for easy combining of matrix operations'''
        if type(other) is SymMat:
            return SymMat(np.dot(self.matrix, other.matrix))
        elif (type(other) is numpy.ndarray) and \
          (other.ndims == 2) and (other.shape[1] > 3):
            return self.map(other)
        else:
            raise ValueError('Cannot perform multiplication on this object')

    def map(self, klist, cols=[1,2,3], inverse=False):
        '''
        Apply this symmetry matrix to a list of k-vectors
        
        Input:
        klist       A 2D Numpy array containing at least i,j,k values of
                    k vectors
        cols        A list of the columns containing the i, j, k 
                    values (default [1,2,3] i.e. 2nd, 3rd and 4th columns)
        inverse     Apply the inverse transform

        Ouptut:
        transformed_klist   A similar vector to input but with k vectors
                    tranformed according to the matrix
        i.e.

        Examples:

        Mapping a klist with identity matrix (i.e. does nothing)

        >>> import numpy as np
        >>> kl = np.array([
        ... [1,0,0,1],
        ... [2,0,0,2],
        ... [3,0,1,0],
        ... [4,0,1,1]])
        >>> sm = SymMat(np.identity(3))
        >>> sm.matrix
        array([[ 1.,  0.,  0.],
               [ 0.,  1.,  0.],
               [ 0.,  0.,  1.]])
        >>> sm.map(kl)
        array([[1, 0, 0, 1],
               [2, 0, 0, 2],
               [3, 0, 1, 0],
               [4, 0, 1, 1]])
        
        Now try with a tau offset ...

        >>> sm.tau_offsets = np.array([0,2,0])
        >>> tkl = sm.map(kl)
        >>> tkl
        array([[1, 0, 2, 1],
               [2, 0, 2, 2],
               [3, 0, 3, 0],
               [4, 0, 3, 1]])


        Now do the inverse ...
        >>> sm.map(tkl, inverse=True)
        array([[1, 0, 0, 1],
               [2, 0, 0, 2],
               [3, 0, 1, 0],
               [4, 0, 1, 1]])
        
        '''
        # Check input
        if klist.ndim != 2:
            raise ValueError('klist must be 2D - actual dimensions: %d' % klist.ndim)
        # Return a copy
        transformed_klist = klist.copy()
        # Extract just the list of kpoints
        kpoints = transformed_klist[:,cols]
        # n.b. Normally use <Operator>.<Vector> to transform but for Numpy to
        # iterate over a list of vectors this has to be reversed and hence the
        # operator transposed
        if inverse == True:
            kpoints = np.dot(kpoints, np.linalg.inv(self.matrix).transpose())
            # Remove the tau offset to put the kpoints back on the high
            # symmetry areas
            kpoints[:,0] = kpoints[:,0] - self.tau_offsets[0]
            kpoints[:,1] = kpoints[:,1] - self.tau_offsets[1]
            kpoints[:,2] = kpoints[:,2] - self.tau_offsets[2]
        else:
            # Add the tau offsets to nudge the kpoints off high symmetry
            # areas
            kpoints[:,0] = kpoints[:,0] + self.tau_offsets[0]
            kpoints[:,1] = kpoints[:,1] + self.tau_offsets[1]
            kpoints[:,2] = kpoints[:,2] + self.tau_offsets[2]
            kpoints = np.dot(kpoints, self.matrix.transpose())
        # Copy back to original klist to match ids, energies, whatever
        transformed_klist[:,cols] = kpoints
        return transformed_klist

if __name__ == '__main__':
    import doctest
    doctest.testmod()
