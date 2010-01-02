'''
SymmetryMatrix.py

An object to store symmetry matrix data
'''

__all__ = ['SymmetryMatrix']

import numpy as np

class SymmetryMatrix(object):
    '''
    An object to store symmetry matrix data
    '''
    def __init__(self, matrix=None, tau_offsets=None, id=None):
        if matrix == None:
            matrix = np.zeros((3,3)) 
        if tau_offsets == None:
            tau_offsets = np.zeros((3)) 
        self.matrix = matrix
        self.tau_offsets = tau_offsets
        self.id = id

