'''Band.py

A class representing an energy band
'''

class Band(object):
    '''An object representing the bands in the calculation
    '''
    def __init__(self, name='', character=''):
        self.name = name
        self.character = character
        self.kpoints = None
        self.energies = None