'''Band.py

A class representing an energy band
'''

class Band(object):
    '''An object representing the bands in the calculation
    '''
    def __init__(self, name='', character=''):
        self.name = name
        self.character = character
        self.data = None
        
    def kpoints(self):
        if self.data is not None:
            return self.data[:,:3] # This slice notation is supremely gay - means columns 1-3
        else:
            return None
        
    def i_vals(self):
        if self.data is not None:
            return self.data[:,0]
        else:
            return None
    def j_vals(self):
        if self.data is not None:
            return self.data[:,1]
        else:
            return None
        
    def k_vals(self):
        if self.data is not None:
            return self.data[:,2]
        else:
            return None

    def energies(self):
        if self.data is not None:
            return self.data[:,3] # This means column 4
        else:
            return None

    # Overwrite the function labels with properties
    kpoints = property(kpoints)
    i_vals = property(i_vals)
    j_vals = property(j_vals)
    k_vals = property(k_vals)
    energies = property(energies)
    