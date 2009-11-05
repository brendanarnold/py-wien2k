'''
StructReader.py

Reads a WIEN2k .struct file into an object
'''

__all__ = ['StructReader']

class StructReader(object):
    '''Reads a WIEN2k .struct file specified in 'filename' parameter into an object

    At present, properties read include,
        lattice_type:   Single letter lattice type specified
        space_group:    Textual space group (i.e. 'I4/mmm')
        units:          Units (i.e. 'bohr')
        a:              Lattice a parameter
        b:              Lattice b parameter
        c:              Lattice c parameter
        alpha:          Lattice alpha angle
        beta:           Lattice beta angle
        gamma:          Lattice gamma angle
    '''
    def __init__(self, filename):
        self.filename = filename
        # Populate values from the struct file
        file_handle = open(self.filename, 'r')
        line_num = 1
        for line in file_handle:
            if line_num == 2:
                self.lattice_type = line[0]
                self.space_group = line.split(':')[-1].strip().split('_')[1]
            if line_num == 3:
                self.units = line.split('=')[-1].strip()
            if line_num == 4:
                self.a, self.b, self.c, self.alpha, self.beta, self.gamma = \
                    [float(num.strip()) for num in line.split(' ') if num.strip() != '']
            line_num = line_num + 1
        # TODO: Rest of parsing
        file_handle.close()