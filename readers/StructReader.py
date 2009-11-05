'''
StructReader.py

Reads a WIEN2k .struct file into an object
'''

__all__ = ['StructReader']

from wien2k.errors import UnexpectedFileFormat

lattice_type_line = 2
units_line = 3
lattice_parameters_line = 4

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
        line_num = 0
        for line in file_handle:
            line_num = line_num + 1
            if line_num == lattice_type_line:
                self.lattice_type = line[0]
                self.space_group = line.split(':')[-1].strip().split('_')[1]
            if line_num == units_line:
                self.units = line.split('=')[-1].strip()
            if line_num == lattice_parameters_line:
                try:
                    self.a, self.b, self.c, self.alpha, self.beta, self.gamma = \
                        [float(num.strip()) for num in line.split(' ') if num.strip() != '']
                except:
                    raise UnexpectedFileFormat('Line idetified as containing lattice parameters could not be parsed (line: %d)' % line_num)
        # TODO: Rest of parsing
        file_handle.close()