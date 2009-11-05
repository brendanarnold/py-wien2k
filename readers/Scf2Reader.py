'''Scf2Reader.py

Reads a WIEN2k .scf2 file into an object
'''

__all__ = ['Scf2Reader']

from wien2k.errors import UnexpectedFileFormat

fermi_energy_line_startswith = ':FER'

class Scf2Reader(object):
    '''Reads a WIEN2k .struct file specified in 'filename' parameter into an object

    At present, properties read include,
        fermi_energy: Fermi energy for this calculation
    '''
    def __init__(self, filename):
        self.filename = filename
        self.fermi_energy = None
        file_handle = open(self.filename, 'r')
        line_num = 0
        for line in file_handle:
            line_num = line_num + 1
            # Takes the last fermi energy (i.e. from last interation) in the file
            if line.strip().startswith(fermi_energy_line_startswith):
                try:
                    self.fermi_energy = float(line.split('=')[-1].strip())
                except:
                    raise UnexpectedFileFormat('Could not parse a Fermi energy from a line identified as containing one (line: %d)' % line_num)
        # TODO: Include other stuff
        file_handle.close()