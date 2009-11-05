'''Scf2Reader.py

Reads a WIEN2k .scf2 file into an object
'''

__all__ = ['Scf2Reader']

class Scf2Reader(object):
    '''Reads a WIEN2k .struct file specified in 'filename' parameter into an object

    At present, properties read include,
        fermi_energy: Fermi energy for this calculation
    '''
    def __init__(self, filename):
        self.filename = filename
        self.fermi_energy = None
        file_handle = open(self.filename, 'r')
        for line in file_handle:
            # Takes the last fermi energy (i.e. from last interation) in the file
            if line.strip().startswith(':FER'):
                self.fermi_energy = float(line.split('=')[-1].strip())
        # TODO: Include other stuff
        file_handle.close()