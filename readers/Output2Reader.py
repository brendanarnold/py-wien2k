'''
Output2Reader.py

A class which reads from the .output2 file which contains band lmits
'''

__all__ = ['Output2Reader']

from wien2k.errors import UnexpectedFileFormat

fmt = {
    # The line preceding the band ranges
    'band_range_header' : 'Bandranges (emin - emax):',
    # The start of the line containing the band energy limits
    'band_line_prefix' : 'band',
}

class Output2Reader(object):
    '''A class which reads from the .output2 file which contains band lmits'''
    def __init__(self, filename):
        self.filename = filename
        self.band_limits = []
        flag = {
            'in_band_ranges' : False
        }
        file_handle = open(self.filename, 'r')
        line_num = 0
        for line in file_handle:
            line_num = line_num + 1
            if line.strip().startswith(fmt['band_range_header']):
                flag['in_band_ranges'] = True
                continue
            if flag['in_band_ranges'] == True:
                if line.strip().startswith(fmt['band_line_prefix']):
                    try:
                        vals = [val.strip() for val in line.split(fmt['band_line_prefix'])[1].split(' ') if val.strip() != '']
                        band_num = int(vals[0])
                        band_min = float(vals[1])
                        band_max = float(vals[2])
                    except:
                        raise UnexpectedFileFormat('Could not parse a band energy limits line from file (line: %d)' % line_num)
                    # Make sure array is large enough to hold band number
                    while len(self.band_limits) < band_num:
                        self.band_limits.append(())
                    self.band_limits[band_num-1] = (band_num, band_min, band_max)
                else:
                    flag['in_band_ranges'] == False
        # TODO: Parse the rest
        file_handle.close()