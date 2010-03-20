__all__ = ['remove_duplicates']

import numpy as np

def remove_duplicates(data, cols=(1,2,3), tolerance=None):
    '''
    Removes duplicate vectors from a list of data points
    Parameters:
        data        An MxN array of N vectors of dimension M 
        cols        An iterable of the columns that must match 
                    in order to constitute a duplicate 
                    (default: 1,2,3 for typical Klist data array) 
        tolerance   An iterable of three tolerances or a single 
                    tolerance for each dimension (default: 0)
    Returns:
        MxI Array   An array of I vectors (minus the 
                    duplicates)
    '''
    # Deal with the parameters
##     if tolerance is not None:
##         try:
##             i_tol, j_tol, k_tol = iter(tolerance)
##         except:
##             i_tol = j_tol = k_tol = tolerance
##     else:
##         tolerance = 0
    # TODO: For now, use a slow Python 'for' loop, try to find a more
    # numponic way later - see: http://stackoverflow.com/questions/2433882/
    cols = np.array(cols)
    sorted_indexes = np.lexsort(tuple([data[:,col] for col in cols]))
    data = data[sorted_indexes]
    unique_kpts = []
    for i in xrange(len(data)):
        if i == 0:
            unique_kpts.append(data[i,:].tolist())    
        else:
            if (data[i, cols] == data[i-1, cols]).all():
                continue
            else:
                unique_kpts.append(data[i,:].tolist())    
    return np.array(unique_kpts)
