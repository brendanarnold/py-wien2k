__all__ = ['remove_duplicates']

import numpy as np

def remove_duplicates(data, cols=[1,2,3], tol=None):
    '''
    Removes duplicate vectors from a list of data points
    Parameters:
        data        An MxN array of N vectors of dimension M 
        cols        An iterable of the columns that must match 
                    in order to constitute a duplicate 
                    (default: 1,2,3 for typical Klist data array) 
        tol         An iterable of three tolerances or a single 
                    tolerance for each dimension. Uses this to round 
                    the values before performing the removal. Note 
                    these values are rounded to the nearest negative
                    power of 10 as they actually represent decimal
                    places to round to
                    (default: None)
                        
    Returns:
        MxI Array   An array of I vectors (minus the 
                    duplicates)

    EXAMPLE:
    
    Remove a duplicate

    >>> import wien2k.utils
    >>> import numpy as np
    >>> vecs = np.array([[1, 0, 0, 0],
    ...     [2, 0, 0, 0],
    ...     [3, 0, 0, 1]])
    >>> remove_duplicates(vecs)
    array([[1, 0, 0, 0],
           [3, 0, 0, 1]])

    Remove a duplicate with a tolerance

    >>> import wien2k.utils
    >>> import numpy as np
    >>> vecs = np.array([[1, 0, 0, 0],
    ...     [2, 0, 0, 0.0000001],
    ...     [3, 0, 0, 1]])
    >>> remove_duplicates(vecs, tol=1e-2)
    array([[ 1.,  0.,  0.,  0.],
           [ 3.,  0.,  0.,  1.]])


    '''
    # Deal with the parameters
    if tol is not None:
        # test to see if already an iterable
        try:
            null = iter(tol)
            tol = np.array(tol)
        except TypeError:
            tol = np.ones_like(cols) * tol
        # Convert to numbers of decimal places
        tol = np.around((np.log10(10*np.ones_like(tol)/tol))).astype(int)

    rnd_data = data.copy()
    # set the tolerances
    if tol is not None:
        for i,col in enumerate(cols):
            rnd_data[:,col] = np.round(rnd_data[:,col], tol[i])

    # TODO: For now, use a slow Python 'for' loop, try to find a more
    # numponic way later - see: http://stackoverflow.com/questions/2433882/
    cols = np.array(cols)
    sorted_indexes = np.lexsort(tuple([data[:,col] for col in cols]))
    rnd_data = rnd_data[sorted_indexes]
    unique_kpts = []
    for i in xrange(len(rnd_data)):
        if i == 0:
            unique_kpts.append(sorted_indexes[i])    
        else:
            if (rnd_data[i, cols] == rnd_data[i-1, cols]).all():
                continue
            else:
                unique_kpts.append(sorted_indexes[i])    
    
    return data[unique_kpts,:]


if __name__ == '__main__':
    import doctest
    from os.path import join, abspath
    doctest.testmod()
    doctest.testfile(abspath(join('..', 'tests', 'remove_duplicates_test.txt')))
