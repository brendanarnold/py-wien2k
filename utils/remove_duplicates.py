__all__ = ['remove_duplicates']

import numpy as np

def remove_duplicates(data, dp_tol=None, cols=None, sort_by=None):
    '''
    Removes duplicate vectors from a list of data points
    Parameters:
        data        An MxN array of N vectors of dimension M 
        cols        An iterable of the columns that must match 
                    in order to constitute a duplicate 
                    (default: [1,2,3] for typical Klist data array) 
        dp_tol      An iterable of three tolerances or a single 
                    tolerance for all dimensions. Uses this to round 
                    the values to specified number of decimal places 
                    before performing the removal. 
                    (default: None)
        sort_by     An iterable of columns to sort by (default: [0])
                        
    Returns:
        MxI Array   An array of I vectors (minus the 
                    duplicates)

    EXAMPLES:
    
    Remove a duplicate

    >>> import wien2k.utils
    >>> import numpy as np
    >>> vecs1 = np.array([[1, 0, 0, 0],
    ...     [2, 0, 0, 0],
    ...     [3, 0, 0, 1]])
    >>> remove_duplicates(vecs1)
    array([[1, 0, 0, 0],
           [3, 0, 0, 1]])

    Remove duplicates with a tolerance

    >>> vecs2 = np.array([[1, 0, 0, 0  ],
    ...     [2, 0, 0, 0.001 ],
    ...     [3, 0, 0, 0.02  ],
    ...     [4, 0, 0, 1     ]])
    >>> remove_duplicates(vecs2, dp_tol=2)
    array([[ 1.  ,  0.  ,  0.  ,  0.  ],
           [ 3.  ,  0.  ,  0.  ,  0.02],
           [ 4.  ,  0.  ,  0.  ,  1.  ]])

    Remove duplicates and sort by k values

    >>> vecs3 = np.array([[1, 0, 0, 0],
    ...     [2, 0, 0, 2],
    ...     [3, 0, 0, 0],
    ...     [4, 0, 0, 1]])
    >>> remove_duplicates(vecs3, sort_by=[3])
    array([[1, 0, 0, 0],
           [4, 0, 0, 1],
           [2, 0, 0, 2]])

    Change the columns that constitute a duplicate

    >>> vecs4 = np.array([[1, 0, 0, 0],
    ...     [2, 0, 0, 2],
    ...     [1, 0, 0, 0],
    ...     [4, 0, 0, 1]])
    >>> remove_duplicates(vecs4, cols=[0])
    array([[1, 0, 0, 0],
           [2, 0, 0, 2],
           [4, 0, 0, 1]])

    '''
    # Deal with the parameters
    if sort_by is None:
        sort_by = [0]
    if cols is None:
        cols = [1,2,3]
    if dp_tol is not None:
        # test to see if already an iterable
        try:
            null = iter(dp_tol)
            tols = np.array(dp_tol)
        except TypeError:
            tols = np.ones_like(cols) * dp_tol
        # Convert to numbers of decimal places
        # Find the 'order' of the axes
    else:
        tols = None

    rnd_data = data.copy()
    # set the tolerances
    if tols is not None:
        for col,tol in zip(cols, tols):
            rnd_data[:,col] = np.around(rnd_data[:,col], decimals=tol)

    # TODO: For now, use a slow Python 'for' loop, try to find a more
    # numponic way later - see: http://stackoverflow.com/questions/2433882/
    sorted_indexes = np.lexsort(tuple([rnd_data[:,col] for col in cols]))
    rnd_data = rnd_data[sorted_indexes]
    unique_kpts = []
    for i in xrange(len(rnd_data)):
        if i == 0:
            unique_kpts.append(i)    
        else:
            if (rnd_data[i, cols] == rnd_data[i-1, cols]).all():
                continue
            else:
                unique_kpts.append(i)    
    
    rnd_data =  rnd_data[unique_kpts]
    # Now sort
    sorted_indexes = np.lexsort(tuple([rnd_data[:,col] for col in sort_by]))
    rnd_data = rnd_data[sorted_indexes]
    return rnd_data
    


if __name__ == '__main__':
    import doctest
    from os.path import join, abspath
    doctest.testmod()
    doctest.testfile(abspath(join('..', 'tests', 'remove_duplicates_test.txt')))
