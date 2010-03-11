__all__ = ['remove_duplicate_kpoints']

def remove_duplicate_kpoints(kpoints, 
  i_col=None, j_col=None, k_col=None, tolerance=None):
    '''
    Removes duplicate k points from a list of data points
    '''
    # Deal with the parameters
    if tolerance is not None:
        try:
            i_tol, j_tol, k_tol = iter(tolerance)
        except:
            i_tol = j_tol = k_tol = tolerance
    else:
        tolerance = 0
   # TODO the rest ...
    
