'''
extract_isoenergy_mesh.py

Pass in a Kmesh object, an energy and an optional precision

Returns a Nx3 array of i, j, k values corresponding to the closest k points the
the iso-energy surface. Also returns the mean energy for the points and the standard deviation

Interpolates using linear approximation
TODO: Interpolates using radial basis functions (See Scipy) and bisection

RETURNS:
  <Nx3 array i,j,k values>
'''

__all__ = ['extract_isoenergy_mesh']

import numpy as np
from scipy.interpolate import Rbf
import sys
import copy
import pdb

def extract_isoenergy_mesh(orig_kmesh, energy, precision=sys.float_info.epsilon, verbose=False, interp_method='linear'):
    '''Returns a list of i, j, k values that map a surface'''
##  pdb.set_trace()
    kmesh = copy.deepcopy(orig_kmesh)
    # Builds a 'marching cube' (http://en.wikipedia.org/wiki/Marching_cubes)
    # map of the Kmesh based on whether the values lie under or over the
    # specified energy level. This can be used to detect the surface
    # 'edge'
    isoenergy_3d_mesh = kmesh.mesh > energy
    if verbose == True:
        print 'Building marching cube indexes ...'
    marching_cube_indexes = \
          np.ma.array(isoenergy_3d_mesh[:-1,:-1,:-1], dtype=int) * 1 \
        + np.ma.array(isoenergy_3d_mesh[:-1,:-1:,1:], dtype=int) * 2 \
        + np.ma.array(isoenergy_3d_mesh[:-1,1:,:-1], dtype=int) * 4  \
        + np.ma.array(isoenergy_3d_mesh[:-1,1:,1:], dtype=int) * 8   \
        + np.ma.array(isoenergy_3d_mesh[1:,:-1,:-1], dtype=int) * 16 \
        + np.ma.array(isoenergy_3d_mesh[1:,:-1,:1:], dtype=int) * 32 \
        + np.ma.array(isoenergy_3d_mesh[1:,1:,:-1], dtype=int) * 64  \
        + np.ma.array(isoenergy_3d_mesh[1:,1:,1:], dtype=int) * 128
    del(isoenergy_3d_mesh)
    # == Interpolate with requested method ==
    # Now we know at which points the surface is at, we can use an
    # interpolation method and (if necessary) a root finding algorithm
    # to find where the energy lies spacially
    if interp_method == 'rbf':
        surface_i_vals, surface_j_vals, surface_k_vals = _interp_rbf(kmesh, \
          marching_cube_indexes, precision, verbose, energy)
    elif interp_method == 'linear':
        surface_i_vals, surface_j_vals, surface_k_vals = _interp_linear(kmesh, \
          marching_cube_indexes, verbose, energy)
    elif interp_method == 'nearest':
        surface_i_vals, surface_j_vals, surface_k_vals = _interp_nearest(kmesh, \
          marching_cube_indexes, verbose)
    # We have the positions of the surface in terms of location within
    # the mesh, we need to convert the indexes to real values.
    if verbose == True:
        print 'Converting back to real co-ordinates ...'
    a_i = kmesh.i_series_offset
    a_j = kmesh.j_series_offset
    a_k = kmesh.k_series_offset
    b_i = kmesh.i_series_spacing
    b_j = kmesh.j_series_spacing
    b_k = kmesh.k_series_spacing
    del(kmesh)
    i_vals = a_i + surface_i_vals * b_i
    j_vals = a_j + surface_j_vals * b_j
    k_vals = a_k + surface_k_vals * b_k
    return (np.column_stack((i_vals, j_vals, k_vals)))


def _interp_nearest(kmesh, marching_cube_indexes, verbose):
    '''
    Gets the points that lie on the edge -not exactly nearest neighbour
    '''
    surface_i_vals, surface_j_vals, surface_k_vals = np.where( \
        (marching_cube_indexes != 255) & \
        (marching_cube_indexes != 0) \
    )
    return (surface_i_vals, surface_j_vals, surface_k_vals)


def _interp_linear(kmesh, marching_cube_indexes, verbose, energy):
    '''
    Use a 1D linear approximation to obtain better values for the surface
    '''
    def _linear_equation(x1, x2, y1, y2, y3):
        '''
        Given two points (x1, y1) and (x2, y2) and a third y3,
        gives the corresponding x3
        '''
        x3 = (y3-y2)*(x2-x1)/((y2-y1).astype(float)) + x2
        return x3
    if verbose == True:
        print 'Using linear interpolation ...'
    surface_i_vals = np.array([])
    surface_j_vals = np.array([])
    surface_k_vals = np.array([])
    if verbose == True:
        print 'Interpolating along i direction ...'
    i_ind, j_ind, k_ind = np.where( \
      (((marching_cube_indexes & 1) + (marching_cube_indexes & 16)) == 1) | \
      (((marching_cube_indexes & 1) + (marching_cube_indexes & 16)) == 16) \
    )
    new_i_vals = _linear_equation(i_ind, i_ind+1, \
      kmesh.mesh[i_ind, j_ind, k_ind], kmesh.mesh[i_ind+1, j_ind, k_ind], energy)
    surface_i_vals = np.append(surface_i_vals, new_i_vals)
    surface_j_vals = np.append(surface_j_vals, j_ind)
    surface_k_vals = np.append(surface_k_vals, k_ind)
    if verbose == True:
        print 'Interpolating along j direction ...'
    i_ind, j_ind, k_ind = np.where( \
      (((marching_cube_indexes & 1) + (marching_cube_indexes & 4)) == 1) | \
      (((marching_cube_indexes & 1) + (marching_cube_indexes & 4)) == 4) \
    )
    new_j_vals = _linear_equation(j_ind, j_ind+1, \
      kmesh.mesh[i_ind, j_ind, k_ind], kmesh.mesh[i_ind, j_ind+1, k_ind], energy)
    surface_i_vals = np.append(surface_i_vals, i_ind)
    surface_j_vals = np.append(surface_j_vals, new_j_vals)
    surface_k_vals = np.append(surface_k_vals, k_ind)
    if verbose == True:
        print 'Interpolating along k direction ...'
    i_ind, j_ind, k_ind = np.where( \
      (((marching_cube_indexes & 1) + (marching_cube_indexes & 2)) == 1) | \
      (((marching_cube_indexes & 1) + (marching_cube_indexes & 2)) == 2) \
    )
    new_k_vals = _linear_equation(k_ind, k_ind+1, \
      kmesh.mesh[i_ind, j_ind, k_ind], kmesh.mesh[i_ind, j_ind, k_ind+1], energy)
    surface_i_vals = np.append(surface_i_vals, i_ind)
    surface_j_vals = np.append(surface_j_vals, j_ind)
    surface_k_vals = np.append(surface_k_vals, new_k_vals)
    return (surface_i_vals, surface_j_vals, surface_k_vals)
    

def _interp_rbf(kmesh, marching_cube_indexes, precision, verbose, energy):
    '''
    Use the radial basis function from Scipy to obtain better values for
    the surface
    TODO: Find out exactly what radial basis functions actually does ...
    TODO: Get this to work on smaller numbers of points with a moving window
    '''
    if verbose == True:
        print 'Generating Radial Basis Function for interpolation ...'
    i_indexes, j_indexes, k_indexes = np.where(kmesh.mesh.mask == False)
    rbf_energy_function = Rbf(i_indexes, j_indexes, k_indexes, kmesh.mesh[i_indexes, j_indexes, k_indexes])
    kmesh.mesh.mask = kmesh.mesh.mask | (marching_cube_indexes == 0) | (marching_cube_indexes == 255)
    surface_i_vals = np.array([])
    surface_j_vals = np.array([])
    surface_k_vals = np.array([])
    # Pick out points where corner 1 is inside surface but corner 16 is not,
    # then bisect 1->16
    if verbose == True:
        print 'Interpolating along i direction ...'
    tmp_i_vals, tmp_j_vals, tmp_k_vals = np.where(((marching_cube_indexes & 1) \
            + (marching_cube_indexes & 16)) == 1)
    tmp_i_vals, tmp_j_vals, tmp_k_vals = _bisect_along_line(rbf_energy_function, \
            tmp_i_vals, tmp_j_vals, tmp_k_vals, direction='i', reverse=False, \
            energy=energy, precision=precision, verbose=verbose)
    surface_i_vals = np.append(surface_i_vals, tmp_i_vals)
    surface_j_vals = np.append(surface_j_vals, tmp_j_vals)
    surface_k_vals = np.append(surface_k_vals, tmp_k_vals)
    # Pick out points where corner 16 is inside surface but corner 1 is not,
    # then bisect 16->1
    tmp_i_vals, tmp_j_vals, tmp_k_vals = np.where(((marching_cube_indexes & 1) \
            + (marching_cube_indexes & 16)) == 16)
    tmp_i_vals, tmp_j_vals, tmp_k_vals = _bisect_along_line(rbf_energy_function, \
            tmp_i_vals, tmp_j_vals, tmp_k_vals, direction='i', reverse=True, \
            energy=energy, precision=precision, verbose=verbose)
    surface_i_vals = np.append(surface_i_vals, tmp_i_vals)
    surface_j_vals = np.append(surface_j_vals, tmp_j_vals)
    surface_k_vals = np.append(surface_k_vals, tmp_k_vals)
    # Pick out points where corner 1 is inside surface but corner 4 is not,
    # then bisect 1->4
    if verbose == True:
        print 'Interpolating along j direction ...'
    tmp_i_vals, tmp_j_vals, tmp_k_vals = np.where(((marching_cube_indexes & 1) \
            + (marching_cube_indexes & 4)) == 1)
    tmp_i_vals, tmp_j_vals, tmp_k_vals = _bisect_along_line(rbf_energy_function, \
            tmp_i_vals, tmp_j_vals, tmp_k_vals, direction='j', reverse=False, \
            energy=energy, precision=precision, verbose=verbose)
    surface_i_vals = np.append(surface_i_vals, tmp_i_vals)
    surface_j_vals = np.append(surface_j_vals, tmp_j_vals)
    surface_k_vals = np.append(surface_k_vals, tmp_k_vals)
    # Pick out points where corner 4 is inside surface but corner 1 is not,
    # then bisect 4->1
    tmp_i_vals, tmp_j_vals, tmp_k_vals = np.where(((marching_cube_indexes & 1) \
            + (marching_cube_indexes & 4)) == 4)
    tmp_i_vals, tmp_j_vals, tmp_k_vals = _bisect_along_line(rbf_energy_function, \
            tmp_i_vals, tmp_j_vals, tmp_k_vals, direction='j', reverse=True, \
            energy=energy, precision=precision, verbose=verbose)
    surface_i_vals = np.append(surface_i_vals, tmp_i_vals)
    surface_j_vals = np.append(surface_j_vals, tmp_j_vals)
    surface_k_vals = np.append(surface_k_vals, tmp_k_vals)
    # Pick out points where corner 1 is inside surface but corner 2 is not,
    # then bisect 1->2
    if verbose == True:
        print 'Interpolating along k direction ...'
    tmp_i_vals, tmp_j_vals, tmp_k_vals = np.where(((marching_cube_indexes & 1) \
            + (marching_cube_indexes & 2)) == 1)
    tmp_i_vals, tmp_j_vals, tmp_k_vals = _bisect_along_line(rbf_energy_function, \
            tmp_i_vals, tmp_j_vals, tmp_k_vals, direction='k', reverse=False, \
            energy=energy, precision=precision, verbose=verbose)
    surface_i_vals = np.append(surface_i_vals, tmp_i_vals)
    surface_j_vals = np.append(surface_j_vals, tmp_j_vals)
    surface_k_vals = np.append(surface_k_vals, tmp_k_vals)
    # Pick out points where corner 2 is inside surface but corner 1 is not,
    # then bisect 2->1
    tmp_i_vals, tmp_j_vals, tmp_k_vals = np.where(((marching_cube_indexes & 1) \
            + (marching_cube_indexes & 2)) == 2)
    tmp_i_vals, tmp_j_vals, tmp_k_vals = _bisect_along_line(rbf_energy_function, \
            tmp_i_vals, tmp_j_vals, tmp_k_vals, direction='k', reverse=True, \
            energy=energy, precision=precision, verbose=verbose)
    surface_i_vals = np.append(surface_i_vals, tmp_i_vals)
    surface_j_vals = np.append(surface_j_vals, tmp_j_vals)
    surface_k_vals = np.append(surface_k_vals, tmp_k_vals)
    return (surface_i_vals, surface_j_vals, surface_k_vals)
    
def _bisect_along_line(fn, i_vals, j_vals, k_vals, direction='i', reverse=False, energy=0.0, precision=sys.float_info.epsilon, verbose=False):
    '''This helper routine adjusts co-ordinates along a certain direction using
    bisection until the lie within precision*2 of the value'''
    # Each line is unit distance so begin with half unit step
    delta = 0.5
    if direction == 'i':
        adjusted_i_vals = np.ma.array(i_vals)
        increments = np.ma.zeros(adjusted_i_vals.shape)
    elif direction == 'j':
        adjusted_j_vals = np.ma.array(j_vals)
        increments = np.ma.zeros(adjusted_j_vals.shape)
    elif direction == 'k':
        adjusted_k_vals = np.ma.array(k_vals)
        increments = np.ma.zeros(adjusted_k_vals.shape)
    while delta > precision:
        if verbose == True:
            print 'Interpolated within %e' % delta
        if direction == 'i':
            energies = fn(adjusted_i_vals + delta, j_vals, k_vals)
        elif direction == 'j':
            energies = fn(i_vals, adjusted_j_vals + delta, k_vals)
        elif direction == 'k':
            energies = fn(i_vals, j_vals, adjusted_k_vals + delta)
        if reverse == False:
            increments.mask = energies <= energy
        else:
            increments.mask = energies >= energy
        increments = increments + delta
        increments.mask = False
        delta = delta/2.0
    if direction == 'i':
        i_vals = adjusted_i_vals.data
    elif direction == 'j':
        j_vals = adjusted_j_vals.data
    elif direction == 'k':
        k_vals = adjusted_k_vals.data
    return (i_vals, j_vals, k_vals)


if __name__ == '__main__':
    pass

