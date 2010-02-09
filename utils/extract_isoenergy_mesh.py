'''
extract_isoenergy_mesh.py

Pass in a Kmesh object, an energy and an optional precision

Returns a Nx3 array of i, j, k values corresponding to the closest k points the
the iso-energy surface. Also returns the mean energy for the points and the standard deviation

Interpolates using radial basis functions (See Scipy) and bisection

RETURNS:
  <Nx3 array i,j,k values>
'''

__all__ = ['extract_isoenergy_mesh']

import numpy as np
from scipy.interpolate import Rbf
import sys

def extract_isoenergy_mesh(kmesh, energy, precision=sys.float_info.epsilon):
    '''Returns a list of i, j, k values that map a surface'''
    isoenergy_3d_mesh = kmesh.mesh > energy
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
    kmesh.mesh = kmesh.mesh[:-1,:-1,:-1]
    kmesh.mesh.mask = marching_cube_indexes | (marching_cube_indexes == 0) | (marching_cube_indexes == 255)
    i_indexes, j_indexes, k_indexes = np.where(kmesh.mesh.mask == False)

    # == Interpolate and bisect to find energy surface ==
    # Use the radial basis function from Scipy to obtain better values
    # TODO: Find out exactly what radial basis functions actually does ...
    rbf_energy_function = Rbf(i_indexes, j_indexes, k_indexes, kmesh.kmesh[i_indexes, j_indexes, k_indexes])
    surface_i_vals = np.array([])
    surface_j_vals = np.array([])
    surface_k_vals = np.array([])
    # Pick out points where corner 1 is inside surface but corner 16 is not,
    # then bisect 1->16
    tmp_i_vals, tmp_j_vals, tmp_k_vals = np.where(((marching_cube_indexes & 1) \
            + (marching_cube_indexes & 16)) == 1)
    tmp_i_vals, tmp_j_vals, tmp_k_vals = _bisect_along_line(rbf_energy_function, \
            tmp_i_vals, tmp_j_vals, tmp_k_vals, direction='i', reverse=False, \
            energy, precision)
    surface_i_vals = np.append(surface_i_vals, tmp_i_vals)
    surface_j_vals = np.append(surface_j_vals, tmp_j_vals)
    surface_k_vals = np.append(surface_k_vals, tmp_k_vals)
    # Pick out points where corner 16 is inside surface but corner 1 is not,
    # then bisect 16->1
    tmp_i_vals, tmp_j_vals, tmp_k_vals = np.where(((marching_cube_indexes & 1) \
            + (marching_cube_indexes & 16)) == 16)
    tmp_i_vals, tmp_j_vals, tmp_k_vals = _bisect_along_line(rbf_energy_function, \
            tmp_i_vals, tmp_j_vals, tmp_k_vals, direction='i', reverse=True, \
            energy, precision)
    surface_i_vals = np.append(surface_i_vals, tmp_i_vals)
    surface_j_vals = np.append(surface_j_vals, tmp_j_vals)
    surface_k_vals = np.append(surface_k_vals, tmp_k_vals)
    # Pick out points where corner 1 is inside surface but corner 4 is not,
    # then bisect 1->4
    tmp_i_vals, tmp_j_vals, tmp_k_vals = np.where(((marching_cube_indexes & 1) \
            + (marching_cube_indexes & 4)) == 1)
    tmp_i_vals, tmp_j_vals, tmp_k_vals = _bisect_along_line(rbf_energy_function, \
            tmp_i_vals, tmp_j_vals, tmp_k_vals, direction='j', reverse=False, \
            energy, precision)
    surface_i_vals = np.append(surface_i_vals, tmp_i_vals)
    surface_j_vals = np.append(surface_j_vals, tmp_j_vals)
    surface_k_vals = np.append(surface_k_vals, tmp_k_vals)
    # Pick out points where corner 4 is inside surface but corner 1 is not,
    # then bisect 4->1
    tmp_i_vals, tmp_j_vals, tmp_k_vals = np.where(((marching_cube_indexes & 1) \
            + (marching_cube_indexes & 4)) == 4)
    tmp_i_vals, tmp_j_vals, tmp_k_vals = _bisect_along_line(rbf_energy_function, \
            tmp_i_vals, tmp_j_vals, tmp_k_vals, direction='j', reverse=True, \
            energy, precision)
    surface_i_vals = np.append(surface_i_vals, tmp_i_vals)
    surface_j_vals = np.append(surface_j_vals, tmp_j_vals)
    surface_k_vals = np.append(surface_k_vals, tmp_k_vals)
    # Pick out points where corner 1 is inside surface but corner 2 is not,
    # then bisect 1->2
    tmp_i_vals, tmp_j_vals, tmp_k_vals = np.where(((marching_cube_indexes & 1) \
            + (marching_cube_indexes & 2)) == 1)
    tmp_i_vals, tmp_j_vals, tmp_k_vals = _bisect_along_line(rbf_energy_function, \
            tmp_i_vals, tmp_j_vals, tmp_k_vals, direction='k', reverse=False, \
            energy, precision)
    surface_i_vals = np.append(surface_i_vals, tmp_i_vals)
    surface_j_vals = np.append(surface_j_vals, tmp_j_vals)
    surface_k_vals = np.append(surface_k_vals, tmp_k_vals)
    # Pick out points where corner 2 is inside surface but corner 1 is not,
    # then bisect 2->1
    tmp_i_vals, tmp_j_vals, tmp_k_vals = np.where(((marching_cube_indexes & 1) \
            + (marching_cube_indexes & 2)) == 2)
    tmp_i_vals, tmp_j_vals, tmp_k_vals = _bisect_along_line(rbf_energy_function, \
            tmp_i_vals, tmp_j_vals, tmp_k_vals, direction='k', reverse=True, \
            energy, precision)
    surface_i_vals = np.append(surface_i_vals, tmp_i_vals)
    surface_j_vals = np.append(surface_j_vals, tmp_j_vals)
    surface_k_vals = np.append(surface_k_vals, tmp_k_vals)

    # == Convert to real co-ordinates ==
    a_i = kmesh.i_series_offset + 0.5*kmesh.i_series_spacing
    a_j = kmesh.j_series_offset + 0.5*kmesh.j_series_spacing
    a_k = kmesh.k_series_offset + 0.5*kmesh.k_series_spacing
    b_i = kmesh.i_series_spacing
    b_j = kmesh.j_series_spacing
    b_k = kmesh.k_series_spacing
    del(kmesh)
    i_vals = a_i + surface_i_vals * b_i
    j_vals = a_j + surface_j_vals * b_j
    k_vals = a_k + surface_k_vals * b_k

    return (np.column_stack((i_vals, j_vals, k_vals)))


    
def _bisect_along_line(fn, i_vals, j_vals, k_vals, direction, reverse, energy, precision):
    '''This helper routine adjusts co-ordinates along a certain direction using
    bisection until the lie within precision*2 of the value'''
    # Each line is unit distance so begin with half unit step
    delta = 0.5
    if direction == 'i':
        adjusted_i_vals = np.ma.array(i_vals)
    elif direction == 'j':
        adjusted_j_vals = np.ma.array(j_vals)
    elif direction == 'k':
        adjusted_k_vals = np.ma.array(k_vals)
    while delta > precision:
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
    import wien2k.Kmesh as Kmesh

