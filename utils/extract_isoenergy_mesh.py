'''
extract_isoenergy_mesh.py

Pass in a Kmesh object and an energy

Returns a Nx3 array of i, j, k values corresponding to the closest k points the
the iso-energy surface. Also returns the mean energy for the points and the standard deviation

No interpolation at present

RETURNS:
  (<Nx3 array i,j,k values>, mean_energy, stadard_deviation_energy)
'''

__all__ = ['extract_isoenergy_mesh']

import numpy as np

def extract_isoenergy_mesh(kmesh, energy):
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
    a_i = kmesh.i_series_offset + 0.5*kmesh.i_series_spacing
    a_j = kmesh.j_series_offset + 0.5*kmesh.j_series_spacing
    a_k = kmesh.k_series_offset + 0.5*kmesh.k_series_spacing
    b_i = kmesh.i_series_spacing
    b_j = kmesh.j_series_spacing
    b_k = kmesh.k_series_spacing
    standard_deviation_energy = kmesh.mesh.ravel().std()
    mean_energy = kmesh.mesh.ravel().mean()
    i_indexes, j_indexes, k_indexes = np.where(kmesh.mesh.mask == False)
    del(kmesh)
    i_vals = a_i + i_indexes * b_i
    j_vals = a_j + j_indexes * b_j
    k_vals = a_k + k_indexes * b_k
    return (np.column_stack((i_vals, j_vals, k_vals)), \
            mean_energy, standard_deviation_energy)

# TODO: Implement marching cubes/Shepard interpolation
##     marching_cube_points = [
##             [], # 0
##             [()], # 1
##             [()], # 2
##             [()], # 4
##             [()], # 8
##             [()], # 16
##             [()], # 32
##             [()], # 64
##             [()], # 128
##     ]

if __name__ == '__main__':
    import wien2k.Kmesh as Kmesh

