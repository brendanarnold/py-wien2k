
__all__ = ['expand_irreducible_brillouin_zone']

import numpy as np


def expand_irreducible_brillouin_zone(symmetry_matrices, k_points):
    tot_k_points = []
    for k_point in k_points:
        tmp_k_points = []
        for symmetry_matrix in symmetry_matrices:
            symmetry_point = list(np.dot(k_point[1:4], symmetry_matrix.matrix) \
                    + symmetry_matrix.tau_offsets)
            tmp_k_points.append(symmetry_point)
        # == Remove duplicates and add the k point id and any further columns if any ==
        tmp_k_points.sort()
        if len(k_point) > 4:
            tmp_k_points = [[k_point[0]] + kp + k_point[4:].tolist() for i, kp in enumerate(tmp_k_points) \
                if i == 0 or kp != tmp_k_points[i-1]]
        else:
            tmp_k_points = [[k_point[0]] + kp for i, kp in enumerate(tmp_k_points) \
                if i == 0 or kp != tmp_k_points[i-1]]
        tot_k_points = tot_k_points + tmp_k_points
    return np.array(tot_k_points, dtype=float)
    
