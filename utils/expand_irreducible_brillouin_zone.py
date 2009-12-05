import numpy as np

def expand_irreducible_brillouin_zone(symmetry_matrices, k_points):
    tot_k_points = []
    for k_point in k_points:
        tmp_k_points = []
        for symmetry_matrix in symmetry_matrices:
            symmetry_point = list(np.dot(k_point[1:4], symmetry_matrix))
            tmp_k_points.append(symmetry_point)
        # Remove duplicates and add the k point id
        tmp_k_points.sort()
        tmp_k_points = [[k_point[0]] + kp for i, kp in enumerate(tmp_k_points) \
            if i == 0 or kp != tmp_k_points[i-1]]
        tot_k_points = tot_k_points + tmp_k_points
    return np.array(tot_k_points, dtype=float)
    