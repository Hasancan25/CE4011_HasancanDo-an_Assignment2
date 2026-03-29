"""
SOLVER: Linear System Solver for Banded Matrices
PURPOSE: Solves Ax = b using Gaussian Elimination optimized for banded storage.
"""

def solve_banded_system(banded_matrix, b_vector):
    n = banded_matrix.n
    bw = banded_matrix.bw
    b = b_vector[:]
    
    # Forward Elimination
    for k in range(n):
        pivot = banded_matrix.get_element(k, k)
        for i in range(k + 1, min(k + bw + 1, n)):
            factor = banded_matrix.get_element(k, i) / pivot
            for j in range(i, min(k + bw + 1, n)):
                # This part modifies the matrix in-place for efficiency
                # Since we don't use NumPy, we access the data directly
                old_val = banded_matrix.get_element(i, j)
                sub_val = factor * banded_matrix.get_element(k, j)
                # Note: This is simplified for the assignment
                banded_matrix.data[i][j-i] = old_val - sub_val
            b[i] -= factor * b[k]

    # Back Substitution
    x = [0.0] * n
    for i in range(n - 1, -1, -1):
        sum_ax = 0.0
        for j in range(i + 1, min(i + bw + 1, n)):
            sum_ax += banded_matrix.get_element(i, j) * x[j]
        x[i] = (b[i] - sum_ax) / banded_matrix.get_element(i, i)
    return x
