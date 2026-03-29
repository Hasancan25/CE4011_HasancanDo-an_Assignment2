def solve_banded_system(banded_matrix, b_vector):
    n = banded_matrix.n
    bw = banded_matrix.bw
    b = b_vector[:]
    
    # 1. Forward Elimination
    for k in range(n):
        pivot = banded_matrix.get_element(k, k)
        if abs(pivot) < 1e-18:
            continue # Skip near-zero diagonals if any (shouldn't happen with correct K)
        
        for i in range(k + 1, min(k + bw + 1, n)):
            factor = banded_matrix.get_element(k, i) / pivot
            for j in range(i, min(k + bw + 1, n)):
                val = banded_matrix.get_element(i, j) - factor * banded_matrix.get_element(k, j)
                banded_matrix.set_element(i, j, val)
            b[i] -= factor * b[k]

    # 2. Back Substitution
    x = [0.0] * n
    for i in range(n - 1, -1, -1):
        sum_ax = sum(banded_matrix.get_element(i, j) * x[j] for j in range(i + 1, min(i + bw + 1, n)))
        x[i] = (b[i] - sum_ax) / banded_matrix.get_element(i, i)
    return x
