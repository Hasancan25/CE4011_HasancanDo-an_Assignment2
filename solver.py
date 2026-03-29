"""
SOLVER MODULE: Linear System Solver for Banded Matrices
AUTHOR: Hasancan Dogan
PURPOSE: Solves Ax = b using Gaussian Elimination optimized for banded symmetric storage.
NOTE: Matches logic from CE 4011 Lecture Notes (Step 12).
"""

def solve_banded_system(banded_matrix, b_vector):
    """
    PURPOSE: Solve the global system [K]{D} = {F} for displacements {D}.
    INPUTS: 
        banded_matrix: A BandedSymmetricMatrix object (Global K).
        b_vector: A list of floats (Global Load Vector F).
    OUTPUTS: 
        x: A list of floats (Displacement Vector D).
    ASSUMPTIONS: Matrix is non-singular and symmetric.
    """
    n = banded_matrix.n
    bw = banded_matrix.bw
    
    # Working copies to avoid modifying original inputs
    # As suggested in Gauss notes [cite: 743, 744]
    b = b_vector[:]
    
    # 1. FORWARD ELIMINATION (Optimized for bandwidth)
    # Only loops through elements within the band [cite: 805, 806]
    for k in range(n):
        pivot = banded_matrix.get_element(k, k)
        
        if abs(pivot) < 1e-15:
            raise ValueError("Matrix is singular or nearly singular.")
        
        # Only eliminate rows within the bandwidth
        for i in range(k + 1, min(k + bw + 1, n)):
            factor = banded_matrix.get_element(k, i) / pivot
            
            # Update subsequent entries in the current row within band
            for j in range(i, min(k + bw + 1, n)):
                new_val = banded_matrix.get_element(i, j) - (factor * banded_matrix.get_element(k, j))
                banded_matrix.assemble(i, j, -banded_matrix.get_element(i, j) + new_val) # Reset and Update
                
            # Update the load vector [cite: 810]
            b[i] = b[i] - (factor * b[k])

    # 2. BACK SUBSTITUTION
    # Starting from the last equation [cite: 813, 1014]
    x = [0.0] * n
    for i in range(n - 1, -1, -1):
        sum_ax = 0.0
        # Only sum elements within the bandwidth
        for j in range(i + 1, min(i + bw + 1, n)):
            sum_ax += banded_matrix.get_element(i, j) * x[j]
        
        x[i] = (b[i] - sum_ax) / banded_matrix.get_element(i, i)
        
    return x
