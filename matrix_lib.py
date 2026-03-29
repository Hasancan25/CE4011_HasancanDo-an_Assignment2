"""
LIBRARY: Matrix Operations Library
AUTHOR: Hasancan Dogan
PURPOSE: Custom matrix library for structural analysis. Handles symmetric banded storage.
NOTE: No external libraries like NumPy are used as per assignment requirements.
"""

class Matrix:
    """
    PURPOSE: Handles general dense matrix operations (used for element-level k).
    ATTRIBUTES: rows (int), cols (int), data (list of lists)
    """
    def __init__(self, rows, cols):
        self.rows = rows
        self.cols = cols
        self.data = [[0.0 for _ in range(cols)] for _ in range(rows)]

    def set_value(self, i, j, value):
        """Sets value at row i, col j."""
        self.data[i][j] = value

    def get_value(self, i, j):
        """Returns value at row i, col j."""
        return self.data[i][j]

class BandedSymmetricMatrix:
    """
    PURPOSE: Efficiently stores global stiffness matrix K using symmetry and bandwidth.
    REQUIREMENT: Minimizes storage by only keeping the diagonal and upper band.
    ATTRIBUTES: n (total equations), bw (half bandwidth)
    """
    def __init__(self, size, half_bandwidth):
        """
        INPUTS: size (total D.O.F), half_bandwidth (max distance from diagonal)
        PURPOSE: Initializes a 2D list for compact storage.
        """
        self.n = size
        self.bw = half_bandwidth
        # Compact storage: rows = size, cols = bw + 1 (diagonal + upper band)
        # data[i][j] stores the element at global index K[i][i+j]
        self.data = [[0.0 for _ in range(self.bw + 1)] for _ in range(self.n)]

    def assemble(self, i, j, value):
        """
        PURPOSE: Assembles a term into the global banded matrix.
        INPUTS: i, j (global indices), value (stiffness term)
        ASSUMPTION: Matrix is symmetric; only upper triangle is stored.
        """
        # Ensure i is the smaller index for upper triangle storage 
        if i > j:
            i, j = j, i
            
        offset = j - i
        if 0 <= offset <= self.bw:
            self.data[i][offset] += value

    def get_element(self, i, j):
        """
        PURPOSE: Retrieves element from compact storage using symmetry.
        """
        if i > j:
            i, j = j, i
        offset = j - i
        if 0 <= offset <= self.bw:
            return self.data[i][offset]
        return 0.0
