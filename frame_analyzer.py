"""
MODULE: Frame Analyzer Logic
AUTHOR: Hasancan Dogan
PURPOSE: Handles Step 7-13 of the frame analysis algorithm.
"""
import math
from matrix_lib import BandedSymmetricMatrix

class FrameAnalyzer:
    def __init__(self, num_node, num_elem, xy, m_props, connectivity, supports, loads):
        """Initialize analyzer with structural data [cite: 153-174]."""
        self.num_node = num_node
        self.num_elem = num_elem
        self.xy = xy
        self.m_props = m_props
        self.connectivity = connectivity
        self.supports = supports
        self.loads = loads
        self.E = []
        self.num_eq = 0

    def label_active_dof(self):
        """Step 7-8: Equation numbering for active degrees of freedom [cite: 176-198]."""
        self.E = [[-1, -1, -1] for _ in range(self.num_node)]
        for sup in self.supports:
            node_id = sup[0] - 1
            self.E[node_id][0] = 0 if sup[1] == 1 else -1
            self.E[node_id][1] = 0 if sup[2] == 1 else -1
            self.E[node_id][2] = 0 if sup[3] == 1 else -1
        
        count = 1
        for i in range(self.num_node):
            for j in range(3):
                if self.E[i][j] == -1:
                    self.E[i][j] = count
                    count += 1
        self.num_eq = count - 1
        return self.num_eq

    def get_element_stiffness_global(self, elem_id):
        """Step 9.1: Form element stiffness matrix in global coordinates [cite: 507-575]."""
        conn = self.connectivity[elem_id]
        n1, n2 = conn[0]-1, conn[1]-1
        A, I, E = self.m_props[conn[2]-1]
        
        dx = self.xy[n2][0] - self.xy[n1][0]
        dy = self.xy[n2][1] - self.xy[n1][1]
        L = math.sqrt(dx**2 + dy**2)
        c, s = dx/L, dy/L
        
        a = E*A/L
        b1, b2, b3, b4 = 12*E*I/L**3, 6*E*I/L**2, 4*E*I/L, 2*E*I/L

        k = [[0.0]*6 for _ in range(6)]
        # Top-left and Bottom-right blocks
        k[0][0] = a*c*c + b1*s*s; k[0][1] = (a-b1)*c*s; k[0][2] = -b2*s
        k[1][1] = a*s*s + b1*c*c; k[1][2] = b2*c
        k[2][2] = b3
        
        k[3][3] = k[0][0];        k[3][4] = k[0][1];   k[3][5] = b2*s
        k[4][4] = k[1][1];        k[4][5] = -b2*c
        k[5][5] = b3

        # Cross blocks
        k[0][3] = -k[0][0];       k[0][4] = -k[0][1];  k[0][5] = -b2*s
        k[1][3] = -k[0][1];       k[1][4] = -k[1][1];  k[1][5] = -b2*c
        k[2][3] = b2*s;           k[2][4] = -b2*c;     k[2][5] = b4
        
        # Mirror for symmetry
        for i in range(6):
            for j in range(i + 1, 6):
                k[j][i] = k[i][j]
        return k

    def assemble_global_matrix(self):
        """Step 9.2-10: Assemble element stiffness into global K [cite: 205-254]."""
        # Global K using Banded Symmetric storage for NumEq
        K_global = BandedSymmetricMatrix(self.num_eq, 9) 
        for i in range(self.num_elem):
            k_elem = self.get_element_stiffness_global(i)
            n1, n2 = self.connectivity[i][0]-1, self.connectivity[i][1]-1
            g = [self.E[n1][0], self.E[n1][1], self.E[n1][2], 
                 self.E[n2][0], self.E[n2][1], self.E[n2][2]]
            for p in range(6):
                for q in range(6):
                    if g[p] > 0 and g[q] > 0:
                        K_global.assemble(g[p]-1, g[q]-1, k_elem[p][q])
        return K_global

    def construct_load_vector(self):
        """Step 11: Construct global load vector F [cite: 255-263]."""
        F = [0.0] * self.num_eq
        for load in self.loads:
            n = load[0]-1
            for q in range(3):
                Q = self.E[n][q]
                if Q > 0: F[Q-1] += load[q+1]
        return F
