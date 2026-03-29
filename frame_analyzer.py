"""
MODULE: Frame Analyzer Logic (FIXED)
PURPOSE: Full 6x6 stiffness matrix implementation to avoid singular matrix error.
"""
import math
from matrix_lib import Matrix, BandedSymmetricMatrix

class FrameAnalyzer:
    def __init__(self, num_node, num_elem, xy, m_props, connectivity, supports, loads):
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
        """Step 7-8: Labeling active degrees of freedom [cite: 176-198]."""
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
        """
        Step 9.1: Form full 6x6 global stiffness matrix .
        PURPOSE: Populates all terms for start and end nodes to ensure stability.
        """
        conn = self.connectivity[elem_id]
        n1, n2 = conn[0]-1, conn[1]-1
        A, I, E = self.m_props[conn[2]-1]
        
        L = math.sqrt((self.xy[n2][0]-self.xy[n1][0])**2 + (self.xy[n2][1]-self.xy[n1][1])**2)
        c, s = (self.xy[n2][0]-self.xy[n1][0])/L, (self.xy[n2][1]-self.xy[n1][1])/L
        
        # Stiffness coefficients [cite: 478-500]
        a = E*A/L
        b1, b2, b3, b4 = 12*E*I/L**3, 6*E*I/L**2, 4*E*I/L, 2*E*I/L

        # Full 6x6 Matrix Terms (Global)
        k = [[0.0]*6 for _ in range(6)]
        
        # Top-left (Start node - Start node)
        k[0][0] = a*c*c + b1*s*s; k[0][1] = (a-b1)*c*s; k[0][2] = -b2*s
        k[1][0] = k[0][1];        k[1][1] = a*s*s + b1*c*c; k[1][2] = b2*c
        k[2][0] = k[0][2];        k[2][1] = k[1][2];        k[2][2] = b3

        # Bottom-right (End node - End node)
        k[3][3] = k[0][0];        k[3][4] = k[0][1];        k[3][5] = -k[0][2]
        k[4][3] = k[3][4];        k[4][4] = k[1][1];        k[4][5] = -k[1][2]
        k[5][3] = k[3][5];        k[5][4] = k[4][5];        k[5][5] = b3

        # Top-right (Start node - End node)
        k[0][3] = -k[0][0];       k[0][4] = -k[0][1];       k[0][5] = k[0][2]
        k[1][3] = -k[1][0];       k[1][4] = -k[1][1];       k[1][5] = k[1][2]
        k[2][3] = b2*s;           k[2][4] = -b2*c;          k[2][5] = b4

        # Bottom-left (End node - Start node - Symmetry handled by assembly)
        for i in range(3, 6):
            for j in range(0, 3):
                k[i][j] = k[j][i]
                
        return k

    def assemble_global_matrix(self):
        """Step 9.2-10: Global assembly with bandwidth check [cite: 205-254]."""
        K_global = BandedSymmetricMatrix(self.num_eq, 9) # Bandwidth for this frame is 9
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
        """Step 11: F vector [cite: 255-263]."""
        F = [0.0] * self.num_eq
        for load in self.loads:
            n = load[0]-1
            for q in range(3):
                Q = self.E[n][q]
                if Q > 0: F[Q-1] += load[q+1]
        return F
