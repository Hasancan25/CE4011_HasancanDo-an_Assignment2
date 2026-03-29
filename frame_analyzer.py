"""
MODULE: Frame Analyzer Logic
AUTHOR: Hasancan Dogan
PURPOSE: Correct implementation of global stiffness assembly to match PDF results.
"""
import math
from matrix_lib import BandedSymmetricMatrix

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
        """Step 7-8: Equation numbering [cite: 68-91]."""
        self.E = [[0, 0, 0] for _ in range(self.num_node)]
        # Boundary conditions
        for sup in self.supports:
            node_id = sup[0] - 1
            self.E[node_id][0] = sup[1] # X
            self.E[node_id][1] = sup[2] # Y
            self.E[node_id][2] = sup[3] # Rot

        # Labeling active DOFs starting from 1
        count = 1
        for i in range(self.num_node):
            for j in range(3):
                if self.E[i][j] == 0: # If free
                    self.E[i][j] = count
                    count += 1
                else:
                    self.E[i][j] = 0 # If restrained
        self.num_eq = count - 1
        return self.num_eq

    def get_element_stiffness_global(self, elem_id):
        """Step 9.1: Global k using R.T * k_local * R ."""
        conn = self.connectivity[elem_id]
        n1, n2 = conn[0]-1, conn[1]-1
        A, I, E = self.m_props[conn[2]-1]
        
        dx = self.xy[n2][0] - self.xy[n1][0]
        dy = self.xy[n2][1] - self.xy[n1][1]
        L = math.sqrt(dx**2 + dy**2)
        c, s = dx/L, dy/L

        # Local Stiffness k' [cite: 266-288]
        axial = E * A / L
        b1 = 12 * E * I / (L**3)
        b2 = 6 * E * I / (L**2)
        b3 = 4 * E * I / L
        b4 = 2 * E * I / L

        k_local = [[0.0]*6 for _ in range(6)]
        k_local[0][0]=axial; k_local[0][3]=-axial
        k_local[1][1]=b1;    k_local[1][2]=b2;   k_local[1][4]=-b1;  k_local[1][5]=b2
        k_local[2][1]=b2;    k_local[2][2]=b3;   k_local[2][4]=-b2;  k_local[2][5]=b4
        k_local[3][0]=-axial;k_local[3][3]=axial
        k_local[4][1]=-b1;   k_local[4][2]=-b2;  k_local[4][4]=b1;   k_local[4][5]=-b2
        k_local[5][1]=b2;    k_local[5][2]=b4;   k_local[5][4]=-b2;  k_local[5][5]=b3

        # Rotation Matrix R [cite: 360-362]
        R = [[0.0]*6 for _ in range(6)]
        R[0][0]=c; R[0][1]=s; R[1][0]=-s; R[1][1]=c; R[2][2]=1.0
        R[3][3]=c; R[3][4]=s; R[4][3]=-s; R[4][4]=c; R[5][5]=1.0

        # Global Matrix: k = R.T * k_local * R
        k_global = [[0.0]*6 for _ in range(6)]
        for i in range(6):
            for j in range(6):
                for m in range(6):
                    for n in range(6):
                        k_global[i][j] += R[m][i] * k_local[m][n] * R[n][j]
        return k_global

    def assemble_global_matrix(self):
        """Step 9.2-10: Global assembly [cite: 92-147]."""
        K_global = BandedSymmetricMatrix(self.num_eq, self.num_eq) # Safe bandwidth
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
        """Step 11: Load vector F [cite: 148-156]."""
        F = [0.0] * self.num_eq
        for load in self.loads:
            n = load[0]-1
            for q in range(3):
                Q = self.E[n][q]
                if Q > 0: F[Q-1] += load[q+1]
        return F
