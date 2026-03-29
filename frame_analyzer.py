import math
from matrix_lib import BandedSymmetricMatrix

class FrameAnalyzer:
    def __init__(self, num_node, num_elem, xy, m_props, connectivity, supports, loads):
        self.num_node, self.num_elem, self.xy = num_node, num_elem, xy
        self.m_props, self.connectivity = m_props, connectivity
        self.supports, self.loads = supports, loads
        self.E = []

    def label_active_dof(self):
        """D.O.F etiketleme (Step 7-8) [cite: 68-91]."""
        self.E = [[0, 0, 0] for _ in range(self.num_node)]
        for sup in self.supports:
            n = sup[0]-1
            self.E[n] = [sup[1], sup[2], sup[3]]
        
        count = 1
        for i in range(self.num_node):
            for j in range(3):
                if self.E[i][j] == 0: # Free
                    self.E[i][j] = count
                    count += 1
                else: self.E[i][j] = 0 # Restrained
        self.num_eq = count - 1
        return self.num_eq

    def get_element_stiffness_global(self, elem_id):
        """Küresel k matrisi hesabı (Step 9.1) ."""
        conn = self.connectivity[elem_id]
        n1, n2 = conn[0]-1, conn[1]-1
        A, I, E = self.m_props[conn[2]-1]
        
        L = math.sqrt((self.xy[n2][0]-self.xy[n1][0])**2 + (self.xy[n2][1]-self.xy[n1][1])**2)
        c, s = (self.xy[n2][0]-self.xy[n1][0])/L, (self.xy[n2][1]-self.xy[n1][1])/L
        
        # Local matrix k' [cite: 266-296]
        ax, b1, b2, b3, b4 = E*A/L, 12*E*I/L**3, 6*E*I/L**2, 4*E*I/L, 2*E*I/L
        kl = [[0.0]*6 for _ in range(6)]
        kl[0][0]=ax; kl[0][3]=-ax; kl[3][0]=-ax; kl[3][3]=ax
        kl[1][1]=b1; kl[1][2]=b2; kl[1][4]=-b1; kl[1][5]=b2
        kl[2][1]=b2; kl[2][2]=b3; kl[2][4]=-b2; kl[2][5]=b4
        kl[4][1]=-b1; kl[4][2]=-b2; kl[4][4]=b1; kl[4][5]=-b2
        kl[5][1]=b2; kl[5][2]=b4; kl[5][4]=-b2; kl[5][5]=b3

        # Rotation matrix R [cite: 360-362]
        R = [[0.0]*6 for _ in range(6)]
        R[0][0]=c; R[0][1]=s; R[1][0]=-s; R[1][1]=c; R[2][2]=1.0
        R[3][3]=c; R[3][4]=s; R[4][3]=-s; R[4][4]=c; R[5][5]=1.0

        # k = R.T * k_local * R
        kg = [[0.0]*6 for _ in range(6)]
        for i in range(6):
            for j in range(6):
                for m in range(6):
                    for n in range(6):
                        kg[i][j] += R[m][i] * kl[m][n] * R[n][j]
        return kg

    def assemble_global_matrix(self):
        """Global assembly (Step 9.2-10) [cite: 92-147]."""
        K = BandedSymmetricMatrix(self.num_eq, self.num_eq)
        for i in range(self.num_elem):
            ke = self.get_element_stiffness_global(i)
            n1, n2 = self.connectivity[i][0]-1, self.connectivity[i][1]-1
            g = [self.E[n1][0], self.E[n1][1], self.E[n1][2], self.E[n2][0], self.E[n2][1], self.E[n2][2]]
            for p in range(6):
                for q in range(6):
                    if g[p] > 0 and g[q] > 0:
                        K.assemble(g[p]-1, g[q]-1, ke[p][q])
        return K

    def construct_load_vector(self):
        """Load vector F (Step 11) [cite: 148-156]."""
        F = [0.0] * self.num_eq
        for l in self.loads:
            n = l[0]-1
            for q in range(3):
                Q = self.E[n][q]
                if Q > 0: F[Q-1] += l[q+1]
        return F
