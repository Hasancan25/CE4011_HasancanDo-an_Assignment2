"""
MODULE: Frame Analyzer Logic (Q2)
AUTHOR: Hasancan Dogan
PURPOSE: Implements structural analysis Steps 7-13. 
         Includes DOF labeling, global matrix assembly, and member end forces.
"""
import math
from matrix_lib import BandedSymmetricMatrix

class FrameAnalyzer:
    def __init__(self, num_node, num_elem, xy, m_props, connectivity, supports, loads):
        """
        PURPOSE: Initialize structural data.
        INPUTS: num_node, num_elem, xy coordinates, material props, connectivity, supports, loads.
        """
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
        """
        STEPS 7-8: Assignment of equation numbers.
        PURPOSE: Creates the E array to map DOFs to equation numbers.
        """
        self.E = [[0, 0, 0] for _ in range(self.num_node)]
        for sup in self.supports:
            n = sup[0] - 1
            self.E[n] = [sup[1], sup[2], sup[3]]
        
        count = 1
        for i in range(self.num_node):
            for j in range(3):
                if self.E[i][j] == 0: # Free
                    self.E[i][j] = count
                    count += 1
                else:
                    self.E[i][j] = 0 # Restrained
        self.num_eq = count - 1
        return self.num_eq

    def get_element_stiffness_global(self, elem_id):
        """
        STEP 9.1: Element stiffness in global coordinates.
        PURPOSE: Calculates $k = R^T \cdot k_{local} \cdot R$
        """
        conn = self.connectivity[elem_id]
        n1, n2 = conn[0]-1, conn[1]-1
        A, I, E_val = self.m_props[conn[2]-1]
        
        dx = self.xy[n2][0] - self.xy[n1][0]
        dy = self.xy[n2][1] - self.xy[n1][1]
        L = math.sqrt(dx**2 + dy**2)
        c, s = dx/L, dy/L

        # Local k'
        ax, b1, b2, b3, b4 = E_val*A/L, 12*E_val*I/L**3, 6*E_val*I/L**2, 4*E_val*I/L, 2*E_val*I/L
        kl = [[0.0]*6 for _ in range(6)]
        kl[0][0]=ax; kl[0][3]=-ax; kl[1][1]=b1; kl[1][2]=b2; kl[1][4]=-b1; kl[1][5]=b2
        kl[2][1]=b2; kl[2][2]=b3; kl[2][4]=-b2; kl[2][5]=b4; kl[3][0]=-ax; kl[3][3]=ax
        kl[4][1]=-b1; kl[4][2]=-b2; kl[4][4]=b1; kl[4][5]=-b2; kl[5][1]=b2; kl[5][2]=b4
        kl[5][4]=-b2; kl[5][5]=b3
        for i in range(6): 
            for j in range(i+1, 6): kl[j][i] = kl[i][j]

        # Rotation Matrix R
        R = [[0.0]*6 for _ in range(6)]
        R[0][0]=c; R[0][1]=s; R[1][0]=-s; R[1][1]=c; R[2][2]=1.0
        R[3][3]=c; R[3][4]=s; R[4][3]=-s; R[4][4]=c; R[5][5]=1.0

        # kg = R.T * kl * R
        kg = [[0.0]*6 for _ in range(6)]
        for i in range(6):
            for j in range(6):
                for m in range(6):
                    for n in range(6):
                        kg[i][j] += R[m][i] * kl[m][n] * R[n][j]
        return kg

    def assemble_global_matrix(self):
        """
        STEPS 9.2-10: Global matrix assembly.
        PURPOSE: Sums element matrices into the banded symmetric global K.
        """
        K = BandedSymmetricMatrix(self.num_eq, self.num_eq)
        for i in range(self.num_elem):
            ke = self.get_element_stiffness_global(i)
            n1, n2 = self.connectivity[i][0]-1, self.connectivity[i][1]-1
            g = [self.E[n1][0], self.E[n1][1], self.E[n1][2], 
                 self.E[n2][0], self.E[n2][1], self.E[n2][2]]
            for p in range(6):
                for q in range(6):
                    if g[p] > 0 and g[q] > 0:
                        K.assemble(g[p]-1, g[q]-1, ke[p][q])
        return K

    def construct_load_vector(self):
        """
        STEP 11: Global load vector F.
        """
        F = [0.0] * self.num_eq
        for load in self.loads:
            n = load[0]-1
            for q in range(3):
                Q = self.E[n][q]
                if Q > 0: F[Q-1] += load[q+1]
        return F

    def calculate_member_forces(self, global_displacements):
        """
        STEP 13: Member end forces in local coordinates.
        PURPOSE: Calculates internal forces (Axial, Shear, Moment).
        """
        print("\n" + "="*60)
        print(f"{'MEMBER END FORCES (LOCAL COORDINATES)':^60}")
        print("="*60)
        
        for i in range(self.num_elem):
            n1, n2 = self.connectivity[i][0]-1, self.connectivity[i][1]-1
            g = [self.E[n1][0], self.E[n1][1], self.E[n1][2], 
                 self.E[n2][0], self.E[n2][1], self.E[n2][2]]
            
            d_glob = [0.0] * 6
            for j in range(6):
                if g[j] > 0: d_glob[j] = global_displacements[g[j]-1]
            
            A, I, E_val = self.m_props[self.connectivity[i][2]-1]
            L = math.sqrt((self.xy[n2][0]-self.xy[n1][0])**2 + (self.xy[n2][1]-self.xy[n1][1])**2)
            c, s = (self.xy[n2][0]-self.xy[n1][0])/L, (self.xy[n2][1]-self.xy[n1][1])/L

            # Rotation Matrix R
            R = [[0.0]*6 for _ in range(6)]
            R[0][0]=c; R[0][1]=s; R[1][0]=-s; R[1][1]=c; R[2][2]=1.0
            R[3][3]=c; R[3][4]=s; R[4][3]=-s; R[4][4]=c; R[5][5]=1.0

            # Local d' = R * d_global
            d_loc = [0.0] * 6
            for row in range(6):
                for col in range(6): d_loc[row] += R[row][col] * d_glob[col]

            # Local stiffness k'
            ax, b1, b2, b3, b4 = E_val*A/L, 12*E_val*I/L**3, 6*E_val*I/L**2, 4*E_val*I/L, 2*E_val*I/L
            kl = [[0.0]*6 for _ in range(6)]
            kl[0][0]=ax; kl[0][3]=-ax; kl[1][1]=b1; kl[1][2]=b2; kl[1][4]=-b1; kl[1][5]=b2
            kl[2][1]=b2; kl[2][2]=b3; kl[2][4]=-b2; kl[2][5]=b4
            kl[3][0]=-ax; kl[3][3]=ax; kl[4][1]=-b1; kl[4][2]=-b2; kl[4][4]=b1; kl[4][5]=-b2
            kl[5][1]=b2; kl[5][2]=b4; kl[5][4]=-b2; kl[5][5]=b3
            for r in range(6):
                for col in range(r+1, 6): kl[col][r] = kl[r][col]

            # f_local = k_local * d_local
            f_loc = [0.0] * 6
            for row in range(6):
                for col in range(6): f_loc[row] += kl[row][col] * d_loc[col]

            print(f"ELEMENT {i+1}:")
            print(f"  Start Node -> Axial: {f_loc[0]:>10.2f} kN, Shear: {f_loc[1]:>10.2f} kN, Moment: {f_loc[2]:>10.2f} kNm")
            print(f"  End Node   -> Axial: {f_loc[3]:>10.2f} kN, Shear: {f_loc[4]:>10.2f} kN, Moment: {f_loc[5]:>10.2f} kNm")
            print("-" * 60)
