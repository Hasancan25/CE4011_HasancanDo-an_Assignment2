def get_element_stiffness_global(self, elem_id):
        """
        Step 9.1: Form element stiffness matrix using k = R.T * k_local * R.
        PURPOSE: Eliminates manual term errors by performing explicit matrix multiplication.
        """
        conn = self.connectivity[elem_id]
        n1, n2 = conn[0]-1, conn[1]-1
        A, I, E = self.m_props[conn[2]-1]
        
        dx = self.xy[n2][0] - self.xy[n1][0]
        dy = self.xy[n2][1] - self.xy[n1][1]
        L = math.sqrt(dx**2 + dy**2)
        c, s = dx/L, dy/L # cos and sin
        
        # 1. Local Stiffness Matrix (k') [cite: 372-410]
        # order: [u1', v1', th1, u2', v2', th2]
        axial = E * A / L
        b1 = 12 * E * I / (L**3)
        b2 = 6 * E * I / (L**2)
        b3 = 4 * E * I / L
        b4 = 2 * E * I / L

        k_local = [[0.0]*6 for _ in range(6)]
        # Top half
        k_local[0][0] = axial; k_local[0][3] = -axial
        k_local[1][1] = b1;    k_local[1][2] = b2;    k_local[1][4] = -b1;   k_local[1][5] = b2
        k_local[2][2] = b3;    k_local[2][4] = -b2;   k_local[2][5] = b4
        k_local[3][3] = axial
        k_local[4][4] = b1;    k_local[4][5] = -b2
        k_local[5][5] = b3
        # Fill symmetry
        for i in range(6):
            for j in range(i+1, 6): k_local[j][i] = k_local[i][j]

        # 2. Rotation Matrix (R) [cite: 362]
        R = [[0.0]*6 for _ in range(6)]
        R[0][0]=c;  R[0][1]=s;  R[1][0]=-s; R[1][1]=c; R[2][2]=1.0
        R[3][3]=c;  R[3][4]=s;  R[4][3]=-s; R[4][4]=c; R[5][5]=1.0

        # 3. Global Stiffness Calculation: k = R.T * k_local * R
        # Manual multiplication for safety (No NumPy)
        # Temp = k_local * R
        temp = [[0.0]*6 for _ in range(6)]
        for i in range(6):
            for j in range(6):
                for m in range(6):
                    temp[i][j] += k_local[i][m] * R[m][j]
        
        # k_global = R.T * Temp
        k_global = [[0.0]*6 for _ in range(6)]
        for i in range(6):
            for j in range(6):
                for m in range(6):
                    # R.T[i][m] is R[m][i]
                    k_global[i][j] += R[m][i] * temp[m][j]
                    
        return k_global
