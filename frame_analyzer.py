def get_element_stiffness_global(self, elem_id):
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
        # Diagonal blocks (Symmetric parts)
        k[0][0] = a*c*c + b1*s*s; k[0][1] = (a-b1)*c*s; k[0][2] = -b2*s
        k[1][1] = a*s*s + b1*c*c; k[1][2] = b2*c
        k[2][2] = b3
        
        k[3][3] = k[0][0];        k[3][4] = k[0][1];   k[3][5] = b2*s
        k[4][4] = k[1][1];        k[4][5] = -b2*c
        k[5][5] = b3

        # Off-diagonal blocks
        k[0][3] = -k[0][0];       k[0][4] = -k[0][1];  k[0][5] = -b2*s
        k[1][3] = -k[0][1];       k[1][4] = -k[1][1];  k[1][5] = -b2*c
        k[2][3] = b2*s;           k[2][4] = -b2*c;     k[2][5] = b4
        
        # Fill symmetry
        for i in range(6):
            for j in range(i + 1, 6):
                k[j][i] = k[i][j]
        return k
