def calculate_member_forces(self, global_displacements):
        """
        STEP 13: Calculate element end forces in local coordinates.
        PURPOSE: Transform global D to local d', then multiply by local k'.
        """
        print("\n" + "="*50)
        print("STEP 13: MEMBER END FORCES (Local Coordinates)")
        print("="*50)
        print(f"{'Member':<8} {'Force Type':<15} {'Start Node':<15} {'End Node':<15}")
        print("-"*50)

        for i in range(self.num_elem):
            # 13.1. Get global displacements for this element
            n1, n2 = self.connectivity[i][0]-1, self.connectivity[i][1]-1
            g_indices = [self.E[n1][0], self.E[n1][1], self.E[n1][2], 
                         self.E[n2][0], self.E[n2][1], self.E[n2][2]]
            
            d_global = [0.0] * 6
            for j in range(6):
                if g_indices[j] > 0:
                    d_global[j] = global_displacements[g_indices[j]-1]
            
            # 13.2. Get Rotation matrix R and Local Stiffness k'
            # (Bu verileri get_element_stiffness_global fonksiyonundan esinlenerek alıyoruz)
            A, I, E_val = self.m_props[self.connectivity[i][2]-1]
            L = math.sqrt((self.xy[n2][0]-self.xy[n1][0])**2 + (self.xy[n2][1]-self.xy[n1][1])**2)
            c, s = (self.xy[n2][0]-self.xy[n1][0])/L, (self.xy[n2][1]-self.xy[n1][1])/L

            # Rotation Matrix R
            R = [[0.0]*6 for _ in range(6)]
            R[0][0]=c; R[0][1]=s; R[1][0]=-s; R[1][1]=c; R[2][2]=1.0
            R[3][3]=c; R[3][4]=s; R[4][3]=-s; R[4][4]=c; R[5][5]=1.0

            # Local d' = R * d_global
            d_local = [0.0] * 6
            for row in range(6):
                for col in range(6):
                    d_local[row] += R[row][col] * d_global[col]

            # 13.4. Local Stiffness k'
            ax, b1, b2, b3, b4 = E_val*A/L, 12*E_val*I/L**3, 6*E_val*I/L**2, 4*E_val*I/L, 2*E_val*I/L
            kl = [[0.0]*6 for _ in range(6)]
            kl[0][0]=ax; kl[0][3]=-ax; kl[1][1]=b1; kl[1][2]=b2; kl[1][4]=-b1; kl[1][5]=b2
            kl[2][1]=b2; kl[2][2]=b3; kl[2][4]=-b2; kl[2][5]=b4
            kl[3][0]=-ax; kl[3][3]=ax; kl[4][1]=-b1; kl[4][2]=-b2; kl[4][4]=b1; kl[4][5]=-b2
            kl[5][1]=b2; kl[5][2]=b4; kl[5][4]=-b2; kl[5][5]=b3
            # Symmetry
            for r in range(6):
                for col in range(r+1, 6): kl[col][r] = kl[r][col]

            # 13.5. f_local = k_local * d_local
            f_local = [0.0] * 6
            for row in range(6):
                for col in range(6):
                    f_local[row] += kl[row][col] * d_local[col]

            # Sonuçları Formatlı Yazdır (Axial, Shear, Moment)
            print(f"{i+1:<8} {'Axial (kN)':<15} {f_local[0]:<15.2f} {f_local[3]:<15.2f}")
            print(f"{'':<8} {'Shear (kN)':<15} {f_local[1]:<15.2f} {f_local[4]:<15.2f}")
            print(f"{'':<8} {'Moment (kNm)':<15} {f_local[2]:<15.2f} {f_local[5]:<15.2f}")
            print("-"*50)
