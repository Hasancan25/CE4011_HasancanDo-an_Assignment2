"""
PROGRAM: Frame Analysis - Final Verification
PURPOSE: Defines structural data and prints results matching PDF Page 6.
"""
from frame_analyzer import FrameAnalyzer
from solver import solve_banded_system

def main():
    # --- INPUT PHASE (Steps 1-6) [cite: 153-174] ---
    num_node = 4
    num_elem = 4
    # Coordinates [cite: 47-52]
    xy_coords = [[0.0, 0.0], [0.0, 3.0], [4.0, 3.0], [4.0, 0.0]]
    # Material properties (E = 2.0e5 MPa) [cite: 39-43]
    material_props = [[0.02, 0.08, 200000.0], [0.01, 0.01, 200000.0]]
    # Connectivity [cite: 54-57]
    connectivity = [[1, 2, 1], [2, 3, 1], [4, 3, 1], [1, 3, 2]]
    # Supports (Node 1: X,Y fixed; Node 4: Y fixed) [cite: 58-62]
    supports = [[1, 1, 1, 0], [4, 0, 1, 0]]
    # Loads (Nodes 2 & 3: 10kN X, -10kN Y) [cite: 63-67]
    loads = [[2, 10.0, -10.0, 0.0], [3, 10.0, -10.0, 0.0]]

    analyzer = FrameAnalyzer(num_node, num_elem, xy_coords, material_props, connectivity, supports, loads)
    num_eq = analyzer.label_active_dof()
    
    print(f"Sistem Serbestlik Derecesi (NumEq): {num_eq}") # Beklenen: 9 [cite: 198]
    
    K_global = analyzer.assemble_global_matrix()
    F_vector = analyzer.construct_load_vector()
    
    # Solve the system [cite: 264-265]
    displacements = solve_banded_system(K_global, F_vector)

    print("\n--- NODAL POINT DISPLACEMENTS (D) ---")
    for i, d in enumerate(displacements):
        print(f"D[{i+1}] = {d:.6e}")

if __name__ == "__main__":
    main()
