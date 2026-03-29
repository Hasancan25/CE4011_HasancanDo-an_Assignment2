"""
PROGRAM: Frame Analysis Software - Final Execution
AUTHOR: Hasancan Dogan
DATE: 29.03.2026
PURPOSE: Performs full linear static analysis of a 2D frame structure (Steps 1-13).
         Matches verification results from CE 4011 assignment documentation.
"""

from frame_analyzer import FrameAnalyzer
from solver import solve_banded_system

def main():
    # --- STEP A: INPUT PHASE (Steps 1-6) ---
    num_node = 4
    num_elem = 4
    
    # Nodal Coordinates (X, Y) [m]
    xy_coords = [
        [0.0, 0.0], # Node 1
        [0.0, 3.0], # Node 2
        [4.0, 3.0], # Node 3
        [4.0, 0.0]  # Node 4
    ]

    # Material Properties [Area (m2), Inertia (m4), Elasticity (MPa)]
    # Note: Using E = 20000.0 to match the numerical results on PDF Page 6.
    material_props = [
        [0.02, 0.08, 20000.0], # Set 1
        [0.01, 0.01, 20000.0]  # Set 2
    ]

    # Element Connectivity [StartNode, EndNode, MaterialSetID]
    connectivity = [
        [1, 2, 1], # Element 1
        [2, 3, 1], # Element 2
        [4, 3, 1], # Element 3
        [1, 3, 2]  # Element 4
    ]

    # Boundary Conditions [NodeID, X_code, Y_code, Z_rot_code] (1:Restrained, 0:Free)
    supports = [
        [1, 1, 1, 0], # Node 1: Fixed X, Y; Free Rotation
        [4, 0, 1, 0]  # Node 4: Fixed Y; Free X and Rotation
    ]

    # Applied Nodal Loads [NodeID, Force_X (kN), Force_Y (kN), Moment_Z (kNm)]
    loads = [
        [2, 10.0, -10.0, 0.0],
        [3, 10.0, -10.0, 0.0]
    ]

    # --- EXECUTION PHASE ---
    # Initialize the structural analyzer
    analyzer = FrameAnalyzer(num_node, num_elem, xy_coords, material_props, connectivity, supports, loads)

    # STEP 7-8: Assignment of equation numbers
    num_eq = analyzer.label_active_dof()
    print(f"\nAnalysis initialized. Total Active Degrees of Freedom: {num_eq}")

    # STEP 9-10: Global matrix assembly
    print("Assembling global banded symmetric stiffness matrix...")
    K_global = analyzer.assemble_global_matrix()

    # STEP 11: Global load vector construction
    F_vector = analyzer.construct_load_vector()

    # STEP 12: Solve the system [K]{D} = {F}
    print("Solving for nodal displacements using Gaussian Elimination...")
    displacements = solve_banded_system(K_global, F_vector)

    # --- OUTPUT PHASE ---
    print("\n" + "="*50)
    print("STEP 12: NODAL POINT DISPLACEMENTS (D)")
    print("="*50)
    print(f"{'DOF ID':<10} {'Displacement Value':<20}")
    print("-" * 35)
    for i, d in enumerate(displacements):
        print(f"D[{i+1:<2}]      {d:+.6e}")

    # STEP 13: Calculate and print member end forces
    analyzer.calculate_member_forces(displacements)

if __name__ == "__main__":
    main()
