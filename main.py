from frame_analyzer import FrameAnalyzer
from solver import solve_banded_system

def main():
    # --- INPUT PHASE (Step 1-6) [cite: 45-67] ---
    num_node = 4
    num_elem = 4
    xy_coords = [[0.0, 0.0], [0.0, 3.0], [4.0, 3.0], [4.0, 0.0]]
    # E değerini 200,000.0 (2.0e5) olarak düzelttik 
    material_props = [[0.02, 0.08, 200000.0], [0.01, 0.01, 200000.0]]
    connectivity = [[1, 2, 1], [2, 3, 1], [4, 3, 1], [1, 3, 2]]
    supports = [[1, 1, 1, 0], [4, 0, 1, 0]]
    loads = [[2, 10.0, -10.0, 0.0], [3, 10.0, -10.0, 0.0]]

    analyzer = FrameAnalyzer(num_node, num_elem, xy_coords, material_props, connectivity, supports, loads)
    num_eq = analyzer.label_active_dof()
    
    print(f"Toplam Serbestlik Derecesi (NumEq): {num_eq}") # Beklenen: 9
    
    K_global = analyzer.assemble_global_matrix()
    F_vector = analyzer.construct_load_vector()
    
    # Denklem çözümü [cite: 157-158]
    displacements = solve_banded_system(K_global, F_vector)

    print("\n--- DOĞRULAMA: YER DEĞİŞTİRMELER (D) ---")
    for i, d in enumerate(displacements):
        print(f"D[{i+1}] = {d:.6e}")

if __name__ == "__main__":
    main()
