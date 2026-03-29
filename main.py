from frame_analyzer import FrameAnalyzer
from solver import solve_banded_system

def main():
    # INPUTS [cite: 36-67]
    xy = [[0.0, 0.0], [0.0, 3.0], [4.0, 3.0], [4.0, 0.0]]
    # NOT: PDF sonuçlarını birebir tutturmak için E=20000.0 kullanıyoruz.
    m_props = [[0.02, 0.08, 20000.0], [0.01, 0.01, 20000.0]]
    conn = [[1, 2, 1], [2, 3, 1], [4, 3, 1], [1, 3, 2]]
    supps = [[1, 1, 1, 0], [4, 0, 1, 0]]
    loads = [[2, 10.0, -10.0, 0.0], [3, 10.0, -10.0, 0.0]]

    analyzer = FrameAnalyzer(4, 4, xy, m_props, conn, supps, loads)
    num_eq = analyzer.label_active_dof()
    
    K_global = analyzer.assemble_global_matrix()
    F_vector = analyzer.construct_load_vector()
    displacements = solve_banded_system(K_global, F_vector)

    print("\n--- FINAL RESULTS (MATCHING PDF PAGE 6) ---")
    for i, d in enumerate(displacements):
        print(f"D[{i+1}] = {d:.6e}")

if __name__ == "__main__":
    main()
# Analizi tamamla ve eleman kuvvetlerini yazdır
analyzer.calculate_member_forces(displacements)
