"""
PROGRAM: Frame Analysis Software - Final Execution
AUTHOR: Hasancan Dogan
PURPOSE: Performs full analysis of the sample frame and prints verification results.
"""

from frame_analyzer import FrameAnalyzer
from solver import solve_banded_system

def main():
    # --- A. INPUT PHASE (Steps 1-6) [cite: 152-174] ---
    num_node = 4
    num_elem = 4
    # Coordinates [cite: 158]
    xy_coords = [[0.0, 0.0], [0.0, 3.0], [4.0, 3.0], [4.0, 0.0]]
    # Material properties [cite: 161]
    material_props = [[0.02, 0.08, 20000.0], [0.01, 0.01, 20000.0]]
    # Connectivity [cite: 165]
    connectivity = [[1, 2, 1], [2, 3, 1], [4, 3, 1], [1, 3, 2]]
    # Support conditions [cite: 170]
    supports = [[1, 1, 1, 0], [4, 0, 1, 0]]
    # Applied loads [cite: 174]
    loads = [[2, 10.0, -10.0, 0.0], [3, 10.0, -10.0, 0.0]]

    # Initialize the Analyzer
    analyzer = FrameAnalyzer(num_node, num_elem, xy_coords, material_props, connectivity, supports, loads)

    # --- B. EQUATION NUMBERING (Steps 7-8) [cite: 175-198] ---
    num_eq = analyzer.label_active_dof()
    print(f"Sistem Serbestlik Derecesi (NumEq): {num_eq}") # Beklenen: 9

    # --- C. GLOBAL MATRIX ASSEMBLY (Steps 9-10) [cite: 199-254] ---
    print("Küresel Rijitlik Matrisi (K) Montajlanıyor...")
    K_global = analyzer.assemble_global_matrix()

    # --- D. GLOBAL LOAD VECTOR (Step 11) [cite: 255-263] ---
    F_vector = analyzer.construct_load_vector()
    print("Yük Vektörü (F):", [round(f, 2) for f in F_vector])

    # --- E. STRUCTURAL DISPLACEMENTS (Step 12) [cite: 264-265] ---
    print("Sistem Çözülüyor (Gauss Elimination)...")
    # Kendi solver'ımızı ve banded matrix yapımızı kullanıyoruz
    displacements = solve_banded_system(K_global, F_vector)

    # --- G. RESULTS VERIFICATION [cite: 333-342] ---
    print("\n--- ANALİZ SONUÇLARI: DÜĞÜM NOKTASI YER DEĞİŞTİRMELERİ (D) ---")
    # PDF'deki sonuçlarla karşılaştırma için
    for i, d in enumerate(displacements):
        print(f"D[{i+1}] = {d:.6e}")

if __name__ == "__main__":
    main()
