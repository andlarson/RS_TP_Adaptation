"""
This file makes it easy to test the offerings of TetGen.
"""

import tetgen



if __name__ == "__main__":
    
    stl_file_path = "/home/andlars/Desktop/RS_TP_Adaptation/experiments/experiments/"\
                    "stl_to_volumetric_mesh/select_ascii_stl_from_ABC_dataset/00000029.stl"

    tgen = tetgen.TetGen(stl_file_path)

    nodes, elem = tgen.tetrahedralize(quality=True, verbose=2)
    tgen.grid.plot(show_edges=True)

    # tgen.plot(show_edges=True)
    
