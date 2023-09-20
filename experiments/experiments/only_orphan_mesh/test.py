import core.machining.machining as mach
import util.geom as geom
import core.part.part as part
import core.tool_pass.tool_pass as tp

# DEBUG
import core.simulation.simulation as sim
import core.abaqus.abaqus_shim as shim

import sys
import os


"""
This script maps an orphan mesh to a geometry in a simple case. 
"""


path_to_cae_with_only_orphan_mesh = "/home/andlars/Desktop/RS_TP_Adaptation/software/script_testing/only_orphan_mesh/only_orphan_mesh.cae" 
mdb = shim.use_mdb(path_to_cae_with_only_orphan_mesh)
sim.orphan_mesh_to_geometry("PART-1", "Model-1", mdb)
save_loc = "/home/andlars/Desktop/RS_TP_Adaptation/software/script_testing/only_orphan_mesh/test.cae"
shim.save_mdb_as(save_loc, mdb)
shim.close_mdb(mdb)

