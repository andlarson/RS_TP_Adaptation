import sys

# DEBUG
# For Abaqus PDB and allows main.py to be run from any directory.
# Imperfect solution, but good enough for now. Every time this file is imported,
#    this line runs. This pollutes the sys.path list.
sys.path.append("/home/andlars/Desktop/RS_TP_Adaptation/src")

import numpy as np

import core.machining.machining as mach
import util.geom as geom
import core.part.part as part
import core.tool_pass.tool_pass as tp
import core.boundary_conditions.boundary_conditions as bc
import core.material_properties.material_properties as mp


if __name__ == "__main__":


    # ----- Specifying the initial geometry -----

    material = mp.ElasticMaterial(.3, 10**(9))
    path_to_cae = "/home/andlars/Desktop/RS_TP_Adaptation/experiments/experiments/test_initial_geometry_cae/test_initial_geometry.cae" 
    abaqus_part = part.AbaqusDefinedPart("an_example_part", path_to_cae, material)


    # ----- Specifying the clamping setup (aka the boundary conditions) -----

    # Clamp on one side of bar.
    v1 = geom.Point3D(np.array([0, 10, 0]))
    v2 = geom.Point3D(np.array([40, 10, 0]))
    v3 = geom.Point3D(np.array([40, 10, 40]))
    v4 = geom.Point3D(np.array([0, 10, 40]))
    clamp_surface_vertices = [v1, v2, v3, v4]
    clamp_surface1 = geom.NGon3D(clamp_surface_vertices)

    # Clamp on other side of bar.
    v1 = geom.Point3D(np.array([0, 10, 400]))
    v2 = geom.Point3D(np.array([0, 10, 360]))
    v3 = geom.Point3D(np.array([40, 10, 360]))
    v4 = geom.Point3D(np.array([40, 10, 400]))
    clamp_surface_vertices = [v1, v2, v3, v4]
    clamp_surface2 = geom.NGon3D(clamp_surface_vertices)

    # Approximate how clamps restrict part movement. 
    BC_settings = bc.DisplacementBCSettings(True, True, True, True, True, True)

    BC1 = bc.SurfaceBC(clamp_surface1, BC_settings)
    BC2 = bc.SurfaceBC(clamp_surface2, BC_settings)

    BCs = [BC1, BC2]


    # ----- Building the top-level machining object -----

    machining_process = mach.MachiningProcess(abaqus_part, BCs)


    # ----- Specifying the stress profile for the first commitment phase -----

    # The stress profile comes from a user subroutine.
    # We associate the stress profile with the part, then when a simulation is
    #    done using the part, the user subroutine is automatically invoked to
    #    imbue the stress profile. The user subroutine is inherently invoked at
    #    simulation runtime, so only setup can be done defore that.
    path_to_subroutine = "/home/andlars/Desktop/RS_TP_Adaptation/src/core/user_subroutines/def_stress-std.o"
    machining_process.record_estimated_stress_profile(path_to_subroutine)

    
    # ----- Specifying the first tool pass -----

    # Tool Pass #1
    v1 = geom.Point3D(np.array([40, 9, 215]))
    v2 = geom.Point3D(np.array([40, 20, 215]))
    v3 = geom.Point3D(np.array([40, 9, 185]))
    v4 = geom.Point3D(np.array([40, 20, 185]))
    v5 = geom.Point3D(np.array([0, 9, 215]))
    v6 = geom.Point3D(np.array([0, 20, 215]))
    v7 = geom.Point3D(np.array([0, 9, 185]))
    v8 = geom.Point3D(np.array([0, 20, 185]))

    tp1_shape = geom.SpecRightRectPrism(v1, v2, v3, v4, v5, v6, v7, v8)

    tp1 = tp.ToolPass(tp1_shape)


    # ----- Specifying the second tool pass -----

    v1 = geom.Point3D(np.array([40, 5, 40]))
    v2 = geom.Point3D(np.array([40, 5, 80]))
    v3 = geom.Point3D(np.array([40, 30, 40]))
    v4 = geom.Point3D(np.array([40, 30, 80]))
    v5 = geom.Point3D(np.array([0, 5, 40]))
    v6 = geom.Point3D(np.array([0, 5, 80]))
    v7 = geom.Point3D(np.array([0, 30, 40]))
    v8 = geom.Point3D(np.array([0, 30, 80]))

    tp2_shape = geom.SpecRightRectPrism(v1, v2, v3, v4, v5, v6, v7, v8)

    tp2 = tp.ToolPass(tp2_shape)

    # ----- Construct the tool pass plan -----

    tool_passes = [tp1, tp2]
    tool_pass_plan = tp.ToolPassPlan(tool_passes)


    # ----- Simulate the potential tool passes ----- 

    machining_process.sim_potential_tool_passes(tool_pass_plan, "potential_tool_passes")





    # ----- Specifying the first tool pass -----

    # Tool Pass #1
    v1 = geom.Point3D(np.array([40, 9, 215]))
    v2 = geom.Point3D(np.array([40, 20, 215]))
    v3 = geom.Point3D(np.array([40, 9, 185]))
    v4 = geom.Point3D(np.array([40, 20, 185]))
    v5 = geom.Point3D(np.array([0, 9, 215]))
    v6 = geom.Point3D(np.array([0, 20, 215]))
    v7 = geom.Point3D(np.array([0, 9, 185]))
    v8 = geom.Point3D(np.array([0, 20, 185]))

    tp1_shape = geom.SpecRightRectPrism(v1, v2, v3, v4, v5, v6, v7, v8)

    tp1 = tp.ToolPass(tp1_shape)


    # ----- Specifying the second tool pass -----

    v1 = geom.Point3D(np.array([40, 5, 40]))
    v2 = geom.Point3D(np.array([40, 5, 80]))
    v3 = geom.Point3D(np.array([40, 30, 40]))
    v4 = geom.Point3D(np.array([40, 30, 80]))
    v5 = geom.Point3D(np.array([0, 5, 40]))
    v6 = geom.Point3D(np.array([0, 5, 80]))
    v7 = geom.Point3D(np.array([0, 30, 40]))
    v8 = geom.Point3D(np.array([0, 30, 80]))

    tp2_shape = geom.SpecRightRectPrism(v1, v2, v3, v4, v5, v6, v7, v8)

    tp2 = tp.ToolPass(tp2_shape)


    # ----- Construct the tool pass plan -----

    tool_passes = [tp1, tp2]
    tool_pass_plan = tp.ToolPassPlan(tool_passes)


    # ----- Commit to some tool passes ----- 

    machining_process.commit_tool_passes(tool_pass_plan, "commit_1")




    # ----- Specifying the first tool pass -----

    # Tool Pass #1
    v1 = geom.Point3D(np.array([15, 7, 100]))
    v2 = geom.Point3D(np.array([45, 7, 100]))
    v3 = geom.Point3D(np.array([15, 7, 150]))
    v4 = geom.Point3D(np.array([45, 7, 150]))
    v5 = geom.Point3D(np.array([15, 30, 100]))
    v6 = geom.Point3D(np.array([45, 30, 100]))
    v7 = geom.Point3D(np.array([15, 30, 150]))
    v8 = geom.Point3D(np.array([45, 30, 150]))

    tp1_shape = geom.SpecRightRectPrism(v1, v2, v3, v4, v5, v6, v7, v8)

    tp1 = tp.ToolPass(tp1_shape)

    # ----- Construct the tool pass plan -----

    tool_passes = [tp1]
    tool_pass_plan = tp.ToolPassPlan(tool_passes)


    # ----- Simulate the potential tool passes ----- 

    machining_process.sim_potential_tool_passes(tool_pass_plan, "after_commit_1")




    # ----- Specifying the first tool pass -----

    # Tool Pass #1
    v1 = geom.Point3D(np.array([10, 5, 300]))
    v2 = geom.Point3D(np.array([10, 5, 350]))
    v3 = geom.Point3D(np.array([30, 5, 300]))
    v4 = geom.Point3D(np.array([30, 5, 350]))
    v5 = geom.Point3D(np.array([10, 40, 300]))
    v6 = geom.Point3D(np.array([10, 40, 350]))
    v7 = geom.Point3D(np.array([30, 40, 300]))
    v8 = geom.Point3D(np.array([30, 40, 350]))

    tp1_shape = geom.SpecRightRectPrism(v1, v2, v3, v4, v5, v6, v7, v8)

    tp1 = tp.ToolPass(tp1_shape)


    # ----- Specifying the second tool pass -----

    # Tool Pass #2
    v1 = geom.Point3D(np.array([15, 7, 100]))
    v2 = geom.Point3D(np.array([45, 7, 100]))
    v3 = geom.Point3D(np.array([15, 7, 150]))
    v4 = geom.Point3D(np.array([45, 7, 150]))
    v5 = geom.Point3D(np.array([15, 30, 100]))
    v6 = geom.Point3D(np.array([45, 30, 100]))
    v7 = geom.Point3D(np.array([15, 30, 150]))
    v8 = geom.Point3D(np.array([45, 30, 150]))

    tp2_shape = geom.SpecRightRectPrism(v1, v2, v3, v4, v5, v6, v7, v8)

    tp2 = tp.ToolPass(tp2_shape)

    # ----- Construct the tool pass plan -----

    tool_passes = [tp1, tp2]
    tool_pass_plan = tp.ToolPassPlan(tool_passes)


    # ----- Commit to some tool passes ----- 

    machining_process.commit_tool_passes(tool_pass_plan, "commit_2")