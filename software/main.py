import core.machining.machining as mach
import util.geom as geom
import core.part.part as part
import core.tool_pass.tool_pass as tp

# DEBUG
import core.simulation.simulation as sim
import core.abaqus.abaqus_shim as shim

import sys
import os


if __name__ == "__main__":

    # ----- Specifying the initial geometry -----
   
    path_to_cae = "/home/andlars/Desktop/RS_TP_Adaptation/software/script_testing/test_initial_geometry/test_initial_geometry.cae" 
    abaqus_part = part.AbaqusDefinedPart("an_example_part", path_to_cae)


    # ----- Specifying the stress profile -----

    # The stress profile comes from a user subroutine.
    # We associate the stress profile with the part, then when a simulation is
    #    done using the part, the user subroutine is automatically invoked to
    #    imbue the stress profile. The user subroutine is inherently invoked at
    #    simulation runtime, so only setup can be done defore that.
    path_to_subroutine = "/home/andlars/Desktop/RS_TP_Adaptation/software/core/user_subroutines/def_stress-std.o"
    abaqus_part.add_stress_profile(path_to_subroutine)


    # ----- Building the top-level machining object -----

    machining_process = mach.MachiningProcess(abaqus_part)

    
    # ----- Specifying the first tool pass -----

    # Tool Pass #1
    v1 = geom.Point3D(40, 9, 215)
    v2 = geom.Point3D(40, 20, 215)
    v3 = geom.Point3D(40, 9, 185)
    v4 = geom.Point3D(40, 20, 185)
    v5 = geom.Point3D(0, 9, 215)
    v6 = geom.Point3D(0, 20, 215)
    v7 = geom.Point3D(0, 9, 185)
    v8 = geom.Point3D(0, 20, 185)

    tp1_shape = geom.SpecRightRectPrism(v1, v2, v3, v4, v5, v6, v7, v8)

    tp1 = tp.ToolPass(tp1_shape)

    tool_passes = [tp1]

    tool_pass_plan = tp.ToolPassPlan(tool_passes)


    # ----- Simulating the first potential tool pass ----- 

    save_loc = "/home/andlars/Desktop/RS_TP_Adaptation/software/script_testing/test_initial_geometry/test_post_tool_pass.cae"
    machining_process.sim_next_potential_tool_pass(save_loc)

