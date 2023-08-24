import core.machining.machining as mach
import util.geom as geom
import core.part.part as part
import core.tool_pass.tool_pass as tp
import core.boundary_conditions.boundary_conditions as bc

import core.abaqus.abaqus_shim as shim

import sys
import os


if __name__ == "__main__":

    # ----- Specifying the initial geometry -----
   
    path_to_cae = "/home/andlars/Desktop/RS_TP_Adaptation/software/script_testing/test_initial_geometry/test_initial_geometry.cae" 
    abaqus_part = part.AbaqusDefinedPart("an_example_part", path_to_cae)


    # ----- Specifying the clamping setup (aka the boundary conditions) -----

    # Clamp on one side of bar.
    v1 = geom.Point3D(0, 10, 0)
    v2 = geom.Point3D(40, 10, 0)
    v3 = geom.Point3D(40, 10, 40)
    v4 = geom.Point3D(0, 10, 40)
    clamp_surface_vertices = [v1, v2, v3, v4]
    clamp_surface1 = geom.NGon3D(clamp_surface_vertices)

    # DEBUG
    mdb = shim.use_mdb(path_to_cae)
    model = mdb.models["Model-1"]
    part = model.parts["Initial_Geometry"]
    bc.partition_face(clamp_surface1, "test", part, "Model-1", mdb)
    shim.close_mdb(mdb)

    """
    # Clamp on other side of bar.
    v1 = geom.Point3D(0, 10, 400)
    v2 = geom.Point3D(0, 10, 360)
    v3 = geom.Point3D(40, 10, 400)
    v4 = geom.Point3D(40, 10, 360)
    clamp_surface_vertices = [v1, v2, v3, v4]
    clamp_surface2 = geom.NGon3D(clamp_surface_vertices)

    # Approximate how clamps restrict part movement. 
    BC_settings = bc.DisplacementBCSettings(True, True, True, True, True, True)

    BC1 = bc.BC(clamp_surface1, BC_settings)
    BC2 = bc.BC(clamp_surface2, BC_settings)

    BCs = [BC1, BC2]


    # ----- Specifying the stress profile -----

    # The stress profile comes from a user subroutine.
    # We associate the stress profile with the part, then when a simulation is
    #    done using the part, the user subroutine is automatically invoked to
    #    imbue the stress profile. The user subroutine is inherently invoked at
    #    simulation runtime, so only setup can be done defore that.
    path_to_subroutine = "/home/andlars/Desktop/RS_TP_Adaptation/software/core/user_subroutines/def_stress-std.o"
    abaqus_part.add_stress_profile(path_to_subroutine)


    # ----- Building the top-level machining object -----

    machining_process = mach.MachiningProcess(abaqus_part, BCs)

    
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


    # ----- Specifying the second tool pass -----

    v1 = geom.Point3D(40, 5, 40)
    v2 = geom.Point3D(40, 5, 80)
    v3 = geom.Point3D(40, 30, 40)
    v4 = geom.Point3D(40, 30, 80)
    v5 = geom.Point3D(0, 5, 40)
    v6 = geom.Point3D(0, 5, 80)
    v7 = geom.Point3D(0, 30, 40)
    v8 = geom.Point3D(0, 30, 80)

    tp2_shape = geom.SpecRightRectPrism(v1, v2, v3, v4, v5, v6, v7, v8)

    tp2 = tp.ToolPass(tp2_shape)


    # ----- Specifying the third tool pass -----

    v1 = geom.Point3D(20, 5, 50)
    v2 = geom.Point3D(30, 5, 50)
    v3 = geom.Point3D(20, 20, 50)
    v4 = geom.Point3D(30, 20, 50)
    v5 = geom.Point3D(20, 5, 150)
    v6 = geom.Point3D(30, 5, 150)
    v7 = geom.Point3D(20, 20, 150)
    v8 = geom.Point3D(30, 30, 150)

    tp3_shape = geom.SpecRightRectPrism(v1, v2, v3, v4, v5, v6, v7, v8)

    tp3 = tp.ToolPass(tp3_shape)


    # ----- Specifying the fourth tool pass -----

    v1 = geom.Point3D(10, 3, 300)
    v2 = geom.Point3D(30, 3, 300)
    v3 = geom.Point3D(10, 3, 350)
    v4 = geom.Point3D(30, 3, 350)
    v5 = geom.Point3D(10, 30, 300)
    v6 = geom.Point3D(30, 30, 300)
    v7 = geom.Point3D(10, 30, 350)
    v8 = geom.Point3D(30, 30, 350)

    tp4_shape = geom.SpecRightRectPrism(v1, v2, v3, v4, v5, v6, v7, v8)

    tp4 = tp.ToolPass(tp4_shape)


    # ----- Construct the tool pass plan -----

    tool_passes = [tp1, tp2, tp3, tp4]
    tool_pass_plan = tp.ToolPassPlan(tool_passes)


    # ----- Simulate the potential tool passes ----- 

    save_loc = "/home/andlars/Desktop/RS_TP_Adaptation/software/script_testing/test_initial_geometry/test_post_tool_pass.cae"
    machining_process.sim_potential_tool_passes(tool_pass_plan, save_loc)
    """
