"""
This file is a testbed for the library in development. If you wish to use the
   library, this file demonstrates how to use it. 
"""

import numpy as np

import src.core.machining.machining as mach
import src.util.geom as geom
import src.core.part as part
import src.core.tool_pass.tool_pass as tp
import src.core.boundary_conditions.boundary_conditions as bc
import src.core.material_properties.material_properties as mp



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

    # ----- First tool pass plan -----
    p1 = geom.Point3D(np.array((-10, 5, 100)))
    p2 = geom.Point3D(np.array((0, 5, 90)))
    p3 = geom.Point3D(np.array((10, 5, 120)))
    p4 = geom.Point3D(np.array((30, 5, 150)))

    path = geom.PlanarCubicC2Spline3D([p1, p2, p3, p4])
    tp1 = tp.ToolPass(path, 1, 10)
    plan = tp.ToolPassPlan([tp1])

    machining_process.sim_potential_tool_passes(plan, "first_test")

    # ----- Second tool pass plan -----
    p1 = geom.Point3D(np.array((5, 5, 300)))
    p2 = geom.Point3D(np.array((5, 5, 200)))
    path = geom.PlanarCubicC2Spline3D([p1, p2])
    tp1 = tp.ToolPass(path, 2, 10)

    p1 = geom.Point3D(np.array((5, 5, 300)))
    p2 = geom.Point3D(np.array((35, 5, 300)))
    path = geom.PlanarCubicC2Spline3D([p1, p2])
    tp2 = tp.ToolPass(path, 2, 10)

    p1 = geom.Point3D(np.array((35, 5, 300)))
    p2 = geom.Point3D(np.array((35, 5, 200)))
    path = geom.PlanarCubicC2Spline3D([p1, p2])
    tp3 = tp.ToolPass(path, 2, 10)

    plan = tp.ToolPassPlan([tp1, tp2, tp3])

    machining_process.sim_potential_tool_passes(plan, "second_test")

    """
    p1 = geom.Point3D(np.array((5, 5, 200)))
    p2 = geom.Point3D(np.array((35, 5, 200)))
    path = geom.PlanarCubicC2Spline3D([p1, p2])
    tp4 = tp.ToolPass(path, 2, 10)

    p1 = geom.Point3D(np.array((5, 5, 200)))
    p2 = geom.Point3D(np.array((35, 5, 300)))
    path = geom.PlanarCubicC2Spline3D([p1, p2])
    tp5 = tp.ToolPass(path, 2, 10)
    """