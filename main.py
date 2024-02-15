"""
This file is a testbed for the library in development. If you wish to use the
   library, this file demonstrates how to use it. 
"""

# Necessary hack for Abaqus PDE debugging.
import sys
sys.path.append("/home/andlars/Desktop/RS_TP_Adaptation")

import pathlib
import os
import shutil

import numpy as np

import src.core.machining.machining as mach
import src.util.geom as geom
import src.core.part.part as part
import src.core.tool_pass.tool_pass as tp
import src.core.boundary_conditions.boundary_conditions as bc
import src.core.material_properties.material_properties as mp

from src.util.debug import *

if __name__ == "__main__":

    try:
        # ----- Desired names and some paths ----- 
        PATH_TO_CAE = "/home/andlars/Desktop/RS_TP_Adaptation/experiments/experiments/test_initial_geometry_cae/test_initial_geometry.cae" 
        TOOL_PASS_ROOT_DIR = "/home/andlars/Desktop/RS_TP_Adaptation/experiments/experiments/test_initial_geometry_cae/"
        TOOL_PASS_STRESS_RECOVERY_DIR = "/home/andlars/Desktop/RS_TP_Adaptation/experiments/experiments/test_initial_geometry_cae/stress_estimation"
        tool_pass_plan_names = ["first_plan", "second_plan", "committed_plan"]
        PART_NAME = "an_example_part"
        PATH_TO_SUBROUTINE = "/home/andlars/Desktop/RS_TP_Adaptation/src/core/user_subroutines/def_stress.cpp"

        # ----- Delete any directories which already exist at the desired paths. -----

        # Delete any directories matching those for tool passes.
        for name in tool_pass_plan_names:
            full_path = os.path.join(TOOL_PASS_ROOT_DIR, name) 
            fs_path = pathlib.Path(full_path)
            if fs_path.exists():
                if fs_path.is_dir():
                    shutil.rmtree(fs_path)
                else:
                    raise RuntimeError("Something exists but it's not a directory!") 

        fs_path = pathlib.Path(TOOL_PASS_STRESS_RECOVERY_DIR)
        if fs_path.exists():
            if fs_path.is_dir():
                shutil.rmtree(fs_path)
            else:
                raise RuntimeError("Something exists but it's not a directory!") 
        os.mkdir(fs_path)

        # ----- Specifying the initial geometry -----
        material = mp.ElasticMaterial(.3, 10**(9))
        abaqus_part = part.AbaqusDefinedPart(PART_NAME, PATH_TO_CAE, material)

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
        machining_process.record_estimated_stress_profile(PATH_TO_SUBROUTINE)

        # ----- First tool pass plan -----
        p1 = geom.Point3D(np.array((-10, 5, 100)))
        p2 = geom.Point3D(np.array((0, 5, 90)))
        p3 = geom.Point3D(np.array((10, 5, 120)))
        p4 = geom.Point3D(np.array((30, 5, 150)))

        path = geom.PlanarCubicC2Spline3D([p1, p2, p3, p4])
        tp1 = tp.ToolPass(path, 1, 10)
        plan = tp.ToolPassPlan([tp1])

        machining_process.sim_potential_tool_passes(plan, tool_pass_plan_names[0], TOOL_PASS_ROOT_DIR)

        # ----- Second tool pass plan -----
        """
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

        plan = tp.ToolPassPlan([tp1])

        machining_process.sim_potential_tool_passes(plan, tool_pass_plan_names[1], TOOL_PASS_ROOT_DIR)
        """

        # ----- Committing to a Plan -----
        machining_process.commit_tool_passes(plan, tool_pass_plan_names[2], TOOL_PASS_ROOT_DIR)

        # ----- Estimating Residual Stress Due to Committed Tool Pass -----
        _ = machining_process.estimate_stress_via_last_tool_pass(TOOL_PASS_STRESS_RECOVERY_DIR)

    except BaseException as e:
        dump_exception()       



