"""
This file is a testbed for the library in development. If you wish to use the
   library, this file demonstrates how to use it. 
"""



# *******
# Necessary *modifications* for Abaqus PDE debugging.

# Make importing work.
import sys
sys.path.append("/home/andlars/Desktop/RS_TP_Adaptation")

# Workaround for bug in Abaqus 2024 that makes .o files produced from the
#     "Abaqus make" utility unusable. 
import os
os.environ["ABA_DISABLE_DSL_USUB_CHECK"] = "1"

# *******


import numpy as np

import src.core.machining.machining as mach
import src.util.geom as geom
import src.core.part.part as part
import src.core.tool_pass.tool_pass as tp
import src.core.boundary_conditions.boundary_conditions as bc
import src.core.material_properties.material_properties as mp

import src.util.general_util as util
from src.util.debug import *



if __name__ == "__main__":

    try:

        # ************
        #   FS Setup
        # ************

        # ----- Main Sim: FS Paths ----- 
        INIT_STL = "/home/andlars/Desktop/RS_TP_Adaptation/experiments/experiments/full_flow/initial_geometry.stl" 
        INIT_CAE_DIR = "/home/andlars/Desktop/RS_TP_Adaptation/experiments/experiments/full_flow/initial_geom/"
        INIT_CAE_NAME = "initial_geometry.cae" 
        TP_DIR = "/home/andlars/Desktop/RS_TP_Adaptation/experiments/experiments/full_flow/stress_estimation/"
        TP_STRESS_ESTIMATION_DIR = "/home/andlars/Desktop/RS_TP_Adaptation/experiments/experiments/full_flow/stress_estimation/recover_stresses/"
        potential_tpp_names = ["first_plan", "second_plan", "third_plan"]
        committed_tpp_names = ["first_committed_plan", "second_committed_plan"]
        STRESS_FIELD_ESTIMATE = "/home/andlars/Desktop/RS_TP_Adaptation/experiments/experiments/full_flow/stress_fields/stress_field_estimate/stress_field_estimate-std.o"

        # ----- Main Sim: Directory Cleanup -----
        util.nuke_and_remake(INIT_CAE_DIR)
        util.nuke_and_remake(TP_DIR)
        util.nuke_and_remake(TP_STRESS_ESTIMATION_DIR)

        # ----- Real Life Sim: FS Paths -----
        INIT_CAE_DIR_REAL_LIFE = "/home/andlars/Desktop/RS_TP_Adaptation/experiments/experiments/full_flow/in_real_life/initial_geom/"
        INIT_CAE_NAME_REAL_LIFE = INIT_CAE_NAME
        PARALLEL_TP_DIR = "/home/andlars/Desktop/RS_TP_Adaptation/experiments/experiments/full_flow/in_real_life/tool_passes/"
        ACTUAL_STRESS_FIELD = "/home/andlars/Desktop/RS_TP_Adaptation/experiments/experiments/full_flow/stress_fields/actual_stress_field/actual_stress_field-std.o"
        real_life_cae_names = committed_tpp_names

        # ----- Real Life Sim: Directory Cleanup -----
        util.nuke_and_remake(INIT_CAE_DIR_REAL_LIFE)
        util.nuke_and_remake(PARALLEL_TP_DIR)


        # **********************
        #   Simulation Setup 
        # **********************
        
        # ----- Both: Material, Boundary Conditions ----- 

        material = mp.ElasticMaterial(.3, 10**(9))

        # Clamp on one side of bar.
        v1 = geom.Point3D(np.array([0, 10, 0]))
        v2 = geom.Point3D(np.array([40, 10, 0]))
        v3 = geom.Point3D(np.array([0, 0, 0]))
        v4 = geom.Point3D(np.array([40, 0, 0]))
        clamp_surface_vertices1 = [v1, v2, v3, v4]
        # clamp_surface1 = geom.NGon3D(clamp_surface_vertices1)

        # Clamp on other side of bar.
        v1 = geom.Point3D(np.array([0, 10, 400]))
        v2 = geom.Point3D(np.array([40, 10, 400]))
        v3 = geom.Point3D(np.array([0, 0, 400]))
        v4 = geom.Point3D(np.array([40, 0, 400]))
        clamp_surface_vertices2 = [v1, v2, v3, v4]
        # clamp_surface2 = geom.NGon3D(clamp_surface_vertices2)

        # Approximate how clamps restrict part movement. 
        BC_settings = bc.DisplacementBCSettings(True, True, True, True, True, True)
        
        # Surface boundary conditions are buggy right now.
        # BC1 = bc.SurfaceBC(clamp_surface1, BC_settings)
        # BC2 = bc.SurfaceBC(clamp_surface2, BC_settings)
        # BCs = [BC1, BC2]

        # Vertex boundary conditions.
        BC1 = bc.VertexBC(clamp_surface_vertices1, BC_settings) 
        BC2 = bc.VertexBC(clamp_surface_vertices2, BC_settings) 
        BCs = [BC1, BC2]

        # ----- Real Life Sim: Top Level Machining Object -----
        real_life_machining = mach.MachiningProcess(INIT_STL, BCs, material, os.path.join(INIT_CAE_DIR_REAL_LIFE, INIT_CAE_NAME))

        # ----- Main Sim: Top Level Machining Object -----
        main_machining = mach.MachiningProcess(INIT_STL, BCs, material, os.path.join(INIT_CAE_DIR, INIT_CAE_NAME))

        # ----- Real Life Sim: Specifying Actual Stress Profile ----- 
        real_life_machining.use_stress_profile(ACTUAL_STRESS_FIELD)

        # ----- Main Sim: Specifying Initial Estimated Stress Profile -----
        main_machining.use_stress_profile(STRESS_FIELD_ESTIMATE)

        
        # ********************************
        #        First Committed TPP 
        # ********************************

        # ----- Main Sim: Simulate First Potential TPP -----
        p1 = geom.Point3D(np.array((-10, 5, 100)))
        p2 = geom.Point3D(np.array((0, 5, 90)))
        p3 = geom.Point3D(np.array((10, 5, 120)))
        p4 = geom.Point3D(np.array((30, 5, 150)))
        path = geom.PlanarCubicC2Spline3D([p1, p2, p3, p4])
        tp1 = tp.ToolPass(path, 1, 10)

        potential_tpp_1 = tp.ToolPassPlan([tp1])

        main_machining.sim_potential_tool_passes(potential_tpp_1, potential_tpp_names[0], TP_DIR)

        # ----- Main Sim: Simulate Second Potential TPP -----
        # p1 = geom.Point3D(np.array((5, 5, 300)))
        # p2 = geom.Point3D(np.array((5, 5, 200)))
        # path = geom.PlanarCubicC2Spline3D([p1, p2])
        # tp1 = tp.ToolPass(path, 2, 10)
        # 
        # p1 = geom.Point3D(np.array((5, 5, 300)))
        # p2 = geom.Point3D(np.array((35, 5, 300)))
        # path = geom.PlanarCubicC2Spline3D([p1, p2])
        # tp2 = tp.ToolPass(path, 2, 10)

        # potential_tpp_2 = tp.ToolPassPlan([tp1, tp2])

        # main_machining.sim_potential_tool_passes(potential_tpp_2, potential_tpp_names[1], TP_DIR)

        # ----- Main Sim: Commit to a TPP -----
        committed_plan = potential_tpp_1
        main_machining.commit_tool_passes(committed_plan, committed_tpp_names[0], TP_DIR)

        # ----- Real Life Sim: Simulate Committed TPP -----

        # TODO: Need to increase mesh density to better simulate real life. 
        real_life_machining.commit_tool_passes(committed_plan, committed_tpp_names[0], PARALLEL_TP_DIR)
        odb_path = os.path.join(os.path.join(PARALLEL_TP_DIR, committed_tpp_names[0]), committed_tpp_names[0] + ".odb")

        # ----- Main Sim: Pass in the results from real life -----
        main_machining.add_real_world_machining_data(odb_path)

        # ----- Main Sim: Estimate Stresses -----
        _ = main_machining.estimate_stress(TP_STRESS_ESTIMATION_DIR)
        
        """
        # ----- Main Sim: Provide Estimate of Whole Stress Field -----
        main_machining.use_stress_profile(STRESS_FIELD_ESTIMATE)

        
        # *******************************
        #   Searching For 2nd Good TPP 
        # *******************************

        # ----- Main Sim: Simulate Potential TPP -----
        p1 = geom.Point3D(np.array((20, 5, 100)))
        p2 = geom.Point3D(np.array((20, 5, 300)))
        path = geom.PlanarCubicC2Spline3D([p1, p2])
        tp3 = tp.ToolPass(path, 3, 10)

        potential_tpp_3 = tp.ToolPassPlan([tp3])

        main_machining.sim_potential_tool_passes(potential_tpp_3, potential_tpp_names[2], TP_DIR)

        # ----- Main Sim: Commit to a TPP -----
        committed_plan = potential_tpp_3
        main_machining.commit_tool_passes(committed_plan, committed_tpp_names[1], TP_DIR) 
        """

    except BaseException as e:
        dump_exception()       
