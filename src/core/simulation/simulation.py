"""
Provides top-level functionality to do simulations of interest in Abaqus.
"""

import os
from typing import Optional, Any
import pathlib
import random
import shutil
import copy

from abaqus import *
from abaqusConstants import * 
import numpy as np

import src.core.abaqus.abaqus_shim as shim
import src.core.tool_pass.tool_pass as tp
import src.core.boundary_conditions.boundary_conditions as bc
import src.core.metadata.metadata as md
import src.core.metadata.naming as naming
import src.core.residual_stress.residual_stress as rs
import src.core.metadata.abaqus_metadata as abq_md
import src.util.geom as geom

from src.util.debug import *        


def sim_tool_pass_plan(tool_pass_plan: tp.ToolPassPlan, name: str, target_dir: str, 
                       commitment_phase_md: md.CommitmentPhaseMetadata, 
                       stress_subroutine: Optional[str] = None) -> None:
    """Simulates consecutive tool passes and saves off the results.
    
       Args:
           tool_pass_plan:      The tool pass plan to simulate. 
           name:                Desired name of the MDB which results from all the tool 
                                    path simulations. Also, the desired name of the
                                    subdirectory containing the simulation artifacts. 
           target_dir:          Absolute path to directory. In this directory, a 
                                    subdirectory is created and named name. The
                                    simulation artifacts (.cae, .odb, .msg, etc.)
                                    are generated in this subdirectory.
                                If a subdirectory with name name already exists
                                    in this directory, no new subdirectory is created.
                                    In this case, this pre-existing subdirectory
                                    is used.
           commitment_phase_md: The metadata associated with the current commitment
                                    phase.
           stress_subroutine:   Abolsute path to stress subroutine. If this argument is not 
                                    specified, then the initial stress state of the 
                                    part will be sourced from the last simulation in 
                                    the previous commitment phase. If this argument 
                                    is specified, it is used instead.
       
       Returns:
           None.

       Raises:
           RuntimeError: The path of the new subdirectory is already occupied,
                             but it's not occupied by a directory!
    """

    # When this function is called in the context of a new commitment phase,
    #     which is not the very first one, no new subdirectory needs to 
    #     be created. One must already exist and it will be populated by the 
    #     .cae file which contains the information propagated from the last 
    #     commitment phase.
    new_dir_path = os.path.join(target_dir, name)
    fs_path = pathlib.Path(new_dir_path)
    if len(commitment_phase_md.potential_tpps) == 0 and \
       not commitment_phase_md.first_commitment_phase:
        assert fs_path.exists() and fs_path.is_dir()
    else:
        os.mkdir(new_dir_path)
            
    # And set the CWD to this directory. Now all simulation artifacts will be placed
    #    in this directory.
    cwd = os.getcwd()
    os.chdir(new_dir_path)
    
    # The MDB is located in this new directory.
    new_mdb_path = os.path.join(new_dir_path, name)
    
    mdb = shim.use_mdb(commitment_phase_md.init_part.path_to_mdb)

    for tool_pass in tool_pass_plan.plan:
        _sim_single_tool_pass(tool_pass, commitment_phase_md, mdb, stress_subroutine)    
        tool_pass_plan.pop()

    shim.save_mdb_as(new_mdb_path, mdb)
    shim.close_mdb(mdb)

    # Don't permanently change the CWD.
    os.chdir(cwd)



def _sim_single_tool_pass(tool_pass: tp.ToolPass, commit_metadata: md.CommitmentPhaseMetadata, 
                          mdb: Any, stress_subroutine: Optional[str] = None) -> None:
    """Simulates a single tool pass. This tool pass may be the first tool pass
           in a tool pass plan, or the nth tool pass in a tool pass plan.

       Args:
           tool_pass:         The tool path to simulate.
           commit_metadata:   The record associated with the search for a good next 
                                  tool pass to do.
           mdb:               Abaqus MDB object. The MDB in which the tool pass
                                  will be simulated. The MDB may already contain
                                  simulated tool passes.
           stress_subroutine: Path to stress subroutine. If this argument is not 
                                  specified, then the initial stress state of the 
                                  part will be sourced from the last simulation in 
                                  the previous commitment phase. If this argument 
                                  is specified, it is used instead.

       Returns:
           None.

       Raises:
           AssertionError: The MDB is not in a particular recognized state.
    """

    # The way that the next tool pass is simulated depends on the content of
    #    the MDB.
    # If the MDB contains a single model which represents an initial
    #    geometry, then we simulate the first tool pass. 
    # If the MDB contains n (where n > 1) models and the last model to be 
    #    added contains an orphan mesh, then we're simulating an nth tool pass
    #    (where n > 1).
    mdb_metadata = commit_metadata.potential_tpp_mdb_metadata[-1]
    last_model_name = mdb_metadata.model_names[-1]

    if shim.check_basic_geom(False, mdb):
        if stress_subroutine is not None:
            _sim_first_tool_pass(tool_pass, commit_metadata, mdb, stress_subroutine)
        else:
            _sim_first_tool_pass(tool_pass, commit_metadata, mdb)  
    elif shim.check_multiple_steps(False, last_model_name, mdb):
        _sim_nth_tool_pass(tool_pass, commit_metadata, mdb)
    else:
        raise AssertionError("Can't figure out how to do next tool pass...")



def _sim_first_tool_pass(tool_pass: tp.ToolPass, commit_metadata: md.CommitmentPhaseMetadata, 
                         mdb: Any, stress_subroutine: Optional[str] = None) -> None:
    """Simulates the first tool pass in an MDB. The MDB is assumed to not contain
           any tool pass simulations.
       
       Args:
           tool_pass:         The tool path to simulate.
           commit_metadata:   The record associated with the search for a good next 
                                  tool pass to do.
           mdb:               Abaqus MDB object. The MDB in which the tool pass
                                  will be simulated. 
           stress_subroutine: Path to stress subroutine. If this argument is not 
                                  specified, then the initial stress state of the 
                                  part will be sourced from the last simulation in 
                                  the previous commitment phase. If this argument 
                                  is specified, it is used instead.

       Returns:
           None.

       Raises:
           RuntimeError: The MDB is not in the expected state (e.g. it does not
                             have the expected content).
    """

    if not shim.check_basic_geom(False, mdb):
        raise RuntimeError("The MDB is not in the expected state.")

    mdb_metadata = commit_metadata.potential_tpp_mdb_metadata[-1]

    names = naming.ModelNames(naming.ModelTypes.FIRST_TOOL_PASS_IN_MDB, mdb_metadata)

    # Do material and section creation.
    shim.create_material(commit_metadata.init_part.material, names.new_model_name, mdb)
    shim.create_section(names.new_model_name, mdb)
    shim.assign_section_to_whole_part(names.pre_tool_pass_part_name, names.new_model_name, mdb)

    _do_boilerplate_tp_sim_ops(tool_pass, names, mdb_metadata, commit_metadata, mdb)

    if stress_subroutine is not None:
        # Since this modifies the input file directly, this should be the last
        #    thing that happens before the job is submitted and runs.
        shim.inp_add_stress_subroutine(names.new_model_name, mdb)
    else:
        # Use the stress state from the last toolpath in the previous commit.
        path_sim_file = commit_metadata.path_last_commit_sim_file

        assert path_sim_file is not None

        # Since this modifies the input file directly, this should be the last
        #    thing that happens before the job is submitted and runs.
        shim.inp_map_stress(path_sim_file, names.new_model_name, mdb)

    job = shim.create_job(names.new_model_name, mdb_metadata, mdb) 

    if stress_subroutine is not None:
        # Associate the user subroutine with the job.
        shim.add_user_subroutine(job, stress_subroutine)

    shim.run_job(job)



def _sim_nth_tool_pass(tool_pass: tp.ToolPass, commit_metadata: md.CommitmentPhaseMetadata, 
                      mdb: Any) -> None:
    """Simulates the nth tool pass in an MDB. The MDB is assumed to already
           contain at least one tool pass simulation.

       This function does not accept a stress subroutine because it uses the
           stress state which resulted from the last tool pass simulation as
           the initial stress state for the current simulation.

       Args:
           tool_pass:         The tool pass to simulate.
           commit_metadata:   The record associated with the search for a good next 
                                  tool pass to do.
           mdb:               Abaqus MDB object. The MDB in which the tool pass
                                  will be simulated. 

       Returns:
           None.

       Raises:
           None.
    """

    mdb_metadata = commit_metadata.potential_tpp_mdb_metadata[-1]

    names = naming.ModelNames(naming.ModelTypes.NTH_TOOL_PASS_IN_MDB, mdb_metadata)
    last_odb_file_name = naming.last_odb_file_name(mdb_metadata)

    # The .sim file contains the stress information which resulted from the last
    #    simulation.
    last_sim_file_name = naming.last_sim_file_name(mdb_metadata)

    # Create the new model with the deformed part in it from the ODB. 
    shim.create_model_and_part_from_odb(names.pre_tool_pass_part_name, names.new_model_name, last_odb_file_name, mdb_metadata, mdb)

    # Do material and section creation.
    shim.create_material(commit_metadata.init_part.material, names.new_model_name, mdb)
    shim.create_section(names.new_model_name, mdb)

    _orphan_mesh_to_geometry(names.pre_tool_pass_part_name, names.new_model_name, mdb)

    # Do the section assignment only after the orphan mesh has been mapped to a
    #    geometry. This cannot be done before mapping to a geometry. 
    shim.assign_section_to_whole_part(names.pre_tool_pass_part_name, names.new_model_name, mdb)

    _do_boilerplate_tp_sim_ops(tool_pass, names, mdb_metadata, commit_metadata, mdb)

    # Map the stress profile which existed at the end of the last simulation
    #    onto the new geometry. 
    # Since this modifies the input file directly, this should be the last
    #    thing that happens before the job is submitted and runs.
    shim.inp_map_stress(last_sim_file_name, names.new_model_name, mdb)

    # Create the job and submit it.
    job = shim.create_job(names.new_model_name, mdb_metadata, mdb) 
    shim.run_job(job)



def _do_boilerplate_tp_sim_ops(tool_pass: tp.ToolPass, 
                               names: naming.ModelNames, 
                               mdb_metadata: abq_md.AbaqusMdbMetadata,
                               commit_metadata: md.CommitmentPhaseMetadata, 
                               mdb: Any) -> None:
    """Does some operations such as instancing, section assignment, etc. which
           are boilerplate (i.e. they need to be done just about any tool pass
           simulation).
       
       Args:
           tool_pass:       The tool pass to simulate.
           names:           The names associated with the target model in the
                                MDB.
           mdb_metadata:    Metadata asscoaited with the the MDB.
           commit_metadata: Metadata associated with the commitment phase.
           mdb:             Abaqus MDB object. The MDB in which the tool pass
                                will be simulated. 

       Returns:
           None.

       Raises:
           None.
    """

    initial_geom_part = mdb.models[names.new_model_name].parts[names.pre_tool_pass_part_name]
    initial_geom_instance = shim.instance_part_into_assembly(names.pre_tool_pass_part_name, initial_geom_part, False, names.new_model_name, mdb)  

    tool_pass_part = shim.build_tool_pass_part(names.tool_pass_part_name, tool_pass, names.new_model_name, mdb_metadata, mdb)
    tool_pass_part_instance = shim.instance_part_into_assembly(names.tool_pass_part_name, tool_pass_part, False, names.new_model_name, mdb)

    cutting_instances = (tool_pass_part_instance, )
    post_tool_pass_part = shim.cut_instances_in_assembly(names.post_tool_pass_part_name, initial_geom_instance, cutting_instances, names.new_model_name, mdb_metadata, mdb)

    shim.assign_only_section_to_part(post_tool_pass_part, names.new_model_name, mdb)

    # Simplify the part geometry using some heuristics.
    # This cannot be done earlier because parts with virtual topologies cannot be used in boolean operations...
    shim.add_virtual_topology(post_tool_pass_part)

    # Instance the post cut geometry.
    # The instance is independent so that meshing can be done in the Assembly
    #    module.
    post_tool_pass_instance = shim.instance_part_into_assembly(names.post_tool_pass_part_name, post_tool_pass_part, False, names.new_model_name, mdb)

    bc.apply_BCs(commit_metadata.BCs, shim.STANDARD_INITIAL_STEP_NAME, post_tool_pass_instance, names.new_model_name, mdb)

    # Then mesh the part in the assembly module.
    shim.mesh(post_tool_pass_instance, 20, names.new_model_name, mdb)

    last_step_name = mdb_metadata.models_metadata[names.new_model_name].step_names[-1]
    shim.create_step(names.equil_step_name, last_step_name, names.new_model_name, mdb_metadata, mdb)



def _simulate_traction_app(traction: geom.Vec3D,
                           point_near_face: geom.Point3D, 
                           model_to_copy: str,
                           mdb_metadata: abq_md.AbaqusMdbMetadata,
                           commit_phase_md: md.CommitmentPhaseMetadata,
                           mdb: Any) -> str:
    """Does some operations such as instancing, section assignment, traction 
           creation, meshing, etc. and then runs the job.

       Implementation Detail:
       By assuming that the model to copy already has a geometry associated with
           it, there is no need to map an orphan mesh to a geometry in this
           function, which yields significant runtime savings.
       
       Args:
           traction:        The traction vector to apply.
           point_near_face: A point which lies on, or very near, the face on which
                                to apply the traction.
           model_to_copy:   The name of the model to copy. This model should
                                define a part geometry.
           mdb_metadata:    Metadata associated with the MDB. 
           commit_phase_md: Metadata associated with the commitment phase.
           mdb:             Abaqus MDB object. The MDB to use.
    
       Returns:
           The name of the job which was run. This is not an absolute path.
    
       Raises:
           None.
    """

    if not shim.check_standard_model(True, model_to_copy, mdb):
        raise RuntimeError("The model to copy is in an unexpected state.")

    names = naming.ModelNames(naming.ModelTypes.TRACTION_APP, mdb_metadata)
    
    shim.copy_model(model_to_copy, names.new_model_name, mdb_metadata, mdb)

    initial_geom_part = mdb.models[names.new_model_name].parts[names.deformed_part_name]

    shim.create_material(commit_phase_md.init_part.material, names.new_model_name, mdb)
    shim.create_section(names.new_model_name, mdb)
    shim.assign_section_to_whole_part(names.deformed_part_name, names.new_model_name, mdb)

    shim.assign_only_section_to_part(initial_geom_part, names.new_model_name, mdb)

    shim.add_virtual_topology(initial_geom_part)

    # The instance is independent so that meshing can be done in the Assembly
    #    module.
    initial_geom_instance = shim.instance_part_into_assembly(names.deformed_part_name, 
                                                             initial_geom_part, 
                                                             False, 
                                                             names.new_model_name, 
                                                             mdb)  

    bc.apply_BCs(commit_phase_md.BCs, shim.STANDARD_INITIAL_STEP_NAME, 
                 initial_geom_instance, names.new_model_name, mdb)
    
    # Add a step for traction application. 
    model_metadata = mdb_metadata.models_metadata[names.new_model_name]
    last_step_name = model_metadata.step_names[-1]
    shim.create_step(names.traction_step_name, last_step_name, names.new_model_name, 
                     mdb_metadata, mdb)
    
    face = shim.get_closest_face([point_near_face], initial_geom_instance, True)

    shim.create_traction(shim.STANDARD_SURFACE_TRACTION_NAME, names.traction_step_name,
                         traction, face, names.new_model_name, mdb) 

    shim.mesh(initial_geom_instance, 20, names.new_model_name, mdb)

    job = shim.create_job(names.new_model_name, mdb_metadata, mdb)

    shim.run_job(job)

    job_name = mdb_metadata.models_metadata[names.new_model_name].job_name

    return job_name 



def create_mdb_from_odb(new_mdb_name: str, new_mdb_path: str, odb_path: str
                       ) -> tuple[Any, abq_md.AbaqusMdbMetadata, naming.ModelNames]:
    """Creates an MDB with a part in it by sourcing the content of an ODB.
           Saves the MDB before returning.
      
       Implementation Detail:
       The MDB is populated with a single part and that part is named in a
           generic way.
        
       Args:
           new_mdb_name: The desired name of the new MDB.
           new_mdb_path: Absolute path to directory. The MDB (.cae file) will be
                             created in this directory.
           odb_path:     Absolute path to ODB (.odb file).
     
       Returns:
           Tuple containing the new Abaqus MDB object, the metadata that accompanies 
               it, and the names associated with it.
    
       Raises:
           None.
    """

    # Create the new MDB.
    mdb = shim.create_mdb(new_mdb_name, new_mdb_path)
    full_path = os.path.join(new_mdb_path, new_mdb_name)
    mdb_metadata = abq_md.AbaqusMdbMetadata(full_path) 
    
    # Generate the names for the stuff in the MDB.
    names = naming.ModelNames(naming.ModelTypes.ODB_TO_MDB, mdb_metadata)

    # Fetch the result of the committed tool pass and use it as the initial state
    #     of the MDB.
    shim.create_part_from_odb(names.part_from_odb_name, names.new_model_name, odb_path, mdb_metadata, mdb)
    _orphan_mesh_to_geometry(names.part_from_odb_name, names.new_model_name, mdb)
    
    # Save the MDB so it is visible in the file system.
    shim.save_mdb(mdb)

    return mdb, mdb_metadata, names



def _orphan_mesh_to_geometry(part_name: str, model_name: str, mdb: Any) -> None:
    """Converts an orphan mesh to a solid geometry by converting all the faces
           of the orphan mesh to regions/surfaces, and then melding all the
           regions/surfaces together to form a solid geometry.
           
       For some context, the result of a simulation is an orphan mesh (i.e. 
           just a bunch of vertices and elements connecting the vertices 
           together). An orphan mesh by itself is not very useful. For example, 
           to do boolean operations between parts, the parts need to have 
           geometries associated with them. This function adds geometric features 
           to a part, giving it a geometry.

       Args:
           part_name: The name of the part which contains the orphan mesh. This
                          part will have the solid geometry associated with it.
           model_name: The name of the model the part is in.
           mdb: Abaqus MDB object.

       Returns:
           None.

       Raises:
           None.
    """

    part = mdb.models[model_name].parts[part_name]

    unique_elem_faces = shim.get_unique_element_faces(part)

    # For each unique face in the mesh, we know the face is on the surface if it is
    #    associated with exactly one element.
    # We need to construct a region for each face on the surface of the part. 
    # Each region is then used to build a geometric face feature associated
    #    with the part.

    cnt = 0
    for elem_face in unique_elem_faces:
        if len(shim.get_mesh_face_elements(elem_face)) == 1:

            face_reg = shim.build_region_with_elem_face(elem_face, part)

            shim.add_face_from_region(face_reg, part)

            # Clear the cache every once in a while for a potential speed up. 
            cnt += 1
            if cnt % 50 == 0:
                part.clearGeometryCache()

    # Note: The convert_shell_to_solid() function is significantly less failure 
    #     prone than add_solid_from_faces().
    shim.convert_shell_to_solid(part)
    
    # Remove any dependency on the orphan mesh.
    # Not doing this causes the orphan mesh to still appear in the Assembly
    #    module, which can make appearence confusing. 
    shim.suppress_feature(shim.STANDARD_ORPHAN_MESH_FEATURE_NAME, part)



def estimate_residual_stresses(stress_recovery_mdb: Any,
                               stress_recovery_mdb_md: abq_md.AbaqusMdbMetadata,
                               target_geometry_mdb: Any,
                               target_geometry_mdb_md: abq_md.AbaqusMdbMetadata,
                               tool_pass: tp.ToolPass,
                               commit_phase_md: md.CommitmentPhaseMetadata,
                               path: str
                              ) -> rs.ConstantResidualStressField:
    """Estimates the residual stresses which existed in a region of material
           which was removed and caused a deformation of the workpiece.
       
       Args:
           stress_recovery_mdb:    Abaqus MDB object. MDB which contains the 
                                       geometry of the part after the committed 
                                       tool pass. Should be based on real world 
                                       data. 
           stress_recovery_mdb_md: Metadata associated with the MDB containing 
                                       the part after the committed tool pass. 
           target_geometry_mdb:    Abaqus MDB object. MDB which contains the 
                                       geometry of the part before the committed 
                                       tool pass. Should be based on real world 
                                       data.
           target_geometry_mdb_md: Metadata associated with MDB containing the 
                                       part before the committed tool pass.
           tool_pass:              The tool pass that removed material.
           commit_phase_md:        The commitment phase during which the tool 
                                       pass was committed to.
           path:                   Absolute path to directory. Any MDBs which 
                                       need to be created in order to recover 
                                       the residual stresses are placed in this 
                                       directory.
    
       Returns:
           The residual stress field in the material that was removed by the
               committed tool pass plan.
    
       Raises:
           None.
    """

    return _recover_constant_residual_stress(stress_recovery_mdb, stress_recovery_mdb_md,
                                             target_geometry_mdb, target_geometry_mdb_md,
                                             tool_pass, commit_phase_md, path)



def _recover_constant_residual_stress(stress_recovery_mdb: Any,
                                      stress_recovery_mdb_md: abq_md.AbaqusMdbMetadata,
                                      target_geometry_mdb: Any,
                                      target_geometry_mdb_md: abq_md.AbaqusMdbMetadata,
                                      tool_pass: tp.ToolPass,
                                      commit_phase_md: md.CommitmentPhaseMetadata,
                                      path: str
                                     ) -> rs.ConstantResidualStressField:
    """Recovers the residual stresses which existed in the region of material 
           which was removed by a committed tool pass. 

       This technique assumes that the residual stress field which existed in
           the region of material was constant before it was removed.
        
       Args:
           stress_recovery_mdb:    Abaqus MDB object. MDB which contains the 
                                       geometry of the part after the committed 
                                       tool pass. Should be based on real world 
                                       data. 
           stress_recovery_mdb_md: Metadata associated with the MDB containing 
                                       the part after the committed tool pass. 
           target_geometry_mdb:    Abaqus MDB object. MDB which contains the 
                                       geometry of the part before the committed 
                                       tool pass. Should be based on real world 
                                       data.
           target_geometry_mdb_md: Metadata associated with MDB containing the 
                                       part before the committed tool pass.
           tool_pass:              The tool pass that removed material.
           commit_phase_md:        The commitment phase during which the tool 
                                       pass was committed to.
           path:                   Absolute path to directory. Any MDBs which 
                                       need to be created in order to recover 
                                       the residual stresses are placed in this 
                                       directory.

       Returns:
           The residual stress field in the material that was removed by the
               committed tool pass plan.
    
       Raises:
           None.
    """
    
    # Step 1:
    # Apply the committed tool pass plan to the real world data from the previous
    #     commitment phase. The resulting geometry is that which would have been
    #     created if there were no residual stresses present.

    # Step 2:
    # Find a point on the trench for the tool pass.
    # This point will be used to locate the trench in the real-world part data. 
    # TODO: If the deformation due to the committed tool pass plan is large,
    #     this point won't actually be close to the trench. This will be an
    #     issue, but for now, I'll assume that the deformations are relatively
    #     small.
    tool_pass = commit_phase_md.committed_tpp[1].plan[0]
    trench_point = _find_trench_point(tool_pass) 

    # Step 3:
    # Apply a test traction and estimate the gradient of the volumetric difference
    #     function.
    start_traction = geom.Vec3D(np.array([1, 2, 3]))
    _, _ = _estimate_gradient_volume(start_traction, trench_point,
                                     stress_recovery_mdb, stress_recovery_mdb_md,
                                     target_geometry_mdb, target_geometry_mdb_md,
                                     commit_phase_md, path)

    # Step 4:
    # Do gradient descent. 

    # Step 5:
    # Save off the MDBs so they appear in the filesystem.
    shim.save_mdb(stress_recovery_mdb)
    shim.save_mdb(target_geometry_mdb)

    # Step 6:
    # Map the "good traction vector" to the components of the stress tensor
    #     assuming a constant state of stress.

    return None
         

     
def _find_trench_point(toolpass: tp.ToolPass) -> geom.Point3D:
    """Finds a point on the trench of a tool pass.

       Args:
           toolpass: The tool pass for which to find the trench.
    
       Returns:
           Point which is on the trench of a tool pass.
    
       Raises:
           None.
    """
    
    # Assumes an orientation of the tool with respect to the path that the tool
    #     follows.
    median_index = int(len(toolpass.path.v_list) / 2)
    point_on_trench = toolpass.path.v_list[median_index]

    return point_on_trench 



def _estimate_gradient_volume(at: geom.Vec3D, 
                              point_near_face: geom.Point3D,
                              deformed_mdb: Any,
                              deformed_mdb_md: abq_md.AbaqusMdbMetadata,
                              target_geometry_mdb: Any,
                              target_geometry_mdb_md: abq_md.AbaqusMdbMetadata,
                              commit_phase_md: md.CommitmentPhaseMetadata,
                              path: str,
                              about: geom.Vec3D = geom.Vec3D(np.array([0, 0, 0]))
                             ) -> tuple[geom.Vec3D, int]:
    """Estimates the gradient of the volume-difference function at some point
           in its domain. 

       The volume-difference function is a function with domain which is the
           set of all possible traction vectors and range which is all positive 
           reals. It is defined for a particular target geometry, a particular
           deformed geometry, and a face. To evaluate the function at a point, 
           you simply apply the traction (the point) to the deformed geometry on 
           the specified face, overlap the resulting geometry in space with the 
           target geometry, and compute the volume of the pieces which are not
           in the intersection.
                   
       Args:
           at:                     The traction vector to apply. This is the location
                                       in the domain at which the gradient is
                                       estimated.
           point_near_face:        Point very near the face on which to apply the
                                       tractions.
           deformed_mdb:           Abaqus MDB object. The MDB which contains the
                                       deformed geometry to which the traction
                                       is applied.
           deformed_mdb_md:        Metadata associated with deformed MDB.
           target_geometry_mdb:    Abaqus MDB object. The MDB which contains the
                                       target geometry (i.e. the ideal geometry
                                       after traction vector application).
           target_geometry_mdb_md: Metadata associated with the target geometry
                                       MDB.
           commit_phase_md:        Metadata associated with the commitment phase
                                       to which the deformed MDB belongs. 
           path:                   Absolute path to a directory where simulation
                                       results will be dumped.
           about:                  When estimating the gradient at a point, it's
                                       necessary to use some other point. This
                                       is the other point. Note that the estimate
                                       of the gradient is generally better when
                                       the other point which is used is close to
                                       the point of interest.
    
       Returns:
           A vector which is the estimate of the gradient and the value of the
               volume difference function evaluated at the point of interest.
    
       Raises:
           None.
    """
    
    # Change the CWD to fix location of simulation results.
    orig_cwd = os.getcwd()
    os.chdir(path)
    
    # Step 1:
    # Simulate the traction application.
    job_name = _simulate_traction_app(at, point_near_face, shim.STANDARD_MODEL_NAME,
                                      deformed_mdb_md, commit_phase_md, deformed_mdb)
    # Restore the CWD.
    os.chdir(orig_cwd)

    path_to_odb = os.path.join(path, job_name)
    
    # Step 2:
    # Make a copy of the model containing the target geometry in the target
    #     geometry MDB.

    # Generate and retrieve names.
    names = naming.ModelNames(naming.ModelTypes.TARGET_GEOM, target_geometry_mdb_md)

    shim.copy_model(names.target_model_name, names.new_model_name, 
                    target_geometry_mdb_md, target_geometry_mdb)
    
    # Step 3:
    # Map the results from the ODB to a new part in the newly copied model.
    #     Then convert the orphan mesh to a geometry.

    # Create the orphan mesh and map it to a geometry.
    shim.create_part_from_odb(names.target_part_name, names.new_model_name,
                              path_to_odb, target_geometry_mdb_md,
                              target_geometry_mdb)
    _orphan_mesh_to_geometry(names.target_part_name, names.new_model_name,
                             target_geometry_mdb)

    shim.assign_section_to_whole_part(names.target_part_name, names.new_model_name,
                                      target_geometry_mdb)

    # Step 4:
    # Instance the parts, do the boolean operations, and find the volume of the
    #     region which is not in the intersection.


    # TODO:
    return geom.Vec3D(np.array([1, 2, 3])), 10


# *****************************************************************************
#                                 DEPRECATED 
# *****************************************************************************
       
def _recover_linear_relationships(points: list[geom.Point3D], 
                                  face_point: geom.Point3D,
                                  mdb_metadata: abq_md.AbaqusMdbMetadata,
                                  commit_phase_md: md.CommitmentPhaseMetadata, 
                                  mdb: Any,
                                  path: str
                                 ) -> list[Any]:
    """DEPRECATED. Was used to recover linear relationships between applied
          tractions and nodal displacements in an effort to find good tractions. 
    
       Applies test tractions to some face of a workpiece to recover the linear 
           relationships which exist between applied traction and the displacement 
           of some points.

       Args:
           points:          The points at which to recover the linear relationships.
           face_point:      A point which is on the face on which tractions will
                                be applied. You might ask: Why not just pass in
                                the face itself? The reason is that an Abaqus Face
                                object is associated with a particular part/assembly
                                in a model. Since this function creates a model
                                for each traction which is applied, the Abaqus Face
                                object would be invalid for these new models.
           mdb_metadata:    Metadata associated with the MDB.
           commit_phase_md: Metadata associated with the commitment phase.
           mdb:             Abaqus MDB object. The MDB containing the part which has the
                                face. This MDB should contain only a basic part geometry
                                (i.e. it should contain a single model, a single part,
                                no non-default steps, etc.)
           path:            Absolute path to the directory where simulation results 
                                (.odb, etc.) will be dumped.
    
       Returns:
           List of 3 x 3 numpy arrays which describe the linear relationships
               between the applied tractions and the displacements of the
               points. 
           For each matrix M, the matrix M relates the quantities: M * t = d. 
               In this equation, t is a traction vector and d is a displacement vector.
           The number of matrices equals the number of points passed to this
               function.
           Note that these linear relationships only hold when tractions are
               applied to the face passed to this function.
    
       Raises:
           None.
    """

    if not shim.check_basic_geom(True, mdb):
        raise RuntimeError("The MDB should contain a basic geometry.")

    # Set the CWD to dictate location of simulation results.
    orig_cwd = os.getcwd()
    os.chdir(path)

    x = geom.Vec3D(np.array([1, 0, 0]))
    y = geom.Vec3D(np.array([0, 1, 0]))
    z = geom.Vec3D(np.array([0, 0, 1]))
    traction_vecs = (x, y, z)
   
    job_paths = []
    for vec in traction_vecs:
        job_name = _simulate_traction_app(vec, face_point, shim.STANDARD_MODEL_NAME, mdb_metadata, commit_phase_md, mdb)
        job_paths.append(os.path.join(path, job_name))

    os.chdir(orig_cwd)
    
    displacements = []
    for p in job_paths:
        odb_path = p + ".odb"
        points_to_nodes = shim.read_displacements(points, odb_path)
        displacements.append(points_to_nodes)

    # Useful for understanding difference between points passed in and the nodes
    #     at which the linear relationships are computed.
    for d in displacements:
        dp("For a dictionary of displacements: ")
        for key in d:
            dp("At point " + str(key) + " the closest node is " + str(d[key][0]) + " and the displacement of this node is " + str(d[key][1])) 
    
    # And compute the linear relationship for each point. 
    linear_relationships = []
    for point in points: 
        # Get the displacement of the point for each traction vector application. 
        disps = []
        for d in displacements:
            disps.append(d[point][1])

        # Use the displacements to compute the linear relationship.
        M = _compute_linear_relationship(traction_vecs, tuple(disps))
        linear_relationships.append(M)

    return linear_relationships 


    
def _compute_linear_relationship(tractions: tuple[geom.Vec3D, geom.Vec3D, geom.Vec3D],
                                 displacements: tuple[geom.Vec3D, geom.Vec3D, geom.Vec3D]
                                ) -> Any:
    """DEPRECATED. Was used to compute the linear relationships between a traction
           applied on a surface and the displacement of some points.
           
       Computes the linear relationship which relates an applied traction to
           the displacement of a point.
        
       Args:
           tractions:     Three linearly independent traction vectors.
           displacements: Three displacments of a single point. The displacements
                              produced by the applied tractions.
    
       Returns:
           Numpy Array object of size 3x3. If this matrix is called M, then it
               can be used to do M * t = d, where t is a traction vector and
               d is a displacement vector.
    
       Raises:
           None.
    """

    # Normalize the tractions.
    tractions = (tractions[0].normalize(), tractions[1].normalize(), tractions[2].normalize())

    # Check that the tractions are linearly independent.
    tractions = np.stack((tractions[0].components(), tractions[1].components(), tractions[2].components()), axis=-1)
    try:
        np.linalg.inv(tractions)
    except BaseException as e:
        raise RuntimeError("Tractions are linearly dependent!") 

    # Solve for the matrix.
    displacements = np.stack((displacements[0].components(), displacements[1].components(), displacements[2].components()), axis=-1)
    res = np.linalg.solve(tractions.T, displacements.T).T

    return res









    


