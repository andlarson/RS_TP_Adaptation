"""
Provides top-level functionality to do simulations of interest in Abaqus.
"""

import os
from typing import Optional, Any
import pathlib

from abaqus import *
from abaqusConstants import * 

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
                       commitment_phase_md: md.CommittedToolPassPlanMetadata, 
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
    if len(commitment_phase_md.simulated_tool_pass_plans) == 0 and \
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
    
    mdb = shim.use_mdb(commitment_phase_md.path_initial_mdb)

    for tool_pass in tool_pass_plan.plan:
        _sim_single_tool_pass(tool_pass, commitment_phase_md, mdb, stress_subroutine)    
        tool_pass_plan.pop()

        # DEBUG
        shim.save_mdb_as(new_mdb_path, mdb)

    shim.save_mdb_as(new_mdb_path, mdb)
    shim.close_mdb(mdb)

    commitment_phase_md.simulated_tool_pass_plans.append((name, tool_pass_plan, new_dir_path))

    # Don't permanently change the CWD.
    os.chdir(cwd)



def _sim_single_tool_pass(tool_pass: tp.ToolPass, commit_metadata: md.CommittedToolPassPlanMetadata, 
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
    mdb_metadata = commit_metadata.per_mdb_metadata[-1] 
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



def _sim_first_tool_pass(tool_pass: tp.ToolPass, commit_metadata: md.CommittedToolPassPlanMetadata, 
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

    mdb_metadata = commit_metadata.per_mdb_metadata[-1] 

    names = naming.ModelNames(mdb_metadata, True)

    # Do material and section creation.
    shim.create_material(commit_metadata.init_part.material, names.new_model_name, mdb)
    shim.create_section(names.new_model_name, mdb)
    shim.assign_section_to_whole_part(names.pre_tool_pass_part_name, names.new_model_name, mdb)

    _do_boilerplate_tp_sim_ops(tool_pass, names, commit_metadata, mdb)

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



def _sim_nth_tool_pass(tool_pass: tp.ToolPass, commit_metadata: md.CommittedToolPassPlanMetadata, 
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

    mdb_metadata = commit_metadata.per_mdb_metadata[-1] 

    names = naming.ModelNames(mdb_metadata, False)
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

    _do_boilerplate_tp_sim_ops(tool_pass, names, commit_metadata, mdb)

    # Map the stress profile which existed at the end of the last simulation
    #    onto the new geometry. 
    # Since this modifies the input file directly, this should be the last
    #    thing that happens before the job is submitted and runs.
    shim.inp_map_stress(last_sim_file_name, names.new_model_name, mdb)

    # Create the job and submit it.
    job = shim.create_job(names.new_model_name, mdb_metadata, mdb) 
    shim.run_job(job)



def _do_boilerplate_tp_sim_ops(tool_pass: tp.ToolPass, names: naming.ModelNames, 
                               commit_metadata: md.CommittedToolPassPlanMetadata, 
                               mdb: Any) -> None:
    """Does some operations such as instancing, section assignment, etc. which
           are boilerplate (i.e. they need to be done just about any tool pass
           simulation).

       Args:
           tool_pass:       The tool pass to simulate.
           names:           The names associated with the MDB.
           commit_metadata: Metadata associated with the commitment phase.
           mdb:             Abaqus MDB object. The MDB in which the tool pass
                                will be simulated. 

       Returns:
           None.

       Raises:
           None.
    """

    # Instance the initial geometry part.
    initial_geom_part = mdb.models[names.new_model_name].parts[names.pre_tool_pass_part_name]
    initial_geom_instance = shim.instance_part_into_assembly(names.pre_tool_pass_part_name, initial_geom_part, False, names.new_model_name, mdb)  

    # Build the next tool pass path as a part.
    tool_pass_part = shim.build_tool_pass_part(names.tool_pass_part_name, tool_pass, names.new_model_name, commit_metadata.per_mdb_metadata[-1], mdb)
    tool_pass_part_instance = shim.instance_part_into_assembly(names.tool_pass_part_name, tool_pass_part, False, names.new_model_name, mdb)

    # Create the post tool pass geometry as a part.
    cutting_instances = (tool_pass_part_instance, )
    post_tool_pass_part = shim.cut_instances_in_assembly(names.post_tool_pass_part_name, initial_geom_instance, cutting_instances, names.new_model_name, commit_metadata.per_mdb_metadata[-1], mdb)

    # Assign a section to the part.
    shim.assign_only_section_to_part(post_tool_pass_part, names.new_model_name, mdb)

    # Simplify the part geometry using some heuristics.
    # This cannot be done earlier because parts with virtual topologies cannot be used in boolean operations...
    shim.add_virtual_topology(post_tool_pass_part)

    # Instance the post cut geometry.
    # The instance is independent so that meshing can be done in the Assembly
    #    module.
    post_tool_pass_instance = shim.instance_part_into_assembly(names.post_tool_pass_part_name, post_tool_pass_part, False, names.new_model_name, mdb)

    # Apply the boundary conditions. 
    bc.apply_BCs(commit_metadata.BCs, shim.STANDARD_INITIAL_STEP_NAME, post_tool_pass_instance, names.new_model_name, mdb)

    # Then mesh the part in the assembly module.
    shim.mesh(post_tool_pass_instance, 20, names.new_model_name, mdb)

    # Add an equilibrium step following the last step on commit_metadata.
    last_step_name = commit_metadata.per_mdb_metadata[-1].models_metadata[names.new_model_name].step_names[-1]
    shim.create_equilibrium_step(names.equil_step_name, last_step_name, names.new_model_name, commit_metadata.per_mdb_metadata[-1], mdb)



def _do_boilerplate_traction_app_sim_ops(traction: geom.Traction,
                                         face: Any,
                                         names: naming.ModelNames, 
                                         commit_metadata: md.CommittedToolPassPlanMetadata,
                                         mdb: Any) -> None:
    """Does some operations such as instancing, section assignment, traction 
           creation, meshing, etc. which are boilerplate (i.e. they need to be 
           done for just about any simulation of applied traction).
       
       When this function returns, the 
        
       Args:
           traction:        The traction vector to apply.
           face:            Abaqus Face object. The face on which to apply the traction.
           names:           The names associated with the MDB.
           commit_metadata: Metadata associated with the commitment phase.
           mdb:             Abaqus MDB object.
    
       Returns:
           None.
    
       Raises:
           None.
    """



def create_mdb_from_odb(new_mdb_name: str, new_mdb_path: str, odb_path: str
                       ) -> tuple[Any, abq_md.AbaqusMdbMetadata, naming.ModelNames]:
    """Creates an MDB with a part in it by sourcing the content of an ODB.
           Saves the MDB before returning.
        
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
    names = naming.ModelNames(mdb_metadata, True)

    # Fetch the result of the committed tool pass and use it as the initial state
    #     of the MDB.
    shim.create_part_from_odb(names.pre_tool_pass_part_name, names.new_model_name, odb_path, mdb_metadata, mdb)
    _orphan_mesh_to_geometry(names.pre_tool_pass_part_name, names.new_model_name, mdb)
    
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

    """DEPRECATED. The convert_shell_to_solid() function is significantly less
           failure prone than add_solid_from_faces().

    # Build the solid feature from the face features.
    try:
        shim.add_solid_from_faces(part)
    except AbaqusException as e:
        dp("Shell to Solid conversion failure!")
        dp("Saving off the MDB as failed_shell_to_solid.cae")
        mdb.saveAs("failed_shell_to_solid.cae")
        dp("The arguments associated with the exception are " + str(e.args))
        raise

    """

    shim.convert_shell_to_solid(part)
    
    # Remove any dependency on the orphan mesh.
    # Not doing this causes the orphan mesh to still appear in the Assembly
    #    module, which can make appearence confusing. 
    shim.suppress_feature(shim.STANDARD_ORPHAN_MESH_FEATURE_NAME, part)



def estimate_residual_stresses(commitment_phase_md: md.CommittedToolPassPlanMetadata,
                               path: str) -> rs.ConstantResidualStressField:
    """Estimates the residual stresses which existed in a region of material
           which was removed and caused a deformation of the workpiece.
       
       TODO: For now, the only source of deformation data is that which results
           from simulations. It will also be necessary to use deformation data
           which is from real life.
       
       Args:
           commitment_phase_md: The commitment phase of interest. This commitment
                                    phase should contain a committed tool pass
                                    plan.
           path:                Absolute path to directory. Any MDBs which need
                                    to be created in order to recover the residual
                                    stresses are placed in this directory.
    
       Returns:
           The residual stress field in the material that was removed by the
               committed tool pass plan.
    
       Raises:
           None.
    """

    return recover_constant_residual_stress(commitment_phase_md, path)




def recover_constant_residual_stress(commitment_phase_md: md.CommittedToolPassPlanMetadata,
                                     path: str) -> rs.ConstantResidualStressField:
    """Runs a technique to recover the residual stresses which existed in the 
           region of material which was removed by the last committed tool pass
           plan.

       This technique assumes that the residual stress field which existed in
           the region of material was constant before it was removed.
        
       Args:
           commitment_phase_md: The commitment phase of interest. This commitment
                                    phase should contain a committed tool pass
                                    plan composed of exactly a single tool pass.
           path:                Absolute path to directory. Any MDBs which need
                                    to be created in order to recover the residual
                                    stresses are placed this directory.
                                This directory should not be reused for multiple
                                    calls of this method.

       Returns:
           The residual stress field in the material that was removed by the
               committed tool pass plan.
    
       Raises:
           None.
    """

    if commitment_phase_md.committed_tool_pass_plan is not None:
        raise RuntimeError("No tool pass plan was committed to, so no stresses \
                            can be recovered yet.")

    if len(commitment_phase_md.committed_tool_pass_plan) != 1:
        raise RuntimeError("If more than one tool pass happened, then this \
                            the method that this function implements cannot \
                            hope to recover the stresses which existed in \
                            the regions of material removal.")
    
    # Step 1:
    # Create a new MDB for running the necessary simulations for this technique.
    STRESS_RECOVERY_MDB_NAME = "recover_residual_stresses.cae"

    new_mdb_path = os.path.join(path, STRESS_RECOVERY_MDB_NAME)
    if pathlib.Path(new_mdb_path).exists():
        raise RuntimeError("The path for the stress recovery MDB is already occupied!")

    # Use the metadata about the commitment phase to retrieve the name
    #     of the ODB which resulted from the committed tool pass plan being
    #     simulated.
    committed_tpp_mdb_md = commitment_phase_md.per_mdb_metadata[-1]
    odb_name = naming.last_odb_file_name(committed_tpp_mdb_md)
    path_committed_tpp = commitment_phase_md.committed_tool_pass_plan_path
    odb_path = os.path.join(path_committed_tpp, odb_name) 

    mdb, mdb_metadata, names = create_mdb_from_odb(path, STRESS_RECOVERY_MDB_NAME, odb_path)

    # Step 2:
    # Apply a traction vector on the trench of the well in a test direction.

    # 
     

    # Step 3:
    # Fetch the pre-deformed state of the workpiece. This will be used to find
    #     the correct tractions to apply.

    # Step 4:
    # Recover the linear relationship between the applied traction vector and
    #     the displacements / positions of the nodes of interest.

    # Step 5:
    # Use the linear relationship to minimize a measure of geometric dissimilarity.

    # Step 6:
    # Map the "good traction vector" to the components of the stress tensor
    #     assuming a constant state of stress.



def recover_linear_relationships(points: list[geom.Point3D], face: Any,
                                 mdb: Any) -> tuple[Any, ...]:
    """Applies test tractions to some face of a workpiece to recover the linear 
           relationships which exist between applied traction and the displacement 
           of some points.

       Implementation Details:
       This function runs all the simulations necessary to recover the linear
           relationships.
    
       Args:
           points: The points at which to recover the linear relationships.
           face:   Abaqus Face object. The face to apply the tractions to.
                       Must be planar. 
           mdb:    Abaqus MDB object. The MDB containing the part which has the
                       face. This MDB should contain only a basic part geometry
                       (i.e. it should contain a single model, a single part,
                       no non-default steps, etc.)
    
       Returns:
           Tuple of 3 x 3 numpy arrays which describe the linear relationships
               between the applied tractions and the displacements of the
               points. 
           For each matrix M, the matrix M relates the quantities: M * t = d. 
               In this equation, t is a traction vector and d is a displacment vector.
           The number of matrices equals the number of points passed to this
               function.
           Note that these linear relationships only hold when tractions are
               applied to the face passed to this function.
    
       Raises:
           None.
    """

    if not shim.check_basic_geom(True, mdb):
        raise RuntimeError("The MDB should contain a basic geometry.")

     

    





    


    






