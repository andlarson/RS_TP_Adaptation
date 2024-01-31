"""
Provides top-level functionality to do simulations of interest in Abaqus.
"""

import os
import time
from typing import Optional, Any
import pathlib
import shutil
import sys

from abaqus import *
from abaqusConstants import * 

import src.core.abaqus.abaqus_shim as shim
import src.core.tool_pass.tool_pass as tp
import src.core.boundary_conditions.boundary_conditions as bc
import src.core.metadata.metadata as md
import src.core.metadata.naming as naming

from src.util.debug import *        


def sim_tool_pass_plan(tool_pass_plan: tp.ToolPassPlan, name: str, target_dir: str, 
                       commit_metadata: md.CommittedToolPassPlanMetadata, 
                       stress_subroutine: Optional[str] = None) -> None:
    """Simulates consecutive tool passes and saves off the results.
    
       Args:
           tool_pass_plan:    The tool pass plan to simulate. 
           name:              Desired name of the MDB which results from all the tool 
                                  path simulations. Also, the desired name of the
                                  subdirectory created by this function.
           target_dir:        Absolute path to directory. In this directory, a 
                                  subdirectory is created and named name. The
                                  simulation artifacts (.cae, .odb, .msg, etc.)
                                  are generated in this subdirectory.
                              If a directory or file with the desired name already
                                  exists in this directory, it is deleted.
           commit_metadata:   The record associated with the search for a good next 
                                  tool pass to do.
           stress_subroutine: Path to stress subroutine. If this argument is not 
                                  specified, then the initial stress state of the 
                                  part will be sourced from the last simulation in 
                                  the previous commitment phase. If this argument 
                                  is specified, it is used instead.
       
       Returns:
           None.

       Raises:
           RuntimeError: An entity which is neither file nor directory exists
                             at the file system location where the subdirectory
                             needs to be created.
    """

    # Create the directory wherein all the simulation artifacts for this sequence
    #    of tool passes will live.
    new_dir_path = os.path.join(target_dir, name)
    new_mdb_path = os.path.join(new_dir_path, name)

    # If a directory or file of the same name already exists, delete it.
    subdir = pathlib.Path(new_dir_path)     
    if subdir.exists():
        if subdir.is_dir():
            dp("Removing directory at path: " + new_dir_path)
            shutil.rmtree(new_dir_path) 
        elif subdir.is_file():
            dp("Removing file at path: " + new_dir_path)
            os.remove(new_dir_path)
        else:
            raise RuntimeError("Something which is neither file system nor file\
                                exists at the desired location of the subdirectory.")

    os.mkdir(new_dir_path)

    # And set the CWD to this directory. Now all simulation artifacts will be placed
    #    in this directory.
    cwd = os.getcwd()
    os.chdir(new_dir_path)

    mdb = shim.use_mdb(commit_metadata.path_initial_mdb)

    for tool_pass in tool_pass_plan.plan:
        _sim_single_tool_pass(tool_pass, commit_metadata, mdb, stress_subroutine)    
        tool_pass_plan.pop()

    shim.save_mdb_as(new_mdb_path, mdb)
    shim.close_mdb(mdb)

    commit_metadata.simulated_tool_pass_plans.append((name, tool_pass_plan))

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

    _do_boilerplate_sim_ops(tool_pass, names, commit_metadata, mdb)

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

    orphan_mesh_to_geometry(names.pre_tool_pass_part_name, names.new_model_name, mdb)

    # Do the section assignment only after the orphan mesh has been mapped to a
    #    geometry. This cannot be done before mapping to a geometry. 
    shim.assign_section_to_whole_part(names.pre_tool_pass_part_name, names.new_model_name, mdb)

    _do_boilerplate_sim_ops(tool_pass, names, commit_metadata, mdb)

    # Map the stress profile which existed at the end of the last simulation
    #    onto the new geometry. 
    # Since this modifies the input file directly, this should be the last
    #    thing that happens before the job is submitted and runs.
    shim.inp_map_stress(last_sim_file_name, names.new_model_name, mdb)

    # Create the job and submit it.
    job = shim.create_job(names.new_model_name, mdb_metadata, mdb) 
    shim.run_job(job)



def _do_boilerplate_sim_ops(tool_pass: tp.ToolPass, names: naming.ModelNames, 
                            commit_metadata: md.CommittedToolPassPlanMetadata, 
                            mdb: Any) -> None:
    """Does some operations such as instancing, section assignment, etc. which
           are boilerplate (e.g. they need to be done just about any tool pass
           simulation).

       Args:
           tool_pass:         The tool pass to simulate.
           names:             The names associated with the MDB.
           commit_metadata:   The record associated with the search for a good next 
                                  tool pass to do.
           mdb:               Abaqus MDB object. The MDB in which the tool pass
                                  will be simulated. 

       Returns:
           None.

       Raises:
           None.
    """

    # Instance the initial geometry part.
    initial_geom_part = mdb.models[names.new_model_name].parts[names.pre_tool_pass_part_name]
    initial_geom_instance = shim.instance_part_into_assembly(names.pre_tool_pass_part_name, initial_geom_part, False, names.new_model_name, mdb)  

    """DEPRECATED. No longer using bounding box technique to do mesh refinement.

    # Build the tool pass bounding box as a part.
    bounding_box_part = shim.build_bounding_box_part(names["bounding_box_name"], tool_pass, names["new_model_name"], commit_metadata.per_mdb_metadata[-1], mdb)
    bounding_box_instance = shim.instance_part_into_assembly(names["bounding_box_name"], bounding_box_part, False, names["new_model_name"], mdb)

    # Cut off the portion of the bounding box that lives outside the part.
    cutting_instances = (initial_geom_instance, )
    bounding_box_excess_part = shim.cut_instances_in_assembly(names["excess_bounding_box_name"], bounding_box_instance, cutting_instances, names["new_model_name"], commit_metadata.per_mdb_metadata[-1], mdb)
    bounding_box_excess_part_instance = shim.instance_part_into_assembly(names["excess_bounding_box_name"], bounding_box_excess_part, False, names["new_model_name"], mdb)

    # Resume (i.e. unsuppress) the instances that were suppressed during the cut operation.
    shim.resume_all_assembly_features(names["new_model_name"], mdb)

    # Merge the initial geometry and the bounding box. 
    merging_instances = (initial_geom_instance, bounding_box_instance)
    merged_part = shim.merge_instances_in_assembly(names["merged_bbox_geom_name"], merging_instances, True, names["new_model_name"], commit_metadata.per_mdb_metadata[-1], mdb)
    merged_instance = shim.instance_part_into_assembly(names["merged_bbox_geom_name"], merged_part, False, names["new_model_name"], mdb)

    # Cut off the portion of the bounding box that is outside the initial geometry. 
    cutting_instances = (bounding_box_excess_part_instance, )
    initial_geom_with_bbox_part = shim.cut_instances_in_assembly(names["init_geom_with_bbox_name"], merged_instance, cutting_instances, names["new_model_name"], commit_metadata.per_mdb_metadata[-1], mdb)
    initial_geom_with_bbox_instance = shim.instance_part_into_assembly(names["init_geom_with_bbox_name"], initial_geom_with_bbox_part, False, names["new_model_name"], mdb)
    """
    
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



def orphan_mesh_to_geometry(part_name: str, model_name: str, mdb: Any) -> None:
    """Converts an orphan mesh to a solid geometry by converting all the faces
           of the orphan mesh to regions/surfaces, and then melding all the
           regions/surfaces together to form a solid geometr
           
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

    # DEBUG
    debug = False

    # DEBUG
    if debug:
        start = time.clock()

    cnt = 0
    for elem_face in unique_elem_faces:
        if len(shim.get_mesh_face_elements(elem_face)) == 1:

            # DEBUG
            if debug:
                t1 = time.clock()

            face_reg = shim.build_region_with_elem_face(elem_face, part)

            # DEBUG
            if debug:
                t2 = time.clock()

            shim.add_face_from_region(face_reg, part)

            # DEBUG
            if debug:
                t4 = time.clock()
                dp("For element face with index " + str(elem_face.label) + " the wall clock time for building the face region was " + str(t2 - t1) + " and the wall clock time for adding the face from the region was " + str(t4 - t2))

            # Clear the cache every once in a while for a potential speed up. 
            cnt += 1
            if cnt % 50 == 0:
                part.clearGeometryCache()

    # DEBUG
    if debug:
        end = time.clock()
        dp("The total wall clock time for building the geometric faces from the mesh faces was " + str(end - start))

    """DEPRECATED. The convert_shell_to_solid() function is significantly less
           failure prone than add_solid_from_faces().

    # DEBUG
    start = time.clock()

    # Build the solid feature from the face features.
    try:
        shim.add_solid_from_faces(part)
    except AbaqusException as e:
        dp("Shell to Solid conversion failure!")
        dp("Saving off the MDB as failed_shell_to_solid.cae")
        mdb.saveAs("failed_shell_to_solid.cae")
        dp("The arguments associated with the exception are " + str(e.args))
        raise

    # DEBUG
    end = time.clock()
    dp("The total wall clock time for building the solid from the faces was " + str(end - start))
    """

    shim.convert_shell_to_solid(part)
    
    # Remove any dependency on the orphan mesh.
    # Not doing this causes the orphan mesh to still appear in the Assembly
    #    module, which can make appearence confusing. 
    shim.suppress_feature(shim.STANDARD_ORPHAN_MESH_FEATURE_NAME, part)



