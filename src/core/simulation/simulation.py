import os

import core.abaqus.abaqus_shim as shim
import core.tool_pass.tool_pass as tp
import core.boundary_conditions.boundary_conditions as bc
import core.metadata.metadata as md
import core.metadata.naming as naming

# DEBUG
import time
from util.debug import *        


# Simulate consecutive tool passes and save off the results.
#
# Notes:
#    This function has two effects on the file system:
#       Creates a subdirectory in the CWD and runs the simulations. The simulation 
#          artifacts are saved in this new subdirectory. 
#       Saves the .cae file which represents the work for the last tool pass in
#          the subdirectory. Note that this .cae file does not contain the result
#          of the very last tool path geometry. 
#
# Arguments:
#    tool_pass_plan - ToolPassPlan object.
#    name           - String.
#                     Name of the MDB which results from all the tool path simulations.
#                     Also, the name of the directory which is created.
#    record         - CommittedToolPassPlanMetadata object.
#                     The record associated with ithe search for a good next tool pass
#                        to do.
# 
# Optional Arguments:
#    stress_subroutine - String.
#                        Path to stress subroutine. If this argument is not specified,
#                           then the initial stress state of the part will be sourced
#                           from the last simulation in the previous commitment
#                           phase. If this argument is specified, it is used
#                           instead.
#
# Returns:
#    None.
def sim_tool_pass_plan(tool_pass_plan, name, commit_metadata, stress_subroutine=None):
# type: (tp.ToolPassPlan, str, md.CommittedToolPassPlanMetadata, Optional[str]) -> None

    # Create the directory wherein all the simulation artifacts for this sequence
    #    of tool passes will live.
    dir_path = os.path.dirname(commit_metadata.path_initial_mdb)
    new_dir_path = os.path.join(dir_path, name)
    new_mdb_path = os.path.join(new_dir_path, name)
    os.mkdir(new_dir_path)

    # And set the CWD to this directory. Now all simulation artifacts will be placed
    #    in this directory.
    cwd = os.getcwd()
    os.chdir(new_dir_path)

    mdb = shim.use_mdb(commit_metadata.path_initial_mdb)

    for tool_pass in tool_pass_plan.plan:
        sim_single_tool_pass(tool_pass, commit_metadata, mdb, stress_subroutine)    
        tool_pass_plan.pop()

    shim.save_mdb_as(new_mdb_path, mdb)
    shim.close_mdb(mdb)

    commit_metadata.simulated_tool_pass_plans.append((name, tool_pass_plan))

    # Don't permanently change the CWD.
    os.chdir(cwd)



def sim_single_tool_pass(tool_pass, commit_metadata, mdb, stress_subroutine=None):
# type: (tp.ToolPass, md.CommittedToolPassPlanMetadata, Any, Optional[str]) -> None

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
            sim_first_tool_pass(tool_pass, commit_metadata, mdb, stress_subroutine)
        else:
            sim_first_tool_pass(tool_pass, commit_metadata, mdb)  
    elif shim.check_multiple_steps(False, last_model_name, mdb):
        sim_nth_tool_pass(tool_pass, commit_metadata, mdb)
    else:
        raise RuntimeError("Can't figure out how to do next tool pass...")



def sim_first_tool_pass(tool_pass, commit_metadata, mdb, stress_subroutine=None):
# type: (tp.ToolPass, md.CommittedToolPassPlanMetadata, Any, Optional[str]) -> None

    mdb_metadata = commit_metadata.per_mdb_metadata[-1] 

    names = naming.new_model_names(mdb_metadata, True)

    # Do material and section creation.
    shim.create_material(commit_metadata.init_part.material, names["new_model_name"], mdb)
    shim.create_section(names["new_model_name"], mdb)
    shim.assign_section_to_whole_part(names["pre_tool_pass_part_name"], names["new_model_name"], mdb)

    do_boilerplate_sim_ops(tool_pass, names, commit_metadata, mdb)

    if stress_subroutine is not None:
        # Since this modifies the input file directly, this should be the last
        #    thing that happens before the job is submitted and runs.
        shim.inp_add_stress_subroutine(names["new_model_name"], mdb)
    else:
        # Use the stress state from the last toolpath in the previous commit.
        path_sim_file = commit_metadata.path_last_commit_sim_file

        # Since this modifies the input file directly, this should be the last
        #    thing that happens before the job is submitted and runs.
        shim.inp_map_stress(path_sim_file, names["new_model_name"], mdb)

    job = shim.create_job(names["new_model_name"], mdb_metadata, mdb) 

    if stress_subroutine is not None:
        # Associate the user subroutine with the job.
        shim.add_user_subroutine(job, stress_subroutine)

    shim.run_job(job)



def sim_nth_tool_pass(tool_pass, commit_metadata, mdb):
# type: (tp.ToolPass, md.CommittedToolPassPlanMetadata, Any) -> Any

    mdb_metadata = commit_metadata.per_mdb_metadata[-1] 

    names = naming.new_model_names(mdb_metadata, False)
    last_odb_file_name = naming.last_odb_file_name(mdb_metadata)

    # The .sim file contains the stress information which resulted from the last
    #    simulation.
    last_sim_file_name = naming.last_sim_file_name(mdb_metadata)

    # Create the new model with the deformed part in it from the ODB. 
    shim.create_model_and_part_from_odb(names["pre_tool_pass_part_name"], names["new_model_name"], last_odb_file_name, mdb_metadata, mdb)

    # Do material and section creation.
    shim.create_material(commit_metadata.init_part.material, names["new_model_name"], mdb)
    shim.create_section(names["new_model_name"], mdb)

    orphan_mesh_to_geometry(names["pre_tool_pass_part_name"], names["new_model_name"], mdb)

    # Do the section assignment only after the orphan mesh has been mapped to a
    #    geometry. This cannot be done before mapping to a geometry. 
    shim.assign_section_to_whole_part(names["pre_tool_pass_part_name"], names["new_model_name"], mdb)

    do_boilerplate_sim_ops(tool_pass, names, commit_metadata, mdb)

    # Map the stress profile which existed at the end of the last simulation
    #    onto the new geometry. 
    # Since this modifies the input file directly, this should be the last
    #    thing that happens before the job is submitted and runs.
    shim.inp_map_stress(last_sim_file_name, names["new_model_name"], mdb)

    # Create the job and submit it.
    job = shim.create_job(names["new_model_name"], mdb_metadata, mdb) 
    shim.run_job(job)



def do_boilerplate_sim_ops(tool_pass, names, commit_metadata, mdb):
# type: (tp.ToolPass, dict[str, str], md.CommittedToolPassPlanMetadata, Any) -> None

    # Instance the initial geometry part.
    initial_geom_part = shim.get_part(names["pre_tool_pass_part_name"], names["new_model_name"], mdb)
    initial_geom_instance = shim.instance_part_into_assembly(names["pre_tool_pass_part_name"], initial_geom_part, False, names["new_model_name"], mdb)  
    
    # Build the next tool pass path as a part.
    tool_pass_part = shim.build_tool_pass_part(names["tool_pass_part_name"], tool_pass, names["new_model_name"], commit_metadata.per_mdb_metadata[-1], mdb)

    # Instance the tool pass part.
    tool_pass_part_instance = shim.instance_part_into_assembly(names["tool_pass_part_name"], tool_pass_part, False, names["new_model_name"], mdb)

    # Build the tool pass bounding box as a part.
    bounding_box_part = shim.build_bounding_box_part(names["bounding_box_name"], tool_pass, names["new_model_name"], commit_metadata.per_mdb_metadata[-1], mdb)

    # Instance the bounding box part into the assembly.
    bounding_box_instance = shim.instance_part_into_assembly(names["bounding_box_name"], bounding_box_part, False, names["new_model_name"], mdb)

    # Cut off the portion of the bounding box that lives outside the part.
    cut_instances = (initial_geom_instance, )
    bounding_box_excess_part = shim.cut_instances_in_assembly(names["excess_bounding_box_name"], bounding_box_instance, names["new_model_name"], commit_metadata.per_mdb_metadata[-1], mdb)

    # Instance the excess into the assembly.


    # Create the post tool pass geometry as a part.
    cut_instances = (tool_pass_part_instance, )
    post_tool_pass_part = shim.cut_instances_in_assembly(names["post_tool_pass_part_name"], initial_geom_instance, cut_instances, names["new_model_name"], commit_metadata.per_mdb_metadata[-1], mdb)

    # Assign a section to the part.
    shim.assign_only_section_to_part(post_tool_pass_part, names["new_model_name"], mdb)

    # Simplify the part geometry using some heuristics.
    shim.add_virtual_topology(post_tool_pass_part)

    # Instance the post cut geometry.
    # The instance is independent so that meshing can be done in the Assembly
    #    module.
    post_tool_pass_instance = shim.instance_part_into_assembly(names["post_tool_pass_part_name"], post_tool_pass_part, False, names["new_model_name"], mdb)

    # Apply the boundary conditions. 
    bc.apply_BCs(commit_metadata.BCs, shim.STANDARD_INITIAL_STEP_NAME, post_tool_pass_instance, names["new_model_name"], mdb)

    # DEBUG
    mdb.saveAs("test.cae")

    # Then mesh the part in the assembly module.
    shim.naive_mesh(post_tool_pass_instance, 20, names["new_model_name"], mdb)

    # Add an equilibrium step following the last step on commit_metadata.
    last_step_name = commit_metadata.per_mdb_metadata[-1].models_metadata[names["new_model_name"]].step_names[-1]
    shim.create_equilibrium_step(names["equil_step_name"], last_step_name, names["new_model_name"], commit_metadata.per_mdb_metadata[-1], mdb)



# This function adds geometric features to an orphan mesh. The result of a
#    simulation is an orphan mesh (i.e. just a bunch of vertices and elements
#    connecting the vertices together). An orphan mesh by itself is not very
#    useful. For example, to do boolean operations between parts, the parts
#    need to have geometries associated with them. This function adds geoemtric
#    features to a part, giving it a geometry.
#
# Notes:
#    None.
#
# Arguments:
#    part_name  - String.
#    model_name - String.
#    mdb        - Abaqus MDB object.
#
# Returns:
#    None.
def orphan_mesh_to_geometry(part_name, model_name, mdb):
# type: (str, str, Any) -> None 

    part = shim.get_part(part_name, model_name, mdb) 

    unique_elem_faces = shim.get_unique_element_faces(part)

    # For each unique face in the mesh, we know the face is on the surface if it is
    #    associated with exactly one element.
    # We need to construct a region for each face on the surface of the part. 
    # Each region is then used to build a geometric face feature associated
    #    with the part.

    # DEBUG
    start = time.clock()

    for elem_face in unique_elem_faces:
        if len(shim.get_mesh_face_elements(elem_face)) == 1:

            # DEBUG
            t1 = time.clock()

            face_reg = shim.build_region_with_elem_face(elem_face, part)

            # DEBUG
            t2 = time.clock()

            # DEBUG
            t3 = time.clock()

            shim.add_face_from_region(face_reg, part)

            # DEBUG
            t4 = time.clock()
            dp("For element face with index " + str(elem_face.label) + " the wall clock time for building the face region was " + str(t2 - t1) + " and the wall clock time for adding the face from the region was " + str(t4 - t3))

    # DEBUG
    end = time.clock()
    dp("The total wall clock time for building the geometric faces from the mesh faces was " + str(end - start))

    # DEBUG
    start = time.clock()

    # Build the solid feature from the face features.
    shim.add_solid_from_faces(part)

    # DEBUG
    end = time.clock()
    dp("The total wall clock time for building the solid from the faces was " + str(end - start))

    # Remove any dependency on the orphan mesh.
    # Not doing this causes the orphan mesh to still appear in the Assembly
    #    module, which can make appearence confusing. 
    shim.suppress_feature(shim.STANDARD_ORPHAN_MESH_FEATURE_NAME, part)