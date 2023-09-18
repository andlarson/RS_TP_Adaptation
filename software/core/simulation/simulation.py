import os

import core.abaqus.abaqus_shim as shim
import core.tool_pass.tool_pass as tp
import core.boundary_conditions.boundary_conditions as bc
import core.metadata.metadata as md
import core.metadata.abaqus_metadata as abq_md

from util.debug import *        



def sim_consecutive_tool_passes(tool_pass_plan, save_path, record):
# type: (tp.ToolPassPlan, str, md.CommittedToolPassMetadata) -> None

    # This dictates where the results of the simulations go. 
    os.chdir(record.working_dir)

    mdb = shim.use_mdb(record.path_to_mdb)

    tool_pass = tool_pass_plan.pop() 
    while tool_pass is not None:

        sim_single_tool_pass(tool_pass, record, mdb)    
        record.simulated_tool_passes.append(tool_pass)

        tool_pass = tool_pass_plan.pop()
    
    shim.save_mdb_as(save_path, mdb)
    shim.close_mdb(mdb)



def sim_single_tool_pass(tool_pass, record, mdb):
# type: (tp.ToolPass, md.CommittedToolPassMetadata, Any) -> None

    # The way that the next tool pass is simulated depends on the content of
    #    the MDB.
    # If the MDB contains a single model which represents an initial
    #    geometry, then it's necessary to simulate the first tool pass. 
    # If the MDB contains n (where n > 1) models and the last model to be 
    #    added contains an orphan mesh, then we must be simulating the n-th 
    #    tool pass.
    last_model_name = record.abaqus_mdb_metadata.model_names[-1]

    if shim.check_init_geom(False, mdb):

        sim_first_tool_pass(tool_pass, record, mdb)  

    elif shim.check_multiple_steps(False, last_model_name, mdb):

        sim_nth_tool_pass(tool_pass, record, mdb)

    else:
        raise RuntimeError("Can't figure out how to do next tool pass...")



# Simulating the very first tool pass in an MDB.
def sim_first_tool_pass(tool_pass, record, mdb):
# type: (tp.ToolPass, md.CommittedToolPassMetadata, Any) -> None

    names = abq_md.SimNames(record)
    do_boilerplate_sim_ops(tool_pass, names, record, mdb)

    # Very important note: This operation effectively imposes the user
    #    subroutine defined stress profile on the part after the material has
    #    been removed! In order for this to make any sense, the user subroutine
    #    defined stress profile must achieve equilbrium!
    # Since this modifies the input file directly, this should be the last
    #    thing that happens before the job is submitted and runs.
    shim.inp_add_stress_subroutine(names.new_model_name, mdb)

    job = shim.create_job(names.post_tool_pass_part_name, names.new_model_name, record, mdb) 

    # Associate the user subroutine with the job.
    shim.add_user_subroutine(job, record.init_part.path_to_stress_subroutine)

    shim.run_job(job)



def sim_nth_tool_pass(tool_pass, record, mdb):
# type: (tp.ToolPass, md.CommittedToolPassMetadata, Any) -> Any

    names = abq_md.SimNames(record) 

    # Create the new model from the output of the last model.
    new_model = shim.create_model_from_odb(names.last_model_odb_file_name, names.new_model_name, record, mdb)

    orphan_mesh_to_geometry(names.pre_tool_pass_part_name, names.new_model_name, mdb)

    do_boilerplate_sim_ops(tool_pass, names, record, mdb)
       
    # The .sim file produced by the last simulation has the same name as the last part
    #    created in the last simulation.
    sim_file_path = names.pre_tool_pass_part_name + ".sim"

    # Map the stress profile which existed at the end of the last simulation
    #    onto the new geometry. 
    # Since this modifies the input file directly, this should be the last
    #    thing that happens before the job is submitted and runs.
    shim.inp_map_stress(sim_file_path, names.new_model_name, mdb)

    # Create the job and submit it.
    job = shim.create_job(names.post_tool_pass_part_name, names.new_model_name, record, mdb) 
    shim.run_job(job)



def do_boilerplate_sim_ops(tool_pass, names, record, mdb):
# type: (tp.ToolPass, SimNames, md.CommittedToolPassMetadata, Any) -> None

    # Instance the initial geometry part.
    initial_geom_part = shim.get_part(names.pre_tool_pass_part_name, names.new_model_name, mdb)
    initial_geom_instance = shim.instance_part_into_assembly(names.pre_tool_pass_part_name, initial_geom_part, True, names.new_model_name, mdb)  
    
    # Build the next tool pass path as a part.
    tool_pass_geom = tool_pass.geom
    tool_pass_part = shim.build_part(names.tool_pass_part_name, tool_pass_geom, names.new_model_name, record, mdb)

    # Instance the tool pass part.
    tool_pass_part_instance = shim.instance_part_into_assembly(names.tool_pass_part_name, tool_pass_part, True, names.new_model_name, mdb)

    # Create the post tool pass geometry as a part.
    cut_instances = (tool_pass_part_instance, )
    post_tool_pass_part = shim.cut_instances_in_assembly(names.post_tool_pass_part_name, initial_geom_instance, cut_instances, names.new_model_name, record, mdb)

    # Assign a section to the part.
    shim.assign_only_section_to_part(post_tool_pass_part, names.new_model_name, mdb)

    # Simplify the part geometry using some heuristics.
    shim.add_virtual_topology(post_tool_pass_part)

    # Instance the post cut geometry.
    # The instance is independent so that meshing can be done in the Assembly
    #    module.
    post_tool_pass_instance = shim.instance_part_into_assembly(names.post_tool_pass_part_name, post_tool_pass_part, False, names.new_model_name, mdb)

    # Apply the boundary conditions. 
    bc.apply_BCs(record.BCs, shim.STANDARD_INITIAL_STEP_NAME, post_tool_pass_instance, names.new_model_name, mdb)

    # Then mesh the part in the assembly module.
    shim.naive_mesh(post_tool_pass_instance, 30, names.new_model_name, mdb)

    # Add an equilibrium step following the last step on record.
    last_step_name = record.abaqus_mdb_metadata.models_metadata[names.new_model_name].step_names[-1]
    shim.create_equilibrium_step(names.equil_step_name, last_step_name, names.new_model_name, record, mdb)



"""
This function adds geometric features to an orphan mesh. The result of a
   simulation is an orphan mesh (i.e. just a bunch of vertices and elements
   connecting the vertices together). An orphan mesh by itself is not very
   useful. For example, to do boolean operations between parts, the parts
   need to have geometries associated with them. This function adds geoemtric
   features to a part, giving it a geometry.
"""
def orphan_mesh_to_geometry(part_name, model_name, mdb):
# type: (str, str, Any) -> None 

    assert(shim.check_orphan_mesh(True, part_name, model_name, mdb))

    part = shim.get_part(part_name, model_name, mdb) 
    unique_elem_faces = shim.get_unique_element_faces(part)

    # For each unique face in the mesh, we know the face is on the surface if it is
    #    associated with exactly one element.
    # We need to construct a region for each face on the surface of the part. 
    # Each region is then used to build a geometric face feature associated
    #    with the part.
    for elem_face in unique_elem_faces:
        if len(shim.get_mesh_face_elements(elem_face)) == 1:
            face_reg = shim.build_region_with_elem_face(elem_face, part)
            face_feature = shim.add_face_from_region(face_reg, part)

    # Build the solid feature from the face features.
    shim.add_solid_from_faces(part)

    # Remove any dependency on the orphan mesh.
    # Not doing this causes the orphan mesh to still appear in the Assembly
    #    module, which can make appearence confusing. 
    shim.suppress_feature(shim.STANDARD_ORPHAN_MESH_FEATURE_NAME, part)


