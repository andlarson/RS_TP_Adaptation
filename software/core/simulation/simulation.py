import core.abaqus.abaqus_shim as shim
import core.tool_pass.tool_pass as tp
import core.abaqus.abaqus_metadata as abq_md
import core.tool_pass.tool_pass_record_keeping as tprk
import core.boundary_conditions.boundary_conditions as bc

import os

# DEBUG
from util.debug import *



def sim_consecutive_tool_passes(tool_pass_plan, save_path, record):
# type: (tp.ToolPassPlan, str, tprk.CommittedToolPassRecord) -> None

    # The results of the simulations go into the same directory as the MDB. 
    os.chdir(record.working_dir)

    mdb = shim.use_mdb(record.path_to_mdb)

    tool_pass = tool_pass_plan.pop() 
    while tool_pass != None:
        sim_single_tool_pass(tool_pass, tool_pass_plan.done_so_far(), record, mdb)    

        # DEBUG
        dp("Done simulating tool pass number " + str(tool_pass_plan.done_so_far()))

        tool_pass = tool_pass_plan.pop()
    
    shim.save_mdb_as(save_path, mdb)
    shim.close_mdb(mdb)



def sim_single_tool_pass(tool_pass, tool_pass_cnt, record, mdb):
# type: (tp.ToolPass, int, tprk.CommittedToolPassRecord, Any) -> None

    """
    The way that the next tool pass is simulated depends on the content of
       the MDB.
    If the MDB contains a single model which represents an initial
       geometry, then it's necessary to simulate the first tool pass. 
    If the MDB contains n (where n > 1) models and the last model to be 
       added contains an orphan mesh, then we must be simulating the n-th 
       tool pass.
    """
    last_model_name = record.abq_metadata.get_last_model_name()
    if shim.check_init_geom(False, mdb):
        
        # If we have an initial geometry, there better be an associated
        #    stress profile.
        path_to_stress_subroutine = record.init_part.path_to_stress_subroutine
        assert(path_to_stress_subroutine != None)

        # DEBUG
        dp("Simulating a first tool pass!")

        sim_first_tool_pass(tool_pass, tool_pass_cnt, record.BCs, record.abq_metadata, path_to_stress_subroutine, mdb)  

    elif shim.check_multiple_steps(False, last_model_name, mdb):

        # DEBUG
        dp("Simulating an nth tool pass!")

        # If there are multiple steps in the last model, it is assumed that
        #    the model was run. Therefore, it's necessary to use the results
        #    of that model to chain the simulations together. 

        last_odb_name = record.abq_metadata.get_last_job_name(last_model_name) + ".odb"
        last_part_name = record.abq_metadata.get_last_part_name(last_model_name).upper()
        sim_nth_tool_pass(tool_pass, tool_pass_cnt, last_part_name, last_odb_name, record.BCs, record.abq_metadata, mdb)

    else:
        raise RuntimeError("Can't figure out how to do next tool pass...")



# This works for the first tool pass in an MDB.
# Because it's the first tool pass, assumptions are made about how things are
#    named.
def sim_first_tool_pass(tool_pass, tool_pass_cnt, BCs, abq_metadata, path_to_stress_subroutine, mdb):
# type: (tp.ToolPass, int, List[bc.BC], abq_md.ABQMetadata, str, Any) -> Any

    model_name = shim.STANDARD_MODEL_NAME
    orig_part_name = shim.STANDARD_INIT_GEOM_PART_NAME
    tool_pass_part_name = shim.STANDARD_TOOL_PASS_PART_PREFIX + str(tool_pass_cnt)
    post_tool_pass_part_name = shim.STANDARD_POST_TOOL_PASS_PART_PREFIX + str(tool_pass_cnt)
    step_cnt = shim.get_step_cnt(model_name, mdb)
    equil_step_name = shim.STANDARD_EQUIL_STEP_PREFIX + str(step_cnt + 1)

    do_boilerplate_sim_ops(orig_part_name, tool_pass_part_name, post_tool_pass_part_name, equil_step_name, model_name, tool_pass, BCs, abq_metadata, mdb) 

    """
    Now, just before creating and submitting the job, modify the input 
       file so that the stress profile is included!
    There is nuance here:
       1) Each model has an input file which mirrors the functionality in
             the model. Further modifications to the model after the input
             file has been (essentially) manually modified can cause things
             to get out of sync and be messed up. This should be the last 
             thing done before the job is created and runs.
       2) We are really superimposing the user subroutine defined stress
             profile on the part which has undergone material removal via
             a boolean operation.
    """
    shim.inp_add_stress_subroutine(model_name, mdb)

    job_name = post_tool_pass_part_name
    job = shim.create_job(job_name, model_name, abq_metadata, mdb) 

    # Associate the user subroutine with the job.
    shim.add_user_subroutine(job, path_to_stress_subroutine)

    shim.run_job(job)

    return mdb 




"""
If the source is a previous simulation, import the results of the simulation 
   into a new model in the working MDB, add the next cut, and go from 
   there.
The importing process isn't trivial. The output of the previous
   simulation is an ODB which contains a representation of a deformed
   mesh. This deformed mesh can be imported into Abaqus and viewed.
   This mesh is an "Orphan Mesh" and has no associated underlying
   part geometry. In order or simulate the next cut via the boolean
   cut procedure, it's necessary to generate a part geometry from the
   underlying mesh.  
It also requires mapping the stress state which existed at the end of
   the previous simulation in as the initial stress state.
"""
def sim_nth_tool_pass(tool_pass, tool_pass_cnt, last_part_name, last_odb_name, BCs, abq_metadata, mdb):
# type: (tp.ToolPass, int, str, str, List[bc.BC], abq_md.ABQMetadata, Any) -> Any

    """
    1) Create a new model, import the orphan mesh, map the orphan mesh to a
          geometry.
    """
    new_model_name = shim.STANDARD_MODEL_NAME_PREFIX + str(tool_pass_cnt) 
    new_model = shim.create_model_from_odb(last_odb_name, new_model_name, abq_metadata, mdb)

    orphan_mesh_to_geometry(last_part_name, new_model_name, mdb)


    """
    2) Do all of the generic simulation preparation.
    """
    tool_pass_part_name = shim.STANDARD_TOOL_PASS_PART_PREFIX + str(tool_pass_cnt)
    post_tool_pass_part_name = shim.STANDARD_POST_TOOL_PASS_PART_PREFIX + str(tool_pass_cnt)
    step_cnt = shim.get_step_cnt(new_model_name, mdb)
    equil_step_name = shim.STANDARD_EQUIL_STEP_PREFIX + str(step_cnt + 1)

    do_boilerplate_sim_ops(last_part_name, tool_pass_part_name, post_tool_pass_part_name, equil_step_name, new_model_name, tool_pass, BCs, abq_metadata, mdb)
       

    """
    3) Map the stress profile which existed at the end of the last simulation
          onto the new geometry. 
       Since this modifies the input file directly, this should be the last
          thing that happens before the job is submitted and runs.
    """
    sim_file_path = last_part_name + ".sim"
    shim.inp_map_stress(sim_file_path, new_model_name, mdb)


    """
    4) Create the job and submit it.
    """
    job_name = post_tool_pass_part_name
    job = shim.create_job(job_name, new_model_name, abq_metadata, mdb) 
    shim.run_job(job)

    return mdb



# Executes a sequence of operations which commonly occur in preparation for running
#    a simulation.
def do_boilerplate_sim_ops(orig_part_name, tool_pass_part_name, post_tool_pass_part_name, equil_step_name, model_name, tool_pass, BCs, abq_metadata, mdb):
# type: (str, str, str, str, Any, List[bc.BC], abq_md.ABQMetadata)

    # Instance the initial geometry part.
    initial_geom_part = shim.get_part(orig_part_name, model_name, mdb)
    initial_geom_instance = shim.instance_part_into_assembly(orig_part_name, initial_geom_part, True, model_name, mdb)  
    
    # Build the next tool pass path as a part.
    tool_pass_geom = tool_pass.geom
    tool_pass_part = shim.build_part(tool_pass_part_name, tool_pass_geom, model_name, abq_metadata, mdb)

    # Instance the tool pass part.
    tool_pass_part_instance = shim.instance_part_into_assembly(tool_pass_part_name, tool_pass_part, True, model_name, mdb)

    # Create the post tool pass geometry as a part.
    cut_instances = (tool_pass_part_instance, )
    post_tool_pass_part = shim.cut_instances_in_assembly(post_tool_pass_part_name, initial_geom_instance, cut_instances, model_name, abq_metadata, mdb)

    # Assign a section to the part.
    shim.assign_only_section_to_part(post_tool_pass_part, model_name, mdb)

    # Simplify the part geometry using some heuristics.
    shim.add_virtual_topology(post_tool_pass_part)

    # Instance the post cut geometry.
    # The instance is independent so that meshing can be done in the Assembly
    #    module.
    post_tool_pass_instance = shim.instance_part_into_assembly(post_tool_pass_part_name, post_tool_pass_part, False, model_name, mdb)

    # Then mesh the part in the assembly module.
    shim.naive_mesh(post_tool_pass_instance, 30, model_name, mdb)

    # Add an equilibrium step following the last step on record.
    last_step = abq_metadata.get_last_step(model_name)
    shim.create_equilibrium_step(equil_step_name, last_step, model_name, abq_metadata, mdb)

    # Apply the boundary conditions. 
    bc.apply_surface_BCs(BCs, equil_step_name, model_name, mdb)



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

    """
    For each unique face in the mesh, we know the face is on the surface if it is
       associated with exactly one element.
    We need to construct a region for each face on the surface of the part. 
    Each region is then used to build a geometric face feature associated
       with the part.
    """
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
