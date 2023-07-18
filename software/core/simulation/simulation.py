import core.abaqus.abaqus_shim as shim
import core.tool_pass.tool_pass as tp
import core.abaqus.abaqus_metadata as abq_md



# This works for the first tool pass in an MDB.
# Because it's the first tool pass, assumptions are made about the names of
#    entities.
def sim_first_tool_pass(tool_pass, tool_pass_cnt, abq_metadata, path_to_stress_subroutine, mdb):
# type: (tp.ToolPass, int, abq_md.ABQMetadata, str, Any) -> Any

    model_name = shim.STANDARD_MODEL_NAME

    # Instance the initial geometry part.
    initial_geom_part = shim.get_part(shim.STANDARD_INIT_GEOM_PART_NAME, mdb)
    initial_geom_instance = shim.instance_part_into_assembly(shim.STANDARD_INIT_GEOM_PART_NAME, initial_geom_part, True, model_name, mdb)  
    
    # Build the next tool pass path as a part.
    tool_pass_geom = tool_pass.geom
    tool_pass_name = shim.STANDARD_TOOL_PASS_PART_PREFIX + str(tool_pass_cnt)
    tool_pass_part = shim.build_part(tool_pass_name, tool_pass_geom, model_name, abq_metadata, mdb)

    # Instance the tool pass part.
    tool_pass_part_instance = shim.instance_part_into_assembly(tool_pass_name, tool_pass_part, True, model_name, mdb)

    # Create the post tool pass geometry as a part.
    post_tool_pass_name = shim.STANDARD_POST_TOOL_PASS_PART_PREFIX + str(tool_pass_cnt)
    cut_instances = (tool_pass_part_instance, )
    post_tool_pass_part = shim.cut_instances_in_assembly(post_tool_pass_name, initial_geom_instance, cut_instances, model_name, mdb)

    # Assign a section to the post tool pass geometry part.
    shim.assign_standard_sec_to_simple_part(post_tool_pass_part)

    # Instance the post cut geometry.
    # The instance is independent so that meshing can be done in the Assembly
    #    module.
    post_tool_pass_instance = shim.instance_part_into_assembly(post_tool_pass_name, post_tool_pass_part, False, model_name, mdb)

    # Then mesh the part in the assembly module.
    shim.naive_mesh(post_tool_pass_instance, 1, model_name, mdb)

    # Add an equilibrium step following the last step on record.
    last_step = abq_metadata.get_last_step(model_name)
    step_cnt = shim.get_step_cnt(model_name, mdb)
    equil_step_name = shim.STANDARD_EQUIL_STEP_PREFIX + str(step_cnt + 1)
    shim.create_equilibrium_step(equil_step_name, last_step, model_name, abq_metadata, mdb)

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
    shim.modify_inp("*Initial Conditions", ("Type=Stress", "User"), "", model_name, mdb)

    job_name = tool_pass_name
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
   underlying mesh. Doing so requires mapping mesh element faces ->
   part surfaces, and then (potntially) using the virtual topology
   toolset to simplify the part geometry (removing all the unnecessary
   surface partitions, but also potentially degrading the part
   representation). 
It also requires mapping the stress values which existed at the end
   of the previous simulation as the stress values which exist at the
   beginning of the current simulation.
"""
def sim_nth_tool_pass(tool_pass, tool_pass_cnt, last_part_name, last_odb_name, abq_metadata, mdb):
# type: (tp.ToolPass, int, str, str, abq_md.ABQMetadata, Any) -> Any

    # 1) Create a new model, import the orphan mesh, map the orphan mesh to an
    #    underlying geometry.
    
    new_model_name = shim.STANDARD_MODEL_PREFIX + str(tool_pass_cnt) 
    new_model = shim.create_model_from_odb(last_odb_name, new_model_name, abq_metadata, mdb)

    convert_orphan_mesh_to_part()

    # Map the resulting stress profile onto the new geometry.

    # ~~~~ Everything below is repeated behavior from previous function. ~~~~~ 
    # Create the tool pass part, do the boolean intersection, assign a section
    #    to the resulting part.

    # Instance the part into the assembly.

    # Mesh the post tool pass part in the assembly module.

    # Create the equilbirium step.

    # Create the job and submit it.



# The orphan mesh is assumed to be the only thing in the model. 
def convert_orphan_mesh_to_part(part_name, model_name, mdb)

    assert(shim.check_orphan_mesh(part_name, model_name, mdb))

    new_part = 
      
     

