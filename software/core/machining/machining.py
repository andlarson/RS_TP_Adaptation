import core.part.part as part
import core.abaqus.abaqus_shim as shim
import core.material_properties as mat_props 
import core.simulation.simulation as sim
import core.tool_pass.tool_pass as tp
import core.abaqus.metadata as md
import storage.storage as storage

import os




class MachiningProcess:
   
    def __init__(self, name, init_part, tool_pass_plan):
    # type: (str, part.Part, tp.ToolPassPlan) -> None
    
        self.name = name
        self.part = init_part
        self.tool_pass_plan = tool_pass_plan 
        self.part_history = part.PartHistory() 

        # A single machining process may contain many model databases (i.e. MDBs).
        #    This gives the most flexibility. We may want to chain simulations
        #    in a single MDB or do a single simulation per MDB to make things as
        #    modular. 
        # This maps MDB paths to metadata objects.
        self.mdb_to_metadata = {}

        if isinstance(init_part, part.UserDefinedPart):
            # TODO: We don't handle this case right now...
            #       This should build an MDB and then append it to the list
            #          of simulations. Behavior should be symmetric to the
            #          else clause.
            pass
        else:
            assert(name == None)
            self.mdb_to_metadata[init_part.path_to_mdb] = md.MDBMetadata(init_part.path_to_mdb)

            # In general, all modifications to MDB metadata should happen in the
            #    abaqus_shim.py file at the same time as modifications are made
            #    to the MDB.
            # This is an exceptional case because the user has passed an MDB,
            #    so it's necessary to sync the metadata with the content of the
            #    MDB.
            self.mdb_to_meta_data[init_part.path_to_mdb].add_model(STANDARD_MODEL_NAME)

            self.working_mdb_path = init_part.path_to_mdb



    def sim_next_tool_pass(self, save_location):
    # type: (str) -> None
       
        mdb = shim.use_mdb(self.working_mdb_path)
        metadata = self.mdb_to_metadata[self.working_mdb_path]

        # The way that the next tool pass is simulated depends on the content of
        #    the MDB.
        # If the MDB contains a single model which represents an initial
        #    geometry, then we must be simulating the first tool pass.
        # If the MDB contains n (where n > 1) models and the last model to be 
        #    added contains an orphan mesh, then we must be simulating the n-th 
        #    tool pass.
        last_model_name = metadata.get_last_model_name()
        last_part_name = metadata.get_last_part_name()
        if shim.check_init_geom(mdb):
              
        elif shim.check_orphan_mesh(last_part_name, last_model_name, mdb):

        else:
            raise RuntimeError("Unable to identify how to continue with the \
                                next tool pass simulation!")

        # If the source is a simple part geometry specification, then build the
        #    model up and run a first simulation.

        # If the source is a previous simulation, import the results of the simulation 
        #    into a new model in the working MDB, add the next cut, and go from 
        #    there.
        # The importing process isn't trivial. The output of the previous
        #    simulation is an ODB which contains a representation of a deformed
        #    mesh. This deformed mesh can be imported into Abaqus and viewed.
        #    This mesh is an "Orphan Mesh" and has no associated underlying
        #    part geometry. In order or simulate the next cut via the boolean
        #    cut procedure, it's necessary to generate a part geometry from the
        #    underlying mesh. Doing so requires mapping mesh element faces ->
        #    part surfaces, and then (potntially) using the virtual topology
        #    toolset to simplify the part geometry (removing all the unnecessary
        #    surface partitions, but also potentially degrading the part
        #    representation). 
        # It also requires mapping the stress values which existed at the end
        #    of the previous simulation as the stress values which exist at the
        #    beginning of the current simulation.

        # We're only ever working with a single MDB, so there is just a single
        #    save location.

        # Instance the initial geometry part.
        # We assume the presence of an initial geometry part. 
        initial_geom_part = shim.get_part(shim.STANDARD_INIT_GEOM_PART_NAME, mdb)
        initial_geom_instance = shim.instance_part_into_assembly(shim.STANDARD_INIT_GEOM_PART_NAME, initial_geom_part, True, model_name, mdb)  

        # Build the next tool pass path as a part.
        tool_pass = self.tool_pass_plan.pop()
        tool_pass_geom = tool_pass.geom
        tool_pass_name = shim.STANDARD_TOOL_PASS_PART_PREFIX + str(self.tool_pass_plan.done_so_far())
        tool_pass_part = shim.build_part(tool_pass_name, tool_pass_geom, model_name, mdb)

        # Instance the tool pass part.
        tool_pass_part_instance = shim.instance_part_into_assembly(tool_pass_name, tool_pass_part, True, model_name, mdb)

        # Create the post tool pass geometry as a part.
        post_tool_pass_name = shim.STANDARD_POST_TOOL_PASS_PART_PREFIX + str(self.tool_pass_plan.done_so_far())
        cut_instances = (tool_pass_part_instance, )
        post_tool_pass_part = shim.cut_instances_in_assembly(post_tool_pass_name, initial_geom_instance, cut_instances, model_name, mdb)

        # Assign a section to the post tool pass geometry part.
        shim.assign_standard_sec_to_simple_part(post_tool_pass_part)

        # Instance the post cut geometry.
        # The instance is independent so that meshing can be done in the Assembly
        #    module.
        post_tool_pass_instance = shim.instance_part_into_assembly(post_tool_pass_name, post_tool_pass_part, False, model_name, mdb)

        # Then mesh the part in the assembly module.
        # Picking a seed density randomly....
        shim.naive_mesh(post_tool_pass_instance, 1, model_name, mdb)

        # Add an equilibrium step following the last step on record.
        last_step = metadata.get_last_step(model_name)
        step_cnt = shim.get_step_cnt(model_name, mdb)
        equil_step_name = shim.STANDARD_EQUIL_STEP_PREFIX + str(step_cnt + 1)
        shim.create_equilibrium_step(equil_step_name, last_step, model_name, metadata, mdb)

        # Now, just before creating and submitting the job, modify the input 
        #    file so that the stress profile is included!
        # There is nuance here:
        #    1) Each model has an input file which mirrors the functionality in
        #          the model. Further modifications to the model after the input
        #          file has been (essentially) manually modified can cause things
        #          to get out of sync and be messed up. This should be the last 
        #          thing done before the job is created and runs.
        #    2) We are really superimposing the user subroutine defined stress
        #          profile on the part which has undergone material removal via
        #          a boolean operation.
        assert(self.part.path_to_stress_subroutine != None)
        shim.modify_inp("*Initial Conditions", ("Type=Stress", "User"), "", model_name, mdb)

        job_name = tool_pass_name 
        job = shim.create_job(job_name, model_name, mdb) 

        # Associate the user subroutine with the job.
        shim.add_user_subroutine(job, self.part.path_to_stress_subroutine)

        shim.run_job(job)

        # Save off the resulting MDB.
        # The files produced by the job go in the working directory.
        shim.save_mdb(save_location, mdb)

        shim.close_mdb(mdb) 



    # Use multiple consecutive simulations to simulate the results of multiple
    #   consecutive tool passes. 
    # TODO: The current implementation chains simulations (1 simulation = 1 model) 
    #    within the working mdb.
    def sim_consecutive_tool_passes(self, cnt):
    # type: (int) -> None

         
        

    # TODO:
    # This does the backward direction (estimate residual stress which caused
    #   a particular deformation).
    def estimate_stress_via_last_tool_pass(self):
    # type: (bool) -> StressProfile
        
        pass




