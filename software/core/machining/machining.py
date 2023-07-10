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
            self.mdb_to_metadata[init_part.path_to_mdb] = md.MDBMetadata(init_part.path_to_mdb, shim.STANDARD_MODEL_NAME)
            self.working_mdb_path = init_part.path_to_mdb

    
    def get_working_model_name(self):
    # type: (None) -> None

        return self.mdb_to_metadata[self.working_mdb_path].get_working_model_name()


    def sim_next_tool_pass(self, save_location):
    # type: (str) -> None
       
        mdb = shim.use_mdb(self.working_mdb_path)
        model_name = self.get_working_model_name()
        metadata = self.mdb_to_metadata(self.working_mdb_path)
       
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

        # Create the post cut geometry as a part.
        post_tool_pass_name = shim.STANDARD_POST_TOOL_PASS_PART_PREFIX + str(self.tool_pass_plan.done_so_far())
        cut_instances = (tool_pass_part_instance, )
        post_tool_pass_part = shim.cut_instances_in_assembly(post_tool_pass_name, initial_geom_instance, cut_instances, model_name, mdb)

        # Instance the post cut geometry.
        # The instance is independent so that meshing can be done in the Assembly
        #    module.
        post_tool_pass_instance = shim.instance_part_into_assembly(post_tool_pass_name, post_tool_pass_part, False, model_name, mdb)

        # Then mesh the part in the assembly module.
        # Picking a pretty fine seed density heuristically....
        shim.naive_mesh(post_tool_pass_instance, 1, model_name, mdb)

        # Add an equilibrium step following the last step on record.
        last_step = self.mdb_to_metadata[self.working_mdb_path].get_last_step(model_name)
        step_cnt = shim.get_step_cnt(model_name, mdb)
        equil_step_name = shim.STANDARD_EQUIL_STEP_PREFIX + str(step_cnt + 1)
        shim.create_equilibrium_step(equil_step_name, last_step, model_name, metadata, mdb)

        # Then create and run a job.
        job_name = tool_pass_name 
        shim.create_and_run_job(job_name, model_name, mdb) 

        # And save off the result.
        # TODO: Need a general solution here. Maybe access the MDB metadata or
        #    use the storage functionality to determine where to save.
        shim.save_mdb(save_location, mdb)

        shim.close_mdb(mdb) 



    # TODO:
    # Use multiple consecutive simulations to simulate the results of multiple
    #   consecutive tool passes. We could also implement a naive-approach which
    #   could simulate multiple tool passes with a single simulation with the
    #   implicit assumption that each tool path, except possibly the last one,
    #   would not cause a large deformation.
    def sim_consecutive_tool_passes(self, cnt):
    # type: (int) -> None

        pass 
       
         
    # TODO:
    # This does the backward direction (estimate residual stress which caused
    #   a particular deformation).
    def estimate_stress_via_last_tool_pass(self, use_real_meas_data):
    # type: (bool) -> StressProfile
        
        pass
    








