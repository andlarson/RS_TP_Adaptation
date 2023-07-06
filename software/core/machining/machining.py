import core.part.part as part
import core.abaqus.abaqus_shim as shim
import core.material_properties as mat_props 
import core.simulation.simulation as sim
import core.tool_pass.tool_pass as tp
import core.abaqus.metadata as md
import storage.storage as storage

import os




class MachiningProcess:
   
    def __init__(self, name, part, tool_pass_plan):
    # type: (str, part.Part, tp.ToolPassPlan) -> None
    
        self.name = name
        self.part = part
        self.tool_pass_plan = tool_pass_plan 
        self.part_history = part.PartHistory()

        # A single machining process may contain many model databases (i.e. MDBs).
        #    This gives the most flexibility. We may want to chain simulations
        #    in a single MDB or do a single simulation per MDB to make things as
        #    straightforward as possible.
        # This maps MDB paths to metadata objects.
        self.metadata = {}

        if isinstance(part, part.UserDefinedPart):
            # TODO: We don't handle this case right now...
            #       This should build an MDB and then append it to the list
            #          of simulations. Behavior should be symmetric to the
            #          else clause.
        else:
            assert(name == None)

            # Make sure that things are as we expect.
            mdb = shim.use_mdb(path_to_mdb)
            shim.verify_init_geom_mdb(mdb)
            shim.close_mdb(mdb)

            self.metadata[part.path_to_mdb] = md.MDBMetadata(part.path_to_mdb, md.STANDARD_MODEL_NAME)
            self.working_mdb_path = part.path_to_mdb

    
    def get_working_model_name(self):
    # type: (None) -> None

        return self.metadata[working_mdb_path].get_working_model_name()


    def sim_next_tool_pass(self):
    # type: (None) -> None
       
        mdb = shim.use_mdb(self.working_mdb_path)
        model_name = get_working_model_name(self)

        # Build the next tool pass path as a part.
        tool_pass_geom = self.tool_pass_plan.pop()
        tool_pass_name = shim.STANDARD_TOOL_PASS_PART_PREFIX + str(self.tool_pass_plan.done_sofar())
        tool_pass_part = shim.build_part(tool_pass_name, tool_pass_geom, model_name, mdb)

        # Instance the tool pass part.
        tool_pass_part_instance = shim.instance_part_into_assembly(tool_pass_part, True, model_name, mdb)

        # Instance the initial geometry part.
        # Assumed to already exist in the MDB.
        # TODO: Is this name convention always followed?
        initial_geom_part = shim.get_part(shim.STANDARD_INIT_GEOM_PART_NAME, mdb)
        initial_geom_part_instance = shim.instance_part_into_assembly(initial_geom_part, True, model_name, mdb)

        # Create the post cut geometry as a part.
        post_tool_pass_name = shim.STANDARD_POST_TOOL_PASS_PART_PREFIX + str(self.tool_passes.cuts_done)
        cut_instances = [tool_pass_part_instance]
        post_tool_pass_part = shim.cut_instances_in_assembly(post_tool_pass_name, initial_geom_part_instance, cut_instances, model_name, mdb)

        # Instance the post cut geometry.
        post_tool_pass_instance = shim.intsnace_part_into_assembly(post_tool_pass_part, True, model_name, mdb)

        # Then mesh the part in the assembly module.
        shim.naive_mesh(post_tool_pass_instance)

        # Add an equilibrium step following the last step on record.
        last_step = self.metadata[self.working_mdb_path].get_last_step()
        step_cnt = shim.get_step_cnt(model_name, mdb)
        equil_step_name = shim.STANDARD_EQUIL_STEP_PREFIX + str(step_cnt + 1)
        shim.create_equilibrium_step(equil_step_name, last_step, model_name, mdb)

        # Then create and run a job.
        job_name = tool_pass_name 
        shim.create_and_run_job(job_name, model_name, mdb) 

        # And save off the result.
        shim.save_mdb(mdb)

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
    








