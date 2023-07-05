import core.part.part as part
import core.abaqus.abaqus_shim as shim
import core.material_properties as mat_props 
import core.simulation.simulation as sim
import core.tool_pass.tool_pass as tp
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
        # Note that we store paths to MDBs and open/close the MDBs themselves
        #    only as necessary.
        self.mdb_paths = []

        if isinstance(part, part.UserDefinedPart):
            # TODO: We don't handle this case right now...
            #       This should build an MDB and then append it to the list
            #          of simulations. Behavior should be symmetric to the
            #          else clause.
        else:
            assert(name == None)
            self.mdb_paths.append(part.path_to_mdb)


    # Simulate the next tool pass in the MDB most recently associated with
    #    this object.
    def sim_next_tool_pass(self):
    # type: (None) -> None
       
        mdb = shim.use_mdb(self.mdb_paths)

        # Build the next tool pass path as a part.
        tool_pass_geom = self.tool_pass_plan.pop()
        tool_pass_name = shim.STANDARD_TOOL_PASS_PART_PREFIX + str(self.tool_pass_plan.done_sofar())
        tool_pass_part = shim.build_part(tool_pass_name, tool_pass_geom, self.mdb)

        # Instance the tool pass part.
        tool_pass_part_instance = shim.instance_part_into_assembly(tool_pass_part, mdb, True)

        # Instance the initial geometry part.
        initial_geom_part = shim.get_part(shim.STANDARD_INIT_GEOM_PART_NAME, self.mdb)
        initial_geom_part_instance = shim.instance_part_into_assembly(initial_geom_part, self.mdb, True)

        # Create the post cut geometry as a part and instance it.
        post_tool_pass_name = shim.STANDARD_POST_TOOL_PASS_PART_PREFIX + str(self.tool_passes.cuts_done)
        cut_instances = [tool_pass_part_instances]
        post_tool_pass_part = shim.cut_instances_in_assembly(post_tool_pass_name, initial_geom_part_instance, cut_instances, self.mdb)

        # Add an equilibrium step following the last step on record.
        last_known_step = self.sim_metadata.get_last_step()
        step_cnt = shim.get_step_cnt(self.mdb)
        equil_step_name = shim.STANDARD_EQUIL_STEP_PREFIX + str(step_cnt + 1)
        shim.create_equilibrium_step(equil_step_name, last_known_step, self.mdb)

        # Then create a job and run the job.


        # And save off the result.
    
    
    # TODO:
    # Use multiple consecutive simulations to simulate the results of multiple
    #   consecutive tool passes. We could also implement a naive-approach which
    #   could simulate multiple tool passes with a single simulation with the
    #   implicit assumption that each tool path, except possibly the last one,
    #   would not cause a large deformation.
    def sim_consecutive_tool_passes(self, cnt):
    # type: (int) -> None

        pass 
       
         

    # This does the backward direction (estimate residual stress which caused
    #   a particular deformation).
    # If real-world data is available, that data can be used.
    def estimate_stress_via_last_tool_pass(self, use_real_meas_data):
    # type: (bool) -> StressProfile

        if use_real_meas_data and self.all_in_sim:
            raise RuntimeError("There is no data from the real world available \ 
                                for this machining process, so the stress \
                                estimation cannot use any real world data.")


    








