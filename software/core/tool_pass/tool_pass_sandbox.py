import tool_pass as tp
import core.part.part as part
import core.abaqus.abaqus_metadata as abq_md
import core.abaqus.abaqus_shim as shim
import core.simulation.simulation as sim

import os.path


"""
Encapsulates all data associated with a particular executed (in real life or
   in simulation) tool pass.
The philosophy is that a single MDB is used for a single executed tool pass.
   It might be the case that multiple potential tool paths need to be simulated.
   In this case, all simulations live in this single MDB. When a single tool
   path is decided upon, the decision is stored in this object.
"""
class SingleToolPass:

    def __init__(self, init_part, tool_pass_plan):
    # type: (part.Part, ToolPassPlan) -> None
        
        self.tool_pass_plan = tool_pass_plan
        self.init_part = init_part

        if isinstance(init_part, part.UserDefinedPart):
            """
            TODO: We don't handle this case right now...
            """
            pass

        else:
            self.path_to_mdb = init_part.path_to_mdb
            self.abq_metadata = abq_md.ABQMetadata(self.path_to_mdb)

            """
            In general, all modifications to MDB metadata should happen in the
               abaqus_shim.py file at the same time as modifications are made
               to the MDB.
            This is an exceptional case because the user has passed an MDB,
               so it's necessary to sync the metadata with the content of the
               MDB.
            """
            self.abq_metadata.add_model(shim.STANDARD_MODEL_NAME)
            self.abq_metadata.add_part_to_model(shim.STANDARD_MODEL_NAME, shim.STANDARD_INIT_GEOM_PART_NAME)

            # Record the directory that the MDB lives in. 
            self.working_dir = os.path.dirname(self.path_to_mdb)
 

        
    def sim_next_tool_pass(self, save_path):
    # type: (str) -> None

        # The results of this simulation should go in the same directory
        #    that the MDB lives in.
        os.chdir(self.working_dir)

        mdb = shim.use_mdb(self.path_to_mdb)

        next_tool_pass = self.tool_pass_plan.pop()
        tool_pass_cnt = self.tool_pass_plan.done_so_far()

        """
        The way that the next tool pass is simulated depends on the content of
           the MDB.
        If the MDB contains a single model which represents an initial
           geometry, then it's necessary to simulate the first tool pass. 
        If the MDB contains n (where n > 1) models and the last model to be 
           added contains an orphan mesh, then we must be simulating the n-th 
           tool pass.
        """
        last_model_name = self.abq_metadata.get_last_model_name()
        if shim.check_init_geom(mdb):
            
            # If we have an initial geometry, there better be an associated
            #    stress profile.
            path_to_stress_subroutine = self.init_part.path_to_stress_subroutine
            assert(path_to_stress_subroutine != None)

            mdb = sim.sim_first_tool_pass(next_tool_pass, tool_pass_cnt, self.abq_metadata, path_to_stress_subroutine, mdb)  

        elif shim.check_multiple_steps(last_model_name, mdb):

            # If there are multiple steps in the last model, it is assumed that
            #    the model was run. Therefore, it's necessary to use the results
            #    of that model to chain the simulations together. 

            # Assumption here. The last job from the last model produced the 
            #    output we care about. Also assume the last part created in
            #    the last model is deformed, and is therefore responsible for
            #    the orphan mesh.
            last_odb_name = self.abq_metadata.get_last_job_name(last_model_name)
            last_part_name = self.abq_metadata.get_last_part_name(last_model_name)
            mdb = sim.sim_nth_tool_pass(next_tool_pass, tool_pass_cnt, last_part_name, last_odb_name, mdb)

        else:
            raise RuntimeError("Can't figure out how to do next tool pass...")

        shim.save_mdb_as(save_path, mdb)
        shim.close_mdb(mdb)



    def sim_consecutive_tool_passes(self, save_path):
    # type: (str) -> None
        


                
         
