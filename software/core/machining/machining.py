import core.part.part as part
import core.abaqus.abaqus_shim as shim
import core.abaqus.abaqus_metadata as abq_md
import core.tool_pass.tool_pass_sandbox as tp_sandbox

import os
import os.path as path



class MachiningProcess:
   
    def __init__(self, init_part):
    # type: (part.Part) -> None
    
        self.part_history = part.PartHistory()
       
        """
        This list contains objects which each encapsulate all the data associated
           with a single tool pass.
        """
        self.tool_pass_records = []


    
    # When a tool pass is committed, it means that it has been decided upon.
    # If the machining process is happening in real life, then there may be
    #    some measurement data collected after the tool pass happened.
    def commit_tool_pass(self, tool_pass, in_real_life, measurement_data):
    # type: (tp.ToolPass, bool, Any) -> None

        pass



    # Simulate tool passes without committing to anything. 
    def sim_potential_tool_passes(self, tool_pass_plan, save_path):
    # type: (tp.ToolPassPlan, str) -> None
       
        # Retrieve the last record on file. 
        last_record = self.tool_pass_records[-1]

        # Simulate the next tool pass and record as such in the record.
        last_record.sim_next_tool_pass(save_path)



    # TODO:
    # This does the backward direction (estimate residual stress which caused
    #   a particular deformation).
    def estimate_stress_via_last_tool_pass(self):
    # type: (bool) -> StressProfile
        
        pass





