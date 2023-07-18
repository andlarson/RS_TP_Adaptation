import core.part.part as part
import core.abaqus.abaqus_shim as shim
import core.abaqus.abaqus_metadata as abq_md
import core.tool_pass.tool_pass_sandbox as tp_sandbox

import os
import os.path as path



class MachiningProcess:
   
    def __init__(self, name, init_part, tool_pass_plan):
    # type: (str, part.Part, tp.ToolPassPlan) -> None
    
        self.name = name
        self.part_history = part.PartHistory()
       
        """
        This list contains objects which each encapsulate all the data associated
           with a single tool pass.
        """
        self.single_tool_pass_records = []

        record = tp_sandbox.SingleToolPass(init_part, tool_pass_plan)
        self.single_tool_pass_records.append(record)  


    """
    The next tool pass has been decided upon.
    Either the next tool pass will be done in real life or it can be simulated
       if there is no real life data available.
    """
    def exec_tool_pass(self, tool_pass, real_life_data):
    # type: (tp.ToolPass, Any) -> None

        pass


    def sim_next_potential_tool_pass(self, save_path):
    # type: (str) -> None
       
        # Retrieve the last record on file. 
        last_record = self.single_tool_pass_records[-1]

        # Simulate the next tool pass and record as such in the record.
        last_record.sim_next_tool_pass(save_path)



    # Use multiple consecutive simulations to simulate the results of multiple
    #   consecutive tool passes. 
    # TODO: The current implementation chains simulations (1 simulation = 1 model) 
    #    within the working mdb.
    def sim_consecutive_tool_passes(self, cnt):
    # type: (int) -> None
    
        pass
         
        

    # TODO:
    # This does the backward direction (estimate residual stress which caused
    #   a particular deformation).
    def estimate_stress_via_last_tool_pass(self):
    # type: (bool) -> StressProfile
        
        pass





