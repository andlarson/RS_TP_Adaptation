import core.part.part as part
import core.abaqus.abaqus_shim as shim
import core.abaqus.abaqus_metadata as abq_md
import core.tool_pass.tool_pass_record_keeping as tprk
import core.simulation.simulation as sim

import os



class MachiningProcess:
   
    def __init__(self, init_part):
    # type: (part.Part) -> None
   
        """
        It is necessary to keep track of things at two levels.
        The part history tracks part changes across committed tool passes.
           No information about attempted tool passes is maintained.
        Each tool pass record maps to a single committed tool pass. The
           record encapsulates all information which led to the decision to
           commit to the tool pass. The primary purpose of the record is to
           supplement the state that Abaqus maintains about an MDB.
        """
        self.part_history = part.PartHistory()
        self.tool_pass_records = []

        self.init_part = init_part


    """ 
    When a tool pass is committed, it means that it has been decided upon.
    If the machining process is happening in real life, then there may be
       some measurement data collected after the tool pass happened.
    """
    def commit_tool_pass(self, tool_pass, in_real_life, measurement_data):
    # type: (tp.ToolPass, bool, Any) -> None

        raise RuntimeError("Not yet supported.") 



    # Simulate tool passes without committing.
    def sim_potential_tool_passes(self, tool_pass_plan, save_path):
    # type: (tp.ToolPassPlan, str) -> None
        
        if len(self.tool_pass_records) == 0:
            self.tool_pass_records.append(tprk.CommittedToolPassRecord(self.init_part))

        last_record = self.tool_pass_records[-1]

        sim.sim_consecutive_tool_passes(tool_pass_plan, save_path, last_record)



    # This does the backward direction (estimate residual stress which caused
    #   a particular deformation).
    def estimate_stress_via_last_tool_pass(self):
    # type: (bool) -> StressProfile
        
        raise RuntimeError("Not yet supported.")





