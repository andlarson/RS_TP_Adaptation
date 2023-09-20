import core.part.part as part
import core.simulation.simulation as sim
import core.boundary_conditions.boundary_conditions as bc
import core.metadata.metadata as md 


class MachiningProcess:
   
    def __init__(self, init_part, boundary_conditions):
    # type: (part.Part, List[bc.BC]) -> None
  
        self.metadata = []
        if isinstance(init_part, part.UserDefinedPart):
            # When this is supported, it's expected that a UserDefinedPart is
            #    first constructed in an MDB so the method of setting up the
            #    metadata is the same as the AbaqusDefinedPart.
            raise RuntimeError("Not yet supported.")

        elif isinstance(init_part, part.AbaqusDefinedPart):
            first_tp_metadata = md.CommittedToolPassMetadata(init_part, init_part.path_to_mdb, boundary_conditions)
            self.metadata.append(first_tp_metadata)

        # The boundary conditions are assumed to exist for the whole machining process.
        self.boundary_conditions = boundary_conditions 



    def commit_tool_pass(self, tool_pass, in_real_life, measurement_data):
    # type: (tp.ToolPass, bool, Any) -> None

        raise RuntimeError("Not yet supported.") 



    # Simulate tool passes without committing.
    def sim_potential_tool_passes(self, tool_pass_plan, save_path):
    # type: (tp.ToolPassPlan, str) -> None
       
        sim.sim_consecutive_tool_passes(tool_pass_plan, save_path, self.metadata[-1])



    def estimate_stress_via_last_tool_pass(self):
    # type: (None) -> None
        
        raise RuntimeError("Not yet supported.")