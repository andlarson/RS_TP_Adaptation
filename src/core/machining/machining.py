import core.part.part as part
import core.simulation.simulation as sim
import core.boundary_conditions.boundary_conditions as bc
import core.metadata.metadata as md 
import core.metadata.abaqus_metadata as abq_md


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

        # The boundary conditions are assumed to be fixed for the whole machining process. 
        self.boundary_conditions = boundary_conditions 



    def commit_tool_pass(self, tool_pass, in_real_life, measurement_data):
    # type: (tp.ToolPass, bool, Any) -> None

        raise RuntimeError("Not yet supported.") 



    # Simulate some potential tool passes and save off the results. 
    #
    # Notes:
    #    The MDB which resulted from the last committed tool pass is the starting
    #       point for these tool path simulations (or the first MDB for this
    #       MachiningProcess object). That MDB must exist in some directory. 
    #       This function will save the resulting simulation artifacts (the .odb, 
    #       .sim, .inp, etc. and the final .cae file) in a subdirectory of the
    #       directory which this MDB exists in. In this way, each call to this
    #       funciton results in a new MDB being created.
    #
    # Arguments:
    #    tool_pass_plan - ToolPassPlan object.
    #    name           - String.
    #                     The name of the subdirectory in which the simulation artifacts
    #                        will be placed. Also, the name of the MDB which results
    #                        from these simulations and lives in the subdirectory. 
    #
    # Returns:
    #    None. 
    def sim_potential_tool_passes(self, tool_pass_plan, save_name):
    # type: (tp.ToolPassPlan, str) -> None

        # The start point for each sequence of simulations is always the MDB which
        #    resulted from the last committed tool pass (or the very initial MDB).
        path_to_mdb = self.metadata[-1].path_initial_mdb

        # A new MDB is created for this sequence of simulations. Therefore, a
        #    new metadata data structure must exist and accompany this new MDB.
        self.metadata[-1].per_mdb_metadata.append(abq_md.AbaqusMdbMetadata(path_to_mdb))

        sim.sim_consecutive_tool_passes(tool_pass_plan, save_name, self.metadata[-1])



    def estimate_stress_via_last_tool_pass(self):
    # type: (None) -> None
        
        raise RuntimeError("Not yet supported.")