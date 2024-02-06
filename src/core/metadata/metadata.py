"""
This module contains functionality which can be used to keep track of and retrieve 
    metadata associated with abstractions which don't natively exist in Abaqus. 
    Abaqus metadata is embedded in this metadata.
"""

import src.core.part.part as part
import src.core.boundary_conditions.boundary_conditions as bc

from src.util.debug import *


class CommittedToolPassPlanMetadata:

    def __init__(self, init_part: part.AbaqusDefinedPart, path_to_mdb: str, 
                 BCs: list[bc.BC]) -> None:
        """Creates a data structure associated with a single committed tool pass plan.

           The general system flow is assumed to be something like:
               1) User simulates a bunch of potential tool passes.
               2) Based on simulation results, user decides which tool pass to 
                      perform in real life.
               3) Repeat steps 1 and 2 until done machining.
           Thus, there are natural phases to the procedure. A phase corresponds
               to all the stuff done between the last tool pass in real life and
               the next tool pass in real life. This data structure records information
               about what occurred during a phase.

           Args:
               init_part:   The part at the beginning of this phase. Note that 
                                this part may have had some material removed 
                                from it already.
               path_to_mdb: The path to the MDB which contains the state of the 
                                part at the beginning of this commitment phase.
               BCs:         The boundary conditions active during this phase.

           Returns:
               None.

           Raises:
               None.
        """
        
        # Flag indicating if this commitment phase is the very first one.
        self.first_commitment_phase = False

        self.init_part = init_part
        self.path_initial_mdb = path_to_mdb
        self.BCs = BCs

        # The tool pass plan which was committed to for this commitment phase.
        # This acts as a flag. If this is not None, then the commitment phase
        #     is done.
        self.committed_tool_pass_plan = None

        # Absolute path to the directory containing the results of the simulated
        #     tool pass plan.
        self.committed_tool_pass_plan_path = None

        # There is a single MDB for each tool pass plan simulated in this
        #     commitment phase. This keeps track of the metadata for each of these
        #     MDBs.
        # Note that a commitment phase includes both potential tool pass plans
        #     and a single committed tool pass plan.
        self.per_mdb_metadata = [] 

        # Optimization: To speed up committing a tool pass plan, a list of all
        #     the tool pass plans which have been simulated is maintained. That
        #     way, it's possible to check if the tool pass plan has already been
        #     simulated before commitment is done. If so, no need to re-simulate.   
        # Each list entry is a Tuple which contains a name, a ToolPassPlan object,
        #     and the absolute path to the directory where the simulation results
        #     are stored.
        self.simulated_tool_pass_plans = [] 

        # Store the path to the .sim file for the last simulation in the last
        #     commitment phase.
        # This is helpful to keep track of because it can be used as the starting
        #     stress profile for the current commitment phase.
        self.path_last_commit_sim_file = None



