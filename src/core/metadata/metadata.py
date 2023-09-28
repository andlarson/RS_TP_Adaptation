"""
This module contains code which can be used to keep track of and retrieve 
   metadata associated with abstractions which don't natively exist in Abaqus. 
   Abaqus metadata is embedded in this metadata.
"""

from util.debug import *


# Data structure associated with a single committed tool pass plan.
# Encapsulates all metadata associated with a single committed tool pass plan. This
#    includes metadata associated with simulations which led to the decision
#    to commit to the tool pass plan.
class CommittedToolPassPlanMetadata:

    def __init__(self, init_part, path_to_mdb, BCs):
    # type: (part.Part, str, List[bc.BC]) -> None

        self.init_part = init_part
        self.path_initial_mdb = path_to_mdb
        self.committed_tool_pass_plan = None
        self.BCs = BCs

        # There is a single MDB for each tool pass plan simulated in this
        #    commitment phase. This keeps track of the metadata for each of these
        #    MDBs.
        self.per_mdb_metadata = [] 

        # Optimization: To speed up committing a tool pass plan, a list of all
        #    the tool pass plans which have been simulated is maintained. That
        #    way, it's possible to check if the tool pass plan has already been
        #    simulated before commitment is done. If so, no need to re-simulate.   
        # Each list entry is a Tuple which contains a name and a ToolPassPlan
        #    object.
        self.simulated_tool_pass_plans = [] 

        # Store the path to the .sim file for the last simulation in the last
        #    commitment phase.
        # This is helpful to keep track of because it can be used as the starting
        #    stress profile for the current commitment phase.
        self.path_last_commit_sim_file = None