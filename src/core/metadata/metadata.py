"""
This module contains code which can be used to keep track of and retrieve 
   metadata associated with abstractions which don't natively exist in Abaqus. 
   Abaqus metadata is embedded in this metadata.
"""

from util.debug import *


# Data structure associated with a single committed tool pass plan.
# Encapsulates all metadata associated with a single committed tool pass plan. This
#    includes all metadata associated with simulations which led to the decision
#    to commit to the tool pass plan (i.e. exploratory tool pass simulations).
class CommittedToolPassPlanMetadata:

    def __init__(self, init_part, path_to_mdb, BCs):
    # type: (part.Part, str, List[bc.BC]) -> None

        self.init_part = init_part
        self.path_initial_mdb = path_to_mdb

        # Could be dictionary, mapping mdb objects to metadata data structures.
        # Doesn't seem to matter - we always work with the most recent MDB.
        self.per_mdb_metadata = [] 

        # Optimization: To speed up committing a tool pass plan, a list of all
        #    the tool pass plans which have been simulated is maintained. That
        #    way, it's possible to check if the tool pass plan has already been
        #    simulated before commitment is done. If so, no need to re-simulate.   
        self.simulated_tool_pass_plans = []

        self.committed_tool_pass_plan = None
        self.BCs = BCs