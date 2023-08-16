import core.tool_pass.tool_pass as tp
import core.part.part as part
import core.abaqus.abaqus_metadata as abq_md
import core.abaqus.abaqus_shim as shim
import core.boundary_conditions.boundary_conditions as bc

import os.path


"""
Encapsulates all data associated with a single committed tool pass.
The philosophy is that a single MDB is used for a single committed tool pass.
   It might be the case that multiple potential tool paths need to be simulated
   before the user commits to a single tool pass. In this case, all simulations 
   live in this single MDB. When a single tool path is decided upon and committed
   to, the decision is recorded in this object.
"""
class CommittedToolPassRecord:

    def __init__(self, init_part, boundary_conditions):
    # type: (part.Part, List[bc.BC]) -> None
        
        self.init_part = init_part
        self.BCs = boundary_conditions

        if isinstance(init_part, part.UserDefinedPart):
           raise RuntimeError("Can't handle user defined parts right now...") 
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
 
        
