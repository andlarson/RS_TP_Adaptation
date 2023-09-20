"""
This module contains code which can be used to keep track of metadata associated
   with abstractions which don't natively exist in Abaqus.
"""

import os

import core.metadata.abaqus_metadata as abq_md
import core.abaqus.abaqus_shim as shim


# Data structure associated with a single committed tool pass.
# Encapsulates all metadata associated with a single committed tool pass. This
#    includes all metadata associated with simulations which led to the decision
#    to commit to a particular tool pass.
class CommittedToolPassMetadata:

    def __init__(self, init_part, path_to_mdb, BCs):
    # type: (part.Part, str, List[bc.BC]) -> None

        # Idea: The metadata does not depend on the type of Part.
        self.init_part = init_part
        self.path_to_mdb = path_to_mdb
        self.working_dir = os.path.dirname(path_to_mdb)
        self.abaqus_mdb_metadata = abq_md.AbaqusMdbMetadata(path_to_mdb)
        self.committed_tool_pass = None
        self.simulated_tool_passes = []
        self.BCs = BCs



class RealWorldMetadata:

    def __init__(self):
    # type: (None) -> None

        raise RuntimeError("Not yet supported.")



# When creating a new simulation, some names need to be established.
# Note that some names are chosen, and others are looked up. Which names are
#    chosen versus looked up depends on the content of the MDB. 
class SimNames:

    def __init__(self, record):
    # type: (CommittedToolPassMetadata) -> None

        tool_pass_cnt = len(record.simulated_tool_passes)
        self.tool_pass_part_name = shim.STANDARD_TOOL_PASS_PART_PREFIX + str(tool_pass_cnt)
        self.post_tool_pass_part_name = shim.STANDARD_POST_TOOL_PASS_PART_PREFIX + str(tool_pass_cnt)
        step_cnt = 1
        self.equil_step_name = shim.STANDARD_EQUIL_STEP_PREFIX + str(step_cnt + 1)
        
        if tool_pass_cnt == 0:

            self.new_model_name = shim.STANDARD_MODEL_NAME
            self.pre_tool_pass_part_name = shim.STANDARD_INIT_GEOM_PART_NAME

        else:

            self.new_model_name = shim.STANDARD_MODEL_NAME_PREFIX + str(tool_pass_cnt + 1) 

            # Extract some names from the last model on record. 
            last_model_name = record.abaqus_mdb_metadata.model_names[-1] 
            # Abaqus makes the new part name all uppercase for whatever reason.
            self.pre_tool_pass_part_name = record.abaqus_mdb_metadata.models_metadata[last_model_name].part_names[-1].upper()     
            self.last_model_odb_file_name = record.abaqus_mdb_metadata.models_metadata[last_model_name].job_name + ".odb"