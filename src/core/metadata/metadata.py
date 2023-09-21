"""
This module contains code which can be used to keep track of metadata associated
   with abstractions which don't natively exist in Abaqus. Abaqus metadata is
   embedded in this metadata.
"""

import core.metadata.abaqus_metadata as abq_md
import core.abaqus.abaqus_shim as shim

from util.debug import *


# Data structure associated with a single committed tool pass.
# Encapsulates all metadata associated with a single committed tool pass. This
#    includes all metadata associated with simulations which led to the decision
#    to commit to a particular tool pass.
class CommittedToolPassMetadata:

    def __init__(self, init_part, path_to_mdb, BCs):
    # type: (part.Part, str, List[bc.BC]) -> None

        self.init_part = init_part
        self.path_initial_mdb = path_to_mdb

        # Could be dictionary, mapping mdb objects to metadata data structures.
        # Doesn't seem to matter - we are always working with the most recent MDB.
        self.per_mdb_metadata = [] 

        self.committed_tool_pass = None
        self.BCs = BCs



class RealWorldMetadata:

    def __init__(self):
    # type: (None) -> None

        raise RuntimeError("Not yet supported.")



# Figure out how to name all the stuff included in a new Abaqus model.
#
# Notes:
#    Some names are chosen based on convention. Some are looked up because the
#       thing already exists in the Abaqus model. For example, when an MDB is
#       created it comes with a single model which has a single part in it. These
#       things already have names.
#
# Arguments:
#    mdb_metadata - AbaqusMdbMetadata object.
#    is_initial   - Boolean.
#                   Should be True when the first model in an MDB is being worked
#                      with.
#                   Should be False when some model beyond the first is being
#                      created.
#
# Returns:
#    Dict. See code for key-val pairings.
def new_model_names(mdb_metadata, is_initial):
# type: (abq_md.AbaqusMdbMetadata, bool) -> dict[str, str]

    names = {}

    model_cnt = len(mdb_metadata.model_names)

    names["tool_pass_part_name"] = shim.STANDARD_TOOL_PASS_PART_PREFIX + str(model_cnt)
    names["post_tool_pass_part_name"] = shim.STANDARD_POST_TOOL_PASS_PART_PREFIX + str(model_cnt)
    step_cnt = 1
    names["equil_step_name"] = shim.STANDARD_EQUIL_STEP_PREFIX + str(step_cnt + 1)
    
    if is_initial:

        names["new_model_name"] = shim.STANDARD_MODEL_NAME
        names["pre_tool_pass_part_name"] = shim.STANDARD_INIT_GEOM_PART_NAME

    else:

        names["new_model_name"] = shim.STANDARD_MODEL_NAME_PREFIX + str(model_cnt + 1) 

        # Figure out the last model in this MDB. 
        last_model_name = mdb_metadata.model_names[-1] 
        # Abaqus makes the new part name all uppercase for whatever reason.
        names["pre_tool_pass_part_name"] = mdb_metadata.models_metadata[last_model_name].part_names[-1].upper()     

    return names



# Figure out the name of the output database file which resulted from the last
#    job running.
#
# Notes:
#    It only makes sense to call this function if a job has already run!
#
# Arguments:
#    mdb_metadata - AbaqusMdbMetadata object.
#
# Returns:
#    String. End with .odb suffix. 
def last_odb_file_name(mdb_metadata):
# type: (abq_md.AbaqusMdbMetadata) -> str

    # Figure out the last model name in this MDB. 
    last_model_name = mdb_metadata.model_names[-1] 

    return mdb_metadata.models_metadata[last_model_name].job_name + ".odb"



# Figure out the name of the sim file which resulted from the last job running.
#
# Notes:
#    It only makes sense to call this function if a job has already run!
#
# Arguments:
#    mdb_metadata - AbaqusMdbMetadata object.
#
# Returns:
#    String. End with .sim suffix. 
def last_sim_file_name(mdb_metadata):
# type: (abq_md.AbaqusMdbMetadata) -> str

    # Figure out the last model name in this MDB. 
    last_model_name = mdb_metadata.model_names[-1] 

    return mdb_metadata.models_metadata[last_model_name].job_name + ".sim"