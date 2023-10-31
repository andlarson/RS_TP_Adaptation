"""
This module contains code which generates and retrieves names.
"""

import core.abaqus.abaqus_shim as shim



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
#                   Should be True when generating names for the first model in
#                      an MDB, and the model already exists.
#                   Should be False when generating names for some model beyond
#                      the first in an MDB, and that model does not yet exist.
#
# Returns:
#    Dict. See code for key-val pairings.
def new_model_names(mdb_metadata, is_initial):
# type: (abq_md.AbaqusMdbMetadata, bool) -> dict[str, str]

    names = {}

    model_cnt = len(mdb_metadata.model_names)

    step_cnt = 1
    names["equil_step_name"] = shim.STANDARD_EQUIL_STEP_PREFIX + str(step_cnt + 1)
    names["bounding_box_name"] = shim.STANDARD_BOUNDING_BOX_PART_NAME
    names["excess_bounding_box_name"] = shim.STANDARD_EXCESS_BOUNDING_BOX_PART_NAME
    names["merged_bbox_geom_name"] = shim.STANDARD_BOUNDING_BOX_INIT_GEOM_PART_NAME
    names["init_geom_with_bbox_name"] = shim.STANDARD_INIT_GEOM_WITH_BBOX_NAME
    
    if is_initial:

        names["post_tool_pass_part_name"] = shim.STANDARD_POST_TOOL_PASS_PART_PREFIX + str(model_cnt)
        names["new_model_name"] = shim.STANDARD_MODEL_NAME
        names["pre_tool_pass_part_name"] = shim.STANDARD_INIT_GEOM_PART_NAME
        names["tool_pass_part_name"] = shim.STANDARD_TOOL_PASS_PART_PREFIX + str(model_cnt)

    else:

        names["post_tool_pass_part_name"] = shim.STANDARD_POST_TOOL_PASS_PART_PREFIX + str(model_cnt + 1)
        names["new_model_name"] = shim.STANDARD_MODEL_NAME_PREFIX + str(model_cnt + 1) 
        names["tool_pass_part_name"] = shim.STANDARD_TOOL_PASS_PART_PREFIX + str(model_cnt + 1)

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