"""
This module contains functionality to generate and retrieves names.
"""
from enum import Enum

import src.core.abaqus.abaqus_shim as shim
import src.core.metadata.abaqus_metadata as abq_md


# Different types of models have different naming schemes.
class ModelTypes(Enum):
    FIRST_TOOL_PASS_IN_MDB = 1
    NTH_TOOL_PASS_IN_MDB = 2
    FIRST_TRACTION_APP_IN_MDB = 3
    NTH_TRACTION_APP_IN_MDB = 4



class ModelNames:

    def __init__(self, 
                 model_type: ModelTypes, 
                 mdb_metadata: abq_md.AbaqusMdbMetadata, 
                ) -> None:
        """Generates new names and determines pre-existing names for stuff in
               a particular model in an MDB.
           An instance is typically created right before accessing, or adding
               to, the content of a new model in an MDB.
            
           Arguments:
              model_type:   The type of the model for which names need to be
                                generated.
                            For example, if the model is the first in an MDB
                                and is going to be used to simulate a tool pass, 
                                then the names of some stuff are set by Abaqus
                                (e.g. the model name) while other names are
                                set by the user (e.g. the first part they create).
              mdb_metadata: The metadata associated with the MDB of interest. 
           
           Returns:
               None.
           
           Raises:
               None.
        """

        model_cnt = len(mdb_metadata.model_names)
        
        if model_type is ModelTypes.FIRST_TOOL_PASS_IN_MDB:
            # Already exists.
            self.new_model_name = shim.STANDARD_MODEL_NAME

            self.pre_tool_pass_part_name = shim.STANDARD_INIT_GEOM_PART_NAME
            self.tool_pass_part_name = shim.STANDARD_TOOL_PASS_PART_PREFIX 
            self.post_tool_pass_part_name = shim.STANDARD_POST_TOOL_PASS_PART_PREFIX 
            self.equil_step_name = shim.STANDARD_EQUIL_STEP_PREFIX 

        elif model_type is ModelTypes.NTH_TOOL_PASS_IN_MDB:
            self.post_tool_pass_part_name = shim.STANDARD_POST_TOOL_PASS_PART_PREFIX 
            self.new_model_name = shim.STANDARD_MODEL_NAME_PREFIX + str(model_cnt + 1) 
            self.tool_pass_part_name = shim.STANDARD_TOOL_PASS_PART_PREFIX 
            self.equil_step_name = shim.STANDARD_EQUIL_STEP_PREFIX 

            # Figure out the last model in this MDB. 
            last_model_name = mdb_metadata.model_names[-1]  
            # Assumes that an ODB was used to create the initial part in this
            #     model, and that ODB was produced by running the last
            #     model in this MDB. When an ODB is imported to create a part,
            #     the name of the part is, by default, in all caps.
            self.pre_tool_pass_part_name = mdb_metadata.models_metadata[last_model_name].part_names[-1].upper()     

        elif model_type is ModelTypes.FIRST_TRACTION_APP_IN_MDB:
            self.new_model_name = shim.STANDARD_MODEL_NAME
            self.deformed_part_name = shim.SOURCED_FROM_ODB_PART_NAME
            self.traction_step_name = shim.STANDARD_TRACTION_STEP_PREFIX 

        elif model_type is ModelTypes.NTH_TRACTION_APP_IN_MDB:
            self.new_model_name = shim.STANDARD_MODEL_NAME + str(model_cnt)
            self.deformed_part_name = shim.SOURCED_FROM_ODB_PART_NAME
            self.traction_step_name = shim.STANDARD_TRACTION_STEP_PREFIX 



def last_odb_file_name(mdb_metadata: abq_md.AbaqusMdbMetadata) -> str:
    """Determines the name of the output database (.odb file) which resulted
           from the last model being run in a particular MDB.
       
       It only makes sense to call this function if a job has already run!
           This function does not actually check if a job was run, it just
           returns the name of the .odb which would be produced if the job
           were run.

       Args:
           mdb_metadata: Metadata associated with this MDB.

       Returns:
           The name of the .odb file associated with the job of the last model
               which ran.

       Raises:
           RuntimeError: No job was run in the last model, so no odb file exists.
    """

    # Figure out the last model name in this MDB. 
    last_model_name = mdb_metadata.model_names[-1] 

    if mdb_metadata.models_metadata[last_model_name].job_name:
        return mdb_metadata.models_metadata[last_model_name].job_name + ".odb"

    raise RuntimeError("No job recorded, so no odb file exists!")



def last_sim_file_name(mdb_metadata):
    """Determines the name of the sim database (.sim file) which resulted
           from the last model being run in a particular MDB.
       
       It only makes sense to call this function if a job has already run!
           This function does not actually check if a job was run, it just
           returns the name of the .sim which would be produced if the job
           were run.

       Args:
           mdb_metadata: Metadata associated with this MDB.

       Returns:
           The name of the .sim file associated with the job of the last
               model which ran.

       Raises:
           RuntimeError: No job was run in the last model, so no sim file exists.
    """

    # Figure out the last model name in this MDB. 
    last_model_name = mdb_metadata.model_names[-1] 

    if mdb_metadata.models_metadata[last_model_name].job_name:
        return mdb_metadata.models_metadata[last_model_name].job_name + ".sim"

    raise RuntimeError("No job recorded, so no sim file exists!")

