"""
This module contains functionality to generate and retrieves names.
"""

import src.core.abaqus.abaqus_shim as shim
import src.core.metadata.abaqus_metadata as abq_md


class ModelNames:

    def __init__(self, mdb_metadata: abq_md.AbaqusMdbMetadata, 
                 is_initial: bool
                ) -> None:
        """Creates the names associated with an MDB. Uses the metadata for the
               MDB to determine some of the names.
        
           Chooses some names based on convention. Looks up some names because 
               the thing already exists in the Abaqus model. For example, when 
               an MDB is created it comes with a single model which has a single 
               part in it. These things already have names.
           
           Arguments:
              mdb_metadata: The metadata associated with the MDB of interest. 
              is_initial:   Should be True when deciding names for the first 
                                model in an MDB, and the model already exists. 
                                Should be False when generating names for some 
                                model beyond the first in an MDB, and that 
                                model does not yet exist.
           
           Returns:
               None.
           
           Raises:
               None.
        """

        self.model_cnt = len(mdb_metadata.model_names)
        
        # By default, a model starts with a single step in it. 
        step_cnt = 1

        self.equil_step_name = shim.STANDARD_EQUIL_STEP_PREFIX + str(step_cnt + 1)
        
        if is_initial:
            self.post_tool_pass_part_name = shim.STANDARD_POST_TOOL_PASS_PART_PREFIX + str(self.model_cnt)
            self.new_model_name = shim.STANDARD_MODEL_NAME
            self.pre_tool_pass_part_name = shim.STANDARD_INIT_GEOM_PART_NAME
            self.tool_pass_part_name = shim.STANDARD_TOOL_PASS_PART_PREFIX + str(self.model_cnt)

        else:
            self.post_tool_pass_part_name = shim.STANDARD_POST_TOOL_PASS_PART_PREFIX + str(self.model_cnt + 1)
            self.new_model_name = shim.STANDARD_MODEL_NAME_PREFIX + str(self.model_cnt + 1) 
            self.tool_pass_part_name = shim.STANDARD_TOOL_PASS_PART_PREFIX + str(self.model_cnt + 1)

            # Figure out the last model in this MDB. 
            last_model_name = mdb_metadata.model_names[-1] 
            # Abaqus makes the new part name all uppercase for whatever reason.
            self.pre_tool_pass_part_name = mdb_metadata.models_metadata[last_model_name].part_names[-1].upper()     



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

