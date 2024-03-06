"""
This module contains functionality to generate and retrieves names.
"""
from enum import Enum

import src.core.abaqus.abaqus_shim as shim
import src.core.metadata.abaqus_metadata as abq_md


# Different types of models have different naming schemes.
class ModelTypes(Enum):
    DEFAULT_MODEL = 0
    FIRST_TOOL_PASS = 1
    NTH_TOOL_PASS = 2
    TRACTION_APP = 3
    EVAL_VOL_DIFF = 4



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
                                (e.g. the model name) while other names 
                                get to be chosen (e.g. when a part is created).
              mdb_metadata: The metadata associated with the MDB of interest. 
           
           Returns:
               None.
           
           Raises:
               None.
        """

        model_cnt = len(mdb_metadata.model_names)
        
        if model_type is ModelTypes.DEFAULT_MODEL:
            # Assumed to already exist at call time. 
            self.new_model_name = shim.STANDARD_MODEL_NAME
            self.part_from_odb_name = shim.STANDARD_INIT_GEOM_PART_NAME

        elif model_type is ModelTypes.FIRST_TOOL_PASS:
            # Assumed to already exist at call time. 
            self.new_model_name = shim.STANDARD_MODEL_NAME
            self.pre_tool_pass_part_name = shim.STANDARD_INIT_GEOM_PART_NAME
            
            # Assumed to be chosen. 
            self.tool_pass_part_name = shim.STANDARD_TOOL_PASS_PART_PREFIX 
            self.post_tool_pass_part_name = shim.STANDARD_POST_TOOL_PASS_PART_PREFIX + str(model_cnt)
            self.deformation_step_name = shim.STANDARD_DEFORMATION_STEP_NAME
            self.equilibrium_step_name = shim.STANDARD_EQUIL_STEP_NAME 

        elif model_type is ModelTypes.NTH_TOOL_PASS:
            # Assumed to already exist at call time. 
            # Assumes that an ODB was used to create the initial part in this
            #     model, and that ODB was produced by running the last
            #     model in this MDB. When an ODB is imported to create a part,
            #     the name of the part is, by default, in all caps.
            last_model_name = mdb_metadata.model_names[-1]  
            self.pre_tool_pass_part_name = mdb_metadata.models_metadata[last_model_name].part_names[-1].upper()     

            # Assumed to be chosen. 
            self.post_tool_pass_part_name = shim.STANDARD_POST_TOOL_PASS_PART_PREFIX + str(model_cnt + 1)
            self.new_model_name = shim.STANDARD_MODEL_NAME_PREFIX + str(model_cnt + 1) 
            self.tool_pass_part_name = shim.STANDARD_TOOL_PASS_PART_PREFIX 
            self.deformation_step_name = shim.STANDARD_DEFORMATION_STEP_NAME

        elif model_type is ModelTypes.TRACTION_APP:
            # It is assumed that a model which includes traction application
            #     will be created by copying another model which contains a simple 
            #     part geometry. In the copy, the part name will be unchanged.
            
            # Assumed to already exist at call time.
            # Assumed to result from copy operation.
            self.deformed_part_name = shim.STANDARD_INIT_GEOM_PART_NAME

            # Assumed to be chosen.
            self.new_model_name = shim.STANDARD_MODEL_NAME_PREFIX + str(model_cnt + 1)
            self.traction_step_name = shim.STANDARD_TRACTION_STEP_NAME
            self.traction_job_name = shim.STANDARD_JOB_PREFIX + str(model_cnt + 1)

        elif model_type is ModelTypes.EVAL_VOL_DIFF:
            # In order to compute the value of the volume difference function, 
            #     it's assumed that a model containing a deformed geometry will need
            #     to be copied and a part which contains the result of applying
            #     a traction will need to be created.
            
            # Assumed to already exist at call time.
            self.post_cut_post_deform_model_name = shim.STANDARD_MODEL_NAME
            self.post_cut_post_deform_part_name = shim.STANDARD_INIT_GEOM_PART_NAME
            
            # Assumed to be chosen.
            self.new_model_name = shim.STANDARD_DEFORMED_MODEL_PREFIX + str(model_cnt + 1)
            self.post_traction_part_name = shim.STANDARD_POST_TRACTION_PART_NAME



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

