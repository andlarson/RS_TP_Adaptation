"""
This module contains code which can be used to keep track of and supplement the
   state saved in Abaqus.
"""

import core.abaqus.abaqus_shim as shim



# Provides information which supplements the state saved in Abaqus.
# Each object of this type is associated with a single MDB.
# This external metadata is necessary because Abaqus' functionality doesn't
#    suit our needs.
# 
# For example, when adding a step to a model in an MDB, the Abaqus API requires 
#    the name of the step which should precede the new step. However,
#    Abaqus only maintains a dictionary of the steps in each model. There is
#    no notion of order in a dictionary, so it's necessary to maintain the 
#    chronological ordering externally.
# 
# We also care about the order of models in an MDB. For example, we might want
#    to simulate 5 consecutive tool passes. Each tool pass will require its
#    own model, and each model is based on the results of the previous model.
#    Therefore we require some notion of model ordering.
# 
# After a metadata object is created, it is updated via the same functions
#    which update the underlying MDB (i.e. the Abaqus shim functions). This
#    ensures that we don't need to remember to keep the metadata in sync with
#    the underlying MDB.
# 
# Also, this metadata only includes additional information which is necessary.
#    Its only purpose is to supplement the state that Abaqus maintains, not
#    replace or duplicate it. For example, it's unnecessary to track the parts
#    associated with a particular model since we don't care about ordering,
#    etc.
class AbaqusMdbMetadata:

    # Construct data structure based on a default MDB.
    #
    # Notes:
    #    Sets up metadata under assumption that mdb is in its default initial
    #       state. 
    #
    # Arguments:
    #    path - String. 
    #           Path to the MDB. 
    # 
    # Returns:
    #    None.
    def __init__(self, path):
    # type: (str) -> None
    
        self.path_to_mdb = path

        # Order-matters list of models. Order mirrors order in the MDB.
        self.model_names = [shim.STANDARD_MODEL_NAME]

        # Each model has a metadata data structure associated with it. 
        self.models_metadata = {shim.STANDARD_MODEL_NAME: AbaqusModelMetadata()}

    
    # Add a new model to the MDB record. This is special because a new model
    #    maps to a new simulation. 
    #
    # Notes:
    #    None.
    # 
    # Arguments:
    #    name - String. Name of the new model.
    #
    # Returns:
    #    None.
    def add_model(self, name):
    # type: (str) -> None

        self.model_names.append(name)
        self.models_metadata[name] = AbaqusModelMetadata()



# Provides per-model metadata to supplement the state saved in Abaqus.
class AbaqusModelMetadata:

    # Construct data structure based on a default model.
    #
    # Notes:
    #    Assumes the model is in its default, unmodified state. 
    #
    # Arguments:
    #    None. 
    # 
    # Returns:
    #    None.
    def __init__(self):
    # type: (None) -> None
        
        # For all of the following lists, order matters. In particular, the
        #    ordering of these lists mirror the ordering in the Abaqus Model.
        self.step_names = [shim.STANDARD_INITIAL_STEP_NAME]
        self.part_names = [shim.STANDARD_INIT_GEOM_PART_NAME]
        self.job_name = None      # Assume 1-to-1 mapping of model to job.