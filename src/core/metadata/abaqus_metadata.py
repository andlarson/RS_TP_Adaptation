"""
This module contains code which can be used to keep track of and supplement the
    state saved in Abaqus.
"""

from typing import Any
import pathlib
import copy

import src.core.abaqus.abaqus_shim as shim


"""
Provides information which supplements the state saved in Abaqus.
Each object of this type is associated with a single MDB.
This external metadata is necessary because Abaqus' functionality doesn't
    suit our needs.

For example, when adding a step to a model in an MDB, the Abaqus API requires 
    the name of the step which should precede the new step. However,
    Abaqus only maintains a dictionary of the steps in each model. There is
    no notion of order in a dictionary, so it's necessary to maintain the 
    chronological ordering externally.

We also care about the order of models in an MDB. For example, we might want
    to simulate 5 consecutive tool passes. Each tool pass will require its
    own model, and each model is based on the results of the previous model.
    Therefore we require some notion of model ordering.

After a metadata object is created, it is updated via the same functions
    which update the underlying MDB (i.e. the Abaqus shim functions). This
    ensures that we don't need to remember to keep the metadata in sync with
    the underlying MDB.

Also, this metadata only includes additional information which is necessary.
    Its only purpose is to supplement the state that Abaqus maintains, not
    replace or duplicate it. For example, it's unnecessary to track the sketches 
    associated with a particular model since we don't care about their ordering,
    etc.
"""
class AbaqusMdbMetadata:

    def __init__(self, path: str) -> None:
        """Creates a data structure to record state information about a single MDB.

           Sets up metadata under assumption that mdb is in its default initial
               state. 

           Args:
               path: Absolute path to the .cae file. 

           Returns:
               None.

           Raises:
               None.
        """
    
        self.path_to_mdb = path

        # Order-matters list of models. Order mirrors order in the MDB.
        self.model_names = [shim.STANDARD_MODEL_NAME]

        self.models_metadata = {shim.STANDARD_MODEL_NAME: AbaqusModelMetadata()}
   
    
    def mdb_dir(self) -> str:
        """Returns absolute path to the directory that the MDB is in."""
        path = pathlib.Path(self.path_to_mdb)
        return str(path.parent)
    

    def add_model(self, name: str) -> None:
        """Adds a model to the data structure.

           Adding a model is special because it maps to a per-model data structure.

           Args:
               name: The name of the model.

           Returns:
               None.

           Raises:
               None.
        """

        self.model_names.append(name)
        self.models_metadata[name] = AbaqusModelMetadata()


    def add_copy(self, to_copy: str, new_name: str) -> None:
        """Adds a model to the data structure and populates its metadata with
               data from another model.
            
           Args:
               to_copy:  The name of the model to copy.
               new_name: The name of the new model.
        
           Returns:
               None.
        
           Raises:
               None.
        """

        self.model_names.append(new_name)
        self.models_metadata[new_name] = copy.deepcopy(self.models_metadata[to_copy])


    def delete_model(self, model_name: str) -> None:
        """Deletes a model from the data structure."""
    
        for i, name in enumerate(self.model_names):
            if name == model_name:
                del self.model_names[i] 
                return

        raise AssertionError("The model to delete doesn't exist!")


    def rename_model(self, old_model_name: str, new_model_name: str) -> None:
        """Renames a model in the data structure."""

        for i, name in enumerate(self.model_names):
            if name == old_model_name:
                self.model_names[i] = new_model_name
                return

        raise AssertionError("The model to rename doesn't exist!")



""" Provides per-model metadata to supplement the state saved in Abaqus."""
class AbaqusModelMetadata:

    def __init__(self) -> None:
        """Creates a data structure to record state information about a single model.

           Assumes that the model is in its default, unmodified state.

           Args:
               None.

           Returns:
               None.

           Raises:
               None.
        """
        
        # For all of the following lists, order matters. In particular, the
        #    ordering of these lists mirror the ordering in the Abaqus Model.
        self.step_names = [shim.STANDARD_INITIAL_STEP_NAME]
        self.part_names = [shim.STANDARD_INIT_GEOM_PART_NAME]
        self.job_name: str | None = None  # Assume 1-to-1 mapping of model to job.

    
    def add_part(self, name: str) -> None:
        """Adds a part to the data structure."""

        self.part_names.append(name)


    def rename_part(self, new_name: str, old_name: str) -> None:
        """Renames a part in the data structure."""
        
        for i, name in enumerate(self.part_names):
            if name == old_name:
                self.part_names[i] = new_name
                return

        raise AssertionError("Couldn't find matching part name.")


    def delete_part(self, name: str) -> None:
        """Deletes a part from the data structure."""

        for i, cur_name in enumerate(self.part_names):
            if name == cur_name:
                del self.part_names[i]
                return

        raise AssertionError("Couldn't find part to delete.")


