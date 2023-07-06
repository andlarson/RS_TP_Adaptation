
"""
Fulfills the need for metadata assocaited with Abaqus simualtions.
"""


"""
Provides a repository of information which supplements the state saved in
   Abaqus.
Each object of this type is associated with a single MDB.
This external metadata is necessary because Abaqus' functionality doesn't
   suit our needs.
For example, when adding a step to a model in an MDB, it's necessary to
   specify the name of the step which should precede the new step. However,
   Abaqus only maintains a dictionary of the steps in each model. There is
   no notion of order in a dictionary, so it's necessary to maintain the 
   chronological ordering externally.

After a metadata object is created, it is updated via the same functions
   which update the underlying MDB (i.e. the Abaqus shim functions). This
   ensures that we don't need to remember to keep the metadata in sync with
   the underlying MDB.

Also, this metadata only includes additional information which is necessary.
   Its only purpose is to supplement the state that Abaqus maintains, not
   replace or duplicate it. For example, it's unnecessary to track the parts
   associated with a particular model since we don't care about ordering,
   etc.
"""
class MDBMetadata:
    
    def __init__(self, path_to_mdb, working_model_name):
    # type: (str, str) -> None
       
        self.path_to_mdb = path_to_mdb
        self.working_model_name = working_model_name

        # The record maps model names to data structures.
        self.record = {}


    def add_model(self, model_name):
    # type: (str) -> None

        # We don't care about the ordering of models.
        self.record[model_name] = []


    def set_working_model(self, model_name):
    # type: (str) -> None

        assert(model_name in self.record)
        self.working_model_name = model_name


    def get_working_model_name(self):
    # type: (None) -> None

        return self.working_model_name

    
    def add_step_to_model(self, model_name, step_name):
    # type: (str, str) -> None

        # Find the model, then add the step to it.
        # We maintain a chronological ordering here.
        self.record[model_name].append(step_name)
        

    def get_steps(self, model_name):
    # type: (str) -> None

        return self.record[model_name]

    
    def get_last_step(self, model_name):
    # type: (str) -> None

        return self.record[model_name][-1]



