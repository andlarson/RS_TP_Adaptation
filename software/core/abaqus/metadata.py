import abaqus_shim as shim

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
   replace or duplicate it. For example, it's unnecessary to track the parts
   associated with a particular model since we don't care about ordering,
   etc.

Things we care about:
    Model sequence in an MDB.
    Part sequence for each model in an MDB.
    Step sequence for each model in an MDB.
"""
class MDBMetadata:
    
    def __init__(self, path_to_mdb):
    # type: (str, str) -> None
       
        self.path_to_mdb = path_to_mdb

        # Keeps track of order in which models are added to the MDB.
        self.models = []

        # The record maps model names to data structures.
        self.record = {}
    

    def add_model(self, model_name):
    # type: (str) -> None

        self.models.append(model_name)

        # Initialize the data structure for this model.
        self.record[model_name] = {}
        self.record[model_name]["step_names"] = []
        self.record[model_name]["part_names"] = []

        # Every model comes with a default initial step.
        self.add_step_to_model(model_name, shim.STANDARD_INITIAL_STEP_NAME)


    def get_last_model_name(self):
    # type: (None) -> str

        return self.models[-1]


    def add_step_to_model(self, model_name, step_name):
    # type: (str, str) -> None

        self.record[model_name]["step_names"].append(step_name)
        

    def get_steps(self, model_name):
    # type: (str) -> None

        return self.record[model_name]["step_names"]

    
    def get_last_step(self, model_name):
    # type: (str) -> None

        return self.record[model_name]["step_names"][-1]


    def add_part_to_model(self, model_name, part_name):
    # type: (str, str) -> None

        self.record[model_name]["part_names"].append(part_name)


    def get_last_part(self, model_name):
    # type: (str) -> None

        return self.record[model_name]["part_names"][-1]


