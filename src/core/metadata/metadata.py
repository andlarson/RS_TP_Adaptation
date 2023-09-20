import os

import core.metadata.abaqus_metadata as abq_md


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


