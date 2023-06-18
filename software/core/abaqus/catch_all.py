
'''
For now, this is a catch-all file which interfaces with the Abaqus scripting
  interface.
This file will probably need to be broken up.
'''

from abaqus import *
from abaqusConstants import * 
import sketch

import os



STANDARD_MODEL_NAME = "Model-1"

# Create a MDB and open it. 
# This function does not automatically save the MDB. However, the path used
#   here is the location of save when a save does happen.
def create_mdb(name, path):
# type: (str, str) -> Any

    if not os.access(path, os.F_OK):
        raise RuntimeError("Path for MDB creation does not exist.") 

    if not name.endswith(".cae"):
        print("The MDB was suffixed with .cae implicitly.")

    return Mdb(path + name)



# Get a MDB object from the path to a .cae file.
def use_mdb(path_to_mdb):
# type: (str) -> Any 
    
    return openMdb(path_to_mdb)



# Build a part in a MDB. 
# TODO: For now, we assume the part is a rectangular prism. 
def build_part(name, rect_prism, mdb):
# type: (str, RectPrism, Any) -> None

    # Construct a sketch and extrude it to create the correct geometry.
    # Then translate the part to the correct location in space.
    
    v1, v2 = rect_prism.get_rect_corners() 
    mdb.models[STANDARD_MODEL_NAME].


# Used to verify that an MDB is in an expected state.
def verify_mdb_content(mdb):    
# type: (Any) -> None

    if len(mdb.models) != 1:
        raise RuntimeError("Too many or too few models in the MDB.")
    if not (STANDARD_MODEL_NAME in mdb.models):
        raise RuntimeError("The single model in the MDB is improperly named.")
    if len(mdb.models[STANDARD_MODEL_NAME].parts) != 1:
        raise RuntimeError("Too many of too few parts in the MDB.")
    if not ("Initial_Geometry" in mdb.models[STANDARD_MODEL_NAME].parts):
        raise RuntimeError("The single part in the single model in the MDB
                            is improperly named.")


