
'''
For now, this is a catch-all file which interfaces with the Abaqus scripting
  interface.
This file will probably need to be broken up.
'''

from abaqus import *
from abaqusConstants import * 

import os



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
def build_part(name, rect_prism, MDB):
# type: (str, RectPrism, Any) -> None

    
    

    


