
'''
For now, this is a catch-all file which interfaces with the Abaqus scripting
  interface.
This file will probably need to be broken up.
'''

from abaqus import *
from abaqusConstants import * 
import sketch

import os


IDENTITY_SEQ_SEQ = [(1, 0, 0), (0, 1, 0), (0, 0, 1)]
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
# TODO: We assume a special right rectangular prism geometry. 
def build_part(name, spec_right_rect_prism, mdb):
# type: (str, SpecRightRectPrism, Any) -> None

    # The part may be floating in space away from the origin.
    # This necessitates re-centering the sketch origin.
    # We want to re-center in terms of x and y. We don't want to recenter in
    #   terms of z. Doing so would make extrusion impossible. Instead, the
    #   smaller z coordinate among the prism's vertices should be used. 
    centroid = spec_right_rect_prism.get_centroid()
    z_offset = spec_right_rect_prism.get_smaller_z()
    transform_matrix = centroid.append((centroid.x, centroid.y, z_offset)) 

    # Create a sketch object.
    v1, v2 = spec_right_rect_prism.get_rect_corners() 
    largest_coord = max(v1.x, v1.y, v2.x, v2.y)
    sketch = mdb.models[STANDARD_MODEL_NAME].ConstrainedSketch(name=name, sheetSize=(2 * largest_coord), transform=transform_matrix)

    # Draw the rectangle accounting for the translated origin.
    corner_1 = (v1.proj_xy().x - centroid.x, v1.proj_xy().y - centroid.y)
    corner_2 = (v2.proj_xy().x - centroid.x, v2.proj_xy().y - centroid.y)
    if not sketch.rectangle(corner_1, corner_2):
        raise RuntimeError("Failed to create rectangle on sketch!")

    # Build the part.
    # It takes the same name as the sketch.
    part = mdb.models[STANDARD_MODEL_NAME].Part(name=name, dimensionality=THREE_D, type=DEFORMABLE_BODY)

    # Extrude the sketch.
    depth = spec_right_rect_prism.get_dims()[2]
    mdb.models[STANDARD_MODEL_NAME].part[name].BaseSolidExtrude(sketch=sketch, depth=depth)



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


