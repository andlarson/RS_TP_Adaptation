
'''
This file is the only file which makes calls into the Abaqus API.
This file will probably need to be broken up.
'''

from abaqus import *
from abaqusConstants import * 

import Assembly
import Step

import os


# There are standard names/prefixes/suffixes which we expect in Abaqus MDBs.
STANDARD_MODEL_NAME = "Model-1"
STANDARD_INIT_GEOM_PART_NAME = "Initial_Geometry"
STANDARD_TOOL_PASS_PART_PREFIX = "Tool_Pass_"
STANDARD_POST_TOOL_PASS_PART_PREFIX = "Post_Tool_Pass_"
STANDARD_INITIAL_STEP_NAME = "Initial"
STANDARD_EQUIL_STEP_PREFIX = "Equilibrium"


# Create a MDB and open it. 
# This function does not automatically save the MDB. However, the path used
#   here is the location of save when a save does happen.
def create_mdb(name, path):
# type: (str, str) -> Any

    if not os.access(path, os.F_OK):
        raise RuntimeError("Path for MDB creation does not exist.") 

    return Mdb(path + name)



# Get a MDB object from the path to a .cae file.
# Treat opening an MDB like opening a file. Only open it in one place at one time
#    and close it as soon as it is no longer needed.
def use_mdb(path_to_mdb):
# type: (str) -> Any
    
    return openMdb(path_to_mdb)



# Close an MDB which is open.
def close_mdb(mdb):
# type: (Any) -> None

    mdb.close()




# Save a MDB object to the location which was specified during its creation.
def save_mdb(mdb):
# type: (Any) -> None

    save(mdb)
    print("The mdb was saved to: " + mdb.pathName)



# Build a part in a MDB. 
# TODO: We assume a special right rectangular prism geometry. 
# TODO: Break this up into modular chunks.
def build_part(name, spec_right_rect_prism, mdb):
# type: (str, SpecRightRectPrism, Any) -> Any

    # The part may be floating in space away from the origin.
    # This necessitates re-centering the sketch origin.
    # We want to re-center in terms of the x-y centroid and the z plane which
    #    is closer to z=0.
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
    # Note that BaseSolidExtrude() adds a feature to the Part and returns the
    #   associated feature object. We need not keep track of or deal with the
    #   feature object. By adding the feature, the part itself has been
    #   updated.
    depth = spec_right_rect_prism.get_dims()[2]
    mdb.models[STANDARD_MODEL_NAME].part[name].BaseSolidExtrude(sketch=sketch, depth=depth)

    return part



# Verify that the mdb the user passed as an initial part geometry really contains
#    what we expect it to.
def verify_init_geom_mdb(mdb):    
# type: (Any) -> None

    assert len(mdb.models) == 1
    assert STANDARD_MODEL_NAME in mdb.models
    assert len(mdb.models[STANDARD_MODEL_NAME].parts) == 1
    assert STANDARD_INIT_GEOM_PART_NAME in mdb.models[STANDARD_MODEL_NAME].parts
    assert len(mdb.jobs) == 0
    assert mdb.models[STANDARD_MODEL_NAME].rootAssembly == None  
    assert len(mdb.models[STANDARD_MODEL_NAME].predefinedFields) != 0
    assert len(mdb.models[STANDARD_MODEL_NAME].loads) == 0
    assert len(mdb.models[STANDARD_MODEL_NAME].materials) != 0
    assert len(mdb.models[STANDARD_MODEL_NAME].steps) == 0 

    

def get_part(name, mdb):
# type: (str, Any) -> Any    
    
    return mdb.models[STANDARD_MODEL_NAME].parts[name]



# Create part instance in the root assembly for the default model.
def instance_part_into_assembly(part, mdb, dependent):
# type: (Any, Any, bool) -> Any 
    
    return mdb.models[STANDARD_MODEL_NAME].rootAssembly.Instance(part, dependent=dependent)
    


# Cut a single part instance via multiple other part instances. 
# The result is a new part in the model and a new part instance in the root
#   assembly.
def cut_instances_in_assembly(name, instance_to_be_cut, cutting_instances, mdb):
# type: (str, Any, list[Any], Any) -> Any

    return mdb.models[STANDARD_MODEL_NAME].rootAssembly.PartFromBooleanCut(name, instance_to_be_cut, cutting_instances)



# Create an equilbirium step after another specified step.
def create_equilibrium_step(name, name_step_to_follow, mdb):
# type: (Any) -> None

    return mdb.models[STANDARD_MODEL_NAME].StaticStep(name, name_step_to_follow)



# Get the number of steps which exist in the standard single model.
def get_step_cnt(mdb):
# type: (Any) -> None

    return len(mdb.models[STANDARD_MODEL_NAME].steps)
        


def create_job():
    pass



def run_job():
    pass



# Enable visualization to some degree.
def enable_visualization():
    pass
