"""
This file is the only file which makes calls into the Abaqus API.
"""

from abaqus import *
from abaqusConstants import * 

import os

import metadata as md

# DEBUG
from util.debug import *



# There are standard names/prefixes/suffixes which we expect in Abaqus MDBs.
# These act as invariants.
# Whenever an MDB is created, these names are used.
# Whenever an MDB is passed to this program, we check that the existing names
#    match these.
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
def save_mdb(save_path, mdb):
# type: (str, Any) -> None

    mdb.saveAs(save_path)    
    dp("The mdb was saved to: " + mdb.pathName)



# Verify that the mdb the user passed as an initial part geometry really contains
#    what we expect it to.
def verify_init_geom_mdb(mdb):    
# type: (Any) -> None

    assert len(mdb.models) == 1
    assert STANDARD_MODEL_NAME in mdb.models
    assert len(mdb.models[STANDARD_MODEL_NAME].parts) == 1
    assert STANDARD_INIT_GEOM_PART_NAME in mdb.models[STANDARD_MODEL_NAME].parts
    assert len(mdb.jobs) == 0
    assert len(mdb.models[STANDARD_MODEL_NAME].loads) == 0
    assert len(mdb.models[STANDARD_MODEL_NAME].materials) != 0

    # No initial stress field expected.
    assert len(mdb.models[STANDARD_MODEL_NAME].rootAssembly.instances) == 0
    assert len(mdb.models[STANDARD_MODEL_NAME].predefinedFields) == 0

    # There is always an initial step.
    assert len(mdb.models[STANDARD_MODEL_NAME].steps) == 1




# TODO: We assume a special right rectangular prism geometry. 
# TODO: Break this up into modular chunks.
def build_part(name, spec_right_rect_prism, model_name, mdb):
# type: (str, SpecRightRectPrism, str, Any) -> Any
   
    # ----- Part Creation -----

    # Create the part in the model. 
    part = mdb.models[model_name].Part(name=name, dimensionality=THREE_D, type=DEFORMABLE_BODY)


    # ----- Sketch Creation -----

    # The part may be floating in space away from the origin.
    # This necessitates re-centering the sketch origin.
    # We want to re-center in terms of the x-y centroid and the z plane which
    #    is closer to z=0.
    centroid = spec_right_rect_prism.get_centroid()
    z_offset = spec_right_rect_prism.get_smaller_z()

    # The datum plane is parallel to the x-y plane and offset by z_offset in
    #    the z direction.
    # TODO: The objects returned are feature objects. To create the plane, it's
    #    necessary retrieve the datum objects themselves. The datum objects
    #    should probably be recorded as model/part metadata.
    dp1_id = mdb.models[model_name].parts[name].DatumPointByCoordinate((0, 0, z_offset)).id
    dp2_id = mdb.models[model_name].parts[name].DatumPointByCoordinate((1, 0, z_offset)).id
    dp3_id = mdb.models[model_name].parts[name].DatumPointByCoordinate((0, 1, z_offset)).id 
    dp1 = mdb.models[model_name].parts[name].datums[dp1_id]
    dp2 = mdb.models[model_name].parts[name].datums[dp2_id]
    dp3 = mdb.models[model_name].parts[name].datums[dp3_id]

    # Create the plane and retrieve the datum object associated with it.
    sketch_plane_id = mdb.models[model_name].parts[name].DatumPlaneByThreePoints(dp1, dp2, dp3).id
    sketch_plane = mdb.models[model_name].parts[name].datums[sketch_plane_id]

    # Create a transform object associated with the part.
    t = mdb.models[model_name].parts[name].MakeSketchTransform(sketchPlane=sketch_plane, 
                                                               origin=(centroid.x, centroid.y, z_offset))

    # Create the sketch using the transform object.
    v1, v2 = spec_right_rect_prism.get_rect_corners() 
    largest_coord = max(v1.x, v1.y, v2.x, v2.y)
    sketch = mdb.models[model_name].ConstrainedSketch(name=name, sheetSize=(2 * largest_coord), transform=t)
    if sketch == None:
        raise RuntimeError("No sketch could be created!")

    # Draw the rectangle accounting for the translated origin.
    corner_1 = (v1.proj_xy().x - centroid.x, v1.proj_xy().y - centroid.y)
    corner_2 = (v2.proj_xy().x - centroid.x, v2.proj_xy().y - centroid.y)

    # Don't try to check the return value. This returns None even on success...
    # I think the documentation is out-of-date...
    sketch.rectangle(corner_1, corner_2)

    # ----- Part Extrusion -----

    # Extrude the sketch.
    # Note that BaseSolidExtrude() adds a feature to the Part and returns the
    #   associated feature object. We need not keep track of or deal with the
    #   feature object. By adding the feature, the part itself has been
    #   updated.
    depth = spec_right_rect_prism.get_dims()[2]
    mdb.models[model_name].parts[name].BaseSolidExtrude(sketch=sketch, depth=depth)
    
    return part



def get_part(name, mdb):
# type: (str, Any) -> Any    
    
    return mdb.models[STANDARD_MODEL_NAME].parts[name]



def find_step_keyword(kwb):
# type: (Any) -> int
    
    for index, string in enumerate(kwb.sieBlocks):
        if "Step" in string:
            return index

     
# TODO:
# A lot of this is bespoke b/c we only modify the input file for adding stress
#    user subroutine right now.

# Add a keyword and any associated data to the current representation of the 
#    input file associated with a model.
def modify_inp(keyword, parameters, data_lines, model_name, mdb):
# type: (str, tuple[str], str, Any) -> None

    if (keyword != "*Initial Conditions") or \
       (parameters[0] != "Type=Stress") or \
       (parameters[1] != "User") or \
       (data_lines != ""):
        raise RuntimeError("No support for that combination of keyword and parameters and data lines!")

    # First we must synchronize the kwb to the current state of the model.
    # This must happen even if there have been no previous modifications to the
    #    input file.
    kwb = mdb.models[model_name].keywordBlock
    kwb.synchVersions()

    idx_step = find_step_keyword(kwb)
    idx_before_step = idx_step - 1
    kwb.insert(idx_before_step, keyword + ", " + parameters[0] + ", " + parameters[1])
    


def get_step_cnt(model_name, mdb):
# type: (str, Any) -> None

    return len(mdb.models[model_name].steps)
 


# Create part instance in the root assembly for the default model.
def instance_part_into_assembly(instance_name, part, dependent, model_name, mdb):
# type: (str, Any, bool, str, Any) -> Any 
    
    return mdb.models[model_name].rootAssembly.Instance(instance_name, part, dependent=dependent)
    


# Cut a single part instance via multiple other part instances. 
# The result is a new part in the model.
# Remember to instance the resulting part!
def cut_instances_in_assembly(name, instance_to_be_cut, cutting_instances, model_name, mdb):
# type: (str, Any, tuple[Any], str, Any) -> Any

    # Beware, the argument list ordering in the documentation for PartFrom
    #    BooleanCut() appears to be incorrect.
    return mdb.models[model_name].rootAssembly.PartFromBooleanCut(name=name, instanceToBeCut=instance_to_be_cut, cuttingInstances=cutting_instances)



# Try a dumb meshing technique. 
def naive_mesh(part_instance, size, model_name, mdb):
# type: (Any, float, str, Any) -> None

    seq = (part_instance, )
    mdb.models[model_name].rootAssembly.seedPartInstance(seq, size)
    mdb.models[model_name].rootAssembly.generateMesh(regions=seq)



# Create an equilbirium step after another specified step.
def create_equilibrium_step(name, name_step_to_follow, model_name, mdb_metadata, mdb):
# type: (str, str, str, md.MDBMetadata, Any) -> None

    # Step creation means that MDB metadata must be updated.
    mdb_metadata.add_step_to_model(model_name, name)

    return mdb.models[model_name].StaticStep(name, name_step_to_follow)



def create_job(job_name, model_name, mdb):
# type: (str, str, Any) -> None 

    return mdb.Job(job_name, model_name)



def add_user_subroutine(job, path_to_subroutine):
# type: (Any, str) -> None

    assert(job.userSubroutine == "")
    job.setValues(userSubroutine=path_to_subroutine)



def run_job(job):
# type: (Any) -> None

    job.submit()
    job.waitForCompletion()
    print_job_messages(job)



# Check the messages produced by an analysis job.
def print_job_messages(job):
# type: (Any) -> None

    for message in job.messages:
        dp("A message of type " + str(message.type) + " was received when the job ran!")
    


# TODO:
# Enable visualization to some degree.
def enable_visualization():
    pass
