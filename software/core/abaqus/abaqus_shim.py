from abaqus import *
from abaqusConstants import * 
import regionToolset             # Unclear why necessary.   
import part                      # Unclear why necessary.

import os

import abaqus_metadata as abq_md
from util.debug import *


"""
This is the only file which should make calls directly into the Abaqus API.
"""


# There are standard names/prefixes/suffixes which we expect and impose in Abaqus 
#    MDBs.
STANDARD_MODEL_NAME = "Model-1"
STANDARD_MODEL_NAME_PREFIX = "Model-"
STANDARD_INIT_GEOM_PART_NAME = "Initial_Geometry"
STANDARD_TOOL_PASS_PART_PREFIX = "Tool_Pass_"
STANDARD_POST_TOOL_PASS_PART_PREFIX = "Post_Tool_Pass_"
STANDARD_INITIAL_STEP_NAME = "Initial"
STANDARD_EQUIL_STEP_PREFIX = "Equilibrium"
STANDARD_SECTION_NAME = "Section-1"
STANDARD_ORPHAN_MESH_FEATURE_NAME = "Orphan mesh-1"



# Create a MDB and open it. 
# This function does not automatically save the MDB. However, the path used
#   here is the location of save when a save (not save as) does happen.
def create_mdb(name, path):
# type: (str, str) -> Any

    return Mdb(path + name)



# Assign the only section in the model to the whole part.
def assign_only_section_to_part(part, model_name, mdb):
# type: (Any) -> Any

    model = mdb.models[model_name]

    # There better be only a single section.
    assert(len(model.sections) == 1)

    # A little hacky...We could also keep track of section names.
    section_name = model.sections.keys()[0]

    full_set = part.Set(name="simple_set", cells=part.cells)
    
    return part.SectionAssignment(full_set, section_name)



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



def save_mdb_as(save_path, mdb):
# type: (str, Any) -> None

    mdb.saveAs(save_path)    
    dp("The mdb was saved to: " + mdb.pathName)



# Check if the mdb contains an initial geometry.
def check_init_geom(should_print, mdb):    
# type: (bool, Any) -> bool

    if len(mdb.models) != 1 or \
       STANDARD_MODEL_NAME not in mdb.models:
        if should_print:
            dp("check_init_geom failure 1")
        return False

    model = mdb.models[STANDARD_MODEL_NAME]

    if len(model.parts) != 1 or \
       STANDARD_INIT_GEOM_PART_NAME not in model.parts:
        if should_print:
            dp("check_init_geom failure 2")
        return False

    part = model.parts[STANDARD_INIT_GEOM_PART_NAME]

    if len(model.materials) != 1 or \
       len(model.sections) != 1 or \
       STANDARD_SECTION_NAME not in model.sections or \
       len(part.sectionAssignments) != 1:
        if should_print:
            dp("check_init_geom failure 3")
        return False
       
    if len(model.rootAssembly.instances) != 0 or \
       len(model.predefinedFields) != 0 or \
       len(model.steps) != 1 or \
       len(mdb.jobs) != 0 or \
       len(model.loads) != 0:
        if should_print:
            dp("check_init_geom failure 4")
        return False

    return True
    


"""
Checks if the model contains multiple steps.
"""
def check_multiple_steps(should_print, model_name, mdb):
# type: (bool, str, Any) -> bool

    if model_name not in mdb.models:
        if should_print:
            dp("check_multiple_steps failure 1")
        return False

    if len(mdb.models[model_name].steps) > 1:
        if should_print:
            dp("check_multiple_steps failure 2")
        return True

    return False



"""
Checks if the model contains only an orphan mesh.
This is the content we expect when a new model is created and the results of
   a previous simulation have been imported via an ODB.
"""
def check_orphan_mesh(should_print, part_name, model_name, mdb):
# type: (bool, str, str, Any) -> bool

    if model_name not in mdb.models or \
       len(mdb.models[model_name].parts) != 1:
        if should_print:
            dp("check_orphan_mesh failure 1")
        return False

    model = mdb.models[model_name]

    if part_name not in model.parts:
        if should_print:
            dp("check_orphan_mesh failure 2")
        return False

    part = mdb.models[model_name].parts[part_name]

    if len(part.features) != 1 or \
       STANDARD_ORPHAN_MESH_FEATURE_NAME not in part.features:
        if should_print:
            dp("check_orphan_mesh failure 3")
        return False

    if len(model.materials) != 1 or \
       len(model.sections) != 1 or \
       len(part.sectionAssignments) != 1 or \
       len(model.rootAssembly.instances) != 1 or \
       len(model.predefinedFields) != 0:
        if should_print:
            dp("check_orphan_mesh failure 4")
        return False

    if len(model.steps) != 1 or \
       len(model.loads) != 0:
        if should_print:
            dp("check_orphan_mesh failure 5")
        return False

    return True



def suppress_feature(name, part):
# type: (str, Any) -> None

    part.suppressFeatures((name, ))



# TODO: We assume a special right rectangular prism geometry. 
# TODO: Break this up into modular chunks.
def build_part(name, spec_right_rect_prism, model_name, abq_metadata, mdb):
# type: (str, SpecRightRectPrism, str, abq_md.ABQMetadata, Any) -> Any
   
    # ----- Part Creation -----

    # Create the part in the model. 
    part = mdb.models[model_name].Part(name=name, dimensionality=THREE_D, type=DEFORMABLE_BODY)

    # ----- Metadata Updates -----

    abq_metadata.add_part_to_model(model_name, name)

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



def get_part(part_name, model_name, mdb):
# type: (str, str, Any) -> Any    
    
    return mdb.models[model_name].parts[part_name]



def find_step_keyword(kwb):
# type: (Any) -> int
    
    for index, string in enumerate(kwb.sieBlocks):
        if "Step" in string:
            return index

    raise RuntimeError("Failed to find step keyword!")



# Modify the input file so that a stress subroutine is included.
# The stress subroutine must also be associated with the job. This is done
#    elsewhere.
def inp_add_stress_subroutine(model_name, mdb):
# type: (str, Any) -> None

    # Synchronize the kwb to the current state of the model.
    # This must happen even if there have been no previous modifications to the
    #    input file.
    kwb = mdb.models[model_name].keywordBlock
    kwb.synchVersions()

    # The initial condition keyword is placed right before the first step in
    #    the input file.
    idx_step = find_step_keyword(kwb)
    idx_before_step = idx_step - 1
    kwb.insert(idx_before_step, "*Initial Conditions, Type=Stress, User")



# Modify the input file so that the stress state of the part is initialized
#    and mapped from the result of some other simulation. 
def inp_map_stress(path_sim_file, model_name, mdb):
# type: () -> None

    """
    There is some subtlety here.
    It seems like there are two ways to import a stress field which was 
       generated from a previous simulation as the initial state of a new
       simulation.
    We can use the "FILE" parameter in conjunction with the "*Initial Conditions"
       keyword. This technique does not offer any control over how mapping is
       done, which concerns me. What happens if the mesh used in the previous
       simulation and the mesh used for the current simulation are wildly
       different? A .odb or a .sim file can be used with this technique.
    Alternatively, we can use the more general "*External Field" keyword in
       conjunction with the "*Initial Conditions" keyword. This technique offers
       control over how mapping is done. It is also capable of accounting for
       translations and rotations of the geometry. A .sim file must be used
       with this technique.
    """

    # Synchronize the kwb to the current state of the model.
    # This must happen even if there have been no previous modifications to the
    #    input file.
    kwb = mdb.models[model_name].keywordBlock
    kwb.synchVersions()
    
    # The initial condition keyword is placed right before the first step in
    #    the input file.
    idx_step = find_step_keyword(kwb)
    idx_before_step = idx_step - 1
    kwb.insert(idx_before_step, "*Initial Conditions, Type=Stress")
    idx_before_external_field = idx_before_step + 1

    """
    Read Keywords > E > External Field to understand the content of the data
       line. 
    The source region is the whole previous model, the target region is the 
       whole new model, default mapping controls are used, and all stress 
       components are included from the tensors.
    """
    external_field_data_line = " , , , , ,S, , "
    kwb.insert(idx_before_external_field, "*External Field, File=" + path_sim_file + 
               "\n" + external_field_data_line)




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
def cut_instances_in_assembly(name, instance_to_be_cut, cutting_instances, model_name, abq_metadata, mdb):
# type: (str, Any, tuple[Any], str, Any, Any) -> Any

    abq_metadata.add_part_to_model(model_name, name)

    # Beware, the argument list ordering in the documentation for PartFrom
    #    BooleanCut() appears to be incorrect.
    return mdb.models[model_name].rootAssembly.PartFromBooleanCut(name=name, instanceToBeCut=instance_to_be_cut, cuttingInstances=cutting_instances)



def naive_mesh(part_instance, size, model_name, mdb):
# type: (Any, float, str, Any) -> None

    # Tetrahedrons seem to be able to mesh the widest variety of geometries. 
    mdb.models[model_name].rootAssembly.setMeshControls(part_instance.cells, elemShape=TET, technique=FREE)

    seq = (part_instance, )
    mdb.models[model_name].rootAssembly.seedPartInstance(seq, size)
    mdb.models[model_name].rootAssembly.generateMesh(regions=seq)

    while mdb.models[model_name].rootAssembly.getUnmeshedRegions() != None:
        size = size - 1

        if size <= 0:
            raise RuntimeError("There is no global element size that allows meshing for this particular mesh technique.")

        dp("An attempt at meshing failed. Decreasing global element size to " + str(size) + " and giving it another go.") 
        mdb.models[model_name].rootAssembly.deleteSeeds(seq)
        mdb.models[model_name].rootAssembly.seedPartInstance(seq, size)
        mdb.models[model_name].rootAssembly.generateMesh(regions=seq)


# Create an equilbirium step after another specified step.
def create_equilibrium_step(name, name_step_to_follow, model_name, abq_metadata, mdb):
# type: (str, str, str, abq_md.ABQMetadata, Any) -> None

    step = mdb.models[model_name].StaticStep(name, name_step_to_follow)
    
    abq_metadata.add_step_to_model(model_name, name)

    return step 



# Note that the name of the job matches the name of the various files produced
#    by the job when it is run (the .odb, .dat, etc. files).
def create_job(job_name, model_name, abq_metadata, mdb):
# type: (str, str, Any, Any) -> None 

    # The "resultsFormat=BOTH" causes Abaqus to generate .odb and .sim files
    #    when the simulation runs. This functionality is currently undocumented
    #    in the Abaqus documentation.
    job = mdb.Job(job_name, model_name, resultsFormat=BOTH)
    
    abq_metadata.add_job_to_model(model_name, job_name)
    
    return job



def add_user_subroutine(job, path_to_subroutine):
# type: (Any, str) -> None

    assert(job.userSubroutine == "")
    job.setValues(userSubroutine=path_to_subroutine)



def run_job(job):
# type: (Any) -> None

    job.submit()
    job.waitForCompletion()

    if job.status == ABORTED:
        raise RuntimeError("The job with name " + job.name + " was aborted!")

    print_job_messages(job)



# Check the messages produced by an analysis job.
def print_job_messages(job):
# type: (Any) -> None

    if len(job.messages) != 0:
        dp("The job with name " + job.name + " ran and some messages were received!")

    for message in job.messages:
        dp("A message of type " + str(message.type) + " was received when the job ran!")
    


def create_model_from_odb(odb_path, model_name, abq_metadata, mdb):
# type: (str, str, abq_md.ABQMetadata, Any) -> None

    model = mdb.ModelFromOdbFile(model_name, odb_path)

    abq_metadata.add_model(model_name)  

    return model



# Returns all the unique element faces which exist in the part. Unique in this 
#    context means that, even if two elements are next to one another and share
#    a face, the shared face will only be listed once in the returned sequence
#    of MeshFace objects.
def get_unique_element_faces(orphan_mesh_part):
# type: (Any) -> Any  

    return orphan_mesh_part.elementFaces



# Get the elements associated with a particular MeshFace object.
# A MeshFace object corresponds to a face which lives in a mesh.
def get_mesh_face_elements(mesh_face):
# type: (Any) -> Any
    
    return mesh_face.getElements()



def build_region_with_face(face, part):
# type: (Any) -> Any

    # Find the face number that the face belongs to on the element.
    face_num = face.face
    assert(face_num in (FACE1, FACE2, FACE3, FACE4, FACE5, FACE6))

    # Get the element that the face belongs to.
    # There better be only a single element. The face should belong to a unique
    #    element.
    # The result is a tuple of MeshElement objects.
    elems = face.getElements()
    assert(len(elems) == 1)

    # The tuple of MeshElement objects cannot be used directly to create a Region 
    #    object. This is because Region creation requires a MeshSequence, not a 
    #    tuple of MeshElement objects.
    # From the examples at Scripting Reference > Python Commands > Region Commands
    #    > Region Object, it appears that the standard way to obtain a
    #    MeshSequence is to do something like elements[3:5] where elements is
    #    the elements repository associated with a particular part. In our case,
    #    we can't do that directly because we only have the MeshElement object,
    #    not the index into the elements repository that the MeshElement object
    #    lives at.
    # This necessitates a more circuitious solution. We get the label of the
    #    MeshElement and then use the label to create a MeshElementArray, which
    #    I believe is a type of MeshSequence (i.e. MeshElementArray is derived
    #    from MeshSequence).

    # Get the element label for the single element.
    # This label is distinct from the index into the elements repository.
    elem_label = elems[0].label

    # Use the element label to get a MeshElementArray, which I believe is a
    #    type of MeshSequence.
    seq = part.elements.sequenceFromLabels((elem_label,))
        
    if face_num == FACE1:
        region = regionToolset.Region(face1Elements=seq)
    elif face_num == FACE2:
        region = regionToolset.Region(face2Elements=seq)
    elif face_num == FACE3:
        region = regionToolset.Region(face3Elements=seq)
    elif face_num == FACE4:
        region = regionToolset.Region(face4Elements=seq)
    elif face_num == FACE5:
        region = regionToolset.Region(face5Elements=seq)
    elif face_num == FACE6:
        region = regionToolset.Region(face6Elements=seq)
    else:
        raise RuntimeError("Failed to build region from face!")

    return region 



# Add a face feature to a part based on a region. 
def add_face_from_region(region, part):
# type: (Any, Any) -> Any
   
    # According to Scripting Reference > Python Commands > Region Commands >
    #    Region Object, whenever a command accepts a named set or surface, it
    #    will also accept a Region object. However, the converse does not seem
    #    to be true. Therefore a Region object really is required.
    return part.FaceFromElementFaces(region)



# Use a sequence of face features to add a solid feature to the part.
# Note that this may fail, and the Abaqus documentation doesn't give any hint
#    to the circumstances under which it may fail. Presumably it will fail if
#    if the face features don't obviously outline a solid.
def add_solid_from_faces(part):
# type: (Any) -> Any 

    feature = part.AddCells(part.faces)
    assert(feature != None)

    return feature 



# This adds a virtual topology feature to the part. All default options are
#    used. The goal of this is to simplify unimportant, small, and redundant 
#    geometric features, resulting in a part geometry which is easier to mesh.
def add_virtual_topology(part):

    # The documentation says that createVirtualTopology() will not issue exceptions.
    # I found that if createVirtualTopology() does not eliminate any features,
    #    then it issues an exception! It's unclear if this is the only exception
    #    which can be generated by createVirtualTopology().
    try:
        feature = part.createVirtualTopology(ignoreRedundantEntities=True)
    except:
       dp("No excess features were removed by creating a virtual topology!") 


