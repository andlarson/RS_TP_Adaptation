import os

import numpy as np
from abaqus import *
from abaqusConstants import * 
import regionToolset             # Unclear why necessary.   

import util.geom as geom
from util.debug import *


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
STANDARD_BC_PREFIX = "Boundary_Condition_"




# *****************************************************************************
#                         Abaqus Information Retrieval
#       These functions retrieve information/state that Abaqus maintains.
# *****************************************************************************


def find_step_keyword(kwb):
# type: (Any) -> int
    
    for index, string in enumerate(kwb.sieBlocks):
        if "Step" in string:
            return index

    raise RuntimeError("Failed to find step keyword!")



# Find the face closest to the centroid of a group of points.
# 
# Notes:
#    This function is a workaround to writing/calling nasty geometric routines. In
#       particular, it was introduced to avoid writing a routine which checks if
#       a three dimensional polygon is entirely inside a three dimensional polygon
#       which may have holes in it.
#    The default search tolerance is used.
#    Throws an exception if no face can be found. 
# 
# Arguments:
#    points - List of Point3D objects.
#    obj   - An Abaqus Part object or Abaqus PartInstance object.
#
# Returns:
#    Abaqus Face object.
def get_closest_face(points, obj):
# type: (List[Point3D], Any) -> Any

    centroid = geom.find_centroid(points)
    face_and_point = obj.faces.getClosest([centroid.components()])

    face = face_and_point[0][0]

    return face
    


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



# Documentation for getCentroid() is slightly incorrect. A tuple of tuple of
#    floats is returned, instead of a tuple of floats.
def get_face_centroid(face):
# type: (Any) -> Tuple[float, float, float]

    return tuple(face.getCentroid()[0])



def get_face_normal(face):
# type: (Any) -> Tuple[float, float, float] 

    return tuple(face.getNormal())



# Get an edge associated with a face. 
# 
# Notes:
# The edge is chosen arbitrarily.
# 
# Arguments:
#    face - Abaqus Face object.
#    obj  - Abaqus Part object or an Abaqus PartInstance 
#           Should be the object to which the face belongs. No check if done 
#              to ensure that the face actually belongs to this object.
#
# Returns:
#    Abaqus Edge object.
def get_any_edge_on_face(face, obj):
# type: (Any, Any) -> Any
    
    edge_id = face.getEdges()[0]
    return obj.edges[edge_id]



def get_part(part_name, model_name, mdb):
# type: (str, str, Any) -> Any    
    
    return mdb.models[model_name].parts[part_name]



def get_assembly(model_name, mdb):
# type: (str, Any) -> Any
    
    return mdb.models[model_name].rootAssembly



# Get the faces associated with any object which has faces.
# 
# Notes:
#    None.
# 
# Arguments:
#    obj - An Abaqus Part object or an Abaqus PartInstance object. 
#
# Returns:
#    An Abaqus FaceArray object containing all the Abaqus Face object.
def get_faces(obj):
# type: (Any) -> Any

    return obj.faces


# Get all the vertices associated with an object.
# 
# Notes:
# The documentation for the pointOn member of the Vertex object is incorrect. It
#    is really a tuple of tuple of floats, not a simple tuple of floats.
# 
# Arguments:
#    obj - An Abaqus Part object or an Abaqus PartInstance object.
#
# Returns:
#    A list of lists of Point3D objects. 
def get_all_vertices(obj):
# type: (Any) -> List[List[Point3D]]

    faces = get_faces(obj) 

    vertices = []
    for face in faces:
        vertices_on_face = []
        vertex_ids = face.getVertices()

        for vertex_id in vertex_ids:
            vertex = geom.Point3D(*obj.vertices[vertex_id].pointOn[0])
            vertices_on_face.append(vertex)

        vertices.append(vertices_on_face)

    return vertices 



# Get all the vertices associated with an object. This function ensures that, on
#    a per-face basis, the vertices are well-ordered. That is, the ordering of
#    the vertices mirrors the connectivity of the vertices.
# 
# Notes:
# The documentation for the pointOn member of the Vertex object is incorrect. It
#    is really a tuple of tuple of floats, not a simple tuple of floats.
# Consider the square face defined by the vertices: (0, 0), (1, 0), (0, 1), (1, 1).
#    The vertices are not well-ordered when the list of vertices is [(0, 0), (1, 0),
#    (0, 1), (1, 1)] because if you draw a line from (0, 0) to (1, 0), a line from
#    (1, 0) to (0, 1), and a line from (0, 1) to (1, 1) you don't get a square.
# There can be multiple groups of vertices on a single face. This is relevant
#    when there are holes within a face. Each group of vertices lives in its own
#    list. The lists are ordered arbitrarily. It might be the case that a list of
#    vertices which outlines a hole in a face is the first list of vertices for
#    that face.
# 
# Arguments:
#    obj - Abaqus Part object or Abaqus PartInstance object.
#
# Returns:
#    A list of lists of lists of vertices. 
def get_all_vertices_ordered(obj):
# type: (Any) -> List[List[List[geom.Point3D]]]

    vertices = []

    faces = get_faces(obj) 
    for face in faces:
        vertices.append([])

        vertex_ids = face.getVertices()
        vertices_on_face = [obj.vertices[vertex_id] for vertex_id in vertex_ids]

        ordered_vertices_on_face = order_vertices(vertices_on_face, face, obj) 

        for vertex_group in ordered_vertices_on_face:
            vertices[-1].append([geom.Point3D(*vertex.pointOn[0]) for vertex in vertex_group])

    return vertices 



# Orders the vertices of a face to reflect the edge connectivity of the face.
#    In the resulting list of Abaqus Vertex objects, the first vertex is connected
#    to the second vertex on the face, the second vertex is connected to the
#    third on the face, etc.
# A face could contain holes in it. In this case, this function groups the vertices
#    into per-boundary groups. Then, all the points which are connected on the
#    face are in the same group. Within a particular group, the vertices are
#    ordered to reflect connectivity. 
#
# Notes:
#    No lists are circular. 
#    This function does not care whether the face is convex or non-convex.
#
# Arguments:
#    vertices - List of Abaqus Vertex objects.
#    face     - Abaqus Face object.
#               The face on which the vertices lives.
#    obj      - Abaqus Part or Abaqus PartInstance.
#               The vertices and face belong to this.
#
# Returns:
#    List of lists of Abaqus Vertex objects. The order of per-boundary lists is
#       arbitrary.
def order_vertices(vertices, face, obj):
# type: (List[Any], Any, Any) -> List[List[Any]]

    per_group_vertices = []
    vertices_used = [False] * len(vertices) 

    while False in vertices_used:

        for i, vertex in enumerate(vertices):
            
            if vertices_used[i] == False:
              
                # Starting from an unvisited vertex, find all the other connected 
                #    vertices on the face.
                vertex_group = traverse_connected_vertices_on_face(vertices[i], face, obj)
                per_group_vertices.append(vertex_group)

                # Mark all the visited vertices as visited.
                for j, vertex in enumerate(vertices):
                    if vertex in vertex_group:
                        vertices_used[j] = True

    return per_group_vertices 



# Traverses a group of vertices via connections between the vertices, starting
#    from any vertex on a face.
#
# Notes:
#    None.
# 
# Arguments:
#    vertex - Abaqus Vertex object.
#    face   - Abaqus Face object.
#             The face on which the vertex lives.
#    obj    - Abaqus Part or Abaqus PartInstance.
#
# Returns:
#    List of Abaqus Vertex objects.
def traverse_connected_vertices_on_face(vertex, face, obj):
# type: (Any, Any, Any) -> List[Any]

    first_vertex = vertex
    prev_vertex = first_vertex 
    first_neighbors = find_neighbor_vertices(first_vertex, face, obj)
    current_vertex = first_neighbors[0]
    
    well_ordered_vertices = [first_vertex]    
    while current_vertex.index != first_vertex.index:
        
        well_ordered_vertices.append(current_vertex)
    
        current_neighbors = find_neighbor_vertices(current_vertex, face, obj)
    
        if current_neighbors[0].index != prev_vertex.index:
            prev_vertex = current_vertex
            current_vertex = current_neighbors[0]
        else:
            prev_vertex = current_vertex
            current_vertex = current_neighbors[1]

    return well_ordered_vertices



# Find the other vertices that a vertex is connected to on a face.
#
# Notes:
#    None.
#
# Arguments:
#    vertex - Abaqus Vertex object.
#    face   - Abaqus Face object.
#             The face on which the vertex lives.
#    obj    - Abaqus Part or Abaqus PartInstance.
#             The vertex and face belong to this.
#
# Returns:
#    List of exactly two Abaqus Vertex objects.
def find_neighbor_vertices(vertex, face, obj):
# type: (Any, Any, Any) -> List[Any, Any]

    face_edge_ids = face.getEdges()
    face_edges = [obj.edges[edge_id] for edge_id in face_edge_ids] 

    # Find the edges on the face connected to the vertex of interest.
    connected_vertices = []
    for edge in face_edges:
        if vertex.index == edge.getVertices()[0]:
            connected_vertices.append(obj.vertices[edge.getVertices()[1]]) 
        elif vertex.index == edge.getVertices()[1]:
            connected_vertices.append(obj.vertices[edge.getVertices()[0]]) 

    assert(len(connected_vertices) == 2)

    return connected_vertices



# Error in documentation of pointOn member.
def get_edge_vertices(edge, part):
# type: (Any, Any) -> Tuple[Tuple[float, float, float], Tuple[float, float, float]]

    vertex_ids = edge.getVertices()
    vertices = [part.vertices[vertex_id].pointOn[0] for vertex_id in vertex_ids]
    return tuple(vertices)
   


# Error in documentation of pointOn member.
def get_face_vertices(face, part):
# type: (Any, Any) -> Tuple[Tuple[float, float, float], ...]

    vertex_ids = face.getVertices()
    vertices = [part.vertices[vertex_id].pointOn[0] for vertex_id in vertex_ids]
    return tuple(vertices)



def get_step_cnt(model_name, mdb):
# type: (str, Any) -> None

    return len(mdb.models[model_name].steps)



def get_BC_cnt(model_name, mdb):
# type: (Any) -> Int 

    return len(mdb.models[model_name].boundaryConditions)




# *****************************************************************************
#                         Abaqus Object Manipulation 
#    These functions manipulate the state/information that Abaqus maintains.
# *****************************************************************************


# Create a MDB and open it. 
# This function does not automatically save the MDB. However, the path used
#   here is the location of save when a save (not save as) does happen.
def create_mdb(name, path):
# type: (str, str) -> Any

    return Mdb(path + name)



# Assign the only section in the model to the whole part.
def assign_only_section_to_part(part, model_name, mdb):
# type: (Any, str, Any) -> Any

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



# Build a matrix which specifies the location and orientation of a sketch in 
#    space.
#
# Notes:
#    It appears that when the MakeSketchTransform() method is called and an
#       Abaqus Face object is passed in, the sketchUpEdge argument is required
#       even though the documentation says otherwise.
#    All other optional arguments for MakeSketchTransform() are left unspecified.
#       These arguments fully specify the spacial orientation of a sketch. In
#       general, it's not important to specify a particular orientation in space.
#       Rather, it's important to know what the orientation in space is.
#    An Abaqus Transform object has a single method called matrix() that returns
#       the matrix which contains the translation and orientation of the sketch
#       in space. Extracting this information can be useful to map points in the
#       global coordinate system into the coordinate system of this sketch.
#
# Arguments:
#    face - Abaqus Face object. 
#           Locates the sketch in space.
#    edge - Abaqus Edge object. 
#           Orients the sketch in space.
#    obj  - Abaqus Part object or Abaqus Assembly object. 
#           The object that the transform is associated with.
#
# Returns:
#    Abaqus Transform object.
def build_sketch_transform_from_face(face, edge, obj):
# type: (Any, Any, Any) -> Any

    return obj.MakeSketchTransform(face, sketchUpEdge=edge)



# Uses the undocumented matrix() method of an Abaqus Transform object to extract
#    the information necessary to orient the coordinate system of the sketch in
#    terms of the global coordinate system.
#
# Notes:
#    None.
#
# Arguments:
#    transform - Abaqus Transform object. 
#
# Returns:
#    CSys3D object.  
def extract_global_csys_to_sketch_csys(transform):
# type: (Any) -> geom.CSys3D 

    rot_and_trans = transform.matrix()
   
    x_axis = np.array(rot_and_trans[0:3])
    y_axis = np.array(rot_and_trans[3:6])
    z_axis = np.array(rot_and_trans[6:9])
    trans = np.array(rot_and_trans[9:12])

    basis = np.stack((x_axis, y_axis, z_axis), axis=0)

    return geom.CSys3D(trans, basis)



# Build and orient a sketch somewhere in the space.
# 
# Notes:
#    The sheetSize argument doesn't matter when the sketch is not being drawn
#       by hand via a UI.
#    Sketches are only associated with whole models.
#
# Arguments:
#    transform   - Abaqus Transform object. 
#                  Specifies the location and orientation of the sketch in space.
#    sketch_name - String. 
#                  Name of the new sketch.
#    model_name  - String.
#    mdb         - Abaqus MDB object.
#
# Returns:
#    Abaqus ConstrainedSketch object.
def build_constrained_sketch(transform, sketch_name, model_name, mdb):
# type: (Any, str, str, Any) -> None

    return mdb.models[model_name].ConstrainedSketch(sketch_name, sheetSize=1, transform=transform)



# TODO: We assume a special right rectangular prism geometry. 
# TODO: Break this up into modular chunks.
def build_part(name, spec_right_rect_prism, model_name, record, mdb):
# type: (str, SpecRightRectPrism, str, md.CommittedToolPassMetadata, Any) -> Any
   
    # ----- Part Creation -----

    # Create the part in the model. 
    part = mdb.models[model_name].Part(name=name, dimensionality=THREE_D, type=DEFORMABLE_BODY)

    # ----- Metadata Updates -----

    record.abaqus_mdb_metadata.models_metadata[model_name].part_names.append(name)

    # ----- Sketch Creation -----

    # The part may be floating in space away from the origin.
    # This necessitates re-centering the sketch origin.
    # We want to re-center in terms of the x-y centroid and the z plane which
    #    is closer to z=0.
    centroid = spec_right_rect_prism.get_centroid()
    z_offset = spec_right_rect_prism.get_smaller_z()

    # The datum plane is parallel to the x-y plane and offset by z_offset in
    #    the z direction.
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
    corner_1 = (v1.proj_xy().x1 - centroid.x, v1.proj_xy().x2 - centroid.y)
    corner_2 = (v2.proj_xy().x1 - centroid.x, v2.proj_xy().x2 - centroid.y)

    # Don't try to check the return value. This returns None even on success...
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
# type: (str, str, Any) -> None

    # There is some subtlety here.
    # It seems like there are two ways to import a stress field which was 
    #    generated from a previous simulation as the initial state of a new
    #    simulation.
    # We can use the "FILE" parameter in conjunction with the "*Initial Conditions"
    #    keyword. This technique does not offer any control over how mapping is
    #    done, which concerns me. What happens if the mesh used in the previous
    #    simulation and the mesh used for the current simulation are wildly
    #    different? A .odb or a .sim file can be used with this technique.
    # Alternatively, we can use the more general "*External Field" keyword in
    #    conjunction with the "*Initial Conditions" keyword. This technique offers
    #    control over how mapping is done. It is also capable of accounting for
    #    translations and rotations of the geometry. A .sim file must be used
    #    with this technique.

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

    # Read Keywords > E > External Field to understand the content of the data
    #    line. 
    # The source region is the whole previous model, the target region is the 
    #    whole new model, default mapping controls are used, and all stress 
    #    components are included from the tensors.
    external_field_data_line = " , , , , ,S, , "
    kwb.insert(idx_before_external_field, "*External Field, File=" + path_sim_file + 
               "\n" + external_field_data_line)



# Create part instance in the root assembly for the default model.
def instance_part_into_assembly(instance_name, part, dependent, model_name, mdb):
# type: (str, Any, bool, str, Any) -> Any 
    
    return mdb.models[model_name].rootAssembly.Instance(instance_name, part, dependent=dependent)
    


# Cut a single part instance via multiple other part instances. 
# The result is a new part in the model.
# Remember to instance the resulting part!
def cut_instances_in_assembly(name, instance_to_be_cut, cutting_instances, model_name, record, mdb):
# type: (str, Any, tuple[Any], str, md.CommittedToolPassMetadata, Any) -> Any

    record.abaqus_mdb_metadata.models_metadata[model_name].part_names.append(name)

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
def create_equilibrium_step(name, name_step_to_follow, model_name, record, mdb):
# type: (str, str, str, md.CommittedToolPassMetadata, Any) -> None

    step = mdb.models[model_name].StaticStep(name, name_step_to_follow)
    
    record.abaqus_mdb_metadata.models_metadata[model_name].step_names.append(name)

    return step 



# Note that the name of the job matches the name of the various files produced
#    by the job when it is run (the .odb, .dat, etc. files).
def create_job(job_name, model_name, record, mdb):
# type: (str, str, md.CommittedToolPassMetadata, Any) -> None 

    # The "resultsFormat=BOTH" causes Abaqus to generate .odb and .sim files
    #    when the simulation runs. This functionality is currently undocumented
    #    in the Abaqus documentation.
    job = mdb.Job(job_name, model_name, resultsFormat=BOTH)
    
    record.abaqus_mdb_metadata.models_metadata[model_name].job_name = job_name
    
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
    


def create_model_from_odb(odb_path, model_name, record, mdb):
# type: (str, str, md.CommittedToolPassMetadata, Any) -> None

    model = mdb.ModelFromOdbFile(model_name, odb_path)

    record.abaqus_mdb_metadata.add_model(model_name)

    return model



# Build a region on a root assembly using the reference points which already
#    exist on the root assembly.
# Note that this function expects a list of reference point features, which it
#    translates to reference point objects.
def build_region_with_ref_points(ref_point_features, model_name, mdb):
# type: (List[Any], str, Any) -> Any

    root_assembly = mdb.models[model_name].rootAssembly

    # Construct a sequence of ReferencePoint objects.
    ref_point_objects = []
    for ref_point_feature in ref_point_features:
        key = ref_point_feature.id
        ref_point_objects.append(root_assembly.referencePoints[key]) 

    return regionToolset.Region(referencePoints=ref_point_objects)



# Build a region out of a face.
#
# Notes:
#    None.
#
# Arguments:
#    face_feature - Abaqus Face object 
#    obj          - Abaqus Part object or Abaqus PartInstance object.
#                   The face which underlies the feature must be a face on this
#                      object.
#
# Returns:
#    Abaqus Region object.
def build_region_with_face(face, obj):
# type: (Any, Any) -> Any

    # Even though the Abaqus Region() constructor says that it takes "A sequence 
    #    of Face objects", it really takes an Abaqus FaceArray object.
    # If an Abaqus FaceArray object is not passed, an exception will be issued
    #    by the function saying that it expected an object of type "GeomSequence".
    #    I believe that FaceArray, EdgeArray, etc. are all subtypes of
    #    GeomSequence.
    # In order to go from a single Abaqus Face object to an Abaqus FaceArray
    #    object with only a single Abaqus Face object in it, the following
    #    sequence of steps is necessary.
    # Also, note that the findAt() method for a FaceArray object always returns 
    #    an Abaqus FaceArray object, even though the documentation says that 
    #    if it only finds a single Face object, then a single Face object is
    #    returned.
    point_on_face = face.pointOn
    face_seq = obj.faces.findAt(point_on_face)
    region = regionToolset.Region(faces=face_seq)

    return region
    


def build_region_with_elem_face(face, part):
# type: (Any, Any) -> Any

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
# type: (Any) -> None

    # The documentation says that createVirtualTopology() will not issue exceptions.
    # I found that if createVirtualTopology() does not eliminate any features,
    #    then it issues an exception! It's unclear if this is the only exception
    #    which can be generated by createVirtualTopology().
    try:
        feature = part.createVirtualTopology(ignoreRedundantEntities=True)
    except:
       dp("No excess features were removed by creating a virtual topology! This does not necessarily indicate that there is a bug.") 



# Add some reference points to a root assembly.
# The reference points can, for example, then be used to build a region.
# There is a bit of subtlety here. To build a region with reference points, a
#    sequence of ReferencePoint objects is needed. This function returns a
#    sequence of Feature objects which are tied to ReferencePoint objects. To
#    translate from the Feature objects to ReferencePoint objects, it's
#    necessary to use the 'id' data member of the Feature object to get the
#    key for the referencePoint repository which is associated with the root
#    assembly.
def add_ref_points(points, model_name, mdb):
# type: (List[geom.Point3D], str, Any) -> List[Any]

    root_assembly = mdb.models[model_name].rootAssembly
    ref_point_features = []

    for point in points:
        ref_point_features.append(root_assembly.ReferencePoint((point.x, point.y, point.z)))

    return ref_point_features



def create_displacement_bc(BC_name, step_name, region, fix_u1, fix_u2, fix_u3, fix_ur1, fix_ur2, fix_ur3, model_name, mdb):
# type: (str, str, Any, bool, bool, bool, bool, bool, bool, str, Any) 

    fix_x = SET if fix_u1 else UNSET
    fix_y = SET if fix_u2 else UNSET
    fix_z = SET if fix_u3 else UNSET
    prevent_x_rot = SET if fix_ur1 else UNSET
    prevent_y_rot = SET if fix_ur2 else UNSET
    prevent_z_rot = SET if fix_ur3 else UNSET

    mdb.models[model_name].DisplacementBC(name=BC_name, createStepName=step_name, region=region, u1=fix_x, u2=fix_y, u3=fix_z, ur1=prevent_x_rot, ur2=prevent_y_rot, ur3=prevent_z_rot)



# Partitions the face of an object with a sketch.
# 
# Notes:
#    The PartitionFaceBySketch() method claims to take a "sequence of Face
#       objects". I believe it should instead take an Abaqus FaceArray object. 
# 
# Arguments:
#    face     - Abaqus Face object. 
#               The face to partition.
#    edge     - Abaqus Edge object. 
#               Specifies the orientation of the sketch on the
#                  face. This edge should have been used to orient the sketch on 
#                  the face when the sketch was constructed.
#    sketch   - Abaqus ConstrainedSketch object.
#
# Optional Arguments:
#    assembly - Abaqus Assembly object.
#    instance - Abaqus PartInstance object.
#    part     - Abaqus Part object.
#
# Optional Arguments Notes:
#    Either assembly and instance must be specified or part must be specified.
#
# Returns:
#    None. The Abaqus Feature object returned by the PartitionFaceBySketch()
#       method seems to be useless.
def partition_face_with_sketch(face, edge, sketch, assembly=None, instance=None, part=None):
# type: (Any, Any, Any, Optional[Any], Optional[Any], Optional[Any]) -> Any

    point = face.pointOn

    if assembly != None and instance != None:
        face_array = instance.faces.findAt(point)
        assembly.PartitionFaceBySketch(face_array, sketch, sketchUpEdge=edge)
    elif part != None:
        face_array = part.faces.findAt(point)
        part.PartitionFaceBySketch(face_array, sketch, sketchUpEdge=edge)
    else:
        raise RuntimeError("Invalid combination of arguments passed.")


