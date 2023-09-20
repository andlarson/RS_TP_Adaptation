from abaqus import *
from abaqusConstants import * 
import regionToolset 
import numpy as np

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


def get_step_keyword(kwb):
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
# type: (List[geom.Point3D], Any) -> Any

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
# type: (Any) -> List[List[geom.Point3D]]

    faces = get_faces(obj) 

    vertices = []
    for face in faces:
        vertices_on_face = []
        vertex_ids = face.getVertices()

        for vertex_id in vertex_ids:
            vertex = geom.Point3D(np.array(obj.vertices[vertex_id].pointOn[0]))
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

        ordered_vertices_on_face = get_face_ordered_vertices(vertices_on_face, face, obj) 

        for vertex_group in ordered_vertices_on_face:
            vertices[-1].append([geom.Point3D(np.array(vertex.pointOn[0])) for vertex in vertex_group])

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
def get_face_ordered_vertices(vertices, face, obj):
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
    first_neighbors = get_neighbor_vertices(first_vertex, face, obj)
    current_vertex = first_neighbors[0]
    
    well_ordered_vertices = [first_vertex]    
    while current_vertex.index != first_vertex.index:
        
        well_ordered_vertices.append(current_vertex)
    
        current_neighbors = get_neighbor_vertices(current_vertex, face, obj)
    
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
def get_neighbor_vertices(vertex, face, obj):
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
# type: (Any, Any) -> Int 

    return len(mdb.models[model_name].boundaryConditions)



# Check if the face and the ngon have matching vertices. 
# 
# Notes:
#    In general, I don't think that this means that they exactly match. The face
#       could be curved but have vertices which nicely match the ngon.
#
# Arguments:
#    ngon - NGon3D object.
#    face - Abaqus Face object.
#    obj  - Abaqus Part or Abaqus PartInstance object.
#
# Returns:
#    Boolean. 
def check_face_ngon_match(ngon, face, obj):

    face_vertices = get_face_vertices(face, obj)
    ngon_vertices = ngon.get_builtin_rep()
    ngon_face_share_vertices = True 
    if len(face_vertices) == len(ngon_vertices):
        for face_vertex in face_vertices:
            if face_vertex not in ngon_vertices:
                ngon_face_share_vertices = False
    else:
        ngon_face_share_vertices = False

    return ngon_face_share_vertices



# Try to find a face on an object with vertices which exactly match those of some
#    ngon. 
# 
# Notes:
#    If no matching face is found, this function returns None. 
#
# Arguments:
#    ngon - NGon3D object.
#    obj - Abaqus Part or Abaqus PartInstance object.
#
# Returns:
#    Abaqus Face object or None. 
def get_matching_face(ngon, obj):
# type: (geom.NGon3D, Any) -> Any
    
    new_face = None
    for face in get_faces(obj):

        vertices_single_face = get_face_vertices(face, obj)
       
        if len(vertices_single_face) == len(ngon.vertices):

            all_match = True
            for vertex in ngon.vertices:
                if vertex.components() not in vertices_single_face:
                    all_match = False

            if all_match:
                new_face = face 
                break
    
    return new_face




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
    


# Checks if the model contains multiple steps.
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



# Checks if the model contains only an orphan mesh.
# This is the content we expect when a new model is created and the results of
#    a previous simulation have been imported via an ODB.
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



# Extract the necessary information from the an Abaqus Transform object to
#    orient the coordinate system of the sketch in the global coordinate system.
#
# Notes:
#    The .matrix() method of the Abaqus Transform object is undocumented. 
#
# Arguments:
#    transform - Abaqus Transform object. 
#
# Returns:
#    CSys3D object.  
def extract_global_csys_to_sketch_csys(transform):
# type: (Any) -> geom.CSys3D 

    rot_and_trans = transform.matrix()
   
    x_axis = geom.Vec3D(np.array(rot_and_trans[0:3]))
    y_axis = geom.Vec3D(np.array(rot_and_trans[3:6]))
    z_axis = geom.Vec3D(np.array(rot_and_trans[6:9]))
    trans = geom.Vec3D(np.array(rot_and_trans[9:12]))

    basis = geom.Basis3D(x_axis, y_axis, z_axis) 

    return geom.CSys3D(trans, basis)



# Build and orient a sketch somewhere in space.
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

    sketch = mdb.models[model_name].ConstrainedSketch(sketch_name, sheetSize=1, transform=transform)

    if sketch == None:
        raise RuntimeError("Failed to build sketch!!")

    return sketch



# Take the first step towards building a part which is a special right rectangular
#    prism: put and orient a sketch somewhere in space and draw on it. 
# 
# Notes:
#    None.
#
# Arguments:
#    part_name             - String.
#                            Name of the part. Also used to name the sketch which is created.
#    spec_right_rect_prism - SpecRightRectangularPrism.
#                            The geometry of the part to create. 
#    model_name            - String.
#    mdb                   - Abaqus MDB object.
#
# Returns:
#    Abaqus ConstrainedSketch object.
def sketch_spec_right_rect_prism(part_name, spec_right_rect_prism, model_name, mdb):
# type: (str, geom.SpecRightRectPrism, str, Any) -> Any

    # The part may be floating in space away from the origin.
    # This necessitates re-centering the sketch origin.
    # We want to re-center in terms of the x-y centroid and the z plane which
    #    is closer to z=0.
    centroid = spec_right_rect_prism.get_centroid()
    z_offset = spec_right_rect_prism.get_smaller_z()

    # The datum plane is parallel to the x-y plane and offset by z_offset in
    #    the z direction.
    dp1_id = mdb.models[model_name].parts[part_name].DatumPointByCoordinate((0, 0, z_offset)).id
    dp2_id = mdb.models[model_name].parts[part_name].DatumPointByCoordinate((1, 0, z_offset)).id
    dp3_id = mdb.models[model_name].parts[part_name].DatumPointByCoordinate((0, 1, z_offset)).id 
    dp1 = mdb.models[model_name].parts[part_name].datums[dp1_id]
    dp2 = mdb.models[model_name].parts[part_name].datums[dp2_id]
    dp3 = mdb.models[model_name].parts[part_name].datums[dp3_id]

    # Create the plane and retrieve the datum object associated with it.
    sketch_plane_id = mdb.models[model_name].parts[part_name].DatumPlaneByThreePoints(dp1, dp2, dp3).id
    sketch_plane = mdb.models[model_name].parts[part_name].datums[sketch_plane_id]

    # Create a transform object associated with the part.
    t = mdb.models[model_name].parts[part_name].MakeSketchTransform(sketchPlane=sketch_plane, 
                                                                    origin=(centroid.rep[0], centroid.rep[1], z_offset))

    # Create the sketch using the transform object.
    sketch = build_constrained_sketch(t, part_name, model_name, mdb)

    # Draw the rectangle accounting for the translated origin.
    v1, v2 = spec_right_rect_prism.get_rect_corners() 
    corner_1 = (v1.proj_xy().rep[0] - centroid.rep[0], v1.proj_xy().rep[1] - centroid.rep[1])
    corner_2 = (v2.proj_xy().rep[0] - centroid.rep[0], v2.proj_xy().rep[1] - centroid.rep[1])

    # Don't try to check the return value. This returns None even on success...
    sketch.rectangle(corner_1, corner_2)

    return sketch



def build_part(name, spec_right_rect_prism, model_name, record, mdb):
# type: (str, geom.SpecRightRectPrism, str, md.CommittedToolPassMetadata, Any) -> Any
   
    # Create the part in the model. 
    part = mdb.models[model_name].Part(name=name, dimensionality=THREE_D, type=DEFORMABLE_BODY)

    # Update metadata for bookkeeping. 
    record.abaqus_mdb_metadata.models_metadata[model_name].part_names.append(name)

    # Construct the sketch and draw on it.
    sketch = sketch_spec_right_rect_prism(name, spec_right_rect_prism, model_name, mdb)

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
    idx_step = get_step_keyword(kwb)
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
    idx_step = get_step_keyword(kwb)
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



# Create an equilibrium step after another specified step.
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
        part.createVirtualTopology(ignoreRedundantEntities=True)
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



def create_displacement_bc(BC_name, step_name, region, settings, model_name, mdb):
# type: (str, str, Any, BC.BCSettings, str, Any)  -> None

    fix_x = SET if settings.fix_x else UNSET
    fix_y = SET if settings.fix_y else UNSET
    fix_z = SET if settings.fix_z else UNSET
    prevent_x_rot = SET if settings.prevent_x_rot else UNSET
    prevent_y_rot = SET if settings.prevent_y_rot else UNSET
    prevent_z_rot = SET if settings.prevent_z_rot else UNSET

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



# Find the face of the object that an ngon lives on.
# 
# Notes:
#    If the ngon does not exactly live on the face of an object, this function
#       may not be able to detect that this is the case. In this case, this function
#       may return some arbitrary face! This function does its best to avoid
#       such a scenario by doing checks when possible.
#
# Arguments:
#    ngon - NGon3D object.
#    obj  - Abaqus Part object or Abaqus PartInstance object.
#
# Returns:
#    Abaqus Face object.
def find_face_ngon_lives_on(ngon, obj):
# type: (geom.NGon3D, Any) -> Any

    """
    DEPRECATED TECHNIQUE!!!
    # Critical to ensure that the vertices are well-ordered on a per-face basis.
    vertices = shim.get_all_vertices_ordered(obj)

    face_ngon_belongs_to = None

    # Find the face that the ngon lives on.
    for idx, vertices_single_face in enumerate(vertices):

        flattened_vertices_single_face = [vertex for vertex_group in vertices_single_face for vertex in vertex_group]

        # Subtle point: If the points which define the vertices of a face are
        #    not on a plane (i.e. there is a curved surface), then this will
        #    almost always fail. This is intended behavior. No functionality for
        #    partitioning a curved face is desired.
        # TODO: Account for situation where the ngon lives inside two nested
        #    faces. Right now, the first face that the ngon lives in is the one
        #    which is identified and used.
        # TODO: Use the getCurvature() method to check that faces with curvature
        #    are not considered.
        if geom.on_plane_of_ngon(flattened_vertices_single_face, ngon):
            if geom.points_in_ngon_3D(ngon.vertices, geom.NGon3D(vertices_single_face)):
                face_ngon_belongs_to = shim.get_faces(obj)[idx]
                break

    if face_ngon_belongs_to == None:
        raise RuntimeError("Could not identify which face the ngon belongs to. \
                            Can't continue to do partitioning.")
    """

    # Check that the vertices of the ngon do live on the plane DEFINED by some face.
    # Note that, in general, even if the vertices of the ngon live on a plane 
    #    DEFINED by some face, it does not necessarily mean that the ngon actually
    #    lives on that face. That face might have holes in it! That face could be
    #    curved!
    on_a_face = False
    for vertices_single_face in get_all_vertices(obj):
        if geom.on_plane_of_ngon(vertices_single_face, ngon):
            on_a_face = True
    
    if not on_a_face:
        raise RuntimeError("Trying to create a partition with an ngon that doesn't \
                           live on any plane defined by the object's faces!")

    # Assume that the face closest to points which make up the ngon is the face
    #    that the ngon lives on.
    face_ngon_belongs_to = get_closest_face(ngon.vertices, obj)

    return face_ngon_belongs_to



# Sketch the points which make up an ngon on a sketch which has some orientation
#    in space. 
# 
# Notes:
#    None. 
#
# Arguments:
#    ngon      - NGon3D object.
#    transform - Abaqus Transform object.
#    sketch    - Abaqus ConstrainedSketch object.
#
# Returns:
#    Abaqus ConstrainedSketch object. 
def sketch_ngon(ngon, transform, sketch):
# type: (geom.NGon3D, Any, Any) -> Any

    # Extract the coordinate system of the sketch. 
    csys = extract_global_csys_to_sketch_csys(transform)

    # Map the points into the coordinate system of the sketch. 
    points = csys.map_into_csys(ngon.vertices)

    # In the coordinate system of the sketch, it better be the case that the 
    #    component of each point normal to the face is zero. 
    for point in points:
        if not geom.float_equals(point.rep[2], 0):
            raise RuntimeError("After mapping the points into the coordinate" + 
                               " system of the sketch, a point had a nonzero" +
                               " component in the axis normal to the sketch!")

    # Now use the points to construct the ngon on the sketch.
    points_circular = points + points[0:1]
    for idx in range(len(points_circular) - 1):
        cur_point = points_circular[idx].proj_xy().components()
        next_point = points_circular[idx + 1].proj_xy().components()
        
        # The Line() method of an Abaqus Sketch object returns None even on
        #    success. Documentation is wrong!
        sketch.Line(cur_point, next_point)

    return sketch



# Partitions the face of an object.
# 
# Notes:
#    Uses distinct part and instance arguments to avoid type checking Abaqus
#       object types. It's important to know the type of what is being partitioned
#       because the PartitionFaceBySketch() method is valid for Abaqus Part and
#       Abaqus Assembly objects.
#    If the polygon defining the partition does not live on a face of the object,
#       this may not be detected in general and bad things can happen.
#    If the ngon exactly matches a face which already exists, no new face is created.
#       The existing face is returned.
#
# Arguments:
#    ngon          - NGon3D object.
#                    Lies in global coordinate system and should be on a face of the
#                       the object which will be partitioned.
#    new_face_name - String.
#    model_name    - String.
#    mdb           - Abaqus MDB object.
#
# Optional Arguments:
#    part          - Abaqus Part object.
#    instance      - Abaqus PartInstance object.
# 
# Optional Arguments Notes:
#    Either part or instance must be specified.
#
# Returns:
#    Abaqus Face object.
def partition_face(ngon, new_face_name, model_name, mdb, part=None, instance=None):
# type: (geom.NGon3D, str, str, Any, Optional[Any], Optional[Any]) -> Any

    # TODO: UGLY, just do type checking with Abaqus object types...
    if (part is None and instance is None) or (part is not None and instance is not None):
        raise RuntimeError("Bad arguments passed in!")
    obj = part if part else instance

    face_ngon_belongs_to = find_face_ngon_lives_on(ngon, obj)

    # If the face the ngon belongs to has vertices which exactly match the
    #    ngon, we don't want to try to partition the face and create a new face.
    #    This will cause partitioning to fail and Abaqus to give up. 
    # TODO: The ngon could have redundant vertices, but such a partitioning will
    #    actually succeed because new vertices will be added to the face.
    # TODO: Assuming that the vertices which make up the Abaqus Face object are
    #    not redundant.
    # If they do share vertices, no partitioning needs to be done at all.
    if check_face_ngon_match(ngon, face_ngon_belongs_to, obj):
        return face_ngon_belongs_to 

    # Construct and orient a sketch which lives on that face.
    edge = get_any_edge_on_face(face_ngon_belongs_to, obj)
    assembly = mdb.models[model_name].rootAssembly
    transform = build_sketch_transform_from_face(face_ngon_belongs_to, edge, assembly)
    sketch = build_constrained_sketch(transform, new_face_name, model_name, mdb)

    # Sketch the ngon.
    sketch = sketch_ngon(ngon, transform, sketch)

    # Do the partitioning by using the sketch.
    if part is not None: 
        partition_face_with_sketch(face_ngon_belongs_to, edge, sketch, part=part)
    else:
        partition_face_with_sketch(face_ngon_belongs_to, edge, sketch, assembly=assembly, instance=instance)

    # Find the face which was just created via partitioning.
    new_face = get_matching_face(ngon, obj)

    if new_face is None:
        raise RuntimeError("Failed to find the face which was just created!")

    return new_face