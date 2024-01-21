"""
This module offers functions which, generally, are wrappers for Abaqus Python
    API calls. It is desirable to wrap the Abaqus Python API calls for a
    couple reasons:
       1) Some of the Abaqus Python API documentation is wrong or misleading.
             It's therefore helpful to centralize clear and correct behavior
             in a single wrapper function.
       2) There are common sequences of code which surround some Abaqus API
              calls. Encapsulating these common sequences reduces code
              duplication.
       3) From an architectural point of view, it's helpful to separate the
              functionality that Abaqus provides from other functionality.

While all of the above is true, it's completely fine to use the Abaqus API 
    outside of this file. For example, it's often necessary to get an Abaqus 
    Model object by using an Abaqus MDB object, etc. The idea is to encapsulate
    complex and commonly used accesses to Abaqus' API into this file.
"""

from typing import Any, Optional

from abaqus import *
from abaqusConstants import * 
import regionToolset 
import part
import odbAccess
import numpy as np

import src.core.material_properties.material_properties as mp
import src.core.tool_pass.tool_pass as tp
import src.util.geom as geom

from src.util.debug import *


# Built-in names in Abaqus.
STANDARD_MODEL_NAME = "Model-1"
STANDARD_INIT_GEOM_PART_NAME = "Initial_Geometry"
STANDARD_INITIAL_STEP_NAME = "Initial"
STANDARD_SECTION_NAME = "Section-1"
STANDARD_MATERIAL_NAME = "Material-1"
STANDARD_JOB_PREFIX = "Job-"
STANDARD_ORPHAN_MESH_FEATURE_NAME = "Orphan mesh-1"

# Custom names.
STANDARD_MODEL_NAME_PREFIX = "Model-"
STANDARD_TOOL_PASS_PART_PREFIX = "Tool_Pass_"
STANDARD_POST_TOOL_PASS_PART_PREFIX = "Post_Tool_Pass_"
STANDARD_EQUIL_STEP_PREFIX = "Equilibrium"
STANDARD_BC_PREFIX = "Boundary_Condition_"
STANDARD_BOUNDING_BOX_PART_NAME = "Bounding_Box"
STANDARD_EXCESS_BOUNDING_BOX_PART_NAME = "Bounding_Box_Excess"
STANDARD_BOUNDING_BOX_INIT_GEOM_PART_NAME = "Bounding_Box_And_Initial_Geometry"
STANDARD_INIT_GEOM_WITH_BBOX_NAME = "Initial_Geometry_With_Bounding_Box"


# *****************************************************************************
#                         Abaqus Information Retrieval
#       These functions retrieve information/state that Abaqus maintains.
# *****************************************************************************


def get_model_cnt(mdb: Any) -> int:
    """Gets the number of models in an Abaqus MDB.
    
       Args:
           mdb: Abaqus MDB object.

       Returns:
           Number of models in the MDB.

       Raises:
           None.
    """

    return len(mdb.models)



def get_only_part_name(model_name: str, mdb: Any) -> str:
    """In a model with only a single part, gets the name of the part.

       Args:
           model_name: The name of the model in the Abaqus MDB.
           mdb:        Abaqus MDB object. MDB in which the model lives.

       Returns:
           Name of the part.

       Raises:
           None. 
    """
    
    if len(mdb.models[model_name].parts.keys() != 1):
        raise ValueError("There should only be a single part in the model.")

    return mdb.models[model_name].parts.keys()[0]



def get_step_keyword(kwb: Any) -> int:
    """Finds the first occurence of the "Step" keyword in an input file. 

       Args:
           kwb: Abaqus KeywordBlock object.

       Returns:
           The index of the "Step" keyword in the Abaqus KeywordBlock object. 

       Raises:
           RuntimeError: The "Step" keyword was not found. 
    """
    
    for index, string in enumerate(kwb.sieBlocks):
        if "Step" in string:
            return index

    raise RuntimeError("Failed to find Step keyword!")



def get_closest_face(points: list[geom.Point3D], obj: Any) -> Any:
    """Gets the face closest to the centroid of some points.

       This function is a workaround to writing/calling nasty geometric routines. 
           In particular, it was introduced to avoid writing a routine which 
           checks if a three dimensional polygon is entirely inside a three 
           dimensional polygon which may have holes in it.
       The default search tolerance is used.
       Note that the Abaqus routines called in this function may throw an
           exception if no face is found.

       Args:
           points: A nonzero number of points. 
           obj:    Abaqus Part object or Abaqus PartInstance object.

       Returns:
           Abaqus Face object.

       Raises:
           None.
    """

    if len(points) == 0:
        raise ValueError("There should be a nonzero number of points.")

    centroid = geom.find_centroid(points)
    face_and_point = obj.faces.getClosest([centroid.components()])

    face = face_and_point[0][0]

    return face
    


def get_unique_element_faces(part: Any) -> Any:
    """Gets the unique element faces in a part.

       Unique in this context means that, even if two elements are next to 
           one another and share a face, the shared face will only be listed 
           once in the returned sequence. 

       If the part is neither meshed nor an orphan, there are no elements
           associated with it so this function is meaningless.

       Args:
           part: Abaqus Part object.

       Returns:
           Abaqus MeshFaceArray object.

       Raises:
           None.
    """

    return part.elementFaces



def get_mesh_face_elements(mesh_face: Any) -> tuple[Any, ...]:
    """Gets the elements associated with a single mesh face.

       Args:
           mesh_face: Abaqus MeshFace object.

       Returns:
           Tuple of Abaqus MeshElement objects. 

       Raises:
           None.
    """
    
    return mesh_face.getElements()



def get_any_edge_on_face(face: Any, obj: Any) -> Any:
    """Gets an arbitrary edge on a face.

       Args:
           face: Abaqus Face object.
           obj:  Abaqus Part object or Abaqus Partinstance object.

       Returns:
           Abaqus Edge object.

       Raises:
           None.
    """
    
    edge_id = face.getEdges()[0]
    return obj.edges[edge_id]



def get_all_vertices(obj: Any) -> list[list[geom.Point3D]]:
    """Gets all the vertices of an object.
       
       The documentation for the pointOn member of the Vertex object is 
           incorrect. It is really a tuple of tuple of floats, not a simple 
           tuple of floats.

       Args:
           obj: Abaqus Part object or Abaqus PartInstance object.

       Returns:
           The vertices of an object. Note that each inner list of vertices
               contains the vertices which are on a single face of the
               object.

       Raises:
           None.
    """

    faces = obj.faces 

    vertices = []
    for face in faces:
        vertices_on_face = []
        vertex_ids = face.getVertices()

        for vertex_id in vertex_ids:
            vertex = geom.Point3D(np.array(obj.vertices[vertex_id].pointOn[0]))
            vertices_on_face.append(vertex)

        vertices.append(vertices_on_face)

    return vertices 



def get_face_vertices(face: Any, part: Optional[Any] = None, 
                      assembly: Optional[Any] = None, 
                      part_instance: Optional[Any] = None
                     ) -> list[tuple[float, float, float], ...]:
    """Gets the vertices of a face.

       Error in documentation of pointOn member.

       Args:
           face:          Abaqus Face object. The face of interest.
           part:          Abaqus Part object. The part to which the face belongs.
           part_instance: Abaqus PartInstance object. The assembly to which the
                              face belongs.
           assembly:      Abaqus Assembly object. The assembly to which the face 
                              and part instance belong.
           
           Note that either the part argument or the part_instance and the
               assembly objects must be passed.

       Returns:
           List of vertices on the face. 

       Raises:
           None.
    """

    vertex_ids = face.getVertices()
    if assembly is not None and part_instance is not None:
        vertices = [part_instance.vertices[vertex_id].pointOn[0] for vertex_id in vertex_ids]
    elif part is not None:
        vertices = [part.vertices[vertex_id].pointOn[0] for vertex_id in vertex_ids]
    else:
        raise RuntimeError("Invalid combination of arguments passed!")

    return vertices 



def get_bc_cnt(model_name: str, mdb: Any) -> int:
    """Gets the number of boundary conditions present in a model."""

    return len(mdb.models[model_name].boundaryConditions)



# *****************************************************************************
#                         Abaqus Object Manipulation 
#    These functions manipulate the state/information that Abaqus maintains.
# *****************************************************************************


def create_mdb(name: str, path: str) -> Any:
    """Creates an MDB and opens it.
    
       Does not automatically save the MDB.

       Args:
           name: Name of the new MDB.
           path: Location of where a save (not a save as!) will occur.

       Returns:
           Abaqus MDB object.

       Raises:
           None.
    """

    return Mdb(path + name)



# Assign the only section in the model to the whole part.
def assign_only_section_to_part(part: Any, model_name: str, mdb: Any) -> None:
    """Assigns the only section in a model to the entirety of a part.

       This function creates a set which encompasses the whole part.

       Args:
           part:       Abaqus Part object.
           model_name: Name of the model.
           mdb:        MDB in which the model lives.

       Returns:
           None.

       Raises:
           None.
    """

    model = mdb.models[model_name]

    if len(model.sections) != 1:
        raise RuntimeError("There isn't exactly one section in the model!")

    section_name = model.sections.keys()[0]

    full_set = part.Set(name="simple_set", cells=part.cells)
    
    part.SectionAssignment(full_set, section_name)
    
    return



def use_mdb(path_to_mdb: str) -> Any:
    """Opens an MDB.
       
       Treat opening an MDB like opening a file. Only open it in one place at 
           one time and close it as soon as it is no longer needed.

       Args:
           path_to_mdb: Path to the MDB of interest.

       Returns:
           Abaqus MDB object.

       Raises:
           None.
    """
    
    return openMdb(path_to_mdb)



def close_mdb(mdb: Any) -> None:
    """Closes an MDB."""

    mdb.close()



def save_mdb_as(save_path: str, mdb: Any) -> None:
    """Saves an MDB somewhere.
     
       Args:
           save_path: Absolute path of desired save location. This should end
                         with a file with suffix ".cae".
           mdb:       Abaqus MDB object.

       Returns:
           None.

       Raises:
           None.
    """

    mdb.saveAs(save_path)    
    dp("The mdb was saved to: " + mdb.pathName)



def check_basic_geom(should_print: bool, mdb: Any) -> bool:
    """Checks, non-exhaustively, that an MDB contains the basic specification of
           a part geometry.

       Does not check every detail of the contents of the MDB.

       Args:
           should_print: Indicates if failures (i.e. the MDB contains stuff
                             which is unexpected) should cause debug information
                             to being printed.
           mdb:          Abaqus MDB object.

       Returns:
           False if the MDB does not contain a basic specification of a part
               geometry, True otherwise.

       Raises:
           None.
    """

    if not check_standard_model_and_part(should_print, mdb):
        return False 

    if not check_steps_loads_assembly(should_print, mdb) or \
       not check_no_material_and_section(should_print, mdb):
        return False

    return True



def check_multiple_steps(should_print: bool, model_name: bool, mdb: Any) -> bool:
    """Checks that a model contains multiple steps."""

    if model_name not in mdb.models:
        if should_print:
            dp("failure 1")
        return False

    if len(mdb.models[model_name].steps) > 1:
        if should_print:
            dp("failure 2")
        return True

    return False



def check_standard_model_and_part(should_print: bool, mdb: Any):
    """Checks that an MDB was just initialized and contains the default post-
           initialization objects."""

    if STANDARD_MODEL_NAME not in mdb.models or \
       len(mdb.models[STANDARD_MODEL_NAME].parts) != 1:
        if should_print:
            dp("failure 1")
        return False

    model = mdb.models[STANDARD_MODEL_NAME]

    if STANDARD_INIT_GEOM_PART_NAME not in model.parts:
        if should_print:
            dp("failure 2")
        return False

    return True



def check_no_material_and_section(should_print: bool, mdb: Any) -> bool:
    """Checks that the default model does not have a material, a section, or a
           section assignment."""

    model = mdb.models[STANDARD_MODEL_NAME]
    part = model.parts[STANDARD_INIT_GEOM_PART_NAME]

    if len(model.materials) != 0 or \
       len(model.sections) != 0 or \
       len(part.sectionAssignments) != 0: 
        if should_print:
            dp("failure 1")
        return False
    
    return True



def check_steps_loads_assembly(should_print: bool, mdb: Any) -> bool:
    """Checks that the default model has exactly one step, no loads, and that
          the assembly has no instances."""

    model = mdb.models[STANDARD_MODEL_NAME]

    if len(model.steps) != 1 or \
       len(model.loads) != 0:
        if should_print:
            dp("failure 1")
        return False

    if len(model.rootAssembly.instances) != 0:
        if should_print:
            dp("failure 2")
        return False

    return True



def suppress_feature(name: str, part: Any) -> None:
    """Suppresses a feature associated with a part.

       Args:
           name: The name of the feature.
           part: Abaqus Part object.

       Returns:
           None.

       Raises:
           None.
    """

    part.suppressFeatures((name, ))



def build_sketch_transform_from_face(face: Any, edge: Any, obj: Any) -> Any:
    """Specifies the location and orientation of a sketch (basically a plane of
           finite size) in space.
     
       It appears that when the MakeSketchTransform() method is called and an
           Abaqus Face object is passed in, the sketchUpEdge argument is required
           even though the documentation says otherwise.  

       Args:
           face: Abaqus Face object. Locates the sketch in space.
           edge: Abaqus Edge object or Abaqus DatumAxis object. Orients the 
                     sketch in space. This edge specifies the 
           obj:  Abaqus Part object or Abaqus Assembly object. The object that 
                     the transform is associated with.

       Returns:
           Abaqus Transform object which specifies the orientation.

           Note that an Abaqus Transform object has a single method called 
               matrix() that returns the matrix which contains the translation 
               and orientation of the sketch in space. Extracting this information 
               can be useful to map points in the global coordinate system into 
               the coordinate system of this sketch.

       Raises:
           None. 
    """

    return obj.MakeSketchTransform(face, sketchUpEdge=edge)



def extract_global_csys_to_sketch_csys(transform: Any) -> geom.CSys3D:
    """Converts the information associated with a transform to a custom format
           for a coordinate system.

       The format of the return value of the matrix() method for an Abaqus
           Transform object is undocumented.

       Args:
           transform: Abaqus Transform object.

       Returns:
           Custom format for a coordinate system.

       Raises:
           None.
    """

    rot_and_trans = transform.matrix()
   
    x_axis = geom.Vec3D(np.array(rot_and_trans[0:3]))
    y_axis = geom.Vec3D(np.array(rot_and_trans[3:6]))
    z_axis = geom.Vec3D(np.array(rot_and_trans[6:9]))
    trans = geom.Vec3D(np.array(rot_and_trans[9:12]))

    basis = geom.Basis3D(x_axis, y_axis, z_axis) 

    return geom.CSys3D(trans, basis)



def build_constrained_sketch(transform: Any, sketch_name: str, model_name: str, 
                             mdb: Any) -> Any:
    """Builds and orients a sketch somewhere in space.

       The sheetSize argument doesn't matter when the sketch is not being drawn
           by hand via a UI.
       Sketches are only associated with whole models.

       Args:
           transform:   Abaqus Transform object. Specifies the location and 
                            orientation of the sketch in space.
           sketch_name: String. Name of the new sketch.
           model_name:  String.
           mdb:         Abaqus MDB object.

       Returns:
           Abaqus ConstrainedSketch object.

       Raises:
           RuntimeError: The sketch failed to be constructed for some reason. An
                             exception is raised because at this level because
                             Abaqus will not raise an exception if the sketch
                             is not created.
    """

    if transform is not None:
        sketch = mdb.models[model_name].ConstrainedSketch(sketch_name, sheetSize=1, transform=transform)
    else:
        sketch = mdb.models[model_name].ConstrainedSketch(sketch_name, sheetSize=1)

    if sketch == None:
        raise RuntimeError("Failed to build sketch!!")

    return sketch



# Build wire feature. 
# 
# Notes:
#    When a wire feature is created in Abaqus and the spline option is chosen to
#       connect the points, by default, a cubic polynomial with continuous first
#       and second derivatives is used. No other splines can be used.
#    Assumes that the spline is not self-connecting.
#
# Arguments:
#    spline     - PlanarCubicC2Spline3D object.
#    part_name  - String.
#    model_name - String.
#    mdb        - Abaqus MDB object.
#
# Returns:
#    Abaqus Feature object associated with the newly created wire. 
def build_wire(spline, part_name, model_name, mdb):
# type: (geom.PlanarCubicC2Spline3D, str, str, Any) -> Any

    part = mdb.models[model_name].parts[part_name]

    # Create datum points for each point in the spline.
    dps = []
    for v in spline.v_list:
        feature = part.DatumPointByCoordinate(v.components())
        dps.append(part.datums[feature.id])

    # Use the datum points to construct the wire feature.
    return part.WireSpline(dps, mergeType=SEPARATE)



# Sweep a tool cross section along a wire to create a solid feature.
#
# Notes:
#    The primary failure mode of this function is self-intersection. If the
#       path of the tool pass has sufficiently sharp turns and the cross section
#       is too large, the resulting solid would have self-intersections. However,
#       if the feature were to have a self-intersection, Abaqus will fail to
#       create it.
# 
# Arguments:
#    tool_pass  - ToolPass object.
#    edge_array - Abaqus EdgeArray object.
#                 The sequence of edges along which to sweep the tool.
#    part_name  - String.
#    model_name - String.
#    mdb        - Abaqus MDB object.
#
# Return:
#    Abaqus Feature object that corresponds to the newly created solid. 
def sweep_tool_along_wire(tool_pass, edge_array, part_name, model_name, mdb):
# type: (tp.ToolPass, Any, str, str, Any) -> Any

    part = mdb.models[model_name].parts[part_name]

    # This datum axis defines the orientation of the tool with respect to the
    #    path that it follows.
    # For now, the datum axis is parallel to the principal y axis. This guarantees
    #    that it is easy to orient the tool so that it sits on top of the path.
    start_point = tool_pass.path.v_list[0].components()
    above_start_point = (start_point[0], start_point[1] + 1, start_point[2])
    feature = part.DatumAxisByTwoPoint(start_point, above_start_point)
    datum_axis = part.datums[feature.id]

    # Create the sketch that will be swept along the wire.
    sketch = build_constrained_sketch(None, part_name, model_name, mdb)

    # Draw the rectangle.
    sketch.rectangle((tool_pass.radius, 0), (-tool_pass.radius, tool_pass.length))

    sweep = part.SolidSweep(path=edge_array, profile=sketch, sketchUpEdge=datum_axis, sketchOrientation=RIGHT)

    return sweep



# Add the caps at both ends of the tool pass.
# 
# Notes:
#    This function adds features to a part which already exists.
#    This function assumes that the cross section of the tool has already been
#       extruded along the path that the tool takes.
#    This function assumes the orientation of the tool pass with respect to the
#       path that the tool takes.
# 
# Arguments:
#    tool_pass  - ToolPass object.
#    part_name  - String.
#    model_name - String.
#    mdb        - Abaqus MDB object.
#
# Returns:
#    None.
def add_tool_pass_caps(tool_pass, part_name, model_name, mdb):
# type: (tp.ToolPass, str, str, Any) -> None

    part = mdb.models[model_name].parts[part_name]

    start_point = tool_pass.path.v_list[0]
    face = part.faces.findAt(coordinates=start_point.components())
    build_cap("CAP_1", tool_pass, face, part_name, model_name, mdb)

    end_point = tool_pass.path.v_list[-1]
    face = part.faces.findAt(coordinates=end_point.components())
    build_cap("CAP_2", tool_pass, face, part_name, model_name, mdb)



# Add a cap to the face of a part.
# 
# Notes:
#    This function assumes a toolpath orientation with respect to the toolpass
#       path.
#    Adds a feature to the part.
# 
# Arguments:
#    name       - String.
#                 Name of the sketch created for the cap.
#    tool_pass  - ToolPass object.
#                 Dictates the geometry of the cap.
#    face       - Abaqus Face object.
#                 The face on which the cap will be created.
#                 This face should be planar.
#    part_name  - String.
#    model_name - String.
#    mdb        - Abaqus MDB object. 
#
# Returns:
#    None. 
def build_cap(name, tool_pass, face, part_name, model_name, mdb):
# type: (str, Any, Any, str, str, Any) -> None

    part = mdb.models[model_name].parts[part_name]

    point = face.getCentroid()[0]
    datum_feature = part.DatumPointByCoordinate(point)
    point_datum = part.datums[datum_feature.id]

    # Assuming orientation of toolpath with respect to path of tool. 
    above_point = (point[0], point[1] + 1, point[2])
    datum_feature = part.DatumPointByCoordinate(above_point)
    above_point_datum = part.datums[datum_feature.id]

    datum_feature = part.DatumAxisByTwoPoint(point_datum, above_point_datum)
    datum_axis = part.datums[datum_feature.id]

    # Position the sketch in space.
    transform = build_sketch_transform_from_face(face, datum_axis, part)

    # Build the sketch, create its axis of rotation, and draw the shape to be
    #    revolved. 
    sketch = build_constrained_sketch(transform, name, model_name, mdb)
    axis_of_revolution = sketch.ConstructionLine((0., 0.), (0., 1.))
    sketch.assignCenterline(axis_of_revolution)
    sketch.rectangle((tool_pass.radius, tool_pass.length/2), (0, -tool_pass.length/2))

    # Do the revolving.
    part.SolidRevolve(sketchPlane=face, sketchPlaneSide=SIDE1, sketchUpEdge=datum_axis, sketch=sketch, angle=180.0, 
                      sketchOrientation=RIGHT, flipRevolveDirection=OFF)



def build_tool_pass_part(name, tool_pass, model_name, mdb_metadata, mdb):
# type: (str, tool_pass.ToolPass, str, abq_md.AbaqusMdbMetadata, Any) -> Any

    # Create the part in the model. 
    part = mdb.models[model_name].Part(name=name, dimensionality=THREE_D, type=DEFORMABLE_BODY)

    # Update metadata for bookkeeping. 
    mdb_metadata.models_metadata[model_name].part_names.append(name)

    # Create the wire feature.
    build_wire(tool_pass.path, name, model_name, mdb)

    # Extract the Abaqus EdgeArray object which the containts the Abaqus Edge
    #    object that corresponds to the wire.
    # This must immediately follow the wire being created.
    edge = mdb.models[model_name].parts[name].edges[-1:]

    # Sweep the cross section of the cylinder along the wire.
    sweep_tool_along_wire(tool_pass, edge, name, model_name, mdb)

    # Add the rounded portions at both ends of the tool pass.
    add_tool_pass_caps(tool_pass, name, model_name, mdb)

    return part



def build_bounding_box_part(name, tool_pass, model_name, mdb_metadata, mdb):
# type: (str, tp.ToolPass, str, abq_md.AbaqusMdbMetadata, Any) -> Any

    # Create the part in the model. 
    part = mdb.models[model_name].Part(name=name, dimensionality=THREE_D, type=DEFORMABLE_BODY)

    # Update metadata for bookkeeping. 
    mdb_metadata.models_metadata[model_name].part_names.append(name)

    build_tool_pass_bounding_box(name, 3, 3, 3, tool_pass, model_name, mdb)

    return part



# Build bounding box around a tool pass. 
# 
# Notes:
#    The part should be created before this function is called.
# 
# Arguments:
#    name       - String.
#                 Part name. 
#    x_excess   - Float.
#                 Amount of excess in the x direction that the user wants for the
#                    bounding boxes. Expressed in units of the global coordinate
#                    system.
#    y_excess   - Float.
#                 Amount of excess in the y direction that the user wants for the
#                    bounding boxes. Expressed in units of the global coordinate
#                    system.
#    z_excess   - Float.
#                 Amount of excess in the z direction that the user wants for the
#                    bounding boxes. Expressed in units of the global coordinate
#                    system.
#    tool_pass  - ToolPass object.
#    model_name - String.
#    mdb        - Abaqus MDB object.
#
# Returns:
#    None 
def build_tool_pass_bounding_box(name, x_excess, y_excess, z_excess, tool_pass, model_name, mdb):
# type: (str, float, float, float, tp.ToolPass, str, Any) -> None

    bounding_box = tp.create_tool_pass_bounding_box(x_excess, y_excess, z_excess, tool_pass)

    build_spec_right_rect_prism_part(name, bounding_box, model_name, mdb)    



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



# Create part instance in the root assembly. 
def instance_part_into_assembly(instance_name, part, dependent, model_name, mdb):
# type: (str, Any, bool, str, Any) -> Any 
    
    return mdb.models[model_name].rootAssembly.Instance(instance_name, part, dependent=dependent)
    


# Cut instances to create a new part. 
#
# Notes:
#    Creates a new Abaqus Part object. 
#    Causes the instances to be suppressed.
#
# Arguments:
#    name               - String.
#                         The name of the resulting part.
#    instance_to_be_cut - Abaqus PartInstance object.
#    cutting_instances  - Tuple of Abaqus PArtInstance objects.
#                         The instances that do the cutting.
#    model_name         - String.
#    mdb_metadata       - CommittedToolPassPlanMetadata object.
#    mdb                - Abaqus MDB object.
#
# Returns:
#    Abaqus Part object. 
def cut_instances_in_assembly(name, instance_to_be_cut, cutting_instances, model_name, mdb_metadata, mdb):
# type: (str, Any, tuple[Any], str, abq_md.AbaqusMdbMetadata, Any) -> Any

    try:
        # Beware, the argument list ordering in the documentation for PartFrom
        #    BooleanCut() appears to be incorrect.
        part = mdb.models[model_name].rootAssembly.PartFromBooleanCut(name=name, instanceToBeCut=instance_to_be_cut, cuttingInstances=cutting_instances)
        mdb_metadata.models_metadata[model_name].part_names.append(name)
    except AbaqusException as e:
        dp("")
        dp("Failed to do cut operation.")
        dp("The exception message is " + str(e.args))
        dp("")
        raise

    return part



# Merge instances to create a new part.
#
# Notes:
#    Creates a new Abaqus Part object. 
#    Natively, this operation does not cause the instances to be merged. To
#       standardize the behavior, this function explicitly suppresses the instances.
#
# Arguments:
#    name               - String.
#                         The name of the resulting part.
#    instance_to_merge  - Tuple of Abaqus PartInstance objects.
#                         The instances to be merged.
#    keep_intersections - Boolean.
#                         Should the intersecting boundaries be kept after the merge?
#    model_name         - String.
#    mdb_metadata       - CommittedToolPassPlanMetadata object.
#    mdb                - Abaqus MDB object.
#
# Returns:
#    Abaqus Part object. 
def merge_instances_in_assembly(name, instances_to_merge, keep_intersections, model_name, mdb_metadata, mdb):

    try:
        part = mdb.models[model_name].rootAssembly.PartFromBooleanMerge(name, instances_to_merge, keepIntersections=keep_intersections)
        mdb_metadata.models_metadata[model_name].part_names.append(name)

        # Suppress the features used in the merge. 
        for instance in instances_to_merge:
            name = (instance.name, )
            mdb.models[model_name].rootAssembly.suppressFeatures(name)

    except AbaqusException as e:
        dp("")
        dp("Failed to do merge operation.")
        dp("The exception message is " + str(e.args))
        dp("")
        raise

    return part



def mesh(part_instance, size, model_name, mdb):
# type: (Any, float, str, Any) -> None

    mdb.models[model_name].rootAssembly.setMeshControls(part_instance.cells, elemShape=TET, technique=FREE)

    seq = (part_instance, )
    mdb.models[model_name].rootAssembly.seedPartInstance(seq, size, deviationFactor=.03, minSizeFactor=.0001, constraint=FREE)
    mdb.models[model_name].rootAssembly.generateMesh(regions=seq)

    while mdb.models[model_name].rootAssembly.getUnmeshedRegions() != None:
        if size > 1:
            size = size - 1
        else:
            size = float(size) / 2
            
            # DEBUG
            mdb.saveAs("/home/andlars/Desktop/RS_TP_Adaptation/experiments/experiments/test_initial_geometry_cae/small_mesh_necessary.cae")

        if size <= .1:
            raise AssertionError("Even with a miniscule global element size of " + str(size) + ", the mesh still failed to be generated!")

        dp("An attempt at meshing failed. Decreasing global element size to " + str(size) + " and giving it another go.") 
        mdb.models[model_name].rootAssembly.deleteSeeds(seq)
        mdb.models[model_name].rootAssembly.seedPartInstance(seq, size, deviationFactor=.03, minSizeFactor=.0001, constraint=FREE)
        mdb.models[model_name].rootAssembly.generateMesh(regions=seq)

    mesh_stats = mdb.models[model_name].rootAssembly.getMeshStats(regions=seq)
    dp("Meshing succeeded. The total number of tetrahedral elements is " + str(mesh_stats.numTetElems))



# Add a material to a model.
def create_material(material, model_name, mdb):
# type: (mp.Material, str, Any) -> None

    if isinstance(material, mp.ElasticMaterial):
        poissons_ratio = material.poissons_ratio
        youngs_modulus = material.youngs_modulus
        material = mdb.models[model_name].Material(STANDARD_MATERIAL_NAME)
        material.Elastic(((youngs_modulus, poissons_ratio), ), type=ISOTROPIC)
    else:
        raise AssertionError("No support for this type of material.")



# Create section based on default material name.
def create_section(model_name, mdb):
# type: (str, Any) -> None

    section_repo = mdb.models[model_name].sections
    assert(len(section_repo) == 0)

    mdb.models[model_name].HomogeneousSolidSection(STANDARD_SECTION_NAME, STANDARD_MATERIAL_NAME)



# Create an equilibrium step after another specified step.
def create_equilibrium_step(name, name_step_to_follow, model_name, mdb_metadata, mdb):
# type: (str, str, str, abq_md.AbaqusMdbMetadata, Any) -> None

    step = mdb.models[model_name].StaticStep(name, name_step_to_follow)
    
    mdb_metadata.models_metadata[model_name].step_names.append(name)

    return step 



# Note that the name of the job matches the name of the various files produced
#    by the job when it is run (the .odb, .dat, etc. files).
def create_job(model_name, mdb_metadata, mdb):
# type: (str, abq_md.AbaqusMdbMetadata, Any) -> None 

    job_name = STANDARD_JOB_PREFIX + str(get_model_cnt(mdb))

    # The "resultsFormat=BOTH" causes Abaqus to generate .odb and .sim files
    #    when the simulation runs. This functionality is currently undocumented
    #    in the Abaqus documentation.
    job = mdb.Job(job_name, model_name, resultsFormat=BOTH)
    
    mdb_metadata.models_metadata[model_name].job_name = job_name
    
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

    # print_job_messages(job)



# Check the messages produced by an analysis job.
def print_job_messages(job):
# type: (Any) -> None

    if len(job.messages) != 0:
        dp("The job with name " + job.name + " ran and some messages were received!")

    for message in job.messages:
        print_message_info(message)
    


# Dump the message information.
# Assumes an Abaqus Message object is passed to it.
def print_message_info(message):
# type: (Any) -> None

    dp("A message of type " + str(message.type) + " was received when the job ran!")
    for key in message.data:
        dp("For the key " + str(key) + " the associated value is " + str(message.data[key]))



# Create a new model with a single deformed part represented by an orphan mesh.  
#    The orphan mesh comes from an ODB.
#
# Notes:
#    Assumes that the deformed part resulted from the last frame of the last
#       step of the simulation which generated the ODB.
#
# Arguments:
#    part_name    - String.
#                   Desired new part name.
#    model_name   - String.
#                   Desired new model name.
#    path_to_odb  - String.
#                   Path to the ODB to use.
#    mdb_metadata - AbaqusMdbMetadata object.
#    mdb          - Abaqus MDB object.
#
# Returns:
#    None. 
def create_model_and_part_from_odb(part_name, model_name, path_to_odb, mdb_metadata, mdb):
# type: (str, str, str, abq_md.AbaqusMdbMetadata, Any) -> None

    model = mdb.Model(model_name)
    odb = odbAccess.openOdb(path=path_to_odb)
    model.PartFromOdb(part_name, odb, shape=DEFORMED)
    odb.close()

    # Do book keeping.
    mdb_metadata.add_model(model_name)
    mdb_metadata.models_metadata[model_name].part_names.append(part_name)



# Create a part as an orphan mesh via an ODB in a pre-existing model.
#
# Notes:
#    Assumes that the deformed part resulted from the last frame of the last
#       step of the simulation which generated the ODB.
#
# Arguments:
#    part_name    - String.
#                   Desired new part name.
#    model_name   - String.
#                   Existing new model name.
#    path_to_odb  - String.
#                   Path to the ODB to use.
#    mdb_metadata - AbaqusMdbMetadata object.
#    mdb          - Abaqus MDB object.
#
# Returns:
#    None. 
def create_part_from_odb(part_name, model_name, path_to_odb, mdb_metadata, mdb):
# type: (str, str, str, abq_md.AbaqusMdbMetadata, Any) -> None

    odb = odbAccess.openOdb(path=path_to_odb)
    mdb.models[model_name].PartFromOdb(part_name, odb, shape=DEFORMED)
    odb.close()

    # Do book keeping.
    mdb_metadata.models_metadata[model_name].part_names.append(part_name)



"""
# ***** DEPRECATED DUE TO ODB ISSUES *****

# Update all ODBs open in the current session. 
#
# Notes:
#    This function is necessary due to, what I think is, a bug in the Abaqus
#       codebase. The bug signature is a message like "ODB is out-of-date. Please
#       close and reopen all ODBs". This bug might happen even when a path to
#       the ODB is passed to a function like materialsFromOdb(). This makes
#       me think that Abaqus has cached a version of the ODB at the specified
#       path, and realizes that the ODB needs to be updated.
#
# Arguments:
#    None.
# 
# Returns:
#    None.
def update_session_odbs():
# type: () -> None

    for i in range(len(session.odbs.keys())):
        name = session.odbs.keys()[i]
        session.odbs[name].close()
        odbAccess.openOdb(name)
"""



"""
# ***** DEPRECATED DUE TO ODB ISSUES *****

# Propagate the material definition from an ODB into a model.
#
# Notes:
#    None.
#
# Arguments:
#    path_to_odb    - String.
#    odb_model_name - String. 
#                     The name of the model in the ODB which contains the material.
#    new_model_name - String.
#                     The name of the model that the material will be imported
#                        into.
#    mdb            - Abaqus MDB object.
#
# Returns:
#    None. 
def create_material_from_odb(path_to_odb, odb_model_name, new_model_name, mdb):
# type: (str, str, str, Any) -> None 

    assert(len(mdb.models[new_model_name].materials) == 0)

    # Note: Don't do the following!!
    # This often causes a segmentation fault because, I believe, it uses a
    #    cached ODB that is out of date and out of our control. Explicitly 
    #    opening and accessing the ODB is preferable.
    # materials = mdb.models[model_name].materialsFromOdb(path_to_odb)

    odb = odbAccess.openOdb(path_to_odb)

    materials = odb.models[odb_model_name].materials
    assert(len(materials) == 1)

    material = materials[materials.keys()[0]]
    youngs_modulus, poissons_ratio = check_material_properties(material)
    
    new_material = mdb.models[new_model_name].Material(STANDARD_MATERIAL_NAME)
    new_material.Elastic(((youngs_modulus), (poissons_ratio)))

    odb.close()
"""



# Ensure that the material has the expected properties.
#
# Notes:
#   An elastic, isotropic material with a Young's modulus and a Poisson's ratio
#      is expected. 
#
# Arguments:
#    material - Abaqus Material object.
#
# Returns:
#    Tuple of two floats.
def check_material_properties(material):
# type: (Any) -> None 

    assert(material.elastic.type == ISOTROPIC)
    assert(len(material.elastic.table) == 2)
    assert(len(material.elastic.table[0]) == 1)
    assert(len(material.elastic.table[1]) == 1)

    return (material.elastic.table[0], material)



"""
# ***** DEPRECATED DUE TO ODB ISSUES *****

# Propagate the section from an ODB into a model. 
#
# Notes:
#    This does not propagate section assignments.
#    Assumes that the material underlying the section already exists in the
#       new model.
#
# Arguments:
#    path_to_odb     - String.
#    odb_model_name  - String. 
#                      The name of the model in the ODB that the section comes
#                         from.
#    new_model_name  - String.
#                      The name of the model that the section will be imported
#                         into.
#    mdb             - Abaqus MDB object.
#
# Returns:
#    None. 
def create_section_from_odb(path_to_odb, odb_model_name, new_model_name, mdb):
# type: (str, str, str, Any) -> None 

    section_repo = mdb.models[new_model_name].sections
    assert(len(section_repo) == 0)

    # Note: Don't do the following!!
    # This often causes a segmentation fault because, I believe, it uses a
    #    cached ODB that is out of date and out of our control. Explicitly 
    #    opening and accessing the ODB is preferable.
    # new_sections = mdb.models[model_name].sectionsFromOdb(path_to_odb)

    odb = odbAccess.openOdb(path_to_odb)
    sections = odb.models[odb_model_name].sections
    assert(len(sections) == 1)

    section = sections[sections.keys()[0]]
    check_section_properties(section)

    mdb.models[new_model_name].HomogeneousSolidSection(STANDARD_SECTION_NAME, STANDARD_MATERIAL_NAME)

    odb.close()
"""



# Check that the section has the expected properties.
#
# Notes:
#    A homogeneous solid section corresponding is expected.
# 
# Arguments:
#    section - Abaqus section object.
# 
# Returns:
#    None.
def check_section_properties(section):
# type: (Any) -> None

    assert(section.name == STANDARD_SECTION_NAME)
    assert(section.material == STANDARD_MATERIAL_NAME)



# Assign a single section to an entire part. 
#
# Notes:
#    Assumes that the section has standard name.
#
# Arguments:
#    part_name    - String.
#    model_name   - String. 
#    mdb          - Abaqus MDB object.
#
# Returns:
#    None. 
def assign_section_to_whole_part(part_name, model_name, mdb):
# type: (str, str, Any) -> None

    assert(len(mdb.models[model_name].parts[part_name].sectionAssignments) == 0)
    assert(len(mdb.models[model_name].sections) == 1)
    
    # Create a set which encompasses all the cells in the part.
    all_cells = mdb.models[model_name].parts[part_name].cells
    set = mdb.models[model_name].parts[part_name].Set("entire_part", cells=all_cells)

    # Assign the section to that whole part.
    mdb.models[model_name].parts[part_name].SectionAssignment(region=set, sectionName=STANDARD_SECTION_NAME)



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
    #    returned. ACTUALLY, I'M NOT SURE IF THIS IS TRUE...
    point_on_face = face.pointOn
    face_seq = obj.faces.findAt(coordinates=point_on_face)
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
    #    will also accept an Abaqus Region object. However, the converse does not seem
    #    to be true. Therefore a Abaqus Region object really is required.
    # Set associateFace=FALSE and stitch=False to speed this up.
    return part.FaceFromElementFaces(region, stitch=False, associateFace=FALSE)



# Use a sequence of face features to add a solid feature to the part.
# Note that this may fail, and the Abaqus documentation doesn't give any hint
#    to the circumstances under which it may fail. Presumably it will fail if
#    if the face features don't obviously outline a solid.
def add_solid_from_faces(part):
# type: (Any) -> Any 

    feature = part.AddCells(part.faces)
    assert(feature != None)

    return feature 



# Remove redundant geometric features to simplify meshing.
# 
# Notes:
#    In general, it can be useful to call createVirtualTopology recursively on
#       a part. Each call may eliminate more geometric features.
#    Note that Abaqus > Abaqus/CAE > Using Toolsets > The Virtual Topology Toolset >
#       Creating virtual topology based on geometric features says that the
#       ignoreRedundantEntities removes edges that separate an otherwise planar
#       or curved surface. 
#    In general, introducing a virtual geometry can change the geometry significantly!
#       Beware!
# 
# Arguments:
#    part - Abaqus Part object.
#
# Returns:
#    None.
def add_virtual_topology(part):
# type: (Any) -> None

    try:

        ARBITRARY_LARGE_THRESHOLD = 1000

        # Very liberal!
        part.createVirtualTopology(ignoreRedundantEntities=True, mergeShortEdges=True, 
                                   shortEdgeThreshold=ARBITRARY_LARGE_THRESHOLD, 
                                   mergeSmallFaces=True, smallFaceAreaThreshold=ARBITRARY_LARGE_THRESHOLD)

        # Somewhat conservative!
        # part.createVirtualTopology(ignoreRedundantEntities=True, mergeSmallAngleFaces=True, smallFaceCornerAngleThreshold=.5)

        # Very conservative!
        # part.createVirtualTopology(ignoreRedundantEntities=True)

    except AbaqusException as e:
        # If a virtual topology does not remove any entities, then an exception
        #    is sometimes issued!
        dp("")
        dp("AbaqusException raised when trying to create a virtual topology!")
        dp("The exception has arguments: " + str(e.args))
        dp("Continuing without the virtual topology...")
        dp("")



# Convert a shell geometry to a solid geometry.
# 
# Notes:
#    Does not use Abaqus' built in shell to solid conversion utility, which isn't
#       very effective.
#    Assumes that the part is a well formed shell geometry and that the result is
#       is a solid. Note the checkGeometry() method is useless.  
# 
# Arguments:
#    part - Abaqus Part object.
#
# Returns:
#    None. 
def convert_shell_to_solid(part):

    LARGE_STITCH_TOLERANCE = 10

    assert len(part.cells) == 0, "Not a shell!"
    part.Stitch(edgeList=part.edges, stitchTolerance=LARGE_STITCH_TOLERANCE)
    assert len(part.cells) == 1, "Not a solid!"



"""
# SADLY checkGeometry() DOES NOT PRINT TO STDOUT OR STDERR...IT PRINTS TO SOME FILE
#    THAT HOOKS INTO ABAQUS' GUI. THERE IS NO WAY OF KNOWING WHAT FILE THIS IS.

# Check that a part is a shell.
# 
# Notes:
#    Even if a part is a shell, there is no guarantee that the shell forms a
#       closed solid.
# 
# Arguments:
#    part - Abaqus Part object.
#
# Returns:
#    Boolean. 
def check_shell_geometry(part):

    old = sys.stdout
    with tempfile.TemporaryFile() as f:
        sys.stdout = f

        part.checkGeometry(detailed=ON)

        info = f.readlines()
        dp(str(info))

    sys.stdout = old
"""



# Add some reference points to a root assembly.
# 
# Notes:
#    Abaqus ReferencePoint objects can, for example, be used to build an Abaqus 
#       Region object. However, this function only returns Abaqus Feature objects. 
#
# Arguments:
#    points     - List of Point3D objects.
#                 Points to add as reference points.
#    model_name - String.
#    mdb        - Abaqus MDB object.
#
# Returns:
#    List of Abaqus Feature objects. 
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
#       objects". I believe it actually takes an Abaqus FaceArray object. 
# 
# Arguments:
#    face     - Abaqus Face object. 
#               The face to partition.
#    edge     - Abaqus Edge object. 
#               Specifies the orientation of the sketch on the face. This edge 
#                  should have been used to orient the sketch on the face when 
#                  the sketch was constructed.
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
#    Abaqus Face object.
def partition_face_with_sketch(face, edge, sketch, assembly=None, instance=None, part=None):
# type: (Any, Any, Any, Optional[Any], Optional[Any], Optional[Any]) -> Any

    point = face.pointOn

    if assembly != None and instance != None:
        face_array = instance.faces.findAt(coordinates=point)
        face_cnt = len(instance.faces)
        feature = assembly.PartitionFaceBySketch(face_array, sketch, sketchUpEdge=edge)
        post_partition_face_cnt = len(instance.faces)

    elif part != None:
        face_array = part.faces.findAt(coordinates=point)
        face_cnt = len(part.faces)
        feature = part.PartitionFaceBySketch(face_array, sketch, sketchUpEdge=edge)
        post_partition_face_cnt = len(part.faces)

    else:
        raise RuntimeError("Invalid combination of arguments passed.")

    if post_partition_face_cnt != face_cnt + 1 or feature is None:
        raise AssertionError("Partitioning didn't work!")

    # Hack! Assumes that the created face is the last face in the array of faces. 
    if assembly != None and instance != None:
        new_face = instance.faces[-1]
    else:
        new_face = part.faces[-1]

    return new_face



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
        raise RuntimeError("Can't determine what face the ngon lives on!")

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
        line = sketch.Line(cur_point, next_point)

    return sketch



# Partitions the face of an object.
# 
# Notes:
#    The ngon must live on a face of the object. 
#    If the ngon exactly matches a face which already exists, no new face is created.
#       The existing face is returned.
#    The technique used to create the partition is delicate. In particular,
#       it assumes that there are no faces which overlap the new partition or are
#       fully contained inside of it. 
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
#    assembly      - Abaqus Assembly object.
#    instance      - Abaqus PartInstance object.
# 
# Optional Arguments Notes:
#    Either part or instance and assembly must be specified.
#
# Returns:
#    Abaqus Face object.
def partition_face(ngon, part=None, assembly=None, instance=None):
# type: (geom.NGon3D, Optional[Any], Optional[Any], Optional[Any]) -> Any

    module = part if part else assembly
    obj = part if part else instance

    face_ngon_belongs_to = find_face_ngon_lives_on(ngon, obj)

    """ DEPRECATED TECHNIQUE
    # If the face the ngon belongs to has vertices which match the ngon, we don't 
    #    want to try to partition the face and create a new face. This will cause 
    #    partitioning to fail and Abaqus to give up. 
    # TODO: The ngon could have redundant vertices, but such a partitioning will
    #    actually succeed because new vertices will be added to the face.
    # TODO: Assuming that the vertices which make up the Abaqus Face object are
    #    not redundant.
    # If they do share vertices, no partitioning needs to be done.
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
        new_face = partition_face_with_sketch(face_ngon_belongs_to, edge, sketch, part=part)
    else:
        new_face = partition_face_with_sketch(face_ngon_belongs_to, edge, sketch, assembly=assembly, instance=instance)
    """

    pre_partition_face_cnt = len(obj.faces)

    # Create a datum point for each vertex.
    datums = []
    for idx in range(len(ngon.vertices)):

        feature = module.DatumPointByCoordinate(ngon.vertices[idx].components())
        id = feature.id
        datums.append(module.datums[id])

    # Partition the face using the datum pairs to create edges.
    # Once all pairs of vertices of the ngon have been created, a new face should
    #    also be created.
    redundant_edge_cnt = 0
    for idx in range(len(datums)):
        if idx == len(datums) - 1:
            datum_pair = (datums[idx], datums[0])
        else:
            datum_pair = (datums[idx], datums[idx + 1])

        try:
            feature = module.PartitionFaceByShortestPath(face_ngon_belongs_to, datum_pair[0], datum_pair[1])
        except AbaqusException as e:
            # When partitioning, if the edge created matches an edge which already
            #    exists on the object, Abaqus will throw an exception. 
            redundant_edge_cnt += 1

            dp("")
            dp("Failed to create an edge during face partitioning. This may not be a bug.")
            dp("The arguments of the exception are " + str(e.args))
            dp("")
            
    post_partition_face_cnt = len(obj.faces)
    if not (post_partition_face_cnt == pre_partition_face_cnt + 1 or redundant_edge_cnt == len(ngon.vertices)):
        raise AssertionError("Partitioning failed horribly!!")

    # Find the face that matches the ngon.
    for face in obj.faces:

        vertex_ids = face.getVertices()
        vertices = [obj.vertices[id].pointOn[0] for id in vertex_ids]
        points = geom.seq_points(vertices)

        if geom.seq_points_equal(ngon.vertices, points):
            return face

    raise AssertionError("Failed to find new face!")



# *****************************************************************************
#                             DEPRECATED Functions. 
# *****************************************************************************


def get_all_vertices_ordered(obj: Any) -> list[list[list[geom.Point3D]]]:
    """DEPRECATED in deference to the get_closest_face() approach.

       Gets all the vertices of an object and ensures that the vertices are
           well-ordered.

       The documentation for the pointOn member of the Vertex object is 
           incorrect. It is really a tuple of tuple of floats, not a simple 
           tuple of floats.

       Consider the square face defined by the vertices: (0, 0), (1, 0), (0, 1), 
           (1, 1). The vertices are not well-ordered when the list of vertices 
           is [(0, 0), (1, 0), (0, 1), (1, 1)] because if you draw a line from 
           (0, 0) to (1, 0), a line from (1, 0) to (0, 1), and a line from (0, 1) 
           to (1, 1) you don't get a square.

       Args:
           obj: Abaqus Part object or Abaqus PartInstance object.

       Returns:
           All the vertices of the object. The outer list contains lists 
               that correspond to a single face. The inner lists contain
               list of vertices.
           The reason for this nesting is that faces can have holes. A single
               face might have an outer edge with some vertices and two
               holes defined by some other vertices. In this case, each of
               these groups of vertices is included in a single list.

       Raises:
           None.
    """

    vertices = []

    faces = obj.faces 
    for face in faces:
        vertices.append([])

        vertex_ids = face.getVertices()
        vertices_on_face = [obj.vertices[vertex_id] for vertex_id in vertex_ids]

        ordered_vertices_on_face = get_face_ordered_vertices(vertices_on_face, face, obj) 

        for vertex_group in ordered_vertices_on_face:
            vertices[-1].append([geom.Point3D(np.array(vertex.pointOn[0])) for vertex in vertex_group])

    return vertices 



def get_face_ordered_vertices(vertices: list[Any], face: Any, obj: Any) -> list[list[Any]]:
    """DEPRECATED in deference to the get_closest_face() approach.
    
       Orders the vertices of a face.

       A well ordered set of vertices reflects the actual connectivity of
          the lines which connect the vertices. For example, if a square is
          built by connecting a line between v1 and v2, another line from v2
          to v3, another line from v3 to v4, and a final line from v4 to v1,
          then the vertices must be ordered like (v1, v2, v3, v4), or any
          shifted version of this, to be considered well ordered.
       
       This function works regardless of the face being convex or non-convex.

       Args:
           vertices: List of Abaqus Vertex objects. The vertices to be ordered.
           face:     Abaqus Face object. The face to which the vertices belong.
           obj:      Abaqus Part object or Abaqus Part Instance object. The
                         object to which the vertices and the face belong.

       Returns:
           Nested lists of Abaqus Vertex objects. The outer list contains lists
               which each correspond to a single boundary. The inner lists
               contain Abaqus vertex objects. The ordering of the inner lists
               is arbitrary.

       Raises:
           None.
    """

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



def traverse_connected_vertices_on_face(vertex: Any, face: Any, obj: Any) -> list[Any]:
    """DEPRECATED in deference to the get_closest_face() approach.
    
       Traverses vertices on a face via the geometrical connections between
           the vertices.

       Args:
           vertex: Abaqus Vertex object. The vertex from which to start the
                       traversal.
           face:   Abaqus Face object. The face on which the vertex lives.
           obj:    Abaqus Part object or Abaqus PartInstance object. The object
                       to which the face and vertex belong.

       Returns:
           List of Abaqus Vertex objects.

       Raises:
           None.
    """

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



def get_neighbor_vertices(vertex, face, obj):
    """DEPRECATED in deference to the get_closest_face() approach.
    
       Gets the vertices connected to a single vertex on a face.

       Args:
           vertex: Abaqus Vertex object. The vertex of interest.
           face:   Abaqus Face object. The face that the vertex of interest
                       lies on.
           obj:    Abaqus Part object or Abaqus PartInstance object. The part
                       that the vertex and face belong to.

       Returns:
           List of two Abaqus Vertex objects.

       Raises:
           None.
    """

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



def check_face_ngon_match(ngon: geom.NGon3D, face: Any, obj: Any) -> bool:
    """DEPRECATED. This function was previously used in the process of partitioning
           a face. In particular, it was used to check if a face of the exact same
           geometry already existed. In the new approach, it does not matter if
           a face of the exact same geometry already exists.
    
       Checks if a face and a polygon have vertices which match exactly.

       Beware. A face and a polygon could have matching vertices but, due to
           curvature, they might look very different.
       Also, this function checks for an EXACT match without taking into
           account small dicrepencies due to floating point issues. 

       Args:
           ngon: NGon3D object. The polygon of interest.
           face: Abaqus Face object.
           obj:  Abaqus Part object or Abaqus PartInstance object. The object
                     to which the face belongs.

       Returns:
           True if a match, false otherwise.

       Raises:
           None.
    """

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



def change_only_part_name(new_name: str, model_name: str, mdb: Any) -> None:
    """DEPRECATED. Not sure where this was ever used.
       
       In a model with only a single part, changes the name of the part.

       Args:
           new_name:   The desired new name of the model.
           model_name: The current name of the model.
           mdb:        Abaqus MDB object. The MDB that the model belongs to.

       Returns:
           None.

       Raises:
           None.
    """

    if len(mdb.models[model_name].parts.keys() != 1):
        raise ValueError("There should only be a single part in the model.")

    old_name = mdb.models[model_name].parts.keys()[0]

    mdb.models[model_name].parts.changeKey(fromName=old_name, toName=new_name)




def check_ready_for_toolpasses(should_print: bool, mdb: Any) -> bool:    
    """DEPRECATED. Other MDB state checker functions are now being used. 
       
       Checks, non-exhaustively, that an MDB is ready to have tool passes
          simulated in it.

       Does not check every detail of the contents of the MDB.

       Args:
           should_print: Indicates if failures (i.e. the MDB contains stuff
                             which is unexpected) should cause debug information
                             to being printed.
           mdb:          Abaqus MDB object.

       Returns:
           False if the MDB is not ready, True if it is.

       Raises:
           None.
    """

    if not check_standard_model_and_part(should_print, mdb):
        return False 

    if not check_material_and_section(should_print, mdb) or \
       not check_steps_loads_assembly(should_print, mdb):
        return False

    return True



def check_material_and_section(should_print: str, mdb: Any) -> bool:
    """DEPRECATED. Helper for another deprecated function.

       Checks that the MDB contains exactly one material, one section, and that
           the section has been assigned."""

    model = mdb.models[STANDARD_MODEL_NAME]
    part = model.parts[STANDARD_INIT_GEOM_PART_NAME]

    if len(model.materials) != 1 or \
       len(model.sections) != 1 or \
       len(part.sectionAssignments) != 1: 
        if should_print:
            dp("failure 1")
        return False
    
    return True



def build_spec_right_rect_prism_part(part_name: str, spec_right_rect_prism: geom.SpecRightRectPrism, 
                                     model_name: str, mdb: Any) -> None:
    """DEPRECRATED. Was used to construct very simple tool pass geometries and
           bounding boxes.
    
       Builds a right rectangular prism with special characteristics.

       Args:
           part_name:             Name of the part. Also used to name the sketch 
                                      which is created.
           spec_right_rect_prism: The geometry of the part to create. 
           model_name:            The model in which to create the part.
           mdb:                   Abaqus MDB object.

       Returns:
           None.

       Raises:
           None.
    """

    # The part may be floating in space away from the origin.
    # This necessitates re-centering the sketch origin.
    # We want to re-center in terms of the x-y centroid and the z plane which
    #    is closer to z=0.
    centroid = spec_right_rect_prism.get_centroid()
    z_offset = spec_right_rect_prism.get_smaller_z()

    part = mdb.models[model_name].parts[part_name]

    # The datum plane is parallel to the x-y plane and offset by z_offset in
    #    the z direction.
    dp1_id = part.DatumPointByCoordinate((0, 0, z_offset)).id
    dp2_id = part.DatumPointByCoordinate((1, 0, z_offset)).id
    dp3_id = part.DatumPointByCoordinate((0, 1, z_offset)).id 
    dp1 = part.datums[dp1_id]
    dp2 = part.datums[dp2_id]
    dp3 = part.datums[dp3_id]

    # The up direction of the sketch is parallel to Abaqus' global y axis and needs
    #    to pass through the plane of the sketch.
    da1_id = part.DatumAxisByTwoPoint(dp1, dp3).id
    da1 = part.datums[da1_id]

    # Create the plane and retrieve the datum object associated with it.
    sketch_plane_id = part.DatumPlaneByThreePoints(dp1, dp2, dp3).id
    sketch_plane = part.datums[sketch_plane_id]

    # Create the transform to orient the sketch in space. 
    t = part.MakeSketchTransform(sketchPlane=sketch_plane, origin=(centroid.rep[0], centroid.rep[1], z_offset), sketchUpEdge=da1)

    # Create the sketch using the transform object.
    sketch = build_constrained_sketch(t, part_name, model_name, mdb)

    # Draw the rectangle accounting for the translated origin.
    v1, v2 = spec_right_rect_prism.get_rect_corners() 
    corner_1 = (v1.proj_xy().rep[0] - centroid.rep[0], v1.proj_xy().rep[1] - centroid.rep[1])
    corner_2 = (v2.proj_xy().rep[0] - centroid.rep[0], v2.proj_xy().rep[1] - centroid.rep[1])

    # Don't try to check the return value. This returns None even on success...
    sketch.rectangle(corner_1, corner_2)

    # Now extrude the sketch to create the solid.
    # ASSUMPTION HERE: We always select SIDE1 because that happens to work for the
    #    special right rectangular prism.  
    depth = spec_right_rect_prism.get_dims()[2] 
    part.SolidExtrude(sketch_plane, SIDE1, da1, sketch, depth=depth)




