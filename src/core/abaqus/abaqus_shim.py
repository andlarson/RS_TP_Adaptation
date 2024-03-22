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
       3) It's often necessary to track the state of the Abaqus MDB at a level
              which Abaqus does not natively provide. For example, Abaqus'
              API maintains lots of repository (aka dictionary) data structures
              about the state of the MDB. For example, there is a repository
              of models which belong to a single MDB, a repository of steps that
              belong to a single model, etc. However, these repository (aka
              dictionary) data structures do not record any information about
              ordering. For example, it's very useful to know the order of the
              steps in a model. It's therefore necessary to maintain custom data
              structures which record information about the state of the MDB.
              These data structures need to be updated when the state of the MDB
              changes. By funneling all uses of the Abaqus API through calls to
              functions in this file, and writing the functions in this file
              so that they require the custom data structures which maintain
              additional state to be passed in, the programmer need not remember
              to update the custom data structures when they interact with
              Abaqus' API.
       4) From an architectural point of view, it's helpful to separate the
              functionality that Abaqus provides from other functionality.

While all of the above is true, it's completely fine to use the Abaqus API 
    outside of this file. For example, it's often necessary to get an Abaqus 
    Model object by using an Abaqus MDB object, etc. The idea is to encapsulate
    complex and commonly used accesses to Abaqus' API into this file. However,
    note (3) above when doing so. Manipulating the state of an MDB outside
    of this file is dangerous. Reading the state of an MDB outside this file
    is fine.
"""

from typing import Any, Optional
import time
import random

from abaqus import *
from abaqusConstants import * 
import regionToolset 
import part
import odbAccess
import numpy as np

import src.core.material_properties.material_properties as mp
import src.core.tool_pass.tool_pass as tp
import src.util.geom as geom
import src.core.metadata.abaqus_metadata as abq_md
import src.core.boundary_conditions.boundary_conditions as bc

from src.util.debug import *


# *****************************************************************************
#                           Naming Stuff in Abaqus 
# These are names that Abaqus uses by default or are standardized custom names. 
# *****************************************************************************

# Names that Abaqus uses by default. 
STANDARD_MODEL_NAME = "Model-1"
STANDARD_INIT_GEOM_PART_NAME = "Initial_Geometry"
STANDARD_INITIAL_STEP_NAME = "Initial"
STANDARD_SECTION_NAME = "Section-1"
STANDARD_MATERIAL_NAME = "Material-1"
STANDARD_JOB_PREFIX = "Job-"
STANDARD_ORPHAN_MESH_FEATURE_NAME = "Orphan mesh-1"

# Custom names. 

# Generic.
STANDARD_MODEL_NAME_PREFIX = "Model-"
STANDARD_BC_PREFIX = "Boundary_Condition_"

# For tool pass simulations.
STANDARD_TOOL_PASS_PART_PREFIX = "Tool_Pass_"
STANDARD_POST_TOOL_PASS_PART_PREFIX = "Post_Tool_Pass_"
STANDARD_DEFORMATION_STEP_NAME = "Deform_Due_To_Material_Removal"
STANDARD_EQUIL_STEP_NAME = "Equilibrate_User_Specified_Stresses"
STANDARD_INIT_GEOM_WITH_BBOX_NAME = "Initial_Geometry_With_Bounding_Box"

# For surface traction application.
STANDARD_SURFACE_TRACTION_NAME = "Surface_Traction"
STANDARD_TRACTION_STEP_NAME = "Traction_Application"

# For gradient descent of volume difference function.
STANDARD_VOLUME_DIFFERENCE_MODEL_PREFIX = "Volume_Difference_"
STANDARD_POST_TRACTION_MODEL_PREFIX = "Post_Traction"
STANDARD_POST_TRACTION_PART_NAME = "Post_Traction"



# *****************************************************************************
#                         Abaqus Information Retrieval
#       These functions retrieve information/state that Abaqus maintains.
# *****************************************************************************


def _get_model_cnt(mdb: Any) -> int:
    """Gets the number of models in an Abaqus MDB.
    
       Args:
           mdb: Abaqus MDB object.

       Returns:
           Number of models in the MDB.

       Raises:
           None.
    """

    return len(mdb.models)



def _get_step_keyword(kwb: Any) -> int:
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



def get_closest_face(points: list[geom.Point3D], obj: Any,
                     want_seq: bool) -> Any:
    """Gets the face closest to the centroid of some points.

       This function is a workaround to writing/calling nasty geometric routines. 
           In particular, it was introduced to avoid writing a routine which 
           checks if a three dimensional polygon is entirely inside a three 
           dimensional polygon which may have holes in it.
       The default search tolerance is used.
       Note that Abaqus may throw an exception if no face is found.

       Args:
           points:   A nonzero number of points. 
           obj:      Abaqus Part object or Abaqus PartInstance object.
           want_seq: Flag which specifies the type of the return. When True, this
                         function will return an Abaqus Face Array object
                         consisting of exactly one face. When False, this
                         will return an Abaqus Face object.

       Returns:
           Abaqus Face object or Abaqus Face Array object, depending on the
               value of want_seq.

       Raises:
           None.
    """

    if len(points) == 0:
        raise ValueError("There should be a nonzero number of points.")

    centroid = geom.find_centroid(points)

    if not want_seq:
        face_and_point = obj.faces.getClosest([centroid.components()])
        face = face_and_point[0][0]
        return face
    else:
        face_seq = obj.faces.findAt(((centroid.components(), )))      
        assert len(face_seq) == 1, "Unexpectedly found more than one face!"
        return face_seq
    


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



def _get_all_vertices(obj: Any) -> list[list[geom.Point3D]]:
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



def _get_face_vertices(face: Any, part: Optional[Any] = None, 
                       assembly: Optional[Any] = None, 
                       part_instance: Optional[Any] = None
                      ) -> list[tuple[float, float, float]]:
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



def compute_part_volume(part: Any) -> int:
    """Computes the volume of a solid part.
                   
       Args:
           part: Abaqus Part object.
    
       Returns:
           The volume of the part.
    
       Raises:
           None.
    """

    return part.queryGeometry["volume"]



def check_simple_assembly(should_print: bool, mdb: Any) -> bool:
    """Checks that an MDB contains a single model with a single part instance
           in its assembly. No naming scheme is assumed."""

    if len(mdb.models) != 1:
        if should_print:
            dp("failure 1")
        return False
    
    model_name = mdb.models.keys()[0]
    if len(mdb.models[model_name].rootAssembly.instances) != 1:
        if should_print:
            dp("failure 2")
        return False

    return True



def check_first_model_simple_assembly(should_print: bool, mdb_md: abq_md.AbaqusMdbMetadata, 
                                      mdb: Any) -> bool:
    """Checks that the first model in an MDB contains a single part instance
           in its assembly. No naming scheme is assumed."""

    first_model_name = mdb_md.model_names[0]

    if len(mdb.models[first_model_name].rootAssembly.instances) != 1:
        if should_print:
            dp("failure 1")
        return False

    return True



def check_multiple_steps(should_print: bool, model_name: str, mdb: Any) -> bool:
    """Checks that a model contains multiple steps."""

    if model_name not in mdb.models:
        if should_print:
            dp("failure 1")
        return False

    if len(mdb.models[model_name].steps) == 1:
        if should_print:
            dp("failure 2")
        return False 

    return True 



def check_multiple_models(should_print: bool, mdb: Any) -> bool:
    """Checks if an MDB contains multiple models."""

    if len(mdb.models) <= 1: 
        if should_print:
            dp("failure 1")
        return False

    return True



def check_job_submissions(should_print: bool, mdb: Any) -> bool:
    """Checks if an MDB has at least one job submission."""

    if len(mdb.jobs) < 1:
        if should_print:
            dp("failure 1")
        return False

    return True



def check_two_part_model(should_print: bool, model_name: str, mdb: Any) -> bool:
    """Checks, non-exhasutively, that a model contains two parts with simple
           geometries."""

    if len(mdb.models[model_name].parts) != 2:
        if should_print:
            dp("failure 1")
        return False

    if not _check_no_material_and_section(should_print, model_name, mdb):
        if should_print:
            dp("failure 2")
        return False

    if not _check_steps_loads_assembly(should_print, model_name, mdb):
        if should_print:
            dp("failure 3")
        return False

    return True



def check_standard_model(should_print: bool, model: str, mdb: Any) -> bool:
    """Checks, non-exhaustively, that a model simply contains only a specification
           of a part geometry."""

    model = mdb.models[model]

    if STANDARD_INIT_GEOM_PART_NAME not in model.parts:
        if should_print:
            dp("failure 1")
        return False

    if len(model.boundaryConditions) != 0 or \
       len(model.steps) != 1 or \
       len(model.sections) != 0 or \
       len(model.materials) != 0 or \
       len(model.rootAssembly.instances) != 0:
        if should_print:
            dp("failure 2")
        return False

    # Check that the model has a geometric volume.
    # If the model does not have a geoemtric volume, it is probably an orphan
    #     mesh.
    if len(model.parts[STANDARD_INIT_GEOM_PART_NAME].cells) == 0:
        if should_print:
            dp("failure 3")
        return False

    return True
    


def _check_standard_model_and_part(should_print: bool, mdb: Any):
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



def _check_no_material_and_section(should_print: bool, model_name: str,
                                   mdb: Any) -> bool:
    """Checks that the model has neither a material nor a section."""

    model = mdb.models[model_name]

    if len(model.materials) != 0 or len(model.sections) != 0: 
        if should_print:
            dp("failure 1")
        return False
    
    return True



def _check_steps_loads_assembly(should_print: bool, model_name: str, 
                                mdb: Any) -> bool:
    """Checks that the model has exactly one step, no loads, and that
          the assembly has no instances."""

    model = mdb.models[model_name]

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



def check_simple_standard_mdb(should_print: bool, mdb: Any) -> bool:
    """Checks, non-exhaustively, that an MDB contains the basic specification of
           a part geometry and follows a standard naming scheme."""

    if not _check_standard_model_and_part(should_print, mdb):
        return False 

    if not _check_steps_loads_assembly(should_print, STANDARD_MODEL_NAME, mdb) or \
       not _check_no_material_and_section(should_print, STANDARD_MODEL_NAME, mdb):
        return False

    return True



def get_only_instance_name(model: Any) -> str:
    """Gets the name of the only part instance in the model. Assumed that the
           model contains only a single part instance."""
    
    assert len(model.rootAssembly.instances) == 1

    return model.rootAssembly.instances.keys()[0] 



# *****************************************************************************
#                         Abaqus Object Manipulation 
#    These functions manipulate the state/information that Abaqus maintains.
# *****************************************************************************



def create_mdb(name: str, path: str) -> Any:
    """Creates an MDB and opens it.
    
       Does not automatically save the MDB.

       Args:
           name: Name of the new MDB.
           path: Absolute path to directory in which MDB lives. 

       Returns:
           Abaqus MDB object.

       Raises:
           None.
    """

    return Mdb(os.path.join(path, name))



def copy_mdb(src: str, dest: str) -> None:
    """Copies an MDB. 
        
       Args:
           src:  Absolute path to MDB (.cae file). The MDB to copy. Assumed 
                     that this MDB is not already open.
           dest: Absolute path to desired save location. Assumed to end with
                     .cae.
    
       Returns:
           None.
    
       Raises:
           None.
    """

    to_copy = use_mdb(src)
    save_mdb_as(dest, to_copy)
    close_mdb(to_copy)



def assign_only_section_to_part(part: Any, model_name: str, mdb: Any) -> None:
    """Assigns the only section in a model to the entirety of a part.

       This function creates a set which encompasses the whole part.

       TODO: Redundant function.

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



def assign_section_to_whole_part(part_name: str, model_name: str, mdb: Any):
    """Assigns a single section to an entire part. 

       Assumes that the section has standard name, the part has no section 
           assignments, and there is exactly one section in the model.
       
       TODO: Redundant function.
       
       Args:
           part_name:  The name of the part to which to assign the section.
           model_name: Name of the model the part and section live in.
           mdb:        Abaqus MDB object.

       Returns:
           None.

       Raises:
           None.
    """

    if len(mdb.models[model_name].parts[part_name].sectionAssignments) != 0:
        raise RuntimeError("The part already has a section assignment!")

    if len(mdb.models[model_name].sections) != 1:
        raise RuntimeError("There is not exactly one section in the model!")
    
    # Create a set which encompasses all the cells in the part.
    all_cells = mdb.models[model_name].parts[part_name].cells
    set = mdb.models[model_name].parts[part_name].Set("entire_part", cells=all_cells)

    # Assign the section to that whole part.
    mdb.models[model_name].parts[part_name].SectionAssignment(region=set, sectionName=STANDARD_SECTION_NAME)



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

    dump_banner("MDB SAVED")
    dp("")
    dp("The MDB was saved to " + mdb.pathName)
    dp("")
    dump_banner_end()



def save_mdb(mdb: Any) -> None:
    """Saves the MDB to the location specified at creation time.
        
       Args:
           mdb: Abaqus MDB object.

       Returns:
           None.

       Raises:
           None.
    """

    mdb.save()

    dump_banner("MDB SAVED")
    dp("")
    dp("The MDB was saved to " + mdb.pathName)
    dp("")
    dump_banner_end()



def suppress_part_feature(name: str, part: Any) -> None:
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



def delete_assembly_feature(name: str, model: Any) -> None:
    """Deletes a feature associated with the assembly of a model.
                   
       Args:
           name: The name of the feature to delete.
           model: Abaqus Model object.
    
       Returns:
           None.
    
       Raises:
           None.
    """

    model.rootAssembly.deleteFeatures((name, )) 



def resume_assembly_feature(name: str, model: Any) -> None:
    """Resumes a feature associated with the assembly of a model. 
                   
       Args:
           name: The name of the feature.
           mdb:  Abaqus Model object.
    
       Returns:
           None.
    
       Raises:
           None.
    """
    
    model.resumeFeatures((name, ))



def _build_sketch_transform_from_face(face: Any, edge: Any, obj: Any) -> Any:
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



def _extract_global_csys_to_sketch_csys(transform: Any) -> geom.CSys3D:
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



def _build_constrained_sketch(transform: Any, sketch_name: str, model_name: str, 
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
                             exception is raised in this function because
                             Abaqus will not raise an exception if the sketch
                             is not created.
    """

    if transform is not None:
        sketch = mdb.models[model_name].ConstrainedSketch(sketch_name, sheetSize=1, transform=transform)
    else:
        sketch = mdb.models[model_name].ConstrainedSketch(sketch_name, sheetSize=1)

    if sketch is None:
        raise RuntimeError("Failed to build sketch!!")

    return sketch



def _build_wire(spline: geom.PlanarCubicC2Spline3D, part_name: str, model_name: str, 
               mdb: Any) -> Any:
    """Builds a wire feature.

       When a wire feature is created in Abaqus and the spline option is chosen to
           connect the points, by default, a cubic polynomial with continuous first
           and second derivatives is used. No other splines can be used.

       The wire cannot be self-intersecting.

       Args:
           spline:     The spline which defines the geometry of the wire.
           part_name:  The name of the pre-existing part which the wire will be 
                           associated with.
           model_name: The name of the pre-existing model which the wire will be 
                           associated with.
           mdb:        Abaqus MDB object.

       Returns:
           Abaqus Feature object associated with the newly created wire. 

       Raises:
           None.
    """

    part = mdb.models[model_name].parts[part_name]

    # Create datum points for each point in the spline.
    dps = []
    for v in spline.v_list:
        feature = part.DatumPointByCoordinate(v.components())
        dps.append(part.datums[feature.id])

    # Use the datum points to construct the wire feature.
    return part.WireSpline(dps, mergeType=SEPARATE)



def _sweep_tool_along_wire(tool_pass: tp.ToolPass, edge_array: Any, part_name: str, 
                          model_name: str, mdb: Any) -> Any:
    """Sweeps a tool cross section along a wire to create a solid feature.
       
       The primary failure mode of this function is self-intersection. If the
           path of the tool pass has sufficiently sharp turns and the cross section
           is too large, the resulting solid would have self-intersections. However,
           if the feature were to have a self-intersection, Abaqus will fail to
           create it.
       
       Currently, this function sweeps the tool along the wire such that the tool
           always sits above the spline of the wire (where above is in the +y
           direction).

       Args:
           tool_pass:  A representation of the tool pass which will be swept 
                           along the edges in the edge_array.
           edge_array: Abaqus EdgeArray object. The sequence of edges along 
                           which to sweep the tool. When a wire feature is created 
                           in a part, an edge is created which mirrors the wire. 
                           This single edge (which may be curved) is then often 
                           used by this function.
           part_name:  The name of the pre-existing part which the tool pass will 
                           be associated with.
           model_name: The name of the pre-exisitng mdoel which the tool pass will 
                           be associated with.
           mdb:        Abaqus MDB object.

       Returns:
           Abaqus Feature object assocaited with the newly created solid.

       Raises:
           None.
    """

    part = mdb.models[model_name].parts[part_name]

    # This datum axis defines the orientation of the tool with respect to the
    #    path that it follows.
    # For now, the datum axis is parallel to the principal y axis. This makes
    #    it easy to orient the tool so that it sits on top of the path.
    start_point = tool_pass.path.v_list[0].components()
    above_start_point = (start_point[0], start_point[1] + 1, start_point[2])
    feature = part.DatumAxisByTwoPoint(start_point, above_start_point)
    datum_axis = part.datums[feature.id]

    # Create the sketch that will be swept along the wire.
    sketch = _build_constrained_sketch(None, part_name, model_name, mdb)

    # Draw the rectangle.
    sketch.rectangle((tool_pass.radius, 0), (-tool_pass.radius, tool_pass.length))

    sweep = part.SolidSweep(path=edge_array, profile=sketch, sketchUpEdge=datum_axis, sketchOrientation=RIGHT)

    return sweep



def _add_tool_pass_caps(tool_pass: tp.ToolPass, part_name: str, model_name: str, 
                       mdb: Any) -> None:
    """Adds caps at both ends of a tool pass.
       
       A cap is just a half cylinder.

       This function adds features to a part which already exists.
       This function assumes that the cross section of the tool has already been
           extruded along the path that the tool takes.
       This function assumes the orientation of the tool pass with respect to the
           path that the tool takes.

       Args:
           tool_pass: The tool pass for which caps need to be added.
           part_name:  The name of the pre-existing part which the caps will 
                           be associated with.
           model_name: The name of the pre-exisitng model which the caps will 
                           be associated with.
           mdb:        Abaqus MDB object.

       Returns:
           None.

       Raises:
           None.
    """

    part = mdb.models[model_name].parts[part_name]

    start_point = tool_pass.path.v_list[0]
    face = part.faces.findAt(coordinates=start_point.components())
    _build_cap("CAP_1", tool_pass, face, part_name, model_name, mdb)

    end_point = tool_pass.path.v_list[-1]
    face = part.faces.findAt(coordinates=end_point.components())
    _build_cap("CAP_2", tool_pass, face, part_name, model_name, mdb)



def _build_cap(name: str, tool_pass: Any, face: Any, part_name: str, model_name: str, 
              mdb: Any) -> None:
    """Adds a cap (aka a half cylinder) to the face of a part.
    
       This function assumes a toolpath orientation with respect to the toolpass
           path.
       Adds a feature to the part.

       Args:
           name:       Name of the sketch created for the cap.
           tool_pass:  Dictates the geometry of the cap.
           face:       Abaqus Face object. The face on which the cap will be 
                           created. This face should be planar. If it is not
                           planar, the behavior of thsi function is undefined.
           part_name:  The name of the pre-existing part which the cap will 
                           be associated with.
           model_name: The name of the pre-exisitng model which the cap will 
                           be associated with.
           mdb:        Abaqus MDB object.

       Returns:
           None.

       Raises:
           None.
    """

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
    transform = _build_sketch_transform_from_face(face, datum_axis, part)

    # Build the sketch, create its axis of rotation, and draw the shape to be
    #    revolved. 
    sketch = _build_constrained_sketch(transform, name, model_name, mdb)
    axis_of_revolution = sketch.ConstructionLine((0., 0.), (0., 1.))
    sketch.assignCenterline(axis_of_revolution)
    sketch.rectangle((tool_pass.radius, tool_pass.length/2), (0, -tool_pass.length/2))

    # Do the revolving.
    part.SolidRevolve(sketchPlane=face, sketchPlaneSide=SIDE1, sketchUpEdge=datum_axis, 
                      sketch=sketch, angle=180.0, sketchOrientation=RIGHT, 
                      flipRevolveDirection=OFF)



def build_tool_pass_part(name: str, tool_pass: tp.ToolPass, model_name: str, 
                         mdb_metadata: abq_md.AbaqusMdbMetadata, mdb: Any
                        ) -> Any:
    """Builds the part associated with a tool pass.

       Args:
           name:         The name of the part this function creates. 
           tool_pass:    The tool pass that this function creates.
           model_name:   The model that the tool pass part is created in.
           mdb_metadata: The metadata associated with the MDB.
           mdb:          Abaqus MDB object.

       Returns:
           Abaqus Part object. 

       Raises:
           None.
    """

    # Create the part in the model. 
    part = mdb.models[model_name].Part(name=name, dimensionality=THREE_D, type=DEFORMABLE_BODY)

    # Update metadata for bookkeeping. 
    mdb_metadata.models_metadata[model_name].part_names.append(name)

    # Create the wire feature.
    _build_wire(tool_pass.path, name, model_name, mdb)

    # Extract the Abaqus EdgeArray object which the containts the Abaqus Edge
    #    object that corresponds to the wire.
    # This must immediately follow the wire being created.
    edge = mdb.models[model_name].parts[name].edges[-1:]

    # Sweep the cross section of the cylinder along the wire.
    _sweep_tool_along_wire(tool_pass, edge, name, model_name, mdb)

    # Add the rounded portions at both ends of the tool pass.
    _add_tool_pass_caps(tool_pass, name, model_name, mdb)

    return part



def inp_add_stress_subroutine(model_name: str, mdb: Any) -> None:
    """Modifies a .inp file so that a stress subroutine is included.

       Note that the stress subroutine must also be associated with the job.
           This is done elsewhere.

       Args:
           model_name: The name of the model for which to use the .inp file.
           mdb:        Abaqus MDB object.

       Returns:
           None.

       Raises:
           None.
    """

    # Synchronize the kwb to the current state of the model.
    # This must happen even if there have been no previous modifications to the
    #    input file.
    kwb = mdb.models[model_name].keywordBlock
    kwb.synchVersions()

    # The initial condition keyword is placed right before the first step in
    #    the input file.
    idx_step = _get_step_keyword(kwb)
    idx_before_step = idx_step - 1
    kwb.insert(idx_before_step, "*Initial Conditions, Type=Stress, User")



def inp_map_stress(path_sim_file: str, model_name: str, mdb: Any) -> None:
    """Modifies the input file so that the stress state of the part is initialized
           and mapped from the result of some other simulation.

       Args:
           path_sim_file: Absolute path to the simulation.
           model_name:    The name of the model with stress state which needs
                              to be modified.
           mdb:           Abaqus MDB object.

       Returns:
           None.

       Raises:
           None.
    """

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
    idx_step = _get_step_keyword(kwb)
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



def instance_part_into_assembly(instance_name: str, part: Any, dependent: bool, 
                                model_name: str, mdb: Any) -> Any:
    """Instances a part into the root assembly of a model.

       Args:
           instance_name: The name of the new instance.
           part:          Abaqus Part object. The part to be instanced.
           dependent:     Specifies if the instance should be dependent or
                              independent.
           model_name:    The name of the model in which the root assembly will
                              be used. Unclear if the part must also come from
                              this same model.
           mdb:           Abaqus MDB object.

       Returns:
           Abaqus PartInstance object.

       Raises:
           None.
    """
    
    return mdb.models[model_name].rootAssembly.Instance(instance_name, part, dependent=dependent)
    


def cut_instances_in_assembly(name: str, instance_to_be_cut: Any, cutting_instances: tuple[Any, ...], 
                              model_name: str, mdb_metadata: abq_md.AbaqusMdbMetadata, 
                              mdb: Any) -> Any:
    """Cuts instances (via boolean removal) to create a new part.

       Note that this function causes the instances used by it to be suppressed.

       Args:
           name:               The desired name of the resulting part.
           instance_to_be_cut: Abaqus PartInstance object. The instance which 
                                   will be cut (aka have material removed from 
                                   it).
           cutting_instances:  Tuple of Abaqus PartInstance objects. The instances 
                                   which will do the cutting. These instances define 
                                   where material should be removed from the 
                                   instance_to_be_cut.
           model_name:         The name of the model to which all the instances 
                                   should belong.
           mdb_metadata:       The metadata associated with this MDB.
           mdb:                Abaqus MDB object.

       Returns:
           Abaqus Part object.

       Raises:
           AbaqusException: This function re-raises exceptions of type 
                                AbaqusException which may be raised when the 
                                cut operation occurs. The most common
                                reason for the cut operation to fail and an 
                                AbaqusException to be raised is that the cutting 
                                instances do not overlap with the instance to be 
                                cut in any way.
    """

    try:
        # Beware, the argument list ordering in the documentation for PartFrom
        #    BooleanCut() appears to be incorrect.
        part = mdb.models[model_name].rootAssembly.PartFromBooleanCut(name=name, instanceToBeCut=instance_to_be_cut, cuttingInstances=cutting_instances)
        mdb_metadata.models_metadata[model_name].part_names.append(name)
    except AbaqusException as e:
        dump_exception()
        raise

    return part



def mesh(part_instance: Any, size: float, model_name: str, mdb: Any) -> None:
    """Meshes an indepedent part instance in the context of the root Assembly.

       Uses a simple algorithm to do meshing: Starts by trying to mesh with
           the passed seed size. If that fails, the seed size is decreased
           and meshing is given another go. The seed size is continually shrunk
           until meshing succeeds or it becomes smaller than some heuristically
           chosen quantity - at which point an exception is raised.

       Right now, tetrehedrons and the free meshing algorithm are being used.
           This combination yields the best chance at success for an arbitrary
           geometry.

       Args:
           part_instance: Abaqus PartInstance object. The instance to be meshed. 
                              Must be an independent instance.
           size:          The desired starting seed size. Approximately the target element 
                              size.
           model_name:    The model in which the part instance lives.
           mdb:           Abaqus MDB object.

       Returns:
           None.

       Raises:
           RuntimeError: If meshing continuously fails until the size becomes
                             smaller then .1, an exception is issued to prevent
                             exceeedingly dense meshes from being used.
    """

    mdb.models[model_name].rootAssembly.setMeshControls(part_instance.cells, elemShape=TET, technique=FREE)

    seq = (part_instance, )
    mdb.models[model_name].rootAssembly.seedPartInstance(seq, size, deviationFactor=.03, minSizeFactor=.0001, constraint=FREE)
    mdb.models[model_name].rootAssembly.generateMesh(regions=seq)

    while mdb.models[model_name].rootAssembly.getUnmeshedRegions() is not None:
        size = size / 2
            
        if size <= 1:
            raise RuntimeError("Even with a miniscule global element size of " + str(size) + ", the mesh still failed to be generated!")

        dp("An attempt at meshing failed. Decreasing global element size to " + str(size) + " and giving it another go.") 
        mdb.models[model_name].rootAssembly.deleteSeeds(seq)
        mdb.models[model_name].rootAssembly.seedPartInstance(seq, size, deviationFactor=.03, minSizeFactor=.0001, constraint=FREE)
        mdb.models[model_name].rootAssembly.generateMesh(regions=seq)

    mesh_stats = mdb.models[model_name].rootAssembly.getMeshStats(regions=seq)

    dump_banner("MESHING SUCCEEDED")
    dp("")
    dp("Meshing succeeded. The total number of tetrahedral elements is " + str(mesh_stats.numTetElems))
    dp("")
    dump_banner_end()



def create_material(material: mp.ElasticMaterial, model_name: str, mdb: Any):
    """Adds a material to a model.

       Right now, only some material types are acceptable. 

       Args:
           material:   The material to be added.
           model_name: The model to which the material should be added.
           mdb:        Abaqus MDB object.

       Returns:
           None.

       Raises:
           RuntimeError: Does not support the particular variety of material
                             passed.
    """

    if isinstance(material, mp.ElasticMaterial):
        poissons_ratio = material.poissons_ratio
        youngs_modulus = material.youngs_modulus
        abq_material = mdb.models[model_name].Material(STANDARD_MATERIAL_NAME)
        abq_material.Elastic(((youngs_modulus, poissons_ratio), ), type=ISOTROPIC)
    else:
        raise RuntimeError("No support for this type of material.")



def create_section(model_name: str, mdb: Any) -> None:
    """Creates section of material.

       Uses the default material name to pick out the material.

       Args:
           model_name: Name of the model in which to create the section.
           mdb:        Abaqus MDB object.

       Returns:
           None.

       Raises:
           RuntimeError: Some sections already exist or there is not exactly
                             one material that already exists.
    """

    section_repo = mdb.models[model_name].sections
    if len(section_repo) != 0:
        raise RuntimeError("Nonzero number of sections already exist!")

    material_repo = mdb.models[model_name].materials
    if len(material_repo) != 1:
        raise RuntimeError("A single material should already exist!")

    mdb.models[model_name].HomogeneousSolidSection(STANDARD_SECTION_NAME, STANDARD_MATERIAL_NAME)



def create_step(name: str, name_step_to_follow: str, model_name: str, 
                mdb_metadata: abq_md.AbaqusMdbMetadata, mdb: Any
               ) -> None:
    """Creates a step after some other step.

       Args:
           name:                Name of the new step.
           name_step_to_follow: Name of the step to follow.
           model_name:          Model in which to create the step.
           mdb_metadata:        Metadata for this MDB.
           mdb:                 Abaqus MDB object.

       Returns:
           None.

       Raises:
           None.
    """

    step = mdb.models[model_name].StaticStep(name, name_step_to_follow)
    
    mdb_metadata.models_metadata[model_name].step_names.append(name)

    return step 



def create_job(job_name: str, model_name: str, mdb_metadata: abq_md.AbaqusMdbMetadata, 
               mdb: Any) -> Any:
    """Creates a job via a standard naming scheme.

       Configures the job to emit a .sim file, in addition to the standard
           .odb file, when it runs.

       Args:
           job_name:     Name of the job. This dictates the names of the files
                             (.odb, .dat, etc.) produced by the job when it
                             runs.
           model_name:   Name of the model in which to create the job.
           mdb_metadata: Metadata for this MDB.
           mdb:          Abaqus MDB object.

       Returns:
           Abaqus ModelJob object.

       Raises:
           None.
    """

    # The "resultsFormat=BOTH" causes Abaqus to generate .odb and .sim files
    #    when the simulation runs. This functionality is currently undocumented
    #    in the Abaqus documentation.
    job = mdb.Job(job_name, model_name, resultsFormat=BOTH)
    
    mdb_metadata.models_metadata[model_name].job_name = job_name
    
    return job



def add_user_subroutine(job: Any, path_to_subroutine: str) -> None:
    """Associates a user subroutine with a job.

       Running a user subroutine requires that a user subroutine is associated
           with a job and that the user subroutine is invoked via a command 
           in the .inp file.

       Args:
           job:                Abaqus ModelJob object. The job to associate 
                                   the user subroutine.
           path_to_subroutine: An absolute path to the object (.o) file which
                                   contains the compiled user subroutine code.

       Returns:
           None.

       Raises:
           RuntimeError: Overwriting a user subroutine file which was previously
                             specified. There should only ever be a single user
                             subroutine file for a particular job!
    """
    
    if job.userSubroutine != "":
        raise RuntimeError("Overwriting a previous user subroutine file...")

    job.setValues(userSubroutine=path_to_subroutine)



def run_job(job: Any) -> None:
    """Runs a job and waits for it to complete.

       By default, prints messages received about the job. Note that when a
           script is run without the GUI, no job messages are ever returned.
           But Abaqus' technique for waiting for job completion relies on
           a completion message being returned by the job. What this means is
           that, when a script is run without the GUI, it is not possible
           to wait for the completion of a job submission. This is dangerous
           because, if a job is not completed and some output from the job
           is accessed, the behavior is undefined.

       Args:
           job: Abaqus ModelJob object. The job to run.

       Returns:
           None.

       Raises:
           RuntimeError: Job was aborted.        
    """

    job.submit()
    job.waitForCompletion()

    _print_job_messages(job)

    if job.status == ABORTED:
        raise RuntimeError("The job with name " + job.name + " was aborted!")



def _print_job_messages(job: Any) -> None:
    """Prints all the messages received during the analysis of a job."""
    
    dump_banner("MESSAGES RECEIVED DURING JOB ANALYSIS")
    dp("")

    if len(job.messages) != 0:
        dp("The job with name " + job.name + " ran and some messages were received!")
    
    for message in job.messages:
        dp("")
        dp("A message of type " + str(message.type) + " was received when the job ran!")
        for key in message.data:
            dp("For the key " + str(key) + " the associated value is " + str(message.data[key]))
    
    dp("")
    dump_banner_end()
    


def copy_model(model_to_copy: str, new_model: str, 
               mdb_metadata: abq_md.AbaqusMdbMetadata, mdb: Any) -> None:
    """Copies a model which already exists.
        
       Args:
           model_to_copy: The name of the model to copy. This model must already
                              exist in the MDB.
           model_name:    The name of the new model.
           mdb_metadata:  Metadata associated with the MDB.
           mdb:           Abaqus MDB object.
    
       Returns:
           None.
    
       Raises:
           None.
    """

    mdb.Model(name=new_model, objectToCopy=mdb.models[model_to_copy])

    # Do book keeping.
    mdb_metadata.add_copy(model_to_copy, new_model)



def create_model_and_part_from_odb(part_name: str, model_name: str, path_to_odb: str,  
                                   mdb_metadata: abq_md.AbaqusMdbMetadata, mdb: Any
                                  ) -> None:
    """Creates a new model with a single deformed part represented by an orphan mesh. 
           The orphan mesh comes from an ODB.

       Assumes that the deformed part resulted from the last frame of the last
           step of the simulation which generated the ODB.

       Args:
           part_name:    Desired new part name.
           model_name:   Desired new model name.
           path_to_odb:  Absolute path to ODB to use.
           mdb_metadata: Metadata for the MDB.
           mdb:          Abaqus MDB object.

       Returns:
           None.

       Raises:
           None.
    """

    model = mdb.Model(model_name)
    odb = odbAccess.openOdb(path=path_to_odb)
    model.PartFromOdb(part_name, odb, shape=DEFORMED)
    odb.close()

    # Do book keeping.
    mdb_metadata.add_model(model_name)
    mdb_metadata.models_metadata[model_name].part_names.append(part_name)



def create_part_from_odb(part_name: str, model_name: str, path_to_odb: str, 
                         mdb_metadata: abq_md.AbaqusMdbMetadata, mdb: Any
                        ) -> None:
    """Creates a part as an orphan mesh via an ODB.

       Assumes that the deformed part resulted from the last frame of the last
           step of the simulation which generated the ODB.

       Args:
           part_name:    Desired new part name.
           model_name:   Existing model name. 
           path_to_odb:  Absolute path to the ODB to use.
           mdb_metadata: Metadata from the MDB. 
           mdb:          Abaqus MDB object.

       Returns:
           None.

       Raises:
           None.
    """

    odb = odbAccess.openOdb(path=path_to_odb)
    mdb.models[model_name].PartFromOdb(part_name, odb, shape=DEFORMED)
    
    odb.close()

    # Do book keeping.
    mdb_metadata.models_metadata[model_name].part_names.append(part_name)



def build_region_with_face(face: Any, obj: Any) -> Any:
    """Builds a region which matches a face.

       Args:
           face: Abaqus Face object. The face with which to build the region.
           obj:  Abaqus Part object or Abaqus PartInstance object. The face must
                     be on this object. 

       Returns:
           Abaqus Region object.

       Raises:
           None.
    """

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
    


def build_region_with_elem_face(face: Any, part: Any) -> Any:
    """Builds a region from the face of an element in a mesh.

       Args:
           face: Abaqus MeshFace object. The face to build the region from.
           part: Abaqus Part object or Abaqus PartInstance object. The part that 
                     the face/mesh belongs to.

       Returns:
           Abaqus Region object.

       Raises:
           None.
    """

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
        assert(False)

    return region 



def add_face_from_region(region: Any, part: Any) -> Any:
    """Adds a face feature to a part based on a region.

       Args:
           region: Abaqus Region object. The region defining the new face.
           part:   Abaqus Part object or Abaqus PartInstance object. The part to
                       which the region belongs.

       Returns:
           Abaqus Feature object. The feature is associated with the new face.

       Raises:
           None.
    """
   
    # According to Scripting Reference > Python Commands > Region Commands >
    #    Region Object, whenever a command accepts a named set or surface, it
    #    will also accept an Abaqus Region object. However, the converse does not seem
    #    to be true. Therefore a Abaqus Region object really is required.

    # Set associateFace=FALSE and stitch=False to speed this up.
    return part.FaceFromElementFaces(region, stitch=False, associateFace=FALSE)

    # An alternative that does stitching. This makes face creation MUCH slower! 
    # return part.FaceFromElementFaces(region, stitch=True, associateFace=FALSE)



def add_solid_from_faces(part: Any) -> Any:
    """Adds a solid feature to a part based on the faces associated with the  
           part.

       Assumes that the part already has faces which define a solid region.

       Note that this may fail, and the Abaqus documentation doesn't give any 
           hint to the circumstances under which it may fail. Presumably it 
           will fail if the face features don't obviously outline a solid.
       
       Usually better to use convert_shell_to_solid() instead of this function.

       Args:
           part: Abaqus Part object.

       Returns:
           Abaqus Feature object. The feature associated with the solid.

       Raises:
           None.
    """

    feature = part.AddCells(part.faces)
    assert(feature is not None)

    return feature 



def add_virtual_topology(part: Any) -> None:
    """Removes redundant geometric features to simplify meshing.

       If virtual topology creation fails, nothing is done.

       In general, it can be useful to call createVirtualTopology recursively on
           a part. Each call may eliminate more geometric features.

       In general, introducing a virtual geometry can change the geometry 
           significantly! Beware! To this end, different policies can be used 
           to create the virtual topology, some more conservative than others.
           Basically, the more stuff (faces, edges, etc.) that is allowed to
           be combined, the more liberal the policy is. The danger of a very
           liberal policy is that the geometry can change significantly.

       Args:
           part: Abaqus Part object.

       Returns:
           None.

       Raises:
           None.
    """

    try:
        # Simplify the geometry as much as possible until no further simplification
        #     is possible. When no further simplifiction is possible, Abaqus will
        #     issue an exception.
        dump_banner("Geometry Simplification Started.")
        dp("")
        cnt = 0
        while True:
            dp("Attempting simplification number " + str(cnt))

            # Very liberal! This takes a long time to execute because Abaqus
            #     considers many possibilities in far-flung locations.
            # The only real constraint is the corner angle tolerance.
            # ARBITRARY_LARGE_NUMBER = 1000
            # ARBITRARY_SMALL_NUMBER = 1
            # part.createVirtualTopology(ignoreRedundantEntities=True, 
            #                            mergeShortEdges=True, shortEdgeThreshold=ARBITRARY_LARGE_NUMBER,
            #                            mergeSmallFaces=True, smallFaceAreaThreshold=ARBITRARY_LARGE_NUMBER,
            #                            mergeSliverFaces=True, faceAspectRatioThreshold=ARBITRARY_SMALL_NUMBER,
            #                            mergeSmallAngleFaces=True, smallFaceCornerAngleThreshold=ARBITRARY_SMALL_NUMBER,
            #                            mergeThinStairFaces=True, thinStairFaceThreshold=ARBITRARY_SMALL_NUMBER,
            #                            cornerAngleTolerance=60)

            # These magic numbers are the defaults that Abaqus uses when it
            #     automatically generates a virtual topology.
            part.createVirtualTopology(ignoreRedundantEntities=True, 
                                       mergeShortEdges=True, shortEdgeThreshold=.15, 
                                       mergeSmallFaces=True, smallFaceAreaThreshold=.11,
                                       mergeSliverFaces=True, faceAspectRatioThreshold=10.0,
                                       mergeSmallAngleFaces=True, smallFaceCornerAngleThreshold=10,
                                       mergeThinStairFaces=True, thinStairFaceThreshold=.03,
                                       cornerAngleTolerance=30)

            # Somewhat conservative!
            # part.createVirtualTopology(ignoreRedundantEntities=True, mergeSmallAngleFaces=True, smallFaceCornerAngleThreshold=.5)

            # Very conservative!
            # part.createVirtualTopology(ignoreRedundantEntities=True)

            dp("Simplification number " + str(cnt) + " succeeded.")
            cnt += 1

    except AbaqusException as e:
        dp("")
        dump_banner_end()
        
        # Dump the exception because no further exception is possible. 
        dump_exception()



def convert_shell_to_solid(part: Any) -> None:
    """Converts a shell geometry to a solid geometry.
       
       Results in a new solid feature being added to the part.

       Does not use Abaqus' built in shell to solid conversion utility, which isn't
           very effective. This is why the add_solid_from_faces() function is
           generally avoided. The stitching method used in this function is
           significantly more robust.

       Assumes that the part is a well formed shell geometry and that the result
           is a solid. Note that the checkGeometry() method is useless because
           the report it generates on the geometry is dumped to a file (which
           is neither stdout nor stderr) that hooks into Abaqus' GUI and is
           otherwise inpossible to locate.

       Args:
           part: Abaqus Part object.

       Returns:
           None.

       Raises:
           None.
    """

    if len(part.cells) != 0:
        raise RuntimeError("The part is not a shell.")

    LARGE_STITCH_TOLERANCE = 10
    part.Stitch(edgeList=part.edges, stitchTolerance=LARGE_STITCH_TOLERANCE)

    if len(part.cells) != 1:
        raise RuntimeError("The shell to solid conversion failed!")



def create_disp_rot_bc(BC_name: str, step_name: str, region: Any, 
                       settings: bc.DisplacementBCSettings, model_name: str, 
                       mdb: Any):
    """Creates displacement and rotational boundary conditions for some region.

       Args:
           BC_name:    Name of the new boundary condition.
           step_name:  Name of the step in which the boundary conditions are created.
                           Note that by default, the boundary conditions propagate
                           to all future steps.
           region:     Abaqus Region object. The region to apply the boundary condition.
           settings:   The boundary conditions for the region.
           model_name: The model in which to apply the boundary conditions.
           mdb:        Abaqus MDB object.

       Returns:
           None.

       Raises:
           None.
    """

    fix_x = SET if settings.fix_x else UNSET
    fix_y = SET if settings.fix_y else UNSET
    fix_z = SET if settings.fix_z else UNSET
    prevent_x_rot = SET if settings.prevent_x_rot else UNSET
    prevent_y_rot = SET if settings.prevent_y_rot else UNSET
    prevent_z_rot = SET if settings.prevent_z_rot else UNSET

    mdb.models[model_name].DisplacementBC(name=BC_name, createStepName=step_name, region=region, u1=fix_x, u2=fix_y, u3=fix_z, ur1=prevent_x_rot, ur2=prevent_y_rot, ur3=prevent_z_rot)



def _find_face_ngon_lives_on(ngon: geom.NGon3D, obj: Any) -> Any:
    """Finds the face of an object that an ngon exists on.
    
       If the ngon does not exactly live on the face of an object, an arbitrary
           face my be chosen! This is dangerous - the face may not be the face
           that the ngon lives on. To prevent this situation, the chosen face is
           compared to the face closest to the ngon. If they don't match, an
           exception is raised.

       Args:
           ngon: The ngon that lives on a face.
           obj:  Abaqus Part object or Abaqus PartInstance object.

       Returns:
           Abaqus Face object. The face that the ngon exists on.

       Raises:
           RuntimeError: The face that the ngon is detected on is not the face
                             closest to the ngon in space.
    """

    # Check that the vertices of the ngon do live on the plane DEFINED by some face.
    # Note that, in general, even if the vertices of the ngon live on a plane 
    #    DEFINED by some face, it does not necessarily mean that the ngon actually
    #    lives on that face. That face might have holes in it! That face could be
    #    curved!
    on_a_face = False
    for vertices_single_face in _get_all_vertices(obj):
        if geom.on_plane_of_ngon(vertices_single_face, ngon):
            on_a_face = True
    
    if not on_a_face:
        raise RuntimeError("Can't determine what face the ngon lives on!")

    # Assume that the face closest to points which make up the ngon is the face
    #    that the ngon lives on.
    face_ngon_belongs_to = get_closest_face(ngon.vertices, obj, False)

    return face_ngon_belongs_to



def partition_face(ngon: geom.NGon3D, part: Optional[Any] = None, 
                   assembly: Optional[Any] = None, instance: Optional[Any] = None
                  ) -> Any:
    """Partitions the face of an object using an ngon which is located on one
           of the faces of the object.

       The ngon must live on a face of the object. 
       If the ngon exactly matches a face which already exists, no new face is created.
           The existing face is returned.
       The technique used to create the partition is delicate. In particular,
           it assumes that there are no faces which overlap the new partition or are
           fully contained inside of it. 

       Args:
           ngon:          Lies in global coordinate system and should be on 
                             a face of the object which will be partitioned.
           new_face_name: The desired name of the new face.
           model_name:    The name of the model.
           mdb:           Abaqus MDB object.
           part:          Abaqus Part object. The part which has a face that 
                              the ngon lies on.
           assembly:      Abaqus Assembly object. The assembly which has the 
                              instance.
           instance:      Abaqus PartInstance object. The instance which has 
                              a face that the ngon lies on.

           Either the part or the assembly and instance arguments must be passed.

       Returns:
           Abaqus Face object.

       Raises:
           AssertionError: The face failed to be created or the face was created
                               but could not be found.
    """

    module = part if part else assembly
    obj = part if part else instance

    face_ngon_belongs_to = _find_face_ngon_lives_on(ngon, obj)

    """ DEPRECATED TECHNIQUE. This is an old technique used to create the
            partition.

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
    transform = _build_sketch_transform_from_face(face_ngon_belongs_to, edge, assembly)
    sketch = _build_constrained_sketch(transform, new_face_name, model_name, mdb)

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
    for v in ngon.get_builtin_rep():
        feature = module.DatumPointByCoordinate(v)
        id_ = feature.id
        datums.append(module.datums[id_])

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
            
            dump_banner("FACE PARTITIONING FAILURE")
            dp("")
            dp("Failed to create an edge during face partitioning. This may not be a bug.")
            dp("The arguments of the exception are " + str(e.args))
            dp("")
            dump_banner_end()
            
    post_partition_face_cnt = len(obj.faces)
    if not ((post_partition_face_cnt == pre_partition_face_cnt + 1) or (redundant_edge_cnt == len(ngon.vertices))):
        raise AssertionError("Partitioning failed horribly!!")

    # Find the face that matches the ngon.
    for face in obj.faces:

        vertex_ids = face.getVertices()
        vertices = [obj.vertices[id].pointOn[0] for id in vertex_ids]
        points = geom.seq_points(vertices)

        if geom.seq_points_equal(ngon.vertices, points):
            return face

    raise AssertionError("Failed to find new face!")



def create_traction(name: str, step_name: str, traction: geom.Vec3D, face_seq: Any, 
                    model_name: str, mdb: Any) -> None:
    """Creates a surface traction on a face. 
        
       Args:
           name:       The name of the surface traction which is created. This
                           name can be used to refer back to the surface traction
                           at a later point in time.
           step_name:  The name of the step in which to apply the traction.
           traction:   The traction vector to apply.
           face_seq:   Abaqus Face Array object consisting of exactly one face. 
                           The face must be a face on a part instance in the root 
                           assembly, not a face from a part. This is because 
                           tractions are created in the load application step, 
                           which operates on the root assembly. This face is 
                           the face on which the traction is applied.
           model_name: The name of the model in which the traction is created.
           mdb:        Abaqus MDB object.
    
       Returns:
           None.
    
       Raises:
           None.
    """

    if len(face_seq) != 1:
        raise RuntimeError("A traction can only be applied to a single face.")

    # The chosen face does not matter.
    # This surface gets a default name.
    surface = mdb.models[model_name].rootAssembly.Surface(side1Faces=face_seq, name=name)

    mdb.models[model_name].SurfaceTraction(name=name, createStepName=step_name,
                                           region=surface, magnitude=traction.len(),
                                           directionVector=((0, 0, 0), traction.components()),
                                           distributionType=UNIFORM, traction=GENERAL)



def merge_instances_in_assembly(name: str, instances_to_merge: tuple[Any, ...], 
                                keep_intersections: bool, model_name: str, 
                                mdb_metadata: abq_md.AbaqusMdbMetadata, mdb: Any
                               ) -> Any:
    """Merges instances to create a new part.

       Explicitly suppresses the instances used to do the merge.
       
       Args:
           name:               The name of the resulting part.
           instance_to_merge:  Tuple of Abaqus PartInstance objects. The instances 
                                   to be merged.
           keep_intersections: Should the intersecting boundaries be kept after 
                                   the merge?
           model_name:         Model in which the instances to be merged exist. 
           mdb_metadata:       Metadata for this MDB.
           mdb:                Abaqus MDB object.

       Returns:
           Abaqus Part object.

       Raises:
           None.
    """

    try:
        part = mdb.models[model_name].rootAssembly.PartFromBooleanMerge(name, instances_to_merge, keepIntersections=keep_intersections)
        mdb_metadata.models_metadata[model_name].part_names.append(name)

        # Suppress the features used in the merge. 
        for instance in instances_to_merge:
            abq_name = (instance.name, )
            mdb.models[model_name].rootAssembly.suppressFeatures(abq_name)

    except AbaqusException as e:
        dump_exception()
        raise

    return part



# *****************************************************************************
#                                 DEPRECATED 
# *****************************************************************************


def _get_all_vertices_ordered(obj: Any) -> list[list[list[geom.Point3D]]]:
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

        ordered_vertices_on_face = _get_face_ordered_vertices(vertices_on_face, face, obj) 

        for vertex_group in ordered_vertices_on_face:
            vertices[-1].append([geom.Point3D(np.array(vertex.pointOn[0])) for vertex in vertex_group])

    return vertices 



def _get_face_ordered_vertices(vertices: list[Any], face: Any, obj: Any) -> list[list[Any]]:
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
            
            if not vertices_used[i]:
              
                # Starting from an unvisited vertex, find all the other connected 
                #    vertices on the face.
                vertex_group = _traverse_connected_vertices_on_face(vertices[i], face, obj)
                per_group_vertices.append(vertex_group)

                # Mark all the visited vertices as visited.
                for j, vertex in enumerate(vertices):
                    if vertex in vertex_group:
                        vertices_used[j] = True

    return per_group_vertices 



def _traverse_connected_vertices_on_face(vertex: Any, face: Any, obj: Any) -> list[Any]:
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
    first_neighbors = _get_neighbor_vertices(first_vertex, face, obj)
    current_vertex = first_neighbors[0]
    
    well_ordered_vertices = [first_vertex]    
    while current_vertex.index != first_vertex.index:
        
        well_ordered_vertices.append(current_vertex)
    
        current_neighbors = _get_neighbor_vertices(current_vertex, face, obj)
    
        if current_neighbors[0].index != prev_vertex.index:
            prev_vertex = current_vertex
            current_vertex = current_neighbors[0]
        else:
            prev_vertex = current_vertex
            current_vertex = current_neighbors[1]

    return well_ordered_vertices



def _get_neighbor_vertices(vertex, face, obj):
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



def _check_face_ngon_match(ngon: geom.NGon3D, face: Any, obj: Any) -> bool:
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

    face_vertices = _get_face_vertices(face, obj)
    ngon_vertices = ngon.get_builtin_rep()
    ngon_face_share_vertices = True 
    if len(face_vertices) == len(ngon_vertices):
        for face_vertex in face_vertices:
            if face_vertex not in ngon_vertices:
                ngon_face_share_vertices = False
    else:
        ngon_face_share_vertices = False

    return ngon_face_share_vertices



def _change_only_part_name(new_name: str, model_name: str, mdb: Any) -> None:
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




def _check_ready_for_toolpasses(should_print: bool, mdb: Any) -> bool:    
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

    if not _check_standard_model_and_part(should_print, mdb):
        return False 

    if not _check_material_and_section(should_print, mdb) or \
       not _check_steps_loads_assembly(should_print, mdb):
        return False

    return True



def _check_material_and_section(should_print: bool, mdb: Any) -> bool:
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



def _build_spec_right_rect_prism_part(part_name: str, spec_right_rect_prism: geom.SpecRightRectPrism, 
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
    sketch = _build_constrained_sketch(t, part_name, model_name, mdb)

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



def _build_bounding_box_part(name: str, tool_pass: tp.ToolPass, model_name: str, 
                            mdb_metadata: abq_md.AbaqusMdbMetadata, mdb: Any
                           ) -> Any:

    """DEPRECATED. Was used to build bounding boxes around parts to introduce
           variable mesh density before variable mesh densities were realized
           from the built-in mesh functionalities.

       Builds a bounding box part around a tool pass. The tolerances around the
           tool pass are chosen heuristically.

       Args:
           name:         The name of the new part.
           tool_pass:    The tool pass around which the bounding box should be 
                             built.
           model_name:   The pre-existing model in which the bounding box part 
                             should be built.
           mdb_metadata: The metadata for this MDB.
           mdb:          Abaqus MDB object.

       Returns:
           Abaqus Part object.

       Raises:
           None.
    """

    # Create the part in the model. 
    part = mdb.models[model_name].Part(name=name, dimensionality=THREE_D, type=DEFORMABLE_BODY)

    # Update metadata for bookkeeping. 
    mdb_metadata.models_metadata[model_name].part_names.append(name)

    _build_tool_pass_bounding_box(name, 3, 3, 3, tool_pass, model_name, mdb)

    return part



def _build_tool_pass_bounding_box(name: str, x_excess: float, y_excess: float,
                                 z_excess: float, tool_pass: tp.ToolPass, model_name: str,
                                 mdb: Any) -> None:
    """DEPRECATED. Was used as a helper function for creating bounding boxes
           around tool passes.
       
       Builds bounding box around a tool pass.

       Args:
           name:       The name of the new part.
           x_excess:   Amount of excess in the x direction that the user wants 
                           for the bounding boxes. Expressed in units of the 
                           global coordinate system.
           y_excess:   Amount of excess in the y direction that the user wants 
                           for the bounding boxes. Expressed in units of the 
                           global coordinate system.
           z_excess:   Amount of excess in the z direction that the user wants 
                           for the bounding boxes. Expressed in units of the 
                           global coordinate system.
           tool_pass:  The tool pass around which to create the bounding box.
           model_name: The name of the model that the part should be created in.
           mdb:        Abaqus MDB object.

       Returns:
           None.

       Raises:
           None.
    """

    bounding_box = tp.create_tool_pass_bounding_box(x_excess, y_excess, z_excess, tool_pass)

    _build_spec_right_rect_prism_part(name, bounding_box, model_name, mdb)    



def _update_session_odbs() -> None:
    """DEPRECATED. Was used in an attempt to mitigate buggy accesses to out-of-
           date ODBs. In particular, was used in an attempt to propagate
           material definitions across simulations.

       Updates all ODBs open in the current session. 

       This function is necessary due to, what I think is, a bug in the Abaqus
           codebase. The bug signature is a message like "ODB is out-of-date. Please
           close and reopen all ODBs". This bug might happen even when a path to
           the ODB is passed to a function like materialsFromOdb(). This makes
           me think that Abaqus has cached a version of the ODB at the specified
           path, and realizes that the ODB needs to be updated.

       Args:
           None.

       Returns:
           None.

       Raises:
           None.
    """

    for i in range(len(session.odbs.keys())):
        name = session.odbs.keys()[i]
        session.odbs[name].close()
        odbAccess.openOdb(name)



def _create_material_from_odb(path_to_odb: str, odb_model_name: str, 
                             new_model_name: str, mdb: Any) -> None:
    """DEPRECATED due to bug. Often when this function was called immediately 
           after a job finished, a "ODB is out-of-date. Please
           close and reopen all ODBs" bug message was emitted. Unclear why this
           happened because all jobs were allowed to completely finish before
           accessing the ODB. It may have been the case that I was running scripts
           without the interactive flag, in which case the job may not have
           finished before an attempt to access the ODB was made.

       Propagates a material definition from an ODB into a model.

       Assumes that there are no materials in the new model.

       Args:
           path_to_odb:    Absolute path to the ODB. 
           odb_model_name: The name of the model in the ODB which contains the 
                               material.
           new_model_name: The name of the model that the material will be imported
                               into.
           mdb:            Abaqus MDB object.

       Returns:
           None.

       Raises:
           RuntimeError: Material failed to be created by this function.
    """

    # Note: Don't do the following!!
    # This often causes a segmentation fault because, I believe, it uses a
    #    cached ODB that is out of date and out of our control. Explicitly 
    #    opening and accessing the ODB is preferable.
    # materials = mdb.models[model_name].materialsFromOdb(path_to_odb)

    odb = odbAccess.openOdb(path_to_odb)
    materials = odb.models[odb_model_name].materials

    if len(materials) != 0:
        raise RuntimeError("There probably shouldn't be any materials!")

    material = materials[materials.keys()[0]]
    youngs_modulus, poissons_ratio = _check_material_properties(material)
    
    new_material = mdb.models[new_model_name].Material(STANDARD_MATERIAL_NAME)
    new_material.Elastic(((youngs_modulus), (poissons_ratio)))

    assert(len(materials) == 1)

    odb.close()



def _check_material_properties(material: Any) -> tuple[float, float]:
    """DEPRECATED. Was a helper function for exporting a material definition
           from an ODB to a new model. When I moved this function to deprecated
           status, I found a bug in it which may have been the reason why
           exporting the material from an ODB was failing.
    
       Checks that a material has some expected properties.

       Args:
           material: Abaqus Material object.

       Returns:
           Tuple containing Young's Modulus and Poisson's Ratio.

       Raises:
           RuntimeError: The material is missing some properties. 
    """

    if material.elastic.type != ISOTROPIC:
        raise RuntimeError("Not isotropic!")

    if len(material.elastic.table) != 2:
        raise RuntimeError("Missing material constants!")

    if len(material.elastic.table[0]) != 1 or len(material.elastic.table[1]) != 1:
        raise RuntimeError("Missing material constants!")

    return (material.elastic.table[0][0], material.elastic.table[1][0])



def _create_section_from_odb(path_to_odb: str, odb_model_name: str, 
                            new_model_name: str, mdb: Any) -> None:
    """DEPRECATED due to the same bug as create_material_from_odb.
       
       Propagates a single section from an ODB into a model. 
    
       Does not propagate section assignments.
       Assumes that the material underlying the section already exists in the
           new model.

       Args:
           path_to_odb:    Absolute path to an ODB. 
           odb_model_name: The name of the model in the ODB that the section comes
                               from.
           new_model_name: The name of the model that the section will be imported
                               into.
           mdb:            Abaqus MDB object.

       Returns:
           None.

       Raises:
           None.
    """

    sections = mdb.models[new_model_name].sections

    if len(sections) != 0:
        raise RuntimeError("There probably shouldn't be any sections!")

    # Note: Don't do the following!!
    # This often causes a segmentation fault because, I believe, it uses a
    #    cached ODB that is out of date and out of our control. Explicitly 
    #    opening and accessing the ODB is preferable.
    # new_sections = mdb.models[model_name].sectionsFromOdb(path_to_odb)

    odb = odbAccess.openOdb(path_to_odb)
    sections = odb.models[odb_model_name].sections

    assert(len(sections) == 1)

    section = sections[sections.keys()[0]]
    _check_section_properties(section)

    mdb.models[new_model_name].HomogeneousSolidSection(STANDARD_SECTION_NAME, STANDARD_MATERIAL_NAME)

    odb.close()



def _check_section_properties(section: Any) -> None:
    """DEPRECATED. Was a helper for the create_section_from_odb() function.
    
       Checks the properties of a section. 

       Args:

       Returns:
           None.

       Raises:
           RuntimeError: Section doesn't have the expected properties.
    """

    if section.name != STANDARD_SECTION_NAME:
        raise RuntimeError("The section doesn't have the normal name!")

    if section.material != STANDARD_MATERIAL_NAME:
        raise RuntimeError("The material of the section doesn't have the normal name!")



def _partition_face_with_sketch(face: Any, edge: Any, sketch: Any, assembly: Optional[Any] = None, 
                                instance: Optional[Any] = None, part: Optional[Any] = None
                               ) -> Any:
    """DEPRECATED. Now, instead of partitioning a face with a sketch, a face
           is partitioned by creating edges on it which define the new face.
    
       Partitions the face of an object with a sktech.

       The PartitionFaceBySketch() method claims to take a "sequence of Face
           objects". I believe it actually takes an Abaqus FaceArray object. 

       Args:
           face:     Abaqus Face object. The face to partition.
           edge:     Abaqus Edge object. Specifies the orientation of the sketch 
                         on the face. This edge should have been used to orient 
                         the sketch on the face when the sketch was constructed.
           sketch:   Abaqus ConstrainedSketch object. The sketch with which to 
                         do the partitioning.
           assembly: Abaqus Assembly object.
           instance: Abaqus PartInstance object.
           part:     Abaqus Part object.

           Note that either the assembly and instance arguments must be passed
               or part argument must be passed.

       Returns:
           Abaqus Face object. The face created by doing the partitioning.

       Raises:
           AssertionError: The partitioning did not produce a new face.
    """

    point = face.pointOn

    if assembly is not None and instance is not None:
        face_array = instance.faces.findAt(coordinates=point)
        face_cnt = len(instance.faces)
        feature = assembly.PartitionFaceBySketch(face_array, sketch, sketchUpEdge=edge)
        post_partition_face_cnt = len(instance.faces)
    elif part is not None:
        face_array = part.faces.findAt(coordinates=point)
        face_cnt = len(part.faces)
        feature = part.PartitionFaceBySketch(face_array, sketch, sketchUpEdge=edge)
        post_partition_face_cnt = len(part.faces)
    else:
        raise RuntimeError("Invalid combination of arguments passed.")

    if (post_partition_face_cnt != face_cnt + 1) or (feature is None):
        raise AssertionError("Partitioning didn't work!")

    # Hack! Assumes that the created face is the last face in the array of faces. 
    if assembly is not None and instance is not None:
        new_face = instance.faces[-1]
    elif part is not None:
        new_face = part.faces[-1]

    return new_face



def _sketch_ngon(ngon: geom.NGon3D, transform: Any, sketch: Any) -> Any:
    """DEPRECATED. Was used as a helper to partition a face.
    
       Sketches the points which make up an ngon on a sketch which has some
           orientation in space.

       Args:
           ngon:      The ngon to sketch.
           transform: Abaqus Transform object. Describes the orientation of
                          the sketch in space.
           sketch:    Abaqus ConstrainedSketch object. The empty sketch

       Returns:
           Abaqus ConstrainedSketch object. The sketch wich has the ngon on
               it.

       Raises:
           AssertionError: Some point is a nonzero component in the axis normal
                               to the sketch.
    """

    # Extract the coordinate system of the sketch. 
    csys = _extract_global_csys_to_sketch_csys(transform)

    # Map the points into the coordinate system of the sketch. 
    points = csys.map_into_csys(ngon.vertices)

    # In the coordinate system of the sketch, it better be the case that the 
    #    component of each point normal to the face is zero. 
    for point in points:
        if not geom.float_equals(point.rep[2], 0):
            raise AssertionError("After mapping the points into the coordinate" + 
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



def _get_any_edge_on_face(face: Any, obj: Any) -> Any:
    """DEPRECATED. Was helper for partition face technique.
    
       Gets an arbitrary edge on a face.

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



def pick_nodes_randomly(cnt: int, path_to_odb: str) -> list[geom.Point3D]:
    """DEPRECATED. Was used to pick nodes at random in the process of recovering
           linear relationships.

       Picks some nodes at random from a mesh via the content of an ODB.

       Assumes that the ODB contains a single non-initial step. The nodes are 
           picked from nodes in the last frame of the step. 
        
       Args:
           cnt:         The number of nodes to pick.
           path_to_odb: Absolute path to .odb file.
    
       Returns:
           A list of the nodes chosen at random.
    
       Raises:
           None.
    """

    odb = odbAccess.openOdb(path=path_to_odb, readOnly=True)

    if len(odb.steps) != 1:
        raise RuntimeError("The ODB does not contain exactly one step!") 

    name = odb.steps.keys()[0]
    
    # Figure out which frame is last.
    step = odb.steps[name]
    largest_increment_number = step.frames[0].incrementNumber
    for frame in step.frames: 
        if frame.incrementNumber > largest_increment_number:
            largest_increment_number = frame.incrementNumber
    last_frame = step.frames[largest_increment_number]

    assert "U" in last_frame.fieldOutputs, "No displacement measurements in frame!"
    last_frame_displacements = last_frame.fieldOutputs["U"]
    
    # For the displacement field, there is one value at each node in the mesh.
    node_cnt = len(last_frame_displacements.values)
    indices = random.choices(range(node_cnt), k=cnt)

    ret = []
    for idx in indices:
        # This is weird, but I'm pretty sure it's necessary.
        # An alternative approach would be to use the .inp file generated by
        #     the job. Using the ODB seems easier.
        nodes = last_frame_displacements.values[idx].instance.nodes
        ret.append(geom.Point3D(np.array(nodes[idx].coordinates))) 

    odb.close()

    return ret



def read_displacements(points: list[geom.Point3D], path_to_odb: str
                      ) -> dict[geom.Point3D, tuple[geom.Point3D, geom.Vec3D]]:
    """DEPRECATED. Was used to read displacements of nodes in the process of
           using linear relationships to recover good tractions to apply.

       Reads an ODB to determine the displacements of some points.

       Assumes that the ODB contains a single non-initial step and a single
           meshed part instance. The displacements are computed for this single
           step.
        
       Args:
           points:      The points of interest. If these points do not match
                            nodes in the mesh of the single part instance, then 
                            only approximate displacements can be found. See 
                            the caveat in the return.
           path_to_odb: Absolute path to a .odb file.
    
       Returns:
           Dictionary which maps each point to a node and a displacement vector.     
               The node is the node in the mesh closest to the point and the
               displacement vector is the displacement of that node.

       Raises:
           None.
    """

    odb = odbAccess.openOdb(path=path_to_odb, readOnly=True)
    
    if len(odb.steps) != 1:
        raise RuntimeError("The ODB does not contain exactly one step!") 

    name = odb.steps.keys()[0]

    # The step potentially contains many frames.
    # Since we care about the displacement of the points across the whole step,
    #     it's necessary to find the first and last frames in the step.
    step = odb.steps[name]
    smallest_increment_number = step.frames[0].incrementNumber
    largest_increment_number = step.frames[0].incrementNumber
    for frame in step.frames: 
        if frame.incrementNumber < smallest_increment_number:
            smallest_increment_number = frame.incrementNumber
        elif frame.incrementNumber > largest_increment_number:
            largest_increment_number = frame.incrementNumber
    first_frame = step.frames[smallest_increment_number]
    last_frame = step.frames[largest_increment_number]

    assert "U" in first_frame.fieldOutputs, "No displacement measurements in frame!"
    assert "U" in last_frame.fieldOutputs, "No displacement measurements in frame!"
    
    # Technically there is a displacement in the first frame of a step.
    #     This is taken into account.
    first_frame_displacements = first_frame.fieldOutputs["U"]
    last_frame_displacements = last_frame.fieldOutputs["U"]

    # Since displacements are only computed at nodal points, figure out which
    #     nodal point each point is closest to.

    # Start by computing the distance from every point to the first node.
    assert len(odb.rootAssembly.instances) == 1, "Not exactly one part instance in the assembly!"
    instance_name = odb.rootAssembly.instances.keys()[0]
    nodes = odb.rootAssembly.instances[instance_name].nodes
    first_node = geom.Point3D(np.array(nodes[0].coordinates))
    closest_nodes = []
    for point in points:
        # Record the index of the closest node and the distance from the closest
        #     node to the point.
        closest_nodes.append((1, geom.distance(first_node, point)))

    # Then, for each node, compare it to the closest node for each point.
    #     If it is closer, record as such.
    for node in nodes:
        node_id = node.label
        node_point = geom.Point3D(np.array(node.coordinates))
        for idx, point in enumerate(points):
            dist = geom.distance(node_point, point)
            if dist < closest_nodes[idx][1]:
                # This node is closer to the point! Record as such!
                closest_nodes[idx] = (node_id, dist)
                 
    ret = {}
    for idx, point in enumerate(points):
        closest_node_id = closest_nodes[idx][0]
        
        # The node closest to the point.
        # Careful! The node ids go from 1 -> N, so the 0th entry in the list
        #     of nodes corresponds to node with id 1.
        closest_node = geom.Point3D(np.array(nodes[closest_node_id - 1].coordinates))

        # Assumes the field displacement ouputs are ordered according to the 
        #     ids of the nodes.
        init_disp = geom.Vec3D(np.array(first_frame_displacements.values[closest_node_id - 1].data))
        final_disp = geom.Vec3D(np.array(last_frame_displacements.values[closest_node_id - 1].data))
        disp = final_disp - init_disp

        ret[point] = (closest_node, disp)
         
    odb.close()

    return ret 
