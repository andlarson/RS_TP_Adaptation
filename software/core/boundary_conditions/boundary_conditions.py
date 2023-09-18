import numpy as np

import util.geom as geom
import core.abaqus.abaqus_shim as shim
from util.debug import *



# TODO: Build this out if necessary.
class BCSettings:

    def __init__(self):
        raise RuntimeError("Not yet available.")



# A displacement boundary condition is defined by restrictions on displacement
#    and rotation in the three spatial directions.
class DisplacementBCSettings(BCSettings):
    
    def __init__(self, fix_x, fix_y, fix_z, prevent_x_rot, prevent_y_rot, prevent_z_rot):
    # type: (bool, bool, bool, bool, bool, bool) -> None

        if fix_x == False and fix_y == False and fix_z == False and prevent_x_rot \
           == False and prevent_y_rot == False and prevent_z_rot == False:
            raise RuntimeError("The boundary condition settings don't restrict anything!")

        self.fix_x = fix_x
        self.fix_y = fix_y
        self.fix_z = fix_z
        self.prevent_x_rot = prevent_x_rot
        self.prevent_y_rot = prevent_y_rot
        self.prevent_z_rot = prevent_z_rot



# A boundary condtion is essentially composed of a region of space and some settings.
# For example, there might be some boundary condtions on a set of points.
# Alternatively, there might a boundary condition on some surface.
class BC:

    def __init__(self, ngon, bc_settings):
    # type: (geom.NGon3D, BCSettings) -> None
    
        self.region = ngon
        self.bc_settings = bc_settings
        


# Apply boundary conditions to an Assembly.
# Note that boundary conditions can only be applied to an Assembly, not to a
#    Part.  
#
# Arguments:
#    BCs           - List of BC objects. 
#    step_name     - String. 
#                    The step at which the boundary conditions are first applied. 
#                       By default, they propagate to future steps. 
#    part_instance - Abaqus PartInstance object. 
#                    In order to apply boundary conditions to regions which don't
#                       already exist, it's necessary to partition pre-existing
#                       features. In general, an Abaqus Part object or an Abaqus 
#                       PartInstance object can be partitioned. This function
#                       creates new regions on parts, when necessary, by
#                       partitioning this single Abaqus PartInstance object.
#                    Since a PartInstance object must be passed, it means that
#                       the Assembly is not empty. In my experience, if you try
#                       to create a boundary condition on a region of an Abaqus
#                       Part object and the Assembly is empty, Abaqus segfaults.
#    model_name    - String.
#    mdb           - Abaqus MDB object.
#
# Returns:
#    None.
def apply_BCs(BCs, step_name, part_instance, model_name, mdb):
# type: (List[BC], str, Any, str, Any) -> None

    for BC in BCs:

        BC_cnt = shim.get_BC_cnt(model_name, mdb)
        BC_name = shim.STANDARD_BC_PREFIX + str(BC_cnt)

        face_feature = partition_face(BC.region, BC_name, model_name, mdb, instance=part_instance) 
        region = shim.build_region_with_face(face_feature, part_instance)

        if isinstance(BC.bc_settings, DisplacementBCSettings):
            
            # TODO: This is horrific.
            shim.create_displacement_bc(BC_name, step_name, region, BC.bc_settings.fix_x, BC.bc_settings.fix_y, BC.bc_settings.fix_z, BC.bc_settings.prevent_x_rot, BC.bc_settings.prevent_y_rot, BC.bc_settings.prevent_z_rot, model_name, mdb)

        else:
            raise RuntimeError("An attempt was made to add an unsupported type of boundary condition to a region!")



# Partitions the face of an object.
# 
# Notes:
#    Uses distinct part and instance arguments to avoid type checking Abaqus
#       object types. It's important to know the type of what is being partitioned
#       because the PartitionFaceBySketch() method is valid for Abaqus Part and
#       Abaqus Assembly objects.
#    Even if the ngon happens to match the geometry of a face which already
#       exists, this function still returns a new logical face. 
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

    if part == None and instance == None:
        raise RuntimeError("Missing required argument!")
    elif part != None:
        obj = part
    else:
        obj = instance

    assembly = mdb.models[model_name].rootAssembly

    """
    Step 1: Find the face that the ngon belongs to.
    """
  
    # Critical to ensure that the vertices are well-ordered on a per-face basis.
    vertices = shim.get_all_vertices_ordered(obj)

    face_ngon_belongs_to = None

    # Find the face that the ngon lives on.
    for idx, vertices_single_face in enumerate(vertices):

        # Subtle point: If the points which define the vertices of a face are
        #    not on a plane (i.e. there is a curved surface), then this will
        #    almost always fail. This is intended behavior. No functionality for
        #    partitioning a curved face is desired.
        # TODO: Account for situation where the ngon lives inside two nested
        #    faces. Right now, the first face that the ngon lives in is the one
        #    which is identified and used.
        # TODO: Use the getCurvature() method to check that faces with curvature
        #    are not considered.
        # TODO: It turns out that using idx works because the ordering of
        #    faces is the same. This could be formalized in some way.
        if geom.points_on_plane_of_ngon(vertices_single_face, ngon):
            if geom.points_in_ngon_3D(ngon.vertices, geom.NGon3D(vertices_single_face)):
                face_ngon_belongs_to = shim.get_faces(obj)[idx]
                break

    if face_ngon_belongs_to == None:
        raise RuntimeError("Could not identify which face the ngon belongs to. \
                            Can't continue to do partitioning.")
   
    # If the face the ngon belongs to has vertices which exactly match the
    #    ngon, we don't want to try to partition the face and create a new face.
    #    This will cause partitioning to fail and Abaqus to give up. 
    # TODO: The ngon could have redundant vertices, but such a partitioning will
    #    actually succeed because new vertices will be added to the face.
    # TODO: Assuming that the vertices which make up the Abaqus Face object are
    #    not redundant.
    face_vertices = shim.get_face_vertices(face_ngon_belongs_to, obj)
    ngon_vertices = ngon.get_builtin_rep()
    ngon_face_share_vertices = True 
    if len(face_vertices) == len(ngon_vertices):
        for face_vertex in face_vertices:
            if face_vertex not in ngon_vertices:
                ngon_face_share_vertices = False

    # If they do share vertices, no partitioning needs to be done at all.
    if ngon_face_share_vertices:
        return face_ngon_belongs_to

    """
    Step 2: Construct and orient a sketch which lives on that face.
    """

    edge = shim.get_any_edge_on_face(face_ngon_belongs_to, obj)
    transform = shim.build_sketch_transform_from_face(face_ngon_belongs_to, edge, assembly)

    sketch = shim.build_constrained_sketch(transform, new_face_name, model_name, mdb)


    """
    Step 3: Draw the ngon on the sketch.
    """

    # Extract the coordinate system of the sketch. 
    csys = shim.extract_global_csys_to_sketch_csys(transform)

    # Map the points into the coordinate system of the sketch. 
    points = csys.from_global_to_new(ngon.vertices)

    # In the coordinate system of the sketch, it better be the case that the 
    #    component of each point normal to the face is zero. 
    for point in points:
        if not geom.float_equals(point.z, 0):
            raise RuntimeError("After mapping the points into the coordinate" + 
                               " system of the sketch, a point had a nonzero" +
                               " component in the axis normal to the sketch!")

    # Now use the points to construct the ngon on the sketch.
    points_circular = points + points[0:1]
    for idx in range(len(points_circular) - 1):
        cur_point = points_circular[idx].proj_xy().get_components()
        next_point = points_circular[idx + 1].proj_xy().get_components()
        
        # The Line() method of an Abaqus Sketch object returns None even on
        #    success. Documentation is wrong!
        sketch.Line(cur_point, next_point)


    """
    Step 4: Partition the face using the sketch. 
    """

    if part != None: 
        shim.partition_face_with_sketch(face_ngon_belongs_to, edge, sketch, part=part)
    else:
        shim.partition_face_with_sketch(face_ngon_belongs_to, edge, sketch, assembly=assembly, instance=instance)


    """
    Step 5: Find the face which was just created.
    """
    
    # Find the face created by the ngon. 
    vertices = shim.get_all_vertices(obj)
    new_face = None
    for idx, vertices_single_face in enumerate(vertices):
       
        # TODO: Awkward conversion necessary...
        vertices_single_face = [v.get_components() for v in vertices_single_face]

        # TODO: It turns out that using idx works because the ordering of
        #    faces is the same. This could be formalized in some way.
        if len(vertices_single_face) == len(ngon.vertices):

            all_match = True
            for vertex in ngon.vertices:
                if vertex.get_components() not in vertices_single_face:
                    all_match = False

            if all_match:
                new_face = shim.get_faces(obj)[idx]
                break

    if new_face == None:
        raise RuntimeError("Failed to find the face which was just created!")

    return new_face
