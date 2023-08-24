import util.geom as geom
import core.abaqus.abaqus_shim as shim
import numpy as np

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
        


# Apply some boundary conditions to a part in a model. 
def apply_BCs(BCs, step_name, part, model_name, mdb):
# type: (List[BC], str, Any, Any) -> None

    for BC in BCs:

        BC_cnt = shim.get_BC_cnt(model_name, mdb)
        BC_name = shim.STANDARD_BC_PREFIX + str(BC_cnt)

        face_feature = partition_face(BC.region, BC_name, part, model_name, mdb) 
        region = shim.build_region_with_face(face_feature, part)

        if isinstance(BC.bc_settings, DisplacementBCSettings):
            # TODO: This is horrific.
            shim.create_displacement_bc(BC_name, step_name, region, BC.bc_settings.fix_x, BC.bc_settings.fix_y, BC.bc_settings.fix_z, BC.bc_settings.prevent_x_rot, BC.bc_settings.prevent_y_rot, BC.bc_settings.prevent_z_rot, model_name, mdb)
        else:
            raise RuntimeError("An attempt was made to add an unsupported type of boundary condition to a region!")



# Given an ngon, this function finds the face on the part that the ngon
#    belongs to. It then partitions that face so the ngon acts as a standalone
#    face.
# The face that the ngon belongs to is determined via a two step process. First,
#    an attempt is made to find a face that is represented by the same plane as
#    the ngon. If such a face exists, it may not be the one the ngon actually
#    lies within. To determine if the ngon lies entirely within the face, a
#    check is done to see if each vertex which composes the ngon lies within,
#    or on the boundary of, the face. If all vertices of the ngon lie on or
#    within the face, then the face is partitioned according to the geometry of
#    the ngon.
# This function returns the Feature object associated with the new face.
def partition_face(ngon, new_face_name, part, model_name, mdb):
# type: (geom.NGon3D, str, Any, Any, Any) -> Any

    """
    Step 1: Find the face of the part that the ngon belongs to.
    """

    # Get the vertices associated with all the faces.
    vertices = shim.get_part_vertices(part)

    face_ngon_belongs_to = None
    
    # Find the face that the ngon lives on.
    for face_idx, vertices_single_face in enumerate(vertices):
        
        # TODO: This is awkward...
        vertices_single_face = [geom.Point3D(*vertex[0]) for vertex in vertices_single_face]

        # Subtle point: If the points which define the vertices of a face are
        #    not on a plane (i.e. there is a curved surface), then this will
        #    almost always fail. This is intended behavior. No functionality for
        #    partitioning a curved face is desired.
        # TODO: To ensure that the face which is found is not a curved surface,
        #    the getCurvature() method could be checked at a couple points.
        if geom.points_on_plane_of_ngon(vertices_single_face, ngon):
            if geom.points_in_ngon_3D(ngon.vertices, geom.NGon3D(vertices_single_face)):
                face_ngon_belongs_to = shim.get_faces(part)[face_idx]
                break

    if face_ngon_belongs_to == None:
        raise RuntimeError("Could not identify which face the ngon belongs to. \
                            Can't continue to do partitioning.")
    
    # DEBUG
    dp("The vertices associated with the face the ngon lives on are " + str(shim.get_face_vertices(face_ngon_belongs_to, part)))


    """
    Step 2: Construct and orient a sketch which lives on that face. Record
               /deduce the orientation of the sketch for future use.
    """
    
    # When the sketch is built, it's sufficient to specify which face that the
    #    sketch should live on and how the sketch should be oriented on that 
    #    face. However, it's also necessary to record geometric information about
    #    the face and the orientation because it will be necessary to translate
    #    points in the global coordinate system into the coordinate system of
    #    the sketch.

    # ----- Building and orienting the sketch in Abaqus -----

    # Pick any edge belonging to the face.
    edge = shim.get_any_edge_on_face(face_ngon_belongs_to, part) 

    # DEBUG
    dp("The vertices which belong to the edge are " + str(shim.get_edge_vertices(edge, part)))

    face_normal = geom.Vec3D(*shim.get_face_normal(face_ngon_belongs_to))
    face_centroid = geom.Point3D(*shim.get_face_centroid(face_ngon_belongs_to))

    # DEBUG
    dp("The normal to the face the ngon lives on is " + str(face_normal.np_arr))
    dp("The centroid of the face the ngon lives on is " + str(face_centroid.get_components()))

    transform = shim.build_sketch_transform(face_ngon_belongs_to, face_centroid, edge, part)

    # TODO: Bad magic number. I don't think that the sketch size actually matters
    #    at all.
    sketch = shim.build_constrained_sketch(transform, new_face_name, 1000, model_name, mdb)

    # ----- Deducing and recording orientation of sketch for change of csys -----

    # A point and a normal define a plane, but it's also necessary to define 
    #    a coordinate system on the plane. This means that it's necessary 
    #    to define which direction on the plane is the positive x axis and 
    #    which direction is the positive y axis. Obviously, a point needs to
    #    be picked as the origin of the coordinate system.
    # When the constrained sketch is built, the arguments passed for
    #    sketchOrientation and sketchUpEdge define the directions of the
    #    positive x axis and the positive y axis. 
    # If an Edge, as opposed to a DatumAxis, is passed to sketchUpEdge 
    #    and Right is specified for sketchOrientation, we cannot immediately 
    #    determine what the directions of the postive y axis and positive x 
    #    axis are without additional information. They depend on the relationship 
    #    between the edge and the face in a non-trivial way. 
    # Passing Right or Left, as opposed to Up or Down, as sketchOrientation 
    #    means that the edge defines the y axis on the sketch. An edge is
    #    composed of two vertices. Two vectors can be constructed with two
    #    vertices (by swapping the order of subtraction). One of the vectors 
    #    points in the positive y direction and the other points in the negative
    #    y direction. To determine which is which, we can do the following:
    #       1) Construct the two vectors defined by the edge's vertices.
    #       2) Construct the vector normal to the face.
    #       3) Compute the cross product between each vector from (1) and the
    #          vector from (2). The result is two vectors which lie in the
    #          plane of the face.
    #       4) Check which vector from (3), when it emanates as a ray from
    #          the midpoint of the edge, intersects the face. Exactly one of
    #          the two vectors from (3) will intersect the face. The
    #          corresponding vector between the edge's vertices points in the
    #          negative y direction. This conclusion assumes that Right is used
    #          as the sketchOrientation.

    # Get the two vertices for this edge.
    # Construct the two vectors in (1).
    vertices = shim.get_edge_vertices(edge, part)  
    vec_1 = geom.Vec3D(*(np.array(vertices[0]) - np.array(vertices[1])))
    vec_2 = geom.Vec3D(*(np.array(vertices[1]) - np.array(vertices[0])))

    # Compute the cross product between each vector from (1) and the vector
    #    normal to the face. The resulting vectors are those produced in (3).
    vec_1_in_plane_face = geom.find_third_orthonormal(vec_1, face_normal)
    vec_2_in_plane_face = geom.find_third_orthonormal(vec_2, face_normal)

    # Compute the midpoint of the edge.
    midpoint = geom.Vec3D(*((np.array(vertices[0]) + np.array(vertices[1]))/2)) 

    # Construct points very near the midpoint in the two potential directions.
    point_1 = geom.Point3D(*(.001 * vec_1_in_plane_face.np_arr + midpoint.np_arr))
    point_2 = geom.Point3D(*(.001 * vec_2_in_plane_face.np_arr + midpoint.np_arr))

    # Check which point lies in the ngon defined by the face.
    # TODO: Write a function/find function which can check if a ray and a 
    #    polygon in 3D intersect. This would generalize the heuristic approach.
    # If the point lies in the ngon defined by the face, the underlying vector 
    #    from (1) points in the direction of the negative y axis. If it does not
    #    , it points in the direction of the positive y axis.
    face_vertices = shim.get_face_vertices(face_ngon_belongs_to, part)
    face_vertices = [geom.Point3D(*vertex) for vertex in face_vertices] # TODO: Awkward!
    face_ngon = geom.NGon3D(face_vertices)
    if geom.point_in_ngon_3D(point_1, face_ngon):
        positive_y_axis = vec_2.get_norm()
    elif geom.point_in_ngon_3D(point_2, face_ngon):
        positive_y_axis = vec_1.get_norm()
    else:
        raise RuntimeError("Nether point on either side of the edge lies on the face.")

    positive_x_axis = geom.find_third_orthonormal(positive_y_axis, face_normal)  
   
    # DEBUG
    dp("The positive x, y, and z axis are in the directions given by: " + str(positive_x_axis.np_arr) + ", " + str(positive_y_axis.np_arr) + ", " + str(face_normal.np_arr))

    """
    Step 3: Draw the ngon on the sketch.
    """

    # TODO: An orientation, defined by the directions of a positive x axis, a
    #    positive y axis, and a positive z axis, should be a user defined type.
    #    An instance of this type should be passed into change_csys.
    # First, map the vertices of the ngon into the coordinate system of the
    #    sketch. 
    points = geom.change_csys(ngon.vertices, face_centroid, positive_x_axis, positive_y_axis, face_normal)

    # DEBUG
    dp("The coordinates of the points in the new coordinate system are: " + str([point.get_components() for point in points]))
    
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

        # Documentation is potentially wrong. Even on success, the Line() method
        #    of a sketch returns None!
        # if sketch.Line(cur_point, next_point) == None:
        #     raise RuntimeError("Failed to create line for polygon!")
        res = sketch.Line(cur_point, next_point)

        # DEBUG
        dp("Attempted to create line between " + str(cur_point) + " and " + str(next_point))
        dp("The result was " + str(res))

    # DEBUG
    shim.save_mdb_as("/home/andlars/Desktop/RS_TP_Adaptation/software/script_testing/test_initial_geometry/test_initial_geometry_intermed.cae", mdb)

    """
    Step 4: Partition the face using the sketch. 
    """
    return shim.partition_face_with_sketch(face_ngon_belongs_to, sketch, part)
    
