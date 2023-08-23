import util.geom as geom
import core.abaqus.abaqus_shim as shim



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
        


# Apply some boundary conditions to the part. 
def apply_BCs(BCs, step_name, part):
# type: (List[BC], str, Any) -> None

    for BC in BCs:

        BC_cnt = shim.get_BC_cnt(part)
        BC_name = shim.STANDARD_BC_PREFIX + str(BC_cnt)

        face_feature = partition_face(BC.region, BC_name, part) 
        region = shim.build_region_with_face(face_feature, part)

        # And apply the boundary condition to the region.
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
def partition_face(ngon, new_face_name, part):
# type: (geom.NGon3D, str, Any) -> Any

    """
    Step 1: Find the face of the part that the ngon belongs to.
    """

    # Get the vertices associated with all the faces.
    vertices = shim.get_vertices(part)

    face_ngon_belongs_to = None
    
    # Find the face that the ngon lives on.
    for face_idx, vertices_single_face in enumerate(vertices):
        
        # TODO: This is awkward...
        vertices_single_face = [geom.Point3D(vertex[0], vertex[1], vertex[2]) for vertex in vertices_single_face]

        # Subtle point: If the points which define the vertices of a face are
        #    not on a plane (i.e. there is a curved surface), then this will
        #    almost always fail. This is intended behavior. No functionality for
        #    partitioning a curved face is desired.
        # TODO: To ensure that the face which is found is not a curved surface,
        #    the getCurvature() method could be checked at a couple points.
        if geom.points_on_plane_of_ngon(vertices_single_face, ngon):
            if geom.points_in_ngon(ngon.vertices, geom.NGon3D(vertices_single_face)):
                face_ngon_belongs_to = shim.get_vertices(part)[face_idx]

    if face_ngon_belongs_to == None:
        raise RuntimeError("Could not identify which face the ngon belongs to. \
                            Can't continue to do partitioning.")
   
    """
    Step 2: Construct a sketch which lives on that face.
    """

    # First, define the plane the sketch lives on in space. Since the sketch 
    #    should live on the target face, its origin should be at the centroid 
    #    of the face and its normal should match the normal of the face. 
    face_centroid = geom.Point3D(*shim.get_face_centroid(face))
    face_normal = geom.Vec3D(*shim.get_face_normal(face))

    # A point and a normal define a plane, but it's also necessary to define 
    #    how drawing will be done on the plane. This means that it's necessary 
    #    to define which direction on the plane is the positive x axis and 
    #    which direction is the positive y axis.
    # These axes are chosen in arbitrary directions. They must be orthogonal
    #    to one another and lie on the plane.
    first_axis = face_normal.get_orthonormal()
    second_axis = geom.find_third_orthonormal(first, face_normal)

    # Put all the axes together to specify the orientation of the sketch in
    #    3D.
    orientation = np.stack((first_axis, second_axis, face_normal), axis=0)
   
    # TODO: Bad magic number. I don't think that the sketch size actually matters
    #    at all.
    sketch = shim.build_constrained_sketch(orientation, face_centroid, new_face_name, 1, part)

    
    """
    Step 3: Draw the ngon on the sketch.
    """

    # First, map the vertices of the ngon into the coordinate system of the
    #    sketch. 
    points = geom.change_csys(ngon.vertices, face_centroid, first_axis, second_axis, face_normal)
    
    # In the coordinate system of the sketch, it better be the case that the 
    #    component of each point normal to the face is zero. 
    for point in points:
        if not geom.float_equals(point.z, 0):
            raise RuntimeError("After mapping the points into the coordinate" + 
                               " system of the sketch, a point had a nonzero" +
                               " component in the axis normal to the sketch!")

    # Now use the points to construct the ngon on the sketch.
    point_cnt = len(points)
    points_circular = points + points[0]
    for idx in range(point_cnt):
        cur_point = points_circular[idx]
        next_point = points_circular[idx + 1]
        if sketch.Line(cur_point, next_point) == None:
            raise RuntimeError("Failed to create line for polygon!")


    """
    Step 4: Partition the face using the sketch. 
    """

    return shim.partition_face_with_sketch(face, sketch, part)
    
