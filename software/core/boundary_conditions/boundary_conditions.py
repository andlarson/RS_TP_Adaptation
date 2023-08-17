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
    
        self.surface = ngon
        self.bc_settings = bc_settings
        


# Apply some surface boundary conditions to the root assembly at a particular 
#    step.
def apply_surface_BCs(BCs, step_name, model_name, mdb):
# type: (List[BC], str, str, Any) -> None

    for BC in BCs:
        
        ref_point_features = shim.add_ref_points(BC.flat_surface.vertices, model_name, mdb)
        region = shim.build_region_with_ref_points(ref_point_features, model_name, mdb) 
        
        # Construct name for the new boundary condition.
        BC_cnt = shim.get_BC_cnt(model_name, mdb)
        BC_name = shim.STANDARD_BC_PREFIX + str(BC_cnt)

        # And apply the boundary condition to the region.
        if isinstance(BC.bc_settings, DisplacementBCSettings):
            shim.create_displacement_bc(BC_name, step_name, region, BC.bc_settings.fix_x, BC.bc_settings.fix_y, BC.bc_settings.fix_z, BC.bc_settings.prevent_x_rot, BC.bc_settings.prevent_y_rot, BC.bc_settings.prevent_z_rot, model_name, mdb)
        else:
            raise RuntimeError("An attempt was made to add an unsupported type of boundary condition to a region!")



# Given an ngon, this function finds the face on the assembly that the ngon
#    belongs to. It then partitions that face so the ngon acts as a standalone
#    face.
# The face that the ngon belongs to is determined by checking if the ngon 
#    belongs to the same plane as the vertices of the face. This does not
#    guarantee that the ngon lies entirely within the face. We assume that
#    it does. If this is not the case, an exception is thrown when the
#    partitioning inevitably fails.
def partition_assembly_face_from_ngon(ngon, model_name, mdb):

    # Get the vertices associated with all the faces on the Assembly.
    vertices = shim.get_vertices_of_faces_on_assembly(model_name, mdb)

    # Find the face that the ngon lives on.
    for vertices_single_face in vertices:
        
        # TODO: This is an awkward type conversion for posterity...
        vertices_single_face = [geom.Point3D(vertex[0], vertex[1], vertex[2]) for vertex in vertices_single_face]

        if geom.points_on_plane_of_ngon(ngon.):
             



