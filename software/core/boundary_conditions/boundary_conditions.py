import util.geom as geom
import abaqus.abaqus_shim as shim



# TODO: Build this out if necessary.
class BCSettings:

    def __init__(self):
        raise RuntimeError("Not yet available.")



# A displacement boundary condition is defined by restrictions on displacement
#    and rotation in the three spatial directions.
class DisplacementBCSettings(BCSettings):
    
    def __init__(self, fix_x, fix_y, fix_z, prevent_x_rot, prevent_y_rot, prevent_z_rot):
    # type: (bool, bool, bool, bool, bool, bool) -> None

        if fix_x == True and fix_y == True and fix_z == True and prevent_x_rot \
           == True and prevent_y_rot == True and prevent_z_rot == True:
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

    # TODO: For now, we assume that all BCs exist on surfaces in space. 
    def __init__(self, ngon, bc_settings):
    # type: (List[geom.NGon3D], List[BCSettings]) -> None
    
        self.surface = ngon
        self.bc_settings = bc_settings
        


# Apply some surface boundary conditions to a part or part instance at a
#    particular step.
def apply_surface_BCs(BCs, part_rep, step_name, model_name, mdb):
# type: (List[BC], Any, str) -> None

    for BC in BCs:
        
        # Add the reference points. 
        ref_points = shim.add_ref_points(BC.surface.vertices, part_instance)

        # Use the reference points to build a region.
        region = build_region_with_ref_points(ref_points) 
        
        # Construct name for the new boundary condition.
        BC_cnt = shim.get_BC_cnt(model_name, mdb)
        bc_name = shim.STANDARD_BC_PREFIX + str(BC_cnt)

        # And apply the boundary condition to the region.
        shim.create_bc_for_region(BC, region, step_name, bc_name, model_name, mdb)




