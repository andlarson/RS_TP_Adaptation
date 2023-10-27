import numpy as np

import util.geom as geom
import core.abaqus.abaqus_shim as shim
from util.debug import *



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



class BC:

    def __init__(self):
    # type: (None) -> None

        raise RuntimeError("Not supported.")



# A boundary condition across some surface.
class SurfaceBC(BC):
        
    def __init__(self, ngon, bc_settings):
    # type: (geom.NGon3D, BCSettings) -> None
    
        self.region = ngon
        self.bc_settings = bc_settings



# A boundary condition on some set of points.
class VertexBC(BC):

    def __init__(self, vertices, bc_settings):
    # type: (List[geom.Point3D], BCSettings) -> None

        self.region = vertices
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
#    model_name    - String.
#    mdb           - Abaqus MDB object.
#
# Notes:
#    In order to apply boundary conditions to regions which don't already exist, 
#       it's necessary to partition pre-existing features. In general, an Abaqus 
#       Part object or an Abaqus PartInstance object can be partitioned. This 
#       function creates new regions on parts, when necessary, by partitioning 
#       this single Abaqus PartInstance object.
#    Since a PartInstance object must be passed, it means that the Assembly is 
#       not empty. In my experience, if you try to create a boundary condition 
#       on a region of an Abaqus Part object and the Assembly is empty, Abaqus 
#       segfaults.
#
# Returns:
#    None.
def apply_BCs(BCs, step_name, part_instance, model_name, mdb):
# type: (List[BC], str, Any, str, Any) -> None

    for BC in BCs:

        BC_cnt = shim.get_BC_cnt(model_name, mdb)
        BC_name = shim.STANDARD_BC_PREFIX + str(BC_cnt)

        if isinstance(BC, SurfaceBC):
            assembly = mdb.models[model_name].rootAssembly
            face = shim.partition_face(BC.region, instance=part_instance, assembly=assembly) 
            region = shim.build_region_with_face(face, part_instance)

        elif isinstance(BC, VertexBC):
            raise RuntimeError("Not yet supported.")

        else:
            raise RuntimeError("Not yet supported.")

        if isinstance(BC.bc_settings, DisplacementBCSettings):
            shim.create_displacement_bc(BC_name, step_name, region, BC.bc_settings, model_name, mdb)

        else:
            raise RuntimeError("Not yet supported.")