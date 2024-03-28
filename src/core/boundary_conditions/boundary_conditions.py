"""
This module offers classes and functions related to boundary conditions.
"""

from typing import Any, TypeVar
from abc import ABCMeta, abstractmethod

import numpy as np

import src.util.geom as geom
import src.core.abaqus.abaqus_shim as shim

from src.util.debug import *



class DisplacementBCSettings():
    
    def __init__(self, fix_x: bool, fix_y: bool, fix_z: bool, prevent_x_rot: bool, 
                 prevent_y_rot: bool, prevent_z_rot: bool) -> None:
        """Creates a displacement boundary condition.
           
           The displacment boundary condition is not associated with an object
               or surface.

           Args:
               fix_x:         Restricts displacement in the x direction.
               fix_y:         Restricts displacement in the y direction.
               fix_z:         Restricts displacement in the z direction.
               prevent_x_rot: Restricts rotation about the x axis.
               prevent_y_rot: Restricts rotation about the y axis.
               prevent_z_rot: Restricts rotation about the z axis.

               Note that at least one argument must be true or nothing is restricted
                   and there is no need for a boundary condition.

           Returns:
               DsiplacementBCSettings object.

           Raises:
               None.
        """

        if not fix_x and not fix_y and not fix_z and \
           not prevent_x_rot and not prevent_y_rot and \
           not prevent_z_rot:
            raise RuntimeError("The boundary condition settings don't restrict anything!")

        self.fix_x = fix_x
        self.fix_y = fix_y
        self.fix_z = fix_z
        self.prevent_x_rot = prevent_x_rot
        self.prevent_y_rot = prevent_y_rot
        self.prevent_z_rot = prevent_z_rot



class SurfaceBC():
        
    def __init__(self, ngon: geom.NGon3D, bc_settings: DisplacementBCSettings) -> None:
        """Creates a boundary condition on some surface.
           
           Args:
               ngon:        The surface on which the boundary condition applies.
               bc_settings: The restrictions for the surface.

           Returns:
               None.

           Raises:
               None.
        """
    
        self.region = ngon
        self.bc_settings = bc_settings



class VertexBC():

    def __init__(self, vertices: list[geom.Point3D], bc_settings: DisplacementBCSettings 
                ) -> None:
        """Creates a boundary condition on some set of points.

           Args:
               vertices:    The points on which to create the boundary conditions.
               bc_settings: The restrictions on the points.

           Returns:
               None.

           Raises:
               None.
        """

        self.region = vertices
        self.bc_settings = bc_settings


BC = TypeVar("BC", SurfaceBC, VertexBC)


def apply_BCs(BCs: list[BC], step_name: str, part_instance: Any, model_name: str, 
              mdb: Any) -> None:
    """Applies boundary conditions to an assembly.

       Note that boundary conditions can only be applied to an Abaqus Assembly, they
           cannot be applied to an Abaqus Part.

       In order to apply boundary conditions to regions which don't already exist, 
           it's necessary to partition pre-existing features. In general, an Abaqus 
           Part object or an Abaqus PartInstance object can be partitioned. This 
           function creates new regions on parts, when necessary, by partitioning 
           this single Abaqus PartInstance object.
       Since a PartInstance object must be passed, it means that the Assembly is 
           not empty. In my experience, if you try to create a boundary condition 
           on a region of an Abaqus Part object and the Assembly is empty, Abaqus 
           segfaults.

       Args:
           BCs:           The boundary conditions to apply.
           step_name:     The first step at which the boundary conditions should be
                              active. By default, the boundary conditions propagate
                              to all future steps.
           part_instance: Abaqus PartInstance object. The instance on which to
                              apply the boundary conditions.
           model_name:    The name of the model.
           mdb:           Abaqus MDB object.

       Returns:
           None.

       Raises:
           None.
    """
    for bc in BCs:
        bc_cnt = shim.get_bc_cnt(model_name, mdb)
        bc_name = shim.STANDARD_BC_PREFIX + str(bc_cnt)

        if isinstance(bc, SurfaceBC):
            assembly = mdb.models[model_name].rootAssembly
            face = shim.partition_face(bc.region, instance=part_instance, assembly=assembly) 
            region = shim.build_region_with_face(face, part_instance)
        elif isinstance(bc, VertexBC):
            assert False, "Not yet supported."
        else:
            assert False, "Not yet supported."

        if isinstance(bc.bc_settings, DisplacementBCSettings):
            shim.create_disp_rot_bc(bc_name, step_name, region, bc.bc_settings, model_name, mdb)
        else:
            assert False, "Not yet supported."




