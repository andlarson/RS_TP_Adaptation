"""
This file contains functionality related to the estimation of residual stress
    fields.
"""

from typing import Any

import numpy as np

import src.util.geom as geom



class StressTensor:

    def __init__(self, sigma_xx: float, sigma_yy: float, sigma_zz: float,
                 sigma_xy: float, sigma_xz: float, sigma_yz: float) -> None:
        """Creates a stress tensor. Represents the state of stress at a single
               point in space.
            
           Args:
               sigma_xx: Component of stress acting on the x plane (i.e. the
                             plane normal to the x axis) in the +x dirrection.
               sigma_yy: Component of stress acting on the y plane (i.e. the
                             plane normal to the y axis) in the +y direction.
               sigma_zz: Component of stress acting on the z plane (i.e. the
                             plane normal to the z axis) in the +z direction.
               sigma_xy: Component of stress acting on the x plane (i.e. the
                             plane normal to the x axis) in the +y direction.
               sigma_xz: Component of stress acting on the x plane (i.e. the
                             plane normal to the x axis) in the +z direction.
               sigma_yz: Component of stress acting on the y plane (i.e. the
                             plane normal to the y axis) in the +z direction.
               
               The units of all arguments are Force/Area.
        
           Returns:
               None.
        
           Raises:
               None.
        """

        self.rep = np.array([[sigma_xx, sigma_xy, sigma_xz], [sigma_xy, sigma_yy, sigma_yz], [sigma_xz, sigma_yz, sigma_zz]])
    


class ConstantResidualStressField:

    def __init__(self, region: geom.PlanarCubicC2Spline3D, tensor: StressTensor) -> None:
        """Creates a residual stress field in a region which is constant. A
               constant residual stress field is defined by the same stress
               tensor at every point.
            
           Args:
               region: The region that the residual stress field exists.
               tensor: The stress at every point in the region.
        
           Returns:
               None.
        
           Raises:
               None.
        """

        self.region = region
        self.tensor = tensor



