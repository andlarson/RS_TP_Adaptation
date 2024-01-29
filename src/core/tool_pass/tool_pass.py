"""
This file contains code related to tool passes.
"""

import copy

import numpy as np

import src.util.geom as geom

from src.util.debug import *



class ToolPass:

    def __init__(self, path: geom.PlanarCubicC2Spline3D, radius: float, length: float) -> None:
        """Creates a cylindrical tool pass which should not be self-intersecting. 
   
           No checks are done for self-intersection upon creation of the this object, 
               but when the underlying object is created in Abaqus, an unrecoverable 
               failure will occur. 
           
           Args:
              path:   Path that the tool movies in space. See notes/toolpath/toolpath_orientation_1.jpg 
                          and notes/toolpath/toolpath_orientation_2.jpg to understand 
                          how the toolpath is oriented with respect to the path.
              radius: Radius of the cylindrical tool.
              length: Length of the cylindrical tool.
           
           Return:
               None.

           Raises:
               None.
        """

        self.path = path
        self.radius = radius
        self.length = length



class ToolPassPlan:
    
    def __init__(self, tool_passes: list[ToolPass]) -> None:
        """Creates a tool pass plan, which is just a sequence of tool passes."""

        self.plan = tool_passes
        self.passes_done = []
        self.passes_todo = copy.deepcopy(tool_passes)


    def add(self, tool_pass: ToolPass) -> None:
        """Adds a tool pass to a tool pass plan."""

        self.passes_todo.append(tool_pass)

    
    def pop(self) -> ToolPass:
        """Returns the first tool pass in the tool pass plan, if there are any
               tool passes in the plan."""

        if len(self.passes_todo) == 0:
            return None

        tool_pass = self.passes_todo.pop(0)
        self.passes_done.append(tool_pass)

        return tool_pass


    def done_so_far(self) -> int:
        """Returns the number of tool passes in this tool pass plan which were
               done thus far."""

        return len(self.passes_done)



def compare_tool_passes(tp1: ToolPass, tp2: ToolPass) -> bool:
    """Checks if two tool passes are exactly the same. 
        
       This function looks for exact matches and, when doing so, it compares
           floating point numbers. When comparing floating point numbers,
           it does not allow any leeway - two floats are equal iff they have
           the exact same representation.
        
       Arguments:
           tp1: A tool pass. 
           tp2: Another tool pass. 
       
       Returns:
           Boolean. 
       
       Raises:
           None.
    """

    if type(tp1.geom) != type(tp2.geom):
        return False

    # Compare types on a case-by-case basis.
    if isinstance(tp1.geom, geom.SpecRightRectPrism):

        if len(tp1.geom.vertices) == len(tp2.geom.vertices):

            for v1 in tp1.geom.vertices:
                if v1.components() not in [v2.components() for v2 in tp2.geom.vertices]:
                    return False
            
            for v2 in tp2.geom.vertices:
                if v2.components() not in [v1.components() for v1 in tp1.geom.vertices]:
                    return False

        else:
            return False

        return True
    
    else:
        raise AssertionError("Not yet supported.")



def compare_tool_pass_plans(tpp_1: ToolPassPlan, tpp_2: ToolPassPlan) -> bool:
    """Checks if two tool pass plans exactly match in order and representation.
    
       This function does not care if the two tool pass plans have been simulated,
           partially simulated, or not simulated at all. It just checks if the tool
           passes which underlie the two plans exactly match in order and 
           representation. 
       
       Args:
           tpp_1: A sequence of tool passes. 
           tpp_2: Another sequence of tool passes. 
       
       Returns:
           Boolean.
       
       Raises:
           None.
    """

    if len(tpp_1.plan) == len(tpp_2.plan):
        for tp1 in tpp_1.plan:

            in_plan_2 = False            
            for tp2 in tpp_2.plan:
                if compare_tool_passes(tp1, tp2):
                    in_plan_2 = True
                    break

            if not in_plan_2:
                return False

        for tp2 in tpp_2.plan:

            in_plan_1 = False            
            for tp1 in tpp_1.plan:
                if compare_tool_passes(tp1, tp2):
                    in_plan_1 = True
                    break

            if not in_plan_1:
                return False

    else:
        return False

    return True



def create_tool_pass_bounding_box(x_excess: float, y_excess: float, z_excess: float, 
                                  tool_pass: ToolPass) -> geom.SpecRightRectPrism:
    """DEPRECATED. No longer using bounding boxes to do local mesh refinement.
       
       Generates bounding box geometry for a tool pass. Only works for tool passes
           defined by a path which lives on a plane in space.
       
       Implementation Details:
       Assumes an orientation between the tool and the tool pass path.
       
       Args:
           x_excess:  Amount of excess in the x direction that the user wants for 
                          the bounding boxes. Expressed in units of the global 
                          coordinate system.
           y_excess:  Amount of excess in the y direction that the user wants for 
                          the bounding boxes. Expressed in units of the global 
                          coordinate system.
           z_excess:  Amount of excess in the z direction that the user wants for 
                          the bounding boxes. Expressed in units of the global 
                          coordinate system.
           tool_pass: The tool pass around which to create the bounding box.
       
       Returns:
           Bounding box geometry.
       
       Raises:
           None.
    """

    if not isinstance(tool_pass.path, geom.PlanarCubicC2Spline3D):
        raise RuntimeError("Tool pass path not defined on a plane.")

    # For this type of path, all the bounding boxes have faces defined by the same
    #    y values. 
    max_y = tool_pass.path.y + tool_pass.length + y_excess
    min_y = tool_pass.path.y - y_excess 

    max_x, min_x, _, _, max_z, min_z = geom.find_extrema(tool_pass.path.v_list)
    max_x = max_x + x_excess
    min_x = min_x - x_excess
    max_z = max_z + z_excess
    min_z = min_z - z_excess

    p1 = geom.Point3D(np.array((max_x, max_y, max_z))) 
    p2 = geom.Point3D(np.array((max_x, max_y, min_z)))
    p3 = geom.Point3D(np.array((max_x, min_y, max_z)))
    p4 = geom.Point3D(np.array((max_x, min_y, min_z)))
    p5 = geom.Point3D(np.array((min_x, min_y, min_z)))
    p6 = geom.Point3D(np.array((min_x, max_y, min_z)))
    p7 = geom.Point3D(np.array((min_x, max_y, max_z)))
    p8 = geom.Point3D(np.array((min_x, min_y, max_z)))
    bounding_box = geom.SpecRightRectPrism(p1, p2, p3, p4, p5, p6, p7, p8)

    return bounding_box
      
