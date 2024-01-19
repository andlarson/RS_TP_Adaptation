"""
This file contains code related to tool passes.
"""

import copy

import numpy as np

import src.util.geom as geom

from src.util.debug import *



class ToolPass:

    # A cylindrical tool pass which is not self-intersecting! No checks are done
    #    for self-intersection upon creation of the this object, but when the
    #    underlying object is created in Abaqus, an unrecoverable failure will
    #    occur. 
    #
    # Notes:
    #    See notes/toolpath/toolpath_orientation_1.jpg and notes/toolpath/toolpath_
    #       orientation_2.jpg to understand how the toolpath is oriented with
    #       respect to the path.
    #    
    # Arguments:
    #    path   - PlanarCubicC2Spline3D object.
    #    radius - Float.
    #             Radius of the cylindrical tool.
    #    length - Float.
    #             Length of the cylindrical tool.
    #
    # Return:
    #    ToolPass object. 
    def __init__(self, path, radius, length):
    # type: (geom.PlanarCubicC2Spline3D, float, float) -> None

        self.path = path
        self.radius = radius
        self.length = length



class ToolPassPlan:
    
    def __init__(self, tool_passes):
    # type: (list[ToolPass]) -> None

        self.plan = tool_passes
        self.passes_done = []
        self.passes_todo = copy.deepcopy(tool_passes)


    def add(self, tool_pass):
    # type: (ToolPass) -> None

        self.passes_todo.append(tool_pass)

    
    def pop(self):
    # type: (None) -> ToolPass 

        if len(self.passes_todo) == 0:
            return None
        else:
            tool_pass = self.passes_todo.pop(0)
            self.passes_done.append(tool_pass)

            return tool_pass


    def done_so_far(self):
    # type: (None) -> int 

        return len(self.passes_done)



# Check if two tool passes have exact same representation. 
#
# Notes:
#    None.
#
# Arguments:
#    tp1 - ToolPass.
#    tp2 - ToolPass.
#
# Returns:
#    Boolean. 
def compare_tool_passes(tp1, tp2):

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
        raise RuntimeError("Not yet supported.")



# Checks if two tool pass plans exactly match in order and representation.
#
# Notes:
#    This function does not care if the two tool pass plans have been simulated,
#       partially simulated, or not simulated at all. It just checks if the tool
#       passes which underlie the two plans exactly match in order and 
#       representation. 
#
# Arguments:
#    tpp_1 - ToolPassPlan.
#    tpp_2 - ToolPassPlan.
#
# Returns:
#    Boolean.
def compare_tool_pass_plans(tpp_1, tpp_2):
# type: (ToolPassPlan, ToolPassPlan) -> bool

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



# Generate bounding box geometry for a tool pass.
# 
# Notes:
#    Computations assume an orientation between the tool and the tool pass path.
# 
# Arguments:
#    x_excess   - Float.
#                 Amount of excess in the x direction that the user wants for the
#                    bounding boxes. Expressed in units of the global coordinate
#                    system.
#    y_excess   - Float.
#                 Amount of excess in the y direction that the user wants for the
#                    bounding boxes. Expressed in units of the global coordinate
#                    system.
#    z_excess   - Float.
#                 Amount of excess in the z direction that the user wants for the
#                    bounding boxes. Expressed in units of the global coordinate
#                    system.
#    tool_pass - ToolPass object.
#
# Returns:
#    SpecRightRectPrism object.
def create_tool_pass_bounding_box(x_excess, y_excess, z_excess, tool_pass):
# type: (float, float, float, ToolPass) -> geom.SpecRightRectPrism

    # For this type of path, all the bounding boxes have faces defined by the same
    #    y values. 
    assert(isinstance(tool_pass.path, geom.PlanarCubicC2Spline3D))
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
      
