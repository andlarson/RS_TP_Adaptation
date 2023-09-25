"""
This file contains code related to tool passes.
"""

import copy

import util.geom as geom
from util.debug import *



class ToolPass:
    
    def __init__(self, tool_pass_geom):
    # type: (geom.SpecRightRectPrism) -> None

        self.geom = tool_pass_geom
    


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