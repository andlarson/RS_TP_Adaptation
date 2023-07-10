import core.machining.machining as mach
import util.geom as geom
import core.part.part as part
import core.tool_pass.tool_pass as tp

import sys

# DEBUG
from util.debug import *


if __name__ == "__main__":

    # Need a part and a tool pass plan.    
   
    # The initial geometry with a stress profile defined in Abaqus.
    path_to_cae = "/home/andlars/Downloads/script_testing/test_initial_geometry.cae"
    abaqus_part = part.AbaqusDefinedPart("an_example_part", path_to_cae)
    
    dp("Initial geometry loaded.")

    # Specifying the tool passes. 

    # Tool Pass #1
    v1 = geom.Point3D(40, 9, 215)
    v2 = geom.Point3D(40, 20, 215)
    v3 = geom.Point3D(40, 9, 185)
    v4 = geom.Point3D(40, 20, 185)
    v5 = geom.Point3D(0, 9, 215)
    v6 = geom.Point3D(0, 20, 215)
    v7 = geom.Point3D(0, 9, 185)
    v8 = geom.Point3D(0, 20, 185)

    tp1_shape = geom.SpecRightRectPrism(v1, v2, v3, v4, v5, v6, v7, v8)

    tp1 = tp.ToolPass(tp1_shape)

    tool_passes = [tp1]

    tool_pass_plan = tp.ToolPassPlan(tool_passes)

    dp("Tool pass has been specified.\n") 

    # Building the machining object.
    machining_process = mach.MachiningProcess(None, abaqus_part, tool_pass_plan)

    # And simulate the next tool pass.
    machining_process.sim_next_tool_pass("/home/andlars/Downloads/script_testing/test_post_tool_pass.cae")





