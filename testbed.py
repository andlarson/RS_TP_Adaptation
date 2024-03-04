"""
This file acts as a testbed for individual components of the library.
"""

import numpy as np

from abaqus import *
from abaqusConstants import *
import odbAccess

from src.util.debug import *


if __name__ == "__main__":
   
    odb = odbAccess.openOdb(path="/home/andlars/Desktop/RS_TP_Adaptation/experiments/experiments/full_flow/stress_estimation/first_plan/first_plan.odb", readOnly=True)

    step_name = odb.steps.keys()[0]

    # The step potentially contains many frames.
    # Since we care about the displacement of the points across the whole step,
    #     it's necessary to find the first and last frames in the step.
    step = odb.steps[step_name]
    smallest_increment_number = step.frames[0].incrementNumber
    largest_increment_number = step.frames[0].incrementNumber
    for frame in step.frames: 
        if frame.incrementNumber < smallest_increment_number:
            smallest_increment_number = frame.incrementNumber
        elif frame.incrementNumber > largest_increment_number:
            largest_increment_number = frame.incrementNumber
    first_frame = step.frames[smallest_increment_number]
    last_frame = step.frames[largest_increment_number]
    
    last_frame_displacements = last_frame.fieldOutputs["U"]
    node_cnt = len(last_frame_displacements.values)
    
    # DEBUG
    dp("The total number of nodes is " + str(node_cnt))

    nodes = last_frame_displacements.values[0].instance.nodes

    for idx in range(node_cnt):
        coords = nodes[idx].coordinates
        if coords[0] < 0 or coords[0] > 40 or coords[1] < 0 or coords[1] > 10 or coords[2] < 0 or coords[2] > 400:
            dp("Anomolous point detected!")
    
