"""
This file defines the interface between the parent process running Abaqus' 
    Python interpreter and a child process running the code which interacts
    with Blender's Python API via a default Python interpreter.
"""


# The messages include: 
#     REMESH 
#     VOLUME_SYMMETRIC_DIFFERENCE 
#
# Each message is accompanied by a single line of data. The data which
#     accompanies each message is:
#     REMESH -> Absolute path to .stl file. The .stl file contains some
#                   complex mesh.
#     VOLUME_SYMMETRIC_DIFFERENCE -> Two absolute paths to .stl files. The
#                                        paths are separated by a single
#                                        space.
messages_parent_process_to_child_process: list[str] = ["BLENDER_REMESH", "BLENDER_VOLUME_SYMMETRIC_DIFFERENCE"]

# The messages include: 
#     REMESH
#     VOLUME_SYMMETRIC_DIFFERENCE 
#
# Each message is accompanied by a single line of data. The data which
#     accompanies each message is:
#     REMESH -> Absolute path to .stl file. The .stl file contains a
#                   simplified mesh.
#     VOLUME_SYMMETRIC_DIFFERENCE -> A float representing the volume. 
messages_child_process_to_parent_process: list[str] = ["BLENDER_REMESH", "BLENDER_VOLUME_SYMMETRIC_DIFFERENCE"]





