"""
This file defines the interface between the parent process running Abaqus' 
    Python interpreter and a child process running the code which interacts
    with Blender's Python API via a default Python interpreter.
"""

import enum

from typing import *


# *****************************************************************************
#                                 Parent Side 
# *****************************************************************************

class MessageParentToChild(enum.Enum):
    BLENDER_REMESH = "BLENDER_REMESH"
    BLENDER_VOLUME_SYMMETRIC_DIFFERENCE = "BLENDER_VOLUME_SYMMETRIC_DIFFERENCE"


# The following functions standardize the construction of messsage data for
#     parent-to-child messaging.
# There is a single function for each possible message.

def parent_to_child_remesh_data(path_to_stl: str) -> str:
    """Constructs message data for remeshing.

       Args:
           path_to_stl: Absolute path to .stl file to remesh.
       
       Returns:
           Formatted message data.
    
       Raises:
           None.
    """

    return path_to_stl 



def parent_to_child_symm_diff_data(path_to_stl1: str, path_to_stl2: str) -> str:
    """Constructs message data for computing the symmetric difference between
           two meshes.

       Args:
           path_to_stl1: Absolute path to .stl file.
           path_to_stl2: Absolute path to .stl file.
       
       Returns:
           Formatted message data. 
    
       Raises:
           None.
    """

    return path_to_stl1 + " " + path_to_stl2



# The following functions standardize the unpacking of message data which arrives
#     at the child.

def unpack_remesh_data_at_child(message_data: str) -> str:
    """Unpacks data sent by parent requesting remeshing.

       Args:
           message_data: Data received by the child. Assumed that this data
                             is unmodified (i.e. it is an unchanged copy of
                             what was received by the child).
       
       Returns:
           Absolute path to .stl file to remesh.
    
       Raises:
           None.
    """

    return message_data



def unpack_symm_diff_data_at_child(message_data: str) -> list[str]:
    """Unpacks data sent by parent requesting a symmetric difference computation.

       Args:
           message_data: Data received by the child. Assumed that this data
                             is unmodified (i.e. it is an unchanged copy of
                             what was received by the child).
       
       Returns:
           Absolute paths to two .stl files for which the symmetric difference
               operation should be done. 
    
       Raises:
           None.
    """
    
    stl_paths = message_data.rsplit(" ")
    assert(len(stl_paths) == 2)
    return stl_paths



# *****************************************************************************
#                                   Child Side 
# *****************************************************************************

class MessageChildToParent(enum.Enum):
    BLENDER_REMESH = "BLENDER_REMESH"
    BLENDER_VOLUME_SYMMETRIC_DIFFERENCE = "BLENDER_VOLUME_SYMMETRIC_DIFFERENCE"



# The following functions standardize the construction of messsage data for
#     child-to-parent messaging.
# There is a single function for each possible message.

def child_to_parent_remesh_data(path_to_stl: str) -> str:
    """Constructs data for response to remshing request.

       Args:
           path_to_stl: Absolute path to .stl file containing the results of
                            remeshing.
       
       Returns:
           Formatted message data.
    
       Raises:
           None.
    """

    return path_to_stl



def child_to_parent_symm_diff_data(volume: float) -> str:
    """Constructs data for response to symmetric difference request.

       Args:
           volume: Volume of symmetric difference.
       
       Returns:
           Formatted message data.
    
       Raises:
           None.
    """
    
    # Beware, this cuts off some digits.
    return str(volume)



# The following functions standardize the unpacking of message data which arrives
#     at the parent.

def unpack_remesh_data_at_parent(message_data: str) -> str:
    """Unpacks data sent by child responding to remesh request. 

       Args:
           message_data: Data received by the parent. Assumed that this data
                             is unmodified.
       
       Returns:
           Absolute path to .stl file containing remeshed geometry.
    
       Raises:
           None.
    """

    return message_data



def unpack_symm_diff_data_at_parent(message_data: str) -> float:
    """Unpacks data sent by child responding to symmetric difference calculation
           request.

       Args:
           message_data: Data received by the parent. Assumed that this data
                             is unmodified.
       
       Returns:
           Volume of symmetric difference.
    
       Raises:
           None.
    """
    
    return float(message_data)
