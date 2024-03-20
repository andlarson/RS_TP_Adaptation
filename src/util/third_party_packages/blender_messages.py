"""
This file defines the interface between the parent process running Abaqus' 
    Python interpreter and a child process running the code which interacts
    with Blender's Python API via a default Python interpreter.
"""

import enum

from typing import *


STANDARD_RESPONSE_MESSAGE = "DONE"


# *****************************************************************************
#                                 Message Types 
# *****************************************************************************

class Messages(enum.Enum):
    """The types of messages which can be sent from parent to child and from
           child to parent.
       Each message type requires 4 functions: One to pack data to be sent to
           the child process, one to pack data to be sent to the parent process, 
           one to unpack data receieved in the parent process, one to unpack
           data received in the child process.
    """
    BLENDER_REMESH = "BLENDER_REMESH"
    BLENDER_VOLUME_SYMMETRIC_DIFFERENCE = "BLENDER_VOLUME_SYMMETRIC_DIFFERENCE"
    BLENDER_REMOVE_MATERIAL = "BLENDER_REMOVE_MATERIAL"



# *****************************************************************************
#            Functions For Message Passing, Grouped By Message Type
# *****************************************************************************


# ****************************** REMESH ***************************************

def parent_to_child_remesh_data(path_to_stl: str, path_to_save_loc: str) -> str:
    """Constructs message data for remeshing.

       Args:
           path_to_stl:      Absolute path to .stl file to remesh.
           path_to_save_loc: Absolute path to desired save location of .stl file
                                 which contains the remeshed gometry.
       
       Returns:
           Formatted message data.
    
       Raises:
           None.
    """

    return path_to_stl + "###" + path_to_save_loc



def unpack_remesh_data_at_child(message_data: str) -> list[str]:
    """Unpacks data sent by parent requesting remeshing.

       Args:
           message_data: Data received by the child. Assumed that this data
                             is unmodified (i.e. it is an unchanged copy of
                             what was received by the child).
       
       Returns:
           Absolute path to .stl file to remesh and absolute path to location
               where the remeshed .stl should be saved.
    
       Raises:
           None.
    """
    
    datum = message_data.rsplit("###")
    assert(len(datum) == 2)
    return datum 



def child_to_parent_remesh_data() -> str:
    """Constructs data for response to remeshing request.

       Args:
           None.
       
       Returns:
           Formatted message data.
    
       Raises:
           None.
    """

    return STANDARD_RESPONSE_MESSAGE 



def unpack_remesh_data_at_parent(message_data: str):
    """Unpacks data sent by child responding to remesh request. 

       Args:
           message_data: Data received by the parent. Assumed that this data
                             is unmodified.
       
       Returns:
           None.
    
       Raises:
           None.
    """
    
    assert message_data == STANDARD_RESPONSE_MESSAGE



# ************************ SYMMETRIC DIFFERENCE ********************************

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

    return path_to_stl1 + "###" + path_to_stl2



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
    
    datum = message_data.rsplit("###")
    assert(len(datum) == 3)
    return datum 



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



# ************************** MATERIAL REMOVAL *********************************

def parent_to_child_remove_material_data(path_to_stl1: str, path_to_stl2: str, 
                                         path_to_resulting_stl: str) -> str:
    """Constructs message data for computing the symmetric difference between
           two meshes.

       Args:
           path_to_stl1:          Absolute path to .stl file. The geometry from 
                                      which to remove material.
           path_to_stl2:          Absolute path to .stl file. The geometry which to
                                      remove.
           path_to_resulting_stl: Absolute path to desired location of resulting
                                      .stl file. This is the location that the
                                      child process will save the resulting .stl
                                      file to.
       
       Returns:
           Formatted message data. 
    
       Raises:
           None.
    """

    return path_to_stl1 + "###" + path_to_stl2 + "###" + path_to_resulting_stl



def unpack_material_removal_data_at_child(message_data: str) -> list[str]:
    """Unpacks data sent by parent requesting a material removal operation.

       Args:
           message_data: Data received by the child. Assumed that this data
                             is unmodified (i.e. it is an unchanged copy of
                             what was received by the child).
       
       Returns:
           Absolute paths to three .stl files. The first describes the base
               geoemtry. The second describes the geometry to remove from the
               base geoemtry. The third is the desired save location of the
               resulting geometry.
    
       Raises:
           None.
    """
    
    datum = message_data.rsplit("###")
    assert(len(datum) == 3)
    return datum 



def child_to_parent_material_removal_data() -> str:
    """Constructs data for response to material removal request.

       Args:
           None.
       
       Returns:
           Formatted message data.
    
       Raises:
           None.
    """

    return STANDARD_RESPONSE_MESSAGE



def unpack_material_removal_data_at_parent(message_data: str):
    """Unpacks data sent by child responding to material removal request. 

       Args:
           message_data: Data received by the parent. Assumed that this data
                             is unmodified.
       
       Returns:
           None.
    
       Raises:
           None.
    """
    
    assert message_data == STANDARD_RESPONSE_MESSAGE



