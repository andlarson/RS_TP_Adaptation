"""
It is expected that the code contained in this file runs on a standard Python
    interpreter. This makes it feasible for this code to call code in third 
    party packages. Calling code in third party packages is nearly impossible
    from code which is running on the custom Python interpreter included with
    Abaqus.

When the code in this file runs, it is running in the context of a child process
    spawned by a parent process (the parent process is code running on Abaqus'
    custom Python interpeter). The communication pattern between the child 
    process and the parent process is defined and discussed in 
    calling_third_parties.py.

At a high level, this process simply polls stdin, waiting for messages to be
    sent to it. When it receives a message of a known type, it does the work
    associated with that message and sends the result of the work back
    to the parent process by writing to stdout. It then continues to poll, 
    waiting for more messages to arrive. It is expected that this process 
    should never exit on its own -- the spawning parent process will kill 
    it when it no longer has any more work for it to do.
"""

import sys
import time
import os
from typing import *

# Resolving imports.
sys.path.append("/home/andlars/Desktop/RS_TP_Adaptation/")

import src.util.third_party_packages.blender_messages as blender_messages
import src.util.third_party_packages.child_process.blender as blender



def process_blender_message(message: tuple[str, str]) -> tuple[str, str]:
    """Process a message intended to induce some work in Blender.
                   
       Args:
           message: Message which triggers some work in Blender. Assumed to
                        have format (message_type, message_data).
    
       Returns:
           Result of the work in Blender. Always takes the form (message type,
               message data).
    
       Raises:
           None.
    """
    
    message_type = message[0]
    message_data = message[1]
    
    if message_type == blender_messages.Messages.BLENDER_REMESH.value:
        # Unpack the received data.
        stl_to_remesh, save_location = blender_messages.unpack_remesh_data_at_child(message_data)

        blender.remesh(stl_to_remesh, save_location) 

        # Pack the response data. 
        packed_data = blender_messages.child_to_parent_remesh_data()

        return (blender_messages.Messages.BLENDER_REMESH.value, packed_data) 
    elif message_type == blender_messages.Messages.BLENDER_VOLUME_SYMMETRIC_DIFFERENCE.value:
        # Unpack the received data.
        geom1, geom2 = blender_messages.unpack_symm_diff_data_at_child(message_data)

        volume = blender.compute_volume_symmetric_difference(geom1, geom2) 

        # Pack the response data.
        packed_data = blender_messages.child_to_parent_symm_diff_data(volume)

        return (blender_messages.Messages.BLENDER_VOLUME_SYMMETRIC_DIFFERENCE.value, packed_data)
    elif message_type == blender_messages.Messages.BLENDER_REMOVE_MATERIAL.value:
        # Unpack the received data.
        base_geom, to_remove, save_loc = blender_messages.unpack_material_removal_data_at_child(message_data)

        blender.remove_material(base_geom, to_remove, save_loc)

        # Pack the response data.
        packed_data = blender_messages.child_to_parent_material_removal_data()

        return (blender_messages.Messages.BLENDER_REMOVE_MATERIAL.value, packed_data)
    else:
        raise AssertionError("Unknown message type.")



def process_message(message: tuple[str, str]) -> tuple[str, str]:
    """Hands off message processing to the appropriate function based on
           the type of message."""
    
    message_type = message[0]
   
    if message_type in [member.value for member in blender_messages.Messages]:
        return process_blender_message(message)
    else:
        raise RuntimeError("Unknown message type.")



if __name__ == "__main__":
    
    while True:
        # Assumes an ordering. 
        message_type = sys.__stdin__.readline().removesuffix("\n")
        message_data = sys.__stdin__.readline().removesuffix("\n")

        ret = process_message((message_type, message_data))
        sys.stdout.writelines([ret[0] + "\n", ret[1] + "\n"])
        sys.stdout.flush()

