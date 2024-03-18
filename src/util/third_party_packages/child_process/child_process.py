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
from typing import *

import src.util.third_party_packages.blender_messages as blender_messages
import src.util.third_party_packages.child_process.blender as blender



def process_blender_message(message: tuple[str, str]) -> tuple[str, str]:
    """Process a message intended to induce some work in Blender.
                   
       Args:
           message: Message which triggers some work in Blender. Assumed to
                        have format (message_type, message_data).
    
       Returns:
           Result of the work in Blender. The interpretation of the result
               depends on the message which was received.
    
       Raises:
           None.
    """
    
    message_type = message[0]
    message_data = message[1]
    
    if message_type == blender_messages.messages_parent_process_to_child_process[0]:
        response_data = blender.remesh(message_data) 
        return (blender_messages.messages_child_process_to_parent_process[0], response_data) 
    elif message_type == blender_messages.messages_parent_process_to_child_process[1]:
        response_data = blender.compute_volume_symmetric_difference(message_data) 
        return (blender_messages.messages_child_process_to_parent_process[1], response_data) 
    else:
        raise AssertionError("Unknown message type.")



def process_message(message: tuple[str, str]) -> tuple[str, str]:
    """Hands off message processing to the appropriate function based on
           the type of message."""
    
    message_type = message[0]

    # Blender.
    if message_type in blender_messages.messages_parent_process_to_child_process:
        return process_blender_message(message)
    else:
        raise RuntimeError("Unknown message type.")



if __name__ == "__main__":

    while True:
    
        # Check if a (message type, message) pair is available.
        # Danger here! If a write from the parent process is not atomic, then
        #     it's possible that only a partial message would be read here.
        lines = sys.stdin.readlines(2)
    
        # If no message is available, sleep for a while.
        if len(lines) != 2:
            time.sleep(3)  
        else:
            ret = process_message(tuple(lines))
            sys.stdout.writelines(ret)
