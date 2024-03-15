"""
This file contains functionality which makes it easier to call and use code written
    by third parties.
The Abaqus Python API can only be accessed through a custom Python interpreter
    included with the installation of Abaqus. This custom interpreter is based
    on the Python 3.10.5 interpreter. Along with the interpreter, the install 
    includes commonly used packages like numpy, scipy, and matplotlib. The install
    does not include pip, so installing other packages is non-trivial. 
There are two approachs others have taken:
    1) Install a package via pip in the context of a separate, default Python
           interpreter. Then copy the resulting files from the site-packages/
           directory into the site-packages/ directory of the custom Python
           interpreter included with Abaqus, including all dependencies. This
           approach seems difficult to maintain and hacky.
    2) In code running on the custom Python interpreter, spawn a child process
           which runs code on a separate, default Python interpreter. The separate,
           default Python interpreter can have packages installed with it via
           its pip functionality. Establish an interprocess communication mechanism
           so that the child process can communicate the results with the
           spawning parent process.
This file includes functionality which enables approach (2).

Relevant documentation includes:
    https://r1132100503382-eu1-3dswym.3dexperience.3ds.com/community/swym:prd:R1132100503382:community:39?content=swym:prd:R1132100503382:qnaquestion:aRq4LlR1SLao3xrmWH_uDA
    Dassault Systems Knowledge Base articles QA00000051470 and QA00000008399.
"""

import subprocess
import copy
import os
import time
from typing import *


class UseThirdPartyPackage:
    
    def __init__(self, path_to_python: str, script: str) -> None:
        """Sets things up to get access to third-party packages. 

           The child-parent interaction must be fully synchronous.
           It is expected that the child process reads data from stdin, processes
               the data from stdin according to its type, passes the result back 
               to the parent, and continues this loop until the parent process 
               kills it. This is a message passing scheme.
           Since raw data is exchanged between the parent and child, it's always
               necessary to specify the data type in any exchange. A general way
               to solve this is to specify an initial communication interaction
               pattern wherein the types are shared between parent and child. The
               quicker and easier way to solve this is restrict the message passing
               pattern to a single simple, hard-coded scheme. This class takes
               the latter approach. 
           To make things simple, this class assumes that every parent-child 
               interaction takes place over the course of two lines. The first 
               line should describe the type of the message and the second 
               line should contain the message content. Since the message content
               is restricted to a single line, a file should be used to exchange
               lots of data. This class makes no assumptions about the types of
               the messages, but it does assume this two-line message passing
               scheme.
           
           TODO: Really, should be using Python's asyncio functionality so that
               the receive_data() function doesn't need a sleep loop.
                       
           Args:
               path_to_python: Absolute path to desired Python interpreter.
               script:         Absolute path to Python .py file to execute with 
                                   the interpreter.
        
           Returns:
               None.
        
           Raises:
               None.
        """

        self.path_to_python: str = path_to_python
        self.script: str = script
        
        # The child process should not use the default search path when looking
        #     for module files.
        self.env: Mapping[str, str] = copy.deepcopy(os.environ)
        del self.env['PYTHONPATH']

        self.sp: subprocess.Popen | None = None
       

    def start_child(self) -> None:
        """Starts the child process.
                       
           Args:
               None.
        
           Returns:
               None.
        
           Raises:
               None.
        """
        
        sp = subprocess.Popen(args=[self.path_to_python, self.script], 
                              env=self.env, stdin=subprocess.PIPE,
                              stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        self.sp = sp


    def kill_child(self) -> None:
        """Kills the child process if it has not already terminated.
                       
           Args:
               None.
        
           Returns:
               None.
        
           Raises:
               None.
        """
        
        assert self.sp is not None

        self.sp.poll()

        if self.sp.returncode is None:
            self.sp.kill() 


    def send_data(self, message_type: str, message_data: str) -> None:
        """Passes message to the child process in standard format. 

           Assumes that the child process has not terminated.

           Args:
               message_type: The type of the message to pass to the child process.
                                 A single line. 
               message_data: The data to pass to the child process. A single line.
        
           Returns:
               None.
        
           Raises:
               None.
        """
        
        assert self.sp is not None
        self.sp.poll()
        assert self.sp.returncode is None 

        assert len(message_type.splitlines()) == 1
        assert len(message_data.splitlines()) == 1
        
        self.sp.stdin.write(message_type)
        self.sp.stdin.write(message_data)


    def receive_data(self, timeout=60) -> tuple[str, str]:
        """Periodically sleeps until receiving a message from the child process. 
           
           Should only be called when the caller knows that a message will be
               sent by the child process. 
                       
           Args:
               timeout: The timeout period in seconds. 
        
           Returns:
               A tuple containing the message type and message data.
        
           Raises:
               AssertionError: Thrown when the timeout is exceeded.
        """
        
        start_time = time.time()

        message_type = self.sp.stdout.readline()
        message_data = self.sp.stdout.readline()
        
        while message_type == "" or message_data == "":
            cur_time = time.time()
            if cur_time - start_time > timeout:
                raise AssertionError("Timeout exceeded!")
            time.sleep(3)            

            message_type = self.sp.stdout.readline()
            message_data = self.sp.stdout.readline()

        return (message_type, message_data)
