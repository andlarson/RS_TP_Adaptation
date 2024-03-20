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

from src.util.debug import *


# TODO: Yuck! Hard-coded!
PATH_TO_STANDARD_INTERPRETER = "/home/andlars/Desktop/executables_etc/bin/python3"
PATH_TO_STANDARD_PARENT_SCRIPT = "/home/andlars/Desktop/RS_TP_Adaptation/src/util/third_party_packages/parent_process/calling_third_parties.py"



class UseThirdPartyPackage:
    
    def __init__(self, path_to_python: str, script: str) -> None:
        """Sets things up to get access to third-party packages. 

           This class is designed so that only synchronous communication is
               possible between the parent and child.
           It is expected that the child process reads data from stdin, processes
               the data from stdin according to its type, passes the result back 
               to the parent via stdout, and continues this loop until the parent 
               process kills it. This is a message passing scheme.
           Since raw data is exchanged between the parent and child, it's always
               necessary to specify the data type in any exchange. A general way
               to solve this is to specify an initial communication interaction
               pattern wherein the types are shared between parent and child. The
               quicker and easier way to solve this is restrict the message passing
               pattern to a single simple, hard-coded scheme. This class takes
               the latter approach. 
           To make things simple, this class ensures that a message sent from the
               parent process to a child process occupies two lines. The first 
               line describes the type of the message and the second 
               line contains the message content. Since the message content
               is restricted to a single line, it might be necessary to pass a
               path to a file to exchange lots of data. This class assumes that
               the child process responds to the parent with a two-line message
               in the same format as the sent message (i.e. a (message type,
               message data) tuple).
           
           TODO: Really, should be using Python's asyncio functionality so that
               the receive_data() function doesn't need a sleeep loop.
                       
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
        #     for modules.
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
                              stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                              text=True)
        self.sp = sp

        # DEBUG
        dp("The child process was just created! It has pid: " + str(self.sp.pid))


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


    def exchange_data(self, message_type: Any, message_data: str, timeout=120) -> tuple[str, str]:
        """Passes message to the child process, waits for a response, and returns
               the response.

           Assumes that the child process has not terminated.

           Args:
               message_type: The type of the message to pass to the child process.
                                 Should be an enum with string value.
               message_data: The data to pass to the child process. A single line.
               timeout:      The maximum amount of time, in seconds, to wait for
                                 the child process to send data back to the
                                 parent process.
        
           Returns:
               Tuple containing the type of message received and the message
                   content received.
        
           Raises:
               AssertionError if the timeout is exceeded.
        """
        
        assert self.sp is not None
        self.sp.poll()
        assert self.sp.returncode is None 

        assert len(message_type.value.splitlines()) == 1
        assert len(message_data.splitlines()) == 1
        
        self.sp.stdin.writelines([message_type.value, message_data])
        
        start_time = time.time()

        received = self.sp.stdout.readlines(2)
        
        while len(received) != 2:
            cur_time = time.time()
            if cur_time - start_time > timeout:
                raise AssertionError("Timeout exceeded!")
            time.sleep(3)            

            received = self.sp.stdout.readlines(2)

        return received[0], received[1] 
