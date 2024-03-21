"""
Testing how to manipulate and redirect stdin, stdout, stderr, etc.
"""

import sys

from abaqus import *
from abaqusConstants import *


if __name__ == "__main__":

    mdb = openMdb("test.cae")

    # These are printed to the terminal.
    sys.__stderr__.write("Test 1\n")
    sys.__stdout__.write("Test 2\n")

    # This is printed to the in-Abaqus GUI.
    # Evidently Abaqus redirects stderr and stdout to a temporary file that the GUI gets
    #    access to.
    sys.stderr.write("Test 3\n")
    sys.stdout.write("Test 4\n")

    mdb.saveAs("post_test.cae")
    mdb.close()

