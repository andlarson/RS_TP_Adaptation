"""
Figuring out the checkGeometry() method.
"""

# Resolve imports.
sys.path.append("/home/andlars/Desktop/RS_TP_Adaptation/src")

import sys

from abaqus import *
from abaqusConstants import *

import core.abaqus.abaqus_shim as shim


if __name__ == "__main__":

    mdb = openMdb("test.cae")

    part = mdb.models["Model-1"].parts["Part-1"]

    shim.check_shell_geometry(part)

    mdb.saveAs("post_test.cae")
    mdb.close()

