"""
Testing the PartFromMesh() function available through Abaqus' Python Scripting Interface. 
"""

from abaqus import *
from abaqusConstants import *

if __name__ == "__main__":

    mdb = openMdb("test.cae")

    new_part = mdb.models["Model-1"].parts["Part-1"].PartFromMesh("from_mesh")

    mdb.saveAs("post_test.cae")

    


