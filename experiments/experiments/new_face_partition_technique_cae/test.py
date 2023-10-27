"""
Testing new partitioning technique.
"""

# Resolve imports.
sys.path.append("/home/andlars/Desktop/RS_TP_Adaptation/src")

import numpy as np 
from abaqus import *
from abaqusConstants import *

import core.abaqus.abaqus_shim as shim
from util.debug import *
import util.geom as geom


if __name__ == "__main__":

    p1 = geom.Point3D(np.array((0, 10, 0)))
    p2 = geom.Point3D(np.array((40, 10, 0)))
    p3 = geom.Point3D(np.array((40, 10, 40)))
    p4 = geom.Point3D(np.array((0, 10, 40)))

    ngon = geom.NGon3D([p1, p2, p3, p4])

    mdb = openMdb("test.cae")

    part = mdb.models["Model-1"].parts["Part-1"]
    assembly = mdb.models["Model-1"].rootAssembly
    instance = assembly.instances["instance-1"]

    dp("Before partitioning, the instance had faces with the following vertices.")
    for face in instance.faces:
        dp(str(shim.get_face_vertices(face, assembly=assembly, part_instance=instance)))

    new_face = shim.partition_face(ngon, None, None, None, assembly=assembly, instance=instance) 

    dp("After partitioning, the instance had faces with the following vertices.")
    for face in instance.faces:
        dp(str(shim.get_face_vertices(face, assembly=assembly, part_instance=instance)))

    vertex_ids = new_face.getVertices()
    for id in vertex_ids:
        dp("The vertex on the face is: " + str(instance.vertices[id].pointOn))

    mdb.saveAs("post_test.cae")
    mdb.close()