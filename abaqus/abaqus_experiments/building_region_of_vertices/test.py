"""
This file tests building a region out of vertices.
"""

import sys

import numpy as np

# Resolve imports.
sys.path.append("/home/andlars/Desktop/RS_TP_Adaptation/")

import src.core.abaqus.abaqus_shim as shim
import src.util.geom as geom



if __name__ == "__main__":

    mdb = shim.use_mdb("/home/andlars/Desktop/RS_TP_Adaptation/experiments/experiments/building_region_of_vertices/region_of_vertices.cae")

    part_instance = mdb.models[shim.STANDARD_MODEL_NAME].rootAssembly.allInstances["Initial_Geometry-1"]

    p1 = geom.Point3D(np.array([0, 10, 0]))
    p2 = geom.Point3D(np.array([0, 0, 0]))

    shim.build_region_with_vertices([p1, p2], part_instance)

    shim.save_mdb_as("/home/andlars/Desktop/RS_TP_Adaptation/experiments/experiments/building_region_of_vertices/result.cae", mdb)
    shim.close_mdb(mdb)
