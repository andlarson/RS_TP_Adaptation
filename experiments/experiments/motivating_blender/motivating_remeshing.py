"""
This script motivates remeshing by attempting to convert a mesh to geometry.
    I want to find a mesh which fails in the conversion attempt (using this script), 
    and then figure out what characteristics of that mesh cause it to fail (using
    something like Blender).
"""

# Resolving imports.
sys.path.append("/home/andlars/Desktop/RS_TP_Adaptation")

import src.core.simulation.simulation as sim

from src.util.debug import *


if __name__ == "__main__":
    
    ODB_PATH = "/home/andlars/Desktop/RS_TP_Adaptation/experiments/experiments/mesh_produced_by_abaqus.odb" 
    NEW_MDB_NAME = "mesh_to_geom_fail.cae"
    NEW_MDB_PATH = "/home/andlars/Desktop/RS_TP_Adaptation/experiments/experiments/motivating_remeshing"
    sim.create_mdb_from_odb(NEW_MDB_NAME, NEW_MDB_PATH, ODB_PATH)
