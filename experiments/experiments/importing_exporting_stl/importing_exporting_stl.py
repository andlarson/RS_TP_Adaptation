"""
This file acts as a testbed for individual components of the library.
"""

import sys

# Resolving imports.
sys.path.append("/home/andlars/Desktop/RS_TP_Adaptation/")

import src.util.third_party_packages.parent_process.calling_third_parties as third_parties
import src.util.third_party_packages.blender_messages as blender_messages
import src.core.abaqus.abaqus_shim as shim
import src.core.simulation.simulation as sim

from src.util.debug import *


if __name__ == "__main__":

    # Test exporting a .stl from an MDB.
    """
    assembly_module = sim.ModuleNames.ASSEMBLY
    MODEL_NAME = shim.STANDARD_MODEL_NAME
    STL_PATH = "/home/andlars/Desktop/RS_TP_Adaptation/experiments/experiments/importing_exporting_stl/from_cae.stl"
    mdb = shim.use_mdb("/home/andlars/Desktop/RS_TP_Adaptation/experiments/experiments/traction_vectors_cae/traction_vectors.cae")
    sim._export_stl_from_mdb(assembly_module, MODEL_NAME, STL_PATH, mdb)
    shim.close_mdb(mdb)
    """

    # Test the conversion of .odb to .stl format.
    """
    ODB_PATH = "/home/andlars/Desktop/RS_TP_Adaptation/experiments/experiments/traction_vectors_cae/Job-1.odb"
    STL_PATH = "/home/andlars/Desktop/RS_TP_Adaptation/experiments/experiments/importing_exporting_stl/from_odb.stl"
    sim._convert_odb_to_stl(ODB_PATH, STL_PATH)
    """

    # Test importing complex .stl file.
    STL_PATH = "/home/andlars/Desktop/RS_TP_Adaptation/experiments/experiments/importing_exporting_stl/cool_ascii_stl/hubble_main_body.stl"
    MDB_NAME = "importing_stl.cae"
    MDB_SAVE_DIR = "/home/andlars/Desktop/RS_TP_Adaptation/experiments/experiments/importing_exporting_stl/"
    mdb = sim._import_stl_to_mdb(STL_PATH, MDB_NAME, MDB_SAVE_DIR)
    shim.save_mdb(mdb)
    shim.close_mdb(mdb)



   
    
