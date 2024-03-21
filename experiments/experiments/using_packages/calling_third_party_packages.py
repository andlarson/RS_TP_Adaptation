"""
This file attempts to cal third party packages by using utilities provided by
    this library.
"""


import sys

# Resolving imports.
sys.path.append("/home/andlars/Desktop/RS_TP_Adaptation/")

import src.core.abaqus.abaqus_shim as shim
import src.util.third_party_packages.blender_messages as bm 
import src.util.third_party_packages.parent_process.calling_third_parties as ctp

from src.util.debug import *



if __name__ == "__main__":

    # Do some stuff using the Abaqus kernel.
    PATH_TO_MDB = "/home/andlars/Desktop/RS_TP_Adaptation/experiments/experiments/using_packages/traction_vectors.cae"
    mdb = shim.use_mdb(PATH_TO_MDB)
    shim.close_mdb(mdb)

    # Now spin up a child process and use some functionality from a third party 
    #     package. 
    PATH_TO_STL = "/home/andlars/Desktop/RS_TP_Adaptation/experiments/experiments/using_packages/voyager_dish.stl" 
    PATH_TO_SAVE_LOC = "/home/andlars/Desktop/RS_TP_Adaptation/experiments/experiments/using_packages/test.stl"
    remesh_message = bm.Messages.BLENDER_REMESH
    remesh_message_data = bm.parent_to_child_remesh_data(PATH_TO_STL, PATH_TO_SAVE_LOC)
    third_party_package = ctp.UseThirdPartyPackage(ctp.PATH_TO_STANDARD_INTERPRETER, ctp.PATH_TO_STANDARD_CHILD_SCRIPT)
    third_party_package.start_child()
    received_message_type, received_message_data = third_party_package.exchange_data(remesh_message, remesh_message_data)
    dp("The received message has type " + received_message_type)
    bm.unpack_remesh_data_at_parent(received_message_data)
    third_party_package.kill_child()
    dp("End of test!")
    


