import os

MDB_DEFAULT_DIR = "MDB_DEFAULT"


# Get the path to the default MDB save area. 
def MDB_save_area():
# type: (None) -> str

    path_to_this_file = os.path.abspath(__file__)
    
    # Manipulate the path so that it points to the MDB save area.
    # This assumes a relationship between the location of this file and the
    #   location of the MDB save area.
    dir_of_this_file = os.path.dirname(path_to_this_file)
    mdb_default_path = dir_of_this_file + MDB_DEFAULT_DIR 

    if not os.path.exists(mdb_default_path):
        raise RuntimeError("Can't find the default save area for MDB files.") 

    return mdb_default_path


