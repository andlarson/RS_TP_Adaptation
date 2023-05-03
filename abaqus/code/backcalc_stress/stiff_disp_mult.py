import io
import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
from pathlib import Path
import subprocess


"""
Note the Abaqus convention:
    Direction 1 = X Direction
    Direction 2 = Y Direction
    Direction 3 = Z Direction
"""



"""
Input:
    path - String.
           Full path of file.

Functionality:
    Counts the total number of lines in a file via the wc program. Note that
      this is a Unix program, so it's expected that the system has it available.
    It's expected that this is faster than anything Python can do.

Return:
    Int. Number of lines in the file.
"""
def file_len(path):
    p = subprocess.Popen(['wc', '-l', path], stdout=subprocess.PIPE, 
                                             stderr=subprocess.PIPE)
    result, err = p.communicate()
    if p.returncode != 0:
        raise IOError(err)
    
    return int(result.strip().split()[0])



"""
Input:
    path - String.
           Full path of file.

Functionality:
    Returns the last line of the file. Avoids reading all lines into memory so
      it works well for very large files.

Return:
    String. Last line of file.
"""
def read_last_line(path):

    with open(path, 'rb') as f:
        try: 
            f.seek(-2, os.SEEK_END)
            while f.read(1) != b'\n':
                f.seek(-2, os.SEEK_CUR)
        except OSError:
            f.seek(0)

        last_line = f.readline().decode()

    return last_line



"""
Input:
    path - String.
           Full path of .mtx file.

Functionality:
    Parses the last line of the .mtx file uses it to return the size.
    It's assumed that the node number in the last line of the .mtx file reflects
      the matrix size.

Return:
    Int. The size of the stiffness matrix.
"""
def get_mtx_size(path):
    
    last_line = read_last_line(path)

    parts = last_line.split()
    
    assert(len(parts) == 3)

    return int(parts[0])



"""
Input:
    path    - String. 
              Path of .mtx file.
    size    - Integer. 
              Matrix assumed to be square.
    nonzero - Integer. 
              Number of nonzero entries in matrix.

Functionality:
    Adds the necessary metadata to the top of the .mtx file so that the file is
      in true .mtx format and can be interpreted by scipy.

Return:
    None.
"""
def add_mtx_metadata(path, size, nonzero):
    return 0



"""
CURRENTLY NOT IN USE 

Input:
    path - String.
           Full path of .mtx file of interest.

Functionality:
    Abaqus produces .mtx files which don't exactly conform with the .mtx
      standard which scipy expects. In particular, the .mtx files produced
      by Abaqus are missing the first two lines of metadata.
    See the discussions at https://math.nist.gov/MatrixMarket/formats.html
      and https://stackoverflow.com/q/52947403 for more info.
    As such, this function modifies a .mtx file to include this metadata.
    This function assumes that the number of lines in the .mtx file is exactly
      the number of nonzero entries in the matrix. It also assumes that the
      last line of the .mtx file contains the size of the matrix.

Return:
    None.
"""
def make_abaqus_mtx_conform(path):

    if not Path(path).exists():
        raise RuntimeError("Bad path! Path does not exist.")

    line_cnt = file_len(path)
    last_line = read_last_line(path)

    # Now write the necessary metadata.
    mtx_file = open(path, 'r+')
    
    # The generic stuff.

 

"""
Input:
    path    - String.
              Full path of .mtx file which contains stiffness matrix of interest.
    name    - String.
              The desired name of the file where the sparse matrix will be saved.
    verbose - Boolean.
              Turns on verbose output.

Functionality:
    Reads the .mtx file into csr_array format and then saves it at a particular
      location in the .npz format.
    Requires that the matrix is truly in .mtx format so that scipy has the
      ability to read it.
    
Return:
    None.
"""
def save_stiff_mat(path, name, verbose):

    if not Path(path).exists():
        raise RuntimeError("Bad path! Path does not exist!")

    if verbose:
        print("Starting to read the mtx file!")

    # Read the matrix in.
    # Not clear if the resulting matrix is in a sparse format or not.
    #   Depends on the data in the .mtx file.
    stiffness_mtx = sp.io.mmread(path)

    # Force it into sparse format.
    # We'll use the lil foramt because it is efficient for construction.
    stiffness_mtx = sp.sparse.lil_array(stiffness_mtx)

    # Now save it into the stiffness matrix storage area.
    sp.sparse.save_npz("./stiffness_matrices/" + name, stiffness_mtx.tocsr())

    if verbose:
        print("Done reading the mtx file!")



"""
Input:
    path     - String. 
               Path to .mtx file.
    save_dir - String.
               Directory to save the stiffness matrix in.
    name     - String.
               Desired name of .npz file.
    node_cnt - Int.
               The number of nodes in the model's mesh.
    verbose  - Boolean.
               Turns on verbose output.

Functionality:
    Parses the .mtx file line-by-line and builds a sparse matrix. As the matrix
      is built, it is in lil_array (scipy) format. It is stored as a .npz file
      in csr_array (scipy) format. The matrix is always saved to usual location.
    This function is offered as an alternative which does the parsing manually.
      This is useful when the file containing the .mtx file is huge and it
      doesn't make sense to add metadata to it to make it an official .mtx file.

Return:
    None.
"""
def parse_and_save_stiff_mat(path, save_dir, name, node_cnt, verbose):
    
    if not Path(path).exists():
        raise RuntimeError("Bad path! Path does not exist!")

    if verbose:
        print("Starting to parse stiffness matrix")

    stiffness_mtx = sp.sparse.lil_array((3 * node_cnt, node_cnt))

    # Read the file line-by-line without bringing the whole thing into memory.
    with open(path) as infile:
        for line in infile:
            
            # Split the line into its component parts.
            parts = line.split()

            # We expect each line of the .mtx file to have 3 components.
            assert(len(parts) == 3)

            # Subtract 1 because the nodes in Abaqus are 0-indexed.
            i = int(parts[0]) - 1
            j = int(parts[1]) - 1
            val = float(parts[2])

            stiffness_mtx[i, j] = val   

    # Now save it into the stiffness matrix storage area.
    sp.sparse.save_npz(save_dir + name, stiffness_mtx.tocsr())

    if verbose:
        print("Done parsing and saving stiffness matrix")



"""
Input:
    path     - String.
               Full path of a .dat file which contains the nodal displacements
                 of interest.
    save_dir - String.
               Directory to save the displacement vector to.
    name     - String.
               Desired name of saved file.
    node_cnt - Int.
               The number of nodes in the model's mesh.
    verbose  - Boolean.
               Turns on verbose output.

Functionality:
    WARNING: This function uses 'U1' to find the beginning of nodal
      displacments of interest, uses an empty line to find the end
      of the nodal displacments. It implicitly assumes that there is
      only a single collection of nodal displacments in the .dat file.
    Parses the nodal displacements. Assumes there are exactly three
      displacments per node.
    The displacements are in the U1, U2, and U3 directions. It's necessary to
      understand what these directions represent!
    Saves a 2D numpy array (1 by N) which lists the nodal displacements in a 
      format like: [node1_U1, node1_U2, node1_U3, node2_U1, node2_U2, node2_U3, 
      ...].

Return:
    None.
"""
def parse_and_save_displacement_vec(path, save_dir, name, node_cnt, verbose):

    if not Path(path).exists():
        raise RuntimeError("Bad path passed to parse_displacement_vec!") 

    if verbose:
        print("Starting to parse the displacement vector!")

    # Create the displacement vector.
    displacement_vec = np.zeros((1, 3 * node_cnt))

    with open(path) as f:

        # Use 'U1' to find the beginning of the nodal displacements.
        while ('U1' not in f.readline()):
            pass
    
        # Hard-coded move to beginning of data of interest.
        f.readline()
        f.readline()

        # Populate the displacement vector!
        cur_line = f.readline()
        while cur_line.strip() != '': 
            nodal_info = cur_line.split()
            idx = (int(nodal_info[0]) - 1) * 3

            displacement_vec[0, idx] = float(nodal_info[1])
            displacement_vec[0, idx + 1] = float(nodal_info[2])
            displacement_vec[0, idx + 2] = float(nodal_info[3])

            cur_line = f.readline()

    np.save(save_dir + name, displacement_vec)

    if verbose:
        print("Done parsing and saving the displacement vector!")



"""
Input:
    path - String. 
           Full path of stored array.
           
Functionality:
    Reads the stored array (can be in .npz or .npy format) and returns it.
      Note that the result can be a sparse array or a dense array.

Return:
    An array of type corresponding to the saved type.
"""
def load_saved_array(path):

    # Check that the file actaully exists!
    # If not, let the user know.
    if not Path(path).exists():
        raise RuntimeError("Bad path! Can't find the specified saved array.")

    # Detect if we're loading a scipy sparse array or a standard numpy dense array.
    substrings = path.split('.')

    # Assuming that the name ends with either '.npy' or '.npz'.
    if substrings[-1] == 'npy':
        saved_array = np.load(path)
    elif substrings[-1] == 'npz':
       saved_array = sp.sparse.load_npz(path).tocsr()
    else:
        raise RuntimeError("Parsing of path to stored array failed!")

    return saved_array 



if __name__ == '__main__':

    # Metadata helps avoid too much parsing.
    NODE_CNT = 148470

    matrix_storage_area = './stiffness_matrices/'
    disp_vec_storage_area = './displacement_vectors/'

    path_to_mtx = '../../shot_peened_bar_inp/CutShotPeen_STIF2.mtx'
    # parse_and_save_stiff_mat(path_to_mtx, matrix_storage_area, 'shot_peen_bar_mtx', NODE_CNT, True)
    stiff_mat = load_saved_array(matrix_storage_area + 'shot_peen_bar_mtx.npz')

    path_to_disp_vec = '../../shot_peened_bar_inp/CutShotPeen.dat'
    # parse_and_save_displacement_vec(path_to_disp_vec, disp_vec_storage_area, 'shot_peen_bar_dispv', NODE_CNT, True)
    disp_vec = load_saved_array(disp_vec_storage_area + 'shot_peen_bar_dispv.npy').T

    # Do the matrix-vector multiplication.
    res = stiff_mat.dot(disp_vec)
    print(res) 

    # DEBUG
    print(res[0:100, 0])

    # Save off the result for future use.
    np.save('./experimental_results/shot_peen_bar_force_vec', res)



