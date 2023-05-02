import io
import sys
import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
from pathlib import Path
import subprocess



"""
Input:
    path - String.
           Full path of file.

Functionality:
    Counts the total number of lines in a file via the wc program. Note that
      this is a Unix program, so it's expected that the system has it available.
    It's expected that this is faster than anything Python can do.

Return:
    None.
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
    path - String. Path of .mtx file.
    size - Integer. Matrix assumed to be square.
    nonzero - Integer. Number of nonzero entries in matrix.

Functionality:
    Adds the necessary metadata to the top of the .mtx file so that the file is
      in true .mtx format and can be interpreted by scipy.

Return:
    None.
"""
def add_mtx_metadata(path, size, nonzero):

    return 0



"""
CURRENTLY NOT IN USE - Right now, the technique I'm using is to modify the
.mtx file directly to make it conform.

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
    path - String.
           Full path of .mtx file which contains stiffness matrix of interest.
    name - String.
           The desired name of the file where the sparse matrix will be saved.

Functionality:
    Reads the .mtx file into a sparse format and then saves it at a particular
      location in the .npz format.
    I believe that the numpy routines will not understand the sparse format
      the matrix is saved in, so loading it should be done via the appropriate
      scipy routines which recognize sparse matrix routines.

Return:
    None.
"""
def save_stiff_mat(path, name):

    if not Path(path).exists():
        raise RuntimeError("Bad path! Path does not exist!")

    # Read the matrix in.
    # Not clear if the resulting matrix is in a sparse format or not.
    #   Depends on the data in the .mtx file.
    stiffness_mtx = sp.io.mmread(path)

    # Force it into sparse format.
    # We'll use the lil foramt because it is efficient for construction.
    stiffness_mtx = sp.sparse.lil_array(stiffness_mtx)

    # Now save it into the stiffness matrix storage area.
    sp.sparse.save_npz("./stiffness_matrices/" + name, stiffness_mtx.tocsr())
            


"""
Input:
    name - String. 
           Name of saved sparse matrix. It is expected that the matrix is
             stored in the directory of sparse matrices.

Functionality:
    Reads the sparse matrix from .npz format into a sparse array and returns
      the sparse array.

Return:
    A sparse array in csr (scipy) format.
"""
def load_stiff_mat(name):

    # Construct the path, assuming that the matrix is stored in the saved
    #   matrix directory.
    path = './stiffness_matrices/' + name + '.npz'

    # Check that the file actaully exists!
    # If not, let the user know.
    if not Path(path).exists():
        raise RuntimeError("Bad name! Name does not exist in the matrix save area!")

    stiffness_mat = sp.sparse.csr_array(sp.sparse.load_npz(path))

    return stiffness_mat 



"""
Input:
    path - String.
           Full path of a .dat file which contains the nodal displacements
             of interest.

Functionality:
    Parses the nodal displacements. Assumes there are exactly three
      displacments per node.
    The displacements are in the U1, U2, and U3 directions. It's necessary to
      understand what these directions represent!

Return:
    A 2D numpy array (1 by N) which lists the nodal displacements in a format 
      like: [node1_U1, node1_U2, node1_U3, node2_U1, node2_U2, node2_U3, ...]
"""
def parse_displacement_vec(path):

    if not Path(path).exists():
        raise RuntimeError("Bad path passed to parse_displacement_vec!") 

    # Figure out how many nodes are in the model. 
    node_cnt = -1 
    with open(path) as f:
        for line in f:
            if 'number of nodes' in line:
                node_cnt = int(line.split()[3])
    assert(node_cnt > 0)

    # Create the displacement vector.
    displacement_vec = np.zeros((1, 3 * node_cnt))

    with open(path) as f:

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

    return displacement_vec
            



if __name__ == '__main__':

    path_to_mtx = '../../fine_mesh_bar_inp/Job-1_STIF1_modified.mtx'
    save_stiff_mat(path_to_mtx, 'fine_mesh_bar') 
    stiff_mat = load_stiff_mat('fine_mesh_bar')

    path_to_disp_vec = '/home/andrew/Desktop/stuff/umich/umich_machining_project/abaqus/fine_mesh_bar_inp/Job-1.dat'
    disp_vec = parse_displacement_vec(path_to_disp_vec).T

    # Do the matrix-vector multiplication.
    print(stiff_mat.dot(disp_vec)) 


