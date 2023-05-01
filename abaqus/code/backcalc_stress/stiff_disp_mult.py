import io
import sys
import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
from pathlib import Path



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
        raise RuntimeError("Path passed to save_stiff_mat does not exist.")

    # Read the matrix in.
    # Not clear if the resulting matrix is in a sparse format or not.
    #   Depends on the data in the .mtx file.
    stiffness_mtx = sp.io.mmread(path)

    # Force it into sparse format.
    # We'll use the lil foramt because it is efficient for construction.
    stiffness_mtx = sp.sparse.lil_array(stiffness_mtx)

    # Now save it into the stiffness matrix storage area.
    sp.sparse.save_npz("./stiffness_matrices/" + name)
            


"""
Input:
    path - String. 
           Full path of saved sitffness matrix which is in sparse format.

Functionality:
    Reads the sparse matrix from .npz format into a sparse array and returns
      the sparse array.

Return:
    A sparse array in csr (scipy) format.
"""
def load_stiff_mat(path):

    # Check that the file actaully exists!
    # If not, let the user know.
    if not Path(path).exists():
        raise RuntimeError("Path passed to load_stiff_mat does not exist.")

    stiffness_mat = sp.sparse.csr_array(load_npz(path))

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

    path_to_mtx = '/home/andrew/Desktop/stuff/umich/umich_machining_project/abaqus/fine_mesh_bar_inp/Job-1_STIF1.mtx'
    save_stiff_mat(path_to_mtx, 'fine_mesh_bar') 
    stiff_mat = load_stiff_mat('fine_mesh_bar')

    path_to_disp_vec = '/home/andrew/Desktop/stuff/umich/umich_machining_project/abaqus/fine_mesh_bar_inp/Job-1.dat'
    disp_vec = parse_displacement_vec(path_to_disp_vec)

    # Do the matrix-vector multiplication.
    print(stiff_mat.dot(disp_vec))

    
