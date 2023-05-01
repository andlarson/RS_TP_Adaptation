import io
import sys
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
np.set_printoptions(threshold=sys.maxsize)





"""
Input:
    file_name: String. Full path of .mtx file which contains stiffness matrix 
                 of interest.
    verbose: Boolean. Turns on verbose output. 

Functionality:
    Parses the file line-by-line and stores the matrix in a sparse format.

Return:
    Stiffness matrix.
    A sparse array in csr format. This is a scipy format one should be careful
      to avoid using numpy functions with it. Numpy doesn't understand anything
      about how sparse structures are stored!
    Note that the matrix entries listed in a .mtx file are 1-indexed so an 
      adjustment was made to 0-indexing.
"""
def parse_stiff_mat(file_name, verbose=False):

    if verbose:
        print("Starting to parse the stiffness matrix...")

    # Opened in binary mode!
    matrix_size = 0
    with open(file_name, 'rb') as f:

        # The stiffness matrix is a sparse square matrix.
        # To find its dimensionality, we need the largest index into it.
        # This largest index seems to fall on the last line of the file.

        # Move to the last byte before the '\n' at the end of the file.
        # Note use of -2 since each char occupies 2 bytes.
        last_byte_idx = f.seek(-2, io.SEEK_END)

        # Move backwards until a newline character is found.
        while f.read(1) != b'\n':
            f.seek(-2, io.SEEK_CUR)

        last_line = f.readline().decode()

        last_line = last_line.split(" ")

        matrix_size = int(last_line[0])

    assert(matrix_size != 0)
    
    if verbose:
        print("It appears that the matrix is square and has size " + str(matrix_size))

    # Build a sparse array via scipy.
    # Note that the format is lil_array.
    # The lil_array type is good for constructing sparse arrays.
    stiffness_matrix = sp.sparse.lil_array((matrix_size, matrix_size))

    with open(file_name, 'r') as f:
        for line in f:
            entries = line.split()

            # The lil_array type supports basic slicing.
            stiffness_matrix[int(entries[0]) - 1, int(entries[1]) - 1] = float(entries[2])

    # Now that we're done building the sparse array, convert it to a sparse
    #   format which is better for arithmetic.
    # Converting to csr is more efficient than to csc because it is row based,
    #   like lil.
    stiffness_matrix = stiffness_matrix.tocsr()

    if verbose:
        print("Done parsing the sparse matrix. It's now stored in a sparse format.")

    return stiffness_matrix





"""
Input:
    file_name: String. Full path of a .dat file which contains the nodal 
                 displacements of interest.
    verbose: Boolean. Dump information while the displacement vector is parsed.

Functionality:
    Parses the nodal displacements. Assumes there are exactly three
      displacments per node.
    The displacements are in the U1, U2, and U3 directions. It's necessary to
      understand what these directions represent!
    Stores the nodal displacements in a sparse array format. Note that the
      displacment vector is not actually sparse, but we don't want non-sparse
      and sprase formats mingling together.

Return:
    A sparse array in scipy csr format. It is 1 x N and lists the nodal 
      displacements in a format like: [node1_U1, node1_U2, node1_U3, node2_U1, 
      node2_U2, node2_U3, ...]
"""
def parse_displacement_vec(file_name, verbose=False):

    if verbose:
        print("Starting to parse the displacement vector...")

    # Figure out how many nodes are in the model. 
    node_cnt = -1 
    with open(file_name) as f:
        for line in f:
            if 'number of nodes' in line:
                node_cnt = int(line.split()[3])
    assert(node_cnt > 0)

    if verbose:
        print("Counted " + str(node_cnt) + " in the model.")

    # Create the displacement vector.
    # Use lil_array format for construction.
    displacement_vec = sp.sparse.lil_array((1, 3 * node_cnt))

    with open(file_name) as f:

        # Find the first instance of U1.
        # If there are many instances of U1 this will cause problems.
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

    if verbose:
        print("Done parsing the displacements of the nodes.")

    # Convert to csr format, which is better for computation.
    displacement_vec = displacement_vec.tocsr()

    return displacement_vec





"""
******************* Need to Update To Account For Sparse Format! **********************
Input:
    node - The node number of interest. Assumed to be 1-indexed, because Abaqus
      natively does 1-indexing.
    dir - The direction of interest. Note that x = 0, y = 1, z = 2.
    stiff_mat - The stiffness matrix. A 2D ndarray.
    disp_vec - The displacement vector. A 2D ndarray.

Functionality:
    Debugging function which can make it clear why the force for a particular
      node in a particular direction is what it is.
    For example, if the force being claculated at node 100 in the x direction
      appears wonky, this function can dump some information to shed light on
      the wonkiness. The force at node 100 in the x direction is exactly the
      dot product between a row in the stiffness matrix and the displacement
      vector.

Return:
    No return. This function prints an information dump.
"""
def debug_force(node_of_interest, dir_of_interest, stiff_mat, disp_vec):
    
    dir_of_interest_str = ""
    if dir_of_interest == 0:
        dir_of_interest_str = "x"
    elif dir_of_interest == 1:
        dir_of_interest_str = "y"
    else:
        dir_of_interest_str = "z"


    # Map the node + direction information to the correct row of the
    #   stiffness matrix and the computed force vector.
    r = (node_of_interest - 1) * 3 + dir_of_interest

    # Dump the nodal displacements (node numbers and directions) which affect
    #    the force being computed for some node in some direction.
    for c in range(stiff_mat.shape[1]):
        
        # For each nonzero entry in that row, dump which node and which direction
        #   affects the force calculation.
        if stiff_mat[r, c] != 0:

            # Add 1 because node numbering in Abaqus starts at 1.
            node_affecting_force = int(c/3) + 1

            # Figure out the direction of displacement affecting the force
            #   calculation.
            direction = c % 3
            dir_str = ""
            if direction == 0:
                dir_str = "x"
            elif direction == 1:
                dir_str = "y"
            else:
                dir_str = "z"

            print("Displacement at node " + str(node_affecting_force) + " in the " + dir_str + " direction affects the force computed for node " + str(node_of_interest) + " in the direction " + str(dir_of_interest_str))
            print("The displacement at node " + str(node_affecting_force) + " in the " + dir_str + " is " + str(disp_vec[0, c]))
            print("The corresponding coefficient in the stiffness matrix is: " + str(stiff_mat[r, c]))
            print("")

    # Print the force in the direction of interest for the node of interest.
    print("The computed force for node " + str(node_of_interest) + " in the direction " + str(dir_of_interest_str) + " is " + str((stiff_mat @ disp_vec.T)[r, 0]))






if __name__ == '__main__':

    dir_of_interest = "./../../shot_peened_bar_inp/"

    stiff_mat = parse_stiff_mat(dir_of_interest + "CutShotPeen_STIF2.mtx", verbose=True)
    disp_vec = parse_displacement_vec(dir_of_interest + "CutShotPeen.dat", verbose=True)

    # Do the matrix-vector multiplication.
    # The result of this computation is a vector which has 3 consecutive entries
    #   for each node. The 3 entries per node describe the force in the x, y and
    #   z directions which resulted in the observed displacements.
    # Remember that nodes are labeled 1, ...., N but, for example, the forces for
    #   node 1 are located in entries [0, 2].
    print(stiff_mat.dot(disp_vec))

    # Contour of the stiffness matrix.
    # plt.imshow(stiff_mat, cmap='hot', interpolation='nearest')
    # plt.show()
    
    
    







