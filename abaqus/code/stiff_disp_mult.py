import io
import sys
import numpy as np
import matplotlib.pyplot as plt
np.set_printoptions(threshold=sys.maxsize)



"""
Input:
    Full path of .mtx file which contains stiffness matrix of interest.

Functionality:
    Parses the file line-by-line and stores the matrix in a sparse format.

Return:
    Stiffness matrix.
    A 2D numpy array of arbitrary size. Note that the matrix entries listed
      in a .mtx file are 1-indexed so an adjustment was made to 0-indexing.
"""
def parse_stiff_mat(file_name):

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

    stiffness_matrix = np.eye(matrix_size)

    with open(file_name, 'r') as f:
        
        for line in f:
            entries = line.split()
            stiffness_matrix[int(entries[0]) - 1, int(entries[1]) - 1] = float(entries[2])

    return stiffness_matrix
            



"""
Input:
    Full path of a .dat file which contains the nodal displacements of interest.

Functionality:
    Parses the nodal displacements. Assumes there are exactly three
      displacments per node.
    The displacements are in the U1, U2, and U3 directions. It's necessary to
      understand what these directions represent!

Return:
    A 2D numpy array (1 by N) which lists the nodal displacements in a format 
      like: [node1_U1, node1_U2, node1_U3, node2_U1, node2_U2, node2_U3, ...]
"""
def parse_displacement_vec(file_name):

    # Figure out how many nodes are in the model. 
    node_cnt = -1 
    with open(file_name) as f:
        for line in f:
            if 'number of nodes' in line:
                node_cnt = int(line.split()[3])
    assert(node_cnt > 0)

    # Create the displacement vector.
    displacement_vec = np.zeros((1, 3 * node_cnt))

    with open(file_name) as f:

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
    stiff_mat = parse_stiff_mat('/home/andrew/Desktop/stuff/umich/umich_machining_project/abaqus/fine_mesh_bar_inp/Job-1_STIF1.mtx')
    disp_vec = parse_displacement_vec('/home/andrew/Desktop/stuff/umich/umich_machining_project/abaqus/fine_mesh_bar_inp/Job-1.dat')

    # Do the matrix-vector multiplication.
    # The result of this computation is a vector which has 3 consecutive entries
    #   for each node. The 3 entries per node describe the force in the x, y and
    #   z directions which resulted in the observed displacements.
    # Remember that nodes are labeled 1, ...., N but, for example, the forces for
    #   node 1 are located in entries [0, 2].
    print(stiff_mat @ disp_vec.T)
    # print((stiff_mat @ disp_vec.T)[1212 * 3: 1212 * 3 + 3])

    # Contour of the stiffness matrix.
    plt.imshow(stiff_mat, cmap='hot', interpolation='nearest')
    plt.show()
    
    
    # -------------------
    # A setup to examine / debug why the force for a particular node in a
    #   particular direction is as it is.
    # -------------------


    """
    # The node we are intested in.
    node_of_interest = 1460 

    # Each node has a force computed for it in each of its degrees of
    #   freedom. Specify which degree of freedom (direction) that we are
    #   interested in.
    # x = 0, y = 1, z = 2
    dir_of_interest = 2

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
    """ 








