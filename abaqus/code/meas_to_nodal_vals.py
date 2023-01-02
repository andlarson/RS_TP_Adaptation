
"""
Input:
    A .inp file defining an Abaqus simulation.
Functionality:
    Parses the portion of the .inp file which lists the node numbers and
      nodal coordinates.
Return:
    Node number to coordinate mappings. 
    List of tuples. Each tuple contains: node number, x coord, y coord, z
      coord.
"""
def dotinp_parse_nodes(file_name):

    nodal_coords = []

    # Section starts with '*Node' on its own line and ends with a line which
    #    starts with '*Element'.
    # This may change... 
    in_sec = False
    with open(file_name, 'r') as f:
        for line in f:
            if r"*Node" in line:
                in_sec = True
            elif r"*Element" in line:
                in_sec = False
            elif in_sec:
                line = line.replace(" ", "")
                line = line.replace(".", "")
                nums = line.split(",")
                nums = [int(num) for num in nums] 

                assert len(nums) == 4

                nodal_coords.append((nums[0], nums[1], nums[2], nums[3]))

    return nodal_coords




"""
Input: 
    Node number to coordinate mappings. 
Functionality: 
    Collects groups of nodes with the same y-coordinates.
Return:
    Groups of nodes at same y-coordinate.
    List of tuples. Each tuple contains a y-coordinate and a list of nodes
      which lie at that y-coordinate.
"""
def nodal_coords_to_groups(nodal_coords):
    y_grouped_nodes = []

    for node in nodal_coords:
        has_group = False
        for group in y_grouped_nodes:
            # Assumption made about y-coordinate listing number in the .inp
            #   file.
            if group[0] == node[3]:
                has_group = True
                group[1].append(node[0])

        if not(has_group):
            y_grouped_nodes.append((node[3], [node[0]]))

    return y_grouped_nodes



     
# """
# Input: 
#     List of tuples. Each tuple contains a y-coordinate and a list of nodes
#       which lie at that y-coordinate.
# Functionality: 
#     Assigns displacements to each node according to its y-coordinate and then
#       creates the displacement vector. Each node has a displacement scalar
#       for each dimension.
# Return:
#     Vector of nodal displacements. 
#     Represented by list of nodal tuples. One tuple for each node. Each tuple 
#       contains: (node number, displacement tuple). The displacement tuple 
#       contains: (x disp, y disp, z disp).
# """
# def groups_to_disps(y_coord_groups):
# 
# 
# 
# 
# 
# """
# Input:
#     Vector of nodal displacements.
# Functionality:
#     Performs the multiplication between the sparse stiffness matrix and the
#       displacement vector.
# Return:
#     Vector of nodal stresses.
#     List of tuples. Each tuple contains a node number and a scalar stress
#       value. The units of the stress values are !!!!!!!!!!!!!!.
# """
# def disp_stiff_mult(nodal_disp_vec):
# 
# 
# 
# 
# """
# Input:
#     Vector of nodal stresses.
#     Node number to coordinate mappings.
# Functionality:
#     Plots the nodal stresses throughout the part body. For now, no assumption
#       is made about the stress symmetry within the body.
# Return:
#     None.
# """
# def plot_stresses_in_3d(nodal_stress_vec, nodal_coords):
# 
# 
# 
# 
# """ TODO: Function to plot stresses varying y alone."""




if __name__ == "__main__":
    nodal_coords = dotinp_parse_nodes("/home/andrew/Desktop/stuff/umich/umich_machining_project/abaqus/fine_mesh_bar/Job-1.inp")
    print(nodal_coords_to_groups(nodal_coords))    



