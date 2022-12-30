
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
def parse_nodes(file_name):




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




"""
Input: 
    List of tuples. Each tuple contains a y-coordinate and a list of nodes
      which lie at that y-coordinate.
Functionality: 
    Assigns displacements to each node according to its y-coordinate and then
      creates the displacement vector. Each node has a displacement scalar
      for each dimension.
Return:
    Vector of nodal displacements. 
    Represented by list of nodal tuples. One tuple for each node. Each tuple 
      contains: (node number, displacement tuple). The displacement tuple 
      contains: (x disp, y disp, z disp).
"""
def groups_to_disps(y_coord_groups):





"""
Input:
    Vector of nodal displacements.
Functionality:
    Performs the multiplication between the sparse stiffness matrix and the
      displacement vector.
Return:
    Vector of nodal stresses.
    List of tuples. Each tuple contains a node number and a scalar stress
      value. The units of the stress values are !!!!!!!!!!!!!!.
"""
def disp_stiff_mult(nodal_disp_vec):




"""
Input:
    Vector of nodal stresses.
    Node number to coordinate mappings.
Functionality:
    Plots the nodal stresses throughout the part body. For now, no assumption
      is made about the stress symmetry within the body.
Return:
    None.
"""
def plot_stresses_in_3d(nodal_stress_vec, nodal_coords):




""" TODO: Function to plot stresses varying y alone.


