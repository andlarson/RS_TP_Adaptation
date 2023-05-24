import numpy as np
import parse_sim_output as pso 
import matplotlib.pyplot as plt
import random


"""
Input:
    coordinates - ndarray.
                  The 3D coordinates to be plotted.
                  Assumed to be 2D ndarray with shape (N, 3).
    vals        - ndarray.
                  The values to be plotted at each point in space.
                  Assumed to be 1D array of shape (3*N).
    direction   - String.
                  Specifies direction of interest.
  
Functionality:
    WARNING: Some model specific modifications are made to the plot for better
      viewing.
    Creates a dynamic plot of the 3D data. The values are differentiated via
      their color.

Return:
    None.
"""
def plot_directional_forces(coordinates, vals, direction):

    assert(coordinates.shape[1] == 3)
    assert(vals.ndim == 1)

    vals = select_dir(direction, vals)

    fig = plt.figure(figsize=(13,9), constrained_layout=True)
    ax = fig.add_subplot(projection='3d')

    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')

    ax.set_title("Illustrating Forces in the " + direction + " Direction")

    # ---- Model Specific Mods ----

    # Stretch the z-axis as it is longest in the model.
    # This is model specific.
    ax.set_box_aspect(aspect = (1,1,2))

    # Rotate the plot by setting the viewing angle.
    # This is only the initial viewing angle.
    ax.view_init(elev=30, azim=-20, roll=90)
    
    # -----------------------------

    img = ax.scatter(coordinates[:,0], coordinates[:,1], coordinates[:,2], c=vals, norm='log') 
    fig.colorbar(img)
    plt.show()



"""
Input:
    direction - String. 
                Used to select the direction of interest. One of 'x', 'y', 'z'.
    vals      - ndarray.
                1D ndarray assumed to be formatted like [point1_valx, point1_valy,
                  point1_valz, point2_valx, ...]
    
Functionality:
    Selects the values for the direction of interest.

Return:
    1D ndarray formatted like [point1_val1, point2_val1, point3_val3, ...]
"""
def select_dir(direction, vals):
    
    assert(direction == 'x' or direction == 'y' or direction == 'z')
    assert(vals.ndim == 1)

    # Map direction to an integer.
    # Guarenteed to be 0, 1, or 2.
    direction_int = int(direction, base=36) - int('x', base=36)

    return vals[direction_int::3]
    


"""
Input:
    coordinates - ndarray.
                  3D coordinates of interest.
                  Assumed to be 2D ndarray with shape (N, 3).
    vals        - ndarray
                  The values at each point in space.
                  Assumed to be 1D array of shape (3*N).
    factor      - Float.
                  The factor to downsample by. A factor of .5 corresponds to
                    using only 50% of points.

Functionality:
    It is assumed that the values are passed in the form [point1_val1, point1_val2,
      point1_val3, point2_val1, ...]         
    Returns a selection of the points and associated values.
    Some rounding is done if the factor does not produce an integer number of
      points.

Return:
   Tuple: (coords, vals).
   The coords entry is a 2D ndarray with shape (M, 3).
   The vals entry is a 1D ndarray with form identical to what was passed in.
"""
def select_and_downsample(coordinates, vals, factor):
    
    assert(coordinates.shape[1] == 3)
    assert(vals.ndim == 1)

    num_points = coordinates.shape[0]
        
    # Randomly select the indices of points to keep.
    # num_points_to_keep = int(num_points * factor) 
    # indices_points_to_keep = random.sample(range(num_points), num_points_to_keep)

    # Don't do random selection, it makes it difficult to see the shape of the
    #   part.
    one_of_every = int(1 / factor)
    indices_points_to_keep = range(0, num_points, one_of_every)
    num_points_to_keep = len(indices_points_to_keep)

    # Now save off the corresponding points and values.
    kept_points = np.zeros((num_points_to_keep, 3))
    kept_vals = np.zeros(3 * num_points_to_keep) 
    for i, idx in enumerate(indices_points_to_keep):
        kept_points[i, :] = coordinates[idx, :] 
        kept_vals[3*i:3*i+3] = vals[3*idx:3*idx+3] 

    return (kept_points, kept_vals)



if __name__ == '__main__':

    nodal_coords_dir = './nodal_coords/' 
    nodal_coords_name = 'shot_peen_bar_nodalcoords.npy'

    disp_vec_dir = './displacement_vectors/' 
    disp_vec0_name = 'shot_peen_bar_dispv0.npy'
    disp_vec1_name = 'shot_peen_bar_dispv1.npy'
   
    stiff_mat_dir = './stiffness_matrices/'
    stiffness_mat_name = 'shot_peen_bar_mtx.npz'

    # Careful, the stiffness matrix is in sparse scipy format while the others 
    #   are in the normal dense numpy format.
    nodal_coords = pso.load_saved_array(nodal_coords_dir + nodal_coords_name)
    disp_vec0 = pso.load_saved_array(disp_vec_dir + disp_vec0_name)
    disp_vec1 = pso.load_saved_array(disp_vec_dir + disp_vec1_name)
    stiffness_matrix = pso.load_saved_array(stiff_mat_dir + stiffness_mat_name)

    # Some sanity checking.
    # These should be true as long as I don't reformat how things are saved.
    assert(disp_vec0.size == disp_vec1.size)
    assert(stiffness_matrix.shape[0] == disp_vec0.size)
    assert(stiffness_matrix.shape[1] == disp_vec0.size)
    assert(3 * nodal_coords.shape[0] == disp_vec0.size)

    # Experimental Setup:
    #   Initial Positions -> nodal_coords
    #   Equilibirum Step  -> disp_vec0
    #   Cutting Step      -> disp_vec1

    # Compute the displacement due to the cut.
    disp_due_to_cut = disp_vec1 - disp_vec0

    # Compute the coordinates of each node before the cut happens.
    pre_cut_nodal_coords = nodal_coords + disp_vec0.reshape((-1, 3))

    # Compute the force vector.
    force_vec = stiffness_matrix.dot(disp_due_to_cut)

    points_to_disp, vals_to_disp = select_and_downsample(pre_cut_nodal_coords, force_vec, .01) 

    plot_directional_forces(points_to_disp, vals_to_disp, 'x')

    


