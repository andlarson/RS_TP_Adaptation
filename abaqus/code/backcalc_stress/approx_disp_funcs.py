import numpy as np
import parse_sim_output as pso 
import scipy as sp



"""
Input:
    disp_vec - 1D numpy array.
               Assumed to be of the form [node1_U1, node1_U2, node1_U3, node2_U1,
                 ...].
    coords   - 2D numpy array.
               Assumed to be of shape (N, 3). 
    f        - Function. 
               The function used as the model for curve fitting. 
               It must be of the form: f(x, a, b, c, ...). The x can be a scalar or
                 an 2D numpy array, and it encapsulates the function arguments.
                 The a, b, c, .... are the model parameters which are optimized
                 over. They must be scalars.

Functionality:
    Constructs three functions which approximate the displacmenets in the x, y, and z
      directions.

Return:
    Triple of floats. The function parameters which achieved a good fit.      
"""
def displacement_curve_fit(disp_vec, coords, f):

    x_disps = disp_vec[::3]
    y_disps = disp_vec[1::3]
    z_disps = disp_vec[2::3]

    x_params = sp.optimize.curve_fit(f, coords.T, x_disps)
    y_params = sp.optimize.curve_fit(f, coords.T, y_disps)
    z_params = sp.optimize.curve_fit(f, coords.T, z_disps)

    return x_params, y_params, z_params



"""
Input:
    data  - 1D numpy array.
            Assumed to be of length 3.
    a,... - Scalars.
            The model parameters.

Functionality:
    Evaluates the function.

Return:
    Scalar. The function output.
"""
def disp_func(data, a, b, c, d, e, f, g, h, i, j)

    x = data[0]
    y = data[1]
    z = data[2]

    p1 = a * x**3 + b * x**2 + c * x**1
    p2 = d * y**3 + e * y**2 + f * y**1
    p3 = g * z**3 + h * z**2 + i * z**1

    return p1 + p2 + p3 + j



if __name__ == '__main__':

    # Get the displacement data for the simulation.
    disp_vec = pso.load_saved_array('./displacement_vecs/shot_peen_bar_dispv1.npy')

    # Get the nodal coordinates before the displacment occurs.
    nodal_coords = pso.load_saved_array('./nodal_coords/shot_peen_bar_nodalcoords.npy')

    x_params, y_params, z_params = displacement_curve_fit(disp_vec, nodal_coords, disp_func)

    # Evaluate the displacement functions for each node in the mesh.
    x_disps = np.array([disp_func(coords, *x_params) for r in range(nodal_coords[])])

    # Visualize each displacement function.
    # The displacement functions can only be expected to be valid over the
    #   region that the part exists.

    # For the shot peened bar, this is generally over the rectangular region
    #   (0, 0, 0) -> (40, 10, 400).
      


