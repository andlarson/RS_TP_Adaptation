import numpy as np
import parse_sim_output as pso 
import scipy as sp
import plotting



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
    Triple of list of floats. The function parameters which achieved a good fit.      
"""
def displacement_curve_fit(disp_vec, coords, f):

    x_disps = disp_vec[::3]
    y_disps = disp_vec[1::3]
    z_disps = disp_vec[2::3]

    x_params = sp.optimize.curve_fit(f, coords.T, x_disps)[0]
    y_params = sp.optimize.curve_fit(f, coords.T, y_disps)[0]
    z_params = sp.optimize.curve_fit(f, coords.T, z_disps)[0]

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
def disp_func(data, a, b, c, d, e, f, g, h, i, j):

    x = data[0]
    y = data[1]
    z = data[2]

    p1 = a * x**3 + b * x**2 + c * x**1
    p2 = d * y**3 + e * y**2 + f * y**1
    p3 = g * z**3 + h * z**2 + i * z**1

    return p1 + p2 + p3 + j



"""
Input:
    data  - 1D numpy array.
            Assumed to be length 3.
    a,... - Floats.
            The model parameters.
            Some may not be used.
    wrt   - String. One of 'x', 'y', 'z'.
            The direction to differentiate the displacement function with
              respect to.

Functionality:
    Evaluates the derivative of the displacement function with respect to
      a particular direction at a particular point.
    WARNING: THIS FUNCTION IS COUPLED TO disp_func().

Return:
    Scalar. The function value.
"""
def der_disp_func(data, a, b, c, d, e, f, g, h, i, j, wrt):
    
    if wrt == 'x':
        return 3 * a * data[0]**2 + 2 * b * data[0] + c
    elif wrt == 'y':
        return 3 * d * data[1]**2 + 2 * e * data[1] + f
    elif wrt == 'z':
        return 3 * g * data[2]**2 + 2 * h * data[2] + i
    else:
       raise RuntimeError("Can't take a derivative wrt this!") 



"""
Input:
    nodal_coords - 2D ndarray of shape (N, 3).
                   The coordinates to evaluate the strain tensors at.
    x_params     - The parameters for the dx displacement function.
    y_params     - The parameters for the dy displacement function.
    z_params     - The parameters for the dz displacement function.

Functionality:
    Computes the approximated strain tensor at each node in space using the
      approximate displacement functions.

Return:
    2D ndarray of shape (N, 6). 
"""
def compute_approx_strain(nodal_coords, x_params, y_params, z_params):
   
    node_cnt = nodal_coords[:, 0].size
    approx_nodal_strains = np.zeros((node_cnt, 6))

    for node_idx in range(node_cnt):
        node = nodal_coords[node_idx, :]

        eps_x = der_disp_func(node, *x_params, 'x')
        eps_y = der_disp_func(node, *y_params, 'y')
        eps_z = der_disp_func(node, *z_params, 'z')
        gamma_xy = der_disp_func(node, *x_params, 'y') + der_disp_func(node, *y_params, 'x')
        gamma_xz = der_disp_func(node, *x_params, 'z') + der_disp_func(node, *z_params, 'x')
        gamma_yz = der_disp_func(node, *y_params, 'z') + der_disp_func(node, *z_params, 'y')

        approx_nodal_strains[node_idx, :] = np.array([eps_x, eps_y, eps_z, gamma_xy, gamma_xz, gamma_yz])

    return approx_nodal_strains



"""
Input:
    nodal_strains - 2D numpy array of shape (N, 6).
    E             - Modulus of elasticity.
    v             - Poisson's ratio.
    
Functionality:
    Maps the nodal strain tensors to nodal stress tensors. Assumes an elastic
      i.e. Hookian relationship.

Return:
    2D ndarray of shape (N, 6). The nodal stresses which caused the observed
      strains.
"""
def map_strains_to_stresses(nodal_strains, E, v):
    
    # The modulus of elasticity in shear i.e. modulus of rigidity.
    G = E / (2 * (1 + v))

    nodal_stresses = np.zeros(nodal_strains.shape)
    node_cnt = nodal_strains[:, 0].size
    for idx in range(node_cnt):
        strains = nodal_strains[idx, :]
        sigma_x = strains[0] * E
        sigma_y = strains[1] * E
        sigma_z = strains[2] * E
        tau_xy = strains[3] * G
        tau_xz = strains[4] * G
        tau_yz = strains[5] * G
        nodal_stresses[idx, :] = np.array([sigma_x, sigma_y, sigma_z, tau_xy, tau_xz, tau_yz])

    return nodal_stresses



if __name__ == '__main__':

    # Get the displacement data due to the cut.
    disp_vec = pso.load_saved_array('./displacement_vectors/shot_peen_bar_dispv1.npy')

    # Get the nodal coordinates before the displacement occurs.
    nodal_coords = pso.load_saved_array('./nodal_coords/shot_peen_bar_nodalcoords.npy')

    x_params, y_params, z_params = displacement_curve_fit(disp_vec, nodal_coords, disp_func)

    # Visualize each displacement function.
    # The displacement functions can only be expected to be valid over the
    #   region that the part exists.
    # For the shot peened bar, this is generally over the rectangular region
    #   (0, 0, 0) -> (40, 10, 400).
    """
    # Evaluate the approximate displacement function at each node in space.
    approx_disps = []
    for node_idx in range(nodal_coords[:, 0].size):
        approx_disps.append(disp_func(nodal_coords[node_idx, :], *x_params)) 
        approx_disps.append(disp_func(nodal_coords[node_idx, :], *y_params)) 
        approx_disps.append(disp_func(nodal_coords[node_idx, :], *z_params)) 
    approx_disps = np.array(approx_disps)

    to_plot = plotting.select_and_downsample(nodal_coords, approx_disps, .01)
    plotting.plot3d(to_plot[0], to_plot[1][0::3], "Displacements In The x Direction")    
    """ 

    # Visualize the components of stress.
    approx_nodal_strains = compute_approx_strain(nodal_coords, x_params, y_params, z_params)
    E = 200000000000
    v = .3
    approx_nodal_stresses = map_strains_to_stresses(approx_nodal_strains, E, v)
    
    # Do manual down sampling and select the component of stress of interest.
    sel_approx_nodal_stresses = approx_nodal_stresses[::100, 2]
    sel_nodal_coords = nodal_coords[::100, :]
    plotting.plot3d(sel_nodal_coords, sel_approx_nodal_stresses, "Approximated 3,3 Component of Stress")

     

