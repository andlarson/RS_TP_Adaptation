import numpy as np
import stiff_disp_mult as sdm


if __name__ == '__main__':

    nodal_coords_dir = './nodal_coords/' 
    force_vector_dir = './force_vectors/'
    nodal_coords_name = 'shot_peen_bar_nodal_coords.npy'
    force_vector_name = 'shot_peen_bar_force_vec.npy'

    nodal_coords = sdm.load_saved_array(nodal_coords_dir + nodal_coords_name)
    force_vector = sdm.load_saved_array(force_vector_dir + force_vector_name)

    # The number of nodes better match.



