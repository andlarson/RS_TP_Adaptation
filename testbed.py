"""
This file acts as a testbed for individual components of the library.
"""

import numpy as np


if __name__ == "__main__":
   
    t1 = np.array((1, 0, 0))
    t2 = np.array((0, 1, 0))
    t3 = np.array((0, 0, 1))

    d1 = np.array((5, 0, 0))
    d2 = np.array((10, 3.4, 0))
    d3 = np.array((0, 7.7, 2.3))

    # Check that the tractions are linearly independent.
    tractions = np.stack((t1, t2, t3), axis=-1)
    try:
        np.linalg.inv(tractions)
    except BaseException as e:
        raise RuntimeError("Tractions are linearly dependent!") 

    # Solve for the matrix.
    displacements = np.stack((d1, d2, d3), axis=-1)
    res = np.linalg.solve(tractions.T, displacements.T).T

    print(res)
    print(res @ t1)
    
