"""
This module contains funcitonlaity related to material properties.
"""



class ElasticMaterial():

    def __init__(self, poissons_ratio: float, youngs_modulus: float):
        """Creates an isotropic elastic material.

           Args:
               poissons_ratio: The Poisson Ratio of the material.
               youngs_modulus: The Young's Modulus of the material.

           Returns:
               None.

           Raises:
               None.
        """

        self.poissons_ratio = poissons_ratio
        self.youngs_modulus = youngs_modulus
