

class Material:

    def __init__(self):
    # type: (None) -> None

        raise RuntimeError("Not yet supported!")        



class ElasticMaterial(Material):

    # Create an isotropic elastic material.
    def __init__(self, poissons_ratio, youngs_modulus):
    # type: (float, float) -> None

        self.poissons_ratio = poissons_ratio
        self.youngs_modulus = youngs_modulus