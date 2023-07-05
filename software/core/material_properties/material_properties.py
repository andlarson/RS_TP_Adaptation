

class ElasticProperties:
    
    def __init__(self, poissons_ratio, youngs_modulus):
    # type: (float, float) -> None
       
       self.poissons_ratio = poissons_ratio
       self.youngs_modulus = youngs_modulus



class MaterialProperties:
    
    def __init__(self, elastic_properties):
    # type: (ElasticProperties) -> None
        
        self.elastic_properties = elastic_properties

