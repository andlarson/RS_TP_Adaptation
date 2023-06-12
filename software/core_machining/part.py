# An implication of using relative importing is that this file cannot be used
#   as a script. Thus any tests must be located elsewhere.
from ..util import geom



class VertexRep:

    def __init__(self, rep):
      
        # TODO: For now, we assume the vertices properly define a part in
        #   Abaqus.
        self.rep = rep



class AbaqusRep:

    def __init__(self, rep):
        
        # TODO: Sanity checking of the passed Abaqus representation. We were
        #   passed an Abaqus-compatible file right?!?!?

        self.rep = rep 



class PartRepresentation:
   
    def __init__(self, rep):
       
        if type(rep) == VertexRep:
            pass 
        elif type(rep) == AbaqusRep:
            pass



class Part:

    def __init__(self, part_rep, initial_stress_profile, material_properties):
        self.part_rep = part_rep
        self.initial_stress_profile = initial_stress_profile
        self.material_properties = material_properties

    def update_part_with_real_data(part_rep): 
        pass








