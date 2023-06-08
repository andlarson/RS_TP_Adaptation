


class part:

    __init__(self, part_rep, initial_stress_profile, material_properties):
        self.part_rep = part_rep
        self.initial_stress_profile = initial_stress_profile
        self.material_properties = material_properties

    update_part_with_real_data(part_rep): 
        pass



class part_representation:
   
    __init__(self, rep):
       
        # May get a vertex representation of the part or an Abaqus based
        #   representation of the part.
        if isinstance(rep, type(vertex_rep)):
            pass
        elif isinstance(rep, type(abaqus_rep)):
            pass
        else:
            raise AssertionError("The part representation must be an instance \
                                  of vertex_rep or abaqus_rep!")



class vertex_rep:

    __init__(self, rep):
        
        # TODO: Sanity checking of the passed vertex representation.

        self.rep = rep



class abaqus_rep:

    __init__(self, rep):
        
        # TODO: Sanity checking of the passed Abaqus representation. We were
        #   passed an Abaqus-compatible file right?!?!?

        self.rep = rep 




