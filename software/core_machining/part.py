

class VertexRep:

    def __init__(self, rep: list[tuple[float, float, float]]) -> None:
        
        # TODO: Sanity checking of the passed vertex representation.

        self.rep = rep



class AbaqusRep:

    # TODO: Type hinting the representation.
    def __init__(self, rep) -> None:
        
        # TODO: Sanity checking of the passed Abaqus representation. We were
        #   passed an Abaqus-compatible file right?!?!?

        self.rep = rep 



class PartRepresentation:
   
    def __init__(self, rep: VertexRep | AbaqusRep) -> None:
       
        if type(rep) == VertexRep:
            pass
        elif type(rep) == AbaqusRep:
            pass



class Part:

    def __init__(self, 
                 part_rep: PartRepresentation, 
                 initial_stress_profile: StressProfile,
                 material_properties: MaterialProperties 
                 ) -> None:
        self.part_rep = part_rep
        self.initial_stress_profile = initial_stress_profile
        self.material_properties = material_properties

    def update_part_with_real_data(part_rep: PartRepresentation) -> None: 
        pass








