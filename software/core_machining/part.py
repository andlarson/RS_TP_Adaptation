# Importing like this works ONLY when the util directory is in the path that 
#   Python searches when looking for modules. 
# In the main.py, the software/ directory is added to the search path. This
#   then guarantees that the util/ directory is located in the search path.
# If you try to run this file as a script via "python3 part.py", it will not
#   work!
import util.geom as geom



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








