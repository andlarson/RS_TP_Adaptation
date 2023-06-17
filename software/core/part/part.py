# Importing like this works ONLY when the util directory is in the path that 
#   Python searches when looking for modules. 
# In the main.py, the software/ directory is added to the search path. This
#   then guarantees that the util/ directory is located in the search path.
# If you try to run this file as a script via "python3 part.py", it will not
#   work!
import util.geom as geom
import core.abaqus.catch_all as catch_all 

import os


VERTEX_REP = "VertexRep"
ABAQUS_MDB = "AbaqusMDB"


class VertexRep:

    def __init__(self, rep):
    # type: (list[tuple[int, int, int]]) -> None
      
        # TODO: For now, we assume the vertices properly define a part in
        #   Abaqus.
        self.rep = rep



class AbaqusMDB:

    def __init__(self, rep):  
    # type: (str) -> None

        if not rep.endswith(".cae"):
            raise RuntimeError("An Abaqus representation is a .cae file.")
        
        if not os.access(rep, os.F_OK):
            raise RuntimeError("Bad path to .cae file passed.")

        mdb = catch_all.use_mdb(rep)
        if len(mdb.models) != 1:
            raise RuntimeError("Too many or too few models in the MDB.")
        if not ("Model-1" in mdb.models):
            raise RuntimeError("The single model in the MDB is improperly named.")
        if len(mdb.models["Model-1"].parts) != 1:
            raise RuntimeError("Too many of too few parts in the MDB.")
        if not ("Initial_Geometry" in mdb.models["Model-1"].parts):
            raise RuntimeError("The single part in the single model in the MDB
                                is improperly named.")

        self.rep = rep 



class PartRepresentation:
   
    def __init__(self, rep):
    # type: (VertexRep | AbaqusMDB) -> None
       
        if type(rep) == VertexRep:
            self.rep = rep 
            self.format = VERTEX_REP 
        elif type(rep) == AbaqusMDB:
            self.rep = rep
            self.format = ABAQUS_MDB 
        else:
            raise RuntimeError("Bad type for the construction of a \
                                PartRepresentation object.")



class PartHistory:
    
    def __init__(self):
    # type: (None) -> None

        # There are really two parallel part geometry histories.
        # One lives in Abaqus-land.
        # The other lives in vertex-land.
        # These representations are co-mingled. When used, type checking can
        #   be done.
        # The append operation ensures that the histories are chronological.
        self.part_geom_history = []

        # There is a single representation format for a stress profile.
        self.stress_profile_history = []


    def append_geom_history(self, name, part_rep):
    # type: (str, PartRepresentation) -> None

        self.part_geom_history.append((name, part_rep))


    def append_stress_profile(self, name, stress_profile):
    # type: (str, StressProfile) -> None

        self.stress_profile_history.append((name, stress_profile))



class Part:

    # An initial stress profile is only needed if the initial part representation
    #   does not include one. This also applies to the material properties.
    def __init__(self, initial_part_rep, initial_stress_profile=None,
                 material_properties=None):
    # type: (PartRepresentation, StressProfile, MaterialProperties) -> None

        if initial_part_rep.format == VERTEX_REP and \
           (initial_stress_profile == None or material_properties == None):
            raise RuntimeError("If you don't create a Part via a .cae file, \
                                then it must have an initial stress profile \
                                and material properties.")

        self.initial_part_rep = initial_part_rep
        self.initial_stress_profile = initial_stress_profile
        self.material_properties = material_properties

        self.part_history = PartHistory()
        self.part_history.append_geom_history("initial", initial_part_rep)
        self.part_history.append_stress_profile("initial", initial_stress_profile)


    # TODO:
    def update_part_with_real_data(self, part_rep): 
        pass








