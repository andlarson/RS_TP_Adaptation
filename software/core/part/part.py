# Importing like this works ONLY when the util directory is in the path that 
#   Python searches when looking for modules. 
# In the main.py, the software/ directory is added to the search path. This
#   then guarantees that the util/ directory is located in the search path.
# If you try to run this file as a script via "python3 part.py", this importing
#   technique will not work.
import util.geom as geom
import core.simulation.simulation as sim
import core.stress.stress as stress
import core.material_properties.material_properties as mat_props 


class PartHistory:
    
    def __init__(self):
    # type: (None) -> None

        # There are really two parallel part geometry histories.
        # One lives in Abaqus-land, the other lives in user-defined land.
        # These representations are co-mingled. When accessed, type checking
        #    must be done.
        self.part_history = []


    def append(self, part):
    # type: (None, Part) -> None
        
        self.part_history.append(part)
        

# TODO:
# This is really an abstract base class and the proper Python infrastructure
#    should be used. No one should ever create an object of this type.
class Part:

    def __init__(self, name):
    # type: (str) -> None

        self.name = name
        

    # TODO:
    def update_part_with_real_data(self): 
    # type: (None) -> None

        pass



class UserDefinedPart(Part):
    
    # TODO:
    def __init__(self, name, part_rep, stress_profile, material_properties):
    # type: (str, geom.SpecRightRectPrism, stress.StressProfile, mat_props.MaterialProperties) -> None

        Part.__init__(self, name)
        pass



# An part built in Abaqus is really both a part and a simulation. This
#    motivates the use of multiple inheritance.
class AbaqusDefinedPart(Part, sim.Simulation):
    
    def __init__(self, name, path):
    # type: (str, str) -> None
    
        Part.__init__(self, name)

        mdb = shim.use_mdb(path)
        shim.verify_init_geom_mdb(mdb)
        shim.close(mdb)

        self.path_to_mdb = path


