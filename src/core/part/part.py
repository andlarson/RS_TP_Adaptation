import util.geom as geom
import core.material_properties.material_properties as mat_props 
import core.abaqus.abaqus_shim as shim



class Part:

    def __init__(self, name):
    # type: (str) -> None

        self.name = name
        self.path_to_stress_subroutine = None
        

    def update_part_with_real_data(self): 
    # type: (None) -> None
       
        raise RuntimeError("Not yet supported.")


    # Associate a stress profile, as defined by a particular subroutine file,
    #    with this part.
    def add_stress_profile(self, path_to_subroutine):
    # type: (str) -> None
         
         self.path_to_stress_subroutine = path_to_subroutine



class UserDefinedPart(Part):
    
    # TODO:
    def __init__(self, name, part_rep, material_properties):
    # type: (str, geom.SpecRightRectPrism, mat_props.MaterialProperties) -> None

        raise RuntimeError("Can't handle user defined parts right now...")



class AbaqusDefinedPart(Part):
    
    def __init__(self, name, path):
    # type: (str, str) -> None
    
        Part.__init__(self, name)

        mdb = shim.use_mdb(path)

        assert(shim.check_init_geom(True, mdb))

        shim.close_mdb(mdb)

        self.path_to_mdb = path