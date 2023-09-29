import util.geom as geom
import core.material_properties.material_properties as mat_props 
import core.abaqus.abaqus_shim as shim



class Part:

    def __init__(self, name, material):
    # type: (str, material_properties.Material) -> None

        self.name = name
        self.path_to_stress_subroutine = None
        self.material = material
        

    def update_part_with_real_data(self): 
    # type: (None) -> None
       
        raise RuntimeError("Not yet supported.")



class UserDefinedPart(Part):
    
    def __init__(self, name, part_rep, material_properties):
    # type: (str, geom.SpecRightRectPrism, mat_props.MaterialProperties) -> None

        raise RuntimeError("Can't handle user defined parts right now...")



class AbaqusDefinedPart(Part):
    
    def __init__(self, name, path, material):
    # type: (str, str, material_properties.Material) -> None
    
        Part.__init__(self, name, material)

        mdb = shim.use_mdb(path)

        assert(shim.check_basic_geom(True, mdb))

        shim.close_mdb(mdb)

        self.path_to_mdb = path