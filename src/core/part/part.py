"""
Contains functionality associated with parts.
"""

import src.util.geom as geom
import src.core.material_properties.material_properties as mat_props 
import src.core.abaqus.abaqus_shim as shim


class _Part:

    def __init__(self, name: str, material: mat_props.Material):
    """Creates a top level part.

       Args:
           name:     The name of the part.
           material: The material the part is made out of.

       Returns:
           None.

       Raises:
           None.
    """

        self.name = name
        self.path_to_stress_subroutine = None
        self.material = material
        

    def update_part_with_real_data(self): 
        raise RuntimeError("Not yet supported.")



class UserDefinedPart(_Part):
    
    def __init__(self, name, part_rep, material_properties):
        raise RuntimeError("Can't handle user defined parts right now...")



class AbaqusDefinedPart(_Part):
    
    def __init__(self, name: str, path: str, material: mat_props.Material) -> None:
        """Creates a part based on one which exists in Abaqus.

           Assumes that the part, as it is defined in Abaqus, is defined in a basic
               way.

           Args:
               name:     The name of the part.
               path:     The path of the .cae file defining the part.
               material: The material that the part is made out of.

           Returns:
               None.

           Raises:
               None.
        """
    
        _Part.__init__(self, name, material)

        mdb = shim.use_mdb(path)

        if not shim.check_basic_geom(True, mdb):
            raise RuntimeError("The part, as it exists in Abaqus, is not basic!")

        shim.close_mdb(mdb)

        self.path_to_mdb = path
