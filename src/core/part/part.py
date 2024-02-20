"""
Contains functionality associated with parts.
"""

from abc import ABCMeta, abstractmethod

import src.util.geom as geom
import src.core.material_properties.material_properties as mp 
import src.core.abaqus.abaqus_shim as shim


class InitialPart():
    
    def __init__(self, name: str, path: str, material: mp.ElasticMaterial) -> None:
        """Creates a part based on one which exists in Abaqus.

           Assumes that the part, as it is defined in Abaqus, is defined in a basic
               way.

           Args:
               name:     The name of the part.
               path:     Absolute path to .cae file defining the part. 
               material: The material that the part is made out of.

           Returns:
               None.

           Raises:
               None.
        """

        mdb = shim.use_mdb(path)

        if not shim.check_basic_geom(True, mdb):
            raise RuntimeError("The part, as it exists in Abaqus, is not basic!")

        shim.close_mdb(mdb)

        self.name = name
        self.path_to_stress_subroutine = None
        self.material = material
        self.path_to_mdb = path





class MinimalPart():
    
    def __init__(self, name: str, path: str) -> None:
        """Creates a part based on one which exists in Abaqus. This part is
               minimally defined - it is just a geometry and has no material,
               stress subroutine, etc. associated with it.
            
           Args:
               name: The name of the part.
               path: Absolute path to .cae file defining the part.
        
           Returns:
               None.
        
           Raises:
               None.
        """

        mdb = shim.use_mdb(path)

        if not shim.check_basic_geom(True, mdb):
            raise RuntimeError("The part, as it exists in Abaqus, is not basic!")

        shim.close_mdb(mdb)

        self.name = name
        self.path_to_mdb = path


