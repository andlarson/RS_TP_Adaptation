"""
Contains functionality associated with parts.
"""

from abc import ABCMeta, abstractmethod

import src.util.geom as geom
import src.core.material_properties.material_properties as mp 
import src.core.abaqus.abaqus_shim as shim



class MinimalPart():
    
    def __init__(self, path: str) -> None:
        """Creates a part based on one which exists in Abaqus. This part is
               minimally defined - it is just a geometry and has no material,
               stress subroutine, etc. associated with it.
            
           Args:
               path: Absolute path to .cae file defining the part.
        
           Returns:
               None.
        
           Raises:
               None.
        """

        mdb = shim.use_mdb(path)

        if not shim.check_simple_standard_mdb(True, mdb):
            raise RuntimeError("The part, as it exists in Abaqus, is not basic!")

        shim.close_mdb(mdb)

        self.path_to_mdb = path


