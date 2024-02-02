"""
Contains functionality associated with parts.
"""

from abc import ABCMeta, abstractmethod

import src.util.geom as geom
import src.core.material_properties.material_properties as mp 
import src.core.abaqus.abaqus_shim as shim


class Part(metaclass=ABCMeta):
    
    @abstractmethod
    def update_part_with_real_data(self) -> None:
        """Updates the representation of a part with data from the real world."""



class UserDefinedPart(Part):
    
    def __init__(self, name, part_rep, material_properties):
        raise RuntimeError("Can't handle user defined parts right now...")

    def update_part_with_real_data(self) -> None:
        """Updates the representation of a part with data from the real world."""
        raise RuntimeError("Not yet implemented.")



class AbaqusDefinedPart(Part):
    
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

        self.name = name
        self.path_to_stress_subroutine = None
        self.material = material

        mdb = shim.use_mdb(path)

        if not shim.check_basic_geom(True, mdb):
            raise RuntimeError("The part, as it exists in Abaqus, is not basic!")

        shim.close_mdb(mdb)

        self.path_to_mdb = path

    def update_part_with_real_data(self) -> None:
        """Updates the representation of a part with data from the real world."""
        raise RuntimeError("Not yet implemented.")
