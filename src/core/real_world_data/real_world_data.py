"""
Provides functionality for dealing with real world data.
"""


class RealWorldDataFromSim:
    
    def __init__(self, path: str):
        """Initializes some real world data produced by a simulation. This sounds
               ludicrous, but the idea is that a simulation with sufficiently high
               mesh density and a hidden, ground truth residual stress field can 
               be considered real world data.
                       
           Args:
               path: Absolute path to raw data resulting from a simulation.
                         Expected to have suffix .odb.
        
           Returns:
               None.
        
           Raises:
               None.
        """
        
        assert(path.endswith(".odb"))
        self.path_to_raw_data: str = path



class RealWorldData:

    def __init__(self, path: str):
        """Initializes some real world data.
                       
           Args:
               path: Absolute path to raw data collected from the real world.
                         Expected to be file with suffix .stl.
        
           Returns:
               None.
        
           Raises:
               None.
        """

        assert(path.endswith(".stl"))
        self.path_to_raw_data: str = path
