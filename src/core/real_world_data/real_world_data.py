"""
Provides functionality for dealing with real world data.
"""



class SimRealWorldData():
    
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



class ProcessedRealWorldData():

    def __init__(self, path: str):
        """Initializes some real world data which has already been converted
               to .stl format.
                       
           Args:
               path: Absolute path to .stl file containing geometry from the
                         real world.
        
           Returns:
               None.
        
           Raises:
               None.
        """

        assert(path.endswith(".stl"))
        self.path_to_raw_data: str = path


