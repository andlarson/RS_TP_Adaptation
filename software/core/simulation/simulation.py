
"""
Abstractions for the Abaqus simulations.
"""


# TODO:
# This is really an abstract base class and the proper Python infrastructure
#    should be used. No one should ever create an object of this type.
class Simulation:

    def __init__(self):
    # type: (None) -> None
  
        self.metadata = SimulationMetadata()



class ForwardDirection(Simulation):
    
    def __init__(self):
    # type: (None) -> None

        pass



class BackwardDirection(Simulation):

    def __init__(self):
    # type: (None) -> None

        pass


