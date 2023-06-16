
"""
Utilities for general geometric tasks.
"""


class Point:

    def __init__(self, x, y, z):
    # type: (float, float, float) -> None

        self.x = x
        self.y = y
        self.z = z



class RectPrism:
    # TODO: This is a bit hacky. We should probably build up from more
    #   fundumental geometric objects such as parallel lines, etc.
   
    # Assumed that v1 -> v4 have the same z-coordiante.
    # Assumed that v5 -> v8 have the same z-coordinate.
    # Assumed that v1, v2, v5, v6 have the same x-coordinate.
    # Assumed that v3, v4, v7, v8 have the same x-coordinate.
    def __init__(self, v1, v2, v3, v4, v5, v6, v7, v8):
    # type: (Point, Point, Point, Point, Point, Point, Point, Point) -> None

        if (not (v1.z == v2.z == v3.z == v4.z)) or \
           (not (v5.z == v6.z == v7.z == v8.z)):
            raise RuntimeError("Improperly constructed rectangular prism!")

        if (not (v1.x == v2.x == v5.x == v6.x)) or \
           (not (v3.x == v4.x == v7.x == v8.x)):

        self.vertices = (v1, v2, v3, v4, v5, v6, v7, v8)
        self.x_len = abs()

    
    def get_dims(self):
        
        
        


