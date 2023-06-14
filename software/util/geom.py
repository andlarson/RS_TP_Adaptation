
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
    # Assumed that v5 -> v8 have the same z-ccordinate.
    def __init__(self, v1, v2, v3, v4, v5, v6, v7, v8):
    # type: (Point, Point, Point, Point, Point, Point, Point, Point) -> None

        if (not v1.z == v2.z == v3.z == v4.z) or (not v5.z == v6.z == v7.z == v8.z):
            raise RuntimeError("Improperly constructed rectangular prism!")

        self.vertices = (v1, v2, v3, v4, v5, v6, v7, v8)
        
        


