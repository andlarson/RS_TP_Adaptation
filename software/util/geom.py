
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
    # type: (tuple[Point, Point, Point, Point, Point, Point, Point, Point]) -> None

        if (not (v1.z == v2.z == v3.z == v4.z)) or \
           (not (v5.z == v6.z == v7.z == v8.z)):
            raise RuntimeError("Improperly constructed rectangular prism!")

        if (not (v1.x == v2.x == v5.x == v6.x)) or \
           (not (v3.x == v4.x == v7.x == v8.x)):
            raise RuntimeError("Improperly constructed rectangular prism!")

        self.vertices = (v1, v2, v3, v4, v5, v6, v7, v8)
        self.x_len = abs()


    # Compute the length, width, and height of rectangular prism, regardless of
    #   its location in space.
    # The x direction is length, y direction is width, and z direction is
    #   height.
    # The return format is (length, width, height).
    def get_dims(self):
    # type: (None) -> tuple(float, float, float)
       
        v1 = self.vertices[0]
        for v in self.vertices[1:]:
            if v1.x != v.x and v1.z == v.z and v1.y == v.y:
                x_length = abs(v1.x - v.x)
            if v1.y != v.y and v1.z == v.z and v1.x == v.x:
                y_width = abs(v1.y - v.y)
            if v1.z != v.z and v1.x == v.x and v1.y == v.y:
                z_height = abs(v1.z - v.z) 

        return (x_length, y_width, z_height)

    
    # Get 2 vertices which have the same z coordinates but differing x and y
    #   coordinates.
    def get_rect_corners(self):
    # type: (None) -> tuple[Point, Point] 

        v1 = self.vertices[0]
        for v in self.vertices[1:]:
            if v1.z == v.z and v1.x != v.x and v1.y != v1.y:
                return (v1, v) 


    # Get the center point of the rectangular prism.
    def get_center_point(self):
    # type: (None) -> Point

        # Find the average x-coordinate among the points which live on the same
        #   z plane.

        


