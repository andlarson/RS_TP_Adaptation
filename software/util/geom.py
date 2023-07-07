
"""
Utilities for general geometric tasks.
"""

import sys

class Point2DXY:
    
    def __init__(self, x, y):
    # type: (float, float) -> None
        
        self.x = x
        self.y = y


class Point3D:

    def __init__(self, x, y, z):
    # type: (float, float, float) -> None

        self.x = x
        self.y = y
        self.z = z

    def proj_xy(self):
    # type: (None) -> tuple[float, float]
        
        return Point2DXY(self.x, self.y) 



# A right rectangular prism is a 3D object which consists of 8 vertices,
#   all right angles, and opposite faces have equal area.
# This is not only a right rectangular prism, but also a right rectangular
#   prism which edges which are parallel to the standard x, y, and z axes.
# TODO: Generalize so that edges need not be parallel to the standard x, y, and
#    z axes.
class SpecRightRectPrism:
   
    def __init__(self, v1, v2, v3, v4, v5, v6, v7, v8):
    # type: (Point3D, Point3D, Point3D, Point3D, Point3D, Point3D, Point3D, Point3D) -> None

        v_list = [v1, v2, v3, v4, v5, v6, v7, v8]

        # Find the groups of 4 vertices which share the same x coordinates, y
        #   coordinates, and z coordinates.
        self.same_x_g1 = [v1]
        self.same_x_g2 = []
        self.same_y_g1 = [v1]
        self.same_y_g2 = []
        self.same_z_g1 = [v1]
        self.same_z_g2 = []
        for v in v_list[1:]:
            if v.x == self.same_x_g1[0].x:
                self.same_x_g1.append(v)
            else:
                self.same_x_g2.append(v)

            if v.y == self.same_y_g1[0].y:
                self.same_y_g1.append(v)
            else:
                self.same_y_g2.append(v)

            if v.z == self.same_z_g1[0].z:
                self.same_z_g1.append(v)
            else:
                self.same_z_g2.append(v)
        
        assert(len(self.same_x_g1) == 4)
        assert(len(self.same_x_g2) == 4)
        assert(len(self.same_y_g1) == 4)
        assert(len(self.same_y_g2) == 4)
        assert(len(self.same_z_g1) == 4)
        assert(len(self.same_z_g2) == 4)

        self.vertices = v_list

   
    def get_smaller_z(self):
        return min(self.same_z_g1[0].z, self.same_z_g2[0].z) 


    # Get the length in the x-direction, the width in the y-direction, and the
    #   height in the z-direction.
    def get_dims(self):
    # type: (None) -> tuple(float, float, float)
       
        x_length = abs(self.same_x_g1[0].x - self.same_x_g2[0].x)
        y_width = abs(self.same_y_g1[0].y - self.same_y_g2[0].y)
        z_height = abs(self.same_z_g1[0].z - self.same_z_g2[0].z)

        return (x_length, y_width, z_height)

    
    # Get 2 vertices which have the same z coordinates but differing x and y
    #   coordinates.
    def get_rect_corners(self):
    # type: (None) -> tuple[Point3D, Point3D] 

        v1 = self.vertices[0]

        for v in self.vertices[1:]:
            if (v1.z == v.z) and (v1.x != v.x) and (v1.y != v.y):
                return v1, v


    def get_centroid(self):
    # type: (None) -> Point3D
        
        avg_x = (self.same_x_g1[0].x + self.same_x_g2[0].x)/2
        avg_y = (self.same_y_g1[0].y + self.same_y_g2[0].y)/2
        avg_z = (self.same_z_g1[0].z + self.same_z_g2[0].z)/2

        return Point3D(avg_x, avg_y, avg_z) 


