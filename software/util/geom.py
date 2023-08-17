"""
Utilities for general geometric tasks.
"""

import numpy as np
import sys

# DEBUG
from debug import *


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


    def get_xyz(self):
    # type: (None) -> tuple[float, float, float]

        return (self.x, self.y, self.z)



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



class NGon3D:

    def __init__(self, points):
    # type: (List[Point3D]) -> None

        if not on_single_plane(points):
            raise RuntimeError("The points don't describe a valid n-gon!")

        # TODO:
        # Technically we can do further checks here. Not only should all the points
        #    lie on a single plane, but they should not be on the interior of the
        #    n-gon.
        self.vertices = points


    def get_plane_coeffs(self):
    # type: (None) -> np.ndarray 

        return get_plane_coeffs(self.vertices[0], self.vertices[1]. self.vertices[2])



def float_equals(a, b):
# type: (float, float) -> bool

    if abs(a - b) <= 10**(-6):
        return True

    return False



# Check if three points are collinear.
def are_collinear(point1, point2, point3):
# type: (Point3D, Point3D, Point3D) -> bool

    vec1x, vec1y, vec1z = point2.x - point1.x, point2.y - point1.x, point2.z - point1.z
    vec2x, vec2y, vec2z = point3.x - point1.x, point3.y - point1.x, point3.z - point1.z

    potential_scale_factor = vec2x / vec1x

    if float_equals(vec1y * potential_scale_factor, vec2y) and \
       float_equals(vec1z * potential_scale_factor, vec2z):
        return True

    return False



# Given some points, find the coefficients a, b, c, and d so that all the points
#    lie on the plane described by ax + by + cz + d = 0. If the points are
#    collinear, an exception is thrown.
def get_plane_coeffs(point1, point2, point3):
# type: (Point3D, Point3D, Point3D) -> np.ndarray

    if are_collinear(point1, point2, point3):
        raise RuntimeError("Trying to find plane coefficients for some points
        which are collinear!")

    # We cannot know if the plane will pass through the origin (making d = 0)
    #    or not apriori. 

    # We first assume that the plane does not pass through the origin and try
    #    to solve the system of equations accordingly. If there is no solution,
    #    we know that the plane must go through the origin. Thus, the system
    #    of equations is modified accordingly and resolved.
    x1, y1, z1 = point1.x, point1.y, point1.z
    x2, y2, z2 = point2.x, point2.y, point2.z
    x3, y3, z3 = point3.x, point3.y, point3.z

    coords = np.array([[x1, y1, z1], [x2, y2, z2], [x3, y3, z3]])
    vals = np.array([1, 1, 1])
    try:
        coeffs = np.concatenate((np.linalg.solve(coords, vals), np.array([1]))
    except np.linalg.LinAlgError:
        vals = np.array([0, 0, 0])
        coeffs = np.concatenate((np.linalg.solve(coords, vals), np.array([0]))

    return coeffs



# Check if the points all lie on a single plane.
# Assumed that the first three points are not collinear.
def on_any_plane(points):
# type: (List[Point3D]) -> bool

    # Find the coefficients of the plane which passes through the points.
    coeffs = get_plane_coeffs(points[0], points[1], points[2])

    # Check if all the points lie on this plane.
    return on_particular_plane(points, coeffs)



# Check if all points lie on a plane described by the coefficients of the plane.
# The coefficent order should be a, b, c, d where the equation of the plane is
#    ax + by + cz = d.
def on_particular_plane(points, coeffs):
# type: (List[Point3D], Seq[Float]) -> bool

    for point in points:
        if not float_equals(coeffs[0] * point.x + coeffs[1] * point.y + coeffs[2] * point.z, coeffs[3]):
            return False

    return True



# Check if a set of points all lie on the same plane as an ngon. 
def points_on_plane_of_ngon(points, ngon):
# type: (List[Point3D], NGon3D) -> bool

    ngon_plane_coeffs = ngon.plane_coeffs()

    return on_particular_plane(points, ngon_plane_coeffs)
    





