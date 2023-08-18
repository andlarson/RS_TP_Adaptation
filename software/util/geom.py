"""
Utilities for general geometric tasks.
"""

import numpy as np
import matplotlib.Path as path
import sys

# DEBUG
from debug import *


UNDEF = "0/0"
INF = "a/0"


class Point2D:
    
    def __init__(self, x1, x2):
    # type: (float, float) -> None
        
        self.x1 = x1
        self.x2 = x2 


    def compute_scale_factors(self):
    # type: (None) -> Tuple[Union[str, float]]

        if self.x1 == 0 and self.x2 == 0:
            t1 = UNDEF
        elif self.x1 != 0 and self.x2 == 0:
            t1 = INF
        else:
            t1 = self.x1 / self.x2

        return (t1, )


    def get_components(self):
    # type: (None) -> Tuple[float, float]

        return (self.x1, self.x2)



class Point3D:

    def __init__(self, x, y, z):
    # type: (float, float, float) -> None

        self.x = x
        self.y = y
        self.z = z


    def proj_xy(self):
    # type: (None) -> Tuple[float, float]
        
        return Point2D(self.x, self.y) 


    def proj_xz(self):
    # type: (None) -> Tuple[float]

        return Point2D(self.x, self.z)


    def get_components(self):
    # type: (None) -> Tuple[float, float, float]

        return (self.x, self.y, self.z)


    def compute_scale_factors(self):
    # type: (None) -> Tuple[Union[float, str], Union[float, str]]

        if self.x == 0 and self.z == 0:
            t1 = UNDEF
        elif self.x != 0 and self.z == 0:
            t1 = INF
        else:
            t1 = self.x / self.z

        if self.y == 0 and self.z == 0:
            t2 = UNDEF
        elif self.y != 0 and self.z == 0:
            t2 = INF
        else:
            t2 = self.y / self.z

        return (t1, t2)



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
    # type: (None) -> Tuple(float, float, float)
       
        x_length = abs(self.same_x_g1[0].x - self.same_x_g2[0].x)
        y_width = abs(self.same_y_g1[0].y - self.same_y_g2[0].y)
        z_height = abs(self.same_z_g1[0].z - self.same_z_g2[0].z)

        return (x_length, y_width, z_height)

    
    # Get 2 vertices which have the same z coordinates but differing x and y
    #   coordinates.
    def get_rect_corners(self):
    # type: (None) -> Tuple[Point3D, Point3D] 

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



class NGon2D:

    # The points must be defined in order. In other words, the ngon will be
    #    constructed under the assumption that a line segment connects the
    #    first point to the second point, another line segment connects the
    #    second point to the third point, etc. 
    def __init__(self, points):
    # type: (List[Point2D]) -> None

        # TODO:
        # We should really check that an ngon is well defined.
        # In particular, that none of the points are on the interior of the ngon
        #    and that there are no redundant points.
        self.vertices = points


    # Get a representation which is only composed of Python built-in types. 
    # Useful for passing to other libraries.
    def get_builtin_rep(self):
    # type: (None) -> List[Tuple[float, float]] 

        return [(vertex.x1, vertex.x2) for vertex in self.vertices]



class NGon3D:

    # The points must be defined in order. In other words, the ngon will be
    #    constructed under the assumption that a line segment connects the
    #    first point to the second point, another line segment connects the
    #    second point to the third point, etc. 
    def __init__(self, points):
    # type: (List[Point3D]) -> None

        if not on_single_plane(points):
            raise RuntimeError("The points don't describe a valid n-gon!")

        # TODO:
        # We should really check that an ngon is well defined.
        # In particular, that none of the points are on the interior of the ngon
        #    and that there are no redundant points.
        self.vertices = points


    def get_plane_coeffs(self):
    # type: (None) -> np.ndarray 

        return get_plane_coeffs(self.vertices[0], self.vertices[1]. self.vertices[2])

   
    # Projects the ngon onto the xy-plane, returning an NGon2D.
    # If the ngon existed on a plane orthogonal to the xy-plane, then all the
    #    points in the NGon2D will be collinear.
    def proj_xy(self):
    # type: (None) -> NGon2D

        proj_points = [proj_xy(vertex) for vertex in self.vertices]
        return NGon2D(proj_points)


    # Projects the ngon onto the xz-plane, returning an NGon2D.
    # If the ngon existed on a plane orthogonal to the xz-plane, then all the
    #    points in the NGon2D will be collinear.
    def proj_xz(self):
    # type: (None) -> NGon2D

        proj_points = [proj_xz(vertex) for vertex in self.vertices]
        return NGon2D(proj_points)


    # Get a representation which is only composed of Python built-in types. 
    # Useful for passing to other libraries.
    def get_builtin_rep(self):
    # type: (None) -> List[Tuple[float, float, float]] 

        return [(vertex.x, vertex.y, vertex.z) for vertex in self.vertices]



def float_equals(a, b):
# type: (float, float) -> bool

    if abs(a - b) <= 10**(-6):
        return True

    return False



# Check if points are collinear.
def are_collinear(points):
# type: (List[Union[Point3D, Point2D]]) -> bool

    assert(len(points) > 2)

    potential_scale_factors = points[0].compute_scale_factors() 

    for point in points: 
        scale_factors = point.compute_scale_factors() 
        for idx in range(len(potential_scale_factors)):
            if not float_equals(scale_factors[idx], potential_scale_factors[idx]):
                return False
            
    return True 



# Given some points, find the coefficients a, b, c, and d so that all the points
#    lie on the plane described by ax + by + cz + d = 0. If the points are
#    collinear, an exception is thrown.
def get_plane_coeffs(point1, point2, point3):
# type: (Point3D, Point3D, Point3D) -> np.ndarray

    if are_collinear([point1, point2, point3]):
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



# Check if all the points lie in the ngon.
def points_in_ngon_3D(points, ngon):
# type: (List[Point3D], NGon3D) -> bool

    for point in points:
        if not point_in_ngon_3D(point, ngon):
            return False

    return True



# Check if a point lies inside an ngon. 
# Assumed that both the point and the ngon live in 3D.
def point_in_ngon_3D(point, ngon):
# type: (Point3D, NGon3D) -> bool

    # If the point doesn't lie on the same plane as the ngon, it certainly
    #    isn't inside of it.
    if not points_on_plane_of_ngon([point], ngon):
        return False

    # If it does lie on the same plane as the ngon, we also need to check that
    #    it lies inside of the ngon.

    # Claim: To reduce this problem to 2D, we can project the points which
    #    make up the ngon and the point onto any plane which is not orthogonal
    #    to the plane that the point/ngon lie on. Then we can apply a 2D
    #    algorithm to check if the new point lies within the new ngon.

    # Project the point and the ngon onto the xy-plane.
    proj_point = point.proj_xy()
    proj_ngon = ngon.proj_xy()

    # If all of the points are collinear, then project onto a different plane.
    if are_collinear(proj_ngon.vertices + [proj_point]):
        proj_point = point.proj_xz()
        proj_ngon = ngon.proj_xz()

    return point_in_ngon_2D(proj_point, proj_ngon)



# Check if a point lies inside an ngon.
# Assumed that both the point and the ngon live in 2D.
# This algorithm assumes that the points which compose the ngon are in-order. 
# If the point lies on one of the edges of the ngon (within float approximation), 
#    it is deemed to be in the ngon.
def point_in_ngon_2D(point, ngon):
# type: (Point2D, NGon2D) -> bool    

    circular_vertices = ngon.vertices + ngon.vertices[0]

    # If the point lies on any edge, it is in the ngon.
    for idx in range(len(circular_vertices) - 1):
         
        y2 = circular_vertices[idx + 1].x2
        y1 = circular_vertices[idx].x2
        x2 = circular_vertices[idx + 1].x1
        x1 = circular_vertices[idx].x1
        slope = (y2 - y1) / (x2 - x1)
        intercept = (y2 * x1 - y1 * x2) / (x1 - x2)
        if float_equals(slope * point.x1 + intercept, point.x2):
            return True
        
    # Otherwise, check if it's in the interior. 
    ngon = path.Path(ngon.get_builtin_rep, closed=True)
    return ngon.contains_point(point.get_components(), ngon)
