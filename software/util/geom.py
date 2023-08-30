import sys
import math

import numpy as np
import numpy.linalg as linalg
import matplotlib.path as path

# DEBUG
from debug import *



class Point2D:
    
    def __init__(self, x1, x2):
    # type: (Union[float, int], Union[float, int]) -> None
        
        self.x1 = float(x1)
        self.x2 = float(x2)


    def compute_scale_factors(self):
    # type: (None) -> Tuple[float]

        try:
            t1 = self.x1 / self.x2
        except ZeroDivisionError:
            if float_equals(self.x1, 0):
                t1 = float("nan")
            elif self.x1 > 0:
                t1 = float("+inf")
            else:
                t1 = float("-inf")

        return (t1, )


    def get_components(self):
    # type: (None) -> Tuple[float, float]

        return (self.x1, self.x2)



class Point3D:

    def __init__(self, x, y, z):
    # type: (Union[float, int], Union[float, int], Union[float, int]) -> None

        self.x = float(x)
        self.y = float(y)
        self.z = float(z)


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
    # type: (None) -> Tuple[float, float]

        try:
            t1 = self.x / self.z
        except ZeroDivisionError:
            if float_equals(self.x, 0):
                t1 = float("nan")
            elif self.x > 0:
                t1 = float("+inf")
            else:
                t1 = float("-inf")

        try:
            t2 = self.y / self.z
        except ZeroDivisionError:
            if float_equals(self.y, 0):
                t2 = float("nan")
            elif self.y > 0:
                t2 = float("+inf")
            else:
                t2 = float("-inf")

        return (t1, t2)



class Vec2D:

    def __init__(self, x, y):
    # type: (Union[float, int], Union[float, int]) -> None

        if float_equals(x, 0) and float_equals(y, 0):
            raise RuntimeError("Attempt to create vector which is meaningless.")

        self.np_arr = np.array([float(x), float(y)])


    # Check if the vector has unit norm.
    def is_unit(self):
    # type: (None) -> bool
        
        if float_equals(self.get_len(), 1):
            return True
        return False


    # Get a normalized version of the vector, without modifying internal
    #    state.
    def normalize(self):
    # type: (None) -> Vec3D

        res = (1 / linalg.norm(self.np_arr)) * self.np_arr
        return Vec2D(*res)


    def get_len(self):
    # type: (None) -> float

        return linalg.norm(self.np_arr)



class Vec3D:

    def __init__(self, x, y, z):
    # type: (Union[float, int], Union[float, int], Union[float, int]) -> None

        if float_equals(x, 0) and float_equals(y, 0) and float_equals(z, 0):
            raise RuntimeError("Attempt to create vector which is meaningless.")

        self.np_arr = np.array([float(x), float(y), float(z)])


    # Check if the vector has unit norm.
    def is_unit(self):
    # type: (None) -> bool
        
        if float_equals(self.get_len(), 1):
            return True
        return False


    # Get a normalized version of the vector, without modifying internal
    #    state.
    def normalize(self):
    # type: (None) -> Vec3D

        res = (1 / linalg.norm(self.np_arr)) * self.np_arr
        return Vec3D(*res)


    def get_len(self):
    # type: (None) -> float

        return linalg.norm(self.np_arr)
 
    
    # Get a vector orthogonal to the vector. Note that there are an infinite
    #    number of vectors orthogonal to a vector. This function returns
    #    one chosen arbitrarily with unit length.
    def get_orthonormal(self):
    # type: (None) -> Vec3D

        if self.np_arr[1] == 0 and self.np_arr[2] == 0:
            other_vec = np.array([0, 1, 0])
        else:
            other_vec = np.array([1, 0, 0])

        orth_vec = Vec3D(*np.cross(other_vec, self.np_arr))

        return orth_vec.get_norm()
        


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
    # type: (None) -> float

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
    # type: (Union[List[Point3D], Tuple[Point3D, ...]]) -> None

        if not on_any_plane(points):
            raise RuntimeError("The points don't describe a valid n-gon!")

        # TODO:
        # We should really check that an ngon is well defined.
        # In particular, that none of the points are on the interior of the ngon
        #    and that there are no redundant points.
        self.vertices = points


    def get_plane_coeffs(self):
    # type: (None) -> np.ndarray 

        return get_plane_coeffs(self.vertices[0], self.vertices[1], self.vertices[2])

   
    # Projects the ngon onto the xy-plane, returning an NGon2D.
    # If the ngon existed on a plane orthogonal to the xy-plane, then all the
    #    points in the NGon2D will be collinear.
    def proj_xy(self):
    # type: (None) -> NGon2D

        proj_points = [vertex.proj_xy() for vertex in self.vertices]
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



class Line3D:

    def __init__(self, p1, p2):
    # type: (Point3D, Point3D) -> None

        self.p1 = p1

        self.line_dir = Vec3D(p1.x - p2.x, p1.y - p2.y, p1.z - p2.z)


    def is_on(self, point):
    # type: (Point3D) -> bool

        x_parameter = robust_float_div(point.x - self.p1.x, self.line_dir.np_arr[0])
        y_parameter = robust_float_div(point.y - self.p1.y, self.line_dir.np_arr[1])
        z_parameter = robust_float_div(point.z - self.p1.z, self.line_dir.np_arr[2])

        return (robust_float_equals(x_parameter, y_parameter) and 
                robust_float_equals(y_parameter, z_parameter))




class Line2D:

    def __init__(self, p1, p2):
    # type: (Point2D, Point2D) -> None

        self.p1 = p1 

        self.line_dir = Vec2D(p1.x - p2.x, p1.y - p2.y)


    def is_on(self, point):
    # type: (Point2D) -> bool

        x_parameter = robust_float_div(point.x1 - self.p1.x1, self.line_dir.np_arr[0])
        y_parameter = robust_float_div(point.x2 - self.p1.x2, self.line_dir.np_arr[1])

        return robust_float_equals(x_parameter, y_parameter)



def float_equals(a, b):
# type: (float, float) -> bool

    if abs(a - b) <= 10**(-6):
        return True

    return False



def robust_float_div(a, b):
# type: (float, float) -> float

    try:
        res = a / b 
    except ZeroDivisionError:
        if a == 0:
            return float("nan")
        elif a > 0:
            return float("+inf")
        else:
            return float("-inf")

    return res



def robust_float_equals(a, b):
# type: (float, float) -> bool

    if math.isnan(a) and math.isnan(b):
        return True
    elif math.isnan(a) != math.isnan(b):
        return False
    # Bad idea in general for floats. The purpose of this comparison is to check
    #    for matching infinities.
    elif a == b:          
        return True  
    elif float_equals(a, b):
        return True

    return False 



# Check if points are collinear.
def are_collinear(points):
# type: (Tuple[Union[Point3D, Point2D], ...]) -> bool

    assert(len(points) > 2)

    if isinstance(points[0], Point3D):
        
        line = Line3D(*points[:2])

    else:
        
        line = Line2D(*points[:2])

    for point in points:
        if not line.is_on(point):
            return False

    return True



# Given some points, find the coefficients a, b, c, and d so that all the points
#    lie on the plane described by ax + by + cz + d = 0. If the points are
#    collinear, an exception is thrown.
def get_plane_coeffs(point1, point2, point3):
# type: (Point3D, Point3D, Point3D) -> np.ndarray

    if are_collinear([point1, point2, point3]):
        raise RuntimeError("Points are collinear!")

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
        coeffs = np.concatenate((np.linalg.solve(coords, vals), np.array([1])))
    except linalg.LinAlgError:
        vals = np.array([0, 0, 0])
        coeffs = np.concatenate((np.linalg.solve(coords, vals), np.array([0])))

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

    ngon_plane_coeffs = ngon.get_plane_coeffs()

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

        # DEBUG
        dp("collinear")

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

    circular_vertices = ngon.vertices + ngon.vertices[0:1]

    # If the point lies on any edge, it is in the ngon.
    for idx in range(len(circular_vertices)):
         
        y2 = circular_vertices[idx + 1].x2
        y1 = circular_vertices[idx].x2
        x2 = circular_vertices[idx + 1].x1
        x1 = circular_vertices[idx].x1

        # If the line connecting the vertices is vertical, do a manual check.
        # TODO: Some of the floating point comparisons are a bit loose here. 
        if float_equals(x1 - x2, 0):
            if float_equals(point.x1, x1) and (min(y1, y2) <= point.x2 <= max(y1, y2)):
                return True
        # Otherwise, construct the equation of the line and check if the point
        #    lies on it.
        else:
            slope = (y2 - y1) / (x2 - x1)
            intercept = (y2 * x1 - y1 * x2) / (x1 - x2)
            if float_equals(slope * point.x1 + intercept, point.x2):
                return True
        
    # Otherwise, check if it's in the interior. 
    ngon = path.Path(ngon.get_builtin_rep, closed=True)
    return ngon.contains_point(point.get_components(), ngon)



class CSys3D:

    # Describe a coordinate system relative to some global coordinate system.
    # 
    # Notes:
    #    It is expected that all coordinate systems are described relative to
    #       a single global coodinate system.
    # 
    # Arguments:
    #    translation - 1D numpy array with 3 floats. Describes the translation
    #                     in space of coordinate system.
    #    basis       - 2D, 3x3 numpy array of floats. Describes the rotation of
    #                     this coordinate system relative to a global system.
    #                     All basis vectors must be unit length and be orthogonal
    #                     to one another.
    #
    # Returns:
    #    None.
    def __init__(self, translation, basis):
    # type: (Any, Any) -> None

        # Check the validity of the basis matrix.
        v1 = Vec3D(*basis[0])
        v2 = Vec3D(*basis[1])
        v3 = Vec3D(*basis[2])
        assert(v1.is_unit())
        assert(v2.is_unit())
        assert(v3.is_unit())
        assert(are_orthogonal(v1, v2))
        assert(are_orthogonal(v2, v3))
        assert(are_orthogonal(v1, v3))

        self.translation = translation
        self.basis = basis 

    
    # Map points into the coordinate system described by this object.
    #
    # Notes:
    #    This function does a generic mapping. It has no knowledge about the
    #       coordinate system that the points already live in.
    #
    # Arguments:
    #    points - List of Point3D objects.
    # 
    # Returns:
    #    List of Point3D objects.
    def from_global_to_new(self, points):
    # type: (List[Point3D]) -> List[Point3D]

        points_new_csys = []
        for point in points:
            point = np.array([point.x, point.y, point.z])
            point_new_csys = np.matmul(self.basis, (point - self.translation))
            point_new_csys = Point3D(*point_new_csys)
            points_new_csys.append(point_new_csys)

        return points_new_csys



# Check if two vectors are orthogonal.
def are_orthogonal(vec1, vec2):
# type: (Vec3D, Vec3D) -> bool

    if float_equals(np.dot(vec1.np_arr, vec2.np_arr), 0):
        return True
    return False



# Finds a unit norm vector orthogonal to two vectors.
# The two vectors must already be orthogonal or no other vector exists. 
def find_third_orthonormal(vec1, vec2):
# type: (Vec3D, Vec3D) -> Vec3D

    sol_vec = Vec3D(*np.cross(vec1.np_arr, vec2.np_arr))

    return sol_vec.get_norm()

    


