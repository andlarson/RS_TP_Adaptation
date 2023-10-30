import itertools

import numpy as np
import numpy.linalg as linalg
import matplotlib.path as path

from debug import *


class Basis3D:

    # By definition, a basis consists of three vectors which must have unit length
    #    and be orthogonal to one another.
    def __init__(self, v1, v2, v3):
    # type: (Vec3D, Vec3D, Vec3D) -> None

        assert(v1.is_unit())
        assert(v2.is_unit())
        assert(v3.is_unit())
        assert(are_orthogonal(v1, v2))
        assert(are_orthogonal(v2, v3))
        assert(are_orthogonal(v1, v3))

        self.rep = np.stack([v1.rep, v2.rep, v3.rep], axis=0) 



class CSys3D:

    # Describe a coordinate system relative to some global coordinate system.
    # 
    # Notes:
    #    It is expected that all coordinate systems are described relative to
    #       a single global coodinate system.
    # 
    # Arguments:
    #    translation - Vec3D. 
    #                  Describes the translation in space of coordinate system.
    #    basis       - Basis. 
    #
    # Returns:
    #    None.
    def __init__(self, translation, basis):
    # type: (Vec3D, Basis3D) -> None

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
    def map_into_csys(self, points):
    # type: (List[Point3D]) -> List[Point3D]

        return [Point3D(np.matmul(self.basis.rep, (point.rep - self.translation.rep))) for point in points]



class Point2D:
    
    def __init__(self, arr):
    # type: (np.ndarray) -> None

        assert(isinstance(arr, np.ndarray))
        assert(arr.size == 2)
        
        arr = np.array(arr, dtype=float)
        self.rep = arr


    def components(self):
    # type: (None) -> Tuple[np.float, np.float]

        return (self.rep[0], self.rep[1])



class Point3D:

    def __init__(self, arr):
    # type: (np.ndarray) -> None

        assert(isinstance(arr, np.ndarray))
        assert(arr.size == 3)
        
        arr = np.array(arr, dtype=float)
        self.rep = arr


    def proj_xy(self):
    # type: (None) -> Point2D

        return Point2D(self.rep[0:2]) 


    def proj_xz(self):
    # type: (None) -> Point2D 

        return Point2D(np.array([self.rep[0], self.rep[2]]))


    def components(self):
    # type: (None) -> Tuple[np.float, np.float, np.float]

        return (self.rep[0], self.rep[1], self.rep[2])



class Vec2D:

    def __init__(self, arr):
    # type: (np.ndarray) -> None

        assert(isinstance(arr, np.ndarray))
        assert(arr.size == 2)
        
        arr = np.array(arr, dtype=float)
        self.rep = arr



    def is_unit(self):
    # type: (None) -> bool

        if float_equals(self.len(), 1):
            return True
        return False


    def normalize(self):
    # type: (None) -> Vec3D

        res = (1 / linalg.norm(self.np_arr)) * self.np_arr
        return Vec2D(res)


    def len(self):
    # type: (None) -> float

        return linalg.norm(self.np_arr)



class Vec3D:

    def __init__(self, arr):
    # type: (np.ndarrar) -> None

        assert(isinstance(arr, np.ndarray))
        assert(arr.size == 3)
        
        arr = np.array(arr, dtype=float)
        self.rep = arr


    # Check if the vector has unit norm.
    def is_unit(self):
    # type: (None) -> bool
        
        if float_equals(self.len(), 1):
            return True
        return False


    # Get a normalized version of the vector, without modifying internal
    #    state.
    def normalize(self):
    # type: (None) -> Vec3D

        res = (1 / linalg.norm(self.rep)) * self.rep
        return Vec3D(res)


    def len(self):
    # type: (None) -> float

        return linalg.norm(self.rep)
 
    
    # Get a vector orthogonal to the vector. Note that there are an infinite
    #    number of vectors orthogonal to a vector. This function returns
    #    one chosen arbitrarily with unit length.
    def get_orthonormal(self):
    # type: (None) -> Vec3D

        if self.rep[1] == 0 and self.arr[2] == 0:
            other_vec = np.array([0, 1, 0])
        else:
            other_vec = np.array([1, 0, 0])

        orth_vec = Vec3D(np.cross(other_vec, self.rep))

        return orth_vec.get_norm()
        


class PlanarCubicC2Spline3D:

    # A planar cubic spline with continuous first and second derivatives in 3D.
    # 
    # Notes:
    #    Abaqus always uses the a cubic spline with continuous first and second
    #       derivatives to contruct wire features. For this reason, this is a very
    #       natural geometry.
    #    For simplicity, we force the spline to be planar for now. In particular,
    #       the points MUST ALL HAVE THE SAME Y COORDINATE!
    #
    # Arguments:
    #    points - List of Point3D objects.
    #             The ordering of the list is the order in which the points will
    #                be connected.
    #
    # Returns:
    #   None. 
    def __init__(self, points):
    # type: (List[Point3D]) -> None

        assert(len(points) > 1)

        y = points[0].rep[1]
        for point in points:
            if not float_equals(y, point.rep[1]):
                raise AssertionError("Points are not planar!")

        self.y = y
        self.v_list = points



# A right rectangular prism is a 3D object which consists of 8 vertices,
#   all right angles, and opposite faces have equal area.
# This is not only a right rectangular prism, but also a right rectangular
#   prism with faces which are parallel to the standard planes of the coordinate
#   system.
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
            if v.rep[0] == self.same_x_g1[0].rep[0]:
                self.same_x_g1.append(v)
            else:
                self.same_x_g2.append(v)

            if v.rep[1] == self.same_y_g1[0].rep[1]:
                self.same_y_g1.append(v)
            else:
                self.same_y_g2.append(v)

            if v.rep[2] == self.same_z_g1[0].rep[2]:
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

        return min(self.same_z_g1[0].rep[2], self.same_z_g2[0].rep[2]) 


    # Get the length in the x-direction, the width in the y-direction, and the
    #   height in the z-direction.
    def get_dims(self):
    # type: (None) -> Tuple(float, float, float)
       
        x_length = abs(self.same_x_g1[0].rep[0] - self.same_x_g2[0].rep[0])
        y_width = abs(self.same_y_g1[0].rep[1] - self.same_y_g2[0].rep[1])
        z_height = abs(self.same_z_g1[0].rep[2] - self.same_z_g2[0].rep[2])

        return (x_length, y_width, z_height)

    
    # Get 2 vertices which have the same z coordinates but differing x and y
    #   coordinates.
    def get_rect_corners(self):
    # type: (None) -> Tuple[Point3D, Point3D] 

        v1 = self.vertices[0]

        for v in self.vertices[1:]:
            if (v1.rep[2] == v.rep[2]) and (v1.rep[0] != v.rep[0]) and (v1.rep[1] != v.rep[1]):
                return v1, v


    def get_centroid(self):
    # type: (None) -> Point3D
        
        avg_x = (self.same_x_g1[0].rep[0] + self.same_x_g2[0].rep[0]) / 2
        avg_y = (self.same_y_g1[0].rep[1] + self.same_y_g2[0].rep[1]) / 2
        avg_z = (self.same_z_g1[0].rep[2] + self.same_z_g2[0].rep[2]) / 2

        return Point3D(np.array([avg_x, avg_y, avg_z])) 



class NGon2D:

    # The points must be defined in order. In other words, the ngon will be
    #    constructed under the assumption that a line segment connects the
    #    first point to the second point, another line segment connects the
    #    second point to the third point, etc. 
    def __init__(self, vertices):
    # type: (List[Point2D]) -> None

        # TODO:
        # We should really check that an ngon is well defined.
        # In particular, that none of the points are on the interior of the ngon
        #    and that there are no redundant points.
        self.vertices = vertices 


    # Get a representation which is only composed of Python built-in types. 
    # Useful for passing to other libraries.
    def get_builtin_rep(self):
    # type: (None) -> List[Tuple[float, float]] 

        return [(vertex.rep[0], vertex.rep[1]) for vertex in self.vertices]



class NGon3D:

    # The points must be defined in order. In other words, the ngon will be
    #    constructed under the assumption that a line segment connects the
    #    first point to the second point, another line segment connects the
    #    second point to the third point, etc. 
    def __init__(self, vertices):
    # type: (List[Point3D]) -> None

        if not on_any_plane(vertices):
            raise RuntimeError("The points don't describe a valid n-gon!")

        # TODO:
        # We should really check that an ngon is well defined.
        # In particular, that none of the points are on the interior of the ngon
        #    and that there are no redundant points.
        self.vertices = vertices


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

        proj_points = [vertex.proj_xz() for vertex in self.vertices]
        return NGon2D(proj_points)


    # Get a representation which is only composed of Python built-in types. 
    # Useful for passing to other libraries.
    def get_builtin_rep(self):
    # type: (None) -> List[Tuple[np.float, np.float, np.float]] 

        return [(vertex.rep[0], vertex.rep[1], vertex.rep[2]) for vertex in self.vertices]



# An infinite line extends arbitrarily in space. Its length is infinite.
class InfiniteLine3D:

    def __init__(self, p1, p2):
    # type: (Point3D, Point3D) -> None

        # Store a parametric representation of the line.
        self.p1 = p1
        self.line_dir = Vec3D(p1.rep - p2.rep)


    def is_on(self, point):
    # type: (Point3D) -> bool

        return point_on_inf_line(point, (self.p1, self.line_dir))



# An infinite line extends arbitrarily in space. Its length is infinite.
class InfiniteLine2D:

    def __init__(self, p1, p2):
    # type: (Point2D, Point2D) -> None

        self.p1 = p1 

        self.line_dir = Vec2D(np.array([p1.rep[0] - p2.rep[0], p1.rep[1] - p2.rep[1]]))


    def is_on(self, point):
    # type: (Point2D) -> bool

        if point_on_inf_line(point, (self.p1, self.line_dir)):
            return True
        return False



# A finite line has finite length.
class FiniteLine2D:

    # The two points define the length of the line.
    def __init__(self, p1, p2):
    # type: (Point2D, Point2D) -> None

        self.p1 = p1 
        self.p2 = p2

        self.line_dir = Vec2D(np.array([p1.rep[0] - p2.rep[0], p1.rep[1] - p2.rep[1]]))


    def is_on(self, point):
    # type: (Point2D) -> bool

        if point_on_inf_line(point, (self.p1, self.line_dir)) and point_in_box(point, (self.p1, self.p2)):
            return True
        return False



# Check if a point is on a line which has infinite length.
# 
# Notes:
#    None.
# 
# Arguments:
#    point          - Point3D or Point2D.
#    parametric_rep - Tuple containing Point3D and Vec3D or tuple containing Point2D
#                        and Vec2D.
#                     The parametric representation of the line.
#
# Returns:
#    Boolean.
def point_on_inf_line(point, parametric_rep):
# type: (Union[Point3D, Point2D], Union[Tuple(Point3D, Vec3D), Tuple(Point2D, Vec2D)]) -> bool

    if isinstance(point, Point3D) and isinstance(parametric_rep[0], Point3D) and isinstance(parametric_rep[1], Vec3D):
        dims = 3
    elif isinstance(point, Point2D) and isinstance(parametric_rep[0], Point2D) and isinstance(parametric_rep[1], Vec2D):
        dims = 2
    else:
        raise RuntimeError("Bad passed arguments!")

    p1, line_dir = parametric_rep

    # Try to solve for parameters which put the point on the line. 
    diff = point.rep - p1.rep
    params = []
    for i in range(dims):
        params.append(robust_float_div(diff[i], line_dir.rep[i]))

    # Zeros can cause problems.
    # Consider a line defined by: (0, 0, 0) + t * <0, -5, 0>. Say we want
    #    to check if (0, 5, 0) lies on this line (it obviously does).
    # Using the equations above, the parameters could be nan, -1, and nan.
    # This implies that checking for parameter equality is insufficient. We
    #    should instead check if any of the parameters gives a valid solution.

    degenerate_cases = (float("+inf"), float("-inf"), float("nan"))

    for param in params:
        if param not in degenerate_cases:
            p = p1.rep + param * line_dir.rep
            res = Point3D(p)

            if identical_points(point, res):
                return True

    return False



# Check if a point lies in a box defined by two other points. 
# 
# Notes:
#    There are a number of degenerate cases. The point could lie exactly on the
#       other two points. Two points in three dimensions might only define a 
#       plane. Two points in two dimensions might only define a line. If the point
#       lies or in the geometric structure (point, line, plane, box) defined by the
#       two points, this function returns True by convention. 
# 
# Arguments:
#    point      - Point2D object or Point3D object. 
#    box_points - Two element tuple of Point2D or Point3D objects.
# 
# Returns:
#    Bool.
def point_in_box(point, box_points):
# type: (Union[Point2D, Point3D], Tuple[Union[Point2D, Point3D], Union[Point2D, Point3D]]) -> bool

    if isinstance(point, Point2D):
        if min(box_points[0].rep[0], box_points[1].rep[0]) <= point.rep[0] <= max(box_points[0].rep[0], box_points[1].rep[0]) and \
           min(box_points[0].rep[1], box_points[1].rep[1]) <= point.rep[1] <= max(box_points[0].rep[1], box_points[1].rep[1]):
            return True
        else:
            return False

    else:
        if min(box_points[0].rep[0], box_points[1].rep[0]) <= point.rep[0] <= max(box_points[1].rep[0], box_points[0].rep[0]) and \
           min(box_points[0].rep[1], box_points[1].rep[1]) <= point.rep[1] <= max(box_points[1].rep[1], box_points[0].rep[1]) and \
           min(box_points[0].rep[2], box_points[1].rep[2]) <= point.rep[2] <= max(box_points[1].rep[2], box_points[0].rep[2]):
            return True
        else:
            return False



# Check if two floats are equal using some fixed epsilon.
# 
# Notes:
#    None.
# 
# Arguments:
#    a - Float.
#    b - Float.
#
# Returns:
#    Bool. 
def float_equals(a, b):
# type: (float, float) -> bool

    if abs(a - b) <= 10**(-6):
        return True

    return False



# Check if two points are identical. 
# 
# Notes:
#    Not checking for being truly identical, only to float precision.
# 
# Arguments:
#    pt1 - Point3D or Point2D.
#    pt2 - Point3D or Point2D.
#
# Returns:
#    Bool. 
def identical_points(pt1, pt2):
# type: (Union[Point3D, Point2D], Union[Point3D, Point2D]) -> bool

    if isinstance(pt1, Point3D):
        if float_equals(pt1.rep[0], pt2.rep[0]) and float_equals(pt1.rep[1], pt2.rep[1]) and \
           float_equals(pt1.rep[2], pt2.rep[2]):
            return True
    else:
        if float_equals(pt1.rep[0], pt2.rep[0]) and float_equals(pt1.rep[1], pt2.rep[1]):
            return True

    return False



# In a sequence of points, find some pair of points which are not identical. 
# 
# Notes:
#    Not checking for being truly identical, only to float precision.
# 
# Arguments:
#    points - A tuple of Point3D or Point2D objects. 
#
# Returns:
#    Bool. 
def find_non_identical_points(points):
# type: (Tuple[Union[Point3D, Point2D], ...]) -> Tuple(Union[Point3D, Point2D], Union[Point3D, Point2D])

    for point_pair in itertools.combinations(points, 2):
        if not identical_points(*point_pair):
            return point_pair

    raise RuntimeError("Failed to find a pair of non identical points!")



def robust_float_div(a, b):
# type: (numbers.Real, numbers.Real) -> float

    # Convert to Python's float.
    a = float(a)
    b = float(b)

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




# Check if points are collinear.
def are_collinear(points):
# type: (Tuple[Union[Point3D, Point2D], ...]) -> bool

    assert(len(points) > 2)

    points_for_line = find_non_identical_points(points)

    if isinstance(points[0], Point3D):
        line = InfiniteLine3D(*points_for_line)
    else:
        line = InfiniteLine2D(*points_for_line)

    for point in points:
        if not line.is_on(point):
            return False

    return True



# Given some points, find the coefficients a, b, c, and d so that all the points
#    lie on the plane described by ax + by + cz = d. If the points are
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
    coords = np.stack([point1.rep, point2.rep, point3.rep], axis=0)
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

    assert(len(points) >= 3)

    # Find the coefficients of the plane which passes through the points.
    coeffs = get_plane_coeffs(points[0], points[1], points[2])

    # Check if all the points lie on this plane.
    return on_particular_plane(points, coeffs)



# Check if all points lie on a plane described by the coefficients of the plane.
# The coefficent order should be a, b, c, d where the equation of the plane is
#    ax + by + cz = d.
def on_particular_plane(points, coeffs):
# type: (List[Point3D], np.array) -> bool

    for point in points:
        if not float_equals(np.dot(coeffs[0:3], point.rep), coeffs[3]):
            return False

    return True



# Check if a set of points all lie on the same plane as an ngon. 
def on_plane_of_ngon(points, ngon):
# type: (List[Point3D], NGon3D) -> bool

    ngon_plane_coeffs = ngon.get_plane_coeffs()

    return on_particular_plane(points, ngon_plane_coeffs)



# Check if two vectors are orthogonal.
def are_orthogonal(vec1, vec2):
# type: (Vec3D, Vec3D) -> bool

    if float_equals(np.dot(vec1.rep, vec2.rep), 0):
        return True
    return False



# Find the centroid of a group of points.
#
# Notes:
#    None.
#
# Arguments:
#    points - Tuple of Point3D objects.
#
# Returns:
#    Point3D object.
def find_centroid(points):
# type: (Tuple[Point3D, ...]) -> Point3D 

    rep = np.array([0, 0, 0], dtype=float)
    for point in points:
        rep += point.rep
    cnt = len(points)
    return Point3D(rep / cnt)



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                                  Deprecated!
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# Check if all the points lie in the ngon.
def points_in_ngon_3D(points, ngon):
# type: (Tuple[Point3D], NGon3D) -> bool

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
    if not on_plane_of_ngon([point], ngon):
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
# If the point lies on one of the edges of the ngon, it is deemed to be in the 
#    ngon.
def point_in_ngon_2D(point, ngon):
# type: (Point2D, NGon2D) -> bool    

    # If the point lies on an edge, it is in the ngon.
    for idx, vertex in enumerate(ngon.vertices):
        
        if idx == len(ngon.vertices) - 1:
            next_vertex = ngon.vertices[0]
        else:
            next_vertex = ngon.vertices[idx + 1]
        
        # WARNING: Relies on ordering according to connectivity. 
        line = FiniteLine2D(vertex, next_vertex)

        if line.is_on(point):
            return True
    
    # Otherwise, check if it's in the interior. 
    ngon = path.Path(ngon.get_builtin_rep(), closed=True)
    return ngon.contains_point(point.components())



# Find the extrema (i.e. maximum and minimum) x values, y values, and z values.
# 
# Notes:
#    None.
# 
# Arguments:
#    points - List of Point3D objects.
#
# Returns:
#    6-tuple of floats.
def find_extrema(points):
# type: (List[Point3D]) -> Tuple[float, float, float, float, float, float]

    assert(len(points) > 0)

    max_x = points[0].rep[0]
    min_x = points[0].rep[0] 
    max_y = points[0].rep[1] 
    min_y = points[0].rep[1] 
    max_z = points[0].rep[2] 
    min_z = points[0].rep[2]

    for point in points:
        x = point.rep[0]
        y = point.rep[1]
        z = point.rep[2]

        max_x = x if x > max_x else max_x 
        min_x = x if x < min_x else min_x

        max_y = y if y > max_y else max_y 
        min_y = y if y < min_y else min_y

        max_z = z if z > max_z else max_z 
        min_z = z if z < min_z else min_z

    return (max_x, min_x, max_y, min_y, max_z, min_z)