"""
This module contains functionality for doing various geometric operations of
    interest.
"""

import itertools
from typing import Any, Sequence, TypeVar, Iterable

import numpy as np
import numpy.linalg as linalg
import matplotlib.path as path

from src.util.debug import *


class Point2D:
    
    def __init__(self, arr: Any) -> None:
        """Creates a point in two dimensions.

           Args:
               arr: Numpy ndarray of length 2. The coordinates of the point.

           Returns:
               None.

           Raises:
               None.
        """

        assert isinstance(arr, np.ndarray)
        assert arr.size == 2
        
        arr = np.array(arr, dtype=float)
        self.rep = arr


    def components(self) -> tuple[Any, Any]:
        """Returns the components of the point."""

        return (self.rep[0], self.rep[1])



class Point3D:

    def __init__(self, arr: Any) -> None:
        """Creates a point in three dimensions.

           Args:
               arr: Numpy ndarray of length 3. The coordinates of the point.

           Returns:
               None.

           Raises:
               None.
        """

        assert isinstance(arr, np.ndarray)
        assert arr.size == 3
        
        arr = np.array(arr, dtype=float)
        self.rep = arr


    def proj_xy(self) -> Point2D:
        """Returns the projection of the point onto the x-y plane."""

        return Point2D(self.rep[0:2]) 


    def proj_xz(self) -> Point2D:
        """Returns the projection of the point onto the x-z plane."""

        return Point2D(np.array([self.rep[0], self.rep[2]]))


    def components(self) -> tuple[Any, Any, Any]:
        """Returns the components of the point."""

        return (self.rep[0], self.rep[1], self.rep[2])

    def __str__(self) -> str:
        """Overloads string conversion."""
        return str(self.components())


Point = TypeVar("Point", Point3D, Point2D)


class Vec2D:

    def __init__(self, arr: Any) -> None:
        """Creates a vector in two dimensions.

           Args:
               arr: Numpy ndarray of length 2. The coordinates of the point.

           Returns:
               None.

           Raises:
               None.
        """

        assert isinstance(arr, np.ndarray)
        assert arr.size == 2
        
        arr = np.array(arr, dtype=float)
        self.rep = arr



    def is_unit(self) -> bool:
        """Checks if the vector has unit norm."""

        if float_equals(self.len(), 1):
            return True
        return False


    def normalize(self) -> "Vec2D":
        """Returns a normalized version of the vector. Does not modify internal
               state."""

        res = (1 / linalg.norm(self.rep)) * self.rep
        return Vec2D(res)


    def len(self) -> Any:
        """Returns the length of the vector."""

        return linalg.norm(self.rep)


    def components(self) -> tuple[float, float]:
        """Returns the components of the vector."""

        return (self.rep[0], self.rep[1])



class Vec3D:

    def __init__(self, arr: Any) -> None:
        """Creates a vector in three dimensions.

           Args:
               arr: Numpy ndarray of length 3. The coordinates of the point.

           Returns:
               None.

           Raises:
               None.
        """

        assert isinstance(arr, np.ndarray)
        assert arr.size == 3
        
        arr = np.array(arr, dtype=float)
        self.rep = arr


    def is_unit(self) -> bool:
        """Checks if the vector has unit norm."""
        
        if float_equals(self.len(), 1):
            return True
        return False


    def normalize(self) -> "Vec3D":
        """Returns a normalized version of the vector. Does not modify internal
               state."""

        res = (1 / linalg.norm(self.rep)) * self.rep
        return Vec3D(res)


    def len(self) -> Any:
        """Returns the length of the vector."""

        return linalg.norm(self.rep)
 
    
    def get_orthonormal(self) -> "Vec3D":
        """Returns a vector orthogonal to the vector. Note that there are an infinite 
               number of vectors orthogonal to a vector. Returns one chosen arbitrarily 
               with unit length."""

        if self.rep[1] == 0 and self.rep[2] == 0:
            other_vec = np.array([0, 1, 0])
        else:
            other_vec = np.array([1, 0, 0])

        orth_vec = Vec3D(np.cross(other_vec, self.rep))

        return orth_vec.normalize()


    def components(self) -> tuple[float, float, float]:
        """Returns the components of the vector."""

        return (self.rep[0], self.rep[1], self.rep[2])

    
    def __sub__(self, other: "Vec3D") -> "Vec3D":
        """Overloads subtraction."""
        op1 = self.components()
        op2 = other.components()
        return Vec3D(np.array((op1[0] - op2[0], op1[1] - op2[1], op1[2] - op2[2])))


    def __str__(self) -> str:
        """Overloading string representation."""
        return str(self.components())
    


Vec = TypeVar("Vec", Vec3D, Vec2D)


class Basis3D:

    def __init__(self, v1: Vec3D, v2: Vec3D, v3: Vec3D) -> None:
        """Creates a basis which consists of three vectors which are of unit
               length and are orthogonal to one another."""

        assert v1.is_unit()
        assert v2.is_unit()
        assert v3.is_unit()
        assert _are_orthogonal(v1, v2)
        assert _are_orthogonal(v2, v3)
        assert _are_orthogonal(v1, v3)

        self.rep = np.stack([v1.rep, v2.rep, v3.rep], axis=0) 



class CSys3D:

    def __init__(self, translation: Vec3D, basis: Basis3D):
        """Describes a coordinate system relative to some global coordinate system.
        
           It is expected that all coordinate systems are described relative to
               a single global coodinate system.
        
           Args:
               translation: Describes the translation in space of coordinate system.
               basis:       Describes the directions of the axes. 
        
           Returns:
               None.

           Raises:
               None.
        """

        self.translation = translation
        self.basis = basis 

    
    def map_into_csys(self, points: list[Point3D]) -> list[Point3D]:
        """Maps points into the coordinate system described by this object.
        
           This function does a generic mapping. It has no knowledge about the
               coordinate system that the points already live in.
        
           Args:
               points: List of Point3D objects.
           
           Returns:
               The points in the new coordinate system.

           Raises:
               None.
        """

        return [Point3D(np.matmul(self.basis.rep, (point.rep - self.translation.rep))) for point in points]



class PlanarCubicC2Spline3D:

    def __init__(self, points: list[Point3D]) -> None:
        """Creates a planar cubic spline with continuous first and second derivatives 
               in 3D.
        
           Abaqus always uses the a cubic spline with continuous first and second
               derivatives to contruct wire features. For this reason, this is a very
               natural geometry.
           For simplicity, we force the spline to be planar for now. In particular,
               the points MUST ALL HAVE THE SAME Y COORDINATE!
        
           Args:
               points: The points which make up the spline. The ordering of the list 
                           should represent the connectivity. 
        
           Returns:
               None. 
           
           Raises:
               None.
        """
        assert len(points) > 1

        y = points[0].rep[1]
        for point in points:
            if not float_equals(y, point.rep[1]):
                raise AssertionError("Points are not planar!")

        self.y = y
        self.v_list = points



class NGon2D:

    def __init__(self, vertices: list[Point2D]) -> None:
        """Constructs a polygon in two dimensions.

           No checks are done to insure that the polygon is well defined.

           Args:
               vertices: The vertices ordered according to connectivity. In
                             other words, it is assumed that a line segment
                             connects the first vertex to the second, the
                             second to the third, ..., and the last vertex
                             to the first.

           Returns:
               None.

           Raises:
               None.
        """

        # TODO: We should really check that an ngon is well defined.
        #       In particular, that none of the points are on the interior of the ngon
        #           and that there are no redundant points.
        self.vertices = vertices 


    def get_builtin_rep(self) -> list[tuple[Any, Any]]:
        """Returns a representation of the polygon composed of only Python
               built-in types. Useful for interacting with other libraries."""

        return [(vertex.rep[0], vertex.rep[1]) for vertex in self.vertices]



class NGon3D:

    def __init__(self, vertices: list[Point3D]) -> None:
        """Constructs a polygon in three dimensions.

           No checks are done to insure that the polygon is well defined.

           Args:
               vertices: The vertices ordered according to connectivity. In
                             other words, it is assumed that a line segment
                             connects the first vertex to the second, the
                             second to the third, ..., and the last vertex
                             to the first.
                          The vertices must all exist on a single plane.

           Returns:
               None.

           Raises:
               None.
        """

        if not _on_any_plane(vertices):
            raise RuntimeError("The points don't describe a valid n-gon!")

        # TODO: We should really check that an ngon is well defined.
        #       In particular, that none of the points are on the interior of 
        #           the ngon and that there are no redundant points.
        self.vertices = vertices
    


    def get_plane_coeffs(self) -> Any:
        """Gets the coefficients of the the plane that the vertices lie on.

           Args:
               None.

           Returns:
               Numpy ndarray of length 3. The coefficients of the plane.

           Raises:
               None.
        """

        return _get_plane_coeffs(self.vertices[0], self.vertices[1], self.vertices[2])

   
    def proj_xy(self) -> NGon2D:
        """Projects the ngon onto the xy-plane, returning a two dimensional polygon.
    
           If the ngon existed on a plane orthogonal to the xy-plane, then all the
               points in the returned polygon will be collinear.

           Args:
               None.

           Returns:
               A representation of the two dimensional polygon resulting from
                   the projection.

           Raises:
               None.
        """

        proj_points = [vertex.proj_xy() for vertex in self.vertices]
        return NGon2D(proj_points)


    def proj_xz(self) -> NGon2D:
        """Projects the ngon onto the xz-plane, returning a two dimensional polygon.
    
           If the ngon existed on a plane orthogonal to the xz-plane, then all the
               points in the returned polygon will be collinear.

           Args:
               None.

           Returns:
               A representation of the two dimensional polygon resulting from
                   the projection.

           Raises:
               None.
        """

        proj_points = [vertex.proj_xz() for vertex in self.vertices]
        return NGon2D(proj_points)


    # Get a representation which is only composed of Python built-in types. 
    # Useful for passing to other libraries.
    def get_builtin_rep(self) -> list[tuple[Any, Any, Any]]:
        """Returns a representation of the polygon composed of only Python
               built-in types. Useful for interacting with other libraries."""

        return [(vertex.rep[0], vertex.rep[1], vertex.rep[2]) for vertex in self.vertices]



class _InfiniteLine:

    def __init__(self, p1: Point, p2: Point) -> None:
        """Creates an infinitely long line. 

           Args:
               p1: One point defining the line.
               p2: Another point defining the line.

           Returns:
               None.

           Raises:
               None.
        """

        # Store a parametric representation of the line.
        self.p1 = p1

        if isinstance(p1, Point3D):
            self.line_dir = Vec3D(p1.rep - p2.rep)
        elif isinstance(p1, Point2D):
            self.line_dir = Vec2D(p1.rep - p2.rep)
        else:
            raise AssertionError("Shouldn't get here.")


    def is_on(self, point: Point) -> bool:
        """Checks if a point is on the line."""

        if type(point) is not type(self.p1):
            raise RuntimeError("Dimensions don't match!")

        return _point_on_inf_line(point, (self.p1, self.line_dir))



def _point_on_inf_line(point: Point, parametric_rep: tuple[Point3D, Vec3D] | tuple[Point2D, Vec2D]) -> bool:
    """Checks if a point is on a line which has infinite length.
    
       Args:
           point:          The point of interest in two dimensions or three dimensions. 
           parametric_rep: A tuple containing a point and direction in two 
                               dimensions or three dimensions. The parametric 
                               representation of the line.
       
       Returns:
           Indication of point being on the line.
       
       Raises:
           None.
    """

    if isinstance(point, Point3D) and isinstance(parametric_rep[0], Point3D) and isinstance(parametric_rep[1], Vec3D):
        dims = 3
    elif isinstance(point, Point2D) and isinstance(parametric_rep[0], Point2D) and isinstance(parametric_rep[1], Vec2D):
        dims = 2
    else:
        raise TypeError("Dimensionality mismatch!")

    p1, line_dir = parametric_rep

    # Try to solve for parameters which put the point on the line. 
    diff = point.rep - p1.rep
    params = []
    for i in range(dims):
        params.append(_robust_float_div(diff[i], line_dir.rep[i]))

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

            if dims == 3:
                res = Point3D(p)
            else:
                res = Point2D(p)

            if _identical_points(point, res):
                return True

    return False



FLOATING_POINT_TOLERANCE = 10**(-6)
def float_equals(a: float, b: float) -> bool:
    """Checks if two floats are equal up to some epsilon."""

    if abs(a - b) <= FLOATING_POINT_TOLERANCE:
        return True

    return False



def seq_points_equal(s1: Sequence[Point], s2: Sequence[Point]) -> bool:
    """Checks if two sequences of points are the same.

       The ordering of the sequences doesn't matter. If, for every point in one
           sequence, there is an equivalnet point in the other sequence, the
           two sequences are considered equal.

       Args:
           s1: Some points.
           s2: Some points.

       Returns:
           Indication of the sequences being equal.

       Raises:
           None.
    """

    if len(s1) != len(s2):
        return False

    found_match = [False] * len(s2)
    for p1 in s1:
        for index, p2 in enumerate(s2):
            if _identical_points(p1, p2) and not found_match[index]:
                found_match[index] = True 
                break

    if False not in found_match:
        return True

    return False
                


def _identical_points(pt1: Point, pt2: Point) -> bool:
    """Checks if two points are identical within float tolerance."""

    if isinstance(pt1, Point3D):
        if float_equals(pt1.rep[0], pt2.rep[0]) and float_equals(pt1.rep[1], pt2.rep[1]) and \
           float_equals(pt1.rep[2], pt2.rep[2]):
            return True
    else:
        if float_equals(pt1.rep[0], pt2.rep[0]) and float_equals(pt1.rep[1], pt2.rep[1]):
            return True

    return False



def _find_non_identical_points(points: Iterable[Point]) -> tuple[Point, Point]:
    """Finds and returns a pair of points which are not identical within float 
           tolerance.

       Args:
           points: A sequence of points. 

       Returns:
           None.

       Raises:
           RuntimeError: A pair of non identical points could not be found. 
    """

    for point_pair in itertools.combinations(points, 2):
        if not _identical_points(*point_pair):
            return point_pair

    raise RuntimeError("Failed to find a pair of non identical points!")



def _robust_float_div(a: float, b: float) -> float:
    """Divides two numbers with the guarantee that any divide by zero exception
           will be caught and the result will be passed to the caller."""

    a = float(a)
    b = float(b)

    try:
        res = a / b 
    except ZeroDivisionError:
        if a == 0:
            res = float("nan")
        elif a > 0:
            res = float("+inf")
        else:
            res = float("-inf")

    return res



def _are_collinear(points: Sequence[Point]) -> bool:
    """Checks if a sequence of points are collinear. 

       Args:
           points: At least two points.

       Returns:
           Indication of collinearity.

       Raises:
           None.
    """

    assert len(points) > 2

    points_for_line = _find_non_identical_points(points)

    line = _InfiniteLine(*points_for_line)

    for point in points:
        if not line.is_on(point):
            return False

    return True



def _get_plane_coeffs(point1: Point3D, point2: Point3D, point3: Point3D) -> Any:
    """Finds the coefficients of the plane passing through three points.

       The plane is assumed to be of the form: ax + by + cz = d.

       Args:
           point1, ..., point3: Points.

       Returns:
           Numpy ndarray of length 4. The coefficients a, b, c, and d.

       Raises:
           RuntimeError: The points are collinear.
    """

    if _are_collinear([point1, point2, point3]):
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



def _on_any_plane(points: Sequence[Point3D]) -> bool:
    """Checks if three or more points lie on the same plane."""

    assert len(points) >= 3

    # Find the coefficients of the plane which passes through the points.
    coeffs = _get_plane_coeffs(points[0], points[1], points[2])

    # Check if all the points lie on this plane.
    return _on_particular_plane(points, coeffs)



def _on_particular_plane(points: Sequence[Point], coeffs: Any) -> bool:
    """Checks if a sequence of points lie on a plane of interest, within float
           tolerance.

       Args:
           points: Some points of interest.
           coeffs: Numpy ndarray of length 4. Assumed to contain a, b, c, d where
                       the plane is described by ax + by + cz = d.

       Returns:
           Indication of the points lying on the specified plane.

       Raises:
           None.
    """

    for point in points:
        if not float_equals(np.dot(coeffs[0:3], point.rep), coeffs[3]):
            return False

    return True



def on_plane_of_ngon(points: Sequence[Point], ngon: NGon3D) -> bool:
    """Checks if some points lie on the plane of a polygon."""

    ngon_plane_coeffs = ngon.get_plane_coeffs()

    return _on_particular_plane(points, ngon_plane_coeffs)



def _are_orthogonal(vec1: Vec, vec2: Vec) -> bool:
    """Checks if two vectors are orthogonal within float tolerance."""

    if float_equals(np.dot(vec1.rep, vec2.rep), 0):
        return True
    return False



def find_centroid(points: Sequence[Point]) -> Point:
    """Finds the centroid of some points."""
    
    rep = points[0].rep.copy()
    for point in points[1:]:
        rep += point.rep
    cnt = len(points)

    if len(rep) == 3:
        return Point3D(rep / cnt)
    elif len(rep) == 2:
        return Point2D(rep / cnt)
    else:
        raise AssertionError("Shouldn't get here.")



def seq_points(points: Sequence[Sequence[float]]) -> Sequence[Point]:
    """Converts sequence of sequences of points to genuine point objects.

       Args:
           points: Sequence of at least one point, where each point is represented
                       by a sequence of floats. Assumed that all the sequences
                       representing points are of the same length.

       Returns:
           Genuine point objects. 

       Raises:
           None.
    """

    if not points:
        raise RuntimeError("Expected at least one point.")

    is_3d = len(points[0]) == 3
    is_2d = len(points[0]) == 2

    if not is_3d and not is_2d:
        raise RuntimeError("Expected 2D or 3D points.")
    
    if is_3d:
        l = [Point3D(np.array(point)) for point in points]
    elif is_2d:
        l = [Point2D(np.array(point)) for point in points]

    return l



def distance(p1: Point, p2: Point) -> float:
    """Computes the distance between two points."""

    p1_rep = p1.components()
    p2_rep = p2.components()

    if len(p1_rep) == 3:
        dx = p1_rep[0] - p2_rep[0]
        dy = p1_rep[1] - p2_rep[1]
        dz = p1_rep[2] - p2_rep[2]
        return (dx**2+dy**2+dz**2)**(.5)
    elif len(p1_rep) == 2:
        dx = p1_rep[0] - p2_rep[0]
        dy = p1_rep[1] - p2_rep[1]
        return (dx**2+dy**2)**(.5)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                                  DEPRECATED 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


def points_in_ngon_3D(points: Iterable[Point3D], ngon: NGon3D) -> bool:
    """DEPRECATED. Was used to determine the which face should be partitioned.
        
       Checks if all the points lie in a polygon. 
    """

    for point in points:
        if not _point_in_ngon_3D(point, ngon):
            return False

    return True



def _point_in_ngon_3D(point: Point3D, ngon: NGon3D) -> bool:
    """DEPRECATED. Was a helper function for checking if many points were in a
           polygon.
       
       Checks if a point lies in a polygon."""

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
    if _are_collinear(proj_ngon.vertices + [proj_point]):
        proj_point = point.proj_xz()
        proj_ngon = ngon.proj_xz()

    return _point_in_ngon_2D(proj_point, proj_ngon)



def _point_in_ngon_2D(point: Point2D, ngon: NGon2D) -> bool:
    """DEPRECATED. Was a helper function for 3D version of this function.
    
       Checks a if point lies inside a polygon. If the point lies on one of the
           edges of the polygon, it is considered inside the polygon."""

    # If the point lies on an edge, it is in the ngon.
    for idx, vertex in enumerate(ngon.vertices):
        
        if idx == len(ngon.vertices) - 1:
            next_vertex = ngon.vertices[0]
        else:
            next_vertex = ngon.vertices[idx + 1]
        
        # WARNING: Relies on ordering according to connectivity. 
        line = _FiniteLine(vertex, next_vertex)

        if line.is_on(point):
            return True
    
    # Otherwise, check if it's in the interior. 
    ngon = path.Path(ngon.get_builtin_rep(), closed=True)   
    return ngon.contains_point(point.components()) 



def find_extrema(points: Sequence[Point3D]) -> tuple[float, float, float, float, float, float]:
    """DEPRECATED. Was a helper function for bounding box creation.

       Find the extrema coordinates of a sequence of at least one points."""

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



class SpecRightRectPrism:
    """DEPRECATED. Was used for tool pass geometry before splines were introduced."""
   
    def __init__(self, v1: Point3D, v2: Point3D, v3: Point3D, v4: Point3D, 
                 v5: Point3D, v6: Point3D, v7: Point3D, v8: Point3D) -> None:
        """Creates a right rectangular prism which has faces which are paralle to
               the standard planes of the coordinate system.

           Args:
               v1, .., v8: The vertices. 

           Returns:
               None.

           Raises:
               None.
        """

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

        assert len(self.same_x_g1) == 4
        assert len(self.same_x_g2) == 4
        assert len(self.same_y_g1) == 4
        assert len(self.same_y_g2) == 4
        assert len(self.same_z_g1) == 4
        assert len(self.same_z_g2) == 4

        self.vertices = v_list

   
    def get_smaller_z(self) -> float:
        """Returns the value of the z coordinate corresponding to the plane of
               prism which has smallest z value."""

        return min(self.same_z_g1[0].rep[2], self.same_z_g2[0].rep[2]) 


    def get_dims(self) -> tuple[float, float, float]:
        """Returns the length in the x-direction, the width in the y-direction, and the
               height in the z-direction."""
       
        x_length = abs(self.same_x_g1[0].rep[0] - self.same_x_g2[0].rep[0])
        y_width = abs(self.same_y_g1[0].rep[1] - self.same_y_g2[0].rep[1])
        z_height = abs(self.same_z_g1[0].rep[2] - self.same_z_g2[0].rep[2])

        return (x_length, y_width, z_height)

    
    def get_rect_corners(self) -> tuple[Point3D, Point3D]:
        """Returns 2 vertices which have the same z coordinates but differing x and y
               coordinates."""

        v1 = self.vertices[0]

        for v in self.vertices[1:]:
            if (v1.rep[2] == v.rep[2]) and (v1.rep[0] != v.rep[0]) and (v1.rep[1] != v.rep[1]):
                return v1, v

        raise AssertionError("Cannot find two vertices with matching z coordinates.")


    def get_centroid(self) -> Point3D:
        """Returns the centroid of the prism.""" 
        
        avg_x = (self.same_x_g1[0].rep[0] + self.same_x_g2[0].rep[0]) / 2
        avg_y = (self.same_y_g1[0].rep[1] + self.same_y_g2[0].rep[1]) / 2
        avg_z = (self.same_z_g1[0].rep[2] + self.same_z_g2[0].rep[2]) / 2

        return Point3D(np.array([avg_x, avg_y, avg_z])) 



class _FiniteLine:

    def __init__(self, p1: Point, p2: Point) -> None:
        """DEPRECATED. Was used to help do point inside polygon checking.

           Creates a line of finite length. 

           The line only exists on and between the two points which define it.

           Args:
               p1: One point defining the line.
               p2: Another point defining the line.

           Returns:
               None.

           Raises:
               None.
        """

        self.p1 = p1 
        self.p2 = p2
        
        if isinstance(p1, Point2D):
            self.line_dir = Vec2D(np.array([p1.rep[0] - p2.rep[0], p1.rep[1] - p2.rep[1]]))
        elif isinstance(p1, Point3D):
            self.line_dir = Vec3D(np.array([p1.rep[0] - p2.rep[0], p1.rep[1] - p2.rep[1], p1.rep[2] - p2.rep[2]]))
        else:
            raise AssertionError("Shouldn't get here.")


    def is_on(self, point: Point) -> bool:
        """Checks if a point is on the line."""
        
        if not isinstance(point, type(self.p1)):
            raise TypeError("Dimensionality mismatch.")
        
        if _point_on_inf_line(point, (self.p1, self.line_dir)) and _point_in_box(point, (self.p1, self.p2)):
            return True
        return False



def _point_in_box(point: Point, box_points: tuple[Point, Point]) -> bool:
    """DEPRECATED. Was used as helper for checking if point lies on line of finite
           length.

       Checks if a point lies in a box defined by two other points.
       
       There are a number of degenerate cases. The point could lie exactly on 
           the two points defining the box. Two points in three dimensions might 
           only define a plane, not a rectangular prism. Two points in two 
           dimensions might only define a line, not a rectangle. 

       If the point lies in the geometric structure (point, line, plane, box) defined 
           by the two points, this function returns True by convention. 

       Args:
           point:      The point of interest.
           box_points: The points defining the box.
                       
           The dimensionality of the point and points defining the box should be the same.

       Returns:
           Indication if the points lie on or in the box.

       Raises:
           None.
    """

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



