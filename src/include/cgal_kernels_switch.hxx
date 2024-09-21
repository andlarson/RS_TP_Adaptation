#pragma once

#include "CGAL/Surface_mesh.h"

#include "CGAL/Simple_cartesian.h"
using Kernel = CGAL::Simple_cartesian<double>;
using Point = Kernel::Point_3;
using SurfaceMesh = CGAL::Surface_mesh<Point>;
using vertex_descriptor = SurfaceMesh::Vertex_index;
