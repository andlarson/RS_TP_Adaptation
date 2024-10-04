#pragma once

#include "CGAL/Surface_mesh.h"
// #include "CGAL/Simple_cartesian.h"
#include "CGAL/Exact_predicates_exact_constructions_kernel.h"

// using Kernel = CGAL::Simple_cartesian<double>;
using Kernel = CGAL::Exact_predicates_exact_constructions_kernel;
using Point = Kernel::Point_3;
using SurfaceMesh = CGAL::Surface_mesh<Point>;
