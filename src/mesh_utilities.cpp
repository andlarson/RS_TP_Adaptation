// Standard library.
#include <filesystem>
#include <iostream>

// Third party.
#include "CGAL/Polygon_mesh_processing/orientation.h"
#include "CGAL/Polygon_mesh_processing/self_intersections.h"

// Library private. 
#include "cgal_kernels_switch_p.hxx"
#include "mesh_utilities_p.hxx"

bool is_watertight(const SurfaceMesh& sm)
{
    if (CGAL::is_closed(sm))
        if (CGAL::Polygon_mesh_processing::does_bound_a_volume(sm))
            return true;
    return false;
}

bool has_self_intersections(const SurfaceMesh& sm)
{
    if (CGAL::Polygon_mesh_processing::does_self_intersect(sm))
        return true;
    return false;
}
