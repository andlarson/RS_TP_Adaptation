/*
    Provides general utilities for volumetric and surfacic meshes. 
*/

// Standard library.
#include <filesystem>
#include <iostream>

// Third party.
#include "CGAL/Surface_mesh.h"
#include "CGAL/boost/graph/IO/polygon_mesh_io.h"
#include "CGAL/boost/graph/IO/STL.h"
// #include "CGAL/IO/STL.h"

// Library private. 
#include "cgal_kernels_switch.hxx"
#include "mesh_utilities.hxx"

/*
    Checks if a .stl file contains a surface mesh that is watertight.
*/
bool is_watertight(const std::filesystem::path& surface_mesh)
{
    SurfaceMesh m;

    // Check for 2-manifoldness. Note that 2-manifoldness != watertightness == 3-manifoldness.
    const bool res {CGAL::IO::read_STL(surface_mesh.string(), m, CGAL::parameters::verbose(true))};
    
    /*
    std::vector<Point> points;
    std::vector<std::vector<std::size_t>> triangles;
    const bool res {CGAL::IO::read_STL(surface_mesh.string(), points, triangles)};
    */

    if (res)
        if (cgal_surface_mesh_is_watertight(m))
            return true;
    return false;
}

/*
    Checks if a CGAL Surface Mesh is watertight.

    Arguments:
        m: A surface mesh that must be 2-manifold.
*/
bool cgal_surface_mesh_is_watertight(const SurfaceMesh& m)
{
    for (const vertex_descriptor& v : m.vertices())
        if (m.is_border(v))
            return false;
    return true;
}
