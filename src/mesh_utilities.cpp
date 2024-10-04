/*
    Provides general utilities for volumetric and surfacic meshes. 
*/

// Standard library.
#include <filesystem>
#include <iostream>

// Third party.
#include "CGAL/Surface_mesh.h"
#include "CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h"
#include "CGAL/Polygon_mesh_processing/orientation.h"
#include "CGAL/Polygon_mesh_processing/self_intersections.h"

// Library private. 
#include "cgal_kernels_switch_p.hxx"
#include "mesh_utilities_p.hxx"

// DEBUG!?
#include <vector>
#include <CGAL/Polygon_mesh_processing/self_intersections.h>
#include <iterator>
#include<CGAL/draw_surface_mesh.h>

bool is_watertight(const SurfaceMesh& sm)
{
    if (CGAL::is_closed(sm))
        if (CGAL::Polygon_mesh_processing::does_bound_a_volume(sm))
            return true;
    return false;
}

bool has_self_intersections(const SurfaceMesh& sm)
{
    // DEBUG!?
    /*
    using FaceDesc = boost::graph_traits<SurfaceMesh>::face_descriptor;
    std::vector<std::pair<FaceDesc, FaceDesc>> faces;
    CGAL::Polygon_mesh_processing::self_intersections(sm, std::back_inserter(faces));
    draw(sm);
    */

    if (CGAL::Polygon_mesh_processing::does_self_intersect(sm))
        return true;
    return false;
}
