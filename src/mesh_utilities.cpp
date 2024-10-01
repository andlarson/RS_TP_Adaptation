/*
    Provides general utilities for volumetric and surfacic meshes. 
*/

// Standard library.
#include <filesystem>
#include <iostream>

// Third party.
#include "CGAL/Surface_mesh.h"
#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>
#include "CGAL/Polygon_mesh_processing/orientation.h"

// Library private. 
#include "cgal_kernels_switch.hxx"
#include "mesh_utilities.hxx"

/*
    Checks if a .off/.stl/.obj/etc. file contains a surface mesh that is close 
        enough to being watertight.
*/
bool is_watertight(const std::filesystem::path& surface_mesh)
{
    SurfaceMesh m;

    const bool res {CGAL::Polygon_mesh_processing::IO::read_polygon_mesh(surface_mesh.string(), m)};
    
    /*
    // DEBUG: Surface meshes that I know are 2-manifold cannot be processed by 
    //     the CGAL::IO::read_polygon_mesh() function (which doesn't do any
    //     automated repairs). When debugging this read_polygon_mesh() function
    //     call, I discovered that the failure occurs because an arbitrary face
    //     cannot be added to the surface mesh. The code below replicates the
    //     stuff that CGAL::IO::read_polygon_mesh() does under the hood and, when
    //     a face cannot be added, it dumps all the successfully added faces to
    //     a .stl file and the bad face to a different .stl file. For the particular
    //     bad face in the particular surface mesh that I was looking at, the
    //     face cannot be added because one of its vertices is not new and the
    //     halfedge associated with that vertex is not a border halfedge. I think
    //     there are workarounds for this problem, but I've opted to use the
    //     CGAL::Polygon_mesh_processing::IO::read_polygon_mesh() function because
    //     it does automated repairs that seem to work (i.e. produce a 2-manifold
    //     surface mesh) and it does not move any points and does not introduce
    //     any points. 
    std::vector<Point> vertices;
    std::vector<std::vector<std::size_t>> triangles;
    const bool res {CGAL::IO::read_STL(surface_mesh.string(), vertices, triangles)};
    
    std::vector<SurfaceMesh::Vertex_index> new_vertex_indices;
    for (const Point& p : vertices)
        new_vertex_indices.push_back(m.add_vertex(p));
    
    for (const auto& tri : triangles)
    {
        std::vector<SurfaceMesh::Vertex_index> tri_new_vertex_indices;
        for (const auto& v : tri)
            tri_new_vertex_indices.push_back(new_vertex_indices[v]);

        assert(tri_new_vertex_indices.size() == 3);

        SurfaceMesh::Face_index f {m.add_face(tri_new_vertex_indices)};
        if (f == SurfaceMesh::null_face())
        {
            m.add_face(tri_new_vertex_indices);

            std::filesystem::path good_stl {"/Users/andrewlarson/Downloads/good_stl.stl"};
            assert(CGAL::IO::write_STL(good_stl, m, CGAL::parameters::use_binary_mode(false)));
            
            SurfaceMesh m_bad;
            std::vector<SurfaceMesh::Vertex_index> bad_new_vertex_indices;
            for (const Point& p : vertices)
                bad_new_vertex_indices.push_back(m_bad.add_vertex(p));
            std::vector<SurfaceMesh::Vertex_index> tri_bad_new_vertex_indices;
            for (const auto& v : tri)
                tri_bad_new_vertex_indices.push_back(bad_new_vertex_indices[v]);
            assert(tri_new_vertex_indices.size() == 3);
            SurfaceMesh::Face_index f {m_bad.add_face(tri_bad_new_vertex_indices)};
            std::filesystem::path bad_stl {"/Users/andrewlarson/Downloads/bad_stl.stl"};
            assert(CGAL::IO::write_STL(bad_stl, m_bad, CGAL::parameters::use_binary_mode(false)));

            std::exit(EXIT_FAILURE);
        }
    }
    */
    
    if (res)
        if (CGAL::is_closed(m))
            if (CGAL::Polygon_mesh_processing::does_bound_a_volume(m))
                return true;
    return false;
}
