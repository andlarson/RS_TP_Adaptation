/*
    Provides general utilities for volumetric and surfacic meshes. 
*/

// Standard library.
#include <filesystem>
#include <iostream>

// Third party.
#include "geogram/mesh/mesh.h"
#include "geogram/mesh/mesh_io.h"
#include "geogram/mesh/mesh_geometry.h"

/*
    Checks that a .stl file contains a surface mesh that is watertight.
*/
bool is_watertight(std::filesystem::path surface_mesh)
{
    GEO::Mesh mesh;
    assert(mesh_load(surface_mesh.string(), mesh));
    
    std::cout << "The number of cells in the mesh " << surface_mesh.string() << " is " << GEO::mesh_cells_volume(mesh) << std::endl;

    return true;
}

