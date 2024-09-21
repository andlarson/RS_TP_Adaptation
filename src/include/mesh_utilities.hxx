#pragma once

// Standard library.
#include <filesystem>

// Library private.
#include "cgal_kernels_switch.hxx"

bool is_watertight(const std::filesystem::path& surface_mesh);
bool cgal_surface_mesh_is_watertight(const SurfaceMesh& m);
