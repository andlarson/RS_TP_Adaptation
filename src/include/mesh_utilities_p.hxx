#pragma once

// Standard library.
#include <filesystem>

// Library private.
#include "cgal_kernels_switch_p.hxx"

bool is_watertight(const SurfaceMesh& sm);
bool has_self_intersections(const SurfaceMesh& sm);
