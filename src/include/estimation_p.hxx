#pragma once

// Standard library.
#include <vector>
#include <filesystem>

// Library public.
#include "stress.hxx"
#include "estimation.hxx"

// Library private.
#include "cgal_kernels_switch_p.hxx"

class RSEstimator::RSEstimator_Impl
{
    std::vector<std::filesystem::path> scan_files;
    std::vector<std::filesystem::path> tool_path_files;

    std::vector<SurfaceMesh> scans;
    std::vector<SurfaceMesh> tool_paths;

public:
    RSEstimator_Impl(const std::vector<std::filesystem::path>& scans,
                     const std::vector<std::filesystem::path>& tool_paths);
    StressTensor estimate(const std::pair<unsigned int, unsigned int>& estimation_interval);
};
