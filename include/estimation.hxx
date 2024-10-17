#pragma once

// Standard library.
#include <vector>
#include <string>
#include <filesystem>
#include <memory>

// Library public.
#include "stress.hxx"

class RSEstimator 
{
    // PImpl idiom. 
    class RSEstimator_Impl;
    std::unique_ptr<RSEstimator_Impl> pimpl;

public:
    RSEstimator(const std::vector<std::filesystem::path>& scans,
                const std::vector<std::filesystem::path>& tool_paths);

    ~RSEstimator();
    RSEstimator(RSEstimator&&) noexcept;
    RSEstimator(const RSEstimator&) = delete;
    RSEstimator& operator=(RSEstimator &&) noexcept;
    RSEstimator& operator=(const RSEstimator &) = delete;

    StressTensor estimate(const std::pair<unsigned int, unsigned int>& estimation_interval);
};
