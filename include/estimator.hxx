#pragma once

// Standard library.
#include <vector>
#include <string>
#include <filesystem>

// Library public.
#include "stress.hxx"

using Interval = std::pair<unsigned int, unsigned int>;
using StressEstimate = std::pair<StressTensor, std::string>;

class ResidualStressEstimator
{
public:
    ResidualStressEstimator(const std::vector<std::filesystem::path>& scans,
                            const std::vector<std::filesystem::path>& tool_paths);
    
    StressEstimate estimate(const Interval& estimation_interval);

private:
    std::vector<std::string> scans; 
    std::vector<std::string> tool_paths; 
};
