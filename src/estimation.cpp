/*
    Top level functionality to estimate residual stress during machining.
*/

// Standard library.
#include <vector>
#include <string>
#include <filesystem>

// Third party.

// Public includes.
#include "estimator.hxx"
#include "stress.hxx"

// Private includes.
#include "include/mesh_utilities.hxx"

/*
    Sets up to do residual stress estimation.

    Arguments:
        scans:                Absolute paths to surface meshes, in .stl format, containing the scan data
                                  captured during machining. Must be watertight and have
                                  correctly-oriented surface normals. The number of scans
                                  and toolpaths must be equal and all scans must be
                                  in the same coordinate system (i.e. must be
                                  aligned).
        toolpaths:            Absolute paths to surface meshes, in .stl format, containing the path that the
                                  machine tool followed during machining. Must be watertight and 
                                  have correctly-oriented surface normals. The number of scans
                                  and toolpaths must be equal and all toolpaths
                                  must be in the same coordinate system as the
                                  scans.
*/
ResidualStressEstimator::ResidualStressEstimator(const std::vector<std::filesystem::path>& scans,
                                                 const std::vector<std::filesystem::path>& tool_paths)
{
    for (const auto& scan : scans)
    {
        assert(std::filesystem::exists(scan));
        assert(is_watertight(scan));
    }

    for (const auto& tool_path : tool_paths)
    {
        assert(std::filesystem::exists(tool_path));
        assert(is_watertight(tool_path));
    }
}

/*
    Estimates residual stress for a single estimation interval.    
    
    Arguments:
        estimation_interval: Interval across which the residual stress estimate
                                 will be produced. If there are n surface meshes
                                 and n tool paths labeled (0, 1, 2, ..., n - 1),
                                 then every continuous interval of length 2 or
                                 greater contained in this range is a valid estimation
                                 interval.
    
    Return:
        A single estimate of residual stress. A residual stress estimate is 
            composed of a single stress tensor, with some known and unknown 
            components, located in some volume of space. 
*/
StressEstimate ResidualStressEstimator::estimate(const Interval& estimation_interval)
{

}
