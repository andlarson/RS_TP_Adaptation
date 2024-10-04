/*
    Top level functionality to estimate residual stress during machining.
*/

// Standard library.
#include <vector>
#include <string>
#include <filesystem>
#include <assert.h>
#include <memory>

// Third party.
#include "CGAL/Polygon_mesh_processing/corefinement.h"
#include "CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h"

// Public includes.
#include "estimation.hxx"
#include "stress.hxx"

// Private includes.
#include "mesh_utilities_p.hxx"
#include "estimation_p.hxx"

// DEBUG!?
#include "CGAL/draw_surface_mesh.h"
#include "CGAL/Polygon_mesh_processing/self_intersections.h"
#include "CGAL/boost/graph/graph_traits_Surface_mesh.h"
#include <iterator>

// *****************************************************************************
//                      RSEstimator: PImpl Forwarding
// *****************************************************************************

RSEstimator::RSEstimator(const std::vector<std::filesystem::path>& scans,
                         const std::vector<std::filesystem::path>& tool_paths)
    : pimpl{std::make_unique<RSEstimator_Impl>(scans, tool_paths)} 
{
}

StressTensor RSEstimator::estimate(const std::pair<unsigned int, unsigned int>& estimation_interval)
{
    return pimpl->estimate(estimation_interval);
}

RSEstimator::~RSEstimator() = default;
RSEstimator::RSEstimator(RSEstimator&&) noexcept = default;
RSEstimator& RSEstimator::operator=(RSEstimator &&) noexcept = default;

// *****************************************************************************
//                  RSEstimator_Impl: PImpl Implementation 
// *****************************************************************************

/*
    Sets up to do residual stress estimation.

    Arguments:
        scans:                Absolute paths to surface meshes, in .stl format,
                                  containing the scan data captured during
                                  machining. Must be watertight and have
                                  correctly-oriented surface normals. 
                              The scans and toolpaths need to exist in the same
                                  coordinate system and have the same units.
                              There must be one more scan than toolpath.
        toolpaths:            Absolute paths to surface meshes, in .stl format, containing 
                                  the path that the machine tool followed
                                  during machining. Must be watertight and have
                                  correctly-oriented surface normals. 
    
    Return:
        None.
*/
RSEstimator::RSEstimator_Impl::RSEstimator_Impl(const std::vector<std::filesystem::path>& scans,
                                                const std::vector<std::filesystem::path>& tool_paths)
{
    // DEBUG!?
    /*
    assert(scans.size() == tool_paths.size() + 1);
    */

    const std::string STL_FILE_EXTENSION {".stl"};
    for (const auto& scan : scans)
    {
        SurfaceMesh sm;
        assert(std::filesystem::exists(scan));
        assert(scan.extension() == STL_FILE_EXTENSION);
        assert(CGAL::Polygon_mesh_processing::IO::read_polygon_mesh(scan.string(), sm, CGAL::parameters::verbose(true)));
        assert(!has_self_intersections(sm));
        assert(is_watertight(sm));
        this->scans.push_back(sm);
    }
    this->scan_files = scans;

    for (const auto& tool_path : tool_paths)
    {
        SurfaceMesh sm;
        assert(std::filesystem::exists(tool_path));
        assert(tool_path.extension() == STL_FILE_EXTENSION);
        assert(CGAL::Polygon_mesh_processing::IO::read_polygon_mesh(tool_path.string(), sm, CGAL::parameters::verbose(true)));
        assert(!has_self_intersections(sm));
        assert(is_watertight(sm));
        this->tool_paths.push_back(sm);
    }
    this->tool_path_files = tool_paths;
}

/*
    Estimates residual stress for a single estimation interval.    
    
    Arguments:
        estimation_interval: Interval across which the residual stress estimate
                                 will be produced.
                             If there are n scans, then the valid intervals are
                                 all subsets of [0, n-1] with cardinality two.
                             The interval is inclusive. An interval (0, 2) yields
                                 an estimate of residual stress using data from
                                 scan 0, scan 2, toolpath 0, and toolpath 1.
    
    Return:
        A single stress tensor, with some unknown components. Assuming the
            estimation interval is of the form [a, b], then the estimate of
            stress is valid in the region of space occupied by toolpaths
            [a, b].
        In other words, the stress estimate is valid for all material removed
            in toolpath a, a+1, a+2, ..., b-1.
*/
StressTensor RSEstimator::RSEstimator_Impl::estimate(const std::pair<unsigned int, unsigned int>& estimation_interval)
{
    // DEBUG!?
    /*
    assert(estimation_interval.first < estimation_interval.second);
    assert(estimation_interval.first >= 0);
    assert(estimation_interval.second < scan_files.size());
    */

    // DEBUG!? 
    SurfaceMesh diff;
    bool res {CGAL::Polygon_mesh_processing::corefine_and_compute_difference(scans[0], tool_paths[0], diff)};
    assert(res);
    draw(diff);
}
