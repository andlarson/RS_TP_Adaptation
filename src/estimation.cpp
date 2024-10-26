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
#include "CGAL/Graphics_scene.h"
#include "CGAL/Qt/Basic_viewer.h"

// *****************************************************************************
//                      RSEstimator: PImpl Forwarding
// *****************************************************************************

RSEstimator::RSEstimator(const std::vector<std::filesystem::path>& scans,
                         const std::vector<std::filesystem::path>& toolpaths)
    : pimpl{std::make_unique<RSEstimator_Impl>(scans, toolpaths)} 
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
    Main object for residual stress estimation.

    Requires:
        (1) The scans are watertight, 2-manifold, and have no self
                intersections.
        (2) There is one more scan than toolpath.
        (3) The toolpaths are watertight, 2-manifold, and have no self
                intersections.

    Assumes:
        (1) The scans and toolpaths exist in the same coordinate system and have
                the same units.
        (2) The scans are properly oriented.
        (3) The toolpaths are properly oriented.

    Note: This function does a small amount of repair and orientation.
        As such, it is not strictly required that the scans and toolpaths
        are watertight, 2-manifold, lack self-intersections, and are properly
        oriented. However, the repair and orientation that this function does
        are neither thorough nor robust. As such, the aforementioned
        requirements and assumptions really should be met.

    Arguments:
        scans:                Absolute paths to surface meshes, in .stl format,
                                  containing the scan data captured during
                                  machining.  
        toolpaths:            Absolute paths to surface meshes, in .stl format, 
                                  containing the path that the machine tool
                                  followed during machining. 
*/
RSEstimator::RSEstimator_Impl::RSEstimator_Impl(const std::vector<std::filesystem::path>& scans,
                                                const std::vector<std::filesystem::path>& toolpaths)
{
    assert(scans.size() == toolpaths.size() + 1);

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

    for (const auto& toolpath : toolpaths)
    {
        SurfaceMesh sm;
        assert(std::filesystem::exists(toolpath));
        assert(toolpath.extension() == STL_FILE_EXTENSION);
        assert(CGAL::Polygon_mesh_processing::IO::read_polygon_mesh(toolpath.string(), sm, CGAL::parameters::verbose(true)));
        assert(!has_self_intersections(sm));
        assert(is_watertight(sm));
        this->toolpaths.push_back(sm);
    }
    this->toolpath_files = toolpaths;
}

/*
    Estimates residual stress for a single estimation interval.    

    Requires:
        (1) If there are n scans, then the interval must be a subset of {0, 1,
                2, 3, ..., n-1} with cardinality 2.

    Assumes:
        (1) The interval is inclusive. An interval [0, 2] yields an estimate of
                residual stress using data from scan 0, scan 2, toolpath 0, and
                toolpath 1.
    
    Arguments:
        estimation_interval: Interval across which the residual stress estimate
                                 will be produced.
    
    Return:
        A single stress tensor, with some unknown components. Assuming the
            estimation interval is of the form {a, b}, then the estimate of
            stress is valid in the region of space occupied by toolpaths
            {a, a+1, a+2, a+3, ..., b-2, b-1, b}.
        In other words, the stress estimate is valid for all material removed
            in toolpath a, a+1, a+2, ..., b-1.
*/
StressTensor RSEstimator::RSEstimator_Impl::estimate(const std::pair<unsigned int, unsigned int>& estimation_interval)
{
    assert(estimation_interval.first < estimation_interval.second);
    assert(estimation_interval.first >= 0);
    assert(estimation_interval.second < scan_files.size());

    // DEBUG!? 
    SurfaceMesh diff {scans[0]};
    bool res {CGAL::Polygon_mesh_processing::corefine_and_compute_difference(diff, toolpaths[0], diff)};
    assert(res);
    res = CGAL::Polygon_mesh_processing::corefine_and_compute_difference(diff, toolpaths[1], diff);
    assert(res);
    res = CGAL::Polygon_mesh_processing::corefine_and_compute_difference(diff, toolpaths[2], diff);
    assert(res);
    res = CGAL::Polygon_mesh_processing::corefine_and_compute_difference(diff, toolpaths[3], diff);
    assert(res);
    
    CGAL::Graphics_scene scene;
    // CGAL::add_to_graphics_scene(scans[0], scene);
    // CGAL::add_to_graphics_scene(toolpaths[0], scene);
    CGAL::add_to_graphics_scene(diff, scene);
    CGAL::draw_graphics_scene(scene);
}
