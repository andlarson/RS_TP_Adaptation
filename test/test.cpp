// Standard library.
#include <filesystem>
#include <string>

// Third party.
#include "estimation.hxx"

using namespace std;

int main()
{
    const filesystem::path scan_path_prefix {string(TEST_DIR) + "/synthetic_scan_data/simple1/scans/"};
    const vector<filesystem::path> machining_scans {
                                                     scan_path_prefix.string() + "pre_machining.stl",
                                                     scan_path_prefix.string() + "inc1_test.stl",
                                                     scan_path_prefix.string() + "inc2_test.stl",
                                                     scan_path_prefix.string() + "inc3_test.stl",
                                                     scan_path_prefix.string() + "inc4_test.stl",
                                                   };
    
    const filesystem::path toolpath_path_prefix {string(TEST_DIR) + "/synthetic_scan_data/simple1/toolpaths/"};
    const vector<filesystem::path> toolpaths {
                                               toolpath_path_prefix.string() + "move0_test.stl",
                                               toolpath_path_prefix.string() + "move1_test.stl",
                                               toolpath_path_prefix.string() + "move2_test.stl",
                                               toolpath_path_prefix.string() + "move3_test.stl",
                                             };

    RSEstimator estimator {machining_scans, toolpaths};
    estimator.estimate({0, 1});

    return EXIT_SUCCESS;
}
