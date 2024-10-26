// Standard library.
#include <iostream>
#include <filesystem>
#include <string>

// Third party.
#include "estimation.hxx"
#include "stress.hxx"

using namespace std;

int main()
{
    // REAL WORKPIECE SCANS + REAL TOOLPATH MOVEMENTS
    /*
    const string scan_path_prefix {"/Users/andrewlarson/Desktop/umich/RS_TP_Adaptation/RS_TP_Adaptation_physical_experiments/bending_bar_experiment/cleaned_scans/translated/"};
    const vector<filesystem::path> machining_scans {
                                                     scan_path_prefix + "pre_machining.stl",
                                                     // scan_path_prefix + "inc1.stl",
                                                     // scan_path_prefix + "inc2.stl",
                                                     // BAD: scan_path_prefix + "inc3.stl"
                                                     // scan_path_prefix + "inc4.stl",
                                                     // BAD: scan_path_prefix + "inc5.stl"
                                                   };
    
    const string toolpath_path_prefix {"/Users/andrewlarson/Desktop/umich/RS_TP_Adaptation/RS_TP_Adaptation_physical_experiments/bending_bar_experiment/gcode/cleaned_toolpath_surface_meshes/"};
    const vector<filesystem::path> toolpaths {
                                               toolpath_path_prefix + "move_0_linear_and_arc_of_circle.stl",
                                               toolpath_path_prefix + "move_1_arc_of_circle.stl",
                                               toolpath_path_prefix + "move_2_arc_of_circle.stl",
                                               toolpath_path_prefix + "move_3_arc_of_circle.stl",
                                               toolpath_path_prefix + "move_4_arc_of_circle.stl"
                                             };
    */

    // FAKE WORKPIECE SCANS + FAKE TOOLPATH MOVEMENTS 
    const string scan_path_prefix {"/Users/andrewlarson/Desktop/umich/RS_TP_Adaptation/RS_TP_Adaptation/test/"};
    const vector<filesystem::path> machining_scans {
                                                     scan_path_prefix + "pre_machining.stl",
                                                     scan_path_prefix + "inc1_test.stl",
                                                     scan_path_prefix + "inc2_test.stl",
                                                     scan_path_prefix + "inc3_test.stl",
                                                     scan_path_prefix + "inc4_test.stl",
                                                   };
    
    const string toolpath_path_prefix {"/Users/andrewlarson/Desktop/umich/RS_TP_Adaptation/RS_TP_Adaptation/test/"};
    const vector<filesystem::path> toolpaths {
                                               toolpath_path_prefix + "move0_test.stl",
                                               toolpath_path_prefix + "move1_test.stl",
                                               toolpath_path_prefix + "move2_test.stl",
                                               toolpath_path_prefix + "move3_test.stl",
                                             };

    RSEstimator estimator {machining_scans, toolpaths};
    estimator.estimate({0, 1});

    return EXIT_SUCCESS;
}
