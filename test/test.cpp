/*
    Tests for residual stress estimation technique.
*/

// Standard library.
#include <iostream>
#include <filesystem>

// Third party.
#include "estimation.hxx"
#include "stress.hxx"

int main()
{
    std::vector<std::filesystem::path> machining_scans {
                                                        "/Users/andrewlarson/Desktop/umich/RS_TP_Adaptation_notes/bending_bar_experiment/cleaned_scans/2manifold_and_watertight/pre_machining.stl",
                                                        // "/Users/andrewlarson/Desktop/umich/RS_TP_Adaptation_notes/bending_bar_experiment/cleaned_scans/2manifold_and_watertight/inc1.stl",
                                                        // "/Users/andrewlarson/Desktop/umich/RS_TP_Adaptation_notes/bending_bar_experiment/cleaned_scans/2manifold_and_watertight/inc2.stl",
                                                        // BAD: "/Users/andrewlarson/Desktop/umich/RS_TP_Adaptation_notes/bending_bar_experiment/cleaned_scans/2manifold_and_watertight/inc3.stl",
                                                        // "/Users/andrewlarson/Desktop/umich/RS_TP_Adaptation_notes/bending_bar_experiment/cleaned_scans/2manifold_and_watertight/inc4.stl",
                                                        // BAD: "/Users/andrewlarson/Desktop/umich/RS_TP_Adaptation_notes/bending_bar_experiment/cleaned_scans/2manifold_and_watertight/inc5.stl"
                                                       };

    std::vector<std::filesystem::path> tool_paths {
                                                   // BAD: "/Users/andrewlarson/Downloads/corner.stl"
                                                   "/Users/andrewlarson/Downloads/straight_line_scaled.stl",
                                                   // "/Users/andrewlarson/Downloads/zigzag.stl", 
                                                   // BAD: "/Users/andrewlarson/Downloads/horseshoe.stl" 
                                                  };

    RSEstimator estimator {machining_scans, tool_paths};
    estimator.estimate({0, 1});

    return EXIT_SUCCESS;
}
