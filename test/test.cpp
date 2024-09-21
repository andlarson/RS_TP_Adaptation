/*
    Tests for residual stress estimation technique.
*/

// Standard library.
#include <iostream>

// Third party.
#include "estimator.hxx"
#include "stress.hxx"

int main()
{
    std::vector<std::filesystem::path> machining_scans {// "/Users/andrewlarson/Downloads/corner.stl", 
                                                        "/Users/andrewlarson/Downloads/straight_line.stl" 
                                                       };

    std::vector<std::filesystem::path> tool_paths {"/Users/andrewlarson/Downloads/zigzag.stl", 
                                                   "/Users/andrewlarson/Downloads/horseshoe.stl" 
                                                  };

    ResidualStressEstimator estimator {machining_scans, tool_paths};

    return EXIT_SUCCESS;
}
