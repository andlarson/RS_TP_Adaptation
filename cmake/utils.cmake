# When using GCC version 8 (and other older compilers), the default C++ standard
#     library doesn't support the filesystem library . This is a known problem,
#     and this solution was taken from:
#     https://discourse.cmake.org/t/correct-way-to-link-std-filesystem-with-gcc-8/4121/6
#     There is a discussion on the CMake forum about this issue here:
#     https://gitlab.kitware.com/cmake/cmake/-/issues/17834   
# This approach is slightly flawed because it only fixes the problem for GCC.
#     This issue also exists on old versions of Clang. For now, this is good
#     enough.
function(GCC8_filesystem_stl_fix)
    link_libraries("$<$<AND:$<CXX_COMPILER_ID:GNU>,$<VERSION_LESS:$<CXX_COMPILER_VERSION>,9.0>>:-lstdc++fs>")
endfunction()
