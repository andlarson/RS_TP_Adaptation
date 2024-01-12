// Included in Abaqus install.
// Includes macros to call C/C++ from Fortran.
#include <omi_for_c.h>

// For some reason the system is unable to find this header. This header
//    will be helpful if we ever want to call Fortran subroutines defined
//    by Abaqus from our C++ code. This header really has nothing to do with
//    GCC, it is included via the Intel Fortran compiler. One potential
//    approach to get this working is to overwrite the compilation options
//    in the Abaqus environment files via a new local environment file.
// See the Intel Fortran compiler documentation for more information.
// Included via Intel Fortran compiler.
// Includes utilities for calling Fortran from C/C++.
// #include <ISO_Fortran_binding.h>

#include <cstddef>
#include <cstdlib>
#include <cstring>


// GCC is the compiler and C++ 14 is being used.


// Get rid of this eventually...
using namespace std;


// Utility subroutines written in Fortran.
// extern "C" void stdb_abqerr_(int *, char *, int *, double *, char **);


// This function is called at particular material points.
// Note that the argument types are not specified in the Abaqus documentation - 
//    I inferrred what they should be.
extern "C" void FOR_NAME(sigini, SIGINI) (double *sigma, double *coords, int *ntens,
                                          int *ncrds, int *noel, int *npt,
                                          int *layer, int *kspt, int *lrebar,
                                          char **names)
{
    // sigma[0] = coords[0] * coords[0];
    sigma[0] = 10^8;

    // Calling a Fortran subroutine.
    // Doing things this way doesn't quite work. It appears as though the lop
    //    argument is properly interpreted but the my_msg argument is not. This
    //    is evident by looking at the .msg file. There are message headings, but
    //    the messages don't go through properly.
    /* 
    int *ints_to_output = NULL;
    double *reals_to_output = NULL;
    char **strings_to_output = NULL; 
    int lop = 1;
    char *my_msg = "HELLO WORLD!!!";

    stdb_abqerr_(&lop, my_msg, ints_to_output, reals_to_output, strings_to_output);
    */

     

}


