// Required for Fortran/C++ interoperability.
#include <omi_for_c.h>

#include <cstddef>



extern "C" void FOR_NAME(std_abqerr, STD_ABGERR);

// This function is called at particular material points.
// Note that the argument types are not specified in the Abaqus documentation - 
//    I inferrred what they should be.
extern "C" void FOR_NAME(sigini, SIGINI) (double *sigma, double *coords, int *ntens,
                                          int *ncrds, int *noel, int *npt,
                                          int *layer, int *kspt, int *lrebar,
                                          char **names)
{
    sigma[0] = 10;

    int *ints_to_output = NULL;
    double *reals_to_output = NULL;
    char **strings_to_output = NULL;
    std_abqerr(1, "test", ints_to_output, reals_to_output, strings_to_output);
}


