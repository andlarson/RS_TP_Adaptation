// Required for Fortran/C++ interoperability.
#include <omi_for_c.h>

#include <cstddef>
#include <cstdlib>
#include <cstring>

// For DEBUG
#include <fstream>
#include <iostream>


// Get rid of this eventually...
using namespace std;


// Utility subroutines written in Fortran.
extern "C" void stdb_abqerr_(int *, char *, int *, double *, char **);

static int testing = 1;

// This function is called at particular material points.
// Note that the argument types are not specified in the Abaqus documentation - 
//    I inferrred what they should be.
extern "C" void FOR_NAME(sigini, SIGINI) (double *sigma, double *coords, int *ntens,
                                          int *ncrds, int *noel, int *npt,
                                          int *layer, int *kspt, int *lrebar,
                                          char **names)
{
    sigma[0] = 10;

    if (testing == 1) {
        ofstream myfile;
        myfile.open("/home/andlars/Downloads/abaqus_script_logging.txt");
        myfile << "Got here!";
        myfile.close();
        testing = 0;
    }

    int *ints_to_output = NULL;
    double *reals_to_output = NULL;
    char **strings_to_output = (char **) 5;

    int *lop = (int *) malloc(sizeof(int));
    *lop = 1;

    char *message = (char *) malloc(15 * sizeof(char));
    char *my_msg = "HELLO WORLD!!!";
    memcpy(message, my_msg, 15);

    stdb_abqerr_(lop, message, ints_to_output, reals_to_output, strings_to_output);

    free(message);
    free(lop);
}


