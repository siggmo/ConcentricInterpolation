#include "concentric_interpolation.h"

/** \brief Concentric Interpolation demo
 *
 * This is an implementation of the example from Section 7.2 of the original
 * article [https://doi.org/10.1007/s10444-019-09726-5]
 *
 * Optimizes the kernel parameter on a subset of the validation data, using
 * relative errors on a spherical subset of the validation data.
 *
 * Then, computes the absolute error on the whole validation data.
 *
 */

int main()
{
    // set up support data
    char FilenameSupportDirections[]    = "../data/directions/D8_N166_from_N83_sym1_s-2.txt";
    char FilenameSupportRadii[]         = "../data/radii/radii_4.txt";
    char FilenameSupportValues[]        = "../data/values/D8_N166from83_s-2_R4.txt";
    ConcentricData SupportData( FilenameSupportDirections,
                                FilenameSupportRadii,
                                FilenameSupportValues );

    // set up the Concentric Interpolation scheme
    ConcentricInterpolation ConcentricInterpolant(SupportData);

    // set up validation data
    char FilenameValidationDirections[] = "../data/directions/D8_N512_from_N256_sym1_s-2.txt";
    char FilenameValidationRadii[]      = "../data/radii/radii_4.txt";
    char FilenameValidationValues[]     = "../data/values/D8_N512from256_s-2_R4.txt";
    ConcentricData ValidationData( FilenameValidationDirections,
                                   FilenameValidationRadii,
                                   FilenameValidationValues );


    // find the optimal gamma value for this validation set
    const double gamma_min = 0.1;
    const double gamma_max = 2.0;
    const int   N_gamma_regular   = 10;
    const int   N_gamma_bisection = 5;
    ConcentricInterpolant.OptimizeGamma(    gamma_min,
                                            gamma_max,
                                            N_gamma_regular,
                                            N_gamma_bisection,
                                            ValidationData,
                                            1,                  // use relative error, see ConcentricInterpolation::Error()
                                            0,                  // first component of values
                                            0,                  // last component of values
                                            2                   // only compute error on radius with this index (zero-based), see ConcentricInterpolation::Error()
                                       );

    // compute the error on the validation set, using the above found gamma
    printf("\nRMS absolute error on whole validation data = %e\n\n",
           ConcentricInterpolant.Error(ValidationData));

    return 0;
}
