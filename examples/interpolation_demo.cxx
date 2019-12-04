#include "concentric_interpolation.h"

/**! \brief Concentric Interpolation demo
 *
 * Optimizes the kernel parameter on a subset of the validation data,
 * and computes the \e absolute error on the whole validation data.
 */

int main()
{
    // set up support data
    char FilenameSupportDirections[]    = "../data/directions/D8_N166_from_N83_sym1_s-2.txt";
    char FilenameSupportRadii[]         = "../data/radii/radii_5.txt";
    char FilenameSupportValues[]        = "../data/values/D8_N166from83_s-2_R5.txt";
    ConcentricData SupportData( FilenameSupportDirections,
                                FilenameSupportRadii,
                                FilenameSupportValues );

    // set up the Concentric Interpolation scheme
    ConcentricInterpolation ConcentricInterpolant(SupportData);

    // set up validation data
    char FilenameValidationDirections[] = "../data/directions/D8_N512_from_N256_sym1_s-2.txt";
    char FilenameValidationRadii[]      = "../data/radii/radii_5.txt";
    char FilenameValidationValues[]     = "../data/values/D8_N512from256_s-2_R5.txt";
    ConcentricData ValidationData( FilenameValidationDirections,
                                   FilenameValidationRadii,
                                   FilenameValidationValues );


    // find the optimal gamma value for this validation set
    const double gamma_min = 0.5;
    const double gamma_max = 2.0;
    const int   N_gamma_regular   = 10;
    const int   N_gamma_bisection = 5;
    ConcentricInterpolant.OptimizeGamma(    gamma_min,
                                            gamma_max,
                                            N_gamma_regular,
                                            N_gamma_bisection,
                                            ValidationData,
                                            4 // ATTENTION optimize gamma only on the sphere with the radius that is indicated by this zero-based number
                                       );

    // compute the error on the validation set, using the above found gamma
    printf("\n\nerror on whole validation data = %e\n", ConcentricInterpolant.Error(ValidationData));

    return 0;
}
