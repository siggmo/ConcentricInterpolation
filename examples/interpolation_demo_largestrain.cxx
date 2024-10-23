#include "concentric_interpolation.h"

/** \brief Example of interpolating a finite strain hyperelastic material law.
 * This complements the research article
 * AUTHORS: Oliver Kunc and Felix Fritzen
 * TITLE  :
 * JOURNAL:
 * YEAR   : 2020
 * DOI    :
 *
 * The Extended Tube Model [Kaliske, Heinrich, https://doi.org/10.5254/1.3538822]
 * is interpolated. The following material parameters were chosen (according to
 * [Marckmann, Verron, https://doi.org/10.5254/1.3547969]) for the isochoric part:
 *
 *  - \f$ G_c = 0.202 \f$ MPa
 *  - \f$ G_e = 0.153 \f$ MPa
 *  - \f$ \beta = 0.178 \f$ MPa
 *  - \f$ \delta = 0.0856 \f$ MPa
 *
 * For the volumetric part, the energy \f$ U(J) = 0.25 \cdot K \cdot ( (J-1)^2 + {\rm log}(J)^2 ) \f$
 * from [Doll, Schweizerhof, https://doi.org/10.1115/1.321146] is used with the
 * bulk modulus
 *  - \f$ K = 1 \f$ MPa.
 *
 * For the setup of the Concentric Interpolation (CI), the following parameters
 * are chosen:
 *  - \f$ N^{\rm supp}_{\rm dir}=200 \f$ (directions from file \e data/directions/D6_N200_sym0_s-2.txt )
 *  - \f$ N^{\rm supp}_{\rm mag}=10\f$ (\f$ t^{(i)}=i/10, i=1,\ldots,10 \f$)
 *  - \f$ J^* = 1.001 \f$
 *
 * The CI interpolates 28 scalar values on the space of Hencky strains, \f$ \mathbf{E} \f$:
 * - (1) hyperelastic energy density \f$ W_{\rm C}({\rm exp}(2\mathbf{E}) \f$
 * - (6) second Piola-Kirchhoff stress \f$ \mathbf{S} \f$
 * - (21) tangent modulus of the second P.K. stress \f$ \boldmath{C} \f$
 *
 * For validation purposes, exact function data is available at the following
 * Concentric Sampling sites:
 *  - \f$ N^{\rm eval}_{\rm dir}=100 \f$ (directions from file \e data/directions/D6_N100_sym0_s-2.txt )
 *  - \f$ N^{\rm eval}_{\rm mag}=15\f$ (\f$ t^{(i)}=i/10, i=1,\ldots,15 \f$)
 *  - \f$ J^* = 1.001 \f$
 *
 * First, an optimization of the kernel parameter \f$ \gamma \f$ is performed on
 * a subset of the validation data, namely on all validation data with a certain
 * magnitude, employing a relative error measure on the stresses alone.
 *
 * Then, an RMS of the absolute stress error is computed on all validation data.
 * This is the final output.
 */

int main(int argc, char* argv[])
{
    // set up support data
    char FilenameSupportDirections[]    = "../data/directions/D6_N200_sym0_s-2.txt";
    char FilenameSupportRadii[]         = "../data/radii/radii_10.txt";
    char FilenameSupportValues[]        = "../data/values/CI_data_ETM_Nsuppdir200_Nsuppmag10_Jmax1.001.txt";
    ConcentricData SupportData( FilenameSupportDirections,
                                FilenameSupportRadii,
                                FilenameSupportValues,
                                28 /* interpolate 28 scalar values separately & simultaneously: energy (1), Stress(6), Stiffness(21) */
                              );

    // set up the Concentric Interpolation scheme
    ConcentricInterpolation ConcentricInterpolant(SupportData);

    // set up validation data
    char FilenameValidationDirections[] = "../data/directions/D6_N100_sym0_s-2.txt";
    char FilenameValidationRadii[]      = "../data/radii/radii_15.txt";
    char FilenameValidationValues[]     = "../data/values/CI_data_ETM_Nsuppdir100_Nsuppmag15_Jmax1.001.txt";
    ConcentricData ValidationData( FilenameValidationDirections,
                                   FilenameValidationRadii,
                                   FilenameValidationValues,
                                   28 /* interpolate 28 scalar values separately & simultaneously: energy (1), Stress(6), Stiffness(21) */
                                 );

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
                                            1,                  // use relative error, see ConcentricInterpolation::Error()
                                            1,                  // first component of values for stress computation, see ConcentricInterpolation::Error()
                                            6,                  // last component of values for stress computation, see ConcentricInterpolation::Error()
                                            11                  // only compute error on radius with this index (zero-based), see ConcentricInterpolation::Error()
                                       );

    // compute the absolute RMS error on the whole validation set, using the above found gamma
    printf("\nRMS absolute stress error on all validation data = %e\n\n",
           ConcentricInterpolant.Error(ValidationData, 2, 1, 6));

    return 0;
}
