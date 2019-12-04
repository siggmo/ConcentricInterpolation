#ifndef __CONCENTRIC_INTERPOLATION_H__
#define __CONCENTRIC_INTERPOLATION_H__

#include "util.h"
#include "data.h"
#include "tangential_interpolation.h"
#include "radial_interpolation.h"
/*
 *  COPYRIGHT NOTES
 *
 *  ConcentricInterpolation
 * Copyright (C) 2018-2019 by Felix Fritzen (fritzen@mechbau.uni-stuttgart.de)
 *                         and Oliver Kunc (kunc@mechbau.uni-stuttgart.de
 * All rights reserved.
 *
 * This source code is licensed under the BSD 3-Clause License found in the
 * LICENSE file in the root directory of this source tree.
 *
 *  This software package is related to the research article
 *
 *     Oliver Kunc and Felix Fritzen: 'Generation of energy-minimizing point
 *                                     sets on spheres and their application in
 *                                     mesh-free interpolation and
 *                                     differentiation'
 *     Advances in Computational Mathematics, Number/Volume, p. XX-YY, 2019
 *     DOI   10.1007/s10444-019-09726-5
 *     URL   dx.doi.org/10.1007/s10444-019-09726-5
 *
 *  The latest version of this software can be obtained through
 *  https://github.com/EMMA-Group/ConcentricInterpolation
 *
 */

/** \brief Main structure for Concentric Interpolation
 *
 * Create an instance of this class and call Setup() to get started.
 *
 * This is basically a wrapper for instances of RadialInterpolation
 * and TangentialInterpolation providing convenient access to and
 * post-processing routines for their results.
 *
 * You may evaluate this concentric interpolant via Evaluate(), or
 * compute the error w.r.t. some data via Error().
 *
 * If you don't know which kernel parameter \f$ \gamma \f$ to use,
 * see the article for some experience or call OptimizeGamma().
 *
 * Long story short: if \f$ N_{\rm dir}^{\rm supp} \geq 200 \f$,
 * try \f$ \gamma = 2.0 \f$, else try \f$ \gamma = 1.0 \f$.
 *
 * \see RadialInterpolation
 * \see TangentialInterpolation
 */

class ConcentricInterpolation
{
private:
    TangentialInterpolation TangentialInterpolant;
    RadialInterpolation     RadialInterpolant;
    Interpolant             LinearInterpolant;      //!< adjust this to your needs, e.g. QuadraticInterpolant or CubicInterpolant

    double *                RI_result;              //!< this is \f$ \underline{\mathcal{R}}(x) \underline{\underline{K}}^{-1} \f$ of equation (CI) of the article, length \p RI_result_length
    double *                TI_result;              //!< this is \f$ \underline{k}(\mathbf{n}) \f$ of equation (CI) of the article, length \p TI_result_length

    size_t                  RI_result_length;       //!< length of \p RI_result = \f$ N_{\rm dir}^{\rm eval} \cdot D_{\rm val}\f$
    size_t                  TI_result_length;       //!< length of \p TI_result = \f$ N_{\rm dir}^{\rm eval} \f$

    bool                    initialized;
    bool                    print_info;             //!< when calling Error(), a message is printed once-in-a-runtime

public:
    ConcentricInterpolation();                              //!< Default constructor initializing empty objects
    ConcentricInterpolation(ConcentricData& a_SetupData);   //!< Constructor calling Setup()

    ~ConcentricInterpolation();

    /** \brief Get started with this function.
    *
    * Then, see Evaluate() or Error(), for instance.
    */
    void Setup(
        const ConcentricData& a_SetupData,          //!< [in] concentric support data
        const double a_gamma = 0.,                  //!< [in,optional] kernel width parameter \f$ \gamma \f$, \see OptimizeGamma()
        const bool a_sym = false                    //!< [in,optional] symmetry flag. if true, use symmetrized kernel functions for the interpolation
         );

    /** \brief Evaluates ConcentricInterpolation on the GeneralCoordinates provided by the input
     * argument and stores the resulting values to the storage therein.
     *
     * \attention This is untested yet.
     *
     * \attention The interpolation's values will overwrite the values that are possibly stored in the
     * input argument.
     *
     * \todo test this
     * \todo add overloaded function Evaluate(const ConcentricData&)
     */
    void Evaluate(
        GeneralData&    a_DataEvaluation
         );

    /** \brief Evaluates ConcentricInterpolation on concentric coordinates and computes error.
     *
     * Computes the error between the results and the values stored within the input argument, and
     * returns a scalar error value.
     *
     * \note This is where you can implement an error function that suits your purpose.
     *
     * \attention if <tt>a_OnlyRadius==-1</tt> (default value), the error may be biased because ConcentricData includes the origin for each direction
     *
     * \todo make user-defined error functionals more accessible
     */
    double Error(
        const ConcentricData&   a_DataValidation,   //!< [in] ConcentricInterpolation will be evaluated at the contained coordinates and the results will be compared against the contained values
        const int               a_OnlyRadius = -1   //!< [in, optional] zero-based index of the only radius that should be considered (such that the error is computed on a sphere). default value -1 means that all radii will be considered
                );

    /** \brief Changes the kernel parameter \f$ \gamma \f$ and processes the
     * interpolation data accordingly.
     *
     * Changes the parameter via TangentialInterpolation::SetGamma() and adjusts
     * the data via RadialInterpolation::MultiplyKernelMatrixAndReorder().
     *
     * User has to decide whether kernel matrix should be re-computed
     * immediately or not. ConcentricInterpolation can only be used when kernel
     * matrix is re-computed after a change of \f$ \gamma \f$, but the user
     * might prefer to perform the costly re-computation after at a later point
     * (e.g. after calling TangentialInterpolation::SetLambda(), ...)
     */
    void SetGamma(  const double a_gamma = 1.0,                 //!< [in] see TangentialInterpolation::SetGamma()
                    const bool a_recompute_kernel_matrix = true,//!< [in] see TangentialInterpolation::SetGamma()
                    const bool a_quiet = false                  //!< [in] see TangentialInterpolation::SetGamma()
    );

    /** \brief Optimizes the Gaussian kernel parameter \f$ \gamma \f$ such that the error implemented
     * in Error() is minimized on the input data set.
     *
     * First, \f$ \gamma \f$ is varied from \f$ \gamma_{\rm min} \f$ to \f$ \gamma_{\rm max} \f$ by
     * \f$ N_{\gamma,\rm regular} \f$ equidistant steps. Then, a bisection algorithm with
     * \f$ N_{\gamma,\rm bisection} \f$ bisection steps is performed on the interval bounded by the
     * best value of \f$ \gamma \f$ from the first step and its best neighbor. The final result is
     * applied via TangentialInterpolation::SetGamma().
     */
    void OptimizeGamma(
        const double    a_gamma_min,            //!< [in] \f$ \gamma_{\rm min} \f$
        const double    a_gamma_max,            //!< [in] \f$ \gamma_{\rm max} \f$
        const int       a_N_gamma_regular,      //!< [in] \f$ N_{\gamma,\rm regular} \f$
        const int       a_N_gamma_bisection,    //!< [in] \f$ N_{\gamma,\rm bisection} \f$
        const ConcentricData&  a_DataValidation,//!< [in] data set with respect to which \f$ \gamma \f$ is optimized
        const int       a_OnlyRadius = -1       //!< [in, optional] zero-based index of the only radius that should be considered. if negative, then all radii will be considered. \see Error()
               );

};

#endif
