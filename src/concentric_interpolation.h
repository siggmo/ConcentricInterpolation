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
 *  Authors: Oliver Kunc and Felix Fritzen
 *  Title  : Generation of energy-minimizing point sets on spheres and their
 *           application in mesh-free interpolation and differentiation
 *  Journal: Advances in Computational Mathematics 45(5-6), pp. 3021-3056
 *  Year   : 2019
 *  URL    : https://doi.org/10.1007/s10444-019-09726-5
 *
 *  The latest version of this software can be obtained through
 *  https://github.com/EMMA-Group/ConcentricInterpolation
 *
 */

#ifndef __CONCENTRIC_INTERPOLATION_H__
#define __CONCENTRIC_INTERPOLATION_H__

#include "util.h"
#include "data.h"
#include "tangential_interpolation.h"
#include "radial_interpolation.h"

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

    int                     D_val;                  //!< dimension of the data, i.e. number of scalar values to be interpolated at once
    bool                    initialized;
    bool                    printed_info_error;     //!< when calling Error(), a message is printed once-in-a-runtime

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

    /** \brief Efficient single evaluation of ConcentricInterpolation at one
     * point, given by direction and radius
     *
     * \attention The user must assure sufficient memory allocation for the arguments
     *
     */
    inline void Evaluate(
        const double *  Direction,
        const double    Radius,
        double *        ResultValues
         )
    {
        RadialInterpolant.Interpolate( Radius, RI_result );
        TangentialInterpolant.KernelVector( 1, Direction, TI_result );
        MatVecMul(RI_result, TI_result, ResultValues, TI_result_length, D_val, true);
    }

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
     * The second argument specifies how the error is computed. If \p a_error_type is 1, then the mean relative error is computed:
     * \f[ {\rm error} = \left( \sum_{i=1}^{N^{\rm eval}_{\rm points}} \frac{\left| \underline{f}^{\rm error}(x_i)-\widetilde{\underline{f}}^{\rm error}(x_i)\right|}{\left|\underline{f}^{\rm error}(x_i)\right|}\right)/N^{\rm eval}_{\rm points} \f] If \p a_error_type is 2, then the rooted mean square is computed:
     * \f[ {\rm error} = \sqrt{\sum_{i=1}^{N^{\rm eval}_{\rm points}} \left| \underline{f}^{\rm error}(x_i) - \widetilde{\underline{f}}^{\rm error}(x_i) \right|^2 \big/ N^{\rm eval}_{\rm points} }\f]
     * Here, \f$N^{\rm eval}_{\rm points}=N^{\rm eval}_{\rm dir}N^{\rm eval}_{\rm mag}\f$ is the number of evaluation points. The vectors \f$ \underline{f}^{\rm error},\widetilde{\underline{f}}^{\rm error}\in\mathbf{R}^{D^{\rm error}_{\rm values}}\f$ either equal the original vectors \f$ \underline{f},\widetilde{\underline{f}}\in\mathbf{R}^{D_{\rm values}\f$ (when \f$ D^{\rm error}_{\rm values} = D_{\rm values} \f$), or consist of a subset of the components (when \f$ D^{\rm error}_{\rm values} < D_{\rm values} \f$). In the latter case
     * \f$ D^{\rm error}_{\rm values} = \f$ <tt> d_val_end - d_val_start + 1</tt>. The variables on the right-hand side are defined by the arguments \p a_d_val_start and \p a_d_val_end, which default to 0 and \p D_val -1, respectively.
     *
     * > Example: in the example file interpolation_demo_largestrain.cxx, the vector \f$ [W, \underline{S}, \underline{C}] \in \mathbf{R}^{1+6+21} \f$ is interpolated. The stress (\f$\underline{S}\f$) error can be computed by choosing \p a_d_val_start=1 and \p a_d_val_end=6.
     *
     * \note This is where you can implement an error function that suits your purpose.
     *
     * The last argument defines whether the error is computed on the whole set of validation coordinates (default, i.e. when \p a_OnlyRadius = -1) or only on the sphere with index 0 \f$\leq\f$ \p a_OnlyRadius \f$\leq N^{\rm eval}_{\rm dir}-1\f$ (zero-based).
     *
     * \attention if <tt>a_OnlyRadius==-1</tt> (default value), the error may be biased because ConcentricData includes the origin for each direction, and the output is always zero at the origin.
     */
    double Error(
        const ConcentricData&   a_DataValidation,   //!< [in] ConcentricInterpolation will be evaluated at the contained coordinates and the results will be compared against the contained values
        const int               a_error_type = 2,   //!< [in] if 1: mean relative error; if 2: rooted mean square of absolute difference
        const int               a_d_val_start = -1, //!< [in] starting index of \f$ \underline{f}^{\rm error} \f$ within \f$ \underline{f} \f$, zero-based. See also main description of Error()
        const int               a_d_val_end = -1,   //!< [in] ending index of \f$ \underline{f}^{\rm error} \f$ within \f$ \underline{f} \f$, zero-based. See also main description of Error()
        const int               a_OnlyRadius = -1   //!< [in, optional] zero-based index of the only radius that should be considered (such that the error is computed on a sphere). default value -1 means that all radii will be considered, except the radius zero (at which ConcentricData must be zero in any case).
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
        const int       a_error_type,           //!< [in, optional] see Error()
        const int       a_d_val_start,          //!< [in, optional] see Error()
        const int       a_d_val_end,            //!< [in, optional] see Error()
        const int       a_OnlyRadius = -1       //!< [in, optional] see Error()
               );

    void                    PrintInfo(FILE*F=stdout);//!< print information about this object, e.g. value of kernel parameter, number of kernels, number of radial support points, etc.

};

#endif
