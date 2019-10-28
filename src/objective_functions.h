#ifndef __OBJECTIVE_FUNCTIONS_H__
#define __OBJECTIVE_FUNCTIONS_H__

/*
 *  COPYRIGHT NOTES
 *
 *  ConcentricInterpolation
 *  Copyright (C) 2018  Felix Fritzen    ( fritzen@mechbau.uni-stuttgart.de )
 *                      and Oliver Kunc  ( kunc@mechbau.uni-stuttgart.de )
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
 *     JOURNAL NAME, Number/Volume, p. XX-YY, 2019
 *     DOI   ...
 *     URL   dx.doi.org/...
 *
 *  The latest version of this software can be obtained through
 *  https://github.com/EMMA-Group/ConcentricInterpolation
 *
 */

//! NOTE: add your very own objective function, e.g. incorporating physics-motivated weights

/** Objective function that treats functions value, gradient and Hessian equally
 *  This is similar to the Sobolev H2 norm if
 *      Dist_local_AbsDiff is used as local distance function
 *      Dist_RMS_global is used as global distance
 *      the evaluation points allow for equal weight numerical quadrature
 *  In this case, the objective function value is (\f$ g(x_p)\f$ denotes the gradient of \f$ f \f$ at \f$ x_p \f$, and \f$ H(x_p)\f$ denotes the Hessian of \f$ f \f$ at \f$ x_p \f$)
 * \f[\left(\frac{\sum_{p=1}^P|f(x_p)-\widetilde{f}(x_p)|}{P}\right)^{1/2} + \left(\frac{\sum_{p=1}^P|| g(x_p)-\widetilde{g}(x_p)||}{P}\right)^{1/2} + \left(\frac{\sum_{p=1}^P|| H(x_p)-\widetilde{H}(x_p)rVert}{P}\right)^{1/2}\f]
 *  which is not quite the actual Sobolev H2 norm, but close.
 */
inline double H2_like_objective(const double dist_global_values, const double dist_global_gradients, const double dist_global_hessians)
{
    return dist_global_values+dist_global_gradients+dist_global_hessians;
}

/** Objective function that considers the gradient and the function value
 *  This is similar to the Sobolev H1 norm if
 *      Dist_local_AbsDiff is used as local distance function
 *      Dist_RMS_global is used as global distance
 *      the evaluation points allow for equal weight numerical quadrature
 *  In this case, the objective function value is (\f$ g(x_p)\f$ denotes the gradient of \f$ f \f$ at \f$ x_p \f$)
 * \f[\left(\frac{\sum_{p=1}^P|f(x_p)-\widetilde{f}(x_p)|}{P}\right)^{1/2} + \left(\frac{\sum_{p=1}^P|| g(x_p)-\widetilde{g}(x_p)||}{P}\right)^{1/2}\f]
 *  which is not quite the actual Sobolev H1 norm, but close.
 */
inline double H1_like_objective(const double dist_global_values, const double dist_global_gradients, const double dist_global_hessians=0)
{
    return dist_global_values + dist_global_gradients;
}

/** Objective function that considers only the global distance of the values, i.e. ignoring gradients and Hessians.
 *  This is the L2 norm if
 *      Dist_local_AbsDiff is used as local distance function
 *      Dist_RMS_global is used as global distance
 *      the evaluation points allow for equal weight numerical quadrature
 *  In this case, the objective function value is
 * \f[\left(\frac{\sum_{p=1}^P|f(x_p)-\widetilde{f}(x_p)|}{P}\right)^{1/2}\f]
 */
inline double L2_like_objective(const double dist_global_values, const double dist_global_gradients=0, const double dist_global_hessians=0)
{
    return dist_global_values;
}

#endif // define __OBJECTIVE_FUNCTIONS_H__
