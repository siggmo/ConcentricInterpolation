#ifndef __DATA_UTIL_H__
#define __DATA_UTIL_H__

/*
 *  ConcentricInterpolation
 *  Copyright (C) 2018  Felix Fritzen    ( felix.fritzen@mechbau.uni-stuttgart.de )
 *                      and Oliver Kunc  ( oliver.kunc@mechbau.uni-stuttgart.de )
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <https://www.gnu.org/licenses/>.
 *  
 *  
 *  For details or if you like this software please refer to LITERATURE which
 *  contains also BIBTEX information.
 *  
 *  The latest version of this software can be obtained through https://github.com/EMMA-Group/ConcentricInterpolation
 *  
 *  
 */

#include "concentric_interpolation.h"
#include "data.h"
#include "distance_functions.h"

/*! \brief 
 * Interpolate data on a set of coordinates [OPTIONAL: compute the gradient and the Hessian, too] */
void InterpolateOnSet(  ConcentricInterpolation * interpolation, /*!< [in] the interpolation scheme that will be evaluated */
                        Data * data_interpolation,              /*!< [in,out] a Data object containing the point coordinates for evaluation of the interpolation scheme, and the memory for the results */
                        const bool evaluate_gradients= true,    /*!< [in] compute gradients from interpolated values? (implies function interpolation) */
                        const bool evaluate_hessians = true     /*!< [in] compute hessians? (implies gradient) */
                     );

/*! \brief Class for evaluation of the distance function
 * 
 *  TODO
 * 
 * \see \GammaOptFunctional
 * \see \OptimizeGamma
 * \see \DistanceFunction
 * \see \ConcentricInterpolation
 */
class TestOnSet{
public:
    TestOnSet();
    TestOnSet(  ConcentricInterpolation * a_interpolation,
                Data *                   a_data_reference,
                DistanceFunction_local * a_distfunc_local,
                DistanceFunction_global* a_distfunc_global,
                double (*a_objective_function)(const double, const double, const double)
             ); //!< create TestOnSet object and run the test for interpolation scheme \c a_interpolation and coordinates and reference values from \c a_data_reference with the local and global distance functions \c a_distfunc_local and \c a_distfunc_global specified and via the objective function \c a_objective_function
    ~TestOnSet();
    double run(   const bool evaluate_gradients = true, /*!< [in] compute gradients from interpolated values? */
                const bool evaluate_hessians  = true    /*!< [in] compute hessians? */
                                 );
    ConcentricInterpolation * interpolation;     /*!< pointer to the interpolation scheme to be tested */
    Data data_interpolation;                    /*!< object containing the results of the interpolation on the reference coordinates */
    Data * data_reference;                      /*!< pointer to the reference data */
    const int D;                                /*!< dimension of the point coordinates */
    const int P;                                /*!< number of points in the test set */
    double * dist_values_local;                 /*!< point-wise distance (error) of interpolated values from reference values */
    double * dist_gradients_local;              /*!< point-wise distance (error) of approximate gradients from reference gradients */
    double * dist_hessians_local;               /*!< point-wise distance (error) of approximate Hessians from reference Hessians */
    
    DistanceFunction_local * distfunc_local; /*!< local distance function */
    DistanceFunction_global* distfunc_global; /*!< global distance function */
    double (*objective_function)(const double, const double, const double); /*!< global objective function */
};

    


/*! \brief Optimize gamma with respect to target data TODO incomplete!!!!
    * 
    * This function searches \c gamma within the given limits such that the specified functional
    * is minimized. The evaluation points are the given target coordinates, i.e.
    * CreateTargetData must have been called before.
    * 
    * The cost functional GammaOptFunctional is minimized w.r.t. gamma. For this, the given 
    * interval of gammas is bisected until either
    *    1) the given number of iterations is reached or
    *    2) gamma doesn't change significantly or
    *    3) the quality of the interpolation doesn't improve significantly anymore,
    * whatever occurs first.
    * 
    * It is assumed that the quality is a convex function within the given gamma interval. This is
    * backed by experience for sufficiently short intervals and if \c gamma_min is not too small.
    * 
    * Typically, we have found that \c gamma_min should be above 0.5. The value of a reasonable
    * \c gamma_max depends on the data, its dimension and on the number of training directions. It
    * might be necessary to call this function multiple times with different \c gamma limits.
    * 
    * \see CreateTargetData
    * \see GammaOptFunctional
    * */
double  OptimizeGamma(  ConcentricInterpolation * a_interpolation,   /*!< [in] the interpolation scheme */
                        Data * data_reference,                      /*!< [in] reference data on which the objective function will be evaluated */
                        DistanceFunction_local * a_distfunc_local,  /*!< [in] pointer to local distance function */
                        DistanceFunction_global* a_distfunc_global, /*!< [in] pionter to global distance function */
                        double (*objective_function)(const double, const double, const double), /*!< [in] pointer to objective function, i.e. the returned gamma must minimize this */
                        const double gamma_min,                     /*!< [in] minimum gamma value*/
                        const double gamma_max,                     /*!< [in] maximum gamma value*/
                        const int num_regular=10,                   /*!< [in] the interval [gamma_min, gamma_max] will contain num_regular+2 equi-spaced points used for global optimization.*/
                        const int num_bisec=4,                      /*!< [in] maximum number of bisectional iterations*/
                        const double bisection_factor=1,            /*!< [in] bisection factor, from open interval (0,2). this moves the next evaluation either towards the current best value (<1) or further away (>1)*/
                        const bool do_gradients=true,               /*!< [in] compute gradients, their errors, and optimize gamma wrt. this*/
                        const bool do_hessians=true                 /*!< [in] compute hessians, their errors, and optimize gamma wrt. this*/
);




void    BuildInterpolationFromTrainingData(
    ConcentricInterpolation & interpolation,     /*!< [in] the interpolation scheme */
    DataTraining * data_training                /*!< [in] the training data on which the interpolation scheme should be built */
);

void    WriteDataGnuplot(
    Data & data,                            //!< [in] reference to the data that is to be written to a text file
    char * output_filename,                 //!< [in] the name of the output text file
    const bool write_coordinates = true,    //!< [in] if true, the coordinates will be written
    const bool write_values      = true,    //!< [in] if true, the function values will be written
    const bool write_gradients   = true,    //!< [in] if true, the gradients will be written
    const bool write_hessians    = true     //!< [in] if true, the Hessians will be written. ATTENTION: this is the full matrix, i.e. D*D entries will be written per point (i.e. per row)
);

#endif // __DATA_UTIL_H__
