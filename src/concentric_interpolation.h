#ifndef __CONCENTRIC_INTERPOLATION_H__

#define __CONCENTRIC_INTERPOLATION_H__

/*
 *  COPYRIGHT NOTES
 * 
 *  ConcentricInterpolation
 *  Copyright (C) 2019  Felix Fritzen    ( fritzen@mechbau.uni-stuttgart.de )
 *                      and Oliver Kunc  ( kunc@mechbau.uni-stuttgart.de )
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
 *  (the full license is distributed together with the software in a file named
 *  LICENSE)
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



#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <iostream>
#include <lapacke.h>
#include <cubic_interpolant.h>
#include <quadratic_interpolant.h>
#include <util.h>


class ConcentricInterpolation;
//! The class ConcentricInterpolation serves as universal interpolator for smooth
//! functions with directional dependency in \f$R^D\f$. It relies on data provided in terms of
//! data points along directors, i.e. \f${f_{ij} = f(r_i, d_j)}\f$ is required with at least three data
//! points at different radii \f$r_i\f$ per normalized direction \f$d_j\f$.
//! 
//! Along the rays a piecewise cubic polynomial is used for the radial
//! interpolation. Between the directions a kernel method with Gaussian kernel that
//! takes as inputs the geodesic distances of the director is applied.
//! 
//! Member functions contain also differentiation operators up to second order.
//! Note that the input points are reproduced _exactly_ and a piecewise cubic polynomial interpolation
//! is returned along the rays provided as inputs. The number of points along different directions
//! is not required to be identical



class ConcentricInterpolation {
private:
    bool        sym;                //!< exploit point symmetry (only one ray of data points must be provided)
    CubicInterpolant    * S;        //!< array of piecewise cubic functions used for radial interpolation
    bool        init;               //!< initialization flag
    int         N_alloc;            //!< dimension of the pre-allocated memory
    int         D;                  //!< dimension of the inputs
    double      gamma;              //!< kernel width parameter of the Gaussian kernel function
    double      lambda;             //!< regression parameter: controls condition number of kernel matrix. usually not needed.
    
    double      radius;             //!< current r, i.e. || x ||
    // vectors:
    double      *x;                 //!< input vector x after normalization, i.e. d = x / r the direction of the input vector
    double      *v, *dv, *ddv,      //!< v, dv, ddv: function values and its first and second derivatives in radial direciton at radius r for all directions
                *Ki_zeta,           //!< Ki_zeta:   inverse kernel matrix times vector containing zeta vector
                *Ki_v,              //!< inverse kernel matrix times vector containing function values at radius r, interpolated by p.w. polynomials
                *Ki_dv,             //!< inverse kernel matrix times vector containing diff. of polynomial data
                *theta, *zeta, *xi, //!< theta_i:   d_i * d, xi: acos(theta), zeta_i: kernel function evaluated at xi_i
                *sin_xi,            //!< sine of xi_i (re-used several times)
                *zeta_tilde,        //!< for the symmetric case: kernel function at (PI-xi_i) (see paper)
                *zeta_star,         //!< for the symmetric case: kernel function at (PI-xi_i) (see paper)
                *dzeta,             //!< derivative of zeta_j w.r.t. xi_j         (sym:  dzeta_star )
                *ddzeta;            //!< second derivative of zeta_j w.r.t. xi_j  (sym: ddzeta_star )
    int         * w_i;              //!< integer working array for linear solver (IPIV in LAPACK)
    bool        * active;           //!< mark training directions which are active or not. active means not too close to the query direction, i.e. 1/sin(xi) is bound. inactive training directions will be assumed to have a minimum angular distance to the query direciton
    // matrix like 2D arrays (in terms of 1D row_major arrays). these are prefixed with "m_" to emphasize their multidimensional nature.
    double      *m_K,               //!< dense kernel matrix
                *m_Kf;              //!< the LDL factorization of the kernel matrix. used for solving systems of the form K * Ki_a = a for Ki_a
    double      *m_J0;              //!< Hessian for very small values of r

    static const double small; // = 1.e-12; //!< small number; used, e.g., in order to prevent division by zero errors
    static const double theta_max; // = 0.9985; //!< required in order to regularize the derivative of zeta!

    void zero_pointers();           //!< initialize all pointers to zero (before doing anything else)

    bool CheckRadius( const double * a_radii, const int a_R ) const; //!< returns true, if the first point is zero; if the data is not sorted then exit with an error
    
public:
    double      *m_X;               //!< training directions (pre-allocated with size N_alloc)
    int         N;                  //!< number of actually provided directions
    QuadraticInterpolant* Sq;       //!< array of piecewise quadratic functions used for radial interpolation
    
    ConcentricInterpolation( const bool a_sym = false /*!< [in] use symmetrization */ );
    /*!< Default constructur; requires call to Allocate before adding data
     * 
     * If \param a_sym is true, then the interpolation data is symmetrized,
     * i.e. I( X ) = I( -X ) is strictly enforced at (almost) no additional computational expense.
     * More precisely the dimension of the kernel matrix is not increased and the evaluation of the
     * interpolation (and of the gradients) involves only few additional operations.
     * 
     * \see Allocate
    */
    void SetSymmetric( const bool a_sym ) {
        sym = a_sym;
        if( init ) RecomputeKernelMatrix( );
    }
    
    ~ConcentricInterpolation();
    //!< clean up and destroy the object
    void    Free(); //!< free allocated memory \see ~ConcentricInterpolation()
    
    void SetJ0( const double * const J0 ); //!< set J0 (i.e. the Hessian for very "small" vectors of X)
    void GetJ0( double * J0 ) const; //!< return J0
    void AutodefineJ0_psi(); //!< set J0 based on the data for PSI
    void AutodefineJ0_dpsi(); //!< set J0 based on the data for dPSI/dr along training directions
    void AutodefineJ0( const double w_psi, const double w_dpsi ); //!< balanced computation of J0 based on the data

    void    Allocate( int a_N_alloc,        /*!< [in] maximum admissible number of directions used for input data */
                      const int a_D     /*!< [in] dimension of the input vector(s) */ );
    //!< \brief pre-allocate memory for (up to) \c a_N_alloc directions in R^\c a_D
    
    
    double  operator() ( const double* a_x /*!< [in] vector at which the interpolation is computed */ )
                                            //!< \brief Evaluate the interpolation function at the vector \c a_x
        { return Interpolate( a_x ); }
        
    //! \brief Copies the support points of the cubic polynomial for the *first* training directions to the arguments
    void    GetTrainingRadii( double * a_training_radii, int & a_num_training_radii ) const;
        
    /*! \brief Add interpoation data along a new direction
     * 
     * Adds a new interpolation direction with discrete data provided at different radii.
     * The data is interpolate using p.w. cubic polynomials in radial direction.
     * Note that the directions should **not** be identical or parallel (in the symmetric case).
     * otherwise, the kernel matrix becomes singular and the program will fail. 
     * */
    void    AddInterpolationDataDF(
        const double * a_radii, //!< [in] radii at which function values are provided
        const double * a_f,     //!< [in] discrete function values along the sampling direction
        const double * a_df,    //!< [in] discrete radial derivatives of the function along the sampling direction
        int a_n,                //!< [in] number of data points along the given direction (a_n >= 3 is required)
        const double * a_X      //!< [in] (unit) vector of dimension D; direction along which the data is provided
    );
    
    /*! \brief Return the kernel parameter
     * 
     * \see SetGamma
     * 
     * */
    inline double GetGamma( ) const
    {
        return gamma;
    }
    
    
    /*! \brief Set the kernel parameter
     * 
     * This also recomputes the kernel matrix and its inverse. If the kernel method has not been
     * initialized yet, it is also initialized.
     * 
     * \see InitializeKernelMethod
     * 
     * */
    void    SetGamma( const double a_gamma /*!< [in] new kernel parameter */ )
    {
        gamma = a_gamma;
        if(init)    { RecomputeKernelMatrix(); }
        else        { InitializeKernelMethod();}
    }

    /*! \brief change the regression parameter (adds some smoothness) */
    void SetLambda( const double a_lambda );
    
    /*! \brief Re-compute the kernel matrix and its inverse
     *
     * This does not re-allocate memory or change any member variables except m_K and m_Ki,
     * in contrast to InitializeKernelMethod.
     *
     * \see InitializeKernelMethod
     * 
     * */
    void    RecomputeKernelMatrix( ); // TODO OK: can/should this be private?
    
    /*! \brief (Re-)Compute the kernel matrix and its inverse
     * 
     * This member function evaluates the kernel matrix and computes its inverse.
     * A call to this function **must** be made before the interpolation function can be
     * called.
     * 
     * */
    void    InitializeKernelMethod( );
    
    
    /*! \brief Interpolate data at vector \c X
     * 
     * If the vector is zero, then the value of the first interpolation direction at radius zero is returned.
     * This function **must** be called before any calls to Gradient() of Hessian().
     * Before evaluating the interpolation function, the kernel parameter must be set via SetGamma()
     * and the kernel method must have been initialized via InitializeKernelMethod().
     * 
     * The computational complexity for one evaluation is proportional to the number of directions.
     * For each direction, a p.w. cubic polynomial interpolation is required. Further, an inner product is
     * be computed and the kernel function is evaluated. The interpolant is then given by an
     * appropriately weighted sum.
     * 
     * \see SetGamma
     * \see InitializeKernelMethod
     * \see Gradient
     * \see Hessian
     * */
    double  Interpolate( const double * a_x /*!< [in] vector \c X*/ );
    
    /*! Evaluation of piecewise cubic polynomials in radial direction
     * 
     * The derivatives are also computed (since computationally inexpensive
     * and required for \c Gradient/\c Hessian anyway).
     * 
     * \see Gradient
     * \see Hessian
     * 
     * */
    void    RadialInterpolation( const double a_radius /*!< [in] radius (i.e. norm of \c X) */ );
    
    /*! \brief Compute the gradient with respect to the current direction
     * 
     * The output is \f[res_i=\frac{\partial I}{\partial x_i}\f] 
     * 
     * \attention The computation requires a previous call to Interpolate.
     * 
     * \see Interpolate
     * \see Hessian
     * \see RadialInterpolation
     * */
    void    Gradient( double * grad /*!< [out] gradient of the interpolation function*/ );
    
    /*! \brief Compute the second gradient (i.e. the Hessian) of the interpolation function
     * 
     * The output consists of a one-dimensional array representing the dense matrix according to
     * \f[
     *      {m\_res}_{ (i-1)d + j } = J_{ij} = J_{ji} = \frac{ \partial^2 I}{ \partial x_i \, \partial x_j }
     * \f]
     * 
     * \attention The computation requires a previous call to Interpolate and Gradient.
     * 
     * \see Interpolate
     * \see Gradient
     * \see RadialInterpolation
     * 
     * */
    void    Hessian( double * hess /*!< [out] dense symmetric matrix stored in C-style \c{row_major} format */);

    /*!\brief Compute the Hessian via finite difference method (based on the gradient)
     * 
     * The output consists of a one-dimensional array representing the matrix according to
     * \f[
     *      {m\_res}_{ (i-1)d + j } = J_{ij} = J_{ji} = \frac{ \partial^2 I}{ \partial x_i \, \partial x_j }
     * \f]
     * 
     * \attention The computation requires a previous call to Interpolate and Gradient.
     * 
     * \attention After this calls to Gradient/Hessian work on the last perturbed value.
     * 
     * \see Interpolate
     * \see Gradient
     * \see RadialInterpolation
     * */
    void    HessianFDM( double * hess, /*!< [out] symmetric matrix stored in C-style \c{row_major} format */
    const double dx /*!< [in] perturbation parameter */ ); 
};



#endif /* __CONCENTRIC_INTERPOLATION_H__ */

