#ifndef _TANGENTIAL_INTERPOLATION_H_
#define _TANGENTIAL_INTERPOLATION_H_

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
#include <util.h>
#include <cblas.h>

/** Tangential Interpolation is the classical interpolation on spheres by means of spherical basis functions.
 * The class TangentialInterpolation provides all functionalities required for the setup and conduction of
 * such a spherical interpolation. The main function is Weights, which returns \f$ \underline{\underline{K}}^{-1}\underline{\zeta} \f$,
 * where \f$ \underline{\underline{K}} \f$ is the kernel matrix and \f$ \underline{\zeta} \f$ is the vector
 * of kernel values.
 *
 * At the moment, the Gaussian kernel \f$ k_i(\underline{N})=\exp(-\gamma\mathrm{acos}(\underline{N}\cdot\underline{N_i})^2) \f$ is employed,
 * but the user is free to implement the kernel of their choice. The spherical supporting points \f$ \{\underline{N_i}\}_{i=1,\dots,N} \$
 * are also called \emph{directions}.
 *
 * For Concentric Interpolation, this is used in combination with Radial Interpolation.
 *
 * \see RadialInterpolation
 */

class TangentialInterpolation;
// class TangentialInterpolationPseudoSym;

class TangentialInterpolation {
private:

    /*! \brief (Re-)Compute the kernel matrix and its inverse
     *
     * This member function evaluates the kernel matrix and computes its inverse.
     * A call to this function **must** be made before the interpolation function can be
     * called.
     *
     * */
    void    InitializeKernelMethod( );

    /*! \brief Re-compute the kernel matrix and its inverse
     *
     * This does not re-allocate memory or change any member variables except m_K and m_Ki,
     * in contrast to InitializeKernelMethod.
     *
     * \see InitializeKernelMethod
     *
     * */
    virtual void    RecomputeKernelMatrix( );

protected:

    bool        sym;                //!< exploit point symmetry (only one ray of data points must be provided)
    bool        init;               //!< initialization flag
    int         N_alloc;            //!< dimension of the pre-allocated memory
    int         N;                  //!< number of actually provided directions
    int         D;                  //!< dimension of the inputs
    double      gamma;              //!< kernel width parameter of the Gaussian kernel function
    double      lambda;             //!< regression parameter: controls condition number of kernel matrix. usually not needed.

    double      *x;                 //!< input vector x after normalization, i.e. d = x / r the direction of the input vector
    double      *theta, *zeta, *xi, //!< theta_i:   d_i * d, xi: acos(theta), zeta_i: kernel function evaluated at xi_i
//                 *sin_xi,            //!< sine of xi_i (re-used several times)
                *zeta_tilde,        //!< for the symmetric case: kernel function at (PI-xi_i) (see paper)
                *zeta_star;         //!< for the symmetric case: kernel function at (PI-xi_i) (see paper)
    int         * w_i;              //!< integer working array for linear solver (IPIV in LAPACK)
    bool        * active;           //!< mark training directions which are active or not. active means not too close to the query direction, i.e. 1/sin(xi) is bound. inactive training directions will be assumed to have a minimum angular distance to the query direciton
    // matrix like 2D arrays (in terms of 1D row_major arrays). these are prefixed with "m_" to emphasize their multidimensional nature.
    double      *m_K,               //!< dense kernel matrix
                *m_Kf,              //!< the LDL factorization of the kernel matrix. used for solving systems of the form K * Ki_a = a for Ki_a
                *m_X;               //!< training directions, i.e. supporting points on the sphere (pre-allocated with size N_alloc)

    static const double small;
    static const double theta_max;

    void zero_pointers();           //!< initialize all pointers to zero (before doing anything else)

public:

    TangentialInterpolation( const bool a_sym = false /*!< [in] use symmetrization */ );
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

    //!< clean up and destroy the object
    ~TangentialInterpolation();

    void    Free(); //!< free allocated memory \see ~ConcentricInterpolation()


    void    Allocate( int a_N_alloc,        /*!< [in] maximum admissible number of directions used for input data */
                      const int a_D     /*!< [in] dimension of the input vector(s) */ );
    //!< \brief pre-allocate memory for (up to) \c a_N_alloc directions in R^\c a_D

      /*! \brief Compute weights for the direction of the vector \c X (not necessarily normalized)
     *
     * If the vector is zero, then the unit vector e_1 is returned.
     * Otherwise, the (symmetric) kernel function is evaluated and
     * K^-1 * zeta is returned in \c W.
     *
     * optionally, zeta, dzeta and ddzeta can be output
     *
     * */
    void Weights(
        double * a_W,       /*!< [out] vector \c W of weights */
        const int n_dirs,   /*!< [in] number of directions in a_x */
        const double * a_x, /*!< [in] n_dirs x dim matrix of directions (row-wise) \c X */
        double * a_zeta,    /*!< [out, optional] vector \c W of zeta values (if not needed set to NULL) */
        double * a_dzeta,   /*!< [out, optional] vector \c W of dzeta values (if not needed set to NULL) */
        double * a_ddzeta   /*!< [out, optional] vector \c W of ddzeta values (if not needed set to NULL) */
                );

//     void Weights(
//         int    n_vec,       /*!< [in]  number of input vector */
//         double * a_W,       /*!< [out] (n x n_vec) matrix \c W of weights */
//         const double * a_x, /*!< [in]  (n_vec x D) matrix containing directions as /direction \c X */
//         double * w_d        /* sufficiently large working array of doubles */ );

    //! short-cut to Weights( ... ) with n_dirs=1 [see above]
    void operator() (
        double * a_W,       /*!< [out] vector \c W of weights */
        const double * a_x, /*!< [in] single direction \c X */
        double * a_zeta,    /*!< [out, optional] vector \c W of zeta values (if not needed set to NULL) */
        double * a_dzeta,   /*!< [out, optional] vector \c W of dzeta values (if not needed set to NULL) */
        double * a_ddzeta   /*!< [out, optional] vector \c W of ddzeta values (if not needed set to NULL) */
                ) { Weights( a_W, 1, a_x, a_zeta, a_dzeta, a_ddzeta ); }

    /*! \brief Add a new direction
     *
     * Note that the directions **must not** be identical (or parallel in the symmetric case).
     * otherwise, the kernel matrix becomes singular and the program will fail without error handling.
     * */
    void    AddDirection(
        const double * a_X      //!< [in] (unit) vector of dimension D; direction along which the data is provided [normalization is carried out internally]
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
     * \see GetGamma
     *
     * */
    void SetGamma( const double a_gamma /*!< [in] new kernel parameter */ );

    /*! \brief Set a regression parameter.
     * This value is added to the diagonal of the kernel matrix, improving the
     * condition number and worsening the accuracy. */
    void SetLambda( const double a_lambda );

}; /* class TangentialInterpolation */


// class TangentialInterpolationPseudoSym : public TangentialInterpolation {
// private:
//     double *m_Kdiff, *m_Kdiff_f;    //!< difference kernel matrix and its factorization (LDL)
//     int         * w_diff_i;         //!< integer working array for linear solver (IPIV in LAPACK)
//     double      * w_s;              //!< zeta_star and symmtric weights
//     void zero_pointers();
// public:
//     TangentialInterpolationPseudoSym();
//     ~TangentialInterpolationPseudoSym();
//     void Free();
//     void Allocate( int a_N_alloc, const int a_D );
//
//     void RecomputeKernelMatrix();
//
//     void Weights(
//         double * a_W,       /*!< [out] vector \c W of weights */
//         const int n_dirs,   /*!< [in] number of vectors contained in \c X */
//         const double * a_x, /*!< [in] vector(s) [row-matrix of size n_dirs x D] \c X */
//         double * a_zeta,    /*!< [out, optional] vector \c W of zeta values (if not needed set to NULL) */
//         double * a_dzeta,   /*!< [out, optional] vector \c W of dzeta values (if not needed set to NULL) */
//         double * a_ddzeta   /*!< [out, optional] vector \c W of ddzeta values (if not needed set to NULL) */
//                 );
// //    void Weights(
// //         int    n_vec,       /*!< [in]  number of input vector */
// //         double * a_W,       /*!< [out] (n x n_vec) matrix \c W of weights */
// //         const double * a_x, /*!< [in]  (n_vec x D) matrix containing directions as /direction \c X */
// //         double * w_d        /* sufficiently large working array of doubles */ );
//
//    const double * const SymWeights() const { return w_s; }
//
// };

#endif /* _TANGENTIAL_INTERPOLATION_H_ */
