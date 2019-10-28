#ifndef _TANGENTIAL_INTERPOLATION_H_
#define _TANGENTIAL_INTERPOLATION_H_

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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <lapacke.h>
#include <util.h>
#include <cblas.h>

class TangentialInterpolation;

/** Tangential Interpolation is the classical interpolation on spheres by means of spherical basis functions.
 * The class TangentialInterpolation provides all functionalities required for the setup and conduction of
 * such a spherical interpolation. The main function is Weights, which returns \f$ \underline{\underline{K}}^{-1}\underline{\zeta} \f$,
 * where \f$ \underline{\underline{K}} \f$ is the kernel matrix and \f$ \underline{\zeta} \f$ is the vector
 * of kernel values.
 *
 * At the moment, the Gaussian kernel \f$ k_i(\underline{N})=\exp(-\gamma\mathrm{acos}(\underline{N}\cdot\underline{N_i})^2) \f$ is employed,
 * but the user is free to implement the kernel of their choice. The spherical supporting points \f$ \{\underline{N_i}\}_{i=1,\dots,N} \subset \mathbf{R}^D \f$
 * are also called \e directions.
 *
 * For Concentric Interpolation, this is used in combination with RadialInterpolation.
 *
 * To get started, call Setup().
 *
 * \see RadialInterpolation
 * \see ConentricInterpolation
 *
 */

class TangentialInterpolation {
private:

    /*! \brief (Re-)Computes the kernel matrix and its inverse
     *
     * Requires prior call to Allocate().
     *
     * This is automatically called during Weights(), if it has not been called
     * before or if any of the following functions was executed since the last
     * time:
     *  - Allocate()
     *  - SetSymmetric()
     *  - SetGamma()
     *  - SetLambda()
     *
     * TODO
     *
     * */
    virtual void    ComputeKernelMatrix( );

protected:

    bool        sym;                //!< exploit point symmetry (only one ray of data points must be provided)
    bool        init;               //!< initialization flag
    int         N_alloc;            //!< dimension of the pre-allocated memory
    int         N;                  //!< number of \e actually provided directions, i.e. support points of the spherical basis functions
    int         D;                  //!< dimension of the Euclidean space in which the directions are defined
    double      gamma;              //!< kernel width parameter of the Gaussian kernel function
    double      lambda;             //!< regression parameter: controls condition number of kernel matrix. usually not needed.

    double      *x;                 //!< input vector x after normalization, i.e. d = x / r the direction of the input vector
    double      *theta, *zeta, *xi, //!< theta_i:   d_i * d, xi: acos(theta), zeta_i: kernel function evaluated at xi_i
//                 *sin_xi,            //!< sine of xi_i (re-used several times)
                *zeta_tilde,        //!< for the symmetric case: kernel function at (PI-xi_i) (see paper)
                *zeta_star;         //!< for the symmetric case: kernel function at (PI-xi_i) (see paper)
    int         * w_i;              //!< permutation array for factorization, this is \p IPIV of LAPACK's function \p dsytrf
    bool        * active;           //!< mark training directions which are active or not. active means not too close to the query direction, i.e. 1/sin(xi) is bound. inactive training directions will be assumed to have a minimum angular distance to the query direciton
    // matrix like 2D arrays (in terms of 1D row_major arrays). these are prefixed with "m_" to emphasize their multidimensional nature.
    double      *m_K,               //!< dense kernel matrix
                *m_Kf,              //!< the LDL factorization of the kernel matrix. used for solving systems of the form K * Ki_a = a for Ki_a
                *m_X;               //!< training directions, i.e. supporting points on the sphere (pre-allocated with size N_alloc)

    static const double small;
    static const double theta_max;

    void zero_pointers();           //!< initialize all pointers to zero (before doing anything else)

public:

    /*! Default constructur, sets
     *  - TODO
     */
    TangentialInterpolation();

    //! Clean up and destroy the object
    ~TangentialInterpolation();

    void    Free(); //!< free allocated memory \see ~ConcentricInterpolation()

    //! calls Allocate(), AddDirection() and SetGamma() in the correct order. SetGamma() may also be called afterwards, the argument \p gamma is optional. SetLambda() is usually not required, but may also be called afterwards.
    void Setup( const int       a_N,        //!< number of \e desired directions, i.e. support points of the spherical basis functions \see TangentialInterpolation::N
                const int       a_D,        //!< dimension of Euclidean space, i.e. the directions are \p N points on \f$ \mathbf{S}^{D-1} \f$
                const double * const * a_directions,   //!< d-th component of n-th direction --> <tt> a_directions[n][d] </tt>
                const bool      a_sym=false,//!< [\e optional] symmetry flag, defaluts to false meaning no symmetry will be considered. See SetSymmetric().
                const double    gamma=0.    //!< [\e optional] kernel width parameter, defaults to zero meaning it has to be set later. Calls SetGamma().
         );

    /*! \brief pre-allocate memory for (up to) \c a_N_alloc directions in R^\c a_D */
    void    Allocate( int a_N_alloc,    /*!< [in] maximum admissible number of directions used for input data */
                      const int a_D     /*!< [in] dimension of the Euclidean space in which the directions are defined */ );

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
        double * a_zeta  =0,/*!< [out, optional] vector \c W of zeta values */
        double * a_dzeta =0,/*!< [out, optional] vector \c W of dzeta values */
        double * a_ddzeta=0 /*!< [out, optional] vector \c W of ddzeta values */
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
        double * a_zeta  =0,/*!< [out, optional] vector \c W of zeta values */
        double * a_dzeta =0,/*!< [out, optional] vector \c W of dzeta values */
        double * a_ddzeta=0 /*!< [out, optional] vector \c W of ddzeta values */
                ) { Weights( a_W, 1, a_x, a_zeta, a_dzeta, a_ddzeta ); }

     /*! \brief Compute and return the kernel vector \f$\underline{k}(\mathbf{n})\f$, see equation (CI) of the paper.
     *
     * The kernel vector is called \p zeta in this code.
     * The derivatives \p dzeta and \p ddzeta are always computed, and can optionally be output.
     *
     * This function does \e not depend on or use the kernel matrix \f$\underline{\underline{K}}\f$ in any way.
     *
     * \attention If an evaluation direction (i.e. row of \c X) has very small norm (e.g. numerical zero), then the evaluation direction is set to \f$ [1, 0, \ldots , 0] \in\mathbf{R}^N \f$
     *
     * */
    void KernelVector(
        const int n_dirs_eval,  /*!< [in] number of evaluation directions in a_x */
        const double * a_x,     /*!< [in] row-major pseudo-matrix (\p n_dirs x \p D) of evaluation directions \c X */
        double * a_zeta,        /*!< [out] row-major pseudo-matrix (\p n_dirs x \p N) containing the \f$ \underline{\zeta} \f$ vectors */
        double * a_dzeta =0,    /*!< [out, optional] TODO \todo */
        double * a_ddzeta=0     /*!< [out, optional] TODO \todo  */
                );

    /*! \brief Add a new direction
     *
     * Note that the directions **must not** be identical in any case, and **must not** be parallel in the symmetric case.
     * Otherwise, the kernel matrix becomes singular and the program will fail without error handling.
     * */
    void    AddDirection(
        const double * a_X      //!< [in] (unit) vector of dimension D; direction along which the data is provided [normalization is carried out internally]
    );

    /*! \brief Returns the the symmetry flag.
     *
     * \see SetSymmetric
     *
     * */
    inline bool GetSymmetric( ) const
    {
        return sym;
    }

    /*! \brief Set the symmetry parameter. If \p a_sym is true, then the interpolation  is symmetrized,
     * i.e. \f$ \widetilde{f}( \underline{N} ) = \widetilde{f}( -\underline{N} ) \f$ is strictly enforced at (almost) no additional computational expense.
     * More precisely the dimension of the kernel matrix is not increased and the evaluation of the
     * interpolation (and of the gradients) involves only few additional operations.
     *
     * \attention
     *  - if the symmetric interpolation should have \f$ N \f$ supporting points on \f$ \mathbf{S}^{D-1} \f$, then
     * <em> only \f$ N/2 \f$ points must be provided</em>
     *  - the other \f$ N/2 \f$ points are assumed to be the \e antipodes of the provided ones, and are considered implicitly
     *  - <em>DO NOT PROVIDE THE ANTIPODES</em>
     *  - \f$ N \f$ must be even
     *
     * \see GetSymmetric
     * \see AddDirection
    */
    void SetSymmetric( const bool a_sym );

    /*! \brief Returns the kernel width parameter.
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
     * \see GetGamma
     *
     * */
    void SetGamma( const double a_gamma /*!< [in] new kernel parameter */ );

    /*! \brief Returns the regression parameter.
     *
     * \see SetLambda
     *
     * */
    inline double GetLambda( ) const
    {
        return lambda;
    }

    /*! \brief Set the regression parameter.
     * This value is added to the diagonal of the kernel matrix, improving the
     * condition number and compromising the accuracy. Usually unnecessary.
     *
     * \see GetLambda
     */
    void SetLambda( const double a_lambda );

    /*! \brief Returns number of directions N. */
    inline int GetN( ) const
    {
        return N;
    }

    /*! \brief return kernel matrix \f$ \underline{\underline{K}} \f$ */
    inline const double * const GetKernelMatrix() const
    {
        return m_K;
    }

    /*! \brief return LDL factorization of kernel matrix \f$ \underline{\underline{K}} \f$ */
    inline const double * GetKernelMatrixFactorization() const
    {
        return m_Kf;
    }

    /*! \brief return permutation indices of the factorization TangentialInterpolation::w_i */
    inline const int * GetPermutation() const
    {
        return w_i;
    }

    /*! \brief return initialization flag */
    inline const bool GetInit() const
    {
        return init;
    }

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
//     void ComputeKernelMatrix();
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
