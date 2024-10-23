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

#ifndef _TANGENTIAL_INTERPOLATION_H_
#define _TANGENTIAL_INTERPOLATION_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <util.h>

class TangentialInterpolation;

/** \brief This is one of the two main ingredients for ConcentricInterpolation, besides RadialInterpolation.
 *
 * Tangential Interpolation is the classical interpolation on spheres by means of spherical basis functions.
 * The class TangentialInterpolation provides all functionalities required for the setup and conduction of
 * such a spherical interpolation. The main function is KernelVector, which returns \f$ \underline{k}(\mathbf{n}) \f$,
 * the vector of kernel values from equation (CI) of the article.
 *
 * At the moment, the Gaussian kernel \f$ k_i(\underline{N})=\exp(-\gamma\mathrm{acos}(\underline{N}\cdot\underline{N_i})^2) \f$ is employed,
 * but the user is free to implement the kernel of their choice. The spherical supporting points \f$ \{\underline{N_i}\}_{i=1,\dots,N_{\rm dir}^{\rm supp}} \subset \mathbf{S}^{D_{\rm inp}-1} \f$
 * are also called \e directions (in \f$ \mathbf{R}^{D_{\rm inp}}\f$).
 *
 * To get started, call Setup().
 *
 * \see RadialInterpolation
 * \see ConcentricInterpolation
 *
 */

class TangentialInterpolation {
private:
    /*! \brief (Re-)Computes the kernel matrix and its inverse
     *
     * Requires prior call to Allocate().
     *
     * This is automatically called during KernelVector(), if it has not been called
     * before or if any of the following functions was executed since the last
     * time:
     *  - Allocate()
     *  - SetSymmetric()
     *  - SetGamma()
     *  - SetLambda()
     *
     * \todo complete documentation of TangentialInterpolation
     *
     * */
    virtual void    ComputeKernelMatrix( );

protected:

    bool        sym;                //!< exploit point symmetry (only one ray of data points must be provided)
    bool        init;               //!< initialization flag
    int         N_dir_alloc;        //!< dimension of the pre-allocated memory
    int         N_dir_supp;         //!< number of \e actually provided directions, i.e. support points of the spherical basis functions
    int         D_inp;              //!< dimension of the Euclidean space in which the directions are defined
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
                *m_X;               //!< training directions, i.e. supporting points on the sphere (pre-allocated with size N_dir_alloc)

    static const double small;
    static const double theta_max;

    void zero_pointers();           //!< initialize all pointers to zero (before doing anything else)

public:

    /*! Default constructur
     *
     *  \todo comment default constructor of TangentialInterpolation
     */
    TangentialInterpolation();

    //! Clean up and destroy the object
    ~TangentialInterpolation();

    void    Free(); //!< free allocated memory \see ~ConcentricInterpolation()

    //! calls Allocate(), AddDirection() and SetGamma() in the correct order. SetGamma() may also be called afterwards, the argument \p gamma is optional. SetLambda() is usually not required, but may also be called afterwards.
    void Setup( const int       a_N_dir_supp,        //!< number of \e desired directions, i.e. support points of the spherical basis functions \see TangentialInterpolation::N_dir_supp
                const int       a_D_inp,        //!< dimension of Euclidean space, i.e. the directions are \p N_dir_supp points on \f$ \mathbf{S}^{D_{\rm inp}-1} \f$
                const double * const * a_directions,   //!< d-th component of n-th direction --> <tt> a_directions[n][d] </tt>
                const bool      a_sym=false,//!< [\e optional] symmetry flag, defaluts to false meaning no symmetry will be considered. See SetSymmetric().
                const double    gamma=0.    //!< [\e optional] kernel width parameter, defaults to zero meaning it has to be set later. Calls SetGamma().
         );

    /*! \brief pre-allocate memory for (up to) \c a_N_alloc directions in R^\c a_D_inp */
    void    Allocate( int a_N_alloc,    /*!< [in] maximum admissible number of directions used for input data */
                      const int a_D_inp     /*!< [in] dimension of the Euclidean space in which the directions are defined */ );

     /*! \brief Compute and return the kernel vector \f$\underline{k}(\mathbf{n})\f$, see equation (CI) of the paper.
     *
     * The kernel vector is called \p zeta in this code.
     * The derivatives \p dzeta and \p ddzeta are always computed, and can optionally be output.
     *
     * This function does \e not depend on or use the kernel matrix \f$\underline{\underline{K}}\f$ in any way.
     *
     * \attention If an evaluation direction (i.e. row of \c X) has very small norm (e.g. numerical zero), then the evaluation direction is set to \f$ [1, 0, \ldots , 0] \in\mathbf{R}^{D_{\rm inp}} \f$
     *
     * */
    void KernelVector(
        const int n_dirs_eval,  /*!< [in] number of evaluation directions \p n_dirs */
        const double * a_x,     /*!< [in] row-major pseudo-matrix (\p n_dirs x \p D_inp) of evaluation directions \c X */
        double * a_zeta,        /*!< [out] row-major pseudo-matrix (\p n_dirs x \p N_dir_supp) containing the \f$ \underline{\zeta} \f$ vectors */
        double * a_dzeta =0,    /*!< [out, optional] \todo implement \f$ {\rm d}\underline{\zeta} \f$*/
        double * a_ddzeta=0     /*!< [out, optional] \todo implement \f$ {\rm d}^2\underline{\zeta} \f$ */
                );

    /*! \brief Add a new direction
     *
     * Note that the directions **must not** be identical in any case, and **must not** be parallel in the symmetric case.
     * Otherwise, the kernel matrix becomes singular and the program will fail without error handling.
     * */
    void    AddDirection(
        const double * a_X      //!< [in] (unit) vector of dimension D_inp; direction along which the data is provided [normalization is carried out internally]
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
     *  - if the symmetric interpolation should have \f$ N_{\rm dir}^{\rm supp} \f$ supporting points on \f$ \mathbf{S}^{D_{\rm inp}-1} \f$, then
     * <em> only \f$ N_{\rm dir}^{\rm supp}/2 \f$ points must be provided</em>
     *  - the other \f$ N_{\rm dir}^{\rm supp}/2 \f$ points are assumed to be the \e antipodes of the provided ones, and are considered implicitly
     *  - <em>DO NOT PROVIDE THE ANTIPODES</em>
     *  - \f$ N_{\rm dir}^{\rm supp} \f$ must be even
     *
     * \see GetSymmetric
     * \see AddDirection
    */
    void SetSymmetric(  const bool a_sym,
                        const bool a_recompute_kernel_matrix = false    /*!< [in] call ComputeKernelMatrix() now? */
    );

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
    void SetGamma(  const double a_gamma,                           /*!< [in] new kernel parameter */
                    const bool a_recompute_kernel_matrix = false,   /*!< [in] call ComputeKernelMatrix() now? */
                    const bool a_quiet = false                      /*!< [in] if true, then don't print notification */
    );

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
    void SetLambda( const double a_lambda,
                    const bool a_recompute_kernel_matrix = false    /*!< [in] call ComputeKernelMatrix() now? */
    );

    /*! \brief Returns number of directions \p N_dir_supp. */
    inline int GetN_dir_supp( ) const
    {
        return N_dir_supp;
    }

    /*! \brief Returns number of spatial dimensions \p D_inp. */
    inline int GetD_inp( ) const
    {
        return D_inp;
    }

    /*! \brief Same as GetD_inp(). */
    inline int GetDimCoords( ) const
    {
        return D_inp;
    }

    /*! \brief return kernel matrix \f$ \underline{\underline{K}} \f$ */
    inline const double * const GetKernelMatrix() const
    {
        UTILITY::assert_msg( init, "ERROR: initialize TangentialInterpolation before calling GetKernelMatrix. Forgot to call ComputeKernelMatrix?\n");
        return m_K;
    }

    /*! \brief return LDL factorization of kernel matrix \f$ \underline{\underline{K}} \f$ */
    inline const double * GetKernelMatrixFactorization() const
    {
        UTILITY::assert_msg( init, "ERROR: initialize TangentialInterpolation before calling GetKernelMatrixFactorization. Forgot to call ComputeKernelMatrix?\n");
        return m_Kf;
    }

    /*! \brief return permutation indices of the factorization TangentialInterpolation::w_i */
    inline const int * GetPermutation() const
    {
        UTILITY::assert_msg( init, "ERROR: initialize TangentialInterpolation before calling GetPermutation. Forgot to call ComputeKernelMatrix?\n");
        return w_i;
    }

    /*! \brief return initialization flag */
    inline const bool GetInit() const
    {
        return init;
    }

}; /* class TangentialInterpolation */


#endif /* _TANGENTIAL_INTERPOLATION_H_ */
