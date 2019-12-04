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

#include <tangential_interpolation.h>
using namespace UTILITY;

/* *************************************************************************************** */
const double TangentialInterpolation::small     = 1.e-12;   //!< small number; used, e.g., in order to prevent division by zero errors
const double TangentialInterpolation::theta_max = 0.99995;  //!< required in order to regularize the derivative of zeta
/* *************************************************************************************** */
TangentialInterpolation::TangentialInterpolation()
{
    FastAcos::Init();
    sym         = false;
    init        = false;
    N_dir_supp  = 0;
    N_dir_alloc = 0;
    lambda      = 0.;
    D_inp       = 0;
    gamma       = -1.;
    zero_pointers();
}
/* *************************************************************************************** */
void TangentialInterpolation::zero_pointers()
{
    // initialize pointers to NULL

    // protected members:
    w_i         = 0;
    m_K         = 0;
    m_Kf        = 0;
    m_X         = 0;
    theta       = 0;
    zeta        = 0;
    zeta_star   = 0;
    zeta_tilde  = 0;
    active      = 0;
//     sin_xi      = 0;
    xi          = 0;
    x           = 0;
    init        = false;
}
/* *************************************************************************************** */
void TangentialInterpolation::Allocate( int a_N_alloc, const int a_D_inp )
{
    assert_msg( a_N_alloc>0, "ERROR in TangentialInterpolation::Allocate: N_dir_alloc > 0 required\n");
    assert_msg( a_D_inp>0, "ERROR in TangentialInterpolation::Allocate: D_inp > 0 required\n");
    Free();

    N_dir_alloc = a_N_alloc;
    D_inp       = a_D_inp;
    N_dir_supp  = 0;

    active  = new bool [ N_dir_alloc ];
    assert_msg( active != 0, "ERROR in TangentialInterpolation::Allocate: could not allocate memory for 'active'\n");

    m_X     = alloc_array( N_dir_alloc*D_inp );
    x       = alloc_array( D_inp );
//     sin_xi  = alloc_array( N_dir_alloc );
    theta   = alloc_array( N_dir_alloc );
    xi      = alloc_array( N_dir_alloc );
    zeta    = alloc_array( N_dir_alloc );

    // make sure the working array has the appropriate size
    if( w_i != 0 ) { delete [] w_i; w_i = 0; }
    w_i     = new int [ 32*N_dir_alloc*N_dir_alloc ]; memset(w_i, 0, sizeof(int)*size_t(32*N_dir_alloc*N_dir_alloc));
    m_K     = new double [ N_dir_alloc*N_dir_alloc ]; memset(m_K, 0, sizeof(double)*size_t(N_dir_alloc*N_dir_alloc));
    m_Kf    = new double [ N_dir_alloc*N_dir_alloc ]; memset(m_Kf,0, sizeof(double)*size_t(N_dir_alloc*N_dir_alloc));
    if( sym ) {
        zeta_tilde  = alloc_array( N_dir_alloc );
        zeta_star   = alloc_array( N_dir_alloc );
    }
}
/* *************************************************************************************** */
TangentialInterpolation::~TangentialInterpolation() {
    Free();
}
/* *************************************************************************************** */
void TangentialInterpolation::Free()
{
    delete [] w_i; w_i = 0;
    init    = false;
    sym     = false;
    gamma   = 0.;
    lambda  = 0.;
    N_dir_alloc = 0;
    N_dir_supp       = 0;
    free_array( &m_X );
    free_array( &m_K );
    free_array( &m_Kf );
    free_array( &x );
    if( active != 0 )
    {
        delete [] active;
        active = 0;
    }
//     free_array( &sin_xi );
    free_array( &theta );
    free_array( &zeta );
    free_array( &zeta_tilde );
    free_array( &zeta_star );
    free_array( &xi );
}
/* *************************************************************************************** */
void TangentialInterpolation::Setup(
            const int       a_N_dir_supp,
            const int       a_D_inp,
            const double * const * a_directions,
            const bool      a_sym,
            const double    a_gamma
        )
{
    Allocate( a_N_dir_supp, a_D_inp );
    for( int i_dir=0; i_dir<a_N_dir_supp; i_dir++)
        AddDirection( a_directions[i_dir] );
    SetSymmetric(a_sym);
    if( a_gamma>0. )
        SetGamma( a_gamma, true );
    printf("# set up Tangential Interpolant for N_dir_supp = %i, D_inp = %i, sym = %i, gamma = %e\n", a_N_dir_supp, a_D_inp, a_sym, a_gamma);
}
/* *************************************************************************************** */
void TangentialInterpolation::AddDirection( const double * a_X )
{
    assert_msg( N_dir_supp < N_dir_alloc, "ERROR in TangentialInterpolation::AddInterpolationData: allocation size too small (trying to add another entry although n=N_dir_alloc)\n");
    init = false; // reset initialization flag since kernel matrix needs to be re-computed!

    const double l = norm( a_X, D_inp );
    for(int d=0; d<D_inp; d++) m_X[N_dir_supp*D_inp+d] = a_X[d] / l; // make sure the length is 1 for the directions
    N_dir_supp++; // increment the counter for the dimension of m_X, i.e. the number of actually provided training directions a.k.a. support points for the spherical basis functions
}
/* *************************************************************************************** */
void TangentialInterpolation::KernelVector(
        const int n_dirs_eval,
        const double * a_x,
        double * a_zeta,
        double * a_dzeta,
        double * a_ddzeta
                )
{
    const bool b_dzeta      = (a_dzeta != 0);
    const bool b_ddzeta     = (a_ddzeta != 0);

    if( !init ) { ComputeKernelMatrix(); }

    int offset = 0;
    for( int i_dir=0; i_dir < n_dirs_eval; i_dir++)
    {

        double radius = norm( a_x + i_dir*D_inp, D_inp );
        if( radius>1e-8 )
            for(int d=0;d<D_inp;d++) x[d] = a_x[d+i_dir*D_inp]/radius; // normalize inputs, i.e. direction of input
        // ATTENTION in the case of a very small amplitude, change input direction to [1,0,...0] and radius to 1.
        else
        {
            radius  = 1.;
            x[0]    = 1.;
            for(int d=1;d<D_inp;d++) x[d] = 0.;
        }

        // compute theta = cos(xi)
        MatVecMul( m_X, x, theta, N_dir_supp, D_inp );

        for(int n=0; n<N_dir_supp; n++) {
            xi[n] = safeAcos(theta[n]);
            if( theta[n] > theta_max ) {
                active[n]   = false;
//                 sin_xi[n]   = 1e-16;
            }
            else if ( theta[n] < -theta_max ) {
                active[n]   = false;
//                 sin_xi[n]   = 1e-16;
            }
            else {
                active[n]   = true;
//                 sin_xi[n]   = sqrt( 1.-theta[n]*theta[n]); // sin(xi[n]);
            }
        }
        // symmetric case:
        if( sym )
        {
            // the intermediate vectors zeta_tilde and zeta_star are re-used in the gradient computation
    #pragma unroll (4)
            for( int n=0; n<N_dir_supp; n++ ) {
                zeta[n]         = exp(- gamma * xi[n]*xi[n] );
                zeta_tilde[n]   = exp(- gamma * (PI-xi[n]) * (PI-xi[n]) );
                zeta_star[n]    = zeta_tilde[n] + zeta[n];
                if( b_dzeta || b_ddzeta )
                {
                    if(active[n])
                    {
                        if(b_dzeta)     a_dzeta[n+offset]    = - 2.*gamma*( zeta[n]*xi[n] - (PI-xi[n])*zeta_tilde[n] );
                        if(b_ddzeta)    a_ddzeta[n+offset]   = -2.*gamma * (
                                zeta[n]         * ( 1. - 2.*gamma*xi[n]*xi[n]  )
                            +   zeta_tilde[n]   * ( 1. - 2.*gamma*(PI-xi[n])*(PI-xi[n]) )
                                                );
                    }
                    else
                    {
                        if(b_dzeta)     a_dzeta[n+offset]    = 0.;
                        if(b_ddzeta)    a_ddzeta[n+offset]   = 0.;
                    }
                }
            }
        }
        else
        {
        // for the non-symmetric case:
    #pragma unroll (4)
            for(int n=0; n<N_dir_supp; n++)
            {
                zeta[n] = exp(- gamma * xi[n]*xi[n] );
                if( b_dzeta || b_ddzeta )
                {
                    if( active[n] )
                    {
                        if(b_dzeta)     a_dzeta[n+offset]    = - 2.*gamma*zeta[n]*xi[n];
                        if(b_ddzeta)    a_ddzeta[n+offset]   = -2.*gamma * zeta[n] * ( 1. - 2.*gamma*xi[n]*xi[n]  );
                    }
                    else
                    {
                        if(b_dzeta)     a_dzeta[n+offset]    = 0.;
                        if(b_ddzeta)    a_ddzeta[n+offset]   = 0.;
                    }
                }
            }
        }

#pragma unroll (4)
        for(int n=0; n<N_dir_supp; n++)
            a_zeta[n+offset]=zeta[n];

        offset += N_dir_supp;
    } // for( int i_dir=0; i_dir < n_dirs_eval; i_dir++)
}
/* *************************************************************************************** */
void TangentialInterpolation::SetSymmetric( const bool a_sym, const bool a_recompute_kernel_matrix )
{
    sym  = a_sym;
    init = false;
    if( a_recompute_kernel_matrix )
        ComputeKernelMatrix();
}
/* *************************************************************************************** */
void TangentialInterpolation::SetGamma( const double a_gamma, const bool a_recompute_kernel_matrix, const bool a_quiet )
{
    assert_msg( fabs(a_gamma)>1e-8, "ERROR in TangentialInterpolation::SetGamma: kernel width parameter must be significantly positive\n");
    gamma = a_gamma;
    init  = false;
    if( !a_quiet )
        printf("# Tangential Interpolation: set gamma to %f\n", gamma);
    if( a_recompute_kernel_matrix )
        ComputeKernelMatrix();
}
/* *************************************************************************************** */
void TangentialInterpolation::SetLambda( const double a_lambda, const bool a_recompute_kernel_matrix )
{
    assert_msg( a_lambda >= 0., "ERROR in TangentialInterpolation::SetLambda: regression parameter must be non-negative\n");
    lambda = a_lambda;
    init   = false;
    printf("# Tangential Interpolation: set lambda to %f\n", lambda);
    if( a_recompute_kernel_matrix )
        ComputeKernelMatrix();
}
/* *************************************************************************************** */
void TangentialInterpolation::ComputeKernelMatrix()
{
    assert_msg( gamma>1e-3, "ERROR in TangentialInterpolation::ComputeKernelMatrix: gamma must be significantly greater than zero!\n");
    double tmp_alpha = 0.;
    for(int n1=0;n1<N_dir_supp;n1++)
    {
        for(int n2=n1;n2<N_dir_supp;n2++)
        {
            tmp_alpha  = safeAcos( VecVecMul( m_X + n1*D_inp, m_X + n2*D_inp, D_inp ) );
            if ( sym )
                    m_K[n1*N_dir_supp+n2]  = exp(-gamma*tmp_alpha*tmp_alpha)  +  exp(-gamma*(PI-tmp_alpha)*(PI-tmp_alpha));
            else
                    m_K[n1*N_dir_supp+n2]  = exp(-gamma*tmp_alpha*tmp_alpha);
            m_K[n2*N_dir_supp+n1] = m_K[n1*N_dir_supp+n2];
        }
        m_K[n1*(N_dir_supp+1)] += lambda; // add regression parameter
    }

    // store LDL factorization of kernel matrix
    Factorize( m_K, m_Kf, w_i, N_dir_supp );

    // ready to go
    init = true;

}
/* *************************************************************************************** */
