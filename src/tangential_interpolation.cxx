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
    N           = 0;
    N_alloc     = 0;
    lambda      = 0.;
    D           = 0;
    gamma       = -1.;
    zero_pointers();
}
/* *************************************************************************************** */
void TangentialInterpolation::zero_pointers()
{
    // initialize pointers to NULL

    // private members:
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
}
/* *************************************************************************************** */
void TangentialInterpolation::Allocate( int a_N_alloc, const int a_D )
{
    assert_msg( a_N_alloc>0, "ERROR in TangentialInterpolation::Allocate: N_alloc > 0 required\n");
    assert_msg( a_D>0, "ERROR in TangentialInterpolation::Allocate: D > 0 required\n");
    Free();

    N_alloc = a_N_alloc;
    D       = a_D;
    N       = 0;

    active  = new bool [ N_alloc ];
    assert_msg( active != 0, "ERROR in TangentialInterpolation::Allocate: could not allocate memory for 'active'\n");

    m_X     = alloc_array( N_alloc*D );
    x       = alloc_array( D );
//     sin_xi  = alloc_array( N_alloc );
    theta   = alloc_array( N_alloc );
    xi      = alloc_array( N_alloc );
    zeta    = alloc_array( N_alloc );

    // make sure the working array has the appropriate size
    if( w_i != 0 ) { delete [] w_i; w_i = 0; }
    w_i     = new int [ 32*N_alloc*N_alloc ];
    m_K     = new double [ N_alloc*N_alloc ];
    m_Kf    = new double [ N_alloc*N_alloc ];
    if( sym ) {
        zeta_tilde  = alloc_array( N_alloc );
        zeta_star   = alloc_array( N_alloc );
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
    N_alloc = 0;
    N       = 0;
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
            const int       a_N,
            const int       a_D,
            const double * const * a_directions,
            const bool      a_sym,
            const double    a_gamma
        )
{
    Allocate( a_N, a_D );
    for( int i_dir=0; i_dir<a_N; i_dir++)
        AddDirection( a_directions[i_dir] );
    SetSymmetric(a_sym);
    if( a_gamma>0 )
    {
        SetGamma( a_gamma );
        ComputeKernelMatrix();
    }
}
/* *************************************************************************************** */
void TangentialInterpolation::AddDirection( const double * a_X )
{
    assert_msg( N < N_alloc, "ERROR in TangentialInterpolation::AddInterpolationData: allocation size too small (trying to add another entry although n=N_alloc)\n");
    init = false; // reset initialization flag since kernel matrix needs to be re-computed!

    const double l = norm( a_X, D );
    for(int d=0; d<D; d++) m_X[N*D+d] = a_X[d] / l; // make sure the length is 1 for the directions
    N++; // increment the counter for the dimension of m_X, i.e. the number of actually provided training directions a.k.a. support points for the spherical basis functions
}
/* *************************************************************************************** */
void TangentialInterpolation::Weights(
        double * a_W,
        const int n_dirs,
        const double * a_x,
        double * a_zeta,
        double * a_dzeta,
        double * a_ddzeta
                )
{
    fprintf(stderr, "----=== DO NOT USE THIS FUNCTION ===----\nunless you are sure that Kinv is not multiplied into RI\n");
    exit(-1);

    const bool b_zeta       = (a_zeta != 0);
    const bool b_dzeta      = (a_dzeta != 0);
    const bool b_ddzeta     = (a_ddzeta != 0);

    if( !init ) { ComputeKernelMatrix(); }

    double * rhs = alloc_array( n_dirs*N );
    int offset = 0;
    for( int i_dir=0; i_dir < n_dirs; i_dir++)
    {

        double radius = norm( a_x + i_dir*D, D );
        // in the case of a very small amplitude, interpolate value of the first pcw cubic function at radius=0.
        for(int d=0;d<D;d++) x[d] = a_x[d+i_dir*D]/radius; // normalize inputs, i.e. direction of input

        // compute theta = cos(xi)
        MatVecMul( m_X, x, theta, N, D );

        for(int n=0; n<N; n++) {
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
            for( int n=0; n<N; n++ ) {
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
            for(int n=0; n<N; n++)
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

        if( b_zeta )
        {
#pragma unroll (4)
            for(int n=0; n<N; n++)
                a_zeta[n+offset]=zeta[n];
        }

        if( sym )
        {
#pragma unroll (4)
            for(int n=0; n<N; n++)
                rhs[n_dirs*n+i_dir]=zeta_star[n];
        }
        else
        {
#pragma unroll (4)
            for(int n=0; n<N; n++)
                rhs[n_dirs*n+i_dir]=zeta[n];
        }
        offset += N;
    }

//     SolveByFactorization( m_Kf, rhs, a_W, w_i, N, 1 );
    SolveByFactorization( m_Kf, rhs, a_W, w_i, N, n_dirs );

    free_array(&rhs);
}/* *************************************************************************************** */
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

    if( !init ) { ComputeKernelMatrix(); } // TODO this is unnecessary. make code consistent

    int offset = 0;
    for( int i_dir=0; i_dir < n_dirs_eval; i_dir++)
    {

        double radius = norm( a_x + i_dir*D, D );
        if( radius>1e-8 )
            for(int d=0;d<D;d++) x[d] = a_x[d+i_dir*D]/radius; // normalize inputs, i.e. direction of input
        // ATTENTION in the case of a very small amplitude, change input direction to [1,0,...0] and radius to 1.
        else
        {
            radius  = 1.;
            x[0]    = 1.;
            for(int d=1;d<D;d++) x[d] = 0.;
        }

        // compute theta = cos(xi)
        MatVecMul( m_X, x, theta, N, D );

        for(int n=0; n<N; n++) {
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
            for( int n=0; n<N; n++ ) {
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
            for(int n=0; n<N; n++)
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
        for(int n=0; n<N; n++)
            a_zeta[n+offset]=zeta[n];

        offset += N;
    } // for( int i_dir=0; i_dir < n_dirs_eval; i_dir++)
}
// /* *************************************************************************************** */
// void TangentialInterpolation::Weights(
//         int    n_vec,       /*!< [in]  number of input vector */
//         double * a_W,       /*!< [out] (n x n_vec) matrix \c W of weights */
//         const double * a_x, /*!< [in]  (n_vec x D) matrix containing directions as /direction \c X */
//         double * w_d        /* sufficiently large working array of doubles */)
// {
//
//     assert_msg( 1==0, "ERROR: Implementation of TangentialInterpolationPseudoSym::Weights() with multiple r.h.s. is yet outstanding\n");
//
//     if( !init ) { ComputeKernelMatrix(); }
//
// //     double * X = w_d;
// //     double * THETA = w_d + n_vec*D;
// //     double * XI = THETA + n_vec;
// //     double * ZETA = XI + n_vec;
// //     for( int i_dir = 0; i_dir < n_vec; i_dir++ )
// //     {
// //         double radius[i_dir] = norm( a_x, D );
// //         for(int d=0;d<D;d++) X[i_dir*D+d] = a_x[i_dir*D+d]/radius; // normalize inputs, i.e. direction of input
// //     }
// //     // compute theta = cos(xi)
// //     MatMatTMul( m_X, X, N, D, n_vec );
// // #pragma unroll 4
// //     for(int n=0; n<N*n_vec; n++) {
// //         XI[n] = safeAcos(THETA[n]);
// //         ZETA[n] = exp( - gamma * XI[n] * XI[n] );
// //         if( sym )
// //             ZETA[n] += exp( - gamma * (PI-XI[n]) * (PI-XI[n]) );
// //     }
// //
// //     if( sym )
// //         SolveByFactorization( m_Kf, zeta_star, a_W, w_i, N, n_vec );
// //     else
// //         SolveByFactorization( m_Kf, zeta, a_W, w_i, N, n_vec );
// }
/* *************************************************************************************** */
void TangentialInterpolation::SetSymmetric( const bool a_sym )
{
    sym  = a_sym;
    printf("$   set symmetry flag to %i\n", sym);
    init = false; // ComputeKernelMatrix() will be called (again)
}
/* *************************************************************************************** */
void TangentialInterpolation::SetGamma( const double a_gamma )
{
    assert_msg( fabs(a_gamma)>1e-8, "ERROR in TangentialInterpolation::SetGamma: kernel width parameter must be significantly positive\n");
    gamma = a_gamma;
    printf("$   set gamma to %lf\n", gamma);
    init  = false; // ComputeKernelMatrix() will be  called (again)
}
/* *************************************************************************************** */
void TangentialInterpolation::SetLambda( const double a_lambda )
{
    assert_msg( a_lambda >= 0., "ERROR in TangentialInterpolation::SetLambda: regression parameter must be non-negative\n");
    lambda = a_lambda;
    printf("$   set lambda to %lf\n", lambda);
    init   = false; // ComputeKernelMatrix() will be  called (again)

}
/* *************************************************************************************** */
void TangentialInterpolation::ComputeKernelMatrix()
{
    assert_msg( gamma>1e-3, "ERROR in TangentialInterpolation::ComputeKernelMatrix: gamma must be significantly greater than zero!\n");
    double tmp_alpha = 0.;
    for(int n1=0;n1<N;n1++)
    {
        for(int n2=n1;n2<N;n2++)
        {
            tmp_alpha  = safeAcos( VecVecMul( m_X + n1*D, m_X + n2*D, D ) );
            if ( sym )
                    m_K[n1*N+n2]  = exp(-gamma*tmp_alpha*tmp_alpha)  +  exp(-gamma*(PI-tmp_alpha)*(PI-tmp_alpha));
            else
                    m_K[n1*N+n2]  = exp(-gamma*tmp_alpha*tmp_alpha);
            m_K[n2*N+n1] = m_K[n1*N+n2];
        }
        m_K[n1*(N+1)] += lambda; // add regression parameter
    }

    // store LDL factorization of kernel matrix
    Factorize( m_K, m_Kf, w_i, N );

    // ready to go
    init = true;

}
/* *************************************************************************************** */
#if 0
TangentialInterpolationPseudoSym::TangentialInterpolationPseudoSym()
{
    zero_pointers();

}

void TangentialInterpolationPseudoSym::zero_pointers()
{
    TangentialInterpolation::zero_pointers();
    m_Kdiff     = 0;
    m_Kdiff_f   = 0;
    w_s         = 0;
    w_diff_i    = 0;
    sym         = true; // enforeces allocation of zeta_star and zeta_tilde!
}

TangentialInterpolationPseudoSym::~TangentialInterpolationPseudoSym()
{
    Free();

}
void TangentialInterpolationPseudoSym::Free()
{
    TangentialInterpolation::Free();
    free_array( &m_Kdiff );
    free_array( &m_Kdiff_f );
    free_array( &w_s );
    delete [] w_diff_i;  w_diff_i = 0;

}

void TangentialInterpolationPseudoSym::Allocate( int a_N_alloc, const int a_D )
{
    Free();
    TangentialInterpolation::Allocate( a_N_alloc, a_D );
    m_Kdiff     = alloc_array( a_N_alloc*a_N_alloc );
    m_Kdiff_f   = alloc_array( a_N_alloc*a_N_alloc );
    w_diff_i    = new int [ 32*a_N_alloc*a_N_alloc];
    w_s         = alloc_array( a_N_alloc );

}

void TangentialInterpolationPseudoSym::RecomputeKernelMatrix()
{
//     printf("\n\n\nKERNEL MATRIX\n\n\n\nDERIVED\n");
//     fflush(stdout);
    double tmp_alpha = 0., A=0., B=0.;
    for(int n1=0;n1<N;n1++)
    {
        for(int n2=n1;n2<N;n2++)
        {
            tmp_alpha  = safeAcos( VecVecMul( m_X + n1*D, m_X + n2*D, D ) );
            A = exp(-gamma*tmp_alpha*tmp_alpha);
            B = exp(-gamma*(PI-tmp_alpha)*(PI-tmp_alpha));
            m_K[n1*N+n2]        = A + B;
            m_Kdiff[n1*N+n2]    = A - B;
            m_K[n2*N+n1]        = m_K[n1*N+n2];
            m_Kdiff[n2*N+n1]    = m_Kdiff[n1*N+n2];
        }
        m_K[n1*(N+1)]       += lambda; // add regression parameter
        m_Kdiff[n1*(N+1)]   += lambda; // add regression parameter
    }

    // store LDL factorization of kernel matrix
    Factorize( m_K, m_Kf, w_i, N );
    Factorize( m_Kdiff, m_Kdiff_f, w_diff_i, N );

}

void TangentialInterpolationPseudoSym::Weights(
        double * a_W,
        const int n_dirs,
        const double * a_x,
        double * a_zeta,
        double * a_dzeta,
        double * a_ddzeta
                )
{
//     printf("INIT? %i \n", init);

    const bool b_zeta       = (a_zeta != 0);
    const bool b_dzeta      = (a_dzeta != 0);
    const bool b_ddzeta     = (a_ddzeta != 0);

    if( !init ) { ComputeKernelMatrix(); }

    double * rhs_a = alloc_array( n_dirs*N ), * rhs_b = alloc_array( n_dirs*N ), * weights = alloc_array(n_dirs*2*N);
//     printf("# initialized\n");
//     fflush(stdout);
    int offset = 0;
    for( int i_dir = 0; i_dir < n_dirs; i_dir++ )
    {

        double radius = norm( a_x + i_dir*D, D );
        // in the case of a very small amplitude, interpolate value of the first pcw cubic function at radius=0.
        for(int d=0;d<D;d++) x[d] = a_x[d+i_dir*D]/radius; // normalize inputs, i.e. direction of input

        // compute theta = cos(xi)
        MatVecMul( m_X, x, theta, N, D );

        for(int n=0; n<N; n++) {
            xi[n] = safeAcos(theta[n]);
            if( theta[n] > theta_max ) {
                active[n]   = false;
                sin_xi[n]   = 1e-16;
            }
            else if ( theta[n] < -theta_max ) {
                active[n]   = false;
                sin_xi[n]   = 1e-16;
            }
            else {
                active[n]   = true;
                sin_xi[n]   = sqrt( 1.-theta[n]*theta[n]); // sin(xi[n]);
            }
        }
    //     printf("# xi ok\n");
    //     fflush(stdout);
        // symmetric case:
        // the intermediate vectors zeta_tilde and zeta_star are re-used in the gradient computation
    #pragma unroll (4)
        for( int n=0; n<N; n++ ) {
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
#pragma unroll (4)
        if( b_zeta )
            for(int n=0; n<N; n++)
                a_zeta[n+offset]=zeta[n];


#pragma unroll (4)
        for(int n=0; n<N; n++)
            zeta[n]-=zeta_tilde[n];
#pragma unroll (4)
        for(int n=0; n<N; n++)
            rhs_a[n*n_dirs+i_dir] = zeta[n];
#pragma unroll (4)
        for(int n=0; n<N; n++)
            rhs_b[n*n_dirs+i_dir] = zeta_star[n];
        offset += N;
    }
//     printf("# zeta, zeta_tilde, zeta_star ok\n");
//     fflush(stdout);


    SolveByFactorization( m_Kf, rhs_a, weights, w_i, N, n_dirs );
    SolveByFactorization( m_Kdiff_f, rhs_b, weights+n_dirs*N, w_i, N, n_dirs );
    for(int i_dir=0;i_dir<n_dirs;i_dir++)
    {
#pragma unroll (4)
        for(int n=0; n<N; n++)
        {
            a_W[i_dir*2*N+n] = 0.5*(weights[n*n_dirs+i_dir]+weights[(N+n)*n_dirs+i_dir]);
            a_W[i_dir*2*N+N+n] = 0.5*(weights[n*n_dirs+i_dir]-weights[(N+n)*n_dirs+i_dir]);
        }
    }
    free_array(&weights);
    free_array(&rhs_a);
    free_array(&rhs_b);

}
// /* *************************************************************************************** */
// void TangentialInterpolationPseudoSym::Weights(
//         int    n_vec,       /*!< [in]  number of input vector */
//         double * a_W,       /*!< [out] (n x n_vec) matrix \c W of weights */
//         const double * a_x, /*!< [in]  (n_vec x D) matrix containing directions as /direction \c X */
//         double * w_d        /* sufficiently large working array of doubles */)
// {
//     assert_msg( 1==0, "ERROR: Implementation of TangentialInterpolationPseudoSym::Weights() with multiple r.h.s. is yet outstanding\n");
//
//     if( !init ) { ComputeKernelMatrix(); }
//
// //     double * X = w_d;
// //     double * THETA = w_d + n_vec*D;
// //     double * XI = THETA + n_vec;
// //     double * ZETA = XI + n_vec;
// //     for( int i_dir = 0; i_dir < n_vec; i_dir++ )
// //     {
// //         double radius[i_dir] = norm( a_x, D );
// //         for(int d=0;d<D;d++) X[i_dir*D+d] = a_x[i_dir*D+d]/radius; // normalize inputs, i.e. direction of input
// //     }
// //     // compute theta = cos(xi)
// //     MatMatTMul( m_X, X, N, D, n_vec );
// // #pragma unroll 4
// //     for(int n=0; n<N*n_vec; n++) {
// //         XI[n] = safeAcos(THETA[n]);
// //         ZETA[n] = exp( - gamma * XI[n] * XI[n] );
// //         if( sym )
// //             ZETA[n] += exp( - gamma * (PI-XI[n]) * (PI-XI[n]) );
// //     }
// //
// //     if( sym )
// //         SolveByFactorization( m_Kf, zeta_star, a_W, w_i, N, n_vec );
// //     else
// //         SolveByFactorization( m_Kf, zeta, a_W, w_i, N, n_vec );
// }
#endif
