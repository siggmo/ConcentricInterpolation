/* *************************************************************************************** */
/*
 *  TangentialInterpolation
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

#include <tangential_interpolation.h>
using namespace UTILITY;

/* *************************************************************************************** */
const double TangentialInterpolation::small = 1.e-12;    //!< small number; used, e.g., in order to prevent division by zero errors
const double TangentialInterpolation::theta_max = 0.99995 /*1.000000*/; //!< required in order to regularize the derivative of zeta!
/* *************************************************************************************** */
TangentialInterpolation::TangentialInterpolation( const bool a_sym )
{
    FastAcos::Init();
    sym         = a_sym; // symmetry flag
    if (sym) printf("$   symmetry flag is TRUE\n");
    init        = false; // set initialization flag to FALSE
    N           = 0;
    N_alloc     = 0;
    lambda      = 0.;
    D           = 0;
    gamma       = 0;
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
    sin_xi      = 0;
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
    sin_xi  = alloc_array( N_alloc );
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
    free_array( &sin_xi );
    free_array( &theta );
    free_array( &zeta );
    free_array( &zeta_tilde );
    free_array( &zeta_star );
    free_array( &xi );
}
/* *************************************************************************************** */
void TangentialInterpolation::AddDirection( const double * a_X )
{
    assert_msg( N < N_alloc, "ERROR in TangentialInterpolation::AddInterpolationData: allocation size too small (trying to add another entry although n=N_alloc)\n");
    init = false; // reset initialization flag since kernel matrix needs to be re-computed!
        
    const double l = norm( a_X, D );
    for(int d=0; d<D; d++) m_X[N*D+d] = a_X[d] / l; // make sure the length is 1 for the directions
    N++; // increment the counter for the dimension of m_X, i.e. the number of training directions
}
/* *************************************************************************************** */
void TangentialInterpolation::Weights(
        double * o_W,       /* [out] vector \c W of weights */
        const int n_dirs,   /* [in]  number of directions in \c X */
        const double * a_x, /* [in]  vector/direction \c X */
        double * o_zeta,    /* [out] vector \c W of zeta values (if not needed set to NULL) */
        double * o_dzeta,   /* [out] vector \c W of dzeta values (if not needed set to NULL) */
        double * o_ddzeta   /* [out] vector \c W of dzeta values (if not needed set to NULL) */
                )
{
    const bool b_zeta       = (o_zeta != 0);
    const bool b_dzeta      = (o_dzeta != 0);
    const bool b_ddzeta     = (o_ddzeta != 0);

    if( !init ) { InitializeKernelMethod(); }

    printf("n_dirs %i, N %i\n", n_dirs, N); fflush(stdout);
    double * rhs = alloc_array( n_dirs*N ), * weights = alloc_array(n_dirs*N);
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
                        if(b_dzeta)     o_dzeta[n+offset]    = - 2.*gamma*( zeta[n]*xi[n] - (PI-xi[n])*zeta_tilde[n] );
                        if(b_ddzeta)    o_ddzeta[n+offset]   = -2.*gamma * (
                                zeta[n]         * ( 1. - 2.*gamma*xi[n]*xi[n]  )
                            +   zeta_tilde[n]   * ( 1. - 2.*gamma*(PI-xi[n])*(PI-xi[n]) ) 
                                                );
                    }
                    else
                    {
                        if(b_dzeta)     o_dzeta[n+offset]    = 0.;
                        if(b_ddzeta)    o_ddzeta[n+offset]   = 0.;
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
                        if(b_dzeta)     o_dzeta[n+offset]    = - 2.*gamma*zeta[n]*xi[n];
                        if(b_ddzeta)    o_ddzeta[n+offset]   = -2.*gamma * zeta[n] * ( 1. - 2.*gamma*xi[n]*xi[n]  );
                    }
                    else
                    {
                        if(b_dzeta)     o_dzeta[n+offset]    = 0.;
                        if(b_ddzeta)    o_ddzeta[n+offset]   = 0.;
                    }
                }
            }
        }
        
        if( b_zeta ) 
        {
#pragma unroll (4)
            for(int n=0; n<N; n++)
                o_zeta[n+offset]=zeta[n];
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

//     SolveByFactorization( m_Kf, rhs, o_W, w_i, N, 1 );
    SolveByFactorization( m_Kf, rhs, o_W, w_i, N, n_dirs );
    
    free_array(&rhs);
    free_array(&weights);
}
// /* *************************************************************************************** */
// void TangentialInterpolation::Weights(
//         int    n_vec,       /*!< [in]  number of input vector */
//         double * o_W,       /*!< [out] (n x n_vec) matrix \c W of weights */
//         const double * a_x, /*!< [in]  (n_vec x D) matrix containing directions as /direction \c X */
//         double * w_d        /* sufficiently large working array of doubles */)
// {
// 
//     assert_msg( 1==0, "ERROR: Implementation of TangentialInterpolationPseudoSym::Weights() with multiple r.h.s. is yet outstanding\n");
// 
//     if( !init ) { InitializeKernelMethod(); }
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
// //         SolveByFactorization( m_Kf, zeta_star, o_W, w_i, N, n_vec );
// //     else
// //         SolveByFactorization( m_Kf, zeta, o_W, w_i, N, n_vec );
// }
/* *************************************************************************************** */
void TangentialInterpolation::InitializeKernelMethod( )
{
//     if( init ) {
//         free_array( &m_K );
//         free_array( &m_Kf );
//         delete [] w_i;      w_i = 0;
//     }
// 
//     m_K = alloc_array( N*N );
//     m_Kf = alloc_array( N*N );
    // compute the kernel matrix
    RecomputeKernelMatrix();
    Factorize( m_K, m_Kf, w_i, N );

    init    = true;
}
/* *************************************************************************************** */
void   TangentialInterpolation::SetGamma( const double a_gamma )
{
    gamma = a_gamma;
    if(init)    { RecomputeKernelMatrix(); }
    else        { InitializeKernelMethod();}
}
/* *************************************************************************************** */
void TangentialInterpolation::SetLambda( const double a_lambda )
{
    assert_msg( a_lambda >= 0., "ERROR in TangentialInterpolation::SetLambda: regression parameter must be non-negative\n");
    lambda = a_lambda;
    printf("$   set lambda to %lf\n", lambda);

}
/* *************************************************************************************** */
void TangentialInterpolation::RecomputeKernelMatrix()
{
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
    

}
/* *************************************************************************************** */
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
        double * o_W,       /* [out] vector \c W of weights */
        const int n_dirs,   /* [in] number of directions in a_x */
        const double * a_x, /* [in] vector/direction \c X */
        double * o_zeta,    /* [out] vector \c W of zeta values (if not needed set to NULL) */
        double * o_dzeta,   /* [out] vector \c W of dzeta values (if not needed set to NULL) */
        double * o_ddzeta   /* [out] vector \c W of dzeta values (if not needed set to NULL) */
                )
{
//     printf("INIT? %i \n", init);
    
    const bool b_zeta       = (o_zeta != 0);
    const bool b_dzeta      = (o_dzeta != 0);
    const bool b_ddzeta     = (o_ddzeta != 0);

    if( !init ) { InitializeKernelMethod(); }

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
                    if(b_dzeta)     o_dzeta[n+offset]    = - 2.*gamma*( zeta[n]*xi[n] - (PI-xi[n])*zeta_tilde[n] );
                    if(b_ddzeta)    o_ddzeta[n+offset]   = -2.*gamma * (
                            zeta[n]         * ( 1. - 2.*gamma*xi[n]*xi[n]  )
                        +   zeta_tilde[n]   * ( 1. - 2.*gamma*(PI-xi[n])*(PI-xi[n]) ) 
                                            );
                }
                else
                {
                    if(b_dzeta)     o_dzeta[n+offset]    = 0.;
                    if(b_ddzeta)    o_ddzeta[n+offset]   = 0.;
                }
            }
        }
#pragma unroll (4)
        if( b_zeta ) 
            for(int n=0; n<N; n++)
                o_zeta[n+offset]=zeta[n];
        
        
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
            o_W[i_dir*2*N+n] = 0.5*(weights[n*n_dirs+i_dir]+weights[(N+n)*n_dirs+i_dir]);
            o_W[i_dir*2*N+N+n] = 0.5*(weights[n*n_dirs+i_dir]-weights[(N+n)*n_dirs+i_dir]);
        }
    }
    free_array(&weights);
    free_array(&rhs_a);
    free_array(&rhs_b);
    
}
// /* *************************************************************************************** */
// void TangentialInterpolationPseudoSym::Weights(
//         int    n_vec,       /*!< [in]  number of input vector */
//         double * o_W,       /*!< [out] (n x n_vec) matrix \c W of weights */
//         const double * a_x, /*!< [in]  (n_vec x D) matrix containing directions as /direction \c X */
//         double * w_d        /* sufficiently large working array of doubles */)
// {
//     assert_msg( 1==0, "ERROR: Implementation of TangentialInterpolationPseudoSym::Weights() with multiple r.h.s. is yet outstanding\n");
// 
//     if( !init ) { InitializeKernelMethod(); }
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
// //         SolveByFactorization( m_Kf, zeta_star, o_W, w_i, N, n_vec );
// //     else
// //         SolveByFactorization( m_Kf, zeta, o_W, w_i, N, n_vec );
// }
