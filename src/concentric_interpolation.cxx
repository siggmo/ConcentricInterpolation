/* *************************************************************************************** */
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

#include <concentric_interpolation.h>
using namespace UTILITY;

#define TESTBDRYDATA
/* *************************************************************************************** */
const double ConcentricInterpolation::small = 1.e-12;    //!< small number; used, e.g., in order to prevent division by zero errors
const double ConcentricInterpolation::theta_max = 0.99985 /*1.000000*/; //!< required in order to regularize the derivative of zeta!
/* *************************************************************************************** */
ConcentricInterpolation::ConcentricInterpolation( const bool a_sym )
{
    FastAcos::Init();
    sym         = a_sym; // symmetry flag
    if (sym) printf("$   symmetry flag is TRUE\n");
    init        = false; // set initialization flag to FALSE
    N           = 0;
    N_alloc     = 0;
    lambda      = 0.;
    D       = 0;
    gamma       = 0;
//  s_type      = 1; // assume: quadratic functions to the left and to the right
    zero_pointers();
}
/* *************************************************************************************** */
void ConcentricInterpolation::zero_pointers()
{
    // initialize pointers to NULL
    
    // private members:
    S           = 0;
    w_i         = 0;
    m_K         = 0;
    m_Kf        = 0;
    m_X         = 0;
    m_J0        = 0;
    v           = 0;
    dv          = 0;
    ddv         = 0;
    Ki_zeta     = 0;
    Ki_v        = 0;
    Ki_dv       = 0;
    theta       = 0;
    zeta        = 0;
    dzeta       = 0;
    ddzeta      = 0;
    zeta_star   = 0;
    zeta_tilde  = 0;
    active      = 0;
    sin_xi      = 0;
    xi          = 0;
    x           = 0;
}
/* *************************************************************************************** */
void ConcentricInterpolation::Allocate( int a_N_alloc, const int a_D )
{
    assert_msg( a_N_alloc>0, "ERROR in ConcentricInterpolation::Allocate: N_alloc > 0 required\n");
    assert_msg( a_D>0, "ERROR in ConcentricInterpolation::Allocate: D > 0 required\n");
    Free();

    N_alloc = a_N_alloc;
    D   = a_D;
    N       = 0;
    S       = new CubicInterpolant [N_alloc];

    active  = new bool [ N_alloc ];
    assert_msg( active != 0, "ERROR in ConcentricInterpolation::Allocate: could not allocate memory for 'active'\n");
    // m_K and m_Kf are only initialized when actually needed!
    m_X     = alloc_array( N_alloc*D );
    m_J0    = alloc_array( D*D );
    x       = alloc_array( D );
    sin_xi  = alloc_array( N_alloc );
    v       = alloc_array( N_alloc );
    dv      = alloc_array( N_alloc );
    ddv     = alloc_array( N_alloc );
    Ki_zeta = alloc_array( N_alloc );
    Ki_v    = alloc_array( N_alloc );
    theta   = alloc_array( N_alloc );
    xi      = alloc_array( N_alloc );
    zeta    = alloc_array( N_alloc );
    dzeta   = alloc_array( N_alloc );
    ddzeta  = alloc_array( N_alloc );
    Ki_dv   = alloc_array( N_alloc );
    if( sym ) {
        zeta_tilde  = alloc_array( N_alloc );
        zeta_star   = alloc_array( N_alloc );
    }
}
/* *************************************************************************************** */
ConcentricInterpolation::~ConcentricInterpolation() {
    Free();
}
/* *************************************************************************************** */
void ConcentricInterpolation::Free()
{
    if( S != 0 ) delete [] S;
    delete [] w_i; w_i = 0;
    init    = false;
    S       = 0;
    N_alloc = 0;
    N       = 0;
    free_array( &m_X );
    free_array( &m_K );
    free_array( &m_Kf );
    free_array( &m_J0 );
    free_array( &x );
    if( active != 0 )
    {
        delete [] active;
        active = 0;
    }
    free_array( &sin_xi );
    free_array( &v );
    free_array( &dv );
    free_array( &ddv );
    free_array( &Ki_zeta );
    free_array( &Ki_v );
    free_array( &Ki_dv );
    free_array( &theta );
    free_array( &zeta );
    free_array( &dzeta );
    free_array( &ddzeta );
    free_array( &zeta_tilde );
    free_array( &zeta_star );
    free_array( &xi );
}
/* *************************************************************************************** */
void    ConcentricInterpolation::GetTrainingRadii( double * a_training_radii, int & a_num_training_radii ) const
{
    S[0].GetSupport( a_training_radii, a_num_training_radii );
}
/* *************************************************************************************** */
void ConcentricInterpolation::AddInterpolationDataDF( const double * a_radii, const double * a_f, const double * a_df, int a_n, const double * a_X )
{
    assert_msg( N < N_alloc, "ERROR in ConcentricInterpolation::AddInterpolationData: allocation size too small (trying to add another entry although n=N_alloc)\n");
    init = false; // reset initialization flag since kernel matrix needs to be re-computed!
        
    S[N].SetData( a_n, a_radii, a_f, a_df );
    
    const double l = norm( a_X, D );
    for(int d=0; d<D; d++) m_X[N*D+d] = a_X[d] / l; // make sure the length is 1 for the directions
    N++; // increment the counter for the dimension of m_X, i.e. the number of training directions
}
/* *************************************************************************************** */
double ConcentricInterpolation::Interpolate( const double * a_x )
{
    if( !init ) { InitializeKernelMethod(); }
    radius = norm( a_x, D );
    // in the case of a very small amplitude, interpolate value of the first pcw cubic function at radius=0.
    if( radius < small ) {

        for(int n=0; n<N; n++) {
            active[n]   = false;
            xi[n]       = 0.;
            sin_xi[n]   = 1e-16;
        }
        
        double f, df, ddf;
        f=S[0].Interpolate( 0., df, ddf);
        SetVector( dv, N, df );
        SetVector( ddv, N, ddf );
        return f;
    }
    RadialInterpolation( radius ); // updates v, dv, ddv
    
    SolveByFactorization(m_Kf, v, Ki_v, w_i, N);
    
    
    for(int d=0;d<D;d++) x[d] = a_x[d]/radius; // normalize inputs, i.e. direction of input

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
            if( active[n] )
            {
                dzeta[n]    = - 2.*gamma*( zeta[n]*xi[n] - (PI-xi[n])*zeta_tilde[n] );
                ddzeta[n]   = -2.*gamma * (
                        zeta[n]         * ( 1. - 2.*gamma*xi[n]*xi[n]  )
                    +   zeta_tilde[n]   * ( 1. - 2.*gamma*(PI-xi[n])*(PI-xi[n]) ) 
                                        );
            }
            else
            {
                dzeta[n]    = 0.;
                ddzeta[n]   = 0.;
            }
        }
        return VecVecMul( Ki_v, zeta_star, N );
    }
    
    // for the non-symmetric case:
#pragma unroll (4)
    for(int n=0; n<N; n++)
    {
        zeta[n] = exp(- gamma * xi[n]*xi[n] );
        if( active[n] )
        {
            dzeta[n]    = - 2.*gamma*zeta[n]*xi[n];
            ddzeta[n]   = -2.*gamma * zeta[n] * ( 1. - 2.*gamma*xi[n]*xi[n]  );
        }
        else
        {
            dzeta[n]    = 0.;
            ddzeta[n]   = 0.;
        }
    }
    return VecVecMul( Ki_v, zeta, N );
}
/* *************************************************************************************** */
void ConcentricInterpolation::Gradient( double * grad )
{
    // ASSERTION: Interpolate has been called before (re-uses zeta and Ki_v)
    
    /* note: instead of locally computing the directions Q_i pointing from the query 
     * direction N to the training directions N_n, and adding them with weights tau
     * to the radial stress, we
     *    1) compute the Q_n without the subtraction of N, and sum over n
     *    2) substract the part of the part to be substracted (i.e. the part of the N_i that is parallel to N)
     *       using only the theta_i --> only one vector summation needed
     */
    
    
    SetVector( grad, D, 0. );
    // if the input has negligible norm, then return 0 (gradient not defined at 0!)
    if( radius < small ) return;

    // FF: adjust equation number
    // see paper (65):
    // sigma = x * (s'^T * Ki * v) - 1/radius ( s^T * Ki )_j dzeta_j Q_j
    
    // recall x from Interpolate(): normalized input, i.e. direction of input to Interpolate()
    double  factor;

    // *********************************************************
    // PART 1: contribution of the tangential derivatives
    for(int n1=0; n1<N; n1++)
    {
        if( active[n1] )
        {
            factor      = -Ki_v[n1]*dzeta[n1]/sin_xi[n1];              // tau_n1 * radius / sin(xi_n1)
            for(int d=0; d<D; d++) grad[d] += factor * m_X[n1*D+d];
        }
    }
    for(int d=0; d<D; d++) grad[d] *= 1./radius;
    // now: grad = - 1/radius  sum_n2,n1  S_n2 Kinv_n2,n1 dzeta_n1  N_n1  / sin(xi_n1)
    //           = sum_n  tau_n / sin(xi_n) * N_n
    //      i.e. the non-theta part of the Q-part of the gradient
    // *********************************************************
    // PART 2: contribution due to the radial derivative
    if( sym )
        SolveByFactorization(m_Kf, zeta_star, Ki_zeta, w_i, N);
    else
        SolveByFactorization(m_Kf, zeta, Ki_zeta, w_i, N);
    factor = VecVecMul( Ki_zeta, dv, N );                                // sigma_0
    // *********************************************************
    // PART 3: add the radial derivative and subtract the part of Q that was not considered in PART 1
    const double proj=VecVecMul( grad, x, D ); // = projection onto the normal space to x, i.e. onto tangential plane
    for(int d=0; d<D; d++) grad[d] += (factor-proj)*x[d];
}
/* *************************************************************************************** */
void ConcentricInterpolation::AutodefineJ0_psi()
{
    AutodefineJ0( 1., 0.);
}
/* *************************************************************************************** */
void ConcentricInterpolation::AutodefineJ0_dpsi()
{
    AutodefineJ0( 0., 1.);
}
/* *************************************************************************************** */
void ConcentricInterpolation::GetJ0( double * J0 ) const
{
    for(int dd=0; dd<D*D; dd++) J0[dd] = m_J0[dd];
}
/* *************************************************************************************** */
void ConcentricInterpolation::AutodefineJ0( const double w_psi, const double w_dpsi )
{
    double  w1      = w_psi/(w_psi+w_dpsi),
            w2      = w_dpsi/(w_psi+w_dpsi);        // normalize weights
    int     D_sym   = (D*(D+1))/2;                  // dimension of symmetric matrices
    double  * m_H   = alloc_array( D_sym * D_sym ); // the Hessian needed for the fitting procedure
    double  * v_rhs = alloc_array( D_sym );         // residual
    double  n_tmp[D_sym];
    const double sqrt2=sqrt(2.);
    // initialize to zero
    for( int d_s=0; d_s<D_sym; d_s++ )          v_rhs[d_s]  = 0.;
    for( int dd_s=0; dd_s<D_sym*D_sym; dd_s++ ) m_H[dd_s]   = 0.;
    
    // assemble Hessian and right hand side
    double v0, dv0, x0;
    for(int n=0; n<N; n++)
    {
        // cast X[n] times X[n] into a vector using the symmetry convention
        int ct=0;
        for( int d=0; d<D; d++)
        {
            n_tmp[ct++] = m_X[n*D+d]*m_X[n*D+d];
            for(int d2=d+1; d2<D; d2++)
                n_tmp[ct++] = sqrt2 * m_X[n*D+d]*m_X[n*D+d2];
        }
        S[n].FirstV( v0, x0 );
        v0 -= S[n].Value(0.); // substract the value at zero --> shift function
        S[n].FirstDV( dv0, x0 );
        const double f = w1 * ( 0.5*v0*x0*x0 ) + w2 * ( dv0*x0 ); // factor for RHS vector
        const double g = w1 * ( 0.25*x0*x0*x0*x0 ) + w2 * ( x0*x0 ); // factor for Hessian
        
        for( int d=0; d<D_sym; d++ ) v_rhs[d] += f*n_tmp[d];
        for( int d=0; d<D_sym; d++ )
        {
            for( int d2=d; d2<D_sym; d2++ )
            {
                m_H[d*D_sym+d2] += g*n_tmp[d]*n_tmp[d2];
            }
        }
    }
    // do not use w_i as it may not yet have been allocated...
    int *i_tmp = new int [ 32*N_alloc*N_alloc ];
    int err = LAPACKE_dsysv( LAPACK_ROW_MAJOR, 'U', D_sym, 1, m_H, D_sym, i_tmp, v_rhs, 1 ); //ATTENTION what to do with this?

    int ct=0;
    for( int d=0;d<D; d++ )
    {
        m_J0[d*D+d] = v_rhs[ct++];
        for( int d2=d+1;d2<D;d2++)
            m_J0[d*D+d2] = m_J0[d2*D+d] = v_rhs[ct++] / sqrt2;
    }
    free_array( &m_H );
    free_array( &v_rhs );
}
/* *************************************************************************************** */
void ConcentricInterpolation::SetJ0( const double * const J0 )
{
    for(int dd=0; dd<D*D; dd++) m_J0[dd] = J0[dd];
}
/* *************************************************************************************** */
void ConcentricInterpolation::Hessian( double * hess )
{
    // if small radius, put dummy Hessian here
    if( radius < small ) {
        for(int dd=0; dd<D*D; dd++)
            hess[dd]=m_J0[dd];
        return;
    }
    
    // initialize quantities to zero
    double  beta=0., alpha=0., dalpha = 0., omega = 0.;
    double Q[D], Qbar[D];

    for(int d=0; d<D; d++)      Qbar[d]     = 0.;
    for(int dd=0; dd<D*D; dd++) hess[dd]   = 0.;

    // compute Ki_dv, i.e. product of inverse kernel matrix and vector of radial derivatives
    SolveByFactorization(m_Kf, dv, Ki_dv, w_i, N);
    
    // set up the upper half of the Hessian
    for(int n=0; n<N; n++)
    {
        if( active[n] )
        {
            const double factor = 1./sin_xi[n];
            for( int d=0; d<D; d++ ) Q[d] = factor * ( m_X[n*D+d] - theta[n]*x[d] );
            
            omega   = Ki_v[n] * dzeta[n]/(radius*radius);
            dalpha  = omega*theta[n]*factor;
            alpha   += dalpha;

            beta    = - dalpha + Ki_v[n] * ddzeta[n]/(radius*radius);
            omega   += - Ki_dv[n] * dzeta[n]/radius;
            // add beta Q dyad Q and increment Qbar <-- +omega Q
            for(int d=0; d<D; d++)
                for(int d2=d; d2<D; d2++)
                    hess[d*D+d2] += beta * Q[d]*Q[d2];

            for(int d=0; d<D; d++) Qbar[d] += omega*Q[d];
        }
        else {
            // this is the only contribution that yields a nonzero part of the Hessian
            // in the case of xi_i  approx. 0
            alpha -= 2.*gamma/(radius*radius) * Ki_v[n];
        }
    }
    alpha += VecVecMul( dv, Ki_zeta, N ) / radius;
    const double mu_bar = VecVecMul( ddv, Ki_zeta, N );
    for(int d=0; d<D; d++)
    {
        hess[d*(D+1)] += alpha;
        for(int d2=d; d2<D; d2++)
        {
            hess[d*D+d2] +=     ( x[d] * Qbar[d2] + Qbar[d]*x[d2] )
                                +   (mu_bar - alpha) * x[d]*x[d2];
        }
    }
    
    // symmetrize
    for(int d=0; d<D; d++)
        for(int d2=0; d2<d; d2++)
            hess[d*D+d2] = hess[d2*D+d];
}
/* *************************************************************************************** */
void ConcentricInterpolation::HessianFDM( double * hess, const double dx )
{
    assert_msg( dx > 0, "ERROR in ConcentricInterpolation::HessianFDM: perturbation parameter must be positiv\n");
    
    for(int dd=0;dd<D*D; dd++) hess[dd] = 0.;
    // compute Hessian via finite difference method
    double x_new[D], x0[D], g0[D];
    Gradient(g0);
    for( int d=0; d<D; d++) x0[d] = radius*x[d];
    for( int i_pert=0; i_pert<D; i_pert++)
    {
        for( int d=0; d<D; d++) x_new[d] = x0[d];
        x_new[i_pert] += dx;
        Interpolate(x_new);
        Gradient(hess + i_pert*D);
        for( int d=0; d<D; d++) hess[i_pert*D+d] -= g0[d];
    }
    for(int d=0;d<D;d++)
        for(int d2=d;d2<D;d2++)
        {
            double sym_val = 0.5*( hess[ IDX2(d,d2,D)  ] + hess[ IDX2(d2,d,D) ] );
            hess[ IDX2(d,d2,D)  ] = hess[ IDX2(d2,d,D) ] = sym_val/dx;
        }
}
/* *************************************************************************************** */
void ConcentricInterpolation::InitializeKernelMethod( )
{
    if( init ) {
        free_array( &m_K );
        free_array( &m_Kf );
        delete [] w_i;      w_i = 0;
    }

    m_K = alloc_array( N*N );
    m_Kf = alloc_array( N*N );
    // compute the kernel matrix
    double tmp_alpha = 0.;
    for(int n=0;n<N;n++)
    {
        for(int n2=n;n2<N;n2++)
        {
            tmp_alpha  = safeAcos( VecVecMul( m_X + n*D, m_X + n2*D, D ) );
            if ( sym )
                    m_K[n*N+n2]  = exp(-gamma*tmp_alpha*tmp_alpha)  +  exp(-gamma*(PI-tmp_alpha)*(PI-tmp_alpha));
            else  
                    m_K[n*N+n2]  = exp(-gamma*tmp_alpha*tmp_alpha);
            m_K[n2*N+n] = m_K[n*N+n2];
        }
        m_K[n*(N+1)] += lambda;
    }

    // make sure the working array has the appropriate size
    if( w_i != 0 ) { delete [] w_i; w_i = 0; }
    w_i = new int [ 32*N_alloc*N_alloc ];
    Factorize( m_K, m_Kf, w_i, N );
    init    = true;
}
/* *************************************************************************************** */
void ConcentricInterpolation::RadialInterpolation( const double a_radius )
{
    // apply pcw cubic function interpolation (and compute also the first and second derivatives at all nodes)
    radius = a_radius;
    for( int n=0; n<N; n++ )    v[n]    = S[n].Interpolate(radius, dv[n], ddv[n]);
}
/* *************************************************************************************** */
bool ConcentricInterpolation::CheckRadius( const double * a_radii, const int a_R ) const
{
    assert_msg( ( a_R >= 2 ), "ERROR in ConcentricInterpolation::CheckRadius: less than 2 data points provided\n");
    bool sc = fabs( a_radii[0] ) < small; // is the first vlaue 0?
    for( int r=0; r<a_R; r++) {
        assert_msg( a_radii[r] > -small, "ERROR in ConcentricInterpolation::CheckRadius: a_radii[r] < 0 detected\n");
        if(r>0) assert_msg( a_radii[r] > a_radii[r-1], "ERROR in ConcentricInterpolation::CheckRadius: a_radii[r+1] <= a_radii[r] detected\n");
    }
    return sc;
}
/* *************************************************************************************** */
void ConcentricInterpolation::SetLambda( const double a_lambda )
{
    assert_msg( a_lambda >= 0., "ERROR in ConcentricInterpolation::SetLambda: regression parameter must be non-negative\n");
    lambda = a_lambda;
    printf("$   set lambda to %lf\n", lambda);

}
/* *************************************************************************************** */
void ConcentricInterpolation::RecomputeKernelMatrix()
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
