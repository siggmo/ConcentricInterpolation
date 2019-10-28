#include <quadratic_interpolant.h>
/*
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
 *  The latest version of this software can be obtained through https://github.com/EMMA-Group/ConcentricInterpolation
 *
 *
 */

QuadraticInterpolant::QuadraticInterpolant() {
    DefaultInit();
}

void QuadraticInterpolant::DefaultInit()
{
    init    = false;
    num_support_r       = 0;
    support_r       = 0;
    poly_a  = 0;
    poly_b  = 0;
    poly_c  = 0;
}
QuadraticInterpolant::~QuadraticInterpolant() {
    Free();
}
QuadraticInterpolant::QuadraticInterpolant( const int a_num_support_r, const double * a_support_r, const double * S )
{
    DefaultInit();
    assert_msg( a_num_support_r >= 3, "Error in QuadraticInterpolant: there must be more than 2 data samples for the construction of the interpolant");

    SetData( a_num_support_r, a_support_r, S );
}
void QuadraticInterpolant::Free()
{
    if(!init) return;
    free_array( &support_r );
    free_array( &poly_a );
    free_array( &poly_b );
    free_array( &poly_c );
    num_support_r = 0;
    init = false;
}
void QuadraticInterpolant::Allocate()
{
    if(init) Free();
    assert_msg( num_support_r >= 2, "Error in QuadraticInterpolant::Allocate: num_support_r must be >= 2\n");
    support_r = alloc_array( num_support_r );
    poly_a = alloc_array( num_support_r );
    poly_b = alloc_array( num_support_r );
    poly_c = alloc_array( num_support_r );
}
void QuadraticInterpolant::SetData( const int a_num_training_r, const double * a_training_r, const double * S )
{
    if( init ) Free();
    bool has_zero = ( fabs( a_training_r[0] ) < 1e-16);
    assert_msg( has_zero, "Error in QuadraticInterpolant::SetData: data must have zero.\n");
    assert_msg( a_num_training_r>=3, "Error in QuadraticInterpolant::SetData: number of training points must be >= 3.\n");
    assert_msg( a_num_training_r%2==1, "Error in QuadraticInterpolant::SetData: number of training points must be odd.\n");

    num_support_r = (a_num_training_r - 1)/2;
//     printf("num_support_r = %i\n", num_support_r);

    Allocate();

    // set support_r to a_training_r with initial 0
    for( int i_supp=0; i_supp<num_support_r; i_supp++) support_r[i_supp] = a_training_r[2*i_supp];
//     printf("support_r = ");
//     for(int r=0;r<num_support_r; r++)
//         printf("%8.3e ", support_r[r]);
//     printf("\n"), fflush(stdout);

    // set coefficients
    for( int i_supp=0; i_supp<num_support_r; i_supp++)
    {
        assert_msg( a_training_r[i_supp]<a_training_r[i_supp+1], "Error in QuadraticInterpolant: a_training_r[i_supp]<a_training_r[i_supp+1] is not guaranteed\n");
        SolveForQuadraticCoefficients(  0., a_training_r[2*i_supp+1]-a_training_r[2*i_supp+0], a_training_r[2*i_supp+2]-a_training_r[2*i_supp+0],
//         SolveForQuadraticCoefficients(  0., 0.25, 0.5,
                                        S[2*i_supp+0], S[2*i_supp+1], S[2*i_supp+2],
                                        poly_a[i_supp], poly_b[i_supp], poly_c[i_supp] );
//         printf("positions:%lf, %lf, %lf\n", 0., a_training_r[2*i_supp+1]-a_training_r[2*i_supp+0], a_training_r[2*i_supp+2]-a_training_r[2*i_supp+0]);
//         printf("supports: %lf, %lf, %lf\n", S[2*i_supp+0], S[2*i_supp+1], S[2*i_supp+2]);
//         printf("coeffs: a = %lf, b = %lf, c = %lf\n", poly_a[i_supp], poly_b[i_supp], poly_c[i_supp]), fflush(stdout);
    }

    // set initialization flag
    init = true;
}

void QuadraticInterpolant::FirstV( double & v0, double & r ) const
{
    assert_msg( num_support_r > 1, "Error in QuadraticInterpolant::FirstV: Polynomial not initialized (or too few data)");
    v0  = poly_a[1];
    r   = support_r[1];
}
void QuadraticInterpolant::FirstDV( double & dv0, double & r ) const
{
    // slope from the left
    dv0 = poly_b[0] + 2.*support_r[1]*poly_c[0];
    r   = support_r[1];
}
double QuadraticInterpolant::Value( const double r ) const {
    int idx = 0;
    bool sc = false;
    while( !sc && idx < num_support_r )
    {
        sc = (r < support_r[idx+1] ) && (r >= support_r[idx]);
        if(!sc) idx++;
    }
    if(!sc) idx = num_support_r-1; // extrapolate last quadratic function
    const double dx=r-support_r[idx];
    return poly_a[idx]+dx*(poly_b[idx]+dx*poly_c[idx]);
}

double QuadraticInterpolant::Interpolate( const double r, double & dS, double & ddS ) const
{
    int idx = 0;
    bool sc = false;
    while( !sc && idx < num_support_r )
    {
        sc = (r < support_r[idx+1] ) && (r >= support_r[idx]);
        if(!sc) idx++;
    }
    if(!sc) idx = num_support_r-1; // extrapolate last quadratic function ATTENTION TODO this was -1, not -2
//     printf("Interpolate: r = %lf, support_r[%i] = %lf\n", r, idx, support_r[idx]); fflush(stdout);
    const double dr=r-support_r[idx];

    dS      = poly_b[idx]+dr*2.*poly_c[idx];
    ddS     = 2.*poly_c[idx];
//     printf("Sq: r = %4.2f, idx = %i, poly = %f\n", r, idx, poly_a[idx]+dr*(poly_b[idx]+dr*poly_c[idx]));
    return poly_a[idx]+dr*(poly_b[idx]+dr*poly_c[idx]);
}

void QuadraticInterpolant::GetSupport( double * a_support_r, int & a_num_support_r ) const
{
    for(int iradius=0; iradius<num_support_r; iradius++)
        a_support_r[iradius] = support_r[iradius];
    a_num_support_r = num_support_r;
}
