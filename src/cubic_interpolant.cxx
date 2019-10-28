#include <cubic_interpolant.h>
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

CubicInterpolant::CubicInterpolant() {
    DefaultInit();
}

void CubicInterpolant::DefaultInit()
{
    init    = false;
    num_support_r       = 0;
    support_r       = 0;
    poly_a  = 0;
    poly_b  = 0;
    poly_c  = 0;
    poly_d  = 0;
}
CubicInterpolant::~CubicInterpolant() {
    Free();
}
CubicInterpolant::CubicInterpolant( const int a_num_support_r, const double * a_support_r, const double * S, const double * dS )
{
    DefaultInit();
    assert_msg( a_num_support_r >= 2, "Error in CubicInterpolant: there must be more than 1 data samples for the construction of the interpolant\num_support_r");

    SetData( a_num_support_r, a_support_r, S, dS );
}
void CubicInterpolant::Free()
{
    if(!init) return;
    free_array( &support_r );
    free_array( &poly_a );
    free_array( &poly_b );
    free_array( &poly_c );
    free_array( &poly_d );
    num_support_r = 0;
    init = false;
}
void CubicInterpolant::Allocate()
{
    if(init) Free();
    assert_msg( num_support_r >= 2, "Error in CubicInterpolant::Allocate: num_support_r must be >= 2\n");
    support_r = alloc_array( num_support_r );
    poly_a = alloc_array( num_support_r );
    poly_b = alloc_array( num_support_r );
    poly_c = alloc_array( num_support_r );
    poly_d = alloc_array( num_support_r );
}
void CubicInterpolant::SetData( const int a_num_support_r, const double * a_support_r, const double * S, const double * dS )
{
    if( init ) Free();
    num_support_r = a_num_support_r;
    bool has_zero = ( fabs( a_support_r[0] ) < 1e-16);
    if(!has_zero) num_support_r++;
//  assert_msg( fabs( a_support_r[0] ) > 1e-16, "Error in CubicInterpolant: the zero value is not contained");

    Allocate();
    if(!has_zero) {
        // add data for zero: zero function value and zero slope at support_r=0:
        poly_a[0] = 0.;
        poly_b[0] = 0.;
        // define coefficients c and d by conditions at the end of the intervall [0, a_support_r[0]]
        poly_d[0] = (dS[0]*a_support_r[0] - 2.*S[0]) / (a_support_r[0]*a_support_r[0]*a_support_r[0]);
        poly_c[0] = S[0]/(a_support_r[0]*a_support_r[0]) - poly_d[0]*a_support_r[0];
    }
    double S1,S2,dS1,dS2,dr;
    const int zero_offset=( has_zero ? 0 : 1 );

    // set support_r to a_support_r with initial 0
    support_r[0] = 0.;
    for( int iradius=0; iradius<a_num_support_r; iradius++) support_r[iradius+zero_offset] = a_support_r[iradius];

    for( int iradius=0; iradius<a_num_support_r-1; iradius++)
    {
        if( iradius > 0 ) assert_msg( a_support_r[iradius]<a_support_r[iradius+1], "Error in CubicInterpolant: a_support_r[iradius]<a_support_r[iradius+1] is not guaranteed\n");
        const int i = iradius + zero_offset;
        S1  = S[iradius];       dS1 = dS[iradius];
        S2  = S[iradius+1]; dS2 = dS[iradius+1];
        dr  = a_support_r[iradius+1]-a_support_r[iradius];

        poly_a[i] = S1;
        poly_b[i] = dS1;
        poly_d[i] = ( (dS2+dS1) * dr + 2.*(S1-S2) )/ (dr*dr*dr);
        poly_c[i] = (S2-S1-dS1*dr-poly_d[i]*dr*dr*dr)/(dr*dr);
    }

    // set initialization flag
    init = true;
}

void CubicInterpolant::FirstV( double & v0, double & r ) const
{
    assert_msg( num_support_r > 1, "Error in CubicInterpolant::FirstV: Polynomial not initialized (or too few data)");
    v0  = poly_a[1];
    r   = support_r[1];
}
void CubicInterpolant::FirstDV( double & dv0, double & r ) const
{
//  // slope from the left and from the right --> mean value
//  dv0 = 0.5 * ( poly_b[1] + poly_b[0] + 2.*support_r[1]*poly_c[0] + 3.*support_r[1]*support_r[1]*poly_d[0] );
    // slope from the left
    dv0 = poly_b[0] + 2.*support_r[1]*poly_c[0] + 3.*support_r[1]*support_r[1]*poly_d[0];
    r   = support_r[1];
}
double CubicInterpolant::Value( const double r ) const {
    int idx = 0;
    bool sc = false;
    while( !sc && idx < num_support_r-1 )
    {
        sc = (r <= support_r[idx+1] ) && (r >= support_r[idx]);
        if(!sc) idx++;
    }
    if(!sc) idx = num_support_r-2; // extrapolate second to last cubic function
    const double dx=r-support_r[idx];
    return poly_a[idx]+dx*(poly_b[idx]+dx*(poly_c[idx]+dx*poly_d[idx]));
}

double CubicInterpolant::Interpolate( const double r, double & dS, double & ddS ) const
{
    int idx = 0;
    bool sc = false;
    while( !sc && idx < num_support_r-1 )
    {
        sc = (r <= support_r[idx+1] ) && (r >= support_r[idx]);
        if(!sc) idx++;
    }
    if(!sc) idx = num_support_r-2; // extrapolate second to last cubic function
    const double dr=r-support_r[idx];

    dS      = poly_b[idx]+dr*(2.*poly_c[idx]+3.0*dr*poly_d[idx]);
    ddS     = 2.*poly_c[idx]+6.*dr*poly_d[idx];
    return poly_a[idx]+dr*(poly_b[idx]+dr*(poly_c[idx]+dr*poly_d[idx]));
}

void CubicInterpolant::GetSupport( double * a_support_r, int & a_num_support_r ) const
{
    for(int iradius=0; iradius<num_support_r; iradius++)
        a_support_r[iradius] = support_r[iradius];
    a_num_support_r = num_support_r;
}
