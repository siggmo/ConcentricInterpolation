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

#include <radial_interpolation.h>

using namespace UTILITY;

Interpolant::Interpolant(const int a_n_term):
    n_term(a_n_term)
{
    DefaultInit();
}

Interpolant::~Interpolant() {
    Free();
}

void Interpolant::DefaultInit() {
    X   = 0;
    coeff = 0;
    n_comp = 0;
    n_X   = 0;
}

void Interpolant::Free() {
    free_array( &X );
    free_array( &coeff );
    n_comp = 0;
    n_X   = 0;
}

void Interpolant::Allocate( const int a_n, const int a_n_comp )
{
    assert_msg((a_n>=2), "ERROR in Interpolant::Allocate: number of support points >= 2 expected\n");
    assert_msg((a_n_comp>=1), "ERROR in Interpolant::Allocate: data-dimension >= 1 expected\n");

    Free();
    n_comp = a_n_comp;
    n_X   = a_n;
    X   = alloc_array(n_X);
    coeff = alloc_array(n_term*n_X*n_comp);
}

void Interpolant::Train(    const int a_n_X,
                            const double * const a_X,
                            const int a_n_comp,
                            const double * const m_data   )
{
    Free();
    Allocate( a_n_X-1, a_n_comp );
    for( int i_point=0; i_point<n_X; i_point++ ) X[i_point] = a_X[i_point];
    // NOTE:    for a_X > x_max --> extrapolate linearly, i_point.e. last point is not stored
    //          (assumed to continue until infinity)

    size_t ct = 0;
    for( int i_point=0;i_point<n_X; i_point++ )
    {
        for(int i_comp=0; i_comp<n_comp; i_comp++)
        {
            coeff[ct++] = m_data[i_point*n_comp+i_comp];
            coeff[ct++] = (m_data[(i_point+1)*n_comp+i_comp] - m_data[i_point*n_comp+i_comp])/(a_X[i_point+1]-a_X[i_point]);
        }
    }
}

void Interpolant::CompData( const int d ) {
    // setting the dimension implies clearing existing data
    Free();
    n_comp = d;
}

void Interpolant::Interpolate( const double a_X, double * v_Data )
{
    int idx = Interval( a_X );
    InterpolationWeights( idx, a_X, weights );
    size_t ct = 2*n_comp*idx;
    for(int i_comp = 0; i_comp<n_comp; i_comp++)
    {
        v_Data[i_comp] = coeff[ct] + coeff[ct+1]*weights[1]; // weights[0] is always 1, see InterpolationWeights(const int, const double, double*)
        ct += 2;
    }
}
void Interpolant::InterpolationWeights( const double a_X, double * o_w) const
{
    int idx = Interval( a_X );
    InterpolationWeights( idx, a_X, o_w );
}

void Interpolant::InterpolationWeights( const int idx, const double a_X, double * o_w ) const
{
    assert_msg( ((idx >=0) || (idx<n_X)), "ERROR in InterpolationWeights(const int, const double, double *): idx out of bound\n");
    o_w[0] = 1.;
    o_w[1] = a_X-X[idx];
}
















/* ************************************************************************ */
/* ************************************************************************ */
/* ************************************************************************ */
/* ************************************************************************ */

InterpolantQuad::InterpolantQuad():
    Interpolant(3)
{
    assert_msg( false, "InterpolantQuad not tested yet. -> Quit\n");
    DefaultInit();
}

InterpolantQuad::~InterpolantQuad()
{
    Free();
}

void InterpolantQuad::Train( int        a_n_X,
                const double * const    a_X,
                int                     a_n_comp,
                const double * const    m_Data        )
{
    Free();
    assert_msg((a_n_X>=3), "ERROR in InterpolantQuad::Train: number of support points >= 3 expected\n");
    assert_msg( int((a_n_X-1)/2)*2 == a_n_X-1, "ERROR in InterpolantQuad::Train(): a_n_X must be in {3, 5, 7, ...}\n");

    Allocate( (a_n_X-1)/2, a_n_comp );
    for( int i_point=0; i_point<n_X; i_point++ ) X[i_point] = a_X[2*i_point];
    // NOTE:    for a_X > x_max --> extrapolate linearly, i.e. last point is not stored
    //          (assumed to continue until infinity)

    size_t ct = 0;
    for( int i_point=0;i_point<n_X; i_point++ )
    {
        const double dx1 = a_X[2*i_point+1]-a_X[2*i_point], dx2 = a_X[2*i_point+2]-a_X[2*i_point];
        for(int i_comp=0; i_comp<n_comp; i_comp++)
        {
            coeff[ct+0] = m_Data[(2*i_point)*n_comp+i_comp];
            coeff[ct+2] = ( m_Data[(2*i_point+1)*n_comp+i_comp]*dx2 - m_Data[(2*i_point+2)*n_comp+i_comp]*dx1 - (dx2-dx1)*coeff[ct+0] ) / ( dx1*dx2*(dx1-dx2));
            coeff[ct+1] = ( m_Data[(2*i_point+1)*n_comp+i_comp] - coeff[ct+0] - dx1*dx1*coeff[ct+2] )/dx1;
            ct += 3; // 3 = n_term
        }
    }
    print_matrix( X, 1, n_X);
    print_matrix( coeff, n_X, 3);
}


void InterpolantQuad::TrainC1( int        a_n_X,
                const double * const    a_X,
                int                     a_n_comp,
                const double * const    m_Data        )
{
    Free();
    assert_msg((a_n_X>=3), "ERROR in InterpolantQuad::TrainC1: number of support points >= 3 expected\n");

    Allocate( a_n_X-2, a_n_comp );
    X[0] = a_X[0];
    for( int i=1; i<n_X; i++ ) {
        X[i] = a_X[i+1];
        printf("X[%2i] = %10.5f\n", i, X[i] );
    }
    // NOTE:    for a_X > x_max --> extrapolate linearly, i.e. last point is not stored
    //          (assumed to continue until infinity)

    size_t ct = 0;
    // first interval ranging from a_X[0] to a_X[2]:
    const double dx1 = a_X[1]-a_X[0], dx2= a_X[2]-a_X[0];
    for(int j=0; j<n_comp; j++)
    {
        coeff[ct+0] = m_Data[j];
        coeff[ct+2] = ( m_Data[n_comp+j]*dx2 - m_Data[2*n_comp+j]*dx1 - (dx2-dx1)*coeff[j] ) / ( dx1*dx2*(dx1-dx2));
        coeff[ct+1] = ( m_Data[n_comp+j] - coeff[ct] - dx1*dx1*coeff[ct+2] )/dx1;
        ct = ct + n_term;
    }

    for( int i_point=1;i_point<n_X; i_point++ )
    {
        const double * y = m_Data + (i_point+1)*n_comp;
        const double * y_next = m_Data + (i_point+2)*n_comp;
        const double dx_old = X[i_point]-X[i_point-1], dx = a_X[i_point+2]-a_X[i_point+1];
        printf("i_point: %2i, dx_old %8.4f, dx %8.4f\n", i_point, dx_old, dx);
        // i_comp --> component of the data
        for(int i_comp=0; i_comp<n_comp; i_comp++)
        {
            // NOTE:
            // i_point corresponds to input site (i_point+1) [first interval used 3 samples]
            // y_i[i_comp] = coeff[ 3*i_point*n_comp + 3*i_comp + 0]
            // dy_i[i_comp] = coeff[ 3*i_point*n_comp + 3*i_comp + 1]
            // dy_i[i_comp] = coeff[ 3*i_point*n_comp + 3*i_comp + 2]
            // offset
            coeff[ct+0] = y[i_comp];
            // dy_i+1 = dy_i + 2 * dx_old * ddy_i
            coeff[ct+1] = coeff[n_comp*(i_comp+(i_point-1)*n_comp) + 1] + 2. * dx_old * coeff[n_comp*(i_comp+(i_point-1)*n_comp) + 2];
            // ddy_i+1 = ( y_i+2 - ( y_i+1 + dy_i+1 * dx ) ) / ( dx^2 )
            coeff[ct+2] = (y_next[i_comp] - (coeff[ct+0] + coeff[ct+1] * dx ) )/( dx * dx );
            ct = ct + n_term;
        }
    }
    print_matrix( X, n_X, 1);
    print_matrix( coeff, n_X, n_term);
}

void InterpolantQuad::InterpolationWeights( const int idx, const double a_X, double * o_w ) const
{
    assert_msg( ((idx >=0) || (idx<n_X)), "ERROR in InterpolationWeights(const int, const double, double *): idx out of bound\n");
    o_w[0] = 1.;
    o_w[1] = a_X-X[idx];
    o_w[2] = o_w[1]*o_w[1];
}
void InterpolantQuad::Interpolate( const double a_X, double * v_Data )
{
    int idx = Interval( a_X );
    InterpolationWeights( idx, a_X, weights );
    size_t ct = n_term*n_comp*idx;
    for(int d = 0; d<n_comp; d++)
    {
        v_Data[d] = coeff[ct] + coeff[ct+1]*weights[1] + coeff[ct+2]*weights[2];
        ct += n_term;
    }
}
void InterpolantQuad::InterpolationWeights( const double a_X, double * o_w) const
{
    int idx = Interval( a_X );
    InterpolationWeights( idx, a_X, o_w );
}
























/* ************************************************************************ */
/* ************************************************************************ */
/* ************************************************************************ */
/* ************************************************************************ */

RadialInterpolation::RadialInterpolation()
{
    DefaultInit();
}
RadialInterpolation::~RadialInterpolation()
{
    Free();
}

void RadialInterpolation::Free()
{
    In = 0;
    free_array( &m_I );
    free_array( &weights );
    DefaultInit();
}

void RadialInterpolation::DefaultInit()
{
    weights = 0;
    n_dir   = 0;
    n_comp  = 0;
    n_term  = 0;
    n_input = 0;
    m_I     = 0;
    In      = 0;
    ordered = false;
}

void RadialInterpolation::Setup(    const double * const *  a_data,
                                    const int a_n_dir,
                                    const int a_n_x,
                                    const int a_dim_data,
                                    const double * X,
                                    Interpolant & a_In,
                                    TangentialInterpolation & a_TI
                               )
{
    assert_msg( a_TI.GetInit(), "ERROR in RadialInterpolation::Setup: TangentialInterpolation must be set up first!\n");
    assert_msg( a_TI.GetN()==a_n_dir, "ERROR in RadialInterpolation::Setup: N of TangentialInterpolation must equal n_dir of RadialInterpolation!\n");
    for( int i_dir=0; i_dir<a_n_dir; i_dir++)
    {
        a_In.Train( a_n_x, X, a_dim_data, a_data[i_dir] );
        if( i_dir==0 )
        {
            SetInterpolant( a_In );
            Allocate( a_n_dir );
        }
        AddData( a_In );
    }
    // multiply inverse kernel matrix from left, i.e. contract with first index
    const size_t w_size = 32*n_dir*n_dir;
    int * temp_w_I = new int[w_size];
    memcpy(temp_w_I, a_TI.GetPermutation(), w_size);
    SolveByFactorization( a_TI.GetKernelMatrixFactorization(), m_I, m_I, temp_w_I, n_dir, n_int*n_comp*n_term );
    printf("WARNING: multiplied with Kinv, but not considered in TI --> DO NOT USE TangentialInterpolation::Weights()\n");
    delete [] temp_w_I, temp_w_I = 0;
    // re-order, swapping first two indices
    ReorderData();
}

void RadialInterpolation::SetInterpolant( Interpolant & a_In )
{
    In = &a_In;
}

void RadialInterpolation::Allocate( const int a_n_dir )
{
    assert_msg( In != 0, "ERROR in RedialInterpolation::Allocate: Interpolate must be set, first!\n");

    Interpolant * In_tmp = In;
    Free();
    In = In_tmp; In_tmp = 0; // recover existing interpolant

    n_dir   = a_n_dir;
    n_comp  = In->CompData();
    n_term  = In->NTerms();
    n_int   = In->NIntervals();
    // allocate memory for interpolation data and for the weights
    m_I     = alloc_array( n_dir * n_term * n_comp * n_int );
    weights       = alloc_array( n_term );
    ordered = false;
}

void RadialInterpolation::AddData( Interpolant & In )
{
    assert_msg( (In.NIntervals() == n_int), "ERROR in RadialInterpolation::AddData: number of intervals is not matching the pre-defined number\n");

    AddData( In.InterpolationData(0) );
}

void RadialInterpolation::AddData( const double * const coeff )
{
    assert_msg( n_input < n_dir, "ERROR in RadialInterpolation::AddData: number of pre-allocated directions exceeded\n");
    memcpy( m_I + n_input*n_term*n_comp*n_int, coeff, sizeof(double)*n_term*n_comp*n_int );
    n_input++;
    ordered = false;
}

void RadialInterpolation::ReorderData()
{
    if( ordered )
    {
        printf("returning from RadialInterpolation::ReorderData() without reordering\n");
        return;
    }

    double * tmp = alloc_array( n_int * n_comp * n_term * n_dir );

    const int old0 = n_term*n_comp*n_int, old1 = n_term*n_comp, old2 = n_term, old3 = 1;
    const int new0 = n_term*n_comp*n_dir, new1 = n_term*n_comp, new2 = n_term, new3 = 1;

    for( int i_input=0; i_input < n_dir; i_input++)
    {
        for( int i_x =0; i_x<n_int; i_x++ )
        {
            for(int d=0; d<n_comp; d++ )
            {
                for(int i_w=0;i_w<n_term;i_w++)
                {
                    tmp[ i_x*new0 + i_input *new1 + d*new2 + i_w * new3 ] =
                    m_I[ i_input *old0 + i_x*old1 + d*old2 + i_w * old3 ];
                }
            }
        }
    }


    free_array( &m_I );

    m_I = tmp;
    tmp = 0;

    ordered = true;
    printf("reordered RadialInterpolation\n");
}

void RadialInterpolation::Interpolate( const double a_x, double * m_out )
{
    assert_msg( ordered, "ERROR in RadialInterpolation::Interpolate: data is not reorderd (call ReorderData() first)\n");
    assert_msg( n_input==n_dir, "ERROR in RadialInterpolation::Interpolate: n_input must equal n_dir\n");
    int idx=In->Interval(a_x);
    In->InterpolationWeights( idx, a_x, weights );
    MatVecMul( m_I + idx * n_comp * n_term * n_dir, weights, m_out, n_comp*n_dir, n_term );
}
