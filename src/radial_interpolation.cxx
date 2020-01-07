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

#include <radial_interpolation.h>

using namespace UTILITY;

Interpolant::Interpolant(const int a_N_term):
    N_term(a_N_term)
{
    DefaultInit();
}

Interpolant::~Interpolant() {
    Free();
}

void Interpolant::DefaultInit() {
    Radii       = 0;
    coeff       = 0;
    D_val       = 0;
    N_rad_supp  = 0;
}

void Interpolant::Free() {
    free_array( &Radii );
    free_array( &coeff );
    D_val       = 0;
    N_rad_supp  = 0;
}

void Interpolant::Allocate( const int a_n, const int a_D_val )
{
    assert_msg((a_n>=2), "ERROR in Interpolant::Allocate: number of support points >= 2 expected\n");
    assert_msg((a_D_val>=1), "ERROR in Interpolant::Allocate: data-dimension >= 1 expected\n");

    Free();
    D_val       = a_D_val;
    N_rad_supp  = a_n;
    Radii       = alloc_array(N_rad_supp);
    coeff       = alloc_array(N_term*N_rad_supp*D_val);
}

void Interpolant::Train(    const int a_N_rad_supp,
                            const double * const a_Radii,
                            const int a_D_val,
                            const double * const * m_data   )
{
    Free();
    Allocate( a_N_rad_supp-1, a_D_val );  // NOTE: input MUST contain zero, but it is not sto
    for( int i_point=0; i_point<N_rad_supp; i_point++ ) Radii[i_point] = a_Radii[i_point];
    // NOTE:    for a_Radii > x_max --> extrapolate linearly, i_point.e. last point is not stored
    //          (assumed to continue until infinity)

    size_t ct = 0;
    for( int i_point=0;i_point<N_rad_supp; i_point++ )
    {
        for(int i_comp=0; i_comp<D_val; i_comp++)
        {
            coeff[ct++] = m_data[i_point][i_comp];
            coeff[ct++] = (m_data[i_point+1][i_comp] - m_data[i_point][i_comp])/(a_Radii[i_point+1]-a_Radii[i_point]);
//             coeff[ct++] = m_data[i_point*D_val+i_comp];
//             coeff[ct++] = (m_data[(i_point+1)*D_val+i_comp] - m_data[i_point*D_val+i_comp])/(a_Radii[i_point+1]-a_Radii[i_point]);
        }
    }
}

void Interpolant::CompData( const int d ) {
    // setting the dimension implies clearing existing data
    Free();
    D_val = d;
}

void Interpolant::Interpolate( const double a_Radii, double * v_Data )
{
    int idx = Interval( a_Radii );
    InterpolationWeights( idx, a_Radii, weights );
    size_t ct = 2*D_val*idx;
    for(int i_comp = 0; i_comp<D_val; i_comp++)
    {
        v_Data[i_comp] = coeff[ct] + coeff[ct+1]*weights[1]; // weights[0] is always 1, see InterpolationWeights(const int, const double, double*)
        ct += 2;
    }
}
void Interpolant::InterpolationWeights( const double a_Radii, double * o_w) const
{
    int idx = Interval( a_Radii );
    InterpolationWeights( idx, a_Radii, o_w );
}

void Interpolant::InterpolationWeights( const int idx, const double a_Radii, double * o_w ) const
{
    assert_msg( ((idx >=0) || (idx<N_rad_supp)), "ERROR in InterpolationWeights(const int, const double, double *): idx out of bound\n");
    o_w[0] = 1.;
    o_w[1] = a_Radii-Radii[idx];
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

void InterpolantQuad::Train( int        a_N_rad_supp,
                const double * const    a_Radii,
                int                     a_D_val,
                const double * const    m_Data        )
{
    Free();
    assert_msg((a_N_rad_supp>=3), "ERROR in InterpolantQuad::Train: number of support points >= 3 expected\n");
    assert_msg( int((a_N_rad_supp-1)/2)*2 == a_N_rad_supp-1, "ERROR in InterpolantQuad::Train(): a_N_rad_supp must be in {3, 5, 7, ...}\n");

    Allocate( (a_N_rad_supp-1)/2, a_D_val );
    for( int i_point=0; i_point<N_rad_supp; i_point++ ) Radii[i_point] = a_Radii[2*i_point];
    // NOTE:    for a_Radii > x_max --> extrapolate linearly, i.e. last point is not stored
    //          (assumed to continue until infinity)

    size_t ct = 0;
    for( int i_point=0;i_point<N_rad_supp; i_point++ )
    {
        const double dx1 = a_Radii[2*i_point+1]-a_Radii[2*i_point], dx2 = a_Radii[2*i_point+2]-a_Radii[2*i_point];
        for(int i_comp=0; i_comp<D_val; i_comp++)
        {
            coeff[ct+0] = m_Data[(2*i_point)*D_val+i_comp];
            coeff[ct+2] = ( m_Data[(2*i_point+1)*D_val+i_comp]*dx2 - m_Data[(2*i_point+2)*D_val+i_comp]*dx1 - (dx2-dx1)*coeff[ct+0] ) / ( dx1*dx2*(dx1-dx2));
            coeff[ct+1] = ( m_Data[(2*i_point+1)*D_val+i_comp] - coeff[ct+0] - dx1*dx1*coeff[ct+2] )/dx1;
            ct += 3; // 3 = N_term
        }
    }
    print_matrix( Radii, 1, N_rad_supp);
    print_matrix( coeff, N_rad_supp, 3);
}


void InterpolantQuad::TrainC1( int      a_N_rad_supp,
                const double * const    a_Radii,
                int                     a_D_val,
                const double * const    m_Data        )
{
    Free();
    assert_msg((a_N_rad_supp>=3), "ERROR in InterpolantQuad::TrainC1: number of support points >= 3 expected\n");

    Allocate( a_N_rad_supp-2, a_D_val );
    Radii[0] = a_Radii[0];
    for( int i=1; i<N_rad_supp; i++ )
        Radii[i] = a_Radii[i+1];
    // NOTE:    for a_Radii > x_max --> extrapolate linearly, i.e. last point is not stored
    //          (assumed to continue until infinity)

    size_t ct = 0;
    // first interval ranging from a_Radii[0] to a_Radii[2]:
    const double dx1 = a_Radii[1]-a_Radii[0], dx2= a_Radii[2]-a_Radii[0];
    for(int j=0; j<D_val; j++)
    {
        coeff[ct+0] = m_Data[j];
        coeff[ct+2] = ( m_Data[D_val+j]*dx2 - m_Data[2*D_val+j]*dx1 - (dx2-dx1)*coeff[j] ) / ( dx1*dx2*(dx1-dx2));
        coeff[ct+1] = ( m_Data[D_val+j] - coeff[ct] - dx1*dx1*coeff[ct+2] )/dx1;
        ct = ct + N_term;
    }

    for( int i_point=1;i_point<N_rad_supp; i_point++ )
    {
        const double * y = m_Data + (i_point+1)*D_val;
        const double * y_next = m_Data + (i_point+2)*D_val;
        const double dx_old = Radii[i_point]-Radii[i_point-1], dx = a_Radii[i_point+2]-a_Radii[i_point+1];
        // i_comp --> component of the data
        for(int i_comp=0; i_comp<D_val; i_comp++)
        {
            // NOTE:
            // i_point corresponds to input site (i_point+1) [first interval used 3 samples]
            // y_i[i_comp] = coeff[ 3*i_point*D_val + 3*i_comp + 0]
            // dy_i[i_comp] = coeff[ 3*i_point*D_val + 3*i_comp + 1]
            // dy_i[i_comp] = coeff[ 3*i_point*D_val + 3*i_comp + 2]
            // offset
            coeff[ct+0] = y[i_comp];
            // dy_i+1 = dy_i + 2 * dx_old * ddy_i
            coeff[ct+1] = coeff[D_val*(i_comp+(i_point-1)*D_val) + 1] + 2. * dx_old * coeff[D_val*(i_comp+(i_point-1)*D_val) + 2];
            // ddy_i+1 = ( y_i+2 - ( y_i+1 + dy_i+1 * dx ) ) / ( dx^2 )
            coeff[ct+2] = (y_next[i_comp] - (coeff[ct+0] + coeff[ct+1] * dx ) )/( dx * dx );
            ct = ct + N_term;
        }
    }
}

void InterpolantQuad::InterpolationWeights( const int idx, const double a_Radii, double * o_w ) const
{
    assert_msg( ((idx >=0) || (idx<N_rad_supp)), "ERROR in InterpolationWeights(const int, const double, double *): idx out of bound\n");
    o_w[0] = 1.;
    o_w[1] = a_Radii-Radii[idx];
    o_w[2] = o_w[1]*o_w[1];
}
void InterpolantQuad::Interpolate( const double a_Radii, double * v_Data )
{
    int idx = Interval( a_Radii );
    InterpolationWeights( idx, a_Radii, weights );
    size_t ct = N_term*D_val*idx;
    for(int d = 0; d<D_val; d++)
    {
        v_Data[d] = coeff[ct] + coeff[ct+1]*weights[1] + coeff[ct+2]*weights[2];
        ct += N_term;
    }
}
void InterpolantQuad::InterpolationWeights( const double a_Radii, double * o_w) const
{
    int idx = Interval( a_Radii );
    InterpolationWeights( idx, a_Radii, o_w );
}
























/* ************************************************************************ */
/* ************************************************************************ */
/* ************************************************************************ */
/* ************************************************************************ */

RadialInterpolation::RadialInterpolation()
{
    initialized = false;
    DefaultInit();
}
RadialInterpolation::~RadialInterpolation()
{
    Free();
}

void RadialInterpolation::Free()
{
    m_I_length = 0;
    free_array( &m_I );
    free_array( &m_I_old );
    In = nullptr; // TODO destroy object pointed to by In
    free_array( &weights );
    initialized = false;
    DefaultInit();
}

void RadialInterpolation::DefaultInit()
{
    assert_msg( !initialized, "ERROR: RadialInterpolation must not be initialized when calling RadialInterpolation::DefaultInit()\n");
    N_dir   = 0;
    D_val   = 0;
    N_term  = 0;
    N_dir_counter   = 0;
    m_I_length      = 0;
    m_I     = nullptr;
    m_I_old = nullptr;
    In      = nullptr;
    weights = nullptr;
    ordered = false;
    initialized     = true;
}

void RadialInterpolation::Setup(    const double * const * const * a_data, //!< \todo rename data to values
                                    const int a_N_dir,      //!< \todo rename N_dir to N
                                    const int a_n_x,        //!< \todo rename x to R
                                    const int a_dim_data,   //!< \todo rename dim_data to DimValues
                                    const double * a_Radii,
                                    Interpolant & a_In,
                                    TangentialInterpolation & a_TI
                               )
{
    // checks
    assert_msg( a_TI.GetInit(), "ERROR in RadialInterpolation::Setup: TangentialInterpolation must be set up first!\n");
    assert_msg( a_TI.GetN_dir_supp()==a_N_dir, "ERROR in RadialInterpolation::Setup: N of TangentialInterpolation must equal N_dir of RadialInterpolation!\n");

    // allocate memory for the one-dimensional interpolants along each direction
    // train these interpolants along each direction
    // and include these interpolants in the scheme
    for( int i_dir=0; i_dir<a_N_dir; i_dir++)
    {
        a_In.Train( a_n_x, a_Radii, a_dim_data, a_data[i_dir] );
        if( i_dir==0 )
        {
            SetInterpolant( a_In );
            Allocate( a_N_dir );
        }
        AddData( a_In );
    }

    // backup the data before multiplying by inverse Kernel matrix and before re-ordering
    free_array( &m_I_old );
    m_I_old = alloc_array( m_I_length );
    memcpy( m_I_old, m_I, sizeof(double)*m_I_length );

    // multiply by inverse kernel matrix from left and re-order
    MultiplyKernelMatrixAndReorder( a_TI );

    printf("# set up Radial Interpolation for N = %i, R = %i, DimValues = %i\n", a_N_dir, a_n_x, a_dim_data);
}

void RadialInterpolation::SetInterpolant( Interpolant & a_In )
{
    In = &a_In;
}

void RadialInterpolation::Allocate( const int a_N_dir )
{
    assert_msg( In != 0, "ERROR in RedialInterpolation::Allocate: Interpolate must be set, first!\n");

    Interpolant * In_tmp = In;
    Free();
    In = In_tmp; In_tmp = 0; // recover existing interpolant

    N_dir   = a_N_dir;
    D_val   = In->CompData();
    N_term  = In->NTerms();
    N_int   = In->NIntervals();
    // length of m_I, for later use
    m_I_length = size_t(N_dir * N_term * D_val * N_int);
    // allocate memory for interpolation data and for the weights
    m_I     = alloc_array( m_I_length );
    weights       = alloc_array( N_term );
    ordered = false;
}

void RadialInterpolation::AddData( Interpolant & In )
{
    assert_msg( (In.NIntervals() == N_int), "ERROR in RadialInterpolation::AddData: number of intervals is not matching the pre-defined number\n");

    AddData( In.InterpolationData(0) );
}

void RadialInterpolation::AddData( const double * const coeff )
{
    assert_msg( N_dir_counter < N_dir, "ERROR in RadialInterpolation::AddData: number of pre-allocated directions exceeded\n");
    memcpy( m_I + N_dir_counter*N_term*D_val*N_int, coeff, sizeof(double)*size_t(N_term*D_val*N_int) );
    N_dir_counter++;
    ordered = false;
}

void RadialInterpolation::MultiplyKernelMatrixAndReorder( TangentialInterpolation & a_TI )
{
    assert_msg( initialized, "ERROR in RadialInterpolation::MultiplyKernelMatrixAndReorder: must be initialized first\n");

    // get the factorization permutation
    const size_t w_size = 32*N_dir*N_dir; // ATTENTION this must equal the length of a_TI::w_i
    int * temp_w_I = new int[w_size];
    memcpy(temp_w_I, a_TI.GetPermutation(), w_size*sizeof(int));

    // restore backed up m_I from m_I_old. this handles the case of a previous call to Setup(),
    // i.e. when the kernel matrix inverse has been contracted with m_I before.
    memcpy(m_I, m_I_old, sizeof(double)*m_I_length);

    // now contract first index of m_I with second index of inverse kernel matrix
    SolveByFactorization( a_TI.GetKernelMatrixFactorization(), m_I, m_I, temp_w_I, N_dir, N_int*D_val*N_term );
    delete [] temp_w_I, temp_w_I = 0;

    // re-order, swapping first two indices of m_I
    ReorderData( true );
}

void RadialInterpolation::ReorderData(const bool enforce_order)
{
    if( ordered && !enforce_order )
    {
        printf("returning from RadialInterpolation::ReorderData() without reordering\n");
        return;
    }

    double * tmp = alloc_array( N_int * D_val * N_term * N_dir );

    const int old0 = N_term*D_val*N_int, old1 = N_term*D_val, old2 = N_term, old3 = 1;
    const int new0 = N_term*D_val*N_dir, new1 = N_term*D_val, new2 = N_term, new3 = 1;

    for( int i_input=0; i_input < N_dir; i_input++)
    {
        for( int i_x =0; i_x<N_int; i_x++ )
        {
            for(int d=0; d<D_val; d++ )
            {
                for(int i_w=0;i_w<N_term;i_w++)
                {
                    tmp[ i_x*new0 + i_input *new1 + d*new2 + i_w * new3 ] =
                    m_I[ i_input *old0 + i_x*old1 + d*old2 + i_w * old3 ];
                }
            }
        }
    }


    free_array( &m_I );

    m_I = tmp;
    tmp = nullptr;

    ordered = true;
}

void RadialInterpolation::Interpolate( const double a_x, double * m_out )
{
    assert_msg( initialized, "ERROR in in RadialInterpolation::Interpolate: RadialInterpolation is not initialized\n");
    assert_msg( ordered, "ERROR in RadialInterpolation::Interpolate: data is not reorderd (call ReorderData() first)\n");
    assert_msg( N_dir_counter==N_dir, "ERROR in RadialInterpolation::Interpolate: N_dir_counter must equal N_dir\n");

    int idx=In->Interval(a_x);
    In->InterpolationWeights( idx, a_x, weights );
    MatVecMul( m_I + idx * D_val * N_term * N_dir, weights, m_out, D_val*N_dir, N_term );
}
