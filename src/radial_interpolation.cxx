#include <radial_interpolation.h>

using namespace UTILITY;

Interpolant::Interpolant() {
    DefaultInit();
}

Interpolant::~Interpolant() {
    Free();
}

void Interpolant::DefaultInit() {
    X   = 0;
    m_D = 0;
    dim = 0;
    n   = 0;
}

void Interpolant::Free() {
    free_array( &X );
    free_array( &m_D );
    dim = 0;
    n   = 0;    
}

void Interpolant::Allocate( const int a_n, const int a_dim )
{
    assert_msg((a_n>=2), "ERROR in Interpolant::Allocate: number of support points >= 2 expected\n");
    assert_msg((a_dim>=1), "ERROR in Interpolant::Allocate: data-dimension >= 1 expected\n");

    Free();
    dim = a_dim;
    n   = a_n;
    X   = alloc_array(n);
    m_D = alloc_array(2*n*dim);
}

void Interpolant::Train(
    int n_x,        const double * const x,
    int a_dim,      const double * const m_data )
{
    Free();
    Allocate( n_x-1, a_dim );
    for( int i=0; i<n; i++ ) X[i] = x[i];
    // NOTE:    for x > x_max --> extrapolate linearly, i.e. last point is not stored
    //          (assumed to continue until infinity)
    
    size_t ct = 0;
    for( int i=0;i<n; i++ )
    {
        for(int j=0; j<dim; j++)
        {
            m_D[ct++] = m_data[i*dim+j];
            m_D[ct++] = (m_data[(i+1)*dim+j] - m_data[i*dim+j])/(x[i+1]-x[i]);
        }
    }
}

void Interpolant::Dim( const int d ) {
    // setting the dimension implies clearing existing data
    Free();
    dim = d;
}

int Interpolant::Dim() const { return dim; }

void Interpolant::InterpolationWeights( const int idx, const double x, double * o_w ) const
{
    assert_msg( ((idx >=0) || (idx<n)), "ERROR in InterpolationWeights(const int, const double, double *): idx out of bound\n");
    o_w[0] = 1.;
    o_w[1] = x-X[idx];
}
void Interpolant::Interpolate( const double x, double * v_Data )
{
    int idx = Interval( x );
    InterpolationWeights( idx, x, w );
    size_t ct = 2*dim*idx;
    for(int d = 0; d<dim; d++)
    {
        v_Data[d] = m_D[ct] + m_D[ct+1]*w[1];
        ct += 2;
    }
}
void Interpolant::InterpolationWeights( const double x, double * o_w) const
{
    int idx = Interval( x );
    InterpolationWeights( idx, x, o_w );
}

/* ************************************************************************ */
/* ************************************************************************ */
/* ************************************************************************ */
/* ************************************************************************ */

InterpolantQuad::InterpolantQuad()
{
    DefaultInit();
    
}InterpolantQuad::~InterpolantQuad()
{
    Free();
}

void InterpolantQuad::Allocate( const int a_n, const int a_dim )
{
    assert_msg((a_dim>=1), "ERROR in InterpolantQuad::Allocate: data-dimension >= 1 expected\n");

    Free();
    dim = a_dim;
    n   = a_n;
    X   = alloc_array(n);
    m_D = alloc_array(3*n*dim);
}

void InterpolantQuad::Train( int        n_x,
                const double * const    x,
                int                     a_dim,
                const double * const    m_Data        )
{
    Free();
    assert_msg((n_x>=3), "ERROR in InterpolantQuad::Train: number of support points >= 3 expected\n");
    assert_msg( int((n_x-1)/2)*2 == n_x-1, "ERROR in InterpolantQuad::Train(): n_x must be in {3, 5, 7, ...}\n");
    
    Allocate( (n_x-1)/2, a_dim );
    for( int i=0; i<n; i++ ) X[i] = x[2*i];
    // NOTE:    for x > x_max --> extrapolate linearly, i.e. last point is not stored
    //          (assumed to continue until infinity)
    
    size_t ct = 0;
    for( int i=0;i<n; i++ )
    {
        const double dx1 = x[2*i+1]-x[2*i], dx2 = x[2*i+2]-x[2*i];
        for(int j=0; j<dim; j++)
        {
            m_D[ct+0] = m_Data[(2*i)*dim+j];
            m_D[ct+2] = ( m_Data[(2*i+1)*dim+j]*dx2 - m_Data[(2*i+2)*dim+j]*dx1 - (dx2-dx1)*m_D[ct+0] ) / ( dx1*dx2*(dx1-dx2));
            m_D[ct+1] = ( m_Data[(2*i+1)*dim+j] - m_D[ct+0] - dx1*dx1*m_D[ct+2] )/dx1;
            ct += 3;
        }
    }
    print_matrix( X, 1, n);
    print_matrix( m_D, n, 3);
}


void InterpolantQuad::TrainC1( int        n_x,
                const double * const    x,
                int                     a_dim,
                const double * const    m_Data        )
{
    Free();
    assert_msg((n_x>=3), "ERROR in InterpolantQuad::TrainC1: number of support points >= 3 expected\n");
    
    Allocate( n_x-2, a_dim );
    X[0] = x[0];
    for( int i=1; i<n; i++ ) {
        X[i] = x[i+1];
        printf("X[%2i] = %10.5f\n", i, X[i] );
    }
    // NOTE:    for x > x_max --> extrapolate linearly, i.e. last point is not stored
    //          (assumed to continue until infinity)
    
    size_t ct = 0;
    // first interval ranging from x[0] to x[2]:
    const double dx1 = x[1]-x[0], dx2= x[2]-x[0];
    for(int j=0; j<dim; j++)
    {
        m_D[ct+0] = m_Data[j];
        m_D[ct+2] = ( m_Data[dim+j]*dx2 - m_Data[2*dim+j]*dx1 - (dx2-dx1)*m_D[j] ) / ( dx1*dx2*(dx1-dx2));
        m_D[ct+1] = ( m_Data[dim+j] - m_D[ct] - dx1*dx1*m_D[ct+2] )/dx1;
        ct = ct + 3;
    }

    for( int i=1;i<n; i++ )
    {
        const double * y = m_Data + (i+1)*dim;
        const double * y_next = m_Data + (i+2)*dim;
        const double dx_old = X[i]-X[i-1], dx = x[i+2]-x[i+1];
        printf("i: %2i, dx_old %8.4f, dx %8.4f\n", i, dx_old, dx);
        // j --> component of the data
        for(int j=0; j<dim; j++)
        {
            // NOTE:
            // i corresponds to input site (i+1) [first interval used 3 samples]
            // y_i[j] = m_D[ 3*i*dim + 3*j + 0]
            // dy_i[j] = m_D[ 3*i*dim + 3*j + 1]
            // dy_i[j] = m_D[ 3*i*dim + 3*j + 2]
            // offset
            m_D[ct+0] = y[j]; 
            // dy_i+1 = dy_i + 2 * dx_old * ddy_i
            m_D[ct+1] = m_D[3*(j+(i-1)*dim) + 1] + 2. * dx_old * m_D[3*(j+(i-1)*dim) + 2];
            // ddy_i+1 = ( y_i+2 - ( y_i+1 + dy_i+1 * dx ) ) / ( dx^2 )
            m_D[ct+2] = (y_next[j] - (m_D[ct+0] + m_D[ct+1] * dx ) )/( dx * dx );
            ct = ct + 3;
        }
    }
    print_matrix( X, n, 1);
    print_matrix( m_D, n, 3);
}

void InterpolantQuad::InterpolationWeights( const int idx, const double x, double * o_w ) const
{
    assert_msg( ((idx >=0) || (idx<n)), "ERROR in InterpolationWeights(const int, const double, double *): idx out of bound\n");
    o_w[0] = 1.;
    o_w[1] = x-X[idx];
    o_w[2] = o_w[1]*o_w[1];
}
void InterpolantQuad::Interpolate( const double x, double * v_Data )
{
    int idx = Interval( x );
    printf("idx: %i, n: %i\n", idx, n);
    InterpolationWeights( idx, x, w );
    size_t ct = 3*dim*idx;
    for(int d = 0; d<dim; d++)
    {
        v_Data[d] = m_D[ct] + m_D[ct+1]*w[1] + m_D[ct+2]*w[2];
        ct += 3;
    }
}
void InterpolantQuad::InterpolationWeights( const double x, double * o_w) const
{
    int idx = Interval( x );
    InterpolationWeights( idx, x, o_w );
}

/* ************************************************************************ */
/* ************************************************************************ */
/* ************************************************************************ */
/* ************************************************************************ */
