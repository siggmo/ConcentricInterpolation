#include <util.h>

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

using namespace UTILITY;

/* *************************************************************************************** */
void    UTILITY::cleanString( char * in, char * out )
{
    // remove comments (after # symbol), replace double whitespaces by whitespaces, ....
    // [ NOTE: , and TAB (\t) are also considered as a white space character! ]
    int i=0, j=0;
    int l=strlen(in);
    if( l == 0 ){ out[0]='\0'; return; }
    
    // transfer eol character to \0
    if( in[l-1] == '\n' ) { in[l-1]='\0'; l--; }
    
    bool sc=true, ws=false;
    for(i=0;(i<l)&&(sc);i++)
    {
        sc = !(in[i] == '#') && !(in[i] == '\n');
        if( sc ) 
        {
            if( in[i] == ' ' || in[i] == ',' || in[i] == '\t' ) {
                if((j!=0) && (!ws)) out[j++] = ' '; // eliminate also leading whitespace characters
                ws = true;
            }
            else {
                ws = false;
                out[j++] = in[i]; // copy character
            }
        }
    }
    if( ws ) j--; // remove trailing whitespace
    out[j] = '\0';
}
/* *************************************************************************************** */
char * UTILITY::ToNextWhitespace( char * in )
{
    char * out = in;
    while( (out[0] != ' ') && (out[0] != '\t') && (out[0] != '\n') && (out[0] != '\0') ) out++;
    if( (out[0] == '\0') || (out[0]=='\n') ) return 0;

    // skip the whitespace:
    out++;
    return out;
}
/* *************************************************************************************** */
double ** UTILITY::ReadMatrix( int *r, int *c, const char * fn, const int NMAX_LINES, const int NMAX_COL )
{
    FILE * F = fopen(fn, "r");
    assert_msg( F != 0, "ERROR opening file in ReadMatrix\n");
    
    char    line[8192], buffer[8192], *dummy;
    char    * c_str = 0;
    bool    sc      = false;
    // skip all comment lines, i.e. lines commencing with # (or whitespaces, i.e. ' ', '\t' followed by #)
    while( !sc && !feof(F) )
    {
        dummy = fgets(line, 8192, F);
        // remove leading white spaces:
        int ndel=0, l=strlen(line);
        while( (line[ndel]==' '||line[ndel]=='\t') && ndel<l) ndel++;
        c_str = line + ndel; // skip the first ndel characters
        sc = ! ( (c_str[0] == '#') || (l ==  ndel) || (c_str[0] =='\n')  || (c_str[0]=='\0') );
    }
    assert_msg( !feof(F), "ReadMatrix reached end of file before finding useful data.\n");
    
    // first step: find the right number of columns from the first row
    cleanString( line, buffer );
    int nrow=0, ncol=0;
    int l=strlen(buffer);
    for(int i=0;i<l;i++) ncol+=(buffer[i]==' ');
    ncol++;
    double *tmp_matrix=new double[NMAX_LINES*NMAX_COL];
    sc = true;
    size_t ct=0;
    while( sc )
    {
        c_str=buffer;
        for(int i=0;(i<ncol) && sc;i++)
        {
            sc = sscanf(c_str, "%lf", tmp_matrix + ct ); // sc = 0 if reading a double fails
            assert_msg( sc, "Error reading matrix from file: not numerical data/input error\n");
            ct++;
            if( i<ncol-1) {
                c_str   = ToNextWhitespace( c_str );
                sc      = (c_str != 0);
            }
        }
        if( sc ) {
            nrow++;
            dummy = fgets( line, 8192, F );
            if( feof(F) ) sc=false;
            else {
                cleanString( line, buffer );
                while( (strlen(buffer)==0) && (!feof(F)) )
                {
                    dummy = fgets(line, 8192, F);
                    cleanString( line, buffer );
                }
                if( buffer[0] == '\0' ) sc=false;
                else {
                    int tmp_ncol=1;
                    for(int i=0;i<strlen(buffer);i++) tmp_ncol+=(buffer[i]==' ');
                    sc=(tmp_ncol==ncol);
                }
            }
        }
    }
    *r  = nrow;
    *c  = ncol;
    
    double ** matrix = alloc_matrix(nrow, ncol);
    for(size_t irow=0; irow<nrow; irow++)
        for(size_t icol=0; icol<ncol; icol++)
            matrix[irow][icol] = tmp_matrix[irow*ncol+icol];

    delete [] tmp_matrix;
    tmp_matrix = 0;

    fclose(F);
    F=0;
    
    return matrix;
}
/* *************************************************************************************** */
void UTILITY::print_matrix( const double * a, const int m, const int n, FILE * F )
{
    assert_msg( m>0, "error in print_matrix: m>0 required\n");
    assert_msg( n>0, "error in print_matrix: n>0 required\n");
    assert_msg( a!=0, "error in print_matrix: data pointer is NULL\n");
    assert_msg( F!=0, "error in print_matrix: file pointer is NULL\n");
    for(int i=0;i<m;i++)
    {
        fprintf(F,"%24.17le", a[ IDX2( i, 0, n ) ] );
        for(int j=1;j<n;j++)
            fprintf(F,"  %24.17le", a[ IDX2( i, j, n ) ] );
        fprintf(F,"\n");
    }
}
/* *************************************************************************************** */
void UTILITY::assert_msg( const bool condition, const char * str, const bool quit, FILE * outputstream ) 
{
    if(condition)
        return;
    fprintf(outputstream, "%s", str );
    fflush( outputstream );
    
    if(quit)
        exit( DEFAULT_ERROR_CODE );
}
/* *************************************************************************************** */
double * UTILITY::alloc_array(const int N )
{
    assert_msg( N>0, "ERROR: array must have positive dimension\n");
    double * a = 0;
    a = new double [N];
    assert_msg( a != 0, "ERROR in alloc_array (NULL pointer returned)\n");
    return a;
}
/* *************************************************************************************** */
void    UTILITY::free_array( double ** a )
{
    delete [] *a;
    *a = 0;
}
/* *************************************************************************************** */
double      UTILITY::distance( const double * a, const double * b, const int dim ) 
{
    if( dim==0 ) return 0.;
    double res = 0.;

#pragma unroll (4)
    for( int i=0; i<dim; i++ ) res += (a[0]-b[0])*(a[0]-b[0]);

    return sqrt(res);
}
/* *************************************************************************************** */
double  UTILITY::norm( const double * a, const int dim ) 
{
    if( dim==0 ) return 0.;
    double res = 0.;

#pragma unroll (4)
    for( int i=0; i<dim; i++ ) res+=a[i]*a[i];

    return sqrt(res);
}
/* *************************************************************************************** */
void    UTILITY::SetVector( double * a, const int N, const double a0 )
{

#pragma unroll (4)
    for(int i=0;i<N;i++) a[i] = a0;
}
/* *************************************************************************************** */
double      UTILITY::MatVecMul( const double * A, const double * x, double * y, const int M, const int N, const bool T ) 
{
    if( T )
    {
        SetVector(y, N, 0. );
        for(int i=0;i<M; i++ )
        {
            #pragma unroll (4)
            for(int j=0;j<N; j++ )
                y[j]+=A[i*N+j]*x[i];
        }
    }
    else
    {
        for(int i=0;i<M;i++)
            y[i] = VecVecMul( A + i*N, x, N );
    }
}
/* *************************************************************************************** */
void    UTILITY::Factorize( const double * A, double * Af, int * w_i, const int N )
{
#pragma unroll (4)
    for(int n=0;n<N*N;n++) Af[n] = A[n];
    
    int info = LAPACKE_dsytrf( LAPACK_ROW_MAJOR, 'L', N, Af, N, w_i );
}
/* *************************************************************************************** */
void    UTILITY::SolveByFactorization( const double * Af, const double * a, double * Ai_a, int * w_i, const int N, const int Nrhs )
{
#pragma unroll (4)
    for(int n=0;n<N;n++) Ai_a[n] = a[n];
    
    int info = LAPACKE_dsytrs( LAPACK_ROW_MAJOR, 'L', N, Nrhs, Af, N, w_i, Ai_a, Nrhs );
}
/* *************************************************************************************** */
double      FastAcos::data[2*NSTEPS+1];
double      FastAcos::xval[2*NSTEPS];
double      FastAcos::mid_data[2*NSTEPS+1];

void FastAcos::Init()
{
    for(int i=0;i<2*NSTEPS;i++)
    {
        FastAcos::xval[i]=-1.+double(i)*FastAcos::h;
        FastAcos::data[i]=safeAcos(FastAcos::xval[i]);
        FastAcos::mid_data[i]=safeAcos(FastAcos::xval[i]+0.5*h);
    }
    FastAcos::data[2*NSTEPS] = 0.;
    FastAcos::mid_data[2*NSTEPS] = 0.;
    init=true;
}

double FastAcos::Eval(const double x )
{
    if(x<=-1.)      return PI;
    if(x>=1.)       return 0.;
    const int       i       = int((x+1.)/h);
    const double    alpha   = (x-xval[i])/h;
    
    // quadratic interpolation function (via Lagrange polynomials
    return ((data[i] * (alpha-0.5) - 2.*alpha* mid_data[i] ) *(alpha-1.0) + alpha*(alpha-0.5)*data[i+1])*2.;
}
/* *************************************************************************************** */
double ** UTILITY::alloc_matrix( const int m, const int n )
{
    if( (m<=0) || (n<=0) )
    {
        double ** A = 0;
        return A;
    }
    int i;
    double ** A = new double * [m];
    A[0] = new double [ size_t(m)*size_t(n) ];
    assert_msg( A[0] != 0 , "ERROR in double ** alloc_matrix() : out of memory error\n");

    for(i=0;i<m;i++) A[i] = A[0] + i*n;
    
    memset(A[0], 0, sizeof(double)*size_t(m)*size_t(n)); /*set to zero*/
    return A;
}
/* *************************************************************************************** */
void UTILITY::free_matrix( double ** &A, const size_t  m )
{
    assert_msg( ( A == 0 ) || ( m >= 0 ) ,"ERROR: Trying to free an array of negative or zero size in free_matrix(double **&, const size_t)\n");

    if( m == 0 ) return;
    if( A != 0 ) {
        for(size_t i=1;i<m;i++) A[i] = 0;
        delete [] A[0]; A[0] = 0;
        delete [] A;    A = 0;
    }
}

/* *************************************************************************************** */
/* auxiliary namespace that can be used to track the progress of a set operation, i.e. the state of TestOnSet or InterpolateOnSet */
namespace UTILITY{
    namespace progress{
        double progress_new;
        double progress_old;
    }
};
/* *************************************************************************************** */
void UTILITY::progress::init()
{
    progress_old = -1;
    progress_new = -1;
}

void UTILITY::progress::display(const int current, const int total)
{
    progress_new = floor( double(current)/double(total) * 10. );
    if( progress_new > progress_old )
    {
        progress_old = progress_new;
        printf("\r$   finished%3i%%", int(progress_old)*10); // NOTE: the percentage output has been tested under Linux with a standard bash shell
        fflush(stdout);
    }
//  if(current==total)
//      printf("\n");
}
