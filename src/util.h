#ifndef _UTIL_H_
#define _UTIL_H_

/*
 *  COPYRIGHT NOTES
 * 
 *  ConcentricInterpolation
 *  Copyright (C) 2019  Felix Fritzen    ( fritzen@mechbau.uni-stuttgart.de )
 *                      and Oliver Kunc  ( kunc@mechbau.uni-stuttgart.de )
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
 *  (the full license is distributed together with the software in a file named
 *  LICENSE)
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

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <math.h>
#include <string.h>
#include <unistd.h>
#include <sys/stat.h> 
#include <lapacke.h>
#include <random>
#include <chrono>

#ifndef PI
#define PI 3.14159265358979
#endif

#define DEFAULT_ERROR_CODE 10

namespace UTILITY {

inline int IDX2( int i, int j, int n ) { return( i*n+j ); }

//! jump to the position after the next whitespace character in in[]
char *      ToNextWhitespace( char * in );
//! remove comments, leading/trailing whitespace characters, replace multiple whitespaces by single whitespaces and replace separators (,) by spaces
void        cleanString( char * in, char * out );
//! read a matrix from a text file (matrix dimensions are unknown)
double **   ReadMatrix(  int *r             /*!<[out] number of read rows*/,
                         int *c             /*!<[out] number of read columns*/,
                         const char * fn,   /*!<[in] name of the file containing the directions*/
                         const int NMAX_LINES=32768, const int NMAX_COL=32 );
//! dump a pseudo-2d-matrix or an array to a stream. if array, set m = 1.
void        print_matrix(   const double * const a/*!<[in] 1-d or pseudo 2-d pointer*/,
                            const int m/*!<[in] number of rows*/,
                            const int n/*!<[in] number of columns*/,
                            FILE * F = stdout/*!<[in] stream pointer*/);
//! simple assert command: prints str to outputstream if condition is false. if quit is true, it then also exits with DEFAULT_ERROR_CODE
void        assert_msg(const bool condition, const char *str, const bool quit=true, FILE * outputstream=stderr );


/* *************************************************************************************** */
/* memory management */
double *    alloc_array( const int N )  ; //!< allocate array of length \c N (N is the number of doubles to allocate memory for)
void        free_array( double ** )  ; //!< savely free previously allocated array

//! \brief function for allocating pseudo 2-dimensional arrays. memory is filled with zeros.
//! \details \see free_matrix
double ** alloc_matrix( const int m, const int n );

//! \brief function for freeing memory allocated by alloc_matrix
//! \details \see alloc_matrix
void free_matrix( double ** &A, const size_t  m /**< [in] first dimension of A */ );
/* *************************************************************************************** */
double      norm( const double * a, const int dim )  ; //!< compute the \f$l^2\f$ vector norm of \c a in R^\c dim
double      distance( const double * a, const double * b, const int dim )  ; //!< compute the \f$l^2\f$ vector distance of vector \c a and \c b

// interface to LAPACKE and BLAS
void        SetVector( double * a, const int N, const double a0 )  ; //!< initialize all components of a vector, i.e. \f$a_i=a_0/f$
inline double VecVecMul( const double * a, const double * b, const int N ) //!< return inner product \f$a\cdot{}b\f$
{
    if( N==0 ) return 0;
    double res = 0.;

#pragma unroll (4)
    for( int i=0; i<N; i++ ) res += a[i]*b[i];

    return res;
}
/* *************************************************************************************** */
double      MatVecMul( const double * A, const double * x, double * y, const int M, const int N, const bool T=false ); //!< on exit, \f$y_i=A_{ij}x_j\f$ or (if \c T==true) \f$y_j=A_{ij}x_i\f$

/* *************************************************************************************** */
/* use LDL solver to compute the product of the inverse kernel matrix and a vector */
//! stores LDL factorization of \c A to \c Af, and corresponding interchange indices to the integer array \c w_i. \c A and \c Af need to be of size \c N*\c N and the size of \c w_i should be at least the same size or bigger
void        Factorize(  const double * A /*!< [in] N*N matrix*/,
                        double * Af /*!< [out] LDL factorization of A */,
                        int * w_i /*!< [out] integer working array of size >= N*N*/,
                        const int N /*!< [in] dimension of A */);
//! solves the system A * Ai_a = a for \c Ai_a by the LDL factorization \c Af ( computed by \see Factorize )
void        SolveByFactorization( const double * Af /*!< [in] LDL factorization of A; \see Factorize */,
                                  const double * a /*!< [in] right hand side vector (or matrix) */,
                                  double * Ai_a /*!< [out] solution of the linear system A*Ai_a=a */,
                                  int * w_i /*!< [in] permutation vector; \see Factorize */,
                                  const int N /*!< [in] dimension of A */,
                                  const int Nrhs = 1 /*!< number of right-hand sides, i.e. the number of a's columns */ );
/* *************************************************************************************** */
//! generate n_dir random directions of dimension dim (stored as row vectors in output)
double ** RandomDirections( const int dim,  //!< [in] dimension of the vectors
                            const int n_dir //!< [in] number of directions to generate
                          );

/* *************************************************************************************** */


inline double safeAcos( const double a ) {
    if( a<=-1. ) return PI;
    else
    {
        if ( a>=1. ) return 0.;
        else return acos(a);
    }
}; //!< safe implementation of \c arccos, i.e. data range is truncated to [-1,1]


inline void SolveForQuadraticCoefficients(   const double r0, const double r1, const double r2, // coordinates
                                        const double s0, const double s1, const double s2,      // right-hand side
                                        double &a, double &b, double &c                         // polynomial coefficients
                                    )
{
    a = (r0*r1*s2)/(r0*r1 - r0*r2 - r1*r2 + r2*r2) - (r1*r2*s0)/(r0*r1 + r0*r2 - r1*r2 - r0*r0) - (r0*r2*s1)/(r0*r1 - r0*r2 + r1*r2 - r1*r1);
    b = (s0*(r1 + r2))/(r0*r1 + r0*r2 - r1*r2 - r0*r0) - (s2*(r0 + r1))/(r0*r1 - r0*r2 - r1*r2 + r2*r2) + (s1*(r0 + r2))/(r0*r1 - r0*r2 + r1*r2 - r1*r1);
    c = s2/(r0*r1 - r0*r2 - r1*r2 + r2*r2) - s0/(r0*r1 + r0*r2 - r1*r2 - r0*r0) - s1/(r0*r1 - r0*r2 + r1*r2 - r1*r1);
}



//! sub-namespace for displaying progress of interpolation loop
namespace progress{
    void init(); /*!< initialize progress display*/
    void display(const int current, const int total); /*!< display current progress in 10 % steps*/
};

}; // namespace UTILITY

namespace FastAcos {
    static const int        NSTEPS  = 2048;
    static const double     h       = 2./double(2*NSTEPS-1);
    extern double           data[2*NSTEPS+1];
    extern double           xval[2*NSTEPS];
    extern double           mid_data[2*NSTEPS+1];
    static bool             init;

    void Init();
    double Eval(const double x);
}; // namespace FastAcos



#endif /* _UTIL_H_ */
