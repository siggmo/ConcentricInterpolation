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

#ifndef _UTIL_H_
#define _UTIL_H_

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <math.h>
#include <string.h>
#include <unistd.h>
#include <sys/stat.h>
#include <random>
#include <chrono>

#ifdef CONINTER_MKL
    #include <mkl_lapacke.h>
    #include <mkl_cblas.h>
#else
    #include <lapacke.h>
    #include <cblas.h>
#endif /* CONINTER_MKL */

#ifndef PI
#define PI 3.14159265358979
#endif

#define DEFAULT_ERROR_CODE 10

//! Handy utility functions.
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
                         const int NMAX_LINES=2000, const int NMAX_COL=1000 ); // ATTENTION: if this function segfaults or fails in some other way, try increasing NMAX_LINES and NMAX_COL.
//! dump a pseudo-2d-matrix or an array to a stream. if array, set m = 1.
void        print_matrix(   const double * const a/*!<[in] 1-d or pseudo 2-d pointer*/,
                            const int m/*!<[in] number of rows*/,
                            const int n/*!<[in] number of columns*/,
                            FILE * F = stdout/*!<[in] stream pointer*/);
//! nicely print the upper triangular part of a symmetric pseudo-2d-matrix
void        print_matrix_sym( const double * const a/*!<[in] pseudo 2-d pointer*/,
                            const int m/*!<[in] number of rows = number of cols of symmetric matrix*/,
                            FILE * F = stdout/*!<[in] stream pointer*/);
//! simple assert command: prints str to outputstream if condition is false. if quit is true, it then also exits with DEFAULT_ERROR_CODE
void        assert_msg(const bool condition, const char *str, const bool quit=true, FILE * outputstream=stderr );


/* *************************************************************************************** */
/* memory management */
double *    alloc_array( const size_t N )  ; //!< allocate array of length \c N (N is the number of doubles to allocate memory for)
void        free_array( double ** )  ; //!< savely free previously allocated array

//! \brief function for allocating pseudo 2-dimensional arrays. memory is filled with zeros.
//! \details \see free_matrix
double ** alloc_matrix( const size_t m, const size_t n );

//! \brief function for allocating pseudo 3-dimensional arrays. memory is filled with zeros.
//! \details \see free_array3
double *** alloc_array3( const size_t m, const size_t n, const size_t o );

//! \brief function for freeing memory allocated by alloc_matrix
//! \details \see alloc_matrix
void free_matrix( double ** &A, const size_t  m /**< [in] first dimension of A */ );

//! \brief function for freeing memory allocated by alloc_array3
//! \details \see alloc_matrix
void free_array3( double *** &A, const size_t  m,  /**< [in] first dimension of A */ const size_t n /**< [in] second dimension of A */ );
/* *************************************************************************************** */
double      norm( const double * a, const int dim )  ; //!< compute the \f$l^2\f$ vector norm of \c a in R^\c dim
double      distance( const double * a, const double * b, const int dim )  ; //!< compute the \f$l^2\f$ vector distance of vector \c a and \c b
void        scale( const double * vec1, const double val, double * vec2, const int dim );//!< scales \p vec1 by \p val and stores result to \p vec2. assumes both \p vec1 and \p vec2 are of length \p dim.

/* *************************************************************************************** */
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
/** \brief on exit, \f$y_i=A_{ij}x_j\f$ or (if transposition flag \c transpose==true) \f$y_j=A_{ij}x_i\f$.
*
* For use in ConcentricInterpolation, transposition should always be true.
*
* \see ConcentricInterpolation::Evaluate()
* ConcentricInterpolation::Error()
*/
inline void MatVecMul( const double * A, const double * x, double * y, const int M, const int N, const bool transpose = false )
{
    if( transpose )
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
/* use LDL solver to compute the product of the inverse kernel matrix and a vector */
//! stores LDL factorization of \c A to \c Af, and corresponding interchange indices to the integer array \c w_i. \c A and \c Af need to be of size \c N*\c N and the size of \c w_i should be at least the same size or bigger
void        Factorize(  const double * A /*!< [in] N*N matrix*/,
                        double * Af      /*!< [out] LDL factorization of A */,
                        int * w_i        /*!< [out] integer working array of size >= N*N, this is \p IPIV in LAPACK's \p dsytrf*/,
                        const int N      /*!< [in] dimension of A */);
//! solves the system \p A \f$\cdot\f$ \p Ai_a \f$=\f$ \p a for \c Ai_a by the LDL factorization \c Af. The factorization must have been computed by Factorize().
void        SolveByFactorization( const double * Af  /*!< [in] LDL factorization of \p A; see Factorize()*/,
                                  const double * a   /*!< [in] right hand side vector (or matrix) */,
                                  double * Ai_a      /*!< [out] solution of the linear system \p A \f$\cdot\f$ \p Ai_a \f$=\f$ \p a */,
                                  int * w_i          /*!< [in] permutation vector; see Factorize() */,
                                  const int N        /*!< [in] dimension of A */,
                                  const int Nrhs = 1 /*!< [in] number of right-hand sides, i.e. the number of a's columns */ );
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
