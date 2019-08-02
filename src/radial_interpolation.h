#ifndef _RADIAL_INTERPOLATION_H_
#define _RADIAL_INTERPOLATION_H_

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
#include <math.h>
#include <util.h>

class Interpolant;
class InterpolantQuad;
class RadialInterpolation;

/** \brief Class for linear interpolation, serving as parent to all other
 *
 */

class Interpolant {
protected:
    double  * X;            //!< interval data
    double  * m_D;          //!< matrix containing the data for the interpolation (here: m_D[2*dim*i + 2*j + 0] --> offset, m_D[2*dim*i + 2*j + 1] --> slope (of component j)
    int     dim, n;         //!< dimension of the data; number of data positions
    void    DefaultInit();  //!< initialize data to (secure) default values
    void    Free();         //!< clear memory
    int     n_int_data;     //!< size of interpolation data (2 for linear, 3 for quadratic, ...)
private:
    double  w[2];           //!< weights for internal use
    void    Allocate( const int a_n, const int a_dim ); //!< allocate memory for x and m_D
public:
    Interpolant();
    ~Interpolant();

    void    Dim( const int d );    //!< set the dimension of the data (implies clearing existing data)
    int     Dim() const;            //!< get the dimension of the data
    //! return the interval in which the point is found
    inline int Interval(const double x) const {
        if( x <= X[0] )
            return 0; /* take left-most interval and use negative linear extrapolation */
        if( x >= X[n-1] )
            return n-1;

        int l=0, r=n-1, m=0;
        while(r-l>3)
        {
            m=int((l+r)/2);
//             printf("x[l] %8.4f   x[m] %8.4f   x[r] %8.4f   x    %8.4f\n", X[l], X[m], X[r], x);
            if( X[m] < x )  l = m;
            else            r = m;
        }
        while( (l<r) && (X[l+1]<=x)) l++;
//         printf("x_- %8.4f  x %8.4f  X_+ %8.4f\n", X[l], x, X[l+1]);
        return l;
    }

    //! initialize the intervals and the data
    void Train( int n_x                         /*! [in] number of positions */ ,
                const double * const x          /*! [in] data positions */ ,
                int a_dim                       /*! [in] dimension of the data */ ,
                const double * const m_Data     /*! [in] data j at point x_i --> m_Data[i*a_dim+j]*/
                );
    //! interpolate data at given position x
    void Interpolate( const double x, double * v_Data);
    //! get interpolation weights (interval index known)
    void InterpolationWeights( const int idx, const double x, double * o_w  ) const;
    //! get interpolation weights
    void InterpolationWeights( const double x, double * o_w ) const;
    //! return the data of the current interval
    const double * const InterpolationData( const int idx ) const;
    const int NIntervals() const;   //!< return number of data intervals
    const int DataSize() const; //!< return size interpolation data (dim * n_int_data)
}; /* class Interpolant */




class InterpolantQuad : public Interpolant {
private:
    double  w[3];           //!< weights for internal use
protected:
    void Allocate(const int a_n, const int a_dim ); //!< allocate memory for x and m_D
    void    DefaultInit();  //!< initialize data to (secure) default values
public:
    InterpolantQuad();
    ~InterpolantQuad();

    //! initialize the intervals and the data
    void Train( int n_x                         /*! [in] number of positions */ ,
                const double * const x          /*! [in] data positions */ ,
                int a_dim                       /*! [in] dimension of the data */ ,
                const double * const m_Data     /*! [in] data j at point x_i --> m_Data[i*a_dim+j]*/
                );

    //! initialize the intervals and the data (using C1 continuous data)
    void TrainC1( int n_x                       /*! [in] number of positions */ ,
                const double * const x          /*! [in] data positions */ ,
                int a_dim                       /*! [in] dimension of the data */ ,
                const double * const m_Data     /*! [in] data j at point x_i --> m_Data[i*a_dim+j]*/
                );

    //! interpolate data at given position x
    void Interpolate( const double x, double * v_Data);
    //! get interpolation weights (interval index known)
    void InterpolationWeights( const int idx, const double x, double * o_w  ) const;
    //! get interpolation weights
    void InterpolationWeights( const double x, double * o_w ) const;
}; /* class InterpolantQuad */





/** RadialInterpolation provides a general one-dimensional interpolation
 * framework for a set of data samples.
 * Ideally, the data samples pop-up in the same scheme for all interpolants,
 * e.g. in the context of Concentric Interpolation: along all directions.
 * Then the interpolation can be replaced by a straight-forward vector-matrix
 * multiplication which is appealing from an algorithmic perspective.
 *
 * A function that interpolates the data on the same grid is available.
 *
 * The general RadialInterpolation class provides piecewise linear, continuous
 * interpolation. Other types can be derived from it rather easily by replacing
 * the FIT and INTERPOLATIONWEIGHTS functions.
 *
 * For Concentric Interpolation, this is used in combination with Tangential
 * Interpolation.
 *
 * \see TangentialInterpolation
 * */

class RadialInterpolation {
private:
    int         n_dir;      //!< number of interpolation directions
    int         dim;        //!< output dimension of RadialInterpolation
    int         n_int;      //!< number of intervals in the interpolation
    int         n_w;        //!< data size / number of int. weights
    int         n_input;    //!< counter of the input data (i.e. how many RI are actually fed into the object)
    double      * m_I;      //!< interpolation data (big matrix of size n_dir x n_interval x dim x n_w)
    Interpolant * In;       //!< reference interpolation function in radial direction (in order to get intervals)
    double      * w;        //!< the weights that are common for all interpolants
    void        DefaultInit();
    bool        ordered;    //!< make sure the data is in the right ordering
protected:
public:
    RadialInterpolation();
    ~RadialInterpolation();
    void Free();

    //! return the dimension of the outputs
    int Dim() const;

    //! return the data size (dim*n_w) for one radial direction
    int DataSize() const;

    //! return the number of weights of the underlying interpolant class
    int NWeights() const;

    //! return number of directions fed into the object
    int NInput() const;

    //! To be done first: set interpolant
    void SetInterpolant( Interpolant * a_In );

    //! allocate interpolation data and weights (required In to be set)
    void Allocate( const int a_n_dir  );

    //! add radial data
    void AddData( Interpolant * In );

    //! add radial data
    void AddData( const double * const m_D );

    //! reorder the data such that Interpolate is indeed efficient
    void ReorderData();


    //! interpolate data
    void Interpolate( const double x /*! [in] x position for interpolation */ ,
                      double * m_out /*! [out] interpolation data m_out[i*dim+j] --> j-th component in direction i */ ) ;

    const double * const Data() const;
};



#endif /* _RADIAL_INTERPOLATION_H_ */
