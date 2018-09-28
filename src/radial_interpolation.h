#ifndef _RADIAL_INTERPOLATION_H_
#define _RADIAL_INTERPOLATION_H_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <util.h>

/* RadialInterpolation provides a general one-dimensional interpolation
 * framework for a set of data samples.
 * Ideally, the data samples pop-up in the same scheme.
 * Then the interpolation can be replaced by a straight-forward vector-matrix
 * multiplication which is appealing from an algorithmic perspective.
 * 
 * A function that interpolations the data on the same grid is available.
 * 
 * The general RadialInterpolation class provides piecewise linear, continuous
 * interpolation. Other types can be derived from it rather easily by replacing the
 * FIT and INTERPOLATIONWEIGHTS functions.
 * */

class Interpolant;
class InterpolantQuad;

class Interpolant {
protected:
    double  * X;            //!< interval data
    double  * m_D;          //!< matrix containing the data for the interpolation (here: m_D[2*dim*i + 2*j + 0] --> offset, m_D[2*dim*i + 2*j + 1] --> slope (of component j)
    int     dim, n;         //!< dimension of the data; number of data positions
    void    DefaultInit();  //!< initialize data to (secure) default values
    void    Free();         //!< clear memory
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
    //! get interpolation weights (interval index known
    void InterpolationWeights( const int idx, const double x, double * o_w  ) const;
    //! get interpolation weights
    void InterpolationWeights( const double x, double * o_w ) const;
}; /* class Interpolant */




class InterpolantQuad : public Interpolant {
private:
    double  w[3];           //!< weights for internal use
protected:
    void Allocate(const int a_n, const int a_dim ); //!< allocate memory for x and m_D
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
    void TrainC1( int n_x                         /*! [in] number of positions */ ,
                const double * const x          /*! [in] data positions */ ,
                int a_dim                       /*! [in] dimension of the data */ ,
                const double * const m_Data     /*! [in] data j at point x_i --> m_Data[i*a_dim+j]*/
                );

    //! interpolate data at given position x
    void Interpolate( const double x, double * v_Data);
    //! get interpolation weights (interval index known
    void InterpolationWeights( const int idx, const double x, double * o_w  ) const;
    //! get interpolation weights
    void InterpolationWeights( const double x, double * o_w ) const;
}; /* class InterpolantQuad */





#endif /* _RADIAL_INTERPOLATION_H_ */
