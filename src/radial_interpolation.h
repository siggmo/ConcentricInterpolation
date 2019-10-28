#ifndef _RADIAL_INTERPOLATION_H_
#define _RADIAL_INTERPOLATION_H_

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


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <util.h>
#include "tangential_interpolation.h"

class Interpolant;
class InterpolantQuad;
class RadialInterpolation;

/** \brief Parent class for one-dimensional interpolation. This performs linear interpolation, children may perform higher-order interpolation.
 *
 * The interpolants have the form \f$ a\cdot 1+b\cdot x+c\cdot x^2+\dots \f$. The "Vandermonde" array \f$ [1,x,x^2,\dots] \f$ is called \e weights (see Interpolant::weights) and the array \f$ [a,b,c,\dots] \f$ is called \e coefficients (see Interpolant::coeff). The number of terms is stored in Interpolant::n_term.
 *
 * For RadialInterpolation, \e one instance of Interpolant is necessary for RadialInterpolation::Setup. Reminder: throughout this implementation of Concentric Interpolation, it is assumed that the same Interpolant is used for every direction.
 *
 * \see InterpolantQuad
 * \see RadialInterpolation
 *
 */

class Interpolant {
protected:
    double  * X;            //!< supporting points of the piecewise polynomial interpolation, i.e. interval boundaries. contains \p n_X entries. TODO: including zero?
    int     n_X;            //!< number of interval boundaries, i.e. length of \p X.
    const int n_term;       //!< number of terms of the piecewise polynomials (e.g. 2 for linear, 3 for quadratic, ...), i.e. the order of the polynomial plus one. This could also be called \p n_coeff. \see Interpolant::coeff, Interpolant::weights
    int     n_comp;         //!< number of components of the data, each treated individually in a one-dimensional manner. I.e. simultaneous, individual interpolation of \p n_comp one-dimensional functions.
    double  * coeff;        /*!< coefficients \f$ [a,b,c,\dots] \f$ of the polynomial interpolant \f$ a\cdot 1+b\cdot x+c\cdot x^2+\dots \f$. There are Interplant::n_term<tt> = n_coeff </tt> coefficients at each of the \p n_X points \p X, \e and for each of the \p n_comp components. Thus, the length of \p coeff is <tt> n_X * n_comp * n_coeff </tt>. \verbatim Access: at the i-th point, for the j-th component, the k-th coefficient is
    coeff[i*n_comp*n_coeff + j*n_coeff + k]. \endverbatim This will be copied to RadialInterpolation upon calling RadialInterpolation::AddData(Interpolant&). \see Interpolant::weights, Train()*/
    void    DefaultInit();  //!< initialize data to (secure) default values
    void    Allocate( const int a_n, const int a_n_comp ); //!< allocate memory for x and coeff
    void    Free();         //!< clear memory
private:
    double  weights[2];     //!< Interpolation weights. The piecewise polynomials have the form \f$ a\cdot 1+b\cdot x+c\cdot x^2+\dots \f$. The array \p weights contains the values \f$ [1,x,x^2,\dots] \f$. \see InterpolationWeights(), Interpolant::coeff
public:
    Interpolant(const int a_n_term = 2);
    ~Interpolant();

    //! print the type of the interpolant
    inline void Type(FILE * outstream = stdout) { fprintf(outstream, "linear interpolant\n"); }
    //! set the number of components of the data (implies clearing existing data)
    void CompData( const int d );
    //! return the number of components of the data
    inline int CompData() const { return n_comp; }
    //! return the data of the current interval
    inline const double * const InterpolationData( const int idx ) const
    { return coeff + idx*n_comp*n_term; }
    //! return \p n_X
    inline int NIntervals() const { return n_X; }
    //! return \p n_term
    inline int NTerms() const { return n_term; }
    //! return the interval in which the point is found
    inline int Interval(const double x) const {
        if( x <= X[0] )
            return 0; /* take left-most interval and use negative linear extrapolation */
        if( x >= X[n_X-1] )
            return n_X-1;

        int l=0, r=n_X-1, m=0;
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

    //! sets the \p n_X supporting points \p X, the number \p n_comp of components, and fills the coefficients \p coeff \see Interpolant::n_X, Interpolant::X, Interpolant::n_comp, Interpolant::coeff
    void Train( const int a_n_X,                /*!< [in] number of positions */
                const double * const a_X,       /*!< [in] data positions */
                const int a_n_comp,             /*!< [in] number of components of the data, see Interpolant::n_comp*/
                const double * const m_Data     /*!< [in] data component \f$ j \f$ at point \f$ x_i \f$ --> <tt> m_Data[i*a_n_comp+j] </tt>*/
                );
    //! interpolate data at given position x
    void Interpolate(   const double x,         /*!< [in] query point (in Concentric Interpolation: radius) */
                        double * v_Data         /*!< [out]  */
                    );
    //! get \p n_term interpolation weights \f$ [1,x,x^2,\dots] \f$ (interval index known)
    void InterpolationWeights( const int idx, const double x, double * o_w  ) const;
    //! get \p n_term  interpolation weights \f$ [1,x,x^2,\dots] \f$
    void InterpolationWeights( const double x, double * o_w ) const;

}; /* class Interpolant */













/** \brief One-dimensional quadratic interpolation. Can be fitted to data either in a C0 manner (calling Train) or in a C1 manner (calling TrainC1).
 *
 * \see Interpolant
 */

class InterpolantQuad : public Interpolant {
private:
    double  weights[3];                               //!< Interpolation weights. The piecewise polynomials have the form \f$ a\cdot 1+b\cdot x+c\cdot x^2+\dots \f$. The array \p weights contains the values \f$ [1,x,x^2,\dots] \f$.
public:
    InterpolantQuad();
    ~InterpolantQuad();

    //! print the type of the interpolant
    inline void Type(FILE * outstream = stdout) { fprintf(outstream, "quadratic interpolant\n"); }
    //! initialize the intervals and the data
    void Train( const int n_x,              /*! [in] number of positions */
                const double * const a_X,   /*! [in] data positions */
                const int a_n_comp,         /*! [in] number of components of the data */
                const double * const m_Data /*! [in] data component \f$ j \f$ at point \f$ x_i \f$ --> <tt> m_Data[i*a_n_comp+j] </tt>*/
                );

    //! initialize the intervals and the data, accounting for the C1 property of the data by considering the slopes at the end points of the intervals
    void TrainC1( const int n_x,            /*! [in] number of positions */
                const double * const a_X,   /*! [in] data positions */
                int a_n_comp,               /*! [in] number of components of the data */
                const double * const m_Data /*! [in] data component \f$ j \f$ at point \f$ x_i \f$ --> <tt> m_Data[i*a_n_comp+j] </tt>*/
                );

    //! interpolate data at given position x
    void Interpolate( const double x, double * v_Data);
    //! get \p n_term = 3 interpolation weights (interval index known)
    void InterpolationWeights( const int idx, const double x, double * o_w  ) const;
    //! get \p n_term = 3 interpolation weights
    void InterpolationWeights( const double x, double * o_w ) const;
}; /* class InterpolantQuad */



















/** RadialInterpolation is one of the two main ingredients for
 * ConcentricInterpolation, besides TangentialInterpolation. It interpolates
 * \p n_comp dimensional data along each of the \p n_dir directions. Each of
 * the \p n_comp data components is interpolated separately by means of a
 * one-dimensional interpolation. To this end, \e one instance of Interpolant is
 * required.
 *
 * To get started, call Setup().
 *
 * It is assumed that \e Concentric \e Data is provided, meaning the \e radii at
 * which the data is provided <em>are the same for all directions</em>. While
 * this assumption is detrimental to generality, it provides significant
 * potential for efficiency optimizations. These optimizations are realized in
 * the current implementation.
 *
 * In essence, the interpolation is replaced by a <em>single vector-matrix
 * multiplication</em>. TODO
 *
 *
 * \see TangentialInterpolation
 * \see ConcentricInterpolation
 * */

class RadialInterpolation {
private:
    int         n_dir;      //!< number of interpolation directions, must equal TangentialInterpolation::N \see RadialInterpolation::n_input
    int         n_comp;     //!< output components of RadialInterpolation
    int         n_int;      //!< number of intervals in the interpolation
    int         n_term;     //!< data size / number of int. weights
    int         n_input;    //!< counter of the input data: how many interpolants (Interpolant) have been fed to Radial Interpolation. In the language of Concentric Interpolation: counts the directions for which one-dimensional interpolants are included. For each direction, one and only one interpolant shall be added. The ordering of the interpolants is the same as the ordering of the directions in TangentialInterpolation. When calling Interpolate(), \p n_input must equal RadialInterpolation::n_dir. \see AddData()
    double      * m_I;      //!< interpolation data: big \e pseudo four-dimensional array of size \p n_dir \f$\times\f$ \p n_int \f$\times\f$ \p n_comp \f$\times\f$ \p n_term \attention This memory is re-ordered after initial population. \see ReorderData() for indexing information
    Interpolant * In;       //!< reference interpolation function in radial direction (in order to get intervals)
    double      * weights;  //!< the weights that are common for all interpolants \see Interpolant::InterpolationWeights()
    void        DefaultInit();
    bool        ordered;    //!< make sure the data is in the right ordering
protected:
public:
    RadialInterpolation();
    ~RadialInterpolation();
    void Free();

    //! return the number of components of the outputs
    inline int CompData() const { return n_comp; }

    //! return the data size (\p n_comp * \p n_term) for one radial direction
    inline int DataSize() const { return n_term*n_comp; }

    //! return the number of terms of the underlying interpolant class \see Interpolant::n_term
    inline int NWeights() const { return n_term; }

    //! return number of directions, which equals number of inputs. \see RadialInterpolation::n_dir, RadialInterpolation::n_input, TangentialInterpolation::GetN()
    inline int GetN() const { return n_input; }

    //! calls SetInterpolant(), Allocate(), Interpolant::Train(), AddData(), and ReorderData() in the correct order
    void Setup( const double * const * a_data,     //!< data component k along direction i at radius j --> <tt> a_data[i][j*a_dim_data+k] </tt>
                const int       a_n_dir,    //!< number of directions
                const int       a_n_x,      //!< number of amplitudes
                const int       a_dim_data, //!< number of dimensions of the data
                const double *  X,          //!< amplitudes
                Interpolant &   a_In,       //!< single interpolant, the same is used for every direction TODO: employ string and automatic creation of the interpolant
                TangentialInterpolation & a_TI //!< TangentialInterpolation
         );

    //! set interpolant for a direction
    void SetInterpolant( Interpolant & a_In );

    //! allocate interpolation data and weights (required In to be set)
    void Allocate( const int a_n_dir  );

    //! add radial data, i.e. the coefficients of the piecewise polynomials along all directions \see Interpolant::coeff, AddData(const double* const)
    void AddData( Interpolant & In );

    //! add radial data, i.e. the coefficients of the piecewise polynomials along all directions \see AddData(Interpolant&)
    void AddData( const double * const coeff );

    /** Re-order the data for efficiency gains. This must
     * be execuded after the final call of AddData(), and is done in Setup().
     *
     * Let \p i_dir, \p i_int, \p i_comp, and \p i_coeff denote the index of the
     * direction, the interval, the component, and the coefficient (i.e. term of
     * the polynomial), respectively. The row-major format is used to store
     * <center> \p n_dir \f$\times\f$ \p n_int \f$\times\f$ \p n_comp \f$\times\f$ \p n_term </center>
     * values in RadialInterpolation::m_I, which is a \e pseudo four-dimensional array. Recall that \p n_term = number of coefficients, as in Interpolant::n_term.
     *
     * Before this function is invoked, the \e pseudo indexing scheme is
     * \verbatim m_I[i_dir][i_int][i_comp][i_coeff]. \endverbatim
     * After this function has been executed, the first two indices are swapped:
     * \verbatim m_I[i_int][i_dir][i_comp][i_coeff]. \endverbatim
     */
    void ReorderData();


    /** \brief interpolate data: return a \e pseudo matrix of size \p n_dir \f$\times\f$ \p n_comp in row-major format, representing \p n_comp column vectors \f$ \underline{\mathcal{R}}(x) \f$. See equation (CI) in the paper.
     *
     * The \e pseudo four-dimensional data array \p RadialInterpolation::m_I has the indexing scheme
     * \verbatim m_I[i_int][i_dir][i_comp][i_coeff] \endverbatim
     * after a call of ReorderData(). The first dimension \p i_int is specified by the amplitude \p a_x via a call of Interpolant::Interval(). The last index is contracted with the weight vector \f$ [1,x,x^2,\dots] \f$ returned by Interpolant::InterpolationWeights().
     */
    void Interpolate( const double a_x /*! [in] x position for interpolation */ ,
                      double * m_out   /*! [out] data component \p j at radius \p x along direction \p i --> <tt> m_out[i*n_comp+j]</tt>. In other words, the output is a \e pseudo matrix of size \p n_dir x \p n_comp */ ) ;

    //! returns an immutable pointer to the (private) data
    const double * const Data() const
    { return m_I; }
};



#endif /* _RADIAL_INTERPOLATION_H_ */
