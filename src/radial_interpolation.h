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

#ifndef _RADIAL_INTERPOLATION_H_
#define _RADIAL_INTERPOLATION_H_

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
 * The interpolants have the form \f$ a+b\cdot x+c\cdot x^2+\dots \f$. The "Vandermonde" array \f$ [1,x,x^2,\dots] \f$ is called \e weights (see Interpolant::weights) and the array \f$ [a,b,c,\dots] \f$ is called \e coefficients (see Interpolant::coeff). The number of terms is stored in Interpolant::N_term.
 *
 * For RadialInterpolation, \e one instance of Interpolant is necessary for RadialInterpolation::Setup. Reminder: throughout this implementation of Concentric Interpolation, it is assumed that the same Interpolant is used for every direction.
 *
 * \see InterpolantQuad
 * \see RadialInterpolation
 *
 */

class Interpolant {
protected:
    double  * Radii;        //!< supporting radii of the piecewise polynomial interpolation, i.e. interval boundaries. contains \p N_rad_supp entries. TODO: including zero?
    int     N_rad_supp;     //!< number of interval boundaries, i.e. length of \p Radii.
    const int N_term;       //!< number of terms of the piecewise polynomials (e.g. 2 for linear, 3 for quadratic, ...), i.e. the order of the polynomial plus one. This could also be called \p N_term. \see Interpolant::coeff, Interpolant::weights
    int     D_val;          //!< number of components of the data, each treated individually in a one-dimensional manner. I.e. simultaneous, individual interpolation of \p D_val one-dimensional functions.
    double  * coeff;        /*!< coefficients \f$ [a,b,c,\dots] \f$ of the polynomial interpolant \f$ a+b\cdot x+c\cdot x^2+\dots \f$. There are Interplant::N_term coefficients at each of the \p N_rad_supp x \p Radii, \e and for each of the \p D_val components. Thus, the length of \p coeff is <tt> N_rad_supp * D_val * N_term </tt>. \verbatim Access: at the i-th point, for the j-th component, the k-th coefficient is
    coeff[i*D_val*N_term + j*N_term + k]. \endverbatim This will be copied to RadialInterpolation upon calling RadialInterpolation::AddData(Interpolant&). \see Interpolant::weights, Train()*/
    void    DefaultInit();  //!< initialize data to (secure) default values
    void    Allocate( const int a_n, const int a_D_val ); //!< allocate memory for x and coeff
    void    Free();         //!< clear memory
private:
    double  weights[2];     //!< Interpolation weights. The piecewise polynomials have the form \f$ a+b\cdot x+c\cdot x^2+\dots \f$. The array \p weights contains the values \f$ [1,x,x^2,\dots] \f$. \see InterpolationWeights(), Interpolant::coeff
public:
    Interpolant(const int a_N_term = 2);
    ~Interpolant();

    //! print the type of the interpolant
    inline void Type(FILE * outstream = stdout) { fprintf(outstream, "linear interpolant\n"); }
    //! set the number of components of the data (implies clearing existing data)
    void CompData( const int d );
    //! return the number of components of the data
    inline int CompData() const { return D_val; }
    //! return the data of the current interval
    inline const double * const InterpolationData( const int idx ) const
    { return coeff + idx*D_val*N_term; }
    //! return \p N_rad_supp
    inline int NIntervals() const { return N_rad_supp; }
    //! return \p N_term
    inline int NTerms() const { return N_term; }
    //! return the interval in which the point is found
    inline int Interval(const double x) const {
        if( x <= Radii[0] )
            return 0; /* take left-most interval and use negative linear extrapolation */
        if( x >= Radii[N_rad_supp-1] )
            return N_rad_supp-1;

        int l=0, r=N_rad_supp-1, m=0;
        while(r-l>3)
        {
            m=int((l+r)/2);
            if( Radii[m] < x )  l = m;
            else            r = m;
        }
        while( (l<r) && (Radii[l+1]<=x)) l++;
        return l;
    }

    //! sets the \p N_rad_supp supporting radii \p Radii, the number \p D_val of components, and fills the coefficients \p coeff \see Interpolant::N_rad_supp, Interpolant::Radii, Interpolant::D_val, Interpolant::coeff
    void Train( const int a_N_rad_supp,            /*!< [in] number of positions */
                const double * const a_Radii,   /*!< [in] data positions */
                const int a_D_val,         /*!< [in] number of components of the data, see Interpolant::D_val*/
                const double * const* m_Data/*!< [in] data component \f$ j \f$ at point \f$ x_i \f$ --> <tt> m_Data[i][j] </tt>*/
                );
    //! interpolate data at given position x
    void Interpolate(   const double x,         /*!< [in] query point (in Concentric Interpolation: radius) */
                        double * v_Data         /*!< [out]  */
                    );
    //! get \p N_term interpolation weights \f$ [1,x,x^2,\dots] \f$ (interval index known)
    void InterpolationWeights( const int idx, const double x, double * o_w  ) const;
    //! get \p N_term  interpolation weights \f$ [1,x,x^2,\dots] \f$
    void InterpolationWeights( const double x, double * o_w ) const;

}; /* class Interpolant */













/** \brief One-dimensional quadratic interpolation. Can be fitted to data either in a C0 manner (calling Train) or in a C1 manner (calling TrainC1).
 *
 * \see Interpolant
 *
 * \attention This has not been tested for the current version.
 */

class InterpolantQuad : public Interpolant {
private:
    double  weights[3];                               //!< Interpolation weights. The piecewise polynomials have the form \f$ a+b\cdot x+c\cdot x^2+\dots \f$. The array \p weights contains the values \f$ [1,x,x^2,\dots] \f$.
public:
    InterpolantQuad();
    ~InterpolantQuad();

    //! print the type of the interpolant
    inline void Type(FILE * outstream = stdout) { fprintf(outstream, "quadratic interpolant\n"); }
    //! initialize the intervals and the data
    void Train( const int n_x,              /*! [in] number of positions */
                const double * const a_Radii,   /*! [in] data positions */
                const int a_D_val,         /*! [in] number of components of the data */
                const double * const m_Data /*! [in] data component \f$ j \f$ at point \f$ x_i \f$ --> <tt> m_Data[i*a_D_val+j] </tt>*/
                );

    //! initialize the intervals and the data, accounting for the C1 property of the data by considering the slopes at the end radii of the intervals
    void TrainC1( const int n_x,            /*! [in] number of positions */
                const double * const a_Radii,   /*! [in] data positions */
                int a_D_val,               /*! [in] number of components of the data */
                const double * const m_Data /*! [in] data component \f$ j \f$ at point \f$ x_i \f$ --> <tt> m_Data[i*a_D_val+j] </tt>*/
                );

    //! interpolate data at given position x
    void Interpolate( const double x, double * v_Data);
    //! get \p N_term = 3 interpolation weights (interval index known)
    void InterpolationWeights( const int idx, const double x, double * o_w  ) const;
    //! get \p N_term = 3 interpolation weights
    void InterpolationWeights( const double x, double * o_w ) const;
}; /* class InterpolantQuad */



















/** \brief This is one of the two main ingredients for
 * ConcentricInterpolation, besides TangentialInterpolation.
 *
 * This interpolates \p D_val dimensional data along each of the \p N_dir
 * directions. Each of the \p D_val data components is interpolated separately
 * by means of a one-dimensional interpolation. To this end, \e one instance of
 * Interpolant is required.
 *
 * To get started, call Setup().
 *
 * It is assumed that ConcentricData is provided, meaning the radii at
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
    int         N_dir;      //!< number of interpolation directions, must equal TangentialInterpolation::N \see RadialInterpolation::N_dir_counter
    int         D_val;     //!< output components of RadialInterpolation
    int         N_int;      //!< number of intervals in the interpolation
    int         N_term;     //!< data size / number of int. weights
    int         N_dir_counter;  //!< counter of the input data: how many interpolants (Interpolant) have been fed to Radial Interpolation. In the language of Concentric Interpolation: counts the directions for which one-dimensional interpolants are included. For each direction, one and only one interpolant shall be added. The ordering of the interpolants is the same as the ordering of the directions in TangentialInterpolation. When calling Interpolate(), \p N_dir_counter must equal RadialInterpolation::N_dir. \see AddData()
    double      * m_I;      //!< support data: big \e pseudo four-dimensional array, initially of size \p N_dir \f$\times\f$ \p N_int \f$\times\f$ \p D_val \f$\times\f$ \p N_term \attention This memory is re-ordered after initial population. \see ReorderData() for indexing information
    double      * m_I_old;  //!< backup of \p m_I before re-ordering \see Setup()
    size_t      m_I_length; //!< length of \p m_I and \p m_I_old
    Interpolant * In;       //!< reference interpolation function in radial direction (in order to get intervals)
    double      * weights;  //!< the weights that are common for all interpolants \see Interpolant::InterpolationWeights()
    void        DefaultInit();
    bool        ordered;    //!< make sure the data is in the right ordering
    bool        initialized;//!< remember if object was initialized
protected:
public:
    RadialInterpolation();
    ~RadialInterpolation();
    void Free();

    //! return the number of components of the outputs
    inline int CompData() const { return D_val; }

    //! same as CompData()
    inline int GetD_val() const { return D_val; }

    //! return the data size (\p D_val * \p N_term) for one radial direction
    inline int DataSize() const { return N_term*D_val; }

    //! return the number of terms of the underlying interpolant class \see Interpolant::N_term
    inline int NWeights() const { return N_term; }

    //! return number of directions, which equals number of inputs. \see RadialInterpolation::N_dir, RadialInterpolation::N_dir_counter, TangentialInterpolation::GetN_dir_supp()
    inline int GetN_dir_supp() const { return N_dir_counter; }

    //! calls SetInterpolant(), Allocate(), Interpolant::Train(), AddData(), and ReorderData() in the correct (this) order
    void Setup( const double * const * const * a_data,  //!< data component \p k along direction \p i at radius \p j --> <tt> a_data[i][j][k] </tt> \todo make const correctness consistent everywhere.
                const int       a_N_dir,    //!< number of directions
                const int       a_n_x,      //!< number of amplitudes
                const int       a_dim_data, //!< number of dimensions of the data
                const double *  a_Radii,        //!< amplitudes
                Interpolant &   a_In,       //!< single interpolant, the same is used for every direction \todo employ user-friendly string identification and automatic creation of the interpolant
                TangentialInterpolation & a_TI //!< TangentialInterpolation
         );
    // TODO: implement overloaded function accepting double*** as first argument

    //! set interpolant for a direction
    void SetInterpolant( Interpolant & a_In );

    //! allocate interpolation data and weights (required In to be set)
    void Allocate( const int a_N_dir  );

    //! add radial data, i.e. the coefficients of the piecewise polynomials along all directions \see Interpolant::coeff, AddData(const double* const)
    void AddData( Interpolant & In );

    //! add radial data, i.e. the coefficients of the piecewise polynomials along all directions \see AddData(Interpolant&)
    void AddData( const double * const coeff );

    /** \brief Multiply inverse kernel matrix from left and re-order
     * interpolation data.
     *
     * Multiply inverse kernel matrix from left to the data backup data w_I_old,
     * store the result to w_I, and re-order using ReorderData() for more
     * efficient memory access during the application (i.e. calls to
     * Interpolate()).
     *
     * This must be called when the kernel parameter \f$ \gamma \f$ was changed,
     * i.e. via TangentialInterpolation::SetGamma(). Always called from Setup().
     */
    void MultiplyKernelMatrixAndReorder( TangentialInterpolation & a_TI );

    /** \brief Re-order the data for efficiency gains.
     *
     * This must be execuded after the inverse of the kernel matrix has been
     * multiplied from the left to \p m_I, meaning that the second index of the
     * matrix is contracted with the first index of \p m_I. After re-ordering,
     * access to \p m_I is very efficient.
     *
     * Recall that the dimensions of RadialInterpolation::m_I are initially
     * \verbatim m_I[N_dir][N_int][D_val][N_term]. \endverbatim
     * and the dimensions of the kernel matrix is \p N_dir \f$ \times \f$ \p N_dir.
     * So, the formal multiplication \f$ \underline{\underline{K}}^{-1} \f$ \p m_I
     * does not change the dimension.
     *
     * After this multiplication, this function must be called in oder to swap
     * the first two indices, such that the dimensions of \p m_I then are
     * \verbatim m_I[i_int][i_dir][i_comp][i_coeff]. \endverbatim
     */
    void ReorderData(const bool enforce_order = false   /**< [in] re-orders even if \p ordered == true*/ );


    /** \brief interpolate data: return a \e pseudo matrix of size \p N_dir \f$\times\f$ \p D_val in row-major format, representing \p D_val column vectors \f$ \underline{\mathcal{R}}(x) \f$. See equation (CI) in the paper.
     *
     * \todo update the following. TODO
     *
     * The \e pseudo four-dimensional data array \p RadialInterpolation::m_I has the indexing scheme
     * \verbatim m_I[i_int][i_dir][i_comp][i_coeff] \endverbatim
     * after a call of ReorderData(). The first dimension \p i_int is specified by the amplitude \p a_x via a call of Interpolant::Interval(). The last index is contracted with the weight vector \f$ [1,x,x^2,\dots] \f$ returned by Interpolant::InterpolationWeights().
     */
    void Interpolate( const double a_x /*! [in] radial position \p x for interpolation */ ,
                      double * m_out   /*! [out] data component \p j at radius \p x along direction \p i --> <tt> m_out[i*D_val+j]</tt>. In other words, the output is a \e pseudo matrix of size \p N_dir x \p D_val */ ) ;

    //! returns an immutable pointer to the (private) data
    const double * const Data() const
    { return m_I; }
};



#endif /* _RADIAL_INTERPOLATION_H_ */
