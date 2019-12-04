#ifndef __DATA_H__
#define __DATA_H__

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
 *     Oliver Kunc and Felix Fritzen: 'Generation of energy-minimizing point
 *                                     sets on spheres and their application in
 *                                     mesh-free interpolation and
 *                                     differentiation'
 *     Advances in Computational Mathematics, Number/Volume, p. XX-YY, 2019
 *     DOI   10.1007/s10444-019-09726-5
 *     URL   dx.doi.org/10.1007/s10444-019-09726-5
 *
 *  The latest version of this software can be obtained through
 *  https://github.com/EMMA-Group/ConcentricInterpolation
 *
 */

#include "util.h"

using namespace UTILITY;

/*! \brief Container for General Data. "General Data" means possibly irregular
 * coordinates together with their corresponding function values.
 *
 * Usage: for evaluation of Concentric Interpolation, error computation, ...
 *
 * The difference to ConcentricData is that, here, the coordinates are not assumed to be concentric. Also, their format is
 * absolute, meaning that they are taken as-is, in contrast to the amplitude-direction-split within ConcentricData.
 *
 * \todo investigate suitability of certain container classes as possible substitutes for the low-level pointer objects in
 * order to maximize user-friendliness.
 */
class GeneralData{
private:
    double**    Coords;     //!< \p N_pts -by-\p D_inp matrix containing the coordinates corresponding to \p Values
    double**    Values;     //!< \p N_pts -by-\p D_val matrix containing the function values corresponding to \p Coords
    int         D_inp;      //!< dimension of the interpolation's input artuments
    int         D_val;      //!< number of scalar values that are interpolated simultaneously
    int         N_pts;      //!< number of data points, i.e. number of coordinates and of values
    bool        initialized;//!< true if memory has been allocated
public:
    //! Constructor initializing nullpointers and zero sizes \p D_inp, \p D_val, \p N_pts.
    GeneralData();
    //! Constructor allocating memory for coordinates and values.
    GeneralData(
        const int a_D_inp,
        const int a_D_val,
        const int a_N_pts
        )
    { initialized = false, Resize(a_D_inp, a_D_val, a_N_pts); }

    ~GeneralData(){ Free(); }                  //!< destructor that also frees the memory that was allocated for the pointers

    inline int  GetD_inp() const { return D_inp; }
    inline int  GetD_val() const { return D_val; }
    inline int  GetN_pts() const        { return N_pts; }
    inline bool IsInitialized() const { return initialized; }

    void Free();                        //!< frees memory and sets nullpointers
    void Resize(
        const int a_D_inp,
        const int a_D_val,
        const int a_N_pts
         );

    inline const double& GetCoord(const int n_point, const int d_inp) const { return Coords[n_point][d_inp]; }  //!< read-only access to single coordinate component of point \p n_point
    inline const double& GetValue(const int n_point, const int d_val) const { return Values[n_point][d_val]; }  //!< read-only access to single value at point \p n_point
    inline const double * GetCoords(const int n_point)    const { return Coords[n_point]; } //!< read-only access to coordinates of point \p n_point
    inline const double * GetValues(const int n_point)    const { return Values[n_point]; } //!< read-only access to values at point \p n_point

    inline double& SetCoord(const int n_point, const int d_inp) { return Coords[n_point][d_inp]; }    //!< write access to coordinates of point \p n_point
    inline double& SetValue(const int n_point, const int d_val) { return Values[n_point][d_val]; }    //!< write access to values at point \p n_point
    inline double* SetCoords(const int n_point)   { return Coords[n_point]; }   //!< write access to coordinates of point \p n_point
    inline double* SetValues(const int n_point)   { return Values[n_point]; }   //!< write access to values at point \p n_point
};












/*! \brief Container for Concentric Data. "Concentric Data" means directions, amplitudes, and the corresponding function values.
*
* Strictly speaking, the data points may be positioned irregularly. No check of
* "concentricity" or whatsoever is performed. But for ConcentricInterpolation,
* it is assumed that the data is indeed concentric.
*
* Usage: for setup of ConcentricInterpolation.
*
* The difference to GeneralData is that, here, the coordinates are given in the direction-amplitude-split notation. This is essential
* for the setup of Concentric Interpolation.
*/
    class ConcentricData{
    private:
        double**    Directions; //!< \p N_dir -by-\p D_inp matrix containing the coordinates corresponding to \p Values
        double*     Radii;      //!< \p N_rad -sized array containing the radii of the support points \attention must contain zero as first entry
        double***   Values;     //!< \p N_dir -by-\p N_rad -by-\p D_val array: access a scalar function value at a radius along a direction by <tt>Values[n_dir][n_rad][d_val]</tt>
        int         D_inp;   //!< dimension of the interpolation's input artuments
        int         D_val;   //!< number of scalar values that are interpolated simultaneously
        int         N_dir;          //!< number of directions
        int         N_rad;          //!< number of radii
        bool        initialized;//!< true if memory has been allocated
    public:
        //! Constructor initializing nullpointers and zero sizes \p D_inp, \p D_val, \p N_dir, \p N_rad.
        ConcentricData();
        //! Constructor allocating memory for coordinates and values.
        ConcentricData(
            const int a_DimDirections,
            const int a_D_val,
            const int a_N_dir,
            const int a_N_rad
            )
        { initialized = false, Resize(a_DimDirections, a_D_val, a_N_dir, a_N_rad); }
        //! Constructor allocating memory and filling that memory with data from textfiles
        ConcentricData(
            char * a_FilenameSupportDirections, //!< name of file containing directions as rows
            char * a_FilenameSupportRadii,      //!< name of file containing radii in one row
            char * a_FilenameSupportValues,     //!< name of file containing the function values in the following format: one row per direction. in each row, there are \p D_val many entries for the first radius, then \p D_val many entries for the second radius, and so on.
            const int a_D_val = 1               //!< number of values per support point
        );

        ~ConcentricData(){ Free(); }                  //!< frees memory and sets nullpointers

        inline int  GetD_inp() const  { return D_inp; }
        inline int  GetD_val() const  { return D_val; }
        inline int  GetN_dir() const  { return N_dir; }
        inline int  GetN_rad() const  { return N_rad; }
        inline bool IsInitialized() const { return initialized; }

        void Free();
        void Resize(
            const int a_D_inp,
            const int a_D_val,
            const int a_N_dir,
            const int a_N_rad
            );

        inline const double& GetDirection(const int n, const int d_inp) const   { return Directions[n][d_inp]; }//!< read-only access to single direction component
        inline const double * const * GetDirections() const                     { return Directions; }          //!< read-only access to all directions
        inline const double& GetRadius(const int r) const                       { return Radii[r]; }            //!< read-only access to <tt>r</tt>-th radius
        inline const double * GetRadii() const                                  { return Radii; }               //!< read-only access to all radii
        inline const double& GetValue(const int n, const int r, const int d_val) const { return Values[n][r][d_val]; }  //!< read-only access to single value component
        inline const double * const * const * GetValues() const                 { return Values; }              //!< read-only access to all values

        inline double& SetDirection(const int n, const int d_inp)               { return Directions[n][d_inp]; }//!< write access to single direction component
        inline double ** SetDirections()                                        { return Directions; }          //!< write access to all directions
        inline double& SetRadius(const int r)                                   { return Radii[r]; }            //!< write access to <tt>r</tt>-th radius
        inline double * SetRadii()                                              { return Radii; }               //!< write access to all radii
        inline double& SetValue(const int n, const int r, const int d_val)      { return Values[n][r][d_val]; } //!< write access to single value component
        inline double *** SetValues()                                           { return Values; }              //!< write access to all values

    };
    #endif
