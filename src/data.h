#ifndef __DATA_H__
#define __DATA_H__

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

#include "util.h"

using namespace UTILITY;

/*! \brief Container class for general data
 * 
 * This class represents data objects \cD containing coordinates, function values,
 * [optional] gradient values and [optional] Hessians
 * 
 * They can be used as inputs to:
 * 
 * \see \InterpolateOnSet
 * \see \TestOnSet
 * \see \OptimizeGamma
 */
class Data{
public:
    Data();                             //!< default constructor, initializing nullpointers and D=P=-1
    Data(const int a_D, const int a_P); //!< constructor that allocates memory for all pointers
    Data(Data& dat, const bool function_data=false);    //!< copy constructor, always copying coordinates, conditionally copying function data. always allocates memory for all pointers, such that all pointers are ready for use
    ~Data();                            //!< destructor that also frees the memory that was allocated for the pointers
    const int D;          //!< dimension of the coordinates (meant to be const because coordinates and funciton data objects are not designed to be resizable)
    const int P;          //!< number of points (meant to be const because coordinates and funciton data objects are not designed to be resizable)

    // NOTE: in the following we define safe access functions that check the indices for compatibility.
    // this is useful especially during experimental stage, but is inefficient. for maximum efficiency,
    // one could declare the pointers public and operate directly on them
    inline double& SafeAccess_coord(const int p, const int d) {
        assert_msg(p<P && d<D, "ERROR in Data access: out of bound\n"); return coordinates[p][d];}
    inline double* SafeAccess_coord(const int p) {
        assert_msg(p<P       , "ERROR in Data access: out of bound\n"); return coordinates[p];}
    inline double& SafeAccess_values(const int p) {
        assert_msg(p<P,        "ERROR in Data access: out of bound\n"); return values[p];}
    inline double& SafeAccess_gradients(const int p, const int d) {
        assert_msg(p<P && d<D, "ERROR in Data access: out of bound\n"); return gradients[p][d];}
    inline double* SafeAccess_gradients(const int p) {
        assert_msg(p<P       , "ERROR in Data access: out of bound\n"); return gradients[p];}
    inline double& SafeAccess_hessians(const int p, const int component) {
        assert_msg(p<P && component<D*D, "ERROR in Data access: out of bound\n"); return hessians[p][component];}
    inline double* SafeAccess_hessians(const int p) {
        assert_msg(p<P       , "ERROR in Data access: out of bound\n"); return hessians[p];}
    
private:
    double **coordinates;   //!< P-by-D matrix containing the coordinates of the points
    double * values;        //!< array of length P containing the function data values
    double **gradients;     //!< P-by-D matrix containing the funciton data gradients
    double **hessians;      //!< P-by-D*D matrix containing the function data Hessians in symmetric notation
    void initialize();      //!< allocates memory to the pointers
};

/*! \brief Container class for training data
 * 
 * Training data differs from general data in that it
 *   - only stores directions and radii instead of point coordinates
 *   - only stores the radial derivatives instead of gradients
 *   - does not store Hessians
 * 
 * It is assumed that there is one single set of radii for all directions.
 * 
 * In the paper, this is referred to as "training data" but, for the sake of readability.
 * 
 * This can be used as inputs to:
 * 
 * \see \AssembleInterpolation
 */
class DataTraining{
public:
    DataTraining();                                             //!< default constructor, initializing nullpointers and negative constants
    DataTraining(const int a_D, const int a_N, const int a_R);  //!< constructor that allocates memory for all pointers
    DataTraining(DataTraining& dat, bool function_data);        //!< copy constructor, always copying coordinates, conditionally copying function data. always allocates memory for all pointers, such that all pointers are ready for use
    ~DataTraining();                                            //!< destructor that also frees the memory that was allocated for the pointers
    const int D;          //!< dimension of the coordinates (meant to be const because coordinates and funciton data objects are not designed to be resizable)
    const int N;          //!< number of directions
    const int R;          //!< number of radii

    // NOTE: in the following we define safe access functions that check the indices for compatibility.
    // this is useful especially during experimental stage, but is inefficient. for maximum efficiency,
    // one could declare the pointers public and operate directly on them
    inline double& SafeAccess_dir(      const int n, const int d) {
        assert_msg(n<N && d<D, "ERROR in DataTraining access: out of bound\n"); return directions[n][d];}
    inline double* SafeAccess_dir(const int n) {
        assert_msg(n<N       , "ERROR in DataTraining access: out of bound\n"); return directions[n];}
    inline double& SafeAccess_values(   const int n, const int r) {
        assert_msg(n<N && r<R, "ERROR in DataTraining access: out of bound\n"); return values[n][r];}
    inline double* SafeAccess_values(const int n) {
        assert_msg(n<N,        "ERROR in DataTraining access: out of bound\n"); return values[n];}
    inline double& SafeAccess_radderiv( const int n, const int r) {
        assert_msg(n<N && r<R, "ERROR in DataTraining access: out of bound\n"); return radialderiv[n][r];}
    inline double* SafeAccess_radderiv(const int n) {
        assert_msg(n<N       , "ERROR in DataTraining access: out of bound\n"); return radialderiv[n];}
    inline double*     Access_radii() {
                                                                                return radii;}
    inline double& SafeAccess_radii(const int r) {
        assert_msg(r<R       , "ERROR in DataTraining access: out of bound\n"); return radii[r];}
    
private:
    double **directions;    //!< N-by-D matrix containing the unit direction vectors
    double **values;        //!< N-by-R matrix containing the function data values at each radius along each direction
    double **radialderiv;   //!< N-by-R matrix containing the radial derivatives at each radius along each direction
    double * radii;         //!< vector of length R containing the radii
    void initialize();      //!< allocates memory to the pointers
};

#endif
