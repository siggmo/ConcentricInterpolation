#ifndef _QUADRATIC_INTERPOLANT_H_
#define _QUADRATIC_INTERPOLANT_H_

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

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <util.h>

using namespace UTILITY;

class QuadraticInterpolant;

/** TODO comment
 * */
class QuadraticInterpolant {
private:
    int     num_support_r;    //!< number of supporting radii, including zero. this is (num_training_radii-1)/2 in the main function
    double  * poly_a, * poly_b, * poly_c, //!< polynomial coefficient containters, f(r) = a + b*r + c*r*r, piecewisely defined coefficients
            * support_r;      //!< supporting radii, including zero. this is a subset (i.e. indices 0, 2, 4, ... excluding the last) of training_radii in the main function \see SetData
    bool    init;             //!< true, if it has been initialized

    void DefaultInit();       //!< default initialization, sets all member variables zero and false
    void Allocate();          //!< allocates memory for member pointers
    void Free();              //!< frees memory allocated for member pointers
    
public:
    QuadraticInterpolant();
    QuadraticInterpolant(   const int a_num_support_r,      //!<[in] number of supporting radii
                        const double * a_support_r, //!<[in] array of supporting radii
                        const double * S   //!<[in] array of values at the supporting radii
                    );
    ~QuadraticInterpolant(); //!< calls Free()
    
    void SetData(       const int a_num_training_r,  //!<[in] number of training radii, possibly not containing 0
                        const double * a_training_r, //!<[in] array of training radii, possibly not containing 0
                        const double * S           //!<[in] array of values at the supporting points
                );                                          //!< sets up the function by computing the polynomial coefficients
    
    // access functions (const-qualified)
    double Value( const double r    //!<[in] radius at which to evaluate the interpolant
                ) const;            //!< returns the value of the piecewise cubic function at r
    double Interpolate( const double r, //!<[in] radius at which to evaluate the interpolant and its first two derivatives
                        double & dS,    //!<[out] first derivative
                        double & ddS    //!<[out] second derivative
                ) const; //!< same as Value() but additionally computes the first and second derivative
    void FirstV( double & v0, double & r ) const;
    void FirstDV( double & dv0, double & r ) const;
    void GetSupport( double * a_support_r, int & a_num_support_r) const;        //!< copies support_r and num_support_r to the respective arguments
};



#endif /* _QUADRATIC_INTERPOLANT_H_ */
