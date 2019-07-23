#ifndef __DISTANCE_FUNCTIONS_H__
#define __DISTANCE_FUNCTIONS_H__

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

#include <cmath>
#include <stdio.h>

const double DistanceFunction_TOL = 1.e-5;

/** \brief Parent class for local distance functions, where <b>local means pointwise</b>
 *  \see objective_functions.h
 */
class DistanceFunction_local
{
public:
    virtual double eval_scalar(const double f1, const double f2) = 0;
    virtual double eval_vector(const double * f1, const double * f2, const int D) = 0;
};

/** \brief Parent class for global distance functions, which <b>acts on the local distance</b>
 *  \see objective_functions.h
 */
class DistanceFunction_global
{
public:
    virtual double eval(const double * f1, const int P ) = 0;
};


/** TODO 
 * \see objective_functions.h
 */
class Dist_local_RelDiff : public DistanceFunction_local
{
    double eval_scalar(const double f1, const double f2);
    double eval_vector(const double * f1, const double * f2, const int D);
};

/** TODO 
 * \see objective_functions.h
 */
class Dist_local_AbsDiff : public DistanceFunction_local
{
    double eval_scalar(const double f1, const double f2);
    double eval_vector(const double * f1, const double * f2, const int D);
};


/** TODO 
 * \see objective_functions.h
 */
class Dist_global_Mean : public DistanceFunction_global
{
    double eval(const double * f, const int P);
};

/** TODO 
 * \see objective_functions.h
 */
class Dist_global_RMS : public DistanceFunction_global
{
    double eval(const double * f, const int P);
};

/** TODO 
 * \see objective_functions.h
 */
class Dist_global_Max : public DistanceFunction_global
{
    double eval(const double * f, const int P);
};

#endif // define __DISTANCE_FUNCTIONS_H__
