#ifndef __DISTANCE_FUNCTIONS_H__
#define __DISTANCE_FUNCTIONS_H__

/*
 *  ConcentricInterpolation
 *  Copyright (C) 2018  Felix Fritzen    ( felix.fritzen@mechbau.uni-stuttgart.de )
 *                      and Oliver Kunc  ( oliver.kunc@mechbau.uni-stuttgart.de )
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
 *  
 *  
 *  For details or if you like this software please refer to LITERATURE which
 *  contains also BIBTEX information.
 *  
 *  The latest version of this software can be obtained through https://github.com/EMMA-Group/ConcentricInterpolation
 *  
 *  
 */

#include <cmath>
#include <stdio.h>

// FF: documentation still missing

const double DistanceFunction_TOL = 1.e-5;

class DistanceFunction_local
{
public:
    virtual double eval_scalar(const double f1, const double f2) = 0;
    virtual double eval_vector(const double * f1, const double * f2, const int D) = 0;
};

class DistanceFunction_global
{
public:
    virtual double eval(const double * f1, const int P ) = 0;
};



class Dist_RelDiff_local : public DistanceFunction_local
{
    double eval_scalar(const double f1, const double f2);
    double eval_vector(const double * f1, const double * f2, const int D);
};

class Dist_AbsDiff_local : public DistanceFunction_local
{
    double eval_scalar(const double f1, const double f2);
    double eval_vector(const double * f1, const double * f2, const int D);
};



class Dist_Mean_global : public DistanceFunction_global
{
    double eval(const double * f, const int P);
};

class Dist_RMS_global : public DistanceFunction_global
{
    double eval(const double * f, const int P);
};

class Dist_Max_global : public DistanceFunction_global
{
    double eval(const double * f, const int P);
};

#endif // define __DISTANCE_FUNCTIONS_H__
