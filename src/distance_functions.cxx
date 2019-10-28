#include "distance_functions.h"

/*
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
 *  The latest version of this software can be obtained through https://github.com/EMMA-Group/ConcentricInterpolation
 *
 *
 */
/* *************************************************************************************** */
double Dist_local_RelDiff::eval_scalar(const double f1, const double f2)
{
    const double f2abs = fabs(f2);
    if(f2abs>DistanceFunction_TOL)
        return fabs(f1-f2)/f2abs;
    else
        return fabs(f1-f2)/DistanceFunction_TOL;
}
double Dist_local_RelDiff::eval_vector(const double * vec1, const double * vec2, const int D)
{
    double value = 0;
    double normvec2 = 0;
    for(int d=0; d<D; d++)
    {
        value += (vec1[d]-vec2[d])*(vec1[d]-vec2[d]);
        normvec2 += vec2[d]*vec2[d];
    }
    if( sqrt(normvec2)>DistanceFunction_TOL )
        return sqrt(value/normvec2);
    else
        return sqrt(value/DistanceFunction_TOL);
}
/* *************************************************************************************** */
double Dist_local_AbsDiff::eval_scalar(const double f1, const double f2)
{
    return fabs(f1-f2);
}

double Dist_local_AbsDiff::eval_vector(const double * vec1, const double * vec2, const int D)
{
    double value = 0;
    for(int d=0; d<D; d++)
        value += (vec1[d]-vec2[d])*(vec1[d]-vec2[d]);
    return sqrt(value);
}
/* *************************************************************************************** */
double Dist_global_Mean::eval(const double * f, const int P)
{
    double value = 0;
    for(int p=0; p<P; p++)
        value += f[p];
    return value/P;
}
double Dist_global_RMS::eval(const double * f, const int P)
{
    double value = 0;
    for(int p=0; p<P; p++)
        value += f[p]*f[p];
    return sqrt(value/P);
}

double Dist_global_Max::eval(const double * f, const int P)
{
    double fmax = fabs(f[0]);
    for(int p=1; p<P; p++)
        if(fabs(f[p])>fmax)
            fmax=fabs(f[p]);
    return fmax;
}


