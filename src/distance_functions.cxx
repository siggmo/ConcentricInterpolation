#include "distance_functions.h"

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


