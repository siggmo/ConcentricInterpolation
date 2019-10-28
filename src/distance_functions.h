#ifndef __DISTANCE_FUNCTIONS_H__
#define __DISTANCE_FUNCTIONS_H__

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
