#include "data.h"
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

Data::Data():
D(-1),
P(-1),
coordinates(0),
values(0),
gradients(0),
hessians(0)
{}

Data::Data(const int a_D, const int a_P):
D(a_D),
P(a_P)
{
    initialize();
}

Data::Data(Data& dat, const bool function_data):
D(dat.D),
P(dat.P)
{
    initialize();

    // always copy point coordinates
    for(int p=0; p<P; p++)
        for(int d=0; d<D; d++)
            coordinates[p][d] = dat.SafeAccess_coord(p,d);

    // conditionally copy function data
    if(function_data)
        for(int p=0; p<P; p++)
        {
            values[p] = dat.values[p];
            for(int d=0; d<D; d++)
                gradients[p][d] = dat.SafeAccess_gradients(p,d);
            for(int component=0; component<D*D; component++)
                hessians[p][component] = dat.SafeAccess_hessians(p,component);
        }
}

Data::~Data()
{
    if(D<=0 && P<=0)
    {}
    else
    {
        free_matrix(coordinates, P);
        delete [] values; values = 0;
        free_matrix(gradients, P);
        free_matrix(hessians, P);
    }
}

void Data::initialize()
{
    assert_msg(D>0 && P>0, "ERROR in Data constructor: D and P must be greater than zero\n");
    coordinates = alloc_matrix(P,D);
    values =      new double[P];
    gradients =   alloc_matrix(P,D);
    hessians =    alloc_matrix(P,D*D);
}


// /////////////////////////////////////////////////////////////////////////////////////////////

DataTraining::DataTraining():
D(-1),
N(-1),
R(-1),
directions(0),
values(0),
radialderiv(0)
{}

DataTraining::DataTraining(const int a_D, const int a_N, const int a_R):
D(a_D),
N(a_N),
R(a_R)
{
    initialize();
}

DataTraining::DataTraining(DataTraining& dat, const bool function_data):
D(dat.D),
N(dat.N),
R(dat.R)
{
    initialize();

    // copy point directions
    for(int n=0; n<N; n++)
        for(int d=0; d<D; d++)
            directions[n][d] = dat.SafeAccess_dir(n,d);

    // copy radii
    for(int r=0; r<R; r++)
        radii[r] = dat.SafeAccess_radii(r);

    // conditionally copy function data
    if(function_data)
        for(int n=0; n<N; n++)
            for(int r=0; r<R; r++)
            {
                values[n][r] = dat.SafeAccess_values(n,r);
                radialderiv[n][r] = dat.SafeAccess_radderiv(n,r);
            }
}

DataTraining::~DataTraining()
{
    if(N>0)
    {
        free_matrix(directions, N);
        free_matrix(values, N);
        free_matrix(radialderiv, N);
    }
    if(R>0)
        delete [] radii, radii=0;
}

void DataTraining::initialize()
{
    assert_msg(D>0 && N>0 && R>0, "ERROR in DataTraining constructor: D, N, and R must be greater than zero\n");
    directions  = alloc_matrix(N,D);
    values      = alloc_matrix(N,R);
    radialderiv = alloc_matrix(N,R);
    radii       = new double[R];
}
