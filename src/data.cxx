#include "data.h"
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
