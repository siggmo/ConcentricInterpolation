#include "create_data.h"
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

Data * CreateData(  double ** coordinates,
                    const int P,
                    const int D,
                    AnalyticalFunction & func
                 )
{
    // create Data object permanently
    Data * newData = new Data(D,P);
    
    // fill coordinates and function values
    for(int p=0; p<P; p++)
    {
        for(int d=0; d<D; d++)
            newData->SafeAccess_coord(p,d) = coordinates[p][d];
        
        func.evaluate(  coordinates[p],
                        newData->SafeAccess_values(p),
                        newData->SafeAccess_gradients(p),
                        newData->SafeAccess_hessians(p)
                     );
    }
    
    // return pointer to Data object
    return newData;
}
/* *************************************************************************************** */
Data * CreateData(  double ** directions,
                    const int N,
                    const int D,
                    double * radii,
                    const int R,
                    AnalyticalFunction & func
                 )
{
    const int P=N*R; // number of points
    
    // create the coordinates as tensor product of directions and radii
    double ** coordinates = alloc_matrix(P,D);
    for(int n=0; n<N; n++)
        for(int r=0; r<R; r++)
            for(int d=0; d<D; d++)
                coordinates[n*R+r][d] = directions[n][d] * radii[r];
    
    // call main CreateData function
    Data * data_pointer = CreateData(coordinates, P, D, func);
    
    // free allocated memory
    free_matrix(coordinates,P);
    
    // return result
    return data_pointer;
}
/* *************************************************************************************** */
Data * CreateData(  char * DIRECTIONS_FN,
                    double * radii,
                    const int R,
                    AnalyticalFunction & func
                 )
{
    printf("$$$ Begin creating data. Coordinates via tensor product of radii and file-read directions.\n$   Function data via analytical function on coordinates.\n");
    
    // read directions
    printf("$   reading target directions from file %s\n", DIRECTIONS_FN);
    
    int N=0;                // number of directions that have been read so far
    int D=0;                // dimension of the read directions
    double ** directions    // N-by-D matrix containing the coordinates of the direction vectors
                          = ReadMatrix( &N, &D, DIRECTIONS_FN);    // read the directions from the file, store N and D and return pointer to matrix

    assert_msg(D==func.GetDim() , "ERROR in CreateData: dimension of directions must be consistent with that of analytical function\n");
    
    // call directions-overloaded function
    Data * newData = CreateData(directions, N, D, radii, R, func);
    
    // free memory for directions
    free_matrix(directions, N);
    
    return newData;
}
/* *************************************************************************************** */
DataTraining * CreateDataTraining(
                    char * DIRECTIONS_FN,
                    double * radii,
                    const int R,
                    AnalyticalFunction & func
                 )
{
    printf("$$$ Begin creating training data. Store directions, radii,\n$   and values and radial derivatives of analytical function.\n");
    
    // read directions
    printf("$   reading training directions from file %s\n", DIRECTIONS_FN);
    
    int N=0;                // number of directions that have been read so far
    int D=0;                // dimension of the read directions
    double ** directions    // N-by-D matrix containing the coordinates of the direction vectors
                          = ReadMatrix( &N, &D, DIRECTIONS_FN);    // read the directions from the file, store N and D and return pointer to matrix
    
    assert_msg(D==func.GetDim() , "ERROR in CreateTrainingData: dimension of directions must be consistent with that of analytical function\n");
    
    // create DataTraining object
    DataTraining * newDataTraining = new DataTraining(D,N,R);
    
    double temp_gradient[D];
    double temp_coordinates[D];
    for(int n=0; n<N; n++)
    {
        // fill directions
        for(int d=0; d<D; d++)
            newDataTraining->SafeAccess_dir(n,d) = directions[n][d];
        
        // fill values and radial derivatives
        for(int r=0; r<R; r++)
        {
            for(int d=0; d<D; d++)
                temp_coordinates[d] = directions[n][d] * radii[r];
            func.evaluate(temp_coordinates, newDataTraining->SafeAccess_values(n,r), temp_gradient);
            newDataTraining->SafeAccess_radderiv(n,r) = VecVecMul(temp_gradient, directions[n], D);
        }
    }
    
    // fill radii
    for(int r=0; r<R; r++)
        newDataTraining->SafeAccess_radii(r) = radii[r];
    
    return newDataTraining;
}
    
/* *************************************************************************************** */
