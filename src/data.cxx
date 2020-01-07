/*
 *  COPYRIGHT NOTES
 *
 *  ConcentricInterpolation
 * Copyright (C) 2018-2019 by Felix Fritzen (fritzen@mechbau.uni-stuttgart.de)
 *                         and Oliver Kunc (kunc@mechbau.uni-stuttgart.de
 * All rights reserved.
 *
 * This source code is licensed under the BSD 3-Clause License found in the
 * LICENSE file in the root directory of this source tree.
 *
 *  This software package is related to the research article
 *
 *  Authors: Oliver Kunc and Felix Fritzen
 *  Title  : Generation of energy-minimizing point sets on spheres and their
 *           application in mesh-free interpolation and differentiation
 *  Journal: Advances in Computational Mathematics 45(5-6), pp. 3021-3056
 *  Year   : 2019
 *  URL    : https://doi.org/10.1007/s10444-019-09726-5
 *
 *  The latest version of this software can be obtained through
 *  https://github.com/EMMA-Group/ConcentricInterpolation
 *
 */

#include "data.h"

GeneralData::GeneralData()
{
    initialized     = false;
    Coords = Values = nullptr;
    D_inp = D_val = N_pts= 0;
}

void GeneralData::Resize(
        const int a_D_inp,
        const int a_D_val,
        const int a_N_pts
         )
{
    Free();

    D_inp  = a_D_inp;
    D_val  = a_D_val;
    N_pts  = a_N_pts;
    Coords = alloc_matrix(N_pts,D_inp);
    Values = alloc_matrix(N_pts,D_val);
    initialized = true;
}

void GeneralData::Free()
{
    if(!initialized)
        return;

    free_matrix(Coords, N_pts);
    free_matrix(Values, N_pts);
    D_inp = D_val = N_pts = 0;
    initialized = false;
}









ConcentricData::ConcentricData()
{
    Directions  = nullptr;
    Radii       = nullptr;
    Values      = nullptr;
    D_inp   = D_val = N_dir = N_rad = 0;
    initialized = false;
}

ConcentricData::ConcentricData(
        char * a_FilenameSupportDirections,
        char * a_FilenameSupportRadii,
        char * a_FilenameSupportValues,
        const int a_D_val
                )
{
    // allocate & read directions
    Directions = ReadMatrix( &N_dir, &D_inp, a_FilenameSupportDirections );

    // allocate & read radii
    int must_be_one = 0;
    double ** temp = ReadMatrix( &must_be_one, &N_rad, a_FilenameSupportRadii );
        assert_msg( must_be_one==1, "ERROR: radius file must contain radii in one row\n");
        assert_msg( N_rad>=2, "ERROR: at least two radii must be provided\n");
    Radii = new double[N_rad];
    for(int r=0; r<N_rad; r++)
        Radii[r] = temp[0][r];
    free_matrix(temp, must_be_one);

    // allocate & read values
    assert_msg( a_D_val>=1, "ERROR: D_val must be at least one\n");
    D_val = a_D_val;
    int must_be_N_dir = 0;
    int must_be_D_val_times_N_rad = 0;
    double ** temp2 = ReadMatrix( &must_be_N_dir, &must_be_D_val_times_N_rad, a_FilenameSupportValues );
    printf("# ConcentridData: reading %s\n", a_FilenameSupportValues);
        assert_msg( must_be_N_dir==N_dir, "ERROR: values file must have same number of rows as directions file\n");
        assert_msg( must_be_D_val_times_N_rad==(D_val*N_rad), "ERROR: values file's columns are not consistent with number of radii and dimension of values\n");
    Values = alloc_array3(N_dir, N_rad, D_val);
    for(int n_dir=0; n_dir<N_dir; n_dir++)
        for(int n_rad=0; n_rad<N_rad; n_rad++)
            for(int d=0; d<D_val; d++)
                Values[n_dir][n_rad][d] = temp2[n_dir][n_rad*D_val+d];
    free_matrix(temp2, must_be_N_dir);

    // ready to go
    initialized = true;
}

void ConcentricData::Free()
{
    if(!initialized)
        return;

    free_matrix(Directions, N_dir);
    delete [] Radii; Radii = 0;
    free_array3(Values, N_dir, N_rad);

    D_inp = D_val = N_dir = N_rad = 0;
    initialized = false;
}

void ConcentricData::Resize(
        const int a_D_inp,
        const int a_D_val,
        const int a_N_dir,
        const int a_N_rad
         )
{
    Free();

    D_inp       = a_D_inp;
    D_val       = a_D_val;
    N_dir       = a_N_dir;
    N_rad       = a_N_rad;
    Directions  = alloc_matrix(N_dir,D_inp);
    Radii       = new double[N_rad];
    Values      = alloc_array3(N_dir,N_rad,D_val);
    initialized = true;
}
