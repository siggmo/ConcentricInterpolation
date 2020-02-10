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

#include "concentric_interpolation.h"

ConcentricInterpolation::ConcentricInterpolation()
{
    RI_result = TI_result = nullptr;
    RI_result_length = TI_result_length = 0;
    initialized = false;
    printed_info_error  = false;
}

ConcentricInterpolation::ConcentricInterpolation(ConcentricData& a_SetupData)
{
    initialized = false;
    printed_info_error  = false;
    Setup(a_SetupData);
}

ConcentricInterpolation::~ConcentricInterpolation()
{
    delete [] RI_result;
    delete [] TI_result;

    RI_result = TI_result = nullptr;
    RI_result_length = TI_result_length = 0;
    initialized = false;
}

void ConcentricInterpolation::Setup(
        const ConcentricData& a_SetupData,
        const double a_gamma,
        const bool a_sym
         )
{
    // initial checks
    assert_msg( !initialized, "ERROR: ConcentricInterpolation can only be set up once atm\n");
    assert_msg( a_SetupData.IsInitialized(), "ERROR: data must be initialized before ConcentricInterpolation can be set up\n");
    double gamma = a_gamma;
    if( gamma<=0.01 )
    {
        gamma = 1.0;
        printf("# ATTENTION: setting gamma to default value %e\n", gamma);
    }

    // first set up TangentialInterpolant
    TangentialInterpolant.Setup(    a_SetupData.GetN_dir(),
                                    a_SetupData.GetD_inp(),
                                    a_SetupData.GetDirections(),
                                    a_sym,
                                    gamma
                               );

    // then set up RadialInterpolant, which requires TangentialInterpolant to be set up
    RadialInterpolant.Setup(    a_SetupData.GetValues(),
                                a_SetupData.GetN_dir(),
                                a_SetupData.GetN_rad(),
                                a_SetupData.GetD_val(),
                                a_SetupData.GetRadii(),
                                LinearInterpolant, // you may choose a different 1D interpolant here
                                TangentialInterpolant
                           );

    RI_result_length = size_t(a_SetupData.GetN_dir() * a_SetupData.GetD_val());
    RI_result = new double[ RI_result_length ];
    TI_result_length = size_t(a_SetupData.GetN_dir());
    TI_result = new double[ TI_result_length ];

    D_val       = a_SetupData.GetD_val();

    initialized = true;
}




void ConcentricInterpolation::Evaluate(
        GeneralData&  a_DataEvaluation
         )
{
    fprintf(stderr, "WARNING: ConcentricInterpolation::Evaluate() has not been tested yet!!\n");
    assert_msg( initialized, "ERROR: ConcentricInterpolation must be initialized before evaluation\n");
    assert_msg( (a_DataEvaluation.GetD_inp()==TangentialInterpolant.GetD_inp())
                    || (a_DataEvaluation.GetD_val()==D_val),
                "ERROR: dimension mismatch in ConcentricInterpolation::Error\n");

    // preparations
    const int N_pts = a_DataEvaluation.GetN_pts();
    const int D_inp = a_DataEvaluation.GetD_inp();
    assert_msg( D_inp==TangentialInterpolant.GetD_inp(), "ERROR: the evaluation data's dimension must match that of TI's input space\n");
    double CurrentRadius;
    double CurrentDirection[D_inp];

    // loop over evaluation points
    for(int n_pt=0; n_pt<N_pts; n_pt++) // n_pt = "number of point"
    {
        // get radius
        CurrentRadius = norm(a_DataEvaluation.GetCoords(n_pt), D_inp);
        // get direction
        scale(a_DataEvaluation.GetCoords(n_pt), CurrentRadius, CurrentDirection, D_inp);

        Evaluate( CurrentDirection, CurrentRadius, a_DataEvaluation.SetValues(n_pt) );
    }

}




double ConcentricInterpolation::Error( const ConcentricData& a_DataValidation, const int a_error_type, const int a_d_val_start, const int a_d_val_end, const int a_OnlyRadius )
{
    assert_msg( initialized, "ERROR: ConcentricInterpolation must be initialized before error computation\n");
    assert_msg( (a_DataValidation.GetD_inp()==TangentialInterpolant.GetD_inp())
                    || (a_DataValidation.GetD_val()==RadialInterpolant.GetD_val()),
                "ERROR: dimension mismatch in ConcentricInterpolation::Error\n");
    assert_msg( a_OnlyRadius<a_DataValidation.GetN_rad(),
                "ERROR: if ConcentricInterpolation::Error is called with too large special radius index\n");
    assert_msg( a_error_type==1 || a_error_type==2,
                "ERROR: only error types 1 and 2 are currently implemented in ConcentricInterpolation::Error\n");
    assert_msg( a_d_val_start>=-1 && a_d_val_end>=-1 && a_d_val_end<D_val && a_d_val_start<=a_d_val_end,
                "ERROR: impossible values of a_d_val_start and/or a_d_val_end in ConcentricInterpolation::Error\n");
    if( !printed_info_error )
    {
        if( a_error_type == 1)
            printf("# [ using mean relative error ]\n");
        else if( a_error_type == 2)
            printf("# [ using RMS absolute error ]\n");
        printed_info_error = true;
    }

    // preparations
    int       N_rad_eval    = a_DataValidation.GetN_rad();
    const int N_dir_eval    = a_DataValidation.GetN_dir();
    const int D_inp         = a_DataValidation.GetD_inp();
    double CurrentRadius;
    double * CurrentDirection   = new double[D_inp];
    double * CurrentValues      = new double[D_val];
    const int d_val_start   = ( a_d_val_start>=0 ? a_d_val_start : 0 );
    const int d_val_end     = ( a_d_val_end>=0 ? a_d_val_end : D_val-1 );
    double error_at_dir     = 0;
    double error_total      = 0;

    // loop validation radii
    for(int n_rad_eval=1; n_rad_eval<N_rad_eval; n_rad_eval++) // n_rad_eval = "number of evaluation radius". ATTENTION: starts at 1, meaning that the data at radius 0 (which has index 0) is NOT considered. that should be 0 anyways.
    {
        // if only one radius should be evaluated, n_rad_eval and N_rad_eval are now modified
        if(a_OnlyRadius>=0)
        {
            n_rad_eval = a_OnlyRadius;
            N_rad_eval = 1;  // number of evaluation radii is 1, which is important for error normalization
        }

        CurrentRadius = a_DataValidation.GetRadius(n_rad_eval);
        assert_msg( CurrentRadius>1.e-8, "WARNING: very small evaluation radius, might compromise relative error measures\n", false);

        // loop validation directions
        for(int n_dir_eval=0; n_dir_eval<N_dir_eval; n_dir_eval++) // n_dir_eval = "number of evaluation direction"
        {
            // fill current direction
            for(int d_inp=0; d_inp<D_inp; d_inp++) // d_inp = "dimension of input"
                CurrentDirection[d_inp] = a_DataValidation.GetDirection(n_dir_eval,d_inp);

            // evaluate at current point
            Evaluate( CurrentDirection, CurrentRadius, CurrentValues );

            // compute difference to validation data at current point
            for(int d_val=d_val_start; d_val<=d_val_end; d_val++) // ATTENTION: loop includes upper bound
                CurrentValues[d_val] -= a_DataValidation.GetValues()[n_dir_eval][n_rad_eval][d_val];

            /*****************************************/
            /* IMPLEMENT YOUR OWN ERROR MEASURE HERE */
            /* AND RIGHT BEFORE THE RETURN STATEMENT */
            /*****************************************/
            if( a_error_type == 1 ) // mean relative error
                error_at_dir = norm(CurrentValues+d_val_start, d_val_end-d_val_start+1)
                                / norm(a_DataValidation.GetValues()[n_dir_eval][n_rad_eval]+d_val_start, d_val_end-d_val_start+1);
            else if ( a_error_type == 2) // rooted mean square
            {
                error_at_dir = norm(CurrentValues+d_val_start, d_val_end-d_val_start+1);
                error_at_dir *= error_at_dir;
            }
            error_total += error_at_dir;
        }
    }

    delete [] CurrentDirection; CurrentDirection = nullptr;
    delete [] CurrentValues;    CurrentValues    = nullptr;

    // finalize total error
    if( a_error_type == 1) // mean relative error
        error_total /= double(N_rad_eval*N_dir_eval);
    else if( a_error_type == 2) // rooted mean square
        error_total = sqrt( error_total / double(N_rad_eval*N_dir_eval) );

    return error_total;
}




void ConcentricInterpolation::OptimizeGamma(
        const double    a_gamma_min,
        const double    a_gamma_max,
        const int       a_N_gamma_regular,
        const int       a_N_gamma_bisection,
        const ConcentricData&     a_DataValidation,
        const int       a_error_type,
        const int       a_d_val_start,
        const int       a_d_val_end,
        const int       a_OnlyRadius
       )
{
    // checks
    assert_msg( initialized, "ERROR: ConcentricInterpolation must be initialized before gamma optimization\n");
    assert_msg(0.<a_gamma_min, "ERROR in OptimizeGamma : a_gamma_min must be greater than zero\n");
    assert_msg(a_gamma_min<a_gamma_max, "ERROR in OptimizeGamma : a_gamma_min must be less than a_gamma_max\n");
    assert_msg(a_N_gamma_regular>=0, "ERROR in OptimizeGamma : num_regular must be non-negative\n");
    assert_msg(a_N_gamma_bisection>=0, "ERROR in OptimizeGamma : a_N_gamma_bisection must be non-negative\n");
    printed_info_error = false;

    // outputs
    printf("\r### Beginning optimization of gamma with the following parameters:\n");
    printf("#   a_gamma_min         = %lf\n", a_gamma_min);
    printf("#   a_gamma_max         = %lf\n", a_gamma_max);
    printf("#   num_gamma_regular   = %i\n", a_N_gamma_regular);
    printf("#   num_gamma_bisec     = %i\n", a_N_gamma_bisection);
    printf("#   error type          : ");
    if( a_error_type == 1 )
        printf("mean relative\n");
    else if( a_error_type == 2 )
        printf("RMS absolute\n");
    printf("#   components of vector of values for error computation: %i ... %i (inclusive, zero-based)\n", a_d_val_start, a_d_val_end);
    if( a_OnlyRadius<0 )
        printf("#   consider all radii of the provided validation data set\n");
    else
        printf("#   consider only radius with index %i of the provided validation data set\n", a_OnlyRadius);

    // error storage
    double ErrorMinimum;
    double Errors[a_N_gamma_regular+2];

    // initialize regular optimization
    double gammas[a_N_gamma_regular + 2];
    double gamma_delta = (a_gamma_max - a_gamma_min)/(a_N_gamma_regular + 1);
    for(int i_gamma=0; i_gamma<a_N_gamma_regular+2; i_gamma++)
        gammas[i_gamma] = a_gamma_min + i_gamma*gamma_delta;
    int index_best_gamma = 0;

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // perform regular optimization: sample the gammas grid from left to right and find its position with minimum error
    printf("# #gamma     gamma          error          new min?\n");
    int i_gamma = 0; // gamma index
    try{        // experience suggests that if there is no exception thrown for a_gamma_min, then there won't be any later on
        SetGamma(gammas[i_gamma], true, true);
        Errors[i_gamma] = Error(a_DataValidation, a_error_type, a_d_val_start, a_d_val_end, a_OnlyRadius);
        ErrorMinimum    = Errors[i_gamma];
        printf("\r  %4i     %f       %5.3le         !\n", i_gamma, gammas[i_gamma], Errors[i_gamma]);  // \r is carriage return for output of progress (called from InterpolateOnSet)
    } catch(int e) {
        fprintf(stderr, "\nEXCEPTION CAUGHT (%i) in OptimizeGamma : a_gamma_min doesn't work, probably too small.\n", e);
        exit(-1);
    }
    for(i_gamma=1; i_gamma<a_N_gamma_regular + 2; i_gamma++)
    {
        SetGamma(gammas[i_gamma], true, true);
        Errors[i_gamma] = Error(a_DataValidation, a_error_type, a_d_val_start, a_d_val_end, a_OnlyRadius);
        printf("\r  %4i     %f       %5.3le", i_gamma, gammas[i_gamma], Errors[i_gamma]);
        if(Errors[i_gamma]<ErrorMinimum)
        {
            printf("         !");
            ErrorMinimum = Errors[i_gamma];
            index_best_gamma = i_gamma;
        }
        printf("\n");
        fflush(stdout);
    }
    double gamma_best = gammas[index_best_gamma];
    printf("#   regular gamma optimization finished! best gamma is %lf with index %i giving error %5.3le\n", gammas[index_best_gamma], index_best_gamma, ErrorMinimum);
    // finish regular optimization by assuring that the best gamma is not on the boundary
    assert_msg(index_best_gamma!=0 && index_best_gamma!= a_N_gamma_regular+1, "ERROR in OptimizeGamma : gamma_best is on boundary. restart with different boundaries.\n");




    /////////////////////////////////////////////////////////////////////////////////////////////////////
    // bisection algorithm

    // initialize bisection algorithm: find the one neighbor of gamma_best within the regular grid
    // that has the lowest error. this is called gamma2
    double gamma2 = -1;
    double Error2;
    if( Errors[index_best_gamma-1] < Errors[index_best_gamma+1] )
    {
        gamma2 = gammas[index_best_gamma-1];
        Error2 = Errors[index_best_gamma-1];
        printf("#   best neighbor is next smaller gamma (%lf) with error %5.3le\n", gamma2, Error2);
    }
    else
    {
        gamma2 = gammas[index_best_gamma+1];
        Error2 = Errors[index_best_gamma+1];
        printf("#   best neighbor is next larger gamma (%lf) with error %5.3le\n", gamma2, Error2);
    }

    // do the bisection, i.e. evaluate within the interval and compare with the outer values
    if(a_N_gamma_bisection>0)
    {
        printf("#   beginning bisection\n");
        printf("# #gamma     gamma          error          new min?\n");
    }
    else
        return;

    double          gamma_mid           = -1;       // the value of the gamma bisecting the current interval
    double          Error_mid           = -1;       // the value of the error at the bisecting gamma
    bool            keep_best_gamma     = false;
    const double    gamma_change_tol    = 1.e-8;    // tolerance for gamma change
    const double    Errors_change_tol   = 1.e-5;    // tolerance for error change
    const double    bisection_factor    = 1.;       // the bisection factor (or step size), may be adjusted within the interval (0,2)
    for(int n_bisec=0; n_bisec<a_N_gamma_bisection; n_bisec++)
    {
        // define bisecting gamma
        gamma_mid = gamma_best + (gamma2-gamma_best)*0.5*bisection_factor;

        // check gamma change
        if( fabs(gamma_mid-gamma_best)<gamma_change_tol )
        {
            printf("\n#   quitting bisection because gamma change is below tolerance (%5.3le).\n", gamma_change_tol);
            SetGamma(gamma_best);
            break;
        }

        // set the bisecting gamma and compute error
        SetGamma(gamma_mid, true, true);
        Error_mid = Error(a_DataValidation, a_error_type, a_d_val_start, a_d_val_end, a_OnlyRadius);
        printf("\r  %4i     %f       %8.6le", n_bisec, gamma_mid, Error_mid); fflush(stdout);

        // check for new minimum error
        if(Error_mid<ErrorMinimum)
        {
            printf("      !");
            Error2 = ErrorMinimum;
            ErrorMinimum = Error_mid;
            gamma2 = gamma_best;
            gamma_best = gamma_mid;
        }
        else
        {
            Error2 = Error_mid;
            gamma2 = gamma_mid;
            keep_best_gamma = true;
        }
        printf("\n");
        if( (Error2-ErrorMinimum)/ErrorMinimum < Errors_change_tol )
        {
            printf("#   quitting bisection due to low relative Errors change (tolerance is %5.3le)\n", Errors_change_tol);
            // before quitting, set gamma to be the current best gamma
            if(keep_best_gamma)
                SetGamma(gamma_best);
            break;
        }
        keep_best_gamma = false;
    }

    printf("#   finished gamma optimization with gamma = %8.6le yielding errors %8.6le\n", gamma_best, ErrorMinimum);
}




void ConcentricInterpolation::SetGamma( const double a_gamma, const bool a_recompute_kernel_matrix, const bool a_quiet )
{
    TangentialInterpolant.SetGamma( a_gamma, a_recompute_kernel_matrix, a_quiet );
    RadialInterpolant.MultiplyKernelMatrixAndReorder( TangentialInterpolant );
}


void ConcentricInterpolation::PrintInfo(FILE * fileptr)
{
    fprintf(fileptr, "###############################################\n");
    fprintf(fileptr, "# no. support directions:  %5i                #\n", TangentialInterpolant.GetN_dir_supp());
    fprintf(fileptr, "# no. support radii:       %5i                #\n", RadialInterpolant.GetN_rad_supp()); // TODO
    fprintf(fileptr, "# kernel parameter:          %17.9e  #\n", TangentialInterpolant.GetGamma());
    fprintf(fileptr, "# no. input dimensions:    %5i                #\n", TangentialInterpolant.GetD_inp());
    fprintf(fileptr, "# no. output components:   %5i                #\n", RadialInterpolant.GetD_val());
    fprintf(fileptr, "###############################################\n");
}
