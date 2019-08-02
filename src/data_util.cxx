#include "data_util.h"

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
void InterpolateOnSet(  ConcentricInterpolation * interpolation,
                        Data * data_interpolation,
                        const bool evaluate_gradients,
                        const bool evaluate_hessians
                     )
{
    assert_msg( interpolation->GetGamma() > 0 , "WARNING in InterpolateOnSet: gamma must be > 0, except for zero radius.\n", false);
    assert_msg( !( evaluate_hessians && !evaluate_gradients ) , "ERROR in InterpolateOnSet: if evaluate_hessians == true, then evaluate_gradients must also be true\n");
    progress::init(); /* optional progress visualization on shell */
    const int P = data_interpolation->P;
    for(int p=0; p<P; p++)
    {
        progress::display(p, P-1); /* optional progress display on shell */
        // interpolate value
        data_interpolation->SafeAccess_values(p) = interpolation->Interpolate( data_interpolation->SafeAccess_coord(p) );
        // get gradient, optional
        if( evaluate_gradients )
            interpolation->Gradient( data_interpolation->SafeAccess_gradients(p) );
        // get hessian, optional
        if( evaluate_hessians )
            interpolation->Hessian( data_interpolation->SafeAccess_hessians(p) );
    }
}
/* *************************************************************************************** */
// default constructor; (should not be used unless you know what you do)
TestOnSet::TestOnSet():
D(-1),
P(-1)
{
    interpolation = 0;
    data_reference = 0;
    dist_local_values = 0;
    dist_local_gradients = 0;
    dist_local_hessians = 0;
    distfunc_local = 0;
    distfunc_global = 0;
    objective_function = 0;
}

TestOnSet::TestOnSet(   ConcentricInterpolation * a_interpolation,
                        Data * a_data_reference,
                        DistanceFunction_local * a_distfunc_local,
                        DistanceFunction_global* a_distfunc_global,
                        double (*a_objective_function)(const double, const double, const double)
                    ):
interpolation(a_interpolation),
data_interpolation(*a_data_reference),  // default copy constructor copies only coordinates
data_reference(a_data_reference),
D(data_interpolation.D),
P(data_interpolation.P),
distfunc_local(a_distfunc_local),
distfunc_global(a_distfunc_global),
objective_function(a_objective_function)
{
    // allocate memory for 
    dist_local_values = new double[P];
    dist_local_gradients = new double[P];
    dist_local_hessians = new double[P];
}

TestOnSet::~TestOnSet()
{
    delete[] dist_local_values; dist_local_values=0;
    delete[] dist_local_gradients; dist_local_gradients=0;
    delete[] dist_local_hessians; dist_local_hessians=0;
}

double TestOnSet::run(  const bool evaluate_gradients,
                        const bool evaluate_hessians
              )
{
    // consistency check !( evaluate_hessians && !evaluate_gradients ) is done in InterpolateOnSet
    
    // run the interpolation on data_interpolation.coordinates
    InterpolateOnSet( interpolation, &data_interpolation, evaluate_gradients, evaluate_hessians);
    
    // compute the point-wise distances between interpolated data and reference data
    for(int p=0; p<P; p++)
    {
        dist_local_values[p] = distfunc_local->eval_scalar( data_interpolation.SafeAccess_values(p), data_reference->SafeAccess_values(p) );
        if(evaluate_gradients)
            dist_local_gradients[p] = distfunc_local->eval_vector( data_interpolation.SafeAccess_gradients(p), data_reference->SafeAccess_gradients(p), D );
        if(evaluate_hessians)
            dist_local_hessians[p] = distfunc_local->eval_vector( data_interpolation.SafeAccess_hessians(p), data_reference->SafeAccess_hessians(p), D*D );
    }
    
    // compute the global distances
    double dist_global_values    = distfunc_global->eval(dist_local_values, P);   /* distance of interpolated values to reference values */
    double disg_global_gradients = 0;
    double dist_global_hessians = 0;
    if(evaluate_gradients)
        disg_global_gradients = distfunc_global->eval(dist_local_gradients, P);/* distance of interpolated gradients to reference gradients */
    if(evaluate_hessians)
        dist_global_hessians  = distfunc_global->eval(dist_local_hessians, P); /* distance of interpolated Hessians to reference Hessians */
    
    // return the value of the objective function
    return objective_function(dist_global_values, disg_global_gradients, dist_global_hessians);
}
/* *************************************************************************************** */
double    OptimizeGamma(ConcentricInterpolation * interpolation,
                        Data * data_reference,
                        DistanceFunction_local * distfunc_local,
                        DistanceFunction_global* distfunc_global,
                        double (*objective_function)(const double, const double, const double),
                        const double gamma_min,
                        const double gamma_max,
                        const int num_regular,
                        const int num_bisec,
                        const double bisection_factor,
                        const bool do_gradients,
                        const bool do_hessians
                       )
{
    assert_msg(0.<gamma_min, "ERROR in OptimizeGamma : gamma_min must be greater than zero\n");
    assert_msg(gamma_min<gamma_max, "ERROR in OptimizeGamma : gamma_min must be less than gamma_max\n");
    assert_msg(num_regular>=0, "ERROR in OptimizeGamma : num_regular must be non-negative\n");
    assert_msg(num_bisec>=0, "ERROR in OptimizeGamma : num_bisec must be non-negative\n");
    assert_msg(bisection_factor>0 && bisection_factor<2, "ERROR in OptimizeGamma : it must hold    0 < bisection_factor < 2  \n");

    
    printf("\r$$$ Beginning optimization of gamma with the following parameters:\n");
    printf("$   gamma_min        = %lf\n", gamma_min);
    printf("$   gamma_max        = %lf\n", gamma_max);
    printf("$   num_regular      = %i\n", num_regular);
    printf("$   num_bisec        = %i\n", num_bisec);
    printf("$   bisection_factor = %lf\n", bisection_factor);
    printf("$   do_gradients     = %i\n", int(do_gradients));
    printf("$   do_hessians      = %i\n", int(do_hessians));
    
    // initialize the TestOnSet object s.t. there is permanent memory for interpolation results
    TestOnSet test_object(interpolation, data_reference, distfunc_local, distfunc_global, objective_function);
    
    // initialize regular optimization
    double * gammas = new double[num_regular + 2];
    double dgamma = (gamma_max - gamma_min)/(num_regular + 1);
    for(int igamma=0; igamma<num_regular+2; igamma++)
        gammas[igamma] = gamma_min + igamma*dgamma;
    int ind_gamma_best = 0;
    double * obj_fcn_value = new double[num_regular + 2];  // obj_fcn_value for each gamma
    double obj_fcn_value_min = 9e99;                        // minimum of all the previous obj_fcn_value
    
    // do regular optimization: sample the gammas grid from left to right and find its position with minimum obj_fcn_value
    printf("$   regular gamma optimization on the interval [%lf, %lf] with %i intermediate points:\n", gammas[0], gammas[num_regular + 1], num_regular);
    printf("$ #gamma     gamma      obj_fcn_value      new min?\n");
    int igamma = 0;                                 // gamma counter
    try{        // experience suggests that if there is no exception thrown for gamma_min, then there won't be any later on, i.e. exception handling is only necessary at gamma_min
        interpolation->SetGamma(gammas[igamma]);
        obj_fcn_value[igamma] = test_object.run( do_gradients, do_hessians );
        obj_fcn_value_min = obj_fcn_value[igamma];
        printf("\r  %4i     %f       %5.3le         !\n", igamma, gammas[igamma], obj_fcn_value[igamma]);  // \r is carriage return for output of progress (called from InterpolateOnSet)
    } catch(int e) {
        fprintf(stderr, "\nEXCEPTION CAUGHT (%i) in OptimizeGamma : gamma_min doesn't work, probably too small.\n", e);
        exit(-1);
    }
    for(igamma=1; igamma<num_regular+2; igamma++)
    {
        interpolation->SetGamma(gammas[igamma]);
        obj_fcn_value[igamma] = test_object.run( do_gradients, do_hessians );
        printf("\r  %4i     %f       %5.3le", igamma, gammas[igamma], obj_fcn_value[igamma]);
        if(obj_fcn_value[igamma]<obj_fcn_value_min)
        {
            printf("         !");
            obj_fcn_value_min = obj_fcn_value[igamma];
            ind_gamma_best = igamma;
        }
        printf("\n");
        fflush(stdout);
    }
    double gamma_best = gammas[ind_gamma_best];
    printf("$   regular gamma optimization finished! best gamma is %lf with index %i giving obj_fcn_value %5.3le\n", gammas[ind_gamma_best], ind_gamma_best, obj_fcn_value_min);
    // finish regular optimization by assuring that the best gamma is not on the boundary
    assert_msg(ind_gamma_best!=0 && ind_gamma_best!= num_regular+1, "ERROR in OptimizeGamma : gamma_best is on boundary. restart with different boundaries.\n");

    
    
    // initialize bisectional algorithm by finding the best value of the neighboring gammas, i.e. gamma_2
    double gamma_2 = -1;
    double obj_fcn_value_2 = 9e99;
    if( obj_fcn_value[ind_gamma_best-1] < obj_fcn_value[ind_gamma_best+1] )
    {
        gamma_2 = gammas[ind_gamma_best-1];
        obj_fcn_value_2 = obj_fcn_value[ind_gamma_best-1];
        printf("$      best neighbor is next smaller gamma (%lf) with obj_fcn_value %5.3le\n", gamma_2, obj_fcn_value_2);
    }
    else
    {
        gamma_2 = gammas[ind_gamma_best+1];
        obj_fcn_value_2 = obj_fcn_value[ind_gamma_best+1];
        printf("$      best neighbor is next larger gamma (%lf) with obj_fcn_value %5.3le\n", gamma_2, obj_fcn_value_2);
    }
    
    // do the bisection, i.e. evaluate within the interval and compare with the outer values
    if(num_bisec>0)
        printf("\n$   beginning bisection\n$ #gamma     gamma      obj_fcn_value      new min?\n");
    double          gamma_mid               = -1;
    double          obj_fcn_value_mid       = 9e99;
    bool            keep_best_gamma         = false;
    const double    gamma_change_tol        = 1.e-8;
    const double    obj_fcn_value_change_tol= 1.e-5;
    for(int it_bisec=0; it_bisec<num_bisec; it_bisec++)
    {
        gamma_mid = gamma_best + (gamma_2-gamma_best)*0.5*bisection_factor;
        if( fabs(gamma_mid-gamma_best)<gamma_change_tol )
        {
            printf("\n$   quitting bisection because gamma change is below tolerance (%5.3le).\n", gamma_change_tol);
            interpolation->SetGamma(gamma_best);
            break;
        }
        interpolation->SetGamma(gamma_mid);
        obj_fcn_value_mid = test_object.run( do_gradients, do_hessians );
        printf("\r  %4i     %f       %5.3le", it_bisec, gamma_mid, obj_fcn_value_mid); fflush(stdout);
        if(obj_fcn_value_mid<obj_fcn_value_min)
        {
            printf("         !");
            obj_fcn_value_2 = obj_fcn_value_min;
            obj_fcn_value_min = obj_fcn_value_mid;
            gamma_2 = gamma_best;
            gamma_best = gamma_mid;
        }
        else// i.e. if obj_fcn_value_mid>=obj_fcn_value_min
        {
            obj_fcn_value_2 = obj_fcn_value_mid;
            gamma_2 = gamma_mid;
            keep_best_gamma = true;
        }
        printf("\n");
        if( (obj_fcn_value_2-obj_fcn_value_min)/obj_fcn_value_min < obj_fcn_value_change_tol )
        {
            printf("$   quitting bisection due to low relative obj_fcn_value change (tolerance is %5.3le)\n", obj_fcn_value_change_tol);
            // before quitting, set gamma to be the current best gamma
            if(keep_best_gamma)
                interpolation->SetGamma(gamma_best);
            break;
        }
        keep_best_gamma = false;
    }
    
    printf("$   finished gamma optimization with gamma = %5.3le yielding obj_fcn_value %5.3le\n", gamma_best, obj_fcn_value_min);
    
    // terminate
    delete [] obj_fcn_value; obj_fcn_value = 0;
    delete [] gammas; gammas = 0;
    
    return gamma_best;
}
/* *************************************************************************************** */
void BuildInterpolationFromTrainingData(
    ConcentricInterpolation & interpolation,
    DataTraining * data_training
     )
{
    const int D = data_training->D;
    const int N = data_training->N;
    const int R = data_training->R;

    // allocate the memory in the interpolation
    interpolation.Allocate(N,D);
    
    // loop over directions: add data along each direction, i.e. build the cubic interpolants
    for(int n=0; n<N; n++)
        interpolation.AddInterpolationDataDF(
            data_training->Access_radii(),
            data_training->SafeAccess_values(n),
            data_training->SafeAccess_radderiv(n),
            R,
            data_training->SafeAccess_dir(n)
        );
}
/* *************************************************************************************** */
void WriteDataGnuplot( Data & data, char * output_filename, const bool write_coordinates, const bool write_values, const bool write_gradients, const bool write_hessians )
{
    assert_msg( data.D>0 && data.P>0,
        "ERROR in WriteDataGnuplot: data must have been initialized\n");
    const int D = data.D;
    const int P = data.P;
    
    // store results into the output file. see the gnuplot script for visualization suggestions.
    printf("$$$ Storing data to  %s\n", output_filename);
    FILE * OUTPUT_FILE = fopen(output_filename, "w");
    assert_msg( OUTPUT_FILE!=0, "ERROR in WriteDataGnuplot: output file could not be opened\n");
    
    // print header (as a comment, i.e. preceeded by #)
    fprintf(OUTPUT_FILE, "#");
    if(write_coordinates)
        for(int d=0; d<D; d++)
                fprintf(OUTPUT_FILE,    "          coord %2i           ", d);
    if(write_values)
        fprintf(OUTPUT_FILE,            "            value            ");
    if(write_gradients)
        for(int d=0; d<D; d++)
            fprintf(OUTPUT_FILE,        "     gradient component %2i   ", d);
    if(write_hessians)
        for(int dd=0; dd<D*D; dd++)
            fprintf(OUTPUT_FILE,        "     Hessian component %3i   ", dd);
    fprintf(OUTPUT_FILE, "\n");
    
    // print data
    for(int p=0; p<data.P; p++)
    {
        if(write_coordinates)
            for(int d=0; d<D; d++)
                    fprintf(OUTPUT_FILE, " %28.21e", data.SafeAccess_coord(p,d));
        if(write_values)
            fprintf(OUTPUT_FILE, " %28.21e", data.SafeAccess_values(p));
        if(write_gradients)
            for(int d=0; d<D; d++)
                fprintf(OUTPUT_FILE, " %28.21e", data.SafeAccess_gradients(p,d));
        if(write_hessians)
            for(int dd=0; dd<D*D; dd++)
                fprintf(OUTPUT_FILE, " %28.21e", data.SafeAccess_hessians(p,dd));
        fprintf(OUTPUT_FILE, "\n");
    }
    fclose(OUTPUT_FILE), OUTPUT_FILE = 0;
}
