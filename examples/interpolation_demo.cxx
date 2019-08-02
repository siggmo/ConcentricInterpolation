#include "concentric_interpolation.h"   // the main program
#include "data_util.h"                  // utility functions concerning data. includes the kernel parameter optimization
#include "objective_functions.h"        // objective functions with respecto to which the kernel parameter may be optimized
#include "quadratic_interpolant.h"      // piecewise quadratic polynomials for radial interpolation. see also "cubic_interpolant.h"


int main(int argc, char* argv[])
{
    // all arguments are passed via the command line.
    // a gamma argument must be given in any case. if gamma search is performed (see below), then the gamma argument from the command line will be ignored.
    assert_msg(argc==7, "USAGE: [this] [filename training directions] [filename training data] [filename validation directions] [filename validtion data] [filename of output file] [gamma]\n");
    
    // symmetry flag
    const bool SYM = false;         // use symmetrized interpolation scheme? only set this true if the taining directions were created with a symmetrized kernel.

    // temporary variables
    int read_rows = -1, read_cols = -1;
    
    /* ******************************************************************************
    * initialize and setup the interpolation
    */
    ConcentricInterpolation interpolation( SYM );
    
    // read directions from text file
    double ** matrix_directions_training = UTILITY::ReadMatrix( &read_rows, &read_cols, argv[1]);
    printf("$$$ Read training directions: %i x %i\n", read_rows, read_cols);
    const int D = read_cols;                // number of space dimensions (in the notation of the paper: d=D-1)
    const int N_training = read_rows;       // number of training directions
    
    // read data from text file. contains N_training x R_training values, i.e. each row is the training values along one training direction
    double ** matrix_data_training = UTILITY::ReadMatrix(&read_rows, &read_cols, argv[2]);
    assert_msg(read_rows==N_training, "Error: number of training directions does not match number of training data directions\n");
    const int R_training = read_cols;       // number of training radii
    const double training_r[] = {0.,    0.2500,    0.5000,    0.7500,    1.0000};   // training radii
    printf("$   Trianing radii: "), print_matrix(training_r,1,R_training);
    interpolation.Allocate(N_training,D);   // allocate memory in the interpolation scheme
    
    // set the supporting points of the interpolation scheme
    printf("$   Set up interpolation along %i directions at %i radii.\n", N_training, R_training);
    for(int i_direction=0; i_direction<N_training; i_direction++)
    {
        // set training data
        interpolation.Sq[i_direction].SetData(R_training, training_r, matrix_data_training[i_direction]);   // the supporting points of the interpolation polynomials are defined

        // set training direction
        const double l = norm( matrix_directions_training[i_direction], D );
        for(int d=0; d<D; d++) interpolation.m_X[interpolation.N*D+d] = matrix_directions_training[i_direction][d]/ l; // make sure the Euclidean norm is 1 for the vectors representing directions
        
        // increment the counter for the dimension of m_X, i.e. the number of training directions
        interpolation.N++;
    }

    
    /* ******************************************************************************
    * specify error measure and objective function
    */

    printf("$$$ Specify error measure\n");
    const bool do_gradients = false;        // consider gradients in the error measure?
    const bool do_hessians = false;         // consider Hessians in the error measure?
    Dist_local_AbsDiff distance_local;
    Dist_global_RMS distance_global;
    double (*objective)(const double, const double, const double) = &L2_like_objective;

    
    /* ******************************************************************************
    * gamma setup. your choices:
    *   (1) EITHER specify options for an automated search
    *       meaning: perform num_regular computations of the error on the gamma interval [gamma_min,gamma_max]
    *       followed by num_bisec bisections with factor bisection_factor of the best found gamma sub-interval
    *       this option ignores the gamma value given by the command line.
    *       if you do not want this option, set search_gamma false.
    */
    const double gamma_min = 0.1, gamma_max = 4;
    const int num_regular = 10;
    const int num_bisec = 4;
    const double bisection_factor = 1.5;    // recommended value is 1.5
    bool search_gamma = true;               // will search for gamma on first non-zero radius
    double gamma = 0./0.;                   // NaN as default, force the user to decide
    
    /*
    *   (2) OR set gamma given as an argument by the command line.
    *       if you do not want this option, comment the following four lines of code (not counting comment lines)
    */
//     gamma = atof(argv[6]);
//     interpolation.SetGamma(gamma);
//     printf("$$$ Set gamma = %lf\n", gamma);
//     search_gamma = false;

    /* ******************************************************************************
     * read validation directions and data
     */
    double ** matrix_directions_validation = UTILITY::ReadMatrix( &read_rows, &read_cols, argv[3]);
    const int N_validation = read_rows;
    assert_msg( read_cols == D, "Error: dimension of validation directions must match that of the training directions\n");
    printf("$$$ Read validation directions: %i x %i\n", N_validation, D);

    double ** matrix_data_validation = UTILITY::ReadMatrix(&read_rows, &read_cols, argv[4]);
    assert_msg( read_rows == N_validation, "Error: number of validation data rows must equal number of validation directions\n");
    const int R_validation = read_cols;
    printf("$   Read validation data: %i x %i\n", N_validation, R_validation);
    
    /* ******************************************************************************
    * run the interpolation and write gnuplot compatible output file
    */
    // header line as comment (beginning with hash)
    FILE * outfile = fopen(argv[5],"w");    // WARNING: overwrites output file. replace "w" by "a" to append instead.
    fprintf(outfile,"# N_training = %i, R_training = %i, N_validation = %i, R_validation = %i, gamma = %lf\n#training directions  : %s\n#training data        : %s\n#validation directions: %s\n#validation data      : %s\n#radius                         L2error\n",
            N_training, R_training, N_validation, R_validation, gamma, argv[1], argv[2], argv[3], argv[4]);
    fclose(outfile);
    
    // create a validation data object data_validation for each validation radius separately, and run the interpolation on each of these radii.
    double r_validation[R_validation], L2errors[R_validation];
    for(int i_radius=0; i_radius<R_validation; i_radius++)
    {
        r_validation[i_radius] = double(i_radius)/double(R_validation-1);
        Data data_validation(D, 1*N_validation);
        int i_point = 0;
        for(int i_direction=0; i_direction<N_validation; i_direction++)
        {
            for(int i_dimension=0; i_dimension<D; i_dimension++)
                data_validation.SafeAccess_coord(i_point,i_dimension) = matrix_directions_validation[i_direction][i_dimension]*r_validation[i_radius];
            data_validation.SafeAccess_values(i_point) = matrix_data_validation[i_direction][i_radius];
            i_point++;
        }
        // perform gamma search if appropriate
        if(search_gamma && r_validation[i_radius]>1e-16)
        {
            search_gamma = false;
            gamma = OptimizeGamma(
                        &interpolation,                                                 
                        &data_validation,                                                   
                        &distance_local,                                               
                        &distance_global,                                              
                        objective,                                                     
                        gamma_min, gamma_max,                                          
                        num_regular, num_bisec,                                        
                        bisection_factor,                                              
                        do_gradients, do_hessians   );
            interpolation.SetGamma(gamma);
            printf("$$$ Set gamma = %lf\n", gamma);
        }
    
        printf(" testing on all validation directions at radius %5.3f (%3i/%3i) ... ", r_validation[i_radius], i_radius, R_validation); fflush(stdout);
        TestOnSet test_object(&interpolation, &data_validation, &distance_local, &distance_global, objective);
        L2errors[i_radius] = test_object.run(do_gradients, do_hessians);
        // save results after each radius
        FILE * outfile = fopen(argv[5],"a");
        fprintf(outfile, "%28.20e   %28.20e\n", r_validation[i_radius], L2errors[i_radius]);
        fclose(outfile);
    }

    
    /* ******************************************************************************
    * close the program, free memory 
    */
    free_matrix(matrix_data_training, N_training);
    free_matrix(matrix_directions_validation, N_validation);
    free_matrix(matrix_data_validation, N_validation);
    
    printf("\n\n$$$ Program finished regularly\n\n\n");
    
    return 0;
}

