#include "concentric_interpolation.h"   // the main program
#include "create_data.h"                // create data sets
#include "data_util.h"                  // utility functions concerning data. includes the kernel parameter optimization
#include "objective_functions.h"        // objective functions with respecto to which the kernel parameter may be optimized

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


int main(int argc, char* argv[])
{
    /* here is the input data.
    *  for convenience, this could in general be read in from a parameter file during runtime.
    *  however this is just a minimum working example.
    *  on Linux, filenames must be the absolute path, i.e. relative to root directory "/"
    */
    
    // filenames
    if(argc!=4)
    {
        fprintf(stderr, "USAGE: [this] [filename training dirs] [filename target dirs] [filename of output file]\n");
        exit(-1);
    }
    
    char    DIRS_FN_TRAINING[512];  // name of textfile containing the training directions, e.g. created by the MATLAB program. each row is one direction.
    char    DIRS_FN_TARGET[512];    // name of textfile containing the target directions, e.g. created by the MATLAB programm. each row is one direction
    char    FN_OUTPUT[512];         // name of the output file which can be loaded by the GNUPLOT script
    sprintf(DIRS_FN_TRAINING, "%s", argv[1]);
    sprintf(DIRS_FN_TARGET, "%s", argv[2]);
    sprintf(FN_OUTPUT, "%s", argv[3]);
    // symmetry flag
    const bool SYM = false;         // use symmetrized interpolation scheme? only set this true if the training directions were created accordingly
    // radii. here, one single radii set is applied to all directions of one direction set. this could be generalized.
    double  radii_training[]        = {1.,10.},                   // training radii, i.e. factors with which the training directions will be multiplied to give the training points
            radii_target[]          = {5.};                  // target radii, i.e. factors with which the target directions will be multiplied to give the target points

    // numbers of radii R
    const int R_training    = sizeof(radii_training)/sizeof(double),
              R_target  = sizeof(radii_target)/sizeof(double);
    
    // outputs concerning radii
    printf("$$$ %02i training radii   = ", R_training);
    if( R_training<=10 )
        UTILITY::print_matrix(radii_training, 1, R_training);
    else
        { printf("(only first 10) "); UTILITY::print_matrix(radii_training, 1, 10);}
    printf("$$$ %02i target radii   = ", R_target);
    if( R_target<=10 )
        UTILITY::print_matrix(radii_target, 1, R_target);
    else
        { printf("(only first 10) "); UTILITY::print_matrix(radii_target, 1, 10);}
        
    /* ******************************************************************************
     * decide which analytically given function should be approximated. the candidates are given in analytical_functions.h
    *  comment all but the desired occurence.
    */
    
    // linear function f = a^T * x  (where a is a constant vector)
    const int D = 3;
    double a[D]; for(int d=0; d<D; d++) a[d]=0.; a[2] = 1.;
    LinearFunction target_function(a,D);
    
    //  quadratic function f = x^T * A * x * 0.5  (where A is a constant matrix)
//     double ** A = alloc_matrix(D,D); // allocate the memory for A, then define its components
//     double a = 2./3.;
//  A[0][0] = 2*a;    A[0][1] = -a;    A[0][2] = -a;    A[0][3] = 0;    A[0][4] = 0;    A[0][5] = 0;
//  A[1][0] = -a;    A[1][1] = 2*a;    A[1][2] = -a;    A[1][3] = 0;    A[1][4] = 0;    A[1][5] = 0;
//  A[2][0] = -a;    A[2][1] = -a;    A[2][2] = 2*a;    A[2][3] = 0;    A[2][4] = 0;    A[2][5] = 0;
//  A[3][0] = 0;    A[3][1] = 0;    A[3][2] = 0;    A[3][3] = 3*a;    A[3][4] = 0;    A[3][5] = 0;
//  A[4][0] = 0;    A[4][1] = 0;    A[4][2] = 0;    A[4][3] = 0;    A[4][4] = 3*a;    A[4][5] = 0;
//  A[5][0] = 0;	A[5][1] = 0;    A[5][2] = 0;    A[5][3] = 0;    A[5][4] = 0;    A[5][5] = 3*a;

    
    //  pseudo plasticity material law (hyperelasticity mimicking von Mises plasticity, only physical for proportional loading)
    //  default constructor contains parameters similar to Aluminium
//     PseudoPlasticity target_function( 75e3, 0.3, 100., 0. );

    
    /* [end of analytical function block] */
    
    
    /* ******************************************************************************
     * create training data and target data
     */
    printf("$$$ creating training data\n");
    DataTraining * data_training = CreateDataTraining(DIRS_FN_TRAINING, radii_training, R_training, target_function);

    printf("$$$ creating target data\n");
    Data * data_target   = CreateData(DIRS_FN_TARGET, radii_target, R_target, target_function);

    
    UTILITY::assert_msg(data_training->D == data_target->D, "ERROR: training data and target data must have the same dimensions\n");
    
    /* [end of data block] */
    
    
    /* ******************************************************************************
    * initialize and setup the interpolation
    */
    printf("$$$ Setup the interpolation\n");
    // initialize
    ConcentricInterpolation interpolation( SYM );
    
    // data setup
    BuildInterpolationFromTrainingData( interpolation, data_training);

    
    /* ******************************************************************************
    * specify error measure and objective function
    */
    const bool do_gradients = true, do_hessians = false;
    Dist_RelDiff_local distance_local;
    Dist_Mean_global distance_global;
//     Dist_AbsDiff_local distance_local;
//     Dist_RMS_global distance_global;
    double (*objective)(const double, const double, const double) = &H1_like_objective;
    
    
    /* ******************************************************************************
    * gamma setup
    *
    *      your choice: either (1) give specific gamma or (2) do gamma optimization to best fit the just created target data. comment the undesired and uncomment the desired.
    * ***
    *      option (1) :
    */
 const double gamma = 2.56156716;
 interpolation.SetGamma(gamma);
 printf("$$$ set gamma = %lf\n", gamma);
    /* ***
     *      option (2):
     */
//     const double gamma_min = 0.5, gamma_max = 4;
//     const int num_regular = 10;
//     const int num_bisec = 8;
//     const double bisection_factor = 1.5; //ATTENTION recommended value is 1.5
//     const double gamma = OptimizeGamma(
//                                         &interpolation,                                                 
//                                         data_target,                                                   
//                                         &distance_local,                                               
//                                         &distance_global,                                              
//                                         objective,                                                     
//                                         gamma_min, gamma_max,                                          
//                                         num_regular, num_bisec,                                        
//                                         bisection_factor,                                              
//                                         do_gradients, do_hessians   );

    
    /* ******************************************************************************
    * run the interpolation on the target points and compare against target data  
    */
    TestOnSet test_object(&interpolation, data_target, &distance_local, &distance_global, objective);
    printf("error = %le\n", test_object.run(do_gradients, do_hessians));
 
 
    
    /* ******************************************************************************
    * output interpolation results to gnuplottable text file
    */
    WriteDataGnuplot( test_object.data_interpolation, FN_OUTPUT);
    
    
    
    
    /* ******************************************************************************
    * close the program, free memory 
    */
    delete data_training;
    delete data_target;
    // if quadratic function has been set up, free matrix A
//     free_matrix(A,D);
    
    printf("\n\n$$$ Program finished regularly\n\n\n");
    
    return 0;
}

