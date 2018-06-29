#ifndef __CREATE_DATA_H__
#define __CREATE_DATA_H__
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

#include "data.h"
#include "analytical_functions.h"

using namespace UTILITY;

/*! \brief Create data object from coordinates and function
 * 
 * Stores the coordinates and the corresponding function data into the return object.
 * 
 * This can be used as a template to create data from other sources, e.g. external programs, Monte Carlo,
 * experiments, etc. Modifications should be straight-forward.
 * 
 * There is an overloaded type accepting directions + radii:
 * 
 * \see \CreateData TODO
 */
Data * CreateData(  double ** coordinates,      /*!< [in] P-by-D matrix containing the coordinates of the points */
                    const int P,                /*!< [in] number of points */
                    const int D,                /*!< [in] number of dimensions */
                    AnalyticalFunction & func   /*!< [in] analytical pre-/user-defined function that will be evaluated at the coordinates */
                 );

/*! \brief Create data object from directions, radii, and function
 * 
 * Stores the tensor product of directions and radii as coordinates,
 * and the corresponding function data into the return object.
 * 
 * The non-overloades version of this function is used:
 * 
 * \see \CreateData TODO
 */
Data * CreateData(  double ** directions,       /*!< [in] N-by-D matrix containing directions */
                    const int N,                /*!< [in] number of directions */
                    const int D,                /*!< [in] number of dimensions */
                    double * radii,             /*!< [in] array of length R containing the radii, same radii for every direction */
                    const int R,                /*!< [in] number of radii */
                    AnalyticalFunction & func   /*!< [in] analytical pre-/user-defined function that will be evaluated at the coordinates */
                 );

/*! \brief Create data object from file-read directions, radii, and function
 * 
 * Stores the tensor product of directions and radii as coordinates,
 * and the corresponding function data into the return object.
 * 
 * The non-overloades version of this function is used:
 * 
 * \see \CreateData TODO
 */
Data * CreateData(  char * DIRECTIONS_FN,       /*!< [in] name of file containing the directions */
                    double * radii,             /*!< [in] array of length R containing the radii, same radii for every direction */
                    const int R,                /*!< [in] number of radii */
                    AnalyticalFunction & func   /*!< [in] analytical pre-/user-defined function that will be evaluated at the coordinates */
                 );

/*! \brief Create training data object from file-read directions, radii, and function
 * 
 * Stores the directions, the radii (same radii for all directions),
 * the corresponding function values and the radial derivative of the
 * function to the DataTrain object pointed to by the return.
 * 
 * \see \CreateData TODO
 */
DataTraining * CreateDataTraining( 
                    char * DIRECTIONS_FN,       /*!< [in] name of file containing the directions */
                    double * radii,             /*!< [in] array of length R containing the radii, same radii for every direction */
                    const int R,                /*!< [in] number of radii */
                    AnalyticalFunction & func   /*!< [in] analytical pre-/user-defined function that will be evaluated at the coordinates */
                 );

#endif // __CREATE_DATA_H__
