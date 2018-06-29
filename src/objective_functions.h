#ifndef __OBJECTIVE_FUNCTIONS_H__
#define __OBJECTIVE_FUNCTIONS_H__

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

//! NOTE: add your very own objective function (incorporating e.g. physics-motivated weights)

//! Objective function that treats functions value, gradient and Hessian equally
inline double H2_like_objective(const double dist_values, const double dist_gradients, const double dist_hessians)
{
    return dist_values+dist_gradients+dist_hessians;
}

//! Objective function that considers the gradient and the function value
inline double H1_like_objective(const double dist_values, const double dist_gradients, const double dist_hessians)
{
    return dist_values + dist_gradients;
}

//! Objective function that considers only the L2 norm of the function
inline double L2_like_objective(const double dist_values, const double dist_gradients, const double dist_hessians)
{
    return dist_values;
}

#endif // define __OBJECTIVE_FUNCTIONS_H__
