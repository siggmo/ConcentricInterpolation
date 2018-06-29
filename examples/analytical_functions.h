#ifndef __ANALYTICAL_FUNCTIONS_H__
#define __ANALYTICAL_FUNCTIONS_H__

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

#include "util.h"

// general function prototype
class AnalyticalFunction{
public:
    AnalyticalFunction(const int a_D):D(a_D){}
    virtual void evaluate(
        const double * point,   //  [IN] evaluation point
        double & value,         // [OUT] resulting value
        double * gradient = 0,  // [OUT] resulting gradient, optional
        double * hessian = 0    // [OUT] resulting hessian, optional
                ) = 0;
                
    virtual void PrintName() = 0;
    int GetDim() const {return D;};

protected:
    const int D;
};




// constant function 
class ConstantFunction : public AnalyticalFunction{
public:
    ConstantFunction(const int a_D, const double a_const_value);
    void PrintName(){printf("ConstantFunction on %i dimensions, value is %le\n", D, const_value);}
    void evaluate(const double * point, double & value, double * gradient = 0, double * hessian = 0);
    

private:
    const double const_value;
};


// linear function f = a^T * x
class LinearFunction : public AnalyticalFunction{
public:
    LinearFunction(double * a_a, const int a_D);
    void PrintName(){printf("LinearFunction on %i dimensions\n", D);}
    void evaluate(const double * point, double & value, double * gradient = 0, double * hessian = 0);
    

private:
    double * const a;
};

// quadratic function f = x^T * A * x * 0.5
class QuadraticFunction : public AnalyticalFunction{
public:
    QuadraticFunction(double ** a_A, const int a_D);
    void PrintName(){printf("QuadraticFunction on %i dimensions\n", D);}
    void evaluate(const double * point, double & value, double * gradient = 0, double * hessian = 0);

private:
    double ** const A;
};

// pseudo plasticity material law (hyperelasticity mimicking von Mises plasticity, only physical for proportional loading)
// input space is the 6-dimensional space of infinitesimal strain tensors eps:
//   eps[0] = eps_xx, eps[1] = eps_yy, eps[2] = eps_zz,
//   eps[3] = eps_xy*sqrt(2), eps[4] = eps_xz*sqrt(2), eps[5] = eps_yz*sqrt(2)
// stress and stiffness are parametrized accordingly
class PseudoPlasticity : public virtual AnalyticalFunction
{
private:
//     const int D = 6;
    const double sqrt23;
    double s0,         //! yield stress [MPa],
           h,          //! hardening modulus [MPa],
           effMod,     //! effective modulus [MPa] for eps > epsc
           twoG,       //! twice the shear modulus [MPa]
           lambda,     //! Lamé's first parameter [MPa]
           K,          //! bulk modulus [MPa]
           E,          //! Young's modulus [MPa]
           nu,         //! Poisson ratio
           epsc;       //! "epsilon critical": plasticity yield limit
    
    //! utility function computing the effective modulus for eps > epsc
    void ComputeEffMod();
    //! utility function computing the yield strain from the yield stress and the shear modulus
    void ComputeEpsc();
    double DeviNorm(const double * x);
public:
    //! constructor with default values similar to Aluminium
    PseudoPlasticity( const double a_E = 75e3, const double a_nu = 0.3, const double a_s0 = 100., const double a_h = 0. );
    
    void PrintName(){printf("PseudoPlasticity on %i dimensions\n", D);}
    
    //! evaluation function
    void evaluate( const double * eps, double & energy, double * stress = 0, double *stiffness = 0 );
    
    //! functions providing access to material parameters
    void SetBulkModulus( const double a_K );
    void SetYieldStress( const double a_s0 );
    void SetHardeningModulus( const double a_h );
    void SetElasticModuli( const double a_E, const double a_nu );  //! computes twoG, K, lambda, epsc and effMod from given E and nu
};


#endif
