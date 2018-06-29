#include "analytical_functions.h"

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

using namespace UTILITY;

ConstantFunction::ConstantFunction(const int a_D, const double a_const_value):
AnalyticalFunction(a_D),
const_value(a_const_value)
{
    printf("$   created "), PrintName();
}

void ConstantFunction::evaluate(const double * point, double & value, double * gradient, double * hessian)
{
    value = const_value;
    if( gradient!=0 )
        memset(gradient, 0, sizeof(double)*D);
    if( hessian!=0 )
        memset(hessian, 0, sizeof(double)*D*D);
    
}

LinearFunction::LinearFunction(double * a_a, const int a_D):
AnalyticalFunction(a_D),
a(a_a)
{
	printf("$   created %i-dimensional LinearFunction object with a = [ ", D);
	for(int d=0; d<D; d++)
		printf("%lf ", a[d]);
	printf("]\n");
}

void LinearFunction::evaluate(const double * point, double & value, double * gradient, double * hessian)
{
	value = VecVecMul(a, point, D);
	if( gradient!=0 )
		for(int d=0; d<D; d++)
			gradient[d] = a[d];
	if( hessian!=0 )
		memset(hessian, 0, sizeof(double)*D*D);
}

QuadraticFunction::QuadraticFunction(double ** a_A, const int a_D):
AnalyticalFunction(a_D),
A(a_A)
{
    printf("$   created "), PrintName();
}

void QuadraticFunction::evaluate(const double * point, double & value, double * gradient, double * hessian)
{
	value = 0.;
    for(int d1=0; d1<D; d1++)
        for(int d2=0; d2<D; d2++)
            value           += point[d1] * A[d1][d2] * point[d2] * 0.5;

    if( gradient != 0 )
        for(int d1=0; d1<D; d1++)
        {
            gradient[d1] = 0.;
            for(int d2=0; d2<D; d2++)
                gradient[d1]    += A[d1][d2] * point[d2];
        }
        
    if( hessian != 0 )
        for(int d1=0; d1<D; d1++)
            for(int d2=0; d2<D; d2++)
                hessian[d1*D+d2] = A[d1][d2];
}


PseudoPlasticity::PseudoPlasticity( const double a_E, const double a_nu, const double a_s0, const double a_h ):
    AnalyticalFunction(6),
	sqrt23(sqrt(2./3.))
{
    printf("$   created 6-dimensional PseudoPlasticity model\n");
	SetElasticModuli( a_E, a_nu );
	SetYieldStress( a_s0 );
	SetHardeningModulus( a_h );
}
void PseudoPlasticity::ComputeEpsc() {
	epsc = s0*sqrt23/twoG;
}
void PseudoPlasticity::ComputeEffMod() {
	effMod = twoG / ( twoG + 2./3. * h ) ;
}
void PseudoPlasticity::SetBulkModulus( const double a_K ) {
	assert_msg(a_K>0., "ERROR: PseudoPlasticity::SetBulkModulus(): Bulk modulus must be positive\n");
	K = a_K;
}
void PseudoPlasticity::SetYieldStress( const double a_s0 ) {
	assert_msg(a_s0>0., "ERROR: PseudoPlasticity::SetYieldStress(): Yield stress must be positive\n");
	s0 = a_s0;
	ComputeEpsc();
}
void PseudoPlasticity::SetHardeningModulus( const double a_h ) {
	assert_msg(a_h>=0., "ERROR: PseudoPlasticity::SetHardeningModulus(): Hardening modulus must be >= 0.\n");
	h		= a_h;
	ComputeEffMod();
}

void PseudoPlasticity::SetElasticModuli( const double a_E, const double a_nu ) {
	assert_msg(a_E>0., "ERROR: PseudoPlasticity::SetElasticModuli(): E must be positive\n");
	assert_msg((a_nu>0.) && (a_nu<0.5), "ERROR: PseudoPlasticity::SetElasticModuli(): nu must be in (0.0; 0.5)\n");
	E		= a_E;
	nu		= a_nu;
	twoG	= E/(1.+nu);
	K		= E/(3.*(1.-2.*nu));
	lambda	= twoG*nu/(1.-2.*nu);
	ComputeEpsc();
	ComputeEffMod();

}

// the energy is now called psi, stress is dpsi, stiffness is ddpsi
void PseudoPlasticity::evaluate( const double * const eps, double & psi, double * dpsi, double * ddpsi )
{
	const double eps_eq = DeviNorm(eps), tr_eps = eps[0]+eps[1]+eps[2];
	psi = 0.5*K*tr_eps*tr_eps;
	if( eps_eq < epsc ) {
		/*elasticity*/
		// energy
		psi += 0.5 * twoG * eps_eq * eps_eq;
		
		// stress
        if( dpsi != 0 )
        {
            for( int i=0;i<6; i++) dpsi[i] = twoG*eps[i];
            dpsi[0] += lambda*tr_eps;
            dpsi[1] += lambda*tr_eps;
            dpsi[2] += lambda*tr_eps;
        }
		
		// stiffness
        if( ddpsi != 0 )
        {
            ddpsi[ 0]=lambda+twoG; ddpsi[ 1]=lambda;      ddpsi[ 2]=lambda;      ddpsi[ 3]=ddpsi[ 4]=ddpsi[ 5]=0.;
            ddpsi[ 6]=lambda;      ddpsi[ 7]=lambda+twoG; ddpsi[ 8]=lambda;      ddpsi[ 9]=ddpsi[10]=ddpsi[11]=0.;
            ddpsi[12]=lambda;      ddpsi[13]=lambda;      ddpsi[14]=lambda+twoG; ddpsi[15]=ddpsi[16]=ddpsi[17]=0.;
            ddpsi[18]=ddpsi[19]=ddpsi[20]=0.;                                    ddpsi[21]=twoG;         ddpsi[22]=ddpsi[23]=0.;
            ddpsi[24]=ddpsi[25]=ddpsi[26]=0.;                                    ddpsi[27]=0.;   ddpsi[28]=twoG; ddpsi[29]=0.;
            ddpsi[30]=ddpsi[31]=ddpsi[32]=0.;                                    ddpsi[33]=ddpsi[34]=0.;         ddpsi[35]=twoG;
        }
        
	} /* endif "elastic" */
	else
	{ /* "plasticity" */
		// energy
		psi += 0.5*twoG*epsc*epsc + ( sqrt23 * s0  + 2./3. * 0.5 * effMod * h * ( eps_eq - epsc ) ) * ( eps_eq - epsc );
		
		const double eps_d[6] = { eps[0] - tr_eps/3., eps[1]-tr_eps/3. ,eps[2]-tr_eps/3., eps[3], eps[4], eps[5] };
		const double f = effMod * ( 2./3.*h + sqrt23 * s0 / eps_eq );
		
		// stress (i.e. first derivative)
        if( dpsi != 0 )
        {
            for( int i=0;i<3; i++) dpsi[i] = f*eps_d[i]  +  K*tr_eps;
            for( int i=3;i<6; i++) dpsi[i] = f*eps_d[i];
        }
		
		// stiffness (i.e. second derivative)
		if( ddpsi != 0 )
        {
            const double lambda_p = K - 1./3. * f;
            ddpsi[ 0]=lambda_p+f; ddpsi[ 1]=lambda_p;   ddpsi[ 2]=lambda_p;      ddpsi[ 3]=ddpsi[ 4]=ddpsi[ 5]=0.;
            ddpsi[ 6]=lambda_p;   ddpsi[ 7]=lambda_p+f; ddpsi[ 8]=lambda_p;      ddpsi[ 9]=ddpsi[10]=ddpsi[11]=0.;
            ddpsi[12]=lambda_p;   ddpsi[13]=lambda_p;   ddpsi[14]=lambda_p+f;    ddpsi[15]=ddpsi[16]=ddpsi[17]=0.;
            ddpsi[18]=ddpsi[19]=ddpsi[20]=0.;                                 ddpsi[21]=f;     ddpsi[22]=ddpsi[23]=0.;
            ddpsi[24]=ddpsi[25]=ddpsi[26]=0.;                                 ddpsi[27]=0.; ddpsi[28]=f; ddpsi[29]=0.;
            ddpsi[30]=ddpsi[31]=ddpsi[32]=0.;                                 ddpsi[33]=ddpsi[34]=0.;     ddpsi[35]=f;
            
            const double g = - effMod *sqrt23*s0/(eps_eq*eps_eq*eps_eq);
            
            for(int i=0;i<6;i++)
                for(int j=0;j<6;j++)
                    ddpsi[i*6+j] += g*eps_d[i]*eps_d[j];
        }
	}
}

double PseudoPlasticity::DeviNorm(const double * x)
{
	const double p = (x[0]+x[1]+x[2])/3.;
	const double q = ( (x[0]-p)*(x[0]-p)+(x[1]-p)*(x[1]-p)+(x[2]-p)*(x[2]-p) )
						+ x[3]*x[3] + x[4]*x[4] + x[5]*x[5];
	return (( q < 0. ) ? 0.:sqrt(q)); /* q may become less than zero due to round off!!! */
	
}
