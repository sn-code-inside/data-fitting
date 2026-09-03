/*****************************************************************
 *
 * File........: functions.c
 * Function....: model functions and their derivatives
 * Author......: Tilo Strutz
 * last changes: 07.05.2008, 26.10.2009, 18.02.2010, 01.01.2011
 *              25.02.2021, 05.11.2025, 11.05.2026
 *
 * LICENCE DETAILS: see software manual
 *	free academic use
 *  cite source as
 * "Strutz, T.: Data Fitting and Uncertainty. 3rd edition,
 *   Springer Nature, 2026"
 *
*****************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
//#include <assert.h>
#include "macros.h"
#include "functions.h"
#include "defines.h"
#include "tinyexpr.h"

#define DEG2RAD	M_PI/180.
#ifndef WIN32
#include <sys/time.h>
#else
#include <time.h>
#define random rand
#endif


extern double x; /* declared in fitting.c */

/* for reasons of compatibility all derivative functions
 * have the same parameter list
 * xxx_deriv( double (*funct)(int,double*,double*),
 *				unsigned int i, unsigned int j, unsigned int M, double *cond, double *a)
 *
 * However, the only function requiring (*funct) is the
 * generic function f_deriv()
 * M is required by flin2_deriv()
 */

 /*---------------------------------------------------------------
	*	fhyperbel()
	*  f(x|a) = a1 / (x+a2) 
	*--------------------------------------------------------------*/
double
fhyperbel( unsigned int i, double* cond, double* a )
{
	if (a[1] + cond[i] > 0.0 && a[1] + cond[i] < TOL)
		return a[0] / TOL ;
	else if (a[1] + cond[i] < 0.0 && a[1] + cond[i] > TOL)
		return -a[0] / TOL ;
	return a[0] / (cond[i] + a[1]);
}

/*---------------------------------------------------------------
 *	fexpBaseX()
 *  f(x|a) = (a1*x)^a2 + a3
 *--------------------------------------------------------------*/
double
fexpBaseX( unsigned int i, double *cond, double *a )
{
	if (a[0] < 0.0)
		a[0] = TOL;
	if (a[1] < 0.0)
		a[1] = TOL;
	return pow( a[0] * cond[i], a[1] ) + a[2];
}
/*---------------------------------------------------------------
 *	fexpBaseX_deriv()
 *  f(x|a) = (a1*x)^a2 + a3
 *--------------------------------------------------------------*/
double
fexpBaseX_deriv( double(*funct)(unsigned int, double*, double*),
	unsigned int i, unsigned int j, unsigned int M, double* cond, double* a )
{
	(void)funct;
	(void)M;
	if (j == 0)
	{
		if (fabs(a[1]) < TOL) /* exponent may not be zero */
			a[1] = TOL;
		if (a[0] < TOL)
			return 0.0;/* zero base is not allowed */
		if (fabs(cond[i]) < TOL) 
			return 0.0;/* zero base is not allowed */
		return a[1] * pow( a[0] * cond[i], a[1] ) /a[0];	/* derivation of a1 */
	}
	else if (j == 1)
	{
		if (fabs( a[1] ) < TOL) /* exponent may not be zero */
			a[1] = TOL;
		if (a[0] < 0.0)
			a[0] = TOL;
		if (cond[i] < TOL) 
			return 0.0;
		return pow( a[0] * cond[i], a[1] ) * log( a[0] * cond[i] );/* derivation of a2 */
	}
	else
		return 1.;/* derivation of a3 */
}

/*---------------------------------------------------------------
 *	fconstant_deriv()
 *--------------------------------------------------------------*/
double
fconstant_deriv( double( *funct )(unsigned int, double*, double*),
	unsigned int i, unsigned int j, unsigned int M, double *cond, double *a )
{
	(void)funct; /* suppresses warnings */
	(void)M;
	(void)a;
	(void)cond;
	(void)i;
	(void)j;
	return 1.0;	/* derivation of a1 */
}

/*---------------------------------------------------------------
 *	flin_deriv()
 *  f(x|a) = a1 + SUM_j a_j * x_j
 *--------------------------------------------------------------*/
double
flin_deriv( double( *funct )(unsigned int, double*, double*),
	unsigned int i, unsigned int j, unsigned int M, double *cond, double *a )
{
	(void)funct;
	(void)M;
	(void)i;
	(void)a;
	/* i ... number of current observation */
	/* j = 0,1,2,...  ... number of parameter a_j*/
	if (j == 0)
		return 1.;	/* derivation of a1 */
	else
		/* all conditions x_j are stored in a single array cond[.] */
		return cond[(M - 1) * i + j - 1];
}

/*---------------------------------------------------------------
 *	flin2_deriv()
 *  f(x|a) = SUM_j a_j * x_j
 *--------------------------------------------------------------*/
double
flin2_deriv( double( *funct )(unsigned int, double*, double*),
	unsigned int i, unsigned int j, unsigned int M, 
	double *cond, double *a )
{
	(void)a;
	(void)funct;
	/* all conditions x_j are stored in a single array cond[.] */
	return cond[M * i + j];
}

/*---------------------------------------------------------------
 *	fcosineNonLin()
 * f(x|a) = a1 + a2 * cos( x - a3)
 *--------------------------------------------------------------*/
double fcosineNonLin( unsigned int i, double* cond, double* a )
{
	return a[0] + a[1] * cos( (cond[i] - a[2]) * DEG2RAD );
}

/*---------------------------------------------------------------
 *	fcosineNonLin_deriv()
 * f(x|a) = a1 + a2 * cos( x - a3)
 *--------------------------------------------------------------*/
double
fcosineNonLin_deriv( double(*funct)(unsigned int, double*, double*),
	unsigned int i, unsigned int j, unsigned int M, double* cond, double* a )
{
	(void)funct;
	(void)M;
	if (j == 0)
		return 1.;	/* derivation of a1 */
	else if (j == 1)
		return cos( (cond[i] - a[2]) * DEG2RAD );
	else
		return a[1] * sin( (cond[i] - a[2]) * DEG2RAD );
}

/*---------------------------------------------------------------
 *	fcosineLin_deriv()
 * f(x|a) = a1 + a2 * cos( x) + a3 * sin( x)
 *--------------------------------------------------------------*/
double
fcosineLin_deriv( double(*funct)(unsigned int, double*, double*),
	unsigned int i, unsigned int j, unsigned int M, double* cond, double* a )
{
	(void)funct;
	(void)M;
	(void)a;
	if (j == 0)
		return 1.;	/* derivation of a1 */
	else if (j == 1)
		return (cos( cond[i] * DEG2RAD ));
	else
		return (sin( cond[i] * DEG2RAD ));
}

/*---------------------------------------------------------------
 *	fcosine()
 * f(x|a) = a1 + a2 * cos( a3*x + a4)
 * a3 = 2*Pi*f
 *--------------------------------------------------------------*/
double
fcosine( unsigned int i, double *cond, double *a )
{
	if (a[3] >= 2*M_PI) a[3] -= 2 * M_PI;
	else if (a[3] < 0.) a[3] += 2 * M_PI;

	return a[0] + a[1] * cos( (a[2] * cond[i] + a[3]) ); //* DEG2RAD
}

/*---------------------------------------------------------------
 *	fcosine_deriv()
 * f(x|a) = a1 + a2 * cos( a3*x + a4)
 *--------------------------------------------------------------*/
double
fcosine_deriv( double( *funct )(unsigned int, double*, double*),
	unsigned int i, unsigned int j, unsigned int M, double *cond, double *a )
{
	(void)funct;
	(void)M;
	if (j == 0)
		return 1.;	/* derivation of a1 */
	else if (j == 1)
		return cos( (a[2] * cond[i] + a[3]) ); // * DEG2RAD
	else if (j == 2)
		return -a[1] * cond[i] * sin( cond[i] * a[2] + a[3]  );// * DEG2RAD
	else
		return -a[1] * sin( (a[2] * cond[i] + a[3])  ); //* DEG2RAD
}

/*---------------------------------------------------------------
 *	fcosine_trend()
 * f(x|a) = a1 + a2 * x + a3 * cos( x - a4)
 *--------------------------------------------------------------*/
double
fcosine_trend( unsigned int i, double *cond, double *a )
{
	return a[0] + a[1] * cond[i] + a[2] * cos( cond[i] - a[3] );
}

/*---------------------------------------------------------------
 *	fcosine_trend_deriv()
 * f(x|a) = a1 + a2 * x + a3 * cos( x - a4)
 *--------------------------------------------------------------*/
double
fcosine_trend_deriv( double( *funct )(unsigned int, double*, double*),
	unsigned int i, unsigned int j, unsigned int M, double *cond, double *a )
{
	(void)funct;
	(void)M;
	if (j == 0)
		return 1.;	/* derivation of a1 */
	else if (j == 1)
		return cond[i];
	else if (j == 2)
		return cos( cond[i] - a[3] );
	else
		return (a[2] * sin( cond[i] - a[3] ));
}

/*---------------------------------------------------------------
 *	ftrigonometric2()
 * f(x|a) = a1 + a2*cos(a3*x-a4) + a5*cos(2*a3*x-a7)
 *--------------------------------------------------------------*/
double
ftrigonometric2( unsigned int i, double *cond, double *a )
{
	return a[0] + a[1] * cos( a[2] * cond[i] - a[3] )
		+ a[4] * cos( 2 * a[2] * cond[i] - a[5] );
}

/*---------------------------------------------------------------
 *	flogarithmic()
 * f(x|a) = a1 + a2 * log( x + a3)
 *--------------------------------------------------------------*/
double
flogarithmic( unsigned  int i, double* cond, double* a )
{
	if (cond[i] + a[2] < TOL)
		return a[0] + a[1] * log( TOL );
	else
		return a[0] + a[1] * log( cond[i] + a[2] );
}

/*---------------------------------------------------------------
 *	flogarithmic_deriv()
 * f(x|a) = a1 + a2 * log( x + a3)
 *--------------------------------------------------------------*/
double
flogarithmic_deriv( double(*funct)(unsigned int, double*, double*),
	unsigned int i, unsigned int j, unsigned int M, double* cond, double* a )
{
	(void)funct;
	(void)M;
	/* correction of values, which are outside the domain of
		 definition	*/
	if (j == 0) return 1.; 
	else if (j == 1)
	{
		if (cond[i] + a[2] < TOL)
			return log( TOL ); /* take a small value */
		else 
			return a[1] / (cond[i] + a[2]);
	}
	else
	{
		if (cond[i] + a[2] < TOL)
			return a[1] / TOL;
		else
			return a[1] / (cond[i] + a[2]);
	}
}

/*---------------------------------------------------------------
 *	fexponential()
 * f(x|a) = a1 + a2 * exp( a3 * x)
 *--------------------------------------------------------------*/
double
fexponential( unsigned int i, double *cond, double *a )
{
	/* arguments >709 will lead to INFINITY */
	/* arguments <-666 will lead to zero */
	return (a[0] + a[1] * exp( max(-666,min(666, a[2] * cond[i]) )));
}

/*---------------------------------------------------------------
 *	fexponential_deriv()
 * f(x|a) = a1 + a2 * exp( a3 * x)
 *--------------------------------------------------------------*/
double
fexponential_deriv( double( *funct )(unsigned int, double*, double*),
	unsigned int i, unsigned int j, unsigned int M, double *cond, double *a )
{
	(void)funct;
	(void)M;
	if (j == 0)
		return 1.;
	else if (j == 1)
		return (exp( a[2] * cond[i] ));
	else
	{
		return (cond[i] * a[1] * exp( max( -666, min( 666, a[2] * cond[i] ) ) ));
	}
}


/*---------------------------------------------------------------
 *	fpolynomial()
 * f(x|a) = sum_{j=1}^M aj * x^(j-1)
 *--------------------------------------------------------------*/
double fpolynomial( unsigned int i, double *cond, double *a )
{
	long j;
	double y, xj;
	/* since M ist not passed as a parameter, we assume maximal
	 * number. oveflous a[j] must be zero !
	 */
	y = a[0];
	xj = 1.;
	for (j = 1; j < M_MAX; j++)
	{
		xj *= cond[i];
		y += a[j] * xj;
	}
	return y;
}

/*---------------------------------------------------------------
 *	fpolynomial_deriv()
 * f(x|a) = sum_{j=1}^M aj * x^(j-1)
 *--------------------------------------------------------------*/
double
fpolynomial_deriv( double( *funct )(unsigned int, double*, double*),
	unsigned int i, unsigned int j, unsigned int M, double *cond, double *a )
{
	(void)funct;
	(void)M;
	(void)a;
	return pow( cond[i], (double)j );
}

/*---------------------------------------------------------------
 *	fgen_laplace()
 * f(x|a) = a1 * exp( -|x|^a2 * a3)
 *--------------------------------------------------------------*/
double
fgen_laplace( unsigned int i, double *cond, double *a )
{
	if (a[1] < 0) /* limitation of parameter space */
		a[1] = 1E-10;
	return a[0] * exp( -a[2] * pow( fabs( cond[i]), a[1]));
}

/*---------------------------------------------------------------
 *	fgen_laplace_deriv()
 * f(x|a) = a1 * exp( -|x|^a2 * a3)
 *--------------------------------------------------------------*/
double
fgen_laplace_deriv( double( *funct )(unsigned int, double*, double*),
	unsigned int i, unsigned int j, unsigned int M, double *cond, double *a )
{
	(void)funct;
	(void)M;
	if (a[1] < 0) /* limitation of parameter space */
		a[1] = 1E-10;
	if (j == 0)
	{
		/*      y = a0 * exp( -|x|^a1 * a2)
		 * dy/da0 =      exp( -a2 * |x|^a1)
		 */
		return exp( -a[2] * pow( fabs( cond[i] ), a[1] ) );
	}
	else if (j == 1)
	{
		/*      y = a0 * exp( -a2 * |x|^a1)
		 * dy/da1 = a0 * exp( -a2 * |x|^a1) * -a2 * |x|^a1 * ln|x|
		 */
		return -a[0] * exp( -a[2] * pow( fabs( cond[i] ), a[1] ) ) *
			a[2] *
			pow( fabs( cond[i] ), a[1] ) * log( fabs( cond[i] ) );
	}
	else
	{
		/*      y = a0 * exp( -a2 * |x|^a1)
		 * dy/da2 = a0 * exp( -a2 * |x|^a1) * - |x|^a1
		 */
		return -a[0] * exp( -a[2] * pow( fabs( cond[i] ), a[1] ) ) *
			pow( fabs( cond[i] ), a[1] );
	}
}

/*---------------------------------------------------------------
 *	fexpon2()
 * f(x|a) = a1 * exp( a2 * x)
 *--------------------------------------------------------------*/
double
fexpon2( unsigned int i, double *cond, double *a )
{
	return (a[0] * exp( a[1] * cond[i] ));
}

/*---------------------------------------------------------------
 *	fexpon2_deriv()
 * f(x|a) = a1 * exp( a2 * x)
 *--------------------------------------------------------------*/
double
fexpon2_deriv( double( *funct )(unsigned int, double*, double*),
	unsigned int i, unsigned int j, unsigned int M, double *cond, double *a )
{
	(void)funct;
	(void)M;
	if (j == 0)
	{
			return (exp( a[1] * cond[i] ));
	}
	else
	{
			return (cond[i] * a[0] * exp( a[1] * cond[i] ));
	}
}

/*---------------------------------------------------------------
 *	fgauss2()
 *  f(x|a) = a1 * exp( a3 * (x-a2)^2) +
 *           a4 * exp( a6 * (x-a5)^2)
 *--------------------------------------------------------------*/
double
fgauss2( unsigned int i, double *cond, double *a )
{
	double tmp1 = cond[i] - a[1];
	double tmp2 = cond[i] - a[4];
	return (a[0] * exp( a[2] * tmp1 * tmp1 ) +
		a[3] * exp( a[5] * tmp2 * tmp2 ));
}

/*---------------------------------------------------------------
 *	fgauss1()
 *  f(x|a) = a1 * exp( a3 * (x-a2)^2)
 *--------------------------------------------------------------*/
double
fgauss1( unsigned int i, double *cond, double *a )
{
	double tmp1 = cond[i] - a[1];
	return a[0] * exp( a[2] * tmp1 * tmp1 );
}

/*---------------------------------------------------------------
 *	fgauss_deriv()
 *  f(x|a) = a1 * exp( a3 * (x-a2)^2)
 *--------------------------------------------------------------*/
double
fgauss_deriv( double( *funct )(unsigned int, double*, double*),
	unsigned int i, unsigned int j, unsigned int M, double *cond, double *a )
{
	(void)funct;
	(void)M;
	double tmp1 = cond[i] - a[1];

	if (a[1] < 0) /* limitation of parameter space */
		a[1] = 1E-10;
	if (a[2] > 0) /* limitation of parameter space */
		a[2] = -1E-10;
	if (j == 0)
		/*     y = a1 * exp( a3 * (x-a2)^2)
		 *dy/da1 =      exp( a3 * (x-a2)^2)
		 */
		return exp( a[2] * tmp1 * tmp1 );
	else  if (j == 1)
		/*     y = a1 * exp( a3 * (x-a2)^2)
		 *dy/da2 =-a1 * exp( a3 * (x-a2)^2) * a3 * 2 * (x-a2)
		 */
		return -a[0] * exp( a[2] * tmp1 * tmp1 ) * a[2] *
		2 * tmp1;
	else
		/*     y = a1 * exp( a3 * (x-a2)^2)
		 *dy/da3 = a1 * exp( a3 * (x-a2)^2) * (x-a2)^2
		 */
		return  a[0] * exp( a[2] * tmp1 * tmp1 ) * tmp1*tmp1;
}

/*---------------------------------------------------------------
 *	fgauss1_deriv2()
 *  f(x|a) = a1 * exp( a3 * (x-a2)^2)
 * see also http://www.mathetools.de/differenzieren/
 * https://www.wolframalpha.com/calculators/derivative-calculator/
 *--------------------------------------------------------------*/
double
fgauss_deriv2( double( *funct )(unsigned int, double*, double*),
	unsigned int i, unsigned int j, unsigned int k, unsigned int M, 
	double *cond, double *a )
{
	(void)funct;
	(void)M;
	double tmp1 = cond[i] - a[1];

	if (a[2] > 0) /* limitation of parameter space */
		a[2] = -1E-10;

	if (j == 0)
	{
		if (k == 0)
			/*     y = a1 * exp( a3 * (x-a2)^2)
			 *dy/da1 =      exp( a3 * (x-a2)^2)
			 *d2y/d2a1 =      0
			 */
			return 0;
		else if (k == 1)
			/*     y = a1 * exp( a3 * (x-a2)^2)
			 *dy/da1 =      exp( a3 * (x-a2)^2)
			 *d2y/da1da2 = -exp( a3 * (x-a2)^2) * a3 * 2 * (x-a2)
			 */
			return -exp( a[2] * tmp1*tmp1 ) * a[2] * 2 * tmp1;
		else
			/*     y = a1 * exp( a3 * (x-a2)^2)
			 *dy/da1 =      exp( a3 * (x-a2)^2)
			 *d2y/da1da3 =  exp( a3 * (x-a2)^2) * (x-a2)^2
			 */
			return  exp( a[2] * tmp1*tmp1 ) * tmp1* tmp1;
	}
	else  if (j == 1)
	{
		if (k == 0)
			/*     y = a1 * exp( a3 * (x-a2)^2)
			 *dy/da2 =-a1 * exp( a3 * (x-a2)^2) * a3 * 2 * (x-a2)
			 *dy/da2da1 =-exp( a3 * (x-a2)^2) * a3 * 2 * (x-a2)
			 */
			return -exp( a[2] * tmp1*tmp1 ) * a[2] * 2 * tmp1;
		else if (k == 1)
			/*     y = a1 * exp( a3 * (x-a2)^2)
			 *dy/da2 =-a1 * exp( a3 * (x-a2)^2) * a3 * 2 * (x-a2)
			 *dy/da2 =-2*a1*a3 * exp( a3 * (x-a2)^2)  *   (x-a2)
			 *dy/d2a2=-2*a1*a3 * exp( a3 * (x-a2)^2)* (-1) +
								2*a1*a3 * exp( a3 * (x-a2)^2)* a3* 2*(x-a2) * (x-a2)
			 *dy/d2a2= 4*a1*a3^2 * exp( a3 * (x-a2)^2)*(x-a2)^2 +
			 *		    	2*a1*a3 * exp( a3 * (x-a2)^2)
			 */
			return 4 * a[0] * a[2] * a[2] * exp( a[2] * tmp1*tmp1 )*tmp1*tmp1 +
			2 * a[0] * a[2] * exp( a[2] * tmp1*tmp1 );
		else
			/*     y = a1 * exp( a3 * (x-a2)^2)
			 *dy/da2 =-2*a1*a3 * exp( a3 * (x-a2)^2)  *   (x-a2)
			 *dy/da2da3 =-2*a1*(x-a2)^3* a3 * exp( a3 * (x-a2)^2) -
			 *						2 * a1* (x-a2)*exp( a3 * (x-a2)^2)
			 */
			return 2 * a[0] * tmp1*tmp1*tmp1 * a[2] * exp( a[2] * tmp1*tmp1 )
			- 2 * a[0] * tmp1 * exp( a[2] * tmp1*tmp1 );
	}
	else
	{
		if (k == 0)
			/*     y = a1 * exp( a3 * (x-a2)^2)
			 *dy/da3 = a1 * exp( a3 * (x-a2)^2) * (x-a2)^2
			 *dy/da3da1 = exp( a3 * (x-a2)^2) * (x-a2)^2
			 */
			return  exp( a[2] * tmp1*tmp1 ) * tmp1* tmp1;
		else if (k == 1)
			/*     y = a1 * exp( a3 * (x-a2)^2)
			 *dy/da3 = a1 * exp( a3 * (x-a2)^2) * (x-a2)^2
			 *dy/da2da3 =-2*a1*(x-a2)^3* a3 * exp( a3 * (x-a2)^2) -
			 *						2 * a1* (x-a2)*exp( a3 * (x-a2)^2)
			 */
			return 2 * a[0] * tmp1*tmp1*tmp1 * a[2] * exp( a[2] * tmp1*tmp1 )
			- 2 * a[0] * tmp1 * exp( a[2] * tmp1*tmp1 );
		else
			/*     y = a1 * exp( a3 * (x-a2)^2)
			 *dy/da3 = a1 * exp( a3 * (x-a2)^2) * (x-a2)^2
			 *dy/d2a3= a1 * exp( a3 * (x-a2)^2) * (x-a2)^4
			 */
			return a[0] * exp( a[2] * tmp1*tmp1 ) * tmp1*tmp1 * tmp1*tmp1;
	}
}

/*---------------------------------------------------------------
 *	fcircleTLS()
 *  f(x|a) = 0 = (sqrt[(x1-a1)^2 + (x2-a2)^2^] - a3)^2
 *--------------------------------------------------------------*/
double
fcircleTLS( unsigned int i, double *cond, double *a )
{
	double tmp1, tmp2, d;
	/* two conditions per measurement */
	tmp1 = cond[2 * i] - a[0];
	tmp2 = cond[2 * i + 1] - a[1];
	d = sqrt( tmp1*tmp1 + tmp2 * tmp2 ) - a[2];
	return d;
}

/*---------------------------------------------------------------
 *	fcircleTLS_deriv()
 *  f(x|a) = 0 = (sqrt[(x1-a1)^2 + (x2-a2)^2^] - a3)^2
 *--------------------------------------------------------------*/
double
fcircleTLS_deriv( double( *funct )(unsigned int, double*, double*),
	unsigned int i, unsigned int j, unsigned int M, double *cond, double *a )
{
	(void)funct;
	(void)M;
	double b, tmp1, tmp2;
	tmp1 = a[0] - cond[2 * i];
	tmp2 = a[1] - cond[2 * i + 1];
	b = sqrt( tmp1*tmp1 + tmp2 * tmp2 );
	if (j == 0)
		return tmp1 / b;
	else  if (j == 1)
		return tmp2 / b;
	else
		return  (-1);
}
/*---------------------------------------------------------------
 *	fcircle()
 *  f(x|a) = 0 = (x1-a1)^2 + (x2-a2)^2 - a3^2
 *--------------------------------------------------------------*/
double
fcircle( unsigned int i, double *cond, double *a )
{
	double tmp1, tmp2;
	/* two conditions per measurement */
	tmp1 = cond[2 * i] - a[0];
	tmp2 = cond[2 * i + 1] - a[1];

	return tmp1 * tmp1 + tmp2 * tmp2 - a[2] * a[2];
}

/*---------------------------------------------------------------
 *	fcircle_deriv()
 * f(x|a) = 0 = (x-a1)^2 + (x2-a2)^2 - a3^2
 *--------------------------------------------------------------*/
double
fcircle_deriv( double( *funct )(unsigned int, double*, double*),
	unsigned int i, unsigned int j, unsigned int M, double *cond, double *a )
{
	(void)funct;
	(void)M;
	if (j == 0)
		return 2 * (a[0] - cond[2 * i]);
	else  if (j == 1)
		return 2 * (a[1] - cond[2 * i + 1]);
	else
		return  -2 * a[2];
}

/*---------------------------------------------------------------
 *	fcirclelin_deriv()
 * f(x|b) = x1^2 + x2^2 = b1*x1 + b2*x2 - b3
 *--------------------------------------------------------------*/
double
fcirclelin_deriv( double( *funct )(unsigned int, double*, double*),
	unsigned int i, unsigned int j, unsigned int M, double *cond, double *a )
{
	(void)funct;
	(void)M;
	(void)a;
	if (j == 0)
		return cond[2 * i];
	else  if (j == 1)
		return cond[2 * i + 1];
	else
		return  -1.;
}

/*---------------------------------------------------------------
 *	frotation()
 * 21... x = f1(u|a) = a1 + cos(a3) * u - sin(a3) * v
 *       y = f2(u|a) = a2 + sin(a3) * u + cos(a3) * v
 *--------------------------------------------------------------*/
double
frotation( unsigned int i, double *cond, double *a )
{
	if (i % 2 == 0)
	{
		/* equation for x */
		return a[0] + cos( a[2] ) * cond[i] - sin( a[2] ) * cond[i + 1];
	}
	else
	{
		/* equation for y */
		return a[1] + sin( a[2] ) * cond[i - 1] + cos( a[2] ) * cond[i];
	}
}
/*---------------------------------------------------------------
 *	frotation_deriv()
 * 21... x = f1(u|a) = a1 + cos(a3) * u - sin(a3) * v
 *       y = f2(u|a) = a2 + sin(a3) * u + cos(a3) * v
 *--------------------------------------------------------------*/
double
frotation_deriv( double( *funct )(unsigned int, double*, double*),
	unsigned int i, unsigned int j, unsigned int M, double *cond, double *a )
{
	(void)funct;
	(void)M;
	if (i % 2 == 0)
	{
		if (j == 0)
			return 1;
		else  if (j == 1)
			return 0;
		else
			return  -cond[i] * sin( a[2] ) - cond[i + 1] * cos( a[2] );
	}
	else
	{
		if (j == 0)
			return 0;
		else  if (j == 1)
			return 1;
		else
			return  cond[i - 1] * cos( a[2] ) - cond[i] * sin( a[2] );
	}
}

/*---------------------------------------------------------------
 *	fpolynom2_deriv()
 * f(x|a) = a1 + a2 * x + a3 * x*x
 *--------------------------------------------------------------*/
double
fpolynom2_deriv( double( *funct ) (unsigned int, double*, double*),
	unsigned int i, unsigned int j, unsigned int M, double *cond, double *a )
{
	(void)funct;
	(void)M;
	(void)a;
	if (j == 0)
		return 1.;
	else if (j == 1)
	{
		return cond[i];
	}
	else
	{
		return cond[i] * cond[i];
	}
}

/*---------------------------------------------------------------
 *	fpolynom3_deriv()
 * f(x|a) = a1 + a2 * x + a3 * x*x + a4 * x*x*x
 *--------------------------------------------------------------*/
double
fpolynom3_deriv( double( *funct ) (unsigned int, double*, double*),
	unsigned int i, unsigned int j, unsigned int M, double *cond, double *a )
{
	(void)funct;
	(void)M;
	(void)a;
	if (j == 0)
		return 1.;
	else if (j == 1)
	{
		return cond[i];
	}
	else if (j == 2)
	{
		return cond[i] * cond[i];
	}
	else
	{
		return cond[i] * cond[i] * cond[i];
	}
}

/*---------------------------------------------------------------
 *	fquadsurface_deriv()
 * f(x|a) = a1 + a2*x1 + a3*x1^2 + a4*x2 + a5*x2^2
 *--------------------------------------------------------------*/
double
fquadsurface_deriv( double( *funct ) (unsigned int, double*, double*),
	unsigned int i, unsigned int j, unsigned int M, double *cond, double *a )
{
	(void)funct;
	(void)M;
	(void)a;
	if (j == 0)
		return 1.;
	else if (j == 1)
	{
		return cond[2 * i];
	}
	else if (j == 2)
	{
		return cond[2 * i] * cond[2 * i];
	}
	else if (j == 3)
	{
		return cond[2 * i + 1];
	}
	else
	{
		return cond[2 * i + 1] * cond[2 * i + 1];
	}
}

/*---------------------------------------------------------------
 *	fNN_3_3()
 *--------------------------------------------------------------*/
double
fNN_3_3( unsigned int i, double *cond, double *a )
{
	double h1, h2, h3;
	double arg1, arg2, arg3;

	/* 1st hidden neuron */
	arg1 = a[0] + a[1] * cond[3 * i]
		+ a[2] * cond[3 * i + 1]
		+ a[3] * cond[3 * i + 2];
	h1 = tanh( arg1 );

	/* 2nd hidden neuron */
	arg2 = a[4] + a[5] * cond[3 * i]
		+ a[6] * cond[3 * i + 1]
		+ a[7] * cond[3 * i + 2];
	h2 = tanh( arg2 );

	/* 3rd hidden neuron */
	arg3 = a[8] + a[9] * cond[3 * i]
		+ a[10] * cond[3 * i + 1]
		+ a[11] * cond[3 * i + 2];
	h3 = tanh( arg3 );

	/* output neuron */
	return a[12] * h1 + a[13] * h2 + a[14] * h3;
}

/*---------------------------------------------------------------
 *	fNN_3_2()
 * f(x|a) =
 *--------------------------------------------------------------*/
double
fNN_3_2( unsigned int i, double *cond, double *a )
{
	double h1, h2;
	double arg1, arg2;

	/* 1st hidden neuron */
	arg1 = a[0] + a[1] * cond[3 * i]
		+ a[2] * cond[3 * i + 1]
		+ a[3] * cond[3 * i + 2];
	if (arg1 < 0) arg1 = 0;
	h1 = 2. / (1. + exp( -arg1 )) - 1;

	/* 2nd hidden neuron */
	arg2 = a[4] + a[5] * cond[3 * i]
		+ a[6] * cond[3 * i + 1]
		+ a[7] * cond[3 * i + 2];
	if (arg2 < 0) arg2 = 0;
	h2 = 2. / (1. + exp( -arg2 )) - 1;

	/* output neuron */
	return a[8] * h1 + a[9] * h2 + a[10];
}

/*---------------------------------------------------------------
 *	fNN_2_2()
 * f(x|a) =
 *--------------------------------------------------------------*/
double
fNN_2_2( unsigned int i, double *cond, double *a )
{
	double h1, h2;
	double arg1, arg2;

	/* 1st hidden neuron */
	arg1 = a[0] + a[1] * cond[2 * i]
		+ a[2] * cond[2 * i + 1];
	if (arg1 < -100)
		h1 = -1;
	else
		h1 = 2. / (1. + exp( -arg1 )) - 1;

	/* 2nd hidden neuron */
	arg2 = a[3] + a[4] * cond[2 * i]
		+ a[5] * cond[2 * i + 1];
	if (arg2 < -100)
		h2 = -1;
	else
		h2 = 2. / (1. + exp( -arg2 )) - 1;

	/* output neuron */
	return a[6] * h1 + a[7] * h2 + a[8];
}

/*---------------------------------------------------------------
 *	fNN_1_2()
 * f(x|a) =
 *--------------------------------------------------------------*/
double
fNN_1_2( unsigned int i, double *cond, double *a )
{
	double h1, h2;
	double arg1, arg2;

	/* 1st hidden neuron */
	arg1 = a[0] + a[1] * cond[i];
	/*h1 = 2. / (1. + exp( -arg1)) - 1;*/
	h1 = tanh( arg1 );

	/* 2nd hidden neuron */
	arg2 = a[2] + a[3] * cond[i];
	/* h2 = 2. / (1. + exp( -arg2)) - 1; */
	h2 = tanh( arg2 );

	/* output neuron */
	return a[4] * h1 + a[5] * h2 + a[6];
}

/*---------------------------------------------------------------
 *	fNN_1_3()
 * f(x|a) =
 *--------------------------------------------------------------*/
double
fNN_1_3( unsigned int i, double *cond, double *a )
{
	double h1, h2, h3;
	double arg;

	/* 1st hidden neuron */
	arg = a[0] + a[1] * cond[i];
	h1 = tanh( arg );
	/* h1 = 2. / (1. + exp( -arg)) - 1; */

	/* 2nd hidden neuron */
	arg = a[2] + a[3] * cond[i];
	h2 = tanh( arg );
	/* h2 = 2. / (1. + exp( -arg)) - 1; */

	/* 3rd hidden neuron */
	arg = a[4] + a[5] * cond[i];
	h3 = tanh( arg );
	/* h3 = 2. / (1. + exp( -arg)) - 1; */

	/* output neuron */
	return a[6] * h1 + a[7] * h2 + a[8] * h2 + a[9];
}

#ifdef USED
/*---------------------------------------------------------------
 *	fclassif2D()
 * f(x|a) = a1 + a2*x1 + a3*x2 + a4*x1^2 + a5*x2^2 + a6*x1*x2 + a7*x1^3 + a8*x2^3
 *--------------------------------------------------------------*/
double
fclassif2D( unsigned int i, double *cond, double *a )
{
	return a[0] + a[1] * cond[2 * i] + a[2] * cond[2 * i + 1] + a[3] * cond[2 * i] * cond[2 * i]
		+ a[4] * cond[2 * i + 1] * cond[2 * i + 1]
		+ a[5] * cond[2 * i] * cond[2 * i + 1]
		+ a[6] * cond[2 * i] * cond[2 * i] * cond[2 * i]
		+ a[7] * cond[2 * i + 1] * cond[2 * i + 1] * cond[2 * i + 1];
}
#endif

/*---------------------------------------------------------------
 *	fclassif2D_deriv()
 * f(x|a) = a1 + a2*x1 + a3*x2 + a4*x1^2 + a5*x2^2 + a6*x1*x2 + a7*x1^3 + a8*x2^3
 *--------------------------------------------------------------*/
double
fclassif2D_deriv( double(*funct) (unsigned int, double*, double*),
	unsigned int i, unsigned int j, unsigned int M, double* cond, double* a )
{
	(void)funct;
	(void)M;
	(void)a;
	if (j == 0) /* a1 */
		return 1.;
	else if (j == 1) /* a2 */
	{
		return cond[2 * i]; /* x1 and x2 are interleaved */
	}
	else if (j == 2) /* a3 */
	{
		return cond[2 * i + 1]; /* take x2 */
	}
	else if (j == 3) /* a4 */
	{
		return cond[2 * i] * cond[2 * i]; /* x1^2 */
	}
	else if (j == 4) /* a5 */
	{
		return cond[2 * i + 1] * cond[2 * i + 1]; /* take x2^2 */
	}
	else if (j == 5) /* a6 */
	{
		return cond[2 * i] * cond[2 * i + 1]; /* x1 * x2 */
	}
	else if (j == 6) /* a7 */
	{
		return cond[2 * i] * cond[2 * i] * cond[2 * i]; /* x1^3 */
	}
	else //if (j == 7) /* a8 */
	{
		return cond[2 * i + 1] * cond[2 * i + 1] * cond[2 * i + 1]; /* x2^3 */
	}
}

/*---------------------------------------------------------------
 *	fExample()
 * f(x|a) = (a1*a1-1.5*a1)*(x-2)*(x-3)*(x-4)*(x-5)/24 - a1*(x-1)*(x-3)*(x-4)*(x-5)/ 6
 *--------------------------------------------------------------*/
double fExample( unsigned int i, double* cond, double* a )
{
	return (a[0] * a[0] - 1.5 * a[0]) * 
		(cond[i] - 2) * (cond[i] - 3) * (cond[i] - 4) * (cond[i] - 5) / 24
		- a[0] * (cond[i] - 1) * (cond[i] - 3) * (cond[i] - 4) * (cond[i] - 5) / 6;
}

/*---------------------------------------------------------------
 *	fExample_deriv()
 * f(x|a) = (a1*a1-1.5*a1)*(x-2)*(x-3)*(x-4)*(x-5)/24 - a1*(x-1)*(x-3)*(x-4)*(x-5)/ 6
 *--------------------------------------------------------------*/
double
fExample_deriv( double(*funct)(unsigned int, double*, double*),
	unsigned int i, unsigned int j, unsigned int M, double* cond, double* a )
{
	(void)funct;
	(void)M;
	(void)j;
	//return a[0]* pow(cond[i], 4)/24 - 7* a[0]* pow( cond[i], 3 ) / 6. + 
	//	213./36. *a[0] * pow( cond[i], 2 ) - 77./6 * a[0]* x + 10* a[0] - 
	//	11./48 *pow( cond[i], 4 ) + 73./24 * pow( cond[i], 3 ) - 
	//	14.2708 *pow( cond[i], 2 ) + 27.4583 * cond[i] - 17.5;	/* derivation of a1 */
	return  (cond[i] - 5) * (cond[i] - 4) * (cond[i] - 3) *
		(a[0] * (cond[i] - 2) - 2.75 * cond[i] + 3.5) / 12.;
}

/*---------------------------------------------------------------
 *	fExample_deriv2()
 * f(x|a) = (a1*a1-1.5*a1)*(x-2)*(x-3)*(x-4)*(x-5)/24 - a1*(x-1)*(x-3)*(x-4)*(x-5)/ 6
 *--------------------------------------------------------------*/
double
fExample_deriv2( double(*funct)(unsigned int, double*, double*),
	unsigned int i, unsigned int j, unsigned int k, unsigned int M, 
	double* cond, double* a )
{
	(void)funct;
	(void)M;
	(void)j;
	(void)k;
	(void)a;
	return  (cond[i] - 5) * (cond[i] - 4) * (cond[i] - 3) *
		      (cond[i] - 2) / 12;
}

/*---------------------------------------------------------------
 *	f_deriv()
 *  numerical derivation, used for several model functions
 *--------------------------------------------------------------*/
double
f_deriv( double(*funct)(unsigned int, double*, double*),
	unsigned int i, unsigned int j, unsigned int M, double* cond, double* a )
{
	(void)M;
	double del, delTotal;
	double tmp1, tmp2, atmp[M_MAX], ajp, ajm;

	if (fabs( a[j] ) > TOL)
		del = thirdRootDBLEPS * max( 1.0, fabs( a[j] ) ); /* or 1+fabs( a[j] ) */
	else
		del = TOL;

	ajp = a[j] + del; /* plus a bit of current parameter value*/
	ajm = a[j] - del; /* minus % */
	/* actual distance considering rounding errors */
	delTotal = ajp - ajm;

	/* create modified parameter vector,
	 * copy maximum number of parameters
	 */
	/* copy all possible parameters, needed for POLYNOMIAL_REG	*/
	memcpy( atmp, a, sizeof( double ) * M_MAX );
	atmp[j] = ajp;
	/* look, what result is at modified position */
	tmp1 = funct( i, cond, atmp );

	atmp[j] = ajm;
	tmp2 = funct( i, cond, atmp );

	if (!isfinite(tmp1) || !isfinite( tmp2))
	{
		/* can happen, if exp() gets too high argument */
		errno = 0; /* reset error */
		return 0.00000; /* avoid nonsense */
	}
	/* compute gradient (tmp1 - tmp2) / (2 * del) */
	return   (tmp1 - tmp2) / delTotal;
}

/*---------------------------------------------------------------
 *	f_deriv2()
 *  numerical derivation, used for several model functions
 *  using central cross-derivative formula
 *--------------------------------------------------------------*/
double
f_deriv2( double(*funct)(unsigned int, double*, double*),
	unsigned int i, unsigned int j, unsigned int k, unsigned int M, double* cond, double* a )
{
	(void)funct;
	(void)M;
	double tmp, tmpPP, tmpPM, tmpMP, tmpMM, atmp[M_MAX];
	double akp, akm, ajp, ajm;
	double delj, delk, delTotalj, delTotalk;

	delj = fourthRootDBLEPS * max( 1.0, fabs( a[j] ) ); /* or 1+fabs( a[j] ) */
	delk = fourthRootDBLEPS * max( 1.0, fabs( a[j] ) ); /* or 1+fabs( a[j] ) */

	if ( j == k)
	{
		ajp = a[j] + delj; /* plus a bit of current parameter value*/
		ajm = a[j] - delj; /* minus % */
		/* actual distance considering rounding errors */
		delTotalj = ajp - ajm;

		tmp = funct( i, cond, a );

		memcpy( atmp, a, sizeof( double ) * M_MAX );
		atmp[j] = ajp;
		tmpPP = funct( i, cond, atmp );

		atmp[j] = ajm;
		tmpMM = funct( i, cond, atmp );

		if (!isfinite( tmp ) ||
			!isfinite( tmpPP ) ||
			!isfinite( tmpMM ))
		{
			/* can happen, if exp() gets too high argument */
			errno = 0; /* reset error */
			return 0.00000; /* avoid nonsense */
		}

		return   (tmpPP  - 2 * tmp + tmpMM) / (delj * delj);
	}
	else
	{
		ajp = a[j] + delj; /* plus a bit of current parameter value*/
		ajm = a[j] - delj; /* minus % */
		akp = a[k] + delk; /* plus a bit of current parameter value*/
		akm = a[k] - delk; /* minus % */
		/* actual distances considering rounding errors */
		delTotalj = ajp - ajm;
		delTotalk = akp - akm;

		/* create modified parameter vector,
		 * copy maximum number of parameters
		 */
		 /* copy all possible parameters, needed for POLYNOMIAL_REG	*/
		memcpy( atmp, a, sizeof( double ) * M_MAX );
		atmp[j] = ajp;
		atmp[k] = akp;
		/* look, what result is at modified position */
		tmpPP = funct( i, cond, atmp );

		atmp[k] = akm;
		tmpPM = funct( i, cond, atmp );

		memcpy( atmp, a, sizeof( double ) * M_MAX );
		atmp[j] = ajm;
		atmp[k] = akp;
		tmpMP = funct( i, cond, atmp );

		atmp[k] = akm;
		tmpMM = funct( i, cond, atmp );

		if (!isfinite( tmpPP ) ||
			!isfinite( tmpPM ) ||
			!isfinite( tmpMP ) ||
			!isfinite( tmpMM ))
		{
			/* can happen, if exp() gets too high argument */
			errno = 0; /* reset error */
			return 0.00000; /* avoid nonsense */
		}
		/* compute mixed second-order derivative */
		double plus = tmpPP + tmpMM;
		double minus = tmpPM + tmpMP;
		return   (plus - minus) / (2 * delTotalj * delTotalk);
	}
}

#ifdef OBSOLETE
/*---------------------------------------------------------------
 *	f_derivTinyExpr()
 *  numerical derivation, used for model functions given as program option
 * in combination with tinyexpr
 *--------------------------------------------------------------*/
double
f_derivTinyExpr_okay( te_expr *teExpr,	unsigned int j, double* a )
{
	double tmp1, tmp2, ajp, ajm, aTmp;
	double del, delTotal;

	/* inspect positions close to current one */
	if (fabs( a[j] ) > TOL)
		del = thirdRootDBLEPS * fmax( 1.0, fabs( a[j] ) ); /* or 1+fabs( a[j] ) */
	else
		del = TOL;

	ajp = a[j] + del; /* plus a bit of current parameter value*/
	ajm = a[j] - del; /* minus % */
	ajp = nextafter( a[j] + del, INFINITY ); /* guarantee difference */
	ajm = nextafter( a[j] - del, -INFINITY );
	delTotal = ajp - ajm;

	/* create modified parameter vector,
	 * copy maximum number of parameters
	 */
	 /* copy parameter	*/
	aTmp = a[j];
	a[j] = ajp;
	/* look, what result is at modified position */
	tmp1 = te_eval( teExpr );
	a[j] = ajm;
	tmp2 = te_eval( teExpr );
	/* copy parameter	back */
	a[j] = aTmp;

	if ( !isfinite(tmp1) || !isfinite( tmp2))
	{
		/* can happen, if exp() gets too high argument */
		errno = 0; /* reset error */
		return 0.00000; /* avoid nonsense */
	}
	/* compute gradient */
	return (tmp1  - tmp2) / (delTotal);
}
#endif

/*---------------------------------------------------------------
 *	f_derivTinyExpr()
 *  numerical derivation, used for model functions given as program option
 * in combination with tinyexpr
 *--------------------------------------------------------------*/
double
f_derivTinyExpr( te_expr* teExpr, unsigned int j, double* a )
{
	double tmpP, tmpM, aTmp;
	double aP, aM;
	double del, delHalf, d1, d2;

	/* inspect positions close to current one */
	del = thirdRootDBLEPS * fmax( 1.0, fabs( a[j] ) ); /* or 1+fabs( a[j] ) */

	delHalf = 0.5 * del;
	/* copy parameter	*/
	aTmp = a[j];

	/* ---- central difference with del ---- */

	/* plus a bit of current parameter value*/
	aP = nextafter( aTmp + del, INFINITY ); /* guarantee difference */
	a[j] = aP;
	tmpP = te_eval( teExpr );

	/* minus a bit of current parameter value*/
	aM = nextafter( aTmp - del, -INFINITY);
	a[j] = aM;
	tmpM = te_eval( teExpr );

	if (!isfinite( tmpP ) || !isfinite( tmpM ))
	{
		/* can happen, if exp() gets too high argument */
	  /* copy parameter	back */
		a[j] = aTmp;
		errno = 0; /* reset error */
		return 0.00000; /* avoid nonsense */
	}

	d1 = (tmpP - tmpM) / (aP - aM);

	/* ---- central difference with del ---- */

	/* plus a bit of current parameter value*/
	aP = nextafter( aTmp + delHalf, INFINITY ); /* guarantee difference */
	a[j] = aP;
	tmpP = te_eval( teExpr );

	/* minus a bit of current parameter value*/
	aM = nextafter( aTmp - delHalf, -INFINITY );
	a[j] = aM;
	tmpM = te_eval( teExpr );

	d2 = (tmpP - tmpM) / (aP - aM);

	if (!isfinite( tmpP ) || !isfinite( tmpM ))
	{
		/* copy parameter	back */
		a[j] = aTmp;
		errno = 0; /* reset error */
		return 0.00000; /* avoid nonsense */
	}

	/* copy parameter	back */
	a[j] = aTmp;

	/* compute gradient with Richardson extrapolation */
	return fma(4.0, d2, - d1) / 3.0;
}

/*---------------------------------------------------------------
 *	f_deriv2TinyExpr()
 *  numerical derivation, used for model functions given as program option
 * in combination with tinyexpr
 *--------------------------------------------------------------*/
double
f_deriv2TinyExpr( te_expr* teExpr,
	unsigned int j, unsigned int k, double* a )
{
	double tmp, tmpPP, tmpPM, tmpMP, tmpMM;
	double akp, akm, ajp, ajm;
	double delj, delk, delTotalj, delTotalk;

	delj = fourthRootDBLEPS * max( 1.0, fabs( a[j] ) ); /* or 1+fabs( a[j] ) */
	delk = fourthRootDBLEPS * max( 1.0, fabs( a[j] ) ); /* or 1+fabs( a[j] ) */

	if (j == k)
	{
		double aTmp;
		ajp = a[j] + delj; /* plus a bit of current parameter value*/
		ajm = a[j] - delj; /* minus % */
		/* actual distance considering rounding errors */
		delTotalj = ajp - ajm;

		//tmp = funct( i, cond, a );
		tmp = te_eval( teExpr ); /* take current parameters */

		 /* copy parameter	*/
		aTmp = a[j]; /* remember a */
		a[j] = ajp;
		//tmpPP = funct( i, cond, atmp );
		tmpPP = te_eval( teExpr );

		a[j] = ajm;
		//tmpMM = funct( i, cond, atmp );
		tmpMM = te_eval( teExpr );

		if (!isfinite( tmp ) ||
			!isfinite( tmpPP ) ||
			!isfinite( tmpMM ))
		{
			/* can happen, if exp() gets too high argument */
			errno = 0; /* reset error */
			return 0.00000; /* avoid nonsense */
		}
		/* copy parameter	back */
		a[j] = aTmp;

		return   (tmpPP - 2 * tmp + tmpMM) / (delj * delj);
	}
	else
	{
		double ajTmp, akTmp;
		ajp = a[j] + delj; /* plus a bit of current parameter value*/
		ajm = a[j] - delj; /* minus % */
		akp = a[k] + delk; /* plus a bit of current parameter value*/
		akm = a[k] - delk; /* minus % */
		/* actual distances considering rounding errors */
		delTotalj = ajp - ajm;
		delTotalk = akp - akm;

		/* create modified parameter vector,
		 * copy maximum number of parameters
		 */
		ajTmp = a[j]; /* remember a */
		a[j] = ajp;
		akTmp = a[k]; /* remember a */
		a[k] = akp;
		/* look, what result is at modified position */
		tmpPP = te_eval( teExpr ); /* take current parameters */

		a[k] = akm;
		tmpPM = te_eval( teExpr ); /* take current parameters */

		a[j] = ajm;
		a[k] = akp;
		tmpMP = te_eval( teExpr ); /* take current parameters */

		a[k] = akm;
		tmpMM = te_eval( teExpr ); /* take current parameters */

		if (!isfinite( tmpPP ) ||
			!isfinite( tmpPM ) ||
			!isfinite( tmpMP ) ||
			!isfinite( tmpMM ))
		{
			/* can happen, if exp() gets too high argument */
			errno = 0; /* reset error */
			return 0.00000; /* avoid nonsense */
		}
		a[j] = ajTmp;
		a[k] = akTmp;
		/* compute mixed second-order derivative */
		double plus = tmpPP + tmpMM;
		double minus = tmpPM + tmpMP;
		/* copy parameter	back */
		return   (plus - minus) / (2 * delTotalj * delTotalk);
	}
}

