
/*****************************************************************
 *
 * File........: init_NIST.c
 * Function....: parameter initialisation for
 *                    different functions
 * Author......: Tilo Strutz
 * last changes: 28.09.2009, 20.01.2010
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
#include<float.h>
#ifndef WIN32
#include <sys/time.h>
#else
#include <time.h>
#define random rand
#endif

#include "functions.h"
#include "defines.h"
#include "matrixUtils.h"
#include "errmsg.h"
#include "prototypes.h" /* for LS_PARAMS */

 /* compare: https://www.itl.nist.gov/div898/strd/nls/nls_main.shtml */

 /*---------------------------------------------------------------
	*	init_NIST_thurber()
	* f(x|a) =(a1 + a2*x + a3*x**2 + a4*x**3) /
	*          (1 + a5*x + a6*x**2 + a7*x**3)
	*--------------------------------------------------------------*/
int
init_NIST_thurber( unsigned int N, double* obs, double* cond,
	double* a, unsigned char* a_flag, FILE* out )
{
	char* rtn = "init_NIST_MGH09";
	/* set 1 */
	//if (!a_flag[0])	a[0] = 1000;
	//if (!a_flag[1])	a[1] = 1000;
	//if (!a_flag[2])	a[2] = 400;
	//if (!a_flag[3])	a[3] = 40;
	//if (!a_flag[4])	a[4] = 0.7;
	//if (!a_flag[5])	a[5] = 0.3;
	//if (!a_flag[6])	a[6] = 0.03;

	/* change model to linear problem */
	/* y = (a1 + a2*x + a3*x**2 + a4*x**3) /
	*      (1 + a5*x + a6*x**2 + a7*x**3)
	*  y + a5*xy + a6*xxy + a7*xxxy = a1 + a2*x + a3*xx + a4*xxx
	*  y = a1 + a2*x + a3*xx + a4*xxx - a5*xy - a6*xxy - a7*xxxy
	**/

	int err = 0;
	unsigned int M = 7;
	double* aInit = vector( M ); /* initial parameters */
	double* xNew = vector( 7 * N );
	double* weights = vector( N );	/* vector for weights		*/
	double** jacob = matrix( N, M );	/* Jacobian */
	double** covar = matrix( M, M );	/* covariance matrix C	*/

	/* first derivatives */
	for (unsigned int i = 0; i < N; i++)
	{
		weights[i] = 1.0;
		jacob[i][0] = 1;
		jacob[i][1] = cond[i];
		jacob[i][2] = cond[i] * cond[i];
		jacob[i][3] = cond[i] * cond[i] * cond[i];
		jacob[i][4] = -obs[i] * cond[i];
		jacob[i][5] = -obs[i] * cond[i] * cond[i];
		jacob[i][6] = -obs[i] * cond[i] * cond[i] * cond[i];
	}

	fprintf( out, "\n# NIST_Thurber initialization using linear fit" );
	err = solveLinear( N, M, obs, weights,
		jacob, covar, aInit, out );
	/* covar remains unused here */
	if (err)
	{
		errmsg( err, rtn, "ups", 0 );
	}

	/* use parameters, if not provided as option */
	for (unsigned int j = 0; j < M; j++)
	{
		if (!a_flag[j])
		{
			a[j] = aInit[j];
		}
	}

	free_vector( &aInit );
	free_vector( &xNew );
	free_vector( &weights );
	free_matrix( &jacob );
	free_matrix( &covar );

	return 0;
}


/*---------------------------------------------------------------
 *	init_NIST_MGH09()
 *	 f(x|a) = a1 * (x**2 + a2*x) / (x*x + a3*x + a4)
 *--------------------------------------------------------------*/
int init_NIST_MGH09( unsigned int N, double* obs, double* cond,
			double* a, unsigned char* a_flag, FILE* out )
{
	char* rtn = "init_NIST_MGH09";
	/* set 1, far from solution */
	//if (!a_flag[0])	a[0] = 25;
	//if (!a_flag[1])	a[1] = 39;
	//if (!a_flag[2])	a[2] = 41.5;
	//if (!a_flag[3])	a[3] = 39;

  /* change model to linear problem */
  /* y = a1 * (x**2 + a2*x) / (x*x + a3*x + a4)
	 * y * (x*x + a3*x + a4) = a1 * x*x + a1*a2*x
	 * y*x*x + a3*y*x + a4*y = a1 * x*x + b2*x    b2 = a1a2
	 * y*x*x = a1 * x*x + b2*x - a3*y*x - a4*y
	 * y' = a1 * x_1 + b2*x_2 + b3*x_3 + b4*x_4  b3 = -a3; b4 = -a4
	 * with y' = y*x*x; x_1=x*x; x_2 = x; x_3 =x*y; x_4 = y
   */
	int err = 0;
	unsigned int M = 4;
	double* aInit = vector( M ); /* initial parameters */
	double* yNew = vector( N );
	double* weights = vector( N );	/* vector for weights		*/
	double** jacob = matrix( N, M );	/* Jacobian */
	double** covar = matrix( M, M );	/* covariance matrix C	*/

	/* first derivatives */
	for (unsigned int i = 0; i < N; i++)
	{
		weights[i] = 1.0;
		yNew[i] = obs[i] * cond[i] * cond[i];
		jacob[i][0] = cond[i] * cond[i];
		jacob[i][1] = cond[i];
		jacob[i][2] = -obs[i] * cond[i];
		jacob[i][3] = -obs[i];
	}

	fprintf( out, "\n# NIST_MGH09 initialization using linear fit" );
	err = solveLinear( N, M, yNew, weights,
		jacob, covar, aInit, out );
	/* covar remains unused here */
	if (err)
	{
		errmsg( err, rtn, "ups", 0 );
	}
	if (aInit[0] != 0.0)
		aInit[1] = aInit[1] / aInit[0]; /* b2 = a1a2 */
	else aInit[1] = 0.0;

	/* use parameters, if not provided as option */
	for (unsigned int j = 0; j < M; j++)
	{
		if (!a_flag[j])
		{
			a[j] = aInit[j];
		}
	}

	free_vector( &aInit );
	free_vector( &yNew );
	free_vector( &weights );
	free_matrix( &jacob );
	free_matrix( &covar );

	return 0;
}

/*---------------------------------------------------------------
 *	init_NIST_Rat42()
 * f(x|a) = a1 / (1 + exp(a2 - a3*x))
 *--------------------------------------------------------------*/
int init_NIST_Rat42( unsigned int N, double* obs, double* cond,
	double* a, unsigned char* a_flag, FILE* out )
{
	char* rtn = "init_NIST_Rat42";
	/* y = a1 / (1 + exp(a2 - a3 * x))
	 * for x -> infinity, y -> a1
	 * a1 = obs[N-1];
	 * 1/y = 1/a1 * (1 + exp(a2 - a3 * x))
	 * a1/y = 1 + exp(a2 - a3 * x)	
	 * a1/y - 1 = exp(a2 - a3 * x)
	 * ln(a1/y -1) = a2 - a3 * x
	 */
	int err = 0, idx;
	unsigned int M = 3;
	double a1;
	double* aInit = vector( M ); /* initial parameters */
	double* yNew = vector( N );
	double* weights = vector( N );	/* vector for weights		*/
	double** jacob = matrix( N, M );	/* Jacobian */
	double** covar = matrix( M, M );	/* covariance matrix C	*/

	a1 = obs[N-1] * 1.1; /* assume sorted conditions, add 10 %  */

	/* first derivatives */
	idx = 0;
	for (unsigned int i = 0; i < N; i++)
	{
		weights[idx] = 1.0;
		if (obs[i] > 0)
		{ /* since aInit[0] = obs[N-1] * 1.1, argument is > 0.0 */
			yNew[idx] = log( a1 / obs[i] - 1.0 );
		}
		jacob[idx][0] = 1;
		jacob[idx][1] = -cond[i];
		idx++;
	}

	N = idx; /* ignore data with y <= 0 */
	fprintf( out, "\n# NIST_Rat42 initialization using linear fit" );
	err = solveLinear( N, 2/*linear M*/, yNew, weights,
		jacob, covar, aInit, out );
	/* covar remains unused here */
	if (err)
	{
		errmsg( err, rtn, "ups", 0 );
	}

	/* use parameters, if not provided as option */
	if (!a_flag[0])
	{/* attention linear fit has only M=2 parameters*/
		a[0] = a1;
	}
	if (!a_flag[1])
	{/* attention linear fit has only M=2 parameters*/
		a[1] = aInit[0];
	}
	if (!a_flag[2])
	{/* attention linear fit has only M=2 parameters*/
		a[2] = aInit[1];
	}

	free_vector( &aInit );
	free_vector( &yNew );
	free_vector( &weights );
	free_matrix( &jacob );
	free_matrix( &covar );

	return 0;
}

/*---------------------------------------------------------------
	*	init_NIST_Rat42()
	* f(x|a) = a1 / (1 + exp(a2 - a3*x))
	*--------------------------------------------------------------*/
int init_NIST_Rat42_( unsigned int N, double* obs, double* cond,
	double* a, unsigned char* a_flag, FILE * out )
{
	(void) out;
	/* set 1 */
	//if (!a_flag[0])	a[0] = 100;
	//if (!a_flag[1])	a[1] = 10;
	//if (!a_flag[2])	a[2] = 0.1;
	double bestDiff;
	unsigned int idxCloseToMiddle;
	double obsMin, obsMax;
	double z1, zm;

	idxCloseToMiddle = 0;
	bestDiff = fabs( cond[N - 1] - cond[idxCloseToMiddle] );
	obsMin = obsMax = obs[0];
	for (unsigned int i = 1; i < N; i++)
	{
		/* search position closest to middle of all conditions */
		double diff; /* compare distances to first(minimum) and last(maximum) condition */
		diff = fabs( (cond[i] - cond[0]) - (cond[N - 1] - cond[i]) );
		/* is more robust point than (N>>1) when observations are not
		 * evenly distributed along condition axis */
		if (bestDiff > diff)
		{ /* better point found */
			idxCloseToMiddle = i;
			bestDiff = diff;
		}
		if (obsMin > obs[i])
		{
			obsMin = obs[i];
			//obsMinIdx = i;
		}
		if (obsMax < obs[i])
		{
			obsMax = obs[i];
			//obsMaxIdx = i;
		}
	}

	if (!a_flag[0])
	{ /* upper asymptote */
		a[0] = obsMax * 1.01;
		/* curve could alse be decreasing, then obsMin should be taken
		see init_NIST_MGH10() for proper handling */
	}

	/* y = a1 / [1 + exp(a2 - a3*x)] */
	/* 1+ exp(a2 - a3*x) = a1 / y */
	/* exp(a2 - a3*x) = a1 / y - 1 */
	/* (a2 - a3*x) = log(a1 / y - 1);   with z = log(a1 / y - 1) */
	/* z = a2 - a3*;   */
	/* a2 = z + a3*x  */
	/*  z1 + a3*x1 = zm + a3 * xm  */
	/*  a3*(x1-xm) = zm - z1 */
	/*  a3 = (zm - z1) / (x1-xm) */

	z1 = log( a[0] / obs[0] - 1 );
	zm = log( a[0] / obs[idxCloseToMiddle] - 1 );
	if (!a_flag[2])
	{
		a[2] = (zm - z1) / (cond[0] - cond[idxCloseToMiddle]);
	}
	if (!a_flag[1])
	{
		a[1] = zm + a[2] * cond[idxCloseToMiddle];
	}
	return 0;
}

/*---------------------------------------------------------------
 *	init_NIST_Rat43()
 * f(x|a) = a1 / [1 + exp(a2 - a3*x)]^(1/a4)
 *
 * Understanding the meaning of each parameter:
 * Parameter	Meaning 					How to detect it in the data
 * --------------------------------------------------------
 * a1 Upper asymptote (maximum of the curve) Maximum value of y
 * a2 Location parameter   			Related to the inflection point
 * a3 Growth rate / steepness   Slope at the inflection
 * a4 Shape parameter  					Controls left–right asymmetry
 *--------------------------------------------------------------*/
int init_NIST_Rat43( unsigned int N, double* obs, double* cond,
	double* a, unsigned char* a_flag, FILE* out )
{
	(void)out;
	/* set 1 */
	//if (!a_flag[0])	a[0] = 100;
	//if (!a_flag[1])	a[1] = 10;
	//if (!a_flag[2])	a[2] = 1;
	//if (!a_flag[3])	a[3] = 1;

	double bestDiff;
	unsigned int idxCloseToMiddle;
	double obsMin, obsMax;
	double z1, zm;

	idxCloseToMiddle = 0;
	bestDiff = fabs( cond[N - 1] - cond[idxCloseToMiddle] );
	obsMin = obsMax = obs[0];
	for (unsigned int i = 1; i < N; i++)
	{
		/* search position closest to middle of all conditions */
		double diff; /* compare distances to first(minimum) and last(maximum) condition */
		diff = fabs( (cond[i] - cond[0]) - (cond[N - 1] - cond[i]) );
		/* is more robust point than (N>>1) when observations are not
		 * evenly distributed along condition axis */
		if (bestDiff > diff)
		{ /* better point found */
			idxCloseToMiddle = i;
			bestDiff = diff;
		}
		if (obsMin > obs[i])
		{
			obsMin = obs[i];
			//obsMinIdx = i;
		}
		if (obsMax < obs[i])
		{
			obsMax = obs[i];
			//obsMaxIdx = i;
		}
	}

	if (!a_flag[0])
	{ /* upper asymptote */
		a[0] = obsMax *1.01;
		/* curve could alse be decreasing, then obsMin should be taken
		see init_NIST_MGH10() for proper handling */
	}

	/* temporarily fix a4​=1 to reduce the model to the standard logistic function*/
	/* y = a1 / [1 + exp(a2 - a3*x)] */
	/* 1+ exp(a2 - a3*x) = a1 / y */
	/* exp(a2 - a3*x) = a1 / y - 1 */
	/* (a2 - a3*x) = log(a1 / y - 1);   with z = log(a1 / y - 1) */
	/* z = a2 - a3*;   */
	/* a2 = z + a3*x  */
	/*  z1 + a3*x1 = zm + a3 * xm  */
	/*  a3*(x1-xm) = zm - z1 */
	/*  a3 = (zm - z1) / (x1-xm) */

	double arg;
	arg = a[0] / (max( TOL, obs[0] )) - 1;
	z1 = log( arg );
	//if (arg > 0.0) z1 = log( arg );
	//else  
	//	z1 = -1e14;

	arg = a[0] / (max( TOL, obs[idxCloseToMiddle] )) - 1;
	if (arg > 0.0) 	zm = log( arg );
	else 
		zm = -1e14;

	if (!a_flag[2])
	{
		if (fabs( cond[0] - cond[idxCloseToMiddle] ) > TOL)
			a[2] = (zm - z1) / (cond[0] - cond[idxCloseToMiddle]);
		else
			a[2] = 100; /* dummy value */
	}
	if (!a_flag[1])
	{
		a[1] = zm + a[2] * cond[idxCloseToMiddle];
	}
	if (!a_flag[3])
	{
		a[3] = 1; /* typical range: 0.5 .. 5 */
	}

	return 0;
}

/*---------------------------------------------------------------
 *	init_NIST_MGH10()
 *	 f(x|a) = a1 * exp( a2 / (x+a3))
 * see: https://www.itl.nist.gov/div898/strd/nls/data/LINKS/s-mgh10.shtml
 *--------------------------------------------------------------*/
int init_NIST_MGH10( unsigned int N, double* obs,
	double* cond, double* a, unsigned char* a_flag, FILE* out )
{
	/* set 1 */
	//if (!a_flag[0])	a[0] = 2;
	//if (!a_flag[1])	a[1] = 400000;
	//if (!a_flag[2])	a[2] = 25000;

	double sign; /* dummy defaut */
	double bestDiff;
	unsigned int idxCloseToMiddle;
	int variant = 0; /* curvature of graph 1..4  = A..D */
	double obsMin, obsMax;

	/* automatic */
	/* y = a1 * exp( a2 / (x+a3)) */
	/* log(y) = log(a1) + a2 / (x+a3) */
	/*        = constant + hyperbel */
	/* detect curvature */
	idxCloseToMiddle = 0;
	bestDiff = fabs( cond[N - 1] - cond[idxCloseToMiddle] );
	obsMin = obsMax = obs[0];
	for (unsigned int i = 1; i < N; i++)
	{
		/* search position closest to middle of all conditions */
		double diff; /* compare distances to first(minimum) and last(maximum) condition */
		diff = fabs( (cond[i] - cond[0]) - (cond[N - 1] - cond[i]) );
		/* is more robust point than (N>>1) when observations are not
		 * evenly distributed along condition axis */
		if (bestDiff > diff)
		{ /* better point found */
			idxCloseToMiddle = i;
			bestDiff = diff;
		}
		if (obsMin > obs[i])
		{
			obsMin = obs[i];
			//obsMinIdx = i;
		}
		if (obsMax < obs[i])
		{
			obsMax = obs[i];
			//obsMaxIdx = i;
		}
	}

	/* estimating offset a1 */
	if (obs[N - 1] < obs[0]) /* decreasing function */
	{
		fprintf( out, "\n# decreasing function" );
		sign = -1;
		if (fabs( obs[0] - obs[idxCloseToMiddle] ) >
			fabs( obs[N - 1] - obs[idxCloseToMiddle] ))
		{
			variant = 1; /* A: °.. */
			//fprintf( out, ", negative a2, positive scaling factor a1" );
		}
		else
		{
			variant = 3; /* C: °°. */
			//fprintf( out, ", positive a2, negative scaling factor a1" );
		}
	}
	else /* increasing function */
	{
		fprintf( out, "\n# increasing function" );
		sign = +1;
		if (fabs( obs[0] - obs[idxCloseToMiddle] ) >
			fabs( obs[N - 1] - obs[idxCloseToMiddle] ))
		{
			variant = 4; /* D: .°° */
			//fprintf( out, ", negative a2, negative scaling factor a1" );
		}
		else
		{
			variant = 2; /* B: ..° */
			//fprintf( out, ", positive a2, positive scaling factor a1" );
		}
	}

	if (!a_flag[0])
	{ /* assume that curve approaches asymptotically a1 */
		if (variant == 1 || variant == 2) a[0] = obsMin;
		else a[0] = obsMax;
	}

	double diff0, diffN, diffC;
	/* log(y) = log(a1) + a2 / (x+a3) */
	/* (log(y) - log(a1)) * ( x+a3) =  a2 */
	/*  (x+a3) =  a2 /  (log(y) - log(a1))  */
	/*  a3 =  a2 /  (log(y) - log(a1)) - x  = m * z + n */
	/* two equations */
	/*  a3 =  a2 /  (log(y1) - log(a1)) - x1   */
	/*  a3 =  a2 /  (log(yN) - log(a1)) - xN   */
	/*  a2 /  (log(y1) - log(a1)) - x1  =  a2 /  (log(yN) - log(a1)) - xN   */
	/*  a2 * (1./(log(y1) - log(a1)) - 1./(log(yN) - log(a1)) ) = x1 - xN   */
	/*  a2  = (x1 - xN) / * (  (1./(log(y1) - log(a1)) - 1./(log(yN) - log(a1))  )   */
	if (obs[0] - a[0] < TOL)
	{  /* a[0] = obs[0]  ==> diff is always zero ==> may not be used */
		diff0 = (log( obs[idxCloseToMiddle] ) - log( a[0] ));
		diffN = log( obs[N - 1] ) - log( a[0] );
		diffC = cond[idxCloseToMiddle] - cond[N - 1];
	}
	else
	{
		diff0 = log( obs[0] ) - log( a[0] );
		diffN = (log( obs[idxCloseToMiddle] ) - log( a[0] ));
		diffC = cond[0] - cond[idxCloseToMiddle];
	}
	if (!a_flag[1])
	{
		a[1] = diffC / (1. / diff0 - 1. / diffN);
	}

	if (!a_flag[2])
	{	/* y = a1 * exp( a2 / (x+a3)) */
		/* (x+a3) =  a2 /  (log(y) - log(a1)) */
		/* a3 =  a2 /  (log(y) - log(a1)) - x */
		double val1;
		if (obs[0] - a[0] < TOL)
		{
			val1 = a[1] / diff0 - cond[idxCloseToMiddle];
			//val2 = a[1] / diffN - cond[N - 1];
		}
		else
		{
			val1 = a[1] / diff0 - cond[0];
			//val2 = a[1] / diffN - cond[idxCloseToMiddle];
		}
		//a[2] = 0.5 * (val1 + val2);
		a[2] = val1; /* both are identical */
	}

	return 0;
}

/*---------------------------------------------------------------
 *	init_NIST_Eckerle4()
 * f(x|a) = a1 / a2 * exp(-0.5*((x -a3)/ a2)^2)
 *--------------------------------------------------------------*/
int
init_NIST_Eckerle4( unsigned int N, double* obs, double* cond,
	double* a, unsigned char* a_flag, FILE* out )
{
	(void)out;
	int  maxpos;
	double maxval;

	/* estimated parameter of a Gaussian bell */

	/* get mean value a3 at maximum point */
	maxval = obs[0];
	maxpos = 0;
	for (unsigned int i = 1; i < N; i++)
	{
		if (maxval < obs[i])
		{
			maxval = obs[i];
			maxpos = i;
		}
	}
	if (!a_flag[2])
	{
		a[2] = cond[maxpos];
	}

	/* get sigma a2	*/
	if (!a_flag[1])
	{
		a[1] = N >> 1; /* if maxpos+1 == N */
		for (unsigned int i = maxpos + 1; i < N; i++)
		{
			if (obs[i] < maxval / 2)
			{
				a[1] = i - maxpos;
				break;
			}
		}
	}
	/* get magnification */
	if (!a_flag[0])
	{
		a[0] = maxval * a[1];
	}

	return 0;
}

/*---------------------------------------------------------------
 *	init_NIST_Bennett5LOG()
 *	 f(x|a) = a1 * (x+a2)^(-1/a3)
 *	 log(-f(x|a) = log( -a1 * (x+a2)^(-1/a3) )
 *	 log(-f(x|a) = log(-a1) + log( (x+a2)^(-1/a3) )
 *	 log(-f(x|a)) = log(-a1) + (-1/a3) * log(x+a2)
 *	 log(-f(x|a)) = log(-a1) - log(x+a2) / a3
 *	 log(-f(x|a)) = b1 - log(x+a2) / a3
 * * https://www.itl.nist.gov/div898/strd/nls/data/bennett5.shtml
 *--------------------------------------------------------------*/
int
init_NIST_Bennett5LOG( unsigned int N, double* obs,
	double* cond, double* a, unsigned char* a_flag, FILE* out )
{
	(void)out;
	/* log(-y) = log(-a1) - log(x+a2) / a3 */
	/* log(-y) = b1 - log(x+a2) / a3 */
	/* ==> observations are already transformed <== */
	double b3;

	/* x + a2 >= TOL */
	/* a2 >= TOL -x */
	if (!a_flag[1])
	{
		/* use smallest condition x */
		/* assume non-negative shift */
		a[1] = TOL - min( 0, cond[0] );
	}

	/* log(-y) = log(-a1) - log(x+a2) / a3 */
	/* log(-y1) = b1      - log(x1 +a2) * b3 ; b1 = log(-a1); b3 = 1/a3  */
	/* log(-yN) = b1      - log(xN +a2) * b3   */
	/*  b1 =  log(-yN)   + log(xN +a2) * b3  =  log(-y1) + log(x1 +a2) * b3 */
	/*  log(xN +a2) * b3 - log(x1 +a2) * b3  =  log(-y1) - log(-yN) */
	/*  b3 = ( log(-y1) - log(-yN) ) / ( log(xN +a2) - log(x1 +a2)) */

	/* observations are already transformed to log(-y) */
	b3 = (obs[0] - obs[N - 1]) /
		(log( cond[N - 1] + a[1] ) - log( cond[0] + a[1] ));
	if (!a_flag[2])
	{
		a[2] = 1. / b3;
	}
	if (!a_flag[0])
	{
		double b1, bN;
		b1 = obs[0] +  log( cond[0] + a[1] ) * b3;
		//bN = obs[N - 1] + log( cond[N - 1] + a[1] ) * b3;
		bN = obs[N>>1] + log( cond[N >>1] + a[1] ) * b3;
		a[0] =  b1;
	}
	return 0;
}

/*---------------------------------------------------------------
 *	init_NIST_Bennett5()
 *	 f(x|a) = a1 * (x+a2)^(-1/a3)
 * * https://www.itl.nist.gov/div898/strd/nls/data/bennett5.shtml
 *--------------------------------------------------------------*/
int
init_NIST_Bennett5( unsigned int N, double* obs,
	double* cond, double* a, unsigned char* a_flag, FILE* out )
{
	(void)out;
	/* set 1 */
//if (!a_flag[0])	a[0] = -2000;
//if (!a_flag[1])	a[1] = 50;
//if (!a_flag[2])	a[2] = 0.8;

/* y = a1 * (x+a2)^(-1/a3) */
/* y = a1 / (x+a2)^(1/a3) */
/* log(y) = log(a1) - 1/a3 * log (x+a2) */
/* z = b1 + b3 * w  , mit b1 =log(a1),  b3 = -1/a3; w = log (x+a2) */
/* x+a2 > 0; smallest x > -a2 ==>  a2 >-x  a2= -min(x) + eps;  */
	double b1, b3;
	int sign = +1;
	double minObs, maxObs;

	if (!a_flag[1])
	{
		/* assume non-negative shift */
		a[1] = max( 0, -cond[0] ) + TOL;
	}

	minObs = maxObs = obs[0];
	for (unsigned int i = 1; i < N; i++)
	{
		if (maxObs < obs[i]) maxObs = obs[i];
		if (minObs > obs[i]) minObs = obs[i];
	}

	if (minObs < 0.)
	{ /* if observations are negativ, then a1 must be negative */
		sign = -1;
		if (maxObs >= 0)
		{
			fprintf( stderr, "\n## observations must be either all positiv or all negativ! ###" );
			fprintf( stderr, "\n## mismatch between data and model function! ###\n" );
		}
	}

	/* log(y1) = b1 + b3 * log(x1 +a2) ; b1 = log(a1); b3 = -1/a3  */
	/* log(yN) = b1 + b3 * log(xN +a2)  */
	/* or  log(-y) = b1 + b3 * log(xN +a2) ; b1 = log(-a1) */
	/*  b1 =  log(yN) - b3 * log(xN +a2)  =  log(ym) - b3 * log(x1 +a2) */
	/*  b3 * log(x1 +a2) - b3 * log(xN +a2)  =  log(y1) - log(yN) */
	/*  b3 = ( log(y1) - log(yN) ) / ( log(x1 +a2) - log(xN +a2)) */

	/* since sign has been determined, we calculate with absolute observations */
	b3 = (log( fabs(obs[0] )) - log( fabs( obs[N - 1]) )) /
		   (log( cond[0] + a[1] ) - log( cond[N - 1] + a[1] ));
	if (!a_flag[2])
	{
		a[2] = -1. / b3;
	}
	if (!a_flag[0])
	{
		b1 = log( fabs( obs[N - 1] ) ) - b3 * log( cond[N - 1] + a[1] );
		a[0] = sign * exp( b1 );
	}

	return 0;
}

/*---------------------------------------------------------------
 *	init_NIST_BoxBOD()
 *  f(x|a) = a1 * (1 - exp( -a2 * x) )
 *--------------------------------------------------------------*/
int init_NIST_BoxBOD( unsigned int N, double* obs, double* cond,
	double* a, unsigned char* a_flag, FILE* out )
{
	(void)out;
	int i, itmp;
	double mean, maxval;

	itmp = min( 5, N);

	/* estimation of a1 = head of function */
	if (!a_flag[0])
	{
		maxval = obs[0];
		for (i = 1; i < itmp; i++)
		{
			if (maxval < obs[i]) maxval = obs[i];
		}
		a[0] = maxval;
	}

	/* estimation of a1*a2 = gradient at head of function */
	if (!a_flag[1])
	{
		mean = 0;
		for (i = 1; i < itmp; i++)
			mean += (obs[i] - obs[i - 1]) / (cond[i] - cond[i - 1]);
		a[1] = mean / (itmp) / a[0];
		/* use moderate value */
		if (a[1] > 2.) a[1] = 2.;
	}

	return 0;
}
