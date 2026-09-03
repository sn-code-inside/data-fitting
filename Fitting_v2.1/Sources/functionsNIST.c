
/*****************************************************************
 *
 * File....: functions.c
 * Function: model functions and their derivatives
 * Author..: Tilo Strutz
 * Date....: 23.09.2007, 05.08.2026
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
#include <float.h>  /* for DBL_MAX */
#include "defines.h"
#include "functions.h"
#include "macros.h"

/*---------------------------------------------------------------
 *	fNIST_Eckerle4()
 * f(x|a) = (a1 / a2) * exp(-0.5*((x -a3)/ a2)^2) 
 *--------------------------------------------------------------*/
double
fNIST_Eckerle4( unsigned int i, double *cond, double *a)
{
	double x, e;
	if (fabs( a[1] ) > TOL)
		x = (cond[i] - a[2]) / a[1];
	else
	{
		if (a[1] > 0.0)
			x = (cond[i] - a[2]) / TOL;
		else
			x = -(cond[i] - a[2]) / TOL;
	}
	e = max( DBL_MIN, exp( -0.5 * x * x ));

	return (a[0] / a[1]) * e;
}

/*---------------------------------------------------------------
 *	fNIST_Eckerle4_deriv()
 * f(x|a) = (a1 / a2) * exp(-0.5*((x -a3)/ a2)^2) 
 *--------------------------------------------------------------*/
double
fNIST_Eckerle4_deriv( double (*funct)(unsigned int,double*,double*), 
					 unsigned int i, unsigned int j, unsigned int M, double *cond, double *a)
{
	(void)funct;
	(void)M;
	double x, e;
	if (fabs(a[1]) > TOL)
		x = (cond[i] - a[2]) / a[1];
	else
	{
		if (a[1] > 0.0)
			x = (cond[i] - a[2]) / TOL;
		else
			x = -(cond[i] - a[2]) / TOL;
	}
	if (j == 0)
	{
		/* y      = a0 / a1 * exp(-0.5*((x -a2)/ a1)^2)	*/
		e = max( DBL_MIN, exp( -0.5 * x * x ) );
		return   e / a[1];
	}
	else if (j == 1)
	{
		e = max( DBL_MIN, exp( -0.5 * x * x ) );
		return a[0] * e / (a[1] * a[1]) * (x * x - 1);
	}
	else
	{
		e = max( DBL_MIN, exp( -0.5 * x * x ) );
		return a[0] * e * (cond[i] - a[2]) /
			(a[1] * a[1] * a[1]);
	}
}

/*---------------------------------------------------------------
 *	fNIST_Rat43()
 * f(x|a) = a1 / [1 + exp(a2 - a3*x)]^(1/a4)
 *--------------------------------------------------------------*/
double
fNIST_Rat43( unsigned int i, double *cond, double *a)
{
	char* rtn = "fNIST_Rat43";
	double x, power;
	x =  exp( a[1] - a[2]*cond[i]);
	if (x == INFINITY)
	{ /* try to continue with something sane */
		errno = 0;
		x = DBL_MAX;
	}
	if (a[3] > 0.0)
		power = pow( 1 + x, 1. / a[3] );
	else
		power =  pow(1 + x, 1./TOL);

	if (power == INFINITY)
	{
		errno = 0;
		power = DBL_MAX;
	}

	return a[0] /power;
}

/*---------------------------------------------------------------
 *	fNIST_Rat43_deriv()
 * f(x|a) = a1 / [1 + exp(a2 - a3*x)]^(1/a4) 
 *--------------------------------------------------------------*/
double
fNIST_Rat43_deriv( double (*funct)(unsigned int,double*,double*), 
					 unsigned int i, unsigned int j, unsigned int M, double *cond, double *a)
{
	(void)funct;
	(void)M;
	double x;
	x =  exp( a[1] - a[2]*cond[i]);
	/*  y = a0 / [1 + exp( a1-a2*x)]^(1/a3) */
	if (j == 0)
		return 1. / pow(1 + x, 1/a[3]);
	else if (j == 1)
		return -a[0] * x * pow( (x+1), -1/a[3]-1 ) / a[3];
	else if (j == 2)
		return a[0]*cond[i] * x * pow( (x+1), -1/a[3]-1 ) / a[3];
	else
		return a[0] * log( x+1) / (pow( (x+1), 1/a[3]) *a[3]*a[3]);
}

/*---------------------------------------------------------------
 *	fNIST_Rat42()
 * f(x|a) = a1 / (1 + exp(a2 - a3*x)) 
 *--------------------------------------------------------------*/
double
fNIST_Rat42( unsigned int i, double *cond, double *a)
{
	double e;
	e = exp( a[1] - a[2] * cond[i] );
	if (e == INFINITY)
	{
		e = 1. / TOL;
		errno = 0;
	}
	return a[0] / (1+e);
}

/*---------------------------------------------------------------
 *	fNIST_Rat42_deriv()
 * f(x|a) = a1 / (1 + exp(a2 - a3*x)) 
 *--------------------------------------------------------------*/
double
fNIST_Rat42_deriv( double (*funct)(unsigned int,double*,double*), 
					 unsigned int i, unsigned int j, unsigned int M, double *cond, double *a)
{
	(void)funct;
	(void)M;
	/*  y = a0 / (1 + exp( a1 - a2*x)) */
	double e;
	e = exp( a[1] - a[2] * cond[i] );
	if (e > 1e100 || e == INFINITY)
	{ /* keep value sane*/
		e = 1e100;
		errno = 0;
	}

	if (j == 0)
		return 1. / (1+e);
	else if (j == 1)
			return -a[0] * e / pow( (1+e), 2. );
	else
			return a[0] * cond[i] * e /  pow( (1+e), 2. );
}

/*---------------------------------------------------------------
 *	fNIST_thurber()
 * f(x|a) =(a1 + a2*x + a3*x**2 + a4*x**3) / 
 *          (1 + a5*x + a6*x**2 + a7*x**3) 
 *--------------------------------------------------------------*/
double
fNIST_thurber( unsigned int i, double *cond, double *a)
{
	double x2, x3;
	x2 = cond[i]*cond[i];
	x3 = x2*cond[i];
	return (a[0] + a[1]*cond[i] + a[2]*x2 + a[3]*x3) / 
		        (1 + a[4]*cond[i] + a[5]*x2 + a[6]*x3);
}

/*---------------------------------------------------------------
 *	fNIST_thurber_deriv()
 *	 f(x|a) = (a1 + a2*x + a3*x**2 + a4*x**3) / 
 *             (1 + a5*x + a6*x**2 + a7*x**3)
 *--------------------------------------------------------------*/
double
fNIST_thurber_deriv( double (*funct)(unsigned int,double*,double*), 
					 unsigned int i, unsigned int j, unsigned int M, double *cond, double *a)
{
	(void)funct;
	(void)M;
	double x2, x3, div;
	x2 = cond[i]*cond[i];
	x3 = x2*cond[i];
	div = (1 + a[4] * cond[i] + a[5] * x2 + a[6] * x3);
	/* divisor may not be arbitrary small */
	if (div > 0.0 && div < TOL) 
		div = TOL;
	if (div < 0.0 && div > -TOL) 
		div = -TOL;
	switch (j)
	{
	case 0:
		/* y      = (a0 + a1*x + a2*x**2 + a3*x**3) /
								 (1 + a4*x + a5*x**2 + a6*x**3)	*/
		return 1. / div;
	case 1:
		return cond[i] / div;
	case 2:
		return      x2 / div;
	case 3:
		return      x3 / div;
	case 4:
		return -cond[i] * (a[0] + a[1] * cond[i] + a[2] * x2 + a[3] * x3) /
		(div * div); // pow( div, 2. );
	case 5:
		return      -x2 * (a[0] + a[1] * cond[i] + a[2] * x2 + a[3] * x3) /
		(div * div); // pow( div, 2. );
	case 6:
		return      -x3 * (a[0] + a[1] * cond[i] + a[2] * x2 + a[3] * x3) /
		(div * div); // pow( div, 2. );
	}
	return 0.;
}

/*---------------------------------------------------------------
 *	fNIST_thurber_deriv2()
 *	 f(x|a) = (a1 + a2*x + a3*x**2 + a4*x**3) /
 *             (1 + a5*x + a6*x**2 + a7*x**3)
 * see also http://www.mathetools.de/differenzieren/
 * https://www.wolframalpha.com/calculators/derivative-calculator/
 *--------------------------------------------------------------*/
double
fNIST_thurber_deriv2( double(*funct)(unsigned int, double*, double*),
	unsigned int i, unsigned int j, unsigned int k, unsigned int M,
	double* cond, double* a )
{
	(void)funct;
	(void)M;
	double x2, x3, div;
	x2 = cond[i] * cond[i];
	x3 = x2 * cond[i];
	div = (1 + a[4] * cond[i] + a[5] * x2 + a[6] * x3);
	if (div > 0.0 && div < TOL)
		div = TOL;
	if (div < 0.0 && div > -TOL)
		div = -TOL;

	if (j > k)
	{ /* start with smaller index */
		int tmp = j;
		j = k; k = tmp;
	}
	switch (j)
	{
		/* y      = (a0 + a1*x + a2*x**2 + a3*x**3) /
								 (1 + a4*x + a5*x**2 + a6*x**3)	*/
	case 0:
		/* y' = 1/ div = 1 / (1 + a4*x + a5*x**2 + a6*x**3) */
		switch (k)
		{
		case 0: 			return 0.;
		case 1: 			return 0.;
		case 2: 			return 0.;
		case 3: 			return 0.;
		case 4: 			return -cond[i] / (div * div);
		case 5: 			return -x2 / (div * div);
		case 6: 			return -x3 / (div * div);
		}
	case 1:
		/* y' = cond[i]/ div = cond[i] / (1 + a4*x + a5*x**2 + a6*x**3) */
		switch (k)
		{
		case 1: 			return 0.;
		case 2: 			return 0.;
		case 3: 			return 0.;
		case 4: 			return -x2 / (div * div);
		case 5: 			return -x3 / (div * div);
		case 6: 			return -x2*x2 / (div * div);
		}
	case 2:
		/* y' = x2/ div = x2 / (1 + a4*x + a5*x**2 + a6*x**3) */
		switch (k)
		{
		case 2: 			return 0.;
		case 3: 			return 0.;
		case 4: 			return -x3 / (div * div);
		case 5: 			return -x2 * x2 / (div * div);
		case 6: 			return -x3 * x2 / (div * div);
		}
	case 3:
		/* y' = x3/ div = x3 / (1 + a4*x + a5*x**2 + a6*x**3) */
		switch (k)
		{
		case 3: 			return 0.;
		case 4: 			return -x2 * x2 / (div * div);
		case 5: 			return -x3 * x2 / (div * div);
		case 6: 			return -x3 * x3 / (div * div);
		}
	case 4:
		/* y' = -x*(a0 + a1*x + a2*x**2 + a3*x**3)/ (1 + a4*x + a5*x**2 + a6*x**3)^2 */
		switch (k)
		{
		case 4: 			return 2 * x2 * (a[0] + a[1] * cond[i] + a[2] * x2 + a[3] * x3) / (div * div * div);
		case 5: 			return 2 * x3 * (a[0] + a[1] * cond[i] + a[2] * x2 + a[3] * x3) / (div * div * div);
		case 6: 			return 2 * x2*x2 * (a[0] + a[1] * cond[i] + a[2] * x2 + a[3] * x3) / (div * div * div);
		}
	case 5:
		switch (k)
		{
		case 5: 			return 2 * x2 * x2 * (a[0] + a[1] * cond[i] + a[2] * x2 + a[3] * x3) / (div * div * div);
		case 6: 			return 2 * x3 * x2 * (a[0] + a[1] * cond[i] + a[2] * x2 + a[3] * x3) / (div * div * div);
		}
	case 6:
		switch (k)
		{
		case 6: 			return 2 * x3 * x3 * (a[0] + a[1] * cond[i] + a[2] * x2 + a[3] * x3) / (div * div * div);
		}
	}
	return 0.;
}

/*---------------------------------------------------------------
 *	fNIST_MGH09()
 *	 f(x|a) =a1 * (x**2 + a2*x) / (x*x + a3*x + a4) 
 *--------------------------------------------------------------*/
double
fNIST_MGH09( unsigned int i, double *cond, double *a)
{
	double x2;  
	x2 = cond[i]*cond[i];
	return a[0] * (x2 + a[1]*cond[i]) / 
		            (x2 + a[2]*cond[i] + a[3]);
}

/*---------------------------------------------------------------
 *	fNIST_MGH09_deriv()
 *	 f(x|a) =a1 * (x**2 + a2*x) / (x*x + a3*x + a4)
 *--------------------------------------------------------------*/
double
fNIST_MGH09_deriv( double (*funct)(unsigned int, double*, double*),
	unsigned int i, unsigned int j, unsigned int M, double* cond, double* a )
{
	(void)funct;
	(void)M;
	double x2;
	x2 = cond[i] * cond[i];
	switch (j)
	{
	/* y      = a0 * (x**2 + a1*x) / (x*x + a2*x + a3)	*/
	/* dy/da0 = (x**2 + a1*x) / (x*x + a2*x + a3)	*/
	case 0: 
		return (x2 + a[1] * cond[i]) / (x2 + a[2] * cond[i] + a[3]);
	case 1:
		return a[0] * cond[i] / (x2 + a[2] * cond[i] + a[3]);
	case 2:
		return -a[0] * (x2 + a[1] * cond[i]) * cond[i] /
		pow( (x2 + a[2] * cond[i] + a[3]), 2. );
	case 3:
		return -a[0] * (x2 + a[1] * cond[i]) /
		pow( (x2 + a[2] * cond[i] + a[3]), 2. );
	}
	return 0.;
}

/*---------------------------------------------------------------
 *	fNIST_MGH09_deriv2()
 *	 f(x|a) =a1 * (x**2 + a2*x) / (x*x + a3*x + a4)
 * https://www.wolframalpha.com/calculators/derivative-calculator/
 *--------------------------------------------------------------*/
double
fNIST_MGH09_deriv2( double (*funct)(unsigned int, double*, double*),
	unsigned int i, unsigned int j, unsigned int k, unsigned int M,
	double* cond, double* a )
{
	(void)funct;
	(void)M;
	double x2, div;
	x2 = cond[i] * cond[i];
	div = x2 + a[2] * cond[i] + a[3];
	if (j > k)
	{ /* start with smaller index */
		unsigned int tmp = j;
		j = k; k = tmp;
	}
	switch (j)
	{
		/* y      = a0 * (x**2 + a1*x) / (x*x + a2*x + a3)	*/
		/* dy/da0 = (x**2 + a1*x) / (x*x + a2*x + a3)	*/
	case 0:
		switch (k)
		{
		case 0: 			return 0.;
		case 1: 			return cond[i] / div;
		case 2: 			return x2      * (a[1]+ cond[i]) / (div * div);
		case 3: 			return cond[i] * (a[1] + cond[i]) / (div * div);
		}
	case 1:
		switch (k)
		{
		case 1: 			return 0.;
		case 2: 			return a[0] * x2 / (div * div);
		case 3: 			return a[0] * cond[i] / (div * div);
		}
	case 2:
		switch (k)
		{
		case 2: 			return 2 * a[0] * x2 * cond[i] * (a[1] + cond[i]) / (div * div * div);
		case 3: 			return 2 * a[0] * x2 * (a[1] + cond[i]) / (div * div * div);
		}
	case 3:
		switch (k)
		{
		case 3: 			return 2 * a[0] * cond[i] * (a[1] + cond[i]) / (div * div * div);
		}
	}
	return 0.;
}

/*---------------------------------------------------------------
 *	fNIST_MGH10()
 *	 f(x|a) = a1 * exp( a2 / (x+a3)) 
 *--------------------------------------------------------------*/
double
fNIST_MGH10( unsigned int i, double *cond, double *a)
{
	return a[0] * exp( a[1] / (cond[i] + a[2]));
}

/*---------------------------------------------------------------
 *	fNIST_MGH10_deriv()
 *	 f(x|a) = a1 * exp( a2 / (x+a3)) 
 *--------------------------------------------------------------*/
double
fNIST_MGH10_deriv( double (*funct)(unsigned int,double*,double*), 
					 unsigned int i, unsigned int j, unsigned int M, double *cond, double *a)
{
	(void)funct;
	(void)M;
	if (j == 0)
		/* y      = a0 * exp( a1 / (x+a2))	*/
		/* dy/da0 =      exp( a1 / (x+a2))	*/
		return exp( a[1] / (cond[i] + a[2]));
	else	if (j == 1)
		/* y      = a0 * exp( a1 / (x+a2))	*/
		/* dy/da1 = a0 * exp( a1 / (x+a2)) / (x+a2)	*/
		return a[0] * exp( a[1] / (cond[i] + a[2])) / (cond[i]+ a[2]);
	else
		/* y      = a0 * exp( a1 / (x+a2))	*/
		/* dy/da2 =-a0*a1 * exp( a1 / (x+a2)) / (x+a2)^2	*/
		return -a[0] * a[1] * exp( a[1] / (cond[i] + a[2])) /
					((cond[i]+ a[2])*(cond[i]+ a[2]));
}

/*---------------------------------------------------------------
 *	fNIST_MGH10_deriv2()
 *	 f(x|a) = a1 * exp( a2 / (x+a3)) 
 *--------------------------------------------------------------*/
double
fNIST_MGH10_deriv2( double (*funct)(unsigned int,double*,double*), 
					 unsigned int i, unsigned int j,unsigned int k, 
	unsigned int M, double *cond, double *a)
{
	(void)funct;
	(void)M;
	/* i ... number of current observation */
  /* j ... derivation with respect to a_j */ 
  /* k ... derivation with respect to a_k */ 
	if (j == 0)
	{
		if ( k == 0)
			/* y      = a0 * exp( a1 / (x+a2))	*/
			/* dy/da0 =      exp( a1 / (x+a2))	*/
			return 0;	
		else if ( k == 1)	
			/* y      = a0 * exp( a1 / (x+a2))	*/
			/* dy/da0 =      exp( a1 / (x+a2))	*/
		  /* d2y/da0da1 =  exp( a1 / (x+a2)) / (x+a2)	*/
			return exp( a[1] / (cond[i] + a[2])) /
						             (cond[i] + a[2]);	
		else							
			/* y      = a0 * exp( a1 / (x+a2))	*/
			/* dy/da0 =      exp( a1 / (x+a2))	*/
		  /* d2y/da0da2 = -a1 * exp( a1 / (x+a2)) / (x+a2)^2	*/
			return  -a[1] * exp( a[1] / (cond[i] + a[2])) / 
						    ((cond[i] + a[2])*(cond[i] + a[2]));
	}
	else	if (j == 1)
	{
		if ( k == 0)
			/* y      = a0 * exp( a1 / (x+a2))	*/
			/* dy/da1 = a0 * exp( a1 / (x+a2)) / (x+a2)	*/
			/* d2y/da1da0 = exp( a1 / (x+a2)) / (x+a2)	*/
			return exp( a[1] / (cond[i] + a[2])) / 
							           (cond[i] + a[2]);
		else if ( k == 1)
			/* y      = a0 * exp( a1 / (x+a2))	*/
		  /* dy/da1 = a0 * exp( a1 / (x+a2)) / (x+a2)	*/
		  /* d2y/d2a1 = a0 * exp( a1 / (x+a2)) / (x+a2)^2	*/
			return a[0] * exp( a[1] / (cond[i] + a[2])) / 
					    ((cond[i] + a[2])*(cond[i] + a[2]));
		else
			/* y      = a0 * exp( a1 / (x+a2))	*/
		  /* dy/da1 = a0 * exp( a1 / (x+a2)) / (x+a2)	*/
		  /* d2y/da1da2 =-a0    * exp( a1 / (x+a2)) / (x+a2)^2 -
			                a0*a1 * exp( a1 / (x+a2)) / (x+a2)^3  */
		  /* d2y/da1da2 =-a0 * exp( a1 / (x+a2)) / (x+a2)^2 *
			                (1 + a1 / (x+a2))  */
			return -a[0] * exp( a[1] / (cond[i] + a[2])) / 
					     ((cond[i] + a[2])*(cond[i] + a[2])) * 
					(1 + a[1]/(cond[i]+ a[2]));
	}
	else
	{
		if ( k == 0)
			/* y      = a0 * exp( a1 / (x+a2))	*/
			/* dy/da2 =-a0*a1 * exp( a1 / (x+a2)) / (x+a2)^2	*/
		  /* d2y/da0da2 = -a1 * exp( a1 / (x+a2)) / (x+a2)^2	*/
			return  -a[1] * exp( a[1] / (cond[i] + a[2])) / 
						    ((cond[i] + a[2])*(cond[i] + a[2]));
		else if ( k == 1)
			/* y      = a0 * exp( a1 / (x+a2))	*/
			/* dy/da2 =-a0*a1 * exp( a1 / (x+a2)) / (x+a2)^2	*/
		  /* d2y/da1da2 =-a0    * exp( a1 / (x+a2)) / (x+a2)^2 -
			                a0*a1 * exp( a1 / (x+a2)) / (x+a2)^3  */
		  /* d2y/da1da2 =-a0 * exp( a1 / (x+a2)) / (x+a2)^2 *
			                (1 + a1 / (x+a2))  */
			return -a[0] * exp( a[1] / (cond[i] + a[2])) / 
					     ((cond[i] + a[2])*(cond[i] + a[2])) * 
					(1 + a[1]/(cond[i]+ a[2]));
		else
			/* y      = a0 * exp( a1 / (x+a2))	*/
			/* dy/da2 =-a0*a1 * exp( a1 / (x+a2)) / (x+a2)^2	*/
		  /* d2y/da2da2 =2*a0*a1  * exp( a1 / (x+a2)) / (x+a2)^3 +
			                a0*a1^2 * exp( a1 / (x+a2)) / (x+a2)^4  */
		  /* d2y/da2da2 = a0*a1  * exp( a1 / (x+a2)) / (x+a2)^3 *
											(2 + a1 / (x+a2) */
			return a[0] * a[1] * exp( a[1] / (cond[i] + a[2])) / 
					((cond[i]+ a[2])*(cond[i]+ a[2])*(cond[i]+ a[2])) * 
					(2 + a[1] /(cond[i]+ a[2]));
	}
}

/*---------------------------------------------------------------
 *	fNIST_Bennett5()
 *	 f(x|a) = a1 * (x+a2)^(-1/a3) 
 *--------------------------------------------------------------*/
double
fNIST_Bennett5( unsigned int i, double *cond, double *a)
{
	if (fabs( a[2]) > TOL )
	{
		if (cond[i] + a[1] > TOL)
			return a[0] * pow( (cond[i] + a[1]), -1. / a[2] );
		else
			return a[0] * pow( TOL, -1. / a[2] );
	}
	else
	{
		if (cond[i] + a[1] > TOL)
			return a[0] * pow( (cond[i] + a[1]), -1. / TOL );
		else
			return a[0] * pow( TOL, -1. / TOL );
	}
}

/*---------------------------------------------------------------
 *	fNIST_Bennett5_deriv()
 *	 f(x|a) = a1 * (x+a2)^(-1/a3) 
 *--------------------------------------------------------------*/
double
fNIST_Bennett5_deriv( double (*funct)(unsigned int,double*,double*), 
					 unsigned int i, unsigned int j, unsigned int M, double *cond, double *a)
{
	(void)funct;
	(void)M;
	double a1, a2;
	a2 = a[2];
  /* avoid division by very small value */ 
	if (fabs( a2 ) < TOL) a2 = SIGN( a2 ) * TOL;

	a1 = a[1] + cond[i];
	if (a1 < TOL)
	{/* argument of pow() must be positive */
		a[1] = -cond[i] + TOL; /* parameter correction */
		a1 = TOL;
	}

	if (j == 0)
		return pow( a1, -1./a2);	/* derivation of a1 */
	else	if (j == 1)
		return a[0] * pow( a1, -1./a2 - 1.)*(-1./a2);
							/* derivation of a2 */
	else
	{
		return a[0] * pow( a1, -1. / a2 ) *
			log( a1 ) / (a2 * a2);
	}
}
/*---------------------------------------------------------------
 *	fNIST_Bennett5_deriv2()
 *	 f(x|a) = a1 * (x+a2)^(-1/a3)
 * https://www.wolframalpha.com/calculators/derivative-calculator/
 *--------------------------------------------------------------*/
double
fNIST_Bennett5_deriv2( double (*funct)(unsigned int, double*, double*),
	unsigned int i, unsigned int j, unsigned int k, unsigned int M, double* cond, double* a )
{
	(void)funct;
	(void)M;
	double a1, a2;
	a2 = a[2];
	/* avoid division by very small value */
	if (fabs( a2 ) < TOL) a2 = SIGN( a2 ) * TOL;

	a1 = a[1] + cond[i];
	if (a1 < TOL)
	{/* argument of pow() must be positive */
		a[1] = -cond[i] + TOL; /* parameter correction */
		a1 = TOL;
	}	

	if (j > k)
	{ /* start with smaller index */
		unsigned int tmp = j;
		j = k; k = tmp;
	}
	switch (j)
	{
		/* y      = a0 * (x+a1)^(-1/a2)	*/
		/* dy/da0 = (x**2 + a1*x) / (x*x + a2*x + a3)	*/
	case 0:
		switch (k)
		{
		case 0: 			return 0.;
		case 1: 			return pow(a1,-1. - 1./a2) / a2;
		case 2: 			return pow( a1, - 1. / a2 ) * log(a1) / (a2*a2);
		}
	case 1:
		switch (k)
		{
		case 1: 			return a[0] * (a2+1)*pow( a1, -2. - 1. / a2 ) / (a2 * a2);
		case 2: 			return a[0] * pow( a1,-1 -1. / a2 ) * (a2 - log( a1 )) / (a2 * a2 * a2);
		}
	case 2:
		switch (k)
		{
		case 2: 			return a[0] * pow( a1, - 1. / a2 ) *  log( a1 )*(log(a1) - 2*a2) / (a2 * a2 * a2 * a2);
		}
	}
	return 0.;
}

/*---------------------------------------------------------------
 * fNIST_Bennett5LOG()
 * log(-f(x|a)) = log(-a1) - log(x+a2) / a3
 * log(-f(x|a)) = b1 - log(x+a2) / a3
 * log(-f(x|a)) = b1 - log(x+a2) * b3
 *--------------------------------------------------------------*/
double
fNIST_Bennett5LOG( unsigned int i, double* cond, double* a )
{
	double arg = cond[i] + a[1];

	/* avoid division by zero */
	if (fabs( a[2] ) < TOL) a[2] = SIGN( a[2] ) * TOL;

	if (arg < TOL)
	{
		/* avoid non-positive argument */
		a[1] = TOL - cond[i]; 
		return a[0] - log( TOL ) * a[2];
	}
	else
		return a[0] - log( arg ) * a[2];
}

/*---------------------------------------------------------------
 * fNIST_Bennett5LOG_deriv()
 * log(-f(x|a)) = log(-a1) - log(x+a2) / a3
 * log(-f(x|a)) = b1 - log(x+a2) * b3
 *	 d/d a1: 1./ a1
 *	 d/d b1: 1.
 *	 d/d a2: -b3 / (a2+x))
 *	 d/d a3: log(x+a2) / a3^2
 *	 d/d b3: -log(x+a2) 
 *--------------------------------------------------------------*/
double
fNIST_Bennett5LOG_deriv( double (*funct)(unsigned int, double*, double*),
	unsigned int i, unsigned int j, unsigned int M, double* cond, double* a )
{
	(void)funct;
	(void)M;
	/* i ... number of current observation */
	if (j == 0)
	{
		return 1. ;	/* derivation of b1 */
	}
	else	if (j == 1)
	{
		double denom =  (a[1] + cond[i]);
		if (fabs( denom ) < TOL) denom = SIGN( denom ) * TOL;
		return -a[2] / denom;	/* derivation of a1 */
	}
	else  /* derivation of a2 */
	{
		double arg = cond[i] + a[1];
		if (arg < TOL) arg = TOL;
		return -log( arg );
	}
}

/*---------------------------------------------------------------
 *	fNIST_BoxBOD()
 *  f(x|a) = a1 * (1 - exp( -a2 * x) )
 *--------------------------------------------------------------*/
double
fNIST_BoxBOD( unsigned int i, double *cond, double *a)
{
	return ( a[0] * (1 - exp( -a[1] * cond[i])) );
}

/*---------------------------------------------------------------
 *	fNIST_BoxBOD_deriv()
 *  f(x|a) = a1 * (1 - exp( -a2 * x) )
 *--------------------------------------------------------------*/
double
fNIST_BoxBOD_deriv( double (*funct)(unsigned int, double*, double*), 
	unsigned int i, unsigned int j, unsigned int M,
										double *cond, double *a)
{
	(void)funct;
	(void)M;
	/*  y = a0 * (1 - exp( -a1 * x) ) */
	if ( a[1] < 0 ) a[1] = 0;
	if (j == 0)
	{
			return (1 - exp( -a[1] *cond[i]));
	}
	else
	{
			return ( a[0] * cond[i] * exp( -a[1]*cond[i]) );
	}
}
