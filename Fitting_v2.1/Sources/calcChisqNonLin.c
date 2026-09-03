/****************************************************************
 *
 * File........:	calChisqNonLin.c
 * Function....:	calculation of chi^2
 * Author......:	Tilo Strutz
 * last changes:	07.01.2026
 *
 * LICENCE DETAILS: see software manual
 *	free academic use
 *  cite source as
 * "Strutz, T.: Data Fitting and Uncertainty. 3rd edition,
 *   Springer Nature, 2026"
 *
 ****************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h> /* for DBL_MAX */
#include <errno.h> /* for evaluation of errno */
#include <assert.h>
#include "defines.h"
#include "tinyexpr.h"

extern double x; /* declared in fitting.c, for tinyExpr */

/*---------------------------------------------------------------
 *	calChisqNonLin()
 *--------------------------------------------------------------*/
double calcChisqNonLin( double (*funct) (unsigned int, double*, double*), 
	unsigned int N, double* obs, double* cond, double* a,
	double* weights, int expression_flag, te_expr* teExpr)
{
	char* rtn = "calcChisqNonLin";
	double chisq = 0., datac, tmp;
	double c = 0.0; /* error compensation */

	assert( obs != NULL && cond != NULL && weights != NULL );
	assert( expression_flag ? teExpr != NULL : funct != NULL );

	for (unsigned int i = 0; i < N; i++)
	{
		errno = 0; /* reset */
		/* get calculated data points dependent on current
		parameters */
		if (expression_flag)
		{
			x = cond[i];
			datac = te_eval( teExpr ); /* uses global variable x */
		}
		else
		{
			datac = funct( i, cond, a );
		}
		tmp = obs[i] - datac;
		/* weighted and squared differences */
		// chisq += weights[i] * tmp * tmp;
		chisq = fma( weights[i] * tmp, tmp, chisq);

		if (errno)
		{
			char str[511];
			snprintf( str, 510, "\n### %s", rtn );
			perror( str );
			errno = 0;
			return DBL_MAX;
		}
		/* Kahan summation, for compensating rounding errors */
		//double tmpSq = tmp * tmp;
		//double y = tmpSq - c;
		//double t = chisq + y;
		//c = (t - chisq) - y;
		//chisq = t;
	}
	
	return chisq;
}