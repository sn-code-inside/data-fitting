/*****************************************************************
 *
 * File........:	special_preparations.c
 * Function....:	solving linear least squares via SVD
 * Author......:	Tilo Strutz
 * last changes:	25.01.2010
 *
 * LICENCE DETAILS: see software manual
 *	free academic use
 *  cite source as
 * "Strutz, T.: Data Fitting and Uncertainty. 2nd edition,
 *   Springer Vieweg, 2016"
 *
 *****************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "macros.h"  /* for M_PI */
#include "defines.h"

 /*---------------------------------------------------------------
	*	special_preparations()
	*
	*--------------------------------------------------------------*/
void
special_preparations( unsigned int N, double* conditions, 
	double * observations, double* weights, unsigned int type, FILE* out )
{

	/* check input data */
	if (type == COSINE_LIN || type == COSINE)
	{
		/* conditions in degree ? */
		double maxConditionValue, absoluteConditionValue;
		maxConditionValue = fabs( conditions[0] );
		for (unsigned int i = 1; i < N; i++)
		{
			absoluteConditionValue = fabs( conditions[i] );
			if (maxConditionValue < absoluteConditionValue)
				maxConditionValue = absoluteConditionValue;
		}
		if (maxConditionValue < 2 * M_PI)
		{
			fprintf( stderr, "\n======================================" );
			fprintf( stderr, "\n== range of degrees is very small!! ==" );
			fprintf( stderr, "\n== Mismatch with radians??          ==" );
			fprintf( stderr, "\n== Please check.\n                  ==" );
			fprintf( stderr, "\n======================================" );
		}
	}

	/* prepare input data */
	if (type == CIRCLE_LIN)
	{
		for (unsigned int i = 0; i < N; i++)
		{
			observations[i] = conditions[2 * i] * conditions[2 * i] +
				conditions[2 * i + 1] * conditions[2 * i + 1];
		}
	}


	/*
	 * pre-processing if necessary
	 */

	 /* linearisation */
	if (type == EXPONENTIAL2_LIN)
	{
		/* ln(f(x|a)) = ln(a1) + (a2 * x) */
		/* estimates for parameters */
		for (unsigned int i = 0; i < N; i++)
		{
			if (observations[i] <= 0.0)
			{
				/* now we have a problem; this observation is invalid */
				weights[i] = 0;
				observations[i] = -9999.;
			}
			else
				observations[i] = log( observations[i] );
		}
	}
	if (type == EXP_BASE_X_LOG)
	{
		/* ln(f(x|a)) = a2 + ln(a1 * x) */
		/* estimates for parameters */
		for (unsigned int i = 0; i < N; i++)
		{
			if (observations[i] <= 0.0)
			{
				/* now we have a problem; this observation is invalid */
				weights[i] = 0;
				observations[i] = -9999.;
			}
			else
				observations[i] = log( observations[i] );
		}
	}
	return;
}


/*---------------------------------------------------------------
 *	special_postprocessing()
 *
 *--------------------------------------------------------------*/
void special_postprocessing(  unsigned int M, double* a, 
	unsigned int type, FILE* out )
{
	/* correction of parameters, before determination of relative
 * uncertainty, because of phase shift
 */
	if (type == COSINE) /* nonlinear cosine model */
	{
		if (a[1] < 0) /* avoid negative amplitude/radius */
		{
			a[1] = -a[1];
			a[2] = a[2] - 180; /* phase shift of 180 degrees */
		}
		fprintf( out, "\n#\n# corrected Parameters\n# " );
		for (unsigned int j = 0; j < M; j++)
		{
			fprintf( out, "a%d=%16.12G, ", j + 1, a[j] );
		}
		fprintf( out, "\n#" );
	}
	else if (type == TRIGONOMETRIC2)
	{
		if (a[3] > 2 * M_PI) a[3] -= 2 * M_PI;
		else if (a[3] < 0) a[3] += 2 * M_PI;
		if (a[5] > 2 * M_PI) a[6] -= 2 * M_PI;
		else if (a[5] < 0) a[6] += 2 * M_PI;
		fprintf( out, "\n#\n# corrected Parameters\n# " );
		for (unsigned int j = 0; j < M; j++)
		{
			fprintf( out, "a%d=%16.12G, ", j + 1, a[j] );
		}
		fprintf( out, "\n#" );
	}
}