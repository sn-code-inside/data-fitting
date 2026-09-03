/*****************************************************************
 *
 * File........:	specialPreparations.c
 * Function....:	pre- and postprocessing of data for selected model functions
 * Author......:	Tilo Strutz
 * last changes:	03.07.2026
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
#include <math.h>  /* for fabs() */

#include "macros.h"  /* for M_PI */
#include "defines.h"

/*---------------------------------------------------------------
 *	specialPreparations()
 *--------------------------------------------------------------*/
void
specialPreparations( unsigned int N, double* conditions, 
	double * observations, double* weights, unsigned int type )
{
	/* check input data */
	if (type == COSINE_LIN)
	{
		/* conditions in degree ? */
		double maxConditionValue, absoluteConditionValue;
		maxConditionValue = fabs( conditions[0] );
		for (unsigned int i = 1; i < N; i++)
		{
			absoluteConditionValue = fabs( conditions[i] );
			if (maxConditionValue < absoluteConditionValue) maxConditionValue = absoluteConditionValue;
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
			observations[i] = conditions[2 * i] * conditions[2 * i] + conditions[2 * i + 1] * conditions[2 * i + 1];
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
				weights[i] = 0.;
				observations[i] = -9999.;
			}
			else
				observations[i] = log( observations[i] );
		}
	}
	else if (type == NIST_BENNETT5_LOG)
	{
		/* ln( -f( x | a ) ) = (-1 / a3) * ln( -a1 * (x + a2) ) */

		/* estimates for parameters */
		for (unsigned int i = 0; i < N; i++)
		{
			if (observations[i] >= 0.0)
			{
				/* now we have a problem; this observation is invalid */
				weights[i] = 0.;
				observations[i] = -9999.;
			}
			else
				observations[i] = log( -observations[i] );
		}
	}
	if (type == HYPERBEL_LIN)
	{
		/* 1/y = (a2 + x) / a1 */
		/* estimates for parameters */
		for (unsigned int i = 0; i < N; i++)
		{
			if (observations[i] > -TOL && observations[i] < TOL)
			{
				/* now we have a problem; this observation is critical */
				weights[i] = 0.;
				observations[i] = TOL;
			}
			else
				observations[i] = 1. / observations[i];
		}
	}
	return;
}


/*---------------------------------------------------------------
 *	specialPostprocessing()
 *--------------------------------------------------------------*/
void specialPostprocessing(  unsigned int M, double* a, 
	unsigned int type, FILE* out )
{
	/* correction of parameters, before determination of relative
 * uncertainty, because of phase shift
 */
	if (type == COSINE_NONLIN) /* nonlinear cosine model */
	{ /* cos argument in degree (not radians) */
		if (a[1] < 0) /* avoid negative amplitude/radius */
		{
			a[1] = -a[1];
			a[2] = a[2] + 180.; /* phase shift of 180 degrees */
		}
	}
	else if (type == COSINE) /* nonlinear cosine model */
	{
		if (a[1] < 0 ) /* avoid negative amplitude/radius */
		{
			a[1] = -a[1];
			a[3] = a[3] - M_PI; /* phase shift of 180 degrees */
			if (a[3] < 0) a[3] += 2 * M_PI;
			fprintf( out, "\n#\n# corrected Parameters\n# ");
			for (unsigned int j = 0; j < M; j++)
			{
				fprintf( out, "a%d=%16.12G, ", j+1, a[j]);
			}
			fprintf( out, "\n#");
		}
		if (a[2] < 0) /* avoid negative frequencies */
		{
			a[2] = -a[2];
			a[3] = -a[3] + 2* M_PI;
		}
		fprintf( out, "\n# estimated period duration: %.3f", (2 * M_PI) / a[2] );
		fprintf( out, "\n# estimated frequency: %.3f\n#\t\t  (%.8f)", a[2] / (2 * M_PI), a[2] / (2 * M_PI) );
	}
	else if (type == COSINE2)
	{
		if (a[3] >  2*M_PI) a[3] -= 2*M_PI;
		else if  (a[3] < 0) a[3] += 2*M_PI;
		if (a[5] >  2*M_PI) a[6] -= 2*M_PI;
		else if  (a[5] < 0) a[6] += 2*M_PI;
		fprintf( out, "\n#\n# corrected Parameters\n# ");
		for (unsigned int j = 0; j < M; j++)
		{
			fprintf( out, "a%d=%16.12G, ", j+1, a[j]);
		}
		fprintf( out, "\n#");
	}

	if (type == NIST_ECKERLE4)
	{/* solve ambiguity */
		fprintf( out, "\n# solve ambiguity: flip signs of a1 and a2" );
		if (a[1] < 0.)
		{
			a[1] = -a[1];
			a[0] = -a[0];
		}
	}
}