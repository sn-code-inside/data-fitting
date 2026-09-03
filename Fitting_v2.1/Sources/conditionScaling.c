/*****************************************************************
 *
 * File........:	conditionScaling.c
 * Function....:	preparation of conditions for selected model functions
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
#include <math.h> /* for fabs() */

#include "defines.h"

 /*---------------------------------------------------------------
	*	conditionScaling()
	*
	*--------------------------------------------------------------*/
double
conditionScaling( unsigned int N, unsigned int condDim, 
	double* conditions, unsigned int type, FILE* out )
{
	/* get maximum absolute value	*/
	double maxConditionValue, absoluteConditionValue;
	double scalingFactor = 1.0;

	switch (type)
	{
	case LINEAR:
	case LINEAR_2:
	case POLYNOM_2NDORD:
	case POLYNOM_3RDORD:
	case POLYNOMIAL:
	case POLYNOMIAL_REG:
	case NIST_THURBER:
	case NIST_RAT42:
	case NIST_RAT43:
	case NIST_ECKERLE4:
		maxConditionValue = fabs( conditions[0] );
		for (unsigned int i = 1; i < N * condDim; i++)
		{
			absoluteConditionValue = fabs( conditions[i] );
			if (maxConditionValue < absoluteConditionValue)
				maxConditionValue = absoluteConditionValue;
		}
		scalingFactor = 1. / maxConditionValue;

		/* do scaling	*/
		for (unsigned int i = 0; i < N * condDim; i++)
		{
			conditions[i] *= scalingFactor;
		}
		fprintf( out, "\n# -c , scaling activated" );
		fprintf( out, ", scale_fac = %f", scalingFactor );
		/* For cond_dim > 1 it would be even better to scale each
		 * condition separately. This would require cond_dim
		 * different scaling factors
		 */
		break;
		/* scaling has no positive effect for these nonlinear
		 * model functions */
	default:
		fprintf( out, "\n# -c , scaling not supported for '-m %d'", type );
		fprintf( stderr,
			"\n#### scaling not supported for '-m %d'###\n", type );
	}
	return scalingFactor;
}


/*---------------------------------------------------------------
 *	condition_unscaling()
 *--------------------------------------------------------------*/
void condition_unscaling( unsigned int N, 	double* conditions, 
	unsigned int condDim, unsigned int M, double* a, double scalingFactor,
	unsigned int type, FILE* out )
{
	/* correction of parameters	*/
	fprintf( out, "\n# undo the scaling" );
	switch (type)
	{
	case LINEAR:
	case POLYNOM_2NDORD:
	case POLYNOM_3RDORD:
	case EXPONENTIAL2:
		for (unsigned int j = 1; j < M; j++)
		{/* a[0] = a1 must be excluded as it is not a factor of x_i */
			a[j] *= scalingFactor;
		}
		break;
	case LINEAR_2:
		for (unsigned int j = 0; j < M; j++)
		{
			a[j] *= scalingFactor;
		}
		break;
	case EXPONENTIAL:
		/* a[0] = a1 and a[1] = a2 must be excluded as it is not a factor of x_i */
		a[2] *= scalingFactor;
		break;
	case POLYNOMIAL:
	case POLYNOMIAL_REG:
		for (unsigned int j = 1; j < M; j++)
		{/* a[0] = a1 must be excluded as it is not a factor of x_i */
			a[j] *= pow( scalingFactor, (double)j );
		}
		break;
		/* Filip = POLYNOMIAL         f(x|a) = a1 + a2 * x + a3 * x^2 + ... a11 * x^10 */
		/* Wampler 1..6 = POLYNOMIAL  f(x|a) = a1 + a2 * x + a3 * x^2 + ... a6 * x^5 */
		/* Longley = LINEAR           f(x|a) = a1 + a2 * x_1 + a3 * x_2 + ... a7 * x_6 */
	case NIST_THURBER: /* does not change results */
		a[1] *= scalingFactor;
		a[4] *= scalingFactor;
		a[2] *= scalingFactor * scalingFactor;
		a[5] *= scalingFactor * scalingFactor;
		a[3] *= scalingFactor * scalingFactor * scalingFactor;
		a[6] *= scalingFactor * scalingFactor * scalingFactor;
		break;
	case NIST_RAT42:
		a[2] *= scalingFactor;
		break;
	case NIST_RAT43:
		a[2] *= scalingFactor;
		break;
	case NIST_ECKERLE4:
		a[0] /= scalingFactor;
		a[1] /= scalingFactor;
		a[2] /= scalingFactor;
		break;
	default:
		fprintf( out, "\n# scaling not supported for '-m %d'", type );
		fprintf( stderr,
			"\n#### scaling not supported for '-m %d'###\n", type );
	}

	/* unscale conditions for output	*/
	for (unsigned int i = 0; i < N * condDim; i++)
	{
		conditions[i] /= scalingFactor;
	}
}