/*****************************************************************
 *
 * File........:	condition_scaling.c
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
#include <math.h> /* for fabs() */

#include "defines.h"

 /*---------------------------------------------------------------
	*	condition_scaling()
	*
	*--------------------------------------------------------------*/
double
condition_scaling( unsigned int N, unsigned int condDim, 
	double* conditions, unsigned int type, FILE* out )
{
	/* get maximum absolute value	*/
	double maxConditionValue, absoluteConditionValue;
	double scalingFactor = 1.0;

	switch (type)
	{
	case LINEAR:
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
		fprintf( out, "\n# scaling activated" );
		fprintf( out, ", scale_fac = %f", scalingFactor );
		/* For cond_dim > 1 it would be even better to scale each
		 * condition separately. This would require cond_dim
		 * different scaling factors
		 */
		break;
	case POLYNOMIAL:
		/* scaling has no positive effect for these nonlinear
		 * model functions */
	case NIST_RAT42:
	case NIST_ECKERLE4:
	case NIST_MGH10:
		maxConditionValue = fabs( conditions[0] );
		for (unsigned int i = 1; i < N; i++)
		{
			absoluteConditionValue = fabs( conditions[i] );
			if (maxConditionValue < absoluteConditionValue)
				maxConditionValue = absoluteConditionValue;
		}
		scalingFactor = 1. / maxConditionValue;

		/* do scaling	*/
		for (unsigned int i = 0; i < N; i++)
		{
			conditions[i] *= scalingFactor;
		}
		fprintf( out, "\n# scaling activated" );
		fprintf( out, ", scale_fac = %f", scalingFactor );
		break;
	default:
		fprintf( out, "\n# scaling not supported for '-m %d'", type );
		fprintf( stderr,
			"\n#### scaling not supported for '-m %d'###\n", type );
	}
	return scalingFactor;
}


/*---------------------------------------------------------------
 *	condition_unscaling()
 *
 *--------------------------------------------------------------*/
void condition_unscaling( unsigned int M, 
	double* a, double scalingFactor, unsigned int type, FILE* out )
{
	switch (type)
	{
		/* scaling has no positive effect for these nonlinear
		 * model functions */
	case NIST_RAT42:
		a[2] /= scalingFactor;
		break;
	case NIST_ECKERLE4:
		a[0] *= scalingFactor;
		a[1] *= scalingFactor;
		a[2] *= scalingFactor;
		break;
	case NIST_MGH10:
		a[1] *= scalingFactor;
		a[2] *= scalingFactor;
		break;
	}
	fprintf( out, "\n# scaled\n# " );
	/* write initial parameters to output */
	for (unsigned int j = 0; j < M; j++)
	{
		fprintf( out, "a%d=%.9f, ", j + 1, a[j] );
	}
}

/*---------------------------------------------------------------
 *	condition_unscaling2()
 *
 *--------------------------------------------------------------*/
void condition_unscaling2( unsigned int N, 	double* conditions, 
	unsigned int condDim, unsigned int M, double* a, double scalingFactor,
	unsigned int type, FILE* out )
{
	/* correction of parameters	*/
	fprintf( out, "\n# undo the scaling" );
	switch (type)
	{
	case LINEAR:
		for (unsigned int j = 1; j < M; j++)
		{
			a[j] *= scalingFactor;
		}
		break;
	case POLYNOMIAL:
		for (unsigned int j = 1; j < M; j++)
		{
			a[j] *= pow( scalingFactor, (double)j );
		}
		break;
	case NIST_RAT42:
		a[2] *= scalingFactor;
		break;
	case NIST_ECKERLE4:
		a[0] /= scalingFactor;
		a[1] /= scalingFactor;
		a[2] /= scalingFactor;
		break;
	case NIST_MGH10:
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
		//if (type == EXPRESSION)
		//{
		//	x = cond[i];
		//	datac[i] = te_eval( teExpr );
		//}
		//else
		//{
		//	datac[i] = funct( i, cond, a );
		//}

	}
}