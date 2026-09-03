/*****************************************************************
 *
 * File........:	computeNormalMatrixPlusVector.c
 * Function....:	fill Jacobian matrix
 * Author......:	Tilo Strutz
 * last changes:	04.08.2026
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
#include <float.h> /* for DBL_MAX */
#include "errmsg.h"
#include "prototypes.h" /* for  IsFiniteNumber() */
#include "tinyexpr.h"

extern double x; /* declared in fitting.c , for tinyExpr */

/*---------------------------------------------------------------
 *	computeNormalMatrix()
 *--------------------------------------------------------------*/
int computeNormalMatrix(
	unsigned int N, unsigned int M,
	double** jacob, double** normal, double* weights )
{
	char* rtn = "computeNormalMatrix";
	int err = 0, flag = 1;;

	for (unsigned int j = 0; j < M; j++)
	{
		for (unsigned int i = 0; i < M; i++)
		{
			normal[j][i] = 0.;
			for (unsigned int n = 0; n < N; n++)
			{
				//normal[j][i] += jacob[n][j] * jacob[n][i] * weights[n];
				normal[j][i] = fma( jacob[n][j], jacob[n][i] * weights[n], normal[j][i] );
			}
			/* overflow test	*/
			if (flag && !IsFiniteNumber( normal[j][i] ))
			{
				err = errmsg( ERR_IS_INFINITE, rtn, "element of normal", 0 );
				if (normal[j][i] == INFINITY) normal[j][i] = DBL_MAX;
				else   normal[j][i] = -DBL_MAX;
				flag = 0; /* show only once */
			}
		}
	}
	return err;
}

/*---------------------------------------------------------------
 *	computeJWy()
 *--------------------------------------------------------------*/
int computeJWy(
	unsigned int N, unsigned int M, double* observations,
	double** jacob, double* weights, double* tmpvec )
{
	int err = 0;

	/* K = J^(T) * W * y */
	for (unsigned int j = 0; j < M; j++)
	{
		tmpvec[j] = 0.;
		for (unsigned int n = 0; n < N; n++)
		{
			// tmpvec[j] += jacob[n][j] * observations[n] * weights[n];
			tmpvec[j] = fma( jacob[n][j], observations[n] * weights[n], tmpvec[j] );
		}
	}
	return err;
}

/*---------------------------------------------------------------
 *	computeJWr()
 *--------------------------------------------------------------*/
int computeJWr(
	double (*funct) (unsigned int, double*, double*),
	unsigned int N, unsigned int M, double* observations,
	double* conditions, double** jacob, double* weights,
	double* a, double* tmpvec, int expression_flag,
	te_expr* teExpr)
{
	int err = 0;
	double residual;
	/* r contains residuals */
	if (expression_flag)
	{
		for (unsigned int j = 0; j < M; j++)
		{
			tmpvec[j] = 0.;
			for (unsigned int i = 0; i < N; i++)
			{
				// residual = observations[i] - funct( i, conditions, a );
				x = conditions[i];
				residual = observations[i] - te_eval( teExpr );
				//tmpvec[j] += jacob[i][j] * residual * weights[i];
				tmpvec[j] = fma( jacob[i][j], residual * weights[i], tmpvec[j] );
			}
		}
	}
	else
	{
		for (unsigned int j = 0; j < M; j++)
		{
			tmpvec[j] = 0.;
			for (unsigned int i = 0; i < N; i++)
			{
				residual = observations[i] - funct( i, conditions, a );
				// tmpvec[j] += jacob[i][j] * residual * weights[i];
				tmpvec[j] = fma( jacob[i][j], residual * weights[i], tmpvec[j] );
			}
		}
	}
	return err;
}