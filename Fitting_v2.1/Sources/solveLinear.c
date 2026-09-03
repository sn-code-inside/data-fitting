/*****************************************************************
 *
 * File........:	solveLinear.c
 * Function....:	solving linear least squares via SVD
 * Author......:	Tilo Strutz
 * last changes:	25.01.2010, 19.06.2026
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
#include "errmsg.h"
#include "matrixUtils.h"
#include "defines.h"
#include "prototypes.h"

 //#define QDEBUG 
 /*---------------------------------------------------------------
	* inverse normal matrix:  (J^T*W*J)^(-1) = V * diag(s^-2) * V^T
	*
	* Assumes s[j] already holds 1/sigma_j (i.e. already inverted,
	* as done in solveLinear for zero-filtered singular values).
	* Computes only the upper triangle and mirrors it, and factors out
	* the elementwise V*w scaling so it's done once per (i,j)
	*--------------------------------------------------------------*/
static void computeInverseNormal( unsigned int M, double** V,
	double* s, double** normal_i )
{
	double* w = vector( M );     /* w[j] = s[j]^2 = 1/sigma_j^2 */
	double** VW = matrix( M, M ); /* VW[i][j] = V[i][j] * w[j]   */

	for (unsigned int j = 0; j < M; j++)
	{
		w[j] = s[j] * s[j];
	}

	for (unsigned int i = 0; i < M; i++)
	{
		for (unsigned int j = 0; j < M; j++)
		{
			VW[i][j] = V[i][j] * w[j];
		}
	}

	for (unsigned int i = 0; i < M; i++)
	{
		for (unsigned int k = i; k < M; k++)
		{
			double sum = 0.0;
			for (unsigned int j = 0; j < M; j++)
				sum = fma( VW[i][j], V[k][j], sum );

			normal_i[i][k] = sum;
			normal_i[k][i] = sum;   /* symmetric, fill both sides */
		}
	}

	free_vector( &w );
	free_matrix( &VW );
}

/*---------------------------------------------------------------
 *	solveLinear()
 *--------------------------------------------------------------*/
int solveLinear( unsigned int N, unsigned int M, double* obs,
	double* weights, double** jacob, double** normal_i, double* a,
	FILE* out )
{
	char* rtn = "solve_lin";
	int err = 0;
	double thresh, smax, * sqrtwe = NULL;
#ifdef DEACTIVE
	double** tmpmat = NULL;	/* temporary matrix */
	double** tmpmat2 = NULL;	/* temporary matrix */
#endif
	double* s = NULL;	/* singular values */
	double** V = NULL;	/* V matrix */
	double** WJ = NULL;	/* W*J matrix */

	if (N < M)
	{
		err = errmsg( ERR_ANY, rtn, "N < M", N );
		return err;
	}
	V = matrix( M, M );	/* V matrix for SVD */
	s = vector( M );	/* singular values for SVD */
	WJ = matrix( N, M );	/* temporary matrix */
	sqrtwe = vector( N );	/* square root of weights */
#ifdef DEACTIVE
	tmpmat = matrix( M, M );	/* temporary matrix */
	tmpmat2 = matrix( M, N );	/* temporary matrix */
#endif

	/* WJ = sqrt(W)*J  */
	for (unsigned int i = 0; i < N; i++)
	{
		if (weights[i] < 0.0)
		{
			err = errmsg( ERR_ANY, rtn, "weights[i] < 0", i );
			goto endfunc;
		}
		sqrtwe[i] = sqrt( weights[i] ); /* store for later use */
		for (unsigned int j = 0; j < M; j++)
		{
			WJ[i][j] = sqrtwe[i] * jacob[i][j];
		}
	}

	/* do the SVD	*/
	err = singvaldec( WJ, N, M, s, V );
	if (err)
	{
		goto endfunc;
	}

	/* check the singular values	*/
	smax = 0.0;
	for (unsigned int j = 0; j < M; j++)
	{
		if (s[j] > smax)	smax = s[j];
	}
	if (smax < TOL_S)
	{
		fprintf( stderr,
			"\n### %s\n###     singular matrix (WJ), smax = %f", rtn, smax );
		fprintf( out,
			"\n### %s\n###     singular matrix (WJ), smax = %f", rtn, smax );

		err = 1;
		for (unsigned int i = 0; i < M; i++)
		{
			fprintf( out, "\n" );
			fprintf( stderr, "\n" );
			for (unsigned int j = 0; j < M; j++)
			{
				fprintf( out, "%10.6f ", WJ[i][j] );
				fprintf( stderr, "%10.6f ", WJ[i][j] );
			}
		}
		fprintf( out, "\n" );
		fprintf( stderr, "\n" );
		goto endfunc;
	}
	else if (smax > 1.e+31)
	{
		fprintf( stderr,
			"\n### %s\n###     degraded matrix (WJ), smax = huge", rtn );
		fprintf( out,
			"\n### %s\n###     degraded matrix (WJ), smax = huge", rtn );
		err = 1;
		goto endfunc;
	}

	//thresh = min( TOL_S * smax, TOL_S);
	thresh = min( DBL_EPSILON * 100 * M * smax, DBL_EPSILON * M * 100 ); /* for smax < 1 adapt to its value, otherwise keep any s > TOL_S */

	fprintf( out, "\n#\n#  singular values (thresh = %.14G)\n#  ",
		thresh );
	for (unsigned int j = 0; j < M; j++)
	{
		fprintf( out, "s%d=%.14G, ", j + 1, s[j] );
	}

	/* singular values inverted in place: s[j] now holds 1/sigma_j, or 0
	 * if sigma_j was below threshold (rank-deficient direction removed) */
	for (unsigned int j = 0; j < M; j++)
	{
		/* <= in case of smax =0 */
		if (s[j] <= thresh)
			s[j] = 0.0;
		else
			s[j] = 1. / s[j];
	}

	double* z = vector( M );
	double* wObs = vector( N );/* wObs[i] = sqrt(w_i) * y_i, one per observation */
	for (unsigned int i = 0; i < N; i++)
	{
		wObs[i] = sqrtwe[i] * obs[i];
	}
	/* z = diag(1/s) * U^T * sqrt(W) * y
	 *
	 * After singvaldec(), WJ has been overwritten in place and now holds
	 * U (N x M, orthonormal columns) instead of the original weighted
	 * Jacobian. The loop below forms U^T * wObs (length M): each z[j] is
	 * the projection of the weighted residual vector onto the j-th left
	 * singular vector, i.e. how much of the data lies along that
	 * principal direction of sqrt(W)*J. Multiplying by s[j] (=1/sigma_j)
	 * then rescales that projection by the inverse of how "important"
	 * (well-conditioned) that direction is - directions with sigma_j
	 * below the threshold were already zeroed out above, so they
	 * contribute nothing here (rank-deficient directions dropped).
	 *
	 * z therefore holds the least-squares solution expressed in the
	 * rotated/orthogonal basis given by V's columns, *before* rotating
	 * back into the original parameter basis.
	 */
	for (unsigned int j = 0; j < M; j++)
	{
		z[j] = 0.0;
		for (unsigned int i = 0; i < N; i++)
		{
			z[j] = fma( WJ[i][j], wObs[i], z[j] );  /* WJ now holds U */
		}
		z[j] *= s[j];   /* s[j] already 1/sigma_j from earlier inversion */
	}
	/* a = V * z		 */
	for (unsigned int j = 0; j < M; j++)
	{
		a[j] = 0.0;
		for (unsigned int k = 0; k < M; k++)
		{
			a[j] = fma( V[j][k], z[k], a[j] );
		}
	}
	free_vector( &wObs );
	free_vector( &z );

	/* normal_i = (J^T * W * J)^(-1) = V * diag(1/sigma_j^2) * V^T
	 *
	 * The inverse of the (weighted) normal matrix, i.e. the parameter
	 * covariance matrix up to a scale factor (multiply by the reduced
	 * chi-square of the fit residuals if an absolute covariance is
	 * needed and weights are only relative). Computed directly from the
	 * SVD factors already in hand (V, s) rather than re-forming J^T*W*J
	 * and inverting it - this is both cheaper and better conditioned,
	 * since it avoids squaring the condition number of sqrt(W)*J.
	 */
	computeInverseNormal( M, V, s, normal_i );
endfunc:
	free_vector( &s );
	free_vector( &sqrtwe );
	free_matrix( &V );
	free_matrix( &WJ );
#ifdef DEACTIVE
	free_matrix( &tmpmat );
	free_matrix( &tmpmat2 );
#endif

	return err;
}
