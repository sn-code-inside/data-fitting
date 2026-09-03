/*****************************************************************
 *
 * File........:	svd_inversion.c
 * Function....:	matrix inversion via SVD
 * Author......:	Tilo Strutz
 * last changes:	18.01.2010
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
#include "functions.h"

/*---------------------------------------------------------------
 *	svdInversion()
 *--------------------------------------------------------------*/
int svdInversion( unsigned int M, double **normal, double **normal_i, FILE *out)
{
	char *rtn = "svdInversion";
	int err = 0;
	double thresh, smax, smin;
	double **tmpmat = NULL;	/* temporary matrix */
	double *s = NULL;	/* singular values */
	double **V = NULL;	/* V matrix */

		V = matrix( M, M);	/* V matrix for SVD */
		s = vector( M);	/* singular values for SVD */
		tmpmat = matrix( M, M);	/* temporary matrix */

		/* sanity check */
		for (unsigned int j = 0; j < M; j++)
		{
			for (unsigned int k = 0; k < M; k++)
			{
				if (normal[j][k] == INFINITY)  normal[j][k] = DBL_MAX;
				if (normal[j][k] == -INFINITY) normal[j][k] = -DBL_MAX;
			}
		}
		/* normal returns U;  */
		err = singvaldec( normal, M, M, s, V);
		if (err)
		{
				// try to operate with current result
				// goto endfunc;
		}

		smax = smin = s[0];
		for (unsigned int j = 1; j < M; j++)
		{
			if (s[j] > smax)					smax = s[j];
			if (s[j] < smin)					smin = s[j];
		}
		if (smax < TOL_S2)
		{
			fprintf( stderr, 
				"\n### %s\n### singular matrix (normal[][]), smax = %e", rtn, smax);
			fprintf( out, 
				"\n### %s\n### singular matrix (normal[][]), smax = %e", rtn,smax);

			err = ERR_IS_SINGULAR;
			for (unsigned int i = 0; i < M; i++)
			{
				fprintf( out, "\n" );
				fprintf( stderr, "\n" );
				for (unsigned int j = 0; j < M; j++)
				{
					fprintf( out, "%10.6f ", normal[i][j] );
					fprintf( stderr, "%10.6f ", normal[i][j] );
				}
			}
			fprintf( out, "\n" );
			fprintf( stderr, "\n" );
			fflush( stderr );
			goto endfunc;
		}
		if (smin > TOL)
		{
			if (smax / smin > 1.e+99)
			{
				fprintf( stderr,
					"\n### %s\n### high condition number: %e / %e = %e", rtn, smax, smin, smax / smin );
				fprintf( out,
					"\n### %s\n### high condition number: %e / %e = %e", rtn, smax, smin, smax / smin );
				for (unsigned int i = 0; i < M; i++)
				{
					fprintf( stderr, "\n" );
					for (unsigned int j = 0; j < M; j++)
					{
						fprintf( stderr, "%10.6f ", normal[i][j] );
					}
				}
				fprintf( stderr, "\n" );
				fflush( stderr );
				err = ERR_HIGH_CONDITIONNUMBER;
				goto endfunc;
			}
		}

		thresh = min( DBL_EPSILON * 100* M * smax, DBL_EPSILON * M*100 ); /* for smax < 1 adapt to its value, otherwise keep any s > TOL_S */
		//thresh = DBL_EPSILON * M * smax; leads to strange behaviour for BENNETT5 and MGH10 !!

		/* invert singular values */
		for (unsigned int j = 0; j < M; j++)
		{
			/* <= in case of smax =0 */
			if (s[j] <= thresh)	s[j] = 0.0;
			else      					s[j] = 1. / s[j];
			/* Tikhonov-filtered */
			// s[j] = s[j] / (s[j] * s[j] + thresh * thresh);  failure for Eckerle start1 -D and Rat43 -Ge
		}

		/* A = U * S * V^T */
		/* invA = V * invS * U^T */
		/* a[i][j] = Sum_k V[i][k] * invS[k] * U[j][k] */
		for (unsigned int i = 0; i < M; i++)
		{
			for (unsigned int j = 0; j < M; j++)
			{
				long double sum = 0.0;
				for (unsigned int k = 0; k < M; k++)
				{
					sum = fmal( (long double)V[i][k] *
						(long double)s[k],
						(long double)normal[j][k], sum );
				}
				normal_i[i][j] = (double)sum;
			}
		}

endfunc:
	free_vector( &s);
	free_matrix( &V);
	free_matrix( &tmpmat);

	return err;
}
