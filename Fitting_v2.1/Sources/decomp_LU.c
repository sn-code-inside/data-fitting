/*****************************************************************
 *
 * File........:	decomp_LU.c
 * Function....:	LU decomposition
 * Author......:	Tilo Strutz
 * last changes:	20.10.2007
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
#include <math.h>
#include "matrixUtils.h"
#include "errmsg.h"
#define LITTLEBIT 1.0e-20;

/*---------------------------------------------------------------
 *	decomp_LU()
 *--------------------------------------------------------------*/
int
decomp_LU( double **normal, unsigned int M, unsigned int *indx, int *s)
{
	char *rtn = "decomp_LU";
	int err = 0;
	unsigned int imax = 0;
	double maxElement, val, sum, tmp;
	double *rowScale = NULL;	/* scaling of each row */

	/* allocate vector */
	rowScale = vector( M);

	*s = 1;
	/* examine input matrix */
	for (unsigned int i = 0; i < M; i++)
	{
		maxElement = 0.0;
		for (unsigned int j = 0; j < M; j++)
		{
			if (( tmp = fabs( normal[i][j])) > maxElement)
				maxElement = tmp;
		}
		if (maxElement == 0.0)
		{
			err = errmsg( ERR_IS_ZERO, rtn, "'maxElement'", 0);
			goto endfunc;
		}
		rowScale[i] = 1.0 / maxElement;
	}
	/* loop over columns of Crout's method */
	for (unsigned int j = 0; j < M; j++)
	{
		for (unsigned int i = 0; i < j; i++)
		{
			sum = normal[i][j];
			for (unsigned int k = 0; k < i; k++)
			{
				//sum -= normal[i][k] * normal[k][j];
				/* use fma() for less rounding operations */
				sum = fma( -normal[i][k], normal[k][j], sum );
			}
			normal[i][j] = sum;
		}
		maxElement = 0.0;
		for (unsigned int i = j; i < M; i++)
		{
			sum = normal[i][j];
			for (unsigned int k = 0; k < j; k++)
			{
				//sum -= normal[i][k] * normal[k][j];
				/* use fma() for less rounding operations */
				sum = fma( -normal[i][k], normal[k][j], sum );
			}
			normal[i][j] = sum;
			/* is new pivot better than current best ? */
			if (( val = rowScale[i] * fabs( sum)) >= maxElement)
			{
				maxElement = val;
				imax = i;
			}
		}
		if (j != imax)
		{
			/* interchange of rows */
			for (unsigned int k = 0; k < M; k++)
			{
				val = normal[imax][k];
				normal[imax][k] = normal[j][k];
				normal[j][k] = val;
			}
			*s = -( *s);
			/* interschange scale factors */
			rowScale[imax] = rowScale[j];
		}
		indx[j] = imax;
		if (normal[j][j] == 0.0)
		{
			err = errmsg( ERR_IS_ZERO, rtn, "'pivot element'", 0 );
			normal[j][j] = LITTLEBIT;   /* still avoid a hard crash downstream */
		}
		if (j != (M - 1))
		{
			/* divide by the pivot element */
			val = 1.0 / (normal[j][j]);
			for (unsigned int i = j + 1; i < M; i++)
				normal[i][j] *= val;
		}
		/* next column in reduction */
	}

endfunc:
	free_vector( &rowScale);
	return err;
}

/*---------------------------------------------------------------
 *	backsub_LU()
 *--------------------------------------------------------------*/
void
backsub_LU( double **lu, unsigned int N, unsigned int *indx, double back[])
{
	int ii, idx;
	double sum;

	ii = -1;
	for (unsigned int i = 0; i < N; i++)
	{
		idx = indx[i];
		sum = back[idx];
		back[idx] = back[i];
		if (ii >= 0)
		{
			for (unsigned int j = ii; j <= i - 1; j++)
			{
				// sum -= lu[i][j] * back[j];
				/* use fma() for less rounding operations */
				sum = fma( -lu[i][j], back[j], sum );
			}
		}
		else if (sum)
			ii = i;
		back[i] = sum;
	}
	for (int i = N - 1; i >= 0; i--)
	{
		sum = back[i];
		for (unsigned int j = i + 1; j < N; j++)
		{
			//sum -= lu[i][j] * back[j];
			/* use fma() for less rounding operations */
			sum = fma( -lu[i][j], back[j], sum );
		}
		back[i] = sum / lu[i][i];
	}
}
