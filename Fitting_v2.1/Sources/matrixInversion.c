/*****************************************************************
 *
 * File........:	matrixInversion.c
 * Function....:	inversion of a square matrix
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
#include <math.h>/* for fabs() */
#include "errmsg.h"
#include "matrixUtils.h"
#include "defines.h"/* for TOL_S2 */
#include "prototypes.h"

/*---------------------------------------------------------------
*	matrixInversion()
*--------------------------------------------------------------*/
int 
matrixInversion( 	unsigned int M, double** normal, 
	double** normal_i, int algoMode, FILE* out )
{
	char* rtn = "matrixInversion";
	int err = 0;

	if (algoMode == 0)
	{
		double det;
		double** cofac = NULL;	/* cofactor matrix for matrix
														inversion */
		cofac = matrix( M, M );	/* cofactor matrix */

		switch (M)
		{
		case 1:
			/* y = a1 */
			/* df/da1 = 1 */
			det = normal[0][0];	/* determinant */
			if (det != 0)
			{
				cofac[0][0] = 1.;
			}
			else
			{
				err = errmsg( ERR_IS_ZERO, rtn, "determinant", 0 );
				goto endfunc;
			}
			break;

		case 2:
			//det = determinant_2x2( normal );	/* determinant */
			coFactor_2x2( normal, cofac );		/* coFactor matrix */
			det = determinant_2x2( normal );	/* determinant */
			break;

		case 3:
			//det = determinant_3x3( normal );	/* determinant */
			coFactor_3x3( normal, cofac );	/* coFactor matrix */
			det = determinant_3x3_from_cofactor( normal, cofac );	/* determinant */
			break;
		case 4:	/* need 4 Parameters */
			det = inverse_4x4( normal, cofac );
			break;
		case 5:	/* need 5 Parameters */
			det = inverse_5x5( normal, cofac );
			break;

		default:
			fprintf( stderr, "\n too much  parameters (%d) ", M );
			fprintf( stderr, "for standard matrix inversion" );
			err = errmsg( ERR_CALL, rtn,
				"check command-line parameters ", 0 );
			goto endfunc;
		}	/* switch */
		if (fabs( det ) > TOL_S2)
		{
			for (unsigned int i = 0; i < M; i++)
			{
				for (unsigned int j = 0; j < M; j++)
				{
					normal_i[i][j] = cofac[i][j] / det;
				}
			}
		}
		else
		{
			err = errmsg( ERR_IS_ZERO, rtn, "cofactor determinant", 0 );
			fprintf( stderr,
				"\n Please, consider to use option '-a 1' ! " );
			fprintf( out,
				"\n# %s\n### determinant is zero ! ", rtn );
			fprintf( out,
				"\n### Please, consider to use option '-a 1' !\n#" );
			goto endfunc;
		}
	}
	else if (algoMode == 1)	/* SVD -------------------------------*/
	{
		err = svdInversion( M, normal, normal_i, out );
	}	 /* algoMode SVD */
	else if (algoMode == 2)	/* LU decomposition */
	{
		unsigned int* indx = NULL;
		int s;
		double* column = NULL;

		indx = uivector( M );
		column = vector( M );

		/* decompose the matrix */
		err = decomp_LU( normal, M, indx, &s );
		if (err)
		{
			if (err == ERR_IS_ZERO)
			{
				fprintf( out, ERR_IS_ZERO_MSG,
					"decomp_LU", "'max_element'" );
				fprintf( out, "\n#  Parameters:  " );
			}
			free_uivector( &indx );
			free_vector( &column );
			goto endfunc;
		}
		/* find inverse by back-substitution of columns */
		for (unsigned int j = 0; j < M; j++)
		{
			for (unsigned int i = 0; i < M; i++)
				column[i] = 0.0;
			column[j] = 1.0;

			backsub_LU( normal, M, indx, column );

			for (unsigned int i = 0; i < M; i++)
				normal_i[i][j] = column[i];
		}
		free_vector( &column );
		free_uivector( &indx );
	}	 /* LU decomp */

endfunc:
	return err;
}