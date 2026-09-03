/*****************************************************************
 *
 * File........:	lsLinear.c
 * Function....:	least squares with alternative matrix inversion
 * Author......:	Tilo Strutz
 * last changes:	05.02.2008, 28.09.2009, 25.01.2010, 05.11.2025,
 * 05.01.2026, 05.08.2026
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
#include <float.h>
#include "errmsg.h"
#include "matrixUtils.h"
#include "defines.h"
#include "macros.h"
#include "ls.h"
#include "prototypes.h"
#include "functions.h"
#include "tinyexpr.h"
#ifndef WIN32
#include <sys/time.h>
#else
#include <time.h>
#define random rand
#endif


extern double x; /* declared in fitting.c */

/*---------------------------------------------------------------
*	lsLinear()
*
*--------------------------------------------------------------*/
int
lsLinear( double (*funct) (unsigned int, double*, double*),
	double (*funct_deriv) (double(*)(unsigned int, double*, double*),
		unsigned int, unsigned int, unsigned int, double*, double*),
	te_expr* teExpr, 
	unsigned int N, unsigned int M, double* observations, double* conditions, double** jacob,
	double* weights,
	double* a,
	int algoMode, LS_PARAMS* ls_params,
	double** normal_i, FILE* out)
{
	char* rtn = "ls";
	int err = 0;
	double** normal = NULL;			/* N = J^(T) * W * J	*/
	//double** normal_i = NULL;		/* inverse of N				*/
	double* tmpvec = NULL;			/* J^(T) * W * r			*/

	fprintf( out,
		"\n#  -- %s - start  ------------------------------", rtn );

	/* get Jacobian matrix  */
	err = getJacobian( ls_params->expression_flag, funct, funct_deriv,
		teExpr,	N, M, conditions, jacob, a,	out );

	if (ls_params->svd && ls_params->linear)
	{
		/* This function exploits the SVD for solving linear problems.
		* Inverting of the normal matrix is not required, which
		* sometimes runs into problems when inverting the normal matrix
		* of difficult model functions as, for example, polynomials of
		* high order. It returns the estimated parameters in a and the
		*  matrix normal_i
		*/
		err = solveLinear( N, M, observations, weights, jacob, normal_i, a, out );
		if (err)
		{
			fprintf( stderr,
				"\n\n### Unable to solve this linear system!" );
			fprintf( stderr, "\n Abort!\n" );
			goto endfunc;
		}
	}
	else /* a = (J^T * W * J)^(-1) * J^T * y */
	{
		/*
		* allocate memory
		*/
		/* normal matrix N = J^(T) * W * J, its inverse */
		normal = matrix( M, M );
		//normal_i = matrix( M, M );

		tmpvec = vector( M );		/* container for J^(T) * W * G	*/

		/*
		* calculate normal matrix 
		* normal = J^(T) * W * J
		*/
		err = computeNormalMatrix( N, M, jacob, normal, weights );
		/* K = J^(T) * W * y */
		err = computeJWy( N, M, observations, jacob, weights, tmpvec );

		/*
		 * inversion of normal matrix
		 * (cofactor method, LU decomposition, or SVD)
		 *
		 */
		err = matrixInversion( M, normal, normal_i, algoMode, out );

		/* final matrix multiplication to get parameter updates */
		for (unsigned int j = 0; j < M; j++)
		{
			a[j] = 0.0;
			for (unsigned int i = 0; i < M; i++)
			{
				//a[j] += normal_i[j][i] * tmpvec[i];		
				a[j] = fma( normal_i[j][i], tmpvec[i], a[j] );
			}
		}

		if (errno)
		{
				char str[511];
				snprintf( str, 510, "\n### %s: ", rtn );
				perror( str );
			fprintf( stderr, "     errno = %d", errno );
			fprintf( out, "\n# Error in computation (%d), ", errno );
			fprintf( out, "see standard output (console)\n" );
			err = errno;
			errno = 0;
			goto endfunc;
		}
	}

endfunc:
	fprintf( out,
		"\n#  -- %s - end  ------------------------------", rtn );

	free_vector( &tmpvec );
	free_matrix( &normal );

	return err;
}
