/*****************************************************************
 *
 * File........:	getJacobian.c
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
#include <math.h>      /* for NAN and INFINITY */
#include "errmsg.h"
#include "tinyexpr.h"
#include "prototypes.h" /* for IsFiniteNumber() */
#include "functions.h" /* for f_derivTinyExpr() */

extern double x; /* declared in fitting.c, for tinyExpr */

/*---------------------------------------------------------------
*	getJacobian()
*
*--------------------------------------------------------------*/
int 
getJacobian( int expression_flag,
	double (*funct) (unsigned int, double*, double*),
	double (*funct_deriv) (double(*)(unsigned int, double*, double*),
		unsigned int, unsigned int, unsigned int, double*, double*),
	te_expr* teExpr,
	unsigned int N, unsigned int M, double* conditions, double** jacob,
	double* a, FILE* out )
{
	char* rtn = "getJacobian";
	int err = 0;

	if (expression_flag)
	{
		for (unsigned int i = 0; i < N; i++)
		{
			x = conditions[i];
			for (unsigned int j = 0; j < M; j++)
			{
				jacob[i][j] = f_derivTinyExpr( teExpr, j, a );
				if (!IsFiniteNumber( jacob[i][j] ))
				{
					fprintf( stderr,
						"\n# %s\n####  Divergence of approximation: %e", rtn, jacob[i][j] );
					fprintf( out,
						"\n# %s\n#  Divergence of approximation: %e", rtn, jacob[i][j] );
					if (errno)
					{
							char str[511];
							snprintf( str, 510, "\n### %s: ", rtn );
							perror( str );
						fprintf( stderr, "     errno = %d", errno );
						fprintf( out, "\terrno = %d", errno );
						fprintf( out, "\tsee standard output (console)\n" );
						err = errno;
						goto endfunc;
					}
					err = ERR_IS_INFINITE;
					goto endfunc;
				}
			}
		}
	}
	else
	{
		/* in-built model function */
		for (unsigned int i = 0; i < N; i++)
		{
			for (unsigned int j = 0; j < M; j++)
			{
				jacob[i][j] = funct_deriv( funct, i, j, M, conditions, a );
				if (errno)
				{
						char str[511];
						snprintf( str, 510, "\n### %s: ", rtn );
						perror( str );
					fprintf( stderr, "     errno = %d", errno );
					fprintf( out, "\terrno = %d", errno );
					fprintf( out, "\tsee standard output (console)\n" );
					err = errno;
					if (!IsFiniteNumber( jacob[i][j] ) || jacob[i][j] == NAN || jacob[i][j] == INFINITY)
					{
						fprintf( stderr,
							"\n# %s\n####  Divergence of approximation at (i,j) = (%d,%d)", rtn, i, j );
						fprintf( out,
							"\n# %s\n#  Divergence of approximation at (i,j) = (%d,%d)", rtn, i, j );
						/* retry*/
						goto endfunc;
					}
				}
			}
		}
	}
	endfunc:
	return err;
}