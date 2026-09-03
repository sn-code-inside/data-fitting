/*****************************************************************
 *
 * File........:	ls_straightline.c
 * Function....:	least-squares solution for straight lines
 * Author......:	Tilo Strutz
 * last changes:	20.10.2007, 23.07.2026
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

/*---------------------------------------------------------------
 *	ls_straightline()
 *--------------------------------------------------------------*/
void ls_straightline( 
	unsigned int N,					/* number of entries */
	double cond[],	/* vector of conditions */
	double obs[],		/* vector of observations */
	double a[]			/* container for parameters to be estimated */
	)
{
	double Sxx, Sx, Sy, Sxy, tmp;

		Sx = Sxx = Sy = Sxy = 0.;
		for (unsigned int i = 0; i < N; i++)
		{
			Sx += cond[i];
			Sy += obs[i];
			Sxx += cond[i] * cond[i];
			Sxy += cond[i] * obs[i];
		}

		// tmp = N * Sxx - Sx * Sx;
		tmp = diffOfProducts( N, Sxx, Sx, Sx );
		if (tmp > 1.0e-13)
		{
			//a[0] = (Sxx * Sy - Sx * Sxy) / tmp;
			a[0] = diffOfProducts( Sxx, Sy, Sx, Sxy ) / tmp;
			//a[1] = (N * Sxy - Sx * Sy) / tmp;
			a[1] = diffOfProducts( N, Sxy, Sx, Sy ) / tmp;
		}
		else
		{
			a[0] = a[1] = 0;
			fprintf( stderr, "\n### ill-conditioned problem in ls_straightline() !!###\n" );
		}
}
