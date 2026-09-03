/*****************************************************************
 *
 * File........:	checkParameterUncertainties.c
 * Function....:	generates output for evaluation
 * Author......:	Tilo Strutz
 * last changes:	08.12.2025
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

//#include "macros.h"  /* for M_PI */
//#include "defines.h"

 /*---------------------------------------------------------------
	*	checkParameterUncertainties()
	*--------------------------------------------------------------*/
void checkParameterUncertainties( unsigned int M, double* a,
	double** covar, FILE* out, char* outname )
{
	/* check the uncertainty in parameters */
	int countUnreliabilities = 0;
	double val;

	fprintf( out, "\n#\n#  resulting uncertainty of parameters \n#" );
	for (unsigned int j = 0; j < M; j++)
	{
		if (covar[j][j] >= 0)
			fprintf( out, " %15.9G", sqrt( covar[j][j] ) );
		else
			fprintf( out, "     ??      " );
	}

	fprintf( out, "\n#" );

	for (unsigned int j = 0; j < M; j++)
	{
		if (covar[j][j] >= 0)
		{
			val = sqrt( covar[j][j] ) * 100. / fabs( a[j] );
			fprintf( out, "%15.5G%%", val );
			if (val > 10) countUnreliabilities++;
		}
		else
		{
			fprintf( out, "    ??      %%" );
		}
	}

	if (countUnreliabilities == 1)
	{
		fprintf( out,
			"\n# One parameter has relative high uncertainty !" );
		fprintf( stderr,
			"\n# One parameter has relative high uncertainty !" );
		fprintf( stderr, "\n# Please inspect file %s !", outname );
	}
	else if (countUnreliabilities > 1)
	{
		fprintf( out,
			"\n# %d parameters have relative high uncertainty !",
			countUnreliabilities );
		fprintf( stderr,
			"\n# %d parameters have relative high uncertainty !",
			countUnreliabilities );
		fprintf( stderr, "\n# Please inspect file %s !", outname );
	}
}