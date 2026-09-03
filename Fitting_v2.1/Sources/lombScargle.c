/****************************************************************
 *
 * File........:	lombScargle.c
 * Function....:	compute periodogram for a frequency grid
 * Author......:	Tilo Strutz
 * last changes:	01.06.2026, 
 *
 * LICENCE DETAILS: see software manual
 *	free academic use
 *  cite source as
 * "Strutz, T.: Data Fitting and Uncertainty. 3rd edition,
 *   Springer Nature, 2026"
 *
 ****************************************************************/
//#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "defines.h" /* for M_PI */

/*---------------------------------------------------------------
 *	lombScargle()
 *--------------------------------------------------------------*/
int lombScargle( unsigned int N, double* obs, double* cond,
	double obsMean, int fVecLen, double* fVec, double* p, double *bestPhi )
{
	int bestIdx;
	double omega;
	double sumCos, sumCos2;
	double sumSin, sumSin2;
	double maxP, phi;

	maxP = 0.0;
	bestIdx = 0;
	for (int idx =0; idx < fVecLen; idx++)
	{
		omega = 2 * M_PI * fVec[idx];
		/* phi calculation */
		sumCos = sumSin = 0.0;
		for (unsigned int i = 0; i < N; i++)
		{
			sumCos += cos( 2 * omega * cond[i] );
			sumSin += sin( 2 * omega * cond[i] );
		}
		phi = atan( sumSin / sumCos ) * 0.5;
		/* power calculation */
		sumCos = sumSin = 0.0;
		sumCos2 = sumSin2 = 0.0;
		for (unsigned int i = 0; i < N; i++)
		{
			double cosTmp = cos( omega * cond[i] - phi );
			double sinTmp = sin( omega * cond[i] - phi );
			sumCos += cosTmp * (obs[i] - obsMean); /* remove mean value */
			sumSin += sinTmp * (obs[i] - obsMean);
			sumCos2 += cosTmp * cosTmp;
			sumSin2 += sinTmp * sinTmp;
		}
		p[idx] = sumCos * sumCos / sumCos2 + sumSin * sumSin / sumSin2;
		if (maxP < p[idx])
		{
			maxP = p[idx];
			bestIdx = idx;
		}
	}

	{/* phase determination */
		sumCos = 0.0,  sumSin = 0.0;
		omega = 2 * M_PI * fVec[bestIdx];
		for (unsigned int i = 0; i < N; i++)
		{
			double cosTmp = cos( omega * cond[i] );
			double sinTmp = sin( omega * cond[i] );
			sumCos = sumCos + cosTmp * (obs[i] - obsMean);
			sumSin = sumSin + sinTmp * (obs[i] - obsMean);
		}
		phi = atan2( -sumSin, sumCos );

		while (phi < 0.) phi += 2 * M_PI;
		while (phi > 2 * M_PI) phi -= 2 * M_PI;
		*bestPhi = phi;
	}

	return bestIdx;
}