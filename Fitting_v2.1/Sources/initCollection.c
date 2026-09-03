/****************************************************************
 *
 * File........:	init_collection.c
 * Function....:	parameter initialisation for
 *                    different functions
 * Author......:	Tilo Strutz
 * last changes:	02.07.2009, 30.09.2009, 08.01.2010, 18.02.2010,
 * 25.02.2021 new: f(x|a) = (a1*x)^a2 + a3 and f(x|a) = a2 * log(a1*x)
 * 11.05.2026 new: f(x|a) = a1 / (x+a2)
 * 31.05.2026 new: init_cosine() frequency with LombScargle only
 * 28.07.2026 errmsg
 * 
 * LICENCE DETAILS: see software manual
 *	free academic use
 *  cite source as
 * "Strutz, T.: Data Fitting and Uncertainty. 3rd edition,
 *   Springer Nature, 2026"
 *
 ****************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h > /* for uint64_t */
#include "functions.h"
#include "macros.h"
#include "defines.h"
#include "prototypes.h"
#include "matrixUtils.h" /* for vector() */
#include "errmsg.h"


#ifdef __linux__
#include <x86intrin.h> // for __rdtsc
#include <unistd.h>
#endif

#ifndef WIN32
#include <sys/time.h>
#else
#include <time.h>
#define random rand
#endif

#define eTOL (10e-8)

#ifdef TIME_MEASUREMENT
#define myfprintf //fprintf  /*  suppressing time-consuming outputs    */
#else
#define myfprintf fprintf  
#endif


/*---------------------------------------------------------------
 *	init_hyperbel()
	*  f(x|a) = a1 / (x+a2)
 *--------------------------------------------------------------*/
int init_hyperbel( unsigned int N, double* obs, double* cond,
	double* a, unsigned char* a_flag, FILE* out )
{
	(void)out;
	(void)obs;
	double minCond;
	if (!a_flag[0])
	{
		a[0] = 1.0; 
	}
	if (!a_flag[1])
	{
		minCond = cond[0];
		for (unsigned int i = 0; i < N; i++)
			if (minCond < cond[i]) minCond = cond[i];
		a[1] = 1. + minCond;
	}
	return 0;
}

/*---------------------------------------------------------------
 *	init_fexpBaseX()
 *  f(x|a) = (a1*x)^a2 + a3
 *--------------------------------------------------------------*/
int init_fexpBaseX( unsigned int N, double* obs, double* cond,
	double* a, unsigned char* a_flag, FILE* out )
{
	(void)out;
	(void)cond;
	double minObs;
	
	if (!a_flag[2])
	{ /* offset is about the smallest observation */
		minObs = obs[0];
		for (unsigned int i = 0; i < N; i++)
			if (minObs > obs[i]) minObs = obs[i];
		a[2] = minObs;
	}
	if (!a_flag[1])
	{ /* exponent */
		if (obs[N - 1] > obs[0])
		{/* increasing function */
			if (fabs(obs[N - 1] - obs[N >> 1]) > fabs( obs[0] - obs[N >> 1] ))
			{/* increasing slope */
				a[1] = 2.0;
			}
			else
			{/* decreasing slope */
				a[1] = 0.5;
			}
		}
		else
		{/* decreasing function */
			a[1] = -0.5;
		}
	}
	if (!a_flag[0])
	{
		a[0] = 1.; /* scale */
	}
	return 0;
}

/*---------------------------------------------------------------
 *	init_polynomial()
 * f(x|a) = sum_{j=1}^M aj * x^(j-1)
 *--------------------------------------------------------------*/
int init_polynomial( unsigned int N, double* obs, double* cond,
	double* a, unsigned char* a_flag, FILE* out )
{
	(void)N;
	(void)out;
	(void)obs;
	(void)cond;
	long j;

	if (!a_flag[0])
		a[0] = ((double)rand() / RAND_MAX - 0.5) * 10.01;
	if (!a_flag[1])
		a[1] = ((double)rand() / RAND_MAX - 0.5) * 10.01;

	/* assume maximal number of parameters */
	for (j = 2; j < M_MAX; j++)
	{
		if (!a_flag[j])
			a[j] = ((double)rand() / RAND_MAX - 0.5) * 10.01;
	}
	return 0;
}

 /*---------------------------------------------------------------
	*	getPhaseDiff()
	*--------------------------------------------------------------*/
double getPhaseDiff( double a, double b )
{
	double diff = fmod( fabs( a - b ), 2 * M_PI );
	if (diff > M_PI)
	{
		diff = 2 * M_PI - diff;
	}
	return diff;
}

/*---------------------------------------------------------------
 *	init_cosineNonLin()
 * f(x|a) = a1 + a2 * cos( x - a3)  (omega = 2*pi)
 *--------------------------------------------------------------*/
int init_cosineNonLin( unsigned int N, double* obs, double* cond,
	double* a, unsigned char* a_flag, FILE* out )
{
	(void)out;
	(void)cond;
	double sumObs;
	double obsMax, obsMin, obsMean;

	/* mean value of observations for a[0] */
	sumObs = obs[0];

	obsMax = obsMin = obs[0];
	for (unsigned int i = 1; i < N; i++)
	{
		sumObs += obs[i]; /* accumulate observation values */
		if (obsMax < obs[i]) { obsMax = obs[i]; }
		if (obsMin > obs[i]) { obsMin = obs[i]; }
	}

	obsMean = sumObs / N; /* mean is best choice for many samples; */

	if (!a_flag[0])  /* if not set on command line */
	{
		/*  few samples can bias one half-wave leading to unfortune parameter */
		a[0] = obsMean;
	}

	if (!a_flag[1])
	{ /* estimation of a2 = radius */
		a[1] = 0.5 * (obsMax - obsMin);
		/* systematic bias towards too high amplitudes, but no problem */
	}
	if (!a_flag[2])
	{	/* estimation of a3 = phase shift */
		a[2] = 0.0;
	}
	return 0;
}

/*---------------------------------------------------------------
	*	init_cosine()
	* f(x|a) = a1 + a2 * cos( a3*x + a4)
	* f(x|a) = a[0] + a[1] * cos( a[2]*x + a[3])
	*--------------------------------------------------------------*/
int
	init_cosine( unsigned int N, double* obs, double* cond,
		double* a, unsigned char* a_flag, FILE * out )
{
	(void)cond;
	double sumObs;
	double obsMax, obsMin, obsMean;
	
	/* mean value of observations for a[0] */
	sumObs = obs[0];

	/* estimation of a2 = amplitude = radius based on range of observations */
	obsMax = obsMin = obs[0];
	for (unsigned int i = 1; i < N; i++)
	{
		sumObs += obs[i]; /* accumulate observation values */
		if (obsMax < obs[i])		{			obsMax = obs[i]; 		}
		if (obsMin > obs[i])		{			obsMin = obs[i]; 		}
	}

	obsMean = sumObs / N; /* mean is best choice for many samples; */
	//obsMean = (obsMax + obsMin)/ (2); /* better for skewed data */

	if (!a_flag[0])/* if not set on command line */
	{ 
		/*  few samples can bias one half-wave leading to unfortune parameter */
		a[0] = obsMean;
	}

	if (!a_flag[1])
	{
		a[1] = 0.5 * (obsMax - obsMin);
		/* systematic bias towards too high amplitudes, but no problem */
	}

	/* estimation of a3 = frequency = 1/ period
	 * based on Lomb-Scargle periodogram
	 */
	int numFrequencies, bestIdx;
	double deltaF, minF, maxF, Ttot;
	double *p=NULL, * fVec=NULL, bestPhi;
	myfprintf( out, "\n# Lomb-Scargle " );

	Ttot = cond[N - 1] - cond[0];
	deltaF = 1. / Ttot; /* no oversampling */ 
	minF = deltaF;  
	maxF = N / (2* Ttot); /* worst case: fs= two samples period (2/1sec => fmax = fs/2= 1/1sec   */
	myfprintf( out, "\n# frequency grid: fMin=%.4f  fDelta=%.4f  fMax=%.4f  ", minF, deltaF, maxF );

	numFrequencies = 1 + (int)ceil( (maxF - minF) / deltaF );
	p = vector( numFrequencies );
	fVec = vector( numFrequencies );
	fVec[0] = minF;
	for (int i = 1; i < numFrequencies; i++)
	{
		fVec[i] = fVec[i - 1] + deltaF;
	}

	bestIdx = lombScargle( N, obs, cond,
		obsMean, numFrequencies, fVec, p, &bestPhi );

	if (!a_flag[2])
	{
		a[2] = 2 * M_PI * fVec[bestIdx];
	}

	if (!a_flag[3])
	{/* phase determination */
		a[3] = bestPhi;
	}

	free_vector( &p );
	free_vector( &fVec );

	return 0;
}

/*---------------------------------------------------------------
 *	init_cosine2()
 * f(x|a) = a1 + a2*cos(a3*x-a4) + a5*cos(2*a3*x-a6)
 *--------------------------------------------------------------*/
int init_cosine2( unsigned int N, double* obs, double* cond,
	double* a, unsigned char* a_flag, FILE* out )
{
	(void)cond;
	double sumObs;
	double obsMax, obsMin, obsMean;
	double *p, *fVec;

	/* mean value of observations for a[0] */
	sumObs = obs[0];

	/* estimation of a2 = amplitude = radius based on range of observations */
	obsMax = obsMin = obs[0];
	for (unsigned int i = 1; i < N; i++)
	{
		sumObs += obs[i]; /* accumulate observation values */
		if (obsMax < obs[i]) { obsMax = obs[i]; }
		if (obsMin > obs[i]) { obsMin = obs[i]; }
	}

	obsMean = sumObs / N; /* mean is best choice for many samples; */
	//obsMean = (obsMax + obsMin)/ (2); /* better for skewed data */

	if (!a_flag[0])/* if not set on command line */
	{
		/*  few samples can bias one half-wave leading to unfortune parameter */
		a[0] = obsMean;
	}

	/* estimation of dominant frequency = 1/ period
	 * based on Lomb-Scargle periodogram
	 * assume that this is zalso th elower frequency
	 */
	int bestIdx, numFrequencies;

	double deltaF, minF, maxF, Ttot, bestPhi;
	myfprintf( out, "\n# Lomb-Scargle " );

	Ttot = cond[N - 1] - cond[0];
	deltaF = 1. / Ttot; /* no oversampling */
	minF = deltaF;
	maxF = N / (2 * Ttot); /* worst case: fs= two samples period (2/1sec => fmax = fs/2= 1/1sec   */
	myfprintf( out, "\n# frequency grid: fMin=%.4f  fDelta=%.4f  fMax=%.4f  ", minF, deltaF, maxF );

	numFrequencies = 1+ (int)ceil((maxF - minF) / deltaF);
	p = vector( numFrequencies);
	fVec = vector( numFrequencies );
	fVec[0] = minF;
	for (int i = 1; i < numFrequencies; i++)
	{
		fVec[i] = fVec[i-1] + deltaF;
	}
	bestIdx = lombScargle( N, obs, cond,
		obsMean, numFrequencies, fVec, p, &bestPhi );

	if (!a_flag[1])
	{ /* default */
		a[1] = 0.5 * (obsMax - obsMin);
		/* systematic bias towards too high amplitudes, but no problem */
	}
	if (!a_flag[4])
	{ /* default */
		a[4] = 0.5 * (obsMax - obsMin);
	}

	if (!a_flag[2])
	{
		a[2] = 2 * M_PI * fVec[bestIdx];
	}

	if (!a_flag[3])
	{/* phase determination */
		a[3] = bestPhi;
	}
	
	if (!a_flag[4])
	{/* phase determination */
		a[5] = 0; /* dummy */
	}

	free_vector( &p );
	free_vector( &fVec );

	return 0;
}

/*---------------------------------------------------------------
 *	init_cosine_trend()
 * 12: f(x|a) = a1 + a2 * x + a3 * cos( x - a4)
 *-------------------------------------------------------------*/
int init_cosine_trend( unsigned int N, double* obs, double* cond,
	double* a, unsigned char* a_flag, FILE* out )
{
	(void)cond;
	(void)out;
	double mean;

	/* mean value for a[0] */
	if (!a_flag[0])
	{
		mean = 0;
		for (unsigned int i = 0; i < N; i++)
			mean += obs[i];
		a[0] = mean / N;
	}

	/* estimation of a3 = linear trend */
	if (!a_flag[1])
	{
		a[1] = 0.;	/* dummy */
	}

	/* estimation of a3 = radius */
	if (!a_flag[2]) /* if not set on command line */
	{
		mean = 0;
		for (unsigned int i = 0; i < N; i++)
			mean += obs[i] * sqrt( 2. );
		a[2] = mean / N;
	}

	/* estimation of a4 = phase shift */
	if (!a_flag[3])
	{
		a[3] = 0.;	/* dummy */
	}

	return 0;
}

/*---------------------------------------------------------------
 *	init_logarithmic()
 * f(x|a) = a1 + a2 * log( a3 + x)
 *--------------------------------------------------------------*/
int init_logarithmic( unsigned int N, double* obs, double* cond,
	double* a, unsigned char* a_flag, FILE* out )
{
	(void)out;
	double obsCloseToZero, condCloseToZero, obsMin, obsMax;
	unsigned int obsMinIdx, obsMaxIdx;

	obsCloseToZero = fabs( obs[0] );
	condCloseToZero = cond[0];
	obsMin = obsMax = obs[0];
	obsMinIdx = obsMaxIdx = 0;
	for (unsigned int i = 0; i < N; i++)
	{
		if (obsCloseToZero > fabs( obs[i] ))
		{
			obsCloseToZero = fabs( obs[i] );
			condCloseToZero = cond[i];
		}
		if (obsMin > obs[i])
		{
			obsMin = obs[i];
			obsMinIdx = i;
		}
		if (obsMax < obs[i])
		{
			obsMax = obs[i];
			obsMaxIdx = i;
		}
	}

	if (!a_flag[2])
	{
		a[2] = -condCloseToZero;
		if (cond[0] + a[2] <= 0)	a[2] = -cond[0] + 10e-5;

		/* 0= a1 + a2*log(x+a3) ==> x+a3 >= 0 horizontal shift */
	}

	{
		int variant = 0;

		if (obsMaxIdx > obsMinIdx)
		{
			/*
				*   ____
				*  /
				* /
				*           */
			if (!a_flag[0])
			{
				a[0] = obsMin;
			}
		}
		else
		{
			/*
				*   \
				*    \
				*     ----
				*           */
			variant = 1;
			if (!a_flag[0])
			{
				a[0] = obsMax;
			}
		}
		if (!a_flag[1])
		{
			double sum = 0.0;
			for (unsigned int i = 0; i < N; i++)
			{
				sum += (obs[i] - a[0]) / log( a[2] + cond[i] );
			}
			if (variant == 0) 		a[1] = sum / N;
			else		a[1] = - sum / N;
		}
	}

	return 0;
}

/*---------------------------------------------------------------
 *	init_exponentialdummy()
 *  f(x|a) = a1 + a2 * exp( a3 * x)
 *--------------------------------------------------------------*/
int init_exponentialdummy( unsigned int N, double* obs, double* cond,
	double* a, unsigned char* a_flag, FILE* out )
{
	int variant = 0; /* curvature of graph 1..4  = A..D */
	unsigned int idxCloseToMiddle;
	double bestDiff;
	double obsMin, obsMax, obsMid, condMid;
	double obsLeft, obsRight;

	idxCloseToMiddle = 0;
	bestDiff = fabs( cond[N - 1] - cond[idxCloseToMiddle] );
	obsMin = obsMax = obs[0];
	for (unsigned int i = 1; i < N; i++)
	{
		/* search position closest to middle of all conditions */
		double diff; /* compare distances to first(minimum) and last(maximum) condition */
		diff = fabs( (cond[i] - cond[0]) - (cond[N - 1] - cond[i]) );
		/* is more robust point than (N>>1) when observations are not
		 * evenly distributed along condition axis */
		if (bestDiff > diff)
		{ /* better point found */
			idxCloseToMiddle = i;
			bestDiff = diff;
		}
		if (obsMin > obs[i])	obsMin = obs[i]; /* take better value */
		if (obsMax < obs[i])	obsMax = obs[i];
	}
	/* increase robustness against ouliers by computing the avarage */
	obsLeft = obs[0]; // (obs[0] + obs[1]) * 0.5;
	obsMid = (obs[idxCloseToMiddle - 1] + obs[idxCloseToMiddle] + obs[idxCloseToMiddle + 1]) / 3.0;
	condMid = cond[idxCloseToMiddle];
	obsRight = obs[N - 1]; // (obs[N - 2] + obs[N - 1]) * 0.5;

	/*------ determine curvature ------------------*/
	if (obsLeft > obsRight) /* decreasing function */
	{
		fprintf( out, "\n# decreasing function" );
		if (fabs( obsLeft - obsMid ) > fabs( obsRight - obsMid ))
		{
			variant = 1; /* A: �.. */
			fprintf( out, ", negative a3, positive scaling factor a2" );
			if (!a_flag[1])
				a[1] = +1; a[2] = -0.5; /* dummy values */
		}
		else
		{
			variant = 3; /* C: 같. */
			fprintf( out, ", positive a3, negative scaling factor a2" );
			if (!a_flag[1])
				a[1] = -1; a[2] = 0.5; /* dummy values */
		}
	}
	else /* increasing function */
	{
		fprintf( out, "\n# increasing function" );
		if (fabs( obsLeft - obsMid ) > fabs( obsRight - obsMid ))
		{
			variant = 4; /* D: .같 */
			fprintf( out, ", negative a3, negative scaling factor a2" );
			if (!a_flag[1])
				a[1] = -1; a[2] = -0.5; /* dummy values */
		}
		else
		{
			variant = 2; /* B: ..� */
			fprintf( out, ", positive a3, positive scaling factor a2" );
			if (!a_flag[1])
				a[1] = +1; a[2] = 0.5; /* dummy values */
		}
	}
	/*------ estimating offset a1 ------------------*/
	/* assume that curve approaches asymptotically a1 */
	if (variant == 1 || variant == 2)
	{
		if (!a_flag[0])
			a[0] = obsMin;
	}
	else
	{
		if (!a_flag[0])
			a[0] = obsMax;
	}
	return 0;
}

/*---------------------------------------------------------------
	*	init_exponentialLinearised()
	*  f(x|a) = a1 + a2 * exp( a3 * x)
	*--------------------------------------------------------------*/
int init_exponentialLinearised( unsigned int N, double* obs, double* cond,
	double* a, unsigned char* a_flag, FILE* out )
{
	char* rtn = "init_exponentialLinearised";
	int err;
	int variant = 0; /* curvature of graph 1..4  = A..D */
	unsigned int idxCloseToMiddle, M = 2;
	double bestDiff, b[2];
	double thresh, smax, smin;
	double** VinvS = NULL;	/* temporary matrix */
	double** invJ = NULL;	/* temporary matrix */
	double* s = NULL;	/* singular values */
	double** V = NULL;	/* V matrix */
	double** jacob = NULL;	/* Jacobian matrix J		*/
	double* lnObs = NULL;	/* vector for weights		*/
	double obsMin, obsMax, obsMid, condMid;
	double obsLeft, obsRight;

	jacob = matrix( N, M );	/* Jacobian */
	lnObs = vector( N );	/* observations */
	V = matrix( M, M );	/* V matrix for SVD */
	s = vector( M );	/* singular values for SVD */
	VinvS = matrix( M, M );	/* temporary matrix */
	invJ = matrix( M, N );	/* temporary matrix */

	idxCloseToMiddle = 0;
	bestDiff = fabs( cond[N - 1] - cond[idxCloseToMiddle] );
	obsMin = obsMax = obs[0];
	for (unsigned int i = 1; i < N; i++)
	{
		/* search position closest to middle of all conditions */
		double diff; /* compare distances to first(minimum) and last(maximum) condition */
		diff = fabs( (cond[i] - cond[0]) - (cond[N - 1] - cond[i]) );
		/* is more robust point than (N>>1) when observations are not
		 * evenly distributed along condition axis */
		if (bestDiff > diff)
		{ /* better point found */
			idxCloseToMiddle = i;
			bestDiff = diff;
		}
		if (obsMin > obs[i])	obsMin = obs[i]; /* take better value */
		if (obsMax < obs[i])	obsMax = obs[i];
	}
	/* increase robustness against ouliers by computing the avarage */
	obsLeft = obs[0]; // (obs[0] + obs[1]) * 0.5;
	obsMid = (obs[idxCloseToMiddle - 1] + obs[idxCloseToMiddle] + obs[idxCloseToMiddle + 1]) / 3.0;
	condMid = cond[idxCloseToMiddle];
	obsRight = obs[N - 1]; // (obs[N - 2] + obs[N - 1]) * 0.5;

	/*------ determine curvature ------------------*/
	if (obsLeft > obsRight) /* decreasing function */
	{
		fprintf( out, "\n# decreasing function" );
		if (fabs( obsLeft - obsMid ) > fabs( obsRight - obsMid ))
		{
			variant = 1; /* A: �.. */
			fprintf( out, ", negative a3, positive scaling factor a2" );
			if (!a_flag[1])
				a[1] = +1; a[2] = -0.5; /* dummy values */
		}
		else
		{
			variant = 3; /* C: 같. */
			fprintf( out, ", positive a3, negative scaling factor a2" );
			if (!a_flag[1])
				a[1] = -1; a[2] = 0.5; /* dummy values */
		}
	}
	else /* increasing function */
	{
		fprintf( out, "\n# increasing function" );
		if (fabs( obsLeft - obsMid ) > fabs( obsRight - obsMid ))
		{
			variant = 4; /* D: .같 */
			fprintf( out, ", negative a3, negative scaling factor a2" );
			if (!a_flag[1])
				a[1] = -1; a[2] = -0.5; /* dummy values */
		}
		else
		{
			variant = 2; /* B: ..� */
			fprintf( out, ", positive a3, positive scaling factor a2" );
			if (!a_flag[1])
				a[1] = +1; a[2] = 0.5; /* dummy values */
		}
	}

	/* assume that curve approaches asymptotically a1 */
	if (variant == 1 || variant == 2)
	{
		if (!a_flag[0])
			a[0] = obsMin;
	}
	else 
	{
		if (!a_flag[0])
			a[0] = obsMax;
	}

	fprintf( out, "\n# linear regression for initial parameters  #####" );

	for (unsigned int i = 0; i < N; i++)
	{
		switch (variant)
		{
		case 1:
		case 2:
			lnObs[i] = log( obs[i] - a[0] + 10e-8 );
			break;
		case 3:
		case 4:
			lnObs[i] = log( a[0] - obs[i] + 10e-8 );
			break;
		}
	}
	for (unsigned int i = 0; i < N; i++)
	{
		for (unsigned int j = 0; j < M; j++)
		{ /* a is dummy parameter needed only for nonlinear problems */
			jacob[i][j] = flin_deriv( NULL, i, j, M, cond, a );
		}
	}
	/* do the SVD	*/
	err = singvaldec( jacob, N, M, s, V ); /* jacob will contain result of U */
	/* check the singular values	*/
	smax = smin = s[0];
	for (unsigned int j = 1; j < M; j++)
	{
		if (s[j] > smax)					smax = s[j];
		if (s[j] < smin)					smin = s[j];
	}
	if (smax < TOL_S)
	{
		fprintf( stderr,
			"\n### %s\n###     singular matrix, smax = %f", rtn, smax );
		fprintf( out,
			"\n### %s\n###     singular matrix, smax = %f", rtn, smax );

		err = ERR_IS_SINGULAR;
		goto endfunc; /* use dummy values */
	}
	if (smin > TOL)
	{
		if (smax / smin > 1.e+99)
		{
			fprintf( stderr,
				"\n### %s\n###     high condition number: %e / %e = %e", rtn, smax, smin, smax / smin );
			fprintf( out,
				"\n### %s\n###     high condition number: %e / %e = %e", rtn, smax, smin, smax / smin );
			for (unsigned int i = 0; i < M; i++)
			{
				fprintf( stderr, "\n" );
				for (unsigned int j = 0; j < M; j++)
				{
					fprintf( stderr, "%10.6f ", jacob[i][j] );
				}
			}
			fprintf( stderr, "\n" );
			fflush( stderr );
			err = ERR_ANY;
			goto endfunc;
		}
	}

	thresh = min( TOL_S * smax, TOL_S );

	/* invert singular values */
	for (unsigned int j = 0; j < M; j++)
	{
		/* <= in case of smax =0 */
		if (s[j] <= thresh) 		s[j] = 0.0;
		else			s[j] = 1. / s[j];
	}

	/* V * [diag(1/s[j])] */
	for (unsigned int i = 0; i < M; i++)
	{
		for (unsigned int j = 0; j < M; j++)
		{
			VinvS[i][j] = V[i][j] * s[j];
		}
	}

	/* multiplication of VinvS with transposed of U  */
	/* result is: inv(J) = (V*inv(S)) * U'	*/
	for (unsigned int i = 0; i < M; i++)
	{
		for (unsigned int j = 0; j < N; j++)
		{
			invJ[i][j] = 0.;
			for (unsigned int n = 0; n < M; n++)
			{
				invJ[i][j] += VinvS[i][n] * jacob[j][n];
			}
		}
	}

	/* compute the parameter vector a = inv(J)*y */
	for (unsigned int j = 0; j < M; j++)
	{
		b[j] = 0.0;
		for (unsigned int i = 0; i < N; i++)
		{
			b[j] += invJ[j][i] * lnObs[i];
		}
	}

	/* remap linear parameters onto nonlinear problem */
	if (!a_flag[2])
		a[2] = b[1];
	if (!a_flag[1])
	{
		a[1] = exp( b[0] );
		if (variant > 2)
			a[1] = -a[1];
	}
endfunc:
	free_vector( &lnObs );
	free_matrix( &jacob );
	free_vector( &s );
	free_matrix( &V );
	free_matrix( &VinvS );
	free_matrix( &invJ );
	fprintf( out, "\n# a1= %8.4f  a2= %8.4f  a3= %8.4f     #####", a[0], a[1], a[2] );
	return err;
}
/*---------------------------------------------------------------
 *	init_exponential()
 *  f(x|a) = a1 + a2 * exp( a3 * x)
 *--------------------------------------------------------------*/
int init_exponential( unsigned int N, double* obs, double* cond,
	double* a, unsigned char* a_flag, FILE* out )
{
	double slope = 0., sign; /* dummy defaut */
	double bestDiff;
	unsigned int idxCloseToZero, idx1CloseToZero, idx2CloseToZero, idxCloseToMiddle;
	double condCloseToZero;
	double obsMin, obsMax, obsMid, condMid;
	double obsLeft, obsRight, obsRef, condRef;
	double b[3]; /* parameters based on alternative approach */
	int variant = 0; /* curvature of graph 1..4  = A..D */
#if defined DUMMY_EXPONENTIAL || _DEBUG
	int aParamsOk_flag = 0, bParamsOk_flag = 0;
#else
	int aParamsOk_flag = 1, bParamsOk_flag = 1;
#endif

	/*-------- derive some needed values ----------------*/
	idxCloseToZero = 0;
	condCloseToZero = fabs( cond[idxCloseToZero] );
	idxCloseToMiddle = 0;
	bestDiff = fabs( cond[N - 1] - cond[idxCloseToMiddle] );
	obsMin = obsMax = obs[0];
	for (unsigned int i = 1; i < N; i++)
	{
		/* search position closest to middle of all conditions */
		double diff; /* compare distances to first(minimum) and last(maximum) condition */
		diff = fabs( (cond[i] - cond[0]) - (cond[N - 1] - cond[i]) );
		/* is more robust point than (N>>1) when observations are not
		 * evenly distributed along condition axis */
		if (bestDiff > diff)
		{ /* better point found */
			idxCloseToMiddle = i;
			bestDiff = diff;
		}
		/* search position closest to x=0 */
		if (condCloseToZero > fabs( cond[i] ))
		{
			condCloseToZero = fabs( cond[i] );
			idxCloseToZero = i;
		}
		if (obsMin > obs[i])	obsMin = obs[i]; /* take better value */
		if (obsMax < obs[i])	obsMax = obs[i];
	}
	/* increase robustness against ouliers by computing the avarage */
	obsLeft = obs[0]; // (obs[0] + obs[1]) * 0.5;
	obsMid = obs[idxCloseToMiddle]; // (obs[idxCloseToMiddle - 1] + obs[idxCloseToMiddle] + obs[idxCloseToMiddle + 1]) / 3.0;
	condMid = cond[idxCloseToMiddle];
	obsRight = obs[N - 1]; // (obs[N - 2] + obs[N - 1]) * 0.5;

	/*------ estimating offset a1 ------------------*/
	if (obsLeft > obsRight) /* decreasing function */
	{
		fprintf( out, "\n# init: decreasing exponential detected" );
		sign = -1;
		if (fabs( obsLeft - obsMid ) > fabs( obsRight - obsMid ))
		{
			variant = 1; /* A: �.. */
			obsRef = obsLeft;
			condRef = cond[0];
			fprintf( out, ", negative a3, positive scaling factor a2" );
		}
		else
		{
			variant = 3; /* C: 같. */
			obsRef = obsRight;
			condRef = cond[N - 1];
			fprintf( out, ", positive a3, negative scaling factor a2" );
		}
	}
	else /* increasing function */
	{
		fprintf( out, "\n# init: increasing exponential detected" );
		sign = +1;
		if (fabs( obsLeft - obsMid ) > fabs( obsRight - obsMid ))
		{
			variant = 4; /* D: .같 */
			obsRef = obsLeft;
			condRef = cond[0];
			fprintf( out, ", negative a3, negative scaling factor a2" );
		}
		else
		{
			variant = 2; /* B: ..� */
			obsRef = obsRight;
			condRef = cond[N - 1];
			fprintf( out, ", positive a3, positive scaling factor a2" );
		}
	}

	/* assume that curve approaches asymptotically a1 */
	if (variant == 1 || variant == 2)
	{
		if (!a_flag[0])	a[0] = obsMin;
		if (fabs( obsMid - obsMin ) < eTOL)
			obsMid = obsMin + eTOL;
		/* this difference is later used and may not be too small */
	}
	else
	{
		if (!a_flag[0])	a[0] = obsMax;
		if (fabs( obsMid - obsMax ) < eTOL)
			obsMid = obsMax - eTOL;
	}
	/* it is not recommended to use a[0] = obs[0] or obs[N-1], because
		* a[1] depends on a[0] and its calculation could lead to wrong
		* values of a[1] (wrong sign)  */
	b[0] = a[0]; /* copy for second approach */

	/*-------- estimating scale factor a2 ------------------*/
	/* find second point with condition close to zero for slope computation */
	if (cond[idxCloseToZero] > 0.)
	{ /* keep order idx1CloseToZero < idx2CloseToZero */
		if (idxCloseToZero > 0)
		{ /* there is a point on the opposite side */
			idx1CloseToZero = idxCloseToZero-1;
			idx2CloseToZero = idxCloseToZero;
		}
		else
		{ /* we must take a point on the same side */
			idx1CloseToZero = idxCloseToZero;
			idx2CloseToZero = idxCloseToZero + 1;
		}
	}
	else
	{
		if (idxCloseToZero < N - 1)
		{ /* there is a point on the opposite side */
			idx1CloseToZero = idxCloseToZero;
			idx2CloseToZero = idxCloseToZero + 1;
		}
		else
		{ /* we must take a point on the same side */
			idx1CloseToZero = idxCloseToZero - 1;
			idx2CloseToZero = idxCloseToZero;
		}
	}

	while (fabs( cond[idx1CloseToZero] - cond[idx2CloseToZero] ) < eTOL &&
		idx2CloseToZero < N - 1 && idx1CloseToZero > 0)
	{ /* too close */
		if (idx2CloseToZero < N - 1) idx2CloseToZero++;
		else if (idx1CloseToZero > 0)idx1CloseToZero--;
		else break;
	}

	/* check for wrong data, e.g. missing # sign for comments */
	if (fabs(cond[idx1CloseToZero] - cond[idx2CloseToZero]) < TOL)
	{
		fprintf( stderr, "\n\n#### suspicious conditions #####\n" );
		fprintf( out, "\n\n#### suspicious conditions #####\n" );
	}
	while (fabs(cond[idx1CloseToZero] - cond[idx2CloseToZero]) < TOL)
	{ /* wrong conditions */
		if (idx1CloseToZero > 0)	idx1CloseToZero--;
		else if (idx2CloseToZero < N - 1) idx2CloseToZero++;
		else
		{/* dummy values */
			cond[idx2CloseToZero] = cond[idx2CloseToZero] + TOL;
			break;
		}
	}
	/* search for point pair where slope matches to decreasing/increasing behaviour  */
	slope = (obs[idx1CloseToZero] - obs[idx2CloseToZero]) / (cond[idx1CloseToZero] - cond[idx2CloseToZero]);
	while ((sign * slope) < 0.)
	{ /* wrong slope, look for a nother point pair, as long other points are available */
		if (idx1CloseToZero > 0)	idx1CloseToZero--;
		if (idx2CloseToZero < N - 1) idx2CloseToZero++;

		slope = (obs[idx1CloseToZero] - obs[idx2CloseToZero]) /
			(cond[idx1CloseToZero] - cond[idx2CloseToZero]);
		if (idx1CloseToZero == 0 && idx2CloseToZero == N - 1)
		{
			fprintf( out, "\n#### suspicious behavior of function                             #####" );
			fprintf( out, "\n#### cannot find a point pair that satisfies the needed slope sign #####" );
			fprintf( stdout, "\n#### suspicious behavior of function                             #####" );
			fprintf( stdout, "\n#### cannot find a point pair that satisfies the needed slope sign #####" );
			break;
		}
	}

	fprintf( out, "\n# slope: %0.3f", slope );

	if (!a_flag[1])
	{
		a[1] = obs[idxCloseToZero] - slope * cond[idxCloseToZero];
		a[1] -= a[0];
	}

	/*-------- estimating exponential factor a3 ------------------*/
	if (!a_flag[2])
	{
		a[2] = slope / a[1];
	}
	/* check a2 for correct sign */
	if (variant == 1 || variant == 2) /* a2 > 0 */
	{
		if (a[1] <= 0.) aParamsOk_flag = 0;
	}
	else /* a2 < 0 */
	{
		if (a[1] >= 0.) aParamsOk_flag = 0;
	}
	/* check a3 for correct sign */
	if (variant == 1 || variant == 4) /* a3 < 0 */
	{
		if (a[2] >= 0.) aParamsOk_flag = 0;
	}
	else /* a3 > 0 */
	{
		if (a[2] <= 0.) aParamsOk_flag = 0;
	}


	/*-------------- alternative estimation of a1, a2 --------------*/
	/* use obsMid, select second point */
	if (variant == 1 || variant == 2) /* (a[0] == obsMin) */
	{
		b[2] = (log( obsMid - b[0] ) - log( obsRef - b[0] ))
			/ (condMid - condRef);

	}
	else /* (a[0] == obsMax) */
	{
		b[2] = (log( b[0] - obsMid ) - log( b[0] - obsRef ))
			/ (condMid - condRef);
	}
	/* yMid = b1 + b2 * exp(b3*xMid) */
	b[1] = (obsMid - b[0]) / exp( b[2] * condMid );
	/* check a2 for correct sign */
	if (variant == 1 || variant == 2) /* b2 > 0 */
	{
		if (b[1] <= 0.) bParamsOk_flag = 0;
	}
	else /* b2 < 0 */
	{
		if (b[1] >= 0.) bParamsOk_flag = 0;
	}
	/* check b3 for correct sign */
	if (variant == 1 || variant == 4) /* b3 < 0 */
	{
		if (b[2] >= 0.) bParamsOk_flag = 0;
	}
	else /* a3 > 0 */
	{
		if (b[2] <= 0.) bParamsOk_flag = 0;
	}

	fprintf( out, "\n# a1= %8.4f  a2= %8.4f  a3= %8.4f     #####", a[0], a[1], a[2] );
	fprintf( out, "\n# b1= %8.4f  b2= %8.4f  b3= %8.4f     #####", b[0], b[1], b[2] );

	if (aParamsOk_flag)
	{
		if (bParamsOk_flag)
		{/* average of both hoping that the result is closer to the true values */
			if (!a_flag[0])			a[0] = 0.5 * (a[0] + b[0]);
			if (!a_flag[1])			a[1] = 0.5 * (a[1] + b[1]);
			if (!a_flag[2])			a[2] = 0.5 * (a[2] + b[2]);
		}
		/* else keep a version, nothing to change */
	}
	else if (bParamsOk_flag)
	{/* take b version */
		if (!a_flag[0]) 			a[0] = b[0];
		if (!a_flag[1]) 			a[1] = b[1];
		if (!a_flag[2]) 			a[2] = b[2];
	}
	else
	{
		fprintf( out, "\n#\n# estimation of initial parameters failed ! #####" );
		fprintf( out, "\n# use dummy initial parameters ! #####" );
		/* keep a1 */;
		if (variant == 1) { if (!a_flag[1]) a[1] = 1; if (!a_flag[2]) a[2] = -0.5; }
		else if (variant == 2) { if (!a_flag[1]) a[1] = 1; if (!a_flag[2]) a[2] = 0.5; }
		else if (variant == 3) { if (!a_flag[1]) a[1] = -1; if (!a_flag[2]) a[2] = 0.5; }
		else if (variant == 4) { if (!a_flag[1]) a[1] = -1; if (!a_flag[2]) a[2] = -0.5; }
#if defined DUMMY_EXPONENTIAL || _DEBUG
		if (!a_flag[0]) a[0] = 0.;
		if (!a_flag[1]) a[1] = 1.;
		if (!a_flag[2]) a[2] = 0.1;
#endif
	}
	fprintf( out, "\n# a1= %8.4f  a2= %8.4f  a3= %8.4f     #####", a[0], a[1], a[2] );

	return 0;
}

/*---------------------------------------------------------------
 *	init_expon2()
 * f(x|a) = a1 * exp( a2 * x)
 *--------------------------------------------------------------*/
int init_expon2( unsigned int N, double* obs, double* cond,
	double* a, unsigned char* a_flag, FILE* out )
{
	(void)out;
	double slope1 = 0, slope2 = 0; /* dummy defaut */
	double condCloseToZero;
	unsigned int idx1CloseToZero;

	/* search position closest to x=0 */
	condCloseToZero = fabs( obs[0] );
	idx1CloseToZero = 0;
	for (unsigned int i = 0; i < N; i++)
	{
		if (condCloseToZero > fabs( cond[i] ))
		{
			condCloseToZero = fabs( cond[i] );
			idx1CloseToZero = i;
		}
	}
	if (cond[N - 1] - cond[idx1CloseToZero] > 10e-5)
	{
		slope1 = (obs[idx1CloseToZero] - obs[N - 1]) / (cond[idx1CloseToZero] - cond[N - 1]);
	}
	if (cond[idx1CloseToZero] - cond[0] > 10e-5)
	{
		slope2 = (obs[idx1CloseToZero] - obs[0]) / (cond[idx1CloseToZero] - cond[0]);
	}
	if (!a_flag[0])
	{
		a[0] = 0.5 * (obs[idx1CloseToZero] - slope1 * cond[idx1CloseToZero]
			+ obs[idx1CloseToZero] - slope2 * cond[idx1CloseToZero]);
	}

	if (!a_flag[1])
	{
		if (slope2 < TOL)		a[1] = 1 / a[0] * slope1;
		else if (slope1 < TOL)		a[1] = 1 / a[0] * slope2;
		else 	a[1] = 0.5 * (1 / a[0] * slope1 + 1 / a[0] * slope2);
	}
	return 0;
}

/*---------------------------------------------------------------
 *	init_gauss()
 * f(x|a) = a1 * exp( a2 * (x-a3)^2) +
 *--------------------------------------------------------------*/
int init_gauss( unsigned int N, double* obs, double* cond,
	double* a, unsigned char* a_flag, FILE* out )
{
	char* rtn = "init_gauss()";
	int err = 0;		 /* return value */
	int i_mean = 0, i_max, i_min;
	double max_val, min_val, condmin = 0., condmax = 0.;
	double mean, var, sum, sigma, tmp;

	/*
	 * get starting point
	 * assuming that one Gaussian is good enough to fit the data
	 */

	 /* get peak of curve */
	max_val = min_val = obs[1];
	i_max = i_min = 1;
	condmax = cond[1];
	condmin = cond[1];
	for (unsigned int i = 2; i < N - 1; i++) /* let 1 sample border */
	{
		if (max_val < obs[i])
		{
			max_val = obs[i];
			i_max = i;	/* peak position index */
			condmax = cond[i];	/* peak position */
		}
		if (min_val > obs[i])
		{
			min_val = obs[i];
			i_min = i;	/* peak position */
			condmin = cond[i];
		}
	}
	if (max_val == min_val)
	{
		fprintf( out, "\n\n Nothing to fit !!" );
		a[0] = 0.;
		a[2] = -50000000.0;
		a[1] = 0.;
		err = errmsg( ERR_DATA_PROBLEMS, rtn, "Nothing to fit", 0 );
		goto endfunc;
	}

	mean = sum = var = 0.;
	/* take only that part which has the highest peak */
	if (fabs( max_val ) > fabs( min_val ))
	{
		/* positive amplitude */
		for (unsigned int i = 0; i < N; i++)
		{
			if (obs[i] > 0.)
			{
				/* mean and variance of condition
				 * observed value is like probability
				 */
				tmp = cond[i] * obs[i];
				mean += tmp;
				var += cond[i] * tmp;
				sum += obs[i];
			}
		}
		if (sum > 0.)
		{
			mean /= sum;	/* average along cond[i] */
			var = var / sum - mean * mean;
		}
	}
	else
	{
		/* negative amplitude */
		for (unsigned int i = 0; i < N; i++)
		{
			if (obs[i] < 0.)
			{
				tmp = -cond[i] * obs[i];
				mean += tmp;
				var += cond[i] * tmp;
				sum -= obs[i];
			}
		}
		if (sum > 0.)
		{
			mean /= sum;
			var = var / sum - mean * mean;
		}
	}

	/* if only one data point, then sigma is zero */
	if (var > 0.)	sigma = sqrt( var );  /* deviation of Gaussian */
	else
		sigma = 0.0000001;


	/* get index of mean position */
	for (unsigned int i = 1; i < N; i++)
	{
		if (cond[i - 1] <= mean && mean <= cond[i])
		{
			i_mean = i;	/* mean position */
			break;
		}
	}

	/* make values more robust by averaging */
	max_val = (max_val + obs[i_max - 1] + obs[i_max + 1]) / 3;
	min_val = (min_val + obs[i_min - 1] + obs[i_min + 1]) / 3;
	if (obs[i_mean] > 0.)
	{
		{
			/* select highest peak, when there are 2 or more */
			if (!a_flag[0])	a[0] = max_val;
			if (!a_flag[1])	a[1] = condmax;
			/* reduce deviation accordingly */
			if (sigma > fabs( mean - condmax ))
				sigma -= fabs( mean - condmax );
		}
	}
	else
	{
		{
			if (!a_flag[0])	a[0] = min_val;
			if (!a_flag[1])	a[1] = condmin;
			if (sigma > fabs( mean - condmin )) sigma -= fabs( mean - condmin );
		}
	}
	/* transcode deviation */
	if (!a_flag[2])	a[2] = -0.5 / (sigma * sigma);

endfunc:
	return err;
}

/*---------------------------------------------------------------
 *	init_gen_laplace()
 * f(x|a) = a1 * exp( -|x|^a2 * a3)
 *--------------------------------------------------------------*/
int init_gen_laplace( unsigned int N, double* obs, double* cond,
	double* a, unsigned char* a_flag, FILE* out )
{
	(void)N;
	(void)out;
	(void)obs;
	(void)cond;
	/* assumes conditions starting close to zero
	 * y(x=0) = a1 * exp( 0) = a1
	 */
	if (!a_flag[0])
	{
		a[0] = obs[0];
	}
	if (!a_flag[1])
	{
		a[1] = 1.0;
	}
	if (!a_flag[2])
	{
		a[2] = 0.8;
	}

	return 0;
}

/*---------------------------------------------------------------
 *	init_circle()
 *  f(x|a) = 0 = (x1-a1)^2 + (x2-a2)^2 - a3^2
 *--------------------------------------------------------------*/
int init_circle( unsigned int N, double* obs, double* cond,
	double* a, unsigned char* a_flag, FILE* out )
{
	(void)obs;
	double sum_x, sum_y, rad2, diff1, diff2;

	fprintf( out, "\n#\n# init_circle()" );

	/*
	 * determine circle centre
	 */

	 /* compute centroids of conditions */
	sum_x = sum_y = 0.;
	/* two conditions */
	for (unsigned int i = 0; i < 2 * N; i += 2)
	{
		sum_x += cond[i];
		sum_y += cond[i + 1];
	}
	sum_x /= (double)N;
	sum_y /= (double)N;

	rad2 = 0;
	for (unsigned int i = 0; i < 2 * N; i += 2)
	{
		diff1 = cond[i] - sum_x;
		diff2 = cond[i + 1] - sum_y;
		rad2 += sqrt( diff1 * diff1 + diff2 * diff2 );
	}
	rad2 = rad2 / (double)N;
	fprintf( out, "\n#\n# mean of condition coordinates" );
	fprintf( out, "\n#   mean(x)= %f", sum_x );
	fprintf( out, "\n#   mean(y)= %f", sum_y );
	fprintf( out, "\n#   radius = %f", rad2 );

	if (!a_flag[0])	a[0] = sum_x;
	if (!a_flag[1])	a[1] = sum_y;
	if (!a_flag[2])	a[2] = rad2;

	fprintf( out,
		"\n#  f(x|a) =0= (x1-%f)**2 + (x2-%f)**2 - %f**2",
		a[0], a[1], a[2] );

	return 0;
}

/*---------------------------------------------------------------
 *	init_rotation()
 * 21... f1(x|a) = a1 + cos(a3) * x1 - sin(a3) * x2
 *       f2(x|a) = a2 + sin(a3) * x1 + cos(a3) * x2
 *-------------------------------------------------------------*/
int init_rotation( unsigned int N, double* obs, double* cond,
	double* a, unsigned char* a_flag, FILE* out )
{
	(void)obs;
	double sum_x, sum_y, sum_u, sum_v;

	fprintf( out, "\n#\n# init_rotation()" );

	/*
	 * determine rough translation
	 */

	 /* compute centroids of conditions and observations */
	sum_x = sum_y = 0;
	sum_u = sum_v = 0;
	/* assume double observations and conditions */
	for (unsigned int i = 0; i < N * 2; i += 2)
	{
		sum_x += obs[i];
		sum_y += obs[i + 1];
		sum_u += cond[i];
		sum_v += cond[i + 1];
	}
	sum_x /= (double)N;
	sum_y /= (double)N;
	sum_u /= (double)N;
	sum_v /= (double)N;
	fprintf( out, "\n#\n# mean of condition coordinates" );
	fprintf( out, "\n#   mean(u)= %f", sum_u );
	fprintf( out, "\n#   mean(v)= %f", sum_v );
	fprintf( out, "\n# mean of observed coordinates" );
	fprintf( out, "\n#   mean(x)= %f", sum_x );
	fprintf( out, "\n#   mean(y)= %f", sum_y );

	if (!a_flag[0])	a[0] = sum_x - sum_u;
	if (!a_flag[1])	a[1] = sum_y - sum_v;
	if (!a_flag[2])	a[2] = 0;		/* assume no rotation */


	fprintf( out,
		"\n#  f1(u,v) = %f + cos(%f) * u - sin(%f) * v",
		a[0], a[2], a[2] );
	fprintf( out,
		"\n#  f2(u,v) = %f + sin(%f) * u + cos(%f) * v",
		a[1], a[2], a[2] );

	return 0;
}

/*---------------------------------------------------------------
 *	init_NN3x3x1()
 * 3x3x1
 *--------------------------------------------------------------*/
int init_NN3x3x1( unsigned int N, double* obs, double* cond,
	double* a, unsigned char* a_flag, FILE* out )
{
	(void)out;
	(void)obs;
	const int  condDim = 3; /* because there are three input units */
	double minval[3]={0}, maxval[3] = { 0 };
#ifndef WIN32
	struct timeval tv;
	struct timezone tz;

	gettimeofday( &tv, &tz );
	srandom( tv.tv_sec );
#else
	/* Seed the random-number generator with current time so that
	 * the numbers will be different every time we run.
	 */
	srand( (unsigned)time( NULL ) );
#endif

	/* give parameters random values */
	for (int j = 0; j < M_MAX; j++)
	{
		if (!a_flag[j])
			a[j] = 2. * (float)random() / (float)RAND_MAX - 1.;
	}

	/* make random numbers in a range that |cond x param| < 5 */
	minval[0] = maxval[0] = cond[0];
	minval[1] = maxval[1] = cond[1];
	minval[2] = maxval[2] = cond[2];
	for (unsigned int i = 1; i < N; i++)
	{
		for (int j = 0; j < condDim; j++)
		{
			if (minval[j] > cond[condDim * i + j]) minval[j] = cond[condDim * i + j];
			if (maxval[j] < cond[condDim * i + j]) maxval[j] = cond[condDim * i + j];
		}
	}
	/* weights from 1st input */
	if (!a_flag[1])	a[1] = a[1] / (maxval[0] - minval[0]);
	if (!a_flag[5])	a[5] = a[5] / (maxval[0] - minval[0]);
	if (!a_flag[9]) a[9] = a[9] / (maxval[0] - minval[0]);

	/* weights from 2nd input */
	if (!a_flag[2])	a[2] = a[2] / (maxval[1] - minval[1]);
	if (!a_flag[6])	a[6] = a[6] / (maxval[1] - minval[1]);
	if (!a_flag[10])a[10] = a[10] / (maxval[1] - minval[1]);

	/* weights from 3rd input */
	if (!a_flag[3])	a[3] = a[3] / (maxval[2] - minval[2]);
	if (!a_flag[7])	a[7] = a[7] / (maxval[2] - minval[2]);
	if (!a_flag[11])a[11] = a[11] / (maxval[2] - minval[2]);

	return 0;
}

/*---------------------------------------------------------------
 *	init_NN1x3x1()
 * 1x3x1
 *--------------------------------------------------------------*/
int init_NN1x3x1( unsigned int N, double* obs, double* cond,
	double* a, unsigned char* a_flag, FILE* out )
{
	(void)out;
	(void)obs;
	double minval, maxval;
#ifndef WIN32
	struct timeval tv;
	struct timezone tz;

	gettimeofday( &tv, &tz );
	srandom( tv.tv_sec );
#else
	/* Seed the random-number generator with current time so that
	 * the numbers will be different every time we run.
	 */
	srand( (unsigned)time( NULL ) );
#endif

	/* give parameters random values */
	for (unsigned int j = 0; j < M_MAX; j++)
	{
		if (!a_flag[j])
		{
			a[j] = 2. * (float)random() / (float)RAND_MAX - 1.;
		}
	}

	/* make random numbers in a range that |cond x param| < 5 */
	minval = maxval = cond[0];
	for (unsigned int i = 1; i < N; i++)
	{
		if (minval > cond[i]) minval = cond[i];
		if (maxval < cond[i]) maxval = cond[i];
	}
	/* weights from 1st input */
	if (!a_flag[1])	a[1] = a[1] / (maxval - minval);
	if (!a_flag[3])	a[3] = a[3] / (maxval - minval);
	if (!a_flag[5])	a[5] = a[5] / (maxval - minval);

	return 0;
}

/*---------------------------------------------------------------
 *	init_NN()
 * 3x...
 *--------------------------------------------------------------*/
int init_NN( unsigned int N, double* obs, double* cond,
	double* a, unsigned char* a_flag, FILE* out )
{
	(void)N;
	(void)out;
	(void)obs;
	(void)cond;
	int j;
#ifndef WIN32
	struct timeval tv;
	struct timezone tz;

	gettimeofday( &tv, &tz );
	srandom( tv.tv_sec );
#else
	/* Seed the random-number generator with current time so that
	 * the numbers will be different every time we run.
	 */
	srand( (unsigned)time( NULL ) );
#endif

	for (j = 0; j < M_MAX; j++)
	{
		if (!a_flag[j])
		{
			a[j] = 1. * (float)random() / (float)RAND_MAX - 0.5;
		}
	}
	return 0;
}

/*---------------------------------------------------------------
 *	init_tinyExp()
 *  for variable model functions (provided by the user)
 *  only non-negativ parameters assumed
 *--------------------------------------------------------------*/
int init_tinyExp( unsigned int N, double* obs, double* cond,
	double* a, unsigned char* a_flag, FILE* out )
{
	(void)N;
	(void)out;
	(void)obs;
	(void)cond;
	int j;
#ifndef WIN32
	struct timeval tv;
	struct timezone tz;

	gettimeofday( &tv, &tz );
	srandom( tv.tv_sec );
#else
	/* Seed the random-number generator with current time so that
	 * the numbers will be different every time we run.
	 */
	srand( (unsigned)time( NULL ) );
#endif

	for (j = 0; j < M_MAX; j++)
	{
		if (!a_flag[j])
		{
			//a[j] = f * (float)random() / (float)RAND_MAX;
			a[j] = 1.;
		}
	}
	return 0;
}

/*---------------------------------------------------------------
 *	init_Example()
 *--------------------------------------------------------------*/
int init_Example( unsigned int N, double* obs, double* cond,
	double* a, unsigned char* a_flag, FILE* out )
{
	(void)N;
	(void)out;
	(void)obs;
	(void)cond;
	if (!a_flag[0])
	{
		a[0] = 2.25;
	}
	return 0;
}