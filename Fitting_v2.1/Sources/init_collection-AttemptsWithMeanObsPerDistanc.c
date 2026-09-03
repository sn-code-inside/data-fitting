/****************************************************************
 *
 * File........:	init_collection.c
 * Function....:	parameter initialisation for
 *                    different functions
 * Author......:	Tilo Strutz
 * last changes:	02.07.2009, 30.09.2009, 08.01.2010, 18.02.2010
 * 25.02.2021 new: f(x|a) = (a1*x)^a2 + a3 and f(x|a) = a2 * log(a1*x)
 *
 * LICENCE DETAILS: see software manual
 *	free academic use
 *  cite source as
 * "Strutz, T.: Data Fitting and Uncertainty. 2nd edition,
 *   Springer Vieweg, 2016"
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
#include "matrix_utils.h" /* for vecto() */


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


#include <stdlib.h>
#include <stdio.h>
#include <math.h>

 /* Rückgabewerte */
#define RESAMPLE_OK                 0
#define RESAMPLE_ERR_ARGS          -1
#define RESAMPLE_ERR_MONOTONICITY  -2
#define RESAMPLE_ERR_DX_TOO_SMALL  -3
#define RESAMPLE_ERR_ALLOC         -4

/* Epsilon zum Schutz gegen numerische Grenzfälle */
#ifndef EPS_DX
#define EPS_DX 1e-12
#endif

/*---------------------------------------------------------------
*	Lineare Interpolation zwischen (x0,y0) und (x1,y1) an Position x
*--------------------------------------------------------------*/
static double lerp( double x, double x0, double y0, double x1, double y1 )
{
	if (x1 == x0)
	{
		/* degenerierter Fall; gibt y0 zurück */
		return y0;
	}
	return y0 + (y1 - y0) * ((x - x0) / (x1 - x0));
}

/*---------------------------------------------------------------
 * Resampling auf gleichmäßiges Raster mit mittlerer Distanz:
 *
 * Eingabe:
 *   x, y  : Datenpunkte (x[i], y[i]), i=0..n-1 (x strikt aufsteigend)
 *   n     : Anzahl Punkte (n >= 2)
 *
 * Ausgabe:
 *   *x_out: neu allokiertes Array gleichmäßig beabstandeter x-Positionen
 *   *y_out: neu allokiertes Array der resampelten y-Werte
 *   *m_out: Anzahl der Ausgabepunkte
 *
 * Rückgabe:
 *   0 bei Erfolg; negativer Fehlercode sonst (siehe #define oben).
*--------------------------------------------------------------*/
int resample_to_mean_dx( const double* x, const double* y, int n,
	double** x_out, double** y_out, int* m_out )
{
	int i;
	double dx_sum, dx_mean, xmin, xmax, span;
	int m;
	double* xo = NULL;
	double* yo = NULL;

	/* Argumentprüfung */
	if (!x || !y || !x_out || !y_out || !m_out || n < 2)
	{
		return RESAMPLE_ERR_ARGS;
	}

	/* Monotonie prüfen und mittlere Distanz berechnen */
	dx_sum = 0.0;
	for (i = 0; i < n - 1; ++i)
	{
		double dxi = x[i + 1] - x[i];
		if (dxi <= 0.0)
		{
			/* x ist nicht strikt aufsteigend */
			return RESAMPLE_ERR_MONOTONICITY;
		}
		dx_sum += dxi;
	}
	dx_mean = dx_sum / (double)(n - 1);

	if (dx_mean <= EPS_DX)
	{
		return RESAMPLE_ERR_DX_TOO_SMALL;
	}

	xmin = x[0];
	xmax = x[n - 1];
	span = xmax - xmin;

	/* Anzahl gleichmäßig beabstandeter Stützstellen bestimmen
		 Mindestens 2 Punkte, wenn span > 0 */
	if (span < EPS_DX)
	{
		/* Alle x praktisch identisch -> gib genau zwei Punkte zurück,
			 identisch zu xmin/xmax */
		m = 2;
	}
	else
	{
		/* floor(...) + 1 stellt sicher, dass der Start enthalten ist;
			 den letzten Punkt klemmen wir später hart auf xmax */
		m = (int)floor( span / dx_mean ) + 1;
		if (m < 2) m = 2; /* Sicherheitsnetz */
	}

	xo = (double*)malloc( (size_t)m * sizeof( double ) );
	yo = (double*)malloc( (size_t)m * sizeof( double ) );
	if (!xo || !yo)
	{
		if (xo) free( xo );
		if (yo) free( yo );
		return RESAMPLE_ERR_ALLOC;
	}

	/* Gleichmäßige x_out erstellen; letzten Punkt auf xmax setzen (clamp) */
	for (i = 0; i < m; ++i)
	{
		xo[i] = xmin + dx_mean * (double)i;
		if (xo[i] > xmax)
		{
			xo[i] = xmax;
		}
	}
	/* Stelle sicher, dass der letzte exakt xmax ist (schöne Abschlussbedingung) */
	xo[m - 1] = xmax;

	/* Durch das Original laufen und linear interpolieren */
	{
		int seg = 0; /* Segmentindex so, dass xo[j] in [x[seg], x[seg+1]] liegt */
		int j;
		for (j = 0; j < m; ++j)
		{
			double xv = xo[j];

			/* seg so verschieben, dass x[seg] <= xv <= x[seg+1] */
			while (seg < n - 2 && xv > x[seg + 1])
			{
				++seg;
			}

			/* Grenzfälle sauber behandeln */
			if (xv <= x[0])
			{
				yo[j] = y[0];
			}
			else if (xv >= x[n - 1])
			{
				yo[j] = y[n - 1];
			}
			else
			{
				/* Jetzt gilt: x[seg] < xv <= x[seg+1] */
				yo[j] = lerp( xv, x[seg], y[seg], x[seg + 1], y[seg + 1] );
			}
		}
	}

	*x_out = xo;
	*y_out = yo;
	*m_out = m;
	return RESAMPLE_OK;
}
/*---------------------------------------------------------------
 *	init_fexpBaseX_log()
*  f(x|a) = a2 * log(a1*x)
 *--------------------------------------------------------------*/
int init_fexpBaseX_log( unsigned int N, double* obs, double* cond,
	double* a, unsigned char* a_flag, FILE* logfile )
{
	double maxCond;
	if (!a_flag[0])
	{
		maxCond = cond[0];
		for (unsigned int i = 0; i < N; i++)
			if (maxCond < cond[i]) maxCond = cond[i];
		a[0] = 1. / maxCond;
	}
	if (!a_flag[1])
	{
		a[1] = 1.1; /* exponent => factor a2 */
	}
	return 0;
}
/*---------------------------------------------------------------
 *	init_fexpBaseX()
 *  f(x|a) = (a1*x)^a2 + a3
 *--------------------------------------------------------------*/
int init_fexpBaseX( unsigned int N, double* obs, double* cond,
	double* a, unsigned char* a_flag, FILE* logfile )
{
	double maxCond, minObs;
	if (!a_flag[0])
	{
		maxCond = cond[0];
		for (unsigned int i = 0; i < N; i++)
			if (maxCond < cond[i]) maxCond = cond[i];
		a[0] = 1. / maxCond;
	}
	if (!a_flag[1])
	{
		a[1] = 1.1; /* exponent */
	}
	if (!a_flag[2])
	{
		minObs = obs[0];
		for (unsigned int i = 0; i < N; i++)
			if (minObs > obs[i]) minObs = obs[i];
		a[2] = minObs;
	}
	return 0;
}

/*---------------------------------------------------------------
 *	init_polynomial()
 * f(x|a) = sum_{j=1}^M aj * x^(j-1)
 *--------------------------------------------------------------*/
int init_polynomial( unsigned int N, double* obs, double* cond,
	double* a, unsigned char* a_flag, FILE* logfile )
{
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
 *	medianInPlace()
 *--------------------------------------------------------------*/
static void medianInPlace( double* obs, unsigned int N )
{
	/* median(a,b,c)=max(min(a,b), min(max(a,b),c)) */
	for (unsigned int i = 1; i < N - 1; i++)
	{
		obs[i] = max( min( obs[i - 1], obs[i] ), min( max( obs[i - 1], obs[i] ), obs[i + 1] ) );
	}
}

/*---------------------------------------------------------------
 *	init_cosine_nonlin()
 * 5: f(x|a) = a1 + a2 * cos( a3*x + a4)
 * 5: f(x|a) = a[0] + a[1] * cos( a[2]*x + a[3])
 *--------------------------------------------------------------*/
int
init_cosine_nonlin( unsigned int N, double* obs, double* cond,
	double* a, unsigned char* a_flag, FILE* logfile )
{
	double sumObs;
	double obsMax, obsMin, obsMaxMid, obsMinMid;
	//double distExtrema;
	double midLow, midHigh;
	double numOfperiods = 0.0;
	unsigned int numCrossings = 0;
	int idxOfMaxObs, idxOfMinObs, idxOfMaxObsMid, idxOfMinObsMid;
	int countMid = 0;
	//int strengthOfSmoothing = 0;
	uint64_t startTime, endTime;

	startTime = __rdtsc();

	/* mean value of observations for a[0] */
	sumObs = obs[0];

	/* estimation of a2 = radius based on range of observations */
	obsMax = obsMin = obs[0];
	idxOfMaxObs = idxOfMinObs = 0;
	obsMaxMid = obsMinMid = obs[N >> 1];
	idxOfMaxObsMid = idxOfMinObsMid = N >> 1;
	midLow = cond[0] + (cond[N - 1] - cond[0]) * 0.333;
	midHigh = cond[0] + (cond[N - 1] - cond[0]) * 0.667;
	for (unsigned int i = 1; i < N; i++)
	{
		sumObs += obs[i]; /* accumulate observation values */
		if (obsMax < obs[i])
		{
			obsMax = obs[i]; idxOfMaxObs = i;
		}
		if (obsMin > obs[i])
		{
			obsMin = obs[i]; idxOfMinObs = i;
		}
		if (cond[i] >= midLow && cond[i] <= midHigh)
		{ /* get extrema in the middle range of the signal, for phase dtermination */
			countMid++;
			if (obsMaxMid < obs[i])
			{
				obsMaxMid = obs[i]; idxOfMaxObsMid = i;
			}
			if (obsMinMid > obs[i])
			{
				obsMinMid = obs[i]; idxOfMinObsMid = i;
			}
		}
	}
	//distExtrema = fabs( cond[idxOfMaxObs] - cond[idxOfMinObs] );
	fprintf( logfile, "\n# global maximum: x,y=%.3f,%.3f", cond[idxOfMaxObs], obs[idxOfMaxObs] );
	fprintf( logfile, "\n# global minimum: x,y=%.3f,%.3f", cond[idxOfMinObs], obs[idxOfMinObs] );
	fprintf( logfile, "\n# middle maximum: x,y=%.3f,%.3f", cond[idxOfMaxObsMid], obs[idxOfMaxObsMid] );
	fprintf( logfile, "\n# middle minimum: x,y=%.3f,%.3f", cond[idxOfMinObsMid], obs[idxOfMinObsMid] );

	if (!a_flag[0])/* if not set on command line */
	{ /* mean is best choice for many samples; */
		/*  few samples can bias one half-wave leading to unfortune parameter */
		//a[0] = (sumObs + obsMax + obsMin) / (N + 2);
		a[0] = sumObs / N;

		//a[0] = (obsMax + obsMin)/ (2); /* better for skewed data */
	}

	if (!a_flag[1])
		a[1] = 0.5 * (obsMax - obsMin);

#ifdef LOMBSCARGLE
	/* filtering of spikes */
	double* obsC = obs;
#else
	double* obsC = vector( N );
	int countSpikes = 0;
	obsC[0] = obs[0];
	for (unsigned int i = 1; i < N - 1; i++)
	{
		if (obs[i - 1] > a[0] && obs[i] < a[0] && obs[i + 1] > a[0])
		{
			/* ol - mean + or -mean > mean - om */
			double minObs = min( obs[i - 1], obs[i + 1] );
			if (minObs - a[0] > a[0] - obs[i])
			{
				obsC[i] = minObs; countSpikes++;
			}
			else obsC[i] = obs[i];
		}
		else if (obs[i - 1] < a[0] && obs[i] > a[0] && obs[i + 1] < a[0])
		{
			/* mean  - ol+ mean -or > om  - mean */
			double maxObs = max( obs[i - 1], obs[i + 1] );
			if (a[0] - maxObs > obs[i] - a[0])
			{
				obsC[i] = maxObs; countSpikes++;
			}
			else obsC[i] = obs[i];
		}
		else obsC[i] = obs[i];
	}
	obsC[N - 1] = obs[N - 1];

	/*################### output for paper ####################*/
	if (1 || obs[0] == 10.9871310751376)
	{
		fprintf( logfile, "\n# cleaned observations (%d spikes) ", countSpikes );
		for (unsigned int i = 0; i < N; i++)
		{
			fprintf( logfile, "\n# %8.5f \t%8.5f", cond[i], obsC[i] );
		}
		fprintf( logfile, "\n# ---------------- " );
	}
	/*################### output for paper end ####################*/
#endif
	/* estimation of a3 = frequency = 1/ period
	 * based on distance between zero crossings
	 */
	if (!a_flag[2])
	{
#ifdef LOMBSCARGLE
		double deltaF, minF, maxF, T, omega;
		double sumCos, sumCos2;
		double sumSin, sumSin2;
		double p, maxP, phi, bestF = -1.;
		fprintf( logfile, "\n# Lomb-Scargle " );

		T = cond[N - 1] - cond[0];
		deltaF = 1. / (5 * T); // min( , (maxF - minF) / 1000 );
		minF = 0.5 / T; /* half of smallest detectable frequency */
		maxF = N / T; /* theoretically: 1/ (2* min(Delta-t)) */
		fprintf( logfile, "\n# f: %.4f : %.4f :%.4f  ", minF, deltaF, maxF );
		maxP = 0.0;
		for (double f = minF; f <= maxF; f += deltaF)
		{
			omega = 2 * M_PI * f;
			/* phi calculation */
			sumCos = sumSin = 0.0;
			for (unsigned int i = 0; i < N; i++)
			{
				sumCos += cos( 2 * omega * cond[i] );
				sumSin += sin( 2 * omega * cond[i] );
			}
			phi = atan( sumSin / sumCos ) * 0.5;
			sumCos = sumSin = 0.0;
			sumCos2 = sumSin2 = 0.0;
			for (unsigned int i = 0; i < N; i++)
			{
				double cosTmp = cos( omega * cond[i] - phi );
				double sinTmp = sin( omega * cond[i] - phi );
				sumCos += cosTmp * (obs[i] - a[0]); /* remove mean value */
				sumSin += sinTmp * (obs[i] - a[0]);
				sumCos2 += cosTmp * cosTmp;
				sumSin2 += sinTmp * sinTmp;
			}
			p = sumCos * sumCos / sumCos2 + sumSin * sumSin / sumSin2;
			//fprintf( logfile, "\n%.4f \t %.4f ", f, p );
			if (maxP < p)
			{
				maxP = p;
				bestF = f;
			}
		}
		a[2] = 2*M_PI*bestF;
		//fprintf( logfile, "\n# best f: %.3f", bestF );
#else
		double reasonableDistance;
#endif

		double sumOfDistances = 0.0, sumObs;
		int above_flag, numObs;
		/* vector of determined distances between mean crosses*/
		double* crossingDistances = NULL;
		/* vector of conditions where curve crosses the mean value */
		double* crossingCond = NULL;
		double* meanObsBetweenX = NULL; /* average observation between two crossings */
		crossingCond = vector( N );
		crossingCond[0] = 0.0;
		crossingDistances = vector( N ); /* worst case: one crossing per sample; also double distances */
		crossingDistances[0] = 0.0;
		meanObsBetweenX = vector( N ); /* worst case: one crossing per sample; also double distances */
		meanObsBetweenX[0] = 0.0;

		// bringt nix medianInPlace( obs, N ); /* smoothing of observations */

		/* count mean-crossings */
		numObs = 1;
		sumObs = obsC[0];
		if (obsC[0] > a[0]) above_flag = 1;
		else above_flag = 0;
		for (unsigned int i = 1; i < N; i++)
		{
			if (above_flag)
			{
				if (obsC[i] < a[0])
				{
					above_flag = 0;
					if (fabs( cond[i] - cond[i - 1] ) < TOL)
					{/* save guard, avoid division by zero */
						crossingCond[numCrossings] = 0.5 * (cond[i] + cond[i - 1]);
					}
					else /* approximate position of crossing */
					{
						double m, n;
						/* |\      |     obs[i-1]*/
						/* | \     |             */
						/* |  \    |             */
						/* |------------- mean --*/
						/* |    \  |             */
						/* |     \ |             */
						/* |      \       obs[i] */
						/* 0   z   cond[i]       */

						/* straight line y = m * x + n */
						/*  y = (obs[i] - obs[i-1]) / (cond[i] - cond[i-1]) * x + n */
						m = (obsC[i] - obsC[i - 1]) / (cond[i] - cond[i - 1]);
						/*  obs[i]   = m * cond[i]   + n   variant 1: n = obs[i]   - m * cond[i] */
						/*  obs[i-1] = m * cond[i-1] + n   variant 2: n = obs[i-1] - m * cond[i-1] */
						n = obsC[i] - m * cond[i];
						/* search for y == mean */
						/*  mean = m * crossCond + n */
						/*  crossCond = (mean - n)/ m  */
						//crossingCond[numCrossings] = (mean - obs[i - 1]) / m + cond[i - 1];
						crossingCond[numCrossings] = (a[0] - n) / m;
					}
					meanObsBetweenX[numCrossings] = sumObs / numObs - a[0];
					numCrossings++;
					sumObs = obsC[i]; numObs=1;
				}
				else
				{
					sumObs += obsC[i]; numObs++;
				}
			}
			else
			{
				if (obsC[i] > a[0])
				{
					above_flag = 1;
					if (fabs( cond[i] - cond[i - 1] ) < TOL)
					{
						crossingCond[numCrossings] = 0.5 * (cond[i] + cond[i - 1]);
					}
					else
					{
						double m, n;
						m = (obsC[i] - obsC[i - 1]) / (cond[i] - cond[i - 1]);
						n = obsC[i] - m * cond[i];
						crossingCond[numCrossings] = (a[0] - n) / m;
					}
					meanObsBetweenX[numCrossings] = a[0] - sumObs / numObs;
					numCrossings++;
					sumObs = obsC[i]; numObs = 1;
				}
				else
				{
					sumObs += obsC[i]; numObs++;
				}
			}
		}
#ifndef LOMBSCARGLE
		free_vector( &obsC );

		fprintf( logfile, "\n# number of detected crossings: %d", numCrossings );
		if (numCrossings > 1)
		{/* enough crossing to estimate distances */
			unsigned long numOfDistances = numCrossings - 1;

			for (unsigned int i = 1; i < numCrossings; i++) /* keep first index unused => compatible to maxObsBetweenX[] */
			{ /* between crossings */
				crossingDistances[i] = crossingCond[i] - crossingCond[i-1];
				sumOfDistances += crossingDistances[i];
			}

			/* check sum of */
			if (numCrossings > 2 && sumOfDistances < 0.333333 * (cond[N - 1] - cond[0]))
			{ /* sanity check */
				fprintf( logfile, "\n# distances:" );
				/* output for debugging */
				for (unsigned int i = 0; i < numOfDistances; i++)
				{ /* between crossings */
					fprintf( logfile, "\n# %3d \t%.6f", i, crossingDistances[i] );
				}
				fprintf( logfile, "\n# sumOfDistances  %.4f", sumOfDistances );
				fprintf( logfile, " is too low (< 0.333 * %.4f)", cond[N - 1] - cond[0] );
				numCrossings = 1; /* assume fraction of a cosine */
				//goto endOfnumCrossings;
			}
			else
			{
				/* sort in ascending order */
				heap_sort_dd( numOfDistances+1, crossingDistances, meanObsBetweenX );

				fprintf( logfile, "\n# distances + max amplitude:" );
				/* output for debugging */
				for (unsigned int i = 1; i <= numOfDistances; i++)
				{ /* between crossings */
					fprintf( logfile, "\n# %3d \t%.6f \t%.6f", i, crossingDistances[i], meanObsBetweenX[i]);
				}

				if (crossingDistances[numOfDistances] < 2. * crossingDistances[1])
				{
					/* no outliers :-) */
					reasonableDistance = sumOfDistances / (numCrossings-1);
					fprintf( logfile, "\n# only good distances" );
				}
				else if (1)
				{ /* find position between clusters */

					/* create histogram, spurious-free bins */
					int goodIdx = 1;
					double binWidth;
					double threshold;
					double referenceDistance;
					double upperLimitDist;
					//double lowerLimitDist;
					double sumOfSpuriousDistances;
					int numOfGoodDistances, flag;// , numOfSpuriousDistances;
					int binIdx, maxCount, maxIdx;// , lastIdx;
					int binCount[8], numCandidateBins, stop;
					int refFirstIdx[8], refLastIdx[8], refBinIdx1 = -1, refBinIdx2;
					unsigned int numOfBins;

					stop = 0;
					numOfBins = 2; /* init */
					do
					{
						for (unsigned int i = 0; i < numOfBins; i++)
						{
							binCount[i] = 0; /* reset */
						}
						binWidth = crossingDistances[numOfDistances] / numOfBins;
						/* slightly more to integrate last distance into last bin*/
						threshold = 1.001 * binWidth;
						binIdx = 0;
						refFirstIdx[binIdx] = 1; /* keep first element of distances array unused */
						fprintf( logfile, "\n# numOfBins: %d, Thresholds:", numOfBins );
						fprintf( logfile, " %f", threshold );
						for (unsigned int i = 1; i <= numOfDistances; i++)
						{
							if (crossingDistances[i] > threshold)
							{
								threshold += binWidth; /* border of next bin */
								binIdx++;/* next bin */
								refFirstIdx[binIdx] = i; /* remember first index in this bin */
								fprintf( logfile, " %f", threshold );
							}
							binCount[binIdx]++;
							refLastIdx[binIdx] = i; /* remember last index in this bin */
						}
						numCandidateBins = 0;
						for (int i = 0; i <= binIdx; i++) /* <=binIdx instead of <numOfBins, since maybe not all bins are used */
						{ /* check whether bin fullfils the requirement */
							if (crossingDistances[refFirstIdx[i]] * 2 > crossingDistances[refLastIdx[i]])
							{
								numCandidateBins++; /* yes */
								if (refBinIdx1 < 0) refBinIdx1 = i; /* if not set yet */
								refBinIdx2 = i; /* update */
							}
						}
						if (numCandidateBins > 1) stop = 1;
						fprintf( logfile, "\t counts:" );
						for (unsigned int i = 0; i < numOfBins; i++)
						{
							fprintf( logfile, " %d,", binCount[i] );
						}
						fprintf( logfile, " stop: %d", stop );
						numOfBins++;
						if (numOfBins > 5) stop = 1;
					} while (!stop);
					if (numOfBins <= 5)
					{
						if (refBinIdx2 == refBinIdx1 + 1)
						{/* adjacent bins */
							if (crossingDistances[refFirstIdx[refBinIdx1]] * 2 > 
								  crossingDistances[refLastIdx[refBinIdx2]])
							{ /* seams to be contiguous, compute the mean */
								double sum = crossingDistances[refFirstIdx[refBinIdx1]];
								for (int i = refFirstIdx[refBinIdx1] + 1; i <= refLastIdx[refBinIdx2]; i++)
								{
									sum += crossingDistances[i];
								}
								referenceDistance = sum / (refLastIdx[refBinIdx2] - refFirstIdx[refBinIdx1] + 1);
								fprintf( logfile, ", contiguous bins, take mean of both" );
							}
							else /* check whether last bin contains long distances */
							{ 
								int num = refLastIdx[refBinIdx2] - refFirstIdx[refBinIdx2] + 1;
								//int refIdx2 = refFirstIdx[refBinIdx2] + (num >> 1);
								int refIdx2 = refLastIdx[refBinIdx2];
								num = refLastIdx[refBinIdx1] - refFirstIdx[refBinIdx1] + 1;
								//int refIdx1 = refFirstIdx[refBinIdx1] + (num >> 1);
								int refIdx1 = refFirstIdx[refBinIdx1];


								double meanObsPerBin[5] = { 0,0,0,0,0 };

								for (int i = 0; i <= binIdx; i++)
								{
									for (int j = refFirstIdx[i]; j <= refLastIdx[i]; j++)
									{
										meanObsPerBin[i] += meanObsBetweenX[j];
									}
									meanObsPerBin[i] /= refLastIdx[i] - refFirstIdx[i] + 1;
								}

								int distanceIsTooLong = crossingDistances[refIdx2] > (2.5 * crossingDistances[refIdx1]);
								int averageObsComparable = meanObsPerBin[refBinIdx2] < 2 * meanObsPerBin[refBinIdx1];
								if (!averageObsComparable)
								{/* only bin2 contains good distances */
									int num = refLastIdx[refBinIdx2] - refFirstIdx[refBinIdx2] + 1;
									int refIdx = refFirstIdx[refBinIdx2] + (num >> 1);
									if (num % 2)
									{
										referenceDistance = crossingDistances[refIdx];
									}
									else
									{
										int refIdx2 = refFirstIdx[refBinIdx2] + ((num - 1) >> 1);
										referenceDistance = 0.5 * (crossingDistances[refIdx] + crossingDistances[refIdx2]);
									}
									fprintf( logfile, "\n# only bin2 has good distances, take median of bin2" );
								}
								else 
								{
									int lastGoodBin;

									if (distanceIsTooLong && (binCount[refBinIdx2] < binCount[refBinIdx1]))
									{ /* last bin contains long distances */
										binIdx--; /* ignore last bin */
									}
									lastGoodBin = binIdx;
									/* look wether there are even more good bins */
									for (int i = binIdx - 1; i >= 0; i--)
									{
										lastGoodBin = i;
										if (i)
										{
											if (meanObsPerBin[i] > 2 * meanObsPerBin[i - 1])
											{
												fprintf( logfile, "\n# amplitude of bin %d too low", i-1 );
												break;
											}
											if (crossingDistances[refLastIdx[i - 1]] > 2 * crossingDistances[refFirstIdx[i - 1]])
											{
												fprintf( logfile, "\n# distances of bin %d too diverse", i - 1 );
												break;
											}
										}
									}
									int num = refLastIdx[binIdx] - refFirstIdx[lastGoodBin] + 1;
									int refIdx = refFirstIdx[lastGoodBin] + (num >> 1);
									if (num % 2)
									{
										referenceDistance = crossingDistances[refIdx];
									}
									else
									{
										int refIdx2 = refFirstIdx[lastGoodBin] + ((num - 1) >> 1);
										referenceDistance = 0.5 * (crossingDistances[refIdx] + crossingDistances[refIdx2]);
									}
									fprintf( logfile, "\n# bins %d to %d contain good distances, take median of all bins", lastGoodBin, refBinIdx2 );
								}

								//{ /* Bin 2 contains possibly no long distances but good ones: take only upper bin */
								//	fprintf( logfile, ", no long distances, take median of upper bin" );
								//	int num = refLastIdx[refBinIdx2] - refFirstIdx[refBinIdx2] + 1;
								//	int refIdx = refFirstIdx[refBinIdx2] + ((refLastIdx[refBinIdx2] - refFirstIdx[refBinIdx2] + 1) >> 1);
								//	if (num % 2)
								//	{
								//		referenceDistance = crossingDistances[refIdx];
								//	}
								//	else
								//	{
								//		int refIdx2 = refFirstIdx[refBinIdx2] + ((refLastIdx[refBinIdx2] - refFirstIdx[refBinIdx2]) >> 1);
								//		referenceDistance = 0.5 * (crossingDistances[refIdx] + crossingDistances[refIdx2]);
								//	}
								//}
								//else /* possibly long distance */
								//{
								//	//int refIdx = refFirstIdx[refBinIdx1] + ((refLastIdx[refBinIdx1] - refFirstIdx[refBinIdx1] + 1) >> 1);
								//	//referenceDistance = crossingDistances[refIdx]; /* take middle+ of lower bins */
								//	//fprintf( logfile, ", take median of middle+ of both bins" );
								//	if (binCount[refBinIdx1] + binCount[refBinIdx2] > 2 || refBinIdx2 == 1 /* only two bins */)
								//	{
								//		int refIdx = refFirstIdx[refBinIdx1] + ((refLastIdx[refBinIdx2] - refFirstIdx[refBinIdx1] + 1) >> 1);
								//		referenceDistance = crossingDistances[refIdx]; /* take middle of both bins */
								//		fprintf( logfile, ", middle+ of both bins" );
								//	}
								//	else
								//	{
								//		int refIdx = refFirstIdx[refBinIdx1];
								//		referenceDistance = crossingDistances[refIdx]; /* take value of lower bin */
								//		fprintf( logfile, ", take value of lower bin" );
								//	}
								//}
							}
						}
						else
						{
							// does this ever can happen???
							// buggy? int refIdx = refFirstIdx[refBinIdx2] + (refLastIdx[refBinIdx2] >> 1);
							int refIdx = refFirstIdx[refBinIdx2] + ((refLastIdx[refBinIdx2] - refFirstIdx[refBinIdx2] + 1) >> 1);
							referenceDistance = crossingDistances[refIdx]; /* take middle of last bin */
							fprintf( logfile, ", middle of last bins" );
							for(int i=0; i< 40;i++)
								fprintf( stderr, "\n ##### middle of last bins ##########" );
						}
					}
					else /* only one good bin */
					{ /* exception handling */
						fprintf( logfile, "\n# exception handling" );
						numOfBins--;
						maxIdx = 1; /* ignore lowest bin assuming that this belongs to spurious distances */
						maxCount = binCount[maxIdx];
						for (unsigned int i = maxIdx + 1; i < numOfBins; i++)
						{
							if (maxCount <= binCount[i])
							{ /* if equal counts, then take higher index */
								maxCount = binCount[i];
								maxIdx = i; /* remember index of bin with maximum count */
							}
						}
						fprintf( logfile, "\t counts:" );
						for (unsigned int i = 0; i < numOfBins; i++)
						{
							fprintf( logfile, " %d,", binCount[i] );
						}
						fprintf( logfile, " maxCount: %d", maxCount );
						int refIdx = refFirstIdx[maxIdx] + ((refLastIdx[maxIdx] - refFirstIdx[maxIdx]) >> 1);
						referenceDistance = crossingDistances[refIdx];
					}

					fprintf( logfile, "\n# referenceDistance: %.4f", referenceDistance );
					//upperLimitDist = 0.55 * referenceDistance;
					//lowerLimitDist = 0.45 * referenceDistance;
					upperLimitDist = 0.5 * referenceDistance;
					//lowerLimitDist = 0.5 * referenceDistance;
					//fprintf( logfile, "\n# lowerLimitDist: %.4f  upperLimitDist: %.4f", lowerLimitDist, upperLimitDist );
					fprintf( logfile, "\n# limitDist: %.4f", upperLimitDist );
					numOfGoodDistances = 0; flag = 1;
//					numOfSpuriousDistances = 0;
					sumOfSpuriousDistances = 0.0;
					for (unsigned int i = 1; i <= numOfDistances; i++)
					{
						if (crossingDistances[i] > upperLimitDist)
						{
							numOfGoodDistances++;
							if (flag) /* only at first occassion */
							{ /* remember index of smallest good distance */
								goodIdx = i; flag = 0;
							}
						}
						else //if (crossingDistances[i] < lowerLimitDist)
						{
//							numOfSpuriousDistances++;
							sumOfSpuriousDistances += crossingDistances[i];
						}
					}
					if (numOfGoodDistances % 2)
					{ /* take median of good distances */
						reasonableDistance = crossingDistances[goodIdx + (numOfGoodDistances >> 1)];
					}
					else
					{/* take average */
						reasonableDistance = 0.5 * (crossingDistances[goodIdx + (numOfGoodDistances >> 1)] +
							crossingDistances[goodIdx + (numOfGoodDistances >> 1) - 1]);
					}
					fprintf( logfile, "\n# reasonableDistance: %.5f", reasonableDistance );
					fprintf( logfile, "\n# numOfGoodDistances: %d ", numOfGoodDistances );
//					fprintf( logfile, "\t numOfSpuriousDistances: %d ", numOfSpuriousDistances );
					fprintf( logfile, "\n# sumOfSpuriousDistances: %f ", sumOfSpuriousDistances );
					fprintf( logfile, "\t sumOfDistances: %f ", sumOfDistances );

					{ /* number of distances in rejection class*/
//						int numReject = numOfDistances - numOfGoodDistances - numOfSpuriousDistances;
						//double correctionTerm = min( reasonableDistance,
						//	sumOfSpuriousDistances / (numOfGoodDistances + (numReject >> 1)) ); /* limit correction to 100% */
						double correctionTerm = min( reasonableDistance,
							sumOfSpuriousDistances / numOfGoodDistances ); /* limit correction to 100% */
						/* good distances can be reduced to 50%, so adding same amount (100%) restores original length */
						reasonableDistance += correctionTerm;
//						fprintf( logfile, "\n# number of ignored distances: %d ", numReject );
						fprintf( logfile, "\n# correctionTerm: %.4f", correctionTerm );
					}
				}
				else
				{
					double maxMeanObs, aveMeanObs, sumOfSpuriousDistances;
					//double reasonableDistance;
					double correctionTerm;
					unsigned int numOfGoodDistances, numOfSpuriousDistances, idx, aveCount;

					numOfSpuriousDistances = 0; /* default initialization */
					aveMeanObs = maxMeanObs = meanObsBetweenX[numOfDistances];
					aveCount = 1;
					for (int i = numOfDistances - 1; i > 0; i--)
					{
						aveMeanObs = (aveMeanObs * aveCount + meanObsBetweenX[i]) / (aveCount + 1);
						aveCount++;
						//if (maxMeanObs < meanObsBetweenX[i]) maxMeanObs = meanObsBetweenX[i];
						//else 
						if (aveMeanObs > 3 * meanObsBetweenX[i])
						{
							numOfSpuriousDistances = i;
							break;
						}
					}
					fprintf( logfile, "\n# aveMeanObs: %.5f", aveMeanObs );
					numOfGoodDistances = numOfDistances - numOfSpuriousDistances;
					idx = numOfSpuriousDistances + ((numOfGoodDistances + 1) >> 1);
					reasonableDistance = crossingDistances[idx+1]; /* index = 0 is not used */
					sumOfSpuriousDistances = 0.0;
					for (int i = 1; i <= (signed)numOfDistances - (signed)numOfGoodDistances; i++)
					{
						sumOfSpuriousDistances += crossingDistances[i];
					}
					//correctionTerm = min( reasonableDistance,sumOfSpuriousDistances / numOfGoodDistances);
					correctionTerm = sumOfSpuriousDistances / numOfGoodDistances;
					fprintf( logfile, "\n# reasonableDistance: %.5f", reasonableDistance );
					fprintf( logfile, "\n# numOfGoodDistances: %d ", numOfGoodDistances );
					fprintf( logfile, "\n# sumOfSpuriousDistances: %f ", sumOfSpuriousDistances );
					//fprintf( logfile, "\t sumOfDistances: %f ", sumOfDistances );
					fprintf( logfile, "\n# correctionTerm: %.4f", correctionTerm );
					reasonableDistance += correctionTerm;
				}
				/* compare with distance between yMax and yMin */
				fprintf( logfile, "\n# corrected reasonableDistance: %.4f", reasonableDistance );
				/* determine max und min positions in left and right half of signal*/
				if (0) 
				{ /* is only reliable when sampling rate is high enough */
					double distMinMax, distMaxMin, distExtremeOdd;
					double distMaxMax, distMinMin, distExtremeEven;
					double relDoubleOdd, relDoubleEven;

					double obsMaxLeft, obsMinLeft;
					double obsMaxRight, obsMinRight;
					double midCond;
					int idxOfMaxObsLeft, idxOfMinObsLeft;
					int idxOfMaxObsRight, idxOfMinObsRight;
					int relIntOdd, relIntEven;

					obsMaxLeft = obsMinLeft = obs[0];
					idxOfMaxObsLeft = idxOfMinObsLeft = 0;
					obsMaxRight = obsMinRight = obs[N-1];
					idxOfMaxObsRight = idxOfMinObsRight = N-1;
					midCond= cond[0] + (cond[N - 1] - cond[0]) * 0.5;
					for (unsigned int i = 1; i < N-1; i++)
					{
						if (cond[i] < midHigh)
						{
							if (obsMaxLeft < obs[i])
							{
								obsMaxLeft = obs[i]; idxOfMaxObsLeft = i;
							}
							if (obsMinLeft > obs[i])
							{
								obsMinLeft = obs[i]; idxOfMinObsLeft = i;
							}
						}
						else
						{
							if (obsMaxRight < obs[i])
							{
								obsMaxRight = obs[i]; idxOfMaxObsRight = i;
							}
							if (obsMinRight > obs[i])
							{
								obsMinRight = obs[i]; idxOfMinObsRight = i;
							}
						}
					}
					distMaxMin = fabs( cond[idxOfMaxObsRight] - cond[idxOfMinObsLeft] );
					distMinMax = fabs( cond[idxOfMinObsRight] - cond[idxOfMaxObsLeft] );
					distExtremeOdd = max( distMaxMin, distMinMax );
					distMaxMax = fabs( cond[idxOfMaxObsRight] - cond[idxOfMaxObsLeft] );
					distMinMin = fabs( cond[idxOfMinObsRight] - cond[idxOfMinObsLeft] );
					distExtremeEven = max( distMaxMax, distMinMin );

					relDoubleOdd = distExtremeOdd / reasonableDistance;
					relDoubleEven = distExtremeEven / reasonableDistance;
					relIntOdd = (int)floor( relDoubleOdd + 0.5 );
					relIntEven = (int)floor( relDoubleEven + 0.5 );
					fprintf( logfile, "\n# distance of MaxMin: %.3f =  %d times %.4f", distExtremeOdd, relIntOdd, distExtremeOdd / relIntOdd );
					fprintf( logfile, "\n# distance of MM    : %.3f =  %d times %.4f", distExtremeEven, relIntEven, distExtremeEven / relIntEven );
					if (relIntOdd % 2 && relIntOdd > 2)
					{
						if (relIntEven % 2 == 0 && relIntEven > 1)
						{
							reasonableDistance = 0.33333 * (reasonableDistance + distExtremeOdd / relIntOdd + distExtremeEven / relIntEven);
						}
						else
						{
							/* odd multiple, seams to be ok */
							reasonableDistance = 0.5 * (reasonableDistance + distExtremeOdd / relIntOdd);
						}
					}
					else if (relIntEven % 2 == 0 && relIntEven > 1)
					{
						/* even multiple, seams to be ok */
						reasonableDistance = 0.5 * (reasonableDistance + distExtremeEven / relIntEven);
					}
					fprintf( logfile, "\n# final reasonableDistance: %.4f", reasonableDistance );
				}

				a[2] = M_PI / reasonableDistance; /* is time period * 2 * Pi; one period == 2 * medianDistance  */
			}
		}
//endOfnumCrossings:
		if (numCrossings <= 1) /* too less crossings for  proper estimation */
		{
			/* get range of conditions (approx. time period) */
			//a[2] = 0.5 * M_PI / (cond[N-1] - cond[0]); /* presume that conditions cover a quarter oscillation  */
			a[2] = M_PI / (cond[N - 1] - cond[0]); /* presume that conditions cover a half oscillation  */
			fprintf( logfile, "\n# cond[N-1] - cond[0]: %.4f", cond[N - 1] - cond[0] );
		}
#endif // not LOMBSCARGLE
		numOfperiods = a[2] / (2 * M_PI) * (cond[N - 1] - cond[0]); /* P = f * T */
		fprintf( logfile, "\n# number of periods: %.3f", numOfperiods );
		fprintf( logfile, "\n# initial frequency: %.3f", a[2] / (2 * M_PI) );
		free_vector( &crossingDistances );
		free_vector( &crossingCond );
	}

	/* estimation of a4 = phase shift */
	if (!a_flag[3])
	{
		//if (1 || numCrossings > 1) /* range of time is large enough */
		{/* take position of highest extremum for phase shift */
			/* assume  (a3*x+a4) = 0 for maximal amplitude */
			if (countMid < 10 || numOfperiods < 3) /* rough threshold */
			{ /* use global extremum, when only few samples are taken */
				fprintf( logfile, "\n# use global extremum (%d samples in the middle, %.2f periods", countMid , numOfperiods );
				if ((obsMax - a[0]) > fabs( a[0] - obsMin ))
				{
					a[3] = -a[2] * cond[idxOfMaxObs];
				}
				else
				{
					a[3] = -a[2] * cond[idxOfMinObs] + M_PI;
				}
			}
			else
			{ /* use extremum in the middle of the sampled signal */
				fprintf( logfile, "\n# use middle extremum" );
				if ((obsMaxMid - a[0]) > fabs( a[0] - obsMinMid ))
				{
					a[3] = -a[2] * cond[idxOfMaxObsMid];
					fprintf( logfile, "\t maximum=%f at x=%f", obsMaxMid, cond[idxOfMaxObsMid] );
					fprintf( logfile, "\t a4=%f ", a[3] );
				}
				else
				{
					a[3] = -a[2] * cond[idxOfMinObsMid] + M_PI;
					fprintf( logfile, "\t minimum=%f at x=%f", obsMinMid, cond[idxOfMinObsMid] );
					fprintf( logfile, "\t a4=%f ", a[3] );
				}
			}
		}

		while (a[3] < 0.) a[3] += 2 * M_PI;
		while (a[3] > 2 * M_PI) a[3] -= 2 * M_PI;
		//fprintf( logfile, "\t a4=%f ", a[3] );
	}

	endTime = __rdtsc();
	fprintf( logfile, "\n# number of cycles: %lld",	endTime - startTime );

	return 0;
}

/*---------------------------------------------------------------
 *	init_cosine_trend()
 * 12: f(x|a) = a1 + a2 * x + a3 * cos( x - a4)
 *--------------------------------------------------------------*/
int
init_cosine_trend( unsigned int N, double* obs, double* cond,
	double* a, unsigned char* a_flag, FILE* logfile )
{
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
 *	init_trigonometric1()
 * f(x|a) = a1 + a2*cos(a3*x-a4)
 *--------------------------------------------------------------*/
int
init_trigonometric1( unsigned int N, double* obs, double* cond,
	double* a, unsigned char* a_flag, FILE* logfile )
{
	double mean;

	/* mean value for a[0] */
	if (!a_flag[0])
	{
		mean = 0;
		for (unsigned int i = 0; i < N; i++)
			mean += obs[i];
		a[0] = mean / N;
	}

	/* estimation of a3 = period */
	{
		double max_x, min_x;
		max_x = min_x = cond[0];
		for (unsigned int i = 1; i < N; i++)
		{
			if (max_x < cond[i]) max_x = cond[i];
			if (min_x > cond[i]) min_x = cond[i];
		}
		if (!a_flag[2])
			a[2] = 2 * 3.141 / (2 * (max_x - min_x));
	}

	/* estimation of a2 = amplitude */
	{
		double max_obs, min_obs;
		max_obs = min_obs = obs[0];
		for (unsigned int i = 1; i < N; i++)
		{
			if (max_obs < obs[i]) max_obs = obs[i];
			if (min_obs > obs[i]) min_obs = obs[i];
		}
		if (!a_flag[1]) /* if not set on command line */
		{
			a[1] = 0.5 * (max_obs - min_obs);
		}
	}

	/* estimation of a4 = phase shift */
	if (!a_flag[3])
	{
		a[3] = 0.;	/* dummy */
	}

	return 0;
}

/*---------------------------------------------------------------
 *	init_trigonometric2()
 * f(x|a) = a1 + a2*cos(a3*x-a4) + a5*cos(2*a3*x-a6)
 *--------------------------------------------------------------*/
int
init_trigonometric2( unsigned int N, double* obs, double* cond,
	double* a, unsigned char* a_flag, FILE* logfile )
{
	double mean;

	/* mean value for a[0] */
	if (!a_flag[0])
	{
		mean = 0;
		for (unsigned int i = 0; i < N; i++)
			mean += obs[i];
		a[0] = mean / N;
	}

	/* estimation of a3 = period */
	{
		double max_x, min_x;
		max_x = min_x = cond[0];
		for (unsigned int i = 1; i < N; i++)
		{
			if (max_x < cond[i]) max_x = cond[i];
			if (min_x > cond[i]) min_x = cond[i];
		}
		if (!a_flag[2])
			a[2] = 2 * 3.141 / (2 * (max_x - min_x));
	}

	/* estimation of a2, a5  = amplitude */
	{
		double max_obs, min_obs;
		max_obs = min_obs = obs[0];
		for (unsigned int i = 1; i < N; i++)
		{
			if (max_obs < obs[i]) max_obs = obs[i];
			if (min_obs > obs[i]) min_obs = obs[i];
		}
		if (!a_flag[1]) /* if not set on command line */
		{
			a[1] = 0.5 * (max_obs - min_obs);
		}
		if (!a_flag[4]) /* if not set on command line */
		{
			a[4] = rand() * 0.5 * (max_obs - min_obs) / RAND_MAX;
		}
	}

	/* estimation of a4,a6 = phase shift */
	if (!a_flag[3])
	{
		a[3] = 0.;	/* dummy */
	}
	if (!a_flag[5])
	{
		a[5] = 0.;	/* dummy */
	}
	return 0;
}

/*---------------------------------------------------------------
 *	init_logarithmic()
 * f(x|a) = a1 * log( x + a2)
 *--------------------------------------------------------------*/
int
init_logarithmic( unsigned int N, double* obs, double* cond,
	double* a, unsigned char* a_flag, FILE* logfile )
{
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
	if (!a_flag[1])
	{
		a[1] = -condCloseToZero;
		if (cond[0] + a[1] <= 0)	a[1] = -cond[0] + 10e-5;

		/* 0= a1*log(x+a2) ==> x+a2 >= 0 horizontal shift */
	}
	//if (obsMaxIdx >= obsMinIdx)
	if (!a_flag[0])
	{
		if (obsMax >= 0. && obsMin < 0.)
		{ /* graph crosses the x-axis, estimate a[1] should be reliable enough */
			double sum = 0.0;
			/* average all possible values of a1*/
			for (unsigned int i = 0; i < N; i++)
			{
				sum += obs[i] / log( cond[i] + a[1] );
			}
			a[0] = sum / N;
		}
		else /* for variants are possible */
		{
			int variant = 0;
			fprintf( logfile, "\n# graph does not cross the x-axis:" );
			if (obsMax < 0.0)
			{
				/*  -----+----- > x
				 *   ____
				 *  /
				 * /
				 *           */
				if (obsMaxIdx > obsMinIdx) 	variant = 1;
				/* +----------- > x
				 *   \
				 *    \
				 *     ----
				 *           */
				else 	variant = 2;
			}
			else /* obsMin is > 0. */
			{
				/*     ____
				 *    /
				 *   /
				 * +---------- > x          */
				if (obsMaxIdx > obsMinIdx) 	variant = 3;
				/*
				 *  \
				 *   \----
				 * ---------+- > x          */
				else 	variant = 4;
			}
			if (variant == 2 || variant == 3)
				/* use slope of function, presume that x >> a2 */
			{
				double slope = (obs[N - 1] - obs[N >> 1]) / (cond[N - 1] - cond[N >> 1]);
				a[0] = slope * cond[N - 1];
				fprintf( logfile, "\n use slope * x_N = %.2f * %.2f", slope, cond[N - 1] );
			}
			else /* cannot presume that x >> a2 */
			{
				/* super coarse estimation  */
				if (variant == 1) a[0] = 1;
				else a[0] = -1;
			}
		}
	}
	return 0;
}

/*---------------------------------------------------------------
 *	init_exponentialdummy()
 *  f(x|a) = a1 + a2 * exp( a3 * x)
 *--------------------------------------------------------------*/
int init_exponentialdummy( unsigned int N, double* obs, double* cond,
	double* a, unsigned char* a_flag, FILE* logfile )
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
		fprintf( logfile, "\n# decreasing function" );
		if (fabs( obsLeft - obsMid ) > fabs( obsRight - obsMid ))
		{
			variant = 1; /* A: °.. */
			fprintf( logfile, ", negative a3, positive scaling factor a2" );
			a[1] = +1; a[2] = -0.5; /* dummy values */
		}
		else
		{
			variant = 3; /* C: °°. */
			fprintf( logfile, ", positive a3, negative scaling factor a2" );
			a[1] = -1; a[2] = 0.5; /* dummy values */
		}
	}
	else /* increasing function */
	{
		fprintf( logfile, "\n# increasing function" );
		if (fabs( obsLeft - obsMid ) > fabs( obsRight - obsMid ))
		{
			variant = 4; /* D: .°° */
			fprintf( logfile, ", negative a3, negative scaling factor a2" );
			a[1] = -1; a[2] = -0.5; /* dummy values */
		}
		else
		{
			variant = 2; /* B: ..° */
			fprintf( logfile, ", positive a3, positive scaling factor a2" );
			a[1] = +1; a[2] = 0.5; /* dummy values */
		}
	}
	/*------ estimating offset a1 ------------------*/
	/* assume that curve approaches asymptotically a1 */
	if (variant == 1 || variant == 2)		a[0] = obsMin;
	else		a[0] = obsMax;
	return 0;
}

/*---------------------------------------------------------------
	*	init_exponentialLinearised()
	*  f(x|a) = a1 + a2 * exp( a3 * x)
	*--------------------------------------------------------------*/
int init_exponentialLinearised( unsigned int N, double* obs, double* cond,
	double* a, unsigned char* a_flag, FILE* logfile )
{
	int err;
	int variant = 0; /* curvature of graph 1..4  = A..D */
	unsigned int idxCloseToMiddle, M = 2;
	double bestDiff, b[2];
	double thresh, smax;
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
		fprintf( logfile, "\n# decreasing function" );
		if (fabs( obsLeft - obsMid ) > fabs( obsRight - obsMid ))
		{
			variant = 1; /* A: °.. */
			fprintf( logfile, ", negative a3, positive scaling factor a2" );
			a[1] = +1; a[2] = -0.5; /* dummy values */
		}
		else
		{
			variant = 3; /* C: °°. */
			fprintf( logfile, ", positive a3, negative scaling factor a2" );
			a[1] = -1; a[2] = 0.5; /* dummy values */
		}
	}
	else /* increasing function */
	{
		fprintf( logfile, "\n# increasing function" );
		if (fabs( obsLeft - obsMid ) > fabs( obsRight - obsMid ))
		{
			variant = 4; /* D: .°° */
			fprintf( logfile, ", negative a3, negative scaling factor a2" );
			a[1] = -1; a[2] = -0.5; /* dummy values */
		}
		else
		{
			variant = 2; /* B: ..° */
			fprintf( logfile, ", positive a3, positive scaling factor a2" );
			a[1] = +1; a[2] = 0.5; /* dummy values */
		}
	}

	/* assume that curve approaches asymptotically a1 */
	if (variant == 1 || variant == 2)		a[0] = obsMin;
	else		a[0] = obsMax;
	fprintf( logfile, "\n# linear regression for initial parameters  #####" );

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
	smax = 0.0;
	for (unsigned int j = 0; j < M; j++)
	{
		if (s[j] > smax)	smax = s[j];
	}
	if (smax < TOL_S)
	{
		fprintf( stderr,
			"\n###\n###     singular matrix, smax = %f", smax );
		fprintf( logfile,
			"\n###\n###     singular matrix, smax = %f", smax );

		err = 1;
		goto endfunc; /* use dummy values */
	}
	else if (smax > 1.e+31)
	{
		fprintf( stderr,
			"\n###\n###     degraded matrix, smax = huge" );
		fprintf( logfile,
			"\n###\n###     degraded matrix, smax = huge" );
		err = 1;
		goto endfunc; /* use dummy values */
	}

	thresh = MIN( TOL_S * smax, TOL_S );

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
	a[2] = b[1];
	a[1] = exp( b[0] );
	if (variant > 2)		a[1] = -a[1];

endfunc:
	free_vector( &lnObs );
	free_matrix( &jacob );
	free_vector( &s );
	free_matrix( &V );
	free_matrix( &VinvS );
	free_matrix( &invJ );
	fprintf( logfile, "\n# a1= %8.4f  a2= %8.4f  a3= %8.4f     #####", a[0], a[1], a[2] );
	return err;
}
/*---------------------------------------------------------------
 *	init_exponential()
 *  f(x|a) = a1 + a2 * exp( a3 * x)
 *--------------------------------------------------------------*/
int
init_exponential( unsigned int N, double* obs, double* cond,
	double* a, unsigned char* a_flag, FILE* logfile )
{
	double slope = 0., sign; /* dummy defaut */
	double bestDiff;
	unsigned int idx1CloseToZero, idx2CloseToZero, idxCloseToMiddle;
	double condCloseToZero1;
	double obsMin, obsMax, obsMid, condMid;
	double obsLeft, obsRight, obsRef, condRef;
	double b[3]; /* parameters based on alternative approach */
	int variant = 0; /* curvature of graph 1..4  = A..D */
	int aParamsOk_flag = 1, bParamsOk_flag = 1;

	/*-------- derive some needed values ----------------*/
	idx1CloseToZero = 0;
	condCloseToZero1 = fabs( cond[idx1CloseToZero] );
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
		if (condCloseToZero1 > fabs( cond[i] ))
		{
			condCloseToZero1 = fabs( cond[i] );
			idx1CloseToZero = i;
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
		fprintf( logfile, "\n# decreasing function" );
		sign = -1;
		if (fabs( obsLeft - obsMid ) > fabs( obsRight - obsMid ))
		{
			variant = 1; /* A: °.. */
			obsRef = obsLeft;
			condRef = cond[0];
			fprintf( logfile, ", negative a3, positive scaling factor a2" );
		}
		else
		{
			variant = 3; /* C: °°. */
			obsRef = obsRight;
			condRef = cond[N - 1];
			fprintf( logfile, ", positive a3, negative scaling factor a2" );
		}
	}
	else /* increasing function */
	{
		fprintf( logfile, "\n# increasing function" );
		sign = +1;
		if (fabs( obsLeft - obsMid ) > fabs( obsRight - obsMid ))
		{
			variant = 4; /* D: .°° */
			obsRef = obsLeft;
			condRef = cond[0];
			fprintf( logfile, ", negative a3, negative scaling factor a2" );
		}
		else
		{
			variant = 2; /* B: ..° */
			obsRef = obsRight;
			condRef = cond[N - 1];
			fprintf( logfile, ", positive a3, positive scaling factor a2" );
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
	/* it is not recommended to use a[0] = obs[0] oder obs[N-1], because
		* a[1] depends on a[0] and its calculation could lead to wrong
		* values of a[1] (wrong sign)  */
	b[0] = a[0]; /* copy for second approach */

	/*-------- estimating scale factor a2 ------------------*/
	/* find second point with condition close to zero for slope computation */
	if (cond[idx1CloseToZero] > 0.)
	{
		if (idx1CloseToZero > 0)
		{ /* there is a point on the opposite side */
			idx2CloseToZero = idx1CloseToZero;
			idx1CloseToZero--;
		}
		else
		{ /* we maust take a point on the same side */
			idx2CloseToZero = idx1CloseToZero + 1;
		}
	}
	else
	{
		if (idx1CloseToZero < N - 1)
		{ /* there is a point on the opposite side */
			idx2CloseToZero = idx1CloseToZero + 1;
		}
		else
		{ /* we must take a point on the same side */
			idx2CloseToZero = idx1CloseToZero;
			idx1CloseToZero--;
		}
	}

	//if (idx1CloseToZero > idx2CloseToZero)
	//{ /* swap indices if necessary */
	//	int tmp = idx1CloseToZero;
	//	idx1CloseToZero = idx2CloseToZero;
	//	idx2CloseToZero = tmp;
	//}

	while (fabs( cond[idx1CloseToZero] - cond[idx2CloseToZero] ) < eTOL &&
		idx2CloseToZero < N - 1 && idx1CloseToZero > 0)
	{ /* too close */
		if (idx2CloseToZero < N - 1) idx2CloseToZero++;
		else if (idx1CloseToZero > 0)idx1CloseToZero--;
		else break;
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
			fprintf( logfile, "\n#### suspicious behavior of function                             #####" );
			fprintf( logfile, "\n#### cannot find a point pair that satisfies the needed slope sign #####" );
			fprintf( stdout, "\n#### suspicious behavior of function                             #####" );
			fprintf( stdout, "\n#### cannot find a point pair that satisfies the needed slope sign #####" );
			break;
		}
	}

	fprintf( logfile, "\n# slope: %0.3f", slope );

	if (!a_flag[1])
	{
		a[1] = obs[idx1CloseToZero] - slope * cond[idx1CloseToZero];
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

		b[1] = (obsMid - b[0]) / exp( b[2] * condMid );
	}
	else /* (a[0] == obsMax) */
	{
		b[2] = (log( b[0] - obsMid ) - log( b[0] - obsRef ))
			/ (condMid - condRef);

		b[1] = -(obsMid - b[0]) / exp( b[2] * condMid );
	}
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

	fprintf( logfile, "\n# a1= %8.4f  a2= %8.4f  a3= %8.4f     #####", a[0], a[1], a[2] );
	fprintf( logfile, "\n# b1= %8.4f  b2= %8.4f  b3= %8.4f     #####", b[0], b[1], b[2] );

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
		fprintf( logfile, "\n#\n# estimation of initial parameters failed ! #####" );
		fprintf( logfile, "\n# use dummy initial parameters ! #####" );
		/* keep a1 */;
		if (variant == 1) { if (!a_flag[1]) a[1] = 1; if (!a_flag[2]) a[2] = -0.5; }
		else if (variant == 2) { if (!a_flag[1]) a[1] = 1; if (!a_flag[2]) a[2] = 0.5; }
		else if (variant == 3) { if (!a_flag[1]) a[1] = -1; if (!a_flag[2]) a[2] = 0.5; }
		else if (variant == 4) { if (!a_flag[1]) a[1] = -1; if (!a_flag[2]) a[2] = -0.5; }
	}
	fprintf( logfile, "\n# a1= %8.4f  a2= %8.4f  a3= %8.4f     #####", a[0], a[1], a[2] );

	return 0;
}

/*---------------------------------------------------------------
 *	init_expon2()
 * f(x|a) = a1 * exp( a2 * x)
 *--------------------------------------------------------------*/
int
init_expon2( unsigned int N, double* obs, double* cond,
	double* a, unsigned char* a_flag, FILE* out )
{
	double slope1 = 0, slope2 = 0; /* dummy defaut */
	double condCloseToZero;
	unsigned int idx1CloseToZero;
	int stop = 0;

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
	//if (obs[0] > obs[N - 1]) /* decaying process */
	//{/* find second position for slope calculation */
	//	int iStep, jStep;
	//	if (idx1CloseToZero > (N >> 1))
	//	{
	//		iStep = 0;
	//		jStep = -1;  /* look to the left */
	//	}
	//	else
	//	{
	//		iStep = +1; /* look to the right */
	//		jStep = 0;
	//	}
	//	i = idx1CloseToZero + iStep; /* right neighbour */
	//	j = idx1CloseToZero + jStep; /* left neighbour */
	//	do
	//	{
	//		if (obs[j] > obs[i] && fabs(cond[j] - cond[i]) > 10e-5)
	//		{
	//			slope = (obs[j] - obs[i]) / (cond[j] - cond[i]);
	//			stop = 1;
	//		}
	//		i += iStep;
	//		j += jStep;
	//	} while (!stop || i == N);
	//}
	//else /* increasing process */
	//{/* find second position for slope calculation */
	//	int iStep, jStep;
	//	if (idx1CloseToZero > (N >> 1))
	//	{
	//		iStep = 0;
	//		jStep = -1;  /* look to the left */
	//	}
	//	else
	//	{
	//		iStep = +1; /* look to the right */
	//		jStep = 0;
	//	}
	//	i = idx1CloseToZero + iStep; /* right neighbour */
	//	j = idx1CloseToZero + jStep; /* left neighbour */
	//	do
	//	{
	//		if (obs[j] < obs[i] && fabs( cond[j] - cond[i] ) > 10e-5)
	//		{
	//			slope = (obs[j] - obs[i]) / (cond[j] - cond[i]);
	//			stop = 1;
	//		}
	//		i += iStep;
	//		j += jStep;
	//	} while (!stop || i == N);
	//}
	/* estimation of a1 = head of function */

	/* estimation of a2 = gradient at head of function
	 * a2 = f'(0)/a1
	 */
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
int
init_gauss( unsigned int N, double* obs, double* cond,
	double* a, unsigned char* a_flag, FILE* out )
{
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
		err = 8;
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
int
init_gen_laplace( unsigned int N, double* obs, double* cond,
	double* a, unsigned char* a_flag, FILE* out )
{
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
 *	init_circlelin()
 *  f(x|a) = 0 = (x1-a1)^2 + (x2-a2)^2 - a3^2
 *--------------------------------------------------------------*/
int
init_circlelin( unsigned  int N, double* obs, double* cond,
	double* a, unsigned char* a_flag, FILE* logfile )
{
	int err = 0;		 /* return value */
	double b1, b2, b3;

	/* get estimates of centre coordinates and radius */
	err = init_circle( N, obs, cond, a, a_flag, logfile );

	/* convert into vector b	*/
	b1 = 2 * a[0];
	b2 = 2 * a[1];
	b3 = a[0] * a[0] + a[1] * a[1] - a[2] * a[2];

	/* put back to a[]	*/
	a[0] = b1;
	a[1] = b2;
	a[2] = b3;

	return err;
}

/*---------------------------------------------------------------
 *	init_circle()
 *  f(x|a) = 0 = (x1-a1)^2 + (x2-a2)^2 - a3^2
 *--------------------------------------------------------------*/
int
init_circle( unsigned int N, double* obs, double* cond,
	double* a, unsigned char* a_flag, FILE* logfile )
{
	int err = 0;		 /* return value */
	double sum_x, sum_y, rad2, diff1, diff2;

	fprintf( logfile, "\n#\n# init_circle()" );

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
	fprintf( logfile, "\n#\n# mean of condition coordinates" );
	fprintf( logfile, "\n#   mean(x)= %f", sum_x );
	fprintf( logfile, "\n#   mean(y)= %f", sum_y );
	fprintf( logfile, "\n#   radius = %f", rad2 );

	if (!a_flag[0])	a[0] = sum_x;
	if (!a_flag[1])	a[1] = sum_y;
	if (!a_flag[2])	a[2] = rad2;

	fprintf( logfile,
		"\n#  f(x|a) =0= (x1-%f)**2 + (x2-%f)**2 - %f**2",
		a[0], a[1], a[2] );

	return err;
}

/*---------------------------------------------------------------
 *	init_rotation()
 * 21... f1(x|a) = a1 + cos(a3) * x1 - sin(a3) * x2
 *       f2(x|a) = a2 + sin(a3) * x1 + cos(a3) * x2
 *-------------------------------------------------------------*/
int
init_rotation( unsigned int N, double* obs, double* cond,
	double* a, unsigned char* a_flag, FILE* logfile )
{
	int err = 0;		 /* return value */
	double sum_x, sum_y, sum_u, sum_v;

	fprintf( logfile, "\n#\n# init_rotation()" );

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
	fprintf( logfile, "\n#\n# mean of condition coordinates" );
	fprintf( logfile, "\n#   mean(u)= %f", sum_u );
	fprintf( logfile, "\n#   mean(v)= %f", sum_v );
	fprintf( logfile, "\n# mean of observed coordinates" );
	fprintf( logfile, "\n#   mean(x)= %f", sum_x );
	fprintf( logfile, "\n#   mean(y)= %f", sum_y );

	if (!a_flag[0])	a[0] = sum_x - sum_u;
	if (!a_flag[1])	a[1] = sum_y - sum_v;
	if (!a_flag[2])	a[2] = 0;		/* assume no rotation */


	fprintf( logfile,
		"\n#  f1(u,v) = %f + cos(%f) * u - sin(%f) * v",
		a[0], a[2], a[2] );
	fprintf( logfile,
		"\n#  f2(u,v) = %f + sin(%f) * u + cos(%f) * v",
		a[1], a[2], a[2] );

	return err;
}

/*---------------------------------------------------------------
 *	init_NN3x3x1()
 * 3x3x1
 *--------------------------------------------------------------*/
int
init_NN3x3x1( unsigned int N, double* obs, double* cond,
	double* a, unsigned char* a_flag, FILE* logfile )
{
	int j;
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
	for (j = 0; j < M_MAX; j++)
	{
		if (!a_flag[j])
			a[j] = 2. * (float)random() / (float)RAND_MAX - 1.;
	}

	/* make random numbers in a range that |cond x param| < 5 */
	minval = maxval = cond[0];
	for (unsigned int i = 1; i < N; i++)
	{
		if (minval > cond[3 * i]) minval = cond[3 * i];
		if (maxval < cond[3 * i]) maxval = cond[3 * i];
	}
	/* weights from 1st input */
	if (!a_flag[1])	a[1] = a[1] / (maxval - minval);
	if (!a_flag[5])	a[5] = a[5] / (maxval - minval);
	if (!a_flag[9]) a[9] = a[9] / (maxval - minval);

	minval = maxval = cond[1];
	for (unsigned int i = 1; i < N; i++)
	{
		if (minval > cond[3 * i + 1]) minval = cond[3 * i + 1];
		if (maxval < cond[3 * i + 1]) maxval = cond[3 * i + 1];
	}
	/* weights from 2nd input */
	if (!a_flag[2])	a[2] = a[2] / (maxval - minval);
	if (!a_flag[6])	a[6] = a[6] / (maxval - minval);
	if (!a_flag[10])a[10] = a[10] / (maxval - minval);

	minval = maxval = cond[2];
	for (unsigned int i = 1; i < N; i++)
	{
		if (minval > cond[3 * i + 2]) minval = cond[3 * i + 2];
		if (maxval < cond[3 * i + 2]) maxval = cond[3 * i + 2];
	}
	/* weights from 3rd input */
	if (!a_flag[3])	a[3] = a[3] / (maxval - minval);
	if (!a_flag[7])	a[7] = a[7] / (maxval - minval);
	if (!a_flag[11])a[11] = a[11] / (maxval - minval);

	return 0;
}

/*---------------------------------------------------------------
 *	init_NN1x3x1()
 * 1x3x1
 *--------------------------------------------------------------*/
int
init_NN1x3x1( unsigned int N, double* obs, double* cond,
	double* a, unsigned char* a_flag, FILE* logfile )
{
	int  j;
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
	for (j = 0; j < M_MAX; j++)
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
int
init_NN( unsigned int N, double* obs, double* cond,
	double* a, unsigned char* a_flag, FILE* logfile )
{
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
 *  for variable model functions
 * * only non-negativ parameters assumed
 *--------------------------------------------------------------*/
int
init_tinyExp( unsigned int N, double* obs, double* cond,
	double* a, unsigned char* a_flag, FILE* logfile )
{
	int j;
	double f = 0.1;
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
			a[j] = f * (float)random() / (float)RAND_MAX;
		}
	}
	return 0;
}
