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
#include "functions.h"
#include "macros.h"
#include "defines.h"
#include "prototypes.h"

#ifndef WIN32
#include <sys/time.h>
#else
#include <time.h>
#define random rand
#endif

 /*---------------------------------------------------------------
	*	init_fexpBaseX_log()
 *  f(x|a) = a2 * log(a1*x)
	*--------------------------------------------------------------*/
int init_fexpBaseX_log( int N, double *obs, double *cond,
	double *a, unsigned char *a_flag, FILE *logfile )
{
	long i;
	double maxCond;
	if (!a_flag[0])
	{
		maxCond = cond[0];
		for (i = 0; i < N; i++)
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
int init_fexpBaseX( int N, double *obs, double *cond,
	double *a, unsigned char *a_flag, FILE *logfile )
{
	long i;
	double maxCond, minObs;
	if (!a_flag[0])
	{
		maxCond = cond[0];
		for (i = 0; i < N; i++)
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
		for (i = 0; i < N; i++)
			if (minObs > obs[i]) minObs = obs[i];
		a[2] = minObs;
	}
	return 0;
}

/*---------------------------------------------------------------
 *	init_polynomial()
 * f(x|a) = sum_{j=1}^M aj * x^(j-1)
 *--------------------------------------------------------------*/
int init_polynomial( int N, double *obs, double *cond,
	double *a, unsigned char *a_flag, FILE *logfile )
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
 *	init_cosine_nonlin()
 * 5: f(x|a) = a1 + a2 * cos( x - a3)  (omega = 2*pi)
 *--------------------------------------------------------------*/
int
init_cosine_nonlin( int N, double *obs, double *cond,
					 double *a, unsigned char *a_flag, FILE *logfile)
{
	int i;
	double mean;
	
	/* mean value for a[0] */
	if (!a_flag[0])
	{
		mean = 0;
		for (i = 0; i < N; i++)
					mean += obs[i];
		a[0] = mean / N;
	}

	/* estimation of a2 = radius */
	{
		double max_obs, min_obs;
		max_obs = min_obs = obs[0];
		for (i = 1; i < N; i++)
		{
			if (max_obs < obs[i]) max_obs = obs[i];
			if (min_obs > obs[i]) min_obs = obs[i];
		}
		if (!a_flag[1]) /* if not set on command line */
		{
			a[1] = 0.5 * ( max_obs - min_obs);
		}
	}

	/* estimation of a3 = phase shift */
	if (!a_flag[2])
	{
		a[2] = 1.;	/* dummy */
	}
	
	return 0;
}

/*---------------------------------------------------------------
 *	init_cosine_trend()
 * 12: f(x|a) = a1 + a2 * x + a3 * cos( x - a4)
 *--------------------------------------------------------------*/
int
init_cosine_trend( int N, double *obs, double *cond,
					 double *a, unsigned char *a_flag, FILE *logfile)
{
	int i;
	double mean;
	
	/* mean value for a[0] */
	if (!a_flag[0])
	{
		mean = 0;
		for (i = 0; i < N; i++)
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
		for (i = 0; i < N; i++)
					mean += obs[i] * sqrt( 2.);
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
init_trigonometric1( int N, double *obs, double *cond,
					 double *a, unsigned char *a_flag, FILE *logfile)
{
	int i;
	double mean;
	
	/* mean value for a[0] */
	if (!a_flag[0])
	{
		mean = 0;
		for (i = 0; i < N; i++)
					mean += obs[i];
		a[0] = mean / N;
	}

	/* estimation of a3 = period */
	{
		double max_x, min_x;
		max_x = min_x = cond[0];
		for (i = 1; i < N; i++)
		{
			if (max_x < cond[i]) max_x = cond[i];
			if (min_x > cond[i]) min_x = cond[i];
		}
		if (!a_flag[2])
			a[2] = 2*3.141 / (2 * ( max_x - min_x));
	}

	/* estimation of a2 = amplitude */
	{
		double max_obs, min_obs;
		max_obs = min_obs = obs[0];
		for (i = 1; i < N; i++)
		{
			if (max_obs < obs[i]) max_obs = obs[i];
			if (min_obs > obs[i]) min_obs = obs[i];
		}
		if (!a_flag[1]) /* if not set on command line */
		{
			a[1] = 0.5 * ( max_obs - min_obs);
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
init_trigonometric2( int N, double *obs, double *cond,
					 double *a, unsigned char *a_flag, FILE *logfile)
{
	int i;
	double mean;
	
	/* mean value for a[0] */
	if (!a_flag[0])
	{
		mean = 0;
		for (i = 0; i < N; i++)
					mean += obs[i];
		a[0] = mean / N;
	}

	/* estimation of a3 = period */
	{
		double max_x, min_x;
		max_x = min_x = cond[0];
		for (i = 1; i < N; i++)
		{
			if (max_x < cond[i]) max_x = cond[i];
			if (min_x > cond[i]) min_x = cond[i];
		}
		if (!a_flag[2])
			a[2] = 2*3.141 / (2 * ( max_x - min_x));
	}

	/* estimation of a2, a5  = amplitude */
	{
		double max_obs, min_obs;
		max_obs = min_obs = obs[0];
		for (i = 1; i < N; i++)
		{
			if (max_obs < obs[i]) max_obs = obs[i];
			if (min_obs > obs[i]) min_obs = obs[i];
		}
		if (!a_flag[1]) /* if not set on command line */
		{
			a[1] = 0.5 * ( max_obs - min_obs);
		}
		if (!a_flag[4]) /* if not set on command line */
		{
			a[4] = rand() * 0.5 * ( max_obs - min_obs);
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
 *  f(x|a) = log( a1 * x) 
 *--------------------------------------------------------------*/
int
init_logarithmic( int N, double *obs, double *cond,
					 double *a, unsigned char *a_flag, FILE *logfile)
{
	if (!a_flag[0])
	{
		unsigned int idxMin = 0, idxMax = 0;
		double minCond, maxCond;
		minCond = maxCond = cond[0];
		for (unsigned int i = 1; i < N; i--)
		{
			if (minCond > cond[i]) { minCond = cond[i]; idxMin = i;	}
			if (maxCond < cond[i]) { maxCond = cond[i]; idxMax = i; }
		}
		if (minCond > 1.)
		{
			if (obs[idxMin] > 0) a[0] = 1.0;
			else a[0] = 0.0001;
		}
		else
		{
			if (minCond > 0.0)
				a[0] = exp( obs[idxMin] ) / minCond;
			else
				a[0] = 1.;
		}
	}
	return 0;
}
/*---------------------------------------------------------------
 *	init_exponential()
 *  f(x|a) = a1 + a2 * exp( a3 * x) 
 *--------------------------------------------------------------*/
int
init_exponential( int N, double *obs, double *cond,
					 double *a, unsigned char *a_flag, FILE *logfile)
{
	int err = 0;		 /* return value */
	int i, itmp;
	double mean;
	
	/* number of conditions to be inspected */
	itmp = MAX( 0, MIN( 5, N));
	
	/* estimation of a1 = tail of graph */
	if (!a_flag[0])
	{
		mean = 0;
		for (i = N - 1; i > N - itmp; i--)
			mean += obs[i];
		a[0] = mean / itmp;
	}

	/* estimation of a2 = head of function */
	if (!a_flag[1])
	{
		mean = 0;
		for (i = 0; i < itmp; i++)
			mean += obs[i];

		/* assumes conditions starting close to zero
		 * y(x=0) = a1 + a2 * exp( a3 * 0) = a1 + a2
		 */
		a[1] = mean / itmp;
	}

	/* estimation of a3 = gradient at head of function */
	if (!a_flag[2])
	{
		mean = 0;
		for (i = 1; i < itmp; i++)
			mean += (obs[i] - obs[i - 1]) / (cond[i] - cond[i - 1]);
		a[2] = mean / (itmp * 0.5);
	}

	/* if not decaying */
	if (a[0] > a[1])
	{
		fprintf( stderr, "\n a1 > a2 !");
		fprintf( stderr, "\n flip signs of a2 and a3 !\n");
		a[1] = -a[1];
		a[2] = -a[2];
	}
		 
	return err;
}

/*---------------------------------------------------------------
 *	init_expon2()
 * f(x|a) = a1 * exp( a2 * x)
 *--------------------------------------------------------------*/
int
init_expon2( int N, double *obs, double *cond,
	double *a, unsigned char *a_flag, FILE *out)
{
	int i, itmp;
	double mean;

	/* number of conditions to be inspected */
	itmp = MAX( 0, MIN( 5, N));

	/* estimation of a1 = head of function */
	if (!a_flag[0])
	{
		mean = 0;
		for (i = 0; i < itmp; i++)
			mean += obs[i];
		
		/* assumes conditions starting close to zero
		 * y(x=0) = a1 * exp( a2 * 0) = a1
		 */
		a[0] = mean / itmp;
	}

	/* estimation of a2 = gradient at head of function
	 * a2 = f'(0)/a1 
	 */
	if (!a_flag[1])
	{
		mean = 0;
		for (i = 1; i <= itmp; i++)
			mean += (obs[i] - obs[i - 1]) / (cond[i] - cond[i - 1]);
		a[1] = mean / (itmp * a[0]);
	}
	return 0;
}

/*---------------------------------------------------------------
 *	init_gauss()
 * f(x|a) = a1 * exp( a2 * (x-a3)^2) + 
 *--------------------------------------------------------------*/
int
init_gauss( int N, double *obs, double *cond,
						double *a, unsigned char *a_flag, FILE *out)
{
	int err = 0;		 /* return value */
	int i;
	int i_mean = 0, i_max, i_min;
	double max_val, min_val, condmin=0., condmax=0.;
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
	for (i = 2; i < N-1; i++) /* let 1 sample border */
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
		fprintf( out, "\n\n Nothing to fit !!");
		a[0] = 0.;
		a[2] = -50000000.0;
		a[1] = 0.;
		err = 8;
		goto endfunc;
	}
	
	mean = sum = var = 0.;
	/* take only that part which has the highest peak */
	if (fabs(max_val) > fabs(min_val))
	{
		/* positive amplitude */
		for (i = 0; i < N; i++)
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
			var = var/sum - mean*mean;
		}
	}
	else
	{
		/* negative amplitude */
		for (i = 0; i < N; i++)
		{
			if (obs[i] < 0.)
			{
				tmp = - cond[i] * obs[i];
				mean += tmp;
				var += cond[i] * tmp;
				sum -= obs[i];
			}
		}	
		if (sum > 0.)
		{
			mean /= sum;
			var = var/sum - mean*mean;
		}
	}

	/* if only one data point, then sigma is zero */
	if (var > 0.)	sigma = sqrt( var);  /* deviation of Gaussian */
	else
	 sigma = 0.0000001;

	
	/* get index of mean position */
	for (i = 1; i < N; i++)
	{
		if (cond[i-1] <= mean && mean <= cond[i])
		{
			i_mean = i;	/* mean position */
			break;
		}
	}
	
	/* make values more robust by averaging */
	max_val = (max_val + obs[i_max-1] + obs[i_max+1]) /3;
	min_val = (min_val + obs[i_min-1] + obs[i_min+1]) /3;
	if (obs[i_mean] > 0.)
	{
		{
			/* select highest peak, when there are 2 or more */
			if (!a_flag[0])	a[0] = max_val;
			if (!a_flag[1])	a[1] = condmax;
			/* reduce deviation accordingly */
			if (sigma > fabs(mean - condmax)) 
					sigma -= fabs(mean - condmax); 
		}
	}
	else
	{
		{
			if (!a_flag[0])	a[0] = min_val;
			if (!a_flag[1])	a[1] = condmin;
			if (sigma > fabs(mean-condmin)) sigma -= fabs(mean-condmin);
		}
	}
	/* transcode deviation */
	if (!a_flag[2])	a[2] = -0.5 / (sigma*sigma);

endfunc:
	return err;
}

/*---------------------------------------------------------------
 *	init_gen_laplace()
 * f(x|a) = a1 * exp( -|x|^a2 * a3)
 *--------------------------------------------------------------*/
int
init_gen_laplace( int N, double *obs, double *cond,
	double *a, unsigned char *a_flag, FILE *out)
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
init_circlelin( int N, double *obs, double *cond,
	double *a, unsigned char *a_flag, FILE *logfile)
{
	int err = 0;		 /* return value */
	double b1, b2, b3;

	/* get estimates of centre coordinates and radius */	
	err = init_circle( N, obs, cond, a, a_flag, logfile);

	/* convert into vector b	*/
	b1 = 2 * a[0];
	b2 = 2 * a[1];
	b3 = a[0]*a[0] + a[1]*a[1] - a[2]*a[2];

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
init_circle( int N, double *obs, double *cond,
	double *a, unsigned char *a_flag, FILE *logfile)
{
	int err = 0;		 /* return value */
	int i;
	double sum_x, sum_y, rad2, diff1, diff2;

	fprintf( logfile, "\n#\n# init_circle()");
	
	/* 
	 * determine circle centre
	 */

	/* compute centroids of conditions */
	sum_x = sum_y = 0.;
	/* two conditions */
	for (i = 0; i < 2*N; i+=2)
	{
		sum_x += cond[i];
		sum_y += cond[i+1];
	}
	sum_x /= (double)N;
	sum_y /= (double)N;
	
	rad2 = 0;
	for (i = 0; i < 2*N; i+=2)
	{
		diff1 = cond[i] - sum_x;
		diff2 = cond[i+1] - sum_y;
		rad2 += sqrt( diff1*diff1 + diff2*diff2);
	}
	rad2 = rad2 / (double)N;
	fprintf( logfile, "\n#\n# mean of condition coordinates");
	fprintf( logfile, "\n#   mean(x)= %f", sum_x);
	fprintf( logfile, "\n#   mean(y)= %f", sum_y);
	fprintf( logfile, "\n#   radius = %f", rad2);

	if (!a_flag[0])	a[0] = sum_x;
	if (!a_flag[1])	a[1] = sum_y;
	if (!a_flag[2])	a[2] = rad2;
	
	fprintf( logfile, 
					 "\n#  f(x|a) =0= (x1-%f)**2 + (x2-%f)**2 - %f**2", 
					 a[0], a[1], a[2]);
		 
	return err;
}

/*---------------------------------------------------------------
 *	init_rotation()
 * 21... f1(x|a) = a1 + cos(a3) * x1 - sin(a3) * x2
 *       f2(x|a) = a2 + sin(a3) * x1 + cos(a3) * x2
 *-------------------------------------------------------------*/
int
init_rotation( int N, double *obs, double *cond,
	double *a, unsigned char *a_flag, FILE *logfile)
{
	int err = 0;		 /* return value */
	int i;
	double sum_x, sum_y, sum_u, sum_v;

	fprintf( logfile, "\n#\n# init_rotation()");
	
	/* 
	 * determine rough translation
	 */

	/* compute centroids of conditions and observations */
	sum_x = sum_y = 0;
	sum_u = sum_v = 0;
	/* assume double observations and conditions */
	for (i = 0; i < N * 2; i+=2)
	{
		sum_x += obs[i];
		sum_y += obs[i+1];
		sum_u += cond[i];
		sum_v += cond[i+1];
	}
	sum_x /= (double)N;
	sum_y /= (double)N;
	sum_u /= (double)N;
	sum_v /= (double)N;
	fprintf( logfile, "\n#\n# mean of condition coordinates");
	fprintf( logfile, "\n#   mean(u)= %f", sum_u);
	fprintf( logfile, "\n#   mean(v)= %f", sum_v);
	fprintf( logfile, "\n# mean of observed coordinates");
	fprintf( logfile, "\n#   mean(x)= %f", sum_x);
	fprintf( logfile, "\n#   mean(y)= %f", sum_y);

	if (!a_flag[0])	a[0] = sum_x - sum_u;
	if (!a_flag[1])	a[1] = sum_y - sum_v;
	if (!a_flag[2])	a[2] = 0;		/* assume no rotation */
	
	
	fprintf( logfile, 
					 "\n#  f1(u,v) = %f + cos(%f) * u - sin(%f) * v", 
					 a[0], a[2], a[2]);
	fprintf( logfile, 
					 "\n#  f2(u,v) = %f + sin(%f) * u + cos(%f) * v", 
					 a[1], a[2], a[2]);
	 
	return err;
}

/*---------------------------------------------------------------
 *	init_NN3x3x1()
 * 3x3x1 
 *--------------------------------------------------------------*/
int
init_NN3x3x1( int N, double *obs, double *cond,
							double *a, unsigned char *a_flag, FILE *logfile)
{
	int i, j;
	double minval, maxval;
#ifndef WIN32
	struct timeval tv;
	struct timezone tz;
	
	gettimeofday( &tv, &tz);
	srandom( tv.tv_sec);
#else
   /* Seed the random-number generator with current time so that
    * the numbers will be different every time we run.
    */
   srand( (unsigned)time( NULL ) );
#endif
  
	/* give parameters random values */
	for ( j = 0; j < M_MAX; j++)
	{
		if (!a_flag[j])	
			a[j] = 2. * (float)random()/ (float)RAND_MAX - 1.;
	}
	
	/* make random numbers in a range that |cond x param| < 5 */
	minval = maxval = cond[0];
	for ( i = 1; i < N; i++)
	{
		if (minval > cond[3*i]) minval = cond[3*i];
		if (maxval < cond[3*i]) maxval = cond[3*i];
	}
	/* weights from 1st input */
	if (!a_flag[1])	a[1] = a[1] / (maxval-minval);
	if (!a_flag[5])	a[5] = a[5] /	(maxval-minval);
	if (!a_flag[9]) a[9] = a[9] / (maxval-minval);
	
	minval = maxval = cond[1];
	for ( i = 1; i < N; i++)
	{
		if (minval > cond[3*i+1]) minval = cond[3*i+1];
		if (maxval < cond[3*i+1]) maxval = cond[3*i+1];
	}
	/* weights from 2nd input */
	if (!a_flag[2])	a[2] = a[2] / (maxval-minval);
	if (!a_flag[6])	a[6] = a[6] / (maxval-minval);
	if (!a_flag[10])a[10]= a[10] / (maxval-minval);
	
	minval = maxval = cond[2];
	for ( i = 1; i < N; i++)
	{
		if (minval > cond[3*i+2]) minval = cond[3*i+2];
		if (maxval < cond[3*i+2]) maxval = cond[3*i+2];
	}
	/* weights from 3rd input */
	if (!a_flag[3])	a[3] = a[3] / (maxval-minval);
	if (!a_flag[7])	a[7] = a[7] / (maxval-minval);
	if (!a_flag[11])a[11]= a[11]/ (maxval-minval);
	
	return 0;
}

/*---------------------------------------------------------------
 *	init_NN1x3x1()
 * 1x3x1 
 *--------------------------------------------------------------*/
int
init_NN1x3x1( int N, double *obs, double *cond,
							double *a, unsigned char *a_flag, FILE *logfile)
{
	int i, j;
	double minval, maxval;
#ifndef WIN32
	struct timeval tv;
	struct timezone tz;
	
	gettimeofday( &tv, &tz);
	srandom( tv.tv_sec);
#else
   /* Seed the random-number generator with current time so that
    * the numbers will be different every time we run.
    */
   srand( (unsigned)time( NULL ) );
#endif
  
	/* give parameters random values */
	for ( j = 0; j < M_MAX; j++)
	{
		if (!a_flag[j])
		{
			a[j] = 2. * (float)random() / (float)RAND_MAX - 1.;
		}
	}
	
	/* make random numbers in a range that |cond x param| < 5 */
	minval = maxval = cond[0];
	for ( i = 1; i < N; i++)
	{
		if (minval > cond[i]) minval = cond[i];
		if (maxval < cond[i]) maxval = cond[i];
	}
	/* weights from 1st input */
	if (!a_flag[1])	a[1] = a[1] / (maxval-minval);
	if (!a_flag[3])	a[3] = a[3] / (maxval-minval);
	if (!a_flag[5])	a[5] = a[5] / (maxval-minval);
	
	return 0;
}

/*---------------------------------------------------------------
 *	init_NN()
 * 3x... 
 *--------------------------------------------------------------*/
int
init_NN( int N, double *obs, double *cond,
				 double *a, unsigned char *a_flag, FILE *logfile)
{
	int j;
#ifndef WIN32
	struct timeval tv;
	struct timezone tz;
	
	gettimeofday( &tv, &tz);
	srandom( tv.tv_sec);
#else
   /* Seed the random-number generator with current time so that
    * the numbers will be different every time we run.
    */
   srand( (unsigned)time( NULL ) );
#endif
  
	for ( j = 0; j < M_MAX; j++)
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
init_tinyExp( int N, double* obs, double* cond,
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
