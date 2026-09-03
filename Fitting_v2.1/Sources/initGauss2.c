/*****************************************************************
 *
 * File........: init_gauss2.c
 * Function....: parameter initialisation for
 *                    data fitting with 2 Gaussians
 * Author......: Tilo Strutz
 * last changes: 02.07.2009
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
#include "matrixUtils.h"
#include "functions.h"

/*---------------------------------------------------------------
 *	init_gauss1()
 * f(x|a) = a1 * exp( (x-a2)^2 * a3)
 *--------------------------------------------------------------*/
int
init_gauss1( unsigned int N, double *obs, double *cond,
						 double *a, unsigned char *a_flag, 
						 int peak_flag, FILE *logfile)
{
	int err = 0;		 /* return value */
	int i_mean = 0, i_max, i_min;
	double max_val, min_val, condmin=0., condmax=0.;
	double mean, var, sum, sigma, tmp;
	
	/*
	 * get starting point 
	 * assuming that one Gaussian is good enough to fit the data
	 */
	
	/* get peak of curve */
	max_val = min_val = obs[0];
	i_max = i_min = 0;
	condmax = cond[0];
	condmin = cond[0];
	for (unsigned int i = 1; i < N; i++) /* let 1 sample border */
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
		fprintf( logfile, "\n\n Nothing to fit !!");
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
			var = var/sum - mean*mean;
		}
	}
	else
	{
		/* negative amplitude */
		for (unsigned int i = 0; i < N; i++)
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
	for (unsigned int i = 1; i < N; i++)
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
		if (peak_flag)
		{
			/* select highest peak, when there are 2 or more */
			if (!a_flag[0])	a[0] = max_val;
			if (!a_flag[1])	a[1] = condmax;
			/* reduce deviation accordingly */
			if (sigma > fabs(mean - condmax)) 
					sigma -= fabs(mean - condmax); 
		}
		else
		{
			/* increase by value dependent on obs at mean position */
			if (!a_flag[0])	
			{
				a[0] = max_val + (max_val - obs[i_mean]) * 0.5;
			}
			if (!a_flag[1])	a[1] = mean;
		}
	}
	else
	{
		if (peak_flag)
		{
			if (!a_flag[0])	a[0] = min_val;
			if (!a_flag[1])	a[1] = condmin;
			if (sigma > fabs(mean-condmin)) sigma -= fabs(mean-condmin);
		}
		else
		{
			if (!a_flag[0])	
			{
				a[0] = max_val + (max_val - obs[i_mean]) * 0.5;
			}
			if (!a_flag[1])	a[1] = mean;
		}
	}
	/* transcode deviation */
	if (!a_flag[2])	a[2] = -0.5 / (sigma*sigma);

endfunc:
	return err;
}

/*---------------------------------------------------------------
 *	init_gauss2()
 * f(x|a) = a1 * exp( a2 * (x-a3)^2) + 
 *          a4 * exp( a5 * (x-a6)^2)
 *--------------------------------------------------------------*/
int
init_gauss2( unsigned int N, double *obs, double *cond,
						double *a, unsigned char *a_flag, FILE *out)
{
	int err = 0;		 /* return value */
	int  M;
	double sigma1, sigma2;
	double *obs_cpy=NULL;
	
	M = 6;	/* fixed number of parameters */
	obs_cpy = vector( N);

	/*
	 * get starting point 
	 * assuming that one Gaussian is good enough to fit the data
	 */
	/* initialises a[0], a[1], a[2]; select highest peak */
	err = init_gauss1( N, obs, cond, a, a_flag, 1, out);
	if (err) goto endfunc;

	sigma2 = sqrt( -0.5 / a[2]); /* get deviation back */

	fprintf( out, "\n# Initial parameter: first Gaussian");
	fprintf( out, "\n# amplitude: %f, mean: %f, deviation: %f", 
						a[0], a[1], sigma2);
	/* 
	 * compute residuals between observation and model
	 */
	for (unsigned int i = 0; i < N; i++)
	{
			obs_cpy[i] = obs[i] - fgauss1( i, cond, a);
	}
	/* prevent overwriting of first three parameters */
	a[3] = a[0];
	a[4] = a[1];
	a[5] = a[2];

	/* initialises a[0], a[1], a[2]; fitting of residual */
	init_gauss1( N, obs_cpy, cond, a, a_flag, 1, out);
	sigma1 = sqrt( -0.5 / a[2]);
	fprintf( out, "\n# Initial parameter: second Gaussian");
	fprintf( out, "\n# amplitude: %f, mean: %f, deviation: %f", 
						a[0], a[1], sigma1);

	fprintf( out, "\n#\n# f(x) = %f * exp( ( x- %f)**2 * %f)",
					 a[0],a[1],a[2]);
	fprintf( out, "+ %f * exp( ( x- %f)**2 * %f)",
					 a[3],a[4],a[5]);
endfunc:
	free_vector( &obs_cpy);
	return err;
}
