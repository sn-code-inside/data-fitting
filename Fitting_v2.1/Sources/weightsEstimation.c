/*****************************************************************
 *
 * File........:	est_weights.c
 * Function....:	estimation of weights
 * Author......:	Tilo Strutz
 * last changes:	03.07.2009
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
#include "errmsg.h"
#include "matrixUtils.h"
#include "macros.h"
#include "prototypes.h"
#include "defines.h"

/*---------------------------------------------------------------
 *	est_weights1()
 *
 *  weights estimation based on deviates, 50% equal weights
 *	no dependence on standard uncertainty of observations
 *--------------------------------------------------------------*/
void est_weights1( unsigned int N, 
	double *deviates, /* absolute deviates ! */
	double *weights, FILE *out)
{
	char *rtn="est_weights1";
	long *idx_dev=NULL;
	double *dev_sort = NULL;
	double lambda_L, lambda_L2, lambda_L2_inv;
	double max_deviate, bound;
	const double kappa_L = 0.05; 
	
		fprintf( out, "\n# -- %s - start ----------------------", rtn);
		/* 
		 * memory allocation
		 */
		idx_dev = lvector( N);
		dev_sort = vector( N);
		
		/* ascending sorting of deviates and indices */
		memcpy( dev_sort, deviates, sizeof(double) * N);
		for (unsigned int i = 0; i < N; i++)
		{
			if (weights[i] == 0.) 
			{
				/* weights can already be set to zero by purpose
				 * in the linearisation process
				 * corresponding deviates must be ignored somehow
				 */
				dev_sort[i] = 0.;
			}
		}
		/* sorting of absolute deviates and of an index array */
		heap_sort_d_( N, dev_sort, idx_dev);

		/* get max. of all absolute deviates */
		max_deviate = dev_sort[N-1];

		fprintf( out, "\n#  Number of observations: %d", N);
		fprintf( out, "\n#  max_deviate: %f", max_deviate);
		
		/* 
		 * check values
		 */
		if (max_deviate == 0.0)
		{
			fprintf( out, "\n#\n#  all deviates are equal to zero!");
			fprintf( out, "\n#  nothing to weight, perfect fit!");
		}
		else
		{
			/*
			 * initialisation of threshold
			 */
			/* keep half of the values with equal weights */
			lambda_L = dev_sort[N/2];
			bound = max_deviate * kappa_L;
			/* ensure minimum value for this threshold */
			if (lambda_L < bound)
			{
				/* value for worst case, stabilising LS	*/
				lambda_L = bound; 
			}
			lambda_L2 = lambda_L * lambda_L;
			fprintf( out, "\n#  lambda_L = %f", lambda_L);
			
			/*
			 * do the weighting 
			 */
		  /* inspect unsorted data */
		  lambda_L2_inv = 1. / lambda_L2;
			for (unsigned int i = 0; i < N; i++)
		  {
				if (weights[i] > 0.)
				{
					/* weights can already be set to zero by purpose
					 * in the linearisation process
					 */
					if (deviates[i] < lambda_L)
					{
						/* fixed weight */
						weights[i] = lambda_L2_inv;
					}
					else
					{
						/* adapted weight */
						weights[i] = 1. / (deviates[i] * deviates[i]); 
					}
				}
			}

		} /* 	if max_deviate > 0 */
		
		/* output information for debugging */
		fprintf( out, "\n#\n#  dev_sort[i]          weights[ idx[i] ]");
		for (unsigned int i = 0; i < N && i < MAX_LINES_W; i++)
		{
			fprintf( out, "\n#%4d  %14.9e %16.8f", i,
					dev_sort[i], weights[ idx_dev[i] ]);
		}

		free_vector( &dev_sort);
		free_lvector( &idx_dev);
  	fprintf( out, "\n# -- %s - end  -------------------------",rtn);
}

/*---------------------------------------------------------------
 *	est_weights2()
 *
 *  weights estimation based on binning
 *--------------------------------------------------------------*/
void est_weights2( unsigned int N, double *cond, double *obs, 
							double *weights, unsigned int obs_per_bin, FILE *out)
{
	char *rtn="est_weights2";
	int n, m0, b, cnt;
	unsigned int i, m;
	int num_bins;
	long *idx=NULL;
	double w, var, diff, a[2], last_cond;
	double *cond_sort=NULL;
	double *bin_cond = NULL;
	double *bin_obs = NULL;
	
		fprintf( out, 
			"\n#  -- %s - start ------------------------------", rtn);
		
		/* vector of conditions in a single bin	*/
		bin_cond = vector( obs_per_bin);
		/* vector of observations in a single bin	*/
		bin_obs = vector( obs_per_bin);
		
		/* determine the number of bins	*/
		num_bins = N / obs_per_bin;
		/* one bin for the rest	*/
		if (N % obs_per_bin) num_bins++;
	
		cond_sort = vector(N); /* array for sorted conditions */
		idx = lvector( N);	   /* vector for sorted indices   */
		
		/* copy all conditions */
		memcpy( cond_sort, cond, sizeof(double)*N);

		/* ascending sorting of conditions and indices */
		heap_sort_d_( N, cond_sort, idx);

		/* initialise last condition value; just for output */
		last_cond = cond_sort[0];

		fprintf( out, "\n#  Number of bins: %d", num_bins);
		fprintf( out, "\n#  bin   #observ.   a1          ");
		fprintf( out, "a2           last_cond   variance");

		n = 0; /* n... already used observations */

		/* for all bins */
		for ( b = 0; b < num_bins; b++)
		{
			if ( N-n < obs_per_bin/ 2)
			{
				/* applies for the last 2 bins:
				 * if too less remaining observation, 
				 * then let bins overlap
				 */
				m0 = n - obs_per_bin/2;
				obs_per_bin = N-m0;
			}
			else m0 = n;

			/*
			 * piece-wise linear approximation per bin
			 */

			/* copy bin related values */
			cnt = 0;
			for (i = 0, m = m0; i < obs_per_bin && m < N; i++, m++)
			{
				bin_cond[i] = cond_sort[m];
				/* bin_cond[i] = cond[ idx[m] ];  same */
				bin_obs[i] = obs[ idx[m] ];
				cnt++;
			}
			/* determine parameters for piece-wise linear fit */
			ls_straightline( cnt, bin_cond, bin_obs, a);

			/* determine uncertainty */
			m = m0;
			var = 0.;
			for (i = 0; i < obs_per_bin && m < N; i++, m++)
			{
				diff = obs[ idx[m] ] - (a[0] + a[1] * bin_cond[i]);			
				var += diff * diff;
			}
			var = var / (double)cnt;
			if (var > 0.)
				w = 1. / var;
			else w = 0.;
			fprintf( out, 
					"\n#  %4d %3d    %12.4e %12.4e %12.4e %12.3e", 
					b, i, a[0], a[1], bin_cond[i-1], var);
			fprintf( stdout, 
				"\n  g%d(x) = (x > %.4f && x < %.4f) ? %.4f +x*  %.4f : 0", 
				b, last_cond, bin_cond[i-1],  a[0], a[1]);

			last_cond = bin_cond[i-1];
			
			/* klone determined bin weight for all corresponding
			 * observations
			 */ 
			for ( i = 0, m = m0; i < obs_per_bin && m < N; i++, m++)
			{
					weights[ idx[m] ] = w;
			}
			n = m; /* update number of already used observations */
		}

		free_vector( &bin_cond);
		free_vector( &bin_obs);
		free_vector( &cond_sort);
		free_lvector( &idx);
  	fprintf( out, 
			"\n#  -- %s - end  ------------------------------", rtn);
}
