/*****************************************************************
 *
 * File....:	outlier_detection.c
 * Function:	outlier detection
 * Author..:	Tilo Strutz
 * Date....:	03.07.2009, 05.01.2026
 *
 * changes:
 * 10.03.2011, 3.4.2011
 * 20.08.2012 implementation of RANSAC, M-score
 * 05.11.2025 tinyExpr
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
#include <assert.h>
#include <math.h>
#include <time.h>
#include "errmsg.h"
#include "matrixUtils.h"
#include "macros.h"
#include "ls.h"
#include "prototypes.h"
#include "erf.h"
#include "defines.h"
#include "tinyexpr.h"
 /* disables output from outlier detection */
/* #define COMPTEST */

/*---------------------------------------------------------------
 *	outlier_detection1()
 *
 *  outlier_detection based on standard deviation of obs.
 *	plus re-weighting
 *--------------------------------------------------------------*/
int
outlier_detection1( unsigned int N,/* number of observations	*/
				double sigma_y,		/* standard uncertainty	*/
				double *deviates,	/* absolute deviates		*/
				double *weights,	/* weights of observations	*/
				double nu,				/* Chauvenet's parameter		*/
	      FILE* out )
{
	char *rtn = "outlier_detection1";
	int err;
	int	count_outlier; /* counts the outliers	*/
	double erfval;
	double lambda_O;
	double kappa_O;

#ifndef COMPTEST
	fprintf( out, 
		"\n#\n# -- %s - start ------------------------------", rtn);
#endif

	if (nu > 0. && nu < 1.0)
	{
		/* set threshold according to Chauvenet's criterion,
		 * if nu is inside correct range
		 */
		err = erfinv( 1 - nu / N, &erfval);
		if (err) return 0;
		kappa_O = sqrt(2.) * erfval;
	}
	else
	{
		kappa_O = 4;	/* use fixed conservative threshold	*/
	}

	/* set breakdown point (threshold)	*/
	lambda_O = sigma_y * kappa_O;

#ifndef COMPTEST
	fprintf( out, "\n#  sigma_y = %f, kappa_O = %f, lambda_O = %f", 
					sigma_y, kappa_O, lambda_O);
#endif
	
	/* compare all deviates with threshold	*/
	count_outlier = 0;
	for (unsigned int i = 0; i < N; i++)
	{
		if (deviates[i] >= lambda_O)	/* take as outlier */
		{
			weights[i] = 0.0;
			count_outlier++;
		}
	}

#ifndef COMPTEST
	fprintf( out, 
		"\n# -- %s - end  ------------------------------", rtn);
#endif

	return count_outlier;
}

/*---------------------------------------------------------------
 *	outlier_detection2()
 *
 *  cluster-based outlier detection (ClubOD)
 *
 *--------------------------------------------------------------*/
int
outlier_detection2( unsigned int N,/* number of observations	*/ 
				double *deviates,	/* absolute deviates		*/
				double *weights,	/* weights of observations	*/
				FILE *out)
{
	char *rtn="outlier_detection2";
	int start_i;	/* range for outlier search	*/
	int stop;	/* flag	*/
	long *idx_dev=NULL;
	unsigned int cntZ;	/* counter for distances being equal to zero	*/
	unsigned int count_outlier=0; /* counts the outliers	*/
 	double relation, relation_max;
	double dist_loc, dist_ave;
	double *dist_glob = NULL;
	double *dist_all = NULL, dist_all_max;
	double *dist_sort = NULL;
	double *dev_sort = NULL;
	double lambda_O;
	double dist_thresh;
	double kappa1=3.; /* kappa to be used	*/

	/* number of observations for determination of kappa_1	*/ 
	const unsigned int NN[19]=
		{  8,  11,  16,    23,  32,    45,   64,   91, 128,  181, 256, 
			362, 512, 724, 1024, 1448, 2048, 2896, 4096
		};
	/* kappa_1 values according to NN[]	*/
	const double ka1[19]=
	{ 7.3, 7.7, 10.1, 11.8, 14.1, 16.7, 20.3, 25.2, 31.5, 39.6,  51.3, 
	  66.6, 86.4, 112, 150, 198.0, 261.0, 351.0, 
	};
	const double kappa2 = 2.0; 
	/* proportion of observations that should belong to the 
	 * bulk of good observations 
	 */
	const double kappa4 = 0.6; 
	
#ifndef COMPTEST
		fprintf( out, "\n# -- %s - start ----------------------", rtn);
#endif

		/* 
		 * memory allocation
		 */
		idx_dev = lvector( N);	/* index array	*/
		dev_sort = vector( N);	/* sorted absolute deviations	*/
		dist_all = vector( N);	/* distances between deviations	*/
		dist_sort = vector( N);	/* sorted distances	*/
		dist_glob = vector( N);	/* global distances	*/

		/* 
		 * determination of kappa_1 in dependent on number of 
		 * observations (threshold for suspicious distances)
		 */    
		if (N > NN[0])
		{
		  if (N >= NN[18]) 
					kappa1 = ka1[18];
			else
			{
				for (unsigned int i = 1; i < 19; i++)
				{
					if ( N <= NN[i])
					{
						double frac;

						frac = (double)(NN[i] - N) / (double)(NN[i] - NN[i-1]);
						/* linear interpolation	*/
						kappa1 = ka1[i] * (1-frac) + frac * ka1[i-1];
						break;
					}
				}
			}
		}
		else kappa1 = ka1[0];

		/* ascending sorting of deviates and indices */
		memcpy( dev_sort, deviates, sizeof(double) * N);
		heap_sort_d_( N, dev_sort, idx_dev);

		/*
		 * determine all distances
		 */
		for (unsigned int i = N-1; i > 0; i--)
		{
			dist_all[i] = dev_sort[i] - dev_sort[i-1];
		}

		/* dummy */
		dist_all[0] = 0.; 

		/* ascending sorting of distances in order to get 
		 * the median distance
		 */
		memcpy( dist_sort, dist_all, sizeof(double) * N);
		heap_sort_d( N, dist_sort);

#ifndef COMPTEST
		fprintf( out, "\n#  Number of observations: %d", N);
		fprintf( out, "\n#  kappa1: %f", kappa1);
		fprintf( out, "\n#  kappa2: %f", kappa2);
#endif

		/* keep major part as good observations */
		start_i = (int) ceil( kappa4 * N);

		/* determination of dist_glob, 
		 * separate value for each deviate
		 */
		{
			double arg, wj, sum_wj;
			/* for all observations, which have to be tested	*/
			for (unsigned int i = start_i; i < N; i++)
			{
				/* compute weighted average as d_glob,
				 * 2* sigma = N, variance = N*N/4
				 */
				arg = -0.5 * 4. / (N * N); 
				sum_wj = 0.;
				dist_glob[i] = 0;

				/* take all distances into account, exclude [0] */
				for (unsigned int j = 1; j < i; j++)
				{
					wj = 	exp( j*j * arg);
					/* dist_all[] sorted according increasing deviates
					 * largest fac2 for distance of highes deviate
					 */
					dist_glob[i] += dist_all[i-j] * wj;
					sum_wj += wj;
				}
				if (sum_wj > 0.)
					dist_glob[i] /= sum_wj;	/* normalise by sum	*/
				else
					dist_glob[i] = 0.0;
			}
		}
		
		if (dist_glob[start_i] == 0.0) /* check for smallest deviate */
		{
#ifndef COMPTEST
			fprintf( out, "\n#\n#  all distances are equal to zero!");
			fprintf( out, "\n#  nothing to weight, perfect fit!");
#endif
		}
		else
		{
			/* initialisation of thresholds */
			/* reset breakdown point */
			lambda_O = 0;
	
			/*
			 * search for suspicious distance 
			 */

#ifndef COMPTEST
			fprintf( out, 
				"\n#\n#   n  deviate    d[n]     d_glob");
			fprintf( out, "   d[n]/d_glob  d_loc  d[n]/d_loc ");
			fprintf( out, 
				"\n# -----------------------------------------------");
			fprintf( out, "--------------------");
#endif

			stop = 0;
			relation_max = 0.;
			dist_all_max = 0;

			/* check only upper portion of sorted deviates	*/
			for (unsigned int i = start_i; i < N && !stop; i++)
			{
				dist_thresh = dist_glob[i] * kappa1;

#ifndef COMPTEST
				fprintf( out, "\n# %3d %9.3e %9.3e", i, dev_sort[i], 
																								dist_all[i]);
				fprintf( out, " %9.2e %9.2e", 
					dist_glob[i], dist_all[i] / dist_glob[i]);
#endif

				/* if global condition is fulfilled, then check also local
				 * one
				 */
				if (dist_all[i] >= dist_thresh)
				{
					double arg, wj, sum_wj;

					/* determination of local distance value,
					 * 12* sigma = N, variance = N*N/144 
					 */
					arg = -0.5 * 144. / (N * N); 
					cntZ = 0;
					sum_wj = 0;
					dist_ave = 0;
					for (unsigned int j = 1; j < i; j++)
					{
						wj = 	exp( j*j * arg);
						/* largest fac2 for closest neighbour */
						dist_ave += dist_all[i-j] * wj;
						sum_wj += wj;
						/* count distances equal to zero	*/
						if (dist_all[i-j] == 0.0) cntZ++;
					}
					if (sum_wj)
						dist_loc = dist_ave / sum_wj;
					else
						dist_loc = 0.0;
					
					/* exception handling, if more than 50% of distances 
					 * are equal to zero
					 */
					if (cntZ <= i/2)
					{
						/* normal calculation	*/
						if (dist_loc > 0.0)
							relation = dist_all[i] / dist_loc;
						else
							relation = 0.0;
					}
					else /* exception handling	*/
					{
						/* current distance is > zero, border case */
						if (dist_all[i] > 0)
						{
							relation = kappa2;
						}
						else
						{
							/* not considered further */
							relation = 0;
						}
					}

#ifndef COMPTEST
					fprintf( out, " %9.3e", dist_loc);
					fprintf( out, " %8.2e", relation);
#endif

					/* check second condition */
					if (relation >= kappa2)
					{
						/* outliers found */

#ifndef COMPTEST
						fprintf( out, " *");
#endif
						/* look for the highest relation (change in distances)*/
						if (relation > relation_max)
						{
#ifndef COMPTEST
							fprintf( out, " ##");
#endif
							/* deviate of current position is taken as 
							 * breakdown point
							 */
							lambda_O = dev_sort[i];
							/* store current releation as maximal one	*/
							relation_max = relation;
							/* store current distance as maximum distance	*/
							dist_all_max = dist_all[i];
						}
						else if (relation == relation_max)
						{
							/* can happen, for instance, if relation was several 
							 * times equal to kappa2
							 */
							/* if equal, then take teh one corresponding to the
							 * larger absolute deviate
							 */
							if (dist_all[i] >= dist_all_max)
							{
								/* take the one with higher distance */
#ifndef COMPTEST
								fprintf( out, " ##");
#endif
							  lambda_O = dev_sort[i];
							  relation_max = relation;
								dist_all_max = dist_all[i];
							}
						}
						/*stop = 1; */
					}
				} /* if (dist_sort[i] > dist_thresh) */
			} /* for i  */

			/* check whether there was at least one relation 
			 * higher than kappa2
			 */
#ifndef COMPTEST
			if (relation_max > 0.0)
			{
				fprintf( out, "\n#\n#  outliers detected !");
				fprintf( out, "\t  lambda_O = %f", lambda_O);
			}
#endif

		  /* inspect unsorted data, set weights of outliers to zero*/
		  count_outlier = 0;
			if (lambda_O > 0)
			{
				for (unsigned int i = 0; i < N; i++)
				{
					if (deviates[i] >= lambda_O) /* outliers */
					{
							weights[i] = 0.;
							count_outlier++;
					}
				}
			}
		} /* 	if (enough observations) */
		
		/* output information for debugging */
#ifndef COMPTEST
		fprintf( out, "\n#\n#   n  deviate      ");
		fprintf( out, "           d[n]         weights[ idx[i] ]");
		for (unsigned int i = 0; i < N && i < MAX_LINES_W; i++)
		{
			fprintf( out, "\n#%4d  %14.9e %14.9e %16.8f", i,
					dev_sort[i], dist_all[i], 
					weights[ idx_dev[i] ]);
		}
#endif

		free_vector( &dist_glob);
		free_vector( &dev_sort);
		free_vector( &dist_sort);
		free_vector( &dist_all);
		free_lvector( &idx_dev);

#ifndef COMPTEST
  	fprintf( out, "\n# -- %s - end  -----------------------",rtn);
#endif

		return count_outlier;
}

/*---------------------------------------------------------------
 *	outlier_detection3()
 *
 *  outlier_detection based on Median Absolute Deviation (MAD)
 *	plus re-weighting
 *--------------------------------------------------------------*/
int
outlier_detection3( unsigned int N,/* number of observations	*/
				double *deviates,	/* absolute deviates		*/
				double *weights,	/* weights of observations	*/
				double nu,				/* Chauvenet's parameter		*/
				FILE *out)
{
	char *rtn = "outlier_detection3";
	int err;
	int	count_outlier; /* counts the outliers	*/
	double erfval;
	double lambda_O;
	double kappa_O;
	double *dev_sort = NULL, median_dev;

#ifndef COMPTEST
	fprintf( out, 
				"\n# -- %s - start ------------------------------", rtn);
#endif

	if (nu > 0. && nu < 1.0)
	{
		/* set threshold according to Chauvenet's criterion,
		 * if nu is inside correct range
		 */
		err = erfinv( 1 - nu / N, &erfval);
		if (err) return 0;
		kappa_O = sqrt(2.) * erfval;
	}
	else
	{
		kappa_O = 4;	/* use fixed conservative threshold	*/
	}
	
	/*
	 * determine median of absolute deviates
	 */
	dev_sort = vector( N);	/* sorted absolute deviations	*/
	/* ascending sorting of deviates and indices */
	memcpy( dev_sort, deviates, sizeof(double) * N);
	heap_sort_d( N, dev_sort);
	median_dev = dev_sort[N/2];

	/* set breakdown point (threshold)	*/
	/* z-score: lambda_O = sigma_y * kappa_O; */
	lambda_O = 1.4826 * median_dev * kappa_O;

#ifndef COMPTEST
	fprintf( out, 
	"\n#  1.4826 * median_dev = %f, kappa_O = %f, lambda_O = %f", 
					1.4826 * median_dev, kappa_O, lambda_O);
#endif

	/* compare all deviates with threshold	*/
	count_outlier = 0;
	for (unsigned int i = 0; i < N; i++)
	{
		if (deviates[i] >= lambda_O)	/* take as outlier */
		{
			weights[i] = 0.0;
			count_outlier++;
		}
	}

	free_vector( &dev_sort);

#ifndef COMPTEST
	fprintf( out, 
		"\n# -- %s - end  ------------------------------", rtn);
#endif

	return count_outlier;
}

/*---------------------------------------------------------------
 *	ransac()
 *
 *  outlier_detection based on RANSAC/MSAC method
 * This is actually not a method for outlier detection, but
 * a method selecting the most suitable modell-parameter vector 
 * out of many hypothesis of parameters.
 * 
 *--------------------------------------------------------------*/
int
ransac( double (*funct) (unsigned int,double*,double*),
	 double (*funct_deriv) (double(*)(unsigned int,double*,double*),
		 unsigned int, unsigned int, unsigned int,double*,double*),
	 double (*funct_deriv2) (double(*)(unsigned int,double*,double*),
		 unsigned int, unsigned int, unsigned int, unsigned int,double*,double*),
	 te_expr* teExpr,	 unsigned int N, unsigned int M, 
	double *obs, double *cond, double **jacob,
	 double *weights, double *a,
	 int algo_mode, LS_PARAMS *lsParams,
	 double chisqTarget, double **covar, FILE *out,
	 double *deviates_abs,
	unsigned int cond_dim,
	 unsigned int obs_dim, int digitsIntermediate )
{
	char *rtn = "ransac outlier detection";
	int err = 0;
	int iter, num_iterat;
	int rn, rm, idx;
	unsigned int im, Ns;
	int iterations=0;
	int num_outlier, cnt_inliers, cnt_inliers_best;
	double rin = 0.5;		/* percentage of inliers	*/
	double *s_a=NULL;		/* parameter of model function		*/
	double s_obj[M_MAX] = { 0 };	/* best parameter of model function	*/
	double s_a_best[M_MAX];	/* best parameter of model function	*/
	double a_backup[M_MAX];	/* required for tinyexpr, since original a must be overwritten	*/
	double costs, min_costs, threshold;
	double **s_jacob=NULL;/* array for subset of Jacobian matrix */
	double *s_obs = NULL; /* array for subset of observations */
	double *s_cond = NULL; /* array for subset of conditions */
	double *s_weights = NULL; /* array for subset of weights */
	/* array for sorted original deviates */
	double *dev_sort = NULL;	
	/* array for deviates based on subset approximation */
	double *deviates = NULL; 
	/* array of deviates based on best subset */
	double *deviates_best = NULL; 
	/* array for indices of selected data points */
	unsigned int *s_idx = NULL; 
	/* array for indices of selected data points of best subset */
	unsigned int *s_idx_best = NULL; 

#ifndef COMPTEST
	fprintf( out, 
				"\n# -- %s - start ------------------------------", rtn);
#endif

	/* set size of subset	*/
	if (lsParams->linear)
	{
		Ns = M; /* original setup of RANSAC:
					 * take minimum number of data required to fit the model
					 */ 
	}
	else
	{
		Ns = M+1; /* one more than necessary in case of nonlinear models 
							 * leads to somewhat higher numerical stability
							 */
	}

	/* 0.999 = probability that at least one subset is outlier-free */
	num_iterat = (int)(log(1-0.999) / log(1-pow(rin, Ns)));
	fprintf( out, "\n# required trials: %d", num_iterat);

	deviates = vector( N * obs_dim);	
	deviates_best = vector( N * obs_dim);	
	s_obs = vector( Ns * obs_dim);
	s_cond = vector( Ns * cond_dim);	
	s_weights = vector( Ns * obs_dim);
	s_idx = uivector( Ns);	
	s_idx_best = uivector( Ns);
	s_jacob = matrix( Ns * obs_dim, M);	/* Jacobian */
	

	/* todo: check whether local parameters s_a are required */
	if (lsParams->expression_flag)
	{
		/* tinyexpr operates only with original parameters !*/
		memcpy( a_backup, a, M_MAX * sizeof( double ) );
		s_a = a; /* point to original parameters */
	}
	else
	{
		s_a = s_obj; /* point to local memory */
	}

	/*
	 * determine threshold based on median of absolute deviates 
	 * based on current approximation
	 * threshold is a very critical parameter
	 *	if it is too low, then to many points are declared as outliers
	 *  if it is too high, then we might miss some outliers
	 */
	dev_sort = vector( N * obs_dim);	/* sorted absolute deviations	*/
	/* ascending sorting of deviates and indices */
	memcpy( dev_sort, deviates_abs, sizeof(double) * N * obs_dim);
	heap_sort_d( N * obs_dim, dev_sort);
	threshold = 3. * dev_sort[N * obs_dim/2];
	fprintf( out, "\n#    threshold: %f", threshold);

	/* initialize and zero out unnecessary parameters;
	 * especially required for POLYNOMIAL_REG
	 */
	for (unsigned int i = 0; i < M_MAX; i++)
	{
		s_a_best[i] = 0.;
	}

	/*
	 * create hypothesis based on randomly chosen subsets of
	 * the entire set of data points
	 */
  srand( (int)time(NULL));   /* Seed the random number generator. */
	idx = 0;
	min_costs = DBL_MAX; /* arbitrary large number */
	cnt_inliers_best = 0;
	for ( iter = 0; iter < num_iterat; iter++)
	{
		/* genetrate sub-set */
		/* The Knuth algorithm. 
		 * This is a very simple algorithm with a complexity of O(N) 
		 * The algorithm works as follows: 
		 * iterate through all numbers from 1 to N and 
		 * select the current number with probability rm / rn, 
		 * where rm is how many numbers we still need to find, 
		 * and rn is how many numbers we still need to iterate through. 
		 */
		im = 0;
		/* create array of indices (LUT)	*/
		fprintf( out, "\n# subset: ");
		for (unsigned int i = 0; i < N && im < Ns; i++)
		{
			rn = N - i;
			rm = Ns - im;
			if (rand() % rn < rm)
			{
				/* Take it */
				s_idx[im++] = i;
				fprintf( out, "%d ", i);
			}
		}
		assert( im == Ns);

		/* copy selected data */
		for (unsigned int i = 0; i < Ns; i++)
		{
			for (unsigned int j = 0; j < obs_dim; j++)
			{
				/* observations	*/
				s_obs[i * obs_dim+j] = obs[s_idx[i] * obs_dim+j];
				/* weights	*/
				s_weights[i * obs_dim+j] = weights[s_idx[i] * obs_dim+j];	
			}
			for (unsigned int j = 0; j < cond_dim; j++)
			{
				/* conditions	*/
				s_cond[i * cond_dim+j] = cond[s_idx[i]*cond_dim + j];	
			}
		}

		if (!lsParams->linear) /* nonlinear model is used */
		{
			if (cnt_inliers_best < 0.8 * N)
			{
				/* too less good data points; < 80%	*/
				for (unsigned int j = 0; j < M; j++)
				{
					s_a[j] = a[j]; /* reset to originally estimated	values	*/
				}
			}
			else /* use parameter set of best model */
			{
				for (unsigned int j = 0; j < M; j++)
				{
					s_a[j] = s_a_best[j]; 
				}
			}

			fprintf( out, "\n# initial Parameters\n# ");
			/* write initial parameters to output */
			for (unsigned int j = 0; j < M; j++)
			{
				fprintf( out, "a%d=%.9f, ", j+1, s_a[j]);
			}
		}
		else /* lsParams->linear */
		{
			//for (unsigned j = 0; j < M; j++)
			//{
			//	s_a[j] = 1.0; /* may not be zero */
			//}
		}

		/* do the least-squares approximation based on subset
		 * do not overwrite original Parameters in a
		 */
		if (lsParams->linear)
		{
			err =
				lsLinear( funct, funct_deriv,
					teExpr, 
					Ns * obs_dim, M,
					s_obs, s_cond, s_jacob, s_weights,
					s_a, algo_mode,
					lsParams, covar, out);
		}
		else
		{
			err =
				lsNonlinear( funct, funct_deriv, funct_deriv2,
					teExpr, Ns * obs_dim, M,
					s_obs, s_cond, s_jacob, s_weights,
					s_a, algo_mode,
					lsParams, chisqTarget, covar, &iterations, out,
					digitsIntermediate );
		}
		/* if lsLinear approximation failed, go to next subset	*/
		if (err) continue;

		fprintf( out, "\n#\n# subset Parameters:  ");
		for (unsigned int j = 0; j < M; j++)
		{
			fprintf( out, "a%d=%16.12G, ", j + 1, s_a[j]);
		}

		/* compute cost function based on all data points*/
		costs = 0; cnt_inliers = 0;
		for (unsigned int i = 0; i < N * obs_dim; i++)
		{
			/* decide between linear and nonlinear models */
			if (lsParams->linear) 
			{
				double model_val;

				/* get calculated data points dependent on current
					 parameters */
				model_val = 0.0;
				for (unsigned int j = 0; j < M; j++)
				{
					model_val += s_a[j] * jacob[i][j];
				}
				deviates[i] = fabs( model_val - obs[i]);
			}
			else /* nonlinear */
			{
				deviates[i] = fabs( funct( i, cond, s_a) - obs[i]);
			}
			/* check, whether data point (from entire set) is with 
			 * in the limits or not
			 */
			if (deviates[i] < threshold)
			{
				costs += deviates[i] * deviates[i];
				cnt_inliers++;	/* count inliers */
			}
			else
			{
				costs += threshold * threshold;
			}
		}
		fprintf( out, "  inliers: %d", cnt_inliers);
		/* evaluate hypothesis	*/
		if (min_costs > costs)
		{
			fprintf( out, "  ####### best subset so far (%f)", costs);
			fflush( out);
			min_costs = costs;
			/* save current subset as best subset	*/
			for (unsigned int i = 0; i < Ns; i++)
			{
				s_idx_best[i] = s_idx[i];
			}
			/* save deviates of model function based on this subset */
			for (unsigned int i = 0; i < N * obs_dim; i++)
			{
				deviates_best[i] = deviates[i];
			}
			for (unsigned int i = 0; i < M; i++)
			{
				s_a_best[i] = s_a[i];
			}
			cnt_inliers_best = cnt_inliers;

			/* adapt required number of subsets and increase value
			 * somewhat because we are not only in an outlier-free set 
			 * but in a set containing the maximum number of inliers
			 */
			rin = (double)cnt_inliers_best / (N * obs_dim);
			rin *= 0.9; /* this increases the number of trials slightly	*/
			num_iterat = (int)(log(1-0.9999) / log(1-pow( rin, Ns)));
			fprintf( out, "\n# required trials set to: %d", num_iterat);
			fprintf( out, "\t already seen: %d", iter+1);
		}
	}	/* for ( iter */

	num_outlier = N * obs_dim - cnt_inliers_best;
	if (num_outlier > 0.4*N)
	{
		fprintf( out, 
			"\n#\n# consensus set is too small: only %d out of %d!", 
				N * obs_dim-num_outlier, N);
		fprintf( out, 
			"\n#\n# reset to full data set with equal weights!");
		num_outlier = 0;
	}
	else
	{
		fprintf( out, "\n# number of trials: %d", iter);
		fprintf( out, "\n# number of outliers: %d\n#", num_outlier);
		/* for nonlinear models: copy best parameters as initial
		 * vector for final approximation (in case of outliers)
		 */
		if (!lsParams->linear && num_outlier) 
		{
			/* use best parameter set as initial parameter vector for final round */
			for (unsigned int j = 0; j < M; j++)
			{
				a[j] = s_a_best[j];
			}
		}
		/* mark outliers */
		for (unsigned int i = 0; i < N * obs_dim; i++)
		{
			if (deviates_best[i] >= threshold)
			{
				weights[i] = 0.;
			}
		}
	}

	if (lsParams->expression_flag)
	{
		/* restore original parameters */
		memcpy( a, a_backup, M_MAX * sizeof( double ) );
		/* perhaps not needed (?) */
	}

	free_vector( &dev_sort);
	free_vector( &deviates);
	free_vector( &deviates_best);
	free_vector( &s_obs);
	free_vector( &s_cond);
	free_vector( &s_weights);
	free_uivector( &s_idx);
	free_uivector( &s_idx_best);

#ifndef COMPTEST
	fprintf( out, 
		"\n# -- %s - end  ------------------------------", rtn);
#endif

	return num_outlier;
}