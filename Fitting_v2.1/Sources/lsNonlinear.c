/*****************************************************************
 *
 * File........:	lsNonlinear.c
 * Function....:	least squares with alternative matrix inversion
 * Author......:	Tilo Strutz
 * last changes:	05.08.2026
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
#include <float.h>
#include "errmsg.h"
#include "matrixUtils.h"
#include "defines.h"
#include "macros.h"
#include "ls.h"
#include "prototypes.h"
#include "functions.h"
#include "tinyexpr.h"
#ifndef WIN32
#include <sys/time.h>
#else
#include <time.h>
#define random rand
#endif

#define NUMTARGETTRIALS 100
#define GRAD_THRESH 1e-3

long ITERAT_MAX;

extern double x; /* declared in fitting.c */

/*---------------------------------------------------------------
*	lsNonlinear()
*
*--------------------------------------------------------------*/
int
lsNonlinear( double (*funct) (unsigned int, double*, double*),
	double (*funct_deriv) (double(*)(unsigned int, double*, double*),
		unsigned int, unsigned int, unsigned int, double*, double*),
	double (*funct_deriv2) (double(*)(unsigned int, double*, double*),
		unsigned int, unsigned int, unsigned int, unsigned int, double*, double*),
	te_expr* teExpr,
	unsigned int N, unsigned int M, double* observations, double* conditions, double** jacob,
	double* weights,
	double* a,
	int algoMode, LS_PARAMS* lsParams,
	double chisqTarget, double** normal_i, int* iterations, FILE* out,
	int digitsIntermediate )
{
	char* rtn = "lsNonlinear";
	int err = 0;
	int out_flag;/* if number of iteration is to high avoid too many outputs */
	int chisqTargetCounter = 0;
	int converged_flag = 0;
	int targedReached_flag = 0;
	int iterLimitReached_flag = 0;
	unsigned int iter_max = ITERAT_MAX;	/* maximum number of iterations */
	double** cofac = NULL;	/* cofactor matrix for matrix
													inversion */
	double** normal = NULL;			/* N = J^(T) * W * J	*/
	double* gradientVec = NULL;			/* J^(T) * W * r			*/
	double* datac = NULL;	/* calculated values based on parameters
												*/
	double* da = NULL;		/* parameter update delta a */
	double residual, deriv_2nd;
	double* last_a = NULL; /* remember last result */
	double* min_a = NULL; /* remember best result */
	double* best_a = NULL; /* remember best of ll trials */
	double* initial_a = NULL; /* remember initial values for randomization */
	double closestToChisqTarget = DBL_MAX;
	double max_diag;	/* maximum value of Njj */
	double muReduceFac = 0.6666666666666666;
	double muIncreasFac = 4.0;/* should not be a multiple of decreasing factor */
	double sumOfSquaredResiduals = 0., sumOfSquaredResiduals0 = 1.;
	double gradientThreshold[M_MAX];
	int gradientSign[M_MAX];
	int initialWatchDogCounter = 10;

	LS_STATE lsStateObj = { 0 }, * lsState;

	lsState = &lsStateObj;

	fprintf( out,
		"\n#  -- %s - start  ------------------------------", rtn );
	lsState->minChisq = lsState->chisq = DBL_MAX;
	lsState->lastChisq = lsState->bestChisq = DBL_MAX;
	lsState->chisqWatchDog = initialWatchDogCounter; /* allow  iterations of equal costs */
	/*
	* allocate memory
	*/

	/* normal matrix N = J^(T) * W * J, its inverse */
	normal = matrix( M, M );
	cofac = matrix( M, M );	/* cofactor matrix */
	gradientVec = vector( M );		/* container for J^(T) * W * G	*/
	datac = vector( N );		/* calculated data using f(x|a) */
	da = vector( M );				/* model parameter update */
	min_a = vector( M );			/* remember best parameter set of current iteration */
	last_a = vector( M );			/* remember last parameter set */
	best_a = vector( M );			/* remember best parameter set of all trials */
	initial_a = vector( M );			/* remember initial parameter set  for randomisation*/

	lsState->daFac = 1.0;
	lsState->iterCnt = 1;

	/* initialise minimum (best) values */
	for (unsigned int j = 0; j < M; j++)
	{
		min_a[j] = initial_a[j] = last_a[j] = best_a[j] = a[j]; /* initial estimtes */
		gradientThreshold[j] = GRAD_THRESH; /* maximum threshold */
		gradientSign[j] = 0; /* sign not known yet */
	}

	if (lsParams->LM)
	{  /* Levenberg-Marquard */
		fprintf( out, "\n#  initial mu  factor: %lf", lsParams->muFac );
		fprintf( out, "\n#  mu decrease factor: %lf", muReduceFac );
		fprintf( out, "\n#  mu increase factor: %lf", muIncreasFac );
	}

	/* compute initial chi squared */
	lsState->chisq = calcChisqNonLin( funct, N, observations, conditions, a, weights, lsParams->expression_flag, teExpr );
	lsState->bestChisq = lsState->lastChisq = lsState->minChisq = lsState->chisq;
	fprintf( out, "\n#\n#  initial chisq=%.*e", OP_DBL_Digs, lsState->minChisq );

	lsParams->uphillAcceptanceThreshold = 0.0; /* never uphill */
	if (!lsParams->LM)
	{
		if (lsParams->GNspecial_flag)
		{
			lsParams->uphillAcceptanceThreshold = lsState->chisq; /* allow 100% increase at the beginning */
			lsState->muHessian = 0.; /* without true Hessian in the beginning */
			fprintf( out, "\n#  GN: uphill activated, damping of steps forcing convergence" );
			if (lsParams->trueH)
				fprintf( out, ", adapted Hessian" );
			/* initial impact of residuals*/
			sumOfSquaredResiduals0 = TOL; /* prevent division by zero */
			for (unsigned int i = 0; i < N; i++)
			{
				/* damp influence of Q by myHessian */
				if (lsParams->expression_flag)
				{
					x = conditions[i]; /* uses global variable x */
					residual = observations[i] - te_eval( teExpr );
				}
				else
				{
					residual = observations[i] - funct( i, conditions, a );
				}
				sumOfSquaredResiduals0 += residual * residual; /* for adaptation of muHessian */
			}
		}
		else
		{
			lsState->muHessian = 1.; /* full Hessian from the beginning */
			fprintf( out, "\n#  GN: uphill deactivated, no damping of steps" );
			if (lsParams->trueH)
				fprintf( out, ", full Hessian" );
		}
	}

	/* iteration of nonlinear least squares */
	do
	{
		/* feedback on console, the more iterations the less output */
		fprintf( stderr, "\r\t %8d", lsState->iterCnt );
		if (lsState->iterCnt < 3000 || (lsState->iterCnt < 10000 && (lsState->iterCnt % 1000) == 0)
			|| (lsState->iterCnt < 100000 && (lsState->iterCnt % 10000) == 0)
			|| (lsState->iterCnt % 100000) == 0)
		{
			out_flag = 1;
		}
		else out_flag = 0;

		if (out_flag)
		{
			fprintf( out, "\n#\n#  Iteration of least squares: %d",
				lsState->iterCnt );
			fflush( out );
		}
		/* get Jacobian matrix  */
		err = getJacobian( lsParams->expression_flag, funct, funct_deriv,
			teExpr, N, M, conditions, jacob, a, out );

		/*
		* calculate normal matrix N
		* N = J^(T) * W * J
		*/
		err = computeNormalMatrix( N, M, jacob, normal, weights );

		/*
		 * compute the negative gradient vector
		 * -g = J^(T) * W * r
		 */
		err = computeJWr( funct, N, M, observations, conditions, jacob,
			weights, a, gradientVec, lsParams->expression_flag, teExpr );

		if (lsParams->modifyVanishingGradients_flag)
		{ /* avoid vanishing gradients between -GRAD_THRESH and GRAD_THRESH */
			for (unsigned int j = 0; j < M; j++)
			{
				if (gradientVec[j] > 0.0)
				{
					if (gradientSign[j] < 0) /* change of sign detected, minimum detected */
						//gradientThreshold[j] *= 0.9; /* reduce threshold, allow smaller grdients */
						gradientThreshold[j] = 0.0; /* stop tweaking */
					if (gradientVec[j] < gradientThreshold[j])
					{
						gradientVec[j] = gradientThreshold[j];
						if (out_flag)
							fprintf( out, "  , g[%d] fixed to %.2e", j, gradientVec[j] );
					}
					gradientSign[j] = 1; /* remember sign */
				}
				if (gradientVec[j] < 0.0)
				{
					if (gradientSign[j] > 0) /* change of sign detected, minimum detected */
						//gradientThreshold[j] *= 0.9; /* reduce threshold, allow smaller grdients */
						gradientThreshold[j] = 0.0; /* stop tweaking */
					if (gradientVec[j] > -gradientThreshold[j])
					{
						gradientVec[j] = -gradientThreshold[j];
						if (out_flag)
							fprintf( out, "  , g[%d] fixed to %.2e", j, gradientVec[j] );
					}
					gradientSign[j] = -1; /* remember sign */
				}
			}
		}

		if (errno)
		{
			
				char str[511];
				snprintf( str, 510, "\n### %s: ", rtn );
				perror( str );
			
			fprintf( stderr, "     errno = %d", errno );
			fprintf( out, "\n Error in computation (%d), ", errno );
			fprintf( out, "see standard output (console)\n" );
			err = errno;
			errno = 0;
			goto endfunc;
		}

		/* add Levenberg-Marquardt term on main diagonal */
		if (lsParams->LM)
		{
			if (lsState->iterCnt == 1)
			{/* in first loop */
				max_diag = 0.001; // instead of 0.
				for (unsigned int j = 0; j < M; j++)
				{
					/* get maximum value on main diagonal		 */
					if (max_diag < normal[j][j])
					{ /* required only onef for initialization of muLM */
						max_diag = normal[j][j];
					}
				}
				lsState->muLM = lsParams->muFac * max_diag; /* initialize only once */
			}

			if (lsParams->DM)
			{/* use diagonal matrix */
				for (unsigned int j = 0; j < M; j++)
				{
					// normal[j][j] += muLM * diag[j]; 
					//normal[j][j] = fma( lsState->muLM, diag[j], normal[j][j] );
					//normal[j][j] += normal[j][j] * lsState->muLM;
					//normal[j][j] = normal[j][j] * (1. + lsState->muLM);
					normal[j][j] = fma( normal[j][j], lsState->muLM, normal[j][j] );
				}
			}
			else
			{/* use unity/identity matrix */
				for (unsigned int j = 0; j < M; j++)
				{
					normal[j][j] += lsState->muLM;
				}
			}
			/* add Q term ==> true Hessian matrix */
			if (lsParams->trueH)
			{
				if (lsParams->expression_flag)
				{
					for (unsigned int i = 0; i < N; i++)
					{
						x = conditions[i];
						/* damp influence of Q by myHessian */
						residual = observations[i] - te_eval( teExpr ); /* uses global variable x */
						for (unsigned int j = 0; j < M; j++)
						{
							for (unsigned int k = 0; k < M; k++)
							{
								deriv_2nd = f_deriv2TinyExpr( teExpr, j, k, a );
								normal[j][k] -= weights[i] * deriv_2nd * residual;
							}
						}
					}
				}
				else
				{
					for (unsigned int i = 0; i < N; i++)
					{
						/* damp influence of Q by myHessian */
						residual = observations[i] - funct( i, conditions, a );
						for (unsigned int j = 0; j < M; j++)
						{
							for (unsigned int k = 0; k < M; k++)
							{
								deriv_2nd = funct_deriv2( funct, i, j, k, M, conditions, a );
								normal[j][k] -= weights[i] * deriv_2nd * residual;
							}
						}
					}
				}
			}
		}
		else /* Gauss-Newton */
		{
			/* add Q term ==> true Hessian matrix */
			/* including damping */
			if (lsParams->trueH)
			{
				sumOfSquaredResiduals = 0.0;
				for (unsigned int i = 0; i < N; i++)
				{
					/* damp influence of Q by myHessian */
					if (lsParams->expression_flag)
					{
						x = conditions[i]; /* uses global variable x */
						residual = observations[i] - te_eval( teExpr );
					}
					else
					{
						residual = observations[i] - funct( i, conditions, a );
					}
					if (lsParams->GNspecial_flag)
					{
						sumOfSquaredResiduals += residual * residual; /* for adaptation of muHessian */
						residual *= lsState->muHessian; /* damping */
					}

					if (lsParams->expression_flag)
					{
						for (unsigned int j = 0; j < M; j++)
						{
							for (unsigned int k = 0; k < M; k++)
							{
								deriv_2nd = f_deriv2TinyExpr( teExpr, j, k, a );
								normal[j][k] -= weights[i] * deriv_2nd * residual;
							}
						}
					}
					else
					{
						for (unsigned int j = 0; j < M; j++)
						{
							for (unsigned int k = 0; k < M; k++)
							{
								deriv_2nd = funct_deriv2( funct, i, j, k, M, conditions, a );
								normal[j][k] -= weights[i] * deriv_2nd * residual;
							}
						}
					}
				} /* for i */
				if (lsParams->GNspecial_flag)
				{
					lsState->muHessian = min( 1, fabs( sumOfSquaredResiduals0 - sumOfSquaredResiduals ) / sumOfSquaredResiduals0 );
					if (out_flag)
						fprintf( out, ",  muHessian:  %.*e", DBL_DECIMAL_DIG, lsState->muHessian );
					if (lsState->iterCnt == 1) sumOfSquaredResiduals0 = sumOfSquaredResiduals;
				}
			} /* if trueH */
		} /* Gauss-Newton */

		/*
		 * inversion of normal matrix
		 * (cofactor method, LU decomposition, or SVD)
		 */
		err = matrixInversion( M, normal, normal_i, algoMode, out );

		if (err)
		{
			/* ignore error if several trials are carried out to meet target cost */
			if (lsParams->chisqTarget_flag)
			{
				/* fake to trigger re-initialisation of parameters*/
				iterLimitReached_flag = 1;
				/* do-while loop will be continued */
			}
			else
				goto endfunc; 
		}
		else /* do the normal job */
		{
			/* final matrix multiplication to get parameter updates */
			for (unsigned int j = 0; j < M; j++)
			{
				da[j] = 0.;
				for (unsigned int i = 0; i < M; i++)
				{
					//da[j] += normal_i[j][i] * gradientVec[i];		
					da[j] = fma( normal_i[j][i], gradientVec[i], da[j] );
				}
				if (lsParams->GNspecial_flag)
				{
					da[j] *= lsState->daFac; /* for Gauss-Newton: make smaller steps if chisq is not decreasing */
				}
			}

			if (errno)
			{
					char str[511];
					snprintf( str, 510, "\n### %s: ", rtn );
					perror( str );
				fprintf( stderr, "     errno = %d", errno );
				fprintf( out, "\n# Error in computation (%d), ", errno );
				fprintf( out, "see standard output (console)\n" );
				err = errno;
				errno = 0;
				goto endfunc;
			}

			lsState->stopCounter = 0;

			if (out_flag)
			{
				fprintf( out, "\n#  Updates:  " );
			}
			for (unsigned int j = 0; j < M; j++)
			{
				if (out_flag)
					fprintf( out, "da%d=%.*e, ", j + 1, digitsIntermediate, da[j] );
#ifdef ADAPTED_PRECISION
				double eps;
				int RAND_MAX_half = RAND_MAX >> 1;
				if (da[j] < 0.0)
				{ /* compute smallest amount that will change a[j] */
					if (a[j] != 0.0)
						eps = a[j] - nextafter( a[j], INFINITY ); /* a minus next larger  value ==> eps < 0 */
					else /* problems with #DEN */
						eps = 0.0; /* da[] won't change*/
					if (da[j] > eps)
					{ /* da is between eps and zero */
					/* guarantee that da will modify the corresponding a */
						int rndValue = rand();
						if (rndValue > RAND_MAX_half)
						{
							da[j] = eps; /* change in 50% only */
							if (out_flag)
								fprintf( out, "corrected da%d=%.*e\n#\t", j + 1, digitsIntermediate, da[j] );
						}
						else
							if (out_flag)		fprintf( out, "\n#\t" );
						stop_flag++; /* if changes are negligible */
					}
				}
				else
				{ /* compute smallest amount that will change a[j] */
					if (a[j] != 0.0)
						eps = a[j] - nextafter( a[j], -INFINITY ); /* a minus next smaller value ==> eps > 0 */
					else /* problems with #DEN */
						eps = 0.0; /* da[] won't change*/
					if (da[j] < eps)
					{/* da is between zero and eps */
	//#ifdef NOT_USEFUL // could prevent other da to improve result
							/* guarantee that da will modify the corresponding a */
						int rndValue = rand();
						if (rndValue > RAND_MAX_half)
						{
							da[j] = eps; /* change in 50% only */
							if (out_flag)
								fprintf( out, "corrected da%d=%.*e\n#\t", j + 1, digitsIntermediate, da[j] );
						}
						else
							//#endif
							if (out_flag)	fprintf( out, "\n#\t" );
						stop_flag++; /* if changes are negligible */
					}
				}
#else
				/* if changes are negligible */
				if (fabs( da[j] ) < CONVERGENCE_THRESH)
					lsState->stopCounter++;
#endif
			} /* for all j */

			/* adjust parameters */
			memcpy( last_a, a, M * sizeof( double ) ); /* remember last parameters */
			for (unsigned int j = 0; j < M; j++)
			{
				if (lsParams->expression_flag && lsParams->positiveParams)
				{
					/* avoid changes of parameter sign */
					double tmp = a[j] + da[j];
					if (a[j] > 0.0 && (tmp < 0.0))
					{
						da[j] = -0.9 * a[j];
						fprintf( out, "\n# corrected: da%d=%.*e, ", j + 1, digitsIntermediate, da[j] );
					}
					if (a[j] < 0.0 && (tmp > 0.0))
					{
						da[j] = -0.9 * a[j];
						fprintf( out, "\n# corrected: da%d=%.*e, ", j + 1, digitsIntermediate, da[j] );
					}
				}
				/* adjust */
				a[j] += da[j];
			} /* for all j */

			/* compute weighted and squared differences chi-squared */
			lsState->lastChisq = lsState->chisq;
			lsState->chisq = calcChisqNonLin( funct, N, observations, conditions, a,
				weights, lsParams->expression_flag, teExpr );
			if (isnan( lsState->chisq ))
			{
				err = errmsg( ERR_IS_NAN, rtn, "chisq", 0 );
				lsState->stopCounter = M;
			}

			/*-----------------------------------------------
			 * check improvement
			 */
			if (lsState->lastChisq > lsState->chisq)
			{
				/* new result is closer to minimum then last */
				if (lsState->minChisq > lsState->chisq)
				{ /* is also best solution so far */
					lsState->minChisq = lsState->chisq;
					/* remember best parameters */
					for (unsigned int j = 0; j < M; j++)
					{
						min_a[j] = a[j];
					}
				}
				if (lsState->muLM > 2. * DBL_EPSILON)
					lsState->muLM *= muReduceFac; /* decrease Lev-Mar parameter */

				/* Gauss-Newton: adapt damping factor for steps */
				if (lsParams->GNspecial_flag && lsState->daFac < 1.)
				{
					lsState->daFac *= 1.01; /* increase by 1 % */
					if (out_flag) fprintf( out, "\t  daFac: %.*e", digitsIntermediate, lsState->daFac );
				}
				if (lsState->stopCounter == M)
				{
					/* can step be successful if all da have been fallen below the threshold? */
					lsState->stopCounter = 0;
					if (out_flag)
						fprintf( out, "\n#  successful step, although all da are below thresholds !" );
				}
				lsState->chisqWatchDog = initialWatchDogCounter; /* reset after successful step */
			}
			else if (lsState->lastChisq == lsState->chisq)
			{
				lsState->chisqWatchDog--;
				if (out_flag) 
					fprintf( out, "\n#   chisqWatchDog: %d", lsState->chisqWatchDog );
				if (lsState->chisqWatchDog <= 0)
				{
					lsState->stopCounter = M;
					fprintf( out, "\n#   unchanged costs despite significant steps" );
					fprintf( out, "\n#   ==> possibly numerical limit reached" );
				}
				/* assume that parameters have improved despite unchanged cost */
				for (unsigned int j = 0; j < M; j++)
				{
					min_a[j] = a[j];
				}
			}
			else /* costs have increased */
			{
				if (out_flag)
					fprintf( out, "\n#  uphill:" );

				if (lsState->chisq - lsState->lastChisq >= lsParams->uphillAcceptanceThreshold)
				{ /* do not accept uphill step */
					/* try different things */
					unsigned int k;

					if (out_flag)
						fprintf( out, " try one-dimensional step ..." );

					k = 0;
					do
					{
						/* reset parameters */
						for (unsigned int j = 0; j < M; j++)
						{
							a[j] = last_a[j];
						}
						a[k] += da[k]; /* try a single 1D modification */
						lsState->chisq = calcChisqNonLin( funct, N, observations, conditions, a, weights, lsParams->expression_flag, teExpr );
						k++;
						if (k >= M) break;
					} while (lsState->lastChisq < lsState->chisq);

					if (isnan( lsState->chisq ))
					{
						err = errmsg( ERR_IS_NAN, rtn, "chisq", 0 );
						lsState->stopCounter = M; // stop iterations
					}

					if (lsState->lastChisq <= lsState->chisq)
					{ /* no improvement for single 1D steps */
						if (out_flag) fprintf( out, " not successfull" );

						if (lsState->chisqWatchDog < initialWatchDogCounter)
						{ /* continue with watchDog if countdown has already started */
							lsState->chisqWatchDog--;
							if (lsState->chisqWatchDog < initialWatchDogCounter)
							{
								lsState->stopCounter = M;
								fprintf( out, "\n#   unchanged costs despite significant steps" );
								fprintf( out, "\n#   ==> possibly numerical limit reached" );
							}
						}
						/* restore old parameters */
						for (unsigned int j = 0; j < M; j++)
						{
							a[j] = last_a[j];
						}
						lsState->chisq = lsState->lastChisq;

						if (lsParams->LM)
						{ /* Levenberg-Marquardt: adapt damping factor */
							if (lsState->muLM < 1.e+32)
							{
								/* increase Lev-Mar parameter */
								lsState->muLM *= muIncreasFac;
								if (out_flag)
								{
									fprintf( out, ", muLM=%.*e", digitsIntermediate, lsState->muLM );
								}
							}
						}
						else
						{
							if (lsParams->GNspecial_flag)
							{
								lsState->daFac *= 0.99; /* try smaller steps for Gauss-Newton */
								if (out_flag) fprintf( out, "\t  daFac: %.*e", digitsIntermediate, lsState->daFac );
							}
						}
					}
					else /* 1D step was successful */
					{
						if (k > 0 && lsState->stopCounter != M) /* save guard for chisq = nan */
						{
							if (out_flag) fprintf( out, " successfull for a%d", k );
							/* new result is closer to minimum */
							if (lsState->minChisq > lsState->chisq)
							{
								lsState->minChisq = lsState->chisq;
								min_a[k - 1] = a[k - 1];
							}

							if (lsParams->LM) /* modifies anyway */
							{ /* Levenberg-Marquardt */
								if (lsState->muLM < 1.e+16)
								{
									/* increase Lev-Mar parameter */
									lsState->muLM *= muIncreasFac;
									if (out_flag)
										fprintf( out, ", muLM=%.*e", digitsIntermediate, lsState->muLM );
								}
							}
							else
							{
								if (lsState->daFac < 1.0) lsState->daFac *= 1.01; /* try larger steps */
								if (out_flag) fprintf( out, "\t  daFac: %.*e", digitsIntermediate, lsState->daFac );
							}
						}
					}
				} // if upfill 
				else
				{ /* step uphill within acceptance threshold */
					if (out_flag)
						fprintf( out, " do not care ..." );
				}
			} /* else *costs have increased */

			if (out_flag)
			{
				fprintf( out, "\n#  Parameters:  " );
				for (unsigned int j = 0; j < M; j++)
				{
					fprintf( out, "a%d=%.*e, ", j + 1, digitsIntermediate, a[j] );
				}

				fprintf( out, "\n#  chisq: %.*e", OP_DBL_Digs, lsState->chisq );
				if (lsParams->LM) /* Lev-Marquardt method */
				{
					fprintf( out, ", muLM: %*e", digitsIntermediate, lsState->muLM );
					//fprintf( out, ", max_diag: %*e", digitsIntermediate, max_diag );
				}
			}
			/* watch convergence on console */
			fprintf( stderr, "  chisq = %.*e", OP_DBL_Digs, lsState->chisq );

			/* set break condition */
			if (lsState->minChisq == 0.0) lsState->stopCounter = M;

			fflush( stdout );
			if (errno)
			{
					char str[511];
					snprintf( str, 510, "\n### %s: ", rtn );
					perror( str );
				fprintf( stderr, "     errno = %d", errno );
				fprintf( out, "\n Error in computation (%d), ", errno );
				fprintf( out, "see standard output (console)\n" );
				err = errno;
				errno = 0;
				goto endfunc;
			}

			/* convergence reached */
			if (lsState->stopCounter == M) /* all adjustments are negligible */
			{
				converged_flag = 1;

				if (lsState->iterCnt > 3000)
				{
					fprintf( out, "\n#\n#  last updates: " );
					for (unsigned int j = 0; j < M; j++)
					{
						fprintf( out, "da%d=%.14G, ", j + 1, da[j] );
					}
				}
				{
					fprintf( out, "\n#\n#  convergence after %d iterations",
						lsState->iterCnt );
					fprintf( stderr, "\n#\n#  convergence after %d iterations",
						lsState->iterCnt );
					//break; /* stop iterative process */
				}
			} /* if (stopCounter == M)	*/
			else
				lsState->iterCnt++;
			if (lsState->iterCnt >= iter_max)
			{
				iterLimitReached_flag = 1;
				/* take best parameter */
				for (unsigned int j = 0; j < M; j++)
				{
					a[j] = min_a[j];
				}
				fprintf( stderr,
					"\n\n     no convergence after %d iterations",
					lsState->iterCnt );
				fprintf( stderr, "\n     chisq: %e", lsState->chisq );
				fprintf( out,
					"\n#\n#     no convergence after %d iterations",
					lsState->iterCnt );
				fprintf( out, "\n#     chisq: %e", lsState->chisq );
				//break;
			}
			if (lsParams->uphillAcceptanceThreshold > 1e-20)
			{
				/* decrease acceptance by 10% */
				//lsParams->uphillAcceptanceThreshold = lsParams->uphillAcceptanceThreshold *= 0.9;
				/* decrease acceptance by 10%, but it may not be larger than chi squared */
				lsParams->uphillAcceptanceThreshold =
					min( lsState->chisq, lsParams->uphillAcceptanceThreshold *= 0.9 );
			}
			else lsParams->uphillAcceptanceThreshold = 0.0; /* avoid zero steps if chisq==last_chisq */
			if (out_flag)
				fprintf( out,
					",  uphillAcceptanceThreshold: %.*e", digitsIntermediate,
					lsParams->uphillAcceptanceThreshold );

		}
		/* ---------
		 * check target for costs
		 */
		if (lsParams->chisqTarget_flag)
		{
			if (lsState->bestChisq > lsState->minChisq)
			{
				lsState->bestChisq = lsState->minChisq;
				for (unsigned int j = 0; j < M; j++)
				{
					best_a[j] = min_a[j];
				}
			}
			if (converged_flag)
			{
				fprintf( out, "\n#\n#  Iteration of least squares: %d, convergence",
					lsState->iterCnt );
				fprintf( out, "\n#  check target (%e) ", chisqTarget );
				/* check, whether maximum target error is reached */
				if (lsState->chisq <= chisqTarget)
				{
					/* yes, we can stop the optimisation */
					fprintf( stderr, " target met, best chisq: %.*e\n", OP_DBL_Digs, lsState->bestChisq );
					fprintf( out, " target met, best chisq: %.*e\n", OP_DBL_Digs, lsState->bestChisq );
					targedReached_flag = 1;
				}
				else
				{
					if (closestToChisqTarget > lsState->bestChisq)
					{
						closestToChisqTarget = lsState->bestChisq;
						fprintf( stderr, " possibly local minimum, best chisq: %.*e\n", OP_DBL_Digs, closestToChisqTarget );
					}
					else
						fprintf( stderr, " possibly local minimum\n" );
				} /* else if (iter_cnt < iter_max-1)	*/
			}
			if (iterLimitReached_flag)
			{
				if (closestToChisqTarget != DBL_MAX)
				{
					lsState->chisq = closestToChisqTarget;
				}
				fprintf( stderr,
					"\n\n  target value for chisq not reached after %d iterations",
					lsState->iterCnt );
				fprintf( stderr, "\n  target chisq: %e", chisqTarget );
				fprintf( stderr, "\n    best chisq: %e\n", lsState->bestChisq );
				fprintf( out,
					"\n#\n#   target value for chisq not reached after %d iterations",
					lsState->iterCnt );
				fprintf( out, "\n#  target chisq: %e", chisqTarget );
				fprintf( out, "\n#    best chisq: %e\n", lsState->bestChisq );
			}

			if (converged_flag || iterLimitReached_flag)
			{
				/* set flags and randomise parameters */
				fprintf( out, "\n# currently  best parameters:  " );
				for (unsigned int j = 0; j < M; j++)
				{
					fprintf( out, "a%d=%.*e, ", j + 1, digitsIntermediate, best_a[j] );
				}

				if (targedReached_flag)
				{
					fprintf( out, "\n# success in trial %d/%d", chisqTargetCounter + 1, NUMTARGETTRIALS );
				}
				else
				{
					/* lets re-initialise the parameters  */
					fprintf( out, "\n#  re-initialise parameters" );
					fprintf( out, "\n# randomised parameters:  " );
					for (unsigned int j = 0; j < M; j++)
					{
						float z;
						/* -1...+1 */
						z = (2.F * (float)random() / (float)RAND_MAX - 1.F);
						a[j] = initial_a[j] * z * 2; /* max. change 200% of initial values */
						fprintf( out, "a%d=%.*e, ", j + 1, digitsIntermediate, a[j] );
					}
					fprintf( out, "\n" );
					chisqTargetCounter++;
					if (chisqTargetCounter < NUMTARGETTRIALS)
					{/* try again */
						fprintf( out, "\n# trial %d/%d", chisqTargetCounter + 1, NUMTARGETTRIALS );
						fprintf( stderr, "\n# trial %d/%d\n", chisqTargetCounter + 1, NUMTARGETTRIALS );
						lsState->stopCounter = 0;
						lsState->iterCnt = 0;
						converged_flag = 0;
						iterLimitReached_flag = 0;
						closestToChisqTarget = lsState->chisq = DBL_MAX; /* will be copied to lastChisq */
					}
				}
			} /* if (converged_flag || iterLimitReached_flag) */
		} /* if (lsParams->chisqTarget)	*/
	} while (!iterLimitReached_flag && !converged_flag && !targedReached_flag);

	/* final output */

	/*
	 * determination of covariance matrix
	 * if  not yet singular
	 */
	if (err != ERR_IS_SINGULAR)
	{
		/* compute normal w/o muLM !!!! */
		for (unsigned int j = 0; j < M; j++)
		{
			for (unsigned int i = 0; i < M; i++)
			{
				normal[j][i] = 0.;
				for (unsigned int n = 0; n < N; n++)
				{
					//normal[j][i] += jacob[n][j] * jacob[n][i] * weights[n];
					normal[j][i] = fma( jacob[n][j], jacob[n][i] * weights[n], normal[j][i] );
				}
			}
		}

		/* inversion
		 * for ill-conditioned problems (e.g. polynomials of high
		 * order)	the inversion of the normal matrix might fail
		 *
		 */
		err = matrixInversion( M, normal, normal_i, algoMode, out );
		/* take best parameter */
		for (unsigned int j = 0; j < M; j++)
		{
			a[j] = min_a[j];
		}
		/* if algorithm converged, then ignor possible errors
		 * this avoids the output of 'failed/x/' */
		if (converged_flag) err = 0;
	}
endfunc:
	*iterations = lsState->iterCnt;
	fprintf( out,
		"\n#  -- %s - end  ------------------------------", rtn );

	free_vector( &initial_a );
	free_vector( &best_a );
	free_vector( &last_a );
	free_vector( &min_a );
	free_vector( &da );
	free_vector( &gradientVec );
	free_vector( &datac );
	free_matrix( &normal );
	free_matrix( &cofac );

	return err;
}
