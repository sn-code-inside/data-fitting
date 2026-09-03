/*****************************************************************
 *
 * File........:	lsNonlinear.c
 * Function....:	least squares with alternative matrix inversion
 * Author......:	Tilo Strutz
 * last changes:	04.08.2026
 *
 * LICENCE DETAILS: see software manual
 *	free academic use
 *  cite source as
 * "Strutz, T.: Data Fitting and Uncertainty. 2nd edition,
 *   Springer Vieweg, 2016"
 *
*****************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <assert.h>
#include "errmsg.h"
#include "matrix_utils.h"
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

long ITERAT_MAX;

extern double del, delFac;/* for numerical derivation, used in functions.c f_deriv() */

extern double x; /* declared in fitting.c */

/*C2TEXfilter_begin*/
//#define QDEBUG
/*C2TEXfilter_end*/

/*---------------------------------------------------------------
*	lsNonlinear()
*
*--------------------------------------------------------------*/
int
lsNonlinear( double (*funct) (int, double*, double*),
	double (*funct_deriv) (double(*)(int, double*, double*),
		int, int, int, double*, double*),
	double (*funct_deriv2) (double(*)(int, double*, double*),
		int, int, int, int, double*, double*),
	te_expr* teExpr, te_variable* teVars,
	int type,
	unsigned int N, unsigned int M, double* observations, double* conditions, double** jacob,
	double* weights,
	double* a,
	int algoMode, LS_PARAMS* lsParams,
	double chisqTarget, double** covar, int* iterations, FILE* out,
	int digitsIntermediate )
{
	char* rtn = "lsNonlinear";
	int err = 0;
	int Nfree;
	int out_flag;/* if number of iteration is to high avoid too many outputs */
	int iter_max = ITERAT_MAX;	/* maximum number of iterations */
	int chisqTargetCounter = 0;
	int converged_flag = 0;
	int targedReached_flag = 0;
	int iterLimitReached_flag = 0;
	double** cofac = NULL;	/* cofactor matrix for matrix
													inversion */
	double** normal = NULL;			/* N = J^(T) * W * J	*/
	//double** normal_i = NULL;		/* inverse of N				*/
	double* tmpvec = NULL;			/* J^(T) * W * r			*/
	double* datac = NULL;	/* calculated values based on parameters
												*/
	double* da = NULL;		/* parameter update deltasq a */
	double* deltasq = NULL;	/* = w * [observations - f(x|a) ] ^2                XXXX  can be a temporary value, can'nt it? */
	double residual, deriv_2nd;
	double* last_a = NULL; /* remember last result */
	double* min_a = NULL; /* remember best result */
	double closestToChisqTarget = INFINITY;
	double max_diag;	/* maximum value of Njj */
	double diag[M_MAX] = { 0 }; /* vector of diagonal of Normal matrix	*/
	double muReduceFac = 0.6666666666666666;
	double muIncreasFac = 4.0;/* should not be a multiple of decreasing factor */

	LS_STATE lsStateObj = { 0 }, * lsState;

	lsState = &lsStateObj;

	fprintf( out,
		"\n#  -- %s - start  ------------------------------", rtn );
	lsState->lastChisq = lsState->minChisq = lsState->chisq = DBL_MAX;
	/*
	* allocate memory
	*/

	/* normal matrix N = J^(T) * W * J, its inverse */
	normal = matrix( M, M );
	//normal_i = matrix( M, M );

	cofac = matrix( M, M );	/* cofactor matrix */
	tmpvec = vector( M );		/* container for J^(T) * W * G	*/
	datac = vector( N );		/* calculated data using f(x|a) */
	deltasq = vector( N );		/* remaining differences	*/
	da = vector( M );				/* model parameter update */
	min_a = vector( M );			/* remember best parameter set */
	last_a = vector( M );			/* remember last parameter set */

	lsState->daFac = 1.0;
	lsState->iterCnt = 1;

	/* initialise minimum (best) values */
	for (unsigned int j = 0; j < M; j++)
	{
		min_a[j] = last_a[j] = a[j]; /* initial estimtes */
	}

	if (lsParams->LM)
	{  /* Levenberg-Marquard */
		fprintf( out, "\n#  initial mu  factor: %lf", lsParams->muFac );
		fprintf( out, "\n#  mu decrease factor: %lf", muReduceFac );
		fprintf( out, "\n#  mu increase factor: %lf", muIncreasFac );
	}

	/* compute initial chi squared */
	lsState->chisq = calcChisqNonLin( funct, N, observations, conditions, a, weights, &Nfree, lsParams->expression_flag, teExpr );
	lsState->lastChisq = lsState->minChisq = lsState->chisq;
	fprintf( out, "\n#\n#  initial chisq=%.*e", OP_DBL_Digs, lsState->minChisq );

	if (lsParams->LM)
		lsParams->uphillAcceptanceThreshold = 0.0; /* never uphill */
	else
		lsParams->uphillAcceptanceThreshold = lsState->chisq; /* allow 100% increase at the beginning */

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

		/*C2TEXfilter_begin*/
#ifdef QDEBUG
		fprintf( out,
			"\n#\n#== Obs ==  Weight == Jacobian =================" );
		for (i = 0; i < N; i++)
		{
			fprintf( out, "\n# %8.1f %8.6f", obs[i], weights[i] );
			for (j = 0; j < M; j++)
			{
				fprintf( out, " %8.5f", jacob[i][j] );
			}
		}
#endif
		/*C2TEXfilter_end*/

		/* get Jacobian matrix  */
		err = getJacobian( lsParams->expression_flag, funct, funct_deriv,
			teExpr, teVars,
			N, M, conditions, jacob, a,	out );

#ifdef OBSOLETE
		/* update Jacobian matrix for nonlinear models */
		/* am Anfang der Schleife statt am Ende */
		if (lsParams->expression_flag)
		{
			for (unsigned int i = 0; i < N; i++)
			{
				x = conditions[i];
				for (unsigned int j = 0; j < M; j++)
				{
					jacob[i][j] = f_derivTinyExpr( teExpr, teVars, j, M, a );
					if (!IsFiniteNumber( jacob[i][j] ))
					{
						fprintf( stderr,
							"\n# %s\n####  Divergence of approximation: %e", rtn, jacob[i][j] );
						fprintf( out,
							"\n# %s\n#  Divergence of approximation: %e", rtn, jacob[i][j] );
						if (errno)
						{
							perror( "\n### " );
							fprintf( stderr, "     errno = %d", errno );
							fprintf( out, "\terrno = %d", errno );
							fprintf( out, "\tsee standard output (console)\n" );
							err = errno;
							goto endfunc;
						}
						err = ERR_IS_INFINITE;
						goto endfunc;
					}
				}
			}
		}
		else
		{
			/* in-built model function */
			for (unsigned int i = 0; i < N; i++)
			{
				for (unsigned int j = 0; j < M; j++)
				{
					jacob[i][j] = funct_deriv( funct, i, j, M, conditions, a );
					if (errno)
					{
						perror( "\n### " );
						fprintf( stderr, "     errno = %d", errno );
						fprintf( out, "\terrno = %d", errno );
						fprintf( out, "\tsee standard output (console)\n" );
						err = errno;
						if (!IsFiniteNumber( jacob[i][j] ) || jacob[i][j] == NAN || jacob[i][j] == INFINITY)
						{
							fprintf( stderr,
								"\n# %s\n####  Divergence of approximation at (i,j) = (%d,%d)", rtn, i, j );
							fprintf( out,
								"\n# %s\n#  Divergence of approximation at (i,j) = (%d,%d)", rtn, i, j );
							/* retry*/
							goto endfunc;
						}
					}
				}
			}
		}
#endif
		/*
		* calculate normal matrix N
		* N = J^(T) * W * J
		*/
		err = computeNormalMatrix( N, M, jacob, normal, weights );

		/* K = J^(T) * W * r */
		err = computeJWr( funct, N, M, observations, conditions, jacob,
			weights, a, tmpvec, lsParams->expression_flag, teExpr, teVars );

#ifdef OBSOLETE
		if (errno)
		{
			perror( "\n### " );
			fprintf( stderr, "     errno = %d", errno );
			fprintf( out, "\n Error in computation (%d), ", errno );
			fprintf( out, "see standard output (console)\n" );
			err = errno;
			errno = 0;
			goto endfunc;
		}
		max_diag = 0.001; // insead of 0.
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
				/* overflow test	*/
				if (!IsFiniteNumber( normal[j][i] ))
				{
					err = errmsg( ERR_IS_INFINITE, rtn, "normal", 0 );
					goto endfunc;
				}
			}
			/* only for nonlinear models of importance:
			 * get maximum value on main diagonal
			 */
			diag[j] = normal[j][j]; /* diag[] is later used */
			if (max_diag < normal[j][j])
			{
				max_diag = normal[j][j];
			}
			if (errno)
			{
				perror( "\n### " );
				fprintf( stderr, "     errno = %d", errno );
				fprintf( out, "\n Error in computation (%d), ", errno );
				fprintf( out, "see standard output (console)\n" );
				err = errno;
				errno = 0;
				goto endfunc;
			}
		}
#endif

		/*C2TEXfilter_begin*/
#ifdef QDEBUG
		fprintf( out, "\n#\n#== normal matrix =================" );
		for (unsigned int i = 0; i < M; i++)
		{
			fprintf( out, "\n# " );
			for (unsigned int j = 0; j < M; j++)
			{
				fprintf( out, " %12.2f", normal[i][j] );
			}
		}
		fflush( out );
#endif
		/*C2TEXfilter_end*/
		if (errno)
		{
			perror( "\n### " );
			fprintf( stderr, "     errno = %d", errno );
			fprintf( out, "\n Error in computation (%d), ", errno );
			fprintf( out, "see standard output (console)\n" );
			err = errno;
			errno = 0;
			goto endfunc;
		}

#ifdef OBSOLETE
		/* tmpvec = J^(T) * W * r */
		/* r contains residuals */
		if (lsParams->expression_flag)
		{
			for (unsigned int j = 0; j < M; j++)
			{
				tmpvec[j] = 0.;
				for (unsigned int i = 0; i < N; i++)
				{
					// residual = observations[i] - funct( i, conditions, a );
					x = conditions[i];
					residual = observations[i] - te_eval( teExpr );
					//tmpvec[j] += jacob[i][j] * residual * weights[i];
					tmpvec[j] = fma( jacob[i][j], residual * weights[i], tmpvec[j] );
				}
			}
		}
		else
		{
			for (unsigned int j = 0; j < M; j++)
			{
				tmpvec[j] = 0.;
				for (unsigned int i = 0; i < N; i++)
				{
					residual = observations[i] - funct( i, conditions, a );
					// tmpvec[j] += jacob[i][j] * residual * weights[i];
					tmpvec[j] = fma( jacob[i][j], residual * weights[i], tmpvec[j] );
				}
			}
		}
		if (errno)
		{
			perror( "\n### " );
			fprintf( stderr, "     errno = %d", errno );
			fprintf( out, "\n Error in computation (%d), ", errno );
			fprintf( out, "see standard output (console)\n" );
			err = errno;
			errno = 0;
			goto endfunc;
		}
#endif


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
		}
		/* add Q term ==> true Hessian matrix */
		if (lsParams->trueH)
		{
			for (unsigned int i = 0; i < N; i++)
			{
				residual = observations[i] - funct( i, conditions, a );
				for (unsigned int j = 0; j < M; j++)
				{
					for (unsigned int k = 0; k < M; k++)
					{
						deriv_2nd = funct_deriv2( funct, j, k, M, i, conditions, a );
						normal[j][k] += deriv_2nd * residual;
					}
				}
			}
		}

		/*
		 * inversion of normal matrix
		 * (cofactor method, LU decomposition, or SVD)
		 *
		 */
		err = matrixInversion( N, M, normal, covar, algoMode, out );

/*C2TEXfilter_begin*/
#ifdef QDEBUG
		fprintf( out,
			"\n#\n#== inverse normal matrix =================" );
		for (unsigned int i = 0; i < M; i++)
		{
			fprintf( out, "\n# " );
			for (unsigned int j = 0; j < M; j++)
			{
				fprintf( out, " %12.2f", covar[i][j] );
			}
		}
		fflush( out );
#endif
		/*C2TEXfilter_end*/


		/* final matrix multiplication to get parameter updates */
		int delModified_flag = 0;
		for (unsigned int j = 0; j < M; j++)
		{
			da[j] = 0.;
			for (unsigned int i = 0; i < M; i++)
			{
				//da[j] += covar[j][i] * tmpvec[i];		
				da[j] = fma( covar[j][i], tmpvec[i], da[j] );
			}
			da[j] *= lsState->daFac; /* for Gauss-Newton: make smaller steps if chisq is not decreasing */

			if (lsParams->delModification_flag && !delModified_flag &&
				fabs( da[j] ) < TOL * 10 && del < 0.001)
			{ /* increase numerical stability */
				delFac = delFac * 10;
				//if (delFac> 0.001) delFac = delFac * 0.9;
				fprintf( out, "\n#  delFac:  %.*e", DBL_DECIMAL_DIG, delFac );
				delModified_flag = 1;
			}
		}

		if (errno)
		{
			perror( "\n### " );
			fprintf( stderr, "     errno = %d", errno );
			fprintf( out, "\n Error in computation (%d), ", errno );
			fprintf( out, "see standard output (console)\n" );
			err = errno;
			errno = 0;
			goto endfunc;
		}

		lsState->stopCounter = 0;

		if (out_flag)
		{
			fprintf( out, "\n#\n#  Iteration of least squares: %d",
				lsState->iterCnt );
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
			if (fabs( da[j] ) < TOL)
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
		Nfree = -(signed)M;	/* reduce by number of parameters */
		lsState->lastChisq = lsState->chisq;
		lsState->chisq = calcChisqNonLin( funct, N, observations, conditions, a, weights, &Nfree, lsParams->expression_flag, teExpr );
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
			if (!lsParams->LM && lsState->daFac < 1.)
			{
				lsState->daFac *= 1.01; /* increase by 1 % */
				if (out_flag) fprintf( out, "\t  daFac: %.*e", digitsIntermediate, lsState->daFac );
			}
			if (lsState->stopCounter == M)
			{
				lsState->stopCounter = 0; /* can step successful if all da have been fallen below the threshold? */
				if (out_flag)
					fprintf( out, "\n#  successful step  although all da are below thresholds !" );
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
					lsState->chisq = calcChisqNonLin( funct, N, observations, conditions, a, weights, &Nfree, lsParams->expression_flag, teExpr );
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
						lsState->daFac *= 0.99; /* try smaller steps for Gauss-Newton */
						if (out_flag) fprintf( out, "\t  daFac: %.*e", digitsIntermediate, lsState->daFac );
						/* attention: resetting to min_a and randomization ist not useful */
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
			perror( "\n### " );
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
			if (lsState->iterCnt > 3000)
			{
				fprintf( out, "\n#\n#  last updates: " );
				for (unsigned int j = 0; j < M; j++)
				{
					fprintf( out, "da%d=%.14G, ", j + 1, da[j] );
				}
			}
			if (lsParams->chisqTarget)
			{
				fprintf( out, "\n#\n#  Iteration of least squares: %d, convergence",
					lsState->iterCnt );
				fprintf( out, "\n#  check target (%e) ", chisqTarget );
				/* check, whether maximum target error is reached */
				if (lsState->chisq <= chisqTarget)
				{
					/* yes, we can stop the optimisation */
					fprintf( stderr, " target met, best chisq: %.*e\n", OP_DBL_Digs, lsState->chisq );
					fprintf( out, " target met, best chisq: %.*e\n", OP_DBL_Digs, lsState->chisq );
					break;
				}
				else 
				{
					if (closestToChisqTarget > lsState->chisq)
					{
						closestToChisqTarget = lsState->chisq;
						fprintf( stderr, " possibly local minimum, best chisq: %.*e\n", OP_DBL_Digs, closestToChisqTarget );
					}
					else
						fprintf( stderr, " possibly local minimum" );


					/* lets re-initialise the parameters  */
					fprintf( out, "\n#  re-initialise parameters" );

					/* set flags and randomise parameters */
					for (unsigned int j = 0; j < M; j++)
					{
						float z;
						/* -1...+1 */
						z = (2.F * (float)random() / (float)RAND_MAX - 1.F);
						a[j] += min_a[j] * z / 2.; /* max. change 50% */
					}
					chisqTargetCounter++;
					if (chisqTargetCounter < 10)
					{/* try again */
						lsState->stopCounter = 0;
						lsState->iterCnt = 0;
					}

					lsState->minChisq = INFINITY; /* restart monitoring costs */
				} /* else if (iter_cnt < iter_max-1)	*/
			} /* if (lsParams->chisqTarget)	*/
			else
			{
				fprintf( out, "\n#\n#  convergence after %d iterations",
					lsState->iterCnt );
				fprintf( stderr, "\n#\n#  convergence after %d iterations",
					lsState->iterCnt );
				break; /* stop iterative process */
			}
		} /* if (stopCounter == M)	*/

		lsState->iterCnt++;
		if (lsState->iterCnt >= iter_max)
		{
			/* take best parameter */
			for (unsigned int j = 0; j < M; j++)
			{
				a[j] = min_a[j];
			}
			if (lsParams->chisqTarget)
			{
				if (closestToChisqTarget != INFINITY)
				{
					lsState->chisq = closestToChisqTarget;
				}
				fprintf( stderr,
					"\n\n  target value for chisq not reached after %d iterations",
					lsState->iterCnt );
				fprintf( stderr, "\n  target chisq: %e", chisqTarget );
				fprintf( stderr, "\n    best chisq: %e", lsState->chisq );
				fprintf( out,
					"\n#\n#   target value for chisq not reached after %d iterations",
					lsState->iterCnt );
				fprintf( out, "\n#  target chisq: %e", chisqTarget );
				fprintf( out, "\n#    best chisq: %e", lsState->chisq );
			}
			else
			{
				fprintf( stderr,
					"\n\n     no convergence after %d iterations",
					lsState->iterCnt );
				fprintf( stderr, "\n     chisq: %e", lsState->chisq );
				fprintf( out,
					"\n#\n#     no convergence after %d iterations",
					lsState->iterCnt );
				fprintf( out, "\n#     chisq: %e", lsState->chisq );
			}
			break;
		}
		if (lsParams->uphillAcceptanceThreshold > 1e-20)
			lsParams->uphillAcceptanceThreshold *= 0.9; /* decrease acceptance by 10% */
		else lsParams->uphillAcceptanceThreshold = 0.0; /* avoid zero steps if chisq==last_chisq */
		if (out_flag)
			fprintf( out,
				",  uphillAcceptanceThreshold: %.*e", digitsIntermediate,
				lsParams->uphillAcceptanceThreshold );
	} while (1);

	/*
	 * determination of covariance matrix
	 * if  not yet singular
	 */
	if (err != ERR_IS_SINGULAR )
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
		err = matrixInversion( N, M, normal, covar, algoMode, out );
		if (err)
		{
			/* take best parameter */
			for (unsigned int j = 0; j < M; j++)
			{
				a[j] = min_a[j];
			}
		}
	}
endfunc:
	*iterations = lsState->iterCnt;
	fprintf( out,
		"\n#  -- %s - end  ------------------------------", rtn );

	free_vector( &last_a );
	free_vector( &min_a );
	free_vector( &da );
	free_vector( &tmpvec );
	free_vector( &datac );
	free_vector( &deltasq );
	free_matrix( &normal );
	//free_matrix( &normal_i );
	free_matrix( &cofac );

	return err;
}
