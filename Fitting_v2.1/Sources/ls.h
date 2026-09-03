/*****************************************************************
 *
 * File........: lsLinear.h
 * Function....: header for lsLinear.c
 * Author......: Tilo Strutz
 * last changes: 04.08.2026
 * LICENCE DETAILS: see software manual
 *	free academic use
 *  cite source as
 * "Strutz, T.: Data Fitting and Uncertainty. 3rd edition,
 *   Springer Nature, 2026"
 *
 *****************************************************************/
#include "tinyexpr.h"

#ifndef LS_H
#define LS_H

typedef struct
{
	int linear;			/* linear model */
	int svd;				/* special computation for linear models */
	int LM;					/* 0 .. Gaus-Newton, 1 .. Levenberg-Marquardt	*/
	int DM;					/* 0 .. identity, 1 .. diag(JWJ) matrix for Levenberg-Marquardt	*/
	int chisqTarget_flag;	/* indicates that 'chisqTarget' was set	*/
	int trueH;				/* use true Hessian matrix	*/
	int positiveParams; /* for -y : keeps parameters positive */
	int expression_flag;/* indicates use of user-provided model function */
	int GNspecial_flag; /* for Gauss-Newton */
	int modifyVanishingGradients_flag; /* avoids g_j == 0.0 */
	double muFac;/* muLM = muFac * max_diag */
	double uphillAcceptanceThreshold;
} LS_PARAMS;

typedef struct
{
	double chisq;			/* current chi squared, cost function */
	double lastChisq;		/* chi squared of previous iteration */
	double minChisq;	/* minimum chi squared in current trial */
	double bestChisq;	/* minimum chi squared of all trials */
	double muLM;			/* damping factor of Levenberg-Marquardt */
	double muHessian;	/* damping factor of Gauss-Newton, when using true Hessian */
	double daFac;			/* damping for Gauss-Newton */
	unsigned int stopCounter;    /* indicates that iteration must be stopped */
	unsigned int iterCnt;    /* counts  iterations */
	int chisqWatchDog; /* watch change in costs counter */
} LS_STATE;

/* least squares routine */
int
lsLinear( double (*funct) (unsigned int, double*, double*),
	double (*funct_deriv) (double(*)(unsigned int, double*, double*),
		unsigned int, unsigned int, unsigned int, double*, double*),
	te_expr* teExpr, 
	unsigned int N, unsigned int M, double* obs, double* cond, double** jacob,
	double* weights, 
	double* a,
	int algo_mode, LS_PARAMS* ls_flag,
	double** covar, FILE* out);

int
lsNonlinear( double (*funct) (unsigned int, double*, double*),
	double (*funct_deriv) (double(*)(unsigned int, double*, double*),
		unsigned int, unsigned int, unsigned int, double*, double*),
	double (*funct_deriv2) (double(*)(unsigned int, double*, double*),
		unsigned int, unsigned int, unsigned int, unsigned int, double*, double*),
	te_expr* teExpr, 
	unsigned int N, unsigned int M, double* obs, double* cond, double** jacob,
	double* weights,
	double* a,
	int algo_mode, LS_PARAMS* lsParams,
	double chisqTarget, double** covar, int* iterations, FILE* out,
	int digitsIntermediate );

double calcChisqNonLin( double (*funct) (unsigned int, double*, double*),
	unsigned int N, double* obs, double* cond, double* a,
	double* weights, int expression_flag, te_expr* teExpr );

int
solveLinear( unsigned int N, unsigned int M, double* obs, double* weights,
	double** jacob, double** covar, double* a,
	FILE* out );


int getJacobian( int expression_flag,
	double (*funct) (unsigned int, double*, double*),
	double (*funct_deriv) (double(*)(unsigned int, double*, double*),
		unsigned int, unsigned int, unsigned int, double*, double*),
	te_expr* teExpr,
	unsigned int N, unsigned int M, double* conditions, double** jacob,
	double* a,
	FILE* out );

int computeNormalMatrix(
	unsigned int N, unsigned int M,
	double** jacob,
	double** normal, double* weights );

int computeJWy(
	unsigned int N, unsigned int M, double* observations,
	double** jacob, double* weights, double* tmpvec );

int computeJWr(
	double (*funct) (unsigned int, double*, double*),
	unsigned int N, unsigned int M, double* observations,
	double* conditions, double** jacob, double* weights,
	double* a, double* tmpvec, int expression_flag,
	te_expr* teExpr);

int matrixInversion(
	unsigned int M,
	double** normal, double** normal_i, int algoMode,
	FILE* out );

#endif