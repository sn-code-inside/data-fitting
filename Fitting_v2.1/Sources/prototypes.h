/*****************************************************************
 *
 * File........:	prototypes.h
 * Function....:	proto typing for different functions
 * Author......:	Tilo Strutz
 * last changes:	27.01.2010, 30.3.2011, 05.11.2025, 07.01.2026
 *                02.03.2026
 * LICENCE DETAILS: see software manual
 *	free academic use
 *  cite source as 
 * "Strutz, T.: Data Fitting and Uncertainty. 3rd edition,
 *   Springer Nature, 2026"
 *
 *****************************************************************/
#ifndef PROTO_H
#define PROTO_H

#include "ls.h"
#include "tinyexpr.h"

int svdInversion( unsigned int N, double **normal, double **normal_i, 
									FILE *out);
int IsFiniteNumber(double x);

/* parsing of command-line parameters */
char* get_nth_field( char *line, int n);
int is_data_line( char *line, unsigned int N);

/* matrix inversion */
int	singvaldec( double **a, unsigned int N, unsigned int M, double w[], double **v);
void	backsub_LU( double **lu, unsigned int N, unsigned int *indx, double back[]);
int	decomp_LU( double **normal, unsigned int M, unsigned int *indx, int *s);
void heap_sort_d_(unsigned long N, double ra[], long idx[]);
void heap_sort_d(unsigned long N, double ra[]);

/* estimation of weights */
void est_weights1( unsigned int N, double *deltasq,
					double *weights, FILE *out);
void est_weights2( unsigned int N, double *cond, double *obs,
							double *weights, unsigned int obs_per_bin, FILE *out);
int outlier_detection1( unsigned int N, double sigma_y, double *deltasq,
					double *weights, double nu, 
     FILE *out);
int outlier_detection2( unsigned int N, double *deltasq,
					double *weights, FILE *out);
int outlier_detection3( unsigned int N, double *deltasq,
					double *weights, double nu, FILE *out);

int
ransac( double (*funct) (unsigned int,double*,double*),
	 double (*funct_deriv) (double(*)(unsigned int,double*,double*),
		 unsigned int, unsigned int, unsigned int,double*,double*),
	 double (*funct_deriv2) (double(*)(unsigned int,double*,double*),
		 unsigned int, unsigned int, unsigned int, unsigned int,double*,double*),
	 te_expr* teExpr, unsigned int N, unsigned int M, 
	double *obs, double *cond, double **jacob,
	 double *weights, double *a,
	 int algo_mode, LS_PARAMS *ls_flag,
	 double chisqTarget, double **covar, FILE *out,
	 double *deviates_abs,
	unsigned int cond_dim,
	 unsigned int obs_dim, int digitsIntermediate );

	 void
ls_straightline( 
	unsigned int N, /* number of entries */
	double cond[], /* vector of conditions */
	double obs[],		/* vector of observations */
	double a[]			/* container for parameters to be estimated */
	);

	 double conditionScaling( unsigned int N, unsigned int condDim, double* conditions,
		 unsigned int type, FILE* out );
	 void condition_unscaling( unsigned int N, double* conditions,
		 unsigned int condDim, unsigned int M, double* a, double scalingFactor,
		 unsigned int type, FILE* out );
	 void specialPreparations( unsigned int N, double* conditions,
		 double* observations, double* weights, unsigned int type);
	 void specialPostprocessing( unsigned int M, double* a,
		 unsigned int type, FILE* out );
	 void checkParameterUncertainties( unsigned int M, double* a,
		 double** covar, FILE* out, char* outname );

#endif
