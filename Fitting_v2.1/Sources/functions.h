/*****************************************************************
 *
 * File........: functions.h
 * Function....: proto typing for functions.c
 * Author......: Tilo Strutz
 * last changes: 25.09.2009, 06.11.2009, 18.02.2010, 03.01.2011
 *		25.02.2021, 05.11.2025, 11.05.2026
 * LICENCE DETAILS: see software manual
 *	free academic use
 *  cite source as 
 * "Strutz, T.: Data Fitting and Uncertainty. 2nd edition,
 *   Springer Vieweg,2026"
 *
 *****************************************************************/

#ifndef FUNCT_H
#define FUNCT_H
#include "tinyexpr.h"

extern double delFac; /* for numerical derivation, declared in lsLinear.c */
extern double fourthRootDBLEPS, thirdRootDBLEPS; /* declared in fitting.c */

/* linear functions */
double fconstant_deriv( double (*funct)(unsigned int,double*,double*), 
							unsigned int i, unsigned int j, unsigned int M, double *cond,	double *a);
double flin_deriv( double (*funct)(unsigned int,double*,double*), 
							unsigned int i, unsigned int j, unsigned int M, double *cond, double *a);
double flin2_deriv( double (*funct)(unsigned int,double*,double*), 
							unsigned int i, unsigned int j, unsigned int M, double *cond, double *a);
double fcosine_deriv( double (*funct)(unsigned int,double*,double*), 
							unsigned int i, unsigned int j, unsigned int M, double *cond,	double *a);
double fpolynom2_deriv( double (*funct)(unsigned int,double*,double*), 
							unsigned int i, unsigned int j, unsigned int M, double *cond, double *a);
double fpolynom3_deriv( double (*funct)(unsigned int,double*,double*), 
							unsigned int i, unsigned int j, unsigned int M, double *cond, double *a);
double fpolynomial_deriv( double (*funct)(unsigned int,double*,double*), 
							unsigned int i, unsigned int j, unsigned int M, double *cond, double *a);
double fquadsurface_deriv( double (*funct)(unsigned int,double*,double*), 
							unsigned int i, unsigned int j, unsigned int M, double *cond, double *a);

/* nonlinear functions */
double fpolynomial( unsigned int i, double *cond, double *a);
int init_polynomial( unsigned int N, double *obs, double *cond,
							double *a, unsigned char *a_flag, FILE *logfile);

double fcosineLin_deriv( double(*funct)(unsigned int, double*, double*),
	unsigned int i, unsigned int j, unsigned int M, double* cond, double* a );
int init_cosineNonLin( unsigned int N, double* obs, double* cond,
	double* a, unsigned char* a_flag, FILE* logfile );

double fcosineNonLin( unsigned int i, double* cond, double* a );
double fcosineNonLin_deriv( double(*funct)(unsigned int, double*, double*),
	unsigned int i, unsigned int j, unsigned int M, double* cond, double* a );

double fcosine( unsigned int i, double *cond, double *a);
double fcosine_deriv( double (*funct)(unsigned int,double*,double*), 
					unsigned int i, unsigned int j, unsigned int M, double *cond, double *a);
int init_cosine( unsigned int N, double *obs, double *cond,
							double *a, unsigned char *a_flag, FILE *logfile);

double fcosine_trend( unsigned int i, double *cond, double *a);
double fcosine_trend_deriv( double (*funct)(unsigned int,double*,double*), 
					unsigned int i, unsigned int j, unsigned int M, double *cond, double *a);
int init_cosine_trend( unsigned int N, double *obs, double *cond,
					 double *a, unsigned char *a_flag, FILE *logfile);


double ftrigonometric2( unsigned int i, double *cond, double *a);
int init_cosine2( unsigned int N, double *obs, double *cond,
					 double *a, unsigned char *a_flag, FILE *logfile);

double flogarithmic( unsigned int i, double *cond, double *a);
double flogarithmic_deriv( double (*funct)(unsigned int,double*,double*), 
						unsigned int i, unsigned int j, unsigned int M, double *cond, double *a);
//double flogarithmic_deriv2( double (*funct)(unsigned int,double*,double*), 
//	unsigned int i, unsigned int j, unsigned int k, unsigned int M, double *cond, double *a);
int init_logarithmic( unsigned int N, double *obs, double *cond,
					 double *a, unsigned char *a_flag, FILE *logfile);

double fexponential( unsigned int i, double *cond, double *a);
double fexponential_deriv( double (*funct)(unsigned int,double*,double*), 
							unsigned int i, unsigned int j, unsigned int M, double *cond, double *a);
int init_exponential( unsigned int N, double *obs, double *cond,
						double *a, unsigned char *a_flag, FILE *out);
int init_exponentialLinearised( unsigned int N, double* obs, double* cond,
	double* a, unsigned char* a_flag, FILE* logfile );
int init_exponentialdummy( unsigned int N, double* obs, double* cond,
	double* a, unsigned char* a_flag, FILE* logfile );

double fexpon2( unsigned int i, double *cond, double *a);
int init_expon2( unsigned int N, double *obs, double *cond,
						double *a, unsigned char *a_flag, FILE *out);
double fgen_laplace( unsigned int i, double *cond, double *a);
double fgen_laplace_deriv( double (*funct)(unsigned int,double*,double*), 
						unsigned int i, unsigned int j, unsigned int M, double *cond, double *a);
int init_gen_laplace( unsigned int N, double *obs, double *cond,
						double *a, unsigned char *a_flag, FILE *out);
double fexpon2_deriv( double (*funct)(unsigned int,double*,double*), 
							unsigned int i, unsigned int j, unsigned int M, double *cond, double *a);

double fgauss2( unsigned int i, double *cond, double *a);
double fgauss1( unsigned int i, double *cond, double *a);

int init_gauss2( unsigned int N, double *obs, double *cond,
						double *a, unsigned char *a_flag, FILE *out);
int init_gauss1( unsigned int N, double *obs, double *cond,double *a,
					unsigned char *a_flag, int peak_flag, FILE *out);

double fgauss_deriv( double (*funct)(unsigned int,double*,double*), 
				unsigned int i, unsigned int j, unsigned int M, double *cond, double *a);
double fgauss_deriv2( double (*funct)(unsigned int,double*,double*), 
	unsigned int i, unsigned int j, unsigned int k, unsigned int M, 
	double *cond, double *a);
int init_gauss( unsigned int N, double *obs, double *cond,double *a, 
					unsigned char *a_flag, FILE *out);
double fhyperbel( unsigned int i, double* cond, double* a );
int init_hyperbel( unsigned int N, double* obs, double* cond, double* a,
	unsigned char* a_flag, FILE* out );

double frotation( unsigned int i, double *cond, double *a);
double
frotation_deriv( double (*funct)(unsigned int,double*,double*), 
				unsigned int i, unsigned int j, unsigned int M, double *cond, double *a);
int init_rotation( unsigned int N, double *obs, double *cond,
							 double *a, unsigned char *a_flag, FILE *logfile);

double fcircleTLS( unsigned int i, double *cond, double *a);
double fcircleTLS_deriv( double (*funct)(unsigned int,double*,double*), 
				unsigned int i, unsigned int j, unsigned int M, double *cond, double *a);
double fcircle( unsigned int i, double *cond, double *a);
double fcircle_deriv( double (*funct)(unsigned int,double*,double*), 
							unsigned int i, unsigned int j, unsigned int M, double *cond, double *a);
int init_circle( unsigned int N, double *obs, double *cond,
								 double *a, unsigned char *a_flag, FILE *logfile);

double fcirclelin_deriv( double (*funct)(unsigned int,double*,double*), 
				unsigned int i, unsigned int j, unsigned int M, double *cond, double *a);

double fNN_3_3( unsigned int i, double *cond, double *a);
int init_NN3x3x1( unsigned int N, double *obs, double *cond,
						 double *a, unsigned char *a_flag, FILE *logfile);
int init_NN( unsigned int N, double *obs, double *cond,
						 double *a, unsigned char *a_flag, FILE *logfile);
double fNN_3_2( unsigned int i, double *cond, double *a);
double fNN_2_2( unsigned int i, double *cond, double *a);
double fNN_1_2( unsigned int i, double *cond, double *a);
double fNN_1_3( unsigned int i, double *cond, double *a);
int init_NN1x3x1( unsigned int N, double *obs, double *cond,
							double *a, unsigned char *a_flag, FILE *logfile);
double fexpBaseX( unsigned int i, double *cond, double *a );
double fexpBaseX_deriv( double( *funct )(unsigned int, double*, double*),
	unsigned int i, unsigned int j, unsigned int M, double *cond, double *a );
int init_fexpBaseX( unsigned int N, double *obs, double *cond,
	double *a, unsigned char *a_flag, FILE *logfile );


double fExample( unsigned int i, double* cond, double* a );
double
fExample_deriv( double(*funct)(unsigned int, double*, double*),
	unsigned int i, unsigned int j, unsigned int M, double* cond, double* a );
double fExample_deriv2( double(*funct)(unsigned int, double*, double*),
	unsigned int i, unsigned int j, unsigned int k, unsigned int M,
	double* cond, double* a );
int init_Example( unsigned int N, double* obs, double* cond,
	double* a, unsigned char* a_flag, FILE* out );

#ifdef USED
double
fclassif2D( unsigned int i, double *cond, double *a );
#endif
double
fclassif2D_deriv( double( *funct ) (unsigned int, double*, double*),
	unsigned int i, unsigned int j, unsigned int M, double *cond, double *a );

/* common numerical derivation function for all 
 * nonlinear problems
 */
double f_deriv( double (*funct)(unsigned int,double*,double*), 
	unsigned int i, unsigned int j, unsigned int M, 
	double *cond, double *a);
double f_deriv2( double(*funct)(unsigned int, double*, double*),
	unsigned int i, unsigned int j, unsigned int k, unsigned int M,
	double* cond, double* a );


/* numerical derivation function for variable model functions
 * tinyepr
 */
double f_derivTinyExpr(te_expr* teExpr, unsigned int j, double* a);
double f_deriv2TinyExpr( te_expr* teExpr,
	unsigned int j, unsigned int k, double* a );
int init_tinyExp( unsigned int N, double* obs, double* cond,
	double* a, unsigned char* a_flag, FILE* logfile);

int lombScargle( unsigned int N, double* obs, double* cond,
	double obsMean, int fVecLen, double* fVec, double* p, double* phi );

#endif
