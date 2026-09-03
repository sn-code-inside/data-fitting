/*****************************************************************
 *
 * File........:	matrixUtils.h
 * Function....:	special functions (prototyping)
 * Author......:	Tilo Strutz
 * last changes: 20.10.2007, 29.3.2011, 16.02.2026
 *	07.08.2026
 * LICENCE DETAILS: see software manual
 *	free academic use
 *  cite source as 
 * "Strutz, T.: Data Fitting and Uncertainty. 3rd edition,
 *   Springer Nature, 2026"
 *
 *****************************************************************/

#ifndef MATRIX_UTILS_H
#define MATRIX_UTILS_H

#include<float.h>

int *ivector( unsigned long N);
long *lvector( unsigned long N);
unsigned int *uivector( unsigned long N);
float *fvector( unsigned long N);
double *vector( unsigned long N);
double **matrix( unsigned long N, unsigned long M);
float **fmatrix( unsigned long N, unsigned long M);

void free_ivector( int *v[]);
void free_lvector( long* v[] );
void free_uivector( unsigned int *v[]);
void free_vector( double *v[]);
void free_matrix( double **m[]);
void free_fmatrix( float **m[]);


double determinant_2x2( double **a);
double determinant_3x3( double **a);
double determinant_3x3_from_cofactor( double** a, double** b );
double inverse_4x4( double **a, double **b);
double inverse_5x5( double **a, double **b);
void coFactor_2x2( double **a, double **b);
void coFactor_3x3( double **a, double **b);

void multmatsq( unsigned int M, double **a, double **b, double **c);
void multmatsqT( unsigned int N, double **a, double **b, double **c);

/*---------------------------------------------------------------
 *	diffOfProducts()
 * fma-compensated difference of products a*b - c*d, see Kahan (2004) /
 * Muller et al., Handbook of Floating-Point Arithmetic
 *--------------------------------------------------------------*/
static inline double
diffOfProducts( double a, double b, double c, double d )
{
	double cd = c * d;
	double err = fma( c, d, -cd );   /* rounding error of c*d */
	double dop = fma( a, b, -cd );   /* a*b - cd, ecakt via (Sterbenz/FMA) */
	return dop - err;
}
/*---------------------------------------------------------------
 *	diffOfProducts()
 * fma-compensated summation of products a*b - c*d
 *--------------------------------------------------------------*/
static inline double
sumOfProducts( double a, double b, double c, double d )
{
	double cd = c * d;
	double err = fma( c, d, -cd );   // rounding error of  c*d
	double sop = fma( a, b, cd );    // a*b + w, gerundet
	return sop + err;
}

#endif
