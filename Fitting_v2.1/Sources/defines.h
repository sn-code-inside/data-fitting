/*************************************************************
 * File........:	defines.h
 * Function....:	definition of globally used values
 * Author......:	Tilo Strutz
 * last changes:	03.07.2009, 05.11.2025, 11.05.2026
 * 01.06.2026 COSINE_TREND removed
 * 19.07.2026 new NIST_BENNETT5_LOG
 * 12.08.2026 new EXAMPLE
 *
 * LICENCE DETAILS: see software manual
 *	free academic use
 *  cite source as 
 * "Strutz, T.: Data Fitting and Uncertainty. 3rd edition,
 *   Springer Vieweg,2026"
 *
 ***********************************************************/

#ifndef DEFINES_H
#define DEFINES_H

#include<float.h>

#define M_PI        3.14159265358979323846       /* pi */
#define CONVERGENCE_THRESH 1.0e-12		/* for convergence criterion */
#define TOL 1.0e-12		/* for  */
#define TOL_S2 1.0e-22	//1.0e-14		/* for test of maximum singular value */
#define TOL_S 1.0e-12		/* for test of singular values */
 /* this small value is required, for example,
     for the Bennett5 data set */

/* for est_weights.c and outlier_detetction.c	*/
#define MAX_LINES_W 500

/* for fitting.c	*/

/* maximal number of conditions per model function */
/* see functions.h	*/
#define MAX_CONDITIONS (M_MAX-1)
/* output of maximal 100 lines of resulting values as feedback */
#define MAX_LINES 100
/* input file may contain maximal 512 characters per line */
#define MAXLINELENGTH 512


/* for functions.c	*/

/* maximal 15 parameters per model function (see f_deriv()) */
#define M_MAX	20

/* for lsLinear.c	*/
extern long ITERAT_MAX;


/* model functions */
#define CONSTANT				0		/* constant value	*/
#define LINEAR					1		/* f(x|a) = a1 + SUM_j a_j * x_(j-1)	*/
#define LINEAR_2				2		/* f(x|a) =      SUM_j a_j * x_(j-1)	*/
#define POLYNOM_2NDORD	3	  /* 2nd order polynomial f(x|a) = a1 + a2 * x + a3 * x^2 */
#define POLYNOM_3RDORD	4	  /* 3rd order polynomial f(x|a) = a1 + a2 * x + a3 * x^2 + a4 * x^3 */
#define POLYNOMIAL			5	  /* multi-order polynomial f(x|a) = a1 + a2 * x + a3 * x^2 + ...*/
#define POLYNOMIAL_REG	6	  /* regularised   f(x|a) = a1 + a2 * x + a3 * x^2 + ...*/
#define COSINE_LIN			7		/* cosine linear: f(x|a) = a1 + a2 * cos( x) + a3 * sin( x) */
#define COSINE_NONLIN		8		/* cosine nonlinear: f(x|a) = a1 + a2 * cos( x - a3) */
#define COSINE					9		/* cosine nonlinear	f(x|a) = a1 + a2 * cos( a3*x + a4) */
#define COSINE2	        10	/* trigonometric polynom, 2st order	*/
#define COSINE_TREND		11	/* cosine with linear trend	*/
#define EXPONENTIAL			12	/* exponential f(x|a) = a1 + a2 * exp( a3 * x)	*/
#define EXPONENTIAL2		13	/* exponential 2 f(x|a) = a2 * exp( a3 * x) 	*/
#define EXPONENTIAL2_LIN	14	/* exponential 2, linearised: ln(f(x|a)) = ln(a2) + a3 * x*/
#define EXP_BASE_X			15	/* f(x|a) = (a1*x)^a2 + a3	*/
#define LOGARITHM				16	/* f(x|a) = a1 + a2 * log( a3 + x) 	*/
#define HYPERBEL			  17	/* Hyperbel  a1/(x+a2)	*/
#define HYPERBEL_LIN	  18	/* Hyperbel  a1/(x+a2), linearised */
#define GAUSSIAN_1			19	/* Gaussian	f(x|a) = a1 * exp( a2 * (x-a3)^2) */
#define GAUSSIAN_2			20	/* superposition of two Gaussians: 
                             * f(x|a) = a1 * exp( a2 * (x-a3)^2) + 
			                       *          a4 * exp( a5 * (x-a6)^2) */
#define GEN_LAPLACE			21	/* generalised Laplacian distribution f(x|a) = a1 * exp( -|x|^a2 * a3) */

#define QUAD_SURFACE		30	/* parabolic 2D surface: f(x|a) = a1 + a2*x1 + a3*x1^2 + a4*x2 +a5*x2^2  */
#define COORD_TRANSF		31	/* rotation + translation	*/
                             /*
                              * f1(x|a) = a1 + cos(a3) * x1 - sin(a3) * x2
                              * f2(x|a) = a2 + sin(a3) * x1 + cos(a3) * x2
                              */
#define CIRCLE					32	/* circle */
                              /*
                               * f(x|a) = 0 = (x1 - a1)^2 + (x2 - a2)^2 - a3*a3
                               */
#define CIRCLE_LIN			33	/* circle, linearised */
                               /*
                                * f(x|a) = 0 = (x1 - a1)^2 + (x2 - a2)^2 - a3*a3
                                * f(x|b) = x1^2 + x2^2 = b1*x1 + b2*x2 - b3
                                * b1 = 2*a1, b2 = 2*a2, b3 = a1^2 + a2^2 - a3^2
                                */
#define CIRCLE_TLS			34	/* circle, total least squares */
                                /*
                                 * f(x|a) = 0 = (sqrt[(x1 - a1)^2 + (x2 - a2)^2] - a3)^2
                                 */
#define CLASSIF_2D			35	/* f(x|a) = a1 + a2*x1 + a3*x2 + a4*x1^2 + a5*x2^2 + a6*x1*x2 + a7*x1^3 + a8*x2^3	*/

#define NN_3x3x1				40	/* NeuralNet 3x3x1 */
#define NN_3x2x1				41	/* NeuralNet 3x2x1 */
#define NN_1x2x1				42	/* NeuralNet 1x2x1 */
#define NN_2x2x1				43	/* NeuralNet 2x2x1 */
#define NN_1x3x1				44	/* NeuralNet 1x3x1 */

/* Filip = POLYNOMIAL         f(x|a) = a1 + a2 * x + a3 * x^2 + ... a11 * x^10 */
/* Wampler 1..6 = POLYNOMIAL  f(x|a) = a1 + a2 * x + a3 * x^2 + ... a6 * x^5 */
/* Longley = LINEAR           f(x|a) = a1 + a2 * x_1 + a3 * x_2 + ... a7 * x_6 */
#define NIST_THURBER		50	/* NIST data set */
#define NIST_MGH09			51	/* NIST data set */
#define NIST_RAT42			52	/* NIST data set */
#define NIST_RAT43			53	/* NIST data set */
#define NIST_ECKERLE4		54	/* NIST data set */
#define NIST_MGH10			55	/* NIST data set */
#define NIST_BENNETT5		56	/* NIST data set */
#define NIST_BOXBOD			57	/* NIST data set	*/
#define NIST_BENNETT5_LOG	58/* NIST data set */
#define UNSPECIFIED   99
#define EXAMPLE   100


#ifdef DBL_DECIMAL_DIG
  #define OP_DBL_Digs (DBL_DECIMAL_DIG)
#else  
  #ifdef DECIMAL_DIG
    #define OP_DBL_Digs (DECIMAL_DIG)
  #else  
    #define OP_DBL_Digs (DBL_DIG + 3)
  #endif
#endif

#endif
