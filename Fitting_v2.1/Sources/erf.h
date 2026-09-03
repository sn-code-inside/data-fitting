/***************************
 * File........:	erf.h
 * Function....:	prototypes for erf.c
 * Author......:	Tilo Strutz
 * last changes:	05.02.2008
 *
 * LICENCE DETAILS: see software manual
 *	free academic use
 *  cite source as 
 * "Strutz, T.: Data Fitting and Uncertainty. 3rd edition,
 *   Springer Nature, 2026"
 *
 ***************************/
double erf(double x);
double erfc(double x);
int erfinv( double y, double *res );
