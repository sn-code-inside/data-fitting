/************************************************************
 * File........:	erf.c
 * Function....:	error functions
 * Author......:	Tilo Strutz
 * last changes:	05.02.2008, 28.07.2026
 *
 * LICENCE DETAILS: see software manual
 *	free academic use
 *  cite source as 
 * "Strutz, T.: Data Fitting and Uncertainty. 3rd edition,
 *   Springer Nature, 2026"
 *
 ***********************************************************/
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "erf.h"
#include "errmsg.h"

static const double rel_error= 1E-12;
/* calculate 12 significant figures
 * you can adjust rel_error to trade off between accuracy and 
 * speed, but don't ask for > 15 figures 
 *(assuming usual 52 bit mantissa in a double)
 */

/*---------------------------------------------------------------
 *	erf()
 *
 *  erf(x) = 2/sqrt(pi)*integral(exp(-t^2),t,0,x)
 *         = 2/sqrt(pi)*[x - x^3/3 + x^5/5*2! - x^7/7*3! + ...]
 *         = 1-erfc(x)
 *
 *--------------------------------------------------------------*/
double erf(double x)
{
	/* 2/sqrt(pi)*/
  static const double two_sqrtpi = 1.128379167095512574; 
  double sum, term, xsqr;
  int j= 1;

    sum = x;
		term = x;
		xsqr = x * x;
	  
		if (fabs(x) > 2.2)
		{
			/*use continued fraction when fabs(x) > 2.2 */
      return 1.0 - erfc(x); 
    }
    do
		{
			term *= xsqr / j;
			sum -= term /(2*j+1);
			j++;
			term *= xsqr / j;
			sum += term /(2*j+1);
			j++;
    } while (fabs(term) / sum > rel_error);
    return two_sqrtpi * sum;
}

/*---------------------------------------------------------------
 *	erfc()
 *
 *  erfc(x) = 2/sqrt(pi)*integral(exp(-t^2),t,x,inf)
 *         = exp(-x^2)/sqrt(pi) * [1/x + (1/2)/x + (2/2)/x + 
 *																(3/2)/x+ (4/2)/x+ ...]
 *         = 1-erf(x)
 *  expression inside [] is a continued fraction 
 *	so '+' means add to denominator only
 *--------------------------------------------------------------*/
double erfc(double x)
{
	/* 1/sqrt(pi)*/
	static const double one_sqrtpi = 0.564189583547756287;
	double a = 1, b;        /* last two convergent numerators	*/
	double c, d;  /* last two convergent denominators	*/
	double q1, q2;     /* last two convergents (a/c and b/d)	*/
	double n = 1.0, t;

		b = c = x;
		d = x * x + 0.5;
		q2 = b / d;

		if (fabs(x) < 2.2)
		{
			return 1.0 - erf(x);   /* use series when fabs(x) < 2.2	*/
		}
		if (x > 0)
		{               /* continued fraction only valid for x>0 */
			return 2.0 - erfc(-x);
		}
		do
		{
			t = a*n + b*x;
			a = b;
			b = t;
			t = c*n + d*x;
			c = d;
			d = t;
			n += 0.5;
			q1 = q2;
			q2 = b / d;
		} while (fabs(q1-q2) / q2 > rel_error);
		return one_sqrtpi * exp(-x*x) * q2;
}

/*---------------------------------------------------------------
 *	erfinv()
 * Approximation of the inverse error function erfinv(y).
 * Input:
 *   y    : Function argument. Valid range: -1 < y < 1
 *   res  : Pointer to result variable
 * Output:
 *   *res : Approximation of erfinv(y)
 * Return value:
 *   0 : Success
 *   3 : Argument outside valid range
 *--------------------------------------------------------------*/
int erfinv( double y, double* res )
{
  /* Coefficients for rational approximations */
  static double a[] = { 0,  0.886226899, -1.645349621,
                            0.914624893, -0.140543331 };
  static double b[] = { 0, -2.118377725,  1.442710462,
                           -0.329097515,  0.012229801 };
  static double c[] = { 0, -1.970840454, -1.624906493,
                            3.429567803,  1.641345311 };
  static double d[] = { 0,  3.543889200,  1.637067800 };

  int err = 0;      /* Error status */
  double x, z;

  /* Check lower domain boundary: y must be greater than -1 */
  if (y > -1.0)
  {
    if (y >= -0.7)
    {
      if (y <= 0.7)
      {
        /*
         * Central region:
         * Polynomial/rational approximation around y = 0.
         */
        z = y * y;
        x = y *
          (((a[4] * z + a[3]) * z + a[2]) * z + a[1]) /
          ((((b[4] * z + b[3]) * z + b[2]) * z + b[1]) * z + 1);
      }
      else if (y < 1.0)
      {
        /* Upper tail:
         * Approximation for values close to +1.
         */
        z = sqrt( -log( (1 - y) / 2 ) );
        x = (((c[4] * z + c[3]) * z + c[2]) * z + c[1]) /
          ((d[2] * z + d[1]) * z + 1);
      }
      else
      {
        /* Error case: y >= 1.0
         * erfinv(1) tends to +infinity and therefore
         * is not handled by this implementation.
         */
        err = ERR_OUT_OF_RANGE;
        x = 0;
      }
    }
    else
    {
      /* Lower tail:
       * Approximation for values close to -1.
       */
      z = sqrt( -log( (1 + y) / 2 ) );
      x = -(((c[4] * z + c[3]) * z + c[2]) * z + c[1]) /
        ((d[2] * z + d[1]) * z + 1);
    }
  }
  else
  {
    /*
     * Error case:  y <= -1.0
     *
     * erfinv(-1) tends to -infinity and therefore
     * is not handled by this implementation.
     */
    err = ERR_OUT_OF_RANGE;
    x = 0;
  }

  *res = x;
  return err;
}