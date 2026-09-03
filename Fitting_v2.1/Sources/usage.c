/*****************************************************************
 *
 * File....: usage.c
 * Function: parameters for Fitting
 * Author..: Tilo Strutz
 * Date... : 07.05.2008, 01.10.2009, 6.11.2009, 08.01.2010
 *           18.02.2010, 10.03.2010, 05.11.2025, 11.05.2026
 * changes:
 * 28.01.2014 new option cw
 * 25.02.2021 new: f(x|a) = (a1*x)^a2 + a3 and f(x|a) = a2 * log(a1*x)
 * 17.09.2022 new: 50 ... f(x|a) = a1 + a2*x1 + a3*x2 + a4*x1^2 + a5*x2^2 + a6*x1*x2 + a7*x1^3 + a8*x2^3
 * 14.11.2025 new options '-n+', '-y', '-z1', '-z2'
 * 19.06.2026 new option '-D'
 * 26.06.2026 option '-cw' changed to '-cu'
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

/* allowed options:  must start with ':'  !!  */
char *OPTIONSTRING = 
{":a:a1.a2.a3.a4.a5.a6.a7.a8.a9.a10.a11.b:c;cc:co:cu:D;e;f;H;i:I:g;G;m:M:n;n+;o:p;s;t:w:x:y.z0:z1:z2:"};
/* .  ==> argument may start with '-' (minus) */
/* :  ==> argument mmust follow */
/* ;  ==> argument not allowed */
char *NEEDEDOPTIONS = { ":i:o:m:" };  /* required options */

char *title = { "Fitting WLS version 2.1 (08/2026)" };

/*---------------------------------------------------------------
 *  usage()
 *---------------------------------------------------------------*/
void
usage( char *name)
{
    fprintf( stderr, "\n\n %s\n", title);
    fprintf( stderr, "\n\
Usage: %s [options]\n\n\
 Legal Options: \n\
    -i %%s ... input data file (compulsory)\n\
    -o %%s ... output file (compulsory)\n\
    -m %%d ... model function (compulsory)\n\
       0 ... constant f(x|a) = a1\n\
       1 ... f(x|a) = a1 + SUM_{j=2}^M a_j * x_(j-1)\n\
       2 ... f(x|a) =      SUM_{j=2}^M a_j * x_(j-1)\n\
       3 ... f(x|a) = a1 + a2 * x + a3 * x^2\n\
       4 ... f(x|a) = a1 + a2 * x + a3 * x^2 + a4 * x^3\n\
       5 ... f(x|a) = SUM_{j=1}^M aj * x^(j-1) (linear)\n\
       6 ... f(x|a) = SUM_{j=1}^M aj * x^(j-1) (nonlinear)\n\
       7 ... f(x|a) = a1 + a2 * cos( x - a3) (= a1 + a2 * cos( x) + a3 * sin( x)) [in degrees]\n\
       8 ... f(x|a) = a1 + a2 * cos( a3*x + a4) [ in radians]\n\
      10 ... f(x|a) = a1 + a2*cos(a3*x-a4) + a5*cos(2*a3*x-a6) [in radians]\n\
      11 ... f(x|a) = a1 + a2 * x + a3 * cos( x - a4) [in radians]\n\
      12 ... f(x|a) = a1 + a2 * exp( a3 * x)\n\
      13 ... f(x|a) = a1 * exp( a2 * x)\n\
      14 ... ln(f(x|a)) = ln(a1) + a2 * x\n\
      15 ... f(x|a) = (a1 * x)^a2 + a3\n\
      16 ... f(x|a) = a1 + a2 * log( a3 + x)\n\
      17 ... f(x|a) = a1 / (x+a2) \n\
      18 ... 1 / f(x|a) = (x+a2) / a1 = b1 + b2 * x   \n\
      19 ... f(x|a) = a1 * exp( a2 * (x-a3)^2) \n\
      20 ... f(x|a) = a1 * exp( a2 * (x-a3)^2) + \n\
                      a4 * exp( a5 * (x-a6)^2)\n\
      21 ... f(x|a) = a1 * exp( -|x|^a2 * a3)\n\
      30 ... f(x|a) = a1 + a2*x1 + a3*x1^2 + a4*x2 +a5*x2^2\n\
      31 ... f1(x|a) = a1 + cos(a3) * x1 - sin(a3) * x2 [in radians]\n\
             f2(x|a) = a2 + sin(a3) * x1 + cos(a3) * x2 [in radians]\n\
      32 ... f(x|a) = 0 = (x1 - a1)^2 + (x2 - a2)^2 - a3*a3 (circle)\n\
      33 ... f(x|a) = x1^2 + x2^2 = a1*x1 + a2*x2 - a3\n\
                               (circle, linear)\n\
      34 ... f(x|a) = 0 = (sqrt[(x1 - a1)^2 + (x2 - a2)^2] - a3)^2\n\
                                          (circle, TLS)\n\
      35 ... f(x|a) = a1 + a2*x1 + a3*x2 + a4*x1^2 + a5*x2^2 + a6*x1*x2 + a7*x1^3 + a8*x2^3\n\
      40 ... neural network 3x3x1, feed forward\n\
      41 ... neural network 3x2x1\n\
      42 ... neural network 1x2x1\n\
      43 ... neural network 2x2x1\n\
      44 ... neural network 1x3x1\n\
      50 ... f(x|a) =(a1 + a2*x + a3*x*x + a4*x*x*x) /\n\
                      (1 + a5*x + a6*x*x + a7*x*x*x) NIST_THURBER\n\
      51 ... f(x|a) = a1 * (x**2 + a2*x) /\n\
                      (x*x + a3*x + a4)              NIST_MGH09\n\
      52 ... f(x|a) = a1 / (1 + exp(a2 - a3*x))         NIST_Rat42\n\
      53 ... f(x|a) = a1 / [1 + exp(a2 - a3*x)]^(1/a4)	NIST_Rat43\n\
      54 ... f(x|a) = (a1 / a2) * exp(-0.5*((x -a3)/ a2)^2) NIST_Eckerle4\n\
      55 ... f(x|a) = a1 * exp( a2 / (x+a3))	NIST_MGH10\n\
      56 ... f(x|a) = a1 * (x+a2)^(-1/a3)     NIST_Bennett5\n\
      57 ... f(x|a) = a1 *(1 - exp( -a2 * x)  NIST_BoxBOD\n\
      58 ... log(-f(x|a)) = log(-a1) + (-1/a3) * log(x+a2)  NIST_Bennett5 transformed\n\
    -a  %%d ... inversion algorithm (default: 1)\n\
           0 - cofactor method\n\
           1 - singular value decomposition\n\
           2 - LU decomposition\n\
    -a[j]  %%f ... provides initial value for  a_j (j=1,2,..,11)\n\
    -b  %%d ... observations per bin, for '-w 2' (default: 50)\n\
    -c     ... enable scaling of conditions \n\
    -cc %%s ... comma-separated list of column(s) containing \n\
                conditions x (default: 1,2,...)\n\
    -co %%d ... column containing observations y (default: 2)\n\
    -cu %%d ... column containing uncertainties\n\
    -e     ... enable uphill walk and damping of step size, only in combination with '-G'\n\
    -D     ... use diagonal elements of normal matrix instead of identity matrix for Levenberg-Marquardt	 \n\
    -f     ... forget weights after outlier removal\n\
    -H     ... enable true Hessian matrix\n\
    -I %%d   ... maximum number of iterations (default: 2000)\n\
    -M %%d   ... number of parameters (for '-m 1|2|5|6' only)\n\
    -n     ... force usage of numerical derivation\n\
    -g     ... gradient manipulation at saddle points	 \n\
    -G     ... use Gauss-Newton instead of Levenberg-Marquardt	 \n\
    -p     ... in conjunction with '-y': keep parameters positive\n\
    -s     ... disable special SVD function for solving linear model\n\
    -t %%f  ... target value for chisq (maximum error)\n\
                default: iteration until convergence\n\
    -w %%d   ... weighting  (default: 0)\n\
           0 ... no weighting \n\
           1 ... based on deviates \n\
           2 ... binning \n\
           3 ... use uncertainties from file (if -cu is provided) \n\
    -x %%d   ... outlier removal (default: 0)\n\
           0 ... no outlier removal\n\
           1 ... z-score + Chauvenet's criterion\n\
           2 ... cluster criterion (ClubOD)\n\
           3 ... M-score + Chauvenet's criterion \n\
           4 ... RANSAC \n\
     -y %%s  ... model function (e.g. \"a0+a1*x+ a2 * sin(a3*x+a4)\" )\n\
     -z0 %%d  ... format string: number of digits for intermediate output (default: 12)\n\
     -z1 %%d  ... format string: total number of digits for final output (default: 12)\n\
     -z2 %%d  ... format string: number of decimal digits  (default: 6)\n\
    \n\
 ", name);
}
