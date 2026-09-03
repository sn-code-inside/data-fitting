/*****************************************************************
 *
 * File........:	macros.h
 * Function....:	general purpose macros
 * Author......:	Tilo Strutz
 * last changes:	25.01.2010
 *
 * LICENCE DETAILS: see software manual
 *	free academic use
 *  cite source as 
 * "Strutz, T.: Data Fitting and Uncertainty. 3rd edition,
 *   Springer Nature, 2026"
 *
 *****************************************************************/

#ifndef MACROS_H
#define MACROS_H

#ifndef max
  #define max(a,b) (((a) > (b)) ? (a) : (b))
#endif
#ifndef min
  #define min(a,b) (((a) < (b)) ? (a) : (b))
#endif
#ifndef SIGN
  #define SIGN(x) (((x) > 0) ? 1 : (((x) < 0) ? -1: 0))
#endif


#endif
