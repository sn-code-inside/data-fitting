/***********************************************************
 *
 * File........:	singvaldec.c
 * Function....:	singular value decomposition
 * Author......:	Tilo Strutz
 * last changes:	20.10.2007, 22.07.2026
 *
 * LICENCE DETAILS: see software manual
 *	free academic use
 *  cite source as
 * "Strutz, T.: Data Fitting and Uncertainty. 3rd edition,
 *   Springer Nature, 2026"
 *
 ***********************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>/* for copysign() */
#include <float.h> /* for DBL_MAX */
#include <assert.h>
#include "matrixUtils.h"
#include "macros.h"
#include "errmsg.h"

//#define COPYSIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))
#define MAX_ITERATIONS 30

#ifdef OBSOLETE_HYPOT_USED
/*---------------------------------------------------------
 * euclid_dist()
 *
 * Compute sqrt(a² + b²) in a numerically stable
 * manner without unnecessary overflow or
 * underflow.
 *
 * Functionally equivalent to hypot(a,b).
*---------------------------------------------------------*/
double
euclid_dist( double a, double b )
{
	double abs_a, abs_b, val, dval;

	abs_a = fabs( a );
	abs_b = fabs( b );
	if (abs_a > abs_b)
	{
		dval = abs_b / abs_a;
		//val = abs_a * sqrt( 1.0 + dval * dval );
		val = abs_a * sqrt( fma(dval, dval, 1.0 ));
		return val;
	}
	else
	{
		if (abs_b == 0.0)		return 0.0;
		else
		{
			dval = abs_a / abs_b;
			//val = abs_b * sqrt( 1.0 + dval * dval );
			val = abs_b * sqrt( fma( dval, dval, 1.0 ));
			return val;
		}
	}
}
#endif

/*---------------------------------------------------------------
 * singvaldec()
 *
 * Singular Value Decomposition (SVD)
 *
 * Decomposes a matrix A into
 *
 *      A = U * W * V^T
 *
 * where
 *
 *      U : orthogonal left singular vectors
 *      W : diagonal matrix of singular values
 *      V : orthogonal right singular vectors
 *
 * On return:
 *
 *      a -> contains U
 *      w -> contains the singular values
 *      v -> contains V
 *
 * Algorithm:
 *
 * 1. Householder transformations reduce A to
 *    bidiagonal form.
 *
 * 2. The accumulated Householder reflections
 *    are used to construct the orthogonal
 *    matrices U and V.
 *
 * 3. An iterative QR algorithm diagonalizes
 *    the bidiagonal matrix and computes the
 *    final singular values.
 *
 * Numerically sensitive operations use fused
 * multiply-add (FMA) evaluations whenever
 * possible to reduce rounding errors.
 * 
 * translation from http://www.pdas.com/programs/fmm.f90
 *--------------------------------------------------------------*/
int
singvaldec( double** a, /* matrix to be decomposed */
	unsigned int N, /* number of lines */
	unsigned int M, /* number of columns	*/
	double w[], /* matrix for sigular values */
	double** v )
{
	char* rtn = "singvaldec";
	int err = 0;
	int l = 0, flag;
	int nm=0;
	double anorm, c, f, g, h, s, scale, x, y, z, * rv1;

	if (a == NULL || w == NULL || v == NULL)
	{
		err = ERR_ANY;
		return err;
	}
	if (M==0 || N < M)
	{
		err = ERR_ANY;
		return err;
	}

	rv1 = vector( M ); /* check inside vector() */

	/*
	 * Phase 1:
	 * Reduce the input matrix to bidiagonal form
	 * using Householder reflections.
	 *
	 * Alternating left and right Householder
	 * transformations eliminate matrix elements
	 * below the diagonal and to the right of the
	 * first superdiagonal.
	 *
	 * The Householder vectors are stored
	 * implicitly and later reused to construct
	 * the orthogonal matrices U and V.
	 */
	g = scale = anorm = 0.0;
	for (int i = 0; i < (signed)M; i++)
	{
		l = i + 1;
		assert( i >= 0 );
		rv1[i] = scale * g;
		g = s = scale = 0.0;
		if (i < (signed)N)
		{
			for (unsigned int k = i; k < N; k++)
			{
				scale += fabs( a[k][i] );
			}
			if (scale ==  INFINITY) scale = DBL_MAX;
			if (scale)
			{
				for (unsigned int k = i; k < N; k++)
				{
					a[k][i] /= scale;
					//s += a[k][i] * a[k][i];
					s = fma( a[k][i], a[k][i], s );
				}
				f = a[i][i];

				if (s < 0 || !isfinite( s ))
				{ /* checks for NAN and INF*/
					err = ERR_SVD_STRANGE_INPUT;
					goto endfunc;
				}
				//assert( s >= 0 );
				g = -copysign( sqrt( s ), f );
				// h = f * g - s;
				h = fma( f, g, - s);
				a[i][i] = f - g;
				for (unsigned int j = l; j < M; j++)
				{
					s = 0.0;
					for (unsigned int k = i; k < N; k++)
					{
						//s += a[k][i] * a[k][j];
						s = fma( a[k][i], a[k][j], s );
					}
					f = s / h;
					for (unsigned int k = i; k < N; k++)
					{
						//a[k][j] += f * a[k][i];
						a[k][j] = fma( f, a[k][i], a[k][j] );
					}
				}
				for (unsigned int k = i; k < N; k++)
				{
					a[k][i] *= scale;
					if (a[k][i] == INFINITY) a[k][i] = DBL_MAX;
					if (a[k][i] == -INFINITY) a[k][i] = -DBL_MAX;
				}
			}
		}
		w[i] = scale * g;
		g = s = scale = 0.0;
		if (i < (signed)N && i != (signed)M - 1)
		{
			for (unsigned int k = l; k < M; k++)
			{
				scale += fabs( a[i][k] );
			}
			if (scale == INFINITY) scale = DBL_MAX;
			if (scale)
			{
				for (unsigned int k = l; k < M; k++)
				{
					a[i][k] /= scale;
					//s += a[i][k] * a[i][k];
					s = fma( a[i][k], a[i][k], s );
				}
				f = a[i][l];
				if (s < 0 || !isfinite( s ))
				{ /* checks for NAN and INF*/
					err = ERR_SVD_STRANGE_INPUT;
					goto endfunc;
				}
				//assert( s >= 0. );
				g = -copysign( sqrt( s ), f );
				//h = f * g - s;
				h = fma( f, g, - s);
				a[i][l] = f - g;
				for (unsigned int k = l; k < M; k++)
				{
					rv1[k] = a[i][k] / h;
				}
				for (unsigned int j = l; j < N; j++)
				{
					s = 0.0;
					for (unsigned int k = l; k < M; k++)
					{
						//s += a[j][k] * a[i][k];
						s = fma( a[j][k], a[i][k], s );
					}
					for (unsigned int k = l; k < M; k++)
					{
						//a[j][k] += s * rv1[k];
						a[j][k] = fma( s, rv1[k], a[j][k] );
					}
				}
				for (unsigned int k = l; k < M; k++)
				{
					a[i][k] *= scale;
					if (a[i][k] == INFINITY) a[i][k] = DBL_MAX;
					if (a[i][k] == -INFINITY) a[i][k] = -DBL_MAX;
				}
			}
		}
		anorm = max( anorm, (fabs( w[i] ) + fabs( rv1[i] )) );
	}

	/*
	 * Phase 2:
	 * Construct the orthogonal matrix V from the
	 * accumulated right Householder reflections.
	 *
	 * The transformations are applied in reverse
	 * order to recover the complete set of right
	 * singular vectors.
	 */
	for ( int i = (signed)M - 1; i >= 0; i--)
	{
		if (i < (signed)M - 1)
		{
			if (g)
			{
				for (unsigned int j = l; j < M; j++)
				{
					if (a[i][l] == 0.0)
					{
						err = errmsg( ERR_IS_ZERO, rtn, "a[i][l]", 0 );
						goto endfunc;
					}
					v[j][i] = (a[i][j] / a[i][l]) / g;
				}
				for (unsigned int j = l; j < M; j++)
				{
					s = 0.0;
					for (unsigned int k = l; k < M; k++)
					{
						//s += a[i][k] * v[k][j];
						s = fma( a[i][k], v[k][j], s );
					}
					for (unsigned int k = l; k < M; k++)
					{
						//v[k][j] += s * v[k][i];
						v[k][j] = fma( s, v[k][i], v[k][j] );
					}
				}
			}
			for (unsigned int j = l; j < M; j++)
				v[i][j] = v[j][i] = 0.0;
		}
		v[i][i] = 1.0;
		assert( i >= 0 );
		g = rv1[i];
		l = i;
	}

	/*
	 * Phase 3:
	 * Construct the orthogonal matrix U from the
	 * accumulated left Householder reflections.
	 *
	 * The Householder information stored in A is
	 * expanded into the explicit matrix of left
	 * singular vectors.
	 *
	 * After this step, matrix A contains U.
	 */
	for (int i = (int)min( N, M ) - 1; i >= 0; i--)
	{
		l = i + 1;
		g = w[i];
		for (unsigned int j = l; j < M; j++)
			a[i][j] = 0.0;
		if (g && a[i][i])
		{
			g = 1.0 / g;
			for (unsigned int j = l; j < M; j++)
			{
				s = 0.0;
				for (unsigned int k = l; k < N; k++)
				{
					//s += a[k][i] * a[k][j];
					s = fma( a[k][i], a[k][j], s );
				}
				f = (s / a[i][i]) * g;
				for (unsigned int k = i; k < N; k++)
				{
					//a[k][j] += f * a[k][i];
					a[k][j] = fma( f, a[k][i], a[k][j] );
				}
			}
			for (unsigned int j = i; j < N; j++)
				a[j][i] *= g;
		}
		else
		{
			for (unsigned int j = i; j < N; j++)
				a[j][i] = 0.0;
		}
		++a[i][i];
	}

	/*
	 * Phase 4:
	 * Diagonalize the bidiagonal matrix using an
	 * iterative QR algorithm.
	 *
	 * Successive Givens rotations eliminate the
	 * remaining off-diagonal elements while
	 * preserving orthogonality.
	 *
	 * The diagonal elements converge to the
	 * singular values and the corresponding
	 * updates are propagated to U and V.
	 */
	for (int k = (signed)M - 1; k >= 0; k--)	/* loop over singular values */
	{
		/* Process singular values from the last column
		 * towards the first and force convergence of
		 * each singular value individually.
		 */
		assert( k >= 0 );
		unsigned int its; /* is checked outside the loop */
		for ( its = 0; its < MAX_ITERATIONS; its++)
		{
			/* QR iteration loop.
			 * Up to 30 iterations are allowed for each
			 * singular value before convergence is
			 * considered to have failed.
			 */
			flag = 1;
			for (l = k; l >= 0; l--)	/* test for splitting */
			{
				/* Determine whether the bidiagonal matrix can
				 * be split at position l by testing whether
				 * the current superdiagonal element has become
				 * negligible relative to the matrix norm.
				 */
				nm = l - 1;	/* note that rv1[0] is always zero */
				if ((double)(fabs( rv1[l] ) + anorm) == anorm)
				{
					flag = 0;
					break;
				}
				assert( nm >= 0 ); /* should never be violated */
				if ((double)(fabs( w[nm] ) + anorm) == anorm)
					break;
			}
			if (flag)
			{
				assert( l >= 0 );
				c = 0.0;	/* cancellation of rv1[l] if l greater than 1 */
				/* Apply a sequence of Givens rotations to
				 * eliminate rv1[l] and restore bidiagonal
				 * structure before continuing the QR step.
				 */
				s = 1.0;
				for ( int i = l; i < k; i++)
				{
					f = s * rv1[i];
					rv1[i] = c * rv1[i];
					if ((double)(fabs( f ) + anorm) == anorm)
						break;
					g = w[i];
					h = hypot( f, g );
					w[i] = h;
					h = 1.0 / h;
					c = g * h;
					s = -f * h;
					for (unsigned int j = 0; j < N; j++)
					{
						y = a[j][nm];
						z = a[j][i];
						//a[j][nm] = y * c + z * s;
						a[j][nm] = sumOfProducts( y, c, z, s );
						//a[j][i] = z * c - y * s;
						a[j][i] = diffOfProducts( z, c, y, s );
					}
				}
			}
			/* Check whether the current singular value
			 * has converged. Singular values are forced
			 * to be non-negative by adjusting the sign of
			 * the corresponding column of V.
			 */
			z = w[k];
			if (l == k)
			{
				if (z < 0.0)	/* singular value is made non-negative */
				{
					w[k] = -z;
					for (unsigned int j = 0; j < M; j++)
						v[j][k] = -v[j][k];
				}
				break;
			}
			/* Compute the implicit Wilkinson shift from
			 * the lower 2x2 submatrix.
			 *
			 * The shift significantly accelerates
			 * convergence of the QR iteration.
			 */
			x = w[l];
			nm = k - 1;
			y = w[nm];
			g = rv1[nm];
			assert( k >= 0 );
			h = rv1[k];
			//f =	((y - z) * (y + z) + (g - h) * (g +	h))
			//	  / (2.0 * h * y);
			double p1 = y - z, p2 = y + z, p3 = g - h, p4 = g + h;
			f = sumOfProducts( p1, p2, p3, p4 ) / (2.0 * h * y);
			if (!isfinite( f ))
			{
				err = errmsg( ERR_IS_INFINITE, rtn, "f", 0 );
				goto endfunc;
			}

			g = hypot( f, 1.0 );

			if (x == 0.0)
			{
				err = errmsg( ERR_IS_ZERO, rtn, "x", 0 );
				goto endfunc;
			}
			f =
				((x - z) * (x + z) +
					h * ((y / (f + copysign( g, f ))) - h)) / x;

			/* Perform one implicit QR step using Givens
			 * rotations.
			 *
			 * The rotations are applied to both U and V,
			 * preserving the factorization
			 *
			 *      A = U * W * V^T
			 *
			 * while driving the off-diagonal elements
			 * towards zero.
			 */
			c = s = 1.0;
			for ( int j = l; j <= nm; j++)
			{
				int i = j + 1;
				assert( i >= 0 );
				g = rv1[i];
				y = w[i];
				h = s * g;
				g = c * g;
				z = hypot( f, h );
				rv1[j] = z;
				c = f / z;
				s = h / z;
				//f = x * c + g * s;
				f = sumOfProducts( x, c, g, s );
				//g = g * c - x * s;
				g = diffOfProducts( g, c, x, s );
				h = y * s;
				y *= c;
				for (unsigned int jj = 0; jj < M; jj++)
				{
					x = v[jj][j];
					z = v[jj][i];
					//v[jj][j] = x * c + z * s;
					v[jj][j] = sumOfProducts( x, c, z, s );
					//v[jj][i] = z * c - x * s;
					v[jj][i] = diffOfProducts( z, c, x, s );
				}
				z = hypot( f, h );
				w[j] = z;
				if (z)	/* rotation can be arbitrary if z is zero */
				{
					z = 1.0 / z;
					c = f * z;
					s = h * z;
				}
				//f = c * g + s * y;
				f = sumOfProducts( c, g, s, y );
				//x = c * y - s * g;
				x = diffOfProducts( c, y, s, g );
				for (unsigned int jj = 0; jj < N; jj++)
				{
					y = a[jj][j];
					z = a[jj][i];
					//a[jj][j] = y * c + z * s;
					a[jj][j] = sumOfProducts( y, c, z, s );
					//a[jj][i] = z * c - y * s;
					a[jj][i] = diffOfProducts( z, c, y, s );
				}
			}
			assert( l >= 0 );
			assert( k >= 0 );
			rv1[l] = 0.0;
			rv1[k] = f;
			w[k] = x;
		}
		if (its == MAX_ITERATIONS)
		{
			fprintf( stderr,
				"\n%s: No convergence after %d iterations(SVD)\n", rtn, MAX_ITERATIONS );
			err = ERR_SVD_NO_CONVERGENCE;
			goto endfunc;
		}

	}

endfunc:
	free_vector( &rv1 );
	return err;
}
