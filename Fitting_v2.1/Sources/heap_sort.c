/*****************************************************************
 *
 * File........:	heap_sort.c
 * Function....:	sorting of values
 * Author......:	Tilo Strutz
 * last changes:	20.10.2007
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
#include <math.h>

/*---------------------------------------------------------------
 *	heap_sort_d()
 * sorting of values in values[0...N-1] in ascending order
 * values[] is replaced on output by sorted values
 *--------------------------------------------------------------*/
void
heap_sort_d( unsigned long N, double values[])
{
  unsigned long i, ir, j, l;
  double rvalues;
	
  if (N < 2)
		return;
 
	l = (N >> 1);
  ir = N - 1;
	
  for (;;)
  {
		if (l > 0)
		{
			l--;
			rvalues = values[l]; 
		}
		else
		{
			rvalues = values[ir];
			values[ir] = values[0];
			ir--;
			if (ir == 0)
			{
				values[0] = rvalues;
				break;
			}
		}
		i = l;
		
		j = l + l + 2;
		while (j <= ir + 1)
		{
			if (j < ir + 1 && values[j - 1] < values[j])
			{
				j++;
			}
			if (rvalues < values[j - 1])
			{
				values[i] = values[j - 1];
				i = j - 1;
				j <<= 1;
			}
			else
			{
				j = ir + 2;	/* terminate while-loop */
			}
		}
		values[i] = rvalues;
  }
}

/*---------------------------------------------------------------
 *	heap_sort_d_()
 * sorting of values in values[0...N-1] in ascending order
 * values[] is replaced on output by sorted values
 *--------------------------------------------------------------*/
void
heap_sort_d_( unsigned long N, double values[], long idx[])
{
  unsigned long i, ir, j, l;
  double rvalues;
	int iidx;
	
  if (N < 2)
		return;
 
	for ( i = 0; i < N; i++) idx[i] = i;

	l = (N >> 1);
  ir = N - 1;
	
  for (;;)
  {
		if (l > 0)
		{
			l--;
			rvalues = values[l]; iidx = idx[l];
		}
		else
		{
			rvalues = values[ir];			iidx = idx[ir];
			values[ir] = values[0];		idx[ir] = idx[0];
			ir--;
			if (ir == 0)
			{
				values[0] = rvalues;		idx[0] = iidx;
				break;
			}
		}
		i = l;
		
		j = l + l + 2;
		while (j <= ir + 1)
		{
			if (j < ir + 1 && values[j - 1] < values[j])
			{
				j++;
			}
			if (rvalues < values[j - 1])
			{
				values[i] = values[j - 1];	idx[i] = idx[j - 1];
				i = j - 1;
				j <<= 1;
			}
			else
			{
				j = ir + 1 + 1;	/* terminate while-loop */
			}
		}
		values[i] = rvalues;		idx[i] = iidx;
  }
}
