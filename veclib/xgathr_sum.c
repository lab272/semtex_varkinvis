/*****************************************************************************
 *  xgathr_sum:  vector gather, with summation:  z[i] += x[y[i]].
 *
 * $Id$
 *****************************************************************************/

#include <cfemdef.h>

#if defined(__uxp__)
#pragma global novrec
#pragma global noalias
#endif

void dgathr_sum (integer n, const double* x, const integer* y, double* z)
{
  register integer i;

  for (i = 0; i < n; i++) z[i] += x[y[i]];
}


void igathr_sum (integer n, const integer* x, const integer* y, integer* z)
{
  register integer i;

  for (i = 0; i < n; i++) z[i] += x[y[i]];
}


void sgathr_sum (integer n, const float* x, const integer* y, float* z)
{
  register integer i;

  for (i = 0; i < n; i++) z[i] += x[y[i]];
}
