/*****************************************************************************
 * y = (double) x[i]
 *****************************************************************************/

#include <cfemdef.h>

#if defined(__uxp__)
#pragma global novrec
#pragma global noalias
#endif

void vdble (int_t n, const float *x, int_t incx, double *y, int_t incy)
{
   int_t i;

  x += (incx<0) ? (-n+1)*incx : 0;
  y += (incy<0) ? (-n+1)*incy : 0;

  for (i = 0; i < n; i++) y[i*incy] = (double) x[i*incx];
}
