/*****************************************************************************
 * xfill:   x[i] = alpha.
 *****************************************************************************/

#include <cfemdef.h>

#if defined(__uxp__)
#pragma global novrec
#pragma global noalias
#endif


void dfill (int_t n, double alpha, double* x, int_t incx)
{
   int_t i;

  x += (incx<0) ? (-n+1)*incx : 0;

  for (i = 0; i < n; i++) x[i*incx] = alpha;
}


void ifill (int_t n, int_t alpha, int_t* x, int_t incx)
{
   int_t i;

  x += (incx<0) ? (-n+1)*incx : 0;

  for (i = 0; i < n; i++) x[i*incx] = alpha;
}


void sfill (int_t n, float alpha, float* x, int_t incx)
{
   int_t i;

  x += (incx<0) ? (-n+1)*incx : 0;

  for (i = 0; i < n; i++) x[i*incx] = alpha;
}
