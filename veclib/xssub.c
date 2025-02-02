/*****************************************************************************
 * xssub:  y[i] = alpha - x[i].
 *****************************************************************************/

#include <cfemdef.h>

#if defined(__uxp__)
#pragma global novrec
#pragma global noalias
#endif

void dssub (int_t n, double alpha, const double* x, int_t incx,
	                                 double* y, int_t incy)
{
   int_t i;

  x += (incx<0) ? (-n+1)*incx : 0;
  y += (incy<0) ? (-n+1)*incy : 0;

  for (i = 0; i < n; i++) y[i*incy] = alpha - x[i*incx];
}


void isub (int_t n, int_t alpha, const int_t* x, int_t incx,
	                               int_t* y, int_t incy)
{
   int_t i;

  x += (incx<0) ? (-n+1)*incx : 0;
  y += (incy<0) ? (-n+1)*incy : 0;

  for (i = 0; i < n; i++) y[i*incy] = alpha - x[i*incx];
}


void sssub (int_t n, float alpha, const float* x, int_t incx,
	                                float* y, int_t incy)
{
   int_t i;

  x += (incx<0) ? (-n+1)*incx : 0;
  y += (incy<0) ? (-n+1)*incy : 0;

  for (i = 0; i < n; i++) y[i*incy] = alpha - x[i*incx];

}
