/*****************************************************************************
 * xvneg:  y[i] = -x[i].
 *****************************************************************************/

#include <cfemdef.h>

#if defined(__uxp__)
#pragma global novrec
#pragma global noalias
#endif

void dvneg (int_t n, const double* x, int_t incx,
                             double* y, int_t incy)
{
   int_t i;

  x += (incx<0) ? (-n+1)*incx : 0;
  y += (incy<0) ? (-n+1)*incy : 0;

  for (i = 0; i < n; i++) y[i*incy] = -x[i*incx];
}


void ivneg (int_t n, const int_t* x, int_t incx,
                             int_t* y, int_t incy)
{
   int_t i;

  x += (incx<0) ? (-n+1)*incx : 0;
  y += (incy<0) ? (-n+1)*incy : 0;

  for (i = 0; i < n; i++) y[i*incy] = -x[i*incx];
}


void svneg (int_t n, const float* x, int_t incx,
                             float* y, int_t incy)
{
   int_t i;
  
  x += (incx<0) ? (-n+1)*incx : 0;
  y += (incy<0) ? (-n+1)*incy : 0;

  for (i = 0; i < n; i++) y[i*incy] = -x[i*incx];
}
