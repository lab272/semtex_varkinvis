/*****************************************************************************
 * xvsub:  z[i] = x[i] - y[i].
 *****************************************************************************/

#include <cfemdef.h>

#if defined(__uxp__)
#pragma global novrec
#pragma global noalias
#endif

void dvsub (int_t n,
	    const double* x, int_t incx,
	    const double* y, int_t incy,
	          double* z, int_t incz)
{
   int_t i;
  
  x += (incx<0) ? (-n+1)*incx : 0;
  y += (incy<0) ? (-n+1)*incy : 0;
  z += (incz<0) ? (-n+1)*incz : 0;

  for (i = 0; i < n; i++) z[i*incz] = x[i*incx] - y[i*incy];
}


void ivsub (int_t n,
	    const int_t* x, int_t incx,
	    const int_t* y, int_t incy,
	          int_t* z, int_t incz)
{
   int_t i;

  x += (incx<0) ? (-n+1)*incx : 0;
  y += (incy<0) ? (-n+1)*incy : 0;
  z += (incz<0) ? (-n+1)*incz : 0;

  for (i = 0; i < n; i++) z[i*incz] = x[i*incx] - y[i*incy];
}


void svsub (int_t n,
	    const float* x, int_t incx,
	    const float* y, int_t incy,
	          float* z, int_t incz)
{
   int_t i;
  
  x += (incx<0) ? (-n+1)*incx : 0;
  y += (incy<0) ? (-n+1)*incy : 0;
  z += (incz<0) ? (-n+1)*incz : 0;

  for (i = 0; i < n; i++) z[i*incz] = x[i*incx] - y[i*incy];
}
