/*****************************************************************************
 * xzero:  x[i] = 0.
 *****************************************************************************/

#include <string.h>
#include <cfemdef.h>

#if defined(__uxp__)
#pragma global novrec
#pragma global noalias
#endif

void dzero (int_t n, double* x, int_t incx)
{
#if defined(DEBUG)
   int_t i;

  x += (incx < 0) ? (-n+1)*incx : 0;

  for (i = 0; i < n; i++) x[i*incx] = 0.0;
#else
  if (incx == 1)
    memset (x, '\0', n * sizeof (double));

  else {
     int_t i;

    x += (incx < 0) ? (-n+1)*incx : 0;

    for (i = 0; i < n; i++) x[i*incx] = 0.0;
  }
#endif
}


void izero (int_t n, int_t* x, int_t incx)
{
  if (incx == 1)
    memset(x, '\0', n * sizeof (int_t));

  else {
     int_t i;

    x += (incx < 0) ? (-n+1)*incx : 0;

    for (i = 0; i < n; i++) x[i*incx] = 0;
  }
}


void szero (int_t n, float* x, int_t incx)
{
  if (incx == 1)
    memset(x, '\0', n * sizeof (float));

  else {
     int_t i;

    x += (incx < 0) ? (-n+1)*incx : 0;

    for (i = 0; i < n; i++) x[i*incx] = 0.0F;
  }
}
