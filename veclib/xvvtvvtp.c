/*****************************************************************************
 * xvvtvvtp:  z[i] = (v[i] * w[i]) + (x[i] * y[i]).
 *
 * $Id$
 *****************************************************************************/

#ifdef __uxp__
#pragma global novrec
#pragma global noalias
#endif


void dvvtvvtp (int n, const double* v, int incv,
        	      const double* w, int incw,
	              const double* x, int incx,
	              const double* y, int incy,
	                    double* z, int incz)
{
  register int i;

  if (incv == 1 && incw == 1 && incx == 1 && incy == 1 && incz == 1) 
   for (i = 0; i < n; i++) 
     z[i] = v[i] * w[i] + x[i] * y[i]; 

  else {

    w += (incw<0) ? (-n+1)*incw : 0;
    x += (incx<0) ? (-n+1)*incx : 0;
    y += (incy<0) ? (-n+1)*incy : 0;
    z += (incz<0) ? (-n+1)*incz : 0;

    for (i = 0; i < n; i++)
      z[i*incz] = v[i*incv] * w[i*incw] + x[i*incx] * y[i*incy];
  }
}


void svvtvvtp (int n, const float* v, int incv,
	              const float* w, int incw,
	              const float* x, int incx,
	              const float* y, int incy,
	                    float* z, int incz)
{
  register int i;

  if (incv == 1 && incw == 1 && incx == 1 && incy == 1 && incz == 1) 
   for (i = 0; i < n; i++) 
     z[i] = v[i] * w[i] + x[i] * y[i]; 

  else {

    w += (incw<0) ? (-n+1)*incw : 0;
    x += (incx<0) ? (-n+1)*incx : 0;
    y += (incy<0) ? (-n+1)*incy : 0;
    z += (incz<0) ? (-n+1)*incz : 0;

    for (i = 0; i < n; i++)
      z[i*incz] = v[i*incv] * w[i*incw] + x[i*incx] * y[i*incy];
  }
}