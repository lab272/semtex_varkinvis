/*****************************************************************************
 * icount:  number of non-zero values in x.
 * imatch:  number of occurences of integer value in x.
 *****************************************************************************/

#include <cfemdef.h>


int_t icount (int_t n, const int_t* x, int_t incx)
{
  int_t i, sum = 0;

  x += (incx<0) ? (-n+1)*incx : 0;

  for (i = 0; i < n; i++ ) sum += (x[i*incx]) ? 1 : 0;

  return sum;
}


int_t imatch (int_t n, const int_t alpha, const int_t* x, int_t incx)
{
  int_t i, sum = 0;

  x += (incx<0) ? (-n+1)*incx : 0;

  for (i = 0; i < n; i++ ) sum += (x[i*incx] == alpha) ? 1 : 0;

  return sum;
}
