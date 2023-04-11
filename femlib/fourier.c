/*****************************************************************************
 * fourier.c
 *
 * Copyright (c) 1999+, Hugh M Blackburn
 *
 * 1D Fourier transform routines for real data fields based on FFTPACK
 * or Temperton FFT routines, or vendor-supplied alternatives.
 * NB: different restrictions may apply to input args depending on
 * selected routine.
 *
 * See the FFT codes themselves for the definition of a discrete
 * Fourier transform.  If input parameter sign = 1 then a positive
 * sign is used in the complex exponential employed for the Fourier
 * coefficients and a 1/2PI normalizing factor is used in IDFT.  In
 * general: the signs and scaling definitions match those adopted in
 * Numerical Recipes.
 *
 *****************************************************************************/

#include <stdio.h>
#include <stdlib.h>

#include <cfemdef.h>
#include <cveclib.h>
#include <cfemlib.h>

void dDFTr (double*     data,
	    const int_t tlen,
	    const int_t ntrn,
	    const int_t sign)
/* ------------------------------------------------------------------------- *
 * Carry out multiple 1D single--complex Fourier transforms of data.
 * Data is to be Fourier transformed in the direction normal to the most
 * rapid traverse through memory, with sucessive points in the transform
 * separated by ntrn.  Data has a total of tlen * ntrn real points.
 *
 * Input parameters:
 * ----------------
 * tlen: number of real data in each transform, product of prime numbers.
 * ntrn: number of transforms to perform (also skip in data).
 * sign: transform direction: +1 ==> r-->c, -1 ==> c-->r.
 *
 * Notes:
 * -----
 *
 * (1) Data are scaled/normalized with 1/tlen when sign is +1, so that
 * the zeroth Fourier mode contains the spatial average value.
 
 * (2) After forward (r-->c) transform, data are ordered so that
 * within each transform, the zeroth mode datum comes first.  The
 * zeroth mode is followed by the real datum from the maximum
 * frequency (Nyquist) mode, after which the real and imaginary parts
 * for each mode alternate.
 *
 * ------------------------------------------------------------------------- */
{
  const char  routine[] = "dDFTr";
  char        err[STR_MAX];


  if (tlen < 2 || !ntrn) return;

#if defined(_SX)  /* -- Use NEC FFT routines. */

  const int_t ntot = tlen * ntrn;
  int_t       *ifax;
  double      *work, *Wtab;

  ifax = ivector (0, 63);
  work = dvector (0, ntot + tlen - 1);
  Wtab = work + ntot;

  rftfax (tlen, ifax, Wtab);

  if (ifax[0] == -99)
    message (routine, "tlen needs prime factors 2, 3, 5", ERROR);

  if (sign == FORWARD) {
    rfft  (data, work, Wtab, ifax, tlen, ntrn, 1.0 / tlen);
    dcopy ((tlen - 2) * ntrn, data +              ntrn, 1, work,            1);
    dcopy (             ntrn, data + (tlen - 1) * ntrn, 1, data + ntrn,     1);
    dcopy ((tlen - 2) * ntrn, work,                     1, data + 2 * ntrn, 1);
  } else {
    dcopy ((tlen - 2) * ntrn, data + 2 * ntrn, 1, work,                     1);
    dcopy (             ntrn, data + ntrn,     1, data + (tlen - 1) * ntrn, 1);
    dcopy ((tlen - 2) * ntrn, work,            1, data +              ntrn, 1);
    rfft  (data, work, Wtab, ifax, tlen, ntrn, -1.0);
  }

  freeIvector (ifax, 0);
  freeDvector (work, 0);

#elif defined(DEBUG_FFT) /* -- Unvectorized FFTPACK routines (for checking). */

  const int_t ntot = tlen * ntrn;
  int_t       i, *ifax;
  double      *work, *ptr, *Wtab;

  ifax = ivector (0, 14);
  work = dvector (0, 3 * tlen - 1);
  Wtab = work + 2 * tlen;
  ptr  = data;

  drffti (tlen, Wtab, ifax);

  switch (sign) {

  case FORWARD:
    for (i = 0; i < ntrn; i++, ptr++) {
      dcopy  (tlen, ptr, ntrn, work, 1);
      drfftf (tlen, work, work + tlen, Wtab, ifax);
      dcopy  (tlen - 2, work + 1, 1, ptr + 2 * ntrn, ntrn);
      ptr[0]    = work[0];
      ptr[ntrn] = work[tlen - 1];
    }
    dscal (ntot, 1.0 / tlen, data, 1);
    break;

  case INVERSE:
    for (i = 0; i < ntrn; i++, ptr++) {
      work[tlen - 1] = ptr[ntrn];
      work[0]        = ptr[0];
      dcopy  (tlen - 2, ptr + 2 * ntrn, ntrn, work + 1, 1);
      drfftb (tlen, work, work + tlen, Wtab, ifax);
      dcopy  (tlen, work, 1, ptr, ntrn);
    }
    break;

  default:
    message (routine, "illegal direction flag", ERROR);
    break;
  }
  
  freeDvector (work, 0);
  freeIvector (ifax, 0);

#else /* -- Temperton FFT routine is default. */

  static int_t  nlast;
  const int_t   ntot = tlen * ntrn;
  int_t         nwork, dum, ip, iq, ir, ipqr2;  
  static double *work;
  double        *Wtab;


  if (ntrn & 1) {
    sprintf (err, "number of transforms (%1d) must be even", ntrn);
    message (routine, err, ERROR);
  }

  dum = tlen;
  prf235 (&dum, &ip, &iq, &ir, &ipqr2);
  
  if (!dum) {
    sprintf (err, "transform length (%1d) needs prime factors 2, 3, 5", tlen);
    message (routine, err, ERROR);
  }

  nwork = ntot + ipqr2;

  if (nwork != nlast) {
    if (work) freeDvector (work, 0);
    work = dvector (0, nwork - 1);      /* -- Only reallocate if needed. */
  }

  Wtab = work + ntot;

  dsetpf (Wtab, tlen, ip, iq, ir);
  dmpfft (data, work, ntrn, tlen, ip, iq, ir, Wtab, sign);
  if (sign == FORWARD) dscal (ntot, 1.0 / tlen, data, 1);

  nlast = nwork;

#endif
}
