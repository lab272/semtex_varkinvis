/*****************************************************************************
 * fourier.c
 *
 * Copyright (c) 1999+, Hugh M Blackburn
 *
 * 1D Fourier transform routines for real data fields based on
 * FFTPACK, Temperton FFT routines (default), or vendor-supplied
 * alternatives.  NB: different restrictions may apply to input args
 * depending on compile-time options/FFT selection.

 * 1. Input data is to be Fourier transformed in the direction normal
 * to the most rapid traverse through memory, with sucessive points in
 * the transform separated by ntrn.  Data has a total of tlen * ntrn
 * real points.
 *
 * 2. If input parameter sign = +1 then the transform is takesn as
 * real-->complex, the outcome is normalised by 1/tlen, and a positive
 * sign is used in the complex exponential employed for the Fourier
 * coefficients. If sign = -1 then transform is taken as
 * complex-->real and a 1/2PI normalising factor is used.  In general:
 * these signs and scaling definitions match those adopted in
 * Numerical Recipes.
 *
 * 3. The number of transforms, ntrn, must be even owing to the way
 * real-complex transforms are implemented in the default FFT (two
 * real vectors treated as a single complex one).
 *
 * 4. The transform length tlen also must be even, owing to the way
 * Fourier data are assumed to be packed (with Nyquist data packed in
 * as the imaginary part of mode 0).  The transform length is assumed
 * to stay fixed during execution.
 *
 *****************************************************************************/

#include <stdio.h>
#include <stdlib.h>

#include <cfemdef.h>
#include <cveclib.h>
#include <cfemlib.h>


#if defined(_SX)  /* -- Use NEC FFT routines if available. */

static int_t*  ifax;
static double* Wtab;


void preFFT (const int_t tlen)
/* -------------------------------------------------------------------------
 * Carry out FFT rule checks, set up FFT data. 
 * ------------------------------------------------------------------------- */
{
  const char   routine[] = "preFFT";
  char         err[STR_MAX];
  static int_t tlenLast;

  if (tlen < 1 || tlen & 1) {
    sprintf (err, "transform length (%1d) must be even, and positive", tlen);
    message (routine, err, ERROR);
  }

  if (tlen != tlenLast) {	/* -- Setup required. */

    if (!ifax) ifax = ivector (0, 63); /* -- Length is always 64. */

    if (Wtab) freeDvector (Wtab, 0); 
    Wtab = dvector (0, tlen - 1);

    rftfax (tlen, ifax, Wtab);
  
    if (ifax[0] == -99) {
      sprintf (err, "transform length (%1d) needs prime factors 2, 3, 5", tlen);
      message (routine, err, ERROR);
    }

    tlenLast = tlen;
  }
}


void dDFTr (double*     data,
	    const int_t tlen,
	    const int_t ntrn,
	    const int_t sign)
/* ------------------------------------------------------------------------- *
 * Carry out multiple 1D real--complex Fourier transforms of data.
 * ------------------------------------------------------------------------- */
{
  const char   routine[] = "dDFTr";
  char         err[STR_MAX];
  static int_t tlenLast, ntotLast;
  const int_t  ntot = tlen * ntrn;

  if (tlen < 2 || !ntrn) return;

  if (ntrn & 1) {
    sprintf (err, "number of transforms (%1d) must be even", ntrn);
    message (routine, err, ERROR);
  }

  if (tlen != tlenLast) preFFT (tlen);

  if (ntot != ntotLast) {
    if (work) freeDvector (work, 0);
    work = dvector (0, ntot - 1);
  }

  if (sign == FORWARD) {
    rfft  (data, work, Wtab, ifax, tlen, ntrn,  1.0 / tlen);
    dcopy ((tlen - 2) * ntrn, data +              ntrn, 1, work,            1);
    dcopy (             ntrn, data + (tlen - 1) * ntrn, 1, data + ntrn,     1);
    dcopy ((tlen - 2) * ntrn, work,                     1, data + 2 * ntrn, 1);
  } else {
    dcopy ((tlen - 2) * ntrn, data + 2 * ntrn, 1, work,                     1);
    dcopy (             ntrn, data + ntrn,     1, data + (tlen - 1) * ntrn, 1);
    dcopy ((tlen - 2) * ntrn, work,            1, data +              ntrn, 1);
    rfft  (data, work, Wtab, ifax, tlen, ntrn, -1.0);
  }

  tlenLast = tlen;
  ntotLast = ntot;
}


#elif defined(DEBUG_FFT) /* -- FFTPACK routines (for checking). */

static double* Wtab;
static int_t*  ifax;


void preFFT (const int_t tlen)
/* -------------------------------------------------------------------------
 * Carry out FFT rule checks, set up FFT data.  FFTPACK real<-->complex
 * FFT routines do not impose any limitations on tlen but we require
 * it to be even owing to how dDFTr assumes complex data are packed.
 * ------------------------------------------------------------------------- */
{
  const char   routine[] = "preFFT";
  char         err[STR_MAX];
  static int_t tlenLast;

  if (tlen < 1 || tlen & 1) {
    sprintf (err, "transform length (%1d) must be even, and positive", tlen);
    message (routine, err, ERROR);
  }

  if (tlen != tlenLast) {	       /* -- Setup required.      */
    if (!ifax) ifax = ivector (0, 14); /* -- Length is always 15. */

    if (Wtab) freeDvector (Wtab, 0); 
    Wtab = dvector (0, tlen - 1);

    drffti (tlen, Wtab, ifax);

    tlenLast = tlen;
  }
}


void dDFTr (double*     data,
	    const int_t tlen,
	    const int_t ntrn,
	    const int_t sign)
/* ------------------------------------------------------------------------- *
 * Carry out multiple 1D real--complex Fourier transforms of data.
 * ------------------------------------------------------------------------- */
{
  const char     routine[] = "dDFTr";
  char           err[STR_MAX];
  static int_t   tlenLast;
  static double* work;
  double*        ptr;
  int_t          i;

  if (tlen < 2 || !ntrn) return;

  if (ntrn & 1) {
    sprintf (err, "number of transforms (%1d) must be even", ntrn);
    message (routine, err, ERROR);
  }

  if (tlen != tlenLast) {
    preFFT (tlen);
    if (work) freeDvector (work, 0);
    work = dvector (0, 2 * tlen - 1);
  }

  if (sign == FORWARD) {
    for (i = 0, ptr = data; i < ntrn; i++, ptr++) {
      dcopy  (tlen, ptr, ntrn, work, 1);
      drfftf (tlen, work, work + tlen, Wtab, ifax);
      dcopy  (tlen - 2, work + 1, 1, ptr + 2 * ntrn, ntrn);
      ptr[0]    = work[0];
      ptr[ntrn] = work[tlen - 1];
    }
    dscal (ntot, 1.0 / tlen, data, 1);
  } else {
    for (i = 0, ptr = data; i < ntrn; i++, ptr++) {
      work[tlen - 1] = ptr[ntrn];
      work[0]        = ptr[0];
      dcopy  (tlen - 2, ptr + 2 * ntrn, ntrn, work + 1, 1);
      drfftb (tlen, work, work + tlen, Wtab, ifax);
      dcopy  (tlen, work, 1, ptr, ntrn);
    }
  }

  tlenLast = tlen;
}


#else /* -- Temperton FFT routine, which is the default. */

static int_t   ip, iq, ir, ipqr2;
static double* Wtab;


void preFFT (const int_t tlen)
/* ------------------------------------------------------------------------- *
 * Carry out FFT rule checks, set up FFT data.
 * ------------------------------------------------------------------------- */
{
  const char   routine[] = "preFFT";
  char         err[STR_MAX];
  static int_t tlenLast;
  int_t        chk;

  if (tlen < 1 || tlen & 1) {
    sprintf (err, "transform length (%1d) must be even, and positive", tlen);
    message (routine, err, ERROR);
  }

  if (tlen != tlenLast) {	/* -- Setup required. */
    chk = tlen;
    prf235 (&chk, &ip, &iq, &ir, &ipqr2);
  
    if (!chk) {
      sprintf (err, "transform length (%1d) needs prime factors 2, 3, 5", tlen);
      message (routine, err, ERROR);
    }

    if (Wtab) freeDvector (Wtab, 0); 
    Wtab = dvector (0, ipqr2 - 1);

    dsetpf (Wtab, tlen, ip, iq, ir);

    tlenLast = tlen;
  }
}


void dDFTr (double*     data,
	    const int_t tlen,
	    const int_t ntrn,
	    const int_t sign)
/* ------------------------------------------------------------------------- *
 * Carry out multiple 1D real--complex Fourier transforms of data.
 * ------------------------------------------------------------------------- */
{
  const char    routine[] = "dDFTr";
  char          err[STR_MAX];
  static int_t  tlenLast, ntotLast;
  const int_t   ntot = tlen * ntrn;
  static double *work;

  if (tlen < 2 || !ntrn) return;

  if (ntrn & 1) {
    sprintf (err, "number of transforms (%1d) must be even", ntrn);
    message (routine, err, ERROR);
  }

  if (tlen != tlenLast) preFFT (tlen);

  if (ntot != ntotLast) {
    if (work) freeDvector (work, 0);
    work = dvector (0, ntot - 1);
  }

  dmpfft (data, work, ntrn, tlen, ip, iq, ir, Wtab, sign);
  if (sign == FORWARD) dscal (ntot, 1.0 / tlen, data, 1);

  tlenLast = tlen;
  ntotLast = ntot;
}

#endif
