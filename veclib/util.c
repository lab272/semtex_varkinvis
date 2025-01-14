/*****************************************************************************
 *                   FILE, STRING, I-O, TIMING UTILITIES
 *****************************************************************************/

/* #include <sys/types.h> */ /* No longer required? */

#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <string.h>

#include <cfemdef.h>
#include <cveclib.h>


int_t   _vecIreg[NVREG];	/* For FORTRAN linkage. */
char    _vecCreg[NVREG];
float   _vecSreg[NVREG];
double  _vecDreg[NVREG];

char     buf[STR_MAX];		/* A string for general use. */


void message (const char *routine, const char *text, int level)
/* ------------------------------------------------------------------------- *
 * A simple error handler.
 * ------------------------------------------------------------------------- */
{
  switch (level) {
  case WARNING:
    fprintf (stderr, "-- WARNING: %s: %s\n", routine, text); 
    break;
  case ERROR:
    fprintf (stderr, "-- ERROR: %s: %s\n", routine, text); 
    break;
  case REMARK:
    fprintf (stderr, "-- REMARK: %s: %s\n", routine, text);
    break;
  default:
    fprintf (stderr, "bad error level in message: %d\n", level);
    exit (EXIT_FAILURE);
    break;
  }

  if (level == ERROR) exit (EXIT_FAILURE);
}


FILE *efopen (const char *file, const char *mode)
/* ------------------------------------------------------------------------- *
 * fopen file, die if can't.
 * ------------------------------------------------------------------------- */
{
  FILE *fp;

  if (fp = fopen (file, mode)) return fp;

  sprintf (buf, "can't open %s mode %s", file, mode);
  message ("efopen", buf, ERROR);
  
  return (FILE*) 0;
}


#if !defined(i860) && !defined(dclock)


double dclock (void)
/* ------------------------------------------------------------------------- *
 * Double-precision timing routine.
 * ------------------------------------------------------------------------- */
{
  static double tps = 1.0 / CLOCKS_PER_SEC;
  return (double) clock() * tps;
}


float sclock (void)
/* ------------------------------------------------------------------------- *
 * Single-precision timing routine.
 * ------------------------------------------------------------------------- */
{
  static float tps = 1.0F / CLOCKS_PER_SEC;
  return (float) clock() * tps;
}

#endif


static void stringSubst2 (char * const       bsrcptr,
			  char * const       bdstptr,
			  const char * const pattstr,
			  const char * const replstr,
			  const int          pattsiz,
			  const int          replsiz)
/* ------------------------------------------------------------------------- *
 * String substitution routine, borrowed from Scotch. This part is hidden.
 * ------------------------------------------------------------------------- */
{
  char *              pattptr;
  int                 pattidx;

  /* Search for the pattern in the remaining source string   */
  pattptr = strstr (bsrcptr, pattstr);

  /* Get length of unchanged part */
  pattidx = (pattptr == NULL) ? (strlen (bsrcptr) + 1): (pattptr - bsrcptr);

  /* If replacement is smaller, pre-move unchanged part */
  if (replsiz < pattsiz)
    memmove (bdstptr, bsrcptr, pattidx * sizeof (char));

  /* If remaining part of string has to be processed */
  if (pattptr != NULL)
    stringSubst2 (pattptr + pattsiz,
		  bdstptr + pattidx + replsiz,
		  pattstr, replstr,
		  pattsiz, replsiz);

  /* If replacement is longer, post-move unchanged part */
  if (replsiz > pattsiz)
    memmove (bdstptr, bsrcptr, pattidx * sizeof (char));
  
  /* If there is something to replace         */
  if (pattptr != NULL)
    /* Write replacement string */
    memcpy (bdstptr + pattidx, replstr, replsiz * sizeof (char));

  return;
}


void stringSubst (char * const       buffptr, /* String to search into */
		  const char * const pattstr, /* Pattern to search for */
		  const char * const replstr) /* Replacement string    */
/* ------------------------------------------------------------------------- *
 * String substitution routine, borrowed from Scotch.  This part is visible.
 * ------------------------------------------------------------------------- */ 
{
  int pattsiz;
  int replsiz;

  pattsiz = strlen (pattstr);
  replsiz = strlen (replstr);

  stringSubst2 (buffptr, buffptr, pattstr, replstr, pattsiz, replsiz);
}


void printDvector (FILE  *fp     ,
		   int_t  width  ,
		   int_t  prec   ,
		   int_t  ntot   ,
		   int_t  inc    ,
		   int_t  nfield , ...)
/* ------------------------------------------------------------------------- *
 * Print up a variable number of dvectors on fp, in columns.
 * ------------------------------------------------------------------------- */
{
  char    routine[] = "printDvector";
  int_t i, j, k;
  double  **u;
  va_list ap;

  u = (double **) calloc (nfield, sizeof (double*));
  va_start (ap, nfield);
  for (i = 0; i < nfield; i++) u[i] = va_arg (ap, double*); 
  va_end (ap);

  k = (inc < 0) ? (-ntot + 1)*inc : 0;

  for (i = 0; i < ntot; i++) {
    for (j = 0; j < nfield; j++)
      if (fprintf (fp, "%*.*f", width, prec, u[j][k]) < 0)
	message (routine, "unable to write to file", ERROR);
    k += inc;
    fprintf(fp, "\n");
  }

  free (u);
}


void printIvector (FILE  *fp     ,
		   int_t    width  ,
		   int_t    ntot   ,
		   int_t    inc    ,
		   int_t    nfield , ...)
/* ------------------------------------------------------------------------- *
 * Print up a variable number of ivectors on fp, in columns.
 * ------------------------------------------------------------------------- */
{
  char    routine[] = "printIvector";
  int_t i, j, k;
  int_t **u;
  va_list ap;

  u = (int_t **) calloc (nfield, sizeof (int_t*));
  va_start (ap, nfield);
  for (i = 0; i < nfield; i++) u[i] = va_arg (ap, int_t*); 
  va_end (ap);

  k = (inc < 0) ? (-ntot + 1)*inc : 0;

  for (i = 0; i < ntot; i++) {
    for (j = 0; j < nfield; j++)
      if (fprintf (fp, "%*d", width, u[j][k]) < 0)
	message (routine, "couldn't write to file", ERROR);
    k += inc;
    fprintf (fp, "\n");
  }

  free (u);
}


void printSvector (FILE  *fp     ,
		   int_t width  ,
		   int_t prec   ,
		   int_t ntot   ,
		   int_t inc    ,
		   int_t nfield , ...)
/* ------------------------------------------------------------------------- *
 * Write (ASCII) a variable number of svectors on fp, in columns.
 * ------------------------------------------------------------------------- */
{
  char    routine[] = "printSvector";
  int_t   i, j, k;
  float   **u;
  va_list ap;

  u = (float **) calloc (nfield, sizeof (float*));
  va_start (ap, nfield);
  for (i = 0; i < nfield; i++) u[i] = va_arg (ap, float*); 
  va_end (ap);

  k = (inc < 0) ? (-ntot + 1)*inc : 0;

  for (i = 0; i < ntot; i++) {
    for (j = 0; j < nfield; j++)
      if (fprintf (fp, "%*.*f", width, prec, u[j][k]) < 0)
	message (routine, "unable to write to file", ERROR);
    k += inc;
    fprintf (fp, "\n");
  }

  free (u);
}
