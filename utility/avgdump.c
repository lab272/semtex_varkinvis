/*****************************************************************************
 * avgdump.c:  Compute (running) averages of field files.
 *
 * Copyright (c) 1999 Hugh Blackburn
 *
 * SYNOPSIS
 * --------
 * Input is two field files.  The first ("old.file") is assumed to
 * contain a running average of previous dumps, while the second
 * ("new.file") is assumed to contain new data that will be added into
 * the running average.  Both files must be binary (but not
 * necessarily matching) format and storage (number of elements,
 * interpolation orders, scalar fields) must conform.  The step number
 * in "old.file" is taken as the number of averages that have
 * contributed to it, and "new.file"'s data is added in with
 * appropriate weight.  Averaging is initialised using the
 * command-line flag "-i", which adds the data in the two files with
 * equal weight.  Only the first dump in each file is dealt with.
 *
 * A new binary file is written to stdout.  The step number is set to
 * reflect the number of averages that have been done to the data.
 *
 * USAGE
 * -----
 * avgdump [options] old.file new.file
 * options:
 * -h ... print this message
 * -i ... initialise averaging
 *
 * $Id$
 *****************************************************************************/

#include <stdarg.h>
#include <stdlib.h>
#include <malloc.h>
#include <stdio.h>
#include <time.h>
#include <string.h>

#include <femdef.h>
#include <veclib.h>

static char  prog[]    = "avgdump";
static char* hdr_fmt[] = {	 /* -- Header output formatting. */
  "%-25s "             "Session\n",
  "%-25s "             "Created\n",
  "%-5d%-5d%-5d%-10d " "Nr, Ns, Nz, Elements\n",
  "%-25d "             "Step\n",
  "%-25.6g "           "Time\n",
  "%-25.6g "           "Time step\n",
  "%-25.6g "           "Kinvis\n",
  "%-25.6g "           "Beta\n",
  "%-25s "             "Fields written\n",
  "%-25s "             "Format\n"
};
typedef struct {		 /* -- Data information structure. */
  char     session [STR_MAX];
  char     creation[STR_MAX];
  int      np               ;
  int      nz               ;
  int      nel              ;
  int      step             ;
  double   time             ;
  double   timestep         ;
  double   kinvis           ;
  double   beta             ;
  char     field [STR_MAX]  ;
  char     format[STR_MAX]  ;
  double** data             ;
} Dump;


static void getargs   (int, char**, FILE**, FILE**, int*);
static void getheader (FILE*, Dump*);
static void getdata   (FILE*, Dump*);
static void runavg    (Dump*, Dump*, int);
static void printup   (FILE*, Dump*);


int main (int    argc,
	  char** argv)
/* ------------------------------------------------------------------------- *
 * Driver routine.
 * ------------------------------------------------------------------------- */
{
  FILE *oldfile = 0, *newfile = 0;
  Dump *olddump = 0, *newdump = 0;
  int  init = 0;

  getargs (argc, argv, &oldfile, &newfile, &init);

  olddump = (Dump*) calloc (1, sizeof (Dump));
  newdump = (Dump*) calloc (1, sizeof (Dump));

  getheader (oldfile, olddump);
  getheader (newfile, newdump);

  if (olddump -> np  != newdump -> np  ||
      olddump -> nz  != newdump -> nz  ||
      olddump -> nel != newdump -> nel)
    message (prog, "structure of files don't match",           ERROR);

  if (!strstr (olddump -> field, newdump -> field))
    message (prog, "average fields don't match dumped fields", ERROR);

  getdata (oldfile, olddump);
  getdata (newfile, newdump);
  runavg  (olddump, newdump, init);
  printup (stdout,  olddump);

  return (EXIT_SUCCESS);
}


static void getargs (int    argc   ,
		     char** argv   ,
		     FILE** oldfile,
		     FILE** newfile,
		     int*   init   )
/* ------------------------------------------------------------------------- *
 * Parse command-line arguments.
 * ------------------------------------------------------------------------- */
{
  char usage[] =
    "usage: avgdump [options] old.file new.file\n"
    "options:\n"
    "  -h ... display this message\n"
    "  -i ... initialise averaging\n";
    
  char err[STR_MAX], c;

  while (--argc && **++argv == '-')
    switch (c = *++argv[0]) {
    case 'h':
      fprintf (stderr, usage); exit (EXIT_SUCCESS); break;
    case 'i':
      *init = 1; break;
    default:
      sprintf (err, "illegal option: %c\n", c);
      message (prog, err, ERROR);
      break;
    }

  if (argc == 2) {
    *oldfile = efopen (argv[0], "r");
    *newfile = efopen (argv[1], "r");
  } else {
    fprintf (stderr, usage); exit (EXIT_FAILURE);
  }  
}


static void getheader (FILE* f,
		       Dump* h)
/* ------------------------------------------------------------------------- *
 * Find header information.
 * ------------------------------------------------------------------------- */
{
  char buf[STR_MAX];

  fgets  (h -> session,  STR_MAX, f);
  fgets  (h -> creation, STR_MAX, f);
  fscanf (f, "%d %*s %d %d", &h -> np, &h -> nz, &h -> nel);
  fgets  (buf, STR_MAX, f);
  fscanf (f, "%d", &h -> step);
  fgets  (buf, STR_MAX, f);
  fscanf (f, "%lf", &h -> time);
  fgets  (buf, STR_MAX, f);
  fscanf (f, "%lf", &h -> timestep);
  fgets  (buf, STR_MAX, f);
  fscanf (f, "%lf", &h -> kinvis);
  fgets  (buf, STR_MAX, f);
  fscanf (f, "%lf", &h -> beta);
  fgets  (buf, STR_MAX, f);
  fscanf (f, "%s", h -> field);
  fgets  (buf, STR_MAX, f);
  fgets  (h -> format, STR_MAX, f);

  if (!strstr (h -> format,      "binary"))
    message (prog, "input field file not in binary format",     ERROR);
  else if (!strstr (h -> format, "endian"))
    message (prog, "input field file in unknown binary format", WARNING);
}


static void getdata (FILE* f,
		     Dump* h)
/* ------------------------------------------------------------------------- *
 * Find the number of fields, allocate storage & do binary read.
 * ------------------------------------------------------------------------- */
{
  char      localfmt[STR_MAX];
  int       i, swab;
  const int npts    = h -> np * h -> np * h -> nz * h -> nel;
  const int nfields = strlen (h -> field);
  const int ntot    = npts * nfields;

  h -> data = dmatrix (0, nfields - 1, 0, npts - 1);

  if (fread (h -> data[0], sizeof (double), ntot, f) != ntot) {
    sprintf (localfmt, "could not read fields: %s", h -> field);
    message (prog, localfmt, ERROR);
  }

  format (localfmt);
  swab = (strstr (h -> format, "big") && strstr (localfmt,    "little")) || 
         (strstr (localfmt,    "big") && strstr (h -> format, "little"));

  if (swab) dbrev (ntot, h -> data[0], 1, h -> data[0], 1);
}


static void runavg (Dump* a,
		    Dump* b,
		    int   init)
/* ------------------------------------------------------------------------- *
 * Running average of a & b, leave in a.  Optionally initialise.
 * ------------------------------------------------------------------------- */
{
  int       i;
  double    fac;
  const int nfields = strlen (a -> field);
  const int npts    = a -> np * a -> np * a -> nz * a -> nel;
 
  if (init) {
    fac = 0.5;
    a -> step = 2;
    for (i = 0; i < nfields; i++)
      dsvvpt (npts, fac, a -> data[i], 1, b -> data[i], 1, a -> data[i], 1);

  } else {
    fac = (double) a -> step;
    a -> step++;
    for (i = 0; i < nfields; i++) {
      dsvtvp (npts, fac, a -> data[i], 1, b -> data[i], 1, a -> data[i], 1);
      dsmul  (npts, 1.0 / (fac + 1.0),    a -> data[i], 1, a -> data[i], 1);
    }
  }
}


static void printup (FILE* f,
		     Dump* h)
/* ------------------------------------------------------------------------- *
 * Output the modified data.
 * ------------------------------------------------------------------------- */
{
  int       i;
  const int ntot = h -> np * h -> np * h -> nz * h -> nel;
  const int nfields = strlen (h -> field);
  char      s1[STR_MAX], s2[STR_MAX];
  time_t    tp = time (0);

  fprintf  (f, h -> session);
  strftime (s2, 25, "%a %b %d %H:%M:%S %Y", localtime (&tp));
  sprintf  (s1, hdr_fmt[1], s2);
  fprintf  (f, s1);
  sprintf  (s1, hdr_fmt[2], h -> np, h -> np, h -> nz, h -> nel);
  fprintf  (f, s1);
  sprintf  (s1, hdr_fmt[3], h -> step);
  fprintf  (f, s1);
  sprintf  (s1, hdr_fmt[4], h -> time);
  fprintf  (f, s1);
  sprintf  (s1, hdr_fmt[5], h -> timestep);
  fprintf  (f, s1);
  sprintf  (s1, hdr_fmt[6], h -> kinvis);
  fprintf  (f, s1);
  sprintf  (s1, hdr_fmt[7], h -> beta);
  fprintf  (f, s1);
  sprintf  (s1, hdr_fmt[8], h -> field);
  fprintf  (f, s1);
  sprintf  (s2, "binary "); format (s2 + strlen (s2));
  sprintf  (s1, hdr_fmt[9], s2);
  fprintf  (f, s1);

  for (i = 0; i < nfields; i++) 
    if (fwrite (h -> data[i], sizeof (double), ntot, f) != ntot)
      message (prog, "error occurred while writing", ERROR);
}
