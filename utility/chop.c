/*****************************************************************************
 * chop: utility to read an ASCII Input file and reproduce a
 * specified number of lines on standard output.
 * 
 * Usage
 * -----
 * chop [-h] [-s startline] [-n number of lines] [-S skip] [file]
 *
 * Synopsis
 * --------
 * The first two command line arguments specify the first line of the input
 * file to reproduce, and the number of subsequent lines.  The third argument
 * gives a skip between lines of output that are reproduced.
 *
 * Can be used as a filter.
 * If number of lines not specified, read through until EOF.
 * Lines are assumed to be BUFSIZ characters long at most.
 *
 * @file utility/chop.c
 * @ingroup group_utility
 *****************************************************************************/
/* Copyright (c) 1990+, Hugh M Blackburn */

#include <stdio.h>
#include <stdlib.h>

static void chop (FILE* fp, int s, int n, int S);


int main (int argc, char *argv[])
/* ------------------------------------------------------------------------- *
 * Wrapper for chop(), which does the work.  Here we do administration.
 * ------------------------------------------------------------------------- */
{
  static char usage[] = 
    "usage: chop [options] [input]\n"
    "  options:\n"
    "  -h         ... display this message\n"
    "  -n <lines> ... reproduce this many lines of file    [Default: to EOF]\n"
    "  -s <line>  ... start at this line number            [Default: 1]\n"
    "  -S <num>   ... skip <num> lines between each output [Default: 1]\n";
  int   i, start = 1, nlines = 0, skip = 1;
  FILE* fp;
  
  while (--argc && **++argv == '-') {
    switch (*++argv[0]) {
    case 'h':
      fprintf (stderr, usage);
      exit    (EXIT_SUCCESS);
      break;
    case 'n':
      if (*++argv[0]) nlines = atoi (*argv);
      else { --argc; nlines = atoi (*++argv); }
      break;
    case 's':
      if (*++argv[0]) start = atoi (*argv);
      else { --argc;  start = atoi (*++argv); }
      break;
    case 'S':
      if (*++argv[0]) skip = atoi (*argv);
      else { --argc; skip = atoi (*++argv); }
      break;
    default:
      fprintf (stderr, "chop: unknown arg %s\n", *argv);
      exit    (EXIT_FAILURE);
    }
  }

  if (argc) {			/* -- Input from list of named files. */
    for (i = 0; i < argc; i++)
      if (!(fp = fopen (argv[i], "r"))) {
	fprintf (stderr, "chop: can't open %s\n",argv[i]);
	fprintf (stderr, usage);
	exit    (EXIT_FAILURE);
      } else {
	chop    (fp, start, nlines, skip);
	fclose  (fp);
      }

  } else {			/* -- Input from stdin. */
    chop  (stdin, start, nlines, skip);
    while ((i = getchar()) != EOF);
  }
  
  return EXIT_SUCCESS;
}


static void chop (FILE* file , 
		  int   start,
		  int   nline,
		  int   skip )
/* ------------------------------------------------------------------------- *
 * This does the real work.
 * ------------------------------------------------------------------------- */
{
  char buf[BUFSIZ], *OK;
  int  n = 0;

  if (start > 1)
    while (start > 1 && (OK = fgets (buf, BUFSIZ, file))) {
      if (!OK) {
	fprintf (stderr, "chop: Reached EOF before start line.\n");
	exit    (EXIT_FAILURE); 
      }
      --start;
    }

  if (nline) {
    while (nline && (OK = fgets (buf, BUFSIZ, file))) {
      if (!OK) {
	fprintf (stderr, "chop: Reached EOF prematurely.\n") ;
	exit    (EXIT_FAILURE);
      }
      if (n % skip == 0) { fputs (buf, stdout);  --nline; }
      ++n;
    }

  } else
    while (fgets (buf, BUFSIZ, file)) {
      if (n % skip == 0) fputs (buf, stdout);
      ++n;
    }
}




