///////////////////////////////////////////////////////////////////////////////
// meshplot.cpp: make a PostScript plot of a semtex/prism mesh from
// data in a file produced by the meshpr utility or similar.  Graphics
// are based on the PSplot package introduced in Numerical Recipes 3e
// Ch. 22 and described in webnote 26, http://wwww.nr.com/webnotes?26.
//
// Copyright (c) 2019 <--> $Date$, Hugh Blackburn
//
// Usage:
// ------
// meshplot [options] [file]
// options:
//   -h        ... display this message
//   -i        ... show element-internal structure
//   -d 'prog' ... call prog to display PostScript output
//   -o 'file' ... write output to named file [Default: stdout]
//
// Files:
// ------
// Input consists of mesh information in ASCII format (produced by
// meshpr), with a single-line header followed by (x,y) mesh locations
// and optionally by z location list (which is ignored). File may be
// read on stdin.
//
// Mesh header is of form:
//
// 9 9 32 206 NR NS NZ NEL
//
// then follows (x,y) locations of each element's 2D mesh points in
// row-major order, ASCII format.  In all local implementations,
// NR=NS=np.  Finally the NZ + 1 z locations of the 3D mesh are
// supplied, if NZ>1 -- these are ignored here.  
//
// $Id$
///////////////////////////////////////////////////////////////////////////////

#include <cstdlib>		// C standard headers.
#include <cstdio>
#include <cmath>
#include <cfloat>

#include <vector>		// C++ standard headers.
#include <iostream>
#include <fstream>

#include <cfemdef.h>		// Semtex headers.
#include <utility.h>

using namespace std;

#include "psplot.h"

static char prog[] = "meshplot";
static void getargs  (int, char**, char*&, char*&, bool&, istream*&);
static void readmesh (istream&, int&, int&, int&, int&,
		      vector<double>&, vector<double>&,
		      double&, double&, double&, double&);


int main (int    argc,
	  char** argv)
// ---------------------------------------------------------------------------
// An A4 page is 595 pts wide and 841 pts high, i.e. centred near 300, 420.
// ---------------------------------------------------------------------------
{
  char           *ofile = NULL, *disprog = NULL;
  istream*       input;
  bool           internal = false;
  int            i, j, k;
  int            nr, ns, nrns, nel, ntot;
  vector<double> x, y;
  double         xmin = FLT_MAX, ymin = FLT_MAX;
  double         xmax = FLT_MIN, ymax = FLT_MIN;
  double         xcen, ycen, xlen, ylen, lmax;

  getargs (argc, argv, ofile, disprog, internal, input);

  readmesh (*input, nr, ns, nel, ntot, x, y, xmin, xmax, ymin, ymax);
  
  nrns = nr * ns;
  xlen = xmax - xmin;
  ylen = ymax - ymin;
  xcen = 0.5 * (xmin + xmax);
  ycen = 0.5 * (ymin + ymax);
  lmax = max (xlen, ylen);

  // -- Have all the info, now make the plot.

  double xlo, xhi, ylo, yhi;	// -- Limits on the page in pts.
                                // -- Maximum dimension 500 pts.

  if (xlen > ylen) {
    xlo = 300. - 250.;             xhi = 300. + 250.;
    ylo = 420. - ylen/xlen * 250.; yhi = 420. + ylen/xlen * 250.;
  } else {
    xlo = 300. - xlen/ylen * 250.; xhi = 300. + xlen/ylen * 250.;
    ylo = 420. - 250.;             yhi = 420. + 250.;
  }
  
  PSpage pg ("mesh.ps");
  
  PSplot plot (pg, xlo, xhi, ylo, yhi);

#if 1
  plot.setlimits (xmin, xmax, ymin, ymax);
#else
  plot.setlimits (xcen - 0.525 * xlen, xcen + 0.525 * xlen,
		  ycen - 0.525 * ylen, ycen + 0.525 * ylen);
#endif
  
  plot.frame();
  plot.autoscales();

  if (internal) {  // -- Element internal structure, draw first.
    plot.setgray      (0.5);
    plot.setlinewidth (0.5);

    for (i = 0; i < nel; i++) {
      for (j = 1; j < (ns - 1); j++) {
	// -- Lines of constant s.
	k = i * nrns + j * nr;
	plot.polyline (nr, &x[k],  1,  &y[k], 1);
      }
      for (j = 1; j < (nr - 1); j++) {
	// -- Lines of constant r.
	k = i * nrns + j;
	plot.polyline (nr, &x[k], -nr, &y[k], -nr);
      }
    }
  }

  // -- Element outlines.

  plot.setgray      (0.);
  plot.setlinewidth (1.);

  for (i = 0; i < nel; i++) {
    // -- Side 1.
    k = i * nrns;
    plot.polyline (nr, &x[k],  1,  &y[k], 1);
    // -- Side 2.
    k = i * nrns + (nr - 1);
    plot.polyline (ns, &x[k],  nr, &y[k], nr);
    // -- Side 3.
    k = i * nrns + (nr - 1) * ns;
    plot.polyline (nr, &x[k], -1,  &y[k], -1);
    // -- Side 4.
    k = i * nrns;
    plot.polyline (nr, &x[k], -nr, &y[k], -nr);
  }

  if (disprog) plot.display (disprog);

  return (EXIT_SUCCESS);
}


static void getargs (int       argc    ,
		     char**    argv    ,
		     char*&    ofile   ,
		     char*&    disprog ,
		     bool&     internal,
		     istream*& input   )
// ---------------------------------------------------------------------------
// Deal with command-line arguments.
// ---------------------------------------------------------------------------
{
  char usage[] = "Usage: meshplot [options] [file]\n"
    "  options:\n"
    "  -h        ... display this message\n"
    "  -i        ... show element-internal mesh\n"
    "  -d <prog> ... call prog to display PostScript output\n"
    "  -o <file> ... write output to named file [Default: stdout]\n";
 
  while (--argc && **++argv == '-')
    switch (*++argv[0]) {
    case 'h':
      cout << usage;
      exit (EXIT_SUCCESS);
      break;
    case 'i':
      internal = true;
      break;
    case 'd':
      --argc;
      disprog = *++argv;
      break;
    case 'o':
      --argc;
      ofile = *++argv;
      break;
    default:
      cerr << usage;
      exit (EXIT_FAILURE);
      break;
    }

  if (argc == 1) {
    input = new ifstream (*argv);
    if (input -> fail()) message (prog, "unable to open geometry file", ERROR);
  } else input = &cin;
}


static void readmesh (istream&        file,
		      int&            nr  ,
		      int&            ns  ,
		      int&            nel ,
		      int&            ntot,
		      vector<double>& x   ,
		      vector<double>& y   ,
		      double&         xmin,
		      double&         xmax,
		      double&         ymin,
		      double&         ymax)
// ---------------------------------------------------------------------------
// File is already open; extract mesh info.
// ---------------------------------------------------------------------------
{
  const char routine[] = "readmesh";
  char  buf[STR_MAX], err[STR_MAX];
  int   i, nz;
 
  file >> nr >> ns >> nz >> nel;
  file.getline (buf, StrMax);

  if (!strstr (buf, "NR NS NZ NEL")) {
    sprintf (err, "mesh header line should include NR NS NZ NEL: %s", buf);
    message (routine, err, ERROR);
  }

  ntot = nr * ns * nel;
  x.resize (ntot);
  y.resize (ntot);

  for (i = 0; i < ntot; i++) {
    file >> x[i] >> y[i];
    xmin = min (xmin, x[i]);
    xmax = max (xmax, x[i]);
    ymin = min (ymin, y[i]);
    ymax = max (ymax, y[i]);
  }

  if (file.fail()) message (prog, "problem reading mesh", ERROR);
}

