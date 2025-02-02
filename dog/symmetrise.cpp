///////////////////////////////////////////////////////////////////////////////
// symmetrise.C: Read in a data file header and its data (no geometric
// information) and enforce reflection symmetry using node lists
// created by flipmap.  See also drive.cpp.
//
// USAGE
// -----
// symmetrise [options] -m mapfile [file]
// options:
// -h       ... print this message.
//
// If file is not present, read from standard input.  Write to
// standard output.  Mapfile contains the node lists generated by
// flipmap, and also gives the one-character generator of the symmetry
// operation, see below:
//
// If generator == 'x' then we have a refection in the y axis, i.e.
//    u(-x,y) = -u(x,y),  v(-x,y) =  v(x,y),   w(-x,y) = w(x,y)
// If generator == 'y' then we have a refection in the x axis, i.e.
//    u(x,-y) =  u(x,y),  v(x,-y) = -v(x,-y),  w(-x,y) = w(x,y)
// If generator == 'd' then we have reflections in both x and y axes, i.e.
//    u(-x,-y) =  -u(x,y),  v(-x,-y) = -v(x,y),  w(-x,-y) = w(x,y)
//
// NB: It seems that the symmetries above could only be correct for a
// 2D flow field (or the real part of a mode).
//
// Copyright (c) 2008+, Hugh M Blackburn
///////////////////////////////////////////////////////////////////////////////

static char RCS[] = "$Id$";

#include <sem.h>
#include <data2df.h>

static char prog[] = "symmetrise";
static void getargs  (int, char**, istream*&, istream*&);
static void loadmap  (Header&, istream&, char&, vector<int_t>&,vector<int_t>&);


int main (int    argc,
	  char** argv)
// ---------------------------------------------------------------------------
// Driver.
// ---------------------------------------------------------------------------
{
  int_t            i, nf;
  istream          *mapping = 0, *input = 0;
  Header           h;
  vector<Data2DF*> u;
  Data2DF*         tmp;
  char             generator = 'x';
  vector<int_t>    positive, negative;
  real_t           mean2d;

  Femlib::init ();
  
  getargs (argc, argv, mapping, input);

  *input >> h;
  loadmap (h, *mapping, generator, positive, negative);

  cout   << h;
  u.resize (h.nFields());
  tmp = new Data2DF (h.nr, h.nz, h.nel, 0);

  for (i = 0; i < h.nFields(); i++) {

    // -- Create u[i].

    u[i] = new Data2DF (h.nr, h.nz, h.nel, h.flds[i]);

    // -- Read in u[i] and byte-swap if necessary.

    *input >> *u[i];
    if (h.swab()) u[i] -> reverse();

    // -- Make symmetric/antisymmetric parts of selected velocity components.

    (*tmp = *u[i]) . reflect2D (positive, negative); 

    if (generator == 'y') {
      if      (u[i] -> getName() == 'u') *u[i] += *tmp;
      else if (u[i] -> getName() == 'v') *u[i] -= *tmp;
      else if (u[i] -> getName() == 'w') *u[i] += *tmp;
    } else if (generator == 'x') {
      if      (u[i] -> getName() == 'u') *u[i] -= *tmp;
      else if (u[i] -> getName() == 'v') *u[i] += *tmp;
      else if (u[i] -> getName() == 'w') *u[i] += *tmp;
    } else if (generator == 'd') {
      if      (u[i] -> getName() == 'u') *u[i] -= *tmp;
      else if (u[i] -> getName() == 'v') *u[i] -= *tmp;
      else if (u[i] -> getName() == 'w') *u[i] += *tmp;

#if 0
      // -- Enforce zero mean value.

      mean2d = Veclib::sum (u[i] -> _ntot, u[i] -> _plane[0], 1)/u[i] -> _ntot;

      Veclib::sadd (u[i] -> _ntot, -mean2d,
		    u[i] -> _plane[0], 1, u[i] -> _plane[0], 1);
#endif
    }

    *u[i] *= 0.5;
    
    // -- Output u[i];

    cout << *u[i];
  }

  return EXIT_SUCCESS;
}


static void getargs (int       argc ,
		     char**    argv ,
		     istream*& mapfl,
		     istream*& input)
// ---------------------------------------------------------------------------
// Deal with command-line arguments.
// ---------------------------------------------------------------------------
{
  char usage[] = "Usage: symmetrise [options] -m mapfile [file]\n"
    "options:\n"
    "-h       ... print this message\n";
  
  if (argc < 3) {
    cerr << usage ;
    exit (EXIT_FAILURE);
  }
    
  while (--argc && **++argv == '-')
    switch (*++argv[0]) {
    case 'h':
      cout << usage;
      exit (EXIT_SUCCESS);
      break;
    case 'm':
      --argc; ++argv;
      mapfl = new ifstream (*argv);
      if (mapfl -> bad()) Veclib::alert (prog, "unable to open map file",ERROR);
      break;
    default:
      cerr << usage;
      exit (EXIT_FAILURE);
      break;
    }

  if (argc == 1) {
    input = new ifstream (*argv);
    if (input -> bad()) Veclib::alert (prog, "unable to open input file",ERROR);
  } else input = &cin;
}


static void loadmap (Header&        headr    ,
		     istream&       file     ,
		     char&          generator,
		     vector<int_t>& positive ,
		     vector<int_t>& negative )
// ---------------------------------------------------------------------------
// Load symmetry mapping information from file.
// ---------------------------------------------------------------------------
{
  const int_t np  = headr.nr;
  const int_t nel = headr.nel;
  char        buf[StrMax], err[StrMax];
  int_t       i, NR, NS, NEL, NMAP;
  
  if (!file) {
    sprintf (err, "cannot find map file %s", buf);
    Veclib::alert (prog, err, ERROR);
  }

  file >> NR >> NS >> NEL >> NEL;
  file.ignore (StrMax, '\n');

  if (NR != np || NS != np || NEL != nel)
    Veclib::alert (prog, "map file doesn't conform with session file", ERROR);
  file >> generator;
  if (!(generator == 'x' || generator == 'y' || generator == 'd'))
    Veclib::alert (prog,
		   "symmetry generator must be either 'x', 'y', or 'd'", ERROR);
  
  file >> NMAP;

  positive.resize (NMAP);
  negative.resize (NMAP);

  for (i = 0; i < NMAP; i++) file >> positive[i] >> negative[i];

  if (!file)
    Veclib::alert (prog, "bad (premature end of?) map file", ERROR);

}


static void mirror (real_t* tgt)
// ---------------------------------------------------------------------------
// Apply RT-flip-map. Note that to avoid holes in the mapping, the
// gather and scatter vectors have to be used in the order shown.
// ---------------------------------------------------------------------------
{
#if 0
  int_t          i, k;
  const int_t    ND = Geometry::nPert();
  const int_t    NP = Geometry::planeSize();
  const int_t    NZ = Geometry::nZ();
  const int_t    NM = positive.size();
  static real_t* tmp;

  if (!tmp) tmp = new real_t [NP];

  // -- First, the reflection.

  for (i = 0; i < ND; i++)
    for (k = 0; k < NZ; k++) {
      Veclib::copy (NP, tgt + (i*NZ+k)*NP, 1, tmp, 1);
      Veclib::gathr_scatr (NM,tmp,&negative[0],&positive[0],tgt + (i*NZ+k)*NP);
    }
  
  // -- Then the sign change.

  if (generator == 'x')		// -- Change sign of 'u'.
    for (k = 0; k < NZ; k++)
      Veclib::neg (NP, tgt + (0*NZ+k)*NP, 1);
  else				// -- Change sign of 'v'.
    for (k = 0; k < NZ; k++)
      Veclib::neg (NP, tgt + (1*NZ+k)*NP, 1);
#endif
}
