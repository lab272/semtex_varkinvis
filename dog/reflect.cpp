///////////////////////////////////////////////////////////////////////////////
// reflect.cpp: reflect a field dump defined on a half-mesh onto a full
// mesh.  Apply sign change to appropriate velocity components, unless
// explicitly supressed by setting revpar == true (command line).
//
// If negation is called for, the velocity component that would cross
// the reflection boundary is not negated, but the other two are.
//
// Built from semtex/utility/data2df_template.cpp.  See also
// dog/symmetrise.cpp and flipmap.cpp
//
// USAGE
// -----
// reflect [options] -m mapfile fieldfile
// options:
// -h ... print this message.
// -r ... reverse parity
//
// If fieldfile is not present, read from standard input.  Write to
// standard output.  Note that the size of the output file will be
// about twice as large as the input file and the associated fields
// will have twice as many elements.  It is assumed that the mapfile
// has a structure appropriate to the task of the reflection mapping
// on the full domain (mapfile produced by dog/flipmap).
//
// A "hidden constraint" with this plan of attack is that the layout
// of the first half of elements used to generate the flip map has to
// agree with the layout used on the starting mesh.  This means that
// the design/layout of the two session files cannot be independent.
//
// Copyright (c) 2010+, Hugh M Blackburn
///////////////////////////////////////////////////////////////////////////////

#include <sem.h>
#include <data2df.h>

static char prog[] = "reflect";
static void getargs  (int, char**, bool&, istream*&, istream*&);
static void loadmap  (Header&, istream&, char&, vector<int_t>&,vector<int_t>&);
static void halfread (istream*, Data2DF*);


int main (int    argc,
	  char** argv)
// ---------------------------------------------------------------------------
// Driver.
// ---------------------------------------------------------------------------
{
  int_t            i, j;
  istream          *input = 0, *mapping = 0;
  Header           h;
  vector<Data2DF*> u;
  Data2DF*         tmp;
  char             generator = 'x';
  vector<int_t>    positive, negative;
  bool             revpar = false;                       // -- Reverse parity.

  Femlib::init ();
  
  getargs (argc, argv, revpar, mapping, input);

  *input >> h;
  h.nel += h.nel;		// -- Double the number of elements.
  loadmap (h, *mapping, generator, positive, negative);

  cout << h;

  u.resize (h.nFields());
  tmp = new Data2DF (h.nr, h.nz, h.nel, 0);

  for (i = 0; i < h.nFields(); i++) {

    // - Create u[i].

    u[i] = new Data2DF (h.nr, h.nz, h.nel, h.flds[i]);

    const int_t ntot = u[i] -> _ntot;

    // -- Read in u[i] and byte-swap if necessary.

    halfread (input, u[i]);	// -- Input file only has 1/2 the elements.
    if (h.swab()) u[i] -> reverse();

    // -- Create in tmp a spatial reflection of what was just read.

    (*tmp = *u[i]) . reflect2D (positive, negative);

    // -- Copy the reflection of the other half of the data field,
    // -- negating appropriate velocity vector components.

    if (generator == 'y')
      if (revpar) {
	if      (u[i] -> getName() == 'u')
	  for (j = 0; j < ntot; j++)
	    u[i]->_data[j] -= (fabs(u[i]->_data[j]) < EPSDP) 
	      ? tmp->_data[j] : 0.0;
	else if (u[i] -> getName() == 'v')
	  for (j = 0; j < ntot; j++)
	    u[i]->_data[j] += (fabs(u[i]->_data[j]) < EPSDP) 
	      ? tmp->_data[j] : 0.0;
	else if (u[i] -> getName() == 'w')
	  for (j = 0; j < ntot; j++)
	    u[i]->_data[j] -= (fabs(u[i]->_data[j]) < EPSDP) 
	      ? tmp->_data[j] : 0.0;
	else if (u[i] -> getName() == 'p')
	  for (j = 0; j < ntot; j++)
	    u[i]->_data[j] += (fabs(u[i]->_data[j]) < EPSDP) 
	      ? tmp->_data[j] : 0.0;
      }	else {
	if      (u[i] -> getName() == 'u')
	  for (j = 0; j < ntot; j++)
	    u[i]->_data[j] += (fabs(u[i]->_data[j]) < EPSDP) 
	      ? tmp->_data[j] : 0.0;
	else if (u[i] -> getName() == 'v')
	  for (j = 0; j < ntot; j++)
	    u[i]->_data[j] -= (fabs(u[i]->_data[j]) < EPSDP) 
	      ? tmp->_data[j] : 0.0;
	else if (u[i] -> getName() == 'w')
	  for (j = 0; j < ntot; j++)
	    u[i]->_data[j] += (fabs(u[i]->_data[j]) < EPSDP) 
	      ? tmp->_data[j] : 0.0;
	else if (u[i] -> getName() == 'p')
	  for (j = 0; j < ntot; j++)
	    u[i]->_data[j] += (fabs(u[i]->_data[j]) < EPSDP) 
	      ? tmp->_data[j] : 0.0;
      }
    else if (generator == 'x')
      if (revpar) {
	if      (u[i] -> getName() == 'u')
	  for (j = 0; j < ntot; j++)
	    u[i]->_data[j] += (fabs(u[i]->_data[j]) < EPSDP) 
	      ? tmp->_data[j] : 0.0;
	else if (u[i] -> getName() == 'v')
	  for (j = 0; j < ntot; j++)
	    u[i]->_data[j] -= (fabs(u[i]->_data[j]) < EPSDP) 
	      ? tmp->_data[j] : 0.0;
	else if (u[i] -> getName() == 'w')
	  for (j = 0; j < ntot; j++)
	    u[i]->_data[j] -= (fabs(u[i]->_data[j]) < EPSDP) 
	      ? tmp->_data[j] : 0.0;
	else if (u[i] -> getName() == 'p')
	  for (j = 0; j < ntot; j++)
	    u[i]->_data[j] += (fabs(u[i]->_data[j]) < EPSDP) 
	      ? tmp->_data[j] : 0.0;
      }	else {
	if      (u[i] -> getName() == 'u')
	  for (j = 0; j < ntot; j++)
	    u[i]->_data[j] -= (fabs(u[i]->_data[j]) < EPSDP) 
	      ? tmp->_data[j] : 0.0;
	else if (u[i] -> getName() == 'v')
	  for (j = 0; j < ntot; j++)
	    u[i]->_data[j] += (fabs(u[i]->_data[j]) < EPSDP) 
	      ? tmp->_data[j] : 0.0;
	else if (u[i] -> getName() == 'w')
	  for (j = 0; j < ntot; j++)
	    u[i]->_data[j] += (fabs(u[i]->_data[j]) < EPSDP) 
	      ? tmp->_data[j] : 0.0;
	else if (u[i] -> getName() == 'p')
	  for (j = 0; j < ntot; j++)
	    u[i]->_data[j] += (fabs(u[i]->_data[j]) < EPSDP) 
	      ? tmp->_data[j] : 0.0;
      }

    // -- Output u[i];

    cout << *u[i];
  }

  return EXIT_SUCCESS;
}


static void getargs (int       argc  ,
		     char**    argv  ,
		     bool&     revpar,
		     istream*& mapfl ,
		     istream*& input )
// ---------------------------------------------------------------------------
// Deal with command-line arguments.  Based on dog/reflect.cpp.
// ---------------------------------------------------------------------------
{
  char usage[] = "Usage: reflect [options] -m mapfile [file]\n"
    "options:\n"
    "-h ... print this message\n"
    "-r ... reverse parity\n";

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
    case 'r':
      revpar = true;
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
    Veclib::alert (prog, "map file doesn't conform with session file",   ERROR);
  file >> generator;
  if (!(generator == 'x' || generator == 'y'))
    Veclib::alert (prog, "symmetry generator must be either 'x' or 'y'", ERROR);
  
  file >> NMAP;

  positive.resize (NMAP);
  negative.resize (NMAP);

  for (i = 0; i < NMAP; i++) file >> positive[i] >> negative[i];

  if (!file)
    Veclib::alert (prog, "bad (premature end of?) map file", ERROR);
}


static void halfread (istream* file,
		      Data2DF* dat )
// ---------------------------------------------------------------------------
// This is like the input operator for Data2DF but we only read in the
// first half of the data -- since the input data are defined for a
// domain which is only half as big as the output. Leave zero
// everything else (as set by constructor).
// ---------------------------------------------------------------------------
{
  int_t       i;
  const int_t ntot = dat -> _np2 * dat -> _nel / 2;

  for (i = 0; i < dat -> _nz; i++)
    file -> read (reinterpret_cast<char*> (dat -> _plane[i]),
		  ntot * sizeof (real_t));
}
