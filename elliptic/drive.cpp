//////////////////////////////////////////////////////////////////////////////
// drive.cpp: compute solution to elliptic problem, optionally compare to
// exact solution (see getoptions(), below).
//
// USAGE:
// -----
// elliptic [options] session
//   options:
//   -h       ... print this message
//   -i       ... use iterative solver
//   -v[v...] ... increase verbosity level
//
// If session.frc is found, use this field file as forcing for the
// elliptic problem, otherwise use the 'forcing' string in the USER
// section; failing that, set forcing to zero.
//
// Author
// ------
// Hugh Blackburn
// Department of Mechanical & Aerospace Engineering
// Monash University
// Vic 3800
// Australia
// hugh.blackburn@monash.edu
//
// References
// ----------
// Blackburn, Lee, Albrecht & Singh (2019) "Semtex: a spectral
// element-Fourier solver for the incompressible Navier-Stokes
// equations in cylindrical or Cartesian coordinates", CPC 245:106804
//
// Copyright (c) 1994+, Hugh M Blackburn
//////////////////////////////////////////////////////////////////////////////

#include <sem.h>

#ifdef MPI
  static char prog[] = "elliptic_mp";
#else
  static char prog[] = "elliptic";
#endif

static void getargs    (int, char**, char*&);
static void getoptions (FEML*, char*&, char*&);
static void preprocess (const char*, FEML*&, Mesh*&, vector<Element*>&,
			BCmgr*&, Domain*&, AuxField*&);
static void getforcing (const char*, const char*, AuxField*);

void Helmholtz (Domain*, AuxField*);


int main (int    argc,
	  char** argv)
// ---------------------------------------------------------------------------
// Driver.
// ---------------------------------------------------------------------------
{
  char             *session, *forcefunc = 0, *exact = 0;
  vector<Element*> elmt;
  FEML*            file;
  Mesh*            mesh;
  BCmgr*           bman;
  Domain*          domain;
  AuxField*        forcefld;

  Femlib::init  ();
  Message::init (&argc, &argv);

  getargs (argc, argv, session);

  preprocess (session, file, mesh, elmt, bman, domain, forcefld);

  getoptions (file, forcefunc, exact);
  getforcing (session, forcefunc, forcefld);

  domain -> restart();

  Helmholtz (domain, forcefld);

  ROOTONLY if (exact) domain -> u[0] -> errors (mesh, exact);

  domain -> dump();

  Message::stop();

  return EXIT_SUCCESS;
}


static void getargs (int    argc   ,
		     char** argv   ,
		     char*& session)
// ---------------------------------------------------------------------------
// Install default parameters and options, parse command-line for optional
// arguments.  Last argument is name of a session file, not dealt with here.
// ---------------------------------------------------------------------------
{
  const char routine[] = "getargs";
  const char usage[] =
    "Usage: %s [options] session\n"
    "  options:\n"
    "  -h       ... print this message\n"
    "  -i       ... use iterative solver\n"
    "  -v[v...] ... increase verbosity level\n";
  char buf[StrMax];
 
  while (--argc  && **++argv == '-')
    switch (*++argv[0]) {
    case 'h':
      ROOTONLY {
	sprintf (buf, usage, prog);
	cout << buf;
      }
      exit (EXIT_SUCCESS);
      break;
    case 'i':
      Femlib::ivalue ("ITERATIVE", static_cast<int_t>(1));
      break;
    case 'v':
      do
	Femlib::ivalue ("VERBOSE", Femlib::ivalue("VERBOSE")+1);
      while (*++argv[0] == 'v');
      break;
    default:
      ROOTONLY { sprintf (buf, usage, prog); cout << buf; }
      exit (EXIT_FAILURE);
      break;
    }
  
  if (argc != 1) Veclib::messg (routine, "no session definition file", ERROR);

  session = *argv;
}


static void preprocess (const char*       session,
			FEML*&            file   ,
			Mesh*&            mesh   ,
			vector<Element*>& elmt   ,
			BCmgr*&           bman   ,
			Domain*&          domain ,
			AuxField*&        forcing)
// ---------------------------------------------------------------------------
// Create objects needed for execution, given the session file name.
// They are listed in order of creation.
// ---------------------------------------------------------------------------
{
  const int_t        verbose = Femlib::ivalue ("VERBOSE");
  Geometry::CoordSys space;
  int_t              i, np, nz, nel;

  // -- Initialise problem and set up mesh geometry.

  VERBOSE cout << "Building mesh ..." << endl;

  file = new FEML (session);
  mesh = new Mesh (file);

  VERBOSE cout << "done" << endl;

  // -- Set up global geometry variables.

  VERBOSE cout << "Setting geometry ... ";

  nel   =  mesh -> nEl();
  np    =  Femlib::ivalue ("N_P");
  nz    =  Femlib::ivalue ("N_Z");
  space = (Femlib::ivalue ("CYLINDRICAL")) ? 
    Geometry::Cylindrical : Geometry::Cartesian;
  
  Geometry::set (np, nz, nel, space);

  VERBOSE cout << "done" << endl;

  // -- Build all the elements.

  VERBOSE cout << "Building elements ... ";

  elmt.resize (nel);
  for (i = 0; i < nel; i++) elmt[i] = new Element (i, np, mesh);

  VERBOSE cout << "done" << endl;

  // -- Build all the boundary condition applicators.

  VERBOSE cout << "Building boundary condition manager ..." << endl;

  bman = new BCmgr (file, elmt);

  VERBOSE cout << "done" << endl;

  // -- Build the solution domain.

  VERBOSE cout << "Building domain ..." << endl;

  domain = new Domain (file, mesh, elmt, bman);

  VERBOSE cout << "done" << endl;

  // -- Build the forcing field.

  VERBOSE cout << "Building forcing ...";

  forcing = new AuxField(new real_t[(size_t)Geometry::nTotProc()],
			 Geometry::nZProc(), elmt, 'f');

  VERBOSE cout << "done" << endl;
}


static void getoptions (FEML*  feml ,
			char*& forcf,
			char*& exact)
// ---------------------------------------------------------------------------
// Try to load forcing function string and exact solution string from USER
// section of FEML file.  The section is not required to be present.
// 
// Expect something in the form:
// <USER>
// forcing 0
// exact   sin(TWOPI*x)*sinh(TWOPI*y)/sinh(TWOPI)
// </USER>
//
// Either or both of the two strings may be absent.
// ---------------------------------------------------------------------------
{
  char routine[] = "options";
  char s[StrMax];

  if (feml -> seek ("USER")) {
    feml -> stream().ignore (StrMax, '\n');

    while (feml -> stream() >> s) {
      if (strcmp (s, "</USER>") == 0) break;

      upperCase (s);
      if (strcmp (s, "FORCING") == 0)
	feml -> stream() >> (forcf = new char [StrMax]);
      else if (strcmp (s, "EXACT") == 0)
	feml -> stream() >> (exact = new char [StrMax]);
    }

    if (strcmp (s, "</USER>") != 0)
      Veclib::messg
	(routine, "couldn't sucessfully close <USER> section", ERROR);
  }
}


static void getforcing (const char* session  , 
			const char* forcefunc,
			AuxField*   forcefld )
// ---------------------------------------------------------------------------
// If file session.frc is found, use the contents (of the first field
// variable) to initialise forcefld, failing that use the string
// forcefunc, otherwise initialise to zero.  Finally, Fourier
// transform.
// ---------------------------------------------------------------------------
{
  const char routine[] = "getforcing";
  char       restartfile[StrMax];
  
  ROOTONLY cout << "-- Forcing                 : ";
  ifstream file (strcat (strcpy (restartfile, session), ".frc"));

  if (file) {
    ROOTONLY {
      cout << "read from file " << restartfile;
      cout.flush();
    }

    // -- Strip header and check the data conforms.

    int_t         np, nz, nel, ntot, nfields;
    int_t         npchk,  nzchk, nelchk, swab = 0;
    char          s[StrMax], f[StrMax];

    if (file.getline(s, StrMax).eof())
      Veclib::messg (routine, "forcing file is empty", ERROR);

    file.getline(s,StrMax).getline(s,StrMax);

    string        ss(s);
    istringstream sss (ss);

    sss >> np    >> np    >> nz    >> nel;

    forcefld -> describe (f);

    sss.clear();
    sss.str (ss = f);
    sss >> npchk >> npchk >> nzchk >> nelchk;

    if (np  != npchk )
      Veclib::messg (routine, "element size mismatch",       ERROR);
    if (nz  != nzchk )
      Veclib::messg (routine, "number of z planes mismatch", ERROR);
    if (nel != nelchk)
      Veclib::messg (routine, "number of elements mismatch", ERROR);
  
    ntot = np * np * nz * nel;
    if (ntot != Geometry::nTot())
      Veclib::messg (routine, "declared sizes mismatch", ERROR);

    file.getline(s,StrMax).getline(s,StrMax);
    file.getline(s,StrMax).getline(s,StrMax);
    file.getline(s,StrMax).getline(s,StrMax);

    nfields = 0; while (isalpha (s[nfields])) nfields++;
    if (!nfields)
      Veclib::messg
	(routine, "no fields declared in forcing", ERROR);

    file.getline (s, StrMax);
    Veclib::describeFormat (f);

    if (!strstr (s, "binary"))
      Veclib::messg
	(routine, "input field file not in binary format", ERROR);
  
    if (!strstr (s, "endian"))
      Veclib::messg
	(routine, "input field file in unknown binary format", WARNING);
    else {
      swab = ((strstr (s, "big") && strstr (f, "little")) ||
	      (strstr (f, "big") && strstr (s, "little")) );
      ROOTONLY {
	if (swab) cout << " (byte-swapping)";
	cout.flush();
      }
    }

    // -- Read data, byteswap if required. 

    file >> *forcefld;
    if (swab) forcefld -> reverse();
    file.close();

    forcefld -> transform (FORWARD);
    ROOTONLY forcefld -> zeroNyquist();

  } else if (forcefunc) {
    ROOTONLY cout << "set to string: " << forcefunc;
    (*forcefld = forcefunc).transform (FORWARD);
    ROOTONLY forcefld -> zeroNyquist();

  } else {
    ROOTONLY cout << "set to zero";
    *forcefld = 0.0;
  }

  ROOTONLY cout << endl;
}
