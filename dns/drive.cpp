//////////////////////////////////////////////////////////////////////////////
// drive.cpp: control spectral element DNS for incompressible flows.
//
// USAGE:
// -----
// dns [options] session
//   options:
//   -h       ... print usage prompt
//   -i       ... use iterative solver for viscous [and pressure] steps
//   -v[v...] ... increase verbosity level
//   -chk     ... turn off checkpoint field dumps [default: selected]
//   -S|C|N   ... regular skew-symm || convective || Stokes advection
//   -f       ... freeze velocity field (to advect scalar, only)
//
// AUTHOR:
// ------
// Hugh M Blackburn
// Department of Mechanical & Aerospace Engineering
// Monash University
// Vic 3800
// Australia
// hugh.blackburn@monash.edu
//
// REFERENCES
// ----------
// See the more extensive list of references in integrate.cpp.  A
// general reference for semtex is (same as [5] in integrate.cpp):
//
// Blackburn, Lee, Albrecht & Singh (2019) "Semtex: a spectral
// element--Fourier solver for the incompressible Navier--Stokes
// equations in cylindrical or Cartesian coordinates", CPC
// 245:106804.
//
// Copyright (c) 1994+, Hugh M Blackburn
//////////////////////////////////////////////////////////////////////////////

#include <dns.h>

#if defined(MPI_EX)
  static char prog[] = "dns_mp";
#else
  static char prog[] = "dns";
#endif

static void getargs    (int, char**, bool&, char*&);
static void preprocess (const char*, FEML*&, Mesh*&, vector<Element*>&,
			BCmgr*&, Domain*&, FieldForce*&);

void integrate (void (*)
		(Domain*, BCmgr*, AuxField**, AuxField**, FieldForce*),
		Domain*, BCmgr*, DNSAnalyser*, FieldForce*);
void AdvectDiffuse (Domain*, BCmgr*, DNSAnalyser*);


int main (int    argc,
	  char** argv)
// ---------------------------------------------------------------------------
// Driver.
// ---------------------------------------------------------------------------
{
#ifdef _GNU_SOURCE
  feenableexcept (FE_OVERFLOW);    // -- Force SIG8 crash on FP overflow.
#endif

  char*            session;
  int              nproc = 1, iproc = 0, npart2d = 1;  
  bool             freeze = false;
  vector<Element*> elmt;
  FEML*            file;
  Mesh*            mesh;
  BCmgr*           bman;
  Domain*          domain;
  DNSAnalyser*     analyst;
  FieldForce*      FF;

  Femlib::init  ();
  Message::init (&argc, &argv, nproc, iproc);

  Femlib::ivalue ("I_PROC", iproc);
  Femlib::ivalue ("N_PROC", nproc);

  getargs (argc, argv, freeze, session);

  preprocess (session, file, mesh, elmt, bman, domain, FF);

  if ((!domain -> hasScalar()) && freeze)
    Veclib::alert (prog, "need scalar declared if velocity is frozen", ERROR);

  analyst = new DNSAnalyser (domain, bman, file);

  domain -> restart ();

  ROOTONLY domain -> report ();
  
  if (freeze) 
    AdvectDiffuse (domain, bman, analyst); // -- Velocity field doesn't evolve.
  else {
    switch (Femlib::ivalue ("ADVECTION")) {
    case 0: integrate (   skewSymmetric,domain,bman,analyst,FF); break;
    case 1: integrate (altSkewSymmetric,domain,bman,analyst,FF); break;
    case 2: integrate (      convective,domain,bman,analyst,FF); break;
    case 3: integrate (     rotational1,domain,bman,analyst,FF); break;
    case 4: integrate (     rotational2,domain,bman,analyst,FF); break;
    case 5: integrate (          Stokes,domain,bman,analyst,FF); break;
    }
  }

  Message::stop ();

  return EXIT_SUCCESS;
}


static void getargs (int    argc   ,
		     char** argv   ,
		     bool&  freeze ,
		     char*& session)
// ---------------------------------------------------------------------------
// Install default parameters and options, parse command-line for optional
// arguments.  Last argument is name of a session file, not dealt with here.
// ---------------------------------------------------------------------------
{
  char       buf[StrMax];
  const char routine[] = "getargs";
  const char usage[]   = "Usage: %s [options] session-file\n"
    "  [options]:\n"
    "  -h       ... print this message\n"
    "  -f       ... freeze velocity field (scalar advection/diffusion only)\n"
    "  -i       ... use iterative solver for viscous steps\n"
    "  -v[v...] ... increase verbosity level\n"
    "  -chk     ... turn off checkpoint field dumps [default: selected]\n"
    "  -S|C|N   ... regular skew-symm || convective || Stokes advection\n";

  Femlib::ivalue ("ADVECTION", 1); // -- Default is alternating skew symmetric.

  while (--argc  && **++argv == '-')
    switch (*++argv[0]) {
    case 'h':
      sprintf (buf, usage, prog);
      cout << buf;
      exit (EXIT_SUCCESS);
      break;
    case 'S': Femlib::ivalue ("ADVECTION", 0); break;
    case 'C': Femlib::ivalue ("ADVECTION", 2); break;
    case 'N': Femlib::ivalue ("ADVECTION", 5); break;
    case 'f':
      freeze = true;
      break;
    case 'i':
      do			// -- Only allowing ITERATIVE=1 (Viscous).
	Femlib::ivalue ("ITERATIVE", 1);
      while (*++argv[0] == 'i');
      break;
    case 'v':
      do
	Femlib::ivalue ("VERBOSE",   Femlib::ivalue ("VERBOSE")   + 1);
      while (*++argv[0] == 'v');
      break;
    case 'c':
      if (strstr ("chk", *argv))     Femlib::ivalue ("CHKPOINT",    0);
      else { fprintf (stdout, usage, prog); exit (EXIT_FAILURE); }
      break;
    default:
      sprintf (buf, usage, prog);
      cout << buf;
      exit (EXIT_FAILURE);
      break;
    }

  if   (argc != 1) Veclib::alert (routine,
				   "no session definition file", ERROR);
  else             session = *argv;
}


static void preprocess (const char*       session,
			FEML*&            file   ,
			Mesh*&            mesh   ,
			vector<Element*>& elmt   ,
			BCmgr*&           bman   ,
			Domain*&          domain ,
			FieldForce*&      FF     )
// ---------------------------------------------------------------------------
// Create objects needed for execution, given the session file name.
// They are listed above in order of creation.
// ---------------------------------------------------------------------------
{
  const char routine[] = "preprocess";
  const int_t        verbose = Femlib::ivalue ("VERBOSE");
  Geometry::CoordSys space;
  int_t              i, np, nz, nel, procid, seed;
  int                npart2d = 1, ipart2d, npartz, ipartz;  

  // -- Initialise problem and set up mesh geometry.

  VERBOSE cout << "Building mesh ..." << endl;

  file = new FEML (session);
  mesh = new Mesh (file);

  Message::grid ((int) npart2d, ipart2d, npartz, ipartz);

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

  // -- If token RANSEED > 0 then initialize the random number
  //    generator based on wall clock time and process ID (i.e. a "truly"
  //    pseudo-random number).  NB: it is important to have done this
  //    here, before any other possible call to random number routines.

  if (Femlib::ivalue ("RANSEED") > 0) {
    procid = Geometry::procID();
    seed   = -abs((procid + 1) * (char) time(NULL));
  } else seed = -1;
  Veclib::ranInit (seed);

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

  // -- Parse Field Force info from FEML file.

  VERBOSE cout << "Building field forcing ..." << endl;

  FF = new FieldForce (domain, file);

  VERBOSE cout << "done" << endl;

  // -- Sanity checks on installed tokens.  Could be more extensive.

  if (Femlib::ivalue ("SVV_MN") > Geometry::nP())
    Veclib::alert (routine, "SVV_MN exceeds N_P",   ERROR);
  if (Femlib::ivalue ("SVV_MZ") > Geometry::nMode())
    Veclib::alert (routine, "SVV_MZ exceeds N_Z/2", ERROR);
}
