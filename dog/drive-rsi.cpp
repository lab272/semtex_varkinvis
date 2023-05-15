///////////////////////////////////////////////////////////////////////////////
// drive-rsi.cpp: compute leading eigenvalues and eigenvectors for
// flow stability analysis based on timestepping the linearised
// (Navier--Stokes) operators [2].
//
// This version (-rsi) optionally implements a real shift-invert
// solution method as described in [4].  Eigensystem solution is
// carried out using ARPACK (which has a shift-invert mode, though
// possibly not needed), and the required inner iterative linear
// system solution is carried out using Bi-CG-Stab (from Templates).
//
// The user-supplied real shift SIGMA is interpreted in the
// (time-stepping) exponentially-mapped space, i.e. a shift of unity
// (1) corresponds to a shift to the origin (0) in the untransformed
// space.  If shift-inverse is not requested, a standard/forward
// solution of A*x=lambda*x is carried out (again using ARPACK) --
// mainly intended as an optional cross-check.
//
// Despite what the ARPACK documentation says, the MXITER flag
// (iparam[2] here) is only an output value (and is the number of
// restarts).  Convergence seems solely determined by tolerance value
// (which BTW has a different meaning than for the "Barkley"
// algorithm: TOL here is usually lower than for the Barkley
// algorithm.)
//
// The base flow can be either steady or periodic in time, two or
// three component, cylindrical or Cartesian, but must be
// two-dimensional.
//
// The eigenpairs computed in the subspace are related to the Ritz
// estimates of those in the original space in a simple way: the
// eigenvalues are related to those of the original system by a simple
// transformation, and the Ritz eigenvectors are related to the
// subspace eigenvectors through a linear transformation (see [1],
// p.175).
//
// USAGE
// -----
// dog-rsi [options] session
//   session: specifies name of semtex session file.
//   options:
//   -h       ... print this message
//   -v       ... set verbose
//   -a       ... alternatively solve adjoint eigensystem
//   -S <num> ... specify (real) spectral shift, file-scope SIGMA
//   -k <num> ... set dimension of eigensystem Krylov subspace to <num>
//   -n <num> ... converge <num> eigenvalue/eigenvector pairs (n <= k)
//   -t <num> ... set eigensolution tolerance to <num>        [Default: 1e-3]
//   -m <num> ... set max number of inner iterations to <num> [Default: 200]
//   -i <num> ... set linear system solver tolerance to <num> [Default:1e-3]
//   -p       ... recompute the pressure eigenvector
// 
// FILES
// -----
// A number of semtex files are employed --
//   session:     semtex session file
//   session.num: computed automatically if not supplied
//   session.bse: base flow, containing N_SLICE field dumps
//   session.rst: (optional) restart file (else, initialise with white noise)
//
// REFERENCES
// ----------
// [1]  Y Saad (1991), "Numerical methods for large eigenvalue problems"
//      Wiley.
// [2]  LS Tuckerman & D Barkley (2000), "Bifurcation analysis for
//      timesteppers", in Numerical Methods for Bifurcation Problems,
//      ed E Doedel & LS Tuckerman, Springer. 453--466.
// [3]  D Barkley, HM Blackburn & SJ Sherwin (2008), "Direct optimal
//      growth analysis for timesteppers", IJNMF V57, 1435--1458.
// [4]  F Gomez, JM Perez, HM Blackburn & V Theofilis (2015) "On the use of
//      matrix-free shift-invert strategies for global flow instability
//      analysis", Aerosp Sci Tech V44, 69--76.
//
// Copyright (c) 2012+, Hugh Blackburn, Jose-Miguel Perez, Francisco Gomez.
///////////////////////////////////////////////////////////////////////////////

#include <stab.h>

static char             prog[] = "dog-rsi";
static char*            session;
static Domain*          domain;
static StabAnalyser*    analyst;
static vector<Element*> elmt;
static FEML*            file;
static Mesh*            mesh;
static BCmgr*           bman;

static int_t            NTOT;	        // -- Size of an eigenvector.
static int_t            MODE = FORWARD;	// -- else INVERSE.
static problem_t        TASK = PRIMAL;  // -- else ADJOINT.
static real_t           SIGMA;	        // -- Real shift.
static ofstream         RUNINFO;	// -- File for progress reports.

static void getargs (int, char**, int_t&, int_t&, real_t&, int_t&,
		     real_t&, bool&, char*&);

static int_t preprocess (const char*, bool&);
static void  EV_init    (real_t*);

// -- Routine for standard forward mode (also used in inverse mode).

static void  LNS_update (const real_t*, real_t*);

// -- Routines for shift-invert mode.

static void  INV_update (const real_t*, real_t*, const int_t&, const real_t&);
       void  matvec     (const real_t&, const real_t*, const real_t&, real_t*);
       void  ident      (real_t*, const real_t*);


int main (int    argc,
	  char** argv)
// ---------------------------------------------------------------------------
// Driver routine for stability analysis code.
//  
// The default mode is FORWARD (direct eigensystem, no shifting); if SIGMA
// is set by getargs, this is changed to INVERSE (i.e. shift-inverse).
//
// The default task is PRIMAL, for the standard eigensystem based
// on the linearised NavSto equations.  Optionally this can be set
// (again in getargs) to ADJOINT, whereby the eigensystem is based on
// the ajoint linearised NavSto equations.  This should deliver the
// same eigenvalues but different (adjoint) eigenvectors.
// ---------------------------------------------------------------------------
{
  int_t  kdim = 3, nvec = 1, nits = 200;
  real_t evtol = 1.0e-3, lstol = 1.0e-4;
  int_t  i, j, k;
  char   buf[StrMax];
  bool   restart = false, pEV = false;
  
  Femlib::initialize (&argc, &argv);

  getargs (argc, argv, kdim, nvec, evtol, nits, lstol, pEV, session);

  // -- Install lookup copies for reporting purposes (domain.cpp).

  Femlib::ivalue ("KRYLOV_KDIM", kdim);
  Femlib::ivalue ("KRYLOV_NVEC", nvec);
  Femlib::ivalue ("KRYLOV_NITS", nits);
  Femlib::value  ("KRYLOV_KTOL", evtol);

  // -- Check parameter values, intialise solver output.

  if (nvec < 1)      Veclib::alert (prog,"param error: NVEC must be > 1",      ERROR);
  if (kdim < nvec+2) Veclib::alert (prog,"param error: KDIM must be >= NVEC+2",ERROR);

  strcat (strcpy (buf, session), ".evl");
  RUNINFO.open (buf, ios::out);

  if (MODE == FORWARD)
    RUNINFO << "-- Solving forward eigensystem via timestepping." << endl;
  else
    RUNINFO << "-- Solving shift-invert eigensystem via timestepping." << endl;

  // -- Set up to run with semtex.
  
  NTOT = preprocess (session, restart);

  // -- Eigensolution by ARPACK.

  int_t       nconv;
  const int_t DONE = 99, lworkl = 3*kdim*kdim + 6*kdim;
  int_t       ido, info, iparam[11], ipntr[14], select[kdim];

  iparam [0] = 1;		// -- Shifting will be handled by ARPACK.
  iparam [1] = 0; 		// -- Not used.
  iparam [2] = nits;		// -- Input: maximum, output: number done.
  iparam [3] = 1;		// -- Blocksize, ARPACK say = 1.
  iparam [4] = 0;		// -- Output, number of converged values.
  iparam [5] = 0;		// -- Not used.
  if (MODE == FORWARD)
    iparam [6] = 1;             // -- Mode: standard, A x = lambda x.
  else
    iparam [6] = 3;             // -- Mode: shift-inverse.
  iparam [7] = 0; 		// -- For user shifts, not used here.
  iparam [8] = 0;		// -- Output, number of Op x operations.
  iparam [9] = 0;		// -- Output, not used here.
  iparam[10] = 0;		// -- Output, number of re-orthog steps.

  // -- Allocate and zero storage.

  vector<real_t> work(3*NTOT + lworkl + NTOT*kdim + NTOT +
		      2*(nvec+1) + 3*kdim + NTOT*(nvec+1), 0.0);

  real_t*        workd  = &work[0];
  real_t*        workl  = workd + 3*NTOT;
  real_t*        v      = workl + lworkl;
  real_t*        resid  = v + NTOT*kdim;
  real_t*        dr     = resid + NTOT;
  real_t*        di     = dr + nvec + 1;
  real_t*        workev = di + nvec + 1;
  real_t*        z      = workev + 3*kdim;

  // -- Either read in a restart, or set random IC. 

  EV_init (resid);

  // -- Set up for ARPACK reverse communication.

  F77NAME(dnaupd) (ido=0, "I", NTOT, "LM", nvec, evtol, resid, kdim, 
		   v, NTOT, iparam, ipntr, workd, workl, lworkl, info=1);

  if (info != 0)
    Veclib::alert (prog, "ARPACK dnaupd initialisation error", ERROR);

  // -- IRAM iteration.
  
  i = 0;
  while (ido != DONE) {

    RUNINFO << "Outer/ARPACK iteration: " << ++i << endl;
    
    if (MODE == FORWARD)
      LNS_update (workd+ipntr[0]-1, workd+ipntr[1]-1);
    else
      INV_update (workd+ipntr[0]-1, workd+ipntr[1]-1, nits, lstol);

    F77NAME(dnaupd) (ido, "I", NTOT, "LM", nvec, evtol, resid, kdim,
		     v, NTOT, iparam, ipntr, workd, workl, lworkl, info);
  }

  if (info < 0)
    Veclib::alert (prog, "ARPACK dnaupd iteration error",       ERROR);
  if (info > 0)
    Veclib::alert (prog, "ARPACK dnaupd exceeded max restarts", ERROR);
  
  RUNINFO << "--" << endl;
  RUNINFO << "Converged " 
	  << iparam[4] << " Ritz eigenvalue(s) within "
          << evtol     << " in " 
	  << iparam[8] << " operations, "
	  << iparam[2] << " restart(s)." << endl;

  // -- Post-process to obtain eigenvalues and Ritz eigenvectors.

  F77NAME(dneupd) (1, "A", select, dr, di, z, NTOT, SIGMA, 0.0, workev,
		   "I", NTOT, "LM", nvec, evtol, resid, kdim,
		   v, NTOT, iparam, ipntr, workd, workl, lworkl, info);

  // -- Print up eigenvalues.

  real_t       re_ev, im_ev, abs_ev, ang_ev, re_Aev, im_Aev;
  const real_t period = Femlib::value ("D_T * N_STEP");

  RUNINFO.precision(4);
  RUNINFO.setf(ios::scientific, ios::floatfield);

  RUNINFO << "EV  Magnitude   Angle       Growth      Frequency" << endl;

  nconv = iparam[4];
  for (j = 0; j < nconv; j++) {
    re_ev  = dr[j];
    im_ev  = di[j];
    abs_ev = hypot (re_ev, im_ev);
    ang_ev = atan2 (im_ev, re_ev);
    re_Aev = log (abs_ev) / period;
    im_Aev = ang_ev       / period;
    RUNINFO << setw(2)  << j
	    << setw(12) << abs_ev
	    << setw(12) << ang_ev
	    << setw(12) << re_Aev
	    << setw(12) << im_Aev
	    << endl;
  }

  // -- Print up eigenvectors.

  for (j = 0; j < nvec; j++) {
    char     msg[StrMax], nom[StrMax];
    real_t*  src = z + j * NTOT;
    ofstream file;
    for (i = 0; i < Geometry::nPert(); i++)
      for (k = 0; k < Geometry::nZ(); k++)
	domain -> u[i] -> setPlane 
	  (k, src + (i*Geometry::nZ() + k)*Geometry::planeSize());
    if (pEV) // -- Generate the pressure by running LNSE.
      switch (TASK) {
      case PRIMAL:  integrate (linAdvect , domain, bman, analyst); break;
      case ADJOINT: integrate (linAdvectT, domain, bman, analyst); break;
      }
    sprintf   (msg, ".eig.%1d", j);
    strcat    (strcpy (nom, domain -> name), msg);
    file.open (nom, ios::out); file << *domain; file.close();
  }

  RUNINFO.close();
  Femlib::finalize();
  return (EXIT_SUCCESS);
}


static void getargs (int      argc   ,
		     char**   argv   ,
		     int_t&   kdim   ,
		     int_t&   neval  ,
		     real_t&  evtol  ,
		     int_t&   maxit  ,		     
		     real_t&  lstol  ,
		     bool&    pEV    ,
		     char*&   session)
// ---------------------------------------------------------------------------
// Parse command-line arguments.
// ---------------------------------------------------------------------------
{
  char usage[] = "dog-rsi [options] session\n"
    "options:\n"
    "-h       ... print this message\n"
    "-v       ... set verbose\n"
    "-a       ... alternatively solve adjoint eigensystem\n"
    "-S <num> ... set (real) spectral shift to num\n"
    "-k <num> ... set dimension of Krylov subspace to num   [Default: 3   ]\n"
    "-n <num> ... compute num eigenvalue/eigenvector pairs  [Default: 1   ]\n"
    "-t <num> ... set eigensolution tolerance to num        [Default: 1e-3]\n"
    "-m <num> ... set max number of inner iterations to num [Default: 200 ]\n"
    "-i <num> ... set linear solver tolerance to num        [Default: 1e-4]\n"
    "-p       ... compute pressure from converged velocity eigenvector\n";

  while (--argc && **++argv == '-')
    switch (*++argv[0]) {
    case 'h':
      cout << usage;
      exit (EXIT_SUCCESS);
      break;
    case 'v':
      Femlib::ivalue ("VERBOSE", Femlib::ivalue("VERBOSE") + 1);
      break;
    case 'a':
      TASK = ADJOINT;
      break;
    case 'S':
      if (*++argv[0]) SIGMA = atof (  *argv);
      else { --argc;  SIGMA = atof (*++argv); }      
      MODE = INVERSE;
      break;      
    case 'i':
      if (*++argv[0]) lstol = atof (  *argv);
      else { --argc;  lstol = atof (*++argv); }
      break;
    case 'm':
      if (*++argv[0]) maxit = atoi (  *argv);
      else { --argc;  maxit = atoi (*++argv); }
      break;      
    case 'p':
      pEV = true;
      break;
    case 'k':
      if (*++argv[0]) kdim = atoi (  *argv);
      else { --argc;  kdim = atoi (*++argv); }
      break;
    case 'n':
      if (*++argv[0]) neval = atoi (  *argv);
      else { --argc;  neval = atoi (*++argv); }
      break;
    case 't':
      if (*++argv[0]) evtol = atof (  *argv);
      else { --argc;  evtol = atof (*++argv); }
      break;
    default:
      cerr << usage;
      exit (EXIT_FAILURE);
      break;
    }

  if   (argc != 1) Veclib::alert (prog, "no session file",   ERROR);
  else             session = *argv;

  // -- Here is a minor hack, installs TASK in parser so it can be
  //    reported in domain.cpp.

  Femlib::ivalue ("TASK", TASK);

  // -- While Fourier temporal interpolation is the default for base
  //    flow reconstruction, we can switch to 4-point (cubic) Lagrange
  //    interpolation by setting token LAGRANGE_INT. Here we install
  //    it in the parser table but set it to be disabled (0).

  Femlib::ivalue ("LAGRANGE_INT", 0);
}


static int_t preprocess (const char* session,
			 bool&       restart)
// ---------------------------------------------------------------------------
// Create objects needed for semtex execution, given the session file name.
//
// Return length of an eigenvector: the amount of storage required for
// a velocity component * number of components.
// ---------------------------------------------------------------------------
{
  const real_t* z;
  int_t         i, np, nel, npert;

  // -- Set default additional tokens.

  Femlib::value  ("BASE_PERIOD", 0.0);
  Femlib::value  ("T_OFFSET",    0.0);
  Femlib::ivalue ("BIG_RESIDS",  0);

  // -- Start up dealing with session file.

  file  = new FEML (session);
  mesh  = new Mesh (file);

  np    = Femlib::ivalue ("N_P");
  nel   = mesh -> nEl();
  npert = file -> attribute ("FIELDS", "NUMBER") - 1;
  Geometry::set (nel, npert);

  elmt.resize (nel);
  for (i = 0; i < nel; i++) elmt[i] = new Element (i, np, mesh);

  bman   = new BCmgr  (file, elmt);
  domain = new Domain (file, mesh, elmt, bman);

  // -- Load restart and base flow data.

  restart = domain -> restart ();
  domain -> loadBase();
  domain -> report  ();

  analyst = new StabAnalyser (domain, file);

  // -- Over-ride any CHKPOINT flag in session file.

  Femlib::value ("CHKPOINT", 1);

  return Geometry::nPert() * Geometry::planeSize() * Geometry::nZ();
}


static void EV_init (real_t* tgt)
// ---------------------------------------------------------------------------
// Load initial vector from domain velocity fields.
// ---------------------------------------------------------------------------
{
  int_t       i, k;
  const int_t ND = Geometry::nPert();
  const int_t NP = Geometry::planeSize();
  const int_t NZ = Geometry::nZ();
    
  for (i = 0; i < ND; i++)
    for (k = 0; k < NZ; k++)
      domain -> u[i] -> getPlane (k, tgt + (i*NZ + k)*NP);
}


static void LNS_update  (const real_t*   src,
			real_t*         tgt)
// ---------------------------------------------------------------------------
// Generate tgt by applying linear operator (here, a linearised
// Navier--Stokes integrator) to src.  Src and tgt could be the same
// vector.
// ---------------------------------------------------------------------------
{
  int_t       i, k;
  const int_t ND = Geometry::nPert();
  const int_t NP = Geometry::planeSize();
  const int_t NZ = Geometry::nZ();

  for (i = 0; i < ND; i++)
    for (k = 0; k < NZ; k++)
      domain -> u[i] -> setPlane (k, src + (i*NZ + k)*NP);

  switch (TASK) {
  case PRIMAL:			// -- Forward in time.
    integrate (linAdvect , domain, bman, analyst); break;

  case ADJOINT:			// -- Backward in time.
    integrate (linAdvectT, domain, bman, analyst); break;

  default:
    Veclib::alert ("EV_update", "Impossible task", ERROR); break;
  }

  for (i = 0; i < ND; i++)
    for (k = 0; k < NZ; k++)
      domain -> u[i] -> getPlane (k, tgt + (i*NZ + k)*NP);
}


static void INV_update (const real_t* src  ,
			real_t*       tgt  ,
			const int_t&  maxit,
			const real_t& lstol)
// ---------------------------------------------------------------------------
// This routine is called in updating Krylov sequence for ARPACK
// dnaupd, when running in shift-invert (INVERSE) mode.
// ---------------------------------------------------------------------------
{
#if 0
  // -- GMRES.
  
  const int_t           RESTRT = 8; // -- HMB hard-coded; GMRES restart dim.
  const int_t           wdim = NTOT * (8 + RESTRT);

  const int_t           ldw2  = RESTRT + 2;
  const int_t           wdim2 = ldw2 * 2 * (RESTRT + 2);
    
  static vector<real_t> work  (wdim );
  static vector<real_t> work2 (wdim2);
  
  Veclib::zero (wdim,  &work [0], 1);
  Veclib::zero (wdim2, &work2[0], 1);

  real_t* x    = &work[0];
  real_t* b    = x + NTOT;
  real_t* lwrk = b + NTOT;

  real_t  tol  = lstol;
  int_t   itn  = maxit;
  int_t   ier  = 0;

  cout.setf (ios::scientific, ios::floatfield); cout.precision (2);

  Veclib::copy (NTOT, src, 1, b, 1);
  Veclib::zero (NTOT, x, 1);
  
  F77NAME(gmres) (NTOT, b, x, RESTRT, lwrk, NTOT, &work2[0], ldw2,
		  itn, tol, matvec, ident, ier);

  RUNINFO << "  Inner/linear solver iterations: " << itn << endl;

  if (ier < 0) {
    RUNINFO.close();
    Veclib::alert (prog, "Error return from iterative solver", ERROR);
  }

  Veclib::copy (NTOT, x, 1, tgt, 1);

#else  // -- BICGSTAB
  const int_t           wdim = 9 * NTOT;
  static vector<real_t> work (wdim);
  
  Veclib::zero (wdim, &work[0], 1);

  real_t* x    = &work[0];
  real_t* b    = x + NTOT;
  real_t* lwrk = b + NTOT;

  real_t  tol  = lstol;
  int_t   itn  = maxit;
  int_t   ier  = 0;

  cout.setf (ios::scientific, ios::floatfield); cout.precision (2);

  Veclib::copy (NTOT, src, 1, b, 1);
  Veclib::zero (NTOT, x, 1);
      
  F77NAME (bicgstab) (NTOT, b, x, lwrk, NTOT, itn, tol, matvec, ident, ier);

  RUNINFO << "  Inner/linear solver iterations: " << itn << endl;

  if (ier < 0) {
    RUNINFO.close();
    Veclib::alert (prog, "Error return from iterative solver", ERROR);
  }

  Veclib::copy (NTOT, x, 1, tgt, 1);

#endif
}


void matvec (const real_t& alpha,
	     const real_t* x    ,
	     const real_t& beta ,
	     real_t*       y    )
// ---------------------------------------------------------------------------
// This routine is used by Templates solvers, as operator MATVEC.
//
// Its basic design intention mimics BLAS DGEMV:
//   y <-- alpha A x + beta y.
//
// But note that here it incorporates a real shift SIGMA via global variable.
//   y <-- alpha (A - SIGMA) x + beta y  
// ---------------------------------------------------------------------------
{
  static vector<real_t> work (NTOT);

  Veclib::copy (NTOT, y, 1, &work[0], 1);

  LNS_update   (x, y);                           // y <-- (NL + L)x
  Blas::axpy   (NTOT, -SIGMA, x, 1, y, 1);       // y <-- y - sigma * x
  Blas::scal   (NTOT, alpha, y, 1);              // build alpha OP x
  Blas::axpy   (NTOT, beta, &work[0], 1, y, 1);  // build alpha OP x + beta y.
}


void ident (real_t*       tgt,
	    const real_t* src)
// ---------------------------------------------------------------------------
// This is a preconditioner routine supplied to BICGSTAB, as PSOLVE.
// As we don't have a preconditioner for the problem, this is an identity. 
// ---------------------------------------------------------------------------
{
  Veclib::copy (NTOT, src, 1, tgt, 1);
}
