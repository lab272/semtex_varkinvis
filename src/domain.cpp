///////////////////////////////////////////////////////////////////////////////
// domain.cpp: implement domain class functions.
//
// A Domain amounts to a collection of Fields whose data are each to
// be computed by solution of an elliptic (Laplace, Poisson,
// Helmholtz) problem.
//
// Copyright (c) 1994+, Hugh M Blackburn
//
///////////////////////////////////////////////////////////////////////////////

#include <sem.h>


Domain::Domain (FEML*             file   ,
		const Mesh*       mesh   ,
		vector<Element*>& element,
		BCmgr*            bcmgr  ) :
// ---------------------------------------------------------------------------
// Construct a new Domain with all user Fields for use by a time
// integrator or elliptic solver (more generally, variables that have
// to satisfy a PDE and BCs).
//
//  By convention, all Fields stored in the Domain have
// single-character lower-case names.  On input, the names of the
// Fields to be created are stored in the string "field", which is
// supplied from input class bcmgr (obtained from the names in the
// FIELDS section of a session file).  At present, a Domain can
// contain *ordered* storage for one vector field (components u v,
// optionally w), an advected scalar field c, and a constraint scalar
// field p, whose gradient keeps the vector field divergence free
// (i.e. for an incompressible flow, p is the reduced pressure;
// pressure divided bu density).  We also allow just a scalar field
// (so that the Domain class can be used with an elliptic
// solver). Hence legal combinations of Fields declared in a session
// file are: c, uvp, uvwp, uvcp, uvwcp (but I/O field dumps are
// allowed other variables and orderings).  We check/enforce these
// constraints within this constructor, and subsequently the ordering
// and lengths of vectors of Field*, BoundarySys* and NumberSys* held
// by the Domain matches that in the string "field".
//
// Build assembly map information and masks to allow for solution of
// global problems (we no longer use enumerate utility to construct
// these).
// 
// No initialisation of Field MatrixSystems occurs here.
//  
// ---------------------------------------------------------------------------
  elmt (element)
{
  const char     routine[] = "Domain::Domain";
  const int_t    verbose   = Femlib::ivalue ("VERBOSE");
  const int_t    nz        = Geometry::nZProc();
  const int_t    ntot      = Geometry::nTotProc();
  const int_t    nboundary = Geometry::nBnode();
  const int_t    nel       = Geometry::nElmt();
  const int_t    next      = Geometry::nExtElmt();
  const int_t    npnp      = Geometry::nTotElmt();
  
  int_t          i, nfield, *gid;
  real_t*        alloc;
  vector<real_t> unity (Geometry::nTotElmt(), 1.0);

  strcpy ((name = new char [strlen (file -> root()) + 1]), file -> root());
  Femlib::value ("t", time = 0.0);
  step = 0;

  strcpy ((field = new char [strlen (bcmgr -> field()) + 1]), bcmgr -> field());
  nfield = strlen (field);

  if ((nfield < 1) || (nfield == 2))
    message (routine, "session must declare 1, 3, 4 or 5 fields", ERROR);

  if ((nfield == 1) && (field[0] != 'c'))
    message (routine, "session with 1 field: must be c", ERROR);

  if ((nfield == 3) && (strcmp (field, "uvp")))
    message (routine, "session with 3 fields: must be uvp", ERROR);

  if ((nfield == 4) && ((strcmp(field, "uvwp") && (strcmp(field, "uvcp")))))
    message (routine, "session with 4 fields: must be uvwp or uvcp", ERROR);

  if ((nfield == 5) && (strcmp (field, "uvwcp")))
    message (routine, "session with 5 fields: must be uvwcp", ERROR);

  if (nfield > 5)
    message (routine, "session has too many fields", ERROR);
  
  VERBOSE cout << routine << ": Domain will contain fields: " << field << endl;

  // -- Now check that constraints on (radial, azimuthal) velocity
  //    components are satisfied in the case that this is a 3D
  //    cylindrical problem: the types (Dirichlet/Neumann) of these
  //    BCs must match, off-axis, owing to the coupling required of
  //    those variables.  Blackburn & Sherwin (2004).

  if (Geometry::cylindrical() && (Geometry::nZ() > 2))
    this -> checkVBCs (file, field);

  VERBOSE cout << "  Building naive assembly mapping and inverse mass ... ";

  // -- N.B. These vectors are different lengths: the mapping vector
  //    is nel*next (i.e. nboundary) long, and the highest number it
  //    contains is (actually one less than, because we used 0-based
  //    indexing) the number of (unique) globally-numbered
  //    element-boundary nodes in the problem, _nglobal.  OTOH the
  //    inverse mass matrix has just _nglobal (<= nboundary) storage
  //    locations.

  _bmapNaive.resize (nboundary);
  mesh -> buildAssemblyMap (Geometry::nP(), &_bmapNaive[0]);
  _nglobal = _bmapNaive[Veclib::imax(_bmapNaive.size(), &_bmapNaive[0], 1)] + 1;
  _imassNaive.resize (_nglobal);
  
  Veclib::zero (_nglobal, &_imassNaive[0], 1);
  for (gid = &_bmapNaive[0], i = 0; i < nel; i++, gid += next)
    element[i] -> bndryDsSum (gid, &unity[0], &_imassNaive[0]);
  Veclib::vrecp (_nglobal, &_imassNaive[0], 1, &_imassNaive[0], 1);
    
  VERBOSE cout << "done" << endl;

  VERBOSE cout << "  Building table of all field numbering schemes ... ";

  this -> makeAssemblyMaps (file, mesh, bcmgr);

  VERBOSE cout << "done" << endl;

  // -- Build boundary system and field for each variable.
  
  VERBOSE cout << "  Building domain boundary systems and numbering... ";

  b.resize (nfield);
  n.resize (nfield);
  for (i = 0; i < nfield; i++) {
    b[i] = new BoundarySys (bcmgr, elmt,  field[i]);
    n[i] = new NumberSys   (_allMappings, field[i]);
  }

  VERBOSE cout << "done" << endl;
  
  VERBOSE cout << "  Building domain fields ... ";

  u   .resize (nfield);
  udat.resize (nfield);

  alloc = new real_t [static_cast<size_t> (nfield * ntot)];
  for (i = 0; i < nfield; i++) {
    udat[i] = alloc + i * ntot;
    u[i]    = new Field (udat[i], b[i], n[i], nz, elmt, field[i]);
  }

  VERBOSE cout << "done" << endl;
}


void Domain::checkVBCs (FEML*       file ,
			const char* field) const
// ---------------------------------------------------------------------------
// For cylindrical 3D fluids problems, the declared boundary condition types
// for velocity fields v & w must be the same for all groups, to allow
// for coupling of these fields (which uncouples the viscous substep).
//
// Check/assert this by running through each group's BCs and checking
// for tag agreement on v & w BCs.
//
// NB: this check isn't required if the problem is 2D2C or 2D3C, just 3D3C.
// ---------------------------------------------------------------------------
{
  const char routine[] = "Domain::checkVBCs";

  if (!Geometry::cylindrical() || (Geometry::nDim() < 3)) return;
  
  if (!file->seek ("BCS"))  return;
  if (!strchr (field, 'u')) return;
  if (!strchr (field, 'v') || !strchr (field, 'w')) return;

  int_t       i, j, id, nbcs;
  char        vtag, wtag, groupc, fieldc, tagc, tag[StrMax], err[StrMax];
  const int_t N (file->attribute ("BCS", "NUMBER"));

  for (i = 0; i < N; i++) {

    while ((groupc = file->stream().peek()) == '#') // -- Skip comments.
      file->stream().ignore (StrMax, '\n');

    file->stream() >> id >> groupc >> nbcs;
    vtag = wtag = '\0';

    for (j = 0; j < nbcs; j++) {

      file->stream() >> tag;
      if (strchr (tag, '<') && strchr (tag, '>') && (strlen (tag) == 3))
	tagc = tag[1];
      else {
	sprintf (err, "unrecognized BC tag format: %s", tag);
	message (routine, err, ERROR);
      }

      file->stream() >> fieldc;
      if      (fieldc == 'v') vtag = tagc;
      else if (fieldc == 'w') wtag = tagc;
      file->stream().ignore (StrMax, '\n');
    }
    
    if (!(vtag && wtag)) {
      sprintf (err, "group %c: BCs for fields 'v' & 'w' both needed", groupc);
      message (routine, err, ERROR);
    }
    if (vtag != wtag) {
      sprintf (err, "group %c, fields 'v' & 'w': BC type mismatch", groupc);
      message (routine, err, ERROR);
    }
  }
}


char Domain::axialTag (FEML* file) const
// ---------------------------------------------------------------------------
// Return the character tag corresponding to "axis" group, if that exists.
// ---------------------------------------------------------------------------
{
  if (!file->seek ("GROUPS")) return '\0';

  int_t       i;
  char        nextc, buf[StrMax];
  const int_t N (file->attribute ("GROUPS", "NUMBER"));
  
  for (i = 0; i < N; i++) {
    while ((nextc = file->stream().peek()) == '#') // -- Skip comments.
      file->stream().ignore (StrMax, '\n');
    file->stream() >> buf >> nextc >> buf;
    if (strstr (buf, "axis")) return nextc;
  }

  return '\0';
}


void Domain::checkAxialBCs (FEML* file,
			    char  atag) const
// ---------------------------------------------------------------------------
// Run through and ensure that for "axis" group, all BCs are of type <A>.
// ---------------------------------------------------------------------------
{
  if (!file->seek ("BCS")) return;

  const char routine[] = "Domain::checkAxialBCs";

  int_t       i, j, id, nbcs;
  char        groupc, fieldc, tagc, tag[StrMax], err[StrMax];
  const int_t N (file->attribute ("BCS", "NUMBER"));

  for (i = 0; i < N; i++) {

    while ((groupc = file->stream().peek()) == '#') // -- Skip comments.
      file->stream().ignore (StrMax, '\n');

    file->stream() >> id >> groupc >> nbcs;

    for (j = 0; j < nbcs; j++) {

      file->stream() >> tag;
      if (strchr (tag, '<') && strchr (tag, '>') && (strlen (tag) == 3))
	tagc = tag[1];
      else {
	sprintf (err, "unrecognized BC tag format: %s", tag);
	message (routine, err, ERROR);
      }

      file->stream() >> fieldc;

      if (groupc == atag && tagc != 'A') {
	sprintf (err, "group '%c': field '%c' needs axis BC", groupc, fieldc);
	message (routine, err, ERROR);
      }
      file->stream().ignore (StrMax, '\n');
    }
  }
}


bool Domain::multiModalBCs (FEML*       file ,
			    BCmgr*      mgr  ,
			    const char* field) const
// ---------------------------------------------------------------------------
// For cylindrical 3D problems where the axis is present in the
// solution domain, we have to cope with axial boundary conditions
// that vary with Fourier mode, see Blackburn & Sherwin JCP 197
// (2004), and these are automatically assigned by the code.  In all
// other cases, the types of BCs are taken to be the same over all
// Fourier modes considered.  This routine returns true or false (the
// latter for the simple case).  Also check that v & w BCs are able to
// be coupled (are everywhere of same type) if required for non-zero
// modes.  Various checks here were imported from (outdated)
// enumerate.cpp.
// ---------------------------------------------------------------------------
{
  if (!Geometry::cylindrical()) return false;

  // -- Checks on declarations of axis BCs (even if not eventually used).
  
  char axistag = this -> axialTag (file);
  if (axistag)   this -> checkAxialBCs (file, axistag);
  
  if (Geometry::nDim() < 3) return false;

  // -- This test is associated with cylindrical coords, 3D3C:

  this -> checkVBCs (file, field);

  // -- That's dealt with the straightforward cases.  Now we require that
  //    a BC tag of type <A> must be present and also that at least
  //    one SURFACE references a GROUP tagged with the string "axis".

  if ((file -> isStringInSection ("BCS", "<A>")) && (mgr  -> nAxis() > 0))
    return true;   // -- Will need BC sets for different Fourier modes.
  else
    return false;  // -- Will not.
}


void Domain::makeAssemblyMaps (FEML*       file,
			       const Mesh* mesh,
			       BCmgr*      mgr )
// ---------------------------------------------------------------------------
// For all the named fields in the problem, set up internal tables of
// global numbering schemes for subsequent retrieval based on supplied
// character name and Fourier mode index.  Each field name will map to
// a 3-long vector of AssemblyMap* (one each for Fourier modes 0, 1, 2+).
//
// Assembly maps are uniquely determined by
//  
// 1. the mask vector which tells us if an element-boundary node has
// an Essential/Dirichlet BC, i.e. the corresponding value is to be
// lifted out of the solution and
//
// 2. by the elliptic problem solution strategy and hence,
// numbering/assembly methodology applied to the remaining unmasked
// nodes.
//
// We assume that the same solution strategy is going to be applied to
// all the elliptic sub-problems, hence uniqueness of the mask vector
// is sufficient.
//
// Regardless of which process we are on, build AssemblyMap's for Fourier modes
// 0, 1, 2, (if they're indicated).  
// ---------------------------------------------------------------------------
{
  const int_t   strat = Femlib::ivalue ("ENUMERATION");
  int_t         i, j, mode;
  char          name;
  bool          found;
  AssemblyMap*  N;
  vector<int_t> mask (Geometry::nBnode());

  if (!this -> multiModalBCs (file, mgr, this -> field)) {
    
    // -- Numbersystems for all Fourier modes are identical.
    
    for (i = 0; i < strlen (this -> field); i++) {
      name = this -> field[i];
      mesh -> buildLiftMask (Geometry::nP(), name, 0, &mask[0]);
      for (found = false, j = 0; !found && j < _allMappings.size(); j++)
	if (found = _allMappings[j] -> willMatch (mask)) {
	  _allMappings[j] -> addTag (name, 0);
	  _allMappings[j] -> addTag (name, 1);
	  _allMappings[j] -> addTag (name, 2);
	  break;
	}
      if (!found) {
	N = new AssemblyMap (Geometry::nP(), Geometry::nElmt(),
			     strat, _bmapNaive, mask, name, 0);
	N -> addTag (name, 1);
	N -> addTag (name, 2);
	
	_allMappings.push_back (N);
      }
    }

  } else {

    // -- Number systems are different for modes 0, 1, 2+ owing to
    //    presence of axial BCs (Blackburn & Sherwin JCP 197 2004).

    for (i = 0; i < strlen(this -> field); i++) {
      name = this -> field[i];
      for (mode = 0; mode < 3; mode++) {
	mesh -> buildLiftMask (Geometry::nP(), name, mode, &mask[0]);
	for (found = false, j = 0; !found && j < _allMappings.size(); j++)
	  if (found = _allMappings[j] -> willMatch (mask)) {
	    _allMappings[j] -> addTag (name, mode);
	    break;
	  }
	if (!found) {
	  N = new AssemblyMap (Geometry::nP(), Geometry::nElmt(),
			       strat, _bmapNaive, mask, name, mode);
	  _allMappings.push_back (N);
	}
      }
    }
  }
}

#if 0

const AssemblyMap* Domain::getNsys (const char  name,
				    const int_t mode) const
// ---------------------------------------------------------------------------
// Automate retrieval of assembly mapping scheme allocated for a particular
// Field name and (global) Fourier mode index.
// ---------------------------------------------------------------------------
{
  return
    _globalNumbering.at (name)
    [clamp (mode, static_cast<int_t>(0), static_cast<int_t>(2))];
}

#endif

void Domain::report ()
// ---------------------------------------------------------------------------
// Print a run-time summary of domain & timestep information on cout.
// ---------------------------------------------------------------------------
{
  const char routine[] = "Domain::report";
  const real_t t   = time;
  const real_t dt  = Femlib:: value ("D_T");
  const real_t lz  = Femlib:: value ("TWOPI / BETA");
  const int_t  ns  = Femlib::ivalue ("N_STEP");
  const int_t  nt  = Femlib::ivalue ("N_TIME");
  const int_t  chk = Femlib::ivalue ("CHKPOINT");
  const int_t  per = Femlib::ivalue ("IO_FLD");

  cout << "-- Coordinate system       : ";
  if (Geometry::cylindrical())
    cout << "cylindrical" << endl;
  else
    cout << "Cartesian" << endl;

  cout << "   Solution fields         : " << field              << endl;

  cout << "   Number of elements      : " << Geometry::nElmt()  << endl;
  cout << "   Number of planes        : " << Geometry::nZ()     << endl;
  cout << "   Number of processors    : " << Geometry::nProc()  << endl;
  if (Geometry::nZ() > 1) cout << "   Periodic length         : " << lz<< endl;
  cout << "   Polynomial order (np-1) : " << Geometry::nP() - 1 << endl;

  cout << "   Time integration order  : " << nt                 << endl;
  cout << "   Start time              : " << t                  << endl;
  cout << "   Finish time             : " << t + ns * dt        << endl;
  cout << "   Time step               : " << dt                 << endl;
  cout << "   Number of steps         : " << ns                 << endl;
  cout << "   Dump interval (steps)   : " << per;
  if (chk) cout << " (checkpoint)";  
  cout << endl;
  cout << "   Nonlinear terms of type : ";
  switch (Femlib::ivalue ("ADVECTION")) {
  case 0: cout << "0 (Skew symmetric)"; break;
  case 1: cout << "1 (Alternating skew symmetric)"; break;
  case 2: cout << "2 (Convective/nonconservative)"; break;
  case 3: cout << "3 (Rotational-1)"; break;
  case 4: cout << "4 (Rotational-2)"; break;
  case 5: cout << "5 (Stokes/none)"; break;
  default: message (routine, "Unknown advection scheme number", ERROR);
  }
  cout << endl;
}


void Domain::restart ()
// ---------------------------------------------------------------------------
// Initialise all Field variables to zero ICs.  Then if a restart file
// "name".rst can be found, use it for input of the data it contains.
//
// Carry out forwards Fourier transformation, zero Nyquist data.
// ---------------------------------------------------------------------------
{
  int_t       i;
  const int_t nF = nField();
  char        restartfile[StrMax];
  
  for (i = 0; i < nF; i++) *u[i] = 0.0;

  ROOTONLY cout << "-- Initial condition       : ";
  ifstream file (strcat (strcpy (restartfile, name), ".rst"));

  if (file) {
    ROOTONLY {
      cout << "read from file " << restartfile;
      cout.flush();
    }
    file >> *this;
    file.close();
    transform (FORWARD);
    ROOTONLY for (i = 0; i < nF; i++) u[i] -> zeroNyquist();
  } else
    ROOTONLY cout << "set to zero";

  ROOTONLY cout << endl;
  
  Femlib::value ("t", time);
  step = 0;
}


void Domain::dump ()
// ---------------------------------------------------------------------------
// Check if a field-file write is required, carry out.
//
// Fields are inverse Fourier transformed prior to dumping in order to
// provide physical space values.
// ---------------------------------------------------------------------------
{
  const bool periodic = !(step %  Femlib::ivalue ("IO_FLD"));
  const bool initial  =   step == Femlib::ivalue ("IO_FLD");
  const bool final    =   step == Femlib::ivalue ("N_STEP");

  if (!(periodic || final)) return;
  ofstream output;

  Femlib::synchronize();

  ROOTONLY {
    const char  routine[] = "Domain::dump";
    char        dumpfl[StrMax], backup[StrMax], command[StrMax];
    const int_t verbose   = Femlib::ivalue ("VERBOSE");
    const int_t chkpoint  = Femlib::ivalue ("CHKPOINT");

    if (chkpoint) {
      if (final) {
	strcat (strcpy (dumpfl, name), ".fld");
	output.open (dumpfl, ios::out);
      } else {
	strcat (strcpy (dumpfl, name), ".chk");
	if (!initial) {
	  strcat  (strcpy (backup, name), ".chk.bak");
	  rename  (dumpfl, backup);
	}
	output.open (dumpfl, ios::out);
      }
    } else {
      strcat (strcpy (dumpfl, name), ".fld");
      if   (initial) output.open (dumpfl, ios::out);
      else           output.open (dumpfl, ios::app);
    }
    
    if (!output) message (routine, "can't open dump file", ERROR);
    if (verbose) message (routine, ": writing field dump", REMARK);
  }

  Femlib::synchronize();
  this -> transform (INVERSE);
  Femlib::synchronize();

  output << *this;

  Femlib::synchronize();
  this -> transform (FORWARD);
  Femlib::synchronize();

  ROOTONLY output.close();
}


void Domain::transform (const int_t sign)
// ---------------------------------------------------------------------------
// Fourier transform all Fields according to sign.
// ---------------------------------------------------------------------------
{
  int_t       i;
  const int_t N = this -> nField ();

  for (i = 0; i < N; i++) {
    if (sign == INVERSE) u[i] -> zeroNyquist(); 
    u[i] -> transform (sign);
    if (sign == FORWARD) u[i] -> zeroNyquist(); 
  }

}


ostream& operator << (ostream& strm,
		      Domain&  D   )
// ---------------------------------------------------------------------------
// Output all Domain field variables on ostream in prism-compatible
// form.  Binary output only.  Note that output is only done on root
// processor.
// ---------------------------------------------------------------------------
{
  int_t             i;
  const int_t       N = D.u.size();
  vector<AuxField*> field (N);

  for (i = 0; i < N; i++) field[i] = D.u[i];

  writeField (strm, D.name, D.step, D.time, field);

  return strm;
}


istream& operator >> (istream& strm,
		      Domain&  D   )
// ---------------------------------------------------------------------------
// Input all Domain field variables from prism-compatible istream.
//
// Only binary storage format is allowed.  Check if conversion to
// native format (IEEE little/big-endian) is required.
//
// Ordering of fields in file is allowed to differ from that in D.
// ---------------------------------------------------------------------------
{
  const char routine[] = "strm>>Domain";
  int_t      i, j, np, nz, nel, ntot, nfields;
  int_t      npchk,  nzchk, nelchk, verb = Femlib::ivalue ("VERBOSE");
  char       s[StrMax], f[StrMax], err[StrMax], fields[StrMax];
  bool       swap = false, found = false;

  if (strm.getline(s, StrMax).eof()) return strm;

  strm.getline(s,StrMax).getline(s,StrMax);
  
  string ss(s);
  istringstream sss (ss);
  sss >> np >> np >> nz >> nel;

  D.u[0] -> describe (f);
  sss.clear();
  sss.str (ss = f);
  sss >> npchk >> npchk >> nzchk >> nelchk;
  
  if (np  != npchk ) message (routine, "element size mismatch",       ERROR);
  if (nz  != nzchk ) message (routine, "number of z planes mismatch", ERROR);
  if (nel != nelchk) message (routine, "number of elements mismatch", ERROR);
  
  ntot = np * np * nz * nel;
  if (ntot != Geometry::nTot())
    message (routine, "declared sizes mismatch", ERROR);

  strm.getline(s,StrMax);
  sss.clear();
  sss.str  (ss = s);
  sss >> D.step;

  strm.getline(s,StrMax);
  sss.clear();
  sss.str  (ss = s);
  sss >> D.time;
  Femlib::value ("t", D.time);

  strm.getline(s,StrMax).getline(s,StrMax);
  strm.getline(s,StrMax).getline(s,StrMax);

  nfields = 0;
  while (isalpha (s[nfields])) {
    fields[nfields] = s[nfields];
    nfields++;
  }
  fields[nfields] = '\0';

  // -- Check if the restart file and the domain have same fields.
  //    (It's OK if they don't, but issue warnings.)

  ROOTONLY {
    if (nfields != strlen (D.field)) {
      sprintf (err, "  : file: %1d fields, Domain: %1d",
	       (int) nfields, (int) strlen(D.field));
      cerr << endl;
      message (routine, err, REMARK);
    }
    for (i = 0; i < nfields; i++) 
      if (!strchr (D.field, fields[i])) {
	sprintf (err, " : field %c not present in Domain (%s)",
		 fields[i], D.field);
	message (routine, err, REMARK);
      }
    for (i = 0; i < strlen (D.field); i++) 
      if (!strchr (fields, D.field[i])) {
	sprintf (err, "  : field %c not present in restart file (%s)",
		 D.field[i], fields);
	message (routine, err, REMARK);
      }
  }

  strm.getline (s, StrMax);
  Veclib::describeFormat (f);

  if (!strstr (s, "binary"))
    message (routine, "input field file not in binary format", ERROR);
  
  if (!strstr (s, "endian"))
    message (routine, "input field file in unknown binary format", WARNING);
  else {
    swap = ((strstr (s, "big") && strstr (f, "little")) ||
	    (strstr (f, "big") && strstr (s, "little")) );
    ROOTONLY {
      if (swap) cout << " (byte-swapping)";
      cout.flush();
    }
  }

  for (j = 0; j < nfields; j++) {
    for (found = false, i = 0; i < nfields; i++)
      if (fields[j] == D.field[i]) { found = true; break; }
    if (found) {    // -- Read in a matching field variable.
      strm >>  *D.u[i];
      if (swap) D.u[i] -> reverse();
    } else ROOTONLY // -- Skip over a field variable not in the domain.
      strm.seekg (ntot*sizeof(real_t), ios::cur);
  }
    
  ROOTONLY if (strm.bad())
    message (routine, "failed reading field file", ERROR);
    
  return strm;
}

