//////////////////////////////////////////////////////////////////////////////
// fieldforce.cpp: add in distributed body forces of various kinds.
//
// This code was originally contributed by Thomas Albrecht.
//
// This implements a somewhat generic body forcing interface. Basically,
// there are two (types of) classes involved here:
//
// 1. the 'FieldForce' class, which provides an interface to nonlinear.cpp;
//
// 2. several forcing subclasses, which actually compute a specific
//    type of forcing. Those subclasses are derived from the base
//    class 'VirtualForce'.
//
// The FieldForce instance holds the storage for the final force that
// is eventually added to the nonlinear term. To compute the final
// force, FieldForce calls all forcing subclasses and sums up their
// contribution to the final force.
//
// As of 2023, FieldForce routines operate only in physical space.
//
// Copyright (c) 2016+, Thomas Albrecht, Hugh M Blackburn
//////////////////////////////////////////////////////////////////////////////


#include <sem.h>
#include <fieldforce.h>
#include <feml.h>

static int_t NCOM; // -- Number of velocity components.


FieldForce::FieldForce (Domain* D   ,
                        FEML*   file) :
  _D (D)
// ---------------------------------------------------------------------------
// This constructor deals with <FORCING> section of FEML file.
// It maintains a list of forcing subclasses, initialized here.
//
// On initialisation, create concrete forcing type subclasses, which
// are derived from abstract class VirtualForce.  Subclasses then read
// their respective data from the session file.
// ---------------------------------------------------------------------------
{
  const char routine[] = "FieldForce::FieldForce";
  const int_t verbose  = Femlib::ivalue ("VERBOSE");
  
  NCOM = _D -> nVelCmpt();

  // -- Check for FORCE section.

  if (!file -> seek ("FORCE")) {
    VERBOSE cout << "FORCE section not found. Disabling forcing." << endl;
    _enabled = false;
    return;
  }
  _enabled = true;

  // -- Initialise data for "Canonical Steady Boussinesq" (CSB) type
  //    buoyancy.  This isn't a standard body force (rather, it
  //    pre-multiplies all nonlinear and forcing terms), but it seems
  //    cleanest to bury it within the FieldForce class.  See
  //    Blackburn, Lopez, Singh and Smits (2021).

  _CSB_enabled  = false;
  _CSB_T_REF    = 0.0;
  _CSB_BETA_T   = 0.0;
  _CSB_no_hydro = 1;

  if (_D -> hasScalar()) {
    VERBOSE cout << "  " << routine <<
      ": Canonical Steady Boussinesq initialisation" << endl;
    
    if (file -> valueFromSection (&_CSB_T_REF, "FORCE", "CSB_T_REF"))
      VERBOSE cout << "    CSB_T_REF = "        << _CSB_T_REF << endl;
    if (file -> valueFromSection (&_CSB_BETA_T, "FORCE", "CSB_BETA_T"))
      VERBOSE cout << "    CSB_BETA_T = "       << _CSB_BETA_T << endl; 
    if (file -> valueFromSection (&_CSB_no_hydro, "FORCE", "CSB_REMOVE_HYDRO"))
      VERBOSE cout << "    CSB_REMOVE_HYDRO = " << _CSB_no_hydro << endl;

    if (fabs(_CSB_BETA_T) > EPSDP) {
      _CSB_enabled = true;
      VERBOSE cout << "    ENABLED" << endl;
    }

    real_t dummy;
    if (_CSB_enabled &&
	file -> valueFromSection (&dummy, "FORCE", "BOUSSINESQ_TREF"))
      Veclib::messg (routine,
		     "cannot select regular Boussinesq and CSB", ERROR);
  }

  // -- Init body force classes.  NB: Sponge must be first in list
  //    because it initialises rather than adds to forcing.

  _classes.push_back (new SpongeForce         (_D, file));
  _classes.push_back (new CoriolisForce       (_D, file));
  _classes.push_back (new ConstForce          (_D, file));
  _classes.push_back (new WhiteNoiseForce     (_D, file));
  _classes.push_back (new SteadyForce         (_D, file));
  _classes.push_back (new ModulatedForce      (_D, file));
  _classes.push_back (new SpatioTemporalForce (_D, file));
  _classes.push_back (new DragForce           (_D, file));
  _classes.push_back (new SFDForce            (_D, file));
  _classes.push_back (new BuoyancyForce       (_D, file));
}


void FieldForce::addPhysical (AuxField*          Ni ,
                              AuxField*          buf,
                              const int_t        com,
                              vector<AuxField*>& U  )
// ---------------------------------------------------------------------------
// When called from nonlinear() to apply forcing, each subclass's
// update method is called to sum up its contribution to the force,
// which is then added to the nonlinear component term Ni.
//
// Input value com is the directional component index: 0 <==> u; 1 <==>
// v; 2 <==> w. Routine is called component-by-component.
//
// Input buf is used by each subclass as a summation buffer for
// the forcing component.
//
// U contains NADV advected fields (here, supplied in physical space).
// ---------------------------------------------------------------------------
{
  if (!_enabled) return;

  // -- Clear summation buffer.

  *buf = 0.0;

  // -- Loop through all the subclasses and make physical space
  //    additions.  Note that these get successively added into buf.
  
  vector<VirtualForce*>::iterator p;
  for (p = _classes.begin(); p != _classes.end(); p++)
    (*p) -> add (buf, com, U);

  // -- Just as for the nonlinear terms themselves, have to multiply
  //    the axial and radial components by radius if in cylindrical
  //    space.  This requires us to call addPhysical() after the
  //    nonlinear terms have been created and likewise multiplied.
  
  if (Geometry::cylindrical() && (com < 2)) buf -> mulY ();
  *Ni += *buf;
}

void FieldForce::subPhysical (AuxField*          Ni ,
                              AuxField*          buf,
                              const int_t        com,
                              vector<AuxField*>& U  )
// ---------------------------------------------------------------------------
// Like addPhysical but subtracts off (a limited number) of (hydrostatic)
// contributions.
// ---------------------------------------------------------------------------
{
  if (!_enabled) return;

  *buf = 0.0;

  vector<VirtualForce*>::iterator p;
  for (p = _classes.begin(); p != _classes.end(); p++)
    (*p) -> subtract (buf, com, U);

  if (Geometry::cylindrical() && (com < 2)) buf -> mulY ();
  *Ni += *buf;
}


void FieldForce::writeAux (vector<AuxField *>& N)
// ---------------------------------------------------------------------------
// some debug stuff.
// ---------------------------------------------------------------------------
{
  int_t      step     = _D -> step;
  const bool periodic = !(step %  Femlib::ivalue ("IO_FLD"));
  const bool initial  =   step == Femlib::ivalue ("IO_FLD");
  const bool final    =   step == Femlib::ivalue ("N_STEP");
  char       s[StrMax];

  if (!(periodic || final)) return;

  for (int i = 0; i < NCOM; i++) N[i]->setName ( 'u' + i );
  ofstream output;
  sprintf(s, "%s.f.%03i.chk", _D->name, _D->step);
  ROOTONLY output.open (s, ios::out);

  writeField(output, _D -> name, _D->step, _D->time, N);
  ROOTONLY output.close();
}


void FieldForce::canonicalSteadyBoussinesq (AuxField*          work ,
					    vector<AuxField*>& Uphys,
					    vector<AuxField*>& N    )
// ---------------------------------------------------------------------------
// See Blackburn, Lopez, Singh & Smits (2021) "On the Boussinesq
// approximation in arbitrarily accelerating frames of reference", JFM
// 924:R1.
//
//  This method re-implements what was previously in nonlinear.cpp
//  within the FieldForce class.
//  ---------------------------------------------------------------------------
{
  if (!_CSB_enabled) return;

  int_t i;

  *work  = _CSB_T_REF;
  *work -= *Uphys[NCOM];
  *work *= _CSB_BETA_T;
  *work += 1.0;
  
  for (i = 0; i < NCOM; i++) *N[i] *= *work;

  if (_CSB_no_hydro) {

    // -- Subtract out hydrostatic pressure for steady forcing.

    for (i = 0; i < NCOM; i++) this -> subPhysical (N[i], work, i, Uphys);

  }
}


// ---------------------------------------------------------------------------
// VirtualForce is the virtual base class for different types of
// forcing subclasses, providing a uniform interface.
// ---------------------------------------------------------------------------


AuxField* VirtualForce::allocAuxField (Domain *D         ,
				       char    type = '0')
{
  const int_t nTotP = Geometry::nTotProc();
  const int_t nzP   = Geometry::nZProc();

  return new AuxField (new real_t [(size_t)nTotP], nzP, D->elmt, type);
}


void VirtualForce::readSteadyFromFile (char*              fname, 
				       vector<AuxField*>& a    )
// ---------------------------------------------------------------------------
// Read steady forcing from file, store in a.
// ---------------------------------------------------------------------------
{
  const char routine[] = "VirtualForce::readSteadyFromFile";
  ifstream   input;

  input.open (fname);
  if (!input) {
    char s[StrMax];
    sprintf (s, "can't open '%s' file", fname);
    Veclib::messg (routine, s, ERROR);
  }
  readField (input, a);
  input.close ();
}


// ===========================================================================
// Concrete forcing subclasses follow.
// ===========================================================================


ConstForce::ConstForce (Domain* D   ,
			FEML*   file)
// ---------------------------------------------------------------------------
// Constructor.  A force per unit mass that is constant in both space
// in time.  (These replace old tokens FFX, FFY, FFZ.)
//
// Note that these forces are equivalant to a steady reference frame
// acceleration (or actually, its negation, since a positive body
// force per unit mass on the RHS of the NSE is the same as a negative
// frame acceleration on the LHS of the NSE).
// ---------------------------------------------------------------------------
{
  const char  routine[] = "ConstForce::ConstForce";
  const int_t verbose   = Femlib::ivalue ("VERBOSE");
  const char* tok[]     = {"CONST_X", "CONST_Y", "CONST_Z"};
  int_t       i;

  VERBOSE cout << "  " << routine << endl;

  _enabled = false;
  _D = D;

  for (i = 0; i < NCOM; i++) {
    _v[i] = 0.;	// -- default
    if (file -> valueFromSection (&_v[i], "FORCE", tok[i])) {
      _enabled = true;
      VERBOSE cout << "    " << tok[i] << " = " << _v[i] << endl;
    }
  }
}


void ConstForce::add (AuxField*          ff ,
		      const int_t        com,
		      vector<AuxField*>& U  )
// ---------------------------------------------------------------------------
// Applicator.  Add in the force everywhere in physical space (a
// cheaper alternative would be to add only to Fourier mode 0, but we
// choose to only operate in physical space).
// ---------------------------------------------------------------------------
{
  if (!_enabled) return;
  
  if (fabs (_v[com]) > EPSDP) *ff += _v[com];
}


void ConstForce::subtract (AuxField*          ff ,
			   const int_t        com,
			   vector<AuxField*>& U  )
// ---------------------------------------------------------------------------
// Applicator.  Add in the force everywhere in physical space.
// ---------------------------------------------------------------------------
{
  if (!_enabled) return;
  
  if (fabs (_v[com]) > EPSDP) *ff -= _v[com];
}


SteadyForce::SteadyForce (Domain* D   ,
			  FEML*   file)
// ---------------------------------------------------------------------------
// Constructor.  A steady, (possibly spatially varying) force per unit
// mass, computed or read in during pre-processing.  To be applied in
// physical space.
// ---------------------------------------------------------------------------
{
  const char  routine[] = "SteadyForce::SteadyForce";
  const int_t nTotP     = Geometry::nTotProc();
  const int_t nzP       = Geometry::nZProc();
  const int_t verbose   = Femlib::ivalue ("VERBOSE");
  const char* tok[]     = {"STEADY_X", "STEADY_Y", "STEADY_Z"};
  char        fname[StrMax];
  int_t       i;

  VERBOSE cout << "  " << routine << endl;
  _enabled = false;
  _D = D;

  _a.resize (NCOM);

  // -- If given a filename, read steady force from file.

  if (file -> valueFromSection (fname, "FORCE", "STEADY_FILE")) {
    _enabled = true;
    for (i = 0; i < NCOM; i++) _a[i]  = allocAuxField(D, 'u' + i);
    VERBOSE cout << "    reading file " << fname << endl;
    readSteadyFromFile (fname, _a);
  } else for (i = 0; i < NCOM; i++) {
      // -- Otherwise, try to read force components from session file.
      //    If found, allocate storage.
      char a[StrMax];
      sprintf (a, "0");	// -- default
      if (file -> valueFromSection (a, "FORCE", tok[i])) {
	_enabled = true;
	*(_a[i]  = allocAuxField(D, 'u' + i)) = a;
	VERBOSE cout << "    " << tok[i] << " = " << a << endl;
      } else _a[i] = NULL;
    }
}


void SteadyForce::add (AuxField*          ff ,
		       const int_t        com,
		       vector<AuxField*>& U  )
// ---------------------------------------------------------------------------
// Applicator.
// ---------------------------------------------------------------------------
{
  if (!_enabled) return;
  
  if (_a[com]) *ff += (*_a[com]);
}


void SteadyForce::subtract (AuxField*          ff ,
			    const int_t        com,
			    vector<AuxField*>& U  )
// ---------------------------------------------------------------------------
// Applicator.
// ---------------------------------------------------------------------------
{
  if (!_enabled) return;
  
  if (_a[com]) *ff -= (*_a[com]);
}


WhiteNoiseForce::WhiteNoiseForce (Domain* D   ,
				  FEML*   file)
// ---------------------------------------------------------------------------
// Constructor. WhiteNoiseForce adds white noise in given
// direction(s), to all locations, every _apply_step'th step
// ---------------------------------------------------------------------------
{
  const char  routine[] = "WhiteNoiseForce::WhiteNoiseForce";
  const char* tok[]     = {"WHITE_EPS_X", "WHITE_EPS_Y", "WHITE_EPS_Z"};
  const int_t verbose   = Femlib::ivalue ("VERBOSE");
  int_t       i;

  VERBOSE cout << "  " << routine << endl;
  _D = D;

  _enabled = false;

  for (i = 0; i < NCOM; i++) {
    _eps[i] = 0.;	// -- default
    if (file -> valueFromSection (&_eps[i], "FORCE", tok[i]))
      VERBOSE cout << "    " << tok[i] << " = " << _eps[i] << endl;
    _enabled = true;
  }
  _mode = -1;
  if (file -> valueFromSection (&_mode, "FORCE", "WHITE_MODE")) {
    if (_mode != -1)
      Veclib::messg (routine,
		     "WHITE_MODE deprecated, physical space only", WARNING);
    VERBOSE cout << "    " << "WHITE_MODE" << " = " << _mode << endl;
  }

  _apply_step = 1;
  if ((file -> valueFromSection (&_apply_step, "FORCE", "WHITE_APPLY_STEP")))
    VERBOSE cout <<  "  Applied every " << _apply_step << ". step." << endl;
}


void WhiteNoiseForce::add (AuxField*          ff ,
			   const int_t        com, 
			   vector<AuxField*>& U  )
// ---------------------------------------------------------------------------
// Applicator.
// ---------------------------------------------------------------------------
{
  if (!_enabled) return;
  
  if ((fabs(_eps[com]) > EPSDP) && (_D->step % _apply_step == 0))
    ff -> perturb(_eps[com], -1);
}


ModulatedForce::ModulatedForce (Domain* D   ,
				FEML*   file)
// ---------------------------------------------------------------------------
// Constructor.
// ---------------------------------------------------------------------------
{
  const char  routine[]   = "ModulatedForce::ModulatedForce";
  const int_t verbose     = Femlib::ivalue ("VERBOSE");
  const char* tok_a[]     = {"MOD_A_X",     "MOD_A_Y",     "MOD_A_Z"};
  const char* tok_alpha[] = {"MOD_ALPHA_X", "MOD_ALPHA_Y", "MOD_ALPHA_Z"};
  char        a[StrMax], fname[StrMax];
  int_t       i;

  VERBOSE cout << "  " << routine << endl;
  _enabled = false;
  _D = D;

  _a.resize (NCOM);

  bool alpha_found = false;
  for (i = 0; i < NCOM; i++) {
    // -- Try to read time-varying function alpha.
    sprintf(_alpha[i], "0");
    if (file -> valueFromSection (_alpha[i], "FORCE", tok_alpha[i])) {
      alpha_found = true;
      VERBOSE cout << "    " << tok_alpha[i] << " = " << _alpha[i] << endl;
    }
  }

  // -- Disable if no alpha found.

  if (!alpha_found) return;

  // -- Spatially-varying function a.
  //    If given a filename, read from file...

  if (file -> valueFromSection (fname, "FORCE", "MOD_A_FILE")) {
    _enabled = true;
    for (i = 0; i < NCOM; i++) _a[i]  = allocAuxField(D, 'u' + i);
    VERBOSE cout << "    reading file " << fname << endl;
    readSteadyFromFile (fname, _a);
  }
  
  // -- ... otherwise try to read force components from session file.
  //    Allocate storage only if needed.
  
  else for (i = 0; i < NCOM; i++) {
    sprintf (a, "0");	// -- defaults
    if (file -> valueFromSection (a, "FORCE", tok_a[i])) {
      _enabled = true;
      _a[i]  = allocAuxField (D, 'u' + i);
      *_a[i] = a;
      VERBOSE cout << "    " << tok_a[i] << " = " << a << endl;
    }
    else _a[i] = NULL;
  }
}


void ModulatedForce::add (AuxField*          ff,
			  const int_t        com, 
			  vector<AuxField*>& U  )
// ---------------------------------------------------------------------------
// Applicator.
// ---------------------------------------------------------------------------
{
  if (!_enabled) return;
  
  const real_t alpha = Femlib::value (_alpha[com]);

  if (_a[com] && fabs(alpha) > EPSDP) ff -> axpy (alpha, *_a[com]);
}


SpatioTemporalForce::SpatioTemporalForce (Domain* D   ,
					  FEML*   file)
// ---------------------------------------------------------------------------
// Constructor.  SpatioTemporalForce -- f = alpha(x, t).
//
// WARNING! We will evaluate alpha by parsing it on the full field EACH
// TIME STEP, which may severely degrade performance.
// ---------------------------------------------------------------------------
{
  const char  routine[]   = "SpatioTemporalForce::SpatioTemporalForce";
  const int_t verbose     = Femlib::ivalue ("VERBOSE");
  const char* tok_alpha[] = {"SPATIOTEMP_ALPHA_X",
			     "SPATIOTEMP_ALPHA_Y",
			     "SPATIOTEMP_ALPHA_Z"};
  char        a[StrMax], fname[StrMax];
  int_t       i;

  VERBOSE cout << "  " << routine << endl;
  _enabled = false;
  _D = D;

  _a.resize (NCOM);

  for (i = 0; i < NCOM; i++) {
    // -- Try to read spatio-temporally-varying function alpha.
    sprintf(_alpha[i], "0");
    if (file -> valueFromSection (_alpha[i], "FORCE", tok_alpha[i])) {
      _enabled = true;
      VERBOSE cout << "    " << tok_alpha[i] << " = " << _alpha[i] << endl;
      _a[i]  = allocAuxField(D, 'u' + i);
    }
    else _a[i] = NULL;
  }
}


void SpatioTemporalForce::add (AuxField*          ff ,
			       const int_t        com,
			       vector<AuxField*>& U  )
// ---------------------------------------------------------------------------
// Applicator.
// ---------------------------------------------------------------------------
{
  if (!_enabled) return;
  
  if (_a[com]) {
    *_a[com] = _alpha[com];
    *ff += *_a[com];
  }
}


SpongeForce::SpongeForce (Domain* D   ,
			  FEML*   file)
// ---------------------------------------------------------------------------
// Constructor.
// ---------------------------------------------------------------------------
{
  const char  routine[] = "SpongeForce::SpongeForce";
  const int_t verbose   = Femlib::ivalue ("VERBOSE");
  const char* tok_ref[] = {"SPONGE_U", "SPONGE_V", "SPONGE_W"};
  char        s[StrMax];
  int_t       i;
  
  _enabled = false;
  _D = D;
  _update = 0;

  // -- Setup and read sponge mask from session file.

  VERBOSE cout << "  " << routine << endl;
  
  _mask = allocAuxField (D, 'u');
  if (!(file -> valueFromSection (s, "FORCE", "SPONGE_M"))) return;

  _enabled = true;

  VERBOSE cout << "    SPONGE_M = " << s << endl;

  // -- Time-depended mask?

  if ((file -> valueFromSection (&_update, "FORCE", "SPONGE_UPDATE"))) {
    VERBOSE cout <<  "    SPONGE_UPDATE = " << _update << endl;
    strcpy (_mask_func, s);
  }
  else
    *_mask = s;

  // -- Read reference velocity from session file.

  _Uref.resize (3);
  for (i = 0; i < NCOM; i++) {
    sprintf (s, "0");	// -- default
    file -> valueFromSection (s, "FORCE", tok_ref[i]);
     _Uref[i] = allocAuxField (_D, 'r' + i);
    *_Uref[i] = s;
  }
}


void SpongeForce::add (AuxField*          ff ,
		       const int_t        com,
		       vector<AuxField*>& U  )
// ---------------------------------------------------------------------------
// Applicator.  Since this overwrites ff as opposed to adding to it,
// SpongeForce must be first in vector of VirtualForce, see FieldForce
// constructor.
// ---------------------------------------------------------------------------
{
  if (!_enabled) return;

  if (_update && ((_D->step % _update) == 0)) *_mask = _mask_func;

  ff -> vvmvt (*_Uref[com], *U[com], *_mask);   // ff = (uref - u) * mask
}


DragForce::DragForce (Domain* D   ,
		      FEML*   file)
// ---------------------------------------------------------------------------
// Constructor.
//
// Drag force: |F| = - m(X) * U/|U| * |U|^2 (i.e., F and U are coplanar)
// componentwise, this reads:
//   Fx = -m(X) |U| u
//   Fy = -m(X) |U| v
//   Fz = -m(X) |U| w
// m(X) is a shape function.
// ---------------------------------------------------------------------------
{
  const char  routine[] = "DragForce::DragForce";
  const int_t verbose   = Femlib::ivalue ("VERBOSE");
  char        s[StrMax];

  _enabled = false;

  // -- Setup and read mask from session file.

  VERBOSE cout << "  " << routine << endl;
  
  if (!(file -> valueFromSection (s, "FORCE", "DRAG_M")))
    return;
  _D = D;
  _mask = allocAuxField (_D, 0);
  _umag = allocAuxField (_D, 0);

  _enabled = true;
  VERBOSE cout << "  DRAG_M = " << s << endl;
  *_mask = s;
}


void DragForce::add (AuxField*          ff ,
		     const int_t        com,
		     vector<AuxField*>& U  )
// ---------------------------------------------------------------------------
// Applicator.
// ---------------------------------------------------------------------------
{
  if (!_enabled) return;

  // -- Compute velocity magnitude.
  
  if (com == 0) {
    _umag  -> mag(U);
    *_umag *= *_mask;
  }

  ff -> timesMinus (*_umag, *U[com]);
}


CoriolisForce::CoriolisForce (Domain* D   ,
			      FEML*   file)
// ---------------------------------------------------------------------------
// Constructor.
// ---------------------------------------------------------------------------
{
  const char  routine[]     = "CoriolisForce::CoriolisForce";
  const int_t verbose       = Femlib::ivalue ("VERBOSE");
  const char* tokOmega[]    = {"CORIOLIS_OMEGA_X",
			       "CORIOLIS_OMEGA_Y",
			       "CORIOLIS_OMEGA_Z"};
  const char* tokDomegaDt[] = {"CORIOLIS_DOMEGA_X_DT",
			       "CORIOLIS_DOMEGA_Y_DT",
			       "CORIOLIS_DOMEGA_Z_DT"};
  char        s[StrMax];
  int_t       i;

  VERBOSE cout << "  " << routine << endl;

  _enabled = false;
  _unsteady = 0; // int_t, bec. there's no bool version of valueFromSection()
  _D = D;
  _o.resize(3);
  _minus_2o.resize(3);

  _minus_o.resize(3);
  _DoDt.resize(3);

  file -> valueFromSection (&_unsteady, "FORCE", "CORIOLIS_UNSTEADY");

  // -- Try to read omega's components from session file.

  for (i = 0; i < 3; i++) {

    sprintf (_omega[i], "0");	// -- Defaults.
    sprintf (_DomegaDt[i], "0");
    
    if (NCOM == 3 || i == 2) {
      if (file -> valueFromSection (_omega[i], "FORCE", tokOmega[i])) {
        _enabled = true;
        VERBOSE cout << "    " << tokOmega[i] << " = " << _omega[i] << endl;
      }

      if (_unsteady && 
	  (file -> valueFromSection (_DomegaDt[i], "FORCE", tokDomegaDt[i])))
        VERBOSE cout << "    " 
		     << tokDomegaDt[i] << " = " << _DomegaDt[i] << endl;
    }
  }

  if (_enabled && !_unsteady) {
    VERBOSE cout << "    Steady Coriolis," 
      " ignoring possible CORIOLIS_DOMEGA_[XYZ]_DT."     << endl;
    VERBOSE cout << "    If desired," 
      " include centrifugal force via STEADY_[XYZ]." << endl;
   }

  if (!_unsteady)
    for (i = 0; i < 3; i++) _minus_2o[i] = -2. * Femlib::value (_omega[i]);

  _a.resize (NCOM);
  for (i = 0; i < NCOM; i++) _a[i] = allocAuxField(D, 'u' + i);
}


void CoriolisForce::add (AuxField*         ff ,
			 const int_t        com,
			 vector<AuxField*>& U  )
// ---------------------------------------------------------------------------
// Applicator.
// ---------------------------------------------------------------------------
{
  const char  routine[] = "CoriolisForce::add";
  const int_t verbose   = Femlib::ivalue ("VERBOSE");
  int_t       i;

  if (!_enabled) return;

  if (com >= NCOM) return; // -- No Coriolis force applied to the scalar field

  if (NCOM == 2 && Geometry::cylindrical())
    Veclib::messg (routine, "2C: cylindrical not implemented yet.", ERROR);
  if (_unsteady) {
    if (com == 0) 
      for (i = 0; i < 3; i++) {
	if (NCOM == 2) i = 2;
	_o[i] = Femlib::value (_omega[i]);

        // -- we'll need several negative values later on.

	_minus_o[i]  = - _o[i];
	_minus_2o[i] = - 2. * _o[i];
	_DoDt[i]     = - Femlib::value (_DomegaDt[i]);
      }

    ff -> crossXPlus (com, _DoDt);            // - DOmega/Dt x X

    // - Omega x Omega x X  (actually, we compute  + Omega x ((- Omega) x X))

    if (com == 0)
      for (i = 0; i < NCOM; i++) {
	*_a[i] = 0.;
	_a[i] -> crossXPlus (i, _minus_o);
      }
    ff -> crossProductPlus (com, _o, _a);
  }

  ff -> crossProductPlus (com, _minus_2o, U); // - 2 Omega x U
}


SFDForce::SFDForce (Domain* D   ,
		    FEML*   file)
// ---------------------------------------------------------------------------
// Constructor.
// 
// Internal storage _a plays the role of qbar in Akervik et al.'s
// description.  See Akervik et al. (2006), Phys Fluids V18.
// ---------------------------------------------------------------------------
{
  const char  routine[] = "SFDForce::SFDForce";
  const int_t verbose   = Femlib::ivalue ("VERBOSE");
  int_t       i;

  VERBOSE cout << "  " << routine << endl;

  _enabled = false;
  _SFD_DELTA = _SFD_CHI = 0.0;

  if (!(file -> valueFromSection (&_SFD_DELTA, "FORCE", "SFD_DELTA") &&
	file -> valueFromSection (&_SFD_CHI,   "FORCE", "SFD_CHI"  ))) return;

  if ((_SFD_DELTA < EPSDP) || (_SFD_CHI < EPSDP)) {
    VERBOSE Veclib::messg (routine,
		    "SFD_DELTA & SFD_CHI must both be positive to set SFD",
		     REMARK);
    return;
  }

  _D = D;
  _enabled = true;

  _a.resize (NCOM);
  for (i = 0; i < NCOM; i++)
    *(_a[i] = allocAuxField (D, 'u' + i)) = 0.0;

  VERBOSE {
    cout << "  SFD_DELTA = " << _SFD_DELTA << endl;
    cout << "  SFD_CHI   = " << _SFD_CHI   << endl;
  }
}


void SFDForce::add (AuxField*          ff ,
		    const int_t        com,
		    vector<AuxField*>& U  )
// ---------------------------------------------------------------------------
// Applicator.  Forcing -SFD_CHI * (u - qbar) is added to the
// (negative of the) nonlinear terms. Then qbar is updated using
// forwards Euler.
// ---------------------------------------------------------------------------
{
  const char   routine[] = "SFDForce::add";
  const int_t  verbose   = Femlib::ivalue ("VERBOSE");
  const real_t dt        = Femlib::value  ("D_T");
  static int_t step      = 0;  // -- Flag restart.

  if (!_enabled) return;

  // -- Restart qbar from previous velocity field.
  
  if (step < NCOM) *_a[com] = *U[com];

  // -- Subtract CHI*(u-ubar) from -nonlinear terms.

  ff -> axpy (-_SFD_CHI, * U[com]);
  ff -> axpy (+_SFD_CHI, *_a[com]);

  // -- Explicit update for ubar = qbar = _a.

  *_a[com] *= 1.0 - dt / _SFD_DELTA;
  _a[com]  -> axpy (dt / _SFD_DELTA, *U[com]);

  step++;
}


BuoyancyForce::BuoyancyForce (Domain* D   ,
			      FEML*   file)
// ---------------------------------------------------------------------------
// Constructor for Boussinesq gradient-type buoyancy terms, all
// associated with (assumed) small variations to the background
// density field. There are three possible additive contributions,
// with buoyancy force per unit mass driven by
//  
// 1. (standard) a uniform acceleration/gravity field;
// 2. (extended) gradient of kinetic energy.
// 3. (extended) centrifugal force associated with steady frame rotation;
//
// The body forces per unit mass are then
//
// rho'/rho_0 [ g + 0.5 grad (|u|^2) + 0.5 grad (|Omega x r|^2) ]
//              ~              ~                    ~     ~  
// We restrict the treatment of frame rotation centrifugal force:
//   in Cartesian,   only allow Omega_z,
//   in cylindrical, only allow Omega_x;
// in other words, we only allow rotation around the coordinate system
// symmetry axis.  Steady rotation is assumed.
//
// As far as gravity is concerned, that is restricted too:
//   in Cartesian,   only allow x and/or y components,
//   in cylindrical, only allow x component,
// in other words, the forcing obeys the symmetry of the coordinate system.
//  
// And for cylindrical coordinates we also only allow to deal with a
// gravitational component aligned with the x axis.
//
// Cylindrical:
// rho'/rho_0 [ g_x + 0.5 grad (|u|^2) + Omega_x^2 y ]
//                               ~
// Cartesian:
// rho'/rho_0 [ g   + 0.5 grad (|u|^2) + Omega_z^2 x ]
//              ~                ~                 ~
// where the x component of the last term is Omega_z^2 x
// while the y component of the last term is Omega_z^2 y
//
// Reference: Blackburn Lopez Singh & Smits JFM 2021.

// As noted above, what is provided here are the "gradient-type"
// buoyancy terms, only.  The more general treatment outlined in the
// reference, where buoyancy forces are associated with all non-local
// accelerative terms, is dealt with in nonlinear.cpp.
// ---------------------------------------------------------------------------
{
  const char  routine[] = "BuoyancyForce::BuoyancyForce";
  const int_t verbose   = Femlib::ivalue ("VERBOSE");

  real_t      gravMag = 0.0, norm = 0.0;

  VERBOSE cout << "  " << routine << endl;
  _enabled = false;
 
  _TREF = _BETAT = _omega = _gx = _gy = 0.0;
  _centrifugal = _kineticgrad = 0;

  // --Return immediately if mininum data for Boussinesq aren't present.
  
  if (!file -> valueFromSection (&_TREF,   "FORCE", "BOUSSINESQ_TREF"))
    return;
  else
    VERBOSE cout << "    BOUSSINESQ_TREF = "    << _TREF << endl;

  if (file -> valueFromSection (&_BETAT,  "FORCE", "BOUSSINESQ_BETAT"))
    VERBOSE cout << "    BOUSSINESQ_BETAT = "   << _BETAT << endl;
  if (file -> valueFromSection (&gravMag, "FORCE", "BOUSSINESQ_GRAVITY"))
    VERBOSE cout << "    BOUSSINESQ_GRAVITY = " << gravMag << endl;
  if ((gravMag < -EPSDP) || (fabs(_BETAT) < -EPSDP))
    Veclib::messg (routine,
		   "gravity and/or expansion coeff. magnitudes < 0", ERROR);

  if (Geometry::cylindrical()) {
    if (file -> valueFromSection (&_gx, "FORCE", "BOUSSINESQ_GX"))
      VERBOSE cout << "    " << "BOUSSINESQ_GX = " << _gx << endl;
  } else {
    if (file -> valueFromSection (&_gx, "FORCE", "BOUSSINESQ_GX"))
      VERBOSE cout << "    " << "BOUSSINESQ_GX = " << _gx << endl;
    if (file -> valueFromSection (&_gy, "FORCE", "BOUSSINESQ_GY"))
      VERBOSE cout << "    " << "BOUSSINESQ_GY = " << _gy << endl;
  }
  
  if ((norm = sqrt(_gx*_gx + _gy*_gy)) < EPSDP)
    Veclib::messg (routine, "no active gravity vector component", REMARK);
  else {
    _gx *= gravMag/norm;
    _gy *= gravMag/norm;
  }

  // -- Check for extension 2, kinetic energy gradient buoyancy.
    
  if (file -> valueFromSection (&_kineticgrad, "FORCE", "BOUSSINESQ_KINETIC")) {
    if (_kineticgrad == 1) {
      VERBOSE cout << "    Boussinesq buoyancy will include grad(KE)" << endl;
    } else {
      _kineticgrad = 0;	// -- Only allowed values are 0 and 1.
    }
  }

  // -- Check for extension 3, centrifugal buoyancy.

  if (file -> valueFromSection (&_centrifugal,"FORCE","BOUSSINESQ_CENTRIF")) {
    if (_centrifugal == 1) {
      VERBOSE cout << "    Boussinesq centrifugal buoyancy enabled" << endl;
      if (Geometry::cylindrical()) {
	if (!file -> valueFromSection (&_omega, "FORCE", "CORIOLIS_OMEGA_X")) {
	  Veclib::messg (routine,
			 "could not find (expected) CORIOLIS_OMEGA_X",ERROR);
	}
      } else {			// -- Cartesian.
	if (!file -> valueFromSection (&_omega, "FORCE", "CORIOLIS_OMEGA_Z")) {
	  Veclib::messg (routine,
			 "could not find (expected) CORIOLIS_OMEGA_Z", ERROR);
	}
      }
    } else {
      _centrifugal = 0;	// -- Only allowed values are 0 and 1.
    }
  }
  
  // -- If we got this far, everything should be OK.
  
  _enabled = true;
  _D = D;
 
  if (_centrifugal || _kineticgrad) {
    _a.resize (3);
    _a[0] = allocAuxField (_D);   // -- Storage for relative density variation.
    _a[1] = allocAuxField (_D);   // -- Storage for quadratic scalar field.  
    _a[2] = allocAuxField (_D);   // -- Workspace.
  } else {
    _a.resize (1);
    _a[0] = allocAuxField (_D);   // -- Storage for relative density variation.
  }
}


void BuoyancyForce::add (AuxField*          ff ,
			 const int_t        com, // Vel compt index.
			 vector<AuxField*>& U  )
// ---------------------------------------------------------------------------
// Applicator for Boussinesq buoyancy.  We exploit the fact that while
// the velocity components are dealt with in order, the physical space
// velocity and scalar data in U remain the same for each call in a
// timestep.
//
// Cylindrical: (rho'/rho_0) [ g_x + 0.5 grad (|u|^2) + Omega_x^2 y ]
//
// Cartesian:   (rho'/rho_0) [ g   + 0.5 grad (|u|^2) + Omega_z^2 x ]
// ---------------------------------------------------------------------------
{
  if (!_enabled) return;

    if (com == 0) {		// -- First time through on this timestep.
    if (_kineticgrad)
      (_a[1] -> innerProduct (U, U, NCOM)) *= 0.5;
    else if (_centrifugal)
      *_a[1] = 0.0;
    *_a[0]  = _TREF;
    *_a[0] -= *U[NCOM];  // -- U[NCOM] contains temperature.
    *_a[0] *= _BETAT;    // -- _a[0] = rho'/rho_0, relative density variation.
  }

  switch (com) {
    
  case 0:
    if (_centrifugal || _kineticgrad) *_a[2] = 0.0;
    if ((!(Geometry::cylindrical())) && _centrifugal) {
      (*_a[2] = _omega / sqrt(2.0)) . mulX();
      *_a[2] *= *_a[2];
    }
    if (_kineticgrad)
      *_a[2] += *_a[1];
    if (_centrifugal || _kineticgrad) {
      _a[2] -> gradient(0);
      ff -> timesPlus (*_a[0], *_a[2]);
    }
    if (fabs (_gx) > EPSDP) ff -> axpy (_gx, *_a[0]);
    break;
    
  case 1:
    if (_centrifugal || _kineticgrad) *_a[2] = 0.0;
    if (_centrifugal) {
      (*_a[2] = _omega / sqrt(2.0)) . mulY();
      *_a[2] *= *_a[2];
    }
    if (_kineticgrad)
      *_a[2] += *_a[1];
    if (_centrifugal || _kineticgrad) {
      _a[2] -> gradient(1);
      ff -> timesPlus (*_a[0], *_a[2]);
    }
    if (fabs (_gy) > EPSDP) ff -> axpy (_gy, *_a[0]);

    break;

  case 2:
    if (_kineticgrad && Geometry::nDim() > 2) {
      (*_a[2] = *_a[1]) . transform(FORWARD).gradient(2).transform(INVERSE);
      ff -> timesPlus (*_a[0], *_a[2]);
    } 
    break;
  }
}
