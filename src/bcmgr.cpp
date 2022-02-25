///////////////////////////////////////////////////////////////////////////////
// bcmgr.cpp: class functions for managing boundary conditions.
//
// Copyright (c) 1994+, Hugh M Blackburn
//
// SYNOPSIS
// --------
// BCmgr manufactures and stores instances of classes derived from the
// Condition virtual base class, which are subsequently used by the
// Boundary class to apply boundary conditions.
//
// INPUT
// -----
// The following sections are optionally present in FEML file: they
// are needed if the <SURFACES> section used any boundary conditions,
// but it is possible they could be absent if all boundaries are
// periodic.  The sections below are given as examples.
//
// <GROUPS NUMBER=6>
// #	tag	name	descriptor
// #    Note that names "open" and "axis" are reserved and must 
// #    be associated with the BC character tags 'O' and 'A' respectively.
// #    The one-character names should be unique.
// 	1	o	open
//      2       a       axis
// 	3	v	velocity
// 	4	w	wall
//      5       h       radiation
//      6       e       exit 
// </GROUPS>
// 
// <BCS NUMBER=6>
// #	tag	group	number, followed by BCs.
//
// #    BCs for the "open" group are internally computed Robin BCs according
// #    to field name, as outlined in ref [3].
// 	1	o	4
// 			<O>	u                       </O>
// 			<O>	v                       </O>
// 			<O>	w                       </O>
// 			<O>	p                       </O>
// #    BCs for the "axis" group are internally assigned according to
//      field name and Fourier mode, see section below and ref [2].
//      2       a       4
//                      <A>     u                       </A>
//                      <A>     v                       </A>
//                      <A>     w                       </A>
//                      <A>     p                       </A>
// #    For all other groups, the user can select arbitrary combinations.
// 	3	v	4
// 			<D>	u = 1.0-4.0*(y-0.5)^2.0 </D>
// 			<D>	v = 0.0                 </D>
// 			<D>	w = 0.0                 </D>
// 			<H>	p                       </H>
// 	4	w	4
// 			<D>	u = 0.0                 </D>
// 			<D>	v = 0.0                 </D>
// 			<D>	w = 0.0                 </D>
// 			<H>	p                       </H>
//      5       h       4
//                      <M>     u = 1.0,2.0             </M>
//                      <M>     v = 1.0,0.5             </M>
//                      <M>     w = 0.5,2.0             </M>
//                      <H>     p                       </H>
//      6       e       4
//                      <N>     u = 0.                  </N>
//                      <N>     v = 0.                  </N>
//                      <N>     w = 0.                  </N>
//                      <D>     p = 0.                  </D>
// </BCS>
//
// <SURFACES NUMBER=6>
// # Surfaces associate a particular element face with a group of BCs (<B>)
// # OR nominate the face as periodic with another (<P>).
// #       tag     elmt    face    type
//         1       1       1       <P>     3       3       </P>
//         2       2       1       <P>     4       3       </P>
//         3       2       2       <B>     o       </B>
//         4       4       2       <B>     c       </B>
//         5       3       4       <B>     v       </B>
//         6       1       4       <B>     v       </B>
// </SURFACES>
//
// GROUPS
// ------

// Groups provide a way of associating sets of boundary conditions for
// different field variables via a single character tag (as used in
// the SURFACES section), and also provide the capability of
// associating a string description with a character tag, so that
// internal routines can gain additional information about the tag if
// required.  Description strings cannot contain whitespace.  The
// character tag for each group must be unique, but is arbitrary.
//
// Three special group strings are "wall", "axis" and "open".
// E.g. "wall" denotes that all velocity components are zero: this
// means that one Condition can gain information about behaviour of
// another Condition which shares the character tag (i.e. that it is
// also a zero-valued essential boundary for velocity components).
// For the "axis" and "open" groups, note that no BC value information
// is required (values are either set or computed internally) and that
// every BC in these groups must be of the same type.  Other than
// these three special instances, the string used for any group is an
// arbitrary choice of the analyst, and the mixture of types of BCs
// for each field variable is also arbitrary.
//
// The string associated with each group also allows user routines to
// access boundary value storage areas by using the "addForGroup" and
// "zeroForGroup" routines, which use this string as a specifier.
// These capabilities are rarely used, but are provided to give direct
// access to encapsulated storage if required.
//
// BCS
// ---
// BCs associate specific BC types and values with a particular group.
// Available types are:
//   <D> Dirichlet/essential;
//   <N> Neumann/natural;
//   <M> Mixed/Robin;
//   <H> Natural pressure BC (no value specified, since it gets computed);
//   <A> Axis BCs for cylindrical coords.  Also, must belong to "axis" group;
//   <O> Computed open BCs, see [3]. Must belong to the "open" group;
//   <I> Computed open BC for u&v, mixed for w, Dirichlet for c. "inlet" group.
//
// The character tags for variables as shown match those used
// internally as Field names, so that the order in which the BCs are
// supplied for each group is arbitrary, however the number must match
// or exceed the number of Field variables in the problem (this means
// that 3D BC specifications can also be used for 2D or 1D problems,
// but not vice-versa).
//
// For boundary value types where a function is supplied after the
// "=", it is necessary that no spaces appear within the function
// string.
//
// VARIABLE NAMES
// --------------
// The (one character) names of field variables are significant, and have
// the following reserved meanings:
// 
// u:  First velocity component.            (Cylindrical: axial     velocity.)
// v:  Second velocity component.           (Cylindrical: radial    velocity.)
// w:  Third velocity component.            (Cylindrical: azimuthal velocity.)
// p:  Pressure divided by density.
// c:  Scalar for transport or elliptic problems.
//
// In semtex, the velocity component w is always associated with the
// Fourier direction.
//
// CYLINDRICAL COORDINATE SYSTEM AND BCS, see ref [2]
// --------------------------------------------------
//
// Type <A> BCs are used for cylindrical coordinate systems when the
// edge of an element touches the axis (N.B. y (r) negative is
// illegal).  The BCs that are supplied for boundary type <A> depend
// on (a) the field in question (u, v, [w], p, [c]) and (b) the index
// of the Fourier mode.
//
// Fields u (axial velocity), p, & c are have the same treatment on
// the axis --- they have "scalar" type BCs there, due to the need for
// the fields to be single-valued at the axis.
//
// Fields v & w (radial and azimuthal velocities) have their BCs set
// after the fields have been coupled to give v~ and w~ (which
// *decouples* the Helmholtz problems for diffusion).  V & w must have
// a set phase relationship in the first Fourier mode resulting in a
// zero essential BC in the first Fourier mode for v~ and a zero
// natural BC in the first Fourier mode for w~.  In order to be able
// to set BCs for the coupled fields, v & w must have matching BC
// kinds on all remaining boundaries.  This is enforced in the
// enumeration phase.
//
// Summary:
//
//         +------------+------------------------------------------+
//         |  (Coupled) |       Axis BC for Fourier Mode Index     |
//         |   Variable |         0             1           2...   |
//         +------------+------------------------------------------+
//         |   u, p, c  |   du/dr = 0       u   = 0      u  = 0    |
//         |      v~    |     v~  = 0       v~  = 0      v~ = 0    |
//         |      w~    |     w~  = 0    dw~/dr = 0      w~ = 0    |
//         +------------+------------------------------------------+
//
// In order to deal with this modal dependence of BCs, three levels of
// Boundary pointers are maintained for cylindrical Fields,
// corresponding to Fourier modes 0, 1, 2 (and higher), even if there
// are no axial BCs to be applied.  These BC pointers are carried by
// all processors.
//
// The names of the numbering schemes that will be used for cylindrical
// 3D problems with axial BCs are (note case sensitivity):
//         +------------+------------------------------------------+
//         | (Coupled~) | Numbering scheme for Fourier mode index  |
//         |  Variable  |         0             1           2...   |
//         +------------+------------------------------------------+
//         |      u     |         u             U           U      |
//         |      v~    |         v             v           v      |
//         |      w~    |         w             W           w      |
//         |      p     |         p             P           P      |
//         |      c     |         c             C           C      |
//         +------------+------------------------------------------+
//
// REFERENCES
// ----------
// [1] Karniadakis, Israeli & Orszag, JCP 97:414-443 (1991)
// [2] Blackburn & Sherwin, JCP 197:759-778 (2004)
// [3] Dong, JCP 302:300-328 (2015)
// [4] Dong, Karniadakis & Chryssostomidis, JCP 261:83-105 (2014)
//
// N.B. Typo in eq. (37) of [3], confirmed by author: the term n x
// \omega should be n . \nabla x \omega.
//
///////////////////////////////////////////////////////////////////////////////

#include <sem.h>


BCmgr::BCmgr (FEML*             file,
	      vector<Element*>& elmt) :
// ---------------------------------------------------------------------------
// This constructor deals with <FIELDS>, <GROUPS>, <BCS> and
// <SURFACES> sections of FEML file, and loads internal tables for
// later use by BCmgr::getCondition.
//
// In addition, it reads in prebuilt numbering schemes from file
// session.num for later retrieval.
// ---------------------------------------------------------------------------
  _axis   (false),
  _open   (false),
  _toggle (false)     // -- Can only ever become true if _open is too.
{
  const char  routine[] = "BCmgr::BCmgr";
  char        buf[StrMax], err[StrMax], tag[StrMax], gat[StrMax];
  char        groupc, fieldc, testc, tagc;
  char*       trailer;
  const char* grpdesc;
  int_t       verbose = Femlib::ivalue ("VERBOSE");
  int_t       i, j, N, id, nbcs;
  Condition*  C;
  CondRecd*   R;

  // -- Load FIELDS.

  if (file -> seek ("FIELDS")) {
    file -> stream().ignore (StrMax, '\n');
    while (file -> stream().peek() == '#') // -- Skip comments.
      file -> stream().ignore (StrMax, '\n');
    i = 0;
    do
      file -> stream() >> buf[i++];
    while (buf[i - 1] != '<' && i < StrMax);
    if (buf[--i] == '<') {
      buf[i] = '\0';
      file -> stream() >> tag;
      if (!(strstr (tag,    "/FIELDS")))
	   message (routine, "FIELDS section not closed", ERROR);
    } else message (routine, "FIELDS section not closed", ERROR);
  } else   message (routine, "FIELDS section not found",  ERROR);

  strcpy ((_fields = new char [strlen (buf) + 1]), buf);

#if 0
  VERBOSE cout << routine << ": Installing numbering systems ... ";

  this -> buildnum (file -> root(), elmt);

  VERBOSE cout << "done" << endl;
#endif
  
  if (!file -> seek ("GROUPS")) {
    if (verbose)
      message (routine, "no GROUPS, assuming no boundary conditions", WARNING);
    return;
  }

  // -- Load GROUPS.

  N = file -> attribute ("GROUPS", "NUMBER");

  VERBOSE cout << "  Searching for " << N << " GROUPS ... ";

  _group   .resize (N);
  _descript.resize (N);
  
  for (i = 0; i < N; i++) {
    while (file -> stream().peek() == '#') // -- Skip comments.
      file -> stream().ignore (StrMax, '\n');
    file -> stream() >> id >> groupc >> buf;
    _group[i] = groupc;
    strcpy ((_descript[i] = new char [strlen (buf) + 1]), buf);
    if (!strcmp (_descript[i], "axis")) _axis = true;
    if (!strcmp (_descript[i], "open")) _open = true;
  }
  
  VERBOSE cout << "done" << endl;

  // -- Load BCS.

  N = file -> attribute ("BCS", "NUMBER");

  VERBOSE cout << "  Searching for " << N << " BCS ... ";

  for (i = 0; i < N; i++) {

    while (file -> stream().peek() == '#') // -- Skip comments.
      file -> stream().ignore (StrMax, '\n');

    file -> stream() >> id >> groupc >> nbcs;

    grpdesc = groupInfo (groupc);	// -- Ensure matching group exists.

    for (j = 0; j < nbcs; j++) {

      // -- Open tag.

      file -> stream() >> tag;
      if (strchr (tag, '<') && strchr (tag, '>') && (strlen (tag) == 3))
	tagc = tag[1];
      else {
	sprintf (err, "unrecognized BC tag: %s", tag);
	message (routine, err, ERROR);
      }

      // -- Decide if this is a value or function BC.

      file -> stream() >> fieldc;
      file -> stream() >> testc;
      if (testc == '=') { file -> stream() >> buf; strtod (buf, &trailer); }
      if (testc == '<') file -> stream().putback (testc);

      // -- Create appropriate derived Condition structure.
      //    We make at least one per BC group+field variable for 
      //    later retrieval.

      // -- The two primary kinds of BC are Dirichlet/Essential and
      //    Neumann/Natural.  Mixed is a linear combination of these.
      //    Those are all BCs defined by a declared value or function
      //    string that can be evaluated.

      // -- The Axis type amounts to either homogenous Dirichlet or
      //    Neumann depending on variable name and Fourier mode index,
      //    as explained in [2].

      // -- The remaining BC types are computed internally but again
      //    they come out to be either Dirichlet or Neumann types.

      // -- High-order pressure (H) is an internally-computed Neumann BC,
      //    explained in [1].

      // -- Open boundaries (O) comprise a set of computed Mixed/Robin
      //    conditions for velocity and pressure, outlined in [3],
      //    supplemented here by homogeneous Neumann (for scalar). Since
      //    they appear together, the code allocates all of them and
      //    decides which one to apply for a given field variable.

      switch (tagc) {

      case 'D': case 'E':	// -- Dirichlet/Essential BC.
	if (testc != '=') {
	  sprintf (err, "expected an '=' in setting field '%c' BC", fieldc);
	  message (routine, err, ERROR);
	}
	if   (*trailer != 0) C = new EssentialFunction (buf);
	else                 C = new Essential         (buf);
	break;

      case 'N':			// -- Neumann/Natural BC.
	if (testc != '=') {
	  sprintf (err, "expected an '=' in setting field '%c' BC", fieldc);
	  message (routine, err, ERROR);
	}
	if   (*trailer != 0) C = new NaturalFunction (buf);
	else                 C = new Natural         (buf);
	break;

      case 'M':			// -- Mixed BC.
	if (testc != '=') {
	  sprintf (err, "expected an '=' in setting field '%c' BC", fieldc);
	  message (routine, err, ERROR);
	}
	if (!(strchr (buf, ';') || strchr (buf, ','))) {
	  sprintf (buf,"can't find multiplier and reference value in: %s",buf);
	  message (routine, buf, ERROR);
	}
	C = new Mixed (buf);
	break;

      case 'A':			// -- Axis BC.

	if (Geometry::system() != Geometry::Cylindrical)
	  message (routine, "axis BCs disallowed in Cartesian coords", ERROR);

	// -- Create two kinds to be retrieved later.
	//    Ensure that the group name is "axis" to aid retrieval.

	if (!strstr (groupInfo (groupc), "axis"))
	  message (routine, "type 'A' BC must belong to group \"axis\"",ERROR);

	strcpy (buf, "0.0");

	C = new Natural (buf);
	
	_cond.insert (_cond.end(), R = new CondRecd);
	R -> grp    = groupc;
	R -> fld    = fieldc;
	R -> bcn    = C;
	strcpy ((R -> value = new char [strlen (buf) + 1]), buf);
	
	C = new Essential (buf);
	break;

      case 'H':		// -- "High Order" computed natural pressure BC.
	if (fieldc != 'p') {
	  sprintf (err, "expected name 'p' with HOPBC, read '%c'", fieldc);
	  message (routine, err, ERROR);
	}
	C = new NaturalCBCp (this);
	break;

      case 'O':			// -- Open BC.

	if (!strstr (groupInfo (groupc), "open"))
	  message(routine,"type 'O' BC must belong to group \"open\"",ERROR);

	if      (fieldc == 'u') C = new MixedCBCu (this);
	else if (fieldc == 'v') C = new MixedCBCv (this);
	else if (fieldc == 'w') C = new MixedCBCw (this);
        else if (fieldc == 'p') C = new MixedCBCp (this);
	else if (fieldc == 'c') {
	  strcpy (buf, "0.0"); C = new Natural (buf);
	}
	else {
	  sprintf (err,"field name '%c'for open BC not in 'uvwpc'", fieldc);
	  message (routine, err, ERROR);
	}
	break;
	
#if 1
      case 'I':			// -- Inlet BC, set value for scalar
	                        //    and with dw/dn + K w = 0.
                               
	if (!strstr (groupInfo (groupc), "inlet"))
	  message(routine,"type 'I' BC must belong to group \"inlet\"",ERROR);

	if      (fieldc == 'u') C = new MixedCBCu   (this);
	else if (fieldc == 'v') C = new MixedCBCv   (this);
	else if (fieldc == 'w') C = new MixedCBCwIn (this);
	
//	else if (fieldc == 'w') {
//	  C = new Essential ("W_INLET");
//	  strcpy (buf, "1e12;W_INLET"); C = new  Mixed (buf);	  
//	}      
        else if (fieldc == 'p') C = new MixedCBCp (this);
       	else if (fieldc == 'c') {
	  // -- Allocate a set scalar value based on a magic token.
	  strcpy (buf, "T_INLET");
	  //strcpy (buf, "LAPSE*((x-X_CUT)+heav(x-X_CUT)*(X_CUT-x))");
	  C = new EssentialFunction (buf);
	}
	else {
	  sprintf (err,"field name '%c' for openS BC not in 'uvwpc'", fieldc);
	  message (routine, err, ERROR);
	}
	break;
#endif	

      default:
	sprintf (err, "unrecognized BC identifier: %c", tagc);
	message (routine, err, ERROR);
	break;
      }

      // -- Close tag.

      file -> stream() >> gat;

      if (strlen (gat) != 4
	  || gat[0] != '<'
	  || gat[1] != '/' 
	  || gat[2] != tagc
	  || gat[3] != '>') {
	sprintf (err, "close tag %s didn't match open tag %s", gat, tag);
	message (routine, err, ERROR);
      }

      // -- Install new Condition record in internal list.

      _cond.insert (_cond.end(), R = new CondRecd);
      R -> grp = groupc;
      R -> fld = fieldc;
      R -> bcn = C;
      strcpy ((R -> value = new char [strlen (buf) + 1]), buf);
    }
  }

  VERBOSE cout << "done" << endl;

  VERBOSE cout << "  Building internal list of BC edges ... ";

  this -> buildsurf (file, elmt);

  VERBOSE cout << "done" << endl;
}


const Condition* BCmgr::getCondition (const char  group,
				      const char  field,
				      const int_t mode )
// ---------------------------------------------------------------------------
// Search internal list of Conditions for appropriate one, such that
// the stored internal variables match input group & field names.
// Return pointer to the matching Condition.
//
// Have to get into more extended search to locate axis BCs.  Input
// variable 'mode' is only used for cylindrical coordinates: it is the
// Fourier mode number.  See file header above and Ref [2] for more
// information.
// ---------------------------------------------------------------------------
{ 
  const char                  routine[] = "BCmgr::getCondition";
  char                        buf[StrMax], err[StrMax], currgrp, currfld;
  CondRecd*                   C;
  vector<CondRecd*>::iterator c;
  
  for (c = _cond.begin(); c != _cond.end(); c++) {
    C       = *c;
    currgrp = C -> grp;
    currfld = C -> fld;

    if (currgrp == group && currfld == field) {

      C -> bcn -> describe (buf);

      if (strstr (buf, "axis")) {

	switch (field) {
	case 'u': case 'p': case 'c':
	  if ((mode >  0 && strstr (buf, "essential")) ||
	      (mode == 0 && strstr (buf, "natural"    ))) return C -> bcn;
	  break;
	case 'v':
	  if               (strstr (buf, "essential"))    return C -> bcn;
	  break;
	case 'w':
	  if ((mode != 1 && strstr (buf, "essential")) ||
	      (mode == 1 && strstr (buf, "natural"    ))) return C -> bcn;
	  break;
	default:
	  sprintf (err, "unrecognised field '%c' on axis", field);
	  message (routine, err, ERROR);
	  break;
	} 

      } else 			// -- Default.
	return C -> bcn;
    }
  }

  sprintf (err, "can't find record for group '%c', field '%c'", group, field);
  message (routine, err, ERROR);
  return 0;
}


const char* BCmgr::groupInfo (const char name) const
// ---------------------------------------------------------------------------
// Given a group name, return pointer to string descriptor.
//
// Force a halt on failure.  
// ---------------------------------------------------------------------------
{
  const char  routine[] = "BCmgr::groupInfo";
  const int_t N = _group.size();
  char        err[StrMax];
  int_t       i;

  for (i = 0; i < N; i++) if (name == _group[i]) return _descript[i];

  sprintf (err, "unknown group: %c", name);
  message (routine, err, ERROR);
  return 0;
}

#if 0

NumberSys* BCmgr::getNumberSys (const char  name,
				const int_t mode)
// ---------------------------------------------------------------------------
// Return the NumberSys corresponding to the input field name and
// Fourier mode.  Mode numbers begin at zero, and are not expressed
// modulo number of modes per process.
// ---------------------------------------------------------------------------
{
  const char  routine[] = "BCmgr::getNumberSys";
  const int_t nsys  = _numsys.size();
  const int_t cmode = clamp (mode,static_cast<int_t>(0),static_cast<int_t>(2));
  char        err[StrMax], selectname = name;
  NumberSys*  N = 0;

  // -- Switch the name for selected cylindrical coordinate modes.

  if (Geometry::cylindrical() && Geometry::nDim() == 3 && _axis)
    if      (name == 'u' && (cmode == 1 || cmode == 2)) selectname = 'U';
    else if (name == 'w' && (cmode == 1))               selectname = 'W';
    else if (name == 'p' && (cmode == 1 || cmode == 2)) selectname = 'P';
    else if (name == 'c' && (cmode == 1 || cmode == 2)) selectname = 'C';

  // -- Now search.
  
  for (int_t i = 0; i < nsys; i++)
    if (strchr (_numsys[i] -> fields(), selectname)) { N = _numsys[i]; break; }

  if (!N) {
    sprintf (err, "can't find scheme for field %c, mode %1d", name, mode);
    message (routine, err, ERROR);
  }
  
  return N;
}


void BCmgr::buildnum (const char*       session,
		      vector<Element*>& elmt   )
// ---------------------------------------------------------------------------
// Private member function.
//
// Retrieve numbering schemes (btog and bmsk values) from file
// "session.num": this is created by running the "enumerate" utility
// on root processor, if session.num (made by the enumerate) does not
// already exist.
//
// The names of fields and their numbering schemes are significant.
// The convention employed is that the fields have lower-case
// single-character names.  Numbering schemes have the same names,
// *except* in the case of cylindrical coordinate systems where the
// domain includes the symmetry axis.  See top of this file for
// mode-related significance for upper-cased names of numbering
// schemes.
//
// After the numbering schemes have been set up, create the
// corresponding inverse mass matrices.
// ---------------------------------------------------------------------------
{
  const char     routine[] = "BCmgr::buildnum";
  const int_t    np   = Geometry::nP();
  const int_t    NP   = Geometry::nPlane();
  const int_t    nel  = Geometry::nElmt();
  const int_t    npnp = Geometry::nTotElmt();
  const int_t    next = Geometry::nExtElmt();
  const int_t    ntot = Geometry::nBnode();
  char           buf[StrMax], err[StrMax], file[StrMax];
  ifstream       num;
  vector<real_t> work (npnp);
  real_t         *mass, *unity = &work[0];
  int_t i, j, nset, nglobal;
  int_t *gid, *q;

  // -- Read in NumberSystems.

  strcat   (strcpy (file, session), ".num");
  num.open (file);

  if (!num) {
    ROOTONLY {
      sprintf (buf, "enumerate -O3 %s > %s", session, file);
      if (system (buf)) {
        sprintf (err, "couldn't open session file %s, or %s", session, file);
        message (routine, err, ERROR);
      }
    }

    Femlib::synchronize();      // -- Ensure creation has completed.
    num.clear ();
    num.open  (file);

    if (!num) {
      sprintf (err, "couldn't find or create number file %s", file);
      message (routine, err, ERROR);
    }
  }

  num >> buf >> buf >> buf >> buf;
  j = strlen (buf);
  for (i = 0; i < j; i++)
    if (!strchr (_fields, tolower (buf[i]))) {
      sprintf (err, "Fields nominated in %s.num (\"%s\") don't match \"%s\"",
	       session, buf, _fields);
      message (routine, err, ERROR);
    }

  num.getline(buf, StrMax).getline(buf, StrMax);

  num >> buf >> nset >> buf >> buf >> buf;

  _numsys.resize (nset);
  for (i = 0; i < nset; i++) {
    num >> buf;
    _numsys[i] = new NumberSys;
    _numsys[i] -> _fields = new char [strlen (buf) + 1];
    strcpy (_numsys[i] -> _fields, buf);
  }

  for (i = 0; i < nset; i++)
    for (j = 0; j < nset; j++) {
      if (i == j) continue;
      if (strpbrk (_numsys[i] -> _fields, _numsys[j] -> _fields)) {
	sprintf (err, "Field name duplication: %s <--> %s",
		 _numsys[i] -> _fields, _numsys[j] -> _fields);
	message (routine, err, ERROR);
      }
    }

  num >> buf >> buf;
  if (strcmp (buf, "NEL")) {
    sprintf (err, "expected \"NEL\", found %s", buf);
    message (routine, err, ERROR);
  }
  num >> buf;
  for (i = 0; i < nset; i++) {
    num >> j;
    if (j != nel) {
      sprintf (err, "mismatch in number of elements: %1d vs. %1d", j, nel);
      message (routine, err, ERROR);
    }
  }

  num >> buf >> buf;
  if (strcmp (buf, "NP_MAX")) {
    sprintf (err, "expected \"NP_MAX\", found %s", buf);
    message (routine, err, ERROR);
  }
  num >> buf;
  for (i = 0; i < nset; i++) {
    num >> j;
    if (j != np) {
      sprintf (err, "mismatch in element order: %1d vs. %1d", j, np);
      message (routine, err, ERROR);
    }
  }

  num >> buf >> buf;
  if (strcmp (buf, "NEXT_MAX")) {
    sprintf (err, "expected \"NEXT_MAX\", found %s", buf);
    message (routine, err, ERROR);
  }
  num >> buf;
  for (i = 0; i < nset; i++) num >> j;

  num >> buf >> buf;
  if (strcmp (buf, "NINT_MAX")) {
    sprintf (err, "expected \"NINT_MAX\", found %s", buf);
    message (routine, err, ERROR);
  }
  num >> buf;
  for (i = 0; i < nset; i++) num >> j;

  num >> buf >> buf;
  if (strcmp (buf, "NTOTAL")) {
    sprintf (err, "expected \"NTOTAL\", found %s", buf);
    message (routine, err, ERROR);
  }
  num >> buf;
  for (i = 0; i < nset; i++) {
    num >> j;
    if (j != NP) {
      sprintf (err, "mismatch in Field storage requirements: %1d vs %1d",j,NP);
      message (routine, err, ERROR);
    }
  }

  num >> buf >> buf;
  if (strcmp (buf, "NBOUNDARY")) {
    sprintf (err, "expected \"NBOUNDARY\", found %s", buf);
    message (routine, err, ERROR);
  }
  num >> buf;
  for (i = 0; i < nset; i++) {
    num >> j;
    if (j != ntot) {
      sprintf (err,"mismatch in number of boundary nodes: %1d vs. %1d",j,ntot);
      message (routine, err, ERROR);
    }
    _numsys[i] -> _btog  = new int_t [ntot];
    _numsys[i] -> _bmask = new int_t [ntot];
  }

  num >> buf >> buf;
  if (strcmp (buf, "NGLOBAL")) {
    sprintf (err, "expected \"NGLOBAL\", found %s", buf);
    message (routine, err, ERROR);
  }
  num >> buf;
  for (i = 0; i < nset; i++) num >> _numsys[i] -> _nglobal;

  num >> buf >> buf;
  if (strcmp (buf, "NSOLVE")) {
    sprintf (err, "expected \"NSOLVE\", found %s", buf);
    message (routine, err, ERROR);
  }
  num >> buf;
  for (i = 0; i < nset; i++) num >> _numsys[i] -> _nsolve;

  num >> buf >> buf;
  if (strcmp (buf, "OPTIMIZATION")) {
    sprintf (err, "expected \"OPTIMIZATION\", found %s", buf);
    message (routine, err, ERROR);
  }
  num >> buf;
  for (i = 0; i < nset; i++) num >> _numsys[i] -> _optlev;

  num >> buf >> buf;
  if (strcmp (buf, "BANDWIDTH")) {
    sprintf (err, "expected \"BANDWIDTH\", found %s", buf);
    message (routine, err, ERROR);
  }
  num >> buf;
  for (i = 0; i < nset; i++) num >> _numsys[i] -> _nbandw;

  num.getline(buf, StrMax).getline(buf, StrMax).getline(buf, StrMax);

  // -- Consistency checks passed.  Now read in node numbers & BC masks.
  //    A bmask value of 1 implies an essential BC will be imposed at node.

  for (i = 0; i < ntot; i++) {
    num >> buf >> buf >> buf;
    for (j = 0; j < nset; j++) 
      num >> _numsys[j] -> _btog[i] >> _numsys[j] -> _bmask [i];
  }

  if (num.bad())
    message (routine, "failed reading to end of node-number file", ERROR);

  // -- Build emasks by inspecting bmasks for each element.
  //    An emask of 1 implies the element has at least one essential BC node.

  for (j = 0; j < nset; j++) {
    NumberSys* N = _numsys[j];
    N -> _emask = new int_t [nel];
    for (q = N -> _bmask, i = 0; i < nel; i++, q += next)
      N -> _emask[i] = Veclib::any (next, q, 1);
  }
#if 0
  // -- Create inverse mass matrices, avoid division by zero on axis.

  for (j = 0; j < nset; j++) {
    
    // -- Create diagonal global mass matrix.

    nglobal = _numsys[j] -> _nglobal;
    mass    = _numsys[j] -> _imass = new real_t [nglobal];
    Veclib::zero (nglobal, mass, 1);
    for (gid = _numsys[j] -> _btog, i = 0; i < nel; i++, gid += next) {
      Veclib::fill (npnp, 1.0, unity, 1);
      elmt[i] -> bndryDsSum (gid, unity, mass);
    }

    // -- Invert.

    for (i = 0; i < nglobal; i++) mass[i] = 1.0 / mass[i];
  }
#endif  
}

#endif

void BCmgr::buildsurf (FEML*             file,
		       vector<Element*>& Elmt)
// ---------------------------------------------------------------------------
// Private member function.
//
// Assemble list of element sides that have boundary conditions
// attached.  Element and side numbers are decremented by one, i.e are
// zero indexed in internal storage.
//
// As a part of internal checking, we want to ensure that the mesh for
// all "axis" group BCs has y=0.
// ---------------------------------------------------------------------------
{
  const char  routine[] = "BCmgr::buildsurf";
  const int_t nsurf = file -> attribute ("SURFACES", "NUMBER");
  char        err[StrMax], tag[StrMax], group;
  int_t       i, t, elmt, side;
  BCtriple*   BCT;
 
  for (i = 0; i < nsurf; i++) {
    while ((file->stream().peek()) == '#') file->stream().ignore(StrMax, '\n');

    file -> stream() >> t >> elmt >> side >> tag;

    if (strcmp (tag, "<B>") == 0) {
      
      file -> stream() >> group;

      BCT = new BCtriple;

      BCT -> group = group;
      BCT -> elmt  = --elmt;
      BCT -> side  = --side;

      _elmtbc.insert (_elmtbc.end(), BCT);

      file -> stream() >> tag;
      if (strcmp (tag, "</B>") != 0) {
	sprintf (err, "Surface %1d: couldn't close tag <B> with %s", t, tag);
	message (routine, err, ERROR);
      }
    } else
      file -> stream().ignore (StrMax, '\n');
  }

  if (Geometry::cylindrical()) {
    const int_t                 np = Geometry::nP();
    vector<real_t>              work (np);
    vector<BCtriple*>::iterator b;

    for (b = _elmtbc.begin(); b != _elmtbc.end(); b++) {
      BCT = *b;
      if (strstr (groupInfo (BCT -> group), "axis")) {
	Elmt[BCT -> elmt] -> sideGetY (BCT -> side, &work[0]);
	for (i = 0; i < np; i++)
	  if (::fabs (work[i]) > EPSDP) {
	    sprintf (err,
		     "elmt: %1d, side: %1d, offset: %1d, "
		     "y value (%g) too large on axis BC",
		     BCT -> elmt + 1,
		     BCT -> side + 1,
		     i, work[i]);
	    message (routine, err, ERROR);
	  }
      }
    }
  }
}


int_t BCmgr::nWall ()
// ---------------------------------------------------------------------------
// Count up the number of surfaces/element edges that have "wall" descriptor.
// ---------------------------------------------------------------------------
{
  vector<BCtriple*>::const_iterator b;
  int_t                             count = 0;
  BCtriple*                         BCT;

  for (b = _elmtbc.begin(); b != _elmtbc.end(); b++) {
    BCT = *b;
    if (strstr (groupInfo (BCT -> group), "wall")) count++;
  }

  return count;
}


int_t BCmgr::nAxis ()
// ---------------------------------------------------------------------------
// Count up the number of surfaces/element edges that have "axis" descriptor.
// ---------------------------------------------------------------------------
{
  vector<BCtriple*>::const_iterator b;
  int_t                             count = 0;
  BCtriple*                         BCT;

  for (b = _elmtbc.begin(); b != _elmtbc.end(); b++) {
    BCT = *b;
    if (strstr (groupInfo (BCT -> group), "axis")) count++;
  }

  return count;
}


int_t BCmgr::nMatching (const char* descript)
// ---------------------------------------------------------------------------
// Count up the number of surfaces/element edges that have given
//  string descriptor.
// ---------------------------------------------------------------------------
{
  vector<BCtriple*>::const_iterator b;
  int_t                             count = 0;
  BCtriple*                         BCT;

  for (b = _elmtbc.begin(); b != _elmtbc.end(); b++) {
    BCT = *b;
    if (strstr (groupInfo (BCT -> group), descript)) count++;
  }

  return count;
}


///////////////////////////////////////////////////////////////////////////////
// Code following implements computed BC types.
///////////////////////////////////////////////////////////////////////////////


void BCmgr::buildComputedBCs (const Field* master)
// ---------------------------------------------------------------------------
// Build class-scope storage structures required for evaluation of various
// kinds of computed Boundary conditions for Navier--Stokes type problems.
//
// There is some wastage as memory is also allocated for locations that will
// not have computed BCs.
// ---------------------------------------------------------------------------
{
  int_t i, j;

  _nLine = master -> _nline;
  _nEdge = master -> _nbound;
  _nZ    = master -> _nz;
  _nP    = Geometry::nP();
  _nTime = Femlib::ivalue ("N_TIME");

  _work = new real_t [static_cast<size_t>
		      (5*sqr(_nP)                  +
		       7*_nP                       +
		       Integration::OrderMax+1     +
		       _nLine                      )];
  
  _fbuf = new real_t [static_cast<size_t>(_nLine * _nZ)];
  _u2   = new real_t [static_cast<size_t>(_nLine * _nZ)]; // Ref [3]
  _unp  = new real_t [static_cast<size_t>(_nLine * _nZ)]; // Ref [3]
  _Enux = new real_t [static_cast<size_t>(_nLine * _nZ)]; // Ref [3]
  _Enuy = new real_t [static_cast<size_t>(_nLine * _nZ)]; // Ref [3]

  // -- The structure of each of the following arrays is
  //    _xx[time_level][z_plane] which evaluates to a real_t* that is
  //    a pointer to _nline real_t storage.  Thus _xx[0] is an
  //    equivalent pointer to Field->_line.

  // -- Consult references regarding the uses of the following storage
  //    used for computed BCs.

  _u     = new real_t** [static_cast<size_t>(_nTime)]; // -- Velocity component.
  _v     = new real_t** [static_cast<size_t>(_nTime)];
  _w     = new real_t** [static_cast<size_t>(_nTime)];

  _uhat  = new real_t** [static_cast<size_t>(_nTime)]; // -- in Fourier.
  _vhat  = new real_t** [static_cast<size_t>(_nTime)];
  _what  = new real_t** [static_cast<size_t>(_nTime)];
  
  _un    = new real_t** [static_cast<size_t>(_nTime)]; // -- u.n.
  _divu  = new real_t** [static_cast<size_t>(_nTime)]; // -- div(u).
  _gradu = new real_t** [static_cast<size_t>(_nTime)]; // -- n.grad(u).n.
  _hopbc = new real_t** [static_cast<size_t>(_nTime)]; // -- grad(p).n.
  _ndudt = new real_t** [static_cast<size_t>(_nTime)]; // -- (du/dt).n.
  
  for (i = 0; i < _nTime; i++) {
    _u    [i] = new real_t* [static_cast<size_t>(_nZ)];
    _v    [i] = new real_t* [static_cast<size_t>(_nZ)];
    _w    [i] = new real_t* [static_cast<size_t>(_nZ)]; 
    _uhat [i] = new real_t* [static_cast<size_t>(_nZ)];
    _vhat [i] = new real_t* [static_cast<size_t>(_nZ)];
    _what [i] = new real_t* [static_cast<size_t>(_nZ)];       
    _un   [i] = new real_t* [static_cast<size_t>(_nZ)];
    _divu [i] = new real_t* [static_cast<size_t>(_nZ)];
    _gradu[i] = new real_t* [static_cast<size_t>(_nZ)];
    _hopbc[i] = new real_t* [static_cast<size_t>(_nZ)];
    _ndudt[i] = new real_t* [static_cast<size_t>(_nZ)];    
    
    _u    [i][0] = new real_t [static_cast<size_t>(_nLine * _nZ)];
    _v    [i][0] = new real_t [static_cast<size_t>(_nLine * _nZ)];
    _w    [i][0] = new real_t [static_cast<size_t>(_nLine * _nZ)];
    _uhat [i][0] = new real_t [static_cast<size_t>(_nLine * _nZ)];
    _vhat [i][0] = new real_t [static_cast<size_t>(_nLine * _nZ)];
    _what [i][0] = new real_t [static_cast<size_t>(_nLine * _nZ)];    
    _un   [i][0] = new real_t [static_cast<size_t>(_nLine * _nZ)];
    _divu [i][0] = new real_t [static_cast<size_t>(_nLine * _nZ)];
    _gradu[i][0] = new real_t [static_cast<size_t>(_nLine * _nZ)];
    _hopbc[i][0] = new real_t [static_cast<size_t>(_nLine * _nZ)];
    _ndudt[i][0] = new real_t [static_cast<size_t>(_nLine * _nZ)];    
    
    Veclib::zero (_nLine * _nZ, _u    [i][0], 1);
    Veclib::zero (_nLine * _nZ, _v    [i][0], 1);
    Veclib::zero (_nLine * _nZ, _w    [i][0], 1);
    Veclib::zero (_nLine * _nZ, _uhat [i][0], 1);
    Veclib::zero (_nLine * _nZ, _vhat [i][0], 1);
    Veclib::zero (_nLine * _nZ, _what [i][0], 1);    
    Veclib::zero (_nLine * _nZ, _un   [i][0], 1);
    Veclib::zero (_nLine * _nZ, _divu [i][0], 1);
    Veclib::zero (_nLine * _nZ, _gradu[i][0], 1);
    Veclib::zero (_nLine * _nZ, _hopbc[i][0], 1);
    Veclib::zero (_nLine * _nZ, _ndudt[i][0], 1);    

    for (j = 1; j < _nZ; j++) {
      _u    [i][j] = _u    [i][0] + j * _nLine;
      _v    [i][j] = _v    [i][0] + j * _nLine;
      _w    [i][j] = _w    [i][0] + j * _nLine;
      _uhat [i][j] = _uhat [i][0] + j * _nLine;
      _vhat [i][j] = _vhat [i][0] + j * _nLine;
      _what [i][j] = _what [i][0] + j * _nLine;      
      _un   [i][j] = _un   [i][0] + j * _nLine;
      _divu [i][j] = _divu [i][0] + j * _nLine;
      _gradu[i][j] = _gradu[i][0] + j * _nLine;
      _hopbc[i][j] = _hopbc[i][0] + j * _nLine;
      _ndudt[i][j] = _ndudt[i][0] + j * _nLine;      
    }
  }
}


void BCmgr::maintainPhysical (const Field*             master,
			      const vector<AuxField*>& Us    ,
			      const int_t              nCom  )
// ---------------------------------------------------------------------------
// Update storage that (on account of requiring products of terms)
// must be computed from data which are currently in physical space.
//
// nCom is the number of velocity components.
//
// Note that we only need to use this function if we have open BCs.  
// ---------------------------------------------------------------------------
{
  if (!_open) return;		   // -- Hence _toggle remains false always.
  
  const vector<Boundary*>& BC = master -> _bsys -> getBCs (0);
  Boundary*                B;
  int_t                    i, j, k, offset, skip;

  // -- Roll up time-ordered internal storage stack, initialise storage.

  rollv (_u, _nTime);
  rollv (_v, _nTime);
  rollv (_w, _nTime);

  Veclib::zero (_nLine * _nZ, _w[0][0], 1);

  // -- Save domain-edge velocity components for later extrapolation.

  for (k = 0; k < _nZ; k++) {
    for (i = 0; i < _nEdge; i++) {
      B      = BC[i];
      offset = B -> dOff ();
      skip   = B -> dSkip();
      j      = i * _nP;
      Veclib::copy (_nP, Us[0] -> _plane[k] + offset, skip, _u[0][k] + j, 1);
      Veclib::copy (_nP, Us[1] -> _plane[k] + offset, skip, _v[0][k] + j, 1);
      if (nCom == 3)
	Veclib::copy (_nP, Us[2] -> _plane[k] + offset, skip, _w[0][k] + j, 1);
    }
  }

  _toggle = true;  // -- Set trigger for later forward Fourier transform.
}


void BCmgr::maintainFourier (const int_t      step   ,
			     const Field*     master ,
			     const AuxField** Us     ,
			     const AuxField** Uf     ,
			     const int_t      nCom   ,
			     const bool       timedep)
// ---------------------------------------------------------------------------
// Update storage for evaluation of internally computed boundary
// conditions, see Refs [1-3].  Storage order for each edge represents
// a CCW traverse of element boundaries.
//
// If the velocity field varies in time on HOPBC field boundaries
// (e.g. due to time-varying BCs) the local fluid acceleration will be
// estimated from input velocity fields by explicit extrapolation if
// timedep is true.  This correction cannot be carried out at the
// first timestep, since the required extrapolation cannot be done.
// If the acceleration is known, (for example, a known reference frame
// acceleration) it is better to leave timedep false, and to use
// BCmgr::accelerate() to add in the accelerative term.  Note also
// that since grad P is dotted with n, the unit outward normal, at a
// later stage, timedep only needs to be set if there are wall-normal
// accelerative terms.  NB: The default value of timedep is true.
//
// Field* master gives a list of egdes with which to traverse storage
// areas (note this assumes equal-order interpolations).
//
// No smoothing is done to high-order spatial derivatives computed here.
//
// In order to be able to deal with both the original computed-Neumann
// HOPBC BCs appropriate for solid walls (as described e.g. in Refs
// [1,2]) and the computed-mixed pressure open BCs described in Ref
// [3], eq. (37), we separately compute n.d(u)/dt and n.[N(u)+ f - \nu
// curl(curl(u))] where N(u) = -u.grad(u), i.e. the nonlinear terms
// tranposed to RHS of Navier-Stokes.  N(u)+f is supplied by input Uf.
// When required for HOPBC comuputation (elsewhere) we subtract off
// n.d(u)/dt.
//
// If in doubt about (checking) the signs of terms: make up dp/dn from NSE.
// ---------------------------------------------------------------------------
{
  const int_t              base  = Geometry::baseMode();
  const int_t              nMode = Geometry::nModeProc();
  const int_t              mLo   = (Geometry::procID() == 0) ? 1 : 0;
  const int_t              Je    = min (step - 1, _nTime);

  const real_t             nu    = Femlib::value ("KINVIS");
  const real_t             invDt = 1.0 / Femlib::value ("D_T");

  const vector<Boundary*>& BC    = master -> _bsys -> getBCs (0);

  const AuxField*          Ux    = Us[0];
  const AuxField*          Uy    = Us[1];
  const AuxField*          Uz    = (nCom == 3) ? Us[2] : 0;
  const AuxField*          Nx    = Uf[0];
  const AuxField*          Ny    = Uf[1];

  Boundary*                B;
  int_t                    i, j, k, m, q, kr, ki, offset, skip;

  real_t                   *UxRe, *UxIm, *UyRe, *UyIm, *UzRe, *UzIm;
  real_t*                  wrk   = _work;
  real_t*                  xr    = wrk + 5*sqr(_nP) + 3*_nP;
  real_t*                  xi    = xr  + _nP;
  real_t*                  yr    = xi  + _nP;
  real_t*                  yi    = yr  + _nP;
  real_t*                  alpha = yi  + _nP;

  real_t                   betaK;

  if (Je) Integration::StifflyStable (Je, alpha);

  // -- Roll up time-ordered internal storage stacks (except _un),
  //    initialise storage.

  rollv (_uhat,  _nTime);
  rollv (_vhat,  _nTime);
  rollv (_what,  _nTime);  
  rollv (_divu,  _nTime);
  rollv (_gradu, _nTime);
  rollv (_hopbc, _nTime);
  rollv (_ndudt, _nTime);
  
  Veclib::zero (_nLine*_nZ, _divu [0][0], 1);  
  Veclib::zero (_nLine*_nZ, _gradu[0][0], 1);  
  Veclib::zero (_nLine*_nZ, _hopbc[0][0], 1);
  Veclib::zero (_nLine*_nZ, _ndudt[0][0], 1);  

  // -- First (while we still have all the previous time levels of u.n
  //    available), make the most recent u.n and compute the
  //    time-derivative contribution n.dot.d(u)/dt to the high order
  //    pressure BCs.  We use the top level of _divu as temporary
  //    storage of u.n.  Then add in the n.[N+f] terms to _hopbc.

  for (i = 0; i < _nEdge; i++) {
    B      = BC[i];
    offset = B -> dOff ();
    skip   = B -> dSkip();
    j      = i * _nP;

    for (k = 0; k < _nZ; k++) {
      ROOTONLY if (k == 1) continue;

      Veclib::copy  (_nP, Ux -> _plane[k] + offset, skip, _uhat[0][k] + j, 1);	
      Veclib::copy  (_nP, Uy -> _plane[k] + offset, skip, _vhat[0][k] + j, 1);
      if (nCom == 3)
	Veclib::copy  (_nP, Uz -> _plane[k] + offset, skip, _what[0][k] + j, 1);
      
      Veclib::vvtvvtp (_nP, _uhat[0][k] + j, 1, B->nx(), 1,
		            _vhat[0][k] + j, 1, B->ny(), 1,
		            _divu[0][k] + j, 1);

      if (Je && timedep) {
	// -- Make estimate of d(u.n)/dt at end of last time step.
	Veclib::smul (_nP, alpha[0], _divu[0][k] + j , 1, _work, 1);
	for (q = 0; q < Je; q++)
	  Blas::axpy (_nP, alpha[q + 1], _un[q][k] + j, 1, _work, 1);
	Veclib::smul (_nP, invDt, _work, 1, _ndudt[0][k] + j, 1);
      }

      // -- Deal with nonlinear terms + forcing.

      Veclib::zero  (_nP, _work, 1);
      Veclib::vvtvp (_nP, Nx -> _plane[k] + offset, skip, B -> nx(), 1,
		     _work, 1, _work , 1);
      Veclib::vvtvp (_nP, Ny -> _plane[k] + offset, skip, B -> ny(), 1,
		     _work, 1, _work, 1);

      // -- For cylindrical coordinates, N_ are radius-premultiplied. Cancel.

      if (Geometry::cylindrical()) B -> divY (_work);

      Veclib::vadd (_nP, _work, 1, _hopbc[0][k] + j, 1, _hopbc[0][k] + j, 1);
    }
  }

  // -- NOW we can roll _un, load it from _divu, and clear _divu.

  rollv (_un, _nTime);
  Veclib::copy (_nLine*_nZ, _divu[0][0], 1, _un[0][0], 1);
  Veclib::zero (_nLine*_nZ, _divu[0][0], 1);  

  // -- Next, all the terms for which gradients are required. For each
  //    element which has a boundary edge. These make contributions to
  //    _divu, _gradu and _hopbc.

  // -- First, deal with HOPBC -\nu*curlCurl(u) terms.  Although this
  //    produces some redundant gradient operations, any later
  //    gradients will be element-edge-only, hence cheap.  And the
  //    Edge:curlCurl() method is well tested.

  for (i = 0; i < _nEdge; i++) {
    B  = BC[i];
    j  = i * _nP;

    ROOTONLY {			    // -- Deal with 2D/zero Fourier mode terms.
      UxRe = Ux -> _plane[0];
      UyRe = Uy -> _plane[0];

      B -> curlCurl (0,UxRe,0,UyRe,0,0,0,xr,0,yr,0,wrk);

      Veclib::svvttvp(_nP,-nu,xr,1,B->nx(),1,_hopbc[0][0]+j,1,_hopbc[0][0]+j,1);
      Veclib::svvttvp(_nP,-nu,yr,1,B->ny(),1,_hopbc[0][0]+j,1,_hopbc[0][0]+j,1);
    }

    for (m = mLo; m < nMode; m++) { // -- Higher modes.
      kr = 2 * m;
      ki = kr + 1;

      UxRe = Ux -> _plane[kr] ;
      UxIm = Ux -> _plane[ki];
      UyRe = Uy -> _plane[kr];
      UyIm = Uy -> _plane[ki];
      UzRe = Uz -> _plane[kr];
      UzIm = Uz -> _plane[ki];

      B -> curlCurl (m+base,UxRe,UxIm,UyRe,UyIm,UzRe,UzIm,xr,xi,yr,yi,wrk);

      Veclib::svvttvp 
	(_nP, -nu, xr,1, B->nx(),1, _hopbc[0][kr]+j,1, _hopbc[0][kr]+j,1);
      Veclib::svvttvp 
	(_nP, -nu, xi,1, B->nx(),1, _hopbc[0][ki]+j,1, _hopbc[0][ki]+j,1);
      Veclib::svvttvp 
	(_nP, -nu, yr,1, B->ny(),1, _hopbc[0][kr]+j,1, _hopbc[0][kr]+j,1);
     Veclib::svvttvp 
	(_nP, -nu, yi,1, B->ny(),1, _hopbc[0][ki]+j,1, _hopbc[0][ki]+j,1);
    }
  }

  if (!_open) return;		   // -- Nothing else to do in this case.

  // -- Remaining gradient-based terms (ones from Ref [3]).

  for (i = 0; i < _nEdge; i++) {
    B      = BC[i];
    offset = B -> dOff ();
    skip   = B -> dSkip();
    j      = i * _nP;

    ROOTONLY {			    // -- Deal with 2D/zero Fourier mode terms.
      UxRe = Ux -> _plane[0];
      UyRe = Uy -> _plane[0];

      B -> sideGrad   (UxRe, xr, yr, wrk);
      Veclib::vadd    (_nP,  xr, 1, _divu[0][0]+j, 1, _divu[0][0]+j, 1);

      Veclib::vvtvvtp (_nP, xr,1,B->nx(),1, yr,1, B->ny(),1, wrk,1);
      Veclib::vvtvp   (_nP,wrk,1,B->nx(),1,_gradu[0][0]+j,1,_gradu[0][0]+j,1);

      B -> sideGrad   (UyRe, xr, yr, wrk);
      Veclib::vadd    (_nP,  yr, 1, _divu[0][0]+j, 1, _divu[0][0]+j, 1);

      Veclib::vvtvvtp (_nP, xr,1,B->nx(),1, yr,1, B->ny(),1, wrk,1);
      Veclib::vvtvp   (_nP,wrk,1,B->ny(),1,_gradu[0][0]+j,1,_gradu[0][0]+j,1);

      if (Geometry::cylindrical()) {
	B -> get  (UyRe, wrk);
	B -> divY (wrk);
	Veclib::vadd (_nP, wrk, 1, _divu[0][0]+j, 1, _divu[0][0]+j, 1);
      }
    }

    for (m = mLo; m < nMode; m++) { // -- Higher modes.
      betaK = m * Femlib::value ("BETA");
      kr    = 2 * m;
      ki    = kr + 1;

      UxRe = Ux -> _plane[kr];
      UyIm = Uy -> _plane[ki];
      UyRe = Uy -> _plane[kr];
      UyIm = Uy -> _plane[ki];
      UzRe = Uz -> _plane[kr];
      UzIm = Uz -> _plane[ki];

      B -> sideGrad   (UxRe, xr, yr, wrk);
      Veclib::vadd    (_nP,  xr, 1, _divu[0][kr]+j, 1, _divu[0][kr]+j, 1);

      Veclib::vvtvvtp (_nP, xr,1, B->nx(),1, yr,1,B->ny(),1, wrk,1);
      Veclib::vvtvp   (_nP,wrk,1,B->nx(),1,_gradu[0][kr]+j,1,_gradu[0][kr]+j,1);

      B -> sideGrad   (UyRe, xr, yr, wrk);
      Veclib::vadd    (_nP,  yr, 1, _divu[0][kr]+j, 1, _divu[0][kr]+j, 1);

      Veclib::vvtvvtp (_nP, xr,1, B->nx(),1, yr,1,B->ny(),1, wrk,1);
      Veclib::vvtvp   (_nP,wrk,1,B->ny(),1,_gradu[0][kr]+j,1,_gradu[0][kr]+j,1);

      if (Geometry::cylindrical()) {
	B -> get     (UyRe, wrk);
	B -> divY    (wrk);
	Veclib::vadd (_nP, wrk, 1, _divu[0][kr]+j, 1, _divu[0][kr]+j, 1);
	B -> get     (UzIm, wrk);
	B -> divY    (wrk);
	Blas::axpy   (_nP, -betaK, wrk, 1, _divu[0][kr]+j, 1);
      } else {
	B -> get     (UzIm, wrk);
	Blas::axpy   (_nP, -betaK, wrk, 1, _divu[0][kr]+j, 1);
      }

      B -> sideGrad   (UxIm, xr, yr, wrk);
      Veclib::vadd    (_nP, xr, 1, _divu[0][ki]+j, 1, _divu[0][ki]+j, 1);

      Veclib::vvtvvtp (_nP, xr,1, B->nx(),1, yr,1,B->ny(),1, wrk,1);
      Veclib::vvtvp   (_nP,wrk,1,B->nx(),1,_gradu[0][ki]+j,1,_gradu[0][ki]+j,1);

      B -> sideGrad   (UyIm, xr, yr, wrk);
      Veclib::vadd    (_nP, yr, 1, _divu[0][ki]+j, 1, _divu[0][ki]+j, 1);

      Veclib::vvtvvtp (_nP, xr,1, B->nx(),1, yr,1,B->ny(),1, wrk,1);
      Veclib::vvtvp   (_nP,wrk,1,B->ny(),1,_gradu[0][ki]+j,1,_gradu[0][ki]+j,1);

      if (Geometry::cylindrical()) {
	B -> get     (UyIm, wrk);
	B -> divY    (wrk);
	Veclib::vadd (_nP, wrk, 1, _divu[0][ki]+j, 1, _divu[0][ki]+j, 1);
	B -> get     (UzRe, wrk);
	B -> divY    (wrk);
	Blas::axpy   (_nP,  betaK, wrk, 1, _divu[0][ki]+j, 1);
      } else {
	B -> get     (UzRe, wrk);
	Blas::axpy   (_nP,  betaK, wrk, 1, _divu[0][ki]+j, 1);
      }
    }
  }
  Blas::scal (_nLine * _nZ, nu, _divu [0][0], 1);
  Blas::scal (_nLine * _nZ, nu, _gradu[0][0], 1);
}


void BCmgr::evaluateCNBCp (const int_t id   ,
			   const int_t plane,
			   const int_t step ,
			   real_t*     tgt  )
// ---------------------------------------------------------------------------
// "CNBCp" is an acronym for "Computed Natural BC pressure".  A.k.a.
// "high-order" pressure BC.  
//  
//  Refer Ref [1].  Load pressure (p) BC value (tgt) with values
//  obtained from HOPBC multi-level storage. Evaluation is confined to
//  a single element edge: parameter id tells us which this
//  corresponds to in our internal storage.
//
// The boundary condition for evaluation is (obtained by
// re-arrangement of n.NSE):
//
//   dP       /                              du  \
//   -- = n . | N(u) + f - \nu*curlCurl(u) - --  |  =  n . grad P.
//   dn   ~   \ ~ ~    ~                ~    dt  /     ~
//
// Grad P.n is estimated at the end of the current timestep using
// explicit extrapolation.  All the terms for this gradient were
// already stored, here we do the extrapolation only.
//
// Note: the storage for d(u.n)/dt is held separately from the other
// terms (which may also be used in computed mixed pressure BC from
// Ref [3], see following).
// ---------------------------------------------------------------------------
{
  if (step < 1) return;		// -- No evaluation during Field creation.

  Veclib::zero (_nP, tgt, 1);

  ROOTONLY if (plane == 1) return;    // -- Skip Nyquist.

  const int_t Je     = min (step, _nTime);
  const int_t offset = id * _nP;
  real_t*     beta   = _work;
  int_t       q;

  Integration::Extrapolation (Je, beta);
  
  for (q = 0; q < Je; q++) {
    Blas::axpy (_nP,  beta[q], _hopbc[q][plane] + offset, 1, tgt, 1);
    Blas::axpy (_nP, -beta[q], _ndudt[q][plane] + offset, 1, tgt, 1);
  }
}


void BCmgr::accelerate (const Vector& a,
			const Field*  u)
// ---------------------------------------------------------------------------
// Add in frame acceleration term on boundaries that correspond to
// essential velocity BCs (a = - du/dt).  Note that the acceleration
// time level should correspond to the time level in the most recently
// updated pressure gradient storage.  Work only takes place on zeroth
// Fourier mode.
//
// Yes, this is a HACK! To be used if the linear frame acceleration is
// externally available from elsehere, e.g. in fluid-structure
// interaction.  Nowadays we'd probably deal with that in FieldForce
// routines (and let it come in via Uf in maintainFourier).
// ---------------------------------------------------------------------------
{
  const vector<Boundary*>& BC = u -> _bsys -> getBCs (0);
  Boundary*                B;
  int_t                    i;

  for (i = 0; i < u -> _nbound; i++) {
    B = BC[i];
    B -> dotInForGroup ("velocity", a, _hopbc[0][0] + i*_nP);
  }
}


void BCmgr::evaluateCMBCp (const Field* master, // Used for list of boundaries.
			   const int_t  id    , // Index of this boundary edge.
			   const int_t  plane , // Index of Fourier plane.
			   const int_t  step  , // Time step.
			   real_t*      tgt   )
// ---------------------------------------------------------------------------
// "CMBCp" is an acronym for "Computed Mixed BC pressure".
//  
// Refer Ref [3] eq. (37).  Load mixed/Robin pressure (p) BC value
// (tgt) with values obtained from multi-level storage. Evaluation is
// confined to a single element edge: parameter id tells us which this
// corresponds to in our internal storage.
//
// But also, we use a call to this method to trigger computation of
// E(n,u*) terms which rely on products made in physical space (and
// subsequently Fourier transformed), which will also later be used
// by the associated mixed velocity BC.  These computations are
// toggled on the first call to this method in the current timestep.
// ---------------------------------------------------------------------------
{
  const int_t Je = min (step, _nTime);

  if (Je < 1) return;		// -- No evaluation during Field creation.

  Veclib::zero (_nP, tgt, 1);

  static const real_t      iUZD   = 1.0 / Femlib::value ("DONG_UODELTA");
  static const real_t      iNUD   = 1.0 / Femlib::value ("KINVIS*DONG_DO");
  //  static const real_t      BIAS   = Femlib::value ("DONG_BIAS");
  const vector<Boundary*>& BC     = master -> _bsys -> getBCs (0);
  const int_t              offset = id * _nP;
  const Boundary*          B;
  real_t*                  beta = _work;
  int_t                    i, j, k, q;

  Integration::Extrapolation (Je, beta);

  if (_toggle) {
		
    // -- First call to this routine in current timestep.  Make
    //    E(n,u*) i.e. Ref [3], eq. (18).  In our system, it is a 2C
    //    (x,y) vector at each point along the boundary.  Its
    //    computation relies on products of fields, which, as usual,
    //    will be evaluated in physical space, then finally Fourier
    //    transformed.  The velocity fields involved, u*, are values
    //    obtained by extrapolation to the end of the timestep.

    const int_t nTot  = _nLine * _nZ;  
    const int_t nZtot = Geometry::nZ();
    const int_t nPR   = Geometry::nProc();
    const int_t nLP   = _nLine / nPR;
    
    Veclib::zero (nTot, _u2,   1);
    Veclib::zero (nTot, _fbuf, 1);
    Veclib::zero (nTot, _unp,  1);

    for (q = 0; q < Je; q++) {
      Blas::axpy (nTot, beta[q], _u[q][0], 1, _u2,   1); // -- u*.
      Blas::axpy (nTot, beta[q], _v[q][0], 1, _fbuf, 1); // -- v*.
    }

    for (i = 0; i < _nEdge; i++) { // -- Make u*.n, leave in _unp.
      B = BC[i];
      j = i * _nP;
      for (k = 0; k < _nZ; k++)
	Veclib::vvtvvtp (_nP,
			 B -> nx(), 1, _u2   + k*_nLine + j, 1,
			 B -> ny(), 1, _fbuf + k*_nLine + j, 1,
			               _unp  + k*_nLine + j, 1);
    }

    // -- Use u*.n to start construction of E(n,u*), x,y components
    //    (what follows is (n.u*)u*). Ref. [3] eq. (18).

    Veclib::vmul (nTot, _unp, 1, _u2,   1, _Enux, 1);
    Veclib::vmul (nTot, _unp, 1, _fbuf, 1, _Enuy, 1);
    
    // -- Complete making u^2*.  Dealing with w* is potentially wasteful, but.

    Veclib::vmul  (nTot, _u2,   1, _u2,   1, _u2, 1);
    Veclib::vvtvp (nTot, _fbuf, 1, _fbuf, 1, _u2, 1, _u2, 1);

    Veclib::zero (nTot, _fbuf, 1);
    for (q = 0; q < Je; q++)
      Blas::axpy (nTot, beta[q], _w[q][0], 1, _fbuf, 1); // -- w*.

    Veclib::vvtvp (nTot, _fbuf, 1, _fbuf, 1, _u2, 1, _u2, 1);

    // -- Save switch function Theta0=0.25(1-tanh(u*.n/DONG_UODELTA)) in _fbuf.
    //    N.B. Extra factor of 0.5 incorporated here (0.25 = 0.5*0.5).

    Veclib::smul  (nTot, iUZD, _unp,  1, _fbuf, 1);
    //    Veclib::sadd  (nTot, BIAS, _unp,  1, _fbuf, 1);
    Veclib::vtanh (nTot, _fbuf, 1,       _fbuf, 1);
    Veclib::ssub  (nTot,  1.0, _fbuf, 1, _fbuf, 1);
    Blas::scal    (nTot,  0.25,          _fbuf, 1);

    // -- Finish making E(n,u*).

    for (i = 0; i < _nEdge; i++) {
      B = BC[i];
      j = i * _nP;
      for (k = 0; k < _nZ; k++) {
	Veclib::vvtvp (_nP, B -> nx(), 1, _u2   + k*_nLine + j, 1,
		       _Enux + k*_nLine + j, 1, _Enux + k*_nLine + j, 1);
	Veclib::vvtvp (_nP, B -> ny(), 1, _u2   + k*_nLine + j, 1,
		       _Enuy + k*_nLine + j, 1, _Enuy + k*_nLine + j, 1);
      }
    }

    // -- Multiply with Theta0 to make E(n,u*) in physical space.

    Veclib::vmul (nTot, _Enux, 1, _fbuf, 1, _Enux, 1);
    Veclib::vmul (nTot, _Enuy, 1, _fbuf, 1, _Enuy, 1);    

    // -- Forward Fourier transform the outcome.

    if (nPR == 1) {
      if (nZtot > 1) {
	if (nZtot == 2) {
	  Veclib::zero (_nLine, _Enux + _nLine, 1);
	  Veclib::zero (_nLine, _Enuy + _nLine, 1);
	} else {
	  Femlib::DFTr (_Enux,  nZtot, _nLine, FORWARD);
	  Femlib::DFTr (_Enuy,  nZtot, _nLine, FORWARD);
	}
      }
    } else {
      Femlib::exchange (_Enux, _nZ,    _nLine, FORWARD);
      Femlib::DFTr     (_Enux,  nZtot,  nLP,   FORWARD);
      Femlib::exchange (_Enux, _nZ,    _nLine, INVERSE);
      
      Femlib::exchange (_Enuy, _nZ,    _nLine, FORWARD);
      Femlib::DFTr     (_Enuy,  nZtot,  nLP,   FORWARD);
      Femlib::exchange (_Enuy, _nZ,    _nLine, INVERSE);
    }

    _toggle = false;		// -- Turn off extrapolation/Fourier.
  }

  ROOTONLY if (plane == 1) return; // -- Skip Nyquist plane.

  // -- Assemble RHS of Ref [3] eq. (37) in tgt (already set to zero).

  B = BC[id];
  Veclib::vvtvvtp (_nP,
		   B -> nx(), 1, _Enux + plane*_nLine + offset, 1,
		   B -> ny(), 1, _Enuy + plane*_nLine + offset, 1,
		   tgt, 1);
  
  for (q = 0; q < Je; q++)
    Blas::axpy (_nP, -beta[q], _gradu[q][plane] + offset, 1, tgt, 1);
  
  Blas::scal (_nP, -iNUD, tgt, 1);
  
  for (q = 0; q < Je; q++)
    Blas::axpy (_nP,  beta[q], _hopbc[q][plane] + offset, 1, tgt, 1);
}


void BCmgr::evaluateCMBCu (const Field* P    , // Pressure field.
			   const int_t  id   , // Index of this boundary.
			   const int_t  k    , // Index of Fourier plane.
			   const int_t  step , // Time step.
			   const char   cmpt , // Name of velocity component.
			   real_t*      tgt  )
// ---------------------------------------------------------------------------
// "CMBCu" is an acronym for "Computed Mixed BC velocity".  Refer Ref [3].
//  
// Load mixed velocity BC value (tgt) with values obtained from
// multi-level storage to make RHS of eq. (38). Evaluation is confined
// to a single element edge: parameter id tells us which this
// corresponds to in our internal storage.  Input parameter cmpt can
// be 'u', 'v', or 'w'.  Operations are almost exactly the same for
// 'u' and 'v'; the only distinction is whether we post-multiply the
// outcomes by unit outward normal components n_x (for 'u') or n_y
// ('v').  But for 'w' , n_z = 0, and we only have to deal with the
// extrapolative terms.
//
// Input parameter field P is the pressure at the end of the current timestep.
// 
// Note that storage areas for Fourier transform of KE fluxes _Enux
// and _Enuy have already been evaluated earlier in the timestep.
// ---------------------------------------------------------------------------
{
  const int_t Je = min (step, _nTime);

  if (Je < 1) return;		// -- No evaluation during Field creation.

  if (!strchr ("uvw", cmpt))
    message ("BCmgr::evaluateCMBCu", 
	     "called with unknown velocity component", ERROR);

  Veclib::zero (_nP, tgt, 1);

  ROOTONLY if (k == 1) return;	   // -- Skip Nyquist plane.

  static const real_t  iMu    = Femlib::value ("1.0/KINVIS");
  static const real_t  DOdt   = Femlib::value ("DONG_DO/D_T");
  const int_t          offset = id * _nP;
  const Boundary*      B      = P -> _bsys -> getBCs (0)[id];
  
  real_t* alpha = _work;
  real_t* beta  = alpha + Je + 1;
  int_t   q;

  Integration::StifflyStable (Je, alpha);
  Integration::Extrapolation (Je, beta);
  
  if (cmpt == 'u') {
    Veclib::copy (_nP, P -> _plane[k] + B -> dOff(), B -> dSkip(), tgt, 1);
    
    for (q = 0; q < Je; q++)
      Blas::axpy (_nP, -beta[q], _divu[q][k] + offset, 1, tgt, 1);
    
    Veclib::vmul  (_nP, B -> nx(), 1, tgt, 1, tgt, 1);
    Veclib::svvpt (_nP, iMu, _Enux + k*_nLine + offset, 1, tgt, 1, tgt, 1);
    
    for (q = 1; q <= Je; q++)
      Blas::axpy (_nP, -DOdt * alpha[q], _uhat[q-1][k] + offset, 1, tgt, 1);

  } else if (cmpt == 'v') {
    Veclib::copy (_nP, P -> _plane[k] + B -> dOff(), B -> dSkip(), tgt, 1);
    
    for (q = 0; q < Je; q++)
      Blas::axpy (_nP, -beta[q], _divu[q][k] + offset, 1, tgt, 1);
    
    Veclib::vmul  (_nP, B -> ny(), 1, tgt, 1, tgt, 1);
    Veclib::svvpt (_nP, iMu, _Enuy + k*_nLine + offset, 1, tgt, 1, tgt, 1);

    for (q = 1; q <= Je; q++)
      Blas::axpy (_nP, -DOdt * alpha[q], _vhat[q-1][k] + offset, 1, tgt, 1);
    
  } else {			// -- 'w'.
    Veclib::zero (_nP, tgt, 1);

    for (q = 1; q <= Je; q++)
      Blas::axpy (_nP, -DOdt * alpha[q], _what[q-1][k] + offset, 1, tgt, 1);
  }
}
