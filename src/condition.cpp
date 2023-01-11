//////////////////////////////////////////////////////////////////////////////
// condition.cpp: functions used to evaluate & apply BCs along element edges.
//
// Copyright (c) 1994+, Hugh M Blackburn
//
// All classes inherit the (semi-abstract) base class Condition.
// Owing to the different behaviour of essential and natural BCs, a
// number of routines take no action; due also to the generality of
// the base class function calls, many routines do not use all
// parameters.
//
// The two basic kinds of conditions are essential (Dirichlet) and
// natural (Neumann).  Essential conditions are set/imposed directly,
// while natural conditions are applied by summation of integral
// terms, owing to the fact that they derive from integration by
// parts.  So function "set" is only used by essential-type conditions
// while function "sum" is only used by natural-type conditions.
//
// The HOPBC condition are a special case of natural conditions used
// for the pressure Poisson equation: they are derived from the
// momentum equations and computed using an extrapolative process
// described in [1].
//
// A further general class of available condition is the Robin or
// mixed type, of form dc/dn + K( c - C ) = 0. See [2]. 

// NB: if you invent a new BC which is of essential/Dirichlet type,
// you also need to edit mesh.cpp/buildLiftMask() so that assemblymap
// sets the global numbering mask to 1 for this type.  Otherwise the
// BC won't be correctly applied.  In this consideration, Robin/mixed
// count as natural/Neumann.
//
// REFERENCES
// ----------
// [1] Karniadakis, Israeli & Orszag (1991), JCP 91(2)
// [2] Blackburn (2001), Comput Chem Eng 25(2/3)
// [3] Dong (2015), JCP 302:300-328 
// [4] Liu, Xie & Dong (2020), IJHFF 151:119355
///////////////////////////////////////////////////////////////////////////////

#include <sem.h>


void Essential::evaluate (const Field*   src    ,
			  const int_t    id     , 
			  const int_t    plane  , 
			  const Element* elmt   ,
			  const int_t    side   , 
			  const int_t    step   ,
			  const bool     Fourier,
			  real_t*        tgt    ) const
// ---------------------------------------------------------------------------
// Load external value storage area tgt with installed constant.
// Physical space.  The long list of inputs (most of which are not
// used here) is required because the same argument list is shared
// across all evaluation prototypes, which have diverse actions.
// E.g. elmt and side are used by EssentialFunction::evaluate().
// ---------------------------------------------------------------------------
{
  if (!Fourier) Veclib::fill (Geometry::nP(), _value, tgt, 1);
}


void Essential::set (const int_t   side,
		     const int_t*  bmap,
		     const real_t* src ,
		     real_t*       tgt ) const
// ---------------------------------------------------------------------------
// Scatter external value storage area src into globally-numbered tgt. 
// ---------------------------------------------------------------------------
{
  const int_t  nm    = Geometry::nP() - 1;
  const int_t* start = bmap;
  
  switch (side) {
  case 1: start += nm;           break;
  case 2: start += nm + nm;      break;
  case 3: start += nm + nm + nm; break;
  default: break;
  }

  Veclib::scatr (nm, src, start, tgt);
  if   (side == 3) tgt[bmap [ 0]] = src[nm];
  else             tgt[start[nm]] = src[nm];  
}


void Essential::describe (char* tgt) const
// ---------------------------------------------------------------------------
// Load descriptive/diagnostic material into tgt.  N.B. The text of
// these descriptors gets used elsewhere (e.g. bsys.cpp)
// ---------------------------------------------------------------------------
{
  sprintf (tgt, "essential:\t%g", _value);
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

EssentialFunction::EssentialFunction (const char* f)
// ---------------------------------------------------------------------------
// Essential condition that sets value by interpreting function
// string, which can be an explicit function of x, y, z & t as well as
// any installed token values.  Otherwise the same as Essential class.
// ---------------------------------------------------------------------------
{
  strcpy ((_function = new char [strlen (f) + 1]), f);
}


void EssentialFunction::evaluate (const Field*   src    ,
				  const int_t    id     , 
				  const int_t    plane  , 
				  const Element* elmt   ,
				  const int_t    side   , 
				  const int_t    step   ,
				  const bool     Fourier,
				  real_t*        tgt    ) const
// ---------------------------------------------------------------------------
// Load external tgt by interpreting function.
// Physical space.
// ---------------------------------------------------------------------------
{
  if (!Fourier) elmt -> sideEval (side, tgt, _function);
}


void EssentialFunction::set (const int_t   side,
			     const int_t*  bmap,
			     const real_t* src ,
			     real_t*       tgt ) const
// ---------------------------------------------------------------------------
// Scatter external value storage area src into globally-numbered tgt
// (as for Essential class).
// ---------------------------------------------------------------------------
{
  const int_t  nm    = Geometry::nP() - 1;
  const int_t* start = bmap;
  
  switch (side) {
  case 1: start += nm;           break;
  case 2: start += nm + nm;      break;
  case 3: start += nm + nm + nm; break;
  default: break;
  }

  Veclib::scatr (nm, src, start, tgt);
  if   (side == 3) tgt[bmap [ 0]] = src[nm];
  else             tgt[start[nm]] = src[nm];  
}


void EssentialFunction::describe (char* tgt) const
// ---------------------------------------------------------------------------
// Load descriptive/diagnostic material into tgt.
// ---------------------------------------------------------------------------
{
  sprintf (tgt, "essential-function:\t%s", _function);
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

void Natural::evaluate (const Field*   src    ,
			const int_t    id     , 
			const int_t    plane  , 
			const Element* elmt   ,
			const int_t    side   , 
			const int_t    step   ,
			const bool     Fourier,
			real_t*        tgt    ) const
// ---------------------------------------------------------------------------
// Load external value storage area tgt with installed constant.
// Physical space.
// ---------------------------------------------------------------------------
{
  if (!Fourier) Veclib::fill (Geometry::nP(), _value, tgt, 1);
}


void Natural::sum (const int_t   side  ,
		   const int_t*  bmap  ,
		   const real_t* src   ,
		   const real_t* weight,
		   real_t*       work  ,
		   real_t*       tgt   ) const
// ---------------------------------------------------------------------------
// Add boundary-integral terms (from src) into globally-numbered tgt.
// Work vector must be np long.
// ---------------------------------------------------------------------------
{ 
  const int_t  np    = Geometry::nP();
  const int_t  nm    = np - 1;
  const int_t* start = bmap;
  
  switch (side) {
  case 1: start += nm;           break;
  case 2: start += nm + nm;      break;
  case 3: start += nm + nm + nm; break;
  default: break;
  }

  Veclib::vmul (np, src, 1, weight, 1, work, 1);

  Veclib::scatr_sum (nm, work, start, tgt);
  if   (side == 3) tgt[bmap [ 0]] += work[nm];
  else             tgt[start[nm]] += work[nm];  
}


void Natural::describe (char* tgt) const
// ---------------------------------------------------------------------------
// Load descriptive/diagnostic material into tgt.
// ---------------------------------------------------------------------------
{
  sprintf (tgt, "natural:\t%g", _value);
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

NaturalFunction::NaturalFunction (const char* f)
// ---------------------------------------------------------------------------
// Natural condition that sets value by interpreting function string,
// which can be an explicit function of x, y, z & t as well as any
// installed token values.  Otherwise the same as Natural class.
// ---------------------------------------------------------------------------
{
  strcpy ((_function = new char [strlen (f) + 1]), f);
}


void NaturalFunction::evaluate (const Field*   src    ,
				const int_t    id     , 
				const int_t    plane  , 
				const Element* elmt   ,
				const int_t    side   , 
				const int_t    step   ,
				const bool     Fourier,
				real_t*        tgt    ) const
// ---------------------------------------------------------------------------
// Load external tgt by interpreting function.
// Physical space.
// ---------------------------------------------------------------------------
{
  if (!Fourier) elmt -> sideEval (side, tgt, _function);
}


void NaturalFunction::sum (const int_t   side  ,
			   const int_t*  bmap  ,
			   const real_t* src   ,
			   const real_t* weight,
			   real_t*       work  ,
			   real_t*       tgt   ) const
// ---------------------------------------------------------------------------
// Add boundary-integral terms into globally-numbered tgt.
// ---------------------------------------------------------------------------
{ 
  const int_t  np    = Geometry::nP();
  const int_t  nm    = np - 1;
  const int_t* start = bmap;
  
  switch (side) {
  case 1: start += nm;           break;
  case 2: start += nm + nm;      break;
  case 3: start += nm + nm + nm; break;
  default: break;
  }

  Veclib::vmul (np, src, 1, weight, 1, work, 1);

  Veclib::scatr_sum (nm, work, start, tgt);
  if   (side == 3) tgt[bmap [ 0]] += work[nm];
  else             tgt[start[nm]] += work[nm];  
}


void NaturalFunction::describe (char* tgt) const
// ---------------------------------------------------------------------------
// Load descriptive/diagnostic material into tgt.
// ---------------------------------------------------------------------------
{
  sprintf (tgt, "natural-function:\t%s", _function);
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

Mixed::Mixed (const char* v)
// ---------------------------------------------------------------------------
// The format for a Mixed BC is: "field = mulvalue;refvalue", where in
// dc/dn + K( c - C ) = 0, K is mulvalue and C is refvalue.  The
// separator must be ';' or ',' (and without white space). Each of the
// two supplied values is expected to evaluate to a real_t constant
// (NB: K and C could possibly be functions of time and/or space, but
// are only evaluated *once* at the beginning of execution).
// ---------------------------------------------------------------------------
{
  const char sep[] = ";,";
  char       buf[StrMax], *tok;

  strcpy (buf, v);
  _K_ = Femlib::value (tok = strtok (buf,  sep));
  _C_ = Femlib::value (tok = strtok (NULL, sep));
}


void Mixed::sum (const int_t   side  ,
		 const int_t*  bmap  ,
		 const real_t* src   ,
		 const real_t* weight,
		 real_t*       work  ,
		 real_t*       tgt   ) const
// ---------------------------------------------------------------------------
// Add boundary-integral terms into globally-numbered tgt.  This is
// used to add K*C to RHS forcing.  Input src is ignored.
// ---------------------------------------------------------------------------
{ 
  const int_t  np    = Geometry::nP();
  const int_t  nm    = np - 1;
  const int_t* start = bmap;
  
  switch (side) {
  case 1: start += nm;           break;
  case 2: start += nm + nm;      break;
  case 3: start += nm + nm + nm; break;
  default: break;
  }

  Veclib::smul (np, _K_ * _C_, weight, 1, work, 1);

  Veclib::scatr_sum (nm, work, start, tgt);
  if   (side == 3) tgt[bmap [ 0]] += work[nm];
  else             tgt[start[nm]] += work[nm];
}


void Mixed::augmentSC (const int_t   side  ,
		       const int_t   nband ,
		       const int_t   nsolve,
		       const int_t*  bmap  ,
		       const real_t* area  ,
		       real_t*       work  ,
		       real_t*       tgt   )  const
// ---------------------------------------------------------------------------
// Add <K, w> terms to (banded LAPACK) matrix tgt.
// ---------------------------------------------------------------------------
{
  const int_t  np    = Geometry::nP();
  const int_t  nm    = np - 1;
  const int_t* start = bmap;
  int_t        i, k;
  
  switch (side) {
  case 1: start += nm;           break;
  case 2: start += nm + nm;      break;
  case 3: start += nm + nm + nm; break;
  default: break;
  }

  Veclib::smul (np, _K_, area, 1, work, 1);

  for (i = 0; i < nm; i++)
    if ((k = start[i]) < nsolve)
      tgt[Lapack::band_addr (k, k, nband)] += work[i];

  i = (side == 3) ? bmap[0] : start[nm];
  if (i < nsolve) tgt[Lapack::band_addr (i, i, nband)] += work[nm];
}


void Mixed::augmentOp (const int_t   side, 
		       const int_t*  bmap,
		       const real_t* area,
		       const real_t* src ,
		       real_t*       tgt ) const
// ---------------------------------------------------------------------------
// This operation is used to augment the element-wise Helmholtz
// operations where there are mixed BCs.  Add in diagonal terms
// <K*src, w> to tgt.  Both src and tgt are globally-numbered vectors.
// ---------------------------------------------------------------------------
{
  const int_t  np    = Geometry::nP();
  const int_t  nm    = np - 1;
  const int_t* start = bmap;
  int_t        i;
  
  switch (side) {
  case 1: start += nm;           break;
  case 2: start += nm + nm;      break;
  case 3: start += nm + nm + nm; break;
  default: break;
  }

  for (i = 0; i < nm; i++)
    tgt[start[i]] += _K_ * area[i] * src[start[i]];

  i = (side == 3) ? bmap[0] : start[nm];
  tgt[i] += _K_ * area[nm] * src[i];
}


void Mixed::augmentDg (const int_t   side, 
		       const int_t*  bmap,
		       const real_t* area,
		       real_t*       tgt ) const
// ---------------------------------------------------------------------------
// This operation is used to augment the element-wise construction of
// the diagonal of the global Helmholtz matrix where there are mixed
// BCs.  Add in diagonal terms <K, w> to globally-numbered tgt.
// ---------------------------------------------------------------------------
{
  const int_t  np    = Geometry::nP();
  const int_t  nm    = np - 1;
  const int_t* start = bmap;
  int_t        i;
  
  switch (side) {
  case 1: start += nm;           break;
  case 2: start += nm + nm;      break;
  case 3: start += nm + nm + nm; break;
  default: break;
  }

  for (i = 0; i < nm; i++)
    tgt[start[i]] += _K_ * area[i];

  i = (side == 3) ? bmap[0] : start[nm];
  tgt[i] += _K_ * area[nm];
}


void Mixed::describe (char* tgt) const
// ---------------------------------------------------------------------------
// Load descriptive/diagnostic material into tgt.
// ---------------------------------------------------------------------------
{
  sprintf (tgt, "mixed:\t%g\t%g", _K_, _C_);
}

//////////////////////////////////////////////////////////////////////////////
// Internally computed BC types follow.
//
// These are all evaluated in Fourier space.
//////////////////////////////////////////////////////////////////////////////

void NaturalCBCp::evaluate (const Field*   src    ,
			    const int_t    id     , 
			    const int_t    plane  , 
			    const Element* elmt   ,
			    const int_t    side   , 
			    const int_t    step   ,
			    const bool     Fourier,
			    real_t*        tgt    ) const
// ---------------------------------------------------------------------------
// Load external tgt via a call to BCmgr to compute terms.
// ---------------------------------------------------------------------------
{
  if (Fourier) _BCmgr -> evaluateCNBCp (id, plane, step, tgt); 
}


void NaturalCBCp::sum (const int_t   side  ,
		       const int_t*  bmap  ,
		       const real_t* src   ,
		       const real_t* weight,
		       real_t*       work  ,
		       real_t*       tgt   ) const
// ---------------------------------------------------------------------------
// Add boundary-integral src terms into globally-numbered tgt.
// ---------------------------------------------------------------------------
{ 
  const int_t  np    = Geometry::nP();
  const int_t  nm    = np - 1;
  const int_t* start = bmap;
  
  switch (side) {
  case 1: start += nm;           break;
  case 2: start += nm + nm;      break;
  case 3: start += nm + nm + nm; break;
  default: break;
  }

  Veclib::vmul (np, src, 1, weight, 1, work, 1);

  Veclib::scatr_sum (nm, work, start, tgt);
  if   (side == 3) tgt[bmap [ 0]] += work[nm];
  else             tgt[start[nm]] += work[nm];  
}


void NaturalCBCp::describe (char* tgt) const
// ---------------------------------------------------------------------------
// Load descriptive/diagnostic material into tgt.
// ---------------------------------------------------------------------------
{
  sprintf (tgt, "computed-natural-pressure, see KIO91");
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

MixedCBCp::MixedCBCp (BCmgr* B) : _BCmgr (B)
// ---------------------------------------------------------------------------
// Computed mixed/Robin pressure BC for open boundaries, Dong (2015) eq. (37).
// For the LHS, our notation is dp/dn + _K_p, where _K_ is set as below:
// ---------------------------------------------------------------------------
{
  _K_ = 1.0 / Femlib::value ("KINVIS*DONG_DO");
}


void MixedCBCp::evaluate (const Field*   src    ,
			  const int_t    id     , 
			  const int_t    plane  , 
			  const Element* elmt   ,
			  const int_t    side   , 
			  const int_t    step   ,
			  const bool     Fourier,
			  real_t*        tgt    ) const
// ---------------------------------------------------------------------------
// Load external tgt via a call to BCmgr to compute terms.
// ---------------------------------------------------------------------------
{
  if (Fourier) _BCmgr -> evaluateCMBCp (src, id, plane, step, tgt); 
}


void MixedCBCp::sum (const int_t   side  ,
		     const int_t*  bmap  ,
		     const real_t* src   ,
		     const real_t* weight,
		     real_t*       work  ,
		     real_t*       tgt   ) const
// ---------------------------------------------------------------------------
// Add boundary-integral src terms into globally-numbered tgt, the RHS
// of Dong (2015) eq. (37).  Input src (that was created by evaluate)
// takes the place of K*C.
// ---------------------------------------------------------------------------
{
  const int_t  np    = Geometry::nP();
  const int_t  nm    = np - 1;
  const int_t* start = bmap;
  
  switch (side) {
  case 1: start += nm;           break;
  case 2: start += nm + nm;      break;
  case 3: start += nm + nm + nm; break;
  default: break;
  }

  Veclib::vmul (np, src, 1, weight, 1, work, 1);

  Veclib::scatr (nm, work, start, tgt);
  if   (side == 3) tgt[bmap [ 0]] = work[nm];
  else             tgt[start[nm]] = work[nm];  
}


void MixedCBCp::augmentSC (const int_t   side  ,
			   const int_t   nband ,
			   const int_t   nsolve,
			   const int_t*  bmap  ,
			   const real_t* area  ,
			   real_t*       work  ,
			   real_t*       tgt   )  const
// ---------------------------------------------------------------------------
// Add <K, w> terms to (banded LAPACK) matrix tgt.
// ---------------------------------------------------------------------------
{
  const int_t  np    = Geometry::nP();
  const int_t  nm    = np - 1;
  const int_t* start = bmap;
  int_t        i, k;
  
  switch (side) {
  case 1: start += nm;           break;
  case 2: start += nm + nm;      break;
  case 3: start += nm + nm + nm; break;
  default: break;
  }

  Veclib::smul (np, _K_, area, 1, work, 1);

  for (i = 0; i < nm; i++)
    if ((k = start[i]) < nsolve)
      tgt[Lapack::band_addr (k, k, nband)] += work[i];

  i = (side == 3) ? bmap[0] : start[nm];
  if (i < nsolve) tgt[Lapack::band_addr (i, i, nband)] += work[nm];
}


void MixedCBCp::augmentOp (const int_t   side, 
			   const int_t*  bmap,
			   const real_t* area,
			   const real_t* src ,
			   real_t*       tgt ) const
// ---------------------------------------------------------------------------
// This operation is used to augment the element-wise Helmholtz
// operations where there are mixed BCs.  Add in diagonal terms
// <K*src, w> to tgt.  Both src and tgt are globally-numbered vectors.
// ---------------------------------------------------------------------------
{
  const int_t  np    = Geometry::nP();
  const int_t  nm    = np - 1;
  const int_t* start = bmap;
  int_t        i;
  
  switch (side) {
  case 1: start += nm;           break;
  case 2: start += nm + nm;      break;
  case 3: start += nm + nm + nm; break;
  default: break;
  }

  for (i = 0; i < nm; i++)
    tgt[start[i]] += _K_ * area[i] * src[start[i]];

  i = (side == 3) ? bmap[0] : start[nm];
  tgt[i] += _K_ * area[nm] * src[i];
}


void MixedCBCp::augmentDg (const int_t   side, 
			   const int_t*  bmap,
			   const real_t* area,
			   real_t*       tgt ) const
// ---------------------------------------------------------------------------
// This operation is used to augment the element-wise construction of
// the diagonal of the global Helmholtz matrix where there are mixed
// BCs.  Add in diagonal terms <K, w> to globally-numbered tgt.
// ---------------------------------------------------------------------------
{
  const int_t  np    = Geometry::nP();
  const int_t  nm    = np - 1;
  const int_t* start = bmap;
  int_t        i;
  
  switch (side) {
  case 1: start += nm;           break;
  case 2: start += nm + nm;      break;
  case 3: start += nm + nm + nm; break;
  default: break;
  }

  for (i = 0; i < nm; i++)
    tgt[start[i]] += _K_ * area[i];

  i = (side == 3) ? bmap[0] : start[nm];
  tgt[i] += _K_ * area[nm];
}


void MixedCBCp::describe (char* tgt) const
// ---------------------------------------------------------------------------
// Load descriptive/diagnostic material into tgt.
// ---------------------------------------------------------------------------
{
  sprintf (tgt, "computed-mixed-pressure, Dong15 eq. (37)");
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

MixedCBCu::MixedCBCu (BCmgr* B) : _BCmgr (B)
// ---------------------------------------------------------------------------
// Computed mixed/Robin velocity BC for open boundaries, Dong (2015) eq. (38).
// For the LHS, our notation is du/dn + _K_u, where _K_ is set as below.
//
// NOTE: there is a potential problem here, since alpha[0] (and hence
// _K_) depends on time-stepping order, so it will likely be incorrect
// for the first time-step (or two).  However, it is fine for setting up 
// matrices for the asymptotic timestep order.  For all other cases we will
// evaluate _K_ based on the current time step.				  
// ---------------------------------------------------------------------------
{
  _alpha = new real_t [Integration::OrderMax];
  Integration::StifflyStable (Femlib::ivalue ("N_TIME"), _alpha);

  _DoDt = Femlib::value ("DONG_DO/D_T");
  _K_   = _alpha[0] * _DoDt;
}


void MixedCBCu::evaluate (const Field*   src    ,
			  const int_t    id     , 
			  const int_t    plane  , 
			  const Element* elmt   ,
			  const int_t    side   , 
			  const int_t    step   ,
			  const bool     Fourier,
			  real_t*        tgt    ) const
// ---------------------------------------------------------------------------
// Load external tgt via a call to BCmgr to compute terms.
// ---------------------------------------------------------------------------
{
  if (Fourier) _BCmgr -> evaluateCMBCu (src, id, plane, step, 'u', tgt); 
}


void MixedCBCu::sum (const int_t   side  ,
		     const int_t*  bmap  ,
		     const real_t* src   ,
		     const real_t* weight,
		     real_t*       work  ,
		     real_t*       tgt   ) const
// ---------------------------------------------------------------------------
// Add boundary-integral src terms into globally-numbered tgt.
// ---------------------------------------------------------------------------
{ 
  const int_t  np    = Geometry::nP();
  const int_t  nm    = np - 1;
  const int_t* start = bmap;
  
  switch (side) {
  case 1: start += nm;           break;
  case 2: start += nm + nm;      break;
  case 3: start += nm + nm + nm; break;
  default: break;
  }

  Veclib::vmul (np, src, 1, weight, 1, work, 1);

  Veclib::scatr_sum (nm, work, start, tgt);
  if   (side == 3) tgt[bmap [ 0]] += work[nm];
  else             tgt[start[nm]] += work[nm];  
}


void MixedCBCu::augmentSC (const int_t   side  ,
			   const int_t   nband ,
			   const int_t   nsolve,
			   const int_t*  bmap  ,
			   const real_t* area  ,
			   real_t*       work  ,
			   real_t*       tgt   )  const
// ---------------------------------------------------------------------------
// Add <K, w> terms to (banded LAPACK) matrix tgt.  This call only
// occurs during construction of direct-solve matrices, prior to
// time-stepping.  Hence it uses the value of _K_ set during call to
// MixedCBCu constructor, which used the alpha vector for the eventual
// steady-state timestepping order.
// ---------------------------------------------------------------------------
{
  const int_t  np    = Geometry::nP();
  const int_t  nm    = np - 1;
  const int_t* start = bmap;
  int_t        i, k;
  
  switch (side) {
  case 1: start += nm;           break;
  case 2: start += nm + nm;      break;
  case 3: start += nm + nm + nm; break;
  default: break;
  }

  Veclib::smul (np, _K_, area, 1, work, 1);

  for (i = 0; i < nm; i++)
    if ((k = start[i]) < nsolve)
      tgt[Lapack::band_addr (k, k, nband)] += work[i];

  i = (side == 3) ? bmap[0] : start[nm];
  if (i < nsolve) tgt[Lapack::band_addr (i, i, nband)] += work[nm];
}


void MixedCBCu::augmentOp (const int_t   side, 
			   const int_t*  bmap,
			   const real_t* area,
			   const real_t* src ,
			   real_t*       tgt ) const
// ---------------------------------------------------------------------------
// This operation is used to augment the element-wise Helmholtz
// operations where there are mixed BCs.  Add in diagonal terms
// <K*src, w> to tgt.  Both src and tgt are globally-numbered vectors.
// ---------------------------------------------------------------------------
{
  const int_t   np    = Geometry::nP();
  const int_t   nm    = np - 1;
  const int_t*  start = bmap;
  int_t         i;

  i = clamp (Femlib::ivalue("STEP"), 1, Femlib::ivalue("N_TIME"));

  Integration::StifflyStable (i, _alpha);
  const real_t Kloc = _alpha[0] * _DoDt;
  
  switch (side) {
  case 1: start += nm;           break;
  case 2: start += nm + nm;      break;
  case 3: start += nm + nm + nm; break;
  default: break;
  }

  for (i = 0; i < nm; i++)
    tgt[start[i]] += Kloc * area[i] * src[start[i]];

  i = (side == 3) ? bmap[0] : start[nm];
  tgt[i] += Kloc * area[nm] * src[i];
}


void MixedCBCu::augmentDg (const int_t   side, 
			   const int_t*  bmap,
			   const real_t* area,
			   real_t*       tgt ) const
// ---------------------------------------------------------------------------
// This operation is used to augment the element-wise construction of
// the diagonal of the global Helmholtz matrix where there are mixed
// BCs.  Add in diagonal terms <K, w> to globally-numbered tgt.
// ---------------------------------------------------------------------------
{
  const int_t  np    = Geometry::nP();
  const int_t  nm    = np - 1;
  const int_t* start = bmap;
  int_t        i;

  i = clamp (Femlib::ivalue("STEP"), 1, Femlib::ivalue("N_TIME"));
  
  Integration::StifflyStable (i, _alpha);
  const real_t Kloc = _alpha[0] * _DoDt;
  
  switch (side) {
  case 1: start += nm;           break;
  case 2: start += nm + nm;      break;
  case 3: start += nm + nm + nm; break;
  default: break;
  }

  for (i = 0; i < nm; i++)
    tgt[start[i]] += Kloc * area[i];

  i = (side == 3) ? bmap[0] : start[nm];
  tgt[i] += Kloc * area[nm];
}


void MixedCBCu::describe (char* tgt) const
// ---------------------------------------------------------------------------
// Load descriptive/diagnostic material into tgt.
// ---------------------------------------------------------------------------
{
  sprintf (tgt, "computed-mixed-velocity (u cmpt.), Dong15 eq. (38)");
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

MixedCBCv::MixedCBCv (BCmgr* B) : _BCmgr (B)
// ---------------------------------------------------------------------------
// Computed mixed/Robin velocity BC for open boundaries, Dong (2015) eq. (38).
// For the LHS, our notation is du/dn + _K_u, where _K_ is set as below.
//
// NOTE: there is a potential problem here, since alpha[0] (and hence
// _K_) depends on time-stepping order, so it will likely be incorrect
// for the first time-step (or two).  However, it is fine for setting up 
// matrices for the asymptotic timestep order.  For all other cases we will
// evaluate _K_ based on the current time step.				  
// ---------------------------------------------------------------------------
{
  _alpha = new real_t [Integration::OrderMax];
  Integration::StifflyStable (Femlib::ivalue ("N_TIME"), _alpha);

  _DoDt = Femlib::value ("DONG_DO/D_T");
  _K_   = _alpha[0] * _DoDt;
}


void MixedCBCv::evaluate (const Field*   src    ,
			  const int_t    id     , 
			  const int_t    plane  , 
			  const Element* elmt   ,
			  const int_t    side   , 
			  const int_t    step   ,
			  const bool     Fourier,
			  real_t*        tgt    ) const
// ---------------------------------------------------------------------------
// Load external tgt via a call to BCmgr to compute terms.
// ---------------------------------------------------------------------------
{
  if (Fourier) _BCmgr -> evaluateCMBCu (src, id, plane, step, 'v', tgt); 
}


void MixedCBCv::sum (const int_t   side  ,
		     const int_t*  bmap  ,
		     const real_t* src   ,
		     const real_t* weight,
		     real_t*       work  ,
		     real_t*       tgt   ) const
// ---------------------------------------------------------------------------
// Add boundary-integral src terms into globally-numbered tgt.
// ---------------------------------------------------------------------------
{ 
  const int_t  np    = Geometry::nP();
  const int_t  nm    = np - 1;
  const int_t* start = bmap;
  
  switch (side) {
  case 1: start += nm;           break;
  case 2: start += nm + nm;      break;
  case 3: start += nm + nm + nm; break;
  default: break;
  }

  Veclib::vmul (np, src, 1, weight, 1, work, 1);

  Veclib::scatr_sum (nm, work, start, tgt);
  if   (side == 3) tgt[bmap [ 0]] += work[nm];
  else             tgt[start[nm]] += work[nm];  
}


void MixedCBCv::augmentSC (const int_t   side  ,
			   const int_t   nband ,
			   const int_t   nsolve,
			   const int_t*  bmap  ,
			   const real_t* area  ,
			   real_t*       work  ,
			   real_t*       tgt   )  const
// ---------------------------------------------------------------------------
// Add <K, w> terms to (banded LAPACK) matrix tgt.
// ---------------------------------------------------------------------------
{
  const int_t  np    = Geometry::nP();
  const int_t  nm    = np - 1;
  const int_t* start = bmap;
  int_t        i, k;
  
  switch (side) {
  case 1: start += nm;           break;
  case 2: start += nm + nm;      break;
  case 3: start += nm + nm + nm; break;
  default: break;
  }

  Veclib::smul (np, _K_, area, 1, work, 1);

  for (i = 0; i < nm; i++)
    if ((k = start[i]) < nsolve)
      tgt[Lapack::band_addr (k, k, nband)] += work[i];

  i = (side == 3) ? bmap[0] : start[nm];
  if (i < nsolve) tgt[Lapack::band_addr (i, i, nband)] += work[nm];
}


void MixedCBCv::augmentOp (const int_t   side, 
			   const int_t*  bmap,
			   const real_t* area,
			   const real_t* src ,
			   real_t*       tgt ) const
// ---------------------------------------------------------------------------
// This operation is used to augment the element-wise Helmholtz
// operations where there are mixed BCs.  Add in diagonal terms
// <K*src, w> to tgt.  Both src and tgt are globally-numbered vectors.
// ---------------------------------------------------------------------------
{
  const int_t   np    = Geometry::nP();
  const int_t   nm    = np - 1;
  const int_t*  start = bmap;
  int_t         i;

  i = clamp (Femlib::ivalue("STEP"), 1, Femlib::ivalue("N_TIME"));

  Integration::StifflyStable (i, _alpha);
  const real_t Kloc = _alpha[0] * _DoDt;
  
  switch (side) {
  case 1: start += nm;           break;
  case 2: start += nm + nm;      break;
  case 3: start += nm + nm + nm; break;
  default: break;
  }

  for (i = 0; i < nm; i++)
    tgt[start[i]] += Kloc * area[i] * src[start[i]];

  i = (side == 3) ? bmap[0] : start[nm];
  tgt[i] += Kloc * area[nm] * src[i];
}


void MixedCBCv::augmentDg (const int_t   side, 
			   const int_t*  bmap,
			   const real_t* area,
			   real_t*       tgt ) const
// ---------------------------------------------------------------------------
// This operation is used to augment the element-wise construction of
// the diagonal of the global Helmholtz matrix where there are mixed
// BCs.  Add in diagonal terms <K, w> to globally-numbered tgt.
// ---------------------------------------------------------------------------
{
  const int_t  np    = Geometry::nP();
  const int_t  nm    = np - 1;
  const int_t* start = bmap;
  int_t        i;

  i = clamp (Femlib::ivalue("STEP"), 1, Femlib::ivalue("N_TIME"));

  Integration::StifflyStable (i, _alpha);
  const real_t Kloc = _alpha[0] * _DoDt;
  
  switch (side) {
  case 1: start += nm;           break;
  case 2: start += nm + nm;      break;
  case 3: start += nm + nm + nm; break;
  default: break;
  }

  for (i = 0; i < nm; i++)
    tgt[start[i]] += Kloc * area[i];

  i = (side == 3) ? bmap[0] : start[nm];
  tgt[i] += Kloc * area[nm];
}


void MixedCBCv::describe (char* tgt) const
// ---------------------------------------------------------------------------
// Load descriptive/diagnostic material into tgt.
// ---------------------------------------------------------------------------
{
  sprintf (tgt, "computed-mixed-velocity (v cmpt.), Dong15 eq. (38)");
}
 
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

MixedCBCw::MixedCBCw (BCmgr* B) : _BCmgr (B)
// ---------------------------------------------------------------------------
// Computed mixed/Robin velocity BC for open boundaries, Dong (2015) eq. (38).
// For the LHS, our notation is du/dn + _K_u, where _K_ is set as below.
//
// NOTE: there is a potential problem here, since alpha[0] (and hence
// _K_) depends on time-stepping order, so it will likely be incorrect
// for the first time-step (or two).  However, it is fine for setting up 
// matrices for the asymptotic timestep order.  For all other cases we will
// evaluate _K_ based on the current time step.				  
// ---------------------------------------------------------------------------
{
  _alpha = new real_t [Integration::OrderMax];
  Integration::StifflyStable (Femlib::ivalue ("N_TIME"), _alpha);

  _DoDt = Femlib::value ("DONG_DO/D_T");
  _K_   = _alpha[0] * _DoDt;
}


void MixedCBCw::evaluate (const Field*   src    ,
			  const int_t    id     , 
			  const int_t    plane  , 
			  const Element* elmt   ,
			  const int_t    side   , 
			  const int_t    step   ,
			  const bool     Fourier,
			  real_t*        tgt    ) const
// ---------------------------------------------------------------------------
// Load external tgt via a call to BCmgr to compute terms.
// ---------------------------------------------------------------------------
{
  if (Fourier) _BCmgr -> evaluateCMBCu (src, id, plane, step, 'w', tgt); 
}


void MixedCBCw::sum (const int_t   side  ,
		     const int_t*  bmap  ,
		     const real_t* src   ,
		     const real_t* weight,
		     real_t*       work  ,
		     real_t*       tgt   ) const
// ---------------------------------------------------------------------------
// Add boundary-integral src terms into globally-numbered tgt.
// ---------------------------------------------------------------------------
{ 
  const int_t  np    = Geometry::nP();
  const int_t  nm    = np - 1;
  const int_t* start = bmap;
  
  switch (side) {
  case 1: start += nm;           break;
  case 2: start += nm + nm;      break;
  case 3: start += nm + nm + nm; break;
  default: break;
  }

  Veclib::vmul (np, src, 1, weight, 1, work, 1);

  Veclib::scatr_sum (nm, work, start, tgt);
  if   (side == 3) tgt[bmap [ 0]] += work[nm];
  else             tgt[start[nm]] += work[nm];  
}


void MixedCBCw::augmentSC (const int_t   side  ,
			   const int_t   nband ,
			   const int_t   nsolve,
			   const int_t*  bmap  ,
			   const real_t* area  ,
			   real_t*       work  ,
			   real_t*       tgt   )  const
// ---------------------------------------------------------------------------
// Add <K, w> terms to (banded LAPACK) matrix tgt.
// ---------------------------------------------------------------------------
{
  const int_t  np    = Geometry::nP();
  const int_t  nm    = np - 1;
  const int_t* start = bmap;
  int_t        i, k;
  
  switch (side) {
  case 1: start += nm;           break;
  case 2: start += nm + nm;      break;
  case 3: start += nm + nm + nm; break;
  default: break;
  }

  Veclib::smul (np, _K_, area, 1, work, 1);

  for (i = 0; i < nm; i++)
    if ((k = start[i]) < nsolve)
      tgt[Lapack::band_addr (k, k, nband)] += work[i];

  i = (side == 3) ? bmap[0] : start[nm];
  if (i < nsolve) tgt[Lapack::band_addr (i, i, nband)] += work[nm];
}


void MixedCBCw::augmentOp (const int_t   side, 
			   const int_t*  bmap,
			   const real_t* area,
			   const real_t* src ,
			   real_t*       tgt ) const
// ---------------------------------------------------------------------------
// This operation is used to augment the element-wise Helmholtz
// operations where there are mixed BCs.  Add in diagonal terms
// <K*src, w> to tgt.  Both src and tgt are globally-numbered vectors.
// ---------------------------------------------------------------------------
{
  const int_t  np    = Geometry::nP();
  const int_t  nm    = np - 1;
  const int_t* start = bmap;
  int_t        i;

  i = clamp (Femlib::ivalue("STEP"), 1, Femlib::ivalue("N_TIME"));

  Integration::StifflyStable (i, _alpha);
  const real_t Kloc = _alpha[0] * _DoDt;
  
  switch (side) {
  case 1: start += nm;           break;
  case 2: start += nm + nm;      break;
  case 3: start += nm + nm + nm; break;
  default: break;
  }

  for (i = 0; i < nm; i++)
    tgt[start[i]] += Kloc * area[i] * src[start[i]];

  i = (side == 3) ? bmap[0] : start[nm];
  tgt[i] += Kloc * area[nm] * src[i];
}


void MixedCBCw::augmentDg (const int_t   side, 
			   const int_t*  bmap,
			   const real_t* area,
			   real_t*       tgt ) const
// ---------------------------------------------------------------------------
// This operation is used to augment the element-wise construction of
// the diagonal of the global Helmholtz matrix where there are mixed
// BCs.  Add in diagonal terms <K, w> to globally-numbered tgt.
// ---------------------------------------------------------------------------
{
  const int_t  np    = Geometry::nP();
  const int_t  nm    = np - 1;
  const int_t* start = bmap;
  int_t        i;

  i = clamp (Femlib::ivalue("STEP"), 1, Femlib::ivalue("N_TIME"));

  Integration::StifflyStable (i, _alpha);
  const real_t Kloc = _alpha[0] * _DoDt;
  
  switch (side) {
  case 1: start += nm;           break;
  case 2: start += nm + nm;      break;
  case 3: start += nm + nm + nm; break;
  default: break;
  }

  for (i = 0; i < nm; i++)
    tgt[start[i]] += Kloc * area[i];

  i = (side == 3) ? bmap[0] : start[nm];
  tgt[i] += Kloc * area[nm];
}


void MixedCBCw::describe (char* tgt) const
// ---------------------------------------------------------------------------
// Load descriptive/diagnostic material into tgt.
// ---------------------------------------------------------------------------
{
  sprintf (tgt, "computed-mixed-velocity (w cmpt.), Dong15 eq. (38)");
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

MixedCBCwIn::MixedCBCwIn (BCmgr* B) : _BCmgr (B)
// ---------------------------------------------------------------------------
// Computed mixed/Robin velocity BC for open boundaries, Dong (2015) eq. (38).
// For the LHS, our notation is du/dn + _K_u, where _K_ is set as below.
//
// This is a cut-down version of MixedCBCw, which does nothing for the
// boundary integral terms.  In effect this means that _C_ = 0 and we have
// 
//   dw
//   --  +  K w = 0.
//   dn 				      
// ---------------------------------------------------------------------------
{
  _alpha = new real_t [Integration::OrderMax];
  Integration::StifflyStable (Femlib::ivalue ("N_TIME"), _alpha);

  _DoDt = Femlib::value ("DONG_DO/D_T");
  _K_   = _alpha[0] * _DoDt;
}


void MixedCBCwIn::augmentSC (const int_t   side  ,
			     const int_t   nband ,
			     const int_t   nsolve,
			     const int_t*  bmap  ,
			     const real_t* area  ,
			     real_t*       work  ,
			     real_t*       tgt   )  const
// ---------------------------------------------------------------------------
// Add <K, w> terms to (banded LAPACK) matrix tgt.
// ---------------------------------------------------------------------------
{
  const int_t  np    = Geometry::nP();
  const int_t  nm    = np - 1;
  const int_t* start = bmap;
  int_t        i, k;
  
  switch (side) {
  case 1: start += nm;           break;
  case 2: start += nm + nm;      break;
  case 3: start += nm + nm + nm; break;
  default: break;
  }

  Veclib::smul (np, _K_, area, 1, work, 1);

  for (i = 0; i < nm; i++)
    if ((k = start[i]) < nsolve)
      tgt[Lapack::band_addr (k, k, nband)] += work[i];

  i = (side == 3) ? bmap[0] : start[nm];
  if (i < nsolve) tgt[Lapack::band_addr (i, i, nband)] += work[nm];
}


void MixedCBCwIn::augmentOp (const int_t   side, 
			     const int_t*  bmap,
			     const real_t* area,
			     const real_t* src ,
			     real_t*       tgt ) const
// ---------------------------------------------------------------------------
// This operation is used to augment the element-wise Helmholtz
// operations where there are mixed BCs.  Add in diagonal terms
// <K*src, w> to tgt.  Both src and tgt are globally-numbered vectors.
// ---------------------------------------------------------------------------
{
  const int_t   np    = Geometry::nP();
  const int_t   nm    = np - 1;
  const int_t*  start = bmap;
  int_t         i;

  i = clamp (Femlib::ivalue("STEP"), 1, Femlib::ivalue("N_TIME"));

  Integration::StifflyStable (i, _alpha);
  const real_t Kloc = _alpha[0] * _DoDt;
  
  switch (side) {
  case 1: start += nm;           break;
  case 2: start += nm + nm;      break;
  case 3: start += nm + nm + nm; break;
  default: break;
  }

  for (i = 0; i < nm; i++)
    tgt[start[i]] += Kloc * area[i] * src[start[i]];

  i = (side == 3) ? bmap[0] : start[nm];
  tgt[i] += Kloc * area[nm] * src[i];
}


void MixedCBCwIn::augmentDg (const int_t   side, 
			     const int_t*  bmap,
			     const real_t* area,
			     real_t*       tgt ) const
// ---------------------------------------------------------------------------
// This operation is used to augment the element-wise construction of
// the diagonal of the global Helmholtz matrix where there are mixed
// BCs.  Add in diagonal terms <K, w> to globally-numbered tgt.
// ---------------------------------------------------------------------------
{
  const int_t  np    = Geometry::nP();
  const int_t  nm    = np - 1;
  const int_t* start = bmap;
  int_t        i;

  i = clamp (Femlib::ivalue("STEP"), 1, Femlib::ivalue("N_TIME"));

  Integration::StifflyStable (i, _alpha);
  const real_t Kloc = _alpha[0] * _DoDt;
  
  switch (side) {
  case 1: start += nm;           break;
  case 2: start += nm + nm;      break;
  case 3: start += nm + nm + nm; break;
  default: break;
  }

  for (i = 0; i < nm; i++)
    tgt[start[i]] += Kloc * area[i];

  i = (side == 3) ? bmap[0] : start[nm];
  tgt[i] += Kloc * area[nm];
}


void MixedCBCwIn::describe (char* tgt) const
// ---------------------------------------------------------------------------
// Load descriptive/diagnostic material into tgt.
// ---------------------------------------------------------------------------
{
  sprintf (tgt, "computed-mixed-velocity (w cmpt.), Dong15 eq. (38), but C=0");
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

MixedCBCc::MixedCBCc (BCmgr* B) : _BCmgr (B)
// ---------------------------------------------------------------------------
// Computed mixed/Robin velocity BC for open boundaries, Dong (2015) eq. (38).
// For the LHS, our notation is du/dn + _K_u, where _K_ is set as below.
//
// NOTE: there is a potential problem here, since alpha[0] (and hence
// _K_) depends on time-stepping order, so it will likely be incorrect
// for the first time-step (or two).  However, it is fine for setting up 
// matrices for the asymptotic timestep order.  For all other cases we will
// evaluate _K_ based on the current time step.				  
// ---------------------------------------------------------------------------
{
  _alpha = new real_t [Integration::OrderMax];
  Integration::StifflyStable (Femlib::ivalue ("N_TIME"), _alpha);

  _DoDt = Femlib::value ("DONG_DO/D_T");
  _K_   = _alpha[0] * _DoDt;
}


void MixedCBCc::evaluate (const Field*   src    ,
			  const int_t    id     , 
			  const int_t    plane  , 
			  const Element* elmt   ,
			  const int_t    side   , 
			  const int_t    step   ,
			  const bool     Fourier,
			  real_t*        tgt    ) const
// ---------------------------------------------------------------------------
// Load external tgt via a call to BCmgr to compute terms equivalent
// to _K_*_C_ in a regular Mixed BC.
// ---------------------------------------------------------------------------
{
  if (Fourier) _BCmgr -> evaluateCMBCc (id, plane, step, tgt);
}


void MixedCBCc::sum (const int_t   side  ,
		     const int_t*  bmap  ,
		     const real_t* src   ,
		     const real_t* weight,
		     real_t*       work  ,
		     real_t*       tgt   ) const
// ---------------------------------------------------------------------------
// Add boundary-integral src terms (computed by evaluate into *its*
// tgt, delivered here as src) into globally-numbered tgt.
// ---------------------------------------------------------------------------
{ 
  const int_t  np    = Geometry::nP();
  const int_t  nm    = np - 1;
  const int_t* start = bmap;
  
  switch (side) {
  case 1: start += nm;           break;
  case 2: start += nm + nm;      break;
  case 3: start += nm + nm + nm; break;
  default: break;
  }

  Veclib::vmul (np, src, 1, weight, 1, work, 1);

  Veclib::scatr_sum (nm, work, start, tgt);
  if   (side == 3) tgt[bmap [ 0]] += work[nm];
  else             tgt[start[nm]] += work[nm];  
}


void MixedCBCc::augmentSC (const int_t   side  ,
			   const int_t   nband ,
			   const int_t   nsolve,
			   const int_t*  bmap  ,
			   const real_t* area  ,
			   real_t*       work  ,
			   real_t*       tgt   )  const
// ---------------------------------------------------------------------------
// Add <K, w> terms to (banded LAPACK) matrix tgt.
// ---------------------------------------------------------------------------
{
  const int_t  np    = Geometry::nP();
  const int_t  nm    = np - 1;
  const int_t* start = bmap;
  int_t        i, k;
  
  switch (side) {
  case 1: start += nm;           break;
  case 2: start += nm + nm;      break;
  case 3: start += nm + nm + nm; break;
  default: break;
  }

  Veclib::smul (np, _K_, area, 1, work, 1);

  for (i = 0; i < nm; i++)
    if ((k = start[i]) < nsolve)
      tgt[Lapack::band_addr (k, k, nband)] += work[i];

  i = (side == 3) ? bmap[0] : start[nm];
  if (i < nsolve) tgt[Lapack::band_addr (i, i, nband)] += work[nm];
}


void MixedCBCc::augmentOp (const int_t   side, 
			   const int_t*  bmap,
			   const real_t* area,
			   const real_t* src ,
			   real_t*       tgt ) const
// ---------------------------------------------------------------------------
// This operation is used to augment the element-wise Helmholtz
// operations where there are mixed BCs.  Add in diagonal terms
// <K*src, w> to tgt.  Both src and tgt are globally-numbered vectors.
// ---------------------------------------------------------------------------
{
  const int_t   np    = Geometry::nP();
  const int_t   nm    = np - 1;
  const int_t*  start = bmap;
  int_t         i;

  i = clamp (Femlib::ivalue("STEP"), 1, Femlib::ivalue("N_TIME"));

  Integration::StifflyStable (i, _alpha);
  const real_t Kloc = _alpha[0] * _DoDt;
  
  switch (side) {
  case 1: start += nm;           break;
  case 2: start += nm + nm;      break;
  case 3: start += nm + nm + nm; break;
  default: break;
  }

  for (i = 0; i < nm; i++)
    tgt[start[i]] += Kloc * area[i] * src[start[i]];

  i = (side == 3) ? bmap[0] : start[nm];
  tgt[i] += Kloc * area[nm] * src[i];
}


void MixedCBCc::augmentDg (const int_t   side, 
			   const int_t*  bmap,
			   const real_t* area,
			   real_t*       tgt ) const
// ---------------------------------------------------------------------------
// This operation is used to augment the element-wise construction of
// the diagonal of the global Helmholtz matrix where there are mixed
// BCs.  Add in diagonal terms <K, w> to globally-numbered tgt.
// ---------------------------------------------------------------------------
{
  const int_t  np    = Geometry::nP();
  const int_t  nm    = np - 1;
  const int_t* start = bmap;
  int_t        i;

  i = clamp (Femlib::ivalue("STEP"), 1, Femlib::ivalue("N_TIME"));

  Integration::StifflyStable (i, _alpha);
  const real_t Kloc = _alpha[0] * _DoDt;
  
  switch (side) {
  case 1: start += nm;           break;
  case 2: start += nm + nm;      break;
  case 3: start += nm + nm + nm; break;
  default: break;
  }

  for (i = 0; i < nm; i++)
    tgt[start[i]] += Kloc * area[i];

  i = (side == 3) ? bmap[0] : start[nm];
  tgt[i] += Kloc * area[nm];
}


void MixedCBCc::describe (char* tgt) const
// ---------------------------------------------------------------------------
// Load descriptive/diagnostic material into tgt.
// ---------------------------------------------------------------------------
{
  sprintf (tgt, "computed-mixed-scalar, Liu, Xie & Dong 2020 eq. (16b)");
}
