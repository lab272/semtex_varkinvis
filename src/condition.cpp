//////////////////////////////////////////////////////////////////////////////
// condition.cpp: functions used to evaluate & apply BCs along element edges.
//
// Copyright (c) 1994+, Hugh M Blackburn
//
// Conditions have two abstract layers, and a terminal/concrete layer.
// This is to allow all boundary conditions to routines to be
// uniformly applied, use dynamic polymorphism to invoke the correct
// termimal method in each case, and minimize code-re-use.
//
// All concrete classes inherit the abstract base class Condition, as
// well as one of the abstract base classes Essential, Natural, or
// Mixed.  Owing to the different behaviour of essential, natural and
// mixed BCs, some virtual routines may take no action in the
// different cases; due also to the generality of the Condition class
// methods, many routines do not use all parameters.
//
// The two most basic and commonly used kinds of conditions are
// essential (Dirichlet) and natural (Neumann).  Essential conditions
// are set/imposed directly, while natural conditions are applied by
// summation of integral terms, owing to the fact that they derive
// from integration by parts.  So function "set" is only used by
// essential-type conditions while function "sum" is only used by
// natural-type conditions.
//
// The HOPBC condition are a special case of natural conditions used
// for the pressure Poisson equation: they are derived from the
// momentum equations and computed using an extrapolative process
// described in [1].
//
// A further general class of available condition is the mixed or
// Robin type, of form dc/dn + K( c - C ) = 0, equiv. dc/dn + Kc = KC.
// For this type, terms on the diagonal of the associated global
// Helmholtz matrix are augmented; the routines "augment**" implement
// these actions during either construction or application of the
// Helmholtz matrix.  See [2].
//
// Boundary condition methods recognise the fact that the actual
// numerical values for each solution Field are held by the sheet/line
// (3D/2D) storage for that class; that storage "wraps around" the
// periphery of the solution Domain.  The "evaluate" methods are used
// to store appropriate values into that Field storage, while the
// "set", "sum", "augment**" methods are used during solution (and, in
// the case of mixed conditions, construction of) global Helmholtz
// problems; these latter methods retrieve their pre-evaluated data
// from Field storage and deal with it as required for the BC type.
//
// The "set" methods are used for essential BCs, and just impose the
// appropriate values into global RHS storage for values that are
// lifted out of the Helmholtz solution, while the "sum" methods are
// used for natural (and mixed) BCs where the numerical values are
// integrated/summed into the RHS of the Helmholtz problem.  The
// "augment" methods are used by mixed BCs to add terms into the
// diagonal of the global Helmholtz matrix; the fact that there are a
// number of different such methods reflects the fact that we might be
// using a direct or iterative solution method for the Helmholtz
// problem.  See [5].
//
// Current Semtex treatment reflects the fact that BCs may need
// re-evaluation at each time step - in fact, it is assumed that they
// could be remade each time step (though for constant values, this
// has no effect). In addition, there are now a number of computed
// mixed BC types, see [3] and [4], in addition to the computed
// natural BCs introduced in [1].
//
// NB: if you invent a new BC which is of essential/Dirichlet type,
// you also need to edit mesh.cpp/buildLiftMask() so that assemblymap
// sets the global numbering mask to 1 for this type.  Otherwise the
// BC won't be correctly applied.  In this consideration, Robin/mixed
// count as natural/Neumann (nothing needs be done in mesh.cpp).
//
// REFERENCES
// ----------
// [1] Karniadakis, Israeli & Orszag (1991), JCP 91(2)
// [2] Blackburn (2001), Comput Chem Eng 25(2/3)
// [3] Dong (2015), JCP 302:300-328 
// [4] Liu, Xie & Dong (2020), IJHFF 151:119355
// [5] Blackburn, Lee, Albrecht & Singh (2019) CPC 254.
///////////////////////////////////////////////////////////////////////////////

#include <sem.h>


///////////////////////////////////////////////////////////////////////////////
// We first give the semi-abstract shared methods: "set" for essential
// and "sum" for natural and mixed types, and various "augment**" for
// mixed.
///////////////////////////////////////////////////////////////////////////////


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


void Mixed::sum (const int_t   side  ,
		 const int_t*  bmap  ,
		 const real_t* src   ,
		 const real_t* weight,
		 real_t*       work  ,
		 real_t*       tgt   ) const
// ---------------------------------------------------------------------------
// Add boundary-integral src terms into globally-numbered tgt.  For a
// mixed BC, src is essentially the equivalent RHS of dc/dn + K*c =
// K*C.  The method is a copy of the Natural one of the same name.
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


///////////////////////////////////////////////////////////////////////////////
// Next, methods for the terminal concrete classes, including constructors.
// Start with Essential BCs.
///////////////////////////////////////////////////////////////////////////////


EssentialConstant::EssentialConstant (const char* v)
  : _value (strtod (v, 0))
// ---------------------------------------------------------------------------
//
// ---------------------------------------------------------------------------
{
  _descriptor = "essential-constant:\t" + to_string (_value);
}


void EssentialConstant::evaluate (const Field*   src    ,
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


EssentialFunction::EssentialFunction (const char* f)
// ---------------------------------------------------------------------------
//
// ---------------------------------------------------------------------------  
{
  strcpy ((_function = new char [strlen (f) + 1]), f);
  _descriptor  = "essential-function:\t";
  _descriptor += _function;
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


///////////////////////////////////////////////////////////////////////////////
// Natural BCs.
///////////////////////////////////////////////////////////////////////////////


NaturalConstant::NaturalConstant (const char* v) : _value (strtod (v, 0))
// ---------------------------------------------------------------------------
//
// ---------------------------------------------------------------------------  
{
  _descriptor = "natural-constant:\t" + to_string (_value);
}


void NaturalConstant::evaluate (const Field*   src    ,
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


NaturalFunction::NaturalFunction (const char* f)
// ---------------------------------------------------------------------------
// Natural condition that sets value by interpreting function string,
// which can be an explicit function of x, y, z & t as well as any
// installed token values.  Otherwise the same as Natural class.
// ---------------------------------------------------------------------------
{
  strcpy ((_function = new char [strlen (f) + 1]), f);
  _descriptor  = "natural-function:\t";
  _descriptor += _function;
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


NaturalComputed::NaturalComputed (BCmgr*     B,
				  const char t)
  : _BCmgr (B), _tag(t)
// ---------------------------------------------------------------------------
//
// ---------------------------------------------------------------------------  
{
  if (_tag != 'p')
    message ("NaturalComputed::NaturalComputed",
	     "tag p (pressure) is the only recognised option", ERROR);
  
  _descriptor =  "natural-computed-pressure, see KIO91";
}


void NaturalComputed::evaluate (const Field*   src    ,
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
  if (Fourier && _tag == 'p') _BCmgr -> evaluateCNBCp (id, plane, step, tgt); 
}


///////////////////////////////////////////////////////////////////////////////
// Mixed BCs.
///////////////////////////////////////////////////////////////////////////////


MixedConstant::MixedConstant (const char* v)
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

  sprintf (buf, "mixed-constant:\t%g\t%g", _K_, _C_);

  _descriptor  = buf;
}


void MixedConstant::evaluate (const Field*   src    ,
			      const int_t    id     , 
			      const int_t    plane  , 
			      const Element* elmt   ,
			      const int_t    side   , 
			      const int_t    step   ,
			      const bool     Fourier,
			      real_t*        tgt    ) const
// ---------------------------------------------------------------------------
// Load external tgt as RHS of dc/dn + K*c = K*C.
// Physical space.
// ---------------------------------------------------------------------------
{
  if (!Fourier)
    Veclib::fill (Geometry::nP(), _K_ * _C_, tgt, 1);
}


//////////////////////////////////////////////////////////////////////////////
// Internally computed mixed BC types follow.
//
// These are all evaluated in Fourier space.
//////////////////////////////////////////////////////////////////////////////


MixedComputed::MixedComputed (BCmgr*     B,
			      const char t)
  : _BCmgr (B), _tag (t)
// ---------------------------------------------------------------------------
// Computed mixed/Robin BCs for open boundaries, Dong (2015) eqs. (37,38)
// and Liu, Xie & Dong (2020) eq (16b).
//
// NOTE: there is a potential problem here for the velocity and scalar
// computations, since alpha[0] (and hence _K_) depends on
// time-stepping order, so it will likely be incorrect for the first
// time-step (or two).  However, it is fine for setting up matrices
// for the asymptotic timestep order (and the BCs are approximate
// anyhow).
// ---------------------------------------------------------------------------
{
  const char routine[] = "MixedComputed::MixedComputed";
  char       buf[STR_MAX];
  real_t     alpha[5];

  Integration::StifflyStable (Femlib::ivalue ("N_TIME"), &alpha[0]);  

  switch (_tag) {
  case 'p': {
    _K_         = 1.0 / Femlib::value ("KINVIS*DONG_DO");
    _descriptor = "mixed-computed-pressure, Dong 2015 eq. (37)";
  }
  case 'u': {
    _K_         = alpha[0] * Femlib::value ("DONG_DO/D_T");
    _descriptor = "mixed-computed-velocity-u, Dong 2015 eq. (38)";    
  }
  case 'v': {
    _K_         = alpha[0] * Femlib::value ("DONG_DO/D_T");
    _descriptor = "mixed-computed-velocity-v, Dong 2015 eq. (38)";    
  }
  case 'w': {
    _K_         = alpha[0] * Femlib::value ("DONG_DO/D_T");
    _descriptor = "mixed-computed-velocity-w, Dong 2015 eq. (38)";    
  }
  case 'W': {
    _K_         = alpha[0] * Femlib::value ("DONG_DO/D_T");
    _descriptor = "mixed-computed-velocity-W, special case with RHS = 0";    
  }
  case 'c': {
    _K_         = alpha[0] * Femlib::value ("DONG_DO/D_T");
    _descriptor = "mixed-computed-scalar, Liu, Xie & Dong 2020 eq. (16b)";
  }
  default: {
    sprintf (buf, "unrecognised constructor tag (%c)", t);
    message (routine, buf, ERROR);
  }
  }
}


void MixedComputed::evaluate (const Field*   src    ,
			      const int_t    id     , 
			      const int_t    plane  , 
			      const Element* elmt   ,
			      const int_t    side   , 
			      const int_t    step   ,
			      const bool     Fourier,
			      real_t*        tgt    ) const
// ---------------------------------------------------------------------------
// Load external tgt via calls to BCmgr methods.
//   
// NB the special case 'W' is a cut-down version of 'w', which does
// nothing for the boundary integral terms.  In effect this means that
// C = 0 and we have dw/dn + K w = 0.  We may do away with this case
// in future.
// ---------------------------------------------------------------------------
{
  if (Fourier)
    switch (_tag) {
      
    case 'p': _BCmgr -> evaluateCMBCp (src, id, plane, step,      tgt);
      
    case 'u': _BCmgr -> evaluateCMBCu (src, id, plane, step, 'u', tgt);
    case 'v': _BCmgr -> evaluateCMBCu (src, id, plane, step, 'v', tgt);
    case 'w': _BCmgr -> evaluateCMBCu (src, id, plane, step, 'w', tgt);

    case 'c': _BCmgr -> evaluateCMBCc (src, id, plane, step,      tgt);

    case 'W': Veclib::fill (Geometry::nP(), 0.0, tgt, 1);

    }
}
