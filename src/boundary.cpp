//////////////////////////////////////////////////////////////////////////////
// boundary.cpp: implement Boundary class functions.
//
// SYNOPSIS
// --------
// Boundaries correspond to domain edges that have boundary conditions
// applied (as opposed to periodic edges).  The ordering of internal
// storage for condition values and geometric factors corresponds to a
// CCW traverse of 2D element edges.
//
// Copyright (c) 1994+, Hugh M Blackburn
///////////////////////////////////////////////////////////////////////////////

#include <sem.h>


void Boundary::evaluate (const Field* P      ,
			 const int_t  plane  ,
			 const int_t  step   ,
			 const bool   Fourier,
			 real_t*      tgt    ) const
// ---------------------------------------------------------------------------
// Load boundary condition storage area with numeric values.
// ---------------------------------------------------------------------------
{
  _bcond -> evaluate (P, _id, plane, _elmt, _side, step, Fourier, tgt);
}


void Boundary::set (const real_t* src,
		    const int_t*  b2g,
		    real_t*       tgt) const
// ---------------------------------------------------------------------------
// Use (boundary condition) values in src to over-ride (set) values
// in globally-numbered tgt.  This will only take place on essential BCs.
//
// b2g is a pointer to the global node numbers for the appropriate
// element's edge nodes.
// ---------------------------------------------------------------------------
{
  _bcond -> set (_side, b2g, src, tgt);
}


void Boundary::sum (const real_t* src,
		    const int_t*  b2g,
		    real_t*       wrk,
		    real_t*       tgt) const
// ---------------------------------------------------------------------------
// Use (boundary condition) values in src to add in the boundary-integral
// terms generated in constructing the weak form of the MWR into globally-
// numbered tgt.  This will only take place on natural BCs.
//
// b2g is a pointer to the global node numbers for the appropriate
// element's edge nodes.  wrk is a work array, np long.
// ---------------------------------------------------------------------------
{
  _bcond -> sum (_side, b2g, src, _area, wrk, tgt);
}


void Boundary::augmentSC (const int_t  nband ,
			  const int_t  nsolve,
			  const int_t* b2g   ,
			  real_t*      work  ,
			  real_t*      H     ) const
// ---------------------------------------------------------------------------
// Add in diagonal terms <K, w> to (banded LAPACK) H on mixed BCs.
// Work array must be np long.  See p. 135 of HMB's Semtex notebook.
// ---------------------------------------------------------------------------
{
  _bcond -> augmentSC (_side, nband, nsolve, b2g+this->bOff(), _area, work, H);
}


void Boundary::augmentOp (const int_t*  b2g,
			  const real_t* src,
			  real_t*       tgt) const
// ---------------------------------------------------------------------------
// This operation is used to augment the element-wise Helmholtz
// operations where there are mixed BCs.  Add in diagonal terms
// <K*src, w> to tgt.  Both src and tgt are globally-numbered vectors.
// ---------------------------------------------------------------------------
{
  _bcond -> augmentOp (_side, b2g + this -> bOff(), _area, src, tgt);
}


void Boundary::augmentDg (const int_t* b2g,
			  real_t*      tgt) const
// ---------------------------------------------------------------------------
// This operation is used to augment the element-wise construction of
// the diagonal of the global Helmholtz matrix where there are mixed
// BCs.  Add in diagonal terms <K, w> to globally-numbered tgt.
// ---------------------------------------------------------------------------
{
  _bcond -> augmentDg (_side, b2g + this -> bOff(), _area, tgt);
}


void Boundary::print () const
// ---------------------------------------------------------------------------
// (Debugging) utility to print internal information.
// ---------------------------------------------------------------------------
{
  char  info[StrMax];
  int_t i;

  cout << "** Boundary id: " << _id + 1 << " -> ";
  cout << _elmt ->  ID() + 1 << "." << _side + 1;
  cout << " (Element id.side)" << endl;
  
  _bcond -> describe (info);

  cout << info << endl;

  cout << "  " << _np << " (number of points along edge)" << endl;
  cout << "         nx             ny             area";
  cout << endl;

#if 1
  cout << setw(15);
  for (i = 0; i < _np; i++)
    cout << _nx[i] << _ny[i] << _area[i] << endl;
#else  
  printVector (cout, "rrr", _np, _nx, _ny, _area);
#endif
}



