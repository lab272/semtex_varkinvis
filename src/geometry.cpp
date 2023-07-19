///////////////////////////////////////////////////////////////////////////////
// geometry.cpp: define geometrical properties for 2D quad X Fourier spaces.
//
// Most routines are inlined in header file geometry.h
//
// This is effectively a singleton class.  We might want to put the
// static 'private' data somewhere else that is only loaded by
// top-level executables in order to prevent possible conflicts if we
// override the definitions in this file elsewhere.
//
// Copyright (c) 1994+, Hugh M Blackburn
///////////////////////////////////////////////////////////////////////////////

#include <cstdio>
#include <iostream>

#include <cfemdef.h>
#include <utility.h>
#include <geometry.h>
#include <veclib.h>
#include <femlib.h>

int_t Geometry::_pid   = UNSET;	// -- Initialise static private data.
int_t Geometry::_nproc = UNSET;
int_t Geometry::_ndim  = UNSET;
int_t Geometry::_np    = UNSET;
int_t Geometry::_nz    = UNSET;
int_t Geometry::_nzp   = UNSET;
int_t Geometry::_nel   = UNSET;
int_t Geometry::_psize = UNSET;
Geometry::CoordSys Geometry::_csys = Geometry::Cartesian;


static int_t roundUp (const int_t n, const int_t a, const int_t b)
// ---------------------------------------------------------------------------
// Return the first integer greater than or equal to n that has both
// factors a and b.  We assume that all inputs are positive.
// ---------------------------------------------------------------------------
{ int_t m = n; while (m%a || m%b) m++; return m; }


void Geometry::set (const int_t    NP,
		    const int_t    NZ,
		    const int_t    NE,
		    const CoordSys CS)
// ---------------------------------------------------------------------------
// Load values of static internal variables.
//
// The number of processors is restricted: it must either be 1 or be
// be less than or equal to the number of planes / 2.  Furthermore,
// the number of planes on a processor must be even, unless NZ == 1
// (because each Fourier mode is taken to have both real and imaginary
// parts).  Hence, NZ is always even if NZ > 1.
//
// NB: the value of _psize (a.k.a. planeSize) is the value of nPlane
// (nel*np*np), but rounded up if necessary such that it remains an
// even number when divided by the number of processes because:
//
// 1. The even number restriction is to simplify the handling of
// Fourier transforms, which is typically based on a real--complex
// transform (done via the method of transform of two real functions
// simultaneously, see e.g. Numerical Recipes or Bendat & Piersol).
//  
// 2. The restriction to be an integer multiple of twice the number of
// processors is to simplify the structure of memory exchanges
// required for Fourier transforms when computing in parallel.
// ---------------------------------------------------------------------------
{
  static char routine[] = "Geometry::set", err[StrMax];

  if (_np != UNSET)
    Veclib::alert (routine, "cannot re-initialise Geometry", ERROR);

  _pid   = Femlib::ivalue ("I_PROC");
  _nproc = Femlib::ivalue ("N_PROC");

  _np   = NP; _nz = NZ; _nel = NE; _csys = CS;
  _nzp  = _nz / _nproc;
  _ndim = (_nz > 2) ? 3 : 2;

  if (_nz > 1 && _nz & 1) {	// -- 3D problems must have NZ even.
    sprintf (err, "N_Z must be even (%1d)", _nz);
    Veclib::alert (routine, err, ERROR);
  }

  if (_nproc > 1) {		// -- Concurrent execution restrictions.

    if (_nz % (2 * _nproc)) {
      sprintf (err, "No. of planes (%1d) per processor (%1d) must be even",
	       _nz, _nproc);
      Veclib::alert (routine, err, ERROR);
    }

    if (_nproc << 1 > _nz) {
      sprintf (err, "No. of processors (%1d) can at most be half N_Z (%1d)",
	       _nproc, _nz);
      Veclib::alert (routine, err, ERROR);
    }
#if 1
    // -- _psize when divided by number of processes must be even.
    _psize = roundUp (nPlane(), 2*_nproc, 2);
#else
    // -- Old version which didn't use roundUp().
    _psize  = nPlane();
    _psize += 2 * _nproc - nPlane() % (2 * _nproc);   
#endif
    
  } else {
    if (_nz > 1)
      _psize = roundUp (nPlane(), 1, 2);
    else
      _psize = nPlane();
  }
}
