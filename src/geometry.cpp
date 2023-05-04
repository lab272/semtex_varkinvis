///////////////////////////////////////////////////////////////////////////////
// geometry.cpp: define geometrical properties for 2D quad X Fourier spaces.
//
// Most routines are inlined in header file geometry.h
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

int_t Geometry::_pid   = 0;	// -- Initialise static private data.
int_t Geometry::_nproc = 0;
int_t Geometry::_ndim  = 0;
int_t Geometry::_np    = 0;
int_t Geometry::_nz    = 0;
int_t Geometry::_nzp   = 0;
int_t Geometry::_nel   = 0;
int_t Geometry::_psize = 0;
Geometry::CoordSys Geometry::_csys = Geometry::Cartesian;


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
// (nel*np*np), but rounded up if necessary to be an even number and
// also an integer multiple of the number of processors because:
//
// 1. The even number restriction is to simplify the handling of
// Fourier transforms, which is typically based on a real--complex
// transform (done via the method of transform of two real functions
// simultaneously, see e.g. Numerical Recipes or Bendat & Piersol).
//  
// 2. The restriction to be an integer multiple of the number of
// processors is to simplify the structure of memory exchanges
// required for Fourier transforms when computing in parallel.
// ---------------------------------------------------------------------------
{
  static char routine[] = "Geometry::set", err[StrMax];

  _pid   = Femlib::ivalue ("I_PROC");
  _nproc = Femlib::ivalue ("N_PROC");

  _np   = NP; _nz = NZ; _nel = NE; _csys = CS;
  _nzp  = _nz / _nproc;
  _ndim = (_nz > 2) ? 3 : 2;

  if (_nz > 1 && _nz & 1) {	// -- 3D problems must have NZ even.
    sprintf (err, "N_Z must be even (%1d)", _nz);
    Veclib::messg (routine, err, ERROR);
  }

  if (_nproc > 1) {		// -- Concurrent execution restrictions.

    if (_nz % (2 * _nproc)) {
      sprintf (err, "No. of planes (%1d) per processor (%1d) must be even",
	       _nz, _nproc);
      Veclib::messg (routine, err, ERROR);
    }

    if (_nproc << 1 > _nz) {
      sprintf (err, "No. of processors (%1d) can at most be half N_Z (%1d)",
	       _nproc, _nz);
      Veclib::messg (routine, err, ERROR);
    }

    _psize  = nPlane();
    _psize += 2 * _nproc - nPlane() % (2 * _nproc);

  } else {

    _psize = nPlane() + (nPlane() % 2);
  }
}
