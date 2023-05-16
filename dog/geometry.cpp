///////////////////////////////////////////////////////////////////////////////
// geometry.cpp: define geometrical properties for 2D quad X Fourier spaces.
//
// Most routines are inlined in header file geometry.h.
//
// NB: we have capitalized names of all static member data so as to
// avoid potential conflicts with those instantiated in
// src/geometry.cpp (which this file is designed to override).
//
// Copyright (c) 1994+, Hugh M Blackburn
///////////////////////////////////////////////////////////////////////////////

#include <cstdio>
#include <iostream>

#include <cfemdef.h>
#include <utility.h>
#include <veclib.h>
#include <femlib.h>
#include <geometry.h>

int_t Geometry::_Pid    = 0;
int_t Geometry::_Nproc  = 0;
int_t Geometry::_Np     = 0;
int_t Geometry::_Nz     = 0;
int_t Geometry::_Nzp    = 0;
int_t Geometry::_Nel    = 0;
int_t Geometry::_Psize  = 0;

int_t Geometry::_Npert  = 0;
int_t Geometry::_Nbase  = 0;
int_t Geometry::_Nslice = 0;
Geometry::CoordSys Geometry::_Csys = Geometry::Cartesian;
Geometry::Category Geometry::_Cat  = Geometry::O2_3D_SYMM;


void Geometry::set (const int_t nel  ,
		    const int_t npert)
// ---------------------------------------------------------------------------
// Load values of static internal variables.  Session file should
// already have been dealt with.  As well as being specific to stability
// analysis, this version is written for serial execution.
// ---------------------------------------------------------------------------
{
  static char routine[] = "Geometry::set";
  char        err[StrMax];

  _Pid    = Femlib::ivalue ("I_PROC");
  _Nproc  = Femlib::ivalue ("N_PROC");

  _Np     = Femlib::ivalue ("N_P");

  _Nbase  = Femlib::ivalue ("N_BASE");
  _Nslice = Femlib::ivalue ("N_SLICE");
  _Csys   = Femlib::ivalue ("CYLINDRICAL") ? 
                     Geometry::Cylindrical : Geometry::Cartesian;
  _Npert  = npert;
  _Nel    = nel;
  _Psize  = nPlane() + (nPlane() % 2);

  _Nz = _Nzp = Femlib::ivalue ("N_Z");

  if      (_Nbase == 2 && _Npert == 2 && _Nz == 1) _Cat = O2_2D;
  else if (_Nbase == 2 && _Npert == 3 && _Nz == 1) _Cat = O2_3D_SYMM;
  else if (_Nbase == 2 && _Npert == 3 && _Nz == 2) _Cat = O2_3D;
  else if (_Nbase == 3 && _Npert == 3 && _Nz == 1) _Cat = SO2_2D;
  else if (_Nbase == 3 && _Npert == 3 && _Nz == 2) _Cat = SO2_3D;
  else {
    sprintf (err, "illegal: N_BASE = %1d, N_PERT = %1d, N_Z = %1d",
	     _Nbase, _Npert, _Nz); message (routine, err, ERROR);
  }

  // -- We can only allow SO2_2D if BETA != 0.

  if ((_Cat == SO2_2D) && (Femlib::value ("BETA") > EPSDP))
    message
      (routine, "SO2_2D needs BETA=0, i.e. z-invariance, too.", ERROR);

  // -- Other sanity checks.

  if (_Nproc  > 1)
    message (routine, "serial execution only",            ERROR);
  if (_Nslice < 1)
    message (routine, "N_SLICE must be set in session",   ERROR);
}


const char* Geometry::symmetry ()
// ---------------------------------------------------------------------------
// Return a string value of the spatial symmetry category of problem.
// ---------------------------------------------------------------------------
{
  switch (_Cat) {
  case O2_2D:      return "O(2), 2D"; break;
  case O2_3D_SYMM: return "O(2), 3D, standing wave"; break;
  case O2_3D:      return "O(2), 3D, travelling wave"; break;
  case SO2_2D:     return "SO(2), 2D"; break;
  case SO2_3D:     return "SO(2), 3D"; break;
  }

  return "Never happen";
}

