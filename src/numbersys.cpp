///////////////////////////////////////////////////////////////////////////////
// numbersys.cpp: NumberSys class functions.
//
// A NumberSys automates the retrieval of global assembly mapping
// infomation for element-edge nodes appropriate to a given Field and
// Fourier mode.
//
// The information to be returned by class functions are the global
// numbering scheme (assembly map) for a given Field and Fourier mode
// number.  There is one NumberSys for each Field, but a possible
// modal dependence for the appropriate BCs, and hence global assembly
// mapping schemes (in fact, this occurs only for 3D cylindrical
// coordinate systems in which the axis appears, but we assume a modal
// dependence is always possibly the case).
//
// On every process, these are constructed for Fourier modes 0, 1,
// 2(+); subsequently, the correct local system is selected based on
// the requested zero-based mode index.
//
// Copyright (c) 2022+, Hugh M Blackburn
///////////////////////////////////////////////////////////////////////////////

#include <sem.h>


NumberSys::NumberSys (const vector<AssemblyMap*>& allMaps,
		      const char                  name   ) : 
// ---------------------------------------------------------------------------
// Construct vectors of modal assembly map systems for all Fourier
// modes of a named field.  Input value "name" is one of the standard
// Field names: "uvwcp". On entry, the numbering schemes for all
// fields and modes are available in input variable allMaps, which is
// created and stored by the Domain class (it is possible for
// different fields and Fourier modes to share the same assembly map).
// ---------------------------------------------------------------------------
  _field_name (name)
{
  const char routine[] = "NumberSys::NumberSys"; 
  char       err[StrMax];
  int_t      mode, j;
  bool       found;

  _map.resize(3);
  for (mode = 0; mode < 3; mode++)
    
    for (found = false, j = 0; j < allMaps.size(); j++)
      if (found = (allMaps[j] -> matchTag (name, mode))) {
	_map[mode] = allMaps[j];
	break;
	
	if (!found) {			// -- Never happen.
	  sprintf (err,
		   "field %c, mode %1d: can't find matching AssemblyMap",
		   name, mode);
	  Veclib::messg (routine, err, ERROR);
	}
      }
}


const AssemblyMap* NumberSys::getMap (const int_t mode) const
// ---------------------------------------------------------------------------
// Return appropriate AssemblyMap, according to Fourier mode.  Input
// mode is the actual global mode index, counting from zero, not
// modulo number of modes on this process.
// ---------------------------------------------------------------------------
{
  return _map [clamp (mode, static_cast<int_t>(0), static_cast<int_t>(2))];
}
