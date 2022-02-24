///////////////////////////////////////////////////////////////////////////////
// numbersys.cpp: NumberSys class functions.
//
// Copyright (c) 2022+, Hugh M Blackburn
//
// The information to be returned by class functions are the global
// numbering scheme for a given Field and Fourier mode number.  There
// is one NumberSys for each Field, but a possible modal dependence
// for the appropriate BCs, and hence global assembly mapping schemes
// (in fact, this occurs only for 3D cylindrical coordinate systems in
// which the axis appears, but we assume a modal dependence is always
// the case).
//
// (On every process, these are constructed for Fourier modes 0, 1, 2;
// subsequently, the right system is selected based on process ID.)
//
///////////////////////////////////////////////////////////////////////////////

#include <sem.h>


NumberSys::NumberSys (const vector<AssemblyMap*>& allMaps,
		      const char                  name   ) : 
// ---------------------------------------------------------------------------
// Construct vectors of modal assembly map systems for all modes.
// Input value "name" is one of the standard Field names: "uvwcp". On
// entry, the numbering schemes for all fields and modes is available
// in input variable allMaps.
//
// ---------------------------------------------------------------------------
  _field_name (name),
{
  const char  routine[] = "NumberSys::NumberSys"; 
  const int_t verbose = Femlib::ivalue ("VERBOSE");
  char        err[StrMax];
  int_t       i, j;
  bool        found;

  _map.resize(3);
  for (i = 0; i < 3; i++) {
    for (found = false, j = 0; j < allMap.size(); j++)
      if (found = (allMap[j] -> matchTag (name, i)))
	_map[i] = allMap[j];

    if (!found) {			// -- Never happen.
      sprintf (err,
	       "field %c, mode %1d: can't find matching AssemblyMap", name, i);
      message (routine, err, ERROR);
    }
  }
}


const AssemblyMap* NumberSys::getMap (const int_t mode) const
// ---------------------------------------------------------------------------
// Return appropriate AssemblyMap, according
// to Fourier mode.  Mode number is the actual number, counting from
// zero, not modulo number of modes on this process.
// ---------------------------------------------------------------------------
{
  return _map [clamp (mode,static_cast<int_t>(0),static_cast<int_t>(2))];
}

#endif
