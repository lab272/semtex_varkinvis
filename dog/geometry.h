#ifndef GEOMETRY_H
#define GEOMETRY_H

#include <cfemdef.h>

class Geometry
// ===========================================================================
// Details of geometric representation used for scalar fields.  Static
// functions make information globally accessible.
//
// In all cases, 2D quad elements are employed, with a possible
// extension by Fourier expansions in the third dimension.  While the
// representation implied by this class is not necessarily conforming,
// the same order of interpolation is used in each element.
// Equal-order interpolation is used in each direction on faces of
// quads.
//
// With the introduction of concurrent execution, the concept of
// Geometry has been extended to include the processor ID, number of
// processors, number of data planes per processor, etc.
//
// This version adds a number of functions that are used in
// eigensystem analysis.
//
// Copyright (c) 1994+, Hugh M Blackburn
// ===========================================================================
{
public:
  enum CoordSys { Cartesian, Cylindrical };
  enum Category { O2_2D, O2_3D, O2_3D_SYMM, SO2_2D, SO2_3D }; // ->README file.

  static void set (const int_t, const int_t);

  static CoordSys    system()       { return _Csys;               }
  static Category    problem()      { return _Cat;                }
  static const char* symmetry();    // -- Return string corresponding to _Cat.
  static bool        cylindrical () { return _Csys == Geometry::Cylindrical; }

  static int_t nP        () { return _Np;                   }
  static int_t nZ        () { return _Nz;                   }
  static int_t nElmt     () { return _Nel;                  }
  static int_t nTotElmt  () { return _Np * _Np;             }
  static int_t nExtElmt  () { return 4 * (_Np - 1);         }
  static int_t nIntElmt  () { return (_Np - 2) * (_Np - 2); }
  static int_t nMode     () { return (_Nz + 1) >> 1;        }
  static int_t nDim      () { return _Npert;                }

  static int_t nPlane    () { return _Nel * nTotElmt();     }
  static int_t nBnode    () { return _Nel * nExtElmt();     }
  static int_t nInode    () { return _Nel * nIntElmt();     }
  static int_t nTot      () { return _Nz  * nPlane();       }
  static int_t planeSize () { return _Psize;                }
  static int_t nTotal    () { return _Nz * _Psize;          }

  static int_t nProc     () { return _Nproc;                }
  static int_t procID    () { return _Pid;                  }
  static int_t nZProc    () { return _Nzp;                  }
  static int_t nZ32      () { return (_Nproc > 1) ? _Nzp : (3 * _Nz) >> 1; }
  static int_t nTotProc  () { return _Nzp * _Psize;         }
  static int_t nModeProc () { return nMode() / _Nproc;      }
  static int_t baseMode  () { return _Pid * nModeProc();    }
  static int_t basePlane () { return _Pid * _Nzp;           }
  static int_t nBlock    () { return _Psize / _Nproc;       }
  
  // -- These are specific to eigensystem analysis:

  static int_t nBase     () { return _Nbase;                }
  static int_t nPert     () { return _Npert;                }
  static int_t nSlice    () { return _Nslice;               }

private:
  static int_t    _Nproc ;     // Number of processors.
  static int_t    _Pid   ;     // ID for this processor, starting at 0.
  static int_t    _ndim  ;     // Number of space dimensions
  static int_t    _Np    ;     // Number of points along element edge.
  static int_t    _Nz    ;     // Number of planes (total).
  static int_t    _Nzp   ;     // Number of planes per processor.
  static int_t    _Nel   ;     // Number of elements.
  static int_t    _Psize ;     // nPlane rounded up to suit restrictions.
  static int_t    _kfund ;     // Wavenumber of first non-zero Fourier mode.
  static CoordSys _Csys  ;     // Coordinate system (Cartesian/cylindrical).
  static Category _Cat   ;     // Problem category.

  // -- These are specific to eigensystem analysis:

  static int_t    _Npert ;     // Number of perturbation velocity components.
  static int_t    _Nbase ;     // Number of base velocity field components.
  static int_t    _Nslice;     // Number of base velocity fields.

};
#endif
