#ifndef NUMBERSYS_H
#define NUMBERSYS_H

class NumberSys
// ===========================================================================
// This class is a holder for Field Element-boundary numbering
// schemes.  Different Fields can potentially have unique
// NumberSyss, or they may share them with other Fields, according
// to distribution of essential BCs.
//
// The numbering system is associated with the Geometry class, but
// augments it by allowing for connectivity and boundary conditions.
//
// Bmask and emask provide a hierarchy of descriptions of element-
// boundary-node essential boundary condition masks:
//
// 1. bmask is a vector of bool, 4*(np-1)*nel (i.e. nbndry) in length,
//    which describe the CCW traverse of element-boundary nodes.
//    Values set to true indicate that the node has an imposed
//    (essential) BC, values set to false indicate that the node
//    either has a natural or no BC.
//
// 2. emask is the next level of the hierarchy, a vector of ints nel
//    long.  Values are set to tue if the corresponding element has
//    any bmask values set to true, otherwise set to false.
//
// _Nglobal specifies the number of global nodes, while _nsolve ( <=
// _nglobal) specifies the number of global nodes which have _bmask
// values of false, i.e. where the nodal values will be solved for
// instead of set explicitly.
//
// _Nglobal and _nsolve also effectively provide a top level in the
// _bmask/_emask hierarchy, since if their difference is non-zero then
// at least one _bmask value (and _emask value) will be true.
// 
// _Btog gives global node numbers to Element-boundary nodes, same
// length as _bmask (i.e. _nbndry).  The optimization level describes the
// scheme used to build _btog.
//
// Inv_mass is the inverse of the values of the global mass matrix,
// but on Element boundaries only.  Used for Field smoothing
// operations after derivative operators (if required).  Length =
// nglobal.
//
// A NumberSys is uniquely identified by the entries in _btog, or
// equivalently the entries of _bmask, together with optimization
// level/global elliptic problem solution strategy.
// ===========================================================================
{
public:
  NumberSys (const int_t, const int_t, const int_t,
	     const vector<int_t>&, const vector<int_t>&);
 ~NumberSys () { };

  bool          willMatch   (const vector<int_t>&) const;

  int_t         nGlobal () const { return _nglobal; }
  int_t         nSolve  () const { return _nsolve;  }
  int_t         nBand   () const { return _nbandw;  }

  const int_t*  bmask   () const { return &_bmask[0];         }
  const int_t*  emask   () const { return &_emask[0];         }
  const int_t*  btog    () const { return &_btog[0];          }
  int_t         fmask   () const { return _nglobal - _nsolve; }
  
#if 0  
  const real_t* imass   () const { return _imass;             }
#endif
  
private:
  int_t   _optlev ;		// Optimization level used for btog.
  int_t   _nglobal;		// Number of unique globally-numbered nodes.
  int_t   _nbndry ;		// Number of element-edge nodes.
  int_t   _nsolve ;		// Number of non-masked global nodes.
  int_t   _nbandw ;		// Bandwidth of btog (includes diagonal).
  int_t   _np     ;		// Number of nodes on an element edge.
  int_t   _nel    ;		// Number of elements;

  vector<int_t> _bmask ;        // 1 for essential-BC nodes, 0 otherwise.
  vector<int_t> _emask ;	// 1 if associated Element has any esstl set.
  vector<int_t> _btog  ;	// Gives numbers to all element-boundary knots.

  int_t sortGid         (vector<int_t>&, vector<int_t>&);
  int_t buildAdjncySC   (vector<int_t>&, vector<int_t>&, const int_t = 0) const;
  int_t globalBandwidth () const;
  int_t bandwidthSC     (const int_t*, const bool*, const int_t) const;
  void  RCMnumbering    ();
  
#if 0  
  real_t* _imass ;		// Inverse of global mass matrix -- to Domain.
#endif
};

#endif
