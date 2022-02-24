#ifndef ASSEMBLYMAP_H
#define ASSEMBLYMAP_H

class AssemblyMap
// ===========================================================================
// This class is a holder for Field Element-boundary numbering
// schemes.  Different Fields can potentially have unique
// AssemblyMaps, or they may share them with other Fields, according
// to distribution of essential BCs.
//
// The numbering system is associated with the Geometry class, but
// augments it by allowing for connectivity and boundary conditions.
//
// _Bmask and _emask provide a hierarchy of descriptions of element-
// boundary-node essential boundary condition masks:
//
// 1. _bmask is a vector of int_t, 4*(np-1)*nel (i.e. nbndry) in length,
//    which describe the CCW traverse of element-boundary nodes.
//    Values set to 1 indicate that the node has an imposed
//    (essential) BC, values set to 0 indicate that the node
//    either has a natural or no BC.
//
// 2. _emask is the next level of the hierarchy, a vector of int_ts,
//    nel long.  Values are set to 1 if the corresponding element has
//    any bmask values set to true, otherwise set to 0.
//
// (On first acquaintance it might seem more logical to use vectors of
// bool type for these masks, but one cannot obtain an unique byte
// address of an individual bool in a vector, which turns out to be an
// issue.  We could perhaps have used a short int, instead, but that
// would break our veclib conventions, where int_t (same as int) is
// the smallest integer type considered.)
//
// _Nglobal specifies the number of global nodes, while _nsolve ( <=
// _nglobal) specifies the number of global nodes which have _bmask
// values of 0, i.e. where the nodal values will be solved for instead
// of set explicitly.
//
// _Nglobal and _nsolve also effectively provide a top level in the
// _bmask/_emask hierarchy, since if their difference is non-zero then
// at least one _bmask value (and _emask value) will be 1.
// 
// _Btog gives global node numbers to Element-boundary nodes, same
// length as _bmask (i.e. _nbndry).  The optimization level describes
// the scheme used to build _btog.
//
// An AssemblyMap may be uniquely identified by the entries in _btog,
// or equivalently the entries of _bmask, together with optimization
// level/global elliptic problem solution strategy.
// ===========================================================================
{
public:
  AssemblyMap (const int_t, const int_t, const int_t,
	       const vector<int_t>&, const vector<int_t>&,
	       const char, const int_t);
 ~AssemblyMap () { };

  bool         willMatch (const vector<int_t>&) const;
  void         addTag    (const char, const int_t);
  bool         matchTag  (const char, const int_t) const;

  int_t        nGlobal () const { return _nglobal; }
  int_t        nSolve  () const { return _nsolve;  }
  int_t        nBand   () const { return _nbandw;  }

  const int_t* bmask () const { return &_bmask[0];         }
  const int_t* emask () const { return &_emask[0];         }
  const int_t* btog  () const { return &_btog[0];          }
  int_t        fmask () const { return _nglobal - _nsolve; }
  
private:
  int_t _optlev ;	  // Optimization level used for btog.
  int_t _nglobal;	  // Number of unique globally-numbered nodes.
  int_t _nbndry ;	  // Number of element-edge nodes.
  int_t _nsolve ;	  // Number of non-masked global nodes.
  int_t _nbandw ;	  // Bandwidth of btog (includes diagonal).
  int_t _np     ;	  // Number of nodes on an element edge.
  int_t _nel    ;	  // Number of elements;

  vector<int_t> _bmask;	  // 1 for essential-BC nodes, 0 otherwise.
  vector<int_t> _emask;	  // 1 if associated Element has any esstl set.
  vector<int_t> _btog ;	  // Assembly mapping for all element-boundary nodes.

  vector<pair<char, int_t>* > _tag; // Field name and mode tag pairs.

  int_t sortGid         (vector<int_t>&, const vector<int_t>&);
  int_t buildAdjncySC   (vector<int_t>&, vector<int_t>&, const int_t = 0) const;
  int_t globalBandwidth () const;
  void  connectivSC     (vector<vector<int_t> >&, const int_t*,
			 const int_t*, const int_t) const;
  int_t bandwidthSC     (const int_t*, const bool*, const int_t) const;
  void  RCMnumbering    ();
};

#endif
