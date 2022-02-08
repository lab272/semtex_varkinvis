///////////////////////////////////////////////////////////////////////////////
// numbersys.cpp: NumberSys class functions.
//
// Copyright (c) 2022+, Hugh M Blackburn
//
//
///////////////////////////////////////////////////////////////////////////////

#include <sem.h>


NumberSys::NumberSys (const int_t          n_p     , // Element-edge N_P value.
		      const int_t          n_el    , // Number of quad elements.
		      const int_t          strat   , // Enumeration strategy.
		      const vector<int_t>& naiveMap, // Un-optimised assembly.
		      const vector<bool>&  liftMask) // Flagged for lifting.
// ---------------------------------------------------------------------------
// Create internal storage for a global assembly numbering scheme for a mesh
// of 2D quadrilateral elements.
//
// Input naiveMap is an un-optimised element-edge assembly mapping
// generated e.g. by Mesh::buildMap() without reference to BC
// information.  It is 4*(n_p-1)*n_el in length.
//
// Input liftMask (of the same length) flags element-edge nodes which
// will be lifted out of eventual elliptic solution unknowns; a true
// value indicates that the element edge data are to be supplied by an
// essential/Dirichlet boundray condition.  This mask vector could be
// generated e.g. by Mesh::buildMask(). The assembly numbers for the
// associated global nodes will be sorted to the end of the numbering
// scheme created here, and the associated values are not to be solved
// for.
//
// Input strat is at present just a flag for optimisation level to be
// applied to RCM bandwidth minimisation renumbering of global nodes.
//
// Much of what is done here was once buried in enumerate.cpp and BCmgr.cpp.  
//  
// ---------------------------------------------------------------------------
{
  const char routine[] = "NumberSys::NumberSys";

  if (naiveMap.size() != liftMask.size())
    message (routine, "sizes of input vectors don't match",    ERROR);
  if (naiveMap.size() != 4*(n_p-1)*n_el)
    message (routine, "input vectors are of improper lengths", ERROR);

  _optlevel = strat;
  _bmask    = liftMask;
  _btog     = naiveMap;

  _nbndry   = naiveMap.size();
  _nglobal  = naiveMap[Veclib::imax (naiveMap.size(), &naiveMap[0], 1)] + 1;
 
  _nsolve = this -> sortGid (_btog, _bmask);
}


static int cmp1 (const void *a,
		 const void *b)
// ---------------------------------------------------------------------------
// Used by qsort.  Compare first element (global node number) of two arrays.
// ---------------------------------------------------------------------------
{ return static_cast<int>
    (static_cast<const int_t *>(a)[0]-static_cast<const int_t *>(b)[0]); }


static int cmp2 (const void *a,
		 const void *b)
// ---------------------------------------------------------------------------
// Used by qsort.  Compare second element (solve mask) of two arrays.
// ---------------------------------------------------------------------------
{ return static_cast<int>
    (static_cast<const int_t *>(a)[1]-static_cast<const int_t *>(b)[1]); }


int_t NumberSys::sortGid (vector<int_t>&      bmap,
			  const vector<bool>& bmsk)
// ---------------------------------------------------------------------------
// Global node numbers get sorted to place essential-BC nodes last:
// this simplifies the later partition of global matrices.
//
// The non-essential type node numbers can be further sorted to
// optimize global matrix bandwidths, but this is not done here.
//
// A globally-numbered table (reOrder) is constructed, each entry of
// which is two ints: the global node number and the essential BC mask
// value (0/1).  This is then partitioned into "unknown" and "known"
// node numbers, with the "known" numbers last in the table.  Each
// partition is then sorted into ascending node number order, and the
// information is used to create a new element boundary-to-global
// numbering scheme in bmap.
//
// Return number of globally-numbered nodes at which solution is not
// set by essential BCs.
//
// This code is essentially lifted straight from enumerate.cpp and is
// some of the oldest in semtex.
// ---------------------------------------------------------------------------
{
  vector<int_t> work (2 * _nbndry + 3 * _nglobal);
  int_t         *mloc, *bsave, *tmp, *reOrder;
  int_t         unknowns;

  mloc    = &work[0];
  bsave   = mloc  + _nbndry
  tmp     = bsave + _nbndry;
  reOrder = tmp   + _nglobal;

  for (int_t i = 0; i < _nbndry, i++) mloc[i] = bmsk[i]; // Type promotion.

  Veclib::copy  (_nbndry,  &bmap[0], 1, bsave, 1);
  Veclib::scatr (_nbndry , mloc, bsave, tmp);
  Veclib::ramp  (_nglobal, 0,    1, reOrder,     2);
  Veclib::copy  (_nglobal, tmp,  1, reOrder + 1, 2);
  
  unknowns = _nglobal - Veclib::count (_nglobal, tmp, 1);

  if (unknowns < nglobal) {

    // -- Partition into "unknown" nodes & "essential BC" nodes.

    qsort (reOrder,              _nglobal,            2*sizeof (int_t), cmp2);

    // -- Sort each partition into ascending node number order.

    qsort (reOrder,               unknowns,           2*sizeof (int_t), cmp1);
    qsort (reOrder + 2*unknowns, _nglobal - unknowns, 2*sizeof (int_t), cmp1);

    // -- Reload new gids.

    Veclib::copy  (_nglobal, reOrder, 2,  tmp, 1);
    Veclib::ramp  (_nglobal, 0,    1, reOrder, 1);
    Veclib::scatr (_nglobal, reOrder, tmp, reOrder + _nglobal);
    Veclib::gathr (_nbndry , reOrder + _nglobal, bsave, &bmap[0]);
  }

  return unknowns;
}



