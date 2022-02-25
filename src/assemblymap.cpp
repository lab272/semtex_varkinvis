///////////////////////////////////////////////////////////////////////////////
// assemblymap.cpp: AssemblyMap class functions.
//
// Copyright (c) 2022+, Hugh M Blackburn
//
///////////////////////////////////////////////////////////////////////////////

#include <sem.h>


AssemblyMap::AssemblyMap (const int_t          n_p     , // Element N_P value.
			  const int_t          n_el    , // Number of elements.
			  const int_t          strat   , // Enumeration strategy.
			  const vector<int_t>& naiveMap, // Un-optimised assembly.
			  const vector<int_t>& liftMask, // Flagged for lifting.
			  const char           name    , // Field name.
			  const int_t          mode    ) // Fourier mode.
  : _np (n_p), _nel (n_el), _optlev (strat)
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
// essential/Dirichlet boundary condition.  This mask vector could be
// generated e.g. by Mesh::buildMask(). The assembly numbers for the
// associated global nodes will be sorted to the end of the numbering
// scheme created here, and the associated values are not to be solved
// for.
//
// Input strat is at present just a flag for optimisation level to be
// applied to RCM bandwidth minimisation renumbering of global nodes.
//
// Inputs name and mode are initial name tags for this object.  If we
// later find that another field name+mode would get the same assembly
// map, then we add further tags to the present object rather than
// build a new one from scratch.
//
// Much of what is done here was once buried in enumerate.cpp and BCmgr.cpp.  
//  
// ---------------------------------------------------------------------------
{
  const char  routine[] = "AssemblyMap::AssemblyMap";
  const int_t next = _nbndry / _nel;
  int_t       i, j;

  if (naiveMap.size() != liftMask.size())
    message (routine, "sizes of input vectors don't match",    ERROR);
  if (naiveMap.size() != 4*(_np-1)*_nel )
    message (routine, "input vectors are of improper lengths", ERROR);

  _tag.resize (0);		// -- Should initially be this in any case...
  
  _btog  = naiveMap;
  _bmask = liftMask;

  // -- _emask: says if any external nodes on an element are essential/Dirichet.
  
  for (i = 0; i < _nel; i++) {
    _emask[i] = 0;
    for (j = 0; j < next; j++)
      if (_bmask[i*next + j]) {
	_emask[i] = 1;
	break;
      }
  }

  _nbndry   = naiveMap.size();
  _nglobal  = naiveMap[Veclib::imax (naiveMap.size(), &naiveMap[0], 1)] + 1;
 
  _nsolve = this -> sortGid (_btog, _bmask);

  if (_optlev > 0) this -> RCMnumbering ();

  _tag.push_back (new pair<char, int_t> (name, mode)); // -- Initial tags map.
}


bool AssemblyMap::willMatch (const vector<int_t>& candidate) const
// ---------------------------------------------------------------------------
// Return true if candidate is the same as internal storage _bmask ---
// in which case the eventual global numbering system _btog for a new
// AssemblyMap would come out to be the same if the construction
// strategy were also the same (which is assumed).
// ---------------------------------------------------------------------------
{
  const char routine[] = "AssemblyMap::willMatch";

  if (candidate.size() != _bmask.size()) {
    message (routine, "mismatched mask lengths", REMARK);
    return false;
  }

  return (Veclib::same (_bmask.size(), &_bmask[0], 1, &candidate[0], 1)) ?
    true : false;
}


void AssemblyMap::addTag (const char  name,
			  const int_t mode)
// ---------------------------------------------------------------------------
// This assembly map gets another name-tag pair added to aid later retrieval.
// ---------------------------------------------------------------------------
{
  _tag.push_back (new pair<char, int_t> (name, mode));
}


bool AssemblyMap::matchTag (const char  name,
			    const int_t mode) const
// ---------------------------------------------------------------------------
// Retrieval: does this AssemblyMap provide data for name and mode?
// ---------------------------------------------------------------------------
{
  for (int i = 0; i < _tag.size(); i++)
    if ((_tag[i] -> first == name) && (_tag[i] -> second == mode))
      return true;

  return false;
}


static inline int cmp1 (const void *a,
			const void *b)
// ---------------------------------------------------------------------------
// Used by qsort.  Compare first element (global node number) of two arrays.
// ---------------------------------------------------------------------------
{ return static_cast<int>
    (static_cast<const int_t *>(a)[0]-static_cast<const int_t *>(b)[0]); }


static inline int cmp2 (const void *a,
			const void *b)
// ---------------------------------------------------------------------------
// Used by qsort.  Compare second element (solve mask) of two arrays.
// ---------------------------------------------------------------------------
{ return static_cast<int>
    (static_cast<const int_t *>(a)[1]-static_cast<const int_t *>(b)[1]); }


int_t AssemblyMap::sortGid (vector<int_t>&       bmap,
			    const vector<int_t>& bmsk)
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
  bsave   = mloc  + _nbndry;
  tmp     = bsave + _nbndry;
  reOrder = tmp   + _nglobal;

  for (int_t i = 0; i < _nbndry; i++) mloc[i] = bmsk[i]; // Type promotion.

  Veclib::copy  (_nbndry,  &bmap[0], 1, bsave, 1);
  Veclib::scatr (_nbndry , mloc, bsave, tmp);
  Veclib::ramp  (_nglobal, 0,    1, reOrder,     2);
  Veclib::copy  (_nglobal, tmp,  1, reOrder + 1, 2);
  
  unknowns = _nglobal - Veclib::count (_nglobal, tmp, 1);

  if (unknowns < _nglobal) {

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


int_t AssemblyMap::buildAdjncySC (vector<int_t>& adjncy,
				  vector<int_t>& xadj  ,
				  const int_t    base  ) // 0/1-based numbering.
  const
// ---------------------------------------------------------------------------
// From internal tables, generate adjacency information in standard
// format (see e.g. George & Liu 1981) for element-boundary nodes
// (i.e. assuming element-level Static Condensation).  Because the
// information is going to be used subsequently for renumbering
// optimisation of global elliptic direct-solution schemes, nodes
// which correspond to Essential/Dirichlet BCs are not included (they
// get lifted out of solution) - this is determined by the associated
// internally-held mask vector.
//
// On entry, it is assumed that the numbering will already have been
// re-ordered by AssemblyMap::sortGid() to place numbers for the lifted
// nodes last/highest.  If no renumbering is called for, that
// re-ordering is sufficient for solution to proceed by either
// iterative/CG or direct/Cholesky methods.
//
// While the internally-held tables consist of 0-based arrays, the
// adjacency information can be returned with either 0- or 1-based
// numberings.
//
// Return the length of the adjacency table.  
//
// This routine is based on earlier ones from enumerate.cpp.
// ---------------------------------------------------------------------------
{
  int_t                  i, j, k, ntab; 
  vector<vector<int_t> > adjncyList (_nsolve);
  const int_t            next = _nbndry / _nel;

  for (k = 0, ntab = 0; k < _nel; k++) {
    this -> connectivSC (adjncyList, &_btog[0]+ntab, &_bmask[0]+ntab, next);
    ntab += next;
  }

  // -- We will return the total length of the adjncy vector:

  for (k = 0, ntab = 0; k < _nsolve; k++) ntab += adjncyList[k].size();

  // -- At this stage we have all the required information; unpack it:

  adjncy.resize (ntab    + 1);
  xadj.resize   (_nsolve + 1);

  for (k = 0, i = 0; i < _nsolve; i++) {
    xadj[i] = k + base;
    for (j = 0; j < adjncyList[i].size(); j++)
      adjncy[k++] = adjncyList[i][j] + base;
  }

  // -- Conventional terminal values (not used?):

  adjncy[k] = 0;
  xadj[i]   = k + base;
  
  return ntab;
}


void AssemblyMap::connectivSC (vector<vector<int_t> >& adjList,
			       const int_t*            bmap   ,
			       const int_t*            mask   ,
			       const int_t             next   )
  const
// ---------------------------------------------------------------------------
// AdjList is an array of int_t vectors, each of which describes the
// global nodes that have connectivity with the the current node,
// i.e. which make a contribution to the weighted-residual integral
// for this node.  This routine fills in the contribution from the
// current element, as determined by bmap (global numbers of element
// boundary nodes) and mask (which indicates Essential BC nodes, to be
// ignored).
//
// For general finite element Helmholtz matrices, all nodes of an
// element are interconnected, while for statically-condensed
// (high-order) elements, only the boundary nodes are considered
// (since internal nodes are not global).  Hence the traverse around
// the next external nodes.
//
// Essential-BC nodes (those with mask == true) are ignored, since
// we're only interested in mimimizing bandwidths of global matrices.
// ---------------------------------------------------------------------------
{
  int_t                   i, j, gidCurr, gidMate;
  bool                    found;
  vector<int_t>::iterator k;
  
  for (i = 0; i < next; i++) {
    if (! mask[i]) {
      gidCurr = bmap[i];
      
      for (j = 0; j < next; j++) {
	if (i != j && ! mask[j]) {
	  for (k = adjList[gidCurr].begin(), gidMate = bmap[j], found = false;
	       !found && k != adjList[gidCurr].end(); k++)
	    found = (*k == gidMate);
	  if (!found) adjList[gidCurr].push_back (gidMate);
	}
      }
    }
  }
}


void AssemblyMap::RCMnumbering ()
// ---------------------------------------------------------------------------
// For now, this is our only global system renumbering methodology,
// Reverse Cuthill--McKee (RCM), see e.g. George & Liu (1981), from
// where the routines used below are taken (also in SPARSEPAK).
//  
// From the initial ordering specified in bndmap, use RCM to generate
// a reduced-bandwidth numbering scheme.  Reload into _btog.
//
// Different optimization levels are allowed:
//
// 0: Do nothing (no renumbering).
// 1: Use FNROOT (trial root = 1) to find a pseudo-peripheral root node,
//    pass result to RCM for Reverse Cuthill McKee reordering.  Default level.
// 2: Use FNROOT to generate pseudo-peripheral nodes, but with trial
//    roots in steps of 5% of _nsolve, up to _nsolve.  Choose root to
//    minimize global bandwidth.
// 3: Do not use FNROOT.  Try all unknown node numbers as trial roots for RCM.
//    Choose root to minimize global bandwidth.
//
// Reference:
//    A. George and J. W-H. Liu
//    Computer Solution of Large Sparse Positive Definite Systems
//    Prentice-Hall (1981)
//
// Note: all the SPARSEPAK routines are in F77, and assume 1-based
// array indices.  Some minor adjustments are made to the numeberings
// in order to accommodate this difference from our standard base-0
// indexing.
// ---------------------------------------------------------------------------
{
  if (!_optlev || !_nsolve) return;

  const int_t verbose = Femlib::ivalue ("VERBOSE");

  VERBOSE cout << "RCM bandwidth optimisation level: " << _optlev ;

  int_t         i, root, nlvl;
  vector<int_t> adjncy, xadj;
  const int_t   tabSize = this -> buildAdjncySC (adjncy, xadj, 1); // Base-1.

  vector<int_t> work(3 * _nsolve + _nglobal + _nbndry);
  int_t         *perm, *mask, *xls, *invperm, *bsave;

  perm    = &work[0];
  mask    = perm    + _nsolve;
  xls     = mask    + _nsolve;
  invperm = xls     + _nsolve;
  bsave   = invperm + _nglobal;
  
  Veclib::copy (_nbndry, &_btog[0], 1, bsave, 1);
  for (i = _nsolve; i < _nglobal; i++) invperm[i] = i; // -- Dirichlet nodes.

  switch (_optlev) {
  case 1: {
    root = 1;
    Veclib::fill   (_nsolve, 1, mask, 1);
    Femlib::fnroot (root, &xadj[0], &adjncy[0], mask, nlvl, xls, perm);
    Femlib::rcm    (root, &xadj[0], &adjncy[0], mask, perm, nlvl, xls);
    break;
  }
  case 2: {
    int_t rtest, BWtest, BWmin = INT_MAX, best, incr = _nsolve / 20;

    VERBOSE cout << ":" << endl;

    for (root = 1; root <= _nsolve; root += incr) {
      rtest = root;
      Veclib::fill   (_nsolve, 1, mask, 1);
      Femlib::fnroot (rtest, &xadj[0], &adjncy[0], mask, nlvl, xls, perm);
      Femlib::rcm    (rtest, &xadj[0], &adjncy[0], mask, perm, nlvl, xls);

      Veclib::sadd (_nsolve, -1, perm, 1, perm, 1); // -- Revert to base-0.
      for (i = 0; i < _nsolve; i++) invperm[perm[i]] = i;
      Veclib::gathr (_nbndry, invperm, bsave, &_btog[0]);

      BWtest = this -> globalBandwidth();
      if (BWtest < BWmin) {
	BWmin = BWtest;
	best  = rtest;
	VERBOSE cout << "root = " << root << ", BW = " << BWmin << endl;
      }
    }

    Veclib::fill (_nsolve, 1, mask, 1);
    Femlib::rcm  (best, &xadj[0], &adjncy[0], mask, perm, nlvl, xls);

    break;
  }
  case 3: {
    int_t BWtest, BWmin = INT_MAX, best;

    VERBOSE cout << ":" << endl;

    for (root = 1; root <= _nsolve; root++) {
      Veclib::fill (_nsolve, 1, mask, 1);
      Femlib::rcm  (root, &xadj[0], &adjncy[0], mask, perm, nlvl, xls);

      Veclib::sadd (_nsolve, -1, perm, 1, perm, 1); // -- Revert to base-0.
      for (i = 0; i < _nsolve; i++) invperm[perm[i]] = i;
      Veclib::gathr (_nbndry, invperm, bsave, &_btog[0]);

      BWtest = this -> globalBandwidth();
      if (BWtest < BWmin) {
	BWmin = BWtest;
	best  = root;
	VERBOSE cout << "root = " << root << ", BW = " << BWmin << endl;
      }
    }

    Veclib::fill (_nsolve, 1, mask, 1);
    Femlib::rcm  (best, &xadj[0], &adjncy[0], mask, perm, nlvl, xls);

    break;
  }
  default:
    break;
  }

  Veclib::sadd (_nsolve, -1, perm, 1, perm, 1); // -- Revert to base-0.
  for (i = 0; i < _nsolve; i++) invperm[perm[i]] = i;
  Veclib::gathr (_nbndry, invperm, bsave, &_btog[0]);

  VERBOSE cout << endl;
}


int_t AssemblyMap::globalBandwidth () const
// --------------------------------------------------------------------------
// Precompute the bandwidth of assembled global matrix (including
// diagonal).
// --------------------------------------------------------------------------
{
  int_t k, noff, nband = 0;
  const int_t    next = _nbndry / _nel;

  for (k = 0, noff = 0; k < _nel; k++) {
    nband = max (this -> bandwidthSC
		 (&_btog[0]+noff, &_bmask[0]+noff, next), nband);
    noff += next;
  }

  ++nband; // -- Diagonal.

  return nband;
}


int_t AssemblyMap::bandwidthSC (const int_t* bmap,
				const int_t* mask,
				const int_t  next) const
// ---------------------------------------------------------------------------
// Find the global equation bandwidth of this element, excluding
// diagonal.
// ---------------------------------------------------------------------------
{
  int_t i, Min = INT_MAX, Max = INT_MIN;

  for (i = 0; i < next; i++) {
    if (!mask[i]) {
      Min = min (bmap[i], Min);
      Max = max (bmap[i], Max);
    }
  }

  return Max - Min;
}

