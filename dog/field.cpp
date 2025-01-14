///////////////////////////////////////////////////////////////////////////////
// field.cpp: derived from AuxField, Field adds boundary conditions,
// global numbering, and the ability to solve Helmholtz problems.
//
// HELMHOLTZ PROBLEMS
// ------------------
// Solve routines provide solution to the discrete form of the Helmholtz eqn
//                      
//                       div grad u - \lambda^2 u = f,
//
// on domain \Omega, subject to essential BCs u = g on \Gamma_g and
// natural BCs \partial u / \partial n = h on \Gamma_h, where the
// boundary \Gamma of \Omega is the union of (non-overlapping)
// \Gamma_g and \Gamma_h and n is the unit outward normal vector on
// \Gamma.  \lambda^2 is called the Helmholtz constant below.
//
// The Galerkin form, using integration by parts with weighting functions w
// which are zero on \Gamma_g, is
//
//            (grad u, grad w) + \lambda^2 (u, w) = - (f, w) + <h, w>
//
// where (a, b) = \int a.b d\Omega is an integration over the domain and
//       <a, b> = \int a.b d\Gamma is an integration over the domain boundary.
//
// The discrete (finite element) equivalent is to solve
//
//                   K.u + \lambda^2 M.u = - M.f + <h, w>
//
// or
//
//                         H.u = - M.f + <h, w>
//
// where K, M and H are respectively (assembled) "stiffness", "mass"
// and Helmholtz matrices.
//
// Some complications arise from dealing with essential boundary
// conditions, since typically the elemental matrices K^e, M^e which
// are assembled to form K and M do not account for the boundary
// requirements on w.  There are a number of ways of dealing with this
// issue: one approach is to partition H as it is formed (here F =
// -M.f + <h, w>):
//
//   +--------+-------------+ /  \     /  \
//   |        |             | |  |     |  |
//   |   Hp   |     Hc      | |u |     |F |   u : nodal values for solution.
//   |        |(constraint) | | s|     | s|    s
//   |        |             | |  |     |  |       (n_solve values)
//   +--------+-------------+ +--+     +--+
//   |        | H_ess: this | |  |  =  |  |
//   |        | partition   | |  |     |  |
//   |    T   | relates to  | |u |     |F |   u : are given essential BCs.
//   |  Hc    | essential   | | g|     | g|    g
//   |        | BC nodes    | |  |     |  |       (n_global - n_solve values)
//   |        | and is not  | |  |     |  |
//   |        | assembled.  | |  |     |  |
//   +--------+-------------+ \  /     \  /
//
// Partition out the sections of the matrix corresponding to the known
// nodal values (essential BCs), and solve instead the constrained
// problem
//
//   +--------+               /  \     /  \     +-------------+ /  \
//   |        |               |  |     |  |     |             | |  |
//   |   Hp   |               |u |     |F |     |     Hc      | |  |
//   |        |               | s|  =  | s|  -  |             | |u |
//   |        |               |  |     |  |     |             | | g|
//   +--------+               \  /     \  /     +-------------+ |  |.
//                                                              |  |
//                                                              |  |
//                                                              \  /
//
// Here n_global is the number of nodes that receive global node
// numbers, typically those on the mesh edges.  N_solve is the number
// of these nodes that have values that must be solved for,
// i.e. n_global minus the number of global nodes situated on
// essential-type boundaries.
//
// An alternative approach (USED HERE) to the constrained problem is to let
//
//                      u = v + g     (v = 0 on \Gamma_g)
//
// and solve instead
//
//                      H.v = - M.f - H.g + <h, w>
//
// (where only the partition Hp is needed for the matrix H in the
// LHS), then afterwards compose u = v + g to get the full solution.
// The advantage to this method is that the constraint partition Hc
// does not need to be assembled or stored.  The operations M.f and
// H.g can be performed on an element-by-element basis.
//
// FIELD NAMES
// -----------
// The (one character) names of field variables are significant, and have
// the following reserved meanings:
// 
// u:  First velocity component.            (Cylindrical: axial     velocity.)
// v:  Second velocity component.           (Cylindrical: radial    velocity.)
// w:  Third velocity component.            (Cylindrical: azimuthal velocity.)
// p:  Pressure divided by density.
// c:  Scalar for transport or elliptic problems.
//
// Copyright (c) 1994+, Hugh M Blackburn
///////////////////////////////////////////////////////////////////////////////

#include <sem.h>


Field::Field (real_t*           M,
	      BoundarySys*      B,
	      NumberSys*        N,
	      const int_t       nz,
	      vector<Element*>& E,
	      const char        C) :
// ---------------------------------------------------------------------------
// Create storage for a new Field from scratch.
// ---------------------------------------------------------------------------
  AuxField (M, nz, E, C),
  _bsys    (B),
  _nsys    (N)
{
  const int_t              np  = Geometry::nP();
  const int_t              npr = Geometry::nProc();
  const int_t              nzb = Geometry::basePlane();
  const vector<Boundary*>& BC  = _bsys -> getBCs (0);
  real_t*                  p;
  int_t                    i, k;

  // -- Allocate storage for boundary data, round up for Fourier
  //    transform, which is never actually done, but we leave it this
  //    way to conform with the basis code, semtex.
  
  _nbound = _bsys -> nSurf();
  _nline  = _nbound * np;
  if   (npr > 1) _nline += 2 * npr - _nline % (2 * npr);
  else           _nline += _nline % 2;

  _line  = new real_t* [static_cast<size_t>(_nz)];
  _sheet = new real_t  [static_cast<size_t>(_nz * _nline)];

  for (k = 0; k < _nz; k++) _line[k] = _sheet + k*_nline;

  Veclib::zero (_nz * _nline, _sheet, 1);

  // -- Set values for boundary data (in Fourier space), and enforce z = 0
  //    (only because we potentially need to give z a value).

  for (k = 0; k < _nz; k++) {
    Femlib::value ("z", 0.0);
    for (p = _line[k], i = 0; i < _nbound; i++, p += np)
      BC[i] -> evaluate (NULL, k, 0, false, p);
  }

  // -- Do NOT Fourier transform boundary data; already in Fourier space.
}


void Field::evaluateBoundaries (const Field* P      ,
				const int_t  step   ,
				const bool   Fourier)
// ---------------------------------------------------------------------------
// Traverse Boundaries and evaluate according to kind.
// Note that for 3D this evaluation is done in Fourier-transformed space.
//
// This routine is only meant to be used for boundary conditions that must
// be re-evaluated at every step, such as high-order pressure BCs.
// ---------------------------------------------------------------------------
{
  const int_t       np = Geometry::nP();
  vector<Boundary*> BC;
  real_t*           p;
  int_t             i, k;

  if (Geometry::nPert() == 2) BC = _bsys -> getBCs (0);
  else                        BC = _bsys -> getBCs (Femlib::ivalue ("BETA"));

  for (k = 0; k < _nz; k++)
    for (p = _line[k], i = 0; i < _nbound; i++, p += np)
      BC[i] -> evaluate (P, k, step, Fourier, p);
}


void Field::evaluateM0Boundaries (const Field* P   ,
				  const int_t  step)
// ---------------------------------------------------------------------------
// Traverse Boundaries and evaluate according to kind, but only for Mode 0.
// ---------------------------------------------------------------------------
{
  const vector<Boundary*>& BC = _bsys -> getBCs (0);
  const int_t              np = Geometry::nP();
  real_t*                  p;
  int_t                    i;

  for (p = _line[0], i = 0; i < _nbound; i++, p += np)
    BC[i] -> evaluate (P, 0, step, false, p);
}


Field& Field::solve (AuxField*        f,
		     const MatrixSys* M, AuxField*        varkinvis)
// ---------------------------------------------------------------------------
// Problem for solution is
//                                          
//                      div grad u - lambda^2 u = f,
//
// which is set up in discrete form as
//
//                       H v = - M f - H g + <h, w>.
//
// This routine creates the RHS vector from the input forcing field f
// and the Field's boundary conditions g (essential) & h (natural).
// Forcing field f's data area is overwritten/destroyed during
// processing.
//
// For DIRECT (Cholesky) solution:
//
//   The RHS vector is constructed with length of the number of
//   element-edge nodes in the problem (n_gid).  The first n_solve
//   values contain forcing terms for the free (non essential-BC)
//   nodes in the problem, derived from the forcing field and the
//   natural BCs "h", while the remaining values get loaded from
//   essential BC values, "g".
//
// For JACPCG (iterative) solution:
//
//   All vectors are ordered with globally-numbered (element-
//   boundary) nodes first, followed by all element-internal nodes.
//   The zeroing operation which occurs after each application of the
//   Helmholtz operator serves to apply the essential BCs, which are
//   zero during the iteration (see file header).
//
//   The notation follows that used in Fig 2.5 of Barrett et al.,
//   "Templates for the Solution of Linear Systems", netlib.
// ---------------------------------------------------------------------------
{
  const char  routine[] = "Field::solve";
  const int_t nel    = Geometry::nElmt();
  const int_t next   = Geometry::nExtElmt();
  const int_t npnp   = Geometry::nTotElmt();
  const int_t ntot   = Geometry::nPlane();
  int_t       i, k;

  const vector<Boundary*>& B       = M -> _BC;
  const AssemblyMap*       A       = M -> _AM;
  real_t                   lambda2 = M -> _HelmholtzConstant;
  real_t                   betak2  = M -> _FourierConstant;
  int_t                    nsolve  = M -> _nsolve;
  int_t                    nglobal = M -> _nglobal;
  int_t                    nzero   = nglobal - nsolve;

  for (k = 0; k < _nz; k++) {
    real_t* forcing = f -> _plane[k];
    real_t* unknown = _plane     [k];
    real_t* bc      = _line      [k];

    switch (M -> _method) {
    
    case DIRECT: {
      const real_t*  H     = const_cast<const real_t*> (M -> _H);
      const real_t** hii   = const_cast<const real_t**>(M -> _hii);
      const real_t** hbi   = const_cast<const real_t**>(M -> _hbi);
      const int_t*   b2g   = const_cast<const int_t*>  (A -> btog());
      int_t          nband = M -> _nband;

      vector<real_t> work (nglobal + 4*npnp);
      real_t         *RHS = &work[0], *tmp = RHS + nglobal;
      int_t          info;
      
      // -- Build RHS = - M f - H g + <h, w>.

      Veclib::zero (nglobal, RHS, 1);
      
      this -> getEssential (bc, RHS, B, A);
      this -> constrain    (forcing, lambda2, varkinvis, betak2,RHS,A,tmp);
      this -> buildRHS     (forcing, bc,RHS,0,hbi,nsolve,nzero,B,A,tmp);
      
      // -- Solve for unknown global-node values (if any).
    
      if (nsolve) Lapack::pbtrs("U",nsolve,nband-1,1,H,nband,RHS,nglobal,info);
      
      // -- Carry out Schur-complement solution for element-internal nodes.
      
      for (i = 0; i < nel; i++, b2g += next, forcing += npnp, unknown += npnp)
	_elmt[i] -> global2localSC (RHS,b2g,forcing,unknown,hbi[i],hii[i],tmp);

      unknown -= ntot;
    
      // -- Scatter-gather essential BC values into plane.

      Veclib::zero (nglobal, RHS, 1);
    
      this -> getEssential (bc, RHS, B,   A);
      this -> setEssential (RHS, unknown, A);
    }
    break;

    case JACPCG: {
      const int_t    StepMax = Femlib::ivalue ("STEP_MAX");
      const int_t    npts    = M -> _npts;
      real_t         alpha, beta, dotp, epsb2, r2, rho1, rho2 = 0.0;
      vector<real_t> work (5*npts + 4*Geometry::nTotElmt());

      const int_t mode =
	(Geometry::problem() == Geometry::O2_2D ||
	 Geometry::problem() == Geometry::SO2_2D ) ? 0 : 1 ;

      real_t* r   = &work[0];
      real_t* p   = r + npts;
      real_t* q   = p + npts;
      real_t* x   = q + npts;
      real_t* z   = x + npts;
      real_t* wrk = z + npts;

      Veclib::zero (nglobal, x, 1);

      this -> getEssential (bc, x, B, A);  
      this -> constrain    (forcing, lambda2,varkinvis,betak2,x,A,wrk);
      this -> buildRHS     (forcing, bc,r,r+nglobal,0,nsolve,nzero,B,A,wrk);

      epsb2  = Femlib::value ("TOL_REL") * sqrt (Blas::dot (npts, r, 1, r, 1));
      epsb2 *= epsb2;

      // -- Build globally-numbered x from element store.

      this -> local2global (unknown, x, A);
  
      // -- Compute first residual using initial guess: r = b - Ax.
      //    And mask to get residual for the zero-BC problem.

      Veclib::zero (nzero, x + nsolve, 1);   
      Veclib::copy (npts,  x, 1, q, 1);
    
      this -> HelmholtzOperator (q, p, lambda2, varkinvis, betak2, mode, wrk);

      Veclib::zero (nzero, p + nsolve, 1);
      Veclib::zero (nzero, r + nsolve, 1);
      Veclib::vsub (npts, r, 1, p, 1, r, 1);

      r2 = Blas::dot (npts, r, 1, r, 1);

      // -- PCG iteration.

      i = 0;
      while (r2 > epsb2 && ++i < StepMax) {
      
	// -- Preconditioner.

	Veclib::vmul (npts, M -> _PC, 1, r, 1, z, 1);
	
	rho1 = Blas::dot (npts, r, 1, z, 1);

	// -- Update search direction.

	if (i == 1)
	  Veclib::copy  (npts,             z, 1, p, 1); // -- p = z.
	else {
	  beta = rho1 / rho2;	
	  Veclib::svtvp (npts, beta, p, 1, z, 1, p, 1); // -- p = z + beta p.
	}

	// -- Matrix-vector product.

	this -> HelmholtzOperator (p, q, lambda2, varkinvis, betak2, mode, wrk);
	Veclib::zero (nzero, q + nsolve, 1);

	// -- Move in conjugate direction.

	dotp  = Blas::dot (npts, p, 1, q, 1);
	alpha = rho1 / dotp;
	Blas::axpy (npts,  alpha, p, 1, x, 1); // -- x += alpha p.
	Blas::axpy (npts, -alpha, q, 1, r, 1); // -- r -= alpha q.
      
	rho2 = rho1;
	r2   = Blas::dot (npts, r, 1, r, 1);
      }
  
      if (i == StepMax) Veclib::alert (routine, "step limit exceeded", WARNING);
    
      // -- Unpack converged vector x, impose current essential BCs.

      this -> global2local (x, unknown, A);
      
      this -> getEssential (bc, x, B,   A);
      this -> setEssential (x, unknown, A);
  
      if (Femlib::ivalue ("VERBOSE") > 1) {
	char s[StrMax];
	sprintf (s, ":%3d iterations, field '%c'", i, _name);
	Veclib::alert (routine, s, REMARK);
      }
    }
    break;
  
    default:
      Veclib::alert
	(routine, "called with a method that isn't implemented", ERROR);
      break;
    }
  }
  return *this; 
}


void Field::constrain (real_t*            force  ,
		       const real_t       lambda2,
               AuxField* VARKINVIS,
 		       const real_t       betak2 ,
		       const real_t*      esstlbc,
		       const AssemblyMap* N      ,
		       real_t*            work   ) const
// ---------------------------------------------------------------------------
// Replace f's data with constrained weak form of forcing: - M f - H g.
// On input, essential BC values (g) have been loaded into globally-numbered
// esstlbc, other values are zero.
//
// Input vector work should be 4*Geometry::nTotElmt() long.
// ---------------------------------------------------------------------------
{
  const int_t       np    = Geometry::nP();
  const int_t       nel   = Geometry::nElmt();
  const int_t       next  = Geometry::nExtElmt();
  const int_t       npnp  = Geometry::nTotElmt();
  const int_t       ntot  = Geometry::nPlane();
  const int_t*      emask = N -> emask();
  const int_t*      btog  = N -> btog();
   Element* E;
  real_t            *u = work, *tmp = work + npnp;
  int_t             i;

  // -- Manufacture -(M f + H g).

  for (i = 0; i < nel; i++, emask++, btog += next, force += npnp) {
    E = _elmt[i];
    E -> weight (force);	// -- f <-- M f.
    if (*emask) {		// -- f <-- M f + H g.
      Veclib::zero      (npnp, u, 1);
      E -> global2local (u, btog, esstlbc, 0);
      int_t offset = E-> ID() * npnp;
      E -> HelmholtzOp  (lambda2, (VARKINVIS->getData())+offset, betak2, u, u, tmp);
      Veclib::vadd      (npnp, force, 1, u, 1, force, 1);
    }
  }

  force -= ntot;
  Veclib::neg (ntot, force, 1);
}


void Field::HelmholtzOperator (const real_t* x      ,
			       real_t*       y      ,
			       const real_t  lambda2,
                   AuxField* VARKINVIS,
			       const real_t  betak2 ,
			       const int_t   mode   ,
			       real_t*       work   ) const
// ---------------------------------------------------------------------------
// Discrete 2D global Helmholtz operator which takes the vector x into
// vector y, including direct stiffness summation.  Vectors x & y have 
// global ordering: that is, with nglobal (element edge nodes, with
// redundancy removed) coming first, followed by nel blocks of element-
// internal nodes.
//
// Vector work must have length 4*Geometry::nTotElmt().
// ---------------------------------------------------------------------------
{
  const int_t        np      = Geometry::nP();
  const int_t        nel     = Geometry::nElmt();
  const int_t        npnp    = Geometry::nTotElmt();
  const int_t        next    = Geometry::nExtElmt();
  const int_t        nint    = Geometry::nIntElmt();
  const int_t        ntot    = Geometry::nPlane();
  const AssemblyMap* AM      = _nsys -> getMap (mode * Femlib::ivalue ("BETA"));
  const int_t*       gid     = AM -> btog();
  const int_t        nglobal = AM -> nGlobal() + Geometry::nInode();
  const real_t*      xint    = x + AM -> nGlobal();
  real_t*            yint    = y + AM -> nGlobal();
  real_t             *P = work, *tmp = work + npnp;
  int_t              i;
  Element*           E;

  Veclib::zero (nglobal, y, 1);

  // -- Add in contributions from mixed BCs while x & y are global vectors.

  if (_bsys -> mixBC()) {
    const vector<Boundary*>& BC   = _bsys -> getBCs (0);
    const int_t*             bmap = AM    -> btog();

    for (i = 0; i < _nbound; i++)
      BC[i] -> augmentOp (bmap, x, y);
  }

  // -- Add in contributions from elemental Helmholtz operations.

  for (i = 0; i < nel; i++, gid += next, xint += nint, yint += nint) {
    E = _elmt[i];
    E -> global2local    (P, gid, x, xint);
    int_t offset = E-> ID() * npnp;
    E -> HelmholtzOp     (lambda2, VARKINVIS -> getData()+offset, betak2, P, P, tmp);
    E -> local2globalSum (P, gid, y, yint);
  }
}


void Field::buildRHS (real_t*                  force ,
		      const real_t*            bc    ,
		      real_t*                  RHS   ,
		      real_t*                  RHSint,
		      const real_t**           hbi   ,
		      const int_t              nsolve,
		      const int_t              nzero ,
		      const vector<Boundary*>& bnd   ,
		      const AssemblyMap*       Assy  ,
		      real_t*                  work  ) const
// ---------------------------------------------------------------------------
// Build RHS for direct or iterative solution.
//
// Iterative solution is flagged by presence of RHSint, a pointer to
// element-internal node storage.  If RHSint is zero, then hbi, a vector
// of pointers to element interior/exterior coupling matrices, must be 
// non-zero.
//
// Compute RHS vector for direct solution of Helmholtz problem as
//
//                      - M f - H g + <h, w>.
// 
// On input, force contains a plane of
//
//                          - M f - H g
//
// in element (row-major) form, and bc contains the line of BC data values
// for this plane of data: only natural BCs are used in formation of <h, w>.
//
// Input vector work should be Geometry::nTotElmt() long.
// ---------------------------------------------------------------------------
{
  const int_t              np      = Geometry::nP();
  const int_t              nel     = Geometry::nElmt();
  const int_t              next    = Geometry::nExtElmt();
  const int_t              nint    = Geometry::nIntElmt();
  const int_t              npnp    = Geometry::nTotElmt();
  const int_t              nglobal = Assy -> nGlobal();
  const int_t*             gid;
   const Boundary* B;
  int_t                    i, boff;

  if   (RHSint) Veclib::zero (nglobal + Geometry::nInode(), RHS, 1);
  else          Veclib::zero (nglobal,                      RHS, 1);

  // -- Add in contribution from forcing f = - M f - H g.

  for (gid = Assy -> btog(), i = 0; i < nel; i++, force += npnp, gid += next) {
    if (RHSint) {
      _elmt[i] -> local2globalSum   (force, gid, RHS, RHSint); RHSint += nint;
    } else
      _elmt[i] -> local2globalSumSC (force, gid, RHS, hbi[i], work);
  }

  // -- Add in <h, w>.

  for (gid = Assy -> btog(), i = 0; i < _nbound; i++, bc += np) {
    B    = bnd[i];
    boff = B -> bOff();

    B -> sum (bc, gid + boff, work, RHS);
  }

  // -- Zero any contribution that <h, w> made to essential BC nodes.

  Veclib::zero (nzero, RHS + nsolve, 1);
}


void Field::local2global (const real_t*      src,
			  real_t*            tgt,
			  const AssemblyMap* A  ) const
// ---------------------------------------------------------------------------
// Load a plane of data (src) into globally-numbered tgt, with element-
// boundary values in the first A -> nGlobal() places, followed by
// element-internal locations in emap ordering.
// ---------------------------------------------------------------------------
{
  const int_t  nel  = Geometry::nElmt();
  const int_t  next = Geometry::nExtElmt();
  const int_t  nint = Geometry::nIntElmt();
  const int_t  npnp = Geometry::nTotElmt();
  const int_t* gid  = A -> btog();
  int_t        i;
  real_t*      internal = tgt + A -> nGlobal();

  for (i = 0; i < nel; i++, src += npnp, gid += next, internal += nint)
    _elmt[i] -> local2global (src, gid, tgt, internal);
}


void Field::global2local (const real_t*      src,
			  real_t*            tgt,
			  const AssemblyMap* A  ) const
// ---------------------------------------------------------------------------
// Load a plane of data (tgt) from src, which has globally-numbered
// (element- boundary) values in the first nglobal places, followed by
// element-internal locations in emap ordering.
// ---------------------------------------------------------------------------
{
  const int_t   nel  = Geometry::nElmt();
  const int_t   next = Geometry::nExtElmt();
  const int_t   nint = Geometry::nIntElmt();
  const int_t   npnp = Geometry::nTotElmt();
  const int_t*  gid  = A -> btog();
  int_t         i;
  const real_t* internal = src + A -> nGlobal();

  for (i = 0; i < nel; i++, tgt += npnp, gid += next, internal += nint)
    _elmt[i] -> global2local (tgt, gid, src, internal);
}


void Field::local2globalSum (const real_t*      src,
			     real_t*            tgt,
			     const AssemblyMap* A  ) const
/// --------------------------------------------------------------------------
/// Direct stiffness sum a plane of data (src) into globally-numbered
/// tgt, with element-boundary values in the first A -> nGlobal()
/// places, followed by element-internal locations in emap ordering.
// ---------------------------------------------------------------------------
{
  const int_t  nel  = Geometry::nElmt();
  const int_t  next = Geometry::nExtElmt();
  const int_t  nint = Geometry::nIntElmt();
  const int_t  npnp = Geometry::nTotElmt();
  const int_t* gid  = A -> btog();
  int_t        i;
  real_t*      internal = tgt + A -> nGlobal();

  for (i = 0; i < nel; i++, src += npnp, gid += next, internal += nint)
    _elmt[i] -> local2globalSum (src, gid, tgt, internal);
}


void Field::getEssential (const real_t*              src,
			  real_t*                    tgt,
			  const vector<Boundary*>&   bnd,
			  const AssemblyMap*         A  ) const
// ---------------------------------------------------------------------------
// On input, src contains a line of BC values for the current data plane.
// Scatter current values of essential BCs into globally-numbered tgt.
//
// The construction of the essential BCs has to account for cases
// where element corners may touch the domain boundary but do not have
// an edge along a boundary.  This is done by working with a
// globally-numbered vector.
// ---------------------------------------------------------------------------
{
  const int_t     np   = Geometry::nP();
  const int_t*    btog = A -> btog();
  const Boundary* B;
  int_t           i, boff;
  
  for (i = 0; i < _nbound; i++, src += np) {
    B    = bnd[i];
    boff = B -> bOff();
  
    B -> set (src, btog + boff, tgt);
  }
}


void Field::setEssential (const real_t*      src,
			  real_t*            tgt,
			  const AssemblyMap* A  )
// ---------------------------------------------------------------------------
// Gather globally-numbered src into essential BC nodes of current
// data plane.
// ---------------------------------------------------------------------------
{
  const int_t  nel  = Geometry::nElmt();
  const int_t  next = Geometry::nExtElmt();
  const int_t  npnp = Geometry::nTotElmt();
  const int_t* emask = A -> emask();
  const int_t* bmask = A -> bmask();
  const int_t* btog  = A -> btog();
  int_t        i;

  for (i = 0; i < nel; i++, bmask += next, btog += next, tgt += npnp)
    if (emask[i]) _elmt[i] -> bndryMask (bmask, tgt, src, btog);
}


void Field::coupleBCs (Field*      v  ,
		       Field*      w  ,
		       const int_t dir)
// ---------------------------------------------------------------------------
// Couples/uncouple boundary condition values for the radial and
// azimuthal velocity fields in cylindrical coordinates, depending on
// indicated direction.  This action is required due to the coupling
// in the viscous terms of the N--S equations in cylindrical coords.
//
// dir == +1
// ---------
//           v~ <-- v + i w
//           w~ <-- v - i w
// dir == -1
// ---------
//           v  <-- 0.5   * (v~ + w~)
//           w  <-- 0.5 i * (w~ - v~)
//
// Since there is no coupling for the viscous terms in the 2D
// equation, do nothing for the zeroth Fourier mode.
// --------------
{
  if (Geometry::problem() == Geometry::O2_2D ||
      Geometry::problem() == Geometry::SO2_2D ) return;

  const char     routine[] = "Field::couple";
  const int_t    nL        =  v -> _nline;
  vector<real_t> work (nL);
  real_t         *Vr, *Vi, *Wr, *Wi, *tp = &work[0];
  
  if (dir == FORWARD) {

    if (v->_nz == 1) {		// -- Half-complex.
      Vr = v -> _line[0];
      Wi = w -> _line[0];

      Veclib::copy (nL, Vr, 1, tp, 1);
      Veclib::vsub (nL, Vr, 1, Wi, 1, Vr, 1);
      Veclib::vadd (nL, Wi, 1, tp, 1, Wi, 1);
    } else {			// -- Full complex.
      Vr = v -> _line[0];
      Vi = v -> _line[1];
      Wr = w -> _line[0];
      Wi = w -> _line[1];

      Veclib::copy (nL, Vr, 1, tp, 1);
      Veclib::vsub (nL, Vr, 1, Wi, 1, Vr, 1);
      Veclib::vadd (nL, Wi, 1, tp, 1, Wi, 1);
      Veclib::copy (nL, Wr, 1, tp, 1);
      Veclib::copy (nL, Wi, 1, Wr, 1);
      Veclib::vsub (nL, Vi, 1, tp, 1, Wi, 1);
      Veclib::vadd (nL, Vi, 1, tp, 1, Vi, 1);
    }

  } else if (dir == INVERSE) {

    if (v->_nz == 1) {		// -- Half-complex.
      Vr = v -> _line[0];
      Wr = w -> _line[0];

      Veclib::copy  (nL,      Vr, 1, tp, 1);
      Veclib::svvpt (nL, 0.5, Vr, 1, Wr, 1, Vr, 1);
      Veclib::svvmt (nL, 0.5, Wr, 1, tp, 1, Wr, 1);
    } else {			// -- Full complex.
      Vr = v -> _line[0];
      Vi = v -> _line[1];
      Wr = w -> _line[0];
      Wi = w -> _line[1];

      Veclib::copy  (nL,      Vr, 1, tp, 1);
      Veclib::svvpt (nL, 0.5, Vr, 1, Wr, 1, Vr, 1);
      Veclib::svvmt (nL, 0.5, Wr, 1, tp, 1, Wr, 1);
      Veclib::copy  (nL,      Wi, 1, tp, 1);
      Veclib::copy  (nL,      Wr, 1, Wi, 1);
      Veclib::svvpt (nL, 0.5, Vi, 1, tp, 1, Wr, 1);
      Veclib::svvmt (nL, 0.5, Vi, 1, tp, 1, Vi, 1);
    }

  } else
    Veclib::alert (routine, "unknown direction given", ERROR);
}


real_t Field::modeConstant (const char   name,
			    const int_t  mode,
			    const real_t beta)
// ---------------------------------------------------------------------------
// For cylindrical coordinates & 3D, the radial and azimuthal fields
// are coupled before solution of the viscous step.  This means that
// the Fourier constant used for solution may vary from that which
// applies to the axial component.
//
// For Field v~, betak -> betak + 1 while for w~, betak -> betak - 1.
//
// For the uncoupled Fields v, w solved for the zeroth Fourier mode,
// the "Fourier" constant in the Helmholtz equations is +/-1.
// ---------------------------------------------------------------------------
{
  const char routine[] = "Field::modeConstant";

  if (Geometry::problem() == Geometry::O2_2D     ||
      Geometry::problem() == Geometry::SO2_2D    ||
      Geometry::system()  == Geometry::Cartesian || 
      name                ==         'c'         ||
      name                ==         'p'         ||
      name                ==         'u'          ) return beta * mode;

  if      (name == 'v') return beta * mode + 1.0;
  else if (name == 'w') return beta * mode - 1.0;
  else Veclib::alert (routine, "unrecognized Field name given", ERROR);

  return -1.0;			// -- Never happen.
}
