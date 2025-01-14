///////////////////////////////////////////////////////////////////////////////
// field.cpp: Derived from AuxField, Field adds boundary conditions,
// global numbering, and the ability to solve Helmholtz problems.
// See also field.h.
//
// Copyright (c) 1994+, Hugh M Blackburn
///////////////////////////////////////////////////////////////////////////////

#include <sem.h>


Field::Field (real_t*           M ,   // -- Data storage allocation.
	      BoundarySys*      B ,   // -- BC info for all Fourier modes.
	      NumberSys*        N ,   // -- Assembly maps for all Fourier modes.
	      const int_t       nz,   // -- No. of planes of data per process.
	      vector<Element*>& E ,   // -- Element operator pointers.
	      const char        C ) : // -- Single-character name (significant).
/// --------------------------------------------------------------------------
/// Create storage for a new Field from scratch.
// ---------------------------------------------------------------------------
  AuxField (M, nz, E, C),
  _bsys    (B),
  _nsys    (N)
{
  const int_t              np  = Geometry::nP();
  const int_t              npr = Geometry::nProc();
  const int_t              nzb = Geometry::basePlane();
  const vector<Boundary*>& BC  = _bsys -> getBCs (0);
  const real_t             dz  = Femlib::value ("TWOPI / BETA / N_Z");
  real_t*                  p;
  int_t                    i, k;

  // -- Allocate storage for boundary data, round up size for Fourier transform.
  
  _nbound = _bsys -> nSurf();
  _nline  = _nbound * np;
  if   (npr > 1) _nline += 2 * npr - _nline % (2 * npr);
  else           _nline += _nline % 2;

  _line  = new real_t* [static_cast<size_t>(_nz)];
  _sheet = new real_t  [static_cast<size_t>(_nz * _nline)];

  for (k = 0; k < _nz; k++) _line[k] = _sheet + k*_nline;

  Veclib::zero (_nz * _nline, _sheet, 1);

  // -- Set values for boundary data (in physical space, so no Field
  //    is needed and first argument is void in line with false final
  //    argument).

  this -> evaluateBoundaries (NULL, 0, false);

  // -- Fourier transform boundary data.

  this -> bTransform (FORWARD);
}


Field& Field::solve (AuxField*             f  ,
		     const ModalMatrixSys* MMS)
/// --------------------------------------------------------------------------
/// This is the central routine for which the Field class was created.
//
//  Problem for solution is
//                                          
//                      div grad u - lambda^2 u = f,
//
//  which is set up in discrete form as
//
//                       H v = - M f - H g + <h, w>.
//                                      c
//
//  See also the description at the top of field.h.
//  
/// This routine creates the RHS vector from the input forcing field f
/// and the Field's boundary conditions g (essential) & h (natural).
/// Forcing field f's data area is overwritten/destroyed during
/// processing.
///
/// Field data storage on input is assumed to be the existing
/// estimate. This is used for iterative solutions, and for
/// convective-type BCs.
///
/// For DIRECT (Cholesky) solution:
///
///   The RHS vector is constructed with length of the number of
///   element-edge nodes in the problem (n_gid).  The first n_solve
///   values contain forcing terms for the free (non essential-BC)
///   nodes in the problem, derived from the forcing field and the
///   natural BCs "h", while the remaining values get loaded from
///   essential BC values, "g".
///
/// For JACPCG (iterative) solution:
///
///   All vectors are ordered with globally-numbered (element-
///   boundary) nodes first, followed by all element-internal nodes.
///   The zeroing operation which occurs after each application of the
///   Helmholtz operator serves to apply the essential BCs, which are
///   zero during the iteration.
///
///   The notation under JACPCG follows that used in Fig 2.5 of Barrett
///   et al., "Templates for the Solution of Linear Systems", netlib.
///   Iteration stops when ||r|| = ||Ax - b|| < TOL_REL^2 * ||b|
///   (Criterion 2 in Barrett et al.).
//   ---------------------------------------------------------------------------
{
  const char  routine[] = "Field::solve";
  const int_t np    = Geometry::nP();
  const int_t nel   = Geometry::nElmt();
  const int_t next  = Geometry::nExtElmt();
  const int_t npnp  = Geometry::nTotElmt();
  const int_t ntot  = Geometry::nPlane();
  const int_t bmode = Geometry::baseMode(); // -- Process's lowest mode number.
  int_t       i, k, pmode, mode;

  for (k = 0; k < _nz; k++) {	// -- Loop over planes of data.
    
    ROOTONLY if (k == 1) continue;	// -- Nyquist plane always zero.

    // -- Select Fourier mode, set local pointers and variables.

    pmode = k >> 1;
    mode  = bmode + pmode;

    const MatrixSys*         M       = (*MMS)[pmode];
    const vector<Boundary*>& B       = M -> _BC;
    const AssemblyMap*       A       = M -> _AM;
    real_t                   lambda2 = M -> _HelmholtzConstant;
    real_t                   betak2  = M -> _FourierConstant;
    int_t                    nsolve  = M -> _nsolve;
    int_t                    nglobal = M -> _nglobal;
    int_t                    nzero   = nglobal - nsolve;

    real_t*                  forcing = f -> _plane[k];
    real_t*                  unknown = _plane     [k];
    real_t*                  bc      = _line      [k];

    switch (M -> _method) {

    case DIRECT: {
      const real_t*  H     = const_cast<const real_t*>  (M -> _H);
      const real_t** hii   = const_cast<const real_t**> (M -> _hii);
      const real_t** hbi   = const_cast<const real_t**> (M -> _hbi);
      const int_t*   b2g   = const_cast<const int_t*>   (A -> btog());
      int_t          nband = M -> _nband;

      vector<real_t> work (nglobal + 4*npnp);
      real_t         *RHS = &work[0], *tmp = RHS + nglobal;
      int_t          info;
      
      // -- Build RHS = - M f - H g + <h, w>.

      Veclib::zero (nglobal, RHS, 1);
      AuxField* temp;
      this -> getEssential (bc, RHS, B, A);
      this -> constrain    (forcing, lambda2, temp, betak2, RHS, A, tmp);
      this -> buildRHS     (forcing, bc, RHS, 0, hbi, nsolve, nzero,B,A,tmp);
      
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
      const int_t    StepMax =  Femlib::ivalue ("STEP_MAX");
      const int_t    npts    = M -> _npts;
      real_t         alpha, beta, dotp, epsb2, r2, rho1, rho2;
#if defined (_VECTOR_ARCH)
      vector<real_t> work (5 * npts + 3 * Geometry::nPlane());
#else
      vector<real_t> work (5 * npts + 4 * Geometry::nTotElmt());
#endif
      real_t* r   = &work[0];
      real_t* p   = r + npts;
      real_t* q   = p + npts;
      real_t* x   = q + npts;
      real_t* z   = x + npts;
      real_t* wrk = z + npts;

      Veclib::zero (nglobal, x, 1);
      AuxField* temp;

      this -> getEssential (bc,x,B,A);
      this -> constrain    (forcing,lambda2,temp,betak2,x,A,wrk);
      this -> buildRHS     (forcing,bc,r,r+nglobal,0,nsolve,nzero,B,A,wrk);

      epsb2  = Femlib::value ("TOL_REL") * sqrt (Blas::dot (npts, r, 1, r, 1));
      epsb2 *= epsb2;

      // -- Build globally-numbered x from element store.

      this -> local2global (unknown, x, A);

      // -- Compute first residual using initial guess: r = b - Ax.
      //    And mask to get residual for the zero-BC problem.

      Veclib::zero (nzero, x + nsolve, 1);   
      Veclib::copy (npts,  x, 1, q, 1);

      this -> HelmholtzOperator (q, p, lambda2, betak2, mode, wrk);

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

	this -> HelmholtzOperator (p, q, lambda2, betak2, mode, wrk);

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
  
      if (static_cast<int_t>(Femlib::value ("VERBOSE")) > 1) {
	char s[StrMax];
	sprintf (s, ":%3d iterations, field '%c'", i, _name);
	Veclib::alert (routine, s, REMARK);
      }
    }
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
/// --------------------------------------------------------------------------
/// Replace f's data with constrained weak form of forcing: - M f - H g.
/// On input, essential BC values (g) have been loaded into globally-numbered
/// esstlbc, other values are zero.
///
/// Input vector work should be 4*Geometry::nTotElmt() long.
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
  int_t    i;
  real_t            *u = work, *tmp = work + npnp;

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
			       const real_t  betak2 ,
			       const int_t   mode   ,
			       real_t*       work   ) const
/// --------------------------------------------------------------------------
/// Discrete 2D global Helmholtz operator which takes the vector x into
/// vector y, including direct stiffness summation.  Vectors x & y have 
/// global ordering: that is, with nglobal (element edge nodes, with
/// redundancy removed) coming first, followed by nel blocks of element-
/// internal nodes.
//
#if defined (_VECTOR_ARCH)
// Vector work must have length 3 * Geometry::nPlane().
#else
// Vector work must have length 4 * Geometry::nTotElmt().
#endif
//
/// Note that we multiply the input value of mode by BETA in order to
/// pick up the correct mode-dependent set of (axis) BCs for
/// cylindrical-coordinate problems.
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
  int_t              i;

  Veclib::zero (nglobal, y, 1);

  // -- Add in contributions from mixed BCs while x & y are global vectors.

  if (_bsys -> mixBC()) {
    const vector<Boundary*>& BC = _bsys -> getBCs (0);

    for (i = 0; i < _nbound; i++) BC[i] -> augmentOp (gid, x, y);
  }

  // -- Add in contributions from elemental Helmholtz operations.

#if defined (_VECTOR_ARCH)
  const real_t *DV, *DT;
  real_t       *P = work, *R = P + ntot, *S = R + ntot;

  Femlib::quadrature (0, 0, &DV, 0  , np, GLJ, JAC_ALFA, JAC_BETA);
  Femlib::quadrature (0, 0, 0  , &DT, np, GLJ, JAC_ALFA, JAC_BETA);

  Veclib::zero (ntot + ntot, R, 1);

  this -> global2local (x, P, AM);

  Femlib::grad2 (P, P, R, S, DV, DT, np, np, nel);

  for (i = 0; i < nel; i++, R += npnp, S += npnp, P += npnp)
    _elmt[i] -> HelmholtzKern (lambda2, betak2, R, S, P, P);
 
  P -= ntot;
  R -= ntot;
  S -= ntot;
  
  Femlib::grad2 (R, S, P, P, DT, DV, np, np, nel);

  this -> local2globalSum (P, y, AM);

#else

  const real_t *xint    = x + AM -> nGlobal();
  real_t       *yint    = y + AM -> nGlobal();
  real_t       *P = work, *tmp = work + npnp;
  Element      *E;

  for (i = 0; i < nel; i++, gid += next, xint += nint, yint += nint) {
    E = _elmt[i];
    E -> global2local    (P, gid, x, xint);
    real_t* temp;
    E -> HelmholtzOp     (lambda2, temp, betak2, P, P, tmp);
    E -> local2globalSum (P, gid, y, yint);
  }

#endif
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
/// --------------------------------------------------------------------------
/// Build RHS for direct or iterative solution.
//
/// Iterative solution is flagged by presence of RHSint, a pointer to
/// element-internal node storage.  If RHSint is zero, then hbi, a vector
/// of pointers to element interior/exterior coupling matrices, must be 
/// non-zero.
//
//  Compute RHS vector for direct solution of Helmholtz problem as
//
//                      - M f - H g + <h, w>.
// 
//  On input (and output), force contains a plane of
//
//                          - M f - H g
//
//  in element (row-major) form, and bc contains the line of BC data
//  values for this plane of data: only natural/mixed (Neumann/Robin)
//  BCs are used in formation of <h, w>.
//
/// Input vector work should be Geometry::nTotElmt() long.
// ---------------------------------------------------------------------------
{
  const int_t     np      = Geometry::nP();
  const int_t     nel     = Geometry::nElmt();
  const int_t     next    = Geometry::nExtElmt();
  const int_t     nint    = Geometry::nIntElmt();
  const int_t     npnp    = Geometry::nTotElmt();
  const int_t     nglobal = Assy -> nGlobal();
  const int_t*    gid;
  const Boundary* B;
  int_t           i, boff;

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
/// --------------------------------------------------------------------------
/// Load a plane of data (src) into globally-numbered tgt, with element-
/// boundary values in the first A -> nGlobal() places, followed by
/// element-internal locations in emap ordering.
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
/// --------------------------------------------------------------------------
/// Load a plane of data (tgt) from src, which has globally-numbered
/// (element- boundary) values in the first nglobal places, followed by
/// element-internal locations in emap ordering.
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


void Field::getEssential (const real_t*            src,
			  real_t*                  tgt,
			  const vector<Boundary*>& bnd,
			  const AssemblyMap*       A  ) const
/// --------------------------------------------------------------------------
/// On input, src contains a line of BC values for the current data plane.
/// Scatter current values of essential BCs into globally-numbered tgt.
///
/// The construction of the essential BCs has to account for cases
/// where element corners may touch the domain boundary but do not have
/// an edge along a boundary.  This is done by working with a
/// globally-numbered vector.
// ---------------------------------------------------------------------------
{
  const int_t     np = Geometry::nP();
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
/// --------------------------------------------------------------------------
/// Gather globally-numbered src into essential BC nodes of current
/// data plane (this is the inverse of "lifting" them out of solution).
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


void Field::bTransform (const int_t sign)
/// --------------------------------------------------------------------------
/// Compute forward or backward 1D-DFT of boundary value storage areas.
///
/// Normalization is carried out on forward transform, such that the zeroth
/// mode's real_t data are the average over the homogeneous direction of the
/// physical space values.  See also comments for AuxField::transform.
///
/// The Nyquist data do not need to be zeroed as these planes are never
/// evolved regardless of BC.
// ---------------------------------------------------------------------------
{
  const int_t nZ  = Geometry::nZ();
  const int_t nPR = Geometry::nProc();
  const int_t nP  = _nline;
  const int_t nPP = _nline / nPR;

  if (nPR == 1) {
    if (nZ > 1)
      if (nZ == 2)
	if   (sign == FORWARD) Veclib::zero (_nline, _line[1], 1);
	else                   Veclib::copy (_nline, _line[0], 1, _line[1], 1);
      else
	Femlib::DFTr (_sheet, nZ, _nline, sign);
  } else {
    Message::exchange (_sheet, _nz, nP, FORWARD);
    Femlib::DFTr      (_sheet, nZ, nPP, sign   );
    Message::exchange (_sheet, _nz, nP, INVERSE);
  }
}


void Field::evaluateBoundaries (const Field* P      ,
				const int_t  step   ,
				const bool   Fourier)
/// --------------------------------------------------------------------------
/// Traverse Boundaries and evaluate according to kind.  Note that for
/// 3D this evaluation is done in Fourier-transformed space if Fourier
/// = true (the default).
///
/// This routine is mainly intended to be used for computed boundary
/// conditions that must be re-evaluated at every step, such as
/// high-order pressure BCs or mixed velocity and scalar fields that
/// explicitly vary in time.
///
/// Note that if argument Fourier is false, Field* P is eventually not
/// used by Condition::evaluate(), so that we could alternatively use a
/// test on existence (non-NULLness) of P to detect intention.
/// ---------------------------------------------------------------------------
{
  const int_t  np    = Geometry::nP();
  const int_t  nz    = Geometry::nZProc();
  const int_t  bmode = Geometry::baseMode();
  const int_t  nzb   = Geometry::basePlane();
  const real_t dz    = Femlib::value ("TWOPI / BETA / N_Z");
  real_t*      p;
  int_t        i, k, mode;

  if (Fourier) {
    for (k = 0; k < nz; k++) {
      mode = bmode + (k >> 1);
      const vector<Boundary*>& BC =
	_bsys -> getBCs (mode*Femlib::ivalue("BETA"));
      for (p = _line[k], i = 0; i < _nbound; i++, p += np)
	BC[i] -> evaluate (P, k, step, true, p);
    }
  } else {
    for (k = 0; k < _nz; k++) {
      Femlib::value ("z", (nzb + k) * dz);
      const vector<Boundary*>& BC = _bsys -> getBCs (0);
      for (p = _line[k], i = 0; i < _nbound; i++, p += np)
	BC[i] -> evaluate (P, k, step, false, p);
    }
  }
}


void Field::evaluateM0Boundaries (const Field* P   ,
				  const int_t  step)
/// --------------------------------------------------------------------------
/// Traverse Boundaries and evaluate according to kind, but only for
/// Mode 0, i.e. on the z-average mode.
/// ---------------------------------------------------------------------------
{
  ROOTONLY {
    const vector<Boundary*>& BC = _bsys -> getBCs (0);
    const int_t              np = Geometry::nP();
    real_t*                  p;
    int_t                    i;

    for (p = _line[0], i = 0; i < _nbound; i++, p += np)
      BC[i] -> evaluate (P, 0, step, false, p);
  }
}


void Field::addToM0Boundaries (const real_t val,
			       const char*  grp)
/// --------------------------------------------------------------------------
/// Add val to 0th Fourier mode's bc storage area on BC group "grp".
// ---------------------------------------------------------------------------
{
  ROOTONLY {
    const vector<Boundary*>& BC = _bsys -> getBCs (0);
    const int_t              np = Geometry::nP();
    real_t*                  p;
    int_t           i;

    for (p = _line[0], i = 0; i < _nbound; i++, p += np)
      BC[i] -> addForGroup (grp, val, p);
  }
}


real_t Field::scalarFlux (const Field* C)
/// --------------------------------------------------------------------------
/// Static member function.
///
/// Compute edge-normal gradient flux of field C on all "wall" group
/// boundaries.
///
/// This only has to be done on the zero (mean) Fourier mode.
// ---------------------------------------------------------------------------
{
  const vector<Boundary*>& BC = C -> _bsys -> getBCs (0);
  vector<real_t>           work(4 * Geometry::nP());
  real_t                   F = 0.0;
  int_t           i;
  
  for (i = 0; i < C -> _nbound; i++)
    F += BC[i] -> scalarFlux ("wall", C -> _data, &work[0]);

  return F;
}


Vector Field::normTraction (const Field* P)
/// --------------------------------------------------------------------------
/// Static member function.
///
/// Compute normal tractive forces on all WALL boundaries, taking P to
/// be the pressure field.
///
/// This only has to be done on the zero (mean) Fourier mode.
///
/// Note that the z component of this traction is always zero since
/// no component of the unit outward normal points in the z direction in
/// 2.5D geometries.  
// ---------------------------------------------------------------------------
{
  const vector<Boundary*>& BC = P -> _bsys -> getBCs (0);
  const int_t              nsurf = P -> _nbound;
  Vector                   secF, F = {0.0, 0.0, 0.0};
  vector<real_t>           work(Geometry::nP());
  int_t           i;
  
  for (i = 0; i < nsurf; i++) {
    secF = BC[i] -> normTraction ("wall", P -> _data, &work[0]);
    F.x += secF.x;
    F.y += secF.y;
  }

  return F;
}


Vector Field::tangTraction (const Field* U,
			    const Field* V,
			    const Field* W)
/// --------------------------------------------------------------------------
/// Static member function.
///
/// Compute (2D) tangential viscous tractive forces on all WALL
/// boundaries, treating U & V as first and second velocity
/// components, respectively.
//
//  Compute viscous tractive forces on wall from
//
//  t_i  = - T_ij  n_j       (minus sign for force exerted BY fluid ON wall),
//
//  where
//
//  T_ij = viscous stress tensor (here in Cartesian coords)
//                          dU_i    dU_j
//       = RHO * KINVIS * ( ----  + ---- ) .
//                          dx_j    dx_i
//
/// This only has to be done on the zero (mean) Fourier mode.
// ---------------------------------------------------------------------------
{
  const vector<Boundary*>& UBC =       U->_bsys->getBCs(0);
  const vector<Boundary*>& WBC = (W) ? W->_bsys->getBCs(0) : (vector<Boundary*>)0;
  const int_t              np     = Geometry::nP();
  const int_t              nbound = U -> _nbound;
  const real_t             mu     = Femlib::value ("RHO * KINVIS");
  Vector                   secF, F= {0.0, 0.0, 0.0};
  vector<real_t>           work(4 * np);
  int_t           i;

  for (i = 0; i < nbound; i++) {
    secF = UBC[i] -> tangTraction ("wall", U->_data, V->_data, &work[0]);
    F.x        -= mu * secF.x;
    F.y        -= mu * secF.y;
    if (W) F.z -= mu * WBC[i] -> torqueFlux ("wall", W->_data, &work[0]);
  }

  return F;
}


void Field::normTractionV (real_t*      fx,
			   real_t*      fy,
			   const Field* P )
/// --------------------------------------------------------------------------
/// Static member function.
///
/// Compute normal tractive forces at each z location on wall
/// boundaries.  Fx & Fy are assumed to contain sufficient (nZProc)
/// storage and to be zero on entry.  P could be in physical space or
/// Fourier transformed.
// ---------------------------------------------------------------------------
{
  const vector<Boundary*>& BC    = P -> _bsys -> getBCs (0);
  const int_t              np    = Geometry::nP();
  const int_t              nz    = Geometry::nZProc();
  const int_t              nsurf = P -> _nbound;
  Vector                   secF;
  vector<real_t>           work(np);
  real_t                   *p;
  int_t                    i, j;
  
  for (j = 0; j < nz; j++) {
    p = P -> _plane[j];
    for (i = 0; i < nsurf; i++) {
      secF = BC[i] -> normTraction ("wall", p, &work[0]);
      fx[j] += secF.x;
      fy[j] += secF.y;
    }
  }
}


void Field::tangTractionV (real_t*      fx,
			   real_t*      fy,
			   real_t*      fz,
			   const Field* U ,
			   const Field* V ,
			   const Field* W )
/// --------------------------------------------------------------------------
/// Static member function.
///
/// Compute tangential tractive forces at each z location on wall
/// boundaries.  Fx, fy, fz assumed to contain sufficient storage and
/// be zero on entry.  U, V, W could be in physical space or Fourier
/// transformed.
// ---------------------------------------------------------------------------
{
  const vector<Boundary*>& UBC =       U->_bsys->getBCs(0);
  const vector<Boundary*>& WBC = (W) ? W->_bsys->getBCs(0) : (vector<Boundary*>)0;
  const int_t              np     = Geometry::nP();
  const int_t              nz     = Geometry::nZProc();
  const int_t              nbound = U -> _nbound;
  const real_t             mu     = Femlib::value ("RHO * KINVIS");
  Vector                   secF;
  vector<real_t>           work(4 * np);
  real_t                   *u, *v, *w;
  int_t           i, j;

  for (j = 0; j < nz; j++) {
    u = U -> _plane[j];
    v = V -> _plane[j];
    w = (W) ? W -> _plane[j] : 0;
    for (i = 0; i < nbound; i++) {
      secF = UBC[i] -> tangTraction ("wall", u, v, &work[0]);
             fx[j] -= mu * secF.x;
             fy[j] -= mu * secF.y;
      if (W) fz[j] -= mu * WBC[i] -> scalarFlux ("wall", w, &work[0]);
    }
  }
}


void Field::traction (real_t*      n, // Normal/pressure
		      real_t*      t, // In-plane tangent/viscous
		      real_t*      s, // Out-of-plane tangent/viscous
		      const int_t  N ,
		      const int_t  M ,
		      const Field* p ,
		      const Field* u ,
		      const Field* v , 
		      const Field* w )
/// --------------------------------------------------------------------------
/// Static member function.
///
/// Compute the pressure and viscous tractions on the "wall" surfaces
/// (the number of which is given as input parameter N). All
/// computations are carried out on Fourier-transformed
/// variables. Input parameter M is the (exchange-padded) length of
/// each variable's wall-tagged storage, per data plane.  //
/// ---------------------------------------------------------------------------
{
  const vector<Boundary*>& UBC    = u -> _bsys -> getBCs(0);
  const int_t              np     = Geometry::nP();
  const int_t              nz     = Geometry::nZProc();
  const int_t              nbound = u -> _nbound;
  const int_t              bmode  = Geometry::baseMode();
  const real_t             mu     = Femlib::value ("RHO * KINVIS");
  const real_t             *ur, *ui, *vr, *vi, *wr, *wi, *pr, *pi;
  real_t                   *nr, *ni, *tr, *ti, *sr, *si;
  int_t                    i, j, k, mode;
  vector<real_t>           work (4 * np);
    
  for (k = 0; k < nz; k += 2) {
    mode = bmode + (k >> 1);

    pr = p -> _plane[k];
    ur = u -> _plane[k];
    vr = v -> _plane[k];
    wr = (w) ? w -> _plane[k] : 0;
    
    pi = (nz > 1) ? p -> _plane[k+1] : 0;
    ui = (nz > 1) ? u -> _plane[k+1] : 0;
    vi = (nz > 1) ? v -> _plane[k+1] : 0;
    wi = (nz > 1) ? w -> _plane[k+1] : 0;

    for (i = 0, j = 0; i < nbound; i++) {

      // -- We loop over all boundaries but only do the work for walls.
      
      if (UBC[i] -> inGroup ("wall")) {

	nr = n + j*np + k*M;
	tr = t + j*np + k*M;
	sr = s + j*np + k*M;

	ni = (nz > 1) ? n + j*np + (k+1)*M : 0;
	ti = (nz > 1) ? t + j*np + (k+1)*M : 0;
	si = (nz > 1) ? s + j*np + (k+1)*M : 0;

	UBC[i] -> traction (mode, mu, pr, pi, ur, ui, vr, vi, wr, wi,
			    nr, ni, tr, ti, sr, si, &work[0]);
	j++;
      }
    }
  }
}


void Field::coupleBCs (Field*      v  ,
		       Field*      w  ,
		       const int_t dir)
/// --------------------------------------------------------------------------
/// Couple/uncouple boundary condition values for the radial and
/// azimuthal velocity fields in cylindrical coordinates, depending on
/// indicated direction.  This action is required due to the coupling
/// in the viscous terms of the;lo N--S equations in cylindrical coords.
///
/// dir == FORWARD
///
/// +         v~ <-- v + i w
/// +         w~ <-- v - i w
///  
/// dir == INVERSE
///
/// +         v  <-- 0.5   * (v~ + w~)
/// +         w  <-- 0.5 i * (w~ - v~)
//
/// Since there is no coupling for the viscous terms in the 2D
/// equation, do nothing for the zeroth Fourier mode.
// ---------------------------------------------------------------------------
{
  if (Geometry::nDim() < 3) return;

  const char     routine[] = "Field::coupleBCs";
  int_t k, Re, Im;
  const int_t    nL    =  v -> _nline;
  const int_t    nMode =  Geometry::nModeProc();
  const int_t    kLo   = (Geometry::procID() == 0) ? 1 : 0;
  vector<real_t> work (nL);
  real_t         *Vr, *Vi, *Wr, *Wi, *tp = &work[0];
  
  if (dir == FORWARD) {

    for (k = kLo; k < nMode; k++) {
      Re = k  + k;
      Im = Re + 1;

      Vr = v -> _line[Re];
      Vi = v -> _line[Im];
      Wr = w -> _line[Re];
      Wi = w -> _line[Im];

      Veclib::copy (nL, Vr, 1, tp, 1);
      Veclib::vsub (nL, Vr, 1, Wi, 1, Vr, 1);
      Veclib::vadd (nL, Wi, 1, tp, 1, Wi, 1);
      Veclib::copy (nL, Wr, 1, tp, 1);
      Veclib::copy (nL, Wi, 1, Wr, 1);
      Veclib::vsub (nL, Vi, 1, tp, 1, Wi, 1);
      Veclib::vadd (nL, Vi, 1, tp, 1, Vi, 1);
    }

  } else if (dir == INVERSE) {

    for (k = kLo; k < nMode; k++) {
      Re = k  + k;
      Im = Re + 1;

      Vr = v -> _line[Re];
      Vi = v -> _line[Im];
      Wr = w -> _line[Re];
      Wi = w -> _line[Im];

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
/// --------------------------------------------------------------------------
/// For cylindrical coordinates with 3D, the radial and azimuthal fields
/// are coupled before solution of the viscous step.  This means that
/// the Fourier constant used for solution may vary from that which
/// applies to the axial component.
///
/// For Field v~, betak -> betak + 1 while for w~, betak -> betak - 1.
///
/// For the uncoupled Fields v, w solved for the zeroth Fourier mode,
/// the "Fourier" constant in the Helmholtz equations is 1.
// ---------------------------------------------------------------------------
{
  if (Geometry::system() == Geometry::Cartesian || 
      name               ==         'c'         ||
      name               ==         'p'         ||
      name               ==         'u'          ) return beta * mode;

  if      (name == 'v') return (mode == 0) ? 1.0 : beta * mode + 1.0;
  else if (name == 'w') return (mode == 0) ? 1.0 : beta * mode - 1.0;
  else Veclib::alert
	 ("Field::modeConstant", "unrecognized Field name given", ERROR);

  return -1.0;			// -- Never happen.
}

   
void Field::overwriteForGroup (const char*     name,
			       const AuxField* src ,
			       AuxField*       tgt )
/// --------------------------------------------------------------------------
/// On a named boundary group (e.g. "axis"), insert field data from src
/// into field data of tgt.  Could use this for inserting values
/// obtained in src via l'Hopital's rule into those in tgt (e.g. on
/// axis!).
// ---------------------------------------------------------------------------
{
  const int_t              nz = Geometry::nZProc();
  const int_t              np = Geometry::nP();
  const vector<Boundary*>& BC = _bsys -> getBCs (0); // -- Mode-invariant.
  real_t                   *a, *b;
  int_t                    i, k, offset, skip;

  for (k = 0; k < nz; k++) {
    a = src -> _plane [k];
    b = tgt -> _plane [k];
    for (i = 0; i < _nbound; i++)
      if (BC[i] -> inGroup (name)) {
	offset = BC[i] -> dOff();
	skip   = BC[i] -> dSkip();
	Veclib::copy (np, a + offset, skip, b + offset, skip);
      }
  }  
}
