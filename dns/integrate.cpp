///////////////////////////////////////////////////////////////////////////////
// integrate.cpp: Unsteady Navier--Stokes solver, using
// "stiffly-stable" time integration [1,2].  Geometries may be 2- or
// 3-dimensional, Cartesian or cylindrical [3,5].  Fourier expansions
// are used in the homogeneous (z) direction.  This file provides
// integrate as a call-back routine; after initialisation, integrate
// may be called repeatedly without reinitialising internal storage.
//
// For cylindrical coordinates (Fourier in azimuth):
//   u <==> axial     velocity,  x <==> axial     coordinate direction,
//   v <==> radial    velocity,  y <==> radial    coordinate direction,
//   w <==> azimuthal velocity,  z <==> azimuthal coordinate direction.
//
// For Cartesian coordinates (Fourier in z):
//   u <==> x-component  velocity
//   v <==> y-component  velocity
//   w <==> z-component  velocity
//
// In either system, the w velocity component is optional for 2D
// (N_Z=1) (i.e. can have 2D2C or 2D3C).  If 3D (N_Z > 1), w should
// appear in session.
//
// Optionally integrate concentration of advected scalar field c.
//
// Copyright (c) 1994+, Hugh M Blackburn
//
// REFERENCES
// ----------
// [1] Karniadakis, Israeli & Orszag (1991) "High-order splitting methods
//     for the incompressible Navier--Stokes equations", JCP 97:414--443
// [2] Guermond & Shen (2003) "Velocity correction projection methods for
//     incompressible flows", SIAM J Numer Anal 41:112-134
// [3] Blackburn & Sherwin (2004) "Formulation of a Galerkin spectral
//     element--Fourier method for three-dimensional incompressible flows
//     in cylindrical geometries", JCP 179:759-778
// [4] Dong (2015) "A convective-like energy-stable open boundary condition
//     for simulations of incompressible flows", JCP 302:300-328.
// [5] Blackburn, Lee, Albrecht & Singh (2019) "Semtex: a spectral
//     element--Fourier solver for the incompressible Navier--Stokes
//     equations in cylindrical or Cartesian coordinates", CPC 245:106804.
// [6] Blackburn, Lopez, Singh & Smits (2021) "On the Boussinesq
//     approximation in arbitrarily accelerating frames of reference",
//     JFM 924:R1.
//w   
///////////////////////////////////////////////////////////////////////////////

#include <dns.h>

// -- This triggers computation of nonlinear terms for diagnostics.

#define NONLIN_DIAGNOSTIC 0


typedef ModalMatrixSys Msys;

// -- File-scope constants and routines:

static int_t NDIM, NCOM, NORD, NADV;
static bool  C3D;

static void   waveProp  (Domain*, const AuxField***, const AuxField***);
static void   setPForce (const AuxField**, AuxField**);
static void   project   (const Domain*, AuxField**, AuxField**);
static Msys** preSolve  (const Domain*);
static void   Solve     (Domain*, const int_t, AuxField*, Msys*);


void integrate (void (*advection) (Domain*    , 
                                   BCmgr*     ,
                                   AuxField** , 
                                   AuxField** ,
                                   FieldForce*),
                Domain*      D ,
                BCmgr*       B ,
                DNSAnalyser* A ,
                FieldForce*  FF)
// ---------------------------------------------------------------------------
// On entry, D contains storage (in the following order!) for:
// -- velocity Fields 'u', 'v' (and 'w' if 2D3C or 3D),
// -- optional scalar Field 'c',
// -- constraint Field 'p'.
//
// Us is multi-level auxillary Field storage for velocities and
// Uf is multi-level auxillary Field storage for nonlinear forcing terms.
// ---------------------------------------------------------------------------
{
  NCOM = D -> nVelCmpt();              // -- Number of velocity components.
  NADV = D -> nAdvect();               // -- Number of advected fields.
  NDIM = Geometry::nDim();	       // -- Number of space dimensions.
  NORD = Femlib::ivalue ("N_TIME");    // -- Time integration order.
  C3D  = Geometry::cylindrical() && NDIM == 3;
  
  int_t              i, j, k;
  const real_t       dt    = Femlib:: value ("D_T");
  const int_t        nStep = Femlib::ivalue ("N_STEP");
  const int_t        nZ    = Geometry::nZProc();
  static Msys**      MMS;
  static AuxField*** Us;
  static AuxField*** Uf;
  Field*             Pressure = D -> u[NADV];

  if (!MMS) {			// -- Initialise static storage.

    // -- Create multi-level storage for velocities (Us) and forcing (Uf).

    const int_t ntot  = Geometry::nTotProc();
    real_t*     alloc = new real_t [static_cast<size_t>(2 * NADV*NORD * ntot)];
    Us                = new AuxField** [static_cast<size_t>(2 * NORD)];
    Uf                = Us + NORD;

    for (k = 0, i = 0; i < NORD; i++) {
      Us[i] = new AuxField* [static_cast<size_t>(2 * NADV)];
      Uf[i] = Us[i] + NADV;
      for (j = 0; j < NADV; j++) {
        Us[i][j] = new AuxField (alloc + k++ * ntot, nZ, D -> elmt);
        Uf[i][j] = new AuxField (alloc + k++ * ntot, nZ, D -> elmt);
      }
    }

    // -- Create global matrix systems.

    MMS = preSolve (D);

    // -- Create multi-level storage for pressure BCS.

    B -> buildComputedBCs (Pressure, D -> hasScalar());

    // -- Apply coupling to radial & azimuthal velocity BCs.

    if (C3D) Field::coupleBCs (D -> u[1], D -> u[2], FORWARD);
  }

  // -- Because we may restart from scratch on each call, zero these:

  *Pressure = 0.0;

  for (i = 0; i < NORD; i++)
    for (j = 0; j < NADV; j++) {
      *Us[i][j] = 0.0;
      *Uf[i][j] = 0.0;
    }

#if NONLIN_DIAGNOSTIC
  
  // -- Process input to generate nonlinear terms and quit (if nonzero).
  
  advection (D, B, Us[0], Uf[0], FF);
  D -> step += 1; D -> time += dt;
  for (i = 0; i < NADV; i++) *D -> u[i] = *Uf[0][i];
  A -> analyse (Us[0], Uf[0]);

#else  // -- Normal timestepping.
  
  // -- The following timestepping loop implements equations (15--18) in [5].

  while (D -> step < nStep) {

    // -- Compute nonlinear terms from previous velocity field.
    //    Add physical space forcing, again at old time level.
    //    Outcomes are left in Uf[0], while the velocity fields
    //    for the last time level are left in Us[0].

    advection (D, B, Us[0], Uf[0], FF);

    D -> step += 1;
    D -> time += dt;
    Femlib::ivalue ("STEP", D -> step);
    Femlib::value  ("t",    D -> time);

    // -- Update high-order pressure BC storage, first in BCmgr, then
    //    in pressure Field BC area.

    B -> maintainFourier (D -> step, Pressure,
			  const_cast<const AuxField**>(Us[0]),
			  const_cast<const AuxField**>(Uf[0]),
			  NCOM, NADV);
    Pressure -> evaluateBoundaries (Pressure, D -> step);

    // -- Complete unconstrained advective substep and compute
    //    pressure, which is left in D -> u[NADV].

    if (Geometry::cylindrical()) { Us[0][0] -> mulY(); Us[0][1] -> mulY(); }

    waveProp (D, const_cast<const AuxField***>(Us),
	         const_cast<const AuxField***>(Uf));
    for (i = 0; i < NADV; i++) AuxField::swapData (D -> u[i], Us[0][i]);

    rollm     (Uf, NORD, NADV);
    setPForce (const_cast<const AuxField**>(Us[0]), Uf[0]);
    Solve     (D, NADV,  Uf[0][0], MMS[NADV]);

    // -- Correct velocities for pressure.

    project   (D, Us[0], Uf[0]);

    // -- Update multilevel velocity storage.

    for (i = 0; i < NADV; i++) *Us[0][i] = *D -> u[i];
    rollm (Us, NORD, NADV);

    // -- Re-evaluate velocity (possibly time-dependent) BCs, some of
    //    which get made in physical space and then Fourier
    //    transformed, while others (e.g. some computed types) may get
    //    directly evaluated in Fourier space.  Whatever the method,
    //    the outcomes end up in the (Fourier-transformed) BC storage
    //    areas for the relevant Field.

    for (i = 0; i < NADV; i++)  {
      D -> u[i] -> evaluateBoundaries (NULL,     D -> step, false);
      D -> u[i] -> bTransform         (FORWARD);
      D -> u[i] -> evaluateBoundaries (Pressure, D -> step, true);
    }
    if (C3D) Field::coupleBCs (D -> u[1], D -> u[2], FORWARD);

    // -- Viscous correction substep to complete computation of
    //    velocity components (and, if relevant, scalar) for this time
    //    step.

    if (C3D) {
      AuxField::couple (Uf [0][1], Uf [0][2], FORWARD);
      AuxField::couple (D -> u[1], D -> u[2], FORWARD);
    }

    
    for (i = 0; i < NADV; i++) Solve (D, i, Uf[0][i], MMS[i]);
    if (C3D) AuxField::couple (D -> u[1], D -> u[2], INVERSE);

    // -- Process results of this step.

    A -> analyse (Us[0], Uf[0]);
  }
#endif  
}


static void waveProp (Domain*           D ,
		      const AuxField*** Us,
		      const AuxField*** Uf)
// ---------------------------------------------------------------------------
// Compute the first substep of stiffly-stable timestepping scheme.
//
// On entry, the most recent velocity fields are in Us, and the most
// recent nonlinear terms in Uf.  The intermediate velocity field u^ is
// computed and left in D's velocity areas.
//
// This is the only routine that makes explicit use of the multi time
// level structure of Us & Uf.
// ---------------------------------------------------------------------------
{
  int_t             i, q;
  vector<AuxField*> H (NADV);	// -- Mnemonic for u^{Hat}.

  for (i = 0; i < NADV; i++) {
     H[i] = D -> u[i];
    *H[i] = 0.0;
  }

  const int_t    Je = min (D -> step, NORD);
  vector<real_t> alpha (Integration::OrderMax + 1);
  vector<real_t> beta  (Integration::OrderMax);

  Integration::StifflyStable (Je, &alpha[0]);
  Integration::Extrapolation (Je, &beta [0]);
  Blas::scal (Je, Femlib::value ("D_T"), &beta[0],  1);

  for (i = 0; i < NADV; i++)
    for (q = 0; q < Je; q++) {
      H[i] -> axpy (-alpha[q + 1], *Us[q][i]);
      H[i] -> axpy ( beta [q]    , *Uf[q][i]);
    }
}


static void setPForce (const AuxField** Us,
		       AuxField**       Uf)
// ---------------------------------------------------------------------------
// On input, intermediate velocity storage u^ is in Us.  Create div u^ / D_T
// in the first dimension of Uf as a forcing field for discrete PPE.
// ---------------------------------------------------------------------------
{
  int_t        i;
  const real_t dt = Femlib::value ("D_T");

  for (i = 0; i < NDIM; i++) (*Uf[i] = *Us[i]) . gradient (i);

  for (i = 1; i < NDIM; i++) *Uf[0] += *Uf[i];

  *Uf[0] /= dt;
}


static void project (const Domain* D ,
		     AuxField**    Us,
		     AuxField**    Uf)
// ---------------------------------------------------------------------------
// On input, new pressure field is stored in D and intermediate velocity
// level u^ is stored in Us.  Constrain velocity field:
//
//                    u^^ = u^ - D_T * grad P,
//
// then scale by -1.0 / (D_T * KINVIS) to create forcing for viscous step
// (this is -1.0 / (D_T  * diffusivity) in the case of a scalar field).
//
// u^^ is left in Uf.
// ---------------------------------------------------------------------------
{
  int_t        i;
  const real_t alpha = -1.0 / Femlib::value ("D_T * KINVIS");
  const real_t beta  =  1.0 / Femlib::value ("KINVIS");
  const real_t Pr    =        Femlib::value ("PRANDTL");

  for (i = 0; i < NADV; i++) {
    Field::swapData (Us[i], Uf[i]);
    if (Geometry::cylindrical() && i >= 2) Uf[i] -> mulY();
    *Uf[i] *= alpha;
  }

  // -- For scalar, use diffusivity instead of viscosity.
  if (NADV > NCOM) *Uf[NCOM] *= Pr;

  for (i = 0; i < NDIM; i++) {
    (*Us[0] = *D -> u[NADV]) . gradient (i);
    if (Geometry::cylindrical() && i <  2) Us[0] -> mulY();
    Uf[i] -> axpy (beta, *Us[0]);
  }
}


static Msys** preSolve (const Domain* D)
// ---------------------------------------------------------------------------
// Set up ModalMatrixSystems for each Field of D.  If iterative solution
// is selected for any Field, the corresponding ModalMatrixSystem pointer
// is set to zero.
//
// ITERATIVE >= 1 selects iterative solver for velocity components,
// ITERATIVE >= 2 selects iterative solver for non-zero pressure Fourier modes.
// ---------------------------------------------------------------------------
{
  const int_t             nmodes = Geometry::nModeProc();
  const int_t             base   = Geometry::baseMode();
  const int_t             itLev  = Femlib::ivalue ("ITERATIVE");
  const real_t            beta   = Femlib:: value ("BETA");
  const vector<Element*>& E = D -> elmt;
  Msys**                  M = new Msys* [static_cast<size_t>(NADV + 1)];
  int_t                   i;

  vector<real_t> alpha (Integration::OrderMax + 1);
  Integration::StifflyStable (NORD, &alpha[0]);
  real_t         lambda2 = alpha[0] / Femlib::value ("D_T * KINVIS");

  // -- Velocity systems.

  for (i = 0; i < NCOM; i++)
    M[i] = new Msys
      (lambda2, D -> VARKINVIS,  beta, base, nmodes, E, D -> b[i], D -> n[i],
       (itLev) ? JACPCG : DIRECT);

  // -- Scalar system.

  if (NADV != NCOM) {
    lambda2 = alpha[0] / Femlib::value ("D_T * KINVIS / PRANDTL");
    M[NCOM] = new Msys
      (lambda2, D -> VARKINVIS, beta, base, nmodes, E, D -> b[NCOM], D -> n[NCOM],
       (itLev < 1)?DIRECT:JACPCG);
  }

  // -- Pressure system.

  if (itLev > 1)
    M[NADV] = new Msys
      (0.0, D -> VARKINVIS,  beta, base, nmodes, E, D -> b[NADV], D -> n[NADV], MIXED);
  else
    M[NADV] = new Msys
      (0.0, D -> VARKINVIS,  beta, base, nmodes, E, D -> b[NADV], D -> n[NADV], DIRECT);

  return M;
}


static void Solve (Domain*     D,
		   const int_t i,
		   AuxField*   F,
		   Msys*       M)
// ---------------------------------------------------------------------------
// Solve Helmholtz problem for D -> u[i], using F as a forcing Field.
// Iterative or direct solver selected on basis of field type, step,
// time order and command-line arguments.
// ---------------------------------------------------------------------------
{
  const int_t step = D -> step;

  if (i < NADV && step < NORD) { // -- We need a temporary matrix system.
    const int_t Je     = min (step, NORD);
    const int_t base   = Geometry::baseMode();
    const int_t nmodes = Geometry::nModeProc();

    vector<real_t> alpha (Je + 1);
    Integration::StifflyStable (Je, &alpha[0]);
    const real_t   lambda2 = (i == NCOM) ? // -- True for scalar diffusion.
      alpha[0] / Femlib::value ("D_T * KINVIS / PRANDTL") :
      alpha[0] / Femlib::value ("D_T * KINVIS");
    const real_t   beta    = Femlib::value ("BETA");

    Msys* tmp = new Msys
      (lambda2, D -> VARKINVIS,  beta, base, nmodes, D -> elmt, D -> b[i], D -> n[i], JACPCG);

    D -> u[i] -> solve (F, tmp);
    delete tmp;

  } else D -> u[i] -> solve (F, M);
}

#undef NONLIN_DIAGNOSTIC
