///////////////////////////////////////////////////////////////////////////////
// nonlinear.cpp
//
// Compute nonlinear (forcing) terms in Navier--Stokes equations: N(u) + f.
//
// Here N(u) represents the nonlinear advection terms in the N--S
// equations transposed to the RHS i.e. -div.grad(u), and f is a
// vector of body force per unit mass (with possible space-time
// dependency).
//
// Copyright (c) 1994+, Hugh M Blackburn
///////////////////////////////////////////////////////////////////////////////

#include <dns.h>

// -- Turn on/off hydrostatic correction for grad(KE) terms. Usually off.

#define _KE_HYDROSTAT 0


void skewSymmetric (Domain*     D ,
		    BCmgr*      B ,
		    AuxField**  Us,
		    AuxField**  Uf,
		    FieldForce* FF)
// ---------------------------------------------------------------------------
// Velocity field data areas of D (which on entry contain velocity
// data from the previous timestep, in Fourier space) and first level
// of Us are swapped (so that subsequently Us stores the old velocity
// data), then the next stage of nonlinear forcing terms N(u) are
// computed from velocity fields and left in the first (i.e. the
// supplied) level of Uf.
//
// Nonlinear terms N(u) in skew-symmetric form are
//                 ~ ~
//           N  = -0.5 ( u . grad u + div uu )
//           ~           ~        ~       ~~
//
// i.e., in Cartesian component form
//
//           N  = -0.5 ( u  d(u ) / dx  + d(u u ) / dx ).
//            i           j    i      j      i j      j
//
// in cylindrical coordinates
//
//           Nx = -0.5 {ud(u)/dx + vd(u)/dy +  d(uu)/dx + d(vu)/dy +
//                 1/y [wd(u)/dz + d(uw)/dz + vu      ]}
//           Ny = -0.5 {ud(v)/dx + vd(v)/dy +  d(uv)/dx + d(vv)/dy +
//                 1/y [wd(v)/dz + d(vw)/dz + vv - 2ww]}
//           Nz = -0.5 {ud(w)/dx + vd(w)/dy +  d(uw)/dx + d(vw)/dy +
//                 1/y [wd(w)/dz + d(ww)/dz + 3wv     ]}
//
// NB: for the cylindrical coordinate formulation we actually here 
// compute y*Nx, y*Ny, Nz, as outlined in Blackburn & Sherwin (2004).
//
// If a scalar, c, is present, also compute the advection term
// -0.5*[u.grad(c)+div(uc)].  This is a straightforward extension to
// the above for Cartesian coordinates (the loop over i is extended by
// 1 and the last component of 'u_i' is c), whereas in cylindrical
// coordinates we have
//
//           Nc = -0.5 {ud(c)/dx + vd(c)/dy + d(uc)/dx + d(vc)/dy +
//                 1/y [wd(c)/dz + d(wc)/dz + cv ]}
//
// Data are transformed to physical space for most of the operations.
// For gradients in the Fourier direction however, the data must be
// transferred back to Fourier space prior to differention in z, then
// trsnsformed back.  The final form of N delivered at the end of the
// routine is in Fourier space.  Note that there is no (longer any)
// provision for dealiasing in the Fourier direction and that all
// storage areas of D->u (including pressure) are overwritten here.
//  
// NB: for the cylindrical coordinate formulation we actually here 
// compute y*Nx, y*Ny, Nz, as outlined in Blackburn & Sherwin (2004).  
// ---------------------------------------------------------------------------
{
  const int_t NDIM = Geometry::nDim();
  const int_t NADV = D -> nAdvect();
  const int_t NCOM = D -> nVelCmpt();

  vector<AuxField*> U (NADV), N (NADV), Uphys (NADV);
  AuxField*         tmp = D -> u[NADV]; // -- Pressure is used for scratch.
  int_t             i, j;

  for (i = 0; i < NADV; i++) {
    Uphys[i] = D -> u[i];
    N[i]     = Uf[i];
    U[i]     = Us[i];
    *N[i]    = 0.0;
    *U[i]    = *Uphys[i];
    Uphys[i] -> transform (INVERSE);
  }

  B -> maintainPhysical (D -> u[0], Uphys, NCOM, NADV);

  if (Geometry::cylindrical()) {

    for (i = 0; i < NADV; i++) {

      if (i == 0) N[0] -> timesMinus (*Uphys[0], *Uphys[1]);
      if (i == 1) N[1] -> timesMinus (*Uphys[1], *Uphys[1]);

      // -- Terms involving azimuthal derivatives and frame components.

      if (NCOM == 3) {
        if (i == 1) {
          tmp -> times (*Uphys[2], *Uphys[2]);
          N[1] -> axpy ( 2.0, *tmp);
        }

        if (i == 2) {
          tmp -> times (*Uphys[2], *Uphys[1]);
          N[2] -> axpy (-3.0, *tmp);
        }

	if (i == 3)
	  N[3] -> timesMinus (*Uphys[3], *Uphys[1]);
	  
        if (NDIM == 3) {
          (*tmp = *U[i]) . gradient (2) . transform (INVERSE);
          N[i] -> timesMinus (*Uphys[2], *tmp);

          tmp -> times (*Uphys[i], *Uphys[2]);
          (*tmp) . transform (FORWARD). gradient (2). transform (INVERSE);
          *N[i] -= *tmp;
        }
      } else if (i == 2 && (NADV > NCOM))
	N[2] -> timesMinus (*Uphys[1], *Uphys[2]);

      if (i >= 2) N[i] -> divY ();
      
      // -- 2D convective derivatives.

      for (j = 0; j < 2; j++) {
        (*tmp = *Uphys[i]) . gradient (j);
        if (i < 2) tmp -> mulY ();
        N[i] -> timesMinus (*Uphys[j], *tmp);
      }

      // -- 2D conservative derivatives.

      for (j = 0; j < 2; j++) {
        (*tmp). times (*Uphys[j], *Uphys[i]) . gradient (j);
        if (i < 2) tmp -> mulY ();
        *N[i] -= *tmp;
      }

      *N[i] *= 0.5;      // -- Average the two forms to get skew-symmetric.

      if (i < NCOM) FF -> addPhysical (N[i], tmp, i, Uphys);
    }

  } else {			// -- Cartesian coordinates.

    for (i = 0; i < NADV; i++) {
      for (j = 0; j < NDIM; j++) {

        // -- Perform n_i -= u_j d(u_i) / dx_j.

        if (j == 2) (*tmp = *U[i]) . gradient (j) . transform (INVERSE);
        else    (*tmp = *Uphys[i]) . gradient (j);
        N[i] -> timesMinus (*Uphys[j], *tmp);

        // -- Perform n_i -= d(u_i u_j) / dx_j.

        tmp -> times (*Uphys[i], *Uphys[j]);
        if (j == 2) tmp -> transform (FORWARD);
        tmp -> gradient (j);
        if (j == 2) tmp -> transform (INVERSE);
        *N[i] -= *tmp;
      }

      *N[i] *= 0.5;      // -- Average the two forms to get skew-symmetric.

      if (i < NCOM) FF -> addPhysical (N[i], tmp, i, Uphys);
    }
  }

#if 1
  // -- Multiply in density variation (1 + rho'/rho_0) for CSB buoyancy.

  if (D -> hasScalar()) FF -> canonicalSteadyBoussinesq (tmp, Uphys, N);

#else
  
  // -- Multiply in density variation (1 + rho'/rho_0) for LMA13 buoyancy.

  if (D -> hasScalar() && (Femlib::value ("LMA_BETA_T") > EPSDP)) {
    *tmp  = Femlib::value ("LMA_T_REF");
    *tmp -= *Uphys[NCOM];
    *tmp *= Femlib::value ("LMA_BETA_T");
    *tmp += 1.0;
    for (i = 0; i < NCOM; i++) *N[i] *= *tmp;

    // -- Subtract out background hydrostatic contributions.

    // -- 1. Frame-acceleration (incl. gravity) terms.
    
    for (i = 0; i < NCOM; i++) FF -> subPhysical (N[i], tmp, i, Uphys);

#if _KE_HYDROSTAT
    // -- 2. Localised centrifugal buoyancy terms.
    
    tmp -> innerProduct (Uphys, Uphys, NCOM) *= 0.5;

    if (Geometry::cylindrical())
      for (i = 0; i < NCOM; i++)
        if (i == 2) {
	  if (NDIM == 3)
	    *N[i] -= (*Uphys[0] = *tmp) .
	      transform (FORWARD). gradient (i) . transform (INVERSE) . divY();
        } else
	  *N[i] -= (*Uphys[0] = *tmp) . gradient (i) . mulY();
    else
      for (i = 0; i < NCOM; i++)
        if (i == 2) {
	  if (NDIM == 3)
	    *N[i] -= (*Uphys[0] = *tmp) .
	      transform (FORWARD). gradient (i) . transform (INVERSE);
        } else
	  *N[i] -= (*Uphys[0] = *tmp) . gradient (i);
#endif
  }
#endif
  
  for (i = 0; i < NADV; i++) {
    N[i] -> transform (FORWARD);
    N[i] -> smooth (D -> nGlobal(), D -> assemblyNaive(), D -> invMassNaive());
  }
}


void altSkewSymmetric (Domain*     D ,
		       BCmgr*      B ,
		       AuxField**  Us,
		       AuxField**  Uf,
		       FieldForce* FF)
// ---------------------------------------------------------------------------
// This is a cheaper approximation to the full skew-symmetric
// computation (see next routine).  The formulation of nonlinear terms
// used here is so-called "alternating skew symmetric" method (first
// documented by Bob Kerr) which uses the non-conservative and
// conservative forms of the nonlinear terms on alternating
// timesteps. This has shown in testing to be almost as robust as the
// full skew symmetric method but costs half as much.
//
// NB: this is now used as the default advection scheme.  As noted, it
// is about as robust but computationally cheaper than full
// skew-symmetric.  However, it has a slight defect, too: owing to
// alternation, it may produce small high-frequency temporal
// oscillations in the solution, of period 2*D_T.  If that is a
// problem, try full skew-symmetric (previous routine) or the
// non-conservative form of the nonlinear terms (next routine).
//
// Here are the two components of skew symmetric form in cylindrical
// coords (cf. what is written above for the full skew symmetric form):
//
// Non-conservative form:
//
//           Nx = -{ud(u)/dx + vd(u)/dy + 1/y [wd(u)/dz     ]}
//           Ny = -{ud(v)/dx + vd(v)/dy + 1/y [wd(v)/dz - ww]}
//           Nz = -{ud(w)/dx + vd(w)/dy + 1/y [wd(w)/dz + wv]}
//           Nc = -{ud(c)/dx + vd(c)/dy + 1/y [wd(c)/dz     ]}
//
// Conservative form:
//
//           Nx = -{d(uu)/dx + d(vu)/dy + 1/y [d(uw)/dz + vu     ]}
//           Ny = -{d(uv)/dx + d(vv)/dy + 1/y [d(vw)/dz + vv - ww]}
//           Nz = -{d(uw)/dx + d(vw)/dy + 1/y [d(ww)/dz + 2wv    ]}
//           Nc = -{d(uc)/dx + d(vc)/dy + 1/y [d(wc)/dz + cv     ]}
//
// NB: for the cylindrical coordinate formulation we actually here 
// compute y*Nx, y*Ny, Nz, as outlined in Blackburn & Sherwin (2004).  
// ---------------------------------------------------------------------------
{
  const int_t NDIM = Geometry::nDim();
  const int_t NADV = D -> nAdvect();
  const int_t NCOM = D -> nVelCmpt();
  
  vector<AuxField*> U (NADV), N (NADV), Uphys (NADV);
  AuxField*         tmp = D -> u[NADV]; // -- Pressure is used for scratch.
  int_t             i, j;
  static int        toggle = 1;            // -- Switch u.grad(u) or div(uu).
  
  for (i = 0; i < NADV; i++) {
    Uphys[i] = D -> u[i];
    N[i]     = Uf[i];
    U[i]     = Us[i];
    *N[i]    = 0.0;
    *U[i]    = *Uphys[i];
    Uphys[i] -> transform (INVERSE);
  }

  B -> maintainPhysical (D -> u[0], Uphys, NCOM, NADV);

  if (Geometry::cylindrical()) {

    if (toggle) { // -- Convective component u.grad(u).

      for (i = 0; i < NADV; i++) {

        // -- Terms involving azimuthal derivatives and frame components.

        if (NCOM == 3) {
          if (i == 1) N[1] -> times      (*Uphys[2], *Uphys[2]);
          if (i == 2) N[2] -> timesMinus (*Uphys[2], *Uphys[1]);

          if (NDIM == 3) {
            (*tmp = *U[i]) . gradient (2) . transform (INVERSE);
            N[i] -> timesMinus (*Uphys[2], *tmp);
          }
	  if (i >= 2) N[i] -> divY ();	  
        }

        // -- 2D convective derivatives.

        for (j = 0; j < 2; j++) {
          (*tmp = *Uphys[i]) . gradient (j);
          if (i <  2) tmp -> mulY ();
          N[i] -> timesMinus (*Uphys[j], *tmp);
        }

        if (i < NCOM) FF -> addPhysical (N[i], tmp, i, Uphys);
      }

    } else { // -- Conservative component div(uu).

      for (i = 0; i < NADV; i++) {

        // -- Terms involving azimuthal derivatives and frame components.

        if (i == 0) N[0] -> timesMinus (*Uphys[0], *Uphys[1]);
        if (i == 1) N[1] -> timesMinus (*Uphys[1], *Uphys[1]);

        if (NCOM == 3) {
          if (i == 1) N[1] -> timesPlus (*Uphys[2], *Uphys[2]);
          if (i == 2) {
            tmp -> times (*Uphys[2], *Uphys[1]);
            N[2] -> axpy (-2., *tmp);
          }
	  if (i == 3) N[3] -> timesMinus (*Uphys[3], *Uphys[1]);
          if (NDIM == 3) {
            tmp -> times (*Uphys[i], *Uphys[2]);
            (*tmp) . transform (FORWARD). gradient (2). transform (INVERSE);
            *N[i] -= *tmp;
          }
	} else if (i == 2 && (NADV > NCOM))
	  N[2] -> timesMinus (*Uphys[1], *Uphys[2]);

	if (i >= 2) N[i] -> divY ();
	
        // -- 2D conservative derivatives.

        for (j = 0; j < 2; j++) {
          (*tmp). times (*Uphys[j], *Uphys[i]) . gradient (j);
          if (i <  2) tmp -> mulY ();
          *N[i] -= *tmp;
        }

        if (i < NCOM) FF -> addPhysical (N[i], tmp, i, Uphys);

      }
    }

  } else {			// -- Cartesian coordinates.

    if (toggle) { // -- Convective component u.grad(u).

      for (i = 0; i < NADV; i++) {
        for (j = 0; j < NDIM; j++) { // -- Perform n_i -= u_j d(u_i) / dx_j.

          if (j == 2) (*tmp = *U[i]) . gradient (j) . transform (INVERSE);
          else    (*tmp = *Uphys[i]) . gradient (j);
          N[i] -> timesMinus (*Uphys[j], *tmp);
        }

        if (i < NCOM) FF -> addPhysical (N[i], tmp, i, Uphys);

      }

    } else { // -- Conservative component div(uu).

      for (i = 0; i < NADV; i++) {
        for (j = 0; j < NDIM; j++) { // -- Perform n_i -= d(u_i u_j) / dx_j.

          tmp -> times (*Uphys[i], *Uphys[j]);
          if (j == 2) tmp -> transform (FORWARD);
          tmp -> gradient (j);
          if (j == 2) tmp -> transform (INVERSE);
          *N[i] -= *tmp;
        }

        if (i < NCOM) FF -> addPhysical (N[i], tmp, i, Uphys);

      }
    }
  }

#if 1
  // -- Multiply in density variation (1 + rho'/rho_0) for CSB buoyancy.

  if (D -> hasScalar()) FF -> canonicalSteadyBoussinesq (tmp, Uphys, N);

#else
  
  // -- Multiply in density variation (1 + rho'/rho_0) for LMA13 buoyancy.

  if (D -> hasScalar() && (Femlib::value ("LMA_BETA_T") > EPSDP)) {
    *tmp  = Femlib::value ("LMA_T_REF");
    *tmp -= *Uphys[NCOM];
    *tmp *= Femlib::value ("LMA_BETA_T");
    *tmp += 1.0;
    for (i = 0; i < NCOM; i++) *N[i] *= *tmp;

    // -- Subtract out background hydrostatic contributions.

    // -- 1. Frame-acceleration (incl. gravity) terms.
    
    for (i = 0; i < NCOM; i++) FF -> subPhysical (N[i], tmp, i, Uphys);

#if _KE_HYDROSTAT
    // -- 2. Localised centrifugal buoyancy terms.
    
    tmp -> innerProduct (Uphys, Uphys, NCOM) *= 0.5;

    if (Geometry::cylindrical())
      for (i = 0; i < NCOM; i++)
        if (i == 2) {
	  if (NDIM == 3)
	    *N[i] -= (*Uphys[0] = *tmp) .
	      transform (FORWARD). gradient (i) . transform (INVERSE) . divY();
        } else
	  *N[i] -= (*Uphys[0] = *tmp) . gradient (i) . mulY();
    else
      for (i = 0; i < NCOM; i++)
        if (i == 2) {
	  if (NDIM == 3)
	    *N[i] -= (*Uphys[0] = *tmp) .
	      transform (FORWARD). gradient (i) . transform (INVERSE);
        } else
	  *N[i] -= (*Uphys[0] = *tmp) . gradient (i);
#endif
  }
#endif
 
  for (i = 0; i < NADV; i++) {
    N[i] -> transform (FORWARD);
    N[i] -> smooth (D -> nGlobal(), D -> assemblyNaive(), D -> invMassNaive());
  }
  
  toggle = 1 - toggle;
}


void convective (Domain*     D ,
		 BCmgr*      B ,
		 AuxField**  Us,
		 AuxField**  Uf,
		 FieldForce* FF)
// ---------------------------------------------------------------------------
// Nonlinear terms N(u) in standard (non-conservative) convective form are
//                 ~ ~
//           N  = - u . grad u
//           ~      ~        ~
//
// This has a similar operation count to altSkewSymmetric() but seems
// somewhat less stable/robust at high Reynolds numbers.  However, it
// *is* smooth in time, so if stability isn't an issue, but temporal
// smoothness and computational cost are, give it a try if you don't
// want to pay for full skew-symmetric.
//
// i.e., in Cartesian component form
//
//           N  = -  u  d(u ) / dx 
//            i       j    i      j
//
// in cylindrical coordinates
//
//           Nx = -{ud(u)/dx + vd(u)/dy + 1/y [wd(u)/dz]}
//           Ny = -{ud(v)/dx + vd(v)/dy + 1/y [wd(v)/dz - ww]}
//           Nz = -{ud(w)/dx + vd(w)/dy + 1/y [wd(w)/dz + wv]}
//           Nc = -{ud(c)/dx + vd(c)/dy + 1/y [wd(c)/dz     ]}  
//
// NB: for the cylindrical coordinate formulation we actually here 
// compute y*Nx, y*Ny, Nz, as outlined in Blackburn & Sherwin (2004).
// ---------------------------------------------------------------------------
{
  const int_t NDIM = Geometry::nDim();
  const int_t NADV = D -> nAdvect();
  const int_t NCOM = D -> nVelCmpt();

  vector<AuxField*> U (NADV), N (NADV), Uphys (NADV);
  AuxField*         tmp = D -> u[NADV]; // -- Pressure is used for scratch.
  int_t             i, j;

  for (i = 0; i < NADV; i++) {
    Uphys[i] = D -> u[i];
    N[i]     = Uf[i];
    U[i]     = Us[i];
    *N[i]    = 0.0;
    *U[i]    = *Uphys[i];
    Uphys[i] -> transform (INVERSE);
  }

  B -> maintainPhysical (D -> u[0], Uphys, NCOM, NADV);

  if (Geometry::cylindrical()) {

    for (i = 0; i < NADV; i++) {

      // -- Terms involving azimuthal derivatives and frame components.

      if (NCOM == 3) {
        if (i == 1) N[1] -> times      (*Uphys[2], *Uphys[2]);
        if (i == 2) N[2] -> timesMinus (*Uphys[2], *Uphys[1]);

        if (NDIM == 3) {
          (*tmp = *U[i]) . gradient (2) . transform (INVERSE);
          N[i] -> timesMinus (*Uphys[2], *tmp);
        }
	if (i >= 2) N[i] -> divY ();	
      }

      // -- 2D convective derivatives.

      for (j = 0; j < 2; j++) {
        (*tmp = *Uphys[i]) . gradient (j);
        if (i < 2) tmp -> mulY ();
        N[i] -> timesMinus (*Uphys[j], *tmp);
      }

      if (i < NCOM) FF -> addPhysical (N[i], tmp, i, Uphys);

    }

  } else {			// -- Cartesian coordinates.

    for (i = 0; i < NADV; i++) {
      for (j = 0; j < NDIM; j++) { // -- Perform n_i -= u_j d(u_i) / dx_j.

        if (j == 2) (*tmp = *U[i]) . gradient (j) . transform (INVERSE);
        else    (*tmp = *Uphys[i]) . gradient (j);
        N[i] -> timesMinus (*Uphys[j], *tmp);
      }

      if (i < NCOM) FF -> addPhysical (N[i], tmp, i, Uphys);
    }
  }

#if 1
  // -- Multiply in density variation (1 + rho'/rho_0) for CSB buoyancy.

  if (D -> hasScalar()) FF -> canonicalSteadyBoussinesq (tmp, Uphys, N);

#else
  
  // -- Multiply in density variation (1 + rho'/rho_0) for LMA13 buoyancy.

  if (D -> hasScalar() && (Femlib::value ("LMA_BETA_T") > EPSDP)) {
    *tmp  = Femlib::value ("LMA_T_REF");
    *tmp -= *Uphys[NCOM];
    *tmp *= Femlib::value ("LMA_BETA_T");
    *tmp += 1.0;
    for (i = 0; i < NCOM; i++) *N[i] *= *tmp;

    // -- Subtract out background hydrostatic contributions.

    // -- 1. Frame-acceleration (incl. gravity) terms.
    
    for (i = 0; i < NCOM; i++) FF -> subPhysical (N[i], tmp, i, Uphys);

#if _KE_HYDROSTAT
    // -- 2. Localised centrifugal buoyancy terms.
    
    tmp -> innerProduct (Uphys, Uphys, NCOM) *= 0.5;

    if (Geometry::cylindrical())
      for (i = 0; i < NCOM; i++)
        if (i == 2) {
	  if (NDIM == 3)
	    *N[i] -= (*Uphys[0] = *tmp) .
	      transform (FORWARD). gradient (i) . transform (INVERSE) . divY();
        } else
	  *N[i] -= (*Uphys[0] = *tmp) . gradient (i) . mulY();
    else
      for (i = 0; i < NCOM; i++)
        if (i == 2) {
	  if (NDIM == 3)
	    *N[i] -= (*Uphys[0] = *tmp) .
	      transform (FORWARD). gradient (i) . transform (INVERSE);
        } else
	  *N[i] -= (*Uphys[0] = *tmp) . gradient (i);
#endif
  }
#endif

  for (i = 0; i < NADV; i++) {
    N[i] -> transform (FORWARD);
    N[i] -> smooth (D -> nGlobal(), D -> assemblyNaive(), D -> invMassNaive());
  }
}


void rotational1 (Domain*     D ,
		  BCmgr*      B ,
		  AuxField**  Us,
		  AuxField**  Uf,
		  FieldForce* FF)
// ---------------------------------------------------------------------------
// Nonlinear terms N(u) in rotational form are
//                 ~ ~
//           N  = - omega x  u = u x omega = u x curl(u)
//           ~      ~        ~   ~   ~       ~        ~
// where omega is vorticity (curl u) and x is cross-product.  We will deal
// with scalar (c) as for convective form.  
// 
// In Cartesian component form
//
//           Nx = v[d(v)/dx - d(u)/dy] + w[d(w)/dx - d(u)/dz] 
//           Ny = w[d(w)/dy - d(v)/dz] + u[d(u)/dy - d(v)/dx]
//           Nz = u[d(u)/dz - d(w)/dx] + v[d(v)/dz - d(w)/dy]
//           Nc = -ud(c)/dx - vd(c)/dy - wd(c)/dz  
//
// in cylindrical coordinates (premultiplied by y on Nx, Ny, terms, see BS04)
//
//           Nx = y{v[d(v)/dx - d(u)/dy] + wd(w)/dx} - wd(u)/dx
//           Ny = y{u[d(u)/dy - d(v)/dx] + wd(w)/dy} - wd(v)/dz + ww
//           Nz = -ud(w)/dx - vd(w)/dy  + 1/y {vd(v)/dz + ud(u)/dz - vw}
//           Nc = -ud(c)/dx - vd(c)/dy  - 1/y {wd(c)/dz}
//
// ---------------------------------------------------------------------------
{
  const char routine[] = "rotational1";
  
  const int_t NDIM = Geometry::nDim();
  const int_t NADV = D -> nAdvect();
  const int_t NCOM = D -> nVelCmpt();
  const bool  D3C3 = (NDIM == 3) && (NCOM == 3);
  const bool  D2C3 = (NDIM == 2) && (NCOM == 3);
  const bool  D2C2 = (NDIM == 2) && (NCOM == 2);
  const bool  scalar = NADV > NCOM;

  vector<AuxField*> U (NADV), N (NADV), Uphys (NADV);
  AuxField*         tmp = D -> u[NADV]; // -- Pressure is used for scratch.
  int_t             i, j;

  for (i = 0; i < NADV; i++) {
    Uphys[i] = D -> u[i];
    N[i]     = Uf[i];
    U[i]     = Us[i];
    *N[i]    = 0.0;
    *U[i]    = *Uphys[i];
    Uphys[i] -> transform (INVERSE);
  }

  B -> maintainPhysical (D -> u[0], Uphys, NCOM, NADV);

  if (Geometry::cylindrical()) {

    // -- First, the momentum equations.

    if (D3C3) {	          // -- 3D3C: two z-derivatives have to be made twice.

      N[2] -> timesMinus (*Uphys[1], *Uphys[2]);
      (*tmp = *U[0]) . gradient (2) . transform (INVERSE);
      N[2] -> timesPlus (*Uphys[0], *tmp);
      (*tmp = *U[1]) . gradient (2) . transform (INVERSE);
      N[2] -> timesPlus (*Uphys[1], *tmp);
      N[2] -> divY();
 
      (*tmp = *Uphys[0]) . gradient (1);
      N[0] -> timesMinus (*Uphys[1], *tmp);
      N[1] -> timesPlus  (*Uphys[0], *tmp);

      (*tmp = *Uphys[1]) . gradient (0);
      N[0] -> timesPlus  (*Uphys[1], *tmp);
      N[1] -> timesMinus (*Uphys[0], *tmp);

      (*tmp = *Uphys[2]) . gradient (0);
      N[0] -> timesPlus  (*Uphys[2], *tmp) . mulY();
      N[2] -> timesMinus (*Uphys[0], *tmp);

      (*tmp = *Uphys[2]) . gradient (1);
      N[1] -> timesPlus  (*Uphys[2], *tmp) . mulY();
      N[2] -> timesMinus (*Uphys[1], *tmp);

      (*tmp = *U[0]) . gradient (2) . transform (INVERSE);
      N[0] -> timesMinus (*Uphys[2], *tmp);
      (*tmp = *U[1]) . gradient (2) . transform (INVERSE);      
      N[1] -> timesMinus (*Uphys[2], *tmp);

      N[1] -> timesPlus (*Uphys[2], *Uphys[2]);
      
    } else if (D2C3) { // -- 2D3C.

      N[2] -> timesMinus (*Uphys[1], *Uphys[2]);
      N[2] -> divY();
 
      (*tmp = *Uphys[0]) . gradient (1);
      N[0] -> timesMinus (*Uphys[1], *tmp);
      N[1] -> timesPlus  (*Uphys[0], *tmp);

      (*tmp = *Uphys[1]) . gradient (0);
      N[0] -> timesPlus  (*Uphys[1], *tmp);
      N[1] -> timesMinus (*Uphys[0], *tmp);

      (*tmp = *Uphys[2]) . gradient (0);
      N[0] -> timesPlus  (*Uphys[2], *tmp) . mulY();
      N[2] -> timesMinus (*Uphys[0], *tmp);

      (*tmp = *Uphys[2]) . gradient (1);
      N[1] -> timesPlus  (*Uphys[2], *tmp) . mulY();
      N[2] -> timesMinus (*Uphys[1], *tmp);

      N[1] -> timesPlus (*Uphys[2], *Uphys[2]);
      
    } else if (D2C2) { // -- 2D2C.

      (*tmp = *Uphys[0]) . gradient (1);
      N[0] -> timesMinus (*Uphys[1], *tmp);
      N[1] -> timesPlus  (*Uphys[0], *tmp);

      (*tmp = *Uphys[1]) . gradient (0);
      N[0] -> timesPlus  (*Uphys[1], *tmp);
      N[1] -> timesMinus (*Uphys[0], *tmp);

      N[0] -> mulY();
      N[1] -> mulY();

      N[1] -> timesPlus (*Uphys[2], *Uphys[2]);

    } else
      Veclib::alert (routine, "Never get here", ERROR);

    for (i = 0; i < NCOM; i++) FF -> addPhysical (N[i], tmp, i, Uphys);

    // -- Address scalar advection separately.
    
    if (scalar) {
      if (D3C3) {  // -- 3D3C.
	(*tmp = *U[NCOM]) . gradient (2) . transform (INVERSE);
	N[NCOM] -> timesMinus (*Uphys[2], *tmp) . divY();
      }           // -- 2D3C and 2D2C.     
      (*tmp = *Uphys[NCOM]) . gradient (0);
      N[NCOM] -> timesMinus (*Uphys[0], *tmp);
      (*tmp = *Uphys[NCOM]) . gradient (1);
      N[NCOM] -> timesMinus (*Uphys[1], *tmp);
    }	

  } else {			// -- Cartesian coordinates.

    if (D3C3) {

      (*tmp = *Uphys[0]) . gradient (1);
      N[0] -> timesMinus (*Uphys[1], *tmp);
      N[1] -> timesPlus  (*Uphys[0], *tmp);

      (*tmp = *U[0]) . gradient (2) . transform (INVERSE);
      N[0] -> timesMinus (*Uphys[2], *tmp);
      N[2] -> timesPlus  (*Uphys[0], *tmp);

      (*tmp = *Uphys[1]) . gradient (0);
      N[0] -> timesPlus  (*Uphys[1], *tmp);
      N[1] -> timesMinus (*Uphys[0], *tmp);

      (*tmp = *U[1]) . gradient (2) . transform (INVERSE);
      N[1] -> timesMinus (*Uphys[2], *tmp);
      N[2] -> timesPlus  (*Uphys[1], *tmp);

      (*tmp = *Uphys[2]) . gradient (0);
      N[0] -> timesPlus  (*Uphys[2], *tmp);
      N[2] -> timesMinus (*Uphys[0], *tmp);

      (*tmp = *Uphys[2]) . gradient (1);
      N[1] -> timesPlus  (*Uphys[2], *tmp);
      N[2] -> timesMinus (*Uphys[1], *tmp);
      
    } else if (D2C3) {

      (*tmp = *Uphys[0]) . gradient (1);
      N[0] -> timesMinus (*Uphys[1], *tmp);
      N[1] -> timesPlus  (*Uphys[0], *tmp);

      (*tmp = *Uphys[1]) . gradient (0);
      N[0] -> timesPlus  (*Uphys[1], *tmp);
      N[1] -> timesMinus (*Uphys[0], *tmp);

      (*tmp = *Uphys[2]) . gradient (0);
      N[0] -> timesPlus  (*Uphys[2], *tmp);
      N[2] -> timesMinus (*Uphys[0], *tmp);

      (*tmp = *Uphys[2]) . gradient (1);
      N[1] -> timesPlus  (*Uphys[2], *tmp);
      N[2] -> timesMinus (*Uphys[1], *tmp);
      
    } else if (D2C2) {

      (*tmp = *Uphys[0]) . gradient (1);
      N[0] -> timesMinus (*Uphys[1], *tmp);
      N[1] -> timesPlus  (*Uphys[0], *tmp);

      (*tmp = *Uphys[1]) . gradient (0);
      N[0] -> timesPlus  (*Uphys[1], *tmp);
      N[1] -> timesMinus (*Uphys[0], *tmp);

    } else
      Veclib::alert (routine, "Never get here", ERROR);

    for (i = 0; i < NCOM; i++) FF -> addPhysical (N[i], tmp, i, Uphys);

    if (scalar) {
      (*tmp = *Uphys[NCOM]) . gradient (0);
      N[NCOM] -> timesMinus (*Uphys[0], *tmp);
      (*tmp = *Uphys[NCOM]) . gradient (1);
      N[NCOM] -> timesMinus (*Uphys[1], *tmp);
      if (D3C3) {
	(*tmp = *U[NCOM]) . gradient (2) . transform (INVERSE);
	N[NCOM] -> timesMinus (*Uphys[2], *tmp);
      }
    }
  }

#if 1
  // -- Multiply in density variation (1 + rho'/rho_0) for CSB buoyancy.

  if (D -> hasScalar()) FF -> canonicalSteadyBoussinesq (tmp, Uphys, N);

#else
  
  // -- Multiply in density variation (1 + rho'/rho_0) for LMA13 buoyancy.

  if (D -> hasScalar() && (Femlib::value ("LMA_BETA_T") > EPSDP)) {
    *tmp  = Femlib::value ("LMA_T_REF");
    *tmp -= *Uphys[NCOM];
    *tmp *= Femlib::value ("LMA_BETA_T");
    *tmp += 1.0;
    for (i = 0; i < NCOM; i++) *N[i] *= *tmp;

    // -- Subtract out background hydrostatic contributions.

    // -- 1. Frame-acceleration (incl. gravity) terms.
    
    for (i = 0; i < NCOM; i++) FF -> subPhysical (N[i], tmp, i, Uphys);

#if _KE_HYDROSTAT
    // -- 2. Localised centrifugal buoyancy terms.
    
    tmp -> innerProduct (Uphys, Uphys, NCOM) *= 0.5;

    if (Geometry::cylindrical())
      for (i = 0; i < NCOM; i++)
        if (i == 2) {
	  if (NDIM == 3)
	    *N[i] -= (*Uphys[0] = *tmp) .
	      transform (FORWARD). gradient (i) . transform (INVERSE) . divY();
        } else
	  *N[i] -= (*Uphys[0] = *tmp) . gradient (i) . mulY();
    else
      for (i = 0; i < NCOM; i++)
        if (i == 2) {
	  if (NDIM == 3)
	    *N[i] -= (*Uphys[0] = *tmp) .
	      transform (FORWARD). gradient (i) . transform (INVERSE);
        } else
	  *N[i] -= (*Uphys[0] = *tmp) . gradient (i);
#endif
  }
#endif
  

#if 0
    // -- In Rot-1 form we have to explicitly add back terms
    //    associated with non-frame centrifugal buoyancy.

    // -- Leaving this code here as a reminder it has been cut out...
    
    tmp -> innerProduct (Uphys, Uphys, NCOM) *= 0.5;

    *Uphys[1]  = Femlib::value ("LMA_T_REF");
    *Uphys[1] -= *Uphys[NCOM];
    *Uphys[1] *= Femlib::value ("LMA_BETA_T"); // -- rho'/rho_0.

    if (Geometry::cylindrical())
      for (i = 0; i < NCOM; i++)
        if (i == 2) {
	  if (NDIM == 3) {
	    (*Uphys[0] = *tmp) .
	      transform (FORWARD) . gradient (i) . transform (INVERSE) . divY();
	    N[i] -> timesMinus (*Uphys[1], *Uphys[0]);
	  }
        } else {
	  (*Uphys[0] = *tmp) . gradient (i) . mulY();
	  N[i] -> timesMinus (*Uphys[1], *Uphys[0]);
	}
    else
      for (i = 0; i < NCOM; i++)
        if (i == 2) {
	  if (NDIM == 3) {
	    (*Uphys[0] = *tmp) .
	      transform (FORWARD). gradient (i) . transform (INVERSE);
	    N[i] -> timesMinus (*Uphys[1], *Uphys[0]);
	  }
        } else {
	  (*Uphys[0] = *tmp) . gradient (i);
	  N[i] -> timesMinus (*Uphys[1], *Uphys[0]);
	}
    // -- 2. Localised centrifugal buoyancy terms (none present in Rot-1 form).
    
#endif

  for (i = 0; i < NADV; i++) {
    N[i] -> transform (FORWARD);
    N[i] -> smooth (D -> nGlobal(), D -> assemblyNaive(), D -> invMassNaive());
  }  
}


void rotational2 (Domain*     D ,
		  BCmgr*      B ,
		  AuxField**  Us,
		  AuxField**  Uf,
		  FieldForce* FF)
// ---------------------------------------------------------------------------
// Nonlinear terms N(u) in rotational2 form are
//                 ~ ~
//           N  = u x omega - grad(q) = u x curl(u) - grad (q)
//           ~    ~   ~                 ~        ~
// where omega is vorticity x is cross-product and q = KE = 0.5 u.u.
// We will deal with scalar (c) as for convective form.         ~ ~
//
// This is a fairly straightforward modification of rotational1 to include
// terms related to gradient of KE.  
// 
// In Cartesian component form
//
//           Nx = v[d(v)/dx - d(u)/dy] + w[d(w)/dx - d(u)/dz] - d(q)/dx
//           Ny = w[d(w)/dy - d(v)/dz] + u[d(u)/dy - d(v)/dx] - d(q)/dy
//           Nz = u[d(u)/dz - d(w)/dx] + v[d(v)/dz - d(w)/dy] - d(q)/dz
//           Nc = -ud(c)/dx - vd(c)/dy - wd(c)/dz  
//
// in cylindrical coordinates (premultiplied by y on Nx, Ny, terms, see BS04)
//
//           Nx = y{v[d(v)/dx - d(u)/dy] + wd(w)/dx - d(q)/dx} - wd(u)/dx
//           Ny = y{u[d(u)/dy - d(v)/dx] + wd(w)/dy - d(q)/dy} - wd(v)/dz + ww
//           Nz = -ud(w)/dx - vd(w)/dy + 1/y {vd(v)/dz + ud(u)/dz -vw -d(q)/dz}
//           Nc = -ud(c)/dx - vd(c)/dy  - 1/y {wd(c)/dz}
//
// ---------------------------------------------------------------------------
{
  const char routine[] = "rotational2";
  
  const int_t NDIM = Geometry::nDim();
  const int_t NADV = D -> nAdvect();
  const int_t NCOM = D -> nVelCmpt();
  const bool  D3C3 = (NDIM == 3) && (NCOM == 3);
  const bool  D2C3 = (NDIM == 2) && (NCOM == 3);
  const bool  D2C2 = (NDIM == 2) && (NCOM == 2);
  const bool  scalar = NADV > NCOM;

  vector<AuxField*> U (NADV), N (NADV), Uphys (NADV);
  AuxField*         tmp = D -> u[NADV]; // -- Pressure is used for scratch.
  int_t             i, j;

  for (i = 0; i < NADV; i++) {
    Uphys[i] = D -> u[i];
    N[i]     = Uf[i];
    U[i]     = Us[i];
    *N[i]    = 0.0;
    *U[i]    = *Uphys[i];
    Uphys[i] -> transform (INVERSE);
  }

  B -> maintainPhysical (D -> u[0], Uphys, NCOM, NADV);

  if (Geometry::cylindrical()) {

    if (D3C3) {		      // -- Two z-derivatives have to be made twice.

      tmp -> innerProduct (Uphys, Uphys, 3) *= -0.5;

      (*N[0] = *tmp) . gradient (0);
      (*N[1] = *tmp) . gradient (1);
      (*N[2] = *tmp) . transform (FORWARD) . gradient (2) . transform (INVERSE);

      N[2] -> timesMinus (*Uphys[1], *Uphys[2]);
      (*tmp = *U[0]) . gradient (2) . transform (INVERSE);
      N[2] -> timesPlus (*Uphys[0], *tmp);
      (*tmp = *U[1]) . gradient (2) . transform (INVERSE);
      N[2] -> timesPlus (*Uphys[1], *tmp);
      N[2] -> divY();
 
      (*tmp = *Uphys[0]) . gradient (1);
      N[0] -> timesMinus (*Uphys[1], *tmp);
      N[1] -> timesPlus  (*Uphys[0], *tmp);

      (*tmp = *Uphys[1]) . gradient (0);
      N[0] -> timesPlus  (*Uphys[1], *tmp);
      N[1] -> timesMinus (*Uphys[0], *tmp);

      (*tmp = *Uphys[2]) . gradient (0);
      N[0] -> timesPlus  (*Uphys[2], *tmp) . mulY();
      N[2] -> timesMinus (*Uphys[0], *tmp);

      (*tmp = *Uphys[2]) . gradient (1);
      N[1] -> timesPlus  (*Uphys[2], *tmp) . mulY();
      N[2] -> timesMinus (*Uphys[1], *tmp);

      (*tmp = *U[0]) . gradient (2) . transform (INVERSE);
      N[0] -> timesMinus (*Uphys[2], *tmp);
      (*tmp = *U[1]) . gradient (2) . transform (INVERSE);      
      N[1] -> timesMinus (*Uphys[2], *tmp);

      N[1] -> timesPlus (*Uphys[2], *Uphys[2]);
      
    } else if (D2C3) {

      tmp -> innerProduct (Uphys, Uphys, 3) *= -0.5;

      (*N[0] = *tmp) . gradient (0);
      (*N[1] = *tmp) . gradient (1);

      N[2] -> timesMinus (*Uphys[1], *Uphys[2]);
      N[2] -> divY();
 
      (*tmp = *Uphys[0]) . gradient (1);
      N[0] -> timesMinus (*Uphys[1], *tmp);
      N[1] -> timesPlus  (*Uphys[0], *tmp);

      (*tmp = *Uphys[1]) . gradient (0);
      N[0] -> timesPlus  (*Uphys[1], *tmp);
      N[1] -> timesMinus (*Uphys[0], *tmp);

      (*tmp = *Uphys[2]) . gradient (0);
      N[0] -> timesPlus  (*Uphys[2], *tmp) . mulY();
      N[2] -> timesMinus (*Uphys[0], *tmp);

      (*tmp = *Uphys[2]) . gradient (1);
      N[1] -> timesPlus  (*Uphys[2], *tmp) . mulY();
      N[2] -> timesMinus (*Uphys[1], *tmp);

      N[1] -> timesPlus (*Uphys[2], *Uphys[2]);
      
    } else if (D2C2) {

      tmp -> innerProduct (Uphys, Uphys, 2) *= -0.5;

      (*N[0] = *tmp) . gradient (0);
      (*N[1] = *tmp) . gradient (1);

      (*tmp = *Uphys[0]) . gradient (1);
      N[0] -> timesMinus (*Uphys[1], *tmp);
      N[1] -> timesPlus  (*Uphys[0], *tmp);

      (*tmp = *Uphys[1]) . gradient (0);
      N[0] -> timesPlus  (*Uphys[1], *tmp);
      N[1] -> timesMinus (*Uphys[0], *tmp);

      N[0] -> mulY();
      N[1] -> mulY();

      N[1] -> timesPlus (*Uphys[2], *Uphys[2]);

    } else
      Veclib::alert (routine, "Never get here", ERROR);

    for (i = 0; i < NCOM; i++) FF -> addPhysical (N[i], tmp, i, Uphys);

    if (scalar) {
      if (D3C3) {
	(*tmp = *U[NCOM]) . gradient (2) . transform (INVERSE);
	N[NCOM] -> timesMinus (*Uphys[2], *tmp) . divY();
      }      
      (*tmp = *Uphys[NCOM]) . gradient (0);
      N[NCOM] -> timesMinus (*Uphys[0], *tmp);
      (*tmp = *Uphys[NCOM]) . gradient (1);
      N[NCOM] -> timesMinus (*Uphys[1], *tmp);
    }	

  } else {			// -- Cartesian coordinates.

    if (D3C3) {

      tmp -> innerProduct (Uphys, Uphys, 3) *= -0.5;

      (*N[0] = *tmp) . gradient (0);
      (*N[1] = *tmp) . gradient (1);
      (*N[2] = *tmp) . transform (FORWARD) . gradient (2) . transform (INVERSE);

      (*tmp = *Uphys[0]) . gradient (1);
      N[0] -> timesMinus (*Uphys[1], *tmp);
      N[1] -> timesPlus  (*Uphys[0], *tmp);

      (*tmp = *U[0]) . gradient (2) . transform (INVERSE);
      N[0] -> timesMinus (*Uphys[2], *tmp);
      N[2] -> timesPlus  (*Uphys[0], *tmp);

      (*tmp = *Uphys[1]) . gradient (0);
      N[0] -> timesPlus  (*Uphys[1], *tmp);
      N[1] -> timesMinus (*Uphys[0], *tmp);

      (*tmp = *U[1]) . gradient (2) . transform (INVERSE);
      N[1] -> timesMinus (*Uphys[2], *tmp);
      N[2] -> timesPlus  (*Uphys[1], *tmp);

      (*tmp = *Uphys[2]) . gradient (0);
      N[0] -> timesPlus  (*Uphys[2], *tmp);
      N[2] -> timesMinus (*Uphys[0], *tmp);

      (*tmp = *Uphys[2]) . gradient (1);
      N[1] -> timesPlus  (*Uphys[2], *tmp);
      N[2] -> timesMinus (*Uphys[1], *tmp);
      
    } else if (D2C3) {

      tmp -> innerProduct (Uphys, Uphys, 3) *= -0.5;

      (*N[0] = *tmp) . gradient (0);
      (*N[1] = *tmp) . gradient (1);

      (*tmp = *Uphys[0]) . gradient (1);
      N[0] -> timesMinus (*Uphys[1], *tmp);
      N[1] -> timesPlus  (*Uphys[0], *tmp);

      (*tmp = *Uphys[1]) . gradient (0);
      N[0] -> timesPlus  (*Uphys[1], *tmp);
      N[1] -> timesMinus (*Uphys[0], *tmp);

      (*tmp = *Uphys[2]) . gradient (0);
      N[0] -> timesPlus  (*Uphys[2], *tmp);
      N[2] -> timesMinus (*Uphys[0], *tmp);

      (*tmp = *Uphys[2]) . gradient (1);
      N[1] -> timesPlus  (*Uphys[2], *tmp);
      N[2] -> timesMinus (*Uphys[1], *tmp);
      
    } else if (D2C2) {

      tmp -> innerProduct (Uphys, Uphys, 2) *= -0.5;

      (*N[0] = *tmp) . gradient (0);
      (*N[1] = *tmp) . gradient (1);

      (*tmp = *Uphys[0]) . gradient (1);
      N[0] -> timesMinus (*Uphys[1], *tmp);
      N[1] -> timesPlus  (*Uphys[0], *tmp);

      (*tmp = *Uphys[1]) . gradient (0);
      N[0] -> timesPlus  (*Uphys[1], *tmp);
      N[1] -> timesMinus (*Uphys[0], *tmp);

    } else
      Veclib::alert (routine, "Never get here", ERROR);

    for (i = 0; i < NCOM; i++) FF -> addPhysical (N[i], tmp, i, Uphys);

    if (scalar) {
      (*tmp = *Uphys[NCOM]) . gradient (0);
      N[NCOM] -> timesMinus (*Uphys[0], *tmp);
      (*tmp = *Uphys[NCOM]) . gradient (1);
      N[NCOM] -> timesMinus (*Uphys[1], *tmp);
      if (D3C3) {
	(*tmp = *U[NCOM]) . gradient (2) . transform (INVERSE);
	N[NCOM] -> timesMinus (*Uphys[2], *tmp);
      }
    }
  }

#if 1
  // -- Multiply in density variation (1 + rho'/rho_0) for CSB buoyancy.

  if (D -> hasScalar()) FF -> canonicalSteadyBoussinesq (tmp, Uphys, N);

#else
  
  // -- Multiply in density variation (1 + rho'/rho_0) for LMA13 buoyancy.

  if (D -> hasScalar() && (Femlib::value ("LMA_BETA_T") > EPSDP)) {
    *tmp  = Femlib::value ("LMA_T_REF");
    *tmp -= *Uphys[NCOM];
    *tmp *= Femlib::value ("LMA_BETA_T");
    *tmp += 1.0;
    for (i = 0; i < NCOM; i++) *N[i] *= *tmp;

    // -- Subtract out background hydrostatic contributions.

    // -- 1. Frame-acceleration (incl. gravity) terms.
    
    for (i = 0; i < NCOM; i++) FF -> subPhysical (N[i], tmp, i, Uphys);

#if _KE_HYDROSTAT
    // -- 2. Localised centrifugal buoyancy terms.
    
    tmp -> innerProduct (Uphys, Uphys, NCOM) *= 0.5;

    if (Geometry::cylindrical())
      for (i = 0; i < NCOM; i++)
        if (i == 2) {
	  if (NDIM == 3)
	    *N[i] -= (*Uphys[0] = *tmp) .
	      transform (FORWARD). gradient (i) . transform (INVERSE) . divY();
        } else
	  *N[i] -= (*Uphys[0] = *tmp) . gradient (i) . mulY();
    else
      for (i = 0; i < NCOM; i++)
        if (i == 2) {
	  if (NDIM == 3)
	    *N[i] -= (*Uphys[0] = *tmp) .
	      transform (FORWARD). gradient (i) . transform (INVERSE);
        } else
	  *N[i] -= (*Uphys[0] = *tmp) . gradient (i);
#endif
  }
#endif

  for (i = 0; i < NADV; i++) {
    N[i] -> transform (FORWARD);
    N[i] -> smooth (D -> nGlobal(), D -> assemblyNaive(), D -> invMassNaive());
  }  
}


void Stokes (Domain*     D ,
	     BCmgr*      B ,
	     AuxField**  Us,
	     AuxField**  Uf,
	     FieldForce* FF)
// ---------------------------------------------------------------------------
// Stokes flow by definition does not have nonlinear terms but we
// still need to zero storage areas and add in body force ff.
// ---------------------------------------------------------------------------
{
  const int_t NDIM = Geometry::nDim();
  const int_t NADV = D -> nAdvect();
  const int_t NCOM = D -> nVelCmpt();

  vector<AuxField*> U (NADV), N (NADV), Uphys (NADV);
  AuxField*         tmp    = D -> u[NADV]; // -- Pressure is used for scratch.
  int_t             i, j;

  for (i = 0; i < NADV; i++) {
    Uphys[i] = D -> u[i];
    N[i]     = Uf[i];
    U[i]     = Us[i];
    *N[i]    = 0.0;
    *U[i]    = *Uphys[i];
    Uphys[i] -> transform (INVERSE);
  }
  
  B -> maintainPhysical (D -> u[0], Uphys, NCOM, NADV);
 
  for (i = 0; i < NCOM; i++) FF -> addPhysical (N[i], tmp, i, Uphys);

#if 1
  // -- Multiply in density variation (1 + rho'/rho_0) for CSB buoyancy.

  if (D -> hasScalar()) FF -> canonicalSteadyBoussinesq (tmp, Uphys, N);

#else
  
  // -- Multiply in density variation (1 + rho'/rho_0) for LMA13 buoyancy.

  if (D -> hasScalar() && (Femlib::value ("LMA_BETA_T") > EPSDP)) {
    *tmp  = Femlib::value ("LMA_T_REF");
    *tmp -= *Uphys[NCOM];
    *tmp *= Femlib::value ("LMA_BETA_T");
    *tmp += 1.0;
    for (i = 0; i < NCOM; i++) *N[i] *= *tmp;

    // -- Subtract out background hydrostatic contributions.

    // -- 1. Frame-acceleration (incl. gravity) terms.
    
    for (i = 0; i < NCOM; i++) FF -> subPhysical (N[i], tmp, i, Uphys);

#if _KE_HYDROSTAT
    // -- 2. Localised centrifugal buoyancy terms.
    
    tmp -> innerProduct (Uphys, Uphys, NCOM) *= 0.5;

    if (Geometry::cylindrical())
      for (i = 0; i < NCOM; i++)
        if (i == 2) {
	  if (NDIM == 3)
	    *N[i] -= (*Uphys[0] = *tmp) .
	      transform (FORWARD). gradient (i) . transform (INVERSE) . divY();
        } else
	  *N[i] -= (*Uphys[0] = *tmp) . gradient (i) . mulY();
    else
      for (i = 0; i < NCOM; i++)
        if (i == 2) {
	  if (NDIM == 3)
	    *N[i] -= (*Uphys[0] = *tmp) .
	      transform (FORWARD). gradient (i) . transform (INVERSE);
        } else
	  *N[i] -= (*Uphys[0] = *tmp) . gradient (i);
#endif
  }
#endif

  for (i = 0; i < NADV; i++) {
    N[i] -> transform (FORWARD);
    N[i] -> smooth (D -> nGlobal(), D -> assemblyNaive(), D -> invMassNaive());
  }  
}

#undef _KE_HYDROSTAT
