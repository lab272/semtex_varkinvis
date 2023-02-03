//////////////////////////////////////////////////////////////////////////////
// helmholtz.cpp: routines to solve elliptic problems in one variable.
//
// Copyright (c) 1994+, Hugh M Blackburn
//
//////////////////////////////////////////////////////////////////////////////

#include <sem.h>


void Helmholtz (Domain*   D,
		AuxField* F)
// ---------------------------------------------------------------------------
// Solve Helmholtz's equation
//                                  2
//               div grad u - LAMBDA  u = f(x,y,z)
//
// subject to BCs.
//
// See
//
// Blackburn & Sherwin (2004) "Formulation of a Galerkin spectral
// element--Fourier method for three-dimensional incompressible flows
// in cylindrical geometries", JCP 179:759-778
//
// for a discussion of methodology. In particular, we need to
// pre-multiply the forcing by radius if solving in cylindrical
// coordinataes since that factor does not form part of the quadrature
// weights for the right-hand-side of elliptic equations which are
// solved in symmetrised form (pre-multuplied by radius).  The factor
// is omitted in order to correctly cope with right-hand-sides which
// are constructed as a divergence.
// ---------------------------------------------------------------------------
{
  const real_t lambda2 = Femlib::value ("LAMBDA2");
  const real_t beta    = Femlib::value ("BETA");
  const int_t  nmodes  = Geometry::nModeProc();
  const int_t  base    = Geometry::baseMode();
  const int_t  nz      = Geometry::nZProc();
  SolverKind   method  = (Femlib::ivalue("ITERATIVE")) ? JACPCG : DIRECT;

  ModalMatrixSys* M = new ModalMatrixSys
    (lambda2, beta, base, nmodes, D -> elmt, D -> b[0], D -> n[0], method);

  if (Geometry::cylindrical()) F -> mulY();

  D -> u[0] -> solve (F, M);

  D -> step++;
}
