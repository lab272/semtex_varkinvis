#ifndef CFEMDEF_H
#define CFEMDEF_H
/*****************************************************************************
 * Common names and typedefs for C or C++ finite-element codes.
 *****************************************************************************/

#if defined(_SX)		/* NEC SX-4.        */
  typedef long int int_t;
  typedef double   real_t;
#else                           /* Everything else. */
  typedef int     int_t;
  typedef double  real_t;
#endif

typedef struct { real_t x, y, z; } Point;
typedef Point                      Vector;

#define UNSET      -1
#define HDR_LEN    351c
#define STR_MAX    2048
#define F77NAME(x) x##_

/* 
   Choose Chebyshev or Legendre polynomials as the underlying basis
   for quadrature points and weights.  Legendre is the default.
*/

#if defined(CHEBYSHEV)
  #define JAC_ALFA -0.5
  #define JAC_BETA -0.5
#else
  #define JAC_ALFA  0.0
  #define JAC_BETA  0.0
#endif

typedef enum quadrature_kind {
  GJ  = 'G',	/* Gauss-Jacobi quadrature.         */
  GLJ = 'L',	/* Gauss-Lobatto-Jacobi quadrature. */
  GRJ = 'R',	/* Gauss-Radau-Jacobi quadrature.   */
  TRZ = 'U'	/* Trapezoidal (uniform mesh).      */
} QuadratureKind;

typedef enum solver_kind {
  DIRECT,	/* Cholesky back-substitution.                           */
  NESTED,       /* Nested dissection.                                    */
  JACPCG,	/* Conjugate gradient, Jacobi (diagonal) preconditioner. */
  MIXED         /* Direct for mode 0, iterative for all others.          */
} SolverKind;

typedef enum transform_kind {
  INVERSE = -1,	/* Basis    --> physical space. */
  FORWARD = +1	/* Physical --> basis    space. */
} TransformKind;

typedef enum exchange_kind {
  FULL,		/* Full exchange (preserves planar and element structure). */
  HALF		/* Half memory exchange across processes (planar only).    */
} ExchangeKind;

#endif
