///////////////////////////////////////////////////////////////////////////////
// message.cpp: all semtex routines that use message-passing, which
// are currently MPI-specific.
//
// If MPI_EX isn't defined durinhg compilation, these are just stubs.
//
// See the following for a description of what the exchange routines do:
//
// Rudman M & Blackburn HM (2006) Direct numerical simulation of
//  turbulent non-Newtonian flow using a spectral element method, Appl
//  Math Mod, V30N11: 1229-1248.
//
// If you are looking for routines named "message" or "Veclib::messg",
// these are meant for exception warnings and are part of Veclib.
//
// Copyright (c) 1996+, Hugh M Blackburn
///////////////////////////////////////////////////////////////////////////////

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>

#if defined(MPI_EX)
#include <mpi.h>
/* -- File-scope) Cartesian communicators: */
static MPI_Comm grid_comm, row_comm, col_comm;
#endif

#include <femlib.h>
#include <message.h>

#if defined(NUMA)
/* Need this forward declaration for a SGI/NUMA-specific routine. */
extern void _fastbcopy(const void *src, void *dest, size_t n);
    #define __MEMCPY(dest, src, n) _fastbcopy(src, dest, n)
#else
    #define __MEMCPY(dest, src, n) memcpy(dest, src, n)
#endif


namespace Message {

  
  void init (int*    argc,
	     char*** argv)
  // -----------------------------------------------------------------------
  // Do whatever is required to initialise message-passing.  Set up
  // associated global variables.
  // -----------------------------------------------------------------------
  {
#if defined(MPI_EX)

    int  n;
    int  n_part, dim_sizes[1], wrap_around[1], reorder = 1, free_coords[1];
    char s[STR_MAX];

    MPI_Init       (argc, argv);
    //    Femlib::init   ();

    MPI_Comm_rank  (MPI_COMM_WORLD,   &n);
    sprintf        (s, "I_PROC = %1d", n);
    Femlib::ivalue (s);

    MPI_Comm_size  (MPI_COMM_WORLD,   &n);
    sprintf        (s, "N_PROC = %1d", n);
    Femlib::ivalue (s);

    // -- Allocate a 2D Cartesian MPI grid.  For now, it will only have a
    //    single column (a single 2D partition), and work the same as our
    //    original parallel-across-modes semtex.

    n_part = Femlib::ivalue ("N_PART");           // -- Should be 1.
    dim_sizes[0]   = n;
    wrap_around[0] = 0;
  
    MPI_Cart_create (MPI_COMM_WORLD,
		     1, dim_sizes, wrap_around, reorder,
		     &grid_comm);

    free_coords[0] = 1;

    MPI_Cart_sub  (grid_comm, free_coords, &col_comm);

    MPI_Comm_rank  (col_comm,         &n);
    sprintf        (s, "I_PROC = %1d", n);
    Femlib::ivalue (s);

#endif
  }


  void stop ()
  // ------------------------------------------------------------------------
  // Shut down message passing.
  // ------------------------------------------------------------------------
  {
#if defined(MPI_EX)

    MPI_Barrier  (MPI_COMM_WORLD);
    MPI_Finalize ();

#endif
  }


  void sync ()
  // ------------------------------------------------------------------------
  // Block until all processes have entered.
  // ------------------------------------------------------------------------
  {
#if defined(MPI_EX)

    MPI_Barrier (MPI_COMM_WORLD);

#endif
  }

#define MAX(a, b) ((a) > (b) ? (a) : (b))


  static int first (int n, const int* x)
  { 
    int i;
    for (i = 0; i < n; i++) if (x[i]) return i;
    return 0;
  }


  void exchange (real_t*     data,
		 const int_t nZ  ,
		 const int_t nP  ,
		 const int_t sign)
  // ------------------------------------------------------------------------
  // Transpose blocks of data across processors.  Data is a double-precision 
  // vector, total length nP*nZ on each processor.
  //
  // The amount of data held by each processor is nP*nZ.  Each nP-sized plane
  // can be split into nB = nP/nProc sized blocks (it is assumed that nP is an
  // int_t multiple of nProc, i.e. that nB is a whole number).  Initially
  // the data are ordered by as nZ nP-sized planes, with memory traversed 
  // fastest moving over the planes, and the block indices vary more rapidly
  // than the z-indices.
  //
  // The aim of the the exchange is to re-order the data so that each
  // processor holds all the nZ data for one of the blocks, i.e. a gather of
  // all the z-information for a block onto each a single processor, which e.g.
  // can be followed by multiple 1D Fourier transformations over each block.
  //
  // First the data are exchanged within a single processor so that (in
  // terms of blocks) the block (rather than the z) indices vary slowest
  // as memory is traversed.  This is the "in-place scatter".  Then a
  // block-transpose of data across processors is carried out using
  // message passing.
  //
  // NB: order of inter- and intra-processor exchanges must be reversed 
  // in order to invert the exchange with a second exchange: this is the use
  // of input variable sign.
  // ------------------------------------------------------------------------
  {
#if defined(MPI_EX)

    int            i, j;
    const int      ip = Femlib::ivalue ("I_PROC");
    const int      np = Femlib::ivalue ("N_PROC");
    const int      nB = nP / np;       // -- Size of intra-processor block.
    const int      NB = nP / nB;       // -- Number of these blocks in a plane.
    const int      NM = nP * nZ / np;  // -- Size of message block.
    const int      dsize   = sizeof (real_t);
    static double  *tmp   = NULL;
    static int     *kmove = NULL, *kpost, lastk;
    static int     *jmove = NULL, *jpost, lastj;
    static int     lastreq = 0;

    if (np == 1) return;

    if (tmp && lastreq != nP * nZ) { free (tmp); tmp = NULL; }
    if (!tmp) { lastreq = nP * nZ; tmp = (double*) malloc (lastreq * dsize); }

    if (sign == 1) {		// -- "Forwards" exchange.

      // -- Intra-processor exchange.

      if (NB == nZ) {		// -- Symmetric exchange.
	for (i = 0; i < nZ; i++)
	  for (j = i; j < nZ; j++) {
	    if (i != j) {
	      __MEMCPY (tmp,                data + (i*nZ+j)*nB, nB * dsize);
	      __MEMCPY (data + (i*nZ+j)*nB, data + (j*nZ+i)*nB, nB * dsize);
	      __MEMCPY (data + (j*nZ+i)*nB, tmp,                nB * dsize);
	    }
	  }

      } else {			// -- Asymmetric exchange.

	int        k, knext, kconf;
	const int  NBnZm = NB * nZ - 1;

	if (kmove && lastk != nZ*NB) { free (kmove); kmove = NULL; }
	if (!kmove) {
	  lastk = nZ * NB;
	  kmove = (int*) malloc (2*lastk * sizeof (int));
	  kpost = kmove + lastk;
	}

	// -- Build scatter indices.
    
	for (i = 0; i < NB; i++)
	  for (j = 0; j < nZ; j++)
	    kpost[j * NB + i] = i * nZ + j;

	for (i = 1; i < NBnZm; i++) kmove[i] = 1;
	kmove[0] = kmove[NBnZm] = 0;

	// -- Do "in-place" scatter.

	while (k = first (nZ*NB, kmove)) {
	  knext = kconf = kpost[k];
	  __MEMCPY (tmp, data + kconf * nB, nB * dsize);
	  do {
	    __MEMCPY (data + knext * nB, data + k * nB, nB * dsize);
	    kmove[k] = 0;
	    for (i = 1; i < NBnZm; i++)
	      if (kpost[i] == k) {
		k     = i;
		knext = kpost[k];
		break;
	      }
	  } while (k != kconf);
	  __MEMCPY (data + knext * nB, tmp, nB * dsize);
	  kmove[k] = 0;
	}
      }

      // -- Inter-processor transpose, with NB blocks size nZ*nB / processor.

      MPI_Alltoall (data, NM, MPI_DOUBLE, tmp, NM, MPI_DOUBLE, col_comm);
      __MEMCPY     (data, tmp, nP * nZ * dsize);


    } else {			// -- "Backwards" exchange.

      MPI_Alltoall (data, NM, MPI_DOUBLE, tmp, NM, MPI_DOUBLE, col_comm);
      __MEMCPY     (data, tmp, nP * nZ * dsize);

      if (NB == nZ) {

	for (i = 0; i < nZ; i++)
	  for (j = i; j < nZ; j++) {
	    if (i != j) {
	      __MEMCPY (tmp,                data + (i*nZ+j)*nB, nB * dsize);
	      __MEMCPY (data + (i*nZ+j)*nB, data + (j*nZ+i)*nB, nB * dsize);
	      __MEMCPY (data + (j*nZ+i)*nB, tmp,                nB * dsize);
	    }
	  }

      } else {

	int        j, jnext, jconf;
	const int  NBnZm = NB * nZ - 1;

	if (jmove && lastj != nZ*NB) { free (jmove); jmove = NULL; }
	if (!jmove) {
	  lastj = nZ * NB;
	  jmove = (int*) malloc (2*lastj * sizeof (int));
	  jpost = jmove + lastj;
	}

	for (i = 0; i < NB; i++)
	  for (j = 0; j < nZ; j++)
	    jpost[i * nZ + j] = j * NB + i;

	for (i = 1; i < NBnZm; i++) jmove[i] = 1;
	jmove[0] = jmove[NBnZm] = 0;

	while (j = first (nZ*NB, jmove)) {
	  jnext = jconf = jpost[j];
	  __MEMCPY (tmp, data + jconf * nB, nB * dsize);
	  do {
	    __MEMCPY (data + jnext * nB, data + j * nB, nB * dsize);
	    jmove[j] = 0;
	    for (i = 1; i < NBnZm; i++)
	      if (jpost[i] == j) {
		j     = i;
		jnext = jpost[j];
		break;
	      }
	  } while (j != jconf);
	  __MEMCPY (data + jnext * nB, tmp, nB * dsize);
	  jmove[j] = 0;
	}
      }
    }
#endif
  }


  void exchange (int_t*      data,
		 const int_t nZ  ,
		 const int_t nP  ,
		 const int_t sign)
  // ------------------------------------------------------------------------
  // Integer exchange.
  // ------------------------------------------------------------------------
  {
#if defined(MPI_EX)

    int            i, j;
    const int      ip = Femlib::ivalue ("I_PROC");
    const int      np = Femlib::ivalue ("N_PROC");
    const int      nB = nP / np;       // -- Size of intra-processor block.
    const int      NB = nP / nB;       // -- Number of these blocks in a plane.
    const int      NM = nP * nZ / np;  // -- Size of message block.
    const int      dsize   = sizeof (int_t);
    static int     *tmp   = NULL;
    static int     *kmove = NULL, *kpost, lastk;
    static int     *jmove = NULL, *jpost, lastj;
    static int     lastreq = 0;

    if (np == 1) return;

    if (tmp && lastreq != nP * nZ) { free (tmp); tmp = NULL; }
    if (!tmp) { lastreq = nP * nZ; tmp = (int*) malloc (lastreq * dsize); }

    if (sign == 1) {		// -- "Forwards" exchange.

      // -- Intra-processor exchange.

      if (NB == nZ) {		// -- Symmetric exchange.
	for (i = 0; i < nZ; i++)
	  for (j = i; j < nZ; j++) {
	    if (i != j) {
	      __MEMCPY (tmp,                data + (i*nZ+j)*nB, nB * dsize);
	      __MEMCPY (data + (i*nZ+j)*nB, data + (j*nZ+i)*nB, nB * dsize);
	      __MEMCPY (data + (j*nZ+i)*nB, tmp,                nB * dsize);
	    }
	  }

      } else {			// -- Asymmetric exchange.

	int        k, knext, kconf;
	const int  NBnZm = NB * nZ - 1;

	if (kmove && lastk != nZ*NB) { free (kmove); kmove = NULL; }
	if (!kmove) {
	  lastk = nZ * NB;
	  kmove = (int*) malloc (2*lastk * sizeof (int));
	  kpost = kmove + lastk;
	}

	// -- Build scatter indices.
    
	for (i = 0; i < NB; i++)
	  for (j = 0; j < nZ; j++)
	    kpost[j * NB + i] = i * nZ + j;

	for (i = 1; i < NBnZm; i++) kmove[i] = 1;
	kmove[0] = kmove[NBnZm] = 0;

	// -- Do "in-place" scatter.

	while (k = first (nZ*NB, kmove)) {
	  knext = kconf = kpost[k];
	  __MEMCPY (tmp, data + kconf * nB, nB * dsize);
	  do {
	    __MEMCPY (data + knext * nB, data + k * nB, nB * dsize);
	    kmove[k] = 0;
	    for (i = 1; i < NBnZm; i++)
	      if (kpost[i] == k) {
		k     = i;
		knext = kpost[k];
		break;
	      }
	  } while (k != kconf);
	  __MEMCPY (data + knext * nB, tmp, nB * dsize);
	  kmove[k] = 0;
	}
      }

      // -- Inter-processor transpose, with NB blocks size nZ*nB / processor.

      MPI_Alltoall (data, NM, MPI_INT, tmp, NM, MPI_INT, col_comm);
      __MEMCPY     (data, tmp, nP * nZ * dsize);

    } else {			// -- "Backwards" exchange.

      MPI_Alltoall (data, NM, MPI_INT, tmp, NM, MPI_INT, col_comm);
      __MEMCPY     (data, tmp, nP * nZ * dsize);

      if (NB == nZ) {

	for (i = 0; i < nZ; i++)
	  for (j = i; j < nZ; j++) {
	    if (i != j) {
	      __MEMCPY (tmp,                data + (i*nZ+j)*nB, nB * dsize);
	      __MEMCPY (data + (i*nZ+j)*nB, data + (j*nZ+i)*nB, nB * dsize);
	      __MEMCPY (data + (j*nZ+i)*nB, tmp,                nB * dsize);
	    }
	  }

      } else {

	int        j, jnext, jconf;
	const int  NBnZm = NB * nZ - 1;

	if (jmove && lastj != nZ*NB) { free (jmove); jmove = NULL; }
	if (!jmove) {
	  lastj = nZ * NB;
	  jmove = (int*) malloc (2*lastj * sizeof (int));
	  jpost = jmove + lastj;
	}

	for (i = 0; i < NB; i++)
	  for (j = 0; j < nZ; j++)
	    jpost[i * nZ + j] = j * NB + i;

	for (i = 1; i < NBnZm; i++) jmove[i] = 1;
	jmove[0] = jmove[NBnZm] = 0;

	while (j = first (nZ*NB, jmove)) {
	  jnext = jconf = jpost[j];
	  __MEMCPY (tmp, data + jconf * nB, nB * dsize);
	  do {
	    __MEMCPY (data + jnext * nB, data + j * nB, nB * dsize);
	    jmove[j] = 0;
	    for (i = 1; i < NBnZm; i++)
	      if (jpost[i] == j) {
		j     = i;
		jnext = jpost[j];
		break;
	      }
	  } while (j != jconf);
	  __MEMCPY (data + jnext * nB, tmp, nB * dsize);
	  jmove[j] = 0;
	}
      }
    }
#endif
  }


  void send (real_t*     data,
	     const int_t N   ,
	     const int_t tgt )
  // ------------------------------------------------------------------------
  // Send data to processor number tgt.
  // ------------------------------------------------------------------------
  {
#if defined(MPI_EX)

    MPI_Send (data, (int) N, MPI_DOUBLE, (int) tgt, 0, col_comm);

#endif
  }


  void recv (real_t*     data,
	     const int_t N   ,
	     const int_t src )
  // ------------------------------------------------------------------------
  // Receive data from processor number src.
  // ------------------------------------------------------------------------
  {
#if defined(MPI_EX)

    MPI_Status status;

    MPI_Recv (data, (int) N, MPI_DOUBLE, (int) src, 0, col_comm, &status);

#endif
  }


  void send (int_t*      data,
	     const int_t N   ,
	     const int_t tgt )
  // ------------------------------------------------------------------------
  // Send data to processor number tgt.
  // ------------------------------------------------------------------------
  {
#if defined(MPI_EX)

    if (sizeof (int_t) == sizeof (int))
      MPI_Send (data, (int) N, MPI_INT,  (int) tgt, 0, col_comm);
    else
      MPI_Send (data, (int) N, MPI_LONG, (int) tgt, 0, col_comm);

#endif
  }


  void recv (int_t*      data,
	     const int_t N   ,
	     const int_t src )
  // ------------------------------------------------------------------------
  // Receive data from processor number src.
  // ------------------------------------------------------------------------
  {
#if defined(MPI_EX)

    MPI_Status status;

    if (sizeof (int_t) == sizeof (int))
      MPI_Recv (data, (int) N, MPI_INT,  (int) src, 0, col_comm, &status);
    else
      MPI_Recv (data, (int) N, MPI_LONG, (int) src, 0, col_comm, &status);

#endif
  }

}
