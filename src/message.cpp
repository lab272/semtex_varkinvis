///////////////////////////////////////////////////////////////////////////////
// message.cpp: all semtex routines that use message-passing, which
// are currently MPI-specific.
//
// If MPI_EX isn't defined during compilation, these are just stubs.
//
// See the following for a description of what the exchange routines do:
//
// Rudman M & Blackburn HM (2006) Direct numerical simulation of
//  turbulent non-Newtonian flow using a spectral element method, Appl
//  Math Mod, V30N11: 1229-1248.
//
// If you are looking for routines named "message()" or
// "Veclib::alert()", these are meant for exception warnings and are
// part of semtex/veclib.
//
// Copyright (c) 1996+, Hugh M Blackburn
///////////////////////////////////////////////////////////////////////////////

#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <ctime>

#include <utility.h>
#include <veclib.h>
#include <message.h>

#if defined(MPI_EX)
#include <mpi.h>
/* -- File-scope Cartesian communicators: */
static MPI_Comm
  grid_comm = MPI_COMM_NULL,
  row_comm  = MPI_COMM_NULL,
  col_comm  = MPI_COMM_NULL;
#endif

#if defined(NUMA)
// -- Need this forward declaration for a SGI/NUMA-specific routine.
extern void _fastbcopy(const void *src, void *dest, size_t n);
    #define __MEMCPY(dest, src, n) _fastbcopy(src, dest, n)
#else
    #define __MEMCPY(dest, src, n) memcpy(dest, src, n)
#endif


namespace Message {

  void init (int*    argc ,
	     char*** argv ,
	     int&    nproc,
	     int&    iproc)
  // ------------------------------------------------------------------------
  // Initialise message-passing, which strips off MPI-associated
  // runtime command-line arguments (and so should happen before we
  // try to deal with our own).  Return the overall number of
  // processes and the identifier of the current process.
  // ------------------------------------------------------------------------
  {
#if defined(MPI_EX)
    
    MPI_Init      (argc, argv);
    MPI_Comm_size (MPI_COMM_WORLD, &nproc);
    MPI_Comm_rank (MPI_COMM_WORLD, &iproc);

#endif
  }    


  void grid (const int& npart2d, // -- Input [1 .. number of elements in mesh]
	     int&       ipart2d, //    Output
	     int&       npartz , //    Output
	     int&       ipartz ) //    Output
  // ------------------------------------------------------------------------
  // Set up a Cartesian grid for domain decomposition and return the
  // coordinates of the current process within that grid.  The idea is
  // that for our 2D x Fourier discretisation, the decomposition can
  // be set up as a 2D Cartesian grid: the first dimension (rows)
  // spans the 2D partitions, while the second dimension (columns)
  // spans modal blocks in the Fourier direction.
  //
  // The total number of processes is a command-line argument to the
  // message-passing manager (e.g. mpirun) and pre-exists at time of
  // entry to this routine, while the requested number of 2D
  // partitions is a command-line argument to the semtex executable it
  // manages (e.g. dns_mp), and is supplied as input argument npart2d.
  //
  // We check that npart2d is a factor of the available number of
  // processes, but further checking will likely need to be carried
  // out elsewhere: e.g. that npart2d doesn't exceed the number of
  // elements and npartz is consistent with the FFT we use in the
  // Fourier direction.
  // -----------------------------------------------------------------------
  {
#if defined(MPI_EX)

    const char routine[] = "Message::grid";

    int ntot, dim_sizes[2], wrap_around[] = {0, 0}, free_coords[2], reorder = 1;
    
    if (npart2d < 1)
      Veclib::alert
    	(routine, "no. of 2D partitions must be at least 1", ERROR);
    
    MPI_Comm_size (MPI_COMM_WORLD, &ntot);

    if (ntot % npart2d)
      Veclib::alert
    	(routine, "no. of 2D partitions must factor no. of processes", ERROR);

#if defined(XXT_EX)
    
    // -- Will allow 2D mesh partitioning (rows) as well as being
    //    parallel across Fourier modes (columns).  But, potentially,
    //    the number of Fourier modes can be only 1. So we may need a
    //    1D Cartesian MPI grid with only a single row/column (i.e. a
    //    line).

    if (ntot == npart2d) { // -- Solution is 2D.

      // -- Make a row communicator.

      dim_sizes[0]   = ntot;
  
      MPI_Cart_create (MPI_COMM_WORLD,
		       1, dim_sizes, wrap_around, reorder, &grid_comm);

      free_coords[0] = 1;

      MPI_Cart_sub   (grid_comm, free_coords, &row_comm);

      MPI_Comm_rank  (row_comm,  &ipart2d);

      npartz = 1;
      ipartz = 0;
      
    } else {		  // -- 3D solution with both 2D and z partitions.

      // -- Create both row and column communicators.

      int grid_rank, coordinates[2];

      dim_sizes[0] = npart2d;
      dim_sizes[1] = npartz = ntot / npart2d;
  
      MPI_Cart_create  (MPI_COMM_WORLD,
			2, dim_sizes, wrap_around, reorder, &grid_comm);

      free_coords[0] = 0;
      free_coords[1] = 1;

      MPI_Cart_sub    (grid_comm, free_coords, &row_comm);

      free_coords[0] = 1;
      free_coords[1] = 0;

      MPI_Cart_sub    (grid_comm, free_coords, &col_comm);

      MPI_Comm_rank   (grid_comm, &grid_rank);
      MPI_Cart_coords (grid_comm,  grid_rank, 2, coordinates)

      ipart2d = coordinates[0];
      ipartz  = coordinates[1];
    }

#else    
    
    // -- Old-style semtex, parallel-across-Fourier modes, with only a
    //    single 2D mesh partition.  Allocate a 1D Cartesian MPI grid
    //    (say it is a column).

    if (npart2d != 1)
      Veclib::alert
    	(routine, "no. of 2D partitions for Fourier-parallel must be 1", ERROR);

    dim_sizes[0]   = ntot;
    wrap_around[0] = 0;
  
    MPI_Cart_create (MPI_COMM_WORLD,
		     1, dim_sizes, wrap_around, reorder, &grid_comm);

    free_coords[0] = 1;

    MPI_Cart_sub   (grid_comm, free_coords, &col_comm);

    MPI_Comm_rank  (col_comm,  &ipartz);

    npartz  = ntot;
    ipart2d = 0;

#endif
    
#else
    // -- Serial execution; supply default return values (not used).
    
    ipart2d = 0;
    ipartz  = 0, npartz = 1;
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

    int np;
    MPI_Comm_size (col_comm, &np);

    int            i, j;
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

    int np;
    MPI_Comm_size (col_comm, &np);

    int            i, j;
    const int      nB = nP / np;       // -- Size of intra-processor block.
    const int      NB = nP / nB;       // -- Number of these blocks in a plane.
    const int      NM = nP * nZ / np;  // -- Size of message block.
    const int      dsize  = sizeof (int_t);
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
