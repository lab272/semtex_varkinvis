#ifndef STAB_H
#define STAB_H

//////////////////////////////////////////////////////////////////////////////
// stab.h: header file for linearised NS stability solver.
//
// $Id$
//////////////////////////////////////////////////////////////////////////////

#include <sem.h>
#include <data2df.h>

typedef enum { PRIMAL, ADJOINT, GROWTH, SHRINK } problem_t;


class StabAnalyser : public Analyser
// ===========================================================================
// Implement step-by-step processing and output control for flow solver.
// ===========================================================================
{
public:
  StabAnalyser (Domain*, FEML*);
  void analyse (AuxField**, AuxField** = 0);

private:
  vector<HistoryPoint*> base_history; // -- Locations, etc. of history points.
  ofstream              bhs_strm;     // -- File for base history points.
};

void integrate  (void (*)(Domain*, BCmgr*, AuxField**, AuxField**),
		 Domain*, BCmgr*, StabAnalyser*);
void linAdvect  (Domain*, BCmgr*, AuxField**, AuxField**);
void linAdvectT (Domain*, BCmgr*, AuxField**, AuxField**);

// -- Fortran interfaces to things not included in semtex.

extern "C" {
#if defined (ARPACK)

// -- Nonsymmetric eigensystem.

  void F77NAME(dnaupd) 		// -- ARPACK reverse-communications interface.
    (int_t&         ido   ,
     const char*    bmat  ,
     const int_t&   n     ,
     const char*    which ,
     const int_t&   nev   ,
     real_t&        tol   ,
     real_t*        resid ,
     const int_t&   ncv   ,
     real_t*        v     ,
     const int_t&   ldv   ,
     int_t*         iparam,
     int_t*         ipntr ,
     real_t*        workd ,
     real_t*        workl ,
     const int_t&   lworkl,
     int_t&         info  );
  void F77NAME(dneupd)		// -- Postprocessing.
    (const int_t&   rvec  ,
     const char*    howmny,
     const int*     select,
     real_t*        dr    ,
     real_t*        di    ,
     real_t*        z     ,
     const int_t&   ldz   ,
     const real_t&  sigmar,
     const real_t&  sigmai,
     real_t*        workev,
     const char*    bmat  ,	// -- Remainder unchanged after dnaupd.
     const int_t&   n     ,
     const char*    which ,
     const int_t&   nev   ,
     real_t&        tol   ,
     real_t*        resid ,
     const int_t&   ncv   ,
     real_t*        v     ,
     const int_t&   ldv   ,
     int_t*         iparam,
     int_t*         ipntr ,
     real_t*        workd ,
     real_t*        workl ,
     const int_t&   lworkl,
     int_t&         info  );

// -- Symmetric eigensystem.

  void F77NAME(dsaupd) 		// -- ARPACK reverse-communications interface.
    (int_t&         ido   ,
     const char*    bmat  ,
     const int_t&   n     ,
     const char*    which ,
     const int_t&   nev   ,
     real_t&        tol   ,
     real_t*        resid ,
     const int_t&   ncv   ,
     real_t*        v     ,
     const int_t&   ldv   ,
     int_t*         iparam,
     int_t*         ipntr ,
     real_t*        workd ,
     real_t*        workl ,
     const int_t&   lworkl,
     int_t&         info  );
  void F77NAME(dseupd)		// -- Postprocessing.
    (const int_t&   rvec  ,
     const char*    howmny,
     const int*     select,
     real_t*        dr    ,
     real_t*        z     ,
     const int_t&   ldz   ,
     const real_t&  sigma ,
     const char*    bmat  ,	// -- Remainder unchanged after dsaupd.
     const int_t&   n     ,
     const char*    which ,
     const int_t&   nev   ,
     real_t&        tol   ,
     real_t*        resid ,
     const int_t&   ncv   ,
     real_t*        v     ,
     const int_t&   ldv   ,
     int_t*         iparam,
     int_t*         ipntr ,
     real_t*        workd ,
     real_t*        workl ,
     const int_t&   lworkl,
     int_t&         info  );     
#endif
}


#endif
