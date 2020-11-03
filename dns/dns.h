#ifndef DNS_H
#define DNS_H
//////////////////////////////////////////////////////////////////////////////
// dns.h: header file for direct numerical simulation solver.
//
// Copyright (c) 1994 <--> $Date$, Hugh Blackburn
//
// $Id$
//////////////////////////////////////////////////////////////////////////////

#include <sem.h>
#include <fieldforce.h>

class DNSAnalyser : public Analyser
// ===========================================================================
// Implement step-by-step processing and output control for flow solver.
// ===========================================================================
{
public:
  DNSAnalyser  (Domain*, BCmgr*, FEML*);
  void analyse (AuxField**, AuxField**);

private:
  ofstream       _flx_strm;

  bool           _wall; 	// -- True if "wall" boundary group exists.
  bool           _wss;          // -- True if wall shear stress requested.
  ofstream       _wss_strm;
  int_t          _nline;
  int_t          _nwall;
  int_t          _npad;

  vector<real_t> _work;
};

void skewSymmetric    (Domain*,BCmgr*,AuxField**,AuxField**,FieldForce*);
void altSkewSymmetric (Domain*,BCmgr*,AuxField**,AuxField**,FieldForce*);
void convective       (Domain*,BCmgr*,AuxField**,AuxField**,FieldForce*);
void rotational1      (Domain*,BCmgr*,AuxField**,AuxField**,FieldForce*);
void rotational2      (Domain*,BCmgr*,AuxField**,AuxField**,FieldForce*);
void Stokes           (Domain*,BCmgr*,AuxField**,AuxField**,FieldForce*);

#endif
