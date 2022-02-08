#ifndef DOMAIN_H
#define DOMAIN_H

class Domain
// ===========================================================================
// Physical domain storage class for Navier--Stokes and elliptic type
// problems.
//
// Since users are expected to program with entities kept at the
// Domain level, most internal storage is exposed to view (like a C
// struct).
//
// Domain fields are written/read from/to file untransformed (in
// physical space) and at interpolation order N_P.
//
// ===========================================================================
{
friend istream& operator >> (istream&, Domain&);
friend ostream& operator << (ostream&, Domain&);
public:
  Domain (const FEML*, Mesh const*, vector<Element*>&, BCmgr const*);

  char*                name;  // Session name.
  char*                field; // List of lower-case character field names.
  int_t                step;  // Runtime step number.
  real_t               time;  // Simulation time.
  vector<Element*>&    elmt;  // Shared for equal-order interpolations.
  vector<real_t*>      udat;  // Data storage area for solution fields.
  vector<Field*>       u   ;  // Solution fields: velocities, scalar, pressure.
  vector<BoundarySys*> b   ;  // Field boundary systems.


  int_t nField     () const { return u.size(); }
  int_t nAdvect    () const { return u.size() - 1; } // No. of advected terms.
  int_t nVelCmpt   () const { return                 // "" velocity components.
                              (hasScalar()) ? u.size() - 2 : u.size() - 1; }

  bool  hasScalar  () const { return strchr (field, 'c'); }
  void  report     ();
  void  restart    ();
  void  dump       ();
  void  transform  (const int_t);
  
  NumberSys const* getNsys (const char, const int_t) const;

private:
  void  checkVBCs        (const FEML*, const char const*) const;
  char  axialTag         (const FEML*) const;
  bool  multiModalBCs    (const FEML*, const BCmgr const*, char const*) const;
  void  makeAssemblyMaps (const FEML*, const BCmgr const*, char const*) const;

  map <const char, NumberSys*[3]> _globalNumbering;
  vector<NumberSys*>              _n;  // Unique numbering schemes;
};

#endif
