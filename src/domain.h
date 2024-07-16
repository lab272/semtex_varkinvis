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
  Domain (FEML*, const Mesh*, vector<Element*>&, BCmgr*);

  char*                name;  // Session name.
  char*                field; // List of lower-case character field names.
  int_t                step;  // Runtime step number.
  real_t               time;  // Simulation time.
  vector<Element*>&    elmt;  // Shared for equal-order interpolations.
  vector<real_t*>      udat;  // Data storage area for solution fields.
  vector<Field*>       u   ;  // Solution fields: velocities, scalar, pressure.
  vector<BoundarySys*> b   ;  // Corresponding boundary systems.
  vector<NumberSys*>   n   ;  // Corresponding numbering systems.

  int_t nField     () const { return u.size(); }
  int_t nAdvect    () const { return u.size() - 1; } // No. of advected terms.
  int_t nVelCmpt   () const { return                 // "" velocity components.
                              (hasScalar()) ? u.size() - 2 : u.size() - 1; }

  bool  hasScalar  () const { return strchr (field, 'c'); }
  void  report     ();
  void  restart    ();
  void  dump       ();
  void  transform  (const int_t);
  
  int_t         nGlobal       () const { return _nglobal;        }
  const int_t*  assemblyNaive () const { return &_bmapNaive[0];  } 
  const real_t* invMassNaive  () const { return &_imassNaive[0]; }
  
  AuxField* VARKINVIS;  // -- Non-scalar kinematic viscosity field.
  real_t* varkinvisdat; // -- Data storage area for kinematic viscosity auxfield.

private:
  void  checkVBCs        (FEML*, const char*)         const;
  void  checkAxialBCs    (FEML*, char)                const;
  char  axialTag         (FEML*)                      const;
  bool  multiModalBCs    (FEML*, BCmgr*, const char*) const;
  void  makeAssemblyMaps (FEML*, const Mesh*, BCmgr*);

  int_t                _nglobal;     // Number of unique element-edge nodes.
  vector<int_t>        _bmapNaive;   // BC-agnostic assembly map.
  vector<real_t>       _imassNaive;  // Corresp. inverse mass matrix, _nglobal.
  vector<AssemblyMap*> _allMappings; // Complete set of domain AssemblyMaps.
};

#endif
