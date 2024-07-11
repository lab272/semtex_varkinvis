#ifndef DOMAIN_H
#define DOMAIN_H


class Domain
// ===========================================================================
// Physical domain storage class.
//
// Since users are expected to program with entities kept at the
// Domain level, all internal storage is exposed to view.
//
// Domain fields are written/read untransformed (in physical space).
// ===========================================================================
{
friend ifstream& operator >> (ifstream&, Domain&);
friend ofstream& operator << (ofstream&, Domain&);

public:
  Domain (FEML*, const Mesh*, vector<Element*>&, BCmgr*);

  char*                name;	// Session name.
  char*                field;	// List of lower-case Field names.
  int_t                step;	// Runtime step number.
  real_t               time;	// Simulation time.
  vector<Element*>&    elmt;	// Shared for equal-order interpolations.
  vector<real_t*>      udat;	// Data storage area for solution fields.
  vector<Field*>       u   ;	// Solution fields: velocities, pressure.
  vector<BoundarySys*> b   ;	// Corresponding boundary systems.
  vector<NumberSys*>   n   ;    // Corresponding numbering systems.  

  int_t nField    () const { return u.size(); }
  void  report    (ostream& stream = cout);
  bool  restart   ();
  void  dump      ();

  // -- Required for base fields and stability analysis.

  vector<AuxField*> U       ; // -- Base velocity fields - no BCs.
  vector<real_t*>   Udat    ; // -- Data storage area for base auxfields.
  vector<real_t*>   baseFlow; // -- Fourier transformed base velocities.
  real_t            period  ; // -- Temporal period of base flow (if relevant).

  // -- Required for the kinvis field.
  
  AuxField* VARKINVIS;  // -- Non-scalar kinematic viscosity field.
  real_t* varkinvisdat; // -- Data storage area for kinematic viscosity auxfield.
  
  void loadBase  ();
  void updateBase();

  int_t         nGlobal       () const { return _nglobal;        }
  const int_t*  assemblyNaive () const { return &_bmapNaive[0];  } 
  const real_t* invMassNaive  () const { return &_imassNaive[0]; }

private:
  char* baseField;    // Upper-case single character base velocity field names.

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
