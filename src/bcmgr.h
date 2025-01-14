#ifndef BCMGR_H
#define BCMGR_H

typedef struct bctriple { char group; int_t elmt; int_t side; } BCtriple;

class BCmgr
// ===========================================================================
// This is a factory / retrieval service for classes derived from
// Condition, and maintains GROUP descriptors.  In addition, it reads
// and returns NumberSys objects from session.num.  Also, it contains
// code for maintenance and evaluation of computed BC types.
//
// REFERENCES
// ----------
// [1] Karniadakis, Israeli & Orszag, JCP 97:414-443 (1991)
// [2] Blackburn & Sherwin, JCP 197:759-778 (2004)
// [3] Dong, JCP 302:300-328 (2015)
// [4] Dong, Karniadakis & Chryssostomidis, JCP 261:83-105 (2014)
// [5] Liu, Xie & Dong, IJHFF 151:119355 (2020)
// ===========================================================================
{
public:
  BCmgr (FEML*, vector<Element*>&);

  const char*        field        () const { return _fields; }
  const char*        groupInfo    (const char) const;
  const Condition*   getCondition (const char, const char, const int_t = 0);
  vector<BCtriple*>& getBCedges   () { return _elmtbc; }
  int_t              nBCedges     () const { return _elmtbc.size(); }
  int_t              nWall        (); // Why not const: OSX compiler bug?
  int_t              nAxis        ();
  int_t              nMatching    (const char*);

  class CondRecd {
  public: 
    char       grp  ;
    char       fld  ;
    Condition* bcn  ;
    char*      value;
  };

  // -- Routines for maintaining and evaluating computed BCs.

  // -- Chicken & egg: buildComputedBCs can't be in class constructor
  //    because we need a Field to do it.  Right?

  void buildComputedBCs (const Field*, const bool = false);

  void maintainFourier  (const int_t, const Field*, const AuxField**,
			 const AuxField**, const int_t, const int_t,
			 const bool = true);  
  void maintainPhysical (const Field*, const vector<AuxField*>&,
			 const int_t, const int_t);
  void evaluateCNBCp    (const int_t, const int_t, const int_t, real_t*);
  void evaluateCMBCp    (const Field*, const int_t, const int_t, 
			 const int_t, real_t*);
  void evaluateCMBCu    (const Field*, const int_t, const int_t, 
			 const int_t, const char, real_t*);
  void evaluateCMBCc    (const Field*, const int_t, const int_t,
			 const int_t, real_t*);
  void accelerate       (const Vector&, const Field*);

private:
  char*              _fields  ; // String containing field names.
  vector<char>       _group   ; // Single-character group tags.
  vector<char*>      _descript; // Group name strings.
  vector<CondRecd*>  _cond    ; // Conditions in storage.
  vector<BCtriple*>  _elmtbc  ; // Group tags for each element-side BC.
  bool               _axis    ; // Session file declared an axis BC group.
  bool               _open    ; // Session file declared an open BC group.

  void buildnum  (const char*, vector<Element*>&);
  void buildsurf (FEML*, vector<Element*>&);

  // -- Storage of past-time values needed for computed BCs:

  int_t     _nLine;     // Same as for Field storage.
  int_t     _nEdge;     // Number of edges with BCs attached.
  int_t     _nZ;        // Same as for Field storage.
  int_t     _nP;        // Geometry::nP(), number of points on element edge.
  int_t     _nTime;     // N_TIME, time stepping order.

  real_t*** _u;         // (Physical) x velocity component.
  real_t*** _v;         // (Physical) y velocity component.
  real_t*** _w;         // (Physical) z velocity component.
  real_t*** _c;         // (Physical) scalar.  

  real_t*** _uhat;      // (Fourier)  x velocity component.
  real_t*** _vhat;      // (Fourier)  y velocity component.
  real_t*** _what;      // (Fourier)  z velocity component.
  real_t*** _chat;      // (Fourier)  scalar.    

  real_t*** _un;	// (Fourier)  normal velocity u.n for d(u.n)/dt.
  real_t*** _divu;	// (Fourier)  KINVIS * div(u).
  real_t*** _gradu;	// (Fourier)  KINVIS * (normal gradient of velocity).n.
  real_t*** _hopbc;	// (Fourier)  normal component of [N(u)+f-curlCurl(u)].
  real_t*** _ndudt;     // (Fourier)  normal component of (partial) du/dt.

  real_t*   _work;      // Computational workspace (scratch).
  real_t*   _Theta;     // Inflow switch function.
  real_t*   _u2;	// (Physical) 2 * kinetic energy (u^2+v^2+w^2)*.
  real_t*   _unp;       // (Physical) u*.n.
  real_t*   _Enux;	// Energy flux traction, x component. See Ref [3].
  real_t*   _Enuy;	// Energy flux traction, y component.
  real_t*   _H;		// Scalar "energy" fluxion.  See Ref [5].

  bool      _toggle;    // Toggle switch for Fourier transform of KE.
  bool      _scalar;	// Problem has a scalar as well as velocity.  
};

#endif
