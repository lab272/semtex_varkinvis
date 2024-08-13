#ifndef FIELD_H
#define FIELD_H
  

class Field : public AuxField
// ===========================================================================
// Field adds boundary conditions and global numbering to AuxField.
// With the boundary conditions comes knowledge of which sections of
// the boundary hold essential (known) nodes, and which contain
// natural (unknown) nodes.
//
// A Field holds global node numbers and solve masks for element
// boundaries: where mesh value is given by an essential BC, solve
// mask is 0, for all other nodes have value 1.
//
// Helmholtz solution routines are also available.
// ===========================================================================
{
friend class BCmgr;

public:
  Field  (real_t*, BoundarySys*, NumberSys*, const int_t,
	  vector<Element*>&, const char);
 ~Field  () { }

  Field& operator = (const AuxField& z) {AuxField::operator=(z); return *this;}
  Field& operator = (const real_t&   z) {AuxField::operator=(z); return *this;}
  Field& operator = (const char*     z) {AuxField::operator=(z); return *this;}

  Field& solve  (AuxField*, const ModalMatrixSys*);
  Field& solve  (AuxField*, const MatrixSys*, AuxField*);

  void evaluateBoundaries    (const Field*, const int_t, const bool = true);
  void evaluateM0Boundaries  (const Field*, const int_t);

  static void   coupleBCs    (Field*, Field*, const int_t);
  static real_t modeConstant (const char, const int_t, const real_t);

private:
  int_t        _nbound;		// Number of boundary edges.
  int_t        _nline ;		// Length of one boundary line.
  real_t*      _sheet ;		// Wrap-around storage for data boundary.
  real_t**     _line  ;		// Single plane's worth of sheet.
  BoundarySys* _bsys  ;		// Boundary system information.
  NumberSys*   _nsys  ; 	// Assembly mapping information.

  void getEssential      (const real_t*, real_t*,
			  const vector<Boundary*>&, const AssemblyMap*) const;
  void setEssential      (const real_t*, real_t*, const AssemblyMap*);
  void local2global      (const real_t*, real_t*, const AssemblyMap*)   const;
  void global2local      (const real_t*, real_t*, const AssemblyMap*)   const;
  void local2globalSum   (const real_t*, real_t*, const AssemblyMap*)   const;

  void constrain         (real_t*, const real_t, AuxField*, const real_t,
			  const real_t*, const AssemblyMap*, real_t*)   const;
  void buildRHS          (real_t*, const real_t*, real_t*, real_t*,
			  const real_t**, const int_t, const int_t,
			  const vector<Boundary*>&,
			  const AssemblyMap*, real_t*)                  const;
  void HelmholtzOperator (const real_t*, real_t*, const real_t, AuxField*,
			  const real_t, const int_t, real_t*)           const;
};

#endif
