#ifndef CONDITION_H
#define CONDITION_H

class Condition
// ===========================================================================
// Abstract base class for boundary condition application. The base
// class is almost made pure virtual so that identical calls may be
// used for all BCs but only the needed ones of for each concrete
// derived class will be applied.
//
// There is a middle level of semi-virtual base classes for boundary
// conditions, which recognise the different ways in which the three
// basic BC types
//
// 1. Essential or Dirichlet
// 2. Natural or Neumann
// 3. Mixed or Robin
//
// are applied in finite element methods.  These inherit from the
// Condition class.
//
// Terminal classes in turn derive from the three semi-abstract base
// classes listed above.
// 
// 1. EssentialConstant : essential BC with constant, supplied, value.
// 2. EssentialFunction : essential BC, value obtained by parsing a function.
// 3. NaturalConstant   : natural BC with constant, supplied, value.
// 4. NaturalFunction   : natural BC, value obtained by parsing a function.
// 5. NaturalComputed   : "high-order" pressure BC, natural, computed value.
// 6. MixedConstant     : transfer coefficient type, 2 supplied values.
// 7. MixedComputed     : Open boundary BCs, mixed/Robin, computed.
//
// Note that for supplied-constant-value BC types, the value is a
// physical-space description of the BC, and is now set once at
// run-time (cannot be reset).  This is not true for those obtained by
// parsing a function, which are re-parsed every timestep (again, in
// physical space).
//
// All terminal condition classes must provide an "evaluate" method:
// this is used to install values in Field class BC storage area.  For
// essential/Dirichlet condition types, this is the field value that
// is lifted out of the solution (imposed).  For natural/Neumann
// condition types, this is the desired normal derivative of the field
// value along the boundary, which is not directly imposed, but
// approximated by an integral.  This evaluation could take place in
// either physical or Fourier space.
//
// For essential/Dirichlet BCs, the method "set" must be defined;
// For natural/Neumann BCs, the method "sum" must be defined, while
// For mixed/Robin BCs, the three methods "augmentSC", "augmentOp", 
// and "augmentDG" must be defined.
//
// Values for computed BC types are manufactured by BCmgr class.
//
// See also condition.cpp, boundary.h, edge.h, bcmgr.cpp mesh.cpp.
// ===========================================================================
{
public:
  virtual void evaluate  (const Field*, const int_t, const int_t,
			  const Element*, const int_t, const int_t,
			  const bool, real_t*)                        const = 0;
  virtual void set       (const int_t, const int_t*,
			  const real_t*, real_t*)                     const = 0;
  virtual void sum       (const int_t, const int_t*, const real_t*,
			  const real_t*,real_t*,real_t*)              const = 0;
  virtual void augmentSC (const int_t,  const int_t, const int_t,
			  const int_t*,const real_t*,real_t*,real_t*) const = 0;
  virtual void augmentOp (const int_t, const int_t*,
			  const real_t*, const real_t*, real_t*)      const = 0;
  virtual void augmentDg (const int_t, const int_t*, 
			  const real_t*, real_t*)                     const = 0;
  
  virtual ~Condition() = default;

  void describe  (char* tgt) const { sprintf (tgt, _descriptor.c_str ()); }

protected:
  string _descriptor; 		// -- Set during construction for all types.
};


//////////////////////////////////////////////////////////////////////////////
// Next we have derived abstract base classes for the three basic
// types of BC: essential, natural, mixed.  Note that not all of the
// Condition class methods get implemented for each type.  A "final"
// declaration indicates that the method will be implemented at this
// level (and is shared by further-derived concrete classes); ones
// ending with "=0" indicate that they are still virtual and will be
// implemented at the next level of derivation; ones ending in "{}"
// indicate that they are not implemented for the class.
//////////////////////////////////////////////////////////////////////////////


class Essential : public Condition
// ===========================================================================
// Virtual base class for Essential/Dirichlet BC applicators.
// ===========================================================================
{
  void set      (const int_t, const int_t*,
		 const real_t*, real_t*)                           const final;

  void evaluate (const Field*, const int_t, const int_t,
		 const Element*, const int_t, const int_t,
		 const bool, real_t*)                              const = 0;
  
  void sum       (const int_t, const int_t*, const real_t*,
		  const real_t*, real_t*, real_t*)                 const {};
  void augmentSC (const int_t, const int_t, const int_t,
		  const int_t*, const real_t*, real_t*, real_t*)   const {};
  void augmentOp (const int_t, const int_t*,
		  const real_t*, const real_t*, real_t*)           const {};
  void augmentDg (const int_t, const int_t*, 
		  const real_t*, real_t*)                          const {};
};


class Natural : public Condition
// ===========================================================================
// Virtual base class for Essential/Neumann BC applicators.
// ===========================================================================
{
  void sum       (const int_t, const int_t*, const real_t*,
		  const real_t*, real_t*, real_t*)                 const final;
  
  void evaluate  (const Field*, const int_t, const int_t,
		  const Element*, const int_t, const int_t,
		  const bool, real_t*)                             const = 0;

  void set       (const int_t, const int_t*,
		  const real_t*, real_t*)                          const {};
  void augmentSC (const int_t,  const int_t, const int_t,
		  const int_t*, const real_t*, real_t*, real_t*)   const {};
  void augmentOp (const int_t, const int_t*,
		  const real_t*, const real_t*, real_t*)           const {};
  void augmentDg (const int_t, const int_t*, 
		  const real_t*, real_t*)                          const {};
}; 


class Mixed : public Condition
// ===========================================================================
// Virtual base class for Mixed/Robin BC applicators.
// ===========================================================================
{
  void sum       (const int_t, const int_t*, const real_t*,
		  const real_t*, real_t*, real_t*)                 const final;
  void augmentSC (const int_t, const int_t, const int_t,
		  const int_t*, const real_t*,
		  real_t*, real_t*)                                const final;
  void augmentOp (const int_t, const int_t*,
		  const real_t*, const real_t*, real_t*)           const final;
  void augmentDg (const int_t, const int_t*, 
		  const real_t*, real_t*)                          const final;

  void evaluate  (const Field*, const int_t, const int_t,
		  const Element*, const int_t, const int_t,
		  const bool, real_t*)                             const = 0;

  void set       (const int_t, const int_t*,
		  const real_t*, real_t*)                          const {};
protected:
  real_t _K_;
};

  
//////////////////////////////////////////////////////////////////////////////
// Now, the terminal concrete classes; these are specialisations of the
// three basic types.  And the only ones with constructors.
//////////////////////////////////////////////////////////////////////////////


class EssentialConstant :  public Essential
// ===========================================================================
// Essential BC applicator.  This one is for plain (constant value)
// Dirichlet/essential BCs, described by a C string.
// ===========================================================================
{
public:
  EssentialConstant (const char* v);

  void evaluate     (const Field*, const int_t, const int_t,
		     const Element*, const int_t, const int_t,
		     const bool, real_t*)                         const final;
private:
  real_t _value;
};


class EssentialFunction : public Essential
// ===========================================================================
// Essential BC applicator for specified function Dirichlet/essential
// BCs.  This uses the parser to evaluate BC values, potentially at
// each time step. The function is of C-string type.
// ===========================================================================
{
public:
  EssentialFunction (const char* f); 

  void evaluate     (const Field*, const int_t, const int_t,
		     const Element*, const int_t, const int_t,
		     const bool, real_t*)                         const final;
private:
  char* _function;
};


// -- At present, there are no implementations of internally computed
//    essential BCs.  However, this isn't the case for natural and
//    mixed BCs, see below. The internal computations are dealt with
//    by BCmgr class (since it can do them more efficiently than on an
//    edge-by-egde basis).


class NaturalConstant : public Natural
// ===========================================================================
// Natural BC applicator.  This one is for plain (constant value)
// Neumann/natural BCs.
// ===========================================================================
{
public:
  NaturalConstant (const char* v);

  void evaluate (const Field*, const int_t, const int_t,
		 const Element*, const int_t, const int_t,
		 const bool, real_t*)                            const final;
private:
  real_t _value;
};


class NaturalFunction : public Natural
// ===========================================================================
// Natural BC applicator for specified function Neumann/natural BCs.
// This uses parser to set the BC values. 
// ===========================================================================
{
public:
  NaturalFunction (const char* f);

  void evaluate   (const Field*, const int_t, const int_t,
		   const Element*, const int_t, const int_t,
		   const bool, real_t*)                          const final;
private:
  char* _function;
};


class NaturalComputed : public Natural
// ===========================================================================
// Computed Neumann BC for pressure, typical of wall boundaries.  This is the
// prototype High-Order Pressure BC (HOPBC).
//
// Karniadakis, Israeli & Orszag JCP 97, (1991).
//
// Just in case we need to expand the evaluation method in future, we add a
// character tag indentifier.  For now, that will be 'p' (for pressure).
// ===========================================================================
{
public:
  NaturalComputed (BCmgr* B, const char t = 'p');

  void evaluate   (const Field*, const int_t, const int_t,
		   const Element*, const int_t, const int_t,
		   const bool, real_t*)                          const final;
private:
  BCmgr* _BCmgr;
  char   _tag;
};


class MixedConstant : public Mixed
// ===========================================================================
// Boundary condition class for mixed (a.k.a. Robin) type BCs of form
//     dc/dn + K(c - C) = 0.
// Mixed BCs affect problem Helmholtz matrices, but only on the diagonal, 
// and element-boundary, terms. For constant type, syntax in session file is
//     <M> c = K, C </M>  or 
//     <M> c = K; C </M> 
// where 'c' is a field name and K and C can be evaluated as constants
// (perhaps using defined TOKENS). White space following the separators
// above (',' and ';') preceding C is optional.
// ===========================================================================
{
public:
  MixedConstant (const char*);
  
  void evaluate (const Field*, const int_t, const int_t,
		 const Element*, const int_t, const int_t,
		 const bool, real_t*)                             const final;
private:
  real_t _C_;		// -- This is "C" above (_K_ is in base class).
};


// -- At present, there are no implementations of MixedFunction BCs.
//    Classes with internally computed Mixed BC ("CBC") types follow.
//    Methods for evaluate are eventually dealt with by the BCmgr
//    class (see bcmgr.cpp).


class MixedComputed : public Mixed
// ===========================================================================
// Computed mixed BC for pressure on open boundaries.
//
// Dong (2015), JCP 302:300-328.
// ===========================================================================
{
public:
  MixedComputed (BCmgr*, const char);
  
  void evaluate (const Field*, const int_t, const int_t,
		 const Element*, const int_t, const int_t,
		 const bool, real_t*)                             const final;
private:
  BCmgr* _BCmgr;
  char   _tag;
};

#endif
