#ifndef CONDITION_H
#define CONDITION_H

class Condition
// ===========================================================================
// Virtual base class for boundary condition application. The base
// class is made pure virtual so that identical calls may be used for
// all BCs but only the needed ones of for each concrete derived class
// will be applied.
//
// Each concrete class is derived from the virtual base class Condition:
// 1. Essential         : essential BC with constant, supplied, value.
// 2. EssentialFunction : essential BC, value obtained by parsing a function.
// 3. Natural           : natural BC with constant, supplied, value.
// 4. NaturalFunction   : natural BC, value obtained by parsing a function.
// 5. Mixed             : transfer coefficient type, 2 supplied values.
// 6. NaturalCBCp       : "high-order" pressure BC, natural, computed value.
// 7. MixedCBCp         : Open boundary pressure BC, mixed/Robin, computed.
// 8. MixedCBCu         : Open boundary u velocity BC, mixed, computed.
// 9. MixedCBCv         : Open boundary v velocity BC, mixed, computed.
// 10. MixedCBCw        : Open boundary w velocity BC, mixed, computed.
// 11. MixedCBCc        : Open boundary scalar c BC, mixed, computed.
//
// Note that for supplied-value BC types, the value is a physical-space
// description of the BC, and is now set once at run-time (cannot be reset).
// This is not true for those obtained by parsing a function, which are
// re-parsed every timestep.
//
// Also, each condition class derived from the base has to define all the
// pure virtual functions listed below (except the destructor) but
// some of these will just be stubs that do nothing for any particular
// type.  Those stubs are indicated in the present header
// file with the function body "{ };".
//
// All condition classes must provide an "evaluate" method: this is
// used to install values in Field class BC storage area.  For
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
// See also boundary.h, edge.h, bcmgr.cpp mesh.cpp.
// ===========================================================================
{
public:
  virtual void evaluate  (const Field*, const int_t, const int_t,
			  const Element*, const int_t, const int_t,
			  const bool, real_t*)                         const=0;
  virtual void set       (const int_t, const int_t*,
			  const real_t*, real_t*)                      const=0;
  virtual void sum       (const int_t, const int_t*,
		          const real_t*,const real_t*,real_t*,real_t*) const=0;
  virtual void augmentSC (const int_t,  const int_t, const int_t,
			  const int_t*,const real_t*,real_t*,real_t*)  const=0;
  virtual void augmentOp (const int_t, const int_t*,
			  const real_t*, const real_t*, real_t*)       const=0;
  virtual void augmentDg (const int_t, const int_t*, 
			  const real_t*, real_t*)                      const=0;
  virtual void describe  (char* tgt)                                   const=0;

  virtual ~Condition() = default;
};


class Essential : public Condition
// ===========================================================================
// Essential BC applicator.  This one is for plain (constant value)
// Dirichlet/essential BCs.
// ===========================================================================
{
public:
  Essential      (const char* v) : _value (strtod (v, 0)) { }
  void evaluate  (const Field*, const int_t, const int_t,
		  const Element*, const int_t, const int_t,
		  const bool, real_t*)                           const override;
  void set       (const int_t, const int_t*,
		  const real_t*, real_t*)                        const override;
  void describe  (char*)                                         const override;

  // -- Other methods are empty stubs for essential BCs.
  
  void sum       (const int_t, const int_t*, const real_t*,
		  const real_t*, real_t*, real_t*)               const { };
  void augmentSC (const int_t, const int_t, const int_t,
		  const int_t*, const real_t*, real_t*, real_t*) const { };
  void augmentOp (const int_t, const int_t*,
		  const real_t*, const real_t*, real_t*)         const { };
  void augmentDg (const int_t, const int_t*, 
		  const real_t*, real_t*)                        const { };
private:
  real_t _value;
};


class EssentialFunction : public Condition
// ===========================================================================
// Essential BC applicator for specified function Dirichlet/essential BCs.
// This uses parser to set the BC values. 
// ===========================================================================
{
public:
  EssentialFunction      (const char*);
  virtual void evaluate  (const Field*,   const int_t, const int_t,
                          const Element*, const int_t, const int_t,
			  const bool, real_t*)                           const;
  virtual void set       (const int_t, const int_t*,
			  const real_t*, real_t*)                        const;
  virtual void sum       (const int_t, const int_t*,
			  const real_t*,const real_t*,real_t*,real_t*)   const
    { };
  virtual void augmentSC (const int_t, const int_t, const int_t,
			  const int_t*, const real_t*, real_t*, real_t*) const
    { };
  virtual void augmentOp (const int_t, const int_t*,
			  const real_t*, const real_t*, real_t*)         const
    { };
  virtual void augmentDg (const int_t, const int_t*, 
			  const real_t*, real_t*)                        const
    { };
  virtual void describe  (char*)                                         const;
private:
  char* _function;
};


class Natural : public Condition
// ===========================================================================
// Natural BC applicator.  This one is for plain (constant value)
// Neumann/natural BCs.
// ===========================================================================
{
public:
  Natural                (const char* v) : _value (strtod (v, 0)) { }
  virtual void evaluate  (const Field*,   const int_t, const int_t,
                          const Element*, const int_t, const int_t,
			  const bool, real_t*)                           const;
  virtual void set       (const int_t, const int_t*,
			  const real_t*, real_t*)                        const
    { };
  virtual void sum       (const int_t, const int_t*,
			  const real_t*,const real_t*,real_t*,real_t*)   const;
  virtual void augmentSC (const int_t, const int_t, const int_t,
			  const int_t*, const real_t*, real_t*, real_t*) const
    { };
  virtual void augmentOp (const int_t, const int_t*,
			  const real_t*, const real_t*, real_t*)         const
    { };
  virtual void augmentDg (const int_t, const int_t*, 
			  const real_t*, real_t*)                        const
    { };
  virtual void describe  (char*)                                         const;
private:
  real_t _value;
};


class NaturalFunction : public Condition
// ===========================================================================
// Natural BC applicator for specified function Neumann/natural BCs.
// This uses parser to set the BC values. 
// ===========================================================================
{
public:
  NaturalFunction        (const char*);
  virtual void evaluate  (const Field*,   const int_t, const int_t,
                          const Element*, const int_t, const int_t,
			  const bool, real_t*)                           const;
  virtual void set       (const int_t, const int_t*,
			  const real_t*, real_t*)                        const
    { };
  virtual void sum       (const int_t, const int_t*,
			  const real_t*,const real_t*,real_t*,real_t*)   const;
  virtual void augmentSC (const int_t, const int_t, const int_t,
			  const int_t*, const real_t*, real_t*, real_t*) const
    { };
  virtual void augmentOp (const int_t, const int_t*,
			  const real_t*, const real_t*, real_t*)         const
    { };
  virtual void augmentDg (const int_t, const int_t*, 
			  const real_t*, real_t*)                        const
    { };
  virtual void describe  (char*)                                         const;
private:
  char* _function;
};


class Mixed : public Condition
// ===========================================================================
// Boundary condition class for mixed (a.k.a. Robin) type BCs of form
//     dc/dn + K(c - C) = 0.
// Mixed BCs affect problem Helmholtz matrices, but only on the diagonal, 
// and element-boundary, terms. Syntax in session file is
//     <M> c = K, C </M>  or 
//     <M> c = K; C </M> 
// where 'c' is a field name and K and C can be evaluated as constants
// (perhaps using defined TOKENS). White space following the separators
// above (',' and ';') preceding C is optional.
// ===========================================================================
{
public:
  Mixed                  (const char*);
  virtual void evaluate  (const Field*,   const int_t, const int_t,
                          const Element*, const int_t, const int_t,
			  const bool, real_t*)                           const
    { };
  virtual void set       (const int_t, const int_t*,
			  const real_t*, real_t*)                        const
    { };
  virtual void sum       (const int_t, const int_t*, const real_t*,
			  const real_t*, real_t*, real_t*)               const;
  virtual void augmentSC (const int_t, const int_t, const int_t,
			  const int_t*, const real_t*, real_t*, real_t*) const;
  virtual void augmentOp (const int_t, const int_t*,
			  const real_t*, const real_t*, real_t*)         const;
  virtual void augmentDg (const int_t, const int_t*, 
			  const real_t*, real_t*)                        const;
  virtual void describe  (char*)                                         const;
private:
  real_t _K_;		// -- This is "K" above.
  real_t _C_;		// -- This is "C" above.
};


// -- Classes with internally computed BC ("CBC") types follow.
//    Evaluate (and other) functions dealt with in bcmgr.cpp.


class NaturalCBCp : public Condition
// ===========================================================================
// Computed Neumann BC for pressure, typical of wall boundaries.  This is the
// prototype High-Order Pressure BC (HOPBC).
//
// Karniadakis, Israeli & Orszag JCP 97, (1991).
// ===========================================================================
{
public:
  NaturalCBCp            (BCmgr* B) : _BCmgr (B) { }
  virtual void evaluate  (const Field*,   const int_t, const int_t,
                          const Element*, const int_t, const int_t,
			  const bool, real_t*)                           const;
  virtual void set       (const int_t, const int_t*,
			  const real_t*, real_t*)                        const
    { };
  virtual void sum       (const int_t, const int_t*,
			  const real_t*,const real_t*,real_t*,real_t*)   const;
  virtual void augmentSC (const int_t, const int_t, const int_t,
			  const int_t*, const real_t*, real_t*, real_t*) const
    { };
  virtual void augmentDg (const int_t, const int_t*, 
			  const real_t*, real_t*)                        const
    { };
  virtual void augmentOp (const int_t, const int_t*,
			  const real_t*, const real_t*, real_t*)         const
    { };
  virtual void describe  (char*)                                         const;
private:
  BCmgr* _BCmgr; 
};


class MixedCBCp : public Condition
// ===========================================================================
// Computed mixed BC for pressure on open boundaries.
//
// Dong (2015), JCP 302:300-328.
// ===========================================================================
{
public:
  MixedCBCp              (BCmgr* B);
  virtual void evaluate  (const Field*,   const int_t, const int_t,
                          const Element*, const int_t, const int_t,
			  const bool, real_t*)                           const;
  virtual void set       (const int_t, const int_t*,
			  const real_t*, real_t*)                        const
    { };
  virtual void sum       (const int_t, const int_t*,
			  const real_t*,const real_t*,real_t*,real_t*)   const;
  virtual void augmentSC (const int_t, const int_t, const int_t,
			  const int_t*, const real_t*, real_t*, real_t*) const;
  virtual void augmentDg (const int_t, const int_t*, 
			  const real_t*, real_t*)                        const;
  virtual void augmentOp (const int_t, const int_t*,
			  const real_t*, const real_t*, real_t*)         const;
  virtual void describe  (char*)                                         const;
private:
  BCmgr*  _BCmgr;
  real_t  _K_;
};


class MixedCBCu : public Condition
// ===========================================================================
// Computed mixed BC for velocity component 'u' on open boundaries.
//
// Dong (2015), JCP 302:300-328.
// ===========================================================================
{
public:
  MixedCBCu              (BCmgr* B);
  virtual void evaluate  (const Field*,   const int_t, const int_t,
                          const Element*, const int_t, const int_t,
			  const bool, real_t*)                           const;
  virtual void set       (const int_t, const int_t*,
			  const real_t*, real_t*)                        const
    { };
  virtual void sum       (const int_t, const int_t*,
			  const real_t*,const real_t*,real_t*,real_t*)   const;
  virtual void augmentSC (const int_t, const int_t, const int_t,
			  const int_t*, const real_t*, real_t*, real_t*) const;
  virtual void augmentDg (const int_t, const int_t*, 
			  const real_t*, real_t*)                        const;
  virtual void augmentOp (const int_t, const int_t*,
			  const real_t*, const real_t*, real_t*)         const;
  virtual void describe  (char*)                                         const;
private:
  BCmgr*  _BCmgr;
  real_t* _alpha;
  real_t  _K_, _DoDt;
};


class MixedCBCv : public Condition
// ===========================================================================
// Computed mixed BC for velocity component 'v' on open boundaries.
//
// Dong (2015), JCP 302:300-328.
// ===========================================================================
{
public:
  MixedCBCv              (BCmgr* B);
  virtual void evaluate  (const Field*,   const int_t, const int_t,
                          const Element*, const int_t, const int_t,
			  const bool, real_t*)                           const;
  virtual void set       (const int_t, const int_t*,
			  const real_t*, real_t*)                        const
    { };
  virtual void sum       (const int_t, const int_t*,
			  const real_t*,const real_t*,real_t*,real_t*)   const;
  virtual void augmentSC (const int_t, const int_t, const int_t,
			  const int_t*, const real_t*, real_t*, real_t*) const;
  virtual void augmentDg (const int_t, const int_t*, 
			  const real_t*, real_t*)                        const;
  virtual void augmentOp (const int_t, const int_t*,
			  const real_t*, const real_t*, real_t*)         const;
  virtual void describe  (char*)                                         const;
private:
  BCmgr*  _BCmgr;
  real_t* _alpha;
  real_t  _K_, _DoDt;
};


class MixedCBCw : public Condition
// ===========================================================================
// Computed mixed BC for velocity component 'w' on open boundaries.
//
// Dong (2015), JCP 302:300-328.
// ===========================================================================
{
public:
  MixedCBCw              (BCmgr* B);
  virtual void evaluate  (const Field*,   const int_t, const int_t,
                          const Element*, const int_t, const int_t,
			  const bool, real_t*)                           const;
  virtual void set       (const int_t, const int_t*,
			  const real_t*, real_t*)                        const
    { };
  virtual void sum       (const int_t, const int_t*,
			  const real_t*,const real_t*,real_t*,real_t*)   const;
  virtual void augmentSC (const int_t, const int_t, const int_t,
			  const int_t*, const real_t*, real_t*, real_t*) const;
  virtual void augmentDg (const int_t, const int_t*, 
			  const real_t*, real_t*)                        const;
  virtual void augmentOp (const int_t, const int_t*,
			  const real_t*, const real_t*, real_t*)         const;
  virtual void describe  (char*)                                         const;
private:
  BCmgr*  _BCmgr;
  real_t* _alpha;
  real_t  _K_, _DoDt;
};


class MixedCBCwIn : public Condition
// ===========================================================================
// Computed mixed BC for velocity component 'w' on inlet boundaries.
// This version omits adding in of boundary integral terms: evaluate and sum
// are both stubs.
//
// Dong (2015), JCP 302:300-328.
// ===========================================================================
{
public:
  MixedCBCwIn            (BCmgr* B);
  virtual void evaluate  (const Field*,   const int_t, const int_t,
                          const Element*, const int_t, const int_t,
			  const bool, real_t*)                           const
    { };
  virtual void set       (const int_t, const int_t*,
			  const real_t*, real_t*)                        const
    { };
  virtual void sum       (const int_t, const int_t*,
			  const real_t*,const real_t*,real_t*,real_t*)   const
    { };
  virtual void augmentSC (const int_t, const int_t, const int_t,
			  const int_t*, const real_t*, real_t*, real_t*) const;
  virtual void augmentDg (const int_t, const int_t*, 
			  const real_t*, real_t*)                        const;
  virtual void augmentOp (const int_t, const int_t*,
			  const real_t*, const real_t*, real_t*)         const;
  virtual void describe  (char*)                                         const;
private:
  BCmgr*  _BCmgr;
  real_t* _alpha;
  real_t  _K_, _DoDt;
};


class MixedCBCc : public Condition
// ===========================================================================
// Computed mixed BC for scalar on open boundaries.
//
// Liu, Xie & Dong (2020), IJHFF 151.
// ===========================================================================
{
public:
  MixedCBCc              (BCmgr* B);
  virtual void evaluate  (const Field*, const int_t, const int_t,
                          const Element*, const int_t, const int_t,
			  const bool, real_t*)                           const;
  virtual void set       (const int_t, const int_t*,
			  const real_t*, real_t*)                        const
    { };
  virtual void sum       (const int_t, const int_t*,
			  const real_t*,const real_t*,real_t*,real_t*)   const;
  virtual void augmentSC (const int_t, const int_t, const int_t,
			  const int_t*, const real_t*, real_t*, real_t*) const;
  virtual void augmentDg (const int_t, const int_t*, 
			  const real_t*, real_t*)                        const;
  virtual void augmentOp (const int_t, const int_t*,
			  const real_t*, const real_t*, real_t*)         const;
  virtual void describe  (char*)                                         const;
private:
  BCmgr*  _BCmgr;
  real_t* _alpha;
  real_t  _K_, _DoDt;
};

#endif
