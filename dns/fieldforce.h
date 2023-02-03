#ifndef FIELDFORCE_H
#define FIELDFORCE_H


class VirtualForce
// ---------------------------------------------------------------------------
// Virtual base class for various types of body forcing.
// ---------------------------------------------------------------------------
{
public:
  void      allocStorage       (Domain*);
  AuxField* allocAuxField      (Domain*, char);
  void      readSteadyFromFile (char*, vector<AuxField*>&);

  virtual void add      (AuxField*, const int_t, vector<AuxField*>&) {};
  virtual void subtract (AuxField*, const int_t, vector<AuxField*>&) {};

protected:
  Domain*           _D;
  bool              _enabled;
  vector<AuxField*> _a;		// -- storage for pre-processed part  
};


class FieldForce
// ---------------------------------------------------------------------------
// Provides external access for applications. See e.g. calls in
// nonlinear.cpp.
// ---------------------------------------------------------------------------
{
public:
  FieldForce       (Domain*, FEML*);
  
  void addPhysical (AuxField*, AuxField*, const int_t, vector<AuxField*>&);
  void subPhysical (AuxField*, AuxField*, const int_t, vector<AuxField*>&);
  void writeAux	   (vector<AuxField*>&);

  void canonicalSteadyBoussinesq (AuxField*,
				  vector<AuxField*>&, vector<AuxField*>&);
  
protected:
  Domain*		_D;
  bool			_enabled;
  vector<VirtualForce*> _classes;   // -- Concrete body forcing classes.
private:
  real_t _CSB_T_REF, _CSB_BETA_T;
  int_t  _CSB_no_hydro;
  bool   _CSB_enabled;
};


class ConstForce : public VirtualForce
// ---------------------------------------------------------------------------
// A force constant in both space in time.
// ---------------------------------------------------------------------------
{
public:
  ConstForce    (Domain*, FEML*);
  void add      (AuxField*, const int_t, vector<AuxField*>&);
  void subtract (AuxField*, const int_t, vector<AuxField*>&);
private:
  real_t _v[3];	// Force components
};


class SteadyForce : public VirtualForce
// ---------------------------------------------------------------------------
// Constant in time but a function of space.
// ---------------------------------------------------------------------------
{
public:
  SteadyForce   (Domain*, FEML*);
  void add      (AuxField*, const int_t, vector<AuxField*>&);
  void subtract (AuxField*, const int_t, vector<AuxField*>&);
};


class WhiteNoiseForce : virtual public VirtualForce
// ---------------------------------------------------------------------------
// Forcing stochastic in space and time.
// ---------------------------------------------------------------------------
{
public:
  WhiteNoiseForce (Domain*, FEML*);
  void add        (AuxField*, const int_t, vector<AuxField*>&);
private:
  real_t _eps[3];
  int_t  _mode;
  int_t  _apply_step; // apply force every _apply_step'th step.
};


class ModulatedForce : virtual public VirtualForce
// ---------------------------------------------------------------------------
// Forcing which is a constant function of space (may be read from
// file) modulated by a function of time.
// ---------------------------------------------------------------------------
{
public:
  ModulatedForce (Domain*, FEML*);
  void add       (AuxField*, const int_t, vector<AuxField*>&);
private:
  char _alpha[3][StrMax]; // -- temporally varying part.
};


class SpatioTemporalForce : virtual public VirtualForce
// ---------------------------------------------------------------------------
// Forcing that is an arbitrary function of space-time.
// ---------------------------------------------------------------------------
{
public:
  SpatioTemporalForce (Domain*, FEML*);
  void add            (AuxField*, const int_t, vector<AuxField*>&);
private:
  char _alpha[3][StrMax]; // -- spatio-temporal-varying part.
};


class SpongeForce : virtual public VirtualForce
// ---------------------------------------------------------------------------
// Forcing penalises difference between velocity field and a given
// function of spatial position.
// ---------------------------------------------------------------------------
{
public:
  SpongeForce (Domain*, FEML*);
  void add    (AuxField*, const int_t, vector<AuxField*>&);
private:
  vector<AuxField*> _Uref;
  AuxField*         _mask;
  char              _mask_func[StrMax]; // -- mask function, f(x,y,z,t)
  int               _update; // mask update frequency
};


class DragForce : virtual public VirtualForce
// ---------------------------------------------------------------------------
// Forcing acts against velocity field according to its magnitude.
// ---------------------------------------------------------------------------
{
public:
  DragForce (Domain*, FEML*);
  void add  (AuxField*, const int_t, vector<AuxField*>&);
private:
  AuxField *_mask;
  AuxField *_umag;
};


class CoriolisForce : virtual public VirtualForce
// ---------------------------------------------------------------------------
// Forcing appropriate to solution in a rotating frame of reference.
// ---------------------------------------------------------------------------
{
public:
  CoriolisForce (Domain*, FEML*);
  void add      (AuxField*, const int_t, vector<AuxField*>&);
private:
  char           _omega[3][StrMax];    // -- angular velocity = f(t) ..
  char           _DomegaDt[3][StrMax]; // -- and its time derivative
  vector<real_t> _o;		       // -- evaluated at current time
  vector<real_t> _minus_o;	       // -- - omega
  vector<real_t> _minus_2o;	       // -- - 2 * omega
  vector<real_t> _DoDt;		       // -- evaluated at current time
  int_t          _unsteady;            // -- 1 if omega is unsteady
};


class SFDForce : virtual public VirtualForce
// ---------------------------------------------------------------------------
// Selective frequency damping: add forcing designed to achieve steady
// solution to NSE.  Physical space.  Parameters SFD_DELTA and SFD_CHI.
//
// Reference: Akervik et al., (2006), Steady solutions to the
// Navier--Stokes equations by selective frequency damping, Phys
// Fluids 18: 068102.
// ---------------------------------------------------------------------------
{
public:
  SFDForce (Domain*, FEML*);
  void add (AuxField*, const int_t, vector<AuxField*>&);
private:
  real_t _SFD_DELTA, _SFD_CHI;
};


class BuoyancyForce : virtual public VirtualForce
// ---------------------------------------------------------------------------
// Thermally driven buoyancy force as derived from the passive scalar field.
// ---------------------------------------------------------------------------
{
public:
  BuoyancyForce (Domain*, FEML*);
  void add      (AuxField*, const int_t, vector<AuxField*>&);
private:
  real_t _TREF, _BETAT, _gx, _gy, _omega;
  int_t  _centrifugal, _kineticgrad;
};

#endif
