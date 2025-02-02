#ifndef BOUNDARY_H
#define BOUNDARY_H

class Boundary : public Edge
// ===========================================================================
// Scalar field element-wise boundary class.  These evaluate and apply
// a (boundary) Condition to a Field on a particular side (Edge) of an
// Element that lies on the periphery of a computational Domain.
// ===========================================================================
{
public:
  Boundary (const int_t id, const char* group, const Condition* bcond,
	    const Element* elmt, const int_t side):
              Edge (group, elmt, side), _id (id), _bcond (bcond) { }

  int_t ID        () const { return _id; }
  void  print     () const;

  void  evaluate  (const Field*,const int_t,const int_t,const bool,
		   real_t*)                                              const;

  // -- Impose essential BCs:
  void  set       (const real_t*,const int_t*,real_t*)                   const;

  // -- Apply natural BCs:
  void  sum       (const real_t*,const int_t*,real_t*,real_t*)           const;

  // -- Apply mixed BCs:
  void  augmentSC (const int_t,const int_t,const int_t*,real_t*,real_t*) const;
  void  augmentOp (const int_t*,const real_t*,real_t*)                   const;
  void  augmentDg (const int_t*,real_t*)                                 const;

private:
  int_t            _id   ;	// Ident number.
  const Condition* _bcond;	// Boundary condition.
};

#endif
