#ifndef BSYS_H
#define BSYS_H

class BoundarySys
// ===========================================================================
// This class automates the retrieval of boundary condition
// applicators (Boundary objects) for a given Field and Fourier
// mode
//
#if 0
// global numbering schemes (NumberSys) and inverse mass matrix
#endif
//
// ===========================================================================
{
public:
  BoundarySys  (BCmgr*, const vector<Element*>&, const char);
  ~BoundarySys () { };

  char                     field () const { return _field_name; }
  int_t                    nSurf () const { return _nbound; }
  bool                     mixBC () const { return _mixed; }
  const vector<Boundary*>& BCs   (const int_t) const;
#if 0  
  const NumberSys*         Nsys  (const int_t) const;
  const real_t*            Imass (const int_t) const;
#endif

private:
  char               _field_name;
  int_t              _nbound    ;  // Number of element edges with BCs.
  bool               _mixed     ;  // Flags presence of mixed BC type.
  vector<Boundary*>* _boundary  ;  // Boundary*'s  for modes 0, 1, 2.
#if 0  
  NumberSys**        _number    ;  // NumberSys*'s for modes 0, 1, 2.
#endif  
};

#endif
