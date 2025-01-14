#ifndef AUXFIELD_H
#define AUXFIELD_H


class AuxField
// ===========================================================================
// Physical field storage class, no BCs.
//
// An AuxField has a storage area for data, physical space mesh, and a
// matching list of elements which can access this storage area on an
// element-wise basis if required.  Functions to operate on fields as
// a whole are included.  AuxFields have no boundary conditions, no
// global numbering system, and no solution routines.
// ===========================================================================
{
friend istream& operator >> (istream&, AuxField&);
friend ostream& operator << (ostream&, AuxField&);
friend class    Field;
friend class    BCmgr;

public:
  AuxField (real_t*, const int_t, vector<Element*>&, const char = 0);

  char name     ()      const { return _name; }
  void describe (char*) const;

  AuxField& operator  = (const real_t);
  AuxField& operator += (const real_t);
  AuxField& operator -= (const real_t);
  AuxField& operator *= (const real_t);
  AuxField& operator /= (const real_t);

  AuxField& operator  = (const AuxField&);
  AuxField& operator += (const AuxField&);
  AuxField& operator -= (const AuxField&);
  AuxField& operator *= (const AuxField&);

  AuxField& operator  = (const char*);
  AuxField& axpy        (const real_t, const AuxField&);

  AuxField& times       (const AuxField&, const AuxField&);
  AuxField& timesPlus   (const AuxField&, const AuxField&);
  AuxField& timesMinus  (const AuxField&, const AuxField&);

  AuxField& addToPlane  (const int_t, const real_t);
  AuxField& getPlane    (const int_t, real_t*);
  AuxField& setPlane    (const int_t, const real_t*);
  AuxField& setPlane    (const int_t, const real_t);

  AuxField& update   (const int_t, const real_t*, const real_t, const real_t); 

  AuxField& gradient (const int_t);
  AuxField& mulY     ();
  AuxField& divY     ();

  real_t mode_L2     (const int_t)         const;
  real_t CFL         (const int_t, int_t&) const;
  real_t integral    ()                    const;
  real_t area        ()                    const;
  
  real_t probe       (const Element*, const real_t, 
		      const real_t, const int_t)    const;
  real_t probe       (const Element*, const real_t,
		      const real_t, const real_t)   const;

  AuxField&   reverse   ();
  AuxField&   smooth    (const int_t, const int_t*, const real_t*);
  
  static void swapData  (AuxField*, AuxField*);
  static void couple    (AuxField*, AuxField*, const int_t);
  
  real_t* getData() const;

protected:
  char              _name ;	// Identification tag.  '\0' by default.
  vector<Element*>& _elmt ;	// Quadrilateral elements.
  int_t             _nz   ;	// number of data planes (per process).
  int_t             _size ;	// _nz * Geometry::planeSize().
  real_t*           _data ;	// 2/3D data area, element x element x plane.
  real_t**          _plane;	// Pointer into data for each 2D frame.

private:
  AuxField& operator /= (const AuxField&) { return *this; }

};


// -- Related routins for I/O that don't use any private data.

void readField  (istream&, vector<AuxField*>&);
void writeField (ostream&, const char*, const int_t, const real_t,
		 vector<AuxField*>&);

#endif
