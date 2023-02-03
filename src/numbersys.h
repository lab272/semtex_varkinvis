#ifndef NUMBERSYS_H
#define NUMBERSYS_H

class NumberSys
// ===========================================================================
// This class automates the retrieval of assembly map applicators for
// a given Field.  It is analogous to BoundarySys but deals with
// AssemblyMaps.
//
// ===========================================================================
{
public:
  NumberSys  (const vector<AssemblyMap*>&, const char);
  ~NumberSys () { };

  char                field  () const { return _field_name; }
  const AssemblyMap*  getMap (const int_t) const;

private:
  char                 _field_name;
  vector<AssemblyMap*> _map       ;  // AssemblyMaps's for modes 0, 1, 2.
};

#endif
