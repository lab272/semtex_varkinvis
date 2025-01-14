//////////////////////////////////////////////////////////////////////////////
// family.cpp
//
// Copyright (c) 2003+, Hugh M Blackburn
//
// This file maintains a family structure for vectors of real. Doing
// it in C++ ensures that we can delete copies of arrays allocated
// using new. See the equivalent family routines in Femlib.
//
// This would be templated if instantiation was more standard
// across compilation regimes (or I could figure out how to make it
// all work), but the minimum requirement is real arrays...so that's
// all for now.
//
// To do: reimplement this with smart pointers.
//
//////////////////////////////////////////////////////////////////////////////

#include <sem.h>

class rvect { public: int_t size; real_t* data; int_t nrep; };

static vector<rvect*> rv;

namespace Family {
static real_t* adopted (const int_t size, const real_t* src)
{
  vector<rvect*>::iterator p;
  bool found;
  for (found = false, p = rv.begin(); p != rv.end(); p++) {
    if ((*p) -> size != size) continue;
    if (found = src == (*p) -> data) break;
    if (found = Veclib::same (size, src, 1, (*p) -> data, 1))
      { ++(*p) -> nrep; break; }
  }
  return found ? (*p) -> data : 0;
}

void abandon (real_t** vect)
{
  vector<rvect*>::iterator p;
  for (p = rv.begin(); p != rv.end(); p++)
    if ((*p) -> data == *vect) {
      if (--(*p) -> nrep == 0) { delete[] (*p) -> data; rv.erase (p); }
      break;
    }
}

void adopt (const int_t size, real_t** vect)
{
  if (!vect || !*vect) return;

  rvect*  S = 0;
  real_t* member;

  if ((member = adopted (size, *vect)) && member != *vect) {
    delete[] *vect; *vect = member;
  } else {
    S = new rvect; S -> size = size; S -> data = *vect; S -> nrep = 1;
    rv.push_back (S);
  }
}
}
