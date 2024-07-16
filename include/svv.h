#ifndef SVV_H
#define SVV_H

#include <cfemdef.h>

// ===========================================================================
//
// ===========================================================================

namespace SVV {

  const real_t* coeffs    (const int_t);
  const real_t* coeffs_z  (const int_t);

  void          operators (const int_t, const real_t**, const real_t**);
}

#endif
