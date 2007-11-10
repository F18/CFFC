/*!\file SourceTermFields.cc
  \brief Source file initializing/implementing member variables/functions of class SourceTermFields. */

/* Include required C++ libraries. */
// None

/* Using std namespace functions */
// None

/* Include CFFC header files */
#include "SourceTermFields.h"

// ==== Member functions ====
SourceFieldBasicType * SourceTermFields::Source = NULL;
SourceTermFields::SourceTermFields(void){ }

SourceTermFields& SourceTermFields::getInstance(void){
  static SourceTermFields inst;
  std::atexit(Destroy);
  return inst;
}

