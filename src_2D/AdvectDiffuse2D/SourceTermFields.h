/*!\file SourceTermFields.h
  \brief Header file defining 2D source term fields. */

#ifndef _SOURCETERM_FIELDS_INCLUDED
#define _SOURCETERM_FIELDS_INCLUDED

/* Include required C++ libraries. */
// None

/* Using std namespace functions */
// None

/* Include CFFC header files */
#include "../Utilities/Utilities.h"
#include "SourceFields.h"

class SourceTermFields{
public:
  static SourceTermFields& getInstance(void);
  
  //  double operator()(double x, double y){return Source->operator()(x,y);}
  static SourceFieldBasicType *Source;

  void SetSource(int field){
    switch (field){
    case 1:
      Source = new ZERO_SourceField;
      break;

    case 2:
      Source = new Linear_SourceField;
      break;

    case 3:
      Source = new Exponential_SourceField;
      break;
    }
  }

  static void Destroy(void){ delete Source;}

protected:
  SourceTermFields();		//!< Private constructor
  SourceTermFields(const SourceTermFields&); //!< Private copy constructor
  SourceTermFields& operator=(const SourceTermFields&); //!< Private assignment operator
  
private:
  //  SourceFieldBasicType *Source;
};

#endif
