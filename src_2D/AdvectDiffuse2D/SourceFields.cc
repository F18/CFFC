/*!\file SourceFields.cc
  \brief Source file initializing/implementing member variables/functions that belong to classes defined in SourceFields.h */

/* Include required C++ libraries. */
// None

/* Using std namespace functions */
// None

/* Include CFFC header files */
#include "SourceFields.h"

SourceFieldBasicType::~SourceFieldBasicType(void){cout << "Destructor SourceFieldBasicType\n";  }
SourceFieldThatRequireIntegration::~SourceFieldThatRequireIntegration(void){
  cout << "Destructor SourceFieldThatRequireIntegration\n"; }
SourceFieldThatDoesNotRequireIntegration::~SourceFieldThatDoesNotRequireIntegration(void){
  cout << "Destructor SourceFieldThatDoesNotRequireIntegration\n";
}
