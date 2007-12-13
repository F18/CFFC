/*!\file ENO_ExecutionMode.h
  \brief Definition of flags that control the execution of ENO high-order reconstruction. */

#ifndef _ENO_EXECUTIONMODE_INCLUDED
#define _ENO_EXECUTIONMODE_INCLUDED

/* Include required C++ libraries. */
// None

/* Using std namespace functions */
// None

/* Include CFFC header files */
#include "../Utilities/Utilities.h"

/*!
 * Class: ENO_Execution_Mode
 *
 * @brief Definition of flags that control the execution of ENO high-order reconstruction
 *
 */
class ENO_Execution_Mode{
  
public:

  // set all flags to default values
  static void SetDefaults(void);

  /* This flag controls whether the ENO algorithm is used in the computation.
     Turn ON if you want to use the high-order ENO scheme.
     Turn OFF if you don't want to use this scheme.
     ATTENTION: Some settings/memory allocations are done only if this flag is ON.
     ----------------------------------------------------------------------------------------  */
  static short USE_ENO_ALGORITHM;

  /* Use the piecewise constant state(s) to solve the Riemann problem at the interface of two cells 
     if any/both of the interface states calculated with the reconstruction interpolant
     has/have non-physical values (e.g. negative pressure or density).
     Turn ON if this feature is desired (default)
     Turn OFF if you don't want to use this feature.
     --------------------------------------------------------------------------------------- */
  static short FORCE_WITH_PIECEWISE_CONSTANT_AT_INTERFACE;

  /* Get the primitive variables directly from the characteristic ones when ENO reconstruction 
     with characteristics variables is used. 
     Turn ON if this feature is desired.
     Turn OFF if it's desired that the primitive variables are calculated based on the conserved
     ones, which are in turn determined from the characteristic ones. (default)
     ---------------------------------------------------------------------------------------- */
  static short USE_PRIMITIVE_VARIABLES;

private:
  ENO_Execution_Mode(void){};

};

#endif
