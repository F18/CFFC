/*!\file CENO_ExecutionMode.h
  \brief Definition of flags that control the execution of CENO high-order reconstruction. */

#ifndef _CENO_EXECUTIONMODE_INCLUDED
#define _CENO_EXECUTIONMODE_INCLUDED

/* Include required C++ libraries. */
// None

/* Using std namespace functions */
// None

/* Include CFFC header files */
#include "../Utilities/Utilities.h"

/*!
 * Class: CENO_Execution_Mode
 *
 * @brief Definition of flags that control the execution of CENO high-order reconstruction
 *
 */
class CENO_Execution_Mode{
  
public:

  // set all flags to default values
  static void SetDefaults(void);
  
  /* Store the pseudo-inverse of the LHS term in the CENO reconstruction for every computational cell.
     Turn ON if you want to run in speed efficient mode. However, the memory requirements will
     increase considerably!!!
     Turn OFF for running in the memory efficient mode.
     ---------------------------------------------------------------------------------------- */
  static short CENO_SPEED_EFFICIENT;


  /* This flag controls whether the order is dropped or not in regions detected as non-smooth.
     Turn ON if you want the CENO scheme to produce monotone solutions in non-smooth regions.
     Turn OFF if you want to use high-order interpolant throughout the whole computational domain.
     ATTENTION: If the flag is OFF the solution might blow up !!! 
     ----------------------------------------------------------------------------------------  */
  static short CENO_DROP_ORDER;  


  /* Use CENO smoothness indicator with padding (set to ON to have Piecewise Linear reconstructions
     in all the cells surrounding a cell with non-smooth reconstruction)
     ---------------------------------------------------------------------------------------- */
  static short CENO_Padding;

  /* Turn ON this flag if the geometric weighting used in the least-squares problem is 1.0/(Distance^2)
     Turn OFF if the geometric weighting is 1.0/fabs(Distance)
     ---------------------------------------------------------------------------------------- */
  static short CENO_SQUARE_GEOM_WEIGHTING;


  /* Turn ON if you want to include the geometric weights used in the reconstruction process
     as part of the computation of the smoothness indicator.
     This works only if the CENO_SUPER_FAST is turned ON (because the geometric weights are saved).
     Turn OFF for treating all the terms in the residual sum (SS_Residual) and regression sum (SS_Regression)
     equally (i.e. each of them has a geoemtric weight equals to 1 in the computation of the 
     smoothness indicator).
     ---------------------------------------------------------------------------------------- */
  static short CENO_CONSIDER_WEIGHTS;

};

#endif
