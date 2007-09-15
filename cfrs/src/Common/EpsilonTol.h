/* EpsilonTol.h: Header file defining some numerical tolerances used 
                 in calculations.                                    */

#ifndef _EPSILONTOL_INCLUDED
#define _EPSILONTOL_INCLUDED

/* Include header files */
#include <cmath>

class CENO_EpsilonTol{

 public:
  static double epsilon;
  static double epsilon_relative;
  static double epsilon_absolute;
  static double epsilon_relative_square;
  static double epsilon_absolute_square;
  static double cross_epsilon; /*  2*epsilon_relative*epsilon_absolute  */

  /* These functions can be used to determine what 
     an acceptable tolerance is around the quantity U. */
  static double ToleranceAroundValue(const double & U);
  static double SquareToleranceAroundValue(const double & U);

 protected:
  CENO_EpsilonTol(){};

};

/* ToleranceAroundValue:
   Returns the absolute DeltaU around the quantity U,
   based on relative and absolute tolerance.*/
inline double CENO_EpsilonTol::ToleranceAroundValue(const double & U){
  return epsilon_absolute + epsilon_relative*fabs(U);
}

/* SquareToleranceAroundValue:
   Returns the square of DeltaU around the quantity U,
   based on relative and absolute tolerance.*/
inline double CENO_EpsilonTol::SquareToleranceAroundValue(const double & U){
  return epsilon_absolute_square + cross_epsilon*fabs(U) + epsilon_relative_square*U*U;
}

#endif
