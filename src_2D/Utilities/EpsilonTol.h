/*!\file EpsilonTol.h
  \brief Header file defining numerical tolerances used in calculations.  */

#ifndef _EPSILONTOL_INCLUDED
#define _EPSILONTOL_INCLUDED

/* Include required C++ libraries. */
#include <cmath>
#include <iostream>

/* Using std namespace functions */
using std::ostream;

/* Include CFFC header files */
// None

/*!
 * \class EpsilonTol
 * 
 * \brief General usage tolerance class.
 *********************************************************/
class EpsilonTol{

 public:
  static double epsilon;	         //!< epsilon value (e.g. used to avoid division by zero)
  static double epsilon_relative;        //!< epsilon applied relative to a reference value
  static double epsilon_absolute;        //!< absolute epsilon value
  static double epsilon_relative_square; //!< square of the relative epsilon
  static double epsilon_absolute_square; //!< square of the absolute epsilon
  static double cross_epsilon;           //!< equal to 2*epsilon_relative*epsilon_absolute

  static double MachineEps;	         //!< machine epsilon

  /* These functions can be used to determine what 
     an acceptable tolerance is around the quantity U. */
  static double ToleranceAroundValue(const double & U);
  static double SquareToleranceAroundValue(const double & U);
  static double getAccuracyBasedOnExactDigits(const int &digits){ return 0.5*pow(10.0,-digits); }

  /* Output operator. */
  static void Print_Tolerances(ostream& os = std::cout);

 protected:
  EpsilonTol(){};

};

/*! 
 * Returns the absolute DeltaU around the quantity U,
 * based on relative and absolute tolerance.
 */
inline double EpsilonTol::ToleranceAroundValue(const double & U){
  return epsilon_absolute + epsilon_relative*fabs(U);
}

/*! 
 * Returns the square of DeltaU around the quantity U,
 * based on relative and absolute tolerance.
 */
inline double EpsilonTol::SquareToleranceAroundValue(const double & U){
  return epsilon_absolute_square + cross_epsilon*fabs(U) + epsilon_relative_square*U*U;
}

/*! 
 * Output the tolerances to the standard output.
 */
inline void EpsilonTol::Print_Tolerances(ostream& os){
  os << "\nGeneral Tolerances:\n"
     << "epsilon=" << epsilon << "\n"
     << "epsilon_relative=" << epsilon_relative << "\n"
     << "epsilon_relative_square=" << epsilon_relative_square << "\n"
     << "epsilon_absolute=" << epsilon_absolute << "\n"
     << "epsilon_absolute_square=" << epsilon_absolute_square << "\n"
     << "cross_epsilon=" << cross_epsilon << "\n"
     << "MachineEpsilon=" << MachineEps << "\n";
}

#endif
