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

class EpsilonTol{

 public:
  static double epsilon;
  static double epsilon_relative;
  static double epsilon_absolute;
  static double epsilon_relative_square;
  static double epsilon_absolute_square;
  static double cross_epsilon; /*  2*epsilon_relative*epsilon_absolute  */

  static double MachineEps;

  /* These functions can be used to determine what 
     an acceptable tolerance is around the quantity U. */
  static double ToleranceAroundValue(const double & U);
  static double SquareToleranceAroundValue(const double & U);

  /* Output operator. */
  static void Print_Tolerances(ostream& os = std::cout);

 protected:
  EpsilonTol(){};

};

/* ToleranceAroundValue:
   Returns the absolute DeltaU around the quantity U,
   based on relative and absolute tolerance.*/
inline double EpsilonTol::ToleranceAroundValue(const double & U){
  return epsilon_absolute + epsilon_relative*fabs(U);
}

/* SquareToleranceAroundValue:
   Returns the square of DeltaU around the quantity U,
   based on relative and absolute tolerance.*/
inline double EpsilonTol::SquareToleranceAroundValue(const double & U){
  return epsilon_absolute_square + cross_epsilon*fabs(U) + epsilon_relative_square*U*U;
}

/* Print_Tolerances:
   Outputs the tolerances to the standard output. */
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


/* Tolerance class for CENO reconstruction */
class CENO_EpsilonTol: public EpsilonTol{

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

  /* Output operator. */
  static void Print_CENO_Tolerances(ostream& os = std::cout);

 protected:
  CENO_EpsilonTol(){};

};

/* Print_CENO_Tolerances:
   Outputs the CENO tolerances to the standard output. */
inline void CENO_EpsilonTol::Print_CENO_Tolerances(ostream& os){
  os << "\nCENO_Tolerances:\n"
     << "epsilon=" << epsilon << "\n"
     << "epsilon_relative=" << epsilon_relative << "\n"
     << "epsilon_relative_square=" << epsilon_relative_square << "\n"
     << "epsilon_absolute=" << epsilon_absolute << "\n"
     << "epsilon_absolute_square=" << epsilon_absolute_square << "\n"
     << "cross_epsilon=" << cross_epsilon << "\n"
     << "MachineEpsilon=" << MachineEps << "\n";
}

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
