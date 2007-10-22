/*!\file CENO_Tolerances.h
  \brief Header file defining numerical tolerances used in calculations with CENO reconstruction algorithm. */

#ifndef _CENO_TOLERANCES_INCLUDED
#define _CENO_TOLERANCES_INCLUDED

/* Include required C++ libraries. */
// None

/* Using std namespace functions */
// None

/* Include CFFC header files */
#include "../Utilities/EpsilonTol.h" // general usage tolerances

/*!
 * \class CENO_EpsilonTol
 * 
 * \brief Tolerance class for CENO reconstruction algorithm.
 *********************************************************/
class CENO_EpsilonTol: public EpsilonTol{

 public:
  static double epsilon;                  //!< limit the maximum value taken by the smoothness indicator 
  static double epsilon_relative;         //!< used in computing tolerance around value as relative variation
  static double epsilon_absolute;         //!< used in computing tolerance around value as absolute variation
  static double epsilon_relative_square;  //!< the square of the relative epsilon
  static double epsilon_absolute_square;  //!< the square of the absolute epsilon
  static double cross_epsilon;            //!< equal to 2*epsilon_relative*epsilon_absolute

  /* These functions can be used to determine what 
     an acceptable tolerance is around the quantity U. */
  static double ToleranceAroundValue(const double & U);
  static double SquareToleranceAroundValue(const double & U);

  /* Output operator. */
  static void Print_CENO_Tolerances(ostream& os = std::cout);

 protected:
  CENO_EpsilonTol(){};

};

/*!
 * Output the CENO tolerances to the standard output.
 */
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

/*! 
 * Returns the allowed variation DeltaU around the quantity U,
 * based on relative and absolute tolerance, for which 
 * the solution can be considered smooth, without computing
 * the smoothness indicator.
 */
inline double CENO_EpsilonTol::ToleranceAroundValue(const double & U){
  return epsilon_absolute + epsilon_relative*fabs(U);
}

/*! 
 * Returns the square of DeltaU around the quantity U,
 * based on relative and absolute tolerance.
 */
inline double CENO_EpsilonTol::SquareToleranceAroundValue(const double & U){
  return epsilon_absolute_square + cross_epsilon*fabs(U) + epsilon_relative_square*U*U;
}

#endif
