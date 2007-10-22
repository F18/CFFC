/*!\file CENO_Tolerances.cc
  \brief Source file defining the values of the numerical tolerances declared in CENO_Tolerances.h */

#include "CENO_Tolerances.h"

// CENO_EpsilonTol class

/*!
 * This value is used in the post reconstruction analysis.
 */
double CENO_EpsilonTol::epsilon = 1.0e-8;

/*!
 * This value is used to determine the maximum variation around a reference value 
 * for which the smoothness indicator is not computed.
 */
double CENO_EpsilonTol::epsilon_relative = 5.0e-3;
double CENO_EpsilonTol::epsilon_relative_square = CENO_EpsilonTol::epsilon_relative * CENO_EpsilonTol::epsilon_relative;

/*!
 * This value is used to determine the maximum variation around a reference value 
 * for which the smoothness indicator is not computed.
 */
double CENO_EpsilonTol::epsilon_absolute = 5.0e-7;
double CENO_EpsilonTol::epsilon_absolute_square = CENO_EpsilonTol::epsilon_absolute * CENO_EpsilonTol::epsilon_absolute;

double CENO_EpsilonTol::cross_epsilon = 2.0 * CENO_EpsilonTol::epsilon_absolute * CENO_EpsilonTol::epsilon_relative;
