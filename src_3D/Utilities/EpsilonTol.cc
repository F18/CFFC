/*!\file EpsilonTol.cc
  \brief Source file defining the values of the numerical tolerances declared in EpsilonTol.h */

#include "EpsilonTol.h"

// EpsilonTol class
double EpsilonTol::epsilon = 1.0e-14;
double EpsilonTol::epsilon_relative = 1.0e-14;
double EpsilonTol::epsilon_relative_square = EpsilonTol::epsilon_relative * EpsilonTol::epsilon_relative;
double EpsilonTol::epsilon_absolute = 1.0e-14;
double EpsilonTol::epsilon_absolute_square = EpsilonTol::epsilon_absolute * EpsilonTol::epsilon_absolute;
double EpsilonTol::cross_epsilon = 2.0 * EpsilonTol::epsilon_absolute * EpsilonTol::epsilon_relative;
double EpsilonTol::MachineEps = 1.0e-14;

// CENO_EpsilonTol class
double CENO_EpsilonTol::epsilon = 1.0e-8;
double CENO_EpsilonTol::epsilon_relative = 5.0e-3;
double CENO_EpsilonTol::epsilon_relative_square = CENO_EpsilonTol::epsilon_relative * CENO_EpsilonTol::epsilon_relative;
double CENO_EpsilonTol::epsilon_absolute = 5.0e-7;
double CENO_EpsilonTol::epsilon_absolute_square = CENO_EpsilonTol::epsilon_absolute * CENO_EpsilonTol::epsilon_absolute;
double CENO_EpsilonTol::cross_epsilon = 2.0 * CENO_EpsilonTol::epsilon_absolute * CENO_EpsilonTol::epsilon_relative;
