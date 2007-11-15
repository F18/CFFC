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
