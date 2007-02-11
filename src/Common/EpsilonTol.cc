/* EpsilonTol.cc: Source file defining the values of numerical tolerances  
                  declared in EpsilonTol.h                               */

#include "EpsilonTol.h"

double CENO_EpsilonTol::epsilon = 1.0e-6;
double CENO_EpsilonTol::epsilon_relative = 1.0e-6;
double CENO_EpsilonTol::epsilon_relative_square = 1.0e-12;
double CENO_EpsilonTol::epsilon_absolute = 1.0e-8;
double CENO_EpsilonTol::epsilon_absolute_square = 1.0e-16;
double CENO_EpsilonTol::cross_epsilon = 2.0e-14;
