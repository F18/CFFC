/////////////////////////////////////////////////////////////////////
///
/// \file Flame2DState.cc
/// 
/// \author Marc R.J. Charest
/// 
/// \brief This header file contains the class definition for 
///        describing the state of the chemically reacting, 
///        multispecies ideal gas mixture.
///
/////////////////////////////////////////////////////////////////////
#include "Flame2DState.h"

/**
 * Initialization of static variables.
 */
#ifndef STATIC_NUMBER_OF_SPECIES
int Flame2D_State :: n = 0;
int Flame2D_State :: ns = 0;
#endif
double Flame2D_State::Mref=0.5;
double Flame2D_State::gravity_z=-9.81;


