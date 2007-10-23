/*!\file CENO_Tolerances.cc
  \brief Source file defining the values of the numerical tolerances declared in CENO_Tolerances.h */

#include "CENO_Tolerances.h"

// CENO_EpsilonTol class

/*!
 * These flags are 'true' if the correspondent default tolerances have been changed, otherwise are 'false'. 
 */
bool CENO_EpsilonTol::ChangedDefault_Epsilon = false;
bool CENO_EpsilonTol::ChangedDefault_EpsilonAbsolute = false;
bool CENO_EpsilonTol::ChangedDefault_EpsilonRelative = false;

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

// Set default epsilon values
double CENO_EpsilonTol::epsilon_default = CENO_EpsilonTol::epsilon;
double CENO_EpsilonTol::epsilon_absolute_default = CENO_EpsilonTol::epsilon_absolute;
double CENO_EpsilonTol::epsilon_relative_default = CENO_EpsilonTol::epsilon_relative;

/*! Print the current execution mode
 *  at the output stream
 * \param [in] out_file the output stream
 */
void CENO_EpsilonTol::Print_Info(std::ostream & out_file){

  // output epsilon
  if (ChangedDefault_Epsilon){
    out_file << "\n     -> Epsilon = " << epsilon 
	     << "  (default value = " << epsilon_default << ")";
  }

  // output absolute epsilon
  if (ChangedDefault_EpsilonAbsolute){
    out_file << "\n     -> Absolute allowed DeltaU = " << epsilon_absolute
	     << "  (default value = " << epsilon_absolute_default << ")";
  }

  // output relative epsilon
  if (ChangedDefault_EpsilonRelative){
    out_file << "\n     -> Relative allowed DeltaU = " << epsilon_relative
	     << "  (default value = " << epsilon_relative_default << ")";
  }
}

/*!
 * Set default tolerance values
 */
void CENO_EpsilonTol::SetDefaults(void){
  epsilon = epsilon_default;
  epsilon_absolute = epsilon_absolute_default;
  epsilon_relative = epsilon_relative_default;

  UpdateDependentTolerances();
}
