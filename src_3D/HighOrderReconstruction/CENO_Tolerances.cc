/*!\file CENO_Tolerances.cc
  \brief Source file defining the values of the numerical tolerances declared in CENO_Tolerances.h */

/* Include required C++ libraries. */
// None

/* Using std namespace functions */
// None

/* Include CFFC header files */
#include "../MPI/MPI.h"
#include "CENO_Tolerances.h"

// CENO_Tolerances class

/*!
 * This value is used to differentiate between smooth and non-smooth solution reconstructions.
 * It is used in the post-reconstruction analysis.
 * Cells with values of smoothness indicator lower than Fit_Tolerance will have the 
 * reconstruction switched to a piecewise limited linear one.
 */
double CENO_Tolerances::Fit_Tolerance = 4000;

/*!
 * This value is used to decrease the level of switching in the proximity of the Fit_Tolerance.
 * It is used in the post-reconstruction analysis.
 * Cells with values of smoothness indicator lower than Fit_Tolerance_Buffer which were previously 
 * detected as non-smooth will have the reconstruction maintained to a piecewise limited linear one.
 */
double CENO_Tolerances::Fit_Tolerance_NonSensitivity = 5;

double CENO_Tolerances::Fit_Tolerance_Buffer = CENO_Tolerances::Fit_Tolerance_NonSensitivity * CENO_Tolerances::Fit_Tolerance;

/*!
 * This value is used in the computation of AMR criteria based on smoothness indicator
 */
double CENO_Tolerances::AMR_Smoothness_Units = 1.0;


/*!
 * This value is used in the post-reconstruction analysis.
 */
double CENO_Tolerances::epsilon = 1.0e-8;

/*!
 * This value is used to determine the maximum variation around a reference value 
 * for which the smoothness indicator is not computed.
 */
double CENO_Tolerances::epsilon_relative = 5.0e-3;
double CENO_Tolerances::epsilon_relative_square = CENO_Tolerances::epsilon_relative * CENO_Tolerances::epsilon_relative;

/*!
 * This value is used to determine the maximum variation around a reference value 
 * for which the smoothness indicator is not computed.
 */
double CENO_Tolerances::epsilon_absolute = 5.0e-7;
double CENO_Tolerances::epsilon_absolute_square = CENO_Tolerances::epsilon_absolute * CENO_Tolerances::epsilon_absolute;

double CENO_Tolerances::cross_epsilon = 2.0 * CENO_Tolerances::epsilon_absolute * CENO_Tolerances::epsilon_relative;

// Set default epsilon values
double CENO_Tolerances::epsilon_default = CENO_Tolerances::epsilon;
double CENO_Tolerances::epsilon_absolute_default = CENO_Tolerances::epsilon_absolute;
double CENO_Tolerances::epsilon_relative_default = CENO_Tolerances::epsilon_relative;
double CENO_Tolerances::Fit_Tolerance_default = CENO_Tolerances::Fit_Tolerance;
double CENO_Tolerances::Fit_Tolerance_NonSensitivity_default = CENO_Tolerances::Fit_Tolerance_NonSensitivity;
double CENO_Tolerances::AMR_Smoothness_Units_default = CENO_Tolerances::AMR_Smoothness_Units;

/*! Print the current execution mode
 *  at the output stream
 * \param [in] out_file the output stream
 */
void CENO_Tolerances::Print_Info(std::ostream & out_file){

  // output Fit_Tolerance
  out_file << "\n     -> Fit Tolerance: " << Fit_Tolerance
	   << "  (default value = " << Fit_Tolerance_default << ")";

  // output AMR_Smoothness_Units
  if (AMR_Smoothness_Units != AMR_Smoothness_Units_default){
    out_file << "\n     -> AMR Smoothness Units: " << AMR_Smoothness_Units
	     << "  (default value = " << AMR_Smoothness_Units_default << ")";
  }

  // output epsilon
  if (epsilon != epsilon_default){
    out_file << "\n     -> Epsilon: " << epsilon 
	     << "  (default value = " << epsilon_default << ")";
  }

  // output absolute epsilon
  if (epsilon_absolute != epsilon_absolute_default){
    out_file << "\n     -> Absolute allowed DeltaU: " << epsilon_absolute
	     << "  (default value = " << epsilon_absolute_default << ")";
  }

  // output relative epsilon
  if (epsilon_relative != epsilon_relative_default){
    out_file << "\n     -> Relative allowed DeltaU: " << epsilon_relative
	     << "  (default value = " << epsilon_relative_default << ")";
  }
}

/*!
 * Set default tolerance values
 */
void CENO_Tolerances::SetDefaults(void){
  epsilon = epsilon_default;
  epsilon_absolute = epsilon_absolute_default;
  epsilon_relative = epsilon_relative_default;
  Fit_Tolerance = Fit_Tolerance_default;
  Fit_Tolerance_NonSensitivity = Fit_Tolerance_NonSensitivity_default;
  AMR_Smoothness_Units = AMR_Smoothness_Units_default;

  UpdateDependentTolerances();
}

/*!
 * Broadcast the CENO_Tolerances variables to all      
 * processors associated with the specified communicator
 * from the specified processor using the MPI broadcast 
 * routine.
 *
 */
void CENO_Tolerances::Broadcast(void){
#ifdef _MPI_VERSION
  
  MPI::COMM_WORLD.Bcast(&epsilon,
			1, 
			MPI::DOUBLE, 0);
  MPI::COMM_WORLD.Bcast(&epsilon_relative,
			1, 
			MPI::DOUBLE, 0);
  MPI::COMM_WORLD.Bcast(&epsilon_absolute,
			1, 
			MPI::DOUBLE, 0);
  MPI::COMM_WORLD.Bcast(&Fit_Tolerance,
			1, 
			MPI::DOUBLE, 0);
  MPI::COMM_WORLD.Bcast(&Fit_Tolerance_NonSensitivity,
			1, 
			MPI::DOUBLE, 0);
  MPI::COMM_WORLD.Bcast(&AMR_Smoothness_Units,
			1, 
			MPI::DOUBLE, 0);

  // Update all dependent tolerances
  UpdateDependentTolerances();
#endif
}
