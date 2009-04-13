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
#include "../CFD/CFD.h"

/*!
 * \class CENO_Tolerances
 * 
 * \brief Tolerance class for CENO reconstruction algorithm.
 *********************************************************/
class CENO_Tolerances: public EpsilonTol{

 public:
  static double epsilon;                  //!< limit the maximum value taken by the smoothness indicator 
  static double epsilon_relative;         //!< used in computing tolerance around value as relative variation
  static double epsilon_absolute;         //!< used in computing tolerance around value as absolute variation
  static double epsilon_relative_square;  //!< the square of the relative epsilon
  static double epsilon_absolute_square;  //!< the square of the absolute epsilon
  static double cross_epsilon;            //!< equal to 2*epsilon_relative*epsilon_absolute
  static double Fit_Tolerance;	          //!< value used to distinguish between smooth and non-smooth solution reconstructions.
  static double AMR_Smoothness_Units;     //!< value used in computing the refinement criterion based on smoothness indicator
  static double Fit_Tolerance_NonSensitivity; /*!< value used to create a buffer zone for switching a non-smooth solution
					           to a smooth one. */
  static double Fit_Tolerance_Buffer;     /*!< set a threshold for switching a previously identified non-smooth interpolant
					       to a smooth one. */

  /* These functions can be used to determine what 
     an acceptable tolerance is around the quantity U. */
  static double ToleranceAroundValue(const double & U);
  static double SquareToleranceAroundValue(const double & U);

  /* Output operator. */
  static void Print_CENO_Tolerances(ostream& os = std::cout);

  /* Set tolerances at runtime. */
  template<class Input_Parameters_Type>
  static void Parse_Next_Input_Control_Parameter(Input_Parameters_Type & IP, int & i_command);

  static void Print_Info(std::ostream & out_file);

  static void SetDefaults(void);

  static void Broadcast(void);

  // Update tolerances that depend on the values of other tolerances
  static void UpdateDependentTolerances(void);

 protected:
  CENO_Tolerances();      //!< Private default constructor
  CENO_Tolerances(const CENO_Tolerances&); //!< Private copy constructor
  CENO_Tolerances& operator=(const CENO_Tolerances&); //!< Private assignment operator

  static double epsilon_default;          //!< this is a copy of epsilon that cannot be modified at runtime
  static double epsilon_relative_default; //!< this is a copy of epsilon_relative that cannot be modified at runtime
  static double epsilon_absolute_default; //!< this is a copy of epsilon_absolute that cannot be modified at runtime
  static double Fit_Tolerance_default;    //!< this is a copy of Fit_Tolerance that cannot be modified at runtime
  static double AMR_Smoothness_Units_default;  //!< this is a copy of AMR_Smoothness_Units that cannot be modified at runtime
  static double Fit_Tolerance_NonSensitivity_default;

};

/*!
 * Output the CENO tolerances to the standard output.
 */
inline void CENO_Tolerances::Print_CENO_Tolerances(ostream& os){
  os << "\nCENO_Tolerances:\n"
     << "epsilon=" << epsilon << "\n"
     << "epsilon_relative=" << epsilon_relative << "\n"
     << "epsilon_relative_square=" << epsilon_relative_square << "\n"
     << "epsilon_absolute=" << epsilon_absolute << "\n"
     << "epsilon_absolute_square=" << epsilon_absolute_square << "\n"
     << "cross_epsilon=" << cross_epsilon << "\n"
     << "MachineEpsilon=" << MachineEps << "\n"
     << "Fit_Tolerance=" << Fit_Tolerance << "\n"
     << "Fit_Tolerance_NonSensitivity=" << Fit_Tolerance_NonSensitivity << "\n"
     << "Fit_Tolerance_Buffer=" << Fit_Tolerance_Buffer << "\n"
     << "AMR_Smoothness_Units=" << AMR_Smoothness_Units << "\n";
}

/*! 
 * Returns the allowed variation DeltaU around the quantity U,
 * based on relative and absolute tolerance, for which 
 * the solution can be considered smooth, without computing
 * the smoothness indicator.
 */
inline double CENO_Tolerances::ToleranceAroundValue(const double & U){
  return epsilon_absolute + epsilon_relative*fabs(U);
}

/*! 
 * Returns the square of DeltaU around the quantity U,
 * based on relative and absolute tolerance.
 */
inline double CENO_Tolerances::SquareToleranceAroundValue(const double & U){
  return epsilon_absolute_square + cross_epsilon*fabs(U) + epsilon_relative_square*U*U;
}

/*! 
 * Update tolerances that depend on the values of other tolerances
 * (e.g. epsilon_absolute_square depends on epsilon_absolute)
 */
inline void CENO_Tolerances::UpdateDependentTolerances(void){
  epsilon_absolute_square = epsilon_absolute * epsilon_absolute;
  epsilon_relative_square = epsilon_relative * epsilon_relative;
  cross_epsilon = 2.0 * epsilon_absolute * epsilon_relative;
  Fit_Tolerance_Buffer = Fit_Tolerance * Fit_Tolerance_NonSensitivity;
}


/*!
 * Parse the input control parameters for 
 * settings related to CENO_Tolerances class
 */
template<class Input_Parameters_Type>
inline void CENO_Tolerances::Parse_Next_Input_Control_Parameter(Input_Parameters_Type & IP,
								int & i_command){

  // Check if the next control parameter has already been identified
  if (i_command != INVALID_INPUT_CODE){
    return;
  }

  char buffer[256];

  // Try to match the next control parameter
  if (strcmp(IP.Next_Control_Parameter, "CENO_Epsilon") == 0) {
    i_command = 0;
    ++IP.Line_Number;
    IP.Input_File >> epsilon;
    IP.Input_File.getline(buffer, sizeof(buffer));

  } else if (strcmp(IP.Next_Control_Parameter, "CENO_Absolute_Epsilon") == 0) {
    i_command = 0;
    ++IP.Line_Number;
    IP.Input_File >> epsilon_absolute;
    IP.Input_File.getline(buffer, sizeof(buffer));

  } else if (strcmp(IP.Next_Control_Parameter, "CENO_Relative_Epsilon") == 0) {
    i_command = 0;
    ++IP.Line_Number;
    IP.Input_File >> epsilon_relative;
    IP.Input_File.getline(buffer, sizeof(buffer));

  } else if (strcmp(IP.Next_Control_Parameter, "CENO_Tolerance") == 0) {
    i_command = 0;
    ++IP.Line_Number;
    IP.Input_File >> Fit_Tolerance;
    IP.Input_File.getline(buffer, sizeof(buffer));

  } else if (strcmp(IP.Next_Control_Parameter, "CENO_NonSensitivity") == 0) {
    i_command = 0;
    ++IP.Line_Number;
    IP.Input_File >> Fit_Tolerance_NonSensitivity;
    IP.Input_File.getline(buffer, sizeof(buffer));

  } else if (strcmp(IP.Next_Control_Parameter, "CENO_AMR_Units") == 0) {
    i_command = 0;
    ++IP.Line_Number;
    IP.Input_File >> AMR_Smoothness_Units;
    IP.Input_File.getline(buffer, sizeof(buffer));
    
  } else {
    i_command = INVALID_INPUT_CODE;
    return;
  } // endif

  // Update all dependent tolerances
  UpdateDependentTolerances();
}

#endif
