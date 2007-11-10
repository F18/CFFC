/*!\file SourceFields.h
  \brief Header file defining 2D source fields. */

#ifndef _SOURCE_FIELDS_INCLUDED
#define _SOURCE_FIELDS_INCLUDED

/* Include required C++ libraries. */
// None

/* Using std namespace functions */
// None

/* Include CFFC header files */
#include "../Math/Math.h"
#include "../Utilities/Utilities.h"
#include "AdvectDiffuse2DInput.h"

/*! 
 * \class SourceFieldBasicType
 *
 * Basic data type for any source field.
 * This is an ADT (abstract data type).
 */
class SourceFieldBasicType{
public:

  //! @name Virtual member functions
  //@{
  //! Declare a pure virtual destructor
  virtual ~SourceFieldBasicType(void) = 0;

  /*! 
   * Calculate the non-linear field solution based on the location and the solution value
   */
  virtual double FieldSoln(const double &x, const double &y, const double & Soln) = 0;

  /*! 
   * Calculate the non-linear field solution based on the average solution value
   */
  virtual double FieldSoln(const double &SolnAvg) = 0;

  //! This function returns true is the solution needs to be integrated over the cell domain and false if it does not
  virtual bool FieldRequireIntegration(void) = 0;

  //! Parse the input control parameters
  virtual void Parse_Next_Input_Control_Parameter(AdvectDiffuse2D_Input_Parameters & IP, int & i_command) = 0;
  //@}
};

/*! 
 * \class SourceFieldThatRequireIntegration
 *
 * Basic data type for source fields that require numerical integration of the field solution.
 */
class SourceFieldThatRequireIntegration: public SourceFieldBasicType{
public:

  //! This field require integration (information used at higher level)
  bool FieldRequireIntegration(void){ return true; }

  //! Declare a pure virtual destructor
  virtual ~SourceFieldThatRequireIntegration(void) = 0;

  //! Throw a runtime error if this function is called for this class
  double FieldSoln(const double &SolnAvg){
    throw runtime_error("SourceFieldThatRequireIntegration::FieldSoln() ERROR! Wrong number of parameters!");
    return 0.0;
  }
};

/*! 
 * \class SourceFieldThatDoesNotRequireIntegration
 *
 * Basic data type for source fields that DON'T require numerical integration of the field solution.
 * The integral of the field solution is based on the average solution (this includes the area dependency)
 */
class SourceFieldThatDoesNotRequireIntegration: public SourceFieldBasicType{
public:

  //! This field DOES NOT require integration (information used at higher level)
  bool FieldRequireIntegration(void){ return false; }

  //! Declare a pure virtual destructor
  virtual ~SourceFieldThatDoesNotRequireIntegration(void) = 0;

  //! Throw a runtime error if this function is called for this class
  virtual double FieldSoln(const double &x, const double &y, const double & Soln){
    throw runtime_error("SourceFieldThatDoesNotRequireIntegration::FieldSoln() ERROR! Wrong number of parameters!");
    return 0.0;
  }
};


/*! 
 * \class ZERO_SourceField
 * 
 * Implements a ZERO value source field.
 */
class ZERO_SourceField: public SourceFieldThatDoesNotRequireIntegration{
public:

  //! Return ZERO solution (no source)
  double FieldSoln(const double &SolnAvg){return 0.0;}

  //! Parse the input control parameters
  void Parse_Next_Input_Control_Parameter(AdvectDiffuse2D_Input_Parameters & IP, int & i_command){};
};

/*! 
 * \class Linear_SourceField
 * 
 * Implements a linear source field.
 * 
 */
class Linear_SourceField: public SourceFieldThatDoesNotRequireIntegration{
public:

  Linear_SourceField(void): tau(1.0){ };

  //! Return linear solution 
  double FieldSoln(const double &SolnAvg);

  //! Parse the input control parameters
  void Parse_Next_Input_Control_Parameter(AdvectDiffuse2D_Input_Parameters & IP, int & i_command);

private:
  double tau;
};

inline double Linear_SourceField::FieldSoln(const double &SolnAvg){
  return SolnAvg/tau;
}

inline void Linear_SourceField::Parse_Next_Input_Control_Parameter(AdvectDiffuse2D_Input_Parameters & IP,
								   int & i_command){

  // Check if the next control parameter has already been identified
  if (i_command != INVALID_INPUT_CODE){
    return;
  }

  char buffer[256];

  // Try to match the next control parameter
  if (strcmp(IP.Next_Control_Parameter, "Source_Linear_Tau_Coeff") == 0) {
    i_command = 0;
    ++IP.Line_Number;
    IP.Input_File >> tau;
    IP.Input_File.getline(buffer, sizeof(buffer));

  } else {
    i_command = INVALID_INPUT_CODE;
    return;
  } // endif
}

/*! 
 * \class Exponential_SourceField
 * 
 * Non-linear exponential source field.
 * 
 */
class Exponential_SourceField: public SourceFieldThatRequireIntegration{
public:
  
  Exponential_SourceField(void): a(1.0), beta(1.0){ };

  //! Return exponential non-linear solution
  double FieldSoln(const double &x, const double &y, const double & Soln);

  //! Parse the input control parameters
  void Parse_Next_Input_Control_Parameter(AdvectDiffuse2D_Input_Parameters & IP, int & i_command);

private:
  double a, beta;
};

inline double Exponential_SourceField::FieldSoln(const double &x, const double &y, const double & Soln){
  return a*exp(beta*Soln);
}

inline void Exponential_SourceField::Parse_Next_Input_Control_Parameter(AdvectDiffuse2D_Input_Parameters & IP,
									int & i_command){

  // Check if the next control parameter has already been identified
  if (i_command != INVALID_INPUT_CODE){
    return;
  }

  char buffer[256];

  // Try to match the next control parameter
  if (strcmp(IP.Next_Control_Parameter, "Source_Exponential_A") == 0) {
    i_command = 0;
    ++IP.Line_Number;
    IP.Input_File >> a;
    IP.Input_File.getline(buffer, sizeof(buffer));

  } else if (strcmp(IP.Next_Control_Parameter, "Source_Exponential_Beta") == 0) {
    i_command = 0;
    ++IP.Line_Number;
    IP.Input_File >> beta;
    IP.Input_File.getline(buffer, sizeof(buffer));

  } else {
    i_command = INVALID_INPUT_CODE;
    return;
  } // endif
}

#endif
