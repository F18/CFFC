/*!\file Tecplot_ExecutionMode.h
  \brief Definition of flags that control the data plotting with Tecplot software. */

#ifndef _TECPLOT_EXECUTIONMODE_INCLUDED
#define _TECPLOT_EXECUTIONMODE_INCLUDED

/* Include required C++ libraries. */
// None

/* Using std namespace functions */
// None

/* Include CFFC header files */
#include "../Utilities/Utilities.h"
#include "../CFD/CFD.h"

/*!
 * \class Tecplot_Execution_Mode
 *
 * @brief Definition of flags that control the data plotting with Tecplot.
 *
 */
class Tecplot_Execution_Mode{
  
public:
  //! Data type defining different levels of output.
  enum Output_Level_Type { Brief,    //!< Minimum level of output
			   Detailed, //!< A more detailed output
			   Full    , //!< Complete level of output
			   Extended  //!< Complete level plus extra information, not necessarily part of the solution
  };
  
  //! Data type defining different levels of accuracy for the output data
  enum Output_Accuracy_Type { SinglePrec, //!< define single precision
			      DoublePrec  //!< define double precision 
  };
  
  static short OUTPUT_LEVEL;
  static short OUTPUT_ACCURACY;

  // set all flags to default values
  static void SetDefaults(void);

  template<class Input_Parameters_Type>
  static void Parse_Next_Input_Control_Parameter(Input_Parameters_Type & IP, int & i_command);

  static void Print_Info(std::ostream & out_file);

  static void Broadcast(void);

  //! Set the designated switch to single accuracy precision
  static void setDoublePrecisionPlotting(void){ OUTPUT_ACCURACY = SinglePrec; }

  //! Return true if the data must be plotted with double accuracy otherwise false
  static bool IsDoublePrecision(void){ return (OUTPUT_ACCURACY == DoublePrec) ? true: false; }
  //! Return true if the data associated with the detailed format should be output otherwise false
  static bool IsDetailedOutputRequired(void){ return ((OUTPUT_LEVEL == Detailed) ||
						      (OUTPUT_LEVEL == Full) || 
						      (OUTPUT_LEVEL == Extended) ) ? true: false;}
  //! Return true if the data associated with the full format should be output otherwise false
  static bool IsFullOutputRequired(void){ return ( (OUTPUT_LEVEL == Full) ||
						   (OUTPUT_LEVEL == Extended) ) ? true: false; }
  //! Return true if the data associated with the extended format should be output otherwise false
  static bool IsExtendedOutputRequired(void){ return (OUTPUT_LEVEL == Extended) ? true: false; }
  
protected:
  Tecplot_Execution_Mode(void);   //!< Private default constructor
  Tecplot_Execution_Mode(const Tecplot_Execution_Mode&); //!< Private copy constructor
  Tecplot_Execution_Mode& operator=(const Tecplot_Execution_Mode&); //!< Private assignment operator

};

//! Parse the input control parameters for 
//  settings related to Tecplot_Execution_Mode class
template<class Input_Parameters_Type> inline
void Tecplot_Execution_Mode::Parse_Next_Input_Control_Parameter(Input_Parameters_Type & IP, int & i_command){

  // Check if the next control parameter has already been identified
  if (i_command != INVALID_INPUT_CODE){
    return;
  }

  char buffer[256];

  // Try to match the next control parameter
  if (strcmp(IP.Next_Control_Parameter, "Tecplot_Data_Precision") == 0) {
    IP.Get_Next_Input_Control_Parameter();
    if ( strcmp(IP.Next_Control_Parameter, "Single") == 0 ){
      OUTPUT_ACCURACY = SinglePrec;
    } else if ( strcmp(IP.Next_Control_Parameter, "Double") == 0 ) {
      OUTPUT_ACCURACY = DoublePrec;
    } else {
      i_command = INVALID_INPUT_VALUE;
      return;
    }
    i_command = 0;

  } else if (strcmp(IP.Next_Control_Parameter, "Tecplot_Output_Format") == 0) {
    IP.Get_Next_Input_Control_Parameter();
    if ( strcmp(IP.Next_Control_Parameter, "Brief") == 0 ){
      OUTPUT_LEVEL = Brief;
    } else if ( strcmp(IP.Next_Control_Parameter, "Detailed") == 0 ) {
      OUTPUT_LEVEL = Detailed;
    } else if ( strcmp(IP.Next_Control_Parameter, "Full") == 0 ) {
      OUTPUT_LEVEL = Full;
    } else if ( strcmp(IP.Next_Control_Parameter, "Extended") == 0 ) {
      OUTPUT_LEVEL = Extended;
    } else {
      i_command = INVALID_INPUT_VALUE;
      return;
    }
    i_command = 0;

  } else {
    i_command = INVALID_INPUT_CODE;
  } // endif

}

#endif // _TECPLOT_EXECUTIONMODE_INCLUDED
