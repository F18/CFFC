/*! \file    Grid3DHighOrderExecutionMode.h
 *  \brief   Initialize the flags that control the generation of grid
 *           level information used in the CENO reconstruction.
 */

#ifndef _GRID3D_HO_EXECUTIONMODE_INCLUDED
#define _GRID3D_HO_EXECUTIONMODE_INCLUDED

/* Include required C++ libraries. */
#include <sstream>

/* Using std namespace functions */
// None

/* Include CFFC header files */
// None

//! \class Grid3D_HO_Execution_Mode
//  -----------------------------------------------------------------------
/*!  
 *  Definition of flags that control the generation of grid level
 *  information used in the CENO reconstruction.
 *
 *///----------------------------------------------------------------------
class Grid3D_HO_Execution_Mode{
 
public:

  static void Print_Info(std::ostream & out_file);
  static void Broadcast(void);

  // Set all flags to default values
  static void SetDefaults(void);

  /*!
   *  This parameter controls the generation of grid level information
   *  required by CENO reconstruction.  Flag is set based on the
   *  Reconstruction_Type parsed in the CFD_Input_Parameters class.  
   *  - Turned ON if CENO reconstruction is used 
   *  - Turned OFF if any other reconstruction type is used (default)
   */
  static short USE_HO_CENO_GRID;

  /*!
   *  This paramter is set equal to the recontruction order (k).  It is used to
   *  set the container sizes for grid level information required by CENO
   *  reconstruction.  
   *  - Default is set to ZERO (ie. piecewise constant).
   */
  static int RECONSTRUCTION_ORDER;



protected:
  Grid3D_HO_Execution_Mode(void);   //!< Private default constructor
  Grid3D_HO_Execution_Mode(const Grid3D_HO_Execution_Mode&); //!< Private copy constructor
  Grid3D_HO_Execution_Mode& operator=(const Grid3D_HO_Execution_Mode&); //!< Private assignment operator

};

//! Parse the input control parameters for 
//  settings related to Grid3D_HO_Execution_Mode class
//template<class Input_Parameters_Type> inline
//void Grid3D_HO_Execution_Mode::Parse_Next_Input_Control_Parameter(Input_Parameters_Type & IP, int & i_command, 
//								  char *code, stringstream &value){
//
//  string value_string;
//
////  // Check if the next control parameter has already been identified
////    if (i_command != INVALID_INPUT_CODE){
////      return;
////    }
//
//  // Try to match the next control parameter
//
//    if (strcmp(code, "Reconstruction_Type") == 0) {
//      i_command = 50;
//      value >> value_string;
//      strcpy(IP.Reconstruction_Type_For_Grid_Info, value_string.c_str());
//      if (strcmp(IP.Reconstruction_Type_For_Grid_Info, "CENO") == 0) {
//	Grid3D_HO_Execution_Mode::USE_HO_CENO_GRID == ON;
//      } else {
//	Grid3D_HO_Execution_Mode::USE_HO_CENO_GRID == OFF;
//      }
//    }
////    else {
////      i_command = INVALID_INPUT_CODE;
////    } // endif
//}
#endif
