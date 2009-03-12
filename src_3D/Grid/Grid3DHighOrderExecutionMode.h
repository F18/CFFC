/*!\file Grid3DHighOrderExecutionMode.h
  \brief Definition of flags that control the generation of grid level information required by CENO high-order reconstruction. */

#ifndef _GRID3D_HO_EXECUTIONMODE_INCLUDED
#define _GRID3D_HO_EXECUTIONMODE_INCLUDED

/* Include required C++ libraries. */
#include <sstream>

/* Using std namespace functions */
// None

/* Include CFFC header files */

/*!%
 * Class: Grid3D_HO_Execution_Mode
 *
 * @brief Definition of flags that control the generation of grid level information required by CENO high-order reconstruction. 
 *
 */
class Grid3D_HO_Execution_Mode{
  
public:

  static void Print_Info(std::ostream & out_file);
  //  static void Broadcast(void);

  // set all flags to default values
  static void SetDefaults(void);

  /*! Generation of grid level information required by CENO reconstruction (high or low order).\n
      Set automatically in CFD_Input_Parameters based on the input Reconstruction_Type: \n
      Turned ON if CENO reconstruction is used \n
      Turn OFF if any other reconstruction type is used (default) \n
      ---------------------------------------------------------------------------------------- */
  static short USE_HO_CENO_GRID;

  /*! Static integer variable which equals the recontruction order. \n
      Used to set container sizes for grid level information required by CENO reocnstruction. \n
      Default is set to 0 (or piecewise constant).
      ---------------------------------------------------------------------------------------- */
  static int RECONSTRUCTION_ORDER_FOR_GRID_INFO;



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
