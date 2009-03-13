/*!\file Grid3DHighOrderExecutionMode.cc
  \brief Initialize the flags that control the generation of grid level information required by CENO high-order reconstruction. */

/* Include required C++ libraries. */
// None

/* Using std namespace functions */
// None

/* Include CFFC header files */
# include "../CFD/CFD.h"
#include "Grid3DHighOrderExecutionMode.h"

short Grid3D_HO_Execution_Mode::USE_HO_CENO_GRID = OFF;
int Grid3D_HO_Execution_Mode::RECONSTRUCTION_ORDER_FOR_GRID_INFO = 0;

//! Set all flags to default values
// add all flag default values to this function
void Grid3D_HO_Execution_Mode::SetDefaults(void){
  
  USE_HO_CENO_GRID = OFF;
  RECONSTRUCTION_ORDER_FOR_GRID_INFO = 0;
}

//! Print the current execution mode
//  at the output stream
// \param [in] out_file the output stream
void Grid3D_HO_Execution_Mode::Print_Info(std::ostream & out_file){

  // output execution mode
  if (USE_HO_CENO_GRID == ON){
    out_file << "\n     -> Grid Execution Mode: " << "High Order Grid Information Generated and Stored";
  }
}

//! Parse the input control parameters for 
//  settings related to Grid3D_HO_Execution_Mode class
//template<class Input_Parameters_Type> 
//void Grid3D_HO_Execution_Mode::Parse_Next_Input_Control_Parameter(Input_Parameters_Type & IP, int & i_command, 
//								  char *code, stringstream &value){
//
//  string value_string;
//
//  // Check if the next control parameter has already been identified
//  if (i_command != INVALID_INPUT_CODE){
//    return;
//  }
//  
//  // Try to match the next control parameter
//  
//  if (strcmp(code, "Reconstruction_Type") == 0) {
//    i_command = 50;
//    value >> value_string;
//    strcpy(IP.Reconstruction_Type_For_Grid_Info, value_string.c_str());
//    if (strcmp(IP.Reconstruction_Type_For_Grid_Info, "CENO") == 0) {
//      Grid3D_HO_Execution_Mode::USE_HO_CENO_GRID == ON;
//    } else {
//      Grid3D_HO_Execution_Mode::USE_HO_CENO_GRID == OFF;
//    }
//  } else {
//    i_command = INVALID_INPUT_CODE;
//  } // endif
//}


/*!
 * Broadcast the Grid3D_HO_Execution_Mode variables to all      
 * processors associated with the specified communicator
 * from the specified processor using the MPI broadcast 
 * routine.
 */
//void Grid3D_HO_Execution_Mode::Broadcast(void){
//#ifdef _MPI_VERSION
//  
//  MPI::COMM_WORLD.Bcast(&USE_CENO_ALGORITHM,
// 			1, 
// 			MPI::SHORT, 0);
//  MPI::COMM_WORLD.Bcast(&USE_PSEUDO_INVERSE,
// 			1, 
// 			MPI::SHORT, 0);
//  MPI::COMM_WORLD.Bcast(&CENO_DROP_ORDER,
// 			1, 
// 			MPI::SHORT, 0);
//  MPI::COMM_WORLD.Bcast(&CENO_PADDING,
// 			1, 
// 			MPI::SHORT, 0);
//  MPI::COMM_WORLD.Bcast(&CENO_APPLY_GEOMETRIC_WEIGHTING,
// 			1, 
// 			MPI::SHORT, 0);
//  MPI::COMM_WORLD.Bcast(&CENO_SQUARE_GEOM_WEIGHTING,
// 			1, 
// 			MPI::SHORT, 0);
//  MPI::COMM_WORLD.Bcast(&CENO_CONSIDER_WEIGHTS,
// 			1, 
// 			MPI::SHORT, 0);
//  MPI::COMM_WORLD.Bcast(&FORCE_WITH_PIECEWISE_CONSTANT_AT_INTERFACE,
// 			1, 
// 			MPI::SHORT, 0);
//  MPI::COMM_WORLD.Bcast(&CENO_RECONSTRUCTION_WITH_MESSAGE_PASSING,
// 			1, 
// 			MPI::SHORT, 0);
//  MPI::COMM_WORLD.Bcast(&CENO_SMOOTHNESS_INDICATOR_COMPUTATION_WITH_ONLY_FIRST_NEIGHBOURS,
// 			1, 
// 			MPI::SHORT, 0);
//  MPI::COMM_WORLD.Bcast(&CENO_CONSTRAINED_RECONSTRUCTION_WITH_ADDITIONAL_APPROXIMATE_CONSTRAINTS,
// 			1, 
// 			MPI::SHORT, 0);
//  MPI::COMM_WORLD.Bcast(&CENO_CONSTRAINED_RECONSTRUCTION_WITH_EXTENDED_BIASED_STENCIL,
// 			1, 
// 			MPI::SHORT, 0);
//  MPI::COMM_WORLD.Bcast(&USE_LAPACK_LEAST_SQUARES,
// 			1, 
// 			MPI::SHORT, 0);
//  MPI::COMM_WORLD.Bcast(&IGNORE_CURVED_BOUNDARIES_FOR_ACCURACY_ASSESSMENT,
// 			1, 
// 			MPI::SHORT, 0);
//  MPI::COMM_WORLD.Bcast(&USE_SMOOTHNESS_INDICATOR_FOR_AMR_CRITERIA,
// 			1, 
// 			MPI::SHORT, 0);
//  MPI::COMM_WORLD.Bcast(&HIGH_ORDER_MESSAGE_PASSING,
// 			1, 
// 			MPI::SHORT, 0);
//  MPI::COMM_WORLD.Bcast(&CENO_VERBOSE,
// 			1, 
// 			MPI::SHORT, 0);
//  MPI::COMM_WORLD.Bcast(&Limiter,
// 			1, 
// 			MPI::SHORT, 0);  
//
//#endif
//}
