/*!\file Euler2DExactSolutions.cc
  \brief Source file initializing/implementing member variables/functions of class Euler2D_ExactSolutions. */

/* Include required C++ libraries. */
// None

/* Using std namespace functions */
// None

/* Include CFFC header files */
#include "Euler2DExactSolutions.h"
#include "Euler2DInput.h"
#include "../CFD/CFD.h"


// ==== Member variables ====
ExactSolutionBasicType * Euler2D_ExactSolutions::ExactSoln = NULL; //!< no associated exact solution
short Euler2D_ExactSolutions::i_Exact_Solution_Type = EULER2D_NO_EXACT_SOLUTION; //!< value for no exact solution

// ==== Member functions ====
/*! Default constructor is used to create the unique object of this class */
Euler2D_ExactSolutions::Euler2D_ExactSolutions(void){ }

/*! Destroy the exact solution object */
void Euler2D_ExactSolutions::DestroyExactSolutionObject(void){
  delete ExactSoln;
  ExactSoln = NULL;		// set to no associated exact solution
  i_Exact_Solution_Type = EULER2D_NO_EXACT_SOLUTION;	// set to value for no associated exact solution
}

/*!
 * Create and return the address of the exact solution.
 * This is a singleton class (only one instance is created).\n
 * DON'T modify this function unless you are familiar with the Meyers singleton pattern.
 */
Euler2D_ExactSolutions& Euler2D_ExactSolutions::getInstance(void){
  static Euler2D_ExactSolutions inst;
  std::atexit(DestroyExactSolutionObject);		// schedule for deallocation at exit
  return inst;
}

/*!
 * Set the exact solution based on the required solution index
 */
void Euler2D_ExactSolutions::SetExactSolution(const short &SolutionIndex){
  // deallocate the exact solution
  DestroyExactSolutionObject();

  // create the proper exact solution and set the index accordingly
  switch (SolutionIndex){
  case EULER2D_NO_EXACT_SOLUTION:
    // Don't do anything. The values were set by DestroyExactSolutionObject() routine.
    break;
  case EULER2D_EXACT_SOLUTION_RINGLEB_FLOW:
    ExactSoln = new Ringleb_Flow_ExactSolution;
    break;
  case EULER2D_EXACT_SOLUTION_ABGRALL_FUNCTION:
    ExactSoln = new Abgrall_Function_ExactSolution;
    break;
  case EULER2D_EXACT_SOLUTION_SINUSOIDAL_FUNCTION:
    ExactSoln = new Sinusoidal_Function_ExactSolution;
    break;
  case EULER2D_EXACT_SOLUTION_COSSIN_FUNCTION:
    ExactSoln = new CosSin_Function_ExactSolution;
    break;
  default:
    throw runtime_error("Euler2D_ExactSolutions::SetExactSolution() ERROR! Unknown exact solution type index.");
  }
  
  // Store exact solution type
  i_Exact_Solution_Type = SolutionIndex;
}

/*!
 * Parse the input control parameters for settings 
 * related to Euler2D_ExactSolutions class.
 * The first entry related to the exact solution in the
 * input file MUST specify the type of the exact solution.
 * Consequent parameters will be parsed with the parser 
 * particular to the specified exact solution type.
 */
void Euler2D_ExactSolutions::Parse_Next_Input_Control_Parameter(Euler2D_Input_Parameters & IP,
								int & i_command){
  // Check if the next control parameter has already been identified
  if (i_command != INVALID_INPUT_CODE){
    return;
  }

  char buffer[256];

  // Try to match the next control parameter
  if (strcmp(IP.Next_Control_Parameter, "Exact_Solution") == 0) {
    IP.Get_Next_Input_Control_Parameter();
    if ( strcmp(IP.Next_Control_Parameter, "Ringleb_Flow") == 0 ){
      SetExactSolution(EULER2D_EXACT_SOLUTION_RINGLEB_FLOW);
    } else if ( strcmp(IP.Next_Control_Parameter, "Abgrall_Function") == 0 ){
      SetExactSolution(EULER2D_EXACT_SOLUTION_ABGRALL_FUNCTION);     
    } else if ( strcmp(IP.Next_Control_Parameter, "Sinusoidal_Function") == 0 ){
      SetExactSolution(EULER2D_EXACT_SOLUTION_SINUSOIDAL_FUNCTION);      
    } else if ( strcmp(IP.Next_Control_Parameter, "CosSine_Function") == 0 ){
      SetExactSolution(EULER2D_EXACT_SOLUTION_COSSIN_FUNCTION);
    } else {
      i_command = INVALID_INPUT_CODE;
      return;
    }
    i_command = 0;

  } else {
    // Continue parsing with the parser of the current exact solution
    if (ExactSoln != NULL){
      ExactSoln->Parse_Next_Input_Control_Parameter(IP,i_command);
    }
    return;
  } // endif
}

/*!
 * Print the relevant parameters of the Euler2D_ExactSolutions class for the 
 * selected exact solution to the provided output stream.
 */
void Euler2D_ExactSolutions::Print_Info(std::ostream & out_file){

  out_file << "\n  -> Problem Exact Solution: " ; 
  if (ExactSoln != NULL){
    // Output field name
    out_file << ExactSoln->whatExactSolution();

    // Output exact solution characteristic parameters
    ExactSoln->Print_Info(out_file);
  } else {
    // There is no associated exact solution
    out_file << "Not specified";
  }

}

/*!
 * Broadcast the Euler2D_ExactSolutions variables to all      
 * processors associated with the specified communicator
 * from the specified processor using the MPI broadcast 
 * routine.
 */
void Euler2D_ExactSolutions::Broadcast(void){
#ifdef _MPI_VERSION
  
  short i_Exact_Solution_Type_Copy;

  // Copy the value of the i_Exact_Solution_Type on the main CPU
  if (CFFC_Primary_MPI_Processor()){
    i_Exact_Solution_Type_Copy = i_Exact_Solution_Type;
  }

  // Broadcast the type of the exact solution
  MPI::COMM_WORLD.Bcast(&i_Exact_Solution_Type_Copy,
			1, 
			MPI::SHORT, 0);

  // Create the same exact solution object on each CPU different than the primary one
  // using the broadcast value.
  if (!CFFC_Primary_MPI_Processor()){
      SetExactSolution(i_Exact_Solution_Type_Copy);
  }

  // Broadcast the characteristic parameters for the exact solution object
  if (IsExactSolutionSet()){
    ExactSoln->Broadcast();
  }
  
#endif
}

