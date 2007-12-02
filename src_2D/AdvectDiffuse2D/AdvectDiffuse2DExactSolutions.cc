/*!\file AdvectDiffuse2DExactSolutions.cc
  \brief Source file initializing/implementing member variables/functions of class AdvectDiffuse2D_ExactSolutions. */

/* Include required C++ libraries. */
// None

/* Using std namespace functions */
// None

/* Include CFFC header files */
#include "AdvectDiffuse2DExactSolutions.h"
#include "AdvectDiffuse2DInput.h"
#include "../CFD/CFD.h"


// ==== Member variables ====
ExactSolutionBasicType * AdvectDiffuse2D_ExactSolutions::ExactSoln = NULL; //!< no associated exact solution
short AdvectDiffuse2D_ExactSolutions::i_Exact_Solution_Type = -1;          //!< value for no associated exact solution

// ==== Member functions ====
/*! Default constructor is used to create the unique object of this class */
AdvectDiffuse2D_ExactSolutions::AdvectDiffuse2D_ExactSolutions(void){ }

/*! Destroy the exact solution object */
void AdvectDiffuse2D_ExactSolutions::DestroyExactSolutionObject(void){
  delete ExactSoln;
  ExactSoln = NULL;		// set to no associated exact solution
  i_Exact_Solution_Type = -1;	// set to value for no associated exact solution
}

/*!
 * Create and return the address of the exact solution.
 * This is a singleton class (only one instance is created).\n
 * DON'T modify this function unless you are familiar with the Meyers singleton pattern.
 */
AdvectDiffuse2D_ExactSolutions& AdvectDiffuse2D_ExactSolutions::getInstance(void){
  static AdvectDiffuse2D_ExactSolutions inst;
  std::atexit(DestroyExactSolutionObject);		// schedule for deallocation at exit
  return inst;
}

/*!
 * Set the exact solution based on the required solution index
 */
void AdvectDiffuse2D_ExactSolutions::SetExactSolution(const short &SolutionIndex){
  // deallocate the exact solution
  DestroyExactSolutionObject();

  // create the proper exact solution and set the index accordingly
  switch (SolutionIndex){
  case AD2D_EXACT_SOLUTION_LAPLACE_I:
    ExactSoln = new Laplace_I_ExactSolution;
    break;
    
  case AD2D_EXACT_SOLUTION_LAPLACE_II:
    ExactSoln = new Laplace_II_ExactSolution;
    break;
    
  case AD2D_EXACT_SOLUTION_LAPLACE_III:
    ExactSoln = new Laplace_III_ExactSolution;
    break;

  case AD2D_EXACT_SOLUTION_LAPLACE_IV:
    ExactSoln = new Laplace_IV_ExactSolution;
    break;

  case AD2D_EXACT_SOLUTION_LAPLACE_V:
    ExactSoln = new Laplace_V_ExactSolution;
    break;

  case AD2D_EXACT_SOLUTION_POISSON_I:
    ExactSoln = new Poisson_I_ExactSolution;
    break;
    
  case AD2D_EXACT_SOLUTION_POISSON_II:
    ExactSoln = new Poisson_II_ExactSolution;
    break;
    
  case AD2D_EXACT_SOLUTION_POISSON_III:
    ExactSoln = new Poisson_III_ExactSolution;
    break;

  case AD2D_EXACT_SOLUTION_POISSON_IV:
    ExactSoln = new Poisson_IV_ExactSolution;
    break;

  case AD2D_EXACT_SOLUTION_POISSON_V:
    ExactSoln = new Poisson_V_ExactSolution;
    break;

  case AD2D_EXACT_SOLUTION_STATIONARY_HEAT_TRANSFER_WITH_LINEAR_SOURCE:
    ExactSoln = new StationaryHeatEqnWithLinearSource_ExactSolution;
    break;

  case AD2D_EXACT_SOLUTION_PURE_CIRCULAR_ADVECTION_AT_CONSTANT_SPIN:
    ExactSoln = new PureCircularAdvectionAtConstantSpin_ExactSolution;
    break;

  case AD2D_EXACT_SOLUTION_ADVECTION_DIFFUSION_IN_RECTANGULAR_CHANNEL:
    ExactSoln = new AdvectionDiffusionInRectangularChannel_ExactSolution;
    break;
  
  default:
    throw runtime_error("AdvectDiffuse2D_ExactSolutions::SetExactSolution() ERROR! Unknown exact solution type index.");
  }
  
  // Store exact solution type
  i_Exact_Solution_Type = SolutionIndex;
}

/*!
 * Parse the input control parameters for settings 
 * related to AdvectDiffuse2D_ExactSolutions class.
 * The first entry related to the exact solution in the
 * input file MUST specify the type of the exact solution.
 * Consequent parameters will be parsed with the parser 
 * particular to the specified exact solution type.
 */
void AdvectDiffuse2D_ExactSolutions::Parse_Next_Input_Control_Parameter(AdvectDiffuse2D_Input_Parameters & IP,
									int & i_command){
  // Check if the next control parameter has already been identified
  if (i_command != INVALID_INPUT_CODE){
    return;
  }

  char buffer[256];

  // Try to match the next control parameter
  if (strcmp(IP.Next_Control_Parameter, "Exact_Solution") == 0) {
    IP.Get_Next_Input_Control_Parameter();
    if ( strcmp(IP.Next_Control_Parameter, "Laplace_I") == 0 ){
      SetExactSolution(AD2D_EXACT_SOLUTION_LAPLACE_I);
    } else if ( strcmp(IP.Next_Control_Parameter, "Laplace_II") == 0 ) {
      SetExactSolution(AD2D_EXACT_SOLUTION_LAPLACE_II);
    } else if ( strcmp(IP.Next_Control_Parameter, "Laplace_III") == 0 ) {
      SetExactSolution(AD2D_EXACT_SOLUTION_LAPLACE_III);
    } else if ( strcmp(IP.Next_Control_Parameter, "Laplace_IV") == 0 ) {
      SetExactSolution(AD2D_EXACT_SOLUTION_LAPLACE_IV);
    } else if ( strcmp(IP.Next_Control_Parameter, "Laplace_V") == 0 ) {
      SetExactSolution(AD2D_EXACT_SOLUTION_LAPLACE_V);
    } else if ( strcmp(IP.Next_Control_Parameter, "Poisson_I") == 0 ){
      SetExactSolution(AD2D_EXACT_SOLUTION_POISSON_I);
    } else if ( strcmp(IP.Next_Control_Parameter, "Poisson_II") == 0 ) {
      SetExactSolution(AD2D_EXACT_SOLUTION_POISSON_II);
    } else if ( strcmp(IP.Next_Control_Parameter, "Poisson_III") == 0 ) {
      SetExactSolution(AD2D_EXACT_SOLUTION_POISSON_III);
    } else if ( strcmp(IP.Next_Control_Parameter, "Poisson_IV") == 0 ) {
      SetExactSolution(AD2D_EXACT_SOLUTION_POISSON_IV);
    } else if ( strcmp(IP.Next_Control_Parameter, "Poisson_V") == 0 ) {
      SetExactSolution(AD2D_EXACT_SOLUTION_POISSON_V);
    } else if ( strcmp(IP.Next_Control_Parameter, "Heat_Diffusion_With_Linear_Source") == 0 ) {
      SetExactSolution(AD2D_EXACT_SOLUTION_STATIONARY_HEAT_TRANSFER_WITH_LINEAR_SOURCE);
    } else if ( strcmp(IP.Next_Control_Parameter, "Pure_Circular_Advection") == 0 ) {
      SetExactSolution(AD2D_EXACT_SOLUTION_PURE_CIRCULAR_ADVECTION_AT_CONSTANT_SPIN);
    } else if ( strcmp(IP.Next_Control_Parameter, "Advection_Diffusion_In_Rectangular_Channel") == 0 ) {
      SetExactSolution(AD2D_EXACT_SOLUTION_ADVECTION_DIFFUSION_IN_RECTANGULAR_CHANNEL);
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
 * Print the relevant parameters of the AdvectDiffuse2D_ExactSolutions class for the 
 * selected exact solution to the provided output stream.
 */
void AdvectDiffuse2D_ExactSolutions::Print_Info(std::ostream & out_file){

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
 * Broadcast the AdvectDiffuse2D_ExactSolutions variables to all      
 * processors associated with the specified communicator
 * from the specified processor using the MPI broadcast 
 * routine.
 */
void AdvectDiffuse2D_ExactSolutions::Broadcast(void){
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
  ExactSoln->Broadcast();
  
#endif
}

