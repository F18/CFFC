/* AMRInput.cc: Definition of AMR_Input_Parameters class member functions. */

/* Include the AMRInput header file. */

#ifndef _AMRINPUT_INCLUDED
#include "AMRInput.h"
#endif // _AMRINPUT_INCLUDED

/* Define member functions. */

/***************************************************************************
 * AMR_Input_Parameters::Broadcast -- Broadcast to all processors.         *
 ***************************************************************************/
void AMR_Input_Parameters::Broadcast(void) {

#ifdef _MPI_VERSION
    // Multi-block solution-adaption and parallel domain decomposition input parameters:
    if (!CFFC_Primary_MPI_Processor()) {
       Number_of_Processors = CFFC_MPI::Number_of_Processors;
    } /* endif */
    MPI::COMM_WORLD.Bcast(&(Number_of_Blocks_Per_Processor),
                          1,
                          MPI::INT, 0);

    // AMR input parameters:
    MPI::COMM_WORLD.Bcast(&(AMR),
                          1,
                          MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(&(AMR_Frequency),
                          1,
                          MPI::INT,0);
    MPI::COMM_WORLD.Bcast(&(Number_of_Initial_Mesh_Refinements),
                          1,
                          MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(&(Number_of_Uniform_Mesh_Refinements),
                          1,
                          MPI::INT,0);
    MPI::COMM_WORLD.Bcast(&(Number_of_Boundary_Mesh_Refinements),
                          1,
                          MPI::INT,0);
    MPI::COMM_WORLD.Bcast(&(Maximum_Refinement_Level),
                          1,
                          MPI::INT,0);
    MPI::COMM_WORLD.Bcast(&(Minimum_Refinement_Level),
                          1,
                          MPI::INT,0);
    MPI::COMM_WORLD.Bcast(&(Threshold_for_Refinement),
                          1,
                          MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(Threshold_for_Coarsening),
                          1,
                          MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(Number_of_Refinement_Criteria),
                          1,
                          MPI::INT,0);

    // Morton ordering input parameters:
    MPI::COMM_WORLD.Bcast(&(Morton),
                          1,
                          MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(&(Morton_Reordering_Frequency),
                          1,
                          MPI::INT,0);
#endif

}

/***************************************************************************
 * AMR_Input_Parameters::Parse_Next_Input_Control_Parameter - Parse input. *
 ***************************************************************************/
int AMR_Input_Parameters::Parse_Next_Input_Control_Parameter(char *code, 
                                                             stringstream &value) {

// Returns:
//  - INVALID_INPUT_VALUE if code is valid but value is invalid
//  - INVALID_INPUT_CODE  if unknown code

  int i_command = INVALID_INPUT_CODE;
  string value_string;

  if (strcmp(code, "Number_of_Blocks_Per_Processor") == 0) {
     i_command = 4001;
     value >> Number_of_Blocks_Per_Processor;
     if (Number_of_Blocks_Per_Processor < 1) i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(code, "AMR") == 0) {
    i_command = 4002;
    value >> value_string;
    if (value_string == "OFF" || 
        value_string == "0") {
       AMR = OFF;
    } else if (value_string == "ON" || 
               value_string == "1") {
       AMR = ON;
    } else {
       i_command = INVALID_INPUT_VALUE;
    } /* endif */

  } else if (strcmp(code, "AMR_Frequency") == 0) {
    i_command = 4003;
    value >> AMR_Frequency;
    if (AMR_Frequency < 0) i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(code, "Number_of_Initial_Mesh_Refinements") == 0) {
    i_command = 4004;
    value >> Number_of_Initial_Mesh_Refinements;
    if (Number_of_Initial_Mesh_Refinements < 0) Number_of_Initial_Mesh_Refinements = 0;

  } else if (strcmp(code,"Number_of_Uniform_Mesh_Refinements") == 0) {
    i_command = 4005;
    value >> Number_of_Uniform_Mesh_Refinements;
    if (Number_of_Uniform_Mesh_Refinements < 0) Number_of_Uniform_Mesh_Refinements = 0;

  } else if (strcmp(code,"Number_of_Boundary_Mesh_Refinements") == 0) {
    i_command = 4006;
    value >> Number_of_Boundary_Mesh_Refinements;
    if (Number_of_Boundary_Mesh_Refinements < 0) Number_of_Boundary_Mesh_Refinements = 0;

  } else if (strcmp(code,"Maximum_Refinement_Level") == 0) {
    i_command = 4007;
    value >> Maximum_Refinement_Level;
    if (Maximum_Refinement_Level < 1) Maximum_Refinement_Level = 1;

  } else if (strcmp(code,"Minimum_Refinement_Level") == 0) {
    i_command = 4008;
    value >> Minimum_Refinement_Level;
    if (Minimum_Refinement_Level < 1) Minimum_Refinement_Level = 1;

  } else if (strcmp(code, "Threshold_for_Refinement") == 0) {
    i_command = 4009;
    value >> Threshold_for_Refinement;
    if (Threshold_for_Refinement <= ZERO ||
        Threshold_for_Refinement > ONE) Threshold_for_Refinement = 0.50;

  } else if (strcmp(code, "Threshold_for_Coarsening") == 0) {
    i_command = 4010;
    value >> Threshold_for_Coarsening;
    if (Threshold_for_Coarsening < ZERO ||
  	Threshold_for_Coarsening >= ONE) Threshold_for_Coarsening = 0.10;

  } else if (strcmp(code,"Number_of_Refinement_Criteria") == 0) {
    i_command = 4011;
    value >> Number_of_Refinement_Criteria;
    if (Number_of_Refinement_Criteria < 1) Number_of_Refinement_Criteria = 1;

  } else if (strcmp (code, "Morton") == 0) {
    i_command = 4012;
    value >> value_string;
    if (value_string == "OFF" || 
        value_string == "0") {
       Morton = OFF;
    } else if (value_string == "ON" || 
               value_string == "1") {
       Morton = ON;
    } else {
       i_command = INVALID_INPUT_VALUE;
    } /* endif */

  } else if (strcmp(code, "Morton_Reordering_Frequency") == 0) {
    i_command = 4013;
    value >> Morton_Reordering_Frequency;
    if (Morton_Reordering_Frequency < 0) i_command = INVALID_INPUT_VALUE;

  } else {
    i_command = INVALID_INPUT_CODE;

  } /* endif */

  return i_command;
  
}

/***************************************************************************
 * AMR_Input_Parameters::Check_Inputs -- Check input values.               *
 ***************************************************************************/
int AMR_Input_Parameters::Check_Inputs(void) {

  // Input parameters are consistent.  Exit successfully.
  return 0;

}

/***************************************************************************
 * AMR_Input_Parameters -- Input-output operators.                         *
 ***************************************************************************/
ostream &operator << (ostream &out_file,
                      const AMR_Input_Parameters &IP) {

  IP.Output(out_file);
  return (out_file);

}

istream &operator >> (istream &in_file,
                      AMR_Input_Parameters &IP) {

  return in_file;

}

void AMR_Input_Parameters::Output(ostream &out_file) const {

  out_file << "\n  -> Number of Processors: " 
           << Number_of_Processors;
  out_file << "\n  -> Number of Blocks Per Processor: " 
           << Number_of_Blocks_Per_Processor;
  if (AMR) out_file << "\n  -> AMR Frequency: "
                    << AMR_Frequency
                    << " steps (iterations)";
  if (Number_of_Initial_Mesh_Refinements > 0)
     out_file << "\n  -> Number of Initial Mesh Refinements : " 
              << Number_of_Initial_Mesh_Refinements;
  if (Number_of_Uniform_Mesh_Refinements > 0)
     out_file << "\n  -> Number of Uniform Mesh Refinements : " 
  	      << Number_of_Uniform_Mesh_Refinements;
  if (Number_of_Boundary_Mesh_Refinements > 0)
    out_file << "\n  -> Number of Boundary Mesh Refinements : " 
  	     << Number_of_Boundary_Mesh_Refinements;
  out_file << "\n  -> Minimum Refinement Level: "
           << Minimum_Refinement_Level;
  out_file << "\n  -> Maximum Refinement Level: "
           << Maximum_Refinement_Level;
  if (Morton) out_file << "\n  -> Morton Re-Ordering Frequency: "
                       << Morton_Reordering_Frequency
                       << " steps (iterations)";

}
