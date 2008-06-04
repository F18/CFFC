/*!\file HighOrder2DInput.h
  \brief Header file defining 2D high-order input parameters. */

#ifndef _HIGHORDER_2D_INPUT_INCLUDED
#define _HIGHORDER_2D_INPUT_INCLUDED

/* Include required C++ libraries. */
// None

/* Using std namespace functions */
using std::vector;

/* Include CFFC header files */
#include "../CFD/CFD.h"
#include "../Utilities/Utilities.h"

// Declare HighOrder2D_MultiBlock class which is a friend class of HighOrder2D_Input
class HighOrder2D_MultiBlock;

/*!
 * \class HighOrder2D_Input
 *
 * @brief Parsing class for high-order 2D input parameters.
 * \nosubgrouping
 */
class HighOrder2D_Input{
public:

  //! @name Field access
  //@{
  static const int & MaximumReconstructionOrder(void) { return MaxReconstructionOrder; }
  //@}

  //! @name Functions for input-output and broadcast
  //@{
  template<class Input_Parameters_Type>
  static void Parse_Next_Input_Control_Parameter(Input_Parameters_Type & IP, int & i_command);
  template<class Input_Parameters_Type>
  static void Set_Final_Parameters(Input_Parameters_Type & IP);
  static void Print_Info(std::ostream & out_file);
  static void Broadcast(void);
  //@}

  //! Set default values
  static void SetDefaults(void);

protected:
  HighOrder2D_Input(void);   //!< Private default constructor
  HighOrder2D_Input(const HighOrder2D_Input&); //!< Private copy constructor
  HighOrder2D_Input& operator=(const HighOrder2D_Input&); //!< Private assignment operator

  static vector<int> OrdersOfReconstruction;
  static int NumberOfAuxiliaryReconstructions;
  static int NumberOfHighOrderReconstructions;
  static int MaxReconstructionOrder;

  // declaration of friend classes
  friend class HighOrder2D_MultiBlock;
};

//! Parse the input control parameters for 
//  settings related to HighOrder2D_Input class
template<class Input_Parameters_Type> inline
void HighOrder2D_Input::Parse_Next_Input_Control_Parameter(Input_Parameters_Type & IP, int & i_command){

  // Check if the next control parameter has already been identified
  if (i_command != INVALID_INPUT_CODE){
    return;
  }

  char buffer[256];

  // Try to match the next control parameter
  if (strcmp(IP.Next_Control_Parameter, "Number_of_Auxiliary_Reconstructions") == 0 || 
      strcmp(IP.Next_Control_Parameter, "Number_Of_Auxiliary_Reconstructions") == 0) {
    i_command = 0;
    ++IP.Line_Number;
    IP.Input_File >> NumberOfAuxiliaryReconstructions;
    IP.Input_File.getline(buffer, sizeof(buffer));

    // Check for positive value
    if (NumberOfAuxiliaryReconstructions < 0){
      i_command = INVALID_INPUT_VALUE;
      return;
    }

    // Allocate memory for reconstruction orders ( The first rec. order is set based on the spatial accuracy)
    OrdersOfReconstruction.assign(1 + NumberOfAuxiliaryReconstructions, -1);

  } else if (strcmp(IP.Next_Control_Parameter, "Auxiliary_Reconstruction_Orders") == 0 ) {
    i_command = 0;
    ++IP.Line_Number;
    
    // Read the orders
    IP.Input_File.setf(ios::skipws);
    for (int i = 1; i< OrdersOfReconstruction.size(); ++i){
      IP.Input_File >> OrdersOfReconstruction[i];
    }

    // Read the rest of the line
    IP.Input_File.getline(buffer, sizeof(buffer));

  } else if (strcmp(IP.Next_Control_Parameter, "Erase_Auxiliary_Reconstruction_Orders") == 0 ) {
    i_command = 0;
    ++IP.Line_Number;

    // Reset all affected parameters
    OrdersOfReconstruction.assign(1,-1);
    NumberOfAuxiliaryReconstructions = 0;
    NumberOfHighOrderReconstructions = 0;
    MaxReconstructionOrder = 0;

  } else {
    i_command = INVALID_INPUT_CODE;
  } // endif

}

/*! 
 * Set properly the final HighOrder2D input parameters
 * based on user input and correlations.
 */
template<class Input_Parameters_Type>
void HighOrder2D_Input::Set_Final_Parameters(Input_Parameters_Type & IP){

  int i;

  // Create temporary array.
  // The first order should always be -1 if high-order is NOT used to advance the solution in time!
  vector<int> Temp(1,-1);

  /* Set the first reconstruction order based on the spatial accuracy required by the user
   * if the solution is calculated with high-order reconstruction.
   * This is the MAIN reconstruction order and is used for advancing the solution in time,
   * if it is not set otherwise during the simulation (e.g. hp-refinement)
   */
  if (IP.i_Reconstruction == RECONSTRUCTION_HIGH_ORDER){
    Temp[0] = OrdersOfReconstruction[0] = IP.ReconstructionOrder();
  }

  // Decide how many high-order reconstructions need to be instantiated.
  // by checking how many orders are greater than -1.
  for (i=0,NumberOfHighOrderReconstructions = 0,NumberOfAuxiliaryReconstructions = 0; i < OrdersOfReconstruction.size(); ++i){
    if (OrdersOfReconstruction[i] > -1){
      ++NumberOfHighOrderReconstructions;
      // add to the Temp if this reconstruction is an auxiliary one
      if (i != 0){
	Temp.push_back(OrdersOfReconstruction[i]);
	++NumberOfAuxiliaryReconstructions;
      }	// endif

      // determine the maximum reconstruction order
      MaxReconstructionOrder = max(MaxReconstructionOrder,OrdersOfReconstruction[i]);
    } // endif
  }// endfor

  // Set final reconstruction orders
  OrdersOfReconstruction = Temp;

}

#endif
