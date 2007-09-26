/* Input.h:  Header file defining general
             solution input parameter class. */

#ifndef _INPUT_INCLUDED
#define _INPUT_INCLUDED

/* Include the CFD_Input_Parameters class header file. */

#ifndef _CFD_INPUT_INCLUDED
#include "CFDInput.h"
#endif // _CFD_INPUT_INCLUDED

/* Define the class. */

/*!
 * Class: Input_Parameters
 *
 * @brief Input parameters for general flow solution types.
 *
 * This class defines and handles the input variables related to
 * general flow solution types.  It is a solution-state templated 
 * class based on CFD_Input_Parameters class (inherited base 
 * class).
 *
 */
template<class SOLN_pSTATE, class SOLN_cSTATE>
  class Input_Parameters: public CFD_Input_Parameters {
  private:
  public:
  //@{ @name Reference solution states:
  //! Primitive solution variables:
  SOLN_pSTATE Wo;
  //! Converved solution variables:
  SOLN_cSTATE Uo;
  //@}

  //@{ @name Constructor and desctructor
  //! Constructor (assign default values).
  Input_Parameters(void) : CFD_Input_Parameters() {
    Set_Reference_Solution_States();
  }

  //! Destructor
  ~Input_Parameters(void){ }
  //@}

  //@{ @name Other Member functions.
  //! Broadcast input parameters to all processors:
  void Broadcast(void);
  //! Parse next input line:
  int Parse_Next_Input_Control_Parameter(void);
  //! Process/parse all input from input file:
  int Process_Input_Control_Parameter_File(char *Input_File_Name_ptr,
                                           int &Command_Flag);
  //! Set reference solution states:
  void Set_Reference_Solution_States(void);
  //! Read reference solution states:
  void Read_Reference_Solution_States(istream &restart_file);
  //! Write reference solution states:
  void Write_Reference_Solution_States(ostream &restart_file);
  //@}

  //@{ @name Input-output operators:
  void Output(ostream &out_file) const;
  void Output_Problem_Type(ostream &out_file) const;
  void Output_Solver_Type(ostream &out_file) const;
  void Output_ICsBCs_Types(ostream &out_file) const;
  void Output_Solution_Type(ostream &out_file) const;
  void Output_Reference_Solution_States(ostream &out_file) const;
  void Output_IO_Types(ostream &out_file) const;
  void Output_GridAMR_Types(ostream &out_file) const;
  //@}

};

/* Input-output operators. */

template<class SOLN_pSTATE, class SOLN_cSTATE>
istream &operator >> (istream &in_file, 
                      Input_Parameters<SOLN_pSTATE, SOLN_cSTATE> &IP);

template<class SOLN_pSTATE, class SOLN_cSTATE>
ostream &operator << (ostream &out_file,
                      const Input_Parameters<SOLN_pSTATE, SOLN_cSTATE> &IP);

/********************************************************
 * Routine: Broadcast                                   *
 *                                                      *
 * Broadcast the input parameters variables to all      *
 * processors involved in the calculation from the      *
 * primary processor using the MPI broadcast routine.   *
 *                                                      *
*********************************************************/
template<class SOLN_pSTATE, class SOLN_cSTATE>
void Input_Parameters<SOLN_pSTATE, SOLN_cSTATE>::Broadcast(void) {

   /* Broadcast using base class, CFD_Input_Parameters. */

   CFD_Input_Parameters::Broadcast();

#ifdef _MPI_VERSION
   /* Update the reference solution states. */

   if (!CFFC_Primary_MPI_Processor()) {
      Set_Reference_Solution_States();
   } /* endif */
#endif

}

/********************************************************
 * Routine: Parse_Next_Input_Control_Parameter          *
 *                                                      *
 * Parses and executes the next input control parameter *
 * from the input file.                                 *
 *                                                      *
 ********************************************************/
template<class SOLN_pSTATE, class SOLN_cSTATE>
int Input_Parameters<SOLN_pSTATE, SOLN_cSTATE>::
                    Parse_Next_Input_Control_Parameter(void) {

   int i_command;

   /* Parse line using base class, CFD_Input_Parameters, protocol. */

   i_command = CFD_Input_Parameters::Parse_Next_Input_Control_Parameter();

   /* Return the parser command type indicator. */

   return (i_command);
    
}

/********************************************************
 * Routine: Process_Input_Control_Parameter_File        *
 *                                                      *
 * Reads, parses, and executes the list of input        *
 * control parameters from the standard input file.     *
 *                                                      *
 ********************************************************/
template<class SOLN_pSTATE, class SOLN_cSTATE>
int Input_Parameters<SOLN_pSTATE,SOLN_cSTATE>::
                    Process_Input_Control_Parameter_File(char *Input_File_Name_ptr,
                                                         int &Command_Flag) {

   int error_flag;
   
   /* Process input control parameter file using base class, 
      CFD_Input_Parameters, routine. */
 
   error_flag = 
     CFD_Input_Parameters::Process_Input_Control_Parameter_File(Input_File_Name_ptr,
                                                                Command_Flag);

   /* Update the reference solution states. */

   if (!error_flag) {
      Set_Reference_Solution_States();
   } /* endif */

   /* Return after processing input file. */

   return (error_flag);

}

///********************************************************
// * Routine: Set_Reference_Solution_States               *
// *                                                      *
// * Assigns default values to the reference solution     *
// * states.                                              *
// *                                                      *
// ********************************************************/
//template<class SOLN_pSTATE, class SOLN_cSTATE>
//void Input_Parameters<SOLN_pSTATE, SOLN_cSTATE>::
//                     Set_Reference_Solution_States(void) {
//
//   cerr << "\n Explicit Specialization of "
//        << "Input_Parameters::Set_Reference_Solution_States "
//        << "in Input.h requried.\n";
//   exit(1);
//
//}
//
///********************************************************
// * Routine: Read_Reference_Solution_States              *
// *                                                      *
// * Read in the reference solution states.               *
// *                                                      *
// ********************************************************/
//template<class SOLN_pSTATE, class SOLN_cSTATE>
//void Input_Parameters<SOLN_pSTATE, SOLN_cSTATE>::
//                     Read_Reference_Solution_States(istream &restart_file) {
//
//   cerr << "\n Explicit Specialization of "
//        << "Input_Parameters::Read_Reference_Solution_States "
//        << "in Input.h requried.\n";
//   exit(1);
//
//}
//
///********************************************************
// * Routine: Write_Reference_Solution_States             *
// *                                                      *
// * Write out the reference solution states.             *
// *                                                      *
// ********************************************************/
//template<class SOLN_pSTATE, class SOLN_cSTATE>
//void Input_Parameters<SOLN_pSTATE, SOLN_cSTATE>::
//                     Write_Reference_Solution_States(ostream &restart_file) {
//
//   cerr << "\n Explicit Specialization of "
//        << "Input_Parameters::Write_Reference_Solution_States "
//        << "in Input.h requried.\n";
//   exit(1);
//
//}

/*************************************************************
 * Input_Parameters -- Input-output operators.               *
 *************************************************************/
template<class SOLN_pSTATE, class SOLN_cSTATE>
ostream &operator << (ostream &out_file,
                      const Input_Parameters<SOLN_pSTATE, SOLN_cSTATE> &IP) {

   IP.Output(out_file);
   return (out_file);

}

template<class SOLN_pSTATE, class SOLN_cSTATE>
istream  &operator >> (istream &in_file,
                       Input_Parameters<SOLN_pSTATE, SOLN_cSTATE> &IP) {

   return (in_file);

}

template<class SOLN_pSTATE, class SOLN_cSTATE>
void Input_Parameters<SOLN_pSTATE, SOLN_cSTATE>::Output(ostream &out_file) const {

  out_file << setprecision(6);
  Output_Problem_Type(out_file);
  Output_Solver_Type(out_file);
  Output_ICsBCs_Types(out_file);
  Output_Solution_Type(out_file);
  Output_Reference_Solution_States(out_file);
  Output_IO_Types(out_file);
  Output_GridAMR_Types(out_file);

}

template<class SOLN_pSTATE, class SOLN_cSTATE>
void Input_Parameters<SOLN_pSTATE, SOLN_cSTATE>::Output_Problem_Type(ostream &out_file) const {

  CFD_Input_Parameters::Output_Problem_Type(out_file);

}

template<class SOLN_pSTATE, class SOLN_cSTATE>
void Input_Parameters<SOLN_pSTATE, SOLN_cSTATE>::Output_Solver_Type(ostream &out_file) const {

  CFD_Input_Parameters::Output_Solver_Type(out_file);

}

template<class SOLN_pSTATE, class SOLN_cSTATE>
void Input_Parameters<SOLN_pSTATE, SOLN_cSTATE>::Output_ICsBCs_Types(ostream &out_file) const {

  CFD_Input_Parameters::Output_ICsBCs_Types(out_file);

}

template<class SOLN_pSTATE, class SOLN_cSTATE>
void Input_Parameters<SOLN_pSTATE, SOLN_cSTATE>::Output_Solution_Type(ostream &out_file) const {

  CFD_Input_Parameters::Output_Solution_Type(out_file);

}

template<class SOLN_pSTATE, class SOLN_cSTATE>
void Input_Parameters<SOLN_pSTATE, SOLN_cSTATE>::Output_Reference_Solution_States(ostream &out_file) const {

  out_file << "\n\n Reference Solution States";
  out_file << "\n  -> Wo: "
           << Wo;
  out_file << "\n  -> Uo: "
           << Uo;

}

template<class SOLN_pSTATE, class SOLN_cSTATE>
void Input_Parameters<SOLN_pSTATE, SOLN_cSTATE>::Output_IO_Types(ostream &out_file) const {

  CFD_Input_Parameters::Output_IO_Types(out_file);

}


template<class SOLN_pSTATE, class SOLN_cSTATE>
void Input_Parameters<SOLN_pSTATE, SOLN_cSTATE>::Output_GridAMR_Types(ostream &out_file) const {

  CFD_Input_Parameters::Output_GridAMR_Types(out_file);

}

#endif // _INPUT_INCLUDED
