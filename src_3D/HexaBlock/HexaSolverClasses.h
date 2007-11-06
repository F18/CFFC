#ifndef _HEXA_SOLVER_CLASSES_INCLUDED
#define _HEXA_SOLVER_CLASSES_INCLUDED

/*! **************************************************************
 * class:  HexaSolver_Data                                       *
 *                                                               *
 * This class contains the no solution specific data, ie grids,  *
 * block lists, etc. required by the Hexa Solver.  Its primary   *
 * function is to serve as a data container to organize and      *
 * simplify data handling within the Hexa Solver functions.      *
 *                                                               *
 * Data Members:                                                 *
 *    Initial_Mesh               -                               *
 *    Local_Adaptive_Block_List  -                               *
 *    Global_Adaptive_Block_List -                               *
 *    Octree                     -                               *
 *                                                               *
 *****************************************************************/ 
class HexaSolver_Data {
  private:
  protected: 
  public:

  Grid3D_Hexa_Multi_Block             Initial_Mesh;  
  AdaptiveBlock3D_List                Local_Adaptive_Block_List; 
  AdaptiveBlock3D_ResourceList        Global_Adaptive_Block_List;
  Octree_DataStructure                Octree;
  
  // cpu time variables
  CPUTime processor_cpu_time, total_cpu_time;  
  
  // Time step / iteration counters
  int number_of_explicit_time_steps;       
  int number_of_implicit_time_steps;

  // Physical Time for Unsteady Calculations
  double Time;
  // Batch mode flag, deactivates STDOUT ie. cout
  int batch_flag;
  // Residaul File Output Stream
  ofstream residual_file;

  int total_number_of_time_steps() { return (number_of_explicit_time_steps+number_of_implicit_time_steps); }

  /*********************************/
  //default constructor(s)
  HexaSolver_Data(void){}
  HexaSolver_Data(int batch): batch_flag(batch) {}

  // Destructor
  ~HexaSolver_Data() {}
 
};

/*! **************************************************************
 * class:  HexaSolver_Solution_Data                              *
 *                                                               *
 * This class contains the solution specific data (based on      *
 * pState & cState) required by the Hexa Solver.  Its primary    *
 * function is to serve as a data container to organize and      *
 * simplify data handling within the Hexa Solver functions.      *
 *                                                               *
 * Data Members:                                                 *
 *    Input                  - Input Parameters                  *
 *    Local_Solution_Blocks  - Local Solution Data               *
 *    command_flag           - Solver Directive flag             *
 *                                                               *
 *****************************************************************/ 
template<typename SOLN_pSTATE, typename SOLN_cSTATE>
class HexaSolver_Solution_Data {
  private:
  protected: 
  public:

  Input_Parameters<SOLN_pSTATE, SOLN_cSTATE>               Input;
  Hexa_Multi_Block<Hexa_Block<SOLN_pSTATE, SOLN_cSTATE> >  Local_Solution_Blocks; 

  // Control Flag 
  int command_flag;

  // Constructor
  HexaSolver_Solution_Data(void) {} 
  // Destructor
  ~HexaSolver_Solution_Data() {}

  int Get_Input_Parameters(char *Input_File_Name, int batch);

};


/*! **************************************************************
 * class:  HexaSolver_Solution_Data                              *
 *****************************************************************/ 
template<typename SOLN_pSTATE, typename SOLN_cSTATE>
int HexaSolver_Solution_Data<SOLN_pSTATE,SOLN_cSTATE>::
Get_Input_Parameters(char *Input_File_Name,int batch_flag){

  int error_flag(0);

  // The primary MPI processor processes the input parameter file.
  if (CFFC_Primary_MPI_Processor()) {
    if (!batch_flag) {
      cout << "\n Reading input data file `"<< Input_File_Name << "'.";
      cout.flush();
    }    
    // Reads Input & sets "command_flag"
    error_flag = Input.Process_Input_Control_Parameter_File(Input_File_Name,
							    command_flag);
  }
 
  // Broadcast input solution parameters to other MPI processors.
  CFFC_Barrier_MPI(); 
  CFFC_Broadcast_MPI(&error_flag, 1);
  if (error_flag != 0) return (error_flag);

  CFFC_Broadcast_MPI(&command_flag, 1);
  if (command_flag == TERMINATE_CODE) return error_flag;

  Input.Broadcast();
  
  return error_flag;
}



#endif //_HEXA_SOLVER_CLASSES_INCLUDED

