/* FASMultigridInput.h:  Header file declaring multigrid input class. */

#ifndef _FASMULTIGRIDINPUT_INCLUDED
#define _FASMULTIGRIDINPUT_INCLUDED

// Include required C++ libraries.

#include <cstdio>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cassert>
#include <cstdlib>
#include <cstring>

using namespace std;

// Include CFD and MPI header files.

#ifndef _CFD_INCLUDED
#include "../CFD/CFD.h"
#endif // _CFD_INCLUDED

#ifndef _MPI_INCLUDED
#include "../MPI/MPI.h"
#endif // _MPI_INCLUDED

#define MULTIGRID_INPUT_PARAMETER_LENGTH 256

/*!
 * Class: Multigrid_Input_Parameters
 *
 * @brief Input Parameters for Multigrid.
 *
 * This class defines and handles the input variables related to the
 * Full Approximation Storage (FAS) multigrid method for solving 
 * nonlinear partial differential equations.  Note that unsteady
 * calculations can be performed using dual-time-stepping.
 *
 */
class Multigrid_Input_Parameters{
 private:
 public:
  //@{ @name Multigrid cycle definitions
  int Levels; //!< Number of levels in multigrid cycle.
  char Cycle_Type[MULTIGRID_INPUT_PARAMETER_LENGTH];
  int i_Cycle;
  int N_Stage;
  int Number_of_Cycles_per_Stage_for_Full_Multigrid;
  int Defect_Correction; //!< Flag for fine to coarse defect correction.
  int Prolong_Using_Injection;
  int First_Order_Coarse_Mesh_Reconstruction;
  int Apply_Coarse_Mesh_Boundary_Conditions;
  int Injection_at_Dirichlet_Boundary_Conditions;
  int Update_Stability_Switch;
  int Maximum_Number_of_Update_Reductions;
  //@}

  //@{ @name Multigrid smoother type and parameters
  char Smoothing_Type[MULTIGRID_INPUT_PARAMETER_LENGTH];
  int i_Smoothing;
  int Number_of_Smooths_on_Finest_Level;
  int Number_of_Pre_Smooths;
  int Number_of_Post_Smooths;
  int Number_of_Smooths_on_Coarsest_Level;
  double Absolute_Convergence_Tolerance, Relative_Convergence_Tolerance, 
         FMG_Absolute_Convergence_Tolerance, FMG_Relative_Convergence_Tolerance;
  //@}

  //@{ @name Implicit residual smoothing control parameters:
  int Residual_Smoothing;
  double Residual_Smoothing_Epsilon;
  int Residual_Smoothing_Gauss_Seidel_Iterations;
  //@}

  //@{ @name Dual-time-stepping parameters
  int i_Dual_Time_Stepping; //!< Dual-time-stepping flag (on or off).
  int Ncycles_Regular_Multigrid;
  int Ncycles_Full_Multigrid;
  char Physical_Time_Integration_Type[MULTIGRID_INPUT_PARAMETER_LENGTH];
  int i_Physical_Time_Integration;
  double Physical_Time_CFL_Number;
  double Dual_Time_Convergence_Residual_Level;
  int i_Dual_Time_Preconditioning;
  //@}

  //@{ @name Output control parameters
  int Write_Output_Cells_Frequency; // set to zero to turn off
  //@}

  //@{ @name Constructor and desctructor
  //! Constructor (assign default values).
  Multigrid_Input_Parameters(void) {
    // Multigrid cycle parameters:
    strcpy(Cycle_Type,"V");
    i_Cycle = MULTIGRID_V_CYCLE;
    Levels = 2;
    Number_of_Cycles_per_Stage_for_Full_Multigrid = 0;
    Defect_Correction = ON;
    Prolong_Using_Injection = OFF;
    First_Order_Coarse_Mesh_Reconstruction = ON;
    Apply_Coarse_Mesh_Boundary_Conditions = ON;
    Injection_at_Dirichlet_Boundary_Conditions = ON;
    Update_Stability_Switch = ON;
    Maximum_Number_of_Update_Reductions = 5;
    // Multigrid smoother parameters:
    strcpy(Smoothing_Type,"Multistage_Optimal_Smoothing");
    i_Smoothing = TIME_STEPPING_MULTISTAGE_OPTIMAL_SMOOTHING;
    N_Stage = 1;
    Number_of_Smooths_on_Finest_Level = 1;
    Number_of_Pre_Smooths = 1;
    Number_of_Post_Smooths = 0;
    Number_of_Smooths_on_Coarsest_Level = 3;
    Absolute_Convergence_Tolerance  = -1.0;
    Relative_Convergence_Tolerance  = -1.0;
    FMG_Absolute_Convergence_Tolerance = -1.0;
    FMG_Relative_Convergence_Tolerance = -1.0;
    //Implicit residual smoothing control parameters:   
    Residual_Smoothing = 0;
    Residual_Smoothing_Epsilon = ZERO;
    Residual_Smoothing_Gauss_Seidel_Iterations = 2;
    // Dual-time-stepping parameters:
    i_Dual_Time_Stepping = OFF;
    Ncycles_Regular_Multigrid = 0;
    Ncycles_Full_Multigrid = 0;
    strcpy(Physical_Time_Integration_Type,"Implicit Euler");
    i_Physical_Time_Integration = TIME_STEPPING_IMPLICIT_EULER;
    Physical_Time_CFL_Number = 0.50;
    Dual_Time_Convergence_Residual_Level = -1000.0;
    i_Dual_Time_Preconditioning = SCALAR_LOCAL_TIME_STEPPING;
    Write_Output_Cells_Frequency = 0;
  }

  //! Destructor
  ~Multigrid_Input_Parameters(void){}
  //@}

  //@{ @name Other Member functions.
  //! Broadcast input parameters to all processors:
  void Broadcast(void);
  //! Parse next input line:
  int Parse_Next_Input_Control_Parameter(char *code, stringstream &value);
  //! Check validity of specified input parameters:
  int Check_Inputs(void);
  //@}

  //@{ @name Input-output operators:
  friend ostream &operator << (ostream &out_file,
		               Multigrid_Input_Parameters &IP);
  friend istream &operator >> (istream &in_file,
			       Multigrid_Input_Parameters &IP);
  void Output(ostream &out_file) const;
  //@}

};

#endif // _FASMULTIGRIDINPUT_INCLUDED
