/**********************************************************************
 * FASMultigrid2DInput.h: Header file declaring FAS Multigrid related *
 *                        inputs and functions.                       *
 **********************************************************************/

#ifndef _FASMULTIGRID2DINPUT_INCLUDED
#define _FASMULTIGRID2DINPUT_INCLUDED

// Include required C++ libraries.

#include <cstdio>
#include <iostream>
#include <iomanip>
#include <fstream>
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

#define INPUT_PARAMETER_LENGTH_MULTIGRID 128

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
  char Cycle_Type[INPUT_PARAMETER_LENGTH_MULTIGRID];
  int i_Cycle;
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
  char Smoothing_Type[INPUT_PARAMETER_LENGTH_MULTIGRID];
  int i_Smoothing;
  int Number_of_Smooths_on_Finest_Level;
  int Number_of_Pre_Smooths;
  int Number_of_Post_Smooths;
  int Number_of_Smooths_on_Coarsest_Level;
  double Absolute_Convergence_Tolerance, Relative_Convergence_Tolerance, 
         FMG_Absolute_Convergence_Tolerance, FMG_Relative_Convergence_Tolerance;
  //@}

  //@{ @name Dual-time-stepping parameters
  int i_Dual_Time_Stepping; //!< Dual-time-stepping flag (on or off).
  int Ncycles_Regular_Multigrid;
  int Ncycles_Full_Multigrid;
  char Physical_Time_Integration_Type[INPUT_PARAMETER_LENGTH_MULTIGRID];
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
    Number_of_Smooths_on_Finest_Level = 1;
    Number_of_Pre_Smooths = 1;
    Number_of_Post_Smooths = 0;
    Number_of_Smooths_on_Coarsest_Level = 3;
    Absolute_Convergence_Tolerance  = -1.0;
    Relative_Convergence_Tolerance  = -1.0;
    FMG_Absolute_Convergence_Tolerance = -1.0;
    FMG_Relative_Convergence_Tolerance = -1.0;
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
  ~Multigrid_Input_Parameters(void) {}
  //@}

  //@{ @name Broadcast functions.

  //! Broadcast input parameters to all processors.
  void Broadcast_Input_Parameters(void) {
#ifdef _MPI_VERSION
    // Multigrid cycle parameters:
    MPI::COMM_WORLD.Bcast(&(Levels),
			  1,
			  MPI::INT,0);
    MPI::COMM_WORLD.Bcast(Cycle_Type,
			  INPUT_PARAMETER_LENGTH_MULTIGRID,
			  MPI::CHAR,0);
    MPI::COMM_WORLD.Bcast(&(i_Cycle),
			  1,
			  MPI::INT,0);
    MPI::COMM_WORLD.Bcast(&(Number_of_Cycles_per_Stage_for_Full_Multigrid),
			  1,
			  MPI::INT,0);
    MPI::COMM_WORLD.Bcast(&(Defect_Correction),
			  1,
			  MPI::INT,0);
    MPI::COMM_WORLD.Bcast(&(Prolong_Using_Injection),
			  1,
			  MPI::INT,0);
    MPI::COMM_WORLD.Bcast(&(First_Order_Coarse_Mesh_Reconstruction),
			  1,
			  MPI::INT,0);
    MPI::COMM_WORLD.Bcast(&(Apply_Coarse_Mesh_Boundary_Conditions),
			  1,
			  MPI::INT,0);
    MPI::COMM_WORLD.Bcast(&(Injection_at_Dirichlet_Boundary_Conditions),
			  1,
			  MPI::INT,0);
    MPI::COMM_WORLD.Bcast(&(Update_Stability_Switch),
			  1,
			  MPI::INT,0);
    MPI::COMM_WORLD.Bcast(&(Maximum_Number_of_Update_Reductions),
			  1,
			  MPI::INT,0);
    // Multigrid smoother parameters:
    MPI::COMM_WORLD.Bcast(Smoothing_Type,
			  INPUT_PARAMETER_LENGTH_MULTIGRID,
			  MPI::CHAR,0);
    MPI::COMM_WORLD.Bcast(&(i_Smoothing),
			  1,
			  MPI::INT,0);
    MPI::COMM_WORLD.Bcast(&(Number_of_Smooths_on_Finest_Level),
			  1,
			  MPI::INT,0);
    MPI::COMM_WORLD.Bcast(&(Number_of_Pre_Smooths),
			  1,
			  MPI::INT,0);
    MPI::COMM_WORLD.Bcast(&(Number_of_Post_Smooths),
			  1,
			  MPI::INT,0);
    MPI::COMM_WORLD.Bcast(&(Number_of_Smooths_on_Coarsest_Level),
			  1,
			  MPI::INT,0);
    MPI::COMM_WORLD.Bcast(&(Absolute_Convergence_Tolerance),
			  1,
			  MPI::DOUBLE,0);
    MPI::COMM_WORLD.Bcast(&(Relative_Convergence_Tolerance),
			  1,
			  MPI::DOUBLE,0);
    MPI::COMM_WORLD.Bcast(&(FMG_Absolute_Convergence_Tolerance),
			  1,
			  MPI::DOUBLE,0);
    MPI::COMM_WORLD.Bcast(&(FMG_Relative_Convergence_Tolerance),
			  1,
			  MPI::DOUBLE,0);
    // Dual-time-stepping parameters:
    MPI::COMM_WORLD.Bcast(&(i_Dual_Time_Stepping),
			  1,
			  MPI::INT,0);
    MPI::COMM_WORLD.Bcast(&(Ncycles_Regular_Multigrid),
			  1,
			  MPI::INT,0);
    MPI::COMM_WORLD.Bcast(&(Ncycles_Full_Multigrid),
			  1,
			  MPI::INT,0);
    MPI::COMM_WORLD.Bcast(Physical_Time_Integration_Type,
			  INPUT_PARAMETER_LENGTH_MULTIGRID,
			  MPI::CHAR,0);
    MPI::COMM_WORLD.Bcast(&(i_Physical_Time_Integration),
			  1,
			  MPI::INT,0);
    MPI::COMM_WORLD.Bcast(&(Physical_Time_CFL_Number),
			  1,
			  MPI::DOUBLE,0);
    MPI::COMM_WORLD.Bcast(&(Dual_Time_Convergence_Residual_Level),
			  1,
			  MPI::DOUBLE,0);
    MPI::COMM_WORLD.Bcast(&(i_Dual_Time_Preconditioning),
			  1,
			  MPI::INT,0);
    MPI::COMM_WORLD.Bcast(&(Write_Output_Cells_Frequency),
			  1,
			  MPI::INT,0);
#endif
  }

  //! Broadcast input parameters to all processors associated with the
  //! specified communicator.
#ifdef _MPI_VERSION
  void Broadcast_Input_Parameters(MPI::Intracomm &Communicator,
				  const int Source_CPU) {
    // Multigrid cycle parameters:
    Communicator.Bcast(&(Levels),
		       1,
		       MPI::INT,Source_CPU);
    Communicator.Bcast(Cycle_Type,
                       INPUT_PARAMETER_LENGTH_MULTIGRID,
                       MPI::CHAR,Source_CPU);
    Communicator.Bcast(&(i_Cycle),
		       1,
		       MPI::INT,Source_CPU);
    Communicator.Bcast(&(Number_of_Cycles_per_Stage_for_Full_Multigrid),
		       1,
		       MPI::INT,Source_CPU);
    Communicator.Bcast(&(Defect_Correction),
		       1,
		       MPI::INT,Source_CPU);
    Communicator.Bcast(&(Prolong_Using_Injection),
		       1,
		       MPI::INT,Source_CPU);
    Communicator.Bcast(&(First_Order_Coarse_Mesh_Reconstruction),
		       1,
		       MPI::INT,Source_CPU);
    Communicator.Bcast(&(Apply_Coarse_Mesh_Boundary_Conditions),
		       1,
		       MPI::INT,Source_CPU);
    Communicator.Bcast(&(Injection_at_Dirichlet_Boundary_Conditions),
		       1,
		       MPI::INT,Source_CPU);
    Communicator.Bcast(&(Update_Stability_Switch),
		       1,
		       MPI::INT,Source_CPU);
    Communicator.Bcast(&(Maximum_Number_of_Update_Reductions),
		       1,
		       MPI::INT,Source_CPU);
    // Multigrid smoother parameters:
    Communicator.Bcast(Smoothing_Type,
                       INPUT_PARAMETER_LENGTH_MULTIGRID,
                       MPI::CHAR,Source_CPU);
    Communicator.Bcast(&(i_Smoothing),
		       1,
		       MPI::INT,Source_CPU);
    Communicator.Bcast(&(Number_of_Smooths_on_Finest_Level),
		       1,
		       MPI::INT,Source_CPU);
    Communicator.Bcast(&(Number_of_Pre_Smooths),
		       1,
		       MPI::INT,Source_CPU);
    Communicator.Bcast(&(Number_of_Post_Smooths),
		       1,
		       MPI::INT,Source_CPU);
    Communicator.Bcast(&(Number_of_Smooths_on_Coarsest_Level),
		       1,
		       MPI::INT,Source_CPU);
    Communicator.Bcast(&(Absolute_Convergence_Tolerance),
		       1,
		       MPI::DOUBLE,Source_CPU);
    Communicator.Bcast(&(Relative_Convergence_Tolerance),
		       1,
		       MPI::DOUBLE,Source_CPU);
    Communicator.Bcast(&(FMG_Absolute_Convergence_Tolerance),
		       1,
		       MPI::DOUBLE,Source_CPU);
    Communicator.Bcast(&(FMG_Relative_Convergence_Tolerance),
		       1,
		       MPI::DOUBLE,Source_CPU);
    // Dual-time-stepping parameters:
    Communicator.Bcast(&(i_Dual_Time_Stepping),
		       1,
		       MPI::INT,Source_CPU);
    Communicator.Bcast(&(Ncycles_Regular_Multigrid),
		       1,
		       MPI::INT,Source_CPU);
    Communicator.Bcast(&(Ncycles_Full_Multigrid),
		       1,
		       MPI::INT,Source_CPU);
    Communicator.Bcast(Physical_Time_Integration_Type,
                       INPUT_PARAMETER_LENGTH_MULTIGRID,
                       MPI::CHAR,Source_CPU);
    Communicator.Bcast(&(i_Physical_Time_Integration),
		       1,
		       MPI::INT,Source_CPU);
    Communicator.Bcast(&(Physical_Time_CFL_Number),
		       1,
		       MPI::DOUBLE,Source_CPU);
    Communicator.Bcast(&(Dual_Time_Convergence_Residual_Level),
		       1,
		       MPI::DOUBLE,Source_CPU);
    Communicator.Bcast(&(i_Dual_Time_Preconditioning),
		       1,
		       MPI::INT,Source_CPU);
    Communicator.Bcast(&(Write_Output_Cells_Frequency),
		       1,
		       MPI::INT,Source_CPU);
  }
#endif

  //@{ @name Input-output operators
  friend ostream &operator << (ostream &out_file,
		               const Multigrid_Input_Parameters &IP);
  friend istream &operator >> (istream &in_file,
			       Multigrid_Input_Parameters &IP);
  //@}

};

/**********************************************************************
 * Multigrid_Input_Parameters -- Input-output operators.              *
 **********************************************************************/
inline ostream &operator << (ostream &out_file,
			     const Multigrid_Input_Parameters &IP) {
  // Multigrid cycle parameters:
  out_file << "\n  -> Number of Levels in Multigrid: "
	   << IP.Levels;
  out_file << "\n  -> Multigrid Cycle Type: "
	   << IP.Cycle_Type;
  if (IP.Number_of_Cycles_per_Stage_for_Full_Multigrid > 0) {
    out_file << "\n  -> Full Multigrid Startup: "
	     << IP.Number_of_Cycles_per_Stage_for_Full_Multigrid
	     << " cycles per stage";
  }
  out_file << "\n  -> Defect correction: "
	   << IP.Defect_Correction;
  out_file << "\n  -> Prolong using injection: "
	   << IP.Prolong_Using_Injection;
  out_file << "\n  -> First Order Coarse Mesh Reconstruction: "
	   << IP.First_Order_Coarse_Mesh_Reconstruction;
  out_file << "\n  -> Apply Coarse Mesh Boundary Conditions: "
	   << IP.Apply_Coarse_Mesh_Boundary_Conditions;
  out_file << "\n  -> Injection at Dirichlet Boundary Conditions: "
	   << IP.Injection_at_Dirichlet_Boundary_Conditions;
  out_file << "\n  -> Update Stability Switch: "
	   << IP.Update_Stability_Switch;
  out_file << "\n  -> Maximum Number of Update Reductions: "
	   << IP.Maximum_Number_of_Update_Reductions;
  // Multigrid smoother parameters:
  out_file << "\n  -> Multigrid Smoothing Type: "
	   << IP.Smoothing_Type;
  out_file << "\n  -> Number of Smoothing Iterations on the Finest Multigrid Level:  "
	   << IP.Number_of_Smooths_on_Finest_Level;
  out_file << "\n  -> Number of Multigrid Pre-Smoothing Iterations: "
	   << IP.Number_of_Pre_Smooths;
  out_file << "\n  -> Number of Multigrid Post-Smoothing Iterations: "
	   << IP.Number_of_Post_Smooths;
  out_file << "\n  -> Number of Smoothing Iterations on the Coarsest Multigrid Level:  "
	   << IP.Number_of_Smooths_on_Coarsest_Level;
  if (IP.Absolute_Convergence_Tolerance > 0) {
    out_file.setf(ios::scientific);
    out_file << "\n  -> Absolute Tolerance for Regular Multigrid: "
	     << IP.Absolute_Convergence_Tolerance;
    out_file.unsetf(ios::scientific);
  }
  if (IP.Relative_Convergence_Tolerance > 0) {
    out_file.setf(ios::scientific);
    out_file << "\n  -> Relative Tolerance for Regular Multigrid: "
	     << IP.Relative_Convergence_Tolerance;
    out_file.unsetf(ios::scientific);
  }
  if (IP.Number_of_Cycles_per_Stage_for_Full_Multigrid > 0 && IP.FMG_Absolute_Convergence_Tolerance > 0) {
    out_file.setf(ios::scientific);
    out_file << "\n  -> Absolute Tolerance for Full Multigrid: "
             << IP.FMG_Absolute_Convergence_Tolerance;
    out_file.unsetf(ios::scientific);
  }
  if (IP.Number_of_Cycles_per_Stage_for_Full_Multigrid > 0 && IP.FMG_Relative_Convergence_Tolerance > 0) {
    out_file.setf(ios::scientific);
    out_file << "\n  -> Relative Tolerance for Full Multigrid: "
	     << IP.FMG_Relative_Convergence_Tolerance;
    out_file.unsetf(ios::scientific);
  }
  // Dual-time-stepping parameters:
  if (IP.i_Dual_Time_Stepping) {
    out_file << "\n  -> Dual-time-stepping time integration type: "
	     << IP.Physical_Time_Integration_Type;
    out_file << "\n  -> Number of regular multigrid cycles: "
	     << IP.Ncycles_Regular_Multigrid;
    out_file << "\n  -> Number of full multigrid cycles: "
	     << IP.Ncycles_Full_Multigrid;
    out_file << "\n  -> Physical-time-stepping CFL number: "
	     << IP.Physical_Time_CFL_Number;
    out_file << "\n  -> Dual-time-stepping convergence residual level: "
	     << IP.Dual_Time_Convergence_Residual_Level;
    out_file << "\n  -> Dual-time-step type: "
	     << IP.i_Dual_Time_Preconditioning;
  }
  if (IP.Write_Output_Cells_Frequency > 0) {
    out_file << "\n  -> Write Output Cells Frequency for Multigrid: "
	     << IP.Write_Output_Cells_Frequency;
  }
  // Output successful.
  return out_file;
}

inline istream &operator >> (istream &in_file,
			     Multigrid_Input_Parameters &IP) {
  return in_file;
}

/**********************************************************************
 * Routine: Check_Input_Parameters                                    *
 *                                                                    *
 * Check Input Parameters as stored in for consistency.               *
 *                                                                    *
 **********************************************************************/
template <class Quad_Soln_Input_Parameters>
int Check_Input_Parameters(Quad_Soln_Input_Parameters &IP) {

  // Multigrid Related Input Checking.
  if (IP.i_Time_Integration == TIME_STEPPING_MULTIGRID) {

    // Make sure the smoothing method is not Multigrid or NKS.
    if (IP.Multigrid_IP.i_Smoothing == TIME_STEPPING_MULTIGRID ||
	IP.Multigrid_IP.i_Smoothing == TIME_STEPPING_NKS) {
      cout << "Check_Input_Parameters: invalid smoothing method" << endl;
      cout.flush();
      return 1;
    }

    // Make sure the number of cylce levels is greater than two if full
    // multigrid has been selected.  Otherwise set the number of full
    // multigrid cylces to zero.
    if (IP.Multigrid_IP.Number_of_Cycles_per_Stage_for_Full_Multigrid > 0 &&
	IP.Multigrid_IP.Levels < 3) {
      IP.Multigrid_IP.Number_of_Cycles_per_Stage_for_Full_Multigrid = 0;
    }

    // Make sure the time accurate flag is set appropriately.
    if (!IP.Multigrid_IP.i_Dual_Time_Stepping) IP.Time_Accurate = 0;
    else IP.Time_Accurate = 1;

  }

  // Input parameters are consistent.  Exit successfully.
  return 0;

}

#endif // _FASMULTIGRID2DINPUT_INCLUDED
