/* FASMultigridInput.cc: Definition of Multigrid_Input_Parameters class member functions. */

/* Include the FASMultigridInput header file. */

#ifndef _FASMULTIGRIDINPUT_INCLUDED
#include "FASMultigridInput.h"
#endif // _FASMULTIGRIDINPUT_INCLUDED

/* Define member functions. */

/***************************************************************************
 * Multigrid_Input_Parameters::Broadcast -- Broadcast to all processors.   *
 ***************************************************************************/
void Multigrid_Input_Parameters::Broadcast(void) {

#ifdef _MPI_VERSION
    // Multigrid cycle parameters:
    MPI::COMM_WORLD.Bcast(&(Levels),
			  1,
			  MPI::INT,0);
    MPI::COMM_WORLD.Bcast(Cycle_Type,
			  MULTIGRID_INPUT_PARAMETER_LENGTH,
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
			  MULTIGRID_INPUT_PARAMETER_LENGTH,
			  MPI::CHAR,0);
    MPI::COMM_WORLD.Bcast(&(i_Smoothing),
			  1,
			  MPI::INT,0);
    MPI::COMM_WORLD.Bcast(&(N_Stage),
                          1,
                          MPI::INT, 0);
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

    //Implicit residual smoothing control parameters:   
    MPI::COMM_WORLD.Bcast(&(Residual_Smoothing),
			  1,
			  MPI::INT,0);
    MPI::COMM_WORLD.Bcast(&(Residual_Smoothing_Epsilon),
			  1,
			  MPI::DOUBLE,0);
    MPI::COMM_WORLD.Bcast(&(Residual_Smoothing_Gauss_Seidel_Iterations),
			  1,
			  MPI::INT,0);

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
			  MULTIGRID_INPUT_PARAMETER_LENGTH,
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

    // Output control parameters:
    MPI::COMM_WORLD.Bcast(&(Write_Output_Cells_Frequency),
			  1,
			  MPI::INT,0);
#endif

}


/*********************************************************************************
 * Multigrid_Input_Parameters::Parse_Next_Input_Control_Parameter - Parse input. *
 *********************************************************************************/
int Multigrid_Input_Parameters::Parse_Next_Input_Control_Parameter(char *code,
                                                                   stringstream &value) {

// Returns:
//  - INVALID_INPUT_VALUE if code is valid but value is invalid
//  - INVALID_INPUT_CODE  if unknown code

  int i_command = INVALID_INPUT_CODE;
  string value_string;

  if (strcmp(code,"Multigrid_Levels") == 0) {
    i_command = 2001;
    value >> Levels;
    if (Levels <= ONE) i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(code,"Multigrid_Cycle_Type") == 0) {
    i_command = 2002;
    value >> value_string;
    strcpy(Cycle_Type, value_string.c_str());
    if (strcmp(Cycle_Type,"V") == 0 ||
  	strcmp(Cycle_Type,"v") == 0) {
      i_Cycle = MULTIGRID_V_CYCLE;
    } else if (strcmp(Cycle_Type,"W") == 0 ||
  	       strcmp(Cycle_Type,"w") == 0) {
      i_Cycle = MULTIGRID_W_CYCLE;
    } else {
      i_command = INVALID_INPUT_VALUE;
    }

  } else if (strcmp(code,"Full_Multigrid") == 0) {
    i_command = 2003;
    value >> Number_of_Cycles_per_Stage_for_Full_Multigrid;
    if (Number_of_Cycles_per_Stage_for_Full_Multigrid < 0) Number_of_Cycles_per_Stage_for_Full_Multigrid = 0;

  } else if (strcmp(code,"Defect_Correction") == 0) {
    i_command = 2004;
    value >> value_string;
    if (value_string == "ON") {
      Defect_Correction = ON;
    } else if (value_string == "OFF") {
      Defect_Correction = OFF;
    } else {
      i_command = INVALID_INPUT_VALUE;
    }

  } else if (strcmp(code,"First_Order_Coarse_Mesh_Reconstruction") == 0) {
    i_command = 2005;
    value >> value_string;
    if (value_string == "ON") {
      First_Order_Coarse_Mesh_Reconstruction = ON;
    } else if (value_string == "OFF") {
      First_Order_Coarse_Mesh_Reconstruction = OFF;
    } else {
      i_command = INVALID_INPUT_VALUE;
    }

  } else if (strcmp(code,"Prolong_Using_Injection") == 0) {
    i_command = 2006;
    value >> value_string;
    if (value_string == "ON") {
      Prolong_Using_Injection = ON;
    } else if (value_string == "OFF") {
      Prolong_Using_Injection = OFF;
    } else {
      i_command = INVALID_INPUT_VALUE;
    }

  } else if (strcmp(code,"Apply_Coarse_Mesh_Boundary_Conditions") == 0) {
    i_command = 2007;
    value >> value_string;
    if (value_string == "ON") {
      Apply_Coarse_Mesh_Boundary_Conditions = ON;
    } else if (value_string == "OFF") {
      Apply_Coarse_Mesh_Boundary_Conditions = OFF;
    } else {
      i_command = INVALID_INPUT_VALUE;
    }

  } else if (strcmp(code,"Injection_at_Dirichlet_Boundary_Conditions") == 0) {
    i_command = 2008;
    value >> value_string;
    if (value_string == "ON") {
      Injection_at_Dirichlet_Boundary_Conditions = ON;
    } else if (value_string == "OFF") {
      Injection_at_Dirichlet_Boundary_Conditions = OFF;
    } else {
      i_command = INVALID_INPUT_VALUE;
    }

  } else if (strcmp(code,"Update_Stability_Switch") == 0) {
    i_command = 2009;
    value >> value_string;
    if (value_string == "ON") {
      Update_Stability_Switch = ON;
    } else if (value_string == "OFF") {
      Update_Stability_Switch = OFF;
    } else {
      i_command = INVALID_INPUT_VALUE;
    }

  } else if (strcmp(code,"Maximum_Number_of_Update_Reductions") == 0) {
    i_command = 2010;
    value >> Maximum_Number_of_Update_Reductions;
    if (Maximum_Number_of_Update_Reductions < 1) i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(code,"Multigrid_Number_of_Smooths_on_Finest_Level") == 0) {
    i_command = 2011;
    value >> Number_of_Smooths_on_Finest_Level;
    if (Number_of_Smooths_on_Finest_Level < 0) i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(code,"Multigrid_Number_of_Pre_Smooths") == 0) {
    i_command = 2012;
    value >> Number_of_Pre_Smooths;
    if (Number_of_Pre_Smooths < 0) i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(code,"Multigrid_Number_of_Post_Smooths") == 0) {
    i_command = 2013;
    value >> Number_of_Post_Smooths;
    if (Number_of_Post_Smooths < 0) i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(code,"Multigrid_Number_of_Smooths_on_Coarsest_Level") == 0) {
    i_command = 2014;
    value >> Number_of_Smooths_on_Coarsest_Level;
    if (Number_of_Smooths_on_Coarsest_Level < 0) i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(code, "Multigrid_Absolute_Convergence_Tolerance") == 0) {
    i_command = 2015;
    value >> Absolute_Convergence_Tolerance;

  } else if (strcmp(code, "Multigrid_Relative_Convergence_Tolerance") == 0) {
    i_command = 2016;
    value >> Relative_Convergence_Tolerance;

  } else if (strcmp(code, "FMG_Absolute_Convergence_Tolerance") == 0) {
    i_command = 2017;
    value >> FMG_Absolute_Convergence_Tolerance;

  } else if (strcmp(code, "FMG_Relative_Convergence_Tolerance") == 0) {
    i_command = 2018;
    value >> FMG_Relative_Convergence_Tolerance;

  } else if (strcmp(code,"Multigrid_Smoothing_Type") == 0) {
    i_command = 2019;
    value >> value_string;
    strcpy(Smoothing_Type, value_string.c_str());
    if (strcmp(Smoothing_Type,"Explicit_Euler") == 0) {
      i_Smoothing = TIME_STEPPING_EXPLICIT_EULER;
      N_Stage = 1;
    } else if (strcmp(Smoothing_Type,"Explicit_Predictor_Corrector") == 0) {
      i_Smoothing = TIME_STEPPING_EXPLICIT_PREDICTOR_CORRECTOR;
      N_Stage = 2;
    } else if (strcmp(Smoothing_Type,"Explicit_Runge_Kutta") == 0) {
      i_Smoothing = TIME_STEPPING_EXPLICIT_RUNGE_KUTTA;
      N_Stage = 5;
    } else if (strcmp(Smoothing_Type,"Multistage_Optimal_Smoothing") == 0) {
      i_Smoothing = TIME_STEPPING_MULTISTAGE_OPTIMAL_SMOOTHING;
      N_Stage = 4;
    } else {
      i_command = INVALID_INPUT_VALUE;
    }

  } else if (strcmp(code,"Ncycles_Regular_Multigrid") == 0) {
    i_command = 2020;
    value >> Ncycles_Regular_Multigrid;
    if (Ncycles_Regular_Multigrid < 0) i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(code,"Ncycles_Full_Multigrid") == 0) {
    i_command = 2021;
    value >> Ncycles_Full_Multigrid;
    if (Ncycles_Full_Multigrid < 0) i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(code,"Physical_Time_Integration_Type") == 0) {
    i_command = 2022;
    value >> value_string;
    strcpy(Physical_Time_Integration_Type, value_string.c_str());
    if (strcmp(Physical_Time_Integration_Type,"Implicit_Euler") == 0) {
      i_Physical_Time_Integration = TIME_STEPPING_IMPLICIT_EULER;
    } else if (strcmp(Physical_Time_Integration_Type,"Implicit_Trapezoidal") == 0) {
      i_Physical_Time_Integration = TIME_STEPPING_IMPLICIT_TRAPEZOIDAL;
    } else if (strcmp(Physical_Time_Integration_Type,"Second_Order_Backwards") == 0) {
      i_Physical_Time_Integration = TIME_STEPPING_IMPLICIT_SECOND_ORDER_BACKWARD;
    } else if (strcmp(Physical_Time_Integration_Type,"Adams_Type") == 0) {
      i_Physical_Time_Integration = TIME_STEPPING_IMPLICIT_ADAMS_TYPE;
    } else if (strcmp(Physical_Time_Integration_Type,"Lees") == 0) {
      i_Physical_Time_Integration = TIME_STEPPING_IMPLICIT_LEES;
    } else if (strcmp(Physical_Time_Integration_Type,"Two_Step_Trapezoidal") == 0) {
      i_Physical_Time_Integration = TIME_STEPPING_IMPLICIT_TWO_STEP_TRAPEZOIDAL;
    } else if (strcmp(Physical_Time_Integration_Type,"A_Contractive") == 0) {
      i_Physical_Time_Integration = TIME_STEPPING_IMPLICIT_A_CONTRACTIVE;
    } else {
      i_command = INVALID_INPUT_VALUE;
    }

  } else if (strcmp(code, "Multigrid_N_Stage") == 0) {
    i_command = 2023;
    value >> N_Stage;
    if (N_Stage < 0) i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(code,"Physical_Time_CFL_Number") == 0) {
    i_command = 2024;
    value >> Physical_Time_CFL_Number;
    if (Physical_Time_CFL_Number <= ZERO) i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(code,"Dual_Time_Convergence_Residual_Level") == 0) {
    i_command = 2025;
    value >> Dual_Time_Convergence_Residual_Level;
    if (Dual_Time_Convergence_Residual_Level < ZERO) i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(code,"Multigrid_Write_Output_Cells_Frequency") == 0) {
    i_command = 2026;
    value >> Write_Output_Cells_Frequency;

  } else {
    i_command = INVALID_INPUT_CODE;

  } /* endif */

  return i_command;
  
}

/***************************************************************************
 * Multigrid_Input_Parameters::Check_Inputs -- Check input values.         *
 ***************************************************************************/
int Multigrid_Input_Parameters::Check_Inputs(void) {

  // Make sure the smoothing method is not Multigrid or NKS.
  if (i_Smoothing == TIME_STEPPING_MULTIGRID ||
      i_Smoothing == TIME_STEPPING_NKS) {
    cout << "Multigrid::Check_Inputs: invalid smoothing method" << endl;
    cout.flush();
    return 1;
  } /* endif */

  // Make sure the number of cylce levels is greater than two if full
  // multigrid has been selected.  Otherwise set the number of full
  // multigrid cylces to zero.
  if (Number_of_Cycles_per_Stage_for_Full_Multigrid > 0 &&
      Levels < 3) {
    Number_of_Cycles_per_Stage_for_Full_Multigrid = 0;
  } /* endif */

  // Make sure the time accurate flag is set appropriately.
  if (!i_Dual_Time_Stepping) {
    //IP.Time_Accurate = 0;
  } else {
    //IP.Time_Accurate = 1;
  } /* endif */

  if (Residual_Smoothing_Epsilon <= ZERO) {
     Residual_Smoothing = 0;
     Residual_Smoothing_Epsilon = ZERO;
  } else {
     Residual_Smoothing = 1;
  } /* endif */
  
  if (Residual_Smoothing_Gauss_Seidel_Iterations < 0) {
     Residual_Smoothing_Gauss_Seidel_Iterations = 0;
  } /* endif */

  // Input parameters are consistent.  Exit successfully.
  return 0;

}

/***************************************************************************
 * Multigrid_Input_Parameters -- Input-output operators.                   *
 ***************************************************************************/
ostream &operator << (ostream &out_file,
                      const Multigrid_Input_Parameters &IP) {

  IP.Output(out_file);
  return (out_file);

}

istream &operator >> (istream &in_file,
                      Multigrid_Input_Parameters &IP) {

  return in_file;

}

void Multigrid_Input_Parameters::Output(ostream &out_file) const {

  // Multigrid cycle parameters:
  out_file << "\n  -> Number of Levels in Multigrid: "
	   << Levels;
  out_file << "\n  -> Multigrid Cycle Type: "
	   << Cycle_Type;
  if (Number_of_Cycles_per_Stage_for_Full_Multigrid > 0) {
    out_file << "\n  -> Full Multigrid Startup: "
	     << Number_of_Cycles_per_Stage_for_Full_Multigrid
	     << " cycles per stage";
  }
  out_file << "\n  -> Defect correction: "
	   << Defect_Correction;
  out_file << "\n  -> Prolong using injection: "
	   << Prolong_Using_Injection;
  out_file << "\n  -> First Order Coarse Mesh Reconstruction: "
	   << First_Order_Coarse_Mesh_Reconstruction;
  out_file << "\n  -> Apply Coarse Mesh Boundary Conditions: "
	   << Apply_Coarse_Mesh_Boundary_Conditions;
  out_file << "\n  -> Injection at Dirichlet Boundary Conditions: "
	   << Injection_at_Dirichlet_Boundary_Conditions;
  out_file << "\n  -> Update Stability Switch: "
	   << Update_Stability_Switch;
  out_file << "\n  -> Maximum Number of Update Reductions: "
	   << Maximum_Number_of_Update_Reductions;

  // Multigrid smoother parameters:
  out_file << "\n  -> Multigrid Smoothing Type: "
	   << Smoothing_Type;
  out_file << "\n  -> Number of Smoothing Iterations on the Finest Multigrid Level:  "
	   << Number_of_Smooths_on_Finest_Level;
  out_file << "\n  -> Number of Multigrid Pre-Smoothing Iterations: "
	   << Number_of_Pre_Smooths;
  out_file << "\n  -> Number of Multigrid Post-Smoothing Iterations: "
	   << Number_of_Post_Smooths;
  out_file << "\n  -> Number of Smoothing Iterations on the Coarsest Multigrid Level:  "
	   << Number_of_Smooths_on_Coarsest_Level;
  if (Absolute_Convergence_Tolerance > 0) {
    out_file.setf(ios::scientific);
    out_file << "\n  -> Absolute Convergence Tolerance for Regular Multigrid: "
	     << Absolute_Convergence_Tolerance;
    out_file.unsetf(ios::scientific);
  }
  if (Relative_Convergence_Tolerance > 0) {
    out_file.setf(ios::scientific);
    out_file << "\n  -> Relative Convergence Tolerance for Regular Multigrid: "
	     << Relative_Convergence_Tolerance;
    out_file.unsetf(ios::scientific);
  }
  if (Number_of_Cycles_per_Stage_for_Full_Multigrid > 0 && 
      FMG_Absolute_Convergence_Tolerance > 0) {
    out_file.setf(ios::scientific);
    out_file << "\n  -> Absolute Convergence Tolerance for Full Multigrid: "
             << FMG_Absolute_Convergence_Tolerance;
    out_file.unsetf(ios::scientific);
  }
  if (Number_of_Cycles_per_Stage_for_Full_Multigrid > 0 && 
      FMG_Relative_Convergence_Tolerance > 0) {
    out_file.setf(ios::scientific);
    out_file << "\n  -> Relative Convergence Tolerance for Full Multigrid: "
	     << FMG_Relative_Convergence_Tolerance;
    out_file.unsetf(ios::scientific);
  }

  // Dual-time-stepping parameters:
  if (i_Dual_Time_Stepping) {
    out_file << "\n  -> Dual-time-stepping time integration type: "
	     << Physical_Time_Integration_Type;
    out_file << "\n  -> Number of regular multigrid cycles: "
	     << Ncycles_Regular_Multigrid;
    out_file << "\n  -> Number of full multigrid cycles: "
	     << Ncycles_Full_Multigrid;
    out_file << "\n  -> Physical-time-stepping CFL number: "
	     << Physical_Time_CFL_Number;
    out_file << "\n  -> Dual-time-stepping convergence residual level: "
	     << Dual_Time_Convergence_Residual_Level;
    out_file << "\n  -> Dual-time-step type: "
	     << i_Dual_Time_Preconditioning;
  }

  // Output control parameters:
  if (Write_Output_Cells_Frequency > 0) {
    out_file << "\n  -> Write output cells frequency for multigrid: "
	     << Write_Output_Cells_Frequency;
  }

}
