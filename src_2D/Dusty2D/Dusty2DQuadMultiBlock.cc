/**********************************************************************
 * Dusty2DQuadMultiBlock.cc: Multi-block versions of subroutines for  *
 *                           2D dusty multi-block quadrilateral mesh  *
 *                           solution classes.                        *
 **********************************************************************/

// Include 2D Dusty quadrilateral mesh solution header file.

#ifndef _DUSTY2D_QUAD_INCLUDED
#include "Dusty2DQuad.h"
#endif // _DUSTY2D_QUAD_INCLUDED

/**********************************************************************
 * Dusty2D_Quad_Block -- Multiple Block External Subroutines.         *
 **********************************************************************/

/**********************************************************************
 * Routine: Allocate                                                  *
 *                                                                    *
 * Allocate memory for 1D array of 2D quadrilateral multi-block       *
 * solution blocks.                                                   *
 *                                                                    *
 **********************************************************************/
Dusty2D_Quad_Block* Allocate(Dusty2D_Quad_Block *Soln_ptr,
                             Dusty2D_Input_Parameters &Input_Parameters) {

  // Allocate memory.
  Soln_ptr = new Dusty2D_Quad_Block[Input_Parameters.Number_of_Blocks_Per_Processor];

  // Assign the number of dusty2D variables.
  Soln_ptr[0].NUM_VAR_DUSTY2D = Input_Parameters.Wo.NUM_VAR_DUSTY2D;

  // Assign the residual variable.
  Soln_ptr[0].residual_variable = Input_Parameters.i_Residual_Variable;

  // Return memory location.
  return Soln_ptr;

}

/**********************************************************************
 * Routine: Deallocate                                                *
 *                                                                    *
 * Deallocate memory for 1D array of 2D quadrilateral multi-block     *
 * solution blocks.                                                   *
 *                                                                    *
 **********************************************************************/
Dusty2D_Quad_Block* Deallocate(Dusty2D_Quad_Block *Soln_ptr,
                               Dusty2D_Input_Parameters &Input_Parameters) {

  // Deallocate memory.
  for (int nb = 0 ; nb < Input_Parameters.Number_of_Blocks_Per_Processor; nb++) {
    if (Soln_ptr[nb].W != NULL) Soln_ptr[nb].deallocate();
  }
  delete []Soln_ptr;
  Soln_ptr = NULL;
  
  // Return memory location.
  return Soln_ptr;

}

/**********************************************************************
 * Routine: CreateInitialSolutionBlocks                               *
 *                                                                    *
 * This routine calls the Create_Initial_Solution_Blocks routine      *
 * found in AMR.h.  This interface code to the AMR.h version has been *
 * added to facilitate grid management.                               *
 *                                                                    *
 **********************************************************************/
Dusty2D_Quad_Block* CreateInitialSolutionBlocks(Grid2D_Quad_Block **InitMeshBlk,
						Dusty2D_Quad_Block *Soln_ptr,
						Dusty2D_Input_Parameters &Input_Parameters,
						QuadTreeBlock_DataStructure &QuadTree,
						AdaptiveBlockResourceList &GlobalSolnBlockList,
						AdaptiveBlock2D_List &LocalSolnBlockList) {

  int error_flag;

  // Create the initial solution blocks.
  Soln_ptr = Create_Initial_Solution_Blocks(InitMeshBlk,
					    Soln_ptr,
					    Input_Parameters,
					    QuadTree,
					    GlobalSolnBlockList,
					    LocalSolnBlockList);

  // Send solution information between neighbouring blocks to complete
  // prescription of initial data.
  error_flag = Send_All_Messages(Soln_ptr,
				 LocalSolnBlockList,
				 NUM_COMP_VECTOR2D,
				 ON);
  if (error_flag) return NULL;

  // Allocate memory for the electrostatic variables.
  if (Input_Parameters.Electrostatic) {
    for (int nb = 0; nb < LocalSolnBlockList.Nblk; nb++) {
      if (LocalSolnBlockList.Block[nb].used == ADAPTIVEBLOCK2D_USED) {
	Soln_ptr[nb].allocate_electrostatic();
      }
    }
  }

  // Determine the distance to the nearest wall distance.
  error_flag = Determine_Wall_Distance(Soln_ptr,
				       QuadTree,
				       LocalSolnBlockList,
				       Input_Parameters);
  if (error_flag) {
    cout << " NavierStokes2D ERROR: During wall distance calculation."
	 << "  Error #" << error_flag << "." << endl;
    return NULL;
  }
  CFDkit_Broadcast_MPI(&error_flag,1);
  if (error_flag) return NULL;

  // Return the solution blocks.
  return Soln_ptr;

}

/**********************************************************************
 * Routine: ICs                                                       *
 *                                                                    *
 * Assigns initial conditions and data to the solution variables of a *
 * 1D array of 2D quadrilateral multi-block solution blocks.          *
 *                                                                    *
 **********************************************************************/
void ICs(Dusty2D_Quad_Block *Soln_ptr,
         AdaptiveBlock2D_List &Soln_Block_List,
         Dusty2D_Input_Parameters &Input_Parameters) {

  int error_flag;
  Dusty2D_pState Wo[5];

  // Define various reference flow states.
  Wo[0] = Input_Parameters.Wo;
  Wo[1] = Input_Parameters.W1;
  Wo[1].p = 1.20*Wo[1].p;
  Wo[2] = Input_Parameters.W2;

  // Assign initial data for each solution block.
  for (int nb = 0; nb < Soln_Block_List.Nblk; nb++) {
    if (Soln_Block_List.Block[nb].used == ADAPTIVEBLOCK2D_USED) {
      // Set particle formulation flag.
      Soln_ptr[nb].Particles = Input_Parameters.Particles;
      // Set phase-interaction flag.
      Soln_ptr[nb].PhaseInteraction = Input_Parameters.PhaseInteraction;
      // Set flow type indicator (inviscid/viscous).
      Soln_ptr[nb].Flow_Type = Input_Parameters.FlowType;
      // Set flow geometry indicator (planar/axisymmetric).
      Soln_ptr[nb].Axisymmetric = Input_Parameters.Axisymmetric;
      // Set electrostatic flow indicator (on or off).
      Soln_ptr[nb].Electrostatic = Input_Parameters.Electrostatic;
      // Wall velocity and temperature.
      Soln_ptr[nb].Vwall = Input_Parameters.Vwall;
      Soln_ptr[nb].Twall = Input_Parameters.Twall;
      // Set initial data.
      ICs(Soln_ptr[nb],Input_Parameters,Wo);
      // Set the turbulence initial conditions if required.
      error_flag = Turbulence_ICs(Soln_ptr[nb],Input_Parameters,Wo);
      //if (error_flag) return error_flag;
    }
  }

}

/**********************************************************************
 * Routine: Copy_Electrostatic_Field_Variables                        *
 *                                                                    *
 * Copy the electrostatic field variables from the potential solution *
 * to the local solution block for a 1D array of 2D quadrilateral     *
 * multi-block solution blocks.                                       *
 *                                                                    *
 **********************************************************************/
int Copy_Electrostatic_Field_Variables(Dusty2D_Quad_Block *Soln_ptr,
				       Electrostatic2D_Quad_Block *ES_Soln_ptr,
				       AdaptiveBlock2D_List &Soln_Block_List) {

  int error_flag;

  // Copy the electrostatic field solution information.
  for (int nb = 0; nb < Soln_Block_List.Nblk; nb++) {
    if (Soln_Block_List.Block[nb].used == ADAPTIVEBLOCK2D_USED) {
      error_flag = Copy_Electrostatic_Field_Variables(Soln_ptr[nb],
						      ES_Soln_ptr[nb]);
      if (error_flag) return error_flag;
    }
  }

  // Electrostatic field solution information successfully copied.
  return 0;

}

/**********************************************************************
 * Routine: BCs                                                       *
 *                                                                    *
 * Apply boundary conditions at boundaries of a 1D array of 2D        *
 * quadrilateral multi-block solution blocks.                         *
 *                                                                    *
 **********************************************************************/
void BCs(Dusty2D_Quad_Block *Soln_ptr,
         AdaptiveBlock2D_List &Soln_Block_List,
	 Dusty2D_Input_Parameters &Input_Parameters) {

  int error_flag;

  // Prescribe boundary data for each solution block.
  for (int nb = 0; nb < Soln_Block_List.Nblk; nb++) {
    if (Soln_Block_List.Block[nb].used == ADAPTIVEBLOCK2D_USED) {
      BCs(Soln_ptr[nb],Input_Parameters);
    }
  }

}

/**********************************************************************
 * Routine: CFL                                                       *
 *                                                                    *
 * Determines the allowable global and local time steps (for explicit *
 * Euler time stepping scheme) for a 1D array of 2D quadrilateral     *
 * multi-block solution blocks according to the                       *
 * Courant-Friedrichs-Lewy condition.                                 *
 *                                                                    *
 **********************************************************************/
double CFL(Dusty2D_Quad_Block *Soln_ptr,
           AdaptiveBlock2D_List &Soln_Block_List,
	   Dusty2D_Input_Parameters &Input_Parameters) {

  double dtMin = MILLION;

  // Determine the allowable time step for each solution block.
  for (int nb = 0; nb < Soln_Block_List.Nblk; nb++) {
    if (Soln_Block_List.Block[nb].used == ADAPTIVEBLOCK2D_USED) {
      dtMin = min(dtMin,CFL(Soln_ptr[nb],Input_Parameters));
    }
  }

  // Return the global time step.
  return dtMin;

}

/**********************************************************************
 * Routine: Set_Global_TimeStep                                       *
 *                                                                    *
 * Assigns global time step to a 1D array of 2D quadrilateral         *
 * multi-block solution blocks for time-accurate calculations.        *
 *                                                                    *
 **********************************************************************/
void Set_Global_TimeStep(Dusty2D_Quad_Block *Soln_ptr,
                         AdaptiveBlock2D_List &Soln_Block_List,
                         const double &Dt_min) {

  for (int nb = 0; nb < Soln_Block_List.Nblk; nb++) {
    if (Soln_Block_List.Block[nb].used == ADAPTIVEBLOCK2D_USED) {
      Set_Global_TimeStep(Soln_ptr[nb],Dt_min);
    }
  }

}

/**********************************************************************
 * Routine: L1_Norm_Residual                                          *
 *                                                                    *
 * Determines the L1-norm of the solution residual for a 1D array of  *
 * 2D quadrilateral multi-block solution blocks.  Useful for          *
 * monitoring convergence of the solution for steady state problems.  *
 *                                                                    *
 **********************************************************************/
double L1_Norm_Residual(Dusty2D_Quad_Block *Soln_ptr,
                        AdaptiveBlock2D_List &Soln_Block_List) {

  double l1_norm = ZERO;

  // Calculate the L1-norm.  Sum the L1-norm for each solution block.
  for (int nb = 0; nb < Soln_Block_List.Nblk; nb++) {
    if (Soln_Block_List.Block[nb].used == ADAPTIVEBLOCK2D_USED) {
      l1_norm += L1_Norm_Residual(Soln_ptr[nb]);
    }
  }

  // Return the L1-norm.
  return l1_norm;

}

/**********************************************************************
 * Routine: L2_Norm_Residual                                          *
 *                                                                    *
 * Determines the L2-norm of the solution residual for a 1D array of  *
 * 2D quadrilateral multi-block solution blocks.  Useful for          *
 * monitoring convergence of the solution for steady state problems.  *
 *                                                                    *
 **********************************************************************/
double L2_Norm_Residual(Dusty2D_Quad_Block *Soln_ptr,
                        AdaptiveBlock2D_List &Soln_Block_List) {

  double l2_norm = ZERO;

  // Sum the square of the L2-norm for each solution block.
  for (int nb = 0; nb < Soln_Block_List.Nblk; nb++) {
    if (Soln_Block_List.Block[nb].used == ADAPTIVEBLOCK2D_USED) {
      l2_norm += sqr(L2_Norm_Residual(Soln_ptr[nb]));
    }
  }

  // Calculate the L2-norm for all blocks.
  l2_norm = sqrt(l2_norm);

  // Return the L2-norm.
  return l2_norm;

}

/**********************************************************************
 * Routine: Max_Norm_Residual                                         *
 *                                                                    *
 * Determines the maximum norm of the solution residual for a 1D      *
 * array of 2D quadrilateral multi-block solution blocks.  Useful for *
 * monitoring convergence of the solution for steady state problems.  *
 *                                                                    *
 **********************************************************************/
double Max_Norm_Residual(Dusty2D_Quad_Block *Soln_ptr,
                         AdaptiveBlock2D_List &Soln_Block_List) {

  double max_norm = ZERO;

  // Find the maximum norm for all solution blocks.
  for (int nb = 0; nb < Soln_Block_List.Nblk; nb++) {
    if (Soln_Block_List.Block[nb].used == ADAPTIVEBLOCK2D_USED) {
      max_norm = max(max_norm,Max_Norm_Residual(Soln_ptr[nb]));
    }
  }

  // Return the maximum norm.
  return max_norm;

}

/**********************************************************************
 * Routine: Evaluate_Limiters                                         *
 *                                                                    *
 * Set conditions to evaluate the limiters for a 1D array of 2D       *
 * quadrilateral multi-block solution blocks.                         *
 *                                                                    *
 **********************************************************************/
void Evaluate_Limiters(Dusty2D_Quad_Block *Soln_ptr,
                       AdaptiveBlock2D_List &Soln_Block_List) {

  for (int nb = 0; nb < Soln_Block_List.Nblk; nb++) {
    if (Soln_Block_List.Block[nb].used == ADAPTIVEBLOCK2D_USED) {
      Soln_ptr[nb].evaluate_limiters();
    }
  }

}

/**********************************************************************
 * Routine: Freeze_Limiters                                           *
 *                                                                    *
 * Set conditions to freeze the limiters for a 1D array of 2D         *
 * quadrilateral multi-block solution blocks.                         *
 *                                                                    *
 **********************************************************************/
void Freeze_Limiters(Dusty2D_Quad_Block *Soln_ptr,
                     AdaptiveBlock2D_List &Soln_Block_List) {

  for (int nb = 0; nb < Soln_Block_List.Nblk; nb++) {
    if (Soln_Block_List.Block[nb].used == ADAPTIVEBLOCK2D_USED) {
      Soln_ptr[nb].freeze_limiters();
    }
  }

}

/**********************************************************************
 * Routine:  Residual_Smoothing                                       *
 *                                                                    *
 * Applies implicit residual smoothing to a 1D array of 2D            *
 * 2D quadrilateral multi-block solution blocks.                      *
 *                                                                    *
 **********************************************************************/
void Residual_Smoothing(Dusty2D_Quad_Block *Soln_ptr,
                        AdaptiveBlock2D_List &Soln_Block_List,
                        Dusty2D_Input_Parameters &Input_Parameters,
                        const int I_Stage) {

  int k_residual;

  switch(Input_Parameters.i_Time_Integration) {
  case TIME_STEPPING_EXPLICIT_EULER :
    k_residual = 0;
    break;
  case TIME_STEPPING_EXPLICIT_PREDICTOR_CORRECTOR :
    k_residual = 0;
    break;
  case TIME_STEPPING_EXPLICIT_RUNGE_KUTTA :
    k_residual = 0;
    if (Input_Parameters.N_Stage == 4) {
      if (I_Stage == 4) {
	k_residual = 0;
      } else {
	k_residual = I_Stage - 1;
      }
    }
    break;
  case TIME_STEPPING_MULTISTAGE_OPTIMAL_SMOOTHING :
    k_residual = 0;
    break;
  default:
    k_residual = 0;
    break;
  }

  for (int nb = 0; nb < Soln_Block_List.Nblk; nb++) {
    if (Soln_Block_List.Block[nb].used == ADAPTIVEBLOCK2D_USED) {
      Residual_Smoothing(Soln_ptr[nb],
			 k_residual,
			 Input_Parameters.Residual_Smoothing_Epsilon, 
			 Input_Parameters.Residual_Smoothing_Gauss_Seidel_Iterations);
    }
  }

}

/**********************************************************************
 * Routine: Apply_Boundary_Flux_Correction                            *
 *                                                                    *
 * Apply flux corrections at boundaries of a 1D array of 2D           *
 * quadrilateral multi-block solution blocks to ensure that the       *
 * scheme is conservative at boundaries with resolution mesh changes. *
 *                                                                    *
 **********************************************************************/
void Apply_Boundary_Flux_Corrections(Dusty2D_Quad_Block *Soln_ptr,
                                     AdaptiveBlock2D_List &Soln_Block_List) {

  // Prescribe boundary data for each solution block.
  for (int nb = 0; nb < Soln_Block_List.Nblk; nb++) {
    if (Soln_Block_List.Block[nb].used == ADAPTIVEBLOCK2D_USED) {
      Apply_Boundary_Flux_Corrections(Soln_ptr[nb],
				      Soln_Block_List.Block[nb].nN,
				      Soln_Block_List.Block[nb].nS,
				      Soln_Block_List.Block[nb].nE,
				      Soln_Block_List.Block[nb].nW);
    }
  }

}

/**********************************************************************
 * Routine: Apply_Boundary_Flux_Corrections_Multistage_Explicit       *
 *                                                                    *
 * Apply flux corrections at boundaries of a 1D array of 2D           *
 * quadrilateral multi-block solution blocks to ensure that the       *
 * scheme is conservative at boundaries with resolution mesh changes. *
 *                                                                    *
 **********************************************************************/
void Apply_Boundary_Flux_Corrections_Multistage_Explicit(Dusty2D_Quad_Block *Soln_ptr,
                                                         AdaptiveBlock2D_List &Soln_Block_List,
                                                         Dusty2D_Input_Parameters &Input_Parameters,
   	                                                 const int I_Stage) {

  // Prescribe boundary data for each solution block.
  for (int nb = 0; nb < Soln_Block_List.Nblk; nb++) {
    if (Soln_Block_List.Block[nb].used == ADAPTIVEBLOCK2D_USED) {
      Apply_Boundary_Flux_Corrections_Multistage_Explicit(Soln_ptr[nb],
							  I_Stage,
							  Input_Parameters,
							  Soln_Block_List.Block[nb].nN,
							  Soln_Block_List.Block[nb].nS,
							  Soln_Block_List.Block[nb].nE,
							  Soln_Block_List.Block[nb].nW);
    }
  }

}

/**********************************************************************
 * Routine: dUdt_Residual_Evaluation                                  *
 *                                                                    *
 * This routine evaluates the stage solution residual for a 1D array  *
 * of 2D quadrilateral multi-block solution blocks.  A variety of     *
 * multistage explicit time integration and upwind finite-volume      *
 * spatial discretization procedures can be used depending on the     *
 * specified input values.                                            *
 *                                                                    *
 **********************************************************************/
int dUdt_Residual_Evaluation(Dusty2D_Quad_Block *Soln_ptr,
			     AdaptiveBlockResourceList &Global_Soln_Block_List,
			     AdaptiveBlock2D_List &Local_Soln_Block_List,
                             Dusty2D_Input_Parameters &Input_Parameters) {

  int error_flag;

  // Compute wall injection speed.
  error_flag = Wall_Injection_Speed(Soln_ptr,
				    Global_Soln_Block_List,
				    Local_Soln_Block_List,
				    Input_Parameters);
  if (error_flag) return error_flag;

  // Compute dimensionless wall distance.
  error_flag = Dimensionless_Wall_Distance(Soln_ptr,
					   Local_Soln_Block_List,
					   Input_Parameters);
  if (error_flag) return error_flag;

  // Compute dimensionless wall injection speed.
  error_flag = Dimensionless_Wall_Injection_Speed(Soln_ptr,
 						  Local_Soln_Block_List,
 						  Input_Parameters);
  if (error_flag) return error_flag;

  // Evaluate the stage solution residual for each solution block.
  for (int nb = 0; nb < Local_Soln_Block_List.Nblk; nb++) {
    if (Local_Soln_Block_List.Block[nb].used == ADAPTIVEBLOCK2D_USED) {
      error_flag = dUdt_Residual_Evaluation(Soln_ptr[nb],
					    Input_Parameters);
      if (error_flag) return error_flag;
    }
  }

  // Residuals for each quadrilateral multi-block solution block
  // successfully calculated.  Return.
  return 0;

}

/**********************************************************************
 * Routine: dUdt_Multistage_Explicit                                  *
 *                                                                    *
 * This routine evaluates the stage solution residual for a 1D array  *
 * of 2D quadrilateral multi-block solution blocks.  A variety of     *
 * multistage explicit time integration and upwind finite-volume      *
 * spatial discretization procedures can be used depending on the     *
 * specified input values.                                            *
 *                                                                    *
 **********************************************************************/
int dUdt_Multistage_Explicit(Dusty2D_Quad_Block *Soln_ptr,
			     AdaptiveBlockResourceList &Global_Soln_Block_List,
			     AdaptiveBlock2D_List &Local_Soln_Block_List,
                             Dusty2D_Input_Parameters &Input_Parameters,
   	                     const int I_Stage) {

  int error_flag;

  // Compute wall injection speed.
  error_flag = Wall_Injection_Speed(Soln_ptr,
				    Global_Soln_Block_List,
				    Local_Soln_Block_List,
				    Input_Parameters);
  if (error_flag) return error_flag;

  // Compute dimensionless wall distance.
  error_flag = Dimensionless_Wall_Distance(Soln_ptr,
					   Local_Soln_Block_List,
					   Input_Parameters);
  if (error_flag) return error_flag;

  // Compute dimensionless wall injection speed.
  error_flag = Dimensionless_Wall_Injection_Speed(Soln_ptr,
 						  Local_Soln_Block_List,
						  Input_Parameters);
  if (error_flag) return error_flag;

  // Evaluate the stage solution residual for each solution block.
  for (int nb = 0; nb < Local_Soln_Block_List.Nblk; nb++) {
    if (Local_Soln_Block_List.Block[nb].used == ADAPTIVEBLOCK2D_USED) {
      error_flag = dUdt_Multistage_Explicit(Soln_ptr[nb],
					    I_Stage,
					    Input_Parameters);
      if (error_flag) return error_flag;
    }
  }

  // Residuals for each quadrilateral multi-block solution block
  // successfully calculated.  Return.
  return 0;

}

/**********************************************************************
 * Routine: Update_Solution_Multistage_Explicit                       *
 *                                                                    *
 * This routine updates the solution for a 1D array of 2D             *
 * quadrilateral multi-block solution blocks. A variety of multistage *
 * explicit time integration and upwind finite-volume spatial         *
 * discretization procedures can be used depending on the specified   *
 * input values.                                                      *
 *                                                                    *
 **********************************************************************/
int Update_Solution_Multistage_Explicit(Dusty2D_Quad_Block *Soln_ptr,
                                        AdaptiveBlock2D_List &Soln_Block_List,
                                        Dusty2D_Input_Parameters &Input_Parameters,
   	                                const int I_Stage) {

  int error_flag;

  // Update the solution for each solution block.
  for (int nb = 0; nb < Soln_Block_List.Nblk; nb++) {
    if (Soln_Block_List.Block[nb].used == ADAPTIVEBLOCK2D_USED) {
      error_flag = Update_Solution_Multistage_Explicit(Soln_ptr[nb],
						       I_Stage,
						       Input_Parameters);
      if (error_flag) return error_flag;
    }
  }

  // Quadrilateral multi-block solution blocks successfully updated.
  return 0;

}

/**********************************************************************
 * Routine: Adaptive_Mesh_Refinement                                  *
 *                                                                    *
 * This routine calls the AMR routine.  Note that this level has been *
 * added to permit grid management.                                   *
 *                                                                    *
 **********************************************************************/
int Adaptive_Mesh_Refinement(Dusty2D_Quad_Block *Soln_ptr,
			     Dusty2D_Input_Parameters &Input_Parameters,
			     QuadTreeBlock_DataStructure &QuadTree,
			     AdaptiveBlockResourceList &GlobalSolnBlockList,
			     AdaptiveBlock2D_List &LocalSolnBlockList,
			     const int Set_New_Refinement_Flags) {

  int error_flag;

  Flag_Blocks_For_Refinement(Soln_ptr,
			     Input_Parameters,
			     QuadTree,
			     GlobalSolnBlockList,
			     LocalSolnBlockList);
  for (int nb = 0; nb < LocalSolnBlockList.Nblk; nb++) {
    if (LocalSolnBlockList.Block[nb].used == ADAPTIVEBLOCK2D_USED) {
      if (Soln_ptr[nb].Grid.Cell[3][3].Xc.x > ZERO &&
	  Soln_ptr[nb].Grid.BCtypeN[3] == BC_WALL_VISCOUS_ISOTHERMAL) {
	LocalSolnBlockList.RefineFlag[nb] = ADAPTIVEBLOCK2D_REFINE;
      }
    }
  }

  // Call AMR.
  error_flag = AMR(Soln_ptr,
		   Input_Parameters,
		   QuadTree,
		   GlobalSolnBlockList,
		   LocalSolnBlockList,
		   OFF,ON);
// 		   ON,ON);
  if (error_flag) return error_flag;

//   Dusty2D_cState Unew;
//   for (int nb = 0; nb < LocalSolnBlockList.Nblk; nb++) {
//     if (LocalSolnBlockList.Block[nb].used == ADAPTIVEBLOCK2D_USED) {
//       for (int j = Soln_ptr[nb].JCl; j <= Soln_ptr[nb].JCu; j++) {
// 	for (int i = Soln_ptr[nb].ICl; i <= Soln_ptr[nb].ICu; i++) {
// 	  if (Soln_ptr[nb].U[i][j].Unphysical_Properties() ||
// 	      Soln_ptr[nb].W[i][j].Unphysical_Properties()) {
// 	    if (!Soln_ptr[nb].U[i-1][j].Unphysical_Properties() && !Soln_ptr[nb].W[i-1][j].Unphysical_Properties()) Unew = Soln_ptr[nb].U[i-1][j];
// 	    else if (!Soln_ptr[nb].U[i-1][j-1].Unphysical_Properties() && !Soln_ptr[nb].W[i-1][j-1].Unphysical_Properties()) Unew = Soln_ptr[nb].U[i-1][j-1];
// 	    else if (!Soln_ptr[nb].U[i  ][j-1].Unphysical_Properties() && !Soln_ptr[nb].W[i  ][j-1].Unphysical_Properties()) Unew = Soln_ptr[nb].U[i  ][j-1];
// 	    else if (!Soln_ptr[nb].U[i+1][j-1].Unphysical_Properties() && !Soln_ptr[nb].W[i+1][j-1].Unphysical_Properties()) Unew = Soln_ptr[nb].U[i+1][j-1];
// 	    else if (!Soln_ptr[nb].U[i+1][j  ].Unphysical_Properties() && !Soln_ptr[nb].W[i+1][j  ].Unphysical_Properties()) Unew = Soln_ptr[nb].U[i+1][j  ];
// 	    else if (!Soln_ptr[nb].U[i+1][j+1].Unphysical_Properties() && !Soln_ptr[nb].W[i+1][j+1].Unphysical_Properties()) Unew = Soln_ptr[nb].U[i+1][j+1];
// 	    else if (!Soln_ptr[nb].U[i  ][j+1].Unphysical_Properties() && !Soln_ptr[nb].W[i  ][j+1].Unphysical_Properties()) Unew = Soln_ptr[nb].U[i  ][j+1];
// 	    else if (!Soln_ptr[nb].U[i-1][j+1].Unphysical_Properties() && !Soln_ptr[nb].W[i-1][j+1].Unphysical_Properties()) Unew = Soln_ptr[nb].U[i-1][j+1];
// 	    else return 1;
// 	    Soln_ptr[nb].U[i][j] = Unew;
// 	    Soln_ptr[nb].W[i][j] = W(Unew);
// 	  }
// 	}
//       }
//     }
//   }
//   error_flag = Send_All_Messages(Soln_ptr,
// 				 LocalSolnBlockList,
// 				 Input_Parameters.Wo.NUM_VAR_DUSTY2D,
// 				 OFF);
//   if (error_flag) {
//     cout << "\n Dusty2D ERROR: Message passing error during Dusty2D "
// 	 << "solution intialization on processor "
// 	 << LocalSolnBlockList.ThisCPU
// 	 << "." << endl;
//   }
//   error_flag = CFDkit_OR_MPI(error_flag);
//   if (error_flag) return error_flag;

//   ICs(Soln_ptr,
//       LocalSolnBlockList,
//       Input_Parameters);

  // Determine the wall distance and wall distance location if required.
  error_flag = Determine_Wall_Distance(Soln_ptr,
				       QuadTree,
				       LocalSolnBlockList,
				       Input_Parameters);
  if (error_flag) return error_flag;

  // AMR successful.
  return 0;

}

/**********************************************************************
 * Routine: Initial_Adaptive_Mesh_Refinement                          *
 *                                                                    *
 * This routine calls the Initial_AMR routine.  Note that this level  * 
 * has been added to facilitate grid management.                      *
 *                                                                    *
 **********************************************************************/
int Initial_Adaptive_Mesh_Refinement(Dusty2D_Quad_Block *Soln_ptr,
				     Dusty2D_Input_Parameters &Input_Parameters,
				     QuadTreeBlock_DataStructure &QuadTree,
				     AdaptiveBlockResourceList &GlobalSolnBlockList,
				     AdaptiveBlock2D_List &LocalSolnBlockList) {

  // Exit immediately if no initial mesh refinement is required.
  if (!Input_Parameters.Number_of_Initial_Mesh_Refinements) return 0;

  int error_flag;

  // Call Initial AMR.
  error_flag = Initial_AMR(Soln_ptr,
			   Input_Parameters,
			   QuadTree,
			   GlobalSolnBlockList,
			   LocalSolnBlockList);
  if (error_flag) return error_flag;

  // Determine the wall distance and wall distance location if required.
  error_flag = Determine_Wall_Distance(Soln_ptr,
				       QuadTree,
				       LocalSolnBlockList,
				       Input_Parameters);
  if (error_flag) return error_flag;

  // Re-apply initial conditions.
  ICs(Soln_ptr,LocalSolnBlockList,Input_Parameters);

  // Initial AMR successful.
  return 0;

}

/**********************************************************************
 * Routine: Uniform_Adaptive_Mesh_Refinement                          *
 *                                                                    *
 * This routine calls the Uniform_AMR routine.  Note that this level  * 
 * has been added to facilitate grid management.                      *
 *                                                                    *
 **********************************************************************/
int Uniform_Adaptive_Mesh_Refinement(Dusty2D_Quad_Block *Soln_ptr,
				     Dusty2D_Input_Parameters &Input_Parameters,
				     QuadTreeBlock_DataStructure &QuadTree,
				     AdaptiveBlockResourceList &GlobalSolnBlockList,
				     AdaptiveBlock2D_List &LocalSolnBlockList) {

  // Exit immediately if no uniform mesh refinement is required.
  if (!Input_Parameters.Number_of_Uniform_Mesh_Refinements) return 0;

  int error_flag;

  // Call Uniform AMR.
  error_flag = Uniform_AMR(Soln_ptr,
			   Input_Parameters,
			   QuadTree,
			   GlobalSolnBlockList,
			   LocalSolnBlockList);
  if (error_flag) return error_flag;

  // Determine the wall distance and wall distance location if required.
  error_flag = Determine_Wall_Distance(Soln_ptr,
				       QuadTree,
				       LocalSolnBlockList,
				       Input_Parameters);
  if (error_flag) return error_flag;

  // Re-apply initial conditions.
  ICs(Soln_ptr,LocalSolnBlockList,Input_Parameters);

  // Initial uniform AMR successful.
  return 0;

  // Initial uniform AMR successful.
  return 0;

}

/**********************************************************************
 * Routine: Boundary_Adaptive_Mesh_Refinement                         *
 *                                                                    *
 * This routine calls the Boundary_AMR routine.  Note that this level * 
 * has been added to facilitate grid management.                      *
 *                                                                    *
 **********************************************************************/
int Boundary_Adaptive_Mesh_Refinement(Dusty2D_Quad_Block *Soln_ptr,
				      Dusty2D_Input_Parameters &Input_Parameters,
				      QuadTreeBlock_DataStructure &QuadTree,
				      AdaptiveBlockResourceList &GlobalSolnBlockList,
				      AdaptiveBlock2D_List &LocalSolnBlockList) {


  // Exit immediately if no boundary mesh refinement is required.
  if (!Input_Parameters.Number_of_Boundary_Mesh_Refinements) return 0;

  int error_flag;

  // Call Boundary AMR.
  error_flag = Boundary_AMR(Soln_ptr,
			    Input_Parameters,
			    QuadTree,
			    GlobalSolnBlockList,
			    LocalSolnBlockList);
  if (error_flag) return error_flag;

  // Determine the wall distance and wall distance location if required.
  error_flag = Determine_Wall_Distance(Soln_ptr,
				       QuadTree,
				       LocalSolnBlockList,
				       Input_Parameters);
  if (error_flag) return error_flag;

  // Re-apply initial conditions.
  ICs(Soln_ptr,LocalSolnBlockList,Input_Parameters);

  // Initial boundary AMR successful.
  return 0;

}

/**********************************************************************
 * Routine: Multi_Velocity_Component_Particle_Phase_Switch            *
 *                                                                    *
 * This routine conducts particle-phase velocity component switching  *
 * for a 1D array of 2D quadrilateral multi-block solution blocks.    *
 *                                                                    *
 **********************************************************************/
int Multi_Velocity_Component_Particle_Phase_Switch(Dusty2D_Quad_Block *Soln_ptr,
						   Dusty2D_Input_Parameters &Input_Parameters,
						   AdaptiveBlock2D_List &LocalSolnBlockList) {

  // Exit immediately if no particle families are being simulated.
  if (Input_Parameters.Particles != PARTICLE_PHASE_EULERIAN_FORMULATION) return 0;

  // Exit immediately if the single-velocity component formulation is
  // being used for the particle-phase.
  if (Input_Parameters.num_cmp_part == PARTICLE2D_SINGLE_VELOCITY_FORMULATION) return 0;

  int error_flag;

  //  Conduct particle-phase component switching on each solution block.
  for (int nb = 0; nb < LocalSolnBlockList.Nblk; nb++) {
    if (LocalSolnBlockList.Block[nb].used == ADAPTIVEBLOCK2D_USED) {
      error_flag = Multi_Velocity_Component_Particle_Phase_Switch(Soln_ptr[nb],
								  Input_Parameters);
      if (error_flag) return error_flag;
    }
  }

  // Ensure boundary conditions are consistent.
  BCs(Soln_ptr,
      LocalSolnBlockList,
      Input_Parameters);

  // Particle-phase velocity component switching computed successfully.
  return 0;

}

/**********************************************************************
 * Routine: Determine_Conservation_Properties                         *
 **********************************************************************/
int Determine_Conservation_Properties(Dusty2D_Quad_Block *Soln_ptr,
				      AdaptiveBlock2D_List &LocalSolnBlockList,
				      Dusty2D_cState &Umass) {

  int error_flag;

  Umass.Vacuum();

  for (int nb = 0; nb < LocalSolnBlockList.Nblk; nb++) {
    if (LocalSolnBlockList.Block[nb].used == ADAPTIVEBLOCK2D_USED) {
      error_flag = Determine_Conservation_Properties(Soln_ptr[nb],Umass);
      if (error_flag) return error_flag;
    }
  }

  return 0;

}
