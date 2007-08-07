/**********************************************************************
 * HighTemp2DQuadMultiBlock.cc: Multi-block versions of           *
 *                                  subroutines for the 2D High-    *
 *                                  Temp multi-block quadrilateral  *
 *                                  mesh solution classes.            *
 **********************************************************************/

#include "HighTemp2DQuad.h"

/**********************************************************************
 * HighTemp2D_Quad_Block -- Multiple Block External Subroutines.  *
 **********************************************************************/

/**********************************************************************
 * Routine: Allocate                                                  *
 *                                                                    *
 * Allocate memory for 1D array of 2D quadrilateral multi-block       *
 * solution blocks.                                                   *
 *                                                                    *
 **********************************************************************/
HighTemp2D_Quad_Block* Allocate(HighTemp2D_Quad_Block *Soln_ptr,
				    HighTemp2D_Input_Parameters &Input_Parameters) {

  // Allocate memory.
  Soln_ptr = new HighTemp2D_Quad_Block[Input_Parameters.Number_of_Blocks_Per_Processor];

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
HighTemp2D_Quad_Block* Deallocate(HighTemp2D_Quad_Block *Soln_ptr,
				      HighTemp2D_Input_Parameters &Input_Parameters) {

  // Deallocate memory.
  for (int nb = 0 ; nb < Input_Parameters.Number_of_Blocks_Per_Processor; nb++) {
    if (Soln_ptr[nb].W != NULL) Soln_ptr[nb].deallocate();
  }
  if (Soln_ptr != NULL) {  delete []Soln_ptr; Soln_ptr = NULL; }
  
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
HighTemp2D_Quad_Block* CreateInitialSolutionBlocks(Grid2D_Quad_Block **InitMeshBlk,
						       HighTemp2D_Quad_Block *Soln_ptr,
						       HighTemp2D_Input_Parameters &Input_Parameters,
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

  // Determine the distance to the nearest wall distance.
  error_flag = Determine_Wall_Distance(Soln_ptr,
				       QuadTree,
				       LocalSolnBlockList,
				       Input_Parameters);
  if (error_flag) {
    cout << " HighTemp2D ERROR: During wall distance calculation."
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
void ICs(HighTemp2D_Quad_Block *Soln_ptr,
         AdaptiveBlock2D_List &Soln_Block_List,
         HighTemp2D_Input_Parameters &Input_Parameters) {

  int error_flag;
  HighTemp2D_pState Wo[5];

  // Define various reference flow states.
  Wo[0] = Input_Parameters.Wo;
  Wo[1] = Input_Parameters.W1;
  Wo[1].p = 1.20*Wo[1].p;
  Wo[2] = Input_Parameters.W2;

  // Assign initial data for each solution block.
  for (int nb = 0; nb < Soln_Block_List.Nblk; nb++) {
    if (Soln_Block_List.Block[nb].used == ADAPTIVEBLOCK2D_USED) {
      // Set flow type indicator (inviscid/viscous).
      Soln_ptr[nb].Flow_Type = Input_Parameters.FlowType;
      // Set flow geometry indicator (planar/axisymmetric).
      Soln_ptr[nb].Axisymmetric = Input_Parameters.Axisymmetric;
      // Sub-cell reconstruction indicator (linear or quadratic).
      //   Soln_ptr[nb].SubCell_Reconstruction = Input_Parameters.i_SubCell_Reconstruction;
      // Wall velocity and temperature.
      Soln_ptr[nb].Vwall = Input_Parameters.Vwall;
      Soln_ptr[nb].Twall = Input_Parameters.Twall;
      // Set initial data.
      ICs(Soln_ptr[nb],Input_Parameters,Wo);
      // Set the turbulence initial conditions if required.
      error_flag = Turbulent_ICs(Soln_ptr[nb],Input_Parameters,Wo);
      //if (error_flag) return error_flag;
    }
  }

}

/**********************************************************************
 * Routine: Reinitialize_Wedge_Initial_Conditions                     *
 *                                                                    *
 * Reassigns the initial conditions before the incident shock-wave    *
 * for wedge shock relfection problems for a 1D array of 2D           *
 * quadrilateral multi-block solution blocks.                         *
 *                                                                    *
 **********************************************************************/
int Reinitialize_Wedge_Initial_Conditions(HighTemp2D_Quad_Block *Soln_ptr,
					  AdaptiveBlock2D_List &Soln_Block_List,
					  HighTemp2D_Input_Parameters &Input_Parameters) {

//   // Reassign initial data for each solution block.
//   for (int nb = 0; nb < Soln_Block_List.Nblk; nb++) {
//     if (Soln_Block_List.Block[nb].used == ADAPTIVEBLOCK2D_USED) {
//       for (int j = Soln_ptr[nb].JCl-Soln_ptr[nb].Nghost; j <= Soln_ptr[nb].JCu+Soln_ptr[nb].Nghost; j++) {
// 	for (int i = Soln_ptr[nb].ICl-Soln_ptr[nb].Nghost; i <= Soln_ptr[nb].ICu+Soln_ptr[nb].Nghost; i++) {
// 	  if (Soln_ptr[nb].Grid.Cell[i][j].Xc.x < ZERO) {
// 	    Soln_ptr[nb].W[i][j].rho = 1.84671;
// 	    Soln_ptr[nb].W[i][j].v.x = 1.84671;
// 	    Soln_ptr[nb].W[i][j].v.y = 0.00119507;
// 	    Soln_ptr[nb].W[i][j].p = 205540.0;
// 	    Soln_ptr[nb].U[i][j] = U(Soln_ptr[nb].W[i][j]);
// 	  }
// 	}
//       }
//     }
//   }

  return 0;

}

/**********************************************************************
 * Routine: BCs                                                       *
 *                                                                    *
 * Apply boundary conditions at boundaries of a 1D array of 2D        *
 * quadrilateral multi-block solution blocks.                         *
 *                                                                    *
 **********************************************************************/
void BCs(HighTemp2D_Quad_Block *Soln_ptr,
         AdaptiveBlock2D_List &Soln_Block_List,
	 const HighTemp2D_Input_Parameters &Input_Parameters) {

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
 * N-S HighTemp time stepping scheme) for a 1D array of 2D           *
 * quadrilateral multi-block solution blocks according to the         *
 * Courant-Friedrichs-Lewy condition.                                 *
 *                                                                    *
 **********************************************************************/
double CFL(HighTemp2D_Quad_Block *Soln_ptr,
           AdaptiveBlock2D_List &Soln_Block_List,
	   HighTemp2D_Input_Parameters &Input_Parameters) {

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
void Set_Global_TimeStep(HighTemp2D_Quad_Block *Soln_ptr,
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
double L1_Norm_Residual(HighTemp2D_Quad_Block *Soln_ptr,
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
double L2_Norm_Residual(HighTemp2D_Quad_Block *Soln_ptr,
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
double Max_Norm_Residual(HighTemp2D_Quad_Block *Soln_ptr,
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
void Evaluate_Limiters(HighTemp2D_Quad_Block *Soln_ptr,
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
void Freeze_Limiters(HighTemp2D_Quad_Block *Soln_ptr,
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
 * quadrilateral multi-block solution blocks.                         *
 *                                                                    *
 **********************************************************************/
void Residual_Smoothing(HighTemp2D_Quad_Block *Soln_ptr,
                        AdaptiveBlock2D_List &Soln_Block_List,
                        HighTemp2D_Input_Parameters &Input_Parameters,
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
void Apply_Boundary_Flux_Corrections(HighTemp2D_Quad_Block *Soln_ptr,
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
void Apply_Boundary_Flux_Corrections_Multistage_Explicit(HighTemp2D_Quad_Block *Soln_ptr,
                                                         AdaptiveBlock2D_List &Soln_Block_List,
                                                         HighTemp2D_Input_Parameters &Input_Parameters,
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
int dUdt_Residual_Evaluation(HighTemp2D_Quad_Block *Soln_ptr,
			     AdaptiveBlockResourceList &Global_Soln_Block_List,
                             AdaptiveBlock2D_List &Local_Soln_Block_List,
                             HighTemp2D_Input_Parameters &Input_Parameters) {

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
int dUdt_Multistage_Explicit(HighTemp2D_Quad_Block *Soln_ptr,
			     AdaptiveBlockResourceList &Global_Soln_Block_List,
                             AdaptiveBlock2D_List &Local_Soln_Block_List,
                             HighTemp2D_Input_Parameters &Input_Parameters,
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
int Update_Solution_Multistage_Explicit(HighTemp2D_Quad_Block *Soln_ptr,
                                        AdaptiveBlock2D_List &Soln_Block_List,
                                        HighTemp2D_Input_Parameters &Input_Parameters,
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

  // Apply turbulence boundary conditions.
  // error_flag = Turbulent_BCs(Soln_ptr,
  //			     Soln_Block_List,
  //			     Input_Parameters);
  // if (error_flag) return error_flag;

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
int Adaptive_Mesh_Refinement(HighTemp2D_Quad_Block *Soln_ptr,
			     HighTemp2D_Input_Parameters &Input_Parameters,
			     QuadTreeBlock_DataStructure &QuadTree,
			     AdaptiveBlockResourceList &GlobalSolnBlockList,
			     AdaptiveBlock2D_List &LocalSolnBlockList) {

  int error_flag;

  // Call AMR.
  error_flag = AMR(Soln_ptr,
		   Input_Parameters,
		   QuadTree,
		   GlobalSolnBlockList,
		   LocalSolnBlockList,
		   ON,ON);
  if (error_flag) return error_flag;

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
int Initial_Adaptive_Mesh_Refinement(HighTemp2D_Quad_Block *Soln_ptr,
				     HighTemp2D_Input_Parameters &Input_Parameters,
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
int Uniform_Adaptive_Mesh_Refinement(HighTemp2D_Quad_Block *Soln_ptr,
				     HighTemp2D_Input_Parameters &Input_Parameters,
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

}

/**********************************************************************
 * Routine: Boundary_Adaptive_Mesh_Refinement                         *
 *                                                                    *
 * This routine calls the Boundary_AMR routine.  Note that this level * 
 * has been added to facilitate grid management.                      *
 *                                                                    *
 **********************************************************************/
int Boundary_Adaptive_Mesh_Refinement(HighTemp2D_Quad_Block *Soln_ptr,
				      HighTemp2D_Input_Parameters &Input_Parameters,
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
 * Routine: Flat_Plate_Adaptive_Mesh_Refinement                       *
 *                                                                    *
 * Forced refinement on a flat plate.                                 *
 *                                                                    *
 **********************************************************************/
int Flat_Plate_Adaptive_Mesh_Refinement(HighTemp2D_Quad_Block *Soln_ptr,
					HighTemp2D_Input_Parameters &Input_Parameters,
					QuadTreeBlock_DataStructure &QuadTree,
					AdaptiveBlockResourceList &GlobalSolnBlockList,
					AdaptiveBlock2D_List &LocalSolnBlockList) {

  // Exit immediately if no flat-plate mesh refinement is required.
  if (!Input_Parameters.Number_of_Flat_Plate_Mesh_Refinements) return 0;

  int error_flag;

  // Otherwise refine the mesh the number of specified times.
  for (int number_of_flat_plate_mesh_refinements = 1;
       number_of_flat_plate_mesh_refinements <= Input_Parameters.Number_of_Flat_Plate_Mesh_Refinements;
       number_of_flat_plate_mesh_refinements++) {

    // Set refinement flags.
    QuadTree.nochangeAll();
    for (int nb = 0; nb < LocalSolnBlockList.Nblk; nb++) {
      LocalSolnBlockList.RefineFlag[nb] = ADAPTIVEBLOCK2D_NOCHANGE;
      if (LocalSolnBlockList.Block[nb].used == ADAPTIVEBLOCK2D_USED) {
	for (int i = Soln_ptr[nb].ICl-Soln_ptr[nb].Nghost; i <= Soln_ptr[nb].ICu+Soln_ptr[nb].Nghost; i++) {
	  if ((Soln_ptr[nb].Grid.BCtypeS[i] == BC_WALL_VISCOUS_HEATFLUX ||
	       Soln_ptr[nb].Grid.BCtypeS[i] == BC_WALL_VISCOUS_ISOTHERMAL) &&
	      Soln_ptr[nb].Grid.Cell[2][2].Xc.x < 0.25*Input_Parameters.Plate_Length) {
	    if (LocalSolnBlockList.Block[nb].info.level < Input_Parameters.Maximum_Refinement_Level) {
	      LocalSolnBlockList.RefineFlag[nb] = ADAPTIVEBLOCK2D_REFINE;
	    }
	  }
	}
      }
    }

    // Set the global refinement list.
    Get_Refinement_List(QuadTree,
                        LocalSolnBlockList);

    // Call AMR.
    error_flag = AMR(Soln_ptr,
		     Input_Parameters,
		     QuadTree,
		     GlobalSolnBlockList,
		     LocalSolnBlockList,
		     OFF,ON);
    if (error_flag) return error_flag;

    // Send all messages.
    error_flag = Send_All_Messages(Soln_ptr,
				   LocalSolnBlockList,
				   Soln_ptr[0].NumVar(),
				   ON);
    if (error_flag) return error_flag;

    // Output the refinement statistics.
    if (CFDkit_Primary_MPI_Processor()) {
      cout << "\n Refinement Level #" << number_of_flat_plate_mesh_refinements
	   << " : Number of Blocks = " << QuadTree.countUsedBlocks()
	   << ", Number of Cells = " << QuadTree.countUsedCells()
	   << ", Refinement Efficiency = " << QuadTree.efficiencyRefinement();
    }

  }

  // Determine the wall distance and wall distance location if required.
  error_flag = Determine_Wall_Distance(Soln_ptr,
				       QuadTree,
				       LocalSolnBlockList,
				       Input_Parameters);
  if (error_flag) return error_flag;

  // Initial conditions.
  ICs(Soln_ptr,LocalSolnBlockList,Input_Parameters);

  // Send all messages.
  error_flag = Send_All_Messages(Soln_ptr,
				 LocalSolnBlockList,
				 Soln_ptr[0].NumVar(),
				 ON);
  if (error_flag) return error_flag;

  // Initial flat-plate AMR successful.
  return 0;

}

/**********************************************************************
 * Routine: Perform_External_Functions_on_Coarse_Grids                *
 *                                                                    *
 * This routine performs and functions required on the coarse grid    *
 * that were already performed on the finest grid.  This routine is   *
 * not included as a templated routine in the multigrid code since    *
 * the functions to be performed are dependent on the system of       *
 * equations that are being solved.  These routines are perfomed on a *
 * 1D array of 2D quadrilateral multi-block solution blocks.          *
 *                                                                    *
 **********************************************************************/
/*
int Perform_External_Functions_on_Coarse_Grids(HighTemp2D_Quad_Block *Soln_ptr_Coarse,
					       HighTemp2D_Quad_Block *Soln_ptr_Fine,
					       HighTemp2D_Input_Parameters &Input_Parameters,
					       QuadTreeBlock_DataStructure &QuadTree,
					       AdaptiveBlockResourceList &GlobalSolnBlockList,
					       AdaptiveBlock2D_List &LocalSolnBlockList) {

  int error_flag;

  // Determine the wall distance and wall distance location if required.
  error_flag = Determine_Wall_Distance(Soln_ptr_Coarse,
				       QuadTree,
				       LocalSolnBlockList,
				       Input_Parameters);
  if (error_flag) return error_flag;

  // All required functions computed successfully.
  return 0;

}
*/
/**********************************************************************
 * Routine: Zero_Residuals_on_Coarse_Grid                             *
 *                                                                    *
 * Zeroes the residuals for the turbulence parameters k and omega of  *
 * a 1D array of 2D quadrilateral multi-block solution blocks.  Used  *
 * by the FAS-multigrid routine so that k and omega are not updated   *
 * on the coarse grid levels.                                         *
 *                                                                    *
 **********************************************************************/
/*
int Zero_Residuals_on_Coarse_Grid(HighTemp2D_Quad_Block *Soln_ptr,
				  AdaptiveBlock2D_List &Soln_Block_List,
				  HighTemp2D_Input_Parameters &Input_Parameters,
				  const int &I_Stage) {

  // Exit immediately if not a turbulent flow.
  if (Input_Parameters.FlowType != FLOWTYPE_TURBULENT_RANS_K_OMEGA) return 0;

  int error_flag;

  // Zero the required residuals on each solution block in use.
  for (int nb = 0; nb < Soln_Block_List.Nblk; nb++) {
    if (Soln_Block_List.Block[nb].used == ADAPTIVEBLOCK2D_USED) {
      error_flag = Zero_Residuals_on_Coarse_Grid(Soln_ptr[nb],
						 Input_Parameters,
						 I_Stage);
      if (error_flag) return error_flag;
    }
  }

  // Residuals for k and omega zeroed.
  return 0;

}
*/
