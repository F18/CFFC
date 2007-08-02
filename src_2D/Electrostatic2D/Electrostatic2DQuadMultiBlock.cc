/**********************************************************************
 * Electrostatic2DQuadMultiBlock.cc: Multi-block versions of          *
 *                                   subroutines for 2D electrostatic *
 *                                   multi-block quadrilateral mesh   *
 *                                   solution classes.                *
 **********************************************************************/

// Include 2D Electrostatic quadrilateral mesh solution header file.

#ifndef _ELECTROSTATIC2D_QUAD_INCLUDED
#include "Electrostatic2DQuad.h"
#endif // _ELECTROSTATIC2D_QUAD_INCLUDED

/**********************************************************************
 * Electrostatic2D_Quad_Block -- Multiple Block External Subroutines. *
 **********************************************************************/

/**********************************************************************
 * Routine: Allocate                                                  *
 *                                                                    *
 * Allocate memory for 1D array of 2D quadrilateral multi-block       *
 * solution blocks.                                                   *
 *                                                                    *
 **********************************************************************/
Electrostatic2D_Quad_Block* Allocate(Electrostatic2D_Quad_Block *Soln_ptr,
				     Electrostatic2D_Input_Parameters &Input_Parameters) {

  // Allocate memory.
  Soln_ptr = new Electrostatic2D_Quad_Block[Input_Parameters.Number_of_Blocks_Per_Processor];

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
Electrostatic2D_Quad_Block* Deallocate(Electrostatic2D_Quad_Block *Soln_ptr,
				       Electrostatic2D_Input_Parameters &Input_Parameters) {

  // Deallocate memory.
  for (int nb = 0 ; nb < Input_Parameters.Number_of_Blocks_Per_Processor; nb++) {
    if (Soln_ptr[nb].U != NULL) Soln_ptr[nb].deallocate();
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
Electrostatic2D_Quad_Block* CreateInitialSolutionBlocks(Grid2D_Quad_Block **InitMeshBlk,
							Electrostatic2D_Quad_Block *Soln_ptr,
							Electrostatic2D_Input_Parameters &Input_Parameters,
							QuadTreeBlock_DataStructure &QuadTree,
							AdaptiveBlockResourceList &GlobalSolnBlockList,
							AdaptiveBlock2D_List &LocalSolnBlockList) {

  // Create the initial solution blocks.
  Soln_ptr = Create_Initial_Solution_Blocks(InitMeshBlk,
					    Soln_ptr,
					    Input_Parameters,
					    QuadTree,
					    GlobalSolnBlockList,
					    LocalSolnBlockList);

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
void ICs(Electrostatic2D_Quad_Block *Soln_ptr,
         AdaptiveBlock2D_List &Soln_Block_List,
         Electrostatic2D_Input_Parameters &Input_Parameters) {

  int error_flag;
  Electrostatic2DState Uo[5];

  // Define various reference flow states.
  Uo[0] = Input_Parameters.Uo;
  Uo[1] = Input_Parameters.U1;
  Uo[2] = Input_Parameters.U2;

  // Assign initial data for each solution block.
  for (int nb = 0; nb < Soln_Block_List.Nblk; nb++) {
    if (Soln_Block_List.Block[nb].used == ADAPTIVEBLOCK2D_USED) {
      // Set flow geometry indicator (planar/axisymmetric).
      Soln_ptr[nb].Axisymmetric = Input_Parameters.Axisymmetric;
      // Set initial data.
      ICs(Soln_ptr[nb],Input_Parameters,Uo);
    }
  }

}

/**********************************************************************
 * Routine: BCs                                                       *
 *                                                                    *
 * Apply boundary conditions at boundaries of a 1D array of 2D        *
 * quadrilateral multi-block solution blocks.                         *
 *                                                                    *
 **********************************************************************/
void BCs(Electrostatic2D_Quad_Block *Soln_ptr,
         AdaptiveBlock2D_List &Soln_Block_List,
	 Electrostatic2D_Input_Parameters &Input_Parameters) {

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
double CFL(Electrostatic2D_Quad_Block *Soln_ptr,
           AdaptiveBlock2D_List &Soln_Block_List,
	   Electrostatic2D_Input_Parameters &Input_Parameters) {

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
void Set_Global_TimeStep(Electrostatic2D_Quad_Block *Soln_ptr,
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
double L1_Norm_Residual(Electrostatic2D_Quad_Block *Soln_ptr,
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
double L2_Norm_Residual(Electrostatic2D_Quad_Block *Soln_ptr,
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
double Max_Norm_Residual(Electrostatic2D_Quad_Block *Soln_ptr,
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
 * Routine:  Residual_Smoothing                                       *
 *                                                                    *
 * Applies implicit residual smoothing to a 1D array of 2D            *
 * 2D quadrilateral multi-block solution blocks.                      *
 *                                                                    *
 **********************************************************************/
void Residual_Smoothing(Electrostatic2D_Quad_Block *Soln_ptr,
                        AdaptiveBlock2D_List &Soln_Block_List,
                        Electrostatic2D_Input_Parameters &Input_Parameters,
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
void Apply_Boundary_Flux_Corrections(Electrostatic2D_Quad_Block *Soln_ptr,
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
void Apply_Boundary_Flux_Corrections_Multistage_Explicit(Electrostatic2D_Quad_Block *Soln_ptr,
                                                         AdaptiveBlock2D_List &Soln_Block_List,
                                                         Electrostatic2D_Input_Parameters &Input_Parameters,
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
int dUdt_Residual_Evaluation(Electrostatic2D_Quad_Block *Soln_ptr,
                             AdaptiveBlock2D_List &Soln_Block_List,
                             Electrostatic2D_Input_Parameters &Input_Parameters) {

  int error_flag;

  // Evaluate the stage solution residual for each solution block.
  for (int nb = 0; nb < Soln_Block_List.Nblk; nb++) {
    if (Soln_Block_List.Block[nb].used == ADAPTIVEBLOCK2D_USED) {
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
int dUdt_Multistage_Explicit(Electrostatic2D_Quad_Block *Soln_ptr,
                             AdaptiveBlock2D_List &Soln_Block_List,
                             Electrostatic2D_Input_Parameters &Input_Parameters,
   	                     const int I_Stage) {

  int error_flag;

  // Evaluate the stage solution residual for each solution block.
  for (int nb = 0; nb < Soln_Block_List.Nblk; nb++) {
    if (Soln_Block_List.Block[nb].used == ADAPTIVEBLOCK2D_USED) {
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
int Update_Solution_Multistage_Explicit(Electrostatic2D_Quad_Block *Soln_ptr,
                                        AdaptiveBlock2D_List &Soln_Block_List,
                                        Electrostatic2D_Input_Parameters &Input_Parameters,
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
 * Routine: Determine_Electric_Field                                  *
 **********************************************************************/
int Determine_Electric_Field(Electrostatic2D_Quad_Block *Soln_ptr,
			     AdaptiveBlock2D_List &Soln_Block_List,
			     Electrostatic2D_Input_Parameters &Input_Parameters) {

  int error_flag;

  for (int nb = 0; nb < Soln_Block_List.Nblk; nb++) {
    if (Soln_Block_List.Block[nb].used == ADAPTIVEBLOCK2D_USED) {
      error_flag = Determine_Electric_Field(Soln_ptr[nb],
					    Input_Parameters);
      if (error_flag) return error_flag;
    }
  }

  return 0;

}
