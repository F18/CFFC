/**********************************************************************
 * LevelSet2DQuadEikonalMultiBlock.cc                                 *
 *                                                                    *
 * Multi-block versions of subroutines for the solution of the 2D     *
 * Eikonal equations for the 2D Level Set multi-block quadrilateral   *
 * mesh solution classes.                                             *
 *                                                                    *
 **********************************************************************/

// Include 2D LevelSet quadrilateral mesh solution header file.

#ifndef _LEVELSET2D_QUAD_INCLUDED
#include "LevelSet2DQuad.h"
#endif // _LEVELSET2D_QUAD_INCLUDED

/**********************************************************************
 * LevelSet2D_Quad_Block -- Eikonal Multiple Block External           *
 *                          Subroutines.                              *
 **********************************************************************/

/**********************************************************************
 * Routine: Eikonal_Error                                             *
 *                                                                    *
 * This routine calculates the local area-weighted error in the       *
 * Eikonal equation solution.                                         *
 *                                                                    *
 **********************************************************************/
int Eikonal_Error(LevelSet2D_Quad_Block *Soln_ptr,
		  LevelSet2D_Input_Parameters &IP,
		  AdaptiveBlock2D_List &Soln_Block_List,
		  double &global_error,
		  double &global_area) {

  int error_flag;
  global_error = ZERO;
  global_area = ZERO;

  // Variables for each block:
  double block_error, block_area;

  // Compute the area-weighted error for each block.
  for (int nb = 0; nb < Soln_Block_List.Nblk; nb++) {
    if (Soln_Block_List.Block[nb].used == ADAPTIVEBLOCK2D_USED) {
      block_error = ZERO;
      block_area = ZERO;
      error_flag = Eikonal_Error(Soln_ptr[nb],
				 IP,
				 block_error,
				 block_area);
      if (error_flag) return error_flag;
      global_error = global_error + block_error;
      global_area = global_area + block_area;
    }
  }

  // Eikonal error calculation is successful.
  return 0;

}

/**********************************************************************
 * Routine: Explicit_Eikonal_Equation                                 *
 *                                                                    *
 * This routine manages the explicit solution of the Eikonal equation *
 * for iteratively forcing the level set function, psi, to be a       *
 * signed distance fuction.                                           *
 *                                                                    *
 **********************************************************************/
int Explicit_Eikonal_Equation(LevelSet2D_Quad_Block *Soln_ptr,
			      LevelSet2D_Input_Parameters &IP,
			      QuadTreeBlock_DataStructure &QuadTree,
			      AdaptiveBlockResourceList &GlobalSolnBlockList,
			      AdaptiveBlock2D_List &Soln_Block_List,
			      const int Redistance) {

  int error_flag;
  int ni = 0;                    // Iteration counter.
  int Number_of_Iterations = 0;  // Total number of iterations.
  double dTime;                  // Global time-step.
  double l2_norm = ONE;          // L2-norm of the residual.

  // Determine the number of redistancing steps.
  Number_of_Iterations = Redistance ? IP.Number_of_Redistance_Iterations :
                         IP.Number_of_Initial_Redistance_Iterations;

  // Exit immediately if iterative solution of the Eikonal equation is
  // not required.
  if (Number_of_Iterations < 1) return 0;

  // Store the solution before starting solution of the Eikonal equation.
  error_flag = Store_Initial_Eikonal_Solution(Soln_ptr,
					      Soln_Block_List);

  error_flag = CFFC_OR_MPI(error_flag);
  if (error_flag) return error_flag;

  // Determine the sign function before starting solution of the
  // Eikonal equation.
  error_flag = Calculate_Sign_Function(Soln_ptr,
				       Soln_Block_List,
				       IP);
  error_flag = CFFC_OR_MPI(error_flag);
  if (error_flag) return error_flag;

  // Conduct the explicit solution of the Eikonal equation for the
  // required number of iterations or until the l2-norm of the residual 
  // of the signed-distance function has been reduced by the desired
  // amount.
  while (l2_norm > IP.Redistance_Tolerance &&
	 ni < Number_of_Iterations) {

    // Determine the global time step.
    dTime = CFL_Eikonal(Soln_ptr,Soln_Block_List);
    // Find global minimum time step for all processors.
    dTime = CFFC_Minimum_MPI(dTime);
    // Set global time step.
    Set_Global_TimeStep(Soln_ptr,Soln_Block_List,dTime);

    // Update solution for next time step using a multistage time
    // stepping scheme.
    for (int i_stage = 1; i_stage <= IP.N_Stage; i_stage++) {

      // Step 1. Exchange solution information between neighbouring blocks.
      error_flag = Send_All_Messages(Soln_ptr,
				     Soln_Block_List,
				     NUM_VAR_LEVELSET2D,
				     OFF);
      error_flag = CFFC_OR_MPI(error_flag);
      if (error_flag) return error_flag;
	
      // Step 2. Apply boundary conditions for stage.
      BCs(Soln_ptr,Soln_Block_List,IP);

      // Step 3. Redetermine the sign function if required.
      if (IP.i_Eikonal_Sign_Function == EIKONAL_SIGN_FUNCTION_DERIVATIVE) {
	error_flag = Calculate_Sign_Function(Soln_ptr,
					     Soln_Block_List,
					     IP);
      }

      // Step 4. Determine solution residuals for stage.
      error_flag = dUdt_Multistage_Eikonal(Soln_ptr,
					   Soln_Block_List,
					   IP,
					   i_stage);
      //error_flag = CFFC_OR_MPI(error_flag);
      //if (error_flag) return error_flag;

      // Step 5. Update solution for stage.
      error_flag = Update_Multistage_Eikonal(Soln_ptr,
					     Soln_Block_List,
					     IP,
					     i_stage);
      //error_flag = CFFC_OR_MPI(error_flag);
      //if (error_flag) return error_flag;

    }

    // Determine l2-norm residual.
    l2_norm = L2_Norm_Residual(Soln_ptr,Soln_Block_List,1);
    l2_norm = sqr(l2_norm);
    l2_norm = CFFC_Summation_MPI(l2_norm);
    l2_norm = sqrt(l2_norm);

    // Update iteration counter.
    ni++;

  }

  // Apply boundary conditions.
  BCs(Soln_ptr,Soln_Block_List,IP);

  // Eikonal equation solved successfully.
  return 0;

}

/**********************************************************************
 * Routine: CFL_Eikonal                                               *
 *                                                                    *
 * Determines the allowable global and local time steps for the       *
 * solution of the Eikonal equation for a 1D array of 2D              *
 * quadrilateral multi-block solution blocks according to the         *
 * Courant-Friedrichs-Lewy condition.                                 *
 *                                                                    *
 **********************************************************************/
double CFL_Eikonal(LevelSet2D_Quad_Block *Soln_ptr,
		   AdaptiveBlock2D_List &Soln_Block_List) {

  double dtMin = MILLION;
  
  // Determine the allowable time step for each solution block.
  for (int nb = 0; nb < Soln_Block_List.Nblk; nb++) {
    if (Soln_Block_List.Block[nb].used == ADAPTIVEBLOCK2D_USED) {
      dtMin = min(dtMin,CFL_Eikonal(Soln_ptr[nb]));
    }
  }
  
  // Return the global time step.
  return dtMin;
  
}

/**********************************************************************
 * Routine: dUdt_Multistage_Eikonal                                   *
 *                                                                    *
 * This routine evaluates the stage solution residual for the Eikonal *
 * equation for a 1D array of 2D quadrilateral multi-block solution   *
 * blocks.  A variety of multistage explicit time integration and an  *
 * upwind finite-volume spatial discretization procedure can is used  *
 * depending on the specified input data.                             *
 *                                                                    *
 **********************************************************************/
int dUdt_Multistage_Eikonal(LevelSet2D_Quad_Block *Soln_ptr,
			    AdaptiveBlock2D_List &Soln_Block_List,
			    LevelSet2D_Input_Parameters &IP,
			    const int I_Stage) {

  int error_flag;

  // Update the solution for each solution block.
  for (int nb = 0; nb < Soln_Block_List.Nblk; nb++) {
    if (Soln_Block_List.Block[nb].used == ADAPTIVEBLOCK2D_USED) {
      error_flag = dUdt_Multistage_Eikonal(Soln_ptr[nb],
					   I_Stage,
					   IP);
      if (error_flag) return error_flag;
    }
  }

  // Residuals for each quadrilateral multi-block solution block 
  // successfully calculated.
  return 0;

}

/**********************************************************************
 * Routine: Update_Multistage_Eikonal                                 *
 *                                                                    *
 * This routine updates the solution for a 1D array of 2D             *
 * quadrilateral multi-block Level Set solution blocks.  Second-order *
 * multistage explicit time integration and a finite-volume spatial   *
 * discretization procedure is used.                                  *
 *                                                                    *
 **********************************************************************/
int Update_Multistage_Eikonal(LevelSet2D_Quad_Block *Soln_ptr,
			      AdaptiveBlock2D_List &Soln_Block_List,
			      LevelSet2D_Input_Parameters &IP,
			      const int I_Stage) {

  int error_flag;

  // Update the solution for each solution block.
  for (int nb = 0; nb < Soln_Block_List.Nblk; nb++) {
    if (Soln_Block_List.Block[nb].used == ADAPTIVEBLOCK2D_USED) {
      error_flag = Update_Multistage_Eikonal(Soln_ptr[nb],
					     I_Stage,
					     IP);
      if (error_flag) return error_flag;
    }
  }

  // Quadrilateral multi-block solution blocks successfully updated.
  return 0;

}
