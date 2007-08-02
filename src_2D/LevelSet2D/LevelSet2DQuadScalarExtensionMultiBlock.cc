/**********************************************************************
 * LevelSet2DQuadScalarExtensionMultiBlock.cc                         *
 *                                                                    *
 * Scalar extension multi-block versions of subroutines for 2D Level  *
 * Set multi-block quadrilateral mesh solution classes.               *
 *                                                                    *
 **********************************************************************/

// Include 2D Level Set quadrilateral mesh solution header file.

#ifndef _LEVELSET2D_QUAD_INCLUDED
#include "LevelSet2DQuad.h"
#endif // _LEVELSET2D_QUAD_INCLUDED

/**********************************************************************
 * LevelSet2D_Quad_Block -- Scalar Extension Multiple Block External  *
 *                          Subroutines.                              *
 **********************************************************************/

/**********************************************************************
 * Routine: Explicit_Scalar_Extension_Equation                        *
 *                                                                    *
 * This routine manages the explicit solution of the scalar (front    *
 * speed) extension equation which extends any scalar defined on an   *
 * arbitrary interface to the rest of the flow domain in rays normal  *
 * to the front.                                                      *
 *                                                                    *
 **********************************************************************/
int Explicit_Scalar_Extension_Equation(LevelSet2D_Quad_Block *Soln_ptr,
				       LevelSet2D_Input_Parameters &Input_Parameters,
				       QuadTreeBlock_DataStructure &QuadTree,
				       AdaptiveBlockResourceList &GlobalSolnBlockList,
				       AdaptiveBlock2D_List &Soln_Block_List) {

  int error_flag;
  int ni = 0;           // Iteration counter.
  double dTime;         // Global time-step.
  double l2_norm = ONE; // L2-norm of the residual.

  // Exit immediately if extension of scalar quantities is not required.
  if (Input_Parameters.Number_of_Scalar_Extension_Iterations < 1) return 0;

  // Determine the sign function before starting solution of the
  // scalar extension equation.
  error_flag = Calculate_Sign_Function(Soln_ptr,
				       Soln_Block_List,
				       Input_Parameters);
  error_flag = CFDkit_OR_MPI(error_flag);
  if (error_flag) return error_flag;

  // Conduct explicit solution of the scalar extension equation for the
  // required number of iterations or until the l2-norm of the residual 
  // of the scalar quanitity has been reduced by the desired amount.
  while (l2_norm > NANO && ni < Input_Parameters.Number_of_Scalar_Extension_Iterations) {

    // Determine global time step.
    dTime = CFL_Scalar_Extension(Soln_ptr,Soln_Block_List);
    // Find global minimum time step for all processors.
    dTime = CFDkit_Minimum_MPI(dTime);
    // Set global time step.
    Set_Global_TimeStep(Soln_ptr,Soln_Block_List,dTime);

    // Update solution for next time step using a multistage time
    // stepping scheme.
    for (int i_stage = 1; i_stage <= Input_Parameters.N_Stage; i_stage++) {

      // Step 1. Exchange solution information between neighbouring blocks.
      error_flag = Send_All_Messages(Soln_ptr,
				     Soln_Block_List,
				     NUM_VAR_LEVELSET2D,
				     OFF);
      error_flag = CFDkit_OR_MPI(error_flag);
      if (error_flag) return error_flag;

      // Step 2. Apply boundary conditions for stage.
      BCs(Soln_ptr,Soln_Block_List,Input_Parameters);

      // Step 3. Determine solution residuals for stage.
      error_flag = dUdt_Multistage_Scalar_Extension(Soln_ptr,
						    Soln_Block_List,
						    Input_Parameters,
						    i_stage);
      error_flag = CFDkit_OR_MPI(error_flag);
      if (error_flag) return error_flag;

      // Step 4. Update solution for stage.
      error_flag = Update_Multistage_Scalar_Extension(Soln_ptr,
						      Soln_Block_List,
						      Input_Parameters,
						      i_stage);
      error_flag = CFDkit_OR_MPI(error_flag);
      if (error_flag) return error_flag;

    }

    // Determine L2-norm residual.
    l2_norm = L2_Norm_Residual(Soln_ptr,Soln_Block_List,2);
    l2_norm = sqr(l2_norm);
    l2_norm = CFDkit_Summation_MPI(l2_norm);
    l2_norm = sqrt(l2_norm);

    // Update iteration counter.
    ni++;

  }

  // Apply boundary conditions.
  BCs(Soln_ptr,Soln_Block_List,Input_Parameters);

  // Scalar extension equation solved successfully.
  return 0;

}

/**********************************************************************
 * Routine: dUdt_Multistage_Scalar_Extension                          *
 *                                                                    *
 * This routine evaluates the stage solution residual for the scalar  *
 * (front speed) extension equation for a 1D array of 2D              *
 * quadrilateral multi-block solution blocks.  A variety of           *
 * multistage explicit time integration and an upwind finite-volume   *
 * spatial discretization procedure can is used depending on the      *
 * specified input data.                                              *
 *                                                                    *
 **********************************************************************/
int dUdt_Multistage_Scalar_Extension(LevelSet2D_Quad_Block *Soln_ptr,
				     AdaptiveBlock2D_List &Soln_Block_List,
				     LevelSet2D_Input_Parameters &Input_Parameters,
				     const int I_Stage) {

  int error_flag;

  // Update the solution for each solution block.
  for (int nb = 0; nb < Soln_Block_List.Nblk; nb++) {
    if (Soln_Block_List.Block[nb].used == ADAPTIVEBLOCK2D_USED) {
      error_flag = dUdt_Multistage_Scalar_Extension(Soln_ptr[nb],
						    I_Stage,
						    Input_Parameters);
      if (error_flag) return error_flag;
    }
  }
  
  // Residuals for each quadrilateral multi-block solution block 
  // successfully calculated.
  return 0;

}

/**********************************************************************
 * Routine: Update_Multistage_Scalar_Extension                        *
 *                                                                    *
 * This routine updates the solution for a 1D array of 2D             *
 * quadrilateral multi-block solution blocks for the scalar (front    *
 * speed) extension equation.  Second-order multistage explicit time  *
 * integration and a finite-volume spatial discretization procedure   *
 * is used.                                                           *
 *                                                                    *
 **********************************************************************/
int Update_Multistage_Scalar_Extension(LevelSet2D_Quad_Block *Soln_ptr,
				       AdaptiveBlock2D_List &Soln_Block_List,
				       LevelSet2D_Input_Parameters &Input_Parameters,
				       const int I_Stage) {

  int error_flag;

  // Update the solution for each solution block.
  for (int nb = 0; nb < Soln_Block_List.Nblk; nb++) {
    if (Soln_Block_List.Block[nb].used == ADAPTIVEBLOCK2D_USED) {
      error_flag = Update_Multistage_Scalar_Extension(Soln_ptr[nb],
						      I_Stage,
						      Input_Parameters);
      if (error_flag) return error_flag;
    }
  }
  
  // Quadrilateral multi-block solution blocks successfully updated.
  return 0;
 
}

/**********************************************************************
 * Routine: CFL_Scalar_Extension                                      *
 *                                                                    *
 * Determines the allowable global and local time steps for the       *
 * solution of the scalar (front speed) extension equation for a 1D   *
 * array of 2D quadrilateral multi-block solution blocks according to *
 * the Courant-Friedrichs-Lewy condition.                             *
 *                                                                    *
 **********************************************************************/
double CFL_Scalar_Extension(LevelSet2D_Quad_Block *Soln_ptr,
			    AdaptiveBlock2D_List &Soln_Block_List) {

  double dtMin = MILLION;
  
  // Determine the allowable time step for each solution block.
  for (int nb = 0; nb < Soln_Block_List.Nblk; nb++)
    if (Soln_Block_List.Block[nb].used == ADAPTIVEBLOCK2D_USED)
      dtMin = min(dtMin,CFL_Scalar_Extension(Soln_ptr[nb]));
  
  // Return the global time step.
  return dtMin;
  
}
