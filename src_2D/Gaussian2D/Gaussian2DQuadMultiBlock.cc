/* Gaussian2DQuadMultiBlock.cc:  Multi-Block Versions of Subroutines for 2D Gaussian
                                 Multi-Block Quadrilateral Mesh 
                                 Solution Classes. */

/* Include 2D Gaussian quadrilateral mesh solution header file. */

#ifndef _GAUSSIAN2D_QUAD_INCLUDED
#include "Gaussian2DQuad.h"
#endif // _GAUSSIAN2D_QUAD_INCLUDED

/**************************************************************************
 * Gaussian2D_Quad_Block -- Multiple Block External Subroutines.          *
 **************************************************************************/

/********************************************************
 * Routine: Allocate                                    *
 *                                                      *
 * Allocate memory for 1D array of 2D quadrilateral     *
 * multi-block solution blocks.                         *
 *                                                      *
 ********************************************************/
Gaussian2D_Quad_Block* Allocate(Gaussian2D_Quad_Block *Soln_ptr,
                               Gaussian2D_Input_Parameters &Input_Parameters) {

    /* Allocate memory. */

    Soln_ptr = new Gaussian2D_Quad_Block[Input_Parameters.Number_of_Blocks_Per_Processor];

    // Assign the residual variable.
    Soln_ptr[0].residual_variable = Input_Parameters.i_Residual_Variable;

    // Set heat transfer flag.
    Soln_ptr[0].Heat_Transfer = Input_Parameters.Heat_Transfer;

    /* Return memory location. */

    return(Soln_ptr);

}

/********************************************************
 * Routine: Deallocate                                  *
 *                                                      *
 * Deallocate memory for 1D array of 2D quadrilateral   *
 * multi-block solution blocks.                         *
 *                                                      *
 ********************************************************/
Gaussian2D_Quad_Block* Deallocate(Gaussian2D_Quad_Block *Soln_ptr,
                                  Gaussian2D_Input_Parameters &Input_Parameters) {

    int i;
 
    /* Deallocate memory. */

    for ( i = 0 ; i <= Input_Parameters.Number_of_Blocks_Per_Processor-1 ; ++i ) {
       if (Soln_ptr[i].W != NULL) Soln_ptr[i].deallocate();
    }  /* endfor */
    delete []Soln_ptr;
    Soln_ptr = NULL;

    /* Return memory location. */

    return(Soln_ptr);

}

/**********************************************************************
 * Routine: CreateInitialSolutionBlocks                               *
 *                                                                    *
 * This routine calls the Create_Initial_Solution_Blocks routine      *
 * found in AMR.h.                                                    *
 *                                                                    *
 **********************************************************************/
Gaussian2D_Quad_Block* CreateInitialSolutionBlocks(Grid2D_Quad_Block **InitMeshBlk,
						   Gaussian2D_Quad_Block *Soln_ptr,
						   Gaussian2D_Input_Parameters &Input_Parameters,
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

  // Return the solution blocks.
  return Soln_ptr;

}

/********************************************************
 * Routine: ICs                                         *
 *                                                      *
 * Assigns initial conditions and data to the           *
 * solution variables of a 1D array of 2D quadrilateral *
 * multi-block solution blocks.                         *
 *                                                      *
 ********************************************************/
void ICs(Gaussian2D_Quad_Block *Soln_ptr,
         AdaptiveBlock2D_List &Soln_Block_List,
         Gaussian2D_Input_Parameters &Input_Parameters) {

    int i;
    Gaussian2D_pState Wo[5];

    /* Assign the gas constants for the gas of interest. */

    Input_Parameters.Wo.setgas(Input_Parameters.Gas_Type);
    Input_Parameters.Uo.setgas(Input_Parameters.Gas_Type);
    Input_Parameters.Uo = U(Input_Parameters.Wo);

    /* Define various reference flow states. */

    Wo[0] = Input_Parameters.Wo;
    Wo[1] = Input_Parameters.W1;
    Wo[2] = Input_Parameters.W2;

    /* Assign initial data for each solution block. */

    for ( i = 0 ; i <= Soln_Block_List.Nblk-1 ; ++i ) {
       if (Soln_Block_List.Block[i].used == ADAPTIVEBLOCK2D_USED) {
	  // Set flow geometry indicator (planar/axisymmetric) flow block.
 	  Soln_ptr[i].Axisymmetric = Input_Parameters.Axisymmetric;

          // Set initial data.
          ICs(Soln_ptr[i], 
              Input_Parameters,
              Wo);
       } /* endif */
    }  /* endfor */

}

/********************************************************
 * Routine: BCs                                         *
 *                                                      *
 * Apply boundary conditions at boundaries of a 1D      *
 * array of 2D quadrilateral multi-block solution       *
 * blocks.                                              *
 *                                                      *
 ********************************************************/
void BCs(Gaussian2D_Quad_Block *Soln_ptr,
         AdaptiveBlock2D_List &Soln_Block_List,
	 Gaussian2D_Input_Parameters &IP) {

    int i;

    /* Prescribe boundary data for each solution block. */

    for ( i = 0 ; i <= Soln_Block_List.Nblk-1 ; ++i ) {
       if (Soln_Block_List.Block[i].used == ADAPTIVEBLOCK2D_USED) {
	 BCs(Soln_ptr[i],IP);
       } /* endif */
    }  /* endfor */

}

/********************************************************
 * Routine: CFL                                         *
 *                                                      *
 * Determines the allowable global and local time steps *
 * (for explicit Euler time stepping scheme) for a 1D   *
 * array of 2D quadrilateral multi-block solution       *
 * blocks according to the Courant-Friedrichs-Lewy      *
 * condition.                                           *
 *                                                      *
 ********************************************************/
double CFL(Gaussian2D_Quad_Block *Soln_ptr,
           AdaptiveBlock2D_List &Soln_Block_List,
           Gaussian2D_Input_Parameters &Input_Parameters) {

    int i;
    double dtMin;

    dtMin = MILLION;

    /* Determine the allowable time step for each solution block. */

    for ( i = 0 ; i <= Soln_Block_List.Nblk-1 ; ++i ) {
       if (Soln_Block_List.Block[i].used == ADAPTIVEBLOCK2D_USED) {
          dtMin = min(dtMin, CFL(Soln_ptr[i]));
       } /* endif */
    }  /* endfor */

    /* Return the global time step. */

    return (dtMin);

}

/********************************************************
 * Routine: Set_Global_TimeStep                         *
 *                                                      *
 * Assigns global time step to a 1D array of 2D         *
 * quadrilateral multi-block solution blocks for        *
 * time-accurate calculations.                          *
 *                                                      *
 ********************************************************/
void Set_Global_TimeStep(Gaussian2D_Quad_Block *Soln_ptr,
                         AdaptiveBlock2D_List &Soln_Block_List,
                         const double &Dt_min) {

    int i;

    for ( i = 0 ; i <= Soln_Block_List.Nblk-1 ; ++i ) {
       if (Soln_Block_List.Block[i].used == ADAPTIVEBLOCK2D_USED) {
          Set_Global_TimeStep(Soln_ptr[i], Dt_min);
       } /* endif */
    }  /* endfor */

}

/********************************************************
 * Routine: L1_Norm_Residual                            *
 *                                                      *
 * Determines the L1-norm of the solution residual for  *
 * a 1D array of 2D quadrilateral multi-block solution  *
 * blocks.  Useful for monitoring convergence of the    *
 * solution for steady state problems.                  *
 *                                                      *
 ********************************************************/
double L1_Norm_Residual(Gaussian2D_Quad_Block *Soln_ptr,
                        AdaptiveBlock2D_List &Soln_Block_List) {

    int i;
    double l1_norm;

    l1_norm = ZERO;

    /* Calculate the L1-norm.
       Sum the L1-norm for each solution block. */

    for ( i = 0 ; i <= Soln_Block_List.Nblk-1 ; ++i ) {
       if (Soln_Block_List.Block[i].used == ADAPTIVEBLOCK2D_USED) {
          l1_norm += L1_Norm_Residual(Soln_ptr[i]);
       } /* endif */
    }  /* endfor */

    /* Return the L1-norm. */

    return (l1_norm);

}

/********************************************************
 * Routine: L2_Norm_Residual                            *
 *                                                      *
 * Determines the L2-norm of the solution residual for  *
 * a 1D array of 2D quadrilateral multi-block solution  *
 * blocks.  Useful for monitoring convergence of the    *
 * solution for steady state problems.                  *
 *                                                      *
 ********************************************************/
double L2_Norm_Residual(Gaussian2D_Quad_Block *Soln_ptr,
                        AdaptiveBlock2D_List &Soln_Block_List) {

    int i;
    double l2_norm;

    l2_norm = ZERO;

    /* Sum the square of the L2-norm for each solution block. */

    for ( i = 0 ; i <= Soln_Block_List.Nblk-1 ; ++i ) {
       if (Soln_Block_List.Block[i].used == ADAPTIVEBLOCK2D_USED) {
          l2_norm += sqr(L2_Norm_Residual(Soln_ptr[i]));
       } /* endif */
    }  /* endfor */

    /* Calculate the L2-norm for all blocks. */

    l2_norm = sqrt(l2_norm);

    /* Return the L2-norm. */

    return (l2_norm);

}

/********************************************************
 * Routine: Max_Norm_Residual                           *
 *                                                      *
 * Determines the maximum norm of the solution residual *
 * for a 1D array of 2D quadrilateral multi-block       *
 * solution blocks.  Useful for monitoring convergence  *
 * of the solution for steady state problems.           *
 *                                                      *
 ********************************************************/
double Max_Norm_Residual(Gaussian2D_Quad_Block *Soln_ptr,
                         AdaptiveBlock2D_List &Soln_Block_List) {

    int i;
    double max_norm;

    max_norm = ZERO;

    /* Find the maximum norm for all solution blocks. */

    for ( i = 0 ; i <= Soln_Block_List.Nblk-1 ; ++i ) {
       if (Soln_Block_List.Block[i].used == ADAPTIVEBLOCK2D_USED) {
          max_norm = max(max_norm, Max_Norm_Residual(Soln_ptr[i]));
       } /* endif */
    }  /* endfor */

    /* Return the maximum norm. */

    return (max_norm);

}

/********************************************************
 * Routine: Evaluate_Limiters                           *
 *                                                      *
 * Set conditions to evaluate the limiters for a        *
 * 1D array of 2D quadrilateral multi-block solution    *
 * blocks.                                              *
 *                                                      *
 ********************************************************/
void Evaluate_Limiters(Gaussian2D_Quad_Block *Soln_ptr,
                       AdaptiveBlock2D_List &Soln_Block_List) {

    int i;
    for ( i = 0 ; i <= Soln_Block_List.Nblk-1 ; ++i ) {
       if (Soln_Block_List.Block[i].used == ADAPTIVEBLOCK2D_USED) Soln_ptr[i].evaluate_limiters();
    }  /* endfor */

}

/********************************************************
 * Routine: Freeze_Limiters                             *
 *                                                      *
 * Set conditions to freeze the limiters for a          *
 * 1D array of 2D quadrilateral multi-block solution    *
 * blocks.                                              *
 *                                                      *
 ********************************************************/
void Freeze_Limiters(Gaussian2D_Quad_Block *Soln_ptr,
                     AdaptiveBlock2D_List &Soln_Block_List) {

    int i;
    for ( i = 0 ; i <= Soln_Block_List.Nblk-1 ; ++i ) {
       if (Soln_Block_List.Block[i].used == ADAPTIVEBLOCK2D_USED) Soln_ptr[i].freeze_limiters();
    }  /* endfor */

}

/****************************************************************
 * Routine: Apply_Boundary_Flux_Correction                      *
 *                                                              *
 * Apply flux corrections at boundaries of a 1D array           *
 * of 2D quadrilateral multi-block solution blocks to           *
 * ensure that the scheme is conservative at boundaries         *
 * with resolution mesh changes.                                *
 *                                                              *
 ****************************************************************/
void Apply_Boundary_Flux_Corrections(Gaussian2D_Quad_Block *Soln_ptr,
                                     AdaptiveBlock2D_List &Soln_Block_List) {

    int i;

    /* Prescribe boundary data for each solution block. */

    for ( i = 0 ; i <= Soln_Block_List.Nblk-1 ; ++i ) {
       if (Soln_Block_List.Block[i].used == ADAPTIVEBLOCK2D_USED) {
          Apply_Boundary_Flux_Corrections(Soln_ptr[i],
                                          Soln_Block_List.Block[i].nN,
                                          Soln_Block_List.Block[i].nS,
                                          Soln_Block_List.Block[i].nE,
                                          Soln_Block_List.Block[i].nW);
       } /* endif */
    }  /* endfor */

}

/****************************************************************
 * Routine: Apply_Boundary_Flux_Corrections_Multistage_Explicit *
 *                                                              *
 * Apply flux corrections at boundaries of a 1D array           *
 * of 2D quadrilateral multi-block solution blocks to           *
 * ensure that the scheme is conservative at boundaries         *
 * with resolution mesh changes.                                *
 *                                                              *
 ****************************************************************/
void Apply_Boundary_Flux_Corrections_Multistage_Explicit(Gaussian2D_Quad_Block *Soln_ptr,
                                                         AdaptiveBlock2D_List &Soln_Block_List,
                                                         Gaussian2D_Input_Parameters &Input_Parameters,
   	                                                 const int I_Stage) {

    int i;

    /* Prescribe boundary data for each solution block. */

    for ( i = 0 ; i <= Soln_Block_List.Nblk-1 ; ++i ) {
       if (Soln_Block_List.Block[i].used == ADAPTIVEBLOCK2D_USED) {
          Apply_Boundary_Flux_Corrections_Multistage_Explicit(Soln_ptr[i],
                                                              I_Stage,
                                                              Input_Parameters.N_Stage,
                                                              Input_Parameters.CFL_Number,
                                                              Input_Parameters.i_Time_Integration,
                                                              Input_Parameters.i_Reconstruction,
                                                              Input_Parameters.i_Limiter, 
                                                              Input_Parameters.i_Flux_Function,
                                                              Soln_Block_List.Block[i].nN,
                                                              Soln_Block_List.Block[i].nS,
                                                              Soln_Block_List.Block[i].nE,
                                                              Soln_Block_List.Block[i].nW);
       } /* endif */
    }  /* endfor */

}

/********************************************************
 * Routine: dUdt_Multistage_Explicit                    *
 *                                                      *
 * This routine evaluates the stage solution residual   *
 * for a 1D array of 2D quadrilateral multi-block       *
 * solution blocks.  A variety of multistage explicit   *
 * time integration and upwind finite-volume spatial    *
 * discretization procedures can be used depending on   *
 * the specified input values.                          *
 *                                                      *
 ********************************************************/
int dUdt_Multistage_Explicit(Gaussian2D_Quad_Block *Soln_ptr,
                             AdaptiveBlock2D_List &Soln_Block_List,
                             Gaussian2D_Input_Parameters &Input_Parameters,
   	                     const int I_Stage) {

    int i, error_flag;

    error_flag = 0;

    /* Evaluate the stage solution residual for each solution block. */

    for ( i = 0 ; i <= Soln_Block_List.Nblk-1 ; ++i ) {
       if (Soln_Block_List.Block[i].used == ADAPTIVEBLOCK2D_USED) {

          error_flag = dUdt_Multistage_Explicit(Soln_ptr[i],
                                                I_Stage,
                                                Input_Parameters);
          if (error_flag) return (error_flag);
       } /* endif */
    }  /* endfor */

    /* Residuals for each quadrilateral multi-block solution block
       successfully calcualted.  Return. */

    return(error_flag);

}

/********************************************************
 * Routine: Update_Solution_Multistage_Explicit         *
 *                                                      *
 * This routine updates the solution for a 1D array of  *
 * 2D quadrilateral multi-block solution blocks.  A     *
 * variety of multistage explicit time integration      *
 * and upwind finite-volume spatial discretization      *
 * procedures can be used depending on the specified    *
 * input values.                                        *
 *                                                      *
 ********************************************************/
int Update_Solution_Multistage_Explicit(Gaussian2D_Quad_Block *Soln_ptr,
                                        AdaptiveBlock2D_List &Soln_Block_List,
                                        Gaussian2D_Input_Parameters &Input_Parameters,
   	                                const int I_Stage) {

    int i, error_flag;

    error_flag = 0;

    /* Update the solution for each solution block. */

    for ( i = 0 ; i <= Soln_Block_List.Nblk-1 ; ++i ) {
       if (Soln_Block_List.Block[i].used == ADAPTIVEBLOCK2D_USED) {
          error_flag = Update_Solution_Multistage_Explicit(Soln_ptr[i],
                                                           I_Stage,
							   Input_Parameters);
          if (error_flag) return (error_flag);
       } /* endif */
    }  /* endfor */

    /* Quadrilateral multi-block solution blocks
       successfully updated.  Return. */

    return(error_flag);

}


/********************************************************
 * Routine: Ramp_up_Reference_Mach_Number               *
 *                                                      *
 * This routine ramps up the reference mach number      *
 *                                                      *
 ********************************************************/

void Ramp_up_Reference_Mach_Number(Gaussian2D_Quad_Block *SolnBlk,
				   AdaptiveBlock2D_List &Soln_Block_List,
				   Gaussian2D_Input_Parameters &Input_Parameters,
				   const int Number_of_Time_Steps){
  int i;

    /* Update the Reference for each solution block. */

    for ( i = 0 ; i <= Soln_Block_List.Nblk-1 ; ++i ) {
       if (Soln_Block_List.Block[i].used == ADAPTIVEBLOCK2D_USED) {
	 Ramp_up_Reference_Mach_Number(SolnBlk[i],
				       Input_Parameters,
				       Number_of_Time_Steps);          
       } /* endif */
    }  /* endfor */

}




