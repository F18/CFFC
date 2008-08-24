/*! \file Euler2DQuad.cc
  @brief Subroutines for 2D Euler Quadrilateral Mesh Solution Class. */

/* Include required C++ libraries. */
// None

/* Using std namespace functions */
// None

/* Include CFFC header files */
#include "Euler2DQuad.h"                // Euler2D_Quad_Block class

/**********************************************************************
 * Euler2D_Quad_Block -- Create storage for the static variables.     *
 **********************************************************************/
// Initialize residual_variable
int Euler2D_Quad_Block::residual_variable = 1;
// Initialize Number_of_Residual_Norms
int Euler2D_Quad_Block::Number_of_Residual_Norms = 1;
// Initialize Flow_Type
int Euler2D_Quad_Block::Flow_Type = FLOWTYPE_INVISCID;
// Initialize RefU
Euler2D_pState Euler2D_Quad_Block::RefU(1.0);
// Initialize ExactSoln
Euler2D_ExactSolutions *Euler2D_Quad_Block::ExactSoln = NULL;

/***********************************************************************
 * Euler2D_Quad_Block -- Single Block Member Functions.                *
 **********************************************************************/

/*****************************************************//**
 * Copy the solution information of quadrilateral solution 
 * block SolnBlk to the current solution block.
 ********************************************************/
Euler2D_Quad_Block & Euler2D_Quad_Block::operator =(const Euler2D_Quad_Block &Soln){
  
  int i, j, k;

  // Handle self-assignment:
  if (this == & Soln) return *this;

  // check if solution block Soln has memory allocated
  if (Soln.U != NULL){
    /* Allocate (re-allocate) memory for the solution
       of the quadrilateral solution block as necessary. */
    allocate(Soln.NCi-2*Soln.Nghost,
	     Soln.NCj-2*Soln.Nghost,
	     Soln.Nghost);

    /* Set the same number of high-order objects
       as that of the rhs block. */
    allocate_HighOrder_Array(Soln.NumberOfHighOrderVariables);

  } else {
    deallocate();
  }

  /* Set the axisymmetric/planar flow indicator. */
  Axisymmetric = Soln.Axisymmetric;
  
  /* Copy the grid. */
  Grid = Soln.Grid;
  
  /* Copy the solution information from Soln. */
  if (Soln.U != NULL) {
    for ( j  = JCl-Nghost ; j <= JCu+Nghost ; ++j ) {
      for ( i = ICl-Nghost ; i <= ICu+Nghost ; ++i ) {
	U[i][j] = Soln.U[i][j];
	W[i][j] = Soln.W[i][j];
	for ( k = 0 ; k <= NUMBER_OF_RESIDUAL_VECTORS_EULER2D-1 ; ++k ) {
	  dUdt[i][j][k] = Soln.dUdt[i][j][k];
	} /* endfor */
	dWdx[i][j] = Soln.dWdx[i][j];
	dWdy[i][j] = Soln.dWdy[i][j];
	phi[i][j] = Soln.phi[i][j];
	Uo[i][j] = Soln.Uo[i][j];
	dt[i][j] = Soln.dt[i][j];
      } /* endfor */
    } /* endfor */

    for (j  = JCl-Nghost ; j <= JCu+Nghost ; ++j ) {
      WoW[j] = Soln.WoW[j];
      WoE[j] = Soln.WoE[j];
    }/* endfor */
    
    for ( i = ICl-Nghost ; i <= ICu+Nghost ; ++i ) {
      WoS[i] = Soln.WoS[i];
      WoN[i] = Soln.WoN[i];
    }/* endfor */

    // allocate memory for high-order boundary conditions.
    allocate_HighOrder_BoundaryConditions();

    for (j  = JCl ; j <= JCu ; ++j ) {
      // Copy West high-order BCs
      if (HO_WoW != NULL){
	HO_WoW[j] = Soln.HO_WoW[j];
      }

      // Copy East high-order BCs
      if (HO_WoE != NULL){
	HO_WoE[j] = Soln.HO_WoE[j];
      }
    }

    for ( i = ICl ; i <= ICu ; ++i ) {
      // Copy South high-order BCs
      if (HO_WoS != NULL){
	HO_WoS[i] = Soln.HO_WoS[i];
      }
      
      // Copy North high-order BCs
      if (HO_WoN != NULL){
	HO_WoN[i] = Soln.HO_WoN[i];
      }
    }
    
    // Copy the high-order objects
    for (k = 1; k <= NumberOfHighOrderVariables; ++k){
      HighOrderVariable(k-1) = Soln.HighOrderVariable(k-1);
    }/* endfor */
    
  }/* endif */

  // Copy boundary reference states
  Ref_State_BC_North = Soln.Ref_State_BC_North;
  Ref_State_BC_South = Soln.Ref_State_BC_South;
  Ref_State_BC_East = Soln.Ref_State_BC_East;
  Ref_State_BC_West = Soln.Ref_State_BC_West;  

  // Reset accuracy assessment flag
  AssessAccuracy.ResetForNewCalculation();

  return *this;
}

/***************************************************//**
 * Assigns boundary condition reference data based on
 * the boundary type for the specified quadrilateral 
 * solution block based on the input parameters.
 * This routine makes the link between user's specifications
 * and the values that are set as boundary reference states.
 ********************************************************/
void Euler2D_Quad_Block::Set_Boundary_Reference_States_Based_On_Input(const Euler2D_Input_Parameters &IP){

  // Set the reference values for boundary reference states
  Set_Reference_Values_For_Boundary_States(IP.Ref_State_BC_North,
					   IP.Ref_State_BC_South,
					   IP.Ref_State_BC_East,
					   IP.Ref_State_BC_West);
}

/***************************************************//**
 * Assigns reference values for the boundary condition 
 * reference states to what the user specified.
 * These values are used in Set_Boundary_Reference_States()
 * routine.
 ********************************************************/
void Euler2D_Quad_Block::Set_Reference_Values_For_Boundary_States(const Euler2D_pState & Ref_North,
								  const Euler2D_pState & Ref_South,
								  const Euler2D_pState & Ref_East,
								  const Euler2D_pState & Ref_West){
  Ref_State_BC_North = Ref_North;
  Ref_State_BC_South = Ref_South;
  Ref_State_BC_East = Ref_East;
  Ref_State_BC_West = Ref_West;
}

