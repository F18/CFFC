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
// Initialize RefW
Euler2D_pState Euler2D_Quad_Block::RefW(1.0);
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
    
  }/* endif */

  // Copy high-order objects
  copy_HighOrder_Objects(Soln);

  // Copy boundary reference states
  Ref_State_BC_North = Soln.Ref_State_BC_North;
  Ref_State_BC_South = Soln.Ref_State_BC_South;
  Ref_State_BC_East = Soln.Ref_State_BC_East;
  Ref_State_BC_West = Soln.Ref_State_BC_West;  

  // Reset accuracy assessment flag
  AssessAccuracy.ResetForNewCalculation();

  return *this;
}

/*!
 * Copy the high-order objects from the SolnBlk.
 */
void Euler2D_Quad_Block::copy_HighOrder_Objects(const Euler2D_Quad_Block &SolnBlk){

  int i, j, k;

  if (SolnBlk.U != NULL){

    /* Set the same number of high-order objects
       as that of the rhs block. */
    allocate_HighOrder_Array(SolnBlk.NumberOfHighOrderVariables);
    
    // Copy the high-order objects
    for (k = 1; k <= NumberOfHighOrderVariables; ++k){
      // set geometry pointer for each high-order object to the current grid
      HighOrderVariable(k-1).SetGeometryPointer(Grid);
      // copy high-order content from the the SolnBlk
      HighOrderVariable(k-1) = SolnBlk.HighOrderVariable(k-1);
    }/* endfor */

    // Copy the high-order boundary conditions.

    // === North BCs
    if (SolnBlk.HO_WoN != NULL){

      if (HO_WoN != NULL){
	// deallocate memory
	delete [] HO_WoN; HO_WoN = NULL;
      }

      // allocate new memory based on the number of the current grid
      HO_WoN = new BC_Type[NCi];

      for (i=0; i<NCi; ++i){
	// allocate BC memory for each constrained Gauss quadrature point
	BC_NorthCell(i).InitializeCauchyBCs(Grid.NumOfConstrainedGaussQuadPoints_North(i,JCu),
					    Grid.BCtypeN[i]);

	// Copy North high-order BCs
	HO_WoN[i] = SolnBlk.HO_WoN[i];
      }
      
    } else if ( HO_WoN != NULL){
      // deallocate memory
      delete [] HO_WoN; HO_WoN = NULL;
    }


    // === South BCs
    if (SolnBlk.HO_WoS != NULL){

      if (HO_WoS != NULL){
	// deallocate memory
	delete [] HO_WoS; HO_WoS = NULL;
      }

      // allocate new memory    
      HO_WoS = new BC_Type[NCi];

      for (i=0; i<NCi; ++i){
	// allocate BC memory for each constrained Gauss quadrature point
	BC_SouthCell(i).InitializeCauchyBCs(Grid.NumOfConstrainedGaussQuadPoints_South(i,JCl),
					    Grid.BCtypeS[i]);
	// Copy South high-order BCs
	HO_WoS[i] = SolnBlk.HO_WoS[i];
      }

    } else if (HO_WoS != NULL){
      // deallocate memory
      delete [] HO_WoS; HO_WoS = NULL;
    }


    // === East BCs
    if (SolnBlk.HO_WoE != NULL){

      if (HO_WoE != NULL){
	// deallocate memory
	delete [] HO_WoE; HO_WoE = NULL;
      }

      // allocate new memory    
      HO_WoE = new BC_Type[NCj];

      for (j=0; j<NCj; ++j){
	// allocate BC memory for each constrained Gauss quadrature point
	BC_EastCell(j).InitializeCauchyBCs(Grid.NumOfConstrainedGaussQuadPoints_East(ICu,j),
					   Grid.BCtypeE[j]);
	// Copy East high-order BCs
	HO_WoE[j] = SolnBlk.HO_WoE[j];
      }

    } else if (HO_WoE != NULL){
      // deallocate memory
      delete [] HO_WoE; HO_WoE = NULL;
    }

    // === West BCs
    if (SolnBlk.HO_WoW != NULL){

      if (HO_WoW != NULL){
	// deallocate memory
	delete [] HO_WoW; HO_WoW = NULL;
      }

      // allocate new memory    
      HO_WoW = new BC_Type[NCj];

      for (j=0; j<NCj; ++j){
	// allocate BC memory for each constrained Gauss quadrature point
	BC_WestCell(j).InitializeCauchyBCs(Grid.NumOfConstrainedGaussQuadPoints_West(ICl,j),
					   Grid.BCtypeW[j]);
	// Copy West high-order BCs
	HO_WoW[j] = SolnBlk.HO_WoW[j];
      }

    } else if (HO_WoW != NULL){
      // deallocate memory
      delete [] HO_WoW; HO_WoW = NULL;
    }
    
  } // endif
  
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

/*****************************************************************************//**
 * Euler2D_Quad_Block::LoadSendBuffer_C2F -- Loads send message buffer for     
 *                                           coarse to fine block message      
 *                                           passing.                          
 *******************************************************************************/
int Euler2D_Quad_Block::LoadSendBuffer_C2F(double *buffer,
					   int &buffer_count,
					   const int buffer_size,
					   const int i_min, 
					   const int i_max,
					   const int i_inc,
					   const int j_min, 
					   const int j_max,
					   const int j_inc,
					   const int face,
					   const int sector) {
  int i, j, k;
  Vector2D dX;
  Euler2D_pState Wfine;
  Euler2D_cState Ufine;

  if (j_inc > 0) {
    if (i_inc > 0) {
      for ( j = j_min; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max); j += j_inc) {
	for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
	  // Perform limited linear least squares reconstruction in cell (i, j_min).
	  SubcellReconstruction(i, j, LIMITER_VENKATAKRISHNAN);
	  // Evaluate SW sub (fine) cell values if required.
	  if (!(face == NORTH && sector == WEST && Nghost%2 && j == j_min) &&
	      !(face == NORTH && sector == EAST && Nghost%2 && (i == i_min || j == j_min)) &&
	      !(face == SOUTH && sector == EAST && Nghost%2 && i == i_min) &&
	      !(face == EAST && sector == NORTH && Nghost%2 && (i == i_min || j == j_min)) &&
	      !(face == EAST && sector == SOUTH && Nghost%2 && i == i_min) &&
	      !(face == WEST && sector == NORTH && Nghost%2 && j == j_min) &&
	      !(face == NORTH_EAST && Nghost%2 && (i == i_min || j == j_min)) &&
	      !(face == NORTH_WEST && Nghost%2 && j == j_min) &&
	      !(face == SOUTH_EAST && Nghost%2 && i == i_min)) {
	    dX = Grid.centroidSW(i,j) - Grid.Cell[i][j].Xc;
	    Wfine = W[i][j] + (phi[i][j]^dWdx[i][j])*dX.x +
                              (phi[i][j]^dWdy[i][j])*dX.y;
	    Ufine = Wfine.U();
	    for ( k = 1 ; k <= NUM_VAR_EULER2D; ++ k) {
	      buffer_count = buffer_count + 1;
	      if (buffer_count >= buffer_size) return(1);
	      buffer[buffer_count] = Ufine[k];
	    } /* endfor */
	  } /* endif */
	  // Evaluate SE sub (fine) cell values if required.
	  if (!(face == NORTH && sector == WEST && Nghost%2 && (i == i_max || j == j_min)) &&
	      !(face == NORTH && sector == EAST && Nghost%2 && j == j_min) &&
	      !(face == SOUTH && sector == WEST && Nghost%2 && i == i_max) &&
	      !(face == EAST && sector == NORTH && Nghost%2 && j == j_min) &&
	      !(face == WEST && sector == NORTH && Nghost%2 && (i == i_max || j == j_min)) &&
	      !(face == WEST && sector == SOUTH && Nghost%2 && i == i_max) &&
	      !(face == NORTH_EAST && Nghost%2 && j == j_min) &&
	      !(face == NORTH_WEST && Nghost%2 && (i == i_max || j == j_min)) &&
	      !(face == SOUTH_WEST && Nghost%2 && i == i_max)) {
	    dX = Grid.centroidSE(i,j) - Grid.Cell[i][j].Xc;
	    Wfine = W[i][j] + (phi[i][j]^dWdx[i][j])*dX.x +
                              (phi[i][j]^dWdy[i][j])*dX.y;
	    Ufine = Wfine.U();
	    for ( k = 1 ; k <= NUM_VAR_EULER2D; ++ k) {
	      buffer_count = buffer_count + 1;
	      if (buffer_count >= buffer_size) return(1);
	      buffer[buffer_count] = Ufine[k];
	    } /* endfor */
	  } /* endif */
	} /* endfor */
	for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
	  // Evaluate NW sub (fine) cell values if required.
	  if (!(face == NORTH && sector == EAST && Nghost%2 && i == i_min) &&
	      !(face == SOUTH && sector == EAST && Nghost%2 && (i == i_min || j == j_max)) &&
	      !(face == SOUTH && sector == WEST && Nghost%2 && j == j_max) &&
	      !(face == EAST && sector == NORTH && Nghost%2 && i == i_min) &&
	      !(face == EAST && sector == SOUTH && Nghost%2 && (i == i_min || j == j_max)) &&
	      !(face == WEST && sector == SOUTH && Nghost%2 && j == j_max) &&
	      !(face == NORTH_EAST && Nghost%2 && i == i_min) &&
	      !(face == SOUTH_EAST && Nghost%2 && (i == i_min || j == j_max)) &&
	      !(face == SOUTH_WEST && Nghost%2 && j == j_max)) {
	    dX = Grid.centroidNW(i,j) - Grid.Cell[i][j].Xc;
	    Wfine = W[i][j] + (phi[i][j]^dWdx[i][j])*dX.x +
                              (phi[i][j]^dWdy[i][j])*dX.y;
	    Ufine = Wfine.U();
	    for ( k = 1 ; k <= NUM_VAR_EULER2D; ++ k) { 
	      buffer_count = buffer_count + 1;
	      if (buffer_count >= buffer_size) return(1);
	      buffer[buffer_count] = Ufine[k];
	    } /* endfor */
	  } /* endif */
	  // Evaluate NE sub (fine) cell values if required.
	  if (!(face == NORTH && sector == WEST && Nghost%2 && i == i_max) &&
	      !(face == SOUTH && sector == EAST && Nghost%2 && j == j_max) &&
	      !(face == SOUTH && sector == WEST && Nghost%2 && (i == i_max || j == j_max)) &&
	      !(face == EAST && sector == SOUTH && Nghost%2 && j == j_max) &&
	      !(face == WEST && sector == NORTH && Nghost%2 && i == i_max) &&
	      !(face == WEST && sector == SOUTH && Nghost%2 && (i == i_max || j == j_max)) &&
	      !(face == NORTH_WEST && Nghost%2 && i == i_max) &&
	      !(face == SOUTH_EAST && Nghost%2 && j == j_max) &&
	      !(face == SOUTH_WEST && Nghost%2 && (i == i_max || j == j_max))) {
	    dX = Grid.centroidNE(i,j) - Grid.Cell[i][j].Xc;
	    Wfine = W[i][j] + (phi[i][j]^dWdx[i][j])*dX.x +
                              (phi[i][j]^dWdy[i][j])*dX.y;
	    Ufine = Wfine.U();
	    for ( k = 1 ; k <= NUM_VAR_EULER2D; ++ k) {
	      buffer_count = buffer_count + 1;
	      if (buffer_count >= buffer_size) return(1);
	      buffer[buffer_count] = Ufine[k];
	    } /* endfor */
	  } /* endif */
	} /* endfor */
      } /* endfor */

      return 0;

    } /* endif */
  } /* endif */

  // Load send message buffer for the coarse-to-fine grid for cases in
  // which one (or both) of the increments is negative.  Only for two
  // ghost cells.

  if (CENO_Execution_Mode::USE_CENO_ALGORITHM &&
      CENO_Execution_Mode::HIGH_ORDER_MESSAGE_PASSING){ // High-order message passing
    throw runtime_error("Euler2D_Quad_Block::LoadSendBuffer_C2F() ERROR! The case in which one (or both) of the incements is negative has not been implemented for high-order message passing!");
  }

  if (j_min == j_max) { // North or south boundary.
     // Four different orderings to consider depending on the value of i_inc & j_inc.
     if (j_inc > 0) {
        if (i_inc > 0) {
 	   return 1;
        } else {
           for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
              // Perform limited linear least squares reconstruction in cell (i, j_min).
              SubcellReconstruction(i, j_min, LIMITER_VENKATAKRISHNAN);
              // Evaluate SE sub (fine) cell values.
              dX = (HALF*(Grid.Node[i][j_min].X+Grid.Node[i+1][j_min].X)+
                    Grid.Node[i+1][j_min].X+
                    Grid.Cell[i][j_min].Xc+
                    HALF*(Grid.Node[i+1][j_min].X+Grid.Node[i+1][j_min+1].X))/FOUR -
                   Grid.Cell[i][j_min].Xc;
              Wfine = W[i][j_min] +
                      (phi[i][j_min]^dWdx[i][j_min])*dX.x +
                      (phi[i][j_min]^dWdy[i][j_min])*dX.y;
              Ufine = Wfine.U();
              for ( k = 1 ; k <= NUM_VAR_EULER2D; ++ k) {
  	         buffer_count = buffer_count + 1;
                 if (buffer_count >= buffer_size) return(1);
                 buffer[buffer_count] = Ufine[k];
              } /* endfor */
              // Evaluate SW sub (fine) cell values.
              dX = (Grid.Node[i][j_min].X+
                    HALF*(Grid.Node[i][j_min].X+Grid.Node[i+1][j_min].X)+
                    HALF*(Grid.Node[i][j_min].X+Grid.Node[i][j_min+1].X)+
                    Grid.Cell[i][j_min].Xc)/FOUR -
                   Grid.Cell[i][j_min].Xc;
              Wfine = W[i][j_min] +
                      (phi[i][j_min]^dWdx[i][j_min])*dX.x +
                      (phi[i][j_min]^dWdy[i][j_min])*dX.y;
              Ufine = Wfine.U();
              for ( k = 1 ; k <= NUM_VAR_EULER2D; ++ k) {
  	         buffer_count = buffer_count + 1;
                 if (buffer_count >= buffer_size) return(1);
                 buffer[buffer_count] = Ufine[k];
              } /* endfor */
           } /* endfor */
           for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
              // Evaluate NE sub (fine) cell values.
              dX = (Grid.Cell[i][j_min].Xc+
                    HALF*(Grid.Node[i+1][j_min].X+Grid.Node[i+1][j_min+1].X)+
                    HALF*(Grid.Node[i][j_min+1].X+Grid.Node[i+1][j_min+1].X)+
                    Grid.Node[i+1][j_min+1].X)/FOUR -
                   Grid.Cell[i][j_min].Xc;
              Wfine = W[i][j_min] +
                      (phi[i][j_min]^dWdx[i][j_min])*dX.x +
                      (phi[i][j_min]^dWdy[i][j_min])*dX.y;
              Ufine = Wfine.U();
              for ( k = 1 ; k <= NUM_VAR_EULER2D; ++ k) {
  	         buffer_count = buffer_count + 1;
                 if (buffer_count >= buffer_size) return(1);
                 buffer[buffer_count] = Ufine[k];
              } /* endfor */
              // Evaluate NW sub (fine) cell values.
              dX = (HALF*(Grid.Node[i][j_min].X+Grid.Node[i][j_min+1].X)+
                    Grid.Cell[i][j_min].Xc+
                    Grid.Node[i][j_min+1].X+
                    HALF*(Grid.Node[i][j_min+1].X+Grid.Node[i+1][j_min+1].X))/FOUR -
                   Grid.Cell[i][j_min].Xc;
              Wfine = W[i][j_min] +
                      (phi[i][j_min]^dWdx[i][j_min])*dX.x +
                      (phi[i][j_min]^dWdy[i][j_min])*dX.y;
              Ufine = Wfine.U();
              for ( k = 1 ; k <= NUM_VAR_EULER2D; ++ k) { 
                 buffer_count = buffer_count + 1;
                 if (buffer_count >= buffer_size) return(1);
                 buffer[buffer_count] = Ufine[k];
              } /* endfor */
           } /* endfor */
        } /* endif */
     } else {
        if (i_inc > 0) {
           for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
              // Perform limited linear least squares reconstruction in cell (i, j_min).
              SubcellReconstruction(i, j_min, LIMITER_VENKATAKRISHNAN);
              // Evaluate NW sub (fine) cell values.
              dX = (HALF*(Grid.Node[i][j_min].X+Grid.Node[i][j_min+1].X)+
                    Grid.Cell[i][j_min].Xc+
                    Grid.Node[i][j_min+1].X+
                    HALF*(Grid.Node[i][j_min+1].X+Grid.Node[i+1][j_min+1].X))/FOUR -
                   Grid.Cell[i][j_min].Xc;
              Wfine = W[i][j_min] +
                      (phi[i][j_min]^dWdx[i][j_min])*dX.x +
                      (phi[i][j_min]^dWdy[i][j_min])*dX.y;
              Ufine = Wfine.U();
              for ( k = 1 ; k <= NUM_VAR_EULER2D; ++ k) { 
                 buffer_count = buffer_count + 1;
                 if (buffer_count >= buffer_size) return(1);
                 buffer[buffer_count] = Ufine[k];
              } /* endfor */
              // Evaluate NE sub (fine) cell values.
              dX = (Grid.Cell[i][j_min].Xc+
                    HALF*(Grid.Node[i+1][j_min].X+Grid.Node[i+1][j_min+1].X)+
                    HALF*(Grid.Node[i][j_min+1].X+Grid.Node[i+1][j_min+1].X)+
                    Grid.Node[i+1][j_min+1].X)/FOUR -
                   Grid.Cell[i][j_min].Xc;
              Wfine = W[i][j_min] +
                      (phi[i][j_min]^dWdx[i][j_min])*dX.x +
                      (phi[i][j_min]^dWdy[i][j_min])*dX.y;
              Ufine = Wfine.U();
              for ( k = 1 ; k <= NUM_VAR_EULER2D; ++ k) {
  	         buffer_count = buffer_count + 1;
                 if (buffer_count >= buffer_size) return(1);
                 buffer[buffer_count] = Ufine[k];
              } /* endfor */
           } /* endfor */
           for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
              // Evaluate SW sub (fine) cell values.
              dX = (Grid.Node[i][j_min].X+
                    HALF*(Grid.Node[i][j_min].X+Grid.Node[i+1][j_min].X)+
                    HALF*(Grid.Node[i][j_min].X+Grid.Node[i][j_min+1].X)+
                    Grid.Cell[i][j_min].Xc)/FOUR -
                   Grid.Cell[i][j_min].Xc;
              Wfine = W[i][j_min] +
                      (phi[i][j_min]^dWdx[i][j_min])*dX.x +
                      (phi[i][j_min]^dWdy[i][j_min])*dX.y;
              Ufine = Wfine.U();
              for ( k = 1 ; k <= NUM_VAR_EULER2D; ++ k) {
  	         buffer_count = buffer_count + 1;
                 if (buffer_count >= buffer_size) return(1);
                 buffer[buffer_count] = Ufine[k];
              } /* endfor */
              // Evaluate SE sub (fine) cell values.
              dX = (HALF*(Grid.Node[i][j_min].X+Grid.Node[i+1][j_min].X)+
                    Grid.Node[i+1][j_min].X+
                    Grid.Cell[i][j_min].Xc+
                    HALF*(Grid.Node[i+1][j_min].X+Grid.Node[i+1][j_min+1].X))/FOUR -
                   Grid.Cell[i][j_min].Xc;
              Wfine = W[i][j_min] +
                      (phi[i][j_min]^dWdx[i][j_min])*dX.x +
                      (phi[i][j_min]^dWdy[i][j_min])*dX.y;
              Ufine = Wfine.U();
              for ( k = 1 ; k <= NUM_VAR_EULER2D; ++ k) {
  	         buffer_count = buffer_count + 1;
                 if (buffer_count >= buffer_size) return(1);
                 buffer[buffer_count] = Ufine[k];
              } /* endfor */
           } /* endfor */
        } else {
           for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
              // Perform limited linear least squares reconstruction in cell (i, j_min).
              SubcellReconstruction(i, j_min, LIMITER_VENKATAKRISHNAN);
              // Evaluate NE sub (fine) cell values.
              dX = (Grid.Cell[i][j_min].Xc+
                    HALF*(Grid.Node[i+1][j_min].X+Grid.Node[i+1][j_min+1].X)+
                    HALF*(Grid.Node[i][j_min+1].X+Grid.Node[i+1][j_min+1].X)+
                    Grid.Node[i+1][j_min+1].X)/FOUR -
                   Grid.Cell[i][j_min].Xc;
              Wfine = W[i][j_min] +
                      (phi[i][j_min]^dWdx[i][j_min])*dX.x +
                      (phi[i][j_min]^dWdy[i][j_min])*dX.y;
              Ufine = Wfine.U();
              for ( k = 1 ; k <= NUM_VAR_EULER2D; ++ k) {
  	         buffer_count = buffer_count + 1;
                 if (buffer_count >= buffer_size) return(1);
                 buffer[buffer_count] = Ufine[k];
              } /* endfor */
              // Evaluate NW sub (fine) cell values.
              dX = (HALF*(Grid.Node[i][j_min].X+Grid.Node[i][j_min+1].X)+
                    Grid.Cell[i][j_min].Xc+
                    Grid.Node[i][j_min+1].X+
                    HALF*(Grid.Node[i][j_min+1].X+Grid.Node[i+1][j_min+1].X))/FOUR -
                   Grid.Cell[i][j_min].Xc;
              Wfine = W[i][j_min] +
                      (phi[i][j_min]^dWdx[i][j_min])*dX.x +
                      (phi[i][j_min]^dWdy[i][j_min])*dX.y;
              Ufine = Wfine.U();
              for ( k = 1 ; k <= NUM_VAR_EULER2D; ++ k) { 
                 buffer_count = buffer_count + 1;
                 if (buffer_count >= buffer_size) return(1);
                 buffer[buffer_count] = Ufine[k];
              } /* endfor */
           } /* endfor */
           for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
              // Evaluate SE sub (fine) cell values.
              dX = (HALF*(Grid.Node[i][j_min].X+Grid.Node[i+1][j_min].X)+
                    Grid.Node[i+1][j_min].X+
                    Grid.Cell[i][j_min].Xc+
                    HALF*(Grid.Node[i+1][j_min].X+Grid.Node[i+1][j_min+1].X))/FOUR -
                   Grid.Cell[i][j_min].Xc;
              Wfine = W[i][j_min] +
                      (phi[i][j_min]^dWdx[i][j_min])*dX.x +
                      (phi[i][j_min]^dWdy[i][j_min])*dX.y;
              Ufine = Wfine.U();
              for ( k = 1 ; k <= NUM_VAR_EULER2D; ++ k) {
  	         buffer_count = buffer_count + 1;
                 if (buffer_count >= buffer_size) return(1);
                 buffer[buffer_count] = Ufine[k];
              } /* endfor */
              // Evaluate SW sub (fine) cell values.
              dX = (Grid.Node[i][j_min].X+
                    HALF*(Grid.Node[i][j_min].X+Grid.Node[i+1][j_min].X)+
                    HALF*(Grid.Node[i][j_min].X+Grid.Node[i][j_min+1].X)+
                    Grid.Cell[i][j_min].Xc)/FOUR -
                   Grid.Cell[i][j_min].Xc;
              Wfine = W[i][j_min] +
                      (phi[i][j_min]^dWdx[i][j_min])*dX.x +
                      (phi[i][j_min]^dWdy[i][j_min])*dX.y;
              Ufine = Wfine.U();
              for ( k = 1 ; k <= NUM_VAR_EULER2D; ++ k) {
  	         buffer_count = buffer_count + 1;
                 if (buffer_count >= buffer_size) return(1);
                 buffer[buffer_count] = Ufine[k];
              } /* endfor */
           } /* endfor */
        } /* endif */
     } /* endif */
  } else { // East or west boundary.
     // Four different orderings to consider depending on the value of i_inc & j_inc.
     if (j_inc > 0) {
        if (i_inc > 0) {
 	   return 1;
        } else {
           for ( j = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc ) {
              // Perform limited linear least squares reconstruction in cell (i_min, j).
              SubcellReconstruction(i_min, j, LIMITER_VENKATAKRISHNAN);
              // Evaluate SE sub (fine) cell values.
              dX = (HALF*(Grid.Node[i_min][j].X+Grid.Node[i_min+1][j].X)+
                    Grid.Node[i_min+1][j].X+
                    Grid.Cell[i_min][j].Xc+
                    HALF*(Grid.Node[i_min+1][j].X+Grid.Node[i_min+1][j+1].X))/FOUR -
                   Grid.Cell[i_min][j].Xc;
              Wfine = W[i_min][j] +
                      (phi[i_min][j]^dWdx[i_min][j])*dX.x +
                      (phi[i_min][j]^dWdy[i_min][j])*dX.y;
              Ufine = Wfine.U();
              for ( k = 1 ; k <= NUM_VAR_EULER2D; ++ k) {
  	         buffer_count = buffer_count + 1;
                 if (buffer_count >= buffer_size) return(1);
                 buffer[buffer_count] = Ufine[k];
              } /* endfor */
              // Evaluate SW sub (fine) cell values.
              dX = (Grid.Node[i_min][j].X+
                    HALF*(Grid.Node[i_min][j].X+Grid.Node[i_min+1][j].X)+
                    HALF*(Grid.Node[i_min][j].X+Grid.Node[i_min][j+1].X)+
                    Grid.Cell[i_min][j].Xc)/FOUR -
                   Grid.Cell[i_min][j].Xc;
              Wfine = W[i_min][j] +
                      (phi[i_min][j]^dWdx[i_min][j])*dX.x +
                      (phi[i_min][j]^dWdy[i_min][j])*dX.y;
              Ufine = Wfine.U();
              for ( k = 1 ; k <= NUM_VAR_EULER2D; ++ k) {
  	         buffer_count = buffer_count + 1;
                 if (buffer_count >= buffer_size) return(1);
                 buffer[buffer_count] = Ufine[k];
              } /* endfor */
              // Evaluate NE sub (fine) cell values.
             dX = (Grid.Cell[i_min][j].Xc+
                    HALF*(Grid.Node[i_min+1][j].X+Grid.Node[i_min+1][j+1].X)+
                    HALF*(Grid.Node[i_min][j+1].X+Grid.Node[i_min+1][j+1].X)+
                    Grid.Node[i_min+1][j+1].X)/FOUR -
                   Grid.Cell[i_min][j].Xc;
              Wfine = W[i_min][j] +
                      (phi[i_min][j]^dWdx[i_min][j])*dX.x +
                      (phi[i_min][j]^dWdy[i_min][j])*dX.y;
              Ufine = Wfine.U();
              for ( k = 1 ; k <= NUM_VAR_EULER2D; ++ k) {
                 buffer_count = buffer_count + 1;
                 if (buffer_count >= buffer_size) return(1);
                 buffer[buffer_count] = Ufine[k];
              } /* endfor */
              // Evaluate NW sub (fine) cell values.
              dX = (HALF*(Grid.Node[i_min][j].X+Grid.Node[i_min][j+1].X)+
                    Grid.Cell[i_min][j].Xc+
                    Grid.Node[i_min][j+1].X+
                    HALF*(Grid.Node[i_min][j+1].X+Grid.Node[i_min+1][j+1].X))/FOUR -
                   Grid.Cell[i_min][j].Xc;
              Wfine = W[i_min][j] +
                      (phi[i_min][j]^dWdx[i_min][j])*dX.x +
                      (phi[i_min][j]^dWdy[i_min][j])*dX.y;
              Ufine = Wfine.U();
              for ( k = 1 ; k <= NUM_VAR_EULER2D; ++ k) { 
                 buffer_count = buffer_count + 1;
                 if (buffer_count >= buffer_size) return(1);
                 buffer[buffer_count] = Ufine[k];
              } /* endfor */
           } /* endfor */
        } /* endif */
     } else {
        if (i_inc > 0) {
           for ( j = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc ) {
              // Perform limited linear least squares reconstruction in cell (i_min, j).
              SubcellReconstruction(i_min, j, LIMITER_VENKATAKRISHNAN);
              // Evaluate NW sub (fine) cell values.
              dX = (HALF*(Grid.Node[i_min][j].X+Grid.Node[i_min][j+1].X)+
                    Grid.Cell[i_min][j].Xc+
                    Grid.Node[i_min][j+1].X+
                    HALF*(Grid.Node[i_min][j+1].X+Grid.Node[i_min+1][j+1].X))/FOUR -
                   Grid.Cell[i_min][j].Xc;
              Wfine = W[i_min][j] +
                      (phi[i_min][j]^dWdx[i_min][j])*dX.x +
                      (phi[i_min][j]^dWdy[i_min][j])*dX.y;
              Ufine = Wfine.U();
              for ( k = 1 ; k <= NUM_VAR_EULER2D; ++ k) { 
                 buffer_count = buffer_count + 1;
                 if (buffer_count >= buffer_size) return(1);
                 buffer[buffer_count] = Ufine[k];
              } /* endfor */
              // Evaluate NE sub (fine) cell values.
              dX = (Grid.Cell[i_min][j].Xc+
                    HALF*(Grid.Node[i_min+1][j].X+Grid.Node[i_min+1][j+1].X)+
                    HALF*(Grid.Node[i_min][j+1].X+Grid.Node[i_min+1][j+1].X)+
                    Grid.Node[i_min+1][j+1].X)/FOUR -
                   Grid.Cell[i_min][j].Xc;
              Wfine = W[i_min][j] +
                      (phi[i_min][j]^dWdx[i_min][j])*dX.x +
                      (phi[i_min][j]^dWdy[i_min][j])*dX.y;
              Ufine = Wfine.U();
              for ( k = 1 ; k <= NUM_VAR_EULER2D; ++ k) {
                 buffer_count = buffer_count + 1;
                 if (buffer_count >= buffer_size) return(1);
                 buffer[buffer_count] = Ufine[k];
              } /* endfor */
              // Evaluate SW sub (fine) cell values.
              dX = (Grid.Node[i_min][j].X+
                    HALF*(Grid.Node[i_min][j].X+Grid.Node[i_min+1][j].X)+
                    HALF*(Grid.Node[i_min][j].X+Grid.Node[i_min][j+1].X)+
                    Grid.Cell[i_min][j].Xc)/FOUR -
                   Grid.Cell[i_min][j].Xc;
              Wfine = W[i_min][j] +
                      (phi[i_min][j]^dWdx[i_min][j])*dX.x +
                      (phi[i_min][j]^dWdy[i_min][j])*dX.y;
              Ufine = Wfine.U();
              for ( k = 1 ; k <= NUM_VAR_EULER2D; ++ k) {
  	         buffer_count = buffer_count + 1;
                 if (buffer_count >= buffer_size) return(1);
                 buffer[buffer_count] = Ufine[k];
              } /* endfor */
              // Evaluate SE sub (fine) cell values.
              dX = (HALF*(Grid.Node[i_min][j].X+Grid.Node[i_min+1][j].X)+
                    Grid.Node[i_min+1][j].X+
                    Grid.Cell[i_min][j].Xc+
                    HALF*(Grid.Node[i_min+1][j].X+Grid.Node[i_min+1][j+1].X))/FOUR -
                   Grid.Cell[i_min][j].Xc;
              Wfine = W[i_min][j] +
                      (phi[i_min][j]^dWdx[i_min][j])*dX.x +
                      (phi[i_min][j]^dWdy[i_min][j])*dX.y;
              Ufine = Wfine.U();
              for ( k = 1 ; k <= NUM_VAR_EULER2D; ++ k) {
  	         buffer_count = buffer_count + 1;
                 if (buffer_count >= buffer_size) return(1);
                 buffer[buffer_count] = Ufine[k];
              } /* endfor */
           } /* endfor */
        } else {
           for ( j = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc ) {
              // Perform limited linear least squares reconstruction in cell (i_min, j).
              SubcellReconstruction(i_min, j, LIMITER_VENKATAKRISHNAN);
              // Evaluate NE sub (fine) cell values.
              dX = (Grid.Cell[i_min][j].Xc+
                    HALF*(Grid.Node[i_min+1][j].X+Grid.Node[i_min+1][j+1].X)+
                    HALF*(Grid.Node[i_min][j+1].X+Grid.Node[i_min+1][j+1].X)+
                    Grid.Node[i_min+1][j+1].X)/FOUR -
                   Grid.Cell[i_min][j].Xc;
              Wfine = W[i_min][j] +
                      (phi[i_min][j]^dWdx[i_min][j])*dX.x +
                      (phi[i_min][j]^dWdy[i_min][j])*dX.y;
              Ufine = Wfine.U();
              for ( k = 1 ; k <= NUM_VAR_EULER2D; ++ k) {
                 buffer_count = buffer_count + 1;
                 if (buffer_count >= buffer_size) return(1);
                 buffer[buffer_count] = Ufine[k];
              } /* endfor */
              // Evaluate NW sub (fine) cell values.
              dX = (HALF*(Grid.Node[i_min][j].X+Grid.Node[i_min][j+1].X)+
                    Grid.Cell[i_min][j].Xc+
                    Grid.Node[i_min][j+1].X+
                    HALF*(Grid.Node[i_min][j+1].X+Grid.Node[i_min+1][j+1].X))/FOUR -
                   Grid.Cell[i_min][j].Xc;
              Wfine = W[i_min][j] +
                      (phi[i_min][j]^dWdx[i_min][j])*dX.x +
                      (phi[i_min][j]^dWdy[i_min][j])*dX.y;
              Ufine = Wfine.U();
              for ( k = 1 ; k <= NUM_VAR_EULER2D; ++ k) { 
                 buffer_count = buffer_count + 1;
                 if (buffer_count >= buffer_size) return(1);
                 buffer[buffer_count] = Ufine[k];
              } /* endfor */
              // Evaluate SE sub (fine) cell values.
              dX = (HALF*(Grid.Node[i_min][j].X+Grid.Node[i_min+1][j].X)+
                    Grid.Node[i_min+1][j].X+
                    Grid.Cell[i_min][j].Xc+
                    HALF*(Grid.Node[i_min+1][j].X+Grid.Node[i_min+1][j+1].X))/FOUR -
                   Grid.Cell[i_min][j].Xc;
              Wfine = W[i_min][j] +
                      (phi[i_min][j]^dWdx[i_min][j])*dX.x +
                      (phi[i_min][j]^dWdy[i_min][j])*dX.y;
              Ufine = Wfine.U();
              for ( k = 1 ; k <= NUM_VAR_EULER2D; ++ k) {
  	         buffer_count = buffer_count + 1;
                 if (buffer_count >= buffer_size) return(1);
                 buffer[buffer_count] = Ufine[k];
              } /* endfor */
              // Evaluate SW sub (fine) cell values.
              dX = (Grid.Node[i_min][j].X+
                    HALF*(Grid.Node[i_min][j].X+Grid.Node[i_min+1][j].X)+
                    HALF*(Grid.Node[i_min][j].X+Grid.Node[i_min][j+1].X)+
                    Grid.Cell[i_min][j].Xc)/FOUR -
                   Grid.Cell[i_min][j].Xc;
              Wfine = W[i_min][j] +
                      (phi[i_min][j]^dWdx[i_min][j])*dX.x +
                      (phi[i_min][j]^dWdy[i_min][j])*dX.y;
              Ufine = Wfine.U();
              for ( k = 1 ; k <= NUM_VAR_EULER2D; ++ k) {
  	         buffer_count = buffer_count + 1;
                 if (buffer_count >= buffer_size) return(1);
                 buffer[buffer_count] = Ufine[k];
              } /* endfor */
           } /* endfor */
        } /* endif */
     } /* endif */
  } /* endif */
  return(0);
}

/*********************************************************//**
 * This routine evaluates the residual for the specified
 * solution block using a 2nd-order limited upwind      
 * finite-volume spatial discretization scheme with     
 * either the Godunov, Roe, Rusanov, HLLE, Linde, or    
 * HLLC flux functions.                                 
 * The residual is stored in dUdt[][][0].               
 *                                                      
 *************************************************************/
int Euler2D_Quad_Block::dUdt_Residual_Evaluation(const Euler2D_Input_Parameters &IP){

  int i, j;
  Vector2D dX;
  Euler2D_pState Wl, Wr, Wi;
  Euler2D_cState Flux;

  /* Perform the linear reconstruction within each cell
     of the computational grid for this stage. */
    
  switch(IP.i_Reconstruction) {
  case RECONSTRUCTION_GREEN_GAUSS :
    Linear_Reconstruction_GreenGauss(*this,
				     IP.i_Limiter);    
    break;
  case RECONSTRUCTION_LEAST_SQUARES :
    Linear_Reconstruction_LeastSquares(*this,
				       IP.i_Limiter);
    break;
  default:
    Linear_Reconstruction_LeastSquares(*this,
				       IP.i_Limiter);
    break;
  } /* endswitch */

    /* Evaluate the time rate of change of the solution
       (i.e., the solution residuals) using a second-order
       limited upwind scheme with a variety of flux functions. */
    
    // Add i-direction (zeta-direction) fluxes.
  for ( j  = JCl-1 ; j <= JCu+1 ; ++j ) {
    dUdt[ICl-1][j][0] = Euler2D_U_VACUUM;
          
    for ( i = ICl-1 ; i <= ICu ; ++i ) {

      dUdt[i+1][j][0] = Euler2D_U_VACUUM;
    
      if ( j > JCl-1 && j < JCu+1 ) {
    
	/* Evaluate the cell interface i-direction fluxes. */
	  
	if (i == ICl-1 && 
	    (Grid.BCtypeW[j] == BC_REFLECTION ||
	     Grid.BCtypeW[j] == BC_CHARACTERISTIC ||
	     Grid.BCtypeW[j] == BC_CHARACTERISTIC_VELOCITY ||
	     Grid.BCtypeW[j] == BC_BURNING_SURFACE ||
	     Grid.BCtypeW[j] == BC_MASS_INJECTION ||
	     Grid.BCtypeW[j] == BC_RINGLEB_FLOW ||
	     Grid.BCtypeW[j] == BC_FROZEN ||
	     Grid.BCtypeW[j] == BC_EXACT_SOLUTION ||
	     Grid.BCtypeW[j] == BC_WALL_INVISCID)) {
	  dX = Grid.xfaceW(i+1, j)-Grid.Cell[i+1][j].Xc;
	  Wr = W[i+1][j] + 
	    (phi[i+1][j]^dWdx[i+1][j])*dX.x +
	    (phi[i+1][j]^dWdy[i+1][j])*dX.y;
	  if (Grid.BCtypeW[j] == BC_REFLECTION) {
	    Wl = Reflect(Wr, Grid.nfaceW(i+1, j));
	  } else if (Grid.BCtypeW[j] == BC_BURNING_SURFACE) {
	    Wl = BurningSurface(Wr, Grid.nfaceW(i+1, j));
	  } else if (Grid.BCtypeW[j] == BC_MASS_INJECTION) {
	    Wl = MassInjection(Wr,Grid.nfaceW(i+1,j),OFF);
	  } else if (Grid.BCtypeW[j] == BC_CHARACTERISTIC) {
	    Wl = BC_Characteristic_Pressure(Wr, 
					    WoW[j], 
					    Grid.nfaceW(i+1, j));
	  } else if (Grid.BCtypeW[j] == BC_CHARACTERISTIC_VELOCITY) {
	    Wl = BC_Characteristic(Wr, 
				   WoW[j], 
				   Grid.nfaceW(i+1, j));
	  } else if (Grid.BCtypeW[j] == BC_RINGLEB_FLOW) {
	    Wl = RinglebFlow(Wl,Grid.xfaceW(i+1,j));
	    //Wl = Reflect(Wr, Grid.nfaceW(i+1, j));
	  } else if (Grid.BCtypeW[j] == BC_FROZEN) {
	    // calculate Wl based on the ghost cell reconstruction
	    Wl = PiecewiseLinearSolutionAtLocation(i,j,
						   Grid.xfaceW(i+1,j));
	  } else if (Grid.BCtypeW[j] == BC_EXACT_SOLUTION) {
	    // calculate Wl using the exact solution at the interface and avoiding the problem to become ill-posed.
	    if (ExactSoln->IsExactSolutionSet()){
	      Wl = BC_Characteristic_Pressure(Wr,
					      ExactSoln->Solution(Grid.xfaceW(i+1,j).x,Grid.xfaceW(i+1,j).y),
					      Grid.nfaceW(i+1, j));
	    } else {
	      throw runtime_error("Euler2D_Quad_Block::dUdt_Residual_Evaluation() ERROR! There is no exact solution set for the Exact_Solution BC.");
	    }
	  } else if (Grid.BCtypeW[j] == BC_WALL_INVISCID) {
	    Wl = Reflect(Wr, Grid.nfaceW(i+1, j));
	    Wi = PiecewiseLinearSolutionAtLocation(i,j,
						  Grid.xfaceW(i+1,j));
	    Wl.d = Wi.d;
	    Wl.p = Wi.p;
	  } /* endif */
	} else if (i == ICu && 
		   (Grid.BCtypeE[j] == BC_REFLECTION ||
		    Grid.BCtypeE[j] == BC_CHARACTERISTIC ||
		    Grid.BCtypeE[j] == BC_CHARACTERISTIC_VELOCITY ||
		    Grid.BCtypeE[j] == BC_BURNING_SURFACE ||
		    Grid.BCtypeE[j] == BC_MASS_INJECTION ||
		    Grid.BCtypeE[j] == BC_RINGLEB_FLOW ||
		    Grid.BCtypeE[j] == BC_FROZEN ||
		    Grid.BCtypeE[j] == BC_EXACT_SOLUTION ||
		    Grid.BCtypeE[j] == BC_WALL_INVISCID)) {
	  dX = Grid.xfaceE(i, j)-Grid.Cell[i][j].Xc;
	  Wl = W[i][j] + 
	    (phi[i][j]^dWdx[i][j])*dX.x +
	    (phi[i][j]^dWdy[i][j])*dX.y;
	  if (Grid.BCtypeE[j] == BC_REFLECTION) {
	    Wr = Reflect(Wl, Grid.nfaceE(i, j));
	  } else if (Grid.BCtypeE[j] == BC_BURNING_SURFACE) {
	    Wr = BurningSurface(Wl, Grid.nfaceE(i, j));
	  } else if (Grid.BCtypeE[j] == BC_MASS_INJECTION) {
	    Wr = MassInjection(Wl,Grid.nfaceE(i,j),OFF);
	  } else if (Grid.BCtypeE[j] == BC_CHARACTERISTIC) {
	    Wr = BC_Characteristic_Pressure(Wl, 
					    WoE[j], 
					    Grid.nfaceE(i, j));
	  } else if (Grid.BCtypeE[j] == BC_CHARACTERISTIC_VELOCITY) {
	    Wr = BC_Characteristic(Wl, 
				   WoE[j], 
				   Grid.nfaceE(i, j));
	  } else if (Grid.BCtypeE[j] == BC_RINGLEB_FLOW) {
	    Wr = RinglebFlow(Wr,Grid.xfaceE(i,j));
	    //Wr = Reflect(Wl, Grid.nfaceE(i, j));
	  } else if (Grid.BCtypeE[j] == BC_FROZEN) {
	    // calculate Wr based on the ghost cell reconstruction
	    Wr = PiecewiseLinearSolutionAtLocation(i+1,j,
						   Grid.xfaceE(i,j));
	  } else if (Grid.BCtypeE[j] == BC_EXACT_SOLUTION) {
	    // calculate Wr using the exact solution at the interface and avoiding the problem to become ill-posed.
	    if (ExactSoln->IsExactSolutionSet()){
	      Wr = BC_Characteristic_Pressure(Wl,
					      ExactSoln->Solution(Grid.xfaceE(i,j).x,Grid.xfaceE(i,j).y),
					      Grid.nfaceE(i, j));
	    } else {
	      throw runtime_error("Euler2D_Quad_Block::dUdt_Residual_Evaluation() ERROR! There is no exact solution set for the Exact_Solution BC.");
	    }	    
	  } else if (Grid.BCtypeE[j] == BC_WALL_INVISCID) {
	    Wr = Reflect(Wl, Grid.nfaceE(i, j));
	    Wi = PiecewiseLinearSolutionAtLocation(i+1,j,
						   Grid.xfaceW(i+1, j));
	    Wl.d = Wi.d;
	    Wl.p = Wi.p;
	  } /* endif */
	} else {            
	  dX = Grid.xfaceE(i, j)-Grid.Cell[i][j].Xc;
	  Wl = W[i][j] + 
	    (phi[i][j]^dWdx[i][j])*dX.x +
	    (phi[i][j]^dWdy[i][j])*dX.y;
	  dX = Grid.xfaceW(i+1, j)-Grid.Cell[i+1][j].Xc;
	  Wr = W[i+1][j] + 
	    (phi[i+1][j]^dWdx[i+1][j])*dX.x +
	    (phi[i+1][j]^dWdy[i+1][j])*dX.y;
	} /* endif */

	switch(IP.i_Flux_Function) {
	case FLUX_FUNCTION_GODUNOV :
	  Flux = FluxGodunov_n(Wl, Wr, Grid.nfaceE(i, j));
	  break;
	case FLUX_FUNCTION_ROE :
	  Flux = FluxRoe_n(Wl, Wr, Grid.nfaceE(i, j));
	  break;
	case FLUX_FUNCTION_RUSANOV :
	  Flux = FluxRusanov_n(Wl, Wr, Grid.nfaceE(i, j));
	  break;
	case FLUX_FUNCTION_HLLE :
	  Flux = FluxHLLE_n(Wl, Wr, Grid.nfaceE(i, j));
	  break;
	case FLUX_FUNCTION_LINDE :
	  Flux = FluxLinde_n(Wl, Wr, Grid.nfaceE(i, j));
	  break;
	case FLUX_FUNCTION_HLLC :
	  Flux = FluxHLLC_n(Wl, Wr, Grid.nfaceE(i, j));
	  break;
	case FLUX_FUNCTION_VANLEER :
	  Flux = FluxVanLeer_n(Wl, Wr, Grid.nfaceE(i,j));
	  break;
	case FLUX_FUNCTION_AUSM :
	  Flux = FluxAUSM_n(Wl, Wr, Grid.nfaceE(i,j));
	  break;
	case FLUX_FUNCTION_AUSMplus :
	  Flux = FluxAUSMplus_n(Wl, Wr, Grid.nfaceE(i,j));
	  break;
	case FLUX_FUNCTION_ROE_PRECON_WS :
	  Flux = FluxRoe_n_Precon_WS(Wl, Wr, Grid.nfaceE(i, j));
	  break;
	case FLUX_FUNCTION_HLLE_PRECON_WS :
	  Flux = FluxHLLE_n_Precon_WS(Wl, Wr, Grid.nfaceE(i, j));
	  break;
	default:
	  Flux = FluxRoe_n(Wl, Wr, Grid.nfaceE(i, j));
	  break;
	} /* endswitch */

	  /* Evaluate cell-averaged solution changes. */
	  
	dUdt[i][j][0] -= 
	  Flux*Grid.lfaceE(i, j)/
	  Grid.Cell[i][j].A;
	dUdt[i+1][j][0] += 
	  Flux*Grid.lfaceW(i+1, j)/
	  Grid.Cell[i+1][j].A;

	/* Include axisymmetric source terms as required. */

	if (Axisymmetric) {
	  dUdt[i][j][0] += 
	    S(W[i][j], Grid.Cell[i][j].Xc);
	} /* endif */

	  /* Save west and east face boundary flux. */
	  
	if (i == ICl-1) {
	  FluxW[j] = -Flux*Grid.lfaceW(i+1, j);
	} else if (i == ICu) {
	  FluxE[j] = Flux*Grid.lfaceE(i, j);
	} /* endif */ 

      } /* endif */
    } /* endfor */
      
    if ( j > JCl-1 && j < JCu+1 ) {
      dUdt[ICl-1][j][0] = Euler2D_U_VACUUM;
      dUdt[ICu+1][j][0] = Euler2D_U_VACUUM;
    } /* endif */
  } /* endfor */
    
    // Add j-direction (eta-direction) fluxes.
  for ( i = ICl ; i <= ICu ; ++i ) {
    for ( j  = JCl-1 ; j <= JCu ; ++j ) {
	
      /* Evaluate the cell interface j-direction fluxes. */
	
      if (j == JCl-1 && 
	  (Grid.BCtypeS[i] == BC_REFLECTION ||
	   Grid.BCtypeS[i] == BC_CHARACTERISTIC ||
	   Grid.BCtypeS[i] == BC_CHARACTERISTIC_VELOCITY ||
	   Grid.BCtypeS[i] == BC_BURNING_SURFACE ||
	   Grid.BCtypeS[i] == BC_MASS_INJECTION ||
	   Grid.BCtypeS[i] == BC_RINGLEB_FLOW ||
	   Grid.BCtypeS[i] == BC_FROZEN ||
	   Grid.BCtypeS[i] == BC_EXACT_SOLUTION ||
	   Grid.BCtypeS[i] == BC_WALL_INVISCID)) {
	dX = Grid.xfaceS(i, j+1)-Grid.Cell[i][j+1].Xc;
	Wr = W[i][j+1] +
	  (phi[i][j+1]^dWdx[i][j+1])*dX.x +
	  (phi[i][j+1]^dWdy[i][j+1])*dX.y;
	if (Grid.BCtypeS[i] == BC_REFLECTION) {
	  Wl = Reflect(Wr, Grid.nfaceS(i, j+1));
	} else if (Grid.BCtypeS[i] == BC_BURNING_SURFACE) {
	  Wl = BurningSurface(Wr, Grid.nfaceS(i, j+1));
	} else if (Grid.BCtypeS[i] == BC_MASS_INJECTION) {
	  Wl = MassInjection(Wr,Grid.nfaceS(i,j+1),OFF);
	} else if (Grid.BCtypeS[i] == BC_CHARACTERISTIC) {
	  Wl = BC_Characteristic_Pressure(Wr, 
					  WoS[i], 
					  Grid.nfaceS(i, j+1));
	} else if (Grid.BCtypeS[i] == BC_CHARACTERISTIC_VELOCITY) {
	  Wl = BC_Characteristic(Wr, 
				 WoS[i], 
				 Grid.nfaceS(i, j+1));
	} else if (Grid.BCtypeS[i] == BC_RINGLEB_FLOW) {
	  //Wl = RinglebFlow(Wl,Grid.xfaceS(i,j+1));
	  Wl = BC_Characteristic_Pressure(Wr,
					  RinglebFlow(Wr,Grid.xfaceS(i,j+1)), 
					  Grid.nfaceS(i, j+1));
	} else if (Grid.BCtypeS[i] == BC_FROZEN) {
	  // calculate Wl based on the ghost cell reconstruction
	  Wl = PiecewiseLinearSolutionAtLocation(i,j,
						 Grid.xfaceS(i, j+1));
	} else if (Grid.BCtypeS[i] == BC_EXACT_SOLUTION) {
	  // calculate Wl using the exact solution at the interface and avoiding the problem to become ill-posed.
	  if (ExactSoln->IsExactSolutionSet()){
	    Wl = BC_Characteristic_Pressure(Wr,
					    ExactSoln->Solution(Grid.xfaceS(i,j+1).x,Grid.xfaceS(i,j+1).y),
					    Grid.nfaceS(i, j+1));
	  } else {
	    throw runtime_error("Euler2D_Quad_Block::dUdt_Residual_Evaluation() ERROR! There is no exact solution set for the Exact_Solution BC.");
	  }
	} else if (Grid.BCtypeS[i] == BC_WALL_INVISCID) {
	  Wl = Reflect(Wr, Grid.nfaceS(i, j+1));
	  Wi = PiecewiseLinearSolutionAtLocation(i,j,
						 Grid.xfaceS(i, j+1));
	  Wl.d = Wi.d;
	  Wl.p = Wi.p;
	} /* endif */
      } else if (j == JCu && 
		 (Grid.BCtypeN[i] == BC_REFLECTION ||
		  Grid.BCtypeN[i] == BC_CHARACTERISTIC ||
		  Grid.BCtypeN[i] == BC_CHARACTERISTIC_VELOCITY ||
		  Grid.BCtypeN[i] == BC_BURNING_SURFACE ||
		  Grid.BCtypeN[i] == BC_RINGLEB_FLOW || 
		  Grid.BCtypeN[i] == BC_FROZEN ||
		  Grid.BCtypeN[i] == BC_EXACT_SOLUTION ||
		  Grid.BCtypeN[i] == BC_WALL_INVISCID )) {
	dX = Grid.xfaceN(i, j)-Grid.Cell[i][j].Xc;
	Wl = W[i][j] + 
	  (phi[i][j]^dWdx[i][j])*dX.x +
	  (phi[i][j]^dWdy[i][j])*dX.y;
	if (Grid.BCtypeN[i] == BC_REFLECTION) {
	  Wr = Reflect(Wl, Grid.nfaceN(i, j));
	} else if (Grid.BCtypeN[i] == BC_BURNING_SURFACE) {
	  Wr = BurningSurface(Wl, Grid.nfaceN(i, j));
	} else if (Grid.BCtypeS[i] == BC_MASS_INJECTION) {
	  Wr = MassInjection(Wr,Grid.nfaceN(i,j),OFF);
	} else if (Grid.BCtypeN[i] == BC_CHARACTERISTIC) {
	  Wr = BC_Characteristic_Pressure(Wl, 
					  WoN[i], 
					  Grid.nfaceN(i, j));
	} else if (Grid.BCtypeN[i] == BC_CHARACTERISTIC_VELOCITY) {
	  Wr = BC_Characteristic(Wl, 
				 WoN[i], 
				 Grid.nfaceN(i, j));
	} else if (Grid.BCtypeN[i] == BC_RINGLEB_FLOW) {
	  //Wr = RinglebFlow(Wr,Grid.xfaceN(i,j));
	  Wr = BC_Characteristic_Pressure(Wl, 
					  RinglebFlow(Wr,Grid.xfaceN(i,j)), 
					  Grid.nfaceN(i, j));
	} else if (Grid.BCtypeN[i] == BC_FROZEN) {
	  // calculate Wr based on the ghost cell reconstruction
	  Wr = PiecewiseLinearSolutionAtLocation(i,j+1,
						 Grid.xfaceN(i,j));
	} else if (Grid.BCtypeN[i] == BC_EXACT_SOLUTION) {
	  // calculate Wr using the exact solution at the interface and avoiding the problem to become ill-posed.
	  if (ExactSoln->IsExactSolutionSet()){
	    Wr = BC_Characteristic_Pressure(Wl,
					    ExactSoln->Solution(Grid.xfaceN(i,j).x,Grid.xfaceN(i,j).y),
					    Grid.nfaceN(i, j));
	  } else {
	    throw runtime_error("Euler2D_Quad_Block::dUdt_Residual_Evaluation() ERROR! There is no exact solution set for the Exact_Solution BC.");
	  }
	} else if (Grid.BCtypeN[i] == BC_WALL_INVISCID) {
	  Wr = Reflect(Wl, Grid.nfaceN(i, j));
	  Wi = PiecewiseLinearSolutionAtLocation(i,j+1,
						 Grid.xfaceN(i,j));
	  Wr.d = Wi.d;
	  Wr.p = Wi.p;
	} /* endif */
      } else {
	dX = Grid.xfaceN(i, j)-Grid.Cell[i][j].Xc;
	Wl = W[i][j] + 
	  (phi[i][j]^dWdx[i][j])*dX.x +
	  (phi[i][j]^dWdy[i][j])*dX.y;
	dX = Grid.xfaceS(i, j+1)-Grid.Cell[i][j+1].Xc;
	Wr = W[i][j+1] +
	  (phi[i][j+1]^dWdx[i][j+1])*dX.x +
	  (phi[i][j+1]^dWdy[i][j+1])*dX.y;
      } /* endif */
	
      switch(IP.i_Flux_Function) {
      case FLUX_FUNCTION_GODUNOV :
	Flux = FluxGodunov_n(Wl, Wr, Grid.nfaceN(i, j));
	break;
      case FLUX_FUNCTION_ROE :
	Flux = FluxRoe_n(Wl, Wr, Grid.nfaceN(i, j));
	break;
      case FLUX_FUNCTION_RUSANOV :
	Flux = FluxRusanov_n(Wl, Wr, Grid.nfaceN(i, j));
	break;
      case FLUX_FUNCTION_HLLE :
	Flux = FluxHLLE_n(Wl, Wr, Grid.nfaceN(i, j));
	break;
      case FLUX_FUNCTION_LINDE :
	Flux = FluxLinde_n(Wl, Wr, Grid.nfaceN(i, j));
	break;
      case FLUX_FUNCTION_HLLC :
	Flux = FluxHLLC_n(Wl, Wr, Grid.nfaceN(i, j));
	break;
      case FLUX_FUNCTION_VANLEER :
	Flux = FluxVanLeer_n(Wl, Wr, Grid.nfaceN(i, j));
	break;
      case FLUX_FUNCTION_AUSM :
	Flux = FluxAUSM_n(Wl, Wr, Grid.nfaceN(i, j));
	break;
      case FLUX_FUNCTION_AUSMplus :
	Flux = FluxAUSMplus_n(Wl, Wr, Grid.nfaceN(i, j));
	break;
      case FLUX_FUNCTION_ROE_PRECON_WS :
	Flux = FluxRoe_n_Precon_WS(Wl, Wr, Grid.nfaceN(i, j));
	break;
      case FLUX_FUNCTION_HLLE_PRECON_WS :
	Flux = FluxHLLE_n_Precon_WS(Wl, Wr, Grid.nfaceN(i, j));
	break;
      default:
	Flux = FluxRoe_n(Wl, Wr, Grid.nfaceN(i, j));
	break;
      } /* endswitch */
	
      /* Evaluate cell-averaged solution changes. */
	
      dUdt[i][j][0] -= 
	Flux*Grid.lfaceN(i, j)/
	Grid.Cell[i][j].A;
      dUdt[i][j+1][0] += 
	Flux*Grid.lfaceS(i, j+1)/
	Grid.Cell[i][j+1].A;

      /* Save south and north face boundary flux. */
	
      if (j == JCl-1) {
	FluxS[i] = -Flux*Grid.lfaceS(i, j+1);
      } else if (j == JCu) {
	FluxN[i] = Flux*Grid.lfaceN(i, j);
      } /* endif */
	
    } /* endfor */
      
    dUdt[i][JCl-1][0] = Euler2D_U_VACUUM;
    dUdt[i][JCu+1][0] = Euler2D_U_VACUUM;
  } /* endfor */
    
    /* Residual successfully evaluated. */
  return 0;

}

/*********************************************************//**
 * This routine determines the solution residuals for a 
 * given stage of a variety of multi-stage explicit     
 * time integration schemes for a given solution block. 
 *                                                      
 ************************************************************/
int Euler2D_Quad_Block::dUdt_Multistage_Explicit(const int &i_stage,
						 const Euler2D_Input_Parameters &IP){
  
  int i, j, k_residual;
  double omega;
  Vector2D dX;
  Euler2D_pState Wl, Wr, Wi;
  Euler2D_cState Flux;

  /* Evaluate the solution residual for stage 
     i_stage of an N stage scheme. */

  /* Evaluate the time step fraction and residual storage location for the stage. */
    
  switch(IP.i_Time_Integration) {
  case TIME_STEPPING_EXPLICIT_EULER :
    omega = Runge_Kutta(i_stage, IP.N_Stage);
    k_residual = 0;
    break;
  case TIME_STEPPING_EXPLICIT_PREDICTOR_CORRECTOR :
    omega = Runge_Kutta(i_stage, IP.N_Stage);
    k_residual = 0;
    break;
  case TIME_STEPPING_EXPLICIT_RUNGE_KUTTA :
    omega = Runge_Kutta(i_stage, IP.N_Stage);
    k_residual = 0;
    if (IP.N_Stage == 4) {
      if (i_stage == 4) {
	k_residual = 0;
      } else {
	k_residual = i_stage - 1;
      } /* endif */
    } /* endif */
    break;
  case TIME_STEPPING_MULTISTAGE_OPTIMAL_SMOOTHING :
    omega = MultiStage_Optimally_Smoothing(i_stage, 
					   IP.N_Stage,
					   IP.i_Limiter);
    k_residual = 0;
    break;
  default:
    omega = Runge_Kutta(i_stage, IP.N_Stage);
    k_residual = 0;
    break;
  } /* endswitch */
    
    /* Perform the linear reconstruction within each cell
       of the computational grid for this stage. */
    
  switch(IP.i_Reconstruction) {
  case RECONSTRUCTION_GREEN_GAUSS :
    Linear_Reconstruction_GreenGauss(*this,
				     IP.i_Limiter);    
    break;
  case RECONSTRUCTION_LEAST_SQUARES :
    Linear_Reconstruction_LeastSquares(*this,
				       IP.i_Limiter);
    break;
  default:
    Linear_Reconstruction_LeastSquares(*this,
				       IP.i_Limiter);
    break;
  } /* endswitch */

    /* Evaluate the time rate of change of the solution
       (i.e., the solution residuals) using a second-order
       limited upwind scheme with a variety of flux functions. */
    
    // Add i-direction (zeta-direction) fluxes.
  for ( j  = JCl-1 ; j <= JCu+1 ; ++j ) {
    if ( i_stage == 1 ) {
      Uo[ICl-1][j] = U[ICl-1][j];
      dUdt[ICl-1][j][k_residual] = Euler2D_U_VACUUM;
    } else {
      dUdt[ICl-1][j][k_residual] = Euler2D_U_VACUUM;
    } /* endif */
    
    for ( i = ICl-1 ; i <= ICu ; ++i ) {
      if ( i_stage == 1 ) {
	Uo[i+1][j] = U[i+1][j];
	dUdt[i+1][j][k_residual] = Euler2D_U_VACUUM;
      } else if ( j > JCl-1 && j < JCu+1 ) {
	switch(IP.i_Time_Integration) {
	case TIME_STEPPING_EXPLICIT_PREDICTOR_CORRECTOR :
	  //dUdt[i+1][j][k_residual] = 
	  //   dUdt[i+1][j][k_residual];
	  break;
	case TIME_STEPPING_EXPLICIT_RUNGE_KUTTA :
	  if (IP.N_Stage == 2) {
	    //dUdt[i+1][j][k_residual] = 
	    //   dUdt[i+1][j][k_residual];
	  } else if (IP.N_Stage == 4 && i_stage == 4) {
	    dUdt[i+1][j][k_residual] = 
	      dUdt[i+1][j][0] + 
	      TWO*dUdt[i+1][j][1] +
	      TWO*dUdt[i+1][j][2];
	  } else {
	    dUdt[i+1][j][k_residual] = Euler2D_U_VACUUM;
	  } /* endif */
	  break;
	case TIME_STEPPING_MULTISTAGE_OPTIMAL_SMOOTHING :
	  dUdt[i+1][j][k_residual] = Euler2D_U_VACUUM;
	  break;
	default:
	  dUdt[i+1][j][k_residual] = Euler2D_U_VACUUM;
	  break;
	} /* endswitch */
      } /* endif */
    
      if ( j > JCl-1 && j < JCu+1 ) {
    
	/* Evaluate the cell interface i-direction fluxes. */
    
	if (i == ICl-1 && 
	    (Grid.BCtypeW[j] == BC_REFLECTION ||
	     Grid.BCtypeW[j] == BC_CHARACTERISTIC ||
	     Grid.BCtypeW[j] == BC_CHARACTERISTIC_VELOCITY ||
	     Grid.BCtypeW[j] == BC_BURNING_SURFACE ||
	     Grid.BCtypeW[j] == BC_MASS_INJECTION ||
	     Grid.BCtypeW[j] == BC_RINGLEB_FLOW ||
	     Grid.BCtypeW[j] == BC_FROZEN ||
	     Grid.BCtypeW[j] == BC_EXACT_SOLUTION ||
	     Grid.BCtypeW[j] == BC_WALL_INVISCID)) {
	  dX = Grid.xfaceW(i+1, j)-Grid.Cell[i+1][j].Xc;
	  Wr = W[i+1][j] + 
	    (phi[i+1][j]^dWdx[i+1][j])*dX.x +
	    (phi[i+1][j]^dWdy[i+1][j])*dX.y;
	  if (Grid.BCtypeW[j] == BC_REFLECTION) {
	    Wl = Reflect(Wr, Grid.nfaceW(i+1, j));
	  } else if (Grid.BCtypeW[j] == BC_BURNING_SURFACE) {
	    Wl = BurningSurface(Wr, Grid.nfaceW(i+1, j));
	  } else if (Grid.BCtypeW[j] == BC_MASS_INJECTION) {
	    Wl = MassInjection(Wr,Grid.nfaceW(i+1,j),OFF);
	  } else if (Grid.BCtypeW[j] == BC_CHARACTERISTIC) {
	    Wl = BC_Characteristic_Pressure(Wr, 
					    WoW[j], 
					    Grid.nfaceW(i+1, j));
	  } else if (Grid.BCtypeW[j] == BC_CHARACTERISTIC_VELOCITY) {
	    Wl = BC_Characteristic(Wr, 
				   WoW[j], 
				   Grid.nfaceW(i+1, j));
	  } else if (Grid.BCtypeW[j] == BC_RINGLEB_FLOW) {
	    Wl = RinglebFlow(Wl,Grid.xfaceW(i+1,j));
	    //Wl = Reflect(Wr, Grid.nfaceW(i+1, j));
	  } else if (Grid.BCtypeW[j] == BC_FROZEN) {
	    // calculate Wl based on the ghost cell reconstruction
	    Wl = PiecewiseLinearSolutionAtLocation(i,j,
						   Grid.xfaceW(i+1,j));
	  } else if (Grid.BCtypeW[j] == BC_EXACT_SOLUTION) {
	    // calculate Wl using the exact solution at the interface and avoiding the problem to become ill-posed.
	    if (ExactSoln->IsExactSolutionSet()){
	      Wl = BC_Characteristic_Pressure(Wr,
					      ExactSoln->Solution(Grid.xfaceW(i+1,j).x,Grid.xfaceW(i+1,j).y),
					      Grid.nfaceW(i+1, j));
	    } else {
	      throw runtime_error("Euler2D_Quad_Block::dUdt_Residual_Evaluation() ERROR! There is no exact solution set for the Exact_Solution BC.");
	    }
	  } else if (Grid.BCtypeW[j] == BC_WALL_INVISCID) {
	    Wl = Reflect(Wr, Grid.nfaceW(i+1, j));
	    Wi = PiecewiseLinearSolutionAtLocation(i,j,
						  Grid.xfaceW(i+1,j));
	    Wl.d = Wi.d;
	    Wl.p = Wi.p;
	  } /* endif */
	} else if (i == ICu && 
		   (Grid.BCtypeE[j] == BC_REFLECTION ||
		    Grid.BCtypeE[j] == BC_CHARACTERISTIC ||
		    Grid.BCtypeE[j] == BC_CHARACTERISTIC_VELOCITY ||
		    Grid.BCtypeE[j] == BC_BURNING_SURFACE ||
		    Grid.BCtypeE[j] == BC_MASS_INJECTION ||
		    Grid.BCtypeE[j] == BC_RINGLEB_FLOW ||
		    Grid.BCtypeE[j] == BC_FROZEN ||
		    Grid.BCtypeE[j] == BC_EXACT_SOLUTION ||
		    Grid.BCtypeE[j] == BC_WALL_INVISCID )) {
	  dX = Grid.xfaceE(i, j)-Grid.Cell[i][j].Xc;
	  Wl = W[i][j] + 
	    (phi[i][j]^dWdx[i][j])*dX.x +
	    (phi[i][j]^dWdy[i][j])*dX.y;
	  if (Grid.BCtypeE[j] == BC_REFLECTION) {
	    Wr = Reflect(Wl, Grid.nfaceE(i, j));
	  } else if (Grid.BCtypeE[j] == BC_BURNING_SURFACE) {
	    Wr = BurningSurface(Wl, Grid.nfaceE(i, j));
	  } else if (Grid.BCtypeE[j] == BC_MASS_INJECTION) {
	    Wr = MassInjection(Wl,Grid.nfaceE(i,j),OFF);
	  } else if (Grid.BCtypeE[j] == BC_CHARACTERISTIC) {
	    Wr = BC_Characteristic_Pressure(Wl, 
					    WoE[j], 
					    Grid.nfaceE(i, j));
	  } else if (Grid.BCtypeE[j] == BC_CHARACTERISTIC_VELOCITY) {
	    Wr = BC_Characteristic(Wl, 
				   WoE[j], 
				   Grid.nfaceE(i, j));
	  } else if (Grid.BCtypeE[j] == BC_RINGLEB_FLOW) {
	    Wr = RinglebFlow(Wr,Grid.xfaceE(i,j));
	    //Wr = Reflect(Wl, Grid.nfaceE(i, j));
	  } else if (Grid.BCtypeE[j] == BC_FROZEN) {
	    // calculate Wr based on the ghost cell reconstruction
	    Wr = PiecewiseLinearSolutionAtLocation(i+1,j,
						   Grid.xfaceE(i,j));
	  } else if (Grid.BCtypeE[j] == BC_EXACT_SOLUTION) {
	    // calculate Wr using the exact solution at the interface and avoiding the problem to become ill-posed.
	    if (ExactSoln->IsExactSolutionSet()){
	      Wr = BC_Characteristic_Pressure(Wl,
					      ExactSoln->Solution(Grid.xfaceE(i,j).x,Grid.xfaceE(i,j).y),
					      Grid.nfaceE(i, j));
	    } else {
	      throw runtime_error("Euler2D_Quad_Block::dUdt_Residual_Evaluation() ERROR! There is no exact solution set for the Exact_Solution BC.");
	    }	    
	  } else if (Grid.BCtypeE[j] == BC_WALL_INVISCID) {
	    Wr = Reflect(Wl, Grid.nfaceE(i, j));
	    Wi = PiecewiseLinearSolutionAtLocation(i+1,j,
						   Grid.xfaceW(i+1, j));
	    Wl.d = Wi.d;
	    Wl.p = Wi.p;
	  } /* endif */
	} else {            
	  dX = Grid.xfaceE(i, j)-Grid.Cell[i][j].Xc;
	  Wl = W[i][j] + 
	    (phi[i][j]^dWdx[i][j])*dX.x +
	    (phi[i][j]^dWdy[i][j])*dX.y;
	  dX = Grid.xfaceW(i+1, j)-Grid.Cell[i+1][j].Xc;
	  Wr = W[i+1][j] + 
	    (phi[i+1][j]^dWdx[i+1][j])*dX.x +
	    (phi[i+1][j]^dWdy[i+1][j])*dX.y;
	} /* endif */

	switch(IP.i_Flux_Function) {
	case FLUX_FUNCTION_GODUNOV :
	  Flux = FluxGodunov_n(Wl, Wr, Grid.nfaceE(i, j));
	  break;
	case FLUX_FUNCTION_ROE :
	  Flux = FluxRoe_n(Wl, Wr, Grid.nfaceE(i, j));
	  break;
	case FLUX_FUNCTION_RUSANOV :
	  Flux = FluxRusanov_n(Wl, Wr, Grid.nfaceE(i, j));
	  break;
	case FLUX_FUNCTION_HLLE :
	  Flux = FluxHLLE_n(Wl, Wr, Grid.nfaceE(i, j));
	  break;
	case FLUX_FUNCTION_LINDE :
	  Flux = FluxLinde_n(Wl, Wr, Grid.nfaceE(i, j));
	  break;
	case FLUX_FUNCTION_HLLC :
	  Flux = FluxHLLC_n(Wl, Wr, Grid.nfaceE(i, j));
	  break;
	case FLUX_FUNCTION_VANLEER :
	  Flux = FluxVanLeer_n(Wl, Wr, Grid.nfaceE(i,j));
	  break;
	case FLUX_FUNCTION_AUSM :
	  Flux = FluxAUSM_n(Wl, Wr, Grid.nfaceE(i,j));
	  break;
	case FLUX_FUNCTION_AUSMplus :
	  Flux = FluxAUSMplus_n(Wl, Wr, Grid.nfaceE(i,j));
	  break;
	case FLUX_FUNCTION_ROE_PRECON_WS :
	  Flux = FluxRoe_n_Precon_WS(Wl, Wr, Grid.nfaceE(i, j));
	  break;
	case FLUX_FUNCTION_HLLE_PRECON_WS :
	  Flux = FluxHLLE_n_Precon_WS(Wl, Wr, Grid.nfaceE(i, j));
	  break;
	default:
	  Flux = FluxRoe_n(Wl, Wr, Grid.nfaceE(i, j));
	  break;
	} /* endswitch */
    
	/* Evaluate cell-averaged solution changes. */
    
	dUdt[i][j][k_residual] -= 
	  (IP.CFL_Number*dt[i][j])*
	  Flux*Grid.lfaceE(i, j)/
	  Grid.Cell[i][j].A;
	dUdt[i+1][j][k_residual] += 
	  (IP.CFL_Number*dt[i+1][j])*
	  Flux*Grid.lfaceW(i+1, j)/
	  Grid.Cell[i+1][j].A;

	/* Include axisymmetric source terms as required. */

	if (Axisymmetric) {
	  dUdt[i][j][k_residual] += 
	    (IP.CFL_Number*dt[i][j])*
	    S(W[i][j], Grid.Cell[i][j].Xc);
	} /* endif */

	/* Save west and east face boundary flux. */

	if (i == ICl-1) {
	  FluxW[j] = -Flux*Grid.lfaceW(i+1, j);
	} else if (i == ICu) {
	  FluxE[j] = Flux*Grid.lfaceE(i, j);
	} /* endif */ 

      } /* endif */
    } /* endfor */
    
    if ( j > JCl-1 && j < JCu+1 ) {
      dUdt[ICl-1][j][k_residual] = Euler2D_U_VACUUM;
      dUdt[ICu+1][j][k_residual] = Euler2D_U_VACUUM;
    } /* endif */
  } /* endfor */
    
    // Add j-direction (eta-direction) fluxes.
  for ( i = ICl ; i <= ICu ; ++i ) {
    for ( j  = JCl-1 ; j <= JCu ; ++j ) {
    
      /* Evaluate the cell interface j-direction fluxes. */
         
      if (j == JCl-1 && 
	  (Grid.BCtypeS[i] == BC_REFLECTION ||
	   Grid.BCtypeS[i] == BC_CHARACTERISTIC ||
	   Grid.BCtypeS[i] == BC_CHARACTERISTIC_VELOCITY ||
	   Grid.BCtypeS[i] == BC_BURNING_SURFACE ||
	   Grid.BCtypeS[i] == BC_MASS_INJECTION ||
	   Grid.BCtypeS[i] == BC_RINGLEB_FLOW ||
	   Grid.BCtypeS[i] == BC_FROZEN ||
	   Grid.BCtypeS[i] == BC_EXACT_SOLUTION ||
	   Grid.BCtypeS[i] == BC_WALL_INVISCID )) {
	dX = Grid.xfaceS(i, j+1)-Grid.Cell[i][j+1].Xc;
	Wr = W[i][j+1] +
	  (phi[i][j+1]^dWdx[i][j+1])*dX.x +
	  (phi[i][j+1]^dWdy[i][j+1])*dX.y;
	if (Grid.BCtypeS[i] == BC_REFLECTION) {
	  Wl = Reflect(Wr, Grid.nfaceS(i, j+1));
	} else if (Grid.BCtypeS[i] == BC_BURNING_SURFACE) {
	  Wl = BurningSurface(Wr, Grid.nfaceS(i, j+1));
	} else if (Grid.BCtypeS[i] == BC_MASS_INJECTION) {
	  Wl = MassInjection(Wr,Grid.nfaceS(i,j+1),OFF);
	} else if (Grid.BCtypeS[i] == BC_CHARACTERISTIC) {
	  Wl = BC_Characteristic_Pressure(Wr, 
					  WoS[i], 
					  Grid.nfaceS(i, j+1));
	} else if (Grid.BCtypeS[i] == BC_CHARACTERISTIC_VELOCITY) {
	  Wl = BC_Characteristic(Wr, 
				 WoS[i], 
				 Grid.nfaceS(i, j+1));
	} else if (Grid.BCtypeS[i] == BC_RINGLEB_FLOW) {
	  //Wl = RinglebFlow(Grid.xfaceS(i,j+1));
	  Wl = BC_Characteristic_Pressure(Wr,
					  RinglebFlow(Wl,Grid.xfaceS(i,j+1)), 
					  Grid.nfaceS(i, j+1));
	} else if (Grid.BCtypeS[i] == BC_FROZEN) {
	  // calculate Wl based on the ghost cell reconstruction
	  Wl = PiecewiseLinearSolutionAtLocation(i,j,
						 Grid.xfaceS(i, j+1));
	} else if (Grid.BCtypeS[i] == BC_EXACT_SOLUTION) {
	  // calculate Wl using the exact solution at the interface and avoiding the problem to become ill-posed.
	  if (ExactSoln->IsExactSolutionSet()){
	    Wl = BC_Characteristic_Pressure(Wr,
					    ExactSoln->Solution(Grid.xfaceS(i,j+1).x,Grid.xfaceS(i,j+1).y),
					    Grid.nfaceS(i, j+1));
	  } else {
	    throw runtime_error("Euler2D_Quad_Block::dUdt_Residual_Evaluation() ERROR! There is no exact solution set for the Exact_Solution BC.");
	  }
	} else if (Grid.BCtypeS[i] == BC_WALL_INVISCID) {
	  Wl = Reflect(Wr, Grid.nfaceS(i, j+1));
	  Wi = PiecewiseLinearSolutionAtLocation(i,j,
						 Grid.xfaceS(i, j+1));
	  Wl.d = Wi.d;
	  Wl.p = Wi.p;
	} /* endif */
      } else if (j == JCu && 
		 (Grid.BCtypeN[i] == BC_REFLECTION ||
		  Grid.BCtypeN[i] == BC_CHARACTERISTIC ||
		  Grid.BCtypeN[i] == BC_CHARACTERISTIC_VELOCITY ||
		  Grid.BCtypeN[i] == BC_BURNING_SURFACE ||
		  Grid.BCtypeN[i] == BC_MASS_INJECTION ||
		  Grid.BCtypeN[i] == BC_RINGLEB_FLOW || 
		  Grid.BCtypeN[i] == BC_FROZEN ||
		  Grid.BCtypeN[i] == BC_EXACT_SOLUTION ||
		  Grid.BCtypeN[i] == BC_WALL_INVISCID )) {
	dX = Grid.xfaceN(i, j)-Grid.Cell[i][j].Xc;
	Wl = W[i][j] + 
	  (phi[i][j]^dWdx[i][j])*dX.x +
	  (phi[i][j]^dWdy[i][j])*dX.y;
	if (Grid.BCtypeN[i] == BC_REFLECTION) {
	  Wr = Reflect(Wl, Grid.nfaceN(i, j));
	} else if (Grid.BCtypeN[i] == BC_BURNING_SURFACE) {
	  Wr = BurningSurface(Wl, Grid.nfaceN(i, j));
	} else if (Grid.BCtypeN[i] == BC_MASS_INJECTION) {
	  Wr = MassInjection(Wl,Grid.nfaceN(i,j),OFF);
	} else if (Grid.BCtypeN[i] == BC_CHARACTERISTIC) {
	  Wr = BC_Characteristic_Pressure(Wl, 
					  WoN[i], 
					  Grid.nfaceN(i, j));
	} else if (Grid.BCtypeN[i] == BC_CHARACTERISTIC_VELOCITY) {
	  Wr = BC_Characteristic(Wl, 
				 WoN[i], 
				 Grid.nfaceN(i, j));
	} else if (Grid.BCtypeN[i] == BC_RINGLEB_FLOW) {
	  //Wr = RinglebFlow(Grid.xfaceN(i,j));
	  Wr = BC_Characteristic_Pressure(Wl, 
					  RinglebFlow(Wr,Grid.xfaceN(i,j)), 
					  Grid.nfaceN(i, j));
	} else if (Grid.BCtypeN[i] == BC_FROZEN) {
	  // calculate Wr based on the ghost cell reconstruction
	  Wr = PiecewiseLinearSolutionAtLocation(i,j+1,
						 Grid.xfaceN(i,j));
	} else if (Grid.BCtypeN[i] == BC_EXACT_SOLUTION) {
	  // calculate Wr using the exact solution at the interface and avoiding the problem to become ill-posed.
	  if (ExactSoln->IsExactSolutionSet()){
	    Wr = BC_Characteristic_Pressure(Wl,
					    ExactSoln->Solution(Grid.xfaceN(i,j).x,Grid.xfaceN(i,j).y),
					    Grid.nfaceN(i, j));
	  } else {
	    throw runtime_error("Euler2D_Quad_Block::dUdt_Residual_Evaluation() ERROR! There is no exact solution set for the Exact_Solution BC.");
	  }
	} else if (Grid.BCtypeN[i] == BC_WALL_INVISCID) {
	  Wr = Reflect(Wl, Grid.nfaceN(i, j));
	  Wi = PiecewiseLinearSolutionAtLocation(i,j+1,
						 Grid.xfaceN(i,j));
	  Wr.d = Wi.d;
	  Wr.p = Wi.p;
	} /* endif */
      } else {
	dX = Grid.xfaceN(i, j)-Grid.Cell[i][j].Xc;
	Wl = W[i][j] + 
	  (phi[i][j]^dWdx[i][j])*dX.x +
	  (phi[i][j]^dWdy[i][j])*dX.y;
	dX = Grid.xfaceS(i, j+1)-Grid.Cell[i][j+1].Xc;
	Wr = W[i][j+1] +
	  (phi[i][j+1]^dWdx[i][j+1])*dX.x +
	  (phi[i][j+1]^dWdy[i][j+1])*dX.y;
      } /* endif */

      switch(IP.i_Flux_Function) {
      case FLUX_FUNCTION_GODUNOV :
	Flux = FluxGodunov_n(Wl, Wr, Grid.nfaceN(i, j));
	break;
      case FLUX_FUNCTION_ROE :
	Flux = FluxRoe_n(Wl, Wr, Grid.nfaceN(i, j));
	break;
      case FLUX_FUNCTION_RUSANOV :
	Flux = FluxRusanov_n(Wl, Wr, Grid.nfaceN(i, j));
	break;
      case FLUX_FUNCTION_HLLE :
	Flux = FluxHLLE_n(Wl, Wr, Grid.nfaceN(i, j));
	break;
      case FLUX_FUNCTION_LINDE :
	Flux = FluxLinde_n(Wl, Wr, Grid.nfaceN(i, j));
	break;
      case FLUX_FUNCTION_HLLC :
	Flux = FluxHLLC_n(Wl, Wr, Grid.nfaceN(i, j));
	break;
      case FLUX_FUNCTION_VANLEER :
	Flux = FluxVanLeer_n(Wl, Wr, Grid.nfaceN(i, j));
	break;
      case FLUX_FUNCTION_AUSM :
	Flux = FluxAUSM_n(Wl, Wr, Grid.nfaceN(i, j));
	break;
      case FLUX_FUNCTION_AUSMplus :
	Flux = FluxAUSMplus_n(Wl, Wr, Grid.nfaceN(i, j));
	break;
      case FLUX_FUNCTION_ROE_PRECON_WS :
	Flux = FluxRoe_n_Precon_WS(Wl, Wr, Grid.nfaceN(i, j));
	break;
      case FLUX_FUNCTION_HLLE_PRECON_WS :
	Flux = FluxHLLE_n_Precon_WS(Wl, Wr, Grid.nfaceN(i, j));
	break;
      default:
	Flux = FluxRoe_n(Wl, Wr, Grid.nfaceN(i, j));
	break;
      } /* endswitch */
    
      /* Evaluate cell-averaged solution changes. */
    
      dUdt[i][j][k_residual] -= 
	(IP.CFL_Number*dt[i][j])*
	Flux*Grid.lfaceN(i, j)/
	Grid.Cell[i][j].A;
      dUdt[i][j+1][k_residual] += 
	(IP.CFL_Number*dt[i][j+1])*
	Flux*Grid.lfaceS(i, j+1)/
	Grid.Cell[i][j+1].A;

      /* Save south and north face boundary flux. */

      if (j == JCl-1) {
	FluxS[i] = -Flux*Grid.lfaceS(i, j+1);
      } else if (j == JCu) {
	FluxN[i] = Flux*Grid.lfaceN(i, j);
      } /* endif */

    } /* endfor */
    
    dUdt[i][JCl-1][k_residual] = Euler2D_U_VACUUM;
    dUdt[i][JCu+1][k_residual] = Euler2D_U_VACUUM;
  } /* endfor */
    
    /* Residual for the stage successfully calculated. */

  return (0);

}
