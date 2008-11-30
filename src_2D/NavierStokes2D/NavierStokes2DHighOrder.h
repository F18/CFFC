/*!\file NavierStokes2DHighOrder.h
  \brief Header file defining specializations of HighOrder2D template class for NavierStokes2D solution states. */

#ifndef _NAVIERSTOKES2D_HIGHORDER_SPECIALIZATION_INCLUDED
#define _NAVIERSTOKES2D_HIGHORDER_SPECIALIZATION_INCLUDED

/* Include required C++ libraries. */
// None

/* Using std namespace functions */
// None

/* Include CFFC header files */
#include "../HighOrderReconstruction/HighOrder2D.h" /* Include high-order class header file */

//////////////////////////////////////////////////////////
// Specialized routines for NavierStokes2D solution states     //
//////////////////////////////////////////////////////////

/*!
 * Get the number of variables in the solution state.
 * To avoid reconstruction of unused variables, use
 * the number of active variables, which depends on the flow type.
 */
template<> inline
int HighOrder2D<NavierStokes2D_pState>::NumberOfVariables(void) const {
  return NavierStokes2D_pState::NumVarActive();
}


/*!
 * Set the constraints equations in the constrained least-square reconstruction.
 * These constraints are called relational because they put in relationship
 * several solution parameters.
 * The starting position for the entries in the LHS and RHS of the linear
 * system are specified by the (StartRow, StartCol) input parameters.   \n                                
 *                                                                          
 * The data passed:
 *  GQP      -- Gauss Quadrature Point Locations (i.e. flux calculation points)
 *  A_GQP    -- The coefficient for the Dirichlet boundary condition at GQP                                                
 *  B_GQP    -- The coefficient for the Neumann boundary condition at GQP                                                
 *  F1(GQP)  -- The value of the Dirichlet boundary condition at GQP  
 *  F2(GQP)  -- The value of the Neumann boundary condition at GQP    
 *
 * \param SolnBlk the quad block for which the solution reconstruction is done.
 * \param iCell i-index of the reconstructed cell
 * \param jCell j-index of the reconstructed cell
 * \param Constraints_Loc GQP array
 * \param Constraints_Normals normal vectors at GQP locations
 * \param Constraints_BCs provide the boundary condition coefficients (i.e. A_GQP, B_GQP, F1, F2)
 * \param BC_Type the type of the boundary condition for which the relational constraint is built
 * \param ParameterIndex related to the indexes of the solution
 * \param A the LHS assemble matrix 
 * \param B the RHS assemble matrix
 *
 */
template<>
template<> inline
void HighOrder2D<NavierStokes2D_pState>::
Generalized_RelationalConstraints_Equations(NavierStokes2D_Quad_Block & SolnBlk,
					    const int &iCell, const int &jCell,
					    Vector2DArray & Constraints_Loc,
					    Vector2DArray & Constraints_Normals,
					    BC_Type_Array & Constraints_BCs,
					    const int & BC_Type,
					    DenseMatrix & A, DenseMatrix & All_U,
					    const IndexType & ParameterIndex,
					    const int &StartRow, const int &StartCol) {  

  int P1, P2, i;
  int Point;
  double PowXC, PowYC;		/* PowXC = DistXi^(P1-1); PowYC = DistYi^(P2-1) */
  double DistXi, DistYi;
  int IndexP1, IndexP2;
  double GenTerm;

  
  if (SolnBlk.Flow_Type == FLOWTYPE_INVISCID){
    // impose relational constraints for inviscid flow

    switch(BC_Type){
    case BC_REFLECTION:		// treat it in the same way as BC_WALL_INVISCID
    case BC_WALL_INVISCID:
      /* compute for each derivative the corresponding entry in the matrix of the linear system for this type of BCs */
      for(Point = 0; Point < Constraints_Loc.size(); ++Point){

	// Determine distance between the current GQP and the centroid of cell (iCell,jCell)
	DistXi = Constraints_Loc[Point].x - XCellCenter(iCell,jCell);
	DistYi = Constraints_Loc[Point].y - YCellCenter(iCell,jCell);

	// Form the LHS  -- build the row of the matrix A associated with the current Point
	for (i=0; i<=CellTaylorDeriv(iCell,jCell).LastElem(); ++i){
	  // build the row of the matrix
	  P1 = CellTaylorDeriv(iCell,jCell,i).P1();  // identify P1
	  P2 = CellTaylorDeriv(iCell,jCell,i).P2();  // identify P2

	  /* Initialize PowXC & PowYC */
	  PowXC = 1.0/DistXi;
	  PowYC = 1.0/DistYi;

	  /* Update PowXC & PowYC */
	  for (IndexP1 = 1; IndexP1 <= P1; ++IndexP1){ PowXC *= DistXi; }
	  for (IndexP2 = 1; IndexP2 <= P2; ++IndexP2){ PowYC *= DistYi; }

	  /* Impose normal velocity to be zero --> V*n =0 */
	  GenTerm = PowXC*PowYC*DistXi*DistYi;

	  /* set the u-velocity entry */
	  A(StartRow + Point,StartCol + i) = GenTerm*(Constraints_Normals[Point].x);

	  /* set the v-velocity entry */
	  A(StartRow + Point,StartCol + NumberOfTaylorDerivatives() + i) = GenTerm*(Constraints_Normals[Point].y);
	}//endfor (i)

	// Form the RHS  -- build the row of the matrix All_U associated with the current GQP
	/* For BC_WALL_INVISCID the value gets set to zero */
	All_U(StartRow + Point,0) = 0.0;

      }//endfor (Point)
    
      break;
    } // endswitch

  } // endif (Flow_Type)
}



#endif	// _NAVIERSTOKES2D_HIGHORDER_SPECIALIZATION_INCLUDED
