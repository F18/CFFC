/**********************************************************************
 * EmbeddedBoundaries2D: Header file declaring 2D embedded boundary   *
 *                       classes and functions.                       *
 **********************************************************************/

#ifndef _EMBEDDEDBOUNDARIES2D_EULER2D_INCLUDED
#define _EMBEDDEDBOUNDARIES2D_EULER2D_INCLUDED

// Include 2D Euler quadrilateral mesh solution header file.

#ifndef _EULER2D_QUAD_INCLUDED
#include "../Euler2D/Euler2DQuad.h"
#endif // _EULER2D_QUAD_INCLUDED

// Include Embedded Boundaries input header file.

#ifndef _EMBEDDEDBOUNDARIES2D_INCLUDED
#include "EmbeddedBoundaries2D.h"
#endif // _EMBEDDEDBOUNDARIES2D_INCLUDED

/**********************************************************************
 * EmbeddedBoundaries2D::Calculate_Refinement_Criteria --             *
 **********************************************************************/
template <> inline void EmbeddedBoundaries2D<Euler2D_cState,
					     Euler2D_pState,
					     Euler2D_Quad_Block,
					     Euler2D_Input_Parameters>::
Calculate_Refinement_Criteria(const int &nb, const int &i, const int &j,
			      double *refinement_criteria) {

  double grad_rho_x, grad_rho_y, grad_rho_abs, grad_rho_criteria,
         div_V, div_V_criteria, curl_V_z, curl_V_abs, curl_V_criteria;
  int refinement_criteria_number;

  // Refinement criteria that can be used:
  // (1) Refinement criteria based on the gradient of the density field;
  // (2) Refinement criteria based on the divergence of the velocity vector;
  // (3) Refinement criteria based on the curl of the velocity vector;

  // Reconstruct the solution within the cell.
  Linear_Least_Squares_Reconstruction(nb,i,j);

  // Evaluate refinement criteria #1 based on the gradient of the 
  // density field.
  if (IP->Refinement_Criteria_Gradient_Density) {
    grad_rho_x = Local_SolnBlk[nb].dWdx[i][j].d;
    grad_rho_y = Local_SolnBlk[nb].dWdy[i][j].d;
    grad_rho_abs = sqrt(sqr(grad_rho_x) + sqr(grad_rho_y));
    grad_rho_criteria = sqrt(Local_SolnBlk[nb].Grid.Cell[i][j].A)*grad_rho_abs/Local_SolnBlk[nb].W[i][j].d;
  } else {
    grad_rho_criteria = ZERO;
  }

  // Evaluate refinement criteria #2 based on the divergence of the
  // velocity vector.
  if (IP->Refinement_Criteria_Divergence_Velocity) {
    div_V = Local_SolnBlk[nb].dWdx[i][j].v.x + Local_SolnBlk[nb].dWdy[i][j].v.y;
    div_V_criteria = sqrt(Local_SolnBlk[nb].Grid.Cell[i][j].A)*fabs(div_V)/Local_SolnBlk[nb].W[i][j].a();
  } else {
    div_V_criteria = ZERO;
  }

  // Evaluate refinement criteria #3 based on the curl of the
  // velocity vector.
  if (IP->Refinement_Criteria_Curl_Velocity) {
    curl_V_z = Local_SolnBlk[nb].dWdx[i][j].v.y - Local_SolnBlk[nb].dWdy[i][j].v.x; 
    curl_V_abs = sqrt(sqr(curl_V_z)); 
    curl_V_criteria = sqrt(Local_SolnBlk[nb].Grid.Cell[i][j].A)*curl_V_abs/Local_SolnBlk[nb].W[i][j].a();
  } else {
    curl_V_criteria = ZERO;
  }

  // Set the refinement criteria.
  refinement_criteria_number = 0;
  if (IP->Refinement_Criteria_Gradient_Density) {
    refinement_criteria[refinement_criteria_number] = max(refinement_criteria[refinement_criteria_number],
							  grad_rho_criteria);
    refinement_criteria_number++;
  }
  if (IP->Refinement_Criteria_Divergence_Velocity) {
    refinement_criteria[refinement_criteria_number] = max(refinement_criteria[refinement_criteria_number],
							  div_V_criteria);
    refinement_criteria_number++;
  }
  if (IP->Refinement_Criteria_Curl_Velocity) {
    refinement_criteria[refinement_criteria_number] = max(refinement_criteria[refinement_criteria_number],
							  curl_V_criteria);
    refinement_criteria_number++;
  }

}

/**********************************************************************
 * EmbeddedBoundaries2D::Output_Flat_Plate_Tecplot --                 *
 **********************************************************************/
template <> inline void EmbeddedBoundaries2D<Euler2D_cState,
					     Euler2D_pState,
					     Euler2D_Quad_Block,
					     Euler2D_Input_Parameters>::
Output_Flat_Plate_Tecplot(const int &nb,
			  const int &Output_Title_Soln,
			  ostream &Out_File_Soln,
			  const int &Output_Title_Skin,
			  ostream &Out_File_Skin,
			  double &l1_norm,
			  double &l2_norm,
			  double &max_norm,
			  double &area,
			  int &numberofcells,
			  double &l1_norm_cf,
			  double &l2_norm_cf,
			  double &max_norm_cf,
			  double &area_cf,
			  int &numberofcells_cf) {
}

/**********************************************************************
 * EmbeddedBoundaries2D::BCs_Interface --                             *
 **********************************************************************/
template <> inline int EmbeddedBoundaries2D<Euler2D_cState,
					    Euler2D_pState,
					    Euler2D_Quad_Block,
					    Euler2D_Input_Parameters>::
BCs_Interface(const int &nb, const double &Time) {

  // Exit immediately if none of the embedded boundaries are present in
  // the current solution block.
  if (!Adjustment_Data[nb].Interface_Present[0]) return 0;

  int error_flag, Ni, neighbour_flag;
  int Interface_BC_Type;
  double length;
  Vector2D V;

  for (int j = Local_SolnBlk[nb].JCl-Local_SolnBlk[nb].Nghost; j <= Local_SolnBlk[nb].JCu+Local_SolnBlk[nb].Nghost; j++) {
    for (int i = Local_SolnBlk[nb].ICl-Local_SolnBlk[nb].Nghost; i <= Local_SolnBlk[nb].ICu+Local_SolnBlk[nb].Nghost; i++) {

      // Only apply embedded interface boundary data if the current cell
      // is not active (internal to one of the embedded interfaces) and
      // one of its neighbours is active.  Note that the inactive cell
      // status number corresponds to the interface number (1..Ni) the
      // cell is internal to.
      if (Mesh[nb].cell_status[i][j] != CELL_STATUS_ACTIVE) {

	// Store the interface number.
	Ni = Mesh[nb].cell_status[i][j];
	if (Ni < 1 || Ni > Interface_Union_List.Ni) return 4501;

	// Reset neighbour flag, primitive variable solution state, and
	// the total face length.
	neighbour_flag = OFF;
	Local_SolnBlk[nb].W[i][j].Vacuum();
	length = ZERO;

	// Apply boundary condition for any active neighbour cells if
	// required.

	// NORTH face.
	if (j < Local_SolnBlk[nb].JCu+Local_SolnBlk[nb].Nghost) {
	  if (Mesh[nb].cell_status[i][j+1] == CELL_STATUS_ACTIVE &&
	      Local_SolnBlk[nb].Grid.lfaceN(i,j) > TOLER) {
	    length += Local_SolnBlk[nb].Grid.lfaceN(i,j);
	    // Determine the boundary condition type.
	    Interface_BC_Type = Interface_Union_List[Ni].Determine_Interface_BC_Type(Local_SolnBlk[nb].Grid.xfaceN(i,j));
	    // Determine the boundary velocity.
	    V = Interface_Union_List[Ni].Determine_Interface_Velocity(Local_SolnBlk[nb].Grid.xfaceN(i,j),Time);
	    // Apply the boundary condition.
	    switch(Interface_BC_Type) {
	    case INTERFACE_BC_REFLECTION :
	      Local_SolnBlk[nb].W[i][j] += Reflect(Local_SolnBlk[nb].W[i][j+1],-Local_SolnBlk[nb].Grid.nfaceN(i,j),V)*Local_SolnBlk[nb].Grid.lfaceN(i,j);
	      break;
	    case INTERFACE_BC_BURNING_SURFACE :
	      Local_SolnBlk[nb].W[i][j] += BurningSurface(Local_SolnBlk[nb].W[i][j+1],-Local_SolnBlk[nb].Grid.nfaceN(i,j))*Local_SolnBlk[nb].Grid.lfaceN(i,j);
	      break;
	    case INTERFACE_BC_RINGLEB :
	      Local_SolnBlk[nb].W[i][j] += RinglebFlowAverageState(Local_SolnBlk[nb].W[i][j],
								   Local_SolnBlk[nb].Grid.nodeSW(i,j).X,
								   Local_SolnBlk[nb].Grid.nodeSE(i,j).X,
								   Local_SolnBlk[nb].Grid.nodeNE(i,j).X,
								   Local_SolnBlk[nb].Grid.nodeNW(i,j).X)*Local_SolnBlk[nb].Grid.lfaceN(i,j);
	      break;
	    };
	    // Update neighbour_flag.
	    neighbour_flag = ON;
	  }
	}

	// SOUTH face.
	if (j > Local_SolnBlk[nb].JCl-Local_SolnBlk[nb].Nghost) {
	  if (Mesh[nb].cell_status[i][j-1] == CELL_STATUS_ACTIVE &&
	      Local_SolnBlk[nb].Grid.lfaceS(i,j) > TOLER) {
	    length += Local_SolnBlk[nb].Grid.lfaceS(i,j);
	    // Determine the boundary condition type.
	    Interface_BC_Type = Interface_Union_List[Ni].Determine_Interface_BC_Type(Local_SolnBlk[nb].Grid.xfaceS(i,j));
	    // Determine the boundary velocity.
	    V = Interface_Union_List[Ni].Determine_Interface_Velocity(Local_SolnBlk[nb].Grid.xfaceS(i,j),Time);
	    // Apply the boundary condition.
	    switch(Interface_BC_Type) {
	    case INTERFACE_BC_REFLECTION :
	      Local_SolnBlk[nb].W[i][j] += Reflect(Local_SolnBlk[nb].W[i][j-1],-Local_SolnBlk[nb].Grid.nfaceS(i,j),V)*Local_SolnBlk[nb].Grid.lfaceS(i,j);
	      break;
	    case INTERFACE_BC_BURNING_SURFACE :
	      Local_SolnBlk[nb].W[i][j] += BurningSurface(Local_SolnBlk[nb].W[i][j-1],-Local_SolnBlk[nb].Grid.nfaceS(i,j))*Local_SolnBlk[nb].Grid.lfaceS(i,j);
	      break;
	    case INTERFACE_BC_RINGLEB :
	      Local_SolnBlk[nb].W[i][j] += RinglebFlowAverageState(Local_SolnBlk[nb].W[i][j],
								   Local_SolnBlk[nb].Grid.nodeSW(i,j).X,
								   Local_SolnBlk[nb].Grid.nodeSE(i,j).X,
								   Local_SolnBlk[nb].Grid.nodeNE(i,j).X,
								   Local_SolnBlk[nb].Grid.nodeNW(i,j).X)*Local_SolnBlk[nb].Grid.lfaceS(i,j);
	      break;
	    default :
	      cout << endl << " -> " << Local_SolnBlk[nb].Grid.Cell[i][j].Xc << " " << Interface_BC_Type;
	      return 777;
	    };
	    // Update neighbour_flag.
	    neighbour_flag = ON;
	  }
	}

	// EAST face.
	if (i < Local_SolnBlk[nb].ICu+Local_SolnBlk[nb].Nghost) {
	  if (Mesh[nb].cell_status[i+1][j] == CELL_STATUS_ACTIVE &&
	      Local_SolnBlk[nb].Grid.lfaceE(i,j) > TOLER) {
	    length += Local_SolnBlk[nb].Grid.lfaceE(i,j);
	    // Determine the boundary condition type.
	    Interface_BC_Type = Interface_Union_List[Ni].Determine_Interface_BC_Type(Local_SolnBlk[nb].Grid.xfaceE(i,j));
	    // Determine the boundary velocity.
	    V = Interface_Union_List[Ni].Determine_Interface_Velocity(Local_SolnBlk[nb].Grid.xfaceE(i,j),Time);
	    // Apply the boundary condition.
	    switch(Interface_BC_Type) {
	    case INTERFACE_BC_REFLECTION :
	      Local_SolnBlk[nb].W[i][j] += Reflect(Local_SolnBlk[nb].W[i+1][j],-Local_SolnBlk[nb].Grid.nfaceE(i,j),V)*Local_SolnBlk[nb].Grid.lfaceE(i,j);
	      break;
	    case INTERFACE_BC_BURNING_SURFACE :
	      Local_SolnBlk[nb].W[i][j] += BurningSurface(Local_SolnBlk[nb].W[i+1][j],-Local_SolnBlk[nb].Grid.nfaceE(i,j))*Local_SolnBlk[nb].Grid.lfaceE(i,j);
	      break;
	    case INTERFACE_BC_RINGLEB :
	      Local_SolnBlk[nb].W[i][j] += RinglebFlowAverageState(Local_SolnBlk[nb].W[i][j],
								   Local_SolnBlk[nb].Grid.nodeSW(i,j).X,
								   Local_SolnBlk[nb].Grid.nodeSE(i,j).X,
								   Local_SolnBlk[nb].Grid.nodeNE(i,j).X,
								   Local_SolnBlk[nb].Grid.nodeNW(i,j).X)*Local_SolnBlk[nb].Grid.lfaceE(i,j);
	      break;
	    };
	    // Update neighbour_flag.
	    neighbour_flag = ON;
	  }
	}

	// WEST face.
	if (i > Local_SolnBlk[nb].ICl-Local_SolnBlk[nb].Nghost) {
	  if (Mesh[nb].cell_status[i-1][j] == CELL_STATUS_ACTIVE &&
	      Local_SolnBlk[nb].Grid.lfaceW(i,j) > TOLER) {
	    length += Local_SolnBlk[nb].Grid.lfaceW(i,j);
	    // Determine the boundary condition type.
	    Interface_BC_Type = Interface_Union_List[Ni].Determine_Interface_BC_Type(Local_SolnBlk[nb].Grid.xfaceW(i,j));
	    // Determine the boundary velocity.
	    V = Interface_Union_List[Ni].Determine_Interface_Velocity(Local_SolnBlk[nb].Grid.xfaceW(i,j),Time);
	    // Apply the boundary condition.
	    switch(Interface_BC_Type) {
	    case INTERFACE_BC_REFLECTION :
	      Local_SolnBlk[nb].W[i][j] += Reflect(Local_SolnBlk[nb].W[i-1][j],-Local_SolnBlk[nb].Grid.nfaceW(i,j),V)*Local_SolnBlk[nb].Grid.lfaceW(i,j);
	      break;
	    case INTERFACE_BC_BURNING_SURFACE :
	      Local_SolnBlk[nb].W[i][j] += BurningSurface(Local_SolnBlk[nb].W[i-1][j],-Local_SolnBlk[nb].Grid.nfaceW(i,j))*Local_SolnBlk[nb].Grid.lfaceW(i,j);
	      break;
	    case INTERFACE_BC_RINGLEB :
	      Local_SolnBlk[nb].W[i][j] += RinglebFlowAverageState(Local_SolnBlk[nb].W[i][j],
								   Local_SolnBlk[nb].Grid.nodeSW(i,j).X,
								   Local_SolnBlk[nb].Grid.nodeSE(i,j).X,
								   Local_SolnBlk[nb].Grid.nodeNE(i,j).X,
								   Local_SolnBlk[nb].Grid.nodeNW(i,j).X)*Local_SolnBlk[nb].Grid.lfaceW(i,j);
	      break;
	    };
	    // Update neighbour_flag.
	    neighbour_flag = ON;
	  }
	}

	// Apply boundary condition for any active distance-2 neighbour
	// cells if required (if any exist and there were no active
	// neighbour cells.
	if (!neighbour_flag) {

	  // NORTH face.
	  if (j <= Local_SolnBlk[nb].JCu) {
	    if (Mesh[nb].cell_status[i][j+2] == CELL_STATUS_ACTIVE &&
		Local_SolnBlk[nb].Grid.lfaceN(i,j+1) > TOLER) {
	      length += Local_SolnBlk[nb].Grid.lfaceN(i,j+1);
	      // Determine the boundary condition type.
	      Interface_BC_Type = Interface_Union_List[Ni].Determine_Interface_BC_Type(Local_SolnBlk[nb].Grid.xfaceN(i,j+1));
	      // Determine the boundary velocity.
	      V = Interface_Union_List[Ni].Determine_Interface_Velocity(Local_SolnBlk[nb].Grid.xfaceN(i,j+1),Time);
	      // Apply the boundary condition.
	      switch(Interface_BC_Type) {
	      case INTERFACE_BC_REFLECTION :
		Local_SolnBlk[nb].W[i][j] += Reflect(Local_SolnBlk[nb].W[i][j+2],-Local_SolnBlk[nb].Grid.nfaceN(i,j+1),V)*Local_SolnBlk[nb].Grid.lfaceN(i,j+1);
		break;
	      case INTERFACE_BC_BURNING_SURFACE :
		Local_SolnBlk[nb].W[i][j] += BurningSurface(Local_SolnBlk[nb].W[i][j+2],-Local_SolnBlk[nb].Grid.nfaceN(i,j+1))*Local_SolnBlk[nb].Grid.lfaceN(i,j+1);
		break;
	      case INTERFACE_BC_RINGLEB :
		Local_SolnBlk[nb].W[i][j] += RinglebFlowAverageState(Local_SolnBlk[nb].W[i][j],
								     Local_SolnBlk[nb].Grid.nodeSW(i,j).X,
								     Local_SolnBlk[nb].Grid.nodeSE(i,j).X,
								     Local_SolnBlk[nb].Grid.nodeNE(i,j).X,
								     Local_SolnBlk[nb].Grid.nodeNW(i,j).X)*Local_SolnBlk[nb].Grid.lfaceN(i,j+1);
		break;
	      };
	    }
	  }

	  // SOUTH face.
	  if (j >= Local_SolnBlk[nb].JCl) {
	    if (Mesh[nb].cell_status[i][j-2] == CELL_STATUS_ACTIVE &&
		Local_SolnBlk[nb].Grid.lfaceS(i,j-1) > TOLER) {
	      length += Local_SolnBlk[nb].Grid.lfaceS(i,j-1);
	      // Determine the boundary condition type.
	      Interface_BC_Type = Interface_Union_List[Ni].Determine_Interface_BC_Type(Local_SolnBlk[nb].Grid.xfaceS(i,j-1));
	      // Determine the boundary velocity.
	      V = Interface_Union_List[Ni].Determine_Interface_Velocity(Local_SolnBlk[nb].Grid.xfaceS(i,j-1),Time);
	      // Apply the boundary condition.
	      switch(Interface_BC_Type) {
	      case INTERFACE_BC_REFLECTION :
		Local_SolnBlk[nb].W[i][j] += Reflect(Local_SolnBlk[nb].W[i][j-2],-Local_SolnBlk[nb].Grid.nfaceS(i,j-1),V)*Local_SolnBlk[nb].Grid.lfaceS(i,j-1);
		break;
	      case INTERFACE_BC_BURNING_SURFACE :
		Local_SolnBlk[nb].W[i][j] += BurningSurface(Local_SolnBlk[nb].W[i][j-2],-Local_SolnBlk[nb].Grid.nfaceS(i,j-1))*Local_SolnBlk[nb].Grid.lfaceS(i,j-1);
		break;
	      case INTERFACE_BC_RINGLEB :
		Local_SolnBlk[nb].W[i][j] += RinglebFlowAverageState(Local_SolnBlk[nb].W[i][j],
								     Local_SolnBlk[nb].Grid.nodeSW(i,j).X,
								     Local_SolnBlk[nb].Grid.nodeSE(i,j).X,
								     Local_SolnBlk[nb].Grid.nodeNE(i,j).X,
								     Local_SolnBlk[nb].Grid.nodeNW(i,j).X)*Local_SolnBlk[nb].Grid.lfaceS(i,j-1);
		break;
	      };
	    }
	  }

	  // EAST face.
	  if (i <= Local_SolnBlk[nb].ICu) {
	    if (Mesh[nb].cell_status[i+2][j] == CELL_STATUS_ACTIVE &&
		Local_SolnBlk[nb].Grid.lfaceE(i+1,j) > TOLER) {
	      length += Local_SolnBlk[nb].Grid.lfaceE(i+1,j);
	      // Determine the boundary condition type.
	      Interface_BC_Type = Interface_Union_List[Ni].Determine_Interface_BC_Type(Local_SolnBlk[nb].Grid.xfaceE(i+1,j));
	      // Determine the boundary velocity.
	      V = Interface_Union_List[Ni].Determine_Interface_Velocity(Local_SolnBlk[nb].Grid.xfaceE(i+1,j),Time);
	      // Apply the boundary condition.
	      switch(Interface_BC_Type) {
	      case INTERFACE_BC_REFLECTION :
		Local_SolnBlk[nb].W[i][j] += Reflect(Local_SolnBlk[nb].W[i+2][j],-Local_SolnBlk[nb].Grid.nfaceE(i+1,j),V)*Local_SolnBlk[nb].Grid.lfaceE(i+1,j);
		break;
	      case INTERFACE_BC_BURNING_SURFACE :
		Local_SolnBlk[nb].W[i][j] += BurningSurface(Local_SolnBlk[nb].W[i+2][j],-Local_SolnBlk[nb].Grid.nfaceE(i+1,j))*Local_SolnBlk[nb].Grid.lfaceE(i+1,j);
		break;
	      case INTERFACE_BC_RINGLEB :
		Local_SolnBlk[nb].W[i][j] += RinglebFlowAverageState(Local_SolnBlk[nb].W[i][j],
								     Local_SolnBlk[nb].Grid.nodeSW(i,j).X,
								     Local_SolnBlk[nb].Grid.nodeSE(i,j).X,
								     Local_SolnBlk[nb].Grid.nodeNE(i,j).X,
								     Local_SolnBlk[nb].Grid.nodeNW(i,j).X)*Local_SolnBlk[nb].Grid.lfaceE(i+1,j);
		break;
	      };
	    }
	  }

	  // WEST face.
	  if (i >= Local_SolnBlk[nb].ICl) {
	    if (Mesh[nb].cell_status[i-2][j] == CELL_STATUS_ACTIVE &&
		Local_SolnBlk[nb].Grid.lfaceW(i-1,j) > TOLER) {
	      length += Local_SolnBlk[nb].Grid.lfaceW(i-1,j);
	      // Determine the boundary condition type.
	      Interface_BC_Type = Interface_Union_List[Ni].Determine_Interface_BC_Type(Local_SolnBlk[nb].Grid.xfaceW(i-1,j));
	      // Determine the boundary velocity.
	      V = Interface_Union_List[Ni].Determine_Interface_Velocity(Local_SolnBlk[nb].Grid.xfaceW(i-1,j),Time);
	      // Apply the boundary condition.
	      switch(Interface_BC_Type) {
	      case INTERFACE_BC_REFLECTION :
		Local_SolnBlk[nb].W[i][j] += Reflect(Local_SolnBlk[nb].W[i-2][j],-Local_SolnBlk[nb].Grid.nfaceW(i-1,j),V)*Local_SolnBlk[nb].Grid.lfaceW(i-1,j);
		break;
	      case INTERFACE_BC_BURNING_SURFACE :
		Local_SolnBlk[nb].W[i][j] += BurningSurface(Local_SolnBlk[nb].W[i-2][j],-Local_SolnBlk[nb].Grid.nfaceW(i-1,j))*Local_SolnBlk[nb].Grid.lfaceW(i-1,j);
		break;
	      case INTERFACE_BC_RINGLEB :
		Local_SolnBlk[nb].W[i][j] += RinglebFlowAverageState(Local_SolnBlk[nb].W[i][j],
								     Local_SolnBlk[nb].Grid.nodeSW(i,j).X,
								     Local_SolnBlk[nb].Grid.nodeSE(i,j).X,
								     Local_SolnBlk[nb].Grid.nodeNE(i,j).X,
								     Local_SolnBlk[nb].Grid.nodeNW(i,j).X)*Local_SolnBlk[nb].Grid.lfaceW(i-1,j);
		break;
	      };
	    }
	  }

	}

	// Determine the final (length-averaged) boundary state.
	if (length > ZERO) Local_SolnBlk[nb].W[i][j] = Local_SolnBlk[nb].W[i][j]/length;
	else Local_SolnBlk[nb].W[i][j] = IP->Wo;

	// Set conserved variable solution state.
	Local_SolnBlk[nb].U[i][j] = U(Local_SolnBlk[nb].W[i][j]);

      }

    }
  }
 
  // Embedded boundary boundary conditions applied successfully.
  return 0;

}

/**********************************************************************
 * EmbeddedBoundaries2D:: CFL --                                      *
 **********************************************************************/
template <> inline double EmbeddedBoundaries2D<Euler2D_cState,
					       Euler2D_pState,
					       Euler2D_Quad_Block,
					       Euler2D_Input_Parameters>::
CFL(const double &Time) {

  double dtMin = MILLION, a, v, V, dt_inv;

  // Determine the allowable time step for each solution block.
  for (int nb = 0; nb < Local_Solution_Block_List->Nblk; nb++) {
    if (Local_Solution_Block_List->Block[nb].used == ADAPTIVEBLOCK2D_USED) {
      for (int j = Local_SolnBlk[nb].JCl-Local_SolnBlk[nb].Nghost; j <= Local_SolnBlk[nb].JCu+Local_SolnBlk[nb].Nghost; j++) {
	for (int i = Local_SolnBlk[nb].ICl-Local_SolnBlk[nb].Nghost; i <= Local_SolnBlk[nb].ICu+Local_SolnBlk[nb].Nghost; i++) {
	  if (i < Local_SolnBlk[nb].ICl || i > Local_SolnBlk[nb].ICu || j < Local_SolnBlk[nb].JCl || j > Local_SolnBlk[nb].JCu) {
	    Local_SolnBlk[nb].dt[i][j] = ZERO;
	  } else if (Mesh[nb].cell_status[i][j] != CELL_STATUS_ACTIVE) {
	    Local_SolnBlk[nb].dt[i][j] = ZERO;
	  } else {

	    // Initialize the time-step.
	    dt_inv = ZERO;

	    // Determine the speed of sound in the current cell
	    // for the inviscid time-step calculation:
	    a = Local_SolnBlk[nb].W[i][j].a();

	    // Determine the NORTH face contribution if required.
	    if (Local_SolnBlk[nb].Grid.lfaceN(i,j) > TOLER) {
	      // Inviscid time-step.
	      if (Mesh[nb].cell_status[i][j+1] == CELL_STATUS_ACTIVE) V = ZERO;
	      else V = Interface_Union_List[Mesh[nb].cell_status[i][j+1]].Determine_Interface_Velocity(Local_SolnBlk[nb].Grid.xfaceN(i,j),Time)*Local_SolnBlk[nb].Grid.nfaceN(i,j);
	      v = Local_SolnBlk[nb].W[i][j].v*Local_SolnBlk[nb].Grid.nfaceN(i,j);
	      dt_inv += (a+fabs(v-V)+fabs(V))*Local_SolnBlk[nb].Grid.lfaceN(i,j);
	    }

	    // Determine the SOUTH face contribution if required.
	    if (Local_SolnBlk[nb].Grid.lfaceS(i,j) > TOLER) {
	      // Inviscid time-step.
	      if (Mesh[nb].cell_status[i][j-1] == CELL_STATUS_ACTIVE) V = ZERO;
	      else V = Interface_Union_List[Mesh[nb].cell_status[i][j-1]].Determine_Interface_Velocity(Local_SolnBlk[nb].Grid.xfaceS(i,j),Time)*Local_SolnBlk[nb].Grid.nfaceS(i,j);
	      v = Local_SolnBlk[nb].W[i][j].v*Local_SolnBlk[nb].Grid.nfaceS(i,j);
	      dt_inv += (a+fabs(v-V)+fabs(V))*Local_SolnBlk[nb].Grid.lfaceS(i,j);
	    }

	    // Determine the EAST face contribution if required.
	    if (Local_SolnBlk[nb].Grid.lfaceE(i,j) > TOLER) {
	      // Inviscid time-step.
	      if (Mesh[nb].cell_status[i+1][j] == CELL_STATUS_ACTIVE) V = ZERO;
	      else V = Interface_Union_List[Mesh[nb].cell_status[i+1][j]].Determine_Interface_Velocity(Local_SolnBlk[nb].Grid.xfaceE(i,j),Time)*Local_SolnBlk[nb].Grid.nfaceE(i,j);
	      v = Local_SolnBlk[nb].W[i][j].v*Local_SolnBlk[nb].Grid.nfaceE(i,j);
	      dt_inv += (a+fabs(v-V)+fabs(V))*Local_SolnBlk[nb].Grid.lfaceE(i,j);
	    }

	    // Determine the WEST face contribution if required.
	    if (Local_SolnBlk[nb].Grid.lfaceW(i,j) > TOLER) {
	      // Inviscid time-step.
	      if (Mesh[nb].cell_status[i-1][j] == CELL_STATUS_ACTIVE) V = ZERO;
	      else V = Interface_Union_List[Mesh[nb].cell_status[i-1][j]].Determine_Interface_Velocity(Local_SolnBlk[nb].Grid.xfaceW(i,j),Time)*Local_SolnBlk[nb].Grid.nfaceW(i,j);
	      v = Local_SolnBlk[nb].W[i][j].v*Local_SolnBlk[nb].Grid.nfaceW(i,j);
	      dt_inv += (a+fabs(v-V)+fabs(V))*Local_SolnBlk[nb].Grid.lfaceW(i,j);
	    }

	    // Determine the minimum time-step for the cell.
	    Local_SolnBlk[nb].dt[i][j] = Local_SolnBlk[nb].Grid.Cell[i][j].A/dt_inv;

	    // Determine the mimimum global time-step.
	    dtMin = min(dtMin,Local_SolnBlk[nb].dt[i][j]);

	  }
	}
      }

      for (int j = Local_SolnBlk[nb].JCl-Local_SolnBlk[nb].Nghost; j <= Local_SolnBlk[nb].JCu+Local_SolnBlk[nb].Nghost; j++) {
	for (int i = Local_SolnBlk[nb].ICl-Local_SolnBlk[nb].Nghost; i <= Local_SolnBlk[nb].ICu+Local_SolnBlk[nb].Nghost; i++) {
	  if ((i < Local_SolnBlk[nb].ICl || i > Local_SolnBlk[nb].ICu ||
	       j < Local_SolnBlk[nb].JCl || j > Local_SolnBlk[nb].JCu) ||
	      Mesh[nb].cell_status[i][j] != CELL_STATUS_ACTIVE) {
	    Local_SolnBlk[nb].dt[i][j] = dtMin;
	  }
	}
      }

    }
  }

  // Return the global time step.
  return dtMin;

}

/**********************************************************************
 * EmbeddedBoundaries2D::dUdt_Residual_Evaluation --                  *
 *                                                                    *
 * This routine evaluates the residual for the specified solution     *
 * block using a 2nd-order limited upwind finite-volume spatial       *
 * discretization scheme.  The residual is stored in dUdt[][][0].     *
 *                                                                    *
 **********************************************************************/
template <> int EmbeddedBoundaries2D<Euler2D_cState,
				     Euler2D_pState,
				     Euler2D_Quad_Block,
				     Euler2D_Input_Parameters>::
dUdt_Residual_Evaluation(const double &Time) {

  Vector2D dX;
  Euler2D_pState Wl, Wr;
  Euler2D_cState Flux;

  int Interface_BC_Type, Ni;
  Vector2D V;

  double q1, q2, k;

  // Evaluate the time rate of change of the solution (i.e., the
  // solution residuals) using a second-order limited upwind scheme
  // with a variety of flux functions.

  // Evaluate the stage solution residual for each solution block.
  for (int nb = 0; nb < Local_Solution_Block_List->Nblk; nb++) {
    if (Local_Solution_Block_List->Block[nb].used == ADAPTIVEBLOCK2D_USED) {

      // Perform the linear reconstruction within each cell of the
      // computational grid for this stage.
      Linear_Least_Squares_Reconstruction(nb);

      // Add i-direction (zeta-direction) fluxes.
      for (int j = Local_SolnBlk[nb].JCl-1; j <= Local_SolnBlk[nb].JCu+1; j++) {
	Local_SolnBlk[nb].dUdt[Local_SolnBlk[nb].ICl-1][j][0].Vacuum();

	for (int i = Local_SolnBlk[nb].ICl-1; i <= Local_SolnBlk[nb].ICu; i++) {

	  Local_SolnBlk[nb].dUdt[i+1][j][0].Vacuum();

	  if (j >= Local_SolnBlk[nb].JCl && j <= Local_SolnBlk[nb].JCu) {

 	    // Evaluate the cell interface i-direction INVISCID fluxes.
 	    if (Local_SolnBlk[nb].Grid.lfaceE(i,j) < TOLER) {

 	      // EAST face of cell (i,j) has zero length.
 	      Wl.Standard_Atmosphere(); Wr.Standard_Atmosphere(); V = Vector2D_ZERO;

 	    } else if (Mesh[nb].cell_status[i  ][j] != CELL_STATUS_ACTIVE &&
 		       Mesh[nb].cell_status[i+1][j] != CELL_STATUS_ACTIVE) {

 	      // EAST face of cell (i,j) is inactive.
 	      Wl.Standard_Atmosphere(); Wr.Standard_Atmosphere(); V = Vector2D_ZERO;

 	    } else if (Mesh[nb].cell_status[i  ][j] != CELL_STATUS_ACTIVE &&
 		       Mesh[nb].cell_status[i+1][j] == CELL_STATUS_ACTIVE) {

	      // WEST face of cell (i+1,j) corresponds to an embedded boundary.
	      dX = Local_SolnBlk[nb].Grid.xfaceW(i+1,j)-Local_SolnBlk[nb].Grid.Cell[i+1][j].Xc;
	      Wr = Local_SolnBlk[nb].W[i+1][j] + (Local_SolnBlk[nb].phi[i+1][j]^Local_SolnBlk[nb].dWdx[i+1][j])*dX.x +
		                                 (Local_SolnBlk[nb].phi[i+1][j]^Local_SolnBlk[nb].dWdy[i+1][j])*dX.y;
	      Ni = Mesh[nb].cell_status[i][j];
	      // Determine the boundary condition type.
	      Interface_BC_Type = Interface_Union_List[Ni].Determine_Interface_BC_Type(Local_SolnBlk[nb].Grid.xfaceW(i+1,j));
	      // Determine the boundary velocity.
	      V = Interface_Union_List[Ni].Determine_Interface_Velocity(Local_SolnBlk[nb].Grid.xfaceW(i+1,j),Time);
	      // Determine the left state by applying the appropriate
	      // boundary condition.
	      if (Interface_BC_Type == INTERFACE_BC_REFLECTION) {
		// WEST face of cell (i+1,j) is a REFLECTION boundary.
		Wl = Reflect(Wr,Local_SolnBlk[nb].Grid.nfaceW(i+1,j),V);
	      } else if (Interface_BC_Type == INTERFACE_BC_BURNING_SURFACE) {
		// WEST face of cell (i+1,j) is a BURNING_SURFACE boundary.
		Wl = BurningSurface(Wr,Local_SolnBlk[nb].Grid.nfaceW(i+1,j));
	      } else if (Interface_BC_Type == INTERFACE_BC_RINGLEB) {
		// WEST face of cell (i+1,j) is a RINGLEB_FLOW boundary.
		Wl = RinglebFlow(Wr,Local_SolnBlk[nb].Grid.Cell[i][j].Xc,q1,k);
		Wl = RinglebFlow(Wr,Local_SolnBlk[nb].Grid.Cell[i+1][j].Xc,q2,k);
		Wl = RinglebFlow(Wr,Local_SolnBlk[nb].Grid.xfaceW(i+1,j));
		if (q1 < IP->Isotach_Line && q2 > IP->Isotach_Line)
		  Wl = BC_Characteristic_Pressure(Wr,Wl,Local_SolnBlk[nb].Grid.nfaceW(i+1,j));
	      }

 	    } else if (Mesh[nb].cell_status[i  ][j] == CELL_STATUS_ACTIVE &&
 		       Mesh[nb].cell_status[i+1][j] != CELL_STATUS_ACTIVE) {

	      // EAST face of cell (i,j) corresponds to an embedded boundary.
	      dX = Local_SolnBlk[nb].Grid.xfaceE(i,j) - Local_SolnBlk[nb].Grid.Cell[i][j].Xc;
	      Wl = Local_SolnBlk[nb].W[i][j] + (Local_SolnBlk[nb].phi[i][j]^Local_SolnBlk[nb].dWdx[i][j])*dX.x +
	                                       (Local_SolnBlk[nb].phi[i][j]^Local_SolnBlk[nb].dWdy[i][j])*dX.y;
	      Ni = Mesh[nb].cell_status[i+1][j];
	      // Determine the boundary condition type.
	      Interface_BC_Type = Interface_Union_List[Ni].Determine_Interface_BC_Type(Local_SolnBlk[nb].Grid.xfaceE(i,j));
	      // Determine the boundary velocity.
	      V = Interface_Union_List[Ni].Determine_Interface_Velocity(Local_SolnBlk[nb].Grid.xfaceE(i,j),Time);
	      // Determine the right state by applying the appropriate
	      // boundary condition.
	      if (Interface_BC_Type == INTERFACE_BC_REFLECTION) {
		// EAST face of cell (i,j) is a REFLECTION boundary.
		Wr = Reflect(Wl,Local_SolnBlk[nb].Grid.nfaceE(i,j),V);
	      } else if (Interface_BC_Type == INTERFACE_BC_BURNING_SURFACE) {
		// EAST face of cell (i,j) is a BURNING_SURFACE boundary.
		Wr = BurningSurface(Wl,Local_SolnBlk[nb].Grid.nfaceE(i,j));
	      } else if (Interface_BC_Type == INTERFACE_BC_RINGLEB) {
		// EAST face of cell (i,j) is a RINGLEB FLOW boundary.
		Wr = RinglebFlow(Wl,Local_SolnBlk[nb].Grid.Cell[i+1][j].Xc,q1,k);
		Wr = RinglebFlow(Wl,Local_SolnBlk[nb].Grid.Cell[i][j].Xc,q2,k);
		Wr = RinglebFlow(Wl,Local_SolnBlk[nb].Grid.xfaceE(i,j));
		if (q1 < IP->Isotach_Line && q2 > IP->Isotach_Line)
		  Wr = BC_Characteristic_Pressure(Wl,Wr,Local_SolnBlk[nb].Grid.nfaceE(i,j));
	      }

 	    } else if (Mesh[nb].cell_status[i  ][j] == CELL_STATUS_ACTIVE &&
 		       Mesh[nb].cell_status[i+1][j] == CELL_STATUS_ACTIVE) {

 	      V = Vector2D_ZERO;

 	      if (i == Local_SolnBlk[nb].ICl-1 && 
 		  (Local_SolnBlk[nb].Grid.BCtypeW[j] == BC_REFLECTION ||
 		   Local_SolnBlk[nb].Grid.BCtypeW[j] == BC_BURNING_SURFACE ||
 		   Local_SolnBlk[nb].Grid.BCtypeW[j] == BC_RINGLEB_FLOW ||
 		   Local_SolnBlk[nb].Grid.BCtypeW[j] == BC_CHARACTERISTIC)) {

		dX = Local_SolnBlk[nb].Grid.xfaceW(i+1,j) - Local_SolnBlk[nb].Grid.Cell[i+1][j].Xc;
		Wr = Local_SolnBlk[nb].W[i+1][j] + (Local_SolnBlk[nb].phi[i+1][j]^Local_SolnBlk[nb].dWdx[i+1][j])*dX.x +
	                                           (Local_SolnBlk[nb].phi[i+1][j]^Local_SolnBlk[nb].dWdy[i+1][j])*dX.y;

		// WEST face of cell (i+1,j) is a normal boundary.
		if (Local_SolnBlk[nb].Grid.BCtypeW[j] == BC_REFLECTION) {
		  // WEST face of cell (i+1,j) is a REFLECTION boundary.
		  Wl = Reflect(Wr,Local_SolnBlk[nb].Grid.nfaceW(i+1,j));
		} else if (Local_SolnBlk[nb].Grid.BCtypeW[j] == BC_BURNING_SURFACE) {
		  // WEST face of cell (i+1,j) is a BURNING_SURFACE boundary.
		  Wl = BurningSurface(Wr,Local_SolnBlk[nb].Grid.nfaceW(i+1,j));
		} else if (Local_SolnBlk[nb].Grid.BCtypeW[j] == BC_RINGLEB_FLOW) {
		  // WEST face of cell (i+1,j) is a RINGLEB_FLOW boundary.
		  Wl = RinglebFlow(Wl,Local_SolnBlk[nb].Grid.xfaceW(i+1,j));
		} else {
		  // WEST face of cell (i+1,j) is a CHARACTERISTIC boundary.
		  Wl = BC_Characteristic_Pressure(Wr,Local_SolnBlk[nb].WoW[j],Local_SolnBlk[nb].Grid.nfaceW(i+1,j));
		}

	      } else if (i == Local_SolnBlk[nb].ICu &&
			 (Local_SolnBlk[nb].Grid.BCtypeE[j] == BC_REFLECTION ||
			  Local_SolnBlk[nb].Grid.BCtypeE[j] == BC_BURNING_SURFACE ||
			  Local_SolnBlk[nb].Grid.BCtypeE[j] == BC_RINGLEB_FLOW ||
			  Local_SolnBlk[nb].Grid.BCtypeE[j] == BC_CHARACTERISTIC)) {

		dX = Local_SolnBlk[nb].Grid.xfaceE(i,j) - Local_SolnBlk[nb].Grid.Cell[i][j].Xc;
		Wl = Local_SolnBlk[nb].W[i][j] + (Local_SolnBlk[nb].phi[i][j]^Local_SolnBlk[nb].dWdx[i][j])*dX.x +
	                                         (Local_SolnBlk[nb].phi[i][j]^Local_SolnBlk[nb].dWdy[i][j])*dX.y;

		// EAST face of cell (i,j) is a normal boundary.
		if (Local_SolnBlk[nb].Grid.BCtypeE[j] == BC_REFLECTION) {
		  // EAST face of cell (i,j) is a REFLECTION boundary.
		  Wr = Reflect(Wl,Local_SolnBlk[nb].Grid.nfaceE(i,j));
		} else if (Local_SolnBlk[nb].Grid.BCtypeE[j] == BC_BURNING_SURFACE) {
		  // EAST face of cell (i,j) is a BURNING_SURFACE boundary.
		  Wr = BurningSurface(Wl,Local_SolnBlk[nb].Grid.nfaceE(i,j));
		} else if (Local_SolnBlk[nb].Grid.BCtypeE[j] == BC_RINGLEB_FLOW) {
		  // EAST face of cell (i,j) is a RINGLEB_FLOW boundary.
		  Wr = RinglebFlow(Wr,Local_SolnBlk[nb].Grid.xfaceE(i,j));
		} else {
		  // EAST face of cell (i,j) is a CHARACTERISTIC boundary.
		  Wr = BC_Characteristic_Pressure(Wl,Local_SolnBlk[nb].WoE[j],Local_SolnBlk[nb].Grid.nfaceE(i,j));
		}

 	      } else {

		// EAST face is either a normal cell or possibly a FIXED, 
		// NONE or EXTRAPOLATION boundary.
		dX = Local_SolnBlk[nb].Grid.xfaceE(i  ,j) - Local_SolnBlk[nb].Grid.Cell[i  ][j].Xc;
		Wl = Local_SolnBlk[nb].W[i  ][j] + (Local_SolnBlk[nb].phi[i  ][j]^Local_SolnBlk[nb].dWdx[i  ][j])*dX.x +
	                                           (Local_SolnBlk[nb].phi[i  ][j]^Local_SolnBlk[nb].dWdy[i  ][j])*dX.y;
		dX = Local_SolnBlk[nb].Grid.xfaceW(i+1,j) - Local_SolnBlk[nb].Grid.Cell[i+1][j].Xc;
		Wr = Local_SolnBlk[nb].W[i+1][j] + (Local_SolnBlk[nb].phi[i+1][j]^Local_SolnBlk[nb].dWdx[i+1][j])*dX.x +
	                                           (Local_SolnBlk[nb].phi[i+1][j]^Local_SolnBlk[nb].dWdy[i+1][j])*dX.y;

 	      }

 	    }

	    // Determine EAST face INVISCID flux.
	    switch(IP->i_Flux_Function) {
	    case FLUX_FUNCTION_GODUNOV :
	      Flux = FluxGodunov_n(Wl,Wr,Local_SolnBlk[nb].Grid.nfaceE(i,j));
	      break;
	    case FLUX_FUNCTION_ROE :
	      Flux = FluxRoe_n(Wl,Wr,Local_SolnBlk[nb].Grid.nfaceE(i,j));
	      break;
	    case FLUX_FUNCTION_RUSANOV :
	      Flux = FluxRusanov_n(Wl,Wr,Local_SolnBlk[nb].Grid.nfaceE(i,j));
	      break;
	    case FLUX_FUNCTION_HLLE :
	      Flux = FluxHLLE_n(Wl,Wr,Local_SolnBlk[nb].Grid.nfaceE(i,j));
	      break;
	    case FLUX_FUNCTION_LINDE :
	      Flux = FluxLinde_n(Wl,Wr,Local_SolnBlk[nb].Grid.nfaceE(i,j));
	      break;
	    case FLUX_FUNCTION_HLLC :
	      Flux = FluxHLLC_n(Wl,Wr,Local_SolnBlk[nb].Grid.nfaceE(i,j));
	      break;
	    case FLUX_FUNCTION_VANLEER :
	      Flux = FluxVanLeer_n(Wl,Wr,Local_SolnBlk[nb].Grid.nfaceE(i,j));
	      break;
	    case FLUX_FUNCTION_AUSM :
	      Flux = FluxAUSM_n(Wl,Wr,Local_SolnBlk[nb].Grid.nfaceE(i,j));
	      break;
	    case FLUX_FUNCTION_AUSMplus :
	      Flux = FluxAUSMplus_n(Wl,Wr,Local_SolnBlk[nb].Grid.nfaceE(i,j));
	      break;
	    case FLUX_FUNCTION_GODUNOV_MB :
	      Flux = FluxGodunov_MB_n(Wl,Wr,V,Local_SolnBlk[nb].Grid.nfaceE(i,j));
	      break;
	    case FLUX_FUNCTION_ROE_MB :
	      Flux = FluxRoe_MB_n(Wl,Wr,V,Local_SolnBlk[nb].Grid.nfaceE(i,j));
	      break;
	    case FLUX_FUNCTION_HLLE_MB :
	      Flux = FluxHLLE_MB_n(Wl,Wr,V,Local_SolnBlk[nb].Grid.nfaceE(i,j));
	      break;
	    case FLUX_FUNCTION_VANLEER_MB :
	      Flux = FluxVanLeer_MB_n(Wl,Wr,V,Local_SolnBlk[nb].Grid.nfaceE(i,j));
	      break;
	    default:
	      Flux = FluxRoe_n(Wl,Wr,Local_SolnBlk[nb].Grid.nfaceE(i,j));
	      break;
	    };

	    // Evaluate cell-averaged solution changes.
	    Local_SolnBlk[nb].dUdt[i  ][j][0] -= Flux*Local_SolnBlk[nb].Grid.lfaceE(i,j)/
	                                         Local_SolnBlk[nb].Grid.Cell[i][j].A;
	    Local_SolnBlk[nb].dUdt[i+1][j][0] += Flux*Local_SolnBlk[nb].Grid.lfaceW(i+1,j)/
	                                         Local_SolnBlk[nb].Grid.Cell[i+1][j].A;

	    // Include all axisymmetric source term if required.
	    if (Local_SolnBlk[nb].Axisymmetric)
	      Local_SolnBlk[nb].dUdt[i][j][0] += Local_SolnBlk[nb].W[i][j].S(Local_SolnBlk[nb].Grid.Cell[i][j].Xc);

	    // Include area change source term.
	    if (Interface_Union_List.Ni)
	      Local_SolnBlk[nb].dUdt[i][j][0] -= Local_SolnBlk[nb].U[i][j]*dAdt(nb,i,j,Time);

	    // Save west and east face boundary flux.
	    if (i == Local_SolnBlk[nb].ICl-1) {
	      Local_SolnBlk[nb].FluxW[j] = -Flux*Local_SolnBlk[nb].Grid.lfaceW(i+1,j);
	    } else if (i == Local_SolnBlk[nb].ICu) {
	      Local_SolnBlk[nb].FluxE[j] =  Flux*Local_SolnBlk[nb].Grid.lfaceE(i,j);
	    }

 	  }
	}

	if (j > Local_SolnBlk[nb].JCl-1 && j < Local_SolnBlk[nb].JCu+1) {
	  Local_SolnBlk[nb].dUdt[Local_SolnBlk[nb].ICl-1][j][0].Vacuum();
	  Local_SolnBlk[nb].dUdt[Local_SolnBlk[nb].ICu+1][j][0].Vacuum();
	}

      }

      // Add j-direction (eta-direction) fluxes.
      for (int i = Local_SolnBlk[nb].ICl; i <= Local_SolnBlk[nb].ICu; i++) {
 	for (int j = Local_SolnBlk[nb].JCl-1; j <= Local_SolnBlk[nb].JCu; j++) {

 	  // Evaluate the cell interface j-direction INVISCID fluxes.
 	  if (Local_SolnBlk[nb].Grid.lfaceN(i,j) < TOLER) {
 	    // NORTH face of cell (i,j) has zero length.
 	    Wl.Standard_Atmosphere(); Wr.Standard_Atmosphere(); V = Vector2D_ZERO;

 	  } else if (Mesh[nb].cell_status[i][j  ] != CELL_STATUS_ACTIVE &&
 		     Mesh[nb].cell_status[i][j+1] != CELL_STATUS_ACTIVE) {
 	    // NORTH face of cell (i,j) is inactive.
 	    Wl.Standard_Atmosphere(); Wr.Standard_Atmosphere(); V = Vector2D_ZERO;

 	  } else if (Mesh[nb].cell_status[i][j  ] != CELL_STATUS_ACTIVE &&
 		     Mesh[nb].cell_status[i][j+1] == CELL_STATUS_ACTIVE) {
	    // SOUTH face of cell (i,j+1) corresponds to an embedded boundary.
	    dX = Local_SolnBlk[nb].Grid.xfaceS(i,j+1) - Local_SolnBlk[nb].Grid.Cell[i][j+1].Xc;
	    Wr = Local_SolnBlk[nb].W[i][j+1] + (Local_SolnBlk[nb].phi[i][j+1]^Local_SolnBlk[nb].dWdx[i][j+1])*dX.x +
	                                       (Local_SolnBlk[nb].phi[i][j+1]^Local_SolnBlk[nb].dWdy[i][j+1])*dX.y;
	    Ni = Mesh[nb].cell_status[i][j];
	    // Determine the boundary condition type.
	    Interface_BC_Type = Interface_Union_List[Ni].Determine_Interface_BC_Type(Local_SolnBlk[nb].Grid.xfaceS(i,j+1));
	    // Determine the boundary velocity.
	    V = Interface_Union_List[Ni].Determine_Interface_Velocity(Local_SolnBlk[nb].Grid.xfaceS(i,j+1),Time);
	    // Determine the left state by applying the appropriate
	    // boundary condition.
	    if (Interface_BC_Type == INTERFACE_BC_REFLECTION) {
	      // SOUTH face of cell (i,j+1) is a REFLECTION boundary.
	      Wl = Reflect(Wr,Local_SolnBlk[nb].Grid.nfaceS(i,j+1),V);
	    } else if (Interface_BC_Type == INTERFACE_BC_BURNING_SURFACE) {
	      // SOUTH face of cell (i,j+1) is a BURNING_SURFACE boundary.
	      Wl = BurningSurface(Wr,Local_SolnBlk[nb].Grid.nfaceS(i,j+1));
	    } else if (Interface_BC_Type == INTERFACE_BC_RINGLEB) {
	      // SOUTH face of cell (i,j+1) is a RINGLEB_FLOW boundary.
	      Wl = RinglebFlow(Wr,Local_SolnBlk[nb].Grid.Cell[i][j].Xc,q1,k);
	      Wl = RinglebFlow(Wr,Local_SolnBlk[nb].Grid.Cell[i][j+1].Xc,q2,k);
	      Wl = RinglebFlow(Wr,Local_SolnBlk[nb].Grid.xfaceS(i,j+1));
	      if (q1 < IP->Isotach_Line && q2 > IP->Isotach_Line)
		Wl = BC_Characteristic_Pressure(Wr,Wl,Local_SolnBlk[nb].Grid.nfaceS(i,j+1));
	    }

 	  } else if (Mesh[nb].cell_status[i][j  ] == CELL_STATUS_ACTIVE &&
 		     Mesh[nb].cell_status[i][j+1] != CELL_STATUS_ACTIVE) {
	    // NORTH face of cell (i,j) corresponds to an embedded boundary.
	    dX = Local_SolnBlk[nb].Grid.xfaceN(i,j) - Local_SolnBlk[nb].Grid.Cell[i][j].Xc;
	    Wl = Local_SolnBlk[nb].W[i][j] + (Local_SolnBlk[nb].phi[i][j]^Local_SolnBlk[nb].dWdx[i][j])*dX.x +
	                                     (Local_SolnBlk[nb].phi[i][j]^Local_SolnBlk[nb].dWdy[i][j])*dX.y;
	    Ni = Mesh[nb].cell_status[i][j+1];
	    // Determine the boundary condition type.
	    Interface_BC_Type = Interface_Union_List[Ni].Determine_Interface_BC_Type(Local_SolnBlk[nb].Grid.xfaceN(i,j));
	    // Determine the boundary velocity.
	    V = Interface_Union_List[Ni].Determine_Interface_Velocity(Local_SolnBlk[nb].Grid.xfaceN(i,j),Time);
	    // Determine the right state by applying the appropriate
	    // boundary condition.
	    if (Interface_BC_Type == INTERFACE_BC_REFLECTION) {
	      // NORTH face of cell (i,j) is a REFLECTION boundary.
	      Wr = Reflect(Wl,Local_SolnBlk[nb].Grid.nfaceN(i,j),V);
	    } else if (Interface_BC_Type == INTERFACE_BC_BURNING_SURFACE) {
	      // NORTH face of cell (i,j) is a BURNING_SURFACE boundary.
	      Wr = BurningSurface(Wl,Local_SolnBlk[nb].Grid.nfaceN(i,j));
	    } else if (Interface_BC_Type == INTERFACE_BC_RINGLEB) {
	      // NORTH face of cell (i,j) is a RINGLEB_FLOW boundary.
	      Wr = RinglebFlow(Wl,Local_SolnBlk[nb].Grid.Cell[i][j+1].Xc,q1,k);
	      Wr = RinglebFlow(Wl,Local_SolnBlk[nb].Grid.Cell[i][j].Xc,q2,k);
	      Wr = RinglebFlow(Wl,Local_SolnBlk[nb].Grid.xfaceN(i,j));
	      if (q1 < IP->Isotach_Line && q2 > IP->Isotach_Line)
		Wr = BC_Characteristic_Pressure(Wl,Wr,Local_SolnBlk[nb].Grid.nfaceN(i,j));
	    }

 	  } else if (Mesh[nb].cell_status[i][j  ] == CELL_STATUS_ACTIVE &&
 		     Mesh[nb].cell_status[i][j+1] == CELL_STATUS_ACTIVE) {

 	    V = Vector2D_ZERO;

 	    if (j == Local_SolnBlk[nb].JCl-1 && 
 		(Local_SolnBlk[nb].Grid.BCtypeS[i] == BC_REFLECTION ||
 		 Local_SolnBlk[nb].Grid.BCtypeS[i] == BC_BURNING_SURFACE ||
 		 Local_SolnBlk[nb].Grid.BCtypeS[i] == BC_RINGLEB_FLOW ||
 		 Local_SolnBlk[nb].Grid.BCtypeS[i] == BC_CHARACTERISTIC)) {

	      dX = Local_SolnBlk[nb].Grid.xfaceS(i,j+1) - Local_SolnBlk[nb].Grid.Cell[i][j+1].Xc;
	      Wr = Local_SolnBlk[nb].W[i][j+1] + (Local_SolnBlk[nb].phi[i][j+1]^Local_SolnBlk[nb].dWdx[i][j+1])*dX.x +
   	                                         (Local_SolnBlk[nb].phi[i][j+1]^Local_SolnBlk[nb].dWdy[i][j+1])*dX.y;
	      // SOUTH face of cell (i,j+1) is a normal boundary.
	      if (Local_SolnBlk[nb].Grid.BCtypeS[i] == BC_REFLECTION) {
		// SOUTH face of cell (i,j+1) is a REFLECTION boundary.
		Wl = Reflect(Wr,Local_SolnBlk[nb].Grid.nfaceS(i,j+1));
	      } else if (Local_SolnBlk[nb].Grid.BCtypeS[i] == BC_BURNING_SURFACE) {
		// SOUTH face of cell (i,j+1) is a BURNING_SURFACE boundary.
		Wl = BurningSurface(Wr,Local_SolnBlk[nb].Grid.nfaceS(i,j+1));
	      } else if (Local_SolnBlk[nb].Grid.BCtypeS[i] == BC_RINGLEB_FLOW) {
		// SOUTH face of cell (i,j+1) is a RINGLEB_FLOW boundary.
		Wl = RinglebFlow(Wl,Local_SolnBlk[nb].Grid.xfaceS(i,j+1));
		Wl = BC_Characteristic_Pressure(Wr,Wl,Local_SolnBlk[nb].Grid.nfaceS(i,j+1));
	      } else if (Local_SolnBlk[nb].Grid.BCtypeS[i] == BC_CHARACTERISTIC) {
		// SOUTH face of cell (i,j+1) is a CHARACTERISTIC boundary.
		Wl = BC_Characteristic_Pressure(Wr,Local_SolnBlk[nb].WoS[i],Local_SolnBlk[nb].Grid.nfaceS(i,j+1));
	      }

	    } else if (j == Local_SolnBlk[nb].JCu && 
		       (Local_SolnBlk[nb].Grid.BCtypeN[i] == BC_REFLECTION ||
			Local_SolnBlk[nb].Grid.BCtypeN[i] == BC_BURNING_SURFACE ||
			Local_SolnBlk[nb].Grid.BCtypeN[i] == BC_RINGLEB_FLOW ||
			Local_SolnBlk[nb].Grid.BCtypeN[i] == BC_CHARACTERISTIC)) {

	      dX = Local_SolnBlk[nb].Grid.xfaceN(i,j) - Local_SolnBlk[nb].Grid.Cell[i][j].Xc;
	      Wl = Local_SolnBlk[nb].W[i][j] + (Local_SolnBlk[nb].phi[i][j]^Local_SolnBlk[nb].dWdx[i][j])*dX.x +
   	                                       (Local_SolnBlk[nb].phi[i][j]^Local_SolnBlk[nb].dWdy[i][j])*dX.y;
	      // NORTH face of cell (i,j) is a normal boundary.
	      if (Local_SolnBlk[nb].Grid.BCtypeN[i] == BC_REFLECTION) {
		// NORTH face of cell (i,j) is a REFLECTION boundary.
		Wr = Reflect(Wl,Local_SolnBlk[nb].Grid.nfaceN(i,j));
	      } else if (Local_SolnBlk[nb].Grid.BCtypeN[i] == BC_BURNING_SURFACE) {
		// NORTH face of cell (i,j) is a BURNING_SURFACE boundary.
		Wr = BurningSurface(Wl,Local_SolnBlk[nb].Grid.nfaceN(i,j));
	      } else if (Local_SolnBlk[nb].Grid.BCtypeN[i] == BC_RINGLEB_FLOW) {
		// NORTH face of cell (i,j) is a RINGLEB_FLOW boundary.
		Wr = RinglebFlow(Wr,Local_SolnBlk[nb].Grid.xfaceN(i,j));
		Wr = BC_Characteristic_Pressure(Wl,Wr,Local_SolnBlk[nb].Grid.nfaceN(i,j));
	      } else {
		// NORTH face of cell (i,j) is a CHARACTERISTIC boundary.
		Wr = BC_Characteristic_Pressure(Wl,Local_SolnBlk[nb].WoN[i],Local_SolnBlk[nb].Grid.nfaceN(i,j));
	      }

	    } else {

	      // NORTH face is either a normal cell or possibly a FIXED, 
	      // NONE or EXTRAPOLATION boundary.
	      dX = Local_SolnBlk[nb].Grid.xfaceN(i,j  ) - Local_SolnBlk[nb].Grid.Cell[i][j  ].Xc;
	      Wl = Local_SolnBlk[nb].W[i][j  ] + (Local_SolnBlk[nb].phi[i][j  ]^Local_SolnBlk[nb].dWdx[i][j  ])*dX.x +
                                       (Local_SolnBlk[nb].phi[i][j  ]^Local_SolnBlk[nb].dWdy[i][j  ])*dX.y;
	      dX = Local_SolnBlk[nb].Grid.xfaceS(i,j+1) - Local_SolnBlk[nb].Grid.Cell[i][j+1].Xc;
	      Wr = Local_SolnBlk[nb].W[i][j+1] + (Local_SolnBlk[nb].phi[i][j+1]^Local_SolnBlk[nb].dWdx[i][j+1])*dX.x +
                                       (Local_SolnBlk[nb].phi[i][j+1]^Local_SolnBlk[nb].dWdy[i][j+1])*dX.y;

	    }

	  }

	  // Determine NORTH face inviscid flux.
	  switch(IP->i_Flux_Function) {
	  case FLUX_FUNCTION_GODUNOV :
	    Flux = FluxGodunov_n(Wl,Wr,Local_SolnBlk[nb].Grid.nfaceN(i,j));
	    break;
	  case FLUX_FUNCTION_ROE :
	    Flux = FluxRoe_n(Wl,Wr,Local_SolnBlk[nb].Grid.nfaceN(i,j));
	    break;
	  case FLUX_FUNCTION_RUSANOV :
	    Flux = FluxRusanov_n(Wl,Wr,Local_SolnBlk[nb].Grid.nfaceN(i,j));
	    break;
	  case FLUX_FUNCTION_HLLE :
	    Flux = FluxHLLE_n(Wl,Wr,Local_SolnBlk[nb].Grid.nfaceN(i,j));
	    break;
	  case FLUX_FUNCTION_LINDE :
	    Flux = FluxLinde_n(Wl,Wr,Local_SolnBlk[nb].Grid.nfaceN(i,j));
	    break;
	  case FLUX_FUNCTION_HLLC :
	    Flux = FluxHLLC_n(Wl,Wr,Local_SolnBlk[nb].Grid.nfaceN(i,j));
	    break;
	  case FLUX_FUNCTION_VANLEER :
	    Flux = FluxVanLeer_n(Wl,Wr,Local_SolnBlk[nb].Grid.nfaceN(i,j));
	    break;
	  case FLUX_FUNCTION_AUSM :
	    Flux = FluxAUSM_n(Wl,Wr,Local_SolnBlk[nb].Grid.nfaceN(i,j));
	    break;
	  case FLUX_FUNCTION_AUSMplus :
	    Flux = FluxAUSMplus_n(Wl,Wr,Local_SolnBlk[nb].Grid.nfaceN(i,j));
	    break;
	  case FLUX_FUNCTION_GODUNOV_MB :
	    Flux = FluxGodunov_MB_n(Wl,Wr,V,Local_SolnBlk[nb].Grid.nfaceN(i,j));
	    break;
	  case FLUX_FUNCTION_ROE_MB :
	    Flux = FluxRoe_MB_n(Wl,Wr,V,Local_SolnBlk[nb].Grid.nfaceN(i,j));
	    break;
	  case FLUX_FUNCTION_HLLE_MB :
	    Flux = FluxHLLE_MB_n(Wl,Wr,V,Local_SolnBlk[nb].Grid.nfaceN(i,j));
	    break;
	  case FLUX_FUNCTION_VANLEER_MB :
	    Flux = FluxVanLeer_MB_n(Wl,Wr,V,Local_SolnBlk[nb].Grid.nfaceN(i,j));
	    break;
	  default:
	    Flux = FluxRoe_n(Wl,Wr,Local_SolnBlk[nb].Grid.nfaceN(i,j));
	    break;
	  };

	  // Evaluate cell-averaged solution changes.
	  Local_SolnBlk[nb].dUdt[i][j  ][0] -= Flux*Local_SolnBlk[nb].Grid.lfaceN(i,j)/
	                                       Local_SolnBlk[nb].Grid.Cell[i][j].A;
	  Local_SolnBlk[nb].dUdt[i][j+1][0] += Flux*Local_SolnBlk[nb].Grid.lfaceS(i,j+1)/
	                                       Local_SolnBlk[nb].Grid.Cell[i][j+1].A;

	  // Save south and north face boundary flux.
	  if (j == Local_SolnBlk[nb].JCl-1) {
	    Local_SolnBlk[nb].FluxS[i] = -Flux*Local_SolnBlk[nb].Grid.lfaceS(i,j+1);
	  } else if (j == Local_SolnBlk[nb].JCu) {
	    Local_SolnBlk[nb].FluxN[i] = Flux*Local_SolnBlk[nb].Grid.lfaceN(i,j);
	  }

	}

	Local_SolnBlk[nb].dUdt[i][Local_SolnBlk[nb].JCl-1][0].Vacuum();
	Local_SolnBlk[nb].dUdt[i][Local_SolnBlk[nb].JCu+1][0].Vacuum();

      }

    }
  }

  // Residuals for each quadrilateral multi-block solution block
  // successfully calculated.
  return 0;

}

/**********************************************************************
 * EmbeddedBoundaries2D::dUdt_Multistage_Explicit --                  *
 *                                                                    *
 * This routine evaluates the stage solution residual for a 1D array  *
 * of 2D quadrilateral multi-block solution blocks.  A variety of     *
 * multistage explicit time integration and upwind finite-volume      *
 * spatial discretization procedures can be used depending on the     *
 * specified input values.                                            *
 *                                                                    *
 **********************************************************************/
template <> int EmbeddedBoundaries2D<Euler2D_cState,
				     Euler2D_pState,
				     Euler2D_Quad_Block,
				     Euler2D_Input_Parameters>::
dUdt_Multistage_Explicit(const int &i_stage,
			 const double &Time) {

  int error_flag, k_residual;
  double omega;
  Vector2D dX;
  Euler2D_pState Wl, Wr;
  Euler2D_cState Flux;

  int Interface_BC_Type, Ni;
  Vector2D V;

  double q1, q2, k;

  // Evaluate the solution residual for stage i_stage of n_stage scheme.

  // Evaluate the time step fraction and residual storage location for
  // the stage.
  switch(IP->i_Time_Integration) {
  case TIME_STEPPING_EXPLICIT_EULER :
    omega = Runge_Kutta(i_stage,IP->N_Stage);
    k_residual = 0;
    break;
  case TIME_STEPPING_EXPLICIT_PREDICTOR_CORRECTOR :
    omega = Runge_Kutta(i_stage,IP->N_Stage);
    k_residual = 0;
    break;
  case TIME_STEPPING_EXPLICIT_RUNGE_KUTTA :
    omega = Runge_Kutta(i_stage,IP->N_Stage);
    k_residual = 0;
    if (IP->N_Stage == 4) {
      if (i_stage == 4) k_residual = 0;
      else k_residual = i_stage - 1;
    }
    break;
  case TIME_STEPPING_MULTISTAGE_OPTIMAL_SMOOTHING :
    omega = MultiStage_Optimally_Smoothing(i_stage,IP->N_Stage,IP->i_Limiter);
    k_residual = 0;
    break;
  default:
    omega = Runge_Kutta(i_stage,IP->N_Stage);
    k_residual = 0;
    break;
  };

  // Evaluate the time rate of change of the solution (i.e., the
  // solution residuals) using a second-order limited upwind scheme
  // with a variety of flux functions.

  // Evaluate the stage solution residual for each solution block.
  for (int nb = 0; nb < Local_Solution_Block_List->Nblk; nb++) {
    if (Local_Solution_Block_List->Block[nb].used == ADAPTIVEBLOCK2D_USED) {

      // Perform the linear reconstruction within each cell of the
      // computational grid for this stage.
      Linear_Least_Squares_Reconstruction(nb);

      // Add i-direction (zeta-direction) fluxes.
      for (int j = Local_SolnBlk[nb].JCl-1; j <= Local_SolnBlk[nb].JCu+1; j++) {
	if (i_stage == 1) {
	  Local_SolnBlk[nb].Uo[Local_SolnBlk[nb].ICl-1][j] = Local_SolnBlk[nb].U[Local_SolnBlk[nb].ICl-1][j];
	  Local_SolnBlk[nb].dUdt[Local_SolnBlk[nb].ICl-1][j][k_residual].Vacuum();
	} else {
	  Local_SolnBlk[nb].dUdt[Local_SolnBlk[nb].ICl-1][j][k_residual].Vacuum();
	}

	for (int i = Local_SolnBlk[nb].ICl-1; i <= Local_SolnBlk[nb].ICu; i++) {

	  if (i_stage == 1) {
	    Local_SolnBlk[nb].Uo[i+1][j] = Local_SolnBlk[nb].U[i+1][j];
	    Local_SolnBlk[nb].dUdt[i+1][j][k_residual].Vacuum();
	  } else if (j > Local_SolnBlk[nb].JCl-1 && j < Local_SolnBlk[nb].JCu+1) {
	    switch(IP->i_Time_Integration) {
	    case TIME_STEPPING_EXPLICIT_PREDICTOR_CORRECTOR :
	      //Local_SolnBlk[nb].dUdt[i+1][j][k_residual] = Local_SolnBlk[nb].dUdt[i+1][j][k_residual];
	      break;
	    case TIME_STEPPING_EXPLICIT_RUNGE_KUTTA :
	      if (IP->N_Stage == 2) {
		//Local_SolnBlk[nb].dUdt[i+1][j][k_residual] = Local_SolnBlk[nb].dUdt[i+1][j][k_residual];
	      } else if (IP->N_Stage == 4 && i_stage == 4) {
		Local_SolnBlk[nb].dUdt[i+1][j][k_residual] = Local_SolnBlk[nb].dUdt[i+1][j][0] + 
                                                             TWO*Local_SolnBlk[nb].dUdt[i+1][j][1] +
                                                             TWO*Local_SolnBlk[nb].dUdt[i+1][j][2];
	      } else {
		Local_SolnBlk[nb].dUdt[i+1][j][k_residual].Vacuum();
	      }
	      break;
	    case TIME_STEPPING_MULTISTAGE_OPTIMAL_SMOOTHING :
	      Local_SolnBlk[nb].dUdt[i+1][j][k_residual].Vacuum();
	      break;
	    default:
	      Local_SolnBlk[nb].dUdt[i+1][j][k_residual].Vacuum();
	      break;
	    };
	  }

	  if (j >= Local_SolnBlk[nb].JCl && j <= Local_SolnBlk[nb].JCu) {

 	    // Evaluate the cell interface i-direction INVISCID fluxes.
 	    if (Local_SolnBlk[nb].Grid.lfaceE(i,j) < TOLER) {
 	      // EAST face of cell (i,j) has zero length.
 	      Wl.Standard_Atmosphere(); Wr.Standard_Atmosphere(); V = Vector2D_ZERO;

 	    } else if (Mesh[nb].cell_status[i  ][j] != CELL_STATUS_ACTIVE &&
 		       Mesh[nb].cell_status[i+1][j] != CELL_STATUS_ACTIVE) {
 	      // EAST face of cell (i,j) is inactive.
 	      Wl.Standard_Atmosphere(); Wr.Standard_Atmosphere(); V = Vector2D_ZERO;

 	    } else if (Mesh[nb].cell_status[i  ][j] != CELL_STATUS_ACTIVE &&
 		       Mesh[nb].cell_status[i+1][j] == CELL_STATUS_ACTIVE) {
	      // WEST face of cell (i+1,j) corresponds to an embedded boundary.
	      dX = Local_SolnBlk[nb].Grid.xfaceW(i+1,j)-Local_SolnBlk[nb].Grid.Cell[i+1][j].Xc;
	      Wr = Local_SolnBlk[nb].W[i+1][j] + (Local_SolnBlk[nb].phi[i+1][j]^Local_SolnBlk[nb].dWdx[i+1][j])*dX.x +
		                                 (Local_SolnBlk[nb].phi[i+1][j]^Local_SolnBlk[nb].dWdy[i+1][j])*dX.y;
	      Ni = Mesh[nb].cell_status[i][j];
	      // Determine the boundary condition type.
	      Interface_BC_Type = Interface_Union_List[Ni].Determine_Interface_BC_Type(Local_SolnBlk[nb].Grid.xfaceW(i+1,j));
	      // Determine the boundary velocity.
	      V = Interface_Union_List[Ni].Determine_Interface_Velocity(Local_SolnBlk[nb].Grid.xfaceW(i+1,j),Time);
	      // Determine the left state by applying the appropriate
	      // boundary condition.
	      if (Interface_BC_Type == INTERFACE_BC_REFLECTION) {
		// WEST face of cell (i+1,j) is a REFLECTION boundary.
		Wl = Reflect(Wr,Local_SolnBlk[nb].Grid.nfaceW(i+1,j),V);
	      } else if (Interface_BC_Type == INTERFACE_BC_BURNING_SURFACE) {
		// WEST face of cell (i+1,j) is a BURNING_SURFACE boundary.
		Wl = BurningSurface(Wr,Local_SolnBlk[nb].Grid.nfaceW(i+1,j));
	      } else if (Interface_BC_Type == INTERFACE_BC_RINGLEB) {
		// WEST face of cell (i+1,j) is a RINGLEB_FLOW boundary.
		Wl = RinglebFlow(Wr,Local_SolnBlk[nb].Grid.Cell[i][j].Xc,q1,k);
		Wl = RinglebFlow(Wr,Local_SolnBlk[nb].Grid.Cell[i+1][j].Xc,q2,k);
		Wl = RinglebFlow(Wr,Local_SolnBlk[nb].Grid.xfaceW(i+1,j));
		if (q1 < IP->Isotach_Line && q2 > IP->Isotach_Line)
		  Wl = BC_Characteristic_Pressure(Wr,Wl,Local_SolnBlk[nb].Grid.nfaceW(i+1,j));
	      }
 	    } else if (Mesh[nb].cell_status[i  ][j] == CELL_STATUS_ACTIVE &&
 		       Mesh[nb].cell_status[i+1][j] != CELL_STATUS_ACTIVE) {
	      // EAST face of cell (i,j) corresponds to an embedded boundary.
	      dX = Local_SolnBlk[nb].Grid.xfaceE(i,j) - Local_SolnBlk[nb].Grid.Cell[i][j].Xc;
	      Wl = Local_SolnBlk[nb].W[i][j] + (Local_SolnBlk[nb].phi[i][j]^Local_SolnBlk[nb].dWdx[i][j])*dX.x +
	                                       (Local_SolnBlk[nb].phi[i][j]^Local_SolnBlk[nb].dWdy[i][j])*dX.y;
	      Ni = Mesh[nb].cell_status[i+1][j];
	      // Determine the boundary condition type.
	      Interface_BC_Type = Interface_Union_List[Ni].Determine_Interface_BC_Type(Local_SolnBlk[nb].Grid.xfaceE(i,j));
	      // Determine the boundary velocity.
	      V = Interface_Union_List[Ni].Determine_Interface_Velocity(Local_SolnBlk[nb].Grid.xfaceE(i,j),Time);
	      // Determine the right state by applying the appropriate
	      // boundary condition.
	      if (Interface_BC_Type == INTERFACE_BC_REFLECTION) {
		// EAST face of cell (i,j) is a REFLECTION boundary.
		Wr = Reflect(Wl,Local_SolnBlk[nb].Grid.nfaceE(i,j),V);
	      } else if (Interface_BC_Type == INTERFACE_BC_BURNING_SURFACE) {
		// EAST face of cell (i,j) is a BURNING_SURFACE boundary.
		Wr = BurningSurface(Wl,Local_SolnBlk[nb].Grid.nfaceE(i,j));
	      } else if (Interface_BC_Type == INTERFACE_BC_RINGLEB) {
		// EAST face of cell (i,j) is a RINGLEB FLOW boundary.
		Wr = RinglebFlow(Wr,Local_SolnBlk[nb].Grid.Cell[i+1][j].Xc,q1,k);
		Wr = RinglebFlow(Wr,Local_SolnBlk[nb].Grid.Cell[i][j].Xc,q2,k);
		Wr = RinglebFlow(Wr,Local_SolnBlk[nb].Grid.xfaceE(i,j));
		if (q1 < IP->Isotach_Line && q2 > IP->Isotach_Line)
		  Wr = BC_Characteristic_Pressure(Wl,Wr,Local_SolnBlk[nb].Grid.nfaceE(i,j));
	      }

 	    } else if (Mesh[nb].cell_status[i  ][j] == CELL_STATUS_ACTIVE &&
 		       Mesh[nb].cell_status[i+1][j] == CELL_STATUS_ACTIVE) {

 	      V = Vector2D_ZERO;

 	      if (i == Local_SolnBlk[nb].ICl-1 && 
 		  (Local_SolnBlk[nb].Grid.BCtypeW[j] == BC_REFLECTION ||
 		   Local_SolnBlk[nb].Grid.BCtypeW[j] == BC_BURNING_SURFACE ||
 		   Local_SolnBlk[nb].Grid.BCtypeW[j] == BC_RINGLEB_FLOW ||
 		   Local_SolnBlk[nb].Grid.BCtypeW[j] == BC_CHARACTERISTIC)) {

		dX = Local_SolnBlk[nb].Grid.xfaceW(i+1,j) - Local_SolnBlk[nb].Grid.Cell[i+1][j].Xc;
		Wr = Local_SolnBlk[nb].W[i+1][j] + (Local_SolnBlk[nb].phi[i+1][j]^Local_SolnBlk[nb].dWdx[i+1][j])*dX.x +
	                                           (Local_SolnBlk[nb].phi[i+1][j]^Local_SolnBlk[nb].dWdy[i+1][j])*dX.y;

		// WEST face of cell (i+1,j) is a normal boundary.
		if (Local_SolnBlk[nb].Grid.BCtypeW[j] == BC_REFLECTION) {
		  // WEST face of cell (i+1,j) is a REFLECTION boundary.
		  Wl = Reflect(Wr,Local_SolnBlk[nb].Grid.nfaceW(i+1,j));
		} else if (Local_SolnBlk[nb].Grid.BCtypeW[j] == BC_BURNING_SURFACE) {
		  // WEST face of cell (i+1,j) is a BURNING_SURFACE boundary.
		  Wl = BurningSurface(Wr,Local_SolnBlk[nb].Grid.nfaceW(i+1,j));
		} else if (Local_SolnBlk[nb].Grid.BCtypeW[j] == BC_RINGLEB_FLOW) {
		  // WEST face of cell (i+1,j) is a RINGLEB_FLOW boundary.
		  Wl = RinglebFlow(Wl,Local_SolnBlk[nb].Grid.xfaceW(i+1,j));
		} else {
		  // WEST face of cell (i+1,j) is a CHARACTERISTIC boundary.
		  Wl = BC_Characteristic_Pressure(Wr,Local_SolnBlk[nb].WoW[j],Local_SolnBlk[nb].Grid.nfaceW(i+1,j));
		}

	      } else if (i == Local_SolnBlk[nb].ICu &&
			 (Local_SolnBlk[nb].Grid.BCtypeE[j] == BC_REFLECTION ||
			  Local_SolnBlk[nb].Grid.BCtypeE[j] == BC_BURNING_SURFACE ||
			  Local_SolnBlk[nb].Grid.BCtypeE[j] == BC_RINGLEB_FLOW ||
			  Local_SolnBlk[nb].Grid.BCtypeE[j] == BC_CHARACTERISTIC)) {

		dX = Local_SolnBlk[nb].Grid.xfaceE(i,j) - Local_SolnBlk[nb].Grid.Cell[i][j].Xc;
		Wl = Local_SolnBlk[nb].W[i][j] + (Local_SolnBlk[nb].phi[i][j]^Local_SolnBlk[nb].dWdx[i][j])*dX.x +
	                                         (Local_SolnBlk[nb].phi[i][j]^Local_SolnBlk[nb].dWdy[i][j])*dX.y;

		// EAST face of cell (i,j) is a normal boundary.
		if (Local_SolnBlk[nb].Grid.BCtypeE[j] == BC_REFLECTION) {
		  // EAST face of cell (i,j) is a REFLECTION boundary.
		  Wr = Reflect(Wl,Local_SolnBlk[nb].Grid.nfaceE(i,j));
		} else if (Local_SolnBlk[nb].Grid.BCtypeE[j] == BC_BURNING_SURFACE) {
		  // EAST face of cell (i,j) is a BURNING_SURFACE boundary.
		  Wr = BurningSurface(Wl,Local_SolnBlk[nb].Grid.nfaceE(i,j));
		} else if (Local_SolnBlk[nb].Grid.BCtypeE[j] == BC_RINGLEB_FLOW) {
		  // EAST face of cell (i,j) is a RINGLEB_FLOW boundary.
		  Wr = RinglebFlow(Wr,Local_SolnBlk[nb].Grid.xfaceE(i,j));
		} else {
		  // EAST face of cell (i,j) is a CHARACTERISTIC boundary.
		  Wr = BC_Characteristic_Pressure(Wl,Local_SolnBlk[nb].WoE[j],Local_SolnBlk[nb].Grid.nfaceE(i,j));
		}

 	      } else {

		// EAST face is either a normal cell or possibly a FIXED, 
		// NONE or EXTRAPOLATION boundary.
		dX = Local_SolnBlk[nb].Grid.xfaceE(i  ,j) - Local_SolnBlk[nb].Grid.Cell[i  ][j].Xc;
		Wl = Local_SolnBlk[nb].W[i  ][j] + (Local_SolnBlk[nb].phi[i  ][j]^Local_SolnBlk[nb].dWdx[i  ][j])*dX.x +
	                                           (Local_SolnBlk[nb].phi[i  ][j]^Local_SolnBlk[nb].dWdy[i  ][j])*dX.y;
		dX = Local_SolnBlk[nb].Grid.xfaceW(i+1,j) - Local_SolnBlk[nb].Grid.Cell[i+1][j].Xc;
		Wr = Local_SolnBlk[nb].W[i+1][j] + (Local_SolnBlk[nb].phi[i+1][j]^Local_SolnBlk[nb].dWdx[i+1][j])*dX.x +
	                                           (Local_SolnBlk[nb].phi[i+1][j]^Local_SolnBlk[nb].dWdy[i+1][j])*dX.y;

 	      }

 	    }

	    // Determine EAST face INVISCID flux.
	    switch(IP->i_Flux_Function) {
	    case FLUX_FUNCTION_GODUNOV :
	      Flux = FluxGodunov_n(Wl,Wr,Local_SolnBlk[nb].Grid.nfaceE(i,j));
	      break;
	    case FLUX_FUNCTION_ROE :
	      Flux = FluxRoe_n(Wl,Wr,Local_SolnBlk[nb].Grid.nfaceE(i,j));
	      break;
	    case FLUX_FUNCTION_RUSANOV :
	      Flux = FluxRusanov_n(Wl,Wr,Local_SolnBlk[nb].Grid.nfaceE(i,j));
	      break;
	    case FLUX_FUNCTION_HLLE :
	      Flux = FluxHLLE_n(Wl,Wr,Local_SolnBlk[nb].Grid.nfaceE(i,j));
	      break;
	    case FLUX_FUNCTION_LINDE :
	      Flux = FluxLinde_n(Wl,Wr,Local_SolnBlk[nb].Grid.nfaceE(i,j));
	      break;
	    case FLUX_FUNCTION_HLLC :
	      Flux = FluxHLLC_n(Wl,Wr,Local_SolnBlk[nb].Grid.nfaceE(i,j));
	      break;
	    case FLUX_FUNCTION_VANLEER :
	      Flux = FluxVanLeer_n(Wl,Wr,Local_SolnBlk[nb].Grid.nfaceE(i,j));
	      break;
	    case FLUX_FUNCTION_AUSM :
	      Flux = FluxAUSM_n(Wl,Wr,Local_SolnBlk[nb].Grid.nfaceE(i,j));
	      break;
	    case FLUX_FUNCTION_AUSMplus :
	      Flux = FluxAUSMplus_n(Wl,Wr,Local_SolnBlk[nb].Grid.nfaceE(i,j));
	      break;
	    case FLUX_FUNCTION_GODUNOV_MB :
	      Flux = FluxGodunov_MB_n(Wl,Wr,V,Local_SolnBlk[nb].Grid.nfaceE(i,j));
	      break;
	    case FLUX_FUNCTION_ROE_MB :
	      Flux = FluxRoe_MB_n(Wl,Wr,V,Local_SolnBlk[nb].Grid.nfaceE(i,j));
	      break;
	    case FLUX_FUNCTION_HLLE_MB :
	      Flux = FluxHLLE_MB_n(Wl,Wr,V,Local_SolnBlk[nb].Grid.nfaceE(i,j));
	      break;
	    case FLUX_FUNCTION_VANLEER_MB :
	      Flux = FluxVanLeer_MB_n(Wl,Wr,V,Local_SolnBlk[nb].Grid.nfaceE(i,j));
	      break;
	    default:
	      Flux = FluxRoe_n(Wl,Wr,Local_SolnBlk[nb].Grid.nfaceE(i,j));
	      break;
	    };

	    // Evaluate cell-averaged solution changes.
	    Local_SolnBlk[nb].dUdt[i  ][j][k_residual] -= (IP->CFL_Number*Local_SolnBlk[nb].dt[i][j])*
                                                          Flux*Local_SolnBlk[nb].Grid.lfaceE(i,j)/
	                                                  Local_SolnBlk[nb].Grid.Cell[i][j].A;
	    Local_SolnBlk[nb].dUdt[i+1][j][k_residual] += (IP->CFL_Number*Local_SolnBlk[nb].dt[i+1][j])*
	                                                  Flux*Local_SolnBlk[nb].Grid.lfaceW(i+1,j)/
	                                                  Local_SolnBlk[nb].Grid.Cell[i+1][j].A;

	    // Include all axisymmetric source term if required.
	    if (Local_SolnBlk[nb].Axisymmetric)
	      Local_SolnBlk[nb].dUdt[i][j][k_residual] += (IP->CFL_Number*Local_SolnBlk[nb].dt[i][j])*
		                                          Local_SolnBlk[nb].W[i][j].S(Local_SolnBlk[nb].Grid.Cell[i][j].Xc);

 	    // Include area change source term.
 	    if (Interface_Union_List.Ni)
 	      Local_SolnBlk[nb].dUdt[i][j][k_residual] -= (IP->CFL_Number*Local_SolnBlk[nb].dt[i][j])*
 		                                          Local_SolnBlk[nb].U[i][j]*dAdt(nb,i,j,Time);

	    // Save west and east face boundary flux.
	    if (i == Local_SolnBlk[nb].ICl-1) {
	      Local_SolnBlk[nb].FluxW[j] = -Flux*Local_SolnBlk[nb].Grid.lfaceW(i+1,j);
	    } else if (i == Local_SolnBlk[nb].ICu) {
	      Local_SolnBlk[nb].FluxE[j] =  Flux*Local_SolnBlk[nb].Grid.lfaceE(i,j);
	    }

 	  }
	}

	if (j > Local_SolnBlk[nb].JCl-1 && j < Local_SolnBlk[nb].JCu+1) {
	  Local_SolnBlk[nb].dUdt[Local_SolnBlk[nb].ICl-1][j][k_residual].Vacuum();
	  Local_SolnBlk[nb].dUdt[Local_SolnBlk[nb].ICu+1][j][k_residual].Vacuum();
	}

      }

      // Add j-direction (eta-direction) fluxes.
      for (int i = Local_SolnBlk[nb].ICl; i <= Local_SolnBlk[nb].ICu; i++) {
 	for (int j = Local_SolnBlk[nb].JCl-1; j <= Local_SolnBlk[nb].JCu; j++) {

 	  // Evaluate the cell interface j-direction INVISCID fluxes.
 	  if (Local_SolnBlk[nb].Grid.lfaceN(i,j) < TOLER) {
 	    // NORTH face of cell (i,j) has zero length.
 	    Wl.Standard_Atmosphere(); Wr.Standard_Atmosphere(); V = Vector2D_ZERO;

 	  } else if (Mesh[nb].cell_status[i][j  ] != CELL_STATUS_ACTIVE &&
 		     Mesh[nb].cell_status[i][j+1] != CELL_STATUS_ACTIVE) {
 	    // NORTH face of cell (i,j) is inactive.
 	    Wl.Standard_Atmosphere(); Wr.Standard_Atmosphere(); V = Vector2D_ZERO;

 	  } else if (Mesh[nb].cell_status[i][j  ] != CELL_STATUS_ACTIVE &&
 		     Mesh[nb].cell_status[i][j+1] == CELL_STATUS_ACTIVE) {
	    // SOUTH face of cell (i,j+1) corresponds to an embedded boundary.
	    dX = Local_SolnBlk[nb].Grid.xfaceS(i,j+1) - Local_SolnBlk[nb].Grid.Cell[i][j+1].Xc;
	    Wr = Local_SolnBlk[nb].W[i][j+1] + (Local_SolnBlk[nb].phi[i][j+1]^Local_SolnBlk[nb].dWdx[i][j+1])*dX.x +
	                                       (Local_SolnBlk[nb].phi[i][j+1]^Local_SolnBlk[nb].dWdy[i][j+1])*dX.y;
	    Ni = Mesh[nb].cell_status[i][j];
	    // Determine the boundary condition type.
	    Interface_BC_Type = Interface_Union_List[Ni].Determine_Interface_BC_Type(Local_SolnBlk[nb].Grid.xfaceS(i,j+1));
	    // Determine the boundary velocity.
	    V = Interface_Union_List[Ni].Determine_Interface_Velocity(Local_SolnBlk[nb].Grid.xfaceS(i,j+1),Time);
	    // Determine the left state by applying the appropriate
	    // boundary condition.
	    if (Interface_BC_Type == INTERFACE_BC_REFLECTION) {
	      // SOUTH face of cell (i,j+1) is a REFLECTION boundary.
	      Wl = Reflect(Wr,Local_SolnBlk[nb].Grid.nfaceS(i,j+1),V);
	    } else if (Interface_BC_Type == INTERFACE_BC_BURNING_SURFACE) {
	      // SOUTH face of cell (i,j+1) is a BURNING_SURFACE boundary.
	      Wl = BurningSurface(Wr,Local_SolnBlk[nb].Grid.nfaceS(i,j+1));
	    } else if (Interface_BC_Type == INTERFACE_BC_RINGLEB) {
	      // SOUTH face of cell (i,j+1) is a RINGLEB_FLOW boundary.
	      Wl = RinglebFlow(Wl,Local_SolnBlk[nb].Grid.Cell[i][j].Xc,q1,k);
	      Wl = RinglebFlow(Wl,Local_SolnBlk[nb].Grid.Cell[i][j+1].Xc,q2,k);
	      Wl = RinglebFlow(Wl,Local_SolnBlk[nb].Grid.xfaceS(i,j+1));
	      if (q1 < IP->Isotach_Line && q2 > IP->Isotach_Line)
		Wl = BC_Characteristic_Pressure(Wr,Wl,Local_SolnBlk[nb].Grid.nfaceS(i,j+1));
	    }

 	  } else if (Mesh[nb].cell_status[i][j  ] == CELL_STATUS_ACTIVE &&
 		     Mesh[nb].cell_status[i][j+1] != CELL_STATUS_ACTIVE) {
	    // NORTH face of cell (i,j) corresponds to an embedded boundary.
	    dX = Local_SolnBlk[nb].Grid.xfaceN(i,j) - Local_SolnBlk[nb].Grid.Cell[i][j].Xc;
	    Wl = Local_SolnBlk[nb].W[i][j] + (Local_SolnBlk[nb].phi[i][j]^Local_SolnBlk[nb].dWdx[i][j])*dX.x +
	                                     (Local_SolnBlk[nb].phi[i][j]^Local_SolnBlk[nb].dWdy[i][j])*dX.y;
	    Ni = Mesh[nb].cell_status[i][j+1];
	    // Determine the boundary condition type.
	    Interface_BC_Type = Interface_Union_List[Ni].Determine_Interface_BC_Type(Local_SolnBlk[nb].Grid.xfaceN(i,j));
	    // Determine the boundary velocity.
	    V = Interface_Union_List[Ni].Determine_Interface_Velocity(Local_SolnBlk[nb].Grid.xfaceN(i,j),Time);
	    // Determine the right state by applying the appropriate
	    // boundary condition.
	    if (Interface_BC_Type == INTERFACE_BC_REFLECTION) {
	      // NORTH face of cell (i,j) is a REFLECTION boundary.
	      Wr = Reflect(Wl,Local_SolnBlk[nb].Grid.nfaceN(i,j),V);
	    } else if (Interface_BC_Type == INTERFACE_BC_BURNING_SURFACE) {
	      // NORTH face of cell (i,j) is a BURNING_SURFACE boundary.
	      Wr = BurningSurface(Wl,Local_SolnBlk[nb].Grid.nfaceN(i,j));
	    } else if (Interface_BC_Type == INTERFACE_BC_RINGLEB) {
	      // NORTH face of cell (i,j) is a RINGLEB_FLOW boundary.
	      Wr = RinglebFlow(Wr,Local_SolnBlk[nb].Grid.Cell[i][j+1].Xc,q1,k);
	      Wr = RinglebFlow(Wr,Local_SolnBlk[nb].Grid.Cell[i][j].Xc,q2,k);
	      Wr = RinglebFlow(Wr,Local_SolnBlk[nb].Grid.xfaceN(i,j));
	      if (q1 < IP->Isotach_Line && q2 > IP->Isotach_Line)
		Wr = BC_Characteristic_Pressure(Wl,Wr,Local_SolnBlk[nb].Grid.nfaceN(i,j));
	    }

 	  } else if (Mesh[nb].cell_status[i][j  ] == CELL_STATUS_ACTIVE &&
 		     Mesh[nb].cell_status[i][j+1] == CELL_STATUS_ACTIVE) {

 	    V = Vector2D_ZERO;

 	    if (j == Local_SolnBlk[nb].JCl-1 && 
 		(Local_SolnBlk[nb].Grid.BCtypeS[i] == BC_REFLECTION ||
 		 Local_SolnBlk[nb].Grid.BCtypeS[i] == BC_BURNING_SURFACE ||
 		 Local_SolnBlk[nb].Grid.BCtypeS[i] == BC_RINGLEB_FLOW ||
 		 Local_SolnBlk[nb].Grid.BCtypeS[i] == BC_CHARACTERISTIC)) {

	      dX = Local_SolnBlk[nb].Grid.xfaceS(i,j+1) - Local_SolnBlk[nb].Grid.Cell[i][j+1].Xc;
	      Wr = Local_SolnBlk[nb].W[i][j+1] + (Local_SolnBlk[nb].phi[i][j+1]^Local_SolnBlk[nb].dWdx[i][j+1])*dX.x +
   	                                         (Local_SolnBlk[nb].phi[i][j+1]^Local_SolnBlk[nb].dWdy[i][j+1])*dX.y;

	      // SOUTH face of cell (i,j+1) is a normal boundary.
	      if (Local_SolnBlk[nb].Grid.BCtypeS[i] == BC_REFLECTION) {
		// SOUTH face of cell (i,j+1) is a REFLECTION boundary.
		Wl = Reflect(Wr,Local_SolnBlk[nb].Grid.nfaceS(i,j+1));
	      } else if (Local_SolnBlk[nb].Grid.BCtypeS[i] == BC_BURNING_SURFACE) {
		// SOUTH face of cell (i,j+1) is a BURNING_SURFACE boundary.
		Wl = BurningSurface(Wr,Local_SolnBlk[nb].Grid.nfaceS(i,j+1));
	      } else if (Local_SolnBlk[nb].Grid.BCtypeS[i] == BC_RINGLEB_FLOW) {
		// SOUTH face of cell (i,j+1) is a RINGLEB_FLOW boundary.
		Wl = RinglebFlow(Wl,Local_SolnBlk[nb].Grid.xfaceS(i,j+1));
		Wl = BC_Characteristic_Pressure(Wr,Wl,Local_SolnBlk[nb].Grid.nfaceS(i,j+1));
	      } else if (Local_SolnBlk[nb].Grid.BCtypeS[i] == BC_CHARACTERISTIC) {
		// SOUTH face of cell (i,j+1) is a CHARACTERISTIC boundary.
		Wl = BC_Characteristic_Pressure(Wr,Local_SolnBlk[nb].WoS[i],Local_SolnBlk[nb].Grid.nfaceS(i,j+1));
	      }

	    } else if (j == Local_SolnBlk[nb].JCu && 
		       (Local_SolnBlk[nb].Grid.BCtypeN[i] == BC_REFLECTION ||
			Local_SolnBlk[nb].Grid.BCtypeN[i] == BC_BURNING_SURFACE ||
			Local_SolnBlk[nb].Grid.BCtypeN[i] == BC_RINGLEB_FLOW ||
			Local_SolnBlk[nb].Grid.BCtypeN[i] == BC_CHARACTERISTIC)) {

	      dX = Local_SolnBlk[nb].Grid.xfaceN(i,j) - Local_SolnBlk[nb].Grid.Cell[i][j].Xc;
	      Wl = Local_SolnBlk[nb].W[i][j] + (Local_SolnBlk[nb].phi[i][j]^Local_SolnBlk[nb].dWdx[i][j])*dX.x +
   	                                       (Local_SolnBlk[nb].phi[i][j]^Local_SolnBlk[nb].dWdy[i][j])*dX.y;

	      // NORTH face of cell (i,j) is a normal boundary.
	      if (Local_SolnBlk[nb].Grid.BCtypeN[i] == BC_REFLECTION) {
		// NORTH face of cell (i,j) is a REFLECTION boundary.
		Wr = Reflect(Wl,Local_SolnBlk[nb].Grid.nfaceN(i,j));
	      } else if (Local_SolnBlk[nb].Grid.BCtypeN[i] == BC_BURNING_SURFACE) {
		// NORTH face of cell (i,j) is a BURNING_SURFACE boundary.
		Wr = BurningSurface(Wl,Local_SolnBlk[nb].Grid.nfaceN(i,j));
	      } else if (Local_SolnBlk[nb].Grid.BCtypeN[i] == BC_RINGLEB_FLOW) {
		// NORTH face of cell (i,j) is a RINGLEB_FLOW boundary.
		Wr = RinglebFlow(Wr,Local_SolnBlk[nb].Grid.xfaceN(i,j));
		Wr = BC_Characteristic_Pressure(Wl,Wr,Local_SolnBlk[nb].Grid.nfaceN(i,j));
	      } else {
		// NORTH face of cell (i,j) is a CHARACTERISTIC boundary.
		Wr = BC_Characteristic_Pressure(Wl,Local_SolnBlk[nb].WoN[i],Local_SolnBlk[nb].Grid.nfaceN(i,j));
	      }

	    } else {

	      // NORTH face is either a normal cell or possibly a FIXED, 
	      // NONE or EXTRAPOLATION boundary.
	      dX = Local_SolnBlk[nb].Grid.xfaceN(i,j  ) - Local_SolnBlk[nb].Grid.Cell[i][j  ].Xc;
	      Wl = Local_SolnBlk[nb].W[i][j  ] + (Local_SolnBlk[nb].phi[i][j  ]^Local_SolnBlk[nb].dWdx[i][j  ])*dX.x +
                                                 (Local_SolnBlk[nb].phi[i][j  ]^Local_SolnBlk[nb].dWdy[i][j  ])*dX.y;
	      dX = Local_SolnBlk[nb].Grid.xfaceS(i,j+1) - Local_SolnBlk[nb].Grid.Cell[i][j+1].Xc;
	      Wr = Local_SolnBlk[nb].W[i][j+1] + (Local_SolnBlk[nb].phi[i][j+1]^Local_SolnBlk[nb].dWdx[i][j+1])*dX.x +
                                                 (Local_SolnBlk[nb].phi[i][j+1]^Local_SolnBlk[nb].dWdy[i][j+1])*dX.y;

	    }

	  }

	  // Determine NORTH face inviscid flux.
	  switch(IP->i_Flux_Function) {
	  case FLUX_FUNCTION_GODUNOV :
	    Flux = FluxGodunov_n(Wl,Wr,Local_SolnBlk[nb].Grid.nfaceN(i,j));
	    break;
	  case FLUX_FUNCTION_ROE :
	    Flux = FluxRoe_n(Wl,Wr,Local_SolnBlk[nb].Grid.nfaceN(i,j));
	    break;
	  case FLUX_FUNCTION_RUSANOV :
	    Flux = FluxRusanov_n(Wl,Wr,Local_SolnBlk[nb].Grid.nfaceN(i,j));
	    break;
	  case FLUX_FUNCTION_HLLE :
	    Flux = FluxHLLE_n(Wl,Wr,Local_SolnBlk[nb].Grid.nfaceN(i,j));
	    break;
	  case FLUX_FUNCTION_LINDE :
	    Flux = FluxLinde_n(Wl,Wr,Local_SolnBlk[nb].Grid.nfaceN(i,j));
	    break;
	  case FLUX_FUNCTION_HLLC :
	    Flux = FluxHLLC_n(Wl,Wr,Local_SolnBlk[nb].Grid.nfaceN(i,j));
	    break;
	  case FLUX_FUNCTION_VANLEER :
	    Flux = FluxVanLeer_n(Wl,Wr,Local_SolnBlk[nb].Grid.nfaceN(i,j));
	    break;
	  case FLUX_FUNCTION_AUSM :
	    Flux = FluxAUSM_n(Wl,Wr,Local_SolnBlk[nb].Grid.nfaceN(i,j));
	    break;
	  case FLUX_FUNCTION_AUSMplus :
	    Flux = FluxAUSMplus_n(Wl,Wr,Local_SolnBlk[nb].Grid.nfaceN(i,j));
	    break;
	  case FLUX_FUNCTION_GODUNOV_MB :
	    Flux = FluxGodunov_MB_n(Wl,Wr,V,Local_SolnBlk[nb].Grid.nfaceN(i,j));
	    break;
	  case FLUX_FUNCTION_ROE_MB :
	    Flux = FluxRoe_MB_n(Wl,Wr,V,Local_SolnBlk[nb].Grid.nfaceN(i,j));
	    break;
	  case FLUX_FUNCTION_HLLE_MB :
	    Flux = FluxHLLE_MB_n(Wl,Wr,V,Local_SolnBlk[nb].Grid.nfaceN(i,j));
	    break;
	  case FLUX_FUNCTION_VANLEER_MB :
	    Flux = FluxVanLeer_MB_n(Wl,Wr,V,Local_SolnBlk[nb].Grid.nfaceN(i,j));
	    break;
	  default:
	    Flux = FluxRoe_n(Wl,Wr,Local_SolnBlk[nb].Grid.nfaceN(i,j));
	    break;
	  };

	  // Evaluate cell-averaged solution changes.
	  Local_SolnBlk[nb].dUdt[i][j  ][k_residual] -= (IP->CFL_Number*Local_SolnBlk[nb].dt[i][j  ])*
                                                        Flux*Local_SolnBlk[nb].Grid.lfaceN(i,j)/
	                                                Local_SolnBlk[nb].Grid.Cell[i][j].A;
	  Local_SolnBlk[nb].dUdt[i][j+1][k_residual] += (IP->CFL_Number*Local_SolnBlk[nb].dt[i][j+1])*
                                                        Flux*Local_SolnBlk[nb].Grid.lfaceS(i,j+1)/
	                                                Local_SolnBlk[nb].Grid.Cell[i][j+1].A;

	  // Save south and north face boundary flux.
	  if (j == Local_SolnBlk[nb].JCl-1) {
	    Local_SolnBlk[nb].FluxS[i] = -Flux*Local_SolnBlk[nb].Grid.lfaceS(i,j+1);
	  } else if (j == Local_SolnBlk[nb].JCu) {
	    Local_SolnBlk[nb].FluxN[i] = Flux*Local_SolnBlk[nb].Grid.lfaceN(i,j);
	  }

	}

	Local_SolnBlk[nb].dUdt[i][Local_SolnBlk[nb].JCl-1][k_residual].Vacuum();
	Local_SolnBlk[nb].dUdt[i][Local_SolnBlk[nb].JCu+1][k_residual].Vacuum();

      }

    }
  }

  // Residuals for each quadrilateral multi-block solution block
  // successfully calculated.
  return 0;

}

#endif // _EMBEDDEDBOUNDARIES2D_EULER2D_INCLUDED
