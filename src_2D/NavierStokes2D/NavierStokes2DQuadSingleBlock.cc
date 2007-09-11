/**********************************************************************
 * NavierStokes2DQuadSingleBlock.cc                                   *
 *                                                                    *
 * Single-block versions of subroutines for 2D Navier-Stokes multi-   *
 * block quadrilateral mesh solution classes.                         *
 *                                                                    *
 **********************************************************************/

// Include 2D NavierStokes quadrilateral mesh solution header file.

#ifndef _NAVIERSTOKES2D_QUAD_INCLUDED
#include "NavierStokes2DQuad.h"
#endif // _NAVIERSTOKES2D_QUAD_INCLUDED

/**********************************************************************
 * NavierStokes2D_Quad_Block -- Single Block External Subroutines.    *
 **********************************************************************/

/**********************************************************************
 * Routine: Broadcast_Solution_Block                                  *
 *                                                                    *
 * Broadcast quadrilateral solution block to all processors involved  *
 * in the calculation from the primary processor using the MPI        *
 * broadcast routine.                                                 *
 *                                                                    *
 **********************************************************************/
void Broadcast_Solution_Block(NavierStokes2D_Quad_Block &SolnBlk) {

#ifdef _MPI_VERSION
  int ni, nj, ng, nr, block_allocated, buffer_size;
  double *buffer;

  // Broadcast the number of cells in each direction.
  if (CFFC_Primary_MPI_Processor()) {
    ni = SolnBlk.NCi;
    nj = SolnBlk.NCj;
    ng = SolnBlk.Nghost;
    nr = SolnBlk.residual_variable;
    if (SolnBlk.U != NULL) {
      block_allocated = 1;
    } else {
      block_allocated = 0;
    }
  }

  MPI::COMM_WORLD.Bcast(&ni,1,MPI::INT,0);
  MPI::COMM_WORLD.Bcast(&nj,1,MPI::INT,0);
  MPI::COMM_WORLD.Bcast(&ng,1,MPI::INT,0);
  MPI::COMM_WORLD.Bcast(&nr,1,MPI::INT,0);
  MPI::COMM_WORLD.Bcast(&block_allocated,1,MPI::INT,0);

  // On non-primary MPI processors, allocate (re-allocate) memory for 
  // the quadrilateral solution block as necessary.
  if (!CFFC_Primary_MPI_Processor()) {
    if (SolnBlk.NCi != ni || SolnBlk.NCj != nj || SolnBlk.Nghost != ng) { 
      if (SolnBlk.U != NULL) SolnBlk.deallocate();
      if (block_allocated) SolnBlk.allocate(ni-2*ng,nj-2*ng,ng);
    }
    // Set the block static variables if they were not previously assigned.
    if (SolnBlk.residual_variable != nr) SolnBlk.residual_variable = nr;
  }

  // Broadcast the viscous flow indicator.
  MPI::COMM_WORLD.Bcast(&(SolnBlk.Flow_Type),1,MPI::INT,0);

  // Broadcast the axisymmetric/planar flow indicator.
  MPI::COMM_WORLD.Bcast(&(SolnBlk.Axisymmetric),1,MPI::INT,0);

  // Broadcast the compressibility effect correction indicator.
  MPI::COMM_WORLD.Bcast(&(SolnBlk.Compressibility_Effect),1,MPI::INT,0);

  // Broadcast the transition model correction indicator.
  MPI::COMM_WORLD.Bcast(&(SolnBlk.Transition_Model),1,MPI::INT,0);

  // Broadcast the variable Prandtl number indicator.
  MPI::COMM_WORLD.Bcast(&(SolnBlk.Variable_Prandtl),1,MPI::INT,0);

  // Broadcast the wall velocity.
  MPI::COMM_WORLD.Bcast(&(SolnBlk.Vwall.x),1,MPI::DOUBLE,0);
  MPI::COMM_WORLD.Bcast(&(SolnBlk.Vwall.y),1,MPI::DOUBLE,0);

  // Broadcast the wall temperature.
  MPI::COMM_WORLD.Bcast(&(SolnBlk.Twall),1,MPI::DOUBLE,0);

  // Broadcast the grid.
  Broadcast_Quad_Block(SolnBlk.Grid);

  // Broadcast the solution state variables.
  if (block_allocated) {
    ni = (SolnBlk.ICu+SolnBlk.Nghost) - (SolnBlk.ICl-SolnBlk.Nghost) + 1;
    nj = (SolnBlk.JCu+SolnBlk.Nghost) - (SolnBlk.JCl-SolnBlk.Nghost) + 1;
    buffer = new double[NUM_VAR_NAVIERSTOKES2D*ni*nj];

    if (CFFC_Primary_MPI_Processor()) {
      buffer_size = 0;
      for (int j = SolnBlk.JCl-SolnBlk.Nghost; j <= SolnBlk.JCu+SolnBlk.Nghost; j++) {
	for (int i = SolnBlk.ICl-SolnBlk.Nghost; i <= SolnBlk.ICu+SolnBlk.Nghost; i++) {
	  for (int k = 0; k < NUM_VAR_NAVIERSTOKES2D; k++) {
	    buffer[buffer_size] = SolnBlk.U[i][j][k+1];
	    buffer_size++;
	  }
	}
      }
    }

    buffer_size = NUM_VAR_NAVIERSTOKES2D*ni*nj;
    MPI::COMM_WORLD.Bcast(buffer,buffer_size,MPI::DOUBLE,0);

    if (!CFFC_Primary_MPI_Processor()) {
      buffer_size = 0;
      for (int j = SolnBlk.JCl-SolnBlk.Nghost; j <= SolnBlk.JCu+SolnBlk.Nghost; j++) {
	for (int i = SolnBlk.ICl-SolnBlk.Nghost; i <= SolnBlk.ICu+SolnBlk.Nghost; i++) {
	  for (int k = 0; k < NUM_VAR_NAVIERSTOKES2D; k++) {
	    SolnBlk.U[i][j][k+1] = buffer[buffer_size];
	    buffer_size++;
	  }
	  SolnBlk.W[i][j] = W(SolnBlk.U[i][j]);
	}
      }
    }

    delete []buffer; buffer = NULL;

    ni = 1;
    nj = (SolnBlk.JCu+SolnBlk.Nghost) - (SolnBlk.JCl-SolnBlk.Nghost) + 1;
    buffer = new double[2*NUM_VAR_NAVIERSTOKES2D*ni*nj];

    if (CFFC_Primary_MPI_Processor()) {
      buffer_size = 0;
      for (int j = SolnBlk.JCl-SolnBlk.Nghost; j <= SolnBlk.JCu+SolnBlk.Nghost; j++) {
	for (int k = 0; k < NUM_VAR_NAVIERSTOKES2D; k++) {
	  buffer[buffer_size] = SolnBlk.WoW[j][k+1];
	  buffer_size++;
	}
	for (int k = 0; k < NUM_VAR_NAVIERSTOKES2D; k++) {
	  buffer[buffer_size] = SolnBlk.WoE[j][k+1];
	  buffer_size++;
	}
      }
    }

    buffer_size = 2*NUM_VAR_NAVIERSTOKES2D*ni*nj;
    MPI::COMM_WORLD.Bcast(buffer,buffer_size,MPI::DOUBLE,0);

    if (!CFFC_Primary_MPI_Processor()) {
      buffer_size = 0;
      for (int j = SolnBlk.JCl-SolnBlk.Nghost; j <= SolnBlk.JCu+SolnBlk.Nghost; j++) {
	for (int k = 0; k < NUM_VAR_NAVIERSTOKES2D; k++) {
	  SolnBlk.WoW[j][k+1] = buffer[buffer_size];
	  buffer_size++;
	}
	for (int k = 0; k < NUM_VAR_NAVIERSTOKES2D; k++) {
	  SolnBlk.WoE[j][k+1] = buffer[buffer_size];
	  buffer_size++;
	}
      }
    }

    delete []buffer; buffer = NULL;

    ni = (SolnBlk.ICu+SolnBlk.Nghost) - (SolnBlk.ICl-SolnBlk.Nghost) + 1;
    nj = 1;
    buffer = new double[2*NUM_VAR_NAVIERSTOKES2D*ni*nj];

    if (CFFC_Primary_MPI_Processor()) {
      buffer_size = 0;
      for (int i = SolnBlk.ICl-SolnBlk.Nghost; i <= SolnBlk.ICu+SolnBlk.Nghost; i++) {
	for (int k = 0; k < NUM_VAR_NAVIERSTOKES2D; k++) {
	  buffer[buffer_size] = SolnBlk.WoS[i][k+1];
	  buffer_size++;
	}
	for (int k = 0; k < NUM_VAR_NAVIERSTOKES2D; k++) {
	  buffer[buffer_size] = SolnBlk.WoN[i][k+1];
	  buffer_size++;
	}
      }
    }

    buffer_size = 2*NUM_VAR_NAVIERSTOKES2D*ni*nj;
    MPI::COMM_WORLD.Bcast(buffer,buffer_size,MPI::DOUBLE,0);

    if (!CFFC_Primary_MPI_Processor()) {
      buffer_size = 0;
      for (int i = SolnBlk.ICl-SolnBlk.Nghost; i <= SolnBlk.ICu+SolnBlk.Nghost; i++) {
	for (int k = 0; k < NUM_VAR_NAVIERSTOKES2D; k++) {
	  SolnBlk.WoS[i][k+1] = buffer[buffer_size];
	  buffer_size++;
	}
	for (int k = 0; k < NUM_VAR_NAVIERSTOKES2D; k++) {
	  SolnBlk.WoN[i][k+1] = buffer[buffer_size];
	  buffer_size++;
	}
      }
    }

    delete []buffer; buffer = NULL;

  }

#endif

}

#ifdef _MPI_VERSION
/**********************************************************************
 * Routine: Broadcast_Solution_Block                                  *
 *                                                                    *
 * Broadcast quadrilateral solution block to all processors           *
 * associated with the specified communicator from the specified      *
 * processor using the MPI broadcast routine.                         *
 *                                                                    *
 **********************************************************************/
void Broadcast_Solution_Block(NavierStokes2D_Quad_Block &SolnBlk,
                              MPI::Intracomm &Communicator,
                              const int Source_CPU) {

  int Source_Rank = 0;
  int ni, nj, ng, nr, block_allocated, buffer_size;
  double *buffer;

  // Broadcast the number of cells in each direction.
  if (CFFC_MPI::This_Processor_Number == Source_CPU) {
    ni = SolnBlk.NCi;
    nj = SolnBlk.NCj;
    ng = SolnBlk.Nghost;
    nr = SolnBlk.residual_variable;
    if (SolnBlk.U != NULL) {
      block_allocated = 1;
    } else {
      block_allocated = 0;
    } 
  }

  Communicator.Bcast(&ni,1,MPI::INT,Source_Rank);
  Communicator.Bcast(&nj,1,MPI::INT,Source_Rank);
  Communicator.Bcast(&ng,1,MPI::INT,Source_Rank);
  Communicator.Bcast(&nr,1,MPI::INT,Source_Rank);
  Communicator.Bcast(&block_allocated,1,MPI::INT,Source_Rank);

  // On non-source MPI processors, allocate (re-allocate) memory for the
  // quadrilateral solution block as necessary.
  if (!(CFFC_MPI::This_Processor_Number == Source_CPU)) {
    if (SolnBlk.NCi != ni || SolnBlk.NCj != nj || SolnBlk.Nghost != ng) { 
      if (SolnBlk.U != NULL) SolnBlk.deallocate();
      if (block_allocated) SolnBlk.allocate(ni-2*ng,nj-2*ng,ng); 
    }
    // Set the block static variables if they were not previously assigned.
    if (SolnBlk.residual_variable != nr) SolnBlk.residual_variable = nr;
  }

  // Broadcast the viscous flow indicator.
  Communicator.Bcast(&(SolnBlk.Flow_Type),1,MPI::INT,Source_Rank);

  // Broadcast the axisymmetric/planar flow indicator.
  Communicator.Bcast(&(SolnBlk.Axisymmetric),1,MPI::INT,Source_Rank);

  // Broadcast the compressibility effect correction indicator.
  Communicator.Bcast(&(SolnBlk.Compressibility_Effect),1,MPI::INT,Source_Rank);

  // Broadcast the transition model correction indicator.
  Communicator.Bcast(&(SolnBlk.Transition_Model),1,MPI::INT,Source_Rank);

  // Broadcast the variable Prandtl number indicator.
  Communicator.Bcast(&(SolnBlk.Variable_Prandtl),1,MPI::INT,Source_Rank);

  // Broadcast the wall velocity.
  Communicator.Bcast(&(SolnBlk.Vwall.x),1,MPI::DOUBLE,Source_Rank);
  Communicator.Bcast(&(SolnBlk.Vwall.y),1,MPI::DOUBLE,Source_Rank);

  // Broadcast the wall temperature.
  Communicator.Bcast(&(SolnBlk.Twall),1,MPI::DOUBLE,Source_Rank);

  // Broadcast the grid.
  Broadcast_Quad_Block(SolnBlk.Grid,Communicator,Source_CPU);

  // Broadcast the solution state variables.
  if (block_allocated) {
    ni = (SolnBlk.ICu+SolnBlk.Nghost) - (SolnBlk.ICl-SolnBlk.Nghost) + 1;
    nj = (SolnBlk.JCu+SolnBlk.Nghost) - (SolnBlk.JCl-SolnBlk.Nghost) + 1;
    buffer = new double[NUM_VAR_NAVIERSTOKES2D*ni*nj];

    if (CFFC_MPI::This_Processor_Number == Source_CPU) {
      buffer_size = 0;
      for (int j = SolnBlk.JCl-SolnBlk.Nghost; j <= SolnBlk.JCu+SolnBlk.Nghost; j++) {
	for (int i = SolnBlk.ICl-SolnBlk.Nghost; i <= SolnBlk.ICu+SolnBlk.Nghost; i++) {
	  for (int k = 0; k < NUM_VAR_NAVIERSTOKES2D; k++) {
	    buffer[buffer_size] = SolnBlk.U[i][j][k+1];
	    buffer_size++;
	  }
	}
      }
    }

    buffer_size = NUM_VAR_NAVIERSTOKES2D*ni*nj;
    Communicator.Bcast(buffer,buffer_size,MPI::DOUBLE,Source_Rank);

    if (!(CFFC_MPI::This_Processor_Number == Source_CPU)) {
      buffer_size = 0;
      for (int j = SolnBlk.JCl-SolnBlk.Nghost; j <= SolnBlk.JCu+SolnBlk.Nghost; j++) {
	for (int i = SolnBlk.ICl-SolnBlk.Nghost; i <= SolnBlk.ICu+SolnBlk.Nghost; i++) {
	  for (int k = 0; k < NUM_VAR_NAVIERSTOKES2D; k++) {
	    SolnBlk.U[i][j][k+1] = buffer[buffer_size];
	    buffer_size++;
	  }
	  SolnBlk.W[i][j] = W(SolnBlk.U[i][j]);
	}
      }
    }

    delete []buffer; buffer = NULL;

    ni = 1;
    nj = (SolnBlk.JCu+SolnBlk.Nghost) - (SolnBlk.JCl-SolnBlk.Nghost) + 1;
    buffer = new double[2*NUM_VAR_NAVIERSTOKES2D*ni*nj];

    if (CFFC_MPI::This_Processor_Number == Source_CPU) {
      buffer_size = 0;
      for (int j = SolnBlk.JCl-SolnBlk.Nghost; j <= SolnBlk.JCu+SolnBlk.Nghost; j++) {
	for (int k = 0; k < NUM_VAR_NAVIERSTOKES2D; k++) {
	  buffer[buffer_size] = SolnBlk.WoW[j][k+1];
	  buffer_size++;
	}
	for (int k = 0; k < NUM_VAR_NAVIERSTOKES2D; k++) {
	  buffer[buffer_size] = SolnBlk.WoW[j][k+1];
	  buffer_size++;
	}
      }
    }

    buffer_size = 2*NUM_VAR_NAVIERSTOKES2D*ni*nj;
    Communicator.Bcast(buffer,buffer_size,MPI::DOUBLE,Source_Rank);

    if (!(CFFC_MPI::This_Processor_Number == Source_CPU)) {
      buffer_size = 0;
      for (int j  = SolnBlk.JCl-SolnBlk.Nghost; j <= SolnBlk.JCu+SolnBlk.Nghost; j++) {
	for (int k = 0; k < NUM_VAR_NAVIERSTOKES2D; k++) {
	  SolnBlk.WoW[j][k+1] = buffer[buffer_size];
	  buffer_size++;
	}
	for (int k = 0; k < NUM_VAR_NAVIERSTOKES2D; k++) {
	  SolnBlk.WoE[j][k+1] = buffer[buffer_size];
	  buffer_size++;
	}
      }
    }

    delete []buffer; buffer = NULL;

    ni = (SolnBlk.ICu+SolnBlk.Nghost) - (SolnBlk.ICl-SolnBlk.Nghost) + 1;
    nj = 1;
    buffer = new double[2*NUM_VAR_NAVIERSTOKES2D*ni*nj];

    if (CFFC_MPI::This_Processor_Number == Source_CPU) {
      buffer_size = 0;
      for (int i  = SolnBlk.ICl-SolnBlk.Nghost; i <= SolnBlk.ICu+SolnBlk.Nghost; i++) {
	for (int k = 0; k < NUM_VAR_NAVIERSTOKES2D; k++) {
	  buffer[buffer_size] = SolnBlk.WoS[i][k+1];
	  buffer_size++;
	}
	for (int k = 0; k < NUM_VAR_NAVIERSTOKES2D; k++) {
	  buffer[buffer_size] = SolnBlk.WoN[i][k+1];
	  buffer_size++;
	}
      }
    }

    buffer_size = 2*NUM_VAR_NAVIERSTOKES2D*ni*nj;
    Communicator.Bcast(buffer,buffer_size,MPI::DOUBLE,Source_Rank);

    if (!(CFFC_MPI::This_Processor_Number == Source_CPU)) {
      buffer_size = 0;
      for (int i  = SolnBlk.ICl-SolnBlk.Nghost; i <= SolnBlk.ICu+SolnBlk.Nghost; i++) {
	for (int k = 0; k < NUM_VAR_NAVIERSTOKES2D; k++) {
	  SolnBlk.WoS[i][k+1] = buffer[buffer_size];
	  buffer_size++;
	}
	for (int k = 0; k < NUM_VAR_NAVIERSTOKES2D; k++) {
	  SolnBlk.WoN[i][k+1] = buffer[buffer_size];
	  buffer_size++;
	}
      }
    }

    delete []buffer; buffer = NULL;

  }

}
#endif

/**********************************************************************
 * Routine: Copy_Solution_Block                                       *
 *                                                                    *
 * Copies the solution information of quadrilateral solution block    *
 * SolnBlk2 to SolnBlk1.                                              *
 *                                                                    *
 **********************************************************************/
void Copy_Solution_Block(NavierStokes2D_Quad_Block &SolnBlk1,
                         NavierStokes2D_Quad_Block &SolnBlk2) {

  // Allocate (re-allocate) memory for the solution of the 
  // quadrilateral solution block SolnBlk1 as necessary.
  if (SolnBlk1.NCi    != SolnBlk2.NCi || 
      SolnBlk1.NCj    != SolnBlk2.NCj || 
      SolnBlk1.Nghost != SolnBlk2.Nghost) {
    if (SolnBlk1.U != NULL) SolnBlk1.deallocate();
    if (SolnBlk2.U != NULL) SolnBlk1.allocate(SolnBlk2.NCi-2*SolnBlk2.Nghost,
					      SolnBlk2.NCj-2*SolnBlk2.Nghost,
					      SolnBlk2.Nghost);
  }

  // Copy the viscous flow indicator.
  SolnBlk1.Flow_Type = SolnBlk2.Flow_Type;

  // Copy the axisymmetric/planar flow indicator.
  SolnBlk1.Axisymmetric = SolnBlk2.Axisymmetric;

  // Copy the compressibility effect correction indicator.
  SolnBlk1.Compressibility_Effect = SolnBlk2.Compressibility_Effect;

  // Copy the transition model correction indicator.
  SolnBlk1.Transition_Model = SolnBlk2.Transition_Model;

  // Copy the variable Prandtl number indicator.
  SolnBlk1.Variable_Prandtl = SolnBlk2.Variable_Prandtl;

  // Copy the wall velocity.
  SolnBlk1.Vwall = SolnBlk2.Vwall;

  // Copy the wall temperature.
  SolnBlk1.Twall = SolnBlk2.Twall;

  // Copy the grid of the second solution block to the first solution 
  // block.
  Copy_Quad_Block(SolnBlk1.Grid,SolnBlk2.Grid);

  // Copy the solution information from SolnBlk2 to SolnBlk1.
  if (SolnBlk2.U != NULL) {
    for (int j = SolnBlk1.JCl-SolnBlk1.Nghost; j <= SolnBlk1.JCu+SolnBlk1.Nghost; j++) {
      for (int i = SolnBlk1.ICl-SolnBlk1.Nghost; i <= SolnBlk1.ICu+SolnBlk1.Nghost; i++) {
	SolnBlk1.U[i][j] = SolnBlk2.U[i][j];
	SolnBlk1.W[i][j] = SolnBlk2.W[i][j];
	for (int k = 0; k < NUMBER_OF_RESIDUAL_VECTORS_NAVIERSTOKES2D; k++)
	  SolnBlk1.dUdt[i][j][k] = SolnBlk2.dUdt[i][j][k];
	SolnBlk1.dWdx[i][j] = SolnBlk2.dWdx[i][j];
	SolnBlk1.dWdy[i][j] = SolnBlk2.dWdy[i][j];
	SolnBlk1.d_dWdx_dW[i][j] = SolnBlk2.d_dWdx_dW[i][j];
	SolnBlk1.d_dWdy_dW[i][j] = SolnBlk2.d_dWdy_dW[i][j];
	SolnBlk1.phi[i][j]  = SolnBlk2.phi[i][j];
	SolnBlk1.Uo[i][j]   = SolnBlk2.Uo[i][j];
	SolnBlk1.dt[i][j]   = SolnBlk2.dt[i][j];
      }
    }

    for (int j = SolnBlk1.JCl-SolnBlk1.Nghost; j <= SolnBlk1.JCu+SolnBlk1.Nghost; j++) {
      SolnBlk1.WoW[j] = SolnBlk2.WoW[j];
      SolnBlk1.WoE[j] = SolnBlk2.WoE[j];
    }

    for (int i = SolnBlk1.ICl-SolnBlk1.Nghost; i <= SolnBlk1.ICu+SolnBlk1.Nghost; i++) {
      SolnBlk1.WoS[i] = SolnBlk2.WoS[i];
      SolnBlk1.WoN[i] = SolnBlk2.WoN[i];
    }
  }

}

/**********************************************************************
 * Routine: Prolong_Solution_Block                                    *
 *                                                                    *
 * Prolongs the solution information of one of the specified sectors  *
 * of the original quadrilateral solution block SolnBlk_Original to   *
 * the refined solution block SolnBlk_Fine.                           *
 *                                                                    *
 **********************************************************************/
int Prolong_Solution_Block(NavierStokes2D_Quad_Block &SolnBlk_Fine,
			   NavierStokes2D_Quad_Block &SolnBlk_Original,
			   const int Sector) {

  int error_flag;
  int i_min, i_max, j_min, j_max, mesh_refinement_permitted;
  int i_fine, j_fine, i_coarse_min, j_coarse_min, coarse_cell_found;
  double distance, area_total_fine;
  Vector2D dX;
  NavierStokes2D_cState NavierStokes2D_U_STDATM, Ucoarse;
  NavierStokes2D_pState NavierStokes2D_W_STDATM;

  // Allocate (re-allocate) memory for the solution of the refined
  // quadrilateral solution block as necessary.
  if ((SolnBlk_Original.NCi-2*SolnBlk_Original.Nghost-2*((SolnBlk_Original.NCi-2*SolnBlk_Original.Nghost)/2) != 0) || 
      (SolnBlk_Original.NCj-2*SolnBlk_Original.Nghost-2*((SolnBlk_Original.NCj-2*SolnBlk_Original.Nghost)/2) != 0) ||
      (SolnBlk_Original.NCi-2*SolnBlk_Original.Nghost < 2*SolnBlk_Original.Nghost) ||
      (SolnBlk_Original.NCj-2*SolnBlk_Original.Nghost < 2*SolnBlk_Original.Nghost) ||
      (SolnBlk_Original.Grid.Node == NULL) ||
      (SolnBlk_Original.U == NULL)) {
    mesh_refinement_permitted = 0;
  } else {
    mesh_refinement_permitted = 1;
    if (SolnBlk_Fine.NCi != SolnBlk_Original.NCi || 
	SolnBlk_Fine.NCj != SolnBlk_Original.NCj ||
	SolnBlk_Fine.Nghost != SolnBlk_Original.Nghost) {
      if (SolnBlk_Fine.U != NULL) SolnBlk_Fine.deallocate();
      if (SolnBlk_Original.U != NULL)
	SolnBlk_Fine.allocate(SolnBlk_Original.NCi-2*SolnBlk_Original.Nghost,
			      SolnBlk_Original.NCj-2*SolnBlk_Original.Nghost,
			      SolnBlk_Original.Nghost);
      // In this case, create the refined mesh.
      Refine_Mesh(SolnBlk_Fine.Grid,SolnBlk_Original.Grid,Sector);
    }
  }

  if (mesh_refinement_permitted) {

    // Copy the viscous flow indicator.
    SolnBlk_Fine.Flow_Type = SolnBlk_Original.Flow_Type;

    // Copy the axisymmetric/planar flow indicator.
    SolnBlk_Fine.Axisymmetric = SolnBlk_Original.Axisymmetric;

    // Copy the compressibility effect correction indicator.
    SolnBlk_Fine.Compressibility_Effect = SolnBlk_Original.Compressibility_Effect;

    // Copy the transition model correction indicator.
    SolnBlk_Fine.Transition_Model = SolnBlk_Original.Transition_Model;
    
    // Copy the variable Prandtl number indicator.
    SolnBlk_Fine.Variable_Prandtl = SolnBlk_Original.Variable_Prandtl;

    // Copy the wall velocity.
    SolnBlk_Fine.Vwall = SolnBlk_Original.Vwall;

    // Copy the wall temperature.
    SolnBlk_Fine.Twall = SolnBlk_Original.Twall;

    // Prolong the solution information from original solution block
    // to the refined solution block.
    switch(Sector) {
    case GRID2D_QUAD_BLOCK_SECTOR_NW :
      i_min = SolnBlk_Original.ICl;
      i_max = SolnBlk_Original.ICl+(SolnBlk_Original.ICu-SolnBlk_Original.ICl-1)/2;
      j_min = SolnBlk_Original.JCl+(SolnBlk_Original.JCu-SolnBlk_Original.JCl+1)/2; 
      j_max = SolnBlk_Original.JCu;
      break;
    case GRID2D_QUAD_BLOCK_SECTOR_NE :
      i_min = SolnBlk_Original.ICl+(SolnBlk_Original.ICu-SolnBlk_Original.ICl+1)/2;
      i_max = SolnBlk_Original.ICu;
      j_min = SolnBlk_Original.JCl+(SolnBlk_Original.JCu-SolnBlk_Original.JCl+1)/2; 
      j_max = SolnBlk_Original.JCu;
      break;
    case GRID2D_QUAD_BLOCK_SECTOR_SE :
      i_min = SolnBlk_Original.ICl+(SolnBlk_Original.ICu-SolnBlk_Original.ICl+1)/2;
      i_max = SolnBlk_Original.ICu;
      j_min = SolnBlk_Original.JCl; 
      j_max = SolnBlk_Original.JCl+(SolnBlk_Original.JCu-SolnBlk_Original.JCl-1)/2;
      break;
    case GRID2D_QUAD_BLOCK_SECTOR_SW :
      i_min = SolnBlk_Original.ICl;
      i_max = SolnBlk_Original.ICl+(SolnBlk_Original.ICu-SolnBlk_Original.ICl-1)/2;
      j_min = SolnBlk_Original.JCl; 
      j_max = SolnBlk_Original.JCl+(SolnBlk_Original.JCu-SolnBlk_Original.JCl-1)/2;
      break;
    default:
      i_min = SolnBlk_Original.ICl;
      i_max = SolnBlk_Original.ICl+(SolnBlk_Original.ICu-SolnBlk_Original.ICl-1)/2;
      j_min = SolnBlk_Original.JCl+(SolnBlk_Original.JCu-SolnBlk_Original.JCl+1)/2; 
      j_max = SolnBlk_Original.JCu;
      break;
    };

    for (int j = j_min; j <= j_max; j++) {
      for (int i = i_min ; i <= i_max; i++) {
	//area_total_fine = SolnBlk_Fine.Grid.Cell[2*(i-i_min)+SolnBlk_Fine.ICl  ]
	//                                      [2*(j-j_min)+SolnBlk_Fine.JCl  ].A+
	//                SolnBlk_Fine.Grid.Cell[2*(i-i_min)+SolnBlk_Fine.ICl+1]
	//                                      [2*(j-j_min)+SolnBlk_Fine.JCl  ].A+
	//                SolnBlk_Fine.Grid.Cell[2*(i-i_min)+SolnBlk_Fine.ICl  ]
	//                                      [2*(j-j_min)+SolnBlk_Fine.JCl+1].A+
	//                SolnBlk_Fine.Grid.Cell[2*(i-i_min)+SolnBlk_Fine.ICl+1]
	//                                      [2*(j-j_min)+SolnBlk_Fine.JCl+1].A;
	SolnBlk_Fine.U[2*(i-i_min)+SolnBlk_Fine.ICl  ][2*(j-j_min)+SolnBlk_Fine.JCl  ] = SolnBlk_Original.U[i][j];
	//= (SolnBlk_Original.Grid.Cell[i][j].A/area_total_fine)*SolnBlk_Original.U[i][j];
	SolnBlk_Fine.W[2*(i-i_min)+SolnBlk_Fine.ICl  ][2*(j-j_min)+SolnBlk_Fine.JCl  ] = W(SolnBlk_Fine.U[2*(i-i_min)+SolnBlk_Fine.ICl  ][2*(j-j_min)+SolnBlk_Fine.JCl  ]);
	SolnBlk_Fine.U[2*(i-i_min)+SolnBlk_Fine.ICl+1][2*(j-j_min)+SolnBlk_Fine.JCl  ] = SolnBlk_Original.U[i][j];
	//= (SolnBlk_Original.Grid.Cell[i][j].A/area_total_fine)*SolnBlk_Original.U[i][j];
	SolnBlk_Fine.W[2*(i-i_min)+SolnBlk_Fine.ICl+1][2*(j-j_min)+SolnBlk_Fine.JCl  ] = W(SolnBlk_Fine.U[2*(i-i_min)+SolnBlk_Fine.ICl+1][2*(j-j_min)+SolnBlk_Fine.JCl  ]);
	SolnBlk_Fine.U[2*(i-i_min)+SolnBlk_Fine.ICl  ][2*(j-j_min)+SolnBlk_Fine.JCl+1] = SolnBlk_Original.U[i][j];
	//= (SolnBlk_Original.Grid.Cell[i][j].A/area_total_fine)*SolnBlk_Original.U[i][j];
	SolnBlk_Fine.W[2*(i-i_min)+SolnBlk_Fine.ICl  ][2*(j-j_min)+SolnBlk_Fine.JCl+1] = W(SolnBlk_Fine.U[2*(i-i_min)+SolnBlk_Fine.ICl  ][2*(j-j_min)+SolnBlk_Fine.JCl+1]);
	SolnBlk_Fine.U[2*(i-i_min)+SolnBlk_Fine.ICl+1][2*(j-j_min)+SolnBlk_Fine.JCl+1] = SolnBlk_Original.U[i][j];
	//= (SolnBlk_Original.Grid.Cell[i][j].A/area_total_fine)*SolnBlk_Original.U[i][j];
	SolnBlk_Fine.W[2*(i-i_min)+SolnBlk_Fine.ICl+1][2*(j-j_min)+SolnBlk_Fine.JCl+1] = W(SolnBlk_Fine.U[2*(i-i_min)+SolnBlk_Fine.ICl+1][2*(j-j_min)+SolnBlk_Fine.JCl+1]);
// 	if (SolnBlk_Fine.W[2*(i-i_min)+SolnBlk_Fine.ICl  ][2*(j-j_min)+SolnBlk_Fine.JCl  ].rho < ZERO ||
// 	    SolnBlk_Fine.W[2*(i-i_min)+SolnBlk_Fine.ICl+1][2*(j-j_min)+SolnBlk_Fine.JCl  ].rho < ZERO ||
// 	    SolnBlk_Fine.W[2*(i-i_min)+SolnBlk_Fine.ICl  ][2*(j-j_min)+SolnBlk_Fine.JCl+1].rho < ZERO ||
// 	    SolnBlk_Fine.W[2*(i-i_min)+SolnBlk_Fine.ICl+1][2*(j-j_min)+SolnBlk_Fine.JCl+1].rho < ZERO) {
// 	  for (int jj = SolnBlk_Original.JCl; jj <= SolnBlk_Original.JCu; jj++) {
// 	    for (int ii = SolnBlk_Original.ICl; ii <= SolnBlk_Original.ICu; ii++) {
// 	      cout << endl << SolnBlk_Original.Grid.Cell[ii][jj].Xc << " " << SolnBlk_Original.W[ii][jj];
// 	    }
// 	  }
// // 	  cout << endl << i << " " << j << " " << SolnBlk_Original.JCu
// // 	       << endl << SolnBlk_Original.Grid.Cell[i][j].Xc
// // 	       << endl << SolnBlk_Original.U[i][j]
// // 	       << endl << SolnBlk_Original.W[i][j]
// // 	       << endl << SolnBlk_Fine.Grid.Cell[2*(i-i_min)+SolnBlk_Fine.ICl  ][2*(j-j_min)+SolnBlk_Fine.JCl  ].Xc
// // 	       << endl << SolnBlk_Fine.U[2*(i-i_min)+SolnBlk_Fine.ICl  ][2*(j-j_min)+SolnBlk_Fine.JCl  ]
// // 	       << endl << SolnBlk_Fine.W[2*(i-i_min)+SolnBlk_Fine.ICl  ][2*(j-j_min)+SolnBlk_Fine.JCl  ]
// // 	       << endl << SolnBlk_Fine.Grid.Cell[2*(i-i_min)+SolnBlk_Fine.ICl+1][2*(j-j_min)+SolnBlk_Fine.JCl  ].Xc
// // 	       << endl << SolnBlk_Fine.U[2*(i-i_min)+SolnBlk_Fine.ICl+1][2*(j-j_min)+SolnBlk_Fine.JCl  ]
// // 	       << endl << SolnBlk_Fine.W[2*(i-i_min)+SolnBlk_Fine.ICl+1][2*(j-j_min)+SolnBlk_Fine.JCl  ]
// // 	       << endl << SolnBlk_Fine.Grid.Cell[2*(i-i_min)+SolnBlk_Fine.ICl  ][2*(j-j_min)+SolnBlk_Fine.JCl+1].Xc
// // 	       << endl << SolnBlk_Fine.U[2*(i-i_min)+SolnBlk_Fine.ICl  ][2*(j-j_min)+SolnBlk_Fine.JCl+1]
// // 	       << endl << SolnBlk_Fine.W[2*(i-i_min)+SolnBlk_Fine.ICl  ][2*(j-j_min)+SolnBlk_Fine.JCl+1]
// // 	       << endl << SolnBlk_Fine.Grid.Cell[2*(i-i_min)+SolnBlk_Fine.ICl+1][2*(j-j_min)+SolnBlk_Fine.JCl+1].Xc
// // 	       << endl << SolnBlk_Fine.U[2*(i-i_min)+SolnBlk_Fine.ICl+1][2*(j-j_min)+SolnBlk_Fine.JCl+1]
// // 	       << endl << SolnBlk_Fine.W[2*(i-i_min)+SolnBlk_Fine.ICl+1][2*(j-j_min)+SolnBlk_Fine.JCl+1];
// 	}
// 	if (SolnBlk_Original.U[i][j].rho < ZERO) {
// 	  if (SolnBlk_Original.U[i-1][j].rho > ZERO && i-1 >= SolnBlk_Original.ICl) Ucoarse = SolnBlk_Original.U[i-1][j];
// 	  else if (SolnBlk_Original.U[i-1][j-1].rho > ZERO && i-1 >= SolnBlk_Original.ICl && j-1 >= SolnBlk_Original.JCl) Ucoarse = SolnBlk_Original.U[i-1][j-1];
// 	  else if (SolnBlk_Original.U[i  ][j-1].rho > ZERO                                && j-1 >= SolnBlk_Original.JCl) Ucoarse = SolnBlk_Original.U[i  ][j-1];
// 	  else if (SolnBlk_Original.U[i+1][j-1].rho > ZERO && i+1 <= SolnBlk_Original.ICu && j-1 >= SolnBlk_Original.JCl) Ucoarse = SolnBlk_Original.U[i+1][j-1];
// 	  else if (SolnBlk_Original.U[i+1][j  ].rho > ZERO && i+1 <= SolnBlk_Original.ICu                               ) Ucoarse = SolnBlk_Original.U[i+1][j  ];
// 	  else if (SolnBlk_Original.U[i+1][j+1].rho > ZERO && i+1 <= SolnBlk_Original.ICu && j+1 <= SolnBlk_Original.JCu) Ucoarse = SolnBlk_Original.U[i+1][j+1];
// 	  else if (SolnBlk_Original.U[i  ][j+1].rho > ZERO                                && j+1 <= SolnBlk_Original.JCu) Ucoarse = SolnBlk_Original.U[i  ][j+1];
// 	  else if (SolnBlk_Original.U[i-1][j+1].rho > ZERO && i-1 >= SolnBlk_Original.ICl && j+1 <= SolnBlk_Original.JCu) Ucoarse = SolnBlk_Original.U[i-1][j+1];
// 	  else return 1;
// 	  SolnBlk_Fine.U[2*(i-i_min)+SolnBlk_Fine.ICl  ][2*(j-j_min)+SolnBlk_Fine.JCl  ] = Ucoarse;
// 	  //= (SolnBlk_Original.Grid.Cell[i][j].A/area_total_fine)*SolnBlk_Original.U[i][j];
// 	  SolnBlk_Fine.W[2*(i-i_min)+SolnBlk_Fine.ICl  ][2*(j-j_min)+SolnBlk_Fine.JCl  ] = W(SolnBlk_Fine.U[2*(i-i_min)+SolnBlk_Fine.ICl  ][2*(j-j_min)+SolnBlk_Fine.JCl  ]);
// 	  SolnBlk_Fine.U[2*(i-i_min)+SolnBlk_Fine.ICl+1][2*(j-j_min)+SolnBlk_Fine.JCl  ] = Ucoarse;
// 	  //= (SolnBlk_Original.Grid.Cell[i][j].A/area_total_fine)*SolnBlk_Original.U[i][j];
// 	  SolnBlk_Fine.W[2*(i-i_min)+SolnBlk_Fine.ICl+1][2*(j-j_min)+SolnBlk_Fine.JCl  ] = W(SolnBlk_Fine.U[2*(i-i_min)+SolnBlk_Fine.ICl+1][2*(j-j_min)+SolnBlk_Fine.JCl  ]);
// 	  SolnBlk_Fine.U[2*(i-i_min)+SolnBlk_Fine.ICl  ][2*(j-j_min)+SolnBlk_Fine.JCl+1] = Ucoarse;
// 	  //= (SolnBlk_Original.Grid.Cell[i][j].A/area_total_fine)*SolnBlk_Original.U[i][j];
// 	  SolnBlk_Fine.W[2*(i-i_min)+SolnBlk_Fine.ICl  ][2*(j-j_min)+SolnBlk_Fine.JCl+1] = W(SolnBlk_Fine.U[2*(i-i_min)+SolnBlk_Fine.ICl  ][2*(j-j_min)+SolnBlk_Fine.JCl+1]);
// 	  SolnBlk_Fine.U[2*(i-i_min)+SolnBlk_Fine.ICl+1][2*(j-j_min)+SolnBlk_Fine.JCl+1] = Ucoarse;
// 	  //= (SolnBlk_Original.Grid.Cell[i][j].A/area_total_fine)*SolnBlk_Original.U[i][j];
// 	  SolnBlk_Fine.W[2*(i-i_min)+SolnBlk_Fine.ICl+1][2*(j-j_min)+SolnBlk_Fine.JCl+1] = W(SolnBlk_Fine.U[2*(i-i_min)+SolnBlk_Fine.ICl+1][2*(j-j_min)+SolnBlk_Fine.JCl+1]);
// 	}

      }
    }

    // Prolong the east and west boundary states.
    for (int j = j_min-SolnBlk_Original.Nghost/2; j <= j_max+SolnBlk_Original.Nghost/2; j++) {
      SolnBlk_Fine.WoW[2*(j-j_min)+SolnBlk_Fine.JCl  ] = SolnBlk_Original.WoW[j];
      SolnBlk_Fine.WoW[2*(j-j_min)+SolnBlk_Fine.JCl+1] = SolnBlk_Original.WoW[j];
      SolnBlk_Fine.WoE[2*(j-j_min)+SolnBlk_Fine.JCl  ] = SolnBlk_Original.WoE[j];
      SolnBlk_Fine.WoE[2*(j-j_min)+SolnBlk_Fine.JCl+1] = SolnBlk_Original.WoE[j];
    }

    // Prolong the north and south boundary states.
    for (int i = i_min-SolnBlk_Original.Nghost/2; i <= i_max+SolnBlk_Original.Nghost/2; i++) {
      SolnBlk_Fine.WoS[2*(i-i_min)+SolnBlk_Fine.ICl  ] = SolnBlk_Original.WoS[i];
      SolnBlk_Fine.WoS[2*(i-i_min)+SolnBlk_Fine.ICl+1] = SolnBlk_Original.WoS[i];
      SolnBlk_Fine.WoN[2*(i-i_min)+SolnBlk_Fine.ICl  ] = SolnBlk_Original.WoN[i];
      SolnBlk_Fine.WoN[2*(i-i_min)+SolnBlk_Fine.ICl+1] = SolnBlk_Original.WoN[i];
    }

  }

  // Prolongation of solution block was successful.
  return 0;

}

/**********************************************************************
 * Routine: Restrict_Solution_Block                                   *
 *                                                                    *
 * Restricts the solution information of the original quadrilateral   *
 * solution block SolnBlk_Original to the coarsened solution block    *
 * SolnBlk_Coarse.                                                    *
 *                                                                    *
 **********************************************************************/
int Restrict_Solution_Block(NavierStokes2D_Quad_Block &SolnBlk_Coarse,
			    NavierStokes2D_Quad_Block &SolnBlk_Original_SW,
			    NavierStokes2D_Quad_Block &SolnBlk_Original_SE,
			    NavierStokes2D_Quad_Block &SolnBlk_Original_NW,
			    NavierStokes2D_Quad_Block &SolnBlk_Original_NE) {

  int error_flag;
  int i_coarse, j_coarse, mesh_coarsening_permitted;
  double A, total_area;
  Polygon Pc, Pf;
  NavierStokes2D_pState NavierStokes2D_W_STDATM; NavierStokes2D_W_STDATM.Standard_Atmosphere();
  NavierStokes2D_cState NavierStokes2D_U_STDATM; NavierStokes2D_U_STDATM.Standard_Atmosphere();
  NavierStokes2D_cState NavierStokes2D_U_VACUUM; NavierStokes2D_U_VACUUM.Vacuum();

  // Allocate memory for the cells and nodes for the coarsened 
  // quadrilateral mesh block.

  if ((SolnBlk_Original_SW.NCi-2*SolnBlk_Original_SW.Nghost-
       2*((SolnBlk_Original_SW.NCi-2*SolnBlk_Original_SW.Nghost)/2) != 0) || 
      (SolnBlk_Original_SW.NCj-2*SolnBlk_Original_SW.Nghost-
       2*((SolnBlk_Original_SW.NCj-2*SolnBlk_Original_SW.Nghost)/2) != 0) ||
      (SolnBlk_Original_SW.NCi-2*SolnBlk_Original_SW.Nghost < 4) ||
      (SolnBlk_Original_SW.NCj-2*SolnBlk_Original_SW.Nghost < 4) ||
      (SolnBlk_Original_SE.NCi != SolnBlk_Original_SW.NCi) ||
      (SolnBlk_Original_SE.NCj != SolnBlk_Original_SW.NCj) ||
      (SolnBlk_Original_NW.NCi != SolnBlk_Original_SW.NCi) ||
      (SolnBlk_Original_NW.NCj != SolnBlk_Original_SW.NCj) ||
      (SolnBlk_Original_NE.NCi != SolnBlk_Original_SW.NCi) ||
      (SolnBlk_Original_NE.NCj != SolnBlk_Original_SW.NCj) ||
      (SolnBlk_Original_SW.Grid.Node == NULL) ||
      (SolnBlk_Original_SE.Grid.Node == NULL) ||
      (SolnBlk_Original_NW.Grid.Node == NULL) ||
      (SolnBlk_Original_NE.Grid.Node == NULL) ||
      (SolnBlk_Original_SW.U == NULL) ||
      (SolnBlk_Original_SE.U == NULL) ||
      (SolnBlk_Original_NW.U == NULL) ||
      (SolnBlk_Original_NE.U == NULL)) {
    mesh_coarsening_permitted = 0;
  } else {
    mesh_coarsening_permitted = 1;
    if (SolnBlk_Coarse.NCi != SolnBlk_Original_SW.NCi || 
	SolnBlk_Coarse.NCj != SolnBlk_Original_SW.NCj ||
	SolnBlk_Coarse.Nghost != SolnBlk_Original_SW.Nghost) {
      if (SolnBlk_Coarse.U != NULL) SolnBlk_Coarse.deallocate();
      if (SolnBlk_Original_SW.U != NULL)
	SolnBlk_Coarse.allocate(SolnBlk_Original_SW.NCi-2*SolnBlk_Original_SW.Nghost,
				SolnBlk_Original_SW.NCj-2*SolnBlk_Original_SW.Nghost,
				SolnBlk_Original_SW.Nghost);
      // In this case, create the coarse mesh.
      Coarsen_Mesh(SolnBlk_Coarse.Grid,
		   SolnBlk_Original_SW.Grid,
		   SolnBlk_Original_SE.Grid,
		   SolnBlk_Original_NW.Grid,
		   SolnBlk_Original_NE.Grid);
    }
  }

  if (mesh_coarsening_permitted) {

    // Copy the viscous flow indicator.
    SolnBlk_Coarse.Flow_Type = SolnBlk_Original_SW.Flow_Type;

    // Copy the axisymmetric/planar flow indicator.
    SolnBlk_Coarse.Axisymmetric = SolnBlk_Original_SW.Axisymmetric;

    // Copy the compressibility effect correction indicator.
    SolnBlk_Coarse.Compressibility_Effect = SolnBlk_Original_SW.Compressibility_Effect;

    // Copy the transition model correction indicator.
    SolnBlk_Coarse.Transition_Model = SolnBlk_Original_SW.Transition_Model;
    
    // Copy the variable Prandtl number indicator.
    SolnBlk_Coarse.Variable_Prandtl = SolnBlk_Original_SW.Variable_Prandtl;

    // Copy the wall velocity.
    SolnBlk_Coarse.Vwall = SolnBlk_Original_SW.Vwall;

    // Copy the wall temperature.
    SolnBlk_Coarse.Twall = SolnBlk_Original_SW.Twall;

    // Restrict the solution information from the four original
    // solution blocks to the newly coarsened solution block.

    // South-West corner fine block:	
    for (int j_fine = SolnBlk_Original_SW.JCl; j_fine <= SolnBlk_Original_SW.JCu; j_fine += 2) {
      for (int i_fine = SolnBlk_Original_SW.ICl; i_fine <= SolnBlk_Original_SW.ICu; i_fine += 2) {

	i_coarse = (i_fine-SolnBlk_Original_SW.ICl)/2 + SolnBlk_Coarse.ICl;
	j_coarse = (j_fine-SolnBlk_Original_SW.JCl)/2 + SolnBlk_Coarse.JCl;

	// Restrict the solution information from the fine cells to the
	// coarse cell using an area-weighted average of the fine cell
	// solution information.
	SolnBlk_Coarse.U[i_coarse][j_coarse] = (SolnBlk_Original_SW.Grid.Cell[i_fine  ][j_fine  ].A*
						SolnBlk_Original_SW.U[i_fine  ][j_fine  ] +
						SolnBlk_Original_SW.Grid.Cell[i_fine+1][j_fine  ].A*
						SolnBlk_Original_SW.U[i_fine+1][j_fine  ] + 
						SolnBlk_Original_SW.Grid.Cell[i_fine  ][j_fine+1].A*
						SolnBlk_Original_SW.U[i_fine  ][j_fine+1] +
						SolnBlk_Original_SW.Grid.Cell[i_fine+1][j_fine+1].A*
						SolnBlk_Original_SW.U[i_fine+1][j_fine+1])/
                                               (SolnBlk_Original_SW.Grid.Cell[i_fine  ][j_fine  ].A +
						SolnBlk_Original_SW.Grid.Cell[i_fine+1][j_fine  ].A +
						SolnBlk_Original_SW.Grid.Cell[i_fine  ][j_fine+1].A +
						SolnBlk_Original_SW.Grid.Cell[i_fine+1][j_fine+1].A);
                                               //SolnBlk_Coarse.Grid.Cell[i_coarse][j_coarse].A;
	SolnBlk_Coarse.W[i_coarse][j_coarse] = W(SolnBlk_Coarse.U[i_coarse][j_coarse]);
// 	if (SolnBlk_Coarse.U[i_coarse][j_coarse].rho < ZERO) {
// 	  cout << endl << SolnBlk_Coarse.Grid.Cell[i_coarse][j_coarse].Xc;
// 	}

      }
    }

    for (int j_fine = SolnBlk_Original_SW.JCl-SolnBlk_Original_SW.Nghost; 
	 j_fine <= SolnBlk_Original_SW.JCu+SolnBlk_Original_SW.Nghost; j_fine += 2) {
      j_coarse = (j_fine-SolnBlk_Original_SW.JCl)/2+SolnBlk_Coarse.JCl;
      SolnBlk_Coarse.WoW[j_coarse] = SolnBlk_Original_SW.WoW[j_fine];
      if (j_fine == SolnBlk_Original_SW.JCl-SolnBlk_Original_SW.Nghost) {
	SolnBlk_Coarse.WoW[j_coarse-1] = SolnBlk_Original_SW.WoW[j_fine];
      }
    }

    for (int i_fine = SolnBlk_Original_SW.ICl-SolnBlk_Original_SW.Nghost; 
	 i_fine <= SolnBlk_Original_SW.ICu+SolnBlk_Original_SW.Nghost; i_fine += 2) {
      i_coarse = (i_fine-SolnBlk_Original_SW.ICl)/2+SolnBlk_Coarse.ICl;
      SolnBlk_Coarse.WoS[i_coarse] = SolnBlk_Original_SW.WoS[i_fine];
      if (i_fine == SolnBlk_Original_SW.ICl-SolnBlk_Original_SW.Nghost) {
	SolnBlk_Coarse.WoS[i_coarse-1] = SolnBlk_Original_SW.WoS[i_fine];
      }
    }

    // South-East corner fine block:
    for (int j_fine = SolnBlk_Original_SE.JCl; j_fine <= SolnBlk_Original_SE.JCu; j_fine += 2) {
      for (int i_fine = SolnBlk_Original_SE.ICl; i_fine <= SolnBlk_Original_SE.ICu; i_fine += 2) {

	i_coarse = (i_fine-SolnBlk_Original_SE.ICl)/2 + (SolnBlk_Coarse.ICu-SolnBlk_Coarse.ICl+1)/2+SolnBlk_Coarse.ICl;
	j_coarse = (j_fine-SolnBlk_Original_SE.JCl)/2 + SolnBlk_Coarse.JCl;

	// Restrict the solution information from the fine cells to the
	// coarse cell using an area-weighted average of the fine cell
	// solution information.
	SolnBlk_Coarse.U[i_coarse][j_coarse] = (SolnBlk_Original_SE.Grid.Cell[i_fine  ][j_fine  ].A*
						SolnBlk_Original_SE.U[i_fine  ][j_fine  ] +
						SolnBlk_Original_SE.Grid.Cell[i_fine+1][j_fine  ].A*
						SolnBlk_Original_SE.U[i_fine+1][j_fine  ] + 
						SolnBlk_Original_SE.Grid.Cell[i_fine  ][j_fine+1].A*
						SolnBlk_Original_SE.U[i_fine  ][j_fine+1] +
						SolnBlk_Original_SE.Grid.Cell[i_fine+1][j_fine+1].A*
						SolnBlk_Original_SE.U[i_fine+1][j_fine+1])/
	                                       (SolnBlk_Original_SE.Grid.Cell[i_fine  ][j_fine  ].A +
						SolnBlk_Original_SE.Grid.Cell[i_fine+1][j_fine  ].A +
						SolnBlk_Original_SE.Grid.Cell[i_fine  ][j_fine+1].A +
						SolnBlk_Original_SE.Grid.Cell[i_fine+1][j_fine+1].A);
                                               //SolnBlk_Coarse.Grid.Cell[i_coarse][j_coarse].A;
	SolnBlk_Coarse.W[i_coarse][j_coarse] = W(SolnBlk_Coarse.U[i_coarse][j_coarse]);

      }
    }

    for (int j_fine = SolnBlk_Original_SE.JCl-SolnBlk_Original_SE.Nghost; 
	 j_fine <= SolnBlk_Original_SE.JCu+SolnBlk_Original_SE.Nghost; j_fine += 2) {
      j_coarse = (j_fine-SolnBlk_Original_SE.JCl)/2+SolnBlk_Coarse.JCl;
      SolnBlk_Coarse.WoE[j_coarse] = SolnBlk_Original_SE.WoE[j_fine];
      if (j_fine == SolnBlk_Original_SE.JCl-SolnBlk_Original_SE.Nghost) {
	SolnBlk_Coarse.WoE[j_coarse-1] = SolnBlk_Original_SE.WoE[j_fine];
      }
    }

    for (int i_fine = SolnBlk_Original_SE.ICl-SolnBlk_Original_SE.Nghost; 
	 i_fine <= SolnBlk_Original_SE.ICu+SolnBlk_Original_SE.Nghost; i_fine += 2) {
      i_coarse = (i_fine-SolnBlk_Original_SE.ICl)/2+(SolnBlk_Coarse.ICu-SolnBlk_Coarse.ICl+1)/2+SolnBlk_Coarse.ICl;
      SolnBlk_Coarse.WoS[i_coarse] = SolnBlk_Original_SE.WoS[i_fine];
      if (i_fine == SolnBlk_Original_SE.ICu+SolnBlk_Original_SE.Nghost) {
	SolnBlk_Coarse.WoS[i_coarse+1] = SolnBlk_Original_SE.WoS[i_fine];
      }
    }

    // North-West corner fine block:
    for (int j_fine = SolnBlk_Original_NW.JCl; j_fine <= SolnBlk_Original_NW.JCu; j_fine += 2) {
      for (int i_fine = SolnBlk_Original_NW.ICl; i_fine <= SolnBlk_Original_NW.ICu; i_fine += 2) {

	i_coarse = (i_fine-SolnBlk_Original_NW.ICl)/2 + SolnBlk_Coarse.ICl;
	j_coarse = (j_fine-SolnBlk_Original_NW.JCl)/2 + (SolnBlk_Coarse.JCu-SolnBlk_Coarse.JCl+1)/2+SolnBlk_Coarse.JCl;

	// Restrict the solution information from the fine cells to the
	// coarse cell using an area-weighted average of the fine cell
	// solution information.
	SolnBlk_Coarse.U[i_coarse][j_coarse] = (SolnBlk_Original_NW.Grid.Cell[i_fine  ][j_fine  ].A*
						SolnBlk_Original_NW.U[i_fine  ][j_fine  ] +
						SolnBlk_Original_NW.Grid.Cell[i_fine+1][j_fine  ].A*
						SolnBlk_Original_NW.U[i_fine+1][j_fine  ] + 
						SolnBlk_Original_NW.Grid.Cell[i_fine  ][j_fine+1].A*
						SolnBlk_Original_NW.U[i_fine  ][j_fine+1] +
						SolnBlk_Original_NW.Grid.Cell[i_fine+1][j_fine+1].A*
						SolnBlk_Original_NW.U[i_fine+1][j_fine+1])/
                                               (SolnBlk_Original_NW.Grid.Cell[i_fine  ][j_fine  ].A +
						SolnBlk_Original_NW.Grid.Cell[i_fine+1][j_fine  ].A +
						SolnBlk_Original_NW.Grid.Cell[i_fine  ][j_fine+1].A +
						SolnBlk_Original_NW.Grid.Cell[i_fine+1][j_fine+1].A);
                                               //SolnBlk_Coarse.Grid.Cell[i_coarse][j_coarse].A;
	SolnBlk_Coarse.W[i_coarse][j_coarse] = W(SolnBlk_Coarse.U[i_coarse][j_coarse]);

      }
    }

    for (int j_fine = SolnBlk_Original_NW.JCl-SolnBlk_Original_NW.Nghost; 
	 j_fine <= SolnBlk_Original_NW.JCu+SolnBlk_Original_NW.Nghost; j_fine += 2) {
      j_coarse = (j_fine-SolnBlk_Original_NW.JCl)/2 + (SolnBlk_Coarse.JCu-SolnBlk_Coarse.JCl+1)/2+SolnBlk_Coarse.JCl;
      SolnBlk_Coarse.WoW[j_coarse] = SolnBlk_Original_NW.WoW[j_fine];
      if (j_fine == SolnBlk_Original_NW.JCu+SolnBlk_Original_NW.Nghost) {
	SolnBlk_Coarse.WoW[j_coarse+1] = SolnBlk_Original_NW.WoW[j_fine];
      }
    }

    for (int i_fine = SolnBlk_Original_NW.ICl-SolnBlk_Original_NW.Nghost; 
	 i_fine <= SolnBlk_Original_NW.ICu+SolnBlk_Original_NW.Nghost; i_fine += 2) {
      i_coarse = (i_fine-SolnBlk_Original_NW.ICl)/2+SolnBlk_Coarse.ICl;
      SolnBlk_Coarse.WoN[i_coarse] = SolnBlk_Original_NW.WoN[i_fine];
      if (i_fine == SolnBlk_Original_NW.ICl-SolnBlk_Original_NW.Nghost) {
	SolnBlk_Coarse.WoN[i_coarse-1] = SolnBlk_Original_NW.WoN[i_fine];
      }
    }

    // North-East corner fine block:
    for (int j_fine = SolnBlk_Original_NE.JCl; j_fine <= SolnBlk_Original_NE.JCu; j_fine += 2) {
      for (int i_fine = SolnBlk_Original_NE.ICl; i_fine <= SolnBlk_Original_NE.ICu; i_fine += 2) {

	i_coarse = (i_fine-SolnBlk_Original_NE.ICl)/2 + (SolnBlk_Coarse.ICu-SolnBlk_Coarse.ICl+1)/2+SolnBlk_Coarse.ICl;
	j_coarse = (j_fine-SolnBlk_Original_NE.JCl)/2 + (SolnBlk_Coarse.JCu-SolnBlk_Coarse.JCl+1)/2+SolnBlk_Coarse.JCl;

	// Restrict the solution information from the fine cells to the
	// coarse cell using an area-weighted average of the fine cell
	// solution information.
	SolnBlk_Coarse.U[i_coarse][j_coarse] = (SolnBlk_Original_NE.Grid.Cell[i_fine  ][j_fine  ].A*
						SolnBlk_Original_NE.U[i_fine  ][j_fine  ] +
						SolnBlk_Original_NE.Grid.Cell[i_fine+1][j_fine  ].A*
						SolnBlk_Original_NE.U[i_fine+1][j_fine  ] + 
						SolnBlk_Original_NE.Grid.Cell[i_fine  ][j_fine+1].A*
						SolnBlk_Original_NE.U[i_fine  ][j_fine+1] +
						SolnBlk_Original_NE.Grid.Cell[i_fine+1][j_fine+1].A*
						SolnBlk_Original_NE.U[i_fine+1][j_fine+1]) /
                                               (SolnBlk_Original_NE.Grid.Cell[i_fine  ][j_fine  ].A +
						SolnBlk_Original_NE.Grid.Cell[i_fine+1][j_fine  ].A +
						SolnBlk_Original_NE.Grid.Cell[i_fine  ][j_fine+1].A +
						SolnBlk_Original_NE.Grid.Cell[i_fine+1][j_fine+1].A);
	                                       //SolnBlk_Coarse.Grid.Cell[i_coarse][j_coarse].A;
	SolnBlk_Coarse.W[i_coarse][j_coarse] = W(SolnBlk_Coarse.U[i_coarse][j_coarse]);

      }
    }

    for (int j_fine = SolnBlk_Original_NE.JCl-SolnBlk_Original_NE.Nghost; 
	 j_fine <= SolnBlk_Original_NE.JCu+SolnBlk_Original_NE.Nghost; j_fine += 2) {
      j_coarse = (j_fine-SolnBlk_Original_NE.JCl)/2 + (SolnBlk_Coarse.JCu-SolnBlk_Coarse.JCl+1)/2+SolnBlk_Coarse.JCl;
      SolnBlk_Coarse.WoE[j_coarse] = SolnBlk_Original_NE.WoE[j_fine];
      if (j_fine == SolnBlk_Original_NE.JCu+SolnBlk_Original_NE.Nghost) {
	SolnBlk_Coarse.WoE[j_coarse+1] = SolnBlk_Original_NE.WoE[j_fine];
      }
    }

    for (int i_fine = SolnBlk_Original_NE.ICl-SolnBlk_Original_NE.Nghost; 
	 i_fine <= SolnBlk_Original_NE.ICu+SolnBlk_Original_NE.Nghost; i_fine += 2) {
      i_coarse = (i_fine-SolnBlk_Original_NE.ICl)/2 + (SolnBlk_Coarse.ICu-SolnBlk_Coarse.ICl+1)/2+SolnBlk_Coarse.ICl;
      SolnBlk_Coarse.WoN[i_coarse] = SolnBlk_Original_NE.WoN[i_fine];
      if (i_fine == SolnBlk_Original_NE.ICu+SolnBlk_Original_NE.Nghost) {
	SolnBlk_Coarse.WoN[i_coarse+1] = SolnBlk_Original_NE.WoN[i_fine];
      }
    }

  }

  // Restriction of solution block was successful.
  return 0;

}

/**********************************************************************
 * Routine: ICs                                                       *
 *                                                                    *
 * Assigns initial conditions and data to the solution variables of   *
 * the specified quadrilateral solution block.                        *
 *                                                                    *
 **********************************************************************/
void ICs(NavierStokes2D_Quad_Block &SolnBlk,
         NavierStokes2D_Input_Parameters &IP,
         NavierStokes2D_pState *Wo) {

  NavierStokes2D_pState Wl, Wr;
  double eta, f, fp, fpp;
  double utau, tauw;
  NavierStokes2D_pState NavierStokes2D_W_STDATM; NavierStokes2D_W_STDATM.Standard_Atmosphere();
  ifstream indata, comparedata; 
  // Assign the initial data for the IVP of interest.
  switch(IP.i_ICs) {
  case IC_CONSTANT :
  case IC_UNIFORM :
    // Set the solution state to the initial state Wo[0].
    for (int j = SolnBlk.JCl-SolnBlk.Nghost; j <= SolnBlk.JCu+SolnBlk.Nghost; j++) {
      for (int i = SolnBlk.ICl-SolnBlk.Nghost; i <= SolnBlk.ICu+SolnBlk.Nghost; i++) {
	SolnBlk.W[i][j] = Wo[0];
 	if (SolnBlk.Grid.Cell[i][j].Xc.x<ZERO){
 	  SolnBlk.W[i][j].p = SolnBlk.W[i][j].p*1.3;
	  SolnBlk.W[i][j].v.y = IP.Mach_Number* SolnBlk.W[i][j].a();
	   SolnBlk.W[i][j].v.x = ZERO; 
	}else {
	  //SolnBlk.W[i][j].p = SolnBlk.W[i][j].p*(1.4-0.4/0.005*SolnBlk.Grid.Cell[i][j].Xc.x);
	  SolnBlk.W[i][j].v.x = ZERO;
	}
	if (SolnBlk.Flow_Type == FLOWTYPE_TURBULENT_RANS_K_OMEGA) {
	  //SolnBlk.W[i][j].k = max(ONE,0.00005*(sqr(SolnBlk.W[i][j].v.x) + sqr(SolnBlk.W[i][j].v.y)));
	  SolnBlk.W[i][j].k = 0.00005*(sqr(SolnBlk.W[i][j].v.x) + sqr(SolnBlk.W[i][j].v.y));
	  SolnBlk.W[i][j].omega = SolnBlk.W[i][j].k/0.000001;
	}
	if (SolnBlk.Variable_Prandtl == ON) {
	  //SolnBlk.W[i][j].ke = 100.0;
	  //SolnBlk.W[i][j].ee = 0.1;
	  SolnBlk.W[i][j].ke = sqr(0.0001*SolnBlk.W[i][j].p/(IP.Wo.gm1*SolnBlk.W[i][j].rho));
	  //For the inital guess of 'ee' we fix the Prandtl number (0.7 in the present case) and calculate all 'ee' accordingly.
	  double CmuFmu = SolnBlk.W[i][j].muT()*SolnBlk.W[i][j].epsilon()/SolnBlk.W[i][j].rho/max(sqr(SolnBlk.W[i][j].k),sqr(TOLER));
	  double Coeff = sqr(0.7*0.14*SolnBlk.W[i][j].f_lambda(SolnBlk.Wall[i][j].ywall)/max(CmuFmu,sqr(TOLER)));
	  SolnBlk.W[i][j].ee = Coeff*SolnBlk.W[i][j].epsilon()* SolnBlk.W[i][j].ke/ max(SolnBlk.W[i][j].k,TOLER);
	}
	SolnBlk.U[i][j] = U(SolnBlk.W[i][j]);
      }
    }
    if (1 == 2 && IP.i_Grid == GRID_ROCKET_MOTOR &&
	SolnBlk.Flow_Type == FLOWTYPE_TURBULENT_RANS_K_OMEGA) {
      double length = IP.Chamber_Length+IP.Nozzle_Length;
      Vector2D twall, that, v;
      for (int j = SolnBlk.JCl-SolnBlk.Nghost; j <= SolnBlk.JCu+SolnBlk.Nghost; j++) {
	for (int i = SolnBlk.ICl-SolnBlk.Nghost; i <= SolnBlk.ICu+SolnBlk.Nghost; i++) {
	  SolnBlk.W[i][j] = TurbulentPipeFlow(IP.Wo,
					      SolnBlk.Grid.Cell[i][j].Xc,
					      IP.dp,
					      length,
					      SolnBlk.Wall[i][j].Xwall.y,
					      IP.Reynolds_Number);
    //    //SolnBlk.W[i][j].p = IP.Wo.p - (2.0/3.0)*SolnBlk.W[i][j].rho*SolnBlk.W[i][j].k;
    //    SolnBlk.W[i][j].p += (2.0/3.0)*SolnBlk.W[i][j].rho*SolnBlk.W[i][j].k;
    //    SolnBlk.W[i][j].omega *= 100.0;
	  if (SolnBlk.Grid.Cell[2][2].Xc.x > 0.0) {
	    v = SolnBlk.W[i][j].v;
	    twall.x = -SolnBlk.Wall[i][j].nwall.y;
	    twall.y =  SolnBlk.Wall[i][j].nwall.x;
	    that = twall*(SolnBlk.Grid.Cell[i][j].Xc.y/SolnBlk.Wall[i][j].Xwall.y) +
	      ihat*(ONE - SolnBlk.Grid.Cell[i][j].Xc.y/SolnBlk.Wall[i][j].Xwall.y);
	    SolnBlk.W[i][j].v.x = v.x*that.x + v.y*that.y;
	    SolnBlk.W[i][j].v.y = v.x*that.y - v.y*that.x;
	    SolnBlk.U[i][j] = U(SolnBlk.W[i][j]);
	  }
	}
      }
    }
    break;
  case IC_SOD_XDIR :
    // Set the initial data for the Sod IVP in x-direction.
    Wl = NavierStokes2D_W_STDATM;
    Wr = NavierStokes2D_pState(DENSITY_STDATM/EIGHT,ZERO,ZERO,PRESSURE_STDATM/TEN);
    for (int j = SolnBlk.JCl-SolnBlk.Nghost; j <= SolnBlk.JCu+SolnBlk.Nghost; j++) {
      for (int i = SolnBlk.ICl-SolnBlk.Nghost; i <= SolnBlk.ICu+SolnBlk.Nghost; i++) {
	if (SolnBlk.Grid.Cell[i][j].Xc.x <= ZERO) {
	  SolnBlk.W[i][j] = Wl;
	} else {
	  SolnBlk.W[i][j] = Wr;
	}
	SolnBlk.U[i][j] = U(SolnBlk.W[i][j]);
      }
    }
    break;
  case IC_SOD_YDIR :
    // Set the initial data for the Sod IVP in y-direction.
    Wl = NavierStokes2D_W_STDATM;
    Wr = NavierStokes2D_pState(DENSITY_STDATM/EIGHT,ZERO,ZERO,PRESSURE_STDATM/TEN);
    for (int j = SolnBlk.JCl-SolnBlk.Nghost; j <= SolnBlk.JCu+SolnBlk.Nghost; j++) {
      for (int i = SolnBlk.ICl-SolnBlk.Nghost; i <= SolnBlk.ICu+SolnBlk.Nghost; i++) {
	if (SolnBlk.Grid.Cell[i][j].Xc.y <= ZERO) {
	  SolnBlk.W[i][j] = Wl;
	} else {
	  SolnBlk.W[i][j] = Wr;
	}
	SolnBlk.U[i][j] = U(SolnBlk.W[i][j]);
      }
    }
    break;
  case IC_GROTH_XDIR :
    // Set the initial data for the Groth IVP in x-direction.
    Wl = NavierStokes2D_pState(4.696,ZERO,ZERO,404.4e03);
    Wr = NavierStokes2D_pState(1.408,ZERO,ZERO,101.1e03);
    for (int j = SolnBlk.JCl-SolnBlk.Nghost; j <= SolnBlk.JCu+SolnBlk.Nghost; j++) {
      for (int i = SolnBlk.ICl-SolnBlk.Nghost; i <= SolnBlk.ICu+SolnBlk.Nghost; i++) {
	if (SolnBlk.Grid.Cell[i][j].Xc.x <= ZERO) {
	  SolnBlk.W[i][j] = Wl;
	} else {
	  SolnBlk.W[i][j] = Wr;
	}
	SolnBlk.U[i][j] = U(SolnBlk.W[i][j]);
      }
    }
    break;
  case IC_GROTH_YDIR :
    // Set the initial data for the Groth IVP in y-direction.
    Wl = NavierStokes2D_pState(4.696,ZERO,ZERO,404.4e03);
    Wr = NavierStokes2D_pState(1.408,ZERO,ZERO,101.1e03);
    for (int j = SolnBlk.JCl-SolnBlk.Nghost; j <= SolnBlk.JCu+SolnBlk.Nghost; j++) {
      for (int i = SolnBlk.ICl-SolnBlk.Nghost; i <= SolnBlk.ICu+SolnBlk.Nghost; i++) {
	if (SolnBlk.Grid.Cell[i][j].Xc.y <= ZERO) {
	  SolnBlk.W[i][j] = Wl;
	} else {
	  SolnBlk.W[i][j] = Wr;
	}
	SolnBlk.U[i][j] = U(SolnBlk.W[i][j]);
      }
    }
    break;
  case IC_EINFELDT_XDIR :
    // Set the initial data for the Einfeldt IVP in x-direction.
    Wl = NavierStokes2D_pState(ONE,-TWO,ZERO,FOUR/TEN);
    Wr = NavierStokes2D_pState(ONE, TWO,ZERO,FOUR/TEN);
    for (int j = SolnBlk.JCl-SolnBlk.Nghost; j <= SolnBlk.JCu+SolnBlk.Nghost; j++) {
      for (int i = SolnBlk.ICl-SolnBlk.Nghost; i <= SolnBlk.ICu+SolnBlk.Nghost; i++) {
	if (SolnBlk.Grid.Cell[i][j].Xc.x <= ZERO) {
	  SolnBlk.W[i][j] = Wl;
	} else {
	  SolnBlk.W[i][j] = Wr;
	}
	SolnBlk.U[i][j] = U(SolnBlk.W[i][j]);
      }
    }
    break;
  case IC_EINFELDT_YDIR :
    // Set the initial data for the Einfeld IVP in y-direction.
    Wl = NavierStokes2D_pState(ONE,ZERO,-TWO,FOUR/TEN);
    Wr = NavierStokes2D_pState(ONE,ZERO, TWO,FOUR/TEN);
    for (int j = SolnBlk.JCl-SolnBlk.Nghost; j <= SolnBlk.JCu+SolnBlk.Nghost; j++) {
      for (int i = SolnBlk.ICl-SolnBlk.Nghost; i <= SolnBlk.ICu+SolnBlk.Nghost; i++) {
	if (SolnBlk.Grid.Cell[i][j].Xc.y <= ZERO) {
	  SolnBlk.W[i][j] = Wl;
	} else {
	  SolnBlk.W[i][j] = Wr;
	}
	SolnBlk.U[i][j] = U(SolnBlk.W[i][j]);
      }
    }
    break;
  case IC_SHOCK_WAVE_XDIR :
    // Set the initial data for moving shock wave propagating in the x-direction.
    Wl = NavierStokes2D_pState(2.281, 164.83, ZERO, 201.17e03);
    Wr = NavierStokes2D_pState(1.408, ZERO, ZERO, 101.1e03);
    //Wl = NavierStokes2D_pState(1.3189655, 57.543675, ZERO, 126.65625e03);
    //Wr = NavierStokes2D_pState(1.125, ZERO, ZERO, 101.325e03);
    //Wl = NavierStokes2D_pState(1.828125, 186.120997, ZERO, 202650.5);
    //Wr = NavierStokes2D_pState(1.125, ZERO, ZERO, 101325.0);
    //Wl = Wo[0];
    //Wr = NavierStokes2D_W_STDATM;
    for (int j = SolnBlk.JCl-SolnBlk.Nghost; j <= SolnBlk.JCu+SolnBlk.Nghost; j++) {
      for (int i = SolnBlk.ICl-SolnBlk.Nghost; i <= SolnBlk.ICu+SolnBlk.Nghost; i++) {
	if (SolnBlk.Grid.Cell[i][j].Xc.x <= IP.Wave_Position.x) {
	  SolnBlk.W[i][j] = Wl;
	} else {
	  SolnBlk.W[i][j] = Wr;
	}
	SolnBlk.U[i][j] = U(SolnBlk.W[i][j]);
      }
    }
    break;
  case IC_SHOCK_WAVE_YDIR :
    // Set the initial data for moving shock wave propagating in the y-direction.
    Wl = NavierStokes2D_pState(2.281,ZERO,164.83,201.17e03);
    Wr = NavierStokes2D_pState(1.408,ZERO,ZERO,  101.10e03);
   for (int j = SolnBlk.JCl-SolnBlk.Nghost; j <= SolnBlk.JCu+SolnBlk.Nghost; j++) {
      for (int i = SolnBlk.ICl-SolnBlk.Nghost; i <= SolnBlk.ICu+SolnBlk.Nghost; i++) {
	if (SolnBlk.Grid.Cell[i][j].Xc.y <= IP.Wave_Position.y) {
	  SolnBlk.W[i][j] = Wl;
	} else {
	  SolnBlk.W[i][j] = Wr;
	}
	SolnBlk.U[i][j] = U(SolnBlk.W[i][j]);
      }
    }
    break;
  case IC_CONTACT_SURFACE_XDIR :
    // Set the initial data for moving contact surface propagating in the x-direction.
    Wl = NavierStokes2D_pState(1.045,200.0,ZERO,PRESSURE_STDATM);
    Wr = NavierStokes2D_pState(3.483,200.0,ZERO,PRESSURE_STDATM);
    for (int j = SolnBlk.JCl-SolnBlk.Nghost; j <= SolnBlk.JCu+SolnBlk.Nghost; j++) {
      for (int i = SolnBlk.ICl-SolnBlk.Nghost; i <= SolnBlk.ICu+SolnBlk.Nghost; i++) {
	if (SolnBlk.Grid.Cell[i][j].Xc.x <= ZERO) {
	  SolnBlk.W[i][j] = Wl;
	} else {
	  SolnBlk.W[i][j] = Wr;
	}
	SolnBlk.U[i][j] = U(SolnBlk.W[i][j]);
      }
    }
    break;
  case IC_CONTACT_SURFACE_YDIR :
    // Set the initial data for moving contact surface propagating in the y-direction.
    Wl = NavierStokes2D_pState(1.045,ZERO,200.0,PRESSURE_STDATM);
    Wr = NavierStokes2D_pState(3.483,ZERO,200.0,PRESSURE_STDATM);
    for (int j = SolnBlk.JCl-SolnBlk.Nghost; j <= SolnBlk.JCu+SolnBlk.Nghost; j++) {
      for (int i = SolnBlk.ICl-SolnBlk.Nghost; i <= SolnBlk.ICu+SolnBlk.Nghost; i++) {
	if (SolnBlk.Grid.Cell[i][j].Xc.y <= ZERO) {
	  SolnBlk.W[i][j] = Wl;
	} else {
	  SolnBlk.W[i][j] = Wr;
	}
	SolnBlk.U[i][j] = U(SolnBlk.W[i][j]);
      }
    }
    break;
  case IC_RAREFACTION_WAVE_XDIR :
    // Set the initial data for moving rarefaction wave propagating in the x-direction.
    Wl = NavierStokes2D_pState(1.598,-383.64,ZERO, 91.88e03);
    Wr = NavierStokes2D_pState(2.787,-216.97,ZERO,200.09e03);
    for (int j = SolnBlk.JCl-SolnBlk.Nghost; j <= SolnBlk.JCu+SolnBlk.Nghost; j++) {
      for (int i = SolnBlk.ICl-SolnBlk.Nghost; i <= SolnBlk.ICu+SolnBlk.Nghost; i++) {
	if (SolnBlk.Grid.Cell[i][j].Xc.x <= ZERO) {
	  SolnBlk.W[i][j] = Wl;
	} else {
	  SolnBlk.W[i][j] = Wr;
	}
	SolnBlk.U[i][j] = U(SolnBlk.W[i][j]);
      }
    }
    break;
  case IC_RAREFACTION_WAVE_YDIR :
    // Set the initial data for moving rarefaction wave propagating in the y-direction.
    Wl = NavierStokes2D_pState(1.598,ZERO,-383.64,91.88e03);
    Wr = NavierStokes2D_pState(2.787,ZERO,-216.97,200.0e03);
    for (int j = SolnBlk.JCl-SolnBlk.Nghost; j <= SolnBlk.JCu+SolnBlk.Nghost; j++) {
      for (int i = SolnBlk.ICl-SolnBlk.Nghost; i <= SolnBlk.ICu+SolnBlk.Nghost; i++) {
	if (SolnBlk.Grid.Cell[i][j].Xc.y <= ZERO) {
	  SolnBlk.W[i][j] = Wl;
	} else {
	  SolnBlk.W[i][j] = Wr;
	}
	SolnBlk.U[i][j] = U(SolnBlk.W[i][j]);
      }
    }
    break;
  case IC_SHOCK_BOX :
    // Set the initial data for the Aki shock-box IVP.
    Wl = NavierStokes2D_W_STDATM;
    Wr = NavierStokes2D_pState(DENSITY_STDATM*FOUR,ZERO,ZERO,PRESSURE_STDATM*FOUR);
    for (int j = SolnBlk.JCl-SolnBlk.Nghost; j <= SolnBlk.JCu+SolnBlk.Nghost; j++) {
      for (int i = SolnBlk.ICl-SolnBlk.Nghost; i <= SolnBlk.ICu+SolnBlk.Nghost; i++) {
	if (SolnBlk.Grid.Cell[i][j].Xc.x <= ZERO &&
	    SolnBlk.Grid.Cell[i][j].Xc.y <= ZERO) {
	  SolnBlk.W[i][j] = Wl;
	} else {
	  SolnBlk.W[i][j] = Wr;
	}
	SolnBlk.U[i][j] = U(SolnBlk.W[i][j]);
      }
    }
    break;
  case IC_HIGH_PRESSURE_RESERVOIR :
    // Set high-pressure reservoir initial data.
    Wl = NavierStokes2D_pState(HUNDRED*DENSITY_STDATM,ZERO,ZERO,HUNDRED*PRESSURE_STDATM);
    Wr = NavierStokes2D_W_STDATM;
    for (int j = SolnBlk.JCl-SolnBlk.Nghost; j <= SolnBlk.JCu+SolnBlk.Nghost; j++) {
      for (int i = SolnBlk.ICl-SolnBlk.Nghost; i <= SolnBlk.ICu+SolnBlk.Nghost; i++) {
	if (SolnBlk.Grid.Cell[i][j].Xc.x <= IP.Wave_Position.x) {
	  SolnBlk.W[i][j] = Wl;
	} else {
	  SolnBlk.W[i][j] = Wr;
	}
	SolnBlk.U[i][j] = U(SolnBlk.W[i][j]);
      }
    }
    break;
  case IC_LOW_PRESSURE_RESERVOIR :
    // Set low-pressure reservoir initial data.
    Wl = NavierStokes2D_W_STDATM;
    Wr = NavierStokes2D_pState(Wo[0].rho,ZERO,ZERO,Wo[0].p);
    for (int j = SolnBlk.JCl-SolnBlk.Nghost; j <= SolnBlk.JCu+SolnBlk.Nghost; j++) {
      for (int i = SolnBlk.ICl-SolnBlk.Nghost; i <= SolnBlk.ICu+SolnBlk.Nghost; i++) {
	if (SolnBlk.Grid.Cell[i][j].Xc.x <= IP.Wave_Position.x) {
	  SolnBlk.W[i][j] = Wl;
	} else {
	  SolnBlk.W[i][j] = Wr;
	}
  	SolnBlk.U[i][j] = U(SolnBlk.W[i][j]);
      }
    }
    break;
  case IC_RIEMANN :
  case IC_RIEMANN_XDIR :
    // Set the Riemann initial data (2-state initial data, left to right).
    Wr = Wo[1];
    Wl = Wo[0];
    for (int j = SolnBlk.JCl-SolnBlk.Nghost; j <= SolnBlk.JCu+SolnBlk.Nghost; j++) {
      for (int i = SolnBlk.ICl-SolnBlk.Nghost; i <= SolnBlk.ICu+SolnBlk.Nghost; i++) {
	if (SolnBlk.Grid.Cell[i][j].Xc.x <= ZERO) {
	  SolnBlk.W[i][j] = Wl;
	} else {
	  SolnBlk.W[i][j] = Wr;
	}
	SolnBlk.U[i][j] = U(SolnBlk.W[i][j]);
      }
    }
    break;
  case IC_RIEMANN_YDIR :
    // Set the Riemann initial data (2-state initial data, top to bottom).
    Wr = Wo[1];
    Wl = Wo[0];
    for (int j = SolnBlk.JCl-SolnBlk.Nghost; j <= SolnBlk.JCu+SolnBlk.Nghost; j++) {
      for (int i = SolnBlk.ICl-SolnBlk.Nghost; i <= SolnBlk.ICu+SolnBlk.Nghost; i++) {
	if (SolnBlk.Grid.Cell[i][j].Xc.y <= ZERO) {
	  SolnBlk.W[i][j] = Wl;
	} else {
	  SolnBlk.W[i][j] = Wr;
	}
	SolnBlk.U[i][j] = U(SolnBlk.W[i][j]);
      }
    }
    break;
  case IC_CYLINDRICAL_EXPLOSION :
    // Set the initial data for the cylindrical explosion IVP (see Toro's book).
    Wl = NavierStokes2D_W_STDATM;
    Wr = NavierStokes2D_pState(EIGHT*DENSITY_STDATM,ZERO,ZERO,TEN*PRESSURE_STDATM);
    for (int j = SolnBlk.JCl-SolnBlk.Nghost; j <= SolnBlk.JCu+SolnBlk.Nghost; j++) {
      for (int i = SolnBlk.ICl-SolnBlk.Nghost; i <= SolnBlk.ICu+SolnBlk.Nghost; i++) {
	if (abs(SolnBlk.Grid.Cell[i][j].Xc) < 0.25) {
	  SolnBlk.W[i][j] = Wr;
	} else {
	  SolnBlk.W[i][j] = Wl;
	}
	SolnBlk.U[i][j] = U(SolnBlk.W[i][j]);
      }
    }
    break;
  case IC_CYLINDRICAL_IMPLOSION :
    // Set the initial data for the cylindrical implosion IVP (see Toro's book).
    Wl = NavierStokes2D_pState(EIGHT*DENSITY_STDATM,ZERO,ZERO,TEN*PRESSURE_STDATM);
    Wr = NavierStokes2D_W_STDATM;
    for (int j = SolnBlk.JCl-SolnBlk.Nghost; j <= SolnBlk.JCu+SolnBlk.Nghost; j++) {
      for (int i = SolnBlk.ICl-SolnBlk.Nghost; i <= SolnBlk.ICu+SolnBlk.Nghost; i++) {
	if (abs(SolnBlk.Grid.Cell[i][j].Xc) < 0.25) {
	  SolnBlk.W[i][j] = Wr;
	} else {
	  SolnBlk.W[i][j] = Wl;
	}
	SolnBlk.U[i][j] = U(SolnBlk.W[i][j]);
      }
    }
    break;
  case IC_SQUARE_WAVE_XDIR :
    // Set the initial data for a moving square wave in the x-direction.
    for (int j = SolnBlk.JCl-SolnBlk.Nghost; j <= SolnBlk.JCu+SolnBlk.Nghost; j++) {
      for (int i = SolnBlk.ICl-SolnBlk.Nghost; i <= SolnBlk.ICu+SolnBlk.Nghost; i++) {
	if (SolnBlk.Grid.Cell[i][j].Xc.x <= IP.Wave_Position.x+HALF*IP.Wave_Width &&
	    SolnBlk.Grid.Cell[i][j].Xc.x >= IP.Wave_Position.x-HALF*IP.Wave_Width) {
	  SolnBlk.W[i][j] = Wo[0];
	} else {
	  SolnBlk.W[i][j] = NavierStokes2D_W_STDATM;
	}
	SolnBlk.W[i][j].rho = Wo[0].rho;
	SolnBlk.W[i][j].v.x = Wo[0].v.x;
	SolnBlk.W[i][j].v.y = Wo[0].v.y;
	SolnBlk.W[i][j].p = Wo[0].p;
	SolnBlk.U[i][j] = U(SolnBlk.W[i][j]);
      }
    }
    break;
  case IC_SQUARE_WAVE_YDIR :
    // Set the initial data for a moving square wave in the y-direction.
    for (int j = SolnBlk.JCl-SolnBlk.Nghost; j <= SolnBlk.JCu+SolnBlk.Nghost; j++) {
      for (int i = SolnBlk.ICl-SolnBlk.Nghost; i <= SolnBlk.ICu+SolnBlk.Nghost; i++) {
	if (SolnBlk.Grid.Cell[i][j].Xc.y <= IP.Wave_Position.y+HALF*IP.Wave_Width &&
	    SolnBlk.Grid.Cell[i][j].Xc.y >= IP.Wave_Position.y-HALF*IP.Wave_Width) {
	  SolnBlk.W[i][j] = Wo[0];
	} else {
	  SolnBlk.W[i][j] = NavierStokes2D_W_STDATM;
	}
	SolnBlk.W[i][j].rho = Wo[0].rho;
	SolnBlk.W[i][j].v.x = Wo[0].v.x;
	SolnBlk.W[i][j].v.y = Wo[0].v.y;
	SolnBlk.W[i][j].p = Wo[0].p;
	SolnBlk.U[i][j] = U(SolnBlk.W[i][j]);
      }
    }
    break;
  case IC_SINE_WAVE_XDIR :
    // Set the initial data for a moving sine wave in the x-direction.
    for (int j = SolnBlk.JCl-SolnBlk.Nghost; j <= SolnBlk.JCu+SolnBlk.Nghost; j++) {
      for (int i = SolnBlk.ICl-SolnBlk.Nghost; i <= SolnBlk.ICu+SolnBlk.Nghost; i++) {
	if (SolnBlk.Grid.Cell[i][j].Xc.x <= IP.Wave_Position.x+HALF*IP.Wave_Width &&
	    SolnBlk.Grid.Cell[i][j].Xc.x >= IP.Wave_Position.x-HALF*IP.Wave_Width) {
 	  SolnBlk.W[i][j] = Wo[0]*sin(PI*(SolnBlk.Grid.Cell[i][j].Xc.x+HALF*IP.Wave_Width)/IP.Wave_Width);
	} else {
	  SolnBlk.W[i][j] = NavierStokes2D_W_STDATM;
	}
	SolnBlk.W[i][j].rho = Wo[0].rho;
	SolnBlk.W[i][j].v.x = Wo[0].v.x;
	SolnBlk.W[i][j].v.y = Wo[0].v.y;
	SolnBlk.W[i][j].p = Wo[0].p;
	SolnBlk.U[i][j] = U(SolnBlk.W[i][j]);
      }
    }
    break;
  case IC_SINE_WAVE_YDIR :
    // Set the initial data for a moving sine wave in the y-direction.
    for (int j = SolnBlk.JCl-SolnBlk.Nghost; j <= SolnBlk.JCu+SolnBlk.Nghost; j++) {
      for (int i = SolnBlk.ICl-SolnBlk.Nghost; i <= SolnBlk.ICu+SolnBlk.Nghost; i++) {
	if (SolnBlk.Grid.Cell[i][j].Xc.y <= IP.Wave_Position.y+HALF*IP.Wave_Width &&
	    SolnBlk.Grid.Cell[i][j].Xc.y >= IP.Wave_Position.y-HALF*IP.Wave_Width) {
	  SolnBlk.W[i][j] = Wo[0]*sin(PI*(SolnBlk.Grid.Cell[i][j].Xc.y+HALF*IP.Wave_Width)/IP.Wave_Width);
	} else {
	  SolnBlk.W[i][j] = NavierStokes2D_W_STDATM;
	}
	SolnBlk.U[i][j] = U(SolnBlk.W[i][j]);
      }
    }
    break;
  case IC_SINE2_WAVE_XDIR :
    // Set the initial data for a moving sine-squared wave in the x-direction.
    for (int j = SolnBlk.JCl-SolnBlk.Nghost; j <= SolnBlk.JCu+SolnBlk.Nghost; j++) {
      for (int i = SolnBlk.ICl-SolnBlk.Nghost; i <= SolnBlk.ICu+SolnBlk.Nghost; i++) {
	if (SolnBlk.Grid.Cell[i][j].Xc.x <= IP.Wave_Position.x+HALF*IP.Wave_Width &&
	    SolnBlk.Grid.Cell[i][j].Xc.x >= IP.Wave_Position.x-HALF*IP.Wave_Width) {
	  SolnBlk.W[i][j] = Wo[0]*sqr(sin(PI*(SolnBlk.Grid.Cell[i][j].Xc.x+HALF*IP.Wave_Width)/IP.Wave_Width));
	} else {
	  SolnBlk.W[i][j] = NavierStokes2D_W_STDATM;
	}
	SolnBlk.U[i][j] = U(SolnBlk.W[i][j]);
      }
    }
    break;
  case IC_SINE2_WAVE_YDIR :
    // Set the initial data for a moving sine-squared wave in the y-direction.
    for (int j = SolnBlk.JCl-SolnBlk.Nghost; j <= SolnBlk.JCu+SolnBlk.Nghost; j++) {
      for (int i = SolnBlk.ICl-SolnBlk.Nghost; i <= SolnBlk.ICu+SolnBlk.Nghost; i++) {
	if (SolnBlk.Grid.Cell[i][j].Xc.y <= IP.Wave_Position.y+HALF*IP.Wave_Width &&
	    SolnBlk.Grid.Cell[i][j].Xc.y >= IP.Wave_Position.y-HALF*IP.Wave_Width) {
	  SolnBlk.W[i][j] = Wo[0]*sqr(sin(PI*(SolnBlk.Grid.Cell[i][j].Xc.y+HALF*IP.Wave_Width)/IP.Wave_Width));
	} else {
	  SolnBlk.W[i][j] = NavierStokes2D_W_STDATM;
	}
	SolnBlk.U[i][j] = U(SolnBlk.W[i][j]);
      }
    }
    break;
  case IC_COMPRESSION_XDIR :
    // Set the initial data for a one-dimensional compression problem in the x-direction.
    Wl = NavierStokes2D_pState(ONE, 200.0,ZERO,PRESSURE_STDATM);
    Wr = NavierStokes2D_pState(ONE,-200.0,ZERO,PRESSURE_STDATM);
    for (int j = SolnBlk.JCl-SolnBlk.Nghost; j <= SolnBlk.JCu+SolnBlk.Nghost; j++) {
      for (int i = SolnBlk.ICl-SolnBlk.Nghost; i <= SolnBlk.ICu+SolnBlk.Nghost; i++) {
	if (SolnBlk.Grid.Cell[i][j].Xc.x <= ZERO) {
	  SolnBlk.W[i][j] = Wl;
	} else {
	  SolnBlk.W[i][j] = Wr;
	}
	SolnBlk.U[i][j] = U(SolnBlk.W[i][j]);
      }
    }
    break;
  case IC_COMPRESSION_YDIR :
    // Set the initial data for a one-dimensional compression problem in the y-direction.
    Wl = NavierStokes2D_pState(ONE,ZERO, 200.0,PRESSURE_STDATM);
    Wr = NavierStokes2D_pState(ONE,ZERO,-200.0,PRESSURE_STDATM);
    for (int j = SolnBlk.JCl-SolnBlk.Nghost; j <= SolnBlk.JCu+SolnBlk.Nghost; j++) {
      for (int i = SolnBlk.ICl-SolnBlk.Nghost; i <= SolnBlk.ICu+SolnBlk.Nghost; i++) {
	if (SolnBlk.Grid.Cell[i][j].Xc.y <= ZERO) {
	  SolnBlk.W[i][j] = Wl;
	} else {
	  SolnBlk.W[i][j] = Wr;
	}
	SolnBlk.U[i][j] = U(SolnBlk.W[i][j]);
      }
    }
    break;
  case IC_EXPANSION_XDIR :
    // Set the initial data for a one-dimensional expansion problem in the x-direction.
    Wl = NavierStokes2D_pState(ONE,-200.0,ZERO,PRESSURE_STDATM);
    Wr = NavierStokes2D_pState(ONE, 200.0,ZERO,PRESSURE_STDATM);
    for (int j = SolnBlk.JCl-SolnBlk.Nghost; j <= SolnBlk.JCu+SolnBlk.Nghost; j++) {
      for (int i = SolnBlk.ICl-SolnBlk.Nghost; i <= SolnBlk.ICu+SolnBlk.Nghost; i++) {
	if (SolnBlk.Grid.Cell[i][j].Xc.x <= ZERO) {
	  SolnBlk.W[i][j] = Wl;
	} else {
	  SolnBlk.W[i][j] = Wr;
	}
	SolnBlk.U[i][j] = U(SolnBlk.W[i][j]);
      }
    }
    break;
  case IC_EXPANSION_YDIR :
    // Set the initial data for a one-dimensional expansion problem in the y-direction.
    Wl = NavierStokes2D_pState(ONE,ZERO,-200.0,PRESSURE_STDATM);
    Wr = NavierStokes2D_pState(ONE,ZERO, 200.0,PRESSURE_STDATM);
    for (int j = SolnBlk.JCl-SolnBlk.Nghost; j <= SolnBlk.JCu+SolnBlk.Nghost; j++) {
      for (int i = SolnBlk.ICl-SolnBlk.Nghost; i <= SolnBlk.ICu+SolnBlk.Nghost; i++) {
	if (SolnBlk.Grid.Cell[i][j].Xc.y <= ZERO) {
	  SolnBlk.W[i][j] = Wl;
	} else {
	  SolnBlk.W[i][j] = Wr;
	}
	SolnBlk.U[i][j] = U(SolnBlk.W[i][j]);
      }
    }
    break;
  case IC_RINGLEB_FLOW :
    // Set the initial data to the analytic solution to Ringleb's problem.
    for (int j = SolnBlk.JCl-SolnBlk.Nghost; j <= SolnBlk.JCu+SolnBlk.Nghost; j++) {
      for (int i = SolnBlk.ICl-SolnBlk.Nghost; i <= SolnBlk.ICu+SolnBlk.Nghost; i++) {
	//SolnBlk.W[i][j] = RinglebFlow(SolnBlk.W[i][j],SolnBlk.Grid.Cell[i][j].Xc);
	SolnBlk.W[i][j] = RinglebFlowAverageState(SolnBlk.W[i][j],
						  SolnBlk.Grid.nodeSW(i,j).X,
						  SolnBlk.Grid.nodeSE(i,j).X,
						  SolnBlk.Grid.nodeNE(i,j).X,
						  SolnBlk.Grid.nodeNW(i,j).X);
	SolnBlk.U[i][j] = U(SolnBlk.W[i][j]);
      }
    }
    break;
  case IC_VISCOUS_CHANNEL_FLOW :
    // Set the initial data to the analytic solution of a laminar 
    // channel flow in the x-direction.
    // -635.00  -423.33  -211.67  0.00  211.67  423.33  635.00
    for (int j = SolnBlk.JCl-SolnBlk.Nghost; j <= SolnBlk.JCu+SolnBlk.Nghost; j++) {
      for (int i = SolnBlk.ICl-SolnBlk.Nghost; i <= SolnBlk.ICu+SolnBlk.Nghost; i++) {
	SolnBlk.W[i][j] = ViscousChannelFlow(SolnBlk.W[i][j],SolnBlk.Grid.Cell[i][j].Xc,
					     SolnBlk.Vwall,IP.dp,IP.Box_Width,IP.Box_Height);
	SolnBlk.U[i][j] = U(SolnBlk.W[i][j]);
      }
    }
    break;
  case IC_VISCOUS_PIPE_FLOW :
    if (SolnBlk.Flow_Type == FLOWTYPE_INVISCID) {
      for (int j = SolnBlk.JCl-SolnBlk.Nghost; j <= SolnBlk.JCu+SolnBlk.Nghost; j++) {
	for (int i = SolnBlk.ICl-SolnBlk.Nghost; i <= SolnBlk.ICu+SolnBlk.Nghost; i++) {
	  SolnBlk.W[i][j] = Wo[0];
	  SolnBlk.U[i][j] = U(SolnBlk.W[i][j]);
	}
      }
    } else if (SolnBlk.Flow_Type == FLOWTYPE_LAMINAR) {
      // Set the initial data to the analytic solution of a laminar pipe
      // flow in the x-direction.
      for (int j = SolnBlk.JCl-SolnBlk.Nghost; j <= SolnBlk.JCu+SolnBlk.Nghost; j++) {
	for (int i = SolnBlk.ICl-SolnBlk.Nghost; i <= SolnBlk.ICu+SolnBlk.Nghost; i++) {
	  SolnBlk.W[i][j] = ViscousPipeFlow(SolnBlk.W[i][j],SolnBlk.Grid.Cell[i][j].Xc,
					    IP.dp,IP.Pipe_Length,IP.Pipe_Radius);
	  SolnBlk.U[i][j] = U(SolnBlk.W[i][j]);
	}
      }
    } else {
      // Set the initial data consistent to a turbulent pipe flow.
      for (int j = SolnBlk.JCl-SolnBlk.Nghost; j <= SolnBlk.JCu+SolnBlk.Nghost; j++) {
	for (int i = SolnBlk.ICl-SolnBlk.Nghost; i <= SolnBlk.ICu+SolnBlk.Nghost; i++) {
	  SolnBlk.W[i][j] = TurbulentPipeFlow(Wo[0],SolnBlk.Grid.Cell[i][j].Xc,
					      IP.dp,IP.Pipe_Length,IP.Pipe_Radius,IP.Reynolds_Number);
	  SolnBlk.U[i][j] = U(SolnBlk.W[i][j]);
	}
      }
    }
    break;
  case IC_VISCOUS_FLAT_PLATE :
    // Set the initial data to the Blasius (exact) solution for the
    // laminar flow over a flat plate in the x-direction.
    for (int j = SolnBlk.JCl-SolnBlk.Nghost; j <= SolnBlk.JCu+SolnBlk.Nghost; j++) {
      for (int i = SolnBlk.ICl-SolnBlk.Nghost; i <= SolnBlk.ICu+SolnBlk.Nghost; i++) {
	if (SolnBlk.Grid.Cell[i][j].Xc.y >= ZERO) SolnBlk.W[i][j] = FlatPlate(Wo[0],SolnBlk.Grid.Cell[i][j].Xc,IP.Plate_Length,eta,f,fp,fpp);
	else SolnBlk.W[i][j] = FlatPlate(Wo[0],Vector2D(SolnBlk.Grid.Cell[i][j].Xc.x,-SolnBlk.Grid.Cell[i][j].Xc.y),IP.Plate_Length,eta,f,fp,fpp);
	if (SolnBlk.Flow_Type == FLOWTYPE_TURBULENT_RANS_K_OMEGA) {
	  SolnBlk.W[i][j].k = 0.05*(sqr(SolnBlk.W[i][j].v.x) + sqr(SolnBlk.W[i][j].v.y));
	  SolnBlk.W[i][j].omega = 100000.0;
	}
	SolnBlk.U[i][j] = U(SolnBlk.W[i][j]);
      }
    }
    break;
  case IC_VISCOUS_STOKES_FLOW :
    for (int j = SolnBlk.JCl-SolnBlk.Nghost; j <= SolnBlk.JCu+SolnBlk.Nghost; j++) {
      for (int i = SolnBlk.ICl-SolnBlk.Nghost; i <= SolnBlk.ICu+SolnBlk.Nghost; i++) {
	SolnBlk.W[i][j] = Wo[0];
	SolnBlk.U[i][j] = U(SolnBlk.W[i][j]);
      }
    }
    break;
  case IC_VISCOUS_DRIVEN_CAVITY_FLOW :
    for (int j = SolnBlk.JCl-SolnBlk.Nghost; j <= SolnBlk.JCu+SolnBlk.Nghost; j++) {
      for (int i = SolnBlk.ICl-SolnBlk.Nghost; i <= SolnBlk.ICu+SolnBlk.Nghost; i++) {
	SolnBlk.W[i][j] = Wo[0];
	SolnBlk.U[i][j] = U(SolnBlk.W[i][j]);
      }
    }
    break;
  case IC_VISCOUS_BACKWARD_FACING_STEP :
    for (int j = SolnBlk.JCl-SolnBlk.Nghost; j <= SolnBlk.JCu+SolnBlk.Nghost; j++) {
      for (int i = SolnBlk.ICl-SolnBlk.Nghost; i <= SolnBlk.ICu+SolnBlk.Nghost; i++) {
// 	SolnBlk.W[i][j] = BackwardFacingStep(Wo[0],SolnBlk.Grid.Cell[i][j].Xc,
// 					     IP.Step_Height,IP.Top_Wall_Deflection,
// 					     IP.Reynolds_Number,IP.Mach_Number);
// 	SolnBlk.U[i][j] = U(SolnBlk.W[i][j]);
      }
    }
    break;
  case IC_JET_FLOW:
    if (IP.Mach_Number > 1.0){
      double sample_Ke, sample_Ee;
      for (int j = SolnBlk.JCl-SolnBlk.Nghost; j <= SolnBlk.JCu+SolnBlk.Nghost; j++) {
	for (int i = SolnBlk.ICl-SolnBlk.Nghost; i <= SolnBlk.ICu+SolnBlk.Nghost; i++) {
	  SolnBlk.W[i][j] = Wo[0];	
	  SolnBlk.W[i][j].rho = Wo[0].p/(300.0*IP.Wo.R);
	  SolnBlk.W[i][j].v.x = ZERO;
	  SolnBlk.W[i][j].v.y = ZERO;
	  
	  if (SolnBlk.Grid.Cell[i][j].Xc.y <IP.Pipe_Radius){
	    SolnBlk.W[i][j].rho = SolnBlk.W[i][j].p/IP.Wo.R/419.4;
	    //SolnBlk.W[i][j].rho = SolnBlk.W[i][j].p/IP.Wo.R/620.4;
	    SolnBlk.W[i][j].v.x = IP.Mach_Number* SolnBlk.W[i][j].a();
	    SolnBlk.W[i][j].v.y = ZERO;
	    SolnBlk.W[i][j].k = sqr(0.1* SolnBlk.W[i][j].v.x);
	    SolnBlk.W[i][j].omega = IP.Step_Height*SolnBlk.W[i][j].k;
	    sample_Ke = sqr(0.001*SolnBlk.W[i][j].p/(IP.Wo.gm1*SolnBlk.W[i][j].rho));
	    double CmuFmu = SolnBlk.W[i][j].muT()*SolnBlk.W[i][j].epsilon()/SolnBlk.W[i][j].rho/max(sqr(SolnBlk.W[i][j].k),sqr(TOLER));
	    double Coeff = sqr(0.7*0.14*SolnBlk.W[i][j].f_lambda(SolnBlk.Wall[i][j].ywall)/max(CmuFmu,sqr(TOLER)));
	    sample_Ee = Coeff*SolnBlk.W[i][j].epsilon()*sample_Ke/ max(SolnBlk.W[i][j].k,TOLER);
	  }// else if (SolnBlk.Grid.Cell[i][j].Xc.y < TWO*IP.Pipe_Radius && SolnBlk.Grid.Cell[i][j].Xc.x > IP.Pipe_Radius/10.0){
// 	    double temper = 419.0-119.0*(SolnBlk.Grid.Cell[i][j].Xc.y-IP.Pipe_Radius)/FOUR/IP.Pipe_Radius;
// 	    SolnBlk.W[i][j].rho = Wo[0].p/(temper*IP.Wo.R);
// 	  }
	  if (SolnBlk.Flow_Type == FLOWTYPE_TURBULENT_RANS_K_OMEGA) {
	    if (SolnBlk.W[i][j].k == ZERO){
	      //Wedge lenth != physical wedge length -> its just used as an input parameter to give intial K level 
	      SolnBlk.W[i][j].k = IP.Wedge_Length;//*exp(-(SolnBlk.Grid.Cell[i][j].Xc.y-IP.Pipe_Radius)/IP.Pipe_Radius/10.0); 
	    }
	    if (SolnBlk.W[i][j].omega == ZERO){
	      //step height != physical step_height -> its just used as an input parameter to give intial Omega level
	      //Omega is assumed to be propotional to K
	      SolnBlk.W[i][j].omega = IP.Wedge_Length;
	    }
	  }
          
	  /********************************************************************************************
           * THE NEXT IF STATEMENT IS NEVER USED. REMOVE IT SOMETIME.                                 *
	   * The initial guess of Ke is taken to be a small fraction of the specific internal energy  *
	   * at a point and the dissipation rate, ee, calcualted so that the intital value of PrT=0.7 *
	   ********************************************************************************************/	
	  if (SolnBlk.Variable_Prandtl == ON){ 
	    SolnBlk.W[i][j].ke = sample_Ke;
	    SolnBlk.W[i][j].ee = sample_Ee;
	    if (SolnBlk.Grid.Cell[i][j].Xc.y >= IP.Pipe_Radius){
	      SolnBlk.W[i][j].ke = sqr(0.0001*SolnBlk.W[i][j].p/(IP.Wo.gm1*SolnBlk.W[i][j].rho)); 
	      double CmuFmu = SolnBlk.W[i][j].muT()*SolnBlk.W[i][j].epsilon()/SolnBlk.W[i][j].rho/max(sqr(SolnBlk.W[i][j].k),sqr(TOLER));
	      double Coeff = sqr(0.7*0.14*SolnBlk.W[i][j].f_lambda(SolnBlk.Wall[i][j].ywall)/max(CmuFmu,sqr(TOLER)));
	      SolnBlk.W[i][j].ee = Coeff*SolnBlk.W[i][j].epsilon()* SolnBlk.W[i][j].ke/max(SolnBlk.W[i][j].k,TOLER);
//  	      double CmuFmu = SolnBlk.W[i][j].muT()*SolnBlk.W[i][j].epsilon()/SolnBlk.W[i][j].rho/max(sqr(SolnBlk.W[i][j].k),sqr(TOLER));
//  	      double Coeff = sqr(CmuFmu/0.14/max(SolnBlk.W[i][j].f_lambda(SolnBlk.Wall[i][j].ywall),TOLER));
//  	      SolnBlk.W[i][j].ke = Coeff* SolnBlk.W[i][j].ee* SolnBlk.W[i][j].k/ max(SolnBlk.W[i][j].epsilon(),TOLER)/0.7/0.7;
	    }
	  }
	  SolnBlk.U[i][j] = U(SolnBlk.W[i][j]);
	}
      }

    }else {
      
    // We are using the exit solution of a pipe stored in a file called "exit.dat"
    // which has 'y' values in the first column and the other 6 primitive variables
    // in the following 6 colums in the ascending order of 'y', making the 
    // total number of columns (variables) to be 7.

//       /*****************************************
//        *Reading solution from the exit.dat file*
//        *****************************************/

//       int numColumn; // the number of columns(variables) in the file pipe exit file.
//       if (SolnBlk.Variable_Prandtl == ON){
// 	numColumn = 9;
//       }else numColumn = 7;
//       int numRows=0; // number of rows in the file
//       double num;
//       if (SolnBlk.Variable_Prandtl == ON){
// 	indata.open("exit_vp.dat"); //opens the file
//       }else{
// 	indata.open("exit.dat"); //opens the file
//       };
//       if(!indata) { // file couldn't be opened
// 	cerr << "\n Error: The exit.dat file could not be opened" << endl;
// 	exit(1);
//       }
      
//       //Lets find how many rows are there in the file
//       int columnCounter=0;
//       indata >> num;
//       while ( !indata.eof() ) { 
// 	columnCounter++;
// 	indata >>num;
// 	if(columnCounter == numColumn){
// 	  numRows++;
// 	  columnCounter = 0;
// 	}  
//       }
//       indata.close();

//       // allocate a double array of numColumn*numRows
//       double** storedColumn;
//       int j,k=0;
//       storedColumn = (double **) malloc(sizeof(double *) * numColumn);
//       for(k=0; k<numColumn; k++)
// 	storedColumn[k] = (double *) malloc(sizeof(double)*numRows);
   
//       // allocate the file data in the double array.
//       columnCounter=0;
//       int currentRow = 0;  
//       if (SolnBlk.Variable_Prandtl == ON){
// 	comparedata.open("exit_vp.dat"); //opens the file
//       }else{
// 	comparedata.open("exit.dat"); //opens the file
//       };
//       if(!comparedata) { // file couldn't be opened
// 	cerr << "\n Error: The exit.dat file could not be opened" << endl;
// 	exit(1);
//       }

//       comparedata >> num;
//       while ( !comparedata.eof() ) { // keep reading until end-of-file
// 	storedColumn[columnCounter][currentRow] = num;
// 	//Increase the downstream pipe pressure equal to the 
// 	//atmospheric pressure.
// // 	if (columnCounter==4){
// // 	  storedColumn[columnCounter][currentRow]=Wo[0].p;
// // 	}
// 	columnCounter++;
// 	if(columnCounter == numColumn){
// 	  columnCounter = 0;
// 	  currentRow ++;
// 	}
// 	comparedata >> num; // sets EOF flag if no value found
//       }
//       comparedata.close();

//       /**************************************
//       Reading done from the exit.dat file
//       ****************************************/
//       //Finding the maximum value of K and Omega in the Pipe
//       double maxK = ZERO, maxOmega = ZERO;
//       for(k=0; k<numRows; k++){
// 	if (maxK < storedColumn[5][k])
// 	  maxK = storedColumn[5][k];
// 	if (maxOmega < storedColumn[6][k])
// 	  maxOmega = storedColumn[6][k];
//       }

//       double interpol; //stores the interpolation factor between two points
//       int foundYval; // a flag to see that whether or not we need to go ahead 
//       for (j = SolnBlk.JCl-SolnBlk.Nghost; j <= SolnBlk.JCu+SolnBlk.Nghost; j++) {
// 	for (int i = SolnBlk.ICl-SolnBlk.Nghost; i <= SolnBlk.ICu+SolnBlk.Nghost; i++) {
// 	  SolnBlk.W[i][j] = Wo[0];
// 	  SolnBlk.W[i][j].rho = Wo[0].p/IP.Wo.R/295.0; 
// 	  SolnBlk.W[i][j].v.x = ZERO;
// 	  SolnBlk.W[i][j].v.y = ZERO;
// 	  //Adjust the freestream pressure is there is a pressure gradient.
// 	  SolnBlk.W[i][j].p = Wo[0].p + SolnBlk.Grid.Cell[i][j].Xc.x/100.0/IP.Pipe_Radius*IP.dp; 
// 	  double dist = SolnBlk.Grid.Cell[i][j].Xc.y/FIVE/IP.Pipe_Radius;
// 	  // If the Y value is greater than pipe radius we simply use the STDATM conditions 
// 	  // otherwise we interpolate the solutions to get the pipe profile at the exit. 
// 	  if (abs(SolnBlk.Grid.Cell[i][j].Xc.y) < IP.Pipe_Radius && SolnBlk.Grid.Cell[i][j].Xc.x<50.0*IP.Pipe_Radius*(1.0-abs(SolnBlk.Grid.Cell[i][j].Xc.y)/IP.Pipe_Radius)){//if 1 begins
// 	    //if (abs(SolnBlk.Grid.Cell[i][j].Xc.y) < IP.Pipe_Radius && abs(SolnBlk.Grid.Cell[i][j].Xc.x) < 20.0*IP.Pipe_Radius){
// 	    //if (abs(SolnBlk.Grid.Cell[i][j].Xc.y) < IP.Pipe_Radius ){
// 	    foundYval = 0; 
// 	    for(k=0;(k<numRows) && (foundYval == 0);k++){//'for loop' begins
// 	      if( (SolnBlk.Grid.Cell[i][j].Xc.y < storedColumn[0][k]) && foundYval == 0 ){//if 2 begins
// 		foundYval = 1;
// 		if (k==0){
// 		  interpol = 0.0;
// 		  k = 1;
// 		}else{ 
// 		  if (storedColumn[0][k] != storedColumn[0][k-1]){
// 		    interpol = (SolnBlk.Grid.Cell[i][j].Xc.y-storedColumn[0][k-1])/(storedColumn[0][k]-storedColumn[0][k-1]);
// 		  }else interpol = 0.0;
// 		  //Linear interpolation is performed in order to match the values of the exit profile and the grid points
// 		  SolnBlk.W[i][j].rho = storedColumn[1][k-1] + interpol*(storedColumn[1][k]-storedColumn[1][k-1]);
// 		  SolnBlk.W[i][j].v.x = storedColumn[2][k-1] + interpol*(storedColumn[2][k]-storedColumn[2][k-1]);
// 		  SolnBlk.W[i][j].v.y = storedColumn[3][k-1] + interpol*(storedColumn[3][k]-storedColumn[3][k-1]);
// 		  SolnBlk.W[i][j].p = storedColumn[4][k-1] + interpol*(storedColumn[4][k]-storedColumn[4][k-1]);
// 		  SolnBlk.W[i][j].k = storedColumn[5][k-1] + interpol*(storedColumn[5][k]-storedColumn[5][k-1]);	      
// 		  SolnBlk.W[i][j].omega = storedColumn[6][k-1] + interpol*(storedColumn[6][k]-storedColumn[6][k-1]);
// 		  if (SolnBlk.Variable_Prandtl == ON){
// 		    SolnBlk.W[i][j].ke = storedColumn[7][k-1] + interpol*(storedColumn[7][k]-storedColumn[7][k-1]);
// 		    SolnBlk.W[i][j].ee = storedColumn[8][k-1] + interpol*(storedColumn[8][k]-storedColumn[8][k-1]);
// 		  }
// 		}   
// 	      }//if 2 ends
// 	    }//for ends
// 	  }//if 1 ends
	  
      double maxK, maxOmega;
      for (int j = SolnBlk.JCl-SolnBlk.Nghost; j <= SolnBlk.JCu+SolnBlk.Nghost; j++) {
	for (int i = SolnBlk.ICl-SolnBlk.Nghost; i <= SolnBlk.ICu+SolnBlk.Nghost; i++) {
	  SolnBlk.W[i][j] = Wo[0];	
	  SolnBlk.W[i][j].rho = Wo[0].p/(288.8*IP.Wo.R);
	  SolnBlk.W[i][j].v.x = ZERO;
	  SolnBlk.W[i][j].v.y = ZERO;
	  
	  if (SolnBlk.Grid.Cell[i][j].Xc.y <IP.Pipe_Radius){
	    SolnBlk.W[i][j].rho = SolnBlk.W[i][j].p/IP.Wo.R/551.2;
	    SolnBlk.W[i][j].v.x = IP.Mach_Number* SolnBlk.W[i][j].a();
	    SolnBlk.W[i][j].v.y = ZERO;
	    SolnBlk.W[i][j].k = sqr(0.1* SolnBlk.W[i][j].v.x);
	    maxK = SolnBlk.W[i][j].k;
	    SolnBlk.W[i][j].omega = IP.Step_Height*SolnBlk.W[i][j].k;
	    maxOmega = SolnBlk.W[i][j].omega;
	  }
	  
	  /********************************************************************************************
	   * For the K-Omega initial guess we are introducing a damping function for both K and omega *
	   * which introduces some turbulent kinetic energy near the shear layer and damps out to     *
	   * extremely small values in the free stream.                                               *
	   ********************************************************************************************/
	  double dist = SolnBlk.Grid.Cell[i][j].Xc.y/TEN/IP.Pipe_Radius;
	  if (SolnBlk.Flow_Type == FLOWTYPE_TURBULENT_RANS_K_OMEGA) {
	    if (SolnBlk.W[i][j].k == ZERO)
	      SolnBlk.W[i][j].k = maxK*exp(-dist);
	    if (SolnBlk.W[i][j].omega == ZERO)
	      SolnBlk.W[i][j].omega = maxOmega*exp(-dist/IP.Step_Height);
	  }

	  /********************************************************************************************
	   * The initial guess of Ke is taken to be a small fraction of the specific internal energy  *
	   * at a point and the dissipation rate, ee, calcualted so that the intital value of PrT=0.7 *
	   ********************************************************************************************/	
	  if (SolnBlk.Variable_Prandtl == ON){ 
	    if (SolnBlk.W[i][j].ke == ZERO)
	      SolnBlk.W[i][j].ke = sqr(0.001*SolnBlk.W[i][j].p/(IP.Wo.gm1*SolnBlk.W[i][j].rho));
	    if (SolnBlk.W[i][j].ee == ZERO){
	      double CmuFmu = SolnBlk.W[i][j].muT()*SolnBlk.W[i][j].epsilon()/SolnBlk.W[i][j].rho/max(sqr(SolnBlk.W[i][j].k),sqr(TOLER));
	      double Coeff = sqr(0.7*0.14*SolnBlk.W[i][j].f_lambda(SolnBlk.Wall[i][j].ywall)/max(CmuFmu,sqr(TOLER)));
	      SolnBlk.W[i][j].ee = Coeff*SolnBlk.W[i][j].epsilon()* SolnBlk.W[i][j].ke/max(SolnBlk.W[i][j].k,TOLER);
	    }
	  }
	  
	  SolnBlk.U[i][j] = U(SolnBlk.W[i][j]);
	}
      }
      //      delete storedColumn, columnCounter, numRows, numColumn, foundYval;
    }
    break;
    
//   case IC_MIXING_LAYER:
//     if (IP.Mach_Number > ONE){
//       for (int j = SolnBlk.JCl-SolnBlk.Nghost; j <= SolnBlk.JCu+SolnBlk.Nghost; j++) {
// 	for (int i = SolnBlk.ICl-SolnBlk.Nghost; i <= SolnBlk.ICu+SolnBlk.Nghost; i++) {
	  
// 	  SolnBlk.W[i][j] = Wo[0];
// // 	  if (SolnBlk.Flow_Type == FLOWTYPE_TURBULENT_RANS_K_OMEGA) 
// // 	    SolnBlk.W[i][j].k = 0.01*(SolnBlk.W[i][j].v.x);
   
// 	  if (SolnBlk.Grid.Cell[i][j].Xc.y >= ZERO) {
// 	    SolnBlk.W[i][j].v.y = ZERO;
// 	  } else{ 
// 	    SolnBlk.W[i][j].rho = SolnBlk.W[i][j].rho * pow( ONE + IP.Wo.gm1/TWO*sqr(IP.Mach_Number), -ONE);
// // 	    SolnBlk.W[i][j].k = Wo[0].rho*SolnBlk.W[i][j].k/SolnBlk.W[i][j].rho;
// 	    SolnBlk.W[i][j].v.x = SolnBlk.W[i][j].v.y;
// 	    SolnBlk.W[i][j].v.y = ZERO;
// 	  }
// // 	  SolnBlk.W[i][j].omega = SolnBlk.W[i][j].rho * SolnBlk.W[i][j].k/ SolnBlk.W[i][j].mu()*100.0;

// 	 /********************************************************************************************
//           * For the K-Omega initial guess we are introducing a damping function for both K and omega *
//           * which introduces some turbulent kinetic energy near the shear layer and damps out to     *
//           * extremely small values in the free stream.                                               *
//           ********************************************************************************************/
// 	  //We are using "Step_Height" as an input variable just to set and test different initial values of Omega
// 	  if (SolnBlk.Flow_Type == FLOWTYPE_TURBULENT_RANS_K_OMEGA) {
// 	    SolnBlk.W[i][j].k = IP.Step_Height*(-90.0*abs(SolnBlk.Grid.Cell[i][j].Xc.y)/IP.Plate_Length+90.0)*pow(cos(atan(1.0)*TWO*SolnBlk.Grid.Cell[i][j].Xc.y/IP.Plate_Length),2.0);
// 	    SolnBlk.W[i][j].omega = IP.Wedge_Length*SolnBlk.W[i][j].k;
// 	    //IP.Wedge_Length*(-90.0*abs(SolnBlk.Grid.Cell[i][j].Xc.y)/IP.Plate_Length*TWO+90.0)*pow(cos(atan(1.0)*TWO*SolnBlk.Grid.Cell[i][j].Xc.y/IP.Plate_Length*TWO),2.0); 
// 	  }
	  
// 	  SolnBlk.U[i][j] = U(SolnBlk.W[i][j]);
// 	}
//       }
//     }else{
//       for (int j = SolnBlk.JCl-SolnBlk.Nghost; j <= SolnBlk.JCu+SolnBlk.Nghost; j++) {
// 	for (int i = SolnBlk.ICl-SolnBlk.Nghost; i <= SolnBlk.ICu+SolnBlk.Nghost; i++) {
// 	  SolnBlk.W[i][j] = Wo[0];
	  
// // 	  if (SolnBlk.Flow_Type == FLOWTYPE_TURBULENT_RANS_K_OMEGA) {
// // 	    SolnBlk.W[i][j].k = 0.01*(SolnBlk.W[i][j].v.x);
// // 	    SolnBlk.W[i][j].omega = SolnBlk.W[i][j].rho * SolnBlk.W[i][j].k/ SolnBlk.W[i][j].mu()*100.0;
// // 	  }
	  
// 	  if (SolnBlk.Grid.Cell[i][j].Xc.y >= ZERO) {
// 	    SolnBlk.W[i][j].v.y = ZERO;
// 	  } else{ 
// 	      SolnBlk.W[i][j].rho = SolnBlk.W[i][j].rho * pow( ONE + IP.Wo.gm1/TWO*sqr(IP.Mach_Number), -ONE);
// 	      // SolnBlk.W[i][j].rho = DENSITY_STDATM * pow( ONE + (IP.Wo.g-ONE)/TWO*sqr(IP.Mach_Number),-IP.Wo.g*(IP.Wo.g-ONE));
// 	      SolnBlk.W[i][j].v.x = SolnBlk.W[i][j].v.y;
// 	      SolnBlk.W[i][j].v.y = ZERO;
// 	  }
	  
// 	 /********************************************************************************************
//           * For the K-Omega initial guess we are introducing a damping function for both K and omega *
//           * which introduces some turbulent kinetic energy near the shear layer and damps out to     *
//           * extremely small values in the free stream.                                               *
//            ********************************************************************************************/
// 	  //We are using "Step_Height" as an input variable just to set and test different initial values of Omega
// 	  if (SolnBlk.Flow_Type == FLOWTYPE_TURBULENT_RANS_K_OMEGA) {
// 	    SolnBlk.W[i][j].k=IP.Step_Height*(-20.0*abs(SolnBlk.Grid.Cell[i][j].Xc.y)/IP.Plate_Length+20.0)*pow(cos(atan(1.0)*TWO*SolnBlk.Grid.Cell[i][j].Xc.y/IP.Plate_Length),2.0);
// 	    SolnBlk.W[i][j].omega = IP.Wedge_Length*SolnBlk.W[i][j].k;

// 	    //IP.Wedge_Length*(-90.0*abs(SolnBlk.Grid.Cell[i][j].Xc.y)/IP.Plate_Length*TWO+90.0)*pow(cos(atan(1.0)*TWO*SolnBlk.Grid.Cell[i][j].Xc.y/IP.Plate_Length*TWO),2.0);
// 	  }

// 	  SolnBlk.U[i][j] = U(SolnBlk.W[i][j]);
// 	}
//       }
//     }
//     break;
  case IC_MIXING_LAYER :
    for (int j = SolnBlk.JCl-SolnBlk.Nghost; j <= SolnBlk.JCu+SolnBlk.Nghost; j++) {
      for (int i = SolnBlk.ICl-SolnBlk.Nghost; i <= SolnBlk.ICu+SolnBlk.Nghost; i++) {
	SolnBlk.W[i][j] = Wo[0];
	double fac = ONE + IP.Wo.gm1/TWO*sqr(IP.Mach_Number);
	double fac1 = ONE + IP.Wo.gm1/TWO*sqr(IP.Mach_Number2);
	double high = 1.0, high_dens = 1.0;
	// This is in order to increase the stagnation pressure and temperature for high mach numbers
	// otherwise the value of temperature were becoming extremely low.
	if (IP.Mach_Number > 2.0){
	  high = 40.0;
	  high_dens = 15.0;
	}

	SolnBlk.W[i][j].p = high*PRESSURE_STDATM * pow(fac,-IP.Wo.g*IP.Wo.gm1i);
	SolnBlk.W[i][j].rho = high_dens*DENSITY_STDATM * pow(fac,-IP.Wo.gm1i) ;
// 	while (SolnBlk.W[i][j].p/IP.Wo.R/SolnBlk.W[i][j].rho < 110.0){
// 	  high = high+0.2;
// 	  high_dens = high_dens-0.1;
// 	  SolnBlk.W[i][j].p = high*PRESSURE_STDATM * pow(fac,-IP.Wo.g*IP.Wo.gm1i);
// 	  SolnBlk.W[i][j].rho = high_dens*DENSITY_STDATM * pow(fac,-IP.Wo.gm1i) ;
// 	}

	double mu = 1.0e-06;
	if (SolnBlk.Grid.Cell[i][j].Xc.y > ZERO) {
	  SolnBlk.W[i][j].v.x = sqrt(IP.Wo.g*SolnBlk.W[i][j].p/SolnBlk.W[i][j].rho)*IP.Mach_Number;
	  SolnBlk.W[i][j].k = 1.2*sqr(0.01*SolnBlk.W[i][j].v.x)/2.0;
	  SolnBlk.W[i][j].omega = 0.09*IP.Wedge_Length*SolnBlk.W[i][j].rho*sqr(SolnBlk.W[i][j].k)/mu;
	}else{
	  SolnBlk.W[i][j].rho = SolnBlk.W[i][j].rho*fac1/fac ;
	  SolnBlk.W[i][j].v.x = sqrt(IP.Wo.g*SolnBlk.W[i][j].p/SolnBlk.W[i][j].rho)*IP.Mach_Number2;
	  SolnBlk.W[i][j].k = 1.2*sqr(0.01*SolnBlk.W[i][j].v.x)/2.0;
	  SolnBlk.W[i][j].omega = 0.09*IP.Wedge_Length*SolnBlk.W[i][j].rho*sqr(SolnBlk.W[i][j].k)/mu;
	}
	//We are using "Step_Height" as an input variable just to set and test different initial values of Omega
// 	if (SolnBlk.Flow_Type == FLOWTYPE_TURBULENT_RANS_K_OMEGA) {
// 	  SolnBlk.W[i][j].k = IP.Step_Height*(-abs(SolnBlk.Grid.Cell[i][j].Xc.y)+2.001)*pow(cos(atan(1.0)*SolnBlk.Grid.Cell[i][j].Xc.y),4.0);
// 	  SolnBlk.W[i][j].omega = IP.Wedge_Length*SolnBlk.W[i][j].k;
	  //IP.Wedge_Length*(-90.0*abs(SolnBlk.Grid.Cell[i][j].Xc.y)/IP.Plate_Length*TWO+90.0)*pow(cos(atan(1.0)*TWO*SolnBlk.Grid.Cell[i][j].Xc.y/IP.Plate_Length*TWO),2.// 0);
// 	}
	  
	SolnBlk.U[i][j] = U(SolnBlk.W[i][j]);
      }
    }
    break;

  case IC_VISCOUS_BRANCHED_DUCT :
    for (int j = SolnBlk.JCl-SolnBlk.Nghost; j <= SolnBlk.JCu+SolnBlk.Nghost; j++) {
      for (int i = SolnBlk.ICl-SolnBlk.Nghost; i <= SolnBlk.ICu+SolnBlk.Nghost; i++) {
	SolnBlk.W[i][j] = Wo[0];
	if (SolnBlk.Grid.Cell[i][j].Xc.x < ZERO && SolnBlk.Grid.Cell[i][j].Xc.y < 0.1143) {
	  SolnBlk.W[i][j].v.x = Wo[0].v.x*(ONE - sqr((TWO*SolnBlk.Grid.Cell[i][j].Xc.y-0.1143)/0.1143));
	  SolnBlk.W[i][j].v.y = ZERO;
	}
	SolnBlk.U[i][j] = U(SolnBlk.W[i][j]);
      }
    }
    break;
  default:
    // Set the solution state to the initial state Wo[0].
    for (int j = SolnBlk.JCl-SolnBlk.Nghost; j <= SolnBlk.JCu+SolnBlk.Nghost; j++) {
      for (int i = SolnBlk.ICl-SolnBlk.Nghost; i <= SolnBlk.ICu+SolnBlk.Nghost; i++) {
	SolnBlk.W[i][j] = Wo[0];
	SolnBlk.U[i][j] = U(SolnBlk.W[i][j]);
      }
    }
    break;
  };

  // Set default values for the boundary condition reference states.
  for (int j = SolnBlk.JCl-SolnBlk.Nghost; j <= SolnBlk.JCu+SolnBlk.Nghost; j++) {
    if (j >= SolnBlk.JCl && j <= SolnBlk.JCu) {
      SolnBlk.WoW[j] = SolnBlk.W[SolnBlk.ICl][j];
      SolnBlk.WoE[j] = SolnBlk.W[SolnBlk.ICu][j];
    } else if (j < SolnBlk.JCl) {
      SolnBlk.WoW[j] = SolnBlk.W[SolnBlk.ICl][SolnBlk.JCl];
      SolnBlk.WoE[j] = SolnBlk.W[SolnBlk.ICu][SolnBlk.JCl];
    } else {
      SolnBlk.WoW[j] = SolnBlk.W[SolnBlk.ICl][SolnBlk.JCu];
      SolnBlk.WoE[j] = SolnBlk.W[SolnBlk.ICu][SolnBlk.JCu];
    }
  }
  for (int i = SolnBlk.ICl-SolnBlk.Nghost; i <= SolnBlk.ICu+SolnBlk.Nghost; i++) {
    if (i >= SolnBlk.ICl && i <= SolnBlk.ICu) {
      SolnBlk.WoS[i] = SolnBlk.W[i][SolnBlk.JCl];
      SolnBlk.WoN[i] = SolnBlk.W[i][SolnBlk.JCu];
    } else if (i < SolnBlk.ICl) {
      SolnBlk.WoS[i] = SolnBlk.W[SolnBlk.ICl][SolnBlk.JCl];
      SolnBlk.WoN[i] = SolnBlk.W[SolnBlk.ICl][SolnBlk.JCu];
    } else {
      SolnBlk.WoS[i] = SolnBlk.W[SolnBlk.ICu][SolnBlk.JCl];
      SolnBlk.WoN[i] = SolnBlk.W[SolnBlk.ICu][SolnBlk.JCu];
    }
  }

}

/**********************************************************************
 * Routine: BCs                                                       *
 *                                                                    *
 * Apply boundary conditions at boundaries of the specified           *
 * quadrilateral solution block.                                      *
 *                                                                    *
 **********************************************************************/
void BCs(NavierStokes2D_Quad_Block &SolnBlk, NavierStokes2D_Input_Parameters &IP) {

  Vector2D dX;
  NavierStokes2D_pState dW, W;

  // WEST and EAST boundary conditions.
  for (int j = SolnBlk.JCl-SolnBlk.Nghost; j <= SolnBlk.JCu+SolnBlk.Nghost; j++) {

    // WEST boundary.
    if ((j >= SolnBlk.JCl && j <= SolnBlk.JCu) ||
	(j < SolnBlk.JCl && 
	 (SolnBlk.Grid.BCtypeS[SolnBlk.ICl-1] == BC_NONE ||
	  SolnBlk.Grid.BCtypeS[SolnBlk.ICl-1] == BC_PERIODIC ||
	  SolnBlk.Grid.BCtypeS[SolnBlk.ICl-1] == BC_CONSTANT_EXTRAPOLATION ||
	  SolnBlk.Grid.BCtypeS[SolnBlk.ICl-1] == BC_LINEAR_EXTRAPOLATION ||
	  SolnBlk.Grid.BCtypeS[SolnBlk.ICl-1] == BC_CHARACTERISTIC ||
	  SolnBlk.Grid.BCtypeS[SolnBlk.ICl-1] == BC_INFLOW_SUBSONIC ||
	  SolnBlk.Grid.BCtypeS[SolnBlk.ICl-1] == BC_OUTFLOW_SUBSONIC ||
	  SolnBlk.Grid.BCtypeS[SolnBlk.ICl-1] == BC_FIXED_PRESSURE ||
	  SolnBlk.Grid.BCtypeS[SolnBlk.ICl-1] == BC_RINGLEB_FLOW)) ||
	(j > SolnBlk.JCu && 
	 (SolnBlk.Grid.BCtypeN[SolnBlk.ICl-1] == BC_NONE ||
	  SolnBlk.Grid.BCtypeN[SolnBlk.ICl-1] == BC_PERIODIC ||
	  SolnBlk.Grid.BCtypeN[SolnBlk.ICl-1] == BC_CONSTANT_EXTRAPOLATION ||
	  SolnBlk.Grid.BCtypeN[SolnBlk.ICl-1] == BC_LINEAR_EXTRAPOLATION ||
	  SolnBlk.Grid.BCtypeN[SolnBlk.ICl-1] == BC_CHARACTERISTIC ||
	  SolnBlk.Grid.BCtypeN[SolnBlk.ICl-1] == BC_INFLOW_SUBSONIC ||
	  SolnBlk.Grid.BCtypeN[SolnBlk.ICl-1] == BC_OUTFLOW_SUBSONIC ||
	  SolnBlk.Grid.BCtypeN[SolnBlk.ICl-1] == BC_FIXED_PRESSURE ||
	  SolnBlk.Grid.BCtypeN[SolnBlk.ICl-1] == BC_RINGLEB_FLOW))) {

      switch(SolnBlk.Grid.BCtypeW[j]) {
      case BC_NONE :
	break;
      case BC_FIXED :
	SolnBlk.W[SolnBlk.ICl-1][j] = SolnBlk.WoW[j];
	SolnBlk.U[SolnBlk.ICl-1][j] = U(SolnBlk.W[SolnBlk.ICl-1][j]);
	SolnBlk.W[SolnBlk.ICl-2][j] = SolnBlk.WoW[j];
	SolnBlk.U[SolnBlk.ICl-2][j] = U(SolnBlk.W[SolnBlk.ICl-2][j]);
	break;
      case BC_CONSTANT_EXTRAPOLATION :
	SolnBlk.W[SolnBlk.ICl-1][j] = SolnBlk.W[SolnBlk.ICl][j];
	SolnBlk.U[SolnBlk.ICl-1][j] = SolnBlk.U[SolnBlk.ICl][j];
	SolnBlk.W[SolnBlk.ICl-2][j] = SolnBlk.W[SolnBlk.ICl][j];
	SolnBlk.U[SolnBlk.ICl-2][j] = SolnBlk.U[SolnBlk.ICl][j];
	break;
      case BC_INFLOW_SUBSONIC :
 	if (IP.i_ICs == IC_VISCOUS_CHANNEL_FLOW) {
	  // Fixed mass flux (rho and v) and linear extrapolatation of p.
 	  dX = SolnBlk.Grid.Cell[SolnBlk.ICl-1][j].Xc - SolnBlk.Grid.Cell[SolnBlk.ICl][j].Xc; 
  	  SolnBlk.W[SolnBlk.ICl-1][j] = SolnBlk.WoW[j];
 	  SolnBlk.W[SolnBlk.ICl-1][j].p = SolnBlk.W[SolnBlk.ICl][j].p + ((SolnBlk.WoE[j].p - SolnBlk.WoW[j].p)/
									 fabs(SolnBlk.Grid.Cell[SolnBlk.ICu][j].Xc.x - 
									      SolnBlk.Grid.Cell[SolnBlk.ICl][j].Xc.x))*dX.x;
	  SolnBlk.W[SolnBlk.ICl-1][j] = ViscousChannelFlowVelocity(SolnBlk.W[SolnBlk.ICl-1][j],SolnBlk.Grid.Cell[SolnBlk.ICl-1][j].Xc,
								   SolnBlk.Vwall,IP.dp,IP.Box_Width,IP.Box_Height);
	  SolnBlk.U[SolnBlk.ICl-1][j] = U(SolnBlk.W[SolnBlk.ICl-1][j]);
 	  dX = SolnBlk.Grid.Cell[SolnBlk.ICl-2][j].Xc - SolnBlk.Grid.Cell[SolnBlk.ICl][j].Xc;
 	  SolnBlk.W[SolnBlk.ICl-2][j] = SolnBlk.WoW[j];
 	  SolnBlk.W[SolnBlk.ICl-2][j].p = SolnBlk.W[SolnBlk.ICl][j].p + ((SolnBlk.WoE[j].p - SolnBlk.WoW[j].p)/
									 fabs(SolnBlk.Grid.Cell[SolnBlk.ICu][j].Xc.x - 
									      SolnBlk.Grid.Cell[SolnBlk.ICl][j].Xc.x))*dX.x;
	  SolnBlk.W[SolnBlk.ICl-2][j] = ViscousChannelFlowVelocity(SolnBlk.W[SolnBlk.ICl-2][j],SolnBlk.Grid.Cell[SolnBlk.ICl-2][j].Xc,
								   SolnBlk.Vwall,IP.dp,IP.Box_Width,IP.Box_Height);
	  SolnBlk.U[SolnBlk.ICl-2][j] = U(SolnBlk.W[SolnBlk.ICl-2][j]);
	} else if (IP.i_Grid == GRID_PIPE && IP.FlowType == FLOWTYPE_LAMINAR) {
	  // Fixed mass flux (rho and v) and linear extrapolatation of p.
 	  dX = SolnBlk.Grid.Cell[SolnBlk.ICl-1][j].Xc - SolnBlk.Grid.Cell[SolnBlk.ICl][j].Xc; 
  	  SolnBlk.W[SolnBlk.ICl-1][j] = SolnBlk.WoW[j];
 	  SolnBlk.W[SolnBlk.ICl-1][j].p = SolnBlk.W[SolnBlk.ICl][j].p + ((SolnBlk.WoE[j].p - SolnBlk.WoW[j].p)/
									 fabs(SolnBlk.Grid.Cell[SolnBlk.ICu][j].Xc.x - 
									      SolnBlk.Grid.Cell[SolnBlk.ICl][j].Xc.x))*dX.x;
	  SolnBlk.U[SolnBlk.ICl-1][j] = U(SolnBlk.W[SolnBlk.ICl-1][j]);
 	  dX = SolnBlk.Grid.Cell[SolnBlk.ICl-2][j].Xc - SolnBlk.Grid.Cell[SolnBlk.ICl][j].Xc;
 	  SolnBlk.W[SolnBlk.ICl-2][j] = SolnBlk.WoW[j];
 	  SolnBlk.W[SolnBlk.ICl-2][j].p = SolnBlk.W[SolnBlk.ICl][j].p + ((SolnBlk.WoE[j].p - SolnBlk.WoW[j].p)/
									 fabs(SolnBlk.Grid.Cell[SolnBlk.ICu][j].Xc.x - 
									      SolnBlk.Grid.Cell[SolnBlk.ICl][j].Xc.x))*dX.x;
	  SolnBlk.U[SolnBlk.ICl-2][j] = U(SolnBlk.W[SolnBlk.ICl-2][j]);
 	} else if (IP.i_Grid == GRID_PIPE && IP.FlowType == FLOWTYPE_TURBULENT_RANS_K_OMEGA) {
	  // Fixed mass flux (rho and v) and linear extrapolatation of p.
 	  dX = SolnBlk.Grid.Cell[SolnBlk.ICl-1][j].Xc - SolnBlk.Grid.Cell[SolnBlk.ICl][j].Xc; 
  	  SolnBlk.W[SolnBlk.ICl-1][j] = SolnBlk.WoW[j];
	  SolnBlk.W[SolnBlk.ICl-1][j].p = SolnBlk.WoW[j].p - (TWO/THREE)*SolnBlk.W[SolnBlk.ICl][j].dk() + dX.x*IP.dp/IP.Pipe_Length;
	  SolnBlk.W[SolnBlk.ICl-1][j].k = SolnBlk.W[SolnBlk.ICl][j].k;
	  SolnBlk.W[SolnBlk.ICl-1][j].omega = SolnBlk.W[SolnBlk.ICl][j].omega;
	  SolnBlk.W[SolnBlk.ICl-1][j].ke = SolnBlk.W[SolnBlk.ICl][j].ke;
	  SolnBlk.W[SolnBlk.ICl-1][j].ee = SolnBlk.W[SolnBlk.ICl][j].ee;	  
	  SolnBlk.U[SolnBlk.ICl-1][j] = U(SolnBlk.W[SolnBlk.ICl-1][j]);
 	  dX = SolnBlk.Grid.Cell[SolnBlk.ICl-2][j].Xc - SolnBlk.Grid.Cell[SolnBlk.ICl][j].Xc;
 	  SolnBlk.W[SolnBlk.ICl-2][j] = SolnBlk.WoW[j];
	  SolnBlk.W[SolnBlk.ICl-2][j].p = SolnBlk.WoW[j].p - (TWO/THREE)*SolnBlk.W[SolnBlk.ICl][j].dk() + dX.x*IP.dp/IP.Pipe_Length;
	  SolnBlk.W[SolnBlk.ICl-2][j].k = SolnBlk.W[SolnBlk.ICl][j].k;
	  SolnBlk.W[SolnBlk.ICl-2][j].omega = SolnBlk.W[SolnBlk.ICl][j].omega;
	  SolnBlk.W[SolnBlk.ICl-2][j].ke = SolnBlk.W[SolnBlk.ICl][j].ke;
	  SolnBlk.W[SolnBlk.ICl-2][j].ee = SolnBlk.W[SolnBlk.ICl][j].ee;	  
	  SolnBlk.U[SolnBlk.ICl-2][j] = U(SolnBlk.W[SolnBlk.ICl-2][j]);
 	} else {
	  //// Constant extrapolation with fixed v.x and u.x.
	  //SolnBlk.W[SolnBlk.ICl-1][j]     = SolnBlk.WoW[j];
	  //SolnBlk.W[SolnBlk.ICl-1][j].v.x = SolnBlk.W[SolnBlk.ICl][j].v.x;
	  ////SolnBlk.W[SolnBlk.ICl-1][j].u.x = SolnBlk.W[SolnBlk.ICl][j].u.x;
	  //SolnBlk.U[SolnBlk.ICl-1][j]     = U(SolnBlk.W[SolnBlk.ICl-1][j]);
	  //SolnBlk.W[SolnBlk.ICl-2][j]     = SolnBlk.WoW[j];
	  //SolnBlk.W[SolnBlk.ICl-2][j].v.x = SolnBlk.W[SolnBlk.ICl][j].v.x;
	  ////SolnBlk.W[SolnBlk.ICl-2][j].u.x = SolnBlk.W[SolnBlk.ICl][j].u.x;
	  //SolnBlk.U[SolnBlk.ICl-2][j]     = U(SolnBlk.W[SolnBlk.ICl-2][j]);
 	  // Fixed rho and v and constant extrapolation of p.
 	  SolnBlk.W[SolnBlk.ICl-1][j]   = SolnBlk.WoW[j];
 	  SolnBlk.W[SolnBlk.ICl-1][j].p = SolnBlk.W[SolnBlk.ICl][j].p;
 	  SolnBlk.U[SolnBlk.ICl-1][j]   = U(SolnBlk.W[SolnBlk.ICl-1][j]);
 	  SolnBlk.W[SolnBlk.ICl-2][j]   = SolnBlk.WoW[j];
 	  SolnBlk.W[SolnBlk.ICl-2][j].p = SolnBlk.W[SolnBlk.ICl][j].p;
 	  SolnBlk.U[SolnBlk.ICl-2][j]   = U(SolnBlk.W[SolnBlk.ICl-2][j]);
 	}
	break;
      case BC_OUTFLOW_SUBSONIC :
	// Constant extrapolation with fixed pressure.
	SolnBlk.W[SolnBlk.ICl-1][j]   = SolnBlk.W[SolnBlk.ICl][j]; 
	SolnBlk.W[SolnBlk.ICl-1][j].p = SolnBlk.WoW[j].p;
	SolnBlk.U[SolnBlk.ICl-1][j]   = U(SolnBlk.W[SolnBlk.ICl-1][j]);
	SolnBlk.W[SolnBlk.ICl-2][j]   = SolnBlk.W[SolnBlk.ICl][j];
	SolnBlk.W[SolnBlk.ICl-2][j].p = SolnBlk.WoW[j].p;
	SolnBlk.U[SolnBlk.ICl-2][j]   = U(SolnBlk.W[SolnBlk.ICl-2][j]);
	break;
      case BC_FIXED_PRESSURE :
	// Constant extrapolation with fixed pressure.
	SolnBlk.W[SolnBlk.ICl-1][j]   = SolnBlk.W[SolnBlk.ICl][j]; 
	SolnBlk.W[SolnBlk.ICl-1][j].p = SolnBlk.WoW[j].p;
	SolnBlk.U[SolnBlk.ICl-1][j]   = U(SolnBlk.W[SolnBlk.ICl-1][j]);
	SolnBlk.W[SolnBlk.ICl-2][j]   = SolnBlk.W[SolnBlk.ICl][j];
	SolnBlk.W[SolnBlk.ICl-2][j].p = SolnBlk.WoW[j].p;
	SolnBlk.U[SolnBlk.ICl-2][j]   = U(SolnBlk.W[SolnBlk.ICl-2][j]);
	break;
      case BC_LINEAR_EXTRAPOLATION :
	Linear_Reconstruction_LeastSquares_2(SolnBlk,SolnBlk.ICl,j,LIMITER_BARTH_JESPERSEN);
	dX = SolnBlk.Grid.Cell[SolnBlk.ICl-1][j].Xc - SolnBlk.Grid.Cell[SolnBlk.ICl][j].Xc;
	SolnBlk.W[SolnBlk.ICl-1][j] = SolnBlk.W[SolnBlk.ICl][j] + 
	                              (SolnBlk.phi[SolnBlk.ICl][j]^SolnBlk.dWdx[SolnBlk.ICl][j])*dX.x +
	                              (SolnBlk.phi[SolnBlk.ICl][j]^SolnBlk.dWdy[SolnBlk.ICl][j])*dX.y;
	SolnBlk.U[SolnBlk.ICl-1][j] = U(SolnBlk.W[SolnBlk.ICl-1][j]);
	dX = SolnBlk.Grid.Cell[SolnBlk.ICl-2][j].Xc - SolnBlk.Grid.Cell[SolnBlk.ICl][j].Xc;
	SolnBlk.W[SolnBlk.ICl-2][j] = SolnBlk.W[SolnBlk.ICl][j] + 
	                              (SolnBlk.phi[SolnBlk.ICl][j]^SolnBlk.dWdx[SolnBlk.ICl][j])*dX.x +
	                              (SolnBlk.phi[SolnBlk.ICl][j]^SolnBlk.dWdy[SolnBlk.ICl][j])*dX.y;
	SolnBlk.U[SolnBlk.ICl-2][j] = U(SolnBlk.W[SolnBlk.ICl-2][j]);
	break;
      case BC_REFLECTION :
	SolnBlk.W[SolnBlk.ICl-1][j] = Reflect(SolnBlk.W[SolnBlk.ICl][j],SolnBlk.Grid.nfaceW(SolnBlk.ICl,j));
	SolnBlk.U[SolnBlk.ICl-1][j] = U(SolnBlk.W[SolnBlk.ICl-1][j]);
	SolnBlk.W[SolnBlk.ICl-2][j] = Reflect(SolnBlk.W[SolnBlk.ICl+1][j],SolnBlk.Grid.nfaceW(SolnBlk.ICl,j));
	SolnBlk.U[SolnBlk.ICl-2][j] = U(SolnBlk.W[SolnBlk.ICl-2][j]);
	break;
      case BC_WALL_VISCOUS_HEATFLUX :
	SolnBlk.W[SolnBlk.ICl-1][j] = WallViscousHeatFlux(SolnBlk.W[SolnBlk.ICl][j],SolnBlk.Grid.nfaceW(SolnBlk.ICl,j));
	SolnBlk.U[SolnBlk.ICl-1][j] = U(SolnBlk.W[SolnBlk.ICl-1][j]);
	SolnBlk.W[SolnBlk.ICl-2][j] = WallViscousHeatFlux(SolnBlk.W[SolnBlk.ICl+1][j],SolnBlk.Grid.nfaceW(SolnBlk.ICl,j));
	SolnBlk.U[SolnBlk.ICl-2][j] = U(SolnBlk.W[SolnBlk.ICl-2][j]);
	break;
      case BC_WALL_VISCOUS_ISOTHERMAL :
	SolnBlk.W[SolnBlk.ICl-1][j] = WallViscousIsothermal(SolnBlk.W[SolnBlk.ICl][j],SolnBlk.Grid.nfaceW(SolnBlk.ICl,j),SolnBlk.Twall);
	SolnBlk.U[SolnBlk.ICl-1][j] = U(SolnBlk.W[SolnBlk.ICl-1][j]);
	SolnBlk.W[SolnBlk.ICl-2][j] = WallViscousIsothermal(SolnBlk.W[SolnBlk.ICl+1][j],SolnBlk.Grid.nfaceW(SolnBlk.ICl,j),SolnBlk.Twall);
	SolnBlk.U[SolnBlk.ICl-2][j] = U(SolnBlk.W[SolnBlk.ICl-2][j]);
	break;
      case BC_MOVING_WALL_HEATFLUX :
	SolnBlk.W[SolnBlk.ICl-1][j] = MovingWallHeatFlux(SolnBlk.W[SolnBlk.ICl][j],SolnBlk.Grid.nfaceW(SolnBlk.ICl,j),SolnBlk.Vwall.x);
	SolnBlk.U[SolnBlk.ICl-1][j] = U(SolnBlk.W[SolnBlk.ICl-1][j]);
	SolnBlk.W[SolnBlk.ICl-2][j] = MovingWallHeatFlux(SolnBlk.W[SolnBlk.ICl+1][j],SolnBlk.Grid.nfaceW(SolnBlk.ICl,j),SolnBlk.Vwall.x);
	SolnBlk.U[SolnBlk.ICl-2][j] = U(SolnBlk.W[SolnBlk.ICl-2][j]);
	break;
      case BC_MOVING_WALL_ISOTHERMAL :
	SolnBlk.W[SolnBlk.ICl-1][j] = MovingWallIsothermal(SolnBlk.W[SolnBlk.ICl][j],SolnBlk.Grid.nfaceW(SolnBlk.ICl,j),SolnBlk.Vwall.x,SolnBlk.Twall);
	SolnBlk.U[SolnBlk.ICl-1][j] = U(SolnBlk.W[SolnBlk.ICl-1][j]);
	SolnBlk.W[SolnBlk.ICl-2][j] = MovingWallIsothermal(SolnBlk.W[SolnBlk.ICl+1][j],SolnBlk.Grid.nfaceW(SolnBlk.ICl,j),SolnBlk.Vwall.x,SolnBlk.Twall);
	SolnBlk.U[SolnBlk.ICl-2][j] = U(SolnBlk.W[SolnBlk.ICl-2][j]);
	break;
      case BC_BURNING_SURFACE :
	Linear_Reconstruction_LeastSquares_2(SolnBlk,SolnBlk.ICl,j,LIMITER_BARTH_JESPERSEN);
	dX = SolnBlk.Grid.Cell[SolnBlk.ICl-1][j].Xc - SolnBlk.Grid.Cell[SolnBlk.ICl][j].Xc;
	dW = SolnBlk.W[SolnBlk.ICl][j] + (SolnBlk.phi[SolnBlk.ICl][j]^SolnBlk.dWdx[SolnBlk.ICl][j])*dX.x +
                                         (SolnBlk.phi[SolnBlk.ICl][j]^SolnBlk.dWdy[SolnBlk.ICl][j])*dX.y;
	SolnBlk.W[SolnBlk.ICl-1][j] = BurningSurface(dW,SolnBlk.Grid.nfaceW(SolnBlk.ICl,j));
	SolnBlk.U[SolnBlk.ICl-1][j] = U(SolnBlk.W[SolnBlk.ICl-1][j]);
	SolnBlk.W[SolnBlk.ICl-2][j] = SolnBlk.W[SolnBlk.ICl-1][j];
	SolnBlk.U[SolnBlk.ICl-2][j] = SolnBlk.U[SolnBlk.ICl-1][j];
	break;
      case BC_MASS_INJECTION :
	//SolnBlk.W[SolnBlk.ICl-1][j] = MassInjection(SolnBlk.W[SolnBlk.ICl][j],SolnBlk.Grid.nfaceW(SolnBlk.ICl,j));
	SolnBlk.W[SolnBlk.ICl-1][j] = MassInjection2(SolnBlk.W[SolnBlk.ICl][j],
						     SolnBlk.Grid.xfaceW(SolnBlk.ICl,j),
						     SolnBlk.Grid.nfaceW(SolnBlk.ICl,j),
						     SolnBlk.Twall);
	SolnBlk.U[SolnBlk.ICl-1][j] = U(SolnBlk.W[SolnBlk.ICl-1][j]);
	SolnBlk.W[SolnBlk.ICl-2][j] = SolnBlk.W[SolnBlk.ICl-1][j];
	SolnBlk.U[SolnBlk.ICl-2][j] = SolnBlk.U[SolnBlk.ICl-1][j];
	break;
      case BC_PERIODIC :
	SolnBlk.W[SolnBlk.ICl-1][j] = SolnBlk.W[SolnBlk.ICu-1][j];
	SolnBlk.U[SolnBlk.ICl-1][j] = SolnBlk.U[SolnBlk.ICu-1][j];
	SolnBlk.W[SolnBlk.ICl-2][j] = SolnBlk.W[SolnBlk.ICu-2][j];
	SolnBlk.U[SolnBlk.ICl-2][j] = SolnBlk.U[SolnBlk.ICu-2][j];
	break;
      case BC_CHARACTERISTIC :
	SolnBlk.W[SolnBlk.ICl-1][j] = BC_Characteristic_Pressure(SolnBlk.W[SolnBlk.ICl][j],
								 SolnBlk.WoW[j],
								 SolnBlk.Grid.nfaceW(SolnBlk.ICl,j));
	SolnBlk.U[SolnBlk.ICl-1][j] = U(SolnBlk.W[SolnBlk.ICl-1][j]);
	SolnBlk.W[SolnBlk.ICl-2][j] = SolnBlk.W[SolnBlk.ICl-1][j];
	SolnBlk.U[SolnBlk.ICl-2][j] = U(SolnBlk.W[SolnBlk.ICl-2][j]);
	break;
      case BC_RINGLEB_FLOW :
// // 	SolnBlk.W[SolnBlk.ICl-1][j] = Reflect(SolnBlk.W[SolnBlk.ICl][j],SolnBlk.Grid.nfaceW(SolnBlk.ICl,j));
// // 	SolnBlk.U[SolnBlk.ICl-1][j] = U(SolnBlk.W[SolnBlk.ICl-1][j]);
// // 	SolnBlk.W[SolnBlk.ICl-2][j] = Reflect(SolnBlk.W[SolnBlk.ICl+1][j],SolnBlk.Grid.nfaceW(SolnBlk.ICl,j));
// // 	SolnBlk.U[SolnBlk.ICl-2][j] = U(SolnBlk.W[SolnBlk.ICl-2][j]);
// 	SolnBlk.W[SolnBlk.ICl-1][j] = RinglebFlow(SolnBlk.W[SolnBlk.ICl][j],SolnBlk.Grid.Cell[SolnBlk.ICl-1][j].Xc);
// 	SolnBlk.U[SolnBlk.ICl-1][j] = U(SolnBlk.W[SolnBlk.ICl-1][j]);
// 	SolnBlk.W[SolnBlk.ICl-2][j] = RinglebFlow(SolnBlk.W[SolnBlk.ICl-2][j],SolnBlk.Grid.Cell[SolnBlk.ICl-2][j].Xc);
// 	SolnBlk.U[SolnBlk.ICl-2][j] = U(SolnBlk.W[SolnBlk.ICl-2][j]);
 	SolnBlk.W[SolnBlk.ICl-1][j] = RinglebFlowAverageState(SolnBlk.W[SolnBlk.ICl-1][j],
							      SolnBlk.Grid.nodeSW(SolnBlk.ICl-1,j).X,SolnBlk.Grid.nodeSE(SolnBlk.ICl-1,j).X,
							      SolnBlk.Grid.nodeNE(SolnBlk.ICl-1,j).X,SolnBlk.Grid.nodeNW(SolnBlk.ICl-1,j).X);
 	SolnBlk.U[SolnBlk.ICl-1][j] = U(SolnBlk.W[SolnBlk.ICl-1][j]);
 	SolnBlk.W[SolnBlk.ICl-2][j] = RinglebFlowAverageState(SolnBlk.W[SolnBlk.ICl-2][j],
							      SolnBlk.Grid.nodeSW(SolnBlk.ICl-2,j).X,SolnBlk.Grid.nodeSE(SolnBlk.ICl-2,j).X,
							      SolnBlk.Grid.nodeNE(SolnBlk.ICl-2,j).X,SolnBlk.Grid.nodeNW(SolnBlk.ICl-2,j).X);
 	SolnBlk.U[SolnBlk.ICl-2][j] = U(SolnBlk.W[SolnBlk.ICl-2][j]);
	break;
      default:
	SolnBlk.W[SolnBlk.ICl-1][j] = SolnBlk.W[SolnBlk.ICl][j];
	SolnBlk.U[SolnBlk.ICl-1][j] = SolnBlk.U[SolnBlk.ICl][j];
	SolnBlk.W[SolnBlk.ICl-2][j] = SolnBlk.W[SolnBlk.ICl][j];
	SolnBlk.U[SolnBlk.ICl-2][j] = SolnBlk.U[SolnBlk.ICl][j];
	break;
      };
    }

    // EAST boundary.
    if ((j >= SolnBlk.JCl && j <= SolnBlk.JCu) ||
	(j < SolnBlk.JCl && 
	 (SolnBlk.Grid.BCtypeS[SolnBlk.ICu+1] == BC_NONE ||
	  SolnBlk.Grid.BCtypeS[SolnBlk.ICu+1] == BC_PERIODIC ||
	  SolnBlk.Grid.BCtypeS[SolnBlk.ICu+1] == BC_CONSTANT_EXTRAPOLATION ||
	  SolnBlk.Grid.BCtypeS[SolnBlk.ICu+1] == BC_LINEAR_EXTRAPOLATION ||
	  SolnBlk.Grid.BCtypeS[SolnBlk.ICu+1] == BC_CHARACTERISTIC ||
	  SolnBlk.Grid.BCtypeS[SolnBlk.ICu+1] == BC_INFLOW_SUBSONIC ||
	  SolnBlk.Grid.BCtypeS[SolnBlk.ICu+1] == BC_OUTFLOW_SUBSONIC ||
	  SolnBlk.Grid.BCtypeS[SolnBlk.ICu+1] == BC_FIXED_PRESSURE ||
	  SolnBlk.Grid.BCtypeS[SolnBlk.ICu+1] == BC_RINGLEB_FLOW)) ||
	(j > SolnBlk.JCu && 
	 (SolnBlk.Grid.BCtypeN[SolnBlk.ICu+1] == BC_NONE ||
	  SolnBlk.Grid.BCtypeN[SolnBlk.ICu+1] == BC_PERIODIC ||
	  SolnBlk.Grid.BCtypeN[SolnBlk.ICu+1] == BC_CONSTANT_EXTRAPOLATION ||
	  SolnBlk.Grid.BCtypeN[SolnBlk.ICu+1] == BC_LINEAR_EXTRAPOLATION ||
	  SolnBlk.Grid.BCtypeN[SolnBlk.ICu+1] == BC_CHARACTERISTIC ||
	  SolnBlk.Grid.BCtypeN[SolnBlk.ICu+1] == BC_INFLOW_SUBSONIC ||
	  SolnBlk.Grid.BCtypeN[SolnBlk.ICu+1] == BC_OUTFLOW_SUBSONIC ||
	  SolnBlk.Grid.BCtypeN[SolnBlk.ICu+1] == BC_FIXED_PRESSURE ||
	  SolnBlk.Grid.BCtypeN[SolnBlk.ICu+1] == BC_RINGLEB_FLOW))) {

      switch(SolnBlk.Grid.BCtypeE[j]) {
      case BC_NONE :
	break;
      case BC_FIXED :
	SolnBlk.W[SolnBlk.ICu+1][j] = SolnBlk.WoE[j];
	SolnBlk.U[SolnBlk.ICu+1][j] = U(SolnBlk.W[SolnBlk.ICu+1][j]);
	SolnBlk.W[SolnBlk.ICu+2][j] = SolnBlk.WoE[j];
	SolnBlk.U[SolnBlk.ICu+2][j] = U(SolnBlk.W[SolnBlk.ICu+2][j]);
	break;
      case BC_CONSTANT_EXTRAPOLATION :
	SolnBlk.W[SolnBlk.ICu+1][j] = SolnBlk.W[SolnBlk.ICu][j];
	SolnBlk.U[SolnBlk.ICu+1][j] = SolnBlk.U[SolnBlk.ICu][j];
	SolnBlk.W[SolnBlk.ICu+2][j] = SolnBlk.W[SolnBlk.ICu][j];
	SolnBlk.U[SolnBlk.ICu+2][j] = SolnBlk.U[SolnBlk.ICu][j];
	break;
      case BC_INFLOW_SUBSONIC :
	// Constant extrapolation with fixed v.x and u.x.
	SolnBlk.W[SolnBlk.ICu+1][j]     = SolnBlk.WoE[j];
	SolnBlk.W[SolnBlk.ICu+1][j].v.x = SolnBlk.W[SolnBlk.ICu][j].v.x;
	//SolnBlk.W[SolnBlk.ICu+1][j].u.x = SolnBlk.W[SolnBlk.ICu][j].u.x;
	SolnBlk.U[SolnBlk.ICu+1][j]     = U(SolnBlk.W[SolnBlk.ICu+1][j]);
	SolnBlk.W[SolnBlk.ICu+2][j]     = SolnBlk.WoE[j];
	SolnBlk.W[SolnBlk.ICu+2][j].v.x = SolnBlk.W[SolnBlk.ICu][j].v.x;
	//SolnBlk.W[SolnBlk.ICu+2][j].u.x = SolnBlk.W[SolnBlk.ICu][j].u.x;
	SolnBlk.U[SolnBlk.ICu+2][j]     = U(SolnBlk.W[SolnBlk.ICu+2][j]);
	break;
      case BC_OUTFLOW_SUBSONIC :
 	if (IP.i_ICs == IC_VISCOUS_CHANNEL_FLOW) {
	  // Constant extrapolation for rho, v.x, and v.y but linear
	  // extrapolation of p.
	  dX = SolnBlk.Grid.Cell[SolnBlk.ICu+1][j].Xc - SolnBlk.Grid.Cell[SolnBlk.ICu][j].Xc; 
	  SolnBlk.W[SolnBlk.ICu+1][j] = SolnBlk.W[SolnBlk.ICu][j];
 	  SolnBlk.W[SolnBlk.ICu+1][j].p = SolnBlk.W[SolnBlk.ICu][j].p + ((SolnBlk.WoE[j].p - SolnBlk.WoW[j].p)/
									 fabs(SolnBlk.Grid.Cell[SolnBlk.ICu][j].Xc.x - 
									      SolnBlk.Grid.Cell[SolnBlk.ICl][j].Xc.x))*dX.x;
	  SolnBlk.W[SolnBlk.ICu+1][j] = ViscousChannelFlowVelocity(SolnBlk.W[SolnBlk.ICu+1][j],SolnBlk.Grid.Cell[SolnBlk.ICu+1][j].Xc,
								   SolnBlk.Vwall,IP.dp,IP.Box_Width,IP.Box_Height);
	  SolnBlk.U[SolnBlk.ICu+1][j] = U(SolnBlk.W[SolnBlk.ICu+1][j]);
	  SolnBlk.W[SolnBlk.ICu+2][j] = SolnBlk.W[SolnBlk.ICu][j];
 	  dX = SolnBlk.Grid.Cell[SolnBlk.ICu+2][j].Xc - SolnBlk.Grid.Cell[SolnBlk.ICu][j].Xc; 
 	  SolnBlk.W[SolnBlk.ICu+2][j].p = SolnBlk.W[SolnBlk.ICu][j].p + ((SolnBlk.WoE[j].p - SolnBlk.WoW[j].p)/
									 fabs(SolnBlk.Grid.Cell[SolnBlk.ICu][j].Xc.x - 
									      SolnBlk.Grid.Cell[SolnBlk.ICl][j].Xc.x))*dX.x;
	  SolnBlk.W[SolnBlk.ICu+2][j] = ViscousChannelFlowVelocity(SolnBlk.W[SolnBlk.ICu+2][j],SolnBlk.Grid.Cell[SolnBlk.ICu+2][j].Xc,
								   SolnBlk.Vwall,IP.dp,IP.Box_Width,IP.Box_Height);
	  SolnBlk.U[SolnBlk.ICu+2][j] = U(SolnBlk.W[SolnBlk.ICu+2][j]);
	} else if (IP.i_Grid == GRID_PIPE && IP.FlowType == FLOWTYPE_LAMINAR) {
	  dX = SolnBlk.Grid.Cell[SolnBlk.ICu+1][j].Xc - SolnBlk.Grid.Cell[SolnBlk.ICu][j].Xc; 
	  SolnBlk.W[SolnBlk.ICu+1][j] = SolnBlk.W[SolnBlk.ICu][j];
 	  SolnBlk.W[SolnBlk.ICu+1][j].p = SolnBlk.W[SolnBlk.ICu][j].p + ((SolnBlk.WoE[j].p - SolnBlk.WoW[j].p)/
									 fabs(SolnBlk.Grid.Cell[SolnBlk.ICu][j].Xc.x - 
									      SolnBlk.Grid.Cell[SolnBlk.ICl][j].Xc.x))*dX.x;
	  SolnBlk.U[SolnBlk.ICu+1][j] = U(SolnBlk.W[SolnBlk.ICu+1][j]);
	  SolnBlk.W[SolnBlk.ICu+2][j] = SolnBlk.W[SolnBlk.ICu][j];
 	  dX = SolnBlk.Grid.Cell[SolnBlk.ICu+2][j].Xc - SolnBlk.Grid.Cell[SolnBlk.ICu][j].Xc; 
 	  SolnBlk.W[SolnBlk.ICu+2][j].p = SolnBlk.W[SolnBlk.ICu][j].p + ((SolnBlk.WoE[j].p - SolnBlk.WoW[j].p)/
									 fabs(SolnBlk.Grid.Cell[SolnBlk.ICu][j].Xc.x - 
									      SolnBlk.Grid.Cell[SolnBlk.ICl][j].Xc.x))*dX.x;
	  SolnBlk.U[SolnBlk.ICu+2][j] = U(SolnBlk.W[SolnBlk.ICu+2][j]);
 	} else if (IP.i_Grid == GRID_PIPE && IP.FlowType == FLOWTYPE_TURBULENT_RANS_K_OMEGA) {
	  // Constant extrapolation for rho, v.x, and v.y but linear
	  // extrapolation of p.
	  dX = SolnBlk.Grid.Cell[SolnBlk.ICu+1][j].Xc - SolnBlk.Grid.Cell[SolnBlk.ICu][j].Xc; 
	  SolnBlk.W[SolnBlk.ICu+1][j] = SolnBlk.W[SolnBlk.ICu][j];
	  SolnBlk.W[SolnBlk.ICu+1][j].p = (SolnBlk.WoE[j].p - (TWO/THREE)*SolnBlk.W[SolnBlk.ICu][j].dk()) + dX.x*IP.dp/IP.Pipe_Length;
	  SolnBlk.U[SolnBlk.ICu+1][j] = U(SolnBlk.W[SolnBlk.ICu+1][j]);
 	  dX = SolnBlk.Grid.Cell[SolnBlk.ICu+2][j].Xc - SolnBlk.Grid.Cell[SolnBlk.ICu][j].Xc; 
	  SolnBlk.W[SolnBlk.ICu+2][j] = SolnBlk.W[SolnBlk.ICu][j];
	  SolnBlk.W[SolnBlk.ICu+2][j].p = (SolnBlk.WoE[j].p - (TWO/THREE)*SolnBlk.W[SolnBlk.ICu][j].dk()) + dX.x*IP.dp/IP.Pipe_Length;
	  SolnBlk.U[SolnBlk.ICu+2][j] = U(SolnBlk.W[SolnBlk.ICu+2][j]);
 	} else {
 	  // Constant extrapolation with fixed pressure.
 	  SolnBlk.W[SolnBlk.ICu+1][j]   = SolnBlk.W[SolnBlk.ICu][j]; 
 	  SolnBlk.W[SolnBlk.ICu+1][j].p = SolnBlk.WoE[j].p;
 	  SolnBlk.U[SolnBlk.ICu+1][j]   = U(SolnBlk.W[SolnBlk.ICu+1][j]);
 	  SolnBlk.W[SolnBlk.ICu+2][j]   = SolnBlk.W[SolnBlk.ICu][j];
 	  SolnBlk.W[SolnBlk.ICu+2][j].p = SolnBlk.WoE[j].p;
 	  SolnBlk.U[SolnBlk.ICu+2][j]   = U(SolnBlk.W[SolnBlk.ICu+2][j]);
 	}
	break;
      case BC_FIXED_PRESSURE :
	// Constant extrapolation with fixed pressure.
	SolnBlk.W[SolnBlk.ICu+1][j]   = SolnBlk.W[SolnBlk.ICu][j]; 
	SolnBlk.W[SolnBlk.ICu+1][j].p = SolnBlk.WoE[j].p;
	SolnBlk.U[SolnBlk.ICu+1][j]   = U(SolnBlk.W[SolnBlk.ICu+1][j]);
	SolnBlk.W[SolnBlk.ICu+2][j]   = SolnBlk.W[SolnBlk.ICu][j];
	SolnBlk.W[SolnBlk.ICu+2][j].p = SolnBlk.WoE[j].p;
	SolnBlk.U[SolnBlk.ICu+2][j]   = U(SolnBlk.W[SolnBlk.ICu+2][j]);

	//Fixed Temperature Outflow (constant extrapoltaion of pressure)
// 	SolnBlk.W[SolnBlk.ICu+1][j]   = SolnBlk.W[SolnBlk.ICu][j]; 
// 	SolnBlk.W[SolnBlk.ICu+1][j].rho = SolnBlk.W[SolnBlk.ICu][j].rho*SolnBlk.W[SolnBlk.ICu][j].T()/(310.0);
// 	SolnBlk.U[SolnBlk.ICu+1][j]   = U(SolnBlk.W[SolnBlk.ICu+1][j]);
// 	SolnBlk.W[SolnBlk.ICu+2][j]   = SolnBlk.W[SolnBlk.ICu][j];
// 	SolnBlk.W[SolnBlk.ICu+2][j].rho = SolnBlk.W[SolnBlk.ICu][j].rho*SolnBlk.W[SolnBlk.ICu][j].T()/(310.0);
// 	SolnBlk.U[SolnBlk.ICu+2][j]   = U(SolnBlk.W[SolnBlk.ICu+2][j]);

	break;	  
      case BC_LINEAR_EXTRAPOLATION :
	Linear_Reconstruction_LeastSquares_2(SolnBlk,SolnBlk.ICu,j,LIMITER_BARTH_JESPERSEN);
	dX = SolnBlk.Grid.Cell[SolnBlk.ICu+1][j].Xc - SolnBlk.Grid.Cell[SolnBlk.ICu][j].Xc;
	SolnBlk.W[SolnBlk.ICu+1][j] = SolnBlk.W[SolnBlk.ICu][j] + 
                                      (SolnBlk.phi[SolnBlk.ICu][j]^SolnBlk.dWdx[SolnBlk.ICu][j])*dX.x +
                                      (SolnBlk.phi[SolnBlk.ICu][j]^SolnBlk.dWdy[SolnBlk.ICu][j])*dX.y;
	SolnBlk.U[SolnBlk.ICu+1][j] = U(SolnBlk.W[SolnBlk.ICu+1][j]);
	dX = SolnBlk.Grid.Cell[SolnBlk.ICu+2][j].Xc - SolnBlk.Grid.Cell[SolnBlk.ICu][j].Xc;
	SolnBlk.W[SolnBlk.ICu+2][j] = SolnBlk.W[SolnBlk.ICu][j] + 
                                      (SolnBlk.phi[SolnBlk.ICu][j]^SolnBlk.dWdx[SolnBlk.ICu][j])*dX.x +
                                      (SolnBlk.phi[SolnBlk.ICu][j]^SolnBlk.dWdy[SolnBlk.ICu][j])*dX.y;
	SolnBlk.U[SolnBlk.ICu+2][j] = U(SolnBlk.W[SolnBlk.ICu+2][j]);
	break;
      case BC_REFLECTION :
	SolnBlk.W[SolnBlk.ICu+1][j] = Reflect(SolnBlk.W[SolnBlk.ICu][j],SolnBlk.Grid.nfaceE(SolnBlk.ICu,j));
	SolnBlk.U[SolnBlk.ICu+1][j] = U(SolnBlk.W[SolnBlk.ICu+1][j]);
	SolnBlk.W[SolnBlk.ICu+2][j] = Reflect(SolnBlk.W[SolnBlk.ICu-1][j],SolnBlk.Grid.nfaceE(SolnBlk.ICu,j));
	SolnBlk.U[SolnBlk.ICu+2][j] = U(SolnBlk.W[SolnBlk.ICu+2][j]);
	break;
      case BC_WALL_VISCOUS_HEATFLUX :
	SolnBlk.W[SolnBlk.ICu+1][j] = WallViscousHeatFlux(SolnBlk.W[SolnBlk.ICu][j],SolnBlk.Grid.nfaceE(SolnBlk.ICu,j));
	SolnBlk.U[SolnBlk.ICu+1][j] = U(SolnBlk.W[SolnBlk.ICu+1][j]);
	SolnBlk.W[SolnBlk.ICu+2][j] = WallViscousHeatFlux(SolnBlk.W[SolnBlk.ICu-1][j],SolnBlk.Grid.nfaceE(SolnBlk.ICu,j));
	SolnBlk.U[SolnBlk.ICu+2][j] = U(SolnBlk.W[SolnBlk.ICu+2][j]);
	break;
      case BC_WALL_VISCOUS_ISOTHERMAL :
	SolnBlk.W[SolnBlk.ICu+1][j] = WallViscousIsothermal(SolnBlk.W[SolnBlk.ICu][j],
							    SolnBlk.Grid.nfaceE(SolnBlk.ICu,j),
							    SolnBlk.Twall);
	SolnBlk.U[SolnBlk.ICu+1][j] = U(SolnBlk.W[SolnBlk.ICu+1][j]);
	SolnBlk.W[SolnBlk.ICu+2][j] = WallViscousIsothermal(SolnBlk.W[SolnBlk.ICu-1][j],
							    SolnBlk.Grid.nfaceE(SolnBlk.ICu,j),
							    SolnBlk.Twall);
	SolnBlk.U[SolnBlk.ICu+2][j] = U(SolnBlk.W[SolnBlk.ICu+2][j]);
	break;
      case BC_MOVING_WALL_HEATFLUX :
	SolnBlk.W[SolnBlk.ICu+1][j] = MovingWallHeatFlux(SolnBlk.W[SolnBlk.ICu][j],
							 SolnBlk.Grid.nfaceE(SolnBlk.ICu,j),
							 SolnBlk.Vwall.x);
	SolnBlk.U[SolnBlk.ICu+1][j] = U(SolnBlk.W[SolnBlk.ICu+1][j]);
	SolnBlk.W[SolnBlk.ICu+2][j] = MovingWallHeatFlux(SolnBlk.W[SolnBlk.ICu-1][j],
							 SolnBlk.Grid.nfaceE(SolnBlk.ICu,j),
							 SolnBlk.Vwall.x);
	SolnBlk.U[SolnBlk.ICu+2][j] = U(SolnBlk.W[SolnBlk.ICu+2][j]);
	break;
      case BC_MOVING_WALL_ISOTHERMAL :
	SolnBlk.W[SolnBlk.ICu+1][j] = MovingWallIsothermal(SolnBlk.W[SolnBlk.ICu][j],
							   SolnBlk.Grid.nfaceE(SolnBlk.ICu,j),
							   SolnBlk.Vwall.x,SolnBlk.Twall);
	SolnBlk.U[SolnBlk.ICu+1][j] = U(SolnBlk.W[SolnBlk.ICu+1][j]);
	SolnBlk.W[SolnBlk.ICu+2][j] = MovingWallIsothermal(SolnBlk.W[SolnBlk.ICu-1][j],
							   SolnBlk.Grid.nfaceE(SolnBlk.ICu,j),
							   SolnBlk.Vwall.x,SolnBlk.Twall);
	SolnBlk.U[SolnBlk.ICu+2][j] = U(SolnBlk.W[SolnBlk.ICu+2][j]);
	break;
      case BC_BURNING_SURFACE :
	Linear_Reconstruction_LeastSquares_2(SolnBlk,SolnBlk.ICu,j,LIMITER_BARTH_JESPERSEN);
	dX = SolnBlk.Grid.Cell[SolnBlk.ICu+1][j].Xc - SolnBlk.Grid.Cell[SolnBlk.ICu][j].Xc;
	dW = SolnBlk.W[SolnBlk.ICu][j] + (SolnBlk.phi[SolnBlk.ICu][j]^SolnBlk.dWdx[SolnBlk.ICu][j])*dX.x +
                                         (SolnBlk.phi[SolnBlk.ICu][j]^SolnBlk.dWdy[SolnBlk.ICu][j])*dX.y;
	SolnBlk.W[SolnBlk.ICu+1][j] = BurningSurface(dW,SolnBlk.Grid.nfaceE(SolnBlk.ICu,j));
	SolnBlk.U[SolnBlk.ICu+1][j] = U(SolnBlk.W[SolnBlk.ICu+1][j]);
	SolnBlk.W[SolnBlk.ICu+2][j] = SolnBlk.W[SolnBlk.ICu+1][j];
	SolnBlk.U[SolnBlk.ICu+2][j] = SolnBlk.U[SolnBlk.ICu+1][j];
	break;
      case BC_MASS_INJECTION :
	//SolnBlk.W[SolnBlk.ICu+1][j] = MassInjection(SolnBlk.W[SolnBlk.ICu][j],SolnBlk.Grid.nfaceE(SolnBlk.ICu,j));
	SolnBlk.W[SolnBlk.ICu+1][j] = MassInjection2(SolnBlk.W[SolnBlk.ICu][j],
						     SolnBlk.Grid.xfaceE(SolnBlk.ICu,j),
						     SolnBlk.Grid.nfaceE(SolnBlk.ICu,j),
						     SolnBlk.Twall);
	SolnBlk.U[SolnBlk.ICu+1][j] = SolnBlk.U[SolnBlk.ICl+1][j];
	SolnBlk.W[SolnBlk.ICu+2][j] = SolnBlk.W[SolnBlk.ICl+2][j];
	SolnBlk.U[SolnBlk.ICu+2][j] = SolnBlk.U[SolnBlk.ICl+2][j];
	break;
      case BC_PERIODIC :
	SolnBlk.W[SolnBlk.ICu+1][j] = SolnBlk.W[SolnBlk.ICl+1][j];
	SolnBlk.U[SolnBlk.ICu+1][j] = SolnBlk.U[SolnBlk.ICl+1][j];
	SolnBlk.W[SolnBlk.ICu+2][j] = SolnBlk.W[SolnBlk.ICl+2][j];
	SolnBlk.U[SolnBlk.ICu+2][j] = SolnBlk.U[SolnBlk.ICl+2][j];
	break;
      case BC_CHARACTERISTIC :
	////SolnBlk.W[SolnBlk.ICu+1][j] = SolnBlk.WoE[j];
	////SolnBlk.U[SolnBlk.ICu+1][j] = U(SolnBlk.W[SolnBlk.ICu+1][j]);
	//SolnBlk.W[SolnBlk.ICu+2][j] = SolnBlk.W[SolnBlk.ICu+1][j];//SolnBlk.WoE[j];
	//SolnBlk.U[SolnBlk.ICu+2][j] = U(SolnBlk.W[SolnBlk.ICu+2][j]);
	SolnBlk.W[SolnBlk.ICu+1][j] = BC_Characteristic_Pressure(SolnBlk.W[SolnBlk.ICu][j],
								 SolnBlk.WoE[j],
								 SolnBlk.Grid.nfaceE(SolnBlk.ICu,j));
	SolnBlk.U[SolnBlk.ICu+1][j] = U(SolnBlk.W[SolnBlk.ICu+1][j]);
	SolnBlk.W[SolnBlk.ICu+2][j] = SolnBlk.W[SolnBlk.ICu+1][j];
	SolnBlk.U[SolnBlk.ICu+2][j] = U(SolnBlk.W[SolnBlk.ICu+2][j]);
	break;
      case BC_RINGLEB_FLOW :
// 	//xxxSolnBlk.W[SolnBlk.ICu+1][j] = Reflect(SolnBlk.W[SolnBlk.ICu][j],SolnBlk.Grid.nfaceE(SolnBlk.ICu,j));
// 	//xxxSolnBlk.U[SolnBlk.ICu+1][j] = U(SolnBlk.W[SolnBlk.ICu+1][j]);
// 	//xxxSolnBlk.W[SolnBlk.ICu+2][j] = Reflect(SolnBlk.W[SolnBlk.ICu-1][j],SolnBlk.Grid.nfaceE(SolnBlk.ICu,j));
// 	//xxxSolnBlk.U[SolnBlk.ICu+2][j] = U(SolnBlk.W[SolnBlk.ICu+2][j]);
// 	SolnBlk.W[SolnBlk.ICu+1][j] = RinglebFlow(SolnBlk.W[SolnBlk.ICu-1][j],SolnBlk.Grid.Cell[SolnBlk.ICu+1][j].Xc);
// 	SolnBlk.U[SolnBlk.ICu+1][j] = U(SolnBlk.W[SolnBlk.ICu+1][j]);
// 	SolnBlk.W[SolnBlk.ICu+2][j] = RinglebFlow(SolnBlk.W[SolnBlk.ICu-2][j],SolnBlk.Grid.Cell[SolnBlk.ICu+2][j].Xc);
// 	SolnBlk.U[SolnBlk.ICu+2][j] = U(SolnBlk.W[SolnBlk.ICu+2][j]);
 	SolnBlk.W[SolnBlk.ICu+1][j] = RinglebFlowAverageState(SolnBlk.W[SolnBlk.ICu+1][j],
							      SolnBlk.Grid.nodeSW(SolnBlk.ICu+1,j).X,SolnBlk.Grid.nodeSE(SolnBlk.ICu+1,j).X,
							      SolnBlk.Grid.nodeNE(SolnBlk.ICu+1,j).X,SolnBlk.Grid.nodeNW(SolnBlk.ICu+1,j).X);
 	SolnBlk.U[SolnBlk.ICu+1][j] = U(SolnBlk.W[SolnBlk.ICu+1][j]);
 	SolnBlk.W[SolnBlk.ICu+2][j] = RinglebFlowAverageState(SolnBlk.W[SolnBlk.ICu+2][j],
							      SolnBlk.Grid.nodeSW(SolnBlk.ICu+2,j).X,SolnBlk.Grid.nodeSE(SolnBlk.ICu+2,j).X,
							      SolnBlk.Grid.nodeNE(SolnBlk.ICu+2,j).X,SolnBlk.Grid.nodeNW(SolnBlk.ICu+2,j).X);
 	SolnBlk.U[SolnBlk.ICu+2][j] = U(SolnBlk.W[SolnBlk.ICu+2][j]);
	break;
      default:
	SolnBlk.W[SolnBlk.ICu+1][j] = SolnBlk.W[SolnBlk.ICu][j];
	SolnBlk.U[SolnBlk.ICu+1][j] = SolnBlk.U[SolnBlk.ICu][j];
	SolnBlk.W[SolnBlk.ICu+2][j] = SolnBlk.W[SolnBlk.ICu][j];
	SolnBlk.U[SolnBlk.ICu+2][j] = SolnBlk.U[SolnBlk.ICu][j];
	break;
      };
    }

  }

  // NORTH and SOUTH boundary conditions.
  for (int i = SolnBlk.ICl-SolnBlk.Nghost; i <= SolnBlk.ICu+SolnBlk.Nghost; i++) {

    // SOUTH boundary.
    switch(SolnBlk.Grid.BCtypeS[i]) {
    case BC_NONE :
      break;
    case BC_FIXED :
      SolnBlk.W[i][SolnBlk.JCl-1] = SolnBlk.WoS[i];
      SolnBlk.U[i][SolnBlk.JCl-1] = U(SolnBlk.W[i][SolnBlk.JCl-1]);
      SolnBlk.W[i][SolnBlk.JCl-2] = SolnBlk.WoS[i];
      SolnBlk.U[i][SolnBlk.JCl-2] = U(SolnBlk.W[i][SolnBlk.JCl-2]);
      break;
    case BC_CONSTANT_EXTRAPOLATION :
      SolnBlk.W[i][SolnBlk.JCl-1] = SolnBlk.W[i][SolnBlk.JCl];
      SolnBlk.U[i][SolnBlk.JCl-1] = SolnBlk.U[i][SolnBlk.JCl];
      SolnBlk.W[i][SolnBlk.JCl-2] = SolnBlk.W[i][SolnBlk.JCl];
      SolnBlk.U[i][SolnBlk.JCl-2] = SolnBlk.U[i][SolnBlk.JCl];
      break;
    case BC_INFLOW_SUBSONIC :
      // Constant extrapolation with fixed v.y and u.y.
      SolnBlk.W[i][SolnBlk.JCl-1]     = SolnBlk.WoS[i];
      SolnBlk.W[i][SolnBlk.JCl-1].v.y = SolnBlk.W[i][SolnBlk.JCu].v.y;
      //SolnBlk.W[i][SolnBlk.JCl-1].u.y = SolnBlk.W[i][SolnBlk.JCu].v.y;
      SolnBlk.U[i][SolnBlk.JCl-1]     = U(SolnBlk.W[i][SolnBlk.JCl-1]);  
      SolnBlk.W[i][SolnBlk.JCl-2]     = SolnBlk.WoS[i];
      SolnBlk.W[i][SolnBlk.JCl-2].v.y = SolnBlk.W[i][SolnBlk.JCu].v.y;
      //SolnBlk.W[i][SolnBlk.JCl-2].u.y = SolnBlk.W[i][SolnBlk.JCu].v.y;
      SolnBlk.U[i][SolnBlk.JCl-2]     = U(SolnBlk.W[i][SolnBlk.JCl-2]);
      break;
    case BC_OUTFLOW_SUBSONIC :
      // Constant extrapolation with fixed pressure.
      SolnBlk.W[i][SolnBlk.JCl-1]   = SolnBlk.W[i][SolnBlk.JCu];
      SolnBlk.W[i][SolnBlk.JCl-1].p = SolnBlk.WoS[i].p;
      SolnBlk.U[i][SolnBlk.JCl-1]   = U(SolnBlk.W[i][SolnBlk.JCl-1]);
      SolnBlk.W[i][SolnBlk.JCl-2]   = SolnBlk.W[i][SolnBlk.JCu];
      SolnBlk.W[i][SolnBlk.JCl-2].p = SolnBlk.WoS[i].p;
      SolnBlk.U[i][SolnBlk.JCl-2]   = U(SolnBlk.W[i][SolnBlk.JCl-2]);
      break;
    case BC_FIXED_PRESSURE :
      // Constant extrapolation with fixed pressure.
      SolnBlk.W[i][SolnBlk.JCl-1]   = SolnBlk.W[i][SolnBlk.JCu];
      SolnBlk.W[i][SolnBlk.JCl-1].p = SolnBlk.WoS[i].p;
      SolnBlk.U[i][SolnBlk.JCl-1]   = U(SolnBlk.W[i][SolnBlk.JCl-1]);
      SolnBlk.W[i][SolnBlk.JCl-2]   = SolnBlk.W[i][SolnBlk.JCu];
      SolnBlk.W[i][SolnBlk.JCl-2].p = SolnBlk.WoS[i].p;
      SolnBlk.U[i][SolnBlk.JCl-2]   = U(SolnBlk.W[i][SolnBlk.JCl-2]);
      break;
    case BC_LINEAR_EXTRAPOLATION :
      Linear_Reconstruction_LeastSquares_2(SolnBlk,i,SolnBlk.JCl,LIMITER_BARTH_JESPERSEN);
      dX = SolnBlk.Grid.Cell[i][SolnBlk.JCl-1].Xc - SolnBlk.Grid.Cell[i][SolnBlk.JCl].Xc;
      SolnBlk.W[i][SolnBlk.JCl-1] = SolnBlk.W[i][SolnBlk.JCl] + 
	                            (SolnBlk.phi[i][SolnBlk.JCl]^SolnBlk.dWdx[i][SolnBlk.JCl])*dX.x +
	                            (SolnBlk.phi[i][SolnBlk.JCl]^SolnBlk.dWdy[i][SolnBlk.JCl])*dX.y;
      SolnBlk.U[i][SolnBlk.JCl-1] = U(SolnBlk.W[i][SolnBlk.JCl-1]);
      dX = SolnBlk.Grid.Cell[i][SolnBlk.JCl-2].Xc - SolnBlk.Grid.Cell[i][SolnBlk.JCl].Xc;
      SolnBlk.W[i][SolnBlk.JCl-2] = SolnBlk.W[i][SolnBlk.JCl] + 
                	            (SolnBlk.phi[i][SolnBlk.JCl]^SolnBlk.dWdx[i][SolnBlk.JCl])*dX.x +
	                            (SolnBlk.phi[i][SolnBlk.JCl]^SolnBlk.dWdy[i][SolnBlk.JCl])*dX.y;
      SolnBlk.U[i][SolnBlk.JCl-2] = U(SolnBlk.W[i][SolnBlk.JCl-2]);
      break;
    case BC_REFLECTION :
      SolnBlk.W[i][SolnBlk.JCl-1] = Reflect(SolnBlk.W[i][SolnBlk.JCl],
					    SolnBlk.Grid.nfaceS(i,SolnBlk.JCl));
      SolnBlk.U[i][SolnBlk.JCl-1] = U(SolnBlk.W[i][SolnBlk.JCl-1]);
      SolnBlk.W[i][SolnBlk.JCl-2] = Reflect(SolnBlk.W[i][SolnBlk.JCl+1],
					    SolnBlk.Grid.nfaceS(i,SolnBlk.JCl));
      if (IP.Axisymmetric && IP.FlowType == FLOWTYPE_TURBULENT_RANS_K_OMEGA) {
	SolnBlk.W[i][SolnBlk.JCl-2].k = SolnBlk.W[i][SolnBlk.JCl].k;
	SolnBlk.W[i][SolnBlk.JCl-2].omega = SolnBlk.W[i][SolnBlk.JCl].omega;
 	SolnBlk.W[i][SolnBlk.JCl-2].ke = SolnBlk.W[i][SolnBlk.JCl].ke;
	SolnBlk.W[i][SolnBlk.JCl-2].ee = SolnBlk.W[i][SolnBlk.JCl].ee;
      }
      SolnBlk.U[i][SolnBlk.JCl-2] = U(SolnBlk.W[i][SolnBlk.JCl-2]);
      break;
    case BC_WALL_VISCOUS_HEATFLUX :
      SolnBlk.W[i][SolnBlk.JCl-1] = WallViscousHeatFlux(SolnBlk.W[i][SolnBlk.JCl],SolnBlk.Grid.nfaceS(i,SolnBlk.JCl));
      SolnBlk.U[i][SolnBlk.JCl-1] = U(SolnBlk.W[i][SolnBlk.JCl-1]);
      SolnBlk.W[i][SolnBlk.JCl-2] = WallViscousHeatFlux(SolnBlk.W[i][SolnBlk.JCl+1],SolnBlk.Grid.nfaceS(i,SolnBlk.JCl));
      SolnBlk.U[i][SolnBlk.JCl-2] = U(SolnBlk.W[i][SolnBlk.JCl-2]);
      break;
    case BC_WALL_VISCOUS_ISOTHERMAL :
      SolnBlk.W[i][SolnBlk.JCl-1] = WallViscousIsothermal(SolnBlk.W[i][SolnBlk.JCl],SolnBlk.Grid.nfaceS(i,SolnBlk.JCl),SolnBlk.Twall);
      SolnBlk.U[i][SolnBlk.JCl-1] = U(SolnBlk.W[i][SolnBlk.JCl-1]);
      SolnBlk.W[i][SolnBlk.JCl-2] = WallViscousIsothermal(SolnBlk.W[i][SolnBlk.JCl+1],SolnBlk.Grid.nfaceS(i,SolnBlk.JCl),SolnBlk.Twall);
      SolnBlk.U[i][SolnBlk.JCl-2] = U(SolnBlk.W[i][SolnBlk.JCl-2]);
      break;
    case BC_MOVING_WALL_HEATFLUX :
      SolnBlk.W[i][SolnBlk.JCl-1] = MovingWallHeatFlux(SolnBlk.W[i][SolnBlk.JCl],-SolnBlk.Grid.nfaceS(i,SolnBlk.JCl),SolnBlk.Vwall.x);
      SolnBlk.U[i][SolnBlk.JCl-1] = U(SolnBlk.W[i][SolnBlk.JCl-1]);
      SolnBlk.W[i][SolnBlk.JCl-2] = MovingWallHeatFlux(SolnBlk.W[i][SolnBlk.JCl+1],-SolnBlk.Grid.nfaceS(i,SolnBlk.JCl),SolnBlk.Vwall.x);
      SolnBlk.U[i][SolnBlk.JCl-2] = U(SolnBlk.W[i][SolnBlk.JCl-2]);
      break;
    case BC_MOVING_WALL_ISOTHERMAL :
      SolnBlk.W[i][SolnBlk.JCl-1] = MovingWallIsothermal(SolnBlk.W[i][SolnBlk.JCl],SolnBlk.Grid.nfaceS(i,SolnBlk.JCl),SolnBlk.Vwall.x,SolnBlk.Twall);
      SolnBlk.U[i][SolnBlk.JCl-1] = U(SolnBlk.W[i][SolnBlk.JCl-1]);
      SolnBlk.W[i][SolnBlk.JCl-2] = MovingWallIsothermal(SolnBlk.W[i][SolnBlk.JCl+1],SolnBlk.Grid.nfaceS(i,SolnBlk.JCl),SolnBlk.Vwall.x,SolnBlk.Twall);
      SolnBlk.U[i][SolnBlk.JCl-2] = U(SolnBlk.W[i][SolnBlk.JCl-2]);
      break;
    case BC_BURNING_SURFACE :
      Linear_Reconstruction_LeastSquares_2(SolnBlk,i,SolnBlk.JCl,LIMITER_BARTH_JESPERSEN);
      dX = SolnBlk.Grid.Cell[i][SolnBlk.JCl-1].Xc - SolnBlk.Grid.Cell[i][SolnBlk.JCl].Xc;
      dW = SolnBlk.W[i][SolnBlk.JCl] + (SolnBlk.phi[i][SolnBlk.JCl]^SolnBlk.dWdx[i][SolnBlk.JCl])*dX.x +
	                               (SolnBlk.phi[i][SolnBlk.JCl]^SolnBlk.dWdy[i][SolnBlk.JCl])*dX.y;
      SolnBlk.W[i][SolnBlk.JCl-1] = BurningSurface(dW,SolnBlk.Grid.nfaceS(i,SolnBlk.JCl));
      SolnBlk.U[i][SolnBlk.JCl-1] = U(SolnBlk.W[i][SolnBlk.JCl-1]);
      SolnBlk.W[i][SolnBlk.JCl-2] = SolnBlk.W[i][SolnBlk.JCl-1];
      SolnBlk.U[i][SolnBlk.JCl-2] = SolnBlk.U[i][SolnBlk.JCl-1];
      break;
    case BC_MASS_INJECTION :
      //SolnBlk.W[i][SolnBlk.JCl-1] = MassInjection(SolnBlk.W[i][SolnBlk.JCl],SolnBlk.Grid.nfaceS(i,SolnBlk.JCl));
      SolnBlk.W[i][SolnBlk.JCl-1] = MassInjection2(SolnBlk.W[i][SolnBlk.JCl],
						   SolnBlk.Grid.xfaceS(i,SolnBlk.JCl),
						   SolnBlk.Grid.nfaceS(i,SolnBlk.JCl),
						   SolnBlk.Twall);
      SolnBlk.U[i][SolnBlk.JCl-1] = SolnBlk.U[i][SolnBlk.JCu-1];
      SolnBlk.W[i][SolnBlk.JCl-2] = SolnBlk.W[i][SolnBlk.JCu-2];
      SolnBlk.U[i][SolnBlk.JCl-2] = SolnBlk.U[i][SolnBlk.JCu-2];
      break;
    case BC_PERIODIC :
      SolnBlk.W[i][SolnBlk.JCl-1] = SolnBlk.W[i][SolnBlk.JCu-1];
      SolnBlk.U[i][SolnBlk.JCl-1] = SolnBlk.U[i][SolnBlk.JCu-1];
      SolnBlk.W[i][SolnBlk.JCl-2] = SolnBlk.W[i][SolnBlk.JCu-2];
      SolnBlk.U[i][SolnBlk.JCl-2] = SolnBlk.U[i][SolnBlk.JCu-2];
      break;
    case BC_CHARACTERISTIC :
      SolnBlk.W[i][SolnBlk.JCl-1] = BC_Characteristic_Pressure(SolnBlk.W[i][SolnBlk.JCl],SolnBlk.WoS[i],SolnBlk.Grid.nfaceS(i,SolnBlk.JCl));
      SolnBlk.U[i][SolnBlk.JCl-1] = U(SolnBlk.W[i][SolnBlk.JCl-1]);
      SolnBlk.W[i][SolnBlk.JCl-2] = SolnBlk.W[i][SolnBlk.JCl-1];
      SolnBlk.U[i][SolnBlk.JCl-2] = U(SolnBlk.W[i][SolnBlk.JCl-2]);
      break;
    case BC_RINGLEB_FLOW :
//       SolnBlk.W[i][SolnBlk.JCl-1] = RinglebFlow(SolnBlk.W[i][SolnBlk.JCl-1],SolnBlk.Grid.Cell[i][SolnBlk.JCl-1].Xc);
//       SolnBlk.U[i][SolnBlk.JCl-1] = U(SolnBlk.W[i][SolnBlk.JCl-1]);
//       SolnBlk.W[i][SolnBlk.JCl-2] = RinglebFlow(SolnBlk.W[i][SolnBlk.JCl-2],SolnBlk.Grid.Cell[i][SolnBlk.JCl-2].Xc);
//       SolnBlk.U[i][SolnBlk.JCl-2] = U(SolnBlk.W[i][SolnBlk.JCl-2]);
      SolnBlk.W[i][SolnBlk.JCl-1] = RinglebFlowAverageState(SolnBlk.W[i][SolnBlk.JCl-1],
							    SolnBlk.Grid.nodeSW(i,SolnBlk.JCl-1).X,SolnBlk.Grid.nodeSE(i,SolnBlk.JCl-1).X,
							    SolnBlk.Grid.nodeNE(i,SolnBlk.JCl-1).X,SolnBlk.Grid.nodeNW(i,SolnBlk.JCl-1).X);
      SolnBlk.U[i][SolnBlk.JCl-1] = U(SolnBlk.W[i][SolnBlk.JCl-1]);
      SolnBlk.W[i][SolnBlk.JCl-2] = RinglebFlowAverageState(SolnBlk.W[i][SolnBlk.JCl-2],
							    SolnBlk.Grid.nodeSW(i,SolnBlk.JCl-2).X,SolnBlk.Grid.nodeSE(i,SolnBlk.JCl-2).X,
							    SolnBlk.Grid.nodeNE(i,SolnBlk.JCl-2).X,SolnBlk.Grid.nodeNW(i,SolnBlk.JCl-2).X);
      SolnBlk.U[i][SolnBlk.JCl-2] = U(SolnBlk.W[i][SolnBlk.JCl-2]);
      break;
    default:
      SolnBlk.W[i][SolnBlk.JCl-1] = SolnBlk.W[i][SolnBlk.JCl];
      SolnBlk.U[i][SolnBlk.JCl-1] = SolnBlk.U[i][SolnBlk.JCl];
      SolnBlk.W[i][SolnBlk.JCl-2] = SolnBlk.W[i][SolnBlk.JCl];
      SolnBlk.U[i][SolnBlk.JCl-2] = SolnBlk.U[i][SolnBlk.JCl];
      break;
    };

    // NORTH boundary.
    switch(SolnBlk.Grid.BCtypeN[i]) {
    case BC_NONE :
      break;
    case BC_FIXED :
      SolnBlk.W[i][SolnBlk.JCu+1] = SolnBlk.WoN[i];
      SolnBlk.U[i][SolnBlk.JCu+1] = U(SolnBlk.W[i][SolnBlk.JCu+1]);
      SolnBlk.W[i][SolnBlk.JCu+2] = SolnBlk.WoN[i];
      SolnBlk.U[i][SolnBlk.JCu+2] = U(SolnBlk.W[i][SolnBlk.JCu+2]);
      break;
    case BC_CONSTANT_EXTRAPOLATION :
      SolnBlk.W[i][SolnBlk.JCu+1] = SolnBlk.W[i][SolnBlk.JCu];
      SolnBlk.U[i][SolnBlk.JCu+1] = SolnBlk.U[i][SolnBlk.JCu];
      SolnBlk.W[i][SolnBlk.JCu+2] = SolnBlk.W[i][SolnBlk.JCu];
      SolnBlk.U[i][SolnBlk.JCu+2] = SolnBlk.U[i][SolnBlk.JCu];
      break;
    case BC_INFLOW_SUBSONIC :
      // All fixed except v.x (u) which is constant extrapolation
      SolnBlk.W[i][SolnBlk.JCu+1]     = SolnBlk.WoN[i];
      SolnBlk.W[i][SolnBlk.JCu+1].v.y = SolnBlk.W[i][SolnBlk.JCu].v.y;
      //SolnBlk.W[i][SolnBlk.JCu+1].u.y = SolnBlk.W[i][SolnBlk.JCu].u.y;
      SolnBlk.U[i][SolnBlk.JCu+1]     = U(SolnBlk.W[i][SolnBlk.JCu+1]);
      SolnBlk.W[i][SolnBlk.JCu+2]     = SolnBlk.WoN[i];
      SolnBlk.W[i][SolnBlk.JCu+2].v.y = SolnBlk.W[i][SolnBlk.JCu].v.y;
      //SolnBlk.W[i][SolnBlk.JCu+2].u.y = SolnBlk.W[i][SolnBlk.JCu].u.y;
      SolnBlk.U[i][SolnBlk.JCu+2]     = U(SolnBlk.W[i][SolnBlk.JCu+2]);
      break;
    case BC_OUTFLOW_SUBSONIC :
      // Constant extrapolation with fixed pressure.
      SolnBlk.W[i][SolnBlk.JCu+1]   = SolnBlk.W[i][SolnBlk.JCu];
      SolnBlk.W[i][SolnBlk.JCu+1].p = SolnBlk.WoN[i].p;
      SolnBlk.U[i][SolnBlk.JCu+1]   = U(SolnBlk.W[i][SolnBlk.JCu+1]);
      SolnBlk.W[i][SolnBlk.JCu+2]   = SolnBlk.W[i][SolnBlk.JCu];
      SolnBlk.W[i][SolnBlk.JCu+2].p = SolnBlk.WoN[i].p;
      SolnBlk.U[i][SolnBlk.JCu+2]   = U(SolnBlk.W[i][SolnBlk.JCu+2]);
      break;
    case BC_FIXED_PRESSURE :
      // Constant extrapolation with fixed pressure.
      SolnBlk.W[i][SolnBlk.JCu+1]   = SolnBlk.W[i][SolnBlk.JCu];
      SolnBlk.W[i][SolnBlk.JCu+1].p = SolnBlk.WoN[i].p;
      SolnBlk.U[i][SolnBlk.JCu+1]   = U(SolnBlk.W[i][SolnBlk.JCu+1]);
      SolnBlk.W[i][SolnBlk.JCu+2]   = SolnBlk.W[i][SolnBlk.JCu];
      SolnBlk.W[i][SolnBlk.JCu+2].p = SolnBlk.WoN[i].p;
      SolnBlk.U[i][SolnBlk.JCu+2]   = U(SolnBlk.W[i][SolnBlk.JCu+2]);

      //Fixed Temperature Outflow  (constant extrapoltaion of pressure)
//       SolnBlk.W[i][SolnBlk.JCu+1]   = SolnBlk.W[i][SolnBlk.JCu]; 
//       SolnBlk.W[i][SolnBlk.JCu+1].rho = SolnBlk.W[i][SolnBlk.JCu].rho*SolnBlk.W[i][SolnBlk.JCu].T()/(SolnBlk.WoN[i].T());
//       SolnBlk.U[i][SolnBlk.JCu+1]   = U(SolnBlk.W[i][SolnBlk.JCu+1]);
//       SolnBlk.W[i][SolnBlk.JCu+2]   = SolnBlk.W[i][SolnBlk.JCu];
//       SolnBlk.W[i][SolnBlk.JCu+2].rho = SolnBlk.W[i][SolnBlk.JCu].rho*SolnBlk.W[i][SolnBlk.JCu].T()/(SolnBlk.WoN[i].T());
//       SolnBlk.U[i][SolnBlk.JCu+2]   = U(SolnBlk.W[i][SolnBlk.JCu+2]);
      break;
    case BC_LINEAR_EXTRAPOLATION :
      Linear_Reconstruction_LeastSquares_2(SolnBlk,i,SolnBlk.JCu,LIMITER_BARTH_JESPERSEN);
      dX = SolnBlk.Grid.Cell[i][SolnBlk.JCu+1].Xc - SolnBlk.Grid.Cell[i][SolnBlk.JCu].Xc;
      SolnBlk.W[i][SolnBlk.JCu+1] = SolnBlk.W[i][SolnBlk.JCu] + 
	                            (SolnBlk.phi[i][SolnBlk.JCu]^SolnBlk.dWdx[i][SolnBlk.JCu])*dX.x +
	                            (SolnBlk.phi[i][SolnBlk.JCu]^SolnBlk.dWdy[i][SolnBlk.JCu])*dX.y;
      SolnBlk.U[i][SolnBlk.JCu+1] = U(SolnBlk.W[i][SolnBlk.JCu+1]);
      dX = SolnBlk.Grid.Cell[i][SolnBlk.JCu+2].Xc - SolnBlk.Grid.Cell[i][SolnBlk.JCu].Xc;
      SolnBlk.W[i][SolnBlk.JCu+2] = SolnBlk.W[i][SolnBlk.JCu] + 
	                            (SolnBlk.phi[i][SolnBlk.JCu]^SolnBlk.dWdx[i][SolnBlk.JCu])*dX.x +
	                            (SolnBlk.phi[i][SolnBlk.JCu]^SolnBlk.dWdy[i][SolnBlk.JCu])*dX.y;
      SolnBlk.U[i][SolnBlk.JCu+2] = U(SolnBlk.W[i][SolnBlk.JCu+2]);
      break;
    case BC_REFLECTION :
      SolnBlk.W[i][SolnBlk.JCu+1] = Reflect(SolnBlk.W[i][SolnBlk.JCu],SolnBlk.Grid.nfaceN(i,SolnBlk.JCu));
      SolnBlk.U[i][SolnBlk.JCu+1] = U(SolnBlk.W[i][SolnBlk.JCu+1]);
      SolnBlk.W[i][SolnBlk.JCu+2] = Reflect(SolnBlk.W[i][SolnBlk.JCu-1],SolnBlk.Grid.nfaceN(i,SolnBlk.JCu));
      SolnBlk.U[i][SolnBlk.JCu+2] = U(SolnBlk.W[i][SolnBlk.JCu+2]);
      break;
    case BC_WALL_VISCOUS_HEATFLUX :
      SolnBlk.W[i][SolnBlk.JCu+1] = WallViscousHeatFlux(SolnBlk.W[i][SolnBlk.JCu],SolnBlk.Grid.nfaceN(i,SolnBlk.JCu));
      SolnBlk.U[i][SolnBlk.JCu+1] = U(SolnBlk.W[i][SolnBlk.JCu+1]);
      SolnBlk.W[i][SolnBlk.JCu+2] = WallViscousHeatFlux(SolnBlk.W[i][SolnBlk.JCu-1],SolnBlk.Grid.nfaceN(i,SolnBlk.JCu));
      SolnBlk.U[i][SolnBlk.JCu+2] = U(SolnBlk.W[i][SolnBlk.JCu+2]);
      break;
    case BC_WALL_VISCOUS_ISOTHERMAL :
      SolnBlk.W[i][SolnBlk.JCu+1] = WallViscousIsothermal(SolnBlk.W[i][SolnBlk.JCu],SolnBlk.Grid.nfaceN(i,SolnBlk.JCu),SolnBlk.Twall);
      SolnBlk.U[i][SolnBlk.JCu+1] = U(SolnBlk.W[i][SolnBlk.JCu+1]);
      SolnBlk.W[i][SolnBlk.JCu+2] = WallViscousIsothermal(SolnBlk.W[i][SolnBlk.JCu-1],SolnBlk.Grid.nfaceN(i,SolnBlk.JCu),SolnBlk.Twall);
      SolnBlk.U[i][SolnBlk.JCu+2] = U(SolnBlk.W[i][SolnBlk.JCu+2]);
      break;
    case BC_MOVING_WALL_HEATFLUX :
      SolnBlk.W[i][SolnBlk.JCu+1] = MovingWallHeatFlux(SolnBlk.W[i][SolnBlk.JCu],SolnBlk.Grid.nfaceN(i,SolnBlk.JCu),SolnBlk.Vwall.x);
      SolnBlk.U[i][SolnBlk.JCu+1] = U(SolnBlk.W[i][SolnBlk.JCu+1]);
      SolnBlk.W[i][SolnBlk.JCu+2] = MovingWallHeatFlux(SolnBlk.W[i][SolnBlk.JCu-1],SolnBlk.Grid.nfaceN(i,SolnBlk.JCu),SolnBlk.Vwall.x);
      SolnBlk.U[i][SolnBlk.JCu+2] = U(SolnBlk.W[i][SolnBlk.JCu+2]);
      break;
    case BC_MOVING_WALL_ISOTHERMAL :
      SolnBlk.W[i][SolnBlk.JCu+1] = MovingWallIsothermal(SolnBlk.W[i][SolnBlk.JCu],SolnBlk.Grid.nfaceN(i,SolnBlk.JCu),SolnBlk.Vwall.x,SolnBlk.Twall);
      SolnBlk.U[i][SolnBlk.JCu+1] = U(SolnBlk.W[i][SolnBlk.JCu+1]);
      SolnBlk.W[i][SolnBlk.JCu+2] = MovingWallIsothermal(SolnBlk.W[i][SolnBlk.JCu-1],SolnBlk.Grid.nfaceN(i,SolnBlk.JCu),SolnBlk.Vwall.x,SolnBlk.Twall);
      SolnBlk.U[i][SolnBlk.JCu+2] = U(SolnBlk.W[i][SolnBlk.JCu+2]);
      break;
    case BC_BURNING_SURFACE :
      Linear_Reconstruction_LeastSquares_2(SolnBlk,i,SolnBlk.JCu,LIMITER_BARTH_JESPERSEN);
      dX = SolnBlk.Grid.Cell[i][SolnBlk.JCu+1].Xc - SolnBlk.Grid.Cell[i][SolnBlk.JCu].Xc;
      dW = SolnBlk.W[i][SolnBlk.JCu] + (SolnBlk.phi[i][SolnBlk.JCu]^SolnBlk.dWdx[i][SolnBlk.JCu])*dX.x +
	                               (SolnBlk.phi[i][SolnBlk.JCu]^SolnBlk.dWdy[i][SolnBlk.JCu])*dX.y;
      SolnBlk.W[i][SolnBlk.JCu+1] = BurningSurface(dW,SolnBlk.Grid.nfaceN(i,SolnBlk.JCu));
      SolnBlk.U[i][SolnBlk.JCu+1] = U(SolnBlk.W[i][SolnBlk.JCu+1]);
      SolnBlk.W[i][SolnBlk.JCu+2] = SolnBlk.W[i][SolnBlk.JCu+1];
      SolnBlk.U[i][SolnBlk.JCu+2] = SolnBlk.U[i][SolnBlk.JCu+1];
      break;
    case BC_MASS_INJECTION :
      //SolnBlk.W[i][SolnBlk.JCu+1] = MassInjection(SolnBlk.W[i][SolnBlk.JCu],SolnBlk.Grid.nfaceN(i,SolnBlk.JCu));
      SolnBlk.W[i][SolnBlk.JCu+1] = MassInjection2(SolnBlk.W[i][SolnBlk.JCu],
						   SolnBlk.Grid.xfaceN(i,SolnBlk.JCu),
						   SolnBlk.Grid.nfaceN(i,SolnBlk.JCu),
						   SolnBlk.Twall);
      SolnBlk.U[i][SolnBlk.JCu+1] = U(SolnBlk.W[i][SolnBlk.JCu+1]);
      SolnBlk.W[i][SolnBlk.JCu+2] = SolnBlk.W[i][SolnBlk.JCu+1];
      SolnBlk.U[i][SolnBlk.JCu+2] = SolnBlk.U[i][SolnBlk.JCu+1];
      break;
    case BC_PERIODIC :
      SolnBlk.W[i][SolnBlk.JCu+1] = SolnBlk.W[i][SolnBlk.JCl+1];
      SolnBlk.U[i][SolnBlk.JCu+1] = SolnBlk.U[i][SolnBlk.JCl+1];
      SolnBlk.W[i][SolnBlk.JCu+2] = SolnBlk.W[i][SolnBlk.JCl+2];
      SolnBlk.U[i][SolnBlk.JCu+2] = SolnBlk.U[i][SolnBlk.JCl+2];
      break;
    case BC_CHARACTERISTIC :
      SolnBlk.W[i][SolnBlk.JCu+1] = BC_Characteristic_Pressure(SolnBlk.W[i][SolnBlk.JCu],SolnBlk.WoN[i],SolnBlk.Grid.nfaceN(i,SolnBlk.JCu));
      SolnBlk.U[i][SolnBlk.JCu+1] = U(SolnBlk.W[i][SolnBlk.JCu+1]);
      SolnBlk.W[i][SolnBlk.JCu+2] = SolnBlk.W[i][SolnBlk.JCu+1];
      SolnBlk.U[i][SolnBlk.JCu+2] = U(SolnBlk.W[i][SolnBlk.JCu+2]);
      break;
    case BC_RINGLEB_FLOW :
//       //SolnBlk.W[i][SolnBlk.JCu+1] = RinglebFlow(SolnBlk.W[i][SolnBlk.JCu+1],SolnBlk.Grid.Cell[i][SolnBlk.JCu+1].Xc);
//       //SolnBlk.U[i][SolnBlk.JCu+1] = U(SolnBlk.W[i][SolnBlk.JCu+1]);
//       //SolnBlk.W[i][SolnBlk.JCu+2] = RinglebFlow(SolnBlk.W[i][SolnBlk.JCu+2],SolnBlk.Grid.Cell[i][SolnBlk.JCu+2].Xc);
//       //SolnBlk.U[i][SolnBlk.JCu+2] = U(SolnBlk.W[i][SolnBlk.JCu+2]);
      SolnBlk.W[i][SolnBlk.JCu+1] = RinglebFlowAverageState(SolnBlk.W[i][SolnBlk.JCu+1],
							    SolnBlk.Grid.nodeSW(i,SolnBlk.JCu+1).X,SolnBlk.Grid.nodeSE(i,SolnBlk.JCu+1).X,
							    SolnBlk.Grid.nodeNE(i,SolnBlk.JCu+1).X,SolnBlk.Grid.nodeNW(i,SolnBlk.JCu+1).X);
      SolnBlk.U[i][SolnBlk.JCu+1] = U(SolnBlk.W[i][SolnBlk.JCu+1]);
      SolnBlk.W[i][SolnBlk.JCu+2] = RinglebFlowAverageState(SolnBlk.W[i][SolnBlk.JCu+2],
							    SolnBlk.Grid.nodeSW(i,SolnBlk.JCu+2).X,SolnBlk.Grid.nodeSE(i,SolnBlk.JCu+2).X,
							    SolnBlk.Grid.nodeNE(i,SolnBlk.JCu+2).X,SolnBlk.Grid.nodeNW(i,SolnBlk.JCu+2).X);
      SolnBlk.U[i][SolnBlk.JCu+2] = U(SolnBlk.W[i][SolnBlk.JCu+2]);
      break;
    default:
      SolnBlk.W[i][SolnBlk.JCu+1] = SolnBlk.W[i][SolnBlk.JCu];
      SolnBlk.U[i][SolnBlk.JCu+1] = SolnBlk.U[i][SolnBlk.JCu];
      SolnBlk.W[i][SolnBlk.JCu+2] = SolnBlk.W[i][SolnBlk.JCu];
      SolnBlk.U[i][SolnBlk.JCu+2] = SolnBlk.U[i][SolnBlk.JCu];
      break;
    };

  }

  // BC fix for corner points with burning surfaces on either side.
  if ((SolnBlk.Grid.BCtypeW[SolnBlk.JCl] == BC_REFLECTION ||
       SolnBlk.Grid.BCtypeW[SolnBlk.JCl] == BC_WALL_VISCOUS_HEATFLUX ||
       SolnBlk.Grid.BCtypeW[SolnBlk.JCl] == BC_WALL_VISCOUS_ISOTHERMAL ||
       SolnBlk.Grid.BCtypeW[SolnBlk.JCl] == BC_MOVING_WALL_HEATFLUX ||
       SolnBlk.Grid.BCtypeW[SolnBlk.JCl] == BC_MOVING_WALL_ISOTHERMAL ||
       SolnBlk.Grid.BCtypeW[SolnBlk.JCl] == BC_BURNING_SURFACE// ||
       //SolnBlk.Grid.BCtypeW[SolnBlk.JCl] == BC_MASS_INJECTION
       ) &&
      (SolnBlk.Grid.BCtypeS[SolnBlk.ICl] == BC_REFLECTION ||
       SolnBlk.Grid.BCtypeS[SolnBlk.ICl] == BC_WALL_VISCOUS_HEATFLUX ||
       SolnBlk.Grid.BCtypeS[SolnBlk.ICl] == BC_WALL_VISCOUS_ISOTHERMAL ||
       SolnBlk.Grid.BCtypeS[SolnBlk.ICl] == BC_MOVING_WALL_HEATFLUX ||
       SolnBlk.Grid.BCtypeS[SolnBlk.ICl] == BC_MOVING_WALL_ISOTHERMAL ||
       SolnBlk.Grid.BCtypeS[SolnBlk.ICl] == BC_BURNING_SURFACE// ||
       //SolnBlk.Grid.BCtypeS[SolnBlk.ICl] == BC_MASS_INJECTION
       )) {
    SolnBlk.W[SolnBlk.ICl-2][SolnBlk.JCl-2] = HALF*(SolnBlk.W[SolnBlk.ICl-2][SolnBlk.JCl  ]+
						    SolnBlk.W[SolnBlk.ICl  ][SolnBlk.JCl-2]);
    SolnBlk.U[SolnBlk.ICl-2][SolnBlk.JCl-2] = U(SolnBlk.W[SolnBlk.ICl-2][SolnBlk.JCl-2]);
    SolnBlk.W[SolnBlk.ICl-1][SolnBlk.JCl-2] = HALF*(SolnBlk.W[SolnBlk.ICl-1][SolnBlk.JCl  ]+
						    SolnBlk.W[SolnBlk.ICl  ][SolnBlk.JCl-2]);
    SolnBlk.U[SolnBlk.ICl-1][SolnBlk.JCl-2] = U(SolnBlk.W[SolnBlk.ICl-1][SolnBlk.JCl-2]);
    SolnBlk.W[SolnBlk.ICl-2][SolnBlk.JCl-1] = HALF*(SolnBlk.W[SolnBlk.ICl-2][SolnBlk.JCl  ]+
						    SolnBlk.W[SolnBlk.ICl  ][SolnBlk.JCl-1]);
    SolnBlk.U[SolnBlk.ICl-2][SolnBlk.JCl-1] = U(SolnBlk.W[SolnBlk.ICl-2][SolnBlk.JCl-1]);
    SolnBlk.W[SolnBlk.ICl-1][SolnBlk.JCl-1] = HALF*(SolnBlk.W[SolnBlk.ICl-1][SolnBlk.JCl  ]+
						    SolnBlk.W[SolnBlk.ICl  ][SolnBlk.JCl-1]);
    SolnBlk.U[SolnBlk.ICl-1][SolnBlk.JCl-1] = U(SolnBlk.W[SolnBlk.ICl-1][SolnBlk.JCl-1]);
  }

  if ((SolnBlk.Grid.BCtypeW[SolnBlk.JCu] == BC_REFLECTION ||
       SolnBlk.Grid.BCtypeW[SolnBlk.JCu] == BC_WALL_VISCOUS_HEATFLUX ||
       SolnBlk.Grid.BCtypeW[SolnBlk.JCu] == BC_WALL_VISCOUS_ISOTHERMAL ||
       SolnBlk.Grid.BCtypeW[SolnBlk.JCu] == BC_MOVING_WALL_HEATFLUX ||
       SolnBlk.Grid.BCtypeW[SolnBlk.JCu] == BC_MOVING_WALL_ISOTHERMAL ||
       SolnBlk.Grid.BCtypeW[SolnBlk.JCu] == BC_BURNING_SURFACE// ||
       //SolnBlk.Grid.BCtypeW[SolnBlk.JCu] == BC_MASS_INJECTION
       ) &&
      (SolnBlk.Grid.BCtypeN[SolnBlk.ICl] == BC_REFLECTION ||
       SolnBlk.Grid.BCtypeN[SolnBlk.ICl] == BC_WALL_VISCOUS_HEATFLUX ||
       SolnBlk.Grid.BCtypeN[SolnBlk.ICl] == BC_WALL_VISCOUS_ISOTHERMAL ||
       SolnBlk.Grid.BCtypeN[SolnBlk.ICl] == BC_MOVING_WALL_HEATFLUX ||
       SolnBlk.Grid.BCtypeN[SolnBlk.ICl] == BC_MOVING_WALL_ISOTHERMAL ||
       SolnBlk.Grid.BCtypeN[SolnBlk.ICl] == BC_BURNING_SURFACE// ||
       //SolnBlk.Grid.BCtypeN[SolnBlk.ICl] == BC_MASS_INJECTION
       )) {
    SolnBlk.W[SolnBlk.ICl-1][SolnBlk.JCu+1] = HALF*(SolnBlk.W[SolnBlk.ICl-1][SolnBlk.JCu]+
						    SolnBlk.W[SolnBlk.ICl][SolnBlk.JCu+1]);
    SolnBlk.U[SolnBlk.ICl-1][SolnBlk.JCu+1] = U(SolnBlk.W[SolnBlk.ICl-1][SolnBlk.JCu+1]);
    SolnBlk.W[SolnBlk.ICl-2][SolnBlk.JCu+1] = HALF*(SolnBlk.W[SolnBlk.ICl-2][SolnBlk.JCu]+
						    SolnBlk.W[SolnBlk.ICl][SolnBlk.JCu+1]);
    SolnBlk.U[SolnBlk.ICl-2][SolnBlk.JCu+1] = U(SolnBlk.W[SolnBlk.ICl-2][SolnBlk.JCu+1]);
    SolnBlk.W[SolnBlk.ICl-1][SolnBlk.JCu+2] = HALF*(SolnBlk.W[SolnBlk.ICl-1][SolnBlk.JCu]+
						    SolnBlk.W[SolnBlk.ICl][SolnBlk.JCu+2]);
    SolnBlk.U[SolnBlk.ICl-1][SolnBlk.JCu+2] = U(SolnBlk.W[SolnBlk.ICl-1][SolnBlk.JCu+2]);
    SolnBlk.W[SolnBlk.ICl-2][SolnBlk.JCu+2] = HALF*(SolnBlk.W[SolnBlk.ICl-2][SolnBlk.JCu]+
						    SolnBlk.W[SolnBlk.ICl][SolnBlk.JCu+2]);
    SolnBlk.U[SolnBlk.ICl-2][SolnBlk.JCu+2] = U(SolnBlk.W[SolnBlk.ICl-2][SolnBlk.JCu+2]);
  }

  if ((SolnBlk.Grid.BCtypeE[SolnBlk.JCl] == BC_REFLECTION ||
       SolnBlk.Grid.BCtypeE[SolnBlk.JCl] == BC_WALL_VISCOUS_HEATFLUX ||
       SolnBlk.Grid.BCtypeE[SolnBlk.JCl] == BC_WALL_VISCOUS_ISOTHERMAL ||
       SolnBlk.Grid.BCtypeE[SolnBlk.JCl] == BC_MOVING_WALL_HEATFLUX ||
       SolnBlk.Grid.BCtypeE[SolnBlk.JCl] == BC_MOVING_WALL_ISOTHERMAL ||
       SolnBlk.Grid.BCtypeE[SolnBlk.JCl] == BC_BURNING_SURFACE// ||
       //SolnBlk.Grid.BCtypeE[SolnBlk.JCl] == BC_MASS_INJECTION
       ) &&
      (SolnBlk.Grid.BCtypeS[SolnBlk.ICu] == BC_REFLECTION ||
       SolnBlk.Grid.BCtypeS[SolnBlk.ICu] == BC_WALL_VISCOUS_HEATFLUX ||
       SolnBlk.Grid.BCtypeS[SolnBlk.ICu] == BC_WALL_VISCOUS_ISOTHERMAL ||
       SolnBlk.Grid.BCtypeS[SolnBlk.ICu] == BC_MOVING_WALL_HEATFLUX ||
       SolnBlk.Grid.BCtypeS[SolnBlk.ICu] == BC_MOVING_WALL_ISOTHERMAL ||
       SolnBlk.Grid.BCtypeS[SolnBlk.ICu] == BC_BURNING_SURFACE// ||
       //SolnBlk.Grid.BCtypeS[SolnBlk.ICu] == BC_MASS_INJECTION
       )) {
    SolnBlk.W[SolnBlk.ICu+1][SolnBlk.JCl-1] = HALF*(SolnBlk.W[SolnBlk.ICu+1][SolnBlk.JCl]+
						    SolnBlk.W[SolnBlk.ICu][SolnBlk.JCl-1]);
    SolnBlk.U[SolnBlk.ICu+1][SolnBlk.JCl-1] = U(SolnBlk.W[SolnBlk.ICu+1][SolnBlk.JCl-1]);
    SolnBlk.W[SolnBlk.ICu+2][SolnBlk.JCl-1] = HALF*(SolnBlk.W[SolnBlk.ICu+2][SolnBlk.JCl  ] +
						    SolnBlk.W[SolnBlk.ICu  ][SolnBlk.JCl-1]);
    SolnBlk.U[SolnBlk.ICu+2][SolnBlk.JCl-1] = U(SolnBlk.W[SolnBlk.ICu+2][SolnBlk.JCl-1]);
    SolnBlk.W[SolnBlk.ICu+1][SolnBlk.JCl-2] = HALF*(SolnBlk.W[SolnBlk.ICu+1][SolnBlk.JCl  ] +
						    SolnBlk.W[SolnBlk.ICu  ][SolnBlk.JCl-2]);
    SolnBlk.U[SolnBlk.ICu+1][SolnBlk.JCl-2] = U(SolnBlk.W[SolnBlk.ICu+1][SolnBlk.JCl-2]);
    SolnBlk.W[SolnBlk.ICu+2][SolnBlk.JCl-2] = HALF*(SolnBlk.W[SolnBlk.ICu+2][SolnBlk.JCl  ] +
						    SolnBlk.W[SolnBlk.ICu  ][SolnBlk.JCl-2]);
    SolnBlk.U[SolnBlk.ICu+2][SolnBlk.JCl-2] = U(SolnBlk.W[SolnBlk.ICu+2][SolnBlk.JCl-2]);
  }

  if ((SolnBlk.Grid.BCtypeE[SolnBlk.JCu] == BC_REFLECTION ||
       SolnBlk.Grid.BCtypeE[SolnBlk.JCu] == BC_WALL_VISCOUS_HEATFLUX ||
       SolnBlk.Grid.BCtypeE[SolnBlk.JCu] == BC_WALL_VISCOUS_ISOTHERMAL ||
       SolnBlk.Grid.BCtypeE[SolnBlk.JCu] == BC_MOVING_WALL_HEATFLUX ||
       SolnBlk.Grid.BCtypeE[SolnBlk.JCu] == BC_MOVING_WALL_ISOTHERMAL ||
       SolnBlk.Grid.BCtypeE[SolnBlk.JCu] == BC_BURNING_SURFACE// ||
       //SolnBlk.Grid.BCtypeE[SolnBlk.JCu] == BC_MASS_INJECTION
       ) &&
      (SolnBlk.Grid.BCtypeN[SolnBlk.ICu] == BC_REFLECTION ||
       SolnBlk.Grid.BCtypeN[SolnBlk.ICu] == BC_WALL_VISCOUS_HEATFLUX ||
       SolnBlk.Grid.BCtypeN[SolnBlk.ICu] == BC_WALL_VISCOUS_ISOTHERMAL ||
       SolnBlk.Grid.BCtypeN[SolnBlk.ICu] == BC_MOVING_WALL_HEATFLUX ||
       SolnBlk.Grid.BCtypeN[SolnBlk.ICu] == BC_MOVING_WALL_ISOTHERMAL ||
       SolnBlk.Grid.BCtypeN[SolnBlk.ICu] == BC_BURNING_SURFACE// ||
       //SolnBlk.Grid.BCtypeN[SolnBlk.ICu] == BC_MASS_INJECTION
       )) {
    SolnBlk.W[SolnBlk.ICu+1][SolnBlk.JCu+1] = HALF*(SolnBlk.W[SolnBlk.ICu+1][SolnBlk.JCu  ]+
						    SolnBlk.W[SolnBlk.ICu  ][SolnBlk.JCu+1]);
    SolnBlk.U[SolnBlk.ICu+1][SolnBlk.JCu+1] = U(SolnBlk.W[SolnBlk.ICu+1][SolnBlk.JCu+1]);
    SolnBlk.W[SolnBlk.ICu+2][SolnBlk.JCu+1] = HALF*(SolnBlk.W[SolnBlk.ICu+2][SolnBlk.JCu  ] +
						    SolnBlk.W[SolnBlk.ICu  ][SolnBlk.JCu+1]);
    SolnBlk.U[SolnBlk.ICu+2][SolnBlk.JCu+1] = U(SolnBlk.W[SolnBlk.ICu+2][SolnBlk.JCu+1]);
    SolnBlk.W[SolnBlk.ICu+1][SolnBlk.JCu+2] = HALF*(SolnBlk.W[SolnBlk.ICu+1][SolnBlk.JCu  ] +
						    SolnBlk.W[SolnBlk.ICu  ][SolnBlk.JCu+2]);
    SolnBlk.U[SolnBlk.ICu+1][SolnBlk.JCu+2] = U(SolnBlk.W[SolnBlk.ICu+1][SolnBlk.JCu+2]);
    SolnBlk.W[SolnBlk.ICu+2][SolnBlk.JCu+2] = HALF*(SolnBlk.W[SolnBlk.ICu+2][SolnBlk.JCu  ] +
						    SolnBlk.W[SolnBlk.ICu  ][SolnBlk.JCu+2]);
    SolnBlk.U[SolnBlk.ICu+2][SolnBlk.JCu+2] = U(SolnBlk.W[SolnBlk.ICu+2][SolnBlk.JCu+2]);
  }

}

/**********************************************************************
 * Routine: CFL                                                       *
 *                                                                    *
 * Determines the allowable global and local time steps (for explicit *
 * Euler time stepping scheme) for the specified quadrilateral        *
 * solution block according to the Courant-Friedrichs-Lewy condition, *
 * the Neumann condition, and the particle relaxation time-scales.    *
 *                                                                    *
 **********************************************************************/
double CFL(NavierStokes2D_Quad_Block &SolnBlk,
	   NavierStokes2D_Input_Parameters &IP) {

  double dtMin = MILLION, d_i, d_j, v_i, v_j, a;
  double nu, dt_vis;

  for (int j = SolnBlk.JCl-SolnBlk.Nghost; j <= SolnBlk.JCu+SolnBlk.Nghost; j++) {
    for (int i = SolnBlk.ICl-SolnBlk.Nghost; i <= SolnBlk.ICu+SolnBlk.Nghost; i++) {
      if (i < SolnBlk.ICl || i > SolnBlk.ICu || j < SolnBlk.JCl || j > SolnBlk.JCu) {
	SolnBlk.dt[i][j] = ZERO;
      } else {

	d_i = TWO*(SolnBlk.Grid.Cell[i][j].A/
		   (SolnBlk.Grid.lfaceE(i,j)+SolnBlk.Grid.lfaceW(i,j)));
	d_j = TWO*(SolnBlk.Grid.Cell[i][j].A/
		   (SolnBlk.Grid.lfaceN(i,j)+SolnBlk.Grid.lfaceS(i,j)));

	v_i = HALF*((SolnBlk.W[i][j].v)*
		    (SolnBlk.Grid.nfaceE(i,j)-SolnBlk.Grid.nfaceW(i,j)));
	v_j = HALF*((SolnBlk.W[i][j].v)*
		    (SolnBlk.Grid.nfaceN(i,j)-SolnBlk.Grid.nfaceS(i,j)));

	// Inviscid dt calculation.
	a = SolnBlk.W[i][j].c();
	SolnBlk.dt[i][j] = min(d_i/(a+fabs(v_i)),d_j/(a+fabs(v_j)));

	// Viscous dt calculation.
	if (SolnBlk.Flow_Type) {  
	  nu = max(SolnBlk.W[i][j].nu(),SolnBlk.W[i][j].nuT());
	  dt_vis = HALF*min(d_i*d_i,d_j*d_j)/nu;
	  SolnBlk.dt[i][j] = min(SolnBlk.dt[i][j],dt_vis);
	}

	dtMin = min(dtMin,SolnBlk.dt[i][j]);

      }
    }
  }

  for (int j = SolnBlk.JCl-SolnBlk.Nghost; j <= SolnBlk.JCu+SolnBlk.Nghost; j++) {
    for (int i = SolnBlk.ICl-SolnBlk.Nghost; i <= SolnBlk.ICu+SolnBlk.Nghost; i++) {
      if (i < SolnBlk.ICl || i > SolnBlk.ICu || j < SolnBlk.JCl || j > SolnBlk.JCu) {
	SolnBlk.dt[i][j] = dtMin;
      }
    }
  }

  // Return the global time step.
  return dtMin;

}

/**********************************************************************
 * Routine: Set_Global_TimeStep                                       *
 *                                                                    *
 * Assigns global time step to specified solution block for           *
 * time-accurate calculations.                                        *
 *                                                                    *
 **********************************************************************/
void Set_Global_TimeStep(NavierStokes2D_Quad_Block &SolnBlk,
                         const double &Dt_min) {

  for (int j = SolnBlk.JCl; j <= SolnBlk.JCu; j++) {
    for (int i = SolnBlk.ICl; i <= SolnBlk.ICu; i++) {
      SolnBlk.dt[i][j] = Dt_min;
    }
  }

}

/**********************************************************************
 * Routine: L1_Norm_Residual                                          *
 *                                                                    *
 * Determines the L1-norm of the solution residual for the specified  *
 * quadrilateral solution block. Useful for monitoring convergence of *
 * the solution for steady state problems.                            *
 *                                                                    *
 **********************************************************************/
double L1_Norm_Residual(NavierStokes2D_Quad_Block &SolnBlk) {

  double l1_norm = ZERO;

  for (int j = SolnBlk.JCl; j <= SolnBlk.JCu; j++) {
    for (int i = SolnBlk.ICl; i <= SolnBlk.ICu; i++) {
      l1_norm += fabs(SolnBlk.dUdt[i][j][0][SolnBlk.residual_variable]);
    }
  }

  return l1_norm;

}

/**********************************************************************
 * Routine: L2_Norm_Residual                                          *
 *                                                                    *
 * Determines the L2-norm of the solution residual for the specified  *
 * quadrilateral solution block. Useful for monitoring convergence of *
 * the solution for steady state problems.                            *
 *                                                                    *
 **********************************************************************/
double L2_Norm_Residual(NavierStokes2D_Quad_Block &SolnBlk) {

  double l2_norm = ZERO;

  for (int j = SolnBlk.JCl; j <= SolnBlk.JCu; j++) {
    for (int i = SolnBlk.ICl; i <= SolnBlk.ICu; i++) {
      l2_norm += sqr(SolnBlk.dUdt[i][j][0][SolnBlk.residual_variable]);
    }
  }

  l2_norm = sqrt(l2_norm);

  return l2_norm;

}

/**********************************************************************
 * Routine: Max_Norm_Residual                                         *
 *                                                                    *
 * Determines the maximum norm of the solution residual for the       *
 * specified quadrilateral solution block.  Useful for monitoring     *
 * convergence of the solution for steady state problems.             *
 *                                                                    *
 **********************************************************************/
double Max_Norm_Residual(NavierStokes2D_Quad_Block &SolnBlk) {

  double max_norm = ZERO;

  for (int j = SolnBlk.JCl; j <= SolnBlk.JCu; j++) {
    for (int i = SolnBlk.ICl; i <= SolnBlk.ICu; i++) {
      max_norm = max(max_norm,fabs(SolnBlk.dUdt[i][j][0][SolnBlk.residual_variable]));
    }
  }

  return max_norm;

}

double L1_Norm_Residual(NavierStokes2D_Quad_Block &SolnBlk, int rv) {
  double l1_norm = ZERO;
  for (int j = SolnBlk.JCl; j <= SolnBlk.JCu; j++) {
    for (int i = SolnBlk.ICl; i <= SolnBlk.ICu; i++) {
      l1_norm += fabs(SolnBlk.dUdt[i][j][0][rv]);
    }
  }
  return l1_norm;
}

double L2_Norm_Residual(NavierStokes2D_Quad_Block &SolnBlk, int rv) {
  double l2_norm = ZERO;
  for (int j = SolnBlk.JCl; j <= SolnBlk.JCu; j++) {
    for (int i = SolnBlk.ICl; i <= SolnBlk.ICu; i++) {
      l2_norm += sqr(SolnBlk.dUdt[i][j][0][rv]);
    }
  }
  l2_norm = sqrt(l2_norm);
  return l2_norm;
}

double Max_Norm_Residual(NavierStokes2D_Quad_Block &SolnBlk, int rv) {
  double max_norm = ZERO;
  for (int j = SolnBlk.JCl; j <= SolnBlk.JCu; j++) {
    for (int i = SolnBlk.ICl; i <= SolnBlk.ICu; i++) {
      max_norm = max(max_norm, fabs(SolnBlk.dUdt[i][j][0][rv]));
    }
  }
  return max_norm;
}

/**********************************************************************
 * Routine: Linear_Reconstruction_GreenGauss                          *
 *                                                                    *
 * Peforms the reconstruction of a limited piecewise linear solution  *
 * state within a given cell (i,j) of the computational mesh for the  *
 * specified quadrilateral solution block.  A Green-Gauss approach is *
 * used in the evaluation of the unlimited solution gradients.        *
 * Several slope limiters may be used.                                *
 *                                                                    *
 **********************************************************************/
void Linear_Reconstruction_GreenGauss(NavierStokes2D_Quad_Block &SolnBlk,
				      const int i,
                                      const int j,
                                      const int Limiter) {

  int error_flag;
  int n_pts, i_index[8], j_index[8];
  double u0Min, u0Max, uQuad[4], phi;
  double DxDx_ave, DxDy_ave, DyDy_ave;
  double l_north, l_south, l_east, l_west;
  Vector2D n_north, n_south, n_east, n_west, dX;
  NavierStokes2D_pState W_nw, W_ne, W_sw, W_se, W_face,
                 DU, DUDx_ave, DUDy_ave;
  NavierStokes2D_pState NavierStokes2D_W_VACUUM; NavierStokes2D_W_VACUUM.Vacuum();

  // Carry out the limited solution reconstruction in each cell of 
  // the computational mesh.

  // Determine the number of neighbouring cells to be used in the
  // reconstruction procedure.  Away from boundaries this 8 neighbours
  // will be used.
  if (i == SolnBlk.ICl-SolnBlk.Nghost || i == SolnBlk.ICu+SolnBlk.Nghost ||
      j == SolnBlk.JCl-SolnBlk.Nghost || j == SolnBlk.JCu+SolnBlk.Nghost) {
    n_pts = 0;
  } else if ((i == SolnBlk.ICl-SolnBlk.Nghost+1) &&
	     (SolnBlk.Grid.BCtypeW[j] != BC_NONE)) {
    if (j == SolnBlk.JCl-SolnBlk.Nghost+1 || j == SolnBlk.JCu+SolnBlk.Nghost-1) {
      n_pts = 0;
    } else if (SolnBlk.Grid.BCtypeW[j] == BC_PERIODIC ||
	       SolnBlk.Grid.BCtypeW[j] == BC_CONSTANT_EXTRAPOLATION ||
	       SolnBlk.Grid.BCtypeW[j] == BC_LINEAR_EXTRAPOLATION ||
	       SolnBlk.Grid.BCtypeW[j] == BC_INFLOW_SUBSONIC ||
	       SolnBlk.Grid.BCtypeW[j] == BC_OUTFLOW_SUBSONIC ||
	       SolnBlk.Grid.BCtypeW[j] == BC_FIXED_PRESSURE ||
	       SolnBlk.Grid.BCtypeW[j] == BC_CHARACTERISTIC) {
      if (j == SolnBlk.JCl) {
	n_pts = 5;
	i_index[0] = i-1; j_index[0] = j  ;
	i_index[1] = i+1; j_index[1] = j  ;
	i_index[2] = i-1; j_index[2] = j+1;
	i_index[3] = i  ; j_index[3] = j+1;
	i_index[4] = i+1; j_index[4] = j+1;
      } else if (j == SolnBlk.JCu) {
	n_pts = 5;
	i_index[0] = i-1; j_index[0] = j-1;
	i_index[1] = i  ; j_index[1] = j-1;
	i_index[2] = i+1; j_index[2] = j-1;
	i_index[3] = i-1; j_index[3] = j  ;
	i_index[4] = i+1; j_index[4] = j  ;
      } else {
	n_pts = 8;
	i_index[0] = i-1; j_index[0] = j-1;
	i_index[1] = i  ; j_index[1] = j-1;
	i_index[2] = i+1; j_index[2] = j-1;
	i_index[3] = i-1; j_index[3] = j  ;
	i_index[4] = i+1; j_index[4] = j  ;
	i_index[5] = i-1; j_index[5] = j+1;
	i_index[6] = i  ; j_index[6] = j+1;
	i_index[7] = i+1; j_index[7] = j+1;
      }
    } else {
      if (j == SolnBlk.JCl) {
	n_pts = 3;
	i_index[0] = i+1; j_index[0] = j  ;
	i_index[1] = i  ; j_index[1] = j+1;
	i_index[2] = i+1; j_index[2] = j+1;
      } else if (j == SolnBlk.JCu) {
	n_pts = 3;
	i_index[0] = i  ; j_index[0] = j-1;
	i_index[1] = i+1; j_index[1] = j-1;
	i_index[2] = i+1; j_index[2] = j  ;
      } else {
	n_pts = 5;
	i_index[0] = i  ; j_index[0] = j-1;
	i_index[1] = i+1; j_index[1] = j-1;
	i_index[2] = i+1; j_index[2] = j  ;
	i_index[3] = i  ; j_index[3] = j+1;
	i_index[4] = i+1; j_index[4] = j+1;
      }
    }
  } else if ((i == SolnBlk.ICu+SolnBlk.Nghost-1) && 
	     (SolnBlk.Grid.BCtypeE[j] != BC_NONE)) {
    if (j == SolnBlk.JCl-SolnBlk.Nghost+1 || j == SolnBlk.JCu+SolnBlk.Nghost-1) {
      n_pts = 0;
    } else if (SolnBlk.Grid.BCtypeE[j] == BC_PERIODIC ||
	       SolnBlk.Grid.BCtypeE[j] == BC_CONSTANT_EXTRAPOLATION ||
	       SolnBlk.Grid.BCtypeE[j] == BC_LINEAR_EXTRAPOLATION ||
	       SolnBlk.Grid.BCtypeE[j] == BC_INFLOW_SUBSONIC || 
	       SolnBlk.Grid.BCtypeE[j] == BC_OUTFLOW_SUBSONIC || 
	       SolnBlk.Grid.BCtypeE[j] == BC_FIXED_PRESSURE ||
	       SolnBlk.Grid.BCtypeE[j] == BC_CHARACTERISTIC) {
      if (j == SolnBlk.JCl) {
	n_pts = 5;
	i_index[0] = i-1; j_index[0] = j  ;
	i_index[1] = i+1; j_index[1] = j  ;
	i_index[2] = i-1; j_index[2] = j+1;
	i_index[3] = i  ; j_index[3] = j+1;
	i_index[4] = i+1; j_index[4] = j+1;
      } else if (j == SolnBlk.JCu) {
	n_pts = 5;
	i_index[0] = i-1; j_index[0] = j-1;
	i_index[1] = i  ; j_index[1] = j-1;
	i_index[2] = i+1; j_index[2] = j-1;
	i_index[3] = i-1; j_index[3] = j  ;
	i_index[4] = i+1; j_index[4] = j  ;
      } else {
	n_pts = 8;
	i_index[0] = i-1; j_index[0] = j-1;
	i_index[1] = i  ; j_index[1] = j-1;
	i_index[2] = i+1; j_index[2] = j-1;
	i_index[3] = i-1; j_index[3] = j  ;
	i_index[4] = i+1; j_index[4] = j  ;
	i_index[5] = i-1; j_index[5] = j+1;
	i_index[6] = i  ; j_index[6] = j+1;
	i_index[7] = i+1; j_index[7] = j+1;
      }
    } else {
      if (j == SolnBlk.JCl) {
	n_pts = 3;
	i_index[0] = i-1; j_index[0] = j  ;
	i_index[1] = i-1; j_index[1] = j+1;
	i_index[2] = i  ; j_index[2] = j+1;
      } else if (j == SolnBlk.JCu) {
	n_pts = 3;
	i_index[0] = i-1; j_index[0] = j-1;
	i_index[1] = i  ; j_index[1] = j-1;
	i_index[2] = i-1; j_index[2] = j  ;
      } else {
	n_pts = 5;
	i_index[0] = i-1; j_index[0] = j-1;
	i_index[1] = i  ; j_index[1] = j-1;
	i_index[2] = i-1; j_index[2] = j  ;
	i_index[3] = i-1; j_index[3] = j+1;
	i_index[4] = i  ; j_index[4] = j+1;
      }
    }
  } else if ((j == SolnBlk.JCl-SolnBlk.Nghost+1) &&
	     (SolnBlk.Grid.BCtypeS[i] != BC_NONE)) {
    if (i == SolnBlk.ICl-SolnBlk.Nghost+1 || i == SolnBlk.ICu+SolnBlk.Nghost-1) {
      n_pts = 0;
    } else if (SolnBlk.Grid.BCtypeS[i] == BC_PERIODIC ||
	       SolnBlk.Grid.BCtypeS[i] == BC_CONSTANT_EXTRAPOLATION ||
	       SolnBlk.Grid.BCtypeS[i] == BC_LINEAR_EXTRAPOLATION ||
	       SolnBlk.Grid.BCtypeS[j] == BC_INFLOW_SUBSONIC || 
	       SolnBlk.Grid.BCtypeS[j] == BC_OUTFLOW_SUBSONIC || 
	       SolnBlk.Grid.BCtypeS[i] == BC_FIXED_PRESSURE ||
	       SolnBlk.Grid.BCtypeS[i] == BC_CHARACTERISTIC) {
      if (i == SolnBlk.ICl) {
	n_pts = 5;
	i_index[0] = i  ; j_index[0] = j-1;
	i_index[1] = i+1; j_index[1] = j-1;
	i_index[2] = i+1; j_index[2] = j  ;
	i_index[3] = i  ; j_index[3] = j+1;
	i_index[4] = i+1; j_index[4] = j+1;
      } else if (i == SolnBlk.ICu) {
	n_pts = 5;
	i_index[0] = i-1; j_index[0] = j-1;
	i_index[1] = i  ; j_index[1] = j-1;
	i_index[2] = i-1; j_index[2] = j  ;
	i_index[3] = i-1; j_index[3] = j+1;
	i_index[4] = i  ; j_index[4] = j+1;
      } else {
	n_pts = 8;
	i_index[0] = i-1; j_index[0] = j-1;
	i_index[1] = i  ; j_index[1] = j-1;
	i_index[2] = i+1; j_index[2] = j-1;
	i_index[3] = i-1; j_index[3] = j  ;
	i_index[4] = i+1; j_index[4] = j  ;
	i_index[5] = i-1; j_index[5] = j+1;
	i_index[6] = i  ; j_index[6] = j+1;
	i_index[7] = i+1; j_index[7] = j+1;
      }
    } else {
      if (i == SolnBlk.ICl) {
	n_pts = 3;
	i_index[0] = i+1; j_index[0] = j  ;
	i_index[1] = i  ; j_index[1] = j+1;
	i_index[2] = i+1; j_index[2] = j+1;
      } else if (i == SolnBlk.ICu) {
	n_pts = 3;
	i_index[0] = i-1; j_index[0] = j  ;
	i_index[1] = i-1; j_index[1] = j+1;
	i_index[2] = i  ; j_index[2] = j+1;
      } else {
	n_pts = 5;
	i_index[0] = i-1; j_index[0] = j  ;
	i_index[1] = i+1; j_index[1] = j  ;
	i_index[2] = i-1; j_index[2] = j+1;
	i_index[3] = i  ; j_index[3] = j+1;
	i_index[4] = i+1; j_index[4] = j+1;
      }
    }
  } else if ((j == SolnBlk.JCu+SolnBlk.Nghost-1) &&
 	     (SolnBlk.Grid.BCtypeN[i] != BC_NONE)) {
    if (i == SolnBlk.ICl-SolnBlk.Nghost+1 || i == SolnBlk.ICu+SolnBlk.Nghost-1) {
      n_pts = 0;
    } else if (SolnBlk.Grid.BCtypeN[i] == BC_PERIODIC ||
	       SolnBlk.Grid.BCtypeN[i] == BC_CONSTANT_EXTRAPOLATION ||
	       SolnBlk.Grid.BCtypeN[i] == BC_LINEAR_EXTRAPOLATION ||
	       SolnBlk.Grid.BCtypeN[j] == BC_INFLOW_SUBSONIC || 
	       SolnBlk.Grid.BCtypeN[j] == BC_OUTFLOW_SUBSONIC || 
	       SolnBlk.Grid.BCtypeN[i] == BC_FIXED_PRESSURE ||
	       SolnBlk.Grid.BCtypeN[i] == BC_CHARACTERISTIC) {
      if (i == SolnBlk.ICl) {
	n_pts = 5;
	i_index[0] = i  ; j_index[0] = j-1;
	i_index[1] = i+1; j_index[1] = j-1;
	i_index[2] = i+1; j_index[2] = j  ;
	i_index[3] = i  ; j_index[3] = j+1;
	i_index[4] = i+1; j_index[4] = j+1;
      } else if (i == SolnBlk.ICu) {
	n_pts = 5;
	i_index[0] = i-1; j_index[0] = j-1;
	i_index[1] = i  ; j_index[1] = j-1;
	i_index[2] = i-1; j_index[2] = j  ;
	i_index[3] = i-1; j_index[3] = j+1;
	i_index[4] = i  ; j_index[4] = j+1;
      } else {
	n_pts = 8;
	i_index[0] = i-1; j_index[0] = j-1;
	i_index[1] = i  ; j_index[1] = j-1;
	i_index[2] = i+1; j_index[2] = j-1;
	i_index[3] = i-1; j_index[3] = j  ;
	i_index[4] = i+1; j_index[4] = j  ;
	i_index[5] = i-1; j_index[5] = j+1;
	i_index[6] = i  ; j_index[6] = j+1;
	i_index[7] = i+1; j_index[7] = j+1;
      }
    } else {
      if (i == SolnBlk.ICl) {
	n_pts = 3;
	i_index[0] = i  ; j_index[0] = j-1;
	i_index[1] = i+1; j_index[1] = j-1;
	i_index[2] = i+1; j_index[2] = j  ;
      } else if (i == SolnBlk.ICu) {
	n_pts = 3;
	i_index[0] = i-1; j_index[0] = j-1;
	i_index[1] = i  ; j_index[1] = j-1;
	i_index[2] = i-1; j_index[2] = j  ;
      } else {
	n_pts = 5;
	i_index[0] = i-1; j_index[0] = j-1;
	i_index[1] = i  ; j_index[1] = j-1;
	i_index[2] = i+1; j_index[2] = j-1;
	i_index[3] = i-1; j_index[3] = j  ;
	i_index[4] = i+1; j_index[4] = j  ;
      }
    }
  } else {
    n_pts = 8;
    i_index[0] = i-1; j_index[0] = j-1;
    i_index[1] = i  ; j_index[1] = j-1;
    i_index[2] = i+1; j_index[2] = j-1;
    i_index[3] = i-1; j_index[3] = j  ;
    i_index[4] = i+1; j_index[4] = j  ;
    i_index[5] = i-1; j_index[5] = j+1;
    i_index[6] = i  ; j_index[6] = j+1;
    i_index[7] = i+1; j_index[7] = j+1;
  }

  // Perform reconstruction.    
  if (n_pts > 0) {
    // If 8 neighbours are used, apply Green-Gauss reconstruction
    if (n_pts == 8) {
      error_flag = Green_Gauss_Integration(SolnBlk.Grid.nodeNW(i,j).X,SolnBlk.WnNW(i,j),
					   SolnBlk.Grid.nodeSW(i,j).X,SolnBlk.WnSW(i,j),
					   SolnBlk.Grid.nodeSE(i,j).X,SolnBlk.WnSE(i,j),
					   SolnBlk.Grid.nodeNE(i,j).X,SolnBlk.WnNE(i,j),
					   SolnBlk.dWdx[i][j],SolnBlk.dWdy[i][j]);
      //if (error_flag) return error_flag;
//       W_nw = SolnBlk.WnNW(i,j);
//       W_ne = SolnBlk.WnNE(i,j);
//       W_sw = SolnBlk.WnSW(i,j);
//       W_se = SolnBlk.WnSE(i,j);

//       l_north = SolnBlk.Grid.lfaceN(i,j);
//       l_south = SolnBlk.Grid.lfaceS(i,j);
//       l_east = SolnBlk.Grid.lfaceE(i,j);
//       l_west = SolnBlk.Grid.lfaceW(i,j);

//       n_north = SolnBlk.Grid.nfaceN(i,j);
//       n_south = SolnBlk.Grid.nfaceS(i,j);
//       n_east = SolnBlk.Grid.nfaceE(i,j);
//       n_west = SolnBlk.Grid.nfaceW(i,j);

//       W_face = HALF*(W_nw+W_ne)*l_north; 
//       SolnBlk.dWdx[i][j] = W_face*n_north.x;
//       SolnBlk.dWdy[i][j] = W_face*n_north.y;

//       W_face = HALF*(W_sw+W_se)*l_south; 
//       SolnBlk.dWdx[i][j] += W_face*n_south.x;
//       SolnBlk.dWdy[i][j] += W_face*n_south.y;

//       W_face = HALF*(W_ne+W_se)*l_east; 
//       SolnBlk.dWdx[i][j] += W_face*n_east.x;
//       SolnBlk.dWdy[i][j] += W_face*n_east.y;

//       W_face = HALF*(W_nw+W_sw)*l_west; 
//       SolnBlk.dWdx[i][j] += W_face*n_west.x;
//       SolnBlk.dWdy[i][j] += W_face*n_west.y;

//       SolnBlk.dWdx[i][j] = SolnBlk.dWdx[i][j]/SolnBlk.Grid.Cell[i][j].A;
//       SolnBlk.dWdy[i][j] = SolnBlk.dWdy[i][j]/SolnBlk.Grid.Cell[i][j].A;

    } else {

      // If <8 neighbours are used, apply least-squares reconstruction
      DUDx_ave = NavierStokes2D_W_VACUUM;
      DUDy_ave = NavierStokes2D_W_VACUUM;
      DxDx_ave = ZERO;
      DxDy_ave = ZERO;
      DyDy_ave = ZERO;

      for (int n2 = 0; n2 < n_pts; n2++) {
	dX = SolnBlk.Grid.Cell[i_index[n2]][j_index[n2]].Xc - SolnBlk.Grid.Cell[i][j].Xc;
	DU = SolnBlk.W[i_index[n2]][j_index[n2]] - SolnBlk.W[i][j];
	DUDx_ave += DU*dX.x;
	DUDy_ave += DU*dX.y;
	DxDx_ave += dX.x*dX.x;
	DxDy_ave += dX.x*dX.y;
	DyDy_ave += dX.y*dX.y;
      }

      DUDx_ave /= double(n_pts);
      DUDy_ave /= double(n_pts);
      DxDx_ave /= double(n_pts);
      DxDy_ave /= double(n_pts);
      DyDy_ave /= double(n_pts);
      SolnBlk.dWdx[i][j] = (DUDx_ave*DyDy_ave-DUDy_ave*DxDy_ave)/
	                   (DxDx_ave*DyDy_ave-DxDy_ave*DxDy_ave);
      SolnBlk.dWdy[i][j] = (DUDy_ave*DxDx_ave-DUDx_ave*DxDy_ave)/
	                   (DxDx_ave*DyDy_ave-DxDy_ave*DxDy_ave);
    }

    // Calculate slope limiters.
    if (!SolnBlk.Freeze_Limiter) {
      for (int n = 1; n <= NUM_VAR_NAVIERSTOKES2D; n++) {
	u0Min = SolnBlk.W[i][j][n];
	u0Max = u0Min;
	for (int n2 = 0; n2 < n_pts; n2++) {
	  u0Min = min(u0Min,SolnBlk.W[i_index[n2]][j_index[n2]][n]);
	  u0Max = max(u0Max,SolnBlk.W[i_index[n2]][j_index[n2]][n]);
	}

	dX = SolnBlk.Grid.xfaceE(i,j) - SolnBlk.Grid.Cell[i][j].Xc;
	uQuad[0] = SolnBlk.W[i][j][n] + SolnBlk.dWdx[i][j][n]*dX.x +
                                        SolnBlk.dWdy[i][j][n]*dX.y;
	dX = SolnBlk.Grid.xfaceW(i,j) - SolnBlk.Grid.Cell[i][j].Xc;
	uQuad[1] = SolnBlk.W[i][j][n] + SolnBlk.dWdx[i][j][n]*dX.x +
                                        SolnBlk.dWdy[i][j][n]*dX.y;
	dX = SolnBlk.Grid.xfaceN(i,j) - SolnBlk.Grid.Cell[i][j].Xc;
	uQuad[2] = SolnBlk.W[i][j][n] + SolnBlk.dWdx[i][j][n]*dX.x +
	                                SolnBlk.dWdy[i][j][n]*dX.y;
	dX = SolnBlk.Grid.xfaceS(i,j) - SolnBlk.Grid.Cell[i][j].Xc;
	uQuad[3] = SolnBlk.W[i][j][n] + SolnBlk.dWdx[i][j][n]*dX.x +
                                        SolnBlk.dWdy[i][j][n]*dX.y;

	switch(Limiter) {
	case LIMITER_ONE :
	  phi = ONE;
	  break;
	case LIMITER_ZERO :
	  phi = ZERO;
	  break;
	case LIMITER_BARTH_JESPERSEN :
	  phi = Limiter_BarthJespersen(uQuad,SolnBlk.W[i][j][n],u0Min,u0Max,4);
	  break;
	case LIMITER_VENKATAKRISHNAN :
	  phi = Limiter_Venkatakrishnan(uQuad,SolnBlk.W[i][j][n],u0Min,u0Max,4);
	  break;
	case LIMITER_VANLEER :
	  phi = Limiter_VanLeer(uQuad,SolnBlk.W[i][j][n],u0Min,u0Max,4);
	  break;
	case LIMITER_VANALBADA :
	  phi = Limiter_VanAlbada(uQuad,SolnBlk.W[i][j][n],u0Min,u0Max,4);
	  break;
	default:
	  phi = Limiter_BarthJespersen(uQuad,SolnBlk.W[i][j][n],u0Min,u0Max,4);
	  break;
	};

	SolnBlk.phi[i][j][n] = phi;

      }
    }

  } else {
    SolnBlk.dWdx[i][j] = NavierStokes2D_W_VACUUM;
    SolnBlk.dWdy[i][j] = NavierStokes2D_W_VACUUM;
    SolnBlk.phi[i][j] = NavierStokes2D_W_VACUUM;
  }

}

/**********************************************************************
 * Routine: Linear_Reconstruction_GreenGauss                          *
 *                                                                    *
 * Peforms the reconstruction of a limited piecewise linear solution  *
 * state within each cell of the computational mesh for the specified *
 * quadrilateral solution block.  A Green-Gauss approach is used in   *
 * the evaluation of the unlimited solution gradients.  Several slope *
 * limiters may be used.                                              *
 *                                                                    *
 **********************************************************************/
void Linear_Reconstruction_GreenGauss(NavierStokes2D_Quad_Block &SolnBlk,
				      const int Limiter) {

  // Carry out the limited solution reconstruction in each cell of
  // the computational mesh.
  for (int j = SolnBlk.JCl-SolnBlk.Nghost+1; j <= SolnBlk.JCu+SolnBlk.Nghost-1; j++) {
    for (int i = SolnBlk.ICl-SolnBlk.Nghost+1; i <= SolnBlk.ICu+SolnBlk.Nghost-1; i++) {
      Linear_Reconstruction_GreenGauss(SolnBlk,i,j,Limiter);
    }
  }

}

/**********************************************************************
 * Routine: Linear_Reconstruction_LeastSquares                        *
 *                                                                    *
 * Peforms the reconstruction of a limited piecewise linear solution  *
 * state within a given cell (i,j) of the computational mesh for the  *
 * specified quadrilateral solution block.  A least squares approach  *
 * is used in the evaluation of the unlimited solution gradients.     *
 * Several slope limiters may be used.                                *
 *                                                                    *
 **********************************************************************/
void Linear_Reconstruction_LeastSquares(NavierStokes2D_Quad_Block &SolnBlk,
				        const int i,
                                        const int j,
                                        const int Limiter) {

  int n_pts, i_index[8], j_index[8], motion_flag = OFF, n_pts_temp;
  double u0Min, u0Max, uQuad[4], phi;
  double DxDx_ave, DxDy_ave, DyDy_ave, Dx_ave, Dy_ave;
  Vector2D dX;
  NavierStokes2D_pState DU, DUDx_ave, DUDy_ave;
  NavierStokes2D_pState NavierStokes2D_W_VACUUM; NavierStokes2D_W_VACUUM.Vacuum();

  // Carry out the limited solution reconstruction in each cell of the
  // computational mesh.

  if (i == SolnBlk.ICl-SolnBlk.Nghost || i == SolnBlk.ICu+SolnBlk.Nghost ||
      j == SolnBlk.JCl-SolnBlk.Nghost || j == SolnBlk.JCu+SolnBlk.Nghost) {
    n_pts = 0;
  } else if ((i == SolnBlk.ICl-SolnBlk.Nghost+1) && 
	     (SolnBlk.Grid.BCtypeW[j] != BC_NONE)) {
    if (j == SolnBlk.JCl-SolnBlk.Nghost+1 || j == SolnBlk.JCu+SolnBlk.Nghost-1) {
      n_pts = 0;
    } else if (SolnBlk.Grid.BCtypeW[j] == BC_PERIODIC ||
	       SolnBlk.Grid.BCtypeW[j] == BC_CONSTANT_EXTRAPOLATION ||
	       SolnBlk.Grid.BCtypeW[j] == BC_LINEAR_EXTRAPOLATION ||
	       SolnBlk.Grid.BCtypeW[j] == BC_INFLOW_SUBSONIC ||
	       SolnBlk.Grid.BCtypeW[j] == BC_OUTFLOW_SUBSONIC ||
	       SolnBlk.Grid.BCtypeW[j] == BC_FIXED_PRESSURE ||
	       SolnBlk.Grid.BCtypeW[j] == BC_CHARACTERISTIC) {
      if (j == SolnBlk.JCl) {
	n_pts = 5;
	i_index[0] = i-1; j_index[0] = j  ;
	i_index[1] = i+1; j_index[1] = j  ;
	i_index[2] = i-1; j_index[2] = j+1;
	i_index[3] = i  ; j_index[3] = j+1;
	i_index[4] = i+1; j_index[4] = j+1;
      } else if (j == SolnBlk.JCu) {
	n_pts = 5;
	i_index[0] = i-1; j_index[0] = j-1;
	i_index[1] = i  ; j_index[1] = j-1;
	i_index[2] = i+1; j_index[2] = j-1;
	i_index[3] = i-1; j_index[3] = j  ;
	i_index[4] = i+1; j_index[4] = j  ;
      } else {
	n_pts = 8;
	i_index[0] = i-1; j_index[0] = j-1;
	i_index[1] = i  ; j_index[1] = j-1;
	i_index[2] = i+1; j_index[2] = j-1;
	i_index[3] = i-1; j_index[3] = j  ;
	i_index[4] = i+1; j_index[4] = j  ;
	i_index[5] = i-1; j_index[5] = j+1;
	i_index[6] = i  ; j_index[6] = j+1;
	i_index[7] = i+1; j_index[7] = j+1;
      }
    } else {
      if (j == SolnBlk.JCl) {
	n_pts = 3;
	i_index[0] = i+1; j_index[0] = j  ;
	i_index[1] = i  ; j_index[1] = j+1;
	i_index[2] = i+1; j_index[2] = j+1;
      } else if (j == SolnBlk.JCu) {
	n_pts = 3;
	i_index[0] = i  ; j_index[0] = j-1;
	i_index[1] = i+1; j_index[1] = j-1;
	i_index[2] = i+1; j_index[2] = j  ;
      } else {
	n_pts = 5;
	i_index[0] = i  ; j_index[0] = j-1;
	i_index[1] = i+1; j_index[1] = j-1;
	i_index[2] = i+1; j_index[2] = j  ;
	i_index[3] = i  ; j_index[3] = j+1;
	i_index[4] = i+1; j_index[4] = j+1;
      }
    }
  } else if ((i == SolnBlk.ICu+SolnBlk.Nghost-1) && 
	     (SolnBlk.Grid.BCtypeE[j] != BC_NONE)) {
    if (j == SolnBlk.JCl-SolnBlk.Nghost+1 || j == SolnBlk.JCu+SolnBlk.Nghost-1) {
      n_pts = 0;
    } else if (SolnBlk.Grid.BCtypeE[j] == BC_PERIODIC ||
	       SolnBlk.Grid.BCtypeE[j] == BC_CONSTANT_EXTRAPOLATION ||
	       SolnBlk.Grid.BCtypeE[j] == BC_LINEAR_EXTRAPOLATION ||
	       SolnBlk.Grid.BCtypeE[j] == BC_INFLOW_SUBSONIC ||
	       SolnBlk.Grid.BCtypeE[j] == BC_OUTFLOW_SUBSONIC ||
	       SolnBlk.Grid.BCtypeE[j] == BC_FIXED_PRESSURE ||
	       SolnBlk.Grid.BCtypeE[j] == BC_CHARACTERISTIC) {
      if (j == SolnBlk.JCl) {
	n_pts = 5;
	i_index[0] = i-1; j_index[0] = j  ;
	i_index[1] = i+1; j_index[1] = j  ;
	i_index[2] = i-1; j_index[2] = j+1;
	i_index[3] = i  ; j_index[3] = j+1;
	i_index[4] = i+1; j_index[4] = j+1;
      } else if (j == SolnBlk.JCu) {
	n_pts = 5;
	i_index[0] = i-1; j_index[0] = j-1;
	i_index[1] = i  ; j_index[1] = j-1;
	i_index[2] = i+1; j_index[2] = j-1;
	i_index[3] = i-1; j_index[3] = j  ;
	i_index[4] = i+1; j_index[4] = j  ;
      } else {
	n_pts = 8;
	i_index[0] = i-1; j_index[0] = j-1;
	i_index[1] = i  ; j_index[1] = j-1;
	i_index[2] = i+1; j_index[2] = j-1;
	i_index[3] = i-1; j_index[3] = j  ;
	i_index[4] = i+1; j_index[4] = j  ;
	i_index[5] = i-1; j_index[5] = j+1;
	i_index[6] = i  ; j_index[6] = j+1;
	i_index[7] = i+1; j_index[7] = j+1;
      }
    } else {
      if (j == SolnBlk.JCl) {
	n_pts = 3;
	i_index[0] = i-1; j_index[0] = j  ;
	i_index[1] = i-1; j_index[1] = j+1;
	i_index[2] = i  ; j_index[2] = j+1;
      } else if (j == SolnBlk.JCu) {
	n_pts = 3;
	i_index[0] = i-1; j_index[0] = j-1;
	i_index[1] = i  ; j_index[1] = j-1;
	i_index[2] = i-1; j_index[2] = j  ;
      } else {
	n_pts = 5;
	i_index[0] = i-1; j_index[0] = j-1;
	i_index[1] = i  ; j_index[1] = j-1;
	i_index[2] = i-1; j_index[2] = j  ;
	i_index[3] = i-1; j_index[3] = j+1;
	i_index[4] = i  ; j_index[4] = j+1;
      }
    }
  } else if ((j == SolnBlk.JCl-SolnBlk.Nghost+1) && 
	     (SolnBlk.Grid.BCtypeS[i] != BC_NONE)) {
    if (i == SolnBlk.ICl-SolnBlk.Nghost+1 || i == SolnBlk.ICu+SolnBlk.Nghost-1) {
      n_pts = 0;
    } else if (SolnBlk.Grid.BCtypeS[i] == BC_PERIODIC ||
	       SolnBlk.Grid.BCtypeS[i] == BC_CONSTANT_EXTRAPOLATION ||
	       SolnBlk.Grid.BCtypeS[i] == BC_LINEAR_EXTRAPOLATION ||
	       SolnBlk.Grid.BCtypeS[i] == BC_INFLOW_SUBSONIC ||
	       SolnBlk.Grid.BCtypeS[i] == BC_OUTFLOW_SUBSONIC ||
	       SolnBlk.Grid.BCtypeS[i] == BC_FIXED_PRESSURE ||
	       SolnBlk.Grid.BCtypeS[i] == BC_CHARACTERISTIC) {
      if (i == SolnBlk.ICl) {
	n_pts = 5;
	i_index[0] = i  ; j_index[0] = j-1;
	i_index[1] = i+1; j_index[1] = j-1;
	i_index[2] = i+1; j_index[2] = j  ;
	i_index[3] = i  ; j_index[3] = j+1;
	i_index[4] = i+1; j_index[4] = j+1;
      } else if (i == SolnBlk.ICu) {
	n_pts = 5;
	i_index[0] = i-1; j_index[0] = j-1;
	i_index[1] = i  ; j_index[1] = j-1;
	i_index[2] = i-1; j_index[2] = j  ;
	i_index[3] = i-1; j_index[3] = j+1;
	i_index[4] = i  ; j_index[4] = j+1;
      } else {
	n_pts = 8;
	i_index[0] = i-1; j_index[0] = j-1;
	i_index[1] = i  ; j_index[1] = j-1;
	i_index[2] = i+1; j_index[2] = j-1;
	i_index[3] = i-1; j_index[3] = j  ;
	i_index[4] = i+1; j_index[4] = j  ;
	i_index[5] = i-1; j_index[5] = j+1;
	i_index[6] = i  ; j_index[6] = j+1;
	i_index[7] = i+1; j_index[7] = j+1;
      }
    } else {
      if (i == SolnBlk.ICl) {
	n_pts = 3;
	i_index[0] = i+1; j_index[0] = j  ;
	i_index[1] = i  ; j_index[1] = j+1;
	i_index[2] = i+1; j_index[2] = j+1;
      } else if (i == SolnBlk.ICu) {
	n_pts = 3;
	i_index[0] = i-1; j_index[0] = j  ;
	i_index[1] = i-1; j_index[1] = j+1;
	i_index[2] = i  ; j_index[2] = j+1;
      } else {
	n_pts = 5;
	i_index[0] = i-1; j_index[0] = j  ;
	i_index[1] = i+1; j_index[1] = j  ;
	i_index[2] = i-1; j_index[2] = j+1;
	i_index[3] = i  ; j_index[3] = j+1;
	i_index[4] = i+1; j_index[4] = j+1;
      }
    }
  } else if ((j == SolnBlk.JCu+SolnBlk.Nghost-1) &&
	     (SolnBlk.Grid.BCtypeN[i] != BC_NONE)) {
    if (i == SolnBlk.ICl-SolnBlk.Nghost+1 || i == SolnBlk.ICu+SolnBlk.Nghost-1) {
      n_pts = 0;
    } else if (SolnBlk.Grid.BCtypeN[i] == BC_PERIODIC ||
	       SolnBlk.Grid.BCtypeN[i] == BC_CONSTANT_EXTRAPOLATION ||
	       SolnBlk.Grid.BCtypeN[i] == BC_LINEAR_EXTRAPOLATION ||
	       SolnBlk.Grid.BCtypeN[i] == BC_INFLOW_SUBSONIC ||
	       SolnBlk.Grid.BCtypeN[i] == BC_OUTFLOW_SUBSONIC ||
	       SolnBlk.Grid.BCtypeN[i] == BC_FIXED_PRESSURE ||
	       SolnBlk.Grid.BCtypeN[i] == BC_CHARACTERISTIC) {
      if (i == SolnBlk.ICl) {
	n_pts = 5;
	i_index[0] = i  ; j_index[0] = j-1;
	i_index[1] = i+1; j_index[1] = j-1;
	i_index[2] = i+1; j_index[2] = j  ;
	i_index[3] = i  ; j_index[3] = j+1;
	i_index[4] = i+1; j_index[4] = j+1;
      } else if (i == SolnBlk.ICu) {
	n_pts = 5;
	i_index[0] = i-1; j_index[0] = j-1;
	i_index[1] = i  ; j_index[1] = j-1;
	i_index[2] = i-1; j_index[2] = j  ;
	i_index[3] = i-1; j_index[3] = j+1;
	i_index[4] = i  ; j_index[4] = j+1;
      } else {
	n_pts = 8;
	i_index[0] = i-1; j_index[0] = j-1;
	i_index[1] = i  ; j_index[1] = j-1;
	i_index[2] = i+1; j_index[2] = j-1;
	i_index[3] = i-1; j_index[3] = j  ;
	i_index[4] = i+1; j_index[4] = j  ;
	i_index[5] = i-1; j_index[5] = j+1;
	i_index[6] = i  ; j_index[6] = j+1;
	i_index[7] = i+1; j_index[7] = j+1;
      }
    } else {
      if (i == SolnBlk.ICl) {
	n_pts = 3;
	i_index[0] = i  ; j_index[0] = j-1;
	i_index[1] = i+1; j_index[1] = j-1;
	i_index[2] = i+1; j_index[2] = j  ;
      } else if (i == SolnBlk.ICu) {
	n_pts = 3;
	i_index[0] = i-1; j_index[0] = j-1;
	i_index[1] = i  ; j_index[1] = j-1;
	i_index[2] = i-1; j_index[2] = j  ;
      } else {
	n_pts = 5;
	i_index[0] = i-1; j_index[0] = j-1;
	i_index[1] = i  ; j_index[1] = j-1;
	i_index[2] = i+1; j_index[2] = j-1;
	i_index[3] = i-1; j_index[3] = j  ;
	i_index[4] = i+1; j_index[4] = j  ;
      }
    }
  } else {
    n_pts = 8;
    i_index[0] = i-1; j_index[0] = j-1;
    i_index[1] = i  ; j_index[1] = j-1;
    i_index[2] = i+1; j_index[2] = j-1;
    i_index[3] = i-1; j_index[3] = j  ;
    i_index[4] = i+1; j_index[4] = j  ;
    i_index[5] = i-1; j_index[5] = j+1;
    i_index[6] = i  ; j_index[6] = j+1;
    i_index[7] = i+1; j_index[7] = j+1;
  }

  // Perform the linear least squares reconstruction.
  if (n_pts > 0) {
    DUDx_ave = NavierStokes2D_W_VACUUM;
    DUDy_ave = NavierStokes2D_W_VACUUM;
    DxDx_ave = ZERO;
    DxDy_ave = ZERO;
    DyDy_ave = ZERO;
    Dx_ave = ZERO;
    Dy_ave = ZERO;
    
    for (int n2 = 0; n2 < n_pts; n2++) {
      dX = SolnBlk.Grid.Cell[i_index[n2]][j_index[n2]].Xc - SolnBlk.Grid.Cell[i][j].Xc;
      DU = SolnBlk.W[i_index[n2]][j_index[n2]] - SolnBlk.W[i][j];
      DUDx_ave += DU*dX.x;
      DUDy_ave += DU*dX.y;
      DxDx_ave += dX.x*dX.x;
      DxDy_ave += dX.x*dX.y;
      DyDy_ave += dX.y*dX.y;
      Dx_ave += dX.x;
      Dy_ave += dX.y;
    }

    DUDx_ave /= double(n_pts);
    DUDy_ave /= double(n_pts);
    DxDx_ave /= double(n_pts);
    DxDy_ave /= double(n_pts);
    DyDy_ave /= double(n_pts);
    Dx_ave /= double(n_pts);
    Dy_ave /= double(n_pts);
    SolnBlk.dWdx[i][j] = (DUDx_ave*DyDy_ave-DUDy_ave*DxDy_ave)/
                         (DxDx_ave*DyDy_ave-DxDy_ave*DxDy_ave);
    SolnBlk.dWdy[i][j] = (DUDy_ave*DxDx_ave-DUDx_ave*DxDy_ave)/
                         (DxDx_ave*DyDy_ave-DxDy_ave*DxDy_ave);

    /************************************************************************
     * Edited by Pradeep: 21st July, 2007                                   *
     * -----------------------------------                                  *
     *             d  / dW \                                                * 
     * d_dWdx_dW = -- | -- |  It returns the gradients of dw/dx w.r.t       *
     *             dW \ dx /  primitive variables for the solution block.   *
     *                                                                      *
     * Derivation in a nutshell:                                            * 
     * -------------------------                                            *
     *         ---                                                          * 
     * dW/dx = \   Coefficients*W_neighbours + Coeff_cell_centre*W_cell_cen * 
     *         /                                                            *  
     *         ---                                                          *  
     *  Therefore, d (dW/dx) =  Coeff_cell_center                           *   
     *             --                                                       *   
     *             dW => [remember, its cell centre W]                      *
     *                                                                      * 
     * For calculating the final expression for the Coeff_cell_center,      *
     * use the Least Square expression for dWdX and express DU=(U_ngbor-U_c)*
     * and the final expression can be derived.                             *
     ************************************************************************/ 
    
    SolnBlk.d_dWdx_dW[i][j][0] = (-Dx_ave*DyDy_ave+Dy_ave*DxDy_ave)/
                                 (DxDx_ave*DyDy_ave-DxDy_ave*DxDy_ave); 
    SolnBlk.d_dWdy_dW[i][j][0] = (-Dy_ave*DxDx_ave+Dx_ave*DxDy_ave)/
                                 (DxDx_ave*DyDy_ave-DxDy_ave*DxDy_ave); 

    // Calculate slope limiters.
    if (!SolnBlk.Freeze_Limiter) {
      for (int n = 1; n <= NUM_VAR_NAVIERSTOKES2D; n++) {
	u0Min = SolnBlk.W[i][j][n];
	u0Max = u0Min;
	for (int n2 = 0; n2 < n_pts; n2++) {
	  u0Min = min(u0Min,SolnBlk.W[i_index[n2]][j_index[n2]][n]);
	  u0Max = max(u0Max,SolnBlk.W[i_index[n2]][j_index[n2]][n]);
	}

	dX = SolnBlk.Grid.xfaceE(i,j) - SolnBlk.Grid.Cell[i][j].Xc;
	uQuad[0] = SolnBlk.W[i][j][n] + SolnBlk.dWdx[i][j][n]*dX.x +
	                                SolnBlk.dWdy[i][j][n]*dX.y;
	dX = SolnBlk.Grid.xfaceW(i,j) - SolnBlk.Grid.Cell[i][j].Xc;
	uQuad[1] = SolnBlk.W[i][j][n] + SolnBlk.dWdx[i][j][n]*dX.x +
	                                SolnBlk.dWdy[i][j][n]*dX.y;
	dX = SolnBlk.Grid.xfaceN(i,j) - SolnBlk.Grid.Cell[i][j].Xc;
	uQuad[2] = SolnBlk.W[i][j][n] + SolnBlk.dWdx[i][j][n]*dX.x +
                                        SolnBlk.dWdy[i][j][n]*dX.y;
	dX = SolnBlk.Grid.xfaceS(i,j) - SolnBlk.Grid.Cell[i][j].Xc;
	uQuad[3] = SolnBlk.W[i][j][n] + SolnBlk.dWdx[i][j][n]*dX.x +
                                        SolnBlk.dWdy[i][j][n]*dX.y;

	switch(Limiter) {
	case LIMITER_ONE :
	  phi = ONE;
	  break;
	case LIMITER_ZERO :
	  phi = ZERO;
	  break;
	case LIMITER_BARTH_JESPERSEN :
	  phi = Limiter_BarthJespersen(uQuad,SolnBlk.W[i][j][n],u0Min,u0Max,4);
	  break;
	case LIMITER_VENKATAKRISHNAN :
	  phi = Limiter_Venkatakrishnan(uQuad,SolnBlk.W[i][j][n],u0Min,u0Max,4);
	  break;
	case LIMITER_VANLEER :
	  phi = Limiter_VanLeer(uQuad,SolnBlk.W[i][j][n],u0Min,u0Max,4);
	  break;
	case LIMITER_VANALBADA :
	  phi = Limiter_VanAlbada(uQuad,SolnBlk.W[i][j][n],u0Min,u0Max,4);
	  break;
	default:
	  phi = Limiter_BarthJespersen(uQuad,SolnBlk.W[i][j][n],u0Min,u0Max,4);
	  break;
	};

	SolnBlk.phi[i][j][n] = phi;

      }
    }

  } else {
    SolnBlk.dWdx[i][j] = NavierStokes2D_W_VACUUM;
    SolnBlk.dWdy[i][j] = NavierStokes2D_W_VACUUM;
    SolnBlk.phi[i][j]  = NavierStokes2D_W_VACUUM;

  }

}

/**********************************************************************
 * Routine: Linear_Reconstruction_LeastSquares_2                      *
 *                                                                    *
 * Peforms the reconstruction of a limited piecewise linear solution  *
 * state within a given cell (i,j) of the computational mesh for the  *
 * specified quadrilateral solution block.  A least squares approach  *
 * is used in the evaluation of the unlimited solution gradients.     *
 * Several slope limiters may be used.                                *
 *                                                                    *
 **********************************************************************/
void Linear_Reconstruction_LeastSquares_2(NavierStokes2D_Quad_Block &SolnBlk,
				          const int i,
                                          const int j,
                                          const int Limiter) {

  int n_pts, i_index[8], j_index[8];
  double u0Min, u0Max, uQuad[4], phi;
  double DxDx_ave, DxDy_ave, DyDy_ave;
  Vector2D dX;
  NavierStokes2D_pState DU, DUDx_ave, DUDy_ave;
  NavierStokes2D_pState NavierStokes2D_W_VACUUM; NavierStokes2D_W_VACUUM.Vacuum();

  // Carry out the limited solution reconstruction in each cell of the
  // computational mesh.

  if (i == SolnBlk.ICl-SolnBlk.Nghost || i == SolnBlk.ICu+SolnBlk.Nghost ||
      j == SolnBlk.JCl-SolnBlk.Nghost || j == SolnBlk.JCu+SolnBlk.Nghost) {
    n_pts = 0;
  } else if ((i == SolnBlk.ICl-SolnBlk.Nghost+1) && 
	     (SolnBlk.Grid.BCtypeW[j] != BC_NONE)) {
    if (j == SolnBlk.JCl-SolnBlk.Nghost+1 || j == SolnBlk.JCu+SolnBlk.Nghost-1) {
      n_pts = 0;
    } else if (SolnBlk.Grid.BCtypeW[j] == BC_PERIODIC ||
	       SolnBlk.Grid.BCtypeW[j] == BC_CONSTANT_EXTRAPOLATION ||
	       SolnBlk.Grid.BCtypeW[j] == BC_LINEAR_EXTRAPOLATION ||
	       SolnBlk.Grid.BCtypeW[j] == BC_INFLOW_SUBSONIC || 
	       SolnBlk.Grid.BCtypeW[j] == BC_OUTFLOW_SUBSONIC || 
	       SolnBlk.Grid.BCtypeW[j] == BC_FIXED_PRESSURE ||
	       SolnBlk.Grid.BCtypeW[j] == BC_CHARACTERISTIC) {
      if (j == SolnBlk.JCl) {
	n_pts = 5;
	i_index[0] = i-1; j_index[0] = j  ;
	i_index[1] = i+1; j_index[1] = j  ;
	i_index[2] = i-1; j_index[2] = j+1;
	i_index[3] = i  ; j_index[3] = j+1;
	i_index[4] = i+1; j_index[4] = j+1;
      } else if (j == SolnBlk.JCu) {
	n_pts = 5;
	i_index[0] = i-1; j_index[0] = j-1;
	i_index[1] = i  ; j_index[1] = j-1;
	i_index[2] = i+1; j_index[2] = j-1;
	i_index[3] = i-1; j_index[3] = j  ;
	i_index[4] = i+1; j_index[4] = j  ;
      } else {
	n_pts = 8;
	i_index[0] = i-1; j_index[0] = j-1;
	i_index[1] = i  ; j_index[1] = j-1;
	i_index[2] = i+1; j_index[2] = j-1;
	i_index[3] = i-1; j_index[3] = j  ;
	i_index[4] = i+1; j_index[4] = j  ;
	i_index[5] = i-1; j_index[5] = j+1;
	i_index[6] = i  ; j_index[6] = j+1;
	i_index[7] = i+1; j_index[7] = j+1;
      }
    } else {
      if (j == SolnBlk.JCl) {
	n_pts = 0;
      } else if (j == SolnBlk.JCu) {
	n_pts = 0;
      } else {
	n_pts = 0;
      }
    }           
  } else if ((i == SolnBlk.ICu+SolnBlk.Nghost-1) && 
	     (SolnBlk.Grid.BCtypeE[j] != BC_NONE)) {
    if (j == SolnBlk.JCl-SolnBlk.Nghost+1 || j == SolnBlk.JCu+SolnBlk.Nghost-1) {
      n_pts = 0;
    } else if (SolnBlk.Grid.BCtypeE[j] == BC_PERIODIC ||
	       SolnBlk.Grid.BCtypeE[j] == BC_CONSTANT_EXTRAPOLATION ||
	       SolnBlk.Grid.BCtypeE[j] == BC_LINEAR_EXTRAPOLATION ||
	       SolnBlk.Grid.BCtypeE[j] == BC_INFLOW_SUBSONIC || 
	       SolnBlk.Grid.BCtypeE[j] == BC_OUTFLOW_SUBSONIC || 
	       SolnBlk.Grid.BCtypeE[j] == BC_FIXED_PRESSURE ||
	       SolnBlk.Grid.BCtypeE[j] == BC_CHARACTERISTIC) {
      if (j == SolnBlk.JCl) {
	n_pts = 5;
	i_index[0] = i-1; j_index[0] = j  ;
	i_index[1] = i+1; j_index[1] = j  ;
	i_index[2] = i-1; j_index[2] = j+1;
	i_index[3] = i  ; j_index[3] = j+1;
	i_index[4] = i+1; j_index[4] = j+1;
      } else if (j == SolnBlk.JCu) {
	n_pts = 5;
	i_index[0] = i-1; j_index[0] = j-1;
	i_index[1] = i  ; j_index[1] = j-1;
	i_index[2] = i+1; j_index[2] = j-1;
	i_index[3] = i-1; j_index[3] = j  ;
	i_index[4] = i+1; j_index[4] = j  ;
      } else {
	n_pts = 8;
	i_index[0] = i-1; j_index[0] = j-1;
	i_index[1] = i  ; j_index[1] = j-1;
	i_index[2] = i+1; j_index[2] = j-1;
	i_index[3] = i-1; j_index[3] = j  ;
	i_index[4] = i+1; j_index[4] = j  ;
	i_index[5] = i-1; j_index[5] = j+1;
	i_index[6] = i  ; j_index[6] = j+1;
	i_index[7] = i+1; j_index[7] = j+1;
      }
    } else {
      if (j == SolnBlk.JCl) {
	n_pts = 0;
      } else if (j == SolnBlk.JCu) {
	n_pts = 0;
      } else {
	n_pts = 0;
      }
    }
  } else if ((j == SolnBlk.JCl-SolnBlk.Nghost+1) && 
	     (SolnBlk.Grid.BCtypeS[i] != BC_NONE)) {
    if (i == SolnBlk.ICl-SolnBlk.Nghost+1 || i == SolnBlk.ICu+SolnBlk.Nghost-1) {
      n_pts = 0;
    } else if (SolnBlk.Grid.BCtypeS[i] == BC_PERIODIC ||
	       SolnBlk.Grid.BCtypeS[i] == BC_CONSTANT_EXTRAPOLATION ||
	       SolnBlk.Grid.BCtypeS[i] == BC_LINEAR_EXTRAPOLATION ||
	       SolnBlk.Grid.BCtypeS[j] == BC_INFLOW_SUBSONIC || 
	       SolnBlk.Grid.BCtypeS[j] == BC_OUTFLOW_SUBSONIC || 
	       SolnBlk.Grid.BCtypeS[i] == BC_FIXED_PRESSURE ||
	       SolnBlk.Grid.BCtypeS[i] == BC_CHARACTERISTIC) {
      if (i == SolnBlk.ICl) {
	n_pts = 5;
	i_index[0] = i  ; j_index[0] = j-1;
	i_index[1] = i+1; j_index[1] = j-1;
	i_index[2] = i+1; j_index[2] = j  ;
	i_index[3] = i  ; j_index[3] = j+1;
	i_index[4] = i+1; j_index[4] = j+1;
      } else if (i == SolnBlk.ICu) {
	n_pts = 5;
	i_index[0] = i-1; j_index[0] = j-1;
	i_index[1] = i  ; j_index[1] = j-1;
	i_index[2] = i-1; j_index[2] = j  ;
	i_index[3] = i-1; j_index[3] = j+1;
	i_index[4] = i  ; j_index[4] = j+1;
      } else {
	n_pts = 8;
	i_index[0] = i-1; j_index[0] = j-1;
	i_index[1] = i  ; j_index[1] = j-1;
	i_index[2] = i+1; j_index[2] = j-1;
	i_index[3] = i-1; j_index[3] = j  ;
	i_index[4] = i+1; j_index[4] = j  ;
	i_index[5] = i-1; j_index[5] = j+1;
	i_index[6] = i  ; j_index[6] = j+1;
	i_index[7] = i+1; j_index[7] = j+1;
      }
    } else {
      if (i == SolnBlk.ICl) {
	n_pts = 0;
      } else if (i == SolnBlk.ICu) {
	n_pts = 0;
      } else {
	n_pts = 0;
      }
    }
  } else if ((j == SolnBlk.JCu+SolnBlk.Nghost-1) && 
	     (SolnBlk.Grid.BCtypeN[i] != BC_NONE)) {
    if (i == SolnBlk.ICl-SolnBlk.Nghost+1 || i == SolnBlk.ICu+SolnBlk.Nghost-1) {
      n_pts = 0;
    } else if (SolnBlk.Grid.BCtypeN[i] == BC_PERIODIC ||
	       SolnBlk.Grid.BCtypeN[i] == BC_CONSTANT_EXTRAPOLATION ||
	       SolnBlk.Grid.BCtypeN[i] == BC_LINEAR_EXTRAPOLATION ||
	       SolnBlk.Grid.BCtypeN[j] == BC_INFLOW_SUBSONIC || 
	       SolnBlk.Grid.BCtypeN[j] == BC_OUTFLOW_SUBSONIC || 
	       SolnBlk.Grid.BCtypeN[i] == BC_FIXED_PRESSURE ||
	       SolnBlk.Grid.BCtypeN[i] == BC_CHARACTERISTIC) {
      if (i == SolnBlk.ICl) {
	n_pts = 5;
	i_index[0] = i  ; j_index[0] = j-1;
	i_index[1] = i+1; j_index[1] = j-1;
	i_index[2] = i+1; j_index[2] = j  ;
	i_index[3] = i  ; j_index[3] = j+1;
	i_index[4] = i+1; j_index[4] = j+1;
      } else if (i == SolnBlk.ICu) {
	n_pts = 5;
	i_index[0] = i-1; j_index[0] = j-1;
	i_index[1] = i  ; j_index[1] = j-1;
	i_index[2] = i-1; j_index[2] = j  ;
	i_index[3] = i-1; j_index[3] = j+1;
	i_index[4] = i  ; j_index[4] = j+1;
      } else {
	n_pts = 8;
	i_index[0] = i-1; j_index[0] = j-1;
	i_index[1] = i  ; j_index[1] = j-1;
	i_index[2] = i+1; j_index[2] = j-1;
	i_index[3] = i-1; j_index[3] = j  ;
	i_index[4] = i+1; j_index[4] = j  ;
	i_index[5] = i-1; j_index[5] = j+1;
	i_index[6] = i  ; j_index[6] = j+1;
	i_index[7] = i+1; j_index[7] = j+1;
      }
    } else {
      if (i == SolnBlk.ICl) {
	n_pts = 0;
      } else if (i == SolnBlk.ICu) {
	n_pts = 0;
      } else {
	n_pts = 0;
      }
    }
  } else if ((i == SolnBlk.ICl && j == SolnBlk.JCl) && 
	     (SolnBlk.Grid.BCtypeW[j] != BC_NONE &&
	      SolnBlk.Grid.BCtypeS[i] != BC_NONE)) {
    n_pts = 3;
    i_index[0] = i+1; j_index[0] = j  ;
    i_index[1] = i  ; j_index[1] = j+1;
    i_index[2] = i+1; j_index[2] = j+1;
  } else if ((i == SolnBlk.ICl && j == SolnBlk.JCu) && 
	     (SolnBlk.Grid.BCtypeW[j] != BC_NONE &&
	      SolnBlk.Grid.BCtypeN[i] != BC_NONE)) {
    n_pts = 3;
    i_index[0] = i  ; j_index[0] = j-1;
    i_index[1] = i+1; j_index[1] = j-1;
    i_index[2] = i+1; j_index[2] = j  ;
  } else if ((i == SolnBlk.ICu && j == SolnBlk.JCl) && 
	     (SolnBlk.Grid.BCtypeE[j] != BC_NONE &&
	      SolnBlk.Grid.BCtypeS[i] != BC_NONE)) {
    n_pts = 3;
    i_index[0] = i-1; j_index[0] = j  ;
    i_index[1] = i-1; j_index[1] = j+1;
    i_index[2] = i  ; j_index[2] = j+1;
  } else if ((i == SolnBlk.ICu && j == SolnBlk.JCu) && 
	     (SolnBlk.Grid.BCtypeE[j] != BC_NONE &&
	      SolnBlk.Grid.BCtypeN[i] != BC_NONE)) {
    n_pts = 3;
    i_index[0] = i-1; j_index[0] = j-1;
    i_index[1] = i  ; j_index[1] = j-1;
    i_index[2] = i-1; j_index[2] = j  ;
  } else if ((i == SolnBlk.ICl) && 
	     (SolnBlk.Grid.BCtypeW[j] != BC_NONE)) {
    n_pts = 5;
    i_index[0] = i  ; j_index[0] = j-1;
    i_index[1] = i+1; j_index[1] = j-1;
    i_index[2] = i+1; j_index[2] = j  ;
    i_index[3] = i  ; j_index[3] = j+1;
    i_index[4] = i+1; j_index[4] = j+1;
  } else if ((j == SolnBlk.JCu) && 
	     (SolnBlk.Grid.BCtypeN[i] != BC_NONE)) {
    n_pts = 5;
    i_index[0] = i-1; j_index[0] = j-1;
    i_index[1] = i  ; j_index[1] = j-1;
    i_index[2] = i+1; j_index[2] = j-1;
    i_index[3] = i-1; j_index[3] = j  ;
    i_index[4] = i+1; j_index[4] = j  ;
  } else if ((i == SolnBlk.ICu) && 
	     (SolnBlk.Grid.BCtypeE[j] != BC_NONE)) {
    n_pts = 5;
    i_index[0] = i-1; j_index[0] = j-1;
    i_index[1] = i  ; j_index[1] = j-1;
    i_index[2] = i-1; j_index[2] = j  ;
    i_index[3] = i-1; j_index[3] = j+1;
    i_index[4] = i  ; j_index[4] = j+1;
    
  } else if ((j == SolnBlk.JCl) && 
	     (SolnBlk.Grid.BCtypeS[i] != BC_NONE)) {
    n_pts = 5;
    i_index[0] = i-1; j_index[0] = j  ;
    i_index[1] = i+1; j_index[1] = j  ;
    i_index[2] = i-1; j_index[2] = j+1;
    i_index[3] = i  ; j_index[3] = j+1;
    i_index[4] = i+1; j_index[4] = j+1;
  } else {
    n_pts = 8;
    i_index[0] = i-1; j_index[0] = j-1;
    i_index[1] = i  ; j_index[1] = j-1;
    i_index[2] = i+1; j_index[2] = j-1;
    i_index[3] = i-1; j_index[3] = j  ;
    i_index[4] = i+1; j_index[4] = j  ;
    i_index[5] = i-1; j_index[5] = j+1;
    i_index[6] = i  ; j_index[6] = j+1;
    i_index[7] = i+1; j_index[7] = j+1;
  }

  if (n_pts > 0) {
    DUDx_ave = NavierStokes2D_W_VACUUM;
    DUDy_ave = NavierStokes2D_W_VACUUM;
    DxDx_ave = ZERO;
    DxDy_ave = ZERO;
    DyDy_ave = ZERO;

    for (int n2 = 0; n2 < n_pts; n2++) {
      dX = SolnBlk.Grid.Cell[i_index[n2]][j_index[n2]].Xc - SolnBlk.Grid.Cell[i][j].Xc;
      DU = SolnBlk.W[i_index[n2]][j_index[n2]] - SolnBlk.W[i][j];
      DUDx_ave += DU*dX.x;
      DUDy_ave += DU*dX.y;
      DxDx_ave += dX.x*dX.x;
      DxDy_ave += dX.x*dX.y;
      DyDy_ave += dX.y*dX.y;
    }

    DUDx_ave /= double(n_pts);
    DUDy_ave /= double(n_pts);
    DxDx_ave /= double(n_pts);
    DxDy_ave /= double(n_pts);
    DyDy_ave /= double(n_pts);
    SolnBlk.dWdx[i][j] = (DUDx_ave*DyDy_ave-DUDy_ave*DxDy_ave)/
                         (DxDx_ave*DyDy_ave-DxDy_ave*DxDy_ave);
    SolnBlk.dWdy[i][j] = (DUDy_ave*DxDx_ave-DUDx_ave*DxDy_ave)/
	                 (DxDx_ave*DyDy_ave-DxDy_ave*DxDy_ave);

    // Calculate slope limiters.
    if (!SolnBlk.Freeze_Limiter) {
      for (int n = 1; n <= NUM_VAR_NAVIERSTOKES2D; n++) {
	u0Min = SolnBlk.W[i][j][n];
	u0Max = u0Min;
	for (int n2 = 0; n2 < n_pts; n2++) {
	  u0Min = min(u0Min,SolnBlk.W[i_index[n2]][j_index[n2]][n]);
	  u0Max = max(u0Max,SolnBlk.W[i_index[n2]][j_index[n2]][n]);
	}

	dX = SolnBlk.Grid.xfaceE(i,j) - SolnBlk.Grid.Cell[i][j].Xc;
	uQuad[0] = SolnBlk.W[i][j][n] + SolnBlk.dWdx[i][j][n]*dX.x +
                                        SolnBlk.dWdy[i][j][n]*dX.y;
	dX = SolnBlk.Grid.xfaceW(i,j) - SolnBlk.Grid.Cell[i][j].Xc;
	uQuad[1] = SolnBlk.W[i][j][n] + SolnBlk.dWdx[i][j][n]*dX.x +
                                        SolnBlk.dWdy[i][j][n]*dX.y;
	dX = SolnBlk.Grid.xfaceN(i,j) - SolnBlk.Grid.Cell[i][j].Xc;
	uQuad[2] = SolnBlk.W[i][j][n] + SolnBlk.dWdx[i][j][n]*dX.x +
                                        SolnBlk.dWdy[i][j][n]*dX.y;
	dX = SolnBlk.Grid.xfaceS(i,j) - SolnBlk.Grid.Cell[i][j].Xc;
	uQuad[3] = SolnBlk.W[i][j][n] + SolnBlk.dWdx[i][j][n]*dX.x +
                                        SolnBlk.dWdy[i][j][n]*dX.y;

	switch(Limiter) {
	case LIMITER_ONE :
	  phi = ONE;
	  break;
	case LIMITER_ZERO :
	  phi = ZERO;
	  break;
	case LIMITER_BARTH_JESPERSEN :
	  phi = Limiter_BarthJespersen(uQuad,SolnBlk.W[i][j][n],u0Min,u0Max,4);
	  break;
	case LIMITER_VENKATAKRISHNAN :
	  phi = Limiter_Venkatakrishnan(uQuad,SolnBlk.W[i][j][n],u0Min,u0Max,4);
	  break;
	case LIMITER_VANLEER :
	  phi = Limiter_VanLeer(uQuad,SolnBlk.W[i][j][n],u0Min,u0Max,4);
	  break;
	case LIMITER_VANALBADA :
	  phi = Limiter_VanAlbada(uQuad,SolnBlk.W[i][j][n],u0Min,u0Max,4);
	  break;
	default:
	  phi = Limiter_BarthJespersen(uQuad,SolnBlk.W[i][j][n],u0Min,u0Max,4);
	  break;
	};
     
	SolnBlk.phi[i][j][n] = phi;

      }
    }

  } else {
    SolnBlk.dWdx[i][j] = NavierStokes2D_W_VACUUM;
    SolnBlk.dWdy[i][j] = NavierStokes2D_W_VACUUM;
    SolnBlk.phi[i][j] = NavierStokes2D_W_VACUUM;

  }

}

/**********************************************************************
 * Routine: Linear_Reconstruction_LeastSquares                        *
 *                                                                    *         
 * Peforms the reconstruction of a limited piecewise linear solution  *
 * state within each cell of the computational mesh of the specified  *
 * quadrilateral solution block.  A least squares approach is used in *
 * the evaluation of the unlimited solution gradients.  Several slope *
 * limiters may be used.                                              *
 *                                                                    *
 **********************************************************************/
void Linear_Reconstruction_LeastSquares(NavierStokes2D_Quad_Block &SolnBlk,
				        const int Limiter) {

  // Carry out the limited solution reconstruction in each cell of the 
  // computational mesh.
  for (int j = SolnBlk.JCl-SolnBlk.Nghost+1; j <= SolnBlk.JCu+SolnBlk.Nghost-1; j++) {
    for (int i = SolnBlk.ICl-SolnBlk.Nghost+1; i <= SolnBlk.ICu+SolnBlk.Nghost-1; i++) {
      Linear_Reconstruction_LeastSquares(SolnBlk,i,j,Limiter);
      //Linear_Reconstruction_LeastSquares_2(SolnBlk,i,j,Limiter);
    }
  }

}

/**********************************************************************
 * Routine: Residual_Smoothing                                        *
 *                                                                    *
 * Applies implicit residual smoothing to solution residual.  Note    *
 * that only residuals of interior cells are smoothed and residuals   *
 * for cells adjacent to boundaries are not smoothed.                 *
 *                                                                    *
 **********************************************************************/
void Residual_Smoothing(NavierStokes2D_Quad_Block &SolnBlk,
                        const int k_residual,
			double &epsilon,
                        const int number_of_Gauss_Seidel_iterations) {

  for (int j = SolnBlk.JCl+1; j <= SolnBlk.JCu-1; j++) {
    for (int i = SolnBlk.ICl+1; i <= SolnBlk.ICu-1; i++) {
      SolnBlk.dUdt[i][j][k_residual+1] = SolnBlk.dUdt[i][j][k_residual];
    }
  }

  for (int n = 0; n < number_of_Gauss_Seidel_iterations; n++) {
    for (int j = SolnBlk.JCl+1; j <= SolnBlk.JCu-1; j++) {
      for (int i = SolnBlk.ICl+1; i <= SolnBlk.ICu-1; i++) {
	SolnBlk.dUdt[i][j][k_residual+1] = (SolnBlk.dUdt[i][j][k_residual] +
					    epsilon*(SolnBlk.dUdt[i  ][j-1][k_residual+1] +
						     SolnBlk.dUdt[i-1][j  ][k_residual+1] +
						     SolnBlk.dUdt[i+1][j  ][k_residual+1] +
						     SolnBlk.dUdt[i  ][j+1][k_residual+1]))/(ONE + FOUR*epsilon);
      }
    }
  }

  for (int j = SolnBlk.JCl+1; j <= SolnBlk.JCu-1; j++) {
    for (int i = SolnBlk.ICl+1; i <= SolnBlk.ICu-1; i++) {
      SolnBlk.dUdt[i][j][k_residual] = SolnBlk.dUdt[i][j][k_residual+1];
    }
  }

}

/**********************************************************************
 * Routine: Calculate_Refinement_Criteria                             *
 *                                                                    *
 * Calculate refinement criteria for the solution block.              *
 *                                                                    *
 **********************************************************************/
void Calculate_Refinement_Criteria(double *refinement_criteria,
				   NavierStokes2D_Input_Parameters &IP,
                                   int &number_refinement_criteria,
                                   NavierStokes2D_Quad_Block &SolnBlk) {

  double grad_rho_x, grad_rho_y, grad_rho_abs, grad_rho_criteria, grad_rho_criteria_max,
         div_V, div_V_criteria, div_V_criteria_max,
         curl_V_z, curl_V_abs, curl_V_criteria, curl_V_criteria_max,
         grad_k_x, grad_k_y, grad_k_abs, grad_k_criteria, grad_k_criteria_max;
  int refinement_criteria_number;

  // Set the number of refinement criteria to be used:
  // (1) Refinement criteria based on the gradient of the density field;
  // (2) Refinement criteria based on the divergence of the velocity vector;
  // (3) Refinement criteria based on the curl of the velocity vector;
  // (4) Refinement criteria based on the gradient of the turbulence kinetic energy field.
  number_refinement_criteria = IP.Number_of_Refinement_Criteria;

  // Initialize the refinement criteria for the block.
  grad_rho_criteria_max = ZERO;
  div_V_criteria_max = ZERO;
  curl_V_criteria_max = ZERO;
  grad_k_criteria_max = ZERO;

  // Calculate the refinement criteria for each cell of the 
  // computational mesh and assign the maximum value for all cells as 
  // the refinement criteria for the solution block.
  for (int j = SolnBlk.JCl-1; j <= SolnBlk.JCu+1; j++) {
    for (int i = SolnBlk.ICl-1; i <= SolnBlk.ICu+1; i++) {

      if ((j == SolnBlk.JCu+1 && SolnBlk.Grid.BCtypeN[i] != BC_NONE) ||
	  (j == SolnBlk.JCl-1 && SolnBlk.Grid.BCtypeS[i] != BC_NONE) ||
	  (i == SolnBlk.ICu+1 && SolnBlk.Grid.BCtypeE[j] != BC_NONE) ||
	  (i == SolnBlk.ICl-1 && SolnBlk.Grid.BCtypeW[j] != BC_NONE)) {

      } else {

	// Reconstruct the solution within the cell.
	//Linear_Reconstruction_GreenGauss(SolnBlk,i,j,LIMITER_UNLIMITED);
	Linear_Reconstruction_LeastSquares(SolnBlk,i,j,LIMITER_UNLIMITED);
	//Linear_Reconstruction_LeastSquares_2(SolnBlk,i,j,LIMITER_UNLIMITED);

	// Evaluate refinement criteria #1 based on the gradient of the 
	// density field.
	if (IP.Refinement_Criteria_Gradient_Density) {
	  grad_rho_x = SolnBlk.dWdx[i][j].rho;
	  grad_rho_y = SolnBlk.dWdy[i][j].rho;
	  grad_rho_abs = sqrt(sqr(grad_rho_x) + sqr(grad_rho_y));
	  grad_rho_criteria = sqrt(SolnBlk.Grid.Cell[i][j].A)*grad_rho_abs/SolnBlk.W[i][j].rho;
	} else {
	  grad_rho_criteria = ONE;
	}
	grad_rho_criteria_max = max(grad_rho_criteria_max,grad_rho_criteria);

	// Evaluate refinement criteria #2 based on the divergence of the
	// velocity vector.
	if (IP.Refinement_Criteria_Divergence_Velocity) {
	  div_V = SolnBlk.dWdx[i][j].v.x + SolnBlk.dWdy[i][j].v.y;
	  div_V_criteria = sqrt(SolnBlk.Grid.Cell[i][j].A)*fabs(div_V)/SolnBlk.W[i][j].a();
	} else {
	  div_V_criteria = ONE;
	}
	div_V_criteria_max = max(div_V_criteria_max,div_V_criteria);

	// Evaluate refinement criteria #3 based on the curl of the
	// velocity vector.
	if (IP.Refinement_Criteria_Curl_Velocity) {
	  curl_V_z = SolnBlk.dWdx[i][j].v.y - SolnBlk.dWdy[i][j].v.x; 
	  curl_V_abs = sqrt(sqr(curl_V_z)); 
	  curl_V_criteria = sqrt(SolnBlk.Grid.Cell[i][j].A)*curl_V_abs/SolnBlk.W[i][j].a();
	} else {
	  curl_V_criteria = ONE;
	}
	curl_V_criteria_max = max(curl_V_criteria_max,curl_V_criteria);

	// Evaluate refinement criteria #1 based on the gradient of the 
	// turbulence kinetic energy field.
	if (IP.Refinement_Criteria_Gradient_Turbulence_Kinetic_Energy) {
	  grad_k_x = SolnBlk.dWdx[i][j].k;
	  grad_k_y = SolnBlk.dWdy[i][j].k;
	  grad_k_abs = sqrt(sqr(grad_k_x) + sqr(grad_k_y));
	  grad_k_criteria = sqrt(SolnBlk.Grid.Cell[i][j].A)*grad_k_abs/SolnBlk.W[i][j].k;
	} else {
	  grad_k_criteria = ONE;
	}
	grad_k_criteria_max = max(grad_k_criteria_max,grad_k_criteria);

      }

    }
  }

  // Return the refinement criteria.
  refinement_criteria_number = 0;
  if (IP.Refinement_Criteria_Gradient_Density) {
    refinement_criteria[refinement_criteria_number] = grad_rho_criteria_max;
    refinement_criteria_number++;
  }
  if (IP.Refinement_Criteria_Divergence_Velocity) {
    refinement_criteria[refinement_criteria_number] = div_V_criteria_max;
    refinement_criteria_number++;
  }
  if (IP.Refinement_Criteria_Curl_Velocity) {
    refinement_criteria[refinement_criteria_number] = curl_V_criteria_max;
    refinement_criteria_number++;
  }
  if (IP.Refinement_Criteria_Gradient_Turbulence_Kinetic_Energy) {
    refinement_criteria[refinement_criteria_number] = grad_k_criteria_max;
    refinement_criteria_number++;
  }

}

/**********************************************************************
 * Routine: Fix_Refined_Block_Boundaries                              *
 *                                                                    *
 * Adjusts the locations of the boundary nodes of a solution block so *
 * that the new node locations match with cell volumes of adjacent    *
 * solution blocks that have lower levels of mesh refinement (i.e.,   *
 * are coarser solution blocks).                                      *
 *                                                                    *
 **********************************************************************/
void Fix_Refined_Block_Boundaries(NavierStokes2D_Quad_Block &SolnBlk,
                                  const int Fix_North_Boundary,
                                  const int Fix_South_Boundary,
                                  const int Fix_East_Boundary,
                                  const int Fix_West_Boundary) {

  double ds_ratio, dl, dr;

  // Adjust the node locations at the north boundary.
  if (Fix_North_Boundary) {
    for (int i = SolnBlk.Grid.INl+1; i <= SolnBlk.Grid.INu-1; i+=2) {
      dl = abs(SolnBlk.Grid.Node[i  ][SolnBlk.Grid.JNu].X - 
	       SolnBlk.Grid.Node[i-1][SolnBlk.Grid.JNu].X);
      dr = abs(SolnBlk.Grid.Node[i+1][SolnBlk.Grid.JNu].X - 
	       SolnBlk.Grid.Node[i  ][SolnBlk.Grid.JNu].X);
      ds_ratio = dl/(dl+dr);
      SolnBlk.Grid.Node[i][SolnBlk.Grid.JNu].X = 
	SolnBlk.Grid.Node[i-1][SolnBlk.Grid.JNu].X +
	ds_ratio*(SolnBlk.Grid.Node[i+1][SolnBlk.Grid.JNu].X-
		  SolnBlk.Grid.Node[i-1][SolnBlk.Grid.JNu].X);
    }
    for (int i = SolnBlk.ICl; i <= SolnBlk.ICu; i++) {
      SolnBlk.U[i][SolnBlk.JCu] = (SolnBlk.Grid.Cell[i][SolnBlk.JCu].A/
				   SolnBlk.Grid.area(i,SolnBlk.JCu))*SolnBlk.U[i][SolnBlk.JCu];
      SolnBlk.W[i][SolnBlk.JCu] = W(SolnBlk.U[i][SolnBlk.JCu]);
    }
  }

  // Adjust the node locations at the south boundary.
  if (Fix_South_Boundary) {
    for (int i = SolnBlk.Grid.INl+1; i <= SolnBlk.Grid.INu-1; i+=2) {
      dl = abs(SolnBlk.Grid.Node[i  ][SolnBlk.Grid.JNl].X - 
	       SolnBlk.Grid.Node[i-1][SolnBlk.Grid.JNl].X);
      dr = abs(SolnBlk.Grid.Node[i+1][SolnBlk.Grid.JNl].X - 
	       SolnBlk.Grid.Node[i  ][SolnBlk.Grid.JNl].X);
      ds_ratio = dl/(dl+dr);
      SolnBlk.Grid.Node[i][SolnBlk.Grid.JNl].X = 
	SolnBlk.Grid.Node[i-1][SolnBlk.Grid.JNl].X +
	ds_ratio*(SolnBlk.Grid.Node[i+1][SolnBlk.Grid.JNl].X-
		  SolnBlk.Grid.Node[i-1][SolnBlk.Grid.JNl].X);
    }
    for (int i = SolnBlk.ICl; i <= SolnBlk.ICu; i++) {
      SolnBlk.U[i][SolnBlk.JCl] = (SolnBlk.Grid.Cell[i][SolnBlk.JCl].A/
				   SolnBlk.Grid.area(i,SolnBlk.JCl))*SolnBlk.U[i][SolnBlk.JCl];
      SolnBlk.W[i][SolnBlk.JCl] = W(SolnBlk.U[i][SolnBlk.JCl]);
    }
  }

  // Adjust the node locations at the east boundary.
  if (Fix_East_Boundary) {
    for (int j = SolnBlk.Grid.JNl+1; j <= SolnBlk.Grid.JNu-1; j+=2) {
      dl = abs(SolnBlk.Grid.Node[SolnBlk.Grid.INu][j  ].X - 
	       SolnBlk.Grid.Node[SolnBlk.Grid.INu][j-1].X);
      dr = abs(SolnBlk.Grid.Node[SolnBlk.Grid.INu][j+1].X - 
	       SolnBlk.Grid.Node[SolnBlk.Grid.INu][j  ].X);
      ds_ratio = dl/(dl+dr);
      SolnBlk.Grid.Node[SolnBlk.Grid.INu][j].X = 
	SolnBlk.Grid.Node[SolnBlk.Grid.INu][j-1].X +
	ds_ratio*(SolnBlk.Grid.Node[SolnBlk.Grid.INu][j+1].X-
		  SolnBlk.Grid.Node[SolnBlk.Grid.INu][j-1].X);
    }
    for (int j = SolnBlk.JCl; j <= SolnBlk.JCu; j++) {
      SolnBlk.U[SolnBlk.ICu][j] = (SolnBlk.Grid.Cell[SolnBlk.ICu][j].A/
				   SolnBlk.Grid.area(SolnBlk.ICu,j))*SolnBlk.U[SolnBlk.ICu][j];
      SolnBlk.W[SolnBlk.ICu][j] = W(SolnBlk.U[SolnBlk.ICu][j]);
    }
  }

  // Adjust the node locations at the west boundary.
  if (Fix_West_Boundary) {
    for (int j = SolnBlk.Grid.JNl+1; j <= SolnBlk.Grid.JNu-1; j+=2) {
      dl = abs(SolnBlk.Grid.Node[SolnBlk.Grid.INl][j  ].X - 
	       SolnBlk.Grid.Node[SolnBlk.Grid.INl][j-1].X);
      dr = abs(SolnBlk.Grid.Node[SolnBlk.Grid.INl][j+1].X - 
	       SolnBlk.Grid.Node[SolnBlk.Grid.INl][j  ].X);
      ds_ratio = dl/(dl+dr);
      SolnBlk.Grid.Node[SolnBlk.Grid.INl][j].X = 
	SolnBlk.Grid.Node[SolnBlk.Grid.INl][j-1].X +
	ds_ratio*(SolnBlk.Grid.Node[SolnBlk.Grid.INl][j+1].X-
		  SolnBlk.Grid.Node[SolnBlk.Grid.INl][j-1].X);
    }
    for (int j = SolnBlk.JCl; j <= SolnBlk.JCu; j++) {
      SolnBlk.U[SolnBlk.ICl][j] = (SolnBlk.Grid.Cell[SolnBlk.ICl][j].A/
				   SolnBlk.Grid.area(SolnBlk.ICl,j))*SolnBlk.U[SolnBlk.ICl][j];
      SolnBlk.W[SolnBlk.ICl][j] = W(SolnBlk.U[SolnBlk.ICl][j]);
    }
  }

  // Reset the boundary condition types at the block boundaries.
  Set_BCs(SolnBlk.Grid);

  // Recompute the exterior nodes for the block quadrilateral mesh.
  Update_Exterior_Nodes(SolnBlk.Grid);

  // Recompute the cells for the block quadrilateral mesh.
  Update_Cells(SolnBlk.Grid);

}

/**********************************************************************
 * Routine: Unfix_Refined_Block_Boundaries                            *
 *                                                                    *
 * Returns the adjusted the locations of the boundary nodes of a      *
 * solution block to their original unmodified positions.             *
 *                                                                    *
 **********************************************************************/
void Unfix_Refined_Block_Boundaries(NavierStokes2D_Quad_Block &SolnBlk) {

  double sp_l, sp_r, sp_m, ds_ratio, dl, dr;
 
  // Return the nodes at the north boundary to their original positions.
  if (SolnBlk.Grid.BndNorthSpline.np != 0) {
    for (int i = SolnBlk.Grid.INl+1; i < SolnBlk.Grid.INu; i += 2) {
      sp_l = getS(SolnBlk.Grid.Node[i-1][SolnBlk.Grid.JNu].X,
		  SolnBlk.Grid.BndNorthSpline);
      sp_r = getS(SolnBlk.Grid.Node[i+1][SolnBlk.Grid.JNu].X,
		  SolnBlk.Grid.BndNorthSpline);
      dl = abs(SolnBlk.Grid.Node[i  ][SolnBlk.Grid.JNu].X -
	       SolnBlk.Grid.Node[i-1][SolnBlk.Grid.JNu].X);
      dr = abs(SolnBlk.Grid.Node[i+1][SolnBlk.Grid.JNu].X -
	       SolnBlk.Grid.Node[i  ][SolnBlk.Grid.JNu].X);
      ds_ratio = dl/(dl+dr);
      sp_m = sp_l + ds_ratio*(sp_r-sp_l);
      SolnBlk.Grid.Node[i][SolnBlk.Grid.JNu].X = Spline(sp_m,SolnBlk.Grid.BndNorthSpline);
    }
    for (int i = SolnBlk.ICl; i <= SolnBlk.ICu; i++) {
      SolnBlk.U[i][SolnBlk.JCu] = (SolnBlk.Grid.Cell[i][SolnBlk.JCu].A/
				   SolnBlk.Grid.area(i,SolnBlk.JCu))*SolnBlk.U[i][SolnBlk.JCu];
      SolnBlk.W[i][SolnBlk.JCu] = W(SolnBlk.U[i][SolnBlk.JCu]);
    }
  }

  // Return the nodes at the south boundary to their original positions.
  if (SolnBlk.Grid.BndSouthSpline.np != 0) {
    for (int i = SolnBlk.Grid.INl+1; i < SolnBlk.Grid.INu; i += 2) {
      sp_l = getS(SolnBlk.Grid.Node[i-1][SolnBlk.Grid.JNl].X,
		  SolnBlk.Grid.BndSouthSpline);
      sp_r = getS(SolnBlk.Grid.Node[i+1][SolnBlk.Grid.JNl].X,
		  SolnBlk.Grid.BndSouthSpline);
      dl = abs(SolnBlk.Grid.Node[i  ][SolnBlk.Grid.JNl].X -
	       SolnBlk.Grid.Node[i-1][SolnBlk.Grid.JNl].X);
      dr = abs(SolnBlk.Grid.Node[i+1][SolnBlk.Grid.JNl].X -
	       SolnBlk.Grid.Node[i  ][SolnBlk.Grid.JNl].X);
      ds_ratio = dl/(dl+dr);
      sp_m = sp_l + ds_ratio*(sp_r-sp_l);
      SolnBlk.Grid.Node[i][SolnBlk.Grid.JNl].X = Spline(sp_m,SolnBlk.Grid.BndSouthSpline);
    }
    for (int i = SolnBlk.ICl; i <= SolnBlk.ICu; i++) {
      SolnBlk.U[i][SolnBlk.JCl] = (SolnBlk.Grid.Cell[i][SolnBlk.JCl].A/
				   SolnBlk.Grid.area(i,SolnBlk.JCl))*SolnBlk.U[i][SolnBlk.JCl];
      SolnBlk.W[i][SolnBlk.JCl] = W(SolnBlk.U[i][SolnBlk.JCl]);
    }
  }

  // Return the nodes at the east boundary to their original positions.
  if (SolnBlk.Grid.BndEastSpline.np != 0) {
    for (int j = SolnBlk.Grid.JNl+1; j < SolnBlk.Grid.JNu; j += 2) {
      sp_l = getS(SolnBlk.Grid.Node[SolnBlk.Grid.INu][j-1].X,
		  SolnBlk.Grid.BndEastSpline);
      sp_r = getS(SolnBlk.Grid.Node[SolnBlk.Grid.INu][j+1].X,
		  SolnBlk.Grid.BndEastSpline);
      dl = abs(SolnBlk.Grid.Node[SolnBlk.Grid.INu][j  ].X -
	       SolnBlk.Grid.Node[SolnBlk.Grid.INu][j-1].X);
      dr = abs(SolnBlk.Grid.Node[SolnBlk.Grid.INu][j+1].X -
	       SolnBlk.Grid.Node[SolnBlk.Grid.INu][j  ].X);
      ds_ratio = dl/(dl+dr);
      sp_m = sp_l + ds_ratio*(sp_r-sp_l);
      SolnBlk.Grid.Node[SolnBlk.Grid.INu][j].X = Spline(sp_m,SolnBlk.Grid.BndEastSpline);
    }
    for (int j = SolnBlk.JCl; j <= SolnBlk.JCu; j++) {
      SolnBlk.U[SolnBlk.ICu][j] = (SolnBlk.Grid.Cell[SolnBlk.ICu][j].A/
				   SolnBlk.Grid.area(SolnBlk.ICu,j))*SolnBlk.U[SolnBlk.ICu][j];
      SolnBlk.W[SolnBlk.ICu][j] = W(SolnBlk.U[SolnBlk.ICu][j]);
    }
  }

  // Return the nodes at the west boundary to their original positions.
  if (SolnBlk.Grid.BndWestSpline.np != 0) {
    for (int j = SolnBlk.Grid.JNl+1; j < SolnBlk.Grid.JNu; j += 2) {
      sp_l = getS(SolnBlk.Grid.Node[SolnBlk.Grid.INl][j-1].X,
		  SolnBlk.Grid.BndWestSpline);
      sp_r = getS(SolnBlk.Grid.Node[SolnBlk.Grid.INl][j+1].X,
		  SolnBlk.Grid.BndWestSpline);
      dl = abs(SolnBlk.Grid.Node[SolnBlk.Grid.INl][j  ].X -
	       SolnBlk.Grid.Node[SolnBlk.Grid.INl][j-1].X);
      dr = abs(SolnBlk.Grid.Node[SolnBlk.Grid.INl][j+1].X -
	       SolnBlk.Grid.Node[SolnBlk.Grid.INl][j  ].X);
      ds_ratio = dl/(dl+dr);
      sp_m = sp_l + ds_ratio*(sp_r-sp_l);
      SolnBlk.Grid.Node[SolnBlk.Grid.INl][j].X = Spline(sp_m,SolnBlk.Grid.BndWestSpline);
    }
    for (int j = SolnBlk.JCl; j <= SolnBlk.JCu; j++) {
      SolnBlk.U[SolnBlk.ICl][j] = (SolnBlk.Grid.Cell[SolnBlk.ICl][j].A/
				   SolnBlk.Grid.area(SolnBlk.ICl,j))*SolnBlk.U[SolnBlk.ICl][j];
      SolnBlk.W[SolnBlk.ICl][j] = W(SolnBlk.U[SolnBlk.ICl][j]);
    }
  }

  // Reset the boundary condition types at the block boundaries.
  Set_BCs(SolnBlk.Grid);

  // Recompute the exterior nodes for the block quadrilateral mesh.
  Update_Exterior_Nodes(SolnBlk.Grid);

  // Recompute the cells for the block quadrilateral mesh.
  Update_Cells(SolnBlk.Grid);

}

/**********************************************************************
 * Routine: Apply_Boundary_Flux_Corrections                           *
 *                                                                    *
 * Apply flux corrections at boundaries of the solution block to      *
 * ensure that the scheme is conservative at boundaries with mesh     *
 * resolution changes.                                                *
 *                                                                    *
 **********************************************************************/
void Apply_Boundary_Flux_Corrections(NavierStokes2D_Quad_Block &SolnBlk,
                                     const int Number_Neighbours_North_Boundary,
                                     const int Number_Neighbours_South_Boundary,
                                     const int Number_Neighbours_East_Boundary,
                                     const int Number_Neighbours_West_Boundary) {

  // Correct the fluxes at the north boundary as required.
  if (Number_Neighbours_North_Boundary == 2)
    for (int i = SolnBlk.ICl; i <= SolnBlk.ICu; i++)
      SolnBlk.dUdt[i][SolnBlk.JCu][0] -=
	SolnBlk.FluxN[i]/SolnBlk.Grid.Cell[i][SolnBlk.JCu].A;

  // Correct the fluxes at the south boundary as required.
  if (Number_Neighbours_South_Boundary == 2)
    for (int i = SolnBlk.ICl; i <= SolnBlk.ICu; i++)
      SolnBlk.dUdt[i][SolnBlk.JCl][0] -=
	SolnBlk.FluxS[i]/SolnBlk.Grid.Cell[i][SolnBlk.JCl].A;

  // Correct the fluxes at the east boundary as required.
  if (Number_Neighbours_East_Boundary == 2)
    for (int j = SolnBlk.JCl; j <= SolnBlk.JCu; j++)
      SolnBlk.dUdt[SolnBlk.ICu][j][0] -=
	SolnBlk.FluxE[j]/SolnBlk.Grid.Cell[SolnBlk.ICu][j].A;

  // Correct the fluxes at the west boundary as required.
  if (Number_Neighbours_West_Boundary == 2)
    for (int j = SolnBlk.JCl; j <= SolnBlk.JCu; j++)
      SolnBlk.dUdt[SolnBlk.ICl][j][0] -=
	SolnBlk.FluxW[j]/SolnBlk.Grid.Cell[SolnBlk.ICl][j].A;

}

/**********************************************************************
 * Routine: Apply_Boundary_Flux_Corrections_Multistage_Explicit       *
 *                                                                    *
 * Apply flux corrections at boundaries of the solution block to      *
 * ensure that the scheme is conservative at boundaries with mesh     *
 * resolution changes.                                                *
 *                                                                    *
 **********************************************************************/
void Apply_Boundary_Flux_Corrections_Multistage_Explicit(NavierStokes2D_Quad_Block &SolnBlk,
                                                         const int i_stage,
                                                         NavierStokes2D_Input_Parameters &IP,
                                                         const int Number_Neighbours_North_Boundary,
                                                         const int Number_Neighbours_South_Boundary,
                                                         const int Number_Neighbours_East_Boundary,
                                                         const int Number_Neighbours_West_Boundary) {

  int k_residual;
  double omega;

  // Evaluate the time step fraction and residual storage location for the stage.
  switch(IP.i_Time_Integration) {
  case TIME_STEPPING_EXPLICIT_EULER :
    omega = Runge_Kutta(i_stage,IP.N_Stage);
    k_residual = 0;
    break;
  case TIME_STEPPING_EXPLICIT_PREDICTOR_CORRECTOR :
    omega = Runge_Kutta(i_stage,IP.N_Stage);
    k_residual = 0;
    break;
  case TIME_STEPPING_EXPLICIT_RUNGE_KUTTA :
    omega = Runge_Kutta(i_stage,IP.N_Stage);
    k_residual = 0;
    if (IP.N_Stage == 4) {
      if (i_stage == 4) {
	k_residual = 0;
      } else {
	k_residual = i_stage - 1;
      }
    }
    break;
  case TIME_STEPPING_MULTISTAGE_OPTIMAL_SMOOTHING :
    omega = MultiStage_Optimally_Smoothing(i_stage,IP.N_Stage,IP.i_Limiter);
    k_residual = 0;
    break;
  default:
    omega = Runge_Kutta(i_stage,IP.N_Stage);
    k_residual = 0;
    break;
  };

  // Correct the fluxes at the north boundary as required.
  if (Number_Neighbours_North_Boundary == 2)
    for (int i = SolnBlk.ICl; i <= SolnBlk.ICu; i++)
      SolnBlk.dUdt[i][SolnBlk.JCu][k_residual] -= 
	(IP.CFL_Number*SolnBlk.dt[i][SolnBlk.JCu])*SolnBlk.FluxN[i]/
	SolnBlk.Grid.Cell[i][SolnBlk.JCu].A;
  
  // Correct the fluxes at the south boundary as required.
  if (Number_Neighbours_South_Boundary == 2)
    for (int i = SolnBlk.ICl; i <= SolnBlk.ICu; i++)
      SolnBlk.dUdt[i][SolnBlk.JCl][k_residual] -= 
	(IP.CFL_Number*SolnBlk.dt[i][SolnBlk.JCl])*SolnBlk.FluxS[i]/
	SolnBlk.Grid.Cell[i][SolnBlk.JCl].A;

  // Correct the fluxes at the east boundary as required.
  if (Number_Neighbours_East_Boundary == 2)
    for (int j = SolnBlk.JCl; j <= SolnBlk.JCu; j++)
      SolnBlk.dUdt[SolnBlk.ICu][j][k_residual] -= 
	(IP.CFL_Number*SolnBlk.dt[SolnBlk.ICu][j])*SolnBlk.FluxE[j]/
	SolnBlk.Grid.Cell[SolnBlk.ICu][j].A;

  // Correct the fluxes at the west boundary as required.
  if (Number_Neighbours_West_Boundary == 2)
    for (int j= SolnBlk.JCl; j <= SolnBlk.JCu; j++)
      SolnBlk.dUdt[SolnBlk.ICl][j][k_residual] -= 
	(IP.CFL_Number*SolnBlk.dt[SolnBlk.ICl][j])*SolnBlk.FluxW[j]/
	SolnBlk.Grid.Cell[SolnBlk.ICl][j].A;

}

/**********************************************************************
 * Routine: dUdt_Residual_Evaluation                                  *
 *                                                                    *
 * This routine evaluates the residual for the specified solution     *
 * block using a 2nd-order limited upwind finite-volume spatial       *
 * discretization scheme with either the Godunov, Roe, Rusanov, HLLE, *
 * HLLL, or HLLC flux functions for the gas-phase and with either the *
 * Saurel or Equilibrium flux functions for the particle-phase.  The  *
 * residual is stored in dUdt[][][0].                                 *
 *                                                                    *
 **********************************************************************/
int dUdt_Residual_Evaluation(NavierStokes2D_Quad_Block &SolnBlk,
			     NavierStokes2D_Input_Parameters &IP) {

  int error_flag;
  Vector2D dX;
  NavierStokes2D_pState Wl, Wr;
  NavierStokes2D_cState Flux;

  NavierStokes2D_pState Wu, Wd, dWdxl, dWdyl, dWdxr, dWdyr;
  NavierStokes2D_pState dWdx, dWdy;
  Vector2D Xl, Xr, Xu, Xd;
  int viscous_bc_flag;

  int flux_function = IP.i_Flux_Function;

  NavierStokes2D_cState NavierStokes2D_U_VACUUM; NavierStokes2D_U_VACUUM.Vacuum();
  NavierStokes2D_pState NavierStokes2D_W_STDATM; NavierStokes2D_W_STDATM.Standard_Atmosphere();

  // Perform the linear reconstruction within each cell of the
  // computational grid for this stage.
  switch(IP.i_Reconstruction) {
  case RECONSTRUCTION_GREEN_GAUSS :
    Linear_Reconstruction_GreenGauss(SolnBlk,IP.i_Limiter);
    break;
  case RECONSTRUCTION_LINEAR_LEAST_SQUARES :
    Linear_Reconstruction_LeastSquares(SolnBlk,IP.i_Limiter);
    break;
  default:
    Linear_Reconstruction_LeastSquares(SolnBlk,IP.i_Limiter);
    break;
  };

  // Evaluate the time rate of change of the solution (i.e., the
  // solution residuals) using a second-order limited upwind scheme
  // with a variety of flux functions.

  // Add i-direction (zeta-direction) fluxes.
  for (int j = SolnBlk.JCl-1; j <= SolnBlk.JCu+1; j++) {

    SolnBlk.dUdt[SolnBlk.ICl-1][j][0] = NavierStokes2D_U_VACUUM;

    for (int i = SolnBlk.ICl-1; i <= SolnBlk.ICu; i++) {

      SolnBlk.dUdt[i+1][j][0] = NavierStokes2D_U_VACUUM;

      if (j >= SolnBlk.JCl && j <= SolnBlk.JCu) {

	if (i == SolnBlk.ICl-1 && 
	    (SolnBlk.Grid.BCtypeW[j] == BC_REFLECTION ||
	     SolnBlk.Grid.BCtypeW[j] == BC_BURNING_SURFACE ||
	     SolnBlk.Grid.BCtypeW[j] == BC_MASS_INJECTION ||
	     SolnBlk.Grid.BCtypeW[j] == BC_WALL_VISCOUS_HEATFLUX ||
	     SolnBlk.Grid.BCtypeW[j] == BC_WALL_VISCOUS_ISOTHERMAL ||
	     SolnBlk.Grid.BCtypeW[j] == BC_MOVING_WALL_HEATFLUX ||
	     SolnBlk.Grid.BCtypeW[j] == BC_MOVING_WALL_ISOTHERMAL ||
	     SolnBlk.Grid.BCtypeW[j] == BC_RINGLEB_FLOW ||
	     SolnBlk.Grid.BCtypeW[j] == BC_CHARACTERISTIC)) {

	  dX = SolnBlk.Grid.xfaceW(i+1,j) - SolnBlk.Grid.Cell[i+1][j].Xc;
	  Wr = SolnBlk.W[i+1][j] + (SolnBlk.phi[i+1][j]^SolnBlk.dWdx[i+1][j])*dX.x +
	                           (SolnBlk.phi[i+1][j]^SolnBlk.dWdy[i+1][j])*dX.y;

	  // WEST face of cell (i+1,j) is a normal boundary.
	  if (SolnBlk.Grid.BCtypeW[j] == BC_REFLECTION) {
	    // WEST face of cell (i+1,j) is a REFLECTION boundary.
	    Wl = Reflect(Wr,SolnBlk.Grid.nfaceW(i+1,j));
	  } else if (SolnBlk.Grid.BCtypeW[j] == BC_BURNING_SURFACE) {
	    // WEST face of cell (i+1,j) is a BURNING_SURFACE boundary.
	    Wl = BurningSurface(Wr,SolnBlk.Grid.nfaceW(i+1,j));
	  } else if (SolnBlk.Grid.BCtypeW[j] == BC_MASS_INJECTION) {
	    // WEST face of cell (i+1,j) is a MASS_INJECTION boundary.
	    //Wl = MassInjection(Wr,SolnBlk.Grid.nfaceW(i+1,j));
	    Wl = MassInjection2(Wr,SolnBlk.Grid.xfaceW(i+1,j),SolnBlk.Grid.nfaceW(i+1,j),SolnBlk.Twall);
	    flux_function = IP.i_Flux_Function;
	    IP.i_Flux_Function = FLUX_FUNCTION_GODUNOV_WRS;
	  } else if (SolnBlk.Grid.BCtypeW[j] == BC_WALL_VISCOUS_HEATFLUX) {
	    // WEST face of cell (i+1,j) is a WALL_VISCOUS_HEATFLUX boundary.
	    Wl = WallViscousHeatFlux(Wr,SolnBlk.Grid.nfaceW(i+1,j));
	  } else if (SolnBlk.Grid.BCtypeW[j] == BC_WALL_VISCOUS_ISOTHERMAL) {
	    // WEST face of cell (i+1,j) is a WALL_VISCOUS_ISOTHERMAL boundary.
	    Wl = WallViscousIsothermal(Wr,SolnBlk.Grid.nfaceW(i+1,j),SolnBlk.Twall);
	  } else if (SolnBlk.Grid.BCtypeW[j] == BC_MOVING_WALL) {
	    // WEST face of cell (i+1,j) is a MOVINGWALL boundary.
	    Wl = MovingWallHeatFlux(Wr,SolnBlk.Grid.nfaceW(i+1,j),SolnBlk.Vwall.x);
	  } else if (SolnBlk.Grid.BCtypeW[j] == BC_MOVING_WALL_ISOTHERMAL) {
	    // WEST face of cell (i+1,j) is a MOVINGWALL_ISOTHERMAL boundary.
	    Wl = MovingWallIsothermal(Wr,SolnBlk.Grid.nfaceW(i+1,j),SolnBlk.Vwall.x,SolnBlk.Twall);
	  } else if (SolnBlk.Grid.BCtypeW[j] == BC_RINGLEB_FLOW) {
	    // WEST face of cell (i+1,j) is a RINGLEB_FLOW boundary.
	    Wl = RinglebFlow(Wl,SolnBlk.Grid.xfaceW(i+1,j));
	  } else {
	    // WEST face of cell (i+1,j) is a CHARACTERISTIC boundary.
	    Wl = BC_Characteristic_Pressure(Wr,SolnBlk.WoW[j],SolnBlk.Grid.nfaceW(i+1,j));
	  }

	} else if (i == SolnBlk.ICu &&
		   (SolnBlk.Grid.BCtypeE[j] == BC_REFLECTION ||
		    SolnBlk.Grid.BCtypeE[j] == BC_BURNING_SURFACE ||
		    SolnBlk.Grid.BCtypeE[j] == BC_MASS_INJECTION ||
		    SolnBlk.Grid.BCtypeE[j] == BC_WALL_VISCOUS_HEATFLUX ||
		    SolnBlk.Grid.BCtypeE[j] == BC_WALL_VISCOUS_ISOTHERMAL ||
		    SolnBlk.Grid.BCtypeE[j] == BC_MOVING_WALL_HEATFLUX ||
		    SolnBlk.Grid.BCtypeE[j] == BC_MOVING_WALL_ISOTHERMAL ||
		    SolnBlk.Grid.BCtypeE[j] == BC_RINGLEB_FLOW ||
		    SolnBlk.Grid.BCtypeE[j] == BC_CHARACTERISTIC)) {

	  dX = SolnBlk.Grid.xfaceE(i,j) - SolnBlk.Grid.Cell[i][j].Xc;
	  Wl = SolnBlk.W[i][j] + (SolnBlk.phi[i][j]^SolnBlk.dWdx[i][j])*dX.x +
	                         (SolnBlk.phi[i][j]^SolnBlk.dWdy[i][j])*dX.y;

	  // EAST face of cell (i,j) is a normal boundary.
	  if (SolnBlk.Grid.BCtypeE[j] == BC_REFLECTION) {
	    // EAST face of cell (i,j) is a REFLECTION boundary.
	    Wr = Reflect(Wl,SolnBlk.Grid.nfaceE(i,j));
	  } else if (SolnBlk.Grid.BCtypeE[j] == BC_BURNING_SURFACE) {
	    // EAST face of cell (i,j) is a BURNING_SURFACE boundary.
	    Wr = BurningSurface(Wl,SolnBlk.Grid.nfaceE(i,j));
	  } else if (SolnBlk.Grid.BCtypeE[j] == BC_MASS_INJECTION) {
	    // EAST face of cell (i,j) is a MASS_INJECTION boundary.
	    //Wr = MassInjection(Wl,SolnBlk.Grid.nfaceE(i,j));
	    Wr = MassInjection2(Wl,SolnBlk.Grid.xfaceE(i,j),SolnBlk.Grid.nfaceE(i,j),SolnBlk.Twall);
	    flux_function = IP.i_Flux_Function;
	    IP.i_Flux_Function = FLUX_FUNCTION_GODUNOV_WRS;
	  } else if (SolnBlk.Grid.BCtypeE[j] == BC_WALL_VISCOUS_HEATFLUX) {
	    // EAST face of cell (i,j) is a WALL_VISCOUS_HEATFLUX boundary.
	    Wr = WallViscousHeatFlux(Wl,SolnBlk.Grid.nfaceE(i,j));
	  } else if (SolnBlk.Grid.BCtypeE[j] == BC_WALL_VISCOUS_ISOTHERMAL) {
	    // EAST face of cell (i,j) is a WALL_VISCOUS_ISOTHERMAL boundary.
	    Wr = WallViscousIsothermal(Wl,SolnBlk.Grid.nfaceE(i,j),SolnBlk.Twall);
	  } else if (SolnBlk.Grid.BCtypeE[j] == BC_MOVING_WALL_HEATFLUX) {
	    // EAST face of cell (i,j) is a MOVINGWALL_HEATFLUX boundary.
	    Wr = MovingWallHeatFlux(Wl,SolnBlk.Grid.nfaceE(i,j),SolnBlk.Vwall.x);
	  } else if (SolnBlk.Grid.BCtypeE[j] == BC_MOVING_WALL_ISOTHERMAL) {
	    // EAST face of cell (i,j) is a MOVINGWALL_ISOTHERMAL boundary.
	    Wr = MovingWallIsothermal(Wl,SolnBlk.Grid.nfaceE(i,j),SolnBlk.Vwall.x,SolnBlk.Twall);
	  } else if (SolnBlk.Grid.BCtypeE[j] == BC_RINGLEB_FLOW) {
	    // EAST face of cell (i,j) is a RINGLEB_FLOW boundary.
	    Wr = RinglebFlow(Wr,SolnBlk.Grid.xfaceE(i,j));
	  } else {
	    // EAST face of cell (i,j) is a CHARACTERISTIC boundary.
	    Wr = BC_Characteristic_Pressure(Wl,SolnBlk.WoE[j],SolnBlk.Grid.nfaceE(i,j));
	  }

	} else {

	  // EAST face is either a normal cell or possibly a FIXED, 
	  // NONE or EXTRAPOLATION boundary.
	  dX = SolnBlk.Grid.xfaceE(i  ,j) - SolnBlk.Grid.Cell[i  ][j].Xc;
	  Wl = SolnBlk.W[i  ][j] + (SolnBlk.phi[i  ][j]^SolnBlk.dWdx[i  ][j])*dX.x +
	                           (SolnBlk.phi[i  ][j]^SolnBlk.dWdy[i  ][j])*dX.y;
	  dX = SolnBlk.Grid.xfaceW(i+1,j) - SolnBlk.Grid.Cell[i+1][j].Xc;
	  Wr = SolnBlk.W[i+1][j] + (SolnBlk.phi[i+1][j]^SolnBlk.dWdx[i+1][j])*dX.x +
	                           (SolnBlk.phi[i+1][j]^SolnBlk.dWdy[i+1][j])*dX.y;

	}

	// Determine EAST face INVISCID flux.
	switch(IP.i_Flux_Function) {
	case FLUX_FUNCTION_GODUNOV :
	  Flux = FluxGodunov_n(Wl,Wr,SolnBlk.Grid.nfaceE(i,j));
	  break;
	case FLUX_FUNCTION_ROE :
	  Flux = FluxRoe_n(Wl,Wr,SolnBlk.Grid.nfaceE(i,j));
	  break;
	case FLUX_FUNCTION_RUSANOV :
	  Flux = FluxRusanov_n(Wl,Wr,SolnBlk.Grid.nfaceE(i,j));
	  break;
	case FLUX_FUNCTION_HLLE :
	  Flux = FluxHLLE_n(Wl,Wr,SolnBlk.Grid.nfaceE(i,j));
	  break;
	case FLUX_FUNCTION_HLLL :
	  Flux = FluxHLLL_n(Wl,Wr,SolnBlk.Grid.nfaceE(i,j));
	  break;
	case FLUX_FUNCTION_HLLC :
	  Flux = FluxHLLC_n(Wl,Wr,SolnBlk.Grid.nfaceE(i,j));
	  break;
	case FLUX_FUNCTION_VANLEER :
	  Flux = FluxVanLeer_n(Wl,Wr,SolnBlk.Grid.nfaceE(i,j));
	  break;
	case FLUX_FUNCTION_AUSM :
	  Flux = FluxAUSM_n(Wl,Wr,SolnBlk.Grid.nfaceE(i,j));
	  break;
	case FLUX_FUNCTION_AUSMplus :
	  Flux = FluxAUSMplus_n(Wl,Wr,SolnBlk.Grid.nfaceE(i,j));
	  break;
	case FLUX_FUNCTION_GODUNOV_WRS :
	  Flux = FluxGodunov_Wrs_n(Wl,Wr,SolnBlk.Grid.nfaceE(i,j));
	  break;
	default:
	  Flux = FluxRoe_n(Wl,Wr,SolnBlk.Grid.nfaceE(i,j));
	  break;
	};

	if (IP.i_Flux_Function == FLUX_FUNCTION_GODUNOV_WRS) {
	  IP.i_Flux_Function = flux_function;
	}

	// Compute the cell centred stress tensor and heat flux vector if required.
 	if (SolnBlk.Flow_Type == FLOWTYPE_TURBULENT_RANS_K_OMEGA ||
	    (SolnBlk.Flow_Type && SolnBlk.Axisymmetric)) {
	  SolnBlk.W[i][j].ComputeViscousTerms(SolnBlk.dWdx[i][j],
					      SolnBlk.dWdy[i][j],
					      SolnBlk.Grid.Cell[i][j].Xc,
					      SolnBlk.Axisymmetric,
					      OFF,
					      SolnBlk.Wall[i][j].ywall,
					      SolnBlk.Wall[i][j].yplus);
	  SolnBlk.U[i][j].tau = SolnBlk.W[i][j].tau;
	  SolnBlk.U[i][j].q = SolnBlk.W[i][j].q;
	}

	// Evaluate the cell interface i-direction VISCOUS flux if necessary.
	if (SolnBlk.Flow_Type) {
	  // Determine the EAST face VISCOUS flux.
	  if (i == SolnBlk.ICl-1 && 
	      (SolnBlk.Grid.BCtypeW[j] == BC_WALL_VISCOUS_HEATFLUX ||
	       SolnBlk.Grid.BCtypeW[j] == BC_WALL_VISCOUS_ISOTHERMAL ||
	       SolnBlk.Grid.BCtypeW[j] == BC_MOVING_WALL_HEATFLUX ||
	       SolnBlk.Grid.BCtypeW[j] == BC_MOVING_WALL_ISOTHERMAL ||
	       SolnBlk.Grid.BCtypeW[j] == BC_BURNING_SURFACE ||
	       SolnBlk.Grid.BCtypeW[j] == BC_MASS_INJECTION)) {
	    // WEST face of cell (i+1,j) is a normal boundary.
	    Xr = SolnBlk.Grid.Cell[i+1][j].Xc; Wr = SolnBlk.W[i+1][j];
	    if (SolnBlk.Grid.BCtypeW[j] == BC_WALL_VISCOUS_HEATFLUX) {
	      // WEST face of cell (i+1,j) is a WALL_VISCOUS_HEATFLUX boundary.
	      viscous_bc_flag = DIAMONDPATH_RIGHT_TRIANGLE_HEATFLUX;
	      //Wu.rho = Wr.rho; Wu.v.x = ZERO; Wu.v.y = ZERO; Wu.p = Wr.p; Wu.k = Wr.k; Wu.omega = Wr.omega;Wu.ke = Wr.ke; Wu.ee = Wr.ee;
	      Wu = HALF*(Wr+SolnBlk.W[i+1][j+1]); Wu.v = Vector2D_ZERO;
	      Wd = HALF*(Wr+SolnBlk.W[i+1][j-1]); Wd.v = Vector2D_ZERO;
	    } else if (SolnBlk.Grid.BCtypeW[j] == BC_WALL_VISCOUS_ISOTHERMAL) {
	      // WEST face of cell (i+1,j) is a WALL_VISCOUS_ISOTHERMAL boundary.
	      viscous_bc_flag = DIAMONDPATH_RIGHT_TRIANGLE_ISOTHERMAL;
	      //Wu.rho = Wr.p/(Wr.R*SolnBlk.Twall); Wu.v.x = ZERO; Wu.v.y = ZERO; Wu.p = Wr.p; Wu.k = Wr.k; Wu.omega = Wr.omega;Wu.ke = Wr.ke; Wu.ee = Wr.ee;
	      Wu = HALF*(Wr+SolnBlk.W[i+1][j+1]); Wu.v = Vector2D_ZERO; Wu.rho = Wu.p/(Wr.R*SolnBlk.Twall);
	      Wd = HALF*(Wr+SolnBlk.W[i+1][j-1]); Wd.v = Vector2D_ZERO; Wd.rho = Wd.p/(Wr.R*SolnBlk.Twall);
	    } else if (SolnBlk.Grid.BCtypeW[j] == BC_MOVING_WALL) {
	      // WEST face of cell (i+1,j) is a MOVINGWALL_HEATFLUX boundary.
	      viscous_bc_flag = DIAMONDPATH_RIGHT_TRIANGLE_HEATFLUX;
	      //Wu.rho = Wr.rho; Wu.v.x = SolnBlk.Vwall.x; Wu.v.y = SolnBlk.Vwall.y; Wu.p = Wr.p; Wu.k = Wr.k; Wu.omega = Wr.omega;Wu.ke = Wr.ke; Wu.ee = Wr.ee;
	      Wu = HALF*(Wr+SolnBlk.W[i+1][j+1]); Wu.v = SolnBlk.Vwall;
	      Wd = HALF*(Wr+SolnBlk.W[i+1][j-1]); Wd.v = SolnBlk.Vwall;
	    } else if (SolnBlk.Grid.BCtypeW[j] == BC_MOVING_WALL_ISOTHERMAL) {
	      // WEST face of cell (i+1,j) is a MOVINGWALL_ISOTHERMAL boundary.
	      viscous_bc_flag = DIAMONDPATH_RIGHT_TRIANGLE_ISOTHERMAL;
	      //Wu.rho = Wr.p/(Wr.R*SolnBlk.Twall); Wu.v.x = SolnBlk.Vwall.x; Wu.v.y = SolnBlk.Vwall.y; Wu.p = Wr.p; Wu.k = Wr.k; Wu.omega = Wr.omega;Wu.ke = Wr.ke; Wu.ee = Wr.ee;
	      Wu = HALF*(Wr+SolnBlk.W[i+1][j+1]); Wu.v = SolnBlk.Vwall; Wu.rho = Wu.p/(Wr.R*SolnBlk.Twall);
	      Wd = HALF*(Wr+SolnBlk.W[i+1][j-1]); Wd.v = SolnBlk.Vwall; Wd.rho = Wd.p/(Wr.R*SolnBlk.Twall);
	    } else if (SolnBlk.Grid.BCtypeW[j] == BC_BURNING_SURFACE) {
	      // WEST face of cell (i+1,j) is a BURNING_SURFACE boundary.
	      viscous_bc_flag = DIAMONDPATH_RIGHT_TRIANGLE_ISOTHERMAL;
	      Wu = BurningSurface(Wr,SolnBlk.Grid.nfaceW(i+1,j));
	      Wd = Wu;
	    } else if (SolnBlk.Grid.BCtypeW[j] == BC_MASS_INJECTION) {
	      // WEST face of cell (i+1,j) is a MASS_INJECTION boundary.
	      viscous_bc_flag = DIAMONDPATH_RIGHT_TRIANGLE_ISOTHERMAL;
	      //Wu = MassInjection(Wr,SolnBlk.Grid.nfaceW(i+1,j));
	      Wu = MassInjection2(Wr,SolnBlk.Grid.xfaceW(i+1,j),SolnBlk.Grid.nfaceW(i+1,j),SolnBlk.Twall);
	      Wd = Wu;
	    }
	    switch(IP.i_Viscous_Reconstruction) {
	    case VISCOUS_RECONSTRUCTION_CARTESIAN :
	    case VISCOUS_RECONSTRUCTION_DIAMOND_PATH :
	      Xu = SolnBlk.Grid.Node[i+1][j+1].X;
	      Xd = SolnBlk.Grid.Node[i+1][j  ].X;
	      break;
	    case VISCOUS_RECONSTRUCTION_ARITHMETIC_AVERAGE :
	    case VISCOUS_RECONSTRUCTION_HYBRID :
	      Xu = SolnBlk.Grid.xfaceW(i+1,j);
	      Xl = Xu; Wl = Wu;
	      dWdxr = SolnBlk.dWdx[i+1][j]; dWdxl = dWdxr;
	      dWdyr = SolnBlk.dWdy[i+1][j]; dWdyl = dWdyr;
	      break;
	    };
	  } else if (i == SolnBlk.ICu &&
		     (SolnBlk.Grid.BCtypeE[j] == BC_WALL_VISCOUS_HEATFLUX ||
		      SolnBlk.Grid.BCtypeE[j] == BC_WALL_VISCOUS_ISOTHERMAL ||
		      SolnBlk.Grid.BCtypeE[j] == BC_MOVING_WALL_HEATFLUX ||
		      SolnBlk.Grid.BCtypeE[j] == BC_MOVING_WALL_ISOTHERMAL ||
		      SolnBlk.Grid.BCtypeE[j] == BC_BURNING_SURFACE ||
		      SolnBlk.Grid.BCtypeE[j] == BC_MASS_INJECTION)) {
	    // EAST face of cell (i,j) is a normal boundary.
	    Xl = SolnBlk.Grid.Cell[i][j].Xc; Wl = SolnBlk.W[i][j];
	    if (SolnBlk.Grid.BCtypeE[j] == BC_WALL_VISCOUS_HEATFLUX) {
	      // EAST face of cell (i,j) is a WALL_VISCOUS_HEATFLUX boundary.
	      viscous_bc_flag = DIAMONDPATH_LEFT_TRIANGLE_HEATFLUX;
	      //Wu.rho = Wl.rho; Wu.v.x = ZERO; Wu.v.y = ZERO; Wu.p = Wl.p; Wu.k = Wl.k; Wu.omega = Wl.omega;Wu.ke = Wl.ke; Wu.ee = Wl.ee;
	      Wu = HALF*(Wl+SolnBlk.W[i][j+1]); Wu.v = Vector2D_ZERO;
	      Wd = HALF*(Wl+SolnBlk.W[i][j-1]); Wd.v = Vector2D_ZERO;
	    } else if (SolnBlk.Grid.BCtypeE[j] == BC_WALL_VISCOUS_ISOTHERMAL) {
	      // EAST face of cell (i,j) is a WALL_VISCOUS_ISOTHERMAL boundary.
	      viscous_bc_flag = DIAMONDPATH_LEFT_TRIANGLE_ISOTHERMAL;
	      //Wu.rho = Wl.p/(Wl.R*SolnBlk.Twall); Wu.v.x = ZERO; Wu.v.y = ZERO; Wu.p = Wl.p; Wu.k = Wl.k; Wu.omega = Wl.omega;Wu.ke = Wl.ke; Wu.ee = Wl.ee;
	      Wu = HALF*(Wl+SolnBlk.W[i][j+1]); Wu.v = Vector2D_ZERO; Wu.rho = Wu.p/(Wl.R*SolnBlk.Twall);
	      Wd = HALF*(Wl+SolnBlk.W[i][j-1]); Wd.v = Vector2D_ZERO; Wd.rho = Wd.p/(Wl.R*SolnBlk.Twall);
	    } else if (SolnBlk.Grid.BCtypeE[j] == BC_MOVING_WALL_HEATFLUX) {
	      // EAST face of cell (i,j) is a MOVINGWALL_HEATFLUX boundary.
	      viscous_bc_flag = DIAMONDPATH_LEFT_TRIANGLE_HEATFLUX;
	      //Wu.rho = Wl.rho; Wu.v.x = SolnBlk.Vwall.x; Wu.v.y = SolnBlk.Vwall.y; Wu.p = Wl.p; Wu.k = Wl.k; Wu.omega = Wl.omega;Wu.ke = Wl.ke; Wu.ee = Wl.ee;
	      Wu = HALF*(Wl+SolnBlk.W[i][j+1]); Wu.v = SolnBlk.Vwall;
	      Wd = HALF*(Wl+SolnBlk.W[i][j-1]); Wd.v = SolnBlk.Vwall;
	    } else if (SolnBlk.Grid.BCtypeE[j] == BC_MOVING_WALL_ISOTHERMAL) {
	      // EAST face of cell (i,j) is a MOVINGWALL_ISOTHERMAL boundary.
	      viscous_bc_flag = DIAMONDPATH_LEFT_TRIANGLE_ISOTHERMAL;
	      //Wu.rho = Wl.p/(Wl.R*SolnBlk.Twall); Wu.v.x = SolnBlk.Vwall.x; Wu.v.y = SolnBlk.Vwall.y; Wu.p = Wl.p; Wu.k = Wl.k; Wu.omega = Wl.omega;Wu.ke = Wl.ke; Wu.ee = Wl.ee;
	      Wu = HALF*(Wl+SolnBlk.W[i][j+1]); Wu.v = SolnBlk.Vwall; Wu.rho = Wu.p/(Wl.R*SolnBlk.Twall);
	      Wd = HALF*(Wl+SolnBlk.W[i][j-1]); Wd.v = SolnBlk.Vwall; Wd.rho = Wd.p/(Wl.R*SolnBlk.Twall);
	    } else if (SolnBlk.Grid.BCtypeE[j] == BC_MOVING_WALL_ISOTHERMAL) {
	      // EAST face of cell (i,j) is a BURNING_SURFACE boundary.
	      viscous_bc_flag = DIAMONDPATH_LEFT_TRIANGLE_ISOTHERMAL;
	      Wu = BurningSurface(Wr,SolnBlk.Grid.nfaceE(i,j));
	      Wd = Wu;
	    } else if (SolnBlk.Grid.BCtypeE[j] == BC_MASS_INJECTION) {
	      // EAST face of cell (i,j) is a MASS_INJECTION boundary.
	      viscous_bc_flag = DIAMONDPATH_LEFT_TRIANGLE_ISOTHERMAL;
	      //Wu = MassInjection(Wr,SolnBlk.Grid.nfaceE(i,j));
	      Wu = MassInjection2(Wr,SolnBlk.Grid.xfaceE(i,j),SolnBlk.Grid.nfaceE(i,j),SolnBlk.Twall);
	      Wd = Wu;
	    }
	    switch(IP.i_Viscous_Reconstruction) {
	    case VISCOUS_RECONSTRUCTION_CARTESIAN :
	    case VISCOUS_RECONSTRUCTION_DIAMOND_PATH :
	      Xu = SolnBlk.Grid.Node[i+1][j+1].X;
	      Xd = SolnBlk.Grid.Node[i+1][j  ].X;
	      break;
	    case VISCOUS_RECONSTRUCTION_ARITHMETIC_AVERAGE :
	    case VISCOUS_RECONSTRUCTION_HYBRID :
	      Xu = SolnBlk.Grid.xfaceE(i,j);
	      Xr = Xu; Wr = Wu;
	      dWdxl = SolnBlk.dWdx[i][j]; dWdxr = dWdxl;
	      dWdyl = SolnBlk.dWdy[i][j]; dWdyr = dWdyl;
	      break;
	    };
	  } else {
	    // EAST face is either a normal cell or possibly a non-
	    // viscous boundary condition.
	    Xl = SolnBlk.Grid.Cell[i  ][j].Xc; Wl = SolnBlk.W[i  ][j];
	    Xr = SolnBlk.Grid.Cell[i+1][j].Xc; Wr = SolnBlk.W[i+1][j];
	    switch(IP.i_Viscous_Reconstruction) {
	    case VISCOUS_RECONSTRUCTION_CARTESIAN :
	    case VISCOUS_RECONSTRUCTION_DIAMOND_PATH :
	      viscous_bc_flag = DIAMONDPATH_QUADRILATERAL_RECONSTRUCTION;
	      Xu = SolnBlk.Grid.Node[i+1][j+1].X; Wu = SolnBlk.WnNE(i,j);
	      Xd = SolnBlk.Grid.Node[i+1][j  ].X; Wd = SolnBlk.WnSE(i,j);
	      break;
	    case VISCOUS_RECONSTRUCTION_ARITHMETIC_AVERAGE :
	    case VISCOUS_RECONSTRUCTION_HYBRID :
	      Xu = SolnBlk.Grid.xfaceE(i,j); Wu = HALF*(Wl + Wr);
	      dWdxl = SolnBlk.dWdx[i][j]; dWdxr = SolnBlk.dWdy[i+1][j];
	      dWdyl = SolnBlk.dWdy[i][j]; dWdyr = SolnBlk.dWdy[i+1][j];
	      break;
	    };
	  }
	  // Compute the EAST face viscous flux.
	  switch(IP.i_Viscous_Reconstruction) {
	  case VISCOUS_RECONSTRUCTION_CARTESIAN :
	  case VISCOUS_RECONSTRUCTION_DIAMOND_PATH :
	    Flux -= ViscousFluxDiamondPath_n(SolnBlk.Grid.xfaceE(i,j),
					     Xl,Wl,Xd,Wd,Xr,Wr,Xu,Wu,
					     SolnBlk.Grid.nfaceE(i,j),
					     SolnBlk.Axisymmetric,
					     viscous_bc_flag,
					     SolnBlk.Wall[i][j].ywall, 
					     SolnBlk.Wall[i][j].yplus,
					     dWdx,dWdy);
      if (IP.Solver_Type == IMPLICIT && SolnBlk.face_grad_arrays_allocated) {
        SolnBlk.dWdx_faceE[i][j] = dWdx;
        SolnBlk.dWdy_faceE[i][j] = dWdy;
      }
	    break;
	  case VISCOUS_RECONSTRUCTION_ARITHMETIC_AVERAGE :
	    Flux -= ViscousFlux_n(Xu,Wu,HALF*(dWdxl+dWdxr),HALF*(dWdyl+dWdyr),
				  SolnBlk.Grid.nfaceE(i,j),SolnBlk.Axisymmetric,OFF,
				  viscous_bc_flag,
				  SolnBlk.Wall[i][j].ywall, 
				  SolnBlk.Wall[i][j].yplus);
      if (IP.Solver_Type == IMPLICIT && SolnBlk.face_grad_arrays_allocated) {
        SolnBlk.dWdx_faceE[i][j] = HALF*(dWdxl+dWdxr);
        SolnBlk.dWdy_faceE[i][j] = HALF*(dWdyl+dWdyr);
      }
	    break;
	  case VISCOUS_RECONSTRUCTION_HYBRID :
 	    Flux -= ViscousFluxHybrid_n(Xu,Wu,Xl,Wl,dWdxl,dWdyl,
 					Xr,Wr,dWdxr,dWdyr,
 					SolnBlk.Grid.nfaceE(i,j),
 					SolnBlk.Axisymmetric,
 					SolnBlk.Wall[i][j].ywall,
 					SolnBlk.Wall[i][j].yplus,
 					dWdx,dWdy);
      if (IP.Solver_Type == IMPLICIT && SolnBlk.face_grad_arrays_allocated) {
        SolnBlk.dWdx_faceE[i][j] = dWdx;
        SolnBlk.dWdy_faceE[i][j] = dWdy;
      }
	    break;
	  };
	}

	// Evaluate cell-averaged solution changes.
	SolnBlk.dUdt[i  ][j][0] -= Flux*SolnBlk.Grid.lfaceE(i,j)/SolnBlk.Grid.Cell[i][j].A;
	SolnBlk.dUdt[i+1][j][0] += Flux*SolnBlk.Grid.lfaceW(i+1,j)/SolnBlk.Grid.Cell[i+1][j].A;

	// Include axisymmetric source terms if required.
	if (SolnBlk.Axisymmetric) {
	  SolnBlk.dUdt[i][j][0] += SolnBlk.W[i][j].Si(SolnBlk.Grid.Cell[i][j].Xc);
	  if (SolnBlk.Flow_Type)
	    SolnBlk.dUdt[i][j][0] += SolnBlk.W[i][j].Sv(SolnBlk.Grid.Cell[i][j].Xc,
							SolnBlk.dWdy[i][j],
							SolnBlk.Wall[i][j].ywall,
							SolnBlk.Wall[i][j].yplus);
	}

	// Include turbulent production and destruction source term.
	if (SolnBlk.Flow_Type == FLOWTYPE_TURBULENT_RANS_K_OMEGA) {
	  SolnBlk.dUdt[i][j][0] += SolnBlk.W[i][j].St(SolnBlk.Grid.Cell[i][j].Xc,
						      SolnBlk.W[i][j],
						      SolnBlk.dWdx[i][j],
						      SolnBlk.dWdy[i][j],
						      SolnBlk.Axisymmetric,
						      SolnBlk.Wall[i][j].ywall,
						      SolnBlk.Wall[i][j].yplus);
	}

	// Save west and east face boundary flux.
	if (i == SolnBlk.ICl-1) {
	  SolnBlk.FluxW[j] = -Flux*SolnBlk.Grid.lfaceW(i+1,j);
	} else if (i == SolnBlk.ICu) {
	  SolnBlk.FluxE[j] =  Flux*SolnBlk.Grid.lfaceE(i,j);
	}

      }
    }

    if (j > SolnBlk.JCl-1 && j < SolnBlk.JCu+1) {
      SolnBlk.dUdt[SolnBlk.ICl-1][j][0] = NavierStokes2D_U_VACUUM;
      SolnBlk.dUdt[SolnBlk.ICu+1][j][0] = NavierStokes2D_U_VACUUM;
    }

  }

  // Add j-direction (eta-direction) fluxes.
  for (int i = SolnBlk.ICl; i <= SolnBlk.ICu; i++) {
    for (int j = SolnBlk.JCl-1; j <= SolnBlk.JCu; j++) {

      // Evaluate the cell interface j-direction INVISCID fluxes.
      if (j == SolnBlk.JCl-1 && 
	  (SolnBlk.Grid.BCtypeS[i] == BC_REFLECTION ||
	   SolnBlk.Grid.BCtypeS[i] == BC_BURNING_SURFACE ||
	   SolnBlk.Grid.BCtypeS[i] == BC_MASS_INJECTION ||
	   SolnBlk.Grid.BCtypeS[i] == BC_WALL_VISCOUS_HEATFLUX ||
	   SolnBlk.Grid.BCtypeS[i] == BC_WALL_VISCOUS_ISOTHERMAL ||
	   SolnBlk.Grid.BCtypeS[i] == BC_MOVING_WALL_HEATFLUX ||
	   SolnBlk.Grid.BCtypeS[i] == BC_MOVING_WALL_ISOTHERMAL ||
	   SolnBlk.Grid.BCtypeS[i] == BC_RINGLEB_FLOW ||
	   SolnBlk.Grid.BCtypeS[i] == BC_CHARACTERISTIC)) {

	dX = SolnBlk.Grid.xfaceS(i,j+1) - SolnBlk.Grid.Cell[i][j+1].Xc;
	Wr = SolnBlk.W[i][j+1] + (SolnBlk.phi[i][j+1]^SolnBlk.dWdx[i][j+1])*dX.x +
                                 (SolnBlk.phi[i][j+1]^SolnBlk.dWdy[i][j+1])*dX.y;

	// SOUTH face of cell (i,j+1) is a normal boundary.
	if (SolnBlk.Grid.BCtypeS[i] == BC_REFLECTION) {
	  // SOUTH face of cell (i,j+1) is a REFLECTION boundary.
	  Wl = Reflect(Wr,SolnBlk.Grid.nfaceS(i,j+1));
	} else if (SolnBlk.Grid.BCtypeS[i] == BC_BURNING_SURFACE) {
	  // SOUTH face of cell (i,j+1) is a BURNING_SURFACE boundary.
	  Wl = BurningSurface(Wr,SolnBlk.Grid.nfaceS(i,j+1));
	} else if (SolnBlk.Grid.BCtypeS[i] == BC_MASS_INJECTION) {
	  // SOUTH face of cell (i,j+1) is a MASS_INJECTION boundary.
	  //Wl = MassInjection(Wr,SolnBlk.Grid.nfaceS(i,j+1));
	  Wl = MassInjection2(Wr,SolnBlk.Grid.xfaceS(i,j+1),SolnBlk.Grid.nfaceS(i,j+1),SolnBlk.Twall);
	  flux_function = IP.i_Flux_Function;
	  IP.i_Flux_Function = FLUX_FUNCTION_GODUNOV_WRS;
	} else if (SolnBlk.Grid.BCtypeS[i] == BC_WALL_VISCOUS_HEATFLUX) {
	  // SOUTH face of cell (i,j+1) is a WALL_VISCOUS_HEATFLUX boundary.
	  Wl = WallViscousHeatFlux(Wr,SolnBlk.Grid.nfaceS(i,j+1));
	} else if (SolnBlk.Grid.BCtypeS[i] == BC_WALL_VISCOUS_ISOTHERMAL) {
	  // SOUTH face of cell (i,j+1) is a WALL_VISCOUS_ISOTHERMAL boundary.
	  Wl = WallViscousIsothermal(Wr,SolnBlk.Grid.nfaceS(i,j+1),SolnBlk.Twall);
	} else if (SolnBlk.Grid.BCtypeS[i] == BC_MOVING_WALL_HEATFLUX) {
	  // SOUTH face of cell (i,j+1) is a MOVINGWALL_HEATFLUX boundary.
	  Wl = MovingWallHeatFlux(Wr,SolnBlk.Grid.nfaceS(i,j+1),SolnBlk.Vwall.x);
	} else if (SolnBlk.Grid.BCtypeS[i] == BC_MOVING_WALL_ISOTHERMAL) {
	  // SOUTH face of cell (i,j+1) is a MOVINGWALL_ISOTHERMAL boundary.
	  Wl = MovingWallIsothermal(Wr,SolnBlk.Grid.nfaceS(i,j+1),SolnBlk.Vwall.x,SolnBlk.Twall);
	} else if (SolnBlk.Grid.BCtypeS[i] == BC_RINGLEB_FLOW) {
	  // SOUTH face of cell (i,j+1) is a RINGLEB_FLOW boundary.
	  Wl = RinglebFlow(Wl,SolnBlk.Grid.xfaceS(i,j+1));
	  Wl = BC_Characteristic_Pressure(Wr,Wl,SolnBlk.Grid.nfaceS(i,j+1));
	} else if (SolnBlk.Grid.BCtypeS[i] == BC_CHARACTERISTIC) {
	  // SOUTH face of cell (i,j+1) is a CHARACTERISTIC boundary.
	  Wl = BC_Characteristic_Pressure(Wr,SolnBlk.WoS[i],SolnBlk.Grid.nfaceS(i,j+1));
	}

      } else if (j == SolnBlk.JCu && 
		 (SolnBlk.Grid.BCtypeN[i] == BC_REFLECTION ||
		  SolnBlk.Grid.BCtypeN[i] == BC_BURNING_SURFACE ||
		  SolnBlk.Grid.BCtypeN[i] == BC_MASS_INJECTION ||
		  SolnBlk.Grid.BCtypeN[i] == BC_WALL_VISCOUS_HEATFLUX ||
		  SolnBlk.Grid.BCtypeN[i] == BC_WALL_VISCOUS_ISOTHERMAL ||
		  SolnBlk.Grid.BCtypeN[i] == BC_MOVING_WALL_HEATFLUX ||
		  SolnBlk.Grid.BCtypeN[i] == BC_MOVING_WALL_ISOTHERMAL ||
		  SolnBlk.Grid.BCtypeN[i] == BC_RINGLEB_FLOW ||
		  SolnBlk.Grid.BCtypeN[i] == BC_CHARACTERISTIC)) {

	dX = SolnBlk.Grid.xfaceN(i,j) - SolnBlk.Grid.Cell[i][j].Xc;
	Wl = SolnBlk.W[i][j] + (SolnBlk.phi[i][j]^SolnBlk.dWdx[i][j])*dX.x +
	                       (SolnBlk.phi[i][j]^SolnBlk.dWdy[i][j])*dX.y;

	// NORTH face of cell (i,j) is a normal boundary.
	if (SolnBlk.Grid.BCtypeN[i] == BC_REFLECTION) {
	  // NORTH face of cell (i,j) is a REFLECTION boundary.
	  Wr = Reflect(Wl,SolnBlk.Grid.nfaceN(i,j));
	} else if (SolnBlk.Grid.BCtypeN[i] == BC_BURNING_SURFACE) {
	  // NORTH face of cell (i,j) is a BURNING_SURFACE boundary.
	  Wr = BurningSurface(Wl,SolnBlk.Grid.nfaceN(i,j));
	} else if (SolnBlk.Grid.BCtypeN[i] == BC_MASS_INJECTION) {
	  // NORTH face of cell (i,j) is a MASS_INJECTION boundary.
	  //Wr = MassInjection(Wl,SolnBlk.Grid.nfaceN(i,j));
	  Wr = MassInjection2(Wl,SolnBlk.Grid.xfaceN(i,j),SolnBlk.Grid.nfaceN(i,j),SolnBlk.Twall);
	  flux_function = IP.i_Flux_Function;
	  IP.i_Flux_Function = FLUX_FUNCTION_GODUNOV_WRS;
	} else if (SolnBlk.Grid.BCtypeN[i] == BC_WALL_VISCOUS_HEATFLUX) {
	  // NORTH face of cell (i,j) is a WALL_VISCOUS_HEATFLUX boundary.
	  Wr = WallViscousHeatFlux(Wl,SolnBlk.Grid.nfaceN(i,j));
	} else if (SolnBlk.Grid.BCtypeN[i] == BC_WALL_VISCOUS_ISOTHERMAL) {
	  // NORTH face of cell (i,j) is a WALL_VISCOUS_ISOTHERMAL boundary.
	  Wr = WallViscousIsothermal(Wl,SolnBlk.Grid.nfaceN(i,j),SolnBlk.Twall);
	} else if (SolnBlk.Grid.BCtypeN[i] == BC_MOVING_WALL_HEATFLUX) {
	  // NORTH face of cell (i,j) is a MOVINGWALL_HEATFLUX boundary.
	  Wr = MovingWallHeatFlux(Wl,SolnBlk.Grid.nfaceN(i,j),SolnBlk.Vwall.x);
	} else if (SolnBlk.Grid.BCtypeN[i] == BC_MOVING_WALL_ISOTHERMAL) {
	  // NORTH face of cell (i,j) is a MOVINGWALL_ISOTHERMAL boundary.
	  Wr = MovingWallIsothermal(Wl,SolnBlk.Grid.nfaceN(i,j),SolnBlk.Vwall.x,SolnBlk.Twall);
	} else if (SolnBlk.Grid.BCtypeN[i] == BC_RINGLEB_FLOW) {
	  // NORTH face of cell (i,j) is a RINGLEB_FLOW boundary.
	  Wr = RinglebFlow(Wr,SolnBlk.Grid.xfaceN(i,j));
	  Wr = BC_Characteristic_Pressure(Wl,Wr,SolnBlk.Grid.nfaceN(i,j));
	} else {
	  // NORTH face of cell (i,j) is a CHARACTERISTIC boundary.
	  Wr = BC_Characteristic_Pressure(Wl,SolnBlk.WoN[i],SolnBlk.Grid.nfaceN(i,j));
	}

      } else {

	// NORTH face is either a normal cell or possibly a FIXED, 
	// NONE or EXTRAPOLATION boundary.
	dX = SolnBlk.Grid.xfaceN(i,j  ) - SolnBlk.Grid.Cell[i][j  ].Xc;
	Wl = SolnBlk.W[i][j  ] + (SolnBlk.phi[i][j  ]^SolnBlk.dWdx[i][j  ])*dX.x +
                                 (SolnBlk.phi[i][j  ]^SolnBlk.dWdy[i][j  ])*dX.y;
	dX = SolnBlk.Grid.xfaceS(i,j+1) - SolnBlk.Grid.Cell[i][j+1].Xc;
	Wr = SolnBlk.W[i][j+1] + (SolnBlk.phi[i][j+1]^SolnBlk.dWdx[i][j+1])*dX.x +
                                 (SolnBlk.phi[i][j+1]^SolnBlk.dWdy[i][j+1])*dX.y;

      }

      // Determine NORTH face inviscid flux.
      switch(IP.i_Flux_Function) {
      case FLUX_FUNCTION_GODUNOV :
	Flux = FluxGodunov_n(Wl,Wr,SolnBlk.Grid.nfaceN(i,j));
	break;
      case FLUX_FUNCTION_ROE :
	Flux = FluxRoe_n(Wl,Wr,SolnBlk.Grid.nfaceN(i,j));
	break;
      case FLUX_FUNCTION_RUSANOV :
	Flux = FluxRusanov_n(Wl,Wr,SolnBlk.Grid.nfaceN(i,j));
	break;
      case FLUX_FUNCTION_HLLE :
	Flux = FluxHLLE_n(Wl,Wr,SolnBlk.Grid.nfaceN(i,j));
	break;
      case FLUX_FUNCTION_HLLL :
	Flux = FluxHLLL_n(Wl,Wr,SolnBlk.Grid.nfaceN(i,j));
	break;
      case FLUX_FUNCTION_HLLC :
	Flux = FluxHLLC_n(Wl,Wr,SolnBlk.Grid.nfaceN(i,j));
	break;
      case FLUX_FUNCTION_VANLEER :
	Flux = FluxVanLeer_n(Wl,Wr,SolnBlk.Grid.nfaceN(i,j));
	break;
      case FLUX_FUNCTION_AUSM :
	Flux = FluxAUSM_n(Wl,Wr,SolnBlk.Grid.nfaceN(i,j));
	break;
      case FLUX_FUNCTION_AUSMplus :
	Flux = FluxAUSMplus_n(Wl,Wr,SolnBlk.Grid.nfaceN(i,j));
	break;
      case FLUX_FUNCTION_GODUNOV_WRS :
	Flux = FluxGodunov_Wrs_n(Wl,Wr,SolnBlk.Grid.nfaceN(i,j));
	break;
      default:
	Flux = FluxRoe_n(Wl,Wr,SolnBlk.Grid.nfaceN(i,j));
	break;
      };

      if (IP.i_Flux_Function == FLUX_FUNCTION_GODUNOV_WRS) {
	IP.i_Flux_Function = flux_function;
      }

      // Evaluate the cell interface j-direction VISCOUS flux if necessary.
      if (SolnBlk.Flow_Type) {
	// Determine the NORTH face VISCOUS flux.
	if (j == SolnBlk.JCl-1 && 
	    (SolnBlk.Grid.BCtypeS[i] == BC_WALL_VISCOUS_HEATFLUX ||
	     SolnBlk.Grid.BCtypeS[i] == BC_WALL_VISCOUS_ISOTHERMAL ||
	     SolnBlk.Grid.BCtypeS[i] == BC_MOVING_WALL_HEATFLUX ||
	     SolnBlk.Grid.BCtypeS[i] == BC_MOVING_WALL_ISOTHERMAL ||
	     SolnBlk.Grid.BCtypeS[i] == BC_BURNING_SURFACE ||
	     SolnBlk.Grid.BCtypeS[i] == BC_MASS_INJECTION)) {
	  // SOUTH face of cell (i,j+1) is a normal boundary.
	  Xr = SolnBlk.Grid.Cell[i][j+1].Xc; Wr = SolnBlk.W[i][j+1];
	  if (SolnBlk.Grid.BCtypeS[i] == BC_WALL_VISCOUS_HEATFLUX) {
	    // SOUTH face of cell (i,j+1) is a WALL_VISCOUS_HEATFLUX boundary.
	    viscous_bc_flag = DIAMONDPATH_RIGHT_TRIANGLE_HEATFLUX;
	    //Wu.rho = Wr.rho; Wu.v.x = ZERO; Wu.v.y = ZERO; Wu.p = Wr.p; Wu.k = Wr.k; Wu.omega = Wr.omega;Wu.ke = Wr.ke; Wu.ee = Wr.ee;
	    Wu = HALF*(Wr+SolnBlk.W[i-1][j+1]); Wu.v = Vector2D_ZERO;
	    Wd = HALF*(Wr+SolnBlk.W[i+1][j+1]); Wd.v = Vector2D_ZERO;
	  } else if (SolnBlk.Grid.BCtypeS[i] == BC_WALL_VISCOUS_ISOTHERMAL) {
	    // SOUTH face of cell (i,j+1) is a WALL_VISCOUS_ISOTHERMAL boundary.
	    viscous_bc_flag = DIAMONDPATH_RIGHT_TRIANGLE_ISOTHERMAL;
	    //Wu.rho = Wr.p/(Wr.R*SolnBlk.Twall); Wu.v.x = ZERO; Wu.v.y = ZERO; Wu.p = Wr.p; Wu.k = Wr.k; Wu.omega = Wr.omega;Wu.ke = Wr.ke; Wu.ee = Wr.ee;
	    Wu = HALF*(Wr+SolnBlk.W[i-1][j+1]); Wu.v = Vector2D_ZERO; Wu.rho = Wu.p/(Wr.R*SolnBlk.Twall);
	    Wd = HALF*(Wr+SolnBlk.W[i+1][j+1]); Wd.v = Vector2D_ZERO; Wd.rho = Wd.p/(Wr.R*SolnBlk.Twall);
	  } else if (SolnBlk.Grid.BCtypeS[i] == BC_MOVING_WALL_HEATFLUX) {
	    // SOUTH face of cell (i,j+1) is a MOVINGWALL_HEATFLUX boundary.
	    viscous_bc_flag = DIAMONDPATH_RIGHT_TRIANGLE_HEATFLUX;
	    //Wu.rho = Wr.rho; Wu.v.x = SolnBlk.Vwall.x; Wu.v.y = SolnBlk.Vwall.y; Wu.p = Wr.p; Wu.k = Wr.k; Wu.omega = Wr.omega;Wu.ke = Wr.ke; Wu.ee = Wr.ee;
	    Wu = HALF*(Wr+SolnBlk.W[i-1][j+1]); Wu.v = SolnBlk.Vwall;
	    Wd = HALF*(Wr+SolnBlk.W[i+1][j+1]); Wd.v = SolnBlk.Vwall;
	  } else if (SolnBlk.Grid.BCtypeS[i] == BC_MOVING_WALL_ISOTHERMAL) {
	    // SOUTH face of cell (i,j+1) is a MOVINGWALL_ISOTHERMAL boundary.
	    viscous_bc_flag = DIAMONDPATH_RIGHT_TRIANGLE_ISOTHERMAL;
	    //Wu.rho = Wr.p/(Wr.R*SolnBlk.Twall); Wu.v.x = SolnBlk.Vwall.x; Wu.v.y = SolnBlk.Vwall.y; Wu.p = Wr.p; Wu.k = Wr.k; Wu.omega = Wr.omega;Wu.ke = Wr.ke; Wu.ee = Wr.ee;
	    Wu = HALF*(Wr+SolnBlk.W[i-1][j+1]); Wu.v = SolnBlk.Vwall; Wu.rho = Wu.p/(Wr.R*SolnBlk.Twall);
	    Wd = HALF*(Wr+SolnBlk.W[i+1][j+1]); Wd.v = SolnBlk.Vwall; Wd.rho = Wd.p/(Wr.R*SolnBlk.Twall);
	  } else if (SolnBlk.Grid.BCtypeS[i] == BC_BURNING_SURFACE) {
	    // SOUTH face of cell (i,j+1) is a BURNING_SURFACE boundary.
	    viscous_bc_flag = DIAMONDPATH_RIGHT_TRIANGLE_ISOTHERMAL;
	    Wu = BurningSurface(Wr,SolnBlk.Grid.nfaceS(i,j+1));
	    Wd = Wu;
	  } else if (SolnBlk.Grid.BCtypeS[i] == BC_MASS_INJECTION) {
	    // SOUTH face of cell (i,j+1) is a MASS_INJECTION boundary.
	    viscous_bc_flag = DIAMONDPATH_RIGHT_TRIANGLE_ISOTHERMAL;
	    //Wu = MassInjection(Wr,SolnBlk.Grid.nfaceS(i,j+1));
	    Wu = MassInjection2(Wr,SolnBlk.Grid.xfaceS(i,j+1),SolnBlk.Grid.nfaceS(i,j+1),SolnBlk.Twall);
	    Wd = Wu;
	  }
	  switch(IP.i_Viscous_Reconstruction) {
	  case VISCOUS_RECONSTRUCTION_CARTESIAN :
	  case VISCOUS_RECONSTRUCTION_DIAMOND_PATH :
	    Xu = SolnBlk.Grid.Node[i  ][j+1].X;
	    Xd = SolnBlk.Grid.Node[i+1][j+1].X;
	    break;
	  case VISCOUS_RECONSTRUCTION_ARITHMETIC_AVERAGE :
	  case VISCOUS_RECONSTRUCTION_HYBRID :
	    Xu = SolnBlk.Grid.xfaceS(i,j+1);
	    Xl = Xu; Wl = Wu;
	    dWdxr = SolnBlk.dWdx[i][j+1]; dWdxl = dWdxr;
	    dWdyr = SolnBlk.dWdy[i][j+1]; dWdyl = dWdyr;
	    break;
	  };
	} else if (j == SolnBlk.JCu && 
		   (SolnBlk.Grid.BCtypeN[i] == BC_WALL_VISCOUS_HEATFLUX ||
		    SolnBlk.Grid.BCtypeN[i] == BC_WALL_VISCOUS_ISOTHERMAL ||
		    SolnBlk.Grid.BCtypeN[i] == BC_MOVING_WALL_HEATFLUX ||
		    SolnBlk.Grid.BCtypeN[i] == BC_MOVING_WALL_ISOTHERMAL ||
		    SolnBlk.Grid.BCtypeN[i] == BC_BURNING_SURFACE ||
		    SolnBlk.Grid.BCtypeN[i] == BC_MASS_INJECTION)) {
	  // NORTH face of cell (i,j) is a normal boundary.
	  Xl = SolnBlk.Grid.Cell[i][j].Xc; Wl = SolnBlk.W[i][j];
	  if (SolnBlk.Grid.BCtypeN[i] == BC_WALL_VISCOUS_HEATFLUX) {
	    // NORTH face of cell (i,j) is a WALL_VISCOUS_HEATFLUX boundary.
	    viscous_bc_flag = DIAMONDPATH_LEFT_TRIANGLE_HEATFLUX;
	    //Wu.rho = Wl.rho; Wu.v.x = ZERO; Wu.v.y = ZERO; Wu.p = Wl.p; Wu.k = Wl.k; Wu.omega = Wl.omega;Wu.ke = Wl.ke; Wu.ee = Wl.ee;
	    Wu = HALF*(Wl+SolnBlk.W[i-1][j]); Wu.v = Vector2D_ZERO;
	    Wd = HALF*(Wl+SolnBlk.W[i+1][j]); Wd.v = Vector2D_ZERO;
	  } else if (SolnBlk.Grid.BCtypeN[i] == BC_WALL_VISCOUS_ISOTHERMAL) {
	    // NORTH face of cell (i,j) is a WALL_VISCOUS_ISOTHERMAL boundary.
	    viscous_bc_flag = DIAMONDPATH_LEFT_TRIANGLE_ISOTHERMAL;
	    //Wu.rho = Wl.p/(Wl.R*SolnBlk.Twall); Wu.v.x = ZERO; Wu.v.y = ZERO; Wu.p = Wl.p; Wu.k = Wl.k; Wu.omega = Wl.omega;Wu.ke = Wl.ke; Wu.ee = Wl.ee;
	    Wu = HALF*(Wl+SolnBlk.W[i-1][j]); Wu.v = Vector2D_ZERO; Wu.rho = Wu.p/(Wl.R*SolnBlk.Twall);
	    Wd = HALF*(Wl+SolnBlk.W[i+1][j]); Wd.v = Vector2D_ZERO; Wd.rho = Wd.p/(Wl.R*SolnBlk.Twall);
	  } else if (SolnBlk.Grid.BCtypeN[i] == BC_MOVING_WALL_HEATFLUX) {
	    // NORTH face of cell (i,j) is a MOVINGWALL_HEATFLUX boundary.
	    viscous_bc_flag = DIAMONDPATH_LEFT_TRIANGLE_HEATFLUX;
	    //Wu.rho = Wl.rho; Wu.v.x = SolnBlk.Vwall.x; Wu.v.y = SolnBlk.Vwall.y; Wu.p = Wl.p; Wu.k = Wl.k; Wu.omega = Wl.omega;Wu.ke = Wl.ke; Wu.ee = Wl.ee;
	    Wu = HALF*(Wl+SolnBlk.W[i-1][j]); Wu.v = SolnBlk.Vwall;
	    Wd = HALF*(Wl+SolnBlk.W[i+1][j]); Wd.v = SolnBlk.Vwall;
	  } else if (SolnBlk.Grid.BCtypeN[i] == BC_MOVING_WALL_ISOTHERMAL) {
	    // NORTH face of cell (i,j) is a MOVINGWALL_ISOTHERMAL boundary.
	    viscous_bc_flag = DIAMONDPATH_LEFT_TRIANGLE_ISOTHERMAL;
	    //Wu.rho = Wl.p/(Wl.R*SolnBlk.Twall); Wu.v.x = SolnBlk.Vwall.x; Wu.v.y = SolnBlk.Vwall.y; Wu.p = Wl.p; Wu.k = Wl.k; Wu.omega = Wl.omega;Wu.ke = Wl.ke; Wu.ee = Wl.ee;
	    Wu = HALF*(Wl+SolnBlk.W[i-1][j]); Wu.v = SolnBlk.Vwall; Wu.rho = Wu.p/(Wl.R*SolnBlk.Twall);
	    Wd = HALF*(Wl+SolnBlk.W[i+1][j]); Wd.v = SolnBlk.Vwall; Wd.rho = Wd.p/(Wl.R*SolnBlk.Twall);
	  } else if (SolnBlk.Grid.BCtypeN[i] == BC_BURNING_SURFACE) {
	    // NORTH face of cell (i,j) is a BURNING_SURFACE boundary.
	    viscous_bc_flag = DIAMONDPATH_LEFT_TRIANGLE_ISOTHERMAL;
	    Wu = BurningSurface(Wl,SolnBlk.Grid.nfaceN(i,j));
	    Wd = Wu;
	  } else if (SolnBlk.Grid.BCtypeN[i] == BC_MASS_INJECTION) {
	    // NORTH face of cell (i,j) is a MASS_INJECTION boundary.
	    viscous_bc_flag = DIAMONDPATH_LEFT_TRIANGLE_ISOTHERMAL;
	    //Wu = MassInjection(Wl,SolnBlk.Grid.nfaceN(i,j));
	    Wu = MassInjection2(Wl,SolnBlk.Grid.xfaceN(i,j),SolnBlk.Grid.nfaceN(i,j),SolnBlk.Twall);
	    Wd = Wu;
	  }
	  switch(IP.i_Viscous_Reconstruction) {
	  case VISCOUS_RECONSTRUCTION_CARTESIAN :
	  case VISCOUS_RECONSTRUCTION_DIAMOND_PATH :
	    Xu = SolnBlk.Grid.Node[i  ][j+1].X;
	    Xd = SolnBlk.Grid.Node[i+1][j+1].X;
	    break;
	  case VISCOUS_RECONSTRUCTION_ARITHMETIC_AVERAGE :
	  case VISCOUS_RECONSTRUCTION_HYBRID :
	    Xu = SolnBlk.Grid.xfaceN(i,j);
	    Xr = Xu; Wr = Wu;
	    dWdxl = SolnBlk.dWdx[i][j]; dWdxr = dWdxl;
	    dWdyl = SolnBlk.dWdy[i][j]; dWdyr = dWdyl;
	    break;
	  };
	} else {
	  // NORTH face is either a normal cell or possibly a non-viscous
	  // boundary condition.
	  Xl = SolnBlk.Grid.Cell[i][j  ].Xc; Wl = SolnBlk.W[i][j  ];
	  Xr = SolnBlk.Grid.Cell[i][j+1].Xc; Wr = SolnBlk.W[i][j+1];
	  switch(IP.i_Viscous_Reconstruction) {
	  case VISCOUS_RECONSTRUCTION_CARTESIAN :
	  case VISCOUS_RECONSTRUCTION_DIAMOND_PATH :
	    viscous_bc_flag = DIAMONDPATH_QUADRILATERAL_RECONSTRUCTION;
	    Xu = SolnBlk.Grid.Node[i  ][j+1].X; Wu = SolnBlk.WnNW(i,j);
	    Xd = SolnBlk.Grid.Node[i+1][j+1].X; Wd = SolnBlk.WnNE(i,j);
	    break;
	  case VISCOUS_RECONSTRUCTION_ARITHMETIC_AVERAGE :
	  case VISCOUS_RECONSTRUCTION_HYBRID :
	    Xu = SolnBlk.Grid.xfaceN(i,j); Wu = HALF*(Wl + Wr);
	    dWdxl = SolnBlk.dWdx[i][j]; dWdxr = SolnBlk.dWdy[i][j+1];
	    dWdyl = SolnBlk.dWdy[i][j]; dWdyr = SolnBlk.dWdy[i][j+1];
	    break;
	  };
	}
	// Compute the NORTH face viscous flux.
	switch(IP.i_Viscous_Reconstruction) {
	case VISCOUS_RECONSTRUCTION_CARTESIAN :
	case VISCOUS_RECONSTRUCTION_DIAMOND_PATH :
	  Flux -= ViscousFluxDiamondPath_n(SolnBlk.Grid.xfaceN(i,j),
					   Xl,Wl,Xd,Wd,Xr,Wr,Xu,Wu,
					   SolnBlk.Grid.nfaceN(i,j),
					   SolnBlk.Axisymmetric,
					   viscous_bc_flag,
					   SolnBlk.Wall[i][j].ywall,
					   SolnBlk.Wall[i][j].yplus,
					   dWdx,dWdy);
      if (IP.Solver_Type == IMPLICIT && SolnBlk.face_grad_arrays_allocated) {
        SolnBlk.dWdx_faceN[i][j] = dWdx;
        SolnBlk.dWdy_faceN[i][j] = dWdy;
      }
	  break;
	case VISCOUS_RECONSTRUCTION_ARITHMETIC_AVERAGE :
	  Flux -= ViscousFlux_n(Xu,Wu,HALF*(dWdxl+dWdxr),HALF*(dWdyl+dWdyr),
				SolnBlk.Grid.nfaceN(i,j),SolnBlk.Axisymmetric,OFF,
				SolnBlk.Wall[i][j].ywall,
				SolnBlk.Wall[i][j].yplus);
      if (IP.Solver_Type == IMPLICIT && SolnBlk.face_grad_arrays_allocated) {
        SolnBlk.dWdx_faceN[i][j] = HALF*(dWdxl+dWdxr);
        SolnBlk.dWdy_faceN[i][j] = HALF*(dWdyl+dWdyr);
      }
	  break;
	case VISCOUS_RECONSTRUCTION_HYBRID :
 	  Flux -= ViscousFluxHybrid_n(Xu,Wu,Xl,Wl,dWdxl,dWdyl,
 				      Xr,Wr,dWdxr,dWdyr,
 				      SolnBlk.Grid.nfaceN(i,j),
 				      SolnBlk.Axisymmetric,
 				      SolnBlk.Wall[i][j].ywall,
 				      SolnBlk.Wall[i][j].yplus,
 				      dWdx,dWdy);
      if (IP.Solver_Type == IMPLICIT && SolnBlk.face_grad_arrays_allocated) {
        SolnBlk.dWdx_faceN[i][j] = dWdx;
        SolnBlk.dWdy_faceN[i][j] = dWdy;
      }
	  break;
	};
      }

      // Evaluate cell-averaged solution changes.
      SolnBlk.dUdt[i][j  ][0] -= Flux*SolnBlk.Grid.lfaceN(i,j)/SolnBlk.Grid.Cell[i][j].A;
      SolnBlk.dUdt[i][j+1][0] += Flux*SolnBlk.Grid.lfaceS(i,j+1)/SolnBlk.Grid.Cell[i][j+1].A;

      // Save south and north face boundary flux.
      if (j == SolnBlk.JCl-1) {
	SolnBlk.FluxS[i] = -Flux*SolnBlk.Grid.lfaceS(i,j+1);
      } else if (j == SolnBlk.JCu) {
	SolnBlk.FluxN[i] = Flux*SolnBlk.Grid.lfaceN(i,j);
      }

    }

    SolnBlk.dUdt[i][SolnBlk.JCl-1][0] = NavierStokes2D_U_VACUUM;
    SolnBlk.dUdt[i][SolnBlk.JCu+1][0] = NavierStokes2D_U_VACUUM;

  }

  // Zero the residuals for the turbulence variables according to
  // the turbulence boundary condition.
  error_flag = Turbulence_Zero_Residual(SolnBlk,1,IP);
  if (error_flag) return error_flag;

  // Residual successfully evaluated.
  return 0;

}

/**********************************************************************
 * Routine: dUdt_Multistage_Explicit                                  *
 *                                                                    *
 * This routine determines the solution residuals for a given stage   *
 * of a variety of multi-stage explicit time integration schemes for  *
 * a given solution block.                                            *
 *                                                                    *
 **********************************************************************/
int dUdt_Multistage_Explicit(NavierStokes2D_Quad_Block &SolnBlk,
                             const int i_stage,
			     NavierStokes2D_Input_Parameters &IP) {

  int error_flag, k_residual;
  double omega;
  Vector2D dX;
  NavierStokes2D_pState Wl, Wr;
  NavierStokes2D_cState Flux;

  NavierStokes2D_pState Wu, Wd, dWdxl, dWdyl, dWdxr, dWdyr;
  Vector2D Xl, Xr, Xu, Xd;
  int viscous_bc_flag;

  int flux_function = IP.i_Flux_Function;

  NavierStokes2D_cState NavierStokes2D_U_VACUUM; NavierStokes2D_U_VACUUM.Vacuum();
  NavierStokes2D_pState NavierStokes2D_W_STDATM; NavierStokes2D_W_STDATM.Standard_Atmosphere();

  // Evaluate the solution residual for stage i_stage of n_stage scheme.

  // Evaluate the time step fraction and residual storage location for
  // the stage.
  switch(IP.i_Time_Integration) {
  case TIME_STEPPING_EXPLICIT_EULER :
    omega = Runge_Kutta(i_stage,IP.N_Stage);
    k_residual = 0;
    break;
  case TIME_STEPPING_EXPLICIT_PREDICTOR_CORRECTOR :
    omega = Runge_Kutta(i_stage,IP.N_Stage);
    k_residual = 0;
    break;
  case TIME_STEPPING_EXPLICIT_RUNGE_KUTTA :
    omega = Runge_Kutta(i_stage,IP.N_Stage);
    k_residual = 0;
    if (IP.N_Stage == 4) {
      if (i_stage == 4) k_residual = 0;
      else k_residual = i_stage - 1;
    }
    break;
  case TIME_STEPPING_MULTISTAGE_OPTIMAL_SMOOTHING :
    omega = MultiStage_Optimally_Smoothing(i_stage,IP.N_Stage,IP.i_Limiter);
    k_residual = 0;
    break;
  default:
    omega = Runge_Kutta(i_stage,IP.N_Stage);
    k_residual = 0;
    break;
  };

  // Perform the linear reconstruction within each cell of the
  // computational grid for this stage.
  switch(IP.i_Reconstruction) {
  case RECONSTRUCTION_GREEN_GAUSS :
    Linear_Reconstruction_GreenGauss(SolnBlk,IP.i_Limiter);
    break;
  case RECONSTRUCTION_LINEAR_LEAST_SQUARES :
    Linear_Reconstruction_LeastSquares(SolnBlk,IP.i_Limiter);
    break;
  default:
    Linear_Reconstruction_LeastSquares(SolnBlk,IP.i_Limiter);
    break;
  };

  // Evaluate the time rate of change of the solution (i.e., the
  // solution residuals) using a second-order limited upwind scheme
  // with a variety of flux functions.

  // Add i-direction (zeta-direction) fluxes.
  for (int j = SolnBlk.JCl-1; j <= SolnBlk.JCu+1; j++) {
    if (i_stage == 1) {
      SolnBlk.Uo[SolnBlk.ICl-1][j] = SolnBlk.U[SolnBlk.ICl-1][j];
      SolnBlk.dUdt[SolnBlk.ICl-1][j][k_residual] = NavierStokes2D_U_VACUUM;
    } else {
      SolnBlk.dUdt[SolnBlk.ICl-1][j][k_residual] = NavierStokes2D_U_VACUUM;
    }

    for (int i = SolnBlk.ICl-1; i <= SolnBlk.ICu; i++) {
      if (i_stage == 1) {
	SolnBlk.Uo[i+1][j] = SolnBlk.U[i+1][j];
	SolnBlk.dUdt[i+1][j][k_residual] = NavierStokes2D_U_VACUUM;
      } else if (j > SolnBlk.JCl-1 && j < SolnBlk.JCu+1) {
	switch(IP.i_Time_Integration) {
	case TIME_STEPPING_EXPLICIT_PREDICTOR_CORRECTOR :
	  //SolnBlk.dUdt[i+1][j][k_residual] = SolnBlk.dUdt[i+1][j][k_residual];
	  break;
	case TIME_STEPPING_EXPLICIT_RUNGE_KUTTA :
	  if (IP.N_Stage == 2) {
	    //SolnBlk.dUdt[i+1][j][k_residual] = SolnBlk.dUdt[i+1][j][k_residual];
	  } else if (IP.N_Stage == 4 && i_stage == 4) {
	    SolnBlk.dUdt[i+1][j][k_residual] = SolnBlk.dUdt[i+1][j][0] + 
                                               TWO*SolnBlk.dUdt[i+1][j][1] +
                                               TWO*SolnBlk.dUdt[i+1][j][2];
	  } else {
	    SolnBlk.dUdt[i+1][j][k_residual] = NavierStokes2D_U_VACUUM;
	  }
	  break;
	case TIME_STEPPING_MULTISTAGE_OPTIMAL_SMOOTHING :
	  SolnBlk.dUdt[i+1][j][k_residual] = NavierStokes2D_U_VACUUM;
	  break;
	default:
	  SolnBlk.dUdt[i+1][j][k_residual] = NavierStokes2D_U_VACUUM;
	  break;
	};
      }

      if (j >= SolnBlk.JCl && j <= SolnBlk.JCu) {

	// Evaluate the cell interface i-direction INVISCID fluxes.
	if (i == SolnBlk.ICl-1 && 
	    (SolnBlk.Grid.BCtypeW[j] == BC_REFLECTION ||
	     SolnBlk.Grid.BCtypeW[j] == BC_BURNING_SURFACE ||
	     SolnBlk.Grid.BCtypeW[j] == BC_MASS_INJECTION ||
	     SolnBlk.Grid.BCtypeW[j] == BC_WALL_VISCOUS_HEATFLUX ||
	     SolnBlk.Grid.BCtypeW[j] == BC_WALL_VISCOUS_ISOTHERMAL ||
	     SolnBlk.Grid.BCtypeW[j] == BC_MOVING_WALL_HEATFLUX ||
	     SolnBlk.Grid.BCtypeW[j] == BC_MOVING_WALL_ISOTHERMAL ||
	     SolnBlk.Grid.BCtypeW[j] == BC_RINGLEB_FLOW ||
	     SolnBlk.Grid.BCtypeW[j] == BC_CHARACTERISTIC)) {

	  dX = SolnBlk.Grid.xfaceW(i+1,j) - SolnBlk.Grid.Cell[i+1][j].Xc;
	  Wr = SolnBlk.W[i+1][j] + (SolnBlk.phi[i+1][j]^SolnBlk.dWdx[i+1][j])*dX.x +
	                           (SolnBlk.phi[i+1][j]^SolnBlk.dWdy[i+1][j])*dX.y;

	  // WEST face of cell (i+1,j) is a normal boundary.
	  if (SolnBlk.Grid.BCtypeW[j] == BC_REFLECTION) {
	    // WEST face of cell (i+1,j) is a REFLECTION boundary.
	    Wl = Reflect(Wr,SolnBlk.Grid.nfaceW(i+1,j));
	  } else if (SolnBlk.Grid.BCtypeW[j] == BC_BURNING_SURFACE) {
	    // WEST face of cell (i+1,j) is a BURNING_SURFACE boundary.
	    Wl = BurningSurface(Wr,SolnBlk.Grid.nfaceW(i+1,j));
	  } else if (SolnBlk.Grid.BCtypeW[j] == BC_MASS_INJECTION) {
	    // WEST face of cell (i+1,j) is a MASS_INJECTION boundary.
	    //Wl = MassInjection(Wr,SolnBlk.Grid.nfaceW(i+1,j));
	    Wl = MassInjection2(Wr,SolnBlk.Grid.xfaceW(i+1,j),SolnBlk.Grid.nfaceW(i+1,j),SolnBlk.Twall);
	    flux_function = IP.i_Flux_Function;
	    IP.i_Flux_Function = FLUX_FUNCTION_GODUNOV_WRS;
	  } else if (SolnBlk.Grid.BCtypeW[j] == BC_WALL_VISCOUS_HEATFLUX) {
	    // WEST face of cell (i+1,j) is a WALL_VISCOUS_HEATFLUX boundary.
	    Wl = WallViscousHeatFlux(Wr,SolnBlk.Grid.nfaceW(i+1,j));
	  } else if (SolnBlk.Grid.BCtypeW[j] == BC_WALL_VISCOUS_ISOTHERMAL) {
	    // WEST face of cell (i+1,j) is a WALL_VISCOUS_ISOTHERMAL boundary.
	    Wl = WallViscousIsothermal(Wr,SolnBlk.Grid.nfaceW(i+1,j),SolnBlk.Twall);
	  } else if (SolnBlk.Grid.BCtypeW[j] == BC_MOVING_WALL) {
	    // WEST face of cell (i+1,j) is a MOVINGWALL_HEATFLUX boundary.
	    Wl = MovingWallHeatFlux(Wr,SolnBlk.Grid.nfaceW(i+1,j),SolnBlk.Vwall.x);
	  } else if (SolnBlk.Grid.BCtypeW[j] == BC_MOVING_WALL_ISOTHERMAL) {
	    // WEST face of cell (i+1,j) is a MOVINGWALL_ISOTHERMAL boundary.
	    Wl = MovingWallIsothermal(Wr,SolnBlk.Grid.nfaceW(i+1,j),SolnBlk.Vwall.x,SolnBlk.Twall);
	  } else if (SolnBlk.Grid.BCtypeW[j] == BC_RINGLEB_FLOW) {
	    // WEST face of cell (i+1,j) is a RINGLEB_FLOW boundary.
	    Wl = RinglebFlow(Wl,SolnBlk.Grid.xfaceW(i+1,j));
	  } else {
	    // WEST face of cell (i+1,j) is a CHARACTERISTIC boundary.
	    Wl = BC_Characteristic_Pressure(Wr,SolnBlk.WoW[j],SolnBlk.Grid.nfaceW(i+1,j));
	  }

	} else if (i == SolnBlk.ICu &&
		   (SolnBlk.Grid.BCtypeE[j] == BC_REFLECTION ||
		    SolnBlk.Grid.BCtypeE[j] == BC_BURNING_SURFACE ||
		    SolnBlk.Grid.BCtypeE[j] == BC_MASS_INJECTION ||
		    SolnBlk.Grid.BCtypeE[j] == BC_WALL_VISCOUS_HEATFLUX ||
		    SolnBlk.Grid.BCtypeE[j] == BC_WALL_VISCOUS_ISOTHERMAL ||
		    SolnBlk.Grid.BCtypeE[j] == BC_MOVING_WALL_HEATFLUX ||
		    SolnBlk.Grid.BCtypeE[j] == BC_MOVING_WALL_ISOTHERMAL ||
		    SolnBlk.Grid.BCtypeE[j] == BC_RINGLEB_FLOW ||
		    SolnBlk.Grid.BCtypeE[j] == BC_CHARACTERISTIC)) {

	  dX = SolnBlk.Grid.xfaceE(i,j) - SolnBlk.Grid.Cell[i][j].Xc;
	  Wl = SolnBlk.W[i][j] + (SolnBlk.phi[i][j]^SolnBlk.dWdx[i][j])*dX.x +
	                         (SolnBlk.phi[i][j]^SolnBlk.dWdy[i][j])*dX.y;

	  // EAST face of cell (i,j) is a normal boundary.
	  if (SolnBlk.Grid.BCtypeE[j] == BC_REFLECTION) {
	    // EAST face of cell (i,j) is a REFLECTION boundary.
	    Wr = Reflect(Wl,SolnBlk.Grid.nfaceE(i,j));
	  } else if (SolnBlk.Grid.BCtypeE[j] == BC_BURNING_SURFACE) {
	    // EAST face of cell (i,j) is a BURNING_SURFACE boundary.
	    Wr = BurningSurface(Wl,SolnBlk.Grid.nfaceE(i,j));
	  } else if (SolnBlk.Grid.BCtypeE[j] == BC_MASS_INJECTION) {
	    // EAST face of cell (i,j) is a MASS_INJECTION boundary.
	    //Wr = MassInjection(Wl,SolnBlk.Grid.nfaceE(i,j));
	    Wr = MassInjection2(Wl,SolnBlk.Grid.xfaceE(i,j),SolnBlk.Grid.nfaceE(i,j),SolnBlk.Twall);
	    flux_function = IP.i_Flux_Function;
	    IP.i_Flux_Function = FLUX_FUNCTION_GODUNOV_WRS;
	  } else if (SolnBlk.Grid.BCtypeE[j] == BC_WALL_VISCOUS_HEATFLUX) {
	    // EAST face of cell (i,j) is a WALL_VISCOUS_HEATFLUX boundary.
	    Wr = WallViscousHeatFlux(Wl,SolnBlk.Grid.nfaceE(i,j));
	  } else if (SolnBlk.Grid.BCtypeE[j] == BC_WALL_VISCOUS_ISOTHERMAL) {
	    // EAST face of cell (i,j) is a WALL_VISCOUS_ISOTHERMAL boundary.
	    Wr = WallViscousIsothermal(Wl,SolnBlk.Grid.nfaceE(i,j),SolnBlk.Twall);
	  } else if (SolnBlk.Grid.BCtypeE[j] == BC_MOVING_WALL_HEATFLUX) {
	    // EAST face of cell (i,j) is a MOVINGWALL_HEATFLUX boundary.
	    Wr = MovingWallHeatFlux(Wl,SolnBlk.Grid.nfaceE(i,j),SolnBlk.Vwall.x);
	  } else if (SolnBlk.Grid.BCtypeE[j] == BC_MOVING_WALL_ISOTHERMAL) {
	    // EAST face of cell (i,j) is a MOVINGWALL_ISOTHERMAL boundary.
	    Wr = MovingWallIsothermal(Wl,SolnBlk.Grid.nfaceE(i,j),SolnBlk.Vwall.x,SolnBlk.Twall);
	  } else if (SolnBlk.Grid.BCtypeE[j] == BC_RINGLEB_FLOW) {
	    // EAST face of cell (i,j) is a RINGLEB_FLOW boundary.
	    Wr = RinglebFlow(Wr,SolnBlk.Grid.xfaceE(i,j));
	  } else {
	    // EAST face of cell (i,j) is a CHARACTERISTIC boundary.
	    Wr = BC_Characteristic_Pressure(Wl,SolnBlk.WoE[j],SolnBlk.Grid.nfaceE(i,j));
	  }

	} else {

	  // EAST face is either a normal cell or possibly a FIXED, 
	  // NONE or EXTRAPOLATION boundary.
	  dX = SolnBlk.Grid.xfaceE(i  ,j) - SolnBlk.Grid.Cell[i  ][j].Xc;
	  Wl = SolnBlk.W[i  ][j] + (SolnBlk.phi[i  ][j]^SolnBlk.dWdx[i  ][j])*dX.x +
	                           (SolnBlk.phi[i  ][j]^SolnBlk.dWdy[i  ][j])*dX.y;
	  dX = SolnBlk.Grid.xfaceW(i+1,j) - SolnBlk.Grid.Cell[i+1][j].Xc;
	  Wr = SolnBlk.W[i+1][j] + (SolnBlk.phi[i+1][j]^SolnBlk.dWdx[i+1][j])*dX.x +
	                           (SolnBlk.phi[i+1][j]^SolnBlk.dWdy[i+1][j])*dX.y;

	}

	// Determine EAST face INVISCID flux.
	switch(IP.i_Flux_Function) {
	case FLUX_FUNCTION_GODUNOV :
	  Flux = FluxGodunov_n(Wl,Wr,SolnBlk.Grid.nfaceE(i,j));
	  break;
	case FLUX_FUNCTION_ROE :
	  Flux = FluxRoe_n(Wl,Wr,SolnBlk.Grid.nfaceE(i,j));
	  break;
	case FLUX_FUNCTION_RUSANOV :
	  Flux = FluxRusanov_n(Wl,Wr,SolnBlk.Grid.nfaceE(i,j));
	  break;
	case FLUX_FUNCTION_HLLE :
	  Flux = FluxHLLE_n(Wl,Wr,SolnBlk.Grid.nfaceE(i,j));
	  break;
	case FLUX_FUNCTION_HLLL :
	  Flux = FluxHLLL_n(Wl,Wr,SolnBlk.Grid.nfaceE(i,j));
	  break;
	case FLUX_FUNCTION_HLLC :
	  Flux = FluxHLLC_n(Wl,Wr,SolnBlk.Grid.nfaceE(i,j));
	  break;
	case FLUX_FUNCTION_VANLEER :
	  Flux = FluxVanLeer_n(Wl,Wr,SolnBlk.Grid.nfaceE(i,j));
	  break;
	case FLUX_FUNCTION_AUSM :
	  Flux = FluxAUSM_n(Wl,Wr,SolnBlk.Grid.nfaceE(i,j));
	  break;
	case FLUX_FUNCTION_AUSMplus :
	  Flux = FluxAUSMplus_n(Wl,Wr,SolnBlk.Grid.nfaceE(i,j));
	  break;
	case FLUX_FUNCTION_GODUNOV_WRS :
	  Flux = FluxGodunov_Wrs_n(Wl,Wr,SolnBlk.Grid.nfaceN(i,j));
	  break;
	default:
	  Flux = FluxRoe_n(Wl,Wr,SolnBlk.Grid.nfaceE(i,j));
	  break;
	};

	if (IP.i_Flux_Function == FLUX_FUNCTION_GODUNOV_WRS) {
	  IP.i_Flux_Function = flux_function;
	}

 	// Compute the cell centred stress tensor and heat flux vector if required.
 	if (SolnBlk.Flow_Type == FLOWTYPE_TURBULENT_RANS_K_OMEGA ||
	    (SolnBlk.Flow_Type && SolnBlk.Axisymmetric)) {
 	  SolnBlk.W[i][j].ComputeViscousTerms(SolnBlk.dWdx[i][j],
 					      SolnBlk.dWdy[i][j],
 					      SolnBlk.Grid.Cell[i][j].Xc,
 					      SolnBlk.Axisymmetric,
 					      OFF,
					      SolnBlk.Wall[i][j].ywall,
					      SolnBlk.Wall[i][j].yplus);
 	  SolnBlk.U[i][j].tau = SolnBlk.W[i][j].tau;
 	  SolnBlk.U[i][j].q = SolnBlk.W[i][j].q;
 	}

	// Evaluate the cell interface i-direction VISCOUS flux if necessary.
	if (SolnBlk.Flow_Type) {
	  // Determine the EAST face VISCOUS flux.
	  if (i == SolnBlk.ICl-1 && 
	      (SolnBlk.Grid.BCtypeW[j] == BC_WALL_VISCOUS_HEATFLUX ||
	       SolnBlk.Grid.BCtypeW[j] == BC_WALL_VISCOUS_ISOTHERMAL ||
	       SolnBlk.Grid.BCtypeW[j] == BC_MOVING_WALL_HEATFLUX ||
	       SolnBlk.Grid.BCtypeW[j] == BC_MOVING_WALL_ISOTHERMAL ||
	       SolnBlk.Grid.BCtypeW[j] == BC_BURNING_SURFACE ||
	       SolnBlk.Grid.BCtypeW[j] == BC_MASS_INJECTION)) {
	    // WEST face of cell (i+1,j) is a normal boundary.
	    Xr = SolnBlk.Grid.Cell[i+1][j].Xc; Wr = SolnBlk.W[i+1][j];
	    if (SolnBlk.Grid.BCtypeW[j] == BC_WALL_VISCOUS_HEATFLUX) {
	      // WEST face of cell (i+1,j) is a WALL_VISCOUS_HEATFLUX boundary.
	      viscous_bc_flag = DIAMONDPATH_RIGHT_TRIANGLE_HEATFLUX;
	      //Wu.rho = Wr.rho; Wu.v.x = ZERO; Wu.v.y = ZERO; Wu.p = Wr.p; Wu.k = Wr.k; Wu.omega = Wr.omega;Wu.ke = Wr.ke; Wu.ee = Wr.ee;
	      Wu = HALF*(Wr+SolnBlk.W[i+1][j+1]); Wu.v = Vector2D_ZERO;
	      Wd = HALF*(Wr+SolnBlk.W[i+1][j-1]); Wd.v = Vector2D_ZERO;
	    } else if (SolnBlk.Grid.BCtypeW[j] == BC_WALL_VISCOUS_ISOTHERMAL) {
	      // WEST face of cell (i+1,j) is a WALL_VISCOUS_ISOTHERMAL boundary.
	      viscous_bc_flag = DIAMONDPATH_RIGHT_TRIANGLE_ISOTHERMAL;
	      //Wu.rho = Wr.p/(Wr.R*SolnBlk.Twall); Wu.v.x = ZERO; Wu.v.y = ZERO; Wu.p = Wr.p; Wu.k = Wr.k; Wu.omega = Wr.omega;Wu.ke = Wr.ke; Wu.ee = Wr.ee;
	      Wu = HALF*(Wr+SolnBlk.W[i+1][j+1]); Wu.v = Vector2D_ZERO; Wu.rho = Wu.p/(Wr.R*SolnBlk.Twall);
	      Wd = HALF*(Wr+SolnBlk.W[i+1][j-1]); Wd.v = Vector2D_ZERO; Wd.rho = Wd.p/(Wr.R*SolnBlk.Twall);
	    } else if (SolnBlk.Grid.BCtypeW[j] == BC_MOVING_WALL) {
	      // WEST face of cell (i+1,j) is a MOVINGWALL_HEATFLUX boundary.
	      viscous_bc_flag = DIAMONDPATH_RIGHT_TRIANGLE_HEATFLUX;
	      //Wu.rho = Wr.rho; Wu.v.x = SolnBlk.Vwall.x; Wu.v.y = SolnBlk.Vwall.y; Wu.p = Wr.p; Wu.k = Wr.k; Wu.omega = Wr.omega;Wu.ke = Wr.ke; Wu.ee = Wr.ee;
	      Wu = HALF*(Wr+SolnBlk.W[i+1][j+1]); Wu.v = SolnBlk.Vwall;
	      Wd = HALF*(Wr+SolnBlk.W[i+1][j-1]); Wd.v = SolnBlk.Vwall;
	    } else if (SolnBlk.Grid.BCtypeW[j] == BC_MOVING_WALL_ISOTHERMAL) {
	      // WEST face of cell (i+1,j) is a MOVINGWALL_ISOTHERMAL boundary.
	      viscous_bc_flag = DIAMONDPATH_RIGHT_TRIANGLE_ISOTHERMAL;
	      //Wu.rho = Wr.p/(Wr.R*SolnBlk.Twall); Wu.v.x = SolnBlk.Vwall.x; Wu.v.y = SolnBlk.Vwall.y; Wu.p = Wr.p; Wu.k = Wr.k; Wu.omega = Wr.omega;Wu.ke = Wr.ke; Wu.ee = Wr.ee;
	      Wu = HALF*(Wr+SolnBlk.W[i+1][j+1]); Wu.v = SolnBlk.Vwall; Wu.rho = Wu.p/(Wr.R*SolnBlk.Twall);
	      Wd = HALF*(Wr+SolnBlk.W[i+1][j-1]); Wd.v = SolnBlk.Vwall; Wd.rho = Wd.p/(Wr.R*SolnBlk.Twall);
	    } else if (SolnBlk.Grid.BCtypeW[j] == BC_BURNING_SURFACE) {
	      // WEST face of cell (i+1,j) is a BURNING_SURFACE boundary.
	      viscous_bc_flag = DIAMONDPATH_RIGHT_TRIANGLE_ISOTHERMAL;
	      Wu = BurningSurface(Wr,SolnBlk.Grid.nfaceW(i+1,j));
	      Wd = Wu;
	    } else if (SolnBlk.Grid.BCtypeW[j] == BC_MASS_INJECTION) {
	      // WEST face of cell (i+1,j) is a MASS_INJECTION boundary.
	      viscous_bc_flag = DIAMONDPATH_RIGHT_TRIANGLE_ISOTHERMAL;
	      //Wu = MassInjection(Wr,SolnBlk.Grid.nfaceW(i+1,j));
	      Wu = MassInjection2(Wr,SolnBlk.Grid.xfaceW(i+1,j),SolnBlk.Grid.nfaceW(i+1,j),SolnBlk.Twall);
	      Wd = Wu;
	    }
	    switch(IP.i_Viscous_Reconstruction) {
	    case VISCOUS_RECONSTRUCTION_CARTESIAN :
	    case VISCOUS_RECONSTRUCTION_DIAMOND_PATH :
	      Xu = SolnBlk.Grid.Node[i+1][j+1].X;
	      Xd = SolnBlk.Grid.Node[i+1][j  ].X; Wd = Wu;
	      break;
	    case VISCOUS_RECONSTRUCTION_ARITHMETIC_AVERAGE :
	    case VISCOUS_RECONSTRUCTION_HYBRID :
	      Xu = SolnBlk.Grid.xfaceW(i+1,j);
	      Xl = Xu; Wl = Wu;
	      dWdxr = SolnBlk.dWdx[i+1][j]; dWdxl = dWdxr;
	      dWdyr = SolnBlk.dWdy[i+1][j]; dWdyl = dWdyr;
	      break;
	    };
	  } else if (i == SolnBlk.ICu &&
		     (SolnBlk.Grid.BCtypeE[j] == BC_WALL_VISCOUS_HEATFLUX ||
		      SolnBlk.Grid.BCtypeE[j] == BC_WALL_VISCOUS_ISOTHERMAL ||
		      SolnBlk.Grid.BCtypeE[j] == BC_MOVING_WALL_HEATFLUX ||
		      SolnBlk.Grid.BCtypeE[j] == BC_MOVING_WALL_ISOTHERMAL ||
		      SolnBlk.Grid.BCtypeE[j] == BC_BURNING_SURFACE ||
		      SolnBlk.Grid.BCtypeE[j] == BC_MASS_INJECTION)) {
	    // EAST face of cell (i,j) is a normal boundary.
	    Xl = SolnBlk.Grid.Cell[i][j].Xc; Wl = SolnBlk.W[i][j];
	    if (SolnBlk.Grid.BCtypeE[j] == BC_WALL_VISCOUS_HEATFLUX) {
	      // EAST face of cell (i,j) is a WALL_VISCOUS_HEATFLUX boundary.
	      viscous_bc_flag = DIAMONDPATH_LEFT_TRIANGLE_HEATFLUX;
	      //Wu.rho = Wl.rho; Wu.v.x = ZERO; Wu.v.y = ZERO; Wu.p = Wl.p; Wu.k = Wl.k; Wu.omega = Wl.omega;Wu.ke = Wl.ke; Wu.ee = Wl.ee;
	      Wu = HALF*(Wl+SolnBlk.W[i][j+1]); Wu.v = Vector2D_ZERO;
	      Wd = HALF*(Wl+SolnBlk.W[i][j-1]); Wd.v = Vector2D_ZERO;
	    } else if (SolnBlk.Grid.BCtypeE[j] == BC_WALL_VISCOUS_ISOTHERMAL) {
	      // EAST face of cell (i,j) is a WALL_VISCOUS_ISOTHERMAL boundary.
	      viscous_bc_flag = DIAMONDPATH_LEFT_TRIANGLE_ISOTHERMAL;
	      //Wu.rho = Wl.p/(Wl.R*SolnBlk.Twall); Wu.v.x = ZERO; Wu.v.y = ZERO; Wu.p = Wl.p; Wu.k = Wl.k; Wu.omega = Wl.omega;Wu.ke = Wl.ke; Wu.ee = Wl.ee;
	      Wu = HALF*(Wl+SolnBlk.W[i][j+1]); Wu.v = Vector2D_ZERO; Wu.rho = Wu.p/(Wl.R*SolnBlk.Twall);
	      Wd = HALF*(Wl+SolnBlk.W[i][j-1]); Wd.v = Vector2D_ZERO; Wd.rho = Wd.p/(Wl.R*SolnBlk.Twall);
	    } else if (SolnBlk.Grid.BCtypeE[j] == BC_MOVING_WALL_HEATFLUX) {
	      // EAST face of cell (i,j) is a MOVINGWALL_HEATFLUX boundary.
	      viscous_bc_flag = DIAMONDPATH_LEFT_TRIANGLE_HEATFLUX;
	      //Wu.rho = Wl.rho; Wu.v.x = SolnBlk.Vwall.x; Wu.v.y = SolnBlk.Vwall.y; Wu.p = Wl.p; Wu.k = Wl.k; Wu.omega = Wl.omega;Wu.ke = Wl.ke; Wu.ee = Wl.ee;
	      Wu = HALF*(Wl+SolnBlk.W[i][j+1]); Wu.v = SolnBlk.Vwall;
	      Wd = HALF*(Wl+SolnBlk.W[i][j-1]); Wd.v = SolnBlk.Vwall;
	    } else if (SolnBlk.Grid.BCtypeE[j] == BC_MOVING_WALL_ISOTHERMAL) {
	      // EAST face of cell (i,j) is a MOVINGWALL_ISOTHERMAL boundary.
	      viscous_bc_flag = DIAMONDPATH_LEFT_TRIANGLE_ISOTHERMAL;
	      //Wu.rho = Wl.p/(Wl.R*SolnBlk.Twall); Wu.v.x = SolnBlk.Vwall.x; Wu.v.y = SolnBlk.Vwall.y; Wu.p = Wl.p; Wu.k = Wl.k; Wu.omega = Wl.omega;Wu.ke = Wl.ke; Wu.ee = Wl.ee;
	      Wu = HALF*(Wl+SolnBlk.W[i][j+1]); Wu.v = SolnBlk.Vwall; Wu.rho = Wu.p/(Wl.R*SolnBlk.Twall);
	      Wd = HALF*(Wl+SolnBlk.W[i][j-1]); Wd.v = SolnBlk.Vwall; Wd.rho = Wd.p/(Wl.R*SolnBlk.Twall);
	    } else if (SolnBlk.Grid.BCtypeE[j] == BC_MOVING_WALL_ISOTHERMAL) {
	      // EAST face of cell (i,j) is a BURNING_SURFACE boundary.
	      viscous_bc_flag = DIAMONDPATH_LEFT_TRIANGLE_ISOTHERMAL;
	      Wu = BurningSurface(Wr,SolnBlk.Grid.nfaceE(i,j));
	      Wd = Wu;
	    } else if (SolnBlk.Grid.BCtypeE[j] == BC_MASS_INJECTION) {
	      // EAST face of cell (i,j) is a MASS_INJECTION boundary.
	      viscous_bc_flag = DIAMONDPATH_LEFT_TRIANGLE_ISOTHERMAL;
	      //Wu = MassInjection(Wr,SolnBlk.Grid.nfaceE(i,j));
	      Wu = MassInjection2(Wr,SolnBlk.Grid.xfaceE(i,j),SolnBlk.Grid.nfaceE(i,j),SolnBlk.Twall);
	      Wd = Wu;
	    }
	    switch(IP.i_Viscous_Reconstruction) {
	    case VISCOUS_RECONSTRUCTION_CARTESIAN :
	    case VISCOUS_RECONSTRUCTION_DIAMOND_PATH :
	      Xu = SolnBlk.Grid.Node[i+1][j+1].X;
	      Xd = SolnBlk.Grid.Node[i+1][j  ].X;
	      break;
	    case VISCOUS_RECONSTRUCTION_ARITHMETIC_AVERAGE :
	    case VISCOUS_RECONSTRUCTION_HYBRID :
	      Xu = SolnBlk.Grid.xfaceE(i,j);
	      Xr = Xu; Wr = Wu;
	      dWdxl = SolnBlk.dWdx[i][j]; dWdxr = dWdxl;
	      dWdyl = SolnBlk.dWdy[i][j]; dWdyr = dWdyl;
	      break;
	    };
	  } else {
	    // EAST face is either a normal cell or possibly a non-
	    // viscous boundary condition.
	    Xl = SolnBlk.Grid.Cell[i  ][j].Xc; Wl = SolnBlk.W[i  ][j];
	    Xr = SolnBlk.Grid.Cell[i+1][j].Xc; Wr = SolnBlk.W[i+1][j];
	    switch(IP.i_Viscous_Reconstruction) {
	    case VISCOUS_RECONSTRUCTION_CARTESIAN :
	    case VISCOUS_RECONSTRUCTION_DIAMOND_PATH :
	      viscous_bc_flag = DIAMONDPATH_QUADRILATERAL_RECONSTRUCTION;
	      Xu = SolnBlk.Grid.Node[i+1][j+1].X; Wu = SolnBlk.WnNE(i,j);
	      Xd = SolnBlk.Grid.Node[i+1][j  ].X; Wd = SolnBlk.WnSE(i,j);
	      break;
	    case VISCOUS_RECONSTRUCTION_ARITHMETIC_AVERAGE :
	    case VISCOUS_RECONSTRUCTION_HYBRID :
	      Xu = SolnBlk.Grid.xfaceE(i,j); Wu = HALF*(Wl + Wr);
	      dWdxl = SolnBlk.dWdx[i][j]; dWdxr = SolnBlk.dWdy[i+1][j];
	      dWdyl = SolnBlk.dWdy[i][j]; dWdyr = SolnBlk.dWdy[i+1][j];
	      break;
	    };
	  }
	  // Compute the EAST face viscous flux.
	  switch(IP.i_Viscous_Reconstruction) {
	  case VISCOUS_RECONSTRUCTION_CARTESIAN :
	  case VISCOUS_RECONSTRUCTION_DIAMOND_PATH :
	    Flux -= ViscousFluxDiamondPath_n(SolnBlk.Grid.xfaceE(i,j),
					     Xl,Wl,Xd,Wd,Xr,Wr,Xu,Wu,
					     SolnBlk.Grid.nfaceE(i,j),
					     SolnBlk.Axisymmetric,
					     viscous_bc_flag,
					     SolnBlk.Wall[i][j].ywall,
					     SolnBlk.Wall[i][j].yplus);
	    break;
	  case VISCOUS_RECONSTRUCTION_ARITHMETIC_AVERAGE :
	    Flux -= ViscousFlux_n(Xu,Wu,HALF*(dWdxl+dWdxr),
				  HALF*(dWdyl+dWdyr),
				  SolnBlk.Grid.nfaceE(i,j),
				  SolnBlk.Axisymmetric,OFF,
				  SolnBlk.Wall[i][j].ywall,
				  SolnBlk.Wall[i][j].yplus);
	    break;
	  case VISCOUS_RECONSTRUCTION_HYBRID :
 	    Flux -= ViscousFluxHybrid_n(Xu,Wu,Xl,Wl,dWdxl,dWdyl,
 					Xr,Wr,dWdxr,dWdyr,
 					SolnBlk.Grid.nfaceE(i,j),
 					SolnBlk.Axisymmetric,
 					SolnBlk.Wall[i][j].ywall,
 					SolnBlk.Wall[i][j].yplus);
	    break;
	  };
	}

	// Evaluate cell-averaged solution changes.
	SolnBlk.dUdt[i  ][j][k_residual] -= (IP.CFL_Number*SolnBlk.dt[i][j])*
                                            Flux*SolnBlk.Grid.lfaceE(i,j)/SolnBlk.Grid.Cell[i][j].A;
	SolnBlk.dUdt[i+1][j][k_residual] += (IP.CFL_Number*SolnBlk.dt[i+1][j])*
	                                    Flux*SolnBlk.Grid.lfaceW(i+1,j)/SolnBlk.Grid.Cell[i+1][j].A;

   	// Include axisymmetric source terms if required.
       	if (SolnBlk.Axisymmetric) {
   	  SolnBlk.dUdt[i][j][k_residual] += (IP.CFL_Number*SolnBlk.dt[i][j])*
 	                                    SolnBlk.W[i][j].Si(SolnBlk.Grid.Cell[i][j].Xc);
	  if (SolnBlk.Flow_Type)
  	    SolnBlk.dUdt[i][j][k_residual] += (IP.CFL_Number*SolnBlk.dt[i][j])*
  	                                      SolnBlk.W[i][j].Sv(SolnBlk.Grid.Cell[i][j].Xc,
								 SolnBlk.dWdy[i][j],
								 SolnBlk.Wall[i][j].ywall,
								 SolnBlk.Wall[i][j].yplus);
 	}

 	// Include turbulent production and destruction source term if required.
 	if (SolnBlk.Flow_Type == FLOWTYPE_TURBULENT_RANS_K_OMEGA) {
 	  SolnBlk.dUdt[i][j][k_residual] += (IP.CFL_Number*SolnBlk.dt[i][j])*
 	                                    SolnBlk.W[i][j].St(SolnBlk.Grid.Cell[i][j].Xc,
							       SolnBlk.W[i][j],
							       SolnBlk.dWdx[i][j],
							       SolnBlk.dWdy[i][j],
 							       SolnBlk.Axisymmetric,
							       SolnBlk.Wall[i][j].ywall,
							       SolnBlk.Wall[i][j].yplus);
	}

	// Save west and east face boundary flux.
 	if (i == SolnBlk.ICl-1) {
 	  SolnBlk.FluxW[j] = -Flux*SolnBlk.Grid.lfaceW(i+1,j);
 	} else if (i == SolnBlk.ICu) {
 	  SolnBlk.FluxE[j] =  Flux*SolnBlk.Grid.lfaceE(i,j);
 	}

      }
    }

    if (j > SolnBlk.JCl-1 && j < SolnBlk.JCu+1) {
      SolnBlk.dUdt[SolnBlk.ICl-1][j][k_residual] = NavierStokes2D_U_VACUUM;
      SolnBlk.dUdt[SolnBlk.ICu+1][j][k_residual] = NavierStokes2D_U_VACUUM;
    }

  }

  // Add j-direction (eta-direction) fluxes.
  for (int i = SolnBlk.ICl; i <= SolnBlk.ICu; i++) {
    for (int j = SolnBlk.JCl-1; j <= SolnBlk.JCu; j++) {

      // Evaluate the cell interface j-direction INVISCID fluxes.
      if (j == SolnBlk.JCl-1 && 
	  (SolnBlk.Grid.BCtypeS[i] == BC_REFLECTION ||
	   SolnBlk.Grid.BCtypeS[i] == BC_BURNING_SURFACE ||
	   SolnBlk.Grid.BCtypeS[i] == BC_MASS_INJECTION ||
	   SolnBlk.Grid.BCtypeS[i] == BC_WALL_VISCOUS_HEATFLUX ||
	   SolnBlk.Grid.BCtypeS[i] == BC_WALL_VISCOUS_ISOTHERMAL ||
	   SolnBlk.Grid.BCtypeS[i] == BC_MOVING_WALL_HEATFLUX ||
	   SolnBlk.Grid.BCtypeS[i] == BC_MOVING_WALL_ISOTHERMAL ||
	   SolnBlk.Grid.BCtypeS[i] == BC_RINGLEB_FLOW ||
	   SolnBlk.Grid.BCtypeS[i] == BC_CHARACTERISTIC)) {

	dX = SolnBlk.Grid.xfaceS(i,j+1) - SolnBlk.Grid.Cell[i][j+1].Xc;
	Wr = SolnBlk.W[i][j+1] + (SolnBlk.phi[i][j+1]^SolnBlk.dWdx[i][j+1])*dX.x +
   	                         (SolnBlk.phi[i][j+1]^SolnBlk.dWdy[i][j+1])*dX.y;

	// SOUTH face of cell (i,j+1) is a normal boundary.
	if (SolnBlk.Grid.BCtypeS[i] == BC_REFLECTION) {
	  // SOUTH face of cell (i,j+1) is a REFLECTION boundary.
	  Wl = Reflect(Wr,SolnBlk.Grid.nfaceS(i,j+1));
	} else if (SolnBlk.Grid.BCtypeS[i] == BC_BURNING_SURFACE) {
	  // SOUTH face of cell (i,j+1) is a BURNING_SURFACE boundary.
	  Wl = BurningSurface(Wr,SolnBlk.Grid.nfaceS(i,j+1));
	} else if (SolnBlk.Grid.BCtypeS[i] == BC_MASS_INJECTION) {
	  // SOUTH face of cell (i,j+1) is a MASS_INJECTION boundary.
	  //Wl = MassInjection(Wr,SolnBlk.Grid.nfaceS(i,j+1));
	  Wl = MassInjection2(Wr,SolnBlk.Grid.xfaceS(i,j+1),SolnBlk.Grid.nfaceS(i,j+1),SolnBlk.Twall);
	  flux_function = IP.i_Flux_Function;
	  IP.i_Flux_Function = FLUX_FUNCTION_GODUNOV_WRS;
	} else if (SolnBlk.Grid.BCtypeS[i] == BC_WALL_VISCOUS_HEATFLUX) {
	  // SOUTH face of cell (i,j+1) is a WALL_VISCOUS_HEATFLUX boundary.
	  Wl = WallViscousHeatFlux(Wr,SolnBlk.Grid.nfaceS(i,j+1));
	} else if (SolnBlk.Grid.BCtypeS[i] == BC_WALL_VISCOUS_ISOTHERMAL) {
	  // SOUTH face of cell (i,j+1) is a WALL_VISCOUS_ISOTHERMAL boundary.
	  Wl = WallViscousIsothermal(Wr,SolnBlk.Grid.nfaceS(i,j+1),SolnBlk.Twall);
	} else if (SolnBlk.Grid.BCtypeS[i] == BC_MOVING_WALL_HEATFLUX) {
	  // SOUTH face of cell (i,j+1) is a MOVINGWALL_HEATFLUX boundary.
	  Wl = MovingWallHeatFlux(Wr,SolnBlk.Grid.nfaceS(i,j+1),SolnBlk.Vwall.x);
	} else if (SolnBlk.Grid.BCtypeS[i] == BC_MOVING_WALL_ISOTHERMAL) {
	  // SOUTH face of cell (i,j+1) is a MOVINGWALL_ISOTHERMAL boundary.
	  Wl = MovingWallIsothermal(Wr,SolnBlk.Grid.nfaceS(i,j+1),SolnBlk.Vwall.x,SolnBlk.Twall);
	} else if (SolnBlk.Grid.BCtypeS[i] == BC_RINGLEB_FLOW) {
	  // SOUTH face of cell (i,j+1) is a RINGLEB_FLOW boundary.
	  Wl = RinglebFlow(Wl,SolnBlk.Grid.xfaceS(i,j+1));
	  Wl = BC_Characteristic_Pressure(Wr,Wl,SolnBlk.Grid.nfaceS(i,j+1));
	} else if (SolnBlk.Grid.BCtypeS[i] == BC_CHARACTERISTIC) {
	  // SOUTH face of cell (i,j+1) is a CHARACTERISTIC boundary.
	  Wl = BC_Characteristic_Pressure(Wr,SolnBlk.WoS[i],SolnBlk.Grid.nfaceS(i,j+1));
	}

      } else if (j == SolnBlk.JCu && 
		 (SolnBlk.Grid.BCtypeN[i] == BC_REFLECTION ||
		  SolnBlk.Grid.BCtypeN[i] == BC_BURNING_SURFACE ||
		  SolnBlk.Grid.BCtypeN[i] == BC_MASS_INJECTION ||
		  SolnBlk.Grid.BCtypeN[i] == BC_WALL_VISCOUS_HEATFLUX ||
		  SolnBlk.Grid.BCtypeN[i] == BC_WALL_VISCOUS_ISOTHERMAL ||
		  SolnBlk.Grid.BCtypeN[i] == BC_MOVING_WALL_HEATFLUX ||
		  SolnBlk.Grid.BCtypeN[i] == BC_MOVING_WALL_ISOTHERMAL ||
		  SolnBlk.Grid.BCtypeN[i] == BC_RINGLEB_FLOW ||
		  SolnBlk.Grid.BCtypeN[i] == BC_CHARACTERISTIC)) {

	dX = SolnBlk.Grid.xfaceN(i,j) - SolnBlk.Grid.Cell[i][j].Xc;
	Wl = SolnBlk.W[i][j] + (SolnBlk.phi[i][j]^SolnBlk.dWdx[i][j])*dX.x +
   	                       (SolnBlk.phi[i][j]^SolnBlk.dWdy[i][j])*dX.y;

	// NORTH face of cell (i,j) is a normal boundary.
	if (SolnBlk.Grid.BCtypeN[i] == BC_REFLECTION) {
	  // NORTH face of cell (i,j) is a REFLECTION boundary.
	  Wr = Reflect(Wl,SolnBlk.Grid.nfaceN(i,j));
	} else if (SolnBlk.Grid.BCtypeN[i] == BC_BURNING_SURFACE) {
	  // NORTH face of cell (i,j) is a BURNING_SURFACE boundary.
	  Wr = BurningSurface(Wl,SolnBlk.Grid.nfaceN(i,j));
	} else if (SolnBlk.Grid.BCtypeN[i] == BC_MASS_INJECTION) {
	  // NORTH face of cell (i,j) is a MASS_INJECTION boundary.
	  //Wr = MassInjection(Wl,SolnBlk.Grid.nfaceN(i,j));
	  Wr = MassInjection2(Wl,SolnBlk.Grid.xfaceN(i,j),SolnBlk.Grid.nfaceN(i,j),SolnBlk.Twall);
	  flux_function = IP.i_Flux_Function;
	  IP.i_Flux_Function = FLUX_FUNCTION_GODUNOV_WRS;
	} else if (SolnBlk.Grid.BCtypeN[i] == BC_WALL_VISCOUS_HEATFLUX) {
	  // NORTH face of cell (i,j) is a WALL_VISCOUS_HEATFLUX boundary.
	  Wr = WallViscousHeatFlux(Wl,SolnBlk.Grid.nfaceN(i,j));
	} else if (SolnBlk.Grid.BCtypeN[i] == BC_WALL_VISCOUS_ISOTHERMAL) {
	  // NORTH face of cell (i,j) is a WALL_VISCOUS_ISOTHERMAL boundary.
	  Wr = WallViscousIsothermal(Wl,SolnBlk.Grid.nfaceN(i,j),SolnBlk.Twall);
	} else if (SolnBlk.Grid.BCtypeN[i] == BC_MOVING_WALL_HEATFLUX) {
	  // NORTH face of cell (i,j) is a MOVINGWALL_HEATFLUX boundary.
	  Wr = MovingWallHeatFlux(Wl,SolnBlk.Grid.nfaceN(i,j),SolnBlk.Vwall.x);
	} else if (SolnBlk.Grid.BCtypeN[i] == BC_MOVING_WALL_ISOTHERMAL) {
	  // NORTH face of cell (i,j) is a MOVINGWALL_ISOTHERMAL boundary.
	  Wr = MovingWallIsothermal(Wl,SolnBlk.Grid.nfaceN(i,j),SolnBlk.Vwall.x,SolnBlk.Twall);
	} else if (SolnBlk.Grid.BCtypeN[i] == BC_RINGLEB_FLOW) {
	  // NORTH face of cell (i,j) is a RINGLEB_FLOW boundary.
	  Wr = RinglebFlow(Wr,SolnBlk.Grid.xfaceN(i,j));
	  Wr = BC_Characteristic_Pressure(Wl,Wr,SolnBlk.Grid.nfaceN(i,j));
	} else {
	  // NORTH face of cell (i,j) is a CHARACTERISTIC boundary.
	  Wr = BC_Characteristic_Pressure(Wl,SolnBlk.WoN[i],SolnBlk.Grid.nfaceN(i,j));
	}

      } else {
	
	// NORTH face is either a normal cell or possibly a FIXED, 
	// NONE or EXTRAPOLATION boundary.
	dX = SolnBlk.Grid.xfaceN(i,j  ) - SolnBlk.Grid.Cell[i][j  ].Xc;
	Wl = SolnBlk.W[i][j  ] + (SolnBlk.phi[i][j  ]^SolnBlk.dWdx[i][j  ])*dX.x +
                                 (SolnBlk.phi[i][j  ]^SolnBlk.dWdy[i][j  ])*dX.y;
	dX = SolnBlk.Grid.xfaceS(i,j+1) - SolnBlk.Grid.Cell[i][j+1].Xc;
	Wr = SolnBlk.W[i][j+1] + (SolnBlk.phi[i][j+1]^SolnBlk.dWdx[i][j+1])*dX.x +
                                 (SolnBlk.phi[i][j+1]^SolnBlk.dWdy[i][j+1])*dX.y;

      }

      // Determine NORTH face inviscid flux.
      switch(IP.i_Flux_Function) {
      case FLUX_FUNCTION_GODUNOV :
	Flux = FluxGodunov_n(Wl,Wr,SolnBlk.Grid.nfaceN(i,j));
	break;
      case FLUX_FUNCTION_ROE :
	Flux = FluxRoe_n(Wl,Wr,SolnBlk.Grid.nfaceN(i,j));
	break;
      case FLUX_FUNCTION_RUSANOV :
	Flux = FluxRusanov_n(Wl,Wr,SolnBlk.Grid.nfaceN(i,j));
	break;
      case FLUX_FUNCTION_HLLE :
	Flux = FluxHLLE_n(Wl,Wr,SolnBlk.Grid.nfaceN(i,j));
	break;
      case FLUX_FUNCTION_HLLL :
	Flux = FluxHLLL_n(Wl,Wr,SolnBlk.Grid.nfaceN(i,j));
	break;
      case FLUX_FUNCTION_HLLC :
	Flux = FluxHLLC_n(Wl,Wr,SolnBlk.Grid.nfaceN(i,j));
	break;
      case FLUX_FUNCTION_VANLEER :
	Flux = FluxVanLeer_n(Wl,Wr,SolnBlk.Grid.nfaceN(i,j));
	break;
      case FLUX_FUNCTION_AUSM :
	Flux = FluxAUSM_n(Wl,Wr,SolnBlk.Grid.nfaceN(i,j));
	break;
      case FLUX_FUNCTION_AUSMplus :
	Flux = FluxAUSMplus_n(Wl,Wr,SolnBlk.Grid.nfaceN(i,j));
	break;
      case FLUX_FUNCTION_GODUNOV_WRS :
	Flux = FluxGodunov_Wrs_n(Wl,Wr,SolnBlk.Grid.nfaceN(i,j));
	break;
      default:
	Flux = FluxRoe_n(Wl,Wr,SolnBlk.Grid.nfaceN(i,j));
	break;
      };

      if (IP.i_Flux_Function == FLUX_FUNCTION_GODUNOV_WRS) {
	IP.i_Flux_Function = flux_function;
      }

      // Evaluate the cell interface j-direction VISCOUS flux if necessary.
      if (SolnBlk.Flow_Type) {
	// Determine the NORTH face VISCOUS flux.
	if (j == SolnBlk.JCl-1 && 
	    (SolnBlk.Grid.BCtypeS[i] == BC_WALL_VISCOUS_HEATFLUX ||
	     SolnBlk.Grid.BCtypeS[i] == BC_WALL_VISCOUS_ISOTHERMAL ||
	     SolnBlk.Grid.BCtypeS[i] == BC_MOVING_WALL_HEATFLUX ||
	     SolnBlk.Grid.BCtypeS[i] == BC_MOVING_WALL_ISOTHERMAL ||
	     SolnBlk.Grid.BCtypeS[i] == BC_BURNING_SURFACE ||
	     SolnBlk.Grid.BCtypeS[i] == BC_MASS_INJECTION)) {
	  // SOUTH face of cell (i,j+1) is a normal boundary.
	  Xr = SolnBlk.Grid.Cell[i][j+1].Xc; Wr = SolnBlk.W[i][j+1];
	  if (SolnBlk.Grid.BCtypeS[i] == BC_WALL_VISCOUS_HEATFLUX) {
	    // SOUTH face of cell (i,j+1) is a WALL_VISCOUS_HEATFLUX boundary.
	    viscous_bc_flag = DIAMONDPATH_RIGHT_TRIANGLE_HEATFLUX;
	    //Wu.rho = Wr.rho; Wu.v.x = ZERO; Wu.v.y = ZERO; Wu.p = Wr.p; Wu.k = Wr.k; Wu.omega = Wr.omega;Wu.ke = Wr.ke; Wu.ee = Wr.ee;
	    Wu = HALF*(Wr+SolnBlk.W[i-1][j+1]); Wu.v = Vector2D_ZERO;
	    Wd = HALF*(Wr+SolnBlk.W[i+1][j+1]); Wd.v = Vector2D_ZERO;
	  } else if (SolnBlk.Grid.BCtypeS[i] == BC_WALL_VISCOUS_ISOTHERMAL) {
	    // SOUTH face of cell (i,j+1) is a WALL_VISCOUS_ISOTHERMAL boundary.
	    viscous_bc_flag = DIAMONDPATH_RIGHT_TRIANGLE_ISOTHERMAL;
	    //Wu.rho = Wr.p/(Wr.R*SolnBlk.Twall); Wu.v.x = ZERO; Wu.v.y = ZERO; Wu.p = Wr.p; Wu.k = Wr.k; Wu.omega = Wr.omega;Wu.ke = Wr.ke; Wu.ee = Wr.ee;
	    Wu = HALF*(Wr+SolnBlk.W[i-1][j+1]); Wu.v = Vector2D_ZERO; Wu.rho = Wu.p/(Wr.R*SolnBlk.Twall);
	    Wd = HALF*(Wr+SolnBlk.W[i+1][j+1]); Wd.v = Vector2D_ZERO; Wd.rho = Wd.p/(Wr.R*SolnBlk.Twall);
	  } else if (SolnBlk.Grid.BCtypeS[i] == BC_MOVING_WALL_HEATFLUX) {
	    // SOUTH face of cell (i,j+1) is a MOVINGWALL_HEATFLUX boundary.
	    viscous_bc_flag = DIAMONDPATH_RIGHT_TRIANGLE_HEATFLUX;
	    //Wu.rho = Wr.rho; Wu.v.x = SolnBlk.Vwall.x; Wu.v.y = SolnBlk.Vwall.y; Wu.p = Wr.p; Wu.k = Wr.k; Wu.omega = Wr.omega;Wu.ke = Wr.ke; Wu.ee = Wr.ee;
	    Wu = HALF*(Wr+SolnBlk.W[i-1][j+1]); Wu.v = SolnBlk.Vwall;
	    Wd = HALF*(Wr+SolnBlk.W[i+1][j+1]); Wd.v = SolnBlk.Vwall;
	  } else if (SolnBlk.Grid.BCtypeS[i] == BC_MOVING_WALL_ISOTHERMAL) {
	    // SOUTH face of cell (i,j+1) is a MOVINGWALL_ISOTHERMAL boundary.
	    viscous_bc_flag = DIAMONDPATH_RIGHT_TRIANGLE_ISOTHERMAL;
	    //Wu.rho = Wr.p/(Wr.R*SolnBlk.Twall); Wu.v.x = SolnBlk.Vwall.x; Wu.v.y = SolnBlk.Vwall.y; Wu.p = Wr.p; Wu.k = Wr.k; Wu.omega = Wr.omega;Wu.ke = Wr.ke; Wu.ee = Wr.ee;
	    Wu = HALF*(Wr+SolnBlk.W[i-1][j+1]); Wu.v = SolnBlk.Vwall; Wu.rho = Wu.p/(Wr.R*SolnBlk.Twall);
	    Wd = HALF*(Wr+SolnBlk.W[i+1][j+1]); Wd.v = SolnBlk.Vwall; Wd.rho = Wd.p/(Wr.R*SolnBlk.Twall);
	  } else if (SolnBlk.Grid.BCtypeS[i] == BC_BURNING_SURFACE) {
	    // SOUTH face of cell (i,j+1) is a BURNING_SURFACE boundary.
	    viscous_bc_flag = DIAMONDPATH_RIGHT_TRIANGLE_ISOTHERMAL;
	    Wu = BurningSurface(Wr,SolnBlk.Grid.nfaceS(i,j+1));
	    Wd = Wu;
	  } else if (SolnBlk.Grid.BCtypeS[i] == BC_MASS_INJECTION) {
	    // SOUTH face of cell (i,j+1) is a MASS_INJECTION boundary.
	    viscous_bc_flag = DIAMONDPATH_RIGHT_TRIANGLE_ISOTHERMAL;
	    //Wu = MassInjection(Wr,SolnBlk.Grid.nfaceS(i,j+1));
	    Wu = MassInjection2(Wr,SolnBlk.Grid.xfaceS(i,j+1),SolnBlk.Grid.nfaceS(i,j+1),SolnBlk.Twall);
	    Wd = Wu;
	  }
	  switch(IP.i_Viscous_Reconstruction) {
	  case VISCOUS_RECONSTRUCTION_CARTESIAN :
	  case VISCOUS_RECONSTRUCTION_DIAMOND_PATH :
	    Xu = SolnBlk.Grid.Node[i  ][j+1].X;
	    Xd = SolnBlk.Grid.Node[i+1][j+1].X;
	    break;
	  case VISCOUS_RECONSTRUCTION_ARITHMETIC_AVERAGE :
	  case VISCOUS_RECONSTRUCTION_HYBRID :
	    Xu = SolnBlk.Grid.xfaceS(i,j+1);
	    Xl = Xu; Wl = Wu;
	    dWdxr = SolnBlk.dWdx[i][j+1]; dWdxl = dWdxr;
	    dWdyr = SolnBlk.dWdy[i][j+1]; dWdyl = dWdyr;
	    break;
	  };
	} else if (j == SolnBlk.JCu && 
		   (SolnBlk.Grid.BCtypeN[i] == BC_WALL_VISCOUS_HEATFLUX ||
		    SolnBlk.Grid.BCtypeN[i] == BC_WALL_VISCOUS_ISOTHERMAL ||
		    SolnBlk.Grid.BCtypeN[i] == BC_MOVING_WALL_HEATFLUX ||
		    SolnBlk.Grid.BCtypeN[i] == BC_MOVING_WALL_ISOTHERMAL ||
		    SolnBlk.Grid.BCtypeN[i] == BC_BURNING_SURFACE ||
		    SolnBlk.Grid.BCtypeN[i] == BC_MASS_INJECTION)) {
	  // NORTH face of cell (i,j) is a normal boundary.
	  Xl = SolnBlk.Grid.Cell[i][j].Xc; Wl = SolnBlk.W[i][j];
	  if (SolnBlk.Grid.BCtypeN[i] == BC_WALL_VISCOUS_HEATFLUX) {
	    // NORTH face of cell (i,j) is a WALL_VISCOUS_HEATFLUX boundary.
	    viscous_bc_flag = DIAMONDPATH_LEFT_TRIANGLE_HEATFLUX;
	    //Wu.rho = Wl.rho; Wu.v.x = ZERO; Wu.v.y = ZERO; Wu.p = Wl.p; Wu.k = Wl.k; Wu.omega = Wl.omega;Wu.ke = Wl.ke; Wu.ee = Wl.ee;
	    Wu = HALF*(Wl+SolnBlk.W[i-1][j]); Wu.v = Vector2D_ZERO;
	    Wd = HALF*(Wl+SolnBlk.W[i+1][j]); Wd.v = Vector2D_ZERO;
	  } else if (SolnBlk.Grid.BCtypeN[i] == BC_WALL_VISCOUS_ISOTHERMAL) {
	    // NORTH face of cell (i,j) is a WALL_VISCOUS_ISOTHERMAL boundary.
	    viscous_bc_flag = DIAMONDPATH_LEFT_TRIANGLE_ISOTHERMAL;
	    //Wu.rho = Wl.p/(Wl.R*SolnBlk.Twall); Wu.v.x = ZERO; Wu.v.y = ZERO; Wu.p = Wl.p; Wu.k = Wl.k; Wu.omega = Wl.omega;Wu.ke = Wl.ke; Wu.ee = Wl.ee;
	    Wu = HALF*(Wl+SolnBlk.W[i-1][j]); Wu.v = Vector2D_ZERO; Wu.rho = Wu.p/(Wl.R*SolnBlk.Twall);
	    Wd = HALF*(Wl+SolnBlk.W[i+1][j]); Wd.v = Vector2D_ZERO; Wd.rho = Wd.p/(Wl.R*SolnBlk.Twall);
	  } else if (SolnBlk.Grid.BCtypeN[i] == BC_MOVING_WALL_HEATFLUX) {
	    // NORTH face of cell (i,j) is a MOVINGWALL_HEATFLUX boundary.
	    viscous_bc_flag = DIAMONDPATH_LEFT_TRIANGLE_HEATFLUX;
	    //Wu.rho = Wl.rho; Wu.v.x = SolnBlk.Vwall.x; Wu.v.y = SolnBlk.Vwall.y; Wu.p = Wl.p; Wu.k = Wl.k; Wu.omega = Wl.omega;Wu.ke = Wl.ke; Wu.ee = Wl.ee;
	    Wu = HALF*(Wl+SolnBlk.W[i-1][j]); Wu.v = SolnBlk.Vwall;
	    Wd = HALF*(Wl+SolnBlk.W[i+1][j]); Wd.v = SolnBlk.Vwall;
	  } else if (SolnBlk.Grid.BCtypeN[i] == BC_MOVING_WALL_ISOTHERMAL) {
	    // NORTH face of cell (i,j) is a MOVINGWALL_ISOTHERMAL boundary.
	    viscous_bc_flag = DIAMONDPATH_LEFT_TRIANGLE_ISOTHERMAL;
	    //Wu.rho = Wl.p/(Wl.R*SolnBlk.Twall); Wu.v.x = SolnBlk.Vwall.x; Wu.v.y = SolnBlk.Vwall.y; Wu.p = Wl.p; Wu.k = Wl.k; Wu.omega = Wl.omega;Wu.ke = Wl.ke; Wu.ee = Wl.ee;
	    Wu = HALF*(Wl+SolnBlk.W[i-1][j]); Wu.v = SolnBlk.Vwall; Wu.rho = Wu.p/(Wl.R*SolnBlk.Twall);
	    Wd = HALF*(Wl+SolnBlk.W[i+1][j]); Wd.v = SolnBlk.Vwall; Wd.rho = Wd.p/(Wl.R*SolnBlk.Twall);
	  } else if (SolnBlk.Grid.BCtypeN[i] == BC_BURNING_SURFACE) {
	    // NORTH face of cell (i,j) is a BURNING_SURFACE boundary.
	    viscous_bc_flag = DIAMONDPATH_LEFT_TRIANGLE_ISOTHERMAL;
	    Wu = BurningSurface(Wl,SolnBlk.Grid.nfaceN(i,j));
	    Wd = Wu;
	  } else if (SolnBlk.Grid.BCtypeN[i] == BC_MASS_INJECTION) {
	    // NORTH face of cell (i,j) is a MASS_INJECTION boundary.
	    viscous_bc_flag = DIAMONDPATH_LEFT_TRIANGLE_ISOTHERMAL;
	    //Wu = MassInjection(Wl,SolnBlk.Grid.nfaceN(i,j));
	    Wu = MassInjection2(Wl,SolnBlk.Grid.xfaceN(i,j),SolnBlk.Grid.nfaceN(i,j),SolnBlk.Twall);
	    Wd = Wu;
	  }
	  switch(IP.i_Viscous_Reconstruction) {
	  case VISCOUS_RECONSTRUCTION_CARTESIAN :
	  case VISCOUS_RECONSTRUCTION_DIAMOND_PATH :
	    Xu = SolnBlk.Grid.Node[i  ][j+1].X;
	    Xd = SolnBlk.Grid.Node[i+1][j+1].X;
	    break;
	  case VISCOUS_RECONSTRUCTION_ARITHMETIC_AVERAGE :
	  case VISCOUS_RECONSTRUCTION_HYBRID :
	    Xu = SolnBlk.Grid.xfaceN(i,j);
	    Xr = Xu; Wr = Wu;
	    dWdxl = SolnBlk.dWdx[i][j]; dWdxr = dWdxl;
	    dWdyl = SolnBlk.dWdy[i][j]; dWdyr = dWdyl;
	    break;
	  };
	} else {
	  // NORTH face is either a normal cell or possibly a non-viscous
	  // boundary condition.
	  Xl = SolnBlk.Grid.Cell[i][j  ].Xc; Wl = SolnBlk.W[i][j  ];
	  Xr = SolnBlk.Grid.Cell[i][j+1].Xc; Wr = SolnBlk.W[i][j+1];
	  switch(IP.i_Viscous_Reconstruction) {
	  case VISCOUS_RECONSTRUCTION_CARTESIAN :
	  case VISCOUS_RECONSTRUCTION_DIAMOND_PATH :
	    viscous_bc_flag = DIAMONDPATH_QUADRILATERAL_RECONSTRUCTION;
	    Xu = SolnBlk.Grid.Node[i  ][j+1].X; Wu = SolnBlk.WnNW(i,j);
	    Xd = SolnBlk.Grid.Node[i+1][j+1].X; Wd = SolnBlk.WnNE(i,j);
	    break;
	  case VISCOUS_RECONSTRUCTION_ARITHMETIC_AVERAGE :
	  case VISCOUS_RECONSTRUCTION_HYBRID :
	    Xu = SolnBlk.Grid.xfaceN(i,j); Wu = HALF*(Wl + Wr);
	    dWdxl = SolnBlk.dWdx[i][j]; dWdxr = SolnBlk.dWdy[i][j+1];
	    dWdyl = SolnBlk.dWdy[i][j]; dWdyr = SolnBlk.dWdy[i][j+1];
	    break;
	  };
	}
	// Compute the NORTH face viscous flux.
	switch(IP.i_Viscous_Reconstruction) {
	case VISCOUS_RECONSTRUCTION_CARTESIAN :
	case VISCOUS_RECONSTRUCTION_DIAMOND_PATH :
	  Flux -= ViscousFluxDiamondPath_n(SolnBlk.Grid.xfaceN(i,j),
					   Xl,Wl,Xd,Wd,Xr,Wr,Xu,Wu,
					   SolnBlk.Grid.nfaceN(i,j),
					   SolnBlk.Axisymmetric,
					   viscous_bc_flag,
					   SolnBlk.Wall[i][j].ywall,
					   SolnBlk.Wall[i][j].yplus);
	  break;
	case VISCOUS_RECONSTRUCTION_ARITHMETIC_AVERAGE :
	  Flux -= ViscousFlux_n(Xu,Wu,HALF*(dWdxl+dWdxr),
				HALF*(dWdyl+dWdyr),
				SolnBlk.Grid.nfaceN(i,j),
				SolnBlk.Axisymmetric,OFF,
				SolnBlk.Wall[i][j].ywall,
				SolnBlk.Wall[i][j].yplus);
	  break;
	case VISCOUS_RECONSTRUCTION_HYBRID :
	  Flux -= ViscousFluxHybrid_n(Xu,Wu,Xl,Wl,dWdxl,dWdyl,
				      Xr,Wr,dWdxr,dWdyr,
				      SolnBlk.Grid.nfaceN(i,j),
				      SolnBlk.Axisymmetric,
				      SolnBlk.Wall[i][j].ywall,
 				      SolnBlk.Wall[i][j].yplus);
	  break;
	};
      }

      // Evaluate cell-averaged solution changes.
      SolnBlk.dUdt[i][j  ][k_residual] -= (IP.CFL_Number*SolnBlk.dt[i][j  ])*
                                          Flux*SolnBlk.Grid.lfaceN(i,j)/SolnBlk.Grid.Cell[i][j].A;
      SolnBlk.dUdt[i][j+1][k_residual] += (IP.CFL_Number*SolnBlk.dt[i][j+1])*
                                          Flux*SolnBlk.Grid.lfaceS(i,j+1)/SolnBlk.Grid.Cell[i][j+1].A;

      // Save south and north face boundary flux.
      if (j == SolnBlk.JCl-1) {
 	SolnBlk.FluxS[i] = -Flux*SolnBlk.Grid.lfaceS(i,j+1);
      } else if (j == SolnBlk.JCu) {
	SolnBlk.FluxN[i] = Flux*SolnBlk.Grid.lfaceN(i,j);
      }

    }

    SolnBlk.dUdt[i][SolnBlk.JCl-1][k_residual] = NavierStokes2D_U_VACUUM;
    SolnBlk.dUdt[i][SolnBlk.JCu+1][k_residual] = NavierStokes2D_U_VACUUM;

  }

  // Zero the residuals for the turbulence variables according to
  // the turbulence boundary condition.
  error_flag = Turbulence_Zero_Residual(SolnBlk,i_stage,IP);
  if (error_flag) return error_flag;

  // Residual for the stage successfully calculated.
  return 0;

}

/**********************************************************************
 * Routine: Update_Solution_Multistage_Explicit                       *
 *                                                                    *
 * This routine updates solution states of the given solution block   *
 * for a variety of multi-stage explicit time integration schemes.    *
 *                                                                    *
 **********************************************************************/
int Update_Solution_Multistage_Explicit(NavierStokes2D_Quad_Block &SolnBlk,
					const int i_stage,
					NavierStokes2D_Input_Parameters &IP) {

  int k_residual;
  double omega, residual_reduction_factor;

  // Memory for linear system solver.
  DenseMatrix dRdU(NUM_VAR_NAVIERSTOKES2D,NUM_VAR_NAVIERSTOKES2D);
  DenseSystemLinEqs LinSys;

  // Allocate memory for linear system solver.
  LinSys.allocate(NUM_VAR_NAVIERSTOKES2D);

  // Perform update of solution variables for stage i_stage of N_Stage
  // scheme.

  // Evaluate the time step fraction and residual storage location for
  // the stage.
  switch(IP.i_Time_Integration) {
  case TIME_STEPPING_EXPLICIT_EULER :
    omega = Runge_Kutta(i_stage,IP.N_Stage);
    k_residual = 0;
    break;
  case TIME_STEPPING_EXPLICIT_PREDICTOR_CORRECTOR :
    omega = Runge_Kutta(i_stage,IP.N_Stage);
    k_residual = 0;
    break;
  case TIME_STEPPING_EXPLICIT_RUNGE_KUTTA :
    omega = Runge_Kutta(i_stage,IP.N_Stage);
    k_residual = 0;
    if (IP.N_Stage == 4) {
      if (i_stage == 4) {
	k_residual = 0;
      } else {
	k_residual = i_stage - 1;
      }
    }
    break;
  case TIME_STEPPING_MULTISTAGE_OPTIMAL_SMOOTHING :
    omega = MultiStage_Optimally_Smoothing(i_stage,IP.N_Stage,IP.i_Limiter);
    k_residual = 0;
    break;
  default:
    omega = Runge_Kutta(i_stage,IP.N_Stage);
    k_residual = 0;
    break;
  };

  // Update solution variables for this stage.
  for (int j = SolnBlk.JCl; j <= SolnBlk.JCu; j++) {
    for (int i = SolnBlk.ICl; i <= SolnBlk.ICu; i++) {

      // Explicit update.
      if (IP.Local_Time_Stepping == GLOBAL_TIME_STEPPING ||
	  IP.Local_Time_Stepping == SCALAR_LOCAL_TIME_STEPPING) {
	SolnBlk.U[i][j] = SolnBlk.Uo[i][j] + omega*SolnBlk.dUdt[i][j][k_residual];
      }

      // Turbulence parameter residual reductions.
      if (IP.Local_Time_Stepping == SCALAR_LOCAL_TIME_STEPPING &&
	  SolnBlk.Flow_Type == FLOWTYPE_TURBULENT_RANS_K_OMEGA) {
	if (SolnBlk.Uo[i][j].dk > ZERO) {
	  if (fabs(SolnBlk.U[i][j].dk-SolnBlk.Uo[i][j].dk)/SolnBlk.Uo[i][j].dk > 0.25) {
	    SolnBlk.U[i][j].dk = SolnBlk.Uo[i][j].dk + 0.249*omega*SolnBlk.dUdt[i][j][k_residual].dk;
	  }
	}
	if (SolnBlk.Uo[i][j].domega > ZERO) {
	  if (fabs(SolnBlk.U[i][j].domega-SolnBlk.Uo[i][j].domega)/SolnBlk.Uo[i][j].domega > 0.25) {
	    SolnBlk.U[i][j].domega = SolnBlk.Uo[i][j].domega + 0.249*omega*SolnBlk.dUdt[i][j][k_residual].domega;
	  }
	}
	if (SolnBlk.Variable_Prandtl == ON) {
	  if (SolnBlk.Uo[i][j].dke > ZERO) {
	    if (fabs(SolnBlk.U[i][j].dke-SolnBlk.Uo[i][j].dke)/SolnBlk.Uo[i][j].dke > 0.25) {
	      SolnBlk.U[i][j].dke = SolnBlk.Uo[i][j].dke + 0.249*omega*SolnBlk.dUdt[i][j][k_residual].dke;
	    }
	  }
	  if (SolnBlk.Uo[i][j].dee > ZERO) {
	    if (fabs(SolnBlk.U[i][j].dee-SolnBlk.Uo[i][j].dee)/SolnBlk.Uo[i][j].dee > 0.25) {
	      SolnBlk.U[i][j].dee = SolnBlk.Uo[i][j].dee + 0.249*omega*SolnBlk.dUdt[i][j][k_residual].dee;
	    }
	  }
	}
      }
      // Point implicit formulation: set up system of equations 
      // and include source Jacobian in the LHS matrix.
      if (IP.Local_Time_Stepping == SEMI_IMPLICIT_LOCAL_TIME_STEPPING) {

	// Evaluate the source Jacobian.
	dRdU.zero();

	// Calculate source Jacobian for axisymmetric terms if required.
	if (SolnBlk.Axisymmetric) {
	  SolnBlk.W[i][j].dSidU(dRdU,SolnBlk.Grid.Cell[i][j].Xc);
	  if (SolnBlk.Flow_Type)
	    SolnBlk.W[i][j].dSvdU(dRdU,
				  SolnBlk.Grid.Cell[i][j].Xc,
				  SolnBlk.dWdy[i][j]);
	}

	// Calculate the source Jacobian for the turbulent source term if required.
	if (SolnBlk.Flow_Type == FLOWTYPE_TURBULENT_RANS_K_OMEGA) {
	  SolnBlk.W[i][j].dStdU(dRdU,
				SolnBlk.Grid.Cell[i][j].Xc,
				SolnBlk.dWdx[i][j],
				SolnBlk.dWdy[i][j],
				SolnBlk.Axisymmetric);
	}

	// Calculate the source Jacobian for the Variable Prandtl dissipation term if required.
	if (SolnBlk.Variable_Prandtl == ON) {
	  SolnBlk.W[i][j].dSvpdU(dRdU,
				 SolnBlk.Grid.Cell[i][j].Xc,
				 SolnBlk.dWdx[i][j],
				 SolnBlk.dWdy[i][j],
				 SolnBlk.Axisymmetric,
				 SolnBlk.d_dWdx_dW[i][j][0],
				 SolnBlk.d_dWdy_dW[i][j][0],
				 SolnBlk.Wall[i][j].ywall,
				 SolnBlk.Wall[i][j].yplus);
	}

	// Include source Jacobian in the LHS matrix.
	LinSys.A.identity();
	LinSys.A -= (omega*IP.CFL_Number*SolnBlk.dt[i][j])*dRdU;

	// Set the explicit residual as the RHS for the point implicit
	// formulation (already contains the CFL number).
	for (int k = 1; k <= NUM_VAR_NAVIERSTOKES2D; k++)
	  LinSys.b(k-1) = omega*SolnBlk.dUdt[i][j][k_residual][k];

	// Solve system of equations using LU decomposition Gaussian 
	// elimination procedure.
	LinSys.solve(LU_DECOMPOSITION);

	// Update the conserved solution variables.
	for (int k = 1; k <= NUM_VAR_NAVIERSTOKES2D; k++)
	  SolnBlk.U[i][j][k] = SolnBlk.Uo[i][j][k] + LinSys.x(k-1);

      }

      // Perform residual reductions (reduce the time-step) if
      // necessary for scalar or semi-implicit local time-stepping.
      if ((IP.Local_Time_Stepping == SCALAR_LOCAL_TIME_STEPPING ||
	   IP.Local_Time_Stepping == SEMI_IMPLICIT_LOCAL_TIME_STEPPING) &&
	  SolnBlk.U[i][j].Unphysical_Properties()) {
	residual_reduction_factor = ONE;
	for (int n_residual_reduction = 0; n_residual_reduction < 10; n_residual_reduction++) {
	  // Reduce the residual by half.
	  residual_reduction_factor = HALF*residual_reduction_factor;
	  SolnBlk.dt[i][j] = residual_reduction_factor*SolnBlk.dt[i][j];
	  SolnBlk.dUdt[i][j][k_residual] = residual_reduction_factor*SolnBlk.dUdt[i][j][k_residual];
	  // Re-perform the solution update.
	  if (IP.Local_Time_Stepping == SCALAR_LOCAL_TIME_STEPPING) {
	    SolnBlk.U[i][j] = SolnBlk.Uo[i][j] + omega*SolnBlk.dUdt[i][j][k_residual];
	  } else if (IP.Local_Time_Stepping == SEMI_IMPLICIT_LOCAL_TIME_STEPPING) {
	    LinSys.A.identity();
	    LinSys.A -= (omega*IP.CFL_Number*SolnBlk.dt[i][j])*dRdU;
	    for (int k = 1; k <= NUM_VAR_NAVIERSTOKES2D; k++)
	      LinSys.b(k-1) = omega*SolnBlk.dUdt[i][j][k_residual][k];
	    LinSys.solve(LU_DECOMPOSITION);
	    for (int k = 1; k <= NUM_VAR_NAVIERSTOKES2D; k++)
	      SolnBlk.U[i][j][k] = SolnBlk.Uo[i][j][k] + LinSys.x(k-1);
	  }
	  if (!SolnBlk.U[i][j].Unphysical_Properties()) break;
	  if (n_residual_reduction == 10) cout << "n_residual_reductions = " << n_residual_reduction
					       << " @ Xc =" << SolnBlk.Grid.Cell[i][j].Xc << endl;
	}
      }

      // Reset unphysical turbulence state properties.
      if (SolnBlk.Flow_Type == FLOWTYPE_TURBULENT_RANS_K_OMEGA) {
	if (SolnBlk.U[i][j].dk < TOLER) SolnBlk.U[i][j].dk = TOLER;
	if (SolnBlk.U[i][j].domega < TOLER) SolnBlk.U[i][j].domega = TOLER;
      }
      // Reset unphysical turbulence state properties.
      if (SolnBlk.Variable_Prandtl == ON) {
	if (SolnBlk.U[i][j].dke < TOLER*TOLER){
	  SolnBlk.U[i][j].dke = TOLER;
	}
	if (SolnBlk.U[i][j].dee < TOLER*TOLER){
	  SolnBlk.U[i][j].dee = TOLER;
	}
      }

      // Check for unphysical state properties.
      if (SolnBlk.U[i][j].Unphysical_Properties()) {
	cout << "\n " << CFFC_Name() 
	     << " NavierStokes2D ERROR: Negative Density and/or Energy: \n"
	     << " cell = (" << i << "," << j << ") "         << endl
	     << " X    = " << SolnBlk.Grid.Cell[i][j].Xc     << endl
	     << " U    = " << SolnBlk.U[i][j]                << endl 
	     << " W    = " << W(SolnBlk.U[i][j])             << endl
	     << " Uo   = " << SolnBlk.Uo[i][j]               << endl
	     << " Wo   = " << W(SolnBlk.Uo[i][j])            << endl
	     << " dUdt = " << SolnBlk.dUdt[i][j][k_residual] << endl;
	cout.flush();
	return i;
      }

      // Update the primitive variable state.
      SolnBlk.W[i][j] = W(SolnBlk.U[i][j]);

    }
  }

  // Deallocate memory for linear system solver.
  LinSys.deallocate();

  // Solution successfully updated.
  return 0;

}
