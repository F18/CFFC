/**********************************************************************
 * HighTemp2DQuadSingleBlock.cc                                       *
 *                                                                    *
 * Single-block versions of subroutines for 2D High-temp multi-       *
 * block quadrilateral mesh solution classes.                         *
 *                                                                    *
 **********************************************************************/

#include "HighTemp2DQuad.h"

/**********************************************************************
 * HighTemp2D_Quad_Block -- Single Block External Subroutines.        *
 **********************************************************************/

/**********************************************************************
 * Routine: Broadcast_Solution_Block                                  *
 *                                                                    *
 * Broadcast quadrilateral solution block to all processors involved  *
 * in the calculation from the primary processor using the MPI        *
 * broadcast routine.                                                 *
 *                                                                    *
 **********************************************************************/
void Broadcast_Solution_Block(HighTemp2D_Quad_Block &SolnBlk) {

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

  // Broadcast the sub-cell reconstruction indicator.
  //MPI::COMM_WORLD.Bcast(&(SolnBlk.SubCell_Reconstruction),1,MPI::INT,0);

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
    buffer = new double[NUM_VAR_HIGHTEMP2D*ni*nj];

    if (CFFC_Primary_MPI_Processor()) {
      buffer_size = 0;
      for (int j = SolnBlk.JCl-SolnBlk.Nghost; j <= SolnBlk.JCu+SolnBlk.Nghost; j++) {
	for (int i = SolnBlk.ICl-SolnBlk.Nghost; i <= SolnBlk.ICu+SolnBlk.Nghost; i++) {
	  for (int k = 0; k < NUM_VAR_HIGHTEMP2D; k++) {
	    buffer[buffer_size] = SolnBlk.U[i][j][k+1];
	    buffer_size++;
	  }
	}
      }
    }

    buffer_size = NUM_VAR_HIGHTEMP2D*ni*nj;
    MPI::COMM_WORLD.Bcast(buffer,buffer_size,MPI::DOUBLE,0);

    if (!CFFC_Primary_MPI_Processor()) {
      buffer_size = 0;
      for (int j = SolnBlk.JCl-SolnBlk.Nghost; j <= SolnBlk.JCu+SolnBlk.Nghost; j++) {
	for (int i = SolnBlk.ICl-SolnBlk.Nghost; i <= SolnBlk.ICu+SolnBlk.Nghost; i++) {
	  for (int k = 0; k < NUM_VAR_HIGHTEMP2D; k++) {
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
    buffer = new double[2*NUM_VAR_HIGHTEMP2D*ni*nj];

    if (CFFC_Primary_MPI_Processor()) {
      buffer_size = 0;
      for (int j = SolnBlk.JCl-SolnBlk.Nghost; j <= SolnBlk.JCu+SolnBlk.Nghost; j++) {
	for (int k = 0; k < NUM_VAR_HIGHTEMP2D; k++) {
	  buffer[buffer_size] = SolnBlk.WoW[j][k+1];
	  buffer_size++;
	}
	for (int k = 0; k < NUM_VAR_HIGHTEMP2D; k++) {
	  buffer[buffer_size] = SolnBlk.WoE[j][k+1];
	  buffer_size++;
	}
      }
    }

    buffer_size = 2*NUM_VAR_HIGHTEMP2D*ni*nj;
    MPI::COMM_WORLD.Bcast(buffer,buffer_size,MPI::DOUBLE,0);

    if (!CFFC_Primary_MPI_Processor()) {
      buffer_size = 0;
      for (int j = SolnBlk.JCl-SolnBlk.Nghost; j <= SolnBlk.JCu+SolnBlk.Nghost; j++) {
	for (int k = 0; k < NUM_VAR_HIGHTEMP2D; k++) {
	  SolnBlk.WoW[j][k+1] = buffer[buffer_size];
	  buffer_size++;
	}
	for (int k = 0; k < NUM_VAR_HIGHTEMP2D; k++) {
	  SolnBlk.WoE[j][k+1] = buffer[buffer_size];
	  buffer_size++;
	}
      }
    }

    delete []buffer; buffer = NULL;

    ni = (SolnBlk.ICu+SolnBlk.Nghost) - (SolnBlk.ICl-SolnBlk.Nghost) + 1;
    nj = 1;
    buffer = new double[2*NUM_VAR_HIGHTEMP2D*ni*nj];

    if (CFFC_Primary_MPI_Processor()) {
      buffer_size = 0;
      for (int i = SolnBlk.ICl-SolnBlk.Nghost; i <= SolnBlk.ICu+SolnBlk.Nghost; i++) {
	for (int k = 0; k < NUM_VAR_HIGHTEMP2D; k++) {
	  buffer[buffer_size] = SolnBlk.WoS[i][k+1];
	  buffer_size++;
	}
	for (int k = 0; k < NUM_VAR_HIGHTEMP2D; k++) {
	  buffer[buffer_size] = SolnBlk.WoN[i][k+1];
	  buffer_size++;
	}
      }
    }

    buffer_size = 2*NUM_VAR_HIGHTEMP2D*ni*nj;
    MPI::COMM_WORLD.Bcast(buffer,buffer_size,MPI::DOUBLE,0);

    if (!CFFC_Primary_MPI_Processor()) {
      buffer_size = 0;
      for (int i = SolnBlk.ICl-SolnBlk.Nghost; i <= SolnBlk.ICu+SolnBlk.Nghost; i++) {
	for (int k = 0; k < NUM_VAR_HIGHTEMP2D; k++) {
	  SolnBlk.WoS[i][k+1] = buffer[buffer_size];
	  buffer_size++;
	}
	for (int k = 0; k < NUM_VAR_HIGHTEMP2D; k++) {
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
void Broadcast_Solution_Block(HighTemp2D_Quad_Block &SolnBlk,
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

  // Broadcast the sub-cell reconstruction indicator.
  //Communicator.Bcast(&(SolnBlk.SubCell_Reconstruction),1,MPI::INT,Source_Rank);

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
    buffer = new double[NUM_VAR_HIGHTEMP2D*ni*nj];

    if (CFFC_MPI::This_Processor_Number == Source_CPU) {
      buffer_size = 0;
      for (int j = SolnBlk.JCl-SolnBlk.Nghost; j <= SolnBlk.JCu+SolnBlk.Nghost; j++) {
	for (int i = SolnBlk.ICl-SolnBlk.Nghost; i <= SolnBlk.ICu+SolnBlk.Nghost; i++) {
	  for (int k = 0; k < NUM_VAR_HIGHTEMP2D; k++) {
	    buffer[buffer_size] = SolnBlk.U[i][j][k+1];
	    buffer_size++;
	  }
	}
      }
    }

    buffer_size = NUM_VAR_HIGHTEMP2D*ni*nj;
    Communicator.Bcast(buffer,buffer_size,MPI::DOUBLE,Source_Rank);

    if (!(CFFC_MPI::This_Processor_Number == Source_CPU)) {
      buffer_size = 0;
      for (int j = SolnBlk.JCl-SolnBlk.Nghost; j <= SolnBlk.JCu+SolnBlk.Nghost; j++) {
	for (int i = SolnBlk.ICl-SolnBlk.Nghost; i <= SolnBlk.ICu+SolnBlk.Nghost; i++) {
	  for (int k = 0; k < NUM_VAR_HIGHTEMP2D; k++) {
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
    buffer = new double[2*NUM_VAR_HIGHTEMP2D*ni*nj];

    if (CFFC_MPI::This_Processor_Number == Source_CPU) {
      buffer_size = 0;
      for (int j = SolnBlk.JCl-SolnBlk.Nghost; j <= SolnBlk.JCu+SolnBlk.Nghost; j++) {
	for (int k = 0; k < NUM_VAR_HIGHTEMP2D; k++) {
	  buffer[buffer_size] = SolnBlk.WoW[j][k+1];
	  buffer_size++;
	}
	for (int k = 0; k < NUM_VAR_HIGHTEMP2D; k++) {
	  buffer[buffer_size] = SolnBlk.WoW[j][k+1];
	  buffer_size++;
	}
      }
    }

    buffer_size = 2*NUM_VAR_HIGHTEMP2D*ni*nj;
    Communicator.Bcast(buffer,buffer_size,MPI::DOUBLE,Source_Rank);

    if (!(CFFC_MPI::This_Processor_Number == Source_CPU)) {
      buffer_size = 0;
      for (int j  = SolnBlk.JCl-SolnBlk.Nghost; j <= SolnBlk.JCu+SolnBlk.Nghost; j++) {
	for (int k = 0; k < NUM_VAR_HIGHTEMP2D; k++) {
	  SolnBlk.WoW[j][k+1] = buffer[buffer_size];
	  buffer_size++;
	}
	for (int k = 0; k < NUM_VAR_HIGHTEMP2D; k++) {
	  SolnBlk.WoE[j][k+1] = buffer[buffer_size];
	  buffer_size++;
	}
      }
    }

    delete []buffer; buffer = NULL;

    ni = (SolnBlk.ICu+SolnBlk.Nghost) - (SolnBlk.ICl-SolnBlk.Nghost) + 1;
    nj = 1;
    buffer = new double[2*NUM_VAR_HIGHTEMP2D*ni*nj];

    if (CFFC_MPI::This_Processor_Number == Source_CPU) {
      buffer_size = 0;
      for (int i  = SolnBlk.ICl-SolnBlk.Nghost; i <= SolnBlk.ICu+SolnBlk.Nghost; i++) {
	for (int k = 0; k < NUM_VAR_HIGHTEMP2D; k++) {
	  buffer[buffer_size] = SolnBlk.WoS[i][k+1];
	  buffer_size++;
	}
	for (int k = 0; k < NUM_VAR_HIGHTEMP2D; k++) {
	  buffer[buffer_size] = SolnBlk.WoN[i][k+1];
	  buffer_size++;
	}
      }
    }

    buffer_size = 2*NUM_VAR_HIGHTEMP2D*ni*nj;
    Communicator.Bcast(buffer,buffer_size,MPI::DOUBLE,Source_Rank);

    if (!(CFFC_MPI::This_Processor_Number == Source_CPU)) {
      buffer_size = 0;
      for (int i  = SolnBlk.ICl-SolnBlk.Nghost; i <= SolnBlk.ICu+SolnBlk.Nghost; i++) {
	for (int k = 0; k < NUM_VAR_HIGHTEMP2D; k++) {
	  SolnBlk.WoS[i][k+1] = buffer[buffer_size];
	  buffer_size++;
	}
	for (int k = 0; k < NUM_VAR_HIGHTEMP2D; k++) {
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
void Copy_Solution_Block(HighTemp2D_Quad_Block &SolnBlk1,
                         HighTemp2D_Quad_Block &SolnBlk2) {

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

  // Copy the sub-cell reconstruction indicator.
  // SolnBlk1.SubCell_Reconstruction = SolnBlk2.SubCell_Reconstruction;

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
	for (int k = 0; k < NUMBER_OF_RESIDUAL_VECTORS_HIGHTEMP2D; k++)
	  SolnBlk1.dUdt[i][j][k] = SolnBlk2.dUdt[i][j][k];
	SolnBlk1.dWdx[i][j] = SolnBlk2.dWdx[i][j];
	SolnBlk1.dWdy[i][j] = SolnBlk2.dWdy[i][j];
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
int Prolong_Solution_Block(HighTemp2D_Quad_Block &SolnBlk_Fine,
			   HighTemp2D_Quad_Block &SolnBlk_Original,
			   const int Sector) {

  int error_flag;
  int i_min, i_max, j_min, j_max, mesh_refinement_permitted;
  int i_fine, j_fine, i_coarse_min, j_coarse_min, coarse_cell_found;
  double distance, area_total_fine;
  Vector2D dX;
  HighTemp2D_cState HighTemp2D_U_STDATM, Ucoarse;
  HighTemp2D_pState HighTemp2D_W_STDATM;

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

    // Copy the sub-cell reconstruction indicator.
    //SolnBlk_Fine.SubCell_Reconstruction = SolnBlk_Original.SubCell_Reconstruction;

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
 
	if (SolnBlk_Original.U[i][j].rho < ZERO) {
	  if (SolnBlk_Original.U[i-1][j].rho > ZERO && i-1 >= SolnBlk_Original.ICl) Ucoarse = SolnBlk_Original.U[i-1][j];
	  else if (SolnBlk_Original.U[i-1][j-1].rho > ZERO && i-1 >= SolnBlk_Original.ICl && j-1 >= SolnBlk_Original.JCl) Ucoarse = SolnBlk_Original.U[i-1][j-1];
	  else if (SolnBlk_Original.U[i  ][j-1].rho > ZERO                                && j-1 >= SolnBlk_Original.JCl) Ucoarse = SolnBlk_Original.U[i  ][j-1];
	  else if (SolnBlk_Original.U[i+1][j-1].rho > ZERO && i+1 <= SolnBlk_Original.ICu && j-1 >= SolnBlk_Original.JCl) Ucoarse = SolnBlk_Original.U[i+1][j-1];
	  else if (SolnBlk_Original.U[i+1][j  ].rho > ZERO && i+1 <= SolnBlk_Original.ICu                               ) Ucoarse = SolnBlk_Original.U[i+1][j  ];
	  else if (SolnBlk_Original.U[i+1][j+1].rho > ZERO && i+1 <= SolnBlk_Original.ICu && j+1 <= SolnBlk_Original.JCu) Ucoarse = SolnBlk_Original.U[i+1][j+1];
	  else if (SolnBlk_Original.U[i  ][j+1].rho > ZERO                                && j+1 <= SolnBlk_Original.JCu) Ucoarse = SolnBlk_Original.U[i  ][j+1];
	  else if (SolnBlk_Original.U[i-1][j+1].rho > ZERO && i-1 >= SolnBlk_Original.ICl && j+1 <= SolnBlk_Original.JCu) Ucoarse = SolnBlk_Original.U[i-1][j+1];
	  else return 1;
	  SolnBlk_Fine.U[2*(i-i_min)+SolnBlk_Fine.ICl  ][2*(j-j_min)+SolnBlk_Fine.JCl  ] = Ucoarse;
	  //= (SolnBlk_Original.Grid.Cell[i][j].A/area_total_fine)*SolnBlk_Original.U[i][j];
	  SolnBlk_Fine.W[2*(i-i_min)+SolnBlk_Fine.ICl  ][2*(j-j_min)+SolnBlk_Fine.JCl  ] = W(SolnBlk_Fine.U[2*(i-i_min)+SolnBlk_Fine.ICl  ][2*(j-j_min)+SolnBlk_Fine.JCl  ]);
	  SolnBlk_Fine.U[2*(i-i_min)+SolnBlk_Fine.ICl+1][2*(j-j_min)+SolnBlk_Fine.JCl  ] = Ucoarse;
	  //= (SolnBlk_Original.Grid.Cell[i][j].A/area_total_fine)*SolnBlk_Original.U[i][j];
	  SolnBlk_Fine.W[2*(i-i_min)+SolnBlk_Fine.ICl+1][2*(j-j_min)+SolnBlk_Fine.JCl  ] = W(SolnBlk_Fine.U[2*(i-i_min)+SolnBlk_Fine.ICl+1][2*(j-j_min)+SolnBlk_Fine.JCl  ]);
	  SolnBlk_Fine.U[2*(i-i_min)+SolnBlk_Fine.ICl  ][2*(j-j_min)+SolnBlk_Fine.JCl+1] = Ucoarse;
	  //= (SolnBlk_Original.Grid.Cell[i][j].A/area_total_fine)*SolnBlk_Original.U[i][j];
	  SolnBlk_Fine.W[2*(i-i_min)+SolnBlk_Fine.ICl  ][2*(j-j_min)+SolnBlk_Fine.JCl+1] = W(SolnBlk_Fine.U[2*(i-i_min)+SolnBlk_Fine.ICl  ][2*(j-j_min)+SolnBlk_Fine.JCl+1]);
	  SolnBlk_Fine.U[2*(i-i_min)+SolnBlk_Fine.ICl+1][2*(j-j_min)+SolnBlk_Fine.JCl+1] = Ucoarse;
	  //= (SolnBlk_Original.Grid.Cell[i][j].A/area_total_fine)*SolnBlk_Original.U[i][j];
	  SolnBlk_Fine.W[2*(i-i_min)+SolnBlk_Fine.ICl+1][2*(j-j_min)+SolnBlk_Fine.JCl+1] = W(SolnBlk_Fine.U[2*(i-i_min)+SolnBlk_Fine.ICl+1][2*(j-j_min)+SolnBlk_Fine.JCl+1]);
	} 
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
int Restrict_Solution_Block(HighTemp2D_Quad_Block &SolnBlk_Coarse,
			    HighTemp2D_Quad_Block &SolnBlk_Original_SW,
			    HighTemp2D_Quad_Block &SolnBlk_Original_SE,
			    HighTemp2D_Quad_Block &SolnBlk_Original_NW,
			    HighTemp2D_Quad_Block &SolnBlk_Original_NE) {

  int error_flag;
  int i_coarse, j_coarse, mesh_coarsening_permitted;
  double A, total_area;
  Polygon Pc, Pf;
  HighTemp2D_pState HighTemp2D_W_STDATM; HighTemp2D_W_STDATM.Standard_Atmosphere();
  HighTemp2D_cState HighTemp2D_U_STDATM; HighTemp2D_U_STDATM.Standard_Atmosphere();
  HighTemp2D_cState HighTemp2D_U_VACUUM; HighTemp2D_U_VACUUM.Vacuum();

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

    // Copy the sub-cell reconstruction indicator.
    // SolnBlk_Coarse.SubCell_Reconstruction = SolnBlk_Original_SW.SubCell_Reconstruction;

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
	if (SolnBlk_Coarse.U[i_coarse][j_coarse].rho < ZERO) {
	  cout << endl << SolnBlk_Coarse.Grid.Cell[i_coarse][j_coarse].Xc;
	}
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
void ICs(HighTemp2D_Quad_Block &SolnBlk,
         HighTemp2D_Input_Parameters &IP,
         HighTemp2D_pState *Wo) {

  HighTemp2D_pState Wl, Wr;
  double eta, f, fp, fpp;
  HighTemp2D_pState HighTemp2D_W_STDATM; HighTemp2D_W_STDATM.Standard_Atmosphere();

  //set reference Mach number
  Wl.Mref = IP.Mach_Number;
  Wr.Mref = IP.Mach_Number;

  // Assign the initial data for the IVP of interest.
  switch(IP.i_ICs) {
  case IC_CONSTANT :
  case IC_UNIFORM :
    // Set the solution state to the initial state Wo[0].
    for (int j = SolnBlk.JCl-SolnBlk.Nghost; j <= SolnBlk.JCu+SolnBlk.Nghost; j++) {
      for (int i = SolnBlk.ICl-SolnBlk.Nghost; i <= SolnBlk.ICu+SolnBlk.Nghost; i++) {
	SolnBlk.W[i][j] = Wo[0];
       	if (SolnBlk.Flow_Type == FLOWTYPE_TURBULENT_RANS_K_OMEGA) {
	  //SolnBlk.W[i][j].k = max(ONE,0.00005*(sqr(SolnBlk.W[i][j].v.x) + sqr(SolnBlk.W[i][j].v.y)));
	  SolnBlk.W[i][j].k = 0.00005*(sqr(SolnBlk.W[i][j].v.x) + sqr(SolnBlk.W[i][j].v.y));
	  SolnBlk.W[i][j].omega = SolnBlk.W[i][j].k/0.000001;
	}
	SolnBlk.U[i][j] = U(SolnBlk.W[i][j]);
      }
    }
    break;
  case IC_SOD_XDIR :
    // Set the initial data for the Sod IVP in x-direction.
    Wl = HighTemp2D_W_STDATM;
    Wr = HighTemp2D_pState(DENSITY_STDATM/EIGHT,ZERO,ZERO,PRESSURE_STDATM/TEN);
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
    Wl = HighTemp2D_W_STDATM;
    Wr = HighTemp2D_pState(DENSITY_STDATM/EIGHT,ZERO,ZERO,PRESSURE_STDATM/TEN);
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
    Wl = HighTemp2D_pState(4.696,ZERO,ZERO,404.4e03);
    Wr = HighTemp2D_pState(1.408,ZERO,ZERO,101.1e03);
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
    Wl = HighTemp2D_pState(4.696,ZERO,ZERO,404.4e03);
    Wr = HighTemp2D_pState(1.408,ZERO,ZERO,101.1e03);
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
    Wl = HighTemp2D_pState(ONE,-TWO,ZERO,FOUR/TEN);
    Wr = HighTemp2D_pState(ONE, TWO,ZERO,FOUR/TEN);
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
    Wl = HighTemp2D_pState(ONE,ZERO,-TWO,FOUR/TEN);
    Wr = HighTemp2D_pState(ONE,ZERO, TWO,FOUR/TEN);
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
    Wl = HighTemp2D_pState(2.281, 164.83, ZERO, 201.17e03);
    Wr = HighTemp2D_pState(1.408, ZERO, ZERO, 101.1e03);
    //Wl = HighTemp2D_pState(1.3189655, 57.543675, ZERO, 126.65625e03);
    //Wr = HighTemp2D_pState(1.125, ZERO, ZERO, 101.325e03);
    //Wl = HighTemp2D_pState(1.828125, 186.120997, ZERO, 202650.5);
    //Wr = HighTemp2D_pState(1.125, ZERO, ZERO, 101325.0);
    //Wl = Wo[0];
    //Wr = HighTemp2D_W_STDATM;
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
    Wl = HighTemp2D_pState(2.281,ZERO,164.83,201.17e03);
    Wr = HighTemp2D_pState(1.408,ZERO,ZERO,  101.10e03);
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
    Wl = HighTemp2D_pState(1.045,200.0,ZERO,PRESSURE_STDATM);
    Wr = HighTemp2D_pState(3.483,200.0,ZERO,PRESSURE_STDATM);
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
    Wl = HighTemp2D_pState(1.045,ZERO,200.0,PRESSURE_STDATM);
    Wr = HighTemp2D_pState(3.483,ZERO,200.0,PRESSURE_STDATM);
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
    Wl = HighTemp2D_pState(1.598,-383.64,ZERO, 91.88e03);
    Wr = HighTemp2D_pState(2.787,-216.97,ZERO,200.09e03);
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
    Wl = HighTemp2D_pState(1.598,ZERO,-383.64,91.88e03);
    Wr = HighTemp2D_pState(2.787,ZERO,-216.97,200.0e03);
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
		{
		  // double pright = PRESSURE_STDATM*50.0;
    Wl = HighTemp2D_W_STDATM;
    // Wr = HighTemp2D_pState(DENSITY_STDATM*FOUR,0.0,0.0,PRESSURE_STDATM*50.0);
    Wr = HighTemp2D_pState(DENSITY_STDATM*FOUR,0.0,0.0,PRESSURE_STDATM*FOUR);
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
		}
    break;
  case IC_HIGH_PRESSURE_RESERVOIR :
    // Set high-pressure reservoir initial data.
    Wl = HighTemp2D_pState(HUNDRED*DENSITY_STDATM,ZERO,ZERO,HUNDRED*PRESSURE_STDATM);
    Wr = HighTemp2D_W_STDATM;
    for (int j = SolnBlk.JCl-SolnBlk.Nghost; j <= SolnBlk.JCu+SolnBlk.Nghost; j++) {
      for (int i = SolnBlk.ICl-SolnBlk.Nghost; i <= SolnBlk.ICu+SolnBlk.Nghost; i++) {
	if (SolnBlk.Grid.Cell[i][j].Xc.x <= ZERO-0.250+0.05+ZERO*ONE) {//-QUARTER) 
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
    Wl = HighTemp2D_W_STDATM;
    Wr = HighTemp2D_pState(Wo[0].rho,ZERO,ZERO,Wo[0].p);
    for (int j = SolnBlk.JCl-SolnBlk.Nghost; j <= SolnBlk.JCu+SolnBlk.Nghost; j++) {
      for (int i = SolnBlk.ICl-SolnBlk.Nghost; i <= SolnBlk.ICu+SolnBlk.Nghost; i++) {
	if (SolnBlk.Grid.Cell[i][j].Xc.x <= ZERO+0.014) { //ONE+ZERO*
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
    Wl = HighTemp2D_W_STDATM;
    Wr = HighTemp2D_pState(EIGHT*DENSITY_STDATM,ZERO,ZERO,TEN*PRESSURE_STDATM);
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
    Wl = HighTemp2D_pState(EIGHT*DENSITY_STDATM,ZERO,ZERO,TEN*PRESSURE_STDATM);
    Wr = HighTemp2D_W_STDATM;
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
	  SolnBlk.W[i][j] = HighTemp2D_W_STDATM;
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
	  SolnBlk.W[i][j] = HighTemp2D_W_STDATM;
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
	  SolnBlk.W[i][j] = HighTemp2D_W_STDATM;
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
	  SolnBlk.W[i][j] = HighTemp2D_W_STDATM;
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
	  SolnBlk.W[i][j] = HighTemp2D_W_STDATM;
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
	  SolnBlk.W[i][j] = HighTemp2D_W_STDATM;
	}
	SolnBlk.U[i][j] = U(SolnBlk.W[i][j]);
      }
    }
    break;
  case IC_COMPRESSION_XDIR :
    // Set the initial data for a one-dimensional compression problem in the x-direction.
    Wl = HighTemp2D_pState(ONE, 200.0,ZERO,PRESSURE_STDATM);
    Wr = HighTemp2D_pState(ONE,-200.0,ZERO,PRESSURE_STDATM);
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
    Wl = HighTemp2D_pState(ONE,ZERO, 200.0,PRESSURE_STDATM);
    Wr = HighTemp2D_pState(ONE,ZERO,-200.0,PRESSURE_STDATM);
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
    Wl = HighTemp2D_pState(ONE,-200.0,ZERO,PRESSURE_STDATM);
    Wr = HighTemp2D_pState(ONE, 200.0,ZERO,PRESSURE_STDATM);
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
    Wl = HighTemp2D_pState(ONE,ZERO,-200.0,PRESSURE_STDATM);
    Wr = HighTemp2D_pState(ONE,ZERO, 200.0,PRESSURE_STDATM);
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
					      IP.dp,IP.Pipe_Length,IP.Pipe_Radius);
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
	  SolnBlk.W[i][j].k = 0.005*(sqr(SolnBlk.W[i][j].v.x) + sqr(SolnBlk.W[i][j].v.y));
	  SolnBlk.W[i][j].omega = sqrt(SolnBlk.W[i][j].k)/(0.000001);
	  //SolnBlk.W[i][j].omega = 100000.0;
	  //SolnBlk.W[i][j].k = 0.1*(SolnBlk.W[i][j].v.x);
	  //SolnBlk.W[i][j].omega = 10000000.0;
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
  case IC_FORWARD_FACING_STEP :
    for (int j = SolnBlk.JCl-SolnBlk.Nghost; j <= SolnBlk.JCu+SolnBlk.Nghost; j++) {
      for (int i = SolnBlk.ICl-SolnBlk.Nghost; i <= SolnBlk.ICu+SolnBlk.Nghost; i++) {
	//Glaister ICs defined in Input.cc file 
	Wl = Wo[0];
	SolnBlk.W[i][j] = Wo[0];
        SolnBlk.U[i][j] = U(SolnBlk.W[i][j]);
      }
    }
    break;
  case IC_BLUNT_BODY :
    for (int j = SolnBlk.JCl-SolnBlk.Nghost; j <= SolnBlk.JCu+SolnBlk.Nghost; j++) {
      for (int i = SolnBlk.ICl-SolnBlk.Nghost; i <= SolnBlk.ICu+SolnBlk.Nghost; i++) {
	//Glaister ICs defined in Input.cc file 
	Wl = Wo[0];
	SolnBlk.W[i][j] = Wo[0];
        SolnBlk.U[i][j] = U(SolnBlk.W[i][j]);
      }
    }
    break;
  case IC_MIXING_LAYER:
    for (int j = SolnBlk.JCl-SolnBlk.Nghost; j <= SolnBlk.JCu+SolnBlk.Nghost; j++) {
      for (int i = SolnBlk.ICl-SolnBlk.Nghost; i <= SolnBlk.ICu+SolnBlk.Nghost; i++) {
	SolnBlk.W[i][j] = Wo[0];    
	if (SolnBlk.Grid.Cell[i][j].Xc.y >= ZERO) {
	  SolnBlk.W[i][j].v.y = ZERO;
	} else{ 
          SolnBlk.W[i][j].rho = DENSITY_STDATM;
          SolnBlk.W[i][j].v.x = SolnBlk.W[i][j].v.y;
          SolnBlk.W[i][j].v.y = ZERO;
	}
	if (SolnBlk.Flow_Type == FLOWTYPE_TURBULENT_RANS_K_OMEGA) {
	  SolnBlk.W[i][j].k = 0.05*(SolnBlk.W[i][j].v.x);
	  SolnBlk.W[i][j].omega = SolnBlk.W[i][j].rho*SolnBlk.W[i][j].k/(0.000001*SolnBlk.W[i][j].mu());
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
void BCs(HighTemp2D_Quad_Block &SolnBlk, const HighTemp2D_Input_Parameters &IP) {

  Vector2D dX;
  HighTemp2D_pState dW, W;

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
	  SolnBlk.U[SolnBlk.ICl-1][j] = U(SolnBlk.W[SolnBlk.ICl-1][j]);
 	  dX = SolnBlk.Grid.Cell[SolnBlk.ICl-2][j].Xc - SolnBlk.Grid.Cell[SolnBlk.ICl][j].Xc;
 	  SolnBlk.W[SolnBlk.ICl-2][j] = SolnBlk.WoW[j];
	  SolnBlk.W[SolnBlk.ICl-2][j].p = SolnBlk.WoW[j].p - (TWO/THREE)*SolnBlk.W[SolnBlk.ICl][j].dk() + dX.x*IP.dp/IP.Pipe_Length;
	  SolnBlk.W[SolnBlk.ICl-2][j].k = SolnBlk.W[SolnBlk.ICl][j].k;
	  SolnBlk.W[SolnBlk.ICl-2][j].omega = SolnBlk.W[SolnBlk.ICl][j].omega;
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
// 	SolnBlk.W[SolnBlk.ICl-1][j] = Reflect(SolnBlk.W[SolnBlk.ICl][j],SolnBlk.Grid.nfaceW(SolnBlk.ICl,j));
// 	SolnBlk.U[SolnBlk.ICl-1][j] = U(SolnBlk.W[SolnBlk.ICl-1][j]);
// 	SolnBlk.W[SolnBlk.ICl-2][j] = Reflect(SolnBlk.W[SolnBlk.ICl+1][j],SolnBlk.Grid.nfaceW(SolnBlk.ICl,j));
// 	SolnBlk.U[SolnBlk.ICl-2][j] = U(SolnBlk.W[SolnBlk.ICl-2][j]);
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
       SolnBlk.Grid.BCtypeW[SolnBlk.JCl] == BC_BURNING_SURFACE) &&
      (SolnBlk.Grid.BCtypeS[SolnBlk.ICl] == BC_REFLECTION ||
       SolnBlk.Grid.BCtypeS[SolnBlk.ICl] == BC_WALL_VISCOUS_HEATFLUX ||
       SolnBlk.Grid.BCtypeS[SolnBlk.ICl] == BC_WALL_VISCOUS_ISOTHERMAL ||
       SolnBlk.Grid.BCtypeS[SolnBlk.ICl] == BC_MOVING_WALL_HEATFLUX ||
       SolnBlk.Grid.BCtypeS[SolnBlk.ICl] == BC_MOVING_WALL_ISOTHERMAL ||
       SolnBlk.Grid.BCtypeS[SolnBlk.ICl] == BC_BURNING_SURFACE)) {
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
       SolnBlk.Grid.BCtypeW[SolnBlk.JCu] == BC_BURNING_SURFACE) &&
      (SolnBlk.Grid.BCtypeN[SolnBlk.ICl] == BC_REFLECTION ||
       SolnBlk.Grid.BCtypeN[SolnBlk.ICl] == BC_WALL_VISCOUS_HEATFLUX ||
       SolnBlk.Grid.BCtypeN[SolnBlk.ICl] == BC_WALL_VISCOUS_ISOTHERMAL ||
       SolnBlk.Grid.BCtypeN[SolnBlk.ICl] == BC_MOVING_WALL_HEATFLUX ||
       SolnBlk.Grid.BCtypeN[SolnBlk.ICl] == BC_MOVING_WALL_ISOTHERMAL ||
       SolnBlk.Grid.BCtypeN[SolnBlk.ICl] == BC_BURNING_SURFACE)) {
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
       SolnBlk.Grid.BCtypeE[SolnBlk.JCl] == BC_BURNING_SURFACE) &&
      (SolnBlk.Grid.BCtypeS[SolnBlk.ICu] == BC_REFLECTION ||
       SolnBlk.Grid.BCtypeS[SolnBlk.ICu] == BC_WALL_VISCOUS_HEATFLUX ||
       SolnBlk.Grid.BCtypeS[SolnBlk.ICu] == BC_WALL_VISCOUS_ISOTHERMAL ||
       SolnBlk.Grid.BCtypeS[SolnBlk.ICu] == BC_MOVING_WALL_HEATFLUX ||
       SolnBlk.Grid.BCtypeS[SolnBlk.ICu] == BC_MOVING_WALL_ISOTHERMAL ||
       SolnBlk.Grid.BCtypeS[SolnBlk.ICu] == BC_BURNING_SURFACE)) {
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
       SolnBlk.Grid.BCtypeE[SolnBlk.JCu] == BC_BURNING_SURFACE) &&
      (SolnBlk.Grid.BCtypeN[SolnBlk.ICu] == BC_REFLECTION ||
       SolnBlk.Grid.BCtypeN[SolnBlk.ICu] == BC_WALL_VISCOUS_HEATFLUX ||
       SolnBlk.Grid.BCtypeN[SolnBlk.ICu] == BC_WALL_VISCOUS_ISOTHERMAL ||
       SolnBlk.Grid.BCtypeN[SolnBlk.ICu] == BC_MOVING_WALL_HEATFLUX ||
       SolnBlk.Grid.BCtypeN[SolnBlk.ICu] == BC_MOVING_WALL_ISOTHERMAL ||
       SolnBlk.Grid.BCtypeN[SolnBlk.ICu] == BC_BURNING_SURFACE)) {
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
double CFL(HighTemp2D_Quad_Block &SolnBlk,
	   HighTemp2D_Input_Parameters &IP) {

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
void Set_Global_TimeStep(HighTemp2D_Quad_Block &SolnBlk,
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
double L1_Norm_Residual(HighTemp2D_Quad_Block &SolnBlk) {

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
double L2_Norm_Residual(HighTemp2D_Quad_Block &SolnBlk) {

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
double Max_Norm_Residual(HighTemp2D_Quad_Block &SolnBlk) {

  double max_norm = ZERO;

  for (int j = SolnBlk.JCl; j <= SolnBlk.JCu; j++) {
    for (int i = SolnBlk.ICl; i <= SolnBlk.ICu; i++) {
      max_norm = max(max_norm,fabs(SolnBlk.dUdt[i][j][0][SolnBlk.residual_variable]));
    }
  }

  return max_norm;

}

/**********************************************************************
 * Routine: L1_Norm_Residual                                          *
 *                                                                    *
 * Determines the L1-norm of the solution residual for the specified  *
 * quadrilateral solution block. Useful for monitoring convergence of *
 * the solution for steady state problems.                            *
 *                                                                    *
 **********************************************************************/
double L1_Norm_Residual(HighTemp2D_Quad_Block &SolnBlk, int rv) {

  double l1_norm = ZERO;

  for (int j = SolnBlk.JCl; j <= SolnBlk.JCu; j++) {
    for (int i = SolnBlk.ICl; i <= SolnBlk.ICu; i++) {
      l1_norm += fabs(SolnBlk.dUdt[i][j][0][rv]);
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
double L2_Norm_Residual(HighTemp2D_Quad_Block &SolnBlk, int rv) {

  double l2_norm = ZERO;

  for (int j = SolnBlk.JCl; j <= SolnBlk.JCu; j++) {
    for (int i = SolnBlk.ICl; i <= SolnBlk.ICu; i++) {
      l2_norm += sqr(SolnBlk.dUdt[i][j][0][rv]);
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
double Max_Norm_Residual(HighTemp2D_Quad_Block &SolnBlk, int rv) {

  double max_norm = ZERO;

  for (int j = SolnBlk.JCl; j <= SolnBlk.JCu; j++) {
    for (int i = SolnBlk.ICl; i <= SolnBlk.ICu; i++) {
      max_norm = max(max_norm,fabs(SolnBlk.dUdt[i][j][0][rv]));
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
void Linear_Reconstruction_GreenGauss(HighTemp2D_Quad_Block &SolnBlk,
				      const int i,
                                      const int j,
                                      const int Limiter) {

  int error_flag;
  int n_pts, i_index[8], j_index[8];
  double u0Min, u0Max, uQuad[4], phi;
  double DxDx_ave, DxDy_ave, DyDy_ave;
  double l_north, l_south, l_east, l_west;
  Vector2D n_north, n_south, n_east, n_west, dX;
  HighTemp2D_pState W_nw, W_ne, W_sw, W_se, W_face,
                 DU, DUDx_ave, DUDy_ave;
  HighTemp2D_pState HighTemp2D_W_VACUUM; HighTemp2D_W_VACUUM.Vacuum();

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
      DUDx_ave = HighTemp2D_W_VACUUM;
      DUDy_ave = HighTemp2D_W_VACUUM;
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
      for (int n = 1; n <= NUM_VAR_HIGHTEMP2D; n++) {
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
    SolnBlk.dWdx[i][j] = HighTemp2D_W_VACUUM;
    SolnBlk.dWdy[i][j] = HighTemp2D_W_VACUUM;
    SolnBlk.phi[i][j] = HighTemp2D_W_VACUUM;
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
void Linear_Reconstruction_GreenGauss(HighTemp2D_Quad_Block &SolnBlk,
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
void Linear_Reconstruction_LeastSquares(HighTemp2D_Quad_Block &SolnBlk,
				        const int i,
                                        const int j,
                                        const int Limiter) {

  int n_pts, i_index[8], j_index[8], motion_flag = OFF, n_pts_temp;
  double u0Min, u0Max, uQuad[4], phi;
  double DxDx_ave, DxDy_ave, DyDy_ave;
  Vector2D dX;
  HighTemp2D_pState DU, DUDx_ave, DUDy_ave;
  HighTemp2D_pState HighTemp2D_W_VACUUM; HighTemp2D_W_VACUUM.Vacuum();

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
    DUDx_ave = HighTemp2D_W_VACUUM;
    DUDy_ave = HighTemp2D_W_VACUUM;
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
      for (int n = 1; n <= NUM_VAR_HIGHTEMP2D; n++) {
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
    SolnBlk.dWdx[i][j] = HighTemp2D_W_VACUUM;
    SolnBlk.dWdy[i][j] = HighTemp2D_W_VACUUM;
    SolnBlk.phi[i][j]  = HighTemp2D_W_VACUUM;

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
void Linear_Reconstruction_LeastSquares_2(HighTemp2D_Quad_Block &SolnBlk,
				          const int i,
                                          const int j,
                                          const int Limiter) {

  int n_pts, i_index[8], j_index[8];
  double u0Min, u0Max, uQuad[4], phi;
  double DxDx_ave, DxDy_ave, DyDy_ave;
  Vector2D dX;
  HighTemp2D_pState DU, DUDx_ave, DUDy_ave;
  HighTemp2D_pState HighTemp2D_W_VACUUM; HighTemp2D_W_VACUUM.Vacuum();

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
    DUDx_ave = HighTemp2D_W_VACUUM;
    DUDy_ave = HighTemp2D_W_VACUUM;
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
      for (int n = 1; n <= NUM_VAR_HIGHTEMP2D; n++) {
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
    SolnBlk.dWdx[i][j] = HighTemp2D_W_VACUUM;
    SolnBlk.dWdy[i][j] = HighTemp2D_W_VACUUM;
    SolnBlk.phi[i][j] = HighTemp2D_W_VACUUM;

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
void Linear_Reconstruction_LeastSquares(HighTemp2D_Quad_Block &SolnBlk,
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
void Residual_Smoothing(HighTemp2D_Quad_Block &SolnBlk,
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
				   HighTemp2D_Input_Parameters &IP,
                                   int &number_refinement_criteria,
                                   HighTemp2D_Quad_Block &SolnBlk) {

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
void Fix_Refined_Block_Boundaries(HighTemp2D_Quad_Block &SolnBlk,
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
void Unfix_Refined_Block_Boundaries(HighTemp2D_Quad_Block &SolnBlk) {

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
void Apply_Boundary_Flux_Corrections(HighTemp2D_Quad_Block &SolnBlk,
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
void Apply_Boundary_Flux_Corrections_Multistage_Explicit(HighTemp2D_Quad_Block &SolnBlk,
                                                         const int i_stage,
                                                         HighTemp2D_Input_Parameters &IP,
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
int dUdt_Residual_Evaluation(HighTemp2D_Quad_Block &SolnBlk,
			     HighTemp2D_Input_Parameters &IP) {

  int error_flag;
  
  HighTemp2D_pState Wl, Wr, Wu, Wd, Wface;
  HighTemp2D_pState dWdx, dWdy, dWdxl, dWdyl, dWdxr, dWdyr;
  HighTemp2D_cState Flux;
  
  Vector2D dX, Xl, Xr, Xu, Xd;

  HighTemp2D_cState HighTemp2D_U_VACUUM; HighTemp2D_U_VACUUM.Vacuum();
  HighTemp2D_pState HighTemp2D_W_STDATM; HighTemp2D_W_STDATM.Standard_Atmosphere();

  HighTemp2D_pState BlankState; BlankState.Vacuum();
  Vector2D BlankVector; BlankVector.zero();

  if (IP.i_Limiter == LIMITER_ZERO) {
     // Zeroing phi does not need to be done on every call to
     // dUdt_Residual_Evaluation(). But the limiter type can change during
     // the simulation so to be safe we simply zero phi everytime.
     for (int i = 0; i < SolnBlk.NCi; i++) {
	for (int j = 0; j < SolnBlk.NCj; j++) {
	    SolnBlk.phi[i][j].Vacuum();
	}
     }
  } else {
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
	}
  } /* endif */

  // Evaluate the time rate of change of the solution (i.e., the
  // solution residuals) using a second-order limited upwind scheme
  // with a variety of flux functions.

  // Add i-direction (zeta-direction) fluxes.
  for (int j = SolnBlk.JCl-1; j <= SolnBlk.JCu+1; j++) {

    SolnBlk.dUdt[SolnBlk.ICl-1][j][0] = HighTemp2D_U_VACUUM;

    for (int i = SolnBlk.ICl-1; i <= SolnBlk.ICu; i++) {

      SolnBlk.dUdt[i+1][j][0] = HighTemp2D_U_VACUUM;

      if (j >= SolnBlk.JCl && j <= SolnBlk.JCu) {

	if (i == SolnBlk.ICl-1 && 
	    (SolnBlk.Grid.BCtypeW[j] == BC_REFLECTION ||
	     SolnBlk.Grid.BCtypeW[j] == BC_BURNING_SURFACE ||
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
	    Wl = RinglebFlow(Wr,SolnBlk.Grid.xfaceW(i+1,j));
	  } else {
	    // WEST face of cell (i+1,j) is a CHARACTERISTIC boundary.
	    Wl = BC_Characteristic_Pressure(Wr,SolnBlk.WoW[j],SolnBlk.Grid.nfaceW(i+1,j));
	  }
	} else if (i == SolnBlk.ICu &&
		   (SolnBlk.Grid.BCtypeE[j] == BC_REFLECTION ||
		    SolnBlk.Grid.BCtypeE[j] == BC_BURNING_SURFACE ||
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
	    Wr = RinglebFlow(Wl,SolnBlk.Grid.xfaceE(i,j));
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
	case FLUX_FUNCTION_GLAISTER :
	  Flux = FluxGlaister_n(Wl,Wr,SolnBlk.Grid.nfaceE(i,j));
	  break;
	case FLUX_FUNCTION_RUSANOV :
	  Flux = FluxRusanov_n(Wl,Wr,SolnBlk.Grid.nfaceE(i,j));
	  break;
	case FLUX_FUNCTION_HLLE :
	  Flux = FluxHLLE_n(Wl,Wr,SolnBlk.Grid.nfaceE(i,j));
	  break;
        case FLUX_FUNCTION_GHLLE :
	  Flux = FluxGHLLE_n(Wl,Wr,SolnBlk.Grid.nfaceE(i,j));
	  break;
	case FLUX_FUNCTION_HLLL :
	  Flux = FluxHLLL_n(Wl,Wr,SolnBlk.Grid.nfaceE(i,j));
	  break;
	case FLUX_FUNCTION_GHLLL :
	  Flux = FluxGHLLL_n(Wl,Wr,SolnBlk.Grid.nfaceE(i,j));
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
	case FLUX_FUNCTION_AUSM_PLUS_UP :
	  Flux = FluxAUSMplusUP_n(Wl,Wr,SolnBlk.Grid.nfaceE(i,j));
	  break;  
	default:
	  //Flux = FluxRoe_n(Wl,Wr,SolnBlk.Grid.nfaceE(i,j));
	  Flux = FluxGlaister_n(Wl,Wr,SolnBlk.Grid.nfaceE(i,j));
	  break;
	};

	// Compute the cell centred stress tensor and heat flux vector if required.
 	if (SolnBlk.Flow_Type == FLOWTYPE_TURBULENT_RANS_K_OMEGA ||
	    (SolnBlk.Flow_Type && SolnBlk.Axisymmetric)) {
	  SolnBlk.W[i][j].ComputeViscousTerms(SolnBlk.dWdx[i][j],
					      SolnBlk.dWdy[i][j],
					      SolnBlk.Grid.Cell[i][j].Xc,
					      SolnBlk.Axisymmetric);
	  SolnBlk.U[i][j].tau = SolnBlk.W[i][j].tau;
	  SolnBlk.U[i][j].q = SolnBlk.W[i][j].q;
	}

	if (SolnBlk.Flow_Type) {
	   // Determine the i-direction viscous flux.
	   //
	   // At the end of this viscous flux calculation we will call
	   // 
	   //   Flux -= ViscousFlux_n(Xface, Wface, dWdx, dWdy, 
	   //                         norm_dir, SolnBlk.Axisymmetric);
	   //
	   // Xface and norm_dir are simply:
	   Vector2D Xface = SolnBlk.Grid.xfaceE(i, j);
	   Vector2D norm_dir = SolnBlk.Grid.nfaceE(i, j);
	   // regardless of the gradient reconstruction method and
	   // regardless of whether we are on the boundary. Axisymmetric
	   // clearly does not change over the entire simulation.
	   // 
	   // Thus the only purpose of the code between here and the call
	   // to "ViscousFlux_n()" is to determine just the other three
	   // variables: Wface (the solution on the cell interface) and
	   // dWd[xy] (the two components of the solution gradients also
	   // on the face).
	   // 
	   // For Wface: next to a wall, Wface is influenced only by the
	   // type of boundary condition (and not by the type of gradient
	   // reconstruction). Away from a wall, Wface is found using a
	   // bilinear interpolation if we are using diamond-path
	   // reconstruction for the gradients and is found as a simple
	   // average otherwise.
	   // 
	   // For dWd[xy]: first consider the code for cells away from a
	   // boundary which shows fairly well how dWd[xy] is calculated.
	   // Then consider the assumptions near a wall for the different
	   // gradient reconstruction methods.
   	   if (i == SolnBlk.ICl-1 && 
	       (SolnBlk.Grid.BCtypeW[j] == BC_WALL_VISCOUS_HEATFLUX ||
		SolnBlk.Grid.BCtypeW[j] == BC_WALL_VISCOUS_ISOTHERMAL ||
		SolnBlk.Grid.BCtypeW[j] == BC_MOVING_WALL_HEATFLUX ||
		SolnBlk.Grid.BCtypeW[j] == BC_MOVING_WALL_ISOTHERMAL ||
		SolnBlk.Grid.BCtypeW[j] == BC_BURNING_SURFACE)) {
	      // WEST face of cell (i+1,j) is a normal boundary.
	      Xr = SolnBlk.Grid.Cell[i+1][j].Xc; Wr = SolnBlk.W[i+1][j];
	      if (SolnBlk.Grid.BCtypeW[j] == BC_BURNING_SURFACE) {
		 Wu = BurningSurface(Wr, norm_dir);
		 Wd = Wu;
	      } else {
 		 Wu = HALF*(Wr+SolnBlk.W[i+1][j+1]);
		 Wd = HALF*(Wr+SolnBlk.W[i+1][j-1]);
  		 switch (SolnBlk.Grid.BCtypeW[j]) {
		   case BC_WALL_VISCOUS_HEATFLUX:
		     Wu.v.zero(); Wd.v.zero();
		     break;
		   case BC_WALL_VISCOUS_ISOTHERMAL:
		     Wu.v.zero(); Wd.v.zero();
		     switch (HighTemp2D_cState::eos_type) {
			case EOS_TGAS:
			  Wu.rho = Tgas_rho(Wu.p, SolnBlk.Twall);
			  Wd.rho = Tgas_rho(Wd.p, SolnBlk.Twall);
			  break;
			case EOS_IDEAL:
			  Wu.rho = Wu.p / HighTemp2D_cState::R / SolnBlk.Twall;
			  Wd.rho = Wd.p / HighTemp2D_cState::R / SolnBlk.Twall;
			  break;
		     }
		     break;
		   case BC_MOVING_WALL_HEATFLUX:
		     Wu.v = SolnBlk.Vwall; Wd.v = SolnBlk.Vwall;
		     break;
		   case BC_MOVING_WALL_ISOTHERMAL:
		     Wu.v = SolnBlk.Vwall; Wd.v = SolnBlk.Vwall; 
		     switch (HighTemp2D_cState::eos_type) {
			case EOS_TGAS:
			  Wu.rho = Tgas_rho(Wu.p, SolnBlk.Twall);
			  Wd.rho = Tgas_rho(Wd.p, SolnBlk.Twall);
			  break;
			case EOS_IDEAL:
			  Wu.rho = Wu.p / HighTemp2D_cState::R / SolnBlk.Twall;
			  Wd.rho = Wd.p / HighTemp2D_cState::R / SolnBlk.Twall;
			  break;
		     }
		     break;
		 }
	      } /* endif */
  	      Wface = HALF * (Wu + Wd);
 	      switch (IP.i_Viscous_Reconstruction) {
		case VISCOUS_RECONSTRUCTION_CARTESIAN:
		case VISCOUS_RECONSTRUCTION_DIAMOND_PATH: {
		  int stencil_flag = DIAMONDPATH_RIGHT_TRIANGLE;
		  Xu = SolnBlk.Grid.Node[i+1][j+1].X;
		  Xd = SolnBlk.Grid.Node[i+1][j  ].X; 
		  DiamondPath_Find_dWdX(dWdx, dWdy,
					BlankVector, BlankState, 
					Xd, Wd,
					Xr, Wr,
					Xu, Wu,
					stencil_flag);
		  }
		  break;
		case VISCOUS_RECONSTRUCTION_HYBRID :
		  dWdxr = SolnBlk.dWdx[i+1][j]; dWdyr = SolnBlk.dWdy[i+1][j];
                  Xl = Xface;
		  Wl = Wface;
		  dWdxl = dWdxr; dWdyl = dWdyr; // Hmm. 
 		  Hybrid_Find_dWdX(dWdx, dWdy, 
			 	   Xl, Wl, dWdxl, dWdyl,
				   Xr, Wr, dWdxr, dWdyr,
				   norm_dir);
		  break;
		case VISCOUS_RECONSTRUCTION_ARITHMETIC_AVERAGE :
		  dWdx = SolnBlk.dWdx[i+1][j]; dWdy = SolnBlk.dWdy[i+1][j];
		  break;
	      }
	   } else if (i == SolnBlk.ICu &&
		      (SolnBlk.Grid.BCtypeE[j] == BC_WALL_VISCOUS_HEATFLUX ||
		       SolnBlk.Grid.BCtypeE[j] == BC_WALL_VISCOUS_ISOTHERMAL ||
		       SolnBlk.Grid.BCtypeE[j] == BC_MOVING_WALL_HEATFLUX ||
		       SolnBlk.Grid.BCtypeE[j] == BC_MOVING_WALL_ISOTHERMAL ||
		       SolnBlk.Grid.BCtypeE[j] == BC_BURNING_SURFACE)) {
  	      // EAST face of cell (i,j) is a normal boundary.
	      Xl = SolnBlk.Grid.Cell[i][j].Xc; Wl = SolnBlk.W[i][j];
	      if (SolnBlk.Grid.BCtypeE[j] == BC_BURNING_SURFACE) {
		 Wu = BurningSurface(Wl, norm_dir);
		 Wd = Wu;
	      } else {
		 Wu = HALF*(Wl+SolnBlk.W[i][j+1]);
		 Wd = HALF*(Wl+SolnBlk.W[i][j-1]);
		 switch (SolnBlk.Grid.BCtypeE[j]) {
					case BC_WALL_VISCOUS_HEATFLUX:
						Wu.v.zero(); Wd.v.zero();
						break;
					case BC_WALL_VISCOUS_ISOTHERMAL:
						Wu.v.zero(); Wd.v.zero();
						switch (HighTemp2D_cState::eos_type) {
							case EOS_TGAS:
								Wu.rho = Tgas_rho(Wu.p, SolnBlk.Twall);
								Wd.rho = Tgas_rho(Wd.p, SolnBlk.Twall);
								break;
							case EOS_IDEAL:
								Wu.rho = Wu.p / HighTemp2D_cState::R / SolnBlk.Twall;
								Wd.rho = Wd.p / HighTemp2D_cState::R / SolnBlk.Twall;
								break;
						}
						break;
					case BC_MOVING_WALL_HEATFLUX:
						Wu.v = SolnBlk.Vwall; Wd.v = SolnBlk.Vwall;
						break;
					case BC_MOVING_WALL_ISOTHERMAL:
						Wu.v = SolnBlk.Vwall; Wd.v = SolnBlk.Vwall;
						switch (HighTemp2D_cState::eos_type) {
							case EOS_TGAS:
								Wu.rho = Tgas_rho(Wu.p, SolnBlk.Twall);
								Wd.rho = Tgas_rho(Wd.p, SolnBlk.Twall);
								break;
							case EOS_IDEAL:
								Wu.rho = Wu.p / HighTemp2D_cState::R / SolnBlk.Twall;
								Wd.rho = Wd.p / HighTemp2D_cState::R / SolnBlk.Twall;
								break;
						}
						break;
				}
			}

			Wface = HALF * (Wu + Wd);

			switch(IP.i_Viscous_Reconstruction) {
				case VISCOUS_RECONSTRUCTION_CARTESIAN :
				case VISCOUS_RECONSTRUCTION_DIAMOND_PATH :
					{
						int stencil_flag = DIAMONDPATH_LEFT_TRIANGLE;
						Xu = SolnBlk.Grid.Node[i+1][j+1].X;
						Xd = SolnBlk.Grid.Node[i+1][j  ].X; 
						DiamondPath_Find_dWdX(dWdx, dWdy,
								Xl, Wl,
								Xd, Wd,
								BlankVector, BlankState,
								Xu, Wu,
								stencil_flag);
					}
					break;
				case VISCOUS_RECONSTRUCTION_HYBRID :
					dWdxl = SolnBlk.dWdx[i][j]; dWdyl = SolnBlk.dWdy[i][j]; 

					Xr = Xface;
					Wr = Wface;
					dWdxr = dWdxl; dWdyr = dWdyl; 

					Hybrid_Find_dWdX(dWdx, dWdy, 
							Xl, Wl, dWdxl, dWdyl,
							Xr, Wr, dWdxr, dWdyr,
							norm_dir);
					break;
				case VISCOUS_RECONSTRUCTION_ARITHMETIC_AVERAGE :
					dWdx = SolnBlk.dWdx[i][j]; dWdy = SolnBlk.dWdy[i][j];
					break;
			}
		} else {
			// EAST face is either a normal cell or possibly 
			// a non-viscous boundary condition.
			Xl = SolnBlk.Grid.Cell[i  ][j].Xc; Wl = SolnBlk.W[i  ][j];
			Xr = SolnBlk.Grid.Cell[i+1][j].Xc; Wr = SolnBlk.W[i+1][j];

			if (IP.i_Viscous_Reconstruction == VISCOUS_RECONSTRUCTION_CARTESIAN ||
					IP.i_Viscous_Reconstruction == VISCOUS_RECONSTRUCTION_DIAMOND_PATH) {

				int stencil_flag = DIAMONDPATH_QUADRILATERAL_RECONSTRUCTION;
				Xu = SolnBlk.Grid.Node[i+1][j+1].X; 
				Xd = SolnBlk.Grid.Node[i+1][j  ].X; 

				// These two (WnNE and WnSE) are what make diamond path so expensive.
				Wu = SolnBlk.WnNE(i,j); 
				Wd = SolnBlk.WnSE(i,j);

				DiamondPath_Find_dWdX(dWdx, dWdy,
						Xl, Wl, Xd, Wd, Xr, Wr, Xu, Wu, 
						stencil_flag);

				// Find Wface:
				error_flag = Bilinear_Interpolation_ZY(Wl, Xl, Wu, Xu, Wr, Xr, Wd, Xd,
						Xface, Wface);
				// and if error_flag?

			} else {
				Wface = HALF*(Wl + Wr);

				dWdxl = SolnBlk.dWdx[i][j]; dWdxr = SolnBlk.dWdx[i+1][j];
				dWdyl = SolnBlk.dWdy[i][j]; dWdyr = SolnBlk.dWdy[i+1][j];
				
				if (IP.i_Viscous_Reconstruction == VISCOUS_RECONSTRUCTION_HYBRID) {
					Hybrid_Find_dWdX(dWdx, dWdy, 
							Xl, Wl, dWdxl, dWdyl,
							Xr, Wr, dWdxr, dWdyr,
							norm_dir);
				} else if (IP.i_Viscous_Reconstruction == VISCOUS_RECONSTRUCTION_ARITHMETIC_AVERAGE) {
					dWdx = HALF*(dWdxl + dWdxr);
					dWdy = HALF*(dWdyl + dWdyr);
				}
			}
		}
		// Compute the i-direction viscous flux.
		Flux -= ViscousFlux_n(Xface, Wface, dWdx, dWdy, norm_dir, SolnBlk.Axisymmetric);

		if (IP.Solver_Type == IMPLICIT && SolnBlk.face_grad_arrays_allocated) {
			SolnBlk.dWdx_faceE[i][j] = dWdx;
			SolnBlk.dWdy_faceE[i][j] = dWdy;
		}
	} // if (SolnBlk.Flow_Type) 

	// Evaluate cell-averaged solution changes.
	SolnBlk.dUdt[i  ][j][0] -= Flux*SolnBlk.Grid.lfaceE(i,j)/SolnBlk.Grid.Cell[i][j].A;
	SolnBlk.dUdt[i+1][j][0] += Flux*SolnBlk.Grid.lfaceW(i+1,j)/SolnBlk.Grid.Cell[i+1][j].A;

	// Include axisymmetric source terms if required.
	if (SolnBlk.Axisymmetric) {
	  SolnBlk.dUdt[i][j][0] += SolnBlk.W[i][j].Si(SolnBlk.Grid.Cell[i][j].Xc);
	  if (SolnBlk.Flow_Type)
	    SolnBlk.dUdt[i][j][0] += SolnBlk.W[i][j].Sv(SolnBlk.Grid.Cell[i][j].Xc,
							SolnBlk.dWdy[i][j]);
	}

	// Include turbulent production and destruction source term.
	if (SolnBlk.Flow_Type == FLOWTYPE_TURBULENT_RANS_K_OMEGA) {
	  SolnBlk.dUdt[i][j][0] += SolnBlk.W[i][j].St(SolnBlk.Grid.Cell[i][j].Xc,
						      SolnBlk.dWdx[i][j],
						      SolnBlk.dWdy[i][j],
						      SolnBlk.Axisymmetric);
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
      SolnBlk.dUdt[SolnBlk.ICl-1][j][0] = HighTemp2D_U_VACUUM;
      SolnBlk.dUdt[SolnBlk.ICu+1][j][0] = HighTemp2D_U_VACUUM;
    }

  }

  // Add j-direction (eta-direction) fluxes.
  for (int i = SolnBlk.ICl; i <= SolnBlk.ICu; i++) {
    for (int j = SolnBlk.JCl-1; j <= SolnBlk.JCu; j++) {

      // Evaluate the cell interface j-direction INVISCID fluxes.
      if (j == SolnBlk.JCl-1 && 
	  (SolnBlk.Grid.BCtypeS[i] == BC_REFLECTION ||
	   SolnBlk.Grid.BCtypeS[i] == BC_BURNING_SURFACE ||
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
	  Wl = RinglebFlow(Wr,SolnBlk.Grid.xfaceS(i,j+1));
	  Wl = BC_Characteristic_Pressure(Wr,Wl,SolnBlk.Grid.nfaceS(i,j+1));
	} 
	else if (SolnBlk.Grid.BCtypeS[i] == BC_CHARACTERISTIC) {
	  // SOUTH face of cell (i,j+1) is a CHARACTERISTIC boundary.
	  Wl = BC_Characteristic_Pressure(Wr,SolnBlk.WoS[i],SolnBlk.Grid.nfaceS(i,j+1));
	  } 

      } else if (j == SolnBlk.JCu && 
		 (SolnBlk.Grid.BCtypeN[i] == BC_REFLECTION ||
		  SolnBlk.Grid.BCtypeN[i] == BC_BURNING_SURFACE ||
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
	  Wr = RinglebFlow(Wl,SolnBlk.Grid.xfaceN(i,j));
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
      case FLUX_FUNCTION_GLAISTER :
	Flux = FluxGlaister_n(Wl,Wr,SolnBlk.Grid.nfaceN(i,j));
	break;
      case FLUX_FUNCTION_RUSANOV :
	Flux = FluxRusanov_n(Wl,Wr,SolnBlk.Grid.nfaceN(i,j));
	break;
      case FLUX_FUNCTION_HLLE :
	Flux = FluxHLLE_n(Wl,Wr,SolnBlk.Grid.nfaceN(i,j));
	break;
      case FLUX_FUNCTION_GHLLE :
	Flux = FluxGHLLE_n(Wl,Wr,SolnBlk.Grid.nfaceN(i,j));
	break;	
      case FLUX_FUNCTION_HLLL :
	Flux = FluxHLLL_n(Wl,Wr,SolnBlk.Grid.nfaceN(i,j));
	break;
      case FLUX_FUNCTION_GHLLL :
	Flux = FluxGHLLL_n(Wl,Wr,SolnBlk.Grid.nfaceN(i,j));
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
      case FLUX_FUNCTION_AUSM_PLUS_UP :
	Flux = FluxAUSMplusUP_n(Wl,Wr,SolnBlk.Grid.nfaceN(i,j));
	break;
     default:
       //Flux = FluxRoe_n(Wl,Wr,SolnBlk.Grid.nfaceN(i,j));
	Flux = FluxGlaister_n(Wl,Wr,SolnBlk.Grid.nfaceN(i,j));
	break;
      };

		if (SolnBlk.Flow_Type) {

			// Determine the j-direction viscous flux.

			// Please rotate your head 90 degrees to the left before
			// attempting to read this code. For example, Xr is the
			// position of the cell to the north of the face in
			// question. This is necessary, for example, for
			// "DiamondPath_Find_dWdX()", which only deals with "left" 
			// and "right" triangles.

			// See the discussion at the start of code for the 
			// i-direction viscous flux for more on the evaluation 
			// of the viscous flux.

			Vector2D Xface = SolnBlk.Grid.xfaceN(i, j);
			Vector2D norm_dir = SolnBlk.Grid.nfaceN(i, j);

			if (j == SolnBlk.JCl-1 && 
					(SolnBlk.Grid.BCtypeS[i] == BC_WALL_VISCOUS_HEATFLUX ||
					 SolnBlk.Grid.BCtypeS[i] == BC_WALL_VISCOUS_ISOTHERMAL ||
					 SolnBlk.Grid.BCtypeS[i] == BC_MOVING_WALL_HEATFLUX ||
					 SolnBlk.Grid.BCtypeS[i] == BC_MOVING_WALL_ISOTHERMAL ||
					 SolnBlk.Grid.BCtypeS[i] == BC_BURNING_SURFACE)) {

				// SOUTH face of cell (i,j+1) is a normal boundary.

				Xr = SolnBlk.Grid.Cell[i][j+1].Xc; Wr = SolnBlk.W[i][j+1];

				if (SolnBlk.Grid.BCtypeS[i] == BC_BURNING_SURFACE) {
					Wu = BurningSurface(Wr, norm_dir);
					Wd = Wu;
				} else {

					Wu = HALF*(Wr+SolnBlk.W[i-1][j+1]);
					Wd = HALF*(Wr+SolnBlk.W[i+1][j+1]);

					switch (SolnBlk.Grid.BCtypeS[i]) {
						case BC_WALL_VISCOUS_HEATFLUX:
							Wu.v.zero(); Wd.v.zero();
							break;
						case BC_WALL_VISCOUS_ISOTHERMAL:
							Wu.v.zero(); Wd.v.zero();
							switch (HighTemp2D_cState::eos_type) {
								case EOS_TGAS:
									Wu.rho = Tgas_rho(Wu.p, SolnBlk.Twall);
									Wd.rho = Tgas_rho(Wd.p, SolnBlk.Twall);
									break;
								case EOS_IDEAL:
									Wu.rho = Wu.p / HighTemp2D_cState::R / SolnBlk.Twall;
									Wd.rho = Wd.p / HighTemp2D_cState::R / SolnBlk.Twall;
									break;
							}
							break;
						case BC_MOVING_WALL_HEATFLUX:
							Wu.v = SolnBlk.Vwall; Wd.v = SolnBlk.Vwall;
							break;
						case BC_MOVING_WALL_ISOTHERMAL:
							Wu.v = SolnBlk.Vwall; Wd.v = SolnBlk.Vwall; 
							switch (HighTemp2D_cState::eos_type) {
								case EOS_TGAS:
									Wu.rho = Tgas_rho(Wu.p, SolnBlk.Twall);
									Wd.rho = Tgas_rho(Wd.p, SolnBlk.Twall);
									break;
								case EOS_IDEAL:
									Wu.rho = Wu.p / HighTemp2D_cState::R / SolnBlk.Twall;
									Wd.rho = Wd.p / HighTemp2D_cState::R / SolnBlk.Twall;
									break;
							}
							break;
					}
				}

				Wface = HALF * (Wu + Wd);

				switch (IP.i_Viscous_Reconstruction) {
					case VISCOUS_RECONSTRUCTION_CARTESIAN :
					case VISCOUS_RECONSTRUCTION_DIAMOND_PATH :
						{
							int stencil_flag = DIAMONDPATH_RIGHT_TRIANGLE;
							Xu = SolnBlk.Grid.Node[i  ][j+1].X;
							Xd = SolnBlk.Grid.Node[i+1][j+1].X; 
							DiamondPath_Find_dWdX(dWdx, dWdy,
									BlankVector, BlankState,
									Xd, Wd,
									Xr, Wr,
									Xu, Wu,
									stencil_flag);
						}
						break;

					case VISCOUS_RECONSTRUCTION_HYBRID :
						dWdxr = SolnBlk.dWdx[i][j+1]; dWdyr = SolnBlk.dWdy[i][j+1]; 

						Xl = Xface;
						Wl = Wface;
						dWdxl = dWdxr; dWdyl = dWdyr; // Hmm.

						Hybrid_Find_dWdX(dWdx, dWdy,
								Xl, Wl, dWdxl, dWdyl,
								Xr, Wr, dWdxr, dWdyr,
								norm_dir);
						break;
					case VISCOUS_RECONSTRUCTION_ARITHMETIC_AVERAGE :
						dWdx = SolnBlk.dWdx[i][j+1]; dWdy = SolnBlk.dWdy[i][j+1];
						break;
				}

			} else if (j == SolnBlk.JCu && 
					(SolnBlk.Grid.BCtypeN[i] == BC_WALL_VISCOUS_HEATFLUX ||
					 SolnBlk.Grid.BCtypeN[i] == BC_WALL_VISCOUS_ISOTHERMAL ||
					 SolnBlk.Grid.BCtypeN[i] == BC_MOVING_WALL_HEATFLUX ||
					 SolnBlk.Grid.BCtypeN[i] == BC_MOVING_WALL_ISOTHERMAL ||
					 SolnBlk.Grid.BCtypeN[i] == BC_BURNING_SURFACE)) {

				// NORTH face of cell (i,j) is a normal boundary.

				Xl = SolnBlk.Grid.Cell[i][j].Xc; Wl = SolnBlk.W[i][j];

				if (SolnBlk.Grid.BCtypeN[i] == BC_BURNING_SURFACE) {
					Wu = BurningSurface(Wl, norm_dir);
					Wd = Wu;
				} else {

					Wu = HALF*(Wl+SolnBlk.W[i-1][j]);
					Wd = HALF*(Wl+SolnBlk.W[i+1][j]);

					switch (SolnBlk.Grid.BCtypeN[i]) {
						case BC_WALL_VISCOUS_HEATFLUX:
							Wu.v.zero(); Wd.v.zero();
							break;
						case BC_WALL_VISCOUS_ISOTHERMAL:
							Wu.v.zero(); Wd.v.zero();
							switch (HighTemp2D_cState::eos_type) {
								case EOS_TGAS:
									Wu.rho = Tgas_rho(Wu.p, SolnBlk.Twall);
									Wd.rho = Tgas_rho(Wd.p, SolnBlk.Twall);
									break;
								case EOS_IDEAL:
									Wu.rho = Wu.p / HighTemp2D_cState::R / SolnBlk.Twall;
									Wd.rho = Wd.p / HighTemp2D_cState::R / SolnBlk.Twall;
									break;
							}
							break;
						case BC_MOVING_WALL_HEATFLUX:
							Wu.v = SolnBlk.Vwall; Wd.v = SolnBlk.Vwall;
							break;
						case BC_MOVING_WALL_ISOTHERMAL:
							Wu.v = SolnBlk.Vwall; Wd.v = SolnBlk.Vwall; 
							switch (HighTemp2D_cState::eos_type) {
								case EOS_TGAS:
									Wu.rho = Tgas_rho(Wu.p, SolnBlk.Twall);
									Wd.rho = Tgas_rho(Wd.p, SolnBlk.Twall);
									break;
								case EOS_IDEAL:
									Wu.rho = Wu.p / HighTemp2D_cState::R / SolnBlk.Twall;
									Wd.rho = Wd.p / HighTemp2D_cState::R / SolnBlk.Twall;
									break;
							}
							break;
					}
				}

				Wface = HALF * (Wu + Wd);

				switch (IP.i_Viscous_Reconstruction) {
					case VISCOUS_RECONSTRUCTION_CARTESIAN :
					case VISCOUS_RECONSTRUCTION_DIAMOND_PATH :
						{
							int stencil_flag = DIAMONDPATH_LEFT_TRIANGLE;
							Xu = SolnBlk.Grid.Node[i  ][j+1].X;
							Xd = SolnBlk.Grid.Node[i+1][j+1].X; 
							DiamondPath_Find_dWdX(dWdx, dWdy,
									Xl, Wl,
									Xd, Wd,
									BlankVector, BlankState,
									Xu, Wu,
									stencil_flag);
						}
						break;
					case VISCOUS_RECONSTRUCTION_HYBRID :
						dWdxl = SolnBlk.dWdx[i][j]; dWdyl = SolnBlk.dWdy[i][j]; 

						Xr = Xface;
						Wr = Wface;
						dWdxr = dWdxl; dWdyr = dWdyl;

						Hybrid_Find_dWdX(dWdx, dWdy,
								Xl, Wl, dWdxl, dWdyl,
								Xr, Wr, dWdxr, dWdyr,
								norm_dir);
						break;
					case VISCOUS_RECONSTRUCTION_ARITHMETIC_AVERAGE :
						dWdx = SolnBlk.dWdx[i][j]; dWdy = SolnBlk.dWdy[i][j];
						break;
				}
			} else {

				// NORTH face is either a normal cell or possibly a non-viscous
				// boundary condition.
				Xl = SolnBlk.Grid.Cell[i][j  ].Xc; Wl = SolnBlk.W[i][j  ];
				Xr = SolnBlk.Grid.Cell[i][j+1].Xc; Wr = SolnBlk.W[i][j+1];

				if (IP.i_Viscous_Reconstruction == VISCOUS_RECONSTRUCTION_CARTESIAN ||
						IP.i_Viscous_Reconstruction == VISCOUS_RECONSTRUCTION_DIAMOND_PATH) {

					int stencil_flag = DIAMONDPATH_QUADRILATERAL_RECONSTRUCTION;
					Xu = SolnBlk.Grid.Node[i  ][j+1].X; 
					Xd = SolnBlk.Grid.Node[i+1][j+1].X; 

					// These two (WnNW and WnNE) are what make diamond path so expensive.
					Wu = SolnBlk.WnNW(i,j); 
					Wd = SolnBlk.WnNE(i,j);

					DiamondPath_Find_dWdX(dWdx, dWdy,
							Xl, Wl, Xd, Wd, Xr, Wr, Xu, Wu,
							stencil_flag);

					// Find Wface:
					error_flag = Bilinear_Interpolation_ZY(Wl, Xl, Wu, Xu, Wr, Xr, Wd, Xd,
							Xface, Wface);
					// and if error_flag?

				} else {
					Wface = HALF*(Wl + Wr);

					dWdxl = SolnBlk.dWdx[i][j]; dWdxr = SolnBlk.dWdx[i][j+1];
					dWdyl = SolnBlk.dWdy[i][j]; dWdyr = SolnBlk.dWdy[i][j+1];

					if (IP.i_Viscous_Reconstruction == VISCOUS_RECONSTRUCTION_HYBRID) {
						Hybrid_Find_dWdX(dWdx, dWdy,
								Xl, Wl, dWdxl, dWdyl,
								Xr, Wr, dWdxr, dWdyr,
								norm_dir);
					} else if (IP.i_Viscous_Reconstruction == VISCOUS_RECONSTRUCTION_ARITHMETIC_AVERAGE) {
						dWdx = HALF*(dWdxl + dWdxr);
						dWdy = HALF*(dWdyl + dWdyr);
					}
				}
			}
			// Compute the j-direction viscous flux.
			Flux -= ViscousFlux_n(Xface, Wface, dWdx, dWdy, norm_dir, SolnBlk.Axisymmetric);

			if (IP.Solver_Type == IMPLICIT && SolnBlk.face_grad_arrays_allocated) {
				SolnBlk.dWdx_faceN[i][j] = dWdx;
				SolnBlk.dWdy_faceN[i][j] = dWdy;
			}
		} // if (SolnBlk.Flow_Type)

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

    SolnBlk.dUdt[i][SolnBlk.JCl-1][0] = HighTemp2D_U_VACUUM;
    SolnBlk.dUdt[i][SolnBlk.JCu+1][0] = HighTemp2D_U_VACUUM;

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
int dUdt_Multistage_Explicit(HighTemp2D_Quad_Block &SolnBlk,
                             const int i_stage,
			     HighTemp2D_Input_Parameters &IP) {

  int error_flag, k_residual;
  double omega;
 
  HighTemp2D_pState Wl, Wr, Wu, Wd, Wface;
  HighTemp2D_pState dWdx, dWdy, dWdxl, dWdyl, dWdxr, dWdyr;
  HighTemp2D_cState Flux;

  Vector2D dX, Xl, Xr, Xu, Xd;

  HighTemp2D_cState HighTemp2D_U_VACUUM; HighTemp2D_U_VACUUM.Vacuum();
  HighTemp2D_pState HighTemp2D_W_STDATM; HighTemp2D_W_STDATM.Standard_Atmosphere();

  HighTemp2D_pState BlankState; BlankState.Vacuum();
  Vector2D BlankVector; BlankVector.zero();
  
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

  if (IP.i_Limiter == LIMITER_ZERO) {
    // Zeroing phi does not need to be done on every call to
    // dUdt_Residual_Evaluation(). But the limiter type can change during
    // the simulation so to be safe we simply zero phi everytime.
    for (int i = 0; i < SolnBlk.NCi; i++) {
	for (int j = 0; j < SolnBlk.NCj; j++) {
	   SolnBlk.phi[i][j].Vacuum();
	}
    }
  } else {
    // Perform the linear reconstruction within each cell of the
    // computational grid for this stage.
    switch(IP.i_Reconstruction) {
      case RECONSTRUCTION_GREEN_GAUSS :
	Linear_Reconstruction_GreenGauss(SolnBlk,IP.i_Limiter);
	break;
      case RECONSTRUCTION_LINEAR_LEAST_SQUARES :
	Linear_Reconstruction_LeastSquares(SolnBlk,IP.i_Limiter);
	break;
    //case RECONSTRUCTION_QUADRATIC_LEAST_SQUARES :
    //  Quadratic_Reconstruction_LeastSquares(SolnBlk,IP.i_Limiter);
    //  break;
      default:
	Linear_Reconstruction_LeastSquares(SolnBlk,IP.i_Limiter);
	break;
    }
  } /* endif */

  // Evaluate the time rate of change of the solution (i.e., the
  // solution residuals) using a second-order limited upwind scheme
  // with a variety of flux functions.

  // Add i-direction (zeta-direction) fluxes.
  for (int j = SolnBlk.JCl-1; j <= SolnBlk.JCu+1; j++) {
    if (i_stage == 1) {
      SolnBlk.Uo[SolnBlk.ICl-1][j] = SolnBlk.U[SolnBlk.ICl-1][j];
      SolnBlk.dUdt[SolnBlk.ICl-1][j][k_residual] = HighTemp2D_U_VACUUM;
    } else {
      SolnBlk.dUdt[SolnBlk.ICl-1][j][k_residual] = HighTemp2D_U_VACUUM;
    }

    for (int i = SolnBlk.ICl-1; i <= SolnBlk.ICu; i++) {
      if (i_stage == 1) {
	SolnBlk.Uo[i+1][j] = SolnBlk.U[i+1][j];
	SolnBlk.dUdt[i+1][j][k_residual] = HighTemp2D_U_VACUUM;
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
	    SolnBlk.dUdt[i+1][j][k_residual] = HighTemp2D_U_VACUUM;
	  }
	  break;
	case TIME_STEPPING_MULTISTAGE_OPTIMAL_SMOOTHING :
	  SolnBlk.dUdt[i+1][j][k_residual] = HighTemp2D_U_VACUUM;
	  break;
	default:
	  SolnBlk.dUdt[i+1][j][k_residual] = HighTemp2D_U_VACUUM;
	  break;
	};
      }

      if (j >= SolnBlk.JCl && j <= SolnBlk.JCu) {

	// Evaluate the cell interface i-direction INVISCID fluxes.
	if (i == SolnBlk.ICl-1 && 
	    (SolnBlk.Grid.BCtypeW[j] == BC_REFLECTION ||
	     SolnBlk.Grid.BCtypeW[j] == BC_BURNING_SURFACE ||
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
	    Wl = RinglebFlow(Wr,SolnBlk.Grid.xfaceW(i+1,j));
	  } else {
	    // WEST face of cell (i+1,j) is a CHARACTERISTIC boundary.
	    Wl = BC_Characteristic_Pressure(Wr,SolnBlk.WoW[j],SolnBlk.Grid.nfaceW(i+1,j));
	  }
	} else if (i == SolnBlk.ICu &&
		   (SolnBlk.Grid.BCtypeE[j] == BC_REFLECTION ||
		    SolnBlk.Grid.BCtypeE[j] == BC_BURNING_SURFACE ||
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
	    Wr = RinglebFlow(Wl,SolnBlk.Grid.xfaceE(i,j));
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
	case FLUX_FUNCTION_GLAISTER :
	  Flux = FluxGlaister_n(Wl,Wr,SolnBlk.Grid.nfaceE(i,j));
	  break;
	case FLUX_FUNCTION_RUSANOV :
	  Flux = FluxRusanov_n(Wl,Wr,SolnBlk.Grid.nfaceE(i,j));
	  break;
	case FLUX_FUNCTION_HLLE :
	  Flux = FluxHLLE_n(Wl,Wr,SolnBlk.Grid.nfaceE(i,j));
	  break;
	case FLUX_FUNCTION_GHLLE :
	  Flux = FluxGHLLE_n(Wl,Wr,SolnBlk.Grid.nfaceE(i,j));
	  break; 
	case FLUX_FUNCTION_HLLL :
	  Flux = FluxHLLL_n(Wl,Wr,SolnBlk.Grid.nfaceE(i,j));
	  break;
	case FLUX_FUNCTION_GHLLL :
	  Flux = FluxGHLLL_n(Wl,Wr,SolnBlk.Grid.nfaceE(i,j));
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
	case FLUX_FUNCTION_AUSM_PLUS_UP :
	  Flux = FluxAUSMplusUP_n(Wl,Wr,SolnBlk.Grid.nfaceE(i,j));
	  break;
	default:
	  //Flux = FluxRoe_n(Wl,Wr,SolnBlk.Grid.nfaceE(i,j));
	  Flux = FluxGlaister_n(Wl,Wr,SolnBlk.Grid.nfaceE(i,j));
	  break;
	};

 	// Compute the cell centred stress tensor and heat flux vector if required.
 	if (SolnBlk.Flow_Type == FLOWTYPE_TURBULENT_RANS_K_OMEGA ||
	    (SolnBlk.Flow_Type && SolnBlk.Axisymmetric)) {
 	  SolnBlk.W[i][j].ComputeViscousTerms(SolnBlk.dWdx[i][j],
 					      SolnBlk.dWdy[i][j],
 					      SolnBlk.Grid.Cell[i][j].Xc,
 					      SolnBlk.Axisymmetric);
 	  SolnBlk.U[i][j].tau = SolnBlk.W[i][j].tau;
 	  SolnBlk.U[i][j].q = SolnBlk.W[i][j].q;
 	}

	if (SolnBlk.Flow_Type) {

		// Determine the i-direction viscous flux.

		// See the discussion at the start of code for the 
		// i-direction viscous flux in "dUdt_Residual_Evaluation()"
		// for more on the evaluation of the viscous flux.

		Vector2D Xface = SolnBlk.Grid.xfaceE(i, j);
		Vector2D norm_dir = SolnBlk.Grid.nfaceE(i, j);

		if (i == SolnBlk.ICl-1 && 
				(SolnBlk.Grid.BCtypeW[j] == BC_WALL_VISCOUS_HEATFLUX ||
				 SolnBlk.Grid.BCtypeW[j] == BC_WALL_VISCOUS_ISOTHERMAL ||
				 SolnBlk.Grid.BCtypeW[j] == BC_MOVING_WALL_HEATFLUX ||
				 SolnBlk.Grid.BCtypeW[j] == BC_MOVING_WALL_ISOTHERMAL ||
				 SolnBlk.Grid.BCtypeW[j] == BC_BURNING_SURFACE)) {

			// WEST face of cell (i+1,j) is a normal boundary.

			Xr = SolnBlk.Grid.Cell[i+1][j].Xc; Wr = SolnBlk.W[i+1][j];

			if (SolnBlk.Grid.BCtypeW[j] == BC_BURNING_SURFACE) {
				Wu = BurningSurface(Wr, norm_dir);
				Wd = Wu;
			} else {

				Wu = HALF*(Wr+SolnBlk.W[i+1][j+1]);
				Wd = HALF*(Wr+SolnBlk.W[i+1][j-1]);

				switch (SolnBlk.Grid.BCtypeW[j]) {
					case BC_WALL_VISCOUS_HEATFLUX:
						Wu.v.zero(); Wd.v.zero();
						break;
					case BC_WALL_VISCOUS_ISOTHERMAL:
						Wu.v.zero(); Wd.v.zero();
						switch (HighTemp2D_cState::eos_type) {
							case EOS_TGAS:
								Wu.rho = Tgas_rho(Wu.p, SolnBlk.Twall);
								Wd.rho = Tgas_rho(Wd.p, SolnBlk.Twall);
								break;
							case EOS_IDEAL:
								Wu.rho = Wu.p / HighTemp2D_cState::R / SolnBlk.Twall;
								Wd.rho = Wd.p / HighTemp2D_cState::R / SolnBlk.Twall;
								break;
						}
						break;
					case BC_MOVING_WALL_HEATFLUX:
						Wu.v = SolnBlk.Vwall; Wd.v = SolnBlk.Vwall;
						break;
					case BC_MOVING_WALL_ISOTHERMAL:
						Wu.v = SolnBlk.Vwall; Wd.v = SolnBlk.Vwall; 
						switch (HighTemp2D_cState::eos_type) {
							case EOS_TGAS:
								Wu.rho = Tgas_rho(Wu.p, SolnBlk.Twall);
								Wd.rho = Tgas_rho(Wd.p, SolnBlk.Twall);
								break;
							case EOS_IDEAL:
								Wu.rho = Wu.p / HighTemp2D_cState::R / SolnBlk.Twall;
								Wd.rho = Wd.p / HighTemp2D_cState::R / SolnBlk.Twall;
								break;
						}
						break;
				}
			}

			Wface = HALF * (Wu + Wd);

			switch (IP.i_Viscous_Reconstruction) {
				case VISCOUS_RECONSTRUCTION_CARTESIAN :
				case VISCOUS_RECONSTRUCTION_DIAMOND_PATH :
					{
						int stencil_flag = DIAMONDPATH_RIGHT_TRIANGLE;
						Xu = SolnBlk.Grid.Node[i+1][j+1].X;
						Xd = SolnBlk.Grid.Node[i+1][j  ].X; 
						DiamondPath_Find_dWdX(dWdx, dWdy,
								BlankVector, BlankState, 
								Xd, Wd,
								Xr, Wr,
								Xu, Wu,
								stencil_flag);
					}
					break;
				case VISCOUS_RECONSTRUCTION_HYBRID :
					dWdxr = SolnBlk.dWdx[i+1][j]; dWdyr = SolnBlk.dWdy[i+1][j];

					Xl = Xface;
					Wl = Wface;
					dWdxl = dWdxr; dWdyl = dWdyr; // Hmm. 

					Hybrid_Find_dWdX(dWdx, dWdy, 
							Xl, Wl, dWdxl, dWdyl,
							Xr, Wr, dWdxr, dWdyr,
							norm_dir);
					break;
				case VISCOUS_RECONSTRUCTION_ARITHMETIC_AVERAGE :
					dWdx = SolnBlk.dWdx[i+1][j]; dWdy = SolnBlk.dWdy[i+1][j];
					break;
			}
		} else if (i == SolnBlk.ICu &&
				(SolnBlk.Grid.BCtypeE[j] == BC_WALL_VISCOUS_HEATFLUX ||
				 SolnBlk.Grid.BCtypeE[j] == BC_WALL_VISCOUS_ISOTHERMAL ||
				 SolnBlk.Grid.BCtypeE[j] == BC_MOVING_WALL_HEATFLUX ||
				 SolnBlk.Grid.BCtypeE[j] == BC_MOVING_WALL_ISOTHERMAL ||
				 SolnBlk.Grid.BCtypeE[j] == BC_BURNING_SURFACE)) {

			// EAST face of cell (i,j) is a normal boundary.

			Xl = SolnBlk.Grid.Cell[i][j].Xc; Wl = SolnBlk.W[i][j];

			if (SolnBlk.Grid.BCtypeE[j] == BC_BURNING_SURFACE) {
				Wu = BurningSurface(Wl, norm_dir);
				Wd = Wu;
			} else {

				Wu = HALF*(Wl+SolnBlk.W[i][j+1]);
				Wd = HALF*(Wl+SolnBlk.W[i][j-1]);

				switch (SolnBlk.Grid.BCtypeE[j]) {
					case BC_WALL_VISCOUS_HEATFLUX:
						Wu.v.zero(); Wd.v.zero();
						break;
					case BC_WALL_VISCOUS_ISOTHERMAL:
						Wu.v.zero(); Wd.v.zero();
						switch (HighTemp2D_cState::eos_type) {
							case EOS_TGAS:
								Wu.rho = Tgas_rho(Wu.p, SolnBlk.Twall);
								Wd.rho = Tgas_rho(Wd.p, SolnBlk.Twall);
								break;
							case EOS_IDEAL:
								Wu.rho = Wu.p / HighTemp2D_cState::R / SolnBlk.Twall;
								Wd.rho = Wd.p / HighTemp2D_cState::R / SolnBlk.Twall;
								break;
						}
						break;
					case BC_MOVING_WALL_HEATFLUX:
						Wu.v = SolnBlk.Vwall; Wd.v = SolnBlk.Vwall;
						break;
					case BC_MOVING_WALL_ISOTHERMAL:
						Wu.v = SolnBlk.Vwall; Wd.v = SolnBlk.Vwall;
						switch (HighTemp2D_cState::eos_type) {
							case EOS_TGAS:
								Wu.rho = Tgas_rho(Wu.p, SolnBlk.Twall);
								Wd.rho = Tgas_rho(Wd.p, SolnBlk.Twall);
								break;
							case EOS_IDEAL:
								Wu.rho = Wu.p / HighTemp2D_cState::R / SolnBlk.Twall;
								Wd.rho = Wd.p / HighTemp2D_cState::R / SolnBlk.Twall;
								break;
						}
						break;
				}
			}

			Wface = HALF * (Wu + Wd);

			switch(IP.i_Viscous_Reconstruction) {
				case VISCOUS_RECONSTRUCTION_CARTESIAN :
				case VISCOUS_RECONSTRUCTION_DIAMOND_PATH :
					{
						int stencil_flag = DIAMONDPATH_LEFT_TRIANGLE;
						Xu = SolnBlk.Grid.Node[i+1][j+1].X;
						Xd = SolnBlk.Grid.Node[i+1][j  ].X; 
						DiamondPath_Find_dWdX(dWdx, dWdy,
								Xl, Wl,
								Xd, Wd,
								BlankVector, BlankState,
								Xu, Wu,
								stencil_flag);
					}
					break;
				case VISCOUS_RECONSTRUCTION_HYBRID :
					dWdxl = SolnBlk.dWdx[i][j]; dWdyl = SolnBlk.dWdy[i][j]; 

					Xr = Xface;
					Wr = Wface;
					dWdxr = dWdxl; dWdyr = dWdyl; 

					Hybrid_Find_dWdX(dWdx, dWdy, 
							Xl, Wl, dWdxl, dWdyl,
							Xr, Wr, dWdxr, dWdyr,
							norm_dir);
					break;
				case VISCOUS_RECONSTRUCTION_ARITHMETIC_AVERAGE :
					dWdx = SolnBlk.dWdx[i][j]; dWdy = SolnBlk.dWdy[i][j];
					break;
			}
		} else {
			// EAST face is either a normal cell or possibly 
			// a non-viscous boundary condition.
			Xl = SolnBlk.Grid.Cell[i  ][j].Xc; Wl = SolnBlk.W[i  ][j];
			Xr = SolnBlk.Grid.Cell[i+1][j].Xc; Wr = SolnBlk.W[i+1][j];

			if (IP.i_Viscous_Reconstruction == VISCOUS_RECONSTRUCTION_CARTESIAN ||
					IP.i_Viscous_Reconstruction == VISCOUS_RECONSTRUCTION_DIAMOND_PATH) {

				int stencil_flag = DIAMONDPATH_QUADRILATERAL_RECONSTRUCTION;
				Xu = SolnBlk.Grid.Node[i+1][j+1].X; 
				Xd = SolnBlk.Grid.Node[i+1][j  ].X; 

				// These two (WnNE and WnSE) are what make diamond path so expensive.
				Wu = SolnBlk.WnNE(i,j); 
				Wd = SolnBlk.WnSE(i,j);

				DiamondPath_Find_dWdX(dWdx, dWdy,
						Xl, Wl, Xd, Wd, Xr, Wr, Xu, Wu, 
						stencil_flag);

				// Find Wface:
				error_flag = Bilinear_Interpolation_ZY(Wl, Xl, Wu, Xu, Wr, Xr, Wd, Xd,
						Xface, Wface);
				// and if error_flag?

			} else {
				Wface = HALF*(Wl + Wr);

				dWdxl = SolnBlk.dWdx[i][j]; dWdxr = SolnBlk.dWdx[i+1][j];
				dWdyl = SolnBlk.dWdy[i][j]; dWdyr = SolnBlk.dWdy[i+1][j];

				if (IP.i_Viscous_Reconstruction == VISCOUS_RECONSTRUCTION_HYBRID) {
					Hybrid_Find_dWdX(dWdx, dWdy, 
							Xl, Wl, dWdxl, dWdyl,
							Xr, Wr, dWdxr, dWdyr,
							norm_dir);
				} else if (IP.i_Viscous_Reconstruction == VISCOUS_RECONSTRUCTION_ARITHMETIC_AVERAGE) {
					dWdx = HALF*(dWdxl + dWdxr);
					dWdy = HALF*(dWdyl + dWdyr);
				}
			}
		}
		// Compute the i-direction viscous flux.
		Flux -= ViscousFlux_n(Xface, Wface, dWdx, dWdy, norm_dir, SolnBlk.Axisymmetric);

		// We need to save dWd[xy] for use in the implicit
		// preconditioner. Question: would we ever call
		// "dUdt_Multistage_Explicit()" instead of
		// "dUdt_Residual_Evaluation()" in the implicit code? 
		// I say no. Therefore no need to save dWd[xy] here. 
		// -- Alistair Wood Monday, December 4th, 2006 

	} // if (SolnBlk.Flow_Type) 

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
								 SolnBlk.dWdy[i][j]);
 	}

 	// Include turbulent production and destruction source term if required.
 	if (SolnBlk.Flow_Type == FLOWTYPE_TURBULENT_RANS_K_OMEGA) {
 	  SolnBlk.dUdt[i][j][k_residual] += (IP.CFL_Number*SolnBlk.dt[i][j])*
 	                                    SolnBlk.W[i][j].St(SolnBlk.Grid.Cell[i][j].Xc,
							       SolnBlk.dWdx[i][j],
							       SolnBlk.dWdy[i][j],
 							       SolnBlk.Axisymmetric);
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
      SolnBlk.dUdt[SolnBlk.ICl-1][j][k_residual] = HighTemp2D_U_VACUUM;
      SolnBlk.dUdt[SolnBlk.ICu+1][j][k_residual] = HighTemp2D_U_VACUUM;
    }

  }

  // Add j-direction (eta-direction) fluxes.
  for (int i = SolnBlk.ICl; i <= SolnBlk.ICu; i++) {
    for (int j = SolnBlk.JCl-1; j <= SolnBlk.JCu; j++) {

      // Evaluate the cell interface j-direction INVISCID fluxes.
      if (j == SolnBlk.JCl-1 && 
	  (SolnBlk.Grid.BCtypeS[i] == BC_REFLECTION ||
	   SolnBlk.Grid.BCtypeS[i] == BC_BURNING_SURFACE ||
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
	}  else if (SolnBlk.Grid.BCtypeS[i] == BC_BURNING_SURFACE) {
	  // SOUTH face of cell (i,j+1) is a BURNING_SURFACE boundary.
	  Wl = BurningSurface(Wr,SolnBlk.Grid.nfaceS(i,j+1));
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
	  Wl = RinglebFlow(Wr,SolnBlk.Grid.xfaceS(i,j+1));
	  Wl = BC_Characteristic_Pressure(Wr,Wl,SolnBlk.Grid.nfaceS(i,j+1));
	}
	  else if (SolnBlk.Grid.BCtypeS[i] == BC_CHARACTERISTIC) {
	  // SOUTH face of cell (i,j+1) is a CHARACTERISTIC boundary.
	  Wl = BC_Characteristic_Pressure(Wr,SolnBlk.WoS[i],SolnBlk.Grid.nfaceS(i,j+1));
	  } 
      } else if (j == SolnBlk.JCu && 
		 (SolnBlk.Grid.BCtypeN[i] == BC_REFLECTION ||
		  SolnBlk.Grid.BCtypeN[i] == BC_BURNING_SURFACE ||
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
	  Wr = RinglebFlow(Wl,SolnBlk.Grid.xfaceN(i,j));
	  Wr = BC_Characteristic_Pressure(Wl,Wr,SolnBlk.Grid.nfaceN(i,j));
	} 
	 else {
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
      case FLUX_FUNCTION_GLAISTER :
	Flux = FluxGlaister_n(Wl,Wr,SolnBlk.Grid.nfaceN(i,j));
	break;
      case FLUX_FUNCTION_RUSANOV :
	Flux = FluxRusanov_n(Wl,Wr,SolnBlk.Grid.nfaceN(i,j));
	break;
      case FLUX_FUNCTION_HLLE :
	Flux = FluxHLLE_n(Wl,Wr,SolnBlk.Grid.nfaceN(i,j));
	break;
      case FLUX_FUNCTION_GHLLE :
	Flux = FluxGHLLE_n(Wl,Wr,SolnBlk.Grid.nfaceN(i,j));
	break;
      case FLUX_FUNCTION_HLLL :
	Flux = FluxHLLL_n(Wl,Wr,SolnBlk.Grid.nfaceN(i,j));
	break;
      case FLUX_FUNCTION_GHLLL :
	Flux = FluxGHLLL_n(Wl,Wr,SolnBlk.Grid.nfaceN(i,j));
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
      case FLUX_FUNCTION_AUSM_PLUS_UP :
	Flux = FluxAUSMplusUP_n(Wl,Wr,SolnBlk.Grid.nfaceN(i,j));
	break;
      default:
	//Flux = FluxRoe_n(Wl,Wr,SolnBlk.Grid.nfaceN(i,j));
        Flux = FluxGlaister_n(Wl,Wr,SolnBlk.Grid.nfaceN(i,j));
	break;
      };

      if (SolnBlk.Flow_Type) {

				// Determine the j-direction viscous flux.

				// Please rotate your head 90 degrees to the left before
				// attempting to read this code. For example, Xr is the
				// position of the cell to the north of the face in
				// question. This is necessary, for example, for
				// "DiamondPath_Find_dWdX()", which only deals with "left" 
				// and "right" triangles.

				// See the discussion at the start of code for the 
				// i-direction viscous flux in "dUdt_Residual_Evaluation()"
				// for more on the evaluation of the viscous flux.

				Vector2D Xface = SolnBlk.Grid.xfaceN(i, j);
				Vector2D norm_dir = SolnBlk.Grid.nfaceN(i, j);

				if (j == SolnBlk.JCl-1 && 
						(SolnBlk.Grid.BCtypeS[i] == BC_WALL_VISCOUS_HEATFLUX ||
						 SolnBlk.Grid.BCtypeS[i] == BC_WALL_VISCOUS_ISOTHERMAL ||
						 SolnBlk.Grid.BCtypeS[i] == BC_MOVING_WALL_HEATFLUX ||
						 SolnBlk.Grid.BCtypeS[i] == BC_MOVING_WALL_ISOTHERMAL ||
						 SolnBlk.Grid.BCtypeS[i] == BC_BURNING_SURFACE)) {

					// SOUTH face of cell (i,j+1) is a normal boundary.

					Xr = SolnBlk.Grid.Cell[i][j+1].Xc; Wr = SolnBlk.W[i][j+1];

					if (SolnBlk.Grid.BCtypeS[i] == BC_BURNING_SURFACE) {
						Wu = BurningSurface(Wr, norm_dir);
						Wd = Wu;
					} else {

						Wu = HALF*(Wr+SolnBlk.W[i-1][j+1]);
						Wd = HALF*(Wr+SolnBlk.W[i+1][j+1]);

						switch (SolnBlk.Grid.BCtypeS[i]) {
							case BC_WALL_VISCOUS_HEATFLUX:
								Wu.v.zero(); Wd.v.zero();
								break;
							case BC_WALL_VISCOUS_ISOTHERMAL:
								Wu.v.zero(); Wd.v.zero();
								switch (HighTemp2D_cState::eos_type) {
									case EOS_TGAS:
										Wu.rho = Tgas_rho(Wu.p, SolnBlk.Twall);
										Wd.rho = Tgas_rho(Wd.p, SolnBlk.Twall);
										break;
									case EOS_IDEAL:
										Wu.rho = Wu.p / HighTemp2D_cState::R / SolnBlk.Twall;
										Wd.rho = Wd.p / HighTemp2D_cState::R / SolnBlk.Twall;
										break;
								}
								break;
							case BC_MOVING_WALL_HEATFLUX:
								Wu.v = SolnBlk.Vwall; Wd.v = SolnBlk.Vwall;
								break;
							case BC_MOVING_WALL_ISOTHERMAL:
								Wu.v = SolnBlk.Vwall; Wd.v = SolnBlk.Vwall; 
								switch (HighTemp2D_cState::eos_type) {
									case EOS_TGAS:
										Wu.rho = Tgas_rho(Wu.p, SolnBlk.Twall);
										Wd.rho = Tgas_rho(Wd.p, SolnBlk.Twall);
										break;
									case EOS_IDEAL:
										Wu.rho = Wu.p / HighTemp2D_cState::R / SolnBlk.Twall;
										Wd.rho = Wd.p / HighTemp2D_cState::R / SolnBlk.Twall;
										break;
								}
								break;
						} 
					}

					Wface = HALF * (Wu + Wd);

					switch (IP.i_Viscous_Reconstruction) {
						case VISCOUS_RECONSTRUCTION_CARTESIAN :
						case VISCOUS_RECONSTRUCTION_DIAMOND_PATH :
							{
								int stencil_flag = DIAMONDPATH_RIGHT_TRIANGLE;
								Xu = SolnBlk.Grid.Node[i  ][j+1].X;
								Xd = SolnBlk.Grid.Node[i+1][j+1].X; 
								DiamondPath_Find_dWdX(dWdx, dWdy,
										BlankVector, BlankState,
										Xd, Wd,
										Xr, Wr,
										Xu, Wu,
										stencil_flag);
							}
							break;

						case VISCOUS_RECONSTRUCTION_HYBRID :
							dWdxr = SolnBlk.dWdx[i][j+1]; dWdyr = SolnBlk.dWdy[i][j+1]; 

							Xl = Xface;
							Wl = Wface;
							dWdxl = dWdxr; dWdyl = dWdyr; // Hmm.

							Hybrid_Find_dWdX(dWdx, dWdy,
									Xl, Wl, dWdxl, dWdyl,
									Xr, Wr, dWdxr, dWdyr,
									norm_dir);
							break;
						case VISCOUS_RECONSTRUCTION_ARITHMETIC_AVERAGE :
							dWdx = SolnBlk.dWdx[i][j+1]; dWdy = SolnBlk.dWdy[i][j+1];
							break;
					}

				} else if (j == SolnBlk.JCu && 
						(SolnBlk.Grid.BCtypeN[i] == BC_WALL_VISCOUS_HEATFLUX ||
						 SolnBlk.Grid.BCtypeN[i] == BC_WALL_VISCOUS_ISOTHERMAL ||
						 SolnBlk.Grid.BCtypeN[i] == BC_MOVING_WALL_HEATFLUX ||
						 SolnBlk.Grid.BCtypeN[i] == BC_MOVING_WALL_ISOTHERMAL ||
						 SolnBlk.Grid.BCtypeN[i] == BC_BURNING_SURFACE)) {

					// NORTH face of cell (i,j) is a normal boundary.

					Xl = SolnBlk.Grid.Cell[i][j].Xc; Wl = SolnBlk.W[i][j];

					if (SolnBlk.Grid.BCtypeN[i] == BC_BURNING_SURFACE) {
						Wu = BurningSurface(Wl, norm_dir);
						Wd = Wu;
					} else {

						Wu = HALF*(Wl+SolnBlk.W[i-1][j]);
						Wd = HALF*(Wl+SolnBlk.W[i+1][j]);

						switch (SolnBlk.Grid.BCtypeN[i]) {
							case BC_WALL_VISCOUS_HEATFLUX:
								Wu.v.zero(); Wd.v.zero();
								break;
							case BC_WALL_VISCOUS_ISOTHERMAL:
								Wu.v.zero(); Wd.v.zero();
								switch (HighTemp2D_cState::eos_type) {
									case EOS_TGAS:
										Wu.rho = Tgas_rho(Wu.p, SolnBlk.Twall);
										Wd.rho = Tgas_rho(Wd.p, SolnBlk.Twall);
										break;
									case EOS_IDEAL:
										Wu.rho = Wu.p / HighTemp2D_cState::R / SolnBlk.Twall;
										Wd.rho = Wd.p / HighTemp2D_cState::R / SolnBlk.Twall;
										break;
								}
								break;
							case BC_MOVING_WALL_HEATFLUX:
								Wu.v = SolnBlk.Vwall; Wd.v = SolnBlk.Vwall;
								break;
							case BC_MOVING_WALL_ISOTHERMAL:
								Wu.v = SolnBlk.Vwall; Wd.v = SolnBlk.Vwall; 
								switch (HighTemp2D_cState::eos_type) {
									case EOS_TGAS:
										Wu.rho = Tgas_rho(Wu.p, SolnBlk.Twall);
										Wd.rho = Tgas_rho(Wd.p, SolnBlk.Twall);
										break;
									case EOS_IDEAL:
										Wu.rho = Wu.p / HighTemp2D_cState::R / SolnBlk.Twall;
										Wd.rho = Wd.p / HighTemp2D_cState::R / SolnBlk.Twall;
										break;
								}
								break;
						}
					}

					Wface = HALF * (Wu + Wd);

					switch (IP.i_Viscous_Reconstruction) {
						case VISCOUS_RECONSTRUCTION_CARTESIAN :
						case VISCOUS_RECONSTRUCTION_DIAMOND_PATH :
							{
								int stencil_flag = DIAMONDPATH_LEFT_TRIANGLE;
								Xu = SolnBlk.Grid.Node[i  ][j+1].X;
								Xd = SolnBlk.Grid.Node[i+1][j+1].X; 
								DiamondPath_Find_dWdX(dWdx, dWdy,
										Xl, Wl,
										Xd, Wd,
										BlankVector, BlankState,
										Xu, Wu,
										stencil_flag);
							}
							break;
						case VISCOUS_RECONSTRUCTION_HYBRID :
							dWdxl = SolnBlk.dWdx[i][j]; dWdyl = SolnBlk.dWdy[i][j]; 

							Xr = Xface;
							Wr = Wface;
							dWdxr = dWdxl; dWdyr = dWdyl;

							Hybrid_Find_dWdX(dWdx, dWdy,
									Xl, Wl, dWdxl, dWdyl,
									Xr, Wr, dWdxr, dWdyr,
									norm_dir);
							break;
						case VISCOUS_RECONSTRUCTION_ARITHMETIC_AVERAGE :
							dWdx = SolnBlk.dWdx[i][j]; dWdy = SolnBlk.dWdy[i][j];
							break;
					}
				} else {

					// NORTH face is either a normal cell or possibly a non-viscous
					// boundary condition.
					Xl = SolnBlk.Grid.Cell[i][j  ].Xc; Wl = SolnBlk.W[i][j  ];
					Xr = SolnBlk.Grid.Cell[i][j+1].Xc; Wr = SolnBlk.W[i][j+1];

					if (IP.i_Viscous_Reconstruction == VISCOUS_RECONSTRUCTION_CARTESIAN ||
							IP.i_Viscous_Reconstruction == VISCOUS_RECONSTRUCTION_DIAMOND_PATH) {

						int stencil_flag = DIAMONDPATH_QUADRILATERAL_RECONSTRUCTION;
						Xu = SolnBlk.Grid.Node[i  ][j+1].X; 
						Xd = SolnBlk.Grid.Node[i+1][j+1].X; 

						// These two (WnNW and WnNE) are what make diamond path so expensive.
						Wu = SolnBlk.WnNW(i,j); 
						Wd = SolnBlk.WnNE(i,j);

						DiamondPath_Find_dWdX(dWdx, dWdy,
								Xl, Wl, Xd, Wd, Xr, Wr, Xu, Wu,
								stencil_flag);

						// Find Wface:
						error_flag = Bilinear_Interpolation_ZY(Wl, Xl, Wu, Xu, Wr, Xr, Wd, Xd,
								Xface, Wface);
						// and if error_flag?

					} else {
						Wface = HALF*(Wl + Wr);

						dWdxl = SolnBlk.dWdx[i][j]; dWdxr = SolnBlk.dWdx[i][j+1];
						dWdyl = SolnBlk.dWdy[i][j]; dWdyr = SolnBlk.dWdy[i][j+1];

						if (IP.i_Viscous_Reconstruction == VISCOUS_RECONSTRUCTION_HYBRID) {
							Hybrid_Find_dWdX(dWdx, dWdy,
									Xl, Wl, dWdxl, dWdyl,
									Xr, Wr, dWdxr, dWdyr,
									norm_dir);
						} else if (IP.i_Viscous_Reconstruction == VISCOUS_RECONSTRUCTION_ARITHMETIC_AVERAGE) {
							dWdx = HALF*(dWdxl + dWdxr);
							dWdy = HALF*(dWdyl + dWdyr);
						}
					}
				}
				// Compute the j-direction viscous flux.
				Flux -= ViscousFlux_n(Xface, Wface, dWdx, dWdy, norm_dir, SolnBlk.Axisymmetric);

				// Again since this is dUdt_Multistage_Explicit() no need to save dWd[xy].

			} // if (SolnBlk.Flow_Type)

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

    SolnBlk.dUdt[i][SolnBlk.JCl-1][k_residual] = HighTemp2D_U_VACUUM;
    SolnBlk.dUdt[i][SolnBlk.JCu+1][k_residual] = HighTemp2D_U_VACUUM;

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
int Update_Solution_Multistage_Explicit(HighTemp2D_Quad_Block &SolnBlk,
					const int i_stage,
					HighTemp2D_Input_Parameters &IP) {

  int k_residual;
  double omega, residual_reduction_factor;

  // Memory for linear system solver.
  DenseMatrix dRdU(NUM_VAR_HIGHTEMP2D,NUM_VAR_HIGHTEMP2D);
  DenseSystemLinEqs LinSys;

  // Allocate memory for linear system solver.
  LinSys.allocate(NUM_VAR_HIGHTEMP2D);

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

	// Include source Jacobian in the LHS matrix.
	LinSys.A.identity();
	LinSys.A -= (omega*IP.CFL_Number*SolnBlk.dt[i][j])*dRdU;

	// Set the explicit residual as the RHS for the point implicit
	// formulation (already contains the CFL number).
	for (int k = 1; k <= NUM_VAR_HIGHTEMP2D; k++)
	  LinSys.b(k-1) = omega*SolnBlk.dUdt[i][j][k_residual][k];

	// Solve system of equations using LU decomposition Gaussian 
	// elimination procedure.
	LinSys.solve(LU_DECOMPOSITION);

	// Update the conserved solution variables.
	for (int k = 1; k <= NUM_VAR_HIGHTEMP2D; k++)
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
	    for (int k = 1; k <= NUM_VAR_HIGHTEMP2D; k++)
	      LinSys.b(k-1) = omega*SolnBlk.dUdt[i][j][k_residual][k];
	    LinSys.solve(LU_DECOMPOSITION);
	    for (int k = 1; k <= NUM_VAR_HIGHTEMP2D; k++)
	      SolnBlk.U[i][j][k] = SolnBlk.Uo[i][j][k] + LinSys.x(k-1);
	  }
	  if (!SolnBlk.U[i][j].Unphysical_Properties()) break;
	  if (n_residual_reduction == 10) cout << "n_residual_reductions = " << n_residual_reduction
					       << " @ Xc =" << SolnBlk.Grid.Cell[i][j].Xc << endl;
	}
      }

      if (SolnBlk.Flow_Type == FLOWTYPE_TURBULENT_RANS_K_OMEGA) {
	if (SolnBlk.U[i][j].dk < TOLER) SolnBlk.U[i][j].dk = TOLER;
	if (SolnBlk.U[i][j].domega < TOLER) SolnBlk.U[i][j].domega = TOLER;
      }

      // Check for unphysical state properties.
      if (SolnBlk.U[i][j].Unphysical_Properties()) {
	cout << "\n " << CFFC_Name() 
	     << " HighTemp2D ERROR: Negative Density and/or Energy: \n"
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

