/**********************************************************************
 * Electrostatic2DQuadSingleBlock.cc: Single-block versions of        *
 *                                    subroutines for 2D              *
 *                                    electrostatic multi-block       *
 *                                    quadrilateral mesh solution     *
 *                                    classes.                        *
 **********************************************************************/

// Include 2D Electrostatic quadrilateral mesh solution header file.

#ifndef _ELECTROSTATIC2D_QUAD_INCLUDED
#include "Electrostatic2DQuad.h"
#endif // _ELECTROSTATIC2D_QUAD_INCLUDED

/**********************************************************************
 * Electrostatic2D_Quad_Block -- Single Block External Subroutines.   *
 **********************************************************************/

/**********************************************************************
 * Routine: Broadcast_Solution_Block                                  *
 *                                                                    *
 * Broadcast quadrilateral solution block to all processors involved  *
 * in the calculation from the primary processor using the MPI        *
 * broadcast routine.                                                 *
 *                                                                    *
 **********************************************************************/
void Broadcast_Solution_Block(Electrostatic2D_Quad_Block &SolnBlk) {

#ifdef _MPI_VERSION
  int ni, nj, ng, block_allocated, buffer_size;
  double *buffer;

  // Broadcast the number of cells in each direction.
  if (CFDkit_Primary_MPI_Processor()) {
    ni = SolnBlk.NCi;
    nj = SolnBlk.NCj;
    ng = SolnBlk.Nghost;
    if (SolnBlk.U != NULL) {
      block_allocated = 1;
    } else {
      block_allocated = 0;
    }
  }

  MPI::COMM_WORLD.Bcast(&ni,1,MPI::INT,0);
  MPI::COMM_WORLD.Bcast(&nj,1,MPI::INT,0);
  MPI::COMM_WORLD.Bcast(&ng,1,MPI::INT,0);
  MPI::COMM_WORLD.Bcast(&block_allocated,1,MPI::INT,0);

  // On non-primary MPI processors, allocate (re-allocate) memory for 
  // the quadrilateral solution block as necessary.
  if (!CFDkit_Primary_MPI_Processor()) {
    if (SolnBlk.NCi != ni || SolnBlk.NCj != nj || SolnBlk.Nghost != ng) { 
      if (SolnBlk.U != NULL) SolnBlk.deallocate();
      if (block_allocated) SolnBlk.allocate(ni-2*ng,nj-2*ng,ng);
    }
  }

  // Broadcast the axisymmetric/planar flow indicator.
  MPI::COMM_WORLD.Bcast(&(SolnBlk.Axisymmetric),1,MPI::INT,0);

  // Broadcast the grid.
  Broadcast_Quad_Block(SolnBlk.Grid);

  // Broadcast the solution state variables.
  if (block_allocated) {
    ni = (SolnBlk.ICu+SolnBlk.Nghost) - (SolnBlk.ICl-SolnBlk.Nghost) + 1;
    nj = (SolnBlk.JCu+SolnBlk.Nghost) - (SolnBlk.JCl-SolnBlk.Nghost) + 1;
    buffer = new double[NUM_VAR_ELECTROSTATIC2D*ni*nj];

    if (CFDkit_Primary_MPI_Processor()) {
      buffer_size = 0;
      for (int j = SolnBlk.JCl-SolnBlk.Nghost; j <= SolnBlk.JCu+SolnBlk.Nghost; j++) {
	for (int i = SolnBlk.ICl-SolnBlk.Nghost; i <= SolnBlk.ICu+SolnBlk.Nghost; i++) {
	  for (int k = 0; k < NUM_VAR_ELECTROSTATIC2D; k++) {
	    buffer[buffer_size] = SolnBlk.U[i][j][k+1];
	    buffer_size++;
	  }
	}
      }
    }

    buffer_size = NUM_VAR_ELECTROSTATIC2D*ni*nj;
    MPI::COMM_WORLD.Bcast(buffer,buffer_size,MPI::DOUBLE,0);

    if (!CFDkit_Primary_MPI_Processor()) {
      buffer_size = 0;
      for (int j = SolnBlk.JCl-SolnBlk.Nghost; j <= SolnBlk.JCu+SolnBlk.Nghost; j++) {
	for (int i = SolnBlk.ICl-SolnBlk.Nghost; i <= SolnBlk.ICu+SolnBlk.Nghost; i++) {
	  for (int k = 0; k < NUM_VAR_ELECTROSTATIC2D; k++) {
	    SolnBlk.U[i][j][k+1] = buffer[buffer_size];
	    buffer_size++;
	  }
	}
      }
    }

    delete []buffer; buffer = NULL;

    ni = 1;
    nj = (SolnBlk.JCu+SolnBlk.Nghost) - (SolnBlk.JCl-SolnBlk.Nghost) + 1;
    buffer = new double[2*NUM_VAR_ELECTROSTATIC2D*ni*nj];

    if (CFDkit_Primary_MPI_Processor()) {
      buffer_size = 0;
      for (int j = SolnBlk.JCl-SolnBlk.Nghost; j <= SolnBlk.JCu+SolnBlk.Nghost; j++) {
	for (int k = 0; k < NUM_VAR_ELECTROSTATIC2D; k++) {
	  buffer[buffer_size] = SolnBlk.UoW[j][k+1];
	  buffer_size++;
	}
	for (int k = 0; k < NUM_VAR_ELECTROSTATIC2D; k++) {
	  buffer[buffer_size] = SolnBlk.UoE[j][k+1];
	  buffer_size++;
	}
      }
    }

    buffer_size = 2*NUM_VAR_ELECTROSTATIC2D*ni*nj;
    MPI::COMM_WORLD.Bcast(buffer,buffer_size,MPI::DOUBLE,0);

    if (!CFDkit_Primary_MPI_Processor()) {
      buffer_size = 0;
      for (int j = SolnBlk.JCl-SolnBlk.Nghost; j <= SolnBlk.JCu+SolnBlk.Nghost; j++) {
	for (int k = 0; k < NUM_VAR_ELECTROSTATIC2D; k++) {
	  SolnBlk.UoW[j][k+1] = buffer[buffer_size];
	  buffer_size++;
	}
	for (int k = 0; k < NUM_VAR_ELECTROSTATIC2D; k++) {
	  SolnBlk.UoE[j][k+1] = buffer[buffer_size];
	  buffer_size++;
	}
      }
    }

    delete []buffer; buffer = NULL;

    ni = (SolnBlk.ICu+SolnBlk.Nghost) - (SolnBlk.ICl-SolnBlk.Nghost) + 1;
    nj = 1;
    buffer = new double[2*NUM_VAR_ELECTROSTATIC2D*ni*nj];

    if (CFDkit_Primary_MPI_Processor()) {
      buffer_size = 0;
      for (int i = SolnBlk.ICl-SolnBlk.Nghost; i <= SolnBlk.ICu+SolnBlk.Nghost; i++) {
	for (int k = 0; k < NUM_VAR_ELECTROSTATIC2D; k++) {
	  buffer[buffer_size] = SolnBlk.UoS[i][k+1];
	  buffer_size++;
	}
	for (int k = 0; k < NUM_VAR_ELECTROSTATIC2D; k++) {
	  buffer[buffer_size] = SolnBlk.UoN[i][k+1];
	  buffer_size++;
	}
      }
    }

    buffer_size = 2*NUM_VAR_ELECTROSTATIC2D*ni*nj;
    MPI::COMM_WORLD.Bcast(buffer,buffer_size,MPI::DOUBLE,0);

    if (!CFDkit_Primary_MPI_Processor()) {
      buffer_size = 0;
      for (int i = SolnBlk.ICl-SolnBlk.Nghost; i <= SolnBlk.ICu+SolnBlk.Nghost; i++) {
	for (int k = 0; k < NUM_VAR_ELECTROSTATIC2D; k++) {
	  SolnBlk.UoS[i][k+1] = buffer[buffer_size];
	  buffer_size++;
	}
	for (int k = 0; k < NUM_VAR_ELECTROSTATIC2D; k++) {
	  SolnBlk.UoN[i][k+1] = buffer[buffer_size];
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
void Broadcast_Solution_Block(Electrostatic2D_Quad_Block &SolnBlk,
                              MPI::Intracomm &Communicator,
                              const int Source_CPU) {

  int Source_Rank = 0;
  int ni, nj, ng, block_allocated, buffer_size;
  double *buffer;

  // Broadcast the number of cells in each direction.
  if (CFDkit_MPI::This_Processor_Number == Source_CPU) {
    ni = SolnBlk.NCi;
    nj = SolnBlk.NCj;
    ng = SolnBlk.Nghost;
    if (SolnBlk.U != NULL) {
      block_allocated = 1;
    } else {
      block_allocated = 0;
    } 
  }

  Communicator.Bcast(&ni,1,MPI::INT,Source_Rank);
  Communicator.Bcast(&nj,1,MPI::INT,Source_Rank);
  Communicator.Bcast(&ng,1,MPI::INT,Source_Rank);
  Communicator.Bcast(&block_allocated,1,MPI::INT,Source_Rank);

  // On non-source MPI processors, allocate (re-allocate) memory for the
  // quadrilateral solution block as necessary.
  if (!(CFDkit_MPI::This_Processor_Number == Source_CPU)) {
    if (SolnBlk.NCi != ni || SolnBlk.NCj != nj || SolnBlk.Nghost != ng) { 
      if (SolnBlk.U != NULL) SolnBlk.deallocate();
      if (block_allocated) SolnBlk.allocate(ni-2*ng,nj-2*ng,ng); 
    }
  }

  // Broadcast the axisymmetric/planar flow indicator.
  Communicator.Bcast(&(SolnBlk.Axisymmetric),1,MPI::INT,Source_Rank);

  // Broadcast the grid.
  Broadcast_Quad_Block(SolnBlk.Grid,Communicator,Source_CPU);

  // Broadcast the solution state variables.
  if (block_allocated) {
    ni = (SolnBlk.ICu+SolnBlk.Nghost) - (SolnBlk.ICl-SolnBlk.Nghost) + 1;
    nj = (SolnBlk.JCu+SolnBlk.Nghost) - (SolnBlk.JCl-SolnBlk.Nghost) + 1;
    buffer = new double[NUM_VAR_ELECTROSTATIC2D*ni*nj];

    if (CFDkit_MPI::This_Processor_Number == Source_CPU) {
      buffer_size = 0;
      for (int j = SolnBlk.JCl-SolnBlk.Nghost; j <= SolnBlk.JCu+SolnBlk.Nghost; j++) {
	for (int i = SolnBlk.ICl-SolnBlk.Nghost; i <= SolnBlk.ICu+SolnBlk.Nghost; i++) {
	  for (int k = 0; k < NUM_VAR_ELECTROSTATIC2D; k++) {
	    buffer[buffer_size] = SolnBlk.U[i][j][k+1];
	    buffer_size++;
	  }
	}
      }
    }

    buffer_size = NUM_VAR_ELECTROSTATIC2D*ni*nj;
    Communicator.Bcast(buffer,buffer_size,MPI::DOUBLE,Source_Rank);

    if (!(CFDkit_MPI::This_Processor_Number == Source_CPU)) {
      buffer_size = 0;
      for (int j = SolnBlk.JCl-SolnBlk.Nghost; j <= SolnBlk.JCu+SolnBlk.Nghost; j++) {
	for (int i = SolnBlk.ICl-SolnBlk.Nghost; i <= SolnBlk.ICu+SolnBlk.Nghost; i++) {
	  for (int k = 0; k < NUM_VAR_ELECTROSTATIC2D; k++) {
	    SolnBlk.U[i][j][k+1] = buffer[buffer_size];
	    buffer_size++;
	  }
	}
      }
    }

    delete []buffer; buffer = NULL;

    ni = 1;
    nj = (SolnBlk.JCu+SolnBlk.Nghost) - (SolnBlk.JCl-SolnBlk.Nghost) + 1;
    buffer = new double[2*NUM_VAR_ELECTROSTATIC2D*ni*nj];

    if (CFDkit_MPI::This_Processor_Number == Source_CPU) {
      buffer_size = 0;
      for (int j = SolnBlk.JCl-SolnBlk.Nghost; j <= SolnBlk.JCu+SolnBlk.Nghost; j++) {
	for (int k = 0; k < NUM_VAR_ELECTROSTATIC2D; k++) {
	  buffer[buffer_size] = SolnBlk.UoW[j][k+1];
	  buffer_size++;
	}
	for (int k = 0; k < NUM_VAR_ELECTROSTATIC2D; k++) {
	  buffer[buffer_size] = SolnBlk.UoE[j][k+1];
	  buffer_size++;
	}
      }
    }

    buffer_size = 2*NUM_VAR_ELECTROSTATIC2D*ni*nj;
    Communicator.Bcast(buffer,buffer_size,MPI::DOUBLE,Source_Rank);

    if (!(CFDkit_MPI::This_Processor_Number == Source_CPU)) {
      buffer_size = 0;
      for (int j  = SolnBlk.JCl-SolnBlk.Nghost; j <= SolnBlk.JCu+SolnBlk.Nghost; j++) {
	for (int k = 0; k < NUM_VAR_ELECTROSTATIC2D; k++) {
	  SolnBlk.UoW[j][k+1] = buffer[buffer_size];
	  buffer_size++;
	}
	for (int k = 0; k < NUM_VAR_ELECTROSTATIC2D; k++) {
	  SolnBlk.UoE[j][k+1] = buffer[buffer_size];
	  buffer_size++;
	}
      }
    }

    delete []buffer; buffer = NULL;

    ni = (SolnBlk.ICu+SolnBlk.Nghost) - (SolnBlk.ICl-SolnBlk.Nghost) + 1;
    nj = 1;
    buffer = new double[2*NUM_VAR_ELECTROSTATIC2D*ni*nj];

    if (CFDkit_MPI::This_Processor_Number == Source_CPU) {
      buffer_size = 0;
      for (int i  = SolnBlk.ICl-SolnBlk.Nghost; i <= SolnBlk.ICu+SolnBlk.Nghost; i++) {
	for (int k = 0; k < NUM_VAR_ELECTROSTATIC2D; k++) {
	  buffer[buffer_size] = SolnBlk.UoS[i][k+1];
	  buffer_size++;
	}
	for (int k = 0; k < NUM_VAR_ELECTROSTATIC2D; k++) {
	  buffer[buffer_size] = SolnBlk.UoN[i][k+1];
	  buffer_size++;
	}
      }
    }

    buffer_size = 2*NUM_VAR_ELECTROSTATIC2D*ni*nj;
    Communicator.Bcast(buffer,buffer_size,MPI::DOUBLE,Source_Rank);

    if (!(CFDkit_MPI::This_Processor_Number == Source_CPU)) {
      buffer_size = 0;
      for (int i  = SolnBlk.ICl-SolnBlk.Nghost; i <= SolnBlk.ICu+SolnBlk.Nghost; i++) {
	for (int k = 0; k < NUM_VAR_ELECTROSTATIC2D; k++) {
	  SolnBlk.UoS[i][k+1] = buffer[buffer_size];
	  buffer_size++;
	}
	for (int k = 0; k < NUM_VAR_ELECTROSTATIC2D; k++) {
	  SolnBlk.UoN[i][k+1] = buffer[buffer_size];
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
void Copy_Solution_Block(Electrostatic2D_Quad_Block &SolnBlk1,
                         Electrostatic2D_Quad_Block &SolnBlk2) {

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

  // Copy the axisymmetric/planar flow indicator.
  SolnBlk1.Axisymmetric = SolnBlk2.Axisymmetric;

  // Copy the grid of the second solution block to the first solution 
  // block.
  Copy_Quad_Block(SolnBlk1.Grid,SolnBlk2.Grid);

  // Copy the solution information from SolnBlk2 to SolnBlk1.
  if (SolnBlk2.U != NULL) {
    for (int j = SolnBlk1.JCl-SolnBlk1.Nghost; j <= SolnBlk1.JCu+SolnBlk1.Nghost; j++) {
      for (int i = SolnBlk1.ICl-SolnBlk1.Nghost; i <= SolnBlk1.ICu+SolnBlk1.Nghost; i++) {
	SolnBlk1.U[i][j] = SolnBlk2.U[i][j];
	SolnBlk1.U[i][j] = SolnBlk2.U[i][j];
	for (int k = 0; k < NUMBER_OF_RESIDUAL_VECTORS_ELECTROSTATIC2D; k++)
	  SolnBlk1.dUdt[i][j][k] = SolnBlk2.dUdt[i][j][k];
	SolnBlk1.dUdx[i][j] = SolnBlk2.dUdx[i][j];
	SolnBlk1.dUdy[i][j] = SolnBlk2.dUdy[i][j];
	SolnBlk1.Uo[i][j]   = SolnBlk2.Uo[i][j];
	SolnBlk1.dt[i][j]   = SolnBlk2.dt[i][j];
      }
    }

    for (int j = SolnBlk1.JCl-SolnBlk1.Nghost; j <= SolnBlk1.JCu+SolnBlk1.Nghost; j++) {
      SolnBlk1.UoW[j] = SolnBlk2.UoW[j];
      SolnBlk1.UoE[j] = SolnBlk2.UoE[j];
    }

    for (int i = SolnBlk1.ICl-SolnBlk1.Nghost; i <= SolnBlk1.ICu+SolnBlk1.Nghost; i++) {
      SolnBlk1.UoS[i] = SolnBlk2.UoS[i];
      SolnBlk1.UoN[i] = SolnBlk2.UoN[i];
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
int Prolong_Solution_Block(Electrostatic2D_Quad_Block &SolnBlk_Fine,
			   Electrostatic2D_Quad_Block &SolnBlk_Original,
			   const int Sector) {

  int error_flag;
  int i_min, i_max, j_min, j_max, mesh_refinement_permitted;
  int i_fine, j_fine, i_coarse_min, j_coarse_min, coarse_cell_found;
  //double area_total_fine;

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

    // Copy the axisymmetric/planar flow indicator.
    SolnBlk_Fine.Axisymmetric = SolnBlk_Original.Axisymmetric;

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

    // Prolong the solution information using injection (area-weighted
    // average commented out).
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
	SolnBlk_Fine.U[2*(i-i_min)+SolnBlk_Fine.ICl+1][2*(j-j_min)+SolnBlk_Fine.JCl  ] = SolnBlk_Original.U[i][j];
	//= (SolnBlk_Original.Grid.Cell[i][j].A/area_total_fine)*SolnBlk_Original.U[i][j];
	SolnBlk_Fine.U[2*(i-i_min)+SolnBlk_Fine.ICl  ][2*(j-j_min)+SolnBlk_Fine.JCl+1] = SolnBlk_Original.U[i][j];
	//= (SolnBlk_Original.Grid.Cell[i][j].A/area_total_fine)*SolnBlk_Original.U[i][j];
	SolnBlk_Fine.U[2*(i-i_min)+SolnBlk_Fine.ICl+1][2*(j-j_min)+SolnBlk_Fine.JCl+1] = SolnBlk_Original.U[i][j];
	//= (SolnBlk_Original.Grid.Cell[i][j].A/area_total_fine)*SolnBlk_Original.U[i][j];
      }
    }

    // Prolong the east and west boundary states.
    for (int j = j_min-SolnBlk_Original.Nghost/2; j <= j_max+SolnBlk_Original.Nghost/2; j++) {
      SolnBlk_Fine.UoW[2*(j-j_min)+SolnBlk_Fine.JCl  ] = SolnBlk_Original.UoW[j];
      SolnBlk_Fine.UoW[2*(j-j_min)+SolnBlk_Fine.JCl+1] = SolnBlk_Original.UoW[j];
      SolnBlk_Fine.UoE[2*(j-j_min)+SolnBlk_Fine.JCl  ] = SolnBlk_Original.UoE[j];
      SolnBlk_Fine.UoE[2*(j-j_min)+SolnBlk_Fine.JCl+1] = SolnBlk_Original.UoE[j];
    }

    // Prolong the north and south boundary states.
    for (int i = i_min-SolnBlk_Original.Nghost/2; i <= i_max+SolnBlk_Original.Nghost/2; i++) {
      SolnBlk_Fine.UoS[2*(i-i_min)+SolnBlk_Fine.ICl  ] = SolnBlk_Original.UoS[i];
      SolnBlk_Fine.UoS[2*(i-i_min)+SolnBlk_Fine.ICl+1] = SolnBlk_Original.UoS[i];
      SolnBlk_Fine.UoN[2*(i-i_min)+SolnBlk_Fine.ICl  ] = SolnBlk_Original.UoN[i];
      SolnBlk_Fine.UoN[2*(i-i_min)+SolnBlk_Fine.ICl+1] = SolnBlk_Original.UoN[i];
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
int Restrict_Solution_Block(Electrostatic2D_Quad_Block &SolnBlk_Coarse,
			    Electrostatic2D_Quad_Block &SolnBlk_Original_SW,
			    Electrostatic2D_Quad_Block &SolnBlk_Original_SE,
			    Electrostatic2D_Quad_Block &SolnBlk_Original_NW,
			    Electrostatic2D_Quad_Block &SolnBlk_Original_NE) {

  int i_coarse, j_coarse, mesh_coarsening_permitted;

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

    // Copy the axisymmetric/planar flow indicator.
    SolnBlk_Coarse.Axisymmetric = SolnBlk_Original_SW.Axisymmetric;

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

      }
    }

    for (int j_fine = SolnBlk_Original_SW.JCl-SolnBlk_Original_SW.Nghost; 
	 j_fine <= SolnBlk_Original_SW.JCu+SolnBlk_Original_SW.Nghost; j_fine += 2) {
      j_coarse = (j_fine-SolnBlk_Original_SW.JCl)/2+SolnBlk_Coarse.JCl;
      SolnBlk_Coarse.UoW[j_coarse] = SolnBlk_Original_SW.UoW[j_fine];
      if (j_fine == SolnBlk_Original_SW.JCl-SolnBlk_Original_SW.Nghost) {
	SolnBlk_Coarse.UoW[j_coarse-1] = SolnBlk_Original_SW.UoW[j_fine];
      }
    }

    for (int i_fine = SolnBlk_Original_SW.ICl-SolnBlk_Original_SW.Nghost; 
	 i_fine <= SolnBlk_Original_SW.ICu+SolnBlk_Original_SW.Nghost; i_fine += 2) {
      i_coarse = (i_fine-SolnBlk_Original_SW.ICl)/2+SolnBlk_Coarse.ICl;
      SolnBlk_Coarse.UoS[i_coarse] = SolnBlk_Original_SW.UoS[i_fine];
      if (i_fine == SolnBlk_Original_SW.ICl-SolnBlk_Original_SW.Nghost) {
	SolnBlk_Coarse.UoS[i_coarse-1] = SolnBlk_Original_SW.UoS[i_fine];
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

      }
    }

    for (int j_fine = SolnBlk_Original_SE.JCl-SolnBlk_Original_SE.Nghost; 
	 j_fine <= SolnBlk_Original_SE.JCu+SolnBlk_Original_SE.Nghost; j_fine += 2) {
      j_coarse = (j_fine-SolnBlk_Original_SE.JCl)/2+SolnBlk_Coarse.JCl;
      SolnBlk_Coarse.UoE[j_coarse] = SolnBlk_Original_SE.UoE[j_fine];
      if (j_fine == SolnBlk_Original_SE.JCl-SolnBlk_Original_SE.Nghost) {
	SolnBlk_Coarse.UoE[j_coarse-1] = SolnBlk_Original_SE.UoE[j_fine];
      }
    }

    for (int i_fine = SolnBlk_Original_SE.ICl-SolnBlk_Original_SE.Nghost; 
	 i_fine <= SolnBlk_Original_SE.ICu+SolnBlk_Original_SE.Nghost; i_fine += 2) {
      i_coarse = (i_fine-SolnBlk_Original_SE.ICl)/2+(SolnBlk_Coarse.ICu-SolnBlk_Coarse.ICl+1)/2+SolnBlk_Coarse.ICl;
      SolnBlk_Coarse.UoS[i_coarse] = SolnBlk_Original_SE.UoS[i_fine];
      if (i_fine == SolnBlk_Original_SE.ICu+SolnBlk_Original_SE.Nghost) {
	SolnBlk_Coarse.UoS[i_coarse+1] = SolnBlk_Original_SE.UoS[i_fine];
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

      }
    }

    for (int j_fine = SolnBlk_Original_NW.JCl-SolnBlk_Original_NW.Nghost; 
	 j_fine <= SolnBlk_Original_NW.JCu+SolnBlk_Original_NW.Nghost; j_fine += 2) {
      j_coarse = (j_fine-SolnBlk_Original_NW.JCl)/2 + (SolnBlk_Coarse.JCu-SolnBlk_Coarse.JCl+1)/2+SolnBlk_Coarse.JCl;
      SolnBlk_Coarse.UoW[j_coarse] = SolnBlk_Original_NW.UoW[j_fine];
      if (j_fine == SolnBlk_Original_NW.JCu+SolnBlk_Original_NW.Nghost) {
	SolnBlk_Coarse.UoW[j_coarse+1] = SolnBlk_Original_NW.UoW[j_fine];
      }
    }

    for (int i_fine = SolnBlk_Original_NW.ICl-SolnBlk_Original_NW.Nghost; 
	 i_fine <= SolnBlk_Original_NW.ICu+SolnBlk_Original_NW.Nghost; i_fine += 2) {
      i_coarse = (i_fine-SolnBlk_Original_NW.ICl)/2+SolnBlk_Coarse.ICl;
      SolnBlk_Coarse.UoN[i_coarse] = SolnBlk_Original_NW.UoN[i_fine];
      if (i_fine == SolnBlk_Original_NW.ICl-SolnBlk_Original_NW.Nghost) {
	SolnBlk_Coarse.UoN[i_coarse-1] = SolnBlk_Original_NW.UoN[i_fine];
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

      }
    }

    for (int j_fine = SolnBlk_Original_NE.JCl-SolnBlk_Original_NE.Nghost; 
	 j_fine <= SolnBlk_Original_NE.JCu+SolnBlk_Original_NE.Nghost; j_fine += 2) {
      j_coarse = (j_fine-SolnBlk_Original_NE.JCl)/2 + (SolnBlk_Coarse.JCu-SolnBlk_Coarse.JCl+1)/2+SolnBlk_Coarse.JCl;
      SolnBlk_Coarse.UoE[j_coarse] = SolnBlk_Original_NE.UoE[j_fine];
      if (j_fine == SolnBlk_Original_NE.JCu+SolnBlk_Original_NE.Nghost) {
	SolnBlk_Coarse.UoE[j_coarse+1] = SolnBlk_Original_NE.UoE[j_fine];
      }
    }

    for (int i_fine = SolnBlk_Original_NE.ICl-SolnBlk_Original_NE.Nghost; 
	 i_fine <= SolnBlk_Original_NE.ICu+SolnBlk_Original_NE.Nghost; i_fine += 2) {
      i_coarse = (i_fine-SolnBlk_Original_NE.ICl)/2 + (SolnBlk_Coarse.ICu-SolnBlk_Coarse.ICl+1)/2+SolnBlk_Coarse.ICl;
      SolnBlk_Coarse.UoN[i_coarse] = SolnBlk_Original_NE.UoN[i_fine];
      if (i_fine == SolnBlk_Original_NE.ICu+SolnBlk_Original_NE.Nghost) {
	SolnBlk_Coarse.UoN[i_coarse+1] = SolnBlk_Original_NE.UoN[i_fine];
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
void ICs(Electrostatic2D_Quad_Block &SolnBlk,
         Electrostatic2D_Input_Parameters &IP,
         Electrostatic2DState *Uo) {

  Electrostatic2DState Ul, Ur;

  // Assign the initial data for the IVP of interest.
  switch(IP.i_ICs) {
  case IC_CONSTANT :
  case IC_UNIFORM :
    // Set the solution state to the initial state Uo[0].
    for (int j = SolnBlk.JCl-SolnBlk.Nghost; j <= SolnBlk.JCu+SolnBlk.Nghost; j++) {
      for (int i = SolnBlk.ICl-SolnBlk.Nghost; i <= SolnBlk.ICu+SolnBlk.Nghost; i++) {
	SolnBlk.U[i][j] = Uo[0];
      }
    }
    break;
  case IC_ELECTROSTATIC_CHANNEL :
    Ul = Electrostatic2DState(ZERO,ZERO,500.0);
    Ur = Electrostatic2DState(ZERO,ZERO, 60.0);
    for (int j = SolnBlk.JCl-SolnBlk.Nghost; j <= SolnBlk.JCu+SolnBlk.Nghost; j++) {
      for (int i = SolnBlk.ICl-SolnBlk.Nghost; i <= SolnBlk.ICu+SolnBlk.Nghost; i++) {
	if (SolnBlk.Grid.Cell[i][j].Xc.y >= HALF*IP.Box_Height) {
	  SolnBlk.U[i][j] = Ul;
	} else {
	  SolnBlk.U[i][j] = Ur;
	}
      }
    }
    break;
  case IC_DESOLVATION_CHAMBER :
    Ul = Electrostatic2DState(ZERO,ZERO,500.0);
    Ur = Electrostatic2DState(ZERO,ZERO, 60.0);
    for (int j = SolnBlk.JCl-SolnBlk.Nghost; j <= SolnBlk.JCu+SolnBlk.Nghost; j++) {
      for (int i = SolnBlk.ICl-SolnBlk.Nghost; i <= SolnBlk.ICu+SolnBlk.Nghost; i++) {
	if (SolnBlk.Grid.Cell[i][j].Xc.x <= 0.013) {
	  SolnBlk.U[i][j] = Ul;
	} else if (SolnBlk.Grid.Cell[i][j].Xc.x >= 0.014) {
	  SolnBlk.U[i][j] = Ur;
	} else {
	  SolnBlk.U[i][j] = Ur - (Ur-Ul)*(0.014-SolnBlk.Grid.Cell[i][j].Xc.x)/0.001;
	}
      }
    }
    break;
  default:
    // Set the solution state to the initial state Uo[0].
    for (int j = SolnBlk.JCl-SolnBlk.Nghost; j <= SolnBlk.JCu+SolnBlk.Nghost; j++) {
      for (int i = SolnBlk.ICl-SolnBlk.Nghost; i <= SolnBlk.ICu+SolnBlk.Nghost; i++) {
	SolnBlk.U[i][j] = Uo[0];
      }
    }
    break;
  };

  // Set default values for the boundary condition reference states.
  for (int j = SolnBlk.JCl-SolnBlk.Nghost; j <= SolnBlk.JCu+SolnBlk.Nghost; j++) {
    if (j >= SolnBlk.JCl && j <= SolnBlk.JCu) {
      SolnBlk.UoW[j] = SolnBlk.U[SolnBlk.ICl][j];
      SolnBlk.UoE[j] = SolnBlk.U[SolnBlk.ICu][j];
    } else if (j < SolnBlk.JCl) {
      SolnBlk.UoW[j] = SolnBlk.U[SolnBlk.ICl][SolnBlk.JCl];
      SolnBlk.UoE[j] = SolnBlk.U[SolnBlk.ICu][SolnBlk.JCl];
    } else {
      SolnBlk.UoW[j] = SolnBlk.U[SolnBlk.ICl][SolnBlk.JCu];
      SolnBlk.UoE[j] = SolnBlk.U[SolnBlk.ICu][SolnBlk.JCu];
    }
  }
  for (int i = SolnBlk.ICl-SolnBlk.Nghost; i <= SolnBlk.ICu+SolnBlk.Nghost; i++) {
    if (i >= SolnBlk.ICl && i <= SolnBlk.ICu) {
      SolnBlk.UoS[i] = SolnBlk.U[i][SolnBlk.JCl];
      SolnBlk.UoN[i] = SolnBlk.U[i][SolnBlk.JCu];
    } else if (i < SolnBlk.ICl) {
      SolnBlk.UoS[i] = SolnBlk.U[SolnBlk.ICl][SolnBlk.JCl];
      SolnBlk.UoN[i] = SolnBlk.U[SolnBlk.ICl][SolnBlk.JCu];
    } else {
      SolnBlk.UoS[i] = SolnBlk.U[SolnBlk.ICu][SolnBlk.JCl];
      SolnBlk.UoN[i] = SolnBlk.U[SolnBlk.ICu][SolnBlk.JCu];
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
void BCs(Electrostatic2D_Quad_Block &SolnBlk, Electrostatic2D_Input_Parameters &IP) {

  Vector2D dX;

  // WEST and EAST boundary conditions.
  for (int j = SolnBlk.JCl-SolnBlk.Nghost; j <= SolnBlk.JCu+SolnBlk.Nghost; j++) {

    // WEST boundary.
    if ((j >= SolnBlk.JCl && j <= SolnBlk.JCu) ||
	(j < SolnBlk.JCl && 
	 (SolnBlk.Grid.BCtypeS[SolnBlk.ICl-1] == BC_NONE ||
	  SolnBlk.Grid.BCtypeS[SolnBlk.ICl-1] == BC_PERIODIC ||
	  SolnBlk.Grid.BCtypeS[SolnBlk.ICl-1] == BC_CONSTANT_EXTRAPOLATION ||
	  SolnBlk.Grid.BCtypeS[SolnBlk.ICl-1] == BC_LINEAR_EXTRAPOLATION)) ||
	(j > SolnBlk.JCu && 
	 (SolnBlk.Grid.BCtypeN[SolnBlk.ICl-1] == BC_NONE ||
	  SolnBlk.Grid.BCtypeN[SolnBlk.ICl-1] == BC_PERIODIC ||
	  SolnBlk.Grid.BCtypeN[SolnBlk.ICl-1] == BC_CONSTANT_EXTRAPOLATION ||
	  SolnBlk.Grid.BCtypeN[SolnBlk.ICl-1] == BC_LINEAR_EXTRAPOLATION))) {

      switch(SolnBlk.Grid.BCtypeW[j]) {
      case BC_NONE :
	break;
      case BC_FIXED :
	SolnBlk.U[SolnBlk.ICl-1][j] = SolnBlk.UoW[j];
	SolnBlk.U[SolnBlk.ICl-2][j] = SolnBlk.UoW[j];
	break;
      case BC_CONSTANT_EXTRAPOLATION :
	SolnBlk.U[SolnBlk.ICl-1][j] = SolnBlk.U[SolnBlk.ICl][j];
	SolnBlk.U[SolnBlk.ICl-2][j] = SolnBlk.U[SolnBlk.ICl][j];
	break;
      case BC_LINEAR_EXTRAPOLATION :
	Linear_Reconstruction_LeastSquares_2(SolnBlk,SolnBlk.ICl,j,LIMITER_ONE);
	dX = SolnBlk.Grid.Cell[SolnBlk.ICl-1][j].Xc - SolnBlk.Grid.Cell[SolnBlk.ICl][j].Xc;
	SolnBlk.U[SolnBlk.ICl-1][j] = SolnBlk.U[SolnBlk.ICl][j] + 
	                              SolnBlk.dUdx[SolnBlk.ICl][j]*dX.x +
	                              SolnBlk.dUdy[SolnBlk.ICl][j]*dX.y;
	dX = SolnBlk.Grid.Cell[SolnBlk.ICl-2][j].Xc - SolnBlk.Grid.Cell[SolnBlk.ICl][j].Xc;
	SolnBlk.U[SolnBlk.ICl-2][j] = SolnBlk.U[SolnBlk.ICl][j] + 
	                              SolnBlk.dUdx[SolnBlk.ICl][j]*dX.x +
	                              SolnBlk.dUdy[SolnBlk.ICl][j]*dX.y;
	break;
      case BC_REFLECTION :
	SolnBlk.U[SolnBlk.ICl-1][j] = Reflect(SolnBlk.U[SolnBlk.ICl][j],
					      SolnBlk.Grid.nfaceW(SolnBlk.ICl,j));
	SolnBlk.U[SolnBlk.ICl-2][j] = Reflect(SolnBlk.U[SolnBlk.ICl+1][j],
					      SolnBlk.Grid.nfaceW(SolnBlk.ICl,j));
	break;
      case BC_WALL_VISCOUS_HEATFLUX :
	SolnBlk.U[SolnBlk.ICl-1][j] = Mirror(SolnBlk.U[SolnBlk.ICl][j],
					     SolnBlk.Grid.nfaceW(SolnBlk.ICl,j));
	SolnBlk.U[SolnBlk.ICl-2][j] = Mirror(SolnBlk.U[SolnBlk.ICl+1][j],
					     SolnBlk.Grid.nfaceW(SolnBlk.ICl,j));
	break;
      case BC_PERIODIC :
	SolnBlk.U[SolnBlk.ICl-1][j] = SolnBlk.U[SolnBlk.ICu-1][j];
	SolnBlk.U[SolnBlk.ICl-2][j] = SolnBlk.U[SolnBlk.ICu-2][j];
	break;
      default:
	SolnBlk.U[SolnBlk.ICl-1][j] = SolnBlk.U[SolnBlk.ICl][j];
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
	  SolnBlk.Grid.BCtypeS[SolnBlk.ICu+1] == BC_LINEAR_EXTRAPOLATION)) ||
	(j > SolnBlk.JCu && 
	 (SolnBlk.Grid.BCtypeN[SolnBlk.ICu+1] == BC_NONE ||
	  SolnBlk.Grid.BCtypeN[SolnBlk.ICu+1] == BC_PERIODIC ||
	  SolnBlk.Grid.BCtypeN[SolnBlk.ICu+1] == BC_CONSTANT_EXTRAPOLATION ||
	  SolnBlk.Grid.BCtypeN[SolnBlk.ICu+1] == BC_LINEAR_EXTRAPOLATION))) {

      switch(SolnBlk.Grid.BCtypeE[j]) {
      case BC_NONE :
	break;
      case BC_FIXED :
	SolnBlk.U[SolnBlk.ICu+1][j] = SolnBlk.UoE[j];
	SolnBlk.U[SolnBlk.ICu+2][j] = SolnBlk.UoE[j];
	break;
      case BC_CONSTANT_EXTRAPOLATION :
	SolnBlk.U[SolnBlk.ICu+1][j] = SolnBlk.U[SolnBlk.ICu][j];
	SolnBlk.U[SolnBlk.ICu+2][j] = SolnBlk.U[SolnBlk.ICu][j];
	break;
      case BC_LINEAR_EXTRAPOLATION :
	Linear_Reconstruction_LeastSquares_2(SolnBlk,SolnBlk.ICu,j,LIMITER_ONE);
	dX = SolnBlk.Grid.Cell[SolnBlk.ICu+1][j].Xc - SolnBlk.Grid.Cell[SolnBlk.ICu][j].Xc;
	SolnBlk.U[SolnBlk.ICu+1][j] = SolnBlk.U[SolnBlk.ICu][j] + 
                                      SolnBlk.dUdx[SolnBlk.ICu][j]*dX.x +
                                      SolnBlk.dUdy[SolnBlk.ICu][j]*dX.y;
	dX = SolnBlk.Grid.Cell[SolnBlk.ICu+2][j].Xc - SolnBlk.Grid.Cell[SolnBlk.ICu][j].Xc;
	SolnBlk.U[SolnBlk.ICu+2][j] = SolnBlk.U[SolnBlk.ICu][j] + 
                                      SolnBlk.dUdx[SolnBlk.ICu][j]*dX.x +
                                      SolnBlk.dUdy[SolnBlk.ICu][j]*dX.y;
	break;
      case BC_REFLECTION :
	SolnBlk.U[SolnBlk.ICu+1][j] = Reflect(SolnBlk.U[SolnBlk.ICu][j],
					      SolnBlk.Grid.nfaceE(SolnBlk.ICu,j));
	SolnBlk.U[SolnBlk.ICu+2][j] = Reflect(SolnBlk.U[SolnBlk.ICu-1][j],
					      SolnBlk.Grid.nfaceE(SolnBlk.ICu,j));
	break;
      case BC_WALL_VISCOUS_HEATFLUX :
	SolnBlk.U[SolnBlk.ICu+1][j] = Mirror(SolnBlk.U[SolnBlk.ICu][j],
					     SolnBlk.Grid.nfaceE(SolnBlk.ICu,j));
	SolnBlk.U[SolnBlk.ICu+2][j] = Mirror(SolnBlk.U[SolnBlk.ICu-1][j],
					     SolnBlk.Grid.nfaceE(SolnBlk.ICu,j));
	break;
      case BC_PERIODIC :
	SolnBlk.U[SolnBlk.ICu+1][j] = SolnBlk.U[SolnBlk.ICl+1][j];
	SolnBlk.U[SolnBlk.ICu+2][j] = SolnBlk.U[SolnBlk.ICl+2][j];
	break;
      default:
	SolnBlk.U[SolnBlk.ICu+1][j] = SolnBlk.U[SolnBlk.ICu][j];
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
      SolnBlk.U[i][SolnBlk.JCl-1] = SolnBlk.UoS[i];
      SolnBlk.U[i][SolnBlk.JCl-2] = SolnBlk.UoS[i];
      break;
    case BC_CONSTANT_EXTRAPOLATION :
      SolnBlk.U[i][SolnBlk.JCl-1] = SolnBlk.U[i][SolnBlk.JCl];
      SolnBlk.U[i][SolnBlk.JCl-2] = SolnBlk.U[i][SolnBlk.JCl];
      break;
    case BC_LINEAR_EXTRAPOLATION :
      Linear_Reconstruction_LeastSquares_2(SolnBlk,i,SolnBlk.JCl,LIMITER_ONE);
      dX = SolnBlk.Grid.Cell[i][SolnBlk.JCl-1].Xc - SolnBlk.Grid.Cell[i][SolnBlk.JCl].Xc;
      SolnBlk.U[i][SolnBlk.JCl-1] = SolnBlk.U[i][SolnBlk.JCl] + 
	                            SolnBlk.dUdx[i][SolnBlk.JCl]*dX.x +
	                            SolnBlk.dUdy[i][SolnBlk.JCl]*dX.y;
      dX = SolnBlk.Grid.Cell[i][SolnBlk.JCl-2].Xc - SolnBlk.Grid.Cell[i][SolnBlk.JCl].Xc;
      SolnBlk.U[i][SolnBlk.JCl-2] = SolnBlk.U[i][SolnBlk.JCl] + 
                	            SolnBlk.dUdx[i][SolnBlk.JCl]*dX.x +
	                            SolnBlk.dUdy[i][SolnBlk.JCl]*dX.y;
      break;
    case BC_REFLECTION :
      SolnBlk.U[i][SolnBlk.JCl-1] = Reflect(SolnBlk.U[i][SolnBlk.JCl],
					    SolnBlk.Grid.nfaceS(i,SolnBlk.JCl));
      SolnBlk.U[i][SolnBlk.JCl-2] = Reflect(SolnBlk.U[i][SolnBlk.JCl+1],
					    SolnBlk.Grid.nfaceS(i,SolnBlk.JCl));
      break;
    case BC_WALL_VISCOUS_HEATFLUX :
      SolnBlk.U[i][SolnBlk.JCl-1] = Mirror(SolnBlk.U[i][SolnBlk.JCl],
					   SolnBlk.Grid.nfaceS(i,SolnBlk.JCl));
      SolnBlk.U[i][SolnBlk.JCl-2] = Mirror(SolnBlk.U[i][SolnBlk.JCl+1],
					   SolnBlk.Grid.nfaceS(i,SolnBlk.JCl));
      break;
    case BC_PERIODIC :
      SolnBlk.U[i][SolnBlk.JCl-1] = SolnBlk.U[i][SolnBlk.JCu-1];
      SolnBlk.U[i][SolnBlk.JCl-2] = SolnBlk.U[i][SolnBlk.JCu-2];
      break;
    default:
      SolnBlk.U[i][SolnBlk.JCl-1] = SolnBlk.U[i][SolnBlk.JCl];
      SolnBlk.U[i][SolnBlk.JCl-2] = SolnBlk.U[i][SolnBlk.JCl];
      break;
    };
    
    // NORTH boundary.
    switch(SolnBlk.Grid.BCtypeN[i]) {
    case BC_NONE :
      break;
    case BC_FIXED :
      SolnBlk.U[i][SolnBlk.JCu+1] = SolnBlk.UoN[i];
      SolnBlk.U[i][SolnBlk.JCu+2] = SolnBlk.UoN[i];
      break;
    case BC_CONSTANT_EXTRAPOLATION :
      SolnBlk.U[i][SolnBlk.JCu+1] = SolnBlk.U[i][SolnBlk.JCu];
      SolnBlk.U[i][SolnBlk.JCu+2] = SolnBlk.U[i][SolnBlk.JCu];
      break;
    case BC_LINEAR_EXTRAPOLATION :
      Linear_Reconstruction_LeastSquares_2(SolnBlk,i,SolnBlk.JCu,LIMITER_BARTH_JESPERSEN);
      dX = SolnBlk.Grid.Cell[i][SolnBlk.JCu+1].Xc - SolnBlk.Grid.Cell[i][SolnBlk.JCu].Xc;
      SolnBlk.U[i][SolnBlk.JCu+1] = SolnBlk.U[i][SolnBlk.JCu] + 
	                            SolnBlk.dUdx[i][SolnBlk.JCu]*dX.x +
	                            SolnBlk.dUdy[i][SolnBlk.JCu]*dX.y;
      dX = SolnBlk.Grid.Cell[i][SolnBlk.JCu+2].Xc - SolnBlk.Grid.Cell[i][SolnBlk.JCu].Xc;
      SolnBlk.U[i][SolnBlk.JCu+2] = SolnBlk.U[i][SolnBlk.JCu] + 
	                            SolnBlk.dUdx[i][SolnBlk.JCu]*dX.x +
	                            SolnBlk.dUdy[i][SolnBlk.JCu]*dX.y;
      break;
    case BC_REFLECTION :
      SolnBlk.U[i][SolnBlk.JCu+1] = Reflect(SolnBlk.U[i][SolnBlk.JCu],
					    SolnBlk.Grid.nfaceN(i,SolnBlk.JCu));
      SolnBlk.U[i][SolnBlk.JCu+2] = Reflect(SolnBlk.U[i][SolnBlk.JCu-1],
					    SolnBlk.Grid.nfaceN(i,SolnBlk.JCu));
      break;
    case BC_WALL_VISCOUS_HEATFLUX :
      SolnBlk.U[i][SolnBlk.JCu+1] = Mirror(SolnBlk.U[i][SolnBlk.JCu],
					   SolnBlk.Grid.nfaceN(i,SolnBlk.JCu));
      SolnBlk.U[i][SolnBlk.JCu+2] = Mirror(SolnBlk.U[i][SolnBlk.JCu-1],
					   SolnBlk.Grid.nfaceN(i,SolnBlk.JCu));
      break;
    case BC_PERIODIC :
      SolnBlk.U[i][SolnBlk.JCu+1] = SolnBlk.U[i][SolnBlk.JCl+1];
      SolnBlk.U[i][SolnBlk.JCu+2] = SolnBlk.U[i][SolnBlk.JCl+2];
      break;
    default:
      SolnBlk.U[i][SolnBlk.JCu+1] = SolnBlk.U[i][SolnBlk.JCu];
      SolnBlk.U[i][SolnBlk.JCu+2] = SolnBlk.U[i][SolnBlk.JCu];
      break;
    };

  }

   // BC fix for corner points with burning surfaces on either side.
   if ((SolnBlk.Grid.BCtypeW[SolnBlk.JCl] == BC_REFLECTION ||
        SolnBlk.Grid.BCtypeW[SolnBlk.JCl] == BC_WALL_VISCOUS_HEATFLUX) &&
       (SolnBlk.Grid.BCtypeS[SolnBlk.ICl] == BC_REFLECTION ||
        SolnBlk.Grid.BCtypeS[SolnBlk.ICl] == BC_WALL_VISCOUS_HEATFLUX)) {
     SolnBlk.U[SolnBlk.ICl-2][SolnBlk.JCl-2] = HALF*(SolnBlk.U[SolnBlk.ICl-2][SolnBlk.JCl  ]+
						     SolnBlk.U[SolnBlk.ICl  ][SolnBlk.JCl-2]);
     SolnBlk.U[SolnBlk.ICl-1][SolnBlk.JCl-2] = HALF*(SolnBlk.U[SolnBlk.ICl-1][SolnBlk.JCl  ]+
						     SolnBlk.U[SolnBlk.ICl  ][SolnBlk.JCl-2]);
     SolnBlk.U[SolnBlk.ICl-2][SolnBlk.JCl-1] = HALF*(SolnBlk.U[SolnBlk.ICl-2][SolnBlk.JCl  ]+
						     SolnBlk.U[SolnBlk.ICl  ][SolnBlk.JCl-1]);
     SolnBlk.U[SolnBlk.ICl-1][SolnBlk.JCl-1] = HALF*(SolnBlk.U[SolnBlk.ICl-1][SolnBlk.JCl  ]+
						     SolnBlk.U[SolnBlk.ICl  ][SolnBlk.JCl-1]);
   }
  
   if ((SolnBlk.Grid.BCtypeW[SolnBlk.JCu] == BC_REFLECTION ||
        SolnBlk.Grid.BCtypeW[SolnBlk.JCu] == BC_WALL_VISCOUS_HEATFLUX) &&
       (SolnBlk.Grid.BCtypeN[SolnBlk.ICl] == BC_REFLECTION ||
        SolnBlk.Grid.BCtypeN[SolnBlk.ICl] == BC_WALL_VISCOUS_HEATFLUX)) {
     SolnBlk.U[SolnBlk.ICl-1][SolnBlk.JCu+1] = HALF*(SolnBlk.U[SolnBlk.ICl-1][SolnBlk.JCu]+
						     SolnBlk.U[SolnBlk.ICl][SolnBlk.JCu+1]);
     SolnBlk.U[SolnBlk.ICl-2][SolnBlk.JCu+1] = HALF*(SolnBlk.U[SolnBlk.ICl-2][SolnBlk.JCu]+
						     SolnBlk.U[SolnBlk.ICl][SolnBlk.JCu+1]);
     SolnBlk.U[SolnBlk.ICl-1][SolnBlk.JCu+2] = HALF*(SolnBlk.U[SolnBlk.ICl-1][SolnBlk.JCu]+
						     SolnBlk.U[SolnBlk.ICl][SolnBlk.JCu+2]);
     SolnBlk.U[SolnBlk.ICl-2][SolnBlk.JCu+2] = HALF*(SolnBlk.U[SolnBlk.ICl-2][SolnBlk.JCu]+
						     SolnBlk.U[SolnBlk.ICl][SolnBlk.JCu+2]);
   }
  
   if ((SolnBlk.Grid.BCtypeE[SolnBlk.JCl] == BC_REFLECTION ||
        SolnBlk.Grid.BCtypeE[SolnBlk.JCl] == BC_WALL_VISCOUS_HEATFLUX) &&
       (SolnBlk.Grid.BCtypeS[SolnBlk.ICu] == BC_REFLECTION ||
        SolnBlk.Grid.BCtypeS[SolnBlk.ICu] == BC_WALL_VISCOUS_HEATFLUX)) {
     SolnBlk.U[SolnBlk.ICu+1][SolnBlk.JCl-1] = HALF*(SolnBlk.U[SolnBlk.ICu+1][SolnBlk.JCl]+
						     SolnBlk.U[SolnBlk.ICu][SolnBlk.JCl-1]);
     SolnBlk.U[SolnBlk.ICu+2][SolnBlk.JCl-1] = HALF*(SolnBlk.U[SolnBlk.ICu+2][SolnBlk.JCl  ] +
						     SolnBlk.U[SolnBlk.ICu  ][SolnBlk.JCl-1]);
     SolnBlk.U[SolnBlk.ICu+1][SolnBlk.JCl-2] = HALF*(SolnBlk.U[SolnBlk.ICu+1][SolnBlk.JCl  ] +
						     SolnBlk.U[SolnBlk.ICu  ][SolnBlk.JCl-2]);
     SolnBlk.U[SolnBlk.ICu+2][SolnBlk.JCl-2] = HALF*(SolnBlk.U[SolnBlk.ICu+2][SolnBlk.JCl  ] +
						     SolnBlk.U[SolnBlk.ICu  ][SolnBlk.JCl-2]);
   }
  
   if ((SolnBlk.Grid.BCtypeE[SolnBlk.JCu] == BC_REFLECTION ||
        SolnBlk.Grid.BCtypeE[SolnBlk.JCu] == BC_WALL_VISCOUS_HEATFLUX) &&
       (SolnBlk.Grid.BCtypeN[SolnBlk.ICu] == BC_REFLECTION ||
        SolnBlk.Grid.BCtypeN[SolnBlk.ICu] == BC_WALL_VISCOUS_HEATFLUX)) {
     SolnBlk.U[SolnBlk.ICu+1][SolnBlk.JCu+1] = HALF*(SolnBlk.U[SolnBlk.ICu+1][SolnBlk.JCu  ]+
						     SolnBlk.U[SolnBlk.ICu  ][SolnBlk.JCu+1]);
     SolnBlk.U[SolnBlk.ICu+2][SolnBlk.JCu+1] = HALF*(SolnBlk.U[SolnBlk.ICu+2][SolnBlk.JCu  ] +
						     SolnBlk.U[SolnBlk.ICu  ][SolnBlk.JCu+1]);
     SolnBlk.U[SolnBlk.ICu+1][SolnBlk.JCu+2] = HALF*(SolnBlk.U[SolnBlk.ICu+1][SolnBlk.JCu  ] +
						     SolnBlk.U[SolnBlk.ICu  ][SolnBlk.JCu+2]);
     SolnBlk.U[SolnBlk.ICu+2][SolnBlk.JCu+2] = HALF*(SolnBlk.U[SolnBlk.ICu+2][SolnBlk.JCu  ] +
						     SolnBlk.U[SolnBlk.ICu  ][SolnBlk.JCu+2]);
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
double CFL(Electrostatic2D_Quad_Block &SolnBlk,
	   Electrostatic2D_Input_Parameters &IP) {

  double dtMin = MILLION, d_i, d_j;

  for (int j = SolnBlk.JCl-SolnBlk.Nghost; j <= SolnBlk.JCu+SolnBlk.Nghost; j++) {
    for (int i = SolnBlk.ICl-SolnBlk.Nghost; i <= SolnBlk.ICu+SolnBlk.Nghost; i++) {
      if (i < SolnBlk.ICl || i > SolnBlk.ICu || j < SolnBlk.JCl || j > SolnBlk.JCu) {
	SolnBlk.dt[i][j] = ZERO;
      } else {

	d_i = TWO*SolnBlk.Grid.Cell[i][j].A/(SolnBlk.Grid.lfaceE(i,j)+SolnBlk.Grid.lfaceW(i,j));
	d_j = TWO*SolnBlk.Grid.Cell[i][j].A/(SolnBlk.Grid.lfaceN(i,j)+SolnBlk.Grid.lfaceS(i,j));

	SolnBlk.dt[i][j] = HALF*min(d_i*d_i,d_j*d_j);
	dtMin = min(dtMin,SolnBlk.dt[i][j]);

      }
    }
  }

  for (int j = SolnBlk.JCl-SolnBlk.Nghost; j <= SolnBlk.JCu+SolnBlk.Nghost; j++)
    for (int i = SolnBlk.ICl-SolnBlk.Nghost; i <= SolnBlk.ICu+SolnBlk.Nghost; i++)
      if (i < SolnBlk.ICl || i > SolnBlk.ICu || j < SolnBlk.JCl || j > SolnBlk.JCu)
	SolnBlk.dt[i][j] = dtMin;

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
void Set_Global_TimeStep(Electrostatic2D_Quad_Block &SolnBlk,
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
double L1_Norm_Residual(Electrostatic2D_Quad_Block &SolnBlk) {

  double l1_norm = ZERO;

  for (int j = SolnBlk.JCl; j <= SolnBlk.JCu; j++) {
    for (int i = SolnBlk.ICl; i <= SolnBlk.ICu; i++) {
      l1_norm += fabs(SolnBlk.dUdt[i][j][0].V);
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
double L2_Norm_Residual(Electrostatic2D_Quad_Block &SolnBlk) {

  double l2_norm = ZERO;

  for (int j = SolnBlk.JCl; j <= SolnBlk.JCu; j++) {
    for (int i = SolnBlk.ICl; i <= SolnBlk.ICu; i++) {
      l2_norm += sqr(SolnBlk.dUdt[i][j][0].V);
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
double Max_Norm_Residual(Electrostatic2D_Quad_Block &SolnBlk) {

  double max_norm = ZERO;

  for (int j = SolnBlk.JCl; j <= SolnBlk.JCu; j++) {
    for (int i = SolnBlk.ICl; i <= SolnBlk.ICu; i++) {
      max_norm = max(max_norm,fabs(SolnBlk.dUdt[i][j][0].V));
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
void Linear_Reconstruction_GreenGauss(Electrostatic2D_Quad_Block &SolnBlk,
				      const int i,
                                      const int j,
                                      const int Limiter) {

  int error_flag;
  int n_pts, i_index[8], j_index[8];
  double u0Min, u0Max, uQuad[4];
  double DxDx_ave, DxDy_ave, DyDy_ave;
  double l_north, l_south, l_east, l_west;
  Vector2D n_north, n_south, n_east, n_west, dX;
  Electrostatic2DState U_nw, U_ne, U_sw, U_se, U_face,
                       DU, DUDx_ave, DUDy_ave;

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
	       SolnBlk.Grid.BCtypeW[j] == BC_LINEAR_EXTRAPOLATION) {
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
	       SolnBlk.Grid.BCtypeE[j] == BC_LINEAR_EXTRAPOLATION) {
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
	       SolnBlk.Grid.BCtypeS[i] == BC_LINEAR_EXTRAPOLATION) {
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
	       SolnBlk.Grid.BCtypeN[i] == BC_LINEAR_EXTRAPOLATION) {
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
      error_flag = Green_Gauss_Integration(SolnBlk.Grid.nodeNW(i,j).X,SolnBlk.UnNW(i,j),
					   SolnBlk.Grid.nodeSW(i,j).X,SolnBlk.UnSW(i,j),
					   SolnBlk.Grid.nodeSE(i,j).X,SolnBlk.UnSE(i,j),
					   SolnBlk.Grid.nodeNE(i,j).X,SolnBlk.UnNE(i,j),
					   SolnBlk.dUdx[i][j],SolnBlk.dUdy[i][j]);
      //if (error_flag) return error_flag;
//       W_nw = SolnBlk.UnNW(i,j);
//       W_ne = SolnBlk.UnNE(i,j);
//       W_sw = SolnBlk.UnSW(i,j);
//       W_se = SolnBlk.UnSE(i,j);

//       l_north = SolnBlk.Grid.lfaceN(i,j);
//       l_south = SolnBlk.Grid.lfaceS(i,j);
//       l_east = SolnBlk.Grid.lfaceE(i,j);
//       l_west = SolnBlk.Grid.lfaceW(i,j);

//       n_north = SolnBlk.Grid.nfaceN(i,j);
//       n_south = SolnBlk.Grid.nfaceS(i,j);
//       n_east = SolnBlk.Grid.nfaceE(i,j);
//       n_west = SolnBlk.Grid.nfaceW(i,j);

//       W_face = HALF*(W_nw+W_ne)*l_north; 
//       SolnBlk.dUdx[i][j] = W_face*n_north.x;
//       SolnBlk.dUdy[i][j] = W_face*n_north.y;

//       W_face = HALF*(W_sw+W_se)*l_south; 
//       SolnBlk.dUdx[i][j] += W_face*n_south.x;
//       SolnBlk.dUdy[i][j] += W_face*n_south.y;

//       W_face = HALF*(W_ne+W_se)*l_east; 
//       SolnBlk.dUdx[i][j] += W_face*n_east.x;
//       SolnBlk.dUdy[i][j] += W_face*n_east.y;

//       W_face = HALF*(W_nw+W_sw)*l_west; 
//       SolnBlk.dUdx[i][j] += W_face*n_west.x;
//       SolnBlk.dUdy[i][j] += W_face*n_west.y;

//       SolnBlk.dUdx[i][j] = SolnBlk.dUdx[i][j]/SolnBlk.Grid.Cell[i][j].A;
//       SolnBlk.dUdy[i][j] = SolnBlk.dUdy[i][j]/SolnBlk.Grid.Cell[i][j].A;

    } else {

      // If <8 neighbours are used, apply least-squares reconstruction
      DUDx_ave = ELECTROSTATIC2D_VACUUM;
      DUDy_ave = ELECTROSTATIC2D_VACUUM;
      DxDx_ave = ZERO;
      DxDy_ave = ZERO;
      DyDy_ave = ZERO;

      for (int n2 = 0; n2 < n_pts; n2++) {
	dX = SolnBlk.Grid.Cell[i_index[n2]][j_index[n2]].Xc - SolnBlk.Grid.Cell[i][j].Xc;
	DU = SolnBlk.U[i_index[n2]][j_index[n2]] - SolnBlk.U[i][j];
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
      SolnBlk.dUdx[i][j] = (DUDx_ave*DyDy_ave-DUDy_ave*DxDy_ave)/
	                   (DxDx_ave*DyDy_ave-DxDy_ave*DxDy_ave);
      SolnBlk.dUdy[i][j] = (DUDy_ave*DxDx_ave-DUDx_ave*DxDy_ave)/
	                   (DxDx_ave*DyDy_ave-DxDy_ave*DxDy_ave);
    }

  } else {
    SolnBlk.dUdx[i][j] = ELECTROSTATIC2D_VACUUM;
    SolnBlk.dUdy[i][j] = ELECTROSTATIC2D_VACUUM;

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
void Linear_Reconstruction_GreenGauss(Electrostatic2D_Quad_Block &SolnBlk,
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
void Linear_Reconstruction_LeastSquares(Electrostatic2D_Quad_Block &SolnBlk,
				        const int i,
                                        const int j,
                                        const int Limiter) {

  int n_pts, i_index[8], j_index[8], motion_flag = OFF, n_pts_temp;
  double u0Min, u0Max, uQuad[4];
  double DxDx_ave, DxDy_ave, DyDy_ave;
  Vector2D dX;
  Electrostatic2DState DU, DUDx_ave, DUDy_ave;

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
	       SolnBlk.Grid.BCtypeW[j] == BC_LINEAR_EXTRAPOLATION) {
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
	       SolnBlk.Grid.BCtypeE[j] == BC_LINEAR_EXTRAPOLATION) {
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
	       SolnBlk.Grid.BCtypeS[i] == BC_LINEAR_EXTRAPOLATION) {
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
	       SolnBlk.Grid.BCtypeN[i] == BC_LINEAR_EXTRAPOLATION) {
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
    DUDx_ave = ELECTROSTATIC2D_VACUUM;
    DUDy_ave = ELECTROSTATIC2D_VACUUM;
    DxDx_ave = ZERO;
    DxDy_ave = ZERO;
    DyDy_ave = ZERO;

    for (int n2 = 0; n2 < n_pts; n2++) {
      dX = SolnBlk.Grid.Cell[i_index[n2]][j_index[n2]].Xc - SolnBlk.Grid.Cell[i][j].Xc;
      DU = SolnBlk.U[i_index[n2]][j_index[n2]] - SolnBlk.U[i][j];
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
    SolnBlk.dUdx[i][j] = (DUDx_ave*DyDy_ave-DUDy_ave*DxDy_ave)/
                         (DxDx_ave*DyDy_ave-DxDy_ave*DxDy_ave);
    SolnBlk.dUdy[i][j] = (DUDy_ave*DxDx_ave-DUDx_ave*DxDy_ave)/
                         (DxDx_ave*DyDy_ave-DxDy_ave*DxDy_ave);

  } else {
    SolnBlk.dUdx[i][j] = ELECTROSTATIC2D_VACUUM;
    SolnBlk.dUdy[i][j] = ELECTROSTATIC2D_VACUUM;

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
void Linear_Reconstruction_LeastSquares_2(Electrostatic2D_Quad_Block &SolnBlk,
				          const int i,
                                          const int j,
                                          const int Limiter) {

  int n_pts, i_index[8], j_index[8];
  double u0Min, u0Max, uQuad[4];
  double DxDx_ave, DxDy_ave, DyDy_ave;
  Vector2D dX;
  Electrostatic2DState DU, DUDx_ave, DUDy_ave;

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
	       SolnBlk.Grid.BCtypeW[j] == BC_LINEAR_EXTRAPOLATION) {
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
	       SolnBlk.Grid.BCtypeE[j] == BC_LINEAR_EXTRAPOLATION) {
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
	       SolnBlk.Grid.BCtypeS[i] == BC_LINEAR_EXTRAPOLATION) {
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
	       SolnBlk.Grid.BCtypeN[i] == BC_LINEAR_EXTRAPOLATION) {
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
    DUDx_ave = ELECTROSTATIC2D_VACUUM;
    DUDy_ave = ELECTROSTATIC2D_VACUUM;
    DxDx_ave = ZERO;
    DxDy_ave = ZERO;
    DyDy_ave = ZERO;

    for (int n2 = 0; n2 < n_pts; n2++) {
      dX = SolnBlk.Grid.Cell[i_index[n2]][j_index[n2]].Xc - SolnBlk.Grid.Cell[i][j].Xc;
      DU = SolnBlk.U[i_index[n2]][j_index[n2]] - SolnBlk.U[i][j];
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
    SolnBlk.dUdx[i][j] = (DUDx_ave*DyDy_ave-DUDy_ave*DxDy_ave)/
                         (DxDx_ave*DyDy_ave-DxDy_ave*DxDy_ave);
    SolnBlk.dUdy[i][j] = (DUDy_ave*DxDx_ave-DUDx_ave*DxDy_ave)/
	                 (DxDx_ave*DyDy_ave-DxDy_ave*DxDy_ave);

  } else {
    SolnBlk.dUdx[i][j] = ELECTROSTATIC2D_VACUUM;
    SolnBlk.dUdy[i][j] = ELECTROSTATIC2D_VACUUM;

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
void Linear_Reconstruction_LeastSquares(Electrostatic2D_Quad_Block &SolnBlk,
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
void Residual_Smoothing(Electrostatic2D_Quad_Block &SolnBlk,
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
				   Electrostatic2D_Input_Parameters &IP,
                                   int &number_refinement_criteria,
                                   Electrostatic2D_Quad_Block &SolnBlk) {

  double grad_V_x, grad_V_y, grad_V_abs, grad_V_criteria, grad_V_criteria_max,
         div_E, div_E_criteria, div_E_criteria_max,
         curl_E_z, curl_E_abs, curl_E_criteria, curl_E_criteria_max;
  int refinement_criteria_number;

  // Set the number of refinement criteria to be used:
  // (1) Refinement criteria based on the gradient of the potential field;
  // (2) Refinement criteria based on the divergence of the electric field vector;
  // (3) Refinement criteria based on the curl of the electric field vector;
  number_refinement_criteria = IP.Number_of_Refinement_Criteria;

  // Initialize the refinement criteria for the block.
  grad_V_criteria_max = ZERO;
  div_E_criteria_max = ZERO;
  curl_E_criteria_max = ZERO;

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
	// potential field.
	if (IP.Refinement_Criteria_Gradient_Potential_Field) {
	  grad_V_x = SolnBlk.dUdx[i][j].V;
	  grad_V_y = SolnBlk.dUdy[i][j].V;
	  grad_V_abs = sqrt(sqr(grad_V_x) + sqr(grad_V_y));
	  grad_V_criteria = sqrt(SolnBlk.Grid.Cell[i][j].A)*grad_V_abs/SolnBlk.U[i][j].V;
	} else {
	  grad_V_criteria = ONE;
	}
	grad_V_criteria_max = max(grad_V_criteria_max,grad_V_criteria);

	// Evaluate refinement criteria #2 based on the divergence of the
	// electric field.
	if (IP.Refinement_Criteria_Divergence_Electric_Field) {
	  div_E = SolnBlk.dUdx[i][j].E.x + SolnBlk.dUdy[i][j].E.y;
	  div_E_criteria = sqrt(SolnBlk.Grid.Cell[i][j].A)*fabs(div_E)/SolnBlk.U[i][j].E.abs();
	} else {
	  div_E_criteria = ONE;
	}
	div_E_criteria_max = max(div_E_criteria_max,div_E_criteria);

	// Evaluate refinement criteria #3 based on the curl of the
	// electric field.
	if (IP.Refinement_Criteria_Curl_Electric_Field) {
	  curl_E_z = SolnBlk.dUdx[i][j].E.y - SolnBlk.dUdy[i][j].E.x; 
	  curl_E_abs = sqrt(sqr(curl_E_z)); 
	  curl_E_criteria = sqrt(SolnBlk.Grid.Cell[i][j].A)*curl_E_abs/SolnBlk.U[i][j].E.abs();
	} else {
	  curl_E_criteria = ONE;
	}
	curl_E_criteria_max = max(curl_E_criteria_max,curl_E_criteria);

      }

    }
  }

  // Return the refinement criteria.
  refinement_criteria_number = 0;
  if (IP.Refinement_Criteria_Gradient_Potential_Field) {
    refinement_criteria[refinement_criteria_number] = grad_V_criteria_max;
    refinement_criteria_number++;
  }
  if (IP.Refinement_Criteria_Divergence_Electric_Field) {
    refinement_criteria[refinement_criteria_number] = div_E_criteria_max;
    refinement_criteria_number++;
  }
  if (IP.Refinement_Criteria_Curl_Electric_Field) {
    refinement_criteria[refinement_criteria_number] = curl_E_criteria_max;
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
void Fix_Refined_Block_Boundaries(Electrostatic2D_Quad_Block &SolnBlk,
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
void Unfix_Refined_Block_Boundaries(Electrostatic2D_Quad_Block &SolnBlk) {

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
void Apply_Boundary_Flux_Corrections(Electrostatic2D_Quad_Block &SolnBlk,
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
void Apply_Boundary_Flux_Corrections_Multistage_Explicit(Electrostatic2D_Quad_Block &SolnBlk,
                                                         const int i_stage,
                                                         Electrostatic2D_Input_Parameters &IP,
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
    omega = MultiStage_Optimally_Smoothing(i_stage,IP.N_Stage,LIMITER_ONE);
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
int dUdt_Residual_Evaluation(Electrostatic2D_Quad_Block &SolnBlk,
			     Electrostatic2D_Input_Parameters &IP) {

  Electrostatic2DState Ul, Ur, Uu, Ud, dUdx, dUdy, Flux;
  Vector2D Xl, Xr, Xu, Xd;
  int diamond_path;

  // Perform the linear reconstruction within each cell of the
  // computational grid for this stage if required.
  if (IP.i_Flux_Reconstruction == VISCOUS_RECONSTRUCTION_ARITHMETIC_AVERAGE || IP.Axisymmetric) {
    switch(IP.i_Reconstruction) {
    case RECONSTRUCTION_GREEN_GAUSS :
      Linear_Reconstruction_GreenGauss(SolnBlk,LIMITER_ONE);
      break;
    case RECONSTRUCTION_LINEAR_LEAST_SQUARES :
      Linear_Reconstruction_LeastSquares(SolnBlk,LIMITER_ONE);
      break;
    };
  }

  // Evaluate the time rate of change of the solution (i.e., the
  // solution residuals) using a second-order limited upwind scheme
  // with a variety of flux functions.

  // Add i-direction (zeta-direction) fluxes.
  for (int j = SolnBlk.JCl-1; j <= SolnBlk.JCu+1; j++) {

    SolnBlk.dUdt[SolnBlk.ICl-1][j][0] = ELECTROSTATIC2D_VACUUM;

    for (int i = SolnBlk.ICl-1; i <= SolnBlk.ICu; i++) {

      SolnBlk.dUdt[i+1][j][0] = ELECTROSTATIC2D_VACUUM;

      if (j >= SolnBlk.JCl && j <= SolnBlk.JCu) {

	// Evaluate the cell interface i-direction flux.

	// Select the appropriate nodes and states.
	Ul = SolnBlk.U[i  ][j]; Xl = SolnBlk.Grid.Cell[i  ][j].Xc;
	Ur = SolnBlk.U[i+1][j]; Xr = SolnBlk.Grid.Cell[i+1][j].Xc;

	switch(IP.i_Flux_Reconstruction) {
	case VISCOUS_RECONSTRUCTION_ARITHMETIC_AVERAGE :
	  Flux = FluxArithmetic_n(SolnBlk.Grid.xfaceE(i,j),
				  Xl,Ul,SolnBlk.dUdx[i  ][j],SolnBlk.dUdy[i  ][j],
				  Xr,Ur,SolnBlk.dUdx[i+1][j],SolnBlk.dUdy[i+1][j],
				  SolnBlk.Grid.nfaceE(i,j));
	  break;
	case VISCOUS_RECONSTRUCTION_DIAMOND_PATH :
	  Xu = SolnBlk.Grid.Node[i+1][j+1].X;
	  Xd = SolnBlk.Grid.Node[i+1][j  ].X;
	  if (i == SolnBlk.ICl-1 && SolnBlk.Grid.BCtypeW[j] == BC_FIXED) {
	    // WEST face of cell (i+1,j) is a FIXED boundary.
	    diamond_path = DIAMONDPATH_RIGHT_TRIANGLE_HEATFLUX;
	    Uu = SolnBlk.UoW[j]; Ud = Uu;
	  } else if (i == SolnBlk.ICu && SolnBlk.Grid.BCtypeE[j] == BC_FIXED) {
	    // EAST face of cell (i,j) is a FIXED boundary.
	    diamond_path = DIAMONDPATH_LEFT_TRIANGLE_HEATFLUX;
	    Uu = SolnBlk.UoE[j]; Ud = Uu;
	  } else {
	    // EAST face is either a normal cell or possibly a non-fixed
	    // boundary condition.
	    diamond_path = DIAMONDPATH_QUADRILATERAL_RECONSTRUCTION;
	    Uu = SolnBlk.UnNE(i,j); Ud = SolnBlk.UnSE(i,j);
	  }
	  Flux = FluxDiamondPath_n(SolnBlk.Grid.xfaceE(i,j),
				   Xl,Ul,Xu,Uu,Xr,Ur,Xd,Ud,
				   SolnBlk.Grid.nfaceE(i,j),
				   diamond_path);
	  break;
	};

	// Evaluate cell-averaged solution changes.
	SolnBlk.dUdt[i  ][j][0] += Flux*SolnBlk.Grid.lfaceE(i,j)/SolnBlk.Grid.Cell[i][j].A;
	SolnBlk.dUdt[i+1][j][0] -= Flux*SolnBlk.Grid.lfaceW(i+1,j)/SolnBlk.Grid.Cell[i+1][j].A;

	// Include axisymmetric source terms if required.
	if (SolnBlk.Axisymmetric) {
	  SolnBlk.dUdt[i][j][0] -= SolnBlk.U[i][j].S(SolnBlk.Grid.Cell[i][j].Xc,
						     SolnBlk.dUdy[i][j]);
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
      SolnBlk.dUdt[SolnBlk.ICl-1][j][0] = ELECTROSTATIC2D_VACUUM;
      SolnBlk.dUdt[SolnBlk.ICu+1][j][0] = ELECTROSTATIC2D_VACUUM;
    }

  }

  // Add j-direction (eta-direction) fluxes.
  for (int i = SolnBlk.ICl; i <= SolnBlk.ICu; i++) {
    for (int j = SolnBlk.JCl-1; j <= SolnBlk.JCu; j++) {

      // Evaluate the cell interface j-direction flux.

      // Select the appropriate nodes and states.
      Ul = SolnBlk.U[i][j  ]; Xl = SolnBlk.Grid.Cell[i][j  ].Xc;
      Ur = SolnBlk.U[i][j+1]; Xr = SolnBlk.Grid.Cell[i][j+1].Xc;

      switch(IP.i_Flux_Reconstruction) {
      case VISCOUS_RECONSTRUCTION_ARITHMETIC_AVERAGE :
	Flux = FluxArithmetic_n(SolnBlk.Grid.xfaceN(i,j),
				Xl,Ul,SolnBlk.dUdx[i][j  ],SolnBlk.dUdy[i][j  ],
				Xr,Ur,SolnBlk.dUdx[i][j+1],SolnBlk.dUdy[i][j+1],
				SolnBlk.Grid.nfaceN(i,j));
	break;
      case VISCOUS_RECONSTRUCTION_DIAMOND_PATH :
	Xu = SolnBlk.Grid.Node[i  ][j+1].X;
	Xd = SolnBlk.Grid.Node[i+1][j+1].X;
	if (j == SolnBlk.JCl-1 && SolnBlk.Grid.BCtypeS[i] == BC_FIXED) {
	  // SOUTH face of cell (i,j+1) is a FIXED boundary.
	  diamond_path = DIAMONDPATH_RIGHT_TRIANGLE_HEATFLUX;
	  Uu = SolnBlk.UoS[i]; Ud = Uu;
	} else if (j == SolnBlk.JCu && SolnBlk.Grid.BCtypeN[i] == BC_FIXED) {
	  // NORTH face of cell (i,j) is a FIXED boundary.
	  diamond_path = DIAMONDPATH_LEFT_TRIANGLE_HEATFLUX;
	  Uu = SolnBlk.UoN[i]; Ud = Uu;
	} else {
	  // NORTH face is either a normal cell or possibly a non-fixed
	  // boundary condition.
	  diamond_path = DIAMONDPATH_QUADRILATERAL_RECONSTRUCTION;
	  Uu = SolnBlk.UnNW(i,j); Ud = SolnBlk.UnNE(i,j);
	}
	Flux = FluxDiamondPath_n(SolnBlk.Grid.xfaceN(i,j),
				 Xl,Ul,Xu,Uu,Xr,Ur,Xd,Ud,
				 SolnBlk.Grid.nfaceN(i,j),
				 diamond_path);
	break;
      };

      // Evaluate cell-averaged solution changes.
      SolnBlk.dUdt[i][j  ][0] += Flux*SolnBlk.Grid.lfaceN(i,j)/SolnBlk.Grid.Cell[i][j].A;
      SolnBlk.dUdt[i][j+1][0] -= Flux*SolnBlk.Grid.lfaceS(i,j+1)/SolnBlk.Grid.Cell[i][j+1].A;

      // Save south and north face boundary flux.
      if (j == SolnBlk.JCl-1) {
	SolnBlk.FluxS[i] = -Flux*SolnBlk.Grid.lfaceS(i,j+1);
      } else if (j == SolnBlk.JCu) {
	SolnBlk.FluxN[i] = Flux*SolnBlk.Grid.lfaceN(i,j);
      }

    }

    SolnBlk.dUdt[i][SolnBlk.JCl-1][0] = ELECTROSTATIC2D_VACUUM;
    SolnBlk.dUdt[i][SolnBlk.JCu+1][0] = ELECTROSTATIC2D_VACUUM;

  }

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
int dUdt_Multistage_Explicit(Electrostatic2D_Quad_Block &SolnBlk,
                             const int i_stage,
			     Electrostatic2D_Input_Parameters &IP) {

  int k_residual, diamond_path;
  double omega;
  Electrostatic2DState Ul, Ur, Uu, Ud, dUdx, dUdy, Flux;
  Vector2D Xl, Xr, Xu, Xd;

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
      if (i_stage == 4) {
	k_residual = 0;
      } else {
	k_residual = i_stage - 1;
      }
    }
    break;
  case TIME_STEPPING_MULTISTAGE_OPTIMAL_SMOOTHING :
    omega = MultiStage_Optimally_Smoothing(i_stage,IP.N_Stage,LIMITER_ONE);
    k_residual = 0;
    break;
  default:
    omega = Runge_Kutta(i_stage,IP.N_Stage);
    k_residual = 0;
    break;
  };

  // Perform the linear reconstruction within each cell of the
  // computational grid for this stage if required.
  if (IP.i_Flux_Reconstruction == VISCOUS_RECONSTRUCTION_ARITHMETIC_AVERAGE || IP.Axisymmetric) {
    switch(IP.i_Reconstruction) {
    case RECONSTRUCTION_GREEN_GAUSS :
      Linear_Reconstruction_GreenGauss(SolnBlk,LIMITER_ONE);
      break;
    case RECONSTRUCTION_LINEAR_LEAST_SQUARES :
      Linear_Reconstruction_LeastSquares(SolnBlk,LIMITER_ONE);
      break;
    };
  }

  // Evaluate the time rate of change of the solution (i.e., the
  // solution residuals) using a second-order limited upwind scheme
  // with a variety of flux functions.

  // Add i-direction (zeta-direction) fluxes.
  for (int j = SolnBlk.JCl-1; j <= SolnBlk.JCu+1; j++) {
    if (i_stage == 1) {
      SolnBlk.Uo[SolnBlk.ICl-1][j] = SolnBlk.U[SolnBlk.ICl-1][j];
      SolnBlk.dUdt[SolnBlk.ICl-1][j][k_residual] = ELECTROSTATIC2D_VACUUM;
    } else {
      SolnBlk.dUdt[SolnBlk.ICl-1][j][k_residual] = ELECTROSTATIC2D_VACUUM;
    }

    for (int i = SolnBlk.ICl-1; i <= SolnBlk.ICu; i++) {
      if (i_stage == 1) {
	SolnBlk.Uo[i+1][j] = SolnBlk.U[i+1][j];
	SolnBlk.dUdt[i+1][j][k_residual] = ELECTROSTATIC2D_VACUUM;
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
	    SolnBlk.dUdt[i+1][j][k_residual] = ELECTROSTATIC2D_VACUUM;
	  }
	  break;
	case TIME_STEPPING_MULTISTAGE_OPTIMAL_SMOOTHING :
	  SolnBlk.dUdt[i+1][j][k_residual] = ELECTROSTATIC2D_VACUUM;
	  break;
	default:
	  SolnBlk.dUdt[i+1][j][k_residual] = ELECTROSTATIC2D_VACUUM;
	  break;
	};
      }

      if (j >= SolnBlk.JCl && j <= SolnBlk.JCu) {

	// Evaluate the cell interface i-direction flux.

	// Select the appropriate nodes and states.
	Ul = SolnBlk.U[i  ][j]; Xl = SolnBlk.Grid.Cell[i  ][j].Xc;
	Ur = SolnBlk.U[i+1][j]; Xr = SolnBlk.Grid.Cell[i+1][j].Xc;

	switch(IP.i_Flux_Reconstruction) {
	case VISCOUS_RECONSTRUCTION_ARITHMETIC_AVERAGE :
	  Flux = FluxArithmetic_n(SolnBlk.Grid.xfaceE(i,j),
				  Xl,Ul,SolnBlk.dUdx[i  ][j],SolnBlk.dUdy[i  ][j],
				  Xr,Ur,SolnBlk.dUdx[i+1][j],SolnBlk.dUdy[i+1][j],
				  SolnBlk.Grid.nfaceE(i,j));
	  break;
	case VISCOUS_RECONSTRUCTION_DIAMOND_PATH :
	  Xu = SolnBlk.Grid.Node[i+1][j+1].X;
	  Xd = SolnBlk.Grid.Node[i+1][j  ].X;
	  if (i == SolnBlk.ICl-1 && SolnBlk.Grid.BCtypeW[j] == BC_FIXED) {
	    // WEST face of cell (i+1,j) is a FIXED boundary.
	    diamond_path = DIAMONDPATH_RIGHT_TRIANGLE_HEATFLUX;
	    Uu = SolnBlk.UoW[j]; Ud = Uu;
	  } else if (i == SolnBlk.ICu && SolnBlk.Grid.BCtypeE[j] == BC_FIXED) {
	    // EAST face of cell (i,j) is a FIXED boundary.
	    diamond_path = DIAMONDPATH_LEFT_TRIANGLE_HEATFLUX;
	    Uu = SolnBlk.UoE[j]; Ud = Uu;
	  } else {
	    // EAST face is either a normal cell or possibly a non-fixed
	    // boundary condition.
	    diamond_path = DIAMONDPATH_QUADRILATERAL_RECONSTRUCTION;
	    Uu = SolnBlk.UnNE(i,j); Ud = SolnBlk.UnSE(i,j);
	  }
	  Flux = FluxDiamondPath_n(SolnBlk.Grid.xfaceE(i,j),
				   Xl,Ul,Xu,Uu,Xr,Ur,Xd,Ud,
				   SolnBlk.Grid.nfaceE(i,j),
				   diamond_path);
	  break;
	};

	// Evaluate cell-averaged solution changes.
	SolnBlk.dUdt[i  ][j][k_residual] += (IP.CFL_Number*SolnBlk.dt[i][j])*
                                            Flux*SolnBlk.Grid.lfaceE(i,j)/SolnBlk.Grid.Cell[i][j].A;
	SolnBlk.dUdt[i+1][j][k_residual] -= (IP.CFL_Number*SolnBlk.dt[i+1][j])*
	                                    Flux*SolnBlk.Grid.lfaceW(i+1,j)/SolnBlk.Grid.Cell[i+1][j].A;

   	// Include axisymmetric source terms if required.
       	if (SolnBlk.Axisymmetric) {
   	  SolnBlk.dUdt[i][j][k_residual] -= (IP.CFL_Number*SolnBlk.dt[i][j])*
  	                                    SolnBlk.U[i][j].S(SolnBlk.Grid.Cell[i][j].Xc,
							      SolnBlk.dUdy[i][j]);
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
      SolnBlk.dUdt[SolnBlk.ICl-1][j][k_residual] = ELECTROSTATIC2D_VACUUM;
      SolnBlk.dUdt[SolnBlk.ICu+1][j][k_residual] = ELECTROSTATIC2D_VACUUM;
    }

  }

  // Add j-direction (eta-direction) fluxes.
  for (int i = SolnBlk.ICl; i <= SolnBlk.ICu; i++) {
    for (int j = SolnBlk.JCl-1; j <= SolnBlk.JCu; j++) {

      // Evaluate the cell interface j-direction flux.

      // Select the appropriate nodes and states.
      Ul = SolnBlk.U[i][j  ]; Xl = SolnBlk.Grid.Cell[i][j  ].Xc;
      Ur = SolnBlk.U[i][j+1]; Xr = SolnBlk.Grid.Cell[i][j+1].Xc;

      switch(IP.i_Flux_Reconstruction) {
      case VISCOUS_RECONSTRUCTION_ARITHMETIC_AVERAGE :
	Flux = FluxArithmetic_n(SolnBlk.Grid.xfaceN(i,j),
				Xl,Ul,SolnBlk.dUdx[i][j  ],SolnBlk.dUdy[i][j  ],
				Xr,Ur,SolnBlk.dUdx[i][j+1],SolnBlk.dUdy[i][j+1],
				SolnBlk.Grid.nfaceN(i,j));
	break;
      case VISCOUS_RECONSTRUCTION_DIAMOND_PATH :
	Xu = SolnBlk.Grid.Node[i  ][j+1].X;
	Xd = SolnBlk.Grid.Node[i+1][j+1].X;
	if (j == SolnBlk.JCl-1 && SolnBlk.Grid.BCtypeS[i] == BC_FIXED) {
	  // SOUTH face of cell (i,j+1) is a FIXED boundary.
	  diamond_path = DIAMONDPATH_RIGHT_TRIANGLE_HEATFLUX;
	  Uu = SolnBlk.UoS[i]; Ud = Uu;
	} else if (j == SolnBlk.JCu && SolnBlk.Grid.BCtypeN[i] == BC_FIXED) {
	  // NORTH face of cell (i,j) is a FIXED boundary.
	  diamond_path = DIAMONDPATH_LEFT_TRIANGLE_HEATFLUX;
	  Uu = SolnBlk.UoN[i]; Ud = Uu;
	} else {
	  // NORTH face is either a normal cell or possibly a non-fixed
	  // boundary condition.
	  diamond_path = DIAMONDPATH_QUADRILATERAL_RECONSTRUCTION;
	  Uu = SolnBlk.UnNW(i,j); Ud = SolnBlk.UnNE(i,j);
	}
	Flux = FluxDiamondPath_n(SolnBlk.Grid.xfaceN(i,j),
				 Xl,Ul,Xu,Uu,Xr,Ur,Xd,Ud,
				 SolnBlk.Grid.nfaceN(i,j),
				 diamond_path);
	break;
      };

      // Evaluate cell-averaged solution changes.
      SolnBlk.dUdt[i][j  ][k_residual] += (IP.CFL_Number*SolnBlk.dt[i][j  ])*
                                          Flux*SolnBlk.Grid.lfaceN(i,j)/SolnBlk.Grid.Cell[i][j].A;
      SolnBlk.dUdt[i][j+1][k_residual] -= (IP.CFL_Number*SolnBlk.dt[i][j+1])*
                                          Flux*SolnBlk.Grid.lfaceS(i,j+1)/SolnBlk.Grid.Cell[i][j+1].A;

      // Save south and north face boundary flux.
      if (j == SolnBlk.JCl-1) {
 	SolnBlk.FluxS[i] = -Flux*SolnBlk.Grid.lfaceS(i,j+1);
      } else if (j == SolnBlk.JCu) {
	SolnBlk.FluxN[i] = Flux*SolnBlk.Grid.lfaceN(i,j);
      }

    }

    SolnBlk.dUdt[i][SolnBlk.JCl-1][k_residual] = ELECTROSTATIC2D_VACUUM;
    SolnBlk.dUdt[i][SolnBlk.JCu+1][k_residual] = ELECTROSTATIC2D_VACUUM;

  }

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
int Update_Solution_Multistage_Explicit(Electrostatic2D_Quad_Block &SolnBlk,
					const int i_stage,
					Electrostatic2D_Input_Parameters &IP) {

  int k_residual;
  double omega;

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
    omega = MultiStage_Optimally_Smoothing(i_stage,IP.N_Stage,LIMITER_ONE);
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

      SolnBlk.U[i][j].V = SolnBlk.Uo[i][j].V + omega*SolnBlk.dUdt[i][j][k_residual].V;

    }
  }

  // Solution successfully updated.
  return 0;

}

/**********************************************************************
 * Routine: Determine_Electric_Field                                  *
 **********************************************************************/
int Determine_Electric_Field(Electrostatic2D_Quad_Block &SolnBlk,
			     Electrostatic2D_Input_Parameters &IP) {

  switch(IP.i_Reconstruction) {
  case RECONSTRUCTION_GREEN_GAUSS :
    Linear_Reconstruction_GreenGauss(SolnBlk,LIMITER_ONE);
    break;
  case RECONSTRUCTION_LINEAR_LEAST_SQUARES :
    Linear_Reconstruction_LeastSquares(SolnBlk,LIMITER_ONE);
    break;
  };

  for (int j = SolnBlk.JCl-SolnBlk.Nghost+1; j <= SolnBlk.JCu+SolnBlk.Nghost-1; j++) {
    for (int i = SolnBlk.ICl-SolnBlk.Nghost+1; i <= SolnBlk.ICu+SolnBlk.Nghost-1; i++) {
      SolnBlk.U[i][j].E.x = -SolnBlk.dUdx[i][j].V;
      SolnBlk.U[i][j].E.y = -SolnBlk.dUdy[i][j].V;
    }
  }

  return 0;

}
