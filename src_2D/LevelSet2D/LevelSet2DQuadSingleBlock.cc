/**********************************************************************
 * LevelSet2DQuadSingleBlock.cc                                       *
 *                                                                    *
 * Single-block versions of subroutines for 2D Level Set multi-block  *
 * quadrilateral mesh solution classes.                               *
 *                                                                    *
 **********************************************************************/

// Include 2D LevelSet quadrilateral mesh solution header file.

#ifndef _LEVELSET2D_QUAD_INCLUDED
#include "LevelSet2DQuad.h"
#endif // _LEVELSET2D_QUAD_INCLUDED

/**********************************************************************
 * LevelSet2D_Quad_Block -- Single Block External Subroutines.        *
 **********************************************************************/

/**********************************************************************
 * Routine: Broadcast_Solution_Block                                  *
 *                                                                    *
 * Broadcast quadrilateral solution block to all processors involved  *
 * in the calculation from the primary processor using the MPI        *
 * broadcast routine.                                                 *
 *                                                                    *
 **********************************************************************/
void Broadcast_Solution_Block(LevelSet2D_Quad_Block &SolnBlk) {

#ifdef _MPI_VERSION
  int ni, nj, ng, block_allocated, buffer_size;
  double *buffer;
  int *i_buffer;

  // Broadcast the number of cells in each direction.
  if (CFFC_Primary_MPI_Processor()) {
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
  if (!CFFC_Primary_MPI_Processor()) {
    if (SolnBlk.NCi != ni || SolnBlk.NCj != nj || SolnBlk.Nghost != ng) { 
      if (SolnBlk.U != NULL) SolnBlk.deallocate();
      if (block_allocated) SolnBlk.allocate(ni-2*ng,nj-2*ng,ng); 
    }
  }

  // Broadcast the grid.
  Broadcast_Quad_Block(SolnBlk.Grid);

  // Broadcast the solution state variables.
  if (block_allocated) {
    ni = (SolnBlk.ICu+SolnBlk.Nghost) - (SolnBlk.ICl-SolnBlk.Nghost) + 1;
    nj = (SolnBlk.JCu+SolnBlk.Nghost) - (SolnBlk.JCl-SolnBlk.Nghost) + 1;
    buffer = new double[NUM_VAR_LEVELSET2D*ni*nj];
    if (CFFC_Primary_MPI_Processor()) {
      buffer_size = 0;
      for (int j = SolnBlk.JCl-SolnBlk.Nghost; j <= SolnBlk.JCu+SolnBlk.Nghost; j++) {
	for (int i = SolnBlk.ICl-SolnBlk.Nghost; i <= SolnBlk.ICu+SolnBlk.Nghost; i++) {
	  buffer[buffer_size  ] = SolnBlk.U[i][j].psi;
	  buffer[buffer_size+1] = SolnBlk.U[i][j].F;
	  buffer[buffer_size+2] = SolnBlk.U[i][j].V.x;
	  buffer[buffer_size+3] = SolnBlk.U[i][j].V.y;
	  buffer_size = buffer_size + NUM_VAR_LEVELSET2D;
	}
      }
    }

    buffer_size = NUM_VAR_LEVELSET2D*ni*nj;
    MPI::COMM_WORLD.Bcast(buffer,buffer_size,MPI::DOUBLE,0);

    if (!CFFC_Primary_MPI_Processor()) {
      buffer_size = 0;
      for (int j = SolnBlk.JCl-SolnBlk.Nghost; j <= SolnBlk.JCu+SolnBlk.Nghost; j++) {
	for (int i = SolnBlk.ICl-SolnBlk.Nghost; i <= SolnBlk.ICu+SolnBlk.Nghost; i++) {
	  SolnBlk.U[i][j].psi = buffer[buffer_size  ];
	  SolnBlk.U[i][j].F   = buffer[buffer_size+1];
	  SolnBlk.U[i][j].V.x = buffer[buffer_size+2];
	  SolnBlk.U[i][j].V.y = buffer[buffer_size+3];
	  buffer_size = buffer_size + NUM_VAR_LEVELSET2D;
	}
      }
    }

    delete []buffer; 
    buffer = NULL;

  }

  // Broadcast interface list.
  Broadcast_Interface_List(SolnBlk.Interface_List);

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
void Broadcast_Solution_Block(LevelSet2D_Quad_Block &SolnBlk,
                              MPI::Intracomm &Communicator,
                              const int Source_CPU) {

  int Source_Rank = 0;
  int ni, nj, ng, block_allocated, buffer_size;
  double *buffer;
  int *i_buffer;
  
  // Broadcast the number of cells in each direction.
  if (CFFC_MPI::This_Processor_Number == Source_CPU) {
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
  if (!(CFFC_MPI::This_Processor_Number == Source_CPU)) {
    if (SolnBlk.NCi != ni || SolnBlk.NCj != nj || SolnBlk.Nghost != ng) { 
      if (SolnBlk.U != NULL) SolnBlk.deallocate();
      if (block_allocated) SolnBlk.allocate(ni-2*ng,nj-2*ng,ng); 
    }
  }

  // Broadcast the grid.
  Broadcast_Quad_Block(SolnBlk.Grid,Communicator,Source_CPU);

  // Broadcast the solution state variables.
  if (block_allocated) {
    ni = (SolnBlk.ICu+SolnBlk.Nghost) - (SolnBlk.ICl-SolnBlk.Nghost) + 1;
    nj = (SolnBlk.JCu+SolnBlk.Nghost) - (SolnBlk.JCl-SolnBlk.Nghost) + 1;
    buffer = new double[NUM_VAR_LEVELSET2D*ni*nj];

    if (CFFC_MPI::This_Processor_Number == Source_CPU) {
      buffer_size = 0;
      for (int j = SolnBlk.JCl-SolnBlk.Nghost; j <= SolnBlk.JCu+SolnBlk.Nghost; j++) {
	for (int i = SolnBlk.ICl-SolnBlk.Nghost; i <= SolnBlk.ICu+SolnBlk.Nghost; i++) {
	  buffer[buffer_size  ] = SolnBlk.U[i][j].psi;
	  buffer[buffer_size+1] = SolnBlk.U[i][j].F;
	  buffer[buffer_size+2] = SolnBlk.U[i][j].V.x;
	  buffer[buffer_size+3] = SolnBlk.U[i][j].V.y;
	  buffer_size = buffer_size + NUM_VAR_LEVELSET2D;
	}
      }
    }

    buffer_size = NUM_VAR_LEVELSET2D*ni*nj;
    Communicator.Bcast(buffer,buffer_size,MPI::DOUBLE,Source_Rank);

    if (!(CFFC_MPI::This_Processor_Number == Source_CPU)) {
      buffer_size = 0;
      for (int j = SolnBlk.JCl-SolnBlk.Nghost; j <= SolnBlk.JCu+SolnBlk.Nghost; j++) {
	for (int i = SolnBlk.ICl-SolnBlk.Nghost; i <= SolnBlk.ICu+SolnBlk.Nghost; i++) {
	  SolnBlk.U[i][j].psi = buffer[buffer_size  ];
	  SolnBlk.U[i][j].F   = buffer[buffer_size+1];
	  SolnBlk.U[i][j].V.x = buffer[buffer_size+2];
	  SolnBlk.U[i][j].V.y = buffer[buffer_size+3];
 	  buffer_size = buffer_size + NUM_VAR_LEVELSET2D;
	}
      }
    }

    delete []buffer; 
    buffer = NULL;

  }

  // Broadcast interface list.
  Broadcast_Interface_List(SolnBlk.Interface_List,Communicator,Source_CPU);

}
#endif

/**********************************************************************
 * Routine: Copy_Solution_Block                                       *
 *                                                                    *
 * Copies the solution information of quadrilateral solution block    *
 * SolnBlk2 to SolnBlk1.                                              *
 *                                                                    *
 **********************************************************************/
void Copy_Solution_Block(LevelSet2D_Quad_Block &SolnBlk1,
                         LevelSet2D_Quad_Block &SolnBlk2) {

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

  // Copy the grid of the second solution block to the first solution 
  // block.
  Copy_Quad_Block(SolnBlk1.Grid,SolnBlk2.Grid);

  // Copy the solution information from SolnBlk2 to SolnBlk1.
  if (SolnBlk2.U != NULL) {
    for (int j = SolnBlk1.JCl-SolnBlk1.Nghost; j <= SolnBlk1.JCu+SolnBlk1.Nghost; j++) {
      for (int i = SolnBlk1.ICl-SolnBlk1.Nghost; i <= SolnBlk1.ICu+SolnBlk1.Nghost; i++) {
	SolnBlk1.U[i][j] = SolnBlk2.U[i][j];
	for (int k = 0; k < NUMBER_OF_RESIDUAL_VECTORS; k++) {
	  SolnBlk1.dUdt[i][j][k] = SolnBlk2.dUdt[i][j][k];
	}
	SolnBlk1.dUdx[i][j] = SolnBlk2.dUdx[i][j];
	SolnBlk1.dUdy[i][j] = SolnBlk2.dUdy[i][j];
	SolnBlk1.phi[i][j]  = SolnBlk2.phi[i][j];
	SolnBlk1.Uo[i][j]   = SolnBlk2.Uo[i][j];
	SolnBlk1.dt[i][j]   = SolnBlk2.dt[i][j];
      }
    }
  }

  // Copy the interface list.
  SolnBlk1.Interface_List.Copy(SolnBlk2.Interface_List);

}

/**********************************************************************
 * Routine: Prolong_Solution_Block                                    *
 *                                                                    *
 * Prolongs the solution information of one of the specified sectors  *
 * of the original quadrilateral solution block SolnBlk_Original to   *
 * the refined solution block SolnBlk_Fine.                           *
 *                                                                    *
 **********************************************************************/
int Prolong_Solution_Block(LevelSet2D_Quad_Block &SolnBlk_Fine,
			   LevelSet2D_Quad_Block &SolnBlk_Original,
			   const int Sector) {

  int error_flag;
  int i_min, i_max, j_min, j_max, mesh_refinement_permitted;
  Vector2D dX;

  // Allocate (re-allocate) memory for the solution of the refined
  // quadrilateral solution block as necessary.
  if ((SolnBlk_Original.NCi-2*SolnBlk_Original.Nghost-2*((SolnBlk_Original.NCi-2*SolnBlk_Original.Nghost)/2) != 0) || 
      (SolnBlk_Original.NCj-2*SolnBlk_Original.Nghost-2*((SolnBlk_Original.NCj-2*SolnBlk_Original.Nghost)/2) != 0) ||
      (SolnBlk_Original.NCi-2*SolnBlk_Original.Nghost < 2*SolnBlk_Original.Nghost) ||
      (SolnBlk_Original.NCj-2*SolnBlk_Original.Nghost < 2*SolnBlk_Original.Nghost) ||
      (SolnBlk_Original.Grid.Node == NULL) ) {
    mesh_refinement_permitted = 0;
  } else {
    mesh_refinement_permitted = 1;
    if (SolnBlk_Fine.NCi != SolnBlk_Original.NCi || 
	SolnBlk_Fine.NCj != SolnBlk_Original.NCj) {
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
      for (int i = i_min; i <= i_max; i++) {

	Linear_Reconstruction_LeastSquares(SolnBlk_Original,i,j,LIMITER_ONE);

	dX = SolnBlk_Fine.Grid.Cell[2*(i-i_min)+SolnBlk_Fine.ICl  ][2*(j-j_min)+SolnBlk_Fine.JCl  ].Xc -
	     SolnBlk_Original.Grid.Cell[i][j].Xc;
	SolnBlk_Fine.U[2*(i-i_min)+SolnBlk_Fine.ICl  ][2*(j-j_min)+SolnBlk_Fine.JCl  ]
	  = SolnBlk_Original.U[i][j] + SolnBlk_Original.dUdx[i][j]*dX.x + SolnBlk_Original.dUdy[i][j]*dX.y;

	dX = SolnBlk_Fine.Grid.Cell[2*(i-i_min)+SolnBlk_Fine.ICl+1][2*(j-j_min)+SolnBlk_Fine.JCl  ].Xc -
	     SolnBlk_Original.Grid.Cell[i][j].Xc;
	SolnBlk_Fine.U[2*(i-i_min)+SolnBlk_Fine.ICl+1][2*(j-j_min)+SolnBlk_Fine.JCl  ]
	  = SolnBlk_Original.U[i][j] + SolnBlk_Original.dUdx[i][j]*dX.x + SolnBlk_Original.dUdy[i][j]*dX.y;

	dX = SolnBlk_Fine.Grid.Cell[2*(i-i_min)+SolnBlk_Fine.ICl  ][2*(j-j_min)+SolnBlk_Fine.JCl+1].Xc -
	     SolnBlk_Original.Grid.Cell[i][j].Xc;
	SolnBlk_Fine.U[2*(i-i_min)+SolnBlk_Fine.ICl  ][2*(j-j_min)+SolnBlk_Fine.JCl+1]
	  = SolnBlk_Original.U[i][j] + SolnBlk_Original.dUdx[i][j]*dX.x + SolnBlk_Original.dUdy[i][j]*dX.y;

	dX = SolnBlk_Fine.Grid.Cell[2*(i-i_min)+SolnBlk_Fine.ICl+1][2*(j-j_min)+SolnBlk_Fine.JCl+1].Xc -
	     SolnBlk_Original.Grid.Cell[i][j].Xc;
	SolnBlk_Fine.U[2*(i-i_min)+SolnBlk_Fine.ICl+1][2*(j-j_min)+SolnBlk_Fine.JCl+1]
	  = SolnBlk_Original.U[i][j] + SolnBlk_Original.dUdx[i][j]*dX.x + SolnBlk_Original.dUdy[i][j]*dX.y;

      }
    }

    // Prolong the interface list.
    SolnBlk_Fine.Interface_List.Copy(SolnBlk_Original.Interface_List);

  }

  // Prolong solution block was successful.
  return 0;

}

/**********************************************************************
 * Routine: Restrict_Solution_Block                                   *
 *                                                                    *
 * Restricts the solution information of four original fine           *
 * quadrilateral solution blocks to the coarse solution block         *
 * SolnBlk_Coarse.                                                    *
 *                                                                    *
 **********************************************************************/
int Restrict_Solution_Block(LevelSet2D_Quad_Block &SolnBlk_Coarse,
			    LevelSet2D_Quad_Block &SolnBlk_Original_SW,
			    LevelSet2D_Quad_Block &SolnBlk_Original_SE,
			    LevelSet2D_Quad_Block &SolnBlk_Original_NW,
			    LevelSet2D_Quad_Block &SolnBlk_Original_NE) {

  int error_flag;
  int i, j, i_coarse, j_coarse, mesh_coarsening_permitted;
 
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

    // Restrict the solution information from the four original solution
    // blocks to the newly coarsened solution block.

    // South-West corner fine block:
    for (int j = SolnBlk_Original_SW.JCl; j <= SolnBlk_Original_SW.JCu; j += 2) {
      for (int i = SolnBlk_Original_SW.ICl; i <= SolnBlk_Original_SW.ICu; i += 2) {
	i_coarse = (i-SolnBlk_Original_SW.ICl)/2 + SolnBlk_Coarse.ICl;
	j_coarse = (j-SolnBlk_Original_SW.JCl)/2 + SolnBlk_Coarse.JCl;
	SolnBlk_Coarse.U[i_coarse][j_coarse] = SolnBlk_Original_SW.UnNE(i,j);
      }
    }

    // South-East corner fine block:
    for (int j = SolnBlk_Original_SE.JCl; j <= SolnBlk_Original_SE.JCu; j += 2) {
      for (int i = SolnBlk_Original_SE.ICl; i <= SolnBlk_Original_SE.ICu; i += 2) {
	i_coarse = (i-SolnBlk_Original_SE.ICl)/2 + 
	           (SolnBlk_Coarse.ICu-SolnBlk_Coarse.ICl+1)/2 + SolnBlk_Coarse.ICl;
	j_coarse = (j-SolnBlk_Original_SE.JCl)/2 + SolnBlk_Coarse.JCl;
	SolnBlk_Coarse.U[i_coarse][j_coarse] = SolnBlk_Original_SE.UnNE(i,j);
      }
    }

    // North-West corner fine block:
    for (int j = SolnBlk_Original_NW.JCl; j <= SolnBlk_Original_NW.JCu; j += 2) {
      for (int i = SolnBlk_Original_NW.ICl; i <= SolnBlk_Original_NW.ICu; i += 2) {
	i_coarse = (i-SolnBlk_Original_NW.ICl)/2 + SolnBlk_Coarse.ICl;
	j_coarse = (j-SolnBlk_Original_NW.JCl)/2 + 
	           (SolnBlk_Coarse.JCu-SolnBlk_Coarse.JCl+1)/2 + SolnBlk_Coarse.JCl;
	SolnBlk_Coarse.U[i_coarse][j_coarse] = SolnBlk_Original_NW.UnNE(i,j);
      }
    }

    // North-east corner fine block:
    for (int j = SolnBlk_Original_NE.JCl; j <= SolnBlk_Original_NE.JCu; j += 2) {
      for (int i = SolnBlk_Original_NE.ICl; i <= SolnBlk_Original_NE.ICu; i += 2) {
	i_coarse = (i-SolnBlk_Original_NE.ICl)/2 + 
	           (SolnBlk_Coarse.ICu-SolnBlk_Coarse.ICl+1)/2 + SolnBlk_Coarse.ICl;
	j_coarse = (j-SolnBlk_Original_NE.JCl)/2 +
	           (SolnBlk_Coarse.JCu-SolnBlk_Coarse.JCl+1)/2 + SolnBlk_Coarse.JCl;
	SolnBlk_Coarse.U[i_coarse][j_coarse] = SolnBlk_Original_NE.UnNE(i,j);
      }
    }

    // Restrict the interface list.
    SolnBlk_Coarse.Interface_List.Copy(SolnBlk_Original_SW.Interface_List);

  }

  // Restrict solution block was successful.
  return 0;

}

/**********************************************************************
 * Routine: Construct_Bulk_Flow_Field                                 *
 *                                                                    *
 * Constructs the bulk flow field on the specified solution block.    *
 *                                                                    *
 **********************************************************************/
int Construct_Bulk_Flow_Field(LevelSet2D_Quad_Block &SolnBlk,
			      LevelSet2D_Input_Parameters &IP) {

  // Construct the bulk flow field.
  if (IP.i_BulkFlowField_Type == INTERFACE_BULKFLOWFIELD_NONE) {
    for (int j = SolnBlk.JCl-SolnBlk.Nghost; j <= SolnBlk.JCu+SolnBlk.Nghost; j++)
      for (int i = SolnBlk.ICl-SolnBlk.Nghost; i <= SolnBlk.ICu+SolnBlk.Nghost; i++)
	SolnBlk.U[i][j].V = Vector2D_ZERO;
  } else if (IP.i_BulkFlowField_Type == INTERFACE_BULKFLOWFIELD_UNIFORM) {
    for (int j = SolnBlk.JCl-SolnBlk.Nghost; j <= SolnBlk.JCu+SolnBlk.Nghost; j++)
      for (int i = SolnBlk.ICl-SolnBlk.Nghost; i <= SolnBlk.ICu+SolnBlk.Nghost; i++) {
	SolnBlk.U[i][j].V = IP.V;
      }
  } else if (IP.i_BulkFlowField_Type == INTERFACE_BULKFLOWFIELD_SWIRL) {
    if (IP.Interface_IP.Component_List[1].Type == INTERFACE_ZALESAK) {
      for (int j = SolnBlk.JCl-SolnBlk.Nghost; j <= SolnBlk.JCu+SolnBlk.Nghost; j++) {
	for (int i = SolnBlk.ICl-SolnBlk.Nghost; i <= SolnBlk.ICu+SolnBlk.Nghost; i++) {
 	  SolnBlk.U[i][j].V.x = PI*(0.50 - SolnBlk.Grid.Cell[i][j].Xc.y);
 	  SolnBlk.U[i][j].V.y = PI*(SolnBlk.Grid.Cell[i][j].Xc.x - 0.50);
	}
      }
    } else {
      for (int j = SolnBlk.JCl-SolnBlk.Nghost; j <= SolnBlk.JCu+SolnBlk.Nghost; j++) {
	for (int i = SolnBlk.ICl-SolnBlk.Nghost; i <= SolnBlk.ICu+SolnBlk.Nghost; i++) {
	  SolnBlk.U[i][j].V.x = -IP.V.x*SolnBlk.Grid.Cell[i][j].Xc.y;
	  SolnBlk.U[i][j].V.y =  IP.V.x*SolnBlk.Grid.Cell[i][j].Xc.x;
	}
      }
    }
  } else {
    return 1;
  }

  // Bulk/convective flow-field successfully created.
  return 0;

}

/**********************************************************************
 * Routine: BCs                                                       *
 *                                                                    *
 * Apply boundary conditions at boundaries of the specified           *
 * quadrilateral solution block.                                      *
 *                                                                    *
 **********************************************************************/
void BCs(LevelSet2D_Quad_Block &SolnBlk) {

  Vector2D dX;

  // WEST and EAST boundary conditions.
  for (int j = SolnBlk.JCl-SolnBlk.Nghost; j <= SolnBlk.JCu+SolnBlk.Nghost; j++) {
    // WEST boundary.
    switch(SolnBlk.Grid.BCtypeW[j]) {
    case BC_NONE :
      break;
    case BC_CONSTANT_EXTRAPOLATION :
      for (int ghost = 1; ghost <= SolnBlk.Nghost; ghost++) {
	SolnBlk.U[SolnBlk.ICl-ghost][j] = SolnBlk.U[SolnBlk.ICl][j];
      }
      break;
    case BC_LINEAR_EXTRAPOLATION :
      //Linear_Reconstruction_LeastSquares_2(SolnBlk,SolnBlk.ICl,j,LIMITER_ONE);
      Linear_Reconstruction_LeastSquares(SolnBlk,SolnBlk.ICl,j,LIMITER_ONE);
      for (int ghost = 1; ghost <= SolnBlk.Nghost; ghost++) {
	dX = SolnBlk.Grid.Cell[SolnBlk.ICl-ghost][j].Xc - SolnBlk.Grid.Cell[SolnBlk.ICl][j].Xc;
	SolnBlk.U[SolnBlk.ICl-ghost][j].psi = SolnBlk.U[SolnBlk.ICl][j].psi +
	                                      ((SolnBlk.phi[SolnBlk.ICl][j].psi*SolnBlk.dUdx[SolnBlk.ICl][j].psi)*dX.x +
					       (SolnBlk.phi[SolnBlk.ICl][j].psi*SolnBlk.dUdy[SolnBlk.ICl][j].psi)*dX.y);
	SolnBlk.U[SolnBlk.ICl-ghost][j].F = SolnBlk.U[SolnBlk.ICl][j].F;
	SolnBlk.U[SolnBlk.ICl-ghost][j].V = SolnBlk.U[SolnBlk.ICl][j].V;
      }
      break;
    default:
      for (int ghost = 1; ghost <= SolnBlk.Nghost; ghost++) {
	SolnBlk.U[SolnBlk.ICl-ghost][j] = SolnBlk.U[SolnBlk.ICl][j];
      }
      break;
    };

    // EAST boundary.
    switch(SolnBlk.Grid.BCtypeE[j]) {
    case BC_NONE :
      break;
    case BC_CONSTANT_EXTRAPOLATION :
      for (int ghost = 1; ghost <= SolnBlk.Nghost; ghost++) {
	SolnBlk.U[SolnBlk.ICu+ghost][j] = SolnBlk.U[SolnBlk.ICu][j];
      }
      break;
    case BC_LINEAR_EXTRAPOLATION :
      //Linear_Reconstruction_LeastSquares_2(SolnBlk,SolnBlk.ICu,j,LIMITER_ONE);
      Linear_Reconstruction_LeastSquares(SolnBlk,SolnBlk.ICu,j,LIMITER_ONE);
      for (int ghost = 1; ghost <= SolnBlk.Nghost; ghost++){
	dX = SolnBlk.Grid.Cell[SolnBlk.ICu+ghost][j].Xc - SolnBlk.Grid.Cell[SolnBlk.ICu][j].Xc;
	SolnBlk.U[SolnBlk.ICu+ghost][j].psi = SolnBlk.U[SolnBlk.ICu][j].psi + 
	  ((SolnBlk.phi[SolnBlk.ICu][j].psi*SolnBlk.dUdx[SolnBlk.ICu][j].psi)*dX.x +
	   (SolnBlk.phi[SolnBlk.ICu][j].psi*SolnBlk.dUdy[SolnBlk.ICu][j].psi)*dX.y);
	SolnBlk.U[SolnBlk.ICu+ghost][j].F = SolnBlk.U[SolnBlk.ICu][j].F;
	SolnBlk.U[SolnBlk.ICu+ghost][j].V = SolnBlk.U[SolnBlk.ICu][j].V;
      }
      break;
    default:
      for (int ghost=1; ghost <= SolnBlk.Nghost; ghost++) {
	SolnBlk.U[SolnBlk.ICu+ghost][j] = SolnBlk.U[SolnBlk.ICu][j];
      }
      break;
    };

  }

  // NORTH and SOUTH boundary conditions.
  for (int i = SolnBlk.ICl-SolnBlk.Nghost; i <= SolnBlk.ICu+SolnBlk.Nghost; i++) {

    // SOUTH boundary.
    switch(SolnBlk.Grid.BCtypeS[i]) {
    case BC_NONE :
      break;
    case BC_CONSTANT_EXTRAPOLATION :
      for (int ghost = 1; ghost <= SolnBlk.Nghost; ghost++) {
	SolnBlk.U[i][SolnBlk.JCl-ghost] = SolnBlk.U[i][SolnBlk.JCl];
      }
      break;
    case BC_LINEAR_EXTRAPOLATION :
      //Linear_Reconstruction_LeastSquares_2(SolnBlk,i,SolnBlk.JCl,LIMITER_ONE);
      Linear_Reconstruction_LeastSquares(SolnBlk,i,SolnBlk.JCl,LIMITER_ONE);
      for (int ghost = 1; ghost <= SolnBlk.Nghost; ghost++) {
	dX = SolnBlk.Grid.Cell[i][SolnBlk.JCl-ghost].Xc - SolnBlk.Grid.Cell[i][SolnBlk.JCl].Xc;
	SolnBlk.U[i][SolnBlk.JCl-ghost].psi = SolnBlk.U[i][SolnBlk.JCl].psi +
	                                      ((SolnBlk.phi[i][SolnBlk.JCl].psi*SolnBlk.dUdx[i][SolnBlk.JCl].psi)*dX.x +
					       (SolnBlk.phi[i][SolnBlk.JCl].psi*SolnBlk.dUdy[i][SolnBlk.JCl].psi)*dX.y);
	SolnBlk.U[i][SolnBlk.JCl-ghost].F = SolnBlk.U[i][SolnBlk.JCl].F;
	SolnBlk.U[i][SolnBlk.JCl-ghost].V = SolnBlk.U[i][SolnBlk.JCl].V;
      }
      break;
    default:
      for (int ghost = 1; ghost <= SolnBlk.Nghost; ghost++) {
	SolnBlk.U[i][SolnBlk.JCl-ghost] = SolnBlk.U[i][SolnBlk.JCl];
      }
      break;
    };

    // NORTH boundary.
    switch(SolnBlk.Grid.BCtypeN[i]) {
    case BC_NONE :
      break;
    case BC_CONSTANT_EXTRAPOLATION :
      for (int ghost = 1; ghost <= SolnBlk.Nghost; ghost++) {
	SolnBlk.U[i][SolnBlk.JCu+ghost] = SolnBlk.U[i][SolnBlk.JCu];
      }
      break;
    case BC_LINEAR_EXTRAPOLATION :
      //Linear_Reconstruction_LeastSquares_2(SolnBlk,i,SolnBlk.JCu,LIMITER_ONE);
      Linear_Reconstruction_LeastSquares(SolnBlk,i,SolnBlk.JCu,LIMITER_ONE);
      for( int ghost = 1; ghost <= SolnBlk.Nghost; ghost++) {
	dX = SolnBlk.Grid.Cell[i][SolnBlk.JCu+ghost].Xc - SolnBlk.Grid.Cell[i][SolnBlk.JCu].Xc;
	SolnBlk.U[i][SolnBlk.JCu+ghost].psi = SolnBlk.U[i][SolnBlk.JCu].psi + 
	                                      ((SolnBlk.phi[i][SolnBlk.JCu].psi*SolnBlk.dUdx[i][SolnBlk.JCu].psi)*dX.x +
					       (SolnBlk.phi[i][SolnBlk.JCu].psi*SolnBlk.dUdy[i][SolnBlk.JCu].psi)*dX.y);
	SolnBlk.U[i][SolnBlk.JCu+ghost].F = SolnBlk.U[i][SolnBlk.JCu].F;
	SolnBlk.U[i][SolnBlk.JCu+ghost].V = SolnBlk.U[i][SolnBlk.JCu].V;
      }
      break;
    default:
      for (int ghost = 1; ghost <= SolnBlk.Nghost; ghost++) {
	SolnBlk.U[i][SolnBlk.JCu+ghost] = SolnBlk.U[i][SolnBlk.JCu];
      }
      break;
    };

  }

}

/**********************************************************************
 * Routine: Set_Global_TimeStep                                       *
 *                                                                    *
 * Assigns global time step to specified solution block for           *
 * time-accurate calculations.                                        *
 *                                                                    *
 **********************************************************************/
void Set_Global_TimeStep(LevelSet2D_Quad_Block &SolnBlk,
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
double L1_Norm_Residual(LevelSet2D_Quad_Block &SolnBlk,
			const int &n) {

  double l1_norm = ZERO;

  for (int j = SolnBlk.JCl; j <= SolnBlk.JCu; j++) {
    for (int i = SolnBlk.ICl; i <= SolnBlk.ICu; i++) {
      l1_norm += fabs(SolnBlk.dUdt[i][j][0][n]);
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
double L2_Norm_Residual(LevelSet2D_Quad_Block &SolnBlk,
			const int &n) {

  double l2_norm = ZERO;

  for (int j = SolnBlk.JCl; j <= SolnBlk.JCu; j++) {
    for (int i = SolnBlk.ICl; i <= SolnBlk.ICu; i++) {
      l2_norm += sqr(fabs(SolnBlk.dUdt[i][j][0][n]));
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
double Max_Norm_Residual(LevelSet2D_Quad_Block &SolnBlk,
			 const int &n) {

  double dUdt, max_norm = ZERO;

  for (int j = SolnBlk.JCl; j <= SolnBlk.JCu; j++) {
    for (int i = SolnBlk.ICl; i <= SolnBlk.ICu; i++) {
      dUdt = fabs(SolnBlk.dUdt[i][j][0][n]);
      max_norm = max(max_norm,dUdt);
    }
  }

  return max_norm;
  
}

/**********************************************************************
 * Routine: Store_Initial_Eikonal_Solution                            *
 *                                                                    *
 * This routines stores a copy of the solution prior to solving the   *
 * Eikonal equation.                                                  *
 *                                                                    *
 **********************************************************************/
int Store_Initial_Eikonal_Solution(LevelSet2D_Quad_Block &SolnBlk) {

  // Store the initial signed-distance function.
  for (int j = SolnBlk.JCl-SolnBlk.Nghost; j <= SolnBlk.JCu+SolnBlk.Nghost; j++) {
    for (int i = SolnBlk.ICl-SolnBlk.Nghost; i <= SolnBlk.ICu+SolnBlk.Nghost; i++) {
      SolnBlk.Uoo[i][j].psi = SolnBlk.U[i][j].psi;
    }
  }

  // Initial signed-distance function successfully stored.
  return 0;
}

/**********************************************************************
 * Routine: Calculate_Sign_Function                                   *
 *                                                                    *
 * This routines calculates and stores the sign function of the       *
 * signed distance function of the specified quadrilateral solution   *
 * block.  Required for the solution of the Eikonal equation for      *
 * redistancing the distance function and the scalar extension        *
 * equation.                                                          *
 *                                                                    *
 **********************************************************************/
int Calculate_Sign_Function(LevelSet2D_Quad_Block &SolnBlk,
			    LevelSet2D_Input_Parameters &IP) {

  double dx, dU;

  dx = max(TOLER,min(fabs(SolnBlk.Grid.Cell[2][2].Xc.x-SolnBlk.Grid.Cell[1][2].Xc.x),
		     fabs(SolnBlk.Grid.Cell[2][2].Xc.y-SolnBlk.Grid.Cell[2][1].Xc.y)));

  switch(IP.i_Eikonal_Sign_Function) {
  case EIKONAL_SIGN_FUNCTION_DISCRETE :
    for (int j = SolnBlk.JCl-SolnBlk.Nghost; j <= SolnBlk.JCu+SolnBlk.Nghost; j++) {
      for (int i = SolnBlk.ICl-SolnBlk.Nghost; i <= SolnBlk.ICu+SolnBlk.Nghost; i++) {
	if (SolnBlk.U[i][j].psi > NANO) {
	  SolnBlk.sign[i][j] = ONE;
	} else if (SolnBlk.U[i][j].psi < -NANO) {
	  SolnBlk.sign[i][j] = -ONE;
	} else {
	  SolnBlk.sign[i][j] = ZERO;
	}
      }
    }
    break;
  case EIKONAL_SIGN_FUNCTION_SMEARED :
    // Sussman, Smereka, Osher. J. comput. Phys. 114, 146-159 (1994).
    for (int j = SolnBlk.JCl-SolnBlk.Nghost; j <= SolnBlk.JCu+SolnBlk.Nghost; j++) {
      for (int i = SolnBlk.ICl-SolnBlk.Nghost; i <= SolnBlk.ICu+SolnBlk.Nghost; i++) {
	SolnBlk.sign[i][j] = SolnBlk.U[i][j].psi/sqrt(sqr(SolnBlk.U[i][j].psi) + sqr(dx));
      }
    }
    break;
  case EIKONAL_SIGN_FUNCTION_SMEARED_MACDONALD :
    // Colin Macdonald and Steven J. Ruuth.
    for (int j = SolnBlk.JCl-SolnBlk.Nghost; j <= SolnBlk.JCu+SolnBlk.Nghost; j++) {
      for (int i = SolnBlk.ICl-SolnBlk.Nghost; i <= SolnBlk.ICu+SolnBlk.Nghost; i++) {
	SolnBlk.sign[i][j] = SolnBlk.U[i][j].psi/sqrt(sqr(SolnBlk.Uo[i][j].psi)+dx);
      }
    }
    break;
  case EIKONAL_SIGN_FUNCTION_DERIVATIVE :
    // Peng, Merriman, Osher, Zhao and Kang
    for (int j = SolnBlk.JCl-SolnBlk.Nghost+1; j <= SolnBlk.JCu+SolnBlk.Nghost-1; j++) {
      for (int i = SolnBlk.ICl-SolnBlk.Nghost+1; i <= SolnBlk.ICu+SolnBlk.Nghost-1; i++) {
	dU = sqrt(sqr((SolnBlk.U[i+1][j].psi-SolnBlk.U[i][j].psi)/dx) + 
		  sqr((SolnBlk.U[i][j+1].psi-SolnBlk.U[i][j].psi)/dx));
	SolnBlk.sign[i][j] = SolnBlk.U[i][j].psi/sqrt(sqr(SolnBlk.U[i][j].psi) + sqr(dU*dx));
      }
    }
    break;
  };

  // Sign function calculated successfully.
  return 0;

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
void Linear_Reconstruction_GreenGauss(LevelSet2D_Quad_Block &SolnBlk,
				      const int i,
                                      const int j,
				      const int Limiter) {

  int error_flag;
  int n_pts, i_index[8], j_index[8];
  double U0Min, U0Max, UQuad[4], phi;
  double DxDx_ave, DxDy_ave, DyDy_ave;
  double l_north, l_south, l_east, l_west;
  Vector2D n_north, n_south, n_east, n_west, dX;
  LevelSet2DState U_nw, U_ne, U_sw, U_se, U_face, 
                  DU, DUDx_ave, DUDy_ave;

  // Carry out the limited solution reconstruction in each cell of the 
  // computational mesh.

  // Determine the number of neighbouring cells to be used in the 
  // reconstruction procedure.  Away from boundaries this 8 neighbours 
  // will be used.
  if (i == SolnBlk.ICl-SolnBlk.Nghost || i == SolnBlk.ICu+SolnBlk.Nghost ||
      j == SolnBlk.JCl-SolnBlk.Nghost || j == SolnBlk.JCu+SolnBlk.Nghost) {
    n_pts = 0;
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

      error_flag = Green_Gauss_Integration(SolnBlk.Grid.nodeNW(i,j).X,SolnBlk.UnNW(i, j),
					   SolnBlk.Grid.nodeSW(i,j).X,SolnBlk.UnSW(i, j),
					   SolnBlk.Grid.nodeSE(i,j).X,SolnBlk.UnSE(i, j),
					   SolnBlk.Grid.nodeNE(i,j).X,SolnBlk.UnNE(i, j),
					   SolnBlk.dUdx[i][j],SolnBlk.dUdy[i][j]);
      //return error_flag;

    // If < 8 neighbours are used, apply least-squares reconstruction
    } else {
      DUDx_ave = LevelSet2D_ZERO;
      DUDy_ave = LevelSet2D_ZERO;
      DxDx_ave = ZERO;
      DxDy_ave = ZERO;
      DyDy_ave = ZERO;

      for (int n = 0; n < n_pts; n++) {
	dX = SolnBlk.Grid.Cell[i_index[n]][j_index[n]].Xc - SolnBlk.Grid.Cell[i][j].Xc;
	DU = SolnBlk.U[i_index[n]][j_index[n]] - SolnBlk.U[i][j];
	DUDx_ave += DU*dX.x;
	DUDy_ave += DU*dX.y;
	DxDx_ave += dX.x*dX.x;
	DxDy_ave += dX.x*dX.y;
	DyDy_ave += dX.y*dX.y;
      }

      DUDx_ave = DUDx_ave/double(n_pts);
      DUDy_ave = DUDy_ave/double(n_pts);
      DxDx_ave = DxDx_ave/double(n_pts);
      DxDy_ave = DxDy_ave/double(n_pts);
      DyDy_ave = DyDy_ave/double(n_pts);
      SolnBlk.dUdx[i][j] = (DUDx_ave*DyDy_ave-DUDy_ave*DxDy_ave)/
	                   (DxDx_ave*DyDy_ave-DxDy_ave*DxDy_ave);
      SolnBlk.dUdy[i][j] = (DUDy_ave*DxDx_ave-DUDx_ave*DxDy_ave)/
	                   (DxDx_ave*DyDy_ave-DxDy_ave*DxDy_ave);
    }

    // Calculate slope limiters.    
    for (int k = 1; k <= NUM_VAR_LEVELSET2D; k++) {
      U0Min = SolnBlk.U[i][j][k];
      U0Max = U0Min;
      for (int n = 0; n < n_pts; n++) {
	U0Min = min(U0Min,SolnBlk.U[i_index[n]][j_index[n]][k]);
	U0Max = max(U0Max,SolnBlk.U[i_index[n]][j_index[n]][k]);
      }

      dX = SolnBlk.Grid.xfaceE(i,j) - SolnBlk.Grid.Cell[i][j].Xc;
      UQuad[0] = SolnBlk.U[i][j][k] + SolnBlk.dUdx[i][j][k]*dX.x +
                                      SolnBlk.dUdy[i][j][k]*dX.y;
      dX = SolnBlk.Grid.xfaceW(i,j) - SolnBlk.Grid.Cell[i][j].Xc;
      UQuad[1] = SolnBlk.U[i][j][k] + SolnBlk.dUdx[i][j][k]*dX.x +
                                      SolnBlk.dUdy[i][j][k]*dX.y;
      dX = SolnBlk.Grid.xfaceN(i,j) - SolnBlk.Grid.Cell[i][j].Xc;
      UQuad[2] = SolnBlk.U[i][j][k] + SolnBlk.dUdx[i][j][k]*dX.x +
                                      SolnBlk.dUdy[i][j][k]*dX.y;
      dX = SolnBlk.Grid.xfaceS(i,j) - SolnBlk.Grid.Cell[i][j].Xc;
      UQuad[3] = SolnBlk.U[i][j][k] + SolnBlk.dUdx[i][j][k]*dX.x +
                                      SolnBlk.dUdy[i][j][k]*dX.y;

      switch(Limiter) {
      case LIMITER_ONE :
	phi = ONE;
	break;
      case LIMITER_ZERO :
	phi = ZERO;
	break;
      case LIMITER_BARTH_JESPERSEN :
	phi = Limiter_BarthJespersen(UQuad,SolnBlk.U[i][j][k],U0Min,U0Max,4);
	break;
      case LIMITER_VENKATAKRISHNAN :
	phi = Limiter_Venkatakrishnan(UQuad,SolnBlk.U[i][j][k],U0Min,U0Max,4);
	break;
      case LIMITER_VANLEER :
	phi = Limiter_VanLeer(UQuad,SolnBlk.U[i][j][k],U0Min,U0Max,4);
	break;
      case LIMITER_VANALBADA :
	phi = Limiter_VanAlbada(UQuad,SolnBlk.U[i][j][k],U0Min,U0Max,4);
	break;
      default:
	phi = Limiter_BarthJespersen(UQuad,SolnBlk.U[i][j][k],U0Min,U0Max,4);
	break;
      };

      SolnBlk.phi[i][j][k] = phi;

    }

  } else {
    SolnBlk.dUdx[i][j] = LevelSet2D_ZERO;
    SolnBlk.dUdy[i][j] = LevelSet2D_ZERO;
    SolnBlk.phi[i][j]  = LevelSet2D_ZERO;
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
void Linear_Reconstruction_GreenGauss(LevelSet2D_Quad_Block &SolnBlk,
				      const int Limiter) {

  // Carry out the limited solution reconstruction in each cell of the
  // computational mesh.
  for (int j = SolnBlk.JCl-1; j <= SolnBlk.JCu+1; j++)
    for (int i = SolnBlk.ICl-1; i <= SolnBlk.ICu+1; i++)
      Linear_Reconstruction_GreenGauss(SolnBlk,i,j,Limiter);

}

/**********************************************************************
 * Routine: Laplacian_Reconstruction_GreenGauss                       *
 *                                                                    *
 * Peforms the reconstruction of the Laplacian of a specified         *
 * solution variable, n, within a given cell (i,j) of the             *
 * computational mesh for the specified Cartesian solution block.  A  *
 * Green-Gauss approach is used in the evaluation of the unlimited    *
 * solution gradients.                                                *
 *                                                                    *
 **********************************************************************/
void Laplacian_Reconstruction_GreenGauss(LevelSet2D_Quad_Block &SolnBlk,
					 const int i,
					 const int j,
					 const int n,
					 double &ddUdd) {

  int error_flag;
  Vector2D dUd_N, dUd_S, dUd_E, dUd_W, ddU;

  // Determine solution gradient on the north face of cell (i,j).
  error_flag = Green_Gauss_Integration(SolnBlk.Grid.nodeNW(i,j).X,SolnBlk.UnNW(i,j)[n],
				       SolnBlk.Grid.Cell[i][j].Xc,SolnBlk.U[i][j][n],
				       SolnBlk.Grid.nodeNE(i,j).X,SolnBlk.UnNE(i,j)[n],
				       SolnBlk.Grid.Cell[i][j+1].Xc,SolnBlk.U[i][j+1][n],
				       dUd_N.x,dUd_N.y);
  //if (error_flag) return error_flag;

  // Determine solution gradient on the south face of cell (i,j).
  error_flag = Green_Gauss_Integration(SolnBlk.Grid.nodeSW(i,j).X,SolnBlk.UnSW(i,j)[n],
				       SolnBlk.Grid.Cell[i][j-1].Xc,SolnBlk.U[i][j-1][n],
				       SolnBlk.Grid.nodeSE(i,j).X,SolnBlk.UnSE(i,j)[n],
				       SolnBlk.Grid.Cell[i][j].Xc,SolnBlk.U[i][j][n],
				       dUd_S.x,dUd_S.y);
  //if (error_flag) return error_flag;

  // Determine solution gradient on the east face of cell (i,j).
  error_flag = Green_Gauss_Integration(SolnBlk.Grid.nodeSE(i,j).X,SolnBlk.UnSE(i,j)[n],
 				       SolnBlk.Grid.Cell[i+1][j].Xc,SolnBlk.U[i+1][j][n],
 				       SolnBlk.Grid.nodeNE(i,j).X,SolnBlk.UnNE(i,j)[n],
 				       SolnBlk.Grid.Cell[i][j].Xc,SolnBlk.U[i][j][n],
 				       dUd_E.x,dUd_E.y);
  //if (error_flag) return error_flag;

  // Determine solution gradient on the west face of cell (i,j).
  error_flag = Green_Gauss_Integration(SolnBlk.Grid.nodeSW(i,j).X,SolnBlk.UnSW(i,j)[n],
 				       SolnBlk.Grid.Cell[i][j].Xc,SolnBlk.U[i][j][n],
 				       SolnBlk.Grid.nodeNW(i,j).X,SolnBlk.UnNW(i,j)[n],
 				       SolnBlk.Grid.Cell[i-1][j].Xc,SolnBlk.U[i-1][j][n],
 				       dUd_W.x,dUd_W.y);
  //if (error_flag) return error_flag;

  // Determine the cell-centred Laplacian.
  error_flag = Green_Gauss_Integration(SolnBlk.Grid.xfaceN(i,j),dUd_N.x,
				       SolnBlk.Grid.xfaceW(i,j),dUd_W.x,
				       SolnBlk.Grid.xfaceS(i,j),dUd_S.x,
				       SolnBlk.Grid.xfaceE(i,j),dUd_E.x,
				       ddU.x,ddU.y);
  //if (error_flag) return error_flag;
  ddUdd = ddU.x;
  error_flag = Green_Gauss_Integration(SolnBlk.Grid.xfaceN(i,j),dUd_N.y,
				       SolnBlk.Grid.xfaceW(i,j),dUd_W.y,
				       SolnBlk.Grid.xfaceS(i,j),dUd_S.y,
				       SolnBlk.Grid.xfaceE(i,j),dUd_E.y,
				       ddU.x,ddU.y);
  //if (error_flag) return error_flag;
  ddUdd += ddU.y;

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
void Linear_Reconstruction_LeastSquares(LevelSet2D_Quad_Block &SolnBlk,
				        const int i,
                                        const int j,
                                        const int Limiter) {

  int n_pts, i_index[8], j_index[8];
  double U0Min, U0Max, UQuad[4], phi;
  double DxDx_ave, DxDy_ave, DyDy_ave;
  Vector2D dX;
  LevelSet2DState DU, DUDx_ave, DUDy_ave;

  // Carry out the limited solution reconstruction in each cell of the 
  // computational mesh.
  if (i == SolnBlk.ICl-SolnBlk.Nghost || i == SolnBlk.ICu+SolnBlk.Nghost ||
      j == SolnBlk.JCl-SolnBlk.Nghost || j == SolnBlk.JCu+SolnBlk.Nghost) {
    n_pts = 0;
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
    DUDx_ave = LevelSet2D_ZERO;
    DUDy_ave = LevelSet2D_ZERO;
    DxDx_ave = ZERO;
    DxDy_ave = ZERO;
    DyDy_ave = ZERO;

    for (int n = 0; n < n_pts; n++) {
      dX = SolnBlk.Grid.Cell[i_index[n]][j_index[n]].Xc - SolnBlk.Grid.Cell[i][j].Xc;
      DU = SolnBlk.U[i_index[n]][j_index[n]] - SolnBlk.U[i][j];
      DUDx_ave += DU*dX.x;
      DUDy_ave += DU*dX.y;
      DxDx_ave += dX.x*dX.x;
      DxDy_ave += dX.x*dX.y;
      DyDy_ave += dX.y*dX.y;
    }

    DUDx_ave = DUDx_ave/double(n_pts);
    DUDy_ave = DUDy_ave/double(n_pts);
    DxDx_ave = DxDx_ave/double(n_pts);
    DxDy_ave = DxDy_ave/double(n_pts);
    DyDy_ave = DyDy_ave/double(n_pts);
    SolnBlk.dUdx[i][j] = (DUDx_ave*DyDy_ave-DUDy_ave*DxDy_ave)/
	                 (DxDx_ave*DyDy_ave-DxDy_ave*DxDy_ave);
    SolnBlk.dUdy[i][j] = (DUDy_ave*DxDx_ave-DUDx_ave*DxDy_ave)/
                         (DxDx_ave*DyDy_ave-DxDy_ave*DxDy_ave);

    // Calculate slope limiters.
    for (int n = 1; n <= NUM_VAR_LEVELSET2D; n++) {

      U0Min = SolnBlk.U[i][j][n];
      U0Max = U0Min;
      for (int n2 = 0; n2 < n_pts; n2++) {
	U0Min = min(U0Min,SolnBlk.U[i_index[n2]][j_index[n2]][n]);
	U0Max = max(U0Max,SolnBlk.U[i_index[n2]][j_index[n2]][n]);
      }

      dX = SolnBlk.Grid.xfaceE(i,j) - SolnBlk.Grid.Cell[i][j].Xc;
      UQuad[0] = SolnBlk.U[i][j][n] + SolnBlk.dUdx[i][j][n]*dX.x +
	                              SolnBlk.dUdy[i][j][n]*dX.y;
      dX = SolnBlk.Grid.xfaceW(i,j) - SolnBlk.Grid.Cell[i][j].Xc;
      UQuad[1] = SolnBlk.U[i][j][n] + SolnBlk.dUdx[i][j][n]*dX.x +
                                      SolnBlk.dUdy[i][j][n]*dX.y;
      dX = SolnBlk.Grid.xfaceN(i,j) - SolnBlk.Grid.Cell[i][j].Xc;
      UQuad[2] = SolnBlk.U[i][j][n] + SolnBlk.dUdx[i][j][n]*dX.x +
                                      SolnBlk.dUdy[i][j][n]*dX.y;
      dX = SolnBlk.Grid.xfaceS(i,j) - SolnBlk.Grid.Cell[i][j].Xc;
      UQuad[3] = SolnBlk.U[i][j][n] + SolnBlk.dUdx[i][j][n]*dX.x +
                                      SolnBlk.dUdy[i][j][n]*dX.y;

      switch(Limiter) {
      case LIMITER_ONE :
	phi = ONE;
	break;
      case LIMITER_ZERO :
	phi = ZERO;
	break;
      case LIMITER_BARTH_JESPERSEN :
	phi = Limiter_BarthJespersen(UQuad,SolnBlk.U[i][j][n],U0Min,U0Max,4);
	break;
      case LIMITER_VENKATAKRISHNAN :
	phi = Limiter_Venkatakrishnan(UQuad,SolnBlk.U[i][j][n],U0Min,U0Max,4);
	break;
      case LIMITER_VANLEER :
	phi = Limiter_VanLeer(UQuad,SolnBlk.U[i][j][n],U0Min,U0Max,4);
	break;
      case LIMITER_VANALBADA :
	phi = Limiter_VanAlbada(UQuad,SolnBlk.U[i][j][n],U0Min,U0Max,4);
	break;
      default:
	phi = Limiter_BarthJespersen(UQuad,SolnBlk.U[i][j][n],U0Min,U0Max,4);
	break;
      };

      SolnBlk.phi[i][j][n] = phi;

    }

  } else {
    SolnBlk.dUdx[i][j] = LevelSet2D_ZERO;
    SolnBlk.dUdy[i][j] = LevelSet2D_ZERO; 
    SolnBlk.phi[i][j]  = LevelSet2D_ZERO;
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
void Linear_Reconstruction_LeastSquares(LevelSet2D_Quad_Block &SolnBlk,
                                        const int Limiter) {

  // Carry out the limited solution reconstruction in each cell of the 
  // computational mesh.
  for (int j = SolnBlk.JCl-SolnBlk.Nghost; j <= SolnBlk.JCu+SolnBlk.Nghost; j++) {
    for (int i = SolnBlk.ICl-SolnBlk.Nghost; i <= SolnBlk.ICu+SolnBlk.Nghost; i++) {
      Linear_Reconstruction_LeastSquares(SolnBlk,i,j,Limiter);
    }
  }

}

/**********************************************************************
 * Routine: dUdx, dUdy, ddUdx2, ddUdy2, dddUdx3, dddUdy3              *
 *                                                                    *
 * Difference tables used for the essentially non-oscillatory         *
 * gradient reconstruction for the specified cell and the specified   *
 * solution variable.                                                 *
 *                                                                    *
 **********************************************************************/
double dUdx(LevelSet2D_Quad_Block &SolnBlk, const int &i, const int &j, const int &n) {
  return (SolnBlk.U[i+1][j][n]-SolnBlk.U[i][j][n])/(SolnBlk.Grid.Cell[i+1][j].Xc.x-SolnBlk.Grid.Cell[i][j].Xc.x);
}

double dUdy(LevelSet2D_Quad_Block &SolnBlk, const int &i, const int &j, const int &n) {
  return (SolnBlk.U[i][j+1][n]-SolnBlk.U[i][j][n])/(SolnBlk.Grid.Cell[i][j+1].Xc.y-SolnBlk.Grid.Cell[i][j].Xc.y);
}

double ddUdx2(LevelSet2D_Quad_Block &SolnBlk, const int &i, const int &j, const int &n) {
  return (dUdx(SolnBlk,i,j,n)-dUdx(SolnBlk,i-1,j,n))/(2.0*(SolnBlk.Grid.xfaceE(i,j).x-SolnBlk.Grid.xfaceW(i,j).x));
}

double ddUdy2(LevelSet2D_Quad_Block &SolnBlk, const int &i, const int &j, const int &n) {
  return (dUdy(SolnBlk,i,j,n)-dUdy(SolnBlk,i,j-1,n))/(2.0*(SolnBlk.Grid.xfaceN(i,j).y-SolnBlk.Grid.xfaceS(i,j).y));
}

double dddUdx3(LevelSet2D_Quad_Block &SolnBlk, const int &i, const int &j, const int &n) {
  return (ddUdx2(SolnBlk,i+1,j,n)-ddUdx2(SolnBlk,i,j,n))/(3.0*(SolnBlk.Grid.Cell[i+1][j].Xc.x-SolnBlk.Grid.Cell[i][j].Xc.x));
}

double dddUdy3(LevelSet2D_Quad_Block &SolnBlk, const int &i, const int &j, const int &n) {
  return (ddUdy2(SolnBlk,i,j+1,n)-ddUdy2(SolnBlk,i,j,n))/(3.0*(SolnBlk.Grid.Cell[i][j+1].Xc.y-SolnBlk.Grid.Cell[i][j].Xc.y));
}

/**********************************************************************
 * Routine: Reconstruction_EssentiallyNonOscillatory                  *
 *                                                                    *
 * Peforms the reconstruction of a limited piecewise linear solution  *
 * variable, n, within a given cell (i,j) of the computational mesh   *
 * for the specified Cartesian solution block.  A Essentially Non-    *
 * Oscillatory approach is used in the evaluation of the unlimited    *
 * solution gradients.                                                *
 *                                                                    *
 **********************************************************************/
void Reconstruction_EssentiallyNonOscillatory(LevelSet2D_Quad_Block &SolnBlk,
					      const int i,
					      const int j,
					      const int n,
					      const int i_Reconstruction) {

  int ism, isp, jsm, jsp;
  double ddUdx2a, ddUdx2b, ddUdy2a, ddUdy2b;
  double dddUdx3a, dddUdx3b, dddUdy3a, dddUdy3b;

  // Conduct linear ENO reconstruction:
  SolnBlk.dUdxm[i][j][n] = dUdx(SolnBlk,i-1,j,n);
  SolnBlk.dUdxp[i][j][n] = dUdx(SolnBlk,i,j,n);
  SolnBlk.dUdym[i][j][n] = dUdy(SolnBlk,i,j-1,n);
  SolnBlk.dUdyp[i][j][n] = dUdy(SolnBlk,i,j,n);

  // Conduct quadratic ENO reconstruction if required:
  if (i_Reconstruction == RECONSTRUCTION_QUADRATIC_ESSENTIALLY_NON_OSCILLATORY ||
      i_Reconstruction == RECONSTRUCTION_CUBIC_ESSENTIALLY_NON_OSCILLATORY) {
    // For dUdxm:
    ddUdx2a = ddUdx2(SolnBlk,i-1,j,n);
    ddUdx2b = ddUdx2(SolnBlk,i,j,n);
    if (fabs(ddUdx2a) <= fabs(ddUdx2b)) { ism = i-2; }
    else { ddUdx2a = ddUdx2b; ism = i-1; }
    SolnBlk.dUdxm[i][j][n] += ddUdx2a*(SolnBlk.Grid.Cell[i][j].Xc.x-SolnBlk.Grid.Cell[i-1][j].Xc.x);
    // For dUdxp:
    ddUdx2a = ddUdx2(SolnBlk,i,j,n);
    ddUdx2b = ddUdx2(SolnBlk,i+1,j,n);
    if (fabs(ddUdx2a) <= fabs(ddUdx2b)) { isp = i-1; }
    else { ddUdx2a = ddUdx2b; isp = i; }
    SolnBlk.dUdxp[i][j][n] -= ddUdx2a*(SolnBlk.Grid.Cell[i+1][j].Xc.x-SolnBlk.Grid.Cell[i][j].Xc.x);
    // For dUdym:
    ddUdy2a = ddUdy2(SolnBlk,i,j-1,n);
    ddUdy2b = ddUdy2(SolnBlk,i,j,n);
    if (fabs(ddUdy2a) <= fabs(ddUdy2b)) { jsm = j-2; }
    else { ddUdy2a = ddUdy2b; jsm = j-1; }
    SolnBlk.dUdym[i][j][n] += ddUdy2a*(SolnBlk.Grid.Cell[i][j].Xc.y-SolnBlk.Grid.Cell[i][j-1].Xc.y);
    // For dUdyp:
    ddUdy2a = ddUdy2(SolnBlk,i,j,n);
    ddUdy2b = ddUdy2(SolnBlk,i,j+1,n);
    if (fabs(ddUdy2a) <= fabs(ddUdy2b)) { jsp = j-1; }
    else { ddUdy2a = ddUdy2b; jsp = j; }
    SolnBlk.dUdyp[i][j][n] -= ddUdy2a*(SolnBlk.Grid.Cell[i][j+1].Xc.y-SolnBlk.Grid.Cell[i][j].Xc.y);
  }

  // Conduct cubic ENO reconstruction if required:
  if (i_Reconstruction == RECONSTRUCTION_CUBIC_ESSENTIALLY_NON_OSCILLATORY) {
    // For dUdxm:
    dddUdx3a = dddUdx3(SolnBlk,ism,j,n);
    dddUdx3b = dddUdx3(SolnBlk,ism+1,j,n);
    if (fabs(dddUdx3a) <= fabs(dddUdx3b)) { }
    else { dddUdx3a = dddUdx3b; }
    SolnBlk.dUdxm[i][j][n] += dddUdx3a*(3.0*sqr(double(i-ism)) - 6.0*double(i-ism) + 2.0)*sqr(SolnBlk.Grid.Cell[i-1+1][j].Xc.x-SolnBlk.Grid.Cell[i-1][j].Xc.x);
    // For dUdxp:
    dddUdx3a = dddUdx3(SolnBlk,isp,j,n);
    dddUdx3b = dddUdx3(SolnBlk,isp+1,j,n);
    if (fabs(dddUdx3a) <= fabs(dddUdx3b)) { }
    else { dddUdx3a = dddUdx3b; }
    SolnBlk.dUdxp[i][j][n] += dddUdx3a*(3.0*sqr(double(i-isp)) - 6.0*double(i-isp) + 2.0)*sqr(SolnBlk.Grid.Cell[i+1][j].Xc.x-SolnBlk.Grid.Cell[i][j].Xc.x);
    // For dUdym:
    dddUdy3a = dddUdy3(SolnBlk,i,jsm,n);
    dddUdy3b = dddUdy3(SolnBlk,i,jsm+1,n);
    if (fabs(dddUdy3a) <= fabs(dddUdy3b)) { }
    else { dddUdy3a = dddUdy3b; }
    SolnBlk.dUdym[i][j][n] += dddUdy3a*(3.0*sqr(double(j-jsm)) - 6.0*double(j-jsm) + 2.0)*sqr(SolnBlk.Grid.Cell[i][j-1+1].Xc.y-SolnBlk.Grid.Cell[i][j-1].Xc.y);
    // For dUdyp:
    dddUdy3a = dddUdy3(SolnBlk,i,jsp,n);
    dddUdy3b = dddUdy3(SolnBlk,j,jsp+1,n);
    if (fabs(dddUdy3a) <= fabs(dddUdy3b)) { }
    else { dddUdy3a = dddUdy3b; }
    SolnBlk.dUdyp[i][j][n] += dddUdy3a*(3.0*sqr(double(j-jsp)) - 6.0*double(j-jsp) + 2.0)*sqr(SolnBlk.Grid.Cell[i][j+1].Xc.y-SolnBlk.Grid.Cell[i][j].Xc.y);
  }

}

/**********************************************************************
 * Routine: Reconstruction_EssentiallyNonOscillatory                  *
 *                                                                    *
 * Peforms the reconstruction of a limited piecewise linear solution  *
 * variable, n, within each cell of the computational mesh for the    *
 * specified Cartesian solution block.  A Essentially Non-Oscillatory *
 * approach is used in the evaluation of the unlimited solution       *
 * gradients.                                                         *
 *                                                                    *
 **********************************************************************/
void Reconstruction_EssentiallyNonOscillatory(LevelSet2D_Quad_Block &SolnBlk,
					      const int n,
					      const int i_Reconstruction) {

  // Carry out the limited solution reconstruction in each cell of the 
  // computational mesh.
  for (int j = SolnBlk.JCl; j <= SolnBlk.JCu; j++) {
    for (int i = SolnBlk.ICl; i <= SolnBlk.ICu; i++) {
      Reconstruction_EssentiallyNonOscillatory(SolnBlk,i,j,n,i_Reconstruction);
    }
  }

}

/**********************************************************************
 * Routine: Reconstruction_WeightedEssentiallyNonOscillatory          *
 *                                                                    *
 * Calculates the 5th order accurate WENO gradient given 5 finite     *
 * difference stencils.                                               *
 *                                                                    *
 **********************************************************************/
double Reconstruction_WeightedEssentiallyNonOscillatory(const double v1,
							const double v2,
							const double v3,
							const double v4,
							const double v5) {
  
  // ENO approximations:
  double eno1, eno2, eno3;
  eno1 = v1/3 - 7*v2/6 + 11*v3/6;
  eno2 = -v2/6 + 5*v3/6 + v4/3;
  eno3 = v3/3 + 5*v4/6 - v5/6;

  // Smoothness Estimates:
  double S1, S2, S3;
  S1 = 13*sqr(v1-2*v2+v3)/12 + sqr(v1-4*v2+3*v3)/4;
  S2 = 13*sqr(v2-2*v3+v4)/12 + sqr(v2-v4)/4;
  S3 = 13*sqr(v3-2*v4+v5)/12 + sqr(3*v3-4*v4+v5)/4;

  // Alpha constants:
  double a1, a2, a3;
  a1 = 0.1/sqr(S1 + 1e-6);
  a2 = 0.6/sqr(S2 + 1e-6);
  a3 = 0.3/sqr(S3 + 1e-6);

  // Weightings:
  double weight1, weight2, weight3;
  weight1 = a1/(a1+a2+a3);
  weight2 = a2/(a1+a2+a3);
  weight3 = a3/(a1+a2+a3);

  // Return the WENO gradient as a weighted combination of the ENO
  // approximations.
  return weight1*eno1 + weight2*eno2 + weight3*eno3;

}

/**********************************************************************
 * Routine: Reconstruction_WeightedEssentiallyNonOscillatory          *
 *                                                                    *
 * Performs the reconstruction of a limited piecewise linear solution *
 * variable, n, within each cell of the computational mesh for the    *
 * specified Cartesian solution block.  A Weighted Essentially        *
 * Non-Oscillatory approach is used in the evaluation of the          *
 * unlimited solution gradients.                                      *
 *                                                                    *
 * Refer to Osher/Fedkiw's book, page 34, for the WENO scheme.        *
 *                                                                    *
 **********************************************************************/
void Reconstruction_WeightedEssentiallyNonOscillatory(LevelSet2D_Quad_Block &SolnBlk,
						      const int i,
						      const int j,
						      const int n) {

  double dxm, dxp, dym, dyp;
  double dxmm, dxpp, dymm, dypp;
  double dxmmm, dxppp, dymmm, dyppp;

  // Determine the distance to the neighbouring cells.
  dxm   = SolnBlk.Grid.Cell[i  ][j  ].Xc.x - SolnBlk.Grid.Cell[i-1][j  ].Xc.x;
  dxp   = SolnBlk.Grid.Cell[i+1][j  ].Xc.x - SolnBlk.Grid.Cell[i  ][j  ].Xc.x;
  dym   = SolnBlk.Grid.Cell[i  ][j  ].Xc.y - SolnBlk.Grid.Cell[i  ][j-1].Xc.y;
  dyp   = SolnBlk.Grid.Cell[i  ][j+1].Xc.y - SolnBlk.Grid.Cell[i  ][j  ].Xc.y;

  dxmm  = SolnBlk.Grid.Cell[i-1][j  ].Xc.x - SolnBlk.Grid.Cell[i-2][j  ].Xc.x;
  dxpp  = SolnBlk.Grid.Cell[i+2][j  ].Xc.x - SolnBlk.Grid.Cell[i+1][j  ].Xc.x;
  dymm  = SolnBlk.Grid.Cell[i  ][j-1].Xc.y - SolnBlk.Grid.Cell[i  ][j-2].Xc.y;
  dypp  = SolnBlk.Grid.Cell[i  ][j+2].Xc.y - SolnBlk.Grid.Cell[i  ][j+1].Xc.y;

  dxmmm = SolnBlk.Grid.Cell[i-2][j  ].Xc.x - SolnBlk.Grid.Cell[i-3][j  ].Xc.x;
  dxppp = SolnBlk.Grid.Cell[i+3][j  ].Xc.x - SolnBlk.Grid.Cell[i+2][j  ].Xc.x;
  dymmm = SolnBlk.Grid.Cell[i  ][j-2].Xc.y - SolnBlk.Grid.Cell[i  ][j-3].Xc.y;
  dyppp = SolnBlk.Grid.Cell[i  ][j+3].Xc.y - SolnBlk.Grid.Cell[i  ][j+2].Xc.y;

  // Calculate the gradient in the x-direction using a left-biased stencil
  double v1im, v2im, v3im, v4im, v5im;
  v1im = (SolnBlk.U[i-2][j  ][n] - SolnBlk.U[i-3][j  ][n])/(SolnBlk.Grid.Cell[i-2][j  ].Xc.x - SolnBlk.Grid.Cell[i-3][j  ].Xc.x);
  v2im = (SolnBlk.U[i-1][j  ][n] - SolnBlk.U[i-2][j  ][n])/(SolnBlk.Grid.Cell[i-1][j  ].Xc.x - SolnBlk.Grid.Cell[i-2][j  ].Xc.x);
  v3im = (SolnBlk.U[i  ][j  ][n] - SolnBlk.U[i-1][j  ][n])/(SolnBlk.Grid.Cell[i  ][j  ].Xc.x - SolnBlk.Grid.Cell[i-1][j  ].Xc.x);
  v4im = (SolnBlk.U[i+1][j  ][n] - SolnBlk.U[i  ][j  ][n])/(SolnBlk.Grid.Cell[i+1][j  ].Xc.x - SolnBlk.Grid.Cell[i  ][j  ].Xc.x);
  v5im = (SolnBlk.U[i+2][j  ][n] - SolnBlk.U[i+1][j  ][n])/(SolnBlk.Grid.Cell[i+2][j  ].Xc.x - SolnBlk.Grid.Cell[i+1][j  ].Xc.x);
  SolnBlk.dUdxm[i][j][n] = Reconstruction_WeightedEssentiallyNonOscillatory(v1im,v2im,v3im,v4im,v5im);

  // Calculate the gradient in the x-direction using a right-biased stencil.
  double v1ip, v2ip, v3ip, v4ip, v5ip;
  v1ip = (SolnBlk.U[i+3][j  ][n] - SolnBlk.U[i+2][j  ][n])/(SolnBlk.Grid.Cell[i+3][j  ].Xc.x - SolnBlk.Grid.Cell[i+2][j  ].Xc.x);
  v2ip = (SolnBlk.U[i+2][j  ][n] - SolnBlk.U[i+1][j  ][n])/(SolnBlk.Grid.Cell[i+2][j  ].Xc.x - SolnBlk.Grid.Cell[i+1][j  ].Xc.x);
  v3ip = (SolnBlk.U[i+1][j  ][n] - SolnBlk.U[i  ][j  ][n])/(SolnBlk.Grid.Cell[i+1][j  ].Xc.x - SolnBlk.Grid.Cell[i  ][j  ].Xc.x);
  v4ip = (SolnBlk.U[i  ][j  ][n] - SolnBlk.U[i-1][j  ][n])/(SolnBlk.Grid.Cell[i  ][j  ].Xc.x - SolnBlk.Grid.Cell[i-1][j  ].Xc.x);
  v5ip = (SolnBlk.U[i-1][j  ][n] - SolnBlk.U[i-2][j  ][n])/(SolnBlk.Grid.Cell[i-1][j  ].Xc.x - SolnBlk.Grid.Cell[i-2][j  ].Xc.x);
  SolnBlk.dUdxp[i][j][n] = Reconstruction_WeightedEssentiallyNonOscillatory(v1ip,v2ip,v3ip,v4ip,v5ip);

  // Calculate gradient in the y-direction using a left-biased stencil.
  double v1jm, v2jm, v3jm, v4jm, v5jm;
  v1jm = (SolnBlk.U[i  ][j-2][n] - SolnBlk.U[i  ][j-3][n])/(SolnBlk.Grid.Cell[i  ][j-2].Xc.y - SolnBlk.Grid.Cell[i  ][j-3].Xc.y);
  v2jm = (SolnBlk.U[i  ][j-1][n] - SolnBlk.U[i  ][j-2][n])/(SolnBlk.Grid.Cell[i  ][j-1].Xc.y - SolnBlk.Grid.Cell[i  ][j-2].Xc.y);
  v3jm = (SolnBlk.U[i  ][j  ][n] - SolnBlk.U[i  ][j-1][n])/(SolnBlk.Grid.Cell[i  ][j  ].Xc.y - SolnBlk.Grid.Cell[i  ][j-1].Xc.y);
  v4jm = (SolnBlk.U[i  ][j+1][n] - SolnBlk.U[i  ][j  ][n])/(SolnBlk.Grid.Cell[i  ][j+1].Xc.y - SolnBlk.Grid.Cell[i  ][j  ].Xc.y);
  v5jm = (SolnBlk.U[i  ][j+2][n] - SolnBlk.U[i  ][j+1][n])/(SolnBlk.Grid.Cell[i  ][j+2].Xc.y - SolnBlk.Grid.Cell[i  ][j+1].Xc.y);
  SolnBlk.dUdym[i][j][n] = Reconstruction_WeightedEssentiallyNonOscillatory(v1jm,v2jm,v3jm,v4jm,v5jm);

  // Calculate gradient in the y-direction using a right-biased stencil.
  double v1jp, v2jp, v3jp, v4jp, v5jp;
  v1jp = (SolnBlk.U[i  ][j+3][n] - SolnBlk.U[i  ][j+2][n])/(SolnBlk.Grid.Cell[i  ][j+3].Xc.y - SolnBlk.Grid.Cell[i  ][j+2].Xc.y);
  v2jp = (SolnBlk.U[i  ][j+2][n] - SolnBlk.U[i  ][j+1][n])/(SolnBlk.Grid.Cell[i  ][j+2].Xc.y - SolnBlk.Grid.Cell[i  ][j+1].Xc.y);
  v3jp = (SolnBlk.U[i  ][j+1][n] - SolnBlk.U[i  ][j  ][n])/(SolnBlk.Grid.Cell[i  ][j+1].Xc.y - SolnBlk.Grid.Cell[i  ][j  ].Xc.y);
  v4jp = (SolnBlk.U[i  ][j  ][n] - SolnBlk.U[i  ][j-1][n])/(SolnBlk.Grid.Cell[i  ][j  ].Xc.y - SolnBlk.Grid.Cell[i  ][j-1].Xc.y);
  v5jp = (SolnBlk.U[i  ][j-1][n] - SolnBlk.U[i  ][j-2][n])/(SolnBlk.Grid.Cell[i  ][j-1].Xc.y - SolnBlk.Grid.Cell[i  ][j-2].Xc.y);
  SolnBlk.dUdyp[i][j][n] = Reconstruction_WeightedEssentiallyNonOscillatory(v1jp,v2jp,v3jp,v4jp,v5jp);

}

/**********************************************************************
 * Routine: Reconstruction_WeightedEssentiallyNonOscillatory          *
 *                                                                    *
 * Peforms the reconstruction of a limited piecewise linear solution  *
 * variable, n, within each cell of the computational mesh for the    *
 * specified Cartesian solution block.  A Weighted Essentially Non-   *
 * Oscillatory approach is used in the evaluation of the unlimited    *
 * solution gradients.                                                *
 *                                                                    *
 **********************************************************************/
void Reconstruction_WeightedEssentiallyNonOscillatory(LevelSet2D_Quad_Block &SolnBlk,
						      const int n) {

  // Carry out the limited solution reconstruction in each cell of the 
  // computational mesh.
  for (int j = SolnBlk.JCl; j <= SolnBlk.JCu; j++) {
    for (int i = SolnBlk.ICl; i <= SolnBlk.ICu; i++) {
      Reconstruction_WeightedEssentiallyNonOscillatory(SolnBlk,i,j,n);
    }
  }

}

/**********************************************************************
 * Routine: Reconstruction_Curvature                                  *
 *                                                                    *
 * Performs the reconstruction of the curvature of the specified      *
 * variable, n, within each cell of the computational mesh for the    *
 * specified Cartesian solution block.                                *
 *                                                                    *
 **********************************************************************/
void Reconstruction_Curvature(LevelSet2D_Quad_Block &SolnBlk,
			      LevelSet2D_Input_Parameters &IP,
			      const int n) {

  /* Carry out the reconstruction of the curvature in each cell of
     the computational mesh. */
  switch(IP.i_Curvature_Scheme) {
  case CURVATURE_SCHEME_LAPLACIAN:
    for (int j = SolnBlk.JCl; j <= SolnBlk.JCu; j++) {
      for (int i = SolnBlk.ICl; i <= SolnBlk.ICu; i++) {
	Reconstruction_Curvature_Laplacian(SolnBlk,i,j,n);
      }
    }
    break;
  case CURVATURE_SCHEME_REGULAR:
    for (int j = SolnBlk.JCl; j <= SolnBlk.JCu; j++) {
      for (int i = SolnBlk.ICl; i <= SolnBlk.ICu; i++) {
	Reconstruction_Curvature_Regular(SolnBlk,i,j,n);
      }
    }
    break;
  default:
    for (int j = SolnBlk.JCl; j <= SolnBlk.JCu; j++) {
      for (int i = SolnBlk.ICl; i <= SolnBlk.ICu; i++) {
	Reconstruction_Curvature_Laplacian(SolnBlk,i,j,n);
      }
    }
    break;
  };

}

/**********************************************************************
 * Routine: Reconstruction_Curvature_Laplacian                        *
 *                                                                    *
 * Performs the reconstruction of the Laplacian of variable n, for    *
 * use in curvature driven motion of the level set function. N must   *
 * close to a signed distance function for accurate approximation.    *
 *                                                                    *
 * See: Osher p.43. Eqn. (4.6).                                        *
 *                                                                    *
 **********************************************************************/
void Reconstruction_Curvature_Laplacian(LevelSet2D_Quad_Block &SolnBlk,
					const int i,
					const int j,
					const int n) {

  double dx, dy, laplacian;

  // Carry out the reconstruction of the Laplacian in each cell.
  dx = SolnBlk.Grid.Cell[i][j].Xc.x - SolnBlk.Grid.Cell[i-1][j].Xc.x;
  dy = SolnBlk.Grid.Cell[i][j].Xc.y - SolnBlk.Grid.Cell[i][j-1].Xc.y;
  laplacian = ZERO;
  
  SolnBlk.ddUdxx[i][j].psi = (SolnBlk.U[i+1][j].psi - TWO*SolnBlk.U[i][j].psi + SolnBlk.U[i-1][j].psi)/sqr(dx);
  SolnBlk.ddUdyy[i][j].psi = (SolnBlk.U[i][j+1].psi - TWO*SolnBlk.U[i][j].psi + SolnBlk.U[i][j-1].psi)/sqr(dy);
  laplacian = SolnBlk.ddUdxx[i][j].psi + SolnBlk.ddUdyy[i][j].psi;
  
  if ( laplacian > (ONE/max(dx,dy)) ) {
    SolnBlk.kappa[i][j].psi = ONE/max(dx,dy);
  } else if ( laplacian < (-ONE/max(dx,dy)) ) {
    SolnBlk.kappa[i][j].psi = -ONE/max(dx,dy);
  } else {
    SolnBlk.kappa[i][j].psi = laplacian;
  }

  // Set gradient magnitude.
  // Note: The magnitude of the gradient is obviously not 1,
  // but for the approximation, we set it to 1.
  SolnBlk.gradMag[i][j].psi = ONE;

}

/**********************************************************************
 * Routine: Reconstruction_Curvature_Regular                          *
 *                                                                    *
 * Performs the reconstruction of the curvature of the level set for  *
 * use in curvature driven motion. The curvature is calculated using  *
 * the regular expression composed of single and double partial       *
 * derivatives in the x and y directions.                             *
 *                                                                    *
 * See: Osher, pp. 12. Eqn (1.8) and accompanying notes.
 *                                                                    *
 **********************************************************************/
void Reconstruction_Curvature_Regular(LevelSet2D_Quad_Block &SolnBlk,
				      const int i,
				      const int j,
				      const int n) {
  
  double phix, phiy;            // centered, first derivatives
  double phixx, phiyy, phixy;   // centered, second derivatives
  double dx, dy;

  dx = fabs(SolnBlk.Grid.Cell[i][j].Xc.x - SolnBlk.Grid.Cell[i+1][j].Xc.x);
  dy = fabs(SolnBlk.Grid.Cell[i][j].Xc.y - SolnBlk.Grid.Cell[i][j+1].Xc.y);

  // Get the first and second deratives.
  phix = HALF*(SolnBlk.U[i+1][j].psi - SolnBlk.U[i-1][j].psi)/dx;
  phiy = HALF*(SolnBlk.U[i][j+1].psi - SolnBlk.U[i][j-1].psi)/dy;
  phixx = (SolnBlk.U[i+1][j].psi - TWO*SolnBlk.U[i][j].psi + SolnBlk.U[i-1][j].psi) / sqr(dx);
  phiyy = (SolnBlk.U[i][j+1].psi - TWO*SolnBlk.U[i][j].psi + SolnBlk.U[i][j-1].psi) / sqr(dy);
  phixy = ((SolnBlk.U[i+1][j+1].psi-SolnBlk.U[i-1][j+1].psi)-(SolnBlk.U[i+1][j-1].psi-SolnBlk.U[i-1][j-1].psi)) / (FOUR*dx*dy);

  // Get the gradient magnitude.
  SolnBlk.gradMag[i][j].psi = sqrt( sqr(phix)+sqr(phiy) );

  // Get the curvature.
  if (SolnBlk.gradMag[i][j].psi < 1e-3) {
    // Catch the divide by zero case.
    SolnBlk.kappa[i][j].psi = ZERO;
  } else {
    // Regular expression.
    SolnBlk.kappa[i][j].psi = (sqr(phix)*phiyy - TWO*phix*phiy*phixy + sqr(phiy)*phixx) / (SolnBlk.gradMag[i][j].psi*SolnBlk.gradMag[i][j].psi*SolnBlk.gradMag[i][j].psi);
    if (SolnBlk.kappa[i][j].psi > (ONE/dx)) {
      SolnBlk.kappa[i][j].psi = ONE/dx;
    } else if (SolnBlk.kappa[i][j].psi < (-ONE/dx)) {
      SolnBlk.kappa[i][j].psi = -ONE/dx;
    }
  }

}

/**********************************************************************
 * Routine: Calculate_Refinement_Criteria                             *
 *                                                                    *
 * Calculate refinement criteria for the solution block.  Called by   *
 * AMR and Initial_AMR functions found in AMR.h.                      *
 *                                                                    *
 **********************************************************************/
void Calculate_Refinement_Criteria(double *refinement_criteria,
				   LevelSet2D_Input_Parameters &IP,
                                   int &number_refinement_criteria,
                                   LevelSet2D_Quad_Block &SolnBlk) {

  double curvature_psi_criteria, curvature_psi_criteria_max, ddpsidd,
         zero_level_set_criteria;
  Vector2D dpsid;
  int refinement_criteria_number;

  // Set the number of refinement criteria to be used (1):
  // (1) Refinement criteria based on the curvature of the level set
  // function (signed distance function).
  // (2) Refinement criteria based on the zero level set.
  number_refinement_criteria = IP.Number_of_Refinement_Criteria;

  // Initialize the refinement criteria for the block.
  curvature_psi_criteria_max = ZERO;
  zero_level_set_criteria = ZERO;

  // Calculate the refinement criteria for each cell of the computational
  // mesh and assign the maximum value for all cells as the refinement 
  // criteria for the solution block.
  for (int j = SolnBlk.JCl-SolnBlk.Nghost+1; j <= SolnBlk.JCu+SolnBlk.Nghost-1; j++) {
    for (int i = SolnBlk.ICl-SolnBlk.Nghost+1; i <= SolnBlk.ICu+SolnBlk.Nghost-1; i++) {
      // Evaluate refinement criteria based on the curvature of the 
      // level set function.
      curvature_psi_criteria = ZERO;
      if (SolnBlk.Un(i  ,j  ).psi*SolnBlk.Un(i+1,j  ).psi < ZERO ||
 	  SolnBlk.Un(i+1,j  ).psi*SolnBlk.Un(i+1,j+1).psi < ZERO ||
 	  SolnBlk.Un(i+1,j+1).psi*SolnBlk.Un(i  ,j+1).psi < ZERO ||
 	  SolnBlk.Un(i  ,j+1).psi*SolnBlk.Un(i  ,j  ).psi < ZERO) {
	// Determine the local gradient of the level set function.
 	Linear_Reconstruction_LeastSquares(SolnBlk,i,j,LIMITER_ONE);
 	dpsid = Vector2D(SolnBlk.dUdx[i][j].psi,SolnBlk.dUdy[i][j].psi);
	// Deteremine the local Laplacian of the signed distance function.
	Laplacian_Reconstruction_GreenGauss(SolnBlk,i,j,1,ddpsidd);
	// Determine the local curvature.
 	curvature_psi_criteria = ddpsidd/dpsid.abs();
	// Evaluate zero level set refinement criteria if required.
	if (IP.Refinement_Criteria_Zero_Level_Set) zero_level_set_criteria = ONE;
      }
      curvature_psi_criteria_max = max(curvature_psi_criteria_max,
				       fabs(curvature_psi_criteria));
    }
  }

  // Return the refinement criteria.
  refinement_criteria_number = 0;
  if (IP.Refinement_Criteria_Curvature) {
    refinement_criteria[refinement_criteria_number] = curvature_psi_criteria_max;
    refinement_criteria_number++;
  }
  if (IP.Refinement_Criteria_Zero_Level_Set) {
    refinement_criteria[refinement_criteria_number] = zero_level_set_criteria;
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
void Fix_Refined_Block_Boundaries(LevelSet2D_Quad_Block &SolnBlk,
                                  const int Fix_North_Boundary,
                                  const int Fix_South_Boundary,
                                  const int Fix_East_Boundary,
                                  const int Fix_West_Boundary) {

  // Adjust the node locations at the north boundary.
  if (Fix_North_Boundary) {
    for (int i = SolnBlk.Grid.INl+1; i <= SolnBlk.Grid.INu-1 ; i+=2) {
      SolnBlk.Grid.Node[i][SolnBlk.Grid.JNu].X = HALF*(SolnBlk.Grid.Node[i+1][SolnBlk.Grid.JNu].X+
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
      SolnBlk.Grid.Node[i][SolnBlk.Grid.JNl].X = 
	HALF*(SolnBlk.Grid.Node[i+1][SolnBlk.Grid.JNl].X+
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
      SolnBlk.Grid.Node[SolnBlk.Grid.INu][j].X = 
	HALF*(SolnBlk.Grid.Node[SolnBlk.Grid.INu][j+1].X+
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
      SolnBlk.Grid.Node[SolnBlk.Grid.INl][j].X = HALF*(SolnBlk.Grid.Node[SolnBlk.Grid.INl][j+1].X+
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
void Unfix_Refined_Block_Boundaries(LevelSet2D_Quad_Block &SolnBlk) {

  double sp_l, sp_r, sp_m, ds_ratio, dl, dr;

  // Return the nodes at the north boundary to their original positions.
  if (SolnBlk.Grid.BndNorthSpline.np != 0) {
    for (int i = SolnBlk.Grid.INl+1; i < SolnBlk.Grid.INu; i += 2) {
      sp_l = getS(SolnBlk.Grid.Node[i-1][SolnBlk.Grid.JNu].X,
		  SolnBlk.Grid.BndNorthSpline);
      sp_r = getS(SolnBlk.Grid.Node[i+1][SolnBlk.Grid.JNu].X,
		  SolnBlk.Grid.BndNorthSpline);
      dl = abs(SolnBlk.Grid.Node[i  ][SolnBlk.Grid.JNu].X - SolnBlk.Grid.Node[i-1][SolnBlk.Grid.JNu].X);
      dr = abs(SolnBlk.Grid.Node[i+1][SolnBlk.Grid.JNu].X - SolnBlk.Grid.Node[i  ][SolnBlk.Grid.JNu].X);
      ds_ratio = dl/(dl+dr);
      sp_m = sp_l + ds_ratio*(sp_r-sp_l);
      SolnBlk.Grid.Node[i][SolnBlk.Grid.JNu].X = Spline(sp_m,SolnBlk.Grid.BndNorthSpline);
    }
    for (int i = SolnBlk.ICl; i <= SolnBlk.ICu; i++)
      SolnBlk.U[i][SolnBlk.JCu] = (SolnBlk.Grid.Cell[i][SolnBlk.JCu].A/
				   SolnBlk.Grid.area(i,SolnBlk.JCu))*SolnBlk.U[i][SolnBlk.JCu];
  }

  // Return the nodes at the south boundary to their original positions.
  if (SolnBlk.Grid.BndSouthSpline.np != 0) {
    for (int i = SolnBlk.Grid.INl+1; i < SolnBlk.Grid.INu; i += 2) {
      sp_l = getS(SolnBlk.Grid.Node[i-1][SolnBlk.Grid.JNl].X,
		  SolnBlk.Grid.BndSouthSpline);
      sp_r = getS(SolnBlk.Grid.Node[i+1][SolnBlk.Grid.JNl].X,
		  SolnBlk.Grid.BndSouthSpline);
      dl = abs(SolnBlk.Grid.Node[i  ][SolnBlk.Grid.JNl].X - SolnBlk.Grid.Node[i-1][SolnBlk.Grid.JNl].X);
      dr = abs(SolnBlk.Grid.Node[i+1][SolnBlk.Grid.JNl].X - SolnBlk.Grid.Node[i  ][SolnBlk.Grid.JNl].X);
      ds_ratio = dl/(dl+dr);
      sp_m = sp_l + ds_ratio*(sp_r-sp_l);
      SolnBlk.Grid.Node[i][SolnBlk.Grid.JNl].X = Spline(sp_m,SolnBlk.Grid.BndSouthSpline);
    }
    for (int i = SolnBlk.ICl; i <= SolnBlk.ICu; i++)
      SolnBlk.U[i][SolnBlk.JCl] = (SolnBlk.Grid.Cell[i][SolnBlk.JCl].A/
				   SolnBlk.Grid.area(i,SolnBlk.JCl))*SolnBlk.U[i][SolnBlk.JCl];
  }

  // Return the nodes at the east boundary to their original positions.
  if (SolnBlk.Grid.BndEastSpline.np != 0) {
    for (int j = SolnBlk.Grid.JNl+1; j < SolnBlk.Grid.JNu; j += 2) {
      sp_l = getS(SolnBlk.Grid.Node[SolnBlk.Grid.INu][j-1].X,
		  SolnBlk.Grid.BndEastSpline);
      sp_r = getS(SolnBlk.Grid.Node[SolnBlk.Grid.INu][j+1].X,
		  SolnBlk.Grid.BndEastSpline);
      dl = abs(SolnBlk.Grid.Node[SolnBlk.Grid.INu][j  ].X - SolnBlk.Grid.Node[SolnBlk.Grid.INu][j-1].X);
      dr = abs(SolnBlk.Grid.Node[SolnBlk.Grid.INu][j+1].X - SolnBlk.Grid.Node[SolnBlk.Grid.INu][j  ].X);
      ds_ratio = dl/(dl+dr);
      sp_m = sp_l + ds_ratio*(sp_r-sp_l);
      SolnBlk.Grid.Node[SolnBlk.Grid.INu][j].X = Spline(sp_m,SolnBlk.Grid.BndEastSpline);
    }
    for (int j = SolnBlk.JCl; j <= SolnBlk.JCu; j++)
      SolnBlk.U[SolnBlk.ICu][j] = (SolnBlk.Grid.Cell[SolnBlk.ICu][j].A/
				   SolnBlk.Grid.area(SolnBlk.ICu,j))*SolnBlk.U[SolnBlk.ICu][j];
  }

  // Return the nodes at the west boundary to their original positions.
  if (SolnBlk.Grid.BndWestSpline.np != 0) {
    for (int j = SolnBlk.Grid.JNl+1; j < SolnBlk.Grid.JNu; j += 2) {
      sp_l = getS(SolnBlk.Grid.Node[SolnBlk.Grid.INl][j-1].X,
		  SolnBlk.Grid.BndWestSpline);
      sp_r = getS(SolnBlk.Grid.Node[SolnBlk.Grid.INl][j+1].X,
		  SolnBlk.Grid.BndWestSpline);
      dl = abs(SolnBlk.Grid.Node[SolnBlk.Grid.INl][j  ].X - SolnBlk.Grid.Node[SolnBlk.Grid.INl][j-1].X);
      dr = abs(SolnBlk.Grid.Node[SolnBlk.Grid.INl][j+1].X - SolnBlk.Grid.Node[SolnBlk.Grid.INl][j  ].X);
      ds_ratio = dl/(dl+dr);
      sp_m = sp_l + ds_ratio*(sp_r-sp_l);
      SolnBlk.Grid.Node[SolnBlk.Grid.INl][j].X = Spline(sp_m,SolnBlk.Grid.BndWestSpline);
    }
    for (int j = SolnBlk.JCl; j <= SolnBlk.JCu; j++)
      SolnBlk.U[SolnBlk.ICl][j] = (SolnBlk.Grid.Cell[SolnBlk.ICl][j].A/
				   SolnBlk.Grid.area(SolnBlk.ICl,j))*SolnBlk.U[SolnBlk.ICl][j];
  }

  // Reset the boundary condition types at the block boundaries.
  Set_BCs(SolnBlk.Grid);

  // Recompute the exterior nodes for the block quadrilateral mesh.
  Update_Exterior_Nodes(SolnBlk.Grid);

  // Recompute the cells for the block quadrilateral mesh.
  Update_Cells(SolnBlk.Grid);

}
