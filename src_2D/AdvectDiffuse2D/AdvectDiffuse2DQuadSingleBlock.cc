/*! \file AdvectDiffuse2DQuadSingleBlock.cc:
  \brief Temporary Single-Block Versions of Subroutines for 2D Advection Diffusion Equation
  Multi-Block Quadrilateral Mesh Solution Classes. */

/* Include required C++ libraries. */
// None

/* Using std namespace functions */
// None

/* Include CFFC header files */
#include "AdvectDiffuse2DQuad.h"    /* Include 2D advection diffusion equation quadrilateral mesh solution header file. */

/**************************************************************************
 * AdvectDiffuse2D_Quad_Block -- Single Block External Subroutines.       *
 **************************************************************************/

/******************************************************//**
 * Routine: Write_Solution_Block                        
 *                                                      
 * Writes the cell centred solution values of the       
 * specified quadrilateral solution block to the        
 * specified output stream for restart purposes.        
 *                                                      
 ********************************************************/
void Write_Solution_Block(AdvectDiffuse2D_Quad_Block &SolnBlk,
	                  ostream &Out_File) {

  Out_File << setprecision(14) << SolnBlk << setprecision(6);

}

/******************************************************//**
 * Routine: Read_Solution_Block                         
 *                                                      
 * Reads the cell centred solution values for the       
 * specified quadrilateral solution block from the      
 * specified input stream as required for restart       
 * purposes.                                            
 *                                                      
 ********************************************************/
void Read_Solution_Block(AdvectDiffuse2D_Quad_Block &SolnBlk,
	                 istream &In_File) {

  In_File >> SolnBlk;

}

/******************************************************//**
 * Routine: Broadcast_Solution_Block                    
 *                                                      
 * Broadcast quadrilateral solution block to all        
 * processors involved in the calculation from the      
 * primary processor using the MPI broadcast routine.   
 *                                                      
 ********************************************************/
void Broadcast_Solution_Block(AdvectDiffuse2D_Quad_Block &SolnBlk) {

#ifdef _MPI_VERSION
  int i, j, ni, nj, ng, block_allocated, buffer_size;
  double *buffer;

  // High-order related variables
  int NumberOfHighOrderVariables;
  vector<int> ReconstructionOrders;

  /* Broadcast the number of cells in each direction. */

  if (CFFC_Primary_MPI_Processor()) {
    ni = SolnBlk.NCi;
    nj = SolnBlk.NCj;
    ng = SolnBlk.Nghost;
    if (SolnBlk.U != NULL) {
      block_allocated = 1;
    } else {
      block_allocated = 0;
    } /* endif */ 

    // High-order variables and their reconstruction order
    NumberOfHighOrderVariables = SolnBlk.NumberOfHighOrderObjects();
    for (i = 0; i < NumberOfHighOrderVariables; ++i){
      ReconstructionOrders.push_back(SolnBlk.HighOrderVariable(i).RecOrder());
    }
  } /* endif */

  MPI::COMM_WORLD.Bcast(&ni, 1, MPI::INT, 0);
  MPI::COMM_WORLD.Bcast(&nj, 1, MPI::INT, 0);
  MPI::COMM_WORLD.Bcast(&ng, 1, MPI::INT, 0);
  MPI::COMM_WORLD.Bcast(&block_allocated, 1, MPI::INT, 0);

  // Broadcast the number of high-order variables and their reconstruction order
  MPI::COMM_WORLD.Bcast(&NumberOfHighOrderVariables, 1, MPI::INT, 0);
  if (!CFFC_Primary_MPI_Processor()) {
    // reserve memory for the reconstruction orders
    ReconstructionOrders.reserve(NumberOfHighOrderVariables);
  }
  for (i = 0; i < NumberOfHighOrderVariables; ++i){
    MPI::COMM_WORLD.Bcast(&ReconstructionOrders[i], 1, MPI::INT, 0);
  }

  /* On non-primary MPI processors, allocate (re-allocate) 
     memory for the quadrilateral solution block as necessary. */

  if (!CFFC_Primary_MPI_Processor()) {
    if (SolnBlk.NCi != ni || SolnBlk.NCj != nj || SolnBlk.Nghost != ng) { 
      if (SolnBlk.U != NULL) SolnBlk.deallocate();
      if (block_allocated) SolnBlk.allocate(ni-2*ng, nj-2*ng, ng); 
    } /* endif */
  } /* endif */

    /* Broadcast the axisymmetric/planar geometry indicator. */

  MPI::COMM_WORLD.Bcast(&(SolnBlk.Axisymmetric), 1, MPI::INT, 0);

  /* Broadcast the grid. */

  Broadcast_Quad_Block(SolnBlk.Grid);

  /* Allocate memory for high-order variables
     on non-primary MPI processors. */
  if (!CFFC_Primary_MPI_Processor()) {
    // allocate memory for high-order variables AFTER grid broadcast!
    SolnBlk.allocate_HighOrder(NumberOfHighOrderVariables,
			       ReconstructionOrders,
			       false); //< only the basics (e.g. no pseudo-inverse calculation)
    // allocate memory for high-order boundary conditions if necessary
    SolnBlk.allocate_HighOrder_BoundaryConditions();
  }


  /* Broadcast the solution state variables. */

  if (block_allocated) {
    ni = (SolnBlk.ICu+SolnBlk.Nghost) - (SolnBlk.ICl-SolnBlk.Nghost) + 1;
    nj = (SolnBlk.JCu+SolnBlk.Nghost) - (SolnBlk.JCl-SolnBlk.Nghost) + 1;
    buffer = new double[ni*nj];

    if (CFFC_Primary_MPI_Processor()) {
      buffer_size = 0;
      for (j  = SolnBlk.JCl-SolnBlk.Nghost ; j <= SolnBlk.JCu+SolnBlk.Nghost ; ++j ) {
	for ( i = SolnBlk.ICl-SolnBlk.Nghost ; i <= SolnBlk.ICu+SolnBlk.Nghost ; ++i ) {
	  buffer[buffer_size  ] = SolnBlk.U[i][j].u;
	  buffer_size += 1;
	} /* endfor */
      } /* endfor */
    } /* endif */

    buffer_size = ni*nj;
    MPI::COMM_WORLD.Bcast(buffer, buffer_size, MPI::DOUBLE, 0);

    if (!CFFC_Primary_MPI_Processor()) {
      buffer_size = 0;
      for (j  = SolnBlk.JCl-SolnBlk.Nghost ; j <= SolnBlk.JCu+SolnBlk.Nghost ; ++j ) {
	for ( i = SolnBlk.ICl-SolnBlk.Nghost ; i <= SolnBlk.ICu+SolnBlk.Nghost ; ++i ) {
	  SolnBlk.U[i][j].u    = buffer[buffer_size];
	  buffer_size += 1;
	} /* endfor */
      } /* endfor */
    } /* endif */

    delete []buffer; 
    buffer = NULL;

    ni = 1;
    nj = (SolnBlk.JCu+SolnBlk.Nghost) - (SolnBlk.JCl-SolnBlk.Nghost) + 1;
    buffer = new double[2*ni*nj];

    if (CFFC_Primary_MPI_Processor()) {
      buffer_size = 0;
      for (j  = SolnBlk.JCl-SolnBlk.Nghost ; j <= SolnBlk.JCu+SolnBlk.Nghost ; ++j ) {
	buffer[buffer_size   ] = SolnBlk.UoW[j].u;
	buffer[buffer_size+ 1] = SolnBlk.UoE[j].u;
	buffer_size += 2;
      } /* endfor */
    } /* endif */

    buffer_size = 2*ni*nj;
    MPI::COMM_WORLD.Bcast(buffer, buffer_size, MPI::DOUBLE, 0);

    if (!CFFC_Primary_MPI_Processor()) {
      buffer_size = 0;
      for (j  = SolnBlk.JCl-SolnBlk.Nghost ; j <= SolnBlk.JCu+SolnBlk.Nghost ; ++j ) {
	SolnBlk.UoW[j].u    = buffer[buffer_size   ];
	SolnBlk.UoE[j].u    = buffer[buffer_size+ 1];
	buffer_size += 2;
      } /* endfor */
    } /* endif */

    delete []buffer; 
    buffer = NULL;

    ni = (SolnBlk.ICu+SolnBlk.Nghost) - (SolnBlk.ICl-SolnBlk.Nghost) + 1;
    nj = 1;
    buffer = new double[2*ni*nj];

    if (CFFC_Primary_MPI_Processor()) {
      buffer_size = 0;
      for (i  = SolnBlk.ICl-SolnBlk.Nghost ; i <= SolnBlk.ICu+SolnBlk.Nghost ; ++i ) {
	buffer[buffer_size   ] = SolnBlk.UoS[i].u;
	buffer[buffer_size+ 1] = SolnBlk.UoN[i].u;
	buffer_size += 2;
      } /* endfor */
    } /* endif */

    buffer_size = 2*ni*nj;
    MPI::COMM_WORLD.Bcast(buffer, buffer_size, MPI::DOUBLE, 0);

    if (!CFFC_Primary_MPI_Processor()) {
      buffer_size = 0;
      for (i  = SolnBlk.ICl-SolnBlk.Nghost ; i <= SolnBlk.ICu+SolnBlk.Nghost ; ++i ) {
	SolnBlk.UoS[i].u    = buffer[buffer_size   ];
	SolnBlk.UoN[i].u    = buffer[buffer_size+ 1];
	buffer_size += 2;
      } /* endfor */
    } /* endif */

    delete []buffer; 
    buffer = NULL;

    /* Broadcast the high-order variables. */
    for (i = 0; i < NumberOfHighOrderVariables; ++i){
      SolnBlk.HighOrderVariable(i).Broadcast_HighOrder_Data(SolnBlk.Grid);
    }

  } /* endif */
#endif

}

#ifdef _MPI_VERSION
/******************************************************//**
 * Routine: Broadcast_Solution_Block                    
 *                                                      
 * Broadcast quadrilateral solution block to all        
 * processors associated with the specified communicator
 * from the specified processor using the MPI broadcast 
 * routine.                                             
 *                                                      
 ********************************************************/
void Broadcast_Solution_Block(AdvectDiffuse2D_Quad_Block &SolnBlk,
                              MPI::Intracomm &Communicator, 
                              const int Source_CPU) {

  int Source_Rank = 0;
  int i, j, ni, nj, ng, block_allocated, buffer_size;
  double *buffer;

  // High-order related variables
  int NumberOfHighOrderVariables;
  vector<int> ReconstructionOrders;

  /* Broadcast the number of cells in each direction. */

  if (CFFC_MPI::This_Processor_Number == Source_CPU) {
    ni = SolnBlk.NCi;
    nj = SolnBlk.NCj;
    ng = SolnBlk.Nghost;
    if (SolnBlk.U != NULL) {
      block_allocated = 1;
    } else {
      block_allocated = 0;
    } /* endif */ 

    // High-order variables and their reconstruction order
    NumberOfHighOrderVariables = SolnBlk.NumberOfHighOrderObjects();
    for (i = 0; i < NumberOfHighOrderVariables; ++i){
      ReconstructionOrders.push_back(SolnBlk.HighOrderVariable(i).RecOrder());
    }
  } /* endif */

  Communicator.Bcast(&ni, 1, MPI::INT, Source_Rank);
  Communicator.Bcast(&nj, 1, MPI::INT, Source_Rank);
  Communicator.Bcast(&ng, 1, MPI::INT, Source_Rank);
  Communicator.Bcast(&block_allocated, 1, MPI::INT, Source_Rank);

  // Broadcast the number of high-order variables and their reconstruction order
  Communicator.Bcast(&NumberOfHighOrderVariables, 1, MPI::INT, Source_Rank);
  if (!(CFFC_MPI::This_Processor_Number == Source_CPU)) {
    // reserve memory for the reconstruction orders
    ReconstructionOrders.reserve(NumberOfHighOrderVariables);
  }
  for (i = 0; i < NumberOfHighOrderVariables; ++i){
    Communicator.Bcast(&ReconstructionOrders[i], 1, MPI::INT, Source_Rank);
  }

  /* On non-source MPI processors, allocate (re-allocate) 
     memory for the quadrilateral solution block as necessary. */

  if (!(CFFC_MPI::This_Processor_Number == Source_CPU)) {
    if (SolnBlk.NCi != ni || SolnBlk.NCj != nj || SolnBlk.Nghost != ng) { 
      if (SolnBlk.U != NULL) SolnBlk.deallocate();
      if (block_allocated) SolnBlk.allocate(ni-2*ng, nj-2*ng, ng); 
    } /* endif */
  } /* endif */

    /* Broadcast the axisymmetric/planar flow indicator. */

  Communicator.Bcast(&(SolnBlk.Axisymmetric), 1, MPI::INT, Source_Rank);

  /* Broadcast the grid. */

  Broadcast_Quad_Block(SolnBlk.Grid, Communicator, Source_CPU);

  /* Allocate memory for high-order variables
     on non-source MPI processors. */
  if (!(CFFC_MPI::This_Processor_Number == Source_CPU)) {
    // allocate memory for high-order variables AFTER grid broadcast!
    SolnBlk.allocate_HighOrder(NumberOfHighOrderVariables,
			       ReconstructionOrders,
			       false); //< only the basics (e.g. no pseudo-inverse calculation)
    // allocate memory for high-order boundary conditions if necessary
    SolnBlk.allocate_HighOrder_BoundaryConditions();
  }

  /* Broadcast the solution state variables. */

  if (block_allocated) {
    ni = (SolnBlk.ICu+SolnBlk.Nghost) - (SolnBlk.ICl-SolnBlk.Nghost) + 1;
    nj = (SolnBlk.JCu+SolnBlk.Nghost) - (SolnBlk.JCl-SolnBlk.Nghost) + 1;
    buffer = new double[ni*nj];

    if (CFFC_MPI::This_Processor_Number == Source_CPU) {
      buffer_size = 0;
      for (j  = SolnBlk.JCl-SolnBlk.Nghost ; j <= SolnBlk.JCu+SolnBlk.Nghost ; ++j ) {
	for ( i = SolnBlk.ICl-SolnBlk.Nghost ; i <= SolnBlk.ICu+SolnBlk.Nghost ; ++i ) {
	  buffer[buffer_size  ] = SolnBlk.U[i][j].u;
	  buffer_size += 1;
	} /* endfor */
      } /* endfor */
    } /* endif */

    buffer_size = ni*nj;
    Communicator.Bcast(buffer, buffer_size, MPI::DOUBLE, Source_Rank);

    if (!(CFFC_MPI::This_Processor_Number == Source_CPU)) {
      buffer_size = 0;
      for (j  = SolnBlk.JCl-SolnBlk.Nghost ; j <= SolnBlk.JCu+SolnBlk.Nghost ; ++j ) {
	for ( i = SolnBlk.ICl-SolnBlk.Nghost ; i <= SolnBlk.ICu+SolnBlk.Nghost ; ++i ) {
	  SolnBlk.U[i][j].u    = buffer[buffer_size];
	  buffer_size += 1;
	} /* endfor */
      } /* endfor */
    } /* endif */

    delete []buffer; 
    buffer = NULL;

    ni = 1;
    nj = (SolnBlk.JCu+SolnBlk.Nghost) - (SolnBlk.JCl-SolnBlk.Nghost) + 1;
    buffer = new double[2*ni*nj];

    if (CFFC_MPI::This_Processor_Number == Source_CPU) {
      buffer_size = 0;
      for (j  = SolnBlk.JCl-SolnBlk.Nghost ; j <= SolnBlk.JCu+SolnBlk.Nghost ; ++j ) {
	buffer[buffer_size   ] = SolnBlk.UoW[j].u;
	buffer[buffer_size+ 1] = SolnBlk.UoE[j].u;
	buffer_size += 2;
      } /* endfor */
    } /* endif */

    buffer_size = 2*ni*nj;
    Communicator.Bcast(buffer, buffer_size, MPI::DOUBLE, Source_Rank);

    if (!(CFFC_MPI::This_Processor_Number == Source_CPU)) {
      buffer_size = 0;
      for (j  = SolnBlk.JCl-SolnBlk.Nghost ; j <= SolnBlk.JCu+SolnBlk.Nghost ; ++j ) {
	SolnBlk.UoW[j].u    = buffer[buffer_size   ];
	SolnBlk.UoE[j].u    = buffer[buffer_size+ 1];
	buffer_size += 2;
      } /* endfor */
    } /* endif */

    delete []buffer; 
    buffer = NULL;

    ni = (SolnBlk.ICu+SolnBlk.Nghost) - (SolnBlk.ICl-SolnBlk.Nghost) + 1;
    nj = 1;
    buffer = new double[2*ni*nj];

    if (CFFC_MPI::This_Processor_Number == Source_CPU) {
      buffer_size = 0;
      for (i  = SolnBlk.ICl-SolnBlk.Nghost ; i <= SolnBlk.ICu+SolnBlk.Nghost ; ++i ) {
	buffer[buffer_size   ] = SolnBlk.UoS[i].u;
	buffer[buffer_size+ 1] = SolnBlk.UoN[i].u;
	buffer_size += 2;
      } /* endfor */
    } /* endif */

    buffer_size = 2*ni*nj;
    Communicator.Bcast(buffer, buffer_size, MPI::DOUBLE, Source_Rank);

    if (!(CFFC_MPI::This_Processor_Number == Source_CPU)) {
      buffer_size = 0;
      for (i  = SolnBlk.ICl-SolnBlk.Nghost ; i <= SolnBlk.ICu+SolnBlk.Nghost ; ++i ) {
	SolnBlk.UoS[i].u    = buffer[buffer_size   ];
	SolnBlk.UoN[i].u    = buffer[buffer_size+ 1];
	buffer_size += 2;
      } /* endfor */
    } /* endif */

    delete []buffer; 
    buffer = NULL;

    /* Broadcast the high-order variables. */
    for (i = 0; i < NumberOfHighOrderVariables; ++i){
      SolnBlk.HighOrderVariable(i).Broadcast_HighOrder_Data(Communicator,
							    Source_CPU,
							    SolnBlk.Grid);
    }

  } /* endif */

}
#endif


/******************************************************//**
 * Routine: Copy_Solution_Block                         
 *                                                      
 * Copies the solution information of quadrilateral     
 * solution block SolnBlk2 to SolnBlk1.                 
 *                                                      
 ********************************************************/
void Copy_Solution_Block(AdvectDiffuse2D_Quad_Block &SolnBlk1,
                         const AdvectDiffuse2D_Quad_Block &SolnBlk2) {
  SolnBlk1.makeCopy(SolnBlk2);
}

/******************************************************//**
 * Routine: Prolong_Solution_Block                      
 *                                                      
 * Prolongs the solution information of one of the      
 * specified sectors of the original quadrilateral      
 * solution block SolnBlk_Original to the refined       
 * solution block SolnBlk_Fine.                         
 *                                                      
 ********************************************************/
int Prolong_Solution_Block(AdvectDiffuse2D_Quad_Block &SolnBlk_Fine,
			   AdvectDiffuse2D_Quad_Block &SolnBlk_Original,
			   const int Sector) {

  int i, j, k, i_min, i_max, j_min, j_max, mesh_refinement_permitted;
  double area_total_fine;

  // High-order related variables
  int SolnBlk_Original_NumberOfHighOrderVariables(SolnBlk_Original.NumberOfHighOrderObjects());
  vector<int> SolnBlk_Original_ReconstructionOrders; 

  /* Allocate (re-allocate) memory for the solution
     of the refined quadrilateral solution block as necessary. */

  if ( (SolnBlk_Original.NCi-2*SolnBlk_Original.Nghost-
	2*((SolnBlk_Original.NCi-2*SolnBlk_Original.Nghost)/2) != 0) || 
       (SolnBlk_Original.NCj-2*SolnBlk_Original.Nghost-
	2*((SolnBlk_Original.NCj-2*SolnBlk_Original.Nghost)/2) != 0) ||
       (SolnBlk_Original.NCi-2*SolnBlk_Original.Nghost < 2*SolnBlk_Original.Nghost) ||
       (SolnBlk_Original.NCj-2*SolnBlk_Original.Nghost < 2*SolnBlk_Original.Nghost) ||
       (SolnBlk_Original.Grid.Node == NULL) ||
       (SolnBlk_Original.U == NULL) ) {
    mesh_refinement_permitted = 0;
  } else {
    mesh_refinement_permitted = 1;
    if (SolnBlk_Fine.NCi != SolnBlk_Original.NCi || 
	SolnBlk_Fine.NCj != SolnBlk_Original.NCj ||
	SolnBlk_Fine.Nghost != SolnBlk_Original.Nghost) {
      if (SolnBlk_Fine.U != NULL) SolnBlk_Fine.deallocate();
      if (SolnBlk_Original.U != NULL) SolnBlk_Fine.allocate(SolnBlk_Original.NCi-2*SolnBlk_Original.Nghost,
							    SolnBlk_Original.NCj-2*SolnBlk_Original.Nghost,
							    SolnBlk_Original.Nghost);
      // In this case, create the refined mesh. /* 
      Refine_Mesh(SolnBlk_Fine.Grid, 
		  SolnBlk_Original.Grid,
		  Sector);
    } /* endif */

    // Allocate high-order objects if necessary
    for (i = 0; i < SolnBlk_Original_NumberOfHighOrderVariables; ++i){
      SolnBlk_Original_ReconstructionOrders.push_back(SolnBlk_Original.HighOrderVariable(i).RecOrder());
    }
    SolnBlk_Fine.allocate_HighOrder(SolnBlk_Original_NumberOfHighOrderVariables,
				    SolnBlk_Original_ReconstructionOrders);
    // allocate memory for high-order boundary conditions if necessary
    SolnBlk_Fine.allocate_HighOrder_BoundaryConditions();

  } /* endif */

  if (mesh_refinement_permitted) {
   
    /* Set the axisymmetric/planar flow indicator for the fine solution block. */

    SolnBlk_Fine.Axisymmetric = SolnBlk_Original.Axisymmetric;

    /* Prolong the solution information from original solution block
       to the refined solution block. */

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
    } /* endswitch */

    int iCell, jCell;		// indexes for the first fine cell

    for ( j  = j_min; j <= j_max; ++j ) {
      for ( i = i_min; i <= i_max; ++i ) {

	// Compute the indexes of the current fine cell (the first out of the 4 ones)
	iCell = 2*(i-i_min)+SolnBlk_Fine.ICl;
	jCell = 2*(j-j_min)+SolnBlk_Fine.JCl;

	// Use direct injection
	// Update fine cell I
	SolnBlk_Fine.U[iCell  ][jCell  ] = SolnBlk_Original.U[i][j];

	// Update fine cell II
	SolnBlk_Fine.U[iCell+1][jCell  ] = SolnBlk_Original.U[i][j];

	// Update fine cell III
	SolnBlk_Fine.U[iCell  ][jCell+1] = SolnBlk_Original.U[i][j];

	// Update fine cell IV
	SolnBlk_Fine.U[iCell+1][jCell+1] = SolnBlk_Original.U[i][j];
      } /* endfor */
    } /* endfor */

    for ( j  = j_min-SolnBlk_Original.Nghost/2; j <= j_max+SolnBlk_Original.Nghost/2 ; ++j ) {
      
      SolnBlk_Fine.UoW[2*(j-j_min)+SolnBlk_Fine.JCl  ] 	= SolnBlk_Original.UoW[j];
      SolnBlk_Fine.UoW[2*(j-j_min)+SolnBlk_Fine.JCl+1]	= SolnBlk_Original.UoW[j];

      SolnBlk_Fine.UoE[2*(j-j_min)+SolnBlk_Fine.JCl  ]	= SolnBlk_Original.UoE[j];
      SolnBlk_Fine.UoE[2*(j-j_min)+SolnBlk_Fine.JCl+1]	= SolnBlk_Original.UoE[j];
    } /* endfor */

    for ( i = i_min-SolnBlk_Original.Nghost/2 ; i <= i_max+SolnBlk_Original.Nghost/2 ; ++i ) {
      SolnBlk_Fine.UoS[2*(i-i_min)+SolnBlk_Fine.ICl  ]	= SolnBlk_Original.UoS[i];
      SolnBlk_Fine.UoS[2*(i-i_min)+SolnBlk_Fine.ICl+1]	= SolnBlk_Original.UoS[i];

      SolnBlk_Fine.UoN[2*(i-i_min)+SolnBlk_Fine.ICl  ]	= SolnBlk_Original.UoN[i];
      SolnBlk_Fine.UoN[2*(i-i_min)+SolnBlk_Fine.ICl+1]	= SolnBlk_Original.UoN[i];
    } /* endfor */

    // Set the reference values for the boundary states to the ones from the Original solution block
    SolnBlk_Fine.Set_Reference_Values_For_Boundary_States(SolnBlk_Original.Ref_State_BC_North,
							  SolnBlk_Original.Ref_State_BC_South,
							  SolnBlk_Original.Ref_State_BC_East,
							  SolnBlk_Original.Ref_State_BC_West);

    // Enforce analytic values for the boundary reference states if defined
    SolnBlk_Fine.Set_Boundary_Reference_States();
    
  } /* endif */
  
  
    // Prolongation of solution block was successful.
  return 0;
}

/******************************************************//**
 * Routine: Restrict_Solution_Block                     
 *                                                      
 * Restricts the solution information of four original  
 * fine quadrilateral solution blocks to the coarse     
 * solution block SolnBlk_Coarse.                       
 *                                                      
 ********************************************************/
int Restrict_Solution_Block(AdvectDiffuse2D_Quad_Block &SolnBlk_Coarse,
			    AdvectDiffuse2D_Quad_Block &SolnBlk_Original_SW,
			    AdvectDiffuse2D_Quad_Block &SolnBlk_Original_SE,
			    AdvectDiffuse2D_Quad_Block &SolnBlk_Original_NW,
			    AdvectDiffuse2D_Quad_Block &SolnBlk_Original_NE) {

  int i, j, i_coarse, j_coarse, mesh_coarsening_permitted;

  // High-order related variables
  int SolnBlk_Original_NumberOfHighOrderVariables(SolnBlk_Original_SW.NumberOfHighOrderObjects());
  vector<int> SolnBlk_Original_ReconstructionOrders; 
 
  /* Allocate memory for the cells and nodes for the 
     coarsened quadrilateral mesh block. */

  if ( (SolnBlk_Original_SW.NCi-2*SolnBlk_Original_SW.Nghost-
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
       (SolnBlk_Original_NE.U == NULL) ) {
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
    } /* endif */

    // Allocate high-order objects if necessary
    for (i = 0; i < SolnBlk_Original_NumberOfHighOrderVariables; ++i){
      SolnBlk_Original_ReconstructionOrders.push_back(SolnBlk_Original_SW.HighOrderVariable(i).RecOrder());
    }
    SolnBlk_Coarse.allocate_HighOrder(SolnBlk_Original_NumberOfHighOrderVariables,
				      SolnBlk_Original_ReconstructionOrders);
    // allocate memory for high-order boundary conditions if necessary
    SolnBlk_Coarse.allocate_HighOrder_BoundaryConditions();

  } /* endif */

  if (mesh_coarsening_permitted) {

    /* Set the axisymmetric/planar flow indicator for the coarse solution block. */

    SolnBlk_Coarse.Axisymmetric = SolnBlk_Original_SW.Axisymmetric;

    /* Restrict the solution information from the four original solution blocks
       to the newly coarsened solution block. */

    // South-West corner fine block:	
    for ( j = SolnBlk_Original_SW.JCl; j <= SolnBlk_Original_SW.JCu ; j += 2 ) {
      for ( i = SolnBlk_Original_SW.ICl ; i <= SolnBlk_Original_SW.ICu ; i += 2 ) {
	i_coarse = (i-SolnBlk_Original_SW.ICl)/2+
	  SolnBlk_Coarse.ICl;
	j_coarse = (j-SolnBlk_Original_SW.JCl)/2+
	  SolnBlk_Coarse.JCl;
	SolnBlk_Coarse.U[i_coarse][j_coarse] = ( (SolnBlk_Original_SW.Grid.Cell[i  ][j  ].A*
						  SolnBlk_Original_SW.U[i  ][j  ] +
						  SolnBlk_Original_SW.Grid.Cell[i+1][j  ].A*
						  SolnBlk_Original_SW.U[i+1][j  ] + 
						  SolnBlk_Original_SW.Grid.Cell[i  ][j+1].A*
						  SolnBlk_Original_SW.U[i  ][j+1] +
						  SolnBlk_Original_SW.Grid.Cell[i+1][j+1].A*
						  SolnBlk_Original_SW.U[i+1][j+1]) /
						 (SolnBlk_Original_SW.Grid.Cell[i  ][j  ].A +
						  SolnBlk_Original_SW.Grid.Cell[i+1][j  ].A +
						  SolnBlk_Original_SW.Grid.Cell[i  ][j+1].A +
						  SolnBlk_Original_SW.Grid.Cell[i+1][j+1].A) );
	//                                      SolnBlk_Coarse.Grid.Cell[i_coarse][j_coarse].A;
      } /* endfor */
    } /* endfor */

    for ( j = SolnBlk_Original_SW.JCl-SolnBlk_Original_SW.Nghost; 
	  j <= SolnBlk_Original_SW.JCu+SolnBlk_Original_SW.Nghost ; j += 2 ) {
      j_coarse = (j-SolnBlk_Original_SW.JCl)/2+SolnBlk_Coarse.JCl;
      SolnBlk_Coarse.UoW[j_coarse] = SolnBlk_Original_SW.UoW[j];
      if (j == SolnBlk_Original_SW.JCl-SolnBlk_Original_SW.Nghost) {
	SolnBlk_Coarse.UoW[j_coarse-1] = SolnBlk_Original_SW.UoW[j];
      } /* endif */
    } /* endfor */

    for ( i = SolnBlk_Original_SW.ICl-SolnBlk_Original_SW.Nghost ; 
	  i <= SolnBlk_Original_SW.ICu+SolnBlk_Original_SW.Nghost ; i += 2) {
      i_coarse = (i-SolnBlk_Original_SW.ICl)/2+SolnBlk_Coarse.ICl;
      SolnBlk_Coarse.UoS[i_coarse] = SolnBlk_Original_SW.UoS[i];
      if (i == SolnBlk_Original_SW.ICl-SolnBlk_Original_SW.Nghost) {
	SolnBlk_Coarse.UoS[i_coarse-1] = SolnBlk_Original_SW.UoS[i];
      } /* endif */
    } /* endfor */
    
      // South-East corner fine block:
    for ( j = SolnBlk_Original_SE.JCl; j <= SolnBlk_Original_SE.JCu ; j += 2 ) {
      for ( i = SolnBlk_Original_SE.ICl ; i <= SolnBlk_Original_SE.ICu ; i += 2 ) {
	i_coarse = (i-SolnBlk_Original_SE.ICl)/2+
	  (SolnBlk_Coarse.ICu-SolnBlk_Coarse.ICl+1)/2+SolnBlk_Coarse.ICl;
	j_coarse = (j-SolnBlk_Original_SE.JCl)/2+
	  SolnBlk_Coarse.JCl;
	SolnBlk_Coarse.U[i_coarse][j_coarse] = ( (SolnBlk_Original_SE.Grid.Cell[i  ][j  ].A*
						  SolnBlk_Original_SE.U[i  ][j  ] +
						  SolnBlk_Original_SE.Grid.Cell[i+1][j  ].A*
						  SolnBlk_Original_SE.U[i+1][j  ] + 
						  SolnBlk_Original_SE.Grid.Cell[i  ][j+1].A*
						  SolnBlk_Original_SE.U[i  ][j+1] +
						  SolnBlk_Original_SE.Grid.Cell[i+1][j+1].A*
						  SolnBlk_Original_SE.U[i+1][j+1]) /
						 (SolnBlk_Original_SE.Grid.Cell[i  ][j  ].A +
						  SolnBlk_Original_SE.Grid.Cell[i+1][j  ].A +
						  SolnBlk_Original_SE.Grid.Cell[i  ][j+1].A +
						  SolnBlk_Original_SE.Grid.Cell[i+1][j+1].A) );
	//                                        SolnBlk_Coarse.Grid.Cell[i_coarse][j_coarse].A;
      } /* endfor */
    } /* endfor */

    for ( j = SolnBlk_Original_SE.JCl-SolnBlk_Original_SE.Nghost; 
	  j <= SolnBlk_Original_SE.JCu+SolnBlk_Original_SE.Nghost ; j += 2 ) {
      j_coarse = (j-SolnBlk_Original_SE.JCl)/2+SolnBlk_Coarse.JCl;
      SolnBlk_Coarse.UoE[j_coarse] = SolnBlk_Original_SE.UoE[j];
      if (j == SolnBlk_Original_SE.JCl-SolnBlk_Original_SE.Nghost) {
	SolnBlk_Coarse.UoE[j_coarse-1] = SolnBlk_Original_SE.UoE[j];
      } /* endif */
    } /* endfor */

    for ( i = SolnBlk_Original_SE.ICl-SolnBlk_Original_SE.Nghost ; 
	  i <= SolnBlk_Original_SE.ICu+SolnBlk_Original_SE.Nghost ; i += 2) {
      i_coarse = (i-SolnBlk_Original_SE.ICl)/2+
	(SolnBlk_Coarse.ICu-SolnBlk_Coarse.ICl+1)/2+SolnBlk_Coarse.ICl;
      SolnBlk_Coarse.UoS[i_coarse] = SolnBlk_Original_SE.UoS[i];
      if (i == SolnBlk_Original_SE.ICu+SolnBlk_Original_SE.Nghost) {
	SolnBlk_Coarse.UoS[i_coarse+1] = SolnBlk_Original_SE.UoS[i];
      } /* endif */
    } /* endfor */

      // North-West corner fine block:
    for ( j = SolnBlk_Original_NW.JCl; j <= SolnBlk_Original_NW.JCu ; j += 2 ) {
      for ( i = SolnBlk_Original_NW.ICl ; i <= SolnBlk_Original_NW.ICu ; i += 2 ) {
	i_coarse = (i-SolnBlk_Original_NW.ICl)/2+
	  SolnBlk_Coarse.ICl;
	j_coarse = (j-SolnBlk_Original_NW.JCl)/2+
	  (SolnBlk_Coarse.JCu-SolnBlk_Coarse.JCl+1)/2+SolnBlk_Coarse.JCl;
	SolnBlk_Coarse.U[i_coarse][j_coarse] = ( (SolnBlk_Original_NW.Grid.Cell[i  ][j  ].A*
						  SolnBlk_Original_NW.U[i  ][j  ] +
						  SolnBlk_Original_NW.Grid.Cell[i+1][j  ].A*
						  SolnBlk_Original_NW.U[i+1][j  ] + 
						  SolnBlk_Original_NW.Grid.Cell[i  ][j+1].A*
						  SolnBlk_Original_NW.U[i  ][j+1] +
						  SolnBlk_Original_NW.Grid.Cell[i+1][j+1].A*
						  SolnBlk_Original_NW.U[i+1][j+1]) /
						 (SolnBlk_Original_NW.Grid.Cell[i  ][j  ].A +
						  SolnBlk_Original_NW.Grid.Cell[i+1][j  ].A +
						  SolnBlk_Original_NW.Grid.Cell[i  ][j+1].A +
						  SolnBlk_Original_NW.Grid.Cell[i+1][j+1].A) );
	//                                        SolnBlk_Coarse.Grid.Cell[i_coarse][j_coarse].A;
      } /* endfor */
    } /* endfor */

    for ( j = SolnBlk_Original_NW.JCl-SolnBlk_Original_NW.Nghost; 
	  j <= SolnBlk_Original_NW.JCu+SolnBlk_Original_NW.Nghost ; j += 2 ) {
      j_coarse = (j-SolnBlk_Original_NW.JCl)/2+
	(SolnBlk_Coarse.JCu-SolnBlk_Coarse.JCl+1)/2+SolnBlk_Coarse.JCl;
      SolnBlk_Coarse.UoW[j_coarse] = SolnBlk_Original_NW.UoW[j];
      if (j == SolnBlk_Original_NW.JCu+SolnBlk_Original_NW.Nghost) {
	SolnBlk_Coarse.UoW[j_coarse+1] = SolnBlk_Original_NW.UoW[j];
      } /* endif */
    } /* endfor */
    
    for ( i = SolnBlk_Original_NW.ICl-SolnBlk_Original_NW.Nghost ; 
	  i <= SolnBlk_Original_NW.ICu+SolnBlk_Original_NW.Nghost ; i += 2) {
      i_coarse = (i-SolnBlk_Original_NW.ICl)/2+SolnBlk_Coarse.ICl;
      SolnBlk_Coarse.UoN[i_coarse] = SolnBlk_Original_NW.UoN[i];
      if (i == SolnBlk_Original_NW.ICl-SolnBlk_Original_NW.Nghost) {
	SolnBlk_Coarse.UoN[i_coarse-1] = SolnBlk_Original_NW.UoN[i];
      } /* endif */
    } /* endfor */

      // North-east corner fine block:
    for ( j = SolnBlk_Original_NE.JCl; j <= SolnBlk_Original_NE.JCu ; j += 2 ) {
      for ( i = SolnBlk_Original_NE.ICl ; i <= SolnBlk_Original_NE.ICu ; i += 2 ) {
	i_coarse = (i-SolnBlk_Original_NE.ICl)/2+
	  (SolnBlk_Coarse.ICu-SolnBlk_Coarse.ICl+1)/2+SolnBlk_Coarse.ICl;
	j_coarse = (j-SolnBlk_Original_NE.JCl)/2+
	  (SolnBlk_Coarse.JCu-SolnBlk_Coarse.JCl+1)/2+SolnBlk_Coarse.JCl;
	SolnBlk_Coarse.U[i_coarse][j_coarse] = ( (SolnBlk_Original_NE.Grid.Cell[i  ][j  ].A*
						  SolnBlk_Original_NE.U[i  ][j  ] +
						  SolnBlk_Original_NE.Grid.Cell[i+1][j  ].A*
						  SolnBlk_Original_NE.U[i+1][j  ] + 
						  SolnBlk_Original_NE.Grid.Cell[i  ][j+1].A*
						  SolnBlk_Original_NE.U[i  ][j+1] +
						  SolnBlk_Original_NE.Grid.Cell[i+1][j+1].A*
						  SolnBlk_Original_NE.U[i+1][j+1]) /
						 (SolnBlk_Original_NE.Grid.Cell[i  ][j  ].A +
						  SolnBlk_Original_NE.Grid.Cell[i+1][j  ].A +
						  SolnBlk_Original_NE.Grid.Cell[i  ][j+1].A +
						  SolnBlk_Original_NE.Grid.Cell[i+1][j+1].A) );
	//                                        SolnBlk_Coarse.Grid.Cell[i_coarse][j_coarse].A;
      } /* endfor */
    } /* endfor */

    for ( j = SolnBlk_Original_NE.JCl-SolnBlk_Original_NE.Nghost; 
	  j <= SolnBlk_Original_NE.JCu+SolnBlk_Original_NE.Nghost ; j += 2 ) {
      j_coarse = (j-SolnBlk_Original_NE.JCl)/2+
	(SolnBlk_Coarse.JCu-SolnBlk_Coarse.JCl+1)/2+SolnBlk_Coarse.JCl;
      SolnBlk_Coarse.UoE[j_coarse] = SolnBlk_Original_NE.UoE[j];
      if (j == SolnBlk_Original_NE.JCu+SolnBlk_Original_NE.Nghost) {
	SolnBlk_Coarse.UoE[j_coarse+1] = SolnBlk_Original_NE.UoE[j];
      } /* endif */
    } /* endfor */

    for ( i = SolnBlk_Original_NE.ICl-SolnBlk_Original_NE.Nghost ; 
	  i <= SolnBlk_Original_NE.ICu+SolnBlk_Original_NE.Nghost ; i += 2) {
      i_coarse = (i-SolnBlk_Original_NE.ICl)/2+
	(SolnBlk_Coarse.ICu-SolnBlk_Coarse.ICl+1)/2+SolnBlk_Coarse.ICl;
      SolnBlk_Coarse.UoN[i_coarse] = SolnBlk_Original_NE.UoN[i];
      if (i == SolnBlk_Original_NE.ICu+SolnBlk_Original_NE.Nghost) {
	SolnBlk_Coarse.UoN[i_coarse+1] = SolnBlk_Original_NE.UoN[i];
      } /* endif */
    } /* endfor */

    // Set the reference values for the boundary states
    // This approach might not always give the proper reference values.
    SolnBlk_Coarse.Set_Reference_Values_For_Boundary_States(SolnBlk_Original_NW.Ref_State_BC_North,
							    SolnBlk_Original_SE.Ref_State_BC_South,
							    SolnBlk_Original_NE.Ref_State_BC_East,
							    SolnBlk_Original_SW.Ref_State_BC_West);
  } /* endif */

  // Restriction of solution block was successful.
  return 0;

}

/******************************************************//**
 * Routine: Output_Tecplot                              
 *                                                      
 * Writes the solution values at the nodes of the       
 * specified quadrilateral solution block to the        
 * specified output stream suitable for plotting with   
 * TECPLOT.                                             
 *                                                      
 ********************************************************/
void Output_Tecplot(AdvectDiffuse2D_Quad_Block &SolnBlk,
                    AdvectDiffuse2D_Input_Parameters &IP,
                    const int Number_of_Time_Steps,
                    const double &Time,
                    const int Block_Number,
                    const int Output_Title,
	            ostream &Out_File) {

  int i, j;
  AdvectDiffuse2D_State U_node;
  Vector2D Node;

  /* Ensure boundary conditions are updated before
     evaluating solution at the nodes. */
    
  BCs(SolnBlk,IP);

  /* Output node solution data. */

  Out_File << setprecision(14);
  if (Output_Title) {
    Out_File << "TITLE = \"" << CFFC_Name() << ": 2D Advection Diffusion Equation Solution, "
	     << "Time Step/Iteration Level = " << Number_of_Time_Steps
	     << ", Time = " << Time
	     << "\"" << "\n"
	     << "VARIABLES = \"x\" \\ \n"
	     << "\"y\" \\ \n"
	     << "\"u\" \\ \n"
	     << "\"Vx\" \\ \n"
	     << "\"Vy\" \\ \n"
	     << "\"k\" \\ \n"
	     << "\"s\" \\ \n";
    if (SolnBlk.ExactSoln->IsExactSolutionSet()){
      Out_File << "\"ExactSoln\" \\ \n";
    }
 
    Out_File << "ZONE T =  \"Block Number = " << Block_Number
	     << "\" \\ \n"
	     << "I = " << SolnBlk.Grid.INu - SolnBlk.Grid.INl + 1 << " \\ \n"
	     << "J = " << SolnBlk.Grid.JNu - SolnBlk.Grid.JNl + 1 << " \\ \n"
	     << "F = POINT \n";
  } else {
    Out_File << "ZONE T =  \"Block Number = " << Block_Number
	     << "\" \\ \n"
	     << "I = " << SolnBlk.Grid.INu - SolnBlk.Grid.INl + 1 << " \\ \n"
	     << "J = " << SolnBlk.Grid.JNu - SolnBlk.Grid.JNl + 1 << " \\ \n"
	     << "F = POINT \n";
  } /* endif */

  for ( j  = SolnBlk.Grid.JNl ; j <= SolnBlk.Grid.JNu ; ++j ) {
    for ( i = SolnBlk.Grid.INl ; i <= SolnBlk.Grid.INu ; ++i ) {
      U_node = SolnBlk.Un(i, j);
      Node = SolnBlk.Grid.Node[i][j].X;
      Out_File << " " << Node << U_node 
	       << " " << SolnBlk.U[i][j].V(Node.x,Node.y)
	       << " " << SolnBlk.U[i][j].k(Node.x,Node.y,U_node[1]) 
	       << " " << source(Node.x,Node.y,U_node);
      if (SolnBlk.ExactSoln->IsExactSolutionSet()){
	Out_File << " " << SolnBlk.ExactSoln->Solution(Node.x,Node.y);
      }
      Out_File << "\n";
      Out_File.unsetf(ios::scientific);
    } /* endfor */
  } /* endfor */
  Out_File << setprecision(6);
    
}

/******************************************************//**
 * Routine: Output_Cells_Tecplot                        
 *                                                      
 * Writes the cell centred solution values of the       
 * specified quadrilateral solution block to the        
 * specified output stream suitable for plotting with   
 * TECPLOT.                                             
 *                                                      
 ********************************************************/
void Output_Cells_Tecplot(AdvectDiffuse2D_Quad_Block &SolnBlk,
                          AdvectDiffuse2D_Input_Parameters &IP,
                          const int Number_of_Time_Steps,
                          const double &Time,
                          const int Block_Number,
                          const int Output_Title,
	                  ostream &Out_File) {

  int i, j;
  Vector2D Node;

  Out_File << setprecision(14);
  if (Output_Title) {
    Out_File << "TITLE = \"" << CFFC_Name() << ": 2D Advection Diffusion Equation Solution, "
	     << "Time Step/Iteration Level = " << Number_of_Time_Steps
	     << ", Time = " << Time
	     << "\"" << "\n"
	     << "VARIABLES = \"x\" \\ \n"
	     << "\"y\" \\ \n"
	     << "\"u\" \\ \n"
	     << "\"Vx\" \\ \n"
	     << "\"Vy\" \\ \n"
	     << "\"k\" \\ \n"
	     << "\"s\" \\ \n";
    if (SolnBlk.ExactSoln->IsExactSolutionSet()){
      Out_File << "\"ExactSoln\" \\ \n";
    }
    Out_File << "ZONE T =  \"Block Number = " << Block_Number
	     << "\" \\ \n"
	     << "I = " << SolnBlk.Grid.ICu - SolnBlk.Grid.ICl + 2*SolnBlk.Nghost + 1 << " \\ \n"
	     << "J = " << SolnBlk.Grid.JCu - SolnBlk.Grid.JCl + 2*SolnBlk.Nghost + 1 << " \\ \n"
	     << "F = POINT \n";
  } else {
    Out_File << "ZONE T =  \"Block Number = " << Block_Number
	     << "\" \\ \n"
	     << "I = " << SolnBlk.Grid.ICu - SolnBlk.Grid.ICl + 2*SolnBlk.Nghost + 1 << " \\ \n"
	     << "J = " << SolnBlk.Grid.JCu - SolnBlk.Grid.JCl + 2*SolnBlk.Nghost + 1 << " \\ \n"
	     << "F = POINT \n";
  } /* endif */

  for ( j  = SolnBlk.JCl-SolnBlk.Nghost ; j <= SolnBlk.JCu+SolnBlk.Nghost ; ++j ) {
    for ( i = SolnBlk.ICl-SolnBlk.Nghost ; i <= SolnBlk.ICu+SolnBlk.Nghost ; ++i ) {
      Node = SolnBlk.Grid.Cell[i][j].Xc;
      Out_File << " " << SolnBlk.Grid.Cell[i][j].Xc << SolnBlk.U[i][j]
	       << " " << SolnBlk.U[i][j].V(Node.x,Node.y)
	       << " " << SolnBlk.U[i][j].k(Node.x,Node.y,SolnBlk.U[i][j][1])
	       << " " << source(Node.x,Node.y,SolnBlk.U[i][j]);
      if (SolnBlk.ExactSoln->IsExactSolutionSet()){
	Out_File << " " << SolnBlk.ExactSoln->Solution(Node.x,Node.y);
      }
      Out_File << "\n";
    } /* endfor */
  } /* endfor */
  Out_File << setprecision(6);
    
}

/******************************************************//**
 * Routine: Output_Nodes_Tecplot                        
 *                                                      
 * Writes the node values of the specified              
 * quadrilateral solution block to the specified output 
 * stream suitable for plotting with TECPLOT.           
 *                                                      
 ********************************************************/
void Output_Nodes_Tecplot(AdvectDiffuse2D_Quad_Block &SolnBlk,
                          const int Number_of_Time_Steps,
                          const double &Time,
                          const int Block_Number,
                          const int Output_Title,
	                  ostream &Out_File) {

  int i, j;

  Out_File << setprecision(14);
  if (Output_Title) {
    Out_File << "TITLE = \"" << CFFC_Name() << ": 2D Advection Diffusion Node Values, "
	     << "Time Step/Iteration Level = " << Number_of_Time_Steps
	     << ", Time = " << Time
	     << "\"" << "\n"
	     << "VARIABLES = \"x\" \\ \n"
	     << "\"y\" \\ \n";
  } /* endif */

  Out_File << "ZONE T =  \"Block Number = " << Block_Number
	   << "\" \\ \n"
	   << "I = " << SolnBlk.Grid.INu - SolnBlk.Grid.INl + 1 + 2*SolnBlk.Nghost << " \\ \n"
	   << "J = " << SolnBlk.Grid.JNu - SolnBlk.Grid.JNl + 1 + 2*SolnBlk.Nghost << " \\ \n"
	   << "F = POINT \\ \n";

  Out_File.setf(ios::scientific);
  for (j = SolnBlk.Grid.JNl - SolnBlk.Nghost; j <= SolnBlk.Grid.JNu + SolnBlk.Nghost; ++j) {
    for (i = SolnBlk.Grid.INl - SolnBlk.Nghost; i <= SolnBlk.Grid.INu + SolnBlk.Nghost; ++i) {
      Out_File << " " << SolnBlk.Grid.Node[i][j].X << endl;
    }
  }
  Out_File.unsetf(ios::scientific);
  Out_File << setprecision(6);

}

/******************************************************//**
 * Routine: ICs                                         
 *                                                      
 * Assigns initial conditions and data to the           
 * solution variables of the specified quadrilateral    
 * solution block.                                      
 *                                                      
 ********************************************************/
void ICs(AdvectDiffuse2D_Quad_Block &SolnBlk,
	 const AdvectDiffuse2D_Input_Parameters &IP,
         AdvectDiffuse2D_State *Uo) {

  int i, j, k;
  double r0;
  Vector2D x0;
  AdvectDiffuse2D_State Ul, Ur;

  /* Assign the initial data for the IVP of interest. */

  switch(IP.i_ICs) {
  case IC_CONSTANT :
  case IC_UNIFORM :
    // Set the solution state to the initial state Uo[1].
    for (j  = SolnBlk.JCl-SolnBlk.Nghost ; j <= SolnBlk.JCu+SolnBlk.Nghost ; ++j ) {
      for ( i = SolnBlk.ICl-SolnBlk.Nghost ; i <= SolnBlk.ICu+SolnBlk.Nghost ; ++i ) {
	SolnBlk.U[i][j] = Uo[1];
      } /* endfor */
    } /* endfor */
    break;
  case IC_RIEMANN :
  case IC_RIEMANN_XDIR :
    // Set Riemann initial data 
    // (2-state initial data, left to right).
    Ur = Uo[0];
    Ul = Uo[1];
    for (j  = SolnBlk.JCl-SolnBlk.Nghost ; j <= SolnBlk.JCu+SolnBlk.Nghost ; ++j ) {
      for ( i = SolnBlk.ICl-SolnBlk.Nghost ; i <= SolnBlk.ICu+SolnBlk.Nghost ; ++i ) {
	if (SolnBlk.Grid.Cell[i][j].Xc.x <= ZERO) {
	  SolnBlk.U[i][j] = Ul;
	} else {
	  SolnBlk.U[i][j] = Ur;	     
	} /* end if */
      } /* endfor */
    } /* endfor */
    break;
  case IC_RIEMANN_YDIR :
    // Set Riemann initial data 
    // (2-state initial data, top to bottom).
    Ur = Uo[0];
    Ul = Uo[1];
    for (j  = SolnBlk.JCl-SolnBlk.Nghost ; j <= SolnBlk.JCu+SolnBlk.Nghost ; ++j ) {
      for ( i = SolnBlk.ICl-SolnBlk.Nghost ; i <= SolnBlk.ICu+SolnBlk.Nghost ; ++i ) {
	if (SolnBlk.Grid.Cell[i][j].Xc.y <= ZERO) {
	  SolnBlk.U[i][j] = Ul;
	} else {
	  SolnBlk.U[i][j] = Ur;	     
	} /* end if */
      } /* endfor */
    } /* endfor */
    break;
  case IC_SQUARE_BOX_IVP :
    // Set initial data for square box IVP with
    // non-zero solution in square region, x in [-0.75, -0.25], y in [-0.75, -0.25].
    Ur = Uo[1];
    Ul = Uo[0];
    for (j  = SolnBlk.JCl-SolnBlk.Nghost ; j <= SolnBlk.JCu+SolnBlk.Nghost ; ++j ) {
      for ( i = SolnBlk.ICl-SolnBlk.Nghost ; i <= SolnBlk.ICu+SolnBlk.Nghost ; ++i ) {
	if (SolnBlk.Grid.Cell[i][j].Xc.x >= -0.75 &&
	    SolnBlk.Grid.Cell[i][j].Xc.x <= -0.25 &&
	    SolnBlk.Grid.Cell[i][j].Xc.y >= -0.75 &&
	    SolnBlk.Grid.Cell[i][j].Xc.y <= -0.25) {
	  SolnBlk.U[i][j] = Ul;
	} else {
	  SolnBlk.U[i][j] = Ur;	     
	} /* end if */
      } /* endfor */
    } /* endfor */
    break;
  case IC_CIRCULAR_BOX_IVP :
    // Set initial data for circular box IVP with
    // non-zero solution in circular region, X-x0 < 0.50 and
    // x0 = [-0.50, -0.50]
    Ur = Uo[1];
    Ul = Uo[0];
    x0 = Vector2D(-0.50,-0.50);
    for (j  = SolnBlk.JCl-SolnBlk.Nghost ; j <= SolnBlk.JCu+SolnBlk.Nghost ; ++j ) {
      for ( i = SolnBlk.ICl-SolnBlk.Nghost ; i <= SolnBlk.ICu+SolnBlk.Nghost ; ++i ) {
	r0 = abs(SolnBlk.Grid.Cell[i][j].Xc-x0);
	if (r0 <= 0.50) {
	  SolnBlk.U[i][j] = Ul;
	} else {
	  SolnBlk.U[i][j] = Ur;	     
	} /* end if */
      } /* endfor */
    } /* endfor */
    break;
  case IC_EXACT_SOLUTION :
    // Set the solution state by calculating the cell average values with integration of the exact solution
    // Use the ExactSoln pointer to access the exact solution
    if (IP.ExactSoln->IsExactSolutionSet()) {
      for (j  = SolnBlk.JCl-SolnBlk.Nghost ; j <= SolnBlk.JCu+SolnBlk.Nghost ; ++j ) {
 	for ( i = SolnBlk.ICl-SolnBlk.Nghost ; i <= SolnBlk.ICu+SolnBlk.Nghost ; ++i ) {
	  Ul.u = 
	    SolnBlk.Grid.Integration.IntegrateFunctionOverCell(i,j,
							       wrapped_member_function(IP.ExactSoln,
										       &AdvectDiffuse2D_ExactSolutions::Solution,
										       Ul.u),
							       wrapped_member_function(IP.ExactSoln,
										       &AdvectDiffuse2D_ExactSolutions::
										       XDependencyIntegrated_Solution,
										       Ul.u),
							       IP.Exact_Integration_Digits,Ul.u)/SolnBlk.Grid.Cell[i][j].A;
	  SolnBlk.U[i][j] = Ul;
	} /* endfor */
      } /* endfor */
    } else {
      // There is no exact solution set for this problem
      throw runtime_error("ICs() ERROR! No exact solution has been set!");
    }
    break;
  case IC_INTERIOR_UNIFORM_GHOSTCELLS_EXACT :
    if (IP.ExactSoln->IsExactSolutionSet()) {
      for (j  = SolnBlk.JCl-SolnBlk.Nghost ; j <= SolnBlk.JCu+SolnBlk.Nghost ; ++j ) {
	for ( i = SolnBlk.ICl-SolnBlk.Nghost ; i <= SolnBlk.ICu+SolnBlk.Nghost ; ++i ) {
	  if (i<SolnBlk.ICl || i>SolnBlk.ICu || j<SolnBlk.JCl || j>SolnBlk.JCu) {
	    // Set the solution state of the ghost cells to average values calculated with
	    // integration of the exact solution
	    Ul.u = 
	      SolnBlk.Grid.Integration.IntegrateFunctionOverCell(i,j,
								 wrapped_member_function(IP.ExactSoln,
											 &AdvectDiffuse2D_ExactSolutions::Solution,
											 Ul.u),
								 IP.Exact_Integration_Digits,Ul.u)/SolnBlk.Grid.Cell[i][j].A;
	    SolnBlk.U[i][j] = Ul;
	  } else {
	    // Set the solution state of the interior cells to the initial state Uo[0].
	    SolnBlk.U[i][j] = Uo[1];
	  }
	} /* endfor */
      } /* endfor */

    } else {
      // There is no exact solution set for this problem
      throw runtime_error("ICs() ERROR! No exact solution has been set!");
    }
    break;

  default:
    // Set the solution state to the initial state Uo[0].
    for (j  = SolnBlk.JCl-SolnBlk.Nghost ; j <= SolnBlk.JCu+SolnBlk.Nghost ; ++j ) {
      for ( i = SolnBlk.ICl-SolnBlk.Nghost ; i <= SolnBlk.ICu+SolnBlk.Nghost ; ++i ) {
	SolnBlk.U[i][j] = Uo[0];
      } /* endfor */
    } /* endfor */
    break;
  } /* endswitch */

  
  /* Assign the reference state data for the Boundary Conditions of interest. */
  SolnBlk.Set_Boundary_Reference_States_Based_On_Input(IP);

}


/******************************************************//**
 * Routine: BCs                                         
 *                                                      
 * Apply boundary conditions at boundaries of the       
 * specified quadrilateral solution block.              
 * 
 ********************************************************/
void BCs(AdvectDiffuse2D_Quad_Block &SolnBlk,
	 const AdvectDiffuse2D_Input_Parameters &IP) {

  int i(0), j(0);
  int ghost;
  double Vn;
  AdvectDiffuse2D_State dUdn;
  double dx_normal;
  

  for ( j = SolnBlk.JCl-SolnBlk.Nghost ; j <= SolnBlk.JCu+SolnBlk.Nghost ; ++j ) {
    // Prescribe West boundary conditions.
    if ( (j >= SolnBlk.JCl && j <= SolnBlk.JCu) ||                                  // <-- affects W boundary cells
	 (j < SolnBlk.JCl && (SolnBlk.Grid.BCtypeS[SolnBlk.ICl-1] == BC_NONE ) ) || // <-- affects SW corner cells
	 (j > SolnBlk.JCu && (SolnBlk.Grid.BCtypeN[SolnBlk.ICl-1] == BC_NONE ) ) )  // <-- affects NW corner cells 
      {														     
	switch(SolnBlk.Grid.BCtypeW[j]) {	
	case BC_FROZEN :
	  // Equivalent to BC_NONE in the sense that it leaves 
	  // the ghost cell values unchanged and 
	  // uses the ghost cell reconstruction to 
	  // compute the inter-cellular state.
	case BC_NONE :
	  break;

	case BC_PERIODIC :
	  for( ghost = 1; ghost <= SolnBlk.Nghost; ++ghost){
	    SolnBlk.U[SolnBlk.ICl-ghost][j] = SolnBlk.U[SolnBlk.ICu-ghost+1][j];
	  }
	  break;

	case BC_INFLOW :	// Inflow BC is treated as a Dirichlet BC. Make sure that UoW has the right inflow data!
	  if (SolnBlk.Grid.BndWestSpline.getFluxCalcMethod() == ReconstructionBasedFlux){
	    for( ghost = 1; ghost <= SolnBlk.Nghost; ++ghost){
	      SolnBlk.U[SolnBlk.ICl-ghost][j] = SolnBlk.UoW[j];
	    }
	    break;
	  }

	case BC_DIRICHLET :	// Use UoW as reference value
	  // Compute the dx in the normal direction between the Gauss point and the first interior cell centroid.
	  dx_normal = ( SolnBlk.Grid.xfaceW(SolnBlk.ICl,j) - 
			SolnBlk.Grid.CellCentroid(SolnBlk.ICl,j) ) * SolnBlk.Grid.nfaceW(SolnBlk.ICl,j);
	  // Estimate the solution gradient in the normal direction
	  dUdn = (SolnBlk.UoW[j] - SolnBlk.U[SolnBlk.ICl][j])/dx_normal;

	  // Fill in the ghost cells.
	  for( ghost = 1; ghost <= SolnBlk.Nghost; ++ghost){
	    dx_normal = ( SolnBlk.Grid.CellCentroid(SolnBlk.ICl-ghost,j) -
			  SolnBlk.Grid.xfaceW(SolnBlk.ICl,j)) * SolnBlk.Grid.nfaceW(SolnBlk.ICl,j);
	    SolnBlk.U[SolnBlk.ICl-ghost][j] = ( SolnBlk.UoW[j] + dUdn * dx_normal );
	  }
	  break;

	case BC_NEUMANN :
	  // The ghost cell values are set based on the normal gradient provided by the user (UoW[j])
	  // and the average solution of the first interior cell.
	  for( ghost = 1; ghost <= SolnBlk.Nghost; ++ghost){
	    dx_normal = ( SolnBlk.Grid.CellCentroid(SolnBlk.ICl-ghost,j) -
			  SolnBlk.Grid.xfaceW(SolnBlk.ICl,j)) * SolnBlk.Grid.nfaceW(SolnBlk.ICl,j);

	    SolnBlk.U[SolnBlk.ICl-ghost][j] = ( SolnBlk.U[SolnBlk.ICl][j] + SolnBlk.UoW[j] * dx_normal );
	  }
	  break;

	case BC_FARFIELD :
	  /* Farfield BC is implemented differently for flows that 
	     enter the domain than for flows that leave the domain.
	     Whether the flow enters of leaves the domain is decided based on
	     the normal component of the velocity at the face midpoint.
	     --> If the flow enters the domain then the reference data is used.
	     --> If the flow leaves the domain then an interior extrapolated value is used. */

	  // Compute the normal velocity
	  Vn = dot(SolnBlk.VelocityAtLocation(SolnBlk.Grid.xfaceW(SolnBlk.ICl,j)),
		   SolnBlk.Grid.nfaceW(SolnBlk.ICl,j));
	  
	  if (Vn <= ZERO){
	    // The flow enters the domain
	    // Use UoW as reference value

	    if (SolnBlk.Grid.BndWestSpline.getFluxCalcMethod() == ReconstructionBasedFlux) {
	      // Fill in the ghost cells.
	      for( ghost = 1; ghost <= SolnBlk.Nghost; ++ghost){
		SolnBlk.U[SolnBlk.ICl-ghost][j] = SolnBlk.UoW[j];
	      }
	      
	    } else {

	      // Compute the dx in the normal direction between the Gauss point and the first interior cell centroid.
	      dx_normal = ( SolnBlk.Grid.xfaceW(SolnBlk.ICl,j) - 
			    SolnBlk.Grid.CellCentroid(SolnBlk.ICl,j) ) * SolnBlk.Grid.nfaceW(SolnBlk.ICl,j);
	      // Estimate the solution gradient in the normal direction
	      dUdn = (SolnBlk.UoW[j] - SolnBlk.U[SolnBlk.ICl][j])/dx_normal;
	      
	      // Fill in the ghost cells.
	      for( ghost = 1; ghost <= SolnBlk.Nghost; ++ghost){
		dx_normal = ( SolnBlk.Grid.CellCentroid(SolnBlk.ICl-ghost,j) -
			      SolnBlk.Grid.xfaceW(SolnBlk.ICl,j)) * SolnBlk.Grid.nfaceW(SolnBlk.ICl,j);
		SolnBlk.U[SolnBlk.ICl-ghost][j] = ( SolnBlk.UoW[j] + dUdn * dx_normal );
	      }
	    }

	  } else {
	    // The flow leaves the domain
	    // Impose constant extrapolation
	    for( ghost = 1; ghost <= SolnBlk.Nghost; ++ghost){
	      SolnBlk.U[SolnBlk.ICl-ghost][j] = SolnBlk.U[SolnBlk.ICl][j];
	    }
	  }
	  break;

	case BC_SYMMETRY_PLANE :
	  throw runtime_error("BCs() ERROR! Symmetry plane BC hasn't been implemented yet!");
	  break;

	case BC_EXTRAPOLATE :
	case BC_LINEAR_EXTRAPOLATION :
	  throw runtime_error("BCs() ERROR! Linear extrapolation BC hasn't been implemented yet!");
	  break;

	case BC_OUTFLOW :
	  // Impose zero derivative in the i-direction at the boundary
	  // This is equivalent with a constant extrapolation
	case BC_CONSTANT_EXTRAPOLATION :
	  for( ghost = 1; ghost <= SolnBlk.Nghost; ++ghost){
	    SolnBlk.U[SolnBlk.ICl-ghost][j] = SolnBlk.U[SolnBlk.ICl][j];
	  }
	  break;

	case BC_EXACT_SOLUTION :
	  // Leave the ghost cell values unchanged.
	  break;

	default:		
	  // Impose constant extrapolation
	  for( ghost = 1; ghost <= SolnBlk.Nghost; ++ghost){
	    SolnBlk.U[SolnBlk.ICl-ghost][j] = SolnBlk.U[SolnBlk.ICl][j];
	  }
	  break;
	} /* endswitch */
      } /* endif */


    // Prescribe East boundary conditions.
    if ( (j >= SolnBlk.JCl && j <= SolnBlk.JCu) ||                                   // <-- affects E boundary cells
	 (j < SolnBlk.JCl && (SolnBlk.Grid.BCtypeS[SolnBlk.ICu+1] == BC_NONE ) ) ||  // <-- affects SE corner cells
	 (j > SolnBlk.JCu && (SolnBlk.Grid.BCtypeN[SolnBlk.ICu+1] == BC_NONE ) ) )   // <-- affects NE corner cells
      {
	switch(SolnBlk.Grid.BCtypeE[j]) {
	case BC_FROZEN :
	  // Equivalent to BC_NONE in the sense that it leaves 
	  // the ghost cell values unchanged and 
	  // uses the ghost cell reconstruction to 
	  // compute the inter-cellular state.
	case BC_NONE :
	  break;
	  
	case BC_PERIODIC :
	  for( ghost = 1; ghost <= SolnBlk.Nghost; ++ghost){
	    SolnBlk.U[SolnBlk.ICu+ghost][j] = SolnBlk.U[SolnBlk.ICl+ghost-1][j];
	  }
	  break;
	  
	case BC_INFLOW :	// Inflow BC is treated as a Dirichlet BC. Make sure that UoE has the right inflow data!
	  if (SolnBlk.Grid.BndEastSpline.getFluxCalcMethod() == ReconstructionBasedFlux){
	    for( ghost = 1; ghost <= SolnBlk.Nghost; ++ghost){
	      SolnBlk.U[SolnBlk.ICu+ghost][j] = SolnBlk.UoE[j];
	    }
	    break;
	  }

	case BC_DIRICHLET :	// Use UoE as reference value
	  // Compute the dx in the normal direction between the Gauss point and the first interior cell centroid.
	  dx_normal = ( SolnBlk.Grid.xfaceE(SolnBlk.ICu,j) - 
			SolnBlk.Grid.CellCentroid(SolnBlk.ICu,j) ) * SolnBlk.Grid.nfaceE(SolnBlk.ICu,j);
	  // Estimate the solution gradient in the normal direction
	  dUdn = (SolnBlk.UoE[j] - SolnBlk.U[SolnBlk.ICu][j])/dx_normal;

	  // Fill in the ghost cells.
	  for( ghost = 1; ghost <= SolnBlk.Nghost; ++ghost){
	    dx_normal = ( SolnBlk.Grid.CellCentroid(SolnBlk.ICu+ghost,j) -
			  SolnBlk.Grid.xfaceW(SolnBlk.ICu,j)) * SolnBlk.Grid.nfaceE(SolnBlk.ICu,j);
	    SolnBlk.U[SolnBlk.ICu+ghost][j] = ( SolnBlk.UoE[j] + dUdn * dx_normal );
	  }
	  break;
	  
	case BC_NEUMANN :
	  // The ghost cell values are set based on the normal gradient provided by the user (UoE[j])
	  // and the average solution of the first interior cell.
	  for( ghost = 1; ghost <= SolnBlk.Nghost; ++ghost){
	    dx_normal = ( SolnBlk.Grid.CellCentroid(SolnBlk.ICu+ghost,j) -
			  SolnBlk.Grid.xfaceW(SolnBlk.ICu,j)) * SolnBlk.Grid.nfaceE(SolnBlk.ICu,j);

	    SolnBlk.U[SolnBlk.ICu+ghost][j] = ( SolnBlk.U[SolnBlk.ICu][j] + SolnBlk.UoE[j] * dx_normal );
	  }
	  break;

	case BC_FARFIELD :
	  /* Farfield BC is implemented differently for flows that 
	     enter the domain than for flows that leave the domain.
	     Whether the flow enters of leaves the domain is decided based on
	     the normal component of the velocity at the face midpoint.
	     --> If the flow enters the domain then the reference data is used.
	     --> If the flow leaves the domain then an interior extrapolated value is used. */

	  // Compute the normal velocity
	  Vn = dot(SolnBlk.VelocityAtLocation(SolnBlk.Grid.xfaceE(SolnBlk.ICu,j)),
		   SolnBlk.Grid.nfaceE(SolnBlk.ICu,j));
	  
	  if (Vn <= ZERO){
	    // The flow enters the domain
	    // Use UoE as reference value

	    if (SolnBlk.Grid.BndEastSpline.getFluxCalcMethod() == ReconstructionBasedFlux) {
	      // Fill in the ghost cells.
	      for( ghost = 1; ghost <= SolnBlk.Nghost; ++ghost){
		SolnBlk.U[SolnBlk.ICu+ghost][j] = SolnBlk.UoE[j];
	      }

	    } else {

	      // Compute the dx in the normal direction between the Gauss point and the first interior cell centroid.
	      dx_normal = ( SolnBlk.Grid.xfaceE(SolnBlk.ICu,j) - 
			    SolnBlk.Grid.CellCentroid(SolnBlk.ICu,j) ) * SolnBlk.Grid.nfaceE(SolnBlk.ICu,j);
	      // Estimate the solution gradient in the normal direction
	      dUdn = (SolnBlk.UoE[j] - SolnBlk.U[SolnBlk.ICu][j])/dx_normal;
	    
	      // Fill in the ghost cells.
	      for( ghost = 1; ghost <= SolnBlk.Nghost; ++ghost){
		dx_normal = ( SolnBlk.Grid.CellCentroid(SolnBlk.ICu+ghost,j) -
			      SolnBlk.Grid.xfaceW(SolnBlk.ICu,j)) * SolnBlk.Grid.nfaceE(SolnBlk.ICu,j);
		SolnBlk.U[SolnBlk.ICu+ghost][j] = ( SolnBlk.UoE[j] + dUdn * dx_normal );
	      }
	    }

	  } else {
	    // The flow leaves the domain
	    // Impose constant extrapolation
	    for( ghost = 1; ghost <= SolnBlk.Nghost; ++ghost){
	      SolnBlk.U[SolnBlk.ICu+ghost][j] = SolnBlk.U[SolnBlk.ICu][j];
	    }
	  }
	  break;

	case BC_SYMMETRY_PLANE :
	  throw runtime_error("BCs() ERROR! Symmetry plane BC hasn't been implemented yet!");
	  break;

	case BC_EXTRAPOLATE :
	case BC_LINEAR_EXTRAPOLATION :
	  throw runtime_error("BCs() ERROR! Linear extrapolation BC hasn't been implemented yet!");
	  break;

	case BC_OUTFLOW :
	  // Impose zero derivative in the i-direction at the boundary
	  // This is equivalent with a constant extrapolation
	case BC_CONSTANT_EXTRAPOLATION :
	  for( ghost = 1; ghost <= SolnBlk.Nghost; ++ghost){
	    SolnBlk.U[SolnBlk.ICu+ghost][j] = SolnBlk.U[SolnBlk.ICu][j];
	  }
	  break;

	case BC_EXACT_SOLUTION :
	  // Leave the ghost cell values unchanged.
	  break;

	default:		
	  // Impose constant extrapolation
	  for( ghost = 1; ghost <= SolnBlk.Nghost; ++ghost){
	    SolnBlk.U[SolnBlk.ICu+ghost][j] = SolnBlk.U[SolnBlk.ICu][j];
	  }
	  break;
	} /* endswitch */
      } /* endif */
  } /* endfor */


  for ( i = SolnBlk.ICl-SolnBlk.Nghost ; i <= SolnBlk.ICu+SolnBlk.Nghost ; ++i ) {
    // Use the South and North BCs for the corner ghost cells

    // Prescribe South boundary conditions.
    switch(SolnBlk.Grid.BCtypeS[i]) {
    case BC_FROZEN :
      // Equivalent to BC_NONE in the sense that it leaves 
      // the ghost cell values unchanged and 
      // uses the ghost cell reconstruction to 
      // compute the inter-cellular state.
    case BC_NONE :
      break;
      
    case BC_PERIODIC :
      for( ghost = 1; ghost <= SolnBlk.Nghost; ++ghost){
	SolnBlk.U[i][SolnBlk.JCl-ghost] = SolnBlk.U[i][SolnBlk.JCu-ghost+1];
      }
      break;
      
    case BC_INFLOW :	// Inflow BC is treated as a Dirichlet BC. Make sure that UoS has the right inflow data!
      if (SolnBlk.Grid.BndSouthSpline.getFluxCalcMethod() == ReconstructionBasedFlux){
	for( ghost = 1; ghost <= SolnBlk.Nghost; ++ghost){
	  SolnBlk.U[i][SolnBlk.JCl-ghost] = SolnBlk.UoS[i];
	}
	break;
      }
      
    case BC_DIRICHLET :	// Use UoS as reference value
      // Compute the dx in the normal direction between the Gauss point and the first interior cell centroid.
      dx_normal = ( SolnBlk.Grid.xfaceS(i,SolnBlk.JCl) - 
		    SolnBlk.Grid.CellCentroid(i,SolnBlk.JCl) ) * SolnBlk.Grid.nfaceS(i,SolnBlk.JCl);
      // Estimate the solution gradient in the normal direction
      dUdn = (SolnBlk.UoS[i] - SolnBlk.U[i][SolnBlk.JCl])/dx_normal;
      
      // Fill in the ghost cells.
      for( ghost = 1; ghost <= SolnBlk.Nghost; ++ghost){
	dx_normal = ( SolnBlk.Grid.CellCentroid(i,SolnBlk.JCl-ghost) -
		      SolnBlk.Grid.xfaceS(i,SolnBlk.JCl)) * SolnBlk.Grid.nfaceS(i,SolnBlk.JCl);
	SolnBlk.U[i][SolnBlk.JCl-ghost] = ( SolnBlk.UoS[i] + dUdn * dx_normal );
      }
      break;
      
    case BC_NEUMANN :
      // The ghost cell values are set based on the normal gradient provided by the user (UoS[i])
      // and the average solution of the first interior cell.
      for( ghost = 1; ghost <= SolnBlk.Nghost; ++ghost){
	dx_normal = ( SolnBlk.Grid.CellCentroid(i,SolnBlk.JCl-ghost) -
		      SolnBlk.Grid.xfaceS(i,SolnBlk.JCl)) * SolnBlk.Grid.nfaceS(i,SolnBlk.JCl);

	SolnBlk.U[i][SolnBlk.JCl-ghost] = ( SolnBlk.U[i][SolnBlk.JCl] + SolnBlk.UoS[i] * dx_normal );
      }
      break;
      
    case BC_FARFIELD :
      /* Farfield BC is implemented differently for flows that 
	 enter the domain than for flows that leave the domain.
	 Whether the flow enters of leaves the domain is decided based on
	 the normal component of the velocity at the face midpoint.
	 --> If the flow enters the domain then the reference data is used.
	 --> If the flow leaves the domain then an interior extrapolated value is used. */

      // Compute the normal velocity
      Vn = dot(SolnBlk.VelocityAtLocation(SolnBlk.Grid.xfaceS(i,SolnBlk.JCl)),
	       SolnBlk.Grid.nfaceS(i,SolnBlk.JCl));
	  
      if (Vn <= ZERO){
	// The flow enters the domain
	// Use UoS as reference value

	if (SolnBlk.Grid.BndSouthSpline.getFluxCalcMethod() == ReconstructionBasedFlux) {
	  // Fill in the ghost cells.
	  for( ghost = 1; ghost <= SolnBlk.Nghost; ++ghost){
	    SolnBlk.U[i][SolnBlk.JCl-ghost] = SolnBlk.UoS[i];
	  }
	  
	} else {

	  // Compute the dx in the normal direction between the Gauss point and the first interior cell centroid.
	  dx_normal = ( SolnBlk.Grid.xfaceS(i,SolnBlk.JCl) - 
			SolnBlk.Grid.CellCentroid(i,SolnBlk.JCl) ) * SolnBlk.Grid.nfaceS(i,SolnBlk.JCl);
	  // Estimate the solution gradient in the normal direction
	  dUdn = (SolnBlk.UoS[i] - SolnBlk.U[i][SolnBlk.JCl])/dx_normal;
      
	  // Fill in the ghost cells.
	  for( ghost = 1; ghost <= SolnBlk.Nghost; ++ghost){
	    dx_normal = ( SolnBlk.Grid.CellCentroid(i,SolnBlk.JCl-ghost) -
			  SolnBlk.Grid.xfaceS(i,SolnBlk.JCl)) * SolnBlk.Grid.nfaceS(i,SolnBlk.JCl);
	    SolnBlk.U[i][SolnBlk.JCl-ghost] = ( SolnBlk.UoS[i] + dUdn * dx_normal );
	  }
	}

      } else {
	// The flow leaves the domain
	// Impose constant extrapolation
	for( ghost = 1; ghost <= SolnBlk.Nghost; ++ghost){
	  SolnBlk.U[i][SolnBlk.JCl-ghost] = SolnBlk.U[i][SolnBlk.JCl];
	}
      }
      break;
      
    case BC_SYMMETRY_PLANE :
      throw runtime_error("BCs() ERROR! Symmetry plane BC hasn't been implemented yet!");
      break;
      
    case BC_EXTRAPOLATE :
    case BC_LINEAR_EXTRAPOLATION :
      throw runtime_error("BCs() ERROR! Linear extrapolation BC hasn't been implemented yet!");
      break;
      
    case BC_OUTFLOW :
      // Impose zero derivative in the j-direction at the boundary
      // This is equivalent with a constant extrapolation
    case BC_CONSTANT_EXTRAPOLATION :
      for( ghost = 1; ghost <= SolnBlk.Nghost; ++ghost){
	SolnBlk.U[i][SolnBlk.JCl-ghost] = SolnBlk.U[i][SolnBlk.JCl];
      }
      break;

    case BC_EXACT_SOLUTION :
      // Leave the ghost cell values unchanged.
      break;

    default:		
      // Impose constant extrapolation
      for( ghost = 1; ghost <= SolnBlk.Nghost; ++ghost){
	SolnBlk.U[i][SolnBlk.JCl-ghost] = SolnBlk.U[i][SolnBlk.JCl];
      }
      break;
    } /* endswitch */

    
    // Prescribe North boundary conditions.
    switch(SolnBlk.Grid.BCtypeN[i]) {
    case BC_FROZEN :
      // Equivalent to BC_NONE in the sense that it leaves 
      // the ghost cell values unchanged and 
      // uses the ghost cell reconstruction to 
      // compute the inter-cellular state.    
    case BC_NONE :
      break;
      
    case BC_PERIODIC :
      for( ghost = 1; ghost <= SolnBlk.Nghost; ++ghost){
	SolnBlk.U[i][SolnBlk.JCu+ghost] = SolnBlk.U[i][SolnBlk.JCl+ghost-1];
      }
      break;
      
    case BC_INFLOW :	// Inflow BC is treated as a Dirichlet BC. Make sure that UoN has the right inflow data!
      if (SolnBlk.Grid.BndNorthSpline.getFluxCalcMethod() == ReconstructionBasedFlux){
	for( ghost = 1; ghost <= SolnBlk.Nghost; ++ghost){
	  SolnBlk.U[i][SolnBlk.JCu+ghost] = SolnBlk.UoN[i];
	}
	break;
      }

    case BC_DIRICHLET :	// Use UoN as reference value
      // Compute the dx in the normal direction between the Gauss point and the first interior cell centroid.
      dx_normal = ( SolnBlk.Grid.xfaceN(i,SolnBlk.JCu) - 
		    SolnBlk.Grid.CellCentroid(i,SolnBlk.JCu) ) * SolnBlk.Grid.nfaceN(i,SolnBlk.JCu);
      // Estimate the solution gradient in the normal direction
      dUdn = (SolnBlk.UoN[i] - SolnBlk.U[i][SolnBlk.JCu])/dx_normal;

      // Fill in the ghost cells.
      for( ghost = 1; ghost <= SolnBlk.Nghost; ++ghost){
	dx_normal = ( SolnBlk.Grid.CellCentroid(i,SolnBlk.JCu+ghost) -
		      SolnBlk.Grid.xfaceN(i,SolnBlk.JCu)) * SolnBlk.Grid.nfaceN(i,SolnBlk.JCu);
	SolnBlk.U[i][SolnBlk.JCu+ghost] = ( SolnBlk.UoN[i] + dUdn * dx_normal );
      }
      break;
      
    case BC_NEUMANN :
      // The ghost cell values are set based on the normal gradient provided by the user (UoN[i])
      // and the average solution of the first interior cell.
      for( ghost = 1; ghost <= SolnBlk.Nghost; ++ghost){
	dx_normal = ( SolnBlk.Grid.CellCentroid(i,SolnBlk.JCu+ghost) -
		      SolnBlk.Grid.xfaceN(i,SolnBlk.JCu)) * SolnBlk.Grid.nfaceN(i,SolnBlk.JCu);

	SolnBlk.U[i][SolnBlk.JCu+ghost] = ( SolnBlk.U[i][SolnBlk.JCu] + SolnBlk.UoN[i] * dx_normal );
      }
      break;
      
    case BC_FARFIELD :
      /* Farfield BC is implemented differently for flows that 
	 enter the domain than for flows that leave the domain.
	 Whether the flow enters of leaves the domain is decided based on
	 the normal component of the velocity at the face midpoint.
	 --> If the flow enters the domain then the reference data is used.
	 --> If the flow leaves the domain then an interior extrapolated value is used. */

      // Compute the normal velocity
      Vn = dot(SolnBlk.VelocityAtLocation(SolnBlk.Grid.xfaceN(i,SolnBlk.JCu)),
	       SolnBlk.Grid.nfaceN(i,SolnBlk.JCu));
	  
      if (Vn <= ZERO){
	// The flow enters the domain
	// Use UoN as reference value

	if (SolnBlk.Grid.BndNorthSpline.getFluxCalcMethod() == ReconstructionBasedFlux) {
	  // Fill in the ghost cells.
	  for( ghost = 1; ghost <= SolnBlk.Nghost; ++ghost){
	    SolnBlk.U[i][SolnBlk.JCu+ghost] = SolnBlk.UoN[i];
	  }
	  
	} else {

	  // Compute the dx in the normal direction between the Gauss point and the first interior cell centroid.
	  dx_normal = ( SolnBlk.Grid.xfaceN(i,SolnBlk.JCu) - 
			SolnBlk.Grid.CellCentroid(i,SolnBlk.JCu) ) * SolnBlk.Grid.nfaceN(i,SolnBlk.JCu);
	  // Estimate the solution gradient in the normal direction
	  dUdn = (SolnBlk.UoN[i] - SolnBlk.U[i][SolnBlk.JCu])/dx_normal;
	      
	  // Fill in the ghost cells.
	  for( ghost = 1; ghost <= SolnBlk.Nghost; ++ghost){
	    dx_normal = ( SolnBlk.Grid.CellCentroid(i,SolnBlk.JCu+ghost) -
			  SolnBlk.Grid.xfaceN(i,SolnBlk.JCu)) * SolnBlk.Grid.nfaceN(i,SolnBlk.JCu);	    
	    SolnBlk.U[i][SolnBlk.JCu+ghost] = ( SolnBlk.UoN[i] + dUdn * dx_normal );
	  }

	}
      } else {
	// The flow leaves the domain
	// Impose constant extrapolation
	for( ghost = 1; ghost <= SolnBlk.Nghost; ++ghost){
	  SolnBlk.U[i][SolnBlk.JCu+ghost] = SolnBlk.U[i][SolnBlk.JCu];
	}
      }
      break;
      
    case BC_SYMMETRY_PLANE :
      throw runtime_error("BCs() ERROR! Symmetry plane BC hasn't been implemented yet!");
      break;
      
    case BC_EXTRAPOLATE :
    case BC_LINEAR_EXTRAPOLATION :
      throw runtime_error("BCs() ERROR! Linear extrapolation BC hasn't been implemented yet!");
      break;
      
    case BC_OUTFLOW :
      // Impose zero derivative in the j-direction at the boundary
      // This is equivalent with a constant extrapolation
    case BC_CONSTANT_EXTRAPOLATION :
      for( ghost = 1; ghost <= SolnBlk.Nghost; ++ghost){
	SolnBlk.U[i][SolnBlk.JCu+ghost] = SolnBlk.U[i][SolnBlk.JCu];
      }
      break;

    case BC_EXACT_SOLUTION :
      // Leave the ghost cell values unchanged.
      break;

    default:		
      // Impose constant extrapolation
      for( ghost = 1; ghost <= SolnBlk.Nghost; ++ghost){
	SolnBlk.U[i][SolnBlk.JCu+ghost] = SolnBlk.U[i][SolnBlk.JCu];
      }
      break;
    } /* endswitch */
  } /* endfor */


  // Impose high-order boundary conditions
  SolnBlk.BCs_HighOrder();

}

/******************************************************//**
 * Routine: CFL                                         
 *                                                      
 * Determines the allowable global and local time steps 
 * (for explicit Euler time stepping scheme) for the    
 * specified quadrilateral solution block according to  
 * the Courant-Friedrichs-Lewy condition for the        
 * advection terms and a semi-impirical criteria for    
 * the diffusion and source terms.                      
 *                                                      
 ********************************************************/
double CFL(AdvectDiffuse2D_Quad_Block &SolnBlk,
           AdvectDiffuse2D_Input_Parameters &Input_Parameters) {

  int i, j;
  double dtMin, d_i, d_j, v_i, v_j, a, dt_cfl, dt_diff, dt_src;

  dtMin = MILLION;

  for ( j  = SolnBlk.JCl ; j <= SolnBlk.JCu ; ++j ) {
    for ( i = SolnBlk.ICl ; i <= SolnBlk.ICu ; ++i ) {
      /* Use different stability criteria to determine
	 the local time step for all interior cells */

      d_i = TWO*(SolnBlk.Grid.Cell[i][j].A/
		 (SolnBlk.Grid.lfaceE(i, j)+SolnBlk.Grid.lfaceW(i, j)));
      d_j = TWO*(SolnBlk.Grid.Cell[i][j].A/
		 (SolnBlk.Grid.lfaceN(i, j)+SolnBlk.Grid.lfaceS(i, j)));

      // Determine stability limit imposed by the advection term
      v_i = HALF*(SolnBlk.VelocityAtCellCentroid(i,j)*
		  (SolnBlk.Grid.nfaceE(i, j)-SolnBlk.Grid.nfaceW(i, j)));
      v_j = HALF*(SolnBlk.VelocityAtCellCentroid(i,j)*
		  (SolnBlk.Grid.nfaceN(i, j)-SolnBlk.Grid.nfaceS(i, j)));

      if (fabs(v_i) > TOLER) {
	dt_cfl = d_i/fabs(v_i);
      } else {
	dt_cfl = MILLION;
      } /* endif */
      if (fabs(v_j) > TOLER) {
	dt_cfl = min(dt_cfl, d_j/fabs(v_j));
      } /* endif */


	// Determine stability limit imposed by the diffusion term
      if (SolnBlk.DiffusionCoeffAtCellCentroid(i,j) > TOLER) {
	dt_diff = HALF*min(sqr(d_i), sqr(d_j))/fabs(SolnBlk.DiffusionCoeffAtCellCentroid(i,j));
      } else {
	dt_diff = MILLION;
      } /* endif */

	// Determine stability limit imposed by the source term
      dt_src = SolnBlk.SourceTermStabilityLimit(i,j);

      // Determine the locally allowed time step
      SolnBlk.dt[i][j] = min(min(dt_cfl, dt_diff), dt_src);

      // Update the global time step
      dtMin = min(dtMin, SolnBlk.dt[i][j]);
    } /* endfor */
  } /* endfor */

  /* Return the global time step. */

  return (dtMin);

}

/******************************************************//**
 * Routine: Set_Global_TimeStep                         
 *                                                      
 * Assigns global time step to specified solution block 
 * for time-accurate calculations.                      
 *                                                      
 ********************************************************/
void Set_Global_TimeStep(AdvectDiffuse2D_Quad_Block &SolnBlk,
                         const double &Dt_min) {

  int i, j;

  for ( j  = SolnBlk.JCl ; j <= SolnBlk.JCu ; ++j ) {
    for ( i = SolnBlk.ICl ; i <= SolnBlk.ICu ; ++i ) {
      SolnBlk.dt[i][j] = Dt_min;
    } /* endfor */
  } /* endfor */

}

/******************************************************//**
 * Routine: L1_Norm_Residual                            
 *                                                      
 * Determines the L1-norm of the solution residual for  
 * the specified quadrilateral solution block.          
 * Useful for monitoring convergence of the solution    
 * for steady state problems.                           
 *                                                      
 ********************************************************/
double L1_Norm_Residual(AdvectDiffuse2D_Quad_Block &SolnBlk) {

  int i, j;
  double l1_norm;

  l1_norm = ZERO;

  for ( j  = SolnBlk.JCl ; j <= SolnBlk.JCu ; ++j ) {
    for ( i = SolnBlk.ICl ; i <= SolnBlk.ICu ; ++i ) {
      l1_norm += fabs(SolnBlk.dUdt[i][j][0][SolnBlk.residual_variable]);
    } /* endfor */
  } /* endfor */

  return (l1_norm);

}

/******************************************************//**
 * Routine: L2_Norm_Residual                            
 *                                                      
 * Determines the L2-norm of the solution residual for  
 * the specified quadrilateral solution block.          
 * Useful for monitoring convergence of the solution    
 * for steady state problems.                           
 *                                                      
 ********************************************************/
double L2_Norm_Residual(AdvectDiffuse2D_Quad_Block &SolnBlk) {

  int i, j;
  double l2_norm;

  l2_norm = ZERO;

  for ( j  = SolnBlk.JCl ; j <= SolnBlk.JCu ; ++j ) {
    for ( i = SolnBlk.ICl ; i <= SolnBlk.ICu ; ++i ) {
      l2_norm += sqr(SolnBlk.dUdt[i][j][0][SolnBlk.residual_variable]);
    } /* endfor */
  } /* endfor */

  l2_norm = sqrt(l2_norm);

  return (l2_norm);

}

/******************************************************//**
 * Routine: Max_Norm_Residual                           
 *                                                      
 * Determines the maximum norm of the solution residual 
 * for the specified quadrilateral solution block.      
 * Useful for monitoring convergence of the solution    
 * for steady state problems.                           
 *                                                      
 ********************************************************/
double Max_Norm_Residual(AdvectDiffuse2D_Quad_Block &SolnBlk) {

  int i, j;
  double max_norm;

  max_norm = ZERO;

  for ( j  = SolnBlk.JCl ; j <= SolnBlk.JCu ; ++j ) {
    for ( i = SolnBlk.ICl ; i <= SolnBlk.ICu ; ++i ) {
      max_norm = max(max_norm, fabs(SolnBlk.dUdt[i][j][0][SolnBlk.residual_variable]));
    } /* endfor */
  } /* endfor */

  return (max_norm);

}

/******************************************************//**
 * Routine: Linear_Reconstruction_GreenGauss            
 *                                                      
 * Performs the reconstruction of a limited piecewise   
 * linear solution state within a given cell (i,j) of   
 * the computational mesh for the specified             
 * quadrilateral solution block.  A Green-Gauss         
 * approach is used in the evaluation of the unlimited  
 * solution gradients.  Several slope limiters may be   
 * used.                                                
 *                                                      
 ********************************************************/
void Linear_Reconstruction_GreenGauss(AdvectDiffuse2D_Quad_Block &SolnBlk,
				      const int i, 
                                      const int j,
                                      const int Limiter) {

  int n, n_pts, i_index[8], j_index[8];
  double u0Min, u0Max, uQuad[4], phi;
  double DxDx_ave, DxDy_ave, DyDy_ave;
  double l_north, l_south, l_east, l_west;
  Vector2D n_north, n_south, n_east, n_west, dX;
  double u_nw, u_ne, u_sw, u_se, u_face, 
    Du, DuDx_ave, DuDy_ave;

  /* Carry out the limited solution reconstruction in
     each cell of the computational mesh. */

  // Determine the number of neighbouring cells to
  // be used in the reconstruction procedure.  Away from
  // boundaries this 8 neighbours will be used.
  if (i <= SolnBlk.ICl-2 || i >= SolnBlk.ICu+2 || j <= SolnBlk.JCl-2 || j >= SolnBlk.JCu+2) {
    n_pts = 0;
  } else if ((i == SolnBlk.ICl-1) && (SolnBlk.Grid.BCtypeW[j] != BC_NONE)) {
    if (j == SolnBlk.JCl-1 || j == SolnBlk.JCu+1) {
      n_pts = 0;
    } else if (SolnBlk.Grid.BCtypeW[j] == BC_PERIODIC ||
	       SolnBlk.Grid.BCtypeW[j] == BC_NEUMANN ||
	       SolnBlk.Grid.BCtypeW[j] == BC_ROBIN) {
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
      } /* endif */
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
      } /* endif */
    } /* endif */           
  } else if ((i == SolnBlk.ICu+1) && (SolnBlk.Grid.BCtypeE[j] != BC_NONE)) {
    if (j == SolnBlk.JCl-1 || j == SolnBlk.JCu+1) {
      n_pts = 0;
    } else if (SolnBlk.Grid.BCtypeE[j] == BC_PERIODIC ||
	       SolnBlk.Grid.BCtypeE[j] == BC_NEUMANN ||
	       SolnBlk.Grid.BCtypeE[j] == BC_ROBIN) {
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
      } /* endif */
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
      } /* endif */
    } /* endif */
  } else if ((j == SolnBlk.JCl-1) && (SolnBlk.Grid.BCtypeS[i] != BC_NONE)) {
    if (i == SolnBlk.ICl-1 || i == SolnBlk.ICu+1) {
      n_pts = 0;
    } else if (SolnBlk.Grid.BCtypeS[i] == BC_PERIODIC ||
	       SolnBlk.Grid.BCtypeS[i] == BC_NEUMANN ||
	       SolnBlk.Grid.BCtypeS[i] == BC_ROBIN) {
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
      } /* endif */
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
      } /* endif */
    } /* endif */
  } else if ((j == SolnBlk.JCu+1) && (SolnBlk.Grid.BCtypeN[i] != BC_NONE)) {
    if (i == SolnBlk.ICl-1 || i == SolnBlk.ICu+1) {
      n_pts = 0;
    } else if (SolnBlk.Grid.BCtypeN[i] == BC_PERIODIC ||
	       SolnBlk.Grid.BCtypeN[i] == BC_NEUMANN ||
	       SolnBlk.Grid.BCtypeN[i] == BC_ROBIN) {
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
      } /* endif */
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
      } /* endif */
    } /* endif */
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
  } /* endif */
    
    // Perform reconstruction.  
  if (n_pts > 0) {
    // If 8 neighbours are used, apply Green-Gauss reconstruction
    if (n_pts == 8) {
      u_nw = SolnBlk.unNW(i, j);
      u_ne = SolnBlk.unNE(i, j);
      u_sw = SolnBlk.unSW(i, j);
      u_se = SolnBlk.unSE(i, j);

      l_north = SolnBlk.Grid.lfaceN(i, j);
      l_south = SolnBlk.Grid.lfaceS(i, j);
      l_east = SolnBlk.Grid.lfaceE(i, j);
      l_west = SolnBlk.Grid.lfaceW(i, j);

      n_north = SolnBlk.Grid.nfaceN(i, j);
      n_south = SolnBlk.Grid.nfaceS(i, j);
      n_east = SolnBlk.Grid.nfaceE(i, j);
      n_west = SolnBlk.Grid.nfaceW(i, j);

      u_face = HALF*(u_nw+u_ne)*l_north; 
      SolnBlk.dUdx[i][j] = AdvectDiffuse2D_State(u_face*n_north.x);
      SolnBlk.dUdy[i][j] = AdvectDiffuse2D_State(u_face*n_north.y);

      u_face = HALF*(u_sw+u_se)*l_south; 
      SolnBlk.dUdx[i][j] += AdvectDiffuse2D_State(u_face*n_south.x);
      SolnBlk.dUdy[i][j] += AdvectDiffuse2D_State(u_face*n_south.y);

      u_face = HALF*(u_ne+u_se)*l_east; 
      SolnBlk.dUdx[i][j] += AdvectDiffuse2D_State(u_face*n_east.x);
      SolnBlk.dUdy[i][j] += AdvectDiffuse2D_State(u_face*n_east.y);

      u_face = HALF*(u_nw+u_sw)*l_west; 
      SolnBlk.dUdx[i][j] += AdvectDiffuse2D_State(u_face*n_west.x);
      SolnBlk.dUdy[i][j] += AdvectDiffuse2D_State(u_face*n_west.y);

      SolnBlk.dUdx[i][j] = SolnBlk.dUdx[i][j]/
	SolnBlk.Grid.Cell[i][j].A;
      SolnBlk.dUdy[i][j] = SolnBlk.dUdy[i][j]/
	SolnBlk.Grid.Cell[i][j].A;

      // If <8 neighbours are used, apply least-squares reconstruction
    } else {
      DuDx_ave = ZERO;
      DuDy_ave = ZERO;
      DxDx_ave = ZERO;
      DxDy_ave = ZERO;
      DyDy_ave = ZERO;
    
      for ( n = 0 ; n <= n_pts-1 ; ++n ) {
	dX = SolnBlk.Grid.Cell[ i_index[n] ][ j_index[n] ].Xc - SolnBlk.Grid.Cell[i][j].Xc;
	Du = SolnBlk.U[ i_index[n] ][ j_index[n] ][1] - SolnBlk.U[i][j][1];
	DuDx_ave += Du*dX.x;
	DuDy_ave += Du*dX.y;
	DxDx_ave += dX.x*dX.x;
	DxDy_ave += dX.x*dX.y;
	DyDy_ave += dX.y*dX.y;
      } /* endfor */
    					    
      DuDx_ave = DuDx_ave/double(n_pts);
      DuDy_ave = DuDy_ave/double(n_pts);
      DxDx_ave = DxDx_ave/double(n_pts);
      DxDy_ave = DxDy_ave/double(n_pts);
      DyDy_ave = DyDy_ave/double(n_pts);
      SolnBlk.dUdx[i][j][1] = ( (DuDx_ave*DyDy_ave-DuDy_ave*DxDy_ave)/
				(DxDx_ave*DyDy_ave-DxDy_ave*DxDy_ave) );
      SolnBlk.dUdy[i][j][1] = ( (DuDy_ave*DxDx_ave-DuDx_ave*DxDy_ave)/
				(DxDx_ave*DyDy_ave-DxDy_ave*DxDy_ave) );
    } /* endif */
    
    // Calculate slope limiter.
    if (!SolnBlk.Freeze_Limiter) {
      u0Min = SolnBlk.U[i][j][1];
      u0Max = u0Min;
      for ( n = 0 ; n <= n_pts-1 ; ++n ) {
	u0Min = min(u0Min, SolnBlk.U[ i_index[n] ][ j_index[n] ][1]);
	u0Max = max(u0Max, SolnBlk.U[ i_index[n] ][ j_index[n] ][1]);
      } /* endfor */
    
      dX = SolnBlk.Grid.xfaceE(i, j)-SolnBlk.Grid.Cell[i][j].Xc;
      uQuad[0] = SolnBlk.U[i][j][1] + SolnBlk.dUdx[i][j][1]*dX.x + SolnBlk.dUdy[i][j][1]*dX.y;
      dX = SolnBlk.Grid.xfaceW(i, j)-SolnBlk.Grid.Cell[i][j].Xc;
      uQuad[1] = SolnBlk.U[i][j][1] + SolnBlk.dUdx[i][j][1]*dX.x + SolnBlk.dUdy[i][j][1]*dX.y;
      dX = SolnBlk.Grid.xfaceN(i, j)-SolnBlk.Grid.Cell[i][j].Xc;
      uQuad[2] = SolnBlk.U[i][j][1] + SolnBlk.dUdx[i][j][1]*dX.x + SolnBlk.dUdy[i][j][1]*dX.y;
      dX = SolnBlk.Grid.xfaceS(i, j)-SolnBlk.Grid.Cell[i][j].Xc;
      uQuad[3] = SolnBlk.U[i][j][1] + SolnBlk.dUdx[i][j][1]*dX.x + SolnBlk.dUdy[i][j][1]*dX.y;
      
      switch(Limiter) {
      case LIMITER_ONE :
	phi = ONE;
	break;
      case LIMITER_ZERO :
	phi = ZERO;
	break;
      case LIMITER_BARTH_JESPERSEN :
	phi = Limiter_BarthJespersen(uQuad, SolnBlk.U[i][j][1], 
				     u0Min, u0Max, 4);
	break;
      case LIMITER_VENKATAKRISHNAN :
	phi = Limiter_Venkatakrishnan(uQuad, SolnBlk.U[i][j][1], 
				      u0Min, u0Max, 4);
	break;
      case LIMITER_VANLEER :
	phi = Limiter_VanLeer(uQuad, SolnBlk.U[i][j][1], 
			      u0Min, u0Max, 4);
	break;
      case LIMITER_VANALBADA :
	phi = Limiter_VanAlbada(uQuad, SolnBlk.U[i][j][1], 
				u0Min, u0Max, 4);
	break;
      default:
	phi = Limiter_BarthJespersen(uQuad, SolnBlk.U[i][j][1], 
				     u0Min, u0Max, 4);
	break;
      } /* endswitch */
    
      SolnBlk.phi[i][j][1] = phi;
    } /* endif */
  } else {
    SolnBlk.dUdx[i][j].Vacuum();
    SolnBlk.dUdy[i][j].Vacuum(); 
    SolnBlk.phi[i][j].Vacuum();
  } /* endif */
    
}

/******************************************************//**
 * Routine: Linear_Reconstruction_GreenGauss            
 *                                                      
 * Performs the reconstruction of a limited piecewise   
 * linear solution state within each cell of the        
 * computational mesh for the specified quadrilateral   
 * solution block.  A Green-Gauss approach is used      
 * in the evaluation of the unlimited solution          
 * gradients.  Several slope limiters may be used.      
 *                                                      
 ********************************************************/
void Linear_Reconstruction_GreenGauss(AdvectDiffuse2D_Quad_Block &SolnBlk,
				      const int Limiter) {

  int i, j;

  /* Carry out the limited solution reconstruction in
     each cell of the computational mesh. */

  for ( j  = SolnBlk.JCl-SolnBlk.Nghost+1 ; j <= SolnBlk.JCu+SolnBlk.Nghost-1 ; ++j ) {
    for ( i = SolnBlk.ICl-SolnBlk.Nghost+1 ; i <= SolnBlk.ICu+SolnBlk.Nghost-1 ; ++i ) {
      Linear_Reconstruction_GreenGauss(SolnBlk, i, j, Limiter);
    } /* endfor */
  } /* endfor */

}

/******************************************************//**
 * Routine: Linear_Reconstruction_LeastSquares          
 *                                                      
 * Performs the reconstruction of a limited piecewise   
 * linear solution state within a given cell (i,j) of   
 * the computational mesh for the specified             
 * quadrilateral solution block.  A least squares       
 * approach is used in the evaluation of the unlimited  
 * solution gradients.  Several slope limiters may be   
 * used.                                                
 *                                                      
 ********************************************************/
void Linear_Reconstruction_LeastSquares(AdvectDiffuse2D_Quad_Block &SolnBlk,
				        const int i, 
                                        const int j,
                                        const int Limiter) {

  int n, n_pts, i_index[8], j_index[8];
  double u0Min, u0Max, uQuad[4], phi;
  double DxDx_ave, DxDy_ave, DyDy_ave;
  Vector2D dX;
  double Du, DuDx_ave, DuDy_ave;

  /* Carry out the limited solution reconstruction in
     each cell of the computational mesh. */

  // Determine the number of neighbouring cells to
  // be used in the reconstruction procedure.  Away from
  // boundaries this 8 neighbours will be used.
  SolnBlk.SetPiecewiseLinearReconstructionStencil(i,j,
						  i_index,j_index,
						  n_pts);

  // Perform reconstruction.
  if (n_pts > 0) {
    DuDx_ave = ZERO;
    DuDy_ave = ZERO;
    DxDx_ave = ZERO;
    DxDy_ave = ZERO;
    DyDy_ave = ZERO;
    
    for ( n = 0 ; n <= n_pts-1 ; ++n ) {
      dX = SolnBlk.Grid.Cell[ i_index[n] ][ j_index[n] ].Xc - SolnBlk.Grid.Cell[i][j].Xc;
      Du = SolnBlk.U[ i_index[n] ][ j_index[n] ][1] - SolnBlk.U[i][j][1];
      DuDx_ave += Du*dX.x;
      DuDy_ave += Du*dX.y;
      DxDx_ave += dX.x*dX.x;
      DxDy_ave += dX.x*dX.y;
      DyDy_ave += dX.y*dX.y;
    } /* endfor */
    					    
    DuDx_ave = DuDx_ave/double(n_pts);
    DuDy_ave = DuDy_ave/double(n_pts);
    DxDx_ave = DxDx_ave/double(n_pts);
    DxDy_ave = DxDy_ave/double(n_pts);
    DyDy_ave = DyDy_ave/double(n_pts);
    SolnBlk.dUdx[i][j][1] = ( (DuDx_ave*DyDy_ave-DuDy_ave*DxDy_ave)/
			      (DxDx_ave*DyDy_ave-DxDy_ave*DxDy_ave) );
    SolnBlk.dUdy[i][j][1] = ( (DuDy_ave*DxDx_ave-DuDx_ave*DxDy_ave)/
			      (DxDx_ave*DyDy_ave-DxDy_ave*DxDy_ave) );
    
    // Calculate slope limiter.
    if (!SolnBlk.Freeze_Limiter) {
      u0Min = SolnBlk.U[i][j][1];
      u0Max = u0Min;
      for ( n = 0 ; n <= n_pts-1 ; ++n ) {
	u0Min = min(u0Min, SolnBlk.U[ i_index[n] ][ j_index[n] ][1]);
	u0Max = max(u0Max, SolnBlk.U[ i_index[n] ][ j_index[n] ][1]);
      } /* endfor */
    
      dX = SolnBlk.Grid.xfaceE(i, j)-SolnBlk.Grid.Cell[i][j].Xc;
      uQuad[0] = SolnBlk.U[i][j][1] + SolnBlk.dUdx[i][j][1]*dX.x + SolnBlk.dUdy[i][j][1]*dX.y ;
      dX = SolnBlk.Grid.xfaceW(i, j)-SolnBlk.Grid.Cell[i][j].Xc;
      uQuad[1] = SolnBlk.U[i][j][1] + SolnBlk.dUdx[i][j][1]*dX.x + SolnBlk.dUdy[i][j][1]*dX.y ;
      dX = SolnBlk.Grid.xfaceN(i, j)-SolnBlk.Grid.Cell[i][j].Xc;
      uQuad[2] = SolnBlk.U[i][j][1] + SolnBlk.dUdx[i][j][1]*dX.x + SolnBlk.dUdy[i][j][1]*dX.y ;
      dX = SolnBlk.Grid.xfaceS(i, j)-SolnBlk.Grid.Cell[i][j].Xc;
      uQuad[3] = SolnBlk.U[i][j][1] + SolnBlk.dUdx[i][j][1]*dX.x + SolnBlk.dUdy[i][j][1]*dX.y ;
    
      switch(Limiter) {
      case LIMITER_ONE :
	phi = ONE;
	break;
      case LIMITER_ZERO :
	phi = ZERO;
	break;
      case LIMITER_BARTH_JESPERSEN :
	phi = Limiter_BarthJespersen(uQuad, SolnBlk.U[i][j][1], 
				     u0Min, u0Max, 4);
	break;
      case LIMITER_VENKATAKRISHNAN :
	phi = Limiter_Venkatakrishnan(uQuad, SolnBlk.U[i][j][1], 
				      u0Min, u0Max, 4);
	break;
      case LIMITER_VANLEER :
	phi = Limiter_VanLeer(uQuad, SolnBlk.U[i][j][1], 
			      u0Min, u0Max, 4);
	break;
      case LIMITER_VANALBADA :
	phi = Limiter_VanAlbada(uQuad, SolnBlk.U[i][j][1], 
				u0Min, u0Max, 4);
	break;
      default:
	phi = Limiter_BarthJespersen(uQuad, SolnBlk.U[i][j][1], 
				     u0Min, u0Max, 4);
	break;
      } /* endswitch */
    
      SolnBlk.phi[i][j] = AdvectDiffuse2D_State(phi);
    } /* endif */
  } else {
    SolnBlk.dUdx[i][j].Vacuum();
    SolnBlk.dUdy[i][j].Vacuum();
    SolnBlk.phi[i][j].Vacuum();
  } /* endif */
    
}

/******************************************************//**
 * Routine: Linear_Reconstruction_LeastSquares          
 *                                                      
 * Performs the reconstruction of a limited piecewise   
 * linear solution state within each cell of the        
 * computational mesh of the specified quadrilateral    
 * solution block.  A least squares approach is         
 * used in the evaluation of the unlimited solution     
 * gradients.  Several slope limiters may be used.      
 *                                                      
 ********************************************************/
void Linear_Reconstruction_LeastSquares(AdvectDiffuse2D_Quad_Block &SolnBlk,
				        const int Limiter) {

  int i, j;

  /* Carry out the limited solution reconstruction in
     each cell of the computational mesh. */

  for ( j  = SolnBlk.JCl-SolnBlk.Nghost+1 ; j <= SolnBlk.JCu+SolnBlk.Nghost-1 ; ++j ) {
    for ( i = SolnBlk.ICl-SolnBlk.Nghost+1 ; i <= SolnBlk.ICu+SolnBlk.Nghost-1 ; ++i ) {
      Linear_Reconstruction_LeastSquares(SolnBlk, i, j, Limiter);
    } /* endfor */
  } /* endfor */

}

/******************************************************//**
 * Routine: Linear_Reconstruction
 *                                                      
 * Performs the reconstruction of a limited piecewise   
 * linear solution state within each cell of the        
 * computational mesh of the specified quadrilateral    
 * solution block.
 ********************************************************/
void Linear_Reconstruction(AdvectDiffuse2D_Quad_Block &SolnBlk,
			   const int &Reconstruction_Type,
			   const int &Limiter) {

  switch(Reconstruction_Type) {
  case RECONSTRUCTION_GREEN_GAUSS :
    Linear_Reconstruction_GreenGauss(SolnBlk,
				     Limiter);    
    break;
  case RECONSTRUCTION_LEAST_SQUARES :
    Linear_Reconstruction_LeastSquares(SolnBlk,
				       Limiter);
    break;
  default:
    Linear_Reconstruction_LeastSquares(SolnBlk,
				       Limiter);
    break;
  } /* endswitch */
}


/******************************************************//**
 * Routine: Residual_Smoothing                          
 *                                                      
 * Applies implicit residual smoothing to solution      
 * residual.  Note that only residuals of interior cells
 * are smoothed and residuals for cells adjacent to     
 * boundaries are not smoothed.                         
 *                                                      
 ********************************************************/
void Residual_Smoothing(AdvectDiffuse2D_Quad_Block &SolnBlk,
                        const int k_residual,
			double &epsilon,
                        const int number_of_Gauss_Seidel_iterations) {

  int i, j, n;

  for ( j  = SolnBlk.JCl+1; j <= SolnBlk.JCu-1; ++j ) {
    for ( i = SolnBlk.ICl+1; i <= SolnBlk.ICu-1; ++i ) {
      SolnBlk.dUdt[i][j][k_residual+1] =  SolnBlk.dUdt[i][j][k_residual];
    } /* endfor */
  } /* endfor */

  for ( n  = 1; n <= number_of_Gauss_Seidel_iterations; ++n ) {

    for ( j  = SolnBlk.JCl+1; j <= SolnBlk.JCu-1; ++j ) {
      for ( i = SolnBlk.ICl+1; i <= SolnBlk.ICu-1; ++i ) {
	SolnBlk.dUdt[i][j][k_residual+1] = (SolnBlk.dUdt[i][j][k_residual]
					    + epsilon*(SolnBlk.dUdt[i  ][j-1][k_residual+1] +
						       SolnBlk.dUdt[i-1][j  ][k_residual+1] +
						       SolnBlk.dUdt[i+1][j  ][k_residual+1] +
						       SolnBlk.dUdt[i  ][j+1][k_residual+1]))/(ONE + FOUR*epsilon);
      } /* endfor */
    } /* endfor */

  } /* endfor */

  for ( j  = SolnBlk.JCl+1; j <= SolnBlk.JCu-1; ++j ) {
    for ( i = SolnBlk.ICl+1; i <= SolnBlk.ICu-1; ++i ) {
      SolnBlk.dUdt[i][j][k_residual] =  SolnBlk.dUdt[i][j][k_residual+1];
    } /* endfor */
  } /* endfor */

}

/******************************************************//**
 * Routine: Calculate_Refinement_Criteria               
 *                                                      
 * Calculate refinement criteria for the solution       
 * block.                                               
 *                                                      
 * \todo Refactor this function!
 ********************************************************/
void Calculate_Refinement_Criteria(double *refinement_criteria,
				   AdvectDiffuse2D_Input_Parameters &IP,
                                   int &number_refinement_criteria,
                                   AdvectDiffuse2D_Quad_Block &SolnBlk) {

  // Calculate refinement criteria based on smoothness indicator
  if (CENO_Execution_Mode::USE_CENO_ALGORITHM && 
      CENO_Execution_Mode::USE_SMOOTHNESS_INDICATOR_FOR_AMR_CRITERIA) {
    return SolnBlk.Calculate_Refinement_Criteria_HighOrder(refinement_criteria,
							   IP,
							   number_refinement_criteria);
  }

  int i, j;

  double grad_u_x, grad_u_y, grad_u_abs, grad_u_criteria, grad_u_criteria_max;

  /* Set the number of refinement criteria to be used (1):
     (1) refinement criteria #1 based on the gradient of the solution. */

  number_refinement_criteria = 1;

  /* Allocate memory for the refinement criteria */
  SolnBlk.Refinement_Criteria().reserve(number_refinement_criteria);

  /* Initialize the refinement criteria for the block. */

  grad_u_criteria_max = ZERO;

  /* Calculate the refinement criteria for each cell of the 
     computational mesh and assign the maximum value for
     all cells as the refinement criteria for the solution 
     block. */

  for ( j  = SolnBlk.JCl-1 ; j <= SolnBlk.JCu+1 ; ++j ) {
    for ( i = SolnBlk.ICl-1 ; i <= SolnBlk.ICu+1 ; ++i ) {
      // Reconstruct the solution within the cell.
      Linear_Reconstruction_LeastSquares(SolnBlk, i, j, LIMITER_UNLIMITED);

      if (SolnBlk.Grid.Cell[i][j].A > ZERO) {
	// Evaluate refinement criteria #1 based on the gradient
	// of the solution.
	grad_u_x = SolnBlk.dUdx[i][j][1];
	grad_u_y = SolnBlk.dUdy[i][j][1];
	grad_u_abs = sqrt(sqr(grad_u_x) + sqr(grad_u_y));
	// Calculate the grad_u_criteria as DeltaU/(1+fabs(U))
	grad_u_criteria = sqrt(SolnBlk.Grid.Cell[i][j].A)*grad_u_abs/(1.0 + fabs(SolnBlk.U[i][j].u));
	grad_u_criteria_max = max(grad_u_criteria_max, grad_u_criteria);
      } /* endif */

    } /* endfor */
  } /* endfor */

  /* Return the refinement criteria. */
  refinement_criteria[0] = grad_u_criteria_max;

  /* Store the refinement_criteria values in the solution block designated variable */
  SolnBlk.Refinement_Criterion(0) = refinement_criteria[0];

}

/******************************************************//**
 * Routine: Fix_Refined_Block_Boundaries                
 *                                                      
 * Adjusts the locations of the boundary nodes of a     
 * solution block so that the new node locations        
 * match with cell volumes of adjacent solution blocks  
 * that have lower levels of mesh refinement (i.e., are 
 * coarser solution blocks).                            
 *                                                      
 ********************************************************/
void Fix_Refined_Block_Boundaries(AdvectDiffuse2D_Quad_Block &SolnBlk,
                                  const int Fix_North_Boundary,
                                  const int Fix_South_Boundary,
                                  const int Fix_East_Boundary,
                                  const int Fix_West_Boundary) {

  int i, j;
  double ds_ratio, dl, dr;
 
  /* Adjust the node locations at the north boundary. */

  if (Fix_North_Boundary) {
    for ( i = SolnBlk.Grid.INl+1 ; i <= SolnBlk.Grid.INu-1 ; i+=2 ) {
      dl = abs(SolnBlk.Grid.Node[i  ][SolnBlk.Grid.JNu].X - 
	       SolnBlk.Grid.Node[i-1][SolnBlk.Grid.JNu].X);
      dr = abs(SolnBlk.Grid.Node[i+1][SolnBlk.Grid.JNu].X - 
	       SolnBlk.Grid.Node[i  ][SolnBlk.Grid.JNu].X);
      ds_ratio = dl/(dl+dr);
      SolnBlk.Grid.Node[i][SolnBlk.Grid.JNu].X = 
	SolnBlk.Grid.Node[i-1][SolnBlk.Grid.JNu].X +
	ds_ratio*(SolnBlk.Grid.Node[i+1][SolnBlk.Grid.JNu].X-
		  SolnBlk.Grid.Node[i-1][SolnBlk.Grid.JNu].X);
    } /* endfor */
    for ( i = SolnBlk.ICl ; i <= SolnBlk.ICu ; ++i ) {
      SolnBlk.U[i][SolnBlk.JCu] = (SolnBlk.Grid.Cell[i][SolnBlk.JCu].A/
				   SolnBlk.Grid.area(i, SolnBlk.JCu))*SolnBlk.U[i][SolnBlk.JCu];
    } /* endfor */
  } /* endif */

    /* Adjust the node locations at the south boundary. */

  if (Fix_South_Boundary) {
    for ( i = SolnBlk.Grid.INl+1 ; i <= SolnBlk.Grid.INu-1 ; i+=2 ) {
      dl = abs(SolnBlk.Grid.Node[i  ][SolnBlk.Grid.JNl].X - 
	       SolnBlk.Grid.Node[i-1][SolnBlk.Grid.JNl].X);
      dr = abs(SolnBlk.Grid.Node[i+1][SolnBlk.Grid.JNl].X - 
	       SolnBlk.Grid.Node[i  ][SolnBlk.Grid.JNl].X);
      ds_ratio = dl/(dl+dr);
      SolnBlk.Grid.Node[i][SolnBlk.Grid.JNl].X = 
	SolnBlk.Grid.Node[i-1][SolnBlk.Grid.JNl].X +
	ds_ratio*(SolnBlk.Grid.Node[i+1][SolnBlk.Grid.JNl].X-
		  SolnBlk.Grid.Node[i-1][SolnBlk.Grid.JNl].X);
    } /* endfor */
    for ( i = SolnBlk.ICl ; i <= SolnBlk.ICu ; ++i ) {
      SolnBlk.U[i][SolnBlk.JCl] = (SolnBlk.Grid.Cell[i][SolnBlk.JCl].A/
				   SolnBlk.Grid.area(i, SolnBlk.JCl))*SolnBlk.U[i][SolnBlk.JCl];
    } /* endfor */
  } /* endif */

    /* Adjust the node locations at the east boundary. */

  if (Fix_East_Boundary) {
    for ( j  = SolnBlk.Grid.JNl+1; j <= SolnBlk.Grid.JNu-1; j+=2 ) {
      dl = abs(SolnBlk.Grid.Node[SolnBlk.Grid.INu][j  ].X - 
	       SolnBlk.Grid.Node[SolnBlk.Grid.INu][j-1].X);
      dr = abs(SolnBlk.Grid.Node[SolnBlk.Grid.INu][j+1].X - 
	       SolnBlk.Grid.Node[SolnBlk.Grid.INu][j  ].X);
      ds_ratio = dl/(dl+dr);
      SolnBlk.Grid.Node[SolnBlk.Grid.INu][j].X = 
	SolnBlk.Grid.Node[SolnBlk.Grid.INu][j-1].X +
	ds_ratio*(SolnBlk.Grid.Node[SolnBlk.Grid.INu][j+1].X-
		  SolnBlk.Grid.Node[SolnBlk.Grid.INu][j-1].X);
    } /* endfor */
    for ( j = SolnBlk.JCl ; j <= SolnBlk.JCu ; ++j ) {
      SolnBlk.U[SolnBlk.ICu][j] = (SolnBlk.Grid.Cell[SolnBlk.ICu][j].A/
				   SolnBlk.Grid.area(SolnBlk.ICu, j))*SolnBlk.U[SolnBlk.ICu][j];
    } /* endfor */
  } /* endif */

    /* Adjust the node locations at the west boundary. */

  if (Fix_West_Boundary) {
    for ( j  = SolnBlk.Grid.JNl+1; j <= SolnBlk.Grid.JNu-1; j+=2 ) {
      dl = abs(SolnBlk.Grid.Node[SolnBlk.Grid.INl][j  ].X - 
	       SolnBlk.Grid.Node[SolnBlk.Grid.INl][j-1].X);
      dr = abs(SolnBlk.Grid.Node[SolnBlk.Grid.INl][j+1].X - 
	       SolnBlk.Grid.Node[SolnBlk.Grid.INl][j  ].X);
      ds_ratio = dl/(dl+dr);
      SolnBlk.Grid.Node[SolnBlk.Grid.INl][j].X = 
	SolnBlk.Grid.Node[SolnBlk.Grid.INl][j-1].X +
	ds_ratio*(SolnBlk.Grid.Node[SolnBlk.Grid.INl][j+1].X-
		  SolnBlk.Grid.Node[SolnBlk.Grid.INl][j-1].X);
    } /* endfor */
    for ( j = SolnBlk.JCl ; j <= SolnBlk.JCu ; ++j ) {
      SolnBlk.U[SolnBlk.ICl][j] = (SolnBlk.Grid.Cell[SolnBlk.ICl][j].A/
				   SolnBlk.Grid.area(SolnBlk.ICl, j))*SolnBlk.U[SolnBlk.ICl][j];
    } /* endfor */
  } /* endif */

  // Update geometric information only if modifications occurred
  if (Fix_North_Boundary || Fix_South_Boundary || 
      Fix_East_Boundary || Fix_West_Boundary ){

    /* Require update of the interior cells geometric properties. */
    SolnBlk.Grid.Schedule_Interior_Mesh_Update();

    /* Reset the boundary condition types at the block boundaries. */
    Set_BCs(SolnBlk.Grid);

    /* Recompute the exterior nodes for the block quadrilateral mesh. */
    Update_Exterior_Nodes(SolnBlk.Grid);

    /* Recompute the cells for the block quadrilateral mesh. */
    Update_Cells(SolnBlk.Grid);
  }

}

/******************************************************//**
 * Routine: Unfix_Refined_Block_Boundaries              
 *                                                      
 * Returns the adjusted the locations of the boundary   
 * nodes of a solution block to their original          
 * unmodified positions.                                
 *                                                      
 ********************************************************/
void Unfix_Refined_Block_Boundaries(AdvectDiffuse2D_Quad_Block &SolnBlk) {

  int i, j;
  double sp_l, sp_r, sp_m, ds_ratio, dl, dr;
  bool ModifiedGrid(false);
 
  /* Return the nodes at the north boundary
     to their original positions. */

  if (SolnBlk.Grid.BndNorthSpline.np != 0) {
    for ( i = SolnBlk.Grid.INl+1 ; i < SolnBlk.Grid.INu ; i += 2 ) {
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
      SolnBlk.Grid.Node[i][SolnBlk.Grid.JNu].X = Spline(sp_m, SolnBlk.Grid.BndNorthSpline);
    } /* endfor */
    for ( i = SolnBlk.ICl ; i <= SolnBlk.ICu ; ++i ) {
      SolnBlk.U[i][SolnBlk.JCu] = (SolnBlk.Grid.Cell[i][SolnBlk.JCu].A/
				   SolnBlk.Grid.area(i, SolnBlk.JCu))*SolnBlk.U[i][SolnBlk.JCu];
    } /* endfor */
    ModifiedGrid = true;
  } /* endif */

    /* Return the nodes at the south boundary
       to their original positions. */

  if (SolnBlk.Grid.BndSouthSpline.np != 0) {
    for ( i = SolnBlk.Grid.INl+1 ; i < SolnBlk.Grid.INu ; i += 2 ) {
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
      SolnBlk.Grid.Node[i][SolnBlk.Grid.JNl].X = Spline(sp_m, SolnBlk.Grid.BndSouthSpline);
    } /* endfor */
    for ( i = SolnBlk.ICl ; i <= SolnBlk.ICu ; ++i ) {
      SolnBlk.U[i][SolnBlk.JCl] = (SolnBlk.Grid.Cell[i][SolnBlk.JCl].A/
				   SolnBlk.Grid.area(i, SolnBlk.JCl))*SolnBlk.U[i][SolnBlk.JCl];
    } /* endfor */
    ModifiedGrid = true;
  } /* endif */

    /* Return the nodes at the east boundary
       to their original positions. */

  if (SolnBlk.Grid.BndEastSpline.np != 0) {
    for (j  = SolnBlk.Grid.JNl+1 ; j < SolnBlk.Grid.JNu ; j += 2 ) {
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
      SolnBlk.Grid.Node[SolnBlk.Grid.INu][j].X = Spline(sp_m, SolnBlk.Grid.BndEastSpline);
    } /* endfor */
    for ( j = SolnBlk.JCl ; j <= SolnBlk.JCu ; ++j ) {
      SolnBlk.U[SolnBlk.ICu][j] = (SolnBlk.Grid.Cell[SolnBlk.ICu][j].A/
				   SolnBlk.Grid.area(SolnBlk.ICu, j))*SolnBlk.U[SolnBlk.ICu][j];
    } /* endfor */
    ModifiedGrid = true;
  } /* endif */

    /* Return the nodes at the west boundary
       to their original positions. */

  if (SolnBlk.Grid.BndWestSpline.np != 0) {
    for (j  = SolnBlk.Grid.JNl+1 ; j < SolnBlk.Grid.JNu ; j += 2 ) {
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
      SolnBlk.Grid.Node[SolnBlk.Grid.INl][j].X = Spline(sp_m, SolnBlk.Grid.BndWestSpline);
    } /* endfor */
    for ( j = SolnBlk.JCl ; j <= SolnBlk.JCu ; ++j ) {
      SolnBlk.U[SolnBlk.ICl][j] = (SolnBlk.Grid.Cell[SolnBlk.ICl][j].A/
				   SolnBlk.Grid.area(SolnBlk.ICl, j))*SolnBlk.U[SolnBlk.ICl][j];
    } /* endfor */
    ModifiedGrid = true;
  } /* endif */

  if (ModifiedGrid){

    /* Require update of the interior cells geometric properties. */
    SolnBlk.Grid.Schedule_Interior_Mesh_Update();

    /* Reset the boundary condition types at the block boundaries. */
    Set_BCs(SolnBlk.Grid);
    
    /* Recompute the exterior nodes for the block quadrilateral mesh. */
    Update_Exterior_Nodes(SolnBlk.Grid);

    /* Recompute the cells for the block quadrilateral mesh. */
    Update_Cells(SolnBlk.Grid);
  }

}

/**************************************************************//**
 * Routine: Apply_Boundary_Flux_Corrections                     
 *                                                              
 * Apply flux corrections at boundaries of the solution         
 * block to ensure that the scheme is conservative at           
 * boundaries with mesh resolution changes.                     
 *                                                              
 ****************************************************************/
void Apply_Boundary_Flux_Corrections(AdvectDiffuse2D_Quad_Block &SolnBlk,
                                     const int Number_Neighbours_North_Boundary,
                                     const int Number_Neighbours_South_Boundary,
                                     const int Number_Neighbours_East_Boundary,
                                     const int Number_Neighbours_West_Boundary) {

  int i, j;
 
  /* Correct the fluxes at the north boundary as required. */

  if (Number_Neighbours_North_Boundary == 2) {
    for ( i = SolnBlk.ICl ; i <= SolnBlk.ICu ; ++i ) {
      SolnBlk.dUdt[i][SolnBlk.JCu][0] -= 
	SolnBlk.FluxN[i]/SolnBlk.Grid.Cell[i][SolnBlk.JCu].A;
    } /* endfor */
  } /* endif */

    /* Correct the fluxes at the south boundary as required. */

  if (Number_Neighbours_South_Boundary == 2) {
    for ( i = SolnBlk.ICl ; i <= SolnBlk.ICu ; ++i ) {
      SolnBlk.dUdt[i][SolnBlk.JCl][0] -= 
	SolnBlk.FluxS[i]/SolnBlk.Grid.Cell[i][SolnBlk.JCl].A;
    } /* endfor */
  } /* endif */

    /* Correct the fluxes at the east boundary as required. */

  if (Number_Neighbours_East_Boundary == 2) {
    for ( j= SolnBlk.JCl ; j <= SolnBlk.JCu ; ++j ) {
      SolnBlk.dUdt[SolnBlk.ICu][j][0] -= 
	SolnBlk.FluxE[j]/SolnBlk.Grid.Cell[SolnBlk.ICu][j].A;
    } /* endfor */
  } /* endif */

    /* Correct the fluxes at the west boundary as required. */

  if (Number_Neighbours_West_Boundary == 2) {
    for ( j= SolnBlk.JCl ; j <= SolnBlk.JCu ; ++j ) {
      SolnBlk.dUdt[SolnBlk.ICl][j][0] -= 
	SolnBlk.FluxW[j]/SolnBlk.Grid.Cell[SolnBlk.ICl][j].A;
    } /* endfor */
  } /* endif */

}

/**************************************************************//**
 * Routine: Apply_Boundary_Flux_Corrections_Multistage_Explicit 
 *                                                              
 * Apply flux corrections at boundaries of the solution         
 * block to ensure that the scheme is conservative at           
 * boundaries with mesh resolution changes.                     
 *                                                              
 ****************************************************************/
void Apply_Boundary_Flux_Corrections_Multistage_Explicit(AdvectDiffuse2D_Quad_Block &SolnBlk,
                                                         const int i_stage,
                                                         const int n_stage,
	                                                 const double &CFL_Number,
                                                         const int Time_Integration_Type,
                                                         const int Reconstruction_Type,
                                                         const int Limiter_Type,
                                                         const int Number_Neighbours_North_Boundary,
                                                         const int Number_Neighbours_South_Boundary,
                                                         const int Number_Neighbours_East_Boundary,
                                                         const int Number_Neighbours_West_Boundary) {

  int i, j, k_residual;
  double omega;

  /* Evaluate the time step fraction and residual storage location for the stage. */
    
  switch(Time_Integration_Type) {
  case TIME_STEPPING_EXPLICIT_EULER :
    omega = Runge_Kutta(i_stage, n_stage);
    k_residual = 0;
    break;
  case TIME_STEPPING_EXPLICIT_PREDICTOR_CORRECTOR :
    omega = Runge_Kutta(i_stage, n_stage);
    k_residual = 0;
    break;
  case TIME_STEPPING_EXPLICIT_RUNGE_KUTTA :
    omega = Runge_Kutta(i_stage, n_stage);
    k_residual = 0;
    if (n_stage == 4) {
      if (i_stage == 4) {
	k_residual = 0;
      } else {
	k_residual = i_stage - 1;
      } /* endif */
    } /* endif */
    break;
  case TIME_STEPPING_MULTISTAGE_OPTIMAL_SMOOTHING :
    omega = MultiStage_Optimally_Smoothing(i_stage, 
					   n_stage,
					   Limiter_Type);
    k_residual = 0;
    break;
  default:
    omega = Runge_Kutta(i_stage, n_stage);
    k_residual = 0;
    break;
  } /* endswitch */

    /* Correct the fluxes at the north boundary as required. */

  if (Number_Neighbours_North_Boundary == 2) {
    for ( i = SolnBlk.ICl ; i <= SolnBlk.ICu ; ++i ) {
	SolnBlk.dUdt[i][SolnBlk.JCu][k_residual] -= ( (CFL_Number*SolnBlk.dt[i][SolnBlk.JCu])*SolnBlk.FluxN[i]/
						      SolnBlk.Grid.Cell[i][SolnBlk.JCu].A );
    } /* endfor */
  } /* endif */

    /* Correct the fluxes at the south boundary as required. */

  if (Number_Neighbours_South_Boundary == 2) {
    for ( i = SolnBlk.ICl ; i <= SolnBlk.ICu ; ++i ) {
	SolnBlk.dUdt[i][SolnBlk.JCl][k_residual] -= ( (CFL_Number*SolnBlk.dt[i][SolnBlk.JCl])*SolnBlk.FluxS[i]/
						      SolnBlk.Grid.Cell[i][SolnBlk.JCl].A );
    } /* endfor */
  } /* endif */

    /* Correct the fluxes at the east boundary as required. */

  if (Number_Neighbours_East_Boundary == 2) {
    for ( j= SolnBlk.JCl ; j <= SolnBlk.JCu ; ++j ) {
	SolnBlk.dUdt[SolnBlk.ICu][j][k_residual] -= ( (CFL_Number*SolnBlk.dt[SolnBlk.ICu][j])*SolnBlk.FluxE[j]/
						      SolnBlk.Grid.Cell[SolnBlk.ICu][j].A );
    } /* endfor */
  } /* endif */

    /* Correct the fluxes at the west boundary as required. */

  if (Number_Neighbours_West_Boundary == 2) {
    for ( j= SolnBlk.JCl ; j <= SolnBlk.JCu ; ++j ) {
	SolnBlk.dUdt[SolnBlk.ICl][j][k_residual] -= ( (CFL_Number*SolnBlk.dt[SolnBlk.ICl][j])*SolnBlk.FluxW[j]/
						      SolnBlk.Grid.Cell[SolnBlk.ICl][j].A );
    } /* endfor */
  } /* endif */

}

/******************************************************//**
 * Routine: dUdt_Residual_Evaluation                    
 *                                                      
 * This routine evaluates the residual for the specified 
 * solution block using a 2nd-order or a high-order (up to 4th-order)
 * finite-volume spatial discretization scheme, whichever
 * is required in the input parameters.
 * The residual is stored in dUdt[][][0].               
 *                                                      
 ********************************************************/
int dUdt_Residual_Evaluation(AdvectDiffuse2D_Quad_Block &SolnBlk,
			     AdvectDiffuse2D_Input_Parameters &IP) {

  if (IP.i_Reconstruction == RECONSTRUCTION_HIGH_ORDER){
    // calculate the high-order residual
    return SolnBlk.dUdt_Residual_Evaluation_HighOrder(IP);
  } else {
    // calculate the 2nd-order residual
    return SolnBlk.dUdt_Residual_Evaluation(IP);
  }
}

/******************************************************//**
 * Routine: dUdt_Multistage_Explicit                    
 *                                                      
 * This routine determines the solution residuals for a 
 * given stage of a variety of multi-stage explicit     
 * time integration schemes for a given solution block. 
 *                                                      
 ********************************************************/
int dUdt_Multistage_Explicit(AdvectDiffuse2D_Quad_Block &SolnBlk,
                             const int i_stage,
                             AdvectDiffuse2D_Input_Parameters &IP) {

  if (IP.i_Reconstruction == RECONSTRUCTION_HIGH_ORDER){
    // calculate the high-order residual
    return SolnBlk.dUdt_Multistage_Explicit_HighOrder(i_stage,IP);
  } else {
    // calculate the 2nd-order residual
    return SolnBlk.dUdt_Multistage_Explicit(i_stage,IP);
  }
}

/******************************************************//**
 * Routine: Update_Solution_Multistage_Explicit         
 *                                                      
 * This routine updates solution states of the given    
 * solution block for a variety of multi-stage explicit 
 * time integration schemes.                            
 *                                                      
 ********************************************************/
int Update_Solution_Multistage_Explicit(AdvectDiffuse2D_Quad_Block &SolnBlk,
                                        const int i_stage,
                                        AdvectDiffuse2D_Input_Parameters &IP) {

  int i, j, k_residual;
  double omega;

  /* Perform update of solution variables for stage 
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
    
  /* Update solution variables for this stage. */
  for ( j  = SolnBlk.JCl ; j <= SolnBlk.JCu ; ++j ) {
    for ( i = SolnBlk.ICl ; i <= SolnBlk.ICu ; ++i ) {

      // Explicit update.
      SolnBlk.U[i][j] = SolnBlk.Uo[i][j] + omega*SolnBlk.dUdt[i][j][k_residual];

    } /* endfor */    
  } /* endfor */

    /* Solution successfully updated. */

  return (0);

}
