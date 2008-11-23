/*!\file Euler2DQuadSingleBlock.cc
  \brief Single-Block Versions of Subroutines for 2D Euler Multi-Block Quadrilateral Mesh Solution Classes. */

/* Include required C++ libraries. */
// None

/* Using std namespace functions */
// None

/* Include CFFC header files */
#include "Euler2DQuad.h" /* Include 2D Euler quadrilateral mesh solution header file. */
#include "Euler2D_ICs.h" /* Include 2D Euler analytically defined initial conditions header file. */


/**************************************************************************
 * Euler2D_Quad_Block -- Single Block External Subroutines.               *
 **************************************************************************/

/********************************************************
 * Routine: Write_Solution_Block                        *
 *                                                      *
 * Writes the cell centred solution values of the       *
 * specified quadrilateral solution block to the        *
 * specified output stream for restart purposes.        *
 *                                                      *
 ********************************************************/
void Write_Solution_Block(Euler2D_Quad_Block &SolnBlk,
	                  ostream &Out_File) {

    Out_File << setprecision(14) << SolnBlk << setprecision(6);

}

/********************************************************
 * Routine: Read_Solution_Block                         *
 *                                                      *
 * Reads the cell centred solution values for the       *
 * specified quadrilateral solution block from the      *
 * specified input stream as required for restart       *
 * purposes.                                            *
 *                                                      *
 ********************************************************/
void Read_Solution_Block(Euler2D_Quad_Block &SolnBlk,
	                 istream &In_File) {

    In_File >> SolnBlk;

}

/********************************************************
 * Routine: Broadcast_Solution_Block                    *
 *                                                      *
 * Broadcast quadrilateral solution block to all        *
 * processors involved in the calculation from the      *
 * primary processor using the MPI broadcast routine.   *
 *                                                      *
 ********************************************************/
void Broadcast_Solution_Block(Euler2D_Quad_Block &SolnBlk) {

#ifdef _MPI_VERSION
    int i, j, ni, nj, ng, nr, block_allocated, buffer_size;
    double *buffer;

    // High-order related variables
    int NumberOfHighOrderVariables;
    vector<int> ReconstructionOrders;

    /* Broadcast the number of cells in each direction. */

    if (CFFC_Primary_MPI_Processor()) {
      ni = SolnBlk.NCi;
      nj = SolnBlk.NCj;
      ng = SolnBlk.Nghost;
      nr = SolnBlk.residual_variable;
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
    MPI::COMM_WORLD.Bcast(&nr,1,MPI::INT,0);
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
       // Set the block static variables if they were not previously assigned.
       if (SolnBlk.residual_variable != nr) SolnBlk.residual_variable = nr;
    } /* endif */

    /* Broadcast the axisymmetric/planar flow indicator. */

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
       buffer = new double[4*ni*nj];

       if (CFFC_Primary_MPI_Processor()) {
          buffer_size = 0;
          for (j  = SolnBlk.JCl-SolnBlk.Nghost ; j <= SolnBlk.JCu+SolnBlk.Nghost ; ++j ) {
              for ( i = SolnBlk.ICl-SolnBlk.Nghost ; i <= SolnBlk.ICu+SolnBlk.Nghost ; ++i ) {
 	          buffer[buffer_size  ] = SolnBlk.U[i][j].d;
 	          buffer[buffer_size+1] = SolnBlk.U[i][j].dv.x;
 	          buffer[buffer_size+2] = SolnBlk.U[i][j].dv.y;
 	          buffer[buffer_size+3] = SolnBlk.U[i][j].E;
                  buffer_size = buffer_size + 4;
              } /* endfor */
          } /* endfor */
       } /* endif */

       buffer_size = 4*ni*nj;
       MPI::COMM_WORLD.Bcast(buffer, buffer_size, MPI::DOUBLE, 0);

       if (!CFFC_Primary_MPI_Processor()) {
          buffer_size = 0;
          for (j  = SolnBlk.JCl-SolnBlk.Nghost ; j <= SolnBlk.JCu+SolnBlk.Nghost ; ++j ) {
              for ( i = SolnBlk.ICl-SolnBlk.Nghost ; i <= SolnBlk.ICu+SolnBlk.Nghost ; ++i ) {
 	          SolnBlk.U[i][j].d    = buffer[buffer_size];
 	          SolnBlk.U[i][j].dv.x = buffer[buffer_size+1];
 	          SolnBlk.U[i][j].dv.y = buffer[buffer_size+2];
 	          SolnBlk.U[i][j].E    = buffer[buffer_size+3];
                  buffer_size = buffer_size + 4;
 	          SolnBlk.W[i][j] = W(SolnBlk.U[i][j]);
              } /* endfor */
           } /* endfor */
       } /* endif */

       delete []buffer; 
       buffer = NULL;

       ni = 1;
       nj = (SolnBlk.JCu+SolnBlk.Nghost) - (SolnBlk.JCl-SolnBlk.Nghost) + 1;
       buffer = new double[8*ni*nj];

       if (CFFC_Primary_MPI_Processor()) {
          buffer_size = 0;
          for (j  = SolnBlk.JCl-SolnBlk.Nghost ; j <= SolnBlk.JCu+SolnBlk.Nghost ; ++j ) {
 	      buffer[buffer_size  ] = SolnBlk.WoW[j].d;
 	      buffer[buffer_size+1] = SolnBlk.WoW[j].v.x;
 	      buffer[buffer_size+2] = SolnBlk.WoW[j].v.y;
 	      buffer[buffer_size+3] = SolnBlk.WoW[j].p;
 	      buffer[buffer_size+4] = SolnBlk.WoE[j].d;
 	      buffer[buffer_size+5] = SolnBlk.WoE[j].v.x;
 	      buffer[buffer_size+6] = SolnBlk.WoE[j].v.y;
 	      buffer[buffer_size+7] = SolnBlk.WoE[j].p;
              buffer_size = buffer_size + 8;
          } /* endfor */
       } /* endif */

       buffer_size = 8*ni*nj;
       MPI::COMM_WORLD.Bcast(buffer, buffer_size, MPI::DOUBLE, 0);

       if (!CFFC_Primary_MPI_Processor()) {
          buffer_size = 0;
          for (j  = SolnBlk.JCl-SolnBlk.Nghost ; j <= SolnBlk.JCu+SolnBlk.Nghost ; ++j ) {
 	      SolnBlk.WoW[j].d   = buffer[buffer_size];
 	      SolnBlk.WoW[j].v.x = buffer[buffer_size+1];
 	      SolnBlk.WoW[j].v.y = buffer[buffer_size+2];
 	      SolnBlk.WoW[j].p   = buffer[buffer_size+3];
 	      SolnBlk.WoE[j].d   = buffer[buffer_size+4];
 	      SolnBlk.WoE[j].v.x = buffer[buffer_size+5];
 	      SolnBlk.WoE[j].v.y = buffer[buffer_size+6];
 	      SolnBlk.WoE[j].p   = buffer[buffer_size+7];
              buffer_size = buffer_size + 8;
          } /* endfor */
       } /* endif */

       delete []buffer; 
       buffer = NULL;

       ni = (SolnBlk.ICu+SolnBlk.Nghost) - (SolnBlk.ICl-SolnBlk.Nghost) + 1;
       nj = 1;
       buffer = new double[8*ni*nj];

       if (CFFC_Primary_MPI_Processor()) {
          buffer_size = 0;
          for (i  = SolnBlk.ICl-SolnBlk.Nghost ; i <= SolnBlk.ICu+SolnBlk.Nghost ; ++i ) {
 	      buffer[buffer_size  ] = SolnBlk.WoS[i].d;
 	      buffer[buffer_size+1] = SolnBlk.WoS[i].v.x;
 	      buffer[buffer_size+2] = SolnBlk.WoS[i].v.y;
 	      buffer[buffer_size+3] = SolnBlk.WoS[i].p;
 	      buffer[buffer_size+4] = SolnBlk.WoN[i].d;
 	      buffer[buffer_size+5] = SolnBlk.WoN[i].v.x;
 	      buffer[buffer_size+6] = SolnBlk.WoN[i].v.y;
 	      buffer[buffer_size+7] = SolnBlk.WoN[i].p;
              buffer_size = buffer_size + 8;
          } /* endfor */
       } /* endif */

       buffer_size = 8*ni*nj;
       MPI::COMM_WORLD.Bcast(buffer, buffer_size, MPI::DOUBLE, 0);

       if (!CFFC_Primary_MPI_Processor()) {
          buffer_size = 0;
          for (i  = SolnBlk.ICl-SolnBlk.Nghost ; i <= SolnBlk.ICu+SolnBlk.Nghost ; ++i ) {
 	      SolnBlk.WoS[i].d   = buffer[buffer_size];
 	      SolnBlk.WoS[i].v.x = buffer[buffer_size+1];
 	      SolnBlk.WoS[i].v.y = buffer[buffer_size+2];
 	      SolnBlk.WoS[i].p   = buffer[buffer_size+3];
 	      SolnBlk.WoN[i].d   = buffer[buffer_size+4];
 	      SolnBlk.WoN[i].v.x = buffer[buffer_size+5];
 	      SolnBlk.WoN[i].v.y = buffer[buffer_size+6];
 	      SolnBlk.WoN[i].p   = buffer[buffer_size+7];
              buffer_size = buffer_size + 8;
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
/********************************************************
 * Routine: Broadcast_Solution_Block                    *
 *                                                      *
 * Broadcast quadrilateral solution block to all        *
 * processors associated with the specified communicator*
 * from the specified processor using the MPI broadcast *
 * routine.                                             *
 *                                                      *
 ********************************************************/
void Broadcast_Solution_Block(Euler2D_Quad_Block &SolnBlk,
                              MPI::Intracomm &Communicator, 
                              const int Source_CPU) {

    int Source_Rank = 0;
    int i, j, ni, nj, ng, nr, block_allocated, buffer_size;
    double *buffer;

    // High-order related variables
    int NumberOfHighOrderVariables;
    vector<int> ReconstructionOrders;

    /* Broadcast the number of cells in each direction. */

    if (CFFC_MPI::This_Processor_Number == Source_CPU) {
      ni = SolnBlk.NCi;
      nj = SolnBlk.NCj;
      ng = SolnBlk.Nghost;
      nr = SolnBlk.residual_variable;
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
    Communicator.Bcast(&nr,1,MPI::INT,Source_Rank);
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
       // Set the block static variables if they were not previously assigned.
       if (SolnBlk.residual_variable != nr) SolnBlk.residual_variable = nr;
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
       buffer = new double[4*ni*nj];

       if (CFFC_MPI::This_Processor_Number == Source_CPU) {
          buffer_size = 0;
          for (j  = SolnBlk.JCl-SolnBlk.Nghost ; j <= SolnBlk.JCu+SolnBlk.Nghost ; ++j ) {
              for ( i = SolnBlk.ICl-SolnBlk.Nghost ; i <= SolnBlk.ICu+SolnBlk.Nghost ; ++i ) {
 	          buffer[buffer_size  ] = SolnBlk.U[i][j].d;
 	          buffer[buffer_size+1] = SolnBlk.U[i][j].dv.x;
 	          buffer[buffer_size+2] = SolnBlk.U[i][j].dv.y;
 	          buffer[buffer_size+3] = SolnBlk.U[i][j].E;
                  buffer_size = buffer_size + 4;
              } /* endfor */
          } /* endfor */
       } /* endif */

       buffer_size = 4*ni*nj;
       Communicator.Bcast(buffer, buffer_size, MPI::DOUBLE, Source_Rank);

       if (!(CFFC_MPI::This_Processor_Number == Source_CPU)) {
          buffer_size = 0;
          for (j  = SolnBlk.JCl-SolnBlk.Nghost ; j <= SolnBlk.JCu+SolnBlk.Nghost ; ++j ) {
              for ( i = SolnBlk.ICl-SolnBlk.Nghost ; i <= SolnBlk.ICu+SolnBlk.Nghost ; ++i ) {
 	          SolnBlk.U[i][j].d    = buffer[buffer_size];
 	          SolnBlk.U[i][j].dv.x = buffer[buffer_size+1];
 	          SolnBlk.U[i][j].dv.y = buffer[buffer_size+2];
 	          SolnBlk.U[i][j].E    = buffer[buffer_size+3];
                  buffer_size = buffer_size + 4;
 	          SolnBlk.W[i][j] = W(SolnBlk.U[i][j]);
              } /* endfor */
           } /* endfor */
       } /* endif */

       delete []buffer; 
       buffer = NULL;

       ni = 1;
       nj = (SolnBlk.JCu+SolnBlk.Nghost) - (SolnBlk.JCl-SolnBlk.Nghost) + 1;
       buffer = new double[8*ni*nj];

       if (CFFC_MPI::This_Processor_Number == Source_CPU) {
          buffer_size = 0;
          for (j  = SolnBlk.JCl-SolnBlk.Nghost ; j <= SolnBlk.JCu+SolnBlk.Nghost ; ++j ) {
 	      buffer[buffer_size  ] = SolnBlk.WoW[j].d;
 	      buffer[buffer_size+1] = SolnBlk.WoW[j].v.x;
 	      buffer[buffer_size+2] = SolnBlk.WoW[j].v.y;
 	      buffer[buffer_size+3] = SolnBlk.WoW[j].p;
 	      buffer[buffer_size+4] = SolnBlk.WoE[j].d;
 	      buffer[buffer_size+5] = SolnBlk.WoE[j].v.x;
 	      buffer[buffer_size+6] = SolnBlk.WoE[j].v.y;
 	      buffer[buffer_size+7] = SolnBlk.WoE[j].p;
              buffer_size = buffer_size + 8;
          } /* endfor */
       } /* endif */

       buffer_size = 8*ni*nj;
       Communicator.Bcast(buffer, buffer_size, MPI::DOUBLE, Source_Rank);

       if (!(CFFC_MPI::This_Processor_Number == Source_CPU)) {
          buffer_size = 0;
          for (j  = SolnBlk.JCl-SolnBlk.Nghost ; j <= SolnBlk.JCu+SolnBlk.Nghost ; ++j ) {
 	      SolnBlk.WoW[j].d   = buffer[buffer_size];
 	      SolnBlk.WoW[j].v.x = buffer[buffer_size+1];
 	      SolnBlk.WoW[j].v.y = buffer[buffer_size+2];
 	      SolnBlk.WoW[j].p   = buffer[buffer_size+3];
 	      SolnBlk.WoE[j].d   = buffer[buffer_size+4];
 	      SolnBlk.WoE[j].v.x = buffer[buffer_size+5];
 	      SolnBlk.WoE[j].v.y = buffer[buffer_size+6];
 	      SolnBlk.WoE[j].p   = buffer[buffer_size+7];
              buffer_size = buffer_size + 8;
          } /* endfor */
       } /* endif */

       delete []buffer; 
       buffer = NULL;

       ni = (SolnBlk.ICu+SolnBlk.Nghost) - (SolnBlk.ICl-SolnBlk.Nghost) + 1;
       nj = 1;
       buffer = new double[8*ni*nj];

       if (CFFC_MPI::This_Processor_Number == Source_CPU) {
          buffer_size = 0;
          for (i  = SolnBlk.ICl-SolnBlk.Nghost ; i <= SolnBlk.ICu+SolnBlk.Nghost ; ++i ) {
 	      buffer[buffer_size  ] = SolnBlk.WoS[i].d;
 	      buffer[buffer_size+1] = SolnBlk.WoS[i].v.x;
 	      buffer[buffer_size+2] = SolnBlk.WoS[i].v.y;
 	      buffer[buffer_size+3] = SolnBlk.WoS[i].p;
 	      buffer[buffer_size+4] = SolnBlk.WoN[i].d;
 	      buffer[buffer_size+5] = SolnBlk.WoN[i].v.x;
 	      buffer[buffer_size+6] = SolnBlk.WoN[i].v.y;
 	      buffer[buffer_size+7] = SolnBlk.WoN[i].p;
              buffer_size = buffer_size + 8;
          } /* endfor */
       } /* endif */

       buffer_size = 8*ni*nj;
       Communicator.Bcast(buffer, buffer_size, MPI::DOUBLE, Source_Rank);

       if (!(CFFC_MPI::This_Processor_Number == Source_CPU)) {
          buffer_size = 0;
          for (i  = SolnBlk.ICl-SolnBlk.Nghost ; i <= SolnBlk.ICu+SolnBlk.Nghost ; ++i ) {
 	      SolnBlk.WoS[i].d   = buffer[buffer_size];
 	      SolnBlk.WoS[i].v.x = buffer[buffer_size+1];
 	      SolnBlk.WoS[i].v.y = buffer[buffer_size+2];
 	      SolnBlk.WoS[i].p   = buffer[buffer_size+3];
 	      SolnBlk.WoN[i].d   = buffer[buffer_size+4];
 	      SolnBlk.WoN[i].v.x = buffer[buffer_size+5];
 	      SolnBlk.WoN[i].v.y = buffer[buffer_size+6];
 	      SolnBlk.WoN[i].p   = buffer[buffer_size+7];
              buffer_size = buffer_size + 8;
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

/********************************************************
 * Routine: Copy_Solution_Block                         *
 *                                                      *
 * Copies the solution information of quadrilateral     *
 * solution block SolnBlk2 to SolnBlk1.                 *
 *                                                      *
 ********************************************************/
void Copy_Solution_Block(Euler2D_Quad_Block &SolnBlk1,
                         Euler2D_Quad_Block &SolnBlk2) {
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
 * \note Add high-order prolongation! Add mechanism to prolong high-order BCs!
 ********************************************************/
int Prolong_Solution_Block(Euler2D_Quad_Block &SolnBlk_Fine,
			   Euler2D_Quad_Block &SolnBlk_Original,
			   const int Sector) {

    int i, j, i_min, i_max, j_min, j_max, mesh_refinement_permitted;
    double area_total_fine;
    Vector2D dX;

    // High-order related variables
    int SolnBlk_Original_NumberOfHighOrderVariables(SolnBlk_Original.NumberOfHighOrderObjects());
    vector<int> SolnBlk_Original_ReconstructionOrders; 

    /* Allocate (re-allocate) memory for the solution
       of the refined quadrilateral solution block as necessary. */

    if ( (SolnBlk_Original.NCi-2*SolnBlk_Original.Nghost-2*((SolnBlk_Original.NCi-2*SolnBlk_Original.Nghost)/2) != 0) || 
         (SolnBlk_Original.NCj-2*SolnBlk_Original.Nghost-2*((SolnBlk_Original.NCj-2*SolnBlk_Original.Nghost)/2) != 0) ||
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
          if (SolnBlk_Original.U != NULL)
	    SolnBlk_Fine.allocate(SolnBlk_Original.NCi-2*SolnBlk_Original.Nghost,
				  SolnBlk_Original.NCj-2*SolnBlk_Original.Nghost,
				  SolnBlk_Original.Nghost);
          // In this case, create the refined mesh.
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

       for ( j  = j_min; j <= j_max ; ++j ) {
	   for ( i = i_min ; i <= i_max ; ++i ) {
               area_total_fine = SolnBlk_Fine.Grid.Cell[2*(i-i_min)+SolnBlk_Fine.ICl  ]
                                                       [2*(j-j_min)+SolnBlk_Fine.JCl  ].A+
                                 SolnBlk_Fine.Grid.Cell[2*(i-i_min)+SolnBlk_Fine.ICl+1]
                                                       [2*(j-j_min)+SolnBlk_Fine.JCl  ].A+
                                 SolnBlk_Fine.Grid.Cell[2*(i-i_min)+SolnBlk_Fine.ICl  ]
                                                       [2*(j-j_min)+SolnBlk_Fine.JCl+1].A+
                                 SolnBlk_Fine.Grid.Cell[2*(i-i_min)+SolnBlk_Fine.ICl+1]
                                                       [2*(j-j_min)+SolnBlk_Fine.JCl+1].A;

     	       SolnBlk_Fine.U[2*(i-i_min)+SolnBlk_Fine.ICl  ]
                              [2*(j-j_min)+SolnBlk_Fine.JCl  ] 
                   = SolnBlk_Original.U[i][j];
                   //= (SolnBlk_Original.Grid.Cell[i][j].A/area_total_fine)*SolnBlk_Original.U[i][j];
     	       SolnBlk_Fine.W[2*(i-i_min)+SolnBlk_Fine.ICl  ]
                              [2*(j-j_min)+SolnBlk_Fine.JCl  ] 
                   = W(SolnBlk_Fine.U[2*(i-i_min)+SolnBlk_Fine.ICl  ]
                                     [2*(j-j_min)+SolnBlk_Fine.JCl  ]);

     	       SolnBlk_Fine.U[2*(i-i_min)+SolnBlk_Fine.ICl+1]
                              [2*(j-j_min)+SolnBlk_Fine.JCl  ] 
                   = SolnBlk_Original.U[i][j];
                   //= (SolnBlk_Original.Grid.Cell[i][j].A/area_total_fine)*SolnBlk_Original.U[i][j];
     	       SolnBlk_Fine.W[2*(i-i_min)+SolnBlk_Fine.ICl+1]
                              [2*(j-j_min)+SolnBlk_Fine.JCl  ] 
                   = W(SolnBlk_Fine.U[2*(i-i_min)+SolnBlk_Fine.ICl+1]
                                     [2*(j-j_min)+SolnBlk_Fine.JCl  ]);

     	       SolnBlk_Fine.U[2*(i-i_min)+SolnBlk_Fine.ICl  ]
                              [2*(j-j_min)+SolnBlk_Fine.JCl+1]
                   = SolnBlk_Original.U[i][j];
                   //= (SolnBlk_Original.Grid.Cell[i][j].A/area_total_fine)*SolnBlk_Original.U[i][j];
    	       SolnBlk_Fine.W[2*(i-i_min)+SolnBlk_Fine.ICl  ]
                              [2*(j-j_min)+SolnBlk_Fine.JCl+1]
                   = W(SolnBlk_Fine.U[2*(i-i_min)+SolnBlk_Fine.ICl  ]
                                     [2*(j-j_min)+SolnBlk_Fine.JCl+1]);

     	       SolnBlk_Fine.U[2*(i-i_min)+SolnBlk_Fine.ICl+1]
                              [2*(j-j_min)+SolnBlk_Fine.JCl+1]
                   = SolnBlk_Original.U[i][j];
                   //= (SolnBlk_Original.Grid.Cell[i][j].A/area_total_fine)*SolnBlk_Original.U[i][j];
    	       SolnBlk_Fine.W[2*(i-i_min)+SolnBlk_Fine.ICl+1]
                              [2*(j-j_min)+SolnBlk_Fine.JCl+1]
                   = W(SolnBlk_Fine.U[2*(i-i_min)+SolnBlk_Fine.ICl+1]
                                     [2*(j-j_min)+SolnBlk_Fine.JCl+1]);

//                SolnBlk_Original.SubcellReconstruction(i, j, LIMITER_VENKATAKRISHNAN);

//                dX = SolnBlk_Fine.Grid.Cell[2*(i-i_min)+SolnBlk_Fine.ICl  ]
//                                           [2*(j-j_min)+SolnBlk_Fine.JCl  ].Xc -
//                     SolnBlk_Original.Grid.Cell[i][j].Xc;
//     	       SolnBlk_Fine.W[2*(i-i_min)+SolnBlk_Fine.ICl  ]
//                              [2*(j-j_min)+SolnBlk_Fine.JCl  ] 
//                   = SolnBlk_Original.W[i][j] +
//                     (SolnBlk_Original.phi[i][j]^SolnBlk_Original.dWdx[i][j])*dX.x +
//                     (SolnBlk_Original.phi[i][j]^SolnBlk_Original.dWdy[i][j])*dX.y;
//     	       SolnBlk_Fine.U[2*(i-i_min)+SolnBlk_Fine.ICl  ]
//                              [2*(j-j_min)+SolnBlk_Fine.JCl  ] 
//                   = U(SolnBlk_Fine.W[2*(i-i_min)+SolnBlk_Fine.ICl  ]
//                                     [2*(j-j_min)+SolnBlk_Fine.JCl  ]);

//                dX = SolnBlk_Fine.Grid.Cell[2*(i-i_min)+SolnBlk_Fine.ICl+1]
//                                           [2*(j-j_min)+SolnBlk_Fine.JCl  ].Xc -
//                     SolnBlk_Original.Grid.Cell[i][j].Xc;
//     	       SolnBlk_Fine.W[2*(i-i_min)+SolnBlk_Fine.ICl+1]
//                              [2*(j-j_min)+SolnBlk_Fine.JCl  ] 
//                   = SolnBlk_Original.W[i][j] +
//                     (SolnBlk_Original.phi[i][j]^SolnBlk_Original.dWdx[i][j])*dX.x +
//                     (SolnBlk_Original.phi[i][j]^SolnBlk_Original.dWdy[i][j])*dX.y;
//     	       SolnBlk_Fine.U[2*(i-i_min)+SolnBlk_Fine.ICl+1]
//                              [2*(j-j_min)+SolnBlk_Fine.JCl  ] 
//                   = U(SolnBlk_Fine.W[2*(i-i_min)+SolnBlk_Fine.ICl+1]
//                                     [2*(j-j_min)+SolnBlk_Fine.JCl  ]);

//                dX = SolnBlk_Fine.Grid.Cell[2*(i-i_min)+SolnBlk_Fine.ICl  ]
//                                           [2*(j-j_min)+SolnBlk_Fine.JCl+1].Xc -
//                     SolnBlk_Original.Grid.Cell[i][j].Xc;
//     	       SolnBlk_Fine.W[2*(i-i_min)+SolnBlk_Fine.ICl  ]
//                              [2*(j-j_min)+SolnBlk_Fine.JCl+1]
//                   = SolnBlk_Original.W[i][j] +
//                     (SolnBlk_Original.phi[i][j]^SolnBlk_Original.dWdx[i][j])*dX.x +
//                     (SolnBlk_Original.phi[i][j]^SolnBlk_Original.dWdy[i][j])*dX.y;
//    	       SolnBlk_Fine.U[2*(i-i_min)+SolnBlk_Fine.ICl  ]
//                              [2*(j-j_min)+SolnBlk_Fine.JCl+1]
//                   = U(SolnBlk_Fine.W[2*(i-i_min)+SolnBlk_Fine.ICl  ]
//                                     [2*(j-j_min)+SolnBlk_Fine.JCl+1]);

//                dX = SolnBlk_Fine.Grid.Cell[2*(i-i_min)+SolnBlk_Fine.ICl+1]
//                                           [2*(j-j_min)+SolnBlk_Fine.JCl+1].Xc -
//                     SolnBlk_Original.Grid.Cell[i][j].Xc;
//     	       SolnBlk_Fine.W[2*(i-i_min)+SolnBlk_Fine.ICl+1]
//                              [2*(j-j_min)+SolnBlk_Fine.JCl+1]
//                   = SolnBlk_Original.W[i][j] +
//                     (SolnBlk_Original.phi[i][j]^SolnBlk_Original.dWdx[i][j])*dX.x +
//                     (SolnBlk_Original.phi[i][j]^SolnBlk_Original.dWdy[i][j])*dX.y;
//    	       SolnBlk_Fine.U[2*(i-i_min)+SolnBlk_Fine.ICl+1]
//                              [2*(j-j_min)+SolnBlk_Fine.JCl+1]
//                   = U(SolnBlk_Fine.W[2*(i-i_min)+SolnBlk_Fine.ICl+1]
//                                     [2*(j-j_min)+SolnBlk_Fine.JCl+1]);
           } /* endfor */
       } /* endfor */

       for ( j  = j_min-SolnBlk_Original.Nghost/2; j <= j_max+SolnBlk_Original.Nghost/2 ; ++j ) {
           SolnBlk_Fine.WoW[2*(j-j_min)+SolnBlk_Fine.JCl  ]
              = SolnBlk_Original.WoW[j];
           SolnBlk_Fine.WoW[2*(j-j_min)+SolnBlk_Fine.JCl+1]
              = SolnBlk_Original.WoW[j];

           SolnBlk_Fine.WoE[2*(j-j_min)+SolnBlk_Fine.JCl  ]
              = SolnBlk_Original.WoE[j];
           SolnBlk_Fine.WoE[2*(j-j_min)+SolnBlk_Fine.JCl+1]
              = SolnBlk_Original.WoE[j];
       } /* endfor */

       for ( i = i_min-SolnBlk_Original.Nghost/2 ; i <= i_max+SolnBlk_Original.Nghost/2 ; ++i ) {
           SolnBlk_Fine.WoS[2*(i-i_min)+SolnBlk_Fine.ICl  ]
              = SolnBlk_Original.WoS[i];
           SolnBlk_Fine.WoS[2*(i-i_min)+SolnBlk_Fine.ICl+1]
              = SolnBlk_Original.WoS[i];

           SolnBlk_Fine.WoN[2*(i-i_min)+SolnBlk_Fine.ICl  ]
              = SolnBlk_Original.WoN[i];
           SolnBlk_Fine.WoN[2*(i-i_min)+SolnBlk_Fine.ICl+1]
              = SolnBlk_Original.WoN[i];
       } /* endfor */

       // Set the reference values for the boundary states to the ones from the Original solution block
       SolnBlk_Fine.Set_Reference_Values_For_Boundary_States(SolnBlk_Original.Ref_State_BC_North,
							     SolnBlk_Original.Ref_State_BC_South,
							     SolnBlk_Original.Ref_State_BC_East,
							     SolnBlk_Original.Ref_State_BC_West);

    } /* endif */

    // Prolongation of solution block was successful.
    return 0;

}

/********************************************************
 * Routine: Restrict_Solution_Block                     *
 *                                                      *
 * Restricts the solution information of four original  *
 * fine quadrilateral solution blocks to the coarse     *
 * solution block SolnBlk_Coarse.                       *
 *                                                      *
 ********************************************************/
int Restrict_Solution_Block(Euler2D_Quad_Block &SolnBlk_Coarse,
			    Euler2D_Quad_Block &SolnBlk_Original_SW,
			    Euler2D_Quad_Block &SolnBlk_Original_SE,
			    Euler2D_Quad_Block &SolnBlk_Original_NW,
			    Euler2D_Quad_Block &SolnBlk_Original_NE) {

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
             SolnBlk_Coarse.U[i_coarse][j_coarse] = (SolnBlk_Original_SW.Grid.Cell[i  ][j  ].A*
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
                                                     SolnBlk_Original_SW.Grid.Cell[i+1][j+1].A);
                                                    //SolnBlk_Coarse.Grid.Cell[i_coarse][j_coarse].A;
             SolnBlk_Coarse.W[i_coarse][j_coarse] = W(SolnBlk_Coarse.U[i_coarse][j_coarse]);
          } /* endfor */
      } /* endfor */

      for ( j = SolnBlk_Original_SW.JCl-SolnBlk_Original_SW.Nghost; 
            j <= SolnBlk_Original_SW.JCu+SolnBlk_Original_SW.Nghost ; j += 2 ) {
      	  j_coarse = (j-SolnBlk_Original_SW.JCl)/2+SolnBlk_Coarse.JCl;
          SolnBlk_Coarse.WoW[j_coarse] = SolnBlk_Original_SW.WoW[j];
          if (j == SolnBlk_Original_SW.JCl-SolnBlk_Original_SW.Nghost) {
             SolnBlk_Coarse.WoW[j_coarse-1] = SolnBlk_Original_SW.WoW[j];
          } /* endif */
      } /* endfor */

      for ( i = SolnBlk_Original_SW.ICl-SolnBlk_Original_SW.Nghost ; 
            i <= SolnBlk_Original_SW.ICu+SolnBlk_Original_SW.Nghost ; i += 2) {
      	  i_coarse = (i-SolnBlk_Original_SW.ICl)/2+SolnBlk_Coarse.ICl;
          SolnBlk_Coarse.WoS[i_coarse] = SolnBlk_Original_SW.WoS[i];
          if (i == SolnBlk_Original_SW.ICl-SolnBlk_Original_SW.Nghost) {
             SolnBlk_Coarse.WoS[i_coarse-1] = SolnBlk_Original_SW.WoS[i];
          } /* endif */
      } /* endfor */

      // South-East corner fine block:
      for ( j = SolnBlk_Original_SE.JCl; j <= SolnBlk_Original_SE.JCu ; j += 2 ) {
      	  for ( i = SolnBlk_Original_SE.ICl ; i <= SolnBlk_Original_SE.ICu ; i += 2 ) {
      	     i_coarse = (i-SolnBlk_Original_SE.ICl)/2+
                        (SolnBlk_Coarse.ICu-SolnBlk_Coarse.ICl+1)/2+SolnBlk_Coarse.ICl;
      	     j_coarse = (j-SolnBlk_Original_SE.JCl)/2+
                        SolnBlk_Coarse.JCl;
             SolnBlk_Coarse.U[i_coarse][j_coarse] = (SolnBlk_Original_SE.Grid.Cell[i  ][j  ].A*
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
                                                     SolnBlk_Original_SE.Grid.Cell[i+1][j+1].A);
                                                    //SolnBlk_Coarse.Grid.Cell[i_coarse][j_coarse].A;
             SolnBlk_Coarse.W[i_coarse][j_coarse] = W(SolnBlk_Coarse.U[i_coarse][j_coarse]);
          } /* endfor */
      } /* endfor */

      for ( j = SolnBlk_Original_SE.JCl-SolnBlk_Original_SE.Nghost; 
            j <= SolnBlk_Original_SE.JCu+SolnBlk_Original_SE.Nghost ; j += 2 ) {
      	  j_coarse = (j-SolnBlk_Original_SE.JCl)/2+SolnBlk_Coarse.JCl;
          SolnBlk_Coarse.WoE[j_coarse] = SolnBlk_Original_SE.WoE[j];
          if (j == SolnBlk_Original_SE.JCl-SolnBlk_Original_SE.Nghost) {
             SolnBlk_Coarse.WoE[j_coarse-1] = SolnBlk_Original_SE.WoE[j];
          } /* endif */
      } /* endfor */

      for ( i = SolnBlk_Original_SE.ICl-SolnBlk_Original_SE.Nghost ; 
            i <= SolnBlk_Original_SE.ICu+SolnBlk_Original_SE.Nghost ; i += 2) {
     	  i_coarse = (i-SolnBlk_Original_SE.ICl)/2+
                     (SolnBlk_Coarse.ICu-SolnBlk_Coarse.ICl+1)/2+SolnBlk_Coarse.ICl;
          SolnBlk_Coarse.WoS[i_coarse] = SolnBlk_Original_SE.WoS[i];
          if (i == SolnBlk_Original_SE.ICu+SolnBlk_Original_SE.Nghost) {
             SolnBlk_Coarse.WoS[i_coarse+1] = SolnBlk_Original_SE.WoS[i];
          } /* endif */
      } /* endfor */

      // North-West corner fine block:
      for ( j = SolnBlk_Original_NW.JCl; j <= SolnBlk_Original_NW.JCu ; j += 2 ) {
      	  for ( i = SolnBlk_Original_NW.ICl ; i <= SolnBlk_Original_NW.ICu ; i += 2 ) {
      	     i_coarse = (i-SolnBlk_Original_NW.ICl)/2+
                        SolnBlk_Coarse.ICl;
      	     j_coarse = (j-SolnBlk_Original_NW.JCl)/2+
                        (SolnBlk_Coarse.JCu-SolnBlk_Coarse.JCl+1)/2+SolnBlk_Coarse.JCl;
             SolnBlk_Coarse.U[i_coarse][j_coarse] = (SolnBlk_Original_NW.Grid.Cell[i  ][j  ].A*
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
                                                     SolnBlk_Original_NW.Grid.Cell[i+1][j+1].A);
                                                     //SolnBlk_Coarse.Grid.Cell[i_coarse][j_coarse].A;
             SolnBlk_Coarse.W[i_coarse][j_coarse] = W(SolnBlk_Coarse.U[i_coarse][j_coarse]);
          } /* endfor */
      } /* endfor */

      for ( j = SolnBlk_Original_NW.JCl-SolnBlk_Original_NW.Nghost; 
            j <= SolnBlk_Original_NW.JCu+SolnBlk_Original_NW.Nghost ; j += 2 ) {
      	  j_coarse = (j-SolnBlk_Original_NW.JCl)/2+
                     (SolnBlk_Coarse.JCu-SolnBlk_Coarse.JCl+1)/2+SolnBlk_Coarse.JCl;
          SolnBlk_Coarse.WoW[j_coarse] = SolnBlk_Original_NW.WoW[j];
          if (j == SolnBlk_Original_NW.JCu+SolnBlk_Original_NW.Nghost) {
             SolnBlk_Coarse.WoW[j_coarse+1] = SolnBlk_Original_NW.WoW[j];
          } /* endif */
      } /* endfor */

      for ( i = SolnBlk_Original_NW.ICl-SolnBlk_Original_NW.Nghost ; 
            i <= SolnBlk_Original_NW.ICu+SolnBlk_Original_NW.Nghost ; i += 2) {
      	  i_coarse = (i-SolnBlk_Original_NW.ICl)/2+SolnBlk_Coarse.ICl;
          SolnBlk_Coarse.WoN[i_coarse] = SolnBlk_Original_NW.WoN[i];
          if (i == SolnBlk_Original_NW.ICl-SolnBlk_Original_NW.Nghost) {
             SolnBlk_Coarse.WoN[i_coarse-1] = SolnBlk_Original_NW.WoN[i];
          } /* endif */
      } /* endfor */

      // North-East corner fine block:
      for ( j = SolnBlk_Original_NE.JCl; j <= SolnBlk_Original_NE.JCu ; j += 2 ) {
      	  for ( i = SolnBlk_Original_NE.ICl ; i <= SolnBlk_Original_NE.ICu ; i += 2 ) {
      	     i_coarse = (i-SolnBlk_Original_NE.ICl)/2+
                        (SolnBlk_Coarse.ICu-SolnBlk_Coarse.ICl+1)/2+SolnBlk_Coarse.ICl;
      	     j_coarse = (j-SolnBlk_Original_NE.JCl)/2+
                        (SolnBlk_Coarse.JCu-SolnBlk_Coarse.JCl+1)/2+SolnBlk_Coarse.JCl;
             SolnBlk_Coarse.U[i_coarse][j_coarse] = (SolnBlk_Original_NE.Grid.Cell[i  ][j  ].A*
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
                                                     SolnBlk_Original_NE.Grid.Cell[i+1][j+1].A);
	                                             //SolnBlk_Coarse.Grid.Cell[i_coarse][j_coarse].A;
             SolnBlk_Coarse.W[i_coarse][j_coarse] = W(SolnBlk_Coarse.U[i_coarse][j_coarse]);
          } /* endfor */
      } /* endfor */

      for ( j = SolnBlk_Original_NE.JCl-SolnBlk_Original_NE.Nghost; 
            j <= SolnBlk_Original_NE.JCu+SolnBlk_Original_NE.Nghost ; j += 2 ) {
      	  j_coarse = (j-SolnBlk_Original_NE.JCl)/2+
                     (SolnBlk_Coarse.JCu-SolnBlk_Coarse.JCl+1)/2+SolnBlk_Coarse.JCl;
          SolnBlk_Coarse.WoE[j_coarse] = SolnBlk_Original_NE.WoE[j];
          if (j == SolnBlk_Original_NE.JCu+SolnBlk_Original_NE.Nghost) {
             SolnBlk_Coarse.WoE[j_coarse+1] = SolnBlk_Original_NE.WoE[j];
          } /* endif */
      } /* endfor */

      for ( i = SolnBlk_Original_NE.ICl-SolnBlk_Original_NE.Nghost ; 
            i <= SolnBlk_Original_NE.ICu+SolnBlk_Original_NE.Nghost ; i += 2) {
      	  i_coarse = (i-SolnBlk_Original_NE.ICl)/2+
                     (SolnBlk_Coarse.ICu-SolnBlk_Coarse.ICl+1)/2+SolnBlk_Coarse.ICl;
          SolnBlk_Coarse.WoN[i_coarse] = SolnBlk_Original_NE.WoN[i];
          if (i == SolnBlk_Original_NE.ICu+SolnBlk_Original_NE.Nghost) {
             SolnBlk_Coarse.WoN[i_coarse+1] = SolnBlk_Original_NE.WoN[i];
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

/********************************************************
 * Routine: Output_Tecplot                              *
 *                                                      *
 * Writes the solution values at the nodes of the       *
 * specified quadrilateral solution block to the        *
 * specified output stream suitable for plotting with   *
 * TECPLOT.                                             *
 *                                                      *
 ********************************************************/
void Output_Tecplot(Euler2D_Quad_Block &SolnBlk,
		    Euler2D_Input_Parameters &IP,
                    const int Number_of_Time_Steps,
                    const double &Time,
                    const int Block_Number,
                    const int Output_Title,
	            ostream &Out_File) {

  int i, j, nRow, nLoop;
  Euler2D_pState W_node;
  Vector2D Node;


  if (Tecplot_Execution_Mode::IsSmoothNodalSolnOutputRequired()){

    // Output the interpolated nodal solution

    /* Ensure boundary conditions are updated before
       evaluating solution at the nodes. */
    
    BCs(SolnBlk,IP);

    /* Output node solution data. */

    Out_File << setprecision(14);
    if (Output_Title) {
      Out_File << "TITLE = \"" << CFFC_Name() << ": 2D Euler Solution, "
	       << "Time Step/Iteration Level = " << Number_of_Time_Steps
	       << ", Time = " << Time
	       << "\"" << "\n"
	       << "VARIABLES = \"x\" \\ \n"
	       << "\"y\" \\ \n"
	       << "\"rho\" \\ \n"
	       << "\"u\" \\ \n"
	       << "\"v\" \\ \n"
	       << "\"p\" \\ \n"
	       << "\"T\" \n"
	       << "\"M\" \n"
	       << "\"H\" \n"
	       << "\"s\" \n"
	       << "ZONE T =  \"Block Number = " << Block_Number
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
	W_node = SolnBlk.Wn(i, j);
	Out_File << " "  << SolnBlk.Grid.Node[i][j].X << W_node;
	Out_File.setf(ios::scientific);
	Out_File << " " << W_node.T() << " " << W_node.v.abs()/W_node.a() 
		 << " " << W_node.H() << " " << W_node.s() << "\n";
	Out_File.unsetf(ios::scientific);
      } /* endfor */
    } /* endfor */
    Out_File << setprecision(6);
    
  } else {

    // Output the discontinuous nodal solution

    /* Output node solution data. */

    Out_File << setprecision(14);
    if (Output_Title) {
      Out_File << "TITLE = \"" << CFFC_Name() << ": 2D Euler Solution, "
	       << "Time Step/Iteration Level = " << Number_of_Time_Steps
	       << ", Time = " << Time
	       << "\"" << "\n"
	       << "VARIABLES = \"x\" \\ \n"
	       << "\"y\" \\ \n"
	       << "\"rho\" \\ \n"
	       << "\"u\" \\ \n"
	       << "\"v\" \\ \n"
	       << "\"p\" \\ \n";
	
      // Add more variables for the Detailed format
      if (Tecplot_Execution_Mode::IsDetailedOutputRequired()){
	Out_File << "\"T\" \n"
		 << "\"M\" \n"
		 << "\"H\" \n"
		 << "\"s\" \n";
      }

      Out_File << "ZONE T =  \"Block Number = " << Block_Number
	       << "\" \\ \n"
	       << "I = " << (SolnBlk.Grid.ICu - SolnBlk.Grid.ICl + 1)*3 << " \\ \n"
	       << "J = " << (SolnBlk.Grid.JCu - SolnBlk.Grid.JCl + 1)*3 << " \\ \n"
	       << "F = POINT \n";
    } else {
      Out_File << "ZONE T =  \"Block Number = " << Block_Number
	       << "\" \\ \n"
	       << "I = " << (SolnBlk.Grid.ICu - SolnBlk.Grid.ICl + 1)*3 << " \\ \n"
	       << "J = " << (SolnBlk.Grid.JCu - SolnBlk.Grid.JCl + 1)*3 << " \\ \n"
	       << "F = POINT \n";
    } /* endif */

    // Set the accuracy properly
    if (Tecplot_Execution_Mode::IsDoublePrecision()){
      Out_File << "DT = (DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE ";

      // Detail format
      if (Tecplot_Execution_Mode::IsDetailedOutputRequired()){
	Out_File << "DOUBLE DOUBLE DOUBLE DOUBLE ";
      }

      // Close line
      Out_File << " ) \n";
    } // endif (DoublePrecision)

    // Output data
    for ( j  = SolnBlk.Grid.JCl ; j <= SolnBlk.Grid.JCu ; ++j ) { // for every j cell
      for ( nRow = 1; nRow <= 3; ++nRow){ // for 3 rows of nodes
	for ( i = SolnBlk.Grid.ICl ; i <= SolnBlk.Grid.ICu ; ++i ) { // for every i cell  
	  for (nLoop = 1; nLoop <= 3; ++nLoop){	// for every node
	    // Get the node location
	    switch(nRow){
	    case 1: // output the 1st row of nodes (i.e. NodeSW(i,j), xfaceS(i,j), NodeSE(i,j))
	      switch(nLoop){
	      case 1:		// output NodeSW(i,j)
		Node = SolnBlk.Grid.nodeSW(i,j).X;
		break;
	      case 2:		// output xfaceS(i,j)
		Node = SolnBlk.Grid.xfaceS(i,j);
		break;
	      case 3:		// output NodeSE(i,j)
		Node = SolnBlk.Grid.nodeSE(i,j).X;
		break;
	      }
	      break;

	    case 2: // output the 2nd row of nodes (i.e. xfaceW(i,j), Grid.CellCentroid(i,j), xfaceE(i,j))
	      switch(nLoop){
	      case 1:		// output xfaceW(i,j)
		Node = SolnBlk.Grid.xfaceW(i,j);
		break;
	      case 2:		// output Grid.CellCentroid(i,j)
		Node = SolnBlk.Grid.CellCentroid(i,j);
		break;
	      case 3:		// output xfaceE(i,j) 
		Node = SolnBlk.Grid.xfaceE(i,j);
		break;
	      }
	      break;

	    case 3: // output the 3rd row of nodes (i.e. NodeNW(i,j), xfaceN(i,j), NodeNE(i,j))
	      switch(nLoop){
	      case 1:		// output NodeNW(i,j)
		Node = SolnBlk.Grid.nodeNW(i,j).X;
		break;
	      case 2:		// output xfaceN(i,j)
		Node = SolnBlk.Grid.xfaceN(i,j);
		break;
	      case 3:		// output NodeNE(i,j)
		Node = SolnBlk.Grid.nodeNE(i,j).X;
		break;
	      }
	      break;
	    } // endswitch

	    // Output Brief format
	    W_node = SolnBlk.PiecewiseLinearSolutionAtLocation(i,j,Node);
	    Out_File << " "  << Node 
		     << " "  << W_node;

	    // Add more variables for the Detailed format
	    if (Tecplot_Execution_Mode::IsDetailedOutputRequired()){ 
	      Out_File.setf(ios::scientific);
	      Out_File << " " << W_node.T() 
		       << " " << W_node.v.abs()/W_node.a() 
		       << " " << W_node.H() 
		       << " " << W_node.s();
	      Out_File.unsetf(ios::scientific);
	    }

	    // Close line
	    Out_File << "\n";
	    Out_File.unsetf(ios::scientific);

	  }
	} /* endfor */
      }
    } /* endfor */
    Out_File << setprecision(6);  

  } // endif (Tecplot_Execution_Mode::IsSmoothNodalSolnOutputRequired())

}

/********************************************************
 * Routine: Output_Cells_Tecplot                        *
 *                                                      *
 * Writes the cell centred solution values of the       *
 * specified quadrilateral solution block to the        *
 * specified output stream suitable for plotting with   *
 * TECPLOT.                                             *
 *                                                      *
 ********************************************************/
void Output_Cells_Tecplot(Euler2D_Quad_Block &SolnBlk,
                          Euler2D_Input_Parameters &IP,
                          const int Number_of_Time_Steps,
                          const double &Time,
                          const int Block_Number,
                          const int Output_Title,
	                  ostream &Out_File) {

    int i, j;

    Out_File << setprecision(14);
    if (Output_Title) {
       Out_File << "TITLE = \"" << CFFC_Name() << ": 2D Euler Solution, "
                << "Time Step/Iteration Level = " << Number_of_Time_Steps
                << ", Time = " << Time
                << "\"" << "\n"
   	        << "VARIABLES = \"x\" \\ \n"
                << "\"y\" \\ \n"
                << "\"rho\" \\ \n"
                << "\"u\" \\ \n"
                << "\"v\" \\ \n"
                << "\"p\" \\ \n"
                << "\"T\" \n"
                << "\"M\" \n"
                << "\"H\" \n"
                << "\"s\" \n"
                << "ZONE T =  \"Block Number = " << Block_Number
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
           Out_File << " "  << SolnBlk.Grid.Cell[i][j].Xc
                    << SolnBlk.W[i][j];
           Out_File.setf(ios::scientific);
           Out_File << " " << SolnBlk.W[i][j].T() << " " << SolnBlk.W[i][j].v.abs()/SolnBlk.W[i][j].a() 
                    << " " << SolnBlk.W[i][j].H() << " " << SolnBlk.W[i][j].s();
	   Out_File << "\n";
           Out_File.unsetf(ios::scientific);
       } /* endfor */
    } /* endfor */
    Out_File << setprecision(6);
    
}

/********************************************************
 * Routine: Output_Nodes_Tecplot                        *
 *                                                      *
 * Writes the node values of the specified              *
 * quadrilateral solution block to the specified output *
 * stream suitable for plotting with TECPLOT.           *
 *                                                      *
 ********************************************************/
void Output_Nodes_Tecplot(Euler2D_Quad_Block &SolnBlk,
                          const int Number_of_Time_Steps,
                          const double &Time,
                          const int Block_Number,
                          const int Output_Title,
	                  ostream &Out_File) {

    int i, j;

    Out_File << setprecision(14);
    if (Output_Title) {
       Out_File << "TITLE = \"" << CFFC_Name() << ": 2D Euler Solution, "
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
    for (int j = SolnBlk.Grid.JNl - SolnBlk.Nghost; j <= SolnBlk.Grid.JNu + SolnBlk.Nghost; j++) {
      for (int i = SolnBlk.Grid.INl - SolnBlk.Nghost; i <= SolnBlk.Grid.INu + SolnBlk.Nghost; i++) {
	Out_File << " " << SolnBlk.Grid.Node[i][j].X << endl;
      }
    }
    Out_File.unsetf(ios::scientific);
    Out_File << setprecision(6);
    
}

/**********************************************************************
 * Routine: Output_Gradients_Tecplot                                  *
 *                                                                    *
 * Writes the cell centred primitive variable state gradients and     *
 * limiters of the specified quadrilateral solution block to the      *
 * specified output stream suitable for plotting with TECPLOT.        *
 *                                                                    *
 **********************************************************************/
void Output_Gradients_Tecplot(Euler2D_Quad_Block &SolnBlk,
			      const int Number_of_Time_Steps,
			      const double &Time,
			      const int Block_Number,
			      const int Output_Title,
			      ostream &Out_File) {

  Out_File << setprecision(14);
  if (Output_Title) {
    Out_File << "TITLE = \"" << CFFC_Name() << ": 2D Euler Solution, "
	     << "Time Step/Iteration Level = " << Number_of_Time_Steps
	     << ", Time = " << Time
	     << "\"" << "\n"
	     << "VARIABLES = \"x\" \\ \n"
 	     << "\"y\" \\ \n"
 	     << "\"drhodx\" \\ \n"
 	     << "\"drhody\" \\ \n"
 	     << "\"dudx\" \\ \n"
 	     << "\"dudy\" \\ \n"
 	     << "\"dvdx\" \\ \n"
 	     << "\"dvdy\" \\ \n"
 	     << "\"dpdx\" \\ \n"
 	     << "\"dpdy\" \\ \n"
 	     << "\"phi_rho\" \\ \n"
 	     << "\"phi_u\" \\ \n"
 	     << "\"phi_v\" \\ \n"
 	     << "\"phi_p\" \\ \n";
  }
  Out_File << "ZONE T =  \"Block Number = " << Block_Number
	   << "\" \\ \n"
	   << "I = " << SolnBlk.Grid.ICu - SolnBlk.Grid.ICl + 3 << " \\ \n"
	   << "J = " << SolnBlk.Grid.JCu - SolnBlk.Grid.JCl + 3 << " \\ \n"
	   << "F = POINT \n";

  for (int j = SolnBlk.JCl-1; j <= SolnBlk.JCu+1; j++) {
    for (int i = SolnBlk.ICl-1; i <= SolnBlk.ICu+1; i++) {
      Linear_Reconstruction_LeastSquares(SolnBlk,i,j,LIMITER_VENKATAKRISHNAN);
      Out_File.setf(ios::scientific);
      Out_File << " " << SolnBlk.Grid.Cell[i][j].Xc
	       << " " << SolnBlk.dWdx[i][j].d
 	       << " " << SolnBlk.dWdy[i][j].d
 	       << " " << SolnBlk.dWdx[i][j].v.x
 	       << " " << SolnBlk.dWdy[i][j].v.x
 	       << " " << SolnBlk.dWdx[i][j].v.y
 	       << " " << SolnBlk.dWdy[i][j].v.y
 	       << " " << SolnBlk.dWdx[i][j].p
 	       << " " << SolnBlk.dWdy[i][j].p
	       << " " << SolnBlk.phi[i][j].d
 	       << " " << SolnBlk.phi[i][j].v.x
 	       << " " << SolnBlk.phi[i][j].v.y
 	       << " " << SolnBlk.phi[i][j].p
	       << endl;
      Out_File.unsetf(ios::scientific);
    }
  }
  Out_File << setprecision(6);

}

/**********************************************************************
 * Routine: Output_Tecplot_Quasi3D                                    *
 *                                                                    *
 * Writes the solution values at the nodes of the specified           *
 * active/positive quadrilateral solution block to the specified      *
 * output stream suitable for plotting with TECPLOT in a quasi-3D     *
 * format.                                                            *
 *                                                                    *
 **********************************************************************/
void Output_Tecplot_Quasi3D(Euler2D_Quad_Block &SolnBlk,
			    Euler2D_Input_Parameters &IP,
			    const int Number_of_Time_Steps,
			    const double &Time,
			    const int Block_Number,
			    const int Output_Title,
			    ostream &Out_File) {

  Euler2D_pState W_node;
  int numberofrotations = 360/15;
  int numberofnodes = ((SolnBlk.Grid.INu - SolnBlk.Grid.INl + 1)*
		       (SolnBlk.Grid.JNu - SolnBlk.Grid.JNl + 1));
  int numberofcells = ((SolnBlk.Grid.ICu - SolnBlk.Grid.ICl + 1)*
		       (SolnBlk.Grid.JCu - SolnBlk.Grid.JCl + 1));

  // Ensure boundary conditions are updated before evaluating
  // solution at the nodes.
  BCs(SolnBlk,IP);

  // Output node solution data.  
  Out_File << setprecision(14);
  if (Output_Title) {
    Out_File << "TITLE = \"" << CFFC_Name() << ": Quasi 3D Euler Solution, "
	     << "Time Step/Iteration Level = " << Number_of_Time_Steps
	     << ", Time = " << Time
	     << "\"" << "\n"
	     << "VARIABLES = \"x\" \\ \n"
	     << "\"y\" \\ \n"
	     << "\"z\" \\ \n"
	     << "\"rho\" \\ \n"
	     << "\"vx\" \\ \n"
	     << "\"vy\" \\ \n"
	     << "\"vz\" \\ \n"
	     << "\"p\" \\ \n"
	     << "\"T\" \\ \n"
	     << "\"M\" \\ \n"
	     << "\"H\" \\ \n"
	     << "\"s\" \n";
  }
  
  for (int r = 0; r < numberofrotations; r++) {
    Out_File << "ZONE T =  \"Block Number = " << Block_Number << r
	     << "\" \\ \n"
	     << "N = " << 2*numberofnodes << " \\ \n"
	     << "E = " << numberofcells << " \\ \n"
	     << "F = FEPOINT \\ \n"
	     << "ET = BRICK \n";
    for (int j = SolnBlk.Grid.JNl; j <= SolnBlk.Grid.JNu; j++) {
      for (int i = SolnBlk.Grid.INl; i <= SolnBlk.Grid.INu; i++) {
	W_node = SolnBlk.Wn(i,j);
	Out_File << " " << SolnBlk.Grid.Node[i][j].X.x
		 << " " << SolnBlk.Grid.Node[i][j].X.y*sin(TWO*PI*double(r)/double(numberofrotations))
		 << " " << SolnBlk.Grid.Node[i][j].X.y*cos(TWO*PI*double(r)/double(numberofrotations))
		 << " " << W_node.d
		 << " " << W_node.v.x
		 << " " << W_node.v.y*sin(TWO*PI*double(r)/double(numberofrotations))
		 << " " << W_node.v.y*cos(TWO*PI*double(r)/double(numberofrotations))
		 << " " << W_node.p
		 << " " << W_node.T()
		 << " " << W_node.v.abs()/W_node.a()
		 << " " << W_node.H() 
		 << " " << W_node.s()
		 << endl;
      }
    }
    if (r < numberofrotations-1) {
      for (int j = SolnBlk.Grid.JNl; j <= SolnBlk.Grid.JNu; j++) {
	for (int i = SolnBlk.Grid.INl; i <= SolnBlk.Grid.INu; i++) {
	  W_node = SolnBlk.Wn(i,j);
	  Out_File << " " << SolnBlk.Grid.Node[i][j].X.x
		   << " " << SolnBlk.Grid.Node[i][j].X.y*sin(TWO*PI*double(r+1)/double(numberofrotations))
		   << " " << SolnBlk.Grid.Node[i][j].X.y*cos(TWO*PI*double(r+1)/double(numberofrotations))
		   << " " << W_node.d
		   << " " << W_node.v.x
		   << " " << W_node.v.y*sin(TWO*PI*double(r+1)/double(numberofrotations))
		   << " " << W_node.v.y*cos(TWO*PI*double(r+1)/double(numberofrotations))
		   << " " << W_node.p
		   << " " << W_node.T()
		   << " " << W_node.v.abs()/W_node.a()
		   << " " << W_node.H() 
		   << " " << W_node.s()
		   << endl;
	}
      }
    } else {
      for (int j = SolnBlk.Grid.JNl; j <= SolnBlk.Grid.JNu; j++) {
	for (int i = SolnBlk.Grid.INl; i <= SolnBlk.Grid.INu; i++) {
	  W_node = SolnBlk.Wn(i,j);
	  Out_File << " " << SolnBlk.Grid.Node[i][j].X.x
		   << " " << SolnBlk.Grid.Node[i][j].X.y*sin(TWO*PI*double(0)/double(numberofrotations))
		   << " " << SolnBlk.Grid.Node[i][j].X.y*cos(TWO*PI*double(0)/double(numberofrotations))
		   << " " << W_node.d
		   << " " << W_node.v.x
		   << " " << W_node.v.y*sin(TWO*PI*double(0)/double(numberofrotations))
		   << " " << W_node.v.y*cos(TWO*PI*double(0)/double(numberofrotations))
		   << " " << W_node.p
		   << " " << W_node.T()
		   << " " << W_node.v.abs()/W_node.a()
		   << " " << W_node.H() 
		   << " " << W_node.s()
		   << endl;
	}
      }
    }
    // Connectivity table.
    for (int j = SolnBlk.Grid.JCl; j <= SolnBlk.Grid.JCu; j++) {
      for (int i = SolnBlk.Grid.ICl; i <= SolnBlk.Grid.ICu; i++) {
	Out_File << (j-2)*(SolnBlk.Grid.INu-1) + i - 1                         << " ";
	Out_File << (j-2)*(SolnBlk.Grid.INu-1) + i                             << " ";
	Out_File << (j-1)*(SolnBlk.Grid.INu-2) + i + 1 + (j-2)                 << " ";
	Out_File << (j-1)*(SolnBlk.Grid.INu-2) + i     + (j-2)                 << " ";
	Out_File << (j-2)*(SolnBlk.Grid.INu-1) + i - 1         + numberofnodes << " ";
	Out_File << (j-2)*(SolnBlk.Grid.INu-1) + i             + numberofnodes << " ";
	Out_File << (j-1)*(SolnBlk.Grid.INu-2) + i + 1 + (j-2) + numberofnodes << " ";
	Out_File << (j-1)*(SolnBlk.Grid.INu-2) + i     + (j-2) + numberofnodes << endl;
      }
    }
  
  }

  Out_File << setprecision(6);

}

/**********************************************************************
 * Routine: Output_Ringleb_Flow_Solution                              *
 *                                                                    *
 *                                                                    *
 **********************************************************************/
void Output_Ringleb_Flow_Solution(Euler2D_Quad_Block &SolnBlk,
				  const int Block_Number,
				  const int Output_Title,
				  ostream &Out_File) {

  Euler2D_pState W;

  // Output node solution data.  
  Out_File << setprecision(14);
  if (Output_Title)
    Out_File << "TITLE = \"" << CFFC_Name() << ": 2D Euler Ringleb Flow Solution "
	     << "\"" << "\n"
	     << "VARIABLES = \"x\" \\ \n"
	     << "\"y\" \\ \n"
	     << "\"rho\" \\ \n"
	     << "\"a\" \\ \n"
	     << "\"rho_exact\" \\ \n"
	     << "\"a_exact\" \\ \n"
	     << "\"rho_error\" \\ \n"
	     << "\"a_error\" \\ \n";
  Out_File << "ZONE T =  \"Block Number = " << Block_Number
	   << "\" \\ \n"
	   << "I = " << SolnBlk.Grid.ICu - SolnBlk.Grid.ICl + 1 << " \\ \n"
	   << "J = " << SolnBlk.Grid.JCu - SolnBlk.Grid.JCl + 1 << " \\ \n"
	   << "F = POINT \\ \n";

  Out_File.setf(ios::scientific);
  for (int j = SolnBlk.Grid.JCl; j <= SolnBlk.Grid.JCu; j++) {
    for (int i = SolnBlk.Grid.ICl; i <= SolnBlk.Grid.ICu; i++) {
      W = RinglebFlow(W,SolnBlk.Grid.Cell[i][j].Xc);
      Out_File << SolnBlk.Grid.Cell[i][j].Xc << " " 
	       << SolnBlk.W[i][j].d << " " 
	       << SolnBlk.W[i][j].a() << " " 
	       << W.d << " " 
	       << W.a()   << " " 
	       << fabs(SolnBlk.W[i][j].d - W.d) << " " 
	       << fabs(SolnBlk.W[i][j].a() - W.a()) << endl;
    } /* endfor */
  } /* endfor */

}

/**********************************************************************
 * Routine: Output_Ringleb_Flow_Error                                 *
 *                                                                    *
 *                                                                    *
 **********************************************************************/
void Output_Ringleb_Flow_Error(Euler2D_Quad_Block &SolnBlk,
			       double &l1_norm,
			       double &l2_norm,
			       double &max_norm) {

  Euler2D_pState W;

   l1_norm = ZERO;
   l2_norm = ZERO;
   max_norm = ZERO;
   for (int j = SolnBlk.Grid.JCl; j <= SolnBlk.Grid.JCu; j++) {
     for (int i = SolnBlk.Grid.ICl; i <= SolnBlk.Grid.ICu; i++) {
       W = RinglebFlow(W,SolnBlk.Grid.Cell[i][j].Xc);
       l1_norm += fabs(W.d - SolnBlk.W[i][j].d);
       l2_norm += sqr(W.d - SolnBlk.W[i][j].d);
       max_norm = max(max_norm,fabs(W.d - SolnBlk.W[i][j].d));
     }
   }
   l1_norm = l1_norm/((SolnBlk.NCi-2*SolnBlk.Nghost)*(SolnBlk.NCj-2*SolnBlk.Nghost));
   l2_norm = sqrt(l2_norm/((SolnBlk.NCi-2*SolnBlk.Nghost)*(SolnBlk.NCj-2*SolnBlk.Nghost)));

//   Vector2D X, X1, X2, X3, X4;
//   double *w;
//   Euler2D_pState **f;
//   double *epsilon;
//   double *eta;

//   double epsilon1 = -ONE, epsilon2 = ONE, epsilon3 = -ONE, epsilon4 = ONE;
//   double eta1 = -ONE, eta2 = -ONE, eta3 =  ONE, eta4 = ONE;
//   double N1, N2, N3, N4;

//   l1_norm = ZERO;
//   l2_norm = ZERO;
//   max_norm = ZERO;
 
//   // Allocate memory.
//   w = new double[3];
//   f = new Euler2D_pState*[3];
//   epsilon = new double[3];
//   eta = new double[3];
//   for (int i = 0; i < 3; i++) f[i] = new Euler2D_pState[3];

//   // Assign values... the dumbest way possible.
//   w[0] = 5.0/9.0;
//   w[1] = 8.0/9.0;
//   w[2] = 5.0/9.0;
//   epsilon[0] = -sqrt(15.0)/5.0;
//   epsilon[1] = 0.0;
//   epsilon[2] = sqrt(15.0)/5.0;
//   eta[0] = -sqrt(15.0)/5.0;
//   eta[1] = 0.0;
//   eta[2] = sqrt(15.0)/5.0;

//   // Get the value of 'f' at each point then perform the integration and finally compute the norm.
//   for (int j = SolnBlk.Grid.JCl; j <= SolnBlk.Grid.JCu; j++) {
//     for (int i = SolnBlk.Grid.ICl; i <= SolnBlk.Grid.ICu; i++) {

//       // Get the value of 'f' at each point.
//       for (int jj = 0; jj < 3; jj++) {
// 	for (int ii = 0; ii < 3; ii++) {

// 	  // Save node values.
// 	  X1 = SolnBlk.Grid.Node[i  ][j  ].X;
// 	  X2 = SolnBlk.Grid.Node[i+1][j  ].X;
// 	  X3 = SolnBlk.Grid.Node[i+1][j+1].X;
// 	  X4 = SolnBlk.Grid.Node[i  ][j+1].X;

// 	  // Set basis functions.
// 	  N1 = 0.25*(1 + epsilon[ii]*epsilon1)*(1 + eta[jj]*eta1);
// 	  N2 = 0.25*(1 + epsilon[ii]*epsilon2)*(1 + eta[jj]*eta2);
// 	  N3 = 0.25*(1 + epsilon[ii]*epsilon3)*(1 + eta[jj]*eta3);
// 	  N4 = 0.25*(1 + epsilon[ii]*epsilon4)*(1 + eta[jj]*eta4);

// 	  // Get point X.
// 	  X = X1*N1 + X2*N2 + X3*N3 + X4*N4;

// 	  // Determine the value of 'f'.
// 	  f[ii][jj] = RinglebFlow(X);
	  
// 	}
//       }

//       // Perform the integration to find the average value of the exact value of rho.
//       W = Euler2D_W_VACUUM;
//       for (int jj = 0; jj < 3; jj++)
// 	for (int ii = 0; ii < 3; ii++)
// 	  W += w[ii]*w[jj]*f[ii][jj]/FOUR;

//       // Compute the norms.
//       l1_norm += fabs(W.d - SolnBlk.W[i][j].d);
//       l2_norm += sqr(W.d - SolnBlk.W[i][j].d);
//       max_norm = max(max_norm,fabs(W.d - SolnBlk.W[i][j].d));

//     }
//   }
//   l1_norm = l1_norm/((SolnBlk.NCi-2*SolnBlk.Nghost)*(SolnBlk.NCj-2*SolnBlk.Nghost));
//   l2_norm = sqrt(l2_norm/((SolnBlk.NCi-2*SolnBlk.Nghost)*(SolnBlk.NCj-2*SolnBlk.Nghost)));

//   // Delete memory.
//   for (int i = 0; i < 3; i++) {
//     delete []f[i]; f[i] = NULL;
//   }
//   delete []w; w = NULL;
//   delete []f; f = NULL;
//   delete []epsilon; epsilon = NULL;
//   delete []eta; eta = NULL;


}

/********************************************************
 * Routine: ICs                                         *
 *                                                      *
 * Assigns initial conditions and data to the           *
 * solution variables of the specified quadrilateral    *
 * solution block.                                      *
 *                                                      *
 ********************************************************/
void ICs(Euler2D_Quad_Block &SolnBlk,
         Euler2D_Input_Parameters &IP,
         Euler2D_pState *Wo) {

    int i, j, k;
    Euler2D_pState Wl, Wr, Wm;

    /* Assign the initial data for the IVP of interest. */

    switch(IP.i_ICs) {
      case IC_CONSTANT :
      case IC_UNIFORM :
        // Set the solution state to the initial state Wo[0].
        for (j  = SolnBlk.JCl-SolnBlk.Nghost ; j <= SolnBlk.JCu+SolnBlk.Nghost ; ++j ) {
            for ( i = SolnBlk.ICl-SolnBlk.Nghost ; i <= SolnBlk.ICu+SolnBlk.Nghost ; ++i ) {
	       SolnBlk.W[i][j] = Wo[0];
               SolnBlk.U[i][j] = U(SolnBlk.W[i][j]);
            } /* endfor */
        } /* endfor */
        break;
      case IC_SOD_XDIR :
        // Set initial data for Sod IVP in x-direction.
        Wl = Euler2D_W_STDATM;
        Wr = Euler2D_pState(DENSITY_STDATM/EIGHT,
          		    ZERO, ZERO,
         		    PRESSURE_STDATM/TEN);
        for (j  = SolnBlk.JCl-SolnBlk.Nghost ; j <= SolnBlk.JCu+SolnBlk.Nghost ; ++j ) {
            for ( i = SolnBlk.ICl-SolnBlk.Nghost ; i <= SolnBlk.ICu+SolnBlk.Nghost ; ++i ) {
               if (SolnBlk.Grid.Cell[i][j].Xc.x <= ZERO) {
                  SolnBlk.W[i][j] = Wl;
               } else {
                  SolnBlk.W[i][j] = Wr;	     
               } /* end if */
               SolnBlk.U[i][j] = U(SolnBlk.W[i][j]);
            } /* endfor */
        } /* endfor */
        break;
      case IC_SOD_YDIR :
        // Set initial data for Sod IVP in y-direction.
        Wl = Euler2D_W_STDATM;
        Wr = Euler2D_pState(DENSITY_STDATM/EIGHT,
          		    ZERO, ZERO,
         		    PRESSURE_STDATM/TEN); 
        for (j  = SolnBlk.JCl-SolnBlk.Nghost ; j <= SolnBlk.JCu+SolnBlk.Nghost ; ++j ) {
            for ( i = SolnBlk.ICl-SolnBlk.Nghost ; i <= SolnBlk.ICu+SolnBlk.Nghost ; ++i ) {
               if (SolnBlk.Grid.Cell[i][j].Xc.y <= ZERO) {
                  SolnBlk.W[i][j] = Wl;
               } else {
                  SolnBlk.W[i][j] = Wr;	
               } /* end if */
               SolnBlk.U[i][j] = U(SolnBlk.W[i][j]);
            } /* endfor */
        } /* endfor */
        break;
      case IC_GROTH_XDIR :
        // Set initial data for Groth IVP in x-direction.
        Wl = Euler2D_pState(4.696, ZERO, ZERO, 404.4e03);
        Wr = Euler2D_pState(1.408, ZERO, ZERO, 101.1e03);
        for (j  = SolnBlk.JCl-SolnBlk.Nghost ; j <= SolnBlk.JCu+SolnBlk.Nghost ; ++j ) {
            for ( i = SolnBlk.ICl-SolnBlk.Nghost ; i <= SolnBlk.ICu+SolnBlk.Nghost ; ++i ) {
               if (SolnBlk.Grid.Cell[i][j].Xc.x <= ZERO) {
                  SolnBlk.W[i][j] = Wl;
               } else {
                  SolnBlk.W[i][j] = Wr;	     
               } /* end if */
               SolnBlk.U[i][j] = U(SolnBlk.W[i][j]);
            } /* endfor */
        } /* endfor */
        break;
      case IC_GROTH_YDIR :
        // Set initial data for Groth IVP in y-direction.
        Wl = Euler2D_pState(4.696, ZERO, ZERO, 404.4e03);
        Wr = Euler2D_pState(1.408, ZERO, ZERO, 101.1e03);
        for (j  = SolnBlk.JCl-SolnBlk.Nghost ; j <= SolnBlk.JCu+SolnBlk.Nghost ; ++j ) {
            for ( i = SolnBlk.ICl-SolnBlk.Nghost ; i <= SolnBlk.ICu+SolnBlk.Nghost ; ++i ) {
               if (SolnBlk.Grid.Cell[i][j].Xc.y <= ZERO) {
                  SolnBlk.W[i][j] = Wl;
               } else {
                  SolnBlk.W[i][j] = Wr;	
               } /* end if */
               SolnBlk.U[i][j] = U(SolnBlk.W[i][j]);
            } /* endfor */
        } /* endfor */
        break;
      case IC_EINFELDT_XDIR :
        // Set initial data for Einfeldt IVP in x-direction.
        Wl = Euler2D_pState(ONE, -TWO, ZERO, FOUR/TEN);
        Wr = Euler2D_pState(ONE, TWO, ZERO, FOUR/TEN);
        for (j  = SolnBlk.JCl-SolnBlk.Nghost ; j <= SolnBlk.JCu+SolnBlk.Nghost ; ++j ) {
            for ( i = SolnBlk.ICl-SolnBlk.Nghost ; i <= SolnBlk.ICu+SolnBlk.Nghost ; ++i ) {
               if (SolnBlk.Grid.Cell[i][j].Xc.x <= ZERO) {
                  SolnBlk.W[i][j] = Wl;
               } else {
                  SolnBlk.W[i][j] = Wr;	
               } /* end if */
               SolnBlk.U[i][j] = U(SolnBlk.W[i][j]);
            } /* endfor */
        } /* endfor */
        break;
      case IC_EINFELDT_YDIR :
        // Set initial data for Einfeld IVP in y-direction.
        Wl = Euler2D_pState(ONE, ZERO, -TWO, FOUR/TEN);
        Wr = Euler2D_pState(ONE, ZERO, TWO, FOUR/TEN);
        for (j  = SolnBlk.JCl-SolnBlk.Nghost ; j <= SolnBlk.JCu+SolnBlk.Nghost ; ++j ) {
            for ( i = SolnBlk.ICl-SolnBlk.Nghost ; i <= SolnBlk.ICu+SolnBlk.Nghost ; ++i ) {
               if (SolnBlk.Grid.Cell[i][j].Xc.y <= ZERO) {
                  SolnBlk.W[i][j] = Wl;
               } else {
                  SolnBlk.W[i][j] = Wr;	
               } /* end if */
               SolnBlk.U[i][j] = U(SolnBlk.W[i][j]);
            } /* endfor */
        } /* endfor */
        break;
      case IC_SHOCK_WAVE_XDIR :
        // Set initial data for moving shock wave propagating in x-direction.
        //Wl = Euler2D_pState(2.281, 164.83, ZERO, 201.17e03);
        //Wr = Euler2D_pState(1.408, ZERO, ZERO, 101.1e03);
        for (j  = SolnBlk.JCl-SolnBlk.Nghost ; j <= SolnBlk.JCu+SolnBlk.Nghost ; ++j ) {
            for ( i = SolnBlk.ICl-SolnBlk.Nghost ; i <= SolnBlk.ICu+SolnBlk.Nghost ; ++i ) {
               if (SolnBlk.Grid.Cell[i][j].Xc.x <= IP.Wave_Position.x) {
		  SolnBlk.W[i][j] = IP.Wo;
               } else {
		  SolnBlk.W[i][j] = Euler2D_W_STDATM;
               } /* end if */
               SolnBlk.U[i][j] = U(SolnBlk.W[i][j]);
            } /* endfor */
        } /* endfor */
        break;
      case IC_SHOCK_WAVE_YDIR :
        // Set initial data for moving shock wave propagating in y-direction.
        Wl = Euler2D_pState(2.281, ZERO, 164.83, 201.17e03);
        Wr = Euler2D_pState(1.408, ZERO, ZERO, 101.1e03);
        for (j  = SolnBlk.JCl-SolnBlk.Nghost ; j <= SolnBlk.JCu+SolnBlk.Nghost ; ++j ) {
            for ( i = SolnBlk.ICl-SolnBlk.Nghost ; i <= SolnBlk.ICu+SolnBlk.Nghost ; ++i ) {
               if (SolnBlk.Grid.Cell[i][j].Xc.y <= IP.Wave_Position.y) {
                  SolnBlk.W[i][j] = Wl;
               } else {
                  SolnBlk.W[i][j] = Wr;	
               } /* end if */
               SolnBlk.U[i][j] = U(SolnBlk.W[i][j]);
            } /* endfor */
        } /* endfor */
        break;
      case IC_CONTACT_SURFACE_XDIR :
        // Set initial data for moving contact surface propagating in x-direction.
        Wl = Euler2D_pState(1.045, 200.00, ZERO, 300.00e03);
        Wr = Euler2D_pState(3.483, 200.00, ZERO, 300.00e03);
        for (j  = SolnBlk.JCl-SolnBlk.Nghost ; j <= SolnBlk.JCu+SolnBlk.Nghost ; ++j ) {
            for ( i = SolnBlk.ICl-SolnBlk.Nghost ; i <= SolnBlk.ICu+SolnBlk.Nghost ; ++i ) {
               if (SolnBlk.Grid.Cell[i][j].Xc.x <= ZERO) {
                  SolnBlk.W[i][j] = Wl;
               } else {
                  SolnBlk.W[i][j] = Wr;	
               } /* end if */
               SolnBlk.U[i][j] = U(SolnBlk.W[i][j]);
            } /* endfor */
        } /* endfor */
        break;
      case IC_CONTACT_SURFACE_YDIR :
        // Set initial data for moving contact surface propagating in y-direction.
        Wl = Euler2D_pState(1.045, ZERO, 200.00, 300.00e03);
        Wr = Euler2D_pState(3.483, ZERO, 200.00, 300.00e03);
        for (j  = SolnBlk.JCl-SolnBlk.Nghost ; j <= SolnBlk.JCu+SolnBlk.Nghost ; ++j ) {
            for ( i = SolnBlk.ICl-SolnBlk.Nghost ; i <= SolnBlk.ICu+SolnBlk.Nghost ; ++i ) {
               if (SolnBlk.Grid.Cell[i][j].Xc.y <= ZERO) {
                  SolnBlk.W[i][j] = Wl;
               } else {
                  SolnBlk.W[i][j] = Wr;	
               } /* end if */
               SolnBlk.U[i][j] = U(SolnBlk.W[i][j]);
            } /* endfor */
        } /* endfor */
        break;
      case IC_RAREFACTION_WAVE_XDIR :
        // Set initial data for moving rarefaction wave propagating in x-direction.
        Wl = Euler2D_pState(1.598, -383.64, ZERO, 91.88e03);
        Wr = Euler2D_pState(2.787, -216.97, ZERO, 200.0e03);
        for (j  = SolnBlk.JCl-SolnBlk.Nghost ; j <= SolnBlk.JCu+SolnBlk.Nghost ; ++j ) {
            for ( i = SolnBlk.ICl-SolnBlk.Nghost ; i <= SolnBlk.ICu+SolnBlk.Nghost ; ++i ) {
               if (SolnBlk.Grid.Cell[i][j].Xc.x <= ZERO) {
                  SolnBlk.W[i][j] = Wl;
               } else {
                  SolnBlk.W[i][j] = Wr;	
               } /* end if */
               SolnBlk.U[i][j] = U(SolnBlk.W[i][j]);
            } /* endfor */
        } /* endfor */
        break;
      case IC_RAREFACTION_WAVE_YDIR :
        // Set initial data for moving rarefaction wave propagating in y-direction.
        Wl = Euler2D_pState(1.598, ZERO, -383.64, 91.88e03);
        Wr = Euler2D_pState(2.787, ZERO, -216.97, 200.0e03);
        for (j  = SolnBlk.JCl-SolnBlk.Nghost ; j <= SolnBlk.JCu+SolnBlk.Nghost ; ++j ) {
            for ( i = SolnBlk.ICl-SolnBlk.Nghost ; i <= SolnBlk.ICu+SolnBlk.Nghost ; ++i ) {
               if (SolnBlk.Grid.Cell[i][j].Xc.y <= ZERO) {
                  SolnBlk.W[i][j] = Wl;
               } else {
                  SolnBlk.W[i][j] = Wr;	
               } /* end if */
               SolnBlk.U[i][j] = U(SolnBlk.W[i][j]);
            } /* endfor */
        } /* endfor */
        break;
      case IC_SHOCK_BOX :
        // Set initial data for Aki shock-box IVP.
        Wl = Euler2D_W_STDATM;
        Wr = Euler2D_pState(DENSITY_STDATM*FOUR,
         		    ZERO, ZERO,
         		    PRESSURE_STDATM*FOUR);
        for (j  = SolnBlk.JCl-SolnBlk.Nghost ; j <= SolnBlk.JCu+SolnBlk.Nghost ; ++j ) {
            for ( i = SolnBlk.ICl-SolnBlk.Nghost ; i <= SolnBlk.ICu+SolnBlk.Nghost ; ++i ) {
               if (SolnBlk.Grid.Cell[i][j].Xc.x <= ZERO &&
                   SolnBlk.Grid.Cell[i][j].Xc.y <= ZERO) {
                  SolnBlk.W[i][j] = Wl;
               } else {
                  SolnBlk.W[i][j] = Wr;	     
               } /* end if */
               SolnBlk.U[i][j] = U(SolnBlk.W[i][j]);
            } /* endfor */
        } /* endfor */
        break;
      case IC_HIGH_PRESSURE_RESERVOIR :
        // Set high-pressure reservoir initial data.
        Wr = Euler2D_W_STDATM;
        Wl = Euler2D_pState(HUNDRED*DENSITY_STDATM,
          		    ZERO, ZERO,
         		    HUNDRED*PRESSURE_STDATM);
        for (j  = SolnBlk.JCl-SolnBlk.Nghost ; j <= SolnBlk.JCu+SolnBlk.Nghost ; ++j ) {
            for ( i = SolnBlk.ICl-SolnBlk.Nghost ; i <= SolnBlk.ICu+SolnBlk.Nghost ; ++i ) {
               if (SolnBlk.Grid.Cell[i][j].Xc.x <= ZERO) {
                  SolnBlk.W[i][j] = Wl;
               } else {
                  SolnBlk.W[i][j] = Wr;	     
               } /* end if */
               SolnBlk.U[i][j] = U(SolnBlk.W[i][j]);
            } /* endfor */
        } /* endfor */
        break;
      case IC_LOW_PRESSURE_RESERVOIR :
        // Set low-pressure reservoir initial data.
        Wl = Euler2D_W_STDATM;
        Wr = Euler2D_pState(DENSITY_STDATM/THOUSAND,
          		    ZERO, ZERO,
         		    PRESSURE_STDATM/THOUSAND);
        for (j  = SolnBlk.JCl-SolnBlk.Nghost ; j <= SolnBlk.JCu+SolnBlk.Nghost ; ++j ) {
            for ( i = SolnBlk.ICl-SolnBlk.Nghost ; i <= SolnBlk.ICu+SolnBlk.Nghost ; ++i ) {
               if (SolnBlk.Grid.Cell[i][j].Xc.x <= IP.Wave_Position.x) {
                  SolnBlk.W[i][j] = Wl;
               } else {
                  SolnBlk.W[i][j] = Wr;	     
               } /* end if */
               SolnBlk.U[i][j] = U(SolnBlk.W[i][j]);
            } /* endfor */
        } /* endfor */
        break;
      case IC_RIEMANN :
      case IC_RIEMANN_XDIR :
        // Set Riemann initial data 
        // (2-state initial data, left to right).
        Wr = Wo[1];
        Wl = Wo[0];
        for (j  = SolnBlk.JCl-SolnBlk.Nghost ; j <= SolnBlk.JCu+SolnBlk.Nghost ; ++j ) {
            for ( i = SolnBlk.ICl-SolnBlk.Nghost ; i <= SolnBlk.ICu+SolnBlk.Nghost ; ++i ) {
               if (SolnBlk.Grid.Cell[i][j].Xc.x <= ZERO) {
                  SolnBlk.W[i][j] = Wl;
               } else {
                  SolnBlk.W[i][j] = Wr;	     
               } /* end if */
               SolnBlk.U[i][j] = U(SolnBlk.W[i][j]);
            } /* endfor */
        } /* endfor */
        break;
      case IC_RIEMANN_YDIR :
        // Set Riemann initial data 
        // (2-state initial data, top to bottom).
        Wr = Wo[1];
        Wl = Wo[0];
        for (j  = SolnBlk.JCl-SolnBlk.Nghost ; j <= SolnBlk.JCu+SolnBlk.Nghost ; ++j ) {
            for ( i = SolnBlk.ICl-SolnBlk.Nghost ; i <= SolnBlk.ICu+SolnBlk.Nghost ; ++i ) {
               if (SolnBlk.Grid.Cell[i][j].Xc.y <= ZERO) {
                  SolnBlk.W[i][j] = Wl;
               } else {
                  SolnBlk.W[i][j] = Wr;	     
               } /* end if */
               SolnBlk.U[i][j] = U(SolnBlk.W[i][j]);
            } /* endfor */
        } /* endfor */
        break;
      case IC_CYLINDRICAL_EXPLOSION :
	Wl = Euler2D_W_STDATM;
	Wr = Euler2D_pState(EIGHT*DENSITY_STDATM,ZERO,ZERO,TEN*PRESSURE_STDATM);
        for (j  = SolnBlk.JCl-SolnBlk.Nghost ; j <= SolnBlk.JCu+SolnBlk.Nghost ; ++j ) {
            for ( i = SolnBlk.ICl-SolnBlk.Nghost ; i <= SolnBlk.ICu+SolnBlk.Nghost ; ++i ) {
	      if (abs(SolnBlk.Grid.Cell[i][j].Xc) < 2.00) {//0.25) {
		  SolnBlk.W[i][j] = Wr;
	       } else {
		  SolnBlk.W[i][j] = Wl;
	       } /* end if */
	       SolnBlk.U[i][j] = U(SolnBlk.W[i][j]);
	    } /* endfor */
	} /* endfor */
	break;
      case IC_CYLINDRICAL_IMPLOSION :
	Wl = Euler2D_pState(EIGHT*DENSITY_STDATM,ZERO,ZERO,TEN*PRESSURE_STDATM);
	Wr = Euler2D_W_STDATM;
        for (j  = SolnBlk.JCl-SolnBlk.Nghost ; j <= SolnBlk.JCu+SolnBlk.Nghost ; ++j ) {
            for ( i = SolnBlk.ICl-SolnBlk.Nghost ; i <= SolnBlk.ICu+SolnBlk.Nghost ; ++i ) {
	       if (abs(SolnBlk.Grid.Cell[i][j].Xc) < 0.25) {
		  SolnBlk.W[i][j] = Wr;
	       } else {
		  SolnBlk.W[i][j] = Wl;
	       } /* end if */
	       SolnBlk.U[i][j] = U(SolnBlk.W[i][j]);
	    } /* endfor */
	} /* endfor */
	break;
      case IC_WEDGE_FLOW :
	Wl.d = FOUR*DENSITY_STDATM;   Wr.d = DENSITY_STDATM;
	Wl.v = Vector2D_ZERO;         Wr.v = Vector2D_ZERO;
	Wl.p = FOUR*PRESSURE_STDATM;  Wr.p = PRESSURE_STDATM;
        for (j  = SolnBlk.JCl-SolnBlk.Nghost ; j <= SolnBlk.JCu+SolnBlk.Nghost ; ++j ) {
            for ( i = SolnBlk.ICl-SolnBlk.Nghost ; i <= SolnBlk.ICu+SolnBlk.Nghost ; ++i ) {
	       if (SolnBlk.Grid.centroid(i,j).x <= -0.75) {
		  SolnBlk.W[i][j] = Wl;
	       } else {
		  SolnBlk.W[i][j] = Wr;	     
	       } /* end if */
	       SolnBlk.U[i][j] = U(SolnBlk.W[i][j]);
	    } /* endfor */
	} /* endfor */
	break;
      case IC_UNSTEADY_BLUNT_BODY :
	Wl.d = FOUR*DENSITY_STDATM;   Wr.d = DENSITY_STDATM;
	Wl.v = Vector2D_ZERO;         Wr.v = Vector2D_ZERO;
	Wl.p = FOUR*PRESSURE_STDATM;  Wr.p = PRESSURE_STDATM;
        for (j  = SolnBlk.JCl-SolnBlk.Nghost ; j <= SolnBlk.JCu+SolnBlk.Nghost ; ++j ) {
            for ( i = SolnBlk.ICl-SolnBlk.Nghost ; i <= SolnBlk.ICu+SolnBlk.Nghost ; ++i ) {
	       if (SolnBlk.Grid.centroid(i,j).x <= -FIVE) {
		  SolnBlk.W[i][j] = Wl;
	       } else {
		  SolnBlk.W[i][j] = Wr;	     
	       } /* end if */
	       SolnBlk.U[i][j] = U(SolnBlk.W[i][j]);
	    } /* endfor */
	} /* endfor */
	break;
      case IC_RINGLEB_FLOW :
        for (j = SolnBlk.JCl-SolnBlk.Nghost; j <= SolnBlk.JCu+SolnBlk.Nghost; j++) {
            for (i = SolnBlk.ICl-SolnBlk.Nghost; i <= SolnBlk.ICu+SolnBlk.Nghost; i++) {
	      //SolnBlk.W[i][j] = RinglebFlow(SolnBlk.W[i][j],SolnBlk.Grid.Cell[i][j].Xc);
	      SolnBlk.W[i][j] = RinglebFlowAverageState(SolnBlk.W[i][j],
							SolnBlk.Grid.nodeSW(i,j).X,
							SolnBlk.Grid.nodeSE(i,j).X,
							SolnBlk.Grid.nodeNE(i,j).X,
							SolnBlk.Grid.nodeNW(i,j).X);
	      SolnBlk.U[i][j] = U(SolnBlk.W[i][j]);
            } /* endfor */
        } /* endfor */
	break;
      case IC_BLAST_WAVE_INTERACTION:
	// defined in the domain '0' to '1'
	Wl = Euler2D_pState(ONE, ZERO, ZERO, 1.0e03);
	Wm = Euler2D_pState(ONE, ZERO, ZERO, 1.0e-02);
	Wr = Euler2D_pState(ONE, ZERO, ZERO, 1.0e02);
	for (j = SolnBlk.JCl-SolnBlk.Nghost; j <= SolnBlk.JCu+SolnBlk.Nghost; j++) {
	  for (i = SolnBlk.ICl-SolnBlk.Nghost; i <= SolnBlk.ICu+SolnBlk.Nghost; i++) {
	    if (SolnBlk.Grid.Cell[i][j].Xc.x <= 0.1) {
	      SolnBlk.W[i][j] = Wl;
	    } else if (SolnBlk.Grid.Cell[i][j].Xc.x <= 0.9){
	      SolnBlk.W[i][j] = Wm;
	    } else {
	      SolnBlk.W[i][j] = Wr;	
	    } /* end if */
	    SolnBlk.U[i][j] = U(SolnBlk.W[i][j]);
	  } /* endfor */
	} /* endfor */
	break;
      case IC_SHOCK_ACOUSTIC_INTERACTION:
	double density;
	Wl = Euler2D_pState(3.857143, 2.629369, ZERO, 10.333333);
	for (j = SolnBlk.JCl-SolnBlk.Nghost; j <= SolnBlk.JCu+SolnBlk.Nghost; j++) {
	  for (i = SolnBlk.ICl-SolnBlk.Nghost; i <= SolnBlk.ICu+SolnBlk.Nghost; i++) {
	    if (SolnBlk.Grid.Cell[i][j].Xc.x <= -4.0) {
	      SolnBlk.W[i][j] = Wl;
	    } else {
	      density = 1 + 0.2*sin(5*SolnBlk.Grid.Cell[i][j].Xc.x);
	      Wr = Euler2D_pState(density, ZERO, ZERO, ONE);
	      SolnBlk.W[i][j] = Wr;	
	    } /* end if */
	    SolnBlk.U[i][j] = U(SolnBlk.W[i][j]);
	  } /* endfor */
	} /* endfor */
	break;
      case IC_PERIODIC_SINX_WAVE :
	Wl.v.x = 100.0;
	Wl.v.y = 0.0;
	Wl.p = PRESSURE_STDATM;
	for (j = SolnBlk.JCl-SolnBlk.Nghost; j <= SolnBlk.JCu+SolnBlk.Nghost; ++j) {
	  for (i = SolnBlk.ICl-SolnBlk.Nghost; i <= SolnBlk.ICu+SolnBlk.Nghost; ++i) {
	    Wl.d = SolnBlk.Grid.Integration.IntegrateFunctionOverCell(i,j,
								      SinVariationInXDir,
								      IP.Exact_Integration_Digits,Wl.d)/SolnBlk.Grid.Cell[i][j].A;
	    SolnBlk.W[i][j] = Wl;
	    SolnBlk.U[i][j] = U(SolnBlk.W[i][j]);
	  }
	}
	break;
      case IC_PERIODIC_SINY_WAVE :
	Wl.v.x = 0.0;
	Wl.v.y = 100.0;
	Wl.p = PRESSURE_STDATM;
        for (j = SolnBlk.JCl-SolnBlk.Nghost; j <= SolnBlk.JCu+SolnBlk.Nghost; j++) {
	  for (i = SolnBlk.ICl-SolnBlk.Nghost; i <= SolnBlk.ICu+SolnBlk.Nghost; i++) {
	    Wl.d = SolnBlk.Grid.Integration.IntegrateFunctionOverCell(i,j,
								      SinVariationInYDir,
								      IP.Exact_Integration_Digits,Wl.d)/SolnBlk.Grid.Cell[i][j].A;
	    SolnBlk.W[i][j] = Wl;
	    SolnBlk.U[i][j] = U(SolnBlk.W[i][j]);
	  }
	}
      break;
      case IC_ABGRALL_FUNCTION :
	Wl.v.x = 0.0;
	Wl.v.y = 0.0;
	Wl.p = PRESSURE_STDATM;
	for (j = SolnBlk.JCl-SolnBlk.Nghost; j <= SolnBlk.JCu+SolnBlk.Nghost; ++j) {
	  for (i = SolnBlk.ICl-SolnBlk.Nghost; i <= SolnBlk.ICu+SolnBlk.Nghost; ++i) {
	    Wl.d = SolnBlk.Grid.Integration.IntegrateFunctionOverCell(i,j,
								      Abgrall_2D_Function,
								      IP.Exact_Integration_Digits,Wl.d)/SolnBlk.Grid.Cell[i][j].A;
	    SolnBlk.W[i][j] = Wl;
	    SolnBlk.U[i][j] = U(SolnBlk.W[i][j]);
	  }
	}
	break;
      case IC_SIN_EXP_X_WAVE :
	Wl.v.x = 0.0;
	Wl.v.y = 0.0;
	Wl.p = PRESSURE_STDATM;
	for (j = SolnBlk.JCl-SolnBlk.Nghost; j <= SolnBlk.JCu+SolnBlk.Nghost; ++j) {
	  for (i = SolnBlk.ICl-SolnBlk.Nghost; i <= SolnBlk.ICu+SolnBlk.Nghost; ++i) {
	    Wl.d = SolnBlk.Grid.Integration.IntegrateFunctionOverCell(i,j,
								      SinExponentialVariationInXDir,
								      IP.Exact_Integration_Digits,Wl.d)/SolnBlk.Grid.Cell[i][j].A;
	    SolnBlk.W[i][j] = Wl;
	    SolnBlk.U[i][j] = U(SolnBlk.W[i][j]);
	  }
	}
	break;
      case IC_SIN_EXP_Y_WAVE :
	Wl.v.x = 0.0;
	Wl.v.y = 0.0;
	Wl.p = PRESSURE_STDATM;
	for (j = SolnBlk.JCl-SolnBlk.Nghost; j <= SolnBlk.JCu+SolnBlk.Nghost; ++j) {
	  for (i = SolnBlk.ICl-SolnBlk.Nghost; i <= SolnBlk.ICu+SolnBlk.Nghost; ++i) {
	    Wl.d = SolnBlk.Grid.Integration.IntegrateFunctionOverCell(i,j,
								      SinExponentialVariationInYDir,
								      IP.Exact_Integration_Digits,Wl.d)/SolnBlk.Grid.Cell[i][j].A;
	    SolnBlk.W[i][j] = Wl;
	    SolnBlk.U[i][j] = U(SolnBlk.W[i][j]);
	  }
	}
	break;
      case IC_SIN_EXP_ROTATED_WAVE :
	Wl.v.x = 0.0;
	Wl.v.y = 0.0;
	Wl.p = PRESSURE_STDATM;
	for (j = SolnBlk.JCl-SolnBlk.Nghost; j <= SolnBlk.JCu+SolnBlk.Nghost; ++j) {
	  for (i = SolnBlk.ICl-SolnBlk.Nghost; i <= SolnBlk.ICu+SolnBlk.Nghost; ++i) {
	    Wl.d = SolnBlk.Grid.Integration.IntegrateFunctionOverCell(i,j,
								      SinExponentialVariationRotated,
								      IP.Exact_Integration_Digits,Wl.d)/SolnBlk.Grid.Cell[i][j].A;
	    SolnBlk.W[i][j] = Wl;
	    SolnBlk.U[i][j] = U(SolnBlk.W[i][j]);
	  }
	}
	break;
      case IC_COSINE_HILL:
	Wl.v.x = 100.0;
	Wl.v.y = 50.0;
	Wl.p = PRESSURE_STDATM;
	for (j = SolnBlk.JCl-SolnBlk.Nghost; j <= SolnBlk.JCu+SolnBlk.Nghost; j++) {
	  for (i = SolnBlk.ICl-SolnBlk.Nghost; i <= SolnBlk.ICu+SolnBlk.Nghost; i++) {
	    Wl.d = SolnBlk.Grid.Integration.IntegrateFunctionOverCell(i,j,
								      CosineHill,
								      IP.Exact_Integration_Digits,Wl.d)/SolnBlk.Grid.Cell[i][j].A;
	    SolnBlk.W[i][j] = Wl;
	    SolnBlk.U[i][j] = U(SolnBlk.W[i][j]);
	  }
	}
	break;
      case IC_HYPER_TANGENT:
	Wl.v.x = Wo[0].v.x;
	Wl.v.y = Wo[0].v.y;
	Wl.p = Wo[0].p;
	for (j = SolnBlk.JCl-SolnBlk.Nghost; j <= SolnBlk.JCu+SolnBlk.Nghost; j++) {
	  for (i = SolnBlk.ICl-SolnBlk.Nghost; i <= SolnBlk.ICu+SolnBlk.Nghost; i++) {
	    Wl.d = SolnBlk.Grid.Integration.IntegrateFunctionOverCell(i,j,
								      Translated_Solutions::HyperTangentIC,
								      IP.Exact_Integration_Digits,Wl.d)/SolnBlk.Grid.Cell[i][j].A;
	    SolnBlk.W[i][j] = Wl;
	    SolnBlk.U[i][j] = U(SolnBlk.W[i][j]);
	  }
	}
	break;
      case IC_PERIODIC_SINX_MULTIWAVE :
	Wl.v.x = 100.0;
	Wl.v.y = 0.0;
	Wl.p = PRESSURE_STDATM;
	for (j = SolnBlk.JCl-SolnBlk.Nghost; j <= SolnBlk.JCu+SolnBlk.Nghost; j++) {
	  for (i = SolnBlk.ICl-SolnBlk.Nghost; i <= SolnBlk.ICu+SolnBlk.Nghost; i++) {
	    Wl.d = SolnBlk.Grid.Integration.IntegrateFunctionOverCell(i,j,
								      MultipleSinVariationInXDir,
								      IP.Exact_Integration_Digits,Wl.d)/SolnBlk.Grid.Cell[i][j].A;
	    SolnBlk.W[i][j] = Wl;
	    SolnBlk.U[i][j] = U(SolnBlk.W[i][j]);
	  }
	}
	break;
      case IC_PERIODIC_SINY_MULTIWAVE :
	Wl.v.x = 0.0;
	Wl.v.y = 100.0;
	Wl.p = PRESSURE_STDATM;
	for (j = SolnBlk.JCl-SolnBlk.Nghost; j <= SolnBlk.JCu+SolnBlk.Nghost; j++) {
	  for (i = SolnBlk.ICl-SolnBlk.Nghost; i <= SolnBlk.ICu+SolnBlk.Nghost; i++) {
	    Wl.d = SolnBlk.Grid.Integration.IntegrateFunctionOverCell(i,j,
								      MultipleSinVariationInYDir,
								      IP.Exact_Integration_Digits,Wl.d)/SolnBlk.Grid.Cell[i][j].A;
	    SolnBlk.W[i][j] = Wl;
	    SolnBlk.U[i][j] = U(SolnBlk.W[i][j]);
	  }
	}
	break;
      case IC_PERIODIC_COMPLEX_MULTIWAVE :
	Wl.v.x = 250.0;
	Wl.v.y = 0.0;
	Wl.p = PRESSURE_STDATM;
	for (j = SolnBlk.JCl-SolnBlk.Nghost; j <= SolnBlk.JCu+SolnBlk.Nghost; j++) {
	  for (i = SolnBlk.ICl-SolnBlk.Nghost; i <= SolnBlk.ICu+SolnBlk.Nghost; i++) {
	    Wl.d = SolnBlk.Grid.Integration.IntegrateFunctionOverCell(i,j,
								      Complex_2D_Waves,
								      IP.Exact_Integration_Digits,Wl.d)/SolnBlk.Grid.Cell[i][j].A;
	    SolnBlk.W[i][j] = Wl;
	    SolnBlk.U[i][j] = U(SolnBlk.W[i][j]);
	  }
	}
	break;
      case IC_EXACT_SOLUTION :
	// Set the solution state by calculating the cell average values with integration of the exact solution
	// Use the ExactSoln pointer to access the exact solution
	if (IP.ExactSoln->IsExactSolutionSet()) {
	  for (j  = SolnBlk.JCl-SolnBlk.Nghost ; j <= SolnBlk.JCu+SolnBlk.Nghost ; ++j ) {
	    for ( i = SolnBlk.ICl-SolnBlk.Nghost ; i <= SolnBlk.ICu+SolnBlk.Nghost ; ++i ) {
	      // Compute the exact average value for each solution parameter.
	      for ( k = 1; k <= SolnBlk.NumVar(); ++k){
		SolnBlk.W[i][j][k] = 
		  SolnBlk.Grid.Integration.
		  IntegrateFunctionOverCell(i,j,
					    wrapped_member_function_one_parameter(IP.ExactSoln,
										  &Euler2D_ExactSolutions::SolutionForParameter,
										  k,
										  Wl[k]),
					    wrapped_member_function_one_parameter(IP.ExactSoln,
										  &Euler2D_ExactSolutions::
										  XDependencyIntegrated_SolutionForParameter,
										  k,
										  Wl[k]),
					    IP.Exact_Integration_Digits,Wl[k])/SolnBlk.Grid.Cell[i][j].A;
	      }
	      SolnBlk.U[i][j] = U(SolnBlk.W[i][j]);
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
		// Compute the exact average value for each solution parameter.
		for ( k = 1; k <= SolnBlk.NumVar(); ++k){
		  SolnBlk.W[i][j][k] = 
		    SolnBlk.Grid.Integration.
		    IntegrateFunctionOverCell(i,j,
					      wrapped_member_function_one_parameter(IP.ExactSoln,
										    &Euler2D_ExactSolutions::SolutionForParameter,
										    k,
										    Wl[k]),
					      wrapped_member_function_one_parameter(IP.ExactSoln,
										    &Euler2D_ExactSolutions::
										    XDependencyIntegrated_SolutionForParameter,
										    k,
										    Wl[k]),
					      IP.Exact_Integration_Digits,Wl[k])/SolnBlk.Grid.Cell[i][j].A;
		}
		SolnBlk.U[i][j] = U(SolnBlk.W[i][j]);
	      } else {
		// Set the solution state of the interior cells to the initial state Uo[0].
		SolnBlk.W[i][j] = Wo[0];
		SolnBlk.U[i][j] = U(SolnBlk.W[i][j]);
	      }
	    } /* endfor */
	  } /* endfor */
	  
	} else {
	  // There is no exact solution set for this problem
	  throw runtime_error("ICs() ERROR! No exact solution has been set!");
	}
	break;

      default:
        // Set the solution state to the initial state Wo[0].
        for (j  = SolnBlk.JCl-SolnBlk.Nghost ; j <= SolnBlk.JCu+SolnBlk.Nghost ; ++j ) {
            for ( i = SolnBlk.ICl-SolnBlk.Nghost ; i <= SolnBlk.ICu+SolnBlk.Nghost ; ++i ) {
	       SolnBlk.W[i][j] = Wo[0];
               SolnBlk.U[i][j] = U(SolnBlk.W[i][j]);
            } /* endfor */
        } /* endfor */
        break;
    } /* endswitch */

    /* Set default values for the boundary condition reference states. */

    for (j = SolnBlk.JCl-SolnBlk.Nghost ; j <= SolnBlk.JCu+SolnBlk.Nghost ; ++j ) {
       if (j >= SolnBlk.JCl && j <= SolnBlk.JCu) {
          SolnBlk.WoW[j] = SolnBlk.W[SolnBlk.ICl][j];
          SolnBlk.WoE[j] = SolnBlk.W[SolnBlk.ICu][j];
       } else if (j < SolnBlk.JCl) {
          SolnBlk.WoW[j] = SolnBlk.W[SolnBlk.ICl][SolnBlk.JCl];
          SolnBlk.WoE[j] = SolnBlk.W[SolnBlk.ICu][SolnBlk.JCl];
       } else {
          SolnBlk.WoW[j] = SolnBlk.W[SolnBlk.ICl][SolnBlk.JCu];
          SolnBlk.WoE[j] = SolnBlk.W[SolnBlk.ICu][SolnBlk.JCu];
       } /* endif */
    } /* endfor */

    for ( i = SolnBlk.ICl-SolnBlk.Nghost ; i <= SolnBlk.ICu+SolnBlk.Nghost ; ++i ) {
       if (i >= SolnBlk.ICl && i <= SolnBlk.ICu) {
          SolnBlk.WoS[i] = SolnBlk.W[i][SolnBlk.JCl];
          SolnBlk.WoN[i] = SolnBlk.W[i][SolnBlk.JCu];
       } else if (i < SolnBlk.ICl) {
          SolnBlk.WoS[i] = SolnBlk.W[SolnBlk.ICl][SolnBlk.JCl];
          SolnBlk.WoN[i] = SolnBlk.W[SolnBlk.ICl][SolnBlk.JCu];
       } else {
          SolnBlk.WoS[i] = SolnBlk.W[SolnBlk.ICu][SolnBlk.JCl];
          SolnBlk.WoN[i] = SolnBlk.W[SolnBlk.ICu][SolnBlk.JCu];
       } /* endif */
    } /* endfor */

}

/******************************************************//**
 * Routine: BCs                                         
 *                                                      
 * Apply boundary conditions at boundaries of the       
 * specified quadrilateral solution block.              
 *                                                      
 * \note Implement Ringleb BC properly!
 ********************************************************/
void BCs(Euler2D_Quad_Block &SolnBlk,
	 Euler2D_Input_Parameters &IP) {

  int i, j;
  int ghost;
  Vector2D dX;
  Euler2D_pState dW, W;

  for ( j = SolnBlk.JCl-SolnBlk.Nghost ; j <= SolnBlk.JCu+SolnBlk.Nghost ; ++j ) {
    // Prescribe West boundary conditions.
    if ( (j >= SolnBlk.JCl && j <= SolnBlk.JCu) ||                                  // <-- affects W boundary cells
	 (j < SolnBlk.JCl && (SolnBlk.Grid.BCtypeS[SolnBlk.ICl-1] == BC_NONE ) ) || // <-- affects SW corner cells
	 (j > SolnBlk.JCu && (SolnBlk.Grid.BCtypeN[SolnBlk.ICl-1] == BC_NONE ) ) )  // <-- affects NW corner cells
      {
        switch(SolnBlk.Grid.BCtypeW[j]) {
	case BC_NONE :
	  break;
	case BC_FIXED :
	  for(ghost = 1; ghost <= SolnBlk.Nghost; ++ghost){
	    SolnBlk.W[SolnBlk.ICl-ghost][j] = SolnBlk.WoW[j]; 
	    SolnBlk.U[SolnBlk.ICl-ghost][j] = U(SolnBlk.W[SolnBlk.ICl-ghost][j]);
	  }
	  break;
	case BC_CONSTANT_EXTRAPOLATION :
	  for(ghost = 1; ghost <= SolnBlk.Nghost; ++ghost){
	    SolnBlk.W[SolnBlk.ICl-ghost][j] = SolnBlk.W[SolnBlk.ICl][j];
	    SolnBlk.U[SolnBlk.ICl-ghost][j] = SolnBlk.U[SolnBlk.ICl][j];
	  }
	  break;
	case BC_LINEAR_EXTRAPOLATION :
	  Linear_Reconstruction_LeastSquares_2(SolnBlk, 
					       SolnBlk.ICl, j, 
					       LIMITER_BARTH_JESPERSEN);
	  for(ghost = 1; ghost <= SolnBlk.Nghost; ++ghost){
	    SolnBlk.W[SolnBlk.ICl-ghost][j] = SolnBlk.PiecewiseLinearSolutionAtLocation(SolnBlk.ICl,j,
											SolnBlk.Grid.Cell[SolnBlk.ICl-ghost][j].Xc);
	    SolnBlk.U[SolnBlk.ICl-ghost][j] = U(SolnBlk.W[SolnBlk.ICl-ghost][j]);
	  }	      
	  break;
	case BC_REFLECTION :
	  for(ghost = 1; ghost <= SolnBlk.Nghost; ++ghost){
	    SolnBlk.W[SolnBlk.ICl-ghost][j] = Reflect(SolnBlk.W[SolnBlk.ICl + ghost-1][j],
						      SolnBlk.Grid.nfaceW(SolnBlk.ICl,j));
	    SolnBlk.U[SolnBlk.ICl-ghost][j] = U(SolnBlk.W[SolnBlk.ICl-ghost][j]);
	  }
	  break;
	case BC_BURNING_SURFACE :
	  Linear_Reconstruction_LeastSquares_2(SolnBlk, 
					       SolnBlk.ICl, j, 
					       LIMITER_BARTH_JESPERSEN);
	  dW = SolnBlk.PiecewiseLinearSolutionAtLocation(SolnBlk.ICl, j,
							 SolnBlk.Grid.Cell[SolnBlk.ICl-1][j].Xc);
	  SolnBlk.W[SolnBlk.ICl-1][j] = BurningSurface(dW,
						       SolnBlk.Grid.nfaceW(SolnBlk.ICl,j));
	  SolnBlk.U[SolnBlk.ICl-1][j] = U(SolnBlk.W[SolnBlk.ICl-1][j]);
	  for(ghost = 2; ghost <= SolnBlk.Nghost; ++ghost){
	    SolnBlk.W[SolnBlk.ICl-ghost][j] = SolnBlk.W[SolnBlk.ICu-1][j];
	    SolnBlk.U[SolnBlk.ICl-ghost][j] = SolnBlk.U[SolnBlk.ICu-1][j];
	  }
	  break;
	case BC_MASS_INJECTION :  
	  // Invalid for arbitrary number of ghost cells
	  SolnBlk.W[SolnBlk.ICl-1][j] = MassInjection(SolnBlk.W[SolnBlk.ICl][j],SolnBlk.Grid.nfaceW(SolnBlk.ICl,j),ON);
	  SolnBlk.U[SolnBlk.ICl-1][j] = U(SolnBlk.W[SolnBlk.ICl-1][j]);
	  SolnBlk.W[SolnBlk.ICl-2][j] = SolnBlk.W[SolnBlk.ICl-1][j];
	  SolnBlk.U[SolnBlk.ICl-2][j] = SolnBlk.U[SolnBlk.ICl-1][j]; 
	  break;
	case BC_PERIODIC :
	  for(ghost = 1; ghost <= SolnBlk.Nghost; ++ghost){
	    SolnBlk.W[SolnBlk.ICl-ghost][j] = SolnBlk.W[SolnBlk.ICu-ghost+1][j];
	    SolnBlk.U[SolnBlk.ICl-ghost][j] = SolnBlk.U[SolnBlk.ICu-ghost+1][j];
	  }
	  break;
	case BC_CHARACTERISTIC :
	  SolnBlk.W[SolnBlk.ICl-1][j] = 
	    BC_Characteristic_Pressure(SolnBlk.W[SolnBlk.ICl][j],
				       SolnBlk.WoW[j], 
				       SolnBlk.Grid.nfaceW(SolnBlk.ICl,j));
	  SolnBlk.U[SolnBlk.ICl-1][j] = U(SolnBlk.W[SolnBlk.ICl-1][j]);
	  for(ghost = 2; ghost <= SolnBlk.Nghost; ++ghost){
	    SolnBlk.W[SolnBlk.ICl-ghost][j] = SolnBlk.W[SolnBlk.ICl-1][j];
	    SolnBlk.U[SolnBlk.ICl-ghost][j] = SolnBlk.U[SolnBlk.ICl-1][j];
	  }
	  break;
	case BC_RINGLEB_FLOW :
	  // Not properly implemented yet!

	  // 	    RinglebFlowFunctionType RinglebFlowFunction = RinglebFlow;
	    
	  // 	    for(ghost = 1; ghost <= SolnBlk.Nghost; ++ghost){ 
	  // 	      //  SolnBlk.W[SolnBlk.ICl-ghost][j] = RinglebFlow(SolnBlk.Grid.Cell[SolnBlk.ICl-ghost][j].Xc);
	    
	  // 	      SolnBlk.W[SolnBlk.ICl-ghost][j]=(SolnBlk.IntegrateOverTheCell(SolnBlk.ICl-ghost,j,
	  // 									    RinglebFlowFunction,8,
	  // 									    SolnBlk.W[SolnBlk.ICl-ghost][j])/
	  // 					       SolnBlk.Grid.area(SolnBlk.ICl-ghost,j));
	  // 	      // 	      SolnBlk.W[SolnBlk.ICl-ghost][j] = Reflect(SolnBlk.W[SolnBlk.ICl+ghost-1][j],
	  // 	      // 							SolnBlk.Grid.nfaceW(SolnBlk.ICl,j));
	  // 	      SolnBlk.U[SolnBlk.ICl-ghost][j] = U(SolnBlk.W[SolnBlk.ICl-ghost][j]);
	  // 	    }

	  //             SolnBlk.W[SolnBlk.ICl-1][j] = RinglebFlow(SolnBlk.W[SolnBlk.ICl-1][j],
	  // 						      SolnBlk.Grid.Cell[SolnBlk.ICl-1][j].Xc);
	  //             SolnBlk.U[SolnBlk.ICl-1][j] = U(SolnBlk.W[SolnBlk.ICl-1][j]);
	  //             SolnBlk.W[SolnBlk.ICl-2][j] = RinglebFlow(SolnBlk.W[SolnBlk.ICl-2][j],
	  // 						      SolnBlk.Grid.Cell[SolnBlk.ICl-2][j].Xc);
	  //             SolnBlk.U[SolnBlk.ICl-2][j] = U(SolnBlk.W[SolnBlk.ICl-2][j]);
	    
	  SolnBlk.W[SolnBlk.ICl-1][j] = Reflect(SolnBlk.W[SolnBlk.ICl][j],
						SolnBlk.Grid.nfaceW(SolnBlk.ICl,j));
	  SolnBlk.U[SolnBlk.ICl-1][j] = U(SolnBlk.W[SolnBlk.ICl-1][j]);
	  SolnBlk.W[SolnBlk.ICl-2][j] = Reflect(SolnBlk.W[SolnBlk.ICl+1][j],
						SolnBlk.Grid.nfaceW(SolnBlk.ICl,j));
	  SolnBlk.U[SolnBlk.ICl-2][j] = U(SolnBlk.W[SolnBlk.ICl-2][j]);
	  break;
	case BC_EXACT_SOLUTION:
	  // Leave the ghost cell values unchanged.
	  break;
	case BC_FROZEN :
	  // Equivalent to BC_NONE in the sense that it leaves 
	  // the ghost cell values unchanged and 
	  // uses the ghost cell reconstruction to 
	  // compute the inter-cellular state.
	  break;
	case BC_CHARACTERISTIC_VELOCITY:
	  SolnBlk.W[SolnBlk.ICl-1][j] = BC_Characteristic(SolnBlk.W[SolnBlk.ICl][j],
							  SolnBlk.WoW[j], 
							  SolnBlk.Grid.nfaceW(SolnBlk.ICl,j));
	  SolnBlk.U[SolnBlk.ICl-1][j] = U(SolnBlk.W[SolnBlk.ICl-1][j]);
	  for(ghost = 2; ghost <= SolnBlk.Nghost; ++ghost){
	    SolnBlk.W[SolnBlk.ICl-ghost][j] = SolnBlk.W[SolnBlk.ICl-1][j];
	    SolnBlk.U[SolnBlk.ICl-ghost][j] = SolnBlk.U[SolnBlk.ICl-1][j];
	  }
	  break;
	case BC_WALL_INVISCID:
	  // Implement inviscid wall boundary condition 
	  // (i.e. no mass flow through the wall)

	  Linear_Reconstruction_LeastSquares_2(SolnBlk, 
					       SolnBlk.ICl, j, 
					       LIMITER_BARTH_JESPERSEN);

	  for(ghost = 1; ghost <= SolnBlk.Nghost; ++ghost){
	    // Implement reflection of velocities
	    SolnBlk.W[SolnBlk.ICl-ghost][j] = Reflect(SolnBlk.W[SolnBlk.ICl + ghost-1][j],
						      SolnBlk.Grid.nfaceW(SolnBlk.ICl,j));
	    // Implement pressure and density based on linear extrapolation
	    W = SolnBlk.PiecewiseLinearSolutionAtLocation(SolnBlk.ICl,j,
							  SolnBlk.Grid.Cell[SolnBlk.ICl-ghost][j].Xc);
	    SolnBlk.W[SolnBlk.ICl-ghost][j].p = W.p;
	    SolnBlk.W[SolnBlk.ICl-ghost][j].d = W.d;
	    // Calculate the conserved variables
	    SolnBlk.U[SolnBlk.ICl-ghost][j] = U(SolnBlk.W[SolnBlk.ICl-ghost][j]);
	  }
	  break;	  

	default:
	  // Impose constant extrapolation
	  for(ghost = 1; ghost <= SolnBlk.Nghost; ++ghost){
	    SolnBlk.W[SolnBlk.ICl-ghost][j] = SolnBlk.W[SolnBlk.ICl][j];
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
	switch(SolnBlk.Grid.BCtypeE[j]) 
	  {
	  case BC_NONE :
	    break;
	  case BC_FIXED :
	    for(ghost = 1; ghost <= SolnBlk.Nghost; ++ghost){
	      SolnBlk.W[SolnBlk.ICu+ghost][j] = SolnBlk.WoE[j];
	      SolnBlk.U[SolnBlk.ICu+ghost][j] = U(SolnBlk.W[SolnBlk.ICu+ghost][j]);
	    }
	    break;
	  case BC_CONSTANT_EXTRAPOLATION :
	    for(ghost = 1; ghost <= SolnBlk.Nghost; ++ghost){
	      SolnBlk.W[SolnBlk.ICu+ghost][j] = SolnBlk.W[SolnBlk.ICu][j];
	      SolnBlk.U[SolnBlk.ICu+ghost][j] = SolnBlk.U[SolnBlk.ICu][j];
	    }
	    break;
	  case BC_LINEAR_EXTRAPOLATION :
	    Linear_Reconstruction_LeastSquares_2(SolnBlk, 
						 SolnBlk.ICu, j, 
						 LIMITER_BARTH_JESPERSEN);
	    for(ghost = 1; ghost <= SolnBlk.Nghost; ++ghost){
	      SolnBlk.W[SolnBlk.ICu+ghost][j] = SolnBlk.PiecewiseLinearSolutionAtLocation(SolnBlk.ICu, j,
											  SolnBlk.Grid.Cell[SolnBlk.ICu+ghost][j].Xc);
	      SolnBlk.U[SolnBlk.ICu+ghost][j] = U(SolnBlk.W[SolnBlk.ICu+ghost][j]);
	    }    
	    break;
	  case BC_REFLECTION :
	    for(ghost = 1; ghost <= SolnBlk.Nghost; ++ghost){
	      SolnBlk.W[SolnBlk.ICu+ghost][j] = Reflect(SolnBlk.W[SolnBlk.ICu-ghost+1][j],
							SolnBlk.Grid.nfaceE(SolnBlk.ICu,j));
	      SolnBlk.U[SolnBlk.ICu+ghost][j] = U(SolnBlk.W[SolnBlk.ICu+ghost][j]);
	    }       
	    break;
	  case BC_BURNING_SURFACE :
	    Linear_Reconstruction_LeastSquares_2(SolnBlk, 
						 SolnBlk.ICu, j, 
						 LIMITER_BARTH_JESPERSEN);
	    dW = SolnBlk.PiecewiseLinearSolutionAtLocation(SolnBlk.ICu, j,
							   SolnBlk.Grid.Cell[SolnBlk.ICu+1][j].Xc);
	    SolnBlk.W[SolnBlk.ICu+1][j] = BurningSurface(dW,
							 SolnBlk.Grid.nfaceE(SolnBlk.ICu,j));
	    SolnBlk.U[SolnBlk.ICu+1][j] = U(SolnBlk.W[SolnBlk.ICu+1][j]);
	    for(ghost = 2; ghost <= SolnBlk.Nghost; ++ghost){
	      SolnBlk.W[SolnBlk.ICl+ghost][j] = SolnBlk.W[SolnBlk.ICu+1][j];
	      SolnBlk.U[SolnBlk.ICl+ghost][j] = SolnBlk.U[SolnBlk.ICu+1][j];
	    }
	    break;
	  case BC_MASS_INJECTION :
	    // Invalid for arbitrary number of ghost cells
	    SolnBlk.W[SolnBlk.ICu+1][j] = MassInjection(SolnBlk.W[SolnBlk.ICu][j],SolnBlk.Grid.nfaceE(SolnBlk.ICu,j),ON);
	    SolnBlk.U[SolnBlk.ICu+1][j] = SolnBlk.U[SolnBlk.ICl+1][j];
	    SolnBlk.W[SolnBlk.ICu+2][j] = SolnBlk.W[SolnBlk.ICl+2][j];
	    SolnBlk.U[SolnBlk.ICu+2][j] = SolnBlk.U[SolnBlk.ICl+2][j];
	    break;
	  case BC_PERIODIC :
	    for(ghost = 1; ghost <= SolnBlk.Nghost; ++ghost){
	      SolnBlk.W[SolnBlk.ICu+ghost][j] = SolnBlk.W[SolnBlk.ICl+ghost-1][j];
	      SolnBlk.U[SolnBlk.ICu+ghost][j] = SolnBlk.U[SolnBlk.ICl+ghost-1][j];
	    }
	    break;
	  case BC_CHARACTERISTIC :
	    SolnBlk.W[SolnBlk.ICu+1][j] = 
	      BC_Characteristic_Pressure(SolnBlk.W[SolnBlk.ICu][j],
					 SolnBlk.WoE[j],
					 SolnBlk.Grid.nfaceE(SolnBlk.ICu,j));
	    SolnBlk.U[SolnBlk.ICu+1][j] = U(SolnBlk.W[SolnBlk.ICu+1][j]);
	    
	    for(ghost = 2; ghost <= SolnBlk.Nghost; ++ghost){
	      SolnBlk.W[SolnBlk.ICu+ghost][j] = SolnBlk.W[SolnBlk.ICu+1][j];
	      SolnBlk.U[SolnBlk.ICu+ghost][j] = SolnBlk.U[SolnBlk.ICu+1][j];
	    }
	    break;
	  case BC_RINGLEB_FLOW :
	    // Not properly implemented yet!
	    
	    // 	    RinglebFlowFunctionType RinglebFlowFunction = RinglebFlow;
	    // 	    for(ghost = 1; ghost <= SolnBlk.Nghost; ++ghost){
	    // 	      //             SolnBlk.W[SolnBlk.ICu+ghost][j] = RinglebFlow(SolnBlk.Grid.Cell[SolnBlk.ICu+ghost][j].Xc);
	    
	    // 	      SolnBlk.W[SolnBlk.ICu+ghost][j] = (SolnBlk.IntegrateOverTheCell(SolnBlk.ICu+ghost,j,
	    // 									     RinglebFlowFunction,8,
	    // 									     SolnBlk.W[SolnBlk.ICu+ghost][j])/
	    // 						 SolnBlk.Grid.area(SolnBlk.ICu+ghost,j));
	    
	    // 	      // 	      SolnBlk.W[SolnBlk.ICu+ghost][j] = Reflect(SolnBlk.W[SolnBlk.ICu-ghost+1][j],
	    // 	      // 						    SolnBlk.Grid.nfaceE(SolnBlk.ICu,j));
	    // 	      SolnBlk.U[SolnBlk.ICu+ghost][j] = U(SolnBlk.W[SolnBlk.ICu+ghost][j]);
	    // 	    }
	  

	    //  SolnBlk.W[SolnBlk.ICu+1][j] = RinglebFlow(SolnBlk.W[SolnBlk.ICu+1][j],SolnBlk.Grid.Cell[SolnBlk.ICu+1][j].Xc);
	    //  SolnBlk.U[SolnBlk.ICu+1][j] = U(SolnBlk.W[SolnBlk.ICu+1][j]);
	    //  SolnBlk.W[SolnBlk.ICu+2][j] = RinglebFlow(SolnBlk.W[SolnBlk.ICu+2][j],SolnBlk.Grid.Cell[SolnBlk.ICu+2][j].Xc);
	    //  SolnBlk.U[SolnBlk.ICu+2][j] = U(SolnBlk.W[SolnBlk.ICu+2][j]);
	    SolnBlk.W[SolnBlk.ICu+1][j] = Reflect(SolnBlk.W[SolnBlk.ICu][j],
						  SolnBlk.Grid.nfaceE(SolnBlk.ICu,j));
	    SolnBlk.U[SolnBlk.ICu+1][j] = U(SolnBlk.W[SolnBlk.ICu+1][j]);
	    SolnBlk.W[SolnBlk.ICu+2][j] = Reflect(SolnBlk.W[SolnBlk.ICu-1][j],
						  SolnBlk.Grid.nfaceE(SolnBlk.ICu,j));
	    SolnBlk.U[SolnBlk.ICu+2][j] = U(SolnBlk.W[SolnBlk.ICu+2][j]);
	    break;

	  case BC_EXACT_SOLUTION:
	    break;
	  case BC_FROZEN :
	    // Equivalent to BC_NONE in the sense that it leaves 
	    // the ghost cell values unchanged and 
	    // uses the ghost cell reconstruction to 
	    // compute the inter-cellular state.
	    break;
	  case BC_CHARACTERISTIC_VELOCITY:
	    SolnBlk.W[SolnBlk.ICu+1][j] = BC_Characteristic(SolnBlk.W[SolnBlk.ICu][j],
							    SolnBlk.WoE[j],
							    SolnBlk.Grid.nfaceE(SolnBlk.ICu,j));
	    SolnBlk.U[SolnBlk.ICu+1][j] = U(SolnBlk.W[SolnBlk.ICu+1][j]);
	    
	    for(ghost = 2; ghost <= SolnBlk.Nghost; ++ghost){
	      SolnBlk.W[SolnBlk.ICu+ghost][j] = SolnBlk.W[SolnBlk.ICu+1][j];
	      SolnBlk.U[SolnBlk.ICu+ghost][j] = SolnBlk.U[SolnBlk.ICu+1][j];
	    }
	    break;
	  case BC_WALL_INVISCID:
	    // Implement inviscid wall boundary condition 
	    // (i.e. no mass flow through the wall)

	    Linear_Reconstruction_LeastSquares_2(SolnBlk, 
						 SolnBlk.ICu, j, 
						 LIMITER_BARTH_JESPERSEN);

	    for(ghost = 1; ghost <= SolnBlk.Nghost; ++ghost){
	      // Implement reflection of velocities
	      SolnBlk.W[SolnBlk.ICu+ghost][j] = Reflect(SolnBlk.W[SolnBlk.ICu - ghost+1][j],
							SolnBlk.Grid.nfaceE(SolnBlk.ICu,j));
	      // Implement pressure and density based on linear extrapolation
	      W = SolnBlk.PiecewiseLinearSolutionAtLocation(SolnBlk.ICu,j,
							    SolnBlk.Grid.Cell[SolnBlk.ICu+ghost][j].Xc);
	      SolnBlk.W[SolnBlk.ICu+ghost][j].p = W.p;
	      SolnBlk.W[SolnBlk.ICu+ghost][j].d = W.d;
	      // Calculate the conserved variables
	      SolnBlk.U[SolnBlk.ICu+ghost][j] = U(SolnBlk.W[SolnBlk.ICu+ghost][j]);
	    }
	    break;	  

	  default:
	    for(ghost = 1; ghost <= SolnBlk.Nghost; ++ghost){
	      SolnBlk.W[SolnBlk.ICu+ghost][j] = SolnBlk.W[SolnBlk.ICu][j];
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
    case BC_NONE :
      break;
    case BC_FIXED :
      for(ghost = 1; ghost <= SolnBlk.Nghost; ++ghost){
	SolnBlk.W[i][SolnBlk.JCl-ghost] = SolnBlk.WoS[i];
	SolnBlk.U[i][SolnBlk.JCl-ghost] = U(SolnBlk.W[i][SolnBlk.JCl-ghost]);
      }
      break;
    case BC_CONSTANT_EXTRAPOLATION :
      for(ghost = 1; ghost <= SolnBlk.Nghost; ++ghost){
	SolnBlk.W[i][SolnBlk.JCl-ghost] = SolnBlk.W[i][SolnBlk.JCl];
	SolnBlk.U[i][SolnBlk.JCl-ghost] = SolnBlk.U[i][SolnBlk.JCl];
      }
      break;
    case BC_LINEAR_EXTRAPOLATION :
      Linear_Reconstruction_LeastSquares_2(SolnBlk, 
					   i, SolnBlk.JCl, 
					   LIMITER_BARTH_JESPERSEN);
      for(ghost = 1; ghost <= SolnBlk.Nghost; ++ghost){
	SolnBlk.W[i][SolnBlk.JCl-ghost] = SolnBlk.PiecewiseLinearSolutionAtLocation(i,SolnBlk.JCl,
										    SolnBlk.Grid.Cell[i][SolnBlk.JCl-ghost].Xc);
	SolnBlk.U[i][SolnBlk.JCl-ghost] = U(SolnBlk.W[i][SolnBlk.JCl-ghost]);
      }
      break;
    case BC_REFLECTION :
      for(ghost = 1; ghost <= SolnBlk.Nghost; ++ghost){
	SolnBlk.W[i][SolnBlk.JCl-ghost] = Reflect(SolnBlk.W[i][SolnBlk.JCl + ghost-1],
						  SolnBlk.Grid.nfaceS(i, SolnBlk.JCl));
	SolnBlk.U[i][SolnBlk.JCl-ghost] = U(SolnBlk.W[i][SolnBlk.JCl-ghost]);    
      }
      break;
    case BC_BURNING_SURFACE :
      Linear_Reconstruction_LeastSquares_2(SolnBlk, 
					   i, SolnBlk.JCl, 
					   LIMITER_BARTH_JESPERSEN);
      dW = SolnBlk.PiecewiseLinearSolutionAtLocation(i,SolnBlk.JCl,
						     SolnBlk.Grid.Cell[i][SolnBlk.JCl-1].Xc);
      SolnBlk.W[i][SolnBlk.JCl-1] = BurningSurface(dW, SolnBlk.Grid.nfaceS(i, SolnBlk.JCl));
      SolnBlk.U[i][SolnBlk.JCl-1] = U(SolnBlk.W[i][SolnBlk.JCl-1]);
      for(ghost = 2; ghost <= SolnBlk.Nghost; ++ghost){
	SolnBlk.W[i][SolnBlk.JCl-ghost] = SolnBlk.W[i][SolnBlk.JCl-1];
	SolnBlk.U[i][SolnBlk.JCl-ghost] = SolnBlk.U[i][SolnBlk.JCl-1];
      }
      break;
    case BC_MASS_INJECTION :
      // Invalid for arbitrary number of ghost cells
      SolnBlk.W[i][SolnBlk.JCl-1] = MassInjection(SolnBlk.W[i][SolnBlk.JCl],SolnBlk.Grid.nfaceS(i,SolnBlk.JCl),ON);
      SolnBlk.U[i][SolnBlk.JCl-1] = SolnBlk.U[i][SolnBlk.JCu-1];
      SolnBlk.W[i][SolnBlk.JCl-2] = SolnBlk.W[i][SolnBlk.JCu-2];
      SolnBlk.U[i][SolnBlk.JCl-2] = SolnBlk.U[i][SolnBlk.JCu-2];
      break;
    case BC_PERIODIC :
      for(ghost = 1; ghost <= SolnBlk.Nghost; ++ghost){
	SolnBlk.W[i][SolnBlk.JCl-ghost] = SolnBlk.W[i][SolnBlk.JCu-ghost+1];
	SolnBlk.U[i][SolnBlk.JCl-ghost] = SolnBlk.U[i][SolnBlk.JCu-ghost+1];
      }
      break;
    case BC_CHARACTERISTIC :
      SolnBlk.W[i][SolnBlk.JCl-1] = BC_Characteristic_Pressure(SolnBlk.W[i][SolnBlk.JCl],
							       SolnBlk.WoS[i],
							       SolnBlk.Grid.nfaceS(i, SolnBlk.JCl));
      SolnBlk.U[i][SolnBlk.JCl-1] = U(SolnBlk.W[i][SolnBlk.JCl-1]); 
      for(ghost = 2; ghost <= SolnBlk.Nghost; ++ghost){
	SolnBlk.W[i][SolnBlk.JCl-ghost] = SolnBlk.W[i][SolnBlk.JCl-1];
	SolnBlk.U[i][SolnBlk.JCl-ghost] = SolnBlk.U[i][SolnBlk.JCl-1];
      }
      break;
    case BC_RINGLEB_FLOW :
      // Not properly implemented yet!

      // 	    RinglebFlowFunctionType RinglebFlowFunction = RinglebFlow;
      // 	    for(ghost = 1; ghost <= SolnBlk.Nghost; ++ghost){
      // 	      //  SolnBlk.W[i][SolnBlk.JCl-ghost] = RinglebFlow(SolnBlk.Grid.Cell[i][SolnBlk.JCl-ghost].Xc);
      
      // 	      SolnBlk.W[i][SolnBlk.JCl-ghost] = (SolnBlk.IntegrateOverTheCell(i,SolnBlk.JCl-ghost,
      // 									     RinglebFlowFunction,8,
      // 									     SolnBlk.W[i][SolnBlk.JCl-ghost])/
      // 						 SolnBlk.Grid.area(i,SolnBlk.JCl-ghost));
      // 	      SolnBlk.U[i][SolnBlk.JCl-ghost] = U(SolnBlk.W[i][SolnBlk.JCl-ghost]);
      // 	    }

      SolnBlk.W[i][SolnBlk.JCl-1] = RinglebFlow(SolnBlk.W[i][SolnBlk.JCl-1],SolnBlk.Grid.Cell[i][SolnBlk.JCl-1].Xc);
      SolnBlk.U[i][SolnBlk.JCl-1] = U(SolnBlk.W[i][SolnBlk.JCl-1]);
      SolnBlk.W[i][SolnBlk.JCl-2] = RinglebFlow(SolnBlk.W[i][SolnBlk.JCl-2],SolnBlk.Grid.Cell[i][SolnBlk.JCl-2].Xc);
      SolnBlk.U[i][SolnBlk.JCl-2] = U(SolnBlk.W[i][SolnBlk.JCl-2]);
      //
      //             SolnBlk.W[i][SolnBlk.JCl-1] = SolnBlk.W[i][SolnBlk.JCl];
      //             SolnBlk.U[i][SolnBlk.JCl-1] = SolnBlk.U[i][SolnBlk.JCl];
      //             SolnBlk.W[i][SolnBlk.JCl-2] = SolnBlk.W[i][SolnBlk.JCl];
      //             SolnBlk.U[i][SolnBlk.JCl-2] = SolnBlk.U[i][SolnBlk.JCl];
      //
      //             Linear_Reconstruction_LeastSquares_2(SolnBlk, 
      //                                                  i, SolnBlk.JCl, 
      //                                                  LIMITER_BARTH_JESPERSEN);
      //             dX = SolnBlk.Grid.Cell[i][SolnBlk.JCl-1].Xc -
      //                  SolnBlk.Grid.Cell[i][SolnBlk.JCl].Xc;
      //             SolnBlk.W[i][SolnBlk.JCl-1] = SolnBlk.W[i][SolnBlk.JCl] + 
      //                (SolnBlk.phi[i][SolnBlk.JCl]^SolnBlk.dWdx[i][SolnBlk.JCl])*dX.x +
      //                (SolnBlk.phi[i][SolnBlk.JCl]^SolnBlk.dWdy[i][SolnBlk.JCl])*dX.y;
      //             SolnBlk.U[i][SolnBlk.JCl-1] = U(SolnBlk.W[i][SolnBlk.JCl-1]);
      //             dX = SolnBlk.Grid.Cell[i][SolnBlk.JCl-2].Xc -
      //                  SolnBlk.Grid.Cell[i][SolnBlk.JCl].Xc;
      //             SolnBlk.W[i][SolnBlk.JCl-2] = SolnBlk.W[i][SolnBlk.JCl] + 
      //                (SolnBlk.phi[i][SolnBlk.JCl]^SolnBlk.dWdx[i][SolnBlk.JCl])*dX.x +
      //                (SolnBlk.phi[i][SolnBlk.JCl]^SolnBlk.dWdy[i][SolnBlk.JCl])*dX.y;
      //             SolnBlk.U[i][SolnBlk.JCl-2] = U(SolnBlk.W[i][SolnBlk.JCl-2]);
      break;
    case BC_EXACT_SOLUTION:
      break;
    case BC_FROZEN :
      // Equivalent to BC_NONE in the sense that it leaves 
      // the ghost cell values unchanged and 
      // uses the ghost cell reconstruction to 
      // compute the inter-cellular state.
      break;
    case BC_CHARACTERISTIC_VELOCITY:
      SolnBlk.W[i][SolnBlk.JCl-1] = BC_Characteristic(SolnBlk.W[i][SolnBlk.JCl],
						      SolnBlk.WoS[i],
						      SolnBlk.Grid.nfaceS(i, SolnBlk.JCl));
      SolnBlk.U[i][SolnBlk.JCl-1] = U(SolnBlk.W[i][SolnBlk.JCl-1]); 
      for(ghost = 2; ghost <= SolnBlk.Nghost; ++ghost){
	SolnBlk.W[i][SolnBlk.JCl-ghost] = SolnBlk.W[i][SolnBlk.JCl-1];
	SolnBlk.U[i][SolnBlk.JCl-ghost] = SolnBlk.U[i][SolnBlk.JCl-1];
      }
      break;
    case BC_WALL_INVISCID:
      // Implement inviscid wall boundary condition 
      // (i.e. no mass flow through the wall)
      
      Linear_Reconstruction_LeastSquares_2(SolnBlk, 
					   i, SolnBlk.JCl, 
					   LIMITER_BARTH_JESPERSEN);
      
      for(ghost = 1; ghost <= SolnBlk.Nghost; ++ghost){
	// Implement reflection of velocities
	SolnBlk.W[i][SolnBlk.JCl-ghost] = Reflect(SolnBlk.W[i][SolnBlk.JCl + ghost-1],
						  SolnBlk.Grid.nfaceS(i, SolnBlk.JCl));
	// Implement pressure and density based on linear extrapolation
	W = SolnBlk.PiecewiseLinearSolutionAtLocation(i,SolnBlk.JCl,
						      SolnBlk.Grid.Cell[i][SolnBlk.JCl-ghost].Xc);
	SolnBlk.W[i][SolnBlk.JCl-ghost].p = W.p;
	SolnBlk.W[i][SolnBlk.JCl-ghost].d = W.d;
	// Calculate the conserved variables
	SolnBlk.U[i][SolnBlk.JCl-ghost] = U(SolnBlk.W[i][SolnBlk.JCl-ghost]);
      }
      break;	  

    default:
      for(ghost = 1; ghost <= SolnBlk.Nghost; ++ghost){
	SolnBlk.W[i][SolnBlk.JCl-ghost] = SolnBlk.W[i][SolnBlk.JCl];
	SolnBlk.U[i][SolnBlk.JCl-ghost] = SolnBlk.U[i][SolnBlk.JCl];
      }
      break;
    } /* endswitch */

    // Prescribe North boundary conditions.
    switch(SolnBlk.Grid.BCtypeN[i]) {
    case BC_NONE :
      break;
    case BC_FIXED :
      for(ghost = 1; ghost <= SolnBlk.Nghost; ++ghost){
	SolnBlk.W[i][SolnBlk.JCu+ghost] = SolnBlk.WoN[i];
	SolnBlk.U[i][SolnBlk.JCu+ghost] = U(SolnBlk.W[i][SolnBlk.JCu+ghost]);
      }
      break;
    case BC_CONSTANT_EXTRAPOLATION :
      for(ghost = 1; ghost <= SolnBlk.Nghost; ++ghost){
	SolnBlk.W[i][SolnBlk.JCu+ghost] = SolnBlk.W[i][SolnBlk.JCu];
	SolnBlk.U[i][SolnBlk.JCu+ghost] = SolnBlk.U[i][SolnBlk.JCu];
      }
      break;
    case BC_LINEAR_EXTRAPOLATION :
      Linear_Reconstruction_LeastSquares_2(SolnBlk, 
					   i, SolnBlk.JCu, 
					   LIMITER_BARTH_JESPERSEN);
      for(ghost = 1; ghost <= SolnBlk.Nghost; ++ghost){
	SolnBlk.W[i][SolnBlk.JCu+ghost] = SolnBlk.PiecewiseLinearSolutionAtLocation(i,SolnBlk.JCu,
										    SolnBlk.Grid.Cell[i][SolnBlk.JCu+ghost].Xc);
	SolnBlk.U[i][SolnBlk.JCu+ghost] = U(SolnBlk.W[i][SolnBlk.JCu+ghost]);
      }
      break;
    case BC_REFLECTION :
      for(ghost = 1; ghost <= SolnBlk.Nghost; ++ghost){
	SolnBlk.W[i][SolnBlk.JCu+ghost] = Reflect(SolnBlk.W[i][SolnBlk.JCu-ghost+1],
						  SolnBlk.Grid.nfaceN(i, SolnBlk.JCu));
	SolnBlk.U[i][SolnBlk.JCu+ghost] = U(SolnBlk.W[i][SolnBlk.JCu+ghost]);
      }
      break;
    case BC_BURNING_SURFACE :
      Linear_Reconstruction_LeastSquares_2(SolnBlk, 
					   i, SolnBlk.JCu, 
					   LIMITER_BARTH_JESPERSEN);
      dW = SolnBlk.PiecewiseLinearSolutionAtLocation(i, SolnBlk.JCu,
						     SolnBlk.Grid.Cell[i][SolnBlk.JCu+1].Xc);
      SolnBlk.W[i][SolnBlk.JCu+1] = BurningSurface(dW, SolnBlk.Grid.nfaceN(i, SolnBlk.JCu));
      SolnBlk.U[i][SolnBlk.JCu+1] = U(SolnBlk.W[i][SolnBlk.JCu+1]);
      for(ghost = 2; ghost <= SolnBlk.Nghost; ++ghost){
	SolnBlk.W[i][SolnBlk.JCu+ghost] = SolnBlk.W[i][SolnBlk.JCu+1];
	SolnBlk.U[i][SolnBlk.JCu+ghost] = SolnBlk.U[i][SolnBlk.JCu+1];
      }
      break;
    case BC_MASS_INJECTION :
      // Invalid for arbitrary number of ghost cells
      SolnBlk.W[i][SolnBlk.JCu+1] = MassInjection(SolnBlk.W[i][SolnBlk.JCu],SolnBlk.Grid.nfaceN(i,SolnBlk.JCu),ON);
      SolnBlk.U[i][SolnBlk.JCu+1] = U(SolnBlk.W[i][SolnBlk.JCu+1]);
      SolnBlk.W[i][SolnBlk.JCu+2] = SolnBlk.W[i][SolnBlk.JCu+1];
      SolnBlk.U[i][SolnBlk.JCu+2] = SolnBlk.U[i][SolnBlk.JCu+1];
      break;
    case BC_PERIODIC :
      for(ghost = 1; ghost <= SolnBlk.Nghost; ++ghost){
	SolnBlk.W[i][SolnBlk.JCu+ghost] = SolnBlk.W[i][SolnBlk.JCl+ghost-1];
	SolnBlk.U[i][SolnBlk.JCu+ghost] = SolnBlk.U[i][SolnBlk.JCl+ghost-1];
      }
      break;
    case BC_CHARACTERISTIC :
      SolnBlk.W[i][SolnBlk.JCu+1] = BC_Characteristic_Pressure(SolnBlk.W[i][SolnBlk.JCu],
							       SolnBlk.WoN[i],
							       SolnBlk.Grid.nfaceN(i, SolnBlk.JCu));
      SolnBlk.U[i][SolnBlk.JCu+1] = U(SolnBlk.W[i][SolnBlk.JCu+1]);
      for(ghost = 2; ghost <= SolnBlk.Nghost; ++ghost){
	SolnBlk.W[i][SolnBlk.JCu+ghost] = SolnBlk.W[i][SolnBlk.JCu+1];
	SolnBlk.U[i][SolnBlk.JCu+ghost] = SolnBlk.U[i][SolnBlk.JCu+1];
      }
      break;
    case BC_RINGLEB_FLOW :
      // Not properly implemented yet!

      // 	    RinglebFlowFunctionType RinglebFlowFunction = RinglebFlow;
      // 	    for(ghost = 1; ghost <= SolnBlk.Nghost; ++ghost){
      // 	      //	      SolnBlk.W[i][SolnBlk.JCu+ghost] = RinglebFlow(SolnBlk.Grid.Cell[i][SolnBlk.JCu+ghost].Xc);
      // 	      SolnBlk.W[i][SolnBlk.JCu+ghost] = (SolnBlk.IntegrateOverTheCell(i,SolnBlk.JCu+ghost,
      // 									      RinglebFlowFunction,8,
      // 									      SolnBlk.W[i][SolnBlk.JCu+ghost])/
      // 						 SolnBlk.Grid.area(i,SolnBlk.JCu+ghost));
      // 	      SolnBlk.U[i][SolnBlk.JCu+ghost] = U(SolnBlk.W[i][SolnBlk.JCu+ghost]);
      // 	    }

      SolnBlk.W[i][SolnBlk.JCu+1] = RinglebFlow(SolnBlk.W[i][SolnBlk.JCu+1],SolnBlk.Grid.Cell[i][SolnBlk.JCu+1].Xc);
      SolnBlk.U[i][SolnBlk.JCu+1] = U(SolnBlk.W[i][SolnBlk.JCu+1]);
      SolnBlk.W[i][SolnBlk.JCu+2] = RinglebFlow(SolnBlk.W[i][SolnBlk.JCu+2],SolnBlk.Grid.Cell[i][SolnBlk.JCu+2].Xc);
      SolnBlk.U[i][SolnBlk.JCu+2] = U(SolnBlk.W[i][SolnBlk.JCu+2]);
      break;
    case BC_EXACT_SOLUTION:
      break;
    case BC_FROZEN :
      // Equivalent to BC_NONE in the sense that it leaves 
      // the ghost cell values unchanged and 
      // uses the ghost cell reconstruction to 
      // compute the inter-cellular state.
      break;
    case BC_CHARACTERISTIC_VELOCITY:
      SolnBlk.W[i][SolnBlk.JCu+1] = BC_Characteristic(SolnBlk.W[i][SolnBlk.JCu],
						      SolnBlk.WoN[i],
						      SolnBlk.Grid.nfaceN(i, SolnBlk.JCu));
      SolnBlk.U[i][SolnBlk.JCu+1] = U(SolnBlk.W[i][SolnBlk.JCu+1]);
      for(ghost = 2; ghost <= SolnBlk.Nghost; ++ghost){
	SolnBlk.W[i][SolnBlk.JCu+ghost] = SolnBlk.W[i][SolnBlk.JCu+1];
	SolnBlk.U[i][SolnBlk.JCu+ghost] = SolnBlk.U[i][SolnBlk.JCu+1];
      }
      break;
    case BC_WALL_INVISCID:
      // Implement inviscid wall boundary condition 
      // (i.e. no mass flow through the wall)
      
      Linear_Reconstruction_LeastSquares_2(SolnBlk, 
					   i, SolnBlk.JCu, 
					   LIMITER_BARTH_JESPERSEN);
      
      for(ghost = 1; ghost <= SolnBlk.Nghost; ++ghost){
	// Implement reflection of velocities
	SolnBlk.W[i][SolnBlk.JCu+ghost] = Reflect(SolnBlk.W[i][SolnBlk.JCu-ghost+1],
						  SolnBlk.Grid.nfaceN(i, SolnBlk.JCu));
	// Implement pressure and density based on linear extrapolation
	W = SolnBlk.PiecewiseLinearSolutionAtLocation(i,SolnBlk.JCu,
						      SolnBlk.Grid.Cell[i][SolnBlk.JCu+ghost].Xc);
	SolnBlk.W[i][SolnBlk.JCu+ghost].p = W.p;
	SolnBlk.W[i][SolnBlk.JCu+ghost].d = W.d;
	// Calculate the conserved variables
	SolnBlk.U[i][SolnBlk.JCu+ghost] = U(SolnBlk.W[i][SolnBlk.JCu+ghost]);
      }
      break;	  

    default:
      for(ghost = 1; ghost <= SolnBlk.Nghost; ++ghost){
	SolnBlk.W[i][SolnBlk.JCu+ghost] = SolnBlk.W[i][SolnBlk.JCu];
	SolnBlk.U[i][SolnBlk.JCu+ghost] = SolnBlk.U[i][SolnBlk.JCu];
      }
      break;
    } /* endswitch */
  } /* endfor */

  /* BC fix for corner points with burning surfaces on either side. */

  if (SolnBlk.Grid.BCtypeW[SolnBlk.JCl] == BC_BURNING_SURFACE &&
      SolnBlk.Grid.BCtypeS[SolnBlk.ICl] == BC_BURNING_SURFACE) {
    SolnBlk.W[SolnBlk.ICl-1][SolnBlk.JCl-1] = HALF*(SolnBlk.W[SolnBlk.ICl-1][SolnBlk.JCl]+
						    SolnBlk.W[SolnBlk.ICl][SolnBlk.JCl-1]);
    SolnBlk.U[SolnBlk.ICl-1][SolnBlk.JCl-1] = U(SolnBlk.W[SolnBlk.ICl-1][SolnBlk.JCl-1]);
  } /* endif */

  if (SolnBlk.Grid.BCtypeW[SolnBlk.JCu] == BC_BURNING_SURFACE &&
      SolnBlk.Grid.BCtypeN[SolnBlk.ICl] == BC_BURNING_SURFACE) {
    SolnBlk.W[SolnBlk.ICl-1][SolnBlk.JCu+1] = HALF*(SolnBlk.W[SolnBlk.ICl-1][SolnBlk.JCu]+
						    SolnBlk.W[SolnBlk.ICl][SolnBlk.JCu+1]);
    SolnBlk.U[SolnBlk.ICl-1][SolnBlk.JCu+1] = U(SolnBlk.W[SolnBlk.ICl-1][SolnBlk.JCu+1]);
  } /* endif */

  if (SolnBlk.Grid.BCtypeE[SolnBlk.JCl] == BC_BURNING_SURFACE &&
      SolnBlk.Grid.BCtypeS[SolnBlk.ICu] == BC_BURNING_SURFACE) {
    SolnBlk.W[SolnBlk.ICu+1][SolnBlk.JCl-1] = HALF*(SolnBlk.W[SolnBlk.ICu+1][SolnBlk.JCl]+
						    SolnBlk.W[SolnBlk.ICu][SolnBlk.JCl-1]);
    SolnBlk.U[SolnBlk.ICu+1][SolnBlk.JCl-1] = U(SolnBlk.W[SolnBlk.ICu+1][SolnBlk.JCl-1]);
  } /* endif */

  if (SolnBlk.Grid.BCtypeE[SolnBlk.JCu] == BC_BURNING_SURFACE &&
      SolnBlk.Grid.BCtypeN[SolnBlk.ICu] == BC_BURNING_SURFACE) {
    SolnBlk.W[SolnBlk.ICu+1][SolnBlk.JCu+1] = HALF*(SolnBlk.W[SolnBlk.ICu+1][SolnBlk.JCu]+
						    SolnBlk.W[SolnBlk.ICu][SolnBlk.JCu+1]);
    SolnBlk.U[SolnBlk.ICu+1][SolnBlk.JCu+1] = U(SolnBlk.W[SolnBlk.ICu+1][SolnBlk.JCu+1]);
  } /* endif */


  // Impose high-order boundary conditions
  SolnBlk.BCs_HighOrder();

}

/********************************************************
 * Routine: CFL                                         *
 *                                                      *
 * Determines the allowable global and local time steps *
 * (for explicit Euler time stepping scheme) for the    *
 * specified quadrilateral solution block according to  *
 * the Courant-Friedrichs-Lewy condition.               *
 *                                                      *
 ********************************************************/
double CFL(Euler2D_Quad_Block &SolnBlk,
           Euler2D_Input_Parameters &Input_Parameters) {

    int i, j;
    double dtMin, d_i, d_j, v_i, v_j, a;
    double mr, aa_i, aa_j;

    dtMin = MILLION;

    for ( j  = SolnBlk.JCl-SolnBlk.Nghost ; j <= SolnBlk.JCu+SolnBlk.Nghost ; ++j ) {
       for ( i = SolnBlk.ICl-SolnBlk.Nghost ; i <= SolnBlk.ICu+SolnBlk.Nghost ; ++i ) {
 	  if (i < SolnBlk.ICl || i > SolnBlk.ICu ||
              j < SolnBlk.JCl || j > SolnBlk.JCu) {
	     SolnBlk.dt[i][j] = ZERO;
          } else {
   	     d_i = TWO*(SolnBlk.Grid.Cell[i][j].A/
                   (SolnBlk.Grid.lfaceE(i, j)+SolnBlk.Grid.lfaceW(i, j)));
             d_j = TWO*(SolnBlk.Grid.Cell[i][j].A/
                   (SolnBlk.Grid.lfaceN(i, j)+SolnBlk.Grid.lfaceS(i, j)));

	     v_i = HALF*(SolnBlk.W[i][j].v*
                   (SolnBlk.Grid.nfaceE(i, j)-SolnBlk.Grid.nfaceW(i, j)));
	     v_j = HALF*(SolnBlk.W[i][j].v*
                   (SolnBlk.Grid.nfaceN(i, j)-SolnBlk.Grid.nfaceS(i, j)));

             a = SolnBlk.W[i][j].a();

             if (Input_Parameters.Local_Time_Stepping != LOW_MACH_NUMBER_WEISS_SMITH_PRECONDITIONER ||
                 ((Input_Parameters.Local_Time_Stepping == LOW_MACH_NUMBER_WEISS_SMITH_PRECONDITIONER) &&
                  (Input_Parameters.i_Flux_Function != FLUX_FUNCTION_ROE_PRECON_WS)) &&
                 ((Input_Parameters.Local_Time_Stepping == LOW_MACH_NUMBER_WEISS_SMITH_PRECONDITIONER) &&
                  (Input_Parameters.i_Flux_Function != FLUX_FUNCTION_HLLE_PRECON_WS))) {
 	        SolnBlk.dt[i][j] = min(d_i/(a+fabs(v_i)), d_j/(a+fabs(v_j)));
             } else {
	        mr = SolnBlk.W[i][j].Mr();

                aa_i = sqrt(QUARTER*sqr(v_i)*(mr*mr*mr*mr-TWO*mr*mr+ONE)+a*a*mr*mr);
                aa_j = sqrt(QUARTER*sqr(v_j)*(mr*mr*mr*mr-TWO*mr*mr+ONE)+a*a*mr*mr);
                SolnBlk.dt[i][j] = min(d_i/(aa_i+fabs(v_i)*HALF*(ONE+mr*mr)), 
                                       d_j/(aa_j+fabs(v_j)*HALF*(ONE+mr*mr)));
             } /* endif */

             dtMin = min(dtMin, SolnBlk.dt[i][j]);

          } /* endif */
       } /* endfor */
    } /* endfor */

    for ( j  = SolnBlk.JCl-SolnBlk.Nghost ; j <= SolnBlk.JCu+SolnBlk.Nghost ; ++j ) {
       for ( i = SolnBlk.ICl-SolnBlk.Nghost ; i <= SolnBlk.ICu+SolnBlk.Nghost ; ++i ) {
 	  if (i < SolnBlk.ICl || i > SolnBlk.ICu ||
              j < SolnBlk.JCl || j > SolnBlk.JCu) {
	     SolnBlk.dt[i][j] = dtMin;
          } /* endif */
       } /* endfor */
    } /* endfor */

    /* Return the global time step. */

    return (dtMin);

}

/********************************************************
 * Routine: Set_Global_TimeStep                         *
 *                                                      *
 * Assigns global time step to specified solution block *
 * for time-accurate calculations.                      *
 *                                                      *
 ********************************************************/
void Set_Global_TimeStep(Euler2D_Quad_Block &SolnBlk,
                         const double &Dt_min) {

    int i, j;

    for ( j  = SolnBlk.JCl ; j <= SolnBlk.JCu ; ++j ) {
       for ( i = SolnBlk.ICl ; i <= SolnBlk.ICu ; ++i ) {
          SolnBlk.dt[i][j] = Dt_min;
       } /* endfor */
    } /* endfor */

}

/********************************************************
 * Routine: L1_Norm_Residual                            *
 *                                                      *
 * Determines the L1-norm of the solution residual for  *
 * the specified quadrilateral solution block.          *
 * Useful for monitoring convergence of the solution    *
 * for steady state problems.                           *
 *                                                      *
 ********************************************************/
double L1_Norm_Residual(Euler2D_Quad_Block &SolnBlk, const int &norm) {

    double l1_norm(ZERO);
    for ( int j  = SolnBlk.JCl ; j <= SolnBlk.JCu ; ++j ) {
       for ( int i = SolnBlk.ICl ; i <= SolnBlk.ICu ; ++i ) {
          l1_norm += fabs(SolnBlk.dUdt[i][j][0][norm]);
       } 
    } 
    return (l1_norm);
}

/********************************************************
 * Routine: L2_Norm_Residual                            *
 *                                                      *
 * Determines the L2-norm of the solution residual for  *
 * the specified quadrilateral solution block.          *
 * Useful for monitoring convergence of the solution    *
 * for steady state problems.                           *
 *                                                      *
 ********************************************************/
double L2_Norm_Residual(Euler2D_Quad_Block &SolnBlk, const int &norm) {


  double l2_norm(ZERO);
  
  for (int j  = SolnBlk.JCl ; j <= SolnBlk.JCu ; ++j ) {
    for (int i = SolnBlk.ICl ; i <= SolnBlk.ICu ; ++i ) {
      l2_norm += sqr(SolnBlk.dUdt[i][j][0][norm]);
    } 
  } 
  l2_norm = sqrt(l2_norm);
  
  return (l2_norm);

}

/********************************************************
 * Routine: Max_Norm_Residual                           *
 *                                                      *
 * Determines the maximum norm of the solution residual *
 * for the specified quadrilateral solution block.      *
 * Useful for monitoring convergence of the solution    *
 * for steady state problems.                           *
 *                                                      *
 ********************************************************/
double Max_Norm_Residual(Euler2D_Quad_Block &SolnBlk, const int &norm) {

  double max_norm(ZERO);

  for (int j  = SolnBlk.JCl ; j <= SolnBlk.JCu ; ++j ) {
    for (int i = SolnBlk.ICl ; i <= SolnBlk.ICu ; ++i ) {
      max_norm = max(max_norm, fabs(SolnBlk.dUdt[i][j][0][norm]));
    } 
  } 
  return (max_norm);

}

/********************************************************
 * Routine: Linear_Reconstruction_GreenGauss            *
 *                                                      *
 * Performs the reconstruction of a limited piecewise   *
 * linear solution state within a given cell (i,j) of   *
 * the computational mesh for the specified             *
 * quadrilateral solution block.  A Green-Gauss         *
 * approach is used in the evaluation of the unlimited  *
 * solution gradients.  Several slope limiters may be   *
 * used.                                                *
 *                                                      *
 ********************************************************/
void Linear_Reconstruction_GreenGauss(Euler2D_Quad_Block &SolnBlk,
				      const int i, 
                                      const int j,
                                      const int Limiter) {

    int n, n2, n_pts, i_index[8], j_index[8];
    double u0Min, u0Max, uQuad[4], phi;
    double DxDx_ave, DxDy_ave, DyDy_ave;
    double l_north, l_south, l_east, l_west;
    Vector2D n_north, n_south, n_east, n_west, dX;
    Euler2D_pState W_nw, W_ne, W_sw, W_se, W_face, 
                   DU, DUDx_ave, DUDy_ave;

    /* Carry out the limited solution reconstruction in
       the specified cell of the computational mesh. */

    // Determine the number of neighbouring cells to
    // be used in the reconstruction procedure.  Away from
    // boundaries this 8 neighbours will be used.
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
    } else if ((i == SolnBlk.ICu+SolnBlk.Nghost-1) && 
               (SolnBlk.Grid.BCtypeE[j] != BC_NONE)) {
      if (j == SolnBlk.JCl-SolnBlk.Nghost+1 || j == SolnBlk.JCu+SolnBlk.Nghost-1) {
         n_pts = 0;
      } else if (SolnBlk.Grid.BCtypeE[j] == BC_PERIODIC ||
                 SolnBlk.Grid.BCtypeE[j] == BC_CONSTANT_EXTRAPOLATION ||
                 SolnBlk.Grid.BCtypeE[j] == BC_LINEAR_EXTRAPOLATION ||
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
    } else if ((j == SolnBlk.JCl-SolnBlk.Nghost+1) && 
               (SolnBlk.Grid.BCtypeS[i] != BC_NONE)) {
      if (i == SolnBlk.ICl-SolnBlk.Nghost+1 || i == SolnBlk.ICu+SolnBlk.Nghost-1) {
         n_pts = 0;
      } else if (SolnBlk.Grid.BCtypeS[i] == BC_PERIODIC ||
                 SolnBlk.Grid.BCtypeS[i] == BC_CONSTANT_EXTRAPOLATION ||
                 SolnBlk.Grid.BCtypeS[i] == BC_LINEAR_EXTRAPOLATION ||
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
    } else if ((j == SolnBlk.JCu+SolnBlk.Nghost-1) && 
               (SolnBlk.Grid.BCtypeN[i] != BC_NONE)) {
      if (i == SolnBlk.ICl-SolnBlk.Nghost+1 || i == SolnBlk.ICu+SolnBlk.Nghost-1) {
         n_pts = 0;
      } else if (SolnBlk.Grid.BCtypeN[i] == BC_PERIODIC ||
                 SolnBlk.Grid.BCtypeN[i] == BC_CONSTANT_EXTRAPOLATION ||
                 SolnBlk.Grid.BCtypeN[i] == BC_LINEAR_EXTRAPOLATION ||
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
//      } else if (((i == SolnBlk.ICl && j == SolnBlk.JCl) && 
//                  ((SolnBlk.Grid.BCtypeW[j] == BC_REFLECTION &&
//                    SolnBlk.Grid.BCtypeS[i] == BC_REFLECTION) ||
//                   (SolnBlk.Grid.BCtypeW[j] == BC_BURNING_SURFACE &&
//                    SolnBlk.Grid.BCtypeS[i] == BC_BURNING_SURFACE))) || 
//                 ((i == SolnBlk.ICl && j == SolnBlk.JCu) && 
//                  ((SolnBlk.Grid.BCtypeW[j] == BC_REFLECTION &&
//                    SolnBlk.Grid.BCtypeN[i] == BC_REFLECTION) ||
//                   (SolnBlk.Grid.BCtypeW[j] == BC_BURNING_SURFACE &&
//                    SolnBlk.Grid.BCtypeN[i] == BC_BURNING_SURFACE))) ||
//                 ((i == SolnBlk.ICu && j == SolnBlk.JCl) && 
//                  ((SolnBlk.Grid.BCtypeE[j] == BC_REFLECTION &&
//                    SolnBlk.Grid.BCtypeS[i] == BC_REFLECTION) ||
//                   (SolnBlk.Grid.BCtypeE[j] == BC_BURNING_SURFACE &&
//                    SolnBlk.Grid.BCtypeS[i] == BC_BURNING_SURFACE))) ||
//                 ((i == SolnBlk.ICu && j == SolnBlk.JCu) && 
//                  ((SolnBlk.Grid.BCtypeE[j] == BC_REFLECTION &&
//                    SolnBlk.Grid.BCtypeN[i] == BC_REFLECTION) ||
//                   (SolnBlk.Grid.BCtypeE[j] == BC_BURNING_SURFACE &&
//                    SolnBlk.Grid.BCtypeN[i] == BC_BURNING_SURFACE)))) {
//        n_pts = 8;
//        i_index[0] = i-1; j_index[0] = j-1;
//        i_index[1] = i  ; j_index[1] = j-1;
//        i_index[2] = i+1; j_index[2] = j-1;
//        i_index[3] = i-1; j_index[3] = j  ;
//        i_index[4] = i+1; j_index[4] = j  ;
//        i_index[5] = i-1; j_index[5] = j+1;
//        i_index[6] = i  ; j_index[6] = j+1;
//        i_index[7] = i+1; j_index[7] = j+1;
//      } else if ((i == SolnBlk.ICl && j == SolnBlk.JCl) && 
//                 (SolnBlk.Grid.BCtypeW[j] != BC_NONE &&
//                  SolnBlk.Grid.BCtypeS[i] != BC_NONE)) {
//        n_pts = 7;
//        i_index[0] = i  ; j_index[0] = j-1;
//        i_index[1] = i+1; j_index[1] = j-1;
//        i_index[2] = i-1; j_index[2] = j  ;
//        i_index[3] = i+1; j_index[3] = j  ;
//        i_index[4] = i-1; j_index[4] = j+1;
//        i_index[5] = i  ; j_index[5] = j+1;
//        i_index[6] = i+1; j_index[6] = j+1;
//      } else if ((i == SolnBlk.ICl && j == SolnBlk.JCu) && 
//                 (SolnBlk.Grid.BCtypeW[j] != BC_NONE &&
//                  SolnBlk.Grid.BCtypeN[i] != BC_NONE)) {
//        n_pts = 7;
//        i_index[0] = i-1; j_index[0] = j-1;
//        i_index[1] = i  ; j_index[1] = j-1;
//        i_index[2] = i+1; j_index[2] = j-1;
//        i_index[3] = i-1; j_index[3] = j  ;
//        i_index[4] = i+1; j_index[4] = j  ;
//        i_index[5] = i  ; j_index[5] = j+1;
//        i_index[6] = i+1; j_index[6] = j+1;
//      } else if ((i == SolnBlk.ICu && j == SolnBlk.JCl) && 
//                 (SolnBlk.Grid.BCtypeE[j] != BC_NONE &&
//                  SolnBlk.Grid.BCtypeS[i] != BC_NONE)) {
//        n_pts = 7;
//        i_index[0] = i-1; j_index[0] = j-1;
//        i_index[1] = i  ; j_index[1] = j-1;
//        i_index[2] = i-1; j_index[2] = j  ;
//        i_index[3] = i+1; j_index[3] = j  ;
//        i_index[4] = i-1; j_index[4] = j+1;
//        i_index[5] = i  ; j_index[5] = j+1;
//        i_index[6] = i+1; j_index[6] = j+1;
//      } else if ((i == SolnBlk.ICu && j == SolnBlk.JCu) && 
//                 (SolnBlk.Grid.BCtypeE[j] != BC_NONE &&
//                  SolnBlk.Grid.BCtypeN[i] != BC_NONE)) {
//        n_pts = 7;
//        i_index[0] = i-1; j_index[0] = j-1;
//        i_index[1] = i  ; j_index[1] = j-1;
//        i_index[2] = i+1; j_index[2] = j-1;
//        i_index[3] = i-1; j_index[3] = j  ;
//        i_index[4] = i+1; j_index[4] = j  ;
//        i_index[5] = i-1; j_index[5] = j+1;
//        i_index[6] = i  ; j_index[6] = j+1;
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
           W_nw = SolnBlk.WnNW(i, j);
           W_ne = SolnBlk.WnNE(i, j);
           W_sw = SolnBlk.WnSW(i, j);
           W_se = SolnBlk.WnSE(i, j);

           l_north = SolnBlk.Grid.lfaceN(i, j);
           l_south = SolnBlk.Grid.lfaceS(i, j);
           l_east = SolnBlk.Grid.lfaceE(i, j);
           l_west = SolnBlk.Grid.lfaceW(i, j);

           n_north = SolnBlk.Grid.nfaceN(i, j);
           n_south = SolnBlk.Grid.nfaceS(i, j);
           n_east = SolnBlk.Grid.nfaceE(i, j);
           n_west = SolnBlk.Grid.nfaceW(i, j);

           W_face = HALF*(W_nw+W_ne)*l_north; 
           SolnBlk.dWdx[i][j] = W_face*n_north.x;
	   SolnBlk.dWdy[i][j] = W_face*n_north.y;

           W_face = HALF*(W_sw+W_se)*l_south; 
           SolnBlk.dWdx[i][j] += W_face*n_south.x;
	   SolnBlk.dWdy[i][j] += W_face*n_south.y;

           W_face = HALF*(W_ne+W_se)*l_east; 
           SolnBlk.dWdx[i][j] += W_face*n_east.x;
	   SolnBlk.dWdy[i][j] += W_face*n_east.y;

           W_face = HALF*(W_nw+W_sw)*l_west; 
           SolnBlk.dWdx[i][j] += W_face*n_west.x;
	   SolnBlk.dWdy[i][j] += W_face*n_west.y;

           SolnBlk.dWdx[i][j] = SolnBlk.dWdx[i][j]/
                                SolnBlk.Grid.Cell[i][j].A;
           SolnBlk.dWdy[i][j] = SolnBlk.dWdy[i][j]/
                                SolnBlk.Grid.Cell[i][j].A;

        // If <8 neighbours are used, apply least-squares reconstruction
        } else {
           DUDx_ave = Euler2D_W_VACUUM;
           DUDy_ave = Euler2D_W_VACUUM;
           DxDx_ave = ZERO;
           DxDy_ave = ZERO;
           DyDy_ave = ZERO;
    
           for ( n2 = 0 ; n2 <= n_pts-1 ; ++n2 ) {
               dX = SolnBlk.Grid.Cell[ i_index[n2] ][ j_index[n2] ].Xc - 
                    SolnBlk.Grid.Cell[i][j].Xc;
               DU = SolnBlk.W[ i_index[n2] ][ j_index[n2] ] - 
                    SolnBlk.W[i][j];
               DUDx_ave += DU*dX.x;
               DUDy_ave += DU*dX.y;
               DxDx_ave += dX.x*dX.x;
               DxDy_ave += dX.x*dX.y;
               DyDy_ave += dX.y*dX.y;
           } /* endfor */
    					    
           DUDx_ave = DUDx_ave/double(n_pts);
           DUDy_ave = DUDy_ave/double(n_pts);
           DxDx_ave = DxDx_ave/double(n_pts);
           DxDy_ave = DxDy_ave/double(n_pts);
           DyDy_ave = DyDy_ave/double(n_pts);
           SolnBlk.dWdx[i][j] = (DUDx_ave*DyDy_ave-DUDy_ave*DxDy_ave)/
                                (DxDx_ave*DyDy_ave-DxDy_ave*DxDy_ave);
           SolnBlk.dWdy[i][j] = (DUDy_ave*DxDx_ave-DUDx_ave*DxDy_ave)/
                                (DxDx_ave*DyDy_ave-DxDy_ave*DxDy_ave);
        } /* endif */

        // Calculate slope limiters.    
	if (!SolnBlk.Freeze_Limiter) {
           for ( n = 1 ; n <= NUM_VAR_EULER2D ; ++n ) {
              u0Min = SolnBlk.W[i][j][n];
              u0Max = u0Min;
              for ( n2 = 0 ; n2 <= n_pts-1 ; ++n2 ) {
                 u0Min = min(u0Min, SolnBlk.W[ i_index[n2] ][ j_index[n2] ][n]);
                 u0Max = max(u0Max, SolnBlk.W[ i_index[n2] ][ j_index[n2] ][n]);
              } /* endfor */
    
              dX = SolnBlk.Grid.xfaceE(i, j)-SolnBlk.Grid.Cell[i][j].Xc;
              uQuad[0] = SolnBlk.W[i][j][n] + 
                         SolnBlk.dWdx[i][j][n]*dX.x +
                         SolnBlk.dWdy[i][j][n]*dX.y ;
              dX = SolnBlk.Grid.xfaceW(i, j)-SolnBlk.Grid.Cell[i][j].Xc;
              uQuad[1] = SolnBlk.W[i][j][n] + 
                         SolnBlk.dWdx[i][j][n]*dX.x +
                         SolnBlk.dWdy[i][j][n]*dX.y ;
              dX = SolnBlk.Grid.xfaceN(i, j)-SolnBlk.Grid.Cell[i][j].Xc;
              uQuad[2] = SolnBlk.W[i][j][n] + 
                         SolnBlk.dWdx[i][j][n]*dX.x +
                         SolnBlk.dWdy[i][j][n]*dX.y ;
              dX = SolnBlk.Grid.xfaceS(i, j)-SolnBlk.Grid.Cell[i][j].Xc;
              uQuad[3] = SolnBlk.W[i][j][n] + 
                         SolnBlk.dWdx[i][j][n]*dX.x +
                         SolnBlk.dWdy[i][j][n]*dX.y ;
    
              switch(Limiter) {
                case LIMITER_ONE :
                  phi = ONE;
                  break;
                case LIMITER_ZERO :
                  phi = ZERO;
                  break;
                case LIMITER_BARTH_JESPERSEN :
                  phi = Limiter_BarthJespersen(uQuad, SolnBlk.W[i][j][n], 
                                               u0Min, u0Max, 4);
                  break;
                case LIMITER_VENKATAKRISHNAN :
                  phi = Limiter_Venkatakrishnan(uQuad, SolnBlk.W[i][j][n], 
                                                u0Min, u0Max, 4);
                  break;
                case LIMITER_VANLEER :
                  phi = Limiter_VanLeer(uQuad, SolnBlk.W[i][j][n], 
                                        u0Min, u0Max, 4);
                  break;
                case LIMITER_VANALBADA :
                  phi = Limiter_VanAlbada(uQuad, SolnBlk.W[i][j][n], 
                                          u0Min, u0Max, 4);
                  break;
                default:
                  phi = Limiter_BarthJespersen(uQuad, SolnBlk.W[i][j][n], 
                                               u0Min, u0Max, 4);
                  break;
              } /* endswitch */

	      SolnBlk.phi[i][j][n] = phi;

           } /* endfor */
        } /* endif */
    } else {
        SolnBlk.dWdx[i][j] = Euler2D_W_VACUUM;
        SolnBlk.dWdy[i][j] = Euler2D_W_VACUUM; 
        SolnBlk.phi[i][j]  = Euler2D_W_VACUUM;
    } /* endif */
    
}

/********************************************************
 * Routine: Linear_Reconstruction_GreenGauss            *
 *                                                      *
 * Performs the reconstruction of a limited piecewise   *
 * linear solution state within each cell of the        *
 * computational mesh for the specified quadrilateral   *
 * solution block.  A Green-Gauss approach is used      *
 * in the evaluation of the unlimited solution          *
 * gradients.  Several slope limiters may be used.      *
 *                                                      *
 ********************************************************/
void Linear_Reconstruction_GreenGauss(Euler2D_Quad_Block &SolnBlk,
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

/********************************************************
 * Routine: Linear_Reconstruction_LeastSquares          *
 *                                                      *
 * Performs the reconstruction of a limited piecewise   *
 * linear solution state within a given cell (i,j) of   *
 * the computational mesh for the specified             *
 * quadrilateral solution block.  A least squares       *
 * approach is used in the evaluation of the unlimited  *
 * solution gradients.  Several slope limiters may be   *
 * used.                                                *
 *                                                      *
 ********************************************************/
void Linear_Reconstruction_LeastSquares(Euler2D_Quad_Block &SolnBlk,
				        const int i, 
                                        const int j,
                                        const int Limiter) {

    int n, n2, n_pts, i_index[8], j_index[8];
    double u0Min, u0Max, uQuad[4], phi;
    double DxDx_ave, DxDy_ave, DyDy_ave;
    Vector2D dX;
    Euler2D_pState DU, DUDx_ave, DUDy_ave;

    /* Carry out the limited solution reconstruction in
       each cell of the computational mesh. */

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
    } else if ((i == SolnBlk.ICu+SolnBlk.Nghost-1) && 
               (SolnBlk.Grid.BCtypeE[j] != BC_NONE)) {
      if (j == SolnBlk.JCl-SolnBlk.Nghost+1 || j == SolnBlk.JCu+SolnBlk.Nghost-1) {
         n_pts = 0;
      } else if (SolnBlk.Grid.BCtypeE[j] == BC_PERIODIC ||
                 SolnBlk.Grid.BCtypeE[j] == BC_CONSTANT_EXTRAPOLATION ||
                 SolnBlk.Grid.BCtypeE[j] == BC_LINEAR_EXTRAPOLATION ||
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
    } else if ((j == SolnBlk.JCl-SolnBlk.Nghost+1) && 
               (SolnBlk.Grid.BCtypeS[i] != BC_NONE)) {
      if (i == SolnBlk.ICl-SolnBlk.Nghost+1 || i == SolnBlk.ICu+SolnBlk.Nghost-1) {
         n_pts = 0;
      } else if (SolnBlk.Grid.BCtypeS[i] == BC_PERIODIC ||
                 SolnBlk.Grid.BCtypeS[i] == BC_CONSTANT_EXTRAPOLATION ||
                 SolnBlk.Grid.BCtypeS[i] == BC_LINEAR_EXTRAPOLATION ||
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
    } else if ((j == SolnBlk.JCu+SolnBlk.Nghost-1) && 
               (SolnBlk.Grid.BCtypeN[i] != BC_NONE)) {
      if (i == SolnBlk.ICl-SolnBlk.Nghost+1 || i == SolnBlk.ICu+SolnBlk.Nghost-1) {
         n_pts = 0;
      } else if (SolnBlk.Grid.BCtypeN[i] == BC_PERIODIC ||
                 SolnBlk.Grid.BCtypeN[i] == BC_CONSTANT_EXTRAPOLATION ||
                 SolnBlk.Grid.BCtypeN[i] == BC_LINEAR_EXTRAPOLATION ||
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
//      } else if (((i == SolnBlk.ICl && j == SolnBlk.JCl) && 
//                  ((SolnBlk.Grid.BCtypeW[j] == BC_REFLECTION &&
//                    SolnBlk.Grid.BCtypeS[i] == BC_REFLECTION) ||
//                   (SolnBlk.Grid.BCtypeW[j] == BC_BURNING_SURFACE &&
//                    SolnBlk.Grid.BCtypeS[i] == BC_BURNING_SURFACE))) || 
//                 ((i == SolnBlk.ICl && j == SolnBlk.JCu) && 
//                  ((SolnBlk.Grid.BCtypeW[j] == BC_REFLECTION &&
//                    SolnBlk.Grid.BCtypeN[i] == BC_REFLECTION) ||
//                   (SolnBlk.Grid.BCtypeW[j] == BC_BURNING_SURFACE &&
//                    SolnBlk.Grid.BCtypeN[i] == BC_BURNING_SURFACE))) ||
//                 ((i == SolnBlk.ICu && j == SolnBlk.JCl) && 
//                  ((SolnBlk.Grid.BCtypeE[j] == BC_REFLECTION &&
//                    SolnBlk.Grid.BCtypeS[i] == BC_REFLECTION) ||
//                   (SolnBlk.Grid.BCtypeE[j] == BC_BURNING_SURFACE &&
//                    SolnBlk.Grid.BCtypeS[i] == BC_BURNING_SURFACE))) ||
//                 ((i == SolnBlk.ICu && j == SolnBlk.JCu) && 
//                  ((SolnBlk.Grid.BCtypeE[j] == BC_REFLECTION &&
//                    SolnBlk.Grid.BCtypeN[i] == BC_REFLECTION) ||
//                   (SolnBlk.Grid.BCtypeE[j] == BC_BURNING_SURFACE &&
//                    SolnBlk.Grid.BCtypeN[i] == BC_BURNING_SURFACE)))) {
//        n_pts = 8;
//        i_index[0] = i-1; j_index[0] = j-1;
//        i_index[1] = i  ; j_index[1] = j-1;
//        i_index[2] = i+1; j_index[2] = j-1;
//        i_index[3] = i-1; j_index[3] = j  ;
//        i_index[4] = i+1; j_index[4] = j  ;
//        i_index[5] = i-1; j_index[5] = j+1;
//        i_index[6] = i  ; j_index[6] = j+1;
//        i_index[7] = i+1; j_index[7] = j+1;
//      } else if ((i == SolnBlk.ICl && j == SolnBlk.JCl) && 
//                 (SolnBlk.Grid.BCtypeW[j] != BC_NONE &&
//                  SolnBlk.Grid.BCtypeS[i] != BC_NONE)) {
//        n_pts = 7;
//        i_index[0] = i  ; j_index[0] = j-1;
//        i_index[1] = i+1; j_index[1] = j-1;
//        i_index[2] = i-1; j_index[2] = j  ;
//        i_index[3] = i+1; j_index[3] = j  ;
//        i_index[4] = i-1; j_index[4] = j+1;
//        i_index[5] = i  ; j_index[5] = j+1;
//        i_index[6] = i+1; j_index[6] = j+1;
//      } else if ((i == SolnBlk.ICl && j == SolnBlk.JCu) && 
//                 (SolnBlk.Grid.BCtypeW[j] != BC_NONE &&
//                  SolnBlk.Grid.BCtypeN[i] != BC_NONE)) {
//        n_pts = 7;
//        i_index[0] = i-1; j_index[0] = j-1;
//        i_index[1] = i  ; j_index[1] = j-1;
//        i_index[2] = i+1; j_index[2] = j-1;
//        i_index[3] = i-1; j_index[3] = j  ;
//        i_index[4] = i+1; j_index[4] = j  ;
//        i_index[5] = i  ; j_index[5] = j+1;
//        i_index[6] = i+1; j_index[6] = j+1;
//      } else if ((i == SolnBlk.ICu && j == SolnBlk.JCl) && 
//                 (SolnBlk.Grid.BCtypeE[j] != BC_NONE &&
//                  SolnBlk.Grid.BCtypeS[i] != BC_NONE)) {
//        n_pts = 7;
//        i_index[0] = i-1; j_index[0] = j-1;
//        i_index[1] = i  ; j_index[1] = j-1;
//        i_index[2] = i-1; j_index[2] = j  ;
//        i_index[3] = i+1; j_index[3] = j  ;
//        i_index[4] = i-1; j_index[4] = j+1;
//        i_index[5] = i  ; j_index[5] = j+1;
//        i_index[6] = i+1; j_index[6] = j+1;
//      } else if ((i == SolnBlk.ICu && j == SolnBlk.JCu) && 
//                 (SolnBlk.Grid.BCtypeE[j] != BC_NONE &&
//                  SolnBlk.Grid.BCtypeN[i] != BC_NONE)) {
//        n_pts = 7;
//        i_index[0] = i-1; j_index[0] = j-1;
//        i_index[1] = i  ; j_index[1] = j-1;
//        i_index[2] = i+1; j_index[2] = j-1;
//        i_index[3] = i-1; j_index[3] = j  ;
//        i_index[4] = i+1; j_index[4] = j  ;
//        i_index[5] = i-1; j_index[5] = j+1;
//        i_index[6] = i  ; j_index[6] = j+1;
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
    
    if (n_pts > 0) {
        DUDx_ave = Euler2D_W_VACUUM;
        DUDy_ave = Euler2D_W_VACUUM;
        DxDx_ave = ZERO;
        DxDy_ave = ZERO;
        DyDy_ave = ZERO;
    
        for ( n2 = 0 ; n2 <= n_pts-1 ; ++n2 ) {
            dX = SolnBlk.Grid.Cell[ i_index[n2] ][ j_index[n2] ].Xc - 
                 SolnBlk.Grid.Cell[i][j].Xc;
            DU = SolnBlk.W[ i_index[n2] ][ j_index[n2] ] - 
                 SolnBlk.W[i][j];
            DUDx_ave += DU*dX.x;
            DUDy_ave += DU*dX.y;
            DxDx_ave += dX.x*dX.x;
            DxDy_ave += dX.x*dX.y;
            DyDy_ave += dX.y*dX.y;
        } /* endfor */
    					    
        DUDx_ave = DUDx_ave/double(n_pts);
        DUDy_ave = DUDy_ave/double(n_pts);
        DxDx_ave = DxDx_ave/double(n_pts);
        DxDy_ave = DxDy_ave/double(n_pts);
        DyDy_ave = DyDy_ave/double(n_pts);
        SolnBlk.dWdx[i][j] = (DUDx_ave*DyDy_ave-DUDy_ave*DxDy_ave)/
                             (DxDx_ave*DyDy_ave-DxDy_ave*DxDy_ave);
        SolnBlk.dWdy[i][j] = (DUDy_ave*DxDx_ave-DUDx_ave*DxDy_ave)/
                             (DxDx_ave*DyDy_ave-DxDy_ave*DxDy_ave);
    
        // Calculate slope limiters. 
	if (!SolnBlk.Freeze_Limiter) {
           for ( n = 1 ; n <= NUM_VAR_EULER2D ; ++n ) {
              u0Min = SolnBlk.W[i][j][n];
              u0Max = u0Min;
              for ( n2 = 0 ; n2 <= n_pts-1 ; ++n2 ) {
                 u0Min = min(u0Min, SolnBlk.W[ i_index[n2] ][ j_index[n2] ][n]);
                 u0Max = max(u0Max, SolnBlk.W[ i_index[n2] ][ j_index[n2] ][n]);
              } /* endfor */
    
              dX = SolnBlk.Grid.xfaceE(i, j)-SolnBlk.Grid.Cell[i][j].Xc;
              uQuad[0] = SolnBlk.W[i][j][n] + 
                         SolnBlk.dWdx[i][j][n]*dX.x +
                         SolnBlk.dWdy[i][j][n]*dX.y ;
              dX = SolnBlk.Grid.xfaceW(i, j)-SolnBlk.Grid.Cell[i][j].Xc;
              uQuad[1] = SolnBlk.W[i][j][n] + 
                         SolnBlk.dWdx[i][j][n]*dX.x +
                         SolnBlk.dWdy[i][j][n]*dX.y ;
              dX = SolnBlk.Grid.xfaceN(i, j)-SolnBlk.Grid.Cell[i][j].Xc;
              uQuad[2] = SolnBlk.W[i][j][n] + 
                         SolnBlk.dWdx[i][j][n]*dX.x +
                         SolnBlk.dWdy[i][j][n]*dX.y ;
              dX = SolnBlk.Grid.xfaceS(i, j)-SolnBlk.Grid.Cell[i][j].Xc;
              uQuad[3] = SolnBlk.W[i][j][n] + 
                         SolnBlk.dWdx[i][j][n]*dX.x +
                         SolnBlk.dWdy[i][j][n]*dX.y ;
    
              switch(Limiter) {
                case LIMITER_ONE :
                  phi = ONE;
                  break;
                case LIMITER_ZERO :
                  phi = ZERO;
                  break;
                case LIMITER_BARTH_JESPERSEN :
                  phi = Limiter_BarthJespersen(uQuad, SolnBlk.W[i][j][n], 
                                               u0Min, u0Max, 4);
                  break;
                case LIMITER_VENKATAKRISHNAN :
                  phi = Limiter_Venkatakrishnan(uQuad, SolnBlk.W[i][j][n], 
                                                u0Min, u0Max, 4);
                  break;
                case LIMITER_VANLEER :
                  phi = Limiter_VanLeer(uQuad, SolnBlk.W[i][j][n], 
                                        u0Min, u0Max, 4);
                  break;
                case LIMITER_VANALBADA :
                  phi = Limiter_VanAlbada(uQuad, SolnBlk.W[i][j][n], 
                                          u0Min, u0Max, 4);
                  break;
                default:
                  phi = Limiter_BarthJespersen(uQuad, SolnBlk.W[i][j][n], 
                                               u0Min, u0Max, 4);
                  break;
              } /* endswitch */

	      SolnBlk.phi[i][j][n] = phi;

           } /* endfor */
        } /* endif */
    } else {
        SolnBlk.dWdx[i][j] = Euler2D_W_VACUUM;
        SolnBlk.dWdy[i][j] = Euler2D_W_VACUUM; 
        SolnBlk.phi[i][j]  = Euler2D_W_VACUUM;
    } /* endif */
    
}

/********************************************************
 * Routine: Linear_Reconstruction_LeastSquares_2        *
 *                                                      *
 * Performs the reconstruction of a limited piecewise   *
 * linear solution state within a given cell (i,j) of   *
 * the computational mesh for the specified             *
 * quadrilateral solution block.  A least squares       *
 * approach is used in the evaluation of the unlimited  *
 * solution gradients.  Several slope limiters may be   *
 * used.                                                *
 *                                                      *
 ********************************************************/
void Linear_Reconstruction_LeastSquares_2(Euler2D_Quad_Block &SolnBlk,
				          const int i, 
                                          const int j,
                                          const int Limiter) {

    int n, n2, n_pts, i_index[8], j_index[8];
    double u0Min, u0Max, uQuad[4], phi;
    double DxDx_ave, DxDy_ave, DyDy_ave;
    Vector2D dX;
    Euler2D_pState DU, DUDx_ave, DUDy_ave;

    /* Carry out the limited solution reconstruction in
       each cell of the computational mesh. */

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
         } /* endif */
      } else {
         if (j == SolnBlk.JCl) {
            n_pts = 0;
         } else if (j == SolnBlk.JCu) {
            n_pts = 0;
         } else {
            n_pts = 0;
         } /* endif */
      } /* endif */           
    } else if ((i == SolnBlk.ICu+SolnBlk.Nghost-1) && 
               (SolnBlk.Grid.BCtypeE[j] != BC_NONE)) {
      if (j == SolnBlk.JCl-SolnBlk.Nghost+1 || j == SolnBlk.JCu+SolnBlk.Nghost-1) {
         n_pts = 0;
      } else if (SolnBlk.Grid.BCtypeE[j] == BC_PERIODIC ||
                 SolnBlk.Grid.BCtypeE[j] == BC_CONSTANT_EXTRAPOLATION ||
                 SolnBlk.Grid.BCtypeE[j] == BC_LINEAR_EXTRAPOLATION ||
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
         } /* endif */
      } else {
         if (j == SolnBlk.JCl) {
            n_pts = 0;
         } else if (j == SolnBlk.JCu) {
            n_pts = 0;
         } else {
            n_pts = 0;
         } /* endif */
      } /* endif */
    } else if ((j == SolnBlk.JCl-SolnBlk.Nghost+1) && 
               (SolnBlk.Grid.BCtypeS[i] != BC_NONE)) {
      if (i == SolnBlk.ICl-SolnBlk.Nghost+1 || i == SolnBlk.ICu+SolnBlk.Nghost-1) {
         n_pts = 0;
      } else if (SolnBlk.Grid.BCtypeS[i] == BC_PERIODIC ||
                 SolnBlk.Grid.BCtypeS[i] == BC_CONSTANT_EXTRAPOLATION ||
                 SolnBlk.Grid.BCtypeS[i] == BC_LINEAR_EXTRAPOLATION ||
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
         } /* endif */
      } else {
         if (i == SolnBlk.ICl) {
            n_pts = 0;
         } else if (i == SolnBlk.ICu) {
            n_pts = 0;
         } else {
            n_pts = 0;
         } /* endif */
      } /* endif */
    } else if ((j == SolnBlk.JCu+SolnBlk.Nghost-1) && 
               (SolnBlk.Grid.BCtypeN[i] != BC_NONE)) {
      if (i == SolnBlk.ICl-SolnBlk.Nghost+1 || i == SolnBlk.ICu+SolnBlk.Nghost-1) {
         n_pts = 0;
      } else if (SolnBlk.Grid.BCtypeN[i] == BC_PERIODIC ||
                 SolnBlk.Grid.BCtypeN[i] == BC_CONSTANT_EXTRAPOLATION ||
                 SolnBlk.Grid.BCtypeN[i] == BC_LINEAR_EXTRAPOLATION ||
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
         } /* endif */
      } else {
         if (i == SolnBlk.ICl) {
            n_pts = 0;
         } else if (i == SolnBlk.ICu) {
            n_pts = 0;
         } else {
            n_pts = 0;
         } /* endif */
      } /* endif */
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
    } /* endif */
    
    if (n_pts > 0) {
        DUDx_ave = Euler2D_W_VACUUM;
        DUDy_ave = Euler2D_W_VACUUM;
        DxDx_ave = ZERO;
        DxDy_ave = ZERO;
        DyDy_ave = ZERO;
    
        for ( n2 = 0 ; n2 <= n_pts-1 ; ++n2 ) {
            dX = SolnBlk.Grid.Cell[ i_index[n2] ][ j_index[n2] ].Xc - 
                 SolnBlk.Grid.Cell[i][j].Xc;
            DU = SolnBlk.W[ i_index[n2] ][ j_index[n2] ] - 
                 SolnBlk.W[i][j];
            DUDx_ave += DU*dX.x;
            DUDy_ave += DU*dX.y;
            DxDx_ave += dX.x*dX.x;
            DxDy_ave += dX.x*dX.y;
            DyDy_ave += dX.y*dX.y;
        } /* endfor */
    					    
        DUDx_ave = DUDx_ave/double(n_pts);
        DUDy_ave = DUDy_ave/double(n_pts);
        DxDx_ave = DxDx_ave/double(n_pts);
        DxDy_ave = DxDy_ave/double(n_pts);
        DyDy_ave = DyDy_ave/double(n_pts);
        SolnBlk.dWdx[i][j] = (DUDx_ave*DyDy_ave-DUDy_ave*DxDy_ave)/
                             (DxDx_ave*DyDy_ave-DxDy_ave*DxDy_ave);
        SolnBlk.dWdy[i][j] = (DUDy_ave*DxDx_ave-DUDx_ave*DxDy_ave)/
                             (DxDx_ave*DyDy_ave-DxDy_ave*DxDy_ave);
    
        // Calculate slope limiters. 
	if (!SolnBlk.Freeze_Limiter) {
           for ( n = 1 ; n <= NUM_VAR_EULER2D ; ++n ) {
              u0Min = SolnBlk.W[i][j][n];
              u0Max = u0Min;
              for ( n2 = 0 ; n2 <= n_pts-1 ; ++n2 ) {
                 u0Min = min(u0Min, SolnBlk.W[ i_index[n2] ][ j_index[n2] ][n]);
                 u0Max = max(u0Max, SolnBlk.W[ i_index[n2] ][ j_index[n2] ][n]);
              } /* endfor */
    
              dX = SolnBlk.Grid.xfaceE(i, j)-SolnBlk.Grid.Cell[i][j].Xc;
              uQuad[0] = SolnBlk.W[i][j][n] + 
                         SolnBlk.dWdx[i][j][n]*dX.x +
                         SolnBlk.dWdy[i][j][n]*dX.y ;
              dX = SolnBlk.Grid.xfaceW(i, j)-SolnBlk.Grid.Cell[i][j].Xc;
              uQuad[1] = SolnBlk.W[i][j][n] + 
                         SolnBlk.dWdx[i][j][n]*dX.x +
                         SolnBlk.dWdy[i][j][n]*dX.y ;
              dX = SolnBlk.Grid.xfaceN(i, j)-SolnBlk.Grid.Cell[i][j].Xc;
              uQuad[2] = SolnBlk.W[i][j][n] + 
                         SolnBlk.dWdx[i][j][n]*dX.x +
                         SolnBlk.dWdy[i][j][n]*dX.y ;
              dX = SolnBlk.Grid.xfaceS(i, j)-SolnBlk.Grid.Cell[i][j].Xc;
              uQuad[3] = SolnBlk.W[i][j][n] + 
                         SolnBlk.dWdx[i][j][n]*dX.x +
                         SolnBlk.dWdy[i][j][n]*dX.y ;
    
              switch(Limiter) {
                case LIMITER_ONE :
                  phi = ONE;
                  break;
                case LIMITER_ZERO :
                  phi = ZERO;
                  break;
                case LIMITER_BARTH_JESPERSEN :
                  phi = Limiter_BarthJespersen(uQuad, SolnBlk.W[i][j][n], 
                                               u0Min, u0Max, 4);
                  break;
                case LIMITER_VENKATAKRISHNAN :
                  phi = Limiter_Venkatakrishnan(uQuad, SolnBlk.W[i][j][n], 
                                                u0Min, u0Max, 4);
                  break;
                case LIMITER_VANLEER :
                  phi = Limiter_VanLeer(uQuad, SolnBlk.W[i][j][n], 
                                        u0Min, u0Max, 4);
                  break;
                case LIMITER_VANALBADA :
                  phi = Limiter_VanAlbada(uQuad, SolnBlk.W[i][j][n], 
                                          u0Min, u0Max, 4);
                  break;
                default:
                  phi = Limiter_BarthJespersen(uQuad, SolnBlk.W[i][j][n], 
                                               u0Min, u0Max, 4);
                  break;
              } /* endswitch */

	      SolnBlk.phi[i][j][n] = phi;

           } /* endfor */
        } /* endif */
    } else {
        SolnBlk.dWdx[i][j] = Euler2D_W_VACUUM;
        SolnBlk.dWdy[i][j] = Euler2D_W_VACUUM; 
        SolnBlk.phi[i][j]  = Euler2D_W_VACUUM;
    } /* endif */
    
}

/********************************************************
 * Routine: Linear_Reconstruction_LeastSquares          *
 *                                                      *
 * Performs the reconstruction of a limited piecewise   *
 * linear solution state within each cell of the        *
 * computational mesh of the specified quadrilateral    *
 * solution block.  A least squares approach is         *
 * used in the evaluation of the unlimited solution     *
 * gradients.  Several slope limiters may be used.      *
 *                                                      *
 ********************************************************/
void Linear_Reconstruction_LeastSquares(Euler2D_Quad_Block &SolnBlk,
				        const int Limiter) {

    int i, j;

    /* Carry out the limited solution reconstruction in
       each cell of the computational mesh. */

    for ( j  = SolnBlk.JCl-SolnBlk.Nghost+1 ; j <= SolnBlk.JCu+SolnBlk.Nghost-1 ; ++j ) {
       for ( i = SolnBlk.ICl-SolnBlk.Nghost+1 ; i <= SolnBlk.ICu+SolnBlk.Nghost-1 ; ++i ) {
	   Linear_Reconstruction_LeastSquares(SolnBlk, i, j, Limiter);
           //Linear_Reconstruction_LeastSquares_2(SolnBlk, i, j, Limiter);
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
void Linear_Reconstruction(Euler2D_Quad_Block &SolnBlk,
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

/********************************************************
 * Routine: Residual_Smoothing                          *
 *                                                      *
 * Applies implicit residual smoothing to solution      *
 * residual.  Note that only residuals of interior cells*
 * are smoothed and residuals for cells adjacent to     *
 * boundaries are not smoothed.                         *
 *                                                      *
 ********************************************************/
void Residual_Smoothing(Euler2D_Quad_Block &SolnBlk,
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

//              SolnBlk.dUdt[i][j][k_residual+1] = (SolnBlk.dUdt[i][j][k_residual]
//                 + epsilon*(SolnBlk.dUdt[i-1][j-1][k_residual+1] +
//                            SolnBlk.dUdt[i  ][j-1][k_residual+1] +
//                            SolnBlk.dUdt[i+1][j-1][k_residual+1] +
//                            SolnBlk.dUdt[i-1][j  ][k_residual+1] +
//                            SolnBlk.dUdt[i+1][j  ][k_residual+1] +
//                            SolnBlk.dUdt[i-1][j+1][k_residual+1] +
//                            SolnBlk.dUdt[i  ][j+1][k_residual+1] +
//                            SolnBlk.dUdt[i+1][j+1][k_residual+1]))/(ONE + EIGHT*epsilon);
          } /* endfor */
       } /* endfor */

    } /* endfor */

    for ( j  = SolnBlk.JCl+1; j <= SolnBlk.JCu-1; ++j ) {
       for ( i = SolnBlk.ICl+1; i <= SolnBlk.ICu-1; ++i ) {
          SolnBlk.dUdt[i][j][k_residual] =  SolnBlk.dUdt[i][j][k_residual+1];
       } /* endfor */
    } /* endfor */

}

/********************************************************
 * Routine: Calculate_Refinement_Criteria               *
 *                                                      *
 * Calculate refinement criteria for the solution       *
 * block.                                               *
 *                                                      *
 ********************************************************/
void Calculate_Refinement_Criteria(double *refinement_criteria,
				   Euler2D_Input_Parameters &IP,
                                   int &number_refinement_criteria,
                                   Euler2D_Quad_Block &SolnBlk) {


    // Calculate refinement criteria based on smoothness indicator
    if (CENO_Execution_Mode::USE_CENO_ALGORITHM && 
	CENO_Execution_Mode::USE_SMOOTHNESS_INDICATOR_FOR_AMR_CRITERIA) {
      return SolnBlk.Calculate_Refinement_Criteria_HighOrder(refinement_criteria,
							     IP,
							     number_refinement_criteria);
    }

    int i, j;

    double grad_rho_x, grad_rho_y, grad_rho_abs, grad_rho_criteria, grad_rho_criteria_max,
           div_V, div_V_criteria, div_V_criteria_max,
           curl_V_z, curl_V_abs, curl_V_criteria, curl_V_criteria_max;
    int refinement_criteria_number;

    /* Set the number of refinement criteria to be used (3):
       (1) refinement criteria #1 based on the gradient of the density field;
       (2) refinement criteria #2 based on the divergence of the velocity vector;
       (3) refinement criteria #3 based on the curl of the velocity vector. */

    number_refinement_criteria = IP.Number_of_Refinement_Criteria;

    /* Initialize the refinement criteria for the block. */

    grad_rho_criteria_max = ZERO;
    div_V_criteria_max = ZERO;
    curl_V_criteria_max = ZERO;

    /* Calculate the refinement criteria for each cell of the 
       computational mesh and assign the maximum value for
       all cells as the refinement criteria for the solution 
       block. */

    for ( j  = SolnBlk.JCl-1 ; j <= SolnBlk.JCu+1 ; ++j ) {
       for ( i = SolnBlk.ICl-1 ; i <= SolnBlk.ICu+1 ; ++i ) {
          // Reconstruct the solution within the cell.
	  //Linear_Reconstruction_GreenGauss(SolnBlk, i, j, LIMITER_UNLIMITED);
	  Linear_Reconstruction_LeastSquares(SolnBlk, i, j, LIMITER_UNLIMITED);
          //Linear_Reconstruction_LeastSquares_2(SolnBlk, i, j, LIMITER_UNLIMITED);

          if (SolnBlk.Grid.Cell[i][j].A > ZERO) {
	    // Evaluate refinement criteria #1 based on the gradient
	    // of the density field.
	    if (IP.Refinement_Criteria_Gradient_Density) {
	      grad_rho_x = SolnBlk.dWdx[i][j].d;
	      grad_rho_y = SolnBlk.dWdy[i][j].d;
	      grad_rho_abs = sqrt(sqr(grad_rho_x) + sqr(grad_rho_y));
	      grad_rho_criteria = sqrt(SolnBlk.Grid.Cell[i][j].A)*grad_rho_abs/SolnBlk.W[i][j].d;
	    } else {
	      grad_rho_criteria = ONE;
	    }
	    grad_rho_criteria_max = max(grad_rho_criteria_max, grad_rho_criteria);

	    // Evaluate refinement criteria #2 based on the divergence
	    // of the velocity vector.
	    if (IP.Refinement_Criteria_Divergence_Velocity) {
	      div_V = SolnBlk.dWdx[i][j].v.x + SolnBlk.dWdy[i][j].v.y;
	      div_V_criteria = sqrt(SolnBlk.Grid.Cell[i][j].A)*fabs(div_V)/SolnBlk.W[i][j].a();
	    } else {
	      div_V_criteria = ONE;
	    }
	    div_V_criteria_max = max(div_V_criteria_max, div_V_criteria);

	    // Evaluate refinement criteria #3 based on the curl
	    // of the velocity vector.
	    if (IP.Refinement_Criteria_Curl_Velocity) {
	      curl_V_z = SolnBlk.dWdx[i][j].v.y - SolnBlk.dWdy[i][j].v.x; 
	      curl_V_abs = sqrt(sqr(curl_V_z));
	      curl_V_criteria = sqrt(SolnBlk.Grid.Cell[i][j].A)*curl_V_abs/SolnBlk.W[i][j].a(); 
	    } else {
	      curl_V_criteria = ONE;
	    }
	    curl_V_criteria_max = max(curl_V_criteria_max, curl_V_criteria);
          } /* endif */

       } /* endfor */
    } /* endfor */

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

}

/********************************************************
 * Routine: Fix_Refined_Block_Boundaries                *
 *                                                      *
 * Adjusts the locations of the boundary nodes of a     *
 * solution block so that the new node locations        *
 * match with cell volumes of adjacent solution blocks  *
 * that have lower levels of mesh refinement (i.e., are *
 * coarser solution blocks).                            *
 *                                                      *
 ********************************************************/
void Fix_Refined_Block_Boundaries(Euler2D_Quad_Block &SolnBlk,
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
	  SolnBlk.W[i][SolnBlk.JCu] = W(SolnBlk.U[i][SolnBlk.JCu]);
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
	  SolnBlk.W[i][SolnBlk.JCl] = W(SolnBlk.U[i][SolnBlk.JCl]);
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
	  SolnBlk.W[SolnBlk.ICu][j] = W(SolnBlk.U[SolnBlk.ICu][j]);
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
	  SolnBlk.W[SolnBlk.ICl][j] = W(SolnBlk.U[SolnBlk.ICl][j]);
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

/********************************************************
 * Routine: Unfix_Refined_Block_Boundaries              *
 *                                                      *
 * Returns the adjusted the locations of the boundary   *
 * nodes of a solution block to their original          *
 * unmodified positions.                                *
 *                                                      *
 ********************************************************/
void Unfix_Refined_Block_Boundaries(Euler2D_Quad_Block &SolnBlk) {

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
	  SolnBlk.W[i][SolnBlk.JCu] = W(SolnBlk.U[i][SolnBlk.JCu]);
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
	  SolnBlk.W[i][SolnBlk.JCl] = W(SolnBlk.U[i][SolnBlk.JCl]);
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
	  SolnBlk.W[SolnBlk.ICu][j] = W(SolnBlk.U[SolnBlk.ICu][j]);
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
	  SolnBlk.W[SolnBlk.ICl][j] = W(SolnBlk.U[SolnBlk.ICl][j]);
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

/****************************************************************
 * Routine: Apply_Boundary_Flux_Corrections                     *
 *                                                              *
 * Apply flux corrections at boundaries of the solution         *
 * block to ensure that the scheme is conservative at           *
 * boundaries with mesh resolution changes.                     *
 *                                                              *
 ****************************************************************/
void Apply_Boundary_Flux_Corrections(Euler2D_Quad_Block &SolnBlk,
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

/****************************************************************
 * Routine: Apply_Boundary_Flux_Corrections_Multistage_Explicit *
 *                                                              *
 * Apply flux corrections at boundaries of the solution         *
 * block to ensure that the scheme is conservative at           *
 * boundaries with mesh resolution changes.                     *
 *                                                              *
 ****************************************************************/
void Apply_Boundary_Flux_Corrections_Multistage_Explicit(Euler2D_Quad_Block &SolnBlk,
                                                         const int i_stage,
                                                         const int n_stage,
	                                                 const double &CFL_Number,
                                                         const int Time_Integration_Type,
                                                         const int Reconstruction_Type,
                                                         const int Limiter_Type,
	                                                 const int Flux_Function_Type,
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
          SolnBlk.dUdt[i][SolnBlk.JCu][k_residual] -= 
             (CFL_Number*SolnBlk.dt[i][SolnBlk.JCu])*SolnBlk.FluxN[i]/
             SolnBlk.Grid.Cell[i][SolnBlk.JCu].A;
       } /* endfor */
    } /* endif */

    /* Correct the fluxes at the south boundary as required. */

    if (Number_Neighbours_South_Boundary == 2) {
       for ( i = SolnBlk.ICl ; i <= SolnBlk.ICu ; ++i ) {
          SolnBlk.dUdt[i][SolnBlk.JCl][k_residual] -= 
             (CFL_Number*SolnBlk.dt[i][SolnBlk.JCl])*SolnBlk.FluxS[i]/
             SolnBlk.Grid.Cell[i][SolnBlk.JCl].A;
       } /* endfor */
    } /* endif */

    /* Correct the fluxes at the east boundary as required. */

    if (Number_Neighbours_East_Boundary == 2) {
       for ( j= SolnBlk.JCl ; j <= SolnBlk.JCu ; ++j ) {
          SolnBlk.dUdt[SolnBlk.ICu][j][k_residual] -= 
             (CFL_Number*SolnBlk.dt[SolnBlk.ICu][j])*SolnBlk.FluxE[j]/
             SolnBlk.Grid.Cell[SolnBlk.ICu][j].A;
       } /* endfor */
    } /* endif */

    /* Correct the fluxes at the west boundary as required. */

    if (Number_Neighbours_West_Boundary == 2) {
       for ( j= SolnBlk.JCl ; j <= SolnBlk.JCu ; ++j ) {
          SolnBlk.dUdt[SolnBlk.ICl][j][k_residual] -= 
             (CFL_Number*SolnBlk.dt[SolnBlk.ICl][j])*SolnBlk.FluxW[j]/
             SolnBlk.Grid.Cell[SolnBlk.ICl][j].A;
       } /* endfor */
    } /* endif */

}

/********************************************************
 * Routine: dUdt_Residual_Evaluation                    *
 *                                                      *
 * This routine evaluates the residual for the specified*
 * solution block using a 2nd-order limited upwind      *
 * finite-volume spatial discretization scheme with     *
 * either the Godunov, Roe, Rusanov, HLLE, Linde, or    *
 * HLLC flux functions.                                 *
 * The residual is stored in dUdt[][][0].               *
 *                                                      *
 ********************************************************/
int dUdt_Residual_Evaluation(Euler2D_Quad_Block &SolnBlk,
			     Euler2D_Input_Parameters &Input_Parameters) {

  if (Input_Parameters.i_Reconstruction == RECONSTRUCTION_HIGH_ORDER){
    // calculate the high-order residual
    return SolnBlk.dUdt_Residual_Evaluation_HighOrder(Input_Parameters);
  } else {
    // calculate the 2nd-order residual
    return SolnBlk.dUdt_Residual_Evaluation(Input_Parameters);
  }

}

/********************************************************
 * Routine: dUdt_Multistage_Explicit                    *
 *                                                      *
 * This routine determines the solution residuals for a *
 * given stage of a variety of multi-stage explicit     *
 * time integration schemes for a given solution block. *
 *                                                      *
 ********************************************************/
int dUdt_Multistage_Explicit(Euler2D_Quad_Block &SolnBlk,
                             const int i_stage,
                             Euler2D_Input_Parameters &Input_Parameters) {

  if (Input_Parameters.i_Reconstruction == RECONSTRUCTION_HIGH_ORDER){
    // calculate the high-order residual
    return SolnBlk.dUdt_Multistage_Explicit_HighOrder(i_stage,Input_Parameters);
  } else {
    // calculate the 2nd-order residual
    return SolnBlk.dUdt_Multistage_Explicit(i_stage,Input_Parameters);
  }
    
}

/*********************************************************
 * Routine: Rotation_Matrix2                             *
 *                                                       *
 * This function returns either the rotation matrix, A,  *
 * or the inverse of A.                                  *
 *                                                       *
 * Note: A_matrix = 1 for returning A.                   *
 *                = 0 for returning inverse of A.        *
 *                                                       *
 *********************************************************/
DenseMatrix Rotation_Matrix2(Vector2D nface, int A_matrix) 
{
  double cos_angle = nface.x; 
  double sin_angle = nface.y;
    
  DenseMatrix mat(4,4);
  mat.identity();

  if (A_matrix) {
    // Rotation Matrix, A 
    mat(1,1) = cos_angle;
    mat(1,2) = sin_angle;
    mat(2,1) = -sin_angle;
    mat(2,2) = cos_angle;
  } else {
    // Rotation Matrix, Inv A 
    mat(1,1) = cos_angle;
    mat(1,2) = -sin_angle;
    mat(2,1) = sin_angle;
    mat(2,2) = cos_angle;
  } /* endif */

  return mat;

} /* End of Rotation_Matrix. */

/********************************************************
 * Routine: Residual_Jacobian                           *
 *                                                      *
 * This routine returns residual Jacobian matrix for the*
 * specified local solution block.                      *
 *                                                      *
 ********************************************************/
DenseMatrix Residual_Jacobian(const int i, 
			      const int j, 
			      Euler2D_Quad_Block &SolnBlk, 
			      const int cell_flag) {
  
  DenseMatrix dRdU(NUM_VAR_EULER2D, NUM_VAR_EULER2D),
              II(NUM_VAR_EULER2D, NUM_VAR_EULER2D), 
              A_N(NUM_VAR_EULER2D, NUM_VAR_EULER2D), 
              AI_N(NUM_VAR_EULER2D, NUM_VAR_EULER2D),
              A_S(NUM_VAR_EULER2D, NUM_VAR_EULER2D), 
              AI_S(NUM_VAR_EULER2D, NUM_VAR_EULER2D),
              A_E(NUM_VAR_EULER2D, NUM_VAR_EULER2D), 
              AI_E(NUM_VAR_EULER2D, NUM_VAR_EULER2D),
              A_W(NUM_VAR_EULER2D, NUM_VAR_EULER2D), 
              AI_W(NUM_VAR_EULER2D, NUM_VAR_EULER2D),
              dFdU_N(NUM_VAR_EULER2D, NUM_VAR_EULER2D),
              dFdU_S(NUM_VAR_EULER2D, NUM_VAR_EULER2D),
              dFdU_E(NUM_VAR_EULER2D, NUM_VAR_EULER2D),
	      dFdU_W(NUM_VAR_EULER2D, NUM_VAR_EULER2D);

  Vector2D lambdas_N, lambdas_S, lambdas_E, lambdas_W;
  Vector2D nface_N, nface_S, nface_E, nface_W;

  double alpha_N, gamma_N;
  double alpha_S, gamma_S;
  double alpha_W, gamma_W;
  double alpha_E, gamma_E;

  if (i < SolnBlk.ICl || i > SolnBlk.ICu ||
      j < SolnBlk.JCl || j > SolnBlk.JCu) {
    // GHOST CELL
    dRdU.zero();
    return dRdU;
    
  } else {
    // NON-GHOST CELL.
    
    dRdU.zero();
    II.identity();
    
    nface_N = SolnBlk.Grid.nfaceN(i, j);
    nface_S = SolnBlk.Grid.nfaceS(i, j);
    nface_E = SolnBlk.Grid.nfaceE(i, j);
    nface_W = SolnBlk.Grid.nfaceW(i, j);
    
    lambdas_N = HLLE_wavespeeds(W(SolnBlk.Uo[i][j]), 
				W(SolnBlk.Uo[i][j+1]), 
				nface_N);
    lambdas_S = HLLE_wavespeeds(W(SolnBlk.Uo[i][j]),
				W(SolnBlk.Uo[i][j-1]), 
				nface_S);
    lambdas_E = HLLE_wavespeeds(W(SolnBlk.Uo[i][j]), 
				W(SolnBlk.Uo[i+1][j]), 
				nface_E);
    lambdas_W = HLLE_wavespeeds(W(SolnBlk.Uo[i][j]),
				W(SolnBlk.Uo[i-1][j]),
				nface_W);

    dFdU_N.zero();
    dFdU(dFdU_N, Rotate(W(SolnBlk.Uo[i][j]), nface_N));
    A_N  = Rotation_Matrix2(nface_N, 1);
    AI_N = Rotation_Matrix2(nface_N, 0);
    if (lambdas_N.x >= ZERO) {
       dRdU = (-SolnBlk.Grid.lfaceN(i, j)/SolnBlk.Grid.area(i, j))*AI_N*dFdU_N*A_N;
    } else if (lambdas_N.y <= ZERO) {
       // Do nothing.
    } else {
       alpha_N = lambdas_N.y/(lambdas_N.y-lambdas_N.x);
       gamma_N = (lambdas_N.x*lambdas_N.y)/(lambdas_N.y-lambdas_N.x);
       dRdU = (-SolnBlk.Grid.lfaceN(i, j)/SolnBlk.Grid.area(i, j))*AI_N*(alpha_N*dFdU_N-gamma_N*II)*A_N;
    } /* endif */

    dFdU_S.zero();
    dFdU(dFdU_S, Rotate(W(SolnBlk.Uo[i][j]), nface_S));
    A_S  = Rotation_Matrix2(nface_S, 1);
    AI_S = Rotation_Matrix2(nface_S, 0);
    if (lambdas_S.x >= ZERO) {
       dRdU += (-SolnBlk.Grid.lfaceS(i, j)/SolnBlk.Grid.area(i, j))*AI_S*dFdU_S*A_S;
    } else if (lambdas_S.y <= ZERO) {
       // Do nothing.
    } else {
       alpha_S = lambdas_S.y/(lambdas_S.y-lambdas_S.x);
       gamma_S = (lambdas_S.x*lambdas_S.y)/(lambdas_S.y-lambdas_S.x);
       dRdU += (-SolnBlk.Grid.lfaceS(i, j)/SolnBlk.Grid.area(i, j))*AI_S*(alpha_S*dFdU_S-gamma_S*II)*A_S;
    } /* endif */

    dFdU_E.zero();
    dFdU(dFdU_E, Rotate(W(SolnBlk.Uo[i][j]), nface_E));
    A_E  = Rotation_Matrix2(nface_E, 1);
    AI_E = Rotation_Matrix2(nface_E, 0);
    if (lambdas_E.x >= ZERO) {
       dRdU += (-SolnBlk.Grid.lfaceE(i, j)/SolnBlk.Grid.area(i, j))*AI_E*dFdU_E*A_E;
    } else if (lambdas_E.y <= ZERO) {
       // Do nothing.
    } else {
       alpha_E = lambdas_E.y/(lambdas_E.y-lambdas_E.x);
       gamma_E = (lambdas_E.x*lambdas_E.y)/(lambdas_E.y-lambdas_E.x);
       dRdU += (-SolnBlk.Grid.lfaceE(i, j)/SolnBlk.Grid.area(i, j))*AI_E*(alpha_E*dFdU_E-gamma_E*II)*A_E;
    } /* endif */

    dFdU_W.zero();
    dFdU(dFdU_W, Rotate(W(SolnBlk.Uo[i][j]), nface_W));
    A_W  = Rotation_Matrix2(nface_W, 1);
    AI_W = Rotation_Matrix2(nface_W, 0);
    if (lambdas_W.x >= ZERO) {
       dRdU += (-SolnBlk.Grid.lfaceW(i, j)/SolnBlk.Grid.area(i, j))*AI_W*dFdU_W*A_W;
    } else if (lambdas_W.y <= ZERO) {
       // Do nothing.
    } else {
       alpha_W = lambdas_W.y/(lambdas_W.y-lambdas_W.x);
       gamma_W = (lambdas_W.x*lambdas_W.y)/(lambdas_W.y-lambdas_W.x);
       dRdU += (-SolnBlk.Grid.lfaceW(i, j)/SolnBlk.Grid.area(i, j))*AI_W*(alpha_W*dFdU_W-gamma_W*II)*A_W;
    } /* endif */
  } /* endif */
    
  return (dRdU);
  
} /* End of Residual_Jacobian. */

/********************************************************
 * Routine: Update_Solution_Multistage_Explicit         *
 *                                                      *
 * This routine updates solution states of the given    *
 * solution block for a variety of multi-stage explicit *
 * time integration schemes.                            *
 *                                                      *
 ********************************************************/
int Update_Solution_Multistage_Explicit(Euler2D_Quad_Block &SolnBlk,
                                        const int i_stage,
                                        Euler2D_Input_Parameters &Input_Parameters) {

    int i, j, k, l, k_residual, n_residual_reduction;
    double omega, residual_reduction_factor;

    // Memory for linear system solver.
    DenseMatrix dRdU(NUM_VAR_EULER2D,NUM_VAR_EULER2D),
                P_inv(NUM_VAR_EULER2D,NUM_VAR_EULER2D);
    DenseSystemLinEqs LinSys;
    Euler2D_cState dU_precon;

    /* Allocate memory for linear system solver. */

    LinSys.allocate(NUM_VAR_EULER2D);

    /* Perform update of solution variables for stage 
       i_stage of an N stage scheme. */

    /* Evaluate the time step fraction and residual storage location for the stage. */
    
    switch(Input_Parameters.i_Time_Integration) {
      case TIME_STEPPING_EXPLICIT_EULER :
        omega = Runge_Kutta(i_stage, Input_Parameters.N_Stage);
        k_residual = 0;
        break;
      case TIME_STEPPING_EXPLICIT_PREDICTOR_CORRECTOR :
        omega = Runge_Kutta(i_stage, Input_Parameters.N_Stage);
        k_residual = 0;
        break;
      case TIME_STEPPING_EXPLICIT_RUNGE_KUTTA :
        omega = Runge_Kutta(i_stage, Input_Parameters.N_Stage);
        k_residual = 0;
	if (Input_Parameters.N_Stage == 4) {
	   if (i_stage == 4) {
	      k_residual = 0;
           } else {
	      k_residual = i_stage - 1;
           } /* endif */
        } /* endif */
        break;
      case TIME_STEPPING_MULTISTAGE_OPTIMAL_SMOOTHING :
        omega = MultiStage_Optimally_Smoothing(i_stage, 
                                               Input_Parameters.N_Stage,
                                               Input_Parameters.i_Limiter);
        k_residual = 0;
        break;
      default:
        omega = Runge_Kutta(i_stage, Input_Parameters.N_Stage);
        k_residual = 0;
        break;
    } /* endswitch */
    
    /* Update solution variables for this stage. */
    
    for ( j  = SolnBlk.JCl ; j <= SolnBlk.JCu ; ++j ) {
       for ( i = SolnBlk.ICl ; i <= SolnBlk.ICu ; ++i ) {

	 if (Input_Parameters.Local_Time_Stepping == GLOBAL_TIME_STEPPING || 
             Input_Parameters.Local_Time_Stepping == SCALAR_LOCAL_TIME_STEPPING) {
           SolnBlk.U[i][j] = SolnBlk.Uo[i][j] + 
                             omega*SolnBlk.dUdt[i][j][k_residual];
         } else if (Input_Parameters.Local_Time_Stepping == LOW_MACH_NUMBER_WEISS_SMITH_PRECONDITIONER) {
	    // Apply Weiss-Smith low-Mach-number preconditioning to the residual.
            dU_precon = Euler2D_U_VACUUM;
            SolnBlk.W[i][j].P_U_WS_inv(P_inv);
            for ( k = 1 ; k <= NUM_VAR_EULER2D ; ++k ) {  
               for ( l = 1 ; l <= NUM_VAR_EULER2D ; ++l ) { 
                  dU_precon[k] += P_inv(k-1,l-1)*omega*SolnBlk.dUdt[i][j][k_residual][l];
               } /* endfor */
            } /* endfor */
            SolnBlk.U[i][j] = SolnBlk.Uo[i][j] + dU_precon;
	 } /* endif */

	 if (Input_Parameters.Local_Time_Stepping == GLOBAL_TIME_STEPPING && 
	     (SolnBlk.U[i][j].d   <= ZERO ||
	      SolnBlk.U[i][j].E   <= ZERO ||
	      SolnBlk.U[i][j].e() <= ZERO)) {
	     cout << "\n " << CFFC_Name() << " Euler2D ERROR: Negative Density and/or Energy: \n"
		  << " cell = (" << i << ", " << j << ") " 
		  << " X = " << SolnBlk.Grid.Cell[i][j].Xc << "\n"
		  << " Uo = " << SolnBlk.Uo[i][j] << "\n"
		  << " U = " << SolnBlk.U[i][j] << "\n"
		  << " dUdt = " << SolnBlk.dUdt[i][j][k_residual] << "\n";
	     return (i);

         } else if ((Input_Parameters.Local_Time_Stepping == SCALAR_LOCAL_TIME_STEPPING || 
                     Input_Parameters.Local_Time_Stepping == LOW_MACH_NUMBER_WEISS_SMITH_PRECONDITIONER) && 
	            (SolnBlk.U[i][j].d   <= ZERO ||
	             SolnBlk.U[i][j].E   <= ZERO ||
	             SolnBlk.U[i][j].e() <= ZERO)) {
	   residual_reduction_factor = ONE;

	   for (n_residual_reduction = 1; n_residual_reduction <= 10; ++n_residual_reduction) {
	     residual_reduction_factor = HALF*residual_reduction_factor;
	     SolnBlk.dt[i][j] = residual_reduction_factor*SolnBlk.dt[i][j];
	     SolnBlk.dUdt[i][j][k_residual] = residual_reduction_factor*SolnBlk.dUdt[i][j][k_residual];

             if (Input_Parameters.Local_Time_Stepping == SCALAR_LOCAL_TIME_STEPPING) {
                SolnBlk.U[i][j] = SolnBlk.Uo[i][j] + 
	                          omega*SolnBlk.dUdt[i][j][k_residual];
             } else if (Input_Parameters.Local_Time_Stepping == LOW_MACH_NUMBER_WEISS_SMITH_PRECONDITIONER) {
                dU_precon = residual_reduction_factor*dU_precon;
                SolnBlk.U[i][j] = SolnBlk.Uo[i][j] + dU_precon;
	     } /* endif */

	     if (SolnBlk.U[i][j].d   > ZERO &&
		 SolnBlk.U[i][j].E   > ZERO &&
		 SolnBlk.U[i][j].e() > ZERO ) {
		 break;
	     } /* endif */
	   } /* end for */

	   if (SolnBlk.U[i][j].d   <= ZERO ||
	       SolnBlk.U[i][j].E   <= ZERO ||
	       SolnBlk.U[i][j].e() <= ZERO ) {
	     cout << "\n " << CFFC_Name() << " Euler2D ERROR: Negative Density and/or Energy: \n"
		  << " cell = (" << i << ", " << j << ") " 
		  << " X = " << SolnBlk.Grid.Cell[i][j].Xc << "\n"
		  << " Uo = " << SolnBlk.Uo[i][j] << "\n"
		  << " U = " << SolnBlk.U[i][j] << "\n"
		  << " dUdt = " << SolnBlk.dUdt[i][j][k_residual] << "\n";
	     return (i);
	   } /* endif */
	 } /* end if */

	 if (Input_Parameters.Local_Time_Stepping == MATRIX_LOCAL_TIME_STEPPING) {

	   /* Evaluate Block Point-Jacobi Preconditioner dFdU */

	   dRdU = Residual_Jacobian(i,
				    j,
				    SolnBlk,
				    0);

	   /* Carry out matrix time stepping approach by setting up
	      linear system of equation for given cell */

	   LinSys.A = (-1.00*SolnBlk.dt[i][j])*dRdU;

	   for (k = 1; k <= NUM_VAR_EULER2D; ++k) {
	     LinSys.b(k-1) = omega*SolnBlk.dUdt[i][j][k_residual][k];
	   } /* endfor */

	   /* Solve system of equations using LU decomposition
	      Gaussian elimination procedure. */

	   LinSys.solve(LU_DECOMPOSITION);

	   /* Update the conserved solution variables. */

	   for (k = 1; k <= NUM_VAR_EULER2D; ++k) {
	     SolnBlk.U[i][j][k] = SolnBlk.Uo[i][j][k] + LinSys.x(k-1);
	   } /* endfor */

	   if (SolnBlk.U[i][j].d   <= ZERO ||
	       SolnBlk.U[i][j].E   <= ZERO ||
	       SolnBlk.U[i][j].e() <= ZERO) {	     
	     residual_reduction_factor = ONE;

	     for (n_residual_reduction = 1; n_residual_reduction <= 10; ++n_residual_reduction) {
	       residual_reduction_factor = HALF*residual_reduction_factor;
               if (n_residual_reduction == 10) residual_reduction_factor = ZERO;
	       SolnBlk.dt[i][j] = residual_reduction_factor*SolnBlk.dt[i][j];
	       SolnBlk.dUdt[i][j][k_residual] = residual_reduction_factor*SolnBlk.dUdt[i][j][k_residual];
	       LinSys.A = (-1.00*SolnBlk.dt[i][j])*dRdU;
	       for (k = 1; k <= NUM_VAR_EULER2D; ++k) {
		 LinSys.b(k-1) = omega*SolnBlk.dUdt[i][j][k_residual][k];
	       } /* endfor */
	       LinSys.solve(LU_DECOMPOSITION);
	       for (k = 1; k <= NUM_VAR_EULER2D; ++k) {
		 SolnBlk.U[i][j][k] = SolnBlk.Uo[i][j][k] + LinSys.x(k-1);
	       } /* endfor */
	       if (SolnBlk.U[i][j].d   > ZERO &&
		   SolnBlk.U[i][j].E   > ZERO &&
		   SolnBlk.U[i][j].e() > ZERO ) break;
	     } /* endfor */
	   } /* endif */

	   if (SolnBlk.U[i][j].d   <= ZERO ||
	       SolnBlk.U[i][j].E   <= ZERO ||
	       SolnBlk.U[i][j].e() <= ZERO ) {	     
	     cout << "\n " << CFFC_Name() << " Euler2D ERROR: Negative Density and/or Energy: \n"
		  << " cell = (" << i << ", " << j << ") " 
		  << " X = " << SolnBlk.Grid.Cell[i][j].Xc << "\n"
		  << " Uo = " << SolnBlk.Uo[i][j] << "\n"
		  << " U = " << SolnBlk.U[i][j] << "\n"
		  << " dUdt = " << SolnBlk.dUdt[i][j][k_residual] << "\n";
	     return (i);
	   } /* endif */

	 } /* endif */

	 SolnBlk.W[i][j] = W(SolnBlk.U[i][j]);

       } /* endfor */
    } /* endfor */

    /* Deallocate memory for linear system solver. */

    LinSys.deallocate();

    /* Solution successfully updated. */

    return (0);   

}

/**********************************************************************
 * Routine: Determine_Maximum_Surface_Pressure                        *
 *                                                                    *
 * This routine determines the maximum pressure on the surface of     *
 * of the wedge by comparing the known maximum to the pressures from  *
 * the current solution block.                                        *
 *                                                                    *
 **********************************************************************/   
int Determine_Maximum_Surface_Pressure(Euler2D_Quad_Block &SolnBlk,
				       double &maximum_pressure) {

  // Exit immediately if the current block is not on the wedge.
  if (SolnBlk.Grid.Cell[2][2].Xc.x < ZERO) return 0;
  if (SolnBlk.Grid.BCtypeS[2] != BC_REFLECTION) return 0;
 
  for (int i = SolnBlk.Grid.INl; i <= SolnBlk.Grid.INu; i++) {
    maximum_pressure = max(maximum_pressure,SolnBlk.Wn(i,SolnBlk.JCl).p);
  }

  // The maximum surface pressure for the current block has been found.
  return 0;

}

/**********************************************************************
 * Routine: Find_Height_of_Mach_Stem                                  *
 *                                                                    *
 * This routine finds the height of the Mach sterm in the Mach        *
 * reflection.                                                        *
 *                                                                    *
 **********************************************************************/   
int Determine_Mach_Stem_Height(Euler2D_Quad_Block &SolnBlk,
			       Euler2D_Input_Parameters &Input_Parameters,
			       const double &ms_pressure,
			       const int Number_of_Time_Steps,
			       const double &Time,
			       const int Block_Number,
			       const int Output_Title,
			       ostream &Out_File) {

  int error_flag, i, j;

  double Lx, Ly, L, lx, ly, l;
  double max_press_on_surface = ZERO;

  int foundIS = 0;  //foundIS = found Incident Shock
  int foundMS = 0;  //foundMS = found Mach Stem

  double right_bottom_cell_y, right_bottom_cell_x;
  double wedge_surface_y;
  double top_y = 0.5;
  double is_x1, is_x2, is_x, is_y1, is_y2, is_y;
  double ms_x1, ms_x2, ms_x, ms_y1, ms_y2, ms_y;
  //double mach_stem_pressure = 187300.0;

  Out_File << setprecision(18);
  if (Output_Title) {
    Out_File << "Location_of_Incident_Shock(IS)_along_y_direction(x,y_are_interpolated)" << "\n";
    Out_File << "j_index " << "x_location_of_IS " << "y_interpolated " << "y_uninterpolated " << "\n";
  } /* endif */

  /* Shock Centre Pressure ********************************
  //151987.5 = average(101325.0, 202650.0) for xi = 0.5
  //113990.625 = average(101325.0, 126656.25) for xi = 0.8
  ********************************************************/
  double Shock_Cntre_Pressure = (Input_Parameters.Wo.p + 101325.0)/TWO;
  //cout << Shock_Cntre_Pressure << endl;

  /* Only consider the blocks that are over to the wedge side  */
  if (SolnBlk.Grid.Cell[2][2].Xc.x >= 0) {
        
    /* Find the location of incident shock at different y */
    if ((SolnBlk.Wn(SolnBlk.Grid.INl, SolnBlk.Grid.JNu).p >= Shock_Cntre_Pressure) && 
	(SolnBlk.Wn(SolnBlk.Grid.INu, SolnBlk.Grid.JNu).p <= Shock_Cntre_Pressure))  {      
      
      //cout << "SolnBlk.Grid.Node[SolnBlk.INl][SolnBlk.JNu].X.x: " << SolnBlk.Grid.Node[SolnBlk.Grid.INl][SolnBlk.Grid.JNu].X.x 
      //   << " SolnBlk.Grid.Node[SolnBlk.INl][SolnBlk.JNu].X.y: " << SolnBlk.Grid.Node[SolnBlk.Grid.INl][SolnBlk.Grid.JNu].X.y 
      //   << endl;
      
      for (j = SolnBlk.Grid.JNu ; j >= SolnBlk.Grid.JNl ; --j) {
	foundIS = 0;
	i = SolnBlk.Grid.INu;
	
	while ((foundIS == 0) && (i != ONE)) {
	  i = i - 1;
	  if (SolnBlk.Wn(i, j).p >= Shock_Cntre_Pressure) {    
	    is_x1 = SolnBlk.Grid.Node[i+1][j].X.x;
	    is_x2 = SolnBlk.Grid.Node[i][j].X.x;
	    is_y1 = SolnBlk.Grid.Node[i+1][j].X.y;
	    is_y2 = SolnBlk.Grid.Node[i][j].X.y;
	    
	    is_x = (Shock_Cntre_Pressure - 101325.0) / (SolnBlk.Wn(i, j).p - 101325.0) * (is_x2 - is_x1) + is_x1;
	    is_y = (Shock_Cntre_Pressure - 101325.0) / (SolnBlk.Wn(i, j).p - 101325.0) * (is_y2 - is_y1) + is_y1;   
	    foundIS = 1; }
	}
	Out_File << j << " "
		     << is_x << " " << is_y << " " << SolnBlk.Grid.Node[i][j].X.y << "\n" << flush;
      }
      Out_File << "\n" << flush;
    }

    /* Find the location of Mach Stem on the wedge surface  */
    right_bottom_cell_y = SolnBlk.Grid.Cell[SolnBlk.ICu][0].Xc.y;
    right_bottom_cell_x = SolnBlk.Grid.Cell[SolnBlk.ICu][0].Xc.x;
    
    wedge_surface_y = right_bottom_cell_x * tan(Input_Parameters.Wedge_Angle*PI/180);
    j=2;
    if (right_bottom_cell_y <= wedge_surface_y) {
      
      if ((SolnBlk.Wn(SolnBlk.Grid.INl, SolnBlk.Grid.JNl).p >= ms_pressure) && 
	  (SolnBlk.Wn(SolnBlk.Grid.INu, SolnBlk.Grid.JNl).p <= ms_pressure))  {

	i = SolnBlk.Grid.INu;
	while (foundMS == 0) {
	  i = i - 1;
	  if (SolnBlk.Wn(i, j).p >= ms_pressure) {
	    ms_x1 = SolnBlk.Grid.Node[i+1][j].X.x;
	    ms_x2 = SolnBlk.Grid.Node[i][j].X.x;
	    ms_y1 = SolnBlk.Grid.Node[i+1][j].X.y;
	    ms_y2 = SolnBlk.Grid.Node[i][j].X.y;
	    
	    ms_x = (ms_pressure - SolnBlk.Wn(i+1, j).p) / (SolnBlk.Wn(i, j).p - SolnBlk.Wn(i+1, j).p) * (ms_x2 - ms_x1) + ms_x1;
	    ms_y = (ms_pressure - SolnBlk.Wn(i+1, j).p) / (SolnBlk.Wn(i, j).p - SolnBlk.Wn(i+1, j).p) * (ms_y2 - ms_y1) + ms_y1;
	    foundMS = 1;
	  }
	}
	
	Out_File << "\n" << flush;
	Out_File << "j: "<< j << "ms_x: " << ms_x << " ms_y: " << ms_y << "\n" << flush;
	
      }
      
      
    }
  }
  
  
}

/**********************************************************************
 * Routine: Output_Wedge_Solution_Distribution_Tecplot                *
 *                                                                    *
 **********************************************************************/
void Output_Wedge_Solution_Distribution_Tecplot(Euler2D_Quad_Block &SolnBlk,
						const double &Wedge_Angle,
						const int Number_of_Time_Steps,
						const double &Time,
						const int Block_Number,
						const int Output_Title,
						ostream &Out_File) {


  if (Output_Title) {
    Out_File << "TITLE = \"" << CFFC_Name() << ": 2D Euler Solution on the surface of the wedge, "
	     << "Time Step/Iteration Level = " << Number_of_Time_Steps
	     << ", Time = " << Time
	     << "\"" << "\n"
	     << "VARIABLES = \"l\" \\ \n"
	     << "\"y\" \\ \n"
	     << "\"x\" \\ \n"
	     << "\"rho\" \\ \n"
	     << "\"u\" \\ \n"
	     << "\"v\" \\ \n"
	     << "\"p\" \\ \n"
	     << "\"T\" \\ \n"
	     << "\"M\" \\ \n"
	     << "\"H\" \\ \n"
	     << "\"s\" \\ \n";
  }

  // Exit immediately if the current block is not on the wedge.
  if (SolnBlk.Grid.Cell[2][2].Xc.x < ZERO) return ;
  if (SolnBlk.Grid.BCtypeS[2] != BC_REFLECTION) return ;

  Euler2D_pState W_node;
 
  Out_File << "ZONE T =  \"Block Number = " << Block_Number
	   << "\" \\ \n"
	   << "I = " << SolnBlk.Grid.INu - SolnBlk.Grid.INl + 1 << " \\ \n"
	   << "J = " << 1 << " \\ \n"
	   << "F = POINT \n";

  Out_File << setprecision(14);
  for (int i = SolnBlk.Grid.INl; i <= SolnBlk.Grid.INu; i++) {
    W_node = SolnBlk.Wn(i,SolnBlk.JCl);
    Out_File << " " << SolnBlk.Grid.Node[i][SolnBlk.JCl].X.x/cos(Wedge_Angle*PI/180);
    Out_File << SolnBlk.Grid.Node[i][SolnBlk.JCl].X << W_node;
    Out_File.setf(ios::scientific);
    Out_File << " " << W_node.T() << " " << W_node.v.abs()/W_node.a() 
	     << " " << W_node.H() << " " << W_node.s() << "\n";
    Out_File.unsetf(ios::scientific);
  } /* endfor */
  Out_File << setprecision(6);  

}
