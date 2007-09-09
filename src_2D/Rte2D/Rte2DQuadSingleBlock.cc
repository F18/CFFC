/***************** Rte2DQuadSingleBlock.cc ***************************
   Single-Block Versions of Subroutines for 2D Rte 
   Multi-Block Quadrilateral Mesh  Solution Classes. 

   NOTES:
          
   - based on Euler2DQuadSingleBlock.cc
**********************************************************************/

/* Include 2D Rte quadrilateral mesh solution header file. */
#include "Rte2DQuad.h"

/**************************************************************************
 * Rte2D_Quad_Block -- Single Block External Subroutines.                 *
 **************************************************************************/

/********************************************************
 * Routine: Write_Solution_Block                        *
 *                                                      *
 * Writes the cell centred solution values of the       *
 * specified quadrilateral solution block to the        *
 * specified output stream for restart purposes.        *
 *                                                      *
 ********************************************************/
void Write_Solution_Block(Rte2D_Quad_Block &SolnBlk,
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
void Read_Solution_Block(Rte2D_Quad_Block &SolnBlk,
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
void Broadcast_Solution_Block(Rte2D_Quad_Block &SolnBlk) {

#ifdef _MPI_VERSION
    int NUM_VAR_RTE2D = SolnBlk.NumVar(); 
    int NUM_VAR_MEDIUM2D = Medium2D_State::NUM_VAR_MEDIUM2D; 
    int NUM_VAR = NUM_VAR_RTE2D + NUM_VAR_MEDIUM2D;
    int i, j, ni, nj, ng, nr, block_allocated, buffer_size;
    double *buffer;

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
    } /* endif */

    MPI::COMM_WORLD.Bcast(&ni, 1, MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(&nj, 1, MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(&ng, 1, MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(&nr,1,MPI::INT,0);
    MPI::COMM_WORLD.Bcast(&block_allocated, 1, MPI::INT, 0);

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
    MPI::COMM_WORLD.Bcast(&(SolnBlk.NorthWallTemp), 1, MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(SolnBlk.SouthWallTemp), 1, MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(SolnBlk.EastWallTemp), 1, MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(SolnBlk.WestWallTemp), 1, MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(SolnBlk.NorthWallEmiss), 1, MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(SolnBlk.SouthWallEmiss), 1, MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(SolnBlk.EastWallEmiss), 1, MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(SolnBlk.WestWallEmiss), 1, MPI::DOUBLE, 0);

    /* Broadcast the grid. */

    Broadcast_Quad_Block(SolnBlk.Grid);

    /* Set the grid 2D to quasi-3D scaling parameters*/

    if (block_allocated && !CFFC_Primary_MPI_Processor()) {
      SolnBlk.ScaleGridTo3D();
    } /* endif */

    /* Broadcast the solution state variables. */

    if (block_allocated) {
       ni = (SolnBlk.ICu+SolnBlk.Nghost) - (SolnBlk.ICl-SolnBlk.Nghost) + 1;
       nj = (SolnBlk.JCu+SolnBlk.Nghost) - (SolnBlk.JCl-SolnBlk.Nghost) + 1;
       buffer = new double[NUM_VAR*ni*nj];

       if (CFFC_Primary_MPI_Processor()) {
          buffer_size = 0;
          for (j  = SolnBlk.JCl-SolnBlk.Nghost ; j <= SolnBlk.JCu+SolnBlk.Nghost ; ++j ) {
              for ( i = SolnBlk.ICl-SolnBlk.Nghost ; i <= SolnBlk.ICu+SolnBlk.Nghost ; ++i ) {

		/* Rte2D Specific */
		for ( int k = 0; k < NUM_VAR_RTE2D; ++k) {	
		  buffer[buffer_size] = SolnBlk.U[i][j][k+1];
		  buffer_size = buffer_size + 1;
		}	      
		for ( int k = 0; k < NUM_VAR_MEDIUM2D; ++k) {	
		  buffer[buffer_size] = SolnBlk.M[i][j][k+1];
		  buffer_size = buffer_size + 1;
		}	      
		/* End Rte2D Specific */

              } /* endfor */
          } /* endfor */
       } /* endif */

       buffer_size = NUM_VAR*ni*nj;
       MPI::COMM_WORLD.Bcast(buffer, buffer_size, MPI::DOUBLE, 0);

       if (!CFFC_Primary_MPI_Processor()) {
          buffer_size = 0;
          for (j  = SolnBlk.JCl-SolnBlk.Nghost ; j <= SolnBlk.JCu+SolnBlk.Nghost ; ++j ) {
              for ( i = SolnBlk.ICl-SolnBlk.Nghost ; i <= SolnBlk.ICu+SolnBlk.Nghost ; ++i ) {
 
		/* Rte2D Specific */
		for ( int k = 0 ; k < NUM_VAR_RTE2D; ++ k) {	
		  SolnBlk.U[i][j][k+1]= buffer[buffer_size];
		  buffer_size = buffer_size + 1;
		}
		for ( int k = 0 ; k < NUM_VAR_MEDIUM2D; ++ k) {	
		  SolnBlk.M[i][j][k+1]= buffer[buffer_size];
		  buffer_size = buffer_size + 1;
		}
		/* End Rte2D Specific */

              } /* endfor */
           } /* endfor */
       } /* endif */

       delete []buffer; 
       buffer = NULL;

       ni = 1;
       nj = (SolnBlk.JCu+SolnBlk.Nghost) - (SolnBlk.JCl-SolnBlk.Nghost) + 1;
       buffer = new double[2*NUM_VAR_RTE2D*ni*nj];

       if (CFFC_Primary_MPI_Processor()) {
          buffer_size = 0;
          for (j  = SolnBlk.JCl-SolnBlk.Nghost ; j <= SolnBlk.JCu+SolnBlk.Nghost ; ++j ) {
 
	    /* Rte2D Specific */
	    for ( int k = 0 ; k < NUM_VAR_RTE2D; ++ k) {	
	      buffer[buffer_size]= SolnBlk.UoW[j][k+1];
	      buffer_size = buffer_size + 1;
	    }
	    for ( int k = 0 ; k < NUM_VAR_RTE2D; ++ k) {
	      buffer[buffer_size]= SolnBlk.UoE[j][k+1];  
	      buffer_size = buffer_size + 1;
	    }
	    /* End Rte2D Specific */

          } /* endfor */
       } /* endif */

       buffer_size = 2*NUM_VAR_RTE2D*ni*nj;
       MPI::COMM_WORLD.Bcast(buffer, buffer_size, MPI::DOUBLE, 0);

       if (!CFFC_Primary_MPI_Processor()) {
          buffer_size = 0;
          for (j  = SolnBlk.JCl-SolnBlk.Nghost ; j <= SolnBlk.JCu+SolnBlk.Nghost ; ++j ) {

	    /* Rte2D Specific */
	    for ( int k = 0 ; k < NUM_VAR_RTE2D; ++ k) {	
	      SolnBlk.UoW[j][k+1] = buffer[buffer_size];
	      buffer_size = buffer_size + 1;
	    }
	    for ( int k = 0 ; k < NUM_VAR_RTE2D; ++ k) {
	      SolnBlk.UoE[j][k+1] = buffer[buffer_size]; 
	      buffer_size = buffer_size + 1;
	    }
	    /* End Rte2D Specific */

          } /* endfor */
       } /* endif */

       delete []buffer; 
       buffer = NULL;

       ni = (SolnBlk.ICu+SolnBlk.Nghost) - (SolnBlk.ICl-SolnBlk.Nghost) + 1;
       nj = 1;
       buffer = new double[2*NUM_VAR_RTE2D*ni*nj];

       if (CFFC_Primary_MPI_Processor()) {
          buffer_size = 0;
          for (i  = SolnBlk.ICl-SolnBlk.Nghost ; i <= SolnBlk.ICu+SolnBlk.Nghost ; ++i ) {

	    /* Rte2D Specific */
	    for ( int k = 0 ; k < NUM_VAR_RTE2D; ++ k) {	
	      buffer[buffer_size]= SolnBlk.UoS[i][k+1];	  
	      buffer_size = buffer_size + 1;
	    }
	    for ( int k = 0 ; k < NUM_VAR_RTE2D; ++ k) {
	      buffer[buffer_size]= SolnBlk.UoN[i][k+1];
	      buffer_size = buffer_size + 1;
	    }
	    /* End Rte2D Specific */
	    
          } /* endfor */
       } /* endif */

       buffer_size = 2*NUM_VAR_RTE2D*ni*nj;
       MPI::COMM_WORLD.Bcast(buffer, buffer_size, MPI::DOUBLE, 0);

       if (!CFFC_Primary_MPI_Processor()) {
          buffer_size = 0;
          for (i  = SolnBlk.ICl-SolnBlk.Nghost ; i <= SolnBlk.ICu+SolnBlk.Nghost ; ++i ) {

	    /* Rte2D Specific */
	    for ( int k = 0 ; k < NUM_VAR_RTE2D; ++ k) {	
	      SolnBlk.UoS[i][k+1] = buffer[buffer_size]; 
	      buffer_size = buffer_size + 1;
	    }
	    for ( int k = 0 ; k < NUM_VAR_RTE2D; ++ k) {
	      SolnBlk.UoN[i][k+1] = buffer[buffer_size];
	     buffer_size = buffer_size + 1;
	    }
	    /* End Rte2D Specific */

          } /* endfor */
       } /* endif */

       delete []buffer; 
       buffer = NULL;
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
void Broadcast_Solution_Block(Rte2D_Quad_Block &SolnBlk,
                              MPI::Intracomm &Communicator, 
                              const int Source_CPU) {

    int NUM_VAR_RTE2D = SolnBlk.NumVar(); 
    int NUM_VAR_MEDIUM2D = Medium2D_State::NUM_VAR_MEDIUM2D; 
    int NUM_VAR = NUM_VAR_RTE2D + NUM_VAR_MEDIUM2D;
    int Source_Rank = 0;
    int i, j, ni, nj, ng, nr, block_allocated, buffer_size;
    double *buffer;

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
    } /* endif */

    Communicator.Bcast(&ni, 1, MPI::INT, Source_Rank);
    Communicator.Bcast(&nj, 1, MPI::INT, Source_Rank);
    Communicator.Bcast(&ng, 1, MPI::INT, Source_Rank);
    Communicator.Bcast(&nr,1,MPI::INT,Source_Rank);
    Communicator.Bcast(&block_allocated, 1, MPI::INT, Source_Rank);

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
    Communicator.Bcast(&(SolnBlk.NorthWallTemp), 1, MPI::DOUBLE, Source_Rank);
    Communicator.Bcast(&(SolnBlk.SouthWallTemp), 1, MPI::DOUBLE, Source_Rank);
    Communicator.Bcast(&(SolnBlk.EastWallTemp), 1, MPI::DOUBLE, Source_Rank);
    Communicator.Bcast(&(SolnBlk.WestWallTemp), 1, MPI::DOUBLE, Source_Rank);
    Communicator.Bcast(&(SolnBlk.NorthWallEmiss), 1, MPI::DOUBLE, Source_Rank);
    Communicator.Bcast(&(SolnBlk.SouthWallEmiss), 1, MPI::DOUBLE, Source_Rank);
    Communicator.Bcast(&(SolnBlk.EastWallEmiss), 1, MPI::DOUBLE, Source_Rank);
    Communicator.Bcast(&(SolnBlk.WestWallEmiss), 1, MPI::DOUBLE, Source_Rank);

    /* Broadcast the grid. */

    Broadcast_Quad_Block(SolnBlk.Grid, Communicator, Source_CPU);

    /* Set the grid 2D to quasi-3D scaling parameters*/

    if (block_allocated && CFFC_MPI::This_Processor_Number != Source_CPU) {
      SolnBlk.ScaleGridTo3D();
    } /* endif */


    /* Broadcast the solution state variables. */

    if (block_allocated) {
       ni = (SolnBlk.ICu+SolnBlk.Nghost) - (SolnBlk.ICl-SolnBlk.Nghost) + 1;
       nj = (SolnBlk.JCu+SolnBlk.Nghost) - (SolnBlk.JCl-SolnBlk.Nghost) + 1;
       buffer = new double[NUM_VAR*ni*nj];

       if (CFFC_MPI::This_Processor_Number == Source_CPU) {
          buffer_size = 0;
          for (j  = SolnBlk.JCl-SolnBlk.Nghost ; j <= SolnBlk.JCu+SolnBlk.Nghost ; ++j ) {
              for ( i = SolnBlk.ICl-SolnBlk.Nghost ; i <= SolnBlk.ICu+SolnBlk.Nghost ; ++i ) {

		/* Rte2D Specific */
		for ( int k=0; k<NUM_VAR_RTE2D; ++k) {	
		  buffer[buffer_size] = SolnBlk.U[i][j][k+1]; 
		  buffer_size = buffer_size + 1;
		}
		for ( int k=0; k<NUM_VAR_MEDIUM2D; ++k) {	
		  buffer[buffer_size] = SolnBlk.M[i][j][k+1]; 
		  buffer_size = buffer_size + 1;
		}
		/* End Rte2D Specific */

              } /* endfor */
          } /* endfor */
       } /* endif */

       buffer_size = NUM_VAR*ni*nj;
       Communicator.Bcast(buffer, buffer_size, MPI::DOUBLE, Source_Rank);

       if (!(CFFC_MPI::This_Processor_Number == Source_CPU)) {
          buffer_size = 0;
          for (j  = SolnBlk.JCl-SolnBlk.Nghost ; j <= SolnBlk.JCu+SolnBlk.Nghost ; ++j ) {
              for ( i = SolnBlk.ICl-SolnBlk.Nghost ; i <= SolnBlk.ICu+SolnBlk.Nghost ; ++i ) {

		/* Rte2D Specific */
		for ( int k=0; k<NUM_VAR_RTE2D; ++k) {	       
		  SolnBlk.U[i][j][k+1]= buffer[buffer_size];
		  buffer_size = buffer_size + 1;
		}
		for ( int k=0; k<NUM_VAR_MEDIUM2D; ++k) {	       
		  SolnBlk.M[i][j][k+1]= buffer[buffer_size];
		  buffer_size = buffer_size + 1;
		}
		/* End Rte2D Specific */

              } /* endfor */
           } /* endfor */
       } /* endif */

       delete []buffer; 
       buffer = NULL;

       ni = 1;
       nj = (SolnBlk.JCu+SolnBlk.Nghost) - (SolnBlk.JCl-SolnBlk.Nghost) + 1;
       buffer = new double[2*NUM_VAR_RTE2D*ni*nj];

       if (CFFC_MPI::This_Processor_Number == Source_CPU) {
          buffer_size = 0;
          for (j  = SolnBlk.JCl-SolnBlk.Nghost ; j <= SolnBlk.JCu+SolnBlk.Nghost ; ++j ) {

	    /* Rte2D Specific */
	    for ( int k=0; k<NUM_VAR_RTE2D; ++k) {	
	      buffer[buffer_size]= SolnBlk.UoW[j][k+1];	
	      buffer_size = buffer_size + 1;
	    }
	    for ( int k=0; k<NUM_VAR_RTE2D; ++k) {	
	      buffer[buffer_size]= SolnBlk.UoE[j][k+1];
	      buffer_size = buffer_size + 1;
	    }
	    /* End Rte2D Specific */

          } /* endfor */
       } /* endif */

       buffer_size = 2*NUM_VAR_RTE2D*ni*nj;
       Communicator.Bcast(buffer, buffer_size, MPI::DOUBLE, Source_Rank);

       if (!(CFFC_MPI::This_Processor_Number == Source_CPU)) {
          buffer_size = 0;
          for (j  = SolnBlk.JCl-SolnBlk.Nghost ; j <= SolnBlk.JCu+SolnBlk.Nghost ; ++j ) {
	    
	    /* Rte2D Specific */
	    for ( int k=0; k<NUM_VAR_RTE2D; ++k) {	
	      SolnBlk.UoW[j][k+1] = buffer[buffer_size];
	      buffer_size = buffer_size + 1;
	    }
	    for ( int k=0; k<NUM_VAR_RTE2D; ++k) {
	      SolnBlk.UoE[j][k+1] = buffer[buffer_size];
	      buffer_size = buffer_size + 1;
	    }
	    /* End Rte2D Specific */

          } /* endfor */
       } /* endif */

       delete []buffer; 
       buffer = NULL;

       ni = (SolnBlk.ICu+SolnBlk.Nghost) - (SolnBlk.ICl-SolnBlk.Nghost) + 1;
       nj = 1;
       buffer = new double[2*NUM_VAR_RTE2D*ni*nj];

       if (CFFC_MPI::This_Processor_Number == Source_CPU) {
          buffer_size = 0;
          for (i  = SolnBlk.ICl-SolnBlk.Nghost ; i <= SolnBlk.ICu+SolnBlk.Nghost ; ++i ) {

	    /* Rte2D Specific */
	    for ( int k=0; k<NUM_VAR_RTE2D; ++k) {    
	      buffer[buffer_size]= SolnBlk.UoS[i][k+1];   
	      buffer_size = buffer_size + 1;
	    }
	    for ( int k=0; k<NUM_VAR_RTE2D; ++k) {    
	      buffer[buffer_size]= SolnBlk.UoN[i][k+1];
	      buffer_size = buffer_size + 1;
	    } 
	    /* End Rte2D Specific */

          } /* endfor */
       } /* endif */

       buffer_size = 2*NUM_VAR_RTE2D*ni*nj;
       Communicator.Bcast(buffer, buffer_size, MPI::DOUBLE, Source_Rank);

       if (!(CFFC_MPI::This_Processor_Number == Source_CPU)) {
          buffer_size = 0;
          for (i  = SolnBlk.ICl-SolnBlk.Nghost ; i <= SolnBlk.ICu+SolnBlk.Nghost ; ++i ) {

	    /* Rte2D Specific */
	    for ( int k=0; k<NUM_VAR_RTE2D; ++k) {	
	      SolnBlk.UoS[i][k+1] = buffer[buffer_size]; 
	      buffer_size = buffer_size + 1;
	    }
	    for ( int k=0; k<NUM_VAR_RTE2D; ++k) {	
	      SolnBlk.UoN[i][k+1] = buffer[buffer_size]; 
	      buffer_size = buffer_size + 1;
	    }
	    /* End Rte2D Specific */

          } /* endfor */
       } /* endif */

       delete []buffer; 
       buffer = NULL;

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
void Copy_Solution_Block(Rte2D_Quad_Block &SolnBlk1,
                         Rte2D_Quad_Block &SolnBlk2) {

    int i, j, k;

    /* Allocate (re-allocate) memory for the solution
       of the quadrilateral solution block SolnBlk1 as necessary. */

    if (SolnBlk1.NCi != SolnBlk2.NCi || 
	SolnBlk1.NCj != SolnBlk2.NCj || 
	SolnBlk1.Nghost != SolnBlk2.Nghost) {
       if (SolnBlk1.U != NULL) SolnBlk1.deallocate();
       if (SolnBlk2.U != NULL) SolnBlk1.allocate(SolnBlk2.NCi-2*SolnBlk2.Nghost,
                                                 SolnBlk2.NCj-2*SolnBlk2.Nghost,
						 SolnBlk2.Nghost);
    } /* endif */

    /* Set the axisymmetric/planar flow indicator. */

    SolnBlk1.Axisymmetric = SolnBlk2.Axisymmetric;
    SolnBlk1.NorthWallTemp  = SolnBlk2.NorthWallTemp; 
    SolnBlk1.SouthWallTemp  = SolnBlk2.SouthWallTemp;  
    SolnBlk1.EastWallTemp   = SolnBlk2.EastWallTemp;  
    SolnBlk1.WestWallTemp   = SolnBlk2.WestWallTemp; 
    SolnBlk1.NorthWallEmiss = SolnBlk2.NorthWallEmiss;
    SolnBlk1.SouthWallEmiss = SolnBlk2.SouthWallEmiss; 
    SolnBlk1.EastWallEmiss  = SolnBlk2.EastWallEmiss;
    SolnBlk1.WestWallEmiss  = SolnBlk2.WestWallEmiss;

    /* Copy the grid of the second solution block
       to the first solution block. */

    Copy_Quad_Block(SolnBlk1.Grid, SolnBlk2.Grid);

    /* Set the grid 2D to quasi-3D scaling parameters*/
      
    SolnBlk1.ScaleGridTo3D();

    /* Copy the solution information from SolnBlk2 to SolnBlk1. */

    if (SolnBlk2.U != NULL) {
       for ( j  = SolnBlk1.JCl-SolnBlk1.Nghost ; j <= SolnBlk1.JCu+SolnBlk1.Nghost ; ++j ) {
          for ( i = SolnBlk1.ICl-SolnBlk1.Nghost ; i <= SolnBlk1.ICu+SolnBlk1.Nghost ; ++i ) {
             SolnBlk1.U[i][j] = SolnBlk2.U[i][j];
             SolnBlk1.M[i][j] = SolnBlk2.M[i][j];
             for ( k = 0 ; k <= NUMBER_OF_RESIDUAL_VECTORS_RTE2D-1 ; ++k ) {
	        SolnBlk1.dUdt[i][j][k] = SolnBlk2.dUdt[i][j][k];
             } /* endfor */
	     SolnBlk1.dUdx[i][j] = SolnBlk2.dUdx[i][j];
	     SolnBlk1.dUdy[i][j] = SolnBlk2.dUdy[i][j];
	     SolnBlk1.phi[i][j] = SolnBlk2.phi[i][j];
	     SolnBlk1.dUdpsi[i][j] = SolnBlk2.dUdpsi[i][j];
	     SolnBlk1.phi_psi[i][j] = SolnBlk2.phi_psi[i][j];
	     SolnBlk1.Uo[i][j] = SolnBlk2.Uo[i][j];
	     SolnBlk1.dt[i][j] = SolnBlk2.dt[i][j];
	     SolnBlk1.Sp[i][j]  = SolnBlk2.Sp[i][j]; 
	     SolnBlk1.SpN[i][j] = SolnBlk2.SpN[i][j]; 
	     SolnBlk1.SpS[i][j] = SolnBlk2.SpS[i][j]; 
	     SolnBlk1.SpE[i][j] = SolnBlk2.SpE[i][j]; 
	     SolnBlk1.SpW[i][j] = SolnBlk2.SpW[i][j]; 
          } /* endfor */
       } /* endfor */

       for (j  = SolnBlk1.JCl-SolnBlk1.Nghost ; j <= SolnBlk1.JCu+SolnBlk1.Nghost ; ++j ) {
	   SolnBlk1.UoW[j] = SolnBlk2.UoW[j];
           SolnBlk1.UoE[j] = SolnBlk2.UoE[j];
       } /* endfor */

       for ( i = SolnBlk1.ICl-SolnBlk1.Nghost ; i <= SolnBlk1.ICu+SolnBlk1.Nghost ; ++i ) {
           SolnBlk1.UoS[i] = SolnBlk2.UoS[i];
           SolnBlk1.UoN[i] = SolnBlk2.UoN[i];
       } /* endfor */
    } /* endif */

}

/********************************************************
 * Routine: Prolong_Solution_Block                      *
 *                                                      *
 * Prolongs the solution information of one of the      *
 * specified sectors of the original quadrilateral      *
 * solution block SolnBlk_Original to the refined       *
 * solution block SolnBlk_Fine.                         *
 *                                                      *
 ********************************************************/
int Prolong_Solution_Block(Rte2D_Quad_Block &SolnBlk_Fine,
		            Rte2D_Quad_Block &SolnBlk_Original,
                            const int Sector) {

    int i, j, i_min, i_max, j_min, j_max, mesh_refinement_permitted;
    double area_total_fine;
    Vector2D dX;

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
    } /* endif */

    if (mesh_refinement_permitted) {
   
    /* Set the axisymmetric/planar flow indicator for the fine solution block and
       copy some relevant information. */

      SolnBlk_Fine.Axisymmetric = SolnBlk_Original.Axisymmetric;
      SolnBlk_Fine.NorthWallTemp = SolnBlk_Original.NorthWallTemp; 
      SolnBlk_Fine.SouthWallTemp = SolnBlk_Original.SouthWallTemp;  
      SolnBlk_Fine.EastWallTemp = SolnBlk_Original.EastWallTemp;  
      SolnBlk_Fine.WestWallTemp = SolnBlk_Original.WestWallTemp; 
      SolnBlk_Fine.NorthWallEmiss = SolnBlk_Original.NorthWallEmiss;
      SolnBlk_Fine.SouthWallEmiss = SolnBlk_Original.SouthWallEmiss; 
      SolnBlk_Fine.EastWallEmiss = SolnBlk_Original.EastWallEmiss;
      SolnBlk_Fine.WestWallEmiss = SolnBlk_Original.WestWallEmiss;

      /* Set the grid 2D to quasi-3D scaling parameters*/
      
      SolnBlk_Fine.ScaleGridTo3D();

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
//                area_total_fine = SolnBlk_Fine.Grid.Cell[2*(i-i_min)+SolnBlk_Fine.ICl  ]
//                                                        [2*(j-j_min)+SolnBlk_Fine.JCl  ].A+
//                                  SolnBlk_Fine.Grid.Cell[2*(i-i_min)+SolnBlk_Fine.ICl+1]
//                                                        [2*(j-j_min)+SolnBlk_Fine.JCl  ].A+
//                                  SolnBlk_Fine.Grid.Cell[2*(i-i_min)+SolnBlk_Fine.ICl  ]
//                                                        [2*(j-j_min)+SolnBlk_Fine.JCl+1].A+
//                                  SolnBlk_Fine.Grid.Cell[2*(i-i_min)+SolnBlk_Fine.ICl+1]
//                                                        [2*(j-j_min)+SolnBlk_Fine.JCl+1].A;

     	       SolnBlk_Fine.U[2*(i-i_min)+SolnBlk_Fine.ICl  ]
                             [2*(j-j_min)+SolnBlk_Fine.JCl  ] 
                   = SolnBlk_Original.U[i][j];
//                    = (SolnBlk_Original.Grid.Cell[i][j].A/area_total_fine)*SolnBlk_Original.U[i][j];

     	       SolnBlk_Fine.U[2*(i-i_min)+SolnBlk_Fine.ICl+1]
                             [2*(j-j_min)+SolnBlk_Fine.JCl  ] 
                   = SolnBlk_Original.U[i][j];
//                    = (SolnBlk_Original.Grid.Cell[i][j].A/area_total_fine)*SolnBlk_Original.U[i][j];

     	       SolnBlk_Fine.U[2*(i-i_min)+SolnBlk_Fine.ICl  ]
                             [2*(j-j_min)+SolnBlk_Fine.JCl+1]
                   = SolnBlk_Original.U[i][j];
// 	       = (SolnBlk_Original.Grid.Cell[i][j].A/area_total_fine)*SolnBlk_Original.U[i][j];

     	       SolnBlk_Fine.U[2*(i-i_min)+SolnBlk_Fine.ICl+1]
                             [2*(j-j_min)+SolnBlk_Fine.JCl+1]
                   = SolnBlk_Original.U[i][j];
//                    = (SolnBlk_Original.Grid.Cell[i][j].A/area_total_fine)*SolnBlk_Original.U[i][j];

           } /* endfor */
       } /* endfor */

       // and the medium state
       PrescribeFields(SolnBlk_Fine);

      // boundary ref states
       for ( j  = j_min-SolnBlk_Original.Nghost/2; j <= j_max+SolnBlk_Original.Nghost/2 ; ++j ) {
           SolnBlk_Fine.UoW[2*(j-j_min)+SolnBlk_Fine.JCl  ]
              = SolnBlk_Original.UoW[j];
           SolnBlk_Fine.UoW[2*(j-j_min)+SolnBlk_Fine.JCl+1]
              = SolnBlk_Original.UoW[j];

           SolnBlk_Fine.UoE[2*(j-j_min)+SolnBlk_Fine.JCl  ]
              = SolnBlk_Original.UoE[j];
           SolnBlk_Fine.UoE[2*(j-j_min)+SolnBlk_Fine.JCl+1]
              = SolnBlk_Original.UoE[j];
       } /* endfor */

       for ( i = i_min-SolnBlk_Original.Nghost/2 ; i <= i_max+SolnBlk_Original.Nghost/2 ; ++i ) {
           SolnBlk_Fine.UoS[2*(i-i_min)+SolnBlk_Fine.ICl  ]
              = SolnBlk_Original.UoS[i];
           SolnBlk_Fine.UoS[2*(i-i_min)+SolnBlk_Fine.ICl+1]
              = SolnBlk_Original.UoS[i];

           SolnBlk_Fine.UoN[2*(i-i_min)+SolnBlk_Fine.ICl  ]
              = SolnBlk_Original.UoN[i];
           SolnBlk_Fine.UoN[2*(i-i_min)+SolnBlk_Fine.ICl+1]
              = SolnBlk_Original.UoN[i];
       } /* endfor */

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
int Restrict_Solution_Block(Rte2D_Quad_Block &SolnBlk_Coarse,
			    Rte2D_Quad_Block &SolnBlk_Original_SW,
			    Rte2D_Quad_Block &SolnBlk_Original_SE,
			    Rte2D_Quad_Block &SolnBlk_Original_NW,
			    Rte2D_Quad_Block &SolnBlk_Original_NE) {

    int i, j, i_coarse, j_coarse, mesh_coarsening_permitted;
 
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
    } /* endif */

    if (mesh_coarsening_permitted) {

    /* Set the axisymmetric/planar flow indicator for the fine solution block and
       copy some relevant information. */

      SolnBlk_Coarse.Axisymmetric = SolnBlk_Original_SW.Axisymmetric;
      SolnBlk_Coarse.NorthWallTemp = SolnBlk_Original_SW.NorthWallTemp; 
      SolnBlk_Coarse.SouthWallTemp = SolnBlk_Original_SW.SouthWallTemp;  
      SolnBlk_Coarse.EastWallTemp = SolnBlk_Original_SW.EastWallTemp;  
      SolnBlk_Coarse.WestWallTemp = SolnBlk_Original_SW.WestWallTemp; 
      SolnBlk_Coarse.NorthWallEmiss = SolnBlk_Original_SW.NorthWallEmiss;
      SolnBlk_Coarse.SouthWallEmiss = SolnBlk_Original_SW.SouthWallEmiss; 
      SolnBlk_Coarse.EastWallEmiss = SolnBlk_Original_SW.EastWallEmiss;
      SolnBlk_Coarse.WestWallEmiss = SolnBlk_Original_SW.WestWallEmiss;

      /* Set the grid 2D to quasi-3D scaling parameters*/
      
      SolnBlk_Coarse.ScaleGridTo3D();

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
						     SolnBlk_Original_SW.Sp[i  ][j  ]*
                                                     SolnBlk_Original_SW.U[i  ][j  ] +
                                                     SolnBlk_Original_SW.Grid.Cell[i+1][j  ].A*
						     SolnBlk_Original_SW.Sp[i+1][j  ]*
                                                     SolnBlk_Original_SW.U[i+1][j  ] + 
                                                     SolnBlk_Original_SW.Grid.Cell[i  ][j+1].A*
						     SolnBlk_Original_SW.Sp[i  ][j+1]*
                                                     SolnBlk_Original_SW.U[i  ][j+1] +
                                                     SolnBlk_Original_SW.Grid.Cell[i+1][j+1].A*
						     SolnBlk_Original_SW.Sp[i+1][j+1]*
                                                     SolnBlk_Original_SW.U[i+1][j+1]) /
                                                    (SolnBlk_Original_SW.Grid.Cell[i  ][j  ].A*
						     SolnBlk_Original_SW.Sp[i  ][j  ] +
						     SolnBlk_Original_SW.Grid.Cell[i+1][j  ].A*
						     SolnBlk_Original_SW.Sp[i+1][j  ] +
                                                     SolnBlk_Original_SW.Grid.Cell[i  ][j+1].A*
						     SolnBlk_Original_SW.Sp[i  ][j+1] +
                                                     SolnBlk_Original_SW.Grid.Cell[i+1][j+1].A*
						     SolnBlk_Original_SW.Sp[i+1][j+1]);
                                                    //SolnBlk_Coarse.Grid.Cell[i_coarse][j_coarse].A;
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
             SolnBlk_Coarse.U[i_coarse][j_coarse] = (SolnBlk_Original_SE.Grid.Cell[i  ][j  ].A*
						     SolnBlk_Original_SE.Sp[i  ][j  ]*
                                                     SolnBlk_Original_SE.U[i  ][j  ] +
                                                     SolnBlk_Original_SE.Grid.Cell[i+1][j  ].A*
						     SolnBlk_Original_SE.Sp[i+1][j  ]*
                                                     SolnBlk_Original_SE.U[i+1][j  ] + 
                                                     SolnBlk_Original_SE.Grid.Cell[i  ][j+1].A*
						     SolnBlk_Original_SE.Sp[i  ][j+1]*
                                                     SolnBlk_Original_SE.U[i  ][j+1] +
                                                     SolnBlk_Original_SE.Grid.Cell[i+1][j+1].A*
						     SolnBlk_Original_SE.Sp[i+1][j+1]*
                                                     SolnBlk_Original_SE.U[i+1][j+1]) /
 	                                            (SolnBlk_Original_SE.Grid.Cell[i  ][j  ].A*
						     SolnBlk_Original_SE.Sp[i  ][j  ] +
						     SolnBlk_Original_SE.Grid.Cell[i+1][j  ].A*
						     SolnBlk_Original_SE.Sp[i+1][j  ] +
                                                     SolnBlk_Original_SE.Grid.Cell[i  ][j+1].A*
						     SolnBlk_Original_SE.Sp[i  ][j+1] +
                                                     SolnBlk_Original_SE.Grid.Cell[i+1][j+1].A*
						     SolnBlk_Original_SE.Sp[i+1][j+1]);
                                                    //SolnBlk_Coarse.Grid.Cell[i_coarse][j_coarse].A;
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
             SolnBlk_Coarse.U[i_coarse][j_coarse] = (SolnBlk_Original_NW.Grid.Cell[i  ][j  ].A*
						     SolnBlk_Original_NW.Sp[i  ][j  ]*
                                                     SolnBlk_Original_NW.U[i  ][j  ] +
                                                     SolnBlk_Original_NW.Grid.Cell[i+1][j  ].A*
						     SolnBlk_Original_NW.Sp[i+1][j  ]*
                                                     SolnBlk_Original_NW.U[i+1][j  ] + 
                                                     SolnBlk_Original_NW.Grid.Cell[i  ][j+1].A*
						     SolnBlk_Original_NW.Sp[i  ][j+1]*
                                                     SolnBlk_Original_NW.U[i  ][j+1] +
                                                     SolnBlk_Original_NW.Grid.Cell[i+1][j+1].A*
						     SolnBlk_Original_NW.Sp[i+1][j+1]*
                                                     SolnBlk_Original_NW.U[i+1][j+1]) /
                                                    (SolnBlk_Original_NW.Grid.Cell[i  ][j  ].A*
						     SolnBlk_Original_NW.Sp[i  ][j  ] +
                                                     SolnBlk_Original_NW.Grid.Cell[i+1][j  ].A*
						     SolnBlk_Original_NW.Sp[i+1][j  ] +
                                                     SolnBlk_Original_NW.Grid.Cell[i  ][j+1].A*
						     SolnBlk_Original_NW.Sp[i  ][j+1] +
                                                     SolnBlk_Original_NW.Grid.Cell[i+1][j+1].A*
						     SolnBlk_Original_NW.Sp[i+1][j+1]);
                                                     //SolnBlk_Coarse.Grid.Cell[i_coarse][j_coarse].A;
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

      // North-East corner fine block:
      for ( j = SolnBlk_Original_NE.JCl; j <= SolnBlk_Original_NE.JCu ; j += 2 ) {
      	  for ( i = SolnBlk_Original_NE.ICl ; i <= SolnBlk_Original_NE.ICu ; i += 2 ) {
      	     i_coarse = (i-SolnBlk_Original_NE.ICl)/2+
                        (SolnBlk_Coarse.ICu-SolnBlk_Coarse.ICl+1)/2+SolnBlk_Coarse.ICl;
      	     j_coarse = (j-SolnBlk_Original_NE.JCl)/2+
                        (SolnBlk_Coarse.JCu-SolnBlk_Coarse.JCl+1)/2+SolnBlk_Coarse.JCl;
             SolnBlk_Coarse.U[i_coarse][j_coarse] = (SolnBlk_Original_NE.Grid.Cell[i  ][j  ].A*
						     SolnBlk_Original_NE.Sp[i  ][j  ]*
                                                     SolnBlk_Original_NE.U[i  ][j  ] +
                                                     SolnBlk_Original_NE.Grid.Cell[i+1][j  ].A*
						     SolnBlk_Original_NE.Sp[i+1][j  ]*
                                                     SolnBlk_Original_NE.U[i+1][j  ] + 
                                                     SolnBlk_Original_NE.Grid.Cell[i  ][j+1].A*
						     SolnBlk_Original_NE.Sp[i  ][j+1]*
                                                     SolnBlk_Original_NE.U[i  ][j+1] +
                                                     SolnBlk_Original_NE.Grid.Cell[i+1][j+1].A*
						     SolnBlk_Original_NE.Sp[i+1][j+1]*
                                                     SolnBlk_Original_NE.U[i+1][j+1]) /
                                                    (SolnBlk_Original_NE.Grid.Cell[i  ][j  ].A*
						     SolnBlk_Original_NE.Sp[i  ][j  ] +
                                                     SolnBlk_Original_NE.Grid.Cell[i+1][j  ].A*
						     SolnBlk_Original_NE.Sp[i+1][j  ] +
                                                     SolnBlk_Original_NE.Grid.Cell[i  ][j+1].A*
						     SolnBlk_Original_NE.Sp[i  ][j+1] +
                                                     SolnBlk_Original_NE.Grid.Cell[i+1][j+1].A*
						     SolnBlk_Original_NE.Sp[i+1][j+1]);
	                                             //SolnBlk_Coarse.Grid.Cell[i_coarse][j_coarse].A;
          } /* endfor */
      } /* endfor */

      // and the medium state
      PrescribeFields(SolnBlk_Coarse);

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
void Output_Tecplot(Rte2D_Quad_Block &SolnBlk,
		    Rte2D_Input_Parameters &IP,
                    const int Number_of_Time_Steps,
                    const double &Time,
                    const int Block_Number,
                    const int Output_Title,
	            ostream &Out_File) {

    int i, j;
    Rte2D_State U_node;
    Medium2D_State M_node;
    Vector2D q_node;

    /* Output node solution data. */

    Out_File << setprecision(14);
    if (Output_Title) {

       Out_File << "TITLE = \"" << CFFC_Name() << ": 2D Rte Solution, "
                << "Time Step/Iteration Level = " << Number_of_Time_Steps
                << ", Time = " << Time
                << "\"" << "\n"
	        << "VARIABLES = \"x\" \\ \n"
                << "\"y\" \\ \n"
		<<"\"G\" \\ \n"
	        <<"\"q.x\" \\ \n"
		<<"\"q.y\" \n";

       // blackbody
       for(int v=0 ;v<SolnBlk.U[0][0].Nband ;v++) 
	 Out_File << "\"Ib_"<<v+1<<"\" \\ \n";

       // apsorption coefficient
       for(int v=0 ;v<SolnBlk.U[0][0].Nband ;v++) 
	 Out_File << "\"kappa_"<<v+1<<"\" \\ \n";

       // scattering coefficient
       for(int v=0 ;v<SolnBlk.U[0][0].Nband ;v++) 
	 Out_File << "\"sigma_"<<v+1<<"\" \\ \n";
	   
       // intensity
       for(int v=0 ;v<SolnBlk.U[0][0].Nband ;v++) 
	 for(int m=0 ;m<SolnBlk.U[0][0].Npolar ;m++)
	   for(int l=0 ;l<SolnBlk.U[0][0].Nazim[m] ;l++)
	     Out_File <<"\"I -> v="<<v+1<<", m="<<m+1<<", l="<<l+1<<"\" \\ \n";

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
	   GetState(M_node, SolnBlk.Grid.Node[i][j].X);
	   q_node = U_node.q(M_node);
           Out_File << " " << SolnBlk.Grid.Node[i][j].X 
		    << " " << U_node.G(M_node)   // /(FOUR*PI*U_node.Ib)
		    << " " << q_node.x // /(PI*U_node.Ib)
		    << " " << q_node.y // /(PI*U_node.Ib)
		    << M_node
		    << U_node
		    << endl;
       } /* endfor */
    } /* endfor */
    Out_File << setprecision(6);
 
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
void Output_Cells_Tecplot(Rte2D_Quad_Block &SolnBlk,
			  Rte2D_Input_Parameters &IP,
                          const int Number_of_Time_Steps,
                          const double &Time,
                          const int Block_Number,
                          const int Output_Title,
	                  ostream &Out_File) {

    int i, j;
    int NUM_VAR_RTE2D = SolnBlk.NumVar(); 
    Vector2D qc;

    Out_File << setprecision(14);
    if (Output_Title) {
       Out_File << "TITLE = \"" << CFFC_Name() << ": 2D Rte Solution, "
                << "Time Step/Iteration Level = " << Number_of_Time_Steps
                << ", Time = " << Time
                << "\"" << "\n"
   	        << "VARIABLES = \"x\" \\ \n"
                << "\"y\" \\ \n"
		<<"\"G\" \\ \n"
	        <<"\"q.x\" \\ \n"
		<<"\"q.y\" \n";

       // Blackbody intensity
       for(int v=0 ;v<SolnBlk.U[0][0].Nband ;v++) 
	 Out_File << "\"Ib_"<<v+1<<"\" \\ \n";

       // absorbsion coefficient
       for(int v=0 ;v<SolnBlk.U[0][0].Nband ;v++) 
	 Out_File << "\"kappa_"<<v+1<<"\" \\ \n";

       // scattering coefficient
       for(int v=0 ;v<SolnBlk.U[0][0].Nband ;v++) 
	 Out_File << "\"sigma_"<<v+1<<"\" \\ \n";

       // intensity
       for(int v=0 ;v<SolnBlk.U[0][0].Nband ;v++) 
	 for(int m=0 ;m<SolnBlk.U[0][0].Npolar ;m++)
	   for(int l=0 ;l<SolnBlk.U[0][0].Nazim[m] ;l++)
	     Out_File <<"\"I -> v="<<v+1<<", m="<<m+1<<", l="<<l+1<<"\" \\ \n";
       
       // zone details
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

	 qc = SolnBlk.U[i][j].q(SolnBlk.M[i][j]);
	 Out_File << " " << SolnBlk.Grid.Cell[i][j].Xc
		  << " " << SolnBlk.U[i][j].G(SolnBlk.M[i][j])   // /(FOUR*PI*SolnBlk.U[i][j].Ib)
		  << " " << qc.x // /(PI*SolnBlk.U[i][j].Ib)
		  << " " << qc.y // /(PI*SolnBlk.U[i][j].Ib)
		  << SolnBlk.M[i][j]
		  << SolnBlk.U[i][j];
	 Out_File << "\n";
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
void Output_Nodes_Tecplot(Rte2D_Quad_Block &SolnBlk,
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
void Output_Gradients_Tecplot(Rte2D_Quad_Block &SolnBlk,
			      const int Number_of_Time_Steps,
			      const double &Time,
			      const int Block_Number,
			      const int Output_Title,
			      ostream &Out_File) {

  int NUM_VAR_RTE2D = SolnBlk.NumVar(); 

  Out_File << setprecision(14);
  if (Output_Title) {
    Out_File << "TITLE = \"" << CFFC_Name() << ": 2D Euler Solution, "
	     << "Time Step/Iteration Level = " << Number_of_Time_Steps
	     << ", Time = " << Time
	     << "\"" << "\n"
	     << "VARIABLES = \"x\" \\ \n"
 	     << "\"y\" \\ \n";

    // dI/dt
    for(int v=0 ;v<SolnBlk.U[0][0].Nband ;v++) 
      for(int m=0 ;m<SolnBlk.U[0][0].Npolar ;m++)
	for(int l=0 ;l<SolnBlk.U[0][0].Nazim[m] ;l++)
	  Out_File <<"\"dI/dt -> v="<<v+1<<", m="<<m+1<<", l="<<l+1<<"\" \\ \n";
    
    // dI/dx
    for(int v=0 ;v<SolnBlk.U[0][0].Nband ;v++) 
      for(int m=0 ;m<SolnBlk.U[0][0].Npolar ;m++)
	for(int l=0 ;l<SolnBlk.U[0][0].Nazim[m] ;l++)
	  Out_File <<"\"dI/dx -> v="<<v+1<<", m="<<m+1<<", l="<<l+1<<"\" \\ \n";
    
    // dI/dy
    for(int v=0 ;v<SolnBlk.U[0][0].Nband ;v++) 
      for(int m=0 ;m<SolnBlk.U[0][0].Npolar ;m++)
	for(int l=0 ;l<SolnBlk.U[0][0].Nazim[m] ;l++)
	  Out_File <<"\"dI/dy -> v="<<v+1<<", m="<<m+1<<", l="<<l+1<<"\" \\ \n";
    
    // phi
    for(int v=0 ;v<SolnBlk.U[0][0].Nband ;v++) 
      for(int m=0 ;m<SolnBlk.U[0][0].Npolar ;m++)
	for(int l=0 ;l<SolnBlk.U[0][0].Nazim[m] ;l++)
	  Out_File <<"\"phi -> v="<<v+1<<", m="<<m+1<<", l="<<l+1<<"\" \\ \n";
    
    // dI/dpsi
    for(int v=0 ;v<SolnBlk.U[0][0].Nband ;v++) 
      for(int m=0 ;m<SolnBlk.U[0][0].Npolar ;m++)
	for(int l=0 ;l<SolnBlk.U[0][0].Nazim[m] ;l++)
	  Out_File <<"\"dI/dpsi -> v="<<v+1<<", m="<<m+1<<", l="<<l+1<<"\" \\ \n";
    
    // phi_psi
    for(int v=0 ;v<SolnBlk.U[0][0].Nband ;v++) 
      for(int m=0 ;m<SolnBlk.U[0][0].Npolar ;m++)
	for(int l=0 ;l<SolnBlk.U[0][0].Nazim[m] ;l++)
	  Out_File <<"\"phi_psi -> v="<<v+1<<", m="<<m+1<<", l="<<l+1<<"\" \\ \n";

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
      Out_File << " " << SolnBlk.Grid.Cell[i][j].Xc;
      for (int k=1; k<=NUM_VAR_RTE2D; k++) Out_File << " " << SolnBlk.dUdt[i][j][0][k];
      for (int k=1; k<=NUM_VAR_RTE2D; k++) Out_File << " " << SolnBlk.dUdx[i][j][k];
      for (int k=1; k<=NUM_VAR_RTE2D; k++) Out_File << " " << SolnBlk.dUdy[i][j][k];
      for (int k=1; k<=NUM_VAR_RTE2D; k++) Out_File << " " << SolnBlk.phi[i][j][k];
      for (int k=1; k<=NUM_VAR_RTE2D; k++) Out_File << " " << SolnBlk.dUdpsi[i][j][k];
      for (int k=1; k<=NUM_VAR_RTE2D; k++) Out_File << " " << SolnBlk.phi_psi[i][j][k];
      Out_File << endl;
      Out_File.unsetf(ios::scientific);
    }
  }
  Out_File << setprecision(6);

}


/********************************************************
 * Routine: ICs                                         *
 *                                                      *
 * Assigns initial conditions and data to the           *
 * solution variables of the specified quadrilateral    *
 * solution block.                                      *
 *                                                      *
 ********************************************************/
void ICs(Rte2D_Quad_Block &SolnBlk,
         Rte2D_Input_Parameters &IP) {

    int i, j, k;
    Rte2D_State Wl, Wr;

    /* Assign the initial data for the IVP of interest. */

    switch(IP.i_ICs) {

      case IC_CONSTANT :
      case IC_UNIFORM :
        // Set the solution state to the initial state Wo[0].
        for (j  = SolnBlk.JCl-SolnBlk.Nghost ; j <= SolnBlk.JCu+SolnBlk.Nghost ; ++j ) {
            for ( i = SolnBlk.ICl-SolnBlk.Nghost ; i <= SolnBlk.ICu+SolnBlk.Nghost ; ++i ) {
	      SolnBlk.U[i][j] = IP.Uo;
            } /* endfor */
        } /* endfor */
        break;

      default :
        for (j  = SolnBlk.JCl-SolnBlk.Nghost ; j <= SolnBlk.JCu+SolnBlk.Nghost ; ++j ) {
            for ( i = SolnBlk.ICl-SolnBlk.Nghost ; i <= SolnBlk.ICu+SolnBlk.Nghost ; ++i ) {
	      SolnBlk.U[i][j] = IP.Uo;
            } /* endfor */
        } /* endfor */
        break;
    } /* endswitch */



    /* Set default values for the boundary condition reference states. */

    for (j = SolnBlk.JCl-SolnBlk.Nghost ; j <= SolnBlk.JCu+SolnBlk.Nghost ; ++j ) {
       if (j >= SolnBlk.JCl && j <= SolnBlk.JCu) {
	 SolnBlk.UoW[j] = SolnBlk.U[SolnBlk.ICl][j];
	 SolnBlk.UoE[j] = SolnBlk.U[SolnBlk.ICu][j];
       } else if (j < SolnBlk.JCl) {
	 SolnBlk.UoW[j] = SolnBlk.U[SolnBlk.ICl][SolnBlk.JCl];
	 SolnBlk.UoE[j] = SolnBlk.U[SolnBlk.ICu][SolnBlk.JCl];
       } else {
	 SolnBlk.UoW[j] = SolnBlk.U[SolnBlk.ICl][SolnBlk.JCu];
	 SolnBlk.UoE[j] = SolnBlk.U[SolnBlk.ICu][SolnBlk.JCu];
       } /* endif */
    } /* endfor */

    for ( i = SolnBlk.ICl-SolnBlk.Nghost ; i <= SolnBlk.ICu+SolnBlk.Nghost ; ++i ) {
       if (i >= SolnBlk.ICl && i <= SolnBlk.ICu) {
	 SolnBlk.UoS[i] = SolnBlk.U[i][SolnBlk.JCl];
	 SolnBlk.UoN[i] = SolnBlk.U[i][SolnBlk.JCu];
       } else if (i < SolnBlk.ICl) {
	 SolnBlk.UoS[i] = SolnBlk.U[SolnBlk.ICl][SolnBlk.JCl];
	 SolnBlk.UoN[i] = SolnBlk.U[SolnBlk.ICl][SolnBlk.JCu];
       } else {
	 SolnBlk.UoS[i] = SolnBlk.U[SolnBlk.ICu][SolnBlk.JCl];
	 SolnBlk.UoN[i] = SolnBlk.U[SolnBlk.ICu][SolnBlk.JCu];
       } /* endif */
    } /* endfor */

}

/********************************************************
 * Routine: PrescribeFields                             *
 *                                                      *
 * Assigns initial conditions and data to the           *
 * solution variables of the specified quadrilateral    *
 * solution block.  This is for the Medium state.       *
 *                                                      *
 ********************************************************/
void PrescribeFields(Rte2D_Quad_Block &SolnBlk) 
{
    for (int j  = SolnBlk.JCl-SolnBlk.Nghost ; j <= SolnBlk.JCu+SolnBlk.Nghost ; ++j ) {
      for (int i = SolnBlk.ICl-SolnBlk.Nghost ; i <= SolnBlk.ICu+SolnBlk.Nghost ; ++i ) {
	SolnBlk.M[i][j].GetState(SolnBlk.Grid.Cell[i][j].Xc);
      } // endfor
    } // endfor 
}

/********************************************************
 * Routine: BCs                                         *
 *                                                      *
 * Apply boundary conditions at boundaries of the       *
 * specified quadrilateral solution block.              *
 *                                                      *
 ********************************************************/
void BCs(Rte2D_Quad_Block &SolnBlk,
	 Rte2D_Input_Parameters &IP) {

    int i, j;
    double dx_norm;
    Vector2D dX;
    Rte2D_State dU, dUdx, Uwall;


    for ( j = SolnBlk.JCl-SolnBlk.Nghost ; j <= SolnBlk.JCu+SolnBlk.Nghost ; ++j ) {
      // Prescribe West boundary conditions.
      if ( (j >= SolnBlk.JCl && j <= SolnBlk.JCu) ||
           (j < SolnBlk.JCl && 
            (SolnBlk.Grid.BCtypeS[SolnBlk.ICl-1] == BC_NONE ||
             SolnBlk.Grid.BCtypeS[SolnBlk.ICl-1] == BC_PERIODIC ||
             SolnBlk.Grid.BCtypeS[SolnBlk.ICl-1] == BC_CONSTANT_EXTRAPOLATION ||
             SolnBlk.Grid.BCtypeS[SolnBlk.ICl-1] == BC_LINEAR_EXTRAPOLATION ||
	     SolnBlk.Grid.BCtypeS[SolnBlk.ICl-1] == BC_GRAY_WALL ) ) ||
           (j > SolnBlk.JCu && 
            (SolnBlk.Grid.BCtypeN[SolnBlk.ICl-1] == BC_NONE ||
             SolnBlk.Grid.BCtypeN[SolnBlk.ICl-1] == BC_PERIODIC ||
             SolnBlk.Grid.BCtypeN[SolnBlk.ICl-1] == BC_CONSTANT_EXTRAPOLATION ||
             SolnBlk.Grid.BCtypeN[SolnBlk.ICl-1] == BC_LINEAR_EXTRAPOLATION ||
	     SolnBlk.Grid.BCtypeN[SolnBlk.ICl-1] == BC_GRAY_WALL  ) ) ) {
        switch(SolnBlk.Grid.BCtypeW[j]) {
          case BC_NONE :
            break;
          case BC_FIXED :
	    for( int ghost = 1; ghost <= SolnBlk.Nghost; ghost++){
	      SolnBlk.U[SolnBlk.ICl-ghost][j] = SolnBlk.UoW[j]; 
	    }
            break;
          case BC_CONSTANT_EXTRAPOLATION :
	    for( int ghost = 1; ghost <= SolnBlk.Nghost; ghost++){
	      SolnBlk.U[SolnBlk.ICl-ghost][j] = SolnBlk.U[SolnBlk.ICl][j];
	    }
            break;
          case BC_LINEAR_EXTRAPOLATION :
            Linear_Reconstruction_LeastSquares_2(SolnBlk, 
                                                 SolnBlk.ICl, j, 
                                                 LIMITER_BARTH_JESPERSEN);
	    for( int ghost = 1; ghost <= SolnBlk.Nghost; ghost++){
	      dX = SolnBlk.Grid.Cell[SolnBlk.ICl-ghost][j].Xc -
		SolnBlk.Grid.Cell[SolnBlk.ICl][j].Xc;
	      SolnBlk.U[SolnBlk.ICl-ghost][j] = SolnBlk.U[SolnBlk.ICl][j] + 
		(SolnBlk.phi[SolnBlk.ICl][j]^SolnBlk.dUdx[SolnBlk.ICl][j])*dX.x +
		(SolnBlk.phi[SolnBlk.ICl][j]^SolnBlk.dUdy[SolnBlk.ICl][j])*dX.y;
	    }	      
            break;
          case BC_REFLECTION :
	    for( int ghost = 1; ghost <= SolnBlk.Nghost; ghost++){
	      SolnBlk.U[SolnBlk.ICl-ghost][j] = Reflect(SolnBlk.U[SolnBlk.ICl+ghost-1][j],
							SolnBlk.Grid.nfaceW(SolnBlk.ICl,j));
	    }
           break;
          case BC_PERIODIC :
	    for( int ghost = 1; ghost <= SolnBlk.Nghost; ghost++){
	      SolnBlk.U[SolnBlk.ICl-ghost][j] = SolnBlk.U[SolnBlk.ICu-ghost][j];
	    }
            break;
          case BC_GRAY_WALL :
            Linear_Reconstruction_LeastSquares_2(SolnBlk, 
                                                 SolnBlk.ICl, j, 
                                                 LIMITER_BARTH_JESPERSEN);
	    dX = SolnBlk.Grid.xfaceW(SolnBlk.ICl, j)-
	         SolnBlk.Grid.Cell[SolnBlk.ICl][j].Xc;
	    Uwall = SolnBlk.U[SolnBlk.ICl][j] + 
		(SolnBlk.phi[SolnBlk.ICl][j]^SolnBlk.dUdx[SolnBlk.ICl][j])*dX.x +
		(SolnBlk.phi[SolnBlk.ICl][j]^SolnBlk.dUdy[SolnBlk.ICl][j])*dX.y;
	    Uwall = Gray_Wall(Uwall, SolnBlk.M[0][0], 
			      SolnBlk.Grid.nfaceW(SolnBlk.ICl, j), 
			      SolnBlk.WestWallTemp, SolnBlk.WestWallEmiss );
	    dx_norm = dX*SolnBlk.Grid.nfaceW(SolnBlk.ICl, j);
	    dU = Uwall-SolnBlk.U[SolnBlk.ICl][j];
	    dUdx = dU/dx_norm;
	    for( int ghost = 1; ghost <= SolnBlk.Nghost; ghost++){
	      dX = SolnBlk.Grid.Cell[SolnBlk.ICl-ghost][j].Xc -
		   SolnBlk.Grid.Cell[SolnBlk.ICl][j].Xc;
	      dx_norm = dX*SolnBlk.Grid.nfaceW(SolnBlk.ICl, j);
	      SolnBlk.U[SolnBlk.ICl-ghost][j] = SolnBlk.U[SolnBlk.ICl][j] + dUdx*dx_norm;
	    }
// 	    SolnBlk.UoW[j] = Gray_Wall(SolnBlk.U[SolnBlk.ICl][j], SolnBlk.Grid.nfaceW(SolnBlk.ICl, j), 
// 				       SolnBlk.WestWallTemp, SolnBlk.WestWallEmiss );
// 	    for( int ghost = 1; ghost <= SolnBlk.Nghost; ghost++){
// 	      SolnBlk.U[SolnBlk.ICl-ghost][j] = SolnBlk.UoW[j]; 
// 	    }
            break;
          default:
	    for( int ghost = 1; ghost <= SolnBlk.Nghost; ghost++){
	      SolnBlk.U[SolnBlk.ICl-ghost][j] = SolnBlk.U[SolnBlk.ICl][j];
	    }
            break;
        } /* endswitch */
      } /* endif */


      // Prescribe East boundary conditions.
      if ( (j >= SolnBlk.JCl && j <= SolnBlk.JCu) ||
           (j < SolnBlk.JCl && 
            (SolnBlk.Grid.BCtypeS[SolnBlk.ICu+1] == BC_NONE ||
             SolnBlk.Grid.BCtypeS[SolnBlk.ICu+1] == BC_PERIODIC ||
             SolnBlk.Grid.BCtypeS[SolnBlk.ICu+1] == BC_CONSTANT_EXTRAPOLATION ||
             SolnBlk.Grid.BCtypeS[SolnBlk.ICu+1] == BC_LINEAR_EXTRAPOLATION ||
	     SolnBlk.Grid.BCtypeS[SolnBlk.ICu+1] == BC_GRAY_WALL ) ) ||
           (j > SolnBlk.JCu && 
            (SolnBlk.Grid.BCtypeN[SolnBlk.ICu+1] == BC_NONE ||
             SolnBlk.Grid.BCtypeN[SolnBlk.ICu+1] == BC_PERIODIC ||
             SolnBlk.Grid.BCtypeN[SolnBlk.ICu+1] == BC_CONSTANT_EXTRAPOLATION ||
             SolnBlk.Grid.BCtypeN[SolnBlk.ICu+1] == BC_LINEAR_EXTRAPOLATION ||
	     SolnBlk.Grid.BCtypeN[SolnBlk.ICu+1] == BC_GRAY_WALL ) ) ) {
        switch(SolnBlk.Grid.BCtypeE[j]) {
          case BC_NONE :
            break;
          case BC_FIXED :
	    for( int ghost = 1; ghost <= SolnBlk.Nghost; ghost++){
	      SolnBlk.U[SolnBlk.ICu+ghost][j] = SolnBlk.UoE[j];
	    }
            break;
          case BC_CONSTANT_EXTRAPOLATION : 
	    for( int ghost = 1; ghost <= SolnBlk.Nghost; ghost++){
	      SolnBlk.U[SolnBlk.ICu+ghost][j] = SolnBlk.U[SolnBlk.ICu][j];
	    }
            break;
          case BC_LINEAR_EXTRAPOLATION :
	      Linear_Reconstruction_LeastSquares_2(SolnBlk, 
						   SolnBlk.ICu, j, 
						   LIMITER_BARTH_JESPERSEN);
	    for( int ghost = 1; ghost <= SolnBlk.Nghost; ghost++){
	      dX = SolnBlk.Grid.Cell[SolnBlk.ICu+ghost][j].Xc -
		SolnBlk.Grid.Cell[SolnBlk.ICu][j].Xc;
	      SolnBlk.U[SolnBlk.ICu+ghost][j] = SolnBlk.U[SolnBlk.ICu][j] + 
		(SolnBlk.phi[SolnBlk.ICu][j]^SolnBlk.dUdx[SolnBlk.ICu][j])*dX.x +
		(SolnBlk.phi[SolnBlk.ICu][j]^SolnBlk.dUdy[SolnBlk.ICu][j])*dX.y;
	    }	    
            break;
          case BC_REFLECTION :
	    for( int ghost = 1; ghost <= SolnBlk.Nghost; ghost++){
	      SolnBlk.U[SolnBlk.ICu+ghost][j] = Reflect(SolnBlk.U[SolnBlk.ICu-ghost+1][j],
							SolnBlk.Grid.nfaceE(SolnBlk.ICu,j));
	    }       
            break;
          case BC_PERIODIC :
	    for( int ghost = 1; ghost <= SolnBlk.Nghost; ghost++){
	      SolnBlk.U[SolnBlk.ICu+ghost][j] = SolnBlk.U[SolnBlk.ICl+ghost][j];
	    }
            break;
          case BC_GRAY_WALL :
            Linear_Reconstruction_LeastSquares_2(SolnBlk, 
                                                 SolnBlk.ICu, j, 
                                                 LIMITER_BARTH_JESPERSEN);
	    dX = SolnBlk.Grid.xfaceE(SolnBlk.ICu, j)-
	         SolnBlk.Grid.Cell[SolnBlk.ICu][j].Xc;
	    Uwall = SolnBlk.U[SolnBlk.ICu][j] + 
		(SolnBlk.phi[SolnBlk.ICu][j]^SolnBlk.dUdx[SolnBlk.ICu][j])*dX.x +
		(SolnBlk.phi[SolnBlk.ICu][j]^SolnBlk.dUdy[SolnBlk.ICu][j])*dX.y;
	    Uwall = Gray_Wall(Uwall, SolnBlk.M[0][0], 
			      SolnBlk.Grid.nfaceE(SolnBlk.ICu, j), 
			      SolnBlk.EastWallTemp, SolnBlk.EastWallEmiss );
	    dx_norm = dX*SolnBlk.Grid.nfaceE(SolnBlk.ICu, j);
	    dU = Uwall-SolnBlk.U[SolnBlk.ICu][j];
	    dUdx = dU/dx_norm;
	    for( int ghost = 1; ghost <= SolnBlk.Nghost; ghost++){
	      dX = SolnBlk.Grid.Cell[SolnBlk.ICu+ghost][j].Xc -
		   SolnBlk.Grid.Cell[SolnBlk.ICu][j].Xc;
	      dx_norm = dX*SolnBlk.Grid.nfaceE(SolnBlk.ICu, j);
	      SolnBlk.U[SolnBlk.ICu+ghost][j] = SolnBlk.U[SolnBlk.ICu][j] + dUdx*dx_norm;
	    }
// 	    SolnBlk.UoE[j] = Gray_Wall(SolnBlk.U[SolnBlk.ICu][j], SolnBlk.Grid.nfaceE(SolnBlk.ICu, j), 
// 				       SolnBlk.EastWallTemp, SolnBlk.EastWallEmiss );
// 	    for( int ghost = 1; ghost <= SolnBlk.Nghost; ghost++){
// 	      SolnBlk.U[SolnBlk.ICu+ghost][j] = SolnBlk.UoE[j];
// 	    }
            break;

	  default:
	    for( int ghost = 1; ghost <= SolnBlk.Nghost; ghost++){
	      SolnBlk.U[SolnBlk.ICu+ghost][j] = SolnBlk.U[SolnBlk.ICu][j];
	    }
            break;
        } /* endswitch */
      } /* endif */

    } /* endfor */

    for ( i = SolnBlk.ICl-SolnBlk.Nghost ; i <= SolnBlk.ICu+SolnBlk.Nghost ; ++i ) {

        // Prescribe South boundary conditions.
        switch(SolnBlk.Grid.BCtypeS[i]) {
          case BC_NONE :
            break;
          case BC_FIXED :
	    for( int ghost = 1; ghost <= SolnBlk.Nghost; ghost++){
	      SolnBlk.U[i][SolnBlk.JCl-ghost] = SolnBlk.UoS[i];
	    }
	    break;
          case BC_CONSTANT_EXTRAPOLATION :
	    for( int ghost = 1; ghost <= SolnBlk.Nghost; ghost++){
	      SolnBlk.U[i][SolnBlk.JCl-ghost] = SolnBlk.U[i][SolnBlk.JCl];
	    }
            break;
          case BC_LINEAR_EXTRAPOLATION :
            Linear_Reconstruction_LeastSquares_2(SolnBlk, 
                                                 i, SolnBlk.JCl, 
                                                 LIMITER_BARTH_JESPERSEN);
	    for( int ghost = 1; ghost <= SolnBlk.Nghost; ghost++){
	      dX = SolnBlk.Grid.Cell[i][SolnBlk.JCl-ghost].Xc -
		SolnBlk.Grid.Cell[i][SolnBlk.JCl].Xc;
	      SolnBlk.U[i][SolnBlk.JCl-ghost] = SolnBlk.U[i][SolnBlk.JCl] + 
		(SolnBlk.phi[i][SolnBlk.JCl]^SolnBlk.dUdx[i][SolnBlk.JCl])*dX.x +
		(SolnBlk.phi[i][SolnBlk.JCl]^SolnBlk.dUdy[i][SolnBlk.JCl])*dX.y;
	    }
	    break;
          case BC_REFLECTION :
	    for( int ghost = 1; ghost <= SolnBlk.Nghost; ghost++){
	      SolnBlk.U[i][SolnBlk.JCl-ghost] = Reflect(SolnBlk.U[i][SolnBlk.JCl + ghost-1],
							SolnBlk.Grid.nfaceS(i, SolnBlk.JCl));
	    }
          break;
          case BC_PERIODIC : 
	    for( int ghost = 1; ghost <= SolnBlk.Nghost; ghost++){
	      SolnBlk.U[i][SolnBlk.JCl-ghost] = SolnBlk.U[i][SolnBlk.JCu-ghost];
	    }
            break;
          case BC_GRAY_WALL :
            Linear_Reconstruction_LeastSquares_2(SolnBlk, 
                                                 i, SolnBlk.JCl, 
                                                 LIMITER_BARTH_JESPERSEN);
	    dX = SolnBlk.Grid.xfaceS(i,SolnBlk.JCl) -
	         SolnBlk.Grid.Cell[i][SolnBlk.JCl].Xc;
	    Uwall = SolnBlk.U[i][SolnBlk.JCl] + 
		(SolnBlk.phi[i][SolnBlk.JCl]^SolnBlk.dUdx[i][SolnBlk.JCl])*dX.x +
		(SolnBlk.phi[i][SolnBlk.JCl]^SolnBlk.dUdy[i][SolnBlk.JCl])*dX.y;
	    Uwall = Gray_Wall(Uwall, SolnBlk.M[0][0], 
			      SolnBlk.Grid.nfaceS(i,SolnBlk.JCl), 
			      SolnBlk.SouthWallTemp, SolnBlk.SouthWallEmiss );
	    dx_norm = dX*SolnBlk.Grid.nfaceS(i,SolnBlk.JCl);
	    dU = Uwall-SolnBlk.U[i][SolnBlk.JCl];
	    dUdx = dU/dx_norm;
	    for( int ghost = 1; ghost <= SolnBlk.Nghost; ghost++){
	      dX = SolnBlk.Grid.Cell[i][SolnBlk.JCl-ghost].Xc -
		   SolnBlk.Grid.Cell[i][SolnBlk.JCl].Xc;
	      dx_norm = dX*SolnBlk.Grid.nfaceS(i, SolnBlk.JCl);
	      SolnBlk.U[i][SolnBlk.JCl-ghost] = SolnBlk.U[i][SolnBlk.JCl] + dUdx*dx_norm;
	    }
// 	    SolnBlk.UoS[i] = Gray_Wall(SolnBlk.U[i][SolnBlk.JCl], SolnBlk.Grid.nfaceS(i,SolnBlk.JCl), 
// 				       SolnBlk.SouthWallTemp, SolnBlk.SouthWallEmiss );
// 	    for( int ghost = 1; ghost <= SolnBlk.Nghost; ghost++){
// 	      SolnBlk.U[i][SolnBlk.JCl-ghost] = SolnBlk.UoS[i];
// 	    }
            break;

          default:
	    for( int ghost = 1; ghost <= SolnBlk.Nghost; ghost++){
	      SolnBlk.U[i][SolnBlk.JCl-ghost] = SolnBlk.U[i][SolnBlk.JCl];
	    }
           break;
        } /* endswitch */

        // Prescribe North boundary conditions.
        switch(SolnBlk.Grid.BCtypeN[i]) {
          case BC_NONE :
            break;
          case BC_FIXED :
	    for( int ghost = 1; ghost <= SolnBlk.Nghost; ghost++){
	      SolnBlk.U[i][SolnBlk.JCu+ghost] = SolnBlk.UoN[i];
	    }
            break;
          case BC_CONSTANT_EXTRAPOLATION :
	    for( int ghost = 1; ghost <= SolnBlk.Nghost; ghost++){
	      SolnBlk.U[i][SolnBlk.JCu+ghost] = SolnBlk.U[i][SolnBlk.JCu];
	    }
            break;
          case BC_LINEAR_EXTRAPOLATION :
            Linear_Reconstruction_LeastSquares_2(SolnBlk, 
                                                 i, SolnBlk.JCu, 
                                                 LIMITER_BARTH_JESPERSEN);
	    for( int ghost = 1; ghost <= SolnBlk.Nghost; ghost++){
	      dX = SolnBlk.Grid.Cell[i][SolnBlk.JCu+ghost].Xc -
		SolnBlk.Grid.Cell[i][SolnBlk.JCu].Xc;
	      SolnBlk.U[i][SolnBlk.JCu+ghost] = SolnBlk.U[i][SolnBlk.JCu] + 
		(SolnBlk.phi[i][SolnBlk.JCu]^SolnBlk.dUdx[i][SolnBlk.JCu])*dX.x +
		(SolnBlk.phi[i][SolnBlk.JCu]^SolnBlk.dUdy[i][SolnBlk.JCu])*dX.y;
	    }
            break;
          case BC_REFLECTION :
	    for( int ghost = 1; ghost <= SolnBlk.Nghost; ghost++){
	      SolnBlk.U[i][SolnBlk.JCu+ghost] = Reflect(SolnBlk.U[i][SolnBlk.JCu-ghost+1],
							SolnBlk.Grid.nfaceN(i, SolnBlk.JCu));
	    }
            break;
          case BC_PERIODIC :
	    for( int ghost = 1; ghost <= SolnBlk.Nghost; ghost++){
	      SolnBlk.U[i][SolnBlk.JCu+ghost] = SolnBlk.U[i][SolnBlk.JCl+ghost];
	    }
            break;
          case BC_GRAY_WALL :
            Linear_Reconstruction_LeastSquares_2(SolnBlk, 
                                                 i, SolnBlk.JCu, 
                                                 LIMITER_BARTH_JESPERSEN);
	    dX = SolnBlk.Grid.xfaceN(i,SolnBlk.JCu) -
	         SolnBlk.Grid.Cell[i][SolnBlk.JCu].Xc;
	    Uwall = SolnBlk.U[i][SolnBlk.JCu] + 
		(SolnBlk.phi[i][SolnBlk.JCu]^SolnBlk.dUdx[i][SolnBlk.JCu])*dX.x +
		(SolnBlk.phi[i][SolnBlk.JCu]^SolnBlk.dUdy[i][SolnBlk.JCu])*dX.y;
	    Uwall = Gray_Wall(Uwall, SolnBlk.M[0][0],
			      SolnBlk.Grid.nfaceN(i,SolnBlk.JCu), 
			      SolnBlk.NorthWallTemp, SolnBlk.NorthWallEmiss );
	    dx_norm = dX*SolnBlk.Grid.nfaceN(i,SolnBlk.JCu);
	    dU = Uwall-SolnBlk.U[i][SolnBlk.JCu];
	    dUdx = dU/dx_norm;
	    for( int ghost = 1; ghost <= SolnBlk.Nghost; ghost++){
	      dX = SolnBlk.Grid.Cell[i][SolnBlk.JCu+ghost].Xc -
		   SolnBlk.Grid.Cell[i][SolnBlk.JCu].Xc;
	      dx_norm = dX*SolnBlk.Grid.nfaceN(i, SolnBlk.JCu);
	      SolnBlk.U[i][SolnBlk.JCu+ghost] = SolnBlk.U[i][SolnBlk.JCu] + dUdx*dx_norm;
	    }
// 	    SolnBlk.UoN[i] = Gray_Wall(SolnBlk.U[i][SolnBlk.JCu], SolnBlk.Grid.nfaceN(i,SolnBlk.JCu), 
// 				       SolnBlk.NorthWallTemp, SolnBlk.NorthWallEmiss );
// 	    for( int ghost = 1; ghost <= SolnBlk.Nghost; ghost++){
// 	      SolnBlk.U[i][SolnBlk.JCu+ghost] = SolnBlk.UoN[i];
// 	    }
            break;

          default:
	    for( int ghost = 1; ghost <= SolnBlk.Nghost; ghost++){
	      SolnBlk.U[i][SolnBlk.JCu+ghost] = SolnBlk.U[i][SolnBlk.JCu];
	    }
            break;
        } /* endswitch */
    } /* endfor */


}


/********************************************************
 * Routine: BCs_Space_March                             *
 *                                                      *
 * Apply boundary conditions at boundaries of the       *
 * specified quadrilateral solution block.              *
 * Used for Space Marching                              *
 *                                                      *
 ********************************************************/
void BCs_Space_March(Rte2D_Quad_Block &SolnBlk,
		     Rte2D_Input_Parameters &IP) {

    int i, j;
    double dx_norm;
    Vector2D dX;
    Rte2D_State dU, dUdx;


    for ( j = SolnBlk.JCl-SolnBlk.Nghost ; j <= SolnBlk.JCu+SolnBlk.Nghost ; ++j ) {
      // Prescribe West boundary conditions.
      if ( (j >= SolnBlk.JCl && j <= SolnBlk.JCu) ||
           (j < SolnBlk.JCl && 
            (SolnBlk.Grid.BCtypeS[SolnBlk.ICl-1] == BC_NONE ||
	     SolnBlk.Grid.BCtypeS[SolnBlk.ICl-1] == BC_GRAY_WALL ) ) ||
           (j > SolnBlk.JCu && 
            (SolnBlk.Grid.BCtypeN[SolnBlk.ICl-1] == BC_NONE ||
	     SolnBlk.Grid.BCtypeN[SolnBlk.ICl-1] == BC_GRAY_WALL  ) ) ) {
        switch(SolnBlk.Grid.BCtypeW[j]) {
          case BC_NONE :
            break;
          case BC_REFLECTION :
	    Reflect_Space_March(SolnBlk.UoW[j], SolnBlk.Grid.nfaceW(SolnBlk.ICl,j));
	    for( int ghost = 1; ghost <= SolnBlk.Nghost; ghost++){
	      SolnBlk.U[SolnBlk.ICl-ghost][j] = Reflect(SolnBlk.U[SolnBlk.ICl+ghost-1][j],
							SolnBlk.Grid.nfaceW(SolnBlk.ICl,j));
	    }
           break;
          case BC_GRAY_WALL :
	    Gray_Wall_Space_March(SolnBlk.UoW[j], SolnBlk.M[0][0],
				  SolnBlk.Grid.nfaceW(SolnBlk.ICl, j), 
				  SolnBlk.WestWallTemp, SolnBlk.WestWallEmiss );
	    for( int ghost = 1; ghost <= SolnBlk.Nghost; ghost++){
	      SolnBlk.U[SolnBlk.ICl-ghost][j] = SolnBlk.UoW[j]; 
	    }
            break;
          default:
            break;
        } /* endswitch */
      } /* endif */


      // Prescribe East boundary conditions.
      if ( (j >= SolnBlk.JCl && j <= SolnBlk.JCu) ||
           (j < SolnBlk.JCl && 
            (SolnBlk.Grid.BCtypeS[SolnBlk.ICu+1] == BC_NONE ||
	     SolnBlk.Grid.BCtypeS[SolnBlk.ICu+1] == BC_GRAY_WALL ) ) ||
           (j > SolnBlk.JCu && 
            (SolnBlk.Grid.BCtypeN[SolnBlk.ICu+1] == BC_NONE ||
	     SolnBlk.Grid.BCtypeN[SolnBlk.ICu+1] == BC_GRAY_WALL ) ) ) {
        switch(SolnBlk.Grid.BCtypeE[j]) {
          case BC_NONE :
            break;
          case BC_REFLECTION :
	    Reflect_Space_March(SolnBlk.UoE[j], SolnBlk.Grid.nfaceE(SolnBlk.ICu,j));
	    for( int ghost = 1; ghost <= SolnBlk.Nghost; ghost++){
	      SolnBlk.U[SolnBlk.ICu+ghost][j] = Reflect(SolnBlk.U[SolnBlk.ICu-ghost+1][j],
							SolnBlk.Grid.nfaceE(SolnBlk.ICu,j));
	    }       
            break;
          case BC_GRAY_WALL :
	    Gray_Wall_Space_March(SolnBlk.UoE[j], SolnBlk.M[0][0],
				  SolnBlk.Grid.nfaceE(SolnBlk.ICu, j), 
				  SolnBlk.EastWallTemp, SolnBlk.EastWallEmiss );
	    for( int ghost = 1; ghost <= SolnBlk.Nghost; ghost++){
	      SolnBlk.U[SolnBlk.ICu+ghost][j] = SolnBlk.UoE[j];
	    }
            break;

	  default:
            break;
        } /* endswitch */
      } /* endif */

    } /* endfor */

    for ( i = SolnBlk.ICl-SolnBlk.Nghost ; i <= SolnBlk.ICu+SolnBlk.Nghost ; ++i ) {

        // Prescribe South boundary conditions.
        switch(SolnBlk.Grid.BCtypeS[i]) {
          case BC_NONE :
            break;
          case BC_REFLECTION :
	    Reflect_Space_March(SolnBlk.UoS[i], SolnBlk.Grid.nfaceS(i, SolnBlk.JCl));
	    for( int ghost = 1; ghost <= SolnBlk.Nghost; ghost++){
	      SolnBlk.U[i][SolnBlk.JCl-ghost] = Reflect(SolnBlk.U[i][SolnBlk.JCl + ghost-1],
							SolnBlk.Grid.nfaceS(i, SolnBlk.JCl));
	    }
          break;
          case BC_GRAY_WALL :
	    Gray_Wall_Space_March(SolnBlk.UoS[i], SolnBlk.M[0][0],
				  SolnBlk.Grid.nfaceS(i,SolnBlk.JCl), 
				  SolnBlk.SouthWallTemp, SolnBlk.SouthWallEmiss );
	    for( int ghost = 1; ghost <= SolnBlk.Nghost; ghost++){
	      SolnBlk.U[i][SolnBlk.JCl-ghost] = SolnBlk.UoS[i];
	    }
            break;

          default:
           break;
        } /* endswitch */

        // Prescribe North boundary conditions.
        switch(SolnBlk.Grid.BCtypeN[i]) {
          case BC_NONE :
            break;
          case BC_REFLECTION :
	    Reflect_Space_March(SolnBlk.UoN[i], SolnBlk.Grid.nfaceN(i, SolnBlk.JCu));

	    for( int ghost = 1; ghost <= SolnBlk.Nghost; ghost++){
	      SolnBlk.U[i][SolnBlk.JCu+ghost] = Reflect(SolnBlk.U[i][SolnBlk.JCu-ghost+1],
							SolnBlk.Grid.nfaceN(i, SolnBlk.JCu));
	    }
            break;
          case BC_GRAY_WALL :
	    Gray_Wall_Space_March(SolnBlk.UoN[i], SolnBlk.M[0][0],
				  SolnBlk.Grid.nfaceN(i,SolnBlk.JCu), 
				  SolnBlk.NorthWallTemp, SolnBlk.NorthWallEmiss );
	    for( int ghost = 1; ghost <= SolnBlk.Nghost; ghost++){
	      SolnBlk.U[i][SolnBlk.JCu+ghost] = SolnBlk.UoN[i];
	    }
            break;

          default:
            break;
        } /* endswitch */
    } /* endfor */


}



/********************************************************
 * Routine: CFL                                         *
 *                                                      *
 * Determines the allowable global and local time steps *
 * (for explicit Rte time stepping scheme) for the    *
 * specified quadrilateral solution block according to  *
 * the Courant-Friedrichs-Lewy condition.               *
 *                                                      *
 ********************************************************/
double CFL(Rte2D_Quad_Block &SolnBlk,
           Rte2D_Input_Parameters &Input_Parameters) {

    int i, j;
    double dtMin, d_i, d_j, v_i, v_j, a;
    double dt_cfl, dt_src, temp, beta;

    dtMin = MILLION;
    dt_cfl = MILLION;
    dt_src = MILLION;
    
    int JCl_overlap = 0;
    int JCu_overlap = 0;
    int ICu_overlap = 0;
    int ICl_overlap = 0;
    
    // Modifications for NKS overlap 
    if(Input_Parameters.NKS_IP.GMRES_Overlap > 0 ){	
      if (SolnBlk.Grid.BCtypeS[SolnBlk.ICl] == BC_NONE)  JCl_overlap = Input_Parameters.NKS_IP.GMRES_Overlap;
      if (SolnBlk.Grid.BCtypeN[SolnBlk.ICu] == BC_NONE)  JCu_overlap = Input_Parameters.NKS_IP.GMRES_Overlap;
      if (SolnBlk.Grid.BCtypeE[SolnBlk.JCu] == BC_NONE)  ICu_overlap = Input_Parameters.NKS_IP.GMRES_Overlap;
      if (SolnBlk.Grid.BCtypeW[SolnBlk.JCl] == BC_NONE)  ICl_overlap = Input_Parameters.NKS_IP.GMRES_Overlap;       
    }

    for ( j  = SolnBlk.JCl-SolnBlk.Nghost ; j <= SolnBlk.JCu+SolnBlk.Nghost ; ++j ) {
       for ( i = SolnBlk.ICl-SolnBlk.Nghost ; i <= SolnBlk.ICu+SolnBlk.Nghost ; ++i ) {
	 if (i < SolnBlk.ICl-ICl_overlap  || i > SolnBlk.ICu+ICu_overlap ||
	     j < SolnBlk.JCl-JCl_overlap  || j > SolnBlk.JCu+JCu_overlap) {
	   SolnBlk.dt[i][j] = ZERO;
	 } else {

	   /************ Scalar Advection Terms ************************/

	   d_i = TWO*(SolnBlk.Grid.Cell[i][j].A*SolnBlk.Sp[i][j]/
		      (SolnBlk.Grid.lfaceE(i, j)*SolnBlk.SpE[i][j]
		       +SolnBlk.Grid.lfaceW(i, j)*SolnBlk.SpW[i][j]));
	   d_j = TWO*(SolnBlk.Grid.Cell[i][j].A*SolnBlk.Sp[i][j]/
		      (SolnBlk.Grid.lfaceN(i, j)*SolnBlk.SpN[i][j]
		       +SolnBlk.Grid.lfaceS(i, j)*SolnBlk.SpS[i][j]));

	   //------------------------------------------------
	   // determine the minimum allowable timestep
	   //------------------------------------------------

	   // 
	   // loop over dirs
	   //
	   for (int m=0; m<SolnBlk.U[i][j].Npolar; m++) 
	     for (int l=0; l<SolnBlk.U[i][j].Nazim[m]; l++) {
	   
	       // convective contribution
	       v_i = HALF*(SolnBlk.U[i][j].mu[m][l]* 
			   (SolnBlk.Grid.nfaceE(i, j).x-SolnBlk.Grid.nfaceW(i, j).x));
	       v_i += HALF*(SolnBlk.U[i][j].eta[m][l]* 
			    (SolnBlk.Grid.nfaceE(i, j).y-SolnBlk.Grid.nfaceW(i, j).y));
	       v_j = HALF*(SolnBlk.U[i][j].mu[m][l]*
			   (SolnBlk.Grid.nfaceN(i, j).x-SolnBlk.Grid.nfaceS(i, j).x));
	       v_j += HALF*(SolnBlk.U[i][j].eta[m][l]*
			    (SolnBlk.Grid.nfaceN(i, j).y-SolnBlk.Grid.nfaceS(i, j).y));
	   
	       if (fabs(v_i) > TOLER) {
		 temp = d_i/fabs(v_i);
	       } else {
		 temp = MILLION;
	       } /* endif */
	       if (fabs(v_j) > TOLER) {
		 temp = min(temp, d_j/fabs(v_j));
	       } /* endif */
	       
	       dt_cfl = min( temp, dt_cfl );
	       
	       // source contribution
	       for (int v=0; v<SolnBlk.U[i][j].Nband; v++ ) {
		 beta = SolnBlk.M[i][j].beta(v);
		 temp = HALF / (beta*SolnBlk.U[i][j].omega[m][l]);
		 dt_src = min( temp, dt_src );
	       } // endfor - bands

	     } // endfor - dirs
	   
	   // determine the limiting case
	   SolnBlk.dt[i][j] = min( dt_cfl, dt_src );
  
	   
	   /************ Global Minimum ********************************/
	   dtMin = min(dtMin, SolnBlk.dt[i][j]);
	   
	 } /* endif */
       } /* endfor */
    } /* endfor */

    for ( j  = SolnBlk.JCl-SolnBlk.Nghost ; j <= SolnBlk.JCu+SolnBlk.Nghost ; ++j ) {
       for ( i = SolnBlk.ICl-SolnBlk.Nghost ; i <= SolnBlk.ICu+SolnBlk.Nghost ; ++i ) {
 	  if (i < SolnBlk.ICl-ICl_overlap  || i > SolnBlk.ICu+ICu_overlap ||
              j < SolnBlk.JCl-JCl_overlap  || j > SolnBlk.JCu+JCu_overlap) {
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
void Set_Global_TimeStep(Rte2D_Quad_Block &SolnBlk,
                         const double &Dt_min) {

  for ( int j  = SolnBlk.JCl-SolnBlk.Nghost ; j <= SolnBlk.JCu+SolnBlk.Nghost ; ++j ) {
    for ( int i = SolnBlk.ICl-SolnBlk.Nghost ; i <= SolnBlk.ICu+SolnBlk.Nghost ; ++i ) {
      SolnBlk.dt[i][j] = Dt_min;
    } /* endfor */
  } /* endfor */
   
  //cout<<"\n Global Time Step used "; 

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
double L1_Norm_Residual(Rte2D_Quad_Block &SolnBlk, const int &norm) {

    int i, j, n;
    double l1_norm;
    int NUM_VAR_RTE2D = SolnBlk.NumVar();

    l1_norm = ZERO;

    for ( j  = SolnBlk.JCl ; j <= SolnBlk.JCu ; ++j ) {
       for ( i = SolnBlk.ICl ; i <= SolnBlk.ICu ; ++i ) {
	 for ( n = 0 ; n < NUM_VAR_RTE2D ; ++n ) 
	   l1_norm += fabs(SolnBlk.dUdt[i][j][0].I[n]);
       } /* endfor */
    } /* endfor */

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
double L2_Norm_Residual(Rte2D_Quad_Block &SolnBlk, const int &norm) {

    int i, j, n;
    double l2_norm;
    int NUM_VAR_RTE2D = SolnBlk.NumVar();

    l2_norm = ZERO;

    for ( j  = SolnBlk.JCl ; j <= SolnBlk.JCu ; ++j ) {
       for ( i = SolnBlk.ICl ; i <= SolnBlk.ICu ; ++i ) {
	 for ( n = 0 ; n < NUM_VAR_RTE2D ; ++n ) 
	   l2_norm += sqr(SolnBlk.dUdt[i][j][0].I[n]);
       } /* endfor */
    } /* endfor */

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
double Max_Norm_Residual(Rte2D_Quad_Block &SolnBlk, const int &norm) {

    int i, j, n;
    double max_norm;
    double temp(ZERO);
    int NUM_VAR_RTE2D = SolnBlk.NumVar();

    max_norm = ZERO;

    for ( j  = SolnBlk.JCl ; j <= SolnBlk.JCu ; ++j ) {
       for ( i = SolnBlk.ICl ; i <= SolnBlk.ICu ; ++i ) {
	 for ( n = 0 ; n < NUM_VAR_RTE2D ; ++n ) {
	   temp += sqr(SolnBlk.dUdt[i][j][0].I[n]);
	 }
	 max_norm = max(max_norm, sqrt(temp));
       } /* endfor */
    } /* endfor */

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
void Linear_Reconstruction_GreenGauss(Rte2D_Quad_Block &SolnBlk,
				      const int i, 
                                      const int j,
                                      const int Limiter) {

    int NUM_VAR_RTE2D = SolnBlk.NumVar(); 
    int n, n2, n_pts, i_index[8], j_index[8];
    double u0Min, u0Max, uQuad[4], phi;
    double DxDx_ave, DxDy_ave, DyDy_ave;
    double l_north, l_south, l_east, l_west;
    Vector2D n_north, n_south, n_east, n_west, dX;
    Rte2D_State U_nw, U_ne, U_sw, U_se, U_face, 
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
		 SolnBlk.Grid.BCtypeW[j] == BC_GRAY_WALL ) {
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
		 SolnBlk.Grid.BCtypeE[j] == BC_GRAY_WALL ) {
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
                 SolnBlk.Grid.BCtypeS[i] == BC_GRAY_WALL) {
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
                 SolnBlk.Grid.BCtypeN[i] == BC_GRAY_WALL ) {
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
	   U_nw = SolnBlk.UnNW(i, j);
           U_ne = SolnBlk.UnNE(i, j);
           U_sw = SolnBlk.UnSW(i, j);
           U_se = SolnBlk.UnSE(i, j);

           l_north = SolnBlk.Grid.lfaceN(i, j);
           l_south = SolnBlk.Grid.lfaceS(i, j);
           l_east = SolnBlk.Grid.lfaceE(i, j);
           l_west = SolnBlk.Grid.lfaceW(i, j);

           n_north = SolnBlk.Grid.nfaceN(i, j);
           n_south = SolnBlk.Grid.nfaceS(i, j);
           n_east = SolnBlk.Grid.nfaceE(i, j);
           n_west = SolnBlk.Grid.nfaceW(i, j);

           U_face = HALF*(U_nw+U_ne)*l_north; 
           SolnBlk.dUdx[i][j] = U_face*n_north.x;
	   SolnBlk.dUdy[i][j] = U_face*n_north.y;

           U_face = HALF*(U_sw+U_se)*l_south; 
           SolnBlk.dUdx[i][j] += U_face*n_south.x;
	   SolnBlk.dUdy[i][j] += U_face*n_south.y;

           U_face = HALF*(U_ne+U_se)*l_east; 
           SolnBlk.dUdx[i][j] += U_face*n_east.x;
	   SolnBlk.dUdy[i][j] += U_face*n_east.y;

           U_face = HALF*(U_nw+U_sw)*l_west; 
           SolnBlk.dUdx[i][j] += U_face*n_west.x;
	   SolnBlk.dUdy[i][j] += U_face*n_west.y;

           SolnBlk.dUdx[i][j] = SolnBlk.dUdx[i][j]/
                                SolnBlk.Grid.Cell[i][j].A;
           SolnBlk.dUdy[i][j] = SolnBlk.dUdy[i][j]/
                                SolnBlk.Grid.Cell[i][j].A;

        // If <8 neighbours are used, apply least-squares reconstruction
        } else {
	   DUDx_ave.Zero();
           DUDy_ave.Zero();
           DxDx_ave = ZERO;
           DxDy_ave = ZERO;
           DyDy_ave = ZERO;
    
           for ( n2 = 0 ; n2 <= n_pts-1 ; ++n2 ) {
               dX = SolnBlk.Grid.Cell[ i_index[n2] ][ j_index[n2] ].Xc - 
                    SolnBlk.Grid.Cell[i][j].Xc;
               DU = SolnBlk.U[ i_index[n2] ][ j_index[n2] ] - 
                    SolnBlk.U[i][j];
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
           SolnBlk.dUdx[i][j] = (DUDx_ave*DyDy_ave-DUDy_ave*DxDy_ave)/
                                (DxDx_ave*DyDy_ave-DxDy_ave*DxDy_ave);
           SolnBlk.dUdy[i][j] = (DUDy_ave*DxDx_ave-DUDx_ave*DxDy_ave)/
                                (DxDx_ave*DyDy_ave-DxDy_ave*DxDy_ave);
        } /* endif */

        // Calculate slope limiters.    
	if (!SolnBlk.Freeze_Limiter) {
           for ( n = 1 ; n <= NUM_VAR_RTE2D ; ++n ) {
              u0Min = SolnBlk.U[i][j][n];
              u0Max = u0Min;
              for ( n2 = 0 ; n2 <= n_pts-1 ; ++n2 ) {
                 u0Min = min(u0Min, SolnBlk.U[ i_index[n2] ][ j_index[n2] ][n]);
                 u0Max = max(u0Max, SolnBlk.U[ i_index[n2] ][ j_index[n2] ][n]);
              } /* endfor */
    
              dX = SolnBlk.Grid.xfaceE(i, j)-SolnBlk.Grid.Cell[i][j].Xc;
              uQuad[0] = SolnBlk.U[i][j][n] + 
                         SolnBlk.dUdx[i][j][n]*dX.x +
                         SolnBlk.dUdy[i][j][n]*dX.y ;
              dX = SolnBlk.Grid.xfaceW(i, j)-SolnBlk.Grid.Cell[i][j].Xc;
              uQuad[1] = SolnBlk.U[i][j][n] + 
                         SolnBlk.dUdx[i][j][n]*dX.x +
                         SolnBlk.dUdy[i][j][n]*dX.y ;
              dX = SolnBlk.Grid.xfaceN(i, j)-SolnBlk.Grid.Cell[i][j].Xc;
              uQuad[2] = SolnBlk.U[i][j][n] + 
                         SolnBlk.dUdx[i][j][n]*dX.x +
                         SolnBlk.dUdy[i][j][n]*dX.y ;
              dX = SolnBlk.Grid.xfaceS(i, j)-SolnBlk.Grid.Cell[i][j].Xc;
              uQuad[3] = SolnBlk.U[i][j][n] + 
                         SolnBlk.dUdx[i][j][n]*dX.x +
                         SolnBlk.dUdy[i][j][n]*dX.y ;
    
              switch(Limiter) {
                case LIMITER_ONE :
                  phi = ONE;
                  break;
                case LIMITER_ZERO :
                  phi = ZERO;
                  break;
                case LIMITER_BARTH_JESPERSEN :
                  phi = Limiter_BarthJespersen(uQuad, SolnBlk.U[i][j][n], 
                                               u0Min, u0Max, 4);
                  break;
                case LIMITER_VENKATAKRISHNAN :
                  phi = Limiter_Venkatakrishnan(uQuad, SolnBlk.U[i][j][n], 
                                                u0Min, u0Max, 4);
                  break;
                case LIMITER_VANLEER :
                  phi = Limiter_VanLeer(uQuad, SolnBlk.U[i][j][n], 
                                        u0Min, u0Max, 4);
                  break;
                case LIMITER_VANALBADA :
                  phi = Limiter_VanAlbada(uQuad, SolnBlk.U[i][j][n], 
                                          u0Min, u0Max, 4);
                  break;
                default:
                  phi = Limiter_BarthJespersen(uQuad, SolnBlk.U[i][j][n], 
                                               u0Min, u0Max, 4);
                  break;
              } /* endswitch */

	      SolnBlk.phi[i][j][n] = phi;

           } /* endfor */
        } /* endif */
    } else {
      SolnBlk.dUdx[i][j].Zero();
        SolnBlk.dUdy[i][j].Zero(); 
        SolnBlk.phi[i][j].Zero();
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
void Linear_Reconstruction_GreenGauss(Rte2D_Quad_Block &SolnBlk,
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
void Linear_Reconstruction_LeastSquares(Rte2D_Quad_Block &SolnBlk,
				        const int i, 
                                        const int j,
                                        const int Limiter) {

    int NUM_VAR_RTE2D = SolnBlk.NumVar(); 
    int n, n2, n_pts, i_index[8], j_index[8];
    double u0Min, u0Max, uQuad[4], phi;
    double DxDx_ave, DxDy_ave, DyDy_ave;
    Vector2D dX;
    Rte2D_State DU, DUDx_ave, DUDy_ave;

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
                 SolnBlk.Grid.BCtypeW[j] == BC_LINEAR_EXTRAPOLATION  ||
                 SolnBlk.Grid.BCtypeW[j] == BC_GRAY_WALL) {
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
                 SolnBlk.Grid.BCtypeE[j] == BC_GRAY_WALL ) {
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
                 SolnBlk.Grid.BCtypeS[i] == BC_LINEAR_EXTRAPOLATION  ||
                 SolnBlk.Grid.BCtypeS[i] == BC_GRAY_WALL ) {
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
                 SolnBlk.Grid.BCtypeN[i] == BC_GRAY_WALL ) {
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
    
    if (n_pts > 0) {
        DUDx_ave.Zero();
        DUDy_ave.Zero();
        DxDx_ave = ZERO;
        DxDy_ave = ZERO;
        DyDy_ave = ZERO;
    
        for ( n2 = 0 ; n2 <= n_pts-1 ; ++n2 ) {
            dX = SolnBlk.Grid.Cell[ i_index[n2] ][ j_index[n2] ].Xc - 
                 SolnBlk.Grid.Cell[i][j].Xc;
            DU = SolnBlk.U[ i_index[n2] ][ j_index[n2] ] - 
                 SolnBlk.U[i][j];
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
        SolnBlk.dUdx[i][j] = (DUDx_ave*DyDy_ave-DUDy_ave*DxDy_ave)/
                             (DxDx_ave*DyDy_ave-DxDy_ave*DxDy_ave);
        SolnBlk.dUdy[i][j] = (DUDy_ave*DxDx_ave-DUDx_ave*DxDy_ave)/
                             (DxDx_ave*DyDy_ave-DxDy_ave*DxDy_ave);
    
        // Calculate slope limiters. 
	if (!SolnBlk.Freeze_Limiter) {
           for ( n = 1 ; n <= NUM_VAR_RTE2D ; ++n ) {
              u0Min = SolnBlk.U[i][j][n];
              u0Max = u0Min;
              for ( n2 = 0 ; n2 <= n_pts-1 ; ++n2 ) {
                 u0Min = min(u0Min, SolnBlk.U[ i_index[n2] ][ j_index[n2] ][n]);
                 u0Max = max(u0Max, SolnBlk.U[ i_index[n2] ][ j_index[n2] ][n]);
              } /* endfor */
    
              dX = SolnBlk.Grid.xfaceE(i, j)-SolnBlk.Grid.Cell[i][j].Xc;
              uQuad[0] = SolnBlk.U[i][j][n] + 
                         SolnBlk.dUdx[i][j][n]*dX.x +
                         SolnBlk.dUdy[i][j][n]*dX.y ;
              dX = SolnBlk.Grid.xfaceW(i, j)-SolnBlk.Grid.Cell[i][j].Xc;
              uQuad[1] = SolnBlk.U[i][j][n] + 
                         SolnBlk.dUdx[i][j][n]*dX.x +
                         SolnBlk.dUdy[i][j][n]*dX.y ;
              dX = SolnBlk.Grid.xfaceN(i, j)-SolnBlk.Grid.Cell[i][j].Xc;
              uQuad[2] = SolnBlk.U[i][j][n] + 
                         SolnBlk.dUdx[i][j][n]*dX.x +
                         SolnBlk.dUdy[i][j][n]*dX.y ;
              dX = SolnBlk.Grid.xfaceS(i, j)-SolnBlk.Grid.Cell[i][j].Xc;
              uQuad[3] = SolnBlk.U[i][j][n] + 
                         SolnBlk.dUdx[i][j][n]*dX.x +
                         SolnBlk.dUdy[i][j][n]*dX.y ;
    
              switch(Limiter) {
                case LIMITER_ONE :
                  phi = ONE;
                  break;
                case LIMITER_ZERO :
                  phi = ZERO;
                  break;
                case LIMITER_BARTH_JESPERSEN :
                  phi = Limiter_BarthJespersen(uQuad, SolnBlk.U[i][j][n], 
                                               u0Min, u0Max, 4);
                  break;
                case LIMITER_VENKATAKRISHNAN :
                  phi = Limiter_Venkatakrishnan(uQuad, SolnBlk.U[i][j][n], 
                                                u0Min, u0Max, 4);
                  break;
                case LIMITER_VANLEER :
                  phi = Limiter_VanLeer(uQuad, SolnBlk.U[i][j][n], 
                                        u0Min, u0Max, 4);
                  break;
                case LIMITER_VANALBADA :
                  phi = Limiter_VanAlbada(uQuad, SolnBlk.U[i][j][n], 
                                          u0Min, u0Max, 4);
                  break;
                default:
                  phi = Limiter_BarthJespersen(uQuad, SolnBlk.U[i][j][n], 
                                               u0Min, u0Max, 4);
                  break;
              } /* endswitch */

	      SolnBlk.phi[i][j][n] = phi;

           } /* endfor */
        } /* endif */
    } else {
        SolnBlk.dUdx[i][j].Zero();
        SolnBlk.dUdy[i][j].Zero(); 
        SolnBlk.phi[i][j].Zero();
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
void Linear_Reconstruction_LeastSquares_2(Rte2D_Quad_Block &SolnBlk,
				          const int i, 
                                          const int j,
                                          const int Limiter) {

    int NUM_VAR_RTE2D = SolnBlk.NumVar(); 
    int n, n2, n_pts, i_index[8], j_index[8];
    double u0Min, u0Max, uQuad[4], phi;
    double DxDx_ave, DxDy_ave, DyDy_ave;
    Vector2D dX;
    Rte2D_State DU, DUDx_ave, DUDy_ave;

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
                 SolnBlk.Grid.BCtypeW[j] == BC_GRAY_WALL ) {
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
                 SolnBlk.Grid.BCtypeE[j] == BC_GRAY_WALL ) {
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
                 SolnBlk.Grid.BCtypeS[i] == BC_GRAY_WALL ) {
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
                 SolnBlk.Grid.BCtypeN[i] == BC_GRAY_WALL ) {
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
        DUDx_ave.Zero();
        DUDy_ave.Zero();
        DxDx_ave = ZERO;
        DxDy_ave = ZERO;
        DyDy_ave = ZERO;
    
        for ( n2 = 0 ; n2 <= n_pts-1 ; ++n2 ) {
            dX = SolnBlk.Grid.Cell[ i_index[n2] ][ j_index[n2] ].Xc - 
                 SolnBlk.Grid.Cell[i][j].Xc;
            DU = SolnBlk.U[ i_index[n2] ][ j_index[n2] ] - 
                 SolnBlk.U[i][j];
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
        SolnBlk.dUdx[i][j] = (DUDx_ave*DyDy_ave-DUDy_ave*DxDy_ave)/
                             (DxDx_ave*DyDy_ave-DxDy_ave*DxDy_ave);
        SolnBlk.dUdy[i][j] = (DUDy_ave*DxDx_ave-DUDx_ave*DxDy_ave)/
                             (DxDx_ave*DyDy_ave-DxDy_ave*DxDy_ave);
    
        // Calculate slope limiters. 
	if (!SolnBlk.Freeze_Limiter) {
           for ( n = 1 ; n <= NUM_VAR_RTE2D ; ++n ) {
              u0Min = SolnBlk.U[i][j][n];
              u0Max = u0Min;
              for ( n2 = 0 ; n2 <= n_pts-1 ; ++n2 ) {
                 u0Min = min(u0Min, SolnBlk.U[ i_index[n2] ][ j_index[n2] ][n]);
                 u0Max = max(u0Max, SolnBlk.U[ i_index[n2] ][ j_index[n2] ][n]);
              } /* endfor */
    
              dX = SolnBlk.Grid.xfaceE(i, j)-SolnBlk.Grid.Cell[i][j].Xc;
              uQuad[0] = SolnBlk.U[i][j][n] + 
                         SolnBlk.dUdx[i][j][n]*dX.x +
                         SolnBlk.dUdy[i][j][n]*dX.y ;
              dX = SolnBlk.Grid.xfaceW(i, j)-SolnBlk.Grid.Cell[i][j].Xc;
              uQuad[1] = SolnBlk.U[i][j][n] + 
                         SolnBlk.dUdx[i][j][n]*dX.x +
                         SolnBlk.dUdy[i][j][n]*dX.y ;
              dX = SolnBlk.Grid.xfaceN(i, j)-SolnBlk.Grid.Cell[i][j].Xc;
              uQuad[2] = SolnBlk.U[i][j][n] + 
                         SolnBlk.dUdx[i][j][n]*dX.x +
                         SolnBlk.dUdy[i][j][n]*dX.y ;
              dX = SolnBlk.Grid.xfaceS(i, j)-SolnBlk.Grid.Cell[i][j].Xc;
              uQuad[3] = SolnBlk.U[i][j][n] + 
                         SolnBlk.dUdx[i][j][n]*dX.x +
                         SolnBlk.dUdy[i][j][n]*dX.y ;
    
              switch(Limiter) {
                case LIMITER_ONE :
                  phi = ONE;
                  break;
                case LIMITER_ZERO :
                  phi = ZERO;
                  break;
                case LIMITER_BARTH_JESPERSEN :
                  phi = Limiter_BarthJespersen(uQuad, SolnBlk.U[i][j][n], 
                                               u0Min, u0Max, 4);
                  break;
                case LIMITER_VENKATAKRISHNAN :
                  phi = Limiter_Venkatakrishnan(uQuad, SolnBlk.U[i][j][n], 
                                                u0Min, u0Max, 4);
                  break;
                case LIMITER_VANLEER :
                  phi = Limiter_VanLeer(uQuad, SolnBlk.U[i][j][n], 
                                        u0Min, u0Max, 4);
                  break;
                case LIMITER_VANALBADA :
                  phi = Limiter_VanAlbada(uQuad, SolnBlk.U[i][j][n], 
                                          u0Min, u0Max, 4);
                  break;
                default:
                  phi = Limiter_BarthJespersen(uQuad, SolnBlk.U[i][j][n], 
                                               u0Min, u0Max, 4);
                  break;
              } /* endswitch */

	      SolnBlk.phi[i][j][n] = phi;

           } /* endfor */
        } /* endif */
    } else {
        SolnBlk.dUdx[i][j].Zero();
        SolnBlk.dUdy[i][j].Zero(); 
        SolnBlk.phi[i][j].Zero();
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
void Linear_Reconstruction_LeastSquares(Rte2D_Quad_Block &SolnBlk,
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


/********************************************************
 * Routine: Linear_Reconstruction_LeastSquares_Angular  *
 *                                                      *
 * Performs the reconstruction of a limited piecewise   *
 * linear solution state in the azimuthal direction     *
 * within each cell of the                              *
 * computational mesh of the specified quadrilateral    *
 * solution block.  A least squares approach is         *
 * used in the evaluation of the unlimited solution     *
 * gradients.  Several slope limiters may be used.      *
 *                                                      *
 ********************************************************/
void Linear_Reconstruction_LeastSquares_Angular(Rte2D_Quad_Block &SolnBlk,
						const int Limiter) {

    int i, j;

    /* Carry out the limited solution reconstruction in
       each cell of the computational mesh. */

    for ( j  = SolnBlk.JCl-SolnBlk.Nghost+1 ; j <= SolnBlk.JCu+SolnBlk.Nghost-1 ; ++j ) {
       for ( i = SolnBlk.ICl-SolnBlk.Nghost+1 ; i <= SolnBlk.ICu+SolnBlk.Nghost-1 ; ++i ) {
	   Linear_Reconstruction_LeastSquares_Angular(SolnBlk, i, j, Limiter);
       } /* endfor */
    } /* endfor */

}

/********************************************************
 * Routine: Linear_Reconstruction_LeastSquares_Angular  *
 *                                                      *
 * Performs the reconstruction of a limited piecewise   *
 * linear solution state in the azimuthal direction     *
 * within each cell of the                              *
 * computational mesh of the specified quadrilateral    *
 * solution block.  A least squares approach is         *
 * used in the evaluation of the unlimited solution     *
 * gradients.  Several slope limiters may be used.      *
 *                                                      *
 ********************************************************/
void Linear_Reconstruction_LeastSquares_Angular(Rte2D_Quad_Block &SolnBlk,
						const int i, 
						const int j,
						const int Limiter) {
  // constants
  const int n_pts = 2;

  // declares
  Rte2D_State &U = SolnBlk.U[0][0]; // temp pointer
  int index[n_pts];
  double u0Min, u0Max, uQuad[n_pts], phi;
  double Dx, DxDx_ave;
  double DU, DUDx_ave;

  //------------------------------------------------
  // Carry out the limited solution reconstruction in
  // each control angle. 
  //------------------------------------------------

  //
  // loop over dirs
  //
  for (int m=0; m<U.Npolar; m++) 
    for (int l=0; l<U.Nazim[m]; l++) {

      // set the neighboring nodes
      if (l==0) {
	index[0] = 0; // symmetry
	index[1] = l+1; 
      } else if  (l==U.Nazim[m]-1) {
	index[0] = l-1;
	index[1] = U.Nazim[m]-1; // symmetry
      } else {
	index[0] = l-1;
	index[1] = l+1; 
      } /* endif */
	  
      //
      // loop over bands
      //
      for (int v=0; v<U.Nband; v++) {

	// zero
	DUDx_ave = ZERO;
	DxDx_ave = ZERO;
	for ( int n = 0 ; n < n_pts ; n++ ) {
	  Dx = U.psi[m][ index[n] ] - U.psi[m][l];
	  DU = SolnBlk.U[i][j].In(v,m,index[n]) - 
	    SolnBlk.U[i][j].In(v,m,l);
	  DUDx_ave += DU*Dx;
	  DxDx_ave += Dx*Dx;
	} /* endfor */
    					    
	DUDx_ave = DUDx_ave/double(n_pts);
	DxDx_ave = DxDx_ave/double(n_pts);

	// compute unlimited gradient
	SolnBlk.dUdpsi[i][j].In(v,m,l) = DUDx_ave/DxDx_ave;

	if (!SolnBlk.Freeze_Limiter) {

	  // determine max and min values of neighbors
	  u0Min = SolnBlk.U[i][j].In(v,m,l);
	  u0Max = u0Min;
	  for ( int n = 0; n<n_pts ; n++ ) {
	    u0Min = min(u0Min, SolnBlk.U[i][j].In(v,m,index[n]));
	    u0Max = max(u0Max, SolnBlk.U[i][j].In(v,m,index[n]));
	  } /* endfor */
	    
	    // compute unlimited face values
	  uQuad[0] = SolnBlk.U[i][j].In(v,m,l)  - HALF*SolnBlk.dUdpsi[i][j].In(v,m,l)*U.delta_psi[m][l];
	  uQuad[1] = SolnBlk.U[i][j].In(v,m,l)  + HALF*SolnBlk.dUdpsi[i][j].In(v,m,l)*U.delta_psi[m][l];
	    
	    
	  // determine the limiter
	  switch(Limiter) {
	  case LIMITER_ONE :
	    phi = ONE;
	    break;
	  case LIMITER_ZERO :
	    phi = ZERO;
	    break;
	  case LIMITER_BARTH_JESPERSEN :
	    phi = Limiter_BarthJespersen(uQuad, SolnBlk.U[i][j].In(v,m,l), 
					 u0Min, u0Max, n_pts);
	    break;
	  case LIMITER_VENKATAKRISHNAN :
	    phi = Limiter_Venkatakrishnan(uQuad, SolnBlk.U[i][j].In(v,m,l), 
					  u0Min, u0Max, n_pts);
	    break;
	  case LIMITER_VANLEER :
	    phi = Limiter_VanLeer(uQuad, SolnBlk.U[i][j].In(v,m,l), 
				  u0Min, u0Max, n_pts);
	    break;
	  case LIMITER_VANALBADA :
	    phi = Limiter_VanAlbada(uQuad, SolnBlk.U[i][j].In(v,m,l), 
				    u0Min, u0Max, n_pts);
	    break;
	  default:
	    phi = Limiter_BarthJespersen(uQuad, SolnBlk.U[i][j].In(v,m,l), 
					 u0Min, u0Max, n_pts);
	    break;
	  } // endswitch

	  // set the value
	  SolnBlk.phi_psi[i][j].In(v,m,l) = phi;

	} // endif - freeze limiter
       
      } // endfor - bands

    } // endfor - dirs

  return;

}


/********************************************************
 * Routine: Linear_Reconstruction_GreenGauss_Angular    *
 *                                                      *
 * Performs the reconstruction of a limited piecewise   *
 * linear solution state in the azimuthal direction     *
 * within each cell of the                              *
 * computational mesh for the specified quadrilateral   *
 * solution block.  A Green-Gauss approach is used      *
 * in the evaluation of the unlimited solution          *
 * gradients.  Several slope limiters may be used.      *
 *                                                      *
 ********************************************************/
void Linear_Reconstruction_GreenGauss_Angular(Rte2D_Quad_Block &SolnBlk,
					      const int Limiter) {

    int i, j;

    /* Carry out the limited solution reconstruction in
       each cell of the computational mesh. */

    for ( j  = SolnBlk.JCl-SolnBlk.Nghost+1 ; j <= SolnBlk.JCu+SolnBlk.Nghost-1 ; ++j ) {
       for ( i = SolnBlk.ICl-SolnBlk.Nghost+1 ; i <= SolnBlk.ICu+SolnBlk.Nghost-1 ; ++i ) {
	   Linear_Reconstruction_GreenGauss_Angular(SolnBlk, i, j, Limiter);
       } /* endfor */
    } /* endfor */

}

/********************************************************
 * Routine: Linear_Reconstruction_GreenGauss_Angular    *
 *                                                      *
 * Performs the reconstruction of a limited piecewise   *
 * linear solution state in the azimuthal direction     *
 * within each cell of the                              *
 * computational mesh for the specified quadrilateral   *
 * solution block.  A Green-Gauss approach is used      *
 * in the evaluation of the unlimited solution          *
 * gradients.  Several slope limiters may be used.      *
 *                                                      *
 ********************************************************/
void Linear_Reconstruction_GreenGauss_Angular(Rte2D_Quad_Block &SolnBlk,
					     const int i,
					     const int j,
					     const int Limiter) {

  // constants
  const int n_pts = 2;

  // declares
  Rte2D_State &U = SolnBlk.U[0][0];  // temp pointer
  int index[n_pts];
  double u0Min, u0Max, uQuad[n_pts], phi;

    
  //------------------------------------------------
  // Carry out the limited solution reconstruction in
  // each control angle. 
  //------------------------------------------------

  //
  // loop over directions
  //
  for (int m=0; m<U.Npolar; m++) 
    for (int l=0; l<U.Nazim[m]; l++) {

      // set the neighboring nodes
      if (l==0) {
	index[0] = 0; // symmetry
	index[1] = l+1; 
      } else if  (l==U.Nazim[m]-1) {
	index[0] = l-1;
	index[1] = U.Nazim[m]-1; // symmetry
      } else {
	index[0] = l-1;
	index[1] = l+1; 
      } /* endif */
      
      //
      // loop over bands
      //
      for (int v=0; v<U.Nband; v++) {
					    
	// compute unlimited gradient
	SolnBlk.dUdpsi[i][j].In(v,m,l) = HALF * 
	  ( SolnBlk.U[i][j].In(v,m,index[1]) - 
	    SolnBlk.U[i][j].In(v,m,index[0]) ) /
	  U.delta_psi[m][ l ];

	if (!SolnBlk.Freeze_Limiter) {

	  // determine max and min values of neighbors
	  u0Min = SolnBlk.U[i][j].In(v,m,l);
	  u0Max = u0Min;
	  for ( int n = 0; n<n_pts ; n++ ) {
	    u0Min = min(u0Min, SolnBlk.U[i][j].In(v,m,index[n]));
	    u0Max = max(u0Max, SolnBlk.U[i][j].In(v,m,index[n]));
	  } /* endfor */
	    
	    // compute unlimited face values
	  uQuad[0] = SolnBlk.U[i][j].In(v,m,l)  - HALF*SolnBlk.dUdpsi[i][j].In(v,m,l)*U.delta_psi[m][l];
	  uQuad[1] = SolnBlk.U[i][j].In(v,m,l)  + HALF*SolnBlk.dUdpsi[i][j].In(v,m,l)*U.delta_psi[m][l];
	    
	    
	  // determine the limiter
	  switch(Limiter) {
	  case LIMITER_ONE :
	    phi = ONE;
	    break;
	  case LIMITER_ZERO :
	    phi = ZERO;
	    break;
	  case LIMITER_BARTH_JESPERSEN :
	    phi = Limiter_BarthJespersen(uQuad, SolnBlk.U[i][j].In(v,m,l), 
					 u0Min, u0Max, n_pts);
	    break;
	  case LIMITER_VENKATAKRISHNAN :
	    phi = Limiter_Venkatakrishnan(uQuad, SolnBlk.U[i][j].In(v,m,l), 
					  u0Min, u0Max, n_pts);
	    break;
	  case LIMITER_VANLEER :
	    phi = Limiter_VanLeer(uQuad, SolnBlk.U[i][j].In(v,m,l), 
				  u0Min, u0Max, n_pts);
	    break;
	  case LIMITER_VANALBADA :
	    phi = Limiter_VanAlbada(uQuad, SolnBlk.U[i][j].In(v,m,l), 
				    u0Min, u0Max, n_pts);
	    break;
	  default:
	    phi = Limiter_BarthJespersen(uQuad, SolnBlk.U[i][j].In(v,m,l), 
					 u0Min, u0Max, n_pts);
	    break;
	  } // endswitch 

	  SolnBlk.phi_psi[i][j].In(v,m,l) = phi;

	} // endif - freeze limiter

      } // endfor - bands
  
    } // endfor - dirs

  return;

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
void Residual_Smoothing(Rte2D_Quad_Block &SolnBlk,
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
				   Rte2D_Input_Parameters &IP,
                                   int &number_refinement_criteria,
                                   Rte2D_Quad_Block &SolnBlk) {

  int NUM_VAR_RTE2D = SolnBlk.NumVar();
  int i, j;
  double grad_G_x, 
    grad_G_y, 
    grad_G_abs, 
    grad_G_criteria, 
    grad_G_criteria_max;
  int refinement_criteria_number;

  /* Set the number of refinement criteria to be used:
     (1) refinement criteria #1 based on the gradient of the intensity field; */
  number_refinement_criteria = IP.Number_of_Refinement_Criteria;

  /* Initialize the refinement criteria for the block. */
    
  grad_G_criteria_max = ZERO;

  /* Calculate the refinement criteria for each cell of the 
     computational mesh and assign the maximum value for
     all cells as the refinement criteria for the solution 
     block. */

  for ( j  = SolnBlk.JCl-1 ; j <= SolnBlk.JCu+1 ; ++j ) {
    for ( i = SolnBlk.ICl-1 ; i <= SolnBlk.ICu+1 ; ++i ) {
      // Reconstruct the solution within the cell.
      //Linear_Reconstruction_GreenGauss(SolnBlk, i, j, LIMITER_UNLIMITED);
      Linear_Reconstruction_LeastSquares_2(SolnBlk, i, j, LIMITER_UNLIMITED);
      //Linear_Reconstruction_LeastSquares_2(SolnBlk, i, j, LIMITER_UNLIMITED);

      if (SolnBlk.Grid.Cell[i][j].A > ZERO) {
	// Evaluate refinement criteria #1 based on the gradient
	// of the intensity field.
	grad_G_x = SolnBlk.dUdx[i][j].G(SolnBlk.M[i][j]);
	grad_G_y = SolnBlk.dUdy[i][j].G(SolnBlk.M[i][j]);
	grad_G_abs = sqrt(sqr(grad_G_x) + sqr(grad_G_y));
	grad_G_criteria = sqrt(SolnBlk.Grid.Cell[i][j].A)*grad_G_abs / 
	                  fabs(SolnBlk.U[i][j].G(SolnBlk.M[i][j]));
	grad_G_criteria_max = max(grad_G_criteria_max, grad_G_criteria);
	      
      } /* endif */

    } /* endfor */
  } /* endfor */


  // Return the refinement criteria.
  refinement_criteria_number = 0;
  refinement_criteria[refinement_criteria_number] = grad_G_criteria_max;
  refinement_criteria_number++;

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
void Fix_Refined_Block_Boundaries(Rte2D_Quad_Block SolnBlk,
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
	  SolnBlk.U[i][SolnBlk.JCu] = 
	    (SolnBlk.Grid.Cell[i][SolnBlk.JCu].A*SolnBlk.Sp[i][SolnBlk.JCu] /
	     (SolnBlk.Grid.area(i, SolnBlk.JCu)*SolnBlk.Sp_c(i, SolnBlk.JCu))) *
	    SolnBlk.U[i][SolnBlk.JCu];
	  SolnBlk.M[i][SolnBlk.JCu].GetState(SolnBlk.Grid.centroid(i,SolnBlk.JCu));
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
	  SolnBlk.U[i][SolnBlk.JCl] = 
	    (SolnBlk.Grid.Cell[i][SolnBlk.JCl].A*SolnBlk.Sp[i][SolnBlk.JCl] /
             (SolnBlk.Grid.area(i, SolnBlk.JCl)*SolnBlk.Sp_c(i, SolnBlk.JCl))) * 
	    SolnBlk.U[i][SolnBlk.JCl];
	  SolnBlk.M[i][SolnBlk.JCl].GetState(SolnBlk.Grid.centroid(i,SolnBlk.JCl));
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
	  SolnBlk.U[SolnBlk.ICu][j] = 
	    (SolnBlk.Grid.Cell[SolnBlk.ICu][j].A*SolnBlk.Sp[SolnBlk.ICu][j] /
             (SolnBlk.Grid.area(SolnBlk.ICu, j)*SolnBlk.Sp_c(SolnBlk.ICu, j))) * 
	    SolnBlk.U[SolnBlk.ICu][j];
	  SolnBlk.M[SolnBlk.ICu][j].GetState(SolnBlk.Grid.centroid(SolnBlk.ICu,j));
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
	  SolnBlk.U[SolnBlk.ICl][j] = 
	    (SolnBlk.Grid.Cell[SolnBlk.ICl][j].A*SolnBlk.Sp[SolnBlk.ICl][j] /
             (SolnBlk.Grid.area(SolnBlk.ICl, j)*SolnBlk.Sp_c(SolnBlk.ICl, j))) * 
	    SolnBlk.U[SolnBlk.ICl][j];
	  SolnBlk.M[SolnBlk.ICl][j].GetState(SolnBlk.Grid.centroid(SolnBlk.ICl,j));
       } /* endfor */
    } /* endif */

    /* Reset the boundary condition types at the block boundaries. */
 
    Set_BCs(SolnBlk.Grid);

    /* Recompute the exterior nodes for the block quadrilateral mesh. */

    Update_Exterior_Nodes(SolnBlk.Grid);

    /* Recompute the cells for the block quadrilateral mesh. */

    Update_Cells(SolnBlk.Grid);

    /* Recompute the grid 2D to quasi-3D scaling parameters*/

    SolnBlk.ScaleGridTo3D();

}

/********************************************************
 * Routine: Unfix_Refined_Block_Boundaries              *
 *                                                      *
 * Returns the adjusted the locations of the boundary   *
 * nodes of a solution block to their original          *
 * unmodified positions.                                *
 *                                                      *
 ********************************************************/
void Unfix_Refined_Block_Boundaries(Rte2D_Quad_Block SolnBlk) {

    int i, j;
    double sp_l, sp_r, sp_m, ds_ratio, dl, dr;
 
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
	  SolnBlk.U[i][SolnBlk.JCu] = 
	    (SolnBlk.Grid.Cell[i][SolnBlk.JCu].A*SolnBlk.Sp[i][SolnBlk.JCu] /
             (SolnBlk.Grid.area(i, SolnBlk.JCu)*SolnBlk.Sp_c(i, SolnBlk.JCu))) *
	    SolnBlk.U[i][SolnBlk.JCu];
	  SolnBlk.M[i][SolnBlk.JCu].GetState(SolnBlk.Grid.centroid(i,SolnBlk.JCu));
       } /* endfor */
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
	  SolnBlk.U[i][SolnBlk.JCl] = 
	    (SolnBlk.Grid.Cell[i][SolnBlk.JCl].A*SolnBlk.Sp[i][SolnBlk.JCl] /
             (SolnBlk.Grid.area(i, SolnBlk.JCl)*SolnBlk.Sp_c(i, SolnBlk.JCl))) *
	    SolnBlk.U[i][SolnBlk.JCl];
	  SolnBlk.M[i][SolnBlk.JCl].GetState(SolnBlk.Grid.centroid(i,SolnBlk.JCl));
       } /* endfor */
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
	  SolnBlk.U[SolnBlk.ICu][j] = 
	    (SolnBlk.Grid.Cell[SolnBlk.ICu][j].A*SolnBlk.Sp[SolnBlk.ICu][j] /
             (SolnBlk.Grid.area(SolnBlk.ICu, j)*SolnBlk.Sp_c(SolnBlk.ICu, j))) *
	    SolnBlk.U[SolnBlk.ICu][j];
	  SolnBlk.M[SolnBlk.ICu][j].GetState(SolnBlk.Grid.centroid(SolnBlk.ICu,j));
       } /* endfor */
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
	  SolnBlk.U[SolnBlk.ICl][j] = 
	    (SolnBlk.Grid.Cell[SolnBlk.ICl][j].A*SolnBlk.Sp[SolnBlk.ICl][j] /
             (SolnBlk.Grid.area(SolnBlk.ICl, j)*SolnBlk.Sp_c(SolnBlk.ICl, j))) *
	    SolnBlk.U[SolnBlk.ICl][j];
	  SolnBlk.M[SolnBlk.ICl][j].GetState(SolnBlk.Grid.centroid(SolnBlk.ICl,j));
       } /* endfor */
    } /* endif */

    /* Reset the boundary condition types at the block boundaries. */
 
    Set_BCs(SolnBlk.Grid);

    /* Recompute the exterior nodes for the block quadrilateral mesh. */

    Update_Exterior_Nodes(SolnBlk.Grid);

    /* Recompute the cells for the block quadrilateral mesh. */

    Update_Cells(SolnBlk.Grid);

    /* Recompute the grid 2D to quasi-3D scaling parameters*/

    SolnBlk.ScaleGridTo3D();


}

/****************************************************************
 * Routine: Apply_Boundary_Flux_Corrections                     *
 *                                                              *
 * Apply flux corrections at boundaries of the solution         *
 * block to ensure that the scheme is conservative at           *
 * boundaries with mesh resolution changes.                     *
 *                                                              *
 ****************************************************************/
void Apply_Boundary_Flux_Corrections(Rte2D_Quad_Block SolnBlk,
                                     const int Number_Neighbours_North_Boundary,
                                     const int Number_Neighbours_South_Boundary,
                                     const int Number_Neighbours_East_Boundary,
                                     const int Number_Neighbours_West_Boundary) {

    int i, j;
 
    /* Correct the fluxes at the north boundary as required. */

    if (Number_Neighbours_North_Boundary == 2) {
       for ( i = SolnBlk.ICl ; i <= SolnBlk.ICu ; ++i ) {
          SolnBlk.dUdt[i][SolnBlk.JCu][0] -= 
	    SolnBlk.FluxN[i]/(SolnBlk.Grid.Cell[i][SolnBlk.JCu].A*SolnBlk.Sp[i][SolnBlk.JCu]);
       } /* endfor */
    } /* endif */

    /* Correct the fluxes at the south boundary as required. */

    if (Number_Neighbours_South_Boundary == 2) {
       for ( i = SolnBlk.ICl ; i <= SolnBlk.ICu ; ++i ) {
          SolnBlk.dUdt[i][SolnBlk.JCl][0] -= 
	    SolnBlk.FluxS[i]/(SolnBlk.Grid.Cell[i][SolnBlk.JCl].A*SolnBlk.Sp[i][SolnBlk.JCl]);
       } /* endfor */
    } /* endif */

    /* Correct the fluxes at the east boundary as required. */

    if (Number_Neighbours_East_Boundary == 2) {
       for ( j= SolnBlk.JCl ; j <= SolnBlk.JCu ; ++j ) {
          SolnBlk.dUdt[SolnBlk.ICu][j][0] -= 
	    SolnBlk.FluxE[j]/(SolnBlk.Grid.Cell[SolnBlk.ICu][j].A*SolnBlk.Sp[SolnBlk.ICu][j]);
       } /* endfor */
    } /* endif */

    /* Correct the fluxes at the west boundary as required. */

    if (Number_Neighbours_West_Boundary == 2) {
       for ( j= SolnBlk.JCl ; j <= SolnBlk.JCu ; ++j ) {
          SolnBlk.dUdt[SolnBlk.ICl][j][0] -= 
	    SolnBlk.FluxW[j]/(SolnBlk.Grid.Cell[SolnBlk.ICl][j].A*SolnBlk.Sp[SolnBlk.ICl][j]);
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
void Apply_Boundary_Flux_Corrections_Multistage_Explicit(Rte2D_Quad_Block SolnBlk,
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
	     (SolnBlk.Grid.Cell[i][SolnBlk.JCu].A*SolnBlk.Sp[i][SolnBlk.JCu]);
       } /* endfor */
    } /* endif */

    /* Correct the fluxes at the south boundary as required. */

    if (Number_Neighbours_South_Boundary == 2) {
       for ( i = SolnBlk.ICl ; i <= SolnBlk.ICu ; ++i ) {
          SolnBlk.dUdt[i][SolnBlk.JCl][k_residual] -= 
             (CFL_Number*SolnBlk.dt[i][SolnBlk.JCl])*SolnBlk.FluxS[i]/
             (SolnBlk.Grid.Cell[i][SolnBlk.JCl].A*SolnBlk.Sp[i][SolnBlk.JCl]);
       } /* endfor */
    } /* endif */

    /* Correct the fluxes at the east boundary as required. */

    if (Number_Neighbours_East_Boundary == 2) {
       for ( j= SolnBlk.JCl ; j <= SolnBlk.JCu ; ++j ) {
          SolnBlk.dUdt[SolnBlk.ICu][j][k_residual] -= 
             (CFL_Number*SolnBlk.dt[SolnBlk.ICu][j])*SolnBlk.FluxE[j]/
             (SolnBlk.Grid.Cell[SolnBlk.ICu][j].A*SolnBlk.Sp[SolnBlk.ICu][j]);
       } /* endfor */
    } /* endif */

    /* Correct the fluxes at the west boundary as required. */

    if (Number_Neighbours_West_Boundary == 2) {
       for ( j= SolnBlk.JCl ; j <= SolnBlk.JCu ; ++j ) {
          SolnBlk.dUdt[SolnBlk.ICl][j][k_residual] -= 
             (CFL_Number*SolnBlk.dt[SolnBlk.ICl][j])*SolnBlk.FluxW[j]/
             (SolnBlk.Grid.Cell[SolnBlk.ICl][j].A*SolnBlk.Sp[SolnBlk.ICl][j]);
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
int dUdt_Residual_Evaluation(Rte2D_Quad_Block &SolnBlk,
			     Rte2D_Input_Parameters &Input_Parameters){


  int JCl_overlap = 0;
  int JCu_overlap = 0;
  int ICu_overlap = 0;
  int ICl_overlap = 0;
  
  // Modifications for NKS overlap 
  if(Input_Parameters.NKS_IP.GMRES_Overlap > 0 ){	
    if (SolnBlk.Grid.BCtypeS[SolnBlk.ICl] == BC_NONE)  JCl_overlap = Input_Parameters.NKS_IP.GMRES_Overlap;
    if (SolnBlk.Grid.BCtypeN[SolnBlk.ICu] == BC_NONE)  JCu_overlap = Input_Parameters.NKS_IP.GMRES_Overlap;
    if (SolnBlk.Grid.BCtypeE[SolnBlk.JCu] == BC_NONE)  ICu_overlap = Input_Parameters.NKS_IP.GMRES_Overlap;
    if (SolnBlk.Grid.BCtypeW[SolnBlk.JCl] == BC_NONE)  ICl_overlap = Input_Parameters.NKS_IP.GMRES_Overlap;    
  }

    int i, j;
    Vector2D dX;
    Rte2D_State Ul, Ur;
    Rte2D_State Flux;

    /* Perform the linear reconstruction within each cell
       of the computational grid for this stage. */
    
    switch(Input_Parameters.i_Reconstruction) {
    case RECONSTRUCTION_GREEN_GAUSS :
      Linear_Reconstruction_GreenGauss(SolnBlk,
                                       Input_Parameters.i_Limiter); 
      if(Input_Parameters.Axisymmetric)
	Linear_Reconstruction_GreenGauss_Angular(SolnBlk, LIMITER_ONE); 
      break;
    case RECONSTRUCTION_LEAST_SQUARES :
      Linear_Reconstruction_LeastSquares(SolnBlk,
					 Input_Parameters.i_Limiter);
      if(Input_Parameters.Axisymmetric)
	Linear_Reconstruction_LeastSquares_Angular(SolnBlk, LIMITER_ONE); 
      break;
    default:
      Linear_Reconstruction_LeastSquares(SolnBlk,
                                         Input_Parameters.i_Limiter);
      break;
    } /* endswitch */

    /* Evaluate the time rate of change of the solution
       (i.e., the solution residuals) using a second-order
       limited upwind scheme with a variety of flux functions. */
    
    // Add i-direction (zeta-direction) fluxes.
    for ( j  = SolnBlk.JCl-1-JCl_overlap ; j <= SolnBlk.JCu+1+JCu_overlap ; ++j ) {

      SolnBlk.dUdt[SolnBlk.ICl-1-ICl_overlap][j][0].Zero();
          
      for ( i = SolnBlk.ICl-1-ICl_overlap; i <= SolnBlk.ICu+ICu_overlap ; ++i ) {

	SolnBlk.dUdt[i+1][j][0].Zero();
    
	if ( j > SolnBlk.JCl-1-JCl_overlap && j < SolnBlk.JCu+1+JCu_overlap ) {
    
	  /* Evaluate the cell interface i-direction fluxes. */
	  
	  if (i == SolnBlk.ICl-1 && 
	      (SolnBlk.Grid.BCtypeW[j] == BC_REFLECTION ||
	       SolnBlk.Grid.BCtypeW[j] == BC_GRAY_WALL)) {

	    dX = SolnBlk.Grid.xfaceW(i+1, j)-SolnBlk.Grid.Cell[i+1][j].Xc;
	    Ur = SolnBlk.U[i+1][j] + 
	      (SolnBlk.phi[i+1][j]^SolnBlk.dUdx[i+1][j])*dX.x +
	      (SolnBlk.phi[i+1][j]^SolnBlk.dUdy[i+1][j])*dX.y;

	    if (SolnBlk.Grid.BCtypeW[j] == BC_REFLECTION) {
	      Ul = Reflect(Ur, SolnBlk.Grid.nfaceW(i+1, j));
	    } else if (SolnBlk.Grid.BCtypeW[j] == BC_GRAY_WALL) {
	      Ul = Gray_Wall(Ur, SolnBlk.M[0][0], SolnBlk.Grid.nfaceW(i+1, j), 
			SolnBlk.WestWallTemp, SolnBlk.WestWallEmiss );
	    }/* endif */

	  } else if (i == SolnBlk.ICu && 
		     (SolnBlk.Grid.BCtypeE[j] == BC_REFLECTION ||
		      SolnBlk.Grid.BCtypeE[j] == BC_GRAY_WALL )) {

	    dX = SolnBlk.Grid.xfaceE(i, j)-SolnBlk.Grid.Cell[i][j].Xc;
	    Ul = SolnBlk.U[i][j] + 
	      (SolnBlk.phi[i][j]^SolnBlk.dUdx[i][j])*dX.x +
	      (SolnBlk.phi[i][j]^SolnBlk.dUdy[i][j])*dX.y;

	    if (SolnBlk.Grid.BCtypeE[j] == BC_REFLECTION) {
	      Ur = Reflect(Ul, SolnBlk.Grid.nfaceE(i, j));
	    } else if (SolnBlk.Grid.BCtypeE[j] == BC_GRAY_WALL) {
	      Ur = Gray_Wall(Ul, SolnBlk.M[0][0], SolnBlk.Grid.nfaceE(i, j), 
			SolnBlk.EastWallTemp, SolnBlk.EastWallEmiss );
	    }/* endif */

	  } else {            
	    dX = SolnBlk.Grid.xfaceE(i, j)-SolnBlk.Grid.Cell[i][j].Xc;
	    Ul = SolnBlk.U[i][j] + 
	      (SolnBlk.phi[i][j]^SolnBlk.dUdx[i][j])*dX.x +
	      (SolnBlk.phi[i][j]^SolnBlk.dUdy[i][j])*dX.y;
	    dX = SolnBlk.Grid.xfaceW(i+1, j)-SolnBlk.Grid.Cell[i+1][j].Xc;
	    Ur = SolnBlk.U[i+1][j] + 
	      (SolnBlk.phi[i+1][j]^SolnBlk.dUdx[i+1][j])*dX.x +
	      (SolnBlk.phi[i+1][j]^SolnBlk.dUdy[i+1][j])*dX.y;
	  } /* endif */

	  Flux = Flux_n(Ul, Ur, SolnBlk.Grid.nfaceE(i, j));

	  /* Evaluate cell-averaged solution changes. */
	  
	  SolnBlk.dUdt[i][j][0] -= 
	    Flux*SolnBlk.Grid.lfaceE(i, j)*SolnBlk.SpE[i][j]/
	    (SolnBlk.Grid.Cell[i][j].A*SolnBlk.Sp[i][j]);
	  SolnBlk.dUdt[i+1][j][0] += 
	    Flux*SolnBlk.Grid.lfaceW(i+1, j)*SolnBlk.SpW[i+1][j]/
	    (SolnBlk.Grid.Cell[i+1][j].A*SolnBlk.Sp[i+1][j]);

	  /* Include general source term */
	  
	  SolnBlk.dUdt[i][j][0] += SolnBlk.U[i][j].S(SolnBlk.M[i][j]);

	  /* Include axisymmetric source terms as required. */

	  if (SolnBlk.Axisymmetric) {
	    SolnBlk.dUdt[i][j][0] += 
	      SolnBlk.U[i][j].Sa( SolnBlk.dUdpsi[i][j], 
				  SolnBlk.phi_psi[i][j],
				  SolnBlk.Axisymmetric) / SolnBlk.Sp[i][j];
	  } /* endif */

	  /* Save west and east face boundary flux. */
	  
	  if (i == SolnBlk.ICl-1) {
	    SolnBlk.FluxW[j] = -Flux*SolnBlk.Grid.lfaceW(i+1, j)*SolnBlk.SpW[i+1][j];
	  } else if (i == SolnBlk.ICu) {
	    SolnBlk.FluxE[j] = Flux*SolnBlk.Grid.lfaceE(i, j)*SolnBlk.SpE[i][j];
	  } /* endif */ 

	} /* endif */
      } /* endfor */
      

      if ( j > SolnBlk.JCl-1-JCl_overlap && j < SolnBlk.JCu+1+JCu_overlap ) {  
	SolnBlk.dUdt[SolnBlk.ICl-1-ICl_overlap][j][0].Zero();      
	SolnBlk.dUdt[SolnBlk.ICu+1+ICu_overlap][j][0].Zero();
      } /* endif */
    } /* endfor */

    // Add j-direction (eta-direction) fluxes.
    for ( i = SolnBlk.ICl-ICl_overlap ; i <= SolnBlk.ICu+ICu_overlap ; ++i ) {
      for ( j  = SolnBlk.JCl-1-JCl_overlap ; j <= SolnBlk.JCu+JCu_overlap ; ++j ) {
	
	/* Evaluate the cell interface j-direction fluxes. */
	
	if (j == SolnBlk.JCl-1 && 
	    (SolnBlk.Grid.BCtypeS[i] == BC_REFLECTION ||
	     SolnBlk.Grid.BCtypeS[i] == BC_GRAY_WALL )) {

	  dX = SolnBlk.Grid.xfaceS(i, j+1)-SolnBlk.Grid.Cell[i][j+1].Xc;
	  Ur = SolnBlk.U[i][j+1] +
	    (SolnBlk.phi[i][j+1]^SolnBlk.dUdx[i][j+1])*dX.x +
	    (SolnBlk.phi[i][j+1]^SolnBlk.dUdy[i][j+1])*dX.y;

	  if (SolnBlk.Grid.BCtypeS[i] == BC_REFLECTION ) {
	    Ul = Reflect(Ur, SolnBlk.Grid.nfaceS(i, j+1));
	  } else if (SolnBlk.Grid.BCtypeS[i] == BC_GRAY_WALL) {
	    Ul = Gray_Wall(Ur, SolnBlk.M[0][0], SolnBlk.Grid.nfaceS(i, j+1), 
		      SolnBlk.SouthWallTemp, SolnBlk.SouthWallEmiss );
	  } /* endif */

	} else if (j == SolnBlk.JCu && 
		   (SolnBlk.Grid.BCtypeN[i] == BC_REFLECTION ||
		    SolnBlk.Grid.BCtypeN[i] == BC_GRAY_WALL )) {

	  dX = SolnBlk.Grid.xfaceN(i, j)-SolnBlk.Grid.Cell[i][j].Xc;
	  Ul = SolnBlk.U[i][j] + 
	    (SolnBlk.phi[i][j]^SolnBlk.dUdx[i][j])*dX.x +
	    (SolnBlk.phi[i][j]^SolnBlk.dUdy[i][j])*dX.y;

	  if (SolnBlk.Grid.BCtypeN[i] == BC_REFLECTION) {
	    Ur = Reflect(Ul, SolnBlk.Grid.nfaceN(i, j));
	  } else if (SolnBlk.Grid.BCtypeN[i] == BC_GRAY_WALL) {
	    Ur = Gray_Wall(Ul, SolnBlk.M[0][0], SolnBlk.Grid.nfaceN(i, j), 
			   SolnBlk.NorthWallTemp, SolnBlk.NorthWallEmiss );
	  } /* endif */

	} else {
	  dX = SolnBlk.Grid.xfaceN(i, j)-SolnBlk.Grid.Cell[i][j].Xc;
	  Ul = SolnBlk.U[i][j] + 
	    (SolnBlk.phi[i][j]^SolnBlk.dUdx[i][j])*dX.x +
	    (SolnBlk.phi[i][j]^SolnBlk.dUdy[i][j])*dX.y;
	  dX = SolnBlk.Grid.xfaceS(i, j+1)-SolnBlk.Grid.Cell[i][j+1].Xc;
	  Ur = SolnBlk.U[i][j+1] +
	    (SolnBlk.phi[i][j+1]^SolnBlk.dUdx[i][j+1])*dX.x +
	    (SolnBlk.phi[i][j+1]^SolnBlk.dUdy[i][j+1])*dX.y;
	} /* endif */
	
	Flux = Flux_n(Ul, Ur, SolnBlk.Grid.nfaceN(i, j));
	
          /* Evaluate cell-averaged solution changes. */
	
	SolnBlk.dUdt[i][j][0] -= 
	  Flux*SolnBlk.Grid.lfaceN(i, j)*SolnBlk.SpN[i][j]/
	  (SolnBlk.Grid.Cell[i][j].A*SolnBlk.Sp[i][j]);
	SolnBlk.dUdt[i][j+1][0] += 
	  Flux*SolnBlk.Grid.lfaceS(i, j+1)*SolnBlk.SpS[i][j+1]/
	  (SolnBlk.Grid.Cell[i][j+1].A*SolnBlk.Sp[i][j+1]);

	/* Save south and north face boundary flux. */
	
	if (j == SolnBlk.JCl-1) {
	  SolnBlk.FluxS[i] = -Flux*SolnBlk.Grid.lfaceS(i, j+1)*SolnBlk.SpS[i][j+1];
	} else if (j == SolnBlk.JCu) {
	  SolnBlk.FluxN[i] = Flux*SolnBlk.Grid.lfaceN(i, j)*SolnBlk.SpN[i][j];
	} /* endif */
	
      } /* endfor */
 
      SolnBlk.dUdt[i][SolnBlk.JCl-1-JCl_overlap][0].Zero();
      SolnBlk.dUdt[i][SolnBlk.JCu+1+JCu_overlap][0].Zero();
    } /* endfor */
    
    /* Residual successfully evaluated. */
    return 0;

}

/********************************************************
 * Routine: dUdt_Multistage_Explicit                    *
 *                                                      *
 * This routine determines the solution residuals for a *
 * given stage of a variety of multi-stage explicit     *
 * time integration schemes for a given solution block. *
 *                                                      *
 ********************************************************/
int dUdt_Multistage_Explicit(Rte2D_Quad_Block &SolnBlk,
                             const int i_stage,
                             Rte2D_Input_Parameters &Input_Parameters) {

    int i, j, k_residual;
    double omega;
    Vector2D dX;
    Rte2D_State Ul, Ur;
    Rte2D_State Flux;

    /* Evaluate the solution residual for stage 
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
    
    /* Perform the linear reconstruction within each cell
       of the computational grid for this stage. */
    
    switch(Input_Parameters.i_Reconstruction) {
    case RECONSTRUCTION_GREEN_GAUSS :
      Linear_Reconstruction_GreenGauss(SolnBlk,
                                       Input_Parameters.i_Limiter);    
      if(Input_Parameters.Axisymmetric)
	Linear_Reconstruction_GreenGauss_Angular(SolnBlk, LIMITER_ONE); 
      break;
    case RECONSTRUCTION_LEAST_SQUARES :
      Linear_Reconstruction_LeastSquares(SolnBlk,
					 Input_Parameters.i_Limiter);
      if(Input_Parameters.Axisymmetric)
	Linear_Reconstruction_LeastSquares_Angular(SolnBlk, LIMITER_ONE); 
      break;
    default:
      Linear_Reconstruction_LeastSquares(SolnBlk,
                                         Input_Parameters.i_Limiter);
      break;
    } /* endswitch */

    /* Evaluate the time rate of change of the solution
       (i.e., the solution residuals) using a second-order
       limited upwind scheme with a variety of flux functions. */
    
    // Add i-direction (zeta-direction) fluxes.
    for ( j  = SolnBlk.JCl-1 ; j <= SolnBlk.JCu+1 ; ++j ) {
       if ( i_stage == 1 ) {
          SolnBlk.Uo[SolnBlk.ICl-1][j] = SolnBlk.U[SolnBlk.ICl-1][j];
          SolnBlk.dUdt[SolnBlk.ICl-1][j][k_residual].Zero();
       } else {
          SolnBlk.dUdt[SolnBlk.ICl-1][j][k_residual].Zero();
       } /* endif */
    
       for ( i = SolnBlk.ICl-1 ; i <= SolnBlk.ICu ; ++i ) {
          if ( i_stage == 1 ) {
              SolnBlk.Uo[i+1][j] = SolnBlk.U[i+1][j];
              SolnBlk.dUdt[i+1][j][k_residual].Zero();
          } else if ( j > SolnBlk.JCl-1 && j < SolnBlk.JCu+1 ) {
              switch(Input_Parameters.i_Time_Integration) {
                case TIME_STEPPING_EXPLICIT_PREDICTOR_CORRECTOR :
                  //SolnBlk.dUdt[i+1][j][k_residual] = 
                  //   SolnBlk.dUdt[i+1][j][k_residual];
                  break;
                case TIME_STEPPING_EXPLICIT_RUNGE_KUTTA :
                  if (Input_Parameters.N_Stage == 2) {
		     //SolnBlk.dUdt[i+1][j][k_residual] = 
		     //   SolnBlk.dUdt[i+1][j][k_residual];
                  } else if (Input_Parameters.N_Stage == 4 && i_stage == 4) {
                     SolnBlk.dUdt[i+1][j][k_residual] = 
                        SolnBlk.dUdt[i+1][j][0] + 
                        TWO*SolnBlk.dUdt[i+1][j][1] +
                        TWO*SolnBlk.dUdt[i+1][j][2];
                  } else {
                     SolnBlk.dUdt[i+1][j][k_residual].Zero();
                  } /* endif */
                  break;
                case TIME_STEPPING_MULTISTAGE_OPTIMAL_SMOOTHING :
                  SolnBlk.dUdt[i+1][j][k_residual].Zero();
                  break;
                default:
                  SolnBlk.dUdt[i+1][j][k_residual].Zero();
                  break;
              } /* endswitch */
          } /* endif */
    
          if ( j > SolnBlk.JCl-1 && j < SolnBlk.JCu+1 ) {
    
             /* Evaluate the cell interface i-direction fluxes. */
    
	     if (i == SolnBlk.ICl-1 && 
                 (SolnBlk.Grid.BCtypeW[j] == BC_REFLECTION ||
		  SolnBlk.Grid.BCtypeW[j] == BC_GRAY_WALL)) {
	       dX = SolnBlk.Grid.xfaceW(i+1, j)-SolnBlk.Grid.Cell[i+1][j].Xc;
	       Ur = SolnBlk.U[i+1][j] + 
		 (SolnBlk.phi[i+1][j]^SolnBlk.dUdx[i+1][j])*dX.x +
		 (SolnBlk.phi[i+1][j]^SolnBlk.dUdy[i+1][j])*dX.y;

	       if (SolnBlk.Grid.BCtypeW[j] == BC_REFLECTION) {
                 Ul = Reflect(Ur, SolnBlk.Grid.nfaceW(i+1, j));
	       } else if (SolnBlk.Grid.BCtypeW[j] == BC_GRAY_WALL) {
		 Ul = Gray_Wall(Ur, SolnBlk.M[0][0], SolnBlk.Grid.nfaceW(i+1, j), 
			   SolnBlk.WestWallTemp, SolnBlk.WestWallEmiss );
	       }/* endif */

             } else if (i == SolnBlk.ICu && 
                        (SolnBlk.Grid.BCtypeE[j] == BC_REFLECTION ||
			 SolnBlk.Grid.BCtypeE[j] == BC_GRAY_WALL)) {

	       dX = SolnBlk.Grid.xfaceE(i, j)-SolnBlk.Grid.Cell[i][j].Xc;
	       Ul = SolnBlk.U[i][j] + 
		 (SolnBlk.phi[i][j]^SolnBlk.dUdx[i][j])*dX.x +
		 (SolnBlk.phi[i][j]^SolnBlk.dUdy[i][j])*dX.y;

	       if (SolnBlk.Grid.BCtypeE[j] == BC_REFLECTION) {
		 Ur = Reflect(Ul, SolnBlk.Grid.nfaceE(i, j));
	       } else if (SolnBlk.Grid.BCtypeE[j] == BC_GRAY_WALL) {
		 Ur = Gray_Wall(Ul, SolnBlk.M[0][0], SolnBlk.Grid.nfaceE(i, j), 
			   SolnBlk.EastWallTemp, SolnBlk.EastWallEmiss );
	       }/* endif */

             } else {            
  	        dX = SolnBlk.Grid.xfaceE(i, j)-SolnBlk.Grid.Cell[i][j].Xc;
                Ul = SolnBlk.U[i][j] + 
                     (SolnBlk.phi[i][j]^SolnBlk.dUdx[i][j])*dX.x +
	             (SolnBlk.phi[i][j]^SolnBlk.dUdy[i][j])*dX.y;
	        dX = SolnBlk.Grid.xfaceW(i+1, j)-SolnBlk.Grid.Cell[i+1][j].Xc;
                Ur = SolnBlk.U[i+1][j] + 
  	             (SolnBlk.phi[i+1][j]^SolnBlk.dUdx[i+1][j])*dX.x +
	             (SolnBlk.phi[i+1][j]^SolnBlk.dUdy[i+1][j])*dX.y;
  	     } /* endif */

	     Flux = Flux_n(Ul, Ur, SolnBlk.Grid.nfaceE(i, j));
    
             /* Evaluate cell-averaged solution changes. */
    
             SolnBlk.dUdt[i][j][k_residual] -= 
                (Input_Parameters.CFL_Number*SolnBlk.dt[i][j])*
                Flux*SolnBlk.Grid.lfaceE(i, j)*SolnBlk.SpE[i][j]/
                (SolnBlk.Grid.Cell[i][j].A*SolnBlk.Sp[i][j]);
             SolnBlk.dUdt[i+1][j][k_residual] += 
                (Input_Parameters.CFL_Number*SolnBlk.dt[i+1][j])*
                Flux*SolnBlk.Grid.lfaceW(i+1, j)*SolnBlk.SpW[i+1][j]/
                (SolnBlk.Grid.Cell[i+1][j].A*SolnBlk.Sp[i+1][j]);

             /* Include general source term. */

	     SolnBlk.dUdt[i][j][k_residual] += 
	       (Input_Parameters.CFL_Number*SolnBlk.dt[i][j])*
	       SolnBlk.U[i][j].S(SolnBlk.M[i][j]);

             /* Include axisymmetric source terms as required. */

	     if (SolnBlk.Axisymmetric) {
               SolnBlk.dUdt[i][j][k_residual] += 
		 (Input_Parameters.CFL_Number*SolnBlk.dt[i][j]) *
		 SolnBlk.U[i][j].Sa( SolnBlk.dUdpsi[i][j], 
				     SolnBlk.phi_psi[i][j],
				     SolnBlk.Axisymmetric ) / SolnBlk.Sp[i][j];
             } /* endif */

             /* Save west and east face boundary flux. */

             if (i == SolnBlk.ICl-1) {
                SolnBlk.FluxW[j] = -Flux*SolnBlk.Grid.lfaceW(i+1, j)*SolnBlk.SpW[i+1][j];
             } else if (i == SolnBlk.ICu) {
                SolnBlk.FluxE[j] = Flux*SolnBlk.Grid.lfaceE(i, j)*SolnBlk.SpE[i][j];
             } /* endif */ 

          } /* endif */
       } /* endfor */
    
       if ( j > SolnBlk.JCl-1 && j < SolnBlk.JCu+1 ) {
          SolnBlk.dUdt[SolnBlk.ICl-1][j][k_residual].Zero();
          SolnBlk.dUdt[SolnBlk.ICu+1][j][k_residual].Zero();
       } /* endif */
    } /* endfor */
    
    // Add j-direction (eta-direction) fluxes.
    for ( i = SolnBlk.ICl ; i <= SolnBlk.ICu ; ++i ) {
       for ( j  = SolnBlk.JCl-1 ; j <= SolnBlk.JCu ; ++j ) {
    
          /* Evaluate the cell interface j-direction fluxes. */
         
	  if (j == SolnBlk.JCl-1 && 
              (SolnBlk.Grid.BCtypeS[i] == BC_REFLECTION ||
	       SolnBlk.Grid.BCtypeS[i] == BC_GRAY_WALL)) {
	     dX = SolnBlk.Grid.xfaceS(i, j+1)-SolnBlk.Grid.Cell[i][j+1].Xc;
	     Ur = SolnBlk.U[i][j+1] +
	          (SolnBlk.phi[i][j+1]^SolnBlk.dUdx[i][j+1])*dX.x +
	          (SolnBlk.phi[i][j+1]^SolnBlk.dUdy[i][j+1])*dX.y;

	     if (SolnBlk.Grid.BCtypeS[i] == BC_REFLECTION ) {
               Ul = Reflect(Ur, SolnBlk.Grid.nfaceS(i, j+1));
	     } else if (SolnBlk.Grid.BCtypeS[i] == BC_GRAY_WALL) {
	       Ul = Gray_Wall(Ur, SolnBlk.M[0][0], SolnBlk.Grid.nfaceS(i, j+1), 
			 SolnBlk.SouthWallTemp, SolnBlk.SouthWallEmiss );
	     } /* endif */

          } else if (j == SolnBlk.JCu && 
                     (SolnBlk.Grid.BCtypeN[i] == BC_REFLECTION ||
		      SolnBlk.Grid.BCtypeN[i] == BC_GRAY_WALL )) {

	    dX = SolnBlk.Grid.xfaceN(i, j)-SolnBlk.Grid.Cell[i][j].Xc;
	    Ul = SolnBlk.U[i][j] + 
	      (SolnBlk.phi[i][j]^SolnBlk.dUdx[i][j])*dX.x +
	      (SolnBlk.phi[i][j]^SolnBlk.dUdy[i][j])*dX.y;

	    if (SolnBlk.Grid.BCtypeN[i] == BC_REFLECTION) {
	      Ur = Reflect(Ul, SolnBlk.Grid.nfaceN(i, j));
	    } else if (SolnBlk.Grid.BCtypeN[i] == BC_GRAY_WALL) {
	      Ur = Gray_Wall(Ul, SolnBlk.M[0][0], SolnBlk.Grid.nfaceN(i, j), 
			SolnBlk.NorthWallTemp, SolnBlk.NorthWallEmiss );
	    } /* endif */

          } else {
  	     dX = SolnBlk.Grid.xfaceN(i, j)-SolnBlk.Grid.Cell[i][j].Xc;
             Ul = SolnBlk.U[i][j] + 
                  (SolnBlk.phi[i][j]^SolnBlk.dUdx[i][j])*dX.x +
	          (SolnBlk.phi[i][j]^SolnBlk.dUdy[i][j])*dX.y;
	     dX = SolnBlk.Grid.xfaceS(i, j+1)-SolnBlk.Grid.Cell[i][j+1].Xc;
             Ur = SolnBlk.U[i][j+1] +
                  (SolnBlk.phi[i][j+1]^SolnBlk.dUdx[i][j+1])*dX.x +
                  (SolnBlk.phi[i][j+1]^SolnBlk.dUdy[i][j+1])*dX.y;
          } /* endif */

	  Flux = Flux_n(Ul, Ur, SolnBlk.Grid.nfaceN(i, j));
    
          /* Evaluate cell-averaged solution changes. */
    
          SolnBlk.dUdt[i][j][k_residual] -= 
             (Input_Parameters.CFL_Number*SolnBlk.dt[i][j])*
             Flux*SolnBlk.Grid.lfaceN(i, j)*SolnBlk.SpN[i][j]/
             (SolnBlk.Grid.Cell[i][j].A*SolnBlk.Sp[i][j]);
          SolnBlk.dUdt[i][j+1][k_residual] += 
             (Input_Parameters.CFL_Number*SolnBlk.dt[i][j+1])*
             Flux*SolnBlk.Grid.lfaceS(i, j+1)*SolnBlk.SpS[i][j+1]/
             (SolnBlk.Grid.Cell[i][j+1].A*SolnBlk.Sp[i][j+1]);

          /* Save south and north face boundary flux. */

          if (j == SolnBlk.JCl-1) {
             SolnBlk.FluxS[i] = -Flux*SolnBlk.Grid.lfaceS(i, j+1)*SolnBlk.SpS[i][j+1];
          } else if (j == SolnBlk.JCu) {
             SolnBlk.FluxN[i] = Flux*SolnBlk.Grid.lfaceN(i, j)*SolnBlk.SpN[i][j];
          } /* endif */

       } /* endfor */
    
       SolnBlk.dUdt[i][SolnBlk.JCl-1][k_residual].Zero();
       SolnBlk.dUdt[i][SolnBlk.JCu+1][k_residual].Zero();
    } /* endfor */
    
    /* Residual for the stage successfully calculated. */

    return (0);
    
}



/********************************************************
 * Routine: dUdt_Space_March_Flux_Eval                  *
 *                                                      *
 * Evaluates the next upwind face values and the cell   *
 * average value for the space marching procedure.      *
 * The upwind and several 2nd-order and high-res spacial*
 * discretization schemes may be used.                  *
 *                                                      *
 * Note: xd and yd are indexes that point in the        *
 * downstream direction.  eg. xd = 1 -> ray traveling   *
 * in +ve x,  xd = -1 -> ray traveling in -ve x.        *
 * Similarily, yd for y-dir.                            *
 *                                                      *
 ********************************************************/
int dUdt_Space_March_Flux_Eval(Rte2D_Quad_Block &SolnBlk,
			       Rte2D_Input_Parameters &Input_Parameters,
			       const int &i,  // x position index
			       const int &j,  // y position index
			       const int &v,  // band index
			       const int &m,  // polar angle index
			       const int &l,  // azim angle index
			       const int &xd, // i+xd -> downstream x index
			       const int &yd, // i+yd -> downstream y index
			       double* &Ix_f, // array of upstream face intensities along x-dir
			       double &Iy_f)  // upstream face intensitie along y-dir
{
  //
  // declares
  //
  Rte2D_State &U = SolnBlk.U[0][0]; // alias
  double V, Ae, Aw, An, As;         // cell face areas, volumes
  double num, denom;                // numerator and denominator
  double Ix_out, Iy_out;            // outgoing face intensity
  double xf, yf;                    // face coordinates                
  double a_N, a_E, a_S, a_W;        // finite volume coefficients
  double corr;                      // high order correction term
  double wx, wwx, wy, wwy,          // weighting factors for differerncing
         wpsi, wwpsi;               //  -> value = w*value_upstream + (1-w)*value_downstream
  double tau_x, tau_y;              // optical thicknesses
  double dx, dy;                    // average face length
  double vx, vy;                    // convective velocity
  double S;                         // source term (for GM scheme)
  double a, b;                      // GM coefficients
  int i_SpaceMarch_Scheme;          // temporary space marching scheme flag

  // variable index
  int index = U.Index[v][m][l];

  //
  // for high res methods, use the upwind scheme at the boundaries
  //
  if ( ( j == SolnBlk.JCl || j == SolnBlk.JCu ||
	 i == SolnBlk.ICl || i == SolnBlk.ICu ) &&
       ( Input_Parameters.i_SpaceMarch_Scheme == SPACE_MARCH_CLAM || 
	 Input_Parameters.i_SpaceMarch_Scheme == SPACE_MARCH_GM ) ) {
    i_SpaceMarch_Scheme = SPACE_MARCH_UPWIND;
  } else {
    i_SpaceMarch_Scheme = Input_Parameters.i_SpaceMarch_Scheme;
  } /* endif */

  
  //
  // For Axisymmetric DOM, set the value of the first special direction 
  // I_(m,1/2) to that of I_(m,1).
  // See Jendoubi, Lee, and Kim, J Thermo & Heat Trans, v7, p213-219. (1993)
  //
  // Also, set the angular discretization coefficients.
  //    value = w*value_upstream + (1-w)*value_downstream
  // For upwind spacial scheme, use upwinding (=1).
  // For higher order spacial schemes, use central differencing (=0.5).
  //
  if (SolnBlk.Axisymmetric && 
      Input_Parameters.i_RTE_Solver == RTE2D_SOLVER_DOM ) {
    
    // initialization
    if (l==0) SolnBlk.U[i][j].I_half[m][l] = SolnBlk.U[i][j].I[index];
    
    // set diff parameter
    if (Input_Parameters.i_SpaceMarch_Scheme==SPACE_MARCH_UPWIND) 
      wpsi=ONE; 
    else 
      wpsi = HALF;
    wwpsi = ( ONE - wpsi ) / wpsi;
      
  } // endif


    // compute the areas and volumes
    V =  SolnBlk.Grid.Cell[i][j].A* SolnBlk.Sp[i][j]; 
    Ae = SolnBlk.Grid.lfaceE(i, j)*SolnBlk.SpE[i][j];
    Aw = SolnBlk.Grid.lfaceW(i, j)*SolnBlk.SpW[i][j];
    An = SolnBlk.Grid.lfaceN(i, j)*SolnBlk.SpN[i][j];
    As = SolnBlk.Grid.lfaceS(i, j)*SolnBlk.SpS[i][j];

    // flux derivatives
    a_E = SolnBlk.U[i][j].dFndU(SolnBlk.Grid.nfaceE(i,j),m,l)*Ae;
    a_W = SolnBlk.U[i][j].dFndU(SolnBlk.Grid.nfaceW(i,j),m,l)*Aw;
    a_N = SolnBlk.U[i][j].dFndU(SolnBlk.Grid.nfaceN(i,j),m,l)*An;
    a_S = SolnBlk.U[i][j].dFndU(SolnBlk.Grid.nfaceS(i,j),m,l)*As;

    //
    // perform differencing to predict downstream face values
    // for higher order correction term
    // and set differencing switch (0.5 for central, )
    // Note: high resolution schemes use the deffered correction method
    // i.e. a high order correction to the upwind solution is 
    // computed based on the solution from the previous iteration.
    // See :
    // S.G. Rubin and P.K. Khosla, "Polynomial interpolation method for
    // viscous flow calculations," in Journal of Computational Physics,
    // vol 27, pp. 153-168. 1982.
    //
    switch (i_SpaceMarch_Scheme) {

    //------------------------------------------------
    // CLAM high resolution differencing scheme. 
    //------------------------------------------------
    // See: B. Van Leer, "Towards the ultimate conservation      
    // difference scheme. II. Monotonicity anc conservation 
    // combined in a second-order scheme," in Journal       
    // Computational physics, v14, 1974.  
    // Note: this scheme uses the deffered correction method
    case SPACE_MARCH_CLAM:
      xf = SolnBlk.Grid.Cell[i][j].Xc.x + 
	   xd*QUARTER*(SolnBlk.Grid.lfaceN(i,j)+SolnBlk.Grid.lfaceS(i,j));
      Ix_out = CLAM(SolnBlk.Uo[i-xd][j].I[index],   // I upstream node
		    SolnBlk.Uo[i   ][j].I[index],   // I centroid node
		    SolnBlk.Uo[i+xd][j].I[index],   // I downstream node
		    SolnBlk.Grid.Cell[i-xd][j].Xc.x, // x upstream node
		    SolnBlk.Grid.Cell[i   ][j].Xc.x, // x centroid node
		    SolnBlk.Grid.Cell[i+xd][j].Xc.x, // x downstream node
		    xf);                             // x downstream face

      yf = SolnBlk.Grid.Cell[i][j].Xc.y + 
	   yd*QUARTER*(SolnBlk.Grid.lfaceE(i,j)+SolnBlk.Grid.lfaceW(i,j));
      Iy_out = CLAM(SolnBlk.Uo[i][j-yd].I[index],   // I upstream node
		    SolnBlk.Uo[i][j   ].I[index],   // I centroid node
		    SolnBlk.Uo[i][j+yd].I[index],   // I downstream node
		    SolnBlk.Grid.Cell[i][j-yd].Xc.y, // y upstream node
		    SolnBlk.Grid.Cell[i][j   ].Xc.y, // y centroid node
		    SolnBlk.Grid.Cell[i][j+yd].Xc.y, // y downstream node
		    yf);                             // y downstream face

      // high order correction term
      corr = ( max(a_E, 0.0) + max(a_W, 0.0) ) * ( SolnBlk.U[i][j].I[index] - Ix_out ) +
             ( max(a_N, 0.0) + max(a_S, 0.0) ) * ( SolnBlk.U[i][j].I[index] - Iy_out );
      
      // solve upwind solution with correction term
      wx = ONE;
      wy = ONE;
      break;

    //------------------------------------------------
    // General Multidimensional high resolution differencing 
    // scheme for uniform grids. 
    //------------------------------------------------
    // See: D. Balsara, "Fast and accurate discrete ordinates methods
    // for multidimensional radiative heat transfer," in Journal       
    // of Quantitative Spectroscopy, v69, 2001.    
    // Note: this scheme uses the deffered correction method
    case SPACE_MARCH_GM:

      dx = HALF * (SolnBlk.Grid.lfaceN(i,j) + SolnBlk.Grid.lfaceS(i,j) );
      dy = HALF * (SolnBlk.Grid.lfaceE(i,j) + SolnBlk.Grid.lfaceW(i,j) );
      vx = fabs(U.mu[m][l]);
      vy = fabs(U.eta[m][l]);
      if ( vx/dx <= vy/dy ) {
	S = SolnBlk.Uo[i+xd][j].S(SolnBlk.M[i+xd][j],v,m,l) - 
	    SolnBlk.Uo[i+xd][j].dSdU(SolnBlk.M[i+xd][j],v,m,l)*SolnBlk.Uo[i+xd][j].I[index];
	a = vx*(SolnBlk.Uo[i+xd][j].I[index]-SolnBlk.Uo[i][j-yd].I[index]) - dx*S;
	b = (vx-vy*dx/dy) * (SolnBlk.Uo[i][j].I[index] - SolnBlk.Uo[i][j-yd].I[index]);
	Ix_out = SolnBlk.Uo[i][j-yd].I[index] + 0.5*vanalbada(a,b,0.10)/vx;
	Iy_out = SolnBlk.Uo[i][j].I[index];

      } else {
	S = SolnBlk.Uo[i][j+yd].S(SolnBlk.M[i][j+yd],v,m,l) - 
	    SolnBlk.Uo[i][j+yd].dSdU(SolnBlk.M[i][j+yd],v,m,l)*SolnBlk.Uo[i][j+yd][index];
	a = vy*(SolnBlk.Uo[i][j+yd].I[index] - SolnBlk.Uo[i-xd][j].I[index]) - dy*S;
	b = (vy-vx*dy/dx )*(SolnBlk.Uo[i][j].I[index] - SolnBlk.Uo[i-xd][j].I[index]);
	Ix_out = SolnBlk.Uo[i][j].I[index];
	Iy_out = SolnBlk.Uo[i-xd][j].I[index] + 0.5*vanalbada(a,b,0.10)/vy;
      }

      // high order correction term
      corr = ( max(a_E, 0.0) + max(a_W, 0.0) ) * ( SolnBlk.U[i][j].I[index] - Ix_out ) +
             ( max(a_N, 0.0) + max(a_S, 0.0) ) * ( SolnBlk.U[i][j].I[index] - Iy_out );
      
      // solve upwind solution with correction term
      wx = ONE;
      wy = ONE;
      break;

    //------------------------------------------------
    // Exponential difference scheme
    //------------------------------------------------
    // See:
    // K.D. Lathrop, "Spatial differencing of the transport 
    // equation: positivity vs. accuracy," in Journal of 
    // computational physics, vol 4, pp. 475-498, 1969.
    case SPACE_MARCH_EXPONENTIAL:
      // cell size
      dx = HALF * (SolnBlk.Grid.lfaceN(i,j) + SolnBlk.Grid.lfaceS(i,j) );
      dy = HALF * (SolnBlk.Grid.lfaceE(i,j) + SolnBlk.Grid.lfaceW(i,j) );
      // optical thickness in each dir
      tau_x = SolnBlk.M[i][j].beta(v) * dx / fabs(U.mu[m][l]/U.omega[m][l]);  // (mu/omega cause weighting worked into mu)
      tau_y = SolnBlk.M[i][j].beta(v) * dy / fabs(U.eta[m][l]/U.omega[m][l]); // (eta/omega cause weighting worked into mu)
      // weights
      wx = ONE / (ONE - exp(-tau_x)) - ONE/tau_x;
      wy = ONE / (ONE - exp(-tau_y)) - ONE/tau_y;
      corr = ZERO;
      break;
      
    //------------------------------------------------
    // Central difference scheme
    //------------------------------------------------
    case SPACE_MARCH_CENTRAL:
      wx = HALF;
      wy = HALF;
      corr = ZERO;
      break;

    //------------------------------------------------
    // Upwind difference scheme
    //------------------------------------------------
    case SPACE_MARCH_UPWIND:
    default:
      wx = ONE;
      wy = ONE;
      corr = ZERO;
      break;

    //------------------------------------------------
    } // end switch
    //------------------------------------------------
      
    
    //
    // compute intensity
    //
    wwx = ( ONE - wx ) / wx;
    wwy = ( ONE - wy ) / wy;
    num = SolnBlk.U[i][j].S(SolnBlk.M[i][j],v,m,l) * V + corr
      + ( max( -a_W, 0.0) + max( a_E, 0.0)*wwx ) * Ix_f[j]
      + ( max( -a_E, 0.0) + max( a_W, 0.0)*wwx ) * Ix_f[j]
      + ( max( -a_S, 0.0) + max( a_N, 0.0)*wwy ) * Iy_f
      + ( max( -a_N, 0.0) + max( a_S, 0.0)*wwy ) * Iy_f;    
    denom = SolnBlk.U[i][j].dSdU(SolnBlk.M[i][j],v,m,l) * V
      + max( a_E, 0.0)/wx
      + max( a_W, 0.0)/wx
      + max( a_N, 0.0)/wy
      + max( a_S, 0.0)/wy;

    //
    // add the axisymmetric source term
    //
    if (SolnBlk.Axisymmetric && 
	Input_Parameters.i_RTE_Solver == RTE2D_SOLVER_FVM) {
      num += SolnBlk.U[i][j].Sa_FVM(v,m,l,SolnBlk.Axisymmetric) * 
	     SolnBlk.Grid.Cell[i][j].A;
      denom += SolnBlk.U[i][j].dSadU_FVM(v,m,l,SolnBlk.Axisymmetric) * 
	       SolnBlk.Grid.Cell[i][j].A;
    } else if (SolnBlk.Axisymmetric && 
	       Input_Parameters.i_RTE_Solver == RTE2D_SOLVER_DOM) {
      num += SolnBlk.U[i][j].Sa_DOM(v,m,l);
      denom += SolnBlk.U[i][j].dSadU_DOM(v,m,l);
    }
    
    SolnBlk.U[i][j].I[index] = num/denom; //max(num/denom, ZERO);


    //
    // store outgoing face value ->  i.e. new incoming face value
    // perform differencing to compute new downstream face values
    //
    switch (i_SpaceMarch_Scheme) {

    //------------------------------------------------
    // CLAM
    //------------------------------------------------
    case SPACE_MARCH_CLAM:
      xf = SolnBlk.Grid.Cell[i][j].Xc.x + 
	xd*QUARTER*(SolnBlk.Grid.lfaceN(i,j)+SolnBlk.Grid.lfaceS(i,j));
      Ix_f[j] = CLAM(SolnBlk.U[i-xd][j].I[index],   // I upstream node
		     SolnBlk.U[i   ][j].I[index],    // I centroid node
		     SolnBlk.U[i+xd][j].I[index],    // I downstream node
		     SolnBlk.Grid.Cell[i-xd][j].Xc.x, // x upstream node
		     SolnBlk.Grid.Cell[i   ][j].Xc.x, // x centroid node
		     SolnBlk.Grid.Cell[i+xd][j].Xc.x, // x downstream node
		     xf);                             // x downstream face
      
      yf = SolnBlk.Grid.Cell[i][j].Xc.y + 
	yd*QUARTER*(SolnBlk.Grid.lfaceE(i,j)+SolnBlk.Grid.lfaceW(i,j));
      Iy_f = CLAM(SolnBlk.U[i][j-yd].I[index],      // I upstream node
		  SolnBlk.U[i][j   ].I[index],    // I centroid node
		  SolnBlk.U[i][j+yd].I[index],    // I downstream node
		  SolnBlk.Grid.Cell[i][j-yd].Xc.y, // y upstream node
		  SolnBlk.Grid.Cell[i][j   ].Xc.y, // y centroid node
		  SolnBlk.Grid.Cell[i][j+yd].Xc.y, // y downstream node
		  yf);                             // y downstream face
      break;

    //------------------------------------------------
    // General Multidimensional
    //------------------------------------------------
    case SPACE_MARCH_GM:
      // vx, vy, dx, dy are already computed
      if ( vx/dx <= vy/dy ) {
	S = SolnBlk.U[i+xd][j].S(SolnBlk.M[i+xd][j],v,m,l) - 
	    SolnBlk.U[i+xd][j].dSdU(SolnBlk.M[i+xd][j],v,m,l)*SolnBlk.U[i+xd][j].I[index];
	a = vx*(SolnBlk.U[i+xd][j].I[index] - SolnBlk.U[i][j-yd].I[index]) - dx*S;
	b = (vx-vy*dx/dy) * (SolnBlk.U[i][j].I[index] - SolnBlk.U[i][j-yd].I[index]);
	Ix_f[j] = SolnBlk.U[i][j-yd].I[index] + 0.5*vanalbada(a,b,0.10)/vx;
	Iy_f = SolnBlk.U[i][j].I[index];

      } else {
	S = SolnBlk.U[i][j+yd].S(SolnBlk.M[i][j+yd],v,m,l) - 
	    SolnBlk.U[i][j+yd].dSdU(SolnBlk.M[i][j+yd],v,m,l)*SolnBlk.U[i][j+yd].I[index];
	a = vy*(SolnBlk.U[i][j+yd].I[index] - SolnBlk.U[i-xd][j].I[index]) - dy*S;
	b = (vy-vx*dy/dx ) * (SolnBlk.U[i][j].I[index] - SolnBlk.U[i-xd][j].I[index]);
	Ix_f[j] = SolnBlk.U[i][j].I[index];
	Iy_f = SolnBlk.U[i-xd][j].I[index] + 0.5*vanalbada(a,b,0.10)/vy;
      }
      break;

    //------------------------------------------------
    // Methods that don't use deffered correction
    //------------------------------------------------
    case SPACE_MARCH_EXPONENTIAL:
    case SPACE_MARCH_CENTRAL:
    case SPACE_MARCH_UPWIND:
    default:
      Ix_f[j] = SolnBlk.U[i][j].I[index]/wx - wwx*Ix_f[j];
      Iy_f = SolnBlk.U[i][j].I[index]/wy - wwy*Iy_f;
      break;

    //------------------------------------------------
    } // end switch
    //------------------------------------------------

    //
    // determine the intensity at the next special direction for DOM
    //
    if (SolnBlk.Axisymmetric && 
	Input_Parameters.i_RTE_Solver == RTE2D_SOLVER_DOM) {
      SolnBlk.U[i][j].I_half[m][l+1] = SolnBlk.U[i][j].I[index]/wpsi - 
	                               wwpsi*SolnBlk.U[i][j].I_half[m][l];
    }

    // Residual successfully evaluated.
    return (0);
    
}

/********************************************************
 * Routine: dUdt_Space_March                            *
 *                                                      *
 * This routine evaluates the residual for the specified*
 * solution block using a space marching procedure.     * 
 * Marching is performed in the direction of the rays.  *
 * The upwind and several 2nd-order and high-res spacial*
 * discretization schemes may be used.                  *
 * The residual is stored in dUdt[][][0].               *
 *                                                      *
 *                                                      *
 *                                                      *
 ********************************************************/
int dUdt_Space_March(Rte2D_Quad_Block &SolnBlk,
		     Rte2D_Input_Parameters &Input_Parameters)
{
  //
  // declares
  //
  int i, j;                               // index counters
  Rte2D_State &U = SolnBlk.U[0][0];       // alias
  double  Iy_f;                           // y face intensity
  double *Ix_f = new double[SolnBlk.NCj]; // x face intensity
  int varindex;                           // variable index
  int xd, yd;                             // i+xd -> downstream x index
                                          // i+yd -> downstream y index

  //
  // For axisymmetric dom, compute ART terms
  //
  if (Input_Parameters.i_RTE_Solver==RTE2D_SOLVER_DOM && 
      Input_Parameters.Axisymmetric) 
    for ( i = SolnBlk.ICl ; i <= SolnBlk.ICu ; ++i )
      for ( j = SolnBlk.JCl ; j <= SolnBlk.JCu ; ++j ) 
	SolnBlk.U[i][j].SetupART_DOM( SolnBlk.Grid.nfaceE(i,j), 
				      SolnBlk.Grid.lfaceE(i,j)*SolnBlk.SpE[i][j],
				      SolnBlk.Grid.nfaceW(i,j), 
				      SolnBlk.Grid.lfaceW(i,j)*SolnBlk.SpW[i][j],
				      SolnBlk.Grid.nfaceN(i,j), 
				      SolnBlk.Grid.lfaceN(i,j)*SolnBlk.SpN[i][j],
				      SolnBlk.Grid.nfaceS(i,j), 
				      SolnBlk.Grid.lfaceS(i,j)*SolnBlk.SpS[i][j] );

  // store the old solution
  for ( i = SolnBlk.ICl ; i <= SolnBlk.ICu ; ++i ) 	      
    for ( j = SolnBlk.JCl ; j <= SolnBlk.JCu ; ++j )    
      SolnBlk.Uo[i][j] = SolnBlk.U[i][j];

  // Evaluate the time rate of change of the solution
  // (i.e., the solution residuals) using a second-order
  // limited upwind scheme with a variety of flux functions.
   	
  //
  // Loop over each band, direction
  // Note: Be careful, you have to sweep the directions 
  // in such a manner that the r-direction cosing is increasing 
  // monotonically.
  // 
  for (int v=0; v<SolnBlk.U[0][0].Nband; v++) 
    for(int m=0 ; m<SolnBlk.U[0][0].Npolar ; m++) 
      for(int l=0 ; l<SolnBlk.U[0][0].Nazim[m] ; l++) {
	  
	// the index
	varindex = U.Index[v][m][l]+1;
	  
	//------------------------------------------------
	// mu<0, eta<0
	//------------------------------------------------
	if ( U.mu[m][l]<ZERO && U.eta[m][l]<ZERO ) {
	    
	  // downstream index directions
	  xd = -1;  yd = -1;

	  // set EAST BCs
	  for ( j = SolnBlk.JCl ; j <= SolnBlk.JCu ; ++j ) 
	    Ix_f[j] = SolnBlk.UoE[j][varindex];

	  //
	  // Loop over the grid
	  //
	  for ( i = SolnBlk.ICu ; i >= SolnBlk.ICl ; --i ) {
	      	      
	    // set NORTH BCs
	    Iy_f = SolnBlk.UoN[i][varindex];
	      
	    for ( j = SolnBlk.JCu ; j >= SolnBlk.JCl ; --j ) {    
	      dUdt_Space_March_Flux_Eval(SolnBlk,Input_Parameters,
					 i, j, v, m, l, xd, yd, Ix_f, Iy_f);
	    } // endfor - j

	    // save SOUTH BCs
	    SolnBlk.UoS[i][varindex] = Iy_f;

	  } // endfor - i

	  // save WEST BCs
	  for ( j = SolnBlk.JCl ; j <= SolnBlk.JCu ; ++j )      
	    SolnBlk.UoW[j][varindex] = Ix_f[j];
    

	//------------------------------------------------
	// mu<0, eta>0 
	//------------------------------------------------
	} else if ( U.mu[m][l]<ZERO && U.eta[m][l]>=ZERO ) {

	  // downstream index directions
	  xd = -1;  yd = 1;

	  // set EAST BCs
	  for ( j = SolnBlk.JCl ; j <= SolnBlk.JCu ; ++j )      
	    Ix_f[j] = SolnBlk.UoE[j][varindex];

	  //
	  // Loop over the grid
	  //
	  for ( i = SolnBlk.ICu ; i >= SolnBlk.ICl ; --i ) {
	      
	    // set SOUTH BCs
	    Iy_f = SolnBlk.UoS[i][varindex];

	    for ( j = SolnBlk.JCl ; j <= SolnBlk.JCu ; ++j ) {  
	      dUdt_Space_March_Flux_Eval(SolnBlk,Input_Parameters,
					 i, j, v, m, l, xd, yd, Ix_f, Iy_f);
	    } // endfor - j

	    // save NORTH BCs
	    SolnBlk.UoN[i][varindex] = Iy_f;
	      
	  } // endfor - i 

	  // save WEST BCs
	  for ( j = SolnBlk.JCl ; j <= SolnBlk.JCu ; ++j )      
	    SolnBlk.UoW[j][varindex] = Ix_f[j];


	//------------------------------------------------
	// mu>0, eta<0 
	//------------------------------------------------
	} else if ( U.mu[m][l]>=ZERO && U.eta[m][l]<ZERO ) {

	  // downstream index directions
	  xd = 1;  yd = -1;

	  // set WEST BCs
	  for ( j = SolnBlk.JCl ; j <= SolnBlk.JCu ; ++j )      
	    Ix_f[j] = SolnBlk.UoW[j][varindex];
	    
	  //
	  // Loop over the grid
	  //
	  for ( i = SolnBlk.ICl ; i <= SolnBlk.ICu ; ++i ) {
	      
	    // set NORTH BCs
	    Iy_f = SolnBlk.UoN[i][varindex];

	    for ( j = SolnBlk.JCu ; j >= SolnBlk.JCl ; --j ) {    

	      dUdt_Space_March_Flux_Eval(SolnBlk, Input_Parameters,
					 i, j, v, m, l, xd, yd, Ix_f, Iy_f);
   
	    } // endfor - j

	    // save SOUTH BCs
	    SolnBlk.UoS[i][varindex] = Iy_f;
	      	    
	  } // endfor - i

	  // save EAST BCs
	  for ( j = SolnBlk.JCl ; j <= SolnBlk.JCu ; ++j )      
	    SolnBlk.UoE[j][varindex] = Ix_f[j];

	//------------------------------------------------
	// mu>0, eta>0
	//------------------------------------------------
	} else if ( U.mu[m][l]>=ZERO && U.eta[m][l]>=ZERO ) {

	  // downstream index directions
	  xd = 1;  yd = 1;

	  // set WEST BCs
	  for ( j = SolnBlk.JCl ; j <= SolnBlk.JCu ; ++j )      
	    Ix_f[j] = SolnBlk.UoW[j][varindex];
	    
	  // Loop over the grid
	  for ( i = SolnBlk.ICl ; i <= SolnBlk.ICu ; ++i ) {

	    // set SOUTH BCs
	    Iy_f = SolnBlk.UoS[i][varindex];

	    for ( j = SolnBlk.JCl ; j <= SolnBlk.JCu ; ++j ) {      
	      dUdt_Space_March_Flux_Eval(SolnBlk,Input_Parameters,
					 i, j, v, m, l, xd, yd, Ix_f, Iy_f);
	    } // endfor - j

	    // save NORTH BCs
	    SolnBlk.UoN[i][varindex] = Iy_f;
	    SolnBlk.FluxN[i][varindex] = SolnBlk.UoN[i].Fn(SolnBlk.Grid.nfaceN(i,SolnBlk.JCu),v,m,l) * 
	      SolnBlk.Grid.lfaceN(i, SolnBlk.JCu)*SolnBlk.SpN[i][SolnBlk.JCu];
	      
	  } // endfor - i 

	  // save EAST BCs
	  for ( j = SolnBlk.JCl ; j <= SolnBlk.JCu ; ++j )      
	    SolnBlk.UoE[j][varindex] = Ix_f[j];

	//------------------------------------------------
	} // endif - mu, eta>0
	//------------------------------------------------
  
      } // endfor - bands, directions
	
  // update solution change
  for ( i = SolnBlk.ICl ; i <= SolnBlk.ICu ; ++i ) 	      
    for ( j = SolnBlk.JCl ; j <= SolnBlk.JCu ; ++j )    
      SolnBlk.dUdt[i][j][0] = SolnBlk.U[i][j] - SolnBlk.Uo[i][j];

  // store the boundary fluxes
  for ( i = SolnBlk.ICl ; i <= SolnBlk.ICu ; ++i ) {
    SolnBlk.FluxS[i] = SolnBlk.UoS[i].Fn(SolnBlk.Grid.nfaceS(i,j)) * 
      SolnBlk.Grid.lfaceS(i, j)*SolnBlk.SpS[i][j];
    SolnBlk.FluxN[i] = SolnBlk.UoN[i].Fn(SolnBlk.Grid.nfaceN(i,j)) * 
      SolnBlk.Grid.lfaceN(i, j)*SolnBlk.SpN[i][j];
  }     
  for ( j = SolnBlk.JCl ; j <= SolnBlk.JCu ; ++j ) {
    SolnBlk.FluxE[j] = SolnBlk.UoE[j].Fn(SolnBlk.Grid.nfaceE(i,j)) * 
      SolnBlk.Grid.lfaceE(i, j)*SolnBlk.SpE[i][j];
    SolnBlk.FluxW[j] = SolnBlk.UoW[j].Fn(SolnBlk.Grid.nfaceW(i,j)) * 
      SolnBlk.Grid.lfaceW(i, j)*SolnBlk.SpW[i][j];
  }
    
  /* Residual successfully evaluated. */
  delete[] Ix_f; Ix_f = NULL;
    
  return (0);

}

/********************************************************
 * Routine: Residual_Jacobian                           *
 *                                                      *
 * This routine returns residual Jacobian matrix for the*
 * specified local solution block.                      *
 *                                                      *
 ********************************************************/
DenseMatrix Residual_Jacobian(const int i, 
			      const int j, 
			      Rte2D_Quad_Block &SolnBlk, 
			      const int cell_flag) {
  
  int NUM_VAR_RTE2D = SolnBlk.NumVar(); 

  DenseMatrix dRdU(NUM_VAR_RTE2D, NUM_VAR_RTE2D),
              dFdU(NUM_VAR_RTE2D, NUM_VAR_RTE2D);

  Vector2D nface_N, nface_S, nface_E, nface_W;

  if (i < SolnBlk.ICl || i > SolnBlk.ICu ||
      j < SolnBlk.JCl || j > SolnBlk.JCu) {
    // GHOST CELL
    dRdU.zero();
    return dRdU;
    
  } else {
    // NON-GHOST CELL.
    
    dRdU.zero();
    
    nface_N = SolnBlk.Grid.nfaceN(i  , j-1);
    nface_S = SolnBlk.Grid.nfaceS(i  , j+1);
    nface_E = SolnBlk.Grid.nfaceE(i-1, j);
    nface_W = SolnBlk.Grid.nfaceW(i+1, j);

    cout << "\n\n***Residual_Jacobian in Rte2DQuadSingleBlock "
	 << "not implemented yet!!!!!\n\n";
    
    dFdU.zero();
    dFndU(dFdU, SolnBlk.Uo[i][j], nface_N);
    dRdU = (SolnBlk.Grid.lfaceN(i,j-1)*SolnBlk.SpN[i][j-1]/(SolnBlk.Grid.area(i, j)*SolnBlk.Sp[i][j]))*dFdU;

    dFdU.zero();
    dFndU(dFdU, SolnBlk.Uo[i][j], nface_S);
    dRdU += (SolnBlk.Grid.lfaceS(i,j+1)*SolnBlk.SpS[i][j+1]/(SolnBlk.Grid.area(i, j)*SolnBlk.Sp[i][j]))*dFdU;

    dFdU.zero();
    dFndU(dFdU, SolnBlk.Uo[i][j], nface_E);
    dRdU += (SolnBlk.Grid.lfaceE(i-1,j)*SolnBlk.SpE[i-1][j]/(SolnBlk.Grid.area(i, j)*SolnBlk.Sp[i][j]))*dFdU;

    dFdU.zero();
    dFndU(dFdU, SolnBlk.Uo[i][j], nface_W);
    dRdU += (SolnBlk.Grid.lfaceW(i+1,j)*SolnBlk.SpW[i+1][j]/(SolnBlk.Grid.area(i, j)*SolnBlk.Sp[i][j]))*dFdU;


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
int Update_Solution_Multistage_Explicit(Rte2D_Quad_Block &SolnBlk,
                                        const int i_stage,
                                        Rte2D_Input_Parameters &Input_Parameters) {

    int NUM_VAR_RTE2D = SolnBlk.NumVar(); 
    int i, j, n, k, l, k_residual, n_residual_reduction;
    double omega, residual_reduction_factor;

    // Memory for linear system solver.
    DenseMatrix dRdU(NUM_VAR_RTE2D,NUM_VAR_RTE2D),
                P_inv(NUM_VAR_RTE2D,NUM_VAR_RTE2D);
    DenseSystemLinEqs LinSys;
    Rte2D_State dU_precon;

    /* Allocate memory for linear system solver. */

    LinSys.allocate(NUM_VAR_RTE2D);

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
    
    /* Perform update of solution variables for stage i_stage of an N stage scheme. */
    /* Update solution variables for this stage. */
    
    for ( j  = SolnBlk.JCl ; j <= SolnBlk.JCu ; ++j ) {
       for ( i = SolnBlk.ICl ; i <= SolnBlk.ICu ; ++i ) {

	 /******************************************************/
	 /********** Fully Explicit ****************************/
	 /******************************************************/
	 if (Input_Parameters.Local_Time_Stepping == GLOBAL_TIME_STEPPING || 
             Input_Parameters.Local_Time_Stepping == SCALAR_LOCAL_TIME_STEPPING) {
           SolnBlk.U[i][j] = SolnBlk.Uo[i][j] + 
                             omega*SolnBlk.dUdt[i][j][k_residual];
//          } else if (Input_Parameters.Local_Time_Stepping == LOW_MACH_NUMBER_WEISS_SMITH_PRECONDITIONER) {
// 	    // Apply Weiss-Smith low-Mach-number preconditioning to the residual.
//             dU_precon.Zero();
//             SolnBlk.W[i][j].P_U_WS_inv(P_inv);
//             for ( k = 1 ; k <= NUM_VAR_EULER2D ; ++k ) {  
//                for ( l = 1 ; l <= NUM_VAR_EULER2D ; ++l ) { 
//                   dU_precon[k] += P_inv(k-1,l-1)*omega*SolnBlk.dUdt[i][j][k_residual][l];
//                } /* endfor */
//             } /* endfor */
//             SolnBlk.U[i][j] = SolnBlk.Uo[i][j] + dU_precon;
	 } /* endif */

	 //Check for unphysical properties  
	 /**********************************************************/
	 /* If unphysical properties and using global timestepping */ 
	 /* stop simulation                                        */
	 /**********************************************************/   
	 if (Input_Parameters.Local_Time_Stepping == GLOBAL_TIME_STEPPING && 
	     SolnBlk.U[i][j].Unphysical_Properties() ) {
	     cout << "\n " << CFFC_Name() << " Rte2D ERROR: Negative Value: \n"
		  << " cell = (" << i << ", " << j << ") " 
		  << " X = " << SolnBlk.Grid.Cell[i][j].Xc << "\n"
		  << " Uo = " << SolnBlk.Uo[i][j] << "\n"
		  << " U = " << SolnBlk.U[i][j] << "\n"
		  << " dUdt = " << SolnBlk.dUdt[i][j][k_residual] << "\n";
	     return (i);

 	/*********************************************************/
	/* If unphysical properties and using local timestepping */ 
	/* try reducing step size                                */
	/*********************************************************/    
        } else if ((Input_Parameters.Local_Time_Stepping == SCALAR_LOCAL_TIME_STEPPING || 
		    Input_Parameters.Local_Time_Stepping == LOW_MACH_NUMBER_WEISS_SMITH_PRECONDITIONER) && 
		   SolnBlk.U[i][j].Unphysical_Properties()) {
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

	     if (!SolnBlk.U[i][j].Unphysical_Properties()) {
		 break;
	     } /* endif */
	   } /* end for */

	   if (SolnBlk.U[i][j].Unphysical_Properties()) {
	     cout << "\n " << CFFC_Name() << " Rte2D ERROR: Negative Intensity: \n"
		  << " cell = (" << i << ", " << j << ") " 
		  << " X = " << SolnBlk.Grid.Cell[i][j].Xc << "\n"
		  << " Uo = " << SolnBlk.Uo[i][j] << "\n"
		  << " U = " << SolnBlk.U[i][j] << "\n"
		  << " dUdt = " << SolnBlk.dUdt[i][j][k_residual] << "\n";
	     return (i);
	   } /* endif */
	 } /* end if */
	 
	 /**************************************************************/
	 /************ SEMI-IMPLICIT AND/OR PRECONDITIONING ************/
	 /**************************************************************/
	 if (Input_Parameters.Local_Time_Stepping == MATRIX_LOCAL_TIME_STEPPING) {
	   
	   /* Evaluate Block Point-Jacobi Preconditioner dFdU */
	   
	   dRdU = Residual_Jacobian(i,
				    j,
				    SolnBlk,
				    0);
	   
	   /* Carry out matrix time stepping approach by setting up
	      linear system of equation for given cell */
	   
	   LinSys.A = (-1.00*SolnBlk.dt[i][j])*dRdU;
	   
	   for (k = 1; k <= NUM_VAR_RTE2D; ++k) {
	     LinSys.b(k-1) = omega*SolnBlk.dUdt[i][j][k_residual][k];
	   } /* endfor */
	   
	   /* Solve system of equations using LU decomposition
	      Gaussian elimination procedure. */
	   
	   LinSys.solve(LU_DECOMPOSITION);
	   
	   /* Update the conserved solution variables. */
	   
	   for (k = 1; k <= NUM_VAR_RTE2D; ++k) {
	     SolnBlk.U[i][j][k] = SolnBlk.Uo[i][j][k] + LinSys.x(k-1);
	   } /* endfor */
	   
	   if (SolnBlk.U[i][j].Unphysical_Properties()) {	     
	     residual_reduction_factor = ONE;
	     
	     for (n_residual_reduction = 1; n_residual_reduction <= 10; ++n_residual_reduction) {
	       residual_reduction_factor = HALF*residual_reduction_factor;
               if (n_residual_reduction = 10) residual_reduction_factor = ZERO;
	       SolnBlk.dt[i][j] = residual_reduction_factor*SolnBlk.dt[i][j];
	       SolnBlk.dUdt[i][j][k_residual] = residual_reduction_factor*SolnBlk.dUdt[i][j][k_residual];
	       LinSys.A = (-1.00*SolnBlk.dt[i][j])*dRdU;
	       for (k = 1; k <= NUM_VAR_RTE2D; ++k) {
		 LinSys.b(k-1) = omega*SolnBlk.dUdt[i][j][k_residual][k];
	       } /* endfor */
	       LinSys.solve(LU_DECOMPOSITION);
	       for (k = 1; k <= NUM_VAR_RTE2D; ++k) {
		 SolnBlk.U[i][j][k] = SolnBlk.Uo[i][j][k] + LinSys.x(k-1);
	       } /* endfor */
	       if (!SolnBlk.U[i][j].Unphysical_Properties()) break;
	     } /* endfor */
	   } /* endif */
	   
	   if (SolnBlk.U[i][j].Unphysical_Properties()) {	     
	     cout << "\n " << CFFC_Name() << " Rte2D ERROR: Negative Intensity: \n"
		  << " cell = (" << i << ", " << j << ") " 
		  << " X = " << SolnBlk.Grid.Cell[i][j].Xc << "\n U = " 
		  << " Uo = " << SolnBlk.Uo[i][j] << "\n"
		  << " U = " << SolnBlk.U[i][j] << "\n"
		  << " dUdt = " << SolnBlk.dUdt[i][j][k_residual] << "\n";
	     return (i);
	   } /* endif */

	 } /* endif */   

         	 
       } /* endfor */    
    } /* endfor */
    
    /* Deallocate memory for linear system solver. */

    LinSys.deallocate();

    /* Solution successfully updated. */
    
    return (0);   
}




/**********************************************************************
 * Routine: Output_Exact                                              *
 *                                                                    *
 * Writes the exact and computed black enclosure rad  solution values *
 * at the nodes of the specified quadrilateral solution block to the  *
 * specified output stream suitable for plotting with TECPLOT.  The   *
 * error norms are also computed.                                     *
 *                                                                    *
 **********************************************************************/
void Output_Exact(Rte2D_Quad_Block &SolnBlk,
		  const int Block_Number,
		  const int Output_Title,
		  ostream &Out_File,
		  double &l1_norm,
		  double &l2_norm,
		  double &max_norm,
		  Rte2D_Input_Parameters &IP) {

  // decalres
  double G, G_e, qx, qx_e, qy, qy_e;
  double tau, c;
  double rpos, zpos;
  double offset;

  // determine some case parameters 
  if (IP.Axisymmetric) {
    c = (IP.Pipe_Length/TWO)/IP.Pipe_Radius;
    tau = IP.AbsorptionCoef*IP.Pipe_Radius;
    offset = IP.Pipe_Length / TWO;
  }

  // Output node solution data header.  
  Out_File << setprecision(14);
  if (Output_Title)
    Out_File << "TITLE = \"" << CFFC_Name() << ": 2D Rte2D Black Enclosure Solution "
	     << "\"" << "\n"
	     << "VARIABLES = \"x\" \\ \n"
	     << "\"y\" \\ \n"
	     << "\"G\" \\ \n"
	     << "\"q.x\" \\ \n"
	     << "\"q.y\" \\ \n"
	     << "\"G_exact\" \\ \n"
	     << "\"q.x_exact\" \\ \n"
	     << "\"q.y_exact\" \\ \n"
	     << "\"G_error\" \\ \n"
	     << "\"q.x_error\" \\ \n"
	     << "\"q.y_error\" \\ \n";
  Out_File << "ZONE T =  \"Block Number = " << Block_Number
	   << "\" \\ \n"
	   << "I = " << SolnBlk.Grid.NCi << " \\ \n"
	   << "J = " << SolnBlk.Grid.NCj << " \\ \n"
	   << "F = POINT \\ \n";

  Out_File.setf(ios::scientific);
 
  //
  // loop through the grid
  //
  for (int j = 0; j < SolnBlk.Grid.NCj; j++) 
    for (int i = 0; i < SolnBlk.Grid.NCi; i++) {

      // values normalized by blackbody intensity
      G = SolnBlk.U[i][j].G(SolnBlk.M[i][j]) / ( FOUR*PI*BlackBody(IP.Temperature) );
      qx = SolnBlk.U[i][j].q(SolnBlk.M[i][j]).x / ( PI*BlackBody(IP.Temperature) );
      qy = SolnBlk.U[i][j].q(SolnBlk.M[i][j]).y / ( PI*BlackBody(IP.Temperature) );

      if (i < SolnBlk.Nghost || i > SolnBlk.ICu ||
	  j < SolnBlk.Nghost || j > SolnBlk.JCu) {
	G_e = G;
	qx_e = qx;
	qy_e = qy;

      } else {

	//
	// compute the exact solution depending upon flow type
	//
	switch (IP.Axisymmetric) {

	//------------------------------------------------
	// PLANAR RADIATION
	//------------------------------------------------
	case PLANAR:
	  RectangularEnclosure( IP.Temperature,
				IP.AbsorptionCoef,
				-IP.Box_Width/TWO,
				 IP.Box_Width/TWO,
				-IP.Box_Height/TWO,
				 IP.Box_Height/TWO,
				SolnBlk.Grid.Cell[i][j].Xc.x,
				SolnBlk.Grid.Cell[i][j].Xc.y,
				G_e, 
				qx_e, 
				qy_e );
	  break;
	  
	//------------------------------------------------
	// AXISYMMETRIC RADIATION (r->x)
	//------------------------------------------------
	case AXISYMMETRIC_X:
	  
	  rpos = SolnBlk.Grid.Cell[i][j].Xc.x / IP.Pipe_Radius;
	  zpos = (SolnBlk.Grid.Cell[i][j].Xc.y-offset) / IP.Pipe_Radius;
	  CylindricalEnclosure( IP.Temperature,
				c,
				tau,
				rpos,
				zpos,
				G_e, 
				qx_e, 
				qy_e );
	  break;

	//------------------------------------------------
	// AXISYMMETRIC RADIATION (r->y)
	//------------------------------------------------
	case AXISYMMETRIC_Y:
	  rpos = SolnBlk.Grid.Cell[i][j].Xc.y / IP.Pipe_Radius;
	  zpos = (SolnBlk.Grid.Cell[i][j].Xc.x-offset) / IP.Pipe_Radius;
	  CylindricalEnclosure( IP.Temperature,
				c,
				tau,
				rpos,
				zpos,
				G_e, 
				qy_e, 
				qx_e );
	  break;

	//------------------------------------------------
	} // endswitch Axisymmetric
	//------------------------------------------------

	// compute the norms
	l1_norm += fabs(G - G_e);
	l2_norm += sqr(G - G_e);
	max_norm = max(max_norm,fabs(G - G_e));

      } // endif

      // Output node solution data
      Out_File << SolnBlk.Grid.Cell[i][j].Xc << " " 
	       << G << " " 
	       << qx << " " 
	       << qy << " " 
	       << G_e << " " 
	       << qx_e << " " 
	       << qy_e << " " 
	       << (G - G_e) << " " 
	       << (qx - qx_e) << " "
	       << (qy - qy_e) << endl;
      
    } // endfor - grid
  
  // the norm
  l2_norm = sqrt(l2_norm);

}
