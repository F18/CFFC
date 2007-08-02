/************** Chem2DQuadSingleBlock.cc ***************************
   Single-Block Versions of Subroutines for 2D Euler 
   Multi-Block Quadrilateral Mesh  Solution Classes. 

   NOTES:
          
   - based on Euler2DQuadSingleBlock.cc
**********************************************************************/

/* Include 2D Euler quadrilateral mesh solution header file. */
#ifndef _CHEM2D_dRdU_INCLUDED
#include "dRdU.h"
#endif // _CHEM2D_dRdU_INCLUDE

/*************************************************************************
* Chem2D_Quad_Block -- Single Block External Subroutines.                *
**************************************************************************/

/********************************************************
 * Routine: Write_Solution_Block                        *
 *                                                      *
 * Writes the cell centred solution values of the       *
 * specified quadrilateral solution block to the        *
 * specified output stream for restart purposes.        *
 *                                                      *
 ********************************************************/
void Write_Solution_Block(Chem2D_Quad_Block &SolnBlk,
	                  ostream &Out_File) {

   Out_File << setprecision(14) << SolnBlk << setprecision(6);

}

/********************************************************
 * Routine: Read_Solution_Block                         *
 *                                                      *
 * Reads the cell centred solution values for the       *
 * specified quadrilateral solution block from the      *
 * specified output stream as required for restart      *
 * purposes.                                            *
 *                                                      *
 ********************************************************/
void Read_Solution_Block(Chem2D_Quad_Block &SolnBlk,
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
void Broadcast_Solution_Block(Chem2D_Quad_Block &SolnBlk) {

#ifdef _MPI_VERSION

    int ni, nj, ng, nr, block_allocated, buffer_size;
    double *buffer;

    int NUM_VAR_CHEM2D = SolnBlk.NumVar(); 

    /* Broadcast the number of cells in each direction. */
    if (CFDkit_Primary_MPI_Processor()) {
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
    MPI::COMM_WORLD.Bcast(&nr, 1, MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(&block_allocated, 1, MPI::INT, 0);

    /* On non-primary MPI processors, allocate (re-allocate) 
       memory for the quadrilateral solution block as necessary. */

    if (!CFDkit_Primary_MPI_Processor()) {
      if (SolnBlk.NCi != ni || SolnBlk.NCj != nj || SolnBlk.Nghost != ng ) { 
	if (SolnBlk.U != NULL) SolnBlk.deallocate();
	if (block_allocated) SolnBlk.allocate(ni-2*ng, nj-2*ng, ng); 
      } /* endif */
      // Set the block static variables if they were not previously assigned.
      if (SolnBlk.residual_variable != nr) SolnBlk.residual_variable = nr;
    } /* endif */

    /* Broadcast the axisymmetric/planar flow, viscous, turbulent, and gravity indicators. */
   
    MPI::COMM_WORLD.Bcast(&(SolnBlk.Flow_Type), 1, MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(&(SolnBlk.Axisymmetric), 1, MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(&(SolnBlk.Gravity), 1, MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(&(SolnBlk.debug_level), 1, MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(&(SolnBlk.Moving_wall_velocity), 1, MPI::DOUBLE, 0);

    /* Broadcast the grid. */

    Broadcast_Quad_Block(SolnBlk.Grid);

    /* Broadcast the solution state variables. */

    if (block_allocated) {
      ni = (SolnBlk.ICu+SolnBlk.Nghost) - (SolnBlk.ICl-SolnBlk.Nghost) + 1;
      nj = (SolnBlk.JCu+SolnBlk.Nghost) - (SolnBlk.JCl-SolnBlk.Nghost) + 1; 
      buffer = new double[NUM_VAR_CHEM2D*ni*nj];
    
      if (CFDkit_Primary_MPI_Processor()) {
	buffer_size = 0;
	for (int j  = SolnBlk.JCl-SolnBlk.Nghost ; j <= SolnBlk.JCu+SolnBlk.Nghost ; ++j ) {
	  for (int i = SolnBlk.ICl-SolnBlk.Nghost ; i <= SolnBlk.ICu+SolnBlk.Nghost ; ++i ) {
	    //Changed for Chem2D
	    for ( int k = 0; k < NUM_VAR_CHEM2D; ++k) {	
	      buffer[buffer_size] = SolnBlk.U[i][j][k+1];
	      buffer_size = buffer_size + 1;
	    }	      
	  } /* endfor */
	} /* endfor */
      } /* endif */
    
       buffer_size = NUM_VAR_CHEM2D*ni*nj;
       MPI::COMM_WORLD.Bcast(buffer, buffer_size, MPI::DOUBLE, 0);

       if (!CFDkit_Primary_MPI_Processor()) {
	 buffer_size = 0;
	 for (int j  = SolnBlk.JCl-SolnBlk.Nghost ; j <= SolnBlk.JCu+SolnBlk.Nghost ; ++j ) {
	   for (int i = SolnBlk.ICl-SolnBlk.Nghost ; i <= SolnBlk.ICu+SolnBlk.Nghost ; ++i ) {
	     for ( int k = 0 ; k < NUM_VAR_CHEM2D; ++ k) {	
	       SolnBlk.U[i][j][k+1]= buffer[buffer_size];
	       buffer_size = buffer_size + 1;
	     }
	     SolnBlk.W[i][j] = W(SolnBlk.U[i][j]);
	   } /* endfor */
	 } /* endfor */
       } /* endif */
       
       delete []buffer; 
       buffer = NULL;

       ni = 1;
       nj = (SolnBlk.JCu+SolnBlk.Nghost) - (SolnBlk.JCl-SolnBlk.Nghost) + 1;
       buffer = new double[2*NUM_VAR_CHEM2D*ni*nj];

       if (CFDkit_Primary_MPI_Processor()) {
	 buffer_size = 0;
         for (int j  = SolnBlk.JCl-SolnBlk.Nghost ; j <= SolnBlk.JCu+SolnBlk.Nghost ; ++j ) {
	   for ( int k = 0 ; k < NUM_VAR_CHEM2D; ++ k) {	
	     buffer[buffer_size]= SolnBlk.WoW[j][k+1];
	     buffer_size = buffer_size + 1;
	   }
	   for ( int k = 0 ; k < NUM_VAR_CHEM2D; ++ k) {
	     buffer[buffer_size]= SolnBlk.WoE[j][k+1];  
	     buffer_size = buffer_size + 1;
	   }
	 } /* endfor */
       } /* endif */

       buffer_size = 2*NUM_VAR_CHEM2D*ni*nj;
       MPI::COMM_WORLD.Bcast(buffer, buffer_size, MPI::DOUBLE, 0);

       if (!CFDkit_Primary_MPI_Processor()) {
	 buffer_size = 0;
	 for (int j  = SolnBlk.JCl-SolnBlk.Nghost ; j <= SolnBlk.JCu+SolnBlk.Nghost ; ++j ) {    
	   for ( int k = 0 ; k < NUM_VAR_CHEM2D; ++ k) {	
	     SolnBlk.WoW[j][k+1] = buffer[buffer_size];
	     buffer_size = buffer_size + 1;
	   }
	   for ( int k = 0 ; k < NUM_VAR_CHEM2D; ++ k) {
	     SolnBlk.WoE[j][k+1] = buffer[buffer_size]; 
	     buffer_size = buffer_size + 1;
	   }
          } /* endfor */
       } /* endif */

       delete []buffer; 
       buffer = NULL;

       ni = (SolnBlk.ICu+SolnBlk.Nghost) - (SolnBlk.ICl-SolnBlk.Nghost) + 1;
       nj = 1;
       buffer = new double[2*NUM_VAR_CHEM2D*ni*nj];

       if (CFDkit_Primary_MPI_Processor()) {
	 buffer_size = 0;
	 for (int i  = SolnBlk.JCl-SolnBlk.Nghost ; i <= SolnBlk.JCu+SolnBlk.Nghost ; ++i ) {   
	   for ( int k = 0 ; k < NUM_VAR_CHEM2D; ++ k) {	
	     buffer[buffer_size]= SolnBlk.WoS[i][k+1];	  
	     buffer_size = buffer_size + 1;
	   }
	   for ( int k = 0 ; k < NUM_VAR_CHEM2D; ++ k) {
	     buffer[buffer_size]= SolnBlk.WoN[i][k+1];
	     buffer_size = buffer_size + 1;
	   }
	 } /* endfor */
       } /* endif */

       buffer_size = 2*NUM_VAR_CHEM2D*ni*nj;
       MPI::COMM_WORLD.Bcast(buffer, buffer_size, MPI::DOUBLE, 0);

       if (!CFDkit_Primary_MPI_Processor()) {
	 buffer_size = 0;
	 for (int i = SolnBlk.ICl-SolnBlk.Nghost ; i <= SolnBlk.ICu+SolnBlk.Nghost ; ++i ) {     
	   for ( int k = 0 ; k < NUM_VAR_CHEM2D; ++ k) {	
	     SolnBlk.WoS[i][k+1] = buffer[buffer_size]; 
	     buffer_size = buffer_size + 1;
	    }
	   for ( int k = 0 ; k < NUM_VAR_CHEM2D; ++ k) {
	     SolnBlk.WoN[i][k+1] = buffer[buffer_size];
	     buffer_size = buffer_size + 1;
	   }

	 } /* endfor */
       } /* endif */
       
       delete []buffer; 
       buffer = NULL;

    } /* endif */
#endif

}

#ifdef _MPI_VERSION
//used in AMR
/********************************************************
* Routine: Broadcast_Solution_Block                    *
*                                                      *
* Broadcast quadrilateral solution block to all        *
* processors associated with the specified communicator*
* from the specified processor using the MPI broadcast *
* routine.                                             *
*                                                      *
********************************************************/
void Broadcast_Solution_Block(Chem2D_Quad_Block &SolnBlk,
                             MPI::Intracomm &Communicator, 
                             const int Source_CPU) {

  int Source_Rank = 0;
  int ni, nj, ng, nr, block_allocated, buffer_size;
  double *buffer;

  int NUM_VAR_CHEM2D = SolnBlk.NumVar();

  /* Broadcast the number of cells in each direction. */

  if (CFDkit_MPI::This_Processor_Number == Source_CPU) {
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

  if (!(CFDkit_MPI::This_Processor_Number == Source_CPU)) {
    if (SolnBlk.NCi != ni || SolnBlk.NCj != nj || SolnBlk.Nghost != ng) { 
      if (SolnBlk.U != NULL) SolnBlk.deallocate();
      if (block_allocated) SolnBlk.allocate(ni-2*ng, nj-2*ng, ng);
    } /* endif */
    // Set the block static variables if they were not previously assigned.
    if (SolnBlk.residual_variable != nr) SolnBlk.residual_variable = nr;
  } /* endif */

    /* Broadcast the axisymmetric/planar flow indicator. */

    Communicator.Bcast(&(SolnBlk.Flow_Type), 1, MPI::INT, Source_Rank);
    Communicator.Bcast(&(SolnBlk.Axisymmetric), 1, MPI::INT, Source_Rank);
    Communicator.Bcast(&(SolnBlk.Gravity), 1, MPI::INT, Source_Rank);
    Communicator.Bcast(&(SolnBlk.debug_level), 1, MPI::INT, Source_Rank);
    Communicator.Bcast(&(SolnBlk.Moving_wall_velocity), 1, MPI::DOUBLE, Source_Rank);

    /* Broadcast the grid. */

    Broadcast_Quad_Block(SolnBlk.Grid, Communicator, Source_CPU);
    
    /* Broadcast the solution state variables. */

    if (block_allocated) {
       ni = (SolnBlk.ICu+SolnBlk.Nghost) - (SolnBlk.ICl-SolnBlk.Nghost) + 1;
       nj = (SolnBlk.JCu+SolnBlk.Nghost) - (SolnBlk.JCl-SolnBlk.Nghost) + 1;
       buffer = new double[NUM_VAR_CHEM2D*ni*nj];

       if (CFDkit_MPI::This_Processor_Number == Source_CPU) {
	 buffer_size = 0;
	 for (int j  = SolnBlk.JCl-SolnBlk.Nghost ; j <= SolnBlk.JCu+SolnBlk.Nghost ; ++j ) {
	   for (int i = SolnBlk.ICl-SolnBlk.Nghost ; i <= SolnBlk.ICu+SolnBlk.Nghost ; ++i ) {
	     //Changed for Chem2D
	     for ( int k=0; k<NUM_VAR_CHEM2D; ++k) {	
	       buffer[buffer_size] = SolnBlk.U[i][j][k+1]; 
	       buffer_size = buffer_size + 1;
		}
	   } /* endfor */
	 } /* endfor */
       } /* endif */

       buffer_size = NUM_VAR_CHEM2D*ni*nj;
       Communicator.Bcast(buffer, buffer_size, MPI::DOUBLE, Source_Rank);

       if (!(CFDkit_MPI::This_Processor_Number == Source_CPU)) {
	 buffer_size = 0;
          for (int j  = SolnBlk.JCl-SolnBlk.Nghost ; j <= SolnBlk.JCu+SolnBlk.Nghost ; ++j ) {
	    for (int i = SolnBlk.ICl-SolnBlk.Nghost ; i <= SolnBlk.ICu+SolnBlk.Nghost ; ++i ) {
	      for ( int k=0; k<NUM_VAR_CHEM2D; ++k) {	       
		SolnBlk.U[i][j][k+1]= buffer[buffer_size];
		buffer_size = buffer_size + 1;
	      }
	      SolnBlk.W[i][j] = W(SolnBlk.U[i][j]);
	    } /* endfor */
	  } /* endfor */
       } /* endif */

       delete []buffer; 
       buffer = NULL;

       ni = 1;
       nj = (SolnBlk.JCu+SolnBlk.Nghost) - (SolnBlk.JCl-SolnBlk.Nghost) + 1;
       buffer = new double[2*NUM_VAR_CHEM2D*ni*nj];

       if (CFDkit_MPI::This_Processor_Number == Source_CPU) {
          buffer_size = 0;
	  for (int j  = SolnBlk.JCl-SolnBlk.Nghost ; j <= SolnBlk.JCu+SolnBlk.Nghost ; ++j ) {
	    for ( int k=0; k<NUM_VAR_CHEM2D; ++k) {	
	      buffer[buffer_size]= SolnBlk.WoW[j][k+1];	
	      buffer_size = buffer_size + 1;
	    }
	    for ( int k=0; k<NUM_VAR_CHEM2D; ++k) {	
	      buffer[buffer_size]= SolnBlk.WoE[j][k+1];
	      buffer_size = buffer_size + 1;
	    }
          } /* endfor */
       } /* endif */

       buffer_size = 2*NUM_VAR_CHEM2D*ni*nj;
       Communicator.Bcast(buffer, buffer_size, MPI::DOUBLE, Source_Rank);

       if (!(CFDkit_MPI::This_Processor_Number == Source_CPU)) {
          buffer_size = 0;
	  for (int j  = SolnBlk.JCl-SolnBlk.Nghost ; j <= SolnBlk.JCu+SolnBlk.Nghost ; ++j ) {    
	    for ( int k=0; k<NUM_VAR_CHEM2D; ++k) {	
	       SolnBlk.WoW[j][k+1] = buffer[buffer_size];
	       buffer_size = buffer_size + 1;
	    }
	    for ( int k=0; k<NUM_VAR_CHEM2D; ++k) {
	      SolnBlk.WoE[j][k+1] = buffer[buffer_size];
	      buffer_size = buffer_size + 1;
	    }
          } /* endfor */
       } /* endif */

       delete []buffer; 
       buffer = NULL;
       ni = (SolnBlk.ICu+SolnBlk.Nghost) - (SolnBlk.ICl-SolnBlk.Nghost) + 1;
       nj = 1;
       buffer = new double[2*NUM_VAR_CHEM2D*ni*nj];

       if (CFDkit_MPI::This_Processor_Number == Source_CPU) {
          buffer_size = 0;  
          for (int i  = SolnBlk.ICl-SolnBlk.Nghost ; i <= SolnBlk.ICu+SolnBlk.Nghost ; ++i ) {
	    for ( int k=0; k<NUM_VAR_CHEM2D; ++k) {    
	      buffer[buffer_size]= SolnBlk.WoS[i][k+1];   
	      buffer_size = buffer_size + 1;
	    }
	    for ( int k=0; k<NUM_VAR_CHEM2D; ++k) {    
	      buffer[buffer_size]= SolnBlk.WoN[i][k+1];
	      buffer_size = buffer_size + 1;
	    } 
          } /* endfor */
       } /* endif */
       buffer_size = 2*NUM_VAR_CHEM2D*ni*nj;
       Communicator.Bcast(buffer, buffer_size, MPI::DOUBLE, Source_Rank);

       if (!(CFDkit_MPI::This_Processor_Number == Source_CPU)) {
	 buffer_size = 0;
	  for (int i  = SolnBlk.ICl-SolnBlk.Nghost ; i <= SolnBlk.ICu+SolnBlk.Nghost ; ++i ) {
	    for ( int k=0; k<NUM_VAR_CHEM2D; ++k) {	
	      SolnBlk.WoS[i][k+1] = buffer[buffer_size]; 
	      buffer_size = buffer_size + 1;
	    }
	    for ( int k=0; k<NUM_VAR_CHEM2D; ++k) {	
	      SolnBlk.WoN[i][k+1] = buffer[buffer_size]; 
	      buffer_size = buffer_size + 1;
	    }
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
void Copy_Solution_Block(Chem2D_Quad_Block &SolnBlk1,
                         Chem2D_Quad_Block &SolnBlk2) {

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

    /* Set the axisymmetric/planar flow, viscous, and gravity indicators. */

    SolnBlk1.Flow_Type = SolnBlk2.Flow_Type;
    SolnBlk1.Axisymmetric = SolnBlk2.Axisymmetric;
    SolnBlk1.Gravity = SolnBlk2.Gravity;
    SolnBlk1.debug_level = SolnBlk2.debug_level;
    SolnBlk1.Moving_wall_velocity = SolnBlk2.Moving_wall_velocity; 

    /* Copy the grid of the second solution block
       to the first solution block. */

    Copy_Quad_Block(SolnBlk1.Grid, SolnBlk2.Grid);

    /* Copy the solution information from SolnBlk2 to SolnBlk1. */

    if (SolnBlk2.U != NULL) {
       for ( j  = SolnBlk1.JCl-SolnBlk1.Nghost ; j <= SolnBlk1.JCu+SolnBlk1.Nghost ; ++j ) {
          for ( i = SolnBlk1.ICl-SolnBlk1.Nghost ; i <= SolnBlk1.ICu+SolnBlk1.Nghost ; ++i ) {
             SolnBlk1.U[i][j] = SolnBlk2.U[i][j];
             SolnBlk1.W[i][j] = SolnBlk2.W[i][j];
             for ( k = 0 ; k <= NUMBER_OF_RESIDUAL_VECTORS_CHEM2D-1 ; ++k ) {
	        SolnBlk1.dUdt[i][j][k] = SolnBlk2.dUdt[i][j][k];
             } /* endfor */
	     SolnBlk1.dWdx[i][j] = SolnBlk2.dWdx[i][j];
	     SolnBlk1.dWdy[i][j] = SolnBlk2.dWdy[i][j];
	     SolnBlk1.phi[i][j] = SolnBlk2.phi[i][j];
	     SolnBlk1.Uo[i][j] = SolnBlk2.Uo[i][j];
	     SolnBlk1.dt[i][j] = SolnBlk2.dt[i][j];
          } /* endfor */
       } /* endfor */

       for (j  = SolnBlk1.JCl-SolnBlk1.Nghost ; j <= SolnBlk1.JCu+SolnBlk1.Nghost ; ++j ) {
	   SolnBlk1.WoW[j] = SolnBlk2.WoW[j];
           SolnBlk1.WoE[j] = SolnBlk2.WoE[j];
       } /* endfor */

       for ( i = SolnBlk1.ICl-SolnBlk1.Nghost ; i <= SolnBlk1.ICu+SolnBlk1.Nghost ; ++i ) {
           SolnBlk1.WoS[i] = SolnBlk2.WoS[i];
           SolnBlk1.WoN[i] = SolnBlk2.WoN[i];
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
int Prolong_Solution_Block(Chem2D_Quad_Block &SolnBlk_Fine,
			   Chem2D_Quad_Block &SolnBlk_Original,
			   const int Sector) {

    int i, j, k, i_min, i_max, j_min, j_max, mesh_refinement_permitted;
    double area_total_fine;

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
          // In this case, create the refined mesh. /* 
          Refine_Mesh(SolnBlk_Fine.Grid, 
                      SolnBlk_Original.Grid,
                      Sector);
       } /* endif */
    } /* endif */

    if (mesh_refinement_permitted) {
   
    /* Set the axisymmetric/planar flow, viscsous, and gravity indicator for the fine solution block. */

       SolnBlk_Fine.Flow_Type = SolnBlk_Original.Flow_Type;
       SolnBlk_Fine.Axisymmetric = SolnBlk_Original.Axisymmetric;
       SolnBlk_Fine.Gravity = SolnBlk_Original.Gravity;
       SolnBlk_Fine.debug_level = SolnBlk_Original.debug_level; 
       SolnBlk_Fine.Moving_wall_velocity = SolnBlk_Original.Moving_wall_velocity; 

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
//                   = (SolnBlk_Original.Grid.Cell[i][j].A/area_total_fine)*SolnBlk_Original.U[i][j];
    	       SolnBlk_Fine.W[2*(i-i_min)+SolnBlk_Fine.ICl  ]
                             [2*(j-j_min)+SolnBlk_Fine.JCl  ] 
                  = W(SolnBlk_Fine.U[2*(i-i_min)+SolnBlk_Fine.ICl  ]
                                    [2*(j-j_min)+SolnBlk_Fine.JCl  ]);

    	       SolnBlk_Fine.U[2*(i-i_min)+SolnBlk_Fine.ICl+1]
                             [2*(j-j_min)+SolnBlk_Fine.JCl  ] 
                  = SolnBlk_Original.U[i][j];
//                   = (SolnBlk_Original.Grid.Cell[i][j].A/area_total_fine)*SolnBlk_Original.U[i][j];
    	       SolnBlk_Fine.W[2*(i-i_min)+SolnBlk_Fine.ICl+1]
                             [2*(j-j_min)+SolnBlk_Fine.JCl  ] 
                  = W(SolnBlk_Fine.U[2*(i-i_min)+SolnBlk_Fine.ICl+1]
                                    [2*(j-j_min)+SolnBlk_Fine.JCl  ]);

    	       SolnBlk_Fine.U[2*(i-i_min)+SolnBlk_Fine.ICl  ]
                             [2*(j-j_min)+SolnBlk_Fine.JCl+1]
                  = SolnBlk_Original.U[i][j];
//                   = (SolnBlk_Original.Grid.Cell[i][j].A/area_total_fine)*SolnBlk_Original.U[i][j];
    	       SolnBlk_Fine.W[2*(i-i_min)+SolnBlk_Fine.ICl  ]
                             [2*(j-j_min)+SolnBlk_Fine.JCl+1]
                  = W(SolnBlk_Fine.U[2*(i-i_min)+SolnBlk_Fine.ICl  ]
                                    [2*(j-j_min)+SolnBlk_Fine.JCl+1]);

    	       SolnBlk_Fine.U[2*(i-i_min)+SolnBlk_Fine.ICl+1]
                             [2*(j-j_min)+SolnBlk_Fine.JCl+1]
                  = SolnBlk_Original.U[i][j];
//                   = (SolnBlk_Original.Grid.Cell[i][j].A/area_total_fine)*SolnBlk_Original.U[i][j];
    	       SolnBlk_Fine.W[2*(i-i_min)+SolnBlk_Fine.ICl+1]
                             [2*(j-j_min)+SolnBlk_Fine.JCl+1]
                  = W(SolnBlk_Fine.U[2*(i-i_min)+SolnBlk_Fine.ICl+1]
                                    [2*(j-j_min)+SolnBlk_Fine.JCl+1]);
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

    } /* endif */

    // Prolongation of solution block.
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
int Restrict_Solution_Block(Chem2D_Quad_Block &SolnBlk_Coarse,
			    Chem2D_Quad_Block &SolnBlk_Original_SW,
			    Chem2D_Quad_Block &SolnBlk_Original_SE,
			    Chem2D_Quad_Block &SolnBlk_Original_NW,
			    Chem2D_Quad_Block &SolnBlk_Original_NE) {

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

       /* Set the axisymmetric/planar flow indicator for the coarse solution block. */

       SolnBlk_Coarse.Flow_Type = SolnBlk_Original_SW.Flow_Type;
       SolnBlk_Coarse.Axisymmetric = SolnBlk_Original_SW.Axisymmetric;
       SolnBlk_Coarse.Gravity = SolnBlk_Original_SW.Gravity;
       SolnBlk_Coarse.debug_level = SolnBlk_Original_SW.debug_level;       
       SolnBlk_Coarse.Moving_wall_velocity = SolnBlk_Original_SW.Moving_wall_velocity; 

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
 //                                                     SolnBlk_Coarse.Grid.Cell[i_coarse][j_coarse].A;
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
//                                                      SolnBlk_Coarse.Grid.Cell[i_coarse][j_coarse].A;
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
//                                                      SolnBlk_Coarse.Grid.Cell[i_coarse][j_coarse].A;
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

      // North-east corner fine block:
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
//                                                      SolnBlk_Coarse.Grid.Cell[i_coarse][j_coarse].A;
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

   } /* endif */

    // Restriction of solution blocks successful.
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
void Output_Tecplot(Chem2D_Quad_Block &SolnBlk,
		    Chem2D_Input_Parameters &Input_Parameters,
                    const int Number_of_Time_Steps,
                    const double &Time,
                    const int Block_Number,
                    const int Output_Title,
	            ostream &Out_File) {

  Chem2D_pState W_node;
   
  /* Cell centered shear and qflux */
  if (SolnBlk.Flow_Type != FLOWTYPE_INVISCID) {
    Viscous_Calculations(SolnBlk);
  }

  /* Ensure boundary conditions are updated before
     evaluating solution at the nodes. */
  BCs(SolnBlk,Input_Parameters);

  /* Output node solution data. */
  Out_File << setprecision(14);
  if (Output_Title) {
 
    Out_File << "TITLE = \"" << CFDkit_Name() << ": 2D Chem2D Solution, "
	     << "Time Step/Iteration Level = " << Number_of_Time_Steps
	     << ", Time = " << Time
	     << "\"" << "\n"
	     << "VARIABLES = \"x\" \\ \n"
	     << "\"y\" \\ \n"
	     << "\"rho\" \\ \n"
	     << "\"u\" \\ \n"
	     << "\"v\" \\ \n"
	     << "\"p\" \\ \n"
	     << "\"k\" \\ \n"
	     << "\"omega\" \\ \n";
    //n species mass fractions names
    for(int i =0 ;i<SolnBlk.W[0][0].ns ;i++){
      Out_File <<"\"c_"<<SolnBlk.W[0][0].specdata[i].Speciesname()<<"\" \\ \n";
    }  
    //   Viscous Terms 
    Out_File << "\"qflux_x\" \\ \n"  
	     << "\"qflux_y\" \\  \n"   
	     << "\"Tau_xx\" \\  \n"  //rr -axisymmetric
	     << "\"Tau_xy\" \\  \n"  //rz
	     << "\"Tau_yy\" \\  \n"  //zz
	     << "\"Tau_zz\" \\  \n"
	     << "\"theta_x\" \\  \n"  
	     << "\"theta_y\" \\  \n"   
	     << "\"lambda_xx\" \\  \n"  //rr -axisymmetric
	     << "\"lambda_xy\" \\  \n"  //rz
	     << "\"lambda_yy\" \\  \n"  //zz
	     << "\"lambda_zz\" \\  \n";   
    //  Calculated values
    Out_File<< "\"T\" \\  \n" 
	    << "\"R\" \\  \n"
	    << "\"M\" \\  \n"
            << "\"viscosity\" \\  \n"
            << "\"thermal conduct\" \\  \n"
	    << "\"Prandtl\" \\  \n"
	    << "\"Prandtl_turb\" \\  \n"
	    << "\"rho*H\"  \\ \n"  
	    <<"\"h\" \\ \n"
	    <<"\"h_s\" \\ \n"
	    <<"\"h_ref\" \\ \n"
	    <<"\"rho*E\" \\ \n"
	    << "\"e\" \\  \n" 
	    << "\"e_s\" \\ \n"
	    << "\"e_ref\" \\ \n"
      // Zone details
	    << "ZONE T =  \"Block Number = " << Block_Number
	    << "\" \\ \n"
	    << "I = " << SolnBlk.Grid.INu - SolnBlk.Grid.INl + 1 << " \\ \n"
	    << "J = " << SolnBlk.Grid.JNu - SolnBlk.Grid.JNl + 1 << " \\ \n"
	    << "F = POINT \n";
  }
  // each block zones 
  else {
 
    Out_File << "ZONE T =  \"Block Number = " << Block_Number
	     << "\" \\ \n"
	     << "I = " << SolnBlk.Grid.INu - SolnBlk.Grid.INl + 1 << " \\ \n"
	     << "J = " << SolnBlk.Grid.JNu - SolnBlk.Grid.JNl + 1 << " \\ \n"
	     << "F = POINT \n";
  } /* endif */

  for ( int j  = SolnBlk.Grid.JNl ; j <= SolnBlk.Grid.JNu ; ++j ) {
    for ( int i = SolnBlk.Grid.INl ; i <= SolnBlk.Grid.INu ; ++i ) {
      W_node = SolnBlk.Wn(i, j);
      //coordinates, cell properties
      Out_File << " " << SolnBlk.Grid.Node[i][j].X<<endl<<W_node;
      //T,M,H,s, and all the rest of the calcuatated parameters
      Out_File.setf(ios::scientific);
      Out_File << " " << W_node.T()<< " " << W_node.Rtot()
	       << " " << W_node.v.abs()/W_node.a() 
	       << " " << W_node.mu() <<" "<< W_node.kappa()
	       << " " << W_node.Prandtl()
	       << " " << W_node.Pr_turb()
	       << " " << W_node.H() <<" "<< W_node.h() <<" "<< W_node.hs()<<" "<< W_node.href()
	       << " " << W_node.E() <<" "<< W_node.e() <<" "<< W_node.es()<<" "<< W_node.eref()<<endl;
	     //   <<" "<< W_node.h() <<" "<< W_node.hs()<<" "<< W_node.href()
// 	       << " " << W_node.E() <<" "<< W_node.e() <<" "<< W_node.es()<<" "<< W_node.eref()<<endl;
  
      Out_File.unsetf(ios::scientific);
      
    } /* endfor */
  } /* endfor */
  Out_File << setprecision(6);
  
  //cout<<" end of Out_Tecplot \n";
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
void Output_Cells_Tecplot(Chem2D_Quad_Block &SolnBlk,
                          const int Number_of_Time_Steps,
                          const double &Time,
                          const int Block_Number,
                          const int Output_Title,
	                  ostream &Out_File) {

    int i, j;

    /* Cell centered shear and qflux */
    if (SolnBlk.Flow_Type != FLOWTYPE_INVISCID) {
      Viscous_Calculations(SolnBlk);
    }

    Out_File << setprecision(14);
    if (Output_Title) {
       Out_File << "TITLE = \"" << CFDkit_Name() << ": 2D Chem2D Solution, "
                << "Time Step/Iteration Level = " << Number_of_Time_Steps
                << ", Time = " << Time
                << "\"" << "\n"
   	        << "VARIABLES = \"x\" \\ \n"
                << "\"y\" \\ \n"
                << "\"rho\" \\ \n"
                << "\"u\" \\ \n"
                << "\"v\" \\ \n"
		<< "\"p\" \\ \n"
		<< "\"k\" \\ \n"
		<< "\"omega\" \\ \n";
       //n species mass fractions names
       for(int i =0 ;i<SolnBlk.W[0][0].ns ;i++){
	 Out_File <<"\"c"<<SolnBlk.W[0][0].specdata[i].Speciesname()<<"\" \\ \n";
       }
       //Viscous Terms 
       Out_File << "\"qflux_x\" \\ \n"  
		<< "\"qflux_y\" \\ \n"   
		<< "\"Tau_xx\" \\ \n"  //rr -axisymmetric
		<< "\"Tau_xy\" \\ \n"  //rz
		<< "\"Tau_yy\" \\ \n"  //zz
		<< "\"Tau_zz\" \\ \n"  //thetatheta
		<< "\"theta_x\" \\ \n"  
		<< "\"theta_y\" \\ \n"   
		<< "\"lambda_xx\" \\ \n"   //rr -axisymmetric
		<< "\"lambda_xy\" \\ \n"   //rz
		<< "\"lambda_yy\" \\ \n"   //zz
		<< "\"lambda_zz\" \\ \n";  //thetatheta
        //Calculated values
       Out_File << "\"T\" \\ \n"
	        << "\"M\" \\ \n"
	        << "\"viscosity\" \\ \n"
	        << "\"thermal conduct\" \\ \n"
	        << "\"Prandtl\" \\ \n";
       //Prandtl, Schmidt, Lewis
       for(int i =0 ;i<SolnBlk.W[0][0].ns;i++){
	 Out_File<<"\"Sc_"<<SolnBlk.W[0][0].specdata[i].Speciesname()<<"\" \\ \n"
		 <<"\"Le_"<<SolnBlk.W[0][0].specdata[i].Speciesname()<<"\" \\ \n"; 
       } 
       
//        //Limiters 
//        for(int i =1 ;i<= SolnBlk.NumVar();i++){
// 	    Out_File<<"\"phi_"<<i<<"\" \\ \n";
//        }       
//        //other stored cell-center data
//        for(int i =0 ;i<SolnBlk.W[0][0].ns;i++){
// 	    Out_File<<"\"gradc.x_"<<SolnBlk.W[0][0].specdata[i].Speciesname()<<"\" \\ \n"
// 		  <<"\"gradc.y_"<<SolnBlk.W[0][0].specdata[i].Speciesname()<<"\" \\ \n"
// 		  <<"\"Dif_coef_"<<SolnBlk.W[0][0].specdata[i].Speciesname()<<"\" \\ \n"; 
//        }       
//        //gradients
//        Out_File<< "\"dx_rho\" \\ \n"
// 	       << "\"dx_u\" \\ \n"
// 	       << "\"dx_v\" \\ \n"
// 	       << "\"dx_p\" \\ \n";
//        //n species mass fractions names
//        for(int i =0 ;i<SolnBlk.W[0][0].ns ;i++){
// 	     Out_File <<"\"dx_c"<<SolnBlk.W[0][0].specdata[i].Speciesname()<<"\" \\ \n";
//        }  
//        Out_File<< "\"dx_q_x\" \n"  
// 	       << "\"dx_q_y\" \n"  
// 	       << "\"dx_Tau_xx\" \n"  //rr -axisymmetric
// 	       << "\"dx_Tau_xy\" \n"  //rz
// 	       << "\"dx_Tau_yy\" \n"  //zz
// 	       << "\"dx_Tau_zz\" \n"; //thetatheta       
//        Out_File<< "\"dy_rho\" \\ \n"
// 	       << "\"dy_u\" \\ \n"
// 	       << "\"dy_v\" \\ \n"
// 	       << "\"dy_p\" \\ \n";
//        //n species mass fractions names
//        for(int i =0 ;i<SolnBlk.W[0][0].ns ;i++){
//  	     Out_File <<"\"dy_c"<<SolnBlk.W[0][0].specdata[i].Speciesname()<<"\" \\ \n";
//        }  
//        Out_File<< "\"dy_q_x\" \n"  
// 	       << "\"dy_q_y\" \n"   
// 	       << "\"dy_Tau_xx\" \n"  //rr -axisymmetric
// 	       << "\"dy_Tau_xy\" \n"  //rz
// 	       << "\"dy_Tau_yy\" \n"  //zz
// 	       << "\"dy_Tau_zz\" \n"; //thetatheta

       // Zone details
       Out_File<< "ZONE T =  \"Block Number = " << Block_Number
	       << "\" \\ \n"
	       << "I = " << SolnBlk.Grid.ICu - SolnBlk.Grid.ICl + 2*SolnBlk.Nghost + 1 << " \\ \n"
	       << "J = " << SolnBlk.Grid.JCu - SolnBlk.Grid.JCl + 2*SolnBlk.Nghost + 1<< " \\ \n"
	       << "F = POINT \n";    
    } else {
      Out_File << "ZONE T =  \"Block Number = " << Block_Number
	       << "\" \\ \n"
                << "I = " << SolnBlk.Grid.ICu - SolnBlk.Grid.ICl + 2*SolnBlk.Nghost + 1 << " \\ \n"
                << "J = " << SolnBlk.Grid.JCu - SolnBlk.Grid.JCl + 2*SolnBlk.Nghost + 1 << " \\ \n"
                << "F = POINT \n";
    } /* endif */

    for ( j  = SolnBlk.JCl-2 ; j <= SolnBlk.JCu+2 ; ++j ) {
       for ( i = SolnBlk.ICl-2 ; i <= SolnBlk.ICu+2 ; ++i ) {
           Out_File << " "  << SolnBlk.Grid.Cell[i][j].Xc
                    << SolnBlk.W[i][j];
           Out_File.setf(ios::scientific);
	   //Temperature
	   Out_File << " " << SolnBlk.W[i][j].T()
		    << " " << SolnBlk.W[i][j].v.abs()/SolnBlk.W[i][j].a() 
		    << " " << SolnBlk.W[i][j].mu()
		    << " " << SolnBlk.W[i][j].kappa()
		    << " " << SolnBlk.W[i][j].Prandtl(); 
	 //   //Prandtl, Schmidt, Lewis   
	   for(int k =0 ;k<SolnBlk.W[0][0].ns ;k++){	  
 	       Out_File<<" "<<SolnBlk.W[i][j].Schmidt_No(k) 
		       <<" "<<SolnBlk.W[i][j].Lewis(k);  
	   }
	  		 
//	   //limiters
// 	   for(int k =1; k<=SolnBlk.NumVar(); k++){
// 	     Out_File <<" "<<SolnBlk.phi[i][j][k];
// 	   }
// 	   //coef..
// 	   for(int k =0 ;k<SolnBlk.W[0][0].ns ;k++){
// 	     Out_File <<SolnBlk.W[i][j].spec[k].gradc<<" "
// 		      <<SolnBlk.W[i][j].spec[k].diffusion_coef;
// 	   }
// 	   //gradients
// 	   Out_File <<SolnBlk.dWdx[i][j]
// 		    <<SolnBlk.dWdy[i][j];
	   Out_File<< "\n";
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
void Output_Nodes_Tecplot(Chem2D_Quad_Block &SolnBlk,
                          const int Number_of_Time_Steps,
                          const double &Time,
                          const int Block_Number,
                          const int Output_Title,
	                  ostream &Out_File) {

    int i, j;

    Out_File << setprecision(14);
    if (Output_Title) {
       Out_File << "TITLE = \"" << CFDkit_Name() << ": 2D Chem2D Solution, "
                << "Time Step/Iteration Level = " << Number_of_Time_Steps
                << ", Time = " << Time
                << "\"" << "\n"
   	        << "VARIABLES = \"x\" \\ \n";
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

/********************************************************
 * Routine: Output_RHS                                  *
 * For testing Jacobians                                *
 *                                                      *
 ********************************************************/
void Output_RHS(Chem2D_Quad_Block &SolnBlk,
		const int Number_of_Time_Steps,
		const double &Time,
		const int Block_Number,
		const int Output_Title,
		ostream &Out_File) {

//   Chem2D_pState W_node;
//   /* Ensure boundary conditions are updated before
//      evaluating solution at the nodes. */
  
//   BCs(SolnBlk);
 
//   //   Output node solution data.
//   Out_File<<"VARIABLES =   r  , rho  , u  , k  , omega  , c1,  c2,  T  "<<endl;
//   Out_File<<"ZONE J=  "<< (SolnBlk.Grid.JCu-1) <<", F=POINT"<<endl;
//   Out_File << setprecision(14);

//   double rad = 0.1;
//   for ( int j  = SolnBlk.Grid.JCl ; j <= SolnBlk.Grid.JCu ; ++j ) {
//     Out_File << " " << SolnBlk.Grid.Cell[6][j].Xc.y/rad <<"  "<<SolnBlk.W[6][j].rho<<"  "<<SolnBlk.W[6][j].v.x
// 	     <<"  "<<SolnBlk.W[6][j].k<<"  "<<SolnBlk.W[6][j].omega<<"  "<<SolnBlk.W[6][j].spec[0].c
// 	     <<"  "<<SolnBlk.W[6][j].spec[1].c<<"  " <<SolnBlk.W[6][j].T()<<endl;
    
//     Out_File.unsetf(ios::scientific);
    
//   } /* endfor */

//   Out_File<<"ZONE J=  "<< (SolnBlk.Grid.JCu-1) <<", F=POINT"<<endl;
//   Out_File << setprecision(14);
//   for ( int j  = SolnBlk.Grid.JCl ; j <= SolnBlk.Grid.JCu ; ++j ) {
//     Out_File << " " << SolnBlk.Grid.Cell[16][j].Xc.y/rad<<"  "<<SolnBlk.W[16][j].rho<<"  "<<SolnBlk.W[16][j].v.x
// 	     <<"  "<<SolnBlk.W[16][j].k<<"  "<<SolnBlk.W[16][j].omega
// 	     <<"  "<<SolnBlk.W[16][j].spec[0].c<<"  "<<SolnBlk.W[16][j].spec[1].c<<"  " <<SolnBlk.W[16][j].T()<<endl;
    
//     Out_File.unsetf(ios::scientific);
    
    
//   } /* endfor */  

//   Out_File<<"ZONE J=  "<< (SolnBlk.Grid.JCu-1) <<", F=POINT"<<endl;
//   Out_File << setprecision(14);
//   for ( int j  = SolnBlk.Grid.JCl ; j <= SolnBlk.Grid.JCu ; ++j ) {
//     Out_File << " " << SolnBlk.Grid.Cell[36][j].Xc.y/rad <<"  "<<SolnBlk.W[36][j].rho<<"  "<<SolnBlk.W[36][j].v.x
// 	     <<"  "<<SolnBlk.W[36][j].k<<"  "<<SolnBlk.W[36][j].omega
// 	     <<"  "<<SolnBlk.W[36][j].spec[0].c<<"  "<<SolnBlk.W[36][j].spec[1].c<<"  " <<SolnBlk.W[36][j].T()<<endl;
    
//     Out_File.unsetf(ios::scientific);
    
  
//   } /* endfor */  


  Out_File << setprecision(6);
  
  //cout<<" end of Out_Tecplot \n";
}
 
/********************************************************
 * Routine: ICs                                         *
 *                                                      *
 * Assigns initial conditions and data to the           *
 * solution variables of the specified quadrilateral    *
 * solution block.                                      *
 *                                                      *
 ********************************************************/
void ICs(Chem2D_Quad_Block &SolnBlk,
	 const int i_ICtype,
	 Chem2D_pState *Wo, Chem2D_Input_Parameters &Input_Parameters) {

    Chem2D_pState Wl, Wr;
    Chem2D_pState Chem2D_W_STDATM;
    Chem2D_cState Chem2D_U_STDATM; //default
    Vector2D dX;

    double delta_pres = ZERO; //should be input parameter but I am lazy     
    double fuel_spacing, fuel_velocity, fuel_temp_inlet, ignition_temp;
    double air_spacing, air_velocity, air_temp_inlet,tube_thickness ;
    double friction_vel = 1.587; //friction velocity for turbulent pipe flow
    double eta, f, fp, fpp;

    /* Assign the initial data for the IVP of interest. */
    switch(i_ICtype) {
    case IC_CONSTANT :
    case IC_UNIFORM :
      // Set the solution state to the initial state Wo[0].
      for (int j  = SolnBlk.JCl-SolnBlk.Nghost ; j <= SolnBlk.JCu+SolnBlk.Nghost ; ++j ) {
	for ( int i = SolnBlk.ICl-SolnBlk.Nghost ; i <= SolnBlk.ICu+SolnBlk.Nghost ; ++i ) {
	  SolnBlk.W[i][j] = Wo[0];
	  SolnBlk.U[i][j] = U(SolnBlk.W[i][j]);
	} 
      }
      break;

    case IC_VISCOUS_PIPE_FLOW :
      delta_pres = 635.54; // -635.00  -423.33  -211.67  0.00  211.67  423.33  635.00
      for ( int j = SolnBlk.JCl-SolnBlk.Nghost; j <= SolnBlk.JCu+SolnBlk.Nghost; j++) {
	for ( int i = SolnBlk.ICl-SolnBlk.Nghost; i <= SolnBlk.ICu+SolnBlk.Nghost; i++) {
          // Apply uniform solution state
	  SolnBlk.W[i][j] = Wo[0];
	  //velocity
	  SolnBlk.W[i][j].v.x = 0.0;
	  
	  if (SolnBlk.Flow_Type == FLOWTYPE_TURBULENT_RANS_K_OMEGA ||
              SolnBlk.Flow_Type == FLOWTYPE_TURBULENT_RANS_K_EPSILON){
	    //turbulent kinetic energy //MARKTHIS XINFENG
	    SolnBlk.W[i][j].k = friction_vel*friction_vel/sqrt(9.0/100.0);
	  	    //specific dissipation rate
	    double radius;
	    if(SolnBlk.Axisymmetric ==1){
	      //pressure 
	      SolnBlk.W[i][j].p = Wo[0].p - ((i-SolnBlk.ICl)*delta_pres/(SolnBlk.ICu-SolnBlk.ICl+1));
	      radius = SolnBlk.Grid.Cell[i][j].Xc.y;
	      if (j < SolnBlk.JCl) {
		SolnBlk.W[i][j].k = TOLER;//~ZERO;
	      }
	      if (j > SolnBlk.JCu){
		SolnBlk.W[i][j].k = TOLER;//~ZERO;
	      }
	    }
	    if(SolnBlk.Axisymmetric ==2){
	      SolnBlk.W[i][j].p = Wo[0].p - ((j-SolnBlk.JCl)*delta_pres/(SolnBlk.JCu-SolnBlk.JCl+1));
	      radius = SolnBlk.Grid.Cell[i][j].Xc.x;
	      if (i < SolnBlk.ICl) {
		SolnBlk.W[i][j].k = TOLER;//~ZERO;
	      }
	      if (i > SolnBlk.ICu){
		SolnBlk.W[i][j].k = TOLER;//~ZERO;
	      }
	    }
	    SolnBlk.W[i][j].omega = friction_vel/
	      (sqrt(9.0/100.0)*0.41*(Input_Parameters.Pipe_Radius-radius)); 
	    if (Input_Parameters.i_Turbulence_BCs == TURBULENT_BC_DIRECT_INTEGRATION &&
		SolnBlk.Wall[i][j].yplus < Input_Parameters.yplus_sublayer){
	      SolnBlk.W[i][j].omega = SolnBlk.W[i][j].omega_sublayer_BC(radius);
	    } //for turbulent flow
	  }
 
          //conservative solution state
	  SolnBlk.U[i][j]   = U(SolnBlk.W[i][j]);
	}
      }
      //omega set-up for ghost cells j direction
      for ( int i = SolnBlk.ICl-SolnBlk.Nghost; i <= SolnBlk.ICu+SolnBlk.Nghost; i++) {
	SolnBlk.W[i][ SolnBlk.JCl-1].omega =  SolnBlk.W[i][ SolnBlk.JCl].omega;
	SolnBlk.W[i][ SolnBlk.JCl-2].omega =  SolnBlk.W[i][ SolnBlk.JCl+1].omega;
	SolnBlk.W[i][ SolnBlk.JCu+1].omega =  SolnBlk.W[i][ SolnBlk.JCu].omega;
	SolnBlk.W[i][ SolnBlk.JCu+2].omega =  SolnBlk.W[i][ SolnBlk.JCu-1].omega;

	SolnBlk.U[i][ SolnBlk.JCl-1]   = U(SolnBlk.W[i][ SolnBlk.JCl-1]);
	SolnBlk.U[i][ SolnBlk.JCl-2]   = U(SolnBlk.W[i][ SolnBlk.JCl-2]);
	SolnBlk.U[i][ SolnBlk.JCu+1]   = U(SolnBlk.W[i][ SolnBlk.JCu+1]);
	SolnBlk.U[i][ SolnBlk.JCu+2]   = U(SolnBlk.W[i][ SolnBlk.JCu+2]);
      }
      
      break;
     
    case IC_SOD_XDIR :
      // Set initial data for Sod IVP in x-direction.
      Wl = Wo[0];
      Wr = Wo[0];
      Wr.rho = Wo[0].rho/EIGHT;
      Wr.p = Wo[0].p/TEN;

      for (int j  = SolnBlk.JCl-SolnBlk.Nghost ; j <= SolnBlk.JCu+SolnBlk.Nghost ; ++j ) {
	for ( int i = SolnBlk.ICl-SolnBlk.Nghost ; i <= SolnBlk.ICu+SolnBlk.Nghost ; ++i ) {
	  if (SolnBlk.Grid.Cell[i][j].Xc.x <= ZERO) {
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
      Wl = Wo[0]; Wr = Wo[0];
      Wl.rho = 2.281; Wl.v.x = 164.83; Wl.v.y = ZERO; Wl.p = 201.17e03;
      Wr.rho = 1.408; Wr.v.x = ZERO; Wr.v.y = ZERO; Wr.p = 101.1e03;
      for (int j  = SolnBlk.JCl-SolnBlk.Nghost ; j <= SolnBlk.JCu+SolnBlk.Nghost ; ++j ) {
	for ( int i = SolnBlk.ICl-SolnBlk.Nghost ; i <= SolnBlk.ICu+SolnBlk.Nghost ; ++i ) {
	  if (SolnBlk.Grid.Cell[i][j].Xc.x <= ZERO) {
	    SolnBlk.W[i][j] = Wl;
	  } else {
	    SolnBlk.W[i][j] = Wr;	
	  } /* end if */
	  SolnBlk.U[i][j] = U(SolnBlk.W[i][j]);
	} /* endfor */
      } /* endfor */
      break;
      
    case IC_PRESSURE_GRADIENT_X :
      delta_pres = 635.54;
      //starts with linear pressure gradient 
      for ( int j  = SolnBlk.JCl-SolnBlk.Nghost ; j <= SolnBlk.JCu+SolnBlk.Nghost ; ++j ) {
	for (int  i = SolnBlk.ICl-SolnBlk.Nghost ; i <= SolnBlk.ICu+SolnBlk.Nghost ; ++i ) {
	  SolnBlk.W[i][j] = Wo[0];
	  SolnBlk.W[i][j].p = Wo[0].p - (i-SolnBlk.ICl-SolnBlk.Nghost)*delta_pres/(SolnBlk.ICu-SolnBlk.ICl);	 
	  //start
	  if( i == SolnBlk.ICl-SolnBlk.Nghost || i == SolnBlk.ICl-SolnBlk.Nghost+1 || i ==SolnBlk.ICl){
	    SolnBlk.W[i][j].p = Wo[0].p;
	  }
	  //end
	  if( i == SolnBlk.ICu+SolnBlk.Nghost || i == SolnBlk.ICu+SolnBlk.Nghost-1 || i ==SolnBlk.ICu){
	    SolnBlk.W[i][j].p = Wo[0].p - delta_pres; 
	  }	   
	  SolnBlk.U[i][j] = U(SolnBlk.W[i][j]);
	} 
      }       
      break; 

    case IC_PRESSURE_GRADIENT_Y :
      delta_pres = 635.54;
      //starts with linear pressure gradient 
      for ( int j  = SolnBlk.JCl-SolnBlk.Nghost ; j <= SolnBlk.JCu+SolnBlk.Nghost ; ++j ) {
	for (int  i = SolnBlk.ICl-SolnBlk.Nghost ; i <= SolnBlk.ICu+SolnBlk.Nghost ; ++i ) {
	  SolnBlk.W[i][j] = Wo[0];
	  SolnBlk.W[i][j].p = Wo[0].p - (j-SolnBlk.JCl-SolnBlk.Nghost)*delta_pres/(SolnBlk.JCu-SolnBlk.JCl);	 
	  //start
	  if( j == SolnBlk.JCl-SolnBlk.Nghost || j == SolnBlk.JCl-SolnBlk.Nghost+1 || j ==SolnBlk.JCl){
	    SolnBlk.W[i][j].p = Wo[0].p;
	  }
	  //end
	  if( j == SolnBlk.JCu+SolnBlk.Nghost || j == SolnBlk.JCu+SolnBlk.Nghost-1 || j ==SolnBlk.JCu){
	    SolnBlk.W[i][j].p = Wo[0].p - delta_pres; 
	  }	   
	  SolnBlk.U[i][j] = U(SolnBlk.W[i][j]);
	} 
      }   

      //Velocity profile
  
      break; 
      
      /***********************************************************************
       * The following three are classic Viscous test cases.                 * 
       *                                                                     *
       * Starting with exact solution, assuming length = 0.2m, height=0.001m *
       ***********************************************************************/
//     case IC_VISCOUS_COUETTE: 
//       //Coutte flow with no pressure gradient        
//       for (int j = SolnBlk.JCl-SolnBlk.Nghost; j <= SolnBlk.JCu+SolnBlk.Nghost; j++) {
// 	for (int i = SolnBlk.ICl-SolnBlk.Nghost; i <= SolnBlk.ICu+SolnBlk.Nghost; i++) {
// 	  SolnBlk.W[i][j] = Wo[0];
// 	  if (j >= SolnBlk.JCl && j <= SolnBlk.JCu){
// 	    SolnBlk.W[i][j].v.x =  SolnBlk.Moving_wall_velocity*
// 	      ((SolnBlk.Grid.Cell[i][j].Xc.y + 0.0005)/0.001);
// 	  }
// 	  SolnBlk.U[i][j] = U(SolnBlk.W[i][j]);	  
// 	}
//       }    
//       break;

    case IC_VISCOUS_COUETTE_PRESSURE_GRADIENT:
      //Couette flow with pressure gradient
      // -635.00  -423.33  -211.67  0.00  211.67  423.33  635.00
      delta_pres = 635.54; //should be Input_Parameter.delta_pres;
      for (int j = SolnBlk.JCl-SolnBlk.Nghost; j <= SolnBlk.JCu+SolnBlk.Nghost; j++) {
	for (int i = SolnBlk.ICl-SolnBlk.Nghost; i <= SolnBlk.ICu+SolnBlk.Nghost; i++) {
	  SolnBlk.W[i][j] = Wo[0];

	  dX.x = fabs(SolnBlk.Grid.Cell[SolnBlk.ICl][j].Xc.x - SolnBlk.Grid.Cell[SolnBlk.ICl-1][j].Xc.x);
	  
	  //start
	  if( i == SolnBlk.ICl-SolnBlk.Nghost){
	    SolnBlk.W[i][j].p = Wo[0].p + delta_pres*5*dX.x*2;	 
	  } else if (i == SolnBlk.ICl-SolnBlk.Nghost+1){
	    SolnBlk.W[i][j].p = Wo[0].p + delta_pres*5*dX.x;	 
	  //end
	  } else if( i == SolnBlk.ICu+SolnBlk.Nghost){
	    SolnBlk.W[i][j].p = Wo[0].p - delta_pres - delta_pres*5*dX.x*2;	 
	  } else if (i == SolnBlk.ICu+SolnBlk.Nghost-1){
	    SolnBlk.W[i][j].p = Wo[0].p - delta_pres - delta_pres*5*dX.x; 
	  } else {	   
	    SolnBlk.W[i][j].p = Wo[0].p - double(i-SolnBlk.ICl)*delta_pres/double(SolnBlk.ICu-SolnBlk.ICl); 
	  }

	  if (j >= SolnBlk.JCl && j <= SolnBlk.JCu){
	    SolnBlk.W[i][j].v.x = (HALF/SolnBlk.W[i][j].mu())*(-delta_pres/0.2)*
	      (pow(SolnBlk.Grid.Cell[i][j].Xc.y,TWO) -(0.001*0.001/4.0))
	      	      + SolnBlk.Moving_wall_velocity*(SolnBlk.Grid.Cell[i][j].Xc.y/0.001 + 0.5); 
	  }
	  SolnBlk.U[i][j]   = U(SolnBlk.W[i][j]);
	}
      }
      break;

    case IC_VISCOUS_FLAT_PLATE :
      for (int j = SolnBlk.JCl-SolnBlk.Nghost; j <= SolnBlk.JCu+SolnBlk.Nghost; j++) {
	for (int i = SolnBlk.ICl-SolnBlk.Nghost; i <= SolnBlk.ICu+SolnBlk.Nghost; i++) {	
	  if (SolnBlk.Grid.Cell[i][j].Xc.x <= ZERO) {
	    SolnBlk.W[i][j] = Wo[0];
	  } else {
	    // Compute the Blasius solution at point Xc.
	    SolnBlk.W[i][j] = FlatPlate(Wo[0],SolnBlk.Grid.Cell[i][j].Xc,eta,f,fp,fpp);
	  }
	  SolnBlk.U[i][j] = U(SolnBlk.W[i][j]);
	}
      }
      break;

    case IC_RINGLEB_FLOW :
      for (int j = SolnBlk.JCl-SolnBlk.Nghost; j <= SolnBlk.JCu+SolnBlk.Nghost; j++) {
	for (int i = SolnBlk.ICl-SolnBlk.Nghost; i <= SolnBlk.ICu+SolnBlk.Nghost; i++) {
	  SolnBlk.W[i][j] = RinglebFlow(SolnBlk.W[i][j],SolnBlk.Grid.Cell[i][j].Xc);
	  SolnBlk.U[i][j] = U(SolnBlk.W[i][j]);
	}
      }
      break;

    case IC_SHOCK_BOX :
        // Set initial data for Aki shock-box IVP.
        Wl = Wo[0];
	Wl.v.x = ZERO; Wl.v.y=ZERO;
        Wr = Wo[0];
	Wr.rho = Wr.rho*FOUR;
	Wr.v.x = ZERO; Wr.v.y = ZERO;
	Wr.p = Wr.p*FOUR;
        for (int j  = SolnBlk.JCl-SolnBlk.Nghost ; j <= SolnBlk.JCu+SolnBlk.Nghost ; ++j ) {
            for (int i = SolnBlk.ICl-SolnBlk.Nghost ; i <= SolnBlk.ICu+SolnBlk.Nghost ; ++i ) {
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
	
     /***************************************************/
     /***************** CHEM2D***************************/
     /***************************************************/     
    case IC_GAS_MIX :
      //set initial data with box half one gas/ half another
      // in this case the opposite mixture rules
      Wl = Wo[0]; Wl.v.y = 0.0; Wl.v.x=1.0;
      Wr = Wo[0]; Wr.v.y = 0.0; Wr.v.x=1.0;
      //Wr.p = Wo[0].p*0.9;

      for(int k=0; k<Wo[0].ns; k++){
	Wl.spec[k] = Wo[0].spec[Wo[0].ns-1-k];
      }

      for (int j  = SolnBlk.JCl-SolnBlk.Nghost ; j <= SolnBlk.JCu+SolnBlk.Nghost ; ++j ) {
	for ( int i = SolnBlk.ICl-SolnBlk.Nghost ; i <= SolnBlk.ICu+SolnBlk.Nghost ; ++i ) {
	  if (SolnBlk.Grid.Cell[i][j].Xc.x <= ZERO){
	    SolnBlk.W[i][j] = Wl;
	  } else {
	    SolnBlk.W[i][j] = Wr;
	  } /* end if */
	  SolnBlk.U[i][j] = U(SolnBlk.W[i][j]);
	} /* endfor */
      } /* endfor */
      break;
      
      /************************************************
       Flame Speed Steady State Calculation  
        
        -> Set initial data with left stoichiometric premixed gas
           and the right hand side being the associated combusted
           products  
      ************************************************/
    case IC_CHEM_1DFLAME:

      Wl = Wo[0]; 
      Wr = Wo[0]; 
 
      //Set downstream products, flamespeeds, Temperatures,
      if(Wo[0].React.reactset_flag == CH4_1STEP){
	//set to phi=1.0 values from CHEMKIN
	Wl.v.x = 0.4101;	
	Wr.v.x = 3.103; 

// 	//phi = 0.6
// 	Wr.spec[0] = ZERO;       //CH4
// 	Wr.spec[1] = 0.0904;     //O2
// 	Wr.spec[2] = 0.0927;     //CO2
// 	Wr.spec[3] = 0.0759;     //H2O 
// 	Wr.rho = Wr.p/(Wr.Rtot()*1650); 

	//phi = 0.8
// 	Wr.spec[0] = ZERO;       //CH4
// 	Wr.spec[1] = 0.0452;     //O2
// 	Wr.spec[2] = 0.1221;     //CO2
// 	Wr.spec[3] = 0.0999;     //H2O 
// 	Wr.rho = Wr.p/(Wr.Rtot()*2000); 

	//phi = 1.0
	Wr.spec[0] = ZERO;       //CH4
	Wr.spec[1] = ZERO;       //O2
	Wr.spec[2] = 0.1511;     //CO2
	Wr.spec[3] = 0.1242;     //H2O 
	Wr.rho = Wr.p/(Wr.Rtot()*2320); //2234
	
//  	//phi = 1.2
// 	Wr.spec[0] = 0.0107;     //CH4
// 	Wr.spec[1] = ZERO;       //O2
// 	Wr.spec[2] = 0.1498;     //CO2
// 	Wr.spec[3] = 0.1227;     //H2O 
// 	Wr.rho = Wr.p/(Wr.Rtot()*2260);

//   	//phi = 1.5
// 	Wr.spec[0] = 0.0267;     //CH4
// 	Wr.spec[1] = ZERO;       //O2
// 	Wr.spec[2] = 0.1474;     //CO2
// 	Wr.spec[3] = 0.1207;     //H2O 
// 	Wr.rho = Wr.p/(Wr.Rtot()*2150);
 
      }	else if (Wo[0].React.reactset_flag == CH4_2STEP){ 
	Wl.v.x = 0.4101;
	Wr.v.x = 3.103;

// 	//phi = 0.6
// 	Wr.spec[0] = ZERO;       //CH4
// 	Wr.spec[1] = 0.0918;     //O2
// 	Wr.spec[2] = 0.0887;     //CO2
// 	Wr.spec[3] = 0.0759;     //H2O  
// 	Wr.spec[4] = 0.0026;     //CO
// 	Wr.rho = Wr.p/(Wr.Rtot()*1650); 

// 	//phi = 0.8
// 	Wr.spec[0] = ZERO;       //CH4
// 	Wr.spec[1] = 0.0472;     //O2
// 	Wr.spec[2] = 0.1166;     //CO2
// 	Wr.spec[3] = 0.0999;     //H2O  
// 	Wr.spec[4] = 0.0035;     //CO
// 	Wr.rho = Wr.p/(Wr.Rtot()*2000); 

	//phi = 1.0
 	Wr.spec[0] = ZERO;       //CH4
 	Wr.spec[1] = 0.005;     //O2
	Wr.spec[2] = 0.1378;     //CO2
	Wr.spec[3] = 0.1237;     //H2O 
	Wr.spec[4] = 0.0088;     //CO
	Wr.rho = Wr.p/(Wr.Rtot()*2250);

//  	//phi = 1.2
// 	Wr.spec[0] = 0.007;     //CH4
// 	Wr.spec[1] = ZERO;       //O2
// 	Wr.spec[2] = 0.1197;     //CO2
// 	Wr.spec[3] = 0.1309;     //H2O 
// 	Wr.spec[4] = 0.0256;     //CO 
// 	Wr.rho = Wr.p/(Wr.Rtot()*2260);

//   	//phi = 1.5
// 	Wr.spec[0] = 0.0232;     //CH4
// 	Wr.spec[1] = ZERO;       //O2
// 	Wr.spec[2] = 0.1188;     //CO2
// 	Wr.spec[3] = 0.1285;     //H2O 
// 	Wr.spec[4] = 0.0243;     //CO 
// 	Wr.rho = Wr.p/(Wr.Rtot()*2150);

      } else if (Wo[0].React.reactset_flag == H2O2_2STEP) {
	Wr.spec[0] = ZERO;       //H2
 	Wr.spec[1] = ZERO;       //O2
	Wr.spec[2] = 0.0049;     //OH
	Wr.spec[3] = 0.2499;     //H2O 

	Wl.v.x = 0.79625;
	Wr.rho = Wr.p/(Wr.Rtot()*3000.0); 
	Wr.v.x = 6.0;
      } else {
	cout<<"\n No 1D_Premixed Flame Initial Conditions for "<<Wo[0].React.Reaction_system; exit(1);
      }
     
      // Set Initial condtions on 1D grid
      for (int j  = SolnBlk.JCl-SolnBlk.Nghost ; j <= SolnBlk.JCu+SolnBlk.Nghost ; ++j ) {
	for ( int i = SolnBlk.ICl-SolnBlk.Nghost ; i <= SolnBlk.ICu+SolnBlk.Nghost ; ++i ) {
	  if (SolnBlk.Grid.Cell[i][j].Xc.x <= 0.01){ //spatial relation so grid independent 
	    SolnBlk.W[i][j] = Wl;  
	  } else {
 	    SolnBlk.W[i][j] = Wr;	     
 	  } /* end if */
	  SolnBlk.U[i][j] = U(SolnBlk.W[i][j]);
	} 
      } 
      break;
      
      /**************************************************************************
       The Core Flame Initial Conditions sets up a laminar diffusion 
       flame in the axisymmetric coordinate system, with the axis of symmetry
       being the left or west boundary y-axis.
       
       The condition's are based on those outlined by Mohammed et al. as well
       as Day and Bell.  
 
      **************************************************************************/
    case IC_CHEM_CORE_FLAME :           
    
      fuel_spacing = 0.002;      //m
      fuel_velocity = 0.70;      //m/s  //0.70
      fuel_temp_inlet = 298.0;   //K 
      tube_thickness = 0.00038;  //m delta
      air_spacing = 0.025;       //m 
      air_velocity = 0.35;       //m/s  0.35
      air_temp_inlet = 298.0;    //K
      ignition_temp = 1300.0;    //K
      
      Wr = Wo[0];
      Wl = Wo[0];
      
      //fuel 65% FUEL & 35% N2
      Wl.spec[0] = 0.5149;   //CH4
      Wl.spec[Wl.ns-1] = 0.4851;   //N2
      for(int q=1; q < Wl.ns-1; q++){
	Wl.spec[q].c =ZERO;
      }
      Wl.rho = Wl.p/(Wl.Rtot()*fuel_temp_inlet);
      Wl.v.y = fuel_velocity;
      Wl.v.x = ZERO;

      //air 21% O2 & 79% N2
      for(int q=0; q < Wr.ns; q++){
	if(q ==1 || q == Wr.ns-1){
	  Wr.spec[1] = 0.232;
	  Wr.spec[Wr.ns-1] = 0.768;
	} else {
	  Wr.spec[q] = ZERO;
	}
      }
      Wr.rho = Wr.p/(Wr.Rtot()*air_temp_inlet);
      Wr.v.zero();

      for (int j  = SolnBlk.JCl-SolnBlk.Nghost ; j <= SolnBlk.JCu+SolnBlk.Nghost ; ++j ) {
	for ( int i = SolnBlk.ICl-SolnBlk.Nghost ; i <= SolnBlk.ICu+SolnBlk.Nghost ; ++i ) {
  	  //region for injected fuel parabolic profile
	  if (SolnBlk.Grid.Cell[i][j].Xc.x <= fuel_spacing ){ 
	    if (SolnBlk.Grid.Cell[i][j].Xc.y <0.006 ){
	      SolnBlk.W[i][j] = Wl;
	      SolnBlk.W[i][j].rho = Wl.p/(Wl.Rtot()*fuel_temp_inlet);  //set BC value to proper inlet temp
	      SolnBlk.W[i][j].v.y = (ONE - pow((SolnBlk.Grid.Cell[i][j].Xc.x/fuel_spacing),TWO))*fuel_velocity;
	      SolnBlk.W[i][j].v.x = ZERO;
	    } else {
	      SolnBlk.W[i][j] = Wr;
	      SolnBlk.W[i][j].rho = Wr.p/(Wr.Rtot()*air_temp_inlet);
	      SolnBlk.W[i][j].v.y = (ONE - pow((SolnBlk.Grid.Cell[i][j].Xc.x/fuel_spacing),TWO))*fuel_velocity;
	      SolnBlk.W[i][j].v.x = ZERO;
	    }
	    //region for injected air
	  } else if (SolnBlk.Grid.Cell[i][j].Xc.x > fuel_spacing+tube_thickness && 
		     SolnBlk.Grid.Cell[i][j].Xc.x <= air_spacing ){		    
	    SolnBlk.W[i][j] = Wr;
	    SolnBlk.W[i][j].rho = Wr.p/(Wr.Rtot()*air_temp_inlet);
	    SolnBlk.W[i][j].v.y = air_velocity;
	    SolnBlk.W[i][j].v.x = ZERO;
	    //region for quiesent air
	  } else {
	    SolnBlk.W[i][j] = Wr;	  
	    SolnBlk.W[i][j].v.zero(); 	    
	  } 

	  //IGNITOR across fuel and air inlets
	  if( SolnBlk.Grid.Cell[i][j].Xc.y < 0.006 && SolnBlk.Grid.Cell[i][j].Xc.y > 0.003){ 
	    if ( SolnBlk.Grid.Cell[i][j].Xc.x <= fuel_spacing && SolnBlk.Grid.Cell[i][j].Xc.y <0.006 && 
	      		 SolnBlk.Grid.Cell[i][j].Xc.x > fuel_spacing*0.25){
	      SolnBlk.W[i][j].rho = Wl.p/(Wl.Rtot()*ignition_temp);
	    } else if (SolnBlk.Grid.Cell[i][j].Xc.x <= fuel_spacing){
	      SolnBlk.W[i][j].rho = Wr.p/(Wr.Rtot()*ignition_temp);
	    } else if (SolnBlk.Grid.Cell[i][j].Xc.x <= air_spacing*0.25){
	      SolnBlk.W[i][j].rho = Wr.p/(Wr.Rtot()*ignition_temp);
	    } else {
	      //left at air
	    }
	  }
	  SolnBlk.U[i][j] = U(SolnBlk.W[i][j]);
	} 
      } 
      break;

      /************************************************
       The Inverse Flame Initial Conditions sets up a
       laminar diffusion flame in the axisymmetric
       coordinate system in the same way as the core
       flame, except the oxidizer is the center annulus, 
       with the fuel coflowed in the outer annulus.
      ************************************************/
 //    case IC_CHEM_INVERSE_FLAME :	
      
//       fuel_spacing = 0.014;      //m
//       fuel_velocity = 0.1076;    //m/s 
//       fuel_temp_inlet = 298.0;   //K 
//       tube_thickness = 0.0007;  //m delta
//       air_spacing = 0.0055;       //m 
//       air_velocity = 0.315;       //m/s  
//       air_temp_inlet = 298.0;    //K
//       ignition_temp = 1200.0;    //K
      
//       Wr = Wo[0]; //oxidizer
//       Wl = Wo[0]; //fuel

//       //fuel 65% FUEL & 35% N2
//       Wl.spec[0] = 1.0; //0.5149;   //CH4
//       Wl.spec[Wl.ns-1] = 0.0; //  0.4851;   //N2
//       for(int q=1; q < Wl.ns-1; q++){
// 	Wl.spec[q].c =ZERO;
//       }
//       Wl.rho = Wl.p/(Wl.Rtot()*fuel_temp_inlet);
//       Wl.v.y = fuel_velocity;
//       Wl.v.x = ZERO;

//       //air 21% O2 & 79% N2
//       for(int q=0; q < Wr.ns; q++){
// 	if(q ==1 || q == Wr.ns-1){
// 	  Wr.spec[1] = 0.232;
// 	  Wr.spec[Wr.ns-1] = 0.768;
// 	} else {
// 	  Wr.spec[q] = ZERO;
// 	}
//       }
//       Wr.rho = Wr.p/(Wr.Rtot()*air_temp_inlet);
//       Wr.v.zero();

//       for (int j  = SolnBlk.JCl-SolnBlk.Nghost ; j <= SolnBlk.JCu+SolnBlk.Nghost ; ++j ) {
// 	for ( int i = SolnBlk.ICl-SolnBlk.Nghost ; i <= SolnBlk.ICu+SolnBlk.Nghost ; ++i ) {
//    	  //region for injected air parabolic profile
// 	  if (SolnBlk.Grid.Cell[i][j].Xc.x <= air_spacing && SolnBlk.Grid.Cell[i][j].Xc.y <0.01) { 
// 	    SolnBlk.W[i][j] = Wr;
// 	    SolnBlk.W[i][j].rho = Wr.p/(Wr.Rtot()*air_temp_inlet);  //set BC value to proper inlet temp
// 	    SolnBlk.W[i][j].v.y = (ONE - pow((SolnBlk.Grid.Cell[i][j].Xc.x/air_spacing),TWO))*air_velocity;
// 	    SolnBlk.W[i][j].v.x = ZERO;
	    
// 	    //gap
// 	  } else if ( SolnBlk.Grid.Cell[i][j].Xc.x <= air_spacing+tube_thickness && SolnBlk.Grid.Cell[i][j].Xc.y <0.01 ){ 
// 	    SolnBlk.W[i][j] = Wr;
// 	    SolnBlk.W[i][j].rho = Wr.p/(Wr.Rtot()*air_temp_inlet);  //set BC value to proper inlet temp
// 	    //SolnBlk.W[i][j] = Wl;
// 	    //SolnBlk.W[i][j].rho = Wl.p/(Wl.Rtot()*fuel_temp_inlet);
// 	    SolnBlk.W[i][j].v.zero();
	
// 	    //region for injected fuel
// 	  } else if (SolnBlk.Grid.Cell[i][j].Xc.x > air_spacing+tube_thickness && 
// 		     SolnBlk.Grid.Cell[i][j].Xc.x <= fuel_spacing
// 		     && SolnBlk.Grid.Cell[i][j].Xc.y < 0.01) {
// 	    SolnBlk.W[i][j] = Wl;
// 	    SolnBlk.W[i][j].rho = Wl.p/(Wl.Rtot()*fuel_temp_inlet);
// 	    SolnBlk.W[i][j].v.y = fuel_velocity;
// 	    SolnBlk.W[i][j].v.x = ZERO;

// 	    //region for quiesent fuel 
// 	  } else if ( SolnBlk.Grid.Cell[i][j].Xc.x > fuel_spacing ){ 
// 	    SolnBlk.W[i][j] = Wl;
// 	    SolnBlk.W[i][j].rho = Wl.p/(Wl.Rtot()*fuel_temp_inlet);
// 	    SolnBlk.W[i][j].v.zero();
// 	  } else {
// 	    SolnBlk.W[i][j] = Wl;
// 	    SolnBlk.W[i][j].rho = Wl.p/(Wl.Rtot()*fuel_temp_inlet);
// 	    SolnBlk.W[i][j].v.y = air_velocity;
// 	    SolnBlk.W[i][j].v.x = ZERO;
// 	  } 

// 	  //IGNITOR across fuel and air inlets
// 	  if( SolnBlk.Grid.Cell[i][j].Xc.y < 0.006 && SolnBlk.Grid.Cell[i][j].Xc.y > 0.003){ 
// 	    //if ( SolnBlk.Grid.Cell[i][j].Xc.x <= air_spacing && SolnBlk.Grid.Cell[i][j].Xc.y <0.004 && 
// 	    //  		 SolnBlk.Grid.Cell[i][j].Xc.x > air_spacing*0.6){
// 	    //  SolnBlk.W[i][j].rho = Wr.p/(Wr.Rtot()*ignition_temp);
// 	    //} else 
// 	    if (SolnBlk.Grid.Cell[i][j].Xc.x <= air_spacing){
// 	      SolnBlk.W[i][j].rho = Wr.p/(Wr.Rtot()*ignition_temp);
// 	    } else if (SolnBlk.Grid.Cell[i][j].Xc.x <= fuel_spacing*0.6 ){
// 	      SolnBlk.W[i][j].rho = Wl.p/(Wl.Rtot()*ignition_temp);
// 	    } else {
// 	      //left at temp
// 	    }
// 	  }
// 	  SolnBlk.U[i][j] = U(SolnBlk.W[i][j]);
// 	} 
//       } 
//       break;
    case IC_CHEM_INVERSE_FLAME :	
      
      fuel_spacing = 0.05;      //m
      fuel_velocity = 0.1076;    //m/s 
      fuel_temp_inlet = 298.0;   //K 
      tube_thickness = 0.0007;  //m delta
      air_spacing = 0.0055;       //m 
      air_velocity = 0.315;       //m/s  
      air_temp_inlet = 298.0;    //K
      ignition_temp = 1200.0;    //K
      
      Wr = Wo[0]; //oxidizer
      Wl = Wo[0]; //fuel

      //fuel 100% FUEL
      Wl.spec[0] = 1.0;   //CH4
      Wl.spec[Wl.ns-1] = 0.0; //  //N2
      for(int q=1; q < Wl.ns-1; q++){
	Wl.spec[q].c =ZERO;
      }
      Wl.rho = Wl.p/(Wl.Rtot()*fuel_temp_inlet);
      Wl.v.y = fuel_velocity;
      Wl.v.x = ZERO;

      //air 21% O2 & 79% N2
      for(int q=0; q < Wr.ns; q++){
	if(q ==1 || q == Wr.ns-1){
	  Wr.spec[1] = 0.232;
	  Wr.spec[Wr.ns-1] = 0.768;
	} else {
	  Wr.spec[q] = ZERO;
	} 
      }
      Wr.rho = Wr.p/(Wr.Rtot()*air_temp_inlet);
      Wr.v.zero();

      for (int j  = SolnBlk.JCl-SolnBlk.Nghost ; j <= SolnBlk.JCu+SolnBlk.Nghost ; ++j ) {
	for ( int i = SolnBlk.ICl-SolnBlk.Nghost ; i <= SolnBlk.ICu+SolnBlk.Nghost ; ++i ) {
   	  //region for injected air parabolic profile
	  if (SolnBlk.Grid.Cell[i][j].Xc.x <= air_spacing && SolnBlk.Grid.Cell[i][j].Xc.y <0.01) { 
	    SolnBlk.W[i][j] = Wr;
	    SolnBlk.W[i][j].rho = Wr.p/(Wr.Rtot()*air_temp_inlet);  //set BC value to proper inlet temp
	    SolnBlk.W[i][j].v.y = air_velocity; //(ONE - pow((SolnBlk.Grid.Cell[i][j].Xc.x/air_spacing),TWO))*air_velocity;
	    SolnBlk.W[i][j].v.x = ZERO;
	    
	    //gap
	  } else if ( SolnBlk.Grid.Cell[i][j].Xc.x <= air_spacing+tube_thickness && SolnBlk.Grid.Cell[i][j].Xc.y <0.005 ){ 
	    SolnBlk.W[i][j] = Wr;
	    SolnBlk.W[i][j].rho = Wr.p/(Wr.Rtot()*air_temp_inlet);  //set BC value to proper inlet temp
	    //SolnBlk.W[i][j] = Wl;
	    //SolnBlk.W[i][j].rho = Wl.p/(Wl.Rtot()*fuel_temp_inlet);
	    SolnBlk.W[i][j].v.zero();
	
	    //region for injected fuel
	  } else if (SolnBlk.Grid.Cell[i][j].Xc.x > air_spacing+tube_thickness && 
		     SolnBlk.Grid.Cell[i][j].Xc.x <= fuel_spacing
		     && SolnBlk.Grid.Cell[i][j].Xc.y < 0.01) {
	    SolnBlk.W[i][j] = Wl;
	    SolnBlk.W[i][j].rho = Wl.p/(Wl.Rtot()*fuel_temp_inlet);
	    SolnBlk.W[i][j].v.y = fuel_velocity;
	    SolnBlk.W[i][j].v.x = ZERO;

	    //region for quiesent fuel 
	  } else if ( SolnBlk.Grid.Cell[i][j].Xc.x > fuel_spacing ){ 
	    SolnBlk.W[i][j] = Wl;
	    SolnBlk.W[i][j].rho = Wl.p/(Wl.Rtot()*fuel_temp_inlet);
	    SolnBlk.W[i][j].v.zero();
	  } else {
	    SolnBlk.W[i][j] = Wl;
	    SolnBlk.W[i][j].rho = Wl.p/(Wl.Rtot()*fuel_temp_inlet);
	    SolnBlk.W[i][j].v.y = air_velocity;
	    SolnBlk.W[i][j].v.x = ZERO;
	  } 

	  //IGNITOR across fuel and air inlets
	  if( SolnBlk.Grid.Cell[i][j].Xc.y < 0.006 && SolnBlk.Grid.Cell[i][j].Xc.y > 0.003){ 
	    if (SolnBlk.Grid.Cell[i][j].Xc.x <= air_spacing){
	      SolnBlk.W[i][j].rho = Wr.p/(Wr.Rtot()*ignition_temp);
	    } else if (SolnBlk.Grid.Cell[i][j].Xc.x <= fuel_spacing*0.3 ){
	      SolnBlk.W[i][j].rho = Wl.p/(Wl.Rtot()*ignition_temp);
	    } else {
	      //left at temp
	    }
	  }
	  SolnBlk.U[i][j] = U(SolnBlk.W[i][j]);
	} 
      } 
      break;

      /******************** DEFAULT *******************/
    default:
      // Set the solution state to the initial state Wo[0].
      for (int j  = SolnBlk.JCl-SolnBlk.Nghost ; j <= SolnBlk.JCu+SolnBlk.Nghost ; ++j ) {
	for ( int i = SolnBlk.ICl-SolnBlk.Nghost ; i <= SolnBlk.ICu+SolnBlk.Nghost ; ++i ) {
	  SolnBlk.W[i][j] = Wo[0];
	  SolnBlk.U[i][j] = U(SolnBlk.W[i][j]);
	} /* endfor */
      } /* endfor */
      break;
    } /* endswitch */

    /***********************************************************************
     ***********************************************************************
     ***********************************************************************/
 
    /* Set the solution residuals, gradients, limiters, and 
       other values to zero. */
    
    for (int j  = SolnBlk.JCl-SolnBlk.Nghost ; j <= SolnBlk.JCu+SolnBlk.Nghost ; ++j ) {
      for ( int i = SolnBlk.ICl-SolnBlk.Nghost ; i <= SolnBlk.ICu+SolnBlk.Nghost ; ++i ) {
	for ( int k = 0 ; k <= NUMBER_OF_RESIDUAL_VECTORS_CHEM2D-1 ; ++k ) {
	  SolnBlk.dUdt[i][j][k].Vacuum();
	} /* endfor */
	SolnBlk.dWdx[i][j].Vacuum();
	SolnBlk.dWdy[i][j].Vacuum();
	SolnBlk.phi[i][j].Vacuum();
	SolnBlk.Uo[i][j].Vacuum();
	SolnBlk.dt[i][j] = ZERO;
      } /* endfor */
    } /* endfor */
    
    /* Set default values for the boundary conditions
       reference states. */

    for (int j = SolnBlk.JCl-SolnBlk.Nghost ; j <= SolnBlk.JCu+SolnBlk.Nghost ; ++j ) {
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
    
    for ( int i = SolnBlk.ICl-SolnBlk.Nghost ; i <= SolnBlk.ICu+SolnBlk.Nghost ; ++i ) {
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

/********************************************************
 * Routine: BCs                                         *
 *                                                      *
 * Apply boundary conditions at boundaries of the       *
 * specified quadrilateral solution block.              *
 *                                                      *
 ********************************************************/
void BCs(Chem2D_Quad_Block &SolnBlk,
	 Chem2D_Input_Parameters &Input_Parameters) {

  int i, j;
  Vector2D dX;
  Chem2D_pState dW, W, W_wall_src;
  double dp = 3177.7; //needs to be input parameter dp/dx or dp/1m

  W_wall_src = Chem2D_pState(DENSITY_STDATM, 
                             Vector2D_ZERO, 
                             DENSITY_STDATM*W_wall_src.Rtot()*TEMPERATURE_STDATM, 
                             W_wall_src.k,  
                             W_wall_src.omega );

  //WEST
  for ( j = SolnBlk.JCl-SolnBlk.Nghost ; j <= SolnBlk.JCu+SolnBlk.Nghost ; ++j ) { 
      if ( (j >= SolnBlk.JCl && j <= SolnBlk.JCu) ||  
	   // These are conditions for the corner ghostcells, its confusing....
           (j < SolnBlk.JCl && 
            (SolnBlk.Grid.BCtypeS[SolnBlk.ICl-1] == BC_NONE ||
             SolnBlk.Grid.BCtypeS[SolnBlk.ICl-1] == BC_PERIODIC ||
             SolnBlk.Grid.BCtypeS[SolnBlk.ICl-1] == BC_CONSTANT_EXTRAPOLATION ||
             SolnBlk.Grid.BCtypeS[SolnBlk.ICl-1] == BC_LINEAR_EXTRAPOLATION ||
             SolnBlk.Grid.BCtypeS[SolnBlk.ICl-1] == BC_CHARACTERISTIC ||
	     SolnBlk.Grid.BCtypeS[SolnBlk.ICl-1] == BC_INFLOW_SUBSONIC ||
	     SolnBlk.Grid.BCtypeS[SolnBlk.ICl-1] == BC_OUTFLOW_SUBSONIC || 
	     SolnBlk.Grid.BCtypeS[SolnBlk.ICl-1] == BC_FLAME_INFLOW ||
	     SolnBlk.Grid.BCtypeS[SolnBlk.ICl-1] == BC_FLAME_OUTFLOW) ) ||
           (j > SolnBlk.JCu && 
            (SolnBlk.Grid.BCtypeN[SolnBlk.ICl-1] == BC_NONE ||
             SolnBlk.Grid.BCtypeN[SolnBlk.ICl-1] == BC_PERIODIC ||
             SolnBlk.Grid.BCtypeN[SolnBlk.ICl-1] == BC_CONSTANT_EXTRAPOLATION ||
             SolnBlk.Grid.BCtypeN[SolnBlk.ICl-1] == BC_LINEAR_EXTRAPOLATION ||
             SolnBlk.Grid.BCtypeN[SolnBlk.ICl-1] == BC_CHARACTERISTIC ||
	     SolnBlk.Grid.BCtypeN[SolnBlk.ICl-1] == BC_INFLOW_SUBSONIC ||
	     SolnBlk.Grid.BCtypeN[SolnBlk.ICl-1] == BC_OUTFLOW_SUBSONIC || 
	     SolnBlk.Grid.BCtypeN[SolnBlk.ICl-1] == BC_FLAME_INFLOW ||
	     SolnBlk.Grid.BCtypeN[SolnBlk.ICl-1] == BC_FLAME_OUTFLOW) ) ) {
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
          case BC_LINEAR_EXTRAPOLATION :
            Linear_Reconstruction_LeastSquares_2(SolnBlk, 
                                                 SolnBlk.ICl, j, 
                                                 LIMITER_BARTH_JESPERSEN);
            dX = SolnBlk.Grid.Cell[SolnBlk.ICl-1][j].Xc -
                 SolnBlk.Grid.Cell[SolnBlk.ICl][j].Xc;
            SolnBlk.W[SolnBlk.ICl-1][j] = SolnBlk.W[SolnBlk.ICl][j] + 
               (SolnBlk.phi[SolnBlk.ICl][j]^SolnBlk.dWdx[SolnBlk.ICl][j])*dX.x +
               (SolnBlk.phi[SolnBlk.ICl][j]^SolnBlk.dWdy[SolnBlk.ICl][j])*dX.y;
            SolnBlk.U[SolnBlk.ICl-1][j] = U(SolnBlk.W[SolnBlk.ICl-1][j]);
            dX = SolnBlk.Grid.Cell[SolnBlk.ICl-2][j].Xc -
                 SolnBlk.Grid.Cell[SolnBlk.ICl][j].Xc;
            SolnBlk.W[SolnBlk.ICl-2][j] = SolnBlk.W[SolnBlk.ICl][j] + 
               (SolnBlk.phi[SolnBlk.ICl][j]^SolnBlk.dWdx[SolnBlk.ICl][j])*dX.x +
               (SolnBlk.phi[SolnBlk.ICl][j]^SolnBlk.dWdy[SolnBlk.ICl][j])*dX.y;
            SolnBlk.U[SolnBlk.ICl-2][j] = U(SolnBlk.W[SolnBlk.ICl-2][j]);
            break;
          case BC_REFLECTION :
            SolnBlk.W[SolnBlk.ICl-1][j] = Reflect(SolnBlk.W[SolnBlk.ICl][j],
	    					  SolnBlk.Grid.nfaceW(SolnBlk.ICl,j));
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
            SolnBlk.W[SolnBlk.ICl-1][j] = 
               BC_Characteristic_Pressure(SolnBlk.W[SolnBlk.ICl][j],
                                          SolnBlk.WoW[j], 
                                          SolnBlk.Grid.nfaceW(SolnBlk.ICl,j));
            SolnBlk.U[SolnBlk.ICl-1][j] = U(SolnBlk.W[SolnBlk.ICl-1][j]);
            SolnBlk.W[SolnBlk.ICl-2][j] = SolnBlk.W[SolnBlk.ICl-1][j];
            SolnBlk.U[SolnBlk.ICl-2][j] = U(SolnBlk.W[SolnBlk.ICl-2][j]);
            break;

	    /* Currently set to Fixed Wall Temperature */
      case BC_FREE_SLIP :
	SolnBlk.W[SolnBlk.ICl-1][j] = Free_Slip(SolnBlk.W[SolnBlk.ICl][j],
						SolnBlk.WoW[j],
						SolnBlk.Grid.nfaceW(SolnBlk.ICl,j),
						FIXED_TEMPERATURE_WALL);
	SolnBlk.U[SolnBlk.ICl-1][j] = U(SolnBlk.W[SolnBlk.ICl-1][j]);
	SolnBlk.W[SolnBlk.ICl-2][j] = Free_Slip(SolnBlk.W[SolnBlk.ICl+1][j],
						SolnBlk.WoW[j],
						SolnBlk.Grid.nfaceW(SolnBlk.ICl,j),
						FIXED_TEMPERATURE_WALL);
	SolnBlk.U[SolnBlk.ICl-2][j] = U(SolnBlk.W[SolnBlk.ICl-2][j]);
	break;
      case BC_NO_SLIP :
	SolnBlk.W[SolnBlk.ICl-1][j] = No_Slip(SolnBlk.W[SolnBlk.ICl][j],
					      SolnBlk.WoW[j],
					      SolnBlk.Grid.nfaceW(SolnBlk.ICl,j),
					      FIXED_TEMPERATURE_WALL);
	SolnBlk.U[SolnBlk.ICl-1][j] = U(SolnBlk.W[SolnBlk.ICl-1][j]);
	SolnBlk.W[SolnBlk.ICl-2][j] = No_Slip(SolnBlk.W[SolnBlk.ICl+1][j],
					      SolnBlk.WoW[j],
					      SolnBlk.Grid.nfaceW(SolnBlk.ICl,j),
					      FIXED_TEMPERATURE_WALL);
	SolnBlk.U[SolnBlk.ICl-2][j] = U(SolnBlk.W[SolnBlk.ICl-2][j]);
	break;
      case BC_MOVING_WALL :	 
	SolnBlk.W[SolnBlk.ICl-1][j] = Moving_Wall(SolnBlk.W[SolnBlk.ICl][j],
						  SolnBlk.WoW[j],
						  SolnBlk.Grid.nfaceW(SolnBlk.ICl,j),
						  SolnBlk.Moving_wall_velocity,
						  FIXED_TEMPERATURE_WALL);
	SolnBlk.U[SolnBlk.ICl-1][j] = U(SolnBlk.W[SolnBlk.ICl-1][j]);
	SolnBlk.W[SolnBlk.ICl-2][j] = Moving_Wall(SolnBlk.W[SolnBlk.ICl+1][j],
						  SolnBlk.WoW[j],
						  SolnBlk.Grid.nfaceW(SolnBlk.ICl,j),
						  SolnBlk.Moving_wall_velocity,
						  FIXED_TEMPERATURE_WALL);
	SolnBlk.U[SolnBlk.ICl-2][j] = U(SolnBlk.W[SolnBlk.ICl-2][j]);
	break;
      case BC_INFLOW_SUBSONIC :
	// v.x (u) is constant extrapolation, p is linearly extrapolated, 
	// everything else is fixed

	dX = SolnBlk.Grid.Cell[SolnBlk.ICl-1][j].Xc - SolnBlk.Grid.Cell[SolnBlk.ICl][j].Xc; 

	SolnBlk.W[SolnBlk.ICl-1][j] = SolnBlk.WoW[j];
	SolnBlk.W[SolnBlk.ICl-1][j].v.x = SolnBlk.W[SolnBlk.ICl][j].v.x;
 	SolnBlk.W[SolnBlk.ICl-1][j].p = SolnBlk.WoW[j].p - dp*dX.x;

	SolnBlk.U[SolnBlk.ICl-1][j] = U(SolnBlk.W[SolnBlk.ICl-1][j]);

	dX = SolnBlk.Grid.Cell[SolnBlk.ICl-2][j].Xc - SolnBlk.Grid.Cell[SolnBlk.ICl-1][j].Xc;

	SolnBlk.W[SolnBlk.ICl-2][j] = SolnBlk.WoW[j];
	SolnBlk.W[SolnBlk.ICl-2][j].v.x = SolnBlk.W[SolnBlk.ICl][j].v.x;
 	SolnBlk.W[SolnBlk.ICl-2][j].p = SolnBlk.W[SolnBlk.ICl-1][j].p - dp*dX.x;
	
	SolnBlk.U[SolnBlk.ICl-2][j] = U(SolnBlk.W[SolnBlk.ICl-2][j]);
	break;
      case BC_OUTFLOW_SUBSONIC :
	// all constant extrapolation except pressure which is fixed.
	SolnBlk.W[SolnBlk.ICl-1][j] = SolnBlk.W[SolnBlk.ICl][j]; 
	SolnBlk.W[SolnBlk.ICl-1][j].p = SolnBlk.WoW[j].p;
	SolnBlk.U[SolnBlk.ICl-1][j] = U(SolnBlk.W[SolnBlk.ICl-1][j]);
	
	SolnBlk.W[SolnBlk.ICl-2][j] = SolnBlk.W[SolnBlk.ICl][j];
	SolnBlk.W[SolnBlk.ICl-2][j].p =  SolnBlk.WoW[j].p;
	SolnBlk.U[SolnBlk.ICl-2][j] = U(SolnBlk.W[SolnBlk.ICl-2][j]);
	break;
	
      case BC_FLAME_INFLOW :
	//fixed except for velocity??? 
	SolnBlk.W[SolnBlk.ICl-1][j] = 
	  BC_Flame_Inflow(SolnBlk.W[SolnBlk.ICl][j],
			  SolnBlk.WoW[j], 
			  SolnBlk.W[SolnBlk.ICu][j],
			  SolnBlk.Grid.nfaceW(SolnBlk.ICl,j)); 
	SolnBlk.U[SolnBlk.ICl-1][j] = U(SolnBlk.W[SolnBlk.ICl-1][j]);
	SolnBlk.W[SolnBlk.ICl-2][j] = SolnBlk.W[SolnBlk.ICl-1][j];
	SolnBlk.U[SolnBlk.ICl-2][j] = U(SolnBlk.W[SolnBlk.ICl-2][j]);
	break;

	case BC_FLAME_OUTFLOW :
	  //isentropic condition for velocity
	  SolnBlk.W[SolnBlk.ICl-1][j] = 
	    BC_Flame_Outflow(SolnBlk.W[SolnBlk.ICl][j],
			     SolnBlk.WoW[j],
			     SolnBlk.W[SolnBlk.ICu][j],
			     SolnBlk.Grid.nfaceW(SolnBlk.ICl,j));
	  SolnBlk.U[SolnBlk.ICl-1][j] = U(SolnBlk.W[SolnBlk.ICl-1][j]);
	  SolnBlk.W[SolnBlk.ICl-2][j] = SolnBlk.W[SolnBlk.ICl-1][j];
	  SolnBlk.U[SolnBlk.ICl-2][j] = U(SolnBlk.W[SolnBlk.ICl-2][j]);
	  break;

	default:
	    SolnBlk.W[SolnBlk.ICl-1][j] = SolnBlk.W[SolnBlk.ICl][j];
            SolnBlk.U[SolnBlk.ICl-1][j] = SolnBlk.U[SolnBlk.ICl][j];
            SolnBlk.W[SolnBlk.ICl-2][j] = SolnBlk.W[SolnBlk.ICl][j];
            SolnBlk.U[SolnBlk.ICl-2][j] = SolnBlk.U[SolnBlk.ICl][j];
            break;
        } /* endswitch WEST*/
      } /* endif WEST*/     

    //EAST
    if ( (j >= SolnBlk.JCl && j <= SolnBlk.JCu) ||
	 (j < SolnBlk.JCl && 
	  (SolnBlk.Grid.BCtypeS[SolnBlk.ICu+1] == BC_NONE ||
	   SolnBlk.Grid.BCtypeS[SolnBlk.ICu+1] == BC_PERIODIC ||
	   SolnBlk.Grid.BCtypeS[SolnBlk.ICu+1] == BC_CONSTANT_EXTRAPOLATION ||
	   SolnBlk.Grid.BCtypeS[SolnBlk.ICu+1] == BC_FIXED_TEMP_WALL ||
	   SolnBlk.Grid.BCtypeS[SolnBlk.ICu+1] == BC_ADIABATIC_WALL ||
	   SolnBlk.Grid.BCtypeS[SolnBlk.ICu+1] == BC_LINEAR_EXTRAPOLATION ||
	   SolnBlk.Grid.BCtypeS[SolnBlk.ICu+1] == BC_CHARACTERISTIC ||
	   SolnBlk.Grid.BCtypeS[SolnBlk.ICu+1] == BC_INFLOW_SUBSONIC ||
	   SolnBlk.Grid.BCtypeS[SolnBlk.ICu+1] == BC_OUTFLOW_SUBSONIC || 
	   SolnBlk.Grid.BCtypeS[SolnBlk.ICu+1] == BC_FLAME_INFLOW ||
	   SolnBlk.Grid.BCtypeS[SolnBlk.ICu+1] == BC_FLAME_OUTFLOW) ) ||
	 (j > SolnBlk.JCu && 
            (SolnBlk.Grid.BCtypeN[SolnBlk.ICu+1] == BC_NONE ||
             SolnBlk.Grid.BCtypeN[SolnBlk.ICu+1] == BC_PERIODIC ||
             SolnBlk.Grid.BCtypeN[SolnBlk.ICu+1] == BC_CONSTANT_EXTRAPOLATION ||
             SolnBlk.Grid.BCtypeN[SolnBlk.ICu+1] == BC_LINEAR_EXTRAPOLATION ||
             SolnBlk.Grid.BCtypeN[SolnBlk.ICu+1] == BC_CHARACTERISTIC ||
	     SolnBlk.Grid.BCtypeN[SolnBlk.ICu+1] == BC_INFLOW_SUBSONIC ||
	     SolnBlk.Grid.BCtypeN[SolnBlk.ICu+1] == BC_OUTFLOW_SUBSONIC || 
	     SolnBlk.Grid.BCtypeN[SolnBlk.ICu+1] == BC_FLAME_INFLOW ||
	     SolnBlk.Grid.BCtypeN[SolnBlk.ICu+1] == BC_FLAME_OUTFLOW) ) ) { 

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
	case BC_FIXED_TEMP_WALL :
	case BC_ADIABATIC_WALL :
	  if (SolnBlk.W[SolnBlk.ICu][j].v*SolnBlk.Grid.nfaceE(SolnBlk.ICu, j) >= ZERO) {
	    SolnBlk.W[SolnBlk.ICu+1][j] = SolnBlk.W[SolnBlk.ICu][j];
	    SolnBlk.U[SolnBlk.ICu+1][j] = SolnBlk.U[SolnBlk.ICu][j];
	    SolnBlk.W[SolnBlk.ICu+2][j] = SolnBlk.W[SolnBlk.ICu][j];
	    SolnBlk.U[SolnBlk.ICu+2][j] = SolnBlk.U[SolnBlk.ICu][j];
	  } else {
	    SolnBlk.W[SolnBlk.ICu+1][j] = W_wall_src;
	    //SolnBlk.W[SolnBlk.ICu+1][j] = SolnBlk.WoE[j];
	    SolnBlk.U[SolnBlk.ICu+1][j] = U(SolnBlk.W[SolnBlk.ICu+1][j]);
	    SolnBlk.W[SolnBlk.ICu+2][j] = W_wall_src;
	    //SolnBlk.W[SolnBlk.ICu+2][j] = SolnBlk.WoE[j];
	    SolnBlk.U[SolnBlk.ICu+2][j] = U(SolnBlk.W[SolnBlk.ICu+2][j]);
	  } /* endif */
	  break;
	case BC_LINEAR_EXTRAPOLATION :
	  Linear_Reconstruction_LeastSquares_2(SolnBlk, 
					       SolnBlk.ICu, j, 
                                                 LIMITER_BARTH_JESPERSEN);
            dX = SolnBlk.Grid.Cell[SolnBlk.ICu+1][j].Xc -
                 SolnBlk.Grid.Cell[SolnBlk.ICu][j].Xc;
            SolnBlk.W[SolnBlk.ICu+1][j] = SolnBlk.W[SolnBlk.ICu][j] + 
               (SolnBlk.phi[SolnBlk.ICu][j]^SolnBlk.dWdx[SolnBlk.ICu][j])*dX.x +
               (SolnBlk.phi[SolnBlk.ICu][j]^SolnBlk.dWdy[SolnBlk.ICu][j])*dX.y;
            SolnBlk.U[SolnBlk.ICu+1][j] = U(SolnBlk.W[SolnBlk.ICu+1][j]);
            dX = SolnBlk.Grid.Cell[SolnBlk.ICu+2][j].Xc -
                 SolnBlk.Grid.Cell[SolnBlk.ICu][j].Xc;
            SolnBlk.W[SolnBlk.ICu+2][j] = SolnBlk.W[SolnBlk.ICu][j] + 
               (SolnBlk.phi[SolnBlk.ICu][j]^SolnBlk.dWdx[SolnBlk.ICu][j])*dX.x +
               (SolnBlk.phi[SolnBlk.ICu][j]^SolnBlk.dWdy[SolnBlk.ICu][j])*dX.y;
            SolnBlk.U[SolnBlk.ICu+2][j] = U(SolnBlk.W[SolnBlk.ICu+2][j]);
            break;
           case BC_REFLECTION :
	     //cout<<" "<<j;
	    SolnBlk.W[SolnBlk.ICu+1][j] = Reflect(SolnBlk.W[SolnBlk.ICu][j],
				       SolnBlk.Grid.nfaceE(SolnBlk.ICu,j));
            SolnBlk.U[SolnBlk.ICu+1][j] = U(SolnBlk.W[SolnBlk.ICu+1][j]);
            SolnBlk.W[SolnBlk.ICu+2][j] = Reflect(SolnBlk.W[SolnBlk.ICu-1][j],
	                                SolnBlk.Grid.nfaceE(SolnBlk.ICu,j));
            SolnBlk.U[SolnBlk.ICu+2][j] = U(SolnBlk.W[SolnBlk.ICu+2][j]);
	    break;
          case BC_PERIODIC :
            SolnBlk.W[SolnBlk.ICu+1][j] = SolnBlk.W[SolnBlk.ICl+1][j];
            SolnBlk.U[SolnBlk.ICu+1][j] = SolnBlk.U[SolnBlk.ICl+1][j];
            SolnBlk.W[SolnBlk.ICu+2][j] = SolnBlk.W[SolnBlk.ICl+2][j];
            SolnBlk.U[SolnBlk.ICu+2][j] = SolnBlk.U[SolnBlk.ICl+2][j];
            break;
          case BC_CHARACTERISTIC :
            SolnBlk.W[SolnBlk.ICu+1][j] = 
               BC_Characteristic_Pressure(SolnBlk.W[SolnBlk.ICu][j],
                                          SolnBlk.WoE[j],
                                          SolnBlk.Grid.nfaceE(SolnBlk.ICu,j));
            SolnBlk.U[SolnBlk.ICu+1][j] = U(SolnBlk.W[SolnBlk.ICu+1][j]);
            SolnBlk.W[SolnBlk.ICu+2][j] = SolnBlk.W[SolnBlk.ICu+1][j];
            SolnBlk.U[SolnBlk.ICu+2][j] = U(SolnBlk.W[SolnBlk.ICu+2][j]);
            break;
	case BC_FREE_SLIP :
	  SolnBlk.W[SolnBlk.ICu+1][j] = Free_Slip(SolnBlk.W[SolnBlk.ICu][j],
						SolnBlk.WoE[j],
						SolnBlk.Grid.nfaceE(SolnBlk.ICu,j),
						FIXED_TEMPERATURE_WALL);
	  SolnBlk.U[SolnBlk.ICu+1][j] = U(SolnBlk.W[SolnBlk.ICu+1][j]);
	  SolnBlk.W[SolnBlk.ICu+2][j] = Free_Slip(SolnBlk.W[SolnBlk.ICu-1][j],
						SolnBlk.WoE[j],
						SolnBlk.Grid.nfaceE(SolnBlk.ICu,j),
						FIXED_TEMPERATURE_WALL);
	  SolnBlk.U[SolnBlk.ICu+2][j] = U(SolnBlk.W[SolnBlk.ICu+2][j]);
	  break;
	case BC_NO_SLIP :
	  SolnBlk.W[SolnBlk.ICu+1][j] = No_Slip(SolnBlk.W[SolnBlk.ICu][j],
						SolnBlk.WoE[j],
						SolnBlk.Grid.nfaceE(SolnBlk.ICu,j),
						FIXED_TEMPERATURE_WALL);
	  SolnBlk.U[SolnBlk.ICu+1][j] = U(SolnBlk.W[SolnBlk.ICu+1][j]);
	  SolnBlk.W[SolnBlk.ICu+2][j] = No_Slip(SolnBlk.W[SolnBlk.ICu-1][j],
						SolnBlk.WoE[j],
						SolnBlk.Grid.nfaceE(SolnBlk.ICu,j),
						FIXED_TEMPERATURE_WALL);
	  SolnBlk.U[SolnBlk.ICu+2][j] = U(SolnBlk.W[SolnBlk.ICu+2][j]);
	  break;
	case BC_MOVING_WALL :
	  SolnBlk.W[SolnBlk.ICu+1][j] = Moving_Wall(SolnBlk.W[SolnBlk.ICu][j],
						    SolnBlk.WoE[j],
						    SolnBlk.Grid.nfaceE(SolnBlk.ICu,j),
						    SolnBlk.Moving_wall_velocity,
						    FIXED_TEMPERATURE_WALL);
	  SolnBlk.U[SolnBlk.ICu+1][j] = U(SolnBlk.W[SolnBlk.ICu+1][j]);
	  SolnBlk.W[SolnBlk.ICu+2][j] = Moving_Wall(SolnBlk.W[SolnBlk.ICu-1][j],
						    SolnBlk.WoE[j],
						    SolnBlk.Grid.nfaceE(SolnBlk.ICu,j),
						    SolnBlk.Moving_wall_velocity,
						    FIXED_TEMPERATURE_WALL);
	  SolnBlk.U[SolnBlk.ICu+2][j] = U(SolnBlk.W[SolnBlk.ICu+2][j]);
	  break;	
	case BC_INFLOW_SUBSONIC :
	  // all fixed except v.x (u) which is constant extrapolation	  
	  SolnBlk.W[SolnBlk.ICu+1][j] = SolnBlk.WoE[j];
	  SolnBlk.W[SolnBlk.ICu+1][j].v.x = SolnBlk.W[SolnBlk.ICu][j].v.x;
	  SolnBlk.U[SolnBlk.ICu+1][j] = U(SolnBlk.W[SolnBlk.ICu+1][j]);

	  SolnBlk.W[SolnBlk.ICu+2][j] = SolnBlk.WoE[j];
	  SolnBlk.W[SolnBlk.ICu+2][j].v.x = SolnBlk.W[SolnBlk.ICu][j].v.x;
	  SolnBlk.U[SolnBlk.ICu+2][j] = U(SolnBlk.W[SolnBlk.ICu+2][j]);
	  break;
	case BC_OUTFLOW_SUBSONIC :
	  // all constant extrapolation except pressure which linearly extrapolated

	  dX = SolnBlk.Grid.Cell[SolnBlk.ICu+1][j].Xc - SolnBlk.Grid.Cell[SolnBlk.ICu][j].Xc; 
	  SolnBlk.W[SolnBlk.ICu+1][j] = SolnBlk.W[SolnBlk.ICu][j]; 
	  SolnBlk.W[SolnBlk.ICu+1][j].p = SolnBlk.WoE[j].p - dp*dX.x; 	  
	  SolnBlk.U[SolnBlk.ICu+1][j] = U(SolnBlk.W[SolnBlk.ICu+1][j]);

	  dX = SolnBlk.Grid.Cell[SolnBlk.ICu+2][j].Xc - SolnBlk.Grid.Cell[SolnBlk.ICu+1][j].Xc; 
	  SolnBlk.W[SolnBlk.ICu+2][j] = SolnBlk.W[SolnBlk.ICu][j];
	  SolnBlk.W[SolnBlk.ICu+2][j].p = SolnBlk.W[SolnBlk.ICu+1][j].p - dp*dX.x; 	    
	  SolnBlk.U[SolnBlk.ICu+2][j] = U(SolnBlk.W[SolnBlk.ICu+2][j]);
	  break;

	case BC_FLAME_INFLOW :	
	  SolnBlk.W[SolnBlk.ICu+1][j] = 
	    BC_Flame_Inflow(SolnBlk.W[SolnBlk.ICu][j],
			    SolnBlk.WoE[j], 
			    SolnBlk.W[SolnBlk.ICl][j],
			    SolnBlk.Grid.nfaceE(SolnBlk.ICu,j)); 

	  SolnBlk.U[SolnBlk.ICu+1][j] = U(SolnBlk.W[SolnBlk.ICu+1][j]);
	  SolnBlk.W[SolnBlk.ICu+2][j] = SolnBlk.W[SolnBlk.ICu+1][j];
	  SolnBlk.U[SolnBlk.ICu+2][j] = U(SolnBlk.W[SolnBlk.ICu+2][j]);	
	  break;

	case BC_FLAME_OUTFLOW : 
	  //isentropic condition for velocity ??
	  SolnBlk.W[SolnBlk.ICu+1][j] = 
	    BC_Flame_Outflow(SolnBlk.W[SolnBlk.ICu][j],
			     SolnBlk.WoE[j],
			     SolnBlk.W[SolnBlk.ICl][j],
			     SolnBlk.Grid.nfaceE(SolnBlk.ICu,j));
	  SolnBlk.U[SolnBlk.ICu+1][j] = U(SolnBlk.W[SolnBlk.ICu+1][j]);
	  SolnBlk.W[SolnBlk.ICu+2][j] = SolnBlk.W[SolnBlk.ICu+1][j];
	  SolnBlk.U[SolnBlk.ICu+2][j] = U(SolnBlk.W[SolnBlk.ICu+2][j]);
	  break;

	default:
            SolnBlk.W[SolnBlk.ICu+1][j] = SolnBlk.W[SolnBlk.ICu][j];
            SolnBlk.U[SolnBlk.ICu+1][j] = SolnBlk.U[SolnBlk.ICu][j];
            SolnBlk.W[SolnBlk.ICu+2][j] = SolnBlk.W[SolnBlk.ICu][j];
            SolnBlk.U[SolnBlk.ICu+2][j] = SolnBlk.U[SolnBlk.ICu][j];
            break;
        } /* endswitch EAST */
      } /* endif EAST*/
    } /*endfor EAST & WEST */
   

    //SOUTH
    for ( i = SolnBlk.ICl-SolnBlk.Nghost ; i <= SolnBlk.ICu+SolnBlk.Nghost ; ++i ) {
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
      case BC_LINEAR_EXTRAPOLATION :
	Linear_Reconstruction_LeastSquares_2(SolnBlk, 
					     i, SolnBlk.JCl, 
					     LIMITER_BARTH_JESPERSEN);
	dX = SolnBlk.Grid.Cell[i][SolnBlk.JCl-1].Xc -
	  SolnBlk.Grid.Cell[i][SolnBlk.JCl].Xc;
	SolnBlk.W[i][SolnBlk.JCl-1] = SolnBlk.W[i][SolnBlk.JCl] + 
	  (SolnBlk.phi[i][SolnBlk.JCl]^SolnBlk.dWdx[i][SolnBlk.JCl])*dX.x +
	  (SolnBlk.phi[i][SolnBlk.JCl]^SolnBlk.dWdy[i][SolnBlk.JCl])*dX.y;
	SolnBlk.U[i][SolnBlk.JCl-1] = U(SolnBlk.W[i][SolnBlk.JCl-1]);
	dX = SolnBlk.Grid.Cell[i][SolnBlk.JCl-2].Xc -
	  SolnBlk.Grid.Cell[i][SolnBlk.JCl].Xc;
	SolnBlk.W[i][SolnBlk.JCl-2] = SolnBlk.W[i][SolnBlk.JCl] + 
	  (SolnBlk.phi[i][SolnBlk.JCl]^SolnBlk.dWdx[i][SolnBlk.JCl])*dX.x +
	  (SolnBlk.phi[i][SolnBlk.JCl]^SolnBlk.dWdy[i][SolnBlk.JCl])*dX.y;
	SolnBlk.U[i][SolnBlk.JCl-2] = U(SolnBlk.W[i][SolnBlk.JCl-2]);
	break;
      case BC_REFLECTION :  
	//cout<<" "<<i;
	SolnBlk.W[i][SolnBlk.JCl-1] = Reflect(SolnBlk.W[i][SolnBlk.JCl],
					      SolnBlk.Grid.nfaceS(i, SolnBlk.JCl));
	SolnBlk.U[i][SolnBlk.JCl-1] = U(SolnBlk.W[i][SolnBlk.JCl-1]);
	SolnBlk.W[i][SolnBlk.JCl-2] = Reflect(SolnBlk.W[i][SolnBlk.JCl+1],
					      SolnBlk.Grid.nfaceS(i, SolnBlk.JCl));
	SolnBlk.U[i][SolnBlk.JCl-2] = U(SolnBlk.W[i][SolnBlk.JCl-2]);
	break;
      case BC_PERIODIC :
	SolnBlk.W[i][SolnBlk.JCl-1] = SolnBlk.W[i][SolnBlk.JCu-1];
	SolnBlk.U[i][SolnBlk.JCl-1] = SolnBlk.U[i][SolnBlk.JCu-1];
	SolnBlk.W[i][SolnBlk.JCl-2] = SolnBlk.W[i][SolnBlk.JCu-2];
	SolnBlk.U[i][SolnBlk.JCl-2] = SolnBlk.U[i][SolnBlk.JCu-2];
	break;
      case BC_CHARACTERISTIC :
	SolnBlk.W[i][SolnBlk.JCl-1] = 
	  BC_Characteristic_Pressure(SolnBlk.W[i][SolnBlk.JCl],
				     SolnBlk.WoS[i],
				     SolnBlk.Grid.nfaceS(i, SolnBlk.JCl));
	SolnBlk.U[i][SolnBlk.JCl-1] = U(SolnBlk.W[i][SolnBlk.JCl-1]);
	SolnBlk.W[i][SolnBlk.JCl-2] = SolnBlk.W[i][SolnBlk.JCl-1];
	SolnBlk.U[i][SolnBlk.JCl-2] = U(SolnBlk.W[i][SolnBlk.JCl-2]);
	break;
      case BC_FREE_SLIP :
	SolnBlk.W[i][SolnBlk.JCl-1] = Free_Slip(SolnBlk.W[i][SolnBlk.JCl],
						SolnBlk.WoS[i],
						SolnBlk.Grid.nfaceS(i, SolnBlk.JCl),
						FIXED_TEMPERATURE_WALL);
	SolnBlk.U[i][SolnBlk.JCl-1] = U(SolnBlk.W[i][SolnBlk.JCl-1]);
	SolnBlk.W[i][SolnBlk.JCl-2] = Free_Slip(SolnBlk.W[i][SolnBlk.JCl+1],
						SolnBlk.WoS[i],
						SolnBlk.Grid.nfaceS(i, SolnBlk.JCl),
						FIXED_TEMPERATURE_WALL);
	SolnBlk.U[i][SolnBlk.JCl-2] = U(SolnBlk.W[i][SolnBlk.JCl-2]);
	break;
      case BC_NO_SLIP :
	
	//		if (i >= SolnBlk.ICl && j <= SolnBlk.ICu){
	SolnBlk.W[i][SolnBlk.JCl-1] = No_Slip(SolnBlk.W[i][SolnBlk.JCl],
					      SolnBlk.WoS[i],
					      SolnBlk.Grid.nfaceS(i, SolnBlk.JCl),
					      FIXED_TEMPERATURE_WALL);
	SolnBlk.U[i][SolnBlk.JCl-1] = U(SolnBlk.W[i][SolnBlk.JCl-1]);
	SolnBlk.W[i][SolnBlk.JCl-2] = No_Slip(SolnBlk.W[i][SolnBlk.JCl+1],
					      SolnBlk.WoS[i],
					      SolnBlk.Grid.nfaceS(i, SolnBlk.JCl),
					      FIXED_TEMPERATURE_WALL);
	SolnBlk.U[i][SolnBlk.JCl-2] = U(SolnBlk.W[i][SolnBlk.JCl-2]);

	//		}

     	break;
      case BC_MOVING_WALL :
	SolnBlk.W[i][SolnBlk.JCl-1] = Moving_Wall(SolnBlk.W[i][SolnBlk.JCl], 
						  SolnBlk.WoS[i],
						  SolnBlk.Grid.nfaceS(i, SolnBlk.JCl),
						  SolnBlk.Moving_wall_velocity,
						  FIXED_TEMPERATURE_WALL);
	SolnBlk.U[i][SolnBlk.JCl-1] = U(SolnBlk.W[i][SolnBlk.JCl-1]);
	SolnBlk.W[i][SolnBlk.JCl-2] = Moving_Wall(SolnBlk.W[i][SolnBlk.JCl+1], 
						  SolnBlk.WoS[i],
						  SolnBlk.Grid.nfaceS(i, SolnBlk.JCl),
						  SolnBlk.Moving_wall_velocity,
						  FIXED_TEMPERATURE_WALL);
	SolnBlk.U[i][SolnBlk.JCl-2] = U(SolnBlk.W[i][SolnBlk.JCl-2]);
	break; 
      case BC_INFLOW_SUBSONIC :
	// all fixed except v.x (u) which is constant extrapolation
	SolnBlk.W[i][SolnBlk.JCl-1] = SolnBlk.WoS[i];
	SolnBlk.W[i][SolnBlk.JCl-1].v.y = SolnBlk.W[i][SolnBlk.JCl].v.y;
	SolnBlk.U[i][SolnBlk.JCl-1] = U(SolnBlk.W[i][SolnBlk.JCl-1]);
	
	SolnBlk.W[i][SolnBlk.JCl-2] = SolnBlk.WoS[i];
	SolnBlk.W[i][SolnBlk.JCl-2].v.y = SolnBlk.W[i][SolnBlk.JCl].v.y;
	SolnBlk.U[i][SolnBlk.JCl-2] = U(SolnBlk.W[i][SolnBlk.JCl-2]);
	break;
      case BC_OUTFLOW_SUBSONIC :
	// all constant extrapolation except pressure which is fixed.
	SolnBlk.W[i][SolnBlk.JCl-1] = SolnBlk.W[i][SolnBlk.JCl];
	SolnBlk.W[i][SolnBlk.JCl-1].p =  SolnBlk.WoS[i].p;
	SolnBlk.U[i][SolnBlk.JCl-1] = U(SolnBlk.W[i][SolnBlk.JCl-1]);
	  
	SolnBlk.W[i][SolnBlk.JCl-2] = SolnBlk.W[i][SolnBlk.JCl];
	SolnBlk.W[i][SolnBlk.JCl-2].p =  SolnBlk.WoS[i].p;
	SolnBlk.U[i][SolnBlk.JCl-2] = U(SolnBlk.W[i][SolnBlk.JCl-2]);
	break;
	
      case BC_FLAME_INFLOW :
	cout<<"\n BC issues ";
	break;
	
      case BC_FLAME_OUTFLOW :
	cout<<"\n BC issues ";
	break;
	
      default:
	SolnBlk.W[i][SolnBlk.JCl-1] = SolnBlk.W[i][SolnBlk.JCl];
	SolnBlk.U[i][SolnBlk.JCl-1] = SolnBlk.U[i][SolnBlk.JCl];
	SolnBlk.W[i][SolnBlk.JCl-2] = SolnBlk.W[i][SolnBlk.JCl];
	SolnBlk.U[i][SolnBlk.JCl-2] = SolnBlk.U[i][SolnBlk.JCl];
	break;
      } /* endswitch South */
      
      //NORTH        
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
      case BC_LINEAR_EXTRAPOLATION :
	Linear_Reconstruction_LeastSquares_2(SolnBlk, 
					     i, SolnBlk.JCu, 
					     LIMITER_BARTH_JESPERSEN);
	dX = SolnBlk.Grid.Cell[i][SolnBlk.JCu+1].Xc -
                 SolnBlk.Grid.Cell[i][SolnBlk.JCu].Xc;
	SolnBlk.W[i][SolnBlk.JCu+1] = SolnBlk.W[i][SolnBlk.JCu] + 
	  (SolnBlk.phi[i][SolnBlk.JCu]^SolnBlk.dWdx[i][SolnBlk.JCu])*dX.x +
	  (SolnBlk.phi[i][SolnBlk.JCu]^SolnBlk.dWdy[i][SolnBlk.JCu])*dX.y;
	SolnBlk.U[i][SolnBlk.JCu+1] = U(SolnBlk.W[i][SolnBlk.JCu+1]);
            dX = SolnBlk.Grid.Cell[i][SolnBlk.JCu+2].Xc -
	      SolnBlk.Grid.Cell[i][SolnBlk.JCu].Xc;
            SolnBlk.W[i][SolnBlk.JCu+2] = SolnBlk.W[i][SolnBlk.JCu] + 
	      (SolnBlk.phi[i][SolnBlk.JCu]^SolnBlk.dWdx[i][SolnBlk.JCu])*dX.x +
               (SolnBlk.phi[i][SolnBlk.JCu]^SolnBlk.dWdy[i][SolnBlk.JCu])*dX.y;
            SolnBlk.U[i][SolnBlk.JCu+2] = U(SolnBlk.W[i][SolnBlk.JCu+2]);
            break;
      case BC_REFLECTION :  
	//cout<<" "<<i;
	SolnBlk.W[i][SolnBlk.JCu+1] = Reflect(SolnBlk.W[i][SolnBlk.JCu],
					      SolnBlk.Grid.nfaceN(i, SolnBlk.JCu));
	SolnBlk.U[i][SolnBlk.JCu+1] = U(SolnBlk.W[i][SolnBlk.JCu+1]);
	SolnBlk.W[i][SolnBlk.JCu+2] = Reflect(SolnBlk.W[i][SolnBlk.JCu-1],
					      SolnBlk.Grid.nfaceN(i, SolnBlk.JCu));
	SolnBlk.U[i][SolnBlk.JCu+2] = U(SolnBlk.W[i][SolnBlk.JCu+2]);
	break;
      case BC_PERIODIC :
	SolnBlk.W[i][SolnBlk.JCu+1] = SolnBlk.W[i][SolnBlk.JCl+1];
	SolnBlk.U[i][SolnBlk.JCu+1] = SolnBlk.U[i][SolnBlk.JCl+1];
	SolnBlk.W[i][SolnBlk.JCu+2] = SolnBlk.W[i][SolnBlk.JCl+2];
	SolnBlk.U[i][SolnBlk.JCu+2] = SolnBlk.U[i][SolnBlk.JCl+2];
	break;
      case BC_CHARACTERISTIC :
	SolnBlk.W[i][SolnBlk.JCu+1] = 
	  BC_Characteristic_Pressure(SolnBlk.W[i][SolnBlk.JCu],
				     SolnBlk.WoN[i],
				     SolnBlk.Grid.nfaceN(i, SolnBlk.JCu));
	SolnBlk.U[i][SolnBlk.JCu+1] = U(SolnBlk.W[i][SolnBlk.JCu+1]);
	SolnBlk.W[i][SolnBlk.JCu+2] = SolnBlk.W[i][SolnBlk.JCu+1];
	SolnBlk.U[i][SolnBlk.JCu+2] = U(SolnBlk.W[i][SolnBlk.JCu+2]);
	break;
      case BC_FREE_SLIP :
	SolnBlk.W[i][SolnBlk.JCu+1] = Free_Slip(SolnBlk.W[i][SolnBlk.JCu],
						SolnBlk.WoN[i],
						SolnBlk.Grid.nfaceN(i, SolnBlk.JCu),
						FIXED_TEMPERATURE_WALL);
	SolnBlk.U[i][SolnBlk.JCu+1] = U(SolnBlk.W[i][SolnBlk.JCu+1]);
	SolnBlk.W[i][SolnBlk.JCu+2] = Free_Slip(SolnBlk.W[i][SolnBlk.JCu-1],
						SolnBlk.WoN[i],
						SolnBlk.Grid.nfaceN(i, SolnBlk.JCu),
						FIXED_TEMPERATURE_WALL);
	SolnBlk.U[i][SolnBlk.JCu+2] = U(SolnBlk.W[i][SolnBlk.JCu+2]);
      case BC_NO_SLIP :
	SolnBlk.W[i][SolnBlk.JCu+1] = No_Slip(SolnBlk.W[i][SolnBlk.JCu],
					      SolnBlk.WoN[i],
					      SolnBlk.Grid.nfaceN(i, SolnBlk.JCu),
					      FIXED_TEMPERATURE_WALL);
	SolnBlk.U[i][SolnBlk.JCu+1] = U(SolnBlk.W[i][SolnBlk.JCu+1]);
	SolnBlk.W[i][SolnBlk.JCu+2] = No_Slip(SolnBlk.W[i][SolnBlk.JCu-1],
					      SolnBlk.WoN[i],
					      SolnBlk.Grid.nfaceN(i, SolnBlk.JCu),
					      FIXED_TEMPERATURE_WALL);
	SolnBlk.U[i][SolnBlk.JCu+2] = U(SolnBlk.W[i][SolnBlk.JCu+2]);
	break;
      case BC_MOVING_WALL :
	//	if (i >= SolnBlk.ICl && j <= SolnBlk.ICu){
	SolnBlk.W[i][SolnBlk.JCu+1] = Moving_Wall(SolnBlk.W[i][SolnBlk.JCu],
						  SolnBlk.WoN[i],
						  SolnBlk.Grid.nfaceN(i, SolnBlk.JCu),
						  SolnBlk.Moving_wall_velocity,
						  FIXED_TEMPERATURE_WALL);
	SolnBlk.U[i][SolnBlk.JCu+1] = U(SolnBlk.W[i][SolnBlk.JCu+1]);
	SolnBlk.W[i][SolnBlk.JCu+2] = Moving_Wall(SolnBlk.W[i][SolnBlk.JCu-1],
						  SolnBlk.WoN[i],
						  SolnBlk.Grid.nfaceN(i, SolnBlk.JCu),
						  SolnBlk.Moving_wall_velocity,
						  FIXED_TEMPERATURE_WALL); 
	SolnBlk.U[i][SolnBlk.JCu+2] = U(SolnBlk.W[i][SolnBlk.JCu+2]);
	//	}
	break; 
      case BC_INFLOW_SUBSONIC :
	// all fixed except v.x (u) which is constant extrapolation
	SolnBlk.W[i][SolnBlk.JCu+1] = SolnBlk.WoN[i];
	SolnBlk.W[i][SolnBlk.JCu+1].v.y = SolnBlk.W[i][SolnBlk.JCu].v.y;
	SolnBlk.U[i][SolnBlk.JCu+1] = U(SolnBlk.W[i][SolnBlk.JCu+1]);
	
	SolnBlk.W[i][SolnBlk.JCu+2] = SolnBlk.WoN[i];
	SolnBlk.W[i][SolnBlk.JCu+2].v.y = SolnBlk.W[i][SolnBlk.JCu].v.y;
	SolnBlk.U[i][SolnBlk.JCu+2] = U(SolnBlk.W[i][SolnBlk.JCu+2]);
	break;
      case BC_OUTFLOW_SUBSONIC :
	// all constant extrapolation except pressure which is fixed.
	SolnBlk.W[i][SolnBlk.JCu+1] = SolnBlk.W[i][SolnBlk.JCu];
	SolnBlk.W[i][SolnBlk.JCu+1].p =  SolnBlk.WoN[i].p;
	SolnBlk.U[i][SolnBlk.JCu+1] = U(SolnBlk.W[i][SolnBlk.JCu+1]);
	
	SolnBlk.W[i][SolnBlk.JCu+2] = SolnBlk.W[i][SolnBlk.JCu];
	SolnBlk.W[i][SolnBlk.JCu+2].p =  SolnBlk.WoN[i].p;
	SolnBlk.U[i][SolnBlk.JCu+2] = U(SolnBlk.W[i][SolnBlk.JCu+2]);
	break;
	
      case BC_FLAME_INFLOW :
	cout<<"\n BC issues ";
	break;
	
      case BC_FLAME_OUTFLOW :
	SolnBlk.W[i][SolnBlk.JCu+1] = 
	  BC_Flame_Outflow(SolnBlk.W[i][SolnBlk.JCu],
			   SolnBlk.WoN[i],
			   SolnBlk.W[i][SolnBlk.JCl],
			   SolnBlk.Grid.nfaceN(i,SolnBlk.JCu));
	SolnBlk.U[i][SolnBlk.JCu+1] = U(SolnBlk.W[i][SolnBlk.JCu+1]);
	SolnBlk.W[i][SolnBlk.JCu+2] = SolnBlk.W[i][SolnBlk.JCu+1];
	SolnBlk.U[i][SolnBlk.JCu+2] = U(SolnBlk.W[i][SolnBlk.JCu+2]);
	break;
	
      default:
	SolnBlk.W[i][SolnBlk.JCu+1] = SolnBlk.W[i][SolnBlk.JCu];
	SolnBlk.U[i][SolnBlk.JCu+1] = SolnBlk.U[i][SolnBlk.JCu];
	SolnBlk.W[i][SolnBlk.JCu+2] = SolnBlk.W[i][SolnBlk.JCu];
	SolnBlk.U[i][SolnBlk.JCu+2] = SolnBlk.U[i][SolnBlk.JCu];
            break;
      } /* endswitch NORTH */ 
    } /* endfor SOUTH & NORTH */
      
}

/********************************************************
 * Routine: CFL                                         *
 *                                                      *
 * Determines the allowable global and local time steps *
 * (for explicit Chem time stepping scheme) for the     *
 * specified quadrilateral solution block according to  *
 * the Courant-Friedrichs-Lewy condition.               *
 *                                                      *
 ********************************************************/
double CFL(Chem2D_Quad_Block &SolnBlk,
           Chem2D_Input_Parameters &Input_Parameters) {

   
  double dtMin, dt_vis, dt_chem, d_i, d_j, v_i, v_j, a, rhomu, rhomut, delta_n;

    dtMin = MILLION;
    dt_vis = MILLION;
    dt_chem = MILLION;
  
    for ( int j  = SolnBlk.JCl-SolnBlk.Nghost ; j <= SolnBlk.JCu+SolnBlk.Nghost ; ++j ) {
       for ( int i = SolnBlk.ICl-SolnBlk.Nghost ; i <= SolnBlk.ICu+SolnBlk.Nghost ; ++i ) {
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

	     /******** Inviscid deltat calculation **********************/   
	     //no preconditioning

	     if(Input_Parameters.Preconditioning == 0){
      	    
	       a = SolnBlk.W[i][j].a();
	   
	       SolnBlk.dt[i][j] = min(d_i/(a+fabs(v_i)), d_j/(a+fabs(v_j)));
	   
	       //Low Mach Number Preconditioning
	     } else if(Input_Parameters.Preconditioning == 1) { 
	 
	       delta_n = min(fabs(d_i),fabs(d_j));
	       SolnBlk.dt[i][j] = min(d_i/SolnBlk.W[i][j].u_plus_aprecon(fabs(v_i),
									 delta_n),
	       			      d_j/SolnBlk.W[i][j].u_plus_aprecon(fabs(v_j),
									 delta_n));	    
	     }

	     /******** Viscous deltat calculation ************************/   
	     if (SolnBlk.Flow_Type != FLOWTYPE_INVISCID && 
                 Input_Parameters.Preconditioning == 0) {  
	       
	       rhomu = SolnBlk.W[i][j].mu()/SolnBlk.W[i][j].rho;
	       if(SolnBlk.Flow_Type == FLOWTYPE_TURBULENT_RANS_K_OMEGA ||
		  SolnBlk. Flow_Type == FLOWTYPE_TURBULENT_RANS_K_EPSILON) {
		 rhomut = SolnBlk.W[i][j].eddy_viscosity()/SolnBlk.W[i][j].rho; 
		 rhomu = max(rhomu, rhomut);
               }
	      
	       dt_vis = min((d_i*d_i)/(TWO*rhomu), (d_j*d_j)/(TWO*rhomu)); 
	       SolnBlk.dt[i][j]  = min(dt_vis, SolnBlk.dt[i][j]);
	     }	  
	     
	     /******** Chemical Source Term deltat calculation ************/   
	     if (SolnBlk.W[i][j].React.reactset_flag != NO_REACTIONS){
	       dt_chem = HALF/SolnBlk.W[i][j].dSwdU_max_diagonal(Input_Parameters.Preconditioning,
                                                                 delta_n);       
	       SolnBlk.dt[i][j] = min(dt_chem, SolnBlk.dt[i][j]);
	     }
	     
	     /************ Global Minimum ********************************/
	     dtMin = min(dtMin, SolnBlk.dt[i][j]);

	    //  cout<<" Viscous "<<dt_chem;  
// 	     cout<<" Chemistry "<<dt_chem;  
// 	     cout<<" Final " <<dtMin<<endl;

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
void Set_Global_TimeStep(Chem2D_Quad_Block &SolnBlk,
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
double L1_Norm_Residual(Chem2D_Quad_Block &SolnBlk) {

    int i, j;
    double l1_norm;

    l1_norm = ZERO;

    for ( j  = SolnBlk.JCl ; j <= SolnBlk.JCu ; ++j ) {
       for ( i = SolnBlk.ICl ; i <= SolnBlk.ICu ; ++i ) {
	 l1_norm += fabs(SolnBlk.dUdt[i][j][0][SolnBlk.residual_variable]);
	 //l1_norm += fabs(SolnBlk.dUdt[i][j][0].rho);
	 //l1_norm += abs(SolnBlk.dUdt[i][j][0].rhov); 
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
double L2_Norm_Residual(Chem2D_Quad_Block &SolnBlk) {

    int i, j;
    double l2_norm;

    l2_norm = ZERO;

    for ( j  = SolnBlk.JCl ; j <= SolnBlk.JCu ; ++j ) {
       for ( i = SolnBlk.ICl ; i <= SolnBlk.ICu ; ++i ) {
	 l2_norm += sqr(SolnBlk.dUdt[i][j][0][SolnBlk.residual_variable]);
	 //l2_norm += sqr(SolnBlk.dUdt[i][j][0].rho);
	 //l2_norm += sqr(SolnBlk.dUdt[i][j][0].rhov); 
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
double Max_Norm_Residual(Chem2D_Quad_Block &SolnBlk) {

    int i, j;
    double max_norm;

    max_norm = ZERO;

    for ( j  = SolnBlk.JCl ; j <= SolnBlk.JCu ; ++j ) {
       for ( i = SolnBlk.ICl ; i <= SolnBlk.ICu ; ++i ) {
	 max_norm = max(max_norm, fabs(SolnBlk.dUdt[i][j][0][SolnBlk.residual_variable]));
	 //max_norm = max(max_norm, fabs(SolnBlk.dUdt[i][j][0].rho));
	 //max_norm = max(max_norm, abs(SolnBlk.dUdt[i][j][0].rhov));
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
void Linear_Reconstruction_GreenGauss(Chem2D_Quad_Block &SolnBlk,
				      const int i, 
                                      const int j,
                                      const int Limiter) {

    int n, n2, n_pts, i_index[8], j_index[8];
    double u0Min, u0Max, uQuad[4], phi;
    double DxDx_ave, DxDy_ave, DyDy_ave;
    double l_north, l_south, l_east, l_west;
    Vector2D n_north, n_south, n_east, n_west, dX;
    Chem2D_pState W_nw, W_ne, W_sw, W_se, W_face, 
                   DU, DUDx_ave, DUDy_ave;
    
    int NUM_VAR_CHEM2D = SolnBlk.NumVar();
  
   /* Carry out the limited solution reconstruction in
       the specified cell of the computational mesh. */
    
    // Determine the number of neighbouring cells to
    // be used in the reconstruction procedure.  Away from
    // boundaries this 8 neighbours will be used.

//     if (i == SolnBlk.ICl-SolnBlk.Nghost || i == SolnBlk.ICu+SolnBlk.Nghost ||
//         j == SolnBlk.JCl-SolnBlk.Nghost || j == SolnBlk.JCu+SolnBlk.Nghost) {
//       n_pts = 0;
//     } else if ((i == SolnBlk.ICl-SolnBlk.Nghost+1) && 
//                (SolnBlk.Grid.BCtypeW[j] != BC_NONE)) {
//       if (j == SolnBlk.JCl-SolnBlk.Nghost+1 || j == SolnBlk.JCu+SolnBlk.Nghost-1) {
//          n_pts = 0;
//       } else if (SolnBlk.Grid.BCtypeW[j] == BC_PERIODIC ||
//                  SolnBlk.Grid.BCtypeW[j] == BC_CONSTANT_EXTRAPOLATION || 
//                  SolnBlk.Grid.BCtypeW[j] == BC_LINEAR_EXTRAPOLATION ||
//                  SolnBlk.Grid.BCtypeW[j] == BC_CHARACTERISTIC ||
// 		 SolnBlk.Grid.BCtypeW[j] == BC_FLAME_OUTFLOW) {
//          if (j == SolnBlk.JCl) {
//             n_pts = 5;
//             i_index[0] = i-1; j_index[0] = j  ;
//             i_index[1] = i+1; j_index[1] = j  ;
//             i_index[2] = i-1; j_index[2] = j+1;
//             i_index[3] = i  ; j_index[3] = j+1;
//             i_index[4] = i+1; j_index[4] = j+1;
//          } else if (j == SolnBlk.JCu) {
//             n_pts = 5;
//             i_index[0] = i-1; j_index[0] = j-1;
//             i_index[1] = i  ; j_index[1] = j-1;
//             i_index[2] = i+1; j_index[2] = j-1;
//             i_index[3] = i-1; j_index[3] = j  ;
//             i_index[4] = i+1; j_index[4] = j  ;
//          } else {
//             n_pts = 8;
//             i_index[0] = i-1; j_index[0] = j-1;
//             i_index[1] = i  ; j_index[1] = j-1;
//             i_index[2] = i+1; j_index[2] = j-1;
//             i_index[3] = i-1; j_index[3] = j  ;
//             i_index[4] = i+1; j_index[4] = j  ;
//             i_index[5] = i-1; j_index[5] = j+1;
//             i_index[6] = i  ; j_index[6] = j+1;
//             i_index[7] = i+1; j_index[7] = j+1;
//          } /* endif */
//       } else {
//          if (j == SolnBlk.JCl) {
//             n_pts = 3;
//             i_index[0] = i+1; j_index[0] = j  ;
//             i_index[1] = i  ; j_index[1] = j+1;
//             i_index[2] = i+1; j_index[2] = j+1;
//          } else if (j == SolnBlk.JCu) {
//             n_pts = 3;
//             i_index[0] = i  ; j_index[0] = j-1;
//             i_index[1] = i+1; j_index[1] = j-1;
//             i_index[2] = i+1; j_index[2] = j  ;
//          } else {
//             n_pts = 5;
//             i_index[0] = i  ; j_index[0] = j-1;
//             i_index[1] = i+1; j_index[1] = j-1;
//             i_index[2] = i+1; j_index[2] = j  ;
//             i_index[3] = i  ; j_index[3] = j+1;
//             i_index[4] = i+1; j_index[4] = j+1;
//          } /* endif */
//       } /* endif */           
//     } else if ((i == SolnBlk.ICu+SolnBlk.Nghost-1) && 
//                (SolnBlk.Grid.BCtypeE[j] != BC_NONE)) {
//       if (j == SolnBlk.JCl-SolnBlk.Nghost+1 || j == SolnBlk.JCu+SolnBlk.Nghost-1) {
//          n_pts = 0;
//       } else if (SolnBlk.Grid.BCtypeE[j] == BC_PERIODIC ||
//                  SolnBlk.Grid.BCtypeE[j] == BC_CONSTANT_EXTRAPOLATION ||
//                  SolnBlk.Grid.BCtypeE[j] == BC_LINEAR_EXTRAPOLATION ||
//                  SolnBlk.Grid.BCtypeE[j] == BC_CHARACTERISTIC ||
// 		 SolnBlk.Grid.BCtypeE[j] == BC_FLAME_OUTFLOW) {
//          if (j == SolnBlk.JCl) {
//             n_pts = 5;
//             i_index[0] = i-1; j_index[0] = j  ;
//             i_index[1] = i+1; j_index[1] = j  ;
//             i_index[2] = i-1; j_index[2] = j+1;
//             i_index[3] = i  ; j_index[3] = j+1;
//             i_index[4] = i+1; j_index[4] = j+1;
//          } else if (j == SolnBlk.JCu) {
//             n_pts = 5;
//             i_index[0] = i-1; j_index[0] = j-1;
//             i_index[1] = i  ; j_index[1] = j-1;
//             i_index[2] = i+1; j_index[2] = j-1;
//             i_index[3] = i-1; j_index[3] = j  ;
//             i_index[4] = i+1; j_index[4] = j  ;
//          } else {
//             n_pts = 8;
//             i_index[0] = i-1; j_index[0] = j-1;
//             i_index[1] = i  ; j_index[1] = j-1;
//             i_index[2] = i+1; j_index[2] = j-1;
//             i_index[3] = i-1; j_index[3] = j  ;
//             i_index[4] = i+1; j_index[4] = j  ;
//             i_index[5] = i-1; j_index[5] = j+1;
//             i_index[6] = i  ; j_index[6] = j+1;
//             i_index[7] = i+1; j_index[7] = j+1;
//          } /* endif */
//       } else {
//          if (j == SolnBlk.JCl) {
//             n_pts = 3;
//             i_index[0] = i-1; j_index[0] = j  ;
//             i_index[1] = i-1; j_index[1] = j+1;
//             i_index[2] = i  ; j_index[2] = j+1;
//          } else if (j == SolnBlk.JCu) {
//             n_pts = 3;
//             i_index[0] = i-1; j_index[0] = j-1;
//             i_index[1] = i  ; j_index[1] = j-1;
//             i_index[2] = i-1; j_index[2] = j  ;
//          } else {
//             n_pts = 5;
//             i_index[0] = i-1; j_index[0] = j-1;
//             i_index[1] = i  ; j_index[1] = j-1;
//             i_index[2] = i-1; j_index[2] = j  ;
//             i_index[3] = i-1; j_index[3] = j+1;
//             i_index[4] = i  ; j_index[4] = j+1;
//          } /* endif */
//       } /* endif */
//     } else if ((j == SolnBlk.JCl-SolnBlk.Nghost+1) && 
//                (SolnBlk.Grid.BCtypeS[i] != BC_NONE)) {
//       if (i == SolnBlk.ICl-SolnBlk.Nghost+1 || i == SolnBlk.ICu+SolnBlk.Nghost-1) {
//          n_pts = 0;
//       } else if (SolnBlk.Grid.BCtypeS[i] == BC_PERIODIC ||
//                  SolnBlk.Grid.BCtypeS[i] == BC_CONSTANT_EXTRAPOLATION ||
//                  SolnBlk.Grid.BCtypeS[i] == BC_LINEAR_EXTRAPOLATION ||
//                  SolnBlk.Grid.BCtypeS[i] == BC_CHARACTERISTIC ||
// 		 SolnBlk.Grid.BCtypeS[i] == BC_FLAME_OUTFLOW) {
//          if (i == SolnBlk.ICl) {
//             n_pts = 5;
//             i_index[0] = i  ; j_index[0] = j-1;
//             i_index[1] = i+1; j_index[1] = j-1;
//             i_index[2] = i+1; j_index[2] = j  ;
//             i_index[3] = i  ; j_index[3] = j+1;
//             i_index[4] = i+1; j_index[4] = j+1;
//          } else if (i == SolnBlk.ICu) {
//             n_pts = 5;
//             i_index[0] = i-1; j_index[0] = j-1;
//             i_index[1] = i  ; j_index[1] = j-1;
//             i_index[2] = i-1; j_index[2] = j  ;
//             i_index[3] = i-1; j_index[3] = j+1;
//             i_index[4] = i  ; j_index[4] = j+1;
//          } else {
//             n_pts = 8;
//             i_index[0] = i-1; j_index[0] = j-1;
//             i_index[1] = i  ; j_index[1] = j-1;
//             i_index[2] = i+1; j_index[2] = j-1;
//             i_index[3] = i-1; j_index[3] = j  ;
//             i_index[4] = i+1; j_index[4] = j  ;
//             i_index[5] = i-1; j_index[5] = j+1;
//             i_index[6] = i  ; j_index[6] = j+1;
//             i_index[7] = i+1; j_index[7] = j+1;
//          } /* endif */
//       } else {
//          if (i == SolnBlk.ICl) {
//             n_pts = 3;
//             i_index[0] = i+1; j_index[0] = j  ;
//             i_index[1] = i  ; j_index[1] = j+1;
//             i_index[2] = i+1; j_index[2] = j+1;
//          } else if (i == SolnBlk.ICu) {
//             n_pts = 3;
//             i_index[0] = i-1; j_index[0] = j  ;
//             i_index[1] = i-1; j_index[1] = j+1;
//             i_index[2] = i  ; j_index[2] = j+1;
//          } else {
//             n_pts = 5;
//             i_index[0] = i-1; j_index[0] = j  ;
//             i_index[1] = i+1; j_index[1] = j  ;
//             i_index[2] = i-1; j_index[2] = j+1;
//             i_index[3] = i  ; j_index[3] = j+1;
//             i_index[4] = i+1; j_index[4] = j+1;
//          } /* endif */
//       } /* endif */
//     } else if ((j == SolnBlk.JCu+SolnBlk.Nghost-1) && 
//                (SolnBlk.Grid.BCtypeN[i] != BC_NONE)) {
//       if (i == SolnBlk.ICl-SolnBlk.Nghost+1 || i == SolnBlk.ICu+SolnBlk.Nghost-1) {
//          n_pts = 0;
//       } else if (SolnBlk.Grid.BCtypeN[i] == BC_PERIODIC ||
//                  SolnBlk.Grid.BCtypeN[i] == BC_CONSTANT_EXTRAPOLATION ||
//                  SolnBlk.Grid.BCtypeN[i] == BC_LINEAR_EXTRAPOLATION ||
//                  SolnBlk.Grid.BCtypeN[i] == BC_CHARACTERISTIC ||
// 		 SolnBlk.Grid.BCtypeN[i] == BC_FLAME_OUTFLOW) {
//          if (i == SolnBlk.ICl) {
//             n_pts = 5;
//             i_index[0] = i  ; j_index[0] = j-1;
//             i_index[1] = i+1; j_index[1] = j-1;
//             i_index[2] = i+1; j_index[2] = j  ;
//             i_index[3] = i  ; j_index[3] = j+1;
//             i_index[4] = i+1; j_index[4] = j+1;
//          } else if (i == SolnBlk.ICu) {
//             n_pts = 5;
//             i_index[0] = i-1; j_index[0] = j-1;
//             i_index[1] = i  ; j_index[1] = j-1;
//             i_index[2] = i-1; j_index[2] = j  ;
//             i_index[3] = i-1; j_index[3] = j+1;
//             i_index[4] = i  ; j_index[4] = j+1;
//          } else {
//             n_pts = 8;
//             i_index[0] = i-1; j_index[0] = j-1;
//             i_index[1] = i  ; j_index[1] = j-1;
//             i_index[2] = i+1; j_index[2] = j-1;
//             i_index[3] = i-1; j_index[3] = j  ;
//             i_index[4] = i+1; j_index[4] = j  ;
//             i_index[5] = i-1; j_index[5] = j+1;
//             i_index[6] = i  ; j_index[6] = j+1;
//             i_index[7] = i+1; j_index[7] = j+1;
//          } /* endif */
//       } else {
//          if (i == SolnBlk.ICl) {
//             n_pts = 3;
//             i_index[0] = i  ; j_index[0] = j-1;
//             i_index[1] = i+1; j_index[1] = j-1;
//             i_index[2] = i+1; j_index[2] = j  ;
//          } else if (i == SolnBlk.ICu) {
//             n_pts = 3;
//             i_index[0] = i-1; j_index[0] = j-1;
//             i_index[1] = i  ; j_index[1] = j-1;
//             i_index[2] = i-1; j_index[2] = j  ;
//          } else {
//             n_pts = 5;
//             i_index[0] = i-1; j_index[0] = j-1;
//             i_index[1] = i  ; j_index[1] = j-1;
//             i_index[2] = i+1; j_index[2] = j-1;
//             i_index[3] = i-1; j_index[3] = j  ;
//             i_index[4] = i+1; j_index[4] = j  ;
//          } /* endif */
//       } /* endif */
// //      } else if (((i == SolnBlk.ICl && j == SolnBlk.JCl) && 
// //                  ((SolnBlk.Grid.BCtypeW[j] == BC_REFLECTION &&
// //                    SolnBlk.Grid.BCtypeS[i] == BC_REFLECTION) ||
// //                   (SolnBlk.Grid.BCtypeW[j] == BC_BURNINGSURFACE &&
// //                    SolnBlk.Grid.BCtypeS[i] == BC_BURNINGSURFACE))) || 
// //                 ((i == SolnBlk.ICl && j == SolnBlk.JCu) && 
// //                  ((SolnBlk.Grid.BCtypeW[j] == BC_REFLECTION &&
// //                    SolnBlk.Grid.BCtypeN[i] == BC_REFLECTION) ||
// //                   (SolnBlk.Grid.BCtypeW[j] == BC_BURNINGSURFACE &&
// //                    SolnBlk.Grid.BCtypeN[i] == BC_BURNINGSURFACE))) ||
// //                 ((i == SolnBlk.ICu && j == SolnBlk.JCl) && 
// //                  ((SolnBlk.Grid.BCtypeE[j] == BC_REFLECTION &&
// //                    SolnBlk.Grid.BCtypeS[i] == BC_REFLECTION) ||
// //                   (SolnBlk.Grid.BCtypeE[j] == BC_BURNINGSURFACE &&
// //                    SolnBlk.Grid.BCtypeS[i] == BC_BURNINGSURFACE))) ||
// //                 ((i == SolnBlk.ICu && j == SolnBlk.JCu) && 
// //                  ((SolnBlk.Grid.BCtypeE[j] == BC_REFLECTION &&
// //                    SolnBlk.Grid.BCtypeN[i] == BC_REFLECTION) ||
// //                   (SolnBlk.Grid.BCtypeE[j] == BC_BURNINGSURFACE &&
// //                    SolnBlk.Grid.BCtypeN[i] == BC_BURNINGSURFACE)))) {
// //        n_pts = 8;
// //        i_index[0] = i-1; j_index[0] = j-1;
// //        i_index[1] = i  ; j_index[1] = j-1;
// //        i_index[2] = i+1; j_index[2] = j-1;
// //        i_index[3] = i-1; j_index[3] = j  ;
// //        i_index[4] = i+1; j_index[4] = j  ;
// //        i_index[5] = i-1; j_index[5] = j+1;
// //        i_index[6] = i  ; j_index[6] = j+1;
// //        i_index[7] = i+1; j_index[7] = j+1;
// //      } else if ((i == SolnBlk.ICl && j == SolnBlk.JCl) && 
// //                 (SolnBlk.Grid.BCtypeW[j] != BC_NONE &&
// //                  SolnBlk.Grid.BCtypeS[i] != BC_NONE)) {
// //        n_pts = 7;
// //        i_index[0] = i  ; j_index[0] = j-1;
// //        i_index[1] = i+1; j_index[1] = j-1;
// //        i_index[2] = i-1; j_index[2] = j  ;
// //        i_index[3] = i+1; j_index[3] = j  ;
// //        i_index[4] = i-1; j_index[4] = j+1;
// //        i_index[5] = i  ; j_index[5] = j+1;
// //        i_index[6] = i+1; j_index[6] = j+1;
// //      } else if ((i == SolnBlk.ICl && j == SolnBlk.JCu) && 
// //                 (SolnBlk.Grid.BCtypeW[j] != BC_NONE &&
// //                  SolnBlk.Grid.BCtypeN[i] != BC_NONE)) {
// //        n_pts = 7;
// //        i_index[0] = i-1; j_index[0] = j-1;
// //        i_index[1] = i  ; j_index[1] = j-1;
// //        i_index[2] = i+1; j_index[2] = j-1;
// //        i_index[3] = i-1; j_index[3] = j  ;
// //        i_index[4] = i+1; j_index[4] = j  ;
// //        i_index[5] = i  ; j_index[5] = j+1;
// //        i_index[6] = i+1; j_index[6] = j+1;
// //      } else if ((i == SolnBlk.ICu && j == SolnBlk.JCl) && 
// //                 (SolnBlk.Grid.BCtypeE[j] != BC_NONE &&
// //                  SolnBlk.Grid.BCtypeS[i] != BC_NONE)) {
// //        n_pts = 7;
// //        i_index[0] = i-1; j_index[0] = j-1;
// //        i_index[1] = i  ; j_index[1] = j-1;
// //        i_index[2] = i-1; j_index[2] = j  ;
// //        i_index[3] = i+1; j_index[3] = j  ;
// //        i_index[4] = i-1; j_index[4] = j+1;
// //        i_index[5] = i  ; j_index[5] = j+1;
// //        i_index[6] = i+1; j_index[6] = j+1;
// //      } else if ((i == SolnBlk.ICu && j == SolnBlk.JCu) && 
// //                 (SolnBlk.Grid.BCtypeE[j] != BC_NONE &&
// //                  SolnBlk.Grid.BCtypeN[i] != BC_NONE)) {
// //        n_pts = 7;
// //        i_index[0] = i-1; j_index[0] = j-1;
// //        i_index[1] = i  ; j_index[1] = j-1;
// //        i_index[2] = i+1; j_index[2] = j-1;
// //        i_index[3] = i-1; j_index[3] = j  ;
// //        i_index[4] = i+1; j_index[4] = j  ;
// //        i_index[5] = i-1; j_index[5] = j+1;
// //        i_index[6] = i  ; j_index[6] = j+1;
//     } else {
//       n_pts = 8;
//       i_index[0] = i-1; j_index[0] = j-1;
//       i_index[1] = i  ; j_index[1] = j-1;
//       i_index[2] = i+1; j_index[2] = j-1;
//       i_index[3] = i-1; j_index[3] = j  ;
//       i_index[4] = i+1; j_index[4] = j  ;
//       i_index[5] = i-1; j_index[5] = j+1;
//       i_index[6] = i  ; j_index[6] = j+1;
//       i_index[7] = i+1; j_index[7] = j+1;
//     } /* endif */
    
/****************************************************************/
    //FOR VISCOUS -> CHANGED TO USE ALL 8
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
    
    /****************************************************************/
    // Perform reconstruction.    
    if (n_pts > 0) {
        // If 8 neighbours are used, apply Green-Gauss reconstruction
        //if (n_pts == 8) {
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
 //        } else {
// 	  DUDx_ave.Vacuum();
// 	  DUDy_ave.Vacuum();
//            DxDx_ave = ZERO;
//            DxDy_ave = ZERO;
//            DyDy_ave = ZERO;
    
//            for ( n2 = 0 ; n2 <= n_pts-1 ; ++n2 ) {
//                dX = SolnBlk.Grid.Cell[ i_index[n2] ][ j_index[n2] ].Xc - 
//                     SolnBlk.Grid.Cell[i][j].Xc;
//                DU = SolnBlk.W[ i_index[n2] ][ j_index[n2] ] - 
//                     SolnBlk.W[i][j];
//                DUDx_ave += DU*dX.x;
//                DUDy_ave += DU*dX.y;
//                DxDx_ave += dX.x*dX.x;
//                DxDy_ave += dX.x*dX.y;
//                DyDy_ave += dX.y*dX.y;
//            } /* endfor */
    					    
//            DUDx_ave = DUDx_ave/double(n_pts);
//            DUDy_ave = DUDy_ave/double(n_pts);
//            DxDx_ave = DxDx_ave/double(n_pts);
//            DxDy_ave = DxDy_ave/double(n_pts);
//            DyDy_ave = DyDy_ave/double(n_pts);
//            SolnBlk.dWdx[i][j] = (DUDx_ave*DyDy_ave-DUDy_ave*DxDy_ave)/
//                                 (DxDx_ave*DyDy_ave-DxDy_ave*DxDy_ave);
//            SolnBlk.dWdy[i][j] = (DUDy_ave*DxDx_ave-DUDx_ave*DxDy_ave)/
//                                 (DxDx_ave*DyDy_ave-DxDy_ave*DxDy_ave);
//         } /* endif */

	// Calculate slope limiters.    
       if (!SolnBlk.Freeze_Limiter) {
	  for ( n = 1 ; n <= NUM_VAR_CHEM2D ; ++n ) {
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
	} // end limiter if
	  
    } else {
      SolnBlk.dWdx[i][j].Vacuum();
      SolnBlk.dWdy[i][j].Vacuum();
      SolnBlk.phi[i][j].Vacuum();
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
void Linear_Reconstruction_GreenGauss(Chem2D_Quad_Block &SolnBlk,
				      const int Limiter) {

    /* Carry out the limited solution reconstruction in
       each cell of the computational mesh. */
    for (int j  = SolnBlk.JCl-SolnBlk.Nghost+1 ; j <= SolnBlk.JCu+SolnBlk.Nghost-1 ; ++j ) {
      for ( int i = SolnBlk.ICl-SolnBlk.Nghost+1 ; i <= SolnBlk.ICu+SolnBlk.Nghost-1 ; ++i ) {
	Linear_Reconstruction_GreenGauss(SolnBlk, i, j, Limiter);
      } /* endfor */
    } /* endfor */

}

// /********************************************************
//  * Routine: Linear_Reconstruction_LeastSquares          *
//  *                                                      *
//  * Performs the reconstruction of a limited piecewise   *
//  * linear solution state within a given cell (i,j) of   *
//  * the computational mesh for the specified             *
//  * quadrilateral solution block.  A least squares       *
//  * approach is used in the evaluation of the unlimited  *
//  * solution gradients.  Several slope limiters may be   *
//  * used.                                                *
//  *                                                      *
//  ********************************************************/
// void Linear_Reconstruction_LeastSquares(Chem2D_Quad_Block &SolnBlk,
// 				        const int i, 
//                                         const int j,
//                                         const int Limiter) {

//     int n, n2, n_pts, i_index[8], j_index[8];
//     double u0Min, u0Max, uQuad[4], phi;
//     double DxDx_ave, DxDy_ave, DyDy_ave;
//     Vector2D dX;
//     Chem2D_pState DU, DUDx_ave, DUDy_ave;
    
//     int NUM_VAR_CHEM2D = SolnBlk.NumVar();

//     /* Carry out the limited solution reconstruction in
//        each cell of the computational mesh. */
    
// //     if (i == SolnBlk.ICl-SolnBlk.Nghost || i == SolnBlk.ICu+SolnBlk.Nghost ||
// //         j == SolnBlk.JCl-SolnBlk.Nghost || j == SolnBlk.JCu+SolnBlk.Nghost) {
// //       n_pts = 0;
// //     } else if ((i == SolnBlk.ICl-SolnBlk.Nghost+1) && 
// //                (SolnBlk.Grid.BCtypeW[j] != BC_NONE)) {
// //       if (j == SolnBlk.JCl-SolnBlk.Nghost+1 || j == SolnBlk.JCu+SolnBlk.Nghost-1) {
// // 	n_pts = 0;
// //       } else if (SolnBlk.Grid.BCtypeW[j] == BC_PERIODIC ||
// //                  SolnBlk.Grid.BCtypeW[j] == BC_CONSTANT_EXTRAPOLATION ||
// //                  SolnBlk.Grid.BCtypeW[j] == BC_LINEAR_EXTRAPOLATION ||
// //                  SolnBlk.Grid.BCtypeW[j] == BC_CHARACTERISTIC ||
// // 		    SolnBlk.Grid.BCtypeW[j] == BC_FLAME_OUTFLOW) {
// //          if (j == SolnBlk.JCl) {
// //             n_pts = 5;
// //             i_index[0] = i-1; j_index[0] = j  ;
// //             i_index[1] = i+1; j_index[1] = j  ;
// //             i_index[2] = i-1; j_index[2] = j+1;
// //             i_index[3] = i  ; j_index[3] = j+1;
// //             i_index[4] = i+1; j_index[4] = j+1;
// //          } else if (j == SolnBlk.JCu) {
// //             n_pts = 5;
// //             i_index[0] = i-1; j_index[0] = j-1;
// //             i_index[1] = i  ; j_index[1] = j-1;
// //             i_index[2] = i+1; j_index[2] = j-1;
// //             i_index[3] = i-1; j_index[3] = j  ;
// //             i_index[4] = i+1; j_index[4] = j  ;
// //          } else {
// //             n_pts = 8;
// //             i_index[0] = i-1; j_index[0] = j-1;
// //             i_index[1] = i  ; j_index[1] = j-1;
// //             i_index[2] = i+1; j_index[2] = j-1;
// //             i_index[3] = i-1; j_index[3] = j  ;
// //             i_index[4] = i+1; j_index[4] = j  ;
// //             i_index[5] = i-1; j_index[5] = j+1;
// //             i_index[6] = i  ; j_index[6] = j+1;
// //             i_index[7] = i+1; j_index[7] = j+1;
// //          } /* endif */
// //       } else {
// //          if (j == SolnBlk.JCl) {
// //             n_pts = 3;
// //             i_index[0] = i+1; j_index[0] = j  ;
// //             i_index[1] = i  ; j_index[1] = j+1;
// //             i_index[2] = i+1; j_index[2] = j+1;
// //          } else if (j == SolnBlk.JCu) {
// //             n_pts = 3;
// //             i_index[0] = i  ; j_index[0] = j-1;
// //             i_index[1] = i+1; j_index[1] = j-1;
// //             i_index[2] = i+1; j_index[2] = j  ;
// //          } else {
// //             n_pts = 5;
// //             i_index[0] = i  ; j_index[0] = j-1;
// //             i_index[1] = i+1; j_index[1] = j-1;
// //             i_index[2] = i+1; j_index[2] = j  ;
// //             i_index[3] = i  ; j_index[3] = j+1;
// //             i_index[4] = i+1; j_index[4] = j+1;
// //          } /* endif */
// //       } /* endif */           
// //     } else if ((i == SolnBlk.ICu+SolnBlk.Nghost-1) && 
// //                (SolnBlk.Grid.BCtypeE[j] != BC_NONE)) {
// //       if (j == SolnBlk.JCl-SolnBlk.Nghost+1 || j == SolnBlk.JCu+SolnBlk.Nghost-1) {
// //          n_pts = 0;
// //       } else if (SolnBlk.Grid.BCtypeE[j] == BC_PERIODIC ||
// //                  SolnBlk.Grid.BCtypeE[j] == BC_CONSTANT_EXTRAPOLATION ||
// //                  SolnBlk.Grid.BCtypeE[j] == BC_LINEAR_EXTRAPOLATION ||
// //                  SolnBlk.Grid.BCtypeE[j] == BC_CHARACTERISTIC ||
// // 		 SolnBlk.Grid.BCtypeE[j] == BC_FLAME_OUTFLOW) {
// //          if (j == SolnBlk.JCl) {
// //             n_pts = 5;
// //             i_index[0] = i-1; j_index[0] = j  ;
// //             i_index[1] = i+1; j_index[1] = j  ;
// //             i_index[2] = i-1; j_index[2] = j+1;
// //             i_index[3] = i  ; j_index[3] = j+1;
// //             i_index[4] = i+1; j_index[4] = j+1;
// //          } else if (j == SolnBlk.JCu) {
// //             n_pts = 5;
// //             i_index[0] = i-1; j_index[0] = j-1;
// //             i_index[1] = i  ; j_index[1] = j-1;
// //             i_index[2] = i+1; j_index[2] = j-1;
// //             i_index[3] = i-1; j_index[3] = j  ;
// //             i_index[4] = i+1; j_index[4] = j  ;
// //          } else {
// //             n_pts = 8;
// //             i_index[0] = i-1; j_index[0] = j-1;
// //             i_index[1] = i  ; j_index[1] = j-1;
// //             i_index[2] = i+1; j_index[2] = j-1;
// //             i_index[3] = i-1; j_index[3] = j  ;
// //             i_index[4] = i+1; j_index[4] = j  ;
// //             i_index[5] = i-1; j_index[5] = j+1;
// //             i_index[6] = i  ; j_index[6] = j+1;
// //             i_index[7] = i+1; j_index[7] = j+1;
// //          } /* endif */
// //       } else {
// //          if (j == SolnBlk.JCl) {
// //             n_pts = 3;
// //             i_index[0] = i-1; j_index[0] = j  ;
// //             i_index[1] = i-1; j_index[1] = j+1;
// //             i_index[2] = i  ; j_index[2] = j+1;
// //          } else if (j == SolnBlk.JCu) {
// //             n_pts = 3;
// //             i_index[0] = i-1; j_index[0] = j-1;
// //             i_index[1] = i  ; j_index[1] = j-1;
// //             i_index[2] = i-1; j_index[2] = j  ;
// //          } else {
// //             n_pts = 5;
// //             i_index[0] = i-1; j_index[0] = j-1;
// //             i_index[1] = i  ; j_index[1] = j-1;
// //             i_index[2] = i-1; j_index[2] = j  ;
// //             i_index[3] = i-1; j_index[3] = j+1;
// //             i_index[4] = i  ; j_index[4] = j+1;
// //          } /* endif */
// //       } /* endif */
// //     } else if ((j == SolnBlk.JCl-SolnBlk.Nghost+1) && 
// //                (SolnBlk.Grid.BCtypeS[i] != BC_NONE)) {
// //       if (i == SolnBlk.ICl-SolnBlk.Nghost+1 || i == SolnBlk.ICu+SolnBlk.Nghost-1) {
// //          n_pts = 0;
// //       } else if (SolnBlk.Grid.BCtypeS[i] == BC_PERIODIC ||
// //                  SolnBlk.Grid.BCtypeS[i] == BC_CONSTANT_EXTRAPOLATION ||
// //                  SolnBlk.Grid.BCtypeS[i] == BC_LINEAR_EXTRAPOLATION ||
// //                  SolnBlk.Grid.BCtypeS[i] == BC_CHARACTERISTIC ||
// // 		 SolnBlk.Grid.BCtypeS[i] == BC_FLAME_OUTFLOW) {
// //          if (i == SolnBlk.ICl) {
// //             n_pts = 5;
// //             i_index[0] = i  ; j_index[0] = j-1;
// //             i_index[1] = i+1; j_index[1] = j-1;
// //             i_index[2] = i+1; j_index[2] = j  ;
// //             i_index[3] = i  ; j_index[3] = j+1;
// //             i_index[4] = i+1; j_index[4] = j+1;
// //          } else if (i == SolnBlk.ICu) {
// //             n_pts = 5;
// //             i_index[0] = i-1; j_index[0] = j-1;
// //             i_index[1] = i  ; j_index[1] = j-1;
// //             i_index[2] = i-1; j_index[2] = j  ;
// //             i_index[3] = i-1; j_index[3] = j+1;
// //             i_index[4] = i  ; j_index[4] = j+1;
// //          } else {
// //             n_pts = 8;
// //             i_index[0] = i-1; j_index[0] = j-1;
// //             i_index[1] = i  ; j_index[1] = j-1;
// //             i_index[2] = i+1; j_index[2] = j-1;
// //             i_index[3] = i-1; j_index[3] = j  ;
// //             i_index[4] = i+1; j_index[4] = j  ;
// //             i_index[5] = i-1; j_index[5] = j+1;
// //             i_index[6] = i  ; j_index[6] = j+1;
// //             i_index[7] = i+1; j_index[7] = j+1;
// //          } /* endif */
// //       } else {
// //          if (i == SolnBlk.ICl) {
// //             n_pts = 3;
// //             i_index[0] = i+1; j_index[0] = j  ;
// //             i_index[1] = i  ; j_index[1] = j+1;
// //             i_index[2] = i+1; j_index[2] = j+1;
// //          } else if (i == SolnBlk.ICu) {
// //             n_pts = 3;
// //             i_index[0] = i-1; j_index[0] = j  ;
// //             i_index[1] = i-1; j_index[1] = j+1;
// //             i_index[2] = i  ; j_index[2] = j+1;
// //          } else {
// //             n_pts = 5;
// //             i_index[0] = i-1; j_index[0] = j  ;
// //             i_index[1] = i+1; j_index[1] = j  ;
// //             i_index[2] = i-1; j_index[2] = j+1;
// //             i_index[3] = i  ; j_index[3] = j+1;
// //             i_index[4] = i+1; j_index[4] = j+1;
// //          } /* endif */
// //       } /* endif */
// //     } else if ((j == SolnBlk.JCu+SolnBlk.Nghost-1) && 
// //                (SolnBlk.Grid.BCtypeN[i] != BC_NONE)) {
// //       if (i == SolnBlk.ICl-SolnBlk.Nghost+1 || i == SolnBlk.ICu+SolnBlk.Nghost-1) {
// //          n_pts = 0;
// //       } else if (SolnBlk.Grid.BCtypeN[i] == BC_PERIODIC ||
// //                  SolnBlk.Grid.BCtypeN[i] == BC_CONSTANT_EXTRAPOLATION ||
// //                  SolnBlk.Grid.BCtypeN[i] == BC_LINEAR_EXTRAPOLATION ||
// //                  SolnBlk.Grid.BCtypeN[i] == BC_CHARACTERISTIC ||
// // 		 SolnBlk.Grid.BCtypeN[i] == BC_FLAME_OUTFLOW) {
// //          if (i == SolnBlk.ICl) {
// //             n_pts = 5;
// //             i_index[0] = i  ; j_index[0] = j-1;
// //             i_index[1] = i+1; j_index[1] = j-1;
// //             i_index[2] = i+1; j_index[2] = j  ;
// //             i_index[3] = i  ; j_index[3] = j+1;
// //             i_index[4] = i+1; j_index[4] = j+1;
// //          } else if (i == SolnBlk.ICu) {
// //             n_pts = 5;
// //             i_index[0] = i-1; j_index[0] = j-1;
// //             i_index[1] = i  ; j_index[1] = j-1;
// //             i_index[2] = i-1; j_index[2] = j  ;
// //             i_index[3] = i-1; j_index[3] = j+1;
// //             i_index[4] = i  ; j_index[4] = j+1;
// //          } else {
// //             n_pts = 8;
// //             i_index[0] = i-1; j_index[0] = j-1;
// //             i_index[1] = i  ; j_index[1] = j-1;
// //             i_index[2] = i+1; j_index[2] = j-1;
// //             i_index[3] = i-1; j_index[3] = j  ;
// //             i_index[4] = i+1; j_index[4] = j  ;
// //             i_index[5] = i-1; j_index[5] = j+1;
// //             i_index[6] = i  ; j_index[6] = j+1;
// //             i_index[7] = i+1; j_index[7] = j+1;
// //          } /* endif */
// //       } else {
// //          if (i == SolnBlk.ICl) {
// //             n_pts = 3;
// //             i_index[0] = i  ; j_index[0] = j-1;
// //             i_index[1] = i+1; j_index[1] = j-1;
// //             i_index[2] = i+1; j_index[2] = j  ;
// //          } else if (i == SolnBlk.ICu) {
// //             n_pts = 3;
// //             i_index[0] = i-1; j_index[0] = j-1;
// //             i_index[1] = i  ; j_index[1] = j-1;
// //             i_index[2] = i-1; j_index[2] = j  ;
// //          } else {
// //             n_pts = 5;
// //             i_index[0] = i-1; j_index[0] = j-1;
// //             i_index[1] = i  ; j_index[1] = j-1;
// //             i_index[2] = i+1; j_index[2] = j-1;
// //             i_index[3] = i-1; j_index[3] = j  ;
// //             i_index[4] = i+1; j_index[4] = j  ;
// //          } /* endif */
// //       } /* endif */
// // //      } else if (((i == SolnBlk.ICl && j == SolnBlk.JCl) && 
// // //                  ((SolnBlk.Grid.BCtypeW[j] == BC_REFLECTION &&
// // //                    SolnBlk.Grid.BCtypeS[i] == BC_REFLECTION) ||
// // //                   (SolnBlk.Grid.BCtypeW[j] == BC_BURNINGSURFACE &&
// // //                    SolnBlk.Grid.BCtypeS[i] == BC_BURNINGSURFACE))) || 
// // //                 ((i == SolnBlk.ICl && j == SolnBlk.JCu) && 
// // //                  ((SolnBlk.Grid.BCtypeW[j] == BC_REFLECTION &&
// // //                    SolnBlk.Grid.BCtypeN[i] == BC_REFLECTION) ||
// // //                   (SolnBlk.Grid.BCtypeW[j] == BC_BURNINGSURFACE &&
// // //                    SolnBlk.Grid.BCtypeN[i] == BC_BURNINGSURFACE))) ||
// // //                 ((i == SolnBlk.ICu && j == SolnBlk.JCl) && 
// // //                  ((SolnBlk.Grid.BCtypeE[j] == BC_REFLECTION &&
// // //                    SolnBlk.Grid.BCtypeS[i] == BC_REFLECTION) ||
// // //                   (SolnBlk.Grid.BCtypeE[j] == BC_BURNINGSURFACE &&
// // //                    SolnBlk.Grid.BCtypeS[i] == BC_BURNINGSURFACE))) ||
// // //                 ((i == SolnBlk.ICu && j == SolnBlk.JCu) && 
// // //                  ((SolnBlk.Grid.BCtypeE[j] == BC_REFLECTION &&
// // //                    SolnBlk.Grid.BCtypeN[i] == BC_REFLECTION) ||
// // //                   (SolnBlk.Grid.BCtypeE[j] == BC_BURNINGSURFACE &&
// // //                    SolnBlk.Grid.BCtypeN[i] == BC_BURNINGSURFACE)))) {
// // //        n_pts = 8;
// // //        i_index[0] = i-1; j_index[0] = j-1;
// // //        i_index[1] = i  ; j_index[1] = j-1;
// // //        i_index[2] = i+1; j_index[2] = j-1;
// // //        i_index[3] = i-1; j_index[3] = j  ;
// // //        i_index[4] = i+1; j_index[4] = j  ;
// // //        i_index[5] = i-1; j_index[5] = j+1;
// // //        i_index[6] = i  ; j_index[6] = j+1;
// // //        i_index[7] = i+1; j_index[7] = j+1;
// // //      } else if ((i == SolnBlk.ICl && j == SolnBlk.JCl) && 
// // //                 (SolnBlk.Grid.BCtypeW[j] != BC_NONE &&
// // //                  SolnBlk.Grid.BCtypeS[i] != BC_NONE)) {
// // //        n_pts = 7;
// // //        i_index[0] = i  ; j_index[0] = j-1;
// // //        i_index[1] = i+1; j_index[1] = j-1;
// // //        i_index[2] = i-1; j_index[2] = j  ;
// // //        i_index[3] = i+1; j_index[3] = j  ;
// // //        i_index[4] = i-1; j_index[4] = j+1;
// // //        i_index[5] = i  ; j_index[5] = j+1;
// // //        i_index[6] = i+1; j_index[6] = j+1;
// // //      } else if ((i == SolnBlk.ICl && j == SolnBlk.JCu) && 
// // //                 (SolnBlk.Grid.BCtypeW[j] != BC_NONE &&
// // //                  SolnBlk.Grid.BCtypeN[i] != BC_NONE)) {
// // //        n_pts = 7;
// // //        i_index[0] = i-1; j_index[0] = j-1;
// // //        i_index[1] = i  ; j_index[1] = j-1;
// // //        i_index[2] = i+1; j_index[2] = j-1;
// // //        i_index[3] = i-1; j_index[3] = j  ;
// // //        i_index[4] = i+1; j_index[4] = j  ;
// // //        i_index[5] = i  ; j_index[5] = j+1;
// // //        i_index[6] = i+1; j_index[6] = j+1;
// // //      } else if ((i == SolnBlk.ICu && j == SolnBlk.JCl) && 
// // //                 (SolnBlk.Grid.BCtypeE[j] != BC_NONE &&
// // //                  SolnBlk.Grid.BCtypeS[i] != BC_NONE)) {
// // //        n_pts = 7;
// // //        i_index[0] = i-1; j_index[0] = j-1;
// // //        i_index[1] = i  ; j_index[1] = j-1;
// // //        i_index[2] = i-1; j_index[2] = j  ;
// // //        i_index[3] = i+1; j_index[3] = j  ;
// // //        i_index[4] = i-1; j_index[4] = j+1;
// // //        i_index[5] = i  ; j_index[5] = j+1;
// // //        i_index[6] = i+1; j_index[6] = j+1;
// // //      } else if ((i == SolnBlk.ICu && j == SolnBlk.JCu) && 
// // //                 (SolnBlk.Grid.BCtypeE[j] != BC_NONE &&
// // //                  SolnBlk.Grid.BCtypeN[i] != BC_NONE)) {
// // //        n_pts = 7;
// // //        i_index[0] = i-1; j_index[0] = j-1;
// // //        i_index[1] = i  ; j_index[1] = j-1;
// // //        i_index[2] = i+1; j_index[2] = j-1;
// // //        i_index[3] = i-1; j_index[3] = j  ;
// // //        i_index[4] = i+1; j_index[4] = j  ;
// // //        i_index[5] = i-1; j_index[5] = j+1;
// // //        i_index[6] = i  ; j_index[6] = j+1;
// //     } else {
// //       n_pts = 8;
// //       i_index[0] = i-1; j_index[0] = j-1;
// //       i_index[1] = i  ; j_index[1] = j-1;
// //       i_index[2] = i+1; j_index[2] = j-1;
// //       i_index[3] = i-1; j_index[3] = j  ;
// //       i_index[4] = i+1; j_index[4] = j  ;
// //       i_index[5] = i-1; j_index[5] = j+1;
// //       i_index[6] = i  ; j_index[6] = j+1;
// //       i_index[7] = i+1; j_index[7] = j+1;
// //     } /* endif */

//     /****************************************************************/
//     //FOR VISCOUS -> CHANGED TO USE ALL 8
//     if (i == SolnBlk.ICl-SolnBlk.Nghost || i == SolnBlk.ICu+SolnBlk.Nghost ||
// 	j == SolnBlk.JCl-SolnBlk.Nghost || j == SolnBlk.JCu+SolnBlk.Nghost) {
//       n_pts = 0;
//     } else {
//       n_pts = 8;
//       i_index[0] = i-1; j_index[0] = j-1;
//       i_index[1] = i  ; j_index[1] = j-1;
//       i_index[2] = i+1; j_index[2] = j-1;
//       i_index[3] = i-1; j_index[3] = j  ;
//       i_index[4] = i+1; j_index[4] = j  ;
//       i_index[5] = i-1; j_index[5] = j+1;
//       i_index[6] = i  ; j_index[6] = j+1;
//       i_index[7] = i+1; j_index[7] = j+1;
//     }    
//     /****************************************************************/
    
//     if (n_pts > 0) {
//         DUDx_ave.Vacuum();
//         DUDy_ave.Vacuum();
//         DxDx_ave = ZERO;
//         DxDy_ave = ZERO;
//         DyDy_ave = ZERO;
    
//         for ( n2 = 0 ; n2 < n_pts ; n2++ ) {
//             dX = SolnBlk.Grid.Cell[ i_index[n2] ][ j_index[n2] ].Xc - 
//                  SolnBlk.Grid.Cell[i][j].Xc;
//             DU = SolnBlk.W[ i_index[n2] ][ j_index[n2] ] - 
//                  SolnBlk.W[i][j];
//             DUDx_ave += DU*dX.x;
//             DUDy_ave += DU*dX.y;
//             DxDx_ave += dX.x*dX.x;
//             DxDy_ave += dX.x*dX.y;
//             DyDy_ave += dX.y*dX.y;

//         } /* endfor */
    					    
//         DUDx_ave = DUDx_ave/double(n_pts);
//         DUDy_ave = DUDy_ave/double(n_pts);
//         DxDx_ave = DxDx_ave/double(n_pts);
//         DxDy_ave = DxDy_ave/double(n_pts);
//         DyDy_ave = DyDy_ave/double(n_pts);
//         SolnBlk.dWdx[i][j] = (DUDx_ave*DyDy_ave-DUDy_ave*DxDy_ave)/
//                              (DxDx_ave*DyDy_ave-DxDy_ave*DxDy_ave);
//         SolnBlk.dWdy[i][j] = (DUDy_ave*DxDx_ave-DUDx_ave*DxDy_ave)/
//                              (DxDx_ave*DyDy_ave-DxDy_ave*DxDy_ave);

//     	if (!SolnBlk.Freeze_Limiter) {
// 	  for ( n = 1 ; n <= NUM_VAR_CHEM2D ; ++n ) {
// 	    u0Min = SolnBlk.W[i][j][n];
// 	    u0Max = u0Min;
// 	    for ( n2 = 0 ; n2 <= n_pts-1 ; ++n2 ) {
//               u0Min = min(u0Min, SolnBlk.W[ i_index[n2] ][ j_index[n2] ][n]);
//               u0Max = max(u0Max, SolnBlk.W[ i_index[n2] ][ j_index[n2] ][n]);
// 	    } /* endfor */
	    
// 	    dX = SolnBlk.Grid.xfaceE(i, j)-SolnBlk.Grid.Cell[i][j].Xc;
// 	    uQuad[0] = SolnBlk.W[i][j][n] + 
// 	      SolnBlk.dWdx[i][j][n]*dX.x +
// 	      SolnBlk.dWdy[i][j][n]*dX.y ;
// 	    dX = SolnBlk.Grid.xfaceW(i, j)-SolnBlk.Grid.Cell[i][j].Xc;
// 	    uQuad[1] = SolnBlk.W[i][j][n] + 
// 	      SolnBlk.dWdx[i][j][n]*dX.x +
// 	      SolnBlk.dWdy[i][j][n]*dX.y ;
// 	    dX = SolnBlk.Grid.xfaceN(i, j)-SolnBlk.Grid.Cell[i][j].Xc;
// 	    uQuad[2] = SolnBlk.W[i][j][n] + 
// 	      SolnBlk.dWdx[i][j][n]*dX.x +
// 	      SolnBlk.dWdy[i][j][n]*dX.y ;
// 	    dX = SolnBlk.Grid.xfaceS(i, j)-SolnBlk.Grid.Cell[i][j].Xc;
// 	    uQuad[3] = SolnBlk.W[i][j][n] + 
// 	      SolnBlk.dWdx[i][j][n]*dX.x +
// 	      SolnBlk.dWdy[i][j][n]*dX.y ;
	    
// 	    switch(Limiter) {
// 	    case LIMITER_ONE :
//                phi = ONE;
//                break;
// 	    case LIMITER_ZERO :
// 	      phi = ZERO;
// 	      break;
// 	    case LIMITER_BARTH_JESPERSEN :
//                phi = Limiter_BarthJespersen(uQuad, SolnBlk.W[i][j][n], 
//                                             u0Min, u0Max, 4);
//                break;
// 	    case LIMITER_VENKATAKRISHNAN :
// 	      phi = Limiter_Venkatakrishnan(uQuad, SolnBlk.W[i][j][n], 
//                                              u0Min, u0Max, 4);
// 	      break;
// 	    case LIMITER_VANLEER :
//                phi = Limiter_VanLeer(uQuad, SolnBlk.W[i][j][n], 
//                                      u0Min, u0Max, 4);
//                break;
// 	    case LIMITER_VANALBADA :
// 	      phi = Limiter_VanAlbada(uQuad, SolnBlk.W[i][j][n], 
// 				      u0Min, u0Max, 4);
// 	      break;
// 	    default:
// 	      phi = Limiter_BarthJespersen(uQuad, SolnBlk.W[i][j][n], 
// 					   u0Min, u0Max, 4);
// 	      break;
// 	    } /* endswitch */
	    
// 	    SolnBlk.phi[i][j][n] = phi;
// 	  } /* endfor */
// 	}//end limiter if
//     } else {
//         SolnBlk.dWdx[i][j].Vacuum();
//         SolnBlk.dWdy[i][j].Vacuum();
//         SolnBlk.phi[i][j].Vacuum(); 
//     } /* endif */
    
// }

void Linear_Reconstruction_LeastSquares(Chem2D_Quad_Block &SolnBlk,
				        const int i, 
                                        const int j,
                                        const int Limiter) {
  //COME

    int n, n2, n_pts, i_index[8], j_index[8];
    double u0Min, u0Max, uQuad[4], phi;
    double DxDx_ave, DxDy_ave, DyDy_ave;
    Vector2D dX;
    Chem2D_pState DU, DUDx_ave, DUDy_ave;
    
    //For the derivatives of gradients purpose
    // i_index_neigbor [8][1]  1 -- North , 8 -- its 8 sorrounding cells
    // i_index_neigbor [8][2]  2 -- South , 8 -- its 8 sorrounding cells
    // i_index_neigbor [8][3]  3 -- South , 8 -- its 8 sorrounding cells
    // i_index_neigbor [8][4]  4 -- South , 8 -- its 8 sorrounding cells

    int i_index_neigbor [8][5], j_index_neigbor[8][5];
    double dxdx_neigbor, dxdy_neigbor, dydy_neigbor;
    Vector2D dX_neigbor;
    Chem2D_pState DW, DWDx_ave, DWDy_ave;

    int num=0;

    int NUM_VAR_CHEM2D = SolnBlk.NumVar();

    /* Carry out the limited solution reconstruction in
       each cell of the computational mesh. */
    
//     if (i == SolnBlk.ICl-SolnBlk.Nghost || i == SolnBlk.ICu+SolnBlk.Nghost ||
//         j == SolnBlk.JCl-SolnBlk.Nghost || j == SolnBlk.JCu+SolnBlk.Nghost) {
//       n_pts = 0;
//     } else if ((i == SolnBlk.ICl-SolnBlk.Nghost+1) && 
//                (SolnBlk.Grid.BCtypeW[j] != BC_NONE)) {
//       if (j == SolnBlk.JCl-SolnBlk.Nghost+1 || j == SolnBlk.JCu+SolnBlk.Nghost-1) {
// 	n_pts = 0;
//       } else if (SolnBlk.Grid.BCtypeW[j] == BC_PERIODIC ||
//                  SolnBlk.Grid.BCtypeW[j] == BC_CONSTANT_EXTRAPOLATION ||
//                  SolnBlk.Grid.BCtypeW[j] == BC_LINEAR_EXTRAPOLATION ||
//                  SolnBlk.Grid.BCtypeW[j] == BC_CHARACTERISTIC ||
// 		    SolnBlk.Grid.BCtypeW[j] == BC_FLAME_OUTFLOW) {
//          if (j == SolnBlk.JCl) {
//             n_pts = 5;
//             i_index[0] = i-1; j_index[0] = j  ;
//             i_index[1] = i+1; j_index[1] = j  ;
//             i_index[2] = i-1; j_index[2] = j+1;
//             i_index[3] = i  ; j_index[3] = j+1;
//             i_index[4] = i+1; j_index[4] = j+1;
//          } else if (j == SolnBlk.JCu) {
//             n_pts = 5;
//             i_index[0] = i-1; j_index[0] = j-1;
//             i_index[1] = i  ; j_index[1] = j-1;
//             i_index[2] = i+1; j_index[2] = j-1;
//             i_index[3] = i-1; j_index[3] = j  ;
//             i_index[4] = i+1; j_index[4] = j  ;
//          } else {
//             n_pts = 8;
//             i_index[0] = i-1; j_index[0] = j-1;
//             i_index[1] = i  ; j_index[1] = j-1;
//             i_index[2] = i+1; j_index[2] = j-1;
//             i_index[3] = i-1; j_index[3] = j  ;
//             i_index[4] = i+1; j_index[4] = j  ;
//             i_index[5] = i-1; j_index[5] = j+1;
//             i_index[6] = i  ; j_index[6] = j+1;
//             i_index[7] = i+1; j_index[7] = j+1;
//          } /* endif */
//       } else {
//          if (j == SolnBlk.JCl) {
//             n_pts = 3;
//             i_index[0] = i+1; j_index[0] = j  ;
//             i_index[1] = i  ; j_index[1] = j+1;
//             i_index[2] = i+1; j_index[2] = j+1;
//          } else if (j == SolnBlk.JCu) {
//             n_pts = 3;
//             i_index[0] = i  ; j_index[0] = j-1;
//             i_index[1] = i+1; j_index[1] = j-1;
//             i_index[2] = i+1; j_index[2] = j  ;
//          } else {
//             n_pts = 5;
//             i_index[0] = i  ; j_index[0] = j-1;
//             i_index[1] = i+1; j_index[1] = j-1;
//             i_index[2] = i+1; j_index[2] = j  ;
//             i_index[3] = i  ; j_index[3] = j+1;
//             i_index[4] = i+1; j_index[4] = j+1;
//          } /* endif */
//       } /* endif */           
//     } else if ((i == SolnBlk.ICu+SolnBlk.Nghost-1) && 
//                (SolnBlk.Grid.BCtypeE[j] != BC_NONE)) {
//       if (j == SolnBlk.JCl-SolnBlk.Nghost+1 || j == SolnBlk.JCu+SolnBlk.Nghost-1) {
//          n_pts = 0;
//       } else if (SolnBlk.Grid.BCtypeE[j] == BC_PERIODIC ||
//                  SolnBlk.Grid.BCtypeE[j] == BC_CONSTANT_EXTRAPOLATION ||
//                  SolnBlk.Grid.BCtypeE[j] == BC_LINEAR_EXTRAPOLATION ||
//                  SolnBlk.Grid.BCtypeE[j] == BC_CHARACTERISTIC ||
// 		 SolnBlk.Grid.BCtypeE[j] == BC_FLAME_OUTFLOW) {
//          if (j == SolnBlk.JCl) {
//             n_pts = 5;
//             i_index[0] = i-1; j_index[0] = j  ;
//             i_index[1] = i+1; j_index[1] = j  ;
//             i_index[2] = i-1; j_index[2] = j+1;
//             i_index[3] = i  ; j_index[3] = j+1;
//             i_index[4] = i+1; j_index[4] = j+1;
//          } else if (j == SolnBlk.JCu) {
//             n_pts = 5;
//             i_index[0] = i-1; j_index[0] = j-1;
//             i_index[1] = i  ; j_index[1] = j-1;
//             i_index[2] = i+1; j_index[2] = j-1;
//             i_index[3] = i-1; j_index[3] = j  ;
//             i_index[4] = i+1; j_index[4] = j  ;
//          } else {
//             n_pts = 8;
//             i_index[0] = i-1; j_index[0] = j-1;
//             i_index[1] = i  ; j_index[1] = j-1;
//             i_index[2] = i+1; j_index[2] = j-1;
//             i_index[3] = i-1; j_index[3] = j  ;
//             i_index[4] = i+1; j_index[4] = j  ;
//             i_index[5] = i-1; j_index[5] = j+1;
//             i_index[6] = i  ; j_index[6] = j+1;
//             i_index[7] = i+1; j_index[7] = j+1;
//          } /* endif */
//       } else {
//          if (j == SolnBlk.JCl) {
//             n_pts = 3;
//             i_index[0] = i-1; j_index[0] = j  ;
//             i_index[1] = i-1; j_index[1] = j+1;
//             i_index[2] = i  ; j_index[2] = j+1;
//          } else if (j == SolnBlk.JCu) {
//             n_pts = 3;
//             i_index[0] = i-1; j_index[0] = j-1;
//             i_index[1] = i  ; j_index[1] = j-1;
//             i_index[2] = i-1; j_index[2] = j  ;
//          } else {
//             n_pts = 5;
//             i_index[0] = i-1; j_index[0] = j-1;
//             i_index[1] = i  ; j_index[1] = j-1;
//             i_index[2] = i-1; j_index[2] = j  ;
//             i_index[3] = i-1; j_index[3] = j+1;
//             i_index[4] = i  ; j_index[4] = j+1;
//          } /* endif */
//       } /* endif */
//     } else if ((j == SolnBlk.JCl-SolnBlk.Nghost+1) && 
//                (SolnBlk.Grid.BCtypeS[i] != BC_NONE)) {
//       if (i == SolnBlk.ICl-SolnBlk.Nghost+1 || i == SolnBlk.ICu+SolnBlk.Nghost-1) {
//          n_pts = 0;
//       } else if (SolnBlk.Grid.BCtypeS[i] == BC_PERIODIC ||
//                  SolnBlk.Grid.BCtypeS[i] == BC_CONSTANT_EXTRAPOLATION ||
//                  SolnBlk.Grid.BCtypeS[i] == BC_LINEAR_EXTRAPOLATION ||
//                  SolnBlk.Grid.BCtypeS[i] == BC_CHARACTERISTIC ||
// 		 SolnBlk.Grid.BCtypeS[i] == BC_FLAME_OUTFLOW) {
//          if (i == SolnBlk.ICl) {
//             n_pts = 5;
//             i_index[0] = i  ; j_index[0] = j-1;
//             i_index[1] = i+1; j_index[1] = j-1;
//             i_index[2] = i+1; j_index[2] = j  ;
//             i_index[3] = i  ; j_index[3] = j+1;
//             i_index[4] = i+1; j_index[4] = j+1;
//          } else if (i == SolnBlk.ICu) {
//             n_pts = 5;
//             i_index[0] = i-1; j_index[0] = j-1;
//             i_index[1] = i  ; j_index[1] = j-1;
//             i_index[2] = i-1; j_index[2] = j  ;
//             i_index[3] = i-1; j_index[3] = j+1;
//             i_index[4] = i  ; j_index[4] = j+1;
//          } else {
//             n_pts = 8;
//             i_index[0] = i-1; j_index[0] = j-1;
//             i_index[1] = i  ; j_index[1] = j-1;
//             i_index[2] = i+1; j_index[2] = j-1;
//             i_index[3] = i-1; j_index[3] = j  ;
//             i_index[4] = i+1; j_index[4] = j  ;
//             i_index[5] = i-1; j_index[5] = j+1;
//             i_index[6] = i  ; j_index[6] = j+1;
//             i_index[7] = i+1; j_index[7] = j+1;
//          } /* endif */
//       } else {
//          if (i == SolnBlk.ICl) {
//             n_pts = 3;
//             i_index[0] = i+1; j_index[0] = j  ;
//             i_index[1] = i  ; j_index[1] = j+1;
//             i_index[2] = i+1; j_index[2] = j+1;
//          } else if (i == SolnBlk.ICu) {
//             n_pts = 3;
//             i_index[0] = i-1; j_index[0] = j  ;
//             i_index[1] = i-1; j_index[1] = j+1;
//             i_index[2] = i  ; j_index[2] = j+1;
//          } else {
//             n_pts = 5;
//             i_index[0] = i-1; j_index[0] = j  ;
//             i_index[1] = i+1; j_index[1] = j  ;
//             i_index[2] = i-1; j_index[2] = j+1;
//             i_index[3] = i  ; j_index[3] = j+1;
//             i_index[4] = i+1; j_index[4] = j+1;
//          } /* endif */
//       } /* endif */
//     } else if ((j == SolnBlk.JCu+SolnBlk.Nghost-1) && 
//                (SolnBlk.Grid.BCtypeN[i] != BC_NONE)) {
//       if (i == SolnBlk.ICl-SolnBlk.Nghost+1 || i == SolnBlk.ICu+SolnBlk.Nghost-1) {
//          n_pts = 0;
//       } else if (SolnBlk.Grid.BCtypeN[i] == BC_PERIODIC ||
//                  SolnBlk.Grid.BCtypeN[i] == BC_CONSTANT_EXTRAPOLATION ||
//                  SolnBlk.Grid.BCtypeN[i] == BC_LINEAR_EXTRAPOLATION ||
//                  SolnBlk.Grid.BCtypeN[i] == BC_CHARACTERISTIC ||
// 		 SolnBlk.Grid.BCtypeN[i] == BC_FLAME_OUTFLOW) {
//          if (i == SolnBlk.ICl) {
//             n_pts = 5;
//             i_index[0] = i  ; j_index[0] = j-1;
//             i_index[1] = i+1; j_index[1] = j-1;
//             i_index[2] = i+1; j_index[2] = j  ;
//             i_index[3] = i  ; j_index[3] = j+1;
//             i_index[4] = i+1; j_index[4] = j+1;
//          } else if (i == SolnBlk.ICu) {
//             n_pts = 5;
//             i_index[0] = i-1; j_index[0] = j-1;
//             i_index[1] = i  ; j_index[1] = j-1;
//             i_index[2] = i-1; j_index[2] = j  ;
//             i_index[3] = i-1; j_index[3] = j+1;
//             i_index[4] = i  ; j_index[4] = j+1;
//          } else {
//             n_pts = 8;
//             i_index[0] = i-1; j_index[0] = j-1;
//             i_index[1] = i  ; j_index[1] = j-1;
//             i_index[2] = i+1; j_index[2] = j-1;
//             i_index[3] = i-1; j_index[3] = j  ;
//             i_index[4] = i+1; j_index[4] = j  ;
//             i_index[5] = i-1; j_index[5] = j+1;
//             i_index[6] = i  ; j_index[6] = j+1;
//             i_index[7] = i+1; j_index[7] = j+1;
//          } /* endif */
//       } else {
//          if (i == SolnBlk.ICl) {
//             n_pts = 3;
//             i_index[0] = i  ; j_index[0] = j-1;
//             i_index[1] = i+1; j_index[1] = j-1;
//             i_index[2] = i+1; j_index[2] = j  ;
//          } else if (i == SolnBlk.ICu) {
//             n_pts = 3;
//             i_index[0] = i-1; j_index[0] = j-1;
//             i_index[1] = i  ; j_index[1] = j-1;
//             i_index[2] = i-1; j_index[2] = j  ;
//          } else {
//             n_pts = 5;
//             i_index[0] = i-1; j_index[0] = j-1;
//             i_index[1] = i  ; j_index[1] = j-1;
//             i_index[2] = i+1; j_index[2] = j-1;
//             i_index[3] = i-1; j_index[3] = j  ;
//             i_index[4] = i+1; j_index[4] = j  ;
//          } /* endif */
//       } /* endif */
// //      } else if (((i == SolnBlk.ICl && j == SolnBlk.JCl) && 
// //                  ((SolnBlk.Grid.BCtypeW[j] == BC_REFLECTION &&
// //                    SolnBlk.Grid.BCtypeS[i] == BC_REFLECTION) ||
// //                   (SolnBlk.Grid.BCtypeW[j] == BC_BURNINGSURFACE &&
// //                    SolnBlk.Grid.BCtypeS[i] == BC_BURNINGSURFACE))) || 
// //                 ((i == SolnBlk.ICl && j == SolnBlk.JCu) && 
// //                  ((SolnBlk.Grid.BCtypeW[j] == BC_REFLECTION &&
// //                    SolnBlk.Grid.BCtypeN[i] == BC_REFLECTION) ||
// //                   (SolnBlk.Grid.BCtypeW[j] == BC_BURNINGSURFACE &&
// //                    SolnBlk.Grid.BCtypeN[i] == BC_BURNINGSURFACE))) ||
// //                 ((i == SolnBlk.ICu && j == SolnBlk.JCl) && 
// //                  ((SolnBlk.Grid.BCtypeE[j] == BC_REFLECTION &&
// //                    SolnBlk.Grid.BCtypeS[i] == BC_REFLECTION) ||
// //                   (SolnBlk.Grid.BCtypeE[j] == BC_BURNINGSURFACE &&
// //                    SolnBlk.Grid.BCtypeS[i] == BC_BURNINGSURFACE))) ||
// //                 ((i == SolnBlk.ICu && j == SolnBlk.JCu) && 
// //                  ((SolnBlk.Grid.BCtypeE[j] == BC_REFLECTION &&
// //                    SolnBlk.Grid.BCtypeN[i] == BC_REFLECTION) ||
// //                   (SolnBlk.Grid.BCtypeE[j] == BC_BURNINGSURFACE &&
// //                    SolnBlk.Grid.BCtypeN[i] == BC_BURNINGSURFACE)))) {
// //        n_pts = 8;
// //        i_index[0] = i-1; j_index[0] = j-1;
// //        i_index[1] = i  ; j_index[1] = j-1;
// //        i_index[2] = i+1; j_index[2] = j-1;
// //        i_index[3] = i-1; j_index[3] = j  ;
// //        i_index[4] = i+1; j_index[4] = j  ;
// //        i_index[5] = i-1; j_index[5] = j+1;
// //        i_index[6] = i  ; j_index[6] = j+1;
// //        i_index[7] = i+1; j_index[7] = j+1;
// //      } else if ((i == SolnBlk.ICl && j == SolnBlk.JCl) && 
// //                 (SolnBlk.Grid.BCtypeW[j] != BC_NONE &&
// //                  SolnBlk.Grid.BCtypeS[i] != BC_NONE)) {
// //        n_pts = 7;
// //        i_index[0] = i  ; j_index[0] = j-1;
// //        i_index[1] = i+1; j_index[1] = j-1;
// //        i_index[2] = i-1; j_index[2] = j  ;
// //        i_index[3] = i+1; j_index[3] = j  ;
// //        i_index[4] = i-1; j_index[4] = j+1;
// //        i_index[5] = i  ; j_index[5] = j+1;
// //        i_index[6] = i+1; j_index[6] = j+1;
// //      } else if ((i == SolnBlk.ICl && j == SolnBlk.JCu) && 
// //                 (SolnBlk.Grid.BCtypeW[j] != BC_NONE &&
// //                  SolnBlk.Grid.BCtypeN[i] != BC_NONE)) {
// //        n_pts = 7;
// //        i_index[0] = i-1; j_index[0] = j-1;
// //        i_index[1] = i  ; j_index[1] = j-1;
// //        i_index[2] = i+1; j_index[2] = j-1;
// //        i_index[3] = i-1; j_index[3] = j  ;
// //        i_index[4] = i+1; j_index[4] = j  ;
// //        i_index[5] = i  ; j_index[5] = j+1;
// //        i_index[6] = i+1; j_index[6] = j+1;
// //      } else if ((i == SolnBlk.ICu && j == SolnBlk.JCl) && 
// //                 (SolnBlk.Grid.BCtypeE[j] != BC_NONE &&
// //                  SolnBlk.Grid.BCtypeS[i] != BC_NONE)) {
// //        n_pts = 7;
// //        i_index[0] = i-1; j_index[0] = j-1;
// //        i_index[1] = i  ; j_index[1] = j-1;
// //        i_index[2] = i-1; j_index[2] = j  ;
// //        i_index[3] = i+1; j_index[3] = j  ;
// //        i_index[4] = i-1; j_index[4] = j+1;
// //        i_index[5] = i  ; j_index[5] = j+1;
// //        i_index[6] = i+1; j_index[6] = j+1;
// //      } else if ((i == SolnBlk.ICu && j == SolnBlk.JCu) && 
// //                 (SolnBlk.Grid.BCtypeE[j] != BC_NONE &&
// //                  SolnBlk.Grid.BCtypeN[i] != BC_NONE)) {
// //        n_pts = 7;
// //        i_index[0] = i-1; j_index[0] = j-1;
// //        i_index[1] = i  ; j_index[1] = j-1;
// //        i_index[2] = i+1; j_index[2] = j-1;
// //        i_index[3] = i-1; j_index[3] = j  ;
// //        i_index[4] = i+1; j_index[4] = j  ;
// //        i_index[5] = i-1; j_index[5] = j+1;
// //        i_index[6] = i  ; j_index[6] = j+1;
//     } else {
//       n_pts = 8;
//       i_index[0] = i-1; j_index[0] = j-1;
//       i_index[1] = i  ; j_index[1] = j-1;
//       i_index[2] = i+1; j_index[2] = j-1;
//       i_index[3] = i-1; j_index[3] = j  ;
//       i_index[4] = i+1; j_index[4] = j  ;
//       i_index[5] = i-1; j_index[5] = j+1;
//       i_index[6] = i  ; j_index[6] = j+1;
//       i_index[7] = i+1; j_index[7] = j+1;
//     } /* endif */

    /****************************************************************/
    //FOR VISCOUS -> CHANGED TO USE ALL 8
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
    /****************************************************************/
    
    if (n_pts > 0) {
      DUDx_ave.Vacuum();
      DUDy_ave.Vacuum();
      DxDx_ave = ZERO;
      DxDy_ave = ZERO;
      DyDy_ave = ZERO;
     
      double Sum_dx =ZERO;
      double Sum_dy =ZERO;
   
      for ( n2 = 0 ; n2 <= n_pts-1 ; ++n2 ) {
	dX = SolnBlk.Grid.Cell[ i_index[n2] ][ j_index[n2] ].Xc - 
	  SolnBlk.Grid.Cell[i][j].Xc;
// 	//************************************
// 	SolnBlk.W[ i_index[n2] ][ j_index[n2] ]= W(SolnBlk.U[ i_index[n2] ][ j_index[n2] ]);
	//***********************************
	DU = SolnBlk.W[ i_index[n2] ][ j_index[n2] ] - SolnBlk.W[i][j];
	DUDx_ave += DU*dX.x;
	DUDy_ave += DU*dX.y;
	DxDx_ave += dX.x*dX.x;
	DxDy_ave += dX.x*dX.y;
	DyDy_ave += dX.y*dX.y;
/**************************************************************************
The construction of point implicit block Jaocbians need  
the information of

   /DW\         /DW\
D |____|     D |____| 
   \Dx /        \Dy /   
 ________     ________     
   DW           DW 
It is convenient to compute them here since there are some grid geometry information 
and related computation is done here during the reconstruction by using least squares
method.
One thing is to know that those derivatives of the primitive parameter gradients are 
different for the center cell and neigbours. Details see my notes.

Xinfeng    Jan. 16, 2004 

***************************************************************************************/
	    Sum_dx += dX.x;
	    Sum_dy += dX.y;
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

	/*************************************************************/
	/* For center cell, the derivatives of gradients*/ 
	Sum_dx =Sum_dx/double(n_pts);
	Sum_dy =Sum_dy/double(n_pts);

	SolnBlk.d_dWdx_dW[i][j][0] =(-Sum_dx*DyDy_ave+Sum_dy*DxDy_ave)/
	  (DxDx_ave*DyDy_ave-DxDy_ave*DxDy_ave); 
	SolnBlk.d_dWdy_dW[i][j][0] =(-Sum_dy*DxDx_ave+Sum_dx*DxDy_ave)/
	  (DxDx_ave*DyDy_ave-DxDy_ave*DxDy_ave); 
	
  
    	if (!SolnBlk.Freeze_Limiter) {
	  for ( n = 1 ; n <= NUM_VAR_CHEM2D ; ++n ) {
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
	}//end limiter if
    } else {
        SolnBlk.dWdx[i][j].Vacuum();
        SolnBlk.dWdy[i][j].Vacuum();
        SolnBlk.phi[i][j].Vacuum(); 
	SolnBlk.d_dWdx_dW[i][j][0] =ZERO;
	SolnBlk.d_dWdy_dW[i][j][0] = ZERO;
    } /* endif */
    // 	/************************************************************/
    // Neigbour cells of (i,j) and the indces of the neibour cell's  8 sorrounding cells
    // 1 --- North (i,j+1)
    // 2 --- South (i, j-1)
    // 3 --- West  (i-1, j)
    // 4 --- East  (i+1, j)
    
    if (i < SolnBlk.ICl || i > SolnBlk.ICu ||j < SolnBlk.JCl || j > SolnBlk.JCu  ) 
      {
	for (int index_neigbor=1; index_neigbor<5; index_neigbor++)
	  {
	    SolnBlk.d_dWdx_dW[i][j][index_neigbor] =ZERO;	// GHOST CELL  MARKTHIS XINFENG
	    SolnBlk.d_dWdy_dW[i][j][index_neigbor] =ZERO;	// GHOST CELL  MARKTHIS XINFENG
	  }
	
      }else {//Catchme

	int i_i = 0;
	int j_j = 0;
	for(int kk_index=1; kk_index<5; kk_index ++)
	  {
	    if (kk_index == 1) { i_i = i; j_j = j+1; }
	    if (kk_index == 2) { i_i = i; j_j = j-1; }
	    if (kk_index == 3) { i_i = i-1; j_j = j; }
	    if (kk_index == 4) { i_i = i+1; j_j = j; }

	    i_index_neigbor[0][kk_index] = i_i-1; j_index_neigbor[0][kk_index] = j_j-1;
	    i_index_neigbor[1][kk_index] = i_i  ; j_index_neigbor[1][kk_index] = j_j-1;
	    i_index_neigbor[2][kk_index] = i_i+1; j_index_neigbor[2][kk_index] = j_j-1;
	    i_index_neigbor[3][kk_index] = i_i-1; j_index_neigbor[3][kk_index] = j_j  ;
	    i_index_neigbor[4][kk_index] = i_i+1; j_index_neigbor[4][kk_index] = j_j  ;
	    i_index_neigbor[5][kk_index] = i_i-1; j_index_neigbor[5][kk_index] = j_j+1;
	    i_index_neigbor[6][kk_index] = i_i  ; j_index_neigbor[6][kk_index] = j_j+1;
	    i_index_neigbor[7][kk_index] = i_i+1; j_index_neigbor[7][kk_index] = j_j+1;

	  }  

	/************************************************************/

	/* For Neigbour cell, the derivatives of gradients*/ 
	/*************************************************************/


	int I_Iid = 0, J_Jid = 0;
	int nn_pts = 8;
	for (int index_neigbor=1; index_neigbor<5; index_neigbor++)
	  {
	    DWDx_ave.Vacuum();
	    DWDy_ave.Vacuum();
	    dxdx_neigbor = ZERO;
	    dxdy_neigbor = ZERO;
	    dydy_neigbor = ZERO;
	    
	    double Sum_dx =ZERO;
	    double Sum_dy =ZERO;
	    
	    if (index_neigbor == 1) { I_Iid = i; J_Jid = j+1; }
	    if (index_neigbor == 2) { I_Iid = i; J_Jid = j-1; }
	    if (index_neigbor == 3) { I_Iid = i-1; J_Jid = j; }
	    if (index_neigbor == 4) { I_Iid = i+1; J_Jid = j; }
   
	
	    for ( int n_n = 0 ; n_n <= nn_pts-1 ; ++n_n ) {
	      
	      int I_cell = i_index_neigbor[n_n][index_neigbor];
	      int J_cell = j_index_neigbor[n_n][index_neigbor];

	      
	      dX_neigbor = SolnBlk.Grid.Cell[I_cell][J_cell].Xc- SolnBlk.Grid.Cell[I_Iid][J_Jid].Xc;
// 	      DW =  SolnBlk.W[I_cell][J_cell] - SolnBlk.W[I_Iid][J_Jid];
// 	      DWDx_ave += DW*dX_neigbor.x;
// 	      DWDy_ave += DW*dX_neigbor.y;
	      dxdx_neigbor += dX_neigbor.x*dX_neigbor.x;
	      dxdy_neigbor += dX_neigbor.x*dX_neigbor.y;
	      dydy_neigbor += dX_neigbor.y*dX_neigbor.y;
	    
	      
	    } /* end inside for */
	  

	    
// 	    DWDx_ave = DWDx_ave/double(nn_pts);
// 	    DWDy_ave = DWDy_ave/double(nn_pts);
	    dxdx_neigbor = dxdx_neigbor/double(nn_pts);
	    dxdy_neigbor = dxdy_neigbor/double(nn_pts);
	    dydy_neigbor = dydy_neigbor/double(nn_pts);
	    
	    
	    Sum_dx = SolnBlk.Grid.Cell[i][j].Xc.x -SolnBlk.Grid.Cell[I_Iid][J_Jid].Xc.x;
	    Sum_dy = SolnBlk.Grid.Cell[i][j].Xc.y -SolnBlk.Grid.Cell[I_Iid][J_Jid].Xc.y;
	    
	    Sum_dx =Sum_dx/double(nn_pts);
	    Sum_dy =Sum_dy/double(nn_pts);
	 	 

	    SolnBlk.d_dWdx_dW[i][j][index_neigbor] =(Sum_dx*dydy_neigbor-Sum_dy*dxdy_neigbor)/
	      (dxdx_neigbor*dydy_neigbor-dxdy_neigbor*dxdy_neigbor); 
	    
	    SolnBlk.d_dWdy_dW[i][j][index_neigbor] =(Sum_dy*dxdx_neigbor-Sum_dx*dxdy_neigbor)/
	      (dxdx_neigbor*dydy_neigbor-dxdy_neigbor*dxdy_neigbor); 

// 	    SolnBlk.d_dWdx_dW[i][j][index_neigbor] =(Sum_dx*dydy_neigbor-Sum_dy*dxdy_neigbor)/
// 	      (dxdx_neigbor*dydy_neigbor-dxdy_neigbor*dxdy_neigbor); 
		  
// 	    SolnBlk.d_dWdy_dW[i][j][index_neigbor] =(Sum_dy*dxdx_neigbor-Sum_dx*dxdy_neigbor)/
// 	      (dxdx_neigbor*dydy_neigbor-dxdy_neigbor*dxdy_neigbor); 
	   //  if(index_neigbor ==1)
// 	      {
// 		cout<<" Noth Sum_dy ="<<Sum_dy<<" Noth Sum_dx ="<<Sum_dx<<"   "<<dxdx_neigbor<<endl;
// 		cout<<"dxdx_neigbor =   "<<dxdx_neigbor<<endl;
// 		cout<<"dydy_neigbor =   "<<dydy_neigbor<<endl;
// 		cout<<"dxdy_neigbor =   "<<dxdy_neigbor<<endl;

// 	      }

	  }
	  
      }
	/************************************************************/


    
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
void Linear_Reconstruction_LeastSquares_2(Chem2D_Quad_Block &SolnBlk,
				          const int i, 
                                          const int j,
                                          const int Limiter) {

  cout<<"\n IN LEAST SQUARES 2!!!!!!!!! \n";

    int n, n2, n_pts, i_index[8], j_index[8];
    double u0Min, u0Max, uQuad[4], phi;
    double DxDx_ave, DxDy_ave, DyDy_ave;
    Vector2D dX;
    Chem2D_pState DU, DUDx_ave, DUDy_ave;
   
    int num=0;
    int NUM_VAR_CHEM2D = SolnBlk.NumVar();

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
                 SolnBlk.Grid.BCtypeW[j] == BC_CHARACTERISTIC ||
		 SolnBlk.Grid.BCtypeW[j] == BC_FLAME_OUTFLOW) {
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
                 SolnBlk.Grid.BCtypeE[j] == BC_CHARACTERISTIC ||
		 SolnBlk.Grid.BCtypeE[j] == BC_FLAME_OUTFLOW) {
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
                 SolnBlk.Grid.BCtypeS[i] == BC_CHARACTERISTIC ||
		 SolnBlk.Grid.BCtypeS[i] == BC_FLAME_OUTFLOW) {
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
                 SolnBlk.Grid.BCtypeN[i] == BC_CHARACTERISTIC ||
		 SolnBlk.Grid.BCtypeN[i] == BC_FLAME_OUTFLOW) {
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
        DUDx_ave.Vacuum();
        DUDy_ave.Vacuum();
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
	  for ( n = 1 ; n <= NUM_VAR_CHEM2D ; ++n ) {
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
	} //end limiter if

    } else {
      SolnBlk.dWdx[i][j].Vacuum();
      SolnBlk.dWdy[i][j].Vacuum();
      SolnBlk.phi[i][j].Vacuum();
    } /* endif */
   
}

//Diamond Path
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
void Linear_Reconstruction_LeastSquares_Diamond(Chem2D_Quad_Block &SolnBlk,
				        const int Limiter) {
 
    /* Carry out the limited solution reconstruction in
       each cell of the computational mesh. */

    for (int j  = SolnBlk.JCl-SolnBlk.Nghost+1 ; j <= SolnBlk.JCu+SolnBlk.Nghost-1 ; ++j ) {
       for (int i = SolnBlk.ICl-SolnBlk.Nghost+1 ; i <= SolnBlk.ICu+SolnBlk.Nghost-1 ; ++i ) {
	   Linear_Reconstruction_LeastSquares_Diamond(SolnBlk, i, j, Limiter);
   	
       } /* endfor */
    } /* endfor */

}
//Diamond path 
void Linear_Reconstruction_LeastSquares_Diamond(Chem2D_Quad_Block &SolnBlk,
						const int i, 
						const int j,
						const int Limiter) {
  int n, n2, n_pts, i_index[8], j_index[8];
  double u0Min, u0Max, uQuad[4], phi;
  double DxDx_ave, DxDy_ave, DyDy_ave;
  double area[4];

  Vector2D deltaX, dX[4];  
  Chem2D_pState DU[4], DUDx_ave, DUDy_ave;
  Chem2D_pState QuadraturePoint_N, QuadraturePoint_E, 
    QuadraturePoint_S, QuadraturePoint_W;  
  Chem2D_pState TopVertex, BottomVertex; 
  
  Chem2D_pState Temp;
  //For the derivatives of gradients purpose
  // i_index_neigbor [8][1]  1 -- North , 8 -- its 8 sorrounding cells
  // i_index_neigbor [8][2]  2 -- South , 8 -- its 8 sorrounding cells
  // i_index_neigbor [8][3]  3 -- South , 8 -- its 8 sorrounding cells
  // i_index_neigbor [8][4]  4 -- South , 8 -- its 8 sorrounding cells
  
  int i_index_neigbor [8][5], j_index_neigbor[8][5];
  double dxdx_neigbor, dxdy_neigbor, dydy_neigbor;
  Vector2D dX_neigbor;
  Chem2D_pState DW, DWDx_ave, DWDy_ave;
  Chem2D_pState W_node;
  
  int num=0;
  int NUM_VAR_CHEM2D = SolnBlk.NumVar();
  
  /****************************************************************/
  //  * A least squares       *
  //  * approach is used in the evaluation of the unlimited  *
  //  * solution gradients on cell faces that is peformed on a diamond path
  // For cell (i,j), in order to reconstruct the gradients on cell faces, 
  // information of four points are needed
  // centeroid (i+1, j) and vertices (i+1, j) and (i+1, j+1)
  // By convention: the quadrature point (mid-edge) (i, j)  
  
  if (i == SolnBlk.ICl-SolnBlk.Nghost || i == SolnBlk.ICu+SolnBlk.Nghost ||
      j == SolnBlk.JCl-SolnBlk.Nghost || j == SolnBlk.JCu+SolnBlk.Nghost) {
    n_pts = 0;
    
    
  } else {
    n_pts = 4;  
    i_index[0] = i  ; j_index[0] = j;
    i_index[1] = i+1; j_index[1] = j  ;
    i_index[2] = i  ; j_index[2] = j+1;
    i_index[3] = i+1; j_index[3] = j+1;
    
    }    
  
  //The following notations are used for viscous Jacobians
  //For this diamond path, d_dWdx_dW[i][j][index] is as follows
  //For d_dWdx_dW[i][j][1] --north face
  //For d_dWdx_dW[i][j][2] --east face
  //For d_dWdx_dW[i][j][3] --west face
  //For d_dWdx_dW[i][j][4] --south face
  
  // Derivatives of Gradients ...
  double dWnNWdWc = ZERO;
  double dWnNEdWc = ZERO;
  double dWnSWdWc = ZERO;
  double dWnSEdWc = ZERO;
  
  string Orient;
  
  Orient = "NW";
  dWnNWdWc = SolnBlk.dWn_dWc(i, j, Orient);
  Orient = "NE";
  dWnNEdWc = SolnBlk.dWn_dWc(i, j, Orient);
  Orient = "SW";
  dWnSWdWc = SolnBlk.dWn_dWc(i, j, Orient);
  Orient = "SE";
  dWnSEdWc = SolnBlk. dWn_dWc(i, j, Orient);
  
  
  if (n_pts > 0) {
    
    /*Formulate the gradients of primitive parameters on the north face of cell (i, j)*/
    DUDx_ave.Vacuum();
    DUDy_ave.Vacuum();
    DxDx_ave = ZERO;
    DxDy_ave = ZERO;
    DyDy_ave = ZERO;
    
   
    TopVertex = SolnBlk.Wn(i,j+1);
    BottomVertex =  SolnBlk.Wn(i+1,j+1);
    QuadraturePoint_N = HALF*( TopVertex + BottomVertex);
      
    //Left state cell (i,j)
    dX[0] = SolnBlk.Grid.Cell[i][j].Xc - SolnBlk.Grid.xfaceN(i,j);
    DU[0] = SolnBlk.W[i][j] - QuadraturePoint_N;
    //rigth state cell (i, j+1)
    dX[1] = SolnBlk.Grid.Cell[i][j+1].Xc - SolnBlk.Grid.xfaceN(i,j);
    DU[1] = SolnBlk.W[i][j+1] - QuadraturePoint_N;
    //top vertex  (i, j+1)
    dX[2] = SolnBlk.Grid.Node[i][j+1].X -  SolnBlk.Grid.xfaceN(i,j);
    DU[2] = TopVertex - QuadraturePoint_N;
    //bottom vertex (i+1,j+1) 
    dX[3] = SolnBlk.Grid.Node[i+1][j+1].X -  SolnBlk.Grid.xfaceN(i,j);
    DU[3] = BottomVertex - QuadraturePoint_N;

    //The calculation of this area is used to weight the gradient at the cell center      
    area[0] = HALF*(SolnBlk.Grid.Node[i+1][j+1].X - SolnBlk.Grid.Node[i][j+1].X )^
      (SolnBlk.Grid.xfaceN(i,j)- SolnBlk.Grid.Cell[i][j].Xc);
    
    
    
    for ( n2 = 0 ; n2 <= n_pts-1 ; ++n2 ) {
      DUDx_ave += DU[n2]*dX[n2].x;
      DUDy_ave += DU[n2]*dX[n2].y;
      DxDx_ave += dX[n2].x*dX[n2].x;
      DxDy_ave += dX[n2].x*dX[n2].y;
      DyDy_ave += dX[n2].y*dX[n2].y;
      
    } /* endfor */
    
    DUDx_ave = DUDx_ave/double(n_pts);
    DUDy_ave = DUDy_ave/double(n_pts);
    DxDx_ave = DxDx_ave/double(n_pts);
    DxDy_ave = DxDy_ave/double(n_pts);
    DyDy_ave = DyDy_ave/double(n_pts);
    
    SolnBlk.dWdx_faceN[i][j] = (DUDx_ave*DyDy_ave-DUDy_ave*DxDy_ave)/
      (DxDx_ave*DyDy_ave-DxDy_ave*DxDy_ave);
    SolnBlk.dWdy_faceN[i][j] = (DUDy_ave*DxDx_ave-DUDx_ave*DxDy_ave)/
      (DxDx_ave*DyDy_ave-DxDy_ave*DxDy_ave);

    
    SolnBlk.d_dWdx_dW[i][j][1] =((dX[0].x+ dX[2].x*dWnNWdWc +  dX[3].x*dWnNEdWc )/double(n_pts)*DyDy_ave 
				 -(dX[0].y+ dX[2].y*dWnNWdWc +  dX[3].y*dWnNEdWc )/double(n_pts)*DxDy_ave)
      /(DxDx_ave*DyDy_ave-DxDy_ave*DxDy_ave) ;
    
    
    SolnBlk.d_dWdy_dW[i][j][1] =((dX[0].y+ dX[2].y*dWnNWdWc + dX[3].y*dWnNEdWc )/double(n_pts)*DxDx_ave
				 -(dX[0].x+ dX[2].x*dWnNWdWc +  dX[3].x*dWnNEdWc )/double(n_pts) *DxDy_ave)/
      (DxDx_ave*DyDy_ave-DxDy_ave*DxDy_ave);
    
    /*Formulate the gradients of primitive parameters on the east face of cell (i, j)*/
    DUDx_ave.Vacuum();
    DUDy_ave.Vacuum();
    DxDx_ave = ZERO;
    DxDy_ave = ZERO;
    DyDy_ave = ZERO;
      
    //needs to assign topvertex and bottomvertex information
    TopVertex = SolnBlk.Wn(i+1,j+1);
    BottomVertex =  SolnBlk.Wn(i+1,j);
    QuadraturePoint_E = HALF*( TopVertex + BottomVertex);
    
    //Left state cell (i,j)
    dX[0] = SolnBlk.Grid.Cell[i][j].Xc - SolnBlk.Grid.xfaceE(i, j);
    DU[0] = SolnBlk.W[i][ j] - QuadraturePoint_E;
    //rigth state cell (i+1, j)
    dX[1] = SolnBlk.Grid.Cell[i+1][j].Xc - SolnBlk.Grid.xfaceE(i,j);
    DU[1] = SolnBlk.W[i+1][ j] - QuadraturePoint_E;
    //top vertex  (i+1, j+1)
    dX[2] = SolnBlk.Grid.Node[i+1][j+1].X - SolnBlk.Grid.xfaceE(i,j);
    DU[2] = TopVertex - QuadraturePoint_E;
    //bottom vertex (i+1,j) 
    dX[3] = SolnBlk.Grid.Node[i+1][j].X - SolnBlk.Grid.xfaceE(i,j);
    DU[3] = BottomVertex - QuadraturePoint_E;
    
      //The calculation of this area is used to weight the gradient at the cell center      
    area[1] = HALF*(SolnBlk.Grid.Node[i+1][j+1].X - SolnBlk.Grid.Node[i+1][j].X )^
      (SolnBlk.Grid.Cell[i][j].Xc - SolnBlk.Grid.xfaceE(i,j));
    
    for ( n2 = 0 ; n2 <= n_pts-1 ; ++n2 ) {
      DUDx_ave += DU[n2]*dX[n2].x;
      DUDy_ave += DU[n2]*dX[n2].y;
      DxDx_ave += dX[n2].x*dX[n2].x;
      DxDy_ave += dX[n2].x*dX[n2].y;
      DyDy_ave += dX[n2].y*dX[n2].y;
      
    } /* endfor */
    
    DUDx_ave = DUDx_ave/double(n_pts);
    DUDy_ave = DUDy_ave/double(n_pts);
    DxDx_ave = DxDx_ave/double(n_pts);
    DxDy_ave = DxDy_ave/double(n_pts);
    DyDy_ave = DyDy_ave/double(n_pts);
    
    SolnBlk.dWdx_faceE[i][j] = (DUDx_ave*DyDy_ave-DUDy_ave*DxDy_ave)/
      (DxDx_ave*DyDy_ave-DxDy_ave*DxDy_ave);
    SolnBlk.dWdy_faceE[i][j] = (DUDy_ave*DxDx_ave-DUDx_ave*DxDy_ave)/
      (DxDx_ave*DyDy_ave-DxDy_ave*DxDy_ave);
    
    
    SolnBlk.d_dWdx_dW[i][j][2] =((dX[0].x+ dX[2].x*dWnNEdWc +  dX[3].x*dWnSEdWc )/double(n_pts)*DyDy_ave 
				 -(dX[0].y+ dX[2].y*dWnNEdWc +  dX[3].y*dWnSEdWc )/double(n_pts)*DxDy_ave)
      /(DxDx_ave*DyDy_ave-DxDy_ave*DxDy_ave) ;
    
    SolnBlk.d_dWdy_dW[i][j][2] =((dX[0].y+ dX[2].y*dWnNEdWc + dX[3].y*dWnSEdWc )/double(n_pts)*DxDx_ave
				 -(dX[0].x+ dX[2].x*dWnNEdWc +  dX[3].x*dWnSEdWc )/double(n_pts)*DxDy_ave)/
      (DxDx_ave*DyDy_ave-DxDy_ave*DxDy_ave);
    
    /*Formulate the gradients of primitive parameters on the west face of cell (i, j)*/
    DUDx_ave.Vacuum();
    DUDy_ave.Vacuum();
    DxDx_ave = ZERO;
    DxDy_ave = ZERO;
    DyDy_ave = ZERO;
  
    TopVertex = SolnBlk.Wn(i,j);
    BottomVertex =  SolnBlk.Wn(i,j+1);
    QuadraturePoint_W = HALF*( TopVertex + BottomVertex);
    
    //Left state cell (i,j)
    dX[0] = SolnBlk.Grid.Cell[i][j].Xc - SolnBlk.Grid.xfaceW(i,j);
    DU[0] = SolnBlk.W[i][ j] - QuadraturePoint_W;
    //rigth state cell (i-1, j)
    dX[1] = SolnBlk.Grid.Cell[i-1][j].Xc - SolnBlk.Grid.xfaceW(i,j);
    DU[1] = SolnBlk.W[i-1][j] - QuadraturePoint_W;
    //top vertex  (i, j)
    dX[2] = SolnBlk.Grid.Node[i][j].X - SolnBlk.Grid.xfaceW(i,j);
    DU[2] = TopVertex - QuadraturePoint_W;
    //bottom vertex (i,j+1) 
    dX[3] = SolnBlk.Grid.Node[i][j+1].X - SolnBlk.Grid.xfaceW(i,j);
    DU[3] = BottomVertex - QuadraturePoint_W;
      
      //The calculation of this area is used to weight the gradient at the cell center      
    area[2] = HALF*(SolnBlk.Grid.Node[i][j+1].X - SolnBlk.Grid.Node[i][j].X )^
      ( SolnBlk.Grid.xfaceW(i, j) - SolnBlk.Grid.Cell[i][j].Xc );
    
    for ( n2 = 0 ; n2 <= n_pts-1 ; ++n2 ) {
      DUDx_ave += DU[n2]*dX[n2].x;
      DUDy_ave += DU[n2]*dX[n2].y;
      DxDx_ave += dX[n2].x*dX[n2].x;
      DxDy_ave += dX[n2].x*dX[n2].y;
      DyDy_ave += dX[n2].y*dX[n2].y;
      
    } /* endfor */
    
    DUDx_ave = DUDx_ave/double(n_pts);
    DUDy_ave = DUDy_ave/double(n_pts);
    DxDx_ave = DxDx_ave/double(n_pts);
    DxDy_ave = DxDy_ave/double(n_pts);
    DyDy_ave = DyDy_ave/double(n_pts);
    
    SolnBlk.dWdx_faceW[i][j] = (DUDx_ave*DyDy_ave-DUDy_ave*DxDy_ave)/
      (DxDx_ave*DyDy_ave-DxDy_ave*DxDy_ave);
    SolnBlk.dWdy_faceW[i][j] = (DUDy_ave*DxDx_ave-DUDx_ave*DxDy_ave)/
      (DxDx_ave*DyDy_ave-DxDy_ave*DxDy_ave);
    
    //Following previous index convenction 3 is for west face
    SolnBlk.d_dWdx_dW[i][j][3] =((dX[0].x+ dX[2].x*dWnSWdWc +  dX[3].x*dWnNWdWc )/double(n_pts)*DyDy_ave 
				 -(dX[0].y+ dX[2].y*dWnSWdWc +  dX[3].y*dWnNWdWc )/double(n_pts)*DxDy_ave)
      /(DxDx_ave*DyDy_ave-DxDy_ave*DxDy_ave) ;
    
    SolnBlk.d_dWdy_dW[i][j][3] =((dX[0].y+ dX[2].y*dWnSWdWc + dX[3].y*dWnNWdWc )/double(n_pts)*DxDx_ave
				 -(dX[0].x+ dX[2].x*dWnSWdWc +  dX[3].x*dWnNWdWc )/double(n_pts)*DxDy_ave)/
      (DxDx_ave*DyDy_ave-DxDy_ave*DxDy_ave);
    
    /*Formulate the gradients of primitive parameters on the south face of cell (i, j)*/
    DUDx_ave.Vacuum();
    DUDy_ave.Vacuum();
    DxDx_ave = ZERO;
    DxDy_ave = ZERO;
    DyDy_ave = ZERO;
    //needs to assign topvertex and bottomvertex information
    //   TopVertex = SolnBlk.WnSE(i,j);
    //       BottomVertex =  SolnBlk.WnSW(i,j);
    TopVertex = SolnBlk.Wn(i+1,j);
    BottomVertex =  SolnBlk.Wn(i,j);
    QuadraturePoint_W = HALF*(TopVertex + BottomVertex);
    
    //Left state cell (i,j)
    dX[0] = SolnBlk.Grid.Cell[i][j].Xc - SolnBlk.Grid.xfaceS(i, j);
    DU[0] = SolnBlk.W[i][j] - QuadraturePoint_S;
    //rigth state cell (i, j-1)
    dX[1] = SolnBlk.Grid.Cell[i][j-1].Xc - SolnBlk.Grid.xfaceS(i, j);
    DU[1] = SolnBlk.W[i][j-1] - QuadraturePoint_S;
    //top vertex  (i+1, j)
    dX[2] = SolnBlk.Grid.Node[i+1][j].X - SolnBlk.Grid.xfaceS(i, j);
    DU[2] = TopVertex - QuadraturePoint_S;
    //bottom vertex (i,j) 
    dX[3] = SolnBlk.Grid.Node[i][j].X - SolnBlk.Grid.xfaceS(i, j);
    DU[3] = BottomVertex - QuadraturePoint_S;
    
    //The calculation of this area is used to weight the gradient at the cell center      
    area[3] = HALF*(SolnBlk.Grid.Node[i+1][j].X - SolnBlk.Grid.Node[i][j].X )^
      (SolnBlk.Grid.Cell[i][j].Xc - SolnBlk.Grid.xfaceS(i,j) );
    
    for ( n2 = 0 ; n2 <= n_pts-1 ; ++n2 ) {
      DUDx_ave += DU[n2]*dX[n2].x;
      DUDy_ave += DU[n2]*dX[n2].y;
      DxDx_ave += dX[n2].x*dX[n2].x;
      DxDy_ave += dX[n2].x*dX[n2].y;
      DyDy_ave += dX[n2].y*dX[n2].y;
      
    } /* endfor */
    
    DUDx_ave = DUDx_ave/double(n_pts);
    DUDy_ave = DUDy_ave/double(n_pts);
    DxDx_ave = DxDx_ave/double(n_pts);
    DxDy_ave = DxDy_ave/double(n_pts);
    DyDy_ave = DyDy_ave/double(n_pts);
    
    SolnBlk.dWdx_faceS[i][j] = (DUDx_ave*DyDy_ave-DUDy_ave*DxDy_ave)/
      (DxDx_ave*DyDy_ave-DxDy_ave*DxDy_ave);
    SolnBlk.dWdy_faceS[i][j] = (DUDy_ave*DxDx_ave-DUDx_ave*DxDy_ave)/
      (DxDx_ave*DyDy_ave-DxDy_ave*DxDy_ave);
    
    //Following previous index convenction 2 is for south face
    SolnBlk.d_dWdx_dW[i][j][4] =((dX[0].x+ dX[2].x*dWnSEdWc +  dX[3].x*dWnSWdWc )/double(n_pts)*DyDy_ave 
				 -(dX[0].y+ dX[2].y*dWnSEdWc +  dX[3].y*dWnSWdWc )/double(n_pts)*DxDy_ave)
      /(DxDx_ave*DyDy_ave-DxDy_ave*DxDy_ave) ;
    
    SolnBlk.d_dWdy_dW[i][j][4] =((dX[0].y+ dX[2].y*dWnSEdWc + dX[3].y*dWnSWdWc )/double(n_pts)*DxDx_ave
				 -(dX[0].x+ dX[2].x*dWnSEdWc +  dX[3].x*dWnSWdWc )/double(n_pts)*DxDy_ave)/
      (DxDx_ave*DyDy_ave-DxDy_ave*DxDy_ave);
      

    //Area weighted gradients at cell centers
    SolnBlk.dWdx[i][j] = SolnBlk.dWdx_faceN[i][j]*area[0]+SolnBlk.dWdx_faceE[i][j]*area[1] +
      SolnBlk.dWdx_faceW[i][j]*area[2] + SolnBlk.dWdx_faceS[i][j]*area[3]; 
    
    SolnBlk.dWdx[i][j] = SolnBlk.dWdx[i][j]/SolnBlk.Grid.Cell[i][j].A; 
    
    SolnBlk.dWdy[i][j] = SolnBlk.dWdy_faceN[i][j]*area[0]+ SolnBlk.dWdy_faceE[i][j]*area[1] +
      SolnBlk.dWdy_faceW[i][j]*area[2] + SolnBlk.dWdy_faceS[i][j]*area[3]; 
    
    SolnBlk.dWdy[i][j] = SolnBlk.dWdy[i][j]/SolnBlk.Grid.Cell[i][j].A; 
    
    //For source Jacobians use
    SolnBlk.d_dWdx_dW[i][j][0] = SolnBlk.d_dWdx_dW[i][j][1]*area[0]+SolnBlk.d_dWdx_dW[i][j][2]*area[1] +
      SolnBlk.d_dWdx_dW[i][j][3]*area[2] + SolnBlk.d_dWdx_dW[i][j][4]*area[3];  
    SolnBlk.d_dWdx_dW[i][j][0] = SolnBlk.d_dWdx_dW[i][j][0]/SolnBlk.Grid.Cell[i][j].A;  
    
    SolnBlk.d_dWdy_dW[i][j][0] = SolnBlk.d_dWdy_dW[i][j][1]*area[0]+SolnBlk.d_dWdy_dW[i][j][2]*area[1] +
      SolnBlk.d_dWdy_dW[i][j][3]*area[2] + SolnBlk.d_dWdy_dW[i][j][4]*area[3];  
    SolnBlk.d_dWdy_dW[i][j][0] = SolnBlk.d_dWdy_dW[i][j][0]/SolnBlk.Grid.Cell[i][j].A;  
    
    
    if (!SolnBlk.Freeze_Limiter) {
      for ( n = 1 ; n <= NUM_VAR_CHEM2D ; ++n ) {
	u0Min = SolnBlk.W[i][j][n];
	u0Max = u0Min;
	for(n2 = 0 ; n2 <= n_pts-1 ; ++n2){
	  Temp =  SolnBlk.Wn(i_index[n2],j_index[n2]);
	  u0Min = min(u0Min, Temp[n] );
	  u0Max = max(u0Max, Temp[n]);
	  
	}
	
	deltaX = SolnBlk.Grid.xfaceE(i, j)-SolnBlk.Grid.Cell[i][j].Xc;
	uQuad[0] = SolnBlk.W[i][j][n] + 
	  SolnBlk.dWdx[i][j][n]*deltaX.x +
	  SolnBlk.dWdy[i][j][n]*deltaX.y ;
	deltaX = SolnBlk.Grid.xfaceW(i, j)-SolnBlk.Grid.Cell[i][j].Xc;
	uQuad[1] = SolnBlk.W[i][j][n] + 
	  SolnBlk.dWdx[i][j][n]*deltaX.x +
	  SolnBlk.dWdy[i][j][n]*deltaX.y ;
	deltaX = SolnBlk.Grid.xfaceN(i, j)-SolnBlk.Grid.Cell[i][j].Xc;
	uQuad[2] = SolnBlk.W[i][j][n] + 
	  SolnBlk.dWdx[i][j][n]*deltaX.x +
	  SolnBlk.dWdy[i][j][n]*deltaX.y ;
	deltaX = SolnBlk.Grid.xfaceS(i, j)-SolnBlk.Grid.Cell[i][j].Xc;
	uQuad[3] = SolnBlk.W[i][j][n] + 
	  SolnBlk.dWdx[i][j][n]*deltaX.x +
	  SolnBlk.dWdy[i][j][n]*deltaX.y ;
	  
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
    }//end limiter if 
  } else {
    SolnBlk.dWdx[i][j].Vacuum();
    SolnBlk.dWdy[i][j].Vacuum();
    SolnBlk.phi[i][j].Vacuum(); 
    SolnBlk.d_dWdx_dW[i][j][0] = ZERO;
    SolnBlk.d_dWdy_dW[i][j][0] = ZERO;
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
void Linear_Reconstruction_LeastSquares(Chem2D_Quad_Block &SolnBlk,
				        const int Limiter) {
 
    /* Carry out the limited solution reconstruction in
       each cell of the computational mesh. */

    for (int j  = SolnBlk.JCl-SolnBlk.Nghost+1 ; j <= SolnBlk.JCu+SolnBlk.Nghost-1 ; ++j ) {
       for (int i = SolnBlk.ICl-SolnBlk.Nghost+1 ; i <= SolnBlk.ICu+SolnBlk.Nghost-1 ; ++i ) {
	   Linear_Reconstruction_LeastSquares(SolnBlk, i, j, Limiter);
           //Linear_Reconstruction_LeastSquares_2(SolnBlk, i, j, Limiter);
       } /* endfor */
    } /* endfor */

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
void Residual_Smoothing(Chem2D_Quad_Block &SolnBlk,
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
 * Routine: Distance_to_Wall                            *
 *                                                      *
 * Determines the distance to the nearest wall for      *
 * given cell centered location.                        *
 *                                                      *
 ********************************************************/
double Distance_to_Wall(Chem2D_Quad_Block &SolnBlk,
                        const Vector2D X_cell) {

    int i, j;
    double y_wall;

    // Initialize y_wall.
    y_wall = 1e70;

    // Check West boundary.
    for ( j = SolnBlk.JCl; j <= SolnBlk.JCu; ++j ) {
       if (SolnBlk.Grid.BCtypeW[j] == BC_WALL_VISCOUS ||
           SolnBlk.Grid.BCtypeW[j] == BC_NO_SLIP  ||
           SolnBlk.Grid.BCtypeW[j] == BC_WALL_VISCOUS_ISOTHERMAL ||
           SolnBlk.Grid.BCtypeW[j] == BC_WALL_VISCOUS_HEATFLUX ||
           SolnBlk.Grid.BCtypeW[j] == BC_ADIABATIC_WALL) {
          y_wall = min(y_wall, abs(SolnBlk.Grid.xfaceW(SolnBlk.ICl,j)-X_cell));

       } /* endif */
    } /* endfor */ 

    // Check East boundary.
    for ( j = SolnBlk.JCl; j <= SolnBlk.JCu; ++j ) {
       if (SolnBlk.Grid.BCtypeE[j] == BC_WALL_VISCOUS ||
           SolnBlk.Grid.BCtypeE[j] == BC_NO_SLIP  ||
           SolnBlk.Grid.BCtypeE[j] == BC_WALL_VISCOUS_ISOTHERMAL ||
           SolnBlk.Grid.BCtypeE[j] == BC_WALL_VISCOUS_HEATFLUX ||
           SolnBlk.Grid.BCtypeE[j] == BC_ADIABATIC_WALL) {
          y_wall = min(y_wall, abs(SolnBlk.Grid.xfaceE(SolnBlk.ICu,j)-X_cell));

       } /* endif */
    } /* endfor */ 

    // Check South boundary.
    for ( i = SolnBlk.ICl; i <= SolnBlk.ICu ; ++i ) {
       if (SolnBlk.Grid.BCtypeS[i] == BC_WALL_VISCOUS ||
           SolnBlk.Grid.BCtypeS[i] == BC_NO_SLIP  ||
           SolnBlk.Grid.BCtypeS[i] == BC_WALL_VISCOUS_ISOTHERMAL ||
           SolnBlk.Grid.BCtypeS[i] == BC_WALL_VISCOUS_HEATFLUX ||
           SolnBlk.Grid.BCtypeS[i] == BC_ADIABATIC_WALL) {
          y_wall = min(y_wall, abs(SolnBlk.Grid.xfaceS(i,SolnBlk.JCl)-X_cell));
       } /* endif */
    } /* endfor */

    // Check North boundary.
    for ( i = SolnBlk.ICl; i <= SolnBlk.ICu ; ++i ) {
       if (SolnBlk.Grid.BCtypeN[i] == BC_WALL_VISCOUS ||
           SolnBlk.Grid.BCtypeN[i] == BC_NO_SLIP  ||
           SolnBlk.Grid.BCtypeN[i] == BC_WALL_VISCOUS_ISOTHERMAL ||
           SolnBlk.Grid.BCtypeN[i] == BC_WALL_VISCOUS_HEATFLUX ||
           SolnBlk.Grid.BCtypeN[i] == BC_ADIABATIC_WALL) {
          y_wall = min(y_wall, abs(SolnBlk.Grid.xfaceN(i,SolnBlk.JCu)-X_cell));
       } /* endif */
    } /* endfor */

    return (y_wall);

}

/********************************************************
 * Routine: Calculate_Refinement_Criteria               *
 *                                                      *
 * Calculate refinement criteria for the solution       *
 * block.                                               *
 *                                                      *
 ********************************************************/
void Calculate_Refinement_Criteria(double *refinement_criteria,
				   Chem2D_Input_Parameters &IP,
                                   int &number_refinement_criteria,
                                   Chem2D_Quad_Block &SolnBlk) {

    int i, j;

    double grad_rho_x, grad_rho_y, grad_rho_abs, grad_rho_criteria, grad_rho_criteria_max,
      div_V, div_V_criteria, div_V_criteria_max,
      curl_V_z, curl_V_abs, curl_V_criteria, curl_V_criteria_max,
      grad_Temp_x, grad_Temp_y, grad_Temp_abs, grad_Temp_criteria, grad_Temp_criteria_max,
      grad_CH4_x, grad_CH4_y, grad_CH4_abs, grad_CH4_criteria, grad_CH4_criteria_max,
      grad_CO2_x, grad_CO2_y, grad_CO2_abs, grad_CO2_criteria, grad_CO2_criteria_max;

    /* Set the number of refinement criteria to be used (3):
       (1) refinement criteria #1 based on the gradient of the density field;
       (2) refinement criteria #2 based on the divergence of the velocity vector;
       (3) refinement criteria #3 based on the curl of the velocity vector. 
       (4) refinement criteria #4 based on the gradient of Temperature
       (5) refinement criteria #5 based on the gradient of CH4 mass fraction
       (6) refinement criteria #6 based on the gradient of CO2 mass fraction
    */
    number_refinement_criteria = 3;

    /* Initialize the refinement criteria for the block. */

    grad_rho_criteria_max = ZERO;
    div_V_criteria_max = ZERO;
    curl_V_criteria_max = ZERO;
    grad_Temp_criteria_max = ZERO;
    grad_CH4_criteria_max = ZERO;      
    grad_CO2_criteria_max = ZERO;

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
             grad_rho_x = SolnBlk.dWdx[i][j].rho;
             grad_rho_y = SolnBlk.dWdy[i][j].rho;
             grad_rho_abs = sqrt(sqr(grad_rho_x) + sqr(grad_rho_y));
             grad_rho_criteria = sqrt(SolnBlk.Grid.Cell[i][j].A)*grad_rho_abs/SolnBlk.W[i][j].rho;
             grad_rho_criteria_max = max(grad_rho_criteria_max, grad_rho_criteria);

 //             // Evaluate refinement criteria #2 based on the divergence
//              // of the velocity vector.
//              div_V = SolnBlk.dWdx[i][j].v.x + SolnBlk.dWdy[i][j].v.y;
//              div_V_criteria = sqrt(SolnBlk.Grid.Cell[i][j].A)*fabs(div_V)/SolnBlk.W[i][j].a();
//              div_V_criteria_max = max(div_V_criteria_max, div_V_criteria);

//              // Evaluate refinement criteria #3 based on the curl
//              // of the velocity vector.
//              curl_V_z = SolnBlk.dWdx[i][j].v.y - SolnBlk.dWdy[i][j].v.x; 
//              curl_V_abs = sqrt(sqr(curl_V_z)); 
//              curl_V_criteria = sqrt(SolnBlk.Grid.Cell[i][j].A)*curl_V_abs/SolnBlk.W[i][j].a(); 
//              curl_V_criteria_max = max(curl_V_criteria_max, curl_V_criteria);

	     // Evaluate refinement criteria #4 based on the gradient
	     // of the Temperature
	     grad_Temp_x = (ONE/( SolnBlk.W[i][j].rho* SolnBlk.W[i][j].Rtot())) * 
	       (SolnBlk.dWdx[i][j].p - (SolnBlk.W[i][j].p/SolnBlk.W[i][j].rho) * SolnBlk.dWdx[i][j].rho);
             grad_Temp_y = (ONE/( SolnBlk.W[i][j].rho* SolnBlk.W[i][j].Rtot())) *
	       (SolnBlk.dWdy[i][j].p - (SolnBlk.W[i][j].p/SolnBlk.W[i][j].rho) * SolnBlk.dWdy[i][j].rho);
             grad_Temp_abs = sqrt(sqr(grad_Temp_x) + sqr(grad_Temp_y));
             grad_Temp_criteria = sqrt(SolnBlk.Grid.Cell[i][j].A)*grad_Temp_abs/SolnBlk.W[i][j].T();
             grad_Temp_criteria_max = max(grad_Temp_criteria_max, grad_Temp_criteria);

// 	     // Evaluate refinement criteria #5 based on the gradient
// 	     // based on the gradient of CH4 mass fraction	     	     
// 	     grad_CH4_x = SolnBlk.dWdx[i][j].spec[0].c;
//              grad_CH4_y = SolnBlk.dWdy[i][j].spec[0].c;
//              grad_CH4_abs = sqrt(sqr(grad_CH4_x) + sqr(grad_CH4_y));
// 	     grad_CH4_criteria = sqrt(SolnBlk.Grid.Cell[i][j].A)*grad_CH4_abs;
//              grad_CH4_criteria_max = max(grad_CH4_criteria_max, grad_CH4_criteria);

	     // Evaluate refinement criteria #6 based on the gradient
	     // based on the gradient of CO2 mass fraction	     	     
	     grad_CO2_x = SolnBlk.dWdx[i][j].spec[2].c;
             grad_CO2_y = SolnBlk.dWdy[i][j].spec[2].c;
             grad_CO2_abs = sqrt(sqr(grad_CO2_x) + sqr(grad_CO2_y));
	     grad_CO2_criteria = sqrt(SolnBlk.Grid.Cell[i][j].A)*grad_CO2_abs;
             grad_CO2_criteria_max = max(grad_CO2_criteria_max, grad_CO2_criteria);

          } /* endif */
       } /* endfor */
    } /* endfor */

    /* Return the refinement criteria. */

    refinement_criteria[0] = grad_rho_criteria_max;
//     refinement_criteria[1] = div_V_criteria_max;
//     refinement_criteria[2] = curl_V_criteria_max;
    refinement_criteria[1] = grad_Temp_criteria_max;
//    refinement_criteria[2] = grad_CH4_criteria_max;
    refinement_criteria[2] = grad_CO2_criteria_max;


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
void Fix_Refined_Block_Boundaries(Chem2D_Quad_Block &SolnBlk,
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

    /* Reset the boundary condition types at the block boundaries. */
 
    Set_BCs(SolnBlk.Grid);

    /* Recompute the exterior nodes for the block quadrilateral mesh. */

    Update_Exterior_Nodes(SolnBlk.Grid);

    /* Recompute the cells for the block quadrilateral mesh. */

    Update_Cells(SolnBlk.Grid);

}

/********************************************************
 * Routine: Unfix_Refined_Block_Boundaries              *
 *                                                      *
 * Returns the adjusted the locations of the boundary   *
 * nodes of a solution block to their original          *
 * unmodified positions.                                *
 *                                                      *
 ********************************************************/
void Unfix_Refined_Block_Boundaries(Chem2D_Quad_Block &SolnBlk) {

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
	  SolnBlk.U[i][SolnBlk.JCu] = (SolnBlk.Grid.Cell[i][SolnBlk.JCu].A/
             SolnBlk.Grid.area(i, SolnBlk.JCu))*SolnBlk.U[i][SolnBlk.JCu];
	  SolnBlk.W[i][SolnBlk.JCu] = W(SolnBlk.U[i][SolnBlk.JCu]);
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
	  SolnBlk.U[i][SolnBlk.JCl] = (SolnBlk.Grid.Cell[i][SolnBlk.JCl].A/
             SolnBlk.Grid.area(i, SolnBlk.JCl))*SolnBlk.U[i][SolnBlk.JCl];
	  SolnBlk.W[i][SolnBlk.JCl] = W(SolnBlk.U[i][SolnBlk.JCl]);
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
	  SolnBlk.U[SolnBlk.ICu][j] = (SolnBlk.Grid.Cell[SolnBlk.ICu][j].A/
             SolnBlk.Grid.area(SolnBlk.ICu, j))*SolnBlk.U[SolnBlk.ICu][j];
	  SolnBlk.W[SolnBlk.ICu][j] = W(SolnBlk.U[SolnBlk.ICu][j]);
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
	  SolnBlk.U[SolnBlk.ICl][j] = (SolnBlk.Grid.Cell[SolnBlk.ICl][j].A/
             SolnBlk.Grid.area(SolnBlk.ICl, j))*SolnBlk.U[SolnBlk.ICl][j];
	  SolnBlk.W[SolnBlk.ICl][j] = W(SolnBlk.U[SolnBlk.ICl][j]);
       } /* endfor */
    } /* endif */

    /* Reset the boundary condition types at the block boundaries. */
 
    Set_BCs(SolnBlk.Grid);

    /* Recompute the exterior nodes for the block quadrilateral mesh. */

    Update_Exterior_Nodes(SolnBlk.Grid);

    /* Recompute the cells for the block quadrilateral mesh. */

    Update_Cells(SolnBlk.Grid);

}

/****************************************************************
 * Routine: Apply_Boundary_Flux_Corrections                     *
 *                                                              *
 * Apply flux corrections at boundaries of the solution         *
 * block to ensure that the scheme is conservative at           *
 * boundaries with mesh resolution changes.                     *
 *                                                              *
 ****************************************************************/
void Apply_Boundary_Flux_Corrections(Chem2D_Quad_Block &SolnBlk,
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
void Apply_Boundary_Flux_Corrections_Multistage_Explicit(Chem2D_Quad_Block &SolnBlk,
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
 *  solution block using a 2nd-ororder limited upwind   *
 * finite-volume spatial discretization scheme with     *
 * either the Godunov, Roe, Rusanov, HLLE, Linde, or    *
 * HLLC flux functions.                                 *
 * The residual is stored in dUdt[][][0].               *
 *                                                      *
 ********************************************************/
int dUdt_Residual_Evaluation(Chem2D_Quad_Block &SolnBlk,
			     Chem2D_Input_Parameters &Input_Parameters) {

  int i, j;
  Vector2D dX;
  Chem2D_pState Wl, Wr;
  Chem2D_cState Flux;

  Chem2D_pState W, W_face, dWdx, dWdy;  
  int NUM_VAR_CHEM2D = SolnBlk.NumVar();
  double delta_n; 

  /* Perform the linear reconstruction within each cell
     of the computational grid for this stage. */
  
  switch(Input_Parameters.i_Reconstruction) {
  case RECONSTRUCTION_GREEN_GAUSS :
    Linear_Reconstruction_GreenGauss(SolnBlk,
				     Input_Parameters.i_Limiter);    
    break;
  case RECONSTRUCTION_LEAST_SQUARES :
    Linear_Reconstruction_LeastSquares(SolnBlk,
				       Input_Parameters.i_Limiter);
   //  // a flag...
//     Linear_Reconstruction_LeastSquares_Diamond(SolnBlk,
// 				       Input_Parameters.i_Limiter);
    break;
  default:
    Linear_Reconstruction_LeastSquares(SolnBlk,
				       Input_Parameters.i_Limiter);
    break;
  }
  
  /********************************************************/
  /********************************************************/
  /* Compute viscous stresses and heat conduction vector if using cell-centered methods. */
  if (SolnBlk.Flow_Type != FLOWTYPE_INVISCID && 
      Input_Parameters.i_Viscous_Flux_Evaluation == VISCOUS_RECONSTRUCTION_ARITHMETIC_AVERAGE) {
     Viscous_Calculations(SolnBlk);
  }

  /********************************************************/
  /********************************************************/
  
    /* Evaluate the time rate of change of the solution
       (i.e., the solution residuals) using a second-order
       limited upwind scheme with a variety of flux functions. */
    
    // Add i-direction (zeta-direction) fluxes.
    for ( j  = SolnBlk.JCl-1 ; j <= SolnBlk.JCu+1 ; ++j ) {
      SolnBlk.dUdt[SolnBlk.ICl-1][j][0].Vacuum();
      
      for ( i = SolnBlk.ICl-1 ; i <= SolnBlk.ICu ; ++i ) {
	SolnBlk.dUdt[i+1][j][0].Vacuum();
	
	if ( j > SolnBlk.JCl-1 && j < SolnBlk.JCu+1 ) {
	  
	  /* Evaluate the cell interface i-direction fluxes. */
	  if (i == SolnBlk.ICl-1 && 
	      (SolnBlk.Grid.BCtypeW[j] == BC_REFLECTION ||
	       SolnBlk.Grid.BCtypeW[j] == BC_CHARACTERISTIC ||
	       SolnBlk.Grid.BCtypeW[j] == BC_FREE_SLIP ||
	       SolnBlk.Grid.BCtypeW[j] == BC_NO_SLIP ||
	       SolnBlk.Grid.BCtypeW[j] == BC_MOVING_WALL ||
	       SolnBlk.Grid.BCtypeW[j] == BC_FLAME_INFLOW ||
	       SolnBlk.Grid.BCtypeW[j] == BC_FLAME_OUTFLOW )) {
	    
	    dX = SolnBlk.Grid.xfaceW(i+1, j)-SolnBlk.Grid.Cell[i+1][j].Xc;	      	      
	    Wr = SolnBlk.W[i+1][j] + 
	      (SolnBlk.phi[i+1][j]^SolnBlk.dWdx[i+1][j])*dX.x +
	      (SolnBlk.phi[i+1][j]^SolnBlk.dWdy[i+1][j])*dX.y;	
	    	    
	    if (SolnBlk.Grid.BCtypeW[j] == BC_REFLECTION) {
	      Wl = Reflect(Wr, SolnBlk.Grid.nfaceW(i+1, j));
	    } else if (SolnBlk.Grid.BCtypeW[j] == BC_FREE_SLIP) {
	      Wl = Free_Slip(Wr,SolnBlk.WoW[j], SolnBlk.Grid.nfaceW(i+1, j),FIXED_TEMPERATURE_WALL);	      
	    } else if (SolnBlk.Grid.BCtypeW[j] == BC_NO_SLIP) {
	      Wl = No_Slip(Wr,SolnBlk.WoW[j], SolnBlk.Grid.nfaceW(i+1, j),FIXED_TEMPERATURE_WALL);
	    } else if (SolnBlk.Grid.BCtypeW[j] == BC_MOVING_WALL) {
	      Wl = Moving_Wall(Wr,SolnBlk.WoW[j], SolnBlk.Grid.nfaceW(i+1, j),
			       SolnBlk.Moving_wall_velocity,FIXED_TEMPERATURE_WALL);

	    } else if (SolnBlk.Grid.BCtypeW[j] == BC_FLAME_INFLOW){
	      Wl = BC_Flame_Inflow(Wr, 
				   SolnBlk.WoW[j],
				   SolnBlk.W[SolnBlk.ICu][j],
				   SolnBlk.Grid.nfaceW(i+1, j));
	      
	    } else if (SolnBlk.Grid.BCtypeW[j] == BC_FLAME_OUTFLOW){
	      Wl = BC_Flame_Outflow(Wr, 
				    SolnBlk.WoW[j], 
				    SolnBlk.W[SolnBlk.ICu][j],
				    SolnBlk.Grid.nfaceW(i+1, j));
	    } else { 
	      Wl = BC_Characteristic_Pressure(Wr, 
					      SolnBlk.WoW[j], 
					      SolnBlk.Grid.nfaceW(i+1, j));
	    } /* endif */

	  } else if (i == SolnBlk.ICu && 
		     (SolnBlk.Grid.BCtypeE[j] == BC_REFLECTION ||
		      SolnBlk.Grid.BCtypeE[j] == BC_CHARACTERISTIC ||
		      SolnBlk.Grid.BCtypeE[j] == BC_FREE_SLIP ||
		      SolnBlk.Grid.BCtypeE[j] == BC_NO_SLIP ||
		      SolnBlk.Grid.BCtypeE[j] == BC_MOVING_WALL || 
		      SolnBlk.Grid.BCtypeE[j] == BC_FLAME_INFLOW ||
		      SolnBlk.Grid.BCtypeE[j] == BC_FLAME_OUTFLOW )) {

	    dX = SolnBlk.Grid.xfaceE(i, j)-SolnBlk.Grid.Cell[i][j].Xc;
	    Wl = SolnBlk.W[i][j] + 
	      (SolnBlk.phi[i][j]^SolnBlk.dWdx[i][j])*dX.x +
	      (SolnBlk.phi[i][j]^SolnBlk.dWdy[i][j])*dX.y;	    
	    if (SolnBlk.Grid.BCtypeE[j] == BC_REFLECTION) {
	      Wr = Reflect(Wl, SolnBlk.Grid.nfaceE(i, j));   
	    } else if (SolnBlk.Grid.BCtypeE[j] == BC_FREE_SLIP) {
	      Wr = Free_Slip(Wl, SolnBlk.WoE[j],SolnBlk.Grid.nfaceE(i, j),FIXED_TEMPERATURE_WALL);
	    } else if (SolnBlk.Grid.BCtypeE[j] == BC_NO_SLIP) {
	      Wr = No_Slip(Wl, SolnBlk.WoE[j],SolnBlk.Grid.nfaceE(i, j),FIXED_TEMPERATURE_WALL);
	    } else if (SolnBlk.Grid.BCtypeE[j] == BC_MOVING_WALL) {
	      Wr = Moving_Wall(Wl,SolnBlk.WoE[j], SolnBlk.Grid.nfaceE(i, j),
			       SolnBlk.Moving_wall_velocity,FIXED_TEMPERATURE_WALL);
	    } else if (SolnBlk.Grid.BCtypeE[j] == BC_FLAME_INFLOW){
	      Wr = BC_Flame_Inflow(Wl, 
				   SolnBlk.WoE[j],
				   SolnBlk.W[SolnBlk.ICl][j],
				   SolnBlk.Grid.nfaceE(i, j));
	    } else if (SolnBlk.Grid.BCtypeE[j] == BC_FLAME_OUTFLOW){
	      Wr = BC_Flame_Outflow(Wl, 
				    SolnBlk.WoE[j],
				    SolnBlk.W[SolnBlk.ICl][j], 
				    SolnBlk.Grid.nfaceE(i, j));
	    } else {
	      Wr = BC_Characteristic_Pressure(Wl, 
					      SolnBlk.WoE[j], 
					      SolnBlk.Grid.nfaceE(i, j));
	    } /* endif */
	  } else {            
	    dX = SolnBlk.Grid.xfaceE(i, j)-SolnBlk.Grid.Cell[i][j].Xc;
	    Wl = SolnBlk.W[i][j] + 
	      (SolnBlk.phi[i][j]^SolnBlk.dWdx[i][j])*dX.x +
	      (SolnBlk.phi[i][j]^SolnBlk.dWdy[i][j])*dX.y;
	    dX = SolnBlk.Grid.xfaceW(i+1, j)-SolnBlk.Grid.Cell[i+1][j].Xc;
	    Wr = SolnBlk.W[i+1][j] + 
	      (SolnBlk.phi[i+1][j]^SolnBlk.dWdx[i+1][j])*dX.x +
	      (SolnBlk.phi[i+1][j]^SolnBlk.dWdy[i+1][j])*dX.y;	       
	  } /* endif */
	  
	  // Spacing for Preconditioner 
	  if(SolnBlk.Flow_Type != FLOWTYPE_INVISCID && 
             Input_Parameters.Preconditioning){
	    delta_n = min( TWO*(SolnBlk.Grid.Cell[i][j].A/
				(SolnBlk.Grid.lfaceE(i, j)+SolnBlk.Grid.lfaceW(i, j))),
			   TWO*(SolnBlk.Grid.Cell[i][j].A/
				(SolnBlk.Grid.lfaceN(i, j)+SolnBlk.Grid.lfaceS(i, j))));
	  } 

	  switch(Input_Parameters.i_Flux_Function) {
	    //                case FLUX_FUNCTION_GODUNOV :
	    //                  Flux = FluxGodunov_n(Wl, Wr, SolnBlk.Grid.nfaceE(i, j));
	    //                  break;
	  case FLUX_FUNCTION_ROE :
	    Flux = FluxRoe_n(Wl, Wr, SolnBlk.Grid.nfaceE(i, j),
			     Input_Parameters.Preconditioning,delta_n);
	    break;
	    //                case FLUX_FUNCTION_RUSANOV :
	    //                  Flux = FluxRusanov_n(Wl, Wr, SolnBlk.Grid.nfaceE(i, j));
	    //                  break;
	  case FLUX_FUNCTION_HLLE :  
	    Flux = FluxHLLE_n(Wl, Wr, SolnBlk.Grid.nfaceE(i, j),Input_Parameters.Preconditioning);
	    break;
	  case FLUX_FUNCTION_LINDE :
	    Flux = FluxLinde_n(Wl, Wr, SolnBlk.Grid.nfaceE(i, j));
	    break;
	    //                case FLUX_FUNCTION_HLLC :
	    //                  Flux = FluxHLLC_n(Wl, Wr, SolnBlk.Grid.nfaceE(i, j));
	    //                  break;
	    //                default:
	    //                  Flux = FluxRoe_n(Wl, Wr, SolnBlk.Grid.nfaceE(i, j));
	    //                  break;
	  } /* endswitch */
    
	  /*****************************************/
	  /***************** VISCOUS ***************/
	  //zeta
	  if(SolnBlk.Flow_Type != FLOWTYPE_INVISCID){ // -ve as on RHS 	  
	    switch(Input_Parameters.i_Viscous_Flux_Evaluation){

	    case VISCOUS_RECONSTRUCTION_ARITHMETIC_AVERAGE :
	      //Arithmetic Mean
	      Flux -= Viscous_FluxArithmetic_n(SolnBlk.U[i][j],
                                               SolnBlk.dWdx[i][j],
                                               SolnBlk.dWdy[i][j],
                                               SolnBlk.U[i+1][j],
                                               SolnBlk.dWdx[i+1][j],
                                               SolnBlk.dWdy[i+1][j],
					       SolnBlk.Grid.nfaceE(i,j));
	      break;
	    case VISCOUS_RECONSTRUCTION_CARTESIAN : 
	      //Reduced Stencil 1D
	      //W = HALF*(SolnBlk.W[i][j] + SolnBlk.W[i+1][j]);

	      W = SolnBlk.W[i+1][j] - (SolnBlk.W[i+1][j] - SolnBlk.W[i][j]) * 
                  (SolnBlk.Grid.Cell[i+1][j].Xc.x-SolnBlk.Grid.xfaceE(i,j).x)
                  /(SolnBlk.Grid.Cell[i+1][j].Xc.x - SolnBlk.Grid.Cell[i][j].Xc.x);
            
	      dWdx = (SolnBlk.W[i+1][j] - SolnBlk.W[i][j])/
		fabs(SolnBlk.Grid.Cell[i][j].Xc.x - SolnBlk.Grid.Cell[i+1][j].Xc.x);
	      dWdy = HALF*(SolnBlk.dWdy[i][j] + SolnBlk.dWdy[i+1][j]);
	      //dWdy.Vacuum(); 
	      // Compute the viscous flux.
	      Flux -= Viscous_Flux_n(W,dWdx,dWdy,
                                     SolnBlk.Axisymmetric,
				     SolnBlk.Grid.xfaceE(i,j),
				     SolnBlk.Grid.nfaceE(i,j));
	      break;	      
	    case VISCOUS_RECONSTRUCTION_DIAMOND_PATH :
	      W_face = HALF*(SolnBlk.Wn(i+1,j+1) +SolnBlk.Wn(i+1,j));
	      dWdx = SolnBlk.dWdx_faceE[i][j];
 	      dWdy = SolnBlk.dWdy_faceE[i][j];
	      Flux -= Viscous_Flux_n(W_face,dWdx,dWdy,
				     SolnBlk.Axisymmetric,
				     SolnBlk.Grid.xfaceE(i,j),
				     SolnBlk.Grid.nfaceE(i,j));
	    
	      break;
	    }	
	  }
	  /*****************************************/
	  /*****************************************/

	  /* Evaluate cell-averaged solution changes. */
          SolnBlk.dUdt[i][j][0] -= 
              Flux*SolnBlk.Grid.lfaceE(i, j)/SolnBlk.Grid.Cell[i][j].A;
           SolnBlk.dUdt[i+1][j][0] += 
              Flux*SolnBlk.Grid.lfaceW(i+1, j)/SolnBlk.Grid.Cell[i+1][j].A;

	  /* Include axisymmetric source terms as required. */
	  if (SolnBlk.Axisymmetric) {
	    SolnBlk.dUdt[i][j][0] += SolnBlk.W[i][j].Sa_inviscid(SolnBlk.Grid.Cell[i][j].Xc,
                                                                 SolnBlk.Axisymmetric);

	    // Include viscous if specified
	    if (SolnBlk.Flow_Type != FLOWTYPE_INVISCID) { //-ve as on RHS 
	      SolnBlk.dUdt[i][j][0] += SolnBlk.W[i][j].Sa_viscous(SolnBlk.dWdx[i][j],
	       						          SolnBlk.dWdy[i][j],
								  SolnBlk.Grid.Cell[i][j].Xc,
								  SolnBlk.Axisymmetric);
	    }
          } /* endif */

          /* Include source terms associated with turbulence model */
          if (SolnBlk.Flow_Type == FLOWTYPE_TURBULENT_RANS_K_OMEGA ||
              SolnBlk.Flow_Type == FLOWTYPE_TURBULENT_RANS_K_EPSILON) {
	      SolnBlk.dUdt[i][j][0] += SolnBlk.W[i][j].S_turbulence_model(SolnBlk.dWdx[i][j],
	       						                  SolnBlk.dWdy[i][j],
								          SolnBlk.Grid.Cell[i][j].Xc,
                                                                          SolnBlk.Axisymmetric);
          } /* endif */

	  /* Include source terms associated with the finite-rate chemistry and 
             turbulence/chemistry interactions */ 
	  if (SolnBlk.W[i][j].React.reactset_flag != NO_REACTIONS) {	  
	    SolnBlk.dUdt[i][j][0] += SolnBlk.W[i][j].Sw(SolnBlk.W[i][j].React.reactset_flag,SolnBlk.W[i][j].flow_type);
	  } /* endif */

          /* Include source terms associated with gravity */
	  if (SolnBlk.Gravity) {	 
	      SolnBlk.dUdt[i][j][0] += SolnBlk.W[i][j].Sg();
          }   
	  /*****************************************/
	  /*****************************************/

	  /* Save west and east face boundary flux. */	
	  if (i == SolnBlk.ICl-1) {
	    SolnBlk.FluxW[j] = -Flux*SolnBlk.Grid.lfaceW(i+1, j);
	  } else if (i == SolnBlk.ICu) {
	    SolnBlk.FluxE[j] = Flux*SolnBlk.Grid.lfaceE(i, j);  
	  } 
	         
	  } /* endif */
       } /* endfor */
      
       if ( j > SolnBlk.JCl-1 && j < SolnBlk.JCu+1 ){
	 SolnBlk.dUdt[SolnBlk.ICu+1][j][0].Vacuum();
	 SolnBlk.dUdt[SolnBlk.ICl-1][j][0].Vacuum();
       }
    } /* endfor */

    
    // Add j-direction (eta-direction) fluxes.
    for ( i = SolnBlk.ICl ; i <= SolnBlk.ICu ; ++i ) {
      for ( j  = SolnBlk.JCl-1 ; j <= SolnBlk.JCu ; ++j ) {
    
	/* Evaluate the cell interface j-direction fluxes. */ 
	if (j == SolnBlk.JCl-1 && 
	    (SolnBlk.Grid.BCtypeS[i] == BC_REFLECTION ||
	     SolnBlk.Grid.BCtypeS[i] == BC_CHARACTERISTIC ||
	     SolnBlk.Grid.BCtypeS[i] == BC_FREE_SLIP ||
	     SolnBlk.Grid.BCtypeS[i] == BC_NO_SLIP ||
	     SolnBlk.Grid.BCtypeS[i] == BC_MOVING_WALL ||
	     SolnBlk.Grid.BCtypeS[i] == BC_FLAME_OUTFLOW )) {

	  dX = SolnBlk.Grid.xfaceS(i, j+1)-SolnBlk.Grid.Cell[i][j+1].Xc;
	  Wr = SolnBlk.W[i][j+1] +
	    (SolnBlk.phi[i][j+1]^SolnBlk.dWdx[i][j+1])*dX.x +
	    (SolnBlk.phi[i][j+1]^SolnBlk.dWdy[i][j+1])*dX.y;
	  if (SolnBlk.Grid.BCtypeS[i] == BC_REFLECTION) {
	    Wl = Reflect(Wr, SolnBlk.Grid.nfaceS(i, j+1)); 
	  } else if (SolnBlk.Grid.BCtypeS[i] == BC_FREE_SLIP) {
	    Wl = Free_Slip(Wr,SolnBlk.WoS[i], SolnBlk.Grid.nfaceS(i, j+1),FIXED_TEMPERATURE_WALL);
	  } else if (SolnBlk.Grid.BCtypeS[i] == BC_NO_SLIP) {
	    Wl = No_Slip(Wr,SolnBlk.WoS[i], SolnBlk.Grid.nfaceS(i, j+1),FIXED_TEMPERATURE_WALL);
	  } else if (SolnBlk.Grid.BCtypeS[i] == BC_MOVING_WALL) {
	    Wl = Moving_Wall(Wr,SolnBlk.WoS[i], SolnBlk.Grid.nfaceS(i, j+1),
			     SolnBlk.Moving_wall_velocity,FIXED_TEMPERATURE_WALL);
	  } else if (SolnBlk.Grid.BCtypeS[i] == BC_FLAME_OUTFLOW){
	    Wl = BC_Flame_Outflow(Wr, 
				  SolnBlk.WoS[i], 
				  SolnBlk.W[i][SolnBlk.JCu],
				  SolnBlk.Grid.nfaceS(i, j+1));
	  } else { 
	    Wl = BC_Characteristic_Pressure(Wr, 
					    SolnBlk.WoS[i], 
					    SolnBlk.Grid.nfaceS(i, j+1));
	  } /* endif */

	} else if (j == SolnBlk.JCu && 
		   (SolnBlk.Grid.BCtypeN[i] == BC_REFLECTION ||
		    SolnBlk.Grid.BCtypeN[i] == BC_CHARACTERISTIC ||
		    SolnBlk.Grid.BCtypeN[i] == BC_FREE_SLIP ||
		    SolnBlk.Grid.BCtypeN[i] == BC_NO_SLIP ||
		    SolnBlk.Grid.BCtypeN[i] == BC_MOVING_WALL ||
		    SolnBlk.Grid.BCtypeN[i] == BC_FLAME_OUTFLOW )) {
	  
	  dX = SolnBlk.Grid.xfaceN(i, j)-SolnBlk.Grid.Cell[i][j].Xc;
	  Wl = SolnBlk.W[i][j] + 
	    (SolnBlk.phi[i][j]^SolnBlk.dWdx[i][j])*dX.x +
	    (SolnBlk.phi[i][j]^SolnBlk.dWdy[i][j])*dX.y;
	  if (SolnBlk.Grid.BCtypeN[i] == BC_REFLECTION) {
	    Wr = Reflect(Wl, SolnBlk.Grid.nfaceN(i, j));	
	  } else if (SolnBlk.Grid.BCtypeN[i] == BC_FREE_SLIP) {
	    Wr = Free_Slip(Wl, SolnBlk.WoN[i], SolnBlk.Grid.nfaceN(i, j),FIXED_TEMPERATURE_WALL);
	  } else if (SolnBlk.Grid.BCtypeN[i] == BC_NO_SLIP) {
	    Wr = No_Slip(Wl, SolnBlk.WoN[i], SolnBlk.Grid.nfaceN(i, j),FIXED_TEMPERATURE_WALL);
	  } else if (SolnBlk.Grid.BCtypeN[i] == BC_MOVING_WALL) {
	    Wr = Moving_Wall(Wl,SolnBlk.WoN[i],  SolnBlk.Grid.nfaceN(i, j),
			     SolnBlk.Moving_wall_velocity,FIXED_TEMPERATURE_WALL);
	  } else if (SolnBlk.Grid.BCtypeN[i] == BC_FLAME_OUTFLOW){
	    Wr = BC_Flame_Outflow(Wl, 
				  SolnBlk.WoN[i], 
				  SolnBlk.W[i][SolnBlk.JCl],
				  SolnBlk.Grid.nfaceW(i, j));
	  } else {
	    Wr = BC_Characteristic_Pressure(Wl, 
					    SolnBlk.WoN[i], 
					    SolnBlk.Grid.nfaceN(i, j));
	  } /* endif */
	  
	} else {
	  dX = SolnBlk.Grid.xfaceN(i, j)-SolnBlk.Grid.Cell[i][j].Xc;
	  Wl = SolnBlk.W[i][j] + 
	    (SolnBlk.phi[i][j]^SolnBlk.dWdx[i][j])*dX.x +
	    (SolnBlk.phi[i][j]^SolnBlk.dWdy[i][j])*dX.y;
	  dX = SolnBlk.Grid.xfaceS(i, j+1)-SolnBlk.Grid.Cell[i][j+1].Xc;
	  Wr = SolnBlk.W[i][j+1] +
	    (SolnBlk.phi[i][j+1]^SolnBlk.dWdx[i][j+1])*dX.x +
	    (SolnBlk.phi[i][j+1]^SolnBlk.dWdy[i][j+1])*dX.y;
	} /* endif */
	  
	// Spacing for Preconditioner 
	if(SolnBlk.Flow_Type != FLOWTYPE_INVISCID && 
           Input_Parameters.Preconditioning){
	  delta_n = min( TWO*(SolnBlk.Grid.Cell[i][j].A/
			      (SolnBlk.Grid.lfaceE(i, j)+SolnBlk.Grid.lfaceW(i, j))),
			 TWO*(SolnBlk.Grid.Cell[i][j].A/
			      (SolnBlk.Grid.lfaceN(i, j)+SolnBlk.Grid.lfaceS(i, j))));
	} 

 	switch(Input_Parameters.i_Flux_Function) {
	  //             case FLUX_FUNCTION_GODUNOV :
	  //               Flux = FluxGodunov_n(Wl, Wr, SolnBlk.Grid.nfaceN(i, j));
	  //               break;
	case FLUX_FUNCTION_ROE :
	  Flux = FluxRoe_n(Wl, Wr, SolnBlk.Grid.nfaceN(i, j),
			   Input_Parameters.Preconditioning,delta_n);
	break;
	  //             case FLUX_FUNCTION_RUSANOV :
	  //               Flux = FluxRusanov_n(Wl, Wr, SolnBlk.Grid.nfaceN(i, j));
	  //               break;
	case FLUX_FUNCTION_HLLE :
	  Flux = FluxHLLE_n(Wl, Wr, SolnBlk.Grid.nfaceN(i, j),
			    Input_Parameters.Preconditioning);
	  break;
	case FLUX_FUNCTION_LINDE :
	  Flux = FluxLinde_n(Wl, Wr, SolnBlk.Grid.nfaceN(i, j));
	  break;
	  //             case FLUX_FUNCTION_HLLC :
	  //               Flux = FluxHLLC_n(Wl, Wr, SolnBlk.Grid.nfaceN(i, j));
	  //               break;
	  //             default:
	  //               Flux = FluxRoe_n(Wl, Wr, SolnBlk.Grid.nfaceN(i, j));
	  //               break;
	} /* endswitch */
	
	/*****************************************/
	/***************** VISCOUS ***************/
	//eta  
	if(SolnBlk.Flow_Type != FLOWTYPE_INVISCID){ // -ve as on RHS    
	  switch(Input_Parameters.i_Viscous_Flux_Evaluation){	    
	  case VISCOUS_RECONSTRUCTION_ARITHMETIC_AVERAGE :
	    //Arithmetic Mean
	    Flux -= Viscous_FluxArithmetic_n(SolnBlk.U[i][j],
                                             SolnBlk.dWdx[i][j],
                                             SolnBlk.dWdy[i][j],
                                             SolnBlk.U[i][j+1],
                                             SolnBlk.dWdx[i][j+1],
                                             SolnBlk.dWdy[i][j+1],
					     SolnBlk.Grid.nfaceN(i,j));
	    break;
	  case VISCOUS_RECONSTRUCTION_CARTESIAN :
	    //Reduced Stencil 1D
	    // W = HALF*(SolnBlk.W[i][j] + SolnBlk.W[i][j+1]);   
	    W = SolnBlk.W[i][j+1] - (SolnBlk.W[i][j+1] - SolnBlk.W[i][j]) * 
	      (SolnBlk.Grid.Cell[i][j+1].Xc.y - SolnBlk.Grid.xfaceN(i,j).y)
	      /(SolnBlk.Grid.Cell[i][j+1].Xc.y - SolnBlk.Grid.Cell[i][j].Xc.y);	    

	    dWdx = HALF*(SolnBlk.dWdx[i][j] + SolnBlk.dWdx[i][j+1]);
	    dWdy = (SolnBlk.W[i][j+1] - SolnBlk.W[i][j])/
	      fabs(SolnBlk.Grid.Cell[i][j+1].Xc.y - SolnBlk.Grid.Cell[i][j].Xc.y);
	  
	    // Compute the viscous flux.
	    Flux -= Viscous_Flux_n(W,dWdx,dWdy,
                                   SolnBlk.Axisymmetric,
				   SolnBlk.Grid.xfaceN(i,j),
				   SolnBlk.Grid.nfaceN(i,j));
	 
	    break;	  
	  case VISCOUS_RECONSTRUCTION_DIAMOND_PATH :
	    W_face = HALF*(SolnBlk.Wn(i,j+1) + SolnBlk.Wn(i+1,j+1));
	    dWdx = SolnBlk.dWdx_faceN[i][j];
	    dWdy = SolnBlk.dWdy_faceN[i][j];
	    Flux -= Viscous_Flux_n(W_face,dWdx,dWdy,
                                   SolnBlk.Axisymmetric,
				   SolnBlk.Grid.xfaceN(i,j),
				   SolnBlk.Grid.nfaceN(i,j));
	    break;
	  }
	}
	/*****************************************/
	/*****************************************/
	
	/* Evaluate cell-averaged solution changes. */
	
	SolnBlk.dUdt[i][j][0] -= 
	  Flux*SolnBlk.Grid.lfaceN(i, j)/
	  SolnBlk.Grid.Cell[i][j].A;
	SolnBlk.dUdt[i][j+1][0] += 
	  Flux*SolnBlk.Grid.lfaceS(i, j+1)/
	  SolnBlk.Grid.Cell[i][j+1].A;
	
        /* Save south and north face boundary flux. */
	// USED for AMR 
	if (j == SolnBlk.JCl-1) {
	  SolnBlk.FluxS[i] = -Flux*SolnBlk.Grid.lfaceS(i, j+1);
	} else if (j == SolnBlk.JCu) {
	  SolnBlk.FluxN[i] = Flux*SolnBlk.Grid.lfaceN(i, j);
	} /* endif */
	

      } /* endfor */

      SolnBlk.dUdt[i][SolnBlk.JCu+1][0].Vacuum();
      SolnBlk.dUdt[i][SolnBlk.JCl-1][0].Vacuum();
    } /* endfor */
    
    // For k-omega turbulence model set residual in laminar sublayer to zero.
    if (SolnBlk.Flow_Type == FLOWTYPE_TURBULENT_RANS_K_OMEGA) {
       for ( i = SolnBlk.ICl ; i <= SolnBlk.ICu ; ++i ) {
          for ( j  = SolnBlk.JCl ; j <= SolnBlk.JCu ; ++j ) {
	     if (SolnBlk.Wall[i][j].yplus <= Input_Parameters.yplus_sublayer) {
	        SolnBlk.dUdt[i][j][0].rhoomega = ZERO;
             } /* endif */
          } /* endfor */
       } /* endfor */
    } /* endif */

    /* Residual successfully calculated. */
    return (0);
    
}

/********************************************************
 * Routine: dUdt_Multistage_Explicit                    *
 *                                                      *
 * This routine determines the solution residuals for a *
 * given stage of a variety of multi-stage explicit     *
 * time integration schemes for a given solution block. *
 *                                                      *
 ********************************************************/ 
int dUdt_Multistage_Explicit(Chem2D_Quad_Block &SolnBlk,
                             const int i_stage,
                             Chem2D_Input_Parameters &Input_Parameters) {

    int i, j, k_residual;
    double omega;
    Vector2D dX;
    Chem2D_pState Wl, Wr;
    Chem2D_cState Flux;

    Chem2D_pState W, W_face, dWdx, dWdy;  
    int NUM_VAR_CHEM2D = SolnBlk.NumVar();
    double delta_n;

    Chem2D_pState Chem2D_W_VACUUM;
    Chem2D_cState Chem2D_U_VACUUM;
    Chem2D_W_VACUUM.Vacuum(); Chem2D_U_VACUUM.Vacuum();

    /* Evaluate the solution residual for stage i_stage of n_stage scheme. */

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
      break;
    case RECONSTRUCTION_LEAST_SQUARES :
      Linear_Reconstruction_LeastSquares(SolnBlk,
					 Input_Parameters.i_Limiter);
      // a flag...
    //   Linear_Reconstruction_LeastSquares_Diamond(SolnBlk,
// 					 Input_Parameters.i_Limiter);
      break;
    default:
      Linear_Reconstruction_LeastSquares(SolnBlk,
                                         Input_Parameters.i_Limiter);
      break;
    } /* endswitch */

   
    /********************************************************/
    /********************************************************/
  
    /* Compute viscous stresses and heat conduction vector if using cell-centered methods. */
    if (SolnBlk.Flow_Type != FLOWTYPE_INVISCID && 
        Input_Parameters.i_Viscous_Flux_Evaluation == VISCOUS_RECONSTRUCTION_ARITHMETIC_AVERAGE) {
       Viscous_Calculations(SolnBlk);
    }
  
    /********************************************************/
    /********************************************************/

    /* Evaluate the time rate of change of the solution
       (i.e., the solution residuals) using a second-order
       limited upwind scheme with a variety of flux functions. */
  
  
    // Add i-direction (zeta-direction) fluxes.
    for ( j  = SolnBlk.JCl-1 ; j <= SolnBlk.JCu+1 ; ++j ) {
      if ( i_stage == 1 ) {
	SolnBlk.Uo[SolnBlk.ICl-1][j] = SolnBlk.U[SolnBlk.ICl-1][j];
	SolnBlk.dUdt[SolnBlk.ICl-1][j][k_residual].Vacuum();
      } else {
	SolnBlk.dUdt[SolnBlk.ICl-1][j][k_residual].Vacuum();
      } /* endif */
      
      for ( i = SolnBlk.ICl-1 ; i <= SolnBlk.ICu ; ++i ) {
	if ( i_stage == 1 ) {
	  SolnBlk.Uo[i+1][j] = SolnBlk.U[i+1][j];
	  SolnBlk.dUdt[i+1][j][k_residual].Vacuum();
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
	      SolnBlk.dUdt[i+1][j][k_residual].Vacuum();
	    } /* endif */
	    break;
	  case TIME_STEPPING_MULTISTAGE_OPTIMAL_SMOOTHING :
	    SolnBlk.dUdt[i+1][j][k_residual].Vacuum(); 
	    break;
	  default:
	    SolnBlk.dUdt[i+1][j][k_residual].Vacuum(); 
	    break;
	  } /* endswitch */
	} /* endif */

  
	if ( j > SolnBlk.JCl-1 && j < SolnBlk.JCu+1 ) {
	  
	  /* Evaluate the cell interface i-direction fluxes. */
	  if (i == SolnBlk.ICl-1 && 
	      (SolnBlk.Grid.BCtypeW[j] == BC_REFLECTION ||
	       SolnBlk.Grid.BCtypeW[j] == BC_CHARACTERISTIC ||
	       SolnBlk.Grid.BCtypeW[j] == BC_FREE_SLIP ||
	       SolnBlk.Grid.BCtypeW[j] == BC_NO_SLIP ||
	       SolnBlk.Grid.BCtypeW[j] == BC_MOVING_WALL ||
	       SolnBlk.Grid.BCtypeW[j] == BC_FLAME_INFLOW ||
	       SolnBlk.Grid.BCtypeW[j] == BC_FLAME_OUTFLOW )) {
	    
	    dX = SolnBlk.Grid.xfaceW(i+1, j)-SolnBlk.Grid.Cell[i+1][j].Xc;	      	      
	    Wr = SolnBlk.W[i+1][j] + 
	      (SolnBlk.phi[i+1][j]^SolnBlk.dWdx[i+1][j])*dX.x +
	      (SolnBlk.phi[i+1][j]^SolnBlk.dWdy[i+1][j])*dX.y;	     
	    
	    if (SolnBlk.Grid.BCtypeW[j] == BC_REFLECTION) {
	      Wl = Reflect(Wr, SolnBlk.Grid.nfaceW(i+1, j));
	    } else if (SolnBlk.Grid.BCtypeW[j] == BC_FREE_SLIP) {
	      Wl = Free_Slip(Wr,SolnBlk.WoW[j], SolnBlk.Grid.nfaceW(i+1, j),FIXED_TEMPERATURE_WALL);	      
	    } else if (SolnBlk.Grid.BCtypeW[j] == BC_NO_SLIP) {
	      Wl = No_Slip(Wr,SolnBlk.WoW[j], SolnBlk.Grid.nfaceW(i+1, j),FIXED_TEMPERATURE_WALL);
	    } else if (SolnBlk.Grid.BCtypeW[j] == BC_MOVING_WALL) {
	      Wl = Moving_Wall(Wr,SolnBlk.WoW[j], SolnBlk.Grid.nfaceW(i+1, j),
			       SolnBlk.Moving_wall_velocity,FIXED_TEMPERATURE_WALL);

	    } else if (SolnBlk.Grid.BCtypeW[j] == BC_FLAME_INFLOW){
	      Wl = BC_Flame_Inflow(Wr, 
				   SolnBlk.WoW[j],
				   SolnBlk.W[SolnBlk.ICu][j],
				   SolnBlk.Grid.nfaceW(i+1, j));
	      
	    } else if (SolnBlk.Grid.BCtypeW[j] == BC_FLAME_OUTFLOW){
	      Wl = BC_Flame_Outflow(Wr, 
				    SolnBlk.WoW[j], 
				    SolnBlk.W[SolnBlk.ICu][j],
				    SolnBlk.Grid.nfaceW(i+1, j));
	    } else { 
	      Wl = BC_Characteristic_Pressure(Wr, 
					      SolnBlk.WoW[j], 
					      SolnBlk.Grid.nfaceW(i+1, j));
	    } /* endif */
	    
	  } else if (i == SolnBlk.ICu && 
		     (SolnBlk.Grid.BCtypeE[j] == BC_REFLECTION ||
		      SolnBlk.Grid.BCtypeE[j] == BC_CHARACTERISTIC ||
		      SolnBlk.Grid.BCtypeE[j] == BC_FREE_SLIP ||
		      SolnBlk.Grid.BCtypeE[j] == BC_NO_SLIP ||
		      SolnBlk.Grid.BCtypeE[j] == BC_MOVING_WALL || 
		      SolnBlk.Grid.BCtypeE[j] == BC_FLAME_INFLOW ||
		      SolnBlk.Grid.BCtypeE[j] == BC_FLAME_OUTFLOW )) {

	    dX = SolnBlk.Grid.xfaceE(i, j)-SolnBlk.Grid.Cell[i][j].Xc;
	    Wl = SolnBlk.W[i][j] + 
	      (SolnBlk.phi[i][j]^SolnBlk.dWdx[i][j])*dX.x +
	      (SolnBlk.phi[i][j]^SolnBlk.dWdy[i][j])*dX.y;	    
	    if (SolnBlk.Grid.BCtypeE[j] == BC_REFLECTION) {
	      Wr = Reflect(Wl, SolnBlk.Grid.nfaceE(i, j));   
	    } else if (SolnBlk.Grid.BCtypeE[j] == BC_FREE_SLIP) {
	      Wr = Free_Slip(Wl, SolnBlk.WoE[j],SolnBlk.Grid.nfaceE(i, j),FIXED_TEMPERATURE_WALL);
	    } else if (SolnBlk.Grid.BCtypeE[j] == BC_NO_SLIP) {
	      Wr = No_Slip(Wl, SolnBlk.WoE[j],SolnBlk.Grid.nfaceE(i, j),FIXED_TEMPERATURE_WALL);
	    } else if (SolnBlk.Grid.BCtypeE[j] == BC_MOVING_WALL) {
	      Wr = Moving_Wall(Wl,SolnBlk.WoE[j], SolnBlk.Grid.nfaceE(i, j),
			       SolnBlk.Moving_wall_velocity,FIXED_TEMPERATURE_WALL);
	    } else if (SolnBlk.Grid.BCtypeE[j] == BC_FLAME_INFLOW){
	      Wr = BC_Flame_Inflow(Wl, 
				   SolnBlk.WoE[j],
				   SolnBlk.W[SolnBlk.ICl][j],
				   SolnBlk.Grid.nfaceE(i, j));
	    } else if (SolnBlk.Grid.BCtypeE[j] == BC_FLAME_OUTFLOW){
	      Wr = BC_Flame_Outflow(Wl, 
				    SolnBlk.WoE[j],
				    SolnBlk.W[SolnBlk.ICl][j],
				    SolnBlk.Grid.nfaceE(i, j));
	    } else {
	      Wr = BC_Characteristic_Pressure(Wl, SolnBlk.WoE[j], 
					      SolnBlk.Grid.nfaceE(i, j));
	    } /* endif */
	  } else {            
	    dX = SolnBlk.Grid.xfaceE(i, j)-SolnBlk.Grid.Cell[i][j].Xc;
	    Wl = SolnBlk.W[i][j] + 
	      (SolnBlk.phi[i][j]^SolnBlk.dWdx[i][j])*dX.x +
	      (SolnBlk.phi[i][j]^SolnBlk.dWdy[i][j])*dX.y;
	    dX = SolnBlk.Grid.xfaceW(i+1, j)-SolnBlk.Grid.Cell[i+1][j].Xc;
	    Wr = SolnBlk.W[i+1][j] + 
	      (SolnBlk.phi[i+1][j]^SolnBlk.dWdx[i+1][j])*dX.x +
	      (SolnBlk.phi[i+1][j]^SolnBlk.dWdy[i+1][j])*dX.y;	       
	  } /* endif */

	  // Spacing for Preconditioner 
	  if(SolnBlk.Flow_Type != FLOWTYPE_INVISCID && 
             Input_Parameters.Preconditioning){
	    delta_n = min( TWO*(SolnBlk.Grid.Cell[i][j].A/
				(SolnBlk.Grid.lfaceE(i, j)+SolnBlk.Grid.lfaceW(i, j))),
			   TWO*(SolnBlk.Grid.Cell[i][j].A/
				(SolnBlk.Grid.lfaceN(i, j)+SolnBlk.Grid.lfaceS(i, j))));
	  }
	
  
	  switch(Input_Parameters.i_Flux_Function) {
	    //                case FLUX_FUNCTION_GODUNOV :
	    //                  Flux = FluxGodunov_n(Wl, Wr, SolnBlk.Grid.nfaceE(i, j));
	    //                  break;
	   
	  case FLUX_FUNCTION_ROE :
	    Flux = FluxRoe_n(Wl, Wr, SolnBlk.Grid.nfaceE(i, j),
			     Input_Parameters.Preconditioning,delta_n);
	    break;
	    //                case FLUX_FUNCTION_RUSANOV :
	    //                  Flux = FluxRusanov_n(Wl, Wr, SolnBlk.Grid.nfaceE(i, j));
	    //                  break;
	  case FLUX_FUNCTION_HLLE :
	    Flux = FluxHLLE_n(Wl, Wr, SolnBlk.Grid.nfaceE(i, j),Input_Parameters.Preconditioning);
	    break;
	  case FLUX_FUNCTION_LINDE :
	    Flux = FluxLinde_n(Wl, Wr, SolnBlk.Grid.nfaceE(i, j));
	    break;
	    //                case FLUX_FUNCTION_HLLC :
	    //                  Flux = FluxHLLC_n(Wl, Wr, SolnBlk.Grid.nfaceE(i, j));
	    //                  break;
	    //                default:
	    //                  Flux = FluxRoe_n(Wl, Wr, SolnBlk.Grid.nfaceE(i, j));
	    //                  break;
	  } /* endswitch */
	  
	
	  /*****************************************/
	  /***************** VISCOUS ***************/
	  //zeta
	  if(SolnBlk.Flow_Type != FLOWTYPE_INVISCID){ // -ve as on RHS 	  
	    switch(Input_Parameters.i_Viscous_Flux_Evaluation){
	    case VISCOUS_RECONSTRUCTION_ARITHMETIC_AVERAGE :
	      //Arithmetic Mean
	      Flux -= Viscous_FluxArithmetic_n(SolnBlk.U[i][j],
                                               SolnBlk.dWdx[i][j],
                                               SolnBlk.dWdy[i][j],
                                               SolnBlk.U[i+1][j],
                                               SolnBlk.dWdx[i+1][j],
                                               SolnBlk.dWdy[i+1][j],
	        		               SolnBlk.Grid.nfaceE(i,j));  // d/dr Frv
	      break;
	    case VISCOUS_RECONSTRUCTION_CARTESIAN : 
	      //Reduced Stencil 1D
	      //W = HALF*(SolnBlk.W[i][j] + SolnBlk.W[i+1][j]);

               W = SolnBlk.W[i+1][j] - (SolnBlk.W[i+1][j] - SolnBlk.W[i][j]) * 
                  (SolnBlk.Grid.Cell[i+1][j].Xc.x-SolnBlk.Grid.xfaceE(i,j).x)
                  /(SolnBlk.Grid.Cell[i+1][j].Xc.x - SolnBlk.Grid.Cell[i][j].Xc.x);

	      dWdx = (SolnBlk.W[i+1][j] - SolnBlk.W[i][j])/
		fabs(SolnBlk.Grid.Cell[i][j].Xc.x - SolnBlk.Grid.Cell[i+1][j].Xc.x);
	      dWdy = HALF*(SolnBlk.dWdy[i][j] + SolnBlk.dWdy[i+1][j]);
	      //dWdy.Vacuum(); 
	      // Compute the viscous flux.
	      Flux -= Viscous_Flux_n(W,dWdx,dWdy,
                                     SolnBlk.Axisymmetric,
				     SolnBlk.Grid.xfaceE(i,j),
				     SolnBlk.Grid.nfaceE(i,j));
	      break;	      
	    case VISCOUS_RECONSTRUCTION_DIAMOND_PATH :
	      //specify the Quadrature point primitive variable solution and
	      W_face = HALF*(SolnBlk.Wn(i+1,j+1) +SolnBlk.Wn(i+1,j));
	      dWdx = SolnBlk.dWdx_faceE[i][j];
	      dWdy = SolnBlk.dWdy_faceE[i][j];
	      Flux -= Viscous_Flux_n(W_face,dWdx,dWdy,
				     SolnBlk.Axisymmetric,
				     SolnBlk.Grid.xfaceE(i,j),
				     SolnBlk.Grid.nfaceE(i,j));
	    
	      break;
	    }	
	  }

	  /*****************************************/
          /* Evaluate cell-averaged solution changes. */
	  SolnBlk.dUdt[i][j][k_residual] -= 
	    (Input_Parameters.CFL_Number*SolnBlk.dt[i][j])*Flux*SolnBlk.Grid.lfaceE(i, j)/SolnBlk.Grid.Cell[i][j].A;

	  SolnBlk.dUdt[i+1][j][k_residual] += 
	    (Input_Parameters.CFL_Number*SolnBlk.dt[i+1][j])*Flux*SolnBlk.Grid.lfaceW(i+1, j)/SolnBlk.Grid.Cell[i+1][j].A;
	  
	   /*****************************************
	    ************* SOURCE TERMS **************
	    *****************************************/

	  /* Include axisymmetric source terms as required. */
	  if (SolnBlk.Axisymmetric) {
	    SolnBlk.dUdt[i][j][k_residual] += 
	      (Input_Parameters.CFL_Number*SolnBlk.dt[i][j])*
              SolnBlk.W[i][j].Sa_inviscid(SolnBlk.Grid.Cell[i][j].Xc,
                                          SolnBlk.Axisymmetric);
	    // Include Viscous if specified	 
	    if (SolnBlk.Flow_Type != FLOWTYPE_INVISCID) { 
	      SolnBlk.dUdt[i][j][k_residual] += 
		(Input_Parameters.CFL_Number*SolnBlk.dt[i][j])*
                SolnBlk.W[i][j].Sa_viscous(SolnBlk.dWdx[i][j],
					   SolnBlk.dWdy[i][j],
					   SolnBlk.Grid.Cell[i][j].Xc,
                                           SolnBlk.Axisymmetric);
	    }
          } /* endif */

          /* Include source terms associated with turbulence model */
          if (SolnBlk.Flow_Type == FLOWTYPE_TURBULENT_RANS_K_OMEGA ||
              SolnBlk.Flow_Type == FLOWTYPE_TURBULENT_RANS_K_EPSILON) {
	    SolnBlk.dUdt[i][j][k_residual] += 
	      (Input_Parameters.CFL_Number*SolnBlk.dt[i][j])*
	      SolnBlk.W[i][j].S_turbulence_model(SolnBlk.dWdx[i][j],
						 SolnBlk.dWdy[i][j],
						 SolnBlk.Grid.Cell[i][j].Xc,
						 SolnBlk.Axisymmetric);
          } /* endif */
	  /* Include source terms associated with the finite-rate chemistry and 
             turbulence/chemistry interactions */ 
	  if (SolnBlk.W[i][j].React.reactset_flag != NO_REACTIONS) {	 
	    //rho*omega_dot
	    //
	    SolnBlk.dUdt[i][j][k_residual] += (Input_Parameters.CFL_Number*SolnBlk.dt[i][j])* 
	      SolnBlk.W[i][j].Sw(SolnBlk.W[i][j].React.reactset_flag,SolnBlk.W[i][j].flow_type);
	    //Used as ignitor (NOT FINISHED YET)
	    //if(SolnBlk.Heat_source){
	    //  SolnBlk.dUdt[i][j][k_residual].E = Heat_source(SolnBlk.Grid.Cell[i][j].Xc);
	    //}     
	  }
	    
          /* Include source terms associated with gravity */
	  if (SolnBlk.Gravity) {	 
	      SolnBlk.dUdt[i][j][k_residual] += (Input_Parameters.CFL_Number*SolnBlk.dt[i][j])*
                                                SolnBlk.W[i][j].Sg();
          } /* endif */

	  /* Save west and east face boundary flux. */
	  // USED for AMR 	
	  if (i == SolnBlk.ICl-1) {
	    SolnBlk.FluxW[j] = -Flux*SolnBlk.Grid.lfaceW(i+1, j);
	  } else if (i == SolnBlk.ICu) {
	    SolnBlk.FluxE[j] = Flux*SolnBlk.Grid.lfaceE(i, j);
	  } 
	
	 //	  cout<<"\n i "<<i<<" "<<j<<" "<<k_residual<<" "<<SolnBlk.dUdt[i][j][k_residual]; 	

	} /* endif */
      } /* endfor */
                
      if ( j > SolnBlk.JCl-1 && j < SolnBlk.JCu+1 ){
	SolnBlk.dUdt[SolnBlk.ICu+1][j][k_residual].Vacuum();
	SolnBlk.dUdt[SolnBlk.ICl-1][j][k_residual].Vacuum();
      }
    }

    // Add j-direction (eta-direction) fluxes.
    for ( i = SolnBlk.ICl ; i <= SolnBlk.ICu ; ++i ) {
      for ( j  = SolnBlk.JCl-1 ; j <= SolnBlk.JCu ; ++j ) {

	/* Evaluate the cell interface j-direction fluxes. */ 
	if (j == SolnBlk.JCl-1 && 
	    (SolnBlk.Grid.BCtypeS[i] == BC_REFLECTION ||
	     SolnBlk.Grid.BCtypeS[i] == BC_CHARACTERISTIC ||
	     SolnBlk.Grid.BCtypeS[i] == BC_FREE_SLIP ||
	     SolnBlk.Grid.BCtypeS[i] == BC_NO_SLIP ||
	     SolnBlk.Grid.BCtypeS[i] == BC_MOVING_WALL ||
	     SolnBlk.Grid.BCtypeS[i] == BC_FLAME_OUTFLOW )) {

	  dX = SolnBlk.Grid.xfaceS(i, j+1)-SolnBlk.Grid.Cell[i][j+1].Xc;
	  Wr = SolnBlk.W[i][j+1] +
	    (SolnBlk.phi[i][j+1]^SolnBlk.dWdx[i][j+1])*dX.x +
	    (SolnBlk.phi[i][j+1]^SolnBlk.dWdy[i][j+1])*dX.y;
	  if (SolnBlk.Grid.BCtypeS[i] == BC_REFLECTION) {
	    Wl = Reflect(Wr, SolnBlk.Grid.nfaceS(i, j+1)); 
	  } else if (SolnBlk.Grid.BCtypeS[i] == BC_FREE_SLIP) {
	    Wl = Free_Slip(Wr,SolnBlk.WoS[i], SolnBlk.Grid.nfaceS(i, j+1),FIXED_TEMPERATURE_WALL);
	  } else if (SolnBlk.Grid.BCtypeS[i] == BC_NO_SLIP) {
	    Wl = No_Slip(Wr,SolnBlk.WoS[i], SolnBlk.Grid.nfaceS(i, j+1),FIXED_TEMPERATURE_WALL);
	  } else if (SolnBlk.Grid.BCtypeS[i] == BC_MOVING_WALL) {
	    Wl = Moving_Wall(Wr,SolnBlk.WoS[i], SolnBlk.Grid.nfaceS(i, j+1),
			     SolnBlk.Moving_wall_velocity,FIXED_TEMPERATURE_WALL);
	  } else if (SolnBlk.Grid.BCtypeS[i] == BC_FLAME_OUTFLOW){
	    Wl = BC_Flame_Outflow(Wr, 
				  SolnBlk.WoS[i], 
				  SolnBlk.W[i][SolnBlk.JCu],
				  SolnBlk.Grid.nfaceS(i, j+1));
	  } else { 
	    Wl = BC_Characteristic_Pressure(Wr, 
					    SolnBlk.WoS[i], 
					    SolnBlk.Grid.nfaceS(i, j+1));
	  } /* endif */

	} else if (j == SolnBlk.JCu && 
		   (SolnBlk.Grid.BCtypeN[i] == BC_REFLECTION ||
		    SolnBlk.Grid.BCtypeN[i] == BC_CHARACTERISTIC ||
		    SolnBlk.Grid.BCtypeN[i] == BC_FREE_SLIP ||
		    SolnBlk.Grid.BCtypeN[i] == BC_NO_SLIP ||
		    SolnBlk.Grid.BCtypeN[i] == BC_MOVING_WALL ||
		    SolnBlk.Grid.BCtypeN[i] == BC_FLAME_OUTFLOW )) {
	  
	  dX = SolnBlk.Grid.xfaceN(i, j)-SolnBlk.Grid.Cell[i][j].Xc;
	  Wl = SolnBlk.W[i][j] + 
	    (SolnBlk.phi[i][j]^SolnBlk.dWdx[i][j])*dX.x +
	    (SolnBlk.phi[i][j]^SolnBlk.dWdy[i][j])*dX.y;
	  if (SolnBlk.Grid.BCtypeN[i] == BC_REFLECTION) {
	    Wr = Reflect(Wl, SolnBlk.Grid.nfaceN(i, j));	
	  } else if (SolnBlk.Grid.BCtypeN[i] == BC_FREE_SLIP) {
	    Wr = Free_Slip(Wl, SolnBlk.WoN[i], SolnBlk.Grid.nfaceN(i, j),FIXED_TEMPERATURE_WALL);
	  } else if (SolnBlk.Grid.BCtypeN[i] == BC_NO_SLIP) {
	    Wr = No_Slip(Wl, SolnBlk.WoN[i], SolnBlk.Grid.nfaceN(i, j),FIXED_TEMPERATURE_WALL);
	  } else if (SolnBlk.Grid.BCtypeN[i] == BC_MOVING_WALL) {
	    Wr = Moving_Wall(Wl,SolnBlk.WoN[i],  SolnBlk.Grid.nfaceN(i, j),
			     SolnBlk.Moving_wall_velocity,FIXED_TEMPERATURE_WALL);
	  } else if (SolnBlk.Grid.BCtypeN[i] == BC_FLAME_OUTFLOW){
	    Wr = BC_Flame_Outflow(Wl, 
				  SolnBlk.WoN[i], 
				  SolnBlk.W[i][SolnBlk.JCl],
				  SolnBlk.Grid.nfaceW(i, j));
	  } else {
	    Wr = BC_Characteristic_Pressure(Wl, 
					    SolnBlk.WoN[i], 
					    SolnBlk.Grid.nfaceN(i, j));
	  } /* endif */
	  
	} else {
	  dX = SolnBlk.Grid.xfaceN(i, j)-SolnBlk.Grid.Cell[i][j].Xc;
	  Wl = SolnBlk.W[i][j] + 
	    (SolnBlk.phi[i][j]^SolnBlk.dWdx[i][j])*dX.x +
	    (SolnBlk.phi[i][j]^SolnBlk.dWdy[i][j])*dX.y;
	  dX = SolnBlk.Grid.xfaceS(i, j+1)-SolnBlk.Grid.Cell[i][j+1].Xc;
	  Wr = SolnBlk.W[i][j+1] +
	    (SolnBlk.phi[i][j+1]^SolnBlk.dWdx[i][j+1])*dX.x +
	    (SolnBlk.phi[i][j+1]^SolnBlk.dWdy[i][j+1])*dX.y;
	} /* endif */
		
	// Spacing for Preconditioner 
	if(SolnBlk.Flow_Type != FLOWTYPE_INVISCID && 
           Input_Parameters.Preconditioning){
	  delta_n = min( TWO*(SolnBlk.Grid.Cell[i][j].A/
			      (SolnBlk.Grid.lfaceE(i, j)+SolnBlk.Grid.lfaceW(i, j))),
			 TWO*(SolnBlk.Grid.Cell[i][j].A/
			      (SolnBlk.Grid.lfaceN(i, j)+SolnBlk.Grid.lfaceS(i, j))));
	}
	switch(Input_Parameters.i_Flux_Function) {
	  //             case FLUX_FUNCTION_GODUNOV :
	  //               Flux = FluxGodunov_n(Wl, Wr, SolnBlk.Grid.nfaceN(i, j));
	  //               break;
	case FLUX_FUNCTION_ROE :

	  Flux = FluxRoe_n(Wl, Wr, SolnBlk.Grid.nfaceN(i, j),
			   Input_Parameters.Preconditioning,delta_n);
	  break;
	  //             case FLUX_FUNCTION_RUSANOV :
	  //               Flux = FluxRusanov_n(Wl, Wr, SolnBlk.Grid.nfaceN(i, j));
	  //               break;
	case FLUX_FUNCTION_HLLE :
	  Flux = FluxHLLE_n(Wl, Wr, SolnBlk.Grid.nfaceN(i, j),Input_Parameters.Preconditioning);
	  break;
	case FLUX_FUNCTION_LINDE :
	  Flux = FluxLinde_n(Wl, Wr, SolnBlk.Grid.nfaceN(i, j));
	  break;
	  //             case FLUX_FUNCTION_HLLC :
	  //               Flux = FluxHLLC_n(Wl, Wr, SolnBlk.Grid.nfaceN(i, j));
	  //               break;
	  //             default:
	  //               Flux = FluxRoe_n(Wl, Wr, SolnBlk.Grid.nfaceN(i, j));
	  //               break;
	} /* endswitch */
	

	/*****************************************/
	/***************** VISCOUS ***************/
	//eta  
	if(SolnBlk.Flow_Type != FLOWTYPE_INVISCID){ // -ve as on RHS    
	  switch(Input_Parameters.i_Viscous_Flux_Evaluation){	    
	  case VISCOUS_RECONSTRUCTION_ARITHMETIC_AVERAGE :
	    //Arithmetic Mean
	    Flux -= Viscous_FluxArithmetic_n(SolnBlk.U[i][j],
                                             SolnBlk.dWdx[i][j],
                                             SolnBlk.dWdy[i][j],
                                             SolnBlk.U[i][j+1],
                                             SolnBlk.dWdx[i][j+1],
                                             SolnBlk.dWdy[i][j+1],
	     			             SolnBlk.Grid.nfaceN(i,j));  // d/dz Fzv
	    break;
	  case VISCOUS_RECONSTRUCTION_CARTESIAN :
	    //Reduced Stencil 1D
	    //W = HALF*(SolnBlk.W[i][j] + SolnBlk.W[i][j+1]);   
	    W = SolnBlk.W[i][j+1] - (SolnBlk.W[i][j+1] - SolnBlk.W[i][j]) * 
                (SolnBlk.Grid.Cell[i][j+1].Xc.y - SolnBlk.Grid.xfaceN(i,j).y)
                /(SolnBlk.Grid.Cell[i][j+1].Xc.y - SolnBlk.Grid.Cell[i][j].Xc.y);
	    
	    dWdx = HALF*(SolnBlk.dWdx[i][j] + SolnBlk.dWdx[i][j+1]);
	    dWdy = (SolnBlk.W[i][j+1] - SolnBlk.W[i][j])/
	      fabs(SolnBlk.Grid.Cell[i][j+1].Xc.y - SolnBlk.Grid.Cell[i][j].Xc.y);
	    // Compute the viscous flux.
	    Flux -= Viscous_Flux_n(W,dWdx,dWdy,
                                   SolnBlk.Axisymmetric,
				   SolnBlk.Grid.xfaceN(i,j),
				   SolnBlk.Grid.nfaceN(i,j));
	 
	    break;	  
	  case VISCOUS_RECONSTRUCTION_DIAMOND_PATH :
	    //specify the Quadrature poiny primitive variable solution and
	    // primitive variable gradients on the north face
	    W_face = HALF*(SolnBlk.Wn(i,j+1) + SolnBlk.Wn(i+1,j+1));
	    dWdx = SolnBlk.dWdx_faceN[i][j];
	    dWdy = SolnBlk.dWdy_faceN[i][j];
	    Flux -= Viscous_Flux_n(W_face,dWdx,dWdy,
                                   SolnBlk.Axisymmetric,
				   SolnBlk.Grid.xfaceN(i,j),
				   SolnBlk.Grid.nfaceN(i,j));
	    break;
	  }
	}
	
	/* Evaluate cell-averaged solution changes. */
	SolnBlk.dUdt[i][j][k_residual] -= 
	  (Input_Parameters.CFL_Number*SolnBlk.dt[i][j])
	  *Flux*SolnBlk.Grid.lfaceN(i,j)/SolnBlk.Grid.Cell[i][j].A;

	SolnBlk.dUdt[i][j+1][k_residual] += 
	  (Input_Parameters.CFL_Number*SolnBlk.dt[i][j+1])
	  *Flux*SolnBlk.Grid.lfaceS(i,j+1)/SolnBlk.Grid.Cell[i][j+1].A;

	/* Save south and north face boundary flux. */
	// USED for AMR 
	if (j == SolnBlk.JCl-1) {
	  SolnBlk.FluxS[i] = -Flux*SolnBlk.Grid.lfaceS(i, j+1); 
	} else if (j == SolnBlk.JCu) {
	  SolnBlk.FluxN[i] = Flux*SolnBlk.Grid.lfaceN(i, j);  
	} /* endif */	

	//	cout<<"\n j "<<i<<" "<<j<<" "<<k_residual<<" "<<SolnBlk.dUdt[i][j][k_residual]; 	

      } /* endfor */

      SolnBlk.dUdt[i][SolnBlk.JCu+1][k_residual].Vacuum();
      SolnBlk.dUdt[i][SolnBlk.JCl-1][k_residual].Vacuum();

    } /* endfor */
  
    // For k-omega turbulence model set residual in laminar sublayer to zero.
    if (SolnBlk.Flow_Type == FLOWTYPE_TURBULENT_RANS_K_OMEGA) {
       for ( i = SolnBlk.ICl ; i <= SolnBlk.ICu ; ++i ) {
          for ( j  = SolnBlk.JCl ; j <= SolnBlk.JCu ; ++j ) {
	     if (SolnBlk.Wall[i][j].yplus <= Input_Parameters.yplus_sublayer) {
	        SolnBlk.dUdt[i][j][k_residual].rhoomega = ZERO;
             } /* endif */
          } /* endfor */
       } /* endfor */
    } /* endif */

    /* Residual successfully calculated. */

    //   cout<<"\n New Block \n";


    return (0);
    
}

/********************************************************
 * Routine: Update_Solution_Multistage_Explicit         *
 *                                                      *
 * This routine updates solution states of the given    *
 * solution block for a variety of multi-stage explicit *
 * time integration schemes.                            *
 *                                                      *
 ********************************************************/
int Update_Solution_Multistage_Explicit(Chem2D_Quad_Block &SolnBlk,
                                        const int i_stage,
                                        Chem2D_Input_Parameters &Input_Parameters) {


  int i, j, k_residual;
  double omega;
  int NUM_VAR_CHEM2D = SolnBlk.NumVar(); 

  double delta_n;

  // Memory for linear system solver. 
  DenseSystemLinEqs LinSys;
  DenseMatrix dSdU(NUM_VAR_CHEM2D-1,NUM_VAR_CHEM2D-1); //Source Jacobian
  DenseMatrix dRdU(NUM_VAR_CHEM2D-1,NUM_VAR_CHEM2D-1); //Residual Jacobian
  DenseMatrix Precon(NUM_VAR_CHEM2D-1,NUM_VAR_CHEM2D-1); // Low Mach number preconditioner
  DenseMatrix Precon_Inv(NUM_VAR_CHEM2D-1,NUM_VAR_CHEM2D-1); // Inverse of low Mach number preconditioner
  
  /* Allocate memory for linear system solver. */
  LinSys.allocate(NUM_VAR_CHEM2D-1);
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
   
    /* Perform update of solution variables for stage i_stage of n_stage scheme. */
    /* Update solution variables for this stage. */
 
  for ( j  = SolnBlk.JCl ; j <= SolnBlk.JCu ; ++j ) {
    for ( i = SolnBlk.ICl ; i <= SolnBlk.ICu ; ++i ) {
      
      /******************************************************/
      /********** Fully Explicit ****************************/
      /******************************************************/
      if (Input_Parameters.Local_Time_Stepping == GLOBAL_TIME_STEPPING || 
	  Input_Parameters.Local_Time_Stepping == SCALAR_LOCAL_TIME_STEPPING ) {
	//Update
	SolnBlk.U[i][j] = SolnBlk.Uo[i][j] + omega*SolnBlk.dUdt[i][j][k_residual];
	//N-1 species
	SolnBlk.U[i][j][NUM_VAR_CHEM2D] = SolnBlk.U[i][j].rho*(ONE - SolnBlk.U[i][j].sum_species());	   
      }
      
      //Check for unphysical properties  // Xinfeng, write valid solution check!
      /**********************************************************/
      /* If unphysical properties and using global timestepping */ 
      /* stop simulation                                        */
      /**********************************************************/   
      if (Input_Parameters.Local_Time_Stepping == GLOBAL_TIME_STEPPING) {
	if(!SolnBlk.U[i][j].Unphysical_Properties_Check(10)) return (i);
	/*********************************************************/
	/* If unphysical properties and using local timestepping */ 
	/* try reducing step size                                */
	/*********************************************************/    
      } else if (Input_Parameters.Local_Time_Stepping == SCALAR_LOCAL_TIME_STEPPING) {
	if( !SolnBlk.U[i][j].Unphysical_Properties_Check(10)) {  
	  
	  if (SolnBlk.debug_level) {
	    cout<<"\n Using local time stepping step reduction in Update_Solution_Multistage_Explicit\n";
	  }
	  
	  double residual_reduction_factor = ONE;
	  for (int n_residual_reduction = 1; n_residual_reduction <= 10; ++n_residual_reduction) {
	    if(SolnBlk.debug_level) {
	      cout<<".."<<n_residual_reduction<<"..";
	    }
	    residual_reduction_factor = HALF*residual_reduction_factor;
	    SolnBlk.dt[i][j] = residual_reduction_factor*SolnBlk.dt[i][j];
	    SolnBlk.dUdt[i][j][k_residual] = residual_reduction_factor*SolnBlk.dUdt[i][j][k_residual];	     
	    
	    //Update
	    SolnBlk.U[i][j] = SolnBlk.Uo[i][j] + omega*SolnBlk.dUdt[i][j][k_residual];
	    //N-1 species
	    SolnBlk.U[i][j][NUM_VAR_CHEM2D] = SolnBlk.U[i][j].rho*(ONE - SolnBlk.U[i][j].sum_species());
	    
	    if(SolnBlk.U[i][j].Unphysical_Properties_Check(n_residual_reduction))  break;
	  } /* endfor */  
	  if (Input_Parameters.Local_Time_Stepping == 1 && 
	      !SolnBlk.U[i][j].Unphysical_Properties_Check(10)
	      ) return(i);  	      	   
	}
      } /* endif */
      
      /**************************************************************/
      /************ SEMI-IMPLICIT AND/OR PRECONDITIONING ************/
      /**************************************************************/
      if (Input_Parameters.Local_Time_Stepping == SEMI_IMPLICIT_LOCAL_TIME_STEPPING || 
	  Input_Parameters.Local_Time_Stepping == LOW_MACH_NUMBER_WEISS_SMITH_PRECONDITIONER || 
	  Input_Parameters.Local_Time_Stepping == SEMI_IMPLICIT_LOW_MACH_NUMBER_PRECONDITIONER || 
	  Input_Parameters.Local_Time_Stepping == MATRIX_LOCAL_TIME_STEPPING) {
	
	dSdU.zero();
	dRdU.zero();

	/************ FORM LHS FOR SEMI-IMPLICIT METHOD *****************/
	if (Input_Parameters.Local_Time_Stepping == SEMI_IMPLICIT_LOCAL_TIME_STEPPING){

	  /* Semi implicit formulation: set up system of equations 
	     and include source Jacobian in the LHS matrix. */
	  
	  dSdU = SemiImplicitBlockJacobi(SolnBlk,i, j);	  	  
	  LinSys.A.identity();
	     
	  //scalar multiplication
	  LinSys.A -= (omega*Input_Parameters.CFL_Number*SolnBlk.dt[i][j])*dSdU;	     
	 
	}
	
	// Spacing for preconditioner and viscous 
	if (SolnBlk.Flow_Type != FLOWTYPE_INVISCID && 
	    (Input_Parameters.Local_Time_Stepping == LOW_MACH_NUMBER_WEISS_SMITH_PRECONDITIONER || 
	     Input_Parameters.Local_Time_Stepping == SEMI_IMPLICIT_LOW_MACH_NUMBER_PRECONDITIONER) ){	     
	  delta_n = min( TWO*(SolnBlk.Grid.Cell[i][j].A/
			      (SolnBlk.Grid.lfaceE(i, j)+SolnBlk.Grid.lfaceW(i, j))),
			 TWO*(SolnBlk.Grid.Cell[i][j].A/
			      (SolnBlk.Grid.lfaceN(i, j)+SolnBlk.Grid.lfaceS(i, j))));
	}
	
	/************ FORM LHS FOR LOW MACH NUMBER PRECONDITIONING ***************/
	if (Input_Parameters.Local_Time_Stepping == LOW_MACH_NUMBER_WEISS_SMITH_PRECONDITIONER ){	   

	  SolnBlk.Uo[i][j].Low_Mach_Number_Preconditioner(Precon,
							  delta_n);	 
	  LinSys.A.zero();
	  // not Multiplied by omega*CFL_Number*SolnBlk.dt[i][j] as dimensionless
	  LinSys.A = Precon;

	} 

	/******* FORM LHS FOR SEMI-IMPLICIT METHOD WITH PRECONDITIONING ********/
	//LHS = P*( I- Pinv*(dt*(dS/dU)))	     
	if (Input_Parameters.Local_Time_Stepping == SEMI_IMPLICIT_LOW_MACH_NUMBER_PRECONDITIONER ){
	 	  
	  /* Semi implicit formulation: set up system of equations 
	     and include source Jacobian in the LHS matrix. */
	  dSdU = SemiImplicitBlockJacobi(SolnBlk,i, j);
	  LinSys.A.identity();
	  
	  //Currently recalculating Pinv and P, expensive !!! should store between	     
	  //stages of time method ( as based on Uo, not U)
	  //Inverse of Preconditioner
	  SolnBlk.U[i][j].Low_Mach_Number_Preconditioner_Inverse(Precon_Inv,
								 delta_n);
	  
	  //I - Pinv*(dt*(dS/dU))
	  LinSys.A -= Precon_Inv*((omega*Input_Parameters.CFL_Number*SolnBlk.dt[i][j])*dSdU);
	  
	  //Preconditioner 
	  SolnBlk.Uo[i][j].Low_Mach_Number_Preconditioner(Precon,
							  delta_n);
	  
	  //LHS = P*( I - Pinv*(dt*(dS/dU)))	     
	  LinSys.A = Precon*LinSys.A;	    

 
	}//end SEMI_IMPLICIT_LOW_MACH_NUMBER_PRECONDITIONER

	/************ FORM LHS FOR MATRIX LOCAL TIME STEPPING BASED ON *****************
	 ************ POINT IMPLICIT BLOCK JACOBI PRECONDITIONER       *****************/
	if (Input_Parameters.Local_Time_Stepping == MATRIX_LOCAL_TIME_STEPPING) {

	  dRdU = PointImplicitBlockJacobi(SolnBlk,Input_Parameters,i, j);

	  //	  LinSys.A  = -( 1.00*SolnBlk.dt[i][j])*dRdU; 
	  //Talk to Prof. Groth about ... 
	  // Same thing as before the follwoing formula works, but not "theoretically" correct ...
	  LinSys.A.identity();
	  LinSys.A -=(omega*Input_Parameters.CFL_Number*SolnBlk.dt[i][j])*dRdU; 

	}

	/******************************************************************
	 ******************* EVALUATE RHS *********************************
	 ******************************************************************/
	// omega*(dt*RHS)
	for (int k = 0; k < NUM_VAR_CHEM2D - 1; k++) {
	  LinSys.b(k) = omega*SolnBlk.dUdt[i][j][k_residual][k+1];
	}
	/* Solve system of equations using LU decomposition Gaussian elimination procedure. */
	LinSys.solve(LU_DECOMPOSITION);
	
	/* Update the conserved solution variables. */
	for (int k = 0; k < NUM_VAR_CHEM2D - 1; k++) {
	  SolnBlk.U[i][j][k+1] = SolnBlk.Uo[i][j][k+1] + LinSys.x(k);
	} 

	/*********************************************************
	     Using N-1 species equations so need to update 
             last speceis using:
	     c_n = 1- sum(1 to N-1) cs
	*********************************************************/
	SolnBlk.U[i][j][NUM_VAR_CHEM2D] = SolnBlk.U[i][j].rho*(ONE - SolnBlk.U[i][j].sum_species());

	/*-----------------------------------------------------------------*/
	/* Apply low-Reynolds number formulations and wall functions for   */
	/* turbulent flows                                       */
  	if (SolnBlk.Flow_Type == FLOWTYPE_TURBULENT_RANS_K_EPSILON ||
	    SolnBlk.Flow_Type == FLOWTYPE_TURBULENT_RANS_K_OMEGA) {
	  if(Input_Parameters.i_Turbulence_BCs == TURBULENT_BC_DIRECT_INTEGRATION) Low_ReynoldsNumber_Formulation(SolnBlk, Input_Parameters, i, j);
	}

	/*********************************************************/
	/* If unphysical properties and using local timestepping */ 
	/* try reducing step size                                */
	/*********************************************************/
	
	if(!SolnBlk.U[i][j].Unphysical_Properties_Check(10)){
	  
	  if (SolnBlk.debug_level) {
	    cout<<"\n Using local time stepping step reduction in Update_Solution_Multistage_Explicit\n";
	  }
	  
	  double residual_reduction_factor = ONE;
	  for (int n_residual_reduction = 1; n_residual_reduction <= 10; ++n_residual_reduction) {
	    if (SolnBlk.debug_level) {
	      cout<<".."<<n_residual_reduction<<"..";
	    }
	    residual_reduction_factor = HALF*residual_reduction_factor;
	    SolnBlk.dt[i][j] = residual_reduction_factor*SolnBlk.dt[i][j];
	    SolnBlk.dUdt[i][j][k_residual] = residual_reduction_factor*SolnBlk.dUdt[i][j][k_residual];
	    
	    if (Input_Parameters.Local_Time_Stepping == SEMI_IMPLICIT_LOCAL_TIME_STEPPING ){
	      LinSys.A.identity();
	      LinSys.A -= (omega*Input_Parameters.CFL_Number*SolnBlk.dt[i][j])*dSdU;
	    }
	    
	    if (Input_Parameters.Local_Time_Stepping == SEMI_IMPLICIT_LOW_MACH_NUMBER_PRECONDITIONER){		 
	      LinSys.A.identity();		 
	      LinSys.A -= Precon_Inv*((residual_reduction_factor*omega*Input_Parameters.CFL_Number*SolnBlk.dt[i][j])*dSdU);
	      //LHS = P*( I - Pinv*(dt*(dS/dU)))	     
	      LinSys.A = Precon*LinSys.A;	     
	    }
	    
	    if (Input_Parameters.Local_Time_Stepping == MATRIX_LOCAL_TIME_STEPPING) {
	      dRdU = PointImplicitBlockJacobi(SolnBlk,Input_Parameters,i, j);
	      //  LinSys.A = -(1.00*SolnBlk.dt[i][j])*dRdU;//Markthis
	      //Talk to Prof. Groth about ... 
	      // Same thing as before the follwoing formula works, but not "theoretically" correct ...
	      LinSys.A.identity();
	      LinSys.A -=(omega*Input_Parameters.CFL_Number*SolnBlk.dt[i][j])*dRdU; 
	      
	    }
	    
	    for (int k = 0; k < NUM_VAR_CHEM2D-1; k++) {
	      LinSys.b(k) = omega*SolnBlk.dUdt[i][j][k_residual][k+1];
	    } 
	    
	    LinSys.solve(LU_DECOMPOSITION);
	    
	    /* Update the conserved Solution */
	    for (int k = 0; k < NUM_VAR_CHEM2D-1; k++) {
	      SolnBlk.U[i][j][k+1] = SolnBlk.Uo[i][j][k+1] + LinSys.x(k);
	    } 
	    
	    SolnBlk.U[i][j][NUM_VAR_CHEM2D] = SolnBlk.U[i][j].rho*(ONE - SolnBlk.U[i][j].sum_species());
	    
	    if (SolnBlk.Flow_Type == FLOWTYPE_TURBULENT_RANS_K_EPSILON ||
		SolnBlk.Flow_Type == FLOWTYPE_TURBULENT_RANS_K_OMEGA) {
	      if(Input_Parameters.i_Turbulence_BCs == TURBULENT_BC_DIRECT_INTEGRATION) Low_ReynoldsNumber_Formulation(SolnBlk, Input_Parameters, i, j);
	    }
	    if (SolnBlk.U[i][j].Unphysical_Properties_Check(n_residual_reduction))  break;
	  }
	} 
	/*********************************************************/
	//Check for unphysical properties
	if(!SolnBlk.U[i][j].Unphysical_Properties_Check(10)) return (i);

      } //end implicit and/or preconditioned formulations 

      //Calculate Primitive values from updated conserved solution
      SolnBlk.W[i][j] = W(SolnBlk.U[i][j]); 

    } //end i
  } //end j //

  LinSys.deallocate();

  /* Solution successfully updated. */
  return (0);

}

/**********************************************************************
 *   Viscous Calculations - Cell Centered                             *
 *                                                                    *
 *       -> U.qflux                                                   *
 *	 -> U.tau                                                     *
 *	 -> U.rhospec[i].diffusion                                    *
 *	 -> U.rhospec[i].gradc                                        *
 *	 -> local grad Temperature                                    *
 *                                                                    *
 **********************************************************************/
void Viscous_Calculations(Chem2D_Quad_Block &SolnBlk) {

  double mu, kappa, r, div_v;
  double Temperature, Rmix, Cp;
  double mu_t, kappa_t, Dm_t;
  Vector2D grad_T;
  Tensor2D strain_rate;

  for (int j = SolnBlk.JCl-1; j <= SolnBlk.JCu+1; j++) {
    for (int i = SolnBlk.ICl-1; i <= SolnBlk.ICu+1; i++) {

      /* Determine the molecular transport properties. */      
      mu = SolnBlk.W[i][j].mu();
      kappa = SolnBlk.W[i][j].kappa();
      Temperature = SolnBlk.W[i][j].T();
      Rmix = SolnBlk.W[i][j].Rtot();
      Cp =  SolnBlk.W[i][j].Cp();

      /* Determine the turbulence model transport properties. */
      if (SolnBlk.Flow_Type == FLOWTYPE_TURBULENT_RANS_K_OMEGA ||
          SolnBlk.Flow_Type == FLOWTYPE_TURBULENT_RANS_K_EPSILON) {
         mu_t = SolnBlk.W[i][j].eddy_viscosity();
         kappa_t =  mu_t*Cp/SolnBlk.W[i][j].Pr_turb();
         Dm_t = SolnBlk.W[i][j].Dm_turb();
      } /* endif */
  
      /***************** Strain rate (+dilatation) **********************/	
      div_v = SolnBlk.dWdx[i][j].v.x + SolnBlk.dWdy[i][j].v.y;
      if (SolnBlk.Axisymmetric == 2) {
         r = SolnBlk.Grid.Cell[i][j].Xc.x;
         div_v += SolnBlk.W[i][j].v.x/r;
      } else if (SolnBlk.Axisymmetric == 1) {
         r = SolnBlk.Grid.Cell[i][j].Xc.y;
         div_v += SolnBlk.W[i][j].v.y/r;
      } /* endif */
      strain_rate.xx = SolnBlk.dWdx[i][j].v.x-div_v/THREE;
      strain_rate.xy = HALF*(SolnBlk.dWdx[i][j].v.y + SolnBlk.dWdy[i][j].v.x);
      strain_rate.yy = SolnBlk.dWdy[i][j].v.y-div_v/THREE;
      if (SolnBlk.Axisymmetric == 0) {
         strain_rate.zz = -(strain_rate.xx + strain_rate.yy); 
      } else if (SolnBlk.Axisymmetric == 2) {
         strain_rate.zz = SolnBlk.W[i][j].v.x/r-div_v/THREE;
      } else if (SolnBlk.Axisymmetric == 1) {
	 strain_rate.zz = SolnBlk.W[i][j].v.y/r-div_v/THREE;
      } /* endif */

      /**************** Temperature gradient ***************************/
      /* Temperature gradients from using the chain rule 
         with the ideal gas law (P=rho*R*T) 
         dT/dx = 1/rho*R *( dP/dx - P/rho * drho/dx) */
      grad_T.x = (ONE/(SolnBlk.W[i][j].rho*Rmix)) * (SolnBlk.dWdx[i][j].p - 
                 (SolnBlk.W[i][j].p/SolnBlk.W[i][j].rho)*SolnBlk.dWdx[i][j].rho);
      grad_T.y = (ONE/(SolnBlk.W[i][j].rho*Rmix)) * (SolnBlk.dWdy[i][j].p - 
                 (SolnBlk.W[i][j].p/SolnBlk.W[i][j].rho)*SolnBlk.dWdy[i][j].rho);

      /*************** Molecular (Laminar) Diffusion of Species ********/
      // for each of the "n" species
      for( int k=0; k<SolnBlk.U[i][j].ns; k++){
	/***************** Diffusion coefficients **********************/
	// using global Schmidt number relation Sc = mu/rho*Ds
	SolnBlk.U[i][j].rhospec[k].diffusion_coef = mu/SolnBlk.W[i][j].Schmidt[k];
	/***************** mass fraction gradients *********************/
	SolnBlk.U[i][j].rhospec[k].gradc.x = SolnBlk.U[i][j].rho * SolnBlk.dWdx[i][j].spec[k].c;
	SolnBlk.U[i][j].rhospec[k].gradc.y = SolnBlk.U[i][j].rho * SolnBlk.dWdy[i][j].spec[k].c;
      }

      /***************** Molecular (Laminar) Stresses ******************/
      SolnBlk.U[i][j].tau = (TWO*mu)*strain_rate;

      /********************** Turbulent Stresses ***********************/
      if (SolnBlk.Flow_Type == FLOWTYPE_TURBULENT_RANS_K_OMEGA ||
          SolnBlk.Flow_Type == FLOWTYPE_TURBULENT_RANS_K_EPSILON) {
         SolnBlk.U[i][j].lambda = (TWO*mu_t)*strain_rate;
	 SolnBlk.U[i][j].lambda.xx -= (TWO/THREE)*SolnBlk.U[i][j].rhok; 
         SolnBlk.U[i][j].lambda.yy -= (TWO/THREE)*SolnBlk.U[i][j].rhok; 
         SolnBlk.U[i][j].lambda.zz -= (TWO/THREE)*SolnBlk.U[i][j].rhok; 
      } /* endif */

      /************* Molecular (Laminar) Heat flux Vector **************/
      /****************** Thermal Conduction ***************************
         q = - kappa * grad(T)                                         */
      SolnBlk.U[i][j].qflux = - kappa*grad_T;
      /****************** Thermal Diffusion ****************************/
      // q -= rho * sum ( hs * Ds *gradcs)  
      SolnBlk.U[i][j].qflux -= SolnBlk.U[i][j].rho*SolnBlk.U[i][j].thermal_diffusion(Temperature);  
      
      /**************** Turbulent Heat flux Vector *********************/
      /****************** Thermal Conduction ***************************
         q = - kappa * grad(T)                                         */
      if (SolnBlk.Flow_Type == FLOWTYPE_TURBULENT_RANS_K_OMEGA ||
          SolnBlk.Flow_Type == FLOWTYPE_TURBULENT_RANS_K_EPSILON) {
        SolnBlk.U[i][j].theta = - kappa_t*grad_T;
      /****************** Thermal Diffusion ****************************/
      // q -= rho * sum ( hs * Ds *gradcs)   
	for( int k=0; k<SolnBlk.U[i][j].ns; k++){
	  SolnBlk.U[i][j].theta -= Dm_t*SolnBlk.U[i][j].rhospec[k].gradc*
	    (SolnBlk.U[i][j].specdata[k].Enthalpy(Temperature)+
	     SolnBlk.U[i][j].specdata[k].Heatofform());
	  

         }
      } /* endif */
      /****** Set Primitive *******/
      SolnBlk.W[i][j] = W(SolnBlk.U[i][j]);

    }     
  }
}


// Low_ReynoldsNumber_Formulation for k-omega model
extern void Low_ReynoldsNumber_Formulation(Chem2D_Quad_Block &SolnBlk, 
                                           Chem2D_Input_Parameters &Input_Parameters, 
                                           int i, int j){

  double y, nu;
  double Pressure_Gradient = ZERO;
  double tau_wall = ZERO;
  double uf =ZERO; //1.587;// sqrt(tau_wall/1.225);// 1.0776;
  double y_plus =ZERO;
  double omega_bc = ZERO;

  Pressure_Gradient = (SolnBlk.W[SolnBlk.ICl-SolnBlk.Nghost][SolnBlk.JCl-SolnBlk.Nghost].p - 
		       SolnBlk.W[SolnBlk.ICu+SolnBlk.Nghost][SolnBlk.JCu+SolnBlk.Nghost].p)/
                       Input_Parameters.Pipe_Length;

  //tau_wall =  -R/(j+1) *Pressure_Gradient;// (j=0 for planar, j=1 for axisymmetric flow) 
  tau_wall = -Input_Parameters.Pipe_Radius*Pressure_Gradient;
  if (SolnBlk.Axisymmetric) {
    tau_wall = HALF*tau_wall;
  }
  tau_wall = fabs(tau_wall);
  uf = sqrt(tau_wall/SolnBlk.W[i][j].rho);
  Vector2D X_cell = SolnBlk.Grid.Cell[i][j].Xc;
  y = Distance_to_Wall(SolnBlk, X_cell);
  nu = SolnBlk.W[i][j].mu()/SolnBlk.W[i][j].rho; //dynamic viscosity
  y_plus = uf*y/nu;

//  cout<<" i = "<<i <<" j= "<<j<<endl;
//  cout<<"Pressure Gradient =   "<<Pressure_Gradient <<endl;
//  cout<<"tau_wall =   "<<tau_wall <<endl;
//  cout<<"uf = "<<uf<<endl;
//  cout<< "y = "<<y<<" nu = "<<nu<<endl;
//  cout<<"y_plus "<<y_plus<<endl;

  if (y_plus <= Input_Parameters.yplus_sublayer) {
    omega_bc =SolnBlk.W[i][j].omega_sublayer_BC(y);
    SolnBlk.U[i][j].rhoomega = SolnBlk.W[i][j].rho*omega_bc;
  }

}
