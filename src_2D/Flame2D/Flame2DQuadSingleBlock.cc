/************** Flame2DQuadSingleBlock.cc ***************************
   Single-Block Versions of Subroutines for 2D Euler 
   Multi-Block Quadrilateral Mesh  Solution Classes. 

   NOTES:
          
   - based on Euler2DQuadSingleBlock.cc
**********************************************************************/

/* Include 2D Euler quadrilateral mesh solution header file. */
#include "Flame2DdRdU.h"

/* Include header file to compute flame jump conditions */
#include "../Reactions/FlameJump.h"

/*************************************************************************
* Flame2D_Quad_Block -- Single Block External Subroutines.                *
**************************************************************************/

/********************************************************
 * Routine: Write_Solution_Block                        *
 *                                                      *
 * Writes the cell centred solution values of the       *
 * specified quadrilateral solution block to the        *
 * specified output stream for restart purposes.        *
 *                                                      *
 ********************************************************/
void Write_Solution_Block(Flame2D_Quad_Block &SolnBlk,
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
void Read_Solution_Block(Flame2D_Quad_Block &SolnBlk,
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
void Broadcast_Solution_Block(Flame2D_Quad_Block &SolnBlk) {

#ifdef _MPI_VERSION

  int ni, nj, ng, nr, nn, block_allocated, buffer_size;
    double *buffer;

    int NUM_VAR_FLAME2D = SolnBlk.NumVar(); 

    /* Broadcast the number of cells in each direction. */
    if (CFFC_Primary_MPI_Processor()) {
      ni = SolnBlk.NCi;
      nj = SolnBlk.NCj;
      ng = SolnBlk.Nghost;
      nr = SolnBlk.residual_variable;
      nn = SolnBlk.Number_of_Residual_Norms;
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
    MPI::COMM_WORLD.Bcast(&nn, 1, MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(&block_allocated, 1, MPI::INT, 0);

    /* On non-primary MPI processors, allocate (re-allocate) 
       memory for the quadrilateral solution block as necessary. */

    if (!CFFC_Primary_MPI_Processor()) {
      if (SolnBlk.NCi != ni || SolnBlk.NCj != nj || SolnBlk.Nghost != ng ) { 
	if (SolnBlk.U != NULL) SolnBlk.deallocate();
	if (block_allocated) SolnBlk.allocate(ni-2*ng, nj-2*ng, ng); 
      }
      // Set the block static variables if they were not previously assigned.
      if (SolnBlk.residual_variable != nr) SolnBlk.residual_variable = nr;
      if (SolnBlk.Number_of_Residual_Norms != nn) SolnBlk.Number_of_Residual_Norms = nn;
    } 

    /* Broadcast the axisymmetric/planar flow, viscous, turbulent, and gravity indicators. */
   
    MPI::COMM_WORLD.Bcast(&(SolnBlk.Flow_Type), 1, MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(&(SolnBlk.Axisymmetric), 1, MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(&(SolnBlk.Gravity), 1, MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(&(SolnBlk.debug_level), 1, MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(&(SolnBlk.Freeze_Limiter), 1, MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(&(SolnBlk.Moving_wall_velocity), 1, MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(SolnBlk.Pressure_gradient), 1, MPI::DOUBLE, 0);
    
    /* Broadcast the grid. */

    Broadcast_Quad_Block(SolnBlk.Grid);

    /* Broadcast the solution state variables. */

    if (block_allocated) {
      ni = (SolnBlk.ICu+SolnBlk.Nghost) - (SolnBlk.ICl-SolnBlk.Nghost) + 1;
      nj = (SolnBlk.JCu+SolnBlk.Nghost) - (SolnBlk.JCl-SolnBlk.Nghost) + 1; 
      buffer = new double[NUM_VAR_FLAME2D*ni*nj];
    
      if (CFFC_Primary_MPI_Processor()) {
	buffer_size = 0;
	for (int j  = SolnBlk.JCl-SolnBlk.Nghost ; j <= SolnBlk.JCu+SolnBlk.Nghost ; ++j ) {
	  for (int i = SolnBlk.ICl-SolnBlk.Nghost ; i <= SolnBlk.ICu+SolnBlk.Nghost ; ++i ) {
	    //Changed for Flame2D
	    for ( int k = 0; k < NUM_VAR_FLAME2D; ++k) {	
	      buffer[buffer_size] = SolnBlk.W[i][j][k+1];
	      buffer_size = buffer_size + 1;
	    }	      
	  } /* endfor */
	} /* endfor */
      } /* endif */
    
       buffer_size = NUM_VAR_FLAME2D*ni*nj;
       MPI::COMM_WORLD.Bcast(buffer, buffer_size, MPI::DOUBLE, 0);

       if (!CFFC_Primary_MPI_Processor()) {
	 buffer_size = 0;
	 for (int j  = SolnBlk.JCl-SolnBlk.Nghost ; j <= SolnBlk.JCu+SolnBlk.Nghost ; ++j ) {
	   for (int i = SolnBlk.ICl-SolnBlk.Nghost ; i <= SolnBlk.ICu+SolnBlk.Nghost ; ++i ) {
	     for ( int k = 0 ; k < NUM_VAR_FLAME2D; ++ k) {	
	       SolnBlk.W[i][j][k+1]= buffer[buffer_size];
	       buffer_size = buffer_size + 1;
	     }
	   } /* endfor */
	 } /* endfor */
       } /* endif */
       
       delete []buffer; 
       buffer = NULL;

       ni = 1;
       nj = (SolnBlk.JCu+SolnBlk.Nghost) - (SolnBlk.JCl-SolnBlk.Nghost) + 1;
       buffer = new double[2*NUM_VAR_FLAME2D*ni*nj];

       if (CFFC_Primary_MPI_Processor()) {
	 buffer_size = 0;
         for (int j  = SolnBlk.JCl-SolnBlk.Nghost ; j <= SolnBlk.JCu+SolnBlk.Nghost ; ++j ) {
	   for ( int k = 0 ; k < NUM_VAR_FLAME2D; ++ k) {	
	     buffer[buffer_size]= SolnBlk.WoW[j][k+1];
	     buffer_size = buffer_size + 1;
	   }
	   for ( int k = 0 ; k < NUM_VAR_FLAME2D; ++ k) {
	     buffer[buffer_size]= SolnBlk.WoE[j][k+1];  
	     buffer_size = buffer_size + 1;
	   }
	 } /* endfor */
       } /* endif */

       buffer_size = 2*NUM_VAR_FLAME2D*ni*nj;
       MPI::COMM_WORLD.Bcast(buffer, buffer_size, MPI::DOUBLE, 0);

       if (!CFFC_Primary_MPI_Processor()) {
	 buffer_size = 0;
	 for (int j  = SolnBlk.JCl-SolnBlk.Nghost ; j <= SolnBlk.JCu+SolnBlk.Nghost ; ++j ) {    
	   for ( int k = 0 ; k < NUM_VAR_FLAME2D; ++ k) {	
	     SolnBlk.WoW[j][k+1] = buffer[buffer_size];
	     buffer_size = buffer_size + 1;
	   }
	   for ( int k = 0 ; k < NUM_VAR_FLAME2D; ++ k) {
	     SolnBlk.WoE[j][k+1] = buffer[buffer_size]; 
	     buffer_size = buffer_size + 1;
	   }
          } /* endfor */
       } /* endif */

       delete []buffer; 
       buffer = NULL;

       ni = (SolnBlk.ICu+SolnBlk.Nghost) - (SolnBlk.ICl-SolnBlk.Nghost) + 1;
       nj = 1;
       buffer = new double[2*NUM_VAR_FLAME2D*ni*nj];

       if (CFFC_Primary_MPI_Processor()) {
	 buffer_size = 0;
	 for (int i  = SolnBlk.JCl-SolnBlk.Nghost ; i <= SolnBlk.JCu+SolnBlk.Nghost ; ++i ) {   
	   for ( int k = 0 ; k < NUM_VAR_FLAME2D; ++ k) {	
	     buffer[buffer_size]= SolnBlk.WoS[i][k+1];	  
	     buffer_size = buffer_size + 1;
	   }
	   for ( int k = 0 ; k < NUM_VAR_FLAME2D; ++ k) {
	     buffer[buffer_size]= SolnBlk.WoN[i][k+1];
	     buffer_size = buffer_size + 1;
	   }
	 } /* endfor */
       } /* endif */

       buffer_size = 2*NUM_VAR_FLAME2D*ni*nj;
       MPI::COMM_WORLD.Bcast(buffer, buffer_size, MPI::DOUBLE, 0);

       if (!CFFC_Primary_MPI_Processor()) {
	 buffer_size = 0;
	 for (int i = SolnBlk.ICl-SolnBlk.Nghost ; i <= SolnBlk.ICu+SolnBlk.Nghost ; ++i ) {     
	   for ( int k = 0 ; k < NUM_VAR_FLAME2D; ++ k) {	
	     SolnBlk.WoS[i][k+1] = buffer[buffer_size]; 
	     buffer_size = buffer_size + 1;
	    }
	   for ( int k = 0 ; k < NUM_VAR_FLAME2D; ++ k) {
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
void Broadcast_Solution_Block(Flame2D_Quad_Block &SolnBlk,
                             MPI::Intracomm &Communicator, 
                             const int Source_CPU) {

  int Source_Rank = 0;
  int ni, nj, ng, nr, nn, block_allocated, buffer_size;
  double *buffer;

  int NUM_VAR_FLAME2D = SolnBlk.NumVar();

  /* Broadcast the number of cells in each direction. */

  if (CFFC_MPI::This_Processor_Number == Source_CPU) {
    ni = SolnBlk.NCi;
    nj = SolnBlk.NCj;
    ng = SolnBlk.Nghost; 
    nr = SolnBlk.residual_variable;
    nn = SolnBlk.Number_of_Residual_Norms;
    if (SolnBlk.U != NULL) {
      block_allocated = 1;
    } else {
      block_allocated = 0;
    } /* endif */ 
  } /* endif */
  
  Communicator.Bcast(&ni, 1, MPI::INT, Source_Rank);
  Communicator.Bcast(&nj, 1, MPI::INT, Source_Rank); 
  Communicator.Bcast(&ng, 1, MPI::INT, Source_Rank);
  Communicator.Bcast(&nr, 1, MPI::INT, Source_Rank);
  Communicator.Bcast(&nn, 1, MPI::INT, Source_Rank);
  Communicator.Bcast(&block_allocated, 1, MPI::INT, Source_Rank);

  /* On non-source MPI processors, allocate (re-allocate) 
     memory for the quadrilateral solution block as necessary. */

  if (!(CFFC_MPI::This_Processor_Number == Source_CPU)) {
    if (SolnBlk.NCi != ni || SolnBlk.NCj != nj || SolnBlk.Nghost != ng) { 
      if (SolnBlk.U != NULL) SolnBlk.deallocate();
      if (block_allocated) SolnBlk.allocate(ni-2*ng, nj-2*ng, ng);
    } /* endif */
    //Set the block static variables if they were not previously assigned.
    if (SolnBlk.residual_variable != nr) SolnBlk.residual_variable = nr;
    if (SolnBlk.Number_of_Residual_Norms != nn) SolnBlk.Number_of_Residual_Norms = nn;
  } /* endif */

    /* Broadcast the axisymmetric/planar flow indicator. */

    Communicator.Bcast(&(SolnBlk.Flow_Type), 1, MPI::INT, Source_Rank);
    Communicator.Bcast(&(SolnBlk.Axisymmetric), 1, MPI::INT, Source_Rank);
    Communicator.Bcast(&(SolnBlk.Gravity), 1, MPI::INT, Source_Rank);
    Communicator.Bcast(&(SolnBlk.debug_level), 1, MPI::INT, Source_Rank);
    Communicator.Bcast(&(SolnBlk.Freeze_Limiter), 1, MPI::INT, Source_Rank);
     Communicator.Bcast(&(SolnBlk.Moving_wall_velocity), 1, MPI::DOUBLE, Source_Rank);
    Communicator.Bcast(&(SolnBlk.Pressure_gradient), 1, MPI::DOUBLE, Source_Rank);
 
    /* Broadcast the grid. */

    Broadcast_Quad_Block(SolnBlk.Grid, Communicator, Source_CPU);
    
    /* Broadcast the solution state variables. */

    if (block_allocated) {
       ni = (SolnBlk.ICu+SolnBlk.Nghost) - (SolnBlk.ICl-SolnBlk.Nghost) + 1;
       nj = (SolnBlk.JCu+SolnBlk.Nghost) - (SolnBlk.JCl-SolnBlk.Nghost) + 1;
       buffer = new double[NUM_VAR_FLAME2D*ni*nj];

       if (CFFC_MPI::This_Processor_Number == Source_CPU) {
	 buffer_size = 0;
	 for (int j  = SolnBlk.JCl-SolnBlk.Nghost ; j <= SolnBlk.JCu+SolnBlk.Nghost ; ++j ) {
	   for (int i = SolnBlk.ICl-SolnBlk.Nghost ; i <= SolnBlk.ICu+SolnBlk.Nghost ; ++i ) {
	     //Changed for Flame2D
	     for ( int k=0; k<NUM_VAR_FLAME2D; ++k) {	
	       buffer[buffer_size] = SolnBlk.W[i][j][k+1]; 
	       buffer_size = buffer_size + 1;
		}
	   } /* endfor */
	 } /* endfor */
       } /* endif */

       buffer_size = NUM_VAR_FLAME2D*ni*nj;
       Communicator.Bcast(buffer, buffer_size, MPI::DOUBLE, Source_Rank);

       if (!(CFFC_MPI::This_Processor_Number == Source_CPU)) {
	 buffer_size = 0;
          for (int j  = SolnBlk.JCl-SolnBlk.Nghost ; j <= SolnBlk.JCu+SolnBlk.Nghost ; ++j ) {
	    for (int i = SolnBlk.ICl-SolnBlk.Nghost ; i <= SolnBlk.ICu+SolnBlk.Nghost ; ++i ) {
	      for ( int k=0; k<NUM_VAR_FLAME2D; ++k) {	       
		SolnBlk.W[i][j][k+1]= buffer[buffer_size];
		buffer_size = buffer_size + 1;
	      }
	    } /* endfor */
	  } /* endfor */
       } /* endif */

       delete []buffer; 
       buffer = NULL;

       ni = 1;
       nj = (SolnBlk.JCu+SolnBlk.Nghost) - (SolnBlk.JCl-SolnBlk.Nghost) + 1;
       buffer = new double[2*NUM_VAR_FLAME2D*ni*nj];

       if (CFFC_MPI::This_Processor_Number == Source_CPU) {
          buffer_size = 0;
	  for (int j  = SolnBlk.JCl-SolnBlk.Nghost ; j <= SolnBlk.JCu+SolnBlk.Nghost ; ++j ) {
	    for ( int k=0; k<NUM_VAR_FLAME2D; ++k) {	
	      buffer[buffer_size]= SolnBlk.WoW[j][k+1];	
	      buffer_size = buffer_size + 1;
	    }
	    for ( int k=0; k<NUM_VAR_FLAME2D; ++k) {	
	      buffer[buffer_size]= SolnBlk.WoE[j][k+1];
	      buffer_size = buffer_size + 1;
	    }
          } /* endfor */
       } /* endif */

       buffer_size = 2*NUM_VAR_FLAME2D*ni*nj;
       Communicator.Bcast(buffer, buffer_size, MPI::DOUBLE, Source_Rank);

       if (!(CFFC_MPI::This_Processor_Number == Source_CPU)) {
          buffer_size = 0;
	  for (int j  = SolnBlk.JCl-SolnBlk.Nghost ; j <= SolnBlk.JCu+SolnBlk.Nghost ; ++j ) {    
	    for ( int k=0; k<NUM_VAR_FLAME2D; ++k) {	
	       SolnBlk.WoW[j][k+1] = buffer[buffer_size];
	       buffer_size = buffer_size + 1;
	    }
	    for ( int k=0; k<NUM_VAR_FLAME2D; ++k) {
	      SolnBlk.WoE[j][k+1] = buffer[buffer_size];
	      buffer_size = buffer_size + 1;
	    }
          } /* endfor */
       } /* endif */

       delete []buffer; 
       buffer = NULL;
       ni = (SolnBlk.ICu+SolnBlk.Nghost) - (SolnBlk.ICl-SolnBlk.Nghost) + 1;
       nj = 1;
       buffer = new double[2*NUM_VAR_FLAME2D*ni*nj];

       if (CFFC_MPI::This_Processor_Number == Source_CPU) {
          buffer_size = 0;  
          for (int i  = SolnBlk.ICl-SolnBlk.Nghost ; i <= SolnBlk.ICu+SolnBlk.Nghost ; ++i ) {
	    for ( int k=0; k<NUM_VAR_FLAME2D; ++k) {    
	      buffer[buffer_size]= SolnBlk.WoS[i][k+1];   
	      buffer_size = buffer_size + 1;
	    }
	    for ( int k=0; k<NUM_VAR_FLAME2D; ++k) {    
	      buffer[buffer_size]= SolnBlk.WoN[i][k+1];
	      buffer_size = buffer_size + 1;
	    } 
          } /* endfor */
       } /* endif */
       buffer_size = 2*NUM_VAR_FLAME2D*ni*nj;
       Communicator.Bcast(buffer, buffer_size, MPI::DOUBLE, Source_Rank);

       if (!(CFFC_MPI::This_Processor_Number == Source_CPU)) {
	 buffer_size = 0;
	  for (int i  = SolnBlk.ICl-SolnBlk.Nghost ; i <= SolnBlk.ICu+SolnBlk.Nghost ; ++i ) {
	    for ( int k=0; k<NUM_VAR_FLAME2D; ++k) {	
	      SolnBlk.WoS[i][k+1] = buffer[buffer_size]; 
	      buffer_size = buffer_size + 1;
	    }
	    for ( int k=0; k<NUM_VAR_FLAME2D; ++k) {	
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
void Copy_Solution_Block(Flame2D_Quad_Block &SolnBlk1,
                         Flame2D_Quad_Block &SolnBlk2) {

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
    SolnBlk1.Freeze_Limiter = SolnBlk2.Freeze_Limiter;
    SolnBlk1.Moving_wall_velocity = SolnBlk2.Moving_wall_velocity; 
    SolnBlk1.Pressure_gradient =  SolnBlk2.Pressure_gradient;

    /* Copy the grid of the second solution block
       to the first solution block. */

    Copy_Quad_Block(SolnBlk1.Grid, SolnBlk2.Grid);

    /* Copy the solution information from SolnBlk2 to SolnBlk1. */

    if (SolnBlk2.U != NULL) {
       for ( j  = SolnBlk1.JCl-SolnBlk1.Nghost ; j <= SolnBlk1.JCu+SolnBlk1.Nghost ; ++j ) {
          for ( i = SolnBlk1.ICl-SolnBlk1.Nghost ; i <= SolnBlk1.ICu+SolnBlk1.Nghost ; ++i ) {
             SolnBlk1.U[i][j] = SolnBlk2.U[i][j];
             SolnBlk1.W[i][j] = SolnBlk2.W[i][j];
             for ( k = 0 ; k <= NUMBER_OF_RESIDUAL_VECTORS_FLAME2D-1 ; ++k ) {
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
int Prolong_Solution_Block(Flame2D_Quad_Block &SolnBlk_Fine,
		            Flame2D_Quad_Block &SolnBlk_Original,
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
       SolnBlk_Fine.Freeze_Limiter = SolnBlk_Original.Freeze_Limiter;
       SolnBlk_Fine.Moving_wall_velocity = SolnBlk_Original.Moving_wall_velocity;        
       SolnBlk_Fine.Pressure_gradient =  SolnBlk_Original.Pressure_gradient;

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

    	       SolnBlk_Fine.W[2*(i-i_min)+SolnBlk_Fine.ICl  ]
                             [2*(j-j_min)+SolnBlk_Fine.JCl  ] 
                  = SolnBlk_Original.W[i][j];

    	       SolnBlk_Fine.W[2*(i-i_min)+SolnBlk_Fine.ICl+1]
                             [2*(j-j_min)+SolnBlk_Fine.JCl  ] 
                  = SolnBlk_Original.W[i][j];

    	       SolnBlk_Fine.W[2*(i-i_min)+SolnBlk_Fine.ICl  ]
                             [2*(j-j_min)+SolnBlk_Fine.JCl+1]
                  = SolnBlk_Original.W[i][j];

    	       SolnBlk_Fine.W[2*(i-i_min)+SolnBlk_Fine.ICl+1]
                             [2*(j-j_min)+SolnBlk_Fine.JCl+1]
                  = SolnBlk_Original.W[i][j];
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
int Restrict_Solution_Block(Flame2D_Quad_Block &SolnBlk_Coarse,
		             Flame2D_Quad_Block &SolnBlk_Original_SW,
                             Flame2D_Quad_Block &SolnBlk_Original_SE,
                             Flame2D_Quad_Block &SolnBlk_Original_NW,
                             Flame2D_Quad_Block &SolnBlk_Original_NE) {

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
       SolnBlk_Coarse.Freeze_Limiter = SolnBlk_Original_SW.Freeze_Limiter;
       SolnBlk_Coarse.Moving_wall_velocity = SolnBlk_Original_SW.Moving_wall_velocity;        
       SolnBlk_Coarse.Pressure_gradient =  SolnBlk_Original_SW.Pressure_gradient;

      /* Restrict the solution information from the four original solution blocks
         to the newly coarsened solution block. */

      // South-West corner fine block:	
      for ( j = SolnBlk_Original_SW.JCl; j <= SolnBlk_Original_SW.JCu ; j += 2 ) {
      	  for ( i = SolnBlk_Original_SW.ICl ; i <= SolnBlk_Original_SW.ICu ; i += 2 ) {
      	     i_coarse = (i-SolnBlk_Original_SW.ICl)/2+
                        SolnBlk_Coarse.ICl;
      	     j_coarse = (j-SolnBlk_Original_SW.JCl)/2+
                        SolnBlk_Coarse.JCl;
             SolnBlk_Coarse.W[i_coarse][j_coarse] = (SolnBlk_Original_SW.Grid.Cell[i  ][j  ].A*
                                                     SolnBlk_Original_SW.W[i  ][j  ] +
                                                     SolnBlk_Original_SW.Grid.Cell[i+1][j  ].A*
                                                     SolnBlk_Original_SW.W[i+1][j  ] + 
                                                     SolnBlk_Original_SW.Grid.Cell[i  ][j+1].A*
                                                     SolnBlk_Original_SW.W[i  ][j+1] +
                                                     SolnBlk_Original_SW.Grid.Cell[i+1][j+1].A*
                                                     SolnBlk_Original_SW.W[i+1][j+1]) /
                                                    (SolnBlk_Original_SW.Grid.Cell[i  ][j  ].A +
                                                     SolnBlk_Original_SW.Grid.Cell[i+1][j  ].A +
                                                     SolnBlk_Original_SW.Grid.Cell[i  ][j+1].A +
                                                     SolnBlk_Original_SW.Grid.Cell[i+1][j+1].A);
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
             SolnBlk_Coarse.W[i_coarse][j_coarse] = (SolnBlk_Original_SE.Grid.Cell[i  ][j  ].A*
                                                     SolnBlk_Original_SE.W[i  ][j  ] +
                                                     SolnBlk_Original_SE.Grid.Cell[i+1][j  ].A*
                                                     SolnBlk_Original_SE.W[i+1][j  ] + 
                                                     SolnBlk_Original_SE.Grid.Cell[i  ][j+1].A*
                                                     SolnBlk_Original_SE.W[i  ][j+1] +
                                                     SolnBlk_Original_SE.Grid.Cell[i+1][j+1].A*
                                                     SolnBlk_Original_SE.W[i+1][j+1]) /
                                                    (SolnBlk_Original_SE.Grid.Cell[i  ][j  ].A +
                                                     SolnBlk_Original_SE.Grid.Cell[i+1][j  ].A +
                                                     SolnBlk_Original_SE.Grid.Cell[i  ][j+1].A +
                                                     SolnBlk_Original_SE.Grid.Cell[i+1][j+1].A);
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
             SolnBlk_Coarse.W[i_coarse][j_coarse] = (SolnBlk_Original_NW.Grid.Cell[i  ][j  ].A*
                                                     SolnBlk_Original_NW.W[i  ][j  ] +
                                                     SolnBlk_Original_NW.Grid.Cell[i+1][j  ].A*
                                                     SolnBlk_Original_NW.W[i+1][j  ] + 
                                                     SolnBlk_Original_NW.Grid.Cell[i  ][j+1].A*
                                                     SolnBlk_Original_NW.W[i  ][j+1] +
                                                     SolnBlk_Original_NW.Grid.Cell[i+1][j+1].A*
                                                     SolnBlk_Original_NW.W[i+1][j+1]) /
                                                    (SolnBlk_Original_NW.Grid.Cell[i  ][j  ].A +
                                                     SolnBlk_Original_NW.Grid.Cell[i+1][j  ].A +
                                                     SolnBlk_Original_NW.Grid.Cell[i  ][j+1].A +
                                                     SolnBlk_Original_NW.Grid.Cell[i+1][j+1].A);
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
             SolnBlk_Coarse.W[i_coarse][j_coarse] = (SolnBlk_Original_NE.Grid.Cell[i  ][j  ].A*
                                                     SolnBlk_Original_NE.W[i  ][j  ] +
                                                     SolnBlk_Original_NE.Grid.Cell[i+1][j  ].A*
                                                     SolnBlk_Original_NE.W[i+1][j  ] + 
                                                     SolnBlk_Original_NE.Grid.Cell[i  ][j+1].A*
                                                     SolnBlk_Original_NE.W[i  ][j+1] +
                                                     SolnBlk_Original_NE.Grid.Cell[i+1][j+1].A*
                                                     SolnBlk_Original_NE.W[i+1][j+1]) /
                                                    (SolnBlk_Original_NE.Grid.Cell[i  ][j  ].A +
                                                     SolnBlk_Original_NE.Grid.Cell[i+1][j  ].A +
                                                     SolnBlk_Original_NE.Grid.Cell[i  ][j+1].A +
                                                     SolnBlk_Original_NE.Grid.Cell[i+1][j+1].A);
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
void Output_Tecplot(Flame2D_Quad_Block &SolnBlk,
		    Flame2D_Input_Parameters &IP,
                    const int Number_of_Time_Steps,
                    const double &Time,
                    const int Block_Number,
                    const int Output_Title,
	            ostream &Out_File) {
   
  /* Cell centered shear and qflux */
  if (SolnBlk.Flow_Type != FLOWTYPE_INVISCID) {
    // compute shear and qflux
    Viscous_Calculations(SolnBlk);
  }

  /* Ensure boundary conditions are updated before
     evaluating solution at the nodes. */
  BCs(SolnBlk,IP);

  /* Output node solution data. */
  Out_File << setprecision(14);
  if (Output_Title) {
 
    Out_File << "TITLE = \"" << CFFC_Name() << ": 2D Flame2D Solution, "
	     << "Time Step/Iteration Level = " << Number_of_Time_Steps
	     << ", Time = " << Time
	     << "\"" << "\n"
	     << "VARIABLES = \"x\" \\ \n"
	     << "\"y\" \\ \n"
	     << "\"rho\" \\ \n"
	     << "\"u\" \\ \n"
	     << "\"v\" \\ \n"
	     << "\"p\" \\ \n";
    //n species mass fractions names
    for(int i =0; i<Mixture::nSpecies(); i++){
      Out_File <<"\"c_"<<Mixture::speciesName(i)<<"\" \\ \n";
    }  
    //  Calculated values
    Out_File<< "\"T\" \\  \n" 
	    << "\"R\" \\  \n"
	    << "\"M\" \\  \n"
            << "\"viscosity\" \\  \n"
            << "\"thermal conduct\" \\  \n"
	    << "\"Prandtl\" \\  \n"
	    << "\"rho*H\"  \\ \n"  
	    <<"\"h\" \\ \n"
	    <<"\"h_s\" \\ \n"
	    <<"\"rho*E\" \\ \n"
	    << "\"e\" \\  \n" 
	    << "\"e_s\" \\ \n"
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
      const Flame2D_pState W_node( SolnBlk.Wn(i, j) );
      //coordinates, cell properties
      Out_File << " " << SolnBlk.Grid.Node[i][j].X<<endl<<W_node;
      //T,M,H,s, and all the rest of the calcuatated parameters
      Out_File.setf(ios::scientific);
      Out_File << " " << W_node.T()<< " " << W_node.Rtot()
	       << " " << W_node.vabs()/W_node.a() 
	       << " " << W_node.mu() <<" "<< W_node.kappa()
	       << " " << W_node.Pr()
	       << " " << W_node.H() <<" "<< W_node.h() <<" "<< W_node.hs()
	       << " " << W_node.E() <<" "<< W_node.e() <<" "<< W_node.es()<<endl;
  
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
void Output_Cells_Tecplot(Flame2D_Quad_Block &SolnBlk,
			  Flame2D_Input_Parameters &IP,
                          const int Number_of_Time_Steps,
                          const double &Time,
                          const int Block_Number,
                          const int Output_Title,
	                  ostream &Out_File) {

    int i, j;

    /* Cell centered shear and qflux */
    if (SolnBlk.Flow_Type != FLOWTYPE_INVISCID) {
      // compute shear and qflux
      Viscous_Calculations(SolnBlk);
    }

    BCs(SolnBlk,IP);

    Out_File << setprecision(14);
    if (Output_Title) {
       Out_File << "TITLE = \"" << CFFC_Name() << ": 2D Flame2D Solution, "
                << "Time Step/Iteration Level = " << Number_of_Time_Steps
                << ", Time = " << Time
                << "\"" << "\n"
   	        << "VARIABLES = \"x\" \\ \n"
                << "\"y\" \\ \n"
                << "\"rho\" \\ \n"
                << "\"u\" \\ \n"
                << "\"v\" \\ \n"
		<< "\"p\" \\ \n";
       //n species mass fractions names
       for(int i =0; i<Mixture::nSpecies(); i++){
	 Out_File <<"\"c"<<Mixture::speciesName(i)<<"\" \\ \n";
       }
        //Calculated values
       Out_File << "\"T\" \\ \n"
	        << "\"M\" \\ \n"
                << "\"R\" \\ \n"
	        << "\"viscosity\" \\ \n"
	        << "\"thermal conduct\" \\ \n"
	        << "\"Prandtl\" \\ \n";

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

    for ( j  = SolnBlk.JCl-SolnBlk.Nghost ; j <= SolnBlk.JCu + SolnBlk.Nghost ; ++j ) {
       for ( i = SolnBlk.ICl-SolnBlk.Nghost ; i <= SolnBlk.ICu + SolnBlk.Nghost ; ++i ) {
           Out_File << " "  << SolnBlk.Grid.Cell[i][j].Xc
                    << SolnBlk.W[i][j];
           Out_File.setf(ios::scientific);
	   //Temperature
	   Out_File << " " << SolnBlk.W[i][j].T()
		    << " " << SolnBlk.W[i][j].vabs()/SolnBlk.W[i][j].a() 
		    << " " << SolnBlk.W[i][j].Rtot()
		    << " " << SolnBlk.W[i][j].mu()
		    << " " << SolnBlk.W[i][j].kappa()
		    << " " << SolnBlk.W[i][j].Pr(); 
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
void Output_Nodes_Tecplot(Flame2D_Quad_Block &SolnBlk,
                          const int Number_of_Time_Steps,
                          const double &Time,
                          const int Block_Number,
                          const int Output_Title,
	                  ostream &Out_File) {

    int i, j;

    Out_File << setprecision(14);
    if (Output_Title) {
       Out_File << "TITLE = \"" << CFFC_Name() << ": 2D Flame2D Solution, "
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
 * Routine: ICs                                         *
 *                                                      *
 * Assigns initial conditions and data to the           *
 * solution variables of the specified quadrilateral    *
 * solution block.                                      *
 *                                                      *
 ********************************************************/
void ICs(Flame2D_Quad_Block &SolnBlk,
	 const int i_ICtype,
	 Flame2D_pState *Wo, Flame2D_Input_Parameters &Input_Parameters) {

    
}




/********************************************************
 * Routine: BCs                                         *
 *                                                      *
 * Apply boundary conditions at boundaries of the       *
 * specified quadrilateral solution block.              *
 *                                                      *
 ********************************************************/
void BCs(Flame2D_Quad_Block &SolnBlk, 
	 Flame2D_Input_Parameters &IP) {


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
double CFL(Flame2D_Quad_Block &SolnBlk,
           Flame2D_Input_Parameters &Input_Parameters) {

   
  double dtMin, dt_vis, dt_chem, d_i, d_j, v_i, v_j, a, rhomu, rhomut, delta_n;
  
  dtMin = MILLION;
  dt_vis = MILLION;
  dt_chem = MILLION;

  // Modifications for NKS overlap 
  int JCl_overlap = 0, JCu_overlap = 0, ICu_overlap = 0, ICl_overlap = 0;
  if(Input_Parameters.NKS_IP.GMRES_Overlap > 0 ){	
    if (SolnBlk.Grid.BCtypeS[SolnBlk.ICl] == BC_NONE)  JCl_overlap = Input_Parameters.NKS_IP.GMRES_Overlap;
    if (SolnBlk.Grid.BCtypeN[SolnBlk.ICu] == BC_NONE)  JCu_overlap = Input_Parameters.NKS_IP.GMRES_Overlap;
    if (SolnBlk.Grid.BCtypeE[SolnBlk.JCu] == BC_NONE)  ICu_overlap = Input_Parameters.NKS_IP.GMRES_Overlap;
    if (SolnBlk.Grid.BCtypeW[SolnBlk.JCl] == BC_NONE)  ICl_overlap = Input_Parameters.NKS_IP.GMRES_Overlap;       
  }
    
  for ( int j  = SolnBlk.JCl-SolnBlk.Nghost ; j <= SolnBlk.JCu+SolnBlk.Nghost ; ++j ) {
    for ( int i = SolnBlk.ICl-SolnBlk.Nghost ; i <= SolnBlk.ICu+SolnBlk.Nghost ; ++i ) {
      
      if (i < SolnBlk.ICl-ICl_overlap  || i > SolnBlk.ICu+ICu_overlap ||
	  j < SolnBlk.JCl-JCl_overlap  || j > SolnBlk.JCu+JCu_overlap) {
	SolnBlk.dt[i][j] = ZERO;
      } else {
	d_i = TWO*(SolnBlk.Grid.Cell[i][j].A/
		   (SolnBlk.Grid.lfaceE(i, j)+SolnBlk.Grid.lfaceW(i, j)));
	d_j = TWO*(SolnBlk.Grid.Cell[i][j].A/
		   (SolnBlk.Grid.lfaceN(i, j)+SolnBlk.Grid.lfaceS(i, j)));
	v_i = HALF*(SolnBlk.W[i][j].v()*
			 (SolnBlk.Grid.nfaceE(i, j)-SolnBlk.Grid.nfaceW(i, j)));
	v_j = HALF*(SolnBlk.W[i][j].v()*
		    (SolnBlk.Grid.nfaceN(i, j)-SolnBlk.Grid.nfaceS(i, j)));
	
	// set constant pointer,  this is needed to force the compiler to pick
	// the constant version of Flame2D_pState::rho()
	const Flame2D_pState &W = SolnBlk.W[i][j];

	/******** Inviscid deltat calculation **********************/   
	//no preconditioning
	
	if(Input_Parameters.Preconditioning == 0){
	  
	  SolnBlk.dt[i][j] = min(d_i/(W.a()+fabs(v_i)), d_j/(a+fabs(v_j)));
	  
	  //Low Mach Number Preconditioning
	} else if(Input_Parameters.Preconditioning == 1) { 
	  
	  delta_n = min(fabs(d_i),fabs(d_j));
	  SolnBlk.dt[i][j] = min(d_i/W.u_plus_aprecon(fabs(v_i),
						      SolnBlk.Flow_Type,
						      delta_n),
				 d_j/W.u_plus_aprecon(fabs(v_j),
						      SolnBlk.Flow_Type,
						      delta_n));	    
	}
	
	/******** Viscous deltat calculation ************************/   
	if (SolnBlk.Flow_Type != FLOWTYPE_INVISCID && 
	    Input_Parameters.Preconditioning == 0) {  
	  
	  rhomu = W.mu()/W.rho();
	  dt_vis = min((d_i*d_i)/(TWO*rhomu), (d_j*d_j)/(TWO*rhomu)); 
	  SolnBlk.dt[i][j]  = min(dt_vis, SolnBlk.dt[i][j]);
	}	  
	
	/******** Chemical Source Term deltat calculation ************/   
// 	dt_chem = HALF/SolnBlk.W[i][j].dSwdU_max_diagonal(Input_Parameters.Preconditioning,
// 							  delta_n,Input_Parameters.Solver_Type);
// 	dt_chem *= Input_Parameters.Source_Term_Multiplyer; // scale source term
// 	SolnBlk.dt[i][j] = min(dt_chem, SolnBlk.dt[i][j]);	  
	
	/************ Global Minimum ********************************/
	dtMin = min(dtMin, SolnBlk.dt[i][j]);
	
	//cout<<i<<" "<<j<<" Viscous "<<dt_vis<<" Chemistry "<<dt_chem<<" Final " <<dtMin<<endl;
	
      } 
    } 
  } 

  // THERE IS PROBABLY A FASTER WAY TO DO THIS, WITHOUT GOING THROUGH ALL CELLS -> CHANGE LATER
  // Set Ghost Cells dt to dtMin  
  for ( int j  = SolnBlk.JCl-SolnBlk.Nghost ; j <= SolnBlk.JCu+SolnBlk.Nghost ; ++j ) {
    for ( int i = SolnBlk.ICl-SolnBlk.Nghost ; i <= SolnBlk.ICu+SolnBlk.Nghost ; ++i ) {
      if (i < SolnBlk.ICl-ICl_overlap  || i > SolnBlk.ICu+ICu_overlap ||
	  j < SolnBlk.JCl-JCl_overlap  || j > SolnBlk.JCu+JCu_overlap) {
	SolnBlk.dt[i][j] = dtMin;
      } 
    } 
  }   
  
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
void Set_Global_TimeStep(Flame2D_Quad_Block &SolnBlk,
                         const double &Dt_min) {

    for (int j = SolnBlk.JCl ; j <= SolnBlk.JCu ; ++j ) {
       for (int i = SolnBlk.ICl ; i <= SolnBlk.ICu ; ++i ) {
          SolnBlk.dt[i][j] = Dt_min;
       } 
    } 

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
double L1_Norm_Residual(Flame2D_Quad_Block &SolnBlk, const int &norm) {
 
  double l1_norm(ZERO);
  for (int j  = SolnBlk.JCl ; j <= SolnBlk.JCu ;j++ ) {
    for (int i = SolnBlk.ICl ; i <= SolnBlk.ICu ;i++ ) {
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
double L2_Norm_Residual(Flame2D_Quad_Block &SolnBlk, const int &norm) {

  double l2_norm(ZERO);
  for (int j  = SolnBlk.JCl ; j <= SolnBlk.JCu ;j++ ) {
    for (int i = SolnBlk.ICl ; i <= SolnBlk.ICu ;i++ ) {
      l2_norm += sqr(SolnBlk.dUdt[i][j][0][norm]);
    } 
  }   
  return (sqrt(l2_norm));

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
double Max_Norm_Residual(Flame2D_Quad_Block &SolnBlk, const int &norm) {
 
  double max_norm(ZERO);
  for (int j  = SolnBlk.JCl ; j <= SolnBlk.JCu ;j++ ) {
    for (int i = SolnBlk.ICl ; i <= SolnBlk.ICu ;i++ ) {
      max_norm = max(max_norm,fabs(SolnBlk.dUdt[i][j][0][norm]));
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
void Linear_Reconstruction_GreenGauss(Flame2D_Quad_Block &SolnBlk,
				      const int i, 
                                      const int j,
                                      const int Limiter) {

    int n, n2, n_pts, i_index[8], j_index[8];
    double u0Min, u0Max, uQuad[4], phi;
    double DxDx_ave, DxDy_ave, DyDy_ave;
    double l_north, l_south, l_east, l_west;
    Vector2D n_north, n_south, n_east, n_west, dX;
    Flame2D_State W_nw, W_ne, W_sw, W_se, W_face;
    
    int NUM_VAR_FLAME2D = SolnBlk.NumVar();
  
   /* Carry out the limited solution reconstruction in
       the specified cell of the computational mesh. */
    
    // Determine the number of neighbouring cells to
    // be used in the reconstruction procedure.  Away from
    // boundaries this 8 neighbours will be used.
    
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

      } /* endif */
	   
	// Calculate slope limiters.    
      if (!SolnBlk.Freeze_Limiter) {
	for ( n = 1 ; n <= NUM_VAR_FLAME2D ; ++n ) {
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
void Linear_Reconstruction_GreenGauss(Flame2D_Quad_Block &SolnBlk,
				      const int Limiter) {

    /* Carry out the limited solution reconstruction in
       each cell of the computational mesh. */
    for (int j  = SolnBlk.JCl-SolnBlk.Nghost+1 ; j <= SolnBlk.JCu+SolnBlk.Nghost-1 ; ++j ) {
      for ( int i = SolnBlk.ICl-SolnBlk.Nghost+1 ; i <= SolnBlk.ICu+SolnBlk.Nghost-1 ; ++i ) {
	Linear_Reconstruction_GreenGauss(SolnBlk, i, j, Limiter);
      } /* endfor */
    } /* endfor */

}



void Linear_Reconstruction_LeastSquares(Flame2D_Quad_Block &SolnBlk,
				        const int i, 
                                        const int j,
                                        const int Limiter) {

    int n, n2, n_pts, i_index[8], j_index[8];
    double u0Min, u0Max, uQuad[4], phi;
    double DxDx_ave, DxDy_ave, DyDy_ave;
    Vector2D dX;
    Flame2D_State DU, DUDx_ave, DUDy_ave;
    
    //For the derivatives of gradients purpose
    // i_index_neigbor [8][1]  1 -- North , 8 -- its 8 sorrounding cells
    // i_index_neigbor [8][2]  2 -- South , 8 -- its 8 sorrounding cells
    // i_index_neigbor [8][3]  3 -- South , 8 -- its 8 sorrounding cells
    // i_index_neigbor [8][4]  4 -- South , 8 -- its 8 sorrounding cells

    int i_index_neigbor[8][5], j_index_neigbor[8][5];
    double dxdx_neigbor, dxdy_neigbor, dydy_neigbor;
    Vector2D dX_neigbor;
    Flame2D_State DW, DWDx_ave, DWDy_ave;

    int num=0;

    int NUM_VAR_FLAME2D = SolnBlk.NumVar();

    /* Carry out the limited solution reconstruction in
       each cell of the computational mesh. */
    
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
     
      for ( n2 = 0 ; n2 <= n_pts-1 ; ++n2 ) {
	dX = SolnBlk.Grid.Cell[ i_index[n2] ][ j_index[n2] ].Xc - 
	  SolnBlk.Grid.Cell[i][j].Xc;
	DU = SolnBlk.W[ i_index[n2] ][ j_index[n2] ] - SolnBlk.W[i][j];
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
	
     	if (!SolnBlk.Freeze_Limiter) {
	  for ( n = 1 ; n <= NUM_VAR_FLAME2D ; ++n ) {
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
void Linear_Reconstruction_LeastSquares(Flame2D_Quad_Block &SolnBlk,
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
void Linear_Reconstruction_LeastSquares_2(Flame2D_Quad_Block &SolnBlk,
				          const int i, 
                                          const int j,
                                          const int Limiter) {

    int n, n2, n_pts, i_index[8], j_index[8];
    double u0Min, u0Max, uQuad[4], phi;
    double DxDx_ave, DxDy_ave, DyDy_ave;
    Vector2D dX;
    Flame2D_State DU, DUDx_ave, DUDy_ave;
   
    int num=0;
    int NUM_VAR_FLAME2D = SolnBlk.NumVar();

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
		 SolnBlk.Grid.BCtypeW[j] == BC_1DFLAME_OUTFLOW) {
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
		 SolnBlk.Grid.BCtypeE[j] == BC_1DFLAME_OUTFLOW) {
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
		 SolnBlk.Grid.BCtypeS[i] == BC_1DFLAME_OUTFLOW) {
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
		 SolnBlk.Grid.BCtypeN[i] == BC_1DFLAME_OUTFLOW) {
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
	  for ( n = 1 ; n <= NUM_VAR_FLAME2D ; ++n ) {
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
void Linear_Reconstruction_LeastSquares_Diamond(Flame2D_Quad_Block &SolnBlk,
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
void Linear_Reconstruction_LeastSquares_Diamond(Flame2D_Quad_Block &SolnBlk,
						const int i, 
						const int j,
						const int Limiter) {
  int n, n2, n_pts, i_index[8], j_index[8];
  int n_neigbour;
  double u0Min, u0Max, uQuad[4], phi;
  double DxDx_ave, DxDy_ave, DyDy_ave;
  double area[4];

  Vector2D deltaX, dX[8];  
  Flame2D_State DU[8], DUDx_ave, DUDy_ave;
  Flame2D_State QuadraturePoint_N, QuadraturePoint_E, 
    QuadraturePoint_S, QuadraturePoint_W;  
  Flame2D_State TopVertex, BottomVertex; 
  
  Flame2D_State Temp;
  //For the derivatives of gradients purpose
  // i_index_neigbor [8][1]  1 -- North , 8 -- its 8 sorrounding cells
  // i_index_neigbor [8][2]  2 -- South , 8 -- its 8 sorrounding cells
  // i_index_neigbor [8][3]  3 -- South , 8 -- its 8 sorrounding cells
  // i_index_neigbor [8][4]  4 -- South , 8 -- its 8 sorrounding cells
  
  int i_index_neigbor [8][5], j_index_neigbor[8][5];
  double dxdx_neigbor, dxdy_neigbor, dydy_neigbor;
  Vector2D dX_neigbor;
  Flame2D_State DW, DWDx_ave, DWDy_ave;
  Flame2D_State W_node;
  
  int num=0;
  int NUM_VAR_FLAME2D = SolnBlk.NumVar();
  
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
   
  n_neigbour = 4;
  QuadraturePoint_N.Vacuum();
  QuadraturePoint_E.Vacuum();
  QuadraturePoint_S.Vacuum();
  QuadraturePoint_W.Vacuum(); 
//   // Derivatives of Gradients ...
//   double dWnNWdWc = ZERO;
//   double dWnNEdWc = ZERO;
//   double dWnSWdWc = ZERO;
//   double dWnSEdWc = ZERO;
//   double dQPdWc = ZERO;
  

//   dWnNWdWc = SolnBlk.dWn_dWc(i, j+1, NW);    //THEST ARE REQUIRED FOR JACOBIANS ONLY, NOT REQUIRED FOR
//   dWnNEdWc = SolnBlk.dWn_dWc(i+1, j+1, NE);    //dWdx and dWdy , however they appear to be pretty cheap
//   dWnSWdWc = SolnBlk.dWn_dWc(i, j, SW);    //affecting performance <0.1%  - SHOULD RECHECK 
//   dWnSEdWc = SolnBlk.dWn_dWc(i+1, j, SE);
     
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
//     dQPdWc = HALF*(dWnNWdWc + dWnNEdWc); // derivative of quadraturepoint w.r.t Wc
    
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
    
    
    
    for ( n2 = 0 ; n2 <= n_neigbour -1 ; ++n2 ) {
      DUDx_ave += DU[n2]*dX[n2].x;
      DUDy_ave += DU[n2]*dX[n2].y;
      DxDx_ave += dX[n2].x*dX[n2].x;
      DxDy_ave += dX[n2].x*dX[n2].y;
      DyDy_ave += dX[n2].y*dX[n2].y;
      
    } /* endfor */
    
    DUDx_ave = DUDx_ave/double(n_neigbour);
    DUDy_ave = DUDy_ave/double(n_neigbour);
    DxDx_ave = DxDx_ave/double(n_neigbour);
    DxDy_ave = DxDy_ave/double(n_neigbour);
    DyDy_ave = DyDy_ave/double(n_neigbour);
    
    SolnBlk.dWdx_faceN[i][j] = (DUDx_ave*DyDy_ave-DUDy_ave*DxDy_ave)/
      (DxDx_ave*DyDy_ave-DxDy_ave*DxDy_ave);
    SolnBlk.dWdy_faceN[i][j] = (DUDy_ave*DxDx_ave-DUDx_ave*DxDy_ave)/
      (DxDx_ave*DyDy_ave-DxDy_ave*DxDy_ave);

    
//     SolnBlk.d_dWdx_dW[i][j][NORTH] =((dX[0].x*(ONE-dQPdWc) - dX[1].x*dQPdWc+ dX[2].x*(dWnNWdWc -dQPdWc ) 
//                                                            + dX[3].x*(dWnNEdWc-dQPdWc) )/double( n_neigbour)*DyDy_ave 
// 				                              - (dX[0].y*(ONE-dQPdWc)- dX[1].y*dQPdWc+dX[2].y*(dWnNWdWc -dQPdWc) 
//                                                            + dX[3].y*(dWnNEdWc-dQPdWc) )/double(n_neigbour)*DxDy_ave)
//                                                             /(DxDx_ave*DyDy_ave-DxDy_ave*DxDy_ave) ;
        
//     SolnBlk.d_dWdy_dW[i][j][NORTH] =((dX[0].y*(ONE-dQPdWc)-dX[1].y*dQPdWc + dX[2].y*(dWnNWdWc -dQPdWc) 
//                                                                           + dX[3].y*(dWnNEdWc-dQPdWc) )/double(n_neigbour)*DxDx_ave
//                                              			     - (dX[0].x*(ONE-dQPdWc)-dX[1].x*dQPdWc+ dX[2].x*(dWnNWdWc -dQPdWc) 
//                                                                           + dX[3].x*(dWnNEdWc-dQPdWc) )/double(n_neigbour) *DxDy_ave)
//                                                                             /(DxDx_ave*DyDy_ave-DxDy_ave*DxDy_ave);
    
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
//     dQPdWc = HALF*(dWnNEdWc+dWnSEdWc);
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
    
    for ( n2 = 0 ; n2 <=  n_neigbour-1 ; ++n2 ) {
      DUDx_ave += DU[n2]*dX[n2].x;
      DUDy_ave += DU[n2]*dX[n2].y;
      DxDx_ave += dX[n2].x*dX[n2].x;
      DxDy_ave += dX[n2].x*dX[n2].y;
      DyDy_ave += dX[n2].y*dX[n2].y;
      
    } /* endfor */
    
    DUDx_ave = DUDx_ave/double(n_neigbour);
    DUDy_ave = DUDy_ave/double(n_neigbour);
    DxDx_ave = DxDx_ave/double(n_neigbour);
    DxDy_ave = DxDy_ave/double(n_neigbour);
    DyDy_ave = DyDy_ave/double(n_neigbour);
    
    SolnBlk.dWdx_faceE[i][j] = (DUDx_ave*DyDy_ave-DUDy_ave*DxDy_ave)/
      (DxDx_ave*DyDy_ave-DxDy_ave*DxDy_ave);
    SolnBlk.dWdy_faceE[i][j] = (DUDy_ave*DxDx_ave-DUDx_ave*DxDy_ave)/
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
//     dQPdWc = HALF*(dWnSWdWc + dWnNWdWc);
    
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
    
    for ( n2 = 0 ; n2 <=  n_neigbour-1 ; ++n2 ) {
      DUDx_ave += DU[n2]*dX[n2].x;
      DUDy_ave += DU[n2]*dX[n2].y;
      DxDx_ave += dX[n2].x*dX[n2].x;
      DxDy_ave += dX[n2].x*dX[n2].y;
      DyDy_ave += dX[n2].y*dX[n2].y;
      
    } /* endfor */
    
    DUDx_ave = DUDx_ave/double(n_neigbour);
    DUDy_ave = DUDy_ave/double(n_neigbour);
    DxDx_ave = DxDx_ave/double(n_neigbour);
    DxDy_ave = DxDy_ave/double(n_neigbour);
    DyDy_ave = DyDy_ave/double(n_neigbour);
    
    SolnBlk.dWdx_faceW[i][j] = (DUDx_ave*DyDy_ave-DUDy_ave*DxDy_ave)/
      (DxDx_ave*DyDy_ave-DxDy_ave*DxDy_ave);
    SolnBlk.dWdy_faceW[i][j] = (DUDy_ave*DxDx_ave-DUDx_ave*DxDy_ave)/
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
    QuadraturePoint_S = HALF*(TopVertex + BottomVertex);
 //    dQPdWc = HALF*(dWnSEdWc + dWnSWdWc);
    
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
    
    for ( n2 = 0 ; n2 <=  n_neigbour-1 ; ++n2 ) {
      DUDx_ave += DU[n2]*dX[n2].x;
      DUDy_ave += DU[n2]*dX[n2].y;
      DxDx_ave += dX[n2].x*dX[n2].x;
      DxDy_ave += dX[n2].x*dX[n2].y;
      DyDy_ave += dX[n2].y*dX[n2].y;
      
    } /* endfor */
    
    DUDx_ave = DUDx_ave/double(n_neigbour);
    DUDy_ave = DUDy_ave/double(n_neigbour);
    DxDx_ave = DxDx_ave/double(n_neigbour);
    DxDy_ave = DxDy_ave/double(n_neigbour);
    DyDy_ave = DyDy_ave/double(n_neigbour);
    
    SolnBlk.dWdx_faceS[i][j] = (DUDx_ave*DyDy_ave-DUDy_ave*DxDy_ave)/
      (DxDx_ave*DyDy_ave-DxDy_ave*DxDy_ave);
    SolnBlk.dWdy_faceS[i][j] = (DUDy_ave*DxDx_ave-DUDx_ave*DxDy_ave)/
      (DxDx_ave*DyDy_ave-DxDy_ave*DxDy_ave);      

    //Area weighted gradients at cell centers
    SolnBlk.dWdx[i][j] = SolnBlk.dWdx_faceN[i][j]*area[0]+SolnBlk.dWdx_faceE[i][j]*area[1] +
      SolnBlk.dWdx_faceW[i][j]*area[2] + SolnBlk.dWdx_faceS[i][j]*area[3]; 
    
    SolnBlk.dWdx[i][j] = SolnBlk.dWdx[i][j]/SolnBlk.Grid.Cell[i][j].A; 
    
    SolnBlk.dWdy[i][j] = SolnBlk.dWdy_faceN[i][j]*area[0]+ SolnBlk.dWdy_faceE[i][j]*area[1] +
      SolnBlk.dWdy_faceW[i][j]*area[2] + SolnBlk.dWdy_faceS[i][j]*area[3]; 
    
    SolnBlk.dWdy[i][j] = SolnBlk.dWdy[i][j]/SolnBlk.Grid.Cell[i][j].A;
      
    
    if (!SolnBlk.Freeze_Limiter) {
      for ( n = 1 ; n <= NUM_VAR_FLAME2D ; ++n ) {
	u0Min = SolnBlk.W[i][j][n];
	u0Max = u0Min;
	for(n2 = 0 ; n2 <= n_pts-1 ; ++n2){
           u0Min = min(u0Min, SolnBlk.W[ i_index[n2] ][ j_index[n2] ][n]);
           u0Max = max(u0Max, SolnBlk.W[ i_index[n2] ][ j_index[n2] ][n]); 
	  
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
  } /* endif */
  
  
    
}

void Linear_Reconstruction_GreenGauss_Diamond(Flame2D_Quad_Block &SolnBlk,
					      const int Limiter) { 
    /* Carry out the limited solution reconstruction in
       each cell of the computational mesh. */
    for (int j  = SolnBlk.JCl-SolnBlk.Nghost+1 ; j <= SolnBlk.JCu+SolnBlk.Nghost-1 ; ++j ) {
       for (int i = SolnBlk.ICl-SolnBlk.Nghost+1 ; i <= SolnBlk.ICu+SolnBlk.Nghost-1 ; ++i ) {
	   Linear_Reconstruction_GreenGauss_Diamond(SolnBlk, i, j, Limiter);	
       } 
    } 

}

void Linear_Reconstruction_GreenGauss_Diamond(Flame2D_Quad_Block &SolnBlk,
                                              const int i, 
                                              const int j,
                                              const int Limiter) {

  int n_pts, i_index[8], j_index[8];
  double area[4], AREA;
  Flame2D_State W_SW, W_SE, W_NW, W_NE, W_average[4];
  Vector2D norm[4]; 
     
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
    
    W_NE = SolnBlk.Wn(i+1,j+1);
    W_NW = SolnBlk.Wn(i,j+1);
    W_SW = SolnBlk.Wn(i,j);
    W_SE = SolnBlk.Wn(i+1,j);
    
    /*************** NORTH ****************************/
    
    /*Formulate the gradients of primitive parameters on the north face of cell (i, j)*/
    // counterclockwise, starting from the cell center (i,j), nodeNE (i, j), 
    // top cell center (i, j+1), NodeNW (i, j);
    W_average[0] = HALF*(SolnBlk.W[i][j]   + W_NE);
    W_average[1] = HALF*(SolnBlk.W[i][j+1] + W_NE);
    W_average[2] = HALF*(SolnBlk.W[i][j+1] + W_NW);
    W_average[3] = HALF*(SolnBlk.W[i][j]   + W_NW);
    
    //  normal vector of the SE side of a diamond 
    norm[0].x = SolnBlk.Grid.nodeNE(i,j).X.y - SolnBlk.Grid.Cell[i][j].Xc.y;
    norm[0].y = SolnBlk.Grid.Cell[i][j].Xc.x - SolnBlk.Grid.nodeNE(i,j).X.x;
    //  normal vector of the NE side of a diamond 
    norm[1].x = SolnBlk.Grid.Cell[i][j+1].Xc.y - SolnBlk.Grid.nodeNE(i,j).X.y;
    norm[1].y = SolnBlk.Grid.nodeNE(i,j).X.x - SolnBlk.Grid.Cell[i][j+1].Xc.x;
    //  normal vector of the NW side of a diamond 
    norm[2].x = SolnBlk.Grid.nodeNW(i,j).X.y - SolnBlk.Grid.Cell[i][j+1].Xc.y;
    norm[2].y = SolnBlk.Grid.Cell[i][j+1].Xc.x - SolnBlk.Grid.nodeNW(i,j).X.x;
    //  normal vector of the SW side of a diamond 
    norm[3].x = SolnBlk.Grid.Cell[i][j].Xc.y - SolnBlk.Grid.nodeNW(i,j).X.y ;
    norm[3].y = SolnBlk.Grid.nodeNW(i,j).X.x - SolnBlk.Grid.Cell[i][j].Xc.x ;
    
    AREA =  HALF*(fabs((SolnBlk.Grid.nodeNE(i,j).X-SolnBlk.Grid.Cell[i][j].Xc)^
		       (SolnBlk.Grid.nodeNW(i,j).X-SolnBlk.Grid.Cell[i][j].Xc)) +
		  fabs((SolnBlk.Grid.nodeNW(i,j).X-SolnBlk.Grid.Cell[i][j+1].Xc)^
		       (SolnBlk.Grid.nodeNE(i,j).X-SolnBlk.Grid.Cell[i][j+1].Xc)));
    
    SolnBlk.dWdx_faceN[i][j] = ( W_average[0]* norm[0].x +W_average[1]* norm[1].x+ W_average[2]* norm[2].x+W_average[3]* norm[3].x)/AREA;
    SolnBlk.dWdy_faceN[i][j] = ( W_average[0]* norm[0].y +W_average[1]* norm[1].y+ W_average[2]* norm[2].y+W_average[3]* norm[3].y)/AREA;  
    
    /*************** EAST ****************************/
    
    /*Formulate the gradients of primitive parameters on the east face of cell (i, j)*/
    // counterclockwise, starting from  nodeSE(i,j), cell (i+1, j), 
    // nodeNE(i,j), cell center (i, j+1)     
    W_average[2] = W_average[0];
    W_average[0] = HALF*(SolnBlk.W[i+1][j]  + W_SE);      
    W_average[1] = HALF*(SolnBlk.W[i+1][j]  + W_NE);     
    W_average[3] = HALF*(SolnBlk.W[i][j]    + W_SE);
    
    //  normal vector of the NW side of a diamond  = - SE of previous
    norm[2] = - norm[0];
    //  normal vector of the SE side of a diamond 
    norm[0].x = SolnBlk.Grid.Cell[i+1][j].Xc.y - SolnBlk.Grid.nodeSE(i,j).X.y;
    norm[0].y = SolnBlk.Grid.nodeSE(i,j).X.x - SolnBlk.Grid.Cell[i+1][j].Xc.x;
    //  normal vector of the NE side of a diamond 
    norm[1].x = SolnBlk.Grid.nodeNE(i,j).X.y -  SolnBlk.Grid.Cell[i+1][j].Xc.y ;
    norm[1].y = SolnBlk.Grid.Cell[i+1][j].Xc.x - SolnBlk.Grid.nodeNE(i,j).X.x;
    //  normal vector of the SW side of a diamond 
    norm[3].x = SolnBlk.Grid.nodeSE(i,j).X.y - SolnBlk.Grid.Cell[i][j].Xc.y;
    norm[3].y = SolnBlk.Grid.Cell[i][j].Xc.x - SolnBlk.Grid.nodeSE(i,j).X.x;
    
    AREA =  HALF*(fabs((SolnBlk.Grid.nodeNE(i,j).X-SolnBlk.Grid.Cell[i+1][j].Xc)^
		       (SolnBlk.Grid.nodeSE(i,j).X-SolnBlk.Grid.Cell[i+1][j].Xc)) +
		  fabs((SolnBlk.Grid.nodeSE(i,j).X-SolnBlk.Grid.Cell[i][j].Xc)^
		       (SolnBlk.Grid.nodeNE(i,j).X-SolnBlk.Grid.Cell[i][j].Xc)));
    
    SolnBlk.dWdx_faceE[i][j] = ( W_average[0]* norm[0].x +W_average[1]* norm[1].x+ W_average[2]* norm[2].x+W_average[3]* norm[3].x)/AREA;
    SolnBlk.dWdy_faceE[i][j] = ( W_average[0]* norm[0].y +W_average[1]* norm[1].y+ W_average[2]* norm[2].y+W_average[3]* norm[3].y)/AREA;  
    
    /*************** SOUTH ****************************/
    
    /*Formulate the gradients of primitive parameters on the south face of cell (i, j)*/
    // counterclockwise, starting from  cell (i, j-1), nodeSE(i,j) 
    // cell(i,j), nodeSW(i,j)     
    W_average[1] = W_average[3];
    W_average[0] = HALF*(W_SE + SolnBlk.W[i][j-1]);
    W_average[2] = HALF*(SolnBlk.W[i][j]  + W_SW);
    W_average[3] = HALF*(SolnBlk.W[i][j-1] + W_SW);
    
    //  normal vector of the NE side of a diamond = -SW of previous
    norm[1] = -norm[3];
    //  normal vector of the SE side of a diamond 
    norm[0].x = SolnBlk.Grid.nodeSE(i,j).X.y - SolnBlk.Grid.Cell[i][j-1].Xc.y;
    norm[0].y = SolnBlk.Grid.Cell[i][j-1].Xc.x - SolnBlk.Grid.nodeSE(i,j).X.x;
    //  normal vector of the NW side of a diamond 
    norm[2].x = SolnBlk.Grid.nodeSW(i,j).X.y - SolnBlk.Grid.Cell[i][j].Xc.y ;
    norm[2].y = SolnBlk.Grid.Cell[i][j].Xc.x - SolnBlk.Grid.nodeSW(i,j).X.x;
    //  normal vector of the SW side of a diamond 
    norm[3].x = SolnBlk.Grid.Cell[i][j-1].Xc.y - SolnBlk.Grid.nodeSW(i,j).X.y;
    norm[3].y = SolnBlk.Grid.nodeSW(i,j).X.x - SolnBlk.Grid.Cell[i][j-1].Xc.x;
    
    AREA =  HALF*(fabs((SolnBlk.Grid.nodeSE(i,j).X-SolnBlk.Grid.Cell[i][j-1].Xc)^
		       (SolnBlk.Grid.nodeSW(i,j).X-SolnBlk.Grid.Cell[i][j-1].Xc)) +
		  fabs((SolnBlk.Grid.nodeSE(i,j).X-SolnBlk.Grid.Cell[i][j].Xc)^
		       (SolnBlk.Grid.nodeSW(i,j).X-SolnBlk.Grid.Cell[i][j].Xc)));
    
    SolnBlk.dWdx_faceS[i][j] = ( W_average[0]* norm[0].x +W_average[1]* norm[1].x+ W_average[2]* norm[2].x+W_average[3]* norm[3].x)/AREA;
    SolnBlk.dWdy_faceS[i][j] = ( W_average[0]* norm[0].y +W_average[1]* norm[1].y+ W_average[2]* norm[2].y+W_average[3]* norm[3].y)/AREA;  
    
    /*************** WEST ****************************/    
    /*Formulate the gradients of primitive parameters on the west face of cell (i, j ) */
    // counterclockwise, starting from  NodeSW(i,j) 
    // cell(i,j), nodeNW(i,j), cell (i-1, j)     
    W_average[0] = W_average[2];
    W_average[1] = HALF*(SolnBlk.W[i][j]   + W_NW);
    W_average[2] = HALF*(SolnBlk.W[i-1][j] + W_NW);
    W_average[3] = HALF*(SolnBlk.W[i-1][j] + W_SW);
    
    //  normal vector of the SE side of a diamond = - NW of previous
    norm[0] = - norm[2];
    //  normal vector of the NE side of a diamond 
    norm[1].x =  SolnBlk.Grid.nodeNW(i,j).X.y - SolnBlk.Grid.Cell[i][j].Xc.y;
    norm[1].y =  SolnBlk.Grid.Cell[i][j].Xc.x - SolnBlk.Grid.nodeNW(i,j).X.x;
    //  normal vector of the NW side of a diamond 
    norm[2].x =  SolnBlk.Grid.Cell[i-1][j].Xc.y - SolnBlk.Grid.nodeNW(i,j).X.y ;
    norm[2].y =  SolnBlk.Grid.nodeNW(i,j).X.x - SolnBlk.Grid.Cell[i-1][j].Xc.x;
    //  normal vector of the SW side of a diamond 
    norm[3].x =  SolnBlk.Grid.nodeSW(i,j).X.y - SolnBlk.Grid.Cell[i-1][j].Xc.y;
    norm[3].y =  SolnBlk.Grid.Cell[i-1][j].Xc.x - SolnBlk.Grid.nodeSW(i,j).X.x;
    
    AREA =  HALF*(fabs((SolnBlk.Grid.nodeNW(i,j).X-SolnBlk.Grid.Cell[i][j].Xc)^
		       (SolnBlk.Grid.nodeSW(i,j).X-SolnBlk.Grid.Cell[i][j].Xc)) +
		  fabs((SolnBlk.Grid.nodeNW(i,j).X-SolnBlk.Grid.Cell[i-1][j].Xc)^
		       (SolnBlk.Grid.nodeSW(i,j).X-SolnBlk.Grid.Cell[i-1][j].Xc)));
        
    SolnBlk.dWdx_faceW[i][j] = ( W_average[0]* norm[0].x +W_average[1]* norm[1].x+ W_average[2]* norm[2].x+W_average[3]* norm[3].x)/AREA;
    SolnBlk.dWdy_faceW[i][j] = ( W_average[0]* norm[0].y +W_average[1]* norm[1].y+ W_average[2]* norm[2].y+W_average[3]* norm[3].y)/AREA;  
    
    /*************** CENTER ****************************/
    
    // area weighted gradients at cell centers, 4 inside triangles
    area[0] = HALF*(SolnBlk.Grid.Node[i+1][j+1].X - SolnBlk.Grid.Node[i][j+1].X )^
      (SolnBlk.Grid.xfaceN(i,j)- SolnBlk.Grid.Cell[i][j].Xc);
    area[1] = HALF*(SolnBlk.Grid.Node[i+1][j+1].X - SolnBlk.Grid.Node[i+1][j].X )^
      (SolnBlk.Grid.Cell[i][j].Xc - SolnBlk.Grid.xfaceE(i,j));
    area[2] = HALF*(SolnBlk.Grid.Node[i][j+1].X - SolnBlk.Grid.Node[i][j].X )^
      ( SolnBlk.Grid.xfaceW(i, j) - SolnBlk.Grid.Cell[i][j].Xc );
    area[3] = HALF*(SolnBlk.Grid.Node[i+1][j].X - SolnBlk.Grid.Node[i][j].X )^
      (SolnBlk.Grid.Cell[i][j].Xc - SolnBlk.Grid.xfaceS(i,j) );
        
    //Reconstructed cell center gradients
    SolnBlk.dWdx[i][j] = SolnBlk.dWdx_faceN[i][j]*area[0]+SolnBlk.dWdx_faceE[i][j]*area[1] +
      SolnBlk.dWdx_faceW[i][j]*area[2] + SolnBlk.dWdx_faceS[i][j]*area[3];       
    SolnBlk.dWdx[i][j] = SolnBlk.dWdx[i][j]/SolnBlk.Grid.Cell[i][j].A; 
    
    SolnBlk.dWdy[i][j] = SolnBlk.dWdy_faceN[i][j]*area[0]+ SolnBlk.dWdy_faceE[i][j]*area[1] +
      SolnBlk.dWdy_faceW[i][j]*area[2] + SolnBlk.dWdy_faceS[i][j]*area[3];     
    SolnBlk.dWdy[i][j] = SolnBlk.dWdy[i][j]/SolnBlk.Grid.Cell[i][j].A; 
    
    /****************************************************/
        
    double u0Min, u0Max, uQuad[4], phi;
    double DxDx_ave, DxDy_ave, DyDy_ave;
    Vector2D deltaX; 
    
    // Calculate slope limiters.  
    if (!SolnBlk.Freeze_Limiter) {
      for ( int n=1 ; n <= SolnBlk.NumVar(); n++ ) {
	u0Min = SolnBlk.W[i][j][n];
	u0Max = u0Min;
	for(int n2 = 0 ; n2 <= n_pts-1 ; ++n2){
	  u0Min = min(u0Min, SolnBlk.W[ i_index[n2] ][ j_index[n2] ][n]);
	  u0Max = max(u0Max, SolnBlk.W[ i_index[n2] ][ j_index[n2] ][n]); 
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
    SolnBlk.dWdx_faceN[i][j].Vacuum();
    SolnBlk.dWdx_faceS[i][j].Vacuum();
    SolnBlk.dWdx_faceE[i][j].Vacuum();
    SolnBlk.dWdx_faceW[i][j].Vacuum();
    SolnBlk.dWdy_faceN[i][j].Vacuum();
    SolnBlk.dWdy_faceS[i][j].Vacuum();
    SolnBlk.dWdy_faceE[i][j].Vacuum();
    SolnBlk.dWdy_faceW[i][j].Vacuum();
    SolnBlk.dWdx[i][j].Vacuum();
    SolnBlk.dWdy[i][j].Vacuum();
    SolnBlk.phi[i][j].Vacuum(); 
  } 
      
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
void Residual_Smoothing(Flame2D_Quad_Block &SolnBlk,
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
				   Flame2D_Input_Parameters &IP, 
                                   int &number_refinement_criteria,				 
                                   Flame2D_Quad_Block &SolnBlk) {

    int i, j;

    double grad_rho_x, grad_rho_y, grad_rho_abs, grad_rho_criteria, grad_rho_criteria_max,
      div_V, div_V_criteria, div_V_criteria_max,
      curl_V_z, curl_V_abs, curl_V_criteria, curl_V_criteria_max,
      grad_Temp_x, grad_Temp_y, grad_Temp_abs, grad_Temp_criteria, grad_Temp_criteria_max,
      grad_CH4_x, grad_CH4_y, grad_CH4_abs, grad_CH4_criteria, grad_CH4_criteria_max,
      grad_CO2_x, grad_CO2_y, grad_CO2_abs, grad_CO2_criteria, grad_CO2_criteria_max;
//       grad_dudy, grad_dudy_criteria, grad_dudy_criteria_max,
//       grad_pressure_x, grad_pressure_y, grad_pressure_abs, grad_pressure_criteria, grad_pressure_criteria_max;

    /* Set the number of refinement criteria to be used (3):
       (1) refinement criteria #1 based on the gradient of the density field;
       (2) refinement criteria #2 based on the divergence of the velocity vector;
       (3) refinement criteria #3 based on the curl of the velocity vector. 
       (4) refinement criteria #4 based on the gradient of Temperature
       (5) refinement criteria #5 based on the gradient of CH4 mass fraction
       (6) refinement criteria #6 based on the gradient of CO2 mass fraction
       (7) refinement criteria #7 based on du/dy (for pipe flow)
       (8) refinement criteria #8 based on the gradient of Pressure
    */

    number_refinement_criteria = IP.Number_of_Refinement_Criteria;
    /* Initialize the refinement criteria for the block. */

    grad_rho_criteria_max = ZERO;
    div_V_criteria_max = ZERO;
    curl_V_criteria_max = ZERO;
    grad_Temp_criteria_max = ZERO;
    grad_CH4_criteria_max = ZERO;      
    grad_CO2_criteria_max = ZERO;
    //grad_pressure_criteria_max = ZERO;

    /* Calculate the refinement criteria for each cell of the 
       computational mesh and assign the maximum value for
       all cells as the refinement criteria for the solution 
       block. */

    for ( j  = SolnBlk.JCl-1 ; j <= SolnBlk.JCu+1 ; ++j ) {
       for ( i = SolnBlk.ICl-1 ; i <= SolnBlk.ICu+1 ; ++i ) {
	 // Reconstruct the solution within the cell.
	 //Linear_Reconstruction_GreenGauss(SolnBlk, i, j, LIMITER_UNLIMITED);
	 //Linear_Reconstruction_LeastSquares(SolnBlk, i, j, LIMITER_UNLIMITED);
	 //Required if initial refinement, ie. 0 iterations 
	 Linear_Reconstruction_LeastSquares_2(SolnBlk, i, j, LIMITER_UNLIMITED);

          if (SolnBlk.Grid.Cell[i][j].A > ZERO) {

	    // set constant pointer,  this is needed to force the compiler to pick
	    // the constant version of Flame2D_pState::rho()
	    const Flame2D_pState &W = SolnBlk.W[i][j];

             // Evaluate refinement criteria #1 based on the gradient
             // of the density field.
             grad_rho_x = SolnBlk.dWdx[i][j].rho();
             grad_rho_y = SolnBlk.dWdy[i][j].rho();
             grad_rho_abs = sqrt(sqr(grad_rho_x) + sqr(grad_rho_y));
             grad_rho_criteria = sqrt(SolnBlk.Grid.Cell[i][j].A)*grad_rho_abs/W.rho();
             grad_rho_criteria_max = max(grad_rho_criteria_max, grad_rho_criteria);

             // Evaluate refinement criteria #2 based on the divergence
             // of the velocity vector.
             div_V = SolnBlk.dWdx[i][j].vx() + SolnBlk.dWdy[i][j].vy();
             div_V_criteria = sqrt(SolnBlk.Grid.Cell[i][j].A)*fabs(div_V)/W.a();
             div_V_criteria_max = max(div_V_criteria_max, div_V_criteria);

             // Evaluate refinement criteria #3 based on the curl
             // of the velocity vector.
             curl_V_z = SolnBlk.dWdx[i][j].vy() - SolnBlk.dWdy[i][j].vx(); 
             curl_V_abs = sqrt(sqr(curl_V_z)); 
             curl_V_criteria = sqrt(SolnBlk.Grid.Cell[i][j].A)*curl_V_abs/W.a(); 
             curl_V_criteria_max = max(curl_V_criteria_max, curl_V_criteria);

// 	     Evaluate refinement criteria #4 based on the gradient
// 	     of the Temperature
	     grad_Temp_x = (ONE/( W.rho()* W.Rtot())) * 
	       (SolnBlk.dWdx[i][j].p() - (W.p()/W.rho()) * SolnBlk.dWdx[i][j].rho());
             grad_Temp_y = (ONE/( W.rho()* W.Rtot())) *
	       (SolnBlk.dWdy[i][j].p() - (W.p()/W.rho()) * SolnBlk.dWdy[i][j].rho());
             grad_Temp_abs = sqrt(sqr(grad_Temp_x) + sqr(grad_Temp_y));
             grad_Temp_criteria = sqrt(SolnBlk.Grid.Cell[i][j].A)*grad_Temp_abs/SolnBlk.W[i][j].T();
             grad_Temp_criteria_max = max(grad_Temp_criteria_max, grad_Temp_criteria);

// 	     Evaluate refinement criteria #5 based on the gradient
// 	     based on the gradient of CH4 mass fraction	     	     
	     grad_CH4_x = SolnBlk.dWdx[i][j].c(0);
             grad_CH4_y = SolnBlk.dWdy[i][j].c(0);
             grad_CH4_abs = sqrt(sqr(grad_CH4_x) + sqr(grad_CH4_y));
	     grad_CH4_criteria = sqrt(SolnBlk.Grid.Cell[i][j].A)*grad_CH4_abs;
             grad_CH4_criteria_max = max(grad_CH4_criteria_max, grad_CH4_criteria);

// 	     Evaluate refinement criteria #6 based on the gradient
// 	     based on the gradient of CO2 mass fraction	     	     
	     grad_CO2_x = SolnBlk.dWdx[i][j].c(2);
             grad_CO2_y = SolnBlk.dWdy[i][j].c(2);
             grad_CO2_abs = sqrt(sqr(grad_CO2_x) + sqr(grad_CO2_y));
	     grad_CO2_criteria = sqrt(SolnBlk.Grid.Cell[i][j].A)*grad_CO2_abs;
             grad_CO2_criteria_max = max(grad_CO2_criteria_max, grad_CO2_criteria);

// 	    // Evaluate refinement criteria #7 based on the gradient dudy
// 	    grad_dudy = SolnBlk.dWdy[i][j].v.x;
// 	    grad_dudy_criteria = sqrt(SolnBlk.Grid.Cell[i][j].A)*fabs(grad_dudy)/SolnBlk.W[i][j].a();
// 	    grad_dudy_criteria_max = max(grad_dudy_criteria_max, grad_dudy_criteria);

// 	     grad_pressure_x = SolnBlk.dWdx[i][j].p;
//              grad_pressure_y = SolnBlk.dWdy[i][j].p;
//              grad_pressure_abs = sqrt(sqr(grad_pressure_x) + sqr(grad_pressure_y));
//              grad_pressure_criteria = sqrt(SolnBlk.Grid.Cell[i][j].A)*grad_pressure_abs/SolnBlk.W[i][j].p;
//              grad_pressure_criteria_max = max(grad_pressure_criteria_max, grad_pressure_criteria);

          } /* endif */
       } /* endfor */
    } /* endfor */

    /* Return the refinement criteria. */
    int  refinement_criteria_number(0);
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
    if (IP.Refinement_Criteria_Gradient_Temperature){
      refinement_criteria[refinement_criteria_number] = grad_Temp_criteria_max;
      refinement_criteria_number++;
    }
    if (IP.Refinement_Criteria_Gradient_CH4){
      refinement_criteria[refinement_criteria_number] = grad_CO2_criteria_max;
      refinement_criteria_number++;
    }
    if (IP.Refinement_Criteria_Gradient_CO2){
      refinement_criteria[refinement_criteria_number] = grad_CO2_criteria_max;
      refinement_criteria_number++;
    }

//     refinement_criteria[0] = grad_dudy_criteria_max;
//     refinement_criteria[1] = grad_pressure_criteria_max;

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
void Fix_Refined_Block_Boundaries(Flame2D_Quad_Block &SolnBlk,
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
void Unfix_Refined_Block_Boundaries(Flame2D_Quad_Block &SolnBlk) {

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
void Apply_Boundary_Flux_Corrections(Flame2D_Quad_Block &SolnBlk,
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
void Apply_Boundary_Flux_Corrections_Multistage_Explicit(Flame2D_Quad_Block &SolnBlk,
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
             (CFL_Number*SolnBlk.dt[i][SolnBlk.JCu]/
	      SolnBlk.Grid.Cell[i][SolnBlk.JCu].A)*SolnBlk.FluxN[i];
       } /* endfor */
    } /* endif */

    /* Correct the fluxes at the south boundary as required. */

    if (Number_Neighbours_South_Boundary == 2) {
       for ( i = SolnBlk.ICl ; i <= SolnBlk.ICu ; ++i ) {
          SolnBlk.dUdt[i][SolnBlk.JCl][k_residual] -= 
             (CFL_Number*SolnBlk.dt[i][SolnBlk.JCl]/
	      SolnBlk.Grid.Cell[i][SolnBlk.JCl].A)*SolnBlk.FluxS[i];
       } /* endfor */
    } /* endif */

    /* Correct the fluxes at the east boundary as required. */

    if (Number_Neighbours_East_Boundary == 2) {
       for ( j= SolnBlk.JCl ; j <= SolnBlk.JCu ; ++j ) {
          SolnBlk.dUdt[SolnBlk.ICu][j][k_residual] -= 
             (CFL_Number*SolnBlk.dt[SolnBlk.ICu][j]/
	      SolnBlk.Grid.Cell[SolnBlk.ICu][j].A)*SolnBlk.FluxE[j];
       } /* endfor */
    } /* endif */

    /* Correct the fluxes at the west boundary as required. */

    if (Number_Neighbours_West_Boundary == 2) {
       for ( j= SolnBlk.JCl ; j <= SolnBlk.JCu ; ++j ) {
          SolnBlk.dUdt[SolnBlk.ICl][j][k_residual] -= 
             (CFL_Number*SolnBlk.dt[SolnBlk.ICl][j]/
	      SolnBlk.Grid.Cell[SolnBlk.ICl][j].A)*SolnBlk.FluxW[j];
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
int dUdt_Residual_Evaluation(Flame2D_Quad_Block &SolnBlk,
			     Flame2D_Input_Parameters &Input_Parameters) {


  Vector2D dX;
  Flame2D_pState Wl, Wr;
  Flame2D_State Flux;
  
  Flame2D_pState W, W_face, dWdx, dWdy;  
  int NUM_VAR_FLAME2D = SolnBlk.NumVar();
  double delta_n; 

  /* Additional variables for dual time stepping. */
  double dTime;    // Physical time step
  if (Input_Parameters.Dual_Time_Stepping) {
    dTime = Input_Parameters.dTime;
  } 

  // Modifications for NKS overlap 
  int JCl_overlap = 0, JCu_overlap = 0, ICu_overlap = 0, ICl_overlap = 0;
  if(Input_Parameters.NKS_IP.GMRES_Overlap > 0 ){	
    if (SolnBlk.Grid.BCtypeS[SolnBlk.ICl] == BC_NONE)  JCl_overlap = Input_Parameters.NKS_IP.GMRES_Overlap;
    if (SolnBlk.Grid.BCtypeN[SolnBlk.ICu] == BC_NONE)  JCu_overlap = Input_Parameters.NKS_IP.GMRES_Overlap;
    if (SolnBlk.Grid.BCtypeE[SolnBlk.JCu] == BC_NONE)  ICu_overlap = Input_Parameters.NKS_IP.GMRES_Overlap;
    if (SolnBlk.Grid.BCtypeW[SolnBlk.JCl] == BC_NONE)  ICl_overlap = Input_Parameters.NKS_IP.GMRES_Overlap;    
  }

  /* Perform the linear reconstruction within each cell of the computational grid. */  
  switch(Input_Parameters.i_Reconstruction) {
  case RECONSTRUCTION_GREEN_GAUSS :            
    if(Input_Parameters.i_Viscous_Flux_Evaluation ==VISCOUS_RECONSTRUCTION_DIAMONDPATH_GREEN_GAUSS){
      Linear_Reconstruction_GreenGauss_Diamond(SolnBlk,Input_Parameters.i_Limiter);
    } else {
      Linear_Reconstruction_GreenGauss(SolnBlk, Input_Parameters.i_Limiter);    
    }
    break;
  case RECONSTRUCTION_LEAST_SQUARES :
    if(Input_Parameters.i_Viscous_Flux_Evaluation ==VISCOUS_RECONSTRUCTION_DIAMONDPATH_LEAST_SQUARES){        
      Linear_Reconstruction_LeastSquares_Diamond(SolnBlk, Input_Parameters.i_Limiter);       
    } else {
      Linear_Reconstruction_LeastSquares(SolnBlk, Input_Parameters.i_Limiter);
    }
    break;
  default:
    Linear_Reconstruction_LeastSquares(SolnBlk, Input_Parameters.i_Limiter);
    break;
  }
  
  /********************************************************/
  /* Compute viscous stresses and heat conduction vector if using cell-centered methods. */
  if (SolnBlk.Flow_Type != FLOWTYPE_INVISCID && 
      Input_Parameters.i_Viscous_Flux_Evaluation == VISCOUS_RECONSTRUCTION_ARITHMETIC) {
    Viscous_Calculations(SolnBlk);
  }  
  // radiation evaluation if using opticallyt htin approx
  if ( Input_Parameters.Radiation == RADIATION_OPTICALLY_THIN )
    Radiation_Source_Eval( SolnBlk, Input_Parameters );
  /********************************************************/
  
  /* Evaluate the time rate of change of the solution
     (i.e., the solution residuals) using a second-order
     limited upwind scheme with a variety of flux functions. */
  
  // Add i-direction (zeta-direction) fluxes.
  for (int j  = SolnBlk.JCl-1-JCl_overlap ; j <= SolnBlk.JCu+1+JCu_overlap ; j++ ) {
    SolnBlk.dUdt[SolnBlk.ICl-1][j][0].Vacuum();
      
    for (int i = SolnBlk.ICl-1-ICl_overlap; i <= SolnBlk.ICu+ICu_overlap ; i++ ) {	 
      SolnBlk.dUdt[i+1][j][0].Vacuum();
      
      //if ( j > SolnBlk.JCl-1-JCl_overlap && j < SolnBlk.JCu+1+JCu_overlap ) { 
      if ( j >= SolnBlk.JCl-JCl_overlap && j <= SolnBlk.JCu+JCu_overlap ) { 
	/* Evaluate the cell interface i-direction fluxes. */
	if (i == SolnBlk.ICl-1 && 
	    (SolnBlk.Grid.BCtypeW[j] == BC_REFLECTION ||
	     SolnBlk.Grid.BCtypeW[j] == BC_CHARACTERISTIC ||
	     SolnBlk.Grid.BCtypeW[j] == BC_FREE_SLIP_ISOTHERMAL ||
	     SolnBlk.Grid.BCtypeW[j] == BC_WALL_VISCOUS_ISOTHERMAL ||
	     SolnBlk.Grid.BCtypeW[j] == BC_MOVING_WALL_ISOTHERMAL ||
	     SolnBlk.Grid.BCtypeW[j] == BC_WALL_VISCOUS_HEATFLUX ||
	     SolnBlk.Grid.BCtypeW[j] == BC_MOVING_WALL_HEATFLUX ||
	     SolnBlk.Grid.BCtypeW[j] == BC_1DFLAME_INFLOW ||
	     SolnBlk.Grid.BCtypeW[j] == BC_1DFLAME_OUTFLOW ||
	     SolnBlk.Grid.BCtypeW[j] == BC_2DFLAME_INFLOW ||
	     SolnBlk.Grid.BCtypeW[j] == BC_2DFLAME_OUTFLOW )) {
	  
	  dX = SolnBlk.Grid.xfaceW(i+1, j)-SolnBlk.Grid.Cell[i+1][j].Xc;	      	      
	  Wr.setW( SolnBlk.W[i+1][j] + 
		   (SolnBlk.phi[i+1][j]^SolnBlk.dWdx[i+1][j])*dX.x +
		   (SolnBlk.phi[i+1][j]^SolnBlk.dWdy[i+1][j])*dX.y );
	  
	  if (SolnBlk.Grid.BCtypeW[j] == BC_REFLECTION) {
	    Wl.Reflect(Wr, SolnBlk.Grid.nfaceW(i+1, j));
	  } else if (SolnBlk.Grid.BCtypeW[j] == BC_FREE_SLIP_ISOTHERMAL) {
	    Wl.Free_Slip(Wr,SolnBlk.WoW[j], SolnBlk.Grid.nfaceW(i+1, j),FIXED_TEMPERATURE_WALL);	      
	  } else if (SolnBlk.Grid.BCtypeW[j] == BC_WALL_VISCOUS_ISOTHERMAL) {
	    Wl.No_Slip(Wr,SolnBlk.WoW[j], SolnBlk.Grid.nfaceW(i+1, j),FIXED_TEMPERATURE_WALL);
	  } else if (SolnBlk.Grid.BCtypeW[j] == BC_MOVING_WALL_ISOTHERMAL) {
	    Wl.Moving_Wall(Wr,SolnBlk.WoW[j], SolnBlk.Grid.nfaceW(i+1, j),
			   SolnBlk.Moving_wall_velocity,FIXED_TEMPERATURE_WALL);
	  } else if (SolnBlk.Grid.BCtypeW[j] == BC_WALL_VISCOUS_HEATFLUX) {
	    Wl.No_Slip(Wr,SolnBlk.WoW[j], SolnBlk.Grid.nfaceW(i+1, j), ADIABATIC_WALL);
	  } else if (SolnBlk.Grid.BCtypeW[j] == BC_MOVING_WALL_HEATFLUX) {
	    Wl.Moving_Wall(Wr,SolnBlk.WoW[j], SolnBlk.Grid.nfaceW(i+1, j),
			   SolnBlk.Moving_wall_velocity,ADIABATIC_WALL );
	  } else if (SolnBlk.Grid.BCtypeW[j] == BC_1DFLAME_INFLOW){
	    Wl.BC_1DFlame_Inflow(Wr, 
				 SolnBlk.WoW[j],
				 SolnBlk.W[SolnBlk.ICu][j],
				 SolnBlk.Grid.nfaceW(i+1, j));	    
	  } else if (SolnBlk.Grid.BCtypeW[j] == BC_1DFLAME_OUTFLOW){
	    Wl.BC_1DFlame_Outflow(Wr, 
				  SolnBlk.WoW[j], 
				  SolnBlk.W[SolnBlk.ICu][j],
				  SolnBlk.Grid.nfaceW(i+1, j));
	  } else if (SolnBlk.Grid.BCtypeW[j] == BC_CHARACTERISTIC) {
	    Wl.BC_Characteristic_Pressure(Wr, 
					  SolnBlk.WoW[j], 
					  SolnBlk.Grid.nfaceW(i+1, j));
	  } else {
	    cerr<< "\n Wrong BC in  dUdt_Residual_Evaluation "<< SolnBlk.Grid.BCtypeW[j]; exit(1);
	  }
	  
	} else if (i == SolnBlk.ICu && 
		   (SolnBlk.Grid.BCtypeE[j] == BC_REFLECTION ||
		    SolnBlk.Grid.BCtypeE[j] == BC_CHARACTERISTIC ||  
		    SolnBlk.Grid.BCtypeE[j] == BC_FREE_SLIP_ISOTHERMAL ||
		    SolnBlk.Grid.BCtypeE[j] == BC_WALL_VISCOUS_ISOTHERMAL ||
		    SolnBlk.Grid.BCtypeE[j] == BC_MOVING_WALL_ISOTHERMAL ||
		    SolnBlk.Grid.BCtypeE[j] == BC_WALL_VISCOUS_HEATFLUX ||
		    SolnBlk.Grid.BCtypeE[j] == BC_MOVING_WALL_HEATFLUX ||
		    SolnBlk.Grid.BCtypeE[j] == BC_1DFLAME_INFLOW ||
		    SolnBlk.Grid.BCtypeE[j] == BC_1DFLAME_OUTFLOW ||
		    SolnBlk.Grid.BCtypeE[j] == BC_2DFLAME_INFLOW ||
		    SolnBlk.Grid.BCtypeE[j] == BC_2DFLAME_OUTFLOW )) {
	  dX = SolnBlk.Grid.xfaceE(i, j)-SolnBlk.Grid.Cell[i][j].Xc;
	  Wl.setW( SolnBlk.W[i][j] + 
		   (SolnBlk.phi[i][j]^SolnBlk.dWdx[i][j])*dX.x +
		   (SolnBlk.phi[i][j]^SolnBlk.dWdy[i][j])*dX.y );
	  if (SolnBlk.Grid.BCtypeE[j] == BC_REFLECTION) {
	    Wr.Reflect(Wl, SolnBlk.Grid.nfaceE(i, j));   
	  } else if (SolnBlk.Grid.BCtypeE[j] == BC_FREE_SLIP_ISOTHERMAL) {
	    Wr.Free_Slip(Wl, SolnBlk.WoE[j],SolnBlk.Grid.nfaceE(i, j),FIXED_TEMPERATURE_WALL);
	  } else if (SolnBlk.Grid.BCtypeE[j] == BC_WALL_VISCOUS_ISOTHERMAL) {
	    Wr.No_Slip(Wl, SolnBlk.WoE[j],SolnBlk.Grid.nfaceE(i, j),FIXED_TEMPERATURE_WALL);
	  } else if (SolnBlk.Grid.BCtypeE[j] == BC_MOVING_WALL_ISOTHERMAL) {
	    Wr.Moving_Wall(Wl,SolnBlk.WoE[j], SolnBlk.Grid.nfaceE(i, j),
			   SolnBlk.Moving_wall_velocity,FIXED_TEMPERATURE_WALL);
	  } else if (SolnBlk.Grid.BCtypeE[j] == BC_WALL_VISCOUS_HEATFLUX) {
	    Wr.No_Slip(Wl, SolnBlk.WoE[j],SolnBlk.Grid.nfaceE(i, j),ADIABATIC_WALL);
	  } else if (SolnBlk.Grid.BCtypeE[j] == BC_MOVING_WALL_HEATFLUX) {
	    Wr.Moving_Wall(Wl,SolnBlk.WoE[j], SolnBlk.Grid.nfaceE(i, j),
			   SolnBlk.Moving_wall_velocity,ADIABATIC_WALL);
	  } else if (SolnBlk.Grid.BCtypeE[j] == BC_1DFLAME_INFLOW){
	    Wr.BC_1DFlame_Inflow(Wl, 
				 SolnBlk.WoE[j],
				 SolnBlk.W[SolnBlk.ICl][j],
				 SolnBlk.Grid.nfaceE(i, j));
	  } else if (SolnBlk.Grid.BCtypeE[j] == BC_1DFLAME_OUTFLOW){
	    Wr.BC_1DFlame_Outflow(Wl, 
				  SolnBlk.WoE[j],
				  SolnBlk.W[SolnBlk.ICl][j], 
				  SolnBlk.Grid.nfaceE(i, j));
	  } else if (SolnBlk.Grid.BCtypeE[j] == BC_CHARACTERISTIC) { 
	    Wr.BC_Characteristic_Pressure(Wl, 
					  SolnBlk.WoE[j], 
					  SolnBlk.Grid.nfaceE(i, j));
	  } else {
	    cerr<< "\n Wrong BC in  dUdt_Residual_Evaluation "<< SolnBlk.Grid.BCtypeE[j]; exit(1);
	  }

	} else {            
	  dX = SolnBlk.Grid.xfaceE(i, j)-SolnBlk.Grid.Cell[i][j].Xc;
	  Wl.setW( SolnBlk.W[i][j] + 
		   (SolnBlk.phi[i][j]^SolnBlk.dWdx[i][j])*dX.x +
		   (SolnBlk.phi[i][j]^SolnBlk.dWdy[i][j])*dX.y );
	  dX = SolnBlk.Grid.xfaceW(i+1, j)-SolnBlk.Grid.Cell[i+1][j].Xc;
	  Wr.setW( SolnBlk.W[i+1][j] + 
		   (SolnBlk.phi[i+1][j]^SolnBlk.dWdx[i+1][j])*dX.x +
		   (SolnBlk.phi[i+1][j]^SolnBlk.dWdy[i+1][j])*dX.y );
	} /* endif */

    
	// Spacing for Preconditioner 
	if(Input_Parameters.Preconditioning){
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
	  Flux.FluxRoe_n(Wl, Wr, SolnBlk.Grid.nfaceE(i, j),
			 Input_Parameters.Preconditioning,SolnBlk.Flow_Type,delta_n);
	  break;
	  //                case FLUX_FUNCTION_RUSANOV :
	  //                  Flux = FluxRusanov_n(Wl, Wr, SolnBlk.Grid.nfaceE(i, j));
	  //                  break;
	case FLUX_FUNCTION_HLLE :  
	  Flux.FluxHLLE_n(Wl, Wr, SolnBlk.Grid.nfaceE(i, j),Input_Parameters.Preconditioning);
	  break;
	case FLUX_FUNCTION_LINDE :
	  Flux.FluxLinde_n(Wl, Wr, SolnBlk.Grid.nfaceE(i, j));
	  break; 
	case FLUX_FUNCTION_AUSM_PLUS_UP :
	  Flux.FluxAUSMplus_up_n(Wl, Wr, SolnBlk.Grid.nfaceE(i, j));
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
	    
	  case VISCOUS_RECONSTRUCTION_ARITHMETIC :
	    //Arithmetic Mean
	    Flux -= Viscous_FluxArithmetic_n(SolnBlk.U[i][j],
					     SolnBlk.dWdx[i][j],
					     SolnBlk.dWdy[i][j],
					     SolnBlk.U[i+1][j],
					     SolnBlk.dWdx[i+1][j],
					     SolnBlk.dWdy[i+1][j],
					     SolnBlk.Grid.nfaceE(i,j));
	    break;
	  case VISCOUS_RECONSTRUCTION_MEAN_GRADIENT :
	    // Determine state at the interface.
	    W = HALF*(SolnBlk.W[i][j] + SolnBlk.W[i+1][j]);
	    // Determine gradient of the primitive variables at the interface.
	    dWdx = HALF*(SolnBlk.dWdx[i][j] + SolnBlk.dWdx[i+1][j]);
	    dWdy = HALF*(SolnBlk.dWdy[i][j] + SolnBlk.dWdy[i+1][j]);
	    // Compute the viscous flux.
	    Flux -= Viscous_Flux_n(W,dWdx,dWdy,
				   SolnBlk.Axisymmetric,
				   SolnBlk.Grid.xfaceE(i,j),
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
	    
	    SolnBlk.dWdx_faceE[i][j]=dWdx;
	    SolnBlk.dWdy_faceE[i][j]=dWdy;
	    SolnBlk.dWdx_faceW[i+1][j]=dWdx;
	    SolnBlk.dWdy_faceW[i+1][j]=dWdy;
	    
	    // Compute the viscous flux.
	    Flux -= Viscous_Flux_n(W,dWdx,dWdy,
				   SolnBlk.Axisymmetric,
				   SolnBlk.Grid.xfaceE(i,j),
				   SolnBlk.Grid.nfaceE(i,j));
	    break;
	    
	  case VISCOUS_RECONSTRUCTION_DIAMONDPATH_LEAST_SQUARES :	      
	  case VISCOUS_RECONSTRUCTION_DIAMONDPATH_GREEN_GAUSS :
	    W_face = HALF*(SolnBlk.Wn(i+1,j+1) +SolnBlk.Wn(i+1,j)); // IS this right in general ??( Wn_NE + Wn_SE ) 
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
		
	  /* Include source terms associated with the finite-rate chemistry and 
             turbulence/chemistry interactions */ 
	if (SolnBlk.W[i][j].React.reactset_flag != NO_REACTIONS) {	  
	  SolnBlk.dUdt[i][j][0] += SolnBlk.W[i][j].Sw(SolnBlk.W[i][j].React.reactset_flag);
	} /* endif */
	
	// 	  // HEAT SOURCE 
	// 	  if(Input_Parameters.Heat_Source > ZERO){	
	// 	    SolnBlk.dUdt[i][j][0].E = (Input_Parameters.Heat_Source - SolnBlk.U[i][j].E);
	// 	  }
		
	/* Include source terms associated with gravity */
	if (SolnBlk.Gravity) {	 
	  SolnBlk.dUdt[i][j][0] += SolnBlk.W[i][j].Sg();
	}   
	
	/* Include source terms associated with radiation */
	SolnBlk.dUdt[i][j][0].E += SolnBlk.Srad[i][j];
	
	/*****************************************/
	/*****************************************/
	
	/* Save west and east face boundary flux. */	
	if (i == SolnBlk.ICl-1) {
	  SolnBlk.FluxW[j] = -Flux*SolnBlk.Grid.lfaceW(i+1, j);
	} else if (i == SolnBlk.ICu) {
	  SolnBlk.FluxE[j] = Flux*SolnBlk.Grid.lfaceE(i, j);  
	} 
	
      } // end if 
    } // end for i
    
    if ( j > SolnBlk.JCl-1-JCl_overlap && j < SolnBlk.JCu+1+JCu_overlap ) {  
      SolnBlk.dUdt[SolnBlk.ICu+1][j][0].Vacuum();
      SolnBlk.dUdt[SolnBlk.ICl-1][j][0].Vacuum();
    }
  } // end for j
      
  
  // Add j-direction (eta-direction) fluxes.
  for ( int i = SolnBlk.ICl-ICl_overlap ; i <= SolnBlk.ICu+ICu_overlap ; ++i ) {
    for ( int j = SolnBlk.JCl-1-JCl_overlap ; j <= SolnBlk.JCu+JCu_overlap ; ++j ) {
      
      /* Evaluate the cell interface j-direction fluxes. */ 
      if (j == SolnBlk.JCl-1 && 
	  (SolnBlk.Grid.BCtypeS[i] == BC_REFLECTION ||
	   SolnBlk.Grid.BCtypeS[i] == BC_CHARACTERISTIC ||
	   SolnBlk.Grid.BCtypeS[i] == BC_FREE_SLIP_ISOTHERMAL ||
	   SolnBlk.Grid.BCtypeS[i] == BC_WALL_VISCOUS_ISOTHERMAL ||
	   SolnBlk.Grid.BCtypeS[i] == BC_MOVING_WALL_ISOTHERMAL ||
	   SolnBlk.Grid.BCtypeS[i] == BC_WALL_VISCOUS_HEATFLUX ||
	   SolnBlk.Grid.BCtypeS[i] == BC_MOVING_WALL_HEATFLUX || 
	   SolnBlk.Grid.BCtypeS[i] == BC_1DFLAME_INFLOW ||
	   SolnBlk.Grid.BCtypeS[i] == BC_1DFLAME_OUTFLOW || 
	   SolnBlk.Grid.BCtypeS[i] == BC_2DFLAME_INFLOW ||
	   SolnBlk.Grid.BCtypeS[i] == BC_2DFLAME_OUTFLOW )) {
	
	  dX = SolnBlk.Grid.xfaceS(i, j+1)-SolnBlk.Grid.Cell[i][j+1].Xc;
	  Wr = SolnBlk.W[i][j+1] +
	    (SolnBlk.phi[i][j+1]^SolnBlk.dWdx[i][j+1])*dX.x +
	    (SolnBlk.phi[i][j+1]^SolnBlk.dWdy[i][j+1])*dX.y;
	  if (SolnBlk.Grid.BCtypeS[i] == BC_REFLECTION) {
	    Wl = Reflect(Wr, SolnBlk.Grid.nfaceS(i, j+1)); 
	  } else if (SolnBlk.Grid.BCtypeS[i] == BC_FREE_SLIP_ISOTHERMAL) {
	    Wl = Free_Slip(Wr,SolnBlk.WoS[i], SolnBlk.Grid.nfaceS(i, j+1),FIXED_TEMPERATURE_WALL);
	  } else if (SolnBlk.Grid.BCtypeS[i] == BC_WALL_VISCOUS_ISOTHERMAL) {
	    Wl = No_Slip(Wr,SolnBlk.WoS[i], SolnBlk.Grid.nfaceS(i, j+1),FIXED_TEMPERATURE_WALL);
	  } else if (SolnBlk.Grid.BCtypeS[i] == BC_MOVING_WALL_ISOTHERMAL) {
	    Wl = Moving_Wall(Wr,SolnBlk.WoS[i], SolnBlk.Grid.nfaceS(i, j+1),
			     SolnBlk.Moving_wall_velocity,FIXED_TEMPERATURE_WALL);
	  } else if (SolnBlk.Grid.BCtypeS[i] == BC_WALL_VISCOUS_HEATFLUX) {
	    Wl = No_Slip(Wr,SolnBlk.WoS[i], SolnBlk.Grid.nfaceS(i, j+1),ADIABATIC_WALL);
	  } else if (SolnBlk.Grid.BCtypeS[i] == BC_MOVING_WALL_HEATFLUX) {
	    Wl = Moving_Wall(Wr,SolnBlk.WoS[i], SolnBlk.Grid.nfaceS(i, j+1),
			     SolnBlk.Moving_wall_velocity,ADIABATIC_WALL); 
	  } else if (SolnBlk.Grid.BCtypeS[i] == BC_1DFLAME_INFLOW){
	    Wl = BC_1DFlame_Inflow(Wr, 
				 SolnBlk.WoS[i], 
				 SolnBlk.W[i][SolnBlk.JCu],
				 SolnBlk.Grid.nfaceS(i, j+1));
	  } else if (SolnBlk.Grid.BCtypeS[i] == BC_2DFLAME_INFLOW){
	    Wl = BC_2DFlame_Inflow(Wr, 
				   SolnBlk.WoS[i], 
				   SolnBlk.Grid.nfaceS(i, j+1));
	  } else if (SolnBlk.Grid.BCtypeS[i] == BC_1DFLAME_OUTFLOW){
	    Wl = BC_1DFlame_Outflow(Wr, 
				    SolnBlk.WoS[i], 
				    SolnBlk.W[i][SolnBlk.JCu],
				    SolnBlk.Grid.nfaceS(i, j+1));
	  } else if(SolnBlk.Grid.BCtypeS[i] == BC_CHARACTERISTIC)  { 
	    Wl = BC_Characteristic_Pressure(Wr, 
					    SolnBlk.WoS[i], 
					    SolnBlk.Grid.nfaceS(i, j+1));
	  } else {
	    cerr<< "\n Wrong BC in  dUdt_Residual_Evaluation "<< SolnBlk.Grid.BCtypeS[i]; exit(1);
	  } 

	} else if (j == SolnBlk.JCu && 
		   (SolnBlk.Grid.BCtypeN[i] == BC_REFLECTION ||
		    SolnBlk.Grid.BCtypeN[i] == BC_CHARACTERISTIC ||
		    SolnBlk.Grid.BCtypeN[i] == BC_FREE_SLIP_ISOTHERMAL || 
		    SolnBlk.Grid.BCtypeN[i] == BC_WALL_VISCOUS_ISOTHERMAL ||
		    SolnBlk.Grid.BCtypeN[i] == BC_MOVING_WALL_ISOTHERMAL ||
		    SolnBlk.Grid.BCtypeN[i] == BC_WALL_VISCOUS_HEATFLUX ||
		    SolnBlk.Grid.BCtypeN[i] == BC_MOVING_WALL_HEATFLUX ||
		    SolnBlk.Grid.BCtypeN[i] == BC_1DFLAME_INFLOW ||
		    SolnBlk.Grid.BCtypeN[i] == BC_1DFLAME_OUTFLOW||
		    SolnBlk.Grid.BCtypeN[i] == BC_2DFLAME_INFLOW ||
		    SolnBlk.Grid.BCtypeN[i] == BC_2DFLAME_OUTFLOW )) {
	  
	  dX = SolnBlk.Grid.xfaceN(i, j)-SolnBlk.Grid.Cell[i][j].Xc;
	  Wl = SolnBlk.W[i][j] + 
	    (SolnBlk.phi[i][j]^SolnBlk.dWdx[i][j])*dX.x +
	    (SolnBlk.phi[i][j]^SolnBlk.dWdy[i][j])*dX.y;
	  if (SolnBlk.Grid.BCtypeN[i] == BC_REFLECTION) {
	    Wr = Reflect(Wl, SolnBlk.Grid.nfaceN(i, j));	
	  } else if (SolnBlk.Grid.BCtypeN[i] == BC_FREE_SLIP_ISOTHERMAL) {
	    Wr = Free_Slip(Wl, SolnBlk.WoN[i], SolnBlk.Grid.nfaceN(i, j),FIXED_TEMPERATURE_WALL);
	  } else if (SolnBlk.Grid.BCtypeN[i] == BC_WALL_VISCOUS_ISOTHERMAL) {
	    Wr = No_Slip(Wl, SolnBlk.WoN[i], SolnBlk.Grid.nfaceN(i, j),FIXED_TEMPERATURE_WALL);
	  } else if (SolnBlk.Grid.BCtypeN[i] == BC_MOVING_WALL_ISOTHERMAL) {
	    Wr = Moving_Wall(Wl,SolnBlk.WoN[i],  SolnBlk.Grid.nfaceN(i, j),
			     SolnBlk.Moving_wall_velocity,FIXED_TEMPERATURE_WALL);
	  } else if (SolnBlk.Grid.BCtypeN[i] == BC_WALL_VISCOUS_HEATFLUX) {
	    Wr = No_Slip(Wl, SolnBlk.WoN[i], SolnBlk.Grid.nfaceN(i, j),ADIABATIC_WALL );
	  } else if (SolnBlk.Grid.BCtypeN[i] == BC_MOVING_WALL_HEATFLUX) {
	    Wr = Moving_Wall(Wl,SolnBlk.WoN[i],  SolnBlk.Grid.nfaceN(i, j),
			     SolnBlk.Moving_wall_velocity,ADIABATIC_WALL ); 
	  } else if (SolnBlk.Grid.BCtypeN[i] == BC_1DFLAME_INFLOW){
	    Wr = BC_1DFlame_Inflow(Wl, 
				 SolnBlk.WoN[i], 
				 SolnBlk.W[i][SolnBlk.JCl],
				 SolnBlk.Grid.nfaceN(i, j));
	  } else if (SolnBlk.Grid.BCtypeN[i] == BC_1DFLAME_OUTFLOW){
	    Wr = BC_1DFlame_Outflow(Wl, 
				  SolnBlk.WoN[i], 
				  SolnBlk.W[i][SolnBlk.JCl],
				  SolnBlk.Grid.nfaceN(i, j));
	  } else if (SolnBlk.Grid.BCtypeN[i] == BC_2DFLAME_OUTFLOW){
	    Wr = BC_2DFlame_Outflow(Wl, 
				    SolnBlk.WoN[i], 
				    SolnBlk.Grid.nfaceN(i, j));	    
	  } else if(SolnBlk.Grid.BCtypeN[i] == BC_CHARACTERISTIC)  { 
	    Wr = BC_Characteristic_Pressure(Wl, 
					    SolnBlk.WoN[i], 
					    SolnBlk.Grid.nfaceN(i, j));
	  } else {
	    cerr<< "\n Wrong BC in  dUdt_Residual_Evaluation "<< SolnBlk.Grid.BCtypeN[i]; exit(1);
	  }  
	  
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
			   Input_Parameters.Preconditioning,SolnBlk.Flow_Type,delta_n);
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
	case FLUX_FUNCTION_AUSM_PLUS_UP :
	  Flux = FluxAUSMplus_up_n(Wl, Wr, SolnBlk.Grid.nfaceN(i, j));
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
	  case VISCOUS_RECONSTRUCTION_ARITHMETIC :
	    //Arithmetic Mean
	    Flux -= Viscous_FluxArithmetic_n(SolnBlk.U[i][j],
                                             SolnBlk.dWdx[i][j],
                                             SolnBlk.dWdy[i][j],
                                             SolnBlk.U[i][j+1],
                                             SolnBlk.dWdx[i][j+1],
                                             SolnBlk.dWdy[i][j+1],
					     SolnBlk.Grid.nfaceN(i,j));
	    break;
	  case VISCOUS_RECONSTRUCTION_MEAN_GRADIENT :
	    // Determine state at the interface.
	    W = HALF*(SolnBlk.W[i][j] + SolnBlk.W[i][j+1]);
	    // Determine gradient of the primitive variables at the interface.
	    dWdx = HALF*(SolnBlk.dWdx[i][j] + SolnBlk.dWdx[i][j+1]);
	    dWdy = HALF*(SolnBlk.dWdy[i][j] + SolnBlk.dWdy[i][j+1]);
	    // Compute the viscous flux.
	    Flux -= Viscous_Flux_n(W,dWdx,dWdy,
                                   SolnBlk.Axisymmetric,
				   SolnBlk.Grid.xfaceN(i,j),
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

	    SolnBlk.dWdx_faceN[i][j]=dWdx;
	    SolnBlk.dWdy_faceN[i][j]=dWdy;
	    SolnBlk.dWdx_faceS[i][j+1]=dWdx;
	    SolnBlk.dWdy_faceS[i][j+1]=dWdy;

	    // Compute the viscous flux.
	    Flux -= Viscous_Flux_n(W,dWdx,dWdy,
                                   SolnBlk.Axisymmetric,
				   SolnBlk.Grid.xfaceN(i,j),
				   SolnBlk.Grid.nfaceN(i,j));
	 
	    break;

	  case VISCOUS_RECONSTRUCTION_DIAMONDPATH_LEAST_SQUARES :	      
	  case VISCOUS_RECONSTRUCTION_DIAMONDPATH_GREEN_GAUSS :	  	  
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
	
	SolnBlk.dUdt[i][j][0] -=  Flux*SolnBlk.Grid.lfaceN(i, j)/
	  SolnBlk.Grid.Cell[i][j].A;
	SolnBlk.dUdt[i][j+1][0] += Flux*SolnBlk.Grid.lfaceS(i, j+1)/
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


    /* Residual successfully calculated. */

//     //HACK FOR 1D FLAME
//     if( SolnBlk.Grid.BCtypeE[SolnBlk.JCu] == BC_FLAME_OUTFLOW 
// 	||  SolnBlk.Grid.BCtypeW[SolnBlk.JCu] == BC_FLAME_INFLOW ){ 
//       SolnBlk.set_v_zero();
//     }

}

/********************************************************
 * Routine: dUdt_Multistage_Explicit                    *
 *                                                      *
 * This routine determines the solution residuals for a *
 * given stage of a variety of multi-stage explicit     *
 * time integration schemes for a given solution block. *
 *                                                      *
 ********************************************************/ 
int dUdt_Multistage_Explicit(Flame2D_Quad_Block &SolnBlk,
                             const int i_stage,
                             Flame2D_Input_Parameters &Input_Parameters) {

    int i, j, k_residual;
    double omega;
    Vector2D dX;
    Flame2D_pState Wl, Wr;
    Flame2D_cState Flux;

    Flame2D_pState W, W_face, dWdx, dWdy;  
    int NUM_VAR_FLAME2D = SolnBlk.NumVar();
    double delta_n;

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
       if(Input_Parameters.i_Viscous_Flux_Evaluation ==VISCOUS_RECONSTRUCTION_DIAMONDPATH_GREEN_GAUSS){
          Linear_Reconstruction_GreenGauss_Diamond(SolnBlk,Input_Parameters.i_Limiter);
       }else{
          Linear_Reconstruction_GreenGauss(SolnBlk,Input_Parameters.i_Limiter);    
       }
       break;
    case RECONSTRUCTION_LEAST_SQUARES :
       if(Input_Parameters.i_Viscous_Flux_Evaluation ==VISCOUS_RECONSTRUCTION_DIAMONDPATH_LEAST_SQUARES){          
          Linear_Reconstruction_LeastSquares_Diamond(SolnBlk, Input_Parameters.i_Limiter);          
       }else{
          Linear_Reconstruction_LeastSquares(SolnBlk, Input_Parameters.i_Limiter);
       }
       break;
    default:
       Linear_Reconstruction_LeastSquares(SolnBlk,Input_Parameters.i_Limiter);
       break;
    } /* endswitch */

    /********************************************************/
    /********************************************************/
  
    /* Compute viscous stresses and heat conduction vector if using cell-centered methods. */
    if (SolnBlk.Flow_Type != FLOWTYPE_INVISCID && 
        Input_Parameters.i_Viscous_Flux_Evaluation == VISCOUS_RECONSTRUCTION_ARITHMETIC) {
       Viscous_Calculations(SolnBlk);
    }
    // radiation evaluation if using opticallyt htin approx
    if ( Input_Parameters.Radiation == RADIATION_OPTICALLY_THIN )
      Radiation_Source_Eval( SolnBlk, Input_Parameters );
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
	     SolnBlk.Grid.BCtypeW[j] == BC_FREE_SLIP_ISOTHERMAL ||
	     SolnBlk.Grid.BCtypeW[j] == BC_WALL_VISCOUS_ISOTHERMAL ||
	     SolnBlk.Grid.BCtypeW[j] == BC_MOVING_WALL_ISOTHERMAL ||
	     SolnBlk.Grid.BCtypeW[j] == BC_WALL_VISCOUS_HEATFLUX ||
	     SolnBlk.Grid.BCtypeW[j] == BC_MOVING_WALL_HEATFLUX ||
	     SolnBlk.Grid.BCtypeW[j] == BC_1DFLAME_INFLOW ||
	     SolnBlk.Grid.BCtypeW[j] == BC_1DFLAME_OUTFLOW ||
	     SolnBlk.Grid.BCtypeW[j] == BC_2DFLAME_INFLOW ||
	     SolnBlk.Grid.BCtypeW[j] == BC_2DFLAME_OUTFLOW )) {
	  
	  dX = SolnBlk.Grid.xfaceW(i+1, j)-SolnBlk.Grid.Cell[i+1][j].Xc;	      	      
	  Wr = SolnBlk.W[i+1][j] + 
	    (SolnBlk.phi[i+1][j]^SolnBlk.dWdx[i+1][j])*dX.x +
	    (SolnBlk.phi[i+1][j]^SolnBlk.dWdy[i+1][j])*dX.y;	
	  
	  if (SolnBlk.Grid.BCtypeW[j] == BC_REFLECTION) {
	    Wl = Reflect(Wr, SolnBlk.Grid.nfaceW(i+1, j));
	  } else if (SolnBlk.Grid.BCtypeW[j] == BC_FREE_SLIP_ISOTHERMAL) {
	      Wl = Free_Slip(Wr,SolnBlk.WoW[j], SolnBlk.Grid.nfaceW(i+1, j),FIXED_TEMPERATURE_WALL);	      
	  } else if (SolnBlk.Grid.BCtypeW[j] == BC_WALL_VISCOUS_ISOTHERMAL) {
	    Wl = No_Slip(Wr,SolnBlk.WoW[j], SolnBlk.Grid.nfaceW(i+1, j),FIXED_TEMPERATURE_WALL);
	  } else if (SolnBlk.Grid.BCtypeW[j] == BC_MOVING_WALL_ISOTHERMAL) {
	    Wl = Moving_Wall(Wr,SolnBlk.WoW[j], SolnBlk.Grid.nfaceW(i+1, j),
			     SolnBlk.Moving_wall_velocity,FIXED_TEMPERATURE_WALL);
	  } else if (SolnBlk.Grid.BCtypeW[j] == BC_WALL_VISCOUS_HEATFLUX) {
	    Wl = No_Slip(Wr,SolnBlk.WoW[j], SolnBlk.Grid.nfaceW(i+1, j), ADIABATIC_WALL);
	  } else if (SolnBlk.Grid.BCtypeW[j] == BC_MOVING_WALL_HEATFLUX) {
	    Wl = Moving_Wall(Wr,SolnBlk.WoW[j], SolnBlk.Grid.nfaceW(i+1, j),
			     SolnBlk.Moving_wall_velocity,ADIABATIC_WALL );
	  } else if (SolnBlk.Grid.BCtypeW[j] == BC_1DFLAME_INFLOW){
	    Wl = BC_1DFlame_Inflow(Wr, 
				   SolnBlk.WoW[j],
				   SolnBlk.W[SolnBlk.ICu][j],
				   SolnBlk.Grid.nfaceW(i+1, j));	    
	  } else if (SolnBlk.Grid.BCtypeW[j] == BC_1DFLAME_OUTFLOW){
	    Wl = BC_1DFlame_Outflow(Wr, 
				    SolnBlk.WoW[j], 
				    SolnBlk.W[SolnBlk.ICu][j],
				    SolnBlk.Grid.nfaceW(i+1, j));
	  } else if (SolnBlk.Grid.BCtypeW[j] == BC_CHARACTERISTIC) {
	    Wl = BC_Characteristic_Pressure(Wr, 
					    SolnBlk.WoW[j], 
					    SolnBlk.Grid.nfaceW(i+1, j));
	  } else {
	    cerr<< "\n Wrong BC in  dUdt_Residual_Evaluation "<< SolnBlk.Grid.BCtypeW[j]; exit(1);
	  }
	  
	} else if (i == SolnBlk.ICu && 
		   (SolnBlk.Grid.BCtypeE[j] == BC_REFLECTION ||
		    SolnBlk.Grid.BCtypeE[j] == BC_CHARACTERISTIC ||  
		    SolnBlk.Grid.BCtypeE[j] == BC_FREE_SLIP_ISOTHERMAL ||
		    SolnBlk.Grid.BCtypeE[j] == BC_WALL_VISCOUS_ISOTHERMAL ||
		    SolnBlk.Grid.BCtypeE[j] == BC_MOVING_WALL_ISOTHERMAL ||
		    SolnBlk.Grid.BCtypeE[j] == BC_WALL_VISCOUS_HEATFLUX ||
		    SolnBlk.Grid.BCtypeE[j] == BC_MOVING_WALL_HEATFLUX ||
		    SolnBlk.Grid.BCtypeE[j] == BC_1DFLAME_INFLOW ||
		    SolnBlk.Grid.BCtypeE[j] == BC_1DFLAME_OUTFLOW ||
		    SolnBlk.Grid.BCtypeE[j] == BC_2DFLAME_INFLOW ||
		    SolnBlk.Grid.BCtypeE[j] == BC_2DFLAME_OUTFLOW )) {
	  dX = SolnBlk.Grid.xfaceE(i, j)-SolnBlk.Grid.Cell[i][j].Xc;
	  Wl = SolnBlk.W[i][j] + 
	    (SolnBlk.phi[i][j]^SolnBlk.dWdx[i][j])*dX.x +
	    (SolnBlk.phi[i][j]^SolnBlk.dWdy[i][j])*dX.y;	    
	  if (SolnBlk.Grid.BCtypeE[j] == BC_REFLECTION) {
	    Wr = Reflect(Wl, SolnBlk.Grid.nfaceE(i, j));   
	  } else if (SolnBlk.Grid.BCtypeE[j] == BC_FREE_SLIP_ISOTHERMAL) {
	    Wr = Free_Slip(Wl, SolnBlk.WoE[j],SolnBlk.Grid.nfaceE(i, j),FIXED_TEMPERATURE_WALL);
	  } else if (SolnBlk.Grid.BCtypeE[j] == BC_WALL_VISCOUS_ISOTHERMAL) {
	    Wr = No_Slip(Wl, SolnBlk.WoE[j],SolnBlk.Grid.nfaceE(i, j),FIXED_TEMPERATURE_WALL);
	  } else if (SolnBlk.Grid.BCtypeE[j] == BC_MOVING_WALL_ISOTHERMAL) {
	    Wr = Moving_Wall(Wl,SolnBlk.WoE[j], SolnBlk.Grid.nfaceE(i, j),
			     SolnBlk.Moving_wall_velocity,FIXED_TEMPERATURE_WALL);
	  } else if (SolnBlk.Grid.BCtypeE[j] == BC_WALL_VISCOUS_HEATFLUX) {
	    Wr = No_Slip(Wl, SolnBlk.WoE[j],SolnBlk.Grid.nfaceE(i, j),ADIABATIC_WALL);
	  } else if (SolnBlk.Grid.BCtypeE[j] == BC_MOVING_WALL_HEATFLUX) {
	    Wr = Moving_Wall(Wl,SolnBlk.WoE[j], SolnBlk.Grid.nfaceE(i, j),
			       SolnBlk.Moving_wall_velocity,ADIABATIC_WALL);
	  } else if (SolnBlk.Grid.BCtypeE[j] == BC_1DFLAME_INFLOW){
	    Wr = BC_1DFlame_Inflow(Wl, 
				   SolnBlk.WoE[j],
				   SolnBlk.W[SolnBlk.ICl][j],
				   SolnBlk.Grid.nfaceE(i, j));
	  } else if (SolnBlk.Grid.BCtypeE[j] == BC_1DFLAME_OUTFLOW){
	    Wr = BC_1DFlame_Outflow(Wl, 
				    SolnBlk.WoE[j],
				    SolnBlk.W[SolnBlk.ICl][j], 
				    SolnBlk.Grid.nfaceE(i, j));
	  } else if (SolnBlk.Grid.BCtypeE[j] == BC_CHARACTERISTIC) { 
	    Wr = BC_Characteristic_Pressure(Wl, 
					    SolnBlk.WoE[j], 
					    SolnBlk.Grid.nfaceE(i, j));
	  } else {
	    cerr<< "\n Wrong BC in  dUdt_Residual_Evaluation "<< SolnBlk.Grid.BCtypeE[j]; exit(1);
	  }

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
	  if(Input_Parameters.Preconditioning){
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
			     Input_Parameters.Preconditioning,SolnBlk.Flow_Type,delta_n);
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
	  case FLUX_FUNCTION_AUSM_PLUS_UP :
	    Flux = FluxAUSMplus_up_n(Wl, Wr, SolnBlk.Grid.nfaceE(i, j));
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
	    case VISCOUS_RECONSTRUCTION_ARITHMETIC :
	      //Arithmetic Mean
	      Flux -= Viscous_FluxArithmetic_n(SolnBlk.U[i][j],
                                               SolnBlk.dWdx[i][j],
                                               SolnBlk.dWdy[i][j],
                                               SolnBlk.U[i+1][j],
                                               SolnBlk.dWdx[i+1][j],
                                               SolnBlk.dWdy[i+1][j],
	        		               SolnBlk.Grid.nfaceE(i,j));  // d/dr Frv
	      break;
	    case VISCOUS_RECONSTRUCTION_MEAN_GRADIENT :
	      // Determine state at the interface.
	      W = HALF*(SolnBlk.W[i][j] + SolnBlk.W[i+1][j]);
	      // Determine gradient of the primitive variables at the interface.
	      dWdx = HALF*(SolnBlk.dWdx[i][j] + SolnBlk.dWdx[i+1][j]);
	      dWdy = HALF*(SolnBlk.dWdy[i][j] + SolnBlk.dWdy[i+1][j]);
	      // Compute the viscous flux.
	      Flux -= Viscous_Flux_n(W,dWdx,dWdy,
                                     SolnBlk.Axisymmetric,
				     SolnBlk.Grid.xfaceE(i,j),
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

	    case VISCOUS_RECONSTRUCTION_DIAMONDPATH_LEAST_SQUARES :
	    case VISCOUS_RECONSTRUCTION_DIAMONDPATH_GREEN_GAUSS :
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

	  /* Include source terms associated with the finite-rate chemistry and 
             turbulence/chemistry interactions */ 
	  if (SolnBlk.W[i][j].React.reactset_flag != NO_REACTIONS) {	 
	    //rho*omega_dot
	    //
	    SolnBlk.dUdt[i][j][k_residual] += (Input_Parameters.CFL_Number*SolnBlk.dt[i][j])* 
	      SolnBlk.W[i][j].Sw(SolnBlk.W[i][j].React.reactset_flag );
	    //Used as ignitor (NOT FINISHED YET)
	    //if(SolnBlk.Heat_source){
	    //  SolnBlk.dUdt[i][j][k_residual].E = Heat_source(SolnBlk.Grid.Cell[i][j].Xc);
	    //}     
	  }

// 	  if(Input_Parameters.Heat_Source > ZERO){	
// 	    SolnBlk.dUdt[i][j][k_residual].E = (Input_Parameters.Heat_Source - SolnBlk.U[i][j].E);
// 	  }     
 
          /* Include source terms associated with gravity */
	  if (SolnBlk.Gravity) {	 
	      SolnBlk.dUdt[i][j][k_residual] += (Input_Parameters.CFL_Number*SolnBlk.dt[i][j])*
                                                SolnBlk.W[i][j].Sg();
          } /* endif */

	  /* Include source terms associated with radiation */
	  SolnBlk.dUdt[i][j][k_residual].E += 
	    (Input_Parameters.CFL_Number*SolnBlk.dt[i][j])*SolnBlk.Srad[i][j];


	  /* Save west and east face boundary flux. */
	  // USED for AMR 	
	  if (i == SolnBlk.ICl-1) {
	    SolnBlk.FluxW[j] = -Flux*SolnBlk.Grid.lfaceW(i+1, j);
	  } else if (i == SolnBlk.ICu) {
	    SolnBlk.FluxE[j] = Flux*SolnBlk.Grid.lfaceE(i, j);
	  } 
	
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
	     SolnBlk.Grid.BCtypeS[i] == BC_FREE_SLIP_ISOTHERMAL ||
	     SolnBlk.Grid.BCtypeS[i] == BC_WALL_VISCOUS_ISOTHERMAL ||
	     SolnBlk.Grid.BCtypeS[i] == BC_MOVING_WALL_ISOTHERMAL ||
	     SolnBlk.Grid.BCtypeS[i] == BC_WALL_VISCOUS_HEATFLUX ||
	     SolnBlk.Grid.BCtypeS[i] == BC_MOVING_WALL_HEATFLUX || 
	     SolnBlk.Grid.BCtypeS[i] == BC_1DFLAME_INFLOW ||
	     SolnBlk.Grid.BCtypeS[i] == BC_1DFLAME_OUTFLOW || 
	     SolnBlk.Grid.BCtypeS[i] == BC_2DFLAME_INFLOW ||
	     SolnBlk.Grid.BCtypeS[i] == BC_2DFLAME_OUTFLOW )) {
	
	  dX = SolnBlk.Grid.xfaceS(i, j+1)-SolnBlk.Grid.Cell[i][j+1].Xc;
	  Wr = SolnBlk.W[i][j+1] +
	    (SolnBlk.phi[i][j+1]^SolnBlk.dWdx[i][j+1])*dX.x +
	    (SolnBlk.phi[i][j+1]^SolnBlk.dWdy[i][j+1])*dX.y;
	  if (SolnBlk.Grid.BCtypeS[i] == BC_REFLECTION) {
	    Wl = Reflect(Wr, SolnBlk.Grid.nfaceS(i, j+1)); 
	  } else if (SolnBlk.Grid.BCtypeS[i] == BC_FREE_SLIP_ISOTHERMAL) {
	    Wl = Free_Slip(Wr,SolnBlk.WoS[i], SolnBlk.Grid.nfaceS(i, j+1),FIXED_TEMPERATURE_WALL);
	  } else if (SolnBlk.Grid.BCtypeS[i] == BC_WALL_VISCOUS_ISOTHERMAL) {
	    Wl = No_Slip(Wr,SolnBlk.WoS[i], SolnBlk.Grid.nfaceS(i, j+1),FIXED_TEMPERATURE_WALL);
	  } else if (SolnBlk.Grid.BCtypeS[i] == BC_MOVING_WALL_ISOTHERMAL) {
	    Wl = Moving_Wall(Wr,SolnBlk.WoS[i], SolnBlk.Grid.nfaceS(i, j+1),
			     SolnBlk.Moving_wall_velocity,FIXED_TEMPERATURE_WALL);
	  } else if (SolnBlk.Grid.BCtypeS[i] == BC_WALL_VISCOUS_HEATFLUX) {
	    Wl = No_Slip(Wr,SolnBlk.WoS[i], SolnBlk.Grid.nfaceS(i, j+1),ADIABATIC_WALL);
	  } else if (SolnBlk.Grid.BCtypeS[i] == BC_MOVING_WALL_HEATFLUX) {
	    Wl = Moving_Wall(Wr,SolnBlk.WoS[i], SolnBlk.Grid.nfaceS(i, j+1),
			     SolnBlk.Moving_wall_velocity,ADIABATIC_WALL); 
	  } else if (SolnBlk.Grid.BCtypeS[i] == BC_1DFLAME_INFLOW){
	    Wl = BC_1DFlame_Inflow(Wr, 
				 SolnBlk.WoS[i], 
				 SolnBlk.W[i][SolnBlk.JCu],
				 SolnBlk.Grid.nfaceS(i, j+1));
	  } else if (SolnBlk.Grid.BCtypeS[i] == BC_2DFLAME_INFLOW){
	    Wl = BC_2DFlame_Inflow(Wr, 
				   SolnBlk.WoS[i], 
				   SolnBlk.Grid.nfaceS(i, j+1));
	  } else if (SolnBlk.Grid.BCtypeS[i] == BC_1DFLAME_OUTFLOW){
	    Wl = BC_1DFlame_Outflow(Wr, 
				    SolnBlk.WoS[i], 
				    SolnBlk.W[i][SolnBlk.JCu],
				    SolnBlk.Grid.nfaceS(i, j+1));
	  } else if(SolnBlk.Grid.BCtypeS[i] == BC_CHARACTERISTIC)  { 
	    Wl = BC_Characteristic_Pressure(Wr, 
					    SolnBlk.WoS[i], 
					    SolnBlk.Grid.nfaceS(i, j+1));
	  } else {
	    cerr<< "\n Wrong BC in  dUdt_Residual_Evaluation "<< SolnBlk.Grid.BCtypeS[i]; exit(1);
	  } 

	} else if (j == SolnBlk.JCu && 
		   (SolnBlk.Grid.BCtypeN[i] == BC_REFLECTION ||
		    SolnBlk.Grid.BCtypeN[i] == BC_CHARACTERISTIC ||
		    SolnBlk.Grid.BCtypeN[i] == BC_FREE_SLIP_ISOTHERMAL || 
		    SolnBlk.Grid.BCtypeN[i] == BC_WALL_VISCOUS_ISOTHERMAL ||
		    SolnBlk.Grid.BCtypeN[i] == BC_MOVING_WALL_ISOTHERMAL ||
		    SolnBlk.Grid.BCtypeN[i] == BC_WALL_VISCOUS_HEATFLUX ||
		    SolnBlk.Grid.BCtypeN[i] == BC_MOVING_WALL_HEATFLUX ||
		    SolnBlk.Grid.BCtypeN[i] == BC_1DFLAME_INFLOW ||
		    SolnBlk.Grid.BCtypeN[i] == BC_1DFLAME_OUTFLOW||
		    SolnBlk.Grid.BCtypeN[i] == BC_2DFLAME_INFLOW ||
		    SolnBlk.Grid.BCtypeN[i] == BC_2DFLAME_OUTFLOW )) {
	  
	  dX = SolnBlk.Grid.xfaceN(i, j)-SolnBlk.Grid.Cell[i][j].Xc;
	  Wl = SolnBlk.W[i][j] + 
	    (SolnBlk.phi[i][j]^SolnBlk.dWdx[i][j])*dX.x +
	    (SolnBlk.phi[i][j]^SolnBlk.dWdy[i][j])*dX.y;
	  if (SolnBlk.Grid.BCtypeN[i] == BC_REFLECTION) {
	    Wr = Reflect(Wl, SolnBlk.Grid.nfaceN(i, j));	
	  } else if (SolnBlk.Grid.BCtypeN[i] == BC_FREE_SLIP_ISOTHERMAL) {
	    Wr = Free_Slip(Wl, SolnBlk.WoN[i], SolnBlk.Grid.nfaceN(i, j),FIXED_TEMPERATURE_WALL);
	  } else if (SolnBlk.Grid.BCtypeN[i] == BC_WALL_VISCOUS_ISOTHERMAL) {
	    Wr = No_Slip(Wl, SolnBlk.WoN[i], SolnBlk.Grid.nfaceN(i, j),FIXED_TEMPERATURE_WALL);
	  } else if (SolnBlk.Grid.BCtypeN[i] == BC_MOVING_WALL_ISOTHERMAL) {
	    Wr = Moving_Wall(Wl,SolnBlk.WoN[i],  SolnBlk.Grid.nfaceN(i, j),
			     SolnBlk.Moving_wall_velocity,FIXED_TEMPERATURE_WALL);
	  } else if (SolnBlk.Grid.BCtypeN[i] == BC_WALL_VISCOUS_HEATFLUX) {
	    Wr = No_Slip(Wl, SolnBlk.WoN[i], SolnBlk.Grid.nfaceN(i, j),ADIABATIC_WALL );
	  } else if (SolnBlk.Grid.BCtypeN[i] == BC_MOVING_WALL_HEATFLUX) {
	    Wr = Moving_Wall(Wl,SolnBlk.WoN[i],  SolnBlk.Grid.nfaceN(i, j),
			     SolnBlk.Moving_wall_velocity,ADIABATIC_WALL ); 
	  } else if (SolnBlk.Grid.BCtypeN[i] == BC_1DFLAME_INFLOW){
	    Wr = BC_1DFlame_Inflow(Wl, 
				 SolnBlk.WoN[i], 
				 SolnBlk.W[i][SolnBlk.JCl],
				 SolnBlk.Grid.nfaceN(i, j));
	  } else if (SolnBlk.Grid.BCtypeN[i] == BC_1DFLAME_OUTFLOW){
	    Wr = BC_1DFlame_Outflow(Wl, 
				  SolnBlk.WoN[i], 
				  SolnBlk.W[i][SolnBlk.JCl],
				  SolnBlk.Grid.nfaceN(i, j));
	  } else if (SolnBlk.Grid.BCtypeN[i] == BC_2DFLAME_OUTFLOW){
	    Wr = BC_2DFlame_Outflow(Wl, 
				    SolnBlk.WoN[i], 
				    SolnBlk.Grid.nfaceN(i, j));	    
	  } else if(SolnBlk.Grid.BCtypeN[i] == BC_CHARACTERISTIC)  { 
	    Wr = BC_Characteristic_Pressure(Wl, 
					    SolnBlk.WoN[i], 
					    SolnBlk.Grid.nfaceN(i, j));
	  } else {
	    cerr<< "\n Wrong BC in  dUdt_Residual_Evaluation "<< SolnBlk.Grid.BCtypeN[i]; exit(1);
	  }  
	  
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
			   Input_Parameters.Preconditioning,SolnBlk.Flow_Type,delta_n);
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
	case FLUX_FUNCTION_AUSM_PLUS_UP :
	  Flux = FluxAUSMplus_up_n(Wl, Wr, SolnBlk.Grid.nfaceN(i, j));
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
	  case VISCOUS_RECONSTRUCTION_ARITHMETIC :
	    //Arithmetic Mean
	    Flux -= Viscous_FluxArithmetic_n(SolnBlk.U[i][j],
                                             SolnBlk.dWdx[i][j],
                                             SolnBlk.dWdy[i][j],
                                             SolnBlk.U[i][j+1],
                                             SolnBlk.dWdx[i][j+1],
                                             SolnBlk.dWdy[i][j+1],
	     			             SolnBlk.Grid.nfaceN(i,j));  // d/dz Fzv
	    break;
	  case VISCOUS_RECONSTRUCTION_MEAN_GRADIENT :
	    // Determine state at the interface.
	    W = HALF*(SolnBlk.W[i][j] + SolnBlk.W[i][j+1]);
	    // Determine gradient of the primitive variables at the interface.
	    dWdx = HALF*(SolnBlk.dWdx[i][j] + SolnBlk.dWdx[i][j+1]);
	    dWdy = HALF*(SolnBlk.dWdy[i][j] + SolnBlk.dWdy[i][j+1]);
	    // Compute the viscous flux.
	    Flux -= Viscous_Flux_n(W,dWdx,dWdy,
                                   SolnBlk.Axisymmetric,
				   SolnBlk.Grid.xfaceN(i,j),
				   SolnBlk.Grid.nfaceN(i,j));
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
	  case VISCOUS_RECONSTRUCTION_DIAMONDPATH_LEAST_SQUARES :
          case VISCOUS_RECONSTRUCTION_DIAMONDPATH_GREEN_GAUSS :
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

	//cout<<"\n j dUdt "<<i<<" "<<j<<" "<<k_residual<<" "<<SolnBlk.dUdt[i][j][k_residual]; 	         
	  
      } /* endfor */

      SolnBlk.dUdt[i][SolnBlk.JCu+1][k_residual].Vacuum();
      SolnBlk.dUdt[i][SolnBlk.JCl-1][k_residual].Vacuum();

    } /* endfor */


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
int Update_Solution_Multistage_Explicit(Flame2D_Quad_Block &SolnBlk,
                                        const int i_stage,
                                        Flame2D_Input_Parameters &Input_Parameters) {


  int k_residual;
  double omega, delta_n;
  int NUM_VAR_FLAME2D = SolnBlk.NumVar(); 

  // Memory for linear system solver. 
  DenseSystemLinEqs LinSys;
  DenseMatrix dSdU(NUM_VAR_FLAME2D-1,NUM_VAR_FLAME2D-1);       // Source Jacobian
  DenseMatrix dRdU(NUM_VAR_FLAME2D-1,NUM_VAR_FLAME2D-1);       // Residual Jacobian
  DenseMatrix Precon(NUM_VAR_FLAME2D-1,NUM_VAR_FLAME2D-1);     // Low Mach number preconditioner
  DenseMatrix Precon_Inv(NUM_VAR_FLAME2D-1,NUM_VAR_FLAME2D-1); // Inverse of low Mach number preconditioner
  
  /* Allocate memory for linear system solver. */
  LinSys.allocate(NUM_VAR_FLAME2D-1);

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
 
  for (int j  = SolnBlk.JCl ; j <= SolnBlk.JCu ; ++j ) {
    for (int i = SolnBlk.ICl ; i <= SolnBlk.ICu ; ++i ) {
      
      /******************************************************/
      /********** Fully Explicit ****************************/
      /******************************************************/
      if (Input_Parameters.Local_Time_Stepping == GLOBAL_TIME_STEPPING || 
	  Input_Parameters.Local_Time_Stepping == SCALAR_LOCAL_TIME_STEPPING ) {
	//Update
	SolnBlk.U[i][j] = SolnBlk.Uo[i][j] + omega*SolnBlk.dUdt[i][j][k_residual];
	//N-1 species
	SolnBlk.U[i][j][NUM_VAR_FLAME2D] = SolnBlk.U[i][j].rho*(ONE - SolnBlk.U[i][j].sum_species());	   
      }
      
      //Check for unphysical properties  
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
	    SolnBlk.U[i][j][NUM_VAR_FLAME2D] = SolnBlk.U[i][j].rho*(ONE - SolnBlk.U[i][j].sum_species());
	    
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
	  Input_Parameters.Local_Time_Stepping == MATRIX_LOCAL_TIME_STEPPING ||
	  Input_Parameters.Local_Time_Stepping == MATRIX_WITH_LOW_MACH_NUMBER_PRECONDITIONER) {
	
	dSdU.zero();  //RESET, MAY NOT BE NECESSARY AS DONE INSIDE FUNCTIONS ??
	dRdU.zero();

	/************ FORM LHS FOR SEMI-IMPLICIT METHOD *****************/
	if (Input_Parameters.Local_Time_Stepping == SEMI_IMPLICIT_LOCAL_TIME_STEPPING){

	  /* Semi implicit formulation: set up system of equations 
	     and include source Jacobian in the LHS matrix. */	  
	  SemiImplicitBlockJacobi(dSdU,SolnBlk,EXPLICIT,i, j);	  	  
	  LinSys.A.identity();
	     
	  //scalar multiplication
	  LinSys.A -= (omega*Input_Parameters.CFL_Number*SolnBlk.dt[i][j])*dSdU;	     
	 
	}
	
	// Spacing for preconditioner and viscous 
	if (SolnBlk.Flow_Type != FLOWTYPE_INVISCID && 
	    (Input_Parameters.Local_Time_Stepping == LOW_MACH_NUMBER_WEISS_SMITH_PRECONDITIONER || 
	     Input_Parameters.Local_Time_Stepping == SEMI_IMPLICIT_LOW_MACH_NUMBER_PRECONDITIONER ||
	     Input_Parameters.Local_Time_Stepping == MATRIX_WITH_LOW_MACH_NUMBER_PRECONDITIONER) ){	     
	  delta_n = min( TWO*(SolnBlk.Grid.Cell[i][j].A/
			      (SolnBlk.Grid.lfaceE(i, j)+SolnBlk.Grid.lfaceW(i, j))),
			 TWO*(SolnBlk.Grid.Cell[i][j].A/
			      (SolnBlk.Grid.lfaceN(i, j)+SolnBlk.Grid.lfaceS(i, j))));
	}
	
	/************ FORM LHS FOR LOW MACH NUMBER PRECONDITIONING ***************/
	if (Input_Parameters.Local_Time_Stepping == LOW_MACH_NUMBER_WEISS_SMITH_PRECONDITIONER ){	   

	  SolnBlk.Uo[i][j].Low_Mach_Number_Preconditioner(Precon,
							  SolnBlk.Flow_Type,
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
	  SemiImplicitBlockJacobi(dSdU,SolnBlk,EXPLICIT,i, j);
	  LinSys.A.identity();
	  
	  //Currently recalculating Pinv and P, should store between stages of time method ( as based on Uo, not U)
	  //Inverse of Preconditioner
	  SolnBlk.Uo[i][j].Low_Mach_Number_Preconditioner_Inverse(Precon_Inv,
								  SolnBlk.Flow_Type,
								  delta_n);
	  
	  //I - Pinv*(dt*(dS/dU))
	  LinSys.A -= Precon_Inv*((omega*Input_Parameters.CFL_Number*SolnBlk.dt[i][j])*dSdU);
	  
	  //Preconditioner 
	  SolnBlk.Uo[i][j].Low_Mach_Number_Preconditioner(Precon,
							  SolnBlk.Flow_Type,
							  delta_n);
	  
	  /******** LOW MACH NUMBER PRECONDITIONER CHECK ******************
	  cout<<"\n Cell "<<i<<" "<<j<<endl;
 	  cout<<"\n PRECON \n"<<Precon<<"\n -1 \n"<<Precon_Inv;
  	  cout<<"\n P*P-1 \n"<<Precon*Precon_Inv;
	  ******** LOW MACH NUMBER PRECONDITIONER CHECK ******************/

	  //LHS = P*( I - Pinv*(dt*(dS/dU)))	     
	  LinSys.A = Precon*LinSys.A;	    

 
	}//end SEMI_IMPLICIT_LOW_MACH_NUMBER_PRECONDITIONER

	/************ FORM LHS FOR MATRIX LOCAL TIME STEPPING BASED ON *****************
	 ************ POINT IMPLICIT BLOCK JACOBI PRECONDITIONER       *****************/
	//LHS = -dt*(dR/dU)
	if (Input_Parameters.Local_Time_Stepping == MATRIX_LOCAL_TIME_STEPPING || 
	    Input_Parameters.Local_Time_Stepping == MATRIX_WITH_LOW_MACH_NUMBER_PRECONDITIONER){

	  PointImplicitBlockJacobi(dRdU,SolnBlk,Input_Parameters,i, j);	 
	  LinSys.A = -SolnBlk.dt[i][j]*dRdU;
	  
	  //  preseve the nonsigularity of matrix A 
	  for(int index=0; index<NUM_VAR_FLAME2D-1; index++){
	    if(fabs(LinSys.A(index, index))<MICRO) LinSys.A(index,index) += ONE;
	  }
	}
	
	/******************************************************************
	 ******************* EVALUATE RHS *********************************
	 ******************************************************************/
	// omega*(dt*RHS)
	for (int k = 0; k < NUM_VAR_FLAME2D - 1; k++) {
	  LinSys.b(k) = omega*SolnBlk.dUdt[i][j][k_residual][k+1];
	}
	/* Solve system of equations using LU decomposition Gaussian elimination procedure. */
	LinSys.solve(LU_DECOMPOSITION);
	
	/* Update the conserved solution variables. */
	for (int k = 0; k < NUM_VAR_FLAME2D - 1; k++) {
	  SolnBlk.U[i][j][k+1] = SolnBlk.Uo[i][j][k+1] + LinSys.x(k);
	} 

	/*********************************************************
	     Using N-1 species equations so need to update 
             last species using:
	     c_n = 1- sum(1 to N-1) cs
	*********************************************************/
	SolnBlk.U[i][j][NUM_VAR_FLAME2D] = SolnBlk.U[i][j].rho*(ONE - SolnBlk.U[i][j].sum_species());

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
	      PointImplicitBlockJacobi(dRdU,SolnBlk,Input_Parameters,i, j);
	      //  LinSys.A = -(1.00*SolnBlk.dt[i][j])*dRdU;//Markthis
	      //Talk to Prof. Groth about ... 
	      // Same thing as before the follwoing formula works, but not "theoretically" correct ...
	      LinSys.A.identity();
	      LinSys.A -=(omega*Input_Parameters.CFL_Number*SolnBlk.dt[i][j])*dRdU; 
	      
	    }
	    
	    for (int k = 0; k < NUM_VAR_FLAME2D-1; k++) {
	      LinSys.b(k) = omega*SolnBlk.dUdt[i][j][k_residual][k+1];
	    } 
	    
	    LinSys.solve(LU_DECOMPOSITION);
	    
	    /* Update the conserved Solution */
	    for (int k = 0; k < NUM_VAR_FLAME2D-1; k++) {
	      SolnBlk.U[i][j][k+1] = SolnBlk.Uo[i][j][k+1] + LinSys.x(k);
	    } 
	    
	    SolnBlk.U[i][j][NUM_VAR_FLAME2D] = SolnBlk.U[i][j].rho*(ONE - SolnBlk.U[i][j].sum_species());
	    
	    if(SolnBlk.U[i][j].Unphysical_Properties_Check(n_residual_reduction))  break;
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
void Viscous_Calculations(Flame2D_Quad_Block &SolnBlk) {

  // FIXME - MRJC
  double mu, kappa, r, div_v;
  double Temperature, Rmix, Cp;
  Vector2D grad_T;
  Vector2D X;
  for (int j = SolnBlk.JCl-1; j <= SolnBlk.JCu+1; j++) {
    for (int i = SolnBlk.ICl-1; i <= SolnBlk.ICu+1; i++) {

      /* Determine the molecular transport properties. */      
      mu = SolnBlk.W[i][j].mu();
      kappa = SolnBlk.W[i][j].kappa();
      Temperature = SolnBlk.W[i][j].T();
      Rmix = SolnBlk.W[i][j].Rtot();
      Cp =  SolnBlk.W[i][j].Cp();
  
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
// 	SolnBlk.U[i][j].rhospec[k].diffusion_coef = mu/SolnBlk.W[i][j].Schmidt[k];
	/***************** mass fraction gradients *********************/
// 	SolnBlk.U[i][j].rhospec[k].gradc.x = SolnBlk.U[i][j].rho * SolnBlk.dWdx[i][j].spec[k].c;
// 	SolnBlk.U[i][j].rhospec[k].gradc.y = SolnBlk.U[i][j].rho * SolnBlk.dWdy[i][j].spec[k].c;
      }
      X = SolnBlk.Grid.Cell[i][j].Xc;
      /***************** Molecular (Laminar) Stresses ******************/
      SolnBlk.W[i][j].Laminar_Stress(SolnBlk.dWdx[i][j], SolnBlk.dWdy[i][j],
				     SolnBlk.Axisymmetric, X);
      SolnBlk.U[i][j].tau =  SolnBlk.W[i][j].tau;

      /************* Molecular (Laminar) Heat flux Vector **************/
      /****************** Thermal Conduction ***************************
         q = - kappa * grad(T)                                         */
      SolnBlk.U[i][j].qflux = - kappa*grad_T;
      /****************** Thermal Diffusion ****************************/
      // q -= rho * sum ( hs * Ds *gradcs)  
      SolnBlk.U[i][j].qflux -= SolnBlk.U[i][j].rho*SolnBlk.U[i][j].thermal_diffusion(Temperature);  
      
      /****** Set Primitive *******/
      SolnBlk.W[i][j] = W(SolnBlk.U[i][j]);

    }     
  }
}




/**********************************************************************
 * Radiation_Source_Eval                                              *
 *                                                                    *
 * Optically thin radiation source term evaluation.  The radiation    *
 * source term is the divergence of the radiative flux vector.        *
 * Here it is evaluated using the optically thin approximation.       *
 *                                                                    *
 **********************************************************************/
void Radiation_Source_Eval( Flame2D_Quad_Block &SolnBlk,
			    Flame2D_Input_Parameters &Input_Parameters ) {

  //
  // declares
  //
  double Temperature, Pressure;
  double xCO, xH2O, xCO2, xO2;
  static const double fsoot = ZERO;  // no soot for now

  // if SNBCK not allocated, allocate/setup a new object
  if (Flame2D_Quad_Block::PlanckMean_data==NULL) {
    Flame2D_Quad_Block::PlanckMean_data = new PlanckMean(Input_Parameters.CFFC_Path);
  } // endif


  //
  // loop over the block
  //
  for (int j = SolnBlk.JCl-1; j <= SolnBlk.JCu+1; j++) {
    for (int i = SolnBlk.ICl-1; i <= SolnBlk.ICu+1; i++) {

      // get the temperature and pressure at this point
      Temperature = SolnBlk.W[i][j].T();
      Pressure = SolnBlk.W[i][j].p;

      // first, get radiating species concentrations
      SolnBlk.W[i][j].MoleFracOfRadSpec( xCO,  xH2O, xCO2, xO2 );

      // compute the source term using optically thin approx
      // using the SNBCK model
      SolnBlk.Srad[i][j] = 
	Flame2D_Quad_Block::PlanckMean_data->RadSourceOptThin( Pressure/PRESSURE_STDATM, //[atm]
							      Temperature,              //[K]
							      xCO,
							      xH2O,
							      xCO2,
							      xO2,
							      fsoot );

    } // endfor i
  } // endfor j

} // end Radiation_Source_Eval()
