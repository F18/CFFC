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

  const int NUM_VAR_FLAME2D = SolnBlk.NumVar(); 

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
	    buffer[buffer_size] = SolnBlk.U[i][j][k+1];
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
	    SolnBlk.U[i][j][k+1]= buffer[buffer_size];
	    buffer_size = buffer_size + 1;
	  }
	  SolnBlk.W[i][j].setU(SolnBlk.U[i][j]);
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

  const int NUM_VAR_FLAME2D = SolnBlk.NumVar();

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
	    buffer[buffer_size] = SolnBlk.U[i][j][k+1]; 
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
	    SolnBlk.U[i][j][k+1]= buffer[buffer_size];
	    buffer_size = buffer_size + 1;
	  }
	  SolnBlk.W[i][j].setU(SolnBlk.U[i][j]);
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

	SolnBlk_Fine.U[2*(i-i_min)+SolnBlk_Fine.ICl  ]
	  [2*(j-j_min)+SolnBlk_Fine.JCl  ] 
	  = SolnBlk_Original.U[i][j];
	SolnBlk_Fine.W[2*(i-i_min)+SolnBlk_Fine.ICl  ]
	  [2*(j-j_min)+SolnBlk_Fine.JCl  ].setU(SolnBlk_Original.U[i][j]);

	SolnBlk_Fine.U[2*(i-i_min)+SolnBlk_Fine.ICl+1]
	  [2*(j-j_min)+SolnBlk_Fine.JCl  ] 
	  = SolnBlk_Original.U[i][j];
	SolnBlk_Fine.W[2*(i-i_min)+SolnBlk_Fine.ICl+1]
	  [2*(j-j_min)+SolnBlk_Fine.JCl  ].setU(SolnBlk_Original.U[i][j]);

	SolnBlk_Fine.U[2*(i-i_min)+SolnBlk_Fine.ICl  ]
	  [2*(j-j_min)+SolnBlk_Fine.JCl+1]
	  = SolnBlk_Original.U[i][j];
	SolnBlk_Fine.W[2*(i-i_min)+SolnBlk_Fine.ICl  ]
	  [2*(j-j_min)+SolnBlk_Fine.JCl+1].setU(SolnBlk_Original.U[i][j]);

	SolnBlk_Fine.U[2*(i-i_min)+SolnBlk_Fine.ICl+1]
	  [2*(j-j_min)+SolnBlk_Fine.JCl+1]
	  = SolnBlk_Original.U[i][j];
	SolnBlk_Fine.W[2*(i-i_min)+SolnBlk_Fine.ICl+1]
	  [2*(j-j_min)+SolnBlk_Fine.JCl+1].setU(SolnBlk_Original.U[i][j]);
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

  const int NUM_VAR_FLAME2D = SolnBlk_Coarse.NumVar();
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
	for ( int k=1; k<=NUM_VAR_FLAME2D; ++k)
	  SolnBlk_Coarse.U[i_coarse][j_coarse][k] = ( (SolnBlk_Original_SW.Grid.Cell[i  ][j  ].A*
						       SolnBlk_Original_SW.U[i  ][j  ][k] +
						       SolnBlk_Original_SW.Grid.Cell[i+1][j  ].A*
						       SolnBlk_Original_SW.U[i+1][j  ][k] + 
						       SolnBlk_Original_SW.Grid.Cell[i  ][j+1].A*
						       SolnBlk_Original_SW.U[i  ][j+1][k] +
						       SolnBlk_Original_SW.Grid.Cell[i+1][j+1].A*
						       SolnBlk_Original_SW.U[i+1][j+1][k]) /
						      (SolnBlk_Original_SW.Grid.Cell[i  ][j  ].A +
						       SolnBlk_Original_SW.Grid.Cell[i+1][j  ].A +
						       SolnBlk_Original_SW.Grid.Cell[i  ][j+1].A +
						       SolnBlk_Original_SW.Grid.Cell[i+1][j+1].A) );
	SolnBlk_Coarse.W[i_coarse][j_coarse].setU(SolnBlk_Coarse.U[i_coarse][j_coarse]);
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
	for ( int k=1; k<=NUM_VAR_FLAME2D; ++k)
	  SolnBlk_Coarse.U[i_coarse][j_coarse][k] = ( (SolnBlk_Original_SE.Grid.Cell[i  ][j  ].A*
						       SolnBlk_Original_SE.U[i  ][j  ][k] +
						       SolnBlk_Original_SE.Grid.Cell[i+1][j  ].A*
						       SolnBlk_Original_SE.U[i+1][j  ][k] + 
						       SolnBlk_Original_SE.Grid.Cell[i  ][j+1].A*
						       SolnBlk_Original_SE.U[i  ][j+1][k] +
						       SolnBlk_Original_SE.Grid.Cell[i+1][j+1].A*
						       SolnBlk_Original_SE.U[i+1][j+1][k]) /
						      (SolnBlk_Original_SE.Grid.Cell[i  ][j  ].A +
						       SolnBlk_Original_SE.Grid.Cell[i+1][j  ].A +
						       SolnBlk_Original_SE.Grid.Cell[i  ][j+1].A +
						       SolnBlk_Original_SE.Grid.Cell[i+1][j+1].A) );
	SolnBlk_Coarse.W[i_coarse][j_coarse].setU(SolnBlk_Coarse.U[i_coarse][j_coarse]);
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
	for ( int k=1; k<=NUM_VAR_FLAME2D; ++k)
	  SolnBlk_Coarse.U[i_coarse][j_coarse][k] = ( (SolnBlk_Original_NW.Grid.Cell[i  ][j  ].A*
						       SolnBlk_Original_NW.U[i  ][j  ][k] +
						       SolnBlk_Original_NW.Grid.Cell[i+1][j  ].A*
						       SolnBlk_Original_NW.U[i+1][j  ][k] + 
						       SolnBlk_Original_NW.Grid.Cell[i  ][j+1].A*
						       SolnBlk_Original_NW.U[i  ][j+1][k] +
						       SolnBlk_Original_NW.Grid.Cell[i+1][j+1].A*
						       SolnBlk_Original_NW.U[i+1][j+1][k]) /
						      (SolnBlk_Original_NW.Grid.Cell[i  ][j  ].A +
						       SolnBlk_Original_NW.Grid.Cell[i+1][j  ].A +
						       SolnBlk_Original_NW.Grid.Cell[i  ][j+1].A +
						       SolnBlk_Original_NW.Grid.Cell[i+1][j+1].A) );
	SolnBlk_Coarse.W[i_coarse][j_coarse].setU(SolnBlk_Coarse.U[i_coarse][j_coarse]);
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
	for ( int k=1; k<=NUM_VAR_FLAME2D; ++k)
	  SolnBlk_Coarse.U[i_coarse][j_coarse][k] = ( (SolnBlk_Original_NE.Grid.Cell[i  ][j  ].A*
						       SolnBlk_Original_NE.U[i  ][j  ][k] +
						       SolnBlk_Original_NE.Grid.Cell[i+1][j  ].A*
						       SolnBlk_Original_NE.U[i+1][j  ][k] + 
						       SolnBlk_Original_NE.Grid.Cell[i  ][j+1].A*
						       SolnBlk_Original_NE.U[i  ][j+1][k] +
						       SolnBlk_Original_NE.Grid.Cell[i+1][j+1].A*
						       SolnBlk_Original_NE.U[i+1][j+1][k]) /
						      (SolnBlk_Original_NE.Grid.Cell[i  ][j  ].A +
						       SolnBlk_Original_NE.Grid.Cell[i+1][j  ].A +
						       SolnBlk_Original_NE.Grid.Cell[i  ][j+1].A +
						       SolnBlk_Original_NE.Grid.Cell[i+1][j+1].A) );
	SolnBlk_Coarse.W[i_coarse][j_coarse].setU(SolnBlk_Coarse.U[i_coarse][j_coarse]);
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
   
  Flame2D_pState W_node;
  Vector2D qflux;
  Tensor2D tau;
  Vector2D Vcorr;

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
    for(int i =0; i<Flame2D_pState::NumSpecies(); i++){
      Out_File <<"\"c_"<<Flame2D_pState::speciesName(i)<<"\" \\ \n";
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
	    << "\"e_s\" \\ \n";
    //viscous quantities
    Out_File << "\"qflux_x\" \\ \n"  
	     << "\"qflux_y\" \\  \n"   
	     << "\"Tau_xx\" \\  \n"  //rr -axisymmetric
	     << "\"Tau_xy\" \\  \n"  //rz
	     << "\"Tau_yy\" \\  \n"  //zz
	     << "\"Tau_zz\" \\  \n"
	     << "\"Vcorr_x\" \\  \n"
	     << "\"Vcorr_y\" \\  \n";
    // Zone details
    Out_File << "ZONE T =  \"Block Number = " << Block_Number
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
      //get nodal values
      W_node = SolnBlk.Wn(i, j);
      // update related transport properties -> kappa and D_i
      W_node.updateTransport();
      // Compute viscous quantities
      SolnBlk.Viscous_Quantities_n( i, j, qflux, tau, Vcorr );
      //coordinates 
      Out_File << " " << SolnBlk.Grid.Node[i][j].X<<endl;
      //cell properties
      Out_File << W_node;
      //T,M,H,s, and all the rest of the calcuatated parameters
      Out_File.setf(ios::scientific);
      Out_File << " " << W_node.T()<< " " << W_node.Rtot()
	       << " " << W_node.vabs()/W_node.a() 
	       << " " << W_node.mu() <<" "<< W_node.kappa()
	       << " " << W_node.Pr()
	       << " " << W_node.H() 
	       << " " << W_node.h() 
	       << " " << W_node.hs()
	       << " " << ((const Flame2D_pState&)W_node).E() 
	       << " " << W_node.e() 
	       << " " << W_node.es();
      // viscous terms
      Out_File << " " << qflux
	       << " " << tau
	       << " " << Vcorr;
      Out_File.unsetf(ios::scientific);
      Out_File << endl;
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
  double rho, p;
  Vector2D qflux;
  Tensor2D tau;
  Vector2D Vcorr;
  Flame2D_State omega;

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
    for(int i =0; i<Flame2D_pState::NumSpecies(); i++){
      Out_File <<"\"c"<<Flame2D_pState::speciesName(i)<<"\" \\ \n";
    }
    // momentum
    Out_File << "\"rho*u\" \\ \n"
	     << "\"rho*v\" \\ \n";
    // checks
    Out_File << "\"Sum Yi\" \\ \n"
      	     << "\"Sum omega\" \\ \n"
	     << "\"Sum Vk_x\" \\ \n"
	     << "\"Sum Vk_y\" \\ \n";
    //Calculated values
    Out_File << "\"T\" \\ \n"
	     << "\"M\" \\ \n"
	     << "\"R\" \\ \n"
	     << "\"viscosity\" \\ \n"
	     << "\"thermal conduct\" \\ \n"
	     << "\"Prandtl\" \\ \n"
	     << "\"rho*H\"  \\ \n"  
	     <<"\"h\" \\ \n"
	     <<"\"h_s\" \\ \n"
	     <<"\"rho*E\" \\ \n"
	     << "\"e\" \\  \n" 
	     << "\"e_s\" \\ \n";
    for(int i =0; i<Flame2D_pState::NumSpecies(); i++){
      Out_File <<"\"D_"<<Flame2D_pState::speciesName(i)<<"\" \\ \n";
    }
    //viscous quantities
    Out_File << "\"qflux_x\" \\ \n"  
	     << "\"qflux_y\" \\  \n"   
	     << "\"Tau_xx\" \\  \n"  //rr -axisymmetric
	     << "\"Tau_xy\" \\  \n"  //rz
	     << "\"Tau_yy\" \\  \n"  //zz
	     << "\"Tau_zz\" \\  \n"
	     << "\"Vcorr_x\" \\ \n"
	     << "\"Vcorr_y\" \\ \n";
    // Reaction Rates
    for(int i =0; i<Flame2D_pState::NumSpecies(); i++){
      Out_File <<"\"omega_c"<<Flame2D_pState::speciesName(i)<<"\" \\ \n";
    }
    // Gradients
    Out_File << "\"dWdx_rho\" \\ \n"
	     << "\"dWdx_u\" \\ \n"
	     << "\"dWdx_v\" \\ \n"
	     << "\"dWdx_p\" \\ \n";
    for(int i =0; i<Flame2D_pState::NumSpecies(); i++){
      Out_File <<"\"dWdx_c"<<Flame2D_pState::speciesName(i)<<"\" \\ \n";
    }
    Out_File << "\"dWdy_rho\" \\ \n"
	     << "\"dWdy_u\" \\ \n"
	     << "\"dWdy_v\" \\ \n"
	     << "\"dWdy_p\" \\ \n";
    for(int i =0; i<Flame2D_pState::NumSpecies(); i++){
      Out_File <<"\"dWdy_c"<<Flame2D_pState::speciesName(i)<<"\" \\ \n";
    }
    Out_File << "\"phi_rho\" \\ \n"
	     << "\"phi_u\" \\ \n"
	     << "\"phi_v\" \\ \n"
	     << "\"phi_p\" \\ \n";
    for(int i =0; i<Flame2D_pState::NumSpecies(); i++){
      Out_File <<"\"phi_c"<<Flame2D_pState::speciesName(i)<<"\" \\ \n";
    }
    // dUdt
    Out_File << "\"dUdt_rho\" \\ \n"
	     << "\"dUdt_rhou\" \\ \n"
	     << "\"dUdt_rhov\" \\ \n"
	     << "\"dUdt_E\" \\ \n";
    for(int i =0; i<Flame2D_pState::NumSpecies(); i++){
      Out_File <<"\"dUdt_rhoc"<<Flame2D_pState::speciesName(i)<<"\" \\ \n";
    }

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

      // update related transport properties -> kappa and D_i
      SolnBlk.W[i][j].updateTransport();
      // compute viscous terms
      SolnBlk.W[i][j].Viscous_Quantities( SolnBlk.dWdx[i][j], SolnBlk.dWdy[i][j], 
					  SolnBlk.Axisymmetric, 
					  SolnBlk.Grid.Cell[i][j].Xc, 
					  qflux, tau, Vcorr );
      // compute reaction rates
      omega.Vacuum();
      SolnBlk.W[i][j].Sw( omega );
      rho = ((const Flame2D_pState&)SolnBlk.W[i][j]).rho();
      p = ((const Flame2D_pState&)SolnBlk.W[i][j]).p();
      // grid location
      Out_File << " "  << SolnBlk.Grid.Cell[i][j].Xc;
      // cell value
      Out_File << SolnBlk.W[i][j];
      // momentum
      Out_File << " " << SolnBlk.U[i][j].rhovx()
	       << " " << SolnBlk.U[i][j].rhovy();
      // checks
      Out_File << " " << SolnBlk.W[i][j].SpecSum()
	       << " " << SolnBlk.W[i][j].OmegaSum()
	       << " " << SolnBlk.W[i][j].DiffSum(SolnBlk.dWdx[i][j].c(), 
						 SolnBlk.dWdy[i][j].c()) + Vcorr;
      // Calculated Values
      Out_File.setf(ios::scientific);
      Out_File << " " << SolnBlk.W[i][j].T()
	       << " " << SolnBlk.W[i][j].vabs()/SolnBlk.W[i][j].a() 
	       << " " << SolnBlk.W[i][j].Rtot()
	       << " " << SolnBlk.W[i][j].mu()
	       << " " << SolnBlk.W[i][j].kappa()
	       << " " << SolnBlk.W[i][j].Pr()
	       << " " << SolnBlk.W[i][j].H() 
	       << " " << SolnBlk.W[i][j].h() 
	       << " " << SolnBlk.W[i][j].hs()
	       << " " << ((const Flame2D_pState&)SolnBlk.W[i][j]).E() 
	       << " " << SolnBlk.W[i][j].e()
	       << " " << SolnBlk.W[i][j].es();
      for(int k=0; k<Flame2D_pState::NumSpecies(); k++){
	Out_File <<" "<<SolnBlk.W[i][j].Diffusion_coef(k);
      }
      // viscous terms
      Out_File << " " << qflux
	       << " " << tau
	       << " " << Vcorr;
      // reaction rates
      for(int k=0; k<Flame2D_pState::NumSpecies(); k++){
	Out_File <<" "<<omega.rhoc(k) / rho;
      }

      Out_File.unsetf(ios::scientific);
      // Gradients
      Out_File << SolnBlk.dWdx[i][j]
	       << SolnBlk.dWdy[i][j]
	       << SolnBlk.phi[i][j];
      // dUdt
      Out_File << SolnBlk.dUdt[i][j][0];
      Out_File << endl;
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
	 const Flame2D_pState *Wo, 
	 const Flame2D_Input_Parameters &Input_Parameters) {


  //
  // Assign the initial data for the IVP of interest.
  //

  //--------------------------------------------------
  // UNIFORM.
  //--------------------------------------------------
  if (i_ICtype == IC_CONSTANT || i_ICtype == IC_UNIFORM) {
	
    // Set the solution state to the initial state Wo[0].
    for (int j  = SolnBlk.JCl-SolnBlk.Nghost ; j <= SolnBlk.JCu+SolnBlk.Nghost ; ++j ) {
      for ( int i = SolnBlk.ICl-SolnBlk.Nghost ; i <= SolnBlk.ICu+SolnBlk.Nghost ; ++i ) {
	SolnBlk.W[i][j] = Wo[0];
	SolnBlk.W[i][j].getU(SolnBlk.U[i][j]);
      } 
    }
    
    //--------------------------------------------------
    // SOD Problem in X-dir
    //--------------------------------------------------
  } else if (i_ICtype == IC_SOD_XDIR) {
 
    Flame2D_State Wl, Wr;
 
    // Set initial data for Sod IVP in x-direction.
    Wl = Wo[0];
    Wr = Wo[0];
    Wr.rho() = Wo[0].rho()/EIGHT;
    Wr.p() = Wo[0].p()/TEN;
    
    for (int j  = SolnBlk.JCl-SolnBlk.Nghost ; j <= SolnBlk.JCu+SolnBlk.Nghost ; ++j ) {
      for ( int i = SolnBlk.ICl-SolnBlk.Nghost ; i <= SolnBlk.ICu+SolnBlk.Nghost ; ++i ) {
	if (SolnBlk.Grid.Cell[i][j].Xc.x <= ZERO) {
	  SolnBlk.W[i][j] = Wl;
	} else {
	  SolnBlk.W[i][j] = Wr;
	} /* end if */
	SolnBlk.W[i][j].getU(SolnBlk.U[i][j]);
      } /* endfor */
    } /* endfor */
    

      //--------------------------------------------------
      // Shock Box
      //--------------------------------------------------
  } else if (i_ICtype == IC_SHOCK_BOX) {

    Flame2D_State Wl, Wr;

    // Set initial data for Aki shock-box IVP.
    Wl = Wo[0];
    Wl.vx() = ZERO; Wl.vy()=ZERO;
    Wr = Wo[0];
    Wr.rho() = Wr.rho()*FOUR;
    Wr.vx() = ZERO; Wr.vy() = ZERO;
    Wr.p() = Wr.p()*FOUR;
    for (int j  = SolnBlk.JCl-SolnBlk.Nghost ; j <= SolnBlk.JCu+SolnBlk.Nghost ; ++j ) {
      for (int i = SolnBlk.ICl-SolnBlk.Nghost ; i <= SolnBlk.ICu+SolnBlk.Nghost ; ++i ) {
	if (SolnBlk.Grid.Cell[i][j].Xc.x <= ZERO &&
	    SolnBlk.Grid.Cell[i][j].Xc.y <= ZERO) {
	  SolnBlk.W[i][j] = Wl;
	} else {
	  SolnBlk.W[i][j] = Wr;	     
	} /* end if */
	SolnBlk.W[i][j].getU(SolnBlk.U[i][j]);
      } /* endfor */
    } /* endfor */
    

      //--------------------------------------------------
      // Viscous Couette flow (with pressure gradient)
      //--------------------------------------------------
  } else if (i_ICtype == IC_VISCOUS_COUETTE) {

    //Couette flow with pressure gradient
    // -635.00  -423.33  -211.67  0.00  211.67  423.33  635.00

    double pl, pr;
    Vector2D dX;

    //total pressure change
    double delta_pres = Input_Parameters.Pressure_Gradient; // 635.54; 
    //grid spacing 
    dX.x = fabs(SolnBlk.Grid.Cell[SolnBlk.ICl][0].Xc.x - SolnBlk.Grid.Cell[SolnBlk.ICl-1][0].Xc.x);
    //block size
    dX.y = fabs(SolnBlk.Grid.Cell[SolnBlk.ICl][0].Xc.x - SolnBlk.Grid.Cell[SolnBlk.ICu][0].Xc.x) + dX.x; 
    //initial pressure for this block
    pr = Wo[0].p() - delta_pres *(( SolnBlk.Grid.Cell[SolnBlk.ICl][0].Xc.x - dX.x/2.0 + 0.1)/0.2);
    //pressure drop per block
    pl = delta_pres * dX.y/0.2;	  
    
    for (int j = SolnBlk.JCl-SolnBlk.Nghost; j <= SolnBlk.JCu+SolnBlk.Nghost; j++) {
      for (int i = SolnBlk.ICl-SolnBlk.Nghost; i <= SolnBlk.ICu+SolnBlk.Nghost; i++) {
	SolnBlk.W[i][j] = Wo[0];
	SolnBlk.W[i][j].updateViscosity();
	
	//pressure profiles
	if( i == SolnBlk.ICl-SolnBlk.Nghost ){
	  SolnBlk.W[i][j].setPressure( pr + pl*dX.x/dX.y*1.5 );
	} else if (i == SolnBlk.ICl-SolnBlk.Nghost+1 ){
	  SolnBlk.W[i][j].setPressure( pr + pl*dX.x/dX.y*0.5 );
	} else if( i == SolnBlk.ICu+SolnBlk.Nghost ){
	  SolnBlk.W[i][j].setPressure( pr - pl - pl*dX.x/dX.y*1.5 );
	} else if(i == SolnBlk.ICu+SolnBlk.Nghost-1 ){
	  SolnBlk.W[i][j].setPressure( pr - pl - pl*dX.x/dX.y*0.5 );
	} else {	   
	  SolnBlk.W[i][j].setPressure( pr - 
				       double(i-SolnBlk.ICl)*pl/double(SolnBlk.ICu-SolnBlk.ICl+1) - 
				       pl*dX.x/dX.y*0.5 );
	}
	
	//velocity profiles
	if (j >= SolnBlk.JCl && j <= SolnBlk.JCu){
	  SolnBlk.W[i][j].updateViscosity();
	  SolnBlk.W[i][j].setVelocityX( (HALF/SolnBlk.W[i][j].mu())*(-delta_pres/0.2)*
					(pow(SolnBlk.Grid.Cell[i][j].Xc.y,TWO) - (0.001*0.001/4.0)) +
					SolnBlk.Moving_wall_velocity*(SolnBlk.Grid.Cell[i][j].Xc.y/0.001 + 0.5) );
	}
	
	//update conserved 
	SolnBlk.W[i][j].getU(SolnBlk.U[i][j]);
      }
    }
    
    //--------------------------------------------------
    // Viscous Blasius Boundary Layer
    //--------------------------------------------------
  } else if (i_ICtype == IC_VISCOUS_FLAT_PLATE) {

    double eta,f,fp,fpp;

    // Set the initial data to the Blasius (exact) solution for the
    // laminar flow over a flat plate in the x-direction.
    for (int j = SolnBlk.JCl-SolnBlk.Nghost; j <= SolnBlk.JCu+SolnBlk.Nghost; j++) {
      for (int i = SolnBlk.ICl-SolnBlk.Nghost; i <= SolnBlk.ICu+SolnBlk.Nghost; i++) {

	if (SolnBlk.Grid.Cell[i][j].Xc.y >= ZERO) 
	  SolnBlk.W[i][j].FlatPlate(Wo[0], SolnBlk.Grid.Cell[i][j].Xc, eta,f,fp,fpp);
	else 
	  SolnBlk.W[i][j].FlatPlate(Wo[0], Vector2D(SolnBlk.Grid.Cell[i][j].Xc.x,
						    -SolnBlk.Grid.Cell[i][j].Xc.y),
				    eta, f, fp, fpp);
	SolnBlk.W[i][j].getU(SolnBlk.U[i][j]);
      }
    }
    
    //--------------------------------------------------
    // Driven Cavity Flow
    //--------------------------------------------------
  } else if (i_ICtype == IC_VISCOUS_DRIVEN_CAVITY_FLOW) {

    for (int j = SolnBlk.JCl-SolnBlk.Nghost; j <= SolnBlk.JCu+SolnBlk.Nghost; j++) {
      for (int i = SolnBlk.ICl-SolnBlk.Nghost; i <= SolnBlk.ICu+SolnBlk.Nghost; i++) {
	SolnBlk.W[i][j] = Wo[0];
	SolnBlk.W[i][j].getU(SolnBlk.U[i][j]);
      }
    }
    

    //--------------------------------------------------
    // 1D premixed flame
    //--------------------------------------------------
  } else if (i_ICtype == IC_1DFLAME) {
    // Set initial data with left stoichiometric premixed gas
    // and the right hand side being the associated combusted
    // products  
    Flame2D_pState Wl = Wo[0]; // unburnt
    Flame2D_pState Wr = Wo[0]; // burnt

    // get equilibrium burned composition
    Wr.equilibrate_HP();
 
    // set laminar flame speed
    Wl.setVelocity( Input_Parameters.flame_speed, 0.0 );

    // compute flame jump conditions for T and v
    Wr.FlameJumpLowMach(Wl);

    // fix constant pressure
    Wr.setState_TPY( Wr.T(), ((const Flame2D_pState&)Wl).p(), ((const Flame2D_pState&)Wr).c() );

    // Set Initial condtions on 1D grid
    for (int j  = SolnBlk.JCl-SolnBlk.Nghost ; j <= SolnBlk.JCu+SolnBlk.Nghost ; ++j ) {
      for ( int i = SolnBlk.ICl-SolnBlk.Nghost ; i <= SolnBlk.ICu+SolnBlk.Nghost ; ++i ) {
	if (SolnBlk.Grid.Cell[i][j].Xc.x <= ZERO){ //spatial relation, grid independent 
	  SolnBlk.W[i][j] = Wl;  
	} else {
	  SolnBlk.W[i][j] = Wr;	     
	} /* end if */
	SolnBlk.W[i][j].getU(SolnBlk.U[i][j]);
      } 
    } 

    //--------------------------------------------------
    // 2D laminar coflow diffusion flame
    //--------------------------------------------------
  } else if (i_ICtype == IC_CORE_FLAME) {
    // The Core Flame Initial Conditions sets up a laminar diffusion 
    // flame in the axisymmetric coordinate system, with the axis of symmetry
    // being the left or west boundary y-axis.
    // The condition's are based on those outlined by Mohammed et al. as well
    // as Day and Bell.  
    double fuel_spacing = 0.002;      //m
    double fuel_velocity = 0.70;      //m/s  //0.70
    double fuel_temp_inlet = 298.0;   //K 
    double tube_thickness =  0.00038; //m delta	
    
    double air_spacing = 0.030;       //m   //0.025
    double air_velocity = 0.35;       //m/s  0.35
    double air_temp_inlet = 298.0;    //K
    double ignition_temp = 1300.0;    //K
    
    // Temporaries
    Flame2D_pState Wf = Wo[0]; // Fuel
    Flame2D_pState Wa = Wo[0]; // Air
    Flame2D_pState Wi = Wo[0]; // ignitor
    const int nsp = Wo[0].NumSpecies();
    double* y = new double[nsp];
    double Press = Wo[0].p();
   
    //fuel 65% FUEL & 35% N2
    for (int i=0; i<nsp; i++) y[i] = 0.0;
    y[Flame2D_pState::speciesIndex("CH4")] = 0.5149;
    y[Flame2D_pState::speciesIndex("N2")]  = 0.4851;
    Wf.setState_TPY( fuel_temp_inlet, Press, y );
    Wf.setVelocity(0.0, 0.0);
    
    //air 21% O2 & 79% N2
    for (int i=0; i<nsp; i++) y[i] = 0.0;
    y[Flame2D_pState::speciesIndex("O2")] = 0.232;
    y[Flame2D_pState::speciesIndex("N2")] = 0.768;
    Wa.setState_TPY( air_temp_inlet, Press, y );
    Wa.setVelocity(0.0, 0.0);

    // ignitor is equilibrium combustion producs
    Wi.equilibrate_HP();

    //
    // Set the inital conditions everywhere
    //
    for (int j  = SolnBlk.JCl-SolnBlk.Nghost ; j <= SolnBlk.JCu+SolnBlk.Nghost ; ++j ) {
      for ( int i = SolnBlk.ICl-SolnBlk.Nghost ; i <= SolnBlk.ICu+SolnBlk.Nghost ; ++i ) {

	//region for injected fuel parabolic profile
	if (SolnBlk.Grid.Cell[i][j].Xc.x <= fuel_spacing ){ 
	  if (SolnBlk.Grid.Cell[i][j].Xc.y <0.006 ){
	    SolnBlk.W[i][j] = Wf;
	    SolnBlk.W[i][j].setVelocityY( (ONE - pow((SolnBlk.Grid.Cell[i][j].Xc.x/fuel_spacing),TWO))*fuel_velocity );
	  } else {
	    SolnBlk.W[i][j] = Wa;
	    SolnBlk.W[i][j].setVelocityY( (ONE - pow((SolnBlk.Grid.Cell[i][j].Xc.x/fuel_spacing),TWO))*fuel_velocity );
	  }
	//region for injected air parabolic profile
	} else if (SolnBlk.Grid.Cell[i][j].Xc.x > fuel_spacing+tube_thickness && 
		   SolnBlk.Grid.Cell[i][j].Xc.x <= 0.05*air_spacing + fuel_spacing+tube_thickness ){		    
	  SolnBlk.W[i][j] = Wa;
	  SolnBlk.W[i][j].setVelocityY( (ONE - pow( (SolnBlk.Grid.Cell[i][j].Xc.x - 
						     fuel_spacing - tube_thickness -
						     0.05*air_spacing) / 
						    (0.05*air_spacing), TWO )) * air_velocity );
	//region for injected air
	} else if (SolnBlk.Grid.Cell[i][j].Xc.x > fuel_spacing+tube_thickness && 
		   SolnBlk.Grid.Cell[i][j].Xc.x <= air_spacing ){	
	  SolnBlk.W[i][j] = Wa;
	  SolnBlk.W[i][j].setVelocityY( air_velocity );
	//region for quiesent air
	} else {
	  SolnBlk.W[i][j] = Wa;	    	  	    	  
	} 
	
	//IGNITOR across fuel and air inlets   //0.006  & 0.003
	if( SolnBlk.Grid.Cell[i][j].Xc.y < 0.006 && SolnBlk.Grid.Cell[i][j].Xc.y > 0.003){   	   
	  if ( SolnBlk.Grid.Cell[i][j].Xc.x <= fuel_spacing && SolnBlk.Grid.Cell[i][j].Xc.y <0.011 && 
	       SolnBlk.Grid.Cell[i][j].Xc.x > fuel_spacing*0.25){ 
	    SolnBlk.W[i][j] = Wi;
	  } else if (SolnBlk.Grid.Cell[i][j].Xc.x <= fuel_spacing){
	    SolnBlk.W[i][j] = Wi;
	  } else if (SolnBlk.Grid.Cell[i][j].Xc.x <= air_spacing*0.25){
	    SolnBlk.W[i][j] = Wi;
	  } else {
	    //left at air
	  }
	}
	
	// set conserved state
	SolnBlk.W[i][j].getU( SolnBlk.U[i][j] );

      } 
    } 

    // clean up
    delete[] y;

    //--------------------------------------------------
    // Error.
    //--------------------------------------------------
  } else {
	
    cerr << "\nFlame2DQuadSingleBlock::ICs() - Error, unknown ICs specified.\n";
    exit(-1);

  } // endswitch
  

    /***********************************************************************
     ***********************************************************************
     ***********************************************************************/
 
    /* Set the solution residuals, gradients, limiters, and  other values to zero. */
    
  for (int j  = SolnBlk.JCl-SolnBlk.Nghost ; j <= SolnBlk.JCu+SolnBlk.Nghost ; ++j ) {
    for ( int i = SolnBlk.ICl-SolnBlk.Nghost ; i <= SolnBlk.ICu+SolnBlk.Nghost ; ++i ) {
      for ( int k = 0 ; k <= NUMBER_OF_RESIDUAL_VECTORS_FLAME2D-1 ; ++k ) {
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
    SolnBlk.WoS[i] = SolnBlk.W[i][SolnBlk.JCl];
    SolnBlk.WoN[i] = SolnBlk.W[i][SolnBlk.JCu];
    //   } else if (i < SolnBlk.ICl) {
    //     SolnBlk.WoS[i] = SolnBlk.W[SolnBlk.ICl][SolnBlk.JCl];       //CAUSES FIXED CORNER BC ISSUSES ????? 
    //     SolnBlk.WoN[i] = SolnBlk.W[SolnBlk.ICl][SolnBlk.JCu];        
    //   } else {
    //     SolnBlk.WoS[i] = SolnBlk.W[SolnBlk.ICu][SolnBlk.JCl];
    //      	SolnBlk.WoN[i] = SolnBlk.W[SolnBlk.ICu][SolnBlk.JCu];
    //   } 
  }

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

  int i, j;
  Vector2D dX;

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
	   SolnBlk.Grid.BCtypeS[SolnBlk.ICl-1] == BC_1DFLAME_INFLOW ||
	   SolnBlk.Grid.BCtypeS[SolnBlk.ICl-1] == BC_1DFLAME_OUTFLOW ||
	   SolnBlk.Grid.BCtypeS[SolnBlk.ICl-1] == BC_2DFLAME_INFLOW ||
	   SolnBlk.Grid.BCtypeS[SolnBlk.ICl-1] == BC_2DFLAME_OUTFLOW) ) ||
	 (j > SolnBlk.JCu && 
	  (SolnBlk.Grid.BCtypeN[SolnBlk.ICl-1] == BC_NONE ||
	   SolnBlk.Grid.BCtypeN[SolnBlk.ICl-1] == BC_PERIODIC ||
	   SolnBlk.Grid.BCtypeN[SolnBlk.ICl-1] == BC_CONSTANT_EXTRAPOLATION ||
	   SolnBlk.Grid.BCtypeN[SolnBlk.ICl-1] == BC_LINEAR_EXTRAPOLATION ||
	   SolnBlk.Grid.BCtypeN[SolnBlk.ICl-1] == BC_CHARACTERISTIC ||
	   SolnBlk.Grid.BCtypeN[SolnBlk.ICl-1] == BC_INFLOW_SUBSONIC ||
	   SolnBlk.Grid.BCtypeN[SolnBlk.ICl-1] == BC_OUTFLOW_SUBSONIC || 
	   SolnBlk.Grid.BCtypeN[SolnBlk.ICl-1] == BC_1DFLAME_INFLOW ||
	   SolnBlk.Grid.BCtypeN[SolnBlk.ICl-1] == BC_1DFLAME_OUTFLOW ||
	   SolnBlk.Grid.BCtypeN[SolnBlk.ICl-1] == BC_2DFLAME_INFLOW ||
	   SolnBlk.Grid.BCtypeN[SolnBlk.ICl-1] == BC_2DFLAME_OUTFLOW) ) ) {
      switch(SolnBlk.Grid.BCtypeW[j]) {	 
      case BC_NONE :
	break;                                                                       
      case BC_FIXED :  
	for( int ghost = 1; ghost <= SolnBlk.Nghost; ghost++){
	  SolnBlk.W[SolnBlk.ICl-ghost][j] = SolnBlk.WoW[j];
	  SolnBlk.W[SolnBlk.ICl-ghost][j].getU(SolnBlk.U[SolnBlk.ICl-ghost][j]);
	}
	break;
      case BC_CONSTANT_EXTRAPOLATION :  
	for( int ghost = 1; ghost <= SolnBlk.Nghost; ghost++){	  
	  SolnBlk.W[SolnBlk.ICl-ghost][j] = SolnBlk.W[SolnBlk.ICl][j];
	  SolnBlk.U[SolnBlk.ICl-ghost][j] = SolnBlk.U[SolnBlk.ICl][j];
	}	
	break;
      case BC_LINEAR_EXTRAPOLATION :
	SolnBlk.Linear_Reconstruction_LeastSquares_2(SolnBlk.ICl, j, 
						     LIMITER_BARTH_JESPERSEN);
	for( int ghost = 1; ghost <= SolnBlk.Nghost; ghost++){
	  dX = SolnBlk.Grid.Cell[SolnBlk.ICl-ghost][j].Xc -
	    SolnBlk.Grid.Cell[SolnBlk.ICl][j].Xc;
	  SolnBlk.W[SolnBlk.ICl-ghost][j].Reconstruct( SolnBlk.W[SolnBlk.ICl][j], 
						       SolnBlk.phi[SolnBlk.ICl][j],
						       SolnBlk.dWdx[SolnBlk.ICl][j], 
						       SolnBlk.dWdy[SolnBlk.ICl][j], dX );
	  SolnBlk.W[SolnBlk.ICl-ghost][j].getU(SolnBlk.U[SolnBlk.ICl-ghost][j]);
	}
	break;
      case BC_REFLECTION : 	
	for( int ghost = 1; ghost <= SolnBlk.Nghost; ghost++){
	  SolnBlk.W[SolnBlk.ICl-ghost][j].Reflect(SolnBlk.W[SolnBlk.ICl + ghost-1][j],
						  SolnBlk.Grid.nfaceW(SolnBlk.ICl,j));
	  SolnBlk.W[SolnBlk.ICl-ghost][j].getU(SolnBlk.U[SolnBlk.ICl-ghost][j]);
	}
	break;
      case BC_PERIODIC : 
	for( int ghost = 1; ghost <= SolnBlk.Nghost; ghost++){
	  SolnBlk.W[SolnBlk.ICl-ghost][j] = SolnBlk.W[SolnBlk.ICu-ghost][j];
	  SolnBlk.U[SolnBlk.ICl-ghost][j] = SolnBlk.U[SolnBlk.ICu-ghost][j];
	}
	break;
      case BC_CHARACTERISTIC :   
	SolnBlk.W[SolnBlk.ICl-1][j].BC_Characteristic_Pressure(SolnBlk.W[SolnBlk.ICl][j],
							       SolnBlk.WoW[j],
							       SolnBlk.Grid.nfaceW(SolnBlk.ICl,j));
	SolnBlk.W[SolnBlk.ICl-1][j].getU(SolnBlk.U[SolnBlk.ICl-1][j]);
	for( int ghost = 2; ghost <= SolnBlk.Nghost; ghost++){
	  SolnBlk.W[SolnBlk.ICl-ghost][j] = SolnBlk.W[SolnBlk.ICl-ghost+1][j];
	  SolnBlk.W[SolnBlk.ICl-ghost][j].getU(SolnBlk.U[SolnBlk.ICl-ghost][j]);
	}
	break;
	
      case BC_FREE_SLIP_ISOTHERMAL :
	for( int ghost = 1; ghost <= SolnBlk.Nghost; ghost++){
	  SolnBlk.W[SolnBlk.ICl-ghost][j].Free_Slip(SolnBlk.W[SolnBlk.ICl + ghost-1][j],
						    SolnBlk.WoW[j],
						    SolnBlk.Grid.nfaceW(SolnBlk.ICl,j),
						    FIXED_TEMPERATURE_WALL);
	  SolnBlk.W[SolnBlk.ICl-ghost][j].getU(SolnBlk.U[SolnBlk.ICl-ghost][j]);
	}
	break;
      case BC_WALL_VISCOUS_ISOTHERMAL :
	for( int ghost = 1; ghost <= SolnBlk.Nghost; ghost++){
	  SolnBlk.W[SolnBlk.ICl-ghost][j].No_Slip(SolnBlk.W[SolnBlk.ICl + ghost-1][j],
						  SolnBlk.WoW[j],
						  SolnBlk.Grid.nfaceW(SolnBlk.ICl,j),
						  FIXED_TEMPERATURE_WALL);
	  SolnBlk.W[SolnBlk.ICl-ghost][j].getU(SolnBlk.U[SolnBlk.ICl-ghost][j]);
	}
	break;
      case BC_MOVING_WALL_ISOTHERMAL :	 	
	for( int ghost = 1; ghost <= SolnBlk.Nghost; ghost++){
	  SolnBlk.W[SolnBlk.ICl-ghost][j].Moving_Wall(SolnBlk.W[SolnBlk.ICl + ghost-1][j],
						      SolnBlk.WoW[j],
						      SolnBlk.Grid.nfaceW(SolnBlk.ICl,j),
						      SolnBlk.Moving_wall_velocity,
						      FIXED_TEMPERATURE_WALL);
	  SolnBlk.W[SolnBlk.ICl-ghost][j].getU(SolnBlk.U[SolnBlk.ICl-ghost][j]);
	}	
	break;
      case BC_WALL_VISCOUS_HEATFLUX :  
	for( int ghost = 1; ghost <= SolnBlk.Nghost; ghost++){
	  SolnBlk.W[SolnBlk.ICl-ghost][j].No_Slip(SolnBlk.W[SolnBlk.ICl + ghost-1][j],
						  SolnBlk.WoW[j],
						  SolnBlk.Grid.nfaceW(SolnBlk.ICl,j),
						  ADIABATIC_WALL);
	  SolnBlk.W[SolnBlk.ICl-ghost][j].getU(SolnBlk.U[SolnBlk.ICl-ghost][j]);
	}
	break;
      case BC_MOVING_WALL_HEATFLUX:	
	for( int ghost = 1; ghost <= SolnBlk.Nghost; ghost++){
	  SolnBlk.W[SolnBlk.ICl-ghost][j].Moving_Wall(SolnBlk.W[SolnBlk.ICl + ghost-1][j],
						      SolnBlk.WoW[j],
						      SolnBlk.Grid.nfaceW(SolnBlk.ICl,j),
						      SolnBlk.Moving_wall_velocity,
						      ADIABATIC_WALL);
	  SolnBlk.W[SolnBlk.ICl-ghost][j].getU(SolnBlk.U[SolnBlk.ICl-ghost][j]);
	}	
	break;

      case BC_INFLOW_SUBSONIC :

	// 	// all fixed except u, which is constant extrapolation
	// 	for( int ghost = 1; ghost <= SolnBlk.Nghost; ghost++){
	// 	  SolnBlk.W[SolnBlk.ICl-ghost][j] = SolnBlk.WoW[j];
	// 	  SolnBlk.W[SolnBlk.ICl-ghost][j].v.x = SolnBlk.W[SolnBlk.ICl][j].v.x;
	// 	  SolnBlk.U[SolnBlk.ICl-ghost][j] = U(SolnBlk.W[SolnBlk.ICl-ghost][j]);
	// 	}	 
	// 	break;

	//BC_INFLOW_SUBSONIC (WEST) IS SPECIALLY SET UP FOR THE VISCOUS FLOW   //CURRENTLY NOT FIXED FOR "N" GHOST CELLS
	//ACCURACY ASSESMENT 
	// all fixed, and dpdx is fixed 
	for( int ghost = 1; ghost <= SolnBlk.Nghost; ghost++){
	  dX = SolnBlk.Grid.Cell[SolnBlk.ICl-ghost][j].Xc - SolnBlk.Grid.Cell[SolnBlk.ICl][j].Xc; 
	  SolnBlk.W[SolnBlk.ICl-ghost][j] = SolnBlk.WoW[j];
	  //block pressure using IC values!
	  double pr = SolnBlk.WoW[j].p() +  //SolnBlk.W[SolnBlk.ICl][j].p +
	    ((SolnBlk.WoW[j].p() - SolnBlk.WoE[j].p())/
	     (SolnBlk.Grid.Cell[SolnBlk.ICl][j].Xc.x - 
	      SolnBlk.Grid.Cell[SolnBlk.ICu][j].Xc.x))*dX.x;
	  SolnBlk.W[SolnBlk.ICl-ghost][j].setPressure( pr );
	  SolnBlk.W[SolnBlk.ICl-ghost][j].getU(SolnBlk.U[SolnBlk.ICl-ghost][j]);
	}
	break;
      case BC_OUTFLOW_SUBSONIC :
	// all constant extrapolation except pressure which is fixed.
	for( int ghost = 1; ghost <= SolnBlk.Nghost; ghost++){
	  SolnBlk.W[SolnBlk.ICl-ghost][j] = SolnBlk.W[SolnBlk.ICl][j];
	  SolnBlk.W[SolnBlk.ICl-ghost][j].setPressure( SolnBlk.WoW[j].p() );
	  SolnBlk.W[SolnBlk.ICl-ghost][j].getU(SolnBlk.U[SolnBlk.ICl-ghost][j]);
	}
	break;	
      case BC_1DFLAME_INFLOW :
	SolnBlk.W[SolnBlk.ICl-1][j].BC_1DFlame_Inflow(SolnBlk.W[SolnBlk.ICl][j],
						      SolnBlk.WoW[j], 
						      SolnBlk.W[SolnBlk.ICu][j],
						      SolnBlk.Grid.nfaceW(SolnBlk.ICl,j)); 
	SolnBlk.W[SolnBlk.ICl-1][j].getU(SolnBlk.U[SolnBlk.ICl-1][j]);
	for( int ghost = 2; ghost <= SolnBlk.Nghost; ghost++){
	  SolnBlk.W[SolnBlk.ICl-ghost][j] = SolnBlk.W[SolnBlk.ICl-ghost+1][j];
	  SolnBlk.W[SolnBlk.ICl-ghost][j].getU(SolnBlk.U[SolnBlk.ICl-ghost][j]);
	}
	break;
      case BC_2DFLAME_INFLOW :
	SolnBlk.W[SolnBlk.ICl-1][j].BC_2DFlame_Inflow(SolnBlk.W[SolnBlk.ICl][j],
						      SolnBlk.WoW[j],
						      SolnBlk.Grid.nfaceW(SolnBlk.ICl,j)); 
	SolnBlk.W[SolnBlk.ICl-1][j].getU(SolnBlk.U[SolnBlk.ICl-1][j]);
	for( int ghost = 2; ghost <= SolnBlk.Nghost; ghost++){
	  SolnBlk.W[SolnBlk.ICl-ghost][j] = SolnBlk.W[SolnBlk.ICl-ghost+1][j];
	  SolnBlk.W[SolnBlk.ICl-ghost][j].getU(SolnBlk.U[SolnBlk.ICl-ghost][j]);
	}
	break;
      case BC_1DFLAME_OUTFLOW :
	//isentropic condition for velocity
	SolnBlk.W[SolnBlk.ICl-1][j].BC_1DFlame_Outflow(SolnBlk.W[SolnBlk.ICl][j],
						       SolnBlk.WoW[j],
						       SolnBlk.W[SolnBlk.ICu][j],
						       SolnBlk.Grid.nfaceW(SolnBlk.ICl,j));
	SolnBlk.W[SolnBlk.ICl-1][j].getU(SolnBlk.U[SolnBlk.ICl-1][j]);
	for( int ghost = 2; ghost <= SolnBlk.Nghost; ghost++){
	  SolnBlk.W[SolnBlk.ICl-ghost][j] = SolnBlk.W[SolnBlk.ICl-ghost+1][j];
	  SolnBlk.W[SolnBlk.ICl-ghost][j].getU(SolnBlk.U[SolnBlk.ICl-ghost][j]);
	}
	break;
      case BC_2DFLAME_OUTFLOW :
	SolnBlk.W[SolnBlk.ICl-1][j].BC_2DFlame_Outflow(SolnBlk.W[SolnBlk.ICl][j],
						       SolnBlk.WoW[j],
						       SolnBlk.Grid.nfaceW(SolnBlk.ICl,j));
	SolnBlk.W[SolnBlk.ICl-1][j].getU(SolnBlk.U[SolnBlk.ICl-1][j]);
	for( int ghost = 2; ghost <= SolnBlk.Nghost; ghost++){
	  SolnBlk.W[SolnBlk.ICl-ghost][j] = SolnBlk.W[SolnBlk.ICl-ghost+1][j];
	  SolnBlk.W[SolnBlk.ICl-ghost][j].getU(SolnBlk.U[SolnBlk.ICl-ghost][j]);
	}
	break;
      default: //BC_CONSTANT_EXTRAPOLATION
	for( int ghost = 1; ghost <= SolnBlk.Nghost; ghost++){
	  SolnBlk.W[SolnBlk.ICl-ghost][j] = SolnBlk.W[SolnBlk.ICl][j];
	  SolnBlk.U[SolnBlk.ICl-ghost][j] = SolnBlk.U[SolnBlk.ICl][j];
	}
	break;
      } /* endswitch WEST*/
    } /* endif WEST*/     
  
      //EAST
    if ( (j >= SolnBlk.JCl && j <= SolnBlk.JCu) ||
	 (j < SolnBlk.JCl && 
	  (SolnBlk.Grid.BCtypeS[SolnBlk.ICu+1] == BC_NONE ||
	   SolnBlk.Grid.BCtypeS[SolnBlk.ICu+1] == BC_PERIODIC ||
	   SolnBlk.Grid.BCtypeS[SolnBlk.ICu+1] == BC_CONSTANT_EXTRAPOLATION ||
	   SolnBlk.Grid.BCtypeS[SolnBlk.ICu+1] == BC_LINEAR_EXTRAPOLATION ||
	   SolnBlk.Grid.BCtypeS[SolnBlk.ICu+1] == BC_CHARACTERISTIC ||
	   SolnBlk.Grid.BCtypeS[SolnBlk.ICu+1] == BC_INFLOW_SUBSONIC ||
	   SolnBlk.Grid.BCtypeS[SolnBlk.ICu+1] == BC_OUTFLOW_SUBSONIC || 
	   SolnBlk.Grid.BCtypeS[SolnBlk.ICu+1] == BC_1DFLAME_INFLOW ||
	   SolnBlk.Grid.BCtypeS[SolnBlk.ICu+1] == BC_1DFLAME_OUTFLOW ||
	   SolnBlk.Grid.BCtypeS[SolnBlk.ICu+1] == BC_2DFLAME_INFLOW ||
	   SolnBlk.Grid.BCtypeS[SolnBlk.ICu+1] == BC_2DFLAME_OUTFLOW) ) ||
	 (j > SolnBlk.JCu && 
	  (SolnBlk.Grid.BCtypeN[SolnBlk.ICu+1] == BC_NONE ||
	   SolnBlk.Grid.BCtypeN[SolnBlk.ICu+1] == BC_PERIODIC ||
	   SolnBlk.Grid.BCtypeN[SolnBlk.ICu+1] == BC_CONSTANT_EXTRAPOLATION ||
	   SolnBlk.Grid.BCtypeN[SolnBlk.ICu+1] == BC_LINEAR_EXTRAPOLATION ||
	   SolnBlk.Grid.BCtypeN[SolnBlk.ICu+1] == BC_CHARACTERISTIC ||
	   SolnBlk.Grid.BCtypeN[SolnBlk.ICu+1] == BC_INFLOW_SUBSONIC ||
	   SolnBlk.Grid.BCtypeN[SolnBlk.ICu+1] == BC_OUTFLOW_SUBSONIC || 
	   SolnBlk.Grid.BCtypeN[SolnBlk.ICu+1] == BC_1DFLAME_INFLOW ||
	   SolnBlk.Grid.BCtypeN[SolnBlk.ICu+1] == BC_1DFLAME_OUTFLOW ||
	   SolnBlk.Grid.BCtypeN[SolnBlk.ICu+1] == BC_2DFLAME_INFLOW ||
	   SolnBlk.Grid.BCtypeN[SolnBlk.ICu+1] == BC_2DFLAME_OUTFLOW) ) ) { 
    
      switch(SolnBlk.Grid.BCtypeE[j]) {
      case BC_NONE :
	break;
      case BC_FIXED : 
	for( int ghost = 1; ghost <= SolnBlk.Nghost; ghost++){
	  SolnBlk.W[SolnBlk.ICu+ghost][j] = SolnBlk.WoE[j];
	  SolnBlk.W[SolnBlk.ICu+ghost][j].getU(SolnBlk.U[SolnBlk.ICu+ghost][j]);
	}
	break;
      case BC_CONSTANT_EXTRAPOLATION :
	for( int ghost = 1; ghost <= SolnBlk.Nghost; ghost++){ 	  
	  SolnBlk.W[SolnBlk.ICu+ghost][j] = SolnBlk.W[SolnBlk.ICu][j];
	  SolnBlk.U[SolnBlk.ICu+ghost][j] = SolnBlk.U[SolnBlk.ICu][j]; 
	}
	break;
      case BC_LINEAR_EXTRAPOLATION :
	SolnBlk.Linear_Reconstruction_LeastSquares_2(SolnBlk.ICu, j, 
						     LIMITER_BARTH_JESPERSEN); 
	for( int ghost = 1; ghost <= SolnBlk.Nghost; ghost++){
	  dX = SolnBlk.Grid.Cell[SolnBlk.ICu+ghost][j].Xc -
	    SolnBlk.Grid.Cell[SolnBlk.ICu][j].Xc;
	  SolnBlk.W[SolnBlk.ICu+ghost][j].Reconstruct( SolnBlk.W[SolnBlk.ICu][j],
						       SolnBlk.phi[SolnBlk.ICu][j],
						       SolnBlk.dWdx[SolnBlk.ICu][j],
						       SolnBlk.dWdy[SolnBlk.ICu][j], dX );
	  SolnBlk.W[SolnBlk.ICu+ghost][j].getU(SolnBlk.U[SolnBlk.ICu+ghost][j]);
	}
	break;
      case BC_REFLECTION :
	for( int ghost = 1; ghost <= SolnBlk.Nghost; ghost++){
	  SolnBlk.W[SolnBlk.ICu+ghost][j].Reflect(SolnBlk.W[SolnBlk.ICu-ghost+1][j],
						  SolnBlk.Grid.nfaceE(SolnBlk.ICu,j));
	  SolnBlk.W[SolnBlk.ICu+ghost][j].getU(SolnBlk.U[SolnBlk.ICu+ghost][j]);
	}
	break;
      case BC_PERIODIC :  
	for( int ghost = 1; ghost <= SolnBlk.Nghost; ghost++){
	  SolnBlk.W[SolnBlk.ICu+ghost][j] = SolnBlk.W[SolnBlk.ICl+ghost][j];
	  SolnBlk.U[SolnBlk.ICu+ghost][j] = SolnBlk.U[SolnBlk.ICl+ghost][j];
	}
	break;
      case BC_CHARACTERISTIC :
	SolnBlk.W[SolnBlk.ICu+1][j].BC_Characteristic_Pressure(SolnBlk.W[SolnBlk.ICu][j],
							       SolnBlk.WoE[j],
							       SolnBlk.Grid.nfaceE(SolnBlk.ICu,j));
	SolnBlk.W[SolnBlk.ICu+1][j].getU(SolnBlk.U[SolnBlk.ICu+1][j]); 
	for( int ghost = 2; ghost <= SolnBlk.Nghost; ghost++){
	  SolnBlk.W[SolnBlk.ICu+ghost][j] = SolnBlk.W[SolnBlk.ICu+ghost-1][j];
	  SolnBlk.W[SolnBlk.ICu+ghost][j].getU(SolnBlk.U[SolnBlk.ICu+ghost][j]);
	}
	break;
      case BC_FREE_SLIP_ISOTHERMAL :
	for( int ghost = 1; ghost <= SolnBlk.Nghost; ghost++){
	  SolnBlk.W[SolnBlk.ICu+ghost][j].Free_Slip(SolnBlk.W[SolnBlk.ICu-ghost+1][j],
						    SolnBlk.WoE[j],
						    SolnBlk.Grid.nfaceE(SolnBlk.ICu,j),
						    FIXED_TEMPERATURE_WALL);
	  SolnBlk.W[SolnBlk.ICu+ghost][j].getU(SolnBlk.U[SolnBlk.ICu+ghost][j]);
	}
	break;
      case BC_WALL_VISCOUS_ISOTHERMAL :
	for( int ghost = 1; ghost <= SolnBlk.Nghost; ghost++){
	  SolnBlk.W[SolnBlk.ICu+ghost][j].No_Slip(SolnBlk.W[SolnBlk.ICu-ghost+1][j],
						  SolnBlk.WoE[j],
						  SolnBlk.Grid.nfaceE(SolnBlk.ICu,j),
						  FIXED_TEMPERATURE_WALL);
	  SolnBlk.W[SolnBlk.ICu+ghost][j].getU(SolnBlk.U[SolnBlk.ICu+ghost][j]);
	}       
	break;
      case BC_MOVING_WALL_ISOTHERMAL :
	for( int ghost = 1; ghost <= SolnBlk.Nghost; ghost++){
	  SolnBlk.W[SolnBlk.ICu+ghost][j].Moving_Wall(SolnBlk.W[SolnBlk.ICu-ghost+1][j],
						      SolnBlk.WoE[j],
						      SolnBlk.Grid.nfaceE(SolnBlk.ICu,j),
						      SolnBlk.Moving_wall_velocity,
						      FIXED_TEMPERATURE_WALL);
	  SolnBlk.W[SolnBlk.ICu+ghost][j].getU(SolnBlk.U[SolnBlk.ICu+ghost][j]);
	}      
	break;	
      case BC_WALL_VISCOUS_HEATFLUX :
	for( int ghost = 1; ghost <= SolnBlk.Nghost; ghost++){
	  SolnBlk.W[SolnBlk.ICu+ghost][j].No_Slip(SolnBlk.W[SolnBlk.ICu-ghost+1][j],
						  SolnBlk.WoE[j],
						  SolnBlk.Grid.nfaceE(SolnBlk.ICu,j),
						  ADIABATIC_WALL);   
	  SolnBlk.W[SolnBlk.ICu+ghost][j].getU(SolnBlk.U[SolnBlk.ICu+ghost][j]);
	}      
	break;
      case BC_MOVING_WALL_HEATFLUX :
	for( int ghost = 1; ghost <= SolnBlk.Nghost; ghost++){
	  SolnBlk.W[SolnBlk.ICu+ghost][j].Moving_Wall(SolnBlk.W[SolnBlk.ICu-ghost+1][j],
						      SolnBlk.WoE[j],
						      SolnBlk.Grid.nfaceE(SolnBlk.ICu,j),
						      SolnBlk.Moving_wall_velocity,
						      ADIABATIC_WALL);  
	  SolnBlk.W[SolnBlk.ICu+ghost][j].getU(SolnBlk.U[SolnBlk.ICu+ghost][j] );
	}     
	break;	
      case BC_INFLOW_SUBSONIC :
	// all fixed except v.x (u) which is constant extrapolation
	for( int ghost = 1; ghost <= SolnBlk.Nghost; ghost++){
	  SolnBlk.W[SolnBlk.ICu+ghost][j] = SolnBlk.WoE[j];
	  
	  SolnBlk.W[SolnBlk.ICu+ghost][j].setVelocityX( ((const Flame2D_pState&)SolnBlk.W[SolnBlk.ICu][j]).vx() );
	  SolnBlk.W[SolnBlk.ICu+ghost][j].getU(SolnBlk.U[SolnBlk.ICu+ghost][j]);
	}	 
	break;
	
      case BC_OUTFLOW_SUBSONIC :
	// 	//all constant except pressure which is FIXED
	// 	for( int ghost = 1; ghost <= SolnBlk.Nghost; ghost++){
	// 	  SolnBlk.W[SolnBlk.ICu+ghost][j] = SolnBlk.W[SolnBlk.ICu][j];
	// 	  SolnBlk.W[SolnBlk.ICu+ghost][j].p = SolnBlk.WoE[j].p;
	// 	  SolnBlk.U[SolnBlk.ICu+ghost][j] =  U(SolnBlk.W[SolnBlk.ICu+ghost][j]);
	// 	}
	// 	break;
  
	//SPECIAL FOR VISOCOUS ACCURACY PIPE/CHANNEL FLOWS
	// all constant extrapolation except pressure which linearly extrapolated
	for( int ghost = 1; ghost <= SolnBlk.Nghost; ghost++){	
	  dX = SolnBlk.Grid.Cell[SolnBlk.ICu+ghost][j].Xc - SolnBlk.Grid.Cell[SolnBlk.ICu][j].Xc; 
	  SolnBlk.W[SolnBlk.ICu+ghost][j] = SolnBlk.W[SolnBlk.ICu][j]; 	
	  //block pressure using IC values!
	  double pr = SolnBlk.WoE[j].p() + //SolnBlk.W[SolnBlk.ICu][j].p +
	    ((SolnBlk.WoW[j].p() - SolnBlk.WoE[j].p())/
	     (SolnBlk.Grid.Cell[SolnBlk.ICl][j].Xc.x - 
	      SolnBlk.Grid.Cell[SolnBlk.ICu][j].Xc.x))*dX.x;
	  SolnBlk.W[SolnBlk.ICu+ghost][j].setPressure( pr );
	  SolnBlk.W[SolnBlk.ICu+ghost][j].getU(SolnBlk.U[SolnBlk.ICu+1][j]);
	}
	break;
	
      case BC_1DFLAME_INFLOW :	
	SolnBlk.W[SolnBlk.ICu+1][j].BC_1DFlame_Inflow(SolnBlk.W[SolnBlk.ICu][j],
						      SolnBlk.WoE[j], 
						      SolnBlk.W[SolnBlk.ICl][j],
						      SolnBlk.Grid.nfaceE(SolnBlk.ICu,j)); 	
	SolnBlk.W[SolnBlk.ICu+1][j].getU(SolnBlk.U[SolnBlk.ICu+1][j]);
	for( int ghost = 2; ghost <= SolnBlk.Nghost; ghost++){
	  SolnBlk.W[SolnBlk.ICu+ghost][j] = SolnBlk.W[SolnBlk.ICu+ghost-1][j];
	  SolnBlk.W[SolnBlk.ICu+ghost][j].getU(SolnBlk.U[SolnBlk.ICu+ghost][j]);
	}
	break;
      case BC_1DFLAME_OUTFLOW : 
	//isentropic condition for velocity ??
	SolnBlk.W[SolnBlk.ICu+1][j].BC_1DFlame_Outflow(SolnBlk.W[SolnBlk.ICu][j],
						       SolnBlk.WoE[j],
						       SolnBlk.W[SolnBlk.ICl][j],
						       SolnBlk.Grid.nfaceE(SolnBlk.ICu,j));
	SolnBlk.W[SolnBlk.ICu+1][j].getU(SolnBlk.U[SolnBlk.ICu+1][j]);
	for( int ghost = 2; ghost <= SolnBlk.Nghost; ghost++){
	  SolnBlk.W[SolnBlk.ICu+ghost][j] = SolnBlk.W[SolnBlk.ICu+ghost-1][j];
	  SolnBlk.W[SolnBlk.ICu+ghost][j].getU(SolnBlk.U[SolnBlk.ICu+ghost][j]);
	}
	break; 
      case BC_2DFLAME_INFLOW :	
	SolnBlk.W[SolnBlk.ICu+1][j].BC_2DFlame_Inflow(SolnBlk.W[SolnBlk.ICu][j],
						      SolnBlk.WoE[j], 		
						      SolnBlk.Grid.nfaceE(SolnBlk.ICu,j)); 	
	SolnBlk.W[SolnBlk.ICu+1][j].getU(SolnBlk.U[SolnBlk.ICu+1][j]);
	for( int ghost = 2; ghost <= SolnBlk.Nghost; ghost++){
	  SolnBlk.W[SolnBlk.ICu+ghost][j] = SolnBlk.W[SolnBlk.ICu+ghost-1][j];
	  SolnBlk.W[SolnBlk.ICu+ghost][j].getU(SolnBlk.U[SolnBlk.ICu+ghost][j] );
	}
	break;	 
      case BC_2DFLAME_OUTFLOW :
	SolnBlk.W[SolnBlk.ICu+1][j].BC_2DFlame_Outflow(SolnBlk.W[SolnBlk.ICu][j],
						       SolnBlk.WoE[j],			   
						       SolnBlk.Grid.nfaceE(SolnBlk.ICu,j));
	SolnBlk.W[SolnBlk.ICu+1][j].getU(SolnBlk.U[SolnBlk.ICu+1][j]);
	for( int ghost = 2; ghost <= SolnBlk.Nghost; ghost++){
	  SolnBlk.W[SolnBlk.ICu+ghost][j] = SolnBlk.W[SolnBlk.ICu+ghost-1][j];
	  SolnBlk.W[SolnBlk.ICu+ghost][j].getU(SolnBlk.U[SolnBlk.ICu+ghost][j]);
	}
	break;
      default:
	for( int ghost = 1; ghost <= SolnBlk.Nghost; ghost++){
	  SolnBlk.W[SolnBlk.ICu+ghost][j] = SolnBlk.W[SolnBlk.ICu][j];
	  SolnBlk.U[SolnBlk.ICu+ghost][j] = SolnBlk.U[SolnBlk.ICu][j];
	}
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
      for( int ghost = 1; ghost <= SolnBlk.Nghost; ghost++){
	SolnBlk.W[i][SolnBlk.JCl-ghost] = SolnBlk.WoS[i];
	SolnBlk.W[i][SolnBlk.JCl-ghost].getU(SolnBlk.U[i][SolnBlk.JCl-ghost]);
      }
      break;
    case BC_CONSTANT_EXTRAPOLATION :
      for( int ghost = 1; ghost <= SolnBlk.Nghost; ghost++){
	SolnBlk.W[i][SolnBlk.JCl-ghost] = SolnBlk.W[i][SolnBlk.JCl];
	SolnBlk.U[i][SolnBlk.JCl-ghost] = SolnBlk.U[i][SolnBlk.JCl]; 
      }
      break;
    case BC_LINEAR_EXTRAPOLATION :
      SolnBlk.Linear_Reconstruction_LeastSquares_2(i, SolnBlk.JCl, 
						   LIMITER_BARTH_JESPERSEN); 
      for( int ghost = 1; ghost <= SolnBlk.Nghost; ghost++){
	dX = SolnBlk.Grid.Cell[i][SolnBlk.JCl-ghost].Xc -
	  SolnBlk.Grid.Cell[i][SolnBlk.JCl].Xc;
	SolnBlk.W[i][SolnBlk.JCl-ghost].Reconstruct( SolnBlk.W[i][SolnBlk.JCl],
						     SolnBlk.phi[i][SolnBlk.JCl],
						     SolnBlk.dWdx[i][SolnBlk.JCl],
						     SolnBlk.dWdy[i][SolnBlk.JCl], dX);
	SolnBlk.W[i][SolnBlk.JCl-ghost].getU(SolnBlk.U[i][SolnBlk.JCl-ghost]);
      }
      break;
    case BC_REFLECTION :  
      for( int ghost = 1; ghost <= SolnBlk.Nghost; ghost++){
	SolnBlk.W[i][SolnBlk.JCl-ghost].Reflect(SolnBlk.W[i][SolnBlk.JCl + ghost-1],      
						SolnBlk.Grid.nfaceS(i, SolnBlk.JCl));
	SolnBlk.W[i][SolnBlk.JCl-ghost].getU(SolnBlk.U[i][SolnBlk.JCl-ghost]);
      }
      break;
    case BC_PERIODIC : 
      for( int ghost = 1; ghost <= SolnBlk.Nghost; ghost++){
	SolnBlk.W[i][SolnBlk.JCl-ghost] = SolnBlk.W[i][SolnBlk.JCu-ghost];
	SolnBlk.U[i][SolnBlk.JCl-ghost] = SolnBlk.U[i][SolnBlk.JCu-ghost];
      }
      break;
    case BC_CHARACTERISTIC :
      SolnBlk.W[i][SolnBlk.JCl-1].BC_Characteristic_Pressure(SolnBlk.W[i][SolnBlk.JCl],
							     SolnBlk.WoS[i],
							     SolnBlk.Grid.nfaceS(i, SolnBlk.JCl));
      SolnBlk.W[i][SolnBlk.JCl-1].getU(SolnBlk.U[i][SolnBlk.JCl-1]);
      for( int ghost = 2; ghost <= SolnBlk.Nghost; ghost++){
	SolnBlk.W[i][SolnBlk.JCl-ghost] = SolnBlk.W[i][SolnBlk.JCl-ghost+1];
	SolnBlk.W[i][SolnBlk.JCl-ghost].getU(SolnBlk.U[i][SolnBlk.JCl-ghost]);
      }
      break;
    case BC_FREE_SLIP_ISOTHERMAL :
      for( int ghost = 1; ghost <= SolnBlk.Nghost; ghost++){
	SolnBlk.W[i][SolnBlk.JCl-ghost].Free_Slip(SolnBlk.W[i][SolnBlk.JCl + ghost-1],     
						  SolnBlk.WoS[i],
						  SolnBlk.Grid.nfaceS(i, SolnBlk.JCl),
						  FIXED_TEMPERATURE_WALL);
	SolnBlk.W[i][SolnBlk.JCl-ghost].getU(SolnBlk.U[i][SolnBlk.JCl-ghost]);
      }
      break;
    case BC_WALL_VISCOUS_ISOTHERMAL :  
      for( int ghost = 1; ghost <= SolnBlk.Nghost; ghost++){
	SolnBlk.W[i][SolnBlk.JCl-ghost].No_Slip(SolnBlk.W[i][SolnBlk.JCl + ghost-1],
						SolnBlk.WoS[i],
						SolnBlk.Grid.nfaceS(i, SolnBlk.JCl),
						FIXED_TEMPERATURE_WALL);
	SolnBlk.W[i][SolnBlk.JCl-ghost].getU(SolnBlk.U[i][SolnBlk.JCl-ghost]);
      }
      break;
    case BC_MOVING_WALL_ISOTHERMAL :
      for( int ghost = 1; ghost <= SolnBlk.Nghost; ghost++){
	SolnBlk.W[i][SolnBlk.JCl-ghost].Moving_Wall(SolnBlk.W[i][SolnBlk.JCl + ghost-1],
						    SolnBlk.WoS[i],
						    SolnBlk.Grid.nfaceS(i, SolnBlk.JCl),
						    SolnBlk.Moving_wall_velocity,
						    FIXED_TEMPERATURE_WALL);
	SolnBlk.W[i][SolnBlk.JCl-ghost].getU(SolnBlk.U[i][SolnBlk.JCl-ghost]);
      }
      break; 
    case BC_WALL_VISCOUS_HEATFLUX : 
      for( int ghost = 1; ghost <= SolnBlk.Nghost; ghost++){
	SolnBlk.W[i][SolnBlk.JCl-ghost].No_Slip(SolnBlk.W[i][SolnBlk.JCl + ghost-1],
						SolnBlk.WoS[i],
						SolnBlk.Grid.nfaceS(i, SolnBlk.JCl),
						ADIABATIC_WALL);  
	SolnBlk.W[i][SolnBlk.JCl-ghost].getU(SolnBlk.U[i][SolnBlk.JCl-ghost]);
      }
      break;
    case BC_MOVING_WALL_HEATFLUX :
      for( int ghost = 1; ghost <= SolnBlk.Nghost; ghost++){
	SolnBlk.W[i][SolnBlk.JCl-ghost].Moving_Wall(SolnBlk.W[i][SolnBlk.JCl + ghost-1],
						    SolnBlk.WoS[i],
						    SolnBlk.Grid.nfaceS(i, SolnBlk.JCl),
						    SolnBlk.Moving_wall_velocity,
						    ADIABATIC_WALL);  
	SolnBlk.W[i][SolnBlk.JCl-ghost].getU(SolnBlk.U[i][SolnBlk.JCl-ghost]);
      }
      break; 
    case BC_INFLOW_SUBSONIC :
      // all fixed except v.x (u) which is constant extrapolation
      for( int ghost = 1; ghost <= SolnBlk.Nghost; ghost++){
	SolnBlk.W[i][SolnBlk.JCl-ghost] = SolnBlk.WoS[i];
	SolnBlk.W[i][SolnBlk.JCl-ghost].setVelocityY( ((const Flame2D_pState&)SolnBlk.W[i][SolnBlk.JCl]).vy() );
	SolnBlk.W[i][SolnBlk.JCl-ghost].getU(SolnBlk.U[i][SolnBlk.JCl-ghost]);
      }
      break;
    case BC_OUTFLOW_SUBSONIC :
      // all constant extrapolation except pressure which is fixed.
      for( int ghost = 1; ghost <= SolnBlk.Nghost; ghost++){
	SolnBlk.W[i][SolnBlk.JCl-ghost] = SolnBlk.W[i][SolnBlk.JCl];
	SolnBlk.W[i][SolnBlk.JCl-ghost].setPressure( SolnBlk.WoS[i].p() );
	SolnBlk.W[i][SolnBlk.JCl-ghost].getU(SolnBlk.U[i][SolnBlk.JCl-ghost]);
      }      
      break;
    
    case BC_1DFLAME_INFLOW :       
      SolnBlk.W[i][SolnBlk.JCl-1].BC_1DFlame_Inflow(SolnBlk.W[i][SolnBlk.JCl],
						    SolnBlk.WoS[i], 
						    SolnBlk.W[i][SolnBlk.JCu],
						    SolnBlk.Grid.nfaceS(i,SolnBlk.JCl)); 
      SolnBlk.W[i][SolnBlk.JCl-1].getU(SolnBlk.U[i][SolnBlk.JCl-1]);
      for( int ghost = 2; ghost <= SolnBlk.Nghost; ghost++){
	SolnBlk.W[i][SolnBlk.JCl-ghost] = SolnBlk.W[i][SolnBlk.JCl-ghost+1];
	SolnBlk.W[i][SolnBlk.JCl-ghost].getU(SolnBlk.U[i][SolnBlk.JCl-ghost]);
      }
      break;  
    case BC_2DFLAME_INFLOW :       
      SolnBlk.W[i][SolnBlk.JCl-1].BC_2DFlame_Inflow(SolnBlk.W[i][SolnBlk.JCl],
						    SolnBlk.WoS[i], 
						    SolnBlk.Grid.nfaceS(i,SolnBlk.JCl)); 
      SolnBlk.W[i][SolnBlk.JCl-1].getU(SolnBlk.U[i][SolnBlk.JCl-1]);
      for( int ghost = 2; ghost <= SolnBlk.Nghost; ghost++){
	SolnBlk.W[i][SolnBlk.JCl-ghost] = SolnBlk.W[i][SolnBlk.JCl-ghost+1];
	SolnBlk.W[i][SolnBlk.JCl-ghost].getU(SolnBlk.U[i][SolnBlk.JCl-ghost]);
      }       
      break;      
    case BC_1DFLAME_OUTFLOW :
      cerr<<"\n BC_FLAME_OUTFLOW South doesn't exist ";
      break;
    case BC_2DFLAME_OUTFLOW :
      cerr<<"\n BC_2DFLAME_OUTFLOW not setup for SOUTH ";
      break;
    default:
      for( int ghost = 1; ghost <= SolnBlk.Nghost; ghost++){
	SolnBlk.W[i][SolnBlk.JCl-ghost] = SolnBlk.W[i][SolnBlk.JCl];
	SolnBlk.U[i][SolnBlk.JCl-ghost] = SolnBlk.U[i][SolnBlk.JCl];
      }
      break;
    } /* endswitch South */

      //NORTH
    switch(SolnBlk.Grid.BCtypeN[i]) {
    case BC_NONE :
      break;
    case BC_FIXED : 
      for( int ghost = 1; ghost <= SolnBlk.Nghost; ghost++){
	SolnBlk.W[i][SolnBlk.JCu+ghost] = SolnBlk.WoN[i];
	SolnBlk.W[i][SolnBlk.JCu+ghost].getU(SolnBlk.U[i][SolnBlk.JCu+ghost]);
      }
      break;
    case BC_CONSTANT_EXTRAPOLATION :      
      for( int ghost = 1; ghost <= SolnBlk.Nghost; ghost++){	
	SolnBlk.W[i][SolnBlk.JCu+ghost] = SolnBlk.W[i][SolnBlk.JCu];
	SolnBlk.U[i][SolnBlk.JCu+ghost] = SolnBlk.U[i][SolnBlk.JCu];
      }
      break;
    case BC_LINEAR_EXTRAPOLATION :
      SolnBlk.Linear_Reconstruction_LeastSquares_2(i, SolnBlk.JCu, 
						   LIMITER_BARTH_JESPERSEN);
      for( int ghost = 1; ghost <= SolnBlk.Nghost; ghost++){
	dX = SolnBlk.Grid.Cell[i][SolnBlk.JCu+ghost].Xc -
	  SolnBlk.Grid.Cell[i][SolnBlk.JCu].Xc;
	SolnBlk.W[i][SolnBlk.JCu+ghost].Reconstruct( SolnBlk.W[i][SolnBlk.JCu],
						     SolnBlk.phi[i][SolnBlk.JCu],
						     SolnBlk.dWdx[i][SolnBlk.JCu],
						     SolnBlk.dWdy[i][SolnBlk.JCu], dX);
	SolnBlk.W[i][SolnBlk.JCu+ghost].getU(SolnBlk.U[i][SolnBlk.JCu+ghost]);
      }
      break;
    case BC_REFLECTION :  
      for( int ghost = 1; ghost <= SolnBlk.Nghost; ghost++){
	SolnBlk.W[i][SolnBlk.JCu+ghost].Reflect(SolnBlk.W[i][SolnBlk.JCu-ghost+1],
						SolnBlk.Grid.nfaceN(i, SolnBlk.JCu));
	SolnBlk.W[i][SolnBlk.JCu+ghost].getU(SolnBlk.U[i][SolnBlk.JCu+ghost]);
      }      
      break;      
    case BC_PERIODIC :
      for( int ghost = 1; ghost <= SolnBlk.Nghost; ghost++){
	SolnBlk.W[i][SolnBlk.JCu+ghost] = SolnBlk.W[i][SolnBlk.JCl+ghost];
	SolnBlk.U[i][SolnBlk.JCu+ghost] = SolnBlk.U[i][SolnBlk.JCl+ghost];
      }
      break;
    case BC_CHARACTERISTIC :
      SolnBlk.W[i][SolnBlk.JCu+1].BC_Characteristic_Pressure(SolnBlk.W[i][SolnBlk.JCu],
							     SolnBlk.WoN[i],
							     SolnBlk.Grid.nfaceN(i, SolnBlk.JCu));
      SolnBlk.W[i][SolnBlk.JCu+1].getU(SolnBlk.U[i][SolnBlk.JCu+1]);
      for( int ghost = 2; ghost <= SolnBlk.Nghost; ghost++){
	SolnBlk.W[i][SolnBlk.JCu+ghost] = SolnBlk.W[i][SolnBlk.JCu+ghost-1];
	SolnBlk.W[i][SolnBlk.JCu+ghost].getU(SolnBlk.U[i][SolnBlk.JCu+ghost]);
      }
      break;
    case BC_FREE_SLIP_ISOTHERMAL :
      for( int ghost = 1; ghost <= SolnBlk.Nghost; ghost++){
	SolnBlk.W[i][SolnBlk.JCu+ghost].Free_Slip(SolnBlk.W[i][SolnBlk.JCu-ghost+1],
						  SolnBlk.WoN[i],
						  SolnBlk.Grid.nfaceN(i, SolnBlk.JCu),
						  FIXED_TEMPERATURE_WALL);
	SolnBlk.W[i][SolnBlk.JCu+ghost].getU(SolnBlk.U[i][SolnBlk.JCu+ghost]);
      }   
    case BC_WALL_VISCOUS_ISOTHERMAL :
      for( int ghost = 1; ghost <= SolnBlk.Nghost; ghost++){
	SolnBlk.W[i][SolnBlk.JCu+ghost].No_Slip(SolnBlk.W[i][SolnBlk.JCu-ghost+1],
						SolnBlk.WoN[i],
						SolnBlk.Grid.nfaceN(i, SolnBlk.JCu),
						FIXED_TEMPERATURE_WALL);
	SolnBlk.W[i][SolnBlk.JCu+ghost].getU(SolnBlk.U[i][SolnBlk.JCu+ghost]);
      }   
      break;
    case BC_MOVING_WALL_ISOTHERMAL :
      for( int ghost = 1; ghost <= SolnBlk.Nghost; ghost++){
	SolnBlk.W[i][SolnBlk.JCu+ghost].Moving_Wall(SolnBlk.W[i][SolnBlk.JCu-ghost+1],
						    SolnBlk.WoN[i],
						    SolnBlk.Grid.nfaceN(i, SolnBlk.JCu),
						    SolnBlk.Moving_wall_velocity,
						    FIXED_TEMPERATURE_WALL);
	SolnBlk.W[i][SolnBlk.JCu+ghost].getU(SolnBlk.U[i][SolnBlk.JCu+ghost]);
      }  
      break; 
    case BC_WALL_VISCOUS_HEATFLUX :
      for( int ghost = 1; ghost <= SolnBlk.Nghost; ghost++){
	SolnBlk.W[i][SolnBlk.JCu+ghost].No_Slip(SolnBlk.W[i][SolnBlk.JCu-ghost+1],
						SolnBlk.WoN[i],
						SolnBlk.Grid.nfaceN(i, SolnBlk.JCu),
						ADIABATIC_WALL);
	SolnBlk.W[i][SolnBlk.JCu+ghost].getU(SolnBlk.U[i][SolnBlk.JCu+ghost]);
      }   
      break;
    case BC_MOVING_WALL_HEATFLUX :
      for( int ghost = 1; ghost <= SolnBlk.Nghost; ghost++){
	SolnBlk.W[i][SolnBlk.JCu+ghost].Moving_Wall(SolnBlk.W[i][SolnBlk.JCu-ghost+1],
						    SolnBlk.WoN[i],
						    SolnBlk.Grid.nfaceN(i, SolnBlk.JCu),
						    SolnBlk.Moving_wall_velocity,
						    ADIABATIC_WALL);
	SolnBlk.W[i][SolnBlk.JCu+ghost].getU(SolnBlk.U[i][SolnBlk.JCu+ghost]);
      }       
      break; 
    case BC_INFLOW_SUBSONIC :
      // all fixed except v.x (u) which is constant extrapolation
      for( int ghost = 1; ghost <= SolnBlk.Nghost; ghost++){
	SolnBlk.W[i][SolnBlk.JCu+ghost] = SolnBlk.WoN[i];
	SolnBlk.W[i][SolnBlk.JCu+ghost].setVelocityY( ((const Flame2D_pState&)SolnBlk.W[i][SolnBlk.JCu]).vy() );
	SolnBlk.W[i][SolnBlk.JCu+ghost].getU(SolnBlk.U[i][SolnBlk.JCu+ghost]);
      }
      break;
    case BC_OUTFLOW_SUBSONIC :
      // all constant extrapolation except pressure which is fixed.
      for( int ghost = 1; ghost <= SolnBlk.Nghost; ghost++){
	SolnBlk.W[i][SolnBlk.JCu+ghost] = SolnBlk.W[i][SolnBlk.JCu];
	SolnBlk.W[i][SolnBlk.JCu+ghost].setPressure( SolnBlk.WoN[i].p() );
	SolnBlk.W[i][SolnBlk.JCu+ghost].getU(SolnBlk.U[i][SolnBlk.JCu+ghost]);
      }
      break;      
    case BC_1DFLAME_INFLOW :
      cerr<<"\n BC_1DFLAME_INFLOW doesn't exist for North BC ";
      break;	
    case BC_2DFLAME_INFLOW :
      cerr<<"\n BC_2DFLAME_INFLOW doesn't exist for North BC ";
      break;	
    case BC_1DFLAME_OUTFLOW :
      SolnBlk.W[i][SolnBlk.JCu+1].BC_1DFlame_Outflow(SolnBlk.W[i][SolnBlk.JCu],
						     SolnBlk.WoN[i],
						     SolnBlk.W[i][SolnBlk.JCl],
						     SolnBlk.Grid.nfaceN(i,SolnBlk.JCu));
      SolnBlk.W[i][SolnBlk.JCu+1].getU(SolnBlk.U[i][SolnBlk.JCu+1]);
      for( int ghost = 2; ghost <= SolnBlk.Nghost; ghost++){
	SolnBlk.W[i][SolnBlk.JCu+ghost] = SolnBlk.W[i][SolnBlk.JCu+ghost-1];
	SolnBlk.W[i][SolnBlk.JCu+ghost].getU(SolnBlk.U[i][SolnBlk.JCu+ghost]);
      }	
      break; 
    case BC_2DFLAME_OUTFLOW :
      SolnBlk.W[i][SolnBlk.JCu+1].BC_2DFlame_Outflow(SolnBlk.W[i][SolnBlk.JCu],
						     SolnBlk.WoN[i],
						     SolnBlk.Grid.nfaceN(i, SolnBlk.JCu));
      SolnBlk.W[i][SolnBlk.JCu+1].getU(SolnBlk.U[i][SolnBlk.JCu+1]);
      for( int ghost = 2; ghost <= SolnBlk.Nghost; ghost++){
	SolnBlk.W[i][SolnBlk.JCu+ghost] = SolnBlk.W[i][SolnBlk.JCu+ghost-1];
	SolnBlk.W[i][SolnBlk.JCu+ghost].getU(SolnBlk.U[i][SolnBlk.JCu+ghost]);
      }
      break;	
    default:
      for( int ghost = 1; ghost <= SolnBlk.Nghost; ghost++){
	SolnBlk.W[i][SolnBlk.JCu+ghost] = SolnBlk.W[i][SolnBlk.JCu];
	SolnBlk.U[i][SolnBlk.JCu+ghost] = SolnBlk.U[i][SolnBlk.JCu];
      }
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
	

	/******** Inviscid deltat calculation **********************/   
	//no preconditioning
	
	if(Input_Parameters.Preconditioning == 0){
	  
	  a = SolnBlk.W[i][j].a();	  
	  SolnBlk.dt[i][j] = min(d_i/(a+fabs(v_i)), d_j/(a+fabs(v_j)));
	  
	  //Low Mach Number Preconditioning
	} else if(Input_Parameters.Preconditioning == 1) { 
	  
	  delta_n = min(fabs(d_i),fabs(d_j));
	  SolnBlk.dt[i][j] = min(d_i/SolnBlk.W[i][j].u_plus_aprecon(fabs(v_i),
								    SolnBlk.Flow_Type,
								    delta_n),
				 d_j/SolnBlk.W[i][j].u_plus_aprecon(fabs(v_j),
								    SolnBlk.Flow_Type,
								    delta_n));	    
	}
	
	/******** Viscous deltat calculation ************************/   
	if (SolnBlk.Flow_Type != FLOWTYPE_INVISCID && 
	    Input_Parameters.Preconditioning == 0) {  
	  
	  SolnBlk.W[i][j].updateViscosity();
	  rhomu = SolnBlk.W[i][j].mu()/((const Flame2D_pState &)SolnBlk.W[i][j]).rho();
	  dt_vis = min((d_i*d_i)/(TWO*rhomu), (d_j*d_j)/(TWO*rhomu)); 
	  SolnBlk.dt[i][j]  = min(dt_vis, SolnBlk.dt[i][j]);
	}	  
	
	/******** Chemical Source Term deltat calculation ************/   
	dt_chem = HALF/SolnBlk.W[i][j].dSwdU_max_diagonal();
	dt_chem *= Input_Parameters.Source_Term_Multiplyer; // scale source term
	SolnBlk.dt[i][j] = min(dt_chem, SolnBlk.dt[i][j]);	  

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

  for (int j = SolnBlk.JCl-SolnBlk.Nghost ; j <= SolnBlk.JCu+SolnBlk.Nghost ; ++j ) {
    for (int i = SolnBlk.ICl-SolnBlk.Nghost ; i <= SolnBlk.ICu+SolnBlk.Nghost ; ++i ) {
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
  double rho, p;
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
      SolnBlk.Linear_Reconstruction_LeastSquares_2(i, j, LIMITER_UNLIMITED);

      if (SolnBlk.Grid.Cell[i][j].A > ZERO) {

	// store rho and p
	rho = ((const Flame2D_pState &)SolnBlk.W[i][j]).rho();
	p = ((const Flame2D_pState &)SolnBlk.W[i][j]).p();

	// Evaluate refinement criteria #1 based on the gradient
	// of the density field.
	grad_rho_x = SolnBlk.dWdx[i][j].rho();
	grad_rho_y = SolnBlk.dWdy[i][j].rho();
	grad_rho_abs = sqrt(sqr(grad_rho_x) + sqr(grad_rho_y));
	grad_rho_criteria = sqrt(SolnBlk.Grid.Cell[i][j].A)*grad_rho_abs/rho;
	grad_rho_criteria_max = max(grad_rho_criteria_max, grad_rho_criteria);

	// Evaluate refinement criteria #2 based on the divergence
	// of the velocity vector.
	div_V = SolnBlk.dWdx[i][j].vx() + SolnBlk.dWdy[i][j].vy();
	div_V_criteria = sqrt(SolnBlk.Grid.Cell[i][j].A)*fabs(div_V)/SolnBlk.W[i][j].a();
	div_V_criteria_max = max(div_V_criteria_max, div_V_criteria);

	// Evaluate refinement criteria #3 based on the curl
	// of the velocity vector.
	curl_V_z = SolnBlk.dWdx[i][j].vy() - SolnBlk.dWdy[i][j].vx(); 
	curl_V_abs = sqrt(sqr(curl_V_z)); 
	curl_V_criteria = sqrt(SolnBlk.Grid.Cell[i][j].A)*curl_V_abs/SolnBlk.W[i][j].a(); 
	curl_V_criteria_max = max(curl_V_criteria_max, curl_V_criteria);

	// 	     Evaluate refinement criteria #4 based on the gradient
	// 	     of the Temperature
	grad_Temp_x = (ONE/( rho* SolnBlk.W[i][j].Rtot())) * 
	  (SolnBlk.dWdx[i][j].p() - (p/rho) * SolnBlk.dWdx[i][j].rho());
	grad_Temp_y = (ONE/( rho* SolnBlk.W[i][j].Rtot())) *
	  (SolnBlk.dWdy[i][j].p() - (p/rho) * SolnBlk.dWdy[i][j].rho());
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
      SolnBlk.U[i][SolnBlk.JCu].set( SolnBlk.U[i][SolnBlk.JCu], 
				     (SolnBlk.Grid.Cell[i][SolnBlk.JCu].A/
				      SolnBlk.Grid.area(i, SolnBlk.JCu)) );
      SolnBlk.W[i][SolnBlk.JCu].setU(SolnBlk.U[i][SolnBlk.JCu]);
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
      SolnBlk.U[i][SolnBlk.JCl].set( SolnBlk.U[i][SolnBlk.JCl],
				     (SolnBlk.Grid.Cell[i][SolnBlk.JCl].A/
				      SolnBlk.Grid.area(i, SolnBlk.JCl)) );
      SolnBlk.W[i][SolnBlk.JCl].setU(SolnBlk.U[i][SolnBlk.JCl]);
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
      SolnBlk.U[SolnBlk.ICu][j].set( SolnBlk.U[SolnBlk.ICu][j],
				     (SolnBlk.Grid.Cell[SolnBlk.ICu][j].A/
				      SolnBlk.Grid.area(SolnBlk.ICu, j)) );
      SolnBlk.W[SolnBlk.ICu][j].setU(SolnBlk.U[SolnBlk.ICu][j]);
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
      SolnBlk.U[SolnBlk.ICl][j].set( SolnBlk.U[SolnBlk.ICl][j],
				     (SolnBlk.Grid.Cell[SolnBlk.ICl][j].A/
				      SolnBlk.Grid.area(SolnBlk.ICl, j)) );
      SolnBlk.W[SolnBlk.ICl][j].setU(SolnBlk.U[SolnBlk.ICl][j]);
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
      SolnBlk.U[i][SolnBlk.JCu].set( SolnBlk.U[i][SolnBlk.JCu],
				     (SolnBlk.Grid.Cell[i][SolnBlk.JCu].A/
				      SolnBlk.Grid.area(i, SolnBlk.JCu)) );
      SolnBlk.W[i][SolnBlk.JCu].setU(SolnBlk.U[i][SolnBlk.JCu]);
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
      SolnBlk.U[i][SolnBlk.JCl].set( SolnBlk.U[i][SolnBlk.JCl],
				     (SolnBlk.Grid.Cell[i][SolnBlk.JCl].A/
				      SolnBlk.Grid.area(i, SolnBlk.JCl)) );
      SolnBlk.W[i][SolnBlk.JCl].setU(SolnBlk.U[i][SolnBlk.JCl]);
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
      SolnBlk.U[SolnBlk.ICu][j].set( SolnBlk.U[SolnBlk.ICu][j],
				     (SolnBlk.Grid.Cell[SolnBlk.ICu][j].A/
				      SolnBlk.Grid.area(SolnBlk.ICu, j)) );
      SolnBlk.W[SolnBlk.ICu][j].setU(SolnBlk.U[SolnBlk.ICu][j]);
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
      SolnBlk.U[SolnBlk.ICl][j].set( SolnBlk.U[SolnBlk.ICl][j],
				     (SolnBlk.Grid.Cell[SolnBlk.ICl][j].A/
				      SolnBlk.Grid.area(SolnBlk.ICl, j)) );
      SolnBlk.W[SolnBlk.ICl][j].setU(SolnBlk.U[SolnBlk.ICl][j]);
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
      SolnBlk.dUdt[i][SolnBlk.JCu][0].add( SolnBlk.FluxN[i], 
					   -1.0/SolnBlk.Grid.Cell[i][SolnBlk.JCu].A );
    } /* endfor */
  } /* endif */

    /* Correct the fluxes at the south boundary as required. */

  if (Number_Neighbours_South_Boundary == 2) {
    for ( i = SolnBlk.ICl ; i <= SolnBlk.ICu ; ++i ) {
      SolnBlk.dUdt[i][SolnBlk.JCl][0].add( SolnBlk.FluxS[i],
					   -1.0/SolnBlk.Grid.Cell[i][SolnBlk.JCl].A );
    } /* endfor */
  } /* endif */

    /* Correct the fluxes at the east boundary as required. */

  if (Number_Neighbours_East_Boundary == 2) {
    for ( j= SolnBlk.JCl ; j <= SolnBlk.JCu ; ++j ) {
      SolnBlk.dUdt[SolnBlk.ICu][j][0].add( SolnBlk.FluxE[j],
					   -1.0/SolnBlk.Grid.Cell[SolnBlk.ICu][j].A );
    } /* endfor */
  } /* endif */

    /* Correct the fluxes at the west boundary as required. */

  if (Number_Neighbours_West_Boundary == 2) {
    for ( j= SolnBlk.JCl ; j <= SolnBlk.JCu ; ++j ) {
      SolnBlk.dUdt[SolnBlk.ICl][j][0].add( SolnBlk.FluxW[j],
					   -1.0/SolnBlk.Grid.Cell[SolnBlk.ICl][j].A );
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
      SolnBlk.dUdt[i][SolnBlk.JCu][k_residual].add( SolnBlk.FluxN[i],
						    -(CFL_Number*SolnBlk.dt[i][SolnBlk.JCu]/
						      SolnBlk.Grid.Cell[i][SolnBlk.JCu].A) );
    } /* endfor */
  } /* endif */

    /* Correct the fluxes at the south boundary as required. */

  if (Number_Neighbours_South_Boundary == 2) {
    for ( i = SolnBlk.ICl ; i <= SolnBlk.ICu ; ++i ) {
      SolnBlk.dUdt[i][SolnBlk.JCl][k_residual].add( SolnBlk.FluxS[i],
						    -(CFL_Number*SolnBlk.dt[i][SolnBlk.JCl]/
						      SolnBlk.Grid.Cell[i][SolnBlk.JCl].A) );
    } /* endfor */
  } /* endif */

    /* Correct the fluxes at the east boundary as required. */

  if (Number_Neighbours_East_Boundary == 2) {
    for ( j= SolnBlk.JCl ; j <= SolnBlk.JCu ; ++j ) {
      SolnBlk.dUdt[SolnBlk.ICu][j][k_residual].add( SolnBlk.FluxE[j],
						    -(CFL_Number*SolnBlk.dt[SolnBlk.ICu][j]/
						      SolnBlk.Grid.Cell[SolnBlk.ICu][j].A) );
    } /* endfor */
  } /* endif */

    /* Correct the fluxes at the west boundary as required. */

  if (Number_Neighbours_West_Boundary == 2) {
    for ( j= SolnBlk.JCl ; j <= SolnBlk.JCu ; ++j ) {
      SolnBlk.dUdt[SolnBlk.ICl][j][k_residual].add( SolnBlk.FluxW[j],
						    -(CFL_Number*SolnBlk.dt[SolnBlk.ICl][j]/
						      SolnBlk.Grid.Cell[SolnBlk.ICl][j].A) );
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


  // declares
  static Vector2D dX;
  static Flame2D_pState Wl, Wr, W_face;
  static Flame2D_State Flux;
  const int NUM_VAR_FLAME2D = SolnBlk.NumVar();
  double delta_n; 

  //-----------------------------------------------------------

  // Modifications for NKS overlap 
  int JCl_overlap = 0, JCu_overlap = 0, ICu_overlap = 0, ICl_overlap = 0;
  if(Input_Parameters.NKS_IP.GMRES_Overlap > 0 ){	
    if (SolnBlk.Grid.BCtypeS[SolnBlk.ICl] == BC_NONE)  JCl_overlap = Input_Parameters.NKS_IP.GMRES_Overlap;
    if (SolnBlk.Grid.BCtypeN[SolnBlk.ICu] == BC_NONE)  JCu_overlap = Input_Parameters.NKS_IP.GMRES_Overlap;
    if (SolnBlk.Grid.BCtypeE[SolnBlk.JCu] == BC_NONE)  ICu_overlap = Input_Parameters.NKS_IP.GMRES_Overlap;
    if (SolnBlk.Grid.BCtypeW[SolnBlk.JCl] == BC_NONE)  ICl_overlap = Input_Parameters.NKS_IP.GMRES_Overlap;    
  }

  //-----------------------------------------------------------

  // Perform the linear reconstruction within each cell of the
  // computational grid for this stage.
  switch(Input_Parameters.i_Reconstruction) {
  case RECONSTRUCTION_GREEN_GAUSS :            
    if(Input_Parameters.i_Viscous_Flux_Evaluation == VISCOUS_RECONSTRUCTION_DIAMONDPATH_GREEN_GAUSS)
      SolnBlk.Linear_Reconstruction_GreenGauss_Diamond(Input_Parameters.i_Limiter);
    else
      SolnBlk.Linear_Reconstruction_GreenGauss(Input_Parameters.i_Limiter);    
    break;
  case RECONSTRUCTION_LEAST_SQUARES :
    if(Input_Parameters.i_Viscous_Flux_Evaluation ==VISCOUS_RECONSTRUCTION_DIAMONDPATH_LEAST_SQUARES)      
      SolnBlk.Linear_Reconstruction_LeastSquares_Diamond(Input_Parameters.i_Limiter);       
    else
      SolnBlk.Linear_Reconstruction_LeastSquares(Input_Parameters.i_Limiter);
    break;
  default:
    SolnBlk.Linear_Reconstruction_LeastSquares(Input_Parameters.i_Limiter);
    break;
  }
  
  //-----------------------------------------------------------

  // radiation evaluation if using optically thin approx
  if ( Input_Parameters.Radiation == RADIATION_OPTICALLY_THIN )
    Radiation_Source_Eval( SolnBlk, Input_Parameters );
  

  // Evaluate the time rate of change of the solution
  // (i.e., the solution residuals) using a second-order
  // limited upwind scheme with a variety of flux functions.
  
  //*****************************************************************
  // Add i-direction (zeta-direction) fluxes.
  //*****************************************************************
  for (int j  = SolnBlk.JCl-1-JCl_overlap ; j <= SolnBlk.JCu+1+JCu_overlap ; j++ ) {
    SolnBlk.dUdt[SolnBlk.ICl-1][j][0].Vacuum();
      
    for (int i = SolnBlk.ICl-1-ICl_overlap; i <= SolnBlk.ICu+ICu_overlap ; i++ ) {	 
      SolnBlk.dUdt[i+1][j][0].Vacuum();
      
      if ( j >= SolnBlk.JCl-JCl_overlap && j <= SolnBlk.JCu+JCu_overlap ) { 


	// Determine left and right states
	SolnBlk.Reconstructed_LeftandRight_States( Wl, Wr, i, j, X_DIRECTION );
    
	// Spacing for Preconditioner 
	if(Input_Parameters.Preconditioning){
	  delta_n = SolnBlk.delta_n(i,j);
	} 

	//-----------------------------------------------------------

	// Determine EAST face INVISCID flux.
	switch(Input_Parameters.i_Flux_Function) {
	case FLUX_FUNCTION_ROE :
	  Flux.FluxRoe_n(Wl, Wr, SolnBlk.Grid.nfaceE(i, j),
			 Input_Parameters.Preconditioning,SolnBlk.Flow_Type,delta_n);
	  break;
	case FLUX_FUNCTION_HLLE :  
	  Flux.FluxHLLE_n(Wl, Wr, SolnBlk.Grid.nfaceE(i, j));
	  break;
	case FLUX_FUNCTION_LINDE :
	  Flux.FluxLinde_n(Wl, Wr, SolnBlk.Grid.nfaceE(i, j));
	  break; 
	case FLUX_FUNCTION_AUSM_PLUS_UP :
	  Flux.FluxAUSMplus_up_n(Wl, Wr, SolnBlk.Grid.nfaceE(i, j));
	  break;
	} // endswitch
	
	  // Determine EAST face VISCOUS flux.
	if(SolnBlk.Flow_Type != FLOWTYPE_INVISCID){ // -ve as on RHS 	  
	  switch(Input_Parameters.i_Viscous_Flux_Evaluation){
	  case VISCOUS_RECONSTRUCTION_DIAMONDPATH_LEAST_SQUARES :	      
	  case VISCOUS_RECONSTRUCTION_DIAMONDPATH_GREEN_GAUSS :
	    // Wn temporarily stored in Wnd -> calculated in reconstruction
	    W_face.Average(SolnBlk.Wnd[i+1][j+1], SolnBlk.Wnd[i+1][j]); // IS this right in general ??( Wn_NE + Wn_SE ) 
	    Flux.Viscous_Flux_n(W_face,
				SolnBlk.dWdx_faceE[i][j],
				SolnBlk.dWdy_faceE[i][j],
				SolnBlk.Axisymmetric,
				SolnBlk.Grid.xfaceE(i,j),
				SolnBlk.Grid.nfaceE(i,j),
				-1.0); // subtract from inviscid flux
	    break;
	  case VISCOUS_RECONSTRUCTION_HYBRID :
	    W_face.Average(SolnBlk.W[i][j], SolnBlk.W[i+1][j]);
	    Flux.Viscous_FluxHybrid_n(W_face,
				      SolnBlk.dWdx_faceE[i][j],
				      SolnBlk.dWdy_faceE[i][j],
				      SolnBlk.Grid.xfaceE(i,j),
				      SolnBlk.W[i][j],
				      SolnBlk.dWdx[i][j],
				      SolnBlk.dWdy[i][j],
				      SolnBlk.Grid.Cell[i][j].Xc,
				      SolnBlk.W[i+1][j],
				      SolnBlk.dWdy[i+1][j],
				      SolnBlk.dWdy[i+1][j],
				      SolnBlk.Grid.Cell[i+1][j].Xc,
				      SolnBlk.Axisymmetric,
				      SolnBlk.Grid.nfaceE(i,j),
				      -1.0); // subtract from inviscid flux
	    SolnBlk.dWdx_faceW[i+1][j] = SolnBlk.dWdx_faceE[i][j];
	    SolnBlk.dWdy_faceW[i+1][j] = SolnBlk.dWdy_faceE[i][j];
	    break;
	  } // endswitch
	} // endif - viscous
	
	  //-----------------------------------------------------------

	  // Evaluate cell-averaged solution changes.
	SolnBlk.dUdt[i][j][0].add(Flux, -SolnBlk.Grid.lfaceE(i, j)/SolnBlk.Grid.Cell[i][j].A);
	SolnBlk.dUdt[i+1][j][0].add(Flux, SolnBlk.Grid.lfaceW(i+1, j)/SolnBlk.Grid.Cell[i+1][j].A);
	
	// Include axisymmetric source terms as required.
	if (SolnBlk.Axisymmetric) {
	  // inviscid
	  SolnBlk.W[i][j].Sa_inviscid(SolnBlk.dUdt[i][j][0], 
				      SolnBlk.Grid.Cell[i][j].Xc,
				      SolnBlk.Axisymmetric);

	  // viscous
	  if (SolnBlk.Flow_Type != FLOWTYPE_INVISCID) { //-ve as on RHS 
	    SolnBlk.W[i][j].Sa_viscous(SolnBlk.dUdt[i][j][0],
				       SolnBlk.dWdx[i][j],
				       SolnBlk.dWdy[i][j],
				       SolnBlk.Grid.Cell[i][j].Xc,
				       SolnBlk.Axisymmetric);
	  }
	}
		
	// Include source terms associated with the finite-rate chemistry
	SolnBlk.W[i][j].Sw(SolnBlk.dUdt[i][j][0]);
	
		
	// Include source terms associated with gravity
	if (SolnBlk.Gravity) {	 
	  SolnBlk.W[i][j].Sg(SolnBlk.dUdt[i][j][0]);
	}   
	
	// Include source terms associated with radiation
	SolnBlk.dUdt[i][j][0].E() += SolnBlk.Srad[i][j];
	
	//-----------------------------------------------------------
	
	// Save west and east face boundary flux.
	if (i == SolnBlk.ICl-1) {
	  SolnBlk.FluxW[j].set(Flux, -SolnBlk.Grid.lfaceW(i+1, j));
	} else if (i == SolnBlk.ICu) {
	  SolnBlk.FluxE[j].set(Flux,  SolnBlk.Grid.lfaceE(i, j));
	} 
	
      } // end if 
    } // end for i
    
    if ( j > SolnBlk.JCl-1-JCl_overlap && j < SolnBlk.JCu+1+JCu_overlap ) {  
      SolnBlk.dUdt[SolnBlk.ICu+1][j][0].Vacuum();
      SolnBlk.dUdt[SolnBlk.ICl-1][j][0].Vacuum();
    }
  } // end for j
      
  
    //*****************************************************************
    // Add j-direction (eta-direction) fluxes.
    //*****************************************************************
  for ( int i = SolnBlk.ICl-ICl_overlap ; i <= SolnBlk.ICu+ICu_overlap ; ++i ) {
    for ( int j = SolnBlk.JCl-1-JCl_overlap ; j <= SolnBlk.JCu+JCu_overlap ; ++j ) {
      
      // Determine left and right states
      SolnBlk.Reconstructed_LeftandRight_States( Wl, Wr, i, j, Y_DIRECTION );

      // Spacing for Preconditioner 
      if(Input_Parameters.Preconditioning){
	delta_n = SolnBlk.delta_n(i,j);
      } 

      //-----------------------------------------------------------

      // Determine NORTH face inviscid flux.
      switch(Input_Parameters.i_Flux_Function) {
      case FLUX_FUNCTION_ROE :
	Flux.FluxRoe_n(Wl, Wr, SolnBlk.Grid.nfaceN(i, j),
		       Input_Parameters.Preconditioning,SolnBlk.Flow_Type,delta_n);
	break;
      case FLUX_FUNCTION_HLLE :
	Flux.FluxHLLE_n(Wl, Wr, SolnBlk.Grid.nfaceN(i, j));
	break;
      case FLUX_FUNCTION_LINDE :
	Flux.FluxLinde_n(Wl, Wr, SolnBlk.Grid.nfaceN(i, j));
	break;
      case FLUX_FUNCTION_AUSM_PLUS_UP :
	Flux.FluxAUSMplus_up_n(Wl, Wr, SolnBlk.Grid.nfaceN(i, j));
	break;
      } // endswitch
	
	// Determine NORTH face viscous flux.
      if(SolnBlk.Flow_Type != FLOWTYPE_INVISCID){ // -ve as on RHS    
	switch(Input_Parameters.i_Viscous_Flux_Evaluation){	    
	case VISCOUS_RECONSTRUCTION_DIAMONDPATH_LEAST_SQUARES :	      
	case VISCOUS_RECONSTRUCTION_DIAMONDPATH_GREEN_GAUSS :	  	  
	  // Wn temporarily stored in Wnd -> calculated in reconstruction
	  W_face.Average(SolnBlk.Wnd[i][j+1], SolnBlk.Wnd[i+1][j+1]);
	  Flux.Viscous_Flux_n(W_face,
			      SolnBlk.dWdx_faceN[i][j],
			      SolnBlk.dWdy_faceN[i][j],
			      SolnBlk.Axisymmetric,
			      SolnBlk.Grid.xfaceN(i,j),
			      SolnBlk.Grid.nfaceN(i,j),
			      -1.0); // subtract from inviscid flux
	  break;
	case VISCOUS_RECONSTRUCTION_HYBRID :
	  W_face.Average(SolnBlk.W[i][j], SolnBlk.W[i][j+1]);
	  Flux.Viscous_FluxHybrid_n(W_face,
				    SolnBlk.dWdx_faceN[i][j],
				    SolnBlk.dWdy_faceN[i][j],
				    SolnBlk.Grid.xfaceN(i,j),
				    SolnBlk.W[i][j],
				    SolnBlk.dWdx[i][j],
				    SolnBlk.dWdy[i][j],
				    SolnBlk.Grid.Cell[i][j].Xc,
				    SolnBlk.W[i][j+1],
				    SolnBlk.dWdy[i][j+1],
				    SolnBlk.dWdy[i][j+1],
				    SolnBlk.Grid.Cell[i][j+1].Xc,
				    SolnBlk.Axisymmetric,
				    SolnBlk.Grid.nfaceN(i,j),
				    -1.0); // subtract from inviscid flux
	  SolnBlk.dWdx_faceS[i][j+1] = SolnBlk.dWdx_faceN[i][j];
	  SolnBlk.dWdy_faceS[i][j+1] = SolnBlk.dWdy_faceN[i][j];
	  break;
	} // endswitch
      }// endif - viscous
	
      //-----------------------------------------------------------
	
      // Evaluate cell-averaged solution changes.
      SolnBlk.dUdt[i][j][0].add( Flux, -SolnBlk.Grid.lfaceN(i, j)/SolnBlk.Grid.Cell[i][j].A );
      SolnBlk.dUdt[i][j+1][0].add( Flux, SolnBlk.Grid.lfaceS(i, j+1)/SolnBlk.Grid.Cell[i][j+1].A );
	
      // Save south and north face boundary flux.
      // USED for AMR 
      if (j == SolnBlk.JCl-1) {
	SolnBlk.FluxS[i].set( Flux, -SolnBlk.Grid.lfaceS(i, j+1) );
      } else if (j == SolnBlk.JCu) {
	SolnBlk.FluxN[i].set( Flux,  SolnBlk.Grid.lfaceN(i, j) );
      } /* endif */
	

    } // endfor - j

      // reset ghost cell solution changes
    SolnBlk.dUdt[i][SolnBlk.JCu+1][0].Vacuum();
    SolnBlk.dUdt[i][SolnBlk.JCl-1][0].Vacuum();

  } //  endfor - i


    /* Residual successfully calculated. */

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
  
  // declares
  int i, j, k_residual;
  double omega;
  static Vector2D dX;
  static Flame2D_pState Wl, Wr, W_face;
  static Flame2D_State Flux;
  const int NUM_VAR_FLAME2D = SolnBlk.NumVar();
  double delta_n;
  
  // Evaluate the solution residual for stage i_stage of n_stage scheme.

  //-----------------------------------------------------------

  // Evaluate the time step fraction and residual 
  // storage location for the stage.
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
      if (i_stage == 4)	k_residual = 0;
      else k_residual = i_stage - 1;
    }
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
  }
    
  //-----------------------------------------------------------

  // Perform the linear reconstruction within each cell of the
  // computational grid for this stage.
  switch(Input_Parameters.i_Reconstruction) {
  case RECONSTRUCTION_GREEN_GAUSS :
    if(Input_Parameters.i_Viscous_Flux_Evaluation ==VISCOUS_RECONSTRUCTION_DIAMONDPATH_GREEN_GAUSS)
      SolnBlk.Linear_Reconstruction_GreenGauss_Diamond(Input_Parameters.i_Limiter);
    else
      SolnBlk.Linear_Reconstruction_GreenGauss(Input_Parameters.i_Limiter);    
    break;
  case RECONSTRUCTION_LEAST_SQUARES :
    if(Input_Parameters.i_Viscous_Flux_Evaluation ==VISCOUS_RECONSTRUCTION_DIAMONDPATH_LEAST_SQUARES)
      SolnBlk.Linear_Reconstruction_LeastSquares_Diamond(Input_Parameters.i_Limiter);          
    else
      SolnBlk.Linear_Reconstruction_LeastSquares(Input_Parameters.i_Limiter);
    break;
  default:
    SolnBlk.Linear_Reconstruction_LeastSquares(Input_Parameters.i_Limiter);
    break;
  }

  //-----------------------------------------------------------

  // radiation evaluation if using opticallyt htin approx
  if ( Input_Parameters.Radiation == RADIATION_OPTICALLY_THIN )
    Radiation_Source_Eval( SolnBlk, Input_Parameters );

  //-----------------------------------------------------------


  // Evaluate the time rate of change of the solution
  // (i.e., the solution residuals) using a second-order
  // limited upwind scheme with a variety of flux functions.
  
  //*****************************************************************
  // Add i-direction (zeta-direction) fluxes.
  //*****************************************************************
  for ( j  = SolnBlk.JCl-1 ; j <= SolnBlk.JCu+1 ; ++j ) {

    // Initialize states depending upon the stage
    if ( i_stage == 1 ) {
      SolnBlk.Uo[SolnBlk.ICl-1][j] = SolnBlk.U[SolnBlk.ICl-1][j];
      SolnBlk.dUdt[SolnBlk.ICl-1][j][k_residual].Vacuum();
    } else {
      SolnBlk.dUdt[SolnBlk.ICl-1][j][k_residual].Vacuum();
    }
      
    for ( i = SolnBlk.ICl-1 ; i <= SolnBlk.ICu ; ++i ) {

      //-----------------------------------------------------------

      // First stage
      if ( i_stage == 1 ) {
	SolnBlk.Uo[i+1][j] = SolnBlk.U[i+1][j];
	SolnBlk.dUdt[i+1][j][k_residual].Vacuum();

	// update dUdt
      } else if ( j > SolnBlk.JCl-1 && j < SolnBlk.JCu+1 ) {
	switch(Input_Parameters.i_Time_Integration) {
	case TIME_STEPPING_EXPLICIT_PREDICTOR_CORRECTOR :
	  //SolnBlk.dUdt[i+1][j][k_residual] = SolnBlk.dUdt[i+1][j][k_residual];
	  break;
	case TIME_STEPPING_EXPLICIT_RUNGE_KUTTA :
	  if (Input_Parameters.N_Stage == 2) {
	    //SolnBlk.dUdt[i+1][j][k_residual] = SolnBlk.dUdt[i+1][j][k_residual];
	  } else if (Input_Parameters.N_Stage == 4 && i_stage == 4) {
	    SolnBlk.dUdt[i+1][j][k_residual] = 
	      SolnBlk.dUdt[i+1][j][0] + 
	      TWO*SolnBlk.dUdt[i+1][j][1] +
	      TWO*SolnBlk.dUdt[i+1][j][2];
	  } else {
	    SolnBlk.dUdt[i+1][j][k_residual].Vacuum();
	  }
	  break;
	case TIME_STEPPING_MULTISTAGE_OPTIMAL_SMOOTHING :
	  SolnBlk.dUdt[i+1][j][k_residual].Vacuum(); 
	  break;
	default:
	  SolnBlk.dUdt[i+1][j][k_residual].Vacuum(); 
	  break;
	} /* endswitch */

      } // endif - stage

	//-----------------------------------------------------------

      if ( j >= SolnBlk.JCl && j <= SolnBlk.JCu ) {
	 
	// Determine left and right states
	SolnBlk.Reconstructed_LeftandRight_States( Wl, Wr, i, j, X_DIRECTION );
	
	// Spacing for Preconditioner 
	if(Input_Parameters.Preconditioning){
	  delta_n = SolnBlk.delta_n(i,j);
	}
	  
	//-----------------------------------------------------------

	// Determine EAST face INVISCID flux.
	switch(Input_Parameters.i_Flux_Function) {
	case FLUX_FUNCTION_ROE :
	  Flux.FluxRoe_n(Wl, Wr, SolnBlk.Grid.nfaceE(i, j),
			 Input_Parameters.Preconditioning,SolnBlk.Flow_Type,delta_n);
	  break;
	case FLUX_FUNCTION_HLLE :
	  Flux.FluxHLLE_n(Wl, Wr, SolnBlk.Grid.nfaceE(i, j));
	  break;
	case FLUX_FUNCTION_LINDE :
	  Flux.FluxLinde_n(Wl, Wr, SolnBlk.Grid.nfaceE(i, j));
	  break;
	case FLUX_FUNCTION_AUSM_PLUS_UP :
	  Flux.FluxAUSMplus_up_n(Wl, Wr, SolnBlk.Grid.nfaceE(i, j));
	  break;
	} // endswitch
	  	
	  // Determine EAST face VISCOUS flux.
	if(SolnBlk.Flow_Type != FLOWTYPE_INVISCID){ // -ve as on RHS 	  
	  switch(Input_Parameters.i_Viscous_Flux_Evaluation){
	  case VISCOUS_RECONSTRUCTION_DIAMONDPATH_LEAST_SQUARES :	      
	  case VISCOUS_RECONSTRUCTION_DIAMONDPATH_GREEN_GAUSS :
	    // Wn temporarily stored in Wnd -> calculated in reconstruction
	    W_face.Average(SolnBlk.Wnd[i+1][j+1], SolnBlk.Wnd[i+1][j]); // IS this right in general ??( Wn_NE + Wn_SE ) 
	    Flux.Viscous_Flux_n(W_face,
				SolnBlk.dWdx_faceE[i][j],
				SolnBlk.dWdy_faceE[i][j],
				SolnBlk.Axisymmetric,
				SolnBlk.Grid.xfaceE(i,j),
				SolnBlk.Grid.nfaceE(i,j),
				-1.0); // subtract from inviscid flux
	    break;
	  case VISCOUS_RECONSTRUCTION_HYBRID :
	    W_face.Average(SolnBlk.W[i][j], SolnBlk.W[i+1][j]);
	    Flux.Viscous_FluxHybrid_n(W_face,
				      SolnBlk.dWdx_faceE[i][j],
				      SolnBlk.dWdy_faceE[i][j],
				      SolnBlk.Grid.xfaceE(i,j),
				      SolnBlk.W[i][j],
				      SolnBlk.dWdx[i][j],
				      SolnBlk.dWdy[i][j],
				      SolnBlk.Grid.Cell[i][j].Xc,
				      SolnBlk.W[i+1][j],
				      SolnBlk.dWdy[i+1][j],
				      SolnBlk.dWdy[i+1][j],
				      SolnBlk.Grid.Cell[i+1][j].Xc,
				      SolnBlk.Axisymmetric,
				      SolnBlk.Grid.nfaceE(i,j),
				      -1.0); // subtract from inviscid flux
	    SolnBlk.dWdx_faceW[i+1][j] = SolnBlk.dWdx_faceE[i][j];
	    SolnBlk.dWdy_faceW[i+1][j] = SolnBlk.dWdy_faceE[i][j];
	    break;
	  } // endswitch
	} // endif - viscous

	  //-----------------------------------------------------------

	  // Evaluate cell-averaged solution changes
	SolnBlk.dUdt[i][j][k_residual].add( Flux, -(Input_Parameters.CFL_Number*
						    SolnBlk.dt[i][j]*SolnBlk.Grid.lfaceE(i, j)/
						    SolnBlk.Grid.Cell[i][j].A) );

	SolnBlk.dUdt[i+1][j][k_residual].add( Flux, (Input_Parameters.CFL_Number*
						     SolnBlk.dt[i+1][j]*SolnBlk.Grid.lfaceW(i+1, j)/
						     SolnBlk.Grid.Cell[i+1][j].A) );
	  
	// Include axisymmetric source terms as required.
	if (SolnBlk.Axisymmetric) {
	  SolnBlk.W[i][j].Sa_inviscid(SolnBlk.dUdt[i][j][k_residual],
				      SolnBlk.Grid.Cell[i][j].Xc,
				      SolnBlk.Axisymmetric,
				      Input_Parameters.CFL_Number*SolnBlk.dt[i][j]);
	  // Include Viscous if specified	 
	  if (SolnBlk.Flow_Type != FLOWTYPE_INVISCID) { 
	    SolnBlk.W[i][j].Sa_viscous(SolnBlk.dUdt[i][j][k_residual],
				       SolnBlk.dWdx[i][j],
				       SolnBlk.dWdy[i][j],
				       SolnBlk.Grid.Cell[i][j].Xc,
				       SolnBlk.Axisymmetric,
				       Input_Parameters.CFL_Number*SolnBlk.dt[i][j]);
	  }
	} // endif

	  // Include source terms associated with the finite-rate chemistry
	  //rho*omega_dot
	SolnBlk.W[i][j].Sw(SolnBlk.dUdt[i][j][k_residual], 
			   Input_Parameters.CFL_Number*SolnBlk.dt[i][j]);
 
	// Include source terms associated with gravity
	if (SolnBlk.Gravity) {	 
	  SolnBlk.W[i][j].Sg(SolnBlk.dUdt[i][j][k_residual],
			     Input_Parameters.CFL_Number*SolnBlk.dt[i][j]);
	}

	// Include source terms associated with radiation
	SolnBlk.dUdt[i][j][k_residual].E() += 
	  (Input_Parameters.CFL_Number*SolnBlk.dt[i][j])*SolnBlk.Srad[i][j];

	//-----------------------------------------------------------


	// Save west and east face boundary flux.
	// USED for AMR 	
	if (i == SolnBlk.ICl-1) {
	  SolnBlk.FluxW[j].set( Flux, -SolnBlk.Grid.lfaceW(i+1, j) );
	} else if (i == SolnBlk.ICu) {
	  SolnBlk.FluxE[j].set( Flux,  SolnBlk.Grid.lfaceE(i, j) );
	} 
	
      } // endif
    } // endfor i
                
    if ( j > SolnBlk.JCl-1 && j < SolnBlk.JCu+1 ){
      SolnBlk.dUdt[SolnBlk.ICu+1][j][k_residual].Vacuum();
      SolnBlk.dUdt[SolnBlk.ICl-1][j][k_residual].Vacuum();
    }
  } // endfor j

    //*****************************************************************
    // Add j-direction (eta-direction) fluxes.
    //*****************************************************************
  for ( i = SolnBlk.ICl ; i <= SolnBlk.ICu ; ++i ) {
    for ( j  = SolnBlk.JCl-1 ; j <= SolnBlk.JCu ; ++j ) {

      // Determine left and right states
      SolnBlk.Reconstructed_LeftandRight_States( Wl, Wr, i, j, Y_DIRECTION );
		
      // Spacing for Preconditioner 
      if(Input_Parameters.Preconditioning){
	delta_n = SolnBlk.delta_n(i,j);
      }
      
      //-----------------------------------------------------------

      // Determine NORTH face inviscid flux.
      switch(Input_Parameters.i_Flux_Function) {
      case FLUX_FUNCTION_ROE :
	Flux.FluxRoe_n(Wl, Wr, SolnBlk.Grid.nfaceN(i, j),
		       Input_Parameters.Preconditioning,SolnBlk.Flow_Type,delta_n);
	break;
      case FLUX_FUNCTION_HLLE :
	Flux.FluxHLLE_n(Wl, Wr, SolnBlk.Grid.nfaceN(i, j));
	break;
      case FLUX_FUNCTION_LINDE :
	Flux.FluxLinde_n(Wl, Wr, SolnBlk.Grid.nfaceN(i, j));
	break;
      case FLUX_FUNCTION_AUSM_PLUS_UP :
	Flux.FluxAUSMplus_up_n(Wl, Wr, SolnBlk.Grid.nfaceN(i, j));
	break;
      } /* endswitch */
	

	// Determine NORTH face viscous flux.
      if(SolnBlk.Flow_Type != FLOWTYPE_INVISCID){ // -ve as on RHS    
	switch(Input_Parameters.i_Viscous_Flux_Evaluation){	    
	case VISCOUS_RECONSTRUCTION_DIAMONDPATH_LEAST_SQUARES :	      
	case VISCOUS_RECONSTRUCTION_DIAMONDPATH_GREEN_GAUSS :	  	  
	  // Wn temporarily stored in Wnd -> calculated in reconstruction
	  W_face.Average(SolnBlk.Wnd[i][j+1], SolnBlk.Wnd[i+1][j+1]);
	  Flux.Viscous_Flux_n(W_face,
			      SolnBlk.dWdx_faceN[i][j],
			      SolnBlk.dWdy_faceN[i][j],
			      SolnBlk.Axisymmetric,
			      SolnBlk.Grid.xfaceN(i,j),
			      SolnBlk.Grid.nfaceN(i,j),
			      -1.0); // subtract from inviscid flux
	  break;
	case VISCOUS_RECONSTRUCTION_HYBRID :
	  W_face.Average(SolnBlk.W[i][j], SolnBlk.W[i][j+1]);
	  Flux.Viscous_FluxHybrid_n(W_face,
				    SolnBlk.dWdx_faceN[i][j],
				    SolnBlk.dWdy_faceN[i][j],
				    SolnBlk.Grid.xfaceN(i,j),
				    SolnBlk.W[i][j],
				    SolnBlk.dWdx[i][j],
				    SolnBlk.dWdy[i][j],
				    SolnBlk.Grid.Cell[i][j].Xc,
				    SolnBlk.W[i][j+1],
				    SolnBlk.dWdy[i][j+1],
				    SolnBlk.dWdy[i][j+1],
				    SolnBlk.Grid.Cell[i][j+1].Xc,
				    SolnBlk.Axisymmetric,
				    SolnBlk.Grid.nfaceN(i,j),
				    -1.0); // subtract from inviscid flux
	  SolnBlk.dWdx_faceS[i][j+1] = SolnBlk.dWdx_faceN[i][j];
	  SolnBlk.dWdy_faceS[i][j+1] = SolnBlk.dWdy_faceN[i][j];
	  break;
	}// endswitch
      }// endif - viscous

      //-----------------------------------------------------------

      // Evaluate cell-averaged solution changes.
      SolnBlk.dUdt[i][j][k_residual].add( Flux, -(Input_Parameters.CFL_Number*
						  SolnBlk.dt[i][j]*SolnBlk.Grid.lfaceN(i,j)/
						  SolnBlk.Grid.Cell[i][j].A) );

      SolnBlk.dUdt[i][j+1][k_residual].add( Flux, (Input_Parameters.CFL_Number*
						   SolnBlk.dt[i][j+1]*
						   SolnBlk.Grid.lfaceS(i,j+1)/
						   SolnBlk.Grid.Cell[i][j+1].A) );

      // Save south and north face boundary flux.
      // USED for AMR 
      if (j == SolnBlk.JCl-1) {
	SolnBlk.FluxS[i].set( Flux, -SolnBlk.Grid.lfaceS(i, j+1) );
      } else if (j == SolnBlk.JCu) {
	SolnBlk.FluxN[i].set( Flux, SolnBlk.Grid.lfaceN(i, j) );
      }
	  
    } // endfor - j

    SolnBlk.dUdt[i][SolnBlk.JCu+1][k_residual].Vacuum();
    SolnBlk.dUdt[i][SolnBlk.JCl-1][k_residual].Vacuum();

  } // endfor - i

    /* Residual successfully calculated. */

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
  const int NUM_VAR_FLAME2D( SolnBlk.NumVar() );
  static Flame2D_pState Wo;
  bool isGoodState;

  // Memory for linear system solver. 
  DenseSystemLinEqs LinSys;
  static DenseMatrix dSdU(NUM_VAR_FLAME2D,NUM_VAR_FLAME2D);            // Source Jacobian
  static DenseMatrix Precon(NUM_VAR_FLAME2D,NUM_VAR_FLAME2D,ZERO);     // Low Mach number preconditioner
  static DenseMatrix Precon_Inv(NUM_VAR_FLAME2D,NUM_VAR_FLAME2D,ZERO); // Inverse of low Mach number preconditioner
  
  /* Allocate memory for linear system solver. */
  LinSys.allocate(NUM_VAR_FLAME2D);

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


      //cout << endl << "Updating cell i=" << setw(6) << i << ", j=" << setw(6) << j << flush;

      /******************************************************/
      /********** Fully Explicit ****************************/
      /******************************************************/
      if (Input_Parameters.Local_Time_Stepping == GLOBAL_TIME_STEPPING || 
	  Input_Parameters.Local_Time_Stepping == SCALAR_LOCAL_TIME_STEPPING ) {

	//Update
	SolnBlk.U[i][j] = SolnBlk.Uo[i][j];
	SolnBlk.U[i][j].add( SolnBlk.dUdt[i][j][k_residual], omega );
	
	//Check for unphysical properties  
	/**********************************************************/
	/* If unphysical properties and using global timestepping */ 
	/* stop simulation                                        */
	/**********************************************************/   
	if (Input_Parameters.Local_Time_Stepping == GLOBAL_TIME_STEPPING) {
	  if(!SolnBlk.U[i][j].isPhysical(10)) return (i);
	  
	  /*********************************************************/
	  /* If unphysical properties and using local timestepping */ 
	  /* try reducing step size                                */
	  /*********************************************************/    
	} else if (Input_Parameters.Local_Time_Stepping == SCALAR_LOCAL_TIME_STEPPING) {
	  
	  isGoodState = SolnBlk.U[i][j].isPhysical(10);
	  if( !isGoodState ) {  
	    
	    if (SolnBlk.debug_level) {
	      cout<<"\n Using local time stepping step reduction in Update_Solution_Multistage_Explicit\n";
	    }
	    
	    double residual_reduction_factor(ONE);
	    for (int n_residual_reduction = 1; n_residual_reduction <= 10; ++n_residual_reduction) {
	      if(SolnBlk.debug_level) {
		cout<<".."<<n_residual_reduction<<"..";
	      }
	      residual_reduction_factor = HALF*residual_reduction_factor;
	      SolnBlk.dt[i][j] = residual_reduction_factor*SolnBlk.dt[i][j];
	      SolnBlk.dUdt[i][j][k_residual].set( SolnBlk.dUdt[i][j][k_residual], residual_reduction_factor );
	      
	      //Update
	      SolnBlk.U[i][j] = SolnBlk.Uo[i][j];
	      SolnBlk.U[i][j].add( SolnBlk.dUdt[i][j][k_residual], omega );
	      
	      isGoodState = SolnBlk.U[i][j].isPhysical(n_residual_reduction);
	      if(isGoodState)  break;
	    } /* endfor */  
	    if (!isGoodState) return(i);
	  }
	} /* endif */
	
	  /**************************************************************/
	  /************ SEMI-IMPLICIT AND/OR PRECONDITIONING ************/
	  /**************************************************************/
      } else if (Input_Parameters.Local_Time_Stepping == SEMI_IMPLICIT_LOCAL_TIME_STEPPING || 
		 Input_Parameters.Local_Time_Stepping == LOW_MACH_NUMBER_WEISS_SMITH_PRECONDITIONER || 
		 Input_Parameters.Local_Time_Stepping == SEMI_IMPLICIT_LOW_MACH_NUMBER_PRECONDITIONER) {
	
	/************ FORM LHS FOR SEMI-IMPLICIT METHOD *****************/
	if (Input_Parameters.Local_Time_Stepping == SEMI_IMPLICIT_LOCAL_TIME_STEPPING){

	  // get initial pState
	  Wo.setU( SolnBlk.Uo[i][j] );

	  // Semi implicit formulation: set up system of equations 
	  // and include source Jacobian in the LHS matrix. 
	  dSdU.zero();
	  SemiImplicitBlockJacobi(dSdU, SolnBlk, Wo, i, j);	  	  
	  LinSys.A.identity();
	     
	  //scalar multiplication
	  double mult( omega*Input_Parameters.CFL_Number*SolnBlk.dt[i][j] );
	  const int n = LinSys.A.get_n();
	  for (int ii=0; ii<n; ii++)
	    for (int jj=0; jj<n; jj++)
	      LinSys.A(ii,jj) -= mult*dSdU(ii,jj);	 
	
	
	  /************ FORM LHS FOR LOW MACH NUMBER PRECONDITIONING ***************/
	} else if (Input_Parameters.Local_Time_Stepping == LOW_MACH_NUMBER_WEISS_SMITH_PRECONDITIONER ){
	  
	  // get initial pState
	  Wo.setU( SolnBlk.Uo[i][j] );

	  // spacing for preconditioner
	  delta_n = SolnBlk.delta_n(i,j);

	  // build preconditioner
	  LinSys.A.zero();
	  Wo.Low_Mach_Number_Preconditioner(LinSys.A,
					    SolnBlk.Flow_Type,
					    delta_n);
	  // not Multiplied by omega*CFL_Number*SolnBlk.dt[i][j] as dimensionless


	  /******* FORM LHS FOR SEMI-IMPLICIT METHOD WITH PRECONDITIONING ********/
	  //LHS = P*( I- Pinv*(dt*(dS/dU)))	     
	} else if (Input_Parameters.Local_Time_Stepping == SEMI_IMPLICIT_LOW_MACH_NUMBER_PRECONDITIONER ){
	 	  
	  // get initial pState
	  Wo.setU( SolnBlk.Uo[i][j] );

	  // Semi implicit formulation: set up system of equations
	  // and include source Jacobian in the LHS matrix.
	  dSdU.zero();
	  SemiImplicitBlockJacobi(dSdU,SolnBlk,Wo,i, j);
	  
	  // spacing for preconditioner
	  delta_n = SolnBlk.delta_n(i,j);

	  //Inverse of Preconditioner
	  // Precon_Inv.zero(); // <- no need, always writing to the same spot
	  Wo.Low_Mach_Number_Preconditioner_Inverse(Precon_Inv,
						    SolnBlk.Flow_Type,
						    delta_n);
	  
	  //I - Pinv*(dt*(dS/dU))
	  LinSys.A.identity();
	  LinSys.A -= Precon_Inv*((omega*Input_Parameters.CFL_Number*SolnBlk.dt[i][j])*dSdU);
	  
	  //Preconditioner 
	  // Precon.zero(); // <- no need, always writing to the same spot
	  Wo.Low_Mach_Number_Preconditioner(Precon,
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
	
	/******************************************************************
	 ******************* EVALUATE RHS *********************************
	 ******************************************************************/
	// omega*(dt*RHS)
	for (int k = 0; k < NUM_VAR_FLAME2D; k++) {
	  LinSys.b(k) = omega*SolnBlk.dUdt[i][j][k_residual][k+1];
	}
	/* Solve system of equations using LU decomposition Gaussian elimination procedure. */
	LinSys.solve(LU_DECOMPOSITION);
	
	/* Update the conserved solution variables. */
	for (int k = 0; k < NUM_VAR_FLAME2D; k++) {
	  SolnBlk.U[i][j][k+1] = SolnBlk.Uo[i][j][k+1] + LinSys.x(k);
	} 


	/*********************************************************/
	/* If unphysical properties and using local timestepping */ 
	/* try reducing step size                                */
	/*********************************************************/
	
	isGoodState = SolnBlk.U[i][j].isPhysical(10);
	if(!isGoodState){
	  
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
	    SolnBlk.dUdt[i][j][k_residual].set( SolnBlk.dUdt[i][j][k_residual], residual_reduction_factor );
	    
	    if (Input_Parameters.Local_Time_Stepping == SEMI_IMPLICIT_LOCAL_TIME_STEPPING ){
	      LinSys.A.identity();
	      LinSys.A -= (omega*Input_Parameters.CFL_Number*SolnBlk.dt[i][j])*dSdU;

	    } else if (Input_Parameters.Local_Time_Stepping == SEMI_IMPLICIT_LOW_MACH_NUMBER_PRECONDITIONER){
	      LinSys.A.identity();		 
	      LinSys.A -= Precon_Inv*((residual_reduction_factor*omega*Input_Parameters.CFL_Number*SolnBlk.dt[i][j])*dSdU);
	      //LHS = P*( I - Pinv*(dt*(dS/dU)))	     
	      LinSys.A = Precon*LinSys.A;	     
	    }

	    for (int k = 0; k < NUM_VAR_FLAME2D; k++) {
	      LinSys.b(k) = omega*SolnBlk.dUdt[i][j][k_residual][k+1];
	    } 
	    
	    LinSys.solve(LU_DECOMPOSITION);
	    
	    /* Update the conserved Solution */
	    for (int k = 0; k < NUM_VAR_FLAME2D; k++) {
	      SolnBlk.U[i][j][k+1] = SolnBlk.Uo[i][j][k+1] + LinSys.x(k);
	    } 
	    	    
	    isGoodState = SolnBlk.U[i][j].isPhysical(n_residual_reduction);
	    if(isGoodState)  break;
	  }
	} 
	/*********************************************************/
	//Check for unphysical properties
	if(!isGoodState) return (i);
	
      } //end implicit and/or preconditioned formulations 
    
	//Calculate Primitive values from updated conserved solution      
      SolnBlk.W[i][j].setU(SolnBlk.U[i][j]);
       
    } //end i
  } //end j //


  LinSys.deallocate();

  /* Solution successfully updated. */
  return (0);
  
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


} // end Radiation_Source_Eval()
