/* Gaussian2DQuadSingleBlock.cc:  Single-Block Versions of Subroutines for 2D Gaussian
                                  Multi-Block Quadrilateral Mesh 
                                  Solution Classes. */

/* Include 2D Gaussian quadrilateral mesh solution header file. */

#ifndef _GAUSSIAN2D_QUAD_INCLUDED
#include "Gaussian2DQuad.h"
#endif // _GAUSSIAN2D_QUAD_INCLUDED

/**************************************************************************
 * Gaussian2D_Quad_Block -- Single Block Static Variables.                *
 **************************************************************************/
int Gaussian2D_Quad_Block::Flow_Type = FLOWTYPE_INVISCID;
int Gaussian2D_Quad_Block::Heat_Transfer = 0;
int Gaussian2D_Quad_Block::residual_variable = 2;  //x-momentum

/**************************************************************************
 * Gaussian2D_Quad_Block -- Single Block External Subroutines.            *
 **************************************************************************/

/********************************************************
 * Routine: Broadcast_Solution_Block                    *
 *                                                      *
 * Broadcast quadrilateral solution block to all        *
 * processors involved in the calculation from the      *
 * primary processor using the MPI broadcast routine.   *
 *                                                      *
 ********************************************************/
void Broadcast_Solution_Block(Gaussian2D_Quad_Block &SolnBlk) {

#ifdef _MPI_VERSION
  int i, j, ni, nj, ng, nr, heat, block_allocated, buffer_size;
  double *buffer;

  /* Broadcast the number of cells in each direction. */

  if (CFFC_Primary_MPI_Processor()) {
    ni = SolnBlk.NCi;
    nj = SolnBlk.NCj;
    ng = SolnBlk.Nghost;
    nr = SolnBlk.residual_variable;
    heat = SolnBlk.Heat_Transfer;
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
  MPI::COMM_WORLD.Bcast(&heat,1,MPI::INT,0);
  MPI::COMM_WORLD.Bcast(&block_allocated, 1, MPI::INT, 0);

  /* On non-primary MPI processors, allocate (re-allocate) 
     memory for the quadrilateral solution block as necessary. */

  if (!CFFC_Primary_MPI_Processor()) {
    if (SolnBlk.NCi != ni || SolnBlk.NCj != nj || SolnBlk.Nghost != ng) { 
      if (SolnBlk.U != NULL) SolnBlk.deallocate();
      if (block_allocated) SolnBlk.allocate(ni-2*ng, nj-2*ng, ng); 
    } /* endif */
    if (SolnBlk.residual_variable != nr) SolnBlk.residual_variable = nr;  //do I need to do this?  I copied it from Jai
    if (SolnBlk.Heat_Transfer != heat) SolnBlk.Heat_Transfer = heat;   
  } /* endif */

    /* Broadcast the axisymmetric/planar flow indicator. */

  MPI::COMM_WORLD.Bcast(&(SolnBlk.Axisymmetric), 1, MPI::INT, 0);

  /* Broadcast the grid. */

  Broadcast_Quad_Block(SolnBlk.Grid);

  /* Broadcast the solution state variables. */

  if (block_allocated) {
    ni = (SolnBlk.ICu+SolnBlk.Nghost) - (SolnBlk.ICl-SolnBlk.Nghost) + 1;
    nj = (SolnBlk.JCu+SolnBlk.Nghost) - (SolnBlk.JCl-SolnBlk.Nghost) + 1;
    buffer = new double[8*ni*nj];

    if (CFFC_Primary_MPI_Processor()) {
      buffer_size = 0;
      for (j  = SolnBlk.JCl-SolnBlk.Nghost ; j <= SolnBlk.JCu+SolnBlk.Nghost ; ++j ) {
	for ( i = SolnBlk.ICl-SolnBlk.Nghost ; i <= SolnBlk.ICu+SolnBlk.Nghost ; ++i ) {
	  buffer[buffer_size  ] = SolnBlk.U[i][j].d;
	  buffer[buffer_size+1] = SolnBlk.U[i][j].dv.x;
	  buffer[buffer_size+2] = SolnBlk.U[i][j].dv.y;
	  buffer[buffer_size+3] = SolnBlk.U[i][j].E.xx;
	  buffer[buffer_size+4] = SolnBlk.U[i][j].E.xy;
	  buffer[buffer_size+5] = SolnBlk.U[i][j].E.yy;
	  buffer[buffer_size+6] = SolnBlk.U[i][j].E.zz;
	  buffer[buffer_size+7] = SolnBlk.U[i][j].erot;
	  buffer_size = buffer_size + 8;
	} /* endfor */
      } /* endfor */
    } /* endif */

    buffer_size = 8*ni*nj;
    MPI::COMM_WORLD.Bcast(buffer, buffer_size, MPI::DOUBLE, 0);

    if (!CFFC_Primary_MPI_Processor()) {
      buffer_size = 0;
      for (j  = SolnBlk.JCl-SolnBlk.Nghost ; j <= SolnBlk.JCu+SolnBlk.Nghost ; ++j ) {
	for ( i = SolnBlk.ICl-SolnBlk.Nghost ; i <= SolnBlk.ICu+SolnBlk.Nghost ; ++i ) {
	  SolnBlk.U[i][j].d    = buffer[buffer_size];
	  SolnBlk.U[i][j].dv.x = buffer[buffer_size+1];
	  SolnBlk.U[i][j].dv.y = buffer[buffer_size+2];
	  SolnBlk.U[i][j].E.xx = buffer[buffer_size+3];
	  SolnBlk.U[i][j].E.xy = buffer[buffer_size+4];
	  SolnBlk.U[i][j].E.yy = buffer[buffer_size+5];
	  SolnBlk.U[i][j].E.zz = buffer[buffer_size+6];
	  SolnBlk.U[i][j].erot = buffer[buffer_size+7];
	  buffer_size = buffer_size + 8;
	  SolnBlk.W[i][j] = W(SolnBlk.U[i][j]);
	} /* endfor */
      } /* endfor */
    } /* endif */

    delete []buffer; 
    buffer = NULL;

    ni = 1;
    nj = (SolnBlk.JCu+SolnBlk.Nghost) - (SolnBlk.JCl-SolnBlk.Nghost) + 1;
    buffer = new double[18*ni*nj];

    if (CFFC_Primary_MPI_Processor()) {
      buffer_size = 0;
      for (j  = SolnBlk.JCl-SolnBlk.Nghost ; j <= SolnBlk.JCu+SolnBlk.Nghost ; ++j ) {
	buffer[buffer_size  ]  = SolnBlk.WoW[j].d;
	buffer[buffer_size+1]  = SolnBlk.WoW[j].v.x;
	buffer[buffer_size+2]  = SolnBlk.WoW[j].v.y;
	buffer[buffer_size+3]  = SolnBlk.WoW[j].p.xx;
	buffer[buffer_size+4]  = SolnBlk.WoW[j].p.xy;
	buffer[buffer_size+5]  = SolnBlk.WoW[j].p.yy;
	buffer[buffer_size+6]  = SolnBlk.WoW[j].p.zz;
	buffer[buffer_size+7]  = SolnBlk.WoW[j].erot;
	buffer[buffer_size+8]  = SolnBlk.WoE[j].d;
	buffer[buffer_size+9]  = SolnBlk.WoE[j].v.x;
	buffer[buffer_size+10] = SolnBlk.WoE[j].v.y;
	buffer[buffer_size+11] = SolnBlk.WoE[j].p.xx;
	buffer[buffer_size+12] = SolnBlk.WoE[j].p.xy;
	buffer[buffer_size+13] = SolnBlk.WoE[j].p.yy;
	buffer[buffer_size+14] = SolnBlk.WoE[j].p.zz;
	buffer[buffer_size+15] = SolnBlk.WoE[j].erot;
	buffer[buffer_size+16] = SolnBlk.oldT_W[j];
	buffer[buffer_size+17] = SolnBlk.oldT_E[j];
	buffer_size = buffer_size + 18;
      } /* endfor */
    } /* endif */

    buffer_size = 18*ni*nj;
    MPI::COMM_WORLD.Bcast(buffer, buffer_size, MPI::DOUBLE, 0);

    if (!CFFC_Primary_MPI_Processor()) {
      buffer_size = 0;
      for (j  = SolnBlk.JCl-SolnBlk.Nghost ; j <= SolnBlk.JCu+SolnBlk.Nghost ; ++j ) {
	SolnBlk.WoW[j].d    = buffer[buffer_size];
	SolnBlk.WoW[j].v.x  = buffer[buffer_size+1];
	SolnBlk.WoW[j].v.y  = buffer[buffer_size+2];
	SolnBlk.WoW[j].p.xx = buffer[buffer_size+3];
	SolnBlk.WoW[j].p.xy = buffer[buffer_size+4];
	SolnBlk.WoW[j].p.yy = buffer[buffer_size+5];
	SolnBlk.WoW[j].p.zz = buffer[buffer_size+6];
	SolnBlk.WoW[j].erot = buffer[buffer_size+7];
	SolnBlk.WoE[j].d    = buffer[buffer_size+8];
	SolnBlk.WoE[j].v.x  = buffer[buffer_size+9];
	SolnBlk.WoE[j].v.y  = buffer[buffer_size+10];
	SolnBlk.WoE[j].p.xx = buffer[buffer_size+11];
	SolnBlk.WoE[j].p.xy = buffer[buffer_size+12];
	SolnBlk.WoE[j].p.yy = buffer[buffer_size+13];
	SolnBlk.WoE[j].p.zz = buffer[buffer_size+14];
	SolnBlk.WoE[j].erot = buffer[buffer_size+15];
	SolnBlk.oldT_W[j]   = buffer[buffer_size+16];
	SolnBlk.oldT_E[j]   = buffer[buffer_size+17];
	buffer_size = buffer_size + 18;
      } /* endfor */
    } /* endif */

    delete []buffer; 
    buffer = NULL;

    ni = (SolnBlk.ICu+SolnBlk.Nghost) - (SolnBlk.ICl-SolnBlk.Nghost) + 1;
    nj = 1;
    buffer = new double[18*ni*nj];

    if (CFFC_Primary_MPI_Processor()) {
      buffer_size = 0;
      for (i  = SolnBlk.ICl-SolnBlk.Nghost ; i <= SolnBlk.ICu+SolnBlk.Nghost ; ++i ) {
	buffer[buffer_size  ]  = SolnBlk.WoS[i].d;
	buffer[buffer_size+1]  = SolnBlk.WoS[i].v.x;
	buffer[buffer_size+2]  = SolnBlk.WoS[i].v.y;
	buffer[buffer_size+3]  = SolnBlk.WoS[i].p.xx;
	buffer[buffer_size+4]  = SolnBlk.WoS[i].p.xy;
	buffer[buffer_size+5]  = SolnBlk.WoS[i].p.yy;
	buffer[buffer_size+6]  = SolnBlk.WoS[i].p.zz;
	buffer[buffer_size+7]  = SolnBlk.WoS[i].erot;
	buffer[buffer_size+8]  = SolnBlk.WoN[i].d;
	buffer[buffer_size+9]  = SolnBlk.WoN[i].v.x;
	buffer[buffer_size+10] = SolnBlk.WoN[i].v.y;
	buffer[buffer_size+11] = SolnBlk.WoN[i].p.xx;
	buffer[buffer_size+12] = SolnBlk.WoN[i].p.xy;
	buffer[buffer_size+13] = SolnBlk.WoN[i].p.yy;
	buffer[buffer_size+14] = SolnBlk.WoN[i].p.zz;
	buffer[buffer_size+15] = SolnBlk.WoN[i].erot;
	buffer[buffer_size+16] = SolnBlk.oldT_S[i];
	buffer[buffer_size+17] = SolnBlk.oldT_N[i];
	buffer_size = buffer_size + 18;
      } /* endfor */
    } /* endif */

    buffer_size = 18*ni*nj;
    MPI::COMM_WORLD.Bcast(buffer, buffer_size, MPI::DOUBLE, 0);

    if (!CFFC_Primary_MPI_Processor()) {
      buffer_size = 0;
      for (i  = SolnBlk.ICl-SolnBlk.Nghost ; i <= SolnBlk.ICu+SolnBlk.Nghost ; ++i ) {
	SolnBlk.WoS[i].d    = buffer[buffer_size];
	SolnBlk.WoS[i].v.x  = buffer[buffer_size+1];
	SolnBlk.WoS[i].v.y  = buffer[buffer_size+2];
	SolnBlk.WoS[i].p.xx = buffer[buffer_size+3];
	SolnBlk.WoS[i].p.xy = buffer[buffer_size+4];
	SolnBlk.WoS[i].p.yy = buffer[buffer_size+5];
	SolnBlk.WoS[i].p.zz = buffer[buffer_size+6];
	SolnBlk.WoS[i].erot = buffer[buffer_size+7];
	SolnBlk.WoN[i].d    = buffer[buffer_size+8];
	SolnBlk.WoN[i].v.x  = buffer[buffer_size+9];
	SolnBlk.WoN[i].v.y  = buffer[buffer_size+10];
	SolnBlk.WoN[i].p.xx = buffer[buffer_size+11];
	SolnBlk.WoN[i].p.xy = buffer[buffer_size+12];
	SolnBlk.WoN[i].p.yy = buffer[buffer_size+13];
	SolnBlk.WoN[i].p.zz = buffer[buffer_size+14];
	SolnBlk.WoN[i].erot = buffer[buffer_size+15];
	SolnBlk.oldT_S[i]   = buffer[buffer_size+16];
	SolnBlk.oldT_N[i]   = buffer[buffer_size+17];
	buffer_size = buffer_size + 18;
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
void Broadcast_Solution_Block(Gaussian2D_Quad_Block &SolnBlk,
                              MPI::Intracomm &Communicator, 
                              const int Source_CPU) {

  int Source_Rank = 0;
  int i, j, ni, nj, ng, nr, heat, block_allocated, buffer_size;
  double *buffer;

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
    nr = SolnBlk.residual_variable;
    heat = SolnBlk.Heat_Transfer;
  } /* endif */

  Communicator.Bcast(&ni, 1, MPI::INT, Source_Rank);
  Communicator.Bcast(&nj, 1, MPI::INT, Source_Rank);
  Communicator.Bcast(&ng, 1, MPI::INT, Source_Rank);
  Communicator.Bcast(&nr,1,MPI::INT,Source_Rank);
  Communicator.Bcast(&heat,1,MPI::INT,Source_Rank);
  Communicator.Bcast(&block_allocated, 1, MPI::INT, Source_Rank);

  /* On non-source MPI processors, allocate (re-allocate) 
     memory for the quadrilateral solution block as necessary. */

  if (!(CFFC_MPI::This_Processor_Number == Source_CPU)) {
    if (SolnBlk.NCi != ni || SolnBlk.NCj != nj || SolnBlk.Nghost != ng) { 
      if (SolnBlk.U != NULL) SolnBlk.deallocate();
      if (block_allocated) SolnBlk.allocate(ni-2*ng, nj-2*ng, ng); 
    } /* endif */
    // Set the block static variables if they were not previously assigned.
    if (SolnBlk.residual_variable != nr) SolnBlk.residual_variable = nr; //do I need to do this?  I copied it from Jai
    if (SolnBlk.Heat_Transfer != heat) SolnBlk.Heat_Transfer = heat;
  } /* endif */

    /* Broadcast the axisymmetric/planar flow indicator. */

  Communicator.Bcast(&(SolnBlk.Axisymmetric), 1, MPI::INT, Source_Rank);

  /* Broadcast the grid. */

  Broadcast_Quad_Block(SolnBlk.Grid, Communicator, Source_CPU);

  /* Broadcast the solution state variables. */

  if (block_allocated) {
    ni = (SolnBlk.ICu+SolnBlk.Nghost) - (SolnBlk.ICl-SolnBlk.Nghost) + 1;
    nj = (SolnBlk.JCu+SolnBlk.Nghost) - (SolnBlk.JCl-SolnBlk.Nghost) + 1;
    buffer = new double[8*ni*nj];

    if (CFFC_MPI::This_Processor_Number == Source_CPU) {
      buffer_size = 0;
      for (j  = SolnBlk.JCl-SolnBlk.Nghost ; j <= SolnBlk.JCu+SolnBlk.Nghost ; ++j ) {
	for ( i = SolnBlk.ICl-SolnBlk.Nghost ; i <= SolnBlk.ICu+SolnBlk.Nghost ; ++i ) {
	  buffer[buffer_size  ] = SolnBlk.U[i][j].d;
	  buffer[buffer_size+1] = SolnBlk.U[i][j].dv.x;
	  buffer[buffer_size+2] = SolnBlk.U[i][j].dv.y;
	  buffer[buffer_size+3] = SolnBlk.U[i][j].E.xx;
	  buffer[buffer_size+4] = SolnBlk.U[i][j].E.xy;
	  buffer[buffer_size+5] = SolnBlk.U[i][j].E.yy;
	  buffer[buffer_size+6] = SolnBlk.U[i][j].E.zz;
	  buffer[buffer_size+7] = SolnBlk.U[i][j].erot;
	  buffer_size = buffer_size + 8;
	} /* endfor */
      } /* endfor */
    } /* endif */

    buffer_size = 8*ni*nj;
    Communicator.Bcast(buffer, buffer_size, MPI::DOUBLE, Source_Rank);

    if (!(CFFC_MPI::This_Processor_Number == Source_CPU)) {
      buffer_size = 0;
      for (j  = SolnBlk.JCl-SolnBlk.Nghost ; j <= SolnBlk.JCu+SolnBlk.Nghost ; ++j ) {
	for ( i = SolnBlk.ICl-SolnBlk.Nghost ; i <= SolnBlk.ICu+SolnBlk.Nghost ; ++i ) {
	  SolnBlk.U[i][j].d    = buffer[buffer_size];
	  SolnBlk.U[i][j].dv.x = buffer[buffer_size+1];
	  SolnBlk.U[i][j].dv.y = buffer[buffer_size+2];
	  SolnBlk.U[i][j].E.xx = buffer[buffer_size+3];
	  SolnBlk.U[i][j].E.xy = buffer[buffer_size+4];
	  SolnBlk.U[i][j].E.yy = buffer[buffer_size+5];
	  SolnBlk.U[i][j].E.zz = buffer[buffer_size+6];
	  SolnBlk.U[i][j].erot = buffer[buffer_size+7];
	  buffer_size = buffer_size + 8;
	  SolnBlk.W[i][j] = W(SolnBlk.U[i][j]);
	} /* endfor */
      } /* endfor */
    } /* endif */

    delete []buffer; 
    buffer = NULL;

    ni = 1;
    nj = (SolnBlk.JCu+SolnBlk.Nghost) - (SolnBlk.JCl-SolnBlk.Nghost) + 1;
    buffer = new double[18*ni*nj];

    if (CFFC_MPI::This_Processor_Number == Source_CPU) {
      buffer_size = 0;
      for (j  = SolnBlk.JCl-SolnBlk.Nghost ; j <= SolnBlk.JCu+SolnBlk.Nghost ; ++j ) {
	buffer[buffer_size  ]  = SolnBlk.WoW[j].d;
	buffer[buffer_size+1]  = SolnBlk.WoW[j].v.x;
	buffer[buffer_size+2]  = SolnBlk.WoW[j].v.y;
	buffer[buffer_size+3]  = SolnBlk.WoW[j].p.xx;
	buffer[buffer_size+4]  = SolnBlk.WoW[j].p.xy;
	buffer[buffer_size+5]  = SolnBlk.WoW[j].p.yy;
	buffer[buffer_size+6]  = SolnBlk.WoW[j].p.zz;
	buffer[buffer_size+7]  = SolnBlk.WoW[j].erot;
	buffer[buffer_size+8]  = SolnBlk.WoE[j].d;
	buffer[buffer_size+9]  = SolnBlk.WoE[j].v.x;
	buffer[buffer_size+10] = SolnBlk.WoE[j].v.y;
	buffer[buffer_size+11] = SolnBlk.WoE[j].p.xx;
	buffer[buffer_size+12] = SolnBlk.WoE[j].p.xy;
	buffer[buffer_size+13] = SolnBlk.WoE[j].p.yy;
	buffer[buffer_size+14] = SolnBlk.WoE[j].p.zz;
	buffer[buffer_size+15] = SolnBlk.WoE[j].erot;
	buffer[buffer_size+16] = SolnBlk.oldT_W[j];
	buffer[buffer_size+17] = SolnBlk.oldT_E[j];
	buffer_size = buffer_size + 18;
      } /* endfor */
    } /* endif */

    buffer_size = 18*ni*nj;
    Communicator.Bcast(buffer, buffer_size, MPI::DOUBLE, Source_Rank);

    if (!(CFFC_MPI::This_Processor_Number == Source_CPU)) {
      buffer_size = 0;
      for (j  = SolnBlk.JCl-SolnBlk.Nghost ; j <= SolnBlk.JCu+SolnBlk.Nghost ; ++j ) {
	SolnBlk.WoW[j].d    = buffer[buffer_size];
	SolnBlk.WoW[j].v.x  = buffer[buffer_size+1];
	SolnBlk.WoW[j].v.y  = buffer[buffer_size+2];
	SolnBlk.WoW[j].p.xx = buffer[buffer_size+3];
	SolnBlk.WoW[j].p.xy = buffer[buffer_size+4];
	SolnBlk.WoW[j].p.yy = buffer[buffer_size+5];
	SolnBlk.WoW[j].p.zz = buffer[buffer_size+6];
	SolnBlk.WoW[j].erot = buffer[buffer_size+7];
	SolnBlk.WoE[j].d    = buffer[buffer_size+8];
	SolnBlk.WoE[j].v.x  = buffer[buffer_size+9];
	SolnBlk.WoE[j].v.y  = buffer[buffer_size+10];
	SolnBlk.WoE[j].p.xx = buffer[buffer_size+11];
	SolnBlk.WoE[j].p.xy = buffer[buffer_size+12];
	SolnBlk.WoE[j].p.yy = buffer[buffer_size+13];
	SolnBlk.WoE[j].p.zz = buffer[buffer_size+14];
	SolnBlk.WoE[j].erot = buffer[buffer_size+15];
	SolnBlk.oldT_W[j]   = buffer[buffer_size+16];
	SolnBlk.oldT_E[j]   = buffer[buffer_size+17];
	buffer_size = buffer_size + 18;
      } /* endfor */
    } /* endif */

    delete []buffer; 
    buffer = NULL;

    ni = (SolnBlk.ICu+SolnBlk.Nghost) - (SolnBlk.ICl-SolnBlk.Nghost) + 1;
    nj = 1;
    buffer = new double[18*ni*nj];

    if (CFFC_MPI::This_Processor_Number == Source_CPU) {
      buffer_size = 0;
      for (i  = SolnBlk.ICl-SolnBlk.Nghost ; i <= SolnBlk.ICu+SolnBlk.Nghost ; ++i ) {
	buffer[buffer_size  ]  = SolnBlk.WoS[i].d;
	buffer[buffer_size+1]  = SolnBlk.WoS[i].v.x;
	buffer[buffer_size+2]  = SolnBlk.WoS[i].v.y;
	buffer[buffer_size+3]  = SolnBlk.WoS[i].p.xx;
	buffer[buffer_size+4]  = SolnBlk.WoS[i].p.xy;
	buffer[buffer_size+5]  = SolnBlk.WoS[i].p.yy;
	buffer[buffer_size+6]  = SolnBlk.WoS[i].p.zz;
	buffer[buffer_size+7]  = SolnBlk.WoS[i].erot;
	buffer[buffer_size+8]  = SolnBlk.WoN[i].d;
	buffer[buffer_size+9]  = SolnBlk.WoN[i].v.x;
	buffer[buffer_size+10] = SolnBlk.WoN[i].v.y;
	buffer[buffer_size+11] = SolnBlk.WoN[i].p.xx;
	buffer[buffer_size+12] = SolnBlk.WoN[i].p.xy;
	buffer[buffer_size+13] = SolnBlk.WoN[i].p.yy;
	buffer[buffer_size+14] = SolnBlk.WoN[i].p.zz;
	buffer[buffer_size+15] = SolnBlk.WoN[i].erot;
	buffer[buffer_size+16] = SolnBlk.oldT_S[i];
	buffer[buffer_size+17] = SolnBlk.oldT_N[i];
	buffer_size = buffer_size + 18;
      } /* endfor */
    } /* endif */

    buffer_size = 18*ni*nj;
    Communicator.Bcast(buffer, buffer_size, MPI::DOUBLE, Source_Rank);

    if (!(CFFC_MPI::This_Processor_Number == Source_CPU)) {
      buffer_size = 0;
      for (i  = SolnBlk.ICl-SolnBlk.Nghost ; i <= SolnBlk.ICu+SolnBlk.Nghost ; ++i ) {
	SolnBlk.WoS[i].d    = buffer[buffer_size];
	SolnBlk.WoS[i].v.x  = buffer[buffer_size+1];
	SolnBlk.WoS[i].v.y  = buffer[buffer_size+2];
	SolnBlk.WoS[i].p.xx = buffer[buffer_size+3];
	SolnBlk.WoS[i].p.xy = buffer[buffer_size+4];
	SolnBlk.WoS[i].p.yy = buffer[buffer_size+5];
	SolnBlk.WoS[i].p.zz = buffer[buffer_size+6];
	SolnBlk.WoS[i].erot = buffer[buffer_size+7];
	SolnBlk.WoN[i].d    = buffer[buffer_size+8];
	SolnBlk.WoN[i].v.x  = buffer[buffer_size+9];
	SolnBlk.WoN[i].v.y  = buffer[buffer_size+10];
	SolnBlk.WoN[i].p.xx = buffer[buffer_size+11];
	SolnBlk.WoN[i].p.xy = buffer[buffer_size+12];
	SolnBlk.WoN[i].p.yy = buffer[buffer_size+13];
	SolnBlk.WoN[i].p.zz = buffer[buffer_size+14];
	SolnBlk.WoN[i].erot = buffer[buffer_size+15];
	SolnBlk.oldT_S[i]   = buffer[buffer_size+16];
	SolnBlk.oldT_N[i]   = buffer[buffer_size+17];
	buffer_size = buffer_size + 18;
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
void Copy_Solution_Block(Gaussian2D_Quad_Block &SolnBlk1,
                         Gaussian2D_Quad_Block &SolnBlk2) {

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

    /* Copy the grid of the second solution block
       to the first solution block. */

    Copy_Quad_Block(SolnBlk1.Grid, SolnBlk2.Grid);

    /* Copy the solution information from SolnBlk2 to SolnBlk1. */

    if (SolnBlk2.U != NULL) {
       for ( j  = SolnBlk1.JCl-SolnBlk1.Nghost ; j <= SolnBlk1.JCu+SolnBlk1.Nghost ; ++j ) {
          for ( i = SolnBlk1.ICl-SolnBlk1.Nghost ; i <= SolnBlk1.ICu+SolnBlk1.Nghost ; ++i ) {
             SolnBlk1.U[i][j] = SolnBlk2.U[i][j];
             SolnBlk1.W[i][j] = SolnBlk2.W[i][j];
             for ( k = 0 ; k <= NUMBER_OF_RESIDUAL_VECTORS_GAUSSIAN2D-1 ; ++k ) {
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
       for (j  = SolnBlk1.JCl-SolnBlk1.Nghost ; j <= SolnBlk1.JCu+SolnBlk1.Nghost ; ++j ) {
	   SolnBlk1.oldT_W[j] = SolnBlk2.oldT_W[j];
           SolnBlk1.oldT_E[j] = SolnBlk2.oldT_E[j];
       } /* endfor */

       for ( i = SolnBlk1.ICl-SolnBlk1.Nghost ; i <= SolnBlk1.ICu+SolnBlk1.Nghost ; ++i ) {
           SolnBlk1.oldT_S[i] = SolnBlk2.oldT_S[i];
           SolnBlk1.oldT_N[i] = SolnBlk2.oldT_N[i];
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
int Prolong_Solution_Block(Gaussian2D_Quad_Block &SolnBlk_Fine,
		            Gaussian2D_Quad_Block &SolnBlk_Original,
                            const int Sector) {

    int i, j, i_min, i_max, j_min, j_max, mesh_refinement_permitted;
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
          // In this case, create the refined mesh.
          Refine_Mesh(SolnBlk_Fine.Grid, 
                      SolnBlk_Original.Grid,
                      Sector);
       } /* endif */
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
		   =  SolnBlk_Original.U[i][j];
	           //= (SolnBlk_Original.Grid.Cell[i][j].A/area_total_fine)*SolnBlk_Original.U[i][j];
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

       for ( j  = j_min-SolnBlk_Original.Nghost/2; j <= j_max+SolnBlk_Original.Nghost/2 ; ++j ) {
           SolnBlk_Fine.oldT_W[2*(j-j_min)+SolnBlk_Fine.JCl  ]
              = SolnBlk_Original.oldT_W[j];
           SolnBlk_Fine.oldT_W[2*(j-j_min)+SolnBlk_Fine.JCl+1]
              = SolnBlk_Original.oldT_W[j];

           SolnBlk_Fine.oldT_E[2*(j-j_min)+SolnBlk_Fine.JCl  ]
              = SolnBlk_Original.oldT_E[j];
           SolnBlk_Fine.oldT_E[2*(j-j_min)+SolnBlk_Fine.JCl+1]
              = SolnBlk_Original.oldT_E[j];
       } /* endfor */

       for ( i = i_min-SolnBlk_Original.Nghost/2 ; i <= i_max+SolnBlk_Original.Nghost/2 ; ++i ) {
           SolnBlk_Fine.oldT_S[2*(i-i_min)+SolnBlk_Fine.ICl  ]
              = SolnBlk_Original.oldT_S[i];
           SolnBlk_Fine.oldT_S[2*(i-i_min)+SolnBlk_Fine.ICl+1]
              = SolnBlk_Original.oldT_S[i];

           SolnBlk_Fine.oldT_N[2*(i-i_min)+SolnBlk_Fine.ICl  ]
              = SolnBlk_Original.oldT_N[i];
           SolnBlk_Fine.oldT_N[2*(i-i_min)+SolnBlk_Fine.ICl+1]
              = SolnBlk_Original.oldT_N[i];
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
int Restrict_Solution_Block(Gaussian2D_Quad_Block &SolnBlk_Coarse,
		             Gaussian2D_Quad_Block &SolnBlk_Original_SW,
                             Gaussian2D_Quad_Block &SolnBlk_Original_SE,
                             Gaussian2D_Quad_Block &SolnBlk_Original_NW,
                             Gaussian2D_Quad_Block &SolnBlk_Original_NE) {

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
            j <= SolnBlk_Original_SW.JCu; j += 2 ) {
      	  j_coarse = (j-SolnBlk_Original_SW.JCl)/2+SolnBlk_Coarse.JCl;
          SolnBlk_Coarse.WoW[j_coarse] = SolnBlk_Original_SW.WoW[j];
          SolnBlk_Coarse.oldT_W[j_coarse] = SolnBlk_Original_SW.oldT_W[j];
          if (j == SolnBlk_Original_SW.JCl-SolnBlk_Original_SW.Nghost) {
             SolnBlk_Coarse.WoW[j_coarse-1] = SolnBlk_Original_SW.WoW[j];
             SolnBlk_Coarse.oldT_W[j_coarse-1] = SolnBlk_Original_SW.oldT_W[j];
          } /* endif */
      } /* endfor */

      for ( i = SolnBlk_Original_SW.ICl-SolnBlk_Original_SW.Nghost ; 
            i <= SolnBlk_Original_SW.ICu; i += 2) {
      	  i_coarse = (i-SolnBlk_Original_SW.ICl)/2+SolnBlk_Coarse.ICl;
          SolnBlk_Coarse.WoS[i_coarse] = SolnBlk_Original_SW.WoS[i];
          SolnBlk_Coarse.oldT_S[i_coarse] = SolnBlk_Original_SW.oldT_S[i];
          if (i == SolnBlk_Original_SW.ICl-SolnBlk_Original_SW.Nghost) {
             SolnBlk_Coarse.WoS[i_coarse-1] = SolnBlk_Original_SW.WoS[i];
             SolnBlk_Coarse.oldT_S[i_coarse-1] = SolnBlk_Original_SW.oldT_S[i];
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
            j <= SolnBlk_Original_SE.JCu; j += 2 ) {
	  j_coarse = (j-SolnBlk_Original_SE.JCl)/2+SolnBlk_Coarse.JCl;
          SolnBlk_Coarse.WoE[j_coarse] = SolnBlk_Original_SE.WoE[j];
          SolnBlk_Coarse.oldT_E[j_coarse] = SolnBlk_Original_SE.oldT_E[j];
          if (j == SolnBlk_Original_SE.JCl-SolnBlk_Original_SE.Nghost) {
             SolnBlk_Coarse.WoE[j_coarse-1] = SolnBlk_Original_SE.WoE[j];
             SolnBlk_Coarse.oldT_E[j_coarse-1] = SolnBlk_Original_SE.oldT_E[j];
          } /* endif */
      } /* endfor */

      for ( i = SolnBlk_Original_SE.ICl;
            i <= SolnBlk_Original_SE.ICu+SolnBlk_Original_SE.Nghost ; i += 2) {
     	  i_coarse = (i-SolnBlk_Original_SE.ICl)/2+
                     (SolnBlk_Coarse.ICu-SolnBlk_Coarse.ICl+1)/2+SolnBlk_Coarse.ICl;
          SolnBlk_Coarse.WoS[i_coarse] = SolnBlk_Original_SE.WoS[i];
          SolnBlk_Coarse.oldT_S[i_coarse] = SolnBlk_Original_SE.oldT_S[i];
          if (i == SolnBlk_Original_SE.ICu+SolnBlk_Original_SE.Nghost-1) {
             SolnBlk_Coarse.WoS[i_coarse+1] = SolnBlk_Original_SE.WoS[i];
             SolnBlk_Coarse.oldT_S[i_coarse+1] = SolnBlk_Original_SE.oldT_S[i];
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

      for ( j = SolnBlk_Original_NW.JCl;
            j <= SolnBlk_Original_NW.JCu+SolnBlk_Original_NW.Nghost ; j += 2 ) {
      	  j_coarse = (j-SolnBlk_Original_NW.JCl)/2+
                     (SolnBlk_Coarse.JCu-SolnBlk_Coarse.JCl+1)/2+SolnBlk_Coarse.JCl;
          SolnBlk_Coarse.WoW[j_coarse] = SolnBlk_Original_NW.WoW[j];
          SolnBlk_Coarse.oldT_W[j_coarse] = SolnBlk_Original_NW.oldT_W[j];
          if (j == SolnBlk_Original_NW.JCu+SolnBlk_Original_NW.Nghost-1) {
             SolnBlk_Coarse.WoW[j_coarse+1] = SolnBlk_Original_NW.WoW[j];
             SolnBlk_Coarse.oldT_W[j_coarse+1] = SolnBlk_Original_NW.oldT_W[j];
          } /* endif */
      } /* endfor */

      for ( i = SolnBlk_Original_NW.ICl-SolnBlk_Original_NW.Nghost ; 
            i <= SolnBlk_Original_NW.ICu; i += 2) {
      	  i_coarse = (i-SolnBlk_Original_NW.ICl)/2+SolnBlk_Coarse.ICl;
          SolnBlk_Coarse.WoN[i_coarse] = SolnBlk_Original_NW.WoN[i];
          SolnBlk_Coarse.oldT_N[i_coarse] = SolnBlk_Original_NW.oldT_N[i];
          if (i == SolnBlk_Original_NW.ICl-SolnBlk_Original_NW.Nghost) {
             SolnBlk_Coarse.WoN[i_coarse-1] = SolnBlk_Original_NW.WoN[i];
             SolnBlk_Coarse.oldT_N[i_coarse-1] = SolnBlk_Original_NW.oldT_N[i];
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
	                                             //SolnBlk_Coarse.Grid.Cell[i_coarse][j_coarse].A;
             SolnBlk_Coarse.W[i_coarse][j_coarse] = W(SolnBlk_Coarse.U[i_coarse][j_coarse]);
          } /* endfor */
      } /* endfor */

      for ( j = SolnBlk_Original_NE.JCl;
            j <= SolnBlk_Original_NE.JCu+SolnBlk_Original_NE.Nghost ; j += 2 ) {
      	  j_coarse = (j-SolnBlk_Original_NE.JCl)/2+
                     (SolnBlk_Coarse.JCu-SolnBlk_Coarse.JCl+1)/2+SolnBlk_Coarse.JCl;
          SolnBlk_Coarse.WoE[j_coarse] = SolnBlk_Original_NE.WoE[j];
          SolnBlk_Coarse.oldT_E[j_coarse] = SolnBlk_Original_NE.oldT_E[j];
          if (j == SolnBlk_Original_NE.JCu+SolnBlk_Original_NE.Nghost-1) {
             SolnBlk_Coarse.WoE[j_coarse+1] = SolnBlk_Original_NE.WoE[j];
             SolnBlk_Coarse.oldT_E[j_coarse+1] = SolnBlk_Original_NE.oldT_E[j];
          } /* endif */
      } /* endfor */

      for ( i = SolnBlk_Original_NE.ICl;
            i <= SolnBlk_Original_NE.ICu+SolnBlk_Original_NE.Nghost ; i += 2) {
      	  i_coarse = (i-SolnBlk_Original_NE.ICl)/2+
                     (SolnBlk_Coarse.ICu-SolnBlk_Coarse.ICl+1)/2+SolnBlk_Coarse.ICl;
          SolnBlk_Coarse.WoN[i_coarse] = SolnBlk_Original_NE.WoN[i];
          SolnBlk_Coarse.oldT_N[i_coarse] = SolnBlk_Original_NE.oldT_N[i];
          if (i == SolnBlk_Original_NE.ICu+SolnBlk_Original_NE.Nghost-1) {
             SolnBlk_Coarse.WoN[i_coarse+1] = SolnBlk_Original_NE.WoN[i];
             SolnBlk_Coarse.oldT_N[i_coarse+1] = SolnBlk_Original_NE.oldT_N[i];
          } /* endif */
      } /* endfor */

    } /* endif */
    return 0;
}

/********************************************************
 * Routine: ICs                                         *
 *                                                      *
 * Assigns initial conditions and data to the           *
 * solution variables of the specified quadrilateral    *
 * solution block.                                      *
 *                                                      *
 ********************************************************/
void ICs(Gaussian2D_Quad_Block &SolnBlk,
	 Gaussian2D_Input_Parameters &Input_Parameters,
         Gaussian2D_pState *Wo) {

    int i, j, k;
    Gaussian2D_pState Wl, Wr;
    double eta, f, fp, fpp;

    /* Assign the initial data for the IVP of interest. */

    switch(Input_Parameters.i_ICs) {
      case IC_CONSTANT :
      case IC_UNIFORM :
        // Set the solution state to the initial state Wo[0].
        for (j  = SolnBlk.JCl-SolnBlk.Nghost ; j <= SolnBlk.JCu+SolnBlk.Nghost ; ++j ) {
            for ( i = SolnBlk.ICl-SolnBlk.Nghost ; i <= SolnBlk.ICu+SolnBlk.Nghost ; ++i ) {
	       SolnBlk.W[i][j] = Wo[0];

	       /////////////////////////////////////////////////////////////////////////////////
	       ////////////////////Vortex Experiment Addition////////////////////////////////////
	       ///////////////////////////////////////////////////////////////////////////////////
	       ////////////////////////////////////////////////////////////////////////////////////
	       /////////////////////////////////////////////////////////////////////////////////////

//	       if(SolnBlk.Grid.Cell[i][j].Xc.y<0
//		  && sqrt(sqr(SolnBlk.Grid.Cell[i][j].Xc.x)+sqr(SolnBlk.Grid.Cell[i][j].Xc.y)) 
//  		     < 3.0*Input_Parameters.Cylinder_Radius) {
//		 SolnBlk.W[i][j].v.x = 0.0;
//		 SolnBlk.W[i][j].v.x = 0.0;
//	       }


               SolnBlk.U[i][j] = U(SolnBlk.W[i][j]);
            } /* endfor */
        } /* endfor */
        break;
      case IC_SOD_XDIR :
        // Set initial data for Sod IVP in x-direction.
        Wl = Gaussian2D_W_STDATM;
        Wr = Gaussian2D_pState(DENSITY_STDATM/EIGHT,
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
        Wl = Gaussian2D_W_STDATM;
        Wr = Gaussian2D_pState(DENSITY_STDATM/EIGHT,
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
        Wl = Gaussian2D_pState(4.696, ZERO, ZERO, 404.4e03);
        Wr = Gaussian2D_pState(1.408, ZERO, ZERO, 101.1e03);
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
        Wl = Gaussian2D_pState(4.696, ZERO, ZERO, 404.4e03);
        Wr = Gaussian2D_pState(1.408, ZERO, ZERO, 101.1e03);
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
        Wl = Gaussian2D_pState(ONE, -TWO, ZERO, FOUR/TEN);
        Wr = Gaussian2D_pState(ONE, TWO, ZERO, FOUR/TEN);
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
        Wl = Gaussian2D_pState(ONE, ZERO, -TWO, FOUR/TEN);
        Wr = Gaussian2D_pState(ONE, ZERO, TWO, FOUR/TEN);
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
        Wl = Gaussian2D_pState(2.281, 164.83, ZERO, 201.17e03);
        Wr = Gaussian2D_pState(1.408, ZERO, ZERO, 101.1e03);
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
      case IC_SHOCK_WAVE_YDIR :
        // Set initial data for moving shock wave propagating in y-direction.
        Wl = Gaussian2D_pState(2.281, ZERO, 164.83, 201.17e03);
        Wr = Gaussian2D_pState(1.408, ZERO, ZERO, 101.1e03);
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
      case IC_SHOCK_STRUCTURE_M1_1:
	if(Wo[0].gas != GAS_A) {
	  cout << "error....don't change gas from \"A\".\n";
	  assert(1==2);
	}
	Wl = Gaussian2D_pState(1.661, 350.7444241833, 0.0, 101325.0);
	Wr = Gaussian2D_pState(1.90955819477, 305.08951614, 0.0, 127922.8125);
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
      case IC_SHOCK_STRUCTURE_M1_2:
	if(Wo[0].gas != GAS_A) {
	  cout << "error....don't change gas from \"A\".\n";
	  assert(1==2);
	}
	Wl = Gaussian2D_pState(1.661, 382.63, 0.0, 101325.0);
	Wr = Gaussian2D_pState(2.15481, 294.944, 0.0, 157054);
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
      case IC_SHOCK_STRUCTURE_M1_3:
	if(Wo[0].gas != GAS_A) {
	  cout << "error....don't change gas from \"A\".\n";
	  assert(1==2);
	}
	Wl = Gaussian2D_pState(1.661, 414.51592216, 0.0, 101325.0);
	Wr = Gaussian2D_pState(2.39410660981, 287.585750733, 0.0, 188717.8125);
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
      case IC_SHOCK_STRUCTURE_M1_5:
	if(Wo[0].gas != GAS_A) {
	  cout << "error....don't change gas from \"A\".\n";
	  assert(1==2);
	}
	Wl = Gaussian2D_pState(1.661, 478.287602499, 0.0, 101325.0);
	Wr = Gaussian2D_pState(2.84742857143, 279.001101458, 0.0, 259645.3125);
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
      case IC_SHOCK_STRUCTURE_M2_0:
	if(Wo[0].gas != GAS_A) {
	  cout << "error....don't change gas from \"A\".\n";
	  assert(1==2);
	}
	Wl = Gaussian2D_pState(1.661, 637.716803332, 0.0, 101325.0);
	Wr = Gaussian2D_pState(3.79657142857, 279.001101458, 0.0, 481293.75);
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
      case IC_SHOCK_STRUCTURE_M9_0:
	if(Wo[0].gas != GAS_A) {
	  cout << "error....don't change gas from \"A\".\n";
	  assert(1==2);
	}
	Wl = Gaussian2D_pState(1.661, 2869.73, 0.0, 101325.0);
	Wr = Gaussian2D_pState(6.40671, 744.005, 0.0, 1.02338e7);
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
      case IC_SHOCK_STRUCTURE_M10_0:
	if(Wo[0].gas != GAS_A) {
	  cout << "error....don't change gas from \"A\".\n";
	  assert(1==2);
	}
	Wl = Gaussian2D_pState(1.661, 3188.58401666, 0.0, 101325.0);
	Wr = Gaussian2D_pState(6.45048543689, 821.06038429, 0.0, 12640293.75);
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
      case IC_CONTACT_SURFACE_XDIR :
        // Set initial data for moving contact surface propagating in x-direction.
        Wl = Gaussian2D_pState(1.045, 200.00, ZERO, 300.00e03);
        Wr = Gaussian2D_pState(3.483, 200.00, ZERO, 300.00e03);
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
        Wl = Gaussian2D_pState(1.045, ZERO, 200.00, 300.00e03);
        Wr = Gaussian2D_pState(3.483, ZERO, 200.00, 300.00e03);
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
        Wl = Gaussian2D_pState(1.598, -383.64, ZERO, 91.88e03);
        Wr = Gaussian2D_pState(2.787, -216.97, ZERO, 200.0e03);
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
        Wl = Gaussian2D_pState(1.598, ZERO, -383.64, 91.88e03);
        Wr = Gaussian2D_pState(2.787, ZERO, -216.97, 200.0e03);
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
        Wl = Gaussian2D_W_STDATM;
        Wr = Gaussian2D_pState(DENSITY_STDATM*FOUR,
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
        Wr = Gaussian2D_W_STDATM;
        Wl = Gaussian2D_pState(HUNDRED*DENSITY_STDATM,
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
        //Wl = Gaussian2D_W_STDATM;
        //Wr = Gaussian2D_pState(DENSITY_STDATM/1520.0,
        //  		    ZERO, ZERO,
        // 		    PRESSURE_STDATM/1520.0);
	Wl = Wo[0];
	Wl.v.x = 0.0;
	Wl.v.y = 0.0;
	Wr = Wl/1520.0;
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
//      case IC_ROCKET_COLD_FLOW_STEADY_STATE :
//        for (j  = SolnBlk.JCl-SolnBlk.Nghost ; j <= SolnBlk.JCu+SolnBlk.Nghost ; ++j ) {
//            for ( i = SolnBlk.ICl-SolnBlk.Nghost ; i <= SolnBlk.ICu+SolnBlk.Nghost ; ++i ) {
//	       if (SolnBlk.Grid.centroid(i,j).x < -TWO) {
//		  SolnBlk.W[i][j].d = HUNDRED*DENSITY_STDATM;
//		  SolnBlk.W[i][j].v = Vector2D_ZERO;
//		  SolnBlk.W[i][j].p.xx = HUNDRED*PRESSURE_STDATM;
//		  SolnBlk.W[i][j].p.xy = ZERO;
//		  SolnBlk.W[i][j].p.yy = HUNDRED*PRESSURE_STDATM;
//		  SolnBlk.W[i][j].p.zz = HUNDRED*PRESSURE_STDATM;
//		  SolnBlk.W[i][j].erot = HUNDRED*PRESSURE_STDATM;
//
//	       } else if (SolnBlk.Grid.centroid(i,j).x >= -TWO &&
//			  SolnBlk.Grid.centroid(i,j).x <= ZERO) {
//		  SolnBlk.W[i][j].d = FIFTY*DENSITY_STDATM;
//		  SolnBlk.W[i][j].v = Vector2D_ZERO;
//		  SolnBlk.W[i][j].p.xx = FIFTY*PRESSURE_STDATM;
//		  SolnBlk.W[i][j].p.xy = ZERO;
//		  SolnBlk.W[i][j].p.yy = FIFTY*PRESSURE_STDATM;
//		  SolnBlk.W[i][j].p.zz = FIFTY*PRESSURE_STDATM;
//		  SolnBlk.W[i][j].erot = FIFTY*PRESSURE_STDATM;
//	       } else {
//		  SolnBlk.W[i][j].d = DENSITY_STDATM;
//		  SolnBlk.W[i][j].v = Vector2D_ZERO;
//		  SolnBlk.W[i][j].p.xx = PRESSURE_STDATM;
//		  SolnBlk.W[i][j].p.xy = ZERO;
//		  SolnBlk.W[i][j].p.yy = PRESSURE_STDATM;
//		  SolnBlk.W[i][j].p.zz = PRESSURE_STDATM;
//		  SolnBlk.W[i][j].erot = PRESSURE_STDATM;
//	       } /* end if */
//	       SolnBlk.U[i][j] = U(SolnBlk.W[i][j]);
//	    } /* endfor */
//	} /* endfor */
//	break;
      case IC_WEDGE_FLOW :
	Wl.d = FOUR*DENSITY_STDATM;      Wr.d = DENSITY_STDATM;
	Wl.v = Vector2D_ZERO;            Wr.v = Vector2D_ZERO;
	Wl.p.xx = FOUR*PRESSURE_STDATM;  Wr.p.xx = PRESSURE_STDATM;
	Wl.p.xy = ZERO;                  Wr.p.xy = ZERO;
	Wl.p.yy = FOUR*PRESSURE_STDATM;  Wr.p.yy = PRESSURE_STDATM;
	Wl.p.zz = FOUR*PRESSURE_STDATM;  Wr.p.zz = PRESSURE_STDATM;
	Wl.erot = FOUR*PRESSURE_STDATM;  Wr.erot = PRESSURE_STDATM;
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
	Wl.d = FOUR*DENSITY_STDATM;      Wr.d = DENSITY_STDATM;
	Wl.v = Vector2D_ZERO;            Wr.v = Vector2D_ZERO;
	Wl.p.xx = FOUR*PRESSURE_STDATM;  Wr.p.xx = PRESSURE_STDATM;
	Wl.p.xy = ZERO;                  Wr.p.xy = ZERO;
	Wl.p.yy = FOUR*PRESSURE_STDATM;  Wr.p.yy = PRESSURE_STDATM;
	Wl.p.zz = FOUR*PRESSURE_STDATM;  Wr.p.zz = PRESSURE_STDATM;
	Wl.erot = FOUR*PRESSURE_STDATM;  Wr.erot = PRESSURE_STDATM;
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
      case IC_PIPE :
        // (2-state initial data, left to right).
        //Wr = Wo[0];
	//Wr.p.xx = Wr.p.xx - Input_Parameters.Pressure_Drop*THOUSAND;
	//Wr.p.yy = Wr.p.yy - Input_Parameters.Pressure_Drop*THOUSAND;
	//Wr.p.zz = Wr.p.zz - Input_Parameters.Pressure_Drop*THOUSAND;
	//Wr.erot = Wr.erot - Input_Parameters.Pressure_Drop*THOUSAND;
        Wl = Wo[0];
        for (j  = SolnBlk.JCl-SolnBlk.Nghost ; j <= SolnBlk.JCu+SolnBlk.Nghost ; ++j ) {
            for ( i = SolnBlk.ICl-SolnBlk.Nghost ; i <= SolnBlk.ICu+SolnBlk.Nghost ; ++i ) {
	      //if (SolnBlk.Grid.Cell[i][j].Xc.x <= Input_Parameters.Pipe_Length/4.0*3.0) {
                  SolnBlk.W[i][j] = Wl;
	       //} else {
  		  //SolnBlk.W[i][j] = Wr;
	       //} /* end if */
               SolnBlk.U[i][j] = U(SolnBlk.W[i][j]);
            } /* endfor */
        } /* endfor */
        break;
      case IC_COUETTE :
        // (2-state initial data, left to right).
        Wr = Wo[0];
	Wr.p.xx = Wr.p.xx - Input_Parameters.Pressure_Drop*THOUSAND;
	Wr.p.yy = Wr.p.yy - Input_Parameters.Pressure_Drop*THOUSAND;
	Wr.p.zz = Wr.p.zz - Input_Parameters.Pressure_Drop*THOUSAND;
	Wr.erot = Wr.erot - Input_Parameters.Pressure_Drop*THOUSAND;
        Wl = Wo[0];
        // Set the solution state to the initial state and adjust velocity & shear
        for (j  = SolnBlk.JCl-SolnBlk.Nghost ; j <= SolnBlk.JCu+SolnBlk.Nghost ; ++j ) {
            for ( i = SolnBlk.ICl-SolnBlk.Nghost ; i <= SolnBlk.ICu+SolnBlk.Nghost ; ++i ) {
               if (SolnBlk.Grid.Cell[i][j].Xc.x <= ZERO) {
                  SolnBlk.W[i][j] = Wl;
               } else {
                  SolnBlk.W[i][j] = Wr;	     
               } /* end if */
	       SolnBlk.W[i][j].v.x = SolnBlk.Grid.Cell[i][j].Xc.y/
                                               (Input_Parameters.Couette_Plate_Separation/2.0)*
		                                Input_Parameters.Couette_Plate_Velocity;
	       SolnBlk.W[i][j].p.xy = -Input_Parameters.Couette_Plate_Velocity*2.0/
		                      Input_Parameters.Couette_Plate_Separation*SolnBlk.W[i][j].viscosity();
               SolnBlk.U[i][j] = U(SolnBlk.W[i][j]);
            } /* endfor */
        } /* endfor */
        break;
      case IC_VISCOUS_FLAT_PLATE:
	for (j = SolnBlk.JCl-SolnBlk.Nghost; j <= SolnBlk.JCu+SolnBlk.Nghost; j++) {
	  for (i = SolnBlk.ICl-SolnBlk.Nghost; i <= SolnBlk.ICu+SolnBlk.Nghost; i++) {
	    if (SolnBlk.Grid.Cell[i][j].Xc.y >= ZERO) {
	      SolnBlk.W[i][j] = FlatPlate(Wo[0],SolnBlk.Grid.Cell[i][j].Xc,eta,f,fp,fpp);
	    } else {
	      SolnBlk.W[i][j] = FlatPlate(Wo[0],Vector2D(SolnBlk.Grid.Cell[i][j].Xc.x,-SolnBlk.Grid.Cell[i][j].Xc.y),eta,f,fp,fpp);
	    }
	    SolnBlk.U[i][j] = U(SolnBlk.W[i][j]);
	  }
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
	  SolnBlk.oldT_W[j] = SolnBlk.WoW[j].T();
	  SolnBlk.oldT_E[j] = SolnBlk.WoE[j].T();
       } else if (j < SolnBlk.JCl) {
          SolnBlk.WoW[j] = SolnBlk.W[SolnBlk.ICl][SolnBlk.JCl];
          SolnBlk.WoE[j] = SolnBlk.W[SolnBlk.ICu][SolnBlk.JCl];
	  SolnBlk.oldT_W[j] = SolnBlk.WoW[j].T();
	  SolnBlk.oldT_E[j] = SolnBlk.WoE[j].T();
       } else {
          SolnBlk.WoW[j] = SolnBlk.W[SolnBlk.ICl][SolnBlk.JCu];
          SolnBlk.WoE[j] = SolnBlk.W[SolnBlk.ICu][SolnBlk.JCu];
	  SolnBlk.oldT_W[j] = SolnBlk.WoW[j].T();
	  SolnBlk.oldT_E[j] = SolnBlk.WoE[j].T();
       } /* endif */

       //Set Wall velocities for Adiabatic walls or viscous isothermal
       if (SolnBlk.Grid.BCtypeW[j] == BC_ADIABATIC_WALL ||
	   SolnBlk.Grid.BCtypeW[j] == BC_WALL_VISCOUS_ISOTHERMAL || 
	   SolnBlk.Grid.BCtypeW[j] == BC_TEMPERATURE_SLIP){
	 SolnBlk.WoW[j].set_temperature_d(Input_Parameters.Temperature_West_BC);
	 SolnBlk.WoW[j].v.y = 0.0;
	 SolnBlk.WoW[j].v.x = 0.0;
	 SolnBlk.oldT_W[j] = SolnBlk.WoW[j].T();
       } /* endif */
       if (SolnBlk.Grid.BCtypeE[j] == BC_ADIABATIC_WALL ||
	   SolnBlk.Grid.BCtypeE[j] == BC_WALL_VISCOUS_ISOTHERMAL ||
	   SolnBlk.Grid.BCtypeE[j] == BC_TEMPERATURE_SLIP){
	 SolnBlk.WoE[j].set_temperature_d(Input_Parameters.Temperature_East_BC);
	 SolnBlk.WoE[j].v.y = 0.0;
	 SolnBlk.WoE[j].v.x = 0.0;
	 SolnBlk.oldT_E[j] = SolnBlk.WoE[j].T();
       } /* endif */

       //Set Pressure drop for Pipe
       if(Input_Parameters.i_ICs == IC_PIPE) {
	 SolnBlk.WoE[j] = Gaussian2D_pState(SolnBlk.WoE[j].d,SolnBlk.WoE[j].v.x,SolnBlk.WoE[j].v.y,
					    SolnBlk.WoE[j].pressure()-Input_Parameters.Pressure_Drop*THOUSAND);
       }

    } /* endfor */

    for ( i = SolnBlk.ICl-SolnBlk.Nghost ; i <= SolnBlk.ICu+SolnBlk.Nghost ; ++i ) {
       if (i >= SolnBlk.ICl && i <= SolnBlk.ICu) {
          SolnBlk.WoS[i] = SolnBlk.W[i][SolnBlk.JCl];
          SolnBlk.WoN[i] = SolnBlk.W[i][SolnBlk.JCu];
          SolnBlk.oldT_S[i] = SolnBlk.WoS[i].T();
          SolnBlk.oldT_N[i] = SolnBlk.WoN[i].T();
       } else if (i < SolnBlk.ICl) {
          SolnBlk.WoS[i] = SolnBlk.W[SolnBlk.ICl][SolnBlk.JCl];
          SolnBlk.WoN[i] = SolnBlk.W[SolnBlk.ICl][SolnBlk.JCu];
          SolnBlk.oldT_S[i] = SolnBlk.WoS[i].T();
          SolnBlk.oldT_N[i] = SolnBlk.WoN[i].T();
       } else {
          SolnBlk.WoS[i] = SolnBlk.W[SolnBlk.ICu][SolnBlk.JCl];
          SolnBlk.WoN[i] = SolnBlk.W[SolnBlk.ICu][SolnBlk.JCu];
          SolnBlk.oldT_S[i] = SolnBlk.WoS[i].T();
          SolnBlk.oldT_N[i] = SolnBlk.WoN[i].T();
       } /* endif */

       //Set Wall velocities for Adiabatic walls
       if (SolnBlk.Grid.BCtypeS[i] == BC_ADIABATIC_WALL ||
	   SolnBlk.Grid.BCtypeS[i] == BC_WALL_VISCOUS_ISOTHERMAL ||
	   SolnBlk.Grid.BCtypeS[i] == BC_TEMPERATURE_SLIP){
	 SolnBlk.WoS[i].set_temperature_d(Input_Parameters.Temperature_South_BC);
	 SolnBlk.WoS[i].v.y = 0.0;
	 if (Input_Parameters.i_Grid == GRID_ADIABATIC_COUETTE) {
	   SolnBlk.WoS[i].v.x =  - Input_Parameters.Couette_Plate_Velocity;
	 } else {
	   SolnBlk.WoS[i].v.x = 0.0;
	 } /* endif */
	 SolnBlk.oldT_S[i] = SolnBlk.WoS[i].T();
       } /* endif */
       if (SolnBlk.Grid.BCtypeN[i] == BC_ADIABATIC_WALL ||
	   SolnBlk.Grid.BCtypeN[i] == BC_WALL_VISCOUS_ISOTHERMAL ||
	   SolnBlk.Grid.BCtypeN[i] == BC_TEMPERATURE_SLIP){
	 SolnBlk.WoN[i].set_temperature_d(Input_Parameters.Temperature_North_BC);
	 SolnBlk.WoN[i].v.y = 0.0;
	 if (Input_Parameters.i_Grid == GRID_ADIABATIC_COUETTE) {
	 SolnBlk.WoN[i].v.x = Input_Parameters.Couette_Plate_Velocity;
	 } else {
	   SolnBlk.WoN[i].v.x = 0.0;
	 } /* endif */
	 SolnBlk.oldT_N[i] = SolnBlk.WoN[i].T();
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
void BCs(Gaussian2D_Quad_Block &SolnBlk,
	 Gaussian2D_Input_Parameters &IP) {

    int i, j;
    Vector2D dX;
    Gaussian2D_pState dW, W;

    for ( j = SolnBlk.JCl-SolnBlk.Nghost ; j <= SolnBlk.JCu+SolnBlk.Nghost ; ++j ) {
      if ( (j >= SolnBlk.JCl && j <= SolnBlk.JCu) ||
           (j < SolnBlk.JCl && 
            (SolnBlk.Grid.BCtypeS[SolnBlk.ICl-1] == BC_NONE ||
             SolnBlk.Grid.BCtypeS[SolnBlk.ICl-1] == BC_PERIODIC ||
             SolnBlk.Grid.BCtypeS[SolnBlk.ICl-1] == BC_CONSTANT_EXTRAPOLATION ||
             SolnBlk.Grid.BCtypeS[SolnBlk.ICl-1] == BC_LINEAR_EXTRAPOLATION ||
             SolnBlk.Grid.BCtypeS[SolnBlk.ICl-1] == BC_COUETTE ||
             SolnBlk.Grid.BCtypeS[SolnBlk.ICl-1] == BC_CHARACTERISTIC ||
             SolnBlk.Grid.BCtypeS[SolnBlk.ICl-1] == BC_CHARACTERISTIC_VELOCITY) ) ||
           (j > SolnBlk.JCu && 
            (SolnBlk.Grid.BCtypeN[SolnBlk.ICl-1] == BC_NONE ||
             SolnBlk.Grid.BCtypeN[SolnBlk.ICl-1] == BC_PERIODIC ||
             SolnBlk.Grid.BCtypeN[SolnBlk.ICl-1] == BC_CONSTANT_EXTRAPOLATION ||
             SolnBlk.Grid.BCtypeN[SolnBlk.ICl-1] == BC_LINEAR_EXTRAPOLATION ||
             SolnBlk.Grid.BCtypeN[SolnBlk.ICl-1] == BC_COUETTE ||
             SolnBlk.Grid.BCtypeN[SolnBlk.ICl-1] == BC_CHARACTERISTIC ||
             SolnBlk.Grid.BCtypeN[SolnBlk.ICl-1] == BC_CHARACTERISTIC_VELOCITY) ) ) {
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
            Linear_Reconstruction_LeastSquares(SolnBlk, 
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
            SolnBlk.W[SolnBlk.ICl-2][j] = Reflect(SolnBlk.W[SolnBlk.ICl+1][j],
						  SolnBlk.Grid.nfaceW(SolnBlk.ICl,j));
            SolnBlk.U[SolnBlk.ICl-2][j] = U(SolnBlk.W[SolnBlk.ICl-2][j]);
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
	  case BC_CHARACTERISTIC_VELOCITY :
            SolnBlk.W[SolnBlk.ICl-1][j] = 
               BC_Characteristic_Velocity(SolnBlk.W[SolnBlk.ICl][j],
                                          SolnBlk.WoW[j], 
                                          SolnBlk.Grid.nfaceW(SolnBlk.ICl,j));
            SolnBlk.U[SolnBlk.ICl-1][j] = U(SolnBlk.W[SolnBlk.ICl-1][j]);
            SolnBlk.W[SolnBlk.ICl-2][j] = SolnBlk.W[SolnBlk.ICl-1][j];
            SolnBlk.U[SolnBlk.ICl-2][j] = U(SolnBlk.W[SolnBlk.ICl-2][j]);
            break;
	  case BC_ADIABATIC_WALL :
	    SolnBlk.W[SolnBlk.ICl-1][j] = Adiabatic_Wall(SolnBlk.W[SolnBlk.ICl][j],
							SolnBlk.WoW[j].v,
							SolnBlk.Grid.nfaceW(SolnBlk.ICl,j));
            SolnBlk.U[SolnBlk.ICl-1][j] = U(SolnBlk.W[SolnBlk.ICl-1][j]);
            SolnBlk.W[SolnBlk.ICl-2][j] = SolnBlk.W[SolnBlk.ICl-1][j];
            SolnBlk.U[SolnBlk.ICl-2][j] = U(SolnBlk.W[SolnBlk.ICl-2][j]);
	    break;
	  case BC_WALL_VISCOUS_ISOTHERMAL :
	    SolnBlk.W[SolnBlk.ICl-1][j] = Isothermal_Wall(SolnBlk.W[SolnBlk.ICl][j],
							  SolnBlk.WoW[j].v, SolnBlk.WoW[j].T(),
							  SolnBlk.Grid.nfaceW(SolnBlk.ICl,j));
            SolnBlk.U[SolnBlk.ICl-1][j] = U(SolnBlk.W[SolnBlk.ICl-1][j]);
            SolnBlk.W[SolnBlk.ICl-2][j] = SolnBlk.W[SolnBlk.ICl-1][j];
            SolnBlk.U[SolnBlk.ICl-2][j] = U(SolnBlk.W[SolnBlk.ICl-2][j]);
	    break;
	  case BC_TEMPERATURE_SLIP :
	    SolnBlk.W[SolnBlk.ICl-1][j] = Isothermal_Wall_Slip_T(SolnBlk.W[SolnBlk.ICl][j],
								 SolnBlk.WoW[j].v, SolnBlk.WoW[j].T(),
								 SolnBlk.oldT_W[j],
								 SolnBlk.phi[SolnBlk.ICl][j]^SolnBlk.dWdx[SolnBlk.ICl][j],
								 SolnBlk.phi[SolnBlk.ICl][j]^SolnBlk.dWdy[SolnBlk.ICl][j],
								 SolnBlk.Grid.nfaceW(SolnBlk.ICl,j));
            SolnBlk.U[SolnBlk.ICl-1][j] = U(SolnBlk.W[SolnBlk.ICl-1][j]);
            SolnBlk.W[SolnBlk.ICl-2][j] = SolnBlk.W[SolnBlk.ICl-1][j];
            SolnBlk.U[SolnBlk.ICl-2][j] = U(SolnBlk.W[SolnBlk.ICl-2][j]);
	    break;
 	  case BC_COUETTE :
            Linear_Reconstruction_LeastSquares(SolnBlk, 
                                                 SolnBlk.ICl, j, 
                                                 LIMITER_BARTH_JESPERSEN);
            dX = SolnBlk.Grid.Cell[SolnBlk.ICl-1][j].Xc -
                 SolnBlk.Grid.Cell[SolnBlk.ICl][j].Xc;
            SolnBlk.W[SolnBlk.ICl-1][j] = BC_Couette(SolnBlk.W[SolnBlk.ICl][j] + 
               (SolnBlk.phi[SolnBlk.ICl][j]^SolnBlk.dWdx[SolnBlk.ICl][j])*dX.x +
               (SolnBlk.phi[SolnBlk.ICl][j]^SolnBlk.dWdy[SolnBlk.ICl][j])*dX.y,SolnBlk.WoW[j]);
            SolnBlk.U[SolnBlk.ICl-1][j] = U(SolnBlk.W[SolnBlk.ICl-1][j]);
            dX = SolnBlk.Grid.Cell[SolnBlk.ICl-2][j].Xc -
                 SolnBlk.Grid.Cell[SolnBlk.ICl][j].Xc;
            SolnBlk.W[SolnBlk.ICl-2][j] = BC_Couette(SolnBlk.W[SolnBlk.ICl][j] + 
               (SolnBlk.phi[SolnBlk.ICl][j]^SolnBlk.dWdx[SolnBlk.ICl][j])*dX.x +
               (SolnBlk.phi[SolnBlk.ICl][j]^SolnBlk.dWdy[SolnBlk.ICl][j])*dX.y,SolnBlk.WoW[j]);
            SolnBlk.U[SolnBlk.ICl-2][j] = U(SolnBlk.W[SolnBlk.ICl-2][j]);
            break;
	  case BC_DEVELOPED_CHANNEL :
	    SolnBlk.W[SolnBlk.ICl-1][j] = BC_Developed_Channel_Flow(SolnBlk.WoW[j],
								    SolnBlk.Grid.nfaceW(SolnBlk.ICl,j),
								    SolnBlk.Grid.Cell[SolnBlk.ICl-1][j].Xc.y);
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
        } /* endswitch */
      } /* endif */

      if ( (j >= SolnBlk.JCl && j <= SolnBlk.JCu) ||
           (j < SolnBlk.JCl && 
            (SolnBlk.Grid.BCtypeS[SolnBlk.ICu+1] == BC_NONE ||
             SolnBlk.Grid.BCtypeS[SolnBlk.ICu+1] == BC_PERIODIC ||
             SolnBlk.Grid.BCtypeS[SolnBlk.ICu+1] == BC_CONSTANT_EXTRAPOLATION ||
             SolnBlk.Grid.BCtypeS[SolnBlk.ICu+1] == BC_LINEAR_EXTRAPOLATION ||
             SolnBlk.Grid.BCtypeS[SolnBlk.ICu+1] == BC_COUETTE ||
             SolnBlk.Grid.BCtypeS[SolnBlk.ICu+1] == BC_CHARACTERISTIC ||
             SolnBlk.Grid.BCtypeS[SolnBlk.ICu+1] == BC_CHARACTERISTIC_VELOCITY) ) ||
           (j > SolnBlk.JCu && 
            (SolnBlk.Grid.BCtypeN[SolnBlk.ICu+1] == BC_NONE ||
             SolnBlk.Grid.BCtypeN[SolnBlk.ICu+1] == BC_PERIODIC ||
             SolnBlk.Grid.BCtypeN[SolnBlk.ICu+1] == BC_CONSTANT_EXTRAPOLATION ||
             SolnBlk.Grid.BCtypeN[SolnBlk.ICu+1] == BC_LINEAR_EXTRAPOLATION ||
             SolnBlk.Grid.BCtypeN[SolnBlk.ICu+1] == BC_COUETTE ||
             SolnBlk.Grid.BCtypeN[SolnBlk.ICu+1] == BC_CHARACTERISTIC ||
             SolnBlk.Grid.BCtypeN[SolnBlk.ICu+1] == BC_CHARACTERISTIC_VELOCITY) ) ) {
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
          case BC_LINEAR_EXTRAPOLATION :
            Linear_Reconstruction_LeastSquares(SolnBlk, 
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
	  case BC_CHARACTERISTIC_VELOCITY :
            SolnBlk.W[SolnBlk.ICu+1][j] = 
               BC_Characteristic_Velocity(SolnBlk.W[SolnBlk.ICu][j],
                                          SolnBlk.WoE[j],
                                          SolnBlk.Grid.nfaceE(SolnBlk.ICu,j));
            SolnBlk.U[SolnBlk.ICu+1][j] = U(SolnBlk.W[SolnBlk.ICu+1][j]);
            SolnBlk.W[SolnBlk.ICu+2][j] = SolnBlk.W[SolnBlk.ICu+1][j];
            SolnBlk.U[SolnBlk.ICu+2][j] = U(SolnBlk.W[SolnBlk.ICu+2][j]);
            break;
	  case BC_ADIABATIC_WALL :
	    SolnBlk.W[SolnBlk.ICu+1][j] = Adiabatic_Wall(SolnBlk.W[SolnBlk.ICu][j],
							SolnBlk.WoE[j].v,
							SolnBlk.Grid.nfaceE(SolnBlk.ICu,j));
            SolnBlk.U[SolnBlk.ICu+1][j] = U(SolnBlk.W[SolnBlk.ICu+1][j]);
            SolnBlk.W[SolnBlk.ICu+2][j] = SolnBlk.W[SolnBlk.ICu+1][j];
            SolnBlk.U[SolnBlk.ICu+2][j] = U(SolnBlk.W[SolnBlk.ICu+2][j]);
	    break;
	  case BC_WALL_VISCOUS_ISOTHERMAL :
	    SolnBlk.W[SolnBlk.ICu+1][j] = Isothermal_Wall(SolnBlk.W[SolnBlk.ICu][j],
							  SolnBlk.WoE[j].v,SolnBlk.WoE[j].T(),
							  SolnBlk.Grid.nfaceE(SolnBlk.ICu,j));
            SolnBlk.U[SolnBlk.ICu+1][j] = U(SolnBlk.W[SolnBlk.ICu+1][j]);
            SolnBlk.W[SolnBlk.ICu+2][j] = SolnBlk.W[SolnBlk.ICu+1][j];
            SolnBlk.U[SolnBlk.ICu+2][j] = U(SolnBlk.W[SolnBlk.ICu+2][j]);
	    break;
	  case BC_TEMPERATURE_SLIP :
	    SolnBlk.W[SolnBlk.ICu+1][j] = Isothermal_Wall_Slip_T(SolnBlk.W[SolnBlk.ICu][j],
								 SolnBlk.WoE[j].v,SolnBlk.WoE[j].T(),
								 SolnBlk.oldT_E[j],
								 SolnBlk.phi[SolnBlk.ICu][j]^SolnBlk.dWdx[SolnBlk.ICu][j],
								 SolnBlk.phi[SolnBlk.ICu][j]^SolnBlk.dWdy[SolnBlk.ICu][j],
								 SolnBlk.Grid.nfaceE(SolnBlk.ICu,j));
            SolnBlk.U[SolnBlk.ICu+1][j] = U(SolnBlk.W[SolnBlk.ICu+1][j]);
            SolnBlk.W[SolnBlk.ICu+2][j] = SolnBlk.W[SolnBlk.ICu+1][j];
            SolnBlk.U[SolnBlk.ICu+2][j] = U(SolnBlk.W[SolnBlk.ICu+2][j]);
	    break;
 	  case BC_COUETTE :
            Linear_Reconstruction_LeastSquares(SolnBlk, 
                                                 SolnBlk.ICu, j, 
                                                 LIMITER_BARTH_JESPERSEN);
            dX = SolnBlk.Grid.Cell[SolnBlk.ICu+1][j].Xc -
                 SolnBlk.Grid.Cell[SolnBlk.ICu][j].Xc;
            SolnBlk.W[SolnBlk.ICu+1][j] = BC_Couette(SolnBlk.W[SolnBlk.ICu][j] + 
               (SolnBlk.phi[SolnBlk.ICu][j]^SolnBlk.dWdx[SolnBlk.ICu][j])*dX.x +
               (SolnBlk.phi[SolnBlk.ICu][j]^SolnBlk.dWdy[SolnBlk.ICu][j])*dX.y,SolnBlk.WoE[j]);
            SolnBlk.U[SolnBlk.ICu+1][j] = U(SolnBlk.W[SolnBlk.ICu+1][j]);
            dX = SolnBlk.Grid.Cell[SolnBlk.ICu+2][j].Xc -
                 SolnBlk.Grid.Cell[SolnBlk.ICu][j].Xc;
            SolnBlk.W[SolnBlk.ICu+2][j] = BC_Couette(SolnBlk.W[SolnBlk.ICu][j] + 
               (SolnBlk.phi[SolnBlk.ICu][j]^SolnBlk.dWdx[SolnBlk.ICu][j])*dX.x +
               (SolnBlk.phi[SolnBlk.ICu][j]^SolnBlk.dWdy[SolnBlk.ICu][j])*dX.y,SolnBlk.WoE[j]);
            SolnBlk.U[SolnBlk.ICu+2][j] = U(SolnBlk.W[SolnBlk.ICu+2][j]);
            break;
	  case BC_DEVELOPED_CHANNEL :
	    SolnBlk.W[SolnBlk.ICu+1][j] = BC_Developed_Channel_Flow(SolnBlk.WoE[j],
								    SolnBlk.Grid.nfaceE(SolnBlk.ICu,j),
								    SolnBlk.Grid.Cell[SolnBlk.ICu+1][j].Xc.y);
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
        } /* endswitch */
      } /* endif */
    } /* endfor */

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
            Linear_Reconstruction_LeastSquares(SolnBlk, 
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
	  case BC_CHARACTERISTIC_VELOCITY :
            SolnBlk.W[i][SolnBlk.JCl-1] = 
               BC_Characteristic_Velocity(SolnBlk.W[i][SolnBlk.JCl],
                                          SolnBlk.WoS[i],
                                          SolnBlk.Grid.nfaceS(i, SolnBlk.JCl));
            SolnBlk.U[i][SolnBlk.JCl-1] = U(SolnBlk.W[i][SolnBlk.JCl-1]);
            SolnBlk.W[i][SolnBlk.JCl-2] = SolnBlk.W[i][SolnBlk.JCl-1];
            SolnBlk.U[i][SolnBlk.JCl-2] = U(SolnBlk.W[i][SolnBlk.JCl-2]);
            break;
	  case BC_ADIABATIC_WALL :
	    SolnBlk.W[i][SolnBlk.JCl-1] = Adiabatic_Wall(SolnBlk.W[i][SolnBlk.JCl],
							SolnBlk.WoS[i].v,
							SolnBlk.Grid.nfaceS(i,SolnBlk.JCl));
            SolnBlk.U[i][SolnBlk.JCl-1] = U(SolnBlk.W[i][SolnBlk.JCl-1]);
            SolnBlk.W[i][SolnBlk.JCl-2] = SolnBlk.W[i][SolnBlk.JCl-1];
            SolnBlk.U[i][SolnBlk.JCl-2] = U(SolnBlk.W[i][SolnBlk.JCl-2]);
	    break;
	  case BC_WALL_VISCOUS_ISOTHERMAL :
	    SolnBlk.W[i][SolnBlk.JCl-1] = Isothermal_Wall(SolnBlk.W[i][SolnBlk.JCl],
							  SolnBlk.WoS[i].v, SolnBlk.WoS[i].T(),
							  SolnBlk.Grid.nfaceS(i,SolnBlk.JCl));
            SolnBlk.U[i][SolnBlk.JCl-1] = U(SolnBlk.W[i][SolnBlk.JCl-1]);
            SolnBlk.W[i][SolnBlk.JCl-2] = SolnBlk.W[i][SolnBlk.JCl-1];
            SolnBlk.U[i][SolnBlk.JCl-2] = U(SolnBlk.W[i][SolnBlk.JCl-2]);
	    break;
	  case BC_TEMPERATURE_SLIP :
	    SolnBlk.W[i][SolnBlk.JCl-1] = Isothermal_Wall_Slip_T(SolnBlk.W[i][SolnBlk.JCl],
								 SolnBlk.WoS[i].v, SolnBlk.WoS[i].T(),
								 SolnBlk.oldT_S[i],
								 SolnBlk.phi[i][SolnBlk.JCl]^SolnBlk.dWdx[i][SolnBlk.JCl],
								 SolnBlk.phi[i][SolnBlk.JCl]^SolnBlk.dWdy[i][SolnBlk.JCl],
								 SolnBlk.Grid.nfaceS(i,SolnBlk.JCl));
            SolnBlk.U[i][SolnBlk.JCl-1] = U(SolnBlk.W[i][SolnBlk.JCl-1]);
            SolnBlk.W[i][SolnBlk.JCl-2] = SolnBlk.W[i][SolnBlk.JCl-1];
            SolnBlk.U[i][SolnBlk.JCl-2] = U(SolnBlk.W[i][SolnBlk.JCl-2]);
	    break;
	  case BC_COUETTE :
            Linear_Reconstruction_LeastSquares(SolnBlk, 
                                                 i, SolnBlk.JCl, 
                                                 LIMITER_BARTH_JESPERSEN);
            dX = SolnBlk.Grid.Cell[i][SolnBlk.JCl-1].Xc -
                 SolnBlk.Grid.Cell[i][SolnBlk.JCl].Xc;
            SolnBlk.W[i][SolnBlk.JCl-1] = BC_Couette(SolnBlk.W[i][SolnBlk.JCl] + 
               (SolnBlk.phi[i][SolnBlk.JCl]^SolnBlk.dWdx[i][SolnBlk.JCl])*dX.x +
               (SolnBlk.phi[i][SolnBlk.JCl]^SolnBlk.dWdy[i][SolnBlk.JCl])*dX.y,	SolnBlk.WoS[i]);
            SolnBlk.U[i][SolnBlk.JCl-1] = U(SolnBlk.W[i][SolnBlk.JCl-1]);
            dX = SolnBlk.Grid.Cell[i][SolnBlk.JCl-2].Xc -
                 SolnBlk.Grid.Cell[i][SolnBlk.JCl].Xc;
            SolnBlk.W[i][SolnBlk.JCl-2] = BC_Couette(SolnBlk.W[i][SolnBlk.JCl] + 
               (SolnBlk.phi[i][SolnBlk.JCl]^SolnBlk.dWdx[i][SolnBlk.JCl])*dX.x +
               (SolnBlk.phi[i][SolnBlk.JCl]^SolnBlk.dWdy[i][SolnBlk.JCl])*dX.y,	SolnBlk.WoS[i]);
            SolnBlk.U[i][SolnBlk.JCl-2] = U(SolnBlk.W[i][SolnBlk.JCl-2]);
            break;
          default:
            SolnBlk.W[i][SolnBlk.JCl-1] = SolnBlk.W[i][SolnBlk.JCl];
            SolnBlk.U[i][SolnBlk.JCl-1] = SolnBlk.U[i][SolnBlk.JCl];
            SolnBlk.W[i][SolnBlk.JCl-2] = SolnBlk.W[i][SolnBlk.JCl];
            SolnBlk.U[i][SolnBlk.JCl-2] = SolnBlk.U[i][SolnBlk.JCl];
            break;
        } /* endswitch */
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
            Linear_Reconstruction_LeastSquares(SolnBlk, 
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
	  case BC_CHARACTERISTIC_VELOCITY :
            SolnBlk.W[i][SolnBlk.JCu+1] = 
               BC_Characteristic_Velocity(SolnBlk.W[i][SolnBlk.JCu],
                                          SolnBlk.WoN[i],
                                          SolnBlk.Grid.nfaceN(i, SolnBlk.JCu));
            SolnBlk.U[i][SolnBlk.JCu+1] = U(SolnBlk.W[i][SolnBlk.JCu+1]);
            SolnBlk.W[i][SolnBlk.JCu+2] = SolnBlk.W[i][SolnBlk.JCu+1];
            SolnBlk.U[i][SolnBlk.JCu+2] = U(SolnBlk.W[i][SolnBlk.JCu+2]);
            break;
	  case BC_ADIABATIC_WALL :
	    SolnBlk.W[i][SolnBlk.JCu+1] = Adiabatic_Wall(SolnBlk.W[i][SolnBlk.JCu],
							SolnBlk.WoN[i].v,
							SolnBlk.Grid.nfaceN(i,SolnBlk.JCu));
            SolnBlk.U[i][SolnBlk.JCu+1] = U(SolnBlk.W[i][SolnBlk.JCu+1]);
            SolnBlk.W[i][SolnBlk.JCu+2] = SolnBlk.W[i][SolnBlk.JCu+1];
            SolnBlk.U[i][SolnBlk.JCu+2] = U(SolnBlk.W[i][SolnBlk.JCu+2]);
	    break;
	  case BC_WALL_VISCOUS_ISOTHERMAL :
	    SolnBlk.W[i][SolnBlk.JCu+1] = Isothermal_Wall(SolnBlk.W[i][SolnBlk.JCu],
							  SolnBlk.WoN[i].v, SolnBlk.WoN[i].T(),
							  SolnBlk.Grid.nfaceN(i,SolnBlk.JCu));
            SolnBlk.U[i][SolnBlk.JCu+1] = U(SolnBlk.W[i][SolnBlk.JCu+1]);
            SolnBlk.W[i][SolnBlk.JCu+2] = SolnBlk.W[i][SolnBlk.JCu+1];
            SolnBlk.U[i][SolnBlk.JCu+2] = U(SolnBlk.W[i][SolnBlk.JCu+2]);
	    break;
	  case BC_TEMPERATURE_SLIP :
	    SolnBlk.W[i][SolnBlk.JCu+1] = Isothermal_Wall_Slip_T(SolnBlk.W[i][SolnBlk.JCu],
								 SolnBlk.WoN[i].v, SolnBlk.WoN[i].T(),
								 SolnBlk.oldT_N[i],
								 SolnBlk.phi[i][SolnBlk.JCu]^SolnBlk.dWdx[i][SolnBlk.JCu],
								 SolnBlk.phi[i][SolnBlk.JCu]^SolnBlk.dWdy[i][SolnBlk.JCu],
								 SolnBlk.Grid.nfaceN(i,SolnBlk.JCu));
            SolnBlk.U[i][SolnBlk.JCu+1] = U(SolnBlk.W[i][SolnBlk.JCu+1]);
            SolnBlk.W[i][SolnBlk.JCu+2] = SolnBlk.W[i][SolnBlk.JCu+1];
            SolnBlk.U[i][SolnBlk.JCu+2] = U(SolnBlk.W[i][SolnBlk.JCu+2]);
	    break;
	  case BC_COUETTE :
            Linear_Reconstruction_LeastSquares(SolnBlk, 
                                                 i, SolnBlk.JCu, 
                                                 LIMITER_BARTH_JESPERSEN);
            dX = SolnBlk.Grid.Cell[i][SolnBlk.JCu+1].Xc -
                 SolnBlk.Grid.Cell[i][SolnBlk.JCu].Xc;
            SolnBlk.W[i][SolnBlk.JCu+1] = BC_Couette(SolnBlk.W[i][SolnBlk.JCu] + 
               (SolnBlk.phi[i][SolnBlk.JCu]^SolnBlk.dWdx[i][SolnBlk.JCu])*dX.x +
               (SolnBlk.phi[i][SolnBlk.JCu]^SolnBlk.dWdy[i][SolnBlk.JCu])*dX.y,	SolnBlk.WoN[i]);
            SolnBlk.U[i][SolnBlk.JCu+1] = U(SolnBlk.W[i][SolnBlk.JCu+1]);
            dX = SolnBlk.Grid.Cell[i][SolnBlk.JCu+2].Xc -
                 SolnBlk.Grid.Cell[i][SolnBlk.JCu].Xc;
            SolnBlk.W[i][SolnBlk.JCu+2] = BC_Couette(SolnBlk.W[i][SolnBlk.JCu] + 
               (SolnBlk.phi[i][SolnBlk.JCu]^SolnBlk.dWdx[i][SolnBlk.JCu])*dX.x +
               (SolnBlk.phi[i][SolnBlk.JCu]^SolnBlk.dWdy[i][SolnBlk.JCu])*dX.y,	SolnBlk.WoN[i]);
            SolnBlk.U[i][SolnBlk.JCu+2] = U(SolnBlk.W[i][SolnBlk.JCu+2]);
            break;
          default:
            SolnBlk.W[i][SolnBlk.JCu+1] = SolnBlk.W[i][SolnBlk.JCu];
            SolnBlk.U[i][SolnBlk.JCu+1] = SolnBlk.U[i][SolnBlk.JCu];
            SolnBlk.W[i][SolnBlk.JCu+2] = SolnBlk.W[i][SolnBlk.JCu];
            SolnBlk.U[i][SolnBlk.JCu+2] = SolnBlk.U[i][SolnBlk.JCu];
            break;
        } /* endswitch */
    } /* endfor */

    //if east and west BCs are periodic, update corner ghost cells.
    //this has to be done after the North and South BCs have been
    //calculated.
    for(j = SolnBlk.JCl-SolnBlk.Nghost ; j < SolnBlk.JCl ; ++j ) {
      if(SolnBlk.Grid.BCtypeW[j] == BC_PERIODIC) {
	SolnBlk.W[SolnBlk.ICl-1][j] = SolnBlk.W[SolnBlk.ICu-1][j];
	SolnBlk.U[SolnBlk.ICl-1][j] = SolnBlk.U[SolnBlk.ICu-1][j];
	SolnBlk.W[SolnBlk.ICl-2][j] = SolnBlk.W[SolnBlk.ICu-2][j];
	SolnBlk.U[SolnBlk.ICl-2][j] = SolnBlk.U[SolnBlk.ICu-2][j];
      }
      if(SolnBlk.Grid.BCtypeE[j] == BC_PERIODIC) {
            SolnBlk.W[SolnBlk.ICu+1][j] = SolnBlk.W[SolnBlk.ICl+1][j];
            SolnBlk.U[SolnBlk.ICu+1][j] = SolnBlk.U[SolnBlk.ICl+1][j];
            SolnBlk.W[SolnBlk.ICu+2][j] = SolnBlk.W[SolnBlk.ICl+2][j];
            SolnBlk.U[SolnBlk.ICu+2][j] = SolnBlk.U[SolnBlk.ICl+2][j];
      }
    }

    for(j = SolnBlk.JCu+1; j <= SolnBlk.JCu+SolnBlk.Nghost ; ++j ) {
      if(SolnBlk.Grid.BCtypeW[j] == BC_PERIODIC) {
	SolnBlk.W[SolnBlk.ICl-1][j] = SolnBlk.W[SolnBlk.ICu-1][j];
	SolnBlk.U[SolnBlk.ICl-1][j] = SolnBlk.U[SolnBlk.ICu-1][j];
	SolnBlk.W[SolnBlk.ICl-2][j] = SolnBlk.W[SolnBlk.ICu-2][j];
	SolnBlk.U[SolnBlk.ICl-2][j] = SolnBlk.U[SolnBlk.ICu-2][j];
      }
      if(SolnBlk.Grid.BCtypeE[j] == BC_PERIODIC) {
            SolnBlk.W[SolnBlk.ICu+1][j] = SolnBlk.W[SolnBlk.ICl+1][j];
            SolnBlk.U[SolnBlk.ICu+1][j] = SolnBlk.U[SolnBlk.ICl+1][j];
            SolnBlk.W[SolnBlk.ICu+2][j] = SolnBlk.W[SolnBlk.ICl+2][j];
            SolnBlk.U[SolnBlk.ICu+2][j] = SolnBlk.U[SolnBlk.ICl+2][j];
      }
    }
    return;
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
double CFL(Gaussian2D_Quad_Block &SolnBlk) {

    int i, j;
    double cos_angleN, cos_angleS, cos_angleE, cos_angleW;
    double sin_angleN, sin_angleS, sin_angleE, sin_angleW;
    double cos_anglei, cos_anglej, sin_anglei, sin_anglej;
    double dtMin, d_i, d_j, v_i, v_j, a_i, a_j, p_ii, p_jj;
    double dt_heat;

    dtMin = MILLION;

    for ( j  = SolnBlk.JCl-SolnBlk.Nghost ; j <= SolnBlk.JCu+SolnBlk.Nghost ; ++j ) {
       for ( i = SolnBlk.ICl-SolnBlk.Nghost ; i <= SolnBlk.ICu+SolnBlk.Nghost ; ++i ) {
 	  if (i < SolnBlk.ICl || i > SolnBlk.ICu ||
              j < SolnBlk.JCl || j > SolnBlk.JCu) {
	     SolnBlk.dt[i][j] = ZERO;
          } else {

 	     cos_angleN = SolnBlk.Grid.nfaceN(i,j).x;
 	     cos_angleS = SolnBlk.Grid.nfaceS(i,j).x;
 	     cos_angleE = SolnBlk.Grid.nfaceE(i,j).x;
 	     cos_angleW = SolnBlk.Grid.nfaceW(i,j).x;

 	     sin_angleN = SolnBlk.Grid.nfaceN(i,j).y;
 	     sin_angleS = SolnBlk.Grid.nfaceS(i,j).y;
 	     sin_angleE = SolnBlk.Grid.nfaceE(i,j).y;
 	     sin_angleW = SolnBlk.Grid.nfaceW(i,j).y;

	     cos_anglei = HALF*(cos_angleE-cos_angleW);
	     cos_anglej = HALF*(cos_angleN-cos_angleS);
	     sin_anglei = HALF*(sin_angleE-sin_angleW);
	     sin_anglej = HALF*(sin_angleN-sin_angleS);

   	     d_i = TWO*(SolnBlk.Grid.Cell[i][j].A/
                   (SolnBlk.Grid.lfaceE(i, j)+SolnBlk.Grid.lfaceW(i, j)));
             d_j = TWO*(SolnBlk.Grid.Cell[i][j].A/
                   (SolnBlk.Grid.lfaceN(i, j)+SolnBlk.Grid.lfaceS(i, j)));

	     v_i = HALF*(SolnBlk.W[i][j].v*
                   (SolnBlk.Grid.nfaceE(i, j)-SolnBlk.Grid.nfaceW(i, j)));
	     v_j = HALF*(SolnBlk.W[i][j].v*
                   (SolnBlk.Grid.nfaceN(i, j)-SolnBlk.Grid.nfaceS(i, j)));

	     p_ii = SolnBlk.W[i][j].p.xx*cos_anglei*cos_anglei
                    +SolnBlk.W[i][j].p.yy*sin_anglei*sin_anglei
                    +2.0*SolnBlk.W[i][j].p.xy*cos_anglei*sin_anglei;
	     p_jj = SolnBlk.W[i][j].p.xx*cos_anglej*cos_anglej
                    +SolnBlk.W[i][j].p.yy*sin_anglej*sin_anglej
                    +2.0*SolnBlk.W[i][j].p.xy*cos_anglej*sin_anglej;

	     // Inviscid dt calculation.
	     a_i = sqrt(3.0*p_ii/SolnBlk.W[i][j].d);
	     a_j = sqrt(3.0*p_jj/SolnBlk.W[i][j].d);
	     SolnBlk.dt[i][j] = min(d_i/(a_i+fabs(v_i)), d_j/(a_j+fabs(v_j)));

#ifdef _GAUSSIAN_HEAT_TRANSFER_
	     // Heat-Transfer-related dt calculation.
	     if(SolnBlk.Heat_Transfer) {
	       dt_heat = HALF*min(sqr(d_i),sqr(d_j))/(SolnBlk.W[i][j].pr*SolnBlk.W[i][j].tt()*SolnBlk.W[i][j].pressure()/SolnBlk.W[i][j].d);
	       SolnBlk.dt[i][j] = min(SolnBlk.dt[i][j],dt_heat);
	     }
#endif

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
void Set_Global_TimeStep(Gaussian2D_Quad_Block &SolnBlk,
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
double L1_Norm_Residual(Gaussian2D_Quad_Block &SolnBlk) {

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

/********************************************************
 * Routine: L2_Norm_Residual                            *
 *                                                      *
 * Determines the L2-norm of the solution residual for  *
 * the specified quadrilateral solution block.          *
 * Useful for monitoring convergence of the solution    *
 * for steady state problems.                           *
 *                                                      *
 ********************************************************/
double L2_Norm_Residual(Gaussian2D_Quad_Block &SolnBlk) {

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

/********************************************************
 * Routine: Max_Norm_Residual                           *
 *                                                      *
 * Determines the maximum norm of the solution residual *
 * for the specified quadrilateral solution block.      *
 * Useful for monitoring convergence of the solution    *
 * for steady state problems.                           *
 *                                                      *
 ********************************************************/
double Max_Norm_Residual(Gaussian2D_Quad_Block &SolnBlk) {

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
void Linear_Reconstruction_GreenGauss(Gaussian2D_Quad_Block &SolnBlk,
				      const int i, 
                                      const int j,
                                      const int Limiter) {

    int n, n2, n_pts, ii, i_index[8], j_index[8];
    double u0Min, u0Max, uQuad[4], phi;
    double DxDx_ave, DxDy_ave, DyDy_ave;
    double l_north, l_south, l_east, l_west;
    double a, b, c, d;
    Vector2D n_north, n_south, n_east, n_west, dX;
    Gaussian2D_pState W_nw, W_ne, W_sw, W_se, W_face, 
                      DU, DUDx_ave, DUDy_ave, W_temp;

    /* Carry out the limited solution reconstruction in
       each cell of the computational mesh. */

      //If on outside ring of ghost cells, no reconstruction

    if (i == SolnBlk.ICl-SolnBlk.Nghost || i == SolnBlk.ICu+SolnBlk.Nghost ||
        j == SolnBlk.JCl-SolnBlk.Nghost || j == SolnBlk.JCu+SolnBlk.Nghost) {
      n_pts = 0;

      //If on inside ring of ghost cells, reconstruction may be required

      /////////////left edge//////////////////

    } else if (i == SolnBlk.ICl-SolnBlk.Nghost+1) {
        //bottom or top left corner (never used to evaluate fluxes, therefore
        //no reconstruction needed
      if(j == SolnBlk.JCl-SolnBlk.Nghost+1 || j == SolnBlk.JCu+SolnBlk.Nghost-1){
	n_pts = 0;
	//left edge, may need reconstruction
        //(only for BC_NONE and BC_PERIODIC)???
      } else if(SolnBlk.Grid.BCtypeW[j] == BC_FIXED ||
                SolnBlk.Grid.BCtypeW[j] == BC_CHARACTERISTIC ||
                SolnBlk.Grid.BCtypeW[j] == BC_LINEAR_EXTRAPOLATION ||
                SolnBlk.Grid.BCtypeW[j] == BC_COUETTE ||
                SolnBlk.Grid.BCtypeW[j] == BC_CONSTANT_EXTRAPOLATION ||
                SolnBlk.Grid.BCtypeW[j] == BC_REFLECTION ||
                SolnBlk.Grid.BCtypeW[j] == BC_WALL_VISCOUS_ISOTHERMAL ||
                SolnBlk.Grid.BCtypeW[j] == BC_TEMPERATURE_SLIP ||
                SolnBlk.Grid.BCtypeW[j] == BC_ADIABATIC_WALL) {
	n_pts = 0;
      } else { //BC's are Periodic or "none"
	  //Bottom left
	if(j == SolnBlk.JCl-SolnBlk.Nghost+2 &&
               (SolnBlk.Grid.BCtypeS[i] == BC_FIXED ||
                SolnBlk.Grid.BCtypeS[i] == BC_CHARACTERISTIC ||
                SolnBlk.Grid.BCtypeS[i] == BC_LINEAR_EXTRAPOLATION ||
                SolnBlk.Grid.BCtypeS[i] == BC_COUETTE              ||
                SolnBlk.Grid.BCtypeS[i] == BC_CONSTANT_EXTRAPOLATION ||
                SolnBlk.Grid.BCtypeS[i] == BC_REFLECTION ||
                SolnBlk.Grid.BCtypeS[i] == BC_WALL_VISCOUS_ISOTHERMAL ||
                SolnBlk.Grid.BCtypeS[i] == BC_TEMPERATURE_SLIP ||
                SolnBlk.Grid.BCtypeS[i] == BC_ADIABATIC_WALL)) {
	  n_pts = 5;
	  i_index[0] = i-1; j_index[0] = j;
	  i_index[1] = i-1; j_index[1] = j+1;
	  i_index[2] = i;   j_index[2] = j+1;
	  i_index[3] = i+1; j_index[3] = j+1;
	  i_index[4] = i+1; j_index[4] = j;
	  //Top Left
	} else if (j == SolnBlk.JCu+SolnBlk.Nghost-2 &&
                       (SolnBlk.Grid.BCtypeN[i] == BC_FIXED ||
                        SolnBlk.Grid.BCtypeN[i] == BC_CHARACTERISTIC ||
                        SolnBlk.Grid.BCtypeN[i] == BC_LINEAR_EXTRAPOLATION ||
                        SolnBlk.Grid.BCtypeN[i] == BC_COUETTE ||
                        SolnBlk.Grid.BCtypeN[i] == BC_CONSTANT_EXTRAPOLATION ||
                        SolnBlk.Grid.BCtypeN[i] == BC_REFLECTION ||
                        SolnBlk.Grid.BCtypeN[i] == BC_WALL_VISCOUS_ISOTHERMAL ||
                        SolnBlk.Grid.BCtypeN[i] == BC_TEMPERATURE_SLIP ||
                        SolnBlk.Grid.BCtypeN[i] == BC_ADIABATIC_WALL)) {
	  n_pts = 5;
	  i_index[0] = i-1; j_index[0] = j;
	  i_index[1] = i-1; j_index[1] = j-1;
	  i_index[2] = i;   j_index[2] = j-1;
	  i_index[3] = i+1; j_index[3] = j-1;
	  i_index[4] = i+1; j_index[4] = j;
	  //middle left
	} else {
	  n_pts = 8;
	  i_index[0] = i-1; j_index[0] = j-1;
	  i_index[1] = i-1; j_index[1] = j;
	  i_index[2] = i-1; j_index[2] = j+1;
	  i_index[3] = i;   j_index[3] = j+1;
	  i_index[4] = i+1; j_index[4] = j+1;
	  i_index[5] = i+1; j_index[5] = j;
	  i_index[6] = i+1; j_index[6] = j-1;
	  i_index[7] = i;   j_index[7] = j-1;
	} //endif
      } //endif

      /////////////Right edge//////////////////

    } else if (i == SolnBlk.ICu+SolnBlk.Nghost-1) {
        //bottom or top left corner (never used to evaluate fluxes, therefore
        //no reconstruction needed
      if(j == SolnBlk.JCl-SolnBlk.Nghost+1 || j == SolnBlk.JCu+SolnBlk.Nghost-1){
	n_pts = 0;
	//right edge, may need reconstruction
        //(only for BC_NONE and BC_PERIODIC)???
      } else if(SolnBlk.Grid.BCtypeE[j] == BC_FIXED ||
                SolnBlk.Grid.BCtypeE[j] == BC_CHARACTERISTIC ||
                SolnBlk.Grid.BCtypeE[j] == BC_LINEAR_EXTRAPOLATION ||
                SolnBlk.Grid.BCtypeE[j] == BC_COUETTE ||
                SolnBlk.Grid.BCtypeE[j] == BC_CONSTANT_EXTRAPOLATION ||
                SolnBlk.Grid.BCtypeE[j] == BC_REFLECTION ||
                SolnBlk.Grid.BCtypeE[j] == BC_WALL_VISCOUS_ISOTHERMAL ||
                SolnBlk.Grid.BCtypeE[j] == BC_TEMPERATURE_SLIP ||
                SolnBlk.Grid.BCtypeE[j] == BC_ADIABATIC_WALL) {
	n_pts = 0;
      } else { //BC's are Periodic or "none"
	  //Bottom right
	if(j == SolnBlk.JCl-SolnBlk.Nghost+2 &&
               (SolnBlk.Grid.BCtypeS[i] == BC_FIXED ||
                SolnBlk.Grid.BCtypeS[i] == BC_CHARACTERISTIC ||
                SolnBlk.Grid.BCtypeS[i] == BC_LINEAR_EXTRAPOLATION ||
                SolnBlk.Grid.BCtypeS[i] == BC_COUETTE ||
                SolnBlk.Grid.BCtypeS[i] == BC_CONSTANT_EXTRAPOLATION ||
                SolnBlk.Grid.BCtypeS[i] == BC_REFLECTION ||
                SolnBlk.Grid.BCtypeS[i] == BC_WALL_VISCOUS_ISOTHERMAL ||
                SolnBlk.Grid.BCtypeS[i] == BC_TEMPERATURE_SLIP ||
                SolnBlk.Grid.BCtypeS[i] == BC_ADIABATIC_WALL)) {
	  n_pts = 5;
	  i_index[0] = i-1; j_index[0] = j;
	  i_index[1] = i-1; j_index[1] = j+1;
	  i_index[2] = i;   j_index[2] = j+1;
	  i_index[3] = i+1; j_index[3] = j+1;
	  i_index[4] = i+1; j_index[4] = j;
	  //Top right
	} else if (j == SolnBlk.JCu+SolnBlk.Nghost-2 &&
                       (SolnBlk.Grid.BCtypeN[i] == BC_FIXED ||
                        SolnBlk.Grid.BCtypeN[i] == BC_CHARACTERISTIC ||
                        SolnBlk.Grid.BCtypeN[i] == BC_LINEAR_EXTRAPOLATION ||
                        SolnBlk.Grid.BCtypeN[i] == BC_COUETTE              ||
                        SolnBlk.Grid.BCtypeN[i] == BC_CONSTANT_EXTRAPOLATION ||
                        SolnBlk.Grid.BCtypeN[i] == BC_REFLECTION ||
                        SolnBlk.Grid.BCtypeN[i] == BC_WALL_VISCOUS_ISOTHERMAL ||
                        SolnBlk.Grid.BCtypeN[i] == BC_TEMPERATURE_SLIP ||
                        SolnBlk.Grid.BCtypeN[i] == BC_ADIABATIC_WALL)) {
	  n_pts = 5;
	  i_index[0] = i-1; j_index[0] = j;
	  i_index[1] = i-1; j_index[1] = j-1;
	  i_index[2] = i;   j_index[2] = j-1;
	  i_index[3] = i+1; j_index[3] = j-1;
	  i_index[4] = i+1; j_index[4] = j;
	  //middle right or top or bottom with none or priodic bc's
	} else {
	  n_pts = 8;
	  i_index[0] = i-1; j_index[0] = j-1;
	  i_index[1] = i-1; j_index[1] = j;
	  i_index[2] = i-1; j_index[2] = j+1;
	  i_index[3] = i;   j_index[3] = j+1;
	  i_index[4] = i+1; j_index[4] = j+1;
	  i_index[5] = i+1; j_index[5] = j;
	  i_index[6] = i+1; j_index[6] = j-1;
	  i_index[7] = i;   j_index[7] = j-1;
	} //endif
      } //endif

      /////////////bottom edge/////////////
      // (corners already taken care of) //

    } else if (j == SolnBlk.JCl-SolnBlk.Nghost+1) {
      if(SolnBlk.Grid.BCtypeS[i] == BC_FIXED ||
         SolnBlk.Grid.BCtypeS[i] == BC_CHARACTERISTIC ||
         SolnBlk.Grid.BCtypeS[i] == BC_LINEAR_EXTRAPOLATION ||
         SolnBlk.Grid.BCtypeS[i] == BC_COUETTE              ||
         SolnBlk.Grid.BCtypeS[i] == BC_CONSTANT_EXTRAPOLATION ||
         SolnBlk.Grid.BCtypeS[i] == BC_REFLECTION ||
         SolnBlk.Grid.BCtypeS[i] == BC_WALL_VISCOUS_ISOTHERMAL ||
         SolnBlk.Grid.BCtypeS[i] == BC_TEMPERATURE_SLIP ||
         SolnBlk.Grid.BCtypeS[i] == BC_ADIABATIC_WALL) {
	n_pts = 0;
      } else { //BC's are Periodic or "none"
	  //Bottom left
	if(i == SolnBlk.ICl-SolnBlk.Nghost+2 &&
               (SolnBlk.Grid.BCtypeW[j] == BC_FIXED ||
                SolnBlk.Grid.BCtypeW[j] == BC_CHARACTERISTIC ||
                SolnBlk.Grid.BCtypeW[j] == BC_LINEAR_EXTRAPOLATION ||
                SolnBlk.Grid.BCtypeW[j] == BC_COUETTE              ||
                SolnBlk.Grid.BCtypeW[j] == BC_CONSTANT_EXTRAPOLATION ||
                SolnBlk.Grid.BCtypeW[j] == BC_REFLECTION ||
                SolnBlk.Grid.BCtypeW[j] == BC_WALL_VISCOUS_ISOTHERMAL ||
                SolnBlk.Grid.BCtypeW[j] == BC_TEMPERATURE_SLIP ||
                SolnBlk.Grid.BCtypeW[j] == BC_ADIABATIC_WALL)) {
	  n_pts = 5;
	  i_index[0] = i;   j_index[0] = j+1;
	  i_index[1] = i+1; j_index[1] = j+1;
	  i_index[2] = i+1; j_index[2] = j;
	  i_index[3] = i+1; j_index[3] = j-1;
	  i_index[4] = i;   j_index[4] = j-1;
	  //Bottom right
	} else if (i == SolnBlk.ICu+SolnBlk.Nghost-2 &&
                       (SolnBlk.Grid.BCtypeE[j] == BC_FIXED ||
                        SolnBlk.Grid.BCtypeE[j] == BC_CHARACTERISTIC ||
                        SolnBlk.Grid.BCtypeE[j] == BC_LINEAR_EXTRAPOLATION ||
                        SolnBlk.Grid.BCtypeE[j] == BC_COUETTE ||
                        SolnBlk.Grid.BCtypeE[j] == BC_CONSTANT_EXTRAPOLATION ||
                        SolnBlk.Grid.BCtypeE[j] == BC_REFLECTION ||
                        SolnBlk.Grid.BCtypeE[j] == BC_WALL_VISCOUS_ISOTHERMAL ||
                        SolnBlk.Grid.BCtypeE[j] == BC_TEMPERATURE_SLIP ||
                        SolnBlk.Grid.BCtypeE[j] == BC_ADIABATIC_WALL)) {
	  n_pts = 5;
	  i_index[0] = i;   j_index[0] = j+1;
	  i_index[1] = i-1; j_index[1] = j+1;
	  i_index[2] = i-1; j_index[2] = j;
	  i_index[3] = i-1; j_index[3] = j-1;
	  i_index[4] = i;   j_index[4] = j-1;
	  //middle bottom or left or right with none or priodic bc's
	} else {
	  n_pts = 8;
	  i_index[0] = i-1; j_index[0] = j-1;
	  i_index[1] = i-1; j_index[1] = j;
	  i_index[2] = i-1; j_index[2] = j+1;
	  i_index[3] = i;   j_index[3] = j+1;
	  i_index[4] = i+1; j_index[4] = j+1;
	  i_index[5] = i+1; j_index[5] = j;
	  i_index[6] = i+1; j_index[6] = j-1;
	  i_index[7] = i;   j_index[7] = j-1;
	} //endif
      } //endif

      /////////////top edge/////////////
      // (corners already taken care of) //

    } else if (j == SolnBlk.JCu+SolnBlk.Nghost-1) {
      if(SolnBlk.Grid.BCtypeN[i] == BC_FIXED ||
         SolnBlk.Grid.BCtypeN[i] == BC_CHARACTERISTIC ||
         SolnBlk.Grid.BCtypeN[i] == BC_LINEAR_EXTRAPOLATION ||
         SolnBlk.Grid.BCtypeN[i] == BC_COUETTE ||
         SolnBlk.Grid.BCtypeN[i] == BC_CONSTANT_EXTRAPOLATION ||
         SolnBlk.Grid.BCtypeN[i] == BC_REFLECTION ||
         SolnBlk.Grid.BCtypeN[i] == BC_WALL_VISCOUS_ISOTHERMAL ||
         SolnBlk.Grid.BCtypeN[i] == BC_TEMPERATURE_SLIP ||
         SolnBlk.Grid.BCtypeN[i] == BC_ADIABATIC_WALL) {
	n_pts = 0;
      } else { //BC's are Periodic or "none"
	  //Top left
	if(i == SolnBlk.ICl-SolnBlk.Nghost+2 &&
               (SolnBlk.Grid.BCtypeW[j] == BC_FIXED ||
                SolnBlk.Grid.BCtypeW[j] == BC_CHARACTERISTIC ||
                SolnBlk.Grid.BCtypeW[j] == BC_LINEAR_EXTRAPOLATION ||
                SolnBlk.Grid.BCtypeW[j] == BC_COUETTE ||
                SolnBlk.Grid.BCtypeW[j] == BC_CONSTANT_EXTRAPOLATION ||
                SolnBlk.Grid.BCtypeW[j] == BC_REFLECTION ||
                SolnBlk.Grid.BCtypeW[j] == BC_WALL_VISCOUS_ISOTHERMAL ||
                SolnBlk.Grid.BCtypeW[j] == BC_TEMPERATURE_SLIP ||
                SolnBlk.Grid.BCtypeW[j] == BC_ADIABATIC_WALL)) {
	  n_pts = 5;
	  i_index[0] = i;   j_index[0] = j+1;
	  i_index[1] = i+1; j_index[1] = j+1;
	  i_index[2] = i+1; j_index[2] = j;
	  i_index[3] = i+1; j_index[3] = j-1;
	  i_index[4] = i;   j_index[4] = j-1;
	  //Top right
	} else if (i == SolnBlk.ICu+SolnBlk.Nghost-2 &&
                       (SolnBlk.Grid.BCtypeE[j] == BC_FIXED ||
                        SolnBlk.Grid.BCtypeE[j] == BC_CHARACTERISTIC ||
                        SolnBlk.Grid.BCtypeE[j] == BC_LINEAR_EXTRAPOLATION ||
                        SolnBlk.Grid.BCtypeE[j] == BC_COUETTE ||
                        SolnBlk.Grid.BCtypeE[j] == BC_CONSTANT_EXTRAPOLATION ||
                        SolnBlk.Grid.BCtypeE[j] == BC_REFLECTION ||
                        SolnBlk.Grid.BCtypeE[j] == BC_WALL_VISCOUS_ISOTHERMAL ||
                        SolnBlk.Grid.BCtypeE[j] == BC_TEMPERATURE_SLIP ||
                        SolnBlk.Grid.BCtypeE[j] == BC_ADIABATIC_WALL)) {
	  n_pts = 5;
	  i_index[0] = i;   j_index[0] = j+1;
	  i_index[1] = i-1; j_index[1] = j+1;
	  i_index[2] = i-1; j_index[2] = j;
	  i_index[3] = i-1; j_index[3] = j-1;
	  i_index[4] = i;   j_index[4] = j-1;
	  //top bottom or left or right with none or priodic bc's
	} else {
	  n_pts = 8;
	  i_index[0] = i-1; j_index[0] = j-1;
	  i_index[1] = i-1; j_index[1] = j;
	  i_index[2] = i-1; j_index[2] = j+1;
	  i_index[3] = i;   j_index[3] = j+1;
	  i_index[4] = i+1; j_index[4] = j+1;
	  i_index[5] = i+1; j_index[5] = j;
	  i_index[6] = i+1; j_index[6] = j-1;
	  i_index[7] = i;   j_index[7] = j-1;
	} //endif
      } //endif
    
      /////////First Cell Inside Computational Domain////////

      //bottom left
    } else if (i == SolnBlk.ICl-SolnBlk.Nghost+2 && 
	       j == SolnBlk.JCl-SolnBlk.Nghost+2) {
      if((SolnBlk.Grid.BCtypeW[j] == BC_CHARACTERISTIC ||
          SolnBlk.Grid.BCtypeW[j] == BC_LINEAR_EXTRAPOLATION ||
          SolnBlk.Grid.BCtypeW[j] == BC_COUETTE ||
          SolnBlk.Grid.BCtypeW[j] == BC_CONSTANT_EXTRAPOLATION ||
          SolnBlk.Grid.BCtypeW[j] == BC_REFLECTION ||
          SolnBlk.Grid.BCtypeW[j] == BC_WALL_VISCOUS_ISOTHERMAL ||
          SolnBlk.Grid.BCtypeW[j] == BC_TEMPERATURE_SLIP ||
          SolnBlk.Grid.BCtypeW[j] == BC_ADIABATIC_WALL)
                           &&
         (SolnBlk.Grid.BCtypeS[i] == BC_CHARACTERISTIC ||
          SolnBlk.Grid.BCtypeS[i] == BC_LINEAR_EXTRAPOLATION ||
          SolnBlk.Grid.BCtypeS[i] == BC_COUETTE ||
          SolnBlk.Grid.BCtypeS[i] == BC_CONSTANT_EXTRAPOLATION ||
          SolnBlk.Grid.BCtypeS[i] == BC_REFLECTION ||
          SolnBlk.Grid.BCtypeS[i] == BC_WALL_VISCOUS_ISOTHERMAL ||
          SolnBlk.Grid.BCtypeS[i] == BC_TEMPERATURE_SLIP ||
          SolnBlk.Grid.BCtypeS[i] == BC_ADIABATIC_WALL)) {
	n_pts = 3;
	i_index[0] = i;   j_index[0] = j+1;
	i_index[1] = i+1; j_index[1] = j+1;
	i_index[2] = i+1; j_index[2] = j;
      } else if (SolnBlk.Grid.BCtypeS[i] == BC_CHARACTERISTIC ||
                 SolnBlk.Grid.BCtypeS[i] == BC_LINEAR_EXTRAPOLATION ||
                 SolnBlk.Grid.BCtypeS[i] == BC_COUETTE ||
                 SolnBlk.Grid.BCtypeS[i] == BC_CONSTANT_EXTRAPOLATION ||
                 SolnBlk.Grid.BCtypeS[i] == BC_REFLECTION ||
                 SolnBlk.Grid.BCtypeS[i] == BC_WALL_VISCOUS_ISOTHERMAL ||
                 SolnBlk.Grid.BCtypeS[i] == BC_TEMPERATURE_SLIP ||
                 SolnBlk.Grid.BCtypeS[i] == BC_ADIABATIC_WALL) {
	n_pts = 5;
	i_index[0] = i-1; j_index[0] = j;
	i_index[1] = i-1; j_index[1] = j+1;
	i_index[2] = i;   j_index[2] = j+1;
	i_index[3] = i+1; j_index[3] = j+1;
	i_index[4] = i+1; j_index[4] = j;
      } else if (SolnBlk.Grid.BCtypeW[j] == BC_CHARACTERISTIC ||
                 SolnBlk.Grid.BCtypeW[j] == BC_LINEAR_EXTRAPOLATION ||
                 SolnBlk.Grid.BCtypeW[j] == BC_COUETTE ||
                 SolnBlk.Grid.BCtypeW[j] == BC_CONSTANT_EXTRAPOLATION ||
                 SolnBlk.Grid.BCtypeW[j] == BC_REFLECTION ||
                 SolnBlk.Grid.BCtypeW[j] == BC_WALL_VISCOUS_ISOTHERMAL ||
                 SolnBlk.Grid.BCtypeW[j] == BC_TEMPERATURE_SLIP ||
                 SolnBlk.Grid.BCtypeW[j] == BC_ADIABATIC_WALL) {
	n_pts = 5;
	i_index[0] = i;   j_index[0] = j+1;
	i_index[1] = i+1; j_index[1] = j+1;
	i_index[2] = i+1; j_index[2] = j;
	i_index[3] = i+1; j_index[3] = j-1;
	i_index[4] = i;   j_index[4] = j-1;
      } else {
        n_pts = 8;
	i_index[0] = i-1; j_index[0] = j-1;
	i_index[1] = i-1; j_index[1] = j;
	i_index[2] = i-1; j_index[2] = j+1;
	i_index[3] = i;   j_index[3] = j+1;
	i_index[4] = i+1; j_index[4] = j+1;
	i_index[5] = i+1; j_index[5] = j;
	i_index[6] = i+1; j_index[6] = j-1;
	i_index[7] = i;   j_index[7] = j-1;
      } //endif

      //top left
    } else if (i == SolnBlk.ICl-SolnBlk.Nghost+2 && 
	       j == SolnBlk.JCu+SolnBlk.Nghost-2) {
      if((SolnBlk.Grid.BCtypeW[j] == BC_CHARACTERISTIC ||
          SolnBlk.Grid.BCtypeW[j] == BC_LINEAR_EXTRAPOLATION ||
          SolnBlk.Grid.BCtypeW[j] == BC_COUETTE ||
          SolnBlk.Grid.BCtypeW[j] == BC_CONSTANT_EXTRAPOLATION ||
          SolnBlk.Grid.BCtypeW[j] == BC_REFLECTION ||
          SolnBlk.Grid.BCtypeW[j] == BC_WALL_VISCOUS_ISOTHERMAL ||
          SolnBlk.Grid.BCtypeW[j] == BC_TEMPERATURE_SLIP ||
          SolnBlk.Grid.BCtypeW[j] == BC_ADIABATIC_WALL)
                           &&
         (SolnBlk.Grid.BCtypeN[i] == BC_CHARACTERISTIC ||
          SolnBlk.Grid.BCtypeN[i] == BC_LINEAR_EXTRAPOLATION ||
          SolnBlk.Grid.BCtypeN[i] == BC_COUETTE ||
          SolnBlk.Grid.BCtypeN[i] == BC_CONSTANT_EXTRAPOLATION ||
          SolnBlk.Grid.BCtypeN[i] == BC_REFLECTION ||
          SolnBlk.Grid.BCtypeN[i] == BC_WALL_VISCOUS_ISOTHERMAL ||
          SolnBlk.Grid.BCtypeN[i] == BC_TEMPERATURE_SLIP ||
          SolnBlk.Grid.BCtypeN[i] == BC_ADIABATIC_WALL)) {
	n_pts = 3;
	i_index[0] = i;   j_index[0] = j-1;
	i_index[1] = i+1; j_index[1] = j-1;
	i_index[2] = i+1; j_index[2] = j;
      } else if (SolnBlk.Grid.BCtypeN[i] == BC_CHARACTERISTIC ||
                 SolnBlk.Grid.BCtypeN[i] == BC_LINEAR_EXTRAPOLATION ||
                 SolnBlk.Grid.BCtypeN[i] == BC_COUETTE ||
                 SolnBlk.Grid.BCtypeN[i] == BC_CONSTANT_EXTRAPOLATION ||
                 SolnBlk.Grid.BCtypeN[i] == BC_REFLECTION ||
                 SolnBlk.Grid.BCtypeN[i] == BC_WALL_VISCOUS_ISOTHERMAL ||
                 SolnBlk.Grid.BCtypeN[i] == BC_TEMPERATURE_SLIP ||
                 SolnBlk.Grid.BCtypeN[i] == BC_ADIABATIC_WALL) {
	n_pts = 5;
	i_index[0] = i-1; j_index[0] = j;
	i_index[1] = i-1; j_index[1] = j-1;
	i_index[2] = i;   j_index[2] = j-1;
	i_index[3] = i+1; j_index[3] = j-1;
	i_index[4] = i+1; j_index[4] = j;
      } else if (SolnBlk.Grid.BCtypeW[j] == BC_CHARACTERISTIC ||
                 SolnBlk.Grid.BCtypeW[j] == BC_LINEAR_EXTRAPOLATION ||
                 SolnBlk.Grid.BCtypeW[j] == BC_COUETTE ||
                 SolnBlk.Grid.BCtypeW[j] == BC_CONSTANT_EXTRAPOLATION ||
                 SolnBlk.Grid.BCtypeW[j] == BC_REFLECTION ||
                 SolnBlk.Grid.BCtypeW[j] == BC_WALL_VISCOUS_ISOTHERMAL ||
                 SolnBlk.Grid.BCtypeW[j] == BC_TEMPERATURE_SLIP ||
                 SolnBlk.Grid.BCtypeW[j] == BC_ADIABATIC_WALL) {
	n_pts = 5;
	i_index[0] = i;   j_index[0] = j+1;
	i_index[1] = i+1; j_index[1] = j+1;
	i_index[2] = i+1; j_index[2] = j;
	i_index[3] = i+1; j_index[3] = j-1;
	i_index[4] = i;   j_index[4] = j-1;
      } else {
        n_pts = 8;
	i_index[0] = i-1; j_index[0] = j-1;
	i_index[1] = i-1; j_index[1] = j;
	i_index[2] = i-1; j_index[2] = j+1;
	i_index[3] = i;   j_index[3] = j+1;
	i_index[4] = i+1; j_index[4] = j+1;
	i_index[5] = i+1; j_index[5] = j;
	i_index[6] = i+1; j_index[6] = j-1;
	i_index[7] = i;   j_index[7] = j-1;
      } //endif

      //top right
    } else if (i == SolnBlk.ICu+SolnBlk.Nghost-2 && 
	       j == SolnBlk.JCu+SolnBlk.Nghost-2) {
      if((SolnBlk.Grid.BCtypeE[j] == BC_CHARACTERISTIC ||
          SolnBlk.Grid.BCtypeE[j] == BC_LINEAR_EXTRAPOLATION ||
          SolnBlk.Grid.BCtypeE[j] == BC_COUETTE ||
          SolnBlk.Grid.BCtypeE[j] == BC_CONSTANT_EXTRAPOLATION ||
          SolnBlk.Grid.BCtypeE[j] == BC_REFLECTION ||
          SolnBlk.Grid.BCtypeE[j] == BC_WALL_VISCOUS_ISOTHERMAL ||
          SolnBlk.Grid.BCtypeE[j] == BC_TEMPERATURE_SLIP ||
          SolnBlk.Grid.BCtypeE[j] == BC_ADIABATIC_WALL)
                           &&
         (SolnBlk.Grid.BCtypeN[i] == BC_CHARACTERISTIC ||
          SolnBlk.Grid.BCtypeN[i] == BC_LINEAR_EXTRAPOLATION ||
          SolnBlk.Grid.BCtypeN[i] == BC_COUETTE ||
          SolnBlk.Grid.BCtypeN[i] == BC_CONSTANT_EXTRAPOLATION ||
          SolnBlk.Grid.BCtypeN[i] == BC_REFLECTION ||
          SolnBlk.Grid.BCtypeN[i] == BC_WALL_VISCOUS_ISOTHERMAL ||
          SolnBlk.Grid.BCtypeN[i] == BC_TEMPERATURE_SLIP ||
          SolnBlk.Grid.BCtypeN[i] == BC_ADIABATIC_WALL)) {
	n_pts = 3;
	i_index[0] = i-1; j_index[0] = j;
	i_index[1] = i-1; j_index[1] = j-1;
	i_index[2] = i;   j_index[2] = j-1;
      } else if (SolnBlk.Grid.BCtypeN[i] == BC_CHARACTERISTIC ||
                 SolnBlk.Grid.BCtypeN[i] == BC_LINEAR_EXTRAPOLATION ||
                 SolnBlk.Grid.BCtypeN[i] == BC_COUETTE ||
                 SolnBlk.Grid.BCtypeN[i] == BC_CONSTANT_EXTRAPOLATION ||
                 SolnBlk.Grid.BCtypeN[i] == BC_REFLECTION ||
                 SolnBlk.Grid.BCtypeN[i] == BC_WALL_VISCOUS_ISOTHERMAL ||
                 SolnBlk.Grid.BCtypeN[i] == BC_TEMPERATURE_SLIP ||
                 SolnBlk.Grid.BCtypeN[i] == BC_ADIABATIC_WALL) {
	n_pts = 5;
	i_index[0] = i-1; j_index[0] = j;
	i_index[1] = i-1; j_index[1] = j-1;
	i_index[2] = i;   j_index[2] = j-1;
	i_index[3] = i+1; j_index[3] = j-1;
	i_index[4] = i+1; j_index[4] = j;
      } else if (SolnBlk.Grid.BCtypeE[j] == BC_CHARACTERISTIC ||
                 SolnBlk.Grid.BCtypeE[j] == BC_LINEAR_EXTRAPOLATION ||
                 SolnBlk.Grid.BCtypeE[j] == BC_COUETTE ||
                 SolnBlk.Grid.BCtypeE[j] == BC_CONSTANT_EXTRAPOLATION ||
                 SolnBlk.Grid.BCtypeE[j] == BC_REFLECTION ||
                 SolnBlk.Grid.BCtypeE[j] == BC_WALL_VISCOUS_ISOTHERMAL ||
                 SolnBlk.Grid.BCtypeE[j] == BC_TEMPERATURE_SLIP ||
                 SolnBlk.Grid.BCtypeE[j] == BC_ADIABATIC_WALL) {
	n_pts = 5;
	i_index[0] = i;   j_index[0] = j+1;
	i_index[1] = i-1; j_index[1] = j+1;
	i_index[2] = i-1; j_index[2] = j;
	i_index[3] = i-1; j_index[3] = j-1;
	i_index[4] = i;   j_index[4] = j-1;
      } else {
        n_pts = 8;
	i_index[0] = i-1; j_index[0] = j-1;
	i_index[1] = i-1; j_index[1] = j;
	i_index[2] = i-1; j_index[2] = j+1;
	i_index[3] = i;   j_index[3] = j+1;
	i_index[4] = i+1; j_index[4] = j+1;
	i_index[5] = i+1; j_index[5] = j;
	i_index[6] = i+1; j_index[6] = j-1;
	i_index[7] = i;   j_index[7] = j-1;
      } //endif

      //Bottom right
    } else if (i == SolnBlk.ICu+SolnBlk.Nghost-2 && 
	       j == SolnBlk.JCl-SolnBlk.Nghost+2) {
      if((SolnBlk.Grid.BCtypeE[j] == BC_CHARACTERISTIC ||
          SolnBlk.Grid.BCtypeE[j] == BC_LINEAR_EXTRAPOLATION ||
          SolnBlk.Grid.BCtypeE[j] == BC_COUETTE ||
          SolnBlk.Grid.BCtypeE[j] == BC_CONSTANT_EXTRAPOLATION ||
          SolnBlk.Grid.BCtypeE[j] == BC_REFLECTION ||
          SolnBlk.Grid.BCtypeE[j] == BC_WALL_VISCOUS_ISOTHERMAL ||
          SolnBlk.Grid.BCtypeE[j] == BC_TEMPERATURE_SLIP ||
          SolnBlk.Grid.BCtypeE[j] == BC_ADIABATIC_WALL)
                           &&
         (SolnBlk.Grid.BCtypeS[i] == BC_CHARACTERISTIC ||
          SolnBlk.Grid.BCtypeS[i] == BC_LINEAR_EXTRAPOLATION ||
          SolnBlk.Grid.BCtypeS[i] == BC_COUETTE ||
          SolnBlk.Grid.BCtypeS[i] == BC_CONSTANT_EXTRAPOLATION ||
          SolnBlk.Grid.BCtypeS[i] == BC_REFLECTION ||
          SolnBlk.Grid.BCtypeS[i] == BC_WALL_VISCOUS_ISOTHERMAL ||
          SolnBlk.Grid.BCtypeS[i] == BC_TEMPERATURE_SLIP ||
          SolnBlk.Grid.BCtypeS[i] == BC_ADIABATIC_WALL)) {
	n_pts = 3;
	i_index[0] = i;   j_index[0] = j+1;
	i_index[1] = i-1; j_index[1] = j+1;
	i_index[2] = i-1; j_index[2] = j;
      } else if (SolnBlk.Grid.BCtypeS[i] == BC_CHARACTERISTIC ||
                 SolnBlk.Grid.BCtypeS[i] == BC_LINEAR_EXTRAPOLATION ||
                 SolnBlk.Grid.BCtypeS[i] == BC_COUETTE ||
                 SolnBlk.Grid.BCtypeS[i] == BC_CONSTANT_EXTRAPOLATION ||
                 SolnBlk.Grid.BCtypeS[i] == BC_REFLECTION ||
                 SolnBlk.Grid.BCtypeS[i] == BC_WALL_VISCOUS_ISOTHERMAL ||
                 SolnBlk.Grid.BCtypeS[i] == BC_TEMPERATURE_SLIP ||
                 SolnBlk.Grid.BCtypeS[i] == BC_ADIABATIC_WALL) {
	n_pts = 5;
	i_index[0] = i-1; j_index[0] = j;
	i_index[1] = i-1; j_index[1] = j+1;
	i_index[2] = i;   j_index[2] = j+1;
	i_index[3] = i+1; j_index[3] = j+1;
	i_index[4] = i+1; j_index[4] = j;
      } else if (SolnBlk.Grid.BCtypeE[j] == BC_CHARACTERISTIC ||
                 SolnBlk.Grid.BCtypeE[j] == BC_LINEAR_EXTRAPOLATION ||
                 SolnBlk.Grid.BCtypeE[j] == BC_COUETTE ||
                 SolnBlk.Grid.BCtypeE[j] == BC_CONSTANT_EXTRAPOLATION ||
                 SolnBlk.Grid.BCtypeE[j] == BC_REFLECTION ||
                 SolnBlk.Grid.BCtypeE[j] == BC_WALL_VISCOUS_ISOTHERMAL ||
                 SolnBlk.Grid.BCtypeE[j] == BC_TEMPERATURE_SLIP ||
                 SolnBlk.Grid.BCtypeE[j] == BC_ADIABATIC_WALL) {
	n_pts = 5;
	i_index[0] = i;   j_index[0] = j+1;
	i_index[1] = i-1; j_index[1] = j+1;
	i_index[2] = i-1; j_index[2] = j;
	i_index[3] = i-1; j_index[3] = j-1;
	i_index[4] = i;   j_index[4] = j-1;
      } else {
        n_pts = 8;
	i_index[0] = i-1; j_index[0] = j-1;
	i_index[1] = i-1; j_index[1] = j;
	i_index[2] = i-1; j_index[2] = j+1;
	i_index[3] = i;   j_index[3] = j+1;
	i_index[4] = i+1; j_index[4] = j+1;
	i_index[5] = i+1; j_index[5] = j;
	i_index[6] = i+1; j_index[6] = j-1;
	i_index[7] = i;   j_index[7] = j-1;
      } //endif

      //left
    } else if (i == SolnBlk.ICl-SolnBlk.Nghost+2) {

      if(SolnBlk.Grid.BCtypeW[j] == BC_CHARACTERISTIC ||
         SolnBlk.Grid.BCtypeW[j] == BC_LINEAR_EXTRAPOLATION ||
         SolnBlk.Grid.BCtypeW[j] == BC_COUETTE ||
         SolnBlk.Grid.BCtypeW[j] == BC_CONSTANT_EXTRAPOLATION ||
         SolnBlk.Grid.BCtypeW[j] == BC_REFLECTION ||
         SolnBlk.Grid.BCtypeW[j] == BC_WALL_VISCOUS_ISOTHERMAL ||
         SolnBlk.Grid.BCtypeW[j] == BC_TEMPERATURE_SLIP ||
         SolnBlk.Grid.BCtypeW[j] == BC_ADIABATIC_WALL) {
	n_pts = 5;
	i_index[0] = i;   j_index[0] = j+1;
	i_index[1] = i+1; j_index[1] = j+1;
	i_index[2] = i+1; j_index[2] = j;
	i_index[3] = i+1; j_index[3] = j-1;
	i_index[4] = i;   j_index[4] = j-1;
      } else {
        n_pts = 8;
	i_index[0] = i-1; j_index[0] = j-1;
	i_index[1] = i-1; j_index[1] = j;
	i_index[2] = i-1; j_index[2] = j+1;
	i_index[3] = i;   j_index[3] = j+1;
	i_index[4] = i+1; j_index[4] = j+1;
	i_index[5] = i+1; j_index[5] = j;
	i_index[6] = i+1; j_index[6] = j-1;
	i_index[7] = i;   j_index[7] = j-1;
      } //endif

      //top
    } else if (j == SolnBlk.JCu+SolnBlk.Nghost-2) {

      if(SolnBlk.Grid.BCtypeN[i] == BC_CHARACTERISTIC ||
         SolnBlk.Grid.BCtypeN[i] == BC_LINEAR_EXTRAPOLATION ||
         SolnBlk.Grid.BCtypeN[i] == BC_COUETTE ||
         SolnBlk.Grid.BCtypeN[i] == BC_CONSTANT_EXTRAPOLATION ||
         SolnBlk.Grid.BCtypeN[i] == BC_REFLECTION ||
         SolnBlk.Grid.BCtypeN[i] == BC_WALL_VISCOUS_ISOTHERMAL ||
         SolnBlk.Grid.BCtypeN[i] == BC_TEMPERATURE_SLIP ||
         SolnBlk.Grid.BCtypeN[i] == BC_ADIABATIC_WALL) {
	n_pts = 5;
	i_index[0] = i-1; j_index[0] = j;
	i_index[1] = i-1; j_index[1] = j-1;
	i_index[2] = i;   j_index[2] = j-1;
	i_index[3] = i+1; j_index[3] = j-1;
	i_index[4] = i+1; j_index[4] = j;
      } else {
        n_pts = 8;
	i_index[0] = i-1; j_index[0] = j-1;
	i_index[1] = i-1; j_index[1] = j;
	i_index[2] = i-1; j_index[2] = j+1;
	i_index[3] = i;   j_index[3] = j+1;
	i_index[4] = i+1; j_index[4] = j+1;
	i_index[5] = i+1; j_index[5] = j;
	i_index[6] = i+1; j_index[6] = j-1;
	i_index[7] = i;   j_index[7] = j-1;
      } //endif

      //Right
    } else if (i == SolnBlk.ICu+SolnBlk.Nghost-2) {

      if(SolnBlk.Grid.BCtypeE[j] == BC_CHARACTERISTIC ||
         SolnBlk.Grid.BCtypeE[j] == BC_LINEAR_EXTRAPOLATION ||
         SolnBlk.Grid.BCtypeE[j] == BC_COUETTE ||
         SolnBlk.Grid.BCtypeE[j] == BC_CONSTANT_EXTRAPOLATION ||
         SolnBlk.Grid.BCtypeE[j] == BC_REFLECTION ||
         SolnBlk.Grid.BCtypeE[j] == BC_WALL_VISCOUS_ISOTHERMAL ||
         SolnBlk.Grid.BCtypeE[j] == BC_TEMPERATURE_SLIP ||
         SolnBlk.Grid.BCtypeE[j] == BC_ADIABATIC_WALL) {
	n_pts = 5;
	i_index[0] = i;   j_index[0] = j+1;
	i_index[1] = i-1; j_index[1] = j+1;
	i_index[2] = i-1; j_index[2] = j;
	i_index[3] = i-1; j_index[3] = j-1;
	i_index[4] = i;   j_index[4] = j-1;
      } else {
        n_pts = 8;
	i_index[0] = i-1; j_index[0] = j-1;
	i_index[1] = i-1; j_index[1] = j;
	i_index[2] = i-1; j_index[2] = j+1;
	i_index[3] = i;   j_index[3] = j+1;
	i_index[4] = i+1; j_index[4] = j+1;
	i_index[5] = i+1; j_index[5] = j;
	i_index[6] = i+1; j_index[6] = j-1;
	i_index[7] = i;   j_index[7] = j-1;
      } //endif

      //Bottom
    } else if (j == SolnBlk.JCl-SolnBlk.Nghost+2) {

      if(SolnBlk.Grid.BCtypeS[i] == BC_CHARACTERISTIC ||
         SolnBlk.Grid.BCtypeS[i] == BC_LINEAR_EXTRAPOLATION ||
         SolnBlk.Grid.BCtypeS[i] == BC_COUETTE ||
         SolnBlk.Grid.BCtypeS[i] == BC_CONSTANT_EXTRAPOLATION ||
         SolnBlk.Grid.BCtypeS[i] == BC_REFLECTION ||
         SolnBlk.Grid.BCtypeS[i] == BC_WALL_VISCOUS_ISOTHERMAL ||
         SolnBlk.Grid.BCtypeS[i] == BC_TEMPERATURE_SLIP ||
         SolnBlk.Grid.BCtypeS[i] == BC_ADIABATIC_WALL) {
	n_pts = 5;
	i_index[0] = i-1; j_index[0] = j;
	i_index[1] = i-1; j_index[1] = j+1;
	i_index[2] = i;   j_index[2] = j+1;
	i_index[3] = i+1; j_index[3] = j+1;
	i_index[4] = i+1; j_index[4] = j;
      } else {
        n_pts = 8;
	i_index[0] = i-1; j_index[0] = j-1;
	i_index[1] = i-1; j_index[1] = j;
	i_index[2] = i-1; j_index[2] = j+1;
	i_index[3] = i;   j_index[3] = j+1;
	i_index[4] = i+1; j_index[4] = j+1;
	i_index[5] = i+1; j_index[5] = j;
	i_index[6] = i+1; j_index[6] = j-1;
	i_index[7] = i;   j_index[7] = j-1;
      } //endif
    
      ///////inside computational domain//////////

    } else {
      n_pts = 8;
      i_index[0] = i-1; j_index[0] = j-1;
      i_index[1] = i-1; j_index[1] = j;
      i_index[2] = i-1; j_index[2] = j+1;
      i_index[3] = i;   j_index[3] = j+1;
      i_index[4] = i+1; j_index[4] = j+1;
      i_index[5] = i+1; j_index[5] = j;
      i_index[6] = i+1; j_index[6] = j-1;
      i_index[7] = i;   j_index[7] = j-1;
    } //endif

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
           DUDx_ave = Gaussian2D_W_VACUUM;
           DUDy_ave = Gaussian2D_W_VACUUM;
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
           for ( n = 1 ; n <= NUM_VAR_GAUSSIAN2D ; ++n ) {
              u0Min = SolnBlk.W[i][j][n];
              u0Max = u0Min;
              for ( n2 = 0 ; n2 <= n_pts-1 ; ++n2 ) {
                 u0Min = min(u0Min, SolnBlk.W[ i_index[n2] ][ j_index[n2] ][n]);
                 u0Max = max(u0Max, SolnBlk.W[ i_index[n2] ][ j_index[n2] ][n]);
              } /* endfor */

	      ii=0;
	      if(i != SolnBlk.ICu || (SolnBlk.Grid.BCtypeE[j] != BC_ADIABATIC_WALL &&
				      SolnBlk.Grid.BCtypeE[j] != BC_WALL_VISCOUS_ISOTHERMAL &&
				      SolnBlk.Grid.BCtypeE[j] != BC_TEMPERATURE_SLIP)) {
		dX = SolnBlk.Grid.xfaceE(i, j)-SolnBlk.Grid.Cell[i][j].Xc;
		uQuad[ii] = SolnBlk.W[i][j][n] +
			    SolnBlk.dWdx[i][j][n]*dX.x +
			    SolnBlk.dWdy[i][j][n]*dX.y;
		ii++;
	      }
	      if(i != SolnBlk.ICl || (SolnBlk.Grid.BCtypeW[j] != BC_ADIABATIC_WALL &&
				      SolnBlk.Grid.BCtypeW[j] != BC_WALL_VISCOUS_ISOTHERMAL &&
				      SolnBlk.Grid.BCtypeW[j] != BC_TEMPERATURE_SLIP)) {
		dX = SolnBlk.Grid.xfaceW(i, j)-SolnBlk.Grid.Cell[i][j].Xc;
		uQuad[ii] = SolnBlk.W[i][j][n] +
                            SolnBlk.dWdx[i][j][n]*dX.x +
                            SolnBlk.dWdy[i][j][n]*dX.y;
		ii++;
	      }
	      if(j != SolnBlk.JCu || (SolnBlk.Grid.BCtypeN[i] != BC_ADIABATIC_WALL &&
				      SolnBlk.Grid.BCtypeN[i] != BC_WALL_VISCOUS_ISOTHERMAL &&
				      SolnBlk.Grid.BCtypeN[i] != BC_TEMPERATURE_SLIP)) {
		dX = SolnBlk.Grid.xfaceN(i, j)-SolnBlk.Grid.Cell[i][j].Xc;
		uQuad[ii] = SolnBlk.W[i][j][n] +
                           SolnBlk.dWdx[i][j][n]*dX.x +
		           SolnBlk.dWdy[i][j][n]*dX.y;
		ii++;
	      }
	      if(j != SolnBlk.JCl || (SolnBlk.Grid.BCtypeS[i] != BC_ADIABATIC_WALL &&
				      SolnBlk.Grid.BCtypeS[i] != BC_WALL_VISCOUS_ISOTHERMAL &&
				      SolnBlk.Grid.BCtypeS[i] != BC_TEMPERATURE_SLIP)) {
		dX = SolnBlk.Grid.xfaceS(i, j)-SolnBlk.Grid.Cell[i][j].Xc;
		uQuad[ii] = SolnBlk.W[i][j][n] +
                           SolnBlk.dWdx[i][j][n]*dX.x +
                           SolnBlk.dWdy[i][j][n]*dX.y;
		ii++;
	      }

	      switch(Limiter) {
	      case LIMITER_ONE :
		phi = ONE;
		break;
	      case LIMITER_ZERO :
		phi = ZERO;
		break;
	      case LIMITER_BARTH_JESPERSEN :
		phi = Limiter_BarthJespersen(uQuad, SolnBlk.W[i][j][n],
					     u0Min, u0Max, ii);
		break;
	      case LIMITER_VENKATAKRISHNAN :
		phi = Limiter_Venkatakrishnan(uQuad, SolnBlk.W[i][j][n],
					      u0Min, u0Max, ii);
		break;
	      case LIMITER_VANLEER :
		phi = Limiter_VanLeer(uQuad, SolnBlk.W[i][j][n],
				      u0Min, u0Max, ii);
		break;
	      case LIMITER_VANALBADA :
		phi = Limiter_VanAlbada(uQuad, SolnBlk.W[i][j][n],
					u0Min, u0Max, ii);
		break;
	      default:
		phi = Limiter_BarthJespersen(uQuad, SolnBlk.W[i][j][n],
					     u0Min, u0Max, ii);
		break;
	      } /* endswitch */
	      SolnBlk.phi[i][j][n] = phi;
           } /* endfor */

        } /* endif */

	//Adjust limiter for validity at Quadrature points
	/*
	dX = SolnBlk.Grid.xfaceE(i, j)-SolnBlk.Grid.Cell[i][j].Xc;
	W_temp = SolnBlk.W[i][j] + 
	         (SolnBlk.phi[i][j]^SolnBlk.dWdx[i][j])*dX.x +
	         (SolnBlk.phi[i][j]^SolnBlk.dWdy[i][j])*dX.y;
	if(W_temp.invalid()) {
	  if(W_temp.DetP() < SolnBlk.W[i][j].DetP()){
	    SolnBlk.phi[i][j] = Gaussian2D_W_VACUUM;
	  } else if (SolnBlk.W[i][j].DetP() > 0.0) {
	    a = (W_temp.p.xx-SolnBlk.W[i][j].p.xx) *
  	        (W_temp.p.yy-SolnBlk.W[i][j].p.yy) -
	        sqr(W_temp.p.xy-SolnBlk.W[i][j].p.xy);
	    b = (W_temp.p.xx-SolnBlk.W[i][j].p.xx)*SolnBlk.W[i][j].p.yy +
	        (W_temp.p.yy-SolnBlk.W[i][j].p.yy)*SolnBlk.W[i][j].p.xx -
	        2.0*(W_temp.p.xy-SolnBlk.W[i][j].p.xy)*SolnBlk.W[i][j].p.xy;
	    c = SolnBlk.W[i][j].p.xx*SolnBlk.W[i][j].p.yy-sqr(SolnBlk.W[i][j].p.xy);
	    d = -b+sqrt(sqr(b)-4.0*a*c)/(2.0*a);
	    cout << "(" << d << ")";
	    SolnBlk.phi[i][j][4] = SolnBlk.phi[i][j][4]*d;
	    SolnBlk.phi[i][j][5] = SolnBlk.phi[i][j][5]*d;
	    SolnBlk.phi[i][j][6] = SolnBlk.phi[i][j][6]*d;
	  }
	}

	dX = SolnBlk.Grid.xfaceW(i, j)-SolnBlk.Grid.Cell[i][j].Xc;
	W_temp = SolnBlk.W[i][j] + 
	         (SolnBlk.phi[i][j]^SolnBlk.dWdx[i][j])*dX.x +
	         (SolnBlk.phi[i][j]^SolnBlk.dWdy[i][j])*dX.y;
	if(W_temp.invalid()) {
	  if(W_temp.DetP() < SolnBlk.W[i][j].DetP()){
	    SolnBlk.phi[i][j] = Gaussian2D_W_VACUUM;
	  } else if (SolnBlk.W[i][j].DetP() > 0.0) {
	    a = (W_temp.p.xx-SolnBlk.W[i][j].p.xx) *
  	        (W_temp.p.yy-SolnBlk.W[i][j].p.yy) -
	        sqr(W_temp.p.xy-SolnBlk.W[i][j].p.xy);
	    b = (W_temp.p.xx-SolnBlk.W[i][j].p.xx)*SolnBlk.W[i][j].p.yy +
	        (W_temp.p.yy-SolnBlk.W[i][j].p.yy)*SolnBlk.W[i][j].p.xx -
	        2.0*(W_temp.p.xy-SolnBlk.W[i][j].p.xy)*SolnBlk.W[i][j].p.xy;
	    c = SolnBlk.W[i][j].p.xx*SolnBlk.W[i][j].p.yy-sqr(SolnBlk.W[i][j].p.xy);
	    d = -b+sqrt(sqr(b)-4.0*a*c)/(2.0*a);
	    cout << "(" << d << ")";
	    SolnBlk.phi[i][j][4] = SolnBlk.phi[i][j][4]*d;
	    SolnBlk.phi[i][j][5] = SolnBlk.phi[i][j][5]*d;
	    SolnBlk.phi[i][j][6] = SolnBlk.phi[i][j][6]*d;
	  }
	}

	dX = SolnBlk.Grid.xfaceN(i, j)-SolnBlk.Grid.Cell[i][j].Xc;
	W_temp = SolnBlk.W[i][j] + 
	         (SolnBlk.phi[i][j]^SolnBlk.dWdx[i][j])*dX.x +
	         (SolnBlk.phi[i][j]^SolnBlk.dWdy[i][j])*dX.y;
	if(W_temp.invalid()) {
	  if(W_temp.DetP() < SolnBlk.W[i][j].DetP()){
	    SolnBlk.phi[i][j] = Gaussian2D_W_VACUUM;
	  } else if (SolnBlk.W[i][j].DetP() > 0.0) {
	    a = (W_temp.p.xx-SolnBlk.W[i][j].p.xx) *
  	        (W_temp.p.yy-SolnBlk.W[i][j].p.yy) -
	        sqr(W_temp.p.xy-SolnBlk.W[i][j].p.xy);
	    b = (W_temp.p.xx-SolnBlk.W[i][j].p.xx)*SolnBlk.W[i][j].p.yy +
	        (W_temp.p.yy-SolnBlk.W[i][j].p.yy)*SolnBlk.W[i][j].p.xx -
	        2.0*(W_temp.p.xy-SolnBlk.W[i][j].p.xy)*SolnBlk.W[i][j].p.xy;
	    c = SolnBlk.W[i][j].p.xx*SolnBlk.W[i][j].p.yy-sqr(SolnBlk.W[i][j].p.xy);
	    d = -b+sqrt(sqr(b)-4.0*a*c)/(2.0*a);
	    cout << "(" << d << ")";
	    SolnBlk.phi[i][j][4] = SolnBlk.phi[i][j][4]*d;
	    SolnBlk.phi[i][j][5] = SolnBlk.phi[i][j][5]*d;
	    SolnBlk.phi[i][j][6] = SolnBlk.phi[i][j][6]*d;
	  }
	}


	dX = SolnBlk.Grid.xfaceS(i, j)-SolnBlk.Grid.Cell[i][j].Xc;
	W_temp = SolnBlk.W[i][j] + 
	         (SolnBlk.phi[i][j]^SolnBlk.dWdx[i][j])*dX.x +
	         (SolnBlk.phi[i][j]^SolnBlk.dWdy[i][j])*dX.y;
	if(W_temp.invalid()) {
	  if(W_temp.DetP() < SolnBlk.W[i][j].DetP()){
	    SolnBlk.phi[i][j] = Gaussian2D_W_VACUUM;
	  } else if (SolnBlk.W[i][j].DetP() > 0.0) {
	    a = (W_temp.p.xx-SolnBlk.W[i][j].p.xx) *
  	        (W_temp.p.yy-SolnBlk.W[i][j].p.yy) -
	        sqr(W_temp.p.xy-SolnBlk.W[i][j].p.xy);
	    b = (W_temp.p.xx-SolnBlk.W[i][j].p.xx)*SolnBlk.W[i][j].p.yy +
	        (W_temp.p.yy-SolnBlk.W[i][j].p.yy)*SolnBlk.W[i][j].p.xx -
	        2.0*(W_temp.p.xy-SolnBlk.W[i][j].p.xy)*SolnBlk.W[i][j].p.xy;
	    c = SolnBlk.W[i][j].p.xx*SolnBlk.W[i][j].p.yy-sqr(SolnBlk.W[i][j].p.xy);
	    d = -b+sqrt(sqr(b)-4.0*a*c)/(2.0*a);
	    cout << "(" << d << ")";
	    SolnBlk.phi[i][j][4] = SolnBlk.phi[i][j][4]*d;
	    SolnBlk.phi[i][j][5] = SolnBlk.phi[i][j][5]*d;
	    SolnBlk.phi[i][j][6] = SolnBlk.phi[i][j][6]*d;
	  }
	}
	*/

	dX = SolnBlk.Grid.xfaceE(i, j)-SolnBlk.Grid.Cell[i][j].Xc;
	while((SolnBlk.W[i][j] + 
	    (SolnBlk.phi[i][j]^SolnBlk.dWdx[i][j])*dX.x +
	    (SolnBlk.phi[i][j]^SolnBlk.dWdy[i][j])*dX.y).invalid()) {
	  cout << "^";
	  SolnBlk.phi[i][j] = SolnBlk.phi[i][j]*0.99;//Gaussian2D_W_VACUUM;
	}
	dX = SolnBlk.Grid.xfaceW(i, j)-SolnBlk.Grid.Cell[i][j].Xc;
	while((SolnBlk.W[i][j] + 
	    (SolnBlk.phi[i][j]^SolnBlk.dWdx[i][j])*dX.x +
	    (SolnBlk.phi[i][j]^SolnBlk.dWdy[i][j])*dX.y).invalid()) {
	  SolnBlk.phi[i][j] = SolnBlk.phi[i][j]*0.99;//Gaussian2D_W_VACUUM;
	  cout << "^";
	}
	dX = SolnBlk.Grid.xfaceN(i, j)-SolnBlk.Grid.Cell[i][j].Xc;
	while((SolnBlk.W[i][j] + 
	    (SolnBlk.phi[i][j]^SolnBlk.dWdx[i][j])*dX.x +
	    (SolnBlk.phi[i][j]^SolnBlk.dWdy[i][j])*dX.y).invalid()) {
	  SolnBlk.phi[i][j] = SolnBlk.phi[i][j]*0.99;//Gaussian2D_W_VACUUM;
	  cout << "^";
	}
	dX = SolnBlk.Grid.xfaceS(i, j)-SolnBlk.Grid.Cell[i][j].Xc;
	while((SolnBlk.W[i][j] + 
	    (SolnBlk.phi[i][j]^SolnBlk.dWdx[i][j])*dX.x +
	    (SolnBlk.phi[i][j]^SolnBlk.dWdy[i][j])*dX.y).invalid()) {
	  SolnBlk.phi[i][j] = SolnBlk.phi[i][j]*0.99;//Gaussian2D_W_VACUUM;
	  cout << "^";
	}
	

	
    } else {
        SolnBlk.dWdx[i][j] = Gaussian2D_W_VACUUM;
        SolnBlk.dWdy[i][j] = Gaussian2D_W_VACUUM; 
        SolnBlk.phi[i][j]  = Gaussian2D_W_VACUUM;
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
void Linear_Reconstruction_GreenGauss(Gaussian2D_Quad_Block &SolnBlk,
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
void Linear_Reconstruction_LeastSquares(Gaussian2D_Quad_Block &SolnBlk,
				        const int i, 
                                        const int j,
                                        const int Limiter) {

    int n, n2, n_pts, ii, i_index[8], j_index[8], reduction_counter(0);
    double u0Min, u0Max, uQuad[4], phi;
    double DxDx_ave, DxDy_ave, DyDy_ave;
    double a, b, c, d, e;
    Vector2D dX;
    Gaussian2D_pState DU, DUDx_ave, DUDy_ave, W_temp;

    /* Carry out the limited solution reconstruction in
       each cell of the computational mesh. */

      //If on outside ring of ghost cells, no reconstruction

    if (i == SolnBlk.ICl-SolnBlk.Nghost || i == SolnBlk.ICu+SolnBlk.Nghost ||
        j == SolnBlk.JCl-SolnBlk.Nghost || j == SolnBlk.JCu+SolnBlk.Nghost) {
      n_pts = 0;

      //If on inside ring of ghost cells, reconstruction may be required

      /////////////left edge//////////////////

//    } else if (i == SolnBlk.ICl-SolnBlk.Nghost+1) {
//        //bottom or top left corner (never used to evaluate fluxes, therefore
//        //no reconstruction needed
//      if(j == SolnBlk.JCl-SolnBlk.Nghost+1 || j == SolnBlk.JCu+SolnBlk.Nghost-1){
//	n_pts = 0;
//	//left edge, may need reconstruction
//        //(only for BC_NONE and BC_PERIODIC)???
//      } else if(SolnBlk.Grid.BCtypeW[j] == BC_FIXED ||
//                SolnBlk.Grid.BCtypeW[j] == BC_CHARACTERISTIC ||
//                SolnBlk.Grid.BCtypeW[j] == BC_LINEAR_EXTRAPOLATION ||
//                SolnBlk.Grid.BCtypeW[j] == BC_COUETTE ||
//                SolnBlk.Grid.BCtypeW[j] == BC_CONSTANT_EXTRAPOLATION ||
//                SolnBlk.Grid.BCtypeW[j] == BC_REFLECTION ||
//                SolnBlk.Grid.BCtypeW[j] == BC_WALL_VISCOUS_ISOTHERMAL ||
//                SolnBlk.Grid.BCtypeW[j] == BC_TEMPERATURE_SLIP ||
//                SolnBlk.Grid.BCtypeW[j] == BC_ADIABATIC_WALL) {
//	n_pts = 0;
//      } else { //BC's are Periodic or "none"
//	  //Bottom left
//	if(j == SolnBlk.JCl-SolnBlk.Nghost+2 &&
//               (SolnBlk.Grid.BCtypeS[i] == BC_FIXED ||
//                SolnBlk.Grid.BCtypeS[i] == BC_CHARACTERISTIC ||
//                SolnBlk.Grid.BCtypeS[i] == BC_LINEAR_EXTRAPOLATION ||
//                SolnBlk.Grid.BCtypeS[i] == BC_COUETTE              ||
//                SolnBlk.Grid.BCtypeS[i] == BC_CONSTANT_EXTRAPOLATION ||
//                SolnBlk.Grid.BCtypeS[i] == BC_REFLECTION ||
//                SolnBlk.Grid.BCtypeS[i] == BC_WALL_VISCOUS_ISOTHERMAL ||
//                SolnBlk.Grid.BCtypeS[i] == BC_TEMPERATURE_SLIP ||
//                SolnBlk.Grid.BCtypeS[i] == BC_ADIABATIC_WALL)) {
//	  n_pts = 5;
//	  i_index[0] = i-1; j_index[0] = j;
//	  i_index[1] = i-1; j_index[1] = j+1;
//	  i_index[2] = i;   j_index[2] = j+1;
//	  i_index[3] = i+1; j_index[3] = j+1;
//	  i_index[4] = i+1; j_index[4] = j;
//	  //Top Left
//	} else if (j == SolnBlk.JCu+SolnBlk.Nghost-2 &&
//                       (SolnBlk.Grid.BCtypeN[i] == BC_FIXED ||
//                        SolnBlk.Grid.BCtypeN[i] == BC_CHARACTERISTIC ||
//                        SolnBlk.Grid.BCtypeN[i] == BC_LINEAR_EXTRAPOLATION ||
//                        SolnBlk.Grid.BCtypeN[i] == BC_COUETTE ||
//                        SolnBlk.Grid.BCtypeN[i] == BC_CONSTANT_EXTRAPOLATION ||
//                        SolnBlk.Grid.BCtypeN[i] == BC_REFLECTION ||
//                        SolnBlk.Grid.BCtypeN[i] == BC_WALL_VISCOUS_ISOTHERMAL ||
//                        SolnBlk.Grid.BCtypeN[i] == BC_TEMPERATURE_SLIP ||
//                        SolnBlk.Grid.BCtypeN[i] == BC_ADIABATIC_WALL)) {
//	  n_pts = 5;
//	  i_index[0] = i-1; j_index[0] = j;
//	  i_index[1] = i-1; j_index[1] = j-1;
//	  i_index[2] = i;   j_index[2] = j-1;
//	  i_index[3] = i+1; j_index[3] = j-1;
//	  i_index[4] = i+1; j_index[4] = j;
//	  //middle left
//	} else {
//	  n_pts = 8;
//	  i_index[0] = i-1; j_index[0] = j-1;
//	  i_index[1] = i-1; j_index[1] = j;
//	  i_index[2] = i-1; j_index[2] = j+1;
//	  i_index[3] = i;   j_index[3] = j+1;
//	  i_index[4] = i+1; j_index[4] = j+1;
//	  i_index[5] = i+1; j_index[5] = j;
//	  i_index[6] = i+1; j_index[6] = j-1;
//	  i_index[7] = i;   j_index[7] = j-1;
//	} //endif
//      } //endif
//
//      /////////////Right edge//////////////////
//
//    } else if (i == SolnBlk.ICu+SolnBlk.Nghost-1) {
//        //bottom or top left corner (never used to evaluate fluxes, therefore
//        //no reconstruction needed
//      if(j == SolnBlk.JCl-SolnBlk.Nghost+1 || j == SolnBlk.JCu+SolnBlk.Nghost-1){
//	n_pts = 0;
//	//right edge, may need reconstruction
//        //(only for BC_NONE and BC_PERIODIC)???
//      } else if(SolnBlk.Grid.BCtypeE[j] == BC_FIXED ||
//                SolnBlk.Grid.BCtypeE[j] == BC_CHARACTERISTIC ||
//                SolnBlk.Grid.BCtypeE[j] == BC_LINEAR_EXTRAPOLATION ||
//                SolnBlk.Grid.BCtypeE[j] == BC_COUETTE ||
//                SolnBlk.Grid.BCtypeE[j] == BC_CONSTANT_EXTRAPOLATION ||
//                SolnBlk.Grid.BCtypeE[j] == BC_REFLECTION ||
//                SolnBlk.Grid.BCtypeE[j] == BC_WALL_VISCOUS_ISOTHERMAL ||
//                SolnBlk.Grid.BCtypeE[j] == BC_TEMPERATURE_SLIP ||
//                SolnBlk.Grid.BCtypeE[j] == BC_ADIABATIC_WALL) {
//	n_pts = 0;
//      } else { //BC's are Periodic or "none"
//	  //Bottom right
//	if(j == SolnBlk.JCl-SolnBlk.Nghost+2 &&
//               (SolnBlk.Grid.BCtypeS[i] == BC_FIXED ||
//                SolnBlk.Grid.BCtypeS[i] == BC_CHARACTERISTIC ||
//                SolnBlk.Grid.BCtypeS[i] == BC_LINEAR_EXTRAPOLATION ||
//                SolnBlk.Grid.BCtypeS[i] == BC_COUETTE ||
//                SolnBlk.Grid.BCtypeS[i] == BC_CONSTANT_EXTRAPOLATION ||
//                SolnBlk.Grid.BCtypeS[i] == BC_REFLECTION ||
//                SolnBlk.Grid.BCtypeS[i] == BC_WALL_VISCOUS_ISOTHERMAL ||
//                SolnBlk.Grid.BCtypeS[i] == BC_TEMPERATURE_SLIP ||
//                SolnBlk.Grid.BCtypeS[i] == BC_ADIABATIC_WALL)) {
//	  n_pts = 5;
//	  i_index[0] = i-1; j_index[0] = j;
//	  i_index[1] = i-1; j_index[1] = j+1;
//	  i_index[2] = i;   j_index[2] = j+1;
//	  i_index[3] = i+1; j_index[3] = j+1;
//	  i_index[4] = i+1; j_index[4] = j;
//	  //Top right
//	} else if (j == SolnBlk.JCu+SolnBlk.Nghost-2 &&
//                       (SolnBlk.Grid.BCtypeN[i] == BC_FIXED ||
//                        SolnBlk.Grid.BCtypeN[i] == BC_CHARACTERISTIC ||
//                        SolnBlk.Grid.BCtypeN[i] == BC_LINEAR_EXTRAPOLATION ||
//                        SolnBlk.Grid.BCtypeN[i] == BC_COUETTE              ||
//                        SolnBlk.Grid.BCtypeN[i] == BC_CONSTANT_EXTRAPOLATION ||
//                        SolnBlk.Grid.BCtypeN[i] == BC_REFLECTION ||
//                        SolnBlk.Grid.BCtypeN[i] == BC_WALL_VISCOUS_ISOTHERMAL ||
//                        SolnBlk.Grid.BCtypeN[i] == BC_TEMPERATURE_SLIP ||
//                        SolnBlk.Grid.BCtypeN[i] == BC_ADIABATIC_WALL)) {
//	  n_pts = 5;
//	  i_index[0] = i-1; j_index[0] = j;
//	  i_index[1] = i-1; j_index[1] = j-1;
//	  i_index[2] = i;   j_index[2] = j-1;
//	  i_index[3] = i+1; j_index[3] = j-1;
//	  i_index[4] = i+1; j_index[4] = j;
//	  //middle right or top or bottom with none or priodic bc's
//	} else {
//	  n_pts = 8;
//	  i_index[0] = i-1; j_index[0] = j-1;
//	  i_index[1] = i-1; j_index[1] = j;
//	  i_index[2] = i-1; j_index[2] = j+1;
//	  i_index[3] = i;   j_index[3] = j+1;
//	  i_index[4] = i+1; j_index[4] = j+1;
//	  i_index[5] = i+1; j_index[5] = j;
//	  i_index[6] = i+1; j_index[6] = j-1;
//	  i_index[7] = i;   j_index[7] = j-1;
//	} //endif
//      } //endif
//
//      /////////////bottom edge/////////////
//      // (corners already taken care of) //
//
//    } else if (j == SolnBlk.JCl-SolnBlk.Nghost+1) {
//      if(SolnBlk.Grid.BCtypeS[i] == BC_FIXED ||
//         SolnBlk.Grid.BCtypeS[i] == BC_CHARACTERISTIC ||
//         SolnBlk.Grid.BCtypeS[i] == BC_LINEAR_EXTRAPOLATION ||
//         SolnBlk.Grid.BCtypeS[i] == BC_COUETTE              ||
//         SolnBlk.Grid.BCtypeS[i] == BC_CONSTANT_EXTRAPOLATION ||
//         SolnBlk.Grid.BCtypeS[i] == BC_REFLECTION ||
//         SolnBlk.Grid.BCtypeS[i] == BC_WALL_VISCOUS_ISOTHERMAL ||
//         SolnBlk.Grid.BCtypeS[i] == BC_TEMPERATURE_SLIP ||
//         SolnBlk.Grid.BCtypeS[i] == BC_ADIABATIC_WALL) {
//	n_pts = 0;
//      } else { //BC's are Periodic or "none"
//	  //Bottom left
//	if(i == SolnBlk.ICl-SolnBlk.Nghost+2 &&
//               (SolnBlk.Grid.BCtypeW[j] == BC_FIXED ||
//                SolnBlk.Grid.BCtypeW[j] == BC_CHARACTERISTIC ||
//                SolnBlk.Grid.BCtypeW[j] == BC_LINEAR_EXTRAPOLATION ||
//                SolnBlk.Grid.BCtypeW[j] == BC_COUETTE              ||
//                SolnBlk.Grid.BCtypeW[j] == BC_CONSTANT_EXTRAPOLATION ||
//                SolnBlk.Grid.BCtypeW[j] == BC_REFLECTION ||
//                SolnBlk.Grid.BCtypeW[j] == BC_WALL_VISCOUS_ISOTHERMAL ||
//                SolnBlk.Grid.BCtypeW[j] == BC_TEMPERATURE_SLIP ||
//                SolnBlk.Grid.BCtypeW[j] == BC_ADIABATIC_WALL)) {
//	  n_pts = 5;
//	  i_index[0] = i;   j_index[0] = j+1;
//	  i_index[1] = i+1; j_index[1] = j+1;
//	  i_index[2] = i+1; j_index[2] = j;
//	  i_index[3] = i+1; j_index[3] = j-1;
//	  i_index[4] = i;   j_index[4] = j-1;
//	  //Bottom right
//	} else if (i == SolnBlk.ICu+SolnBlk.Nghost-2 &&
//                       (SolnBlk.Grid.BCtypeE[j] == BC_FIXED ||
//                        SolnBlk.Grid.BCtypeE[j] == BC_CHARACTERISTIC ||
//                        SolnBlk.Grid.BCtypeE[j] == BC_LINEAR_EXTRAPOLATION ||
//                        SolnBlk.Grid.BCtypeE[j] == BC_COUETTE ||
//                        SolnBlk.Grid.BCtypeE[j] == BC_CONSTANT_EXTRAPOLATION ||
//                        SolnBlk.Grid.BCtypeE[j] == BC_REFLECTION ||
//                        SolnBlk.Grid.BCtypeE[j] == BC_WALL_VISCOUS_ISOTHERMAL ||
//                        SolnBlk.Grid.BCtypeE[j] == BC_TEMPERATURE_SLIP ||
//                        SolnBlk.Grid.BCtypeE[j] == BC_ADIABATIC_WALL)) {
//	  n_pts = 5;
//	  i_index[0] = i;   j_index[0] = j+1;
//	  i_index[1] = i-1; j_index[1] = j+1;
//	  i_index[2] = i-1; j_index[2] = j;
//	  i_index[3] = i-1; j_index[3] = j-1;
//	  i_index[4] = i;   j_index[4] = j-1;
//	  //middle bottom or left or right with none or priodic bc's
//	} else {
//	  n_pts = 8;
//	  i_index[0] = i-1; j_index[0] = j-1;
//	  i_index[1] = i-1; j_index[1] = j;
//	  i_index[2] = i-1; j_index[2] = j+1;
//	  i_index[3] = i;   j_index[3] = j+1;
//	  i_index[4] = i+1; j_index[4] = j+1;
//	  i_index[5] = i+1; j_index[5] = j;
//	  i_index[6] = i+1; j_index[6] = j-1;
//	  i_index[7] = i;   j_index[7] = j-1;
//	} //endif
//      } //endif
//
//      /////////////top edge/////////////
//      // (corners already taken care of) //
//
//    } else if (j == SolnBlk.JCu+SolnBlk.Nghost-1) {
//      if(SolnBlk.Grid.BCtypeN[i] == BC_FIXED ||
//         SolnBlk.Grid.BCtypeN[i] == BC_CHARACTERISTIC ||
//         SolnBlk.Grid.BCtypeN[i] == BC_LINEAR_EXTRAPOLATION ||
//         SolnBlk.Grid.BCtypeN[i] == BC_COUETTE ||
//         SolnBlk.Grid.BCtypeN[i] == BC_CONSTANT_EXTRAPOLATION ||
//         SolnBlk.Grid.BCtypeN[i] == BC_REFLECTION ||
//         SolnBlk.Grid.BCtypeN[i] == BC_WALL_VISCOUS_ISOTHERMAL ||
//         SolnBlk.Grid.BCtypeN[i] == BC_TEMPERATURE_SLIP ||
//         SolnBlk.Grid.BCtypeN[i] == BC_ADIABATIC_WALL) {
//	n_pts = 0;
//      } else { //BC's are Periodic or "none"
//	  //Top left
//	if(i == SolnBlk.ICl-SolnBlk.Nghost+2 &&
//               (SolnBlk.Grid.BCtypeW[j] == BC_FIXED ||
//                SolnBlk.Grid.BCtypeW[j] == BC_CHARACTERISTIC ||
//                SolnBlk.Grid.BCtypeW[j] == BC_LINEAR_EXTRAPOLATION ||
//                SolnBlk.Grid.BCtypeW[j] == BC_COUETTE ||
//                SolnBlk.Grid.BCtypeW[j] == BC_CONSTANT_EXTRAPOLATION ||
//                SolnBlk.Grid.BCtypeW[j] == BC_REFLECTION ||
//                SolnBlk.Grid.BCtypeW[j] == BC_WALL_VISCOUS_ISOTHERMAL ||
//                SolnBlk.Grid.BCtypeW[j] == BC_TEMPERATURE_SLIP ||
//                SolnBlk.Grid.BCtypeW[j] == BC_ADIABATIC_WALL)) {
//	  n_pts = 5;
//	  i_index[0] = i;   j_index[0] = j+1;
//	  i_index[1] = i+1; j_index[1] = j+1;
//	  i_index[2] = i+1; j_index[2] = j;
//	  i_index[3] = i+1; j_index[3] = j-1;
//	  i_index[4] = i;   j_index[4] = j-1;
//	  //Top right
//	} else if (i == SolnBlk.ICu+SolnBlk.Nghost-2 &&
//                       (SolnBlk.Grid.BCtypeE[j] == BC_FIXED ||
//                        SolnBlk.Grid.BCtypeE[j] == BC_CHARACTERISTIC ||
//                        SolnBlk.Grid.BCtypeE[j] == BC_LINEAR_EXTRAPOLATION ||
//                        SolnBlk.Grid.BCtypeE[j] == BC_COUETTE ||
//                        SolnBlk.Grid.BCtypeE[j] == BC_CONSTANT_EXTRAPOLATION ||
//                        SolnBlk.Grid.BCtypeE[j] == BC_REFLECTION ||
//                        SolnBlk.Grid.BCtypeE[j] == BC_WALL_VISCOUS_ISOTHERMAL ||
//                        SolnBlk.Grid.BCtypeE[j] == BC_TEMPERATURE_SLIP ||
//                        SolnBlk.Grid.BCtypeE[j] == BC_ADIABATIC_WALL)) {
//	  n_pts = 5;
//	  i_index[0] = i;   j_index[0] = j+1;
//	  i_index[1] = i-1; j_index[1] = j+1;
//	  i_index[2] = i-1; j_index[2] = j;
//	  i_index[3] = i-1; j_index[3] = j-1;
//	  i_index[4] = i;   j_index[4] = j-1;
//	  //top bottom or left or right with none or priodic bc's
//	} else {
//	  n_pts = 8;
//	  i_index[0] = i-1; j_index[0] = j-1;
//	  i_index[1] = i-1; j_index[1] = j;
//	  i_index[2] = i-1; j_index[2] = j+1;
//	  i_index[3] = i;   j_index[3] = j+1;
//	  i_index[4] = i+1; j_index[4] = j+1;
//	  i_index[5] = i+1; j_index[5] = j;
//	  i_index[6] = i+1; j_index[6] = j-1;
//	  i_index[7] = i;   j_index[7] = j-1;
//	} //endif
//      } //endif
//    
//      /////////First Cell Inside Computational Domain////////
//
//      //bottom left
//    } else if (i == SolnBlk.ICl-SolnBlk.Nghost+2 && 
//	       j == SolnBlk.JCl-SolnBlk.Nghost+2) {
//      if((SolnBlk.Grid.BCtypeW[j] == BC_CHARACTERISTIC ||
//          SolnBlk.Grid.BCtypeW[j] == BC_LINEAR_EXTRAPOLATION ||
//          SolnBlk.Grid.BCtypeW[j] == BC_COUETTE ||
//          SolnBlk.Grid.BCtypeW[j] == BC_CONSTANT_EXTRAPOLATION ||
//          SolnBlk.Grid.BCtypeW[j] == BC_REFLECTION ||
//          SolnBlk.Grid.BCtypeW[j] == BC_WALL_VISCOUS_ISOTHERMAL ||
//          SolnBlk.Grid.BCtypeW[j] == BC_TEMPERATURE_SLIP ||
//          SolnBlk.Grid.BCtypeW[j] == BC_ADIABATIC_WALL)
//                           &&
//         (SolnBlk.Grid.BCtypeS[i] == BC_CHARACTERISTIC ||
//          SolnBlk.Grid.BCtypeS[i] == BC_LINEAR_EXTRAPOLATION ||
//          SolnBlk.Grid.BCtypeS[i] == BC_COUETTE ||
//          SolnBlk.Grid.BCtypeS[i] == BC_CONSTANT_EXTRAPOLATION ||
//          SolnBlk.Grid.BCtypeS[i] == BC_REFLECTION ||
//          SolnBlk.Grid.BCtypeS[i] == BC_WALL_VISCOUS_ISOTHERMAL ||
//          SolnBlk.Grid.BCtypeS[i] == BC_TEMPERATURE_SLIP ||
//          SolnBlk.Grid.BCtypeS[i] == BC_ADIABATIC_WALL)) {
//	n_pts = 3;
//	i_index[0] = i;   j_index[0] = j+1;
//	i_index[1] = i+1; j_index[1] = j+1;
//	i_index[2] = i+1; j_index[2] = j;
//      } else if (SolnBlk.Grid.BCtypeS[i] == BC_CHARACTERISTIC ||
//                 SolnBlk.Grid.BCtypeS[i] == BC_LINEAR_EXTRAPOLATION ||
//                 SolnBlk.Grid.BCtypeS[i] == BC_COUETTE ||
//                 SolnBlk.Grid.BCtypeS[i] == BC_CONSTANT_EXTRAPOLATION ||
//                 SolnBlk.Grid.BCtypeS[i] == BC_REFLECTION ||
//                 SolnBlk.Grid.BCtypeS[i] == BC_WALL_VISCOUS_ISOTHERMAL ||
//                 SolnBlk.Grid.BCtypeS[i] == BC_TEMPERATURE_SLIP ||
//                 SolnBlk.Grid.BCtypeS[i] == BC_ADIABATIC_WALL) {
//	n_pts = 5;
//	i_index[0] = i-1; j_index[0] = j;
//	i_index[1] = i-1; j_index[1] = j+1;
//	i_index[2] = i;   j_index[2] = j+1;
//	i_index[3] = i+1; j_index[3] = j+1;
//	i_index[4] = i+1; j_index[4] = j;
//      } else if (SolnBlk.Grid.BCtypeW[j] == BC_CHARACTERISTIC ||
//                 SolnBlk.Grid.BCtypeW[j] == BC_LINEAR_EXTRAPOLATION ||
//                 SolnBlk.Grid.BCtypeW[j] == BC_COUETTE ||
//                 SolnBlk.Grid.BCtypeW[j] == BC_CONSTANT_EXTRAPOLATION ||
//                 SolnBlk.Grid.BCtypeW[j] == BC_REFLECTION ||
//                 SolnBlk.Grid.BCtypeW[j] == BC_WALL_VISCOUS_ISOTHERMAL ||
//                 SolnBlk.Grid.BCtypeW[j] == BC_TEMPERATURE_SLIP ||
//                 SolnBlk.Grid.BCtypeW[j] == BC_ADIABATIC_WALL) {
//	n_pts = 5;
//	i_index[0] = i;   j_index[0] = j+1;
//	i_index[1] = i+1; j_index[1] = j+1;
//	i_index[2] = i+1; j_index[2] = j;
//	i_index[3] = i+1; j_index[3] = j-1;
//	i_index[4] = i;   j_index[4] = j-1;
//      } else {
//        n_pts = 8;
//	i_index[0] = i-1; j_index[0] = j-1;
//	i_index[1] = i-1; j_index[1] = j;
//	i_index[2] = i-1; j_index[2] = j+1;
//	i_index[3] = i;   j_index[3] = j+1;
//	i_index[4] = i+1; j_index[4] = j+1;
//	i_index[5] = i+1; j_index[5] = j;
//	i_index[6] = i+1; j_index[6] = j-1;
//	i_index[7] = i;   j_index[7] = j-1;
//      } //endif
//
//      //top left
//    } else if (i == SolnBlk.ICl-SolnBlk.Nghost+2 && 
//	       j == SolnBlk.JCu+SolnBlk.Nghost-2) {
//      if((SolnBlk.Grid.BCtypeW[j] == BC_CHARACTERISTIC ||
//          SolnBlk.Grid.BCtypeW[j] == BC_LINEAR_EXTRAPOLATION ||
//          SolnBlk.Grid.BCtypeW[j] == BC_COUETTE ||
//          SolnBlk.Grid.BCtypeW[j] == BC_CONSTANT_EXTRAPOLATION ||
//          SolnBlk.Grid.BCtypeW[j] == BC_REFLECTION ||
//          SolnBlk.Grid.BCtypeW[j] == BC_WALL_VISCOUS_ISOTHERMAL ||
//          SolnBlk.Grid.BCtypeW[j] == BC_TEMPERATURE_SLIP ||
//          SolnBlk.Grid.BCtypeW[j] == BC_ADIABATIC_WALL)
//                           &&
//         (SolnBlk.Grid.BCtypeN[i] == BC_CHARACTERISTIC ||
//          SolnBlk.Grid.BCtypeN[i] == BC_LINEAR_EXTRAPOLATION ||
//          SolnBlk.Grid.BCtypeN[i] == BC_COUETTE ||
//          SolnBlk.Grid.BCtypeN[i] == BC_CONSTANT_EXTRAPOLATION ||
//          SolnBlk.Grid.BCtypeN[i] == BC_REFLECTION ||
//          SolnBlk.Grid.BCtypeN[i] == BC_WALL_VISCOUS_ISOTHERMAL ||
//          SolnBlk.Grid.BCtypeN[i] == BC_TEMPERATURE_SLIP ||
//          SolnBlk.Grid.BCtypeN[i] == BC_ADIABATIC_WALL)) {
//	n_pts = 3;
//	i_index[0] = i;   j_index[0] = j-1;
//	i_index[1] = i+1; j_index[1] = j-1;
//	i_index[2] = i+1; j_index[2] = j;
//      } else if (SolnBlk.Grid.BCtypeN[i] == BC_CHARACTERISTIC ||
//                 SolnBlk.Grid.BCtypeN[i] == BC_LINEAR_EXTRAPOLATION ||
//                 SolnBlk.Grid.BCtypeN[i] == BC_COUETTE ||
//                 SolnBlk.Grid.BCtypeN[i] == BC_CONSTANT_EXTRAPOLATION ||
//                 SolnBlk.Grid.BCtypeN[i] == BC_REFLECTION ||
//                 SolnBlk.Grid.BCtypeN[i] == BC_WALL_VISCOUS_ISOTHERMAL ||
//                 SolnBlk.Grid.BCtypeN[i] == BC_TEMPERATURE_SLIP ||
//                 SolnBlk.Grid.BCtypeN[i] == BC_ADIABATIC_WALL) {
//	n_pts = 5;
//	i_index[0] = i-1; j_index[0] = j;
//	i_index[1] = i-1; j_index[1] = j-1;
//	i_index[2] = i;   j_index[2] = j-1;
//	i_index[3] = i+1; j_index[3] = j-1;
//	i_index[4] = i+1; j_index[4] = j;
//      } else if (SolnBlk.Grid.BCtypeW[j] == BC_CHARACTERISTIC ||
//                 SolnBlk.Grid.BCtypeW[j] == BC_LINEAR_EXTRAPOLATION ||
//                 SolnBlk.Grid.BCtypeW[j] == BC_COUETTE ||
//                 SolnBlk.Grid.BCtypeW[j] == BC_CONSTANT_EXTRAPOLATION ||
//                 SolnBlk.Grid.BCtypeW[j] == BC_REFLECTION ||
//                 SolnBlk.Grid.BCtypeW[j] == BC_WALL_VISCOUS_ISOTHERMAL ||
//                 SolnBlk.Grid.BCtypeW[j] == BC_TEMPERATURE_SLIP ||
//                 SolnBlk.Grid.BCtypeW[j] == BC_ADIABATIC_WALL) {
//	n_pts = 5;
//	i_index[0] = i;   j_index[0] = j+1;
//	i_index[1] = i+1; j_index[1] = j+1;
//	i_index[2] = i+1; j_index[2] = j;
//	i_index[3] = i+1; j_index[3] = j-1;
//	i_index[4] = i;   j_index[4] = j-1;
//      } else {
//        n_pts = 8;
//	i_index[0] = i-1; j_index[0] = j-1;
//	i_index[1] = i-1; j_index[1] = j;
//	i_index[2] = i-1; j_index[2] = j+1;
//	i_index[3] = i;   j_index[3] = j+1;
//	i_index[4] = i+1; j_index[4] = j+1;
//	i_index[5] = i+1; j_index[5] = j;
//	i_index[6] = i+1; j_index[6] = j-1;
//	i_index[7] = i;   j_index[7] = j-1;
//      } //endif
//
//      //top right
//    } else if (i == SolnBlk.ICu+SolnBlk.Nghost-2 && 
//	       j == SolnBlk.JCu+SolnBlk.Nghost-2) {
//      if((SolnBlk.Grid.BCtypeE[j] == BC_CHARACTERISTIC ||
//          SolnBlk.Grid.BCtypeE[j] == BC_LINEAR_EXTRAPOLATION ||
//          SolnBlk.Grid.BCtypeE[j] == BC_COUETTE ||
//          SolnBlk.Grid.BCtypeE[j] == BC_CONSTANT_EXTRAPOLATION ||
//          SolnBlk.Grid.BCtypeE[j] == BC_REFLECTION ||
//          SolnBlk.Grid.BCtypeE[j] == BC_WALL_VISCOUS_ISOTHERMAL ||
//          SolnBlk.Grid.BCtypeE[j] == BC_TEMPERATURE_SLIP ||
//          SolnBlk.Grid.BCtypeE[j] == BC_ADIABATIC_WALL)
//                           &&
//         (SolnBlk.Grid.BCtypeN[i] == BC_CHARACTERISTIC ||
//          SolnBlk.Grid.BCtypeN[i] == BC_LINEAR_EXTRAPOLATION ||
//          SolnBlk.Grid.BCtypeN[i] == BC_COUETTE ||
//          SolnBlk.Grid.BCtypeN[i] == BC_CONSTANT_EXTRAPOLATION ||
//          SolnBlk.Grid.BCtypeN[i] == BC_REFLECTION ||
//          SolnBlk.Grid.BCtypeN[i] == BC_WALL_VISCOUS_ISOTHERMAL ||
//          SolnBlk.Grid.BCtypeN[i] == BC_TEMPERATURE_SLIP ||
//          SolnBlk.Grid.BCtypeN[i] == BC_ADIABATIC_WALL)) {
//	n_pts = 3;
//	i_index[0] = i-1; j_index[0] = j;
//	i_index[1] = i-1; j_index[1] = j-1;
//	i_index[2] = i;   j_index[2] = j-1;
//      } else if (SolnBlk.Grid.BCtypeN[i] == BC_CHARACTERISTIC ||
//                 SolnBlk.Grid.BCtypeN[i] == BC_LINEAR_EXTRAPOLATION ||
//                 SolnBlk.Grid.BCtypeN[i] == BC_COUETTE ||
//                 SolnBlk.Grid.BCtypeN[i] == BC_CONSTANT_EXTRAPOLATION ||
//                 SolnBlk.Grid.BCtypeN[i] == BC_REFLECTION ||
//                 SolnBlk.Grid.BCtypeN[i] == BC_WALL_VISCOUS_ISOTHERMAL ||
//                 SolnBlk.Grid.BCtypeN[i] == BC_TEMPERATURE_SLIP ||
//                 SolnBlk.Grid.BCtypeN[i] == BC_ADIABATIC_WALL) {
//	n_pts = 5;
//	i_index[0] = i-1; j_index[0] = j;
//	i_index[1] = i-1; j_index[1] = j-1;
//	i_index[2] = i;   j_index[2] = j-1;
//	i_index[3] = i+1; j_index[3] = j-1;
//	i_index[4] = i+1; j_index[4] = j;
//      } else if (SolnBlk.Grid.BCtypeE[j] == BC_CHARACTERISTIC ||
//                 SolnBlk.Grid.BCtypeE[j] == BC_LINEAR_EXTRAPOLATION ||
//                 SolnBlk.Grid.BCtypeE[j] == BC_COUETTE ||
//                 SolnBlk.Grid.BCtypeE[j] == BC_CONSTANT_EXTRAPOLATION ||
//                 SolnBlk.Grid.BCtypeE[j] == BC_REFLECTION ||
//                 SolnBlk.Grid.BCtypeE[j] == BC_WALL_VISCOUS_ISOTHERMAL ||
//                 SolnBlk.Grid.BCtypeE[j] == BC_TEMPERATURE_SLIP ||
//                 SolnBlk.Grid.BCtypeE[j] == BC_ADIABATIC_WALL) {
//	n_pts = 5;
//	i_index[0] = i;   j_index[0] = j+1;
//	i_index[1] = i-1; j_index[1] = j+1;
//	i_index[2] = i-1; j_index[2] = j;
//	i_index[3] = i-1; j_index[3] = j-1;
//	i_index[4] = i;   j_index[4] = j-1;
//      } else {
//        n_pts = 8;
//	i_index[0] = i-1; j_index[0] = j-1;
//	i_index[1] = i-1; j_index[1] = j;
//	i_index[2] = i-1; j_index[2] = j+1;
//	i_index[3] = i;   j_index[3] = j+1;
//	i_index[4] = i+1; j_index[4] = j+1;
//	i_index[5] = i+1; j_index[5] = j;
//	i_index[6] = i+1; j_index[6] = j-1;
//	i_index[7] = i;   j_index[7] = j-1;
//      } //endif
//
//      //Bottom right
//    } else if (i == SolnBlk.ICu+SolnBlk.Nghost-2 && 
//	       j == SolnBlk.JCl-SolnBlk.Nghost+2) {
//      if((SolnBlk.Grid.BCtypeE[j] == BC_CHARACTERISTIC ||
//          SolnBlk.Grid.BCtypeE[j] == BC_LINEAR_EXTRAPOLATION ||
//          SolnBlk.Grid.BCtypeE[j] == BC_COUETTE ||
//          SolnBlk.Grid.BCtypeE[j] == BC_CONSTANT_EXTRAPOLATION ||
//          SolnBlk.Grid.BCtypeE[j] == BC_REFLECTION ||
//	  SolnBlk.Grid.BCtypeE[j] == BC_WALL_VISCOUS_ISOTHERMAL ||
//	  SolnBlk.Grid.BCtypeE[j] == BC_TEMPERATURE_SLIP ||
//          SolnBlk.Grid.BCtypeE[j] == BC_ADIABATIC_WALL)
//                           &&
//         (SolnBlk.Grid.BCtypeS[i] == BC_CHARACTERISTIC ||
//          SolnBlk.Grid.BCtypeS[i] == BC_LINEAR_EXTRAPOLATION ||
//          SolnBlk.Grid.BCtypeS[i] == BC_COUETTE ||
//          SolnBlk.Grid.BCtypeS[i] == BC_CONSTANT_EXTRAPOLATION ||
//          SolnBlk.Grid.BCtypeS[i] == BC_REFLECTION ||
//          SolnBlk.Grid.BCtypeS[i] == BC_WALL_VISCOUS_ISOTHERMAL ||
//          SolnBlk.Grid.BCtypeS[i] == BC_TEMPERATURE_SLIP ||
//          SolnBlk.Grid.BCtypeS[i] == BC_ADIABATIC_WALL)) {
//	n_pts = 3;
//	i_index[0] = i;   j_index[0] = j+1;
//	i_index[1] = i-1; j_index[1] = j+1;
//	i_index[2] = i-1; j_index[2] = j;
//      } else if (SolnBlk.Grid.BCtypeS[i] == BC_CHARACTERISTIC ||
//                 SolnBlk.Grid.BCtypeS[i] == BC_LINEAR_EXTRAPOLATION ||
//                 SolnBlk.Grid.BCtypeS[i] == BC_COUETTE ||
//                 SolnBlk.Grid.BCtypeS[i] == BC_CONSTANT_EXTRAPOLATION ||
//                 SolnBlk.Grid.BCtypeS[i] == BC_REFLECTION ||
//                 SolnBlk.Grid.BCtypeS[i] == BC_WALL_VISCOUS_ISOTHERMAL ||
//                 SolnBlk.Grid.BCtypeS[i] == BC_TEMPERATURE_SLIP ||
//                 SolnBlk.Grid.BCtypeS[i] == BC_ADIABATIC_WALL) {
//	n_pts = 5;
//	i_index[0] = i-1; j_index[0] = j;
//	i_index[1] = i-1; j_index[1] = j+1;
//	i_index[2] = i;   j_index[2] = j+1;
//	i_index[3] = i+1; j_index[3] = j+1;
//	i_index[4] = i+1; j_index[4] = j;
//      } else if (SolnBlk.Grid.BCtypeE[j] == BC_CHARACTERISTIC ||
//                 SolnBlk.Grid.BCtypeE[j] == BC_LINEAR_EXTRAPOLATION ||
//                 SolnBlk.Grid.BCtypeE[j] == BC_COUETTE ||
//                 SolnBlk.Grid.BCtypeE[j] == BC_CONSTANT_EXTRAPOLATION ||
//                 SolnBlk.Grid.BCtypeE[j] == BC_REFLECTION ||
//                 SolnBlk.Grid.BCtypeE[j] == BC_WALL_VISCOUS_ISOTHERMAL ||
//                 SolnBlk.Grid.BCtypeE[j] == BC_TEMPERATURE_SLIP ||
//                 SolnBlk.Grid.BCtypeE[j] == BC_ADIABATIC_WALL) {
//	n_pts = 5;
//	i_index[0] = i;   j_index[0] = j+1;
//	i_index[1] = i-1; j_index[1] = j+1;
//	i_index[2] = i-1; j_index[2] = j;
//	i_index[3] = i-1; j_index[3] = j-1;
//	i_index[4] = i;   j_index[4] = j-1;
//      } else {
//        n_pts = 8;
//	i_index[0] = i-1; j_index[0] = j-1;
//	i_index[1] = i-1; j_index[1] = j;
//	i_index[2] = i-1; j_index[2] = j+1;
//	i_index[3] = i;   j_index[3] = j+1;
//	i_index[4] = i+1; j_index[4] = j+1;
//	i_index[5] = i+1; j_index[5] = j;
//	i_index[6] = i+1; j_index[6] = j-1;
//	i_index[7] = i;   j_index[7] = j-1;
//      } //endif
//
//      //left
//    } else if (i == SolnBlk.ICl-SolnBlk.Nghost+2) {
//
//      if(SolnBlk.Grid.BCtypeW[j] == BC_CHARACTERISTIC ||
//         SolnBlk.Grid.BCtypeW[j] == BC_LINEAR_EXTRAPOLATION ||
//         SolnBlk.Grid.BCtypeW[j] == BC_COUETTE ||
//         SolnBlk.Grid.BCtypeW[j] == BC_CONSTANT_EXTRAPOLATION ||
//         SolnBlk.Grid.BCtypeW[j] == BC_REFLECTION ||
//	 SolnBlk.Grid.BCtypeW[j] == BC_WALL_VISCOUS_ISOTHERMAL ||
//	 SolnBlk.Grid.BCtypeW[j] == BC_TEMPERATURE_SLIP ||
//         SolnBlk.Grid.BCtypeW[j] == BC_ADIABATIC_WALL) {
//	n_pts = 5;
//	i_index[0] = i;   j_index[0] = j+1;
//	i_index[1] = i+1; j_index[1] = j+1;
//	i_index[2] = i+1; j_index[2] = j;
//	i_index[3] = i+1; j_index[3] = j-1;
//	i_index[4] = i;   j_index[4] = j-1;
//      } else {
//        n_pts = 8;
//	i_index[0] = i-1; j_index[0] = j-1;
//	i_index[1] = i-1; j_index[1] = j;
//	i_index[2] = i-1; j_index[2] = j+1;
//	i_index[3] = i;   j_index[3] = j+1;
//	i_index[4] = i+1; j_index[4] = j+1;
//	i_index[5] = i+1; j_index[5] = j;
//	i_index[6] = i+1; j_index[6] = j-1;
//	i_index[7] = i;   j_index[7] = j-1;
//      } //endif
//
//      //top
//    } else if (j == SolnBlk.JCu+SolnBlk.Nghost-2) {
//
//      if(SolnBlk.Grid.BCtypeN[i] == BC_CHARACTERISTIC ||
//         SolnBlk.Grid.BCtypeN[i] == BC_LINEAR_EXTRAPOLATION ||
//         SolnBlk.Grid.BCtypeN[i] == BC_COUETTE ||
//         SolnBlk.Grid.BCtypeN[i] == BC_CONSTANT_EXTRAPOLATION ||
//         SolnBlk.Grid.BCtypeN[i] == BC_REFLECTION ||
//         SolnBlk.Grid.BCtypeN[i] == BC_WALL_VISCOUS_ISOTHERMAL ||
//         SolnBlk.Grid.BCtypeN[i] == BC_TEMPERATURE_SLIP ||
//         SolnBlk.Grid.BCtypeN[i] == BC_ADIABATIC_WALL) {
//	n_pts = 5;
//	i_index[0] = i-1; j_index[0] = j;
//	i_index[1] = i-1; j_index[1] = j-1;
//	i_index[2] = i;   j_index[2] = j-1;
//	i_index[3] = i+1; j_index[3] = j-1;
//	i_index[4] = i+1; j_index[4] = j;
//      } else {
//        n_pts = 8;
//	i_index[0] = i-1; j_index[0] = j-1;
//	i_index[1] = i-1; j_index[1] = j;
//	i_index[2] = i-1; j_index[2] = j+1;
//	i_index[3] = i;   j_index[3] = j+1;
//	i_index[4] = i+1; j_index[4] = j+1;
//	i_index[5] = i+1; j_index[5] = j;
//	i_index[6] = i+1; j_index[6] = j-1;
//	i_index[7] = i;   j_index[7] = j-1;
//      } //endif
//
//      //Right
//    } else if (i == SolnBlk.ICu+SolnBlk.Nghost-2) {
//
//      if(SolnBlk.Grid.BCtypeE[j] == BC_CHARACTERISTIC ||
//         SolnBlk.Grid.BCtypeE[j] == BC_LINEAR_EXTRAPOLATION ||
//         SolnBlk.Grid.BCtypeE[j] == BC_COUETTE ||
//         SolnBlk.Grid.BCtypeE[j] == BC_CONSTANT_EXTRAPOLATION ||
//         SolnBlk.Grid.BCtypeE[j] == BC_REFLECTION ||
//         SolnBlk.Grid.BCtypeE[j] == BC_WALL_VISCOUS_ISOTHERMAL ||
//         SolnBlk.Grid.BCtypeE[j] == BC_TEMPERATURE_SLIP ||
//         SolnBlk.Grid.BCtypeE[j] == BC_ADIABATIC_WALL) {
//	n_pts = 5;
//	i_index[0] = i;   j_index[0] = j+1;
//	i_index[1] = i-1; j_index[1] = j+1;
//	i_index[2] = i-1; j_index[2] = j;
//	i_index[3] = i-1; j_index[3] = j-1;
//	i_index[4] = i;   j_index[4] = j-1;
//      } else {
//        n_pts = 8;
//	i_index[0] = i-1; j_index[0] = j-1;
//	i_index[1] = i-1; j_index[1] = j;
//	i_index[2] = i-1; j_index[2] = j+1;
//	i_index[3] = i;   j_index[3] = j+1;
//	i_index[4] = i+1; j_index[4] = j+1;
//	i_index[5] = i+1; j_index[5] = j;
//	i_index[6] = i+1; j_index[6] = j-1;
//	i_index[7] = i;   j_index[7] = j-1;
//      } //endif
//
//      //Bottom
//    } else if (j == SolnBlk.JCl-SolnBlk.Nghost+2) {
//
//      if(SolnBlk.Grid.BCtypeS[i] == BC_CHARACTERISTIC ||
//         SolnBlk.Grid.BCtypeS[i] == BC_LINEAR_EXTRAPOLATION ||
//         SolnBlk.Grid.BCtypeS[i] == BC_COUETTE ||
//         SolnBlk.Grid.BCtypeS[i] == BC_CONSTANT_EXTRAPOLATION ||
//         SolnBlk.Grid.BCtypeS[i] == BC_REFLECTION ||
//         SolnBlk.Grid.BCtypeS[i] == BC_WALL_VISCOUS_ISOTHERMAL ||
//         SolnBlk.Grid.BCtypeS[i] == BC_TEMPERATURE_SLIP ||
//         SolnBlk.Grid.BCtypeS[i] == BC_ADIABATIC_WALL) {
//	n_pts = 5;
//	i_index[0] = i-1; j_index[0] = j;
//	i_index[1] = i-1; j_index[1] = j+1;
//	i_index[2] = i;   j_index[2] = j+1;
//	i_index[3] = i+1; j_index[3] = j+1;
//	i_index[4] = i+1; j_index[4] = j;
//      } else {
//        n_pts = 8;
//	i_index[0] = i-1; j_index[0] = j-1;
//	i_index[1] = i-1; j_index[1] = j;
//	i_index[2] = i-1; j_index[2] = j+1;
//	i_index[3] = i;   j_index[3] = j+1;
//	i_index[4] = i+1; j_index[4] = j+1;
//	i_index[5] = i+1; j_index[5] = j;
//	i_index[6] = i+1; j_index[6] = j-1;
//	i_index[7] = i;   j_index[7] = j-1;
//      } //endif
//    
//      ///////inside computational domain//////////
//
    } else {
      n_pts = 8;
      i_index[0] = i-1; j_index[0] = j-1;
      i_index[1] = i-1; j_index[1] = j;
      i_index[2] = i-1; j_index[2] = j+1;
      i_index[3] = i;   j_index[3] = j+1;
      i_index[4] = i+1; j_index[4] = j+1;
      i_index[5] = i+1; j_index[5] = j;
      i_index[6] = i+1; j_index[6] = j-1;
      i_index[7] = i;   j_index[7] = j-1;
    } //endif

    if (n_pts > 0) {
        DUDx_ave = Gaussian2D_W_VACUUM;
        DUDy_ave = Gaussian2D_W_VACUUM;
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
           for ( n = 1 ; n <= NUM_VAR_GAUSSIAN2D ; ++n ) {
              u0Min = SolnBlk.W[i][j][n];
              u0Max = u0Min;
              for ( n2 = 0 ; n2 <= n_pts-1 ; ++n2 ) {
                 u0Min = min(u0Min, SolnBlk.W[ i_index[n2] ][ j_index[n2] ][n]);
                 u0Max = max(u0Max, SolnBlk.W[ i_index[n2] ][ j_index[n2] ][n]);
              } /* endfor */

	      ii=0;
	      if(i != SolnBlk.ICu || (SolnBlk.Grid.BCtypeE[j] != BC_ADIABATIC_WALL &&
				      SolnBlk.Grid.BCtypeE[j] != BC_WALL_VISCOUS_ISOTHERMAL &&
				      SolnBlk.Grid.BCtypeE[j] != BC_TEMPERATURE_SLIP)) {
		dX = SolnBlk.Grid.xfaceE(i, j)-SolnBlk.Grid.Cell[i][j].Xc;
		uQuad[ii] = SolnBlk.W[i][j][n] +
			    SolnBlk.dWdx[i][j][n]*dX.x +
			    SolnBlk.dWdy[i][j][n]*dX.y;
		ii++;
	      }
	      if(i != SolnBlk.ICl || (SolnBlk.Grid.BCtypeW[j] != BC_ADIABATIC_WALL &&
				      SolnBlk.Grid.BCtypeW[j] != BC_WALL_VISCOUS_ISOTHERMAL &&
				      SolnBlk.Grid.BCtypeW[j] != BC_TEMPERATURE_SLIP)) {
		dX = SolnBlk.Grid.xfaceW(i, j)-SolnBlk.Grid.Cell[i][j].Xc;
		uQuad[ii] = SolnBlk.W[i][j][n] +
                            SolnBlk.dWdx[i][j][n]*dX.x +
                            SolnBlk.dWdy[i][j][n]*dX.y;
		ii++;
	      }
	      if(j != SolnBlk.JCu || (SolnBlk.Grid.BCtypeN[i] != BC_ADIABATIC_WALL &&
				      SolnBlk.Grid.BCtypeN[i] != BC_WALL_VISCOUS_ISOTHERMAL &&
				      SolnBlk.Grid.BCtypeN[i] != BC_TEMPERATURE_SLIP)) {
		dX = SolnBlk.Grid.xfaceN(i, j)-SolnBlk.Grid.Cell[i][j].Xc;
		uQuad[ii] = SolnBlk.W[i][j][n] +
                           SolnBlk.dWdx[i][j][n]*dX.x +
		           SolnBlk.dWdy[i][j][n]*dX.y;
		ii++;
	      }
	      if(j != SolnBlk.JCl || (SolnBlk.Grid.BCtypeS[i] != BC_ADIABATIC_WALL &&
				      SolnBlk.Grid.BCtypeS[i] != BC_WALL_VISCOUS_ISOTHERMAL &&
				      SolnBlk.Grid.BCtypeS[i] != BC_TEMPERATURE_SLIP)) {
		dX = SolnBlk.Grid.xfaceS(i, j)-SolnBlk.Grid.Cell[i][j].Xc;
		uQuad[ii] = SolnBlk.W[i][j][n] +
                           SolnBlk.dWdx[i][j][n]*dX.x +
                           SolnBlk.dWdy[i][j][n]*dX.y;
		ii++;
	      }

	      switch(Limiter) {
	      case LIMITER_ONE :
		phi = ONE;
		break;
	      case LIMITER_ZERO :
		phi = ZERO;
		break;
	      case LIMITER_BARTH_JESPERSEN :
		phi = Limiter_BarthJespersen(uQuad, SolnBlk.W[i][j][n],
					     u0Min, u0Max, ii);
		break;
	      case LIMITER_VENKATAKRISHNAN :
		phi = Limiter_Venkatakrishnan(uQuad, SolnBlk.W[i][j][n],
					      u0Min, u0Max, ii);
		break;
	      case LIMITER_VANLEER :
		phi = Limiter_VanLeer(uQuad, SolnBlk.W[i][j][n],
				      u0Min, u0Max, ii);
		break;
	      case LIMITER_VANALBADA :
		phi = Limiter_VanAlbada(uQuad, SolnBlk.W[i][j][n],
					u0Min, u0Max, ii);
		break;
	      default:
		phi = Limiter_BarthJespersen(uQuad, SolnBlk.W[i][j][n],
					     u0Min, u0Max, ii);
		break;
	      } /* endswitch */
	      SolnBlk.phi[i][j][n] = phi;
           } /* endfor */

        } /* endif */

	//Check for validity at Quadrature points
	/*	
	dX = SolnBlk.Grid.xfaceE(i, j)-SolnBlk.Grid.Cell[i][j].Xc;
	W_temp = SolnBlk.W[i][j] + 
	         (SolnBlk.phi[i][j]^SolnBlk.dWdx[i][j])*dX.x +
	         (SolnBlk.phi[i][j]^SolnBlk.dWdy[i][j])*dX.y;
	if(W_temp.invalid()) {
	  if(W_temp.DetP() < SolnBlk.W[i][j].DetP() && SolnBlk.W[i][j].DetP() < 0.0){
	    SolnBlk.phi[i][j][4] = 0.0;
	    SolnBlk.phi[i][j][5] = 0.0;
	    SolnBlk.phi[i][j][6] = 0.0;
	    //cout << "<" << W_temp.DetP() << "," << SolnBlk.W[i][j].DetP() << ">";
	  } else if (SolnBlk.W[i][j].DetP() > 0.0) {
	    a = (W_temp.p.xx-SolnBlk.W[i][j].p.xx) *
  	        (W_temp.p.yy-SolnBlk.W[i][j].p.yy) -
	        sqr(W_temp.p.xy-SolnBlk.W[i][j].p.xy);
	    if(a<TOLER) {
	      SolnBlk.phi[i][j][4] = 0.0;
	      SolnBlk.phi[i][j][5] = 0.0;
	      SolnBlk.phi[i][j][6] = 0.0;
	      return;
	    }
	    b = (W_temp.p.xx-SolnBlk.W[i][j].p.xx)*SolnBlk.W[i][j].p.yy +
	        (W_temp.p.yy-SolnBlk.W[i][j].p.yy)*SolnBlk.W[i][j].p.xx -
	        2.0*(W_temp.p.xy-SolnBlk.W[i][j].p.xy)*SolnBlk.W[i][j].p.xy;
	    c = SolnBlk.W[i][j].p.xx*SolnBlk.W[i][j].p.yy-sqr(SolnBlk.W[i][j].p.xy);
	    //cout << "a[" << a << "," << b << "," << c <<"]";
	    //cout.flush();
	    d = (-b+sqrt(sqr(b)-4.0*a*c))/(2.0*a);
	    e = (-b-sqrt(sqr(b)-4.0*a*c))/(2.0*a);
	    if(d >= 0.0 && d <= 1.0) {
	      //cout << "d(" << d << ")" << endl;
	      SolnBlk.phi[i][j][4] = SolnBlk.phi[i][j][4]*d;
	      SolnBlk.phi[i][j][5] = SolnBlk.phi[i][j][5]*d;
	      SolnBlk.phi[i][j][6] = SolnBlk.phi[i][j][6]*d;
	    } else if(e >= 0.0 && e <= 1.0) {
	      //cout << "e(" << e << ")" << endl;
	      SolnBlk.phi[i][j][4] = SolnBlk.phi[i][j][4]*e;
	      SolnBlk.phi[i][j][5] = SolnBlk.phi[i][j][5]*e;
	      SolnBlk.phi[i][j][6] = SolnBlk.phi[i][j][6]*e;
	    } else {
	      cout << "*******(" << d << "," << e << ")*********" << endl;
	      cout << "a[" << a << "," << b << "," << c <<"]";
	      cout.flush();
	    }
	  }
	}

	dX = SolnBlk.Grid.xfaceW(i, j)-SolnBlk.Grid.Cell[i][j].Xc;
	W_temp = SolnBlk.W[i][j] + 
	         (SolnBlk.phi[i][j]^SolnBlk.dWdx[i][j])*dX.x +
	         (SolnBlk.phi[i][j]^SolnBlk.dWdy[i][j])*dX.y;
	if(W_temp.invalid()) {
	  if(W_temp.DetP() < SolnBlk.W[i][j].DetP() && SolnBlk.W[i][j].DetP() < 0.0){
	    SolnBlk.phi[i][j][4] = 0.0;
	    SolnBlk.phi[i][j][5] = 0.0;
	    SolnBlk.phi[i][j][6] = 0.0;
	    //cout << "<" << W_temp.DetP() << "," << SolnBlk.W[i][j].DetP() << ">";
	  } else if (SolnBlk.W[i][j].DetP() > 0.0) {
	    //cout << SolnBlk.W[i][j] << endl << W_temp << endl << SolnBlk.phi[i][j] << endl;
	    //cout.flush();
	    a = (W_temp.p.xx-SolnBlk.W[i][j].p.xx) *
  	        (W_temp.p.yy-SolnBlk.W[i][j].p.yy) -
	        sqr(W_temp.p.xy-SolnBlk.W[i][j].p.xy);
	    //cout << "a=" << a << endl;
	    //cout.flush();
	    if(a<TOLER) {
	      SolnBlk.phi[i][j][4] = 0.0;
	      SolnBlk.phi[i][j][5] = 0.0;
	      SolnBlk.phi[i][j][6] = 0.0;
	      return;
	    }
	    b = (W_temp.p.xx-SolnBlk.W[i][j].p.xx)*SolnBlk.W[i][j].p.yy +
	        (W_temp.p.yy-SolnBlk.W[i][j].p.yy)*SolnBlk.W[i][j].p.xx -
	        2.0*(W_temp.p.xy-SolnBlk.W[i][j].p.xy)*SolnBlk.W[i][j].p.xy;
	    //cout << "b=" << b << endl;
	    //cout.flush();
	    c = SolnBlk.W[i][j].p.xx*SolnBlk.W[i][j].p.yy-sqr(SolnBlk.W[i][j].p.xy);
	    //cout << "b[" << a << "," << b << "," << c <<"]";
	    //cout.flush();
	    d = (-b+sqrt(sqr(b)-4.0*a*c))/(2.0*a);
	    e = (-b-sqrt(sqr(b)-4.0*a*c))/(2.0*a);
	    if(d >= 0.0 && d <= 1.0) {
	      //cout << "d(" << d << ")" << endl;
	      SolnBlk.phi[i][j][4] = SolnBlk.phi[i][j][4]*d;
	      SolnBlk.phi[i][j][5] = SolnBlk.phi[i][j][5]*d;
	      SolnBlk.phi[i][j][6] = SolnBlk.phi[i][j][6]*d;
	    } else if(e >= 0.0 && e <= 1.0) {
	      //cout << "e(" << e << ")" << endl;
	      SolnBlk.phi[i][j][4] = SolnBlk.phi[i][j][4]*e;
	      SolnBlk.phi[i][j][5] = SolnBlk.phi[i][j][5]*e;
	      SolnBlk.phi[i][j][6] = SolnBlk.phi[i][j][6]*e;
	    } else {
	      cout << "*******(" << d << "," << e << ")*********" << endl;
	      cout << "b[" << a << "," << b << "," << c <<"]";
	      cout.flush();
	    }
	  }
	}

	dX = SolnBlk.Grid.xfaceN(i, j)-SolnBlk.Grid.Cell[i][j].Xc;
	W_temp = SolnBlk.W[i][j] + 
	         (SolnBlk.phi[i][j]^SolnBlk.dWdx[i][j])*dX.x +
	         (SolnBlk.phi[i][j]^SolnBlk.dWdy[i][j])*dX.y;
	if(W_temp.invalid()) {
	  if(W_temp.DetP() < SolnBlk.W[i][j].DetP() && SolnBlk.W[i][j].DetP() < 0.0){
	    SolnBlk.phi[i][j][4] = 0.0;
	    SolnBlk.phi[i][j][5] = 0.0;
	    SolnBlk.phi[i][j][6] = 0.0;
	    //cout << "<" << W_temp.DetP() << "," << SolnBlk.W[i][j].DetP() << ">";
	  } else if (SolnBlk.W[i][j].DetP() > 0.0) {
	    a = (W_temp.p.xx-SolnBlk.W[i][j].p.xx) *
  	        (W_temp.p.yy-SolnBlk.W[i][j].p.yy) -
	        sqr(W_temp.p.xy-SolnBlk.W[i][j].p.xy);
	    if(a<TOLER) {
	      SolnBlk.phi[i][j][4] = 0.0;
	      SolnBlk.phi[i][j][5] = 0.0;
	      SolnBlk.phi[i][j][6] = 0.0;
	      return;
	    }
	    b = (W_temp.p.xx-SolnBlk.W[i][j].p.xx)*SolnBlk.W[i][j].p.yy +
	        (W_temp.p.yy-SolnBlk.W[i][j].p.yy)*SolnBlk.W[i][j].p.xx -
	        2.0*(W_temp.p.xy-SolnBlk.W[i][j].p.xy)*SolnBlk.W[i][j].p.xy;
	    c = SolnBlk.W[i][j].p.xx*SolnBlk.W[i][j].p.yy-sqr(SolnBlk.W[i][j].p.xy);
	    //cout << "c[" << a << "," << b << "," << c <<"]";
	    //cout.flush();
	    d = (-b+sqrt(sqr(b)-4.0*a*c))/(2.0*a);
	    e = (-b-sqrt(sqr(b)-4.0*a*c))/(2.0*a);
	    if(d >= 0.0 && d <= 1.0) {
	      //cout << "d(" << d << ")" << endl;
	      SolnBlk.phi[i][j][4] = SolnBlk.phi[i][j][4]*d;
	      SolnBlk.phi[i][j][5] = SolnBlk.phi[i][j][5]*d;
	      SolnBlk.phi[i][j][6] = SolnBlk.phi[i][j][6]*d;
	    } else if(e >= 0.0 && e <= 1.0) {
	      //cout << "e(" << e << ")" << endl;
	      SolnBlk.phi[i][j][4] = SolnBlk.phi[i][j][4]*e;
	      SolnBlk.phi[i][j][5] = SolnBlk.phi[i][j][5]*e;
	      SolnBlk.phi[i][j][6] = SolnBlk.phi[i][j][6]*e;
	    } else {
	      cout << "*******(" << d << "," << e << ")*********" << endl;
	      cout << "c[" << a << "," << b << "," << c <<"]";
	      cout.flush();
	    }
	  }
	}


	dX = SolnBlk.Grid.xfaceS(i, j)-SolnBlk.Grid.Cell[i][j].Xc;
	W_temp = SolnBlk.W[i][j] + 
	         (SolnBlk.phi[i][j]^SolnBlk.dWdx[i][j])*dX.x +
	         (SolnBlk.phi[i][j]^SolnBlk.dWdy[i][j])*dX.y;
	if(W_temp.invalid()) {
	  if(W_temp.DetP() < SolnBlk.W[i][j].DetP() && SolnBlk.W[i][j].DetP() < 0.0){
	    SolnBlk.phi[i][j][4] = 0.0;
	    SolnBlk.phi[i][j][5] = 0.0;
	    SolnBlk.phi[i][j][6] = 0.0;
	    //cout << "<" << W_temp.DetP() << "," << SolnBlk.W[i][j].DetP() << ">";
	  } else if (SolnBlk.W[i][j].DetP() > 0.0) {
	    a = (W_temp.p.xx-SolnBlk.W[i][j].p.xx) *
  	        (W_temp.p.yy-SolnBlk.W[i][j].p.yy) -
	        sqr(W_temp.p.xy-SolnBlk.W[i][j].p.xy);
	    if(a<TOLER) {
	      SolnBlk.phi[i][j][4] = 0.0;
	      SolnBlk.phi[i][j][5] = 0.0;
	      SolnBlk.phi[i][j][6] = 0.0;
	      return;
	    }
	    b = (W_temp.p.xx-SolnBlk.W[i][j].p.xx)*SolnBlk.W[i][j].p.yy +
	        (W_temp.p.yy-SolnBlk.W[i][j].p.yy)*SolnBlk.W[i][j].p.xx -
	        2.0*(W_temp.p.xy-SolnBlk.W[i][j].p.xy)*SolnBlk.W[i][j].p.xy;
	    c = SolnBlk.W[i][j].p.xx*SolnBlk.W[i][j].p.yy-sqr(SolnBlk.W[i][j].p.xy);
	    //cout << "d[" << a << "," << b << "," << c <<"]";
	    //cout.flush();
	    d = (-b+sqrt(sqr(b)-4.0*a*c))/(2.0*a);
	    e = (-b-sqrt(sqr(b)-4.0*a*c))/(2.0*a);
	    if(d >= 0.0 && d <= 1.0) {
	      //cout << "d(" << d << ")" << endl;
	      SolnBlk.phi[i][j][4] = SolnBlk.phi[i][j][4]*d;
	      SolnBlk.phi[i][j][5] = SolnBlk.phi[i][j][5]*d;
	      SolnBlk.phi[i][j][6] = SolnBlk.phi[i][j][6]*d;
	    } else if(e >= 0.0 && e <= 1.0) {
	      //cout << "e(" << e << ")" << endl;
	      SolnBlk.phi[i][j][4] = SolnBlk.phi[i][j][4]*e;
	      SolnBlk.phi[i][j][5] = SolnBlk.phi[i][j][5]*e;
	      SolnBlk.phi[i][j][6] = SolnBlk.phi[i][j][6]*e;
	    } else {
	      cout << "*******(" << d << "," << e << ")*********" << endl;
	      cout << "d[" << a << "," << b << "," << c <<"]";
	      cout.flush();
	    }
	  }
	}
	*/
	/*
	dX = SolnBlk.Grid.xfaceE(i, j)-SolnBlk.Grid.Cell[i][j].Xc;
	while((SolnBlk.W[i][j] + 
	    (SolnBlk.phi[i][j]^SolnBlk.dWdx[i][j])*dX.x +
	    (SolnBlk.phi[i][j]^SolnBlk.dWdy[i][j])*dX.y).invalid()) {
	  reduction_counter++;
	  SolnBlk.phi[i][j] = SolnBlk.phi[i][j]*0.99;//Gaussian2D_W_VACUUM;
	}
	dX = SolnBlk.Grid.xfaceW(i, j)-SolnBlk.Grid.Cell[i][j].Xc;
	while((SolnBlk.W[i][j] + 
	    (SolnBlk.phi[i][j]^SolnBlk.dWdx[i][j])*dX.x +
	    (SolnBlk.phi[i][j]^SolnBlk.dWdy[i][j])*dX.y).invalid()) {
	  SolnBlk.phi[i][j] = SolnBlk.phi[i][j]*0.99;//Gaussian2D_W_VACUUM;
	  reduction_counter++;
	}
	dX = SolnBlk.Grid.xfaceN(i, j)-SolnBlk.Grid.Cell[i][j].Xc;
	while((SolnBlk.W[i][j] + 
	    (SolnBlk.phi[i][j]^SolnBlk.dWdx[i][j])*dX.x +
	    (SolnBlk.phi[i][j]^SolnBlk.dWdy[i][j])*dX.y).invalid()) {
	  SolnBlk.phi[i][j] = SolnBlk.phi[i][j]*0.99;//Gaussian2D_W_VACUUM;
	  reduction_counter++;
	}
	dX = SolnBlk.Grid.xfaceS(i, j)-SolnBlk.Grid.Cell[i][j].Xc;
	while((SolnBlk.W[i][j] + 
	    (SolnBlk.phi[i][j]^SolnBlk.dWdx[i][j])*dX.x +
	    (SolnBlk.phi[i][j]^SolnBlk.dWdy[i][j])*dX.y).invalid()) {
	  SolnBlk.phi[i][j] = SolnBlk.phi[i][j]*0.99;//Gaussian2D_W_VACUUM;
	  reduction_counter++;
	}
	*/
    } else {
        SolnBlk.dWdx[i][j] = Gaussian2D_W_VACUUM;
        SolnBlk.dWdy[i][j] = Gaussian2D_W_VACUUM; 
        SolnBlk.phi[i][j]  = Gaussian2D_W_VACUUM;
    } /* endif */

    if(reduction_counter > 0) {
      cout << "^";
      reduction_counter = 0;
    }

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
void Linear_Reconstruction_LeastSquares(Gaussian2D_Quad_Block &SolnBlk,
				        const int Limiter) {

    int i, j;

    /* Carry out the limited solution reconstruction in
       each cell of the computational mesh. */

    for ( j  = SolnBlk.JCl-SolnBlk.Nghost+1 ; j <= SolnBlk.JCu+SolnBlk.Nghost-1 ; ++j ) {
       for ( i = SolnBlk.ICl-SolnBlk.Nghost+1 ; i <= SolnBlk.ICu+SolnBlk.Nghost-1 ; ++i ) {
	   Linear_Reconstruction_LeastSquares(SolnBlk, i, j, Limiter);
       } /* endfor */
    } /* endfor */
    //cout << "DONE" << endl;cout.flush();
}

/********************************************************
 * Routine: Calculate_Refinement_Criteria               *
 *                                                      *
 * Calculate refinement criteria for the solution       *
 * block.                                               *
 *                                                      *
 ********************************************************/
void Calculate_Refinement_Criteria(double *refinement_criteria,
				   Gaussian2D_Input_Parameters &IP,
                                   int &number_refinement_criteria,
                                   Gaussian2D_Quad_Block &SolnBlk) {

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
	      div_V_criteria = sqrt(SolnBlk.Grid.Cell[i][j].A)*fabs(div_V)/SolnBlk.W[i][j].sound();
	    } else {
	      div_V_criteria = ONE;
	    }
	    div_V_criteria_max = max(div_V_criteria_max, div_V_criteria);

             // Evaluate refinement criteria #3 based on the curl
             // of the velocity vector.
	    if (IP.Refinement_Criteria_Curl_Velocity) {
	      curl_V_z = SolnBlk.dWdx[i][j].v.y - SolnBlk.dWdy[i][j].v.x; 
	      curl_V_abs = sqrt(sqr(curl_V_z)); 
	      curl_V_criteria = sqrt(SolnBlk.Grid.Cell[i][j].A)*curl_V_abs/SolnBlk.W[i][j].sound();
	      //curl_V_criteria = curl_V_abs;
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
void Fix_Refined_Block_Boundaries(Gaussian2D_Quad_Block &SolnBlk,
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
void Unfix_Refined_Block_Boundaries(Gaussian2D_Quad_Block &SolnBlk) {

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
void Apply_Boundary_Flux_Corrections(Gaussian2D_Quad_Block &SolnBlk,
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
void Apply_Boundary_Flux_Corrections_Multistage_Explicit(Gaussian2D_Quad_Block &SolnBlk,
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
 *  solution block using a 2nd-order limited upwind     *
 * finite-volume spatial discretization scheme with     *
 * either the Godunov, Roe, Rusanov, HLLE, Linde, or    *
 * HLLC flux functions.                                 *
 * The residual is stored in dUdt[][][0].               *
 *                                                      *
 ********************************************************/
void dUdt_Residual_Evaluation(Gaussian2D_Quad_Block &SolnBlk,
                              Gaussian2D_Input_Parameters &Input_Parameters) {

    int i, j;
    Vector2D dX;
    Gaussian2D_pState Wl, Wr;
    Gaussian2D_cState Flux;

    Gaussian2D_pState Wu, Wd, dWdxl, dWdyl, dWdxr, dWdyr;
    Vector2D Xl, Xr, Xu, Xd;
    int elliptic_bc_flag;

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
      SolnBlk.dUdt[SolnBlk.ICl-1][j][0] = Gaussian2D_U_VACUUM;
          
      for ( i = SolnBlk.ICl-1 ; i <= SolnBlk.ICu ; ++i ) {

	SolnBlk.dUdt[i+1][j][0] = Gaussian2D_U_VACUUM;
    
	if ( j > SolnBlk.JCl-1 && j < SolnBlk.JCu+1 ) {
    
	  /* Evaluate the cell interface i-direction fluxes. */
	  
	  if (i == SolnBlk.ICl-1 && 
	      (SolnBlk.Grid.BCtypeW[j] == BC_REFLECTION ||
	       SolnBlk.Grid.BCtypeW[j] == BC_CHARACTERISTIC ||
	       SolnBlk.Grid.BCtypeW[j] == BC_CHARACTERISTIC_VELOCITY ||
	       SolnBlk.Grid.BCtypeW[j] == BC_WALL_VISCOUS_ISOTHERMAL ||
	       SolnBlk.Grid.BCtypeW[j] == BC_TEMPERATURE_SLIP ||
	       SolnBlk.Grid.BCtypeW[j] == BC_ADIABATIC_WALL)) {
	    dX = SolnBlk.Grid.xfaceW(i+1, j)-SolnBlk.Grid.Cell[i+1][j].Xc;
	    Wr = SolnBlk.W[i+1][j] + 
	      (SolnBlk.phi[i+1][j]^SolnBlk.dWdx[i+1][j])*dX.x +
	      (SolnBlk.phi[i+1][j]^SolnBlk.dWdy[i+1][j])*dX.y;
	    if (SolnBlk.Grid.BCtypeW[j] == BC_REFLECTION) {
	      Wl = Reflect(Wr, SolnBlk.Grid.nfaceW(i+1, j));
	    } else if(SolnBlk.Grid.BCtypeW[j] == BC_ADIABATIC_WALL) {
	      Wl = Adiabatic_Wall(Wr, SolnBlk.WoW[j].v,
                                  SolnBlk.Grid.nfaceW(i+1, j));
	    } else if(SolnBlk.Grid.BCtypeW[j] == BC_WALL_VISCOUS_ISOTHERMAL) {
	      Wl = Isothermal_Wall(Wr, SolnBlk.WoW[j].v,SolnBlk.WoW[j].T(), SolnBlk.Grid.nfaceW(i+1, j));
	    } else if(SolnBlk.Grid.BCtypeW[j] == BC_TEMPERATURE_SLIP) {
	      Wl = Isothermal_Wall_Slip_T(Wr, SolnBlk.WoW[j].v,SolnBlk.WoW[j].T(),
					  SolnBlk.oldT_W[j],
					  SolnBlk.phi[i+1][j]^SolnBlk.dWdx[i+1][j],
					  SolnBlk.phi[i+1][j]^SolnBlk.dWdy[i+1][j],
					  SolnBlk.Grid.nfaceW(i+1, j));
	    } else if(SolnBlk.Grid.BCtypeW[j] == BC_CHARACTERISTIC_VELOCITY) {
	      Wl = BC_Characteristic_Velocity(Wr, 
					      SolnBlk.WoW[j], 
					      SolnBlk.Grid.nfaceW(i+1, j));
	    } else {
	      Wl = BC_Characteristic_Pressure(Wr, 
					      SolnBlk.WoW[j], 
					      SolnBlk.Grid.nfaceW(i+1, j));
	    } /* endif */
	  } else if (i == SolnBlk.ICu && 
		     (SolnBlk.Grid.BCtypeE[j] == BC_REFLECTION ||
		      SolnBlk.Grid.BCtypeE[j] == BC_CHARACTERISTIC ||
		      SolnBlk.Grid.BCtypeE[j] == BC_CHARACTERISTIC_VELOCITY ||
		      SolnBlk.Grid.BCtypeE[j] == BC_WALL_VISCOUS_ISOTHERMAL ||
		      SolnBlk.Grid.BCtypeE[j] == BC_TEMPERATURE_SLIP ||
		      SolnBlk.Grid.BCtypeE[j] == BC_ADIABATIC_WALL)) {
	    dX = SolnBlk.Grid.xfaceE(i, j)-SolnBlk.Grid.Cell[i][j].Xc;
	    Wl = SolnBlk.W[i][j] + 
	      (SolnBlk.phi[i][j]^SolnBlk.dWdx[i][j])*dX.x +
	      (SolnBlk.phi[i][j]^SolnBlk.dWdy[i][j])*dX.y;
	    if (SolnBlk.Grid.BCtypeE[j] == BC_REFLECTION) {
	      Wr = Reflect(Wl, SolnBlk.Grid.nfaceE(i, j));
	    } else if (SolnBlk.Grid.BCtypeE[j] == BC_ADIABATIC_WALL) {
	      Wr = Adiabatic_Wall(Wl, SolnBlk.WoE[j].v, 
	      			     SolnBlk.Grid.nfaceE(i, j));
	    } else if (SolnBlk.Grid.BCtypeE[j] == BC_WALL_VISCOUS_ISOTHERMAL) {
	      Wr = Isothermal_Wall(Wl, SolnBlk.WoE[j].v, SolnBlk.WoE[j].T(), SolnBlk.Grid.nfaceE(i, j)); 
	    } else if (SolnBlk.Grid.BCtypeE[j] == BC_TEMPERATURE_SLIP) {
	      Wr = Isothermal_Wall_Slip_T(Wl, SolnBlk.WoE[j].v, SolnBlk.WoE[j].T(),
					  SolnBlk.oldT_E[j],
					  SolnBlk.phi[i][j]^SolnBlk.dWdx[i][j],
					  SolnBlk.phi[i][j]^SolnBlk.dWdy[i][j],
					  SolnBlk.Grid.nfaceE(i, j));
	    } else if (SolnBlk.Grid.BCtypeE[j] == BC_CHARACTERISTIC_VELOCITY) {
	      Wr = BC_Characteristic_Velocity(Wl, 
					      SolnBlk.WoE[j], 
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

	  switch(Input_Parameters.i_Flux_Function) {
	    //case FLUX_FUNCTION_GODUNOV :
	    //Flux = FluxGodunov_n(Wl, Wr, SolnBlk.Grid.nfaceE(i, j));
	    //break;
	  case FLUX_FUNCTION_ROE :
	    Flux = FluxRoe_n(Wl, Wr, SolnBlk.Grid.nfaceE(i, j));
	    break;
	    //case FLUX_FUNCTION_RUSANOV :
	    //Flux = FluxRusanov_n(Wl, Wr, SolnBlk.Grid.nfaceE(i, j));
	    //break;
	  case FLUX_FUNCTION_HLLE :
	    Flux = FluxHLLE_n(Wl, Wr, SolnBlk.Grid.nfaceE(i, j));
	    break;
	    //case FLUX_FUNCTION_LINDE :
	    //Flux = FluxLinde_n(Wl, Wr, SolnBlk.Grid.nfaceE(i, j));
	    //break;
	    //case FLUX_FUNCTION_HLLC :
	    //Flux = FluxHLLC_n(Wl, Wr, SolnBlk.Grid.nfaceE(i, j));
	    //break;
	  default:
	    Flux = FluxRoe_n(Wl, Wr, SolnBlk.Grid.nfaceE(i, j));
	    break;
	  } /* endswitch */

#ifdef _GAUSSIAN_HEAT_TRANSFER_
	  // Evaluate the cell interface i-direction ELLIPTIC flux if necessary.
	  if (SolnBlk.Heat_Transfer) {
	    // Determine the EAST face ELLIPTIC flux.
	    if (i == SolnBlk.ICl-1 && 
		(SolnBlk.Grid.BCtypeW[j] == BC_ADIABATIC_WALL ||
		 SolnBlk.Grid.BCtypeW[j] == BC_WALL_VISCOUS_ISOTHERMAL ||
		 SolnBlk.Grid.BCtypeW[j] == BC_TEMPERATURE_SLIP)) {
	      // WEST face of cell (i+1,j) is a normal boundary.
	      Xr = SolnBlk.Grid.Cell[i+1][j].Xc; Wr = SolnBlk.W[i+1][j];
	      if (SolnBlk.Grid.BCtypeW[j] == BC_ADIABATIC_WALL) {
		// WEST face of cell (i+1,j) is an ADIABATIC_WALL boundary.
		elliptic_bc_flag = DIAMONDPATH_RIGHT_TRIANGLE_HEATFLUX;
		Wu = Knudsen_Layer_Adiabatic(Wr, SolnBlk.WoW[j].v,
					     SolnBlk.Grid.nfaceW(i+1, j));
	      }else if (SolnBlk.Grid.BCtypeW[j] == BC_WALL_VISCOUS_ISOTHERMAL) {
		// WEST face of cell (i+1,j) is an ISOTHERMAL_WALL boundary.
		elliptic_bc_flag = DIAMONDPATH_RIGHT_TRIANGLE_HEATFLUX;
		Wu = Knudsen_Layer_Isothermal(Wr, SolnBlk.WoW[j].v, SolnBlk.WoW[j].T(),
					      SolnBlk.Grid.nfaceW(i+1, j));
	      }else if (SolnBlk.Grid.BCtypeW[j] == BC_TEMPERATURE_SLIP) {
		// WEST face of cell (i+1,j) is a TEMPERATURE_SLIP boundary.
		elliptic_bc_flag = DIAMONDPATH_RIGHT_TRIANGLE_HEATFLUX;
		Wu = Knudsen_Layer_Isothermal_Slip_T(Wr, SolnBlk.WoW[j].v, SolnBlk.WoW[j].T(),
						     SolnBlk.oldT_W[j],
						     SolnBlk.phi[i+1][j]^SolnBlk.dWdx[i+1][j],
						     SolnBlk.phi[i+1][j]^SolnBlk.dWdy[i+1][j],
						     SolnBlk.Grid.nfaceW(i+1, j));
	      } else {cout << "Error bad BC for elliptic part" << endl;}
	      switch(Input_Parameters.i_Heat_Reconstruction) {
	      case ELLIPTIC_RECONSTRUCTION_CARTESIAN :
	      case ELLIPTIC_RECONSTRUCTION_DIAMOND_PATH :
		Xu = SolnBlk.Grid.Node[i+1][j+1].X;
		Xd = SolnBlk.Grid.Node[i+1][j  ].X; Wd = Wu;
		break;
	      case ELLIPTIC_RECONSTRUCTION_ARITHMETIC_AVERAGE :
	      case ELLIPTIC_RECONSTRUCTION_HYBRID :
		Xu = SolnBlk.Grid.xfaceW(i+1,j);
		Xl = Xu; Wl = Wu;
		dWdxr = SolnBlk.dWdx[i+1][j]; dWdxl = dWdxr;
		dWdyr = SolnBlk.dWdy[i+1][j]; dWdyl = dWdyr;
		break;
	      };
	    } else if (i == SolnBlk.ICu &&
		       (SolnBlk.Grid.BCtypeE[j] == BC_ADIABATIC_WALL ||
			SolnBlk.Grid.BCtypeE[j] == BC_WALL_VISCOUS_ISOTHERMAL ||
			SolnBlk.Grid.BCtypeE[j] == BC_TEMPERATURE_SLIP)) {
	      // EAST face of cell (i,j) is a normal boundary.
	      Xl = SolnBlk.Grid.Cell[i][j].Xc; Wl = SolnBlk.W[i][j];
	      if (SolnBlk.Grid.BCtypeE[j] == BC_ADIABATIC_WALL) {
		// EAST face of cell (i,j) is an ADIABATIC_WALL boundary.
		elliptic_bc_flag = DIAMONDPATH_LEFT_TRIANGLE_HEATFLUX;
		Wu = Knudsen_Layer_Adiabatic(Wl, SolnBlk.WoE[j].v, 
					     SolnBlk.Grid.nfaceE(i, j));
	      } else if (SolnBlk.Grid.BCtypeE[j] == BC_WALL_VISCOUS_ISOTHERMAL) {
		// EAST face of cell (i,j) is an ISOTHERMAL_WALL boundary.
		elliptic_bc_flag = DIAMONDPATH_LEFT_TRIANGLE_HEATFLUX;
		Wu = Knudsen_Layer_Isothermal(Wl, SolnBlk.WoE[j].v, SolnBlk.WoE[j].T(),
					      SolnBlk.Grid.nfaceE(i, j));
	      } else if (SolnBlk.Grid.BCtypeE[j] == BC_TEMPERATURE_SLIP) {
		// EAST face of cell (i,j) is a TEMPERATURE_SLIP boundary.
		elliptic_bc_flag = DIAMONDPATH_LEFT_TRIANGLE_HEATFLUX;
		Wu = Knudsen_Layer_Isothermal_Slip_T(Wl, SolnBlk.WoE[j].v, SolnBlk.WoE[j].T(),
						     SolnBlk.oldT_E[j],
						     SolnBlk.phi[i][j]^SolnBlk.dWdx[i][j],
						     SolnBlk.phi[i][j]^SolnBlk.dWdy[i][j],
						     SolnBlk.Grid.nfaceE(i, j));
	      } else {cout << "Error bad BC for elliptic part" << endl;}

	      switch(Input_Parameters.i_Heat_Reconstruction) {
	      case ELLIPTIC_RECONSTRUCTION_CARTESIAN :
	      case ELLIPTIC_RECONSTRUCTION_DIAMOND_PATH :
		Xu = SolnBlk.Grid.Node[i+1][j+1].X;
		Xd = SolnBlk.Grid.Node[i+1][j  ].X; Wd = Wu;
		break;
	      case ELLIPTIC_RECONSTRUCTION_ARITHMETIC_AVERAGE :
	      case ELLIPTIC_RECONSTRUCTION_HYBRID :
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
	      switch(Input_Parameters.i_Heat_Reconstruction) {
	      case ELLIPTIC_RECONSTRUCTION_CARTESIAN :
	      case ELLIPTIC_RECONSTRUCTION_DIAMOND_PATH :
		elliptic_bc_flag = DIAMONDPATH_QUADRILATERAL_RECONSTRUCTION;
		Xu = SolnBlk.Grid.Node[i+1][j+1].X; Wu = SolnBlk.WnNE(i,j);
		Xd = SolnBlk.Grid.Node[i+1][j  ].X; Wd = SolnBlk.WnSE(i,j);
		break;
	      case ELLIPTIC_RECONSTRUCTION_ARITHMETIC_AVERAGE :
	      case ELLIPTIC_RECONSTRUCTION_HYBRID :
		Xu = SolnBlk.Grid.xfaceE(i,j); Wu = HALF*(Wl + Wr);
		dWdxl = SolnBlk.dWdx[i][j]; dWdxr = SolnBlk.dWdy[i+1][j];
		dWdyl = SolnBlk.dWdy[i][j]; dWdyr = SolnBlk.dWdy[i+1][j];
		break;
	      };
	    }
	    // Compute the EAST face elliptic flux.
	    switch(Input_Parameters.i_Heat_Reconstruction) {
	    case ELLIPTIC_RECONSTRUCTION_CARTESIAN :
	    case ELLIPTIC_RECONSTRUCTION_DIAMOND_PATH :
	      Flux -= HeatFluxDiamondPath_n(SolnBlk.Grid.xfaceE(i,j),
					       Xl,Wl,Xd,Wd,Xr,Wr,Xu,Wu,
					       SolnBlk.Grid.nfaceE(i,j),
					       SolnBlk.Axisymmetric,
					       elliptic_bc_flag);
	      break;
	    case ELLIPTIC_RECONSTRUCTION_ARITHMETIC_AVERAGE :
	      Flux -= HeatFlux_n(Xu,Wu,HALF*(dWdxl+dWdxr),HALF*(dWdyl+dWdyr),
				    SolnBlk.Grid.nfaceE(i,j),SolnBlk.Axisymmetric);
	      break;
	    case ELLIPTIC_RECONSTRUCTION_HYBRID :
	      Flux -= HeatFluxHybrid_n(Xu,Wu,Xl,Wl,dWdxl,dWdyl,Xr,Wr,dWdxr,dWdyr,
					  SolnBlk.Grid.nfaceE(i,j),SolnBlk.Axisymmetric);
	      break;
	    };
	  }
#endif
	  /* Evaluate cell-averaged solution changes. */
	  
	  SolnBlk.dUdt[i][j][0] -= 
	    Flux*SolnBlk.Grid.lfaceE(i, j)/
	    SolnBlk.Grid.Cell[i][j].A;
	  SolnBlk.dUdt[i+1][j][0] += 
	    Flux*SolnBlk.Grid.lfaceW(i+1, j)/
	    SolnBlk.Grid.Cell[i+1][j].A;

	  /* Include axisymmetric source terms as required. */

	  if (SolnBlk.Axisymmetric) {
	    SolnBlk.dUdt[i][j][0] += 
	      S(SolnBlk.W[i][j], SolnBlk.Grid.Cell[i][j].Xc);
	  } /* endif */

	  /* Save west and east face boundary flux. */
	  
	  if (i == SolnBlk.ICl-1) {
	    SolnBlk.FluxW[j] = -Flux*SolnBlk.Grid.lfaceW(i+1, j);
	  } else if (i == SolnBlk.ICu) {
	    SolnBlk.FluxE[j] = Flux*SolnBlk.Grid.lfaceE(i, j);
	  } /* endif */ 

	} /* endif */
      } /* endfor */
      
      if ( j > SolnBlk.JCl-1 && j < SolnBlk.JCu+1 ) {
	SolnBlk.dUdt[SolnBlk.ICl-1][j][0] = Gaussian2D_U_VACUUM;
	SolnBlk.dUdt[SolnBlk.ICu+1][j][0] = Gaussian2D_U_VACUUM;
      } /* endif */
    } /* endfor */
    
    // Add j-direction (eta-direction) fluxes.
    for ( i = SolnBlk.ICl ; i <= SolnBlk.ICu ; ++i ) {
      for ( j  = SolnBlk.JCl-1 ; j <= SolnBlk.JCu ; ++j ) {
	
	/* Evaluate the cell interface j-direction fluxes. */
	
	if (j == SolnBlk.JCl-1 && 
	    (SolnBlk.Grid.BCtypeS[i] == BC_REFLECTION ||
	     SolnBlk.Grid.BCtypeS[i] == BC_CHARACTERISTIC ||
	     SolnBlk.Grid.BCtypeS[i] == BC_CHARACTERISTIC_VELOCITY ||
	     SolnBlk.Grid.BCtypeS[i] == BC_WALL_VISCOUS_ISOTHERMAL ||
	     SolnBlk.Grid.BCtypeS[i] == BC_TEMPERATURE_SLIP ||
	     SolnBlk.Grid.BCtypeS[i] == BC_ADIABATIC_WALL )) {
	  dX = SolnBlk.Grid.xfaceS(i, j+1)-SolnBlk.Grid.Cell[i][j+1].Xc;
	  Wr = SolnBlk.W[i][j+1] +
	    (SolnBlk.phi[i][j+1]^SolnBlk.dWdx[i][j+1])*dX.x +
	    (SolnBlk.phi[i][j+1]^SolnBlk.dWdy[i][j+1])*dX.y;
	  if (SolnBlk.Grid.BCtypeS[i] == BC_REFLECTION) {
	    Wl = Reflect(Wr, SolnBlk.Grid.nfaceS(i, j+1));
	  } else if (SolnBlk.Grid.BCtypeS[i] == BC_ADIABATIC_WALL) {
	    Wl = Adiabatic_Wall(Wr, SolnBlk.WoS[i].v, 
	    			       SolnBlk.Grid.nfaceS(i, j+1));
	  } else if (SolnBlk.Grid.BCtypeS[i] == BC_WALL_VISCOUS_ISOTHERMAL) {
	    Wl = Isothermal_Wall(Wr, SolnBlk.WoS[i].v, SolnBlk.WoS[i].T(), SolnBlk.Grid.nfaceS(i, j+1));
	  } else if (SolnBlk.Grid.BCtypeS[i] == BC_TEMPERATURE_SLIP) {
	    Wl = Isothermal_Wall_Slip_T(Wr, SolnBlk.WoS[i].v, SolnBlk.WoS[i].T(), 
					SolnBlk.oldT_S[i],
					SolnBlk.phi[i][j+1]^SolnBlk.dWdx[i][j+1],
					SolnBlk.phi[i][j+1]^SolnBlk.dWdy[i][j+1],
					SolnBlk.Grid.nfaceS(i, j+1));
	  } else if (SolnBlk.Grid.BCtypeS[i] == BC_CHARACTERISTIC_VELOCITY) {
	    Wl = BC_Characteristic_Velocity(Wr, 
					    SolnBlk.WoS[i], 
					    SolnBlk.Grid.nfaceS(i, j+1));
	  } else {
	    Wl = BC_Characteristic_Pressure(Wr, 
					    SolnBlk.WoS[i], 
					    SolnBlk.Grid.nfaceS(i, j+1));
	  } /* endif */
	} else if (j == SolnBlk.JCu && 
		   (SolnBlk.Grid.BCtypeN[i] == BC_REFLECTION ||
		    SolnBlk.Grid.BCtypeN[i] == BC_CHARACTERISTIC ||
		    SolnBlk.Grid.BCtypeN[i] == BC_CHARACTERISTIC_VELOCITY ||
		    SolnBlk.Grid.BCtypeN[i] == BC_WALL_VISCOUS_ISOTHERMAL ||
		    SolnBlk.Grid.BCtypeN[i] == BC_TEMPERATURE_SLIP ||
		    SolnBlk.Grid.BCtypeN[i] == BC_ADIABATIC_WALL )) {
	  dX = SolnBlk.Grid.xfaceN(i, j)-SolnBlk.Grid.Cell[i][j].Xc;
	  Wl = SolnBlk.W[i][j] + 
	    (SolnBlk.phi[i][j]^SolnBlk.dWdx[i][j])*dX.x +
	    (SolnBlk.phi[i][j]^SolnBlk.dWdy[i][j])*dX.y;
	  if (SolnBlk.Grid.BCtypeN[i] == BC_REFLECTION) {
	    Wr = Reflect(Wl, SolnBlk.Grid.nfaceN(i, j));
	  } else if (SolnBlk.Grid.BCtypeN[i] == BC_ADIABATIC_WALL) {
	    Wr = Adiabatic_Wall(Wl, SolnBlk.WoN[i].v, 
	    			       SolnBlk.Grid.nfaceN(i, j));
	  } else if (SolnBlk.Grid.BCtypeN[i] == BC_WALL_VISCOUS_ISOTHERMAL) {
	    Wr = Isothermal_Wall(Wl, SolnBlk.WoN[i].v, SolnBlk.WoN[i].T(), SolnBlk.Grid.nfaceN(i, j));
	  } else if (SolnBlk.Grid.BCtypeN[i] == BC_TEMPERATURE_SLIP) {
	    Wr = Isothermal_Wall_Slip_T(Wl, SolnBlk.WoN[i].v, SolnBlk.WoN[i].T(),
					SolnBlk.oldT_N[i],
					SolnBlk.phi[i][j]^SolnBlk.dWdx[i][j],
					SolnBlk.phi[i][j]^SolnBlk.dWdy[i][j],
					SolnBlk.Grid.nfaceN(i, j));
	  } else if (SolnBlk.Grid.BCtypeN[i] == BC_CHARACTERISTIC_VELOCITY) {
	    Wr = BC_Characteristic_Velocity(Wl, 
	    				    SolnBlk.WoN[i], 
					    SolnBlk.Grid.nfaceN(i, j));
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
	
	switch(Input_Parameters.i_Flux_Function) {
	  //case FLUX_FUNCTION_GODUNOV :
	  //Flux = FluxGodunov_n(Wl, Wr, SolnBlk.Grid.nfaceN(i, j));
	  //break;
	case FLUX_FUNCTION_ROE :
	  Flux = FluxRoe_n(Wl, Wr, SolnBlk.Grid.nfaceN(i, j));
	  break;
	  //case FLUX_FUNCTION_RUSANOV :
	  //Flux = FluxRusanov_n(Wl, Wr, SolnBlk.Grid.nfaceN(i, j));
	  //break;
	case FLUX_FUNCTION_HLLE :
	  Flux = FluxHLLE_n(Wl, Wr, SolnBlk.Grid.nfaceN(i, j));
	  break;
	  //case FLUX_FUNCTION_LINDE :
	  //Flux = FluxLinde_n(Wl, Wr, SolnBlk.Grid.nfaceN(i, j));
	  //break;
	  //case FLUX_FUNCTION_HLLC :
	  //Flux = FluxHLLC_n(Wl, Wr, SolnBlk.Grid.nfaceN(i, j));
	  //break;
	default:
	  Flux = FluxRoe_n(Wl, Wr, SolnBlk.Grid.nfaceN(i, j));
	  break;
	} /* endswitch */

#ifdef _GAUSSIAN_HEAT_TRANSFER_
	// Evaluate the cell interface j-direction ELLIPTIC flux if necessary.
	if (SolnBlk.Heat_Transfer) {
	  // Determine the NORTH face ELLIPTIC flux.
	  if (j == SolnBlk.JCl-1 && 
	      (SolnBlk.Grid.BCtypeS[i] == BC_ADIABATIC_WALL  ||
	       SolnBlk.Grid.BCtypeS[i] == BC_WALL_VISCOUS_ISOTHERMAL ||
	       SolnBlk.Grid.BCtypeS[i] == BC_TEMPERATURE_SLIP)) {
	    // SOUTH face of cell (i,j+1) is a normal boundary.
	    Xr = SolnBlk.Grid.Cell[i][j+1].Xc; Wr = SolnBlk.W[i][j+1];
	    if (SolnBlk.Grid.BCtypeS[i] == BC_ADIABATIC_WALL) {
	      // SOUTH face of cell (i,j+1) is an ADIABATIC_WALL boundary.
	      elliptic_bc_flag = DIAMONDPATH_RIGHT_TRIANGLE_HEATFLUX;
	      Wu = Knudsen_Layer_Adiabatic(Wr, SolnBlk.WoS[i].v, 
					   SolnBlk.Grid.nfaceS(i, j+1));
	    } else if (SolnBlk.Grid.BCtypeS[i] == BC_WALL_VISCOUS_ISOTHERMAL) {
	      // SOUTH face of cell (i,j+1) is an ISOTHERMAL_WALL boundary.
	      elliptic_bc_flag = DIAMONDPATH_RIGHT_TRIANGLE_HEATFLUX;
	      Wu = Knudsen_Layer_Isothermal(Wr, SolnBlk.WoS[i].v, SolnBlk.WoS[i].T(),
					    SolnBlk.Grid.nfaceS(i, j+1));
	    } else if (SolnBlk.Grid.BCtypeS[i] == BC_TEMPERATURE_SLIP) {
	      // SOUTH face of cell (i,j+1) is a TEMPERATURE_SLIP boundary.
	      elliptic_bc_flag = DIAMONDPATH_RIGHT_TRIANGLE_HEATFLUX;
	      Wu = Knudsen_Layer_Isothermal_Slip_T(Wr, SolnBlk.WoS[i].v, SolnBlk.WoS[i].T(),
						   SolnBlk.oldT_S[i],
						   SolnBlk.phi[i][j+1]^SolnBlk.dWdx[i][j+1],
						   SolnBlk.phi[i][j+1]^SolnBlk.dWdy[i][j+1],
						   SolnBlk.Grid.nfaceS(i, j+1));
	    } else {cout << "Error bad BC for elliptic part" << endl;}
	    switch(Input_Parameters.i_Heat_Reconstruction) {
	    case ELLIPTIC_RECONSTRUCTION_CARTESIAN :
	    case ELLIPTIC_RECONSTRUCTION_DIAMOND_PATH :
	      Xu = SolnBlk.Grid.Node[i  ][j+1].X;
	      Xd = SolnBlk.Grid.Node[i+1][j+1].X; Wd = Wu;
	      break;
	    case ELLIPTIC_RECONSTRUCTION_ARITHMETIC_AVERAGE :
	    case ELLIPTIC_RECONSTRUCTION_HYBRID :
	      Xu = SolnBlk.Grid.xfaceS(i,j+1);
	      Xl = Xu; Wl = Wu;
	      dWdxr = SolnBlk.dWdx[i][j+1]; dWdxl = dWdxr;
	      dWdyr = SolnBlk.dWdy[i][j+1]; dWdyl = dWdyr;
	      break;
	    };
	  } else if (j == SolnBlk.JCu && 
		     (SolnBlk.Grid.BCtypeN[i] == BC_ADIABATIC_WALL ||
		      SolnBlk.Grid.BCtypeN[i] == BC_WALL_VISCOUS_ISOTHERMAL ||
		      SolnBlk.Grid.BCtypeN[i] == BC_TEMPERATURE_SLIP)) {
	    // NORTH face of cell (i,j) is a normal boundary.
	    Xl = SolnBlk.Grid.Cell[i][j].Xc; Wl = SolnBlk.W[i][j];
	    if (SolnBlk.Grid.BCtypeN[i] == BC_ADIABATIC_WALL) {
	      // NORTH face of cell (i,j) is an ADIABATIC_WALL boundary.
	      elliptic_bc_flag = DIAMONDPATH_LEFT_TRIANGLE_HEATFLUX;
	      Wu = Knudsen_Layer_Adiabatic(Wl, SolnBlk.WoN[i].v, 
					   SolnBlk.Grid.nfaceN(i, j));
	    } else if (SolnBlk.Grid.BCtypeN[i] == BC_WALL_VISCOUS_ISOTHERMAL) {
	      // NORTH face of cell (i,j) is an ISOTHERMAL_WALL boundary.
	      elliptic_bc_flag = DIAMONDPATH_LEFT_TRIANGLE_HEATFLUX;
	      Wu = Knudsen_Layer_Isothermal(Wl, SolnBlk.WoN[i].v, SolnBlk.WoN[i].T(),
					    SolnBlk.Grid.nfaceN(i, j));
	    } else if (SolnBlk.Grid.BCtypeN[i] == BC_TEMPERATURE_SLIP) {
	      // NORTH face of cell (i,j) is a TEMPERATURE_SLIP boundary.
	      elliptic_bc_flag = DIAMONDPATH_LEFT_TRIANGLE_HEATFLUX;
	      Wu = Knudsen_Layer_Isothermal_Slip_T(Wl, SolnBlk.WoN[i].v, SolnBlk.WoN[i].T(),
						   SolnBlk.oldT_N[i],
						   SolnBlk.phi[i][j]^SolnBlk.dWdx[i][j],
						   SolnBlk.phi[i][j]^SolnBlk.dWdy[i][j],
						   SolnBlk.Grid.nfaceN(i, j));
	    } else {cout << "Error bad BC for elliptic part." << endl;}
	    switch(Input_Parameters.i_Heat_Reconstruction) {
	    case ELLIPTIC_RECONSTRUCTION_CARTESIAN :
	    case ELLIPTIC_RECONSTRUCTION_DIAMOND_PATH :
	      Xu = SolnBlk.Grid.Node[i  ][j+1].X;
	      Xd = SolnBlk.Grid.Node[i+1][j+1].X; Wd = Wu;
	      break;
	    case ELLIPTIC_RECONSTRUCTION_ARITHMETIC_AVERAGE :
	    case ELLIPTIC_RECONSTRUCTION_HYBRID :
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
	    switch(Input_Parameters.i_Heat_Reconstruction) {
	    case ELLIPTIC_RECONSTRUCTION_CARTESIAN :
	    case ELLIPTIC_RECONSTRUCTION_DIAMOND_PATH :
	      elliptic_bc_flag = DIAMONDPATH_QUADRILATERAL_RECONSTRUCTION;
	      Xu = SolnBlk.Grid.Node[i  ][j+1].X; Wu = SolnBlk.WnNW(i,j);
	      Xd = SolnBlk.Grid.Node[i+1][j+1].X; Wd = SolnBlk.WnNE(i,j);
	      break;
	    case ELLIPTIC_RECONSTRUCTION_ARITHMETIC_AVERAGE :
	    case ELLIPTIC_RECONSTRUCTION_HYBRID :
	      Xu = SolnBlk.Grid.xfaceN(i,j); Wu = HALF*(Wl + Wr);
	      dWdxl = SolnBlk.dWdx[i][j]; dWdxr = SolnBlk.dWdy[i][j+1];
	      dWdyl = SolnBlk.dWdy[i][j]; dWdyr = SolnBlk.dWdy[i][j+1];
	      break;
	    };
	  }
	  // Compute the NORTH face elliptic flux.
	  switch(Input_Parameters.i_Heat_Reconstruction) {
	  case ELLIPTIC_RECONSTRUCTION_CARTESIAN :
	  case ELLIPTIC_RECONSTRUCTION_DIAMOND_PATH :
	    Flux -= HeatFluxDiamondPath_n(SolnBlk.Grid.xfaceN(i,j),
					  Xl,Wl,Xd,Wd,Xr,Wr,Xu,Wu,
					  SolnBlk.Grid.nfaceN(i,j),
					  SolnBlk.Axisymmetric,
					  elliptic_bc_flag);
	    break;
	  case ELLIPTIC_RECONSTRUCTION_ARITHMETIC_AVERAGE :
	    Flux -= HeatFlux_n(Xu,Wu,HALF*(dWdxl+dWdxr),HALF*(dWdyl+dWdyr),
			       SolnBlk.Grid.nfaceN(i,j),SolnBlk.Axisymmetric);
	    break;
	  case ELLIPTIC_RECONSTRUCTION_HYBRID :
	    Flux -= HeatFluxHybrid_n(Xu,Wu,Xl,Wl,dWdxl,dWdyl,Xr,Wr,dWdxr,dWdyr,
				     SolnBlk.Grid.nfaceN(i,j),SolnBlk.Axisymmetric);
	    break;
	  };
	}
#endif	
          /* Evaluate cell-averaged solution changes. */
	
	SolnBlk.dUdt[i][j][0] -= 
	  Flux*SolnBlk.Grid.lfaceN(i, j)/
	  SolnBlk.Grid.Cell[i][j].A;
	SolnBlk.dUdt[i][j+1][0] += 
	  Flux*SolnBlk.Grid.lfaceS(i, j+1)/
	  SolnBlk.Grid.Cell[i][j+1].A;

	/* Save south and north face boundary flux. */
	
	if (j == SolnBlk.JCl-1) {
	  SolnBlk.FluxS[i] = -Flux*SolnBlk.Grid.lfaceS(i, j+1);
	} else if (j == SolnBlk.JCu) {
	  SolnBlk.FluxN[i] = Flux*SolnBlk.Grid.lfaceN(i, j);
	} /* endif */
	
      } /* endfor */
      
      SolnBlk.dUdt[i][SolnBlk.JCl-1][0] = Gaussian2D_U_VACUUM;
      SolnBlk.dUdt[i][SolnBlk.JCu+1][0] = Gaussian2D_U_VACUUM;
    } /* endfor */
    
    /* Residual successfully evaluated. */

}

/********************************************************
 * Routine: dUdt_Multistage_Explicit                    *
 *                                                      *
 * This routine determines the solution residuals for a *
 * given stage of a variety of multi-stage explicit     *
 * time integration schemes for a given solution block. *
 *                                                      *
 ********************************************************/
int dUdt_Multistage_Explicit(Gaussian2D_Quad_Block &SolnBlk,
                             const int i_stage,
			     Gaussian2D_Input_Parameters &Input_Parameters ) {

    int i, j, k_residual;
    double omega;
    Vector2D dX;
    Gaussian2D_pState Wl, Wr;
    Gaussian2D_cState Flux;

    Gaussian2D_pState Wu, Wd, dWdxl, dWdyl, dWdxr, dWdyr, W_temp;
    Vector2D Xl, Xr, Xu, Xd;
    int elliptic_bc_flag;

    /* Evaluate the solution residual for stage 
       i_stage of n_stage scheme. */

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
          SolnBlk.dUdt[SolnBlk.ICl-1][j][k_residual] = Gaussian2D_U_VACUUM;
       } else {
          SolnBlk.dUdt[SolnBlk.ICl-1][j][k_residual] = Gaussian2D_U_VACUUM;
       } /* endif */
    
       for ( i = SolnBlk.ICl-1 ; i <= SolnBlk.ICu ; ++i ) {
          if ( i_stage == 1 ) {
              SolnBlk.Uo[i+1][j] = SolnBlk.U[i+1][j];
              SolnBlk.dUdt[i+1][j][k_residual] = Gaussian2D_U_VACUUM;
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
                     SolnBlk.dUdt[i+1][j][k_residual] = Gaussian2D_U_VACUUM;
                  } /* endif */
                  break;
                case TIME_STEPPING_MULTISTAGE_OPTIMAL_SMOOTHING :
                  SolnBlk.dUdt[i+1][j][k_residual] = Gaussian2D_U_VACUUM;
                  break;
                default:
                  SolnBlk.dUdt[i+1][j][k_residual] = Gaussian2D_U_VACUUM;
                  break;
              } /* endswitch */
          } /* endif */
    
          if ( j > SolnBlk.JCl-1 && j < SolnBlk.JCu+1 ) {
    
             /* Evaluate the cell interface i-direction fluxes. */
    
	     if (i == SolnBlk.ICl-1 && 
                 (SolnBlk.Grid.BCtypeW[j] == BC_REFLECTION ||
                  SolnBlk.Grid.BCtypeW[j] == BC_CHARACTERISTIC ||
                  SolnBlk.Grid.BCtypeW[j] == BC_CHARACTERISTIC_VELOCITY ||
                  SolnBlk.Grid.BCtypeW[j] == BC_WALL_VISCOUS_ISOTHERMAL ||
                  SolnBlk.Grid.BCtypeW[j] == BC_TEMPERATURE_SLIP ||
                  SolnBlk.Grid.BCtypeW[j] == BC_ADIABATIC_WALL )) {
	       dX = SolnBlk.Grid.xfaceW(i+1, j)-SolnBlk.Grid.Cell[i+1][j].Xc;
	       Wr = SolnBlk.W[i+1][j] + 
		 (SolnBlk.phi[i+1][j]^SolnBlk.dWdx[i+1][j])*dX.x +
		 (SolnBlk.phi[i+1][j]^SolnBlk.dWdy[i+1][j])*dX.y;
	       if (SolnBlk.Grid.BCtypeW[j] == BC_REFLECTION) {
                 Wl = Reflect(Wr, SolnBlk.Grid.nfaceW(i+1, j));
	       } else if (SolnBlk.Grid.BCtypeW[j] == BC_ADIABATIC_WALL) {
		 Wl = Adiabatic_Wall(Wr, SolnBlk.WoW[j].v, 
		 	                 SolnBlk.Grid.nfaceW(i+1, j));
	       } else if (SolnBlk.Grid.BCtypeW[j] == BC_WALL_VISCOUS_ISOTHERMAL) {
		 Wl = Isothermal_Wall(Wr, SolnBlk.WoW[j].v, SolnBlk.WoW[j].T(), SolnBlk.Grid.nfaceW(i+1, j));
	       } else if (SolnBlk.Grid.BCtypeW[j] == BC_TEMPERATURE_SLIP) {
		 Wl = Isothermal_Wall_Slip_T(Wr, SolnBlk.WoW[j].v, SolnBlk.WoW[j].T(),
					     SolnBlk.oldT_W[j],
					     SolnBlk.phi[i+1][j]^SolnBlk.dWdx[i+1][j],
					     SolnBlk.phi[i+1][j]^SolnBlk.dWdy[i+1][j],
					     SolnBlk.Grid.nfaceW(i+1, j));
	       } else if (SolnBlk.Grid.BCtypeW[j] == BC_CHARACTERISTIC_VELOCITY) {
		 Wl = BC_Characteristic_Velocity(Wr, 
						 SolnBlk.WoW[j], 
						 SolnBlk.Grid.nfaceW(i+1, j));
	       } else {
		 Wl = BC_Characteristic_Pressure(Wr, 
						 SolnBlk.WoW[j], 
						 SolnBlk.Grid.nfaceW(i+1, j));
                } /* endif */
             } else if (i == SolnBlk.ICu && 
                        (SolnBlk.Grid.BCtypeE[j] == BC_REFLECTION ||
                         SolnBlk.Grid.BCtypeE[j] == BC_CHARACTERISTIC ||
                         SolnBlk.Grid.BCtypeE[j] == BC_CHARACTERISTIC_VELOCITY ||
                         SolnBlk.Grid.BCtypeE[j] == BC_WALL_VISCOUS_ISOTHERMAL ||
                         SolnBlk.Grid.BCtypeE[j] == BC_TEMPERATURE_SLIP ||
                         SolnBlk.Grid.BCtypeE[j] == BC_ADIABATIC_WALL )) {
	       dX = SolnBlk.Grid.xfaceE(i, j)-SolnBlk.Grid.Cell[i][j].Xc;
	       Wl = SolnBlk.W[i][j] + 
		 (SolnBlk.phi[i][j]^SolnBlk.dWdx[i][j])*dX.x +
		 (SolnBlk.phi[i][j]^SolnBlk.dWdy[i][j])*dX.y;
	       if (SolnBlk.Grid.BCtypeE[j] == BC_REFLECTION) {
		 Wr = Reflect(Wl, SolnBlk.Grid.nfaceE(i, j));
               } else if (SolnBlk.Grid.BCtypeE[j] == BC_ADIABATIC_WALL) {
                 Wr = Adiabatic_Wall(Wl, SolnBlk.WoE[j].v, 
                                         SolnBlk.Grid.nfaceE(i, j));
               } else if (SolnBlk.Grid.BCtypeE[j] == BC_WALL_VISCOUS_ISOTHERMAL) {
                 Wr = Isothermal_Wall(Wl, SolnBlk.WoE[j].v, SolnBlk.WoE[j].T(), SolnBlk.Grid.nfaceE(i, j));
               } else if (SolnBlk.Grid.BCtypeE[j] == BC_TEMPERATURE_SLIP) {
                 Wr = Isothermal_Wall_Slip_T(Wl, SolnBlk.WoE[j].v, SolnBlk.WoE[j].T(),
					     SolnBlk.oldT_E[j],
					     SolnBlk.phi[i][j]^SolnBlk.dWdx[i][j],
					     SolnBlk.phi[i][j]^SolnBlk.dWdy[i][j],
					     SolnBlk.Grid.nfaceE(i, j));
               } else if (SolnBlk.Grid.BCtypeE[j] == BC_CHARACTERISTIC_VELOCITY) {
                 Wr = BC_Characteristic_Velocity(Wl, 
                                                 SolnBlk.WoE[j], 
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

             switch(Input_Parameters.i_Flux_Function) {
               //case FLUX_FUNCTION_GODUNOV :
	       //Flux = FluxGodunov_n(Wl, Wr, SolnBlk.Grid.nfaceE(i, j));
	       //break;
               case FLUX_FUNCTION_ROE :
                 Flux = FluxRoe_n(Wl, Wr, SolnBlk.Grid.nfaceE(i, j));
                 break;
		 //case FLUX_FUNCTION_RUSANOV :
                 //Flux = FluxRusanov_n(Wl, Wr, SolnBlk.Grid.nfaceE(i, j));
                 //break;
	       case FLUX_FUNCTION_HLLE :
                 Flux = FluxHLLE_n(Wl, Wr, SolnBlk.Grid.nfaceE(i, j));
                 break;
		 //case FLUX_FUNCTION_LINDE :
                 //Flux = FluxLinde_n(Wl, Wr, SolnBlk.Grid.nfaceE(i, j));
                 //break;
		 //case FLUX_FUNCTION_HLLC :
                 //Flux = FluxHLLC_n(Wl, Wr, SolnBlk.Grid.nfaceE(i, j));
                 //break;
               default:
                 Flux = FluxRoe_n(Wl, Wr, SolnBlk.Grid.nfaceE(i, j));
                 break;
             } /* endswitch */

#ifdef _GAUSSIAN_HEAT_TRANSFER_
	     // Evaluate the cell interface i-direction ELLIPTIC flux if necessary.
	     if (SolnBlk.Heat_Transfer) {
	       // Determine the EAST face ELLIPTIC flux.
	       if (i == SolnBlk.ICl-1 && 
		   (SolnBlk.Grid.BCtypeW[j] == BC_ADIABATIC_WALL ||
		    SolnBlk.Grid.BCtypeW[j] == BC_WALL_VISCOUS_ISOTHERMAL ||
		    SolnBlk.Grid.BCtypeW[j] == BC_TEMPERATURE_SLIP)) {
		 // WEST face of cell (i+1,j) is a normal boundary.
		 Xr = SolnBlk.Grid.Cell[i+1][j].Xc; Wr = SolnBlk.W[i+1][j];
		 if (SolnBlk.Grid.BCtypeW[j] == BC_ADIABATIC_WALL) {
		   // WEST face of cell (i+1,j) is an ADIABATIC_WALL boundary.
		   elliptic_bc_flag = DIAMONDPATH_RIGHT_TRIANGLE_HEATFLUX;
		   Wu = Knudsen_Layer_Adiabatic(Wr, SolnBlk.WoW[j].v,
						SolnBlk.Grid.nfaceW(i+1, j));
		 } else if (SolnBlk.Grid.BCtypeW[j] == BC_WALL_VISCOUS_ISOTHERMAL) {
		   // WEST face of cell (i+1,j) is an ISOTHERMAL_WALL boundary.
		   elliptic_bc_flag = DIAMONDPATH_RIGHT_TRIANGLE_HEATFLUX;
		   Wu = Knudsen_Layer_Isothermal(Wr, SolnBlk.WoW[j].v, SolnBlk.WoW[j].T(),
						 SolnBlk.Grid.nfaceW(i+1, j));
		 } else if (SolnBlk.Grid.BCtypeW[j] == BC_TEMPERATURE_SLIP) {
		   // WEST face of cell (i+1,j) is a TEMPERATURE_SLIP boundary.
		   elliptic_bc_flag = DIAMONDPATH_RIGHT_TRIANGLE_HEATFLUX;
		   Wu = Knudsen_Layer_Isothermal_Slip_T(Wr, SolnBlk.WoW[j].v, SolnBlk.WoW[j].T(),
							SolnBlk.oldT_W[j],
							SolnBlk.phi[i+1][j]^SolnBlk.dWdx[i+1][j],
							SolnBlk.phi[i+1][j]^SolnBlk.dWdy[i+1][j],
							SolnBlk.Grid.nfaceW(i+1, j));
		 } else {cout << "Error bad BC for elliptic part" << endl;}
		 switch(Input_Parameters.i_Heat_Reconstruction) {
		 case ELLIPTIC_RECONSTRUCTION_CARTESIAN :
		 case ELLIPTIC_RECONSTRUCTION_DIAMOND_PATH :
		   Xu = SolnBlk.Grid.Node[i+1][j+1].X;
		   Xd = SolnBlk.Grid.Node[i+1][j  ].X; Wd = Wu;
		   break;
		 case ELLIPTIC_RECONSTRUCTION_ARITHMETIC_AVERAGE :
		 case ELLIPTIC_RECONSTRUCTION_HYBRID :
		   Xu = SolnBlk.Grid.xfaceW(i+1,j);
		   Xl = Xu; Wl = Wu;
		   dWdxr = SolnBlk.dWdx[i+1][j]; dWdxl = dWdxr;
		   dWdyr = SolnBlk.dWdy[i+1][j]; dWdyl = dWdyr;
		   break;
		 };
	       } else if (i == SolnBlk.ICu &&
			  (SolnBlk.Grid.BCtypeE[j] == BC_ADIABATIC_WALL ||
			   SolnBlk.Grid.BCtypeE[j] == BC_WALL_VISCOUS_ISOTHERMAL ||
			   SolnBlk.Grid.BCtypeE[j] == BC_TEMPERATURE_SLIP)) {
		 // EAST face of cell (i,j) is a normal boundary.
		 Xl = SolnBlk.Grid.Cell[i][j].Xc; Wl = SolnBlk.W[i][j];
		 if (SolnBlk.Grid.BCtypeE[j] == BC_ADIABATIC_WALL) {
		   // EAST face of cell (i,j) is an ADIABATIC_WALL boundary.
		   elliptic_bc_flag = DIAMONDPATH_LEFT_TRIANGLE_HEATFLUX;
		   Wu = Knudsen_Layer_Adiabatic(Wl, SolnBlk.WoE[j].v, 
						SolnBlk.Grid.nfaceE(i, j));
		 } else if (SolnBlk.Grid.BCtypeE[j] == BC_WALL_VISCOUS_ISOTHERMAL) {
		   // EAST face of cell (i,j) is an ISOTHERMAL_WALL boundary.
		   elliptic_bc_flag = DIAMONDPATH_LEFT_TRIANGLE_HEATFLUX;
		   Wu = Knudsen_Layer_Isothermal(Wl, SolnBlk.WoE[j].v, SolnBlk.WoE[j].T(),
						 SolnBlk.Grid.nfaceE(i, j));
		 } else if (SolnBlk.Grid.BCtypeE[j] == BC_TEMPERATURE_SLIP) {
		   // EAST face of cell (i,j) is a TEMPERATURE_SLIP boundary.
		   elliptic_bc_flag = DIAMONDPATH_LEFT_TRIANGLE_HEATFLUX;
		   Wu = Knudsen_Layer_Isothermal_Slip_T(Wl, SolnBlk.WoE[j].v, SolnBlk.WoE[j].T(),
							SolnBlk.oldT_E[j],
							SolnBlk.phi[i][j]^SolnBlk.dWdx[i][j],
							SolnBlk.phi[i][j]^SolnBlk.dWdy[i][j],
							SolnBlk.Grid.nfaceE(i, j));
		 } else {cout << "Error bad BC for elliptic part" << endl;}
		 switch(Input_Parameters.i_Heat_Reconstruction) {
		 case ELLIPTIC_RECONSTRUCTION_CARTESIAN :
		 case ELLIPTIC_RECONSTRUCTION_DIAMOND_PATH :
		   Xu = SolnBlk.Grid.Node[i+1][j+1].X;
		   Xd = SolnBlk.Grid.Node[i+1][j  ].X; Wd = Wu;
		   break;
		 case ELLIPTIC_RECONSTRUCTION_ARITHMETIC_AVERAGE :
		 case ELLIPTIC_RECONSTRUCTION_HYBRID :
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
		 switch(Input_Parameters.i_Heat_Reconstruction) {
		 case ELLIPTIC_RECONSTRUCTION_CARTESIAN :
		 case ELLIPTIC_RECONSTRUCTION_DIAMOND_PATH :
		   elliptic_bc_flag = DIAMONDPATH_QUADRILATERAL_RECONSTRUCTION;
		   Xu = SolnBlk.Grid.Node[i+1][j+1].X; Wu = SolnBlk.WnNE(i,j);
		   Xd = SolnBlk.Grid.Node[i+1][j  ].X; Wd = SolnBlk.WnSE(i,j);
		   break;
		 case ELLIPTIC_RECONSTRUCTION_ARITHMETIC_AVERAGE :
		 case ELLIPTIC_RECONSTRUCTION_HYBRID :
		   Xu = SolnBlk.Grid.xfaceE(i,j); Wu = HALF*(Wl + Wr);
		   dWdxl = SolnBlk.dWdx[i][j]; dWdxr = SolnBlk.dWdy[i+1][j];
		   dWdyl = SolnBlk.dWdy[i][j]; dWdyr = SolnBlk.dWdy[i+1][j];
		   break;
		 };
	       }
	       // Compute the EAST face elliptic flux.
	       switch(Input_Parameters.i_Heat_Reconstruction) {
	       case ELLIPTIC_RECONSTRUCTION_CARTESIAN :
	       case ELLIPTIC_RECONSTRUCTION_DIAMOND_PATH :
		 Flux -= HeatFluxDiamondPath_n(SolnBlk.Grid.xfaceE(i,j),
					       Xl,Wl,Xd,Wd,Xr,Wr,Xu,Wu,
					       SolnBlk.Grid.nfaceE(i,j),
					       SolnBlk.Axisymmetric,
					       elliptic_bc_flag);
		 break;
	       case ELLIPTIC_RECONSTRUCTION_ARITHMETIC_AVERAGE :
		 Flux -= HeatFlux_n(Xu,Wu,HALF*(dWdxl+dWdxr),HALF*(dWdyl+dWdyr),
				    SolnBlk.Grid.nfaceE(i,j),SolnBlk.Axisymmetric);
		 break;
	       case ELLIPTIC_RECONSTRUCTION_HYBRID :
		 Flux -= HeatFluxHybrid_n(Xu,Wu,Xl,Wl,dWdxl,dWdyl,Xr,Wr,dWdxr,dWdyr,
					  SolnBlk.Grid.nfaceE(i,j),SolnBlk.Axisymmetric);
		 break;
	       };
	     }
#endif
             /* Evaluate cell-averaged solution changes. */

             SolnBlk.dUdt[i][j][k_residual] -= 
                (Input_Parameters.CFL_Number*SolnBlk.dt[i][j])*
                Flux*SolnBlk.Grid.lfaceE(i, j)/
                SolnBlk.Grid.Cell[i][j].A;
             SolnBlk.dUdt[i+1][j][k_residual] += 
                (Input_Parameters.CFL_Number*SolnBlk.dt[i+1][j])*
                Flux*SolnBlk.Grid.lfaceW(i+1, j)/
                SolnBlk.Grid.Cell[i+1][j].A;

             /* Include axisymmetric source terms as required. */

	     if (SolnBlk.Axisymmetric) {
               SolnBlk.dUdt[i][j][k_residual] += 
                  (Input_Parameters.CFL_Number*SolnBlk.dt[i][j])*
	          S(SolnBlk.W[i][j], SolnBlk.Grid.Cell[i][j].Xc);
             } /* endif */

             /* Save west and east face boundary flux. */

             if (i == SolnBlk.ICl-1) {
                SolnBlk.FluxW[j] = -Flux*SolnBlk.Grid.lfaceW(i+1, j);
             } else if (i == SolnBlk.ICu) {
                SolnBlk.FluxE[j] = Flux*SolnBlk.Grid.lfaceE(i, j);
             } /* endif */ 

          } /* endif */
       } /* endfor */
    
       if ( j > SolnBlk.JCl-1 && j < SolnBlk.JCu+1 ) {
          SolnBlk.dUdt[SolnBlk.ICl-1][j][k_residual] = Gaussian2D_U_VACUUM;
          SolnBlk.dUdt[SolnBlk.ICu+1][j][k_residual] = Gaussian2D_U_VACUUM;
       } /* endif */
    } /* endfor */
    
    // Add j-direction (eta-direction) fluxes.
    for ( i = SolnBlk.ICl ; i <= SolnBlk.ICu ; ++i ) {
       for ( j  = SolnBlk.JCl-1 ; j <= SolnBlk.JCu ; ++j ) {
    
          /* Evaluate the cell interface j-direction fluxes. */
         
	  if (j == SolnBlk.JCl-1 && 
              (SolnBlk.Grid.BCtypeS[i] == BC_REFLECTION ||
               SolnBlk.Grid.BCtypeS[i] == BC_CHARACTERISTIC ||
               SolnBlk.Grid.BCtypeS[i] == BC_CHARACTERISTIC_VELOCITY ||
               SolnBlk.Grid.BCtypeS[i] == BC_WALL_VISCOUS_ISOTHERMAL ||
               SolnBlk.Grid.BCtypeS[i] == BC_TEMPERATURE_SLIP ||
               SolnBlk.Grid.BCtypeS[i] == BC_ADIABATIC_WALL )) {
	     dX = SolnBlk.Grid.xfaceS(i, j+1)-SolnBlk.Grid.Cell[i][j+1].Xc;
	     Wr = SolnBlk.W[i][j+1] +
	          (SolnBlk.phi[i][j+1]^SolnBlk.dWdx[i][j+1])*dX.x +
	          (SolnBlk.phi[i][j+1]^SolnBlk.dWdy[i][j+1])*dX.y;
             if (SolnBlk.Grid.BCtypeS[i] == BC_REFLECTION) {
               Wl = Reflect(Wr, SolnBlk.Grid.nfaceS(i, j+1));
             } else if (SolnBlk.Grid.BCtypeS[i] == BC_ADIABATIC_WALL) {
	       Wl = Adiabatic_Wall(Wr, SolnBlk.WoS[i].v, 
				   SolnBlk.Grid.nfaceS(i, j+1));
             } else if (SolnBlk.Grid.BCtypeS[i] == BC_WALL_VISCOUS_ISOTHERMAL) {
	       Wl = Isothermal_Wall(Wr, SolnBlk.WoS[i].v, SolnBlk.WoS[i].T(), SolnBlk.Grid.nfaceS(i, j+1));
             } else if (SolnBlk.Grid.BCtypeS[i] == BC_TEMPERATURE_SLIP) {
	       Wl = Isothermal_Wall_Slip_T(Wr, SolnBlk.WoS[i].v, SolnBlk.WoS[i].T(),
					   SolnBlk.oldT_S[i],
					   SolnBlk.phi[i][j+1]^SolnBlk.dWdx[i][j+1],
					   SolnBlk.phi[i][j+1]^SolnBlk.dWdy[i][j+1],
					   SolnBlk.Grid.nfaceS(i, j+1));
             } else if (SolnBlk.Grid.BCtypeS[i] == BC_CHARACTERISTIC_VELOCITY) {
               Wl = BC_Characteristic_Velocity(Wr, 
                                               SolnBlk.WoS[i], 
                                               SolnBlk.Grid.nfaceS(i, j+1));
	     } else {
               Wl = BC_Characteristic_Pressure(Wr, 
                                               SolnBlk.WoS[i], 
                                               SolnBlk.Grid.nfaceS(i, j+1));
             } /* endif */
          } else if (j == SolnBlk.JCu && 
                     (SolnBlk.Grid.BCtypeN[i] == BC_REFLECTION ||
                      SolnBlk.Grid.BCtypeN[i] == BC_CHARACTERISTIC ||
                      SolnBlk.Grid.BCtypeN[i] == BC_CHARACTERISTIC_VELOCITY ||
                      SolnBlk.Grid.BCtypeN[i] == BC_WALL_VISCOUS_ISOTHERMAL ||
                      SolnBlk.Grid.BCtypeN[i] == BC_TEMPERATURE_SLIP ||
                      SolnBlk.Grid.BCtypeN[i] == BC_ADIABATIC_WALL )) {
	    dX = SolnBlk.Grid.xfaceN(i, j)-SolnBlk.Grid.Cell[i][j].Xc;
	    Wl = SolnBlk.W[i][j] + 
	      (SolnBlk.phi[i][j]^SolnBlk.dWdx[i][j])*dX.x +
	      (SolnBlk.phi[i][j]^SolnBlk.dWdy[i][j])*dX.y;
	    if (SolnBlk.Grid.BCtypeN[i] == BC_REFLECTION) {
	      Wr = Reflect(Wl, SolnBlk.Grid.nfaceN(i, j));
	    } else if (SolnBlk.Grid.BCtypeN[i] == BC_ADIABATIC_WALL) {
	      Wr = Adiabatic_Wall(Wl, SolnBlk.WoN[i].v,
	      			      SolnBlk.Grid.nfaceN(i, j));
	    } else if (SolnBlk.Grid.BCtypeN[i] == BC_WALL_VISCOUS_ISOTHERMAL) {
	      Wr = Isothermal_Wall(Wl, SolnBlk.WoN[i].v, SolnBlk.WoN[i].T(), SolnBlk.Grid.nfaceN(i, j));
	    } else if (SolnBlk.Grid.BCtypeN[i] == BC_TEMPERATURE_SLIP) {
	      Wr = Isothermal_Wall_Slip_T(Wl, SolnBlk.WoN[i].v, SolnBlk.WoN[i].T(),
					  SolnBlk.oldT_N[i],
					  SolnBlk.phi[i][j]^SolnBlk.dWdx[i][j],
					  SolnBlk.phi[i][j]^SolnBlk.dWdy[i][j],
					  SolnBlk.Grid.nfaceN(i, j));
	    } else if (SolnBlk.Grid.BCtypeN[i] == BC_CHARACTERISTIC_VELOCITY) {
	      Wr = BC_Characteristic_Velocity(Wl, 
					      SolnBlk.WoN[i], 
					      SolnBlk.Grid.nfaceN(i, j));
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

          switch(Input_Parameters.i_Flux_Function) {
            //case FLUX_FUNCTION_GODUNOV :
	    //Flux = FluxGodunov_n(Wl, Wr, SolnBlk.Grid.nfaceN(i, j));
	    //break;
            case FLUX_FUNCTION_ROE :
              Flux = FluxRoe_n(Wl, Wr, SolnBlk.Grid.nfaceN(i, j));
              break;
	      //case FLUX_FUNCTION_RUSANOV :
              //Flux = FluxRusanov_n(Wl, Wr, SolnBlk.Grid.nfaceN(i, j));
              //break;
	    case FLUX_FUNCTION_HLLE :
              Flux = FluxHLLE_n(Wl, Wr, SolnBlk.Grid.nfaceN(i, j));
              break;
	      //case FLUX_FUNCTION_LINDE :
              //Flux = FluxLinde_n(Wl, Wr, SolnBlk.Grid.nfaceN(i, j));
              //break;
	      //case FLUX_FUNCTION_HLLC :
              //Flux = FluxHLLC_n(Wl, Wr, SolnBlk.Grid.nfaceN(i, j));
              //break;
            default:
              Flux = FluxRoe_n(Wl, Wr, SolnBlk.Grid.nfaceN(i, j));
              break;
	  } /* endswitch */

#ifdef _GAUSSIAN_HEAT_TRANSFER_
	  // Evaluate the cell interface j-direction ELLIPTIC flux if necessary.
	  if (SolnBlk.Heat_Transfer) {
	    // Determine the NORTH face ELLIPTIC flux.
	    if (j == SolnBlk.JCl-1 && 
		(SolnBlk.Grid.BCtypeS[i] == BC_ADIABATIC_WALL ||
		 SolnBlk.Grid.BCtypeS[i] == BC_WALL_VISCOUS_ISOTHERMAL ||
		 SolnBlk.Grid.BCtypeS[i] == BC_TEMPERATURE_SLIP)) {
	      // SOUTH face of cell (i,j+1) is a normal boundary.
	      Xr = SolnBlk.Grid.Cell[i][j+1].Xc; Wr = SolnBlk.W[i][j+1];
	      dX = SolnBlk.Grid.xfaceS(i, j+1)-Xr;
	      W_temp = Wr + (SolnBlk.phi[i][j+1]^SolnBlk.dWdx[i][j+1])*dX.x +
	                    (SolnBlk.phi[i][j+1]^SolnBlk.dWdy[i][j+1])*dX.y;
	      if (SolnBlk.Grid.BCtypeS[i] == BC_ADIABATIC_WALL) {
		// SOUTH face of cell (i,j+1) is an ADIABATIC_WALL boundary.
		elliptic_bc_flag = DIAMONDPATH_RIGHT_TRIANGLE_HEATFLUX;
		Wu = Knudsen_Layer_Adiabatic(W_temp, SolnBlk.WoS[i].v, 
					     SolnBlk.Grid.nfaceS(i, j+1));
	      } else if (SolnBlk.Grid.BCtypeS[i] == BC_WALL_VISCOUS_ISOTHERMAL) {
		// SOUTH face of cell (i,j+1) is an ISOTHERMAL_WALL boundary.
		elliptic_bc_flag = DIAMONDPATH_RIGHT_TRIANGLE_HEATFLUX;
		Wu = Knudsen_Layer_Isothermal(W_temp, SolnBlk.WoS[i].v, SolnBlk.WoS[i].T(),
					      SolnBlk.Grid.nfaceS(i, j+1));
	      } else if (SolnBlk.Grid.BCtypeS[i] == BC_TEMPERATURE_SLIP) {
		// SOUTH face of cell (i,j+1) is a TEMPERATURE_SLIP boundary.
		elliptic_bc_flag = DIAMONDPATH_RIGHT_TRIANGLE_HEATFLUX;
		Wu = Knudsen_Layer_Isothermal_Slip_T(W_temp, SolnBlk.WoS[i].v, SolnBlk.WoS[i].T(),
						     SolnBlk.oldT_S[i],
						     SolnBlk.phi[i][j+1]^SolnBlk.dWdx[i][j+1],
						     SolnBlk.phi[i][j+1]^SolnBlk.dWdy[i][j+1],
						     SolnBlk.Grid.nfaceS(i, j+1));
	      } else {cout << "Error bad BC for elliptic part." << endl;}
	      switch(Input_Parameters.i_Heat_Reconstruction) {
	      case ELLIPTIC_RECONSTRUCTION_CARTESIAN :
	      case ELLIPTIC_RECONSTRUCTION_DIAMOND_PATH :
		Xu = SolnBlk.Grid.Node[i  ][j+1].X;
		Xd = SolnBlk.Grid.Node[i+1][j+1].X; Wd = Wu;
		break;
	      case ELLIPTIC_RECONSTRUCTION_ARITHMETIC_AVERAGE :
	      case ELLIPTIC_RECONSTRUCTION_HYBRID :
		Xu = SolnBlk.Grid.xfaceS(i,j+1);
		Xl = Xu; Wl = Wu;
		dWdxr = SolnBlk.dWdx[i][j+1]; dWdxl = dWdxr;
		dWdyr = SolnBlk.dWdy[i][j+1]; dWdyl = dWdyr;
		break;
	      };
	    } else if (j == SolnBlk.JCu && 
		       (SolnBlk.Grid.BCtypeN[i] == BC_ADIABATIC_WALL ||
			SolnBlk.Grid.BCtypeN[i] == BC_WALL_VISCOUS_ISOTHERMAL ||
			SolnBlk.Grid.BCtypeN[i] == BC_TEMPERATURE_SLIP)) {
	      // NORTH face of cell (i,j) is a normal boundary.
	      Xl = SolnBlk.Grid.Cell[i][j].Xc; Wl = SolnBlk.W[i][j];
	      dX = SolnBlk.Grid.xfaceN(i, j)- Xl;
	      W_temp = Wl + (SolnBlk.phi[i][j]^SolnBlk.dWdx[i][j])*dX.x +
		            (SolnBlk.phi[i][j]^SolnBlk.dWdy[i][j])*dX.y;
	      if (SolnBlk.Grid.BCtypeN[i] == BC_ADIABATIC_WALL) {
		// NORTH face of cell (i,j) is an ADIABATIC_WALL boundary.
		elliptic_bc_flag = DIAMONDPATH_LEFT_TRIANGLE_HEATFLUX;
		Wu = Knudsen_Layer_Adiabatic(W_temp, SolnBlk.WoN[i].v, 
					     SolnBlk.Grid.nfaceN(i, j));
	      } else if (SolnBlk.Grid.BCtypeN[i] == BC_WALL_VISCOUS_ISOTHERMAL) {
		// NORTH face of cell (i,j) is an ISOTHERMAL_WALL boundary.
		elliptic_bc_flag = DIAMONDPATH_LEFT_TRIANGLE_HEATFLUX;
		Wu = Knudsen_Layer_Isothermal(W_temp, SolnBlk.WoN[i].v, SolnBlk.WoN[i].T(),
					      SolnBlk.Grid.nfaceN(i, j));
	      } else if (SolnBlk.Grid.BCtypeN[i] == BC_TEMPERATURE_SLIP) {
		// NORTH face of cell (i,j) is a TEMPERATURE_SLIP boundary.
		elliptic_bc_flag = DIAMONDPATH_LEFT_TRIANGLE_HEATFLUX;
		Wu = Knudsen_Layer_Isothermal_Slip_T(W_temp, SolnBlk.WoN[i].v, SolnBlk.WoN[i].T(),
						     SolnBlk.oldT_N[i],
						     SolnBlk.phi[i][j]^SolnBlk.dWdx[i][j],
						     SolnBlk.phi[i][j]^SolnBlk.dWdy[i][j],
						     SolnBlk.Grid.nfaceN(i, j));
	      } else {cout << "Error bad BC for elliptic part." << endl;}
	      switch(Input_Parameters.i_Heat_Reconstruction) {
	      case ELLIPTIC_RECONSTRUCTION_CARTESIAN :
	      case ELLIPTIC_RECONSTRUCTION_DIAMOND_PATH :
		Xu = SolnBlk.Grid.Node[i  ][j+1].X;
		Xd = SolnBlk.Grid.Node[i+1][j+1].X; Wd = Wu;
		break;
	      case ELLIPTIC_RECONSTRUCTION_ARITHMETIC_AVERAGE :
	      case ELLIPTIC_RECONSTRUCTION_HYBRID :
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
	      switch(Input_Parameters.i_Heat_Reconstruction) {
	      case ELLIPTIC_RECONSTRUCTION_CARTESIAN :
	      case ELLIPTIC_RECONSTRUCTION_DIAMOND_PATH :
		elliptic_bc_flag = DIAMONDPATH_QUADRILATERAL_RECONSTRUCTION;
		Xu = SolnBlk.Grid.Node[i  ][j+1].X; Wu = SolnBlk.WnNW(i,j);
		Xd = SolnBlk.Grid.Node[i+1][j+1].X; Wd = SolnBlk.WnNE(i,j);
		break;
	      case ELLIPTIC_RECONSTRUCTION_ARITHMETIC_AVERAGE :
	      case ELLIPTIC_RECONSTRUCTION_HYBRID :
		Xu = SolnBlk.Grid.xfaceN(i,j); Wu = HALF*(Wl + Wr);
		dWdxl = SolnBlk.dWdx[i][j]; dWdxr = SolnBlk.dWdy[i][j+1];
		dWdyl = SolnBlk.dWdy[i][j]; dWdyr = SolnBlk.dWdy[i][j+1];
		break;
	      };
	    }
	    // Compute the NORTH face elliptic flux.
	    switch(Input_Parameters.i_Heat_Reconstruction) {
	    case ELLIPTIC_RECONSTRUCTION_CARTESIAN :
	    case ELLIPTIC_RECONSTRUCTION_DIAMOND_PATH :
	      Flux -= HeatFluxDiamondPath_n(SolnBlk.Grid.xfaceN(i,j),
					    Xl,Wl,Xd,Wd,Xr,Wr,Xu,Wu,
					    SolnBlk.Grid.nfaceN(i,j),
					    SolnBlk.Axisymmetric,
					    elliptic_bc_flag);
	      break;
	    case ELLIPTIC_RECONSTRUCTION_ARITHMETIC_AVERAGE :
	      Flux -= HeatFlux_n(Xu,Wu,HALF*(dWdxl+dWdxr),HALF*(dWdyl+dWdyr),
				 SolnBlk.Grid.nfaceN(i,j),SolnBlk.Axisymmetric);
	      break;
	    case ELLIPTIC_RECONSTRUCTION_HYBRID :
	      Flux -= HeatFluxHybrid_n(Xu,Wu,Xl,Wl,dWdxl,dWdyl,Xr,Wr,dWdxr,dWdyr,
				       SolnBlk.Grid.nfaceN(i,j),SolnBlk.Axisymmetric);
	      break;
	    };
	  }
#endif

          /* Evaluate cell-averaged solution changes. */
    
          SolnBlk.dUdt[i][j][k_residual] -= 
             (Input_Parameters.CFL_Number*SolnBlk.dt[i][j])*
             Flux*SolnBlk.Grid.lfaceN(i, j)/
             SolnBlk.Grid.Cell[i][j].A;
          SolnBlk.dUdt[i][j+1][k_residual] += 
             (Input_Parameters.CFL_Number*SolnBlk.dt[i][j+1])*
             Flux*SolnBlk.Grid.lfaceS(i, j+1)/
             SolnBlk.Grid.Cell[i][j+1].A;

          /* Save south and north face boundary flux. */

          if (j == SolnBlk.JCl-1) {
             SolnBlk.FluxS[i] = -Flux*SolnBlk.Grid.lfaceS(i, j+1);
          } else if (j == SolnBlk.JCu) {
             SolnBlk.FluxN[i] = Flux*SolnBlk.Grid.lfaceN(i, j);
          } /* endif */

       } /* endfor */
    
       SolnBlk.dUdt[i][SolnBlk.JCl-1][k_residual] = Gaussian2D_U_VACUUM;
       SolnBlk.dUdt[i][SolnBlk.JCu+1][k_residual] = Gaussian2D_U_VACUUM;
    } /* endfor */
    
    /* Residual for the stage successfully calculated. */

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
int Update_Solution_Multistage_Explicit(Gaussian2D_Quad_Block &SolnBlk,
                                        const int i_stage,
					Gaussian2D_Input_Parameters &Input_Parameters) {

    int i, j, k, k_residual, n_residual_reduction;
    double omega, residual_reduction_factor;

    /* Perform update of solution variables for stage 
       i_stage of n_stage scheme. */

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

	 SolnBlk.U[i][j] = SolnBlk.Uo[i][j] + omega*SolnBlk.dUdt[i][j][k_residual];

	 SolnBlk.W[i][j] = W(SolnBlk.U[i][j]);
	 SolnBlk.W[i][j].relax(Input_Parameters.CFL_Number*SolnBlk.dt[i][j],
                               i_stage,W(SolnBlk.Uo[i][j]));

	 SolnBlk.U[i][j] = U(SolnBlk.W[i][j]);	 
	 
	 if (Input_Parameters.Local_Time_Stepping == GLOBAL_TIME_STEPPING && 
	     (SolnBlk.U[i][j].invalid())){

	     //SolnBlk.W[i][j].make_valid();
	     //SolnBlk.U[i][j] = SolnBlk.W[i][j];
	     //cout << "(" << SolnBlk.Grid.Cell[i][j].Xc << ")";

	     /*
	     cout << "\n " << CFFC_Name() << " Gaussian2D ERROR: Invalid State: \n"
		  << " cell = (" << i << ", " << j << ") " 
		  << " X = " << SolnBlk.Grid.Cell[i][j].Xc << "\n Wo = "
		  << W(SolnBlk.Uo[i][j]) << "\n DetP for Wo = "
		  << W(SolnBlk.Uo[i][j]).DetP() << "\n W = " 
		  << W(SolnBlk.U[i][j]) << "\n DetP = "
		  << W(SolnBlk.U[i][j]).DetP() << "\n U = "
		  << SolnBlk.U[i][j] << "\n dUdt = " 
		  << SolnBlk.dUdt[i][j][k_residual] << "\n";

	     return (i);
	     */
	 } else if ((Input_Parameters.Local_Time_Stepping == SCALAR_LOCAL_TIME_STEPPING) &&
	            (SolnBlk.U[i][j].invalid())){

	   residual_reduction_factor = ONE;
	   cout << "*";
	   for (n_residual_reduction = 1; n_residual_reduction <= 15; ++n_residual_reduction) {
	     
	     residual_reduction_factor = HALF*residual_reduction_factor;
	     //SolnBlk.dt[i][j] = residual_reduction_factor*SolnBlk.dt[i][j];
	     //SolnBlk.dUdt[i][j][k_residual] = residual_reduction_factor*SolnBlk.dUdt[i][j][k_residual];
	     
	     SolnBlk.U[i][j] = SolnBlk.Uo[i][j] + residual_reduction_factor*omega*SolnBlk.dUdt[i][j][k_residual];
	     
	     SolnBlk.W[i][j] = W(SolnBlk.U[i][j]);
	     SolnBlk.W[i][j].relax(Input_Parameters.CFL_Number*residual_reduction_factor*SolnBlk.dt[i][j],
                                   i_stage,W(SolnBlk.Uo[i][j]));
	     SolnBlk.U[i][j] = U(SolnBlk.W[i][j]);	 
	     
	     /*
	     SolnBlk.dt[i][j] = HALF*SolnBlk.dt[i][j];
	     SolnBlk.dUdt[i][j][k_residual] = HALF*SolnBlk.dUdt[i][j][k_residual];
	     
	     SolnBlk.U[i][j] = SolnBlk.Uo[i][j] + omega*SolnBlk.dUdt[i][j][k_residual];
	     
	     SolnBlk.W[i][j] = W(SolnBlk.U[i][j]);
	     SolnBlk.W[i][j].relax(Input_Parameters.CFL_Number*SolnBlk.dt[i][j],
                                   i_stage,W(SolnBlk.Uo[i][j]));
	     SolnBlk.U[i][j] = U(SolnBlk.W[i][j]);	 
	     */
	     //	     cout << "*" << n_residual_reduction << "*";
	     if (!SolnBlk.U[i][j].invalid()){
		 break;
	     } /* endif */
	   } /* end for */
	   
	   if (SolnBlk.U[i][j].invalid()){

	     cout << "(" << SolnBlk.Grid.Cell[i][j].Xc << ")";
	     SolnBlk.U[i][j] = SolnBlk.Uo[i][j];
	     SolnBlk.W[i][j] = W(SolnBlk.U[i][j]);
	     //SolnBlk.W[i][j].make_valid();
	     //SolnBlk.U[i][j] = SolnBlk.W[i][j];
	   }
	   
	     /*
	     cout << "\n " << CFFC_Name() << " Gaussian2D ERROR: Invalid State: \n"
		  << " cell = (" << i << ", " << j << ") " 
		  << " X = " << SolnBlk.Grid.Cell[i][j].Xc << "\n Wo = "
		  << W(SolnBlk.Uo[i][j]) << "\n DetP for Wo = "
		  << W(SolnBlk.Uo[i][j]).DetP() << "\n W = " 
		  << W(SolnBlk.U[i][j]) << "\n DetP = "
		  << W(SolnBlk.U[i][j]).DetP() << "\n U = "
		  << SolnBlk.U[i][j] << "\n dUdt = " 
		  << SolnBlk.dUdt[i][j][k_residual] << "\n";

	     return (i);
	     */
	   //} /* endif */
	 } /* end if */	 
       } /* endfor */    
    } /* endfor */
    
    /* Solution successfully updated. */
    
    return (0);   
}


/********************************************************
 * Routine: Ramp_up_Reference_Mach_Number               *
 *                                                      *
 * This routine ramps up the reference mach number      *
 *                                                      *
 ********************************************************/

void Ramp_up_Reference_Mach_Number(Gaussian2D_Quad_Block &SolnBlk,
				   Gaussian2D_Input_Parameters &Input_Parameters,
				   const int Number_of_Time_Steps) {

  int i,j;

  Gaussian2D_pState Wo_Ramped(Input_Parameters.Wo);

  Wo_Ramped.v.x = Wo_Ramped.v.x+Input_Parameters.Ramp_by_Mach_Number*
    ((double)(Number_of_Time_Steps+1)/(double)Input_Parameters.Number_of_Time_Steps_to_Ramp)*
    Wo_Ramped.sound()*cos(TWO*PI*Input_Parameters.Flow_Angle/360.0);

  Wo_Ramped.v.y = Wo_Ramped.v.y+Input_Parameters.Ramp_by_Mach_Number*
    ((double)(Number_of_Time_Steps+1)/(double)Input_Parameters.Number_of_Time_Steps_to_Ramp)*
    Wo_Ramped.sound()*sin(TWO*PI*Input_Parameters.Flow_Angle/360.0);

//  for (j = SolnBlk.JCl-SolnBlk.Nghost ; j <= SolnBlk.JCu+SolnBlk.Nghost ; ++j ) {
//    if (j >= SolnBlk.JCl && j <= SolnBlk.JCu) {
//      SolnBlk.WoW[j] = Wo_Ramped;
//      SolnBlk.WoE[j] = Wo_Ramped;
//    } else if (j < SolnBlk.JCl) {
//      SolnBlk.WoW[j] = Wo_Ramped;
//      SolnBlk.WoE[j] = Wo_Ramped;
//    } else {
//      SolnBlk.WoW[j] = Wo_Ramped;
//      SolnBlk.WoE[j] = Wo_Ramped;
//    } /* endif */
//
//    //Set Wall velocities for Adiabatic walls
//    if (SolnBlk.Grid.BCtypeW[j] == BC_ADIABATIC_WALL){
//      SolnBlk.WoW[j].v.y = 0.0;
//      SolnBlk.WoW[j].v.x = 0.0;
//    } /* endif */
//    if (SolnBlk.Grid.BCtypeE[j] == BC_ADIABATIC_WALL){
//      SolnBlk.WoE[j].v.y = 0.0;
//      SolnBlk.WoE[j].v.x = 0.0;
//    } /* endif */
//  } /* endfor */


  for ( i = SolnBlk.ICl-SolnBlk.Nghost ; i <= SolnBlk.ICu+SolnBlk.Nghost ; ++i ) {
    if (i >= SolnBlk.ICl && i <= SolnBlk.ICu) {
      //SolnBlk.WoS[i] = Wo_Ramped;
      SolnBlk.WoN[i] = Wo_Ramped;
    } else if (i < SolnBlk.ICl) {
      //SolnBlk.WoS[i] = Wo_Ramped;
      SolnBlk.WoN[i] = Wo_Ramped;
    } else {
      //SolnBlk.WoS[i] = Wo_Ramped;
      SolnBlk.WoN[i] = Wo_Ramped;
    } /* endif */
    
//    //Set Wall velocities for Adiabatic walls
//    if (SolnBlk.Grid.BCtypeS[i] == BC_ADIABATIC_WALL){
//      SolnBlk.WoS[i].v.y = 0.0;
//      if (Input_Parameters.i_Grid == GRID_ADIABATIC_COUETTE) {
//	SolnBlk.WoS[i].v.x =  - Input_Parameters.Couette_Plate_Velocity;
//      } else {
//	SolnBlk.WoS[i].v.x = 0.0;
//      } /* endif */
//    } /* endif */
    if (SolnBlk.Grid.BCtypeN[i] == BC_ADIABATIC_WALL){
      SolnBlk.WoN[i].v.y = 0.0;
      if (Input_Parameters.i_Grid == GRID_ADIABATIC_COUETTE) {
	SolnBlk.WoN[i].v.x = Input_Parameters.Couette_Plate_Velocity;
      } else {
	SolnBlk.WoN[i].v.x = 0.0;
      } /* endif */
    } /* endif */
  } /* endfor */

}

/********************************************************
 * Routine: Count_south_adiabatic_wall_cells            *
 *                                                      *
 *  Used for free-molecular stuff                       *
 *                                                      *
 ********************************************************/
extern int Count_south_adiabatic_wall_cells(Gaussian2D_Quad_Block &SolnBlk) {

  int i, num(0);

  for(i = SolnBlk.Grid.ICl; i <= SolnBlk.Grid.ICu; i++) {

    if(SolnBlk.Grid.BCtypeS[i] == BC_ADIABATIC_WALL){
      num++;
    }

  }

  return num;
}
