/* AdvectDiffuse2DQuadSingleBlock.cc:  Single-Block Versions of Subroutines for
                                       2D Advection Diffusion Equation
                                       Multi-Block Quadrilateral Mesh 
                                       Solution Classes. */

/* Include 2D advection diffusion equation quadrilateral mesh solution header file. */

#ifndef _ADVECTDIFFUSE2D_QUAD_INCLUDED
#include "AdvectDiffuse2DQuad.h"
#endif // _ADVECTDIFFUSE2D_QUAD_INCLUDED

/**************************************************************************
 * AdvectDiffuse2D_Quad_Block -- Single Block External Subroutines.       *
 **************************************************************************/

/********************************************************
 * Routine: Write_Solution_Block                        *
 *                                                      *
 * Writes the cell centred solution values of the       *
 * specified quadrilateral solution block to the        *
 * specified output stream for restart purposes.        *
 *                                                      *
 ********************************************************/
void Write_Solution_Block(AdvectDiffuse2D_Quad_Block &SolnBlk,
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
void Read_Solution_Block(AdvectDiffuse2D_Quad_Block &SolnBlk,
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
void Broadcast_Solution_Block(AdvectDiffuse2D_Quad_Block &SolnBlk) {

#ifdef _MPI_VERSION
    int i, j, ni, nj, ng, block_allocated, buffer_size;
    double *buffer;

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
    } /* endif */

    MPI::COMM_WORLD.Bcast(&ni, 1, MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(&nj, 1, MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(&ng, 1, MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(&block_allocated, 1, MPI::INT, 0);

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

    /* Broadcast the solution state variables. */

    if (block_allocated) {
       ni = (SolnBlk.ICu+SolnBlk.Nghost) - (SolnBlk.ICl-SolnBlk.Nghost) + 1;
       nj = (SolnBlk.JCu+SolnBlk.Nghost) - (SolnBlk.JCl-SolnBlk.Nghost) + 1;
       buffer = new double[7*ni*nj];

       if (CFFC_Primary_MPI_Processor()) {
          buffer_size = 0;
          for (j  = SolnBlk.JCl-SolnBlk.Nghost ; j <= SolnBlk.JCu+SolnBlk.Nghost ; ++j ) {
              for ( i = SolnBlk.ICl-SolnBlk.Nghost ; i <= SolnBlk.ICu+SolnBlk.Nghost ; ++i ) {
 	          buffer[buffer_size  ] = SolnBlk.U[i][j].u;
 	          buffer[buffer_size+1] = SolnBlk.U[i][j].Fd.x;
 	          buffer[buffer_size+2] = SolnBlk.U[i][j].Fd.y;
 	          buffer[buffer_size+3] = SolnBlk.U[i][j].V.x;
 	          buffer[buffer_size+4] = SolnBlk.U[i][j].V.y;
 	          buffer[buffer_size+5] = SolnBlk.U[i][j].k;
 	          buffer[buffer_size+6] = SolnBlk.U[i][j].T;
                  buffer_size = buffer_size + 7;
              } /* endfor */
          } /* endfor */
       } /* endif */

       buffer_size = 7*ni*nj;
       MPI::COMM_WORLD.Bcast(buffer, buffer_size, MPI::DOUBLE, 0);

       if (!CFFC_Primary_MPI_Processor()) {
          buffer_size = 0;
          for (j  = SolnBlk.JCl-SolnBlk.Nghost ; j <= SolnBlk.JCu+SolnBlk.Nghost ; ++j ) {
              for ( i = SolnBlk.ICl-SolnBlk.Nghost ; i <= SolnBlk.ICu+SolnBlk.Nghost ; ++i ) {
 	          SolnBlk.U[i][j].u    = buffer[buffer_size];
 	          SolnBlk.U[i][j].Fd.x = buffer[buffer_size+1];
 	          SolnBlk.U[i][j].Fd.y = buffer[buffer_size+2];
 	          SolnBlk.U[i][j].V.x  = buffer[buffer_size+3];
 	          SolnBlk.U[i][j].V.y  = buffer[buffer_size+4];
 	          SolnBlk.U[i][j].k    = buffer[buffer_size+5];
 	          SolnBlk.U[i][j].T    = buffer[buffer_size+6];
                  buffer_size = buffer_size + 7;
              } /* endfor */
           } /* endfor */
       } /* endif */

       delete []buffer; 
       buffer = NULL;

       ni = 1;
       nj = (SolnBlk.JCu+SolnBlk.Nghost) - (SolnBlk.JCl-SolnBlk.Nghost) + 1;
       buffer = new double[14*ni*nj];

       if (CFFC_Primary_MPI_Processor()) {
          buffer_size = 0;
          for (j  = SolnBlk.JCl-SolnBlk.Nghost ; j <= SolnBlk.JCu+SolnBlk.Nghost ; ++j ) {
 	      buffer[buffer_size   ] = SolnBlk.UoW[j].u;
 	      buffer[buffer_size+ 1] = SolnBlk.UoW[j].Fd.x;
 	      buffer[buffer_size+ 2] = SolnBlk.UoW[j].Fd.y;
 	      buffer[buffer_size+ 3] = SolnBlk.UoW[j].V.x;
 	      buffer[buffer_size+ 4] = SolnBlk.UoW[j].V.y;
 	      buffer[buffer_size+ 5] = SolnBlk.UoW[j].k;
 	      buffer[buffer_size+ 6] = SolnBlk.UoW[j].T;
 	      buffer[buffer_size+ 7] = SolnBlk.UoE[j].u;
 	      buffer[buffer_size+ 8] = SolnBlk.UoE[j].Fd.x;
 	      buffer[buffer_size+ 9] = SolnBlk.UoE[j].Fd.y;
 	      buffer[buffer_size+10] = SolnBlk.UoE[j].V.x;
 	      buffer[buffer_size+11] = SolnBlk.UoE[j].V.y;
 	      buffer[buffer_size+12] = SolnBlk.UoE[j].k;
 	      buffer[buffer_size+13] = SolnBlk.UoE[j].T;
              buffer_size = buffer_size + 14;
          } /* endfor */
       } /* endif */

       buffer_size = 14*ni*nj;
       MPI::COMM_WORLD.Bcast(buffer, buffer_size, MPI::DOUBLE, 0);

       if (!CFFC_Primary_MPI_Processor()) {
          buffer_size = 0;
          for (j  = SolnBlk.JCl-SolnBlk.Nghost ; j <= SolnBlk.JCu+SolnBlk.Nghost ; ++j ) {
 	      SolnBlk.UoW[j].u    = buffer[buffer_size   ];
 	      SolnBlk.UoW[j].Fd.x = buffer[buffer_size+ 1];
 	      SolnBlk.UoW[j].Fd.y = buffer[buffer_size+ 2];
 	      SolnBlk.UoW[j].V.x  = buffer[buffer_size+ 3];
 	      SolnBlk.UoW[j].V.y  = buffer[buffer_size+ 4];
 	      SolnBlk.UoW[j].k    = buffer[buffer_size+ 5];
 	      SolnBlk.UoW[j].T    = buffer[buffer_size+ 6];
 	      SolnBlk.UoE[j].u    = buffer[buffer_size+ 7];
 	      SolnBlk.UoE[j].Fd.x = buffer[buffer_size+ 8];
 	      SolnBlk.UoE[j].Fd.y = buffer[buffer_size+ 9];
	      SolnBlk.UoE[j].V.x  = buffer[buffer_size+10];
 	      SolnBlk.UoE[j].V.y  = buffer[buffer_size+11];
 	      SolnBlk.UoE[j].k    = buffer[buffer_size+12];
 	      SolnBlk.UoE[j].T    = buffer[buffer_size+13];
              buffer_size = buffer_size + 14;
          } /* endfor */
       } /* endif */

       delete []buffer; 
       buffer = NULL;

       ni = (SolnBlk.ICu+SolnBlk.Nghost) - (SolnBlk.ICl-SolnBlk.Nghost) + 1;
       nj = 1;
       buffer = new double[14*ni*nj];

       if (CFFC_Primary_MPI_Processor()) {
          buffer_size = 0;
          for (i  = SolnBlk.ICl-SolnBlk.Nghost ; i <= SolnBlk.ICu+SolnBlk.Nghost ; ++i ) {
 	      buffer[buffer_size   ] = SolnBlk.UoS[i].u;
 	      buffer[buffer_size+ 1] = SolnBlk.UoS[i].Fd.x;
 	      buffer[buffer_size+ 2] = SolnBlk.UoS[i].Fd.y;
 	      buffer[buffer_size+ 3] = SolnBlk.UoS[i].V.x;
 	      buffer[buffer_size+ 4] = SolnBlk.UoS[i].V.y;
 	      buffer[buffer_size+ 5] = SolnBlk.UoS[i].k;
 	      buffer[buffer_size+ 6] = SolnBlk.UoS[i].T;
 	      buffer[buffer_size+ 7] = SolnBlk.UoN[i].u;
 	      buffer[buffer_size+ 8] = SolnBlk.UoN[i].Fd.x;
 	      buffer[buffer_size+ 9] = SolnBlk.UoN[i].Fd.y;
 	      buffer[buffer_size+10] = SolnBlk.UoN[i].V.x;
 	      buffer[buffer_size+11] = SolnBlk.UoN[i].V.y;
 	      buffer[buffer_size+12] = SolnBlk.UoN[i].k;
 	      buffer[buffer_size+13] = SolnBlk.UoN[i].T;
              buffer_size = buffer_size + 14;
          } /* endfor */
       } /* endif */

       buffer_size = 14*ni*nj;
       MPI::COMM_WORLD.Bcast(buffer, buffer_size, MPI::DOUBLE, 0);

       if (!CFFC_Primary_MPI_Processor()) {
          buffer_size = 0;
          for (i  = SolnBlk.ICl-SolnBlk.Nghost ; i <= SolnBlk.ICu+SolnBlk.Nghost ; ++i ) {
 	      SolnBlk.UoS[i].u    = buffer[buffer_size   ];
 	      SolnBlk.UoS[i].Fd.x = buffer[buffer_size+ 1];
 	      SolnBlk.UoS[i].Fd.y = buffer[buffer_size+ 2];
 	      SolnBlk.UoS[i].V.y  = buffer[buffer_size+ 3];
 	      SolnBlk.UoS[i].V.x  = buffer[buffer_size+ 4];
 	      SolnBlk.UoS[i].k    = buffer[buffer_size+ 5];
 	      SolnBlk.UoS[i].T    = buffer[buffer_size+ 6];
 	      SolnBlk.UoN[i].u    = buffer[buffer_size+ 7];
 	      SolnBlk.UoN[i].Fd.x = buffer[buffer_size+ 8];
 	      SolnBlk.UoN[i].Fd.y = buffer[buffer_size+ 9];
 	      SolnBlk.UoN[i].V.x  = buffer[buffer_size+10];
 	      SolnBlk.UoN[i].V.y  = buffer[buffer_size+11];
 	      SolnBlk.UoN[i].k    = buffer[buffer_size+12];
 	      SolnBlk.UoN[i].T    = buffer[buffer_size+13];
              buffer_size = buffer_size + 14;
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
void Broadcast_Solution_Block(AdvectDiffuse2D_Quad_Block &SolnBlk,
                              MPI::Intracomm &Communicator, 
                              const int Source_CPU) {

    int Source_Rank = 0;
    int i, j, ni, nj, ng, block_allocated, buffer_size;
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
    } /* endif */

    Communicator.Bcast(&ni, 1, MPI::INT, Source_Rank);
    Communicator.Bcast(&nj, 1, MPI::INT, Source_Rank);
    Communicator.Bcast(&ng, 1, MPI::INT, Source_Rank);
    Communicator.Bcast(&block_allocated, 1, MPI::INT, Source_Rank);

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

    /* Broadcast the solution state variables. */

    if (block_allocated) {
       ni = (SolnBlk.ICu+SolnBlk.Nghost) - (SolnBlk.ICl-SolnBlk.Nghost) + 1;
       nj = (SolnBlk.JCu+SolnBlk.Nghost) - (SolnBlk.JCl-SolnBlk.Nghost) + 1;
       buffer = new double[7*ni*nj];

       if (CFFC_MPI::This_Processor_Number == Source_CPU) {
          buffer_size = 0;
          for (j  = SolnBlk.JCl-SolnBlk.Nghost ; j <= SolnBlk.JCu+SolnBlk.Nghost ; ++j ) {
              for ( i = SolnBlk.ICl-SolnBlk.Nghost ; i <= SolnBlk.ICu+SolnBlk.Nghost ; ++i ) {
 	          buffer[buffer_size  ] = SolnBlk.U[i][j].u;
 	          buffer[buffer_size+1] = SolnBlk.U[i][j].Fd.x;
 	          buffer[buffer_size+2] = SolnBlk.U[i][j].Fd.y;
 	          buffer[buffer_size+3] = SolnBlk.U[i][j].V.x;
 	          buffer[buffer_size+4] = SolnBlk.U[i][j].V.y;
 	          buffer[buffer_size+5] = SolnBlk.U[i][j].k;
 	          buffer[buffer_size+6] = SolnBlk.U[i][j].T;
                  buffer_size = buffer_size + 7;
              } /* endfor */
          } /* endfor */
       } /* endif */

       buffer_size = 7*ni*nj;
       Communicator.Bcast(buffer, buffer_size, MPI::DOUBLE, Source_Rank);

       if (!(CFFC_MPI::This_Processor_Number == Source_CPU)) {
          buffer_size = 0;
          for (j  = SolnBlk.JCl-SolnBlk.Nghost ; j <= SolnBlk.JCu+SolnBlk.Nghost ; ++j ) {
              for ( i = SolnBlk.ICl-SolnBlk.Nghost ; i <= SolnBlk.ICu+SolnBlk.Nghost ; ++i ) {
                  SolnBlk.U[i][j].u    = buffer[buffer_size];
 	          SolnBlk.U[i][j].Fd.x = buffer[buffer_size+1];
 	          SolnBlk.U[i][j].Fd.y = buffer[buffer_size+2];
 	          SolnBlk.U[i][j].V.x  = buffer[buffer_size+3];
 	          SolnBlk.U[i][j].V.y  = buffer[buffer_size+4];
 	          SolnBlk.U[i][j].k    = buffer[buffer_size+5];
 	          SolnBlk.U[i][j].T    = buffer[buffer_size+6];
                  buffer_size = buffer_size + 7;
              } /* endfor */
           } /* endfor */
       } /* endif */

       delete []buffer; 
       buffer = NULL;

       ni = 1;
       nj = (SolnBlk.JCu+SolnBlk.Nghost) - (SolnBlk.JCl-SolnBlk.Nghost) + 1;
       buffer = new double[14*ni*nj];

       if (CFFC_MPI::This_Processor_Number == Source_CPU) {
          buffer_size = 0;
          for (j  = SolnBlk.JCl-SolnBlk.Nghost ; j <= SolnBlk.JCu+SolnBlk.Nghost ; ++j ) {
	      buffer[buffer_size   ] = SolnBlk.UoW[j].u;
 	      buffer[buffer_size+ 1] = SolnBlk.UoW[j].Fd.x;
 	      buffer[buffer_size+ 2] = SolnBlk.UoW[j].Fd.y;
 	      buffer[buffer_size+ 3] = SolnBlk.UoW[j].V.x;
 	      buffer[buffer_size+ 4] = SolnBlk.UoW[j].V.y;
 	      buffer[buffer_size+ 5] = SolnBlk.UoW[j].k;
 	      buffer[buffer_size+ 6] = SolnBlk.UoW[j].T;
 	      buffer[buffer_size+ 7] = SolnBlk.UoE[j].u;
 	      buffer[buffer_size+ 8] = SolnBlk.UoE[j].Fd.x;
 	      buffer[buffer_size+ 9] = SolnBlk.UoE[j].Fd.y;
 	      buffer[buffer_size+10] = SolnBlk.UoE[j].V.x;
 	      buffer[buffer_size+11] = SolnBlk.UoE[j].V.y;
 	      buffer[buffer_size+12] = SolnBlk.UoE[j].k;
 	      buffer[buffer_size+13] = SolnBlk.UoE[j].T;
              buffer_size = buffer_size + 14;
          } /* endfor */
       } /* endif */

       buffer_size = 14*ni*nj;
       Communicator.Bcast(buffer, buffer_size, MPI::DOUBLE, Source_Rank);

       if (!(CFFC_MPI::This_Processor_Number == Source_CPU)) {
          buffer_size = 0;
          for (j  = SolnBlk.JCl-SolnBlk.Nghost ; j <= SolnBlk.JCu+SolnBlk.Nghost ; ++j ) {
 	      SolnBlk.UoW[j].u    = buffer[buffer_size   ];
 	      SolnBlk.UoW[j].Fd.x = buffer[buffer_size+ 1];
 	      SolnBlk.UoW[j].Fd.y = buffer[buffer_size+ 2];
 	      SolnBlk.UoW[j].V.x  = buffer[buffer_size+ 3];
 	      SolnBlk.UoW[j].V.y  = buffer[buffer_size+ 4];
 	      SolnBlk.UoW[j].k    = buffer[buffer_size+ 5];
 	      SolnBlk.UoW[j].T    = buffer[buffer_size+ 6];
 	      SolnBlk.UoE[j].u    = buffer[buffer_size+ 7];
 	      SolnBlk.UoE[j].Fd.x = buffer[buffer_size+ 8];
 	      SolnBlk.UoE[j].Fd.y = buffer[buffer_size+ 9];
	      SolnBlk.UoE[j].V.x  = buffer[buffer_size+10];
 	      SolnBlk.UoE[j].V.y  = buffer[buffer_size+11];
 	      SolnBlk.UoE[j].k    = buffer[buffer_size+12];
 	      SolnBlk.UoE[j].T    = buffer[buffer_size+13];
              buffer_size = buffer_size + 14;
          } /* endfor */
       } /* endif */

       delete []buffer; 
       buffer = NULL;

       ni = (SolnBlk.ICu+SolnBlk.Nghost) - (SolnBlk.ICl-SolnBlk.Nghost) + 1;
       nj = 1;
       buffer = new double[14*ni*nj];

       if (CFFC_MPI::This_Processor_Number == Source_CPU) {
          buffer_size = 0;
          for (i  = SolnBlk.ICl-SolnBlk.Nghost ; i <= SolnBlk.ICu+SolnBlk.Nghost ; ++i ) {
 	      buffer[buffer_size   ] = SolnBlk.UoS[i].u;
 	      buffer[buffer_size+ 1] = SolnBlk.UoS[i].Fd.x;
 	      buffer[buffer_size+ 2] = SolnBlk.UoS[i].Fd.y;
 	      buffer[buffer_size+ 3] = SolnBlk.UoS[i].V.x;
 	      buffer[buffer_size+ 4] = SolnBlk.UoS[i].V.y;
 	      buffer[buffer_size+ 5] = SolnBlk.UoS[i].k;
 	      buffer[buffer_size+ 6] = SolnBlk.UoS[i].T;
 	      buffer[buffer_size+ 7] = SolnBlk.UoN[i].u;
 	      buffer[buffer_size+ 8] = SolnBlk.UoN[i].Fd.x;
 	      buffer[buffer_size+ 9] = SolnBlk.UoN[i].Fd.y;
 	      buffer[buffer_size+10] = SolnBlk.UoN[i].V.x;
 	      buffer[buffer_size+11] = SolnBlk.UoN[i].V.y;
 	      buffer[buffer_size+12] = SolnBlk.UoN[i].k;
 	      buffer[buffer_size+13] = SolnBlk.UoN[i].T;
              buffer_size = buffer_size + 14;
          } /* endfor */
       } /* endif */

       buffer_size = 14*ni*nj;
       Communicator.Bcast(buffer, buffer_size, MPI::DOUBLE, Source_Rank);

       if (!(CFFC_MPI::This_Processor_Number == Source_CPU)) {
          buffer_size = 0;
          for (i  = SolnBlk.ICl-SolnBlk.Nghost ; i <= SolnBlk.ICu+SolnBlk.Nghost ; ++i ) {
 	      SolnBlk.UoS[i].u    = buffer[buffer_size   ];
 	      SolnBlk.UoS[i].Fd.x = buffer[buffer_size+ 1];
 	      SolnBlk.UoS[i].Fd.y = buffer[buffer_size+ 2];
 	      SolnBlk.UoS[i].V.y  = buffer[buffer_size+ 3];
 	      SolnBlk.UoS[i].V.x  = buffer[buffer_size+ 4];
 	      SolnBlk.UoS[i].k    = buffer[buffer_size+ 5];
 	      SolnBlk.UoS[i].T    = buffer[buffer_size+ 6];
 	      SolnBlk.UoN[i].u    = buffer[buffer_size+ 7];
 	      SolnBlk.UoN[i].Fd.x = buffer[buffer_size+ 8];
 	      SolnBlk.UoN[i].Fd.y = buffer[buffer_size+ 9];
 	      SolnBlk.UoN[i].V.x  = buffer[buffer_size+10];
 	      SolnBlk.UoN[i].V.y  = buffer[buffer_size+11];
 	      SolnBlk.UoN[i].k    = buffer[buffer_size+12];
 	      SolnBlk.UoN[i].T    = buffer[buffer_size+13];
              buffer_size = buffer_size + 14;
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
void Copy_Solution_Block(AdvectDiffuse2D_Quad_Block &SolnBlk1,
                         AdvectDiffuse2D_Quad_Block &SolnBlk2) {

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
             for ( k = 0 ; k <= NUMBER_OF_RESIDUAL_VECTORS_ADVECTDIFFUSE2D-1 ; ++k ) {
	        SolnBlk1.dudt[i][j][k] = SolnBlk2.dudt[i][j][k];
             } /* endfor */
	     SolnBlk1.dudx[i][j] = SolnBlk2.dudx[i][j];
	     SolnBlk1.dudy[i][j] = SolnBlk2.dudy[i][j];
	     SolnBlk1.phi[i][j] = SolnBlk2.phi[i][j];
	     SolnBlk1.uo[i][j] = SolnBlk2.uo[i][j];
	     SolnBlk1.dt[i][j] = SolnBlk2.dt[i][j];
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
int Prolong_Solution_Block(AdvectDiffuse2D_Quad_Block &SolnBlk_Fine,
			   AdvectDiffuse2D_Quad_Block &SolnBlk_Original,
			   const int Sector) {

    int i, j, k, i_min, i_max, j_min, j_max, mesh_refinement_permitted;
    double area_total_fine;

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

       // Note the extra range (-+1) for j_min, jmax, i_min, and i_max so as to set
       // V, k, and T on fine ghost cells.
       for ( j  = j_min-1; j <= j_max+1 ; ++j ) {
	   for ( i = i_min-1 ; i <= i_max+1 ; ++i ) {
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
    	       SolnBlk_Fine.U[2*(i-i_min)+SolnBlk_Fine.ICl  ]
		             [2*(j-j_min)+SolnBlk_Fine.JCl  ].V 
                  = SolnBlk_Original.U[i][j].V;
    	       SolnBlk_Fine.U[2*(i-i_min)+SolnBlk_Fine.ICl  ]
		             [2*(j-j_min)+SolnBlk_Fine.JCl  ].k 
                  = SolnBlk_Original.U[i][j].k;
    	       SolnBlk_Fine.U[2*(i-i_min)+SolnBlk_Fine.ICl  ]
		             [2*(j-j_min)+SolnBlk_Fine.JCl  ].T 
                  = SolnBlk_Original.U[i][j].T;

    	       SolnBlk_Fine.U[2*(i-i_min)+SolnBlk_Fine.ICl+1]
                             [2*(j-j_min)+SolnBlk_Fine.JCl  ] 
                  = SolnBlk_Original.U[i][j];
 //                  = (SolnBlk_Original.Grid.Cell[i][j].A/area_total_fine)*SolnBlk_Original.U[i][j];
    	       SolnBlk_Fine.U[2*(i-i_min)+SolnBlk_Fine.ICl+1]
                             [2*(j-j_min)+SolnBlk_Fine.JCl  ].V 
                  = SolnBlk_Original.U[i][j].V;
   	       SolnBlk_Fine.U[2*(i-i_min)+SolnBlk_Fine.ICl+1]
                             [2*(j-j_min)+SolnBlk_Fine.JCl  ].k
                  = SolnBlk_Original.U[i][j].k;
   	       SolnBlk_Fine.U[2*(i-i_min)+SolnBlk_Fine.ICl+1]
                             [2*(j-j_min)+SolnBlk_Fine.JCl  ].T
                  = SolnBlk_Original.U[i][j].T;

    	       SolnBlk_Fine.U[2*(i-i_min)+SolnBlk_Fine.ICl  ]
                             [2*(j-j_min)+SolnBlk_Fine.JCl+1]
                  = SolnBlk_Original.U[i][j];
 //                  = (SolnBlk_Original.Grid.Cell[i][j].A/area_total_fine)*SolnBlk_Original.U[i][j];
    	       SolnBlk_Fine.U[2*(i-i_min)+SolnBlk_Fine.ICl  ]
                             [2*(j-j_min)+SolnBlk_Fine.JCl+1].V
                  = SolnBlk_Original.U[i][j].V;
    	       SolnBlk_Fine.U[2*(i-i_min)+SolnBlk_Fine.ICl  ]
                             [2*(j-j_min)+SolnBlk_Fine.JCl+1].k
                  = SolnBlk_Original.U[i][j].k;
    	       SolnBlk_Fine.U[2*(i-i_min)+SolnBlk_Fine.ICl  ]
                             [2*(j-j_min)+SolnBlk_Fine.JCl+1].T
                  = SolnBlk_Original.U[i][j].T;

    	       SolnBlk_Fine.U[2*(i-i_min)+SolnBlk_Fine.ICl+1]
                             [2*(j-j_min)+SolnBlk_Fine.JCl+1]
                  = SolnBlk_Original.U[i][j];
//                   = (SolnBlk_Original.Grid.Cell[i][j].A/area_total_fine)*SolnBlk_Original.U[i][j];
    	       SolnBlk_Fine.U[2*(i-i_min)+SolnBlk_Fine.ICl+1]
                             [2*(j-j_min)+SolnBlk_Fine.JCl+1].V
                  = SolnBlk_Original.U[i][j].V;
    	       SolnBlk_Fine.U[2*(i-i_min)+SolnBlk_Fine.ICl+1]
                             [2*(j-j_min)+SolnBlk_Fine.JCl+1].k
                  = SolnBlk_Original.U[i][j].k;
    	       SolnBlk_Fine.U[2*(i-i_min)+SolnBlk_Fine.ICl+1]
                             [2*(j-j_min)+SolnBlk_Fine.JCl+1].T
                  = SolnBlk_Original.U[i][j].T;
           } /* endfor */
       } /* endfor */

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
int Restrict_Solution_Block(AdvectDiffuse2D_Quad_Block &SolnBlk_Coarse,
			    AdvectDiffuse2D_Quad_Block &SolnBlk_Original_SW,
			    AdvectDiffuse2D_Quad_Block &SolnBlk_Original_SE,
			    AdvectDiffuse2D_Quad_Block &SolnBlk_Original_NW,
			    AdvectDiffuse2D_Quad_Block &SolnBlk_Original_NE) {

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
//                                                      SolnBlk_Coarse.Grid.Cell[i_coarse][j_coarse].A;
             SolnBlk_Coarse.U[i_coarse][j_coarse].V = SolnBlk_Original_SW.Un(i+1,j+1).V;
             SolnBlk_Coarse.U[i_coarse][j_coarse].k = SolnBlk_Original_SW.Un(i+1,j+1).k;
             SolnBlk_Coarse.U[i_coarse][j_coarse].T = SolnBlk_Original_SW.Un(i+1,j+1).T;
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
             SolnBlk_Coarse.U[i_coarse][j_coarse].V = SolnBlk_Original_SE.Un(i+1,j+1).V;
             SolnBlk_Coarse.U[i_coarse][j_coarse].k = SolnBlk_Original_SE.Un(i+1,j+1).k;
             SolnBlk_Coarse.U[i_coarse][j_coarse].T = SolnBlk_Original_SE.Un(i+1,j+1).T;
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
             SolnBlk_Coarse.U[i_coarse][j_coarse].V = SolnBlk_Original_NW.Un(i+1,j+1).V;
             SolnBlk_Coarse.U[i_coarse][j_coarse].k = SolnBlk_Original_NW.Un(i+1,j+1).k;
             SolnBlk_Coarse.U[i_coarse][j_coarse].T = SolnBlk_Original_NW.Un(i+1,j+1).T;
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
             SolnBlk_Coarse.U[i_coarse][j_coarse].V = SolnBlk_Original_NE.Un(i+1,j+1).V;
             SolnBlk_Coarse.U[i_coarse][j_coarse].k = SolnBlk_Original_NE.Un(i+1,j+1).k;
             SolnBlk_Coarse.U[i_coarse][j_coarse].T = SolnBlk_Original_NE.Un(i+1,j+1).T;
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
void Output_Tecplot(AdvectDiffuse2D_Quad_Block &SolnBlk,
                    const int Number_of_Time_Steps,
                    const double &Time,
                    const int Block_Number,
                    const int Output_Title,
	            ostream &Out_File) {

    int i, j;
    AdvectDiffuse2D_State U_node;

    /* Ensure boundary conditions are updated before
       evaluating solution at the nodes. */
    
    BCs(SolnBlk);

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
                << "\"Fx\" \\ \n"
                << "\"Fy\" \\ \n"
                << "\"Vx\" \n"
                << "\"Vy\" \n"
                << "\"k\" \n"
                << "\"T\" \n"
		<< "\"ue\" \n"
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
	   U_node = SolnBlk.Un(i, j);
           Out_File << " " << SolnBlk.Grid.Node[i][j].X << U_node 
		    << " " << SolnBlk.uen(i,j) << "\n";
           Out_File.unsetf(ios::scientific);
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
void Output_Cells_Tecplot(AdvectDiffuse2D_Quad_Block &SolnBlk,
                          const int Number_of_Time_Steps,
                          const double &Time,
                          const int Block_Number,
                          const int Output_Title,
	                  ostream &Out_File) {

    int i, j;

    Out_File << setprecision(14);
    if (Output_Title) {
       Out_File << "TITLE = \"" << CFFC_Name() << ": 2D Advection Diffusion Equation Solution, "
                << "Time Step/Iteration Level = " << Number_of_Time_Steps
                << ", Time = " << Time
                << "\"" << "\n"
   	        << "VARIABLES = \"x\" \\ \n"
                << "\"y\" \\ \n"
                << "\"u\" \\ \n"
                << "\"Fx\" \\ \n"
                << "\"Fy\" \\ \n"
                << "\"Vx\" \\ \n"
                << "\"Vy\" \\ \n"
                << "\"k\" \n"
                << "\"T\" \n"
		<< "\"ue\" \n"
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
           Out_File << " " 
		    << SolnBlk.Grid.Cell[i][j].Xc
                    << SolnBlk.U[i][j] << " " 
		    << SolnBlk.ue[i][j] << "\n";
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
void ICs(AdvectDiffuse2D_Quad_Block &SolnBlk,
	 const int i_ICtype,
         AdvectDiffuse2D_State *Uo) {

    int i, j, k;
    double r0;
    Vector2D x0;
    AdvectDiffuse2D_State Ul, Ur;

    /* Assign the initial data for the IVP of interest. */

    switch(i_ICtype) {
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
      default:
        // Set the solution state to the initial state Uo[0].
        for (j  = SolnBlk.JCl-SolnBlk.Nghost ; j <= SolnBlk.JCu+SolnBlk.Nghost ; ++j ) {
            for ( i = SolnBlk.ICl-SolnBlk.Nghost ; i <= SolnBlk.ICu+SolnBlk.Nghost ; ++i ) {
	       SolnBlk.U[i][j] = Uo[0];
            } /* endfor */
        } /* endfor */
        break;
    } /* endswitch */
}

/********************************************************
 * Routine: Set_Boundary_Ref_State                      *
 *                                                      *
 * Assigns specific boundary conditions and data to the *
 * solution variables of the specified quadrilateral    *
 * solution block.                                      *
 *                                                      *
 * Note:  Override Default Boundary Conditions here     *
 *                                                      *
 ********************************************************/
void Set_Boundary_Ref_State(AdvectDiffuse2D_Quad_Block &SolnBlk,
			    const int GlobalSolnBlkNum,
			    const QuadTreeBlock_DataStructure &QuadTree,
			    const AdvectDiffuse2D_Input_Parameters &Input_Parameters){
    int i, j, k;

    /* Set default values for the boundary condition reference states. */

    for (j  = SolnBlk.JCl-SolnBlk.Nghost ; j <= SolnBlk.JCu+SolnBlk.Nghost ; ++j ) {

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

    /* Override Default Boundary Conditions Here */

    /* First Ascertain whether Current Solution Block is indeed
       a Root Block.  If it is, store its indices. */
    int I = -99;
    int J = -99;
    for (i = 0; i < QuadTree.NRi; i++){
      for (j = 0; j < QuadTree.NRj; j++){
	if (QuadTree.Roots[i][j].block.gblknum == GlobalSolnBlkNum)
	  {
	    I = i;
	    J = j;
	  } /* endif */
      } /* endfor */
    } /* endfor */

    /* Set Desired Boundary Conditions if this is indeed a Root Block */
    if (I != -99 && J != -99)
      {
	 /* Determine whether problem is Circular Advection or 1D-Diffusion */

	if (Input_Parameters.Kappa != ZERO && 
	    Input_Parameters.i_Grid == GRID_SQUARE &&
	    Input_Parameters.i_Velocity_Field == VELOCITY_FIELD_ZERO) // 1D-Diffusion
	  {
	    /***** X-Direction *****/
	    if (Input_Parameters.i_ICs == IC_RIEMANN_XDIR ||
		Input_Parameters.i_ICs == IC_RIEMANN)
	      {
		/* Set West Boundary if Block is on the Western Edge of Domain
		   It is assumed that West Boundary Blocks 
		   have Root Block indices (0,J) */
		
		if (I == 0)
		  {
		    for (j = SolnBlk.JCl-SolnBlk.Nghost ; j <= SolnBlk.JCu+SolnBlk.Nghost ; ++j ) {
		      SolnBlk.UoW[j].u = ZERO;
		    } /* endfor */
		  }

		/* Set East Boundary if Block is on the Eastern Edge of Domain
		   It is assumed that East Boundary Blocks 
		   have Root Block indices (NRi-1,J) */

		if (I == QuadTree.NRi-1)
		  {
		    for (j = SolnBlk.JCl-SolnBlk.Nghost ; j <= SolnBlk.JCu+SolnBlk.Nghost ; ++j ) {
		      SolnBlk.UoE[j].u = ONE;
		    } /* endfor */
		  }
	      } /* endif */
	    
	    /***** Y-Direction *****/
	    if (Input_Parameters.i_ICs == IC_RIEMANN_YDIR)
	      {
		/* Set South Boundary if Block is on the Southern Edge of 
		   Domain.  It is assumed that South Boundary Blocks 
		   have Root Block indices (I,0) */
		
		if (J == 0)
		  {
		    for (i = SolnBlk.ICl-SolnBlk.Nghost ; i <= SolnBlk.ICu+SolnBlk.Nghost ; ++i ) {
		      SolnBlk.UoS[i].u = ZERO;
		    } /* endfor */
		  }

		/* Set North Boundary if Block is on the Northern Edge of 
		   Domain.  It is assumed that North Boundary Blocks 
		   have Root Block indices (I,NRj-1) */

		if (J == QuadTree.NRj-1)
		  {
		    for (i = SolnBlk.ICl-SolnBlk.Nghost ; i <= SolnBlk.ICu+SolnBlk.Nghost ; ++i ) {
		      SolnBlk.UoN[i].u = ONE;
		    } /* endfor */
		  } /* endif */
	      } /* endif */ 
	  }
	else if (Input_Parameters.Kappa == ZERO &&
		 Input_Parameters.i_Velocity_Field == VELOCITY_FIELD_ROTATING &&
		 Input_Parameters.i_Grid == GRID_SQUARE) // Circular Advection
	  {
	    /* Determine if BCs for current root block needs to be set. */
	    double SW_Cell_x = SolnBlk.Grid.Cell[SolnBlk.ICl][SolnBlk.JCl].Xc.x;
	    double SW_Cell_y = SolnBlk.Grid.Cell[SolnBlk.ICl][SolnBlk.JCl].Xc.y;
	    double Blk_y = SolnBlk.Grid.Cell[SolnBlk.ICl][SolnBlk.JCu].Xc.y - SW_Cell_y;
	    if (SW_Cell_x >= ZERO && SW_Cell_y >= ZERO && SW_Cell_y <= HALF*Blk_y)
	      {
		double x_edge = 0.2;
		double r;
		
		for (i  = SolnBlk.ICl-SolnBlk.Nghost ; i <= SolnBlk.ICu+SolnBlk.Nghost ; ++i ) {
		  r = sqrt(sqr(SolnBlk.Grid.Cell[i][SolnBlk.JCl].Xc.x) + sqr(SolnBlk.Grid.Cell[i][SolnBlk.JCl].Xc.y));
		  if ((r <= x_edge) || (r >= 1-x_edge))
		    {
		      SolnBlk.UoS[i].u = ZERO;
		    }
		  else if ((r > x_edge) && (r < 1-x_edge)) 
		    {
		      SolnBlk.UoS[i].u = sqr(sin(PI*(r-x_edge)/(1-2*x_edge)));
		    } /* endif */
		} /* endfor */
	      } /* endif */
	  } /* end if */
	else if (Input_Parameters.Kappa == ZERO &&
		 Input_Parameters.i_Velocity_Field == VELOCITY_FIELD_ROTATING &&
		 Input_Parameters.i_Grid == GRID_CIRCULAR_CYLINDER) // Circular Advection
	  {
	    /* Determine if BCs for current root block needs to be set. */
	    if (I == QuadTree.NRi-1)
	      {
		double x_edge = 0.2;
		double r;
		
		for (j = SolnBlk.JCl-SolnBlk.Nghost ; j <= SolnBlk.JCu+SolnBlk.Nghost ; ++j ) {
		  r = sqrt(sqr(SolnBlk.Grid.Cell[SolnBlk.ICu][j].Xc.x) + sqr(SolnBlk.Grid.Cell[SolnBlk.ICu][j].Xc.y));
		  if ((r <= x_edge) || (r >= 1-x_edge))
		    {
		      SolnBlk.UoE[j].u = ZERO;
		    }
		  else if ((r > x_edge) && (r < 1-x_edge)) 
		    {
		      SolnBlk.UoE[j].u = sqr(sin(PI*(r-x_edge)/(1-2*x_edge)));
		    } /* endif */
		} /* endfor */
	      } /* endif */
	  } /* endif */
      } /* endif */
}

/********************************************************
 * Routine: Set_Advection_Velocity_Field                *
 *                                                      *
 * Computes and assigns the advection velocity field    *
 * for the specified quadrilateral solution block.      *
 *                                                      *
 ********************************************************/
void Set_Advection_Velocity_Field(AdvectDiffuse2D_Quad_Block &SolnBlk,
	                          const int i_Velocity_Field,
                                  const double &Vx,
                                  const double &Vy) {

    int i, j;
    double omega;

    /* Assign values to the components of the advective velocity
       depending on the velocity field type. */

    switch(i_Velocity_Field) {
      case VELOCITY_FIELD_ZERO :
        // Set the advection velocity field to zero.
        for (j  = SolnBlk.JCl-SolnBlk.Nghost ; j <= SolnBlk.JCu+SolnBlk.Nghost ; ++j ) {
            for ( i = SolnBlk.ICl-SolnBlk.Nghost ; i <= SolnBlk.ICu+SolnBlk.Nghost ; ++i ) {
	       SolnBlk.U[i][j].V = Vector2D_ZERO;
            } /* endfor */
        } /* endfor */
        break;
      case VELOCITY_FIELD_UNIFORM :
        // Set the advection velocity field to uniform or constant.
        for (j  = SolnBlk.JCl-SolnBlk.Nghost ; j <= SolnBlk.JCu+SolnBlk.Nghost ; ++j ) {
            for ( i = SolnBlk.ICl-SolnBlk.Nghost ; i <= SolnBlk.ICu+SolnBlk.Nghost ; ++i ) {
	       SolnBlk.U[i][j].V = Vector2D(Vx, Vy);
            } /* endfor */
        } /* endfor */
        break;
      case VELOCITY_FIELD_ROTATING :
        // Set the advection velocity field to a pure rotation.
        omega = Vx;
        for (j  = SolnBlk.JCl-SolnBlk.Nghost ; j <= SolnBlk.JCu+SolnBlk.Nghost ; ++j ) {
            for ( i = SolnBlk.ICl-SolnBlk.Nghost ; i <= SolnBlk.ICu+SolnBlk.Nghost ; ++i ) {
	       SolnBlk.U[i][j].V.x = -omega*SolnBlk.Grid.Cell[i][j].Xc.y;
	       SolnBlk.U[i][j].V.y = omega*SolnBlk.Grid.Cell[i][j].Xc.x;
            } /* endfor */
        } /* endfor */
        break;
      default:
        // Set the advection velocity field to uniform or constant.
        for (j  = SolnBlk.JCl-SolnBlk.Nghost ; j <= SolnBlk.JCu+SolnBlk.Nghost ; ++j ) {
            for ( i = SolnBlk.ICl-SolnBlk.Nghost ; i <= SolnBlk.ICu+SolnBlk.Nghost ; ++i ) {
	       SolnBlk.U[i][j].V = Vector2D(Vx, Vy);
            } /* endfor */
        } /* endfor */
        break;
    } /* endswitch */

    /* Set the advective velocity for the for the boundary conditions reference states. */

    for (j  = SolnBlk.JCl-SolnBlk.Nghost ; j <= SolnBlk.JCu+SolnBlk.Nghost ; ++j ) {
       SolnBlk.UoW[j].V = SolnBlk.U[SolnBlk.ICl-1][j].V;
       SolnBlk.UoE[j].V = SolnBlk.U[SolnBlk.ICu+1][j].V;  
    } /* endfor */

    for ( i = SolnBlk.ICl-SolnBlk.Nghost ; i <= SolnBlk.ICu+SolnBlk.Nghost ; ++i ) {
       SolnBlk.UoS[i].V = SolnBlk.U[i][SolnBlk.JCl-1].V;
       SolnBlk.UoN[i].V = SolnBlk.U[i][SolnBlk.JCu+1].V;
    } /* endfor */

}

/********************************************************
 * Routine: Set_Analytical_Solution                     *
 *                                                      *
 * Computes and assigns the analytical solution         *
 * for the specified quadrilateral solution block.      *
 *                                                      *
 ********************************************************/
void Set_Analytical_Solution(AdvectDiffuse2D_Quad_Block &SolnBlk,
			     const AdvectDiffuse2D_Input_Parameters &Input_Parameters) {
    int i, j;

    /* Determine whether problem is Circular Advection or 1D-Diffusion */

    if (Input_Parameters.Kappa != ZERO && 
	Input_Parameters.i_Grid == GRID_SQUARE &&
	Input_Parameters.i_Velocity_Field == VELOCITY_FIELD_ZERO) // 1D-Diffusion
      {
	/* Calculate lambda = -1/(kappa*tau) */
	double lambda = -1/(SolnBlk.U[SolnBlk.ICl][SolnBlk.JCl].k*SolnBlk.U[SolnBlk.ICl][SolnBlk.JCl].T);

	/* Loop through every cell */
	for (j  = SolnBlk.JCl-SolnBlk.Nghost ; j <= SolnBlk.JCu+SolnBlk.Nghost ; ++j ) {
	  for ( i = SolnBlk.ICl-SolnBlk.Nghost ; i <= SolnBlk.ICu+SolnBlk.Nghost ; ++i ) {
	    /***** X-Direction *****/
	    if (Input_Parameters.i_ICs == IC_RIEMANN_XDIR ||
		Input_Parameters.i_ICs == IC_RIEMANN)
	      {
		/* BCs are of the form 
		   u(xa) = ua
		   u(xb) = ub */
		double xa = 1.0;
		double ua = 1.0;
		double xb = -1.0;
		double ub = 0.0;
		
		if (lambda > ZERO)
		  {
		    double sqld,detA,c1,c2;
		    sqld = sqrt(lambda);
		    detA = cos(sqld*xa)*sin(sqld*xb)-cos(sqld*xb)*sin(sqld*xa);
		    c1 = ( ua*sin(sqld*xb)-ub*sin(sqld*xa))/detA;
		    c2 = (-ua*cos(sqld*xb)+ub*cos(sqld*xa))/detA;
		    SolnBlk.ue[i][j] = c1*cos(sqld*SolnBlk.Grid.Cell[i][j].Xc.x)+c2*sin(sqld*SolnBlk.Grid.Cell[i][j].Xc.x);
		  }
		else if (lambda < ZERO)
		  {
		    double sqs,detA,c3,c4;
		    sqs = sqrt(-lambda);
		    detA = cosh(sqs*xa)*sinh(sqs*xb)-cosh(sqs*xb)*sinh(sqs*xa);
		    c3 = ( ua*sinh(sqs*xb)-ub*sinh(sqs*xa))/detA;
		    c4 = (-ua*cosh(sqs*xb)+ub*cosh(sqs*xa))/detA;
		    SolnBlk.ue[i][j] = c3*cosh(sqs*SolnBlk.Grid.Cell[i][j].Xc.x)+c4*sinh(sqs*SolnBlk.Grid.Cell[i][j].Xc.x);
		  }
		else /* (lambda == ZERO) */
		  {
		    double detA,c1,c2;
		    detA = xb-xa;
		    c1 = (ua*xb-ub*xa)/detA;
		    c2 = (-ua+ub)/detA;
		    SolnBlk.ue[i][j] = c1+c2*SolnBlk.Grid.Cell[i][j].Xc.x;
		  } /* endif */
	      } /* endif */
	    
	    /***** Y-Direction *****/
	    if (Input_Parameters.i_ICs == IC_RIEMANN_YDIR)
	      {
		/* BCs are of the form 
		   u(ya) = ua
		   u(yb) = ub */
		double ya = 1.0;
		double ua = 1.0;
		double yb = -1.0;
		double ub = 0.0;
		
		if (lambda > ZERO)
		  {
		    double sqld,detA,c1,c2;
		    sqld = sqrt(lambda);
		    detA = cos(sqld*ya)*sin(sqld*yb)-cos(sqld*yb)*sin(sqld*ya);
		    c1 = (ua*sin(sqld*yb)-ub*sin(sqld*ya))/detA;
		    c2 = (-ua*cos(sqld*yb)+ub*cos(sqld*ya))/detA;
		    SolnBlk.ue[i][j] = c1*cos(sqld*SolnBlk.Grid.Cell[i][j].Xc.y)+c2*sin(sqld*SolnBlk.Grid.Cell[i][j].Xc.y);
		  }
		else if (lambda < ZERO)
		  {
		    double sqs,detA,c3,c4;
		    sqs = sqrt(-lambda);
		    detA = cosh(sqs*ya)*sinh(sqs*yb)-cosh(sqs*yb)*sinh(sqs*ya);
		    c3 = (ua*sinh(sqs*yb)-ub*sinh(sqs*ya))/detA;
		    c4 = (-ua*cosh(sqs*yb)+ub*cosh(sqs*ya))/detA;
		    SolnBlk.ue[i][j] = c3*cosh(sqs*SolnBlk.Grid.Cell[i][j].Xc.y)+c4*sinh(sqs*SolnBlk.Grid.Cell[i][j].Xc.y);
		  }
		else /* (lambda == ZERO) */
		  {
		    double detA,c1,c2;
		    detA = yb-ya;
		    c1 = (ua*yb-ub*ya)/detA;
		    c2 = (-ua+ub)/detA;
		    SolnBlk.ue[i][j] = c1+c2*SolnBlk.Grid.Cell[i][j].Xc.y;
		  } /* endif */
	      } /* endif */
	  } /* endfor */
	} /* endfor */
      }
    else if (Input_Parameters.Kappa == ZERO &&
	     Input_Parameters.i_Velocity_Field == VELOCITY_FIELD_ROTATING) // Circular Advection
      {
	/* This variable controls how narrow the sin^2 IC is, and must be set
	   to the same as in the function "Set_Boundary_Ref_State" */
	double x_edge = 0.2;
	double r;
	for (j = SolnBlk.JCl-SolnBlk.Nghost ; j <= SolnBlk.JCu+SolnBlk.Nghost ; ++j ) {
            for ( i = SolnBlk.ICl-SolnBlk.Nghost ; i <= SolnBlk.ICu+SolnBlk.Nghost ; ++i ) {
	       r = sqrt(sqr(SolnBlk.Grid.Cell[i][j].Xc.x) + sqr(SolnBlk.Grid.Cell[i][j].Xc.y));
	       if ((r <= x_edge) || (r >= 1-x_edge))
		 {
		   SolnBlk.ue[i][j] = ZERO;
		 }
	       else if ((r > x_edge) && (r < 1-x_edge))
		 {
		   SolnBlk.ue[i][j] = sqr(sin(PI*(r-x_edge)/(1-2*x_edge)));
		 } /* endif */
            } /* endfor */
        } /* endfor */
      }
    else
      {
	/* No Analytical Solution - Set **ue to zero */
	for (j = SolnBlk.JCl-SolnBlk.Nghost ; j <= SolnBlk.JCu+SolnBlk.Nghost ; ++j ) {
            for ( i = SolnBlk.ICl-SolnBlk.Nghost ; i <= SolnBlk.ICu+SolnBlk.Nghost ; ++i ) {
	      SolnBlk.ue[i][j] = ZERO;
	    } /* endfor */
        } /* endfor */
      } /* endif */
}

/********************************************************
 * Routine: BCs                                         *
 *                                                      *
 * Apply boundary conditions at boundaries of the       *
 * specified quadrilateral solution block.              *
 *                                                      *
 ********************************************************/
void BCs(AdvectDiffuse2D_Quad_Block &SolnBlk) {

    int i, j;
    double dx_norm, du, dudx;
    Vector2D dX;

    for ( j = SolnBlk.JCl-SolnBlk.Nghost ; j <= SolnBlk.JCu+SolnBlk.Nghost ; ++j ) {
      if ( (j >= SolnBlk.JCl && j <= SolnBlk.JCu) ||
           (j < SolnBlk.JCl && 
            (SolnBlk.Grid.BCtypeS[SolnBlk.ICl-1] == BC_NONE ||
             SolnBlk.Grid.BCtypeS[SolnBlk.ICl-1] == BC_PERIODIC ||
             SolnBlk.Grid.BCtypeS[SolnBlk.ICl-1] == BC_DIRICHLET ||
             SolnBlk.Grid.BCtypeS[SolnBlk.ICl-1] == BC_NEUMANN ||
             SolnBlk.Grid.BCtypeS[SolnBlk.ICl-1] == BC_ROBIN) ) ||
           (j > SolnBlk.JCu && 
            (SolnBlk.Grid.BCtypeN[SolnBlk.ICl-1] == BC_NONE ||
             SolnBlk.Grid.BCtypeN[SolnBlk.ICl-1] == BC_PERIODIC ||
             SolnBlk.Grid.BCtypeN[SolnBlk.ICl-1] == BC_DIRICHLET ||
             SolnBlk.Grid.BCtypeN[SolnBlk.ICl-1] == BC_NEUMANN ||
             SolnBlk.Grid.BCtypeN[SolnBlk.ICl-1] == BC_ROBIN) ) ) {
        switch(SolnBlk.Grid.BCtypeW[j]) {
          case BC_NONE :
            break;
          case BC_DIRICHLET :
            SolnBlk.U[SolnBlk.ICl-1][j] = SolnBlk.UoW[j];
            SolnBlk.U[SolnBlk.ICl-2][j] = SolnBlk.UoW[j];
            if (SolnBlk.UoW[j].k > TOLER) {
               dX = SolnBlk.Grid.xfaceW(SolnBlk.ICl, j)-
                    SolnBlk.Grid.Cell[SolnBlk.ICl][j].Xc;
               dx_norm = dX*SolnBlk.Grid.nfaceW(SolnBlk.ICl, j);
               du = SolnBlk.UoW[j].u-SolnBlk.U[SolnBlk.ICl][j].u;
               dudx = du/dx_norm;
               dX = SolnBlk.Grid.Cell[SolnBlk.ICl-1][j].Xc -
                    SolnBlk.Grid.Cell[SolnBlk.ICl][j].Xc;
               dx_norm = dX*SolnBlk.Grid.nfaceW(SolnBlk.ICl, j);
               SolnBlk.U[SolnBlk.ICl-1][j].u = SolnBlk.U[SolnBlk.ICl][j].u + dudx*dx_norm;
               dX = SolnBlk.Grid.Cell[SolnBlk.ICl-2][j].Xc -
                    SolnBlk.Grid.Cell[SolnBlk.ICl][j].Xc;
               dx_norm = dX*SolnBlk.Grid.nfaceW(SolnBlk.ICl, j);
               SolnBlk.U[SolnBlk.ICl-2][j].u = SolnBlk.U[SolnBlk.ICl][j].u + dudx*dx_norm;
            } /* endif */
            break;
          case BC_NEUMANN :
            SolnBlk.U[SolnBlk.ICl-1][j] = SolnBlk.U[SolnBlk.ICl][j];
            SolnBlk.U[SolnBlk.ICl-2][j] = SolnBlk.U[SolnBlk.ICl][j];
            break;
          case BC_ROBIN :

            break;
          case BC_PERIODIC :
            SolnBlk.U[SolnBlk.ICl-1][j] = SolnBlk.U[SolnBlk.ICu-1][j];
            SolnBlk.U[SolnBlk.ICl-2][j] = SolnBlk.U[SolnBlk.ICu-2][j];
            break;
          default:
            SolnBlk.U[SolnBlk.ICl-1][j] = SolnBlk.UoW[j];
            SolnBlk.U[SolnBlk.ICl-2][j] = SolnBlk.UoW[j];
            break;
        } /* endswitch */
      } /* endif */

      if ( (j >= SolnBlk.JCl && j <= SolnBlk.JCu) ||
           (j < SolnBlk.JCl && 
            (SolnBlk.Grid.BCtypeS[SolnBlk.ICu+1] == BC_NONE ||
             SolnBlk.Grid.BCtypeS[SolnBlk.ICu+1] == BC_PERIODIC ||
             SolnBlk.Grid.BCtypeS[SolnBlk.ICu+1] == BC_DIRICHLET ||
             SolnBlk.Grid.BCtypeS[SolnBlk.ICu+1] == BC_NEUMANN ||
             SolnBlk.Grid.BCtypeS[SolnBlk.ICu+1] == BC_ROBIN) ) ||
           (j > SolnBlk.JCu && 
            (SolnBlk.Grid.BCtypeN[SolnBlk.ICu+1] == BC_NONE ||
             SolnBlk.Grid.BCtypeN[SolnBlk.ICu+1] == BC_PERIODIC ||
             SolnBlk.Grid.BCtypeN[SolnBlk.ICu+1] == BC_DIRICHLET ||
             SolnBlk.Grid.BCtypeN[SolnBlk.ICu+1] == BC_NEUMANN ||
             SolnBlk.Grid.BCtypeN[SolnBlk.ICu+1] == BC_ROBIN) ) ) {
        switch(SolnBlk.Grid.BCtypeE[j]) {
          case BC_NONE :
            break;
          case BC_DIRICHLET :
            SolnBlk.U[SolnBlk.ICu+1][j] = SolnBlk.UoE[j];
            SolnBlk.U[SolnBlk.ICu+2][j] = SolnBlk.UoE[j];
            if (SolnBlk.UoE[j].k > TOLER) {
               dX = SolnBlk.Grid.xfaceE(SolnBlk.ICu, j)-
                    SolnBlk.Grid.Cell[SolnBlk.ICu][j].Xc;
               dx_norm = dX*SolnBlk.Grid.nfaceE(SolnBlk.ICu, j);
               du = SolnBlk.UoE[j].u-SolnBlk.U[SolnBlk.ICu][j].u;
               dudx = du/dx_norm;
               dX = SolnBlk.Grid.Cell[SolnBlk.ICu+1][j].Xc -
                    SolnBlk.Grid.Cell[SolnBlk.ICu][j].Xc;
               dx_norm = dX*SolnBlk.Grid.nfaceE(SolnBlk.ICu, j);
               SolnBlk.U[SolnBlk.ICu+1][j].u = SolnBlk.U[SolnBlk.ICu][j].u + dudx*dx_norm;
               dX = SolnBlk.Grid.Cell[SolnBlk.ICu+2][j].Xc -
                    SolnBlk.Grid.Cell[SolnBlk.ICu][j].Xc;
               dx_norm = dX*SolnBlk.Grid.nfaceE(SolnBlk.ICu, j);
               SolnBlk.U[SolnBlk.ICu+2][j].u = SolnBlk.U[SolnBlk.ICu][j].u + dudx*dx_norm;
            } /* endif */
            break;
          case BC_NEUMANN :
            SolnBlk.U[SolnBlk.ICu+1][j] = SolnBlk.U[SolnBlk.ICu][j];
            SolnBlk.U[SolnBlk.ICu+2][j] = SolnBlk.U[SolnBlk.ICu][j];
            break;
          case BC_ROBIN :

            break;
          case BC_PERIODIC :
            SolnBlk.U[SolnBlk.ICu+1][j] = SolnBlk.U[SolnBlk.ICl+1][j];
            SolnBlk.U[SolnBlk.ICu+2][j] = SolnBlk.U[SolnBlk.ICl+2][j];
            break;
          default:
            SolnBlk.U[SolnBlk.ICu+1][j] = SolnBlk.UoE[j];
            SolnBlk.U[SolnBlk.ICu+2][j] = SolnBlk.UoE[j];
            break;
        } /* endswitch */
      } /* endif */
    } /* endfor */

    for ( i = SolnBlk.ICl-SolnBlk.Nghost ; i <= SolnBlk.ICu+SolnBlk.Nghost ; ++i ) {
      if ( (i >= SolnBlk.ICl && i <= SolnBlk.ICu) ||
           (i < SolnBlk.ICl && 
            (SolnBlk.Grid.BCtypeW[SolnBlk.JCl-1] == BC_NONE ||
             SolnBlk.Grid.BCtypeW[SolnBlk.JCl-1] == BC_PERIODIC ||
             SolnBlk.Grid.BCtypeW[SolnBlk.JCl-1] == BC_NEUMANN ||
             SolnBlk.Grid.BCtypeW[SolnBlk.JCl-1] == BC_DIRICHLET||
             SolnBlk.Grid.BCtypeW[SolnBlk.JCl-1] == BC_ROBIN) ) ||
           (i > SolnBlk.ICu && 
            (SolnBlk.Grid.BCtypeE[SolnBlk.JCl-1] == BC_NONE ||
             SolnBlk.Grid.BCtypeE[SolnBlk.JCl-1] == BC_PERIODIC ||
             SolnBlk.Grid.BCtypeE[SolnBlk.JCl-1] == BC_DIRICHLET ||
             SolnBlk.Grid.BCtypeE[SolnBlk.JCl-1] == BC_NEUMANN ||
             SolnBlk.Grid.BCtypeE[SolnBlk.JCl-1] == BC_ROBIN) ) ) {
        switch(SolnBlk.Grid.BCtypeS[i]) {
          case BC_NONE :
            break;
          case BC_DIRICHLET :
            SolnBlk.U[i][SolnBlk.JCl-1] = SolnBlk.UoS[i];
            SolnBlk.U[i][SolnBlk.JCl-2] = SolnBlk.UoS[i];
            if (SolnBlk.UoS[i].k > TOLER) {
               dX = SolnBlk.Grid.xfaceS(i, SolnBlk.JCl)-
                    SolnBlk.Grid.Cell[i][SolnBlk.JCl].Xc;
               dx_norm = dX*SolnBlk.Grid.nfaceS(i, SolnBlk.JCl);
               du = SolnBlk.UoS[i].u-SolnBlk.U[i][SolnBlk.JCl].u;
               dudx = du/dx_norm;
               dX = SolnBlk.Grid.Cell[i][SolnBlk.JCl-1].Xc -
                    SolnBlk.Grid.Cell[i][SolnBlk.JCl].Xc;
               dx_norm = dX*SolnBlk.Grid.nfaceS(i, SolnBlk.JCl);
               SolnBlk.U[i][SolnBlk.JCl-1].u = SolnBlk.U[i][SolnBlk.JCl].u + dudx*dx_norm;
               dX = SolnBlk.Grid.Cell[i][SolnBlk.JCl-2].Xc -
                    SolnBlk.Grid.Cell[i][SolnBlk.JCl].Xc;
               dx_norm = dX*SolnBlk.Grid.nfaceS(i, SolnBlk.JCl);
               SolnBlk.U[i][SolnBlk.JCl-2].u = SolnBlk.U[i][SolnBlk.JCl].u + dudx*dx_norm;
            } /* endif */
            break;
          case BC_NEUMANN :
            SolnBlk.U[i][SolnBlk.JCl-1] = SolnBlk.U[i][SolnBlk.JCl];
            SolnBlk.U[i][SolnBlk.JCl-2] = SolnBlk.U[i][SolnBlk.JCl];
            break;
          case BC_ROBIN :

            break;
          case BC_PERIODIC :
            SolnBlk.U[i][SolnBlk.JCl-1] = SolnBlk.U[i][SolnBlk.JCu-1];
            SolnBlk.U[i][SolnBlk.JCl-2] = SolnBlk.U[i][SolnBlk.JCu-2];
            break;
          default:
            SolnBlk.U[i][SolnBlk.JCl-1] = SolnBlk.UoS[i];
            SolnBlk.U[i][SolnBlk.JCl-2] = SolnBlk.UoS[i];
            break;
        } /* endswitch */
      } /* endif */

      if ( (i >= SolnBlk.ICl && i <= SolnBlk.ICu) ||
           (i < SolnBlk.ICl && 
            (SolnBlk.Grid.BCtypeW[SolnBlk.JCu+1] == BC_NONE ||
             SolnBlk.Grid.BCtypeW[SolnBlk.JCu+1] == BC_PERIODIC ||
             SolnBlk.Grid.BCtypeW[SolnBlk.JCu+1] == BC_NEUMANN ||
             SolnBlk.Grid.BCtypeW[SolnBlk.JCu+1] == BC_DIRICHLET||
             SolnBlk.Grid.BCtypeW[SolnBlk.JCu+1] == BC_ROBIN) ) ||
           (i > SolnBlk.ICu && 
            (SolnBlk.Grid.BCtypeE[SolnBlk.JCu+1] == BC_NONE ||
             SolnBlk.Grid.BCtypeE[SolnBlk.JCu+1] == BC_PERIODIC ||
             SolnBlk.Grid.BCtypeE[SolnBlk.JCu+1] == BC_DIRICHLET ||
             SolnBlk.Grid.BCtypeE[SolnBlk.JCu+1] == BC_NEUMANN ||
             SolnBlk.Grid.BCtypeE[SolnBlk.JCu+1] == BC_ROBIN) ) ) {
        switch(SolnBlk.Grid.BCtypeN[i]) {
          case BC_NONE :
            break;
          case BC_DIRICHLET :
            SolnBlk.U[i][SolnBlk.JCu+1] = SolnBlk.UoN[i];
            SolnBlk.U[i][SolnBlk.JCu+2] = SolnBlk.UoN[i];
            if (SolnBlk.UoN[i].k > TOLER) {
               dX = SolnBlk.Grid.xfaceN(i, SolnBlk.JCu)-
                    SolnBlk.Grid.Cell[i][SolnBlk.JCu].Xc;
               dx_norm = dX*SolnBlk.Grid.nfaceN(i, SolnBlk.JCu);
               du = SolnBlk.UoN[i].u-SolnBlk.U[i][SolnBlk.JCu].u;
               dudx = du/dx_norm;
               dX = SolnBlk.Grid.Cell[i][SolnBlk.JCu+1].Xc -
                    SolnBlk.Grid.Cell[i][SolnBlk.JCu].Xc;
               dx_norm = dX*SolnBlk.Grid.nfaceN(i, SolnBlk.JCu);
               SolnBlk.U[i][SolnBlk.JCu+1].u = SolnBlk.U[i][SolnBlk.JCu].u + dudx*dx_norm;
               dX = SolnBlk.Grid.Cell[i][SolnBlk.JCu+2].Xc -
                    SolnBlk.Grid.Cell[i][SolnBlk.JCu].Xc;
               dx_norm = dX*SolnBlk.Grid.nfaceN(i, SolnBlk.JCu);
               SolnBlk.U[i][SolnBlk.JCu+2].u = SolnBlk.U[i][SolnBlk.JCu].u + dudx*dx_norm;
            } /* endif */
            break;
          case BC_NEUMANN :
            SolnBlk.U[i][SolnBlk.JCu+1] = SolnBlk.U[i][SolnBlk.JCu];
            SolnBlk.U[i][SolnBlk.JCu+2] = SolnBlk.U[i][SolnBlk.JCu];
            break;
          case BC_ROBIN :
            break;
          case BC_PERIODIC :
            SolnBlk.U[i][SolnBlk.JCu+1] = SolnBlk.U[i][SolnBlk.JCl+1];
            SolnBlk.U[i][SolnBlk.JCu+2] = SolnBlk.U[i][SolnBlk.JCl+2];
            break;
          default:
            SolnBlk.U[i][SolnBlk.JCu+1] = SolnBlk.UoN[i];
            SolnBlk.U[i][SolnBlk.JCu+2] = SolnBlk.UoN[i];
            break;
        } /* endswitch */
      } /* endif */
    } /* endfor */

}

/********************************************************
 * Routine: CFL                                         *
 *                                                      *
 * Determines the allowable global and local time steps *
 * (for explicit Euler time stepping scheme) for the    *
 * specified quadrilateral solution block according to  *
 * the Courant-Friedrichs-Lewy condition for the        *
 * advection terms and a semi-impirical criteria for    *
 * the diffusion and source terms.                      *
 *                                                      *
 ********************************************************/
double CFL(AdvectDiffuse2D_Quad_Block &SolnBlk,
           AdvectDiffuse2D_Input_Parameters &Input_Parameters) {

    int i, j;
    double dtMin, d_i, d_j, v_i, v_j, a, dt_cfl, dt_diff, dt_src;

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

	     v_i = HALF*(SolnBlk.U[i][j].V*
                   (SolnBlk.Grid.nfaceE(i, j)-SolnBlk.Grid.nfaceW(i, j)));
	     v_j = HALF*(SolnBlk.U[i][j].V*
                   (SolnBlk.Grid.nfaceN(i, j)-SolnBlk.Grid.nfaceS(i, j)));

             if (fabs(v_i) > TOLER) {
	       dt_cfl = d_i/fabs(v_i);
             } else {
               dt_cfl = MILLION;
             } /* endif */
	     if (fabs(v_j) > TOLER) {
                dt_cfl = min(dt_cfl, d_j/fabs(v_j));
             } /* endif */

	     if (SolnBlk.U[i][j].k > TOLER) {
                dt_diff = HALF*min(sqr(d_i), sqr(d_j))/SolnBlk.U[i][j].k;
             } else {
	        dt_diff = MILLION;
             } /* endif */

             dt_src = HALF*SolnBlk.U[i][j].T;

	     SolnBlk.dt[i][j] = min(min(dt_cfl, dt_diff), dt_src);

             dtMin = min(dtMin, SolnBlk.dt[i][j]);
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
void Set_Global_TimeStep(AdvectDiffuse2D_Quad_Block &SolnBlk,
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
double L1_Norm_Residual(AdvectDiffuse2D_Quad_Block &SolnBlk) {

    int i, j;
    double l1_norm;

    l1_norm = ZERO;

    for ( j  = SolnBlk.JCl ; j <= SolnBlk.JCu ; ++j ) {
       for ( i = SolnBlk.ICl ; i <= SolnBlk.ICu ; ++i ) {
          l1_norm += fabs(SolnBlk.dudt[i][j][0]);
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
double L2_Norm_Residual(AdvectDiffuse2D_Quad_Block &SolnBlk) {

    int i, j;
    double l2_norm;

    l2_norm = ZERO;

    for ( j  = SolnBlk.JCl ; j <= SolnBlk.JCu ; ++j ) {
       for ( i = SolnBlk.ICl ; i <= SolnBlk.ICu ; ++i ) {
          l2_norm += sqr(SolnBlk.dudt[i][j][0]);
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
double Max_Norm_Residual(AdvectDiffuse2D_Quad_Block &SolnBlk) {

    int i, j;
    double max_norm;

    max_norm = ZERO;

    for ( j  = SolnBlk.JCl ; j <= SolnBlk.JCu ; ++j ) {
       for ( i = SolnBlk.ICl ; i <= SolnBlk.ICu ; ++i ) {
          max_norm = max(max_norm, fabs(SolnBlk.dudt[i][j][0]));
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
    if (i == SolnBlk.ICl-SolnBlk.Nghost || i == SolnBlk.ICu+SolnBlk.Nghost ||
        j == SolnBlk.JCl-SolnBlk.Nghost || j == SolnBlk.JCu+SolnBlk.Nghost) {
      n_pts = 0;
    } else if ((i == SolnBlk.ICl-SolnBlk.Nghost+1) && 
               (SolnBlk.Grid.BCtypeW[j] != BC_NONE)) {
      if (j == SolnBlk.JCl-SolnBlk.Nghost+1 || j == SolnBlk.JCu+SolnBlk.Nghost-1) {
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
    } else if ((i == SolnBlk.ICu+SolnBlk.Nghost-1) && 
               (SolnBlk.Grid.BCtypeE[j] != BC_NONE)) {
      if (j == SolnBlk.JCl-SolnBlk.Nghost+1 || j == SolnBlk.JCu+SolnBlk.Nghost-1) {
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
    } else if ((j == SolnBlk.JCl-SolnBlk.Nghost+1) && 
               (SolnBlk.Grid.BCtypeS[i] != BC_NONE)) {
      if (i == SolnBlk.ICl-SolnBlk.Nghost+1 || i == SolnBlk.ICu+SolnBlk.Nghost-1) {
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
    } else if ((j == SolnBlk.JCu+SolnBlk.Nghost-1) && 
               (SolnBlk.Grid.BCtypeN[i] != BC_NONE)) {
      if (i == SolnBlk.ICl-SolnBlk.Nghost+1 || i == SolnBlk.ICu+SolnBlk.Nghost-1) {
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
//     } else if ((i == SolnBlk.ICl && j == SolnBlk.JCl) && 
//                (SolnBlk.Grid.BCtypeW[j] != BC_NONE &&
//                 SolnBlk.Grid.BCtypeS[i] != BC_NONE)) {
//       n_pts = 7;
//       i_index[0] = i  ; j_index[0] = j-1;
//       i_index[1] = i+1; j_index[1] = j-1;
//       i_index[2] = i-1; j_index[2] = j  ;
//       i_index[3] = i+1; j_index[3] = j  ;
//       i_index[4] = i-1; j_index[4] = j+1;
//       i_index[5] = i  ; j_index[5] = j+1;
//       i_index[6] = i+1; j_index[6] = j+1;
//     } else if ((i == SolnBlk.ICl && j == SolnBlk.JCu) && 
//                (SolnBlk.Grid.BCtypeW[j] != BC_NONE &&
//                 SolnBlk.Grid.BCtypeN[i] != BC_NONE)) {
//       n_pts = 7;
//       i_index[0] = i-1; j_index[0] = j-1;
//       i_index[1] = i  ; j_index[1] = j-1;
//       i_index[2] = i+1; j_index[2] = j-1;
//       i_index[3] = i-1; j_index[3] = j  ;
//       i_index[4] = i+1; j_index[4] = j  ;
//       i_index[5] = i  ; j_index[5] = j+1;
//       i_index[6] = i+1; j_index[6] = j+1;
//     } else if ((i == SolnBlk.ICu && j == SolnBlk.JCl) && 
//                (SolnBlk.Grid.BCtypeE[j] != BC_NONE &&
//                 SolnBlk.Grid.BCtypeS[i] != BC_NONE)) {
//       n_pts = 7;
//       i_index[0] = i-1; j_index[0] = j-1;
//       i_index[1] = i  ; j_index[1] = j-1;
//       i_index[2] = i-1; j_index[2] = j  ;
//       i_index[3] = i+1; j_index[3] = j  ;
//       i_index[4] = i-1; j_index[4] = j+1;
//       i_index[5] = i  ; j_index[5] = j+1;
//       i_index[6] = i+1; j_index[6] = j+1;
//     } else if ((i == SolnBlk.ICu && j == SolnBlk.JCu) && 
//                (SolnBlk.Grid.BCtypeE[j] != BC_NONE &&
//                 SolnBlk.Grid.BCtypeN[i] != BC_NONE)) {
//       n_pts = 7;
//       i_index[0] = i-1; j_index[0] = j-1;
//       i_index[1] = i  ; j_index[1] = j-1;
//       i_index[2] = i+1; j_index[2] = j-1;
//       i_index[3] = i-1; j_index[3] = j  ;
//       i_index[4] = i+1; j_index[4] = j  ;
//       i_index[5] = i-1; j_index[5] = j+1;
//       i_index[6] = i  ; j_index[6] = j+1;
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
           SolnBlk.dudx[i][j] = u_face*n_north.x;
	   SolnBlk.dudy[i][j] = u_face*n_north.y;

           u_face = HALF*(u_sw+u_se)*l_south; 
           SolnBlk.dudx[i][j] += u_face*n_south.x;
	   SolnBlk.dudy[i][j] += u_face*n_south.y;

           u_face = HALF*(u_ne+u_se)*l_east; 
           SolnBlk.dudx[i][j] += u_face*n_east.x;
	   SolnBlk.dudy[i][j] += u_face*n_east.y;

           u_face = HALF*(u_nw+u_sw)*l_west; 
           SolnBlk.dudx[i][j] += u_face*n_west.x;
	   SolnBlk.dudy[i][j] += u_face*n_west.y;

           SolnBlk.dudx[i][j] = SolnBlk.dudx[i][j]/
                                SolnBlk.Grid.Cell[i][j].A;
           SolnBlk.dudy[i][j] = SolnBlk.dudy[i][j]/
                                SolnBlk.Grid.Cell[i][j].A;

        // If <8 neighbours are used, apply least-squares reconstruction
        } else {
           DuDx_ave = ZERO;
           DuDy_ave = ZERO;
           DxDx_ave = ZERO;
           DxDy_ave = ZERO;
           DyDy_ave = ZERO;
    
           for ( n = 0 ; n <= n_pts-1 ; ++n ) {
               dX = SolnBlk.Grid.Cell[ i_index[n] ][ j_index[n] ].Xc - 
                    SolnBlk.Grid.Cell[i][j].Xc;
               Du = SolnBlk.U[ i_index[n] ][ j_index[n] ].u - 
                    SolnBlk.U[i][j].u;
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
           SolnBlk.dudx[i][j] = (DuDx_ave*DyDy_ave-DuDy_ave*DxDy_ave)/
                                (DxDx_ave*DyDy_ave-DxDy_ave*DxDy_ave);
           SolnBlk.dudy[i][j] = (DuDy_ave*DxDx_ave-DuDx_ave*DxDy_ave)/
                                (DxDx_ave*DyDy_ave-DxDy_ave*DxDy_ave);
        } /* endif */
    
        // Calculate slope limiter.
	if (!SolnBlk.Freeze_Limiter) {
           u0Min = SolnBlk.U[i][j].u;
           u0Max = u0Min;
           for ( n = 0 ; n <= n_pts-1 ; ++n ) {
              u0Min = min(u0Min, SolnBlk.U[ i_index[n] ][ j_index[n] ].u);
              u0Max = max(u0Max, SolnBlk.U[ i_index[n] ][ j_index[n] ].u);
           } /* endfor */
    
           dX = SolnBlk.Grid.xfaceE(i, j)-SolnBlk.Grid.Cell[i][j].Xc;
           uQuad[0] = SolnBlk.U[i][j].u + 
                      SolnBlk.dudx[i][j]*dX.x +
                      SolnBlk.dudy[i][j]*dX.y ;
           dX = SolnBlk.Grid.xfaceW(i, j)-SolnBlk.Grid.Cell[i][j].Xc;
           uQuad[1] = SolnBlk.U[i][j].u + 
                      SolnBlk.dudx[i][j]*dX.x +
                      SolnBlk.dudy[i][j]*dX.y ;
           dX = SolnBlk.Grid.xfaceN(i, j)-SolnBlk.Grid.Cell[i][j].Xc;
           uQuad[2] = SolnBlk.U[i][j].u  + 
                      SolnBlk.dudx[i][j]*dX.x +
                      SolnBlk.dudy[i][j]*dX.y ;
           dX = SolnBlk.Grid.xfaceS(i, j)-SolnBlk.Grid.Cell[i][j].Xc;
           uQuad[3] = SolnBlk.U[i][j].u + 
                      SolnBlk.dudx[i][j]*dX.x +
                      SolnBlk.dudy[i][j]*dX.y ;
    
           switch(Limiter) {
             case LIMITER_ONE :
               phi = ONE;
               break;
             case LIMITER_ZERO :
               phi = ZERO;
               break;
             case LIMITER_BARTH_JESPERSEN :
               phi = Limiter_BarthJespersen(uQuad, SolnBlk.U[i][j].u, 
                                            u0Min, u0Max, 4);
               break;
             case LIMITER_VENKATAKRISHNAN :
               phi = Limiter_Venkatakrishnan(uQuad, SolnBlk.U[i][j].u, 
                                             u0Min, u0Max, 4);
               break;
             case LIMITER_VANLEER :
               phi = Limiter_VanLeer(uQuad, SolnBlk.U[i][j].u, 
                                     u0Min, u0Max, 4);
               break;
             case LIMITER_VANALBADA :
               phi = Limiter_VanAlbada(uQuad, SolnBlk.U[i][j].u, 
                                       u0Min, u0Max, 4);
               break;
             default:
               phi = Limiter_BarthJespersen(uQuad, SolnBlk.U[i][j].u, 
                                            u0Min, u0Max, 4);
               break;
           } /* endswitch */
    
           SolnBlk.phi[i][j] = phi;
        } /* endif */
    } else {
        SolnBlk.dudx[i][j] = ZERO;
        SolnBlk.dudy[i][j] = ZERO; 
        SolnBlk.phi[i][j]  = ZERO;
    } /* endif */
    
}

/********************************************************
 * Routine: Linear_Reconstruction_GreenGauss2           *
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
void Linear_Reconstruction_GreenGauss2(AdvectDiffuse2D_Quad_Block &SolnBlk,
				       const int i, 
                                       const int j,
                                       const int Limiter) {

    int n, n_pts, i_index[8], j_index[8];
    double u0Min, u0Max, uQuad[4], phi;
    double DxDx_ave, DxDy_ave, DyDy_ave;
    double l_north, l_south, l_east, l_west;
    Vector2D n_north, n_south, n_east, n_west, dX;
    double u_nw, u_ne, u_sw, u_se, u_face, Du, DuDx_ave, DuDy_ave;

    /* Carry out the limited solution reconstruction in
       each cell of the computational mesh. */

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
                 SolnBlk.Grid.BCtypeW[j] == BC_NEUMANN ||
                 SolnBlk.Grid.BCtypeW[j] == BC_ROBIN) {
         n_pts = 8;
         i_index[0] = i-1; j_index[0] = j-1;
         i_index[1] = i  ; j_index[1] = j-1;
         i_index[2] = i+1; j_index[2] = j-1;
         i_index[3] = i-1; j_index[3] = j  ;
         i_index[4] = i+1; j_index[4] = j  ;
         i_index[5] = i-1; j_index[5] = j+1;
         i_index[6] = i  ; j_index[6] = j+1;
         i_index[7] = i+1; j_index[7] = j+1;
//           if (j == SolnBlk.JCl) {
//              n_pts = 5;
//              i_index[0] = i-1; j_index[0] = j  ;
//              i_index[1] = i+1; j_index[1] = j  ;
//              i_index[2] = i-1; j_index[2] = j+1;
//              i_index[3] = i  ; j_index[3] = j+1;
//              i_index[4] = i+1; j_index[4] = j+1;
//           } else if (j == SolnBlk.JCu) {
//              n_pts = 5;
//              i_index[0] = i-1; j_index[0] = j-1;
//              i_index[1] = i  ; j_index[1] = j-1;
//              i_index[2] = i+1; j_index[2] = j-1;
//              i_index[3] = i-1; j_index[3] = j  ;
//              i_index[4] = i+1; j_index[4] = j  ;
//           } else {
//              n_pts = 8;
//              i_index[0] = i-1; j_index[0] = j-1;
//              i_index[1] = i  ; j_index[1] = j-1;
//              i_index[2] = i+1; j_index[2] = j-1;
//              i_index[3] = i-1; j_index[3] = j  ;
//              i_index[4] = i+1; j_index[4] = j  ;
//              i_index[5] = i-1; j_index[5] = j+1;
//              i_index[6] = i  ; j_index[6] = j+1;
//              i_index[7] = i+1; j_index[7] = j+1;
//           } /* endif */
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
//           if (j == SolnBlk.JCl) {
//              n_pts = 3;
//              i_index[0] = i+1; j_index[0] = j  ;
//              i_index[1] = i  ; j_index[1] = j+1;
//              i_index[2] = i+1; j_index[2] = j+1;
//           } else if (j == SolnBlk.JCu) {
//              n_pts = 3;
//              i_index[0] = i  ; j_index[0] = j-1;
//              i_index[1] = i+1; j_index[1] = j-1;
//              i_index[2] = i+1; j_index[2] = j  ;
//           } else {
//              n_pts = 5;
//              i_index[0] = i  ; j_index[0] = j-1;
//              i_index[1] = i+1; j_index[1] = j-1;
//              i_index[2] = i+1; j_index[2] = j  ;
//              i_index[3] = i  ; j_index[3] = j+1;
//              i_index[4] = i+1; j_index[4] = j+1;
//           } /* endif */
      } /* endif */           
    } else if ((i == SolnBlk.ICu+SolnBlk.Nghost-1) && 
               (SolnBlk.Grid.BCtypeE[j] != BC_NONE)) {
      if (j == SolnBlk.JCl-SolnBlk.Nghost+1 || j == SolnBlk.JCu+SolnBlk.Nghost-1) {
         n_pts = 0;
      } else if (SolnBlk.Grid.BCtypeE[j] == BC_PERIODIC ||
                 SolnBlk.Grid.BCtypeE[j] == BC_NEUMANN ||
                 SolnBlk.Grid.BCtypeE[j] == BC_ROBIN) {
         n_pts = 8;
         i_index[0] = i-1; j_index[0] = j-1;
         i_index[1] = i  ; j_index[1] = j-1;
         i_index[2] = i+1; j_index[2] = j-1;
         i_index[3] = i-1; j_index[3] = j  ;
         i_index[4] = i+1; j_index[4] = j  ;
         i_index[5] = i-1; j_index[5] = j+1;
         i_index[6] = i  ; j_index[6] = j+1;
         i_index[7] = i+1; j_index[7] = j+1;
//           if (j == SolnBlk.JCl) {
//              n_pts = 5;
//              i_index[0] = i-1; j_index[0] = j  ;
//              i_index[1] = i+1; j_index[1] = j  ;
//              i_index[2] = i-1; j_index[2] = j+1;
//              i_index[3] = i  ; j_index[3] = j+1;
//              i_index[4] = i+1; j_index[4] = j+1;
//           } else if (j == SolnBlk.JCu) {
//              n_pts = 5;
//              i_index[0] = i-1; j_index[0] = j-1;
//              i_index[1] = i  ; j_index[1] = j-1;
//              i_index[2] = i+1; j_index[2] = j-1;
//              i_index[3] = i-1; j_index[3] = j  ;
//              i_index[4] = i+1; j_index[4] = j  ;
//           } else {
//              n_pts = 8;
//              i_index[0] = i-1; j_index[0] = j-1;
//              i_index[1] = i  ; j_index[1] = j-1;
//              i_index[2] = i+1; j_index[2] = j-1;
//              i_index[3] = i-1; j_index[3] = j  ;
//              i_index[4] = i+1; j_index[4] = j  ;
//              i_index[5] = i-1; j_index[5] = j+1;
//              i_index[6] = i  ; j_index[6] = j+1;
//              i_index[7] = i+1; j_index[7] = j+1;
//           } /* endif */
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
//           if (j == SolnBlk.JCl) {
//              n_pts = 3;
//              i_index[0] = i-1; j_index[0] = j  ;
//              i_index[1] = i-1; j_index[1] = j+1;
//              i_index[2] = i  ; j_index[2] = j+1;
//           } else if (j == SolnBlk.JCu) {
//              n_pts = 3;
//              i_index[0] = i-1; j_index[0] = j-1;
//              i_index[1] = i  ; j_index[1] = j-1;
//              i_index[2] = i-1; j_index[2] = j  ;
//           } else {
//              n_pts = 5;
//              i_index[0] = i-1; j_index[0] = j-1;
//              i_index[1] = i  ; j_index[1] = j-1;
//              i_index[2] = i-1; j_index[2] = j  ;
//              i_index[3] = i-1; j_index[3] = j+1;
//              i_index[4] = i  ; j_index[4] = j+1;
//           } /* endif */
      } /* endif */
    } else if ((j == SolnBlk.JCl-SolnBlk.Nghost+1) && 
               (SolnBlk.Grid.BCtypeS[i] != BC_NONE)) {
      if (i == SolnBlk.ICl-SolnBlk.Nghost+1 || i == SolnBlk.ICu+SolnBlk.Nghost-1) {
         n_pts = 0;
      } else if (SolnBlk.Grid.BCtypeS[i] == BC_PERIODIC ||
                 SolnBlk.Grid.BCtypeS[i] == BC_NEUMANN ||
                 SolnBlk.Grid.BCtypeS[i] == BC_ROBIN) {
         n_pts = 8;
         i_index[0] = i-1; j_index[0] = j-1;
         i_index[1] = i  ; j_index[1] = j-1;
         i_index[2] = i+1; j_index[2] = j-1;
         i_index[3] = i-1; j_index[3] = j  ;
         i_index[4] = i+1; j_index[4] = j  ;
         i_index[5] = i-1; j_index[5] = j+1;
         i_index[6] = i  ; j_index[6] = j+1;
         i_index[7] = i+1; j_index[7] = j+1;
//           if (i == SolnBlk.ICl) {
//              n_pts = 5;
//              i_index[0] = i  ; j_index[0] = j-1;
//              i_index[1] = i+1; j_index[1] = j-1;
//              i_index[2] = i+1; j_index[2] = j  ;
//              i_index[3] = i  ; j_index[3] = j+1;
//              i_index[4] = i+1; j_index[4] = j+1;
//           } else if (i == SolnBlk.ICu) {
//              n_pts = 5;
//              i_index[0] = i-1; j_index[0] = j-1;
//              i_index[1] = i  ; j_index[1] = j-1;
//              i_index[2] = i-1; j_index[2] = j  ;
//              i_index[3] = i-1; j_index[3] = j+1;
//              i_index[4] = i  ; j_index[4] = j+1;
//           } else {
//              n_pts = 8;
//              i_index[0] = i-1; j_index[0] = j-1;
//              i_index[1] = i  ; j_index[1] = j-1;
//              i_index[2] = i+1; j_index[2] = j-1;
//              i_index[3] = i-1; j_index[3] = j  ;
//              i_index[4] = i+1; j_index[4] = j  ;
//              i_index[5] = i-1; j_index[5] = j+1;
//              i_index[6] = i  ; j_index[6] = j+1;
//              i_index[7] = i+1; j_index[7] = j+1;
//           } /* endif */
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
//           if (i == SolnBlk.ICl) {
//              n_pts = 3;
//              i_index[0] = i+1; j_index[0] = j  ;
//              i_index[1] = i  ; j_index[1] = j+1;
//              i_index[2] = i+1; j_index[2] = j+1;
//           } else if (i == SolnBlk.ICu) {
//              n_pts = 3;
//              i_index[0] = i-1; j_index[0] = j  ;
//              i_index[1] = i-1; j_index[1] = j+1;
//              i_index[2] = i  ; j_index[2] = j+1;
//           } else {
//              n_pts = 5;
//              i_index[0] = i-1; j_index[0] = j  ;
//              i_index[1] = i+1; j_index[1] = j  ;
//              i_index[2] = i-1; j_index[2] = j+1;
//              i_index[3] = i  ; j_index[3] = j+1;
//              i_index[4] = i+1; j_index[4] = j+1;
//           } /* endif */
      } /* endif */
    } else if ((j == SolnBlk.JCu+SolnBlk.Nghost-1) && 
               (SolnBlk.Grid.BCtypeN[i] != BC_NONE)) {
      if (i == SolnBlk.ICl-SolnBlk.Nghost+1 || i == SolnBlk.ICu+SolnBlk.Nghost-1) {
         n_pts = 0;
      } else if (SolnBlk.Grid.BCtypeN[i] == BC_PERIODIC ||
                 SolnBlk.Grid.BCtypeN[i] == BC_NEUMANN ||
                 SolnBlk.Grid.BCtypeN[i] == BC_ROBIN) {
         n_pts = 8;
         i_index[0] = i-1; j_index[0] = j-1;
         i_index[1] = i  ; j_index[1] = j-1;
         i_index[2] = i+1; j_index[2] = j-1;
         i_index[3] = i-1; j_index[3] = j  ;
         i_index[4] = i+1; j_index[4] = j  ;
         i_index[5] = i-1; j_index[5] = j+1;
         i_index[6] = i  ; j_index[6] = j+1;
         i_index[7] = i+1; j_index[7] = j+1;
//           if (i == SolnBlk.ICl) {
//              n_pts = 5;
//              i_index[0] = i  ; j_index[0] = j-1;
//              i_index[1] = i+1; j_index[1] = j-1;
//              i_index[2] = i+1; j_index[2] = j  ;
//              i_index[3] = i  ; j_index[3] = j+1;
//              i_index[4] = i+1; j_index[4] = j+1;
//           } else if (i == SolnBlk.ICu) {
//              n_pts = 5;
//              i_index[0] = i-1; j_index[0] = j-1;
//              i_index[1] = i  ; j_index[1] = j-1;
//              i_index[2] = i-1; j_index[2] = j  ;
//              i_index[3] = i-1; j_index[3] = j+1;
//              i_index[4] = i  ; j_index[4] = j+1;
//           } else {
//              n_pts = 8;
//              i_index[0] = i-1; j_index[0] = j-1;
//              i_index[1] = i  ; j_index[1] = j-1;
//              i_index[2] = i+1; j_index[2] = j-1;
//              i_index[3] = i-1; j_index[3] = j  ;
//              i_index[4] = i+1; j_index[4] = j  ;
//              i_index[5] = i-1; j_index[5] = j+1;
//              i_index[6] = i  ; j_index[6] = j+1;
//              i_index[7] = i+1; j_index[7] = j+1;
//           } /* endif */
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
//           if (i == SolnBlk.ICl) {
//              n_pts = 3;
//              i_index[0] = i  ; j_index[0] = j-1;
//              i_index[1] = i+1; j_index[1] = j-1;
//              i_index[2] = i+1; j_index[2] = j  ;
//           } else if (i == SolnBlk.ICu) {
//              n_pts = 3;
//              i_index[0] = i-1; j_index[0] = j-1;
//              i_index[1] = i  ; j_index[1] = j-1;
//              i_index[2] = i-1; j_index[2] = j  ;
//           } else {
//              n_pts = 5;
//              i_index[0] = i-1; j_index[0] = j-1;
//              i_index[1] = i  ; j_index[1] = j-1;
//              i_index[2] = i+1; j_index[2] = j-1;
//              i_index[3] = i-1; j_index[3] = j  ;
//              i_index[4] = i+1; j_index[4] = j  ;
//           } /* endif */
      } /* endif */
//     } else if ((i == SolnBlk.ICl && j == SolnBlk.JCl) && 
//                (SolnBlk.Grid.BCtypeW[j] != BC_NONE &&
//                 SolnBlk.Grid.BCtypeS[i] != BC_NONE)) {
//       n_pts = 7;
//       i_index[0] = i  ; j_index[0] = j-1;
//       i_index[1] = i+1; j_index[1] = j-1;
//       i_index[2] = i-1; j_index[2] = j  ;
//       i_index[3] = i+1; j_index[3] = j  ;
//       i_index[4] = i-1; j_index[4] = j+1;
//       i_index[5] = i  ; j_index[5] = j+1;
//       i_index[6] = i+1; j_index[6] = j+1;
//     } else if ((i == SolnBlk.ICl && j == SolnBlk.JCu) && 
//                (SolnBlk.Grid.BCtypeW[j] != BC_NONE &&
//                 SolnBlk.Grid.BCtypeN[i] != BC_NONE)) {
//       n_pts = 7;
//       i_index[0] = i-1; j_index[0] = j-1;
//       i_index[1] = i  ; j_index[1] = j-1;
//       i_index[2] = i+1; j_index[2] = j-1;
//       i_index[3] = i-1; j_index[3] = j  ;
//       i_index[4] = i+1; j_index[4] = j  ;
//       i_index[5] = i  ; j_index[5] = j+1;
//       i_index[6] = i+1; j_index[6] = j+1;
//     } else if ((i == SolnBlk.ICu && j == SolnBlk.JCl) && 
//                (SolnBlk.Grid.BCtypeE[j] != BC_NONE &&
//                 SolnBlk.Grid.BCtypeS[i] != BC_NONE)) {
//       n_pts = 7;
//       i_index[0] = i-1; j_index[0] = j-1;
//       i_index[1] = i  ; j_index[1] = j-1;
//       i_index[2] = i-1; j_index[2] = j  ;
//       i_index[3] = i+1; j_index[3] = j  ;
//       i_index[4] = i-1; j_index[4] = j+1;
//       i_index[5] = i  ; j_index[5] = j+1;
//       i_index[6] = i+1; j_index[6] = j+1;
//     } else if ((i == SolnBlk.ICu && j == SolnBlk.JCu) && 
//                (SolnBlk.Grid.BCtypeE[j] != BC_NONE &&
//                 SolnBlk.Grid.BCtypeN[i] != BC_NONE)) {
//       n_pts = 7;
//       i_index[0] = i-1; j_index[0] = j-1;
//       i_index[1] = i  ; j_index[1] = j-1;
//       i_index[2] = i+1; j_index[2] = j-1;
//       i_index[3] = i-1; j_index[3] = j  ;
//       i_index[4] = i+1; j_index[4] = j  ;
//       i_index[5] = i-1; j_index[5] = j+1;
//       i_index[6] = i  ; j_index[6] = j+1;
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

	  // USE NODE AVERAGING. 
//             u_nw = SolnBlk.unNW(i, j);
//             u_ne = SolnBlk.unNE(i, j);
//             u_sw = SolnBlk.unSW(i, j);
//             u_se = SolnBlk.unSE(i, j);

//             l_north = SolnBlk.Grid.lfaceN(i, j);
//             l_south = SolnBlk.Grid.lfaceS(i, j);
//             l_east = SolnBlk.Grid.lfaceE(i, j);
//             l_west = SolnBlk.Grid.lfaceW(i, j);

//             n_north = SolnBlk.Grid.nfaceN(i, j);
//             n_south = SolnBlk.Grid.nfaceS(i, j);
//             n_east = SolnBlk.Grid.nfaceE(i, j);
//             n_west = SolnBlk.Grid.nfaceW(i, j);

//             u_face = HALF*(u_nw+u_ne)*l_north; 
//             SolnBlk.dudx[i][j] = u_face*n_north.x;
//  	   SolnBlk.dudy[i][j] = u_face*n_north.y;

//             u_face = HALF*(u_sw+u_se)*l_south; 
//             SolnBlk.dudx[i][j] += u_face*n_south.x;
//  	   SolnBlk.dudy[i][j] += u_face*n_south.y;

//             u_face = HALF*(u_ne+u_se)*l_east; 
//             SolnBlk.dudx[i][j] += u_face*n_east.x;
//  	   SolnBlk.dudy[i][j] += u_face*n_east.y;

//             u_face = HALF*(u_nw+u_sw)*l_west; 
//             SolnBlk.dudx[i][j] += u_face*n_west.x;
//  	   SolnBlk.dudy[i][j] += u_face*n_west.y;

	  // USE CELL AVERAGING. (SAME AS THE ONE DONE ON PAPERS)
           l_north = SolnBlk.Grid.lfaceN(i, j);
           l_south = SolnBlk.Grid.lfaceS(i, j);
           l_east = SolnBlk.Grid.lfaceE(i, j);
           l_west = SolnBlk.Grid.lfaceW(i, j);

           n_north = SolnBlk.Grid.nfaceN(i, j);
           n_south = SolnBlk.Grid.nfaceS(i, j);
           n_east = SolnBlk.Grid.nfaceE(i, j);
           n_west = SolnBlk.Grid.nfaceW(i, j);

           u_face = HALF*(SolnBlk.U[i][j].u+SolnBlk.U[i][j+1].u)*l_north; 
           SolnBlk.dudx[i][j] = u_face*n_north.x;
	   SolnBlk.dudy[i][j] = u_face*n_north.y;

           u_face = HALF*(SolnBlk.U[i][j].u+SolnBlk.U[i][j-1].u)*l_south; 
           SolnBlk.dudx[i][j] += u_face*n_south.x;
	   SolnBlk.dudy[i][j] += u_face*n_south.y;

           u_face = HALF*(SolnBlk.U[i][j].u+SolnBlk.U[i+1][j].u)*l_east; 
           SolnBlk.dudx[i][j] += u_face*n_east.x;
	   SolnBlk.dudy[i][j] += u_face*n_east.y;

           u_face = HALF*(SolnBlk.U[i][j].u+SolnBlk.U[i-1][j].u)*l_west; 
           SolnBlk.dudx[i][j] += u_face*n_west.x;
	   SolnBlk.dudy[i][j] += u_face*n_west.y;

	   SolnBlk.dudx[i][j] = SolnBlk.dudx[i][j]/
	     SolnBlk.Grid.Cell[i][j].A;
	   SolnBlk.dudy[i][j] = SolnBlk.dudy[i][j]/
	     SolnBlk.Grid.Cell[i][j].A;

        // If <8 neighbours are used, apply least-squares reconstruction
        } else {
	  
           DuDx_ave = ZERO;
           DuDy_ave = ZERO;
           DxDx_ave = ZERO;
           DxDy_ave = ZERO;
           DyDy_ave = ZERO;
    
           for ( n = 0 ; n <= n_pts-1 ; ++n ) {
               dX = SolnBlk.Grid.Cell[ i_index[n] ][ j_index[n] ].Xc - 
                    SolnBlk.Grid.Cell[i][j].Xc;
               Du = SolnBlk.U[ i_index[n] ][ j_index[n] ].u - 
                    SolnBlk.U[i][j].u;
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
           SolnBlk.dudx[i][j] = (DuDx_ave*DyDy_ave-DuDy_ave*DxDy_ave)/
                                (DxDx_ave*DyDy_ave-DxDy_ave*DxDy_ave);
           SolnBlk.dudy[i][j] = (DuDy_ave*DxDx_ave-DuDx_ave*DxDy_ave)/
                                (DxDx_ave*DyDy_ave-DxDy_ave*DxDy_ave);
        } /* endif */
    
        // Calculate slope limiter.
	if (!SolnBlk.Freeze_Limiter) {
           u0Min = SolnBlk.U[i][j].u;
           u0Max = u0Min;
           for ( n = 0 ; n <= n_pts-1 ; ++n ) {
              u0Min = min(u0Min, SolnBlk.U[ i_index[n] ][ j_index[n] ].u);
              u0Max = max(u0Max, SolnBlk.U[ i_index[n] ][ j_index[n] ].u);
           } /* endfor */
    
           dX = SolnBlk.Grid.xfaceE(i, j)-SolnBlk.Grid.Cell[i][j].Xc;
           uQuad[0] = SolnBlk.U[i][j].u + 
                      SolnBlk.dudx[i][j]*dX.x +
                      SolnBlk.dudy[i][j]*dX.y ;
           dX = SolnBlk.Grid.xfaceW(i, j)-SolnBlk.Grid.Cell[i][j].Xc;
           uQuad[1] = SolnBlk.U[i][j].u + 
                      SolnBlk.dudx[i][j]*dX.x +
                      SolnBlk.dudy[i][j]*dX.y ;
           dX = SolnBlk.Grid.xfaceN(i, j)-SolnBlk.Grid.Cell[i][j].Xc;
           uQuad[2] = SolnBlk.U[i][j].u  + 
                      SolnBlk.dudx[i][j]*dX.x +
                      SolnBlk.dudy[i][j]*dX.y ;
           dX = SolnBlk.Grid.xfaceS(i, j)-SolnBlk.Grid.Cell[i][j].Xc;
           uQuad[3] = SolnBlk.U[i][j].u + 
                      SolnBlk.dudx[i][j]*dX.x +
                      SolnBlk.dudy[i][j]*dX.y ;
    
           switch(Limiter) {
             case LIMITER_ONE :
               phi = ONE;
               break;
             case LIMITER_ZERO :
               phi = ZERO;
               break;
             case LIMITER_BARTH_JESPERSEN :
               phi = Limiter_BarthJespersen(uQuad, SolnBlk.U[i][j].u, 
                                            u0Min, u0Max, 4);
               break;
             case LIMITER_VENKATAKRISHNAN :
               phi = Limiter_Venkatakrishnan(uQuad, SolnBlk.U[i][j].u, 
                                             u0Min, u0Max, 4);
               break;
             case LIMITER_VANLEER :
               phi = Limiter_VanLeer(uQuad, SolnBlk.U[i][j].u, 
                                     u0Min, u0Max, 4);
               break;
             case LIMITER_VANALBADA :
               phi = Limiter_VanAlbada(uQuad, SolnBlk.U[i][j].u, 
                                       u0Min, u0Max, 4);
               break;
             default:
               phi = Limiter_BarthJespersen(uQuad, SolnBlk.U[i][j].u, 
                                            u0Min, u0Max, 4);
               break;
           } /* endswitch */
    
           SolnBlk.phi[i][j] = phi;
        } /* endif */
    } else {
        SolnBlk.dudx[i][j] = ZERO;
        SolnBlk.dudy[i][j] = ZERO; 
        SolnBlk.phi[i][j]  = ZERO;
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
void Linear_Reconstruction_GreenGauss(AdvectDiffuse2D_Quad_Block &SolnBlk,
				      const int Limiter) {

    int i, j;

    /* Carry out the limited solution reconstruction in
       each cell of the computational mesh. */

    for ( j  = SolnBlk.JCl-SolnBlk.Nghost+1 ; j <= SolnBlk.JCu+SolnBlk.Nghost-1 ; ++j ) {
       for ( i = SolnBlk.ICl-SolnBlk.Nghost+1 ; i <= SolnBlk.ICu+SolnBlk.Nghost-1 ; ++i ) {
	   //Linear_Reconstruction_GreenGauss(SolnBlk, i, j, Limiter);
	   Linear_Reconstruction_GreenGauss2(SolnBlk, i, j, Limiter);
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
    if (i == SolnBlk.ICl-SolnBlk.Nghost || i == SolnBlk.ICu+SolnBlk.Nghost ||
        j == SolnBlk.JCl-SolnBlk.Nghost || j == SolnBlk.JCu+SolnBlk.Nghost) {
      n_pts = 0;
    } else if ((i == SolnBlk.ICl-SolnBlk.Nghost+1) && 
               (SolnBlk.Grid.BCtypeW[j] != BC_NONE)) {
      if (j == SolnBlk.JCl-SolnBlk.Nghost+1 || j == SolnBlk.JCu+SolnBlk.Nghost-1) {
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
    } else if ((i == SolnBlk.ICu+SolnBlk.Nghost-1) && 
               (SolnBlk.Grid.BCtypeE[j] != BC_NONE)) {
      if (j == SolnBlk.JCl-SolnBlk.Nghost+1 || j == SolnBlk.JCu+SolnBlk.Nghost-1) {
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
    } else if ((j == SolnBlk.JCl-SolnBlk.Nghost+1) && 
               (SolnBlk.Grid.BCtypeS[i] != BC_NONE)) {
      if (i == SolnBlk.ICl-SolnBlk.Nghost+1 || i == SolnBlk.ICu+SolnBlk.Nghost-1) {
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
    } else if ((j == SolnBlk.JCu+SolnBlk.Nghost-1) && 
               (SolnBlk.Grid.BCtypeN[i] != BC_NONE)) {
      if (i == SolnBlk.ICl-SolnBlk.Nghost+1 || i == SolnBlk.ICu+SolnBlk.Nghost-1) {
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
//     } else if ((i == SolnBlk.ICl && j == SolnBlk.JCl) && 
//                (SolnBlk.Grid.BCtypeW[j] != BC_NONE &&
//                 SolnBlk.Grid.BCtypeS[i] != BC_NONE)) {
//       n_pts = 7;
//       i_index[0] = i  ; j_index[0] = j-1;
//       i_index[1] = i+1; j_index[1] = j-1;
//       i_index[2] = i-1; j_index[2] = j  ;
//       i_index[3] = i+1; j_index[3] = j  ;
//       i_index[4] = i-1; j_index[4] = j+1;
//       i_index[5] = i  ; j_index[5] = j+1;
//       i_index[6] = i+1; j_index[6] = j+1;
//     } else if ((i == SolnBlk.ICl && j == SolnBlk.JCu) && 
//                (SolnBlk.Grid.BCtypeW[j] != BC_NONE &&
//                 SolnBlk.Grid.BCtypeN[i] != BC_NONE)) {
//       n_pts = 7;
//       i_index[0] = i-1; j_index[0] = j-1;
//       i_index[1] = i  ; j_index[1] = j-1;
//       i_index[2] = i+1; j_index[2] = j-1;
//       i_index[3] = i-1; j_index[3] = j  ;
//       i_index[4] = i+1; j_index[4] = j  ;
//       i_index[5] = i  ; j_index[5] = j+1;
//       i_index[6] = i+1; j_index[6] = j+1;
//     } else if ((i == SolnBlk.ICu && j == SolnBlk.JCl) && 
//                (SolnBlk.Grid.BCtypeE[j] != BC_NONE &&
//                 SolnBlk.Grid.BCtypeS[i] != BC_NONE)) {
//       n_pts = 7;
//       i_index[0] = i-1; j_index[0] = j-1;
//       i_index[1] = i  ; j_index[1] = j-1;
//       i_index[2] = i-1; j_index[2] = j  ;
//       i_index[3] = i+1; j_index[3] = j  ;
//       i_index[4] = i-1; j_index[4] = j+1;
//       i_index[5] = i  ; j_index[5] = j+1;
//       i_index[6] = i+1; j_index[6] = j+1;
//     } else if ((i == SolnBlk.ICu && j == SolnBlk.JCu) && 
//                (SolnBlk.Grid.BCtypeE[j] != BC_NONE &&
//                 SolnBlk.Grid.BCtypeN[i] != BC_NONE)) {
//       n_pts = 7;
//       i_index[0] = i-1; j_index[0] = j-1;
//       i_index[1] = i  ; j_index[1] = j-1;
//       i_index[2] = i+1; j_index[2] = j-1;
//       i_index[3] = i-1; j_index[3] = j  ;
//       i_index[4] = i+1; j_index[4] = j  ;
//       i_index[5] = i-1; j_index[5] = j+1;
//       i_index[6] = i  ; j_index[6] = j+1;
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
        DuDx_ave = ZERO;
        DuDy_ave = ZERO;
        DxDx_ave = ZERO;
        DxDy_ave = ZERO;
        DyDy_ave = ZERO;
    
        for ( n = 0 ; n <= n_pts-1 ; ++n ) {
            dX = SolnBlk.Grid.Cell[ i_index[n] ][ j_index[n] ].Xc - 
                 SolnBlk.Grid.Cell[i][j].Xc;
            Du = SolnBlk.U[ i_index[n] ][ j_index[n] ].u - 
                 SolnBlk.U[i][j].u;
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
        SolnBlk.dudx[i][j] = (DuDx_ave*DyDy_ave-DuDy_ave*DxDy_ave)/
                             (DxDx_ave*DyDy_ave-DxDy_ave*DxDy_ave);
        SolnBlk.dudy[i][j] = (DuDy_ave*DxDx_ave-DuDx_ave*DxDy_ave)/
                             (DxDx_ave*DyDy_ave-DxDy_ave*DxDy_ave);
    
        // Calculate slope limiter.
	if (!SolnBlk.Freeze_Limiter) {
           u0Min = SolnBlk.U[i][j].u;
           u0Max = u0Min;
           for ( n = 0 ; n <= n_pts-1 ; ++n ) {
              u0Min = min(u0Min, SolnBlk.U[ i_index[n] ][ j_index[n] ].u);
              u0Max = max(u0Max, SolnBlk.U[ i_index[n] ][ j_index[n] ].u);
           } /* endfor */
    
           dX = SolnBlk.Grid.xfaceE(i, j)-SolnBlk.Grid.Cell[i][j].Xc;
           uQuad[0] = SolnBlk.U[i][j].u + 
                      SolnBlk.dudx[i][j]*dX.x +
                      SolnBlk.dudy[i][j]*dX.y ;
           dX = SolnBlk.Grid.xfaceW(i, j)-SolnBlk.Grid.Cell[i][j].Xc;
           uQuad[1] = SolnBlk.U[i][j].u + 
                      SolnBlk.dudx[i][j]*dX.x +
                      SolnBlk.dudy[i][j]*dX.y ;
           dX = SolnBlk.Grid.xfaceN(i, j)-SolnBlk.Grid.Cell[i][j].Xc;
           uQuad[2] = SolnBlk.U[i][j].u + 
                      SolnBlk.dudx[i][j]*dX.x +
                      SolnBlk.dudy[i][j]*dX.y ;
           dX = SolnBlk.Grid.xfaceS(i, j)-SolnBlk.Grid.Cell[i][j].Xc;
           uQuad[3] = SolnBlk.U[i][j].u + 
                      SolnBlk.dudx[i][j]*dX.x +
                      SolnBlk.dudy[i][j]*dX.y ;
    
           switch(Limiter) {
             case LIMITER_ONE :
               phi = ONE;
               break;
             case LIMITER_ZERO :
               phi = ZERO;
               break;
             case LIMITER_BARTH_JESPERSEN :
               phi = Limiter_BarthJespersen(uQuad, SolnBlk.U[i][j].u, 
                                            u0Min, u0Max, 4);
               break;
             case LIMITER_VENKATAKRISHNAN :
               phi = Limiter_Venkatakrishnan(uQuad, SolnBlk.U[i][j].u, 
                                             u0Min, u0Max, 4);
               break;
             case LIMITER_VANLEER :
               phi = Limiter_VanLeer(uQuad, SolnBlk.U[i][j].u, 
                                     u0Min, u0Max, 4);
               break;
             case LIMITER_VANALBADA :
               phi = Limiter_VanAlbada(uQuad, SolnBlk.U[i][j].u, 
                                       u0Min, u0Max, 4);
               break;
             default:
               phi = Limiter_BarthJespersen(uQuad, SolnBlk.U[i][j].u, 
                                            u0Min, u0Max, 4);
               break;
           } /* endswitch */
    
           SolnBlk.phi[i][j] = phi;
        } /* endif */
    } else {
        SolnBlk.dudx[i][j] = ZERO;
        SolnBlk.dudy[i][j] = ZERO; 
        SolnBlk.phi[i][j]  = ZERO;
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

/********************************************************
 * Routine: Diffusive_Flux                              *
 *                                                      *
 * Evaluates the diffusive flux for each cell of the    *
 * computational mesh for the specified quadrilateral   *
 * solution block.                                      *
 *                                                      *
 ********************************************************/
void Diffusive_Flux(AdvectDiffuse2D_Quad_Block &SolnBlk) {

    int i, j;

    for ( j  = SolnBlk.JCl-SolnBlk.Nghost+1; j <= SolnBlk.JCu+SolnBlk.Nghost-1; ++j ) {
       for ( i = SolnBlk.ICl-SolnBlk.Nghost+1; i <= SolnBlk.ICu+SolnBlk.Nghost-1; ++i ) {
	   SolnBlk.evalDiffusiveFlux(i, j);
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
void Residual_Smoothing(AdvectDiffuse2D_Quad_Block &SolnBlk,
                        const int k_residual,
			double &epsilon,
                        const int number_of_Gauss_Seidel_iterations) {

    int i, j, n;

    for ( j  = SolnBlk.JCl+1; j <= SolnBlk.JCu-1; ++j ) {
       for ( i = SolnBlk.ICl+1; i <= SolnBlk.ICu-1; ++i ) {
          SolnBlk.dudt[i][j][k_residual+1] =  SolnBlk.dudt[i][j][k_residual];
       } /* endfor */
    } /* endfor */

    for ( n  = 1; n <= number_of_Gauss_Seidel_iterations; ++n ) {

       for ( j  = SolnBlk.JCl+1; j <= SolnBlk.JCu-1; ++j ) {
          for ( i = SolnBlk.ICl+1; i <= SolnBlk.ICu-1; ++i ) {
             SolnBlk.dudt[i][j][k_residual+1] = (SolnBlk.dudt[i][j][k_residual]
                + epsilon*(SolnBlk.dudt[i  ][j-1][k_residual+1] +
                           SolnBlk.dudt[i-1][j  ][k_residual+1] +
                           SolnBlk.dudt[i+1][j  ][k_residual+1] +
                           SolnBlk.dudt[i  ][j+1][k_residual+1]))/(ONE + FOUR*epsilon);

//              SolnBlk.dudt[i][j][k_residual+1] = (SolnBlk.dudt[i][j][k_residual]
//                 + epsilon*(SolnBlk.dudt[i-1][j-1][k_residual+1] +
//                            SolnBlk.dudt[i  ][j-1][k_residual+1] +
//                            SolnBlk.dudt[i+1][j-1][k_residual+1] +
//                            SolnBlk.dudt[i-1][j  ][k_residual+1] +
//                            SolnBlk.dudt[i+1][j  ][k_residual+1] +
//                            SolnBlk.dudt[i-1][j+1][k_residual+1] +
//                            SolnBlk.dudt[i  ][j+1][k_residual+1] +
//                            SolnBlk.dudt[i+1][j+1][k_residual+1]))/(ONE + EIGHT*epsilon);
          } /* endfor */
       } /* endfor */

    } /* endfor */

    for ( j  = SolnBlk.JCl+1; j <= SolnBlk.JCu-1; ++j ) {
       for ( i = SolnBlk.ICl+1; i <= SolnBlk.ICu-1; ++i ) {
          SolnBlk.dudt[i][j][k_residual] =  SolnBlk.dudt[i][j][k_residual+1];
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
				   AdvectDiffuse2D_Input_Parameters &IP,
                                   int &number_refinement_criteria,
                                   AdvectDiffuse2D_Quad_Block &SolnBlk) {

    int i, j;

    double grad_u_x, grad_u_y, grad_u_abs, grad_u_criteria, grad_u_criteria_max;

    /* Set the number of refinement criteria to be used (1):
       (1) refinement criteria #1 based on the gradient of the solution. */

    number_refinement_criteria = 1;

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
             grad_u_x = SolnBlk.dudx[i][j];
             grad_u_y = SolnBlk.dudy[i][j];
             grad_u_abs = sqrt(sqr(grad_u_x) + sqr(grad_u_y));
             grad_u_criteria = sqrt(SolnBlk.Grid.Cell[i][j].A)*grad_u_abs/SolnBlk.U[i][j].u;
             grad_u_criteria_max = max(grad_u_criteria_max, grad_u_criteria);
          } /* endif */

       } /* endfor */
    } /* endfor */

    /* Return the refinement criteria. */

    refinement_criteria[0] = grad_u_criteria_max;

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
void Fix_Refined_Block_Boundaries(AdvectDiffuse2D_Quad_Block SolnBlk,
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
void Unfix_Refined_Block_Boundaries(AdvectDiffuse2D_Quad_Block SolnBlk) {

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
void Apply_Boundary_Flux_Corrections(AdvectDiffuse2D_Quad_Block SolnBlk,
                                     const int Number_Neighbours_North_Boundary,
                                     const int Number_Neighbours_South_Boundary,
                                     const int Number_Neighbours_East_Boundary,
                                     const int Number_Neighbours_West_Boundary) {

    int i, j;
 
    /* Correct the fluxes at the north boundary as required. */

    if (Number_Neighbours_North_Boundary == 2) {
       for ( i = SolnBlk.ICl ; i <= SolnBlk.ICu ; ++i ) {
          SolnBlk.dudt[i][SolnBlk.JCu][0] -= 
             SolnBlk.FluxN[i]/SolnBlk.Grid.Cell[i][SolnBlk.JCu].A;
       } /* endfor */
    } /* endif */

    /* Correct the fluxes at the south boundary as required. */

    if (Number_Neighbours_South_Boundary == 2) {
       for ( i = SolnBlk.ICl ; i <= SolnBlk.ICu ; ++i ) {
          SolnBlk.dudt[i][SolnBlk.JCl][0] -= 
             SolnBlk.FluxS[i]/SolnBlk.Grid.Cell[i][SolnBlk.JCl].A;
       } /* endfor */
    } /* endif */

    /* Correct the fluxes at the east boundary as required. */

    if (Number_Neighbours_East_Boundary == 2) {
       for ( j= SolnBlk.JCl ; j <= SolnBlk.JCu ; ++j ) {
          SolnBlk.dudt[SolnBlk.ICu][j][0] -= 
             SolnBlk.FluxE[j]/SolnBlk.Grid.Cell[SolnBlk.ICu][j].A;
       } /* endfor */
    } /* endif */

    /* Correct the fluxes at the west boundary as required. */

    if (Number_Neighbours_West_Boundary == 2) {
       for ( j= SolnBlk.JCl ; j <= SolnBlk.JCu ; ++j ) {
          SolnBlk.dudt[SolnBlk.ICl][j][0] -= 
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
void Apply_Boundary_Flux_Corrections_Multistage_Explicit(AdvectDiffuse2D_Quad_Block SolnBlk,
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
          SolnBlk.dudt[i][SolnBlk.JCu][k_residual] -= 
             (CFL_Number*SolnBlk.dt[i][SolnBlk.JCu])*SolnBlk.FluxN[i]/
             SolnBlk.Grid.Cell[i][SolnBlk.JCu].A;
       } /* endfor */
    } /* endif */

    /* Correct the fluxes at the south boundary as required. */

    if (Number_Neighbours_South_Boundary == 2) {
       for ( i = SolnBlk.ICl ; i <= SolnBlk.ICu ; ++i ) {
          SolnBlk.dudt[i][SolnBlk.JCl][k_residual] -= 
             (CFL_Number*SolnBlk.dt[i][SolnBlk.JCl])*SolnBlk.FluxS[i]/
             SolnBlk.Grid.Cell[i][SolnBlk.JCl].A;
       } /* endfor */
    } /* endif */

    /* Correct the fluxes at the east boundary as required. */

    if (Number_Neighbours_East_Boundary == 2) {
       for ( j= SolnBlk.JCl ; j <= SolnBlk.JCu ; ++j ) {
          SolnBlk.dudt[SolnBlk.ICu][j][k_residual] -= 
             (CFL_Number*SolnBlk.dt[SolnBlk.ICu][j])*SolnBlk.FluxE[j]/
             SolnBlk.Grid.Cell[SolnBlk.ICu][j].A;
       } /* endfor */
    } /* endif */

    /* Correct the fluxes at the west boundary as required. */

    if (Number_Neighbours_West_Boundary == 2) {
       for ( j= SolnBlk.JCl ; j <= SolnBlk.JCu ; ++j ) {
          SolnBlk.dudt[SolnBlk.ICl][j][k_residual] -= 
             (CFL_Number*SolnBlk.dt[SolnBlk.ICl][j])*SolnBlk.FluxW[j]/
             SolnBlk.Grid.Cell[SolnBlk.ICl][j].A;
       } /* endfor */
    } /* endif */

}

/********************************************************
 * Routine: dUdt_Residual_Evaluation                    *
 *                                                      *
 * This routine evaluate the residual for the specified *
 * solution block using a 2nd-ororder limited upwind    *
 * finite-volume spatial discretization scheme for the  *
 * convective flux coupled with a centrally-weighted    *
 * finite-volume discretization for the diffused flux.  *
 * The residual is stored in dUdt[][][0].               *
 *                                                      *
 ********************************************************/
int dUdt_Residual_Evaluation(AdvectDiffuse2D_Quad_Block &SolnBlk,
			     AdvectDiffuse2D_Input_Parameters &Input_Parameters) {

    int i, j;
    Vector2D dX;
    double ul, ur, flux;
    AdvectDiffuse2D_State Ul, Ur;

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

    /* Evaluate the diffusive fluxes in each cell
       of the computational grid for this stage. */

    Diffusive_Flux(SolnBlk);

    /* Evaluate the time rate of change of the solution
       (i.e., the solution residuals) using a second-order
       limited upwind finite-volume scheme for the convective 
       fluxes and a second-order centrally-weighted finite 
       volume discretization for the diffusive fluxes. */
    
    // Add i-direction (zeta-direction) fluxes.
    for ( j  = SolnBlk.JCl-1 ; j <= SolnBlk.JCu+1 ; ++j ) {
      SolnBlk.dudt[SolnBlk.ICl-1][j][0] = ZERO;
          
       for ( i = SolnBlk.ICl-1 ; i <= SolnBlk.ICu ; ++i ) {
	 SolnBlk.dudt[i+1][j][0] = ZERO;
	 
	 if ( j > SolnBlk.JCl-1 && j < SolnBlk.JCu+1 ) {
    
	   /* Evaluate the cell interface i-direction fluxes. */
	   
	   dX = SolnBlk.Grid.xfaceE(i, j)-SolnBlk.Grid.Cell[i][j].Xc;
	   ul = SolnBlk.U[i][j].u + 
	     (SolnBlk.phi[i][j]*SolnBlk.dudx[i][j])*dX.x +
	     (SolnBlk.phi[i][j]*SolnBlk.dudy[i][j])*dX.y;
	   dX = SolnBlk.Grid.xfaceW(i+1, j)-SolnBlk.Grid.Cell[i+1][j].Xc;
	   ur = SolnBlk.U[i+1][j].u + 
	     (SolnBlk.phi[i+1][j]*SolnBlk.dudx[i+1][j])*dX.x +
	     (SolnBlk.phi[i+1][j]*SolnBlk.dudy[i+1][j])*dX.y;
	   
	   Ul = SolnBlk.U[i][j]; 
	   Ul.u = ul;
	   Ur = SolnBlk.U[i+1][j]; 
	   Ur.u = ur;
	   
	   if ((i == SolnBlk.ICl-1) && 
	       (SolnBlk.Grid.BCtypeW[j] == BC_DIRICHLET)) {
	     Ul.Fd = Ur.Fd;
	   } else if ((i == SolnBlk.ICu) && 
		      (SolnBlk.Grid.BCtypeE[j] == BC_DIRICHLET)) {
	     Ur.Fd = Ul.Fd;
	   } /* endif */
	   
	   flux = Flux_n(Ul, Ur, SolnBlk.Grid.nfaceE(i, j));
	   
	   /* Evaluate cell-averaged solution changes. */
	   
	   SolnBlk.dudt[i][j][0] -= 
	     flux*SolnBlk.Grid.lfaceE(i, j)/
	     SolnBlk.Grid.Cell[i][j].A;
	   SolnBlk.dudt[i+1][j][0] +=
	     flux*SolnBlk.Grid.lfaceW(i+1, j)/
	     SolnBlk.Grid.Cell[i+1][j].A;
	   
           /* Include regular source terms. */

           SolnBlk.dudt[i][j][0] += s(SolnBlk.U[i][j]);

           /* Include axisymmetric source terms as required. */

	   if (SolnBlk.Axisymmetric) {
               SolnBlk.dudt[i][j][0] += 
	        s_axi(SolnBlk.U[i][j], SolnBlk.Grid.Cell[i][j].Xc);
           } /* endif */

           /* Save west and east face boundary flux. */

           if (i == SolnBlk.ICl-1) {
              SolnBlk.FluxW[j] = -flux*SolnBlk.Grid.lfaceW(i+1, j);
           } else if (i == SolnBlk.ICu) {
              SolnBlk.FluxE[j] = flux*SolnBlk.Grid.lfaceE(i, j);
           } /* endif */

         } /* endif */
       } /* endfor */
    
       if ( j > SolnBlk.JCl-1 && j < SolnBlk.JCu+1 ) {
          SolnBlk.dudt[SolnBlk.ICl-1][j][0] = ZERO;
          SolnBlk.dudt[SolnBlk.ICu+1][j][0] = ZERO;
       } /* endif */
    } /* endfor */
    
    // Add j-direction (eta-direction) fluxes.
    for ( i = SolnBlk.ICl ; i <= SolnBlk.ICu ; ++i ) {
       for ( j  = SolnBlk.JCl-1 ; j <= SolnBlk.JCu ; ++j ) {
    
          /* Evaluate the cell interface j-direction fluxes. */
         
  	  dX = SolnBlk.Grid.xfaceN(i, j)-SolnBlk.Grid.Cell[i][j].Xc;
          ul = SolnBlk.U[i][j].u + 
               (SolnBlk.phi[i][j]*SolnBlk.dudx[i][j])*dX.x +
               (SolnBlk.phi[i][j]*SolnBlk.dudy[i][j])*dX.y;
	  dX = SolnBlk.Grid.xfaceS(i, j+1)-SolnBlk.Grid.Cell[i][j+1].Xc;
          ur = SolnBlk.U[i][j+1].u +
               (SolnBlk.phi[i][j+1]*SolnBlk.dudx[i][j+1])*dX.x +
               (SolnBlk.phi[i][j+1]*SolnBlk.dudy[i][j+1])*dX.y;

          Ul = SolnBlk.U[i][j]; 
          Ul.u = ul;
          Ur = SolnBlk.U[i][j+1]; 
          Ur.u = ur;

	  if ((j == SolnBlk.JCl-1) && 
              (SolnBlk.Grid.BCtypeS[i] == BC_DIRICHLET)) {
             Ul.Fd = Ur.Fd;
          } else if ((j == SolnBlk.JCu) && 
                     (SolnBlk.Grid.BCtypeN[i] == BC_DIRICHLET)) {
             Ur.Fd = Ul.Fd;
  	  } /* endif */

          flux = Flux_n(Ul, Ur, SolnBlk.Grid.nfaceN(i, j));
    
          /* Evaluate cell-averaged solution changes. */
    
          SolnBlk.dudt[i][j][0] -=
             flux*SolnBlk.Grid.lfaceN(i, j)/
             SolnBlk.Grid.Cell[i][j].A;
          SolnBlk.dudt[i][j+1][0] += 
             flux*SolnBlk.Grid.lfaceS(i, j+1)/
             SolnBlk.Grid.Cell[i][j+1].A;
          
          /* Save south and north face boundary flux. */

          if (j == SolnBlk.JCl-1) {
             SolnBlk.FluxS[i] = -flux*SolnBlk.Grid.lfaceS(i, j+1);
          } else if (j == SolnBlk.JCu) {
             SolnBlk.FluxN[i] = flux*SolnBlk.Grid.lfaceN(i, j);
          } /* endif */

       } /* endfor */

       SolnBlk.dudt[i][SolnBlk.JCl-1][0] = ZERO;
       SolnBlk.dudt[i][SolnBlk.JCu+1][0] = ZERO;
    } /* endfor */
    
    /* residual evaluation successful. */
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
int dUdt_Multistage_Explicit(AdvectDiffuse2D_Quad_Block &SolnBlk,
                             const int i_stage,
                             AdvectDiffuse2D_Input_Parameters &Input_Parameters) {

    int i, j, k_residual;
    double omega;
    Vector2D dX;
    double ul, ur, flux;
    AdvectDiffuse2D_State Ul, Ur;

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

    /* Evaluate the diffusive fluxes in each cell
       of the computational grid for this stage. */

    Diffusive_Flux(SolnBlk);

    /* Evaluate the time rate of change of the solution
       (i.e., the solution residuals) using a second-order
       limited upwind finite-volume scheme for the convective 
       fluxes and a second-order centrally-weighted finite 
       volume discretization for the diffusive fluxes. */
    
    // Add i-direction (zeta-direction) fluxes.
    for ( j  = SolnBlk.JCl-1 ; j <= SolnBlk.JCu+1 ; ++j ) {
       if ( i_stage == 1 ) {
          SolnBlk.uo[SolnBlk.ICl-1][j] = SolnBlk.U[SolnBlk.ICl-1][j].u;
          SolnBlk.dudt[SolnBlk.ICl-1][j][k_residual] = ZERO;
       } else {
          SolnBlk.dudt[SolnBlk.ICl-1][j][k_residual] = ZERO;
       } /* endif */
    
       for ( i = SolnBlk.ICl-1 ; i <= SolnBlk.ICu ; ++i ) {
          if ( i_stage == 1 ) {
              SolnBlk.uo[i+1][j] = SolnBlk.U[i+1][j].u;
              SolnBlk.dudt[i+1][j][k_residual] = ZERO;
          } else if ( j > SolnBlk.JCl-1 && j < SolnBlk.JCu+1 ) {
              switch(Input_Parameters.i_Time_Integration) {
                case TIME_STEPPING_EXPLICIT_PREDICTOR_CORRECTOR :
                  //SolnBlk.dudt[i+1][j][k_residual] = 
                  //   SolnBlk.dudt[i+1][j][k_residual];
                  break;
                case TIME_STEPPING_EXPLICIT_RUNGE_KUTTA :
                  if (Input_Parameters.N_Stage == 2) {
		     //SolnBlk.dudt[i+1][j][k_residual] = 
		     //   SolnBlk.dudt[i+1][j][k_residual];
                  } else if (Input_Parameters.N_Stage == 4 && i_stage == 4) {
                     SolnBlk.dudt[i+1][j][k_residual] = 
                        SolnBlk.dudt[i+1][j][0] + 
                        TWO*SolnBlk.dudt[i+1][j][1] +
                        TWO*SolnBlk.dudt[i+1][j][2];
                  } else {
                     SolnBlk.dudt[i+1][j][k_residual] = ZERO;
                  } /* endif */
                  break;
                case TIME_STEPPING_MULTISTAGE_OPTIMAL_SMOOTHING :
                  SolnBlk.dudt[i+1][j][k_residual] = ZERO;
                  break;
                default:
                  SolnBlk.dudt[i+1][j][k_residual] = ZERO;
                  break;
              } /* endswitch */
          } /* endif */
    
          if ( j > SolnBlk.JCl-1 && j < SolnBlk.JCu+1 ) {
    
             /* Evaluate the cell interface i-direction fluxes. */
    
  	     dX = SolnBlk.Grid.xfaceE(i, j)-SolnBlk.Grid.Cell[i][j].Xc;
             ul = SolnBlk.U[i][j].u + 
                  (SolnBlk.phi[i][j]*SolnBlk.dudx[i][j])*dX.x +
                  (SolnBlk.phi[i][j]*SolnBlk.dudy[i][j])*dX.y;
	     dX = SolnBlk.Grid.xfaceW(i+1, j)-SolnBlk.Grid.Cell[i+1][j].Xc;
             ur = SolnBlk.U[i+1][j].u + 
                  (SolnBlk.phi[i+1][j]*SolnBlk.dudx[i+1][j])*dX.x +
                  (SolnBlk.phi[i+1][j]*SolnBlk.dudy[i+1][j])*dX.y;

             Ul = SolnBlk.U[i][j]; 
             Ul.u = ul;
             Ur = SolnBlk.U[i+1][j]; 
             Ur.u = ur;

	     if ((i == SolnBlk.ICl-1) && 
                 (SolnBlk.Grid.BCtypeW[j] == BC_DIRICHLET)) {
                Ul.Fd = Ur.Fd;
             } else if ((i == SolnBlk.ICu) && 
                        (SolnBlk.Grid.BCtypeE[j] == BC_DIRICHLET)) {
                Ur.Fd = Ul.Fd;
  	     } /* endif */

             flux = Flux_n(Ul, Ur, SolnBlk.Grid.nfaceE(i, j));
    
             /* Evaluate cell-averaged solution changes. */
    
             SolnBlk.dudt[i][j][k_residual] -= 
                (Input_Parameters.CFL_Number*SolnBlk.dt[i][j])*
                flux*SolnBlk.Grid.lfaceE(i, j)/
                SolnBlk.Grid.Cell[i][j].A;
             SolnBlk.dudt[i+1][j][k_residual] +=
                (Input_Parameters.CFL_Number*SolnBlk.dt[i+1][j])*
                flux*SolnBlk.Grid.lfaceW(i+1, j)/
                SolnBlk.Grid.Cell[i+1][j].A;

             /* Include regular source terms. */

             SolnBlk.dudt[i][j][k_residual] += 
                (Input_Parameters.CFL_Number*SolnBlk.dt[i][j])*s(SolnBlk.U[i][j]);

             /* Include axisymmetric source terms as required. */

	     if (SolnBlk.Axisymmetric) {
               SolnBlk.dudt[i][j][k_residual] += 
                  (Input_Parameters.CFL_Number*SolnBlk.dt[i][j])*
	          s_axi(SolnBlk.U[i][j], SolnBlk.Grid.Cell[i][j].Xc);
             } /* endif */

             /* Save west and east face boundary flux. */

             if (i == SolnBlk.ICl-1) {
                SolnBlk.FluxW[j] = -flux*SolnBlk.Grid.lfaceW(i+1, j);
             } else if (i == SolnBlk.ICu) {
                SolnBlk.FluxE[j] = flux*SolnBlk.Grid.lfaceE(i, j);
             } /* endif */

          } /* endif */
       } /* endfor */
    
       if ( j > SolnBlk.JCl-1 && j < SolnBlk.JCu+1 ) {
          SolnBlk.dudt[SolnBlk.ICl-1][j][k_residual] = ZERO;
          SolnBlk.dudt[SolnBlk.ICu+1][j][k_residual] = ZERO;
       } /* endif */
    } /* endfor */
    
    // Add j-direction (eta-direction) fluxes.
    for ( i = SolnBlk.ICl ; i <= SolnBlk.ICu ; ++i ) {
       for ( j  = SolnBlk.JCl-1 ; j <= SolnBlk.JCu ; ++j ) {
    
          /* Evaluate the cell interface j-direction fluxes. */
         
  	  dX = SolnBlk.Grid.xfaceN(i, j)-SolnBlk.Grid.Cell[i][j].Xc;
          ul = SolnBlk.U[i][j].u + 
               (SolnBlk.phi[i][j]*SolnBlk.dudx[i][j])*dX.x +
               (SolnBlk.phi[i][j]*SolnBlk.dudy[i][j])*dX.y;
	  dX = SolnBlk.Grid.xfaceS(i, j+1)-SolnBlk.Grid.Cell[i][j+1].Xc;
          ur = SolnBlk.U[i][j+1].u +
               (SolnBlk.phi[i][j+1]*SolnBlk.dudx[i][j+1])*dX.x +
               (SolnBlk.phi[i][j+1]*SolnBlk.dudy[i][j+1])*dX.y;

          Ul = SolnBlk.U[i][j]; 
          Ul.u = ul;
          Ur = SolnBlk.U[i][j+1]; 
          Ur.u = ur;

	  if ((j == SolnBlk.JCl-1) && 
              (SolnBlk.Grid.BCtypeS[i] == BC_DIRICHLET)) {
             Ul.Fd = Ur.Fd;
          } else if ((j == SolnBlk.JCu) && 
                     (SolnBlk.Grid.BCtypeN[i] == BC_DIRICHLET)) {
             Ur.Fd = Ul.Fd;
  	  } /* endif */

          flux = Flux_n(Ul, Ur, SolnBlk.Grid.nfaceN(i, j));
    
          /* Evaluate cell-averaged solution changes. */
    
          SolnBlk.dudt[i][j][k_residual] -=
             (Input_Parameters.CFL_Number*SolnBlk.dt[i][j])*
             flux*SolnBlk.Grid.lfaceN(i, j)/
             SolnBlk.Grid.Cell[i][j].A;
          SolnBlk.dudt[i][j+1][k_residual] += 
             (Input_Parameters.CFL_Number*SolnBlk.dt[i][j+1])*
             flux*SolnBlk.Grid.lfaceS(i, j+1)/
             SolnBlk.Grid.Cell[i][j+1].A;
          
          /* Save south and north face boundary flux. */

          if (j == SolnBlk.JCl-1) {
             SolnBlk.FluxS[i] = -flux*SolnBlk.Grid.lfaceS(i, j+1);
          } else if (j == SolnBlk.JCu) {
             SolnBlk.FluxN[i] = flux*SolnBlk.Grid.lfaceN(i, j);
          } /* endif */

       } /* endfor */

       SolnBlk.dudt[i][SolnBlk.JCl-1][k_residual] = ZERO;
       SolnBlk.dudt[i][SolnBlk.JCu+1][k_residual] = ZERO;
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
int Update_Solution_Multistage_Explicit(AdvectDiffuse2D_Quad_Block &SolnBlk,
                                        const int i_stage,
                                        AdvectDiffuse2D_Input_Parameters &Input_Parameters) {

    int i, j, k_residual;
    double omega;

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
          SolnBlk.U[i][j].u = SolnBlk.uo[i][j] + 
                              omega*SolnBlk.dudt[i][j][k_residual];
       } /* endfor */    
    } /* endfor */

    /* Solution successfully updated. */

    return (0);
    
}
