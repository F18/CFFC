/* Ion5Moment2DQuadSingleBlock.cc:  Single-Block Versions of Subroutines for
                                    2D 5-Moment Ion Transport Model
                                    Multi-Block Quadrilateral Mesh 
                                    Solution Classes. */

/* Include 2D 5-moment ion transport model quadrilateral mesh solution header file. */

#ifndef _ION5MOMENT2D_QUAD_INCLUDED
#include "Ion5Moment2DQuad.h"
#endif // _ION5MOMENT2D_QUAD_INCLUDED

/**************************************************************************
 * Ion5Moment2D_Quad_Block -- Single Block External Subroutines.          *
 **************************************************************************/

/********************************************************
 * Routine: Write_Solution_Block                        *
 *                                                      *
 * Writes the cell centred solution values of the       *
 * specified quadrilateral solution block to the        *
 * specified output stream for restart purposes.        *
 *                                                      *
 ********************************************************/
void Write_Solution_Block(Ion5Moment2D_Quad_Block &SolnBlk,
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
void Read_Solution_Block(Ion5Moment2D_Quad_Block &SolnBlk,
	                 istream &In_File) {

    In_File >> SolnBlk;

}

/********************************************************
 * Routine: Read_Neutral_Gas_Solution_Block             *
 *                                                      *
 * Reads the cell centred neutral gas flow solution     *
 * values for the specified quadrilateral solution      *
 * block from the specified input stream.  The netural  *
 * flow solution file is expected to be in ASCII        *
 * TECPLOT point data format.                           *
 *                                                      *
 ********************************************************/
void Read_Neutral_Gas_Solution_Block(Ion5Moment2D_Quad_Block &SolnBlk,
	                             istream &In_File) {

    int i, j;
    double XX, YY, ZZ, PP, TT, UU, VV, WW, RR, MM;

    In_File.setf(ios::skipws);

    /* Read and assign the neutral gas solution states. */

    for (j  = SolnBlk.JCl; j <= SolnBlk.JCu; ++j ) {
        for ( i = SolnBlk.ICl ; i <= SolnBlk.ICu ; ++i ) {
            In_File >> XX >> YY >> ZZ >> PP >> TT >> UU >> VV >> WW >> RR >> MM;
            SolnBlk.Wneut[i][j].d = RR;
            SolnBlk.Wneut[i][j].v.x = UU;
            SolnBlk.Wneut[i][j].v.y = VV;
            SolnBlk.Wneut[i][j].p = PP;

            // Set the ion velocity and temperature to the neutral gas values.
            SolnBlk.W[i][j].v = SolnBlk.Wneut[i][j].v;
            SolnBlk.W[i][j].p = SolnBlk.W[i][j].d*SolnBlk.W[i][j].R*SolnBlk.Wneut[i][j].T();
            SolnBlk.U[i][j] = U(SolnBlk.W[i][j]);
        } /* endfor */
    } /* endfor */

    In_File.unsetf(ios::skipws);

    /* Reset values for the ion boundary condition reference states. */

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
       SolnBlk.WoW[j].v = Vector2D_ZERO;
       SolnBlk.WoE[j].v = Vector2D_ZERO;
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
       SolnBlk.WoS[i].v = Vector2D_ZERO;
       SolnBlk.WoN[i].v = Vector2D_ZERO;
    } /* endfor */

}

/********************************************************
 * Routine: Read_Electric_Field_Solution_Block          *
 *                                                      *
 * Reads the cell centred electric field solution       *
 * values for the specified quadrilateral solution      *
 * block from the specified input stream.  The electric *
 * field solution file is expected to be in ASCII       *
 * TECPLOT point data format.                           *
 *                                                      *
 ********************************************************/
void Read_Electric_Field_Solution_Block(Ion5Moment2D_Quad_Block &SolnBlk,
	                                istream &In_File,
                                        const int Add_Initial_and_Solution_File_Electric_Fields) {

    int i, j;
    double XX, YY, ZZ, EX, EY, EZ, VV;

    In_File.setf(ios::skipws);

    /* Read and assign the electric field solution values. */

    for (j  = SolnBlk.JCl; j <= SolnBlk.JCu; ++j ) {
        for ( i = SolnBlk.ICl ; i <= SolnBlk.ICu ; ++i ) {
            In_File >> XX >> YY >> ZZ >> EX >> EY >> EZ >> VV;
            if (Add_Initial_and_Solution_File_Electric_Fields) {
               SolnBlk.E[i][j].x = SolnBlk.E[i][j].x + EX;
               SolnBlk.E[i][j].y = SolnBlk.E[i][j].y + EY;
               SolnBlk.V[i][j] = SolnBlk.V[i][j] + VV;
            } else {
               SolnBlk.E[i][j].x = EX;
               SolnBlk.E[i][j].y = EY;
               SolnBlk.V[i][j] = VV;
            } /* endif */
        } /* endfor */
    } /* endfor */

    In_File.unsetf(ios::skipws);

}

/********************************************************
 * Routine: Broadcast_Solution_Block                    *
 *                                                      *
 * Broadcast quadrilateral solution block to all        *
 * processors involved in the calculation from the      *
 * primary processor using the MPI broadcast routine.   *
 *                                                      *
 ********************************************************/
void Broadcast_Solution_Block(Ion5Moment2D_Quad_Block &SolnBlk) {

#ifdef _MPI_VERSION
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
    MPI::COMM_WORLD.Bcast(&nr, 1, MPI::INT, 0);
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

    /* Broadcast the grid. */

    Broadcast_Quad_Block(SolnBlk.Grid);

    /* Broadcast the solution state variables. */

    if (block_allocated) {
       ni = (SolnBlk.ICu+SolnBlk.Nghost) - (SolnBlk.ICl-SolnBlk.Nghost) + 1;
       nj = (SolnBlk.JCu+SolnBlk.Nghost) - (SolnBlk.JCl-SolnBlk.Nghost) + 1;
       buffer = new double[11*ni*nj];

       if (CFFC_Primary_MPI_Processor()) {
          buffer_size = 0;
          for (j  = SolnBlk.JCl-SolnBlk.Nghost ; j <= SolnBlk.JCu+SolnBlk.Nghost ; ++j ) {
              for ( i = SolnBlk.ICl-SolnBlk.Nghost ; i <= SolnBlk.ICu+SolnBlk.Nghost ; ++i ) {
 	          buffer[buffer_size   ] = SolnBlk.U[i][j].d;
 	          buffer[buffer_size+ 1] = SolnBlk.U[i][j].dv.x;
 	          buffer[buffer_size+ 2] = SolnBlk.U[i][j].dv.y;
 	          buffer[buffer_size+ 3] = SolnBlk.U[i][j].E;
 	          buffer[buffer_size+ 4] = SolnBlk.Wneut[i][j].d;
 	          buffer[buffer_size+ 5] = SolnBlk.Wneut[i][j].v.x;
 	          buffer[buffer_size+ 6] = SolnBlk.Wneut[i][j].v.y;
 	          buffer[buffer_size+ 7] = SolnBlk.Wneut[i][j].p;
 	          buffer[buffer_size+ 8] = SolnBlk.E[i][j].x;
 	          buffer[buffer_size+ 9] = SolnBlk.E[i][j].y;
 	          buffer[buffer_size+10] = SolnBlk.V[i][j];
                  buffer_size = buffer_size + 11;
              } /* endfor */
          } /* endfor */
       } /* endif */

       buffer_size = 11*ni*nj;
       MPI::COMM_WORLD.Bcast(buffer, buffer_size, MPI::DOUBLE, 0);

       if (!CFFC_Primary_MPI_Processor()) {
          buffer_size = 0;
          for (j  = SolnBlk.JCl-SolnBlk.Nghost ; j <= SolnBlk.JCu+SolnBlk.Nghost ; ++j ) {
              for ( i = SolnBlk.ICl-SolnBlk.Nghost ; i <= SolnBlk.ICu+SolnBlk.Nghost ; ++i ) {
 	          SolnBlk.U[i][j].d       = buffer[buffer_size];
 	          SolnBlk.U[i][j].dv.x    = buffer[buffer_size+1];
 	          SolnBlk.U[i][j].dv.y    = buffer[buffer_size+2];
 	          SolnBlk.U[i][j].E       = buffer[buffer_size+3];
 	          SolnBlk.Wneut[i][j].d   = buffer[buffer_size+ 4];
 	          SolnBlk.Wneut[i][j].v.x = buffer[buffer_size+ 5];
 	          SolnBlk.Wneut[i][j].v.y = buffer[buffer_size+ 6];
 	          SolnBlk.Wneut[i][j].p   = buffer[buffer_size+ 7];
 	          SolnBlk.E[i][j].x       = buffer[buffer_size+ 8];
 	          SolnBlk.E[i][j].y       = buffer[buffer_size+ 9];
 	          SolnBlk.V[i][j]         = buffer[buffer_size+10];
                  buffer_size = buffer_size + 11;
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
void Broadcast_Solution_Block(Ion5Moment2D_Quad_Block &SolnBlk,
                              MPI::Intracomm &Communicator, 
                              const int Source_CPU) {

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
    Communicator.Bcast(&nr, 1, MPI::INT, Source_Rank);
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

    /* Broadcast the grid. */

    Broadcast_Quad_Block(SolnBlk.Grid, Communicator, Source_CPU);

    /* Broadcast the solution state variables. */

    if (block_allocated) {
       ni = (SolnBlk.ICu+SolnBlk.Nghost) - (SolnBlk.ICl-SolnBlk.Nghost) + 1;
       nj = (SolnBlk.JCu+SolnBlk.Nghost) - (SolnBlk.JCl-SolnBlk.Nghost) + 1;
       buffer = new double[11*ni*nj];

       if (CFFC_MPI::This_Processor_Number == Source_CPU) {
          buffer_size = 0;
          for (j  = SolnBlk.JCl-SolnBlk.Nghost ; j <= SolnBlk.JCu+SolnBlk.Nghost ; ++j ) {
              for ( i = SolnBlk.ICl-SolnBlk.Nghost ; i <= SolnBlk.ICu+SolnBlk.Nghost ; ++i ) {
 	          buffer[buffer_size   ] = SolnBlk.U[i][j].d;
 	          buffer[buffer_size+ 1] = SolnBlk.U[i][j].dv.x;
 	          buffer[buffer_size+ 2] = SolnBlk.U[i][j].dv.y;
 	          buffer[buffer_size+ 3] = SolnBlk.U[i][j].E;
 	          buffer[buffer_size+ 4] = SolnBlk.Wneut[i][j].d;
 	          buffer[buffer_size+ 5] = SolnBlk.Wneut[i][j].v.x;
 	          buffer[buffer_size+ 6] = SolnBlk.Wneut[i][j].v.y;
 	          buffer[buffer_size+ 7] = SolnBlk.Wneut[i][j].p;
 	          buffer[buffer_size+ 8] = SolnBlk.E[i][j].x;
 	          buffer[buffer_size+ 9] = SolnBlk.E[i][j].y;
 	          buffer[buffer_size+10] = SolnBlk.V[i][j];
                  buffer_size = buffer_size + 11;
              } /* endfor */
          } /* endfor */
       } /* endif */

       buffer_size = 11*ni*nj;
       Communicator.Bcast(buffer, buffer_size, MPI::DOUBLE, Source_Rank);

       if (!(CFFC_MPI::This_Processor_Number == Source_CPU)) {
          buffer_size = 0;
          for (j  = SolnBlk.JCl-SolnBlk.Nghost ; j <= SolnBlk.JCu+SolnBlk.Nghost ; ++j ) {
              for ( i = SolnBlk.ICl-SolnBlk.Nghost ; i <= SolnBlk.ICu+SolnBlk.Nghost ; ++i ) {
 	          SolnBlk.U[i][j].d       = buffer[buffer_size];
 	          SolnBlk.U[i][j].dv.x    = buffer[buffer_size+1];
 	          SolnBlk.U[i][j].dv.y    = buffer[buffer_size+2];
 	          SolnBlk.U[i][j].E       = buffer[buffer_size+3];
 	          SolnBlk.Wneut[i][j].d   = buffer[buffer_size+ 4];
 	          SolnBlk.Wneut[i][j].v.x = buffer[buffer_size+ 5];
 	          SolnBlk.Wneut[i][j].v.y = buffer[buffer_size+ 6];
 	          SolnBlk.Wneut[i][j].p   = buffer[buffer_size+ 7];
 	          SolnBlk.E[i][j].x       = buffer[buffer_size+ 8];
 	          SolnBlk.E[i][j].y       = buffer[buffer_size+ 9];
 	          SolnBlk.V[i][j]         = buffer[buffer_size+10];
                  buffer_size = buffer_size + 11;
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
void Copy_Solution_Block(Ion5Moment2D_Quad_Block &SolnBlk1,
                         Ion5Moment2D_Quad_Block &SolnBlk2) {

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
             SolnBlk1.Wneut[i][j] = SolnBlk2.Wneut[i][j];
             SolnBlk1.E[i][j] = SolnBlk2.E[i][j];
             SolnBlk1.V[i][j] = SolnBlk2.V[i][j];
             for ( k = 0 ; k <= NUMBER_OF_RESIDUAL_VECTORS_ION5MOMENT2D-1 ; ++k ) {
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
int Prolong_Solution_Block(Ion5Moment2D_Quad_Block &SolnBlk_Fine,
			   Ion5Moment2D_Quad_Block &SolnBlk_Original,
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
    	       SolnBlk_Fine.Wneut[2*(i-i_min)+SolnBlk_Fine.ICl  ]
                                 [2*(j-j_min)+SolnBlk_Fine.JCl  ] 
                  = SolnBlk_Original.Wneut[i][j];
//                   = W((SolnBlk_Original.Grid.Cell[i][j].A/area_total_fine)*
//                       (SolnBlk_Original.Wneut[i][j].U()));
    	       SolnBlk_Fine.E[2*(i-i_min)+SolnBlk_Fine.ICl  ]
                             [2*(j-j_min)+SolnBlk_Fine.JCl  ] 
                  = SolnBlk_Original.E[i][j];
    	       SolnBlk_Fine.V[2*(i-i_min)+SolnBlk_Fine.ICl  ]
                             [2*(j-j_min)+SolnBlk_Fine.JCl  ] 
                  = SolnBlk_Original.V[i][j];

    	       SolnBlk_Fine.U[2*(i-i_min)+SolnBlk_Fine.ICl+1]
                             [2*(j-j_min)+SolnBlk_Fine.JCl  ] 
                  = SolnBlk_Original.U[i][j];
//                   = (SolnBlk_Original.Grid.Cell[i][j].A/area_total_fine)*SolnBlk_Original.U[i][j];
    	       SolnBlk_Fine.W[2*(i-i_min)+SolnBlk_Fine.ICl+1]
                             [2*(j-j_min)+SolnBlk_Fine.JCl  ] 
                  = W(SolnBlk_Fine.U[2*(i-i_min)+SolnBlk_Fine.ICl+1]
                                    [2*(j-j_min)+SolnBlk_Fine.JCl  ]);
    	       SolnBlk_Fine.Wneut[2*(i-i_min)+SolnBlk_Fine.ICl+1]
                                 [2*(j-j_min)+SolnBlk_Fine.JCl  ] 
                  = SolnBlk_Original.Wneut[i][j];
//                   = W((SolnBlk_Original.Grid.Cell[i][j].A/area_total_fine)*
//                       (SolnBlk_Original.Wneut[i][j].U()));
    	       SolnBlk_Fine.E[2*(i-i_min)+SolnBlk_Fine.ICl+1]
                             [2*(j-j_min)+SolnBlk_Fine.JCl  ] 
                  = SolnBlk_Original.E[i][j];
    	       SolnBlk_Fine.V[2*(i-i_min)+SolnBlk_Fine.ICl+1]
                             [2*(j-j_min)+SolnBlk_Fine.JCl  ] 
                  = SolnBlk_Original.V[i][j];

    	       SolnBlk_Fine.U[2*(i-i_min)+SolnBlk_Fine.ICl  ]
                             [2*(j-j_min)+SolnBlk_Fine.JCl+1]
                  = SolnBlk_Original.U[i][j];
//                   = (SolnBlk_Original.Grid.Cell[i][j].A/area_total_fine)*SolnBlk_Original.U[i][j];
    	       SolnBlk_Fine.W[2*(i-i_min)+SolnBlk_Fine.ICl  ]
                             [2*(j-j_min)+SolnBlk_Fine.JCl+1]
                  = W(SolnBlk_Fine.U[2*(i-i_min)+SolnBlk_Fine.ICl  ]
                                    [2*(j-j_min)+SolnBlk_Fine.JCl+1]);
    	       SolnBlk_Fine.Wneut[2*(i-i_min)+SolnBlk_Fine.ICl  ]
                                 [2*(j-j_min)+SolnBlk_Fine.JCl+1]
                  = SolnBlk_Original.Wneut[i][j];
//                   = W((SolnBlk_Original.Grid.Cell[i][j].A/area_total_fine)*
//                       (SolnBlk_Original.Wneut[i][j].U()));
    	       SolnBlk_Fine.E[2*(i-i_min)+SolnBlk_Fine.ICl  ]
                             [2*(j-j_min)+SolnBlk_Fine.JCl+1]
                  = SolnBlk_Original.E[i][j];
    	       SolnBlk_Fine.V[2*(i-i_min)+SolnBlk_Fine.ICl  ]
                             [2*(j-j_min)+SolnBlk_Fine.JCl+1]
                  = SolnBlk_Original.V[i][j];

    	       SolnBlk_Fine.U[2*(i-i_min)+SolnBlk_Fine.ICl+1]
                             [2*(j-j_min)+SolnBlk_Fine.JCl+1]
                  = (SolnBlk_Original.Grid.Cell[i][j].A/area_total_fine)*SolnBlk_Original.U[i][j];
//                   = (SolnBlk_Original.Grid.Cell[i][j].A/area_total_fine)*SolnBlk_Original.U[i][j];
    	       SolnBlk_Fine.W[2*(i-i_min)+SolnBlk_Fine.ICl+1]
                             [2*(j-j_min)+SolnBlk_Fine.JCl+1]
                  = W(SolnBlk_Fine.U[2*(i-i_min)+SolnBlk_Fine.ICl+1]
                                    [2*(j-j_min)+SolnBlk_Fine.JCl+1]);
    	       SolnBlk_Fine.Wneut[2*(i-i_min)+SolnBlk_Fine.ICl+1]
                                 [2*(j-j_min)+SolnBlk_Fine.JCl+1]
                  = SolnBlk_Original.Wneut[i][j];
//                   = W((SolnBlk_Original.Grid.Cell[i][j].A/area_total_fine)*
//                       (SolnBlk_Original.Wneut[i][j].U()));
    	       SolnBlk_Fine.E[2*(i-i_min)+SolnBlk_Fine.ICl+1]
                             [2*(j-j_min)+SolnBlk_Fine.JCl+1]
                  = SolnBlk_Original.E[i][j];
    	       SolnBlk_Fine.V[2*(i-i_min)+SolnBlk_Fine.ICl+1]
                             [2*(j-j_min)+SolnBlk_Fine.JCl+1]
                  = SolnBlk_Original.V[i][j];
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
int Restrict_Solution_Block(Ion5Moment2D_Quad_Block &SolnBlk_Coarse,
			    Ion5Moment2D_Quad_Block &SolnBlk_Original_SW,
			    Ion5Moment2D_Quad_Block &SolnBlk_Original_SE,
			    Ion5Moment2D_Quad_Block &SolnBlk_Original_NW,
			    Ion5Moment2D_Quad_Block &SolnBlk_Original_NE) {

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
             SolnBlk_Coarse.W[i_coarse][j_coarse] = W(SolnBlk_Coarse.U[i_coarse][j_coarse]);
             SolnBlk_Coarse.Wneut[i_coarse][j_coarse] = 
                                                  W((SolnBlk_Original_SW.Grid.Cell[i  ][j  ].A*
                                                     SolnBlk_Original_SW.Wneut[i  ][j  ].U() +
                                                     SolnBlk_Original_SW.Grid.Cell[i+1][j  ].A*
                                                     SolnBlk_Original_SW.Wneut[i+1][j  ].U() + 
                                                     SolnBlk_Original_SW.Grid.Cell[i  ][j+1].A*
                                                     SolnBlk_Original_SW.Wneut[i  ][j+1].U() +
                                                     SolnBlk_Original_SW.Grid.Cell[i+1][j+1].A*
                                                     SolnBlk_Original_SW.Wneut[i+1][j+1].U()) /
                                                    (SolnBlk_Original_SW.Grid.Cell[i  ][j  ].A +
                                                     SolnBlk_Original_SW.Grid.Cell[i+1][j  ].A +
                                                     SolnBlk_Original_SW.Grid.Cell[i  ][j+1].A +
                                                     SolnBlk_Original_SW.Grid.Cell[i+1][j+1].A));
//                                                      SolnBlk_Coarse.Grid.Cell[i_coarse][j_coarse].A);
             SolnBlk_Coarse.E[i_coarse][j_coarse] = SolnBlk_Original_SW.En(i+1,j+1);
             SolnBlk_Coarse.V[i_coarse][j_coarse] = SolnBlk_Original_SW.Vn(i+1,j+1);
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
             SolnBlk_Coarse.Wneut[i_coarse][j_coarse] = 
                                                  W((SolnBlk_Original_SE.Grid.Cell[i  ][j  ].A*
                                                     SolnBlk_Original_SE.Wneut[i  ][j  ].U() +
                                                     SolnBlk_Original_SE.Grid.Cell[i+1][j  ].A*
                                                     SolnBlk_Original_SE.Wneut[i+1][j  ].U() + 
                                                     SolnBlk_Original_SE.Grid.Cell[i  ][j+1].A*
                                                     SolnBlk_Original_SE.Wneut[i  ][j+1].U() +
                                                     SolnBlk_Original_SE.Grid.Cell[i+1][j+1].A*
                                                     SolnBlk_Original_SE.Wneut[i+1][j+1].U()) /
                                                    (SolnBlk_Original_SE.Grid.Cell[i  ][j  ].A +
                                                     SolnBlk_Original_SE.Grid.Cell[i+1][j  ].A +
                                                     SolnBlk_Original_SE.Grid.Cell[i  ][j+1].A +
                                                     SolnBlk_Original_SE.Grid.Cell[i+1][j+1].A));
//                                                      SolnBlk_Coarse.Grid.Cell[i_coarse][j_coarse].A);
             SolnBlk_Coarse.E[i_coarse][j_coarse] = SolnBlk_Original_SE.En(i+1,j+1);
             SolnBlk_Coarse.V[i_coarse][j_coarse] = SolnBlk_Original_SE.Vn(i+1,j+1);
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
             SolnBlk_Coarse.Wneut[i_coarse][j_coarse] = 
                                                  W((SolnBlk_Original_NW.Grid.Cell[i  ][j  ].A*
                                                     SolnBlk_Original_NW.Wneut[i  ][j  ].U() +
                                                     SolnBlk_Original_NW.Grid.Cell[i+1][j  ].A*
                                                     SolnBlk_Original_NW.Wneut[i+1][j  ].U() + 
                                                     SolnBlk_Original_NW.Grid.Cell[i  ][j+1].A*
                                                     SolnBlk_Original_NW.Wneut[i  ][j+1].U() +
                                                     SolnBlk_Original_NW.Grid.Cell[i+1][j+1].A*
                                                     SolnBlk_Original_NW.Wneut[i+1][j+1].U()) /
                                                    (SolnBlk_Original_NW.Grid.Cell[i  ][j  ].A +
                                                     SolnBlk_Original_NW.Grid.Cell[i+1][j  ].A +
                                                     SolnBlk_Original_NW.Grid.Cell[i  ][j+1].A +
                                                     SolnBlk_Original_NW.Grid.Cell[i+1][j+1].A));
//                                                      SolnBlk_Coarse.Grid.Cell[i_coarse][j_coarse].A);
             SolnBlk_Coarse.E[i_coarse][j_coarse] = SolnBlk_Original_NW.En(i+1,j+1);
             SolnBlk_Coarse.V[i_coarse][j_coarse] = SolnBlk_Original_NW.Vn(i+1,j+1);
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
             SolnBlk_Coarse.Wneut[i_coarse][j_coarse] = 
                                                  W((SolnBlk_Original_NE.Grid.Cell[i  ][j  ].A*
                                                     SolnBlk_Original_NE.Wneut[i  ][j  ].U() +
                                                     SolnBlk_Original_NE.Grid.Cell[i+1][j  ].A*
                                                     SolnBlk_Original_NE.Wneut[i+1][j  ].U() + 
                                                     SolnBlk_Original_NE.Grid.Cell[i  ][j+1].A*
                                                     SolnBlk_Original_NE.Wneut[i  ][j+1].U() +
                                                     SolnBlk_Original_NE.Grid.Cell[i+1][j+1].A*
                                                     SolnBlk_Original_NE.Wneut[i+1][j+1].U()) /
                                                    (SolnBlk_Original_NE.Grid.Cell[i  ][j  ].A +
                                                     SolnBlk_Original_NE.Grid.Cell[i+1][j  ].A +
                                                     SolnBlk_Original_NE.Grid.Cell[i  ][j+1].A +
                                                     SolnBlk_Original_NE.Grid.Cell[i+1][j+1].A));
//                                                      SolnBlk_Coarse.Grid.Cell[i_coarse][j_coarse].A);
             SolnBlk_Coarse.E[i_coarse][j_coarse] = SolnBlk_Original_NE.En(i+1,j+1);
             SolnBlk_Coarse.V[i_coarse][j_coarse] = SolnBlk_Original_NE.Vn(i+1,j+1);
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
void Output_Tecplot(Ion5Moment2D_Quad_Block &SolnBlk,
                    Ion5Moment2D_Input_Parameters &Input_Parameters,
                    const int Number_of_Time_Steps,
                    const double &Time,
                    const int Block_Number,
                    const int Output_Title,
	            ostream &Out_File) {

    int i, j;
    Ion5Moment2D_pState W_node;
    Euler2D_pState Wneut_node;
    Vector2D E_node, J_node;
    double V_node;

    /* Ensure boundary conditions are updated before
       evaluating solution at the nodes. */
    
    BCs(SolnBlk);

    /* Output node solution data. */

    Out_File << setprecision(14);
    if (Output_Title) {
       Out_File << "TITLE = \"" << CFFC_Name() << ": 2D 5-Moment Ion Transport Solution, "
                << "Time Step/Iteration Level = " << Number_of_Time_Steps
                << ", Time = " << Time
                << "\"" << "\n"
	        << "VARIABLES = \"x\" \\ \n"
                << "\"y\" \\ \n"
                << "\"ni\" \\ \n"
                << "\"ui\" \\ \n"
                << "\"vi\" \\ \n"
                << "\"pi\" \\ \n"
                << "\"Ti\" \n"
                << "\"Mi\" \n"
                << "\"rhon\" \\ \n"
                << "\"un\" \\ \n"
                << "\"vn\" \\ \n"
                << "\"pn\" \\ \n"
                << "\"Tn\" \n"
                << "\"Mn\" \n"
                << "\"Ex\" \\ \n"
                << "\"Ey\" \\ \n"
                << "\"V\" \\ \n"
                << "\"Jx\" \\ \n"
                << "\"Jy\" \\ \n"
                << "\"J\" \\ \n"
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
	   Wneut_node = SolnBlk.Wneutn(i, j);
	   E_node = SolnBlk.En(i, j);
	   V_node = SolnBlk.Vn(i, j);
           J_node = (W_node.q*W_node.n())*W_node.v;
           Out_File << " " << SolnBlk.Grid.Node[i][j].X;
           Out_File.setf(ios::scientific);
           Out_File << " " << W_node.n() 
                    << " " << W_node.v.x 
                    << " " << W_node.v.y 
                    << " " << W_node.p
                    << " " << W_node.T()
                    << " " << W_node.M();
           Out_File.unsetf(ios::scientific);
           Out_File << Wneut_node;
           Out_File.setf(ios::scientific);
           Out_File << " " << Wneut_node.T()
                    << " " << Wneut_node.M();
           Out_File.unsetf(ios::scientific);
           Out_File << E_node;
           Out_File.setf(ios::scientific);
           Out_File << " " << V_node;
           Out_File.unsetf(ios::scientific);
           Out_File << J_node;
           Out_File.setf(ios::scientific);
           Out_File << " " << abs(J_node) << "\n";
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
void Output_Cells_Tecplot(Ion5Moment2D_Quad_Block &SolnBlk,
                          Ion5Moment2D_Input_Parameters &Input_Parameters,
                          const int Number_of_Time_Steps,
                          const double &Time,
                          const int Block_Number,
                          const int Output_Title,
	                  ostream &Out_File) {

    int i, j;

    Out_File << setprecision(14);
    if (Output_Title) {
       Out_File << "TITLE = \"" << CFFC_Name() << ": 2D 5-Moment Ion Transport Solution, "
                << "Time Step/Iteration Level = " << Number_of_Time_Steps
                << ", Time = " << Time
                << "\"" << "\n"
   	        << "VARIABLES = \"x\" \\ \n"
                << "\"y\" \\ \n"
                << "\"rhoi\" \\ \n"
                << "\"ui\" \\ \n"
                << "\"vi\" \\ \n"
                << "\"pi\" \\ \n"
                << "\"Ti\" \n"
                << "\"rhon\" \\ \n"
                << "\"un\" \\ \n"
                << "\"vn\" \\ \n"
                << "\"pn\" \\ \n"
                << "\"Tn\" \n"
                << "\"Ex\" \\ \n"
                << "\"Ey\" \\ \n"
                << "\"V\" \\ \n"
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
           Out_File << " " << SolnBlk.W[i][j].T();
           Out_File.unsetf(ios::scientific);
           Out_File << SolnBlk.Wneut[i][j];
           Out_File.setf(ios::scientific);
           Out_File << " " << SolnBlk.Wneut[i][j].T();
           Out_File.unsetf(ios::scientific);
           Out_File << SolnBlk.E[i][j];
           Out_File.setf(ios::scientific);
           Out_File << " " << SolnBlk.V[i][j] << "\n";
           Out_File.unsetf(ios::scientific);
       } /* endfor */
    } /* endfor */
    Out_File << setprecision(6);
    
}

/********************************************************
 * Routine: ICs                                         *
 *                                                      *
 * Assigns initial conditions and data to the ion       *
 * solution variables of the specified quadrilateral    *
 * solution block.                                      *
 *                                                      *
 ********************************************************/
void ICs(Ion5Moment2D_Quad_Block &SolnBlk,
	 const int i_ICtype,
         Ion5Moment2D_pState *Wo) {

    int i, j, k;
    Ion5Moment2D_pState Wl, Wr;

    /* Assign initial data for the ion flow IVP of interest. */

    switch(i_ICtype) {
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
        Wl = Ion5Moment2D_W_STDATM;
        Wr = Ion5Moment2D_pState(DENSITY_STDATM/EIGHT,
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
        Wl = Ion5Moment2D_W_STDATM;
        Wr = Ion5Moment2D_pState(DENSITY_STDATM/EIGHT,
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
        Wl = Ion5Moment2D_pState(4.696, ZERO, ZERO, 404.4e03);
        Wr = Ion5Moment2D_pState(1.408, ZERO, ZERO, 101.1e03);
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
        Wl = Ion5Moment2D_pState(4.696, ZERO, ZERO, 404.4e03);
        Wr = Ion5Moment2D_pState(1.408, ZERO, ZERO, 101.1e03);
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
        Wl = Ion5Moment2D_pState(ONE, -TWO, ZERO, FOUR/TEN);
        Wr = Ion5Moment2D_pState(ONE, TWO, ZERO, FOUR/TEN);
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
        Wl = Ion5Moment2D_pState(ONE, ZERO, -TWO, FOUR/TEN);
        Wr = Ion5Moment2D_pState(ONE, ZERO, TWO, FOUR/TEN);
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
        Wl = Ion5Moment2D_pState(2.281, 164.83, ZERO, 201.17e03);
        Wr = Ion5Moment2D_pState(1.408, ZERO, ZERO, 101.1e03);
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
        Wl = Ion5Moment2D_pState(2.281, ZERO, 164.83, 201.17e03);
        Wr = Ion5Moment2D_pState(1.408, ZERO, ZERO, 101.1e03);
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
      case IC_CONTACT_SURFACE_XDIR :
        // Set initial data for moving contact surface propagating in x-direction.
        Wl = Ion5Moment2D_pState(1.045, 200.00, ZERO, 300.00e03);
        Wr = Ion5Moment2D_pState(3.483, 200.00, ZERO, 300.00e03);
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
        Wl = Ion5Moment2D_pState(1.045, ZERO, 200.00, 300.00e03);
        Wr = Ion5Moment2D_pState(3.483, ZERO, 200.00, 300.00e03);
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
        Wl = Ion5Moment2D_pState(1.598, -383.64, ZERO, 91.88e03);
        Wr = Ion5Moment2D_pState(2.787, -216.97, ZERO, 200.0e03);
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
        Wl = Ion5Moment2D_pState(1.598, ZERO, -383.64, 91.88e03);
        Wr = Ion5Moment2D_pState(2.787, ZERO, -216.97, 200.0e03);
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
        Wl = Ion5Moment2D_W_STDATM;
        Wr = Ion5Moment2D_pState(DENSITY_STDATM*FOUR,
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
        Wr = Wo[0];
        Wl = Ion5Moment2D_pState(HUNDRED*Wo[0].d,
          		         ZERO, ZERO,
         		         HUNDRED*Wo[0].p);
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
        Wl = Wo[0];
        Wr = Ion5Moment2D_pState(Wo[0].d/(HUNDRED*THOUSAND),
          		         ZERO, ZERO,
         		         Wo[0].p/(HUNDRED*THOUSAND));
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

    /* Set the solution residuals, gradients, limiters, and 
       other values to zero. */

    for (j  = SolnBlk.JCl-SolnBlk.Nghost ; j <= SolnBlk.JCu+SolnBlk.Nghost ; ++j ) {
       for ( i = SolnBlk.ICl-SolnBlk.Nghost ; i <= SolnBlk.ICu+SolnBlk.Nghost ; ++i ) {
          for ( k = 0 ; k <= NUMBER_OF_RESIDUAL_VECTORS_ION5MOMENT2D-1 ; ++k ) {
	     SolnBlk.dUdt[i][j][k] = Ion5Moment2D_U_VACUUM;
          } /* endfor */
	  SolnBlk.dWdx[i][j] = Ion5Moment2D_W_VACUUM;
	  SolnBlk.dWdy[i][j] = Ion5Moment2D_W_VACUUM;
	  SolnBlk.phi[i][j] = Ion5Moment2D_W_VACUUM;
	  SolnBlk.Uo[i][j] = Ion5Moment2D_U_VACUUM;
	  SolnBlk.dt[i][j] = ZERO;
       } /* endfor */
    } /* endfor */

    /* Set default values for the ion boundary condition reference states. */

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

/********************************************************
 * Routine: ICs_Neutral_Gas                             *
 *                                                      *
 * Assigns initial conditions and data to the neutral   *
 * gas solution variables of the specified quadrilateral*
 * solution block.                                      *
 *                                                      *
 ********************************************************/
void ICs_Neutral_Gas(Ion5Moment2D_Quad_Block &SolnBlk,
	             const int i_ICtype,
                     Euler2D_pState &Wno) {

    int i, j;

    /* Assign initial data for the neutral gas flow. */

    for (j  = SolnBlk.JCl-SolnBlk.Nghost ; j <= SolnBlk.JCu+SolnBlk.Nghost ; ++j ) {
       for ( i = SolnBlk.ICl-SolnBlk.Nghost ; i <= SolnBlk.ICu+SolnBlk.Nghost ; ++i ) {
          SolnBlk.Wneut[i][j] = Wno;
       } /* endfor */
    } /* endfor */

}

/********************************************************
 * Routine: ICs_Electric_Field                          *
 *                                                      *
 * Assigns initial conditions and data to the electric  *
 * field solution variables of the specified            *
 * quadrilateral solution block.                        *
 *                                                      *
 ********************************************************/
void ICs_Electric_Field(Ion5Moment2D_Quad_Block &SolnBlk,
                        const int i_Electric_Field,
 	                const double &Electric_Field_Strength,
 	                const double &Electric_Field_Angle) {

    int i, j;
    double q, voltage_rf_peak_to_peak, frequency_rf, 
           v_pseudo_potential_rfrods, r_pseudo_potential,
           x0_pseudo_potential, x1_pseudo_potential, 
           x2_pseudo_potential, x3_pseudo_potential, 
           xx;

    /* Set various constant for electric field calculation. */

    voltage_rf_peak_to_peak = 120.00;
    frequency_rf = 6.28e06;
    r_pseudo_potential = 3.0e-03;
    q = TWO*SolnBlk.W[SolnBlk.JCl][SolnBlk.ICl].q*voltage_rf_peak_to_peak/
        (SolnBlk.W[SolnBlk.JCl][SolnBlk.ICl].m*sqr(frequency_rf*r_pseudo_potential));
    q = 0.3;
    x0_pseudo_potential = 0.00;
    x1_pseudo_potential = 6.5e-03;
    x2_pseudo_potential = 60.0e-03;
    x3_pseudo_potential = 66.5e-03;

    /* Assign initial data for the electric field. */

    switch(i_Electric_Field) {
      case IC_ELECTRIC_FIELD_QUADRUPOLE :
        // Psuedo-potential field associated with quadrapole time-varying fields.
        v_pseudo_potential_rfrods = sqr(q*frequency_rf)*SolnBlk.W[SolnBlk.JCl][SolnBlk.ICl].m/
                                    (16.00*SolnBlk.W[SolnBlk.JCl][SolnBlk.ICl].q);
        for (j  = SolnBlk.JCl-SolnBlk.Nghost ; j <= SolnBlk.JCu+SolnBlk.Nghost ; ++j ) {
          for ( i = SolnBlk.ICl-SolnBlk.Nghost ; i <= SolnBlk.ICu+SolnBlk.Nghost ; ++i ) {
             if (SolnBlk.Grid.Cell[i][j].Xc.x <= x0_pseudo_potential) {
                SolnBlk.V[i][j] = ZERO;
                SolnBlk.E[i][j] = Vector2D_ZERO;
             } else if (SolnBlk.Grid.Cell[i][j].Xc.x <= x1_pseudo_potential) {
                xx = (SolnBlk.Grid.Cell[i][j].Xc.x-x0_pseudo_potential)/
                     (x1_pseudo_potential-x0_pseudo_potential);
                if (SolnBlk.Grid.Cell[i][j].Xc.y > r_pseudo_potential) {
                   SolnBlk.V[i][j] = v_pseudo_potential_rfrods*xx*sqr(r_pseudo_potential);
                   SolnBlk.E[i][j].x = -v_pseudo_potential_rfrods*sqr(r_pseudo_potential)/
                                       (x1_pseudo_potential-x0_pseudo_potential);
                   SolnBlk.E[i][j].y = ZERO;
                } else {
                   SolnBlk.V[i][j] = v_pseudo_potential_rfrods*xx*sqr(SolnBlk.Grid.Cell[i][j].Xc.y);
                   SolnBlk.E[i][j].x = -v_pseudo_potential_rfrods*sqr(SolnBlk.Grid.Cell[i][j].Xc.y)/
                                       (x1_pseudo_potential-x0_pseudo_potential);
                   SolnBlk.E[i][j].y = -TWO*v_pseudo_potential_rfrods*xx*SolnBlk.Grid.Cell[i][j].Xc.y;
                } /* endif */
             } else if (SolnBlk.Grid.Cell[i][j].Xc.x <= x2_pseudo_potential) {
	        if (SolnBlk.Grid.Cell[i][j].Xc.y > r_pseudo_potential) {
                   SolnBlk.V[i][j] = v_pseudo_potential_rfrods*sqr(r_pseudo_potential);
                   SolnBlk.E[i][j] = Vector2D_ZERO;
                } else {
                   SolnBlk.V[i][j] = v_pseudo_potential_rfrods*sqr(SolnBlk.Grid.Cell[i][j].Xc.y);
                   SolnBlk.E[i][j].x = ZERO;
                   SolnBlk.E[i][j].y = -TWO*v_pseudo_potential_rfrods*SolnBlk.Grid.Cell[i][j].Xc.y;
                } /* endif */
             } else if (SolnBlk.Grid.Cell[i][j].Xc.x <= x3_pseudo_potential) {
	        xx = (SolnBlk.Grid.Cell[i][j].Xc.x-x2_pseudo_potential)/
                     (x3_pseudo_potential-x2_pseudo_potential);
	        if (SolnBlk.Grid.Cell[i][j].Xc.y > r_pseudo_potential) {
                   SolnBlk.V[i][j] = v_pseudo_potential_rfrods*(ONE-xx)*sqr(r_pseudo_potential);
                   SolnBlk.E[i][j].x = v_pseudo_potential_rfrods*sqr(r_pseudo_potential)/
                                       (x3_pseudo_potential-x2_pseudo_potential);
                   SolnBlk.E[i][j].y = ZERO;
                } else {
                   SolnBlk.V[i][j] = v_pseudo_potential_rfrods*(ONE-xx)*sqr(SolnBlk.Grid.Cell[i][j].Xc.y);
                   SolnBlk.E[i][j].x = v_pseudo_potential_rfrods*sqr(SolnBlk.Grid.Cell[i][j].Xc.y)/
                                       (x3_pseudo_potential-x2_pseudo_potential);
                   SolnBlk.E[i][j].y = -TWO*v_pseudo_potential_rfrods*(ONE-xx)*SolnBlk.Grid.Cell[i][j].Xc.y;
                } /* endif */
             } else {
                SolnBlk.V[i][j] = ZERO;
                SolnBlk.E[i][j] = Vector2D_ZERO;
             } /* endif */
          } /* endfor */
        } /* endfor */
        break;

      case IC_ELECTRIC_FIELD_OCTAPOLE :
        // Psuedo-potential field associated with octapole time-varying fields.
        v_pseudo_potential_rfrods = sqr(q*frequency_rf)*SolnBlk.W[SolnBlk.JCl][SolnBlk.ICl].m/
                                    (FOUR*SolnBlk.W[SolnBlk.JCl][SolnBlk.ICl].q);
        for (j  = SolnBlk.JCl-SolnBlk.Nghost ; j <= SolnBlk.JCu+SolnBlk.Nghost ; ++j ) {
          for ( i = SolnBlk.ICl-SolnBlk.Nghost ; i <= SolnBlk.ICu+SolnBlk.Nghost ; ++i ) {
             if (SolnBlk.Grid.Cell[i][j].Xc.x <= x0_pseudo_potential) {
                SolnBlk.V[i][j] = ZERO;
                SolnBlk.E[i][j] = Vector2D_ZERO;
             } else if (SolnBlk.Grid.Cell[i][j].Xc.x <= x1_pseudo_potential) {
	        xx = (SolnBlk.Grid.Cell[i][j].Xc.x-x0_pseudo_potential)/
                     (x1_pseudo_potential-x0_pseudo_potential);
                if (SolnBlk.Grid.Cell[i][j].Xc.y > r_pseudo_potential) {
                   SolnBlk.V[i][j] = v_pseudo_potential_rfrods*xx*pow(r_pseudo_potential, SIX);
                   SolnBlk.E[i][j].x = -v_pseudo_potential_rfrods*pow(r_pseudo_potential, SIX)/
                                       (x1_pseudo_potential-x0_pseudo_potential);
                   SolnBlk.E[i][j].y = ZERO;
                } else {
                   SolnBlk.V[i][j] = v_pseudo_potential_rfrods*xx*pow(SolnBlk.Grid.Cell[i][j].Xc.y, SIX);
                   SolnBlk.E[i][j].x = -v_pseudo_potential_rfrods*pow(SolnBlk.Grid.Cell[i][j].Xc.y, SIX)/
                                       (x1_pseudo_potential-x0_pseudo_potential);
                   SolnBlk.E[i][j].y = -SIX*v_pseudo_potential_rfrods*xx*pow(SolnBlk.Grid.Cell[i][j].Xc.y, FIVE);
                } /* endif */
             } else if (SolnBlk.Grid.Cell[i][j].Xc.x <= x2_pseudo_potential) {
	        if (SolnBlk.Grid.Cell[i][j].Xc.y > r_pseudo_potential) {
                   SolnBlk.V[i][j] = v_pseudo_potential_rfrods*pow(r_pseudo_potential, SIX);
                   SolnBlk.E[i][j] = Vector2D_ZERO;
                } else {
                   SolnBlk.V[i][j] = v_pseudo_potential_rfrods*pow(SolnBlk.Grid.Cell[i][j].Xc.y, SIX);
                   SolnBlk.E[i][j].x = ZERO;
                   SolnBlk.E[i][j].y = -SIX*v_pseudo_potential_rfrods*pow(SolnBlk.Grid.Cell[i][j].Xc.y, FIVE);
                } /* endif */
             } else if (SolnBlk.Grid.Cell[i][j].Xc.x <= x3_pseudo_potential) {
	        xx = (SolnBlk.Grid.Cell[i][j].Xc.x-x2_pseudo_potential)/
                     (x3_pseudo_potential-x2_pseudo_potential);
	        if (SolnBlk.Grid.Cell[i][j].Xc.y > r_pseudo_potential) {
                   SolnBlk.V[i][j] = v_pseudo_potential_rfrods*(ONE-xx)*pow(r_pseudo_potential, SIX);
                   SolnBlk.E[i][j].x = v_pseudo_potential_rfrods*pow(r_pseudo_potential, SIX)/
                                       (x3_pseudo_potential-x2_pseudo_potential);
                   SolnBlk.E[i][j].y = ZERO;
                } else {
                   SolnBlk.V[i][j] = v_pseudo_potential_rfrods*(ONE-xx)*pow(SolnBlk.Grid.Cell[i][j].Xc.y, SIX);
                   SolnBlk.E[i][j].x = v_pseudo_potential_rfrods*pow(SolnBlk.Grid.Cell[i][j].Xc.y, SIX)/
                                       (x3_pseudo_potential-x2_pseudo_potential);
                   SolnBlk.E[i][j].y = -SIX*v_pseudo_potential_rfrods*(ONE-xx)*pow(SolnBlk.Grid.Cell[i][j].Xc.y, SIX);
                } /* endif */
             } else {
                SolnBlk.V[i][j] = ZERO;
                SolnBlk.E[i][j] = Vector2D_ZERO;
             } /* endif */
          } /* endfor */
        } /* endfor */
        break;

      case IC_ELECTRIC_FIELD_UNIFORM :
      default:
        // Uniform electric field.
        for (j  = SolnBlk.JCl-SolnBlk.Nghost ; j <= SolnBlk.JCu+SolnBlk.Nghost ; ++j ) {
          for ( i = SolnBlk.ICl-SolnBlk.Nghost ; i <= SolnBlk.ICu+SolnBlk.Nghost ; ++i ) {
	     //if (SolnBlk.Grid.Cell[i][j].Xc.x <= -0.001) {
             //   SolnBlk.E[i][j] = Vector2D_ZERO;
             //} else {
                SolnBlk.E[i][j].x = Electric_Field_Strength*cos(TWO*PI*Electric_Field_Angle/360.00);
                SolnBlk.E[i][j].y = Electric_Field_Strength*sin(TWO*PI*Electric_Field_Angle/360.00);
             //} /* endif */
             SolnBlk.V[i][j] = ZERO;
          } /* endfor */
        } /* endfor */
        break;

    } /* endswitch */

}

/********************************************************
 * Routine: BCs                                         *
 *                                                      *
 * Apply ion boundary conditions at boundaries of the   *
 * specified quadrilateral solution block.              *
 *                                                      *
 ********************************************************/
void BCs(Ion5Moment2D_Quad_Block &SolnBlk) {

    int i, j;
    Vector2D dX;
    Ion5Moment2D_pState dW, W, W_wall_src;

    W_wall_src = Ion5Moment2D_pState(1.0e08*W_wall_src.m, Vector2D_ZERO, 
                                     1.0e08*W_wall_src.m*W_wall_src.R*TEMPERATURE_STDATM);

    for ( j = SolnBlk.JCl-SolnBlk.Nghost ; j <= SolnBlk.JCu+SolnBlk.Nghost ; ++j ) {
      if ( (j >= SolnBlk.JCl && j <= SolnBlk.JCu) ||
           (j < SolnBlk.JCl && 
            (SolnBlk.Grid.BCtypeS[SolnBlk.ICl-1] == BC_NONE ||
             SolnBlk.Grid.BCtypeS[SolnBlk.ICl-1] == BC_PERIODIC ||
             SolnBlk.Grid.BCtypeS[SolnBlk.ICl-1] == BC_CONSTANT_EXTRAPOLATION ||
             SolnBlk.Grid.BCtypeS[SolnBlk.ICl-1] == BC_FIXED_TEMP_WALL ||
             SolnBlk.Grid.BCtypeS[SolnBlk.ICl-1] == BC_ADIABATIC_WALL ||
             SolnBlk.Grid.BCtypeS[SolnBlk.ICl-1] == BC_LINEAR_EXTRAPOLATION ||
             SolnBlk.Grid.BCtypeS[SolnBlk.ICl-1] == BC_CHARACTERISTIC) ) ||
           (j > SolnBlk.JCu && 
            (SolnBlk.Grid.BCtypeN[SolnBlk.ICl-1] == BC_NONE ||
             SolnBlk.Grid.BCtypeN[SolnBlk.ICl-1] == BC_PERIODIC ||
             SolnBlk.Grid.BCtypeN[SolnBlk.ICl-1] == BC_CONSTANT_EXTRAPOLATION ||
             SolnBlk.Grid.BCtypeN[SolnBlk.ICl-1] == BC_FIXED_TEMP_WALL ||
             SolnBlk.Grid.BCtypeN[SolnBlk.ICl-1] == BC_ADIABATIC_WALL ||
             SolnBlk.Grid.BCtypeN[SolnBlk.ICl-1] == BC_LINEAR_EXTRAPOLATION ||
             SolnBlk.Grid.BCtypeN[SolnBlk.ICl-1] == BC_CHARACTERISTIC) ) ) {
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
          case BC_FIXED_TEMP_WALL :
          case BC_ADIABATIC_WALL :
            if (SolnBlk.W[SolnBlk.ICl][j].v*SolnBlk.Grid.nfaceW(SolnBlk.ICl, j) >= ZERO) {
               SolnBlk.W[SolnBlk.ICl-1][j] = SolnBlk.W[SolnBlk.ICl][j];
               SolnBlk.U[SolnBlk.ICl-1][j] = SolnBlk.U[SolnBlk.ICl][j];
               SolnBlk.W[SolnBlk.ICl-2][j] = SolnBlk.W[SolnBlk.ICl][j];
               SolnBlk.U[SolnBlk.ICl-2][j] = SolnBlk.U[SolnBlk.ICl][j];
            } else {
               SolnBlk.W[SolnBlk.ICl-1][j] = W_wall_src;
               //SolnBlk.W[SolnBlk.ICl-1][j] = SolnBlk.WoW[j];
               SolnBlk.U[SolnBlk.ICl-1][j] = U(SolnBlk.W[SolnBlk.ICl-1][j]);
               SolnBlk.W[SolnBlk.ICl-2][j] = W_wall_src;
               //SolnBlk.W[SolnBlk.ICl-2][j] = SolnBlk.WoW[j];
               SolnBlk.U[SolnBlk.ICl-2][j] = U(SolnBlk.W[SolnBlk.ICl-2][j]);
            } /* endif */
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
             SolnBlk.Grid.BCtypeS[SolnBlk.ICu+1] == BC_FIXED_TEMP_WALL ||
             SolnBlk.Grid.BCtypeS[SolnBlk.ICu+1] == BC_ADIABATIC_WALL ||
             SolnBlk.Grid.BCtypeS[SolnBlk.ICu+1] == BC_LINEAR_EXTRAPOLATION ||
             SolnBlk.Grid.BCtypeS[SolnBlk.ICu+1] == BC_CHARACTERISTIC) ) ||
           (j > SolnBlk.JCu && 
            (SolnBlk.Grid.BCtypeN[SolnBlk.ICu+1] == BC_NONE ||
             SolnBlk.Grid.BCtypeN[SolnBlk.ICu+1] == BC_PERIODIC ||
             SolnBlk.Grid.BCtypeN[SolnBlk.ICu+1] == BC_CONSTANT_EXTRAPOLATION ||
             SolnBlk.Grid.BCtypeN[SolnBlk.ICu+1] == BC_FIXED_TEMP_WALL ||
             SolnBlk.Grid.BCtypeN[SolnBlk.ICu+1] == BC_ADIABATIC_WALL ||
             SolnBlk.Grid.BCtypeN[SolnBlk.ICu+1] == BC_LINEAR_EXTRAPOLATION ||
             SolnBlk.Grid.BCtypeN[SolnBlk.ICu+1] == BC_CHARACTERISTIC) ) ) {
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
          case BC_FIXED_TEMP_WALL :
          case BC_ADIABATIC_WALL :
            if (SolnBlk.W[i][SolnBlk.JCl].v*SolnBlk.Grid.nfaceS(i, SolnBlk.JCl) >= ZERO) {
               SolnBlk.W[i][SolnBlk.JCl-1] = SolnBlk.W[i][SolnBlk.JCl];
               SolnBlk.U[i][SolnBlk.JCl-1] = SolnBlk.U[i][SolnBlk.JCl];
               SolnBlk.W[i][SolnBlk.JCl-2] = SolnBlk.W[i][SolnBlk.JCl];
               SolnBlk.U[i][SolnBlk.JCl-2] = SolnBlk.U[i][SolnBlk.JCl];
            } else {
               SolnBlk.W[i][SolnBlk.JCl-1] = W_wall_src;
               //SolnBlk.W[i][SolnBlk.JCl-1] = SolnBlk.WoS[i];
               SolnBlk.U[i][SolnBlk.JCl-1] = U(SolnBlk.W[i][SolnBlk.JCl-1]);
               SolnBlk.W[i][SolnBlk.JCl-2] = W_wall_src;
               //SolnBlk.W[i][SolnBlk.JCl-2] = SolnBlk.WoS[i];
               SolnBlk.U[i][SolnBlk.JCl-2] = U(SolnBlk.W[i][SolnBlk.JCl-2]);
            } /* endif */
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
          case BC_FIXED_TEMP_WALL :
          case BC_ADIABATIC_WALL :
            if (SolnBlk.W[i][SolnBlk.JCu].v*SolnBlk.Grid.nfaceN(i, SolnBlk.JCu) >= ZERO) {
               SolnBlk.W[i][SolnBlk.JCu+1] = SolnBlk.W[i][SolnBlk.JCu];
               SolnBlk.U[i][SolnBlk.JCu+1] = SolnBlk.U[i][SolnBlk.JCu];
               SolnBlk.W[i][SolnBlk.JCu+2] = SolnBlk.W[i][SolnBlk.JCu];
               SolnBlk.U[i][SolnBlk.JCu+2] = SolnBlk.U[i][SolnBlk.JCu];
            } else {
               SolnBlk.W[i][SolnBlk.JCu+1] = W_wall_src;
               //SolnBlk.W[i][SolnBlk.JCu+1] = SolnBlk.WoN[i];
               SolnBlk.U[i][SolnBlk.JCu+1] = U(SolnBlk.W[i][SolnBlk.JCu+1]);
               SolnBlk.W[i][SolnBlk.JCu+2] = W_wall_src;
               //SolnBlk.W[i][SolnBlk.JCu+2] = SolnBlk.WoN[i];
               SolnBlk.U[i][SolnBlk.JCu+2] = U(SolnBlk.W[i][SolnBlk.JCu+2]);
            } /* endif */
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
          default:
            SolnBlk.W[i][SolnBlk.JCu+1] = SolnBlk.W[i][SolnBlk.JCu];
            SolnBlk.U[i][SolnBlk.JCu+1] = SolnBlk.U[i][SolnBlk.JCu];
            SolnBlk.W[i][SolnBlk.JCu+2] = SolnBlk.W[i][SolnBlk.JCu];
            SolnBlk.U[i][SolnBlk.JCu+2] = SolnBlk.U[i][SolnBlk.JCu];
            break;
        } /* endswitch */
    } /* endfor */

}

/********************************************************
 * Routine: BCs_Neutral_Gas                             *
 *                                                      *
 * Apply boundary conditions for the neutral gas flow   *
 * at boundaries of the specified quadrilateral         *
 * solution block.                                      *
 *                                                      *
 ********************************************************/
void BCs_Neutral_Gas(Ion5Moment2D_Quad_Block &SolnBlk) {

    int i, j;

    for ( j = SolnBlk.JCl-SolnBlk.Nghost ; j <= SolnBlk.JCu+SolnBlk.Nghost ; ++j ) {
      if ( (j >= SolnBlk.JCl && j <= SolnBlk.JCu) ||
           (j < SolnBlk.JCl && 
            (SolnBlk.Grid.BCtypeS[SolnBlk.ICl-1] == BC_NONE ||
             SolnBlk.Grid.BCtypeS[SolnBlk.ICl-1] == BC_PERIODIC ||
             SolnBlk.Grid.BCtypeS[SolnBlk.ICl-1] == BC_CONSTANT_EXTRAPOLATION ||
             SolnBlk.Grid.BCtypeS[SolnBlk.ICl-1] == BC_REFLECTION ||
             SolnBlk.Grid.BCtypeS[SolnBlk.ICl-1] == BC_FIXED_TEMP_WALL ||
             SolnBlk.Grid.BCtypeS[SolnBlk.ICl-1] == BC_ADIABATIC_WALL ||
             SolnBlk.Grid.BCtypeS[SolnBlk.ICl-1] == BC_LINEAR_EXTRAPOLATION ||
             SolnBlk.Grid.BCtypeS[SolnBlk.ICl-1] == BC_CHARACTERISTIC) ) ||
           (j > SolnBlk.JCu && 
            (SolnBlk.Grid.BCtypeN[SolnBlk.ICl-1] == BC_NONE ||
             SolnBlk.Grid.BCtypeN[SolnBlk.ICl-1] == BC_PERIODIC ||
             SolnBlk.Grid.BCtypeN[SolnBlk.ICl-1] == BC_CONSTANT_EXTRAPOLATION ||
             SolnBlk.Grid.BCtypeN[SolnBlk.ICl-1] == BC_REFLECTION ||
             SolnBlk.Grid.BCtypeN[SolnBlk.ICl-1] == BC_FIXED_TEMP_WALL ||
             SolnBlk.Grid.BCtypeN[SolnBlk.ICl-1] == BC_ADIABATIC_WALL ||
             SolnBlk.Grid.BCtypeN[SolnBlk.ICl-1] == BC_LINEAR_EXTRAPOLATION ||
             SolnBlk.Grid.BCtypeN[SolnBlk.ICl-1] == BC_CHARACTERISTIC) ) ) {
        switch(SolnBlk.Grid.BCtypeW[j]) {
          case BC_NONE :
          case BC_FIXED :
          case BC_CONSTANT_EXTRAPOLATION :
          case BC_LINEAR_EXTRAPOLATION :
          case BC_CHARACTERISTIC :
            SolnBlk.Wneut[SolnBlk.ICl-1][j] = SolnBlk.Wneut[SolnBlk.ICl][j];
            SolnBlk.Wneut[SolnBlk.ICl-2][j] = SolnBlk.Wneut[SolnBlk.ICl][j];
            break;
          case BC_REFLECTION :
            SolnBlk.Wneut[SolnBlk.ICl-1][j] = Reflect(SolnBlk.Wneut[SolnBlk.ICl][j],
				 		      SolnBlk.Grid.nfaceW(SolnBlk.ICl,j));
            SolnBlk.Wneut[SolnBlk.ICl-2][j] = Reflect(SolnBlk.Wneut[SolnBlk.ICl+1][j],
				 		      SolnBlk.Grid.nfaceW(SolnBlk.ICl,j));
            break;
          case BC_FIXED_TEMP_WALL :
          case BC_ADIABATIC_WALL :
            SolnBlk.Wneut[SolnBlk.ICl-1][j] = NoSlip(SolnBlk.Wneut[SolnBlk.ICl][j],
				 		     SolnBlk.Grid.nfaceW(SolnBlk.ICl,j));
            SolnBlk.Wneut[SolnBlk.ICl-2][j] = NoSlip(SolnBlk.Wneut[SolnBlk.ICl+1][j],
				 		     SolnBlk.Grid.nfaceW(SolnBlk.ICl,j));
            break;
          case BC_PERIODIC :
            SolnBlk.Wneut[SolnBlk.ICl-1][j] = SolnBlk.Wneut[SolnBlk.ICu-1][j];
            SolnBlk.Wneut[SolnBlk.ICl-2][j] = SolnBlk.Wneut[SolnBlk.ICu-2][j];
            break;
          default:
            SolnBlk.Wneut[SolnBlk.ICl-1][j] = SolnBlk.Wneut[SolnBlk.ICl][j];
            SolnBlk.Wneut[SolnBlk.ICl-2][j] = SolnBlk.Wneut[SolnBlk.ICl][j];
            break;
        } /* endswitch */
      } /* endif */

      if ( (j >= SolnBlk.JCl && j <= SolnBlk.JCu) ||
           (j < SolnBlk.JCl && 
            (SolnBlk.Grid.BCtypeS[SolnBlk.ICu+1] == BC_NONE ||
             SolnBlk.Grid.BCtypeS[SolnBlk.ICu+1] == BC_PERIODIC ||
             SolnBlk.Grid.BCtypeS[SolnBlk.ICu+1] == BC_CONSTANT_EXTRAPOLATION ||
             SolnBlk.Grid.BCtypeS[SolnBlk.ICu+1] == BC_REFLECTION ||
             SolnBlk.Grid.BCtypeS[SolnBlk.ICu+1] == BC_FIXED_TEMP_WALL ||
             SolnBlk.Grid.BCtypeS[SolnBlk.ICu+1] == BC_ADIABATIC_WALL ||
             SolnBlk.Grid.BCtypeS[SolnBlk.ICu+1] == BC_LINEAR_EXTRAPOLATION ||
             SolnBlk.Grid.BCtypeS[SolnBlk.ICu+1] == BC_CHARACTERISTIC) ) ||
           (j > SolnBlk.JCu && 
            (SolnBlk.Grid.BCtypeN[SolnBlk.ICu+1] == BC_NONE ||
             SolnBlk.Grid.BCtypeN[SolnBlk.ICu+1] == BC_PERIODIC ||
             SolnBlk.Grid.BCtypeN[SolnBlk.ICu+1] == BC_CONSTANT_EXTRAPOLATION ||
             SolnBlk.Grid.BCtypeN[SolnBlk.ICu+1] == BC_REFLECTION ||
             SolnBlk.Grid.BCtypeN[SolnBlk.ICu+1] == BC_FIXED_TEMP_WALL ||
             SolnBlk.Grid.BCtypeN[SolnBlk.ICu+1] == BC_ADIABATIC_WALL ||
             SolnBlk.Grid.BCtypeN[SolnBlk.ICu+1] == BC_LINEAR_EXTRAPOLATION ||
             SolnBlk.Grid.BCtypeN[SolnBlk.ICu+1] == BC_CHARACTERISTIC) ) ) {
        switch(SolnBlk.Grid.BCtypeE[j]) {
          case BC_NONE :
          case BC_FIXED :
          case BC_CONSTANT_EXTRAPOLATION :
          case BC_LINEAR_EXTRAPOLATION :
          case BC_CHARACTERISTIC :
            SolnBlk.Wneut[SolnBlk.ICu+1][j] = SolnBlk.Wneut[SolnBlk.ICu][j];
            SolnBlk.Wneut[SolnBlk.ICu+2][j] = SolnBlk.Wneut[SolnBlk.ICu][j];
            break;
          case BC_REFLECTION :
            SolnBlk.Wneut[SolnBlk.ICu+1][j] = Reflect(SolnBlk.Wneut[SolnBlk.ICu][j],
                                                      SolnBlk.Grid.nfaceE(SolnBlk.ICu,j));
            SolnBlk.Wneut[SolnBlk.ICu+2][j] = Reflect(SolnBlk.Wneut[SolnBlk.ICu-1][j],
                                                      SolnBlk.Grid.nfaceE(SolnBlk.ICu,j));
            break;
          case BC_FIXED_TEMP_WALL :
          case BC_ADIABATIC_WALL :
            SolnBlk.Wneut[SolnBlk.ICu+1][j] = NoSlip(SolnBlk.Wneut[SolnBlk.ICu][j],
                                                     SolnBlk.Grid.nfaceE(SolnBlk.ICu,j));
            SolnBlk.Wneut[SolnBlk.ICu+2][j] = NoSlip(SolnBlk.Wneut[SolnBlk.ICu-1][j],
                                                     SolnBlk.Grid.nfaceE(SolnBlk.ICu,j));
            break;
          case BC_PERIODIC :
            SolnBlk.Wneut[SolnBlk.ICu+1][j] = SolnBlk.Wneut[SolnBlk.ICl+1][j];
            SolnBlk.Wneut[SolnBlk.ICu+2][j] = SolnBlk.Wneut[SolnBlk.ICl+2][j];
            break;
          default:
            SolnBlk.Wneut[SolnBlk.ICu+1][j] = SolnBlk.Wneut[SolnBlk.ICu][j];
            SolnBlk.Wneut[SolnBlk.ICu+2][j] = SolnBlk.Wneut[SolnBlk.ICu][j];
            break;
        } /* endswitch */
      } /* endif */
    } /* endfor */

    for ( i = SolnBlk.ICl-SolnBlk.Nghost ; i <= SolnBlk.ICu+SolnBlk.Nghost ; ++i ) {
        switch(SolnBlk.Grid.BCtypeS[i]) {
          case BC_NONE :
          case BC_FIXED :
          case BC_CONSTANT_EXTRAPOLATION :
          case BC_LINEAR_EXTRAPOLATION :
          case BC_CHARACTERISTIC :
            SolnBlk.Wneut[i][SolnBlk.JCl-1] = SolnBlk.Wneut[i][SolnBlk.JCl];
            SolnBlk.Wneut[i][SolnBlk.JCl-2] = SolnBlk.Wneut[i][SolnBlk.JCl];
            break;
          case BC_REFLECTION :
            SolnBlk.Wneut[i][SolnBlk.JCl-1] = Reflect(SolnBlk.Wneut[i][SolnBlk.JCl],
                                                      SolnBlk.Grid.nfaceS(i, SolnBlk.JCl));
            SolnBlk.Wneut[i][SolnBlk.JCl-2] = Reflect(SolnBlk.Wneut[i][SolnBlk.JCl+1],
                                                      SolnBlk.Grid.nfaceS(i, SolnBlk.JCl));
            break;
          case BC_FIXED_TEMP_WALL :
          case BC_ADIABATIC_WALL :
            SolnBlk.Wneut[i][SolnBlk.JCl-1] = NoSlip(SolnBlk.Wneut[i][SolnBlk.JCl],
                                                     SolnBlk.Grid.nfaceS(i, SolnBlk.JCl));
            SolnBlk.Wneut[i][SolnBlk.JCl-2] = NoSlip(SolnBlk.Wneut[i][SolnBlk.JCl+1],
                                                     SolnBlk.Grid.nfaceS(i, SolnBlk.JCl));
            break;
          case BC_PERIODIC :
            SolnBlk.Wneut[i][SolnBlk.JCl-1] = SolnBlk.Wneut[i][SolnBlk.JCu-1];
            SolnBlk.Wneut[i][SolnBlk.JCl-2] = SolnBlk.Wneut[i][SolnBlk.JCu-2];
            break;
          default:
            SolnBlk.Wneut[i][SolnBlk.JCl-1] = SolnBlk.Wneut[i][SolnBlk.JCl];
            SolnBlk.Wneut[i][SolnBlk.JCl-2] = SolnBlk.Wneut[i][SolnBlk.JCl];
            break;
        } /* endswitch */

        switch(SolnBlk.Grid.BCtypeN[i]) {
          case BC_NONE :
          case BC_FIXED :
          case BC_CONSTANT_EXTRAPOLATION :
          case BC_LINEAR_EXTRAPOLATION :
          case BC_CHARACTERISTIC :
            SolnBlk.Wneut[i][SolnBlk.JCu+1] = SolnBlk.Wneut[i][SolnBlk.JCu];
            SolnBlk.Wneut[i][SolnBlk.JCu+2] = SolnBlk.Wneut[i][SolnBlk.JCu];
            break;
          case BC_REFLECTION :
            SolnBlk.Wneut[i][SolnBlk.JCu+1] = Reflect(SolnBlk.Wneut[i][SolnBlk.JCu],
                                                      SolnBlk.Grid.nfaceN(i, SolnBlk.JCu));
            SolnBlk.Wneut[i][SolnBlk.JCu+2] = Reflect(SolnBlk.Wneut[i][SolnBlk.JCu-1],
                                                      SolnBlk.Grid.nfaceN(i, SolnBlk.JCu));
            break;
          case BC_FIXED_TEMP_WALL :
          case BC_ADIABATIC_WALL :
            SolnBlk.Wneut[i][SolnBlk.JCu+1] = NoSlip(SolnBlk.Wneut[i][SolnBlk.JCu],
                                                     SolnBlk.Grid.nfaceN(i, SolnBlk.JCu));
            SolnBlk.Wneut[i][SolnBlk.JCu+2] = NoSlip(SolnBlk.Wneut[i][SolnBlk.JCu-1],
                                                     SolnBlk.Grid.nfaceN(i, SolnBlk.JCu));
            break;
          case BC_PERIODIC :
            SolnBlk.Wneut[i][SolnBlk.JCu+1] = SolnBlk.Wneut[i][SolnBlk.JCl+1];
            SolnBlk.Wneut[i][SolnBlk.JCu+2] = SolnBlk.Wneut[i][SolnBlk.JCl+2];
            break;
          default:
            SolnBlk.Wneut[i][SolnBlk.JCu+1] = SolnBlk.Wneut[i][SolnBlk.JCu];
            SolnBlk.Wneut[i][SolnBlk.JCu+2] = SolnBlk.Wneut[i][SolnBlk.JCu];
            break;
        } /* endswitch */
    } /* endfor */

}

/********************************************************
 * Routine: BCs_Electric_Field                          *
 *                                                      *
 * Apply boundary conditions for the electric field     *
 * at boundaries of the specified quadrilateral         *
 * solution block.                                      *
 *                                                      *
 ********************************************************/
void BCs_Electric_Field(Ion5Moment2D_Quad_Block &SolnBlk) {

    int i, j;

    for ( j = SolnBlk.JCl-SolnBlk.Nghost ; j <= SolnBlk.JCu+SolnBlk.Nghost ; ++j ) {
      if ( (j >= SolnBlk.JCl && j <= SolnBlk.JCu) ||
           (j < SolnBlk.JCl && 
            (SolnBlk.Grid.BCtypeS[SolnBlk.ICl-1] == BC_NONE ||
             SolnBlk.Grid.BCtypeS[SolnBlk.ICl-1] == BC_PERIODIC ||
             SolnBlk.Grid.BCtypeS[SolnBlk.ICl-1] == BC_CONSTANT_EXTRAPOLATION ||
             SolnBlk.Grid.BCtypeS[SolnBlk.ICl-1] == BC_REFLECTION ||
             SolnBlk.Grid.BCtypeS[SolnBlk.ICl-1] == BC_FIXED_TEMP_WALL ||
             SolnBlk.Grid.BCtypeS[SolnBlk.ICl-1] == BC_ADIABATIC_WALL ||
             SolnBlk.Grid.BCtypeS[SolnBlk.ICl-1] == BC_LINEAR_EXTRAPOLATION ||
             SolnBlk.Grid.BCtypeS[SolnBlk.ICl-1] == BC_CHARACTERISTIC) ) ||
           (j > SolnBlk.JCu && 
            (SolnBlk.Grid.BCtypeN[SolnBlk.ICl-1] == BC_NONE ||
             SolnBlk.Grid.BCtypeN[SolnBlk.ICl-1] == BC_PERIODIC ||
             SolnBlk.Grid.BCtypeN[SolnBlk.ICl-1] == BC_CONSTANT_EXTRAPOLATION ||
             SolnBlk.Grid.BCtypeN[SolnBlk.ICl-1] == BC_REFLECTION ||
             SolnBlk.Grid.BCtypeN[SolnBlk.ICl-1] == BC_FIXED_TEMP_WALL ||
             SolnBlk.Grid.BCtypeN[SolnBlk.ICl-1] == BC_ADIABATIC_WALL ||
             SolnBlk.Grid.BCtypeN[SolnBlk.ICl-1] == BC_LINEAR_EXTRAPOLATION ||
             SolnBlk.Grid.BCtypeN[SolnBlk.ICl-1] == BC_CHARACTERISTIC) ) ) {
        switch(SolnBlk.Grid.BCtypeW[j]) {
          case BC_NONE :
          case BC_FIXED :
          case BC_CONSTANT_EXTRAPOLATION :
          case BC_LINEAR_EXTRAPOLATION :
          case BC_CHARACTERISTIC :
          case BC_REFLECTION :
          case BC_FIXED_TEMP_WALL :
          case BC_ADIABATIC_WALL :
            SolnBlk.E[SolnBlk.ICl-1][j] = SolnBlk.E[SolnBlk.ICl][j];
            SolnBlk.E[SolnBlk.ICl-2][j] = SolnBlk.E[SolnBlk.ICl][j];
            SolnBlk.V[SolnBlk.ICl-1][j] = SolnBlk.V[SolnBlk.ICl][j];
            SolnBlk.V[SolnBlk.ICl-2][j] = SolnBlk.V[SolnBlk.ICl][j];
            break;
          case BC_PERIODIC :
            SolnBlk.E[SolnBlk.ICl-1][j] = SolnBlk.E[SolnBlk.ICu-1][j];
            SolnBlk.E[SolnBlk.ICl-2][j] = SolnBlk.E[SolnBlk.ICu-2][j];
            SolnBlk.V[SolnBlk.ICl-1][j] = SolnBlk.V[SolnBlk.ICu-1][j];
            SolnBlk.V[SolnBlk.ICl-2][j] = SolnBlk.V[SolnBlk.ICu-2][j];
            break;
          default:
            SolnBlk.E[SolnBlk.ICl-1][j] = SolnBlk.E[SolnBlk.ICl][j];
            SolnBlk.E[SolnBlk.ICl-2][j] = SolnBlk.E[SolnBlk.ICl][j];
            SolnBlk.V[SolnBlk.ICl-1][j] = SolnBlk.V[SolnBlk.ICl][j];
            SolnBlk.V[SolnBlk.ICl-2][j] = SolnBlk.V[SolnBlk.ICl][j];
            break;
        } /* endswitch */
      } /* endif */

      if ( (j >= SolnBlk.JCl && j <= SolnBlk.JCu) ||
           (j < SolnBlk.JCl && 
            (SolnBlk.Grid.BCtypeS[SolnBlk.ICu+1] == BC_NONE ||
             SolnBlk.Grid.BCtypeS[SolnBlk.ICu+1] == BC_PERIODIC ||
             SolnBlk.Grid.BCtypeS[SolnBlk.ICu+1] == BC_CONSTANT_EXTRAPOLATION ||
             SolnBlk.Grid.BCtypeS[SolnBlk.ICu+1] == BC_REFLECTION ||
             SolnBlk.Grid.BCtypeS[SolnBlk.ICu+1] == BC_FIXED_TEMP_WALL ||
             SolnBlk.Grid.BCtypeS[SolnBlk.ICu+1] == BC_ADIABATIC_WALL ||
             SolnBlk.Grid.BCtypeS[SolnBlk.ICu+1] == BC_LINEAR_EXTRAPOLATION ||
             SolnBlk.Grid.BCtypeS[SolnBlk.ICu+1] == BC_CHARACTERISTIC) ) ||
           (j > SolnBlk.JCu && 
            (SolnBlk.Grid.BCtypeN[SolnBlk.ICu+1] == BC_NONE ||
             SolnBlk.Grid.BCtypeN[SolnBlk.ICu+1] == BC_PERIODIC ||
             SolnBlk.Grid.BCtypeN[SolnBlk.ICu+1] == BC_CONSTANT_EXTRAPOLATION ||
             SolnBlk.Grid.BCtypeN[SolnBlk.ICu+1] == BC_REFLECTION ||
             SolnBlk.Grid.BCtypeN[SolnBlk.ICu+1] == BC_FIXED_TEMP_WALL ||
             SolnBlk.Grid.BCtypeN[SolnBlk.ICu+1] == BC_ADIABATIC_WALL ||
             SolnBlk.Grid.BCtypeN[SolnBlk.ICu+1] == BC_LINEAR_EXTRAPOLATION ||
             SolnBlk.Grid.BCtypeN[SolnBlk.ICu+1] == BC_CHARACTERISTIC) ) ) {
        switch(SolnBlk.Grid.BCtypeE[j]) {
          case BC_NONE :
          case BC_FIXED :
          case BC_CONSTANT_EXTRAPOLATION :
          case BC_LINEAR_EXTRAPOLATION :
          case BC_CHARACTERISTIC :
          case BC_REFLECTION :
          case BC_FIXED_TEMP_WALL :
          case BC_ADIABATIC_WALL :
            SolnBlk.E[SolnBlk.ICu+1][j] = SolnBlk.E[SolnBlk.ICu][j];
            SolnBlk.E[SolnBlk.ICu+2][j] = SolnBlk.E[SolnBlk.ICu][j];
            SolnBlk.V[SolnBlk.ICu+1][j] = SolnBlk.V[SolnBlk.ICu][j];
            SolnBlk.V[SolnBlk.ICu+2][j] = SolnBlk.V[SolnBlk.ICu][j];
            break;
          case BC_PERIODIC :
            SolnBlk.E[SolnBlk.ICu+1][j] = SolnBlk.E[SolnBlk.ICl+1][j];
            SolnBlk.E[SolnBlk.ICu+2][j] = SolnBlk.E[SolnBlk.ICl+2][j];
            SolnBlk.V[SolnBlk.ICu+1][j] = SolnBlk.V[SolnBlk.ICl+1][j];
            SolnBlk.V[SolnBlk.ICu+2][j] = SolnBlk.V[SolnBlk.ICl+2][j];
            break;
          default:
            SolnBlk.E[SolnBlk.ICu+1][j] = SolnBlk.E[SolnBlk.ICu][j];
            SolnBlk.E[SolnBlk.ICu+2][j] = SolnBlk.E[SolnBlk.ICu][j];
            SolnBlk.V[SolnBlk.ICu+1][j] = SolnBlk.V[SolnBlk.ICu][j];
            SolnBlk.V[SolnBlk.ICu+2][j] = SolnBlk.V[SolnBlk.ICu][j];
            break;
        } /* endswitch */
      } /* endif */
    } /* endfor */

    for ( i = SolnBlk.ICl-SolnBlk.Nghost ; i <= SolnBlk.ICu+SolnBlk.Nghost ; ++i ) {
        switch(SolnBlk.Grid.BCtypeS[i]) {
          case BC_NONE :
          case BC_FIXED :
          case BC_CONSTANT_EXTRAPOLATION :
          case BC_LINEAR_EXTRAPOLATION :
          case BC_CHARACTERISTIC :
          case BC_REFLECTION :
          case BC_FIXED_TEMP_WALL :
          case BC_ADIABATIC_WALL :
            SolnBlk.E[i][SolnBlk.JCl-1] = SolnBlk.E[i][SolnBlk.JCl];
            SolnBlk.E[i][SolnBlk.JCl-2] = SolnBlk.E[i][SolnBlk.JCl];
            SolnBlk.V[i][SolnBlk.JCl-1] = SolnBlk.V[i][SolnBlk.JCl];
            SolnBlk.V[i][SolnBlk.JCl-2] = SolnBlk.V[i][SolnBlk.JCl];
            break;
          case BC_PERIODIC :
            SolnBlk.E[i][SolnBlk.JCl-1] = SolnBlk.E[i][SolnBlk.JCu-1];
            SolnBlk.E[i][SolnBlk.JCl-2] = SolnBlk.E[i][SolnBlk.JCu-2];
            SolnBlk.V[i][SolnBlk.JCl-1] = SolnBlk.V[i][SolnBlk.JCu-1];
            SolnBlk.V[i][SolnBlk.JCl-2] = SolnBlk.V[i][SolnBlk.JCu-2];
            break;
          default:
            SolnBlk.E[i][SolnBlk.JCl-1] = SolnBlk.E[i][SolnBlk.JCl];
            SolnBlk.E[i][SolnBlk.JCl-2] = SolnBlk.E[i][SolnBlk.JCl];
            SolnBlk.V[i][SolnBlk.JCl-1] = SolnBlk.V[i][SolnBlk.JCl];
            SolnBlk.V[i][SolnBlk.JCl-2] = SolnBlk.V[i][SolnBlk.JCl];
            break;
        } /* endswitch */

        switch(SolnBlk.Grid.BCtypeN[i]) {
          case BC_NONE :
          case BC_FIXED :
          case BC_CONSTANT_EXTRAPOLATION :
          case BC_LINEAR_EXTRAPOLATION :
          case BC_CHARACTERISTIC :
          case BC_REFLECTION :
          case BC_FIXED_TEMP_WALL :
          case BC_ADIABATIC_WALL :
            SolnBlk.E[i][SolnBlk.JCu+1] = SolnBlk.E[i][SolnBlk.JCu];
            SolnBlk.E[i][SolnBlk.JCu+2] = SolnBlk.E[i][SolnBlk.JCu];
            SolnBlk.V[i][SolnBlk.JCu+1] = SolnBlk.V[i][SolnBlk.JCu];
            SolnBlk.V[i][SolnBlk.JCu+2] = SolnBlk.V[i][SolnBlk.JCu];
            break;
          case BC_PERIODIC :
            SolnBlk.E[i][SolnBlk.JCu+1] = SolnBlk.E[i][SolnBlk.JCl+1];
            SolnBlk.E[i][SolnBlk.JCu+2] = SolnBlk.E[i][SolnBlk.JCl+2];
            SolnBlk.V[i][SolnBlk.JCu+1] = SolnBlk.V[i][SolnBlk.JCl+1];
            SolnBlk.V[i][SolnBlk.JCu+2] = SolnBlk.V[i][SolnBlk.JCl+2];
            break;
          default:
            SolnBlk.E[i][SolnBlk.JCu+1] = SolnBlk.E[i][SolnBlk.JCu];
            SolnBlk.E[i][SolnBlk.JCu+2] = SolnBlk.E[i][SolnBlk.JCu];
            SolnBlk.V[i][SolnBlk.JCu+1] = SolnBlk.V[i][SolnBlk.JCu];
            SolnBlk.V[i][SolnBlk.JCu+2] = SolnBlk.V[i][SolnBlk.JCu];
            break;
        } /* endswitch */
    } /* endfor */

}

/********************************************************
 * Routine: CFL                                         *
 *                                                      *
 * Determines the allowable global and local time steps *
 * (for explicit Euler time stepping scheme) for the    *
 * specified quadrilateral solution block according to  *
 * the Courant-Friedrichs-Lewy condition for the        *
 * wave motion terms and a semi-impirical criteria for  *
 * the ion-neutral collisions source terms.             *
 *                                                      *
 ********************************************************/
double CFL(Ion5Moment2D_Quad_Block &SolnBlk,
           Ion5Moment2D_Input_Parameters &Input_Parameters) {

    int i, j;
    double dtMin, d_i, d_j, v_i, v_j, a, dt_cfl, dt_src;

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
  
	     dt_cfl = min(d_i/(a+fabs(v_i)), d_j/(a+fabs(v_j)));

             dt_src = 1.00e05/max(SolnBlk.W[i][j].nu(SolnBlk.Wneut[i][j]), TOLER*TOLER);
  
	     SolnBlk.dt[i][j] = min(dt_cfl, dt_src);

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
void Set_Global_TimeStep(Ion5Moment2D_Quad_Block &SolnBlk,
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
double L1_Norm_Residual(Ion5Moment2D_Quad_Block &SolnBlk) {

    int i, j;
    double l1_norm;

    l1_norm = ZERO;

    for ( j  = SolnBlk.JCl ; j <= SolnBlk.JCu ; ++j ) {
       for ( i = SolnBlk.ICl ; i <= SolnBlk.ICu ; ++i ) {
          l1_norm += abs(SolnBlk.dUdt[i][j][0].dv);
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
double L2_Norm_Residual(Ion5Moment2D_Quad_Block &SolnBlk) {

    int i, j;
    double l2_norm;

    l2_norm = ZERO;

    for ( j  = SolnBlk.JCl ; j <= SolnBlk.JCu ; ++j ) {
       for ( i = SolnBlk.ICl ; i <= SolnBlk.ICu ; ++i ) {
          l2_norm += sqr(SolnBlk.dUdt[i][j][0].dv);
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
double Max_Norm_Residual(Ion5Moment2D_Quad_Block &SolnBlk) {

    int i, j;
    double max_norm;

    max_norm = ZERO;

    for ( j  = SolnBlk.JCl ; j <= SolnBlk.JCu ; ++j ) {
       for ( i = SolnBlk.ICl ; i <= SolnBlk.ICu ; ++i ) {
          max_norm = max(max_norm, abs(SolnBlk.dUdt[i][j][0].dv));
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
void Linear_Reconstruction_GreenGauss(Ion5Moment2D_Quad_Block &SolnBlk,
				      const int i, 
                                      const int j,
                                      const int Limiter) {

    int n, n2, n_pts, i_index[8], j_index[8];
    double u0Min, u0Max, uQuad[4], phi;
    double DxDx_ave, DxDy_ave, DyDy_ave;
    double l_north, l_south, l_east, l_west;
    Vector2D n_north, n_south, n_east, n_west, dX;
    Ion5Moment2D_pState W_nw, W_ne, W_sw, W_se, W_face, 
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
           DUDx_ave = Ion5Moment2D_W_VACUUM;
           DUDy_ave = Ion5Moment2D_W_VACUUM;
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
        for ( n = 1 ; n <= NUM_VAR_ION5MOMENT2D ; ++n ) {
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

    } else {
        SolnBlk.dWdx[i][j] = Ion5Moment2D_W_VACUUM;
        SolnBlk.dWdy[i][j] = Ion5Moment2D_W_VACUUM; 
        SolnBlk.phi[i][j]  = Ion5Moment2D_W_VACUUM;
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
void Linear_Reconstruction_GreenGauss(Ion5Moment2D_Quad_Block &SolnBlk,
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
void Linear_Reconstruction_LeastSquares(Ion5Moment2D_Quad_Block &SolnBlk,
				        const int i, 
                                        const int j,
                                        const int Limiter) {

    int n, n2, n_pts, i_index[8], j_index[8];
    double u0Min, u0Max, uQuad[4], phi;
    double DxDx_ave, DxDy_ave, DyDy_ave;
    double Tratio, TratioMin, TratioMax;
    Vector2D dX;
    Ion5Moment2D_pState DU, DUDx_ave, DUDy_ave;

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
        DUDx_ave = Ion5Moment2D_W_VACUUM;
        DUDy_ave = Ion5Moment2D_W_VACUUM;
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
    
        for ( n = 1 ; n <= NUM_VAR_ION5MOMENT2D ; ++n ) {
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

        TratioMax = ONE;
        TratioMin = ONE;
        for ( n2 = 0 ; n2 <= n_pts-1 ; ++n2 ) {
            Tratio = SolnBlk.W[ i_index[n2] ][ j_index[n2] ].T()/SolnBlk.W[i][j].T();
            TratioMin = min(TratioMin, Tratio);
            TratioMax = max(TratioMax, Tratio);
        } /* endfor */
        if (TratioMin < 0.005 ||
            TratioMax > 200.00) {
	   //cout << i << " " << j << " " << TratioMin << " " << TratioMax << "\n";
	   if (TratioMax > 200.00) phi = max(TOLER, exp(-(TratioMax-200.00)/25.00));
           if (TratioMin < 0.005) phi = max(TOLER, exp(-(ONE/TratioMin-200.00)/25.00));
           SolnBlk.phi[i][j]  = phi*SolnBlk.phi[i][j];
        } /* endif */
    } else {
        SolnBlk.dWdx[i][j] = Ion5Moment2D_W_VACUUM;
        SolnBlk.dWdy[i][j] = Ion5Moment2D_W_VACUUM; 
        SolnBlk.phi[i][j]  = Ion5Moment2D_W_VACUUM;
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
void Linear_Reconstruction_LeastSquares_2(Ion5Moment2D_Quad_Block &SolnBlk,
				          const int i, 
                                          const int j,
                                          const int Limiter) {

    int n, n2, n_pts, i_index[8], j_index[8];
    double u0Min, u0Max, uQuad[4], phi;
    double DxDx_ave, DxDy_ave, DyDy_ave;
    double Tratio, TratioMin, TratioMax;
    Vector2D dX;
    Ion5Moment2D_pState DU, DUDx_ave, DUDy_ave;

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
        DUDx_ave = Ion5Moment2D_W_VACUUM;
        DUDy_ave = Ion5Moment2D_W_VACUUM;
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
    
        for ( n = 1 ; n <= NUM_VAR_ION5MOMENT2D ; ++n ) {
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

        TratioMax = ONE;
        TratioMin = ONE;
        for ( n2 = 0 ; n2 <= n_pts-1 ; ++n2 ) {
            Tratio = SolnBlk.W[ i_index[n2] ][ j_index[n2] ].T()/SolnBlk.W[i][j].T();
            TratioMin = min(TratioMin, Tratio);
            TratioMax = max(TratioMax, Tratio);
        } /* endfor */
        if (TratioMin < 0.005 ||
            TratioMax > 200.00) {
	   //cout << i << " " << j << " " << TratioMin << " " << TratioMax << "\n";
	   if (TratioMax > 200.00) phi = max(TOLER, exp(-(TratioMax-200.00)/25.00));
           if (TratioMin < 0.005) phi = max(TOLER, exp(-(ONE/TratioMin-200.00)/25.00));
           SolnBlk.phi[i][j]  = phi*SolnBlk.phi[i][j];
        } /* endif */
    } else {
        SolnBlk.dWdx[i][j] = Ion5Moment2D_W_VACUUM;
        SolnBlk.dWdy[i][j] = Ion5Moment2D_W_VACUUM; 
        SolnBlk.phi[i][j]  = Ion5Moment2D_W_VACUUM;
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
void Linear_Reconstruction_LeastSquares(Ion5Moment2D_Quad_Block &SolnBlk,
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
 * Routine: Residual_Smoothing                          *
 *                                                      *
 * Applies implicit residual smoothing to solution      *
 * residual.  Note that only residuals of interior cells*
 * are smoothed and residuals for cells adjacent to     *
 * boundaries are not smoothed.                         *
 *                                                      *
 ********************************************************/
void Residual_Smoothing(Ion5Moment2D_Quad_Block &SolnBlk,
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
				   Ion5Moment2D_Input_Parameters &IP,
                                   int &number_refinement_criteria,
                                   Ion5Moment2D_Quad_Block &SolnBlk) {

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

    /* Return the refinement criteria. */
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
void Fix_Refined_Block_Boundaries(Ion5Moment2D_Quad_Block SolnBlk,
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

	  SolnBlk.Wneut[i][SolnBlk.JCu] = W((SolnBlk.Grid.Cell[i][SolnBlk.JCu].A/
             SolnBlk.Grid.area(i, SolnBlk.JCu))*SolnBlk.Wneut[i][SolnBlk.JCu].U());
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

	  SolnBlk.Wneut[i][SolnBlk.JCl] = W((SolnBlk.Grid.Cell[i][SolnBlk.JCl].A/
             SolnBlk.Grid.area(i, SolnBlk.JCl))*SolnBlk.Wneut[i][SolnBlk.JCl].U());
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

	  SolnBlk.Wneut[SolnBlk.ICu][j] = W((SolnBlk.Grid.Cell[SolnBlk.ICu][j].A/
             SolnBlk.Grid.area(SolnBlk.ICu, j))*SolnBlk.Wneut[SolnBlk.ICu][j].U());
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

	  SolnBlk.Wneut[SolnBlk.ICl][j] = W((SolnBlk.Grid.Cell[SolnBlk.ICl][j].A/
             SolnBlk.Grid.area(SolnBlk.ICl, j))*SolnBlk.Wneut[SolnBlk.ICl][j].U());
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
void Unfix_Refined_Block_Boundaries(Ion5Moment2D_Quad_Block SolnBlk) {

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

	  SolnBlk.Wneut[i][SolnBlk.JCu] = W((SolnBlk.Grid.Cell[i][SolnBlk.JCu].A/
             SolnBlk.Grid.area(i, SolnBlk.JCu))*SolnBlk.Wneut[i][SolnBlk.JCu].U());
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

	  SolnBlk.Wneut[i][SolnBlk.JCl] = W((SolnBlk.Grid.Cell[i][SolnBlk.JCl].A/
             SolnBlk.Grid.area(i, SolnBlk.JCl))*SolnBlk.Wneut[i][SolnBlk.JCl].U());
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

	  SolnBlk.Wneut[SolnBlk.ICu][j] = W((SolnBlk.Grid.Cell[SolnBlk.ICu][j].A/
             SolnBlk.Grid.area(SolnBlk.ICu, j))*SolnBlk.Wneut[SolnBlk.ICu][j].U());
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

	  SolnBlk.Wneut[SolnBlk.ICl][j] = W((SolnBlk.Grid.Cell[SolnBlk.ICl][j].A/
             SolnBlk.Grid.area(SolnBlk.ICl, j))*SolnBlk.Wneut[SolnBlk.ICl][j].U());
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
void Apply_Boundary_Flux_Corrections(Ion5Moment2D_Quad_Block SolnBlk,
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
void Apply_Boundary_Flux_Corrections_Multistage_Explicit(Ion5Moment2D_Quad_Block SolnBlk,
                                                         const int i_stage,
                                                         const int n_stage,
	                                                 const double &CFL_Number,
                                                         const int Time_Integration_Type,
                                                         const int Local_Time_Stepping,
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
 * solution block using a 2nd-ororder limited upwind    *
 * finite-volume spatial discretization scheme with     *
 * either the Godunov, Roe, Rusanov, HLLE, Linde, or    *
 * HLLC flux functions.                                 *
 * The residual is stored in dUdt[][][0].               *
 *                                                      *
 ********************************************************/
int dUdt_Residual_Evaluation(Ion5Moment2D_Quad_Block &SolnBlk,
			     Ion5Moment2D_Input_Parameters &Input_Parameters) {

    int i, j;
    Vector2D dX;
    Ion5Moment2D_pState Wl, Wr;
    Ion5Moment2D_cState Flux;

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
      SolnBlk.dUdt[SolnBlk.ICl-1][j][0] = Ion5Moment2D_U_VACUUM;
    
       for ( i = SolnBlk.ICl-1 ; i <= SolnBlk.ICu ; ++i ) {
	 SolnBlk.dUdt[i+1][j][0] = Ion5Moment2D_U_VACUUM;
    
	 if ( j > SolnBlk.JCl-1 && j < SolnBlk.JCu+1 ) {
    
             /* Evaluate the cell interface i-direction fluxes. */
    
	     if (i == SolnBlk.ICl-1 && 
                 (SolnBlk.Grid.BCtypeW[j] == BC_REFLECTION ||
                  SolnBlk.Grid.BCtypeW[j] == BC_CHARACTERISTIC)) {
	       dX = SolnBlk.Grid.xfaceW(i+1, j)-SolnBlk.Grid.Cell[i+1][j].Xc;
	       Wr = SolnBlk.W[i+1][j] + 
		 (SolnBlk.phi[i+1][j]^SolnBlk.dWdx[i+1][j])*dX.x +
		 (SolnBlk.phi[i+1][j]^SolnBlk.dWdy[i+1][j])*dX.y;
               if (Wr.d <= ZERO || Wr.p <= ZERO) Wr = SolnBlk.W[i+1][j];
	       if (SolnBlk.Grid.BCtypeW[j] == BC_REFLECTION) {
                 Wl = Reflect(Wr, SolnBlk.Grid.nfaceW(i+1, j));
	       } else {
		 Wl = BC_Characteristic_Pressure(Wr, 
						 SolnBlk.WoW[j], 
						 SolnBlk.Grid.nfaceW(i+1, j));
                } /* endif */
             } else if (i == SolnBlk.ICu && 
                        (SolnBlk.Grid.BCtypeE[j] == BC_REFLECTION ||
                         SolnBlk.Grid.BCtypeE[j] == BC_CHARACTERISTIC)) {
	       dX = SolnBlk.Grid.xfaceE(i, j)-SolnBlk.Grid.Cell[i][j].Xc;
	       Wl = SolnBlk.W[i][j] + 
		 (SolnBlk.phi[i][j]^SolnBlk.dWdx[i][j])*dX.x +
		 (SolnBlk.phi[i][j]^SolnBlk.dWdy[i][j])*dX.y;
               if (Wl.d <= ZERO || Wl.p <= ZERO) Wl = SolnBlk.W[i][j];
	       if (SolnBlk.Grid.BCtypeE[j] == BC_REFLECTION) {
		 Wr = Reflect(Wl, SolnBlk.Grid.nfaceE(i, j));
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
                if (Wl.d <= ZERO || Wl.p <= ZERO) Wl = SolnBlk.W[i][j];
	        dX = SolnBlk.Grid.xfaceW(i+1, j)-SolnBlk.Grid.Cell[i+1][j].Xc;
                Wr = SolnBlk.W[i+1][j] + 
  	             (SolnBlk.phi[i+1][j]^SolnBlk.dWdx[i+1][j])*dX.x +
	             (SolnBlk.phi[i+1][j]^SolnBlk.dWdy[i+1][j])*dX.y;
                if (Wr.d <= ZERO || Wr.p <= ZERO) Wr = SolnBlk.W[i+1][j];
  	     } /* endif */

             switch(Input_Parameters.i_Flux_Function) {
               case FLUX_FUNCTION_GODUNOV :
                 Flux = FluxGodunov_n(Wl, Wr, SolnBlk.Grid.nfaceE(i, j));
                 break;
               case FLUX_FUNCTION_ROE :
                 Flux = FluxRoe_n(Wl, Wr, SolnBlk.Grid.nfaceE(i, j));
                 break;
               case FLUX_FUNCTION_RUSANOV :
                 Flux = FluxRusanov_n(Wl, Wr, SolnBlk.Grid.nfaceE(i, j));
                 break;
               case FLUX_FUNCTION_HLLE :
                 Flux = FluxHLLE_n(Wl, Wr, SolnBlk.Grid.nfaceE(i, j));
                 break;
               case FLUX_FUNCTION_LINDE :
                 Flux = FluxLinde_n(Wl, Wr, SolnBlk.Grid.nfaceE(i, j));
                 break;
               case FLUX_FUNCTION_HLLC :
                 Flux = FluxHLLC_n(Wl, Wr, SolnBlk.Grid.nfaceE(i, j));
                 break;
               default:
                 Flux = FluxRoe_n(Wl, Wr, SolnBlk.Grid.nfaceE(i, j));
                 break;
             } /* endswitch */
    
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
	          Sa(SolnBlk.W[i][j], SolnBlk.Grid.Cell[i][j].Xc);
             } /* endif */

             /* Include source terms associated with the 
                electric field. */

             SolnBlk.dUdt[i][j][0] += 
	       Se(SolnBlk.W[i][j], SolnBlk.E[i][j]);

             /* Include source terms associated with the 
                ion-neutral collision processes. */

             SolnBlk.dUdt[i][j][0] += 
	       Sc(SolnBlk.W[i][j], SolnBlk.Wneut[i][j]);

             /* Save west and east face boundary flux. */

             if (i == SolnBlk.ICl-1) {
                SolnBlk.FluxW[j] = -Flux*SolnBlk.Grid.lfaceW(i+1, j);
             } else if (i == SolnBlk.ICu) {
                SolnBlk.FluxE[j] = Flux*SolnBlk.Grid.lfaceE(i, j);
             } /* endif */

          } /* endif */
       } /* endfor */
    
       if ( j > SolnBlk.JCl-1 && j < SolnBlk.JCu+1 ) {
          SolnBlk.dUdt[SolnBlk.ICl-1][j][0] = Ion5Moment2D_U_VACUUM;
          SolnBlk.dUdt[SolnBlk.ICu+1][j][0] = Ion5Moment2D_U_VACUUM;
       } /* endif */
    } /* endfor */
    
    // Add j-direction (eta-direction) fluxes.
    for ( i = SolnBlk.ICl ; i <= SolnBlk.ICu ; ++i ) {
       for ( j  = SolnBlk.JCl-1 ; j <= SolnBlk.JCu ; ++j ) {
    
          /* Evaluate the cell interface j-direction fluxes. */
         
	  if (j == SolnBlk.JCl-1 && 
              (SolnBlk.Grid.BCtypeS[i] == BC_REFLECTION ||
               SolnBlk.Grid.BCtypeS[i] == BC_CHARACTERISTIC)) {
	     dX = SolnBlk.Grid.xfaceS(i, j+1)-SolnBlk.Grid.Cell[i][j+1].Xc;
	     Wr = SolnBlk.W[i][j+1] +
	          (SolnBlk.phi[i][j+1]^SolnBlk.dWdx[i][j+1])*dX.x +
	          (SolnBlk.phi[i][j+1]^SolnBlk.dWdy[i][j+1])*dX.y;
             if (Wr.d <= ZERO || Wr.p <= ZERO) Wr = SolnBlk.W[i][j+1];
             if (SolnBlk.Grid.BCtypeS[i] == BC_REFLECTION) {
               Wl = Reflect(Wr, SolnBlk.Grid.nfaceS(i, j+1));
             } else {
               Wl = BC_Characteristic_Pressure(Wr, 
                                               SolnBlk.WoS[i], 
                                               SolnBlk.Grid.nfaceS(i, j+1));
             } /* endif */
          } else if (j == SolnBlk.JCu && 
                     (SolnBlk.Grid.BCtypeN[i] == BC_REFLECTION ||
                      SolnBlk.Grid.BCtypeN[i] == BC_CHARACTERISTIC)) {
	    dX = SolnBlk.Grid.xfaceN(i, j)-SolnBlk.Grid.Cell[i][j].Xc;
	    Wl = SolnBlk.W[i][j] + 
	      (SolnBlk.phi[i][j]^SolnBlk.dWdx[i][j])*dX.x +
	      (SolnBlk.phi[i][j]^SolnBlk.dWdy[i][j])*dX.y;
            if (Wl.d <= ZERO || Wl.p <= ZERO) Wl = SolnBlk.W[i][j];
	    if (SolnBlk.Grid.BCtypeN[i] == BC_REFLECTION) {
	      Wr = Reflect(Wl, SolnBlk.Grid.nfaceN(i, j));
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
             if (Wl.d <= ZERO || Wl.p <= ZERO) Wl = SolnBlk.W[i][j];
	     dX = SolnBlk.Grid.xfaceS(i, j+1)-SolnBlk.Grid.Cell[i][j+1].Xc;
             Wr = SolnBlk.W[i][j+1] +
                  (SolnBlk.phi[i][j+1]^SolnBlk.dWdx[i][j+1])*dX.x +
                  (SolnBlk.phi[i][j+1]^SolnBlk.dWdy[i][j+1])*dX.y;
             if (Wr.d <= ZERO || Wr.p <= ZERO) Wr = SolnBlk.W[i][j+1];
          } /* endif */

          switch(Input_Parameters.i_Flux_Function) {
            case FLUX_FUNCTION_GODUNOV :
              Flux = FluxGodunov_n(Wl, Wr, SolnBlk.Grid.nfaceN(i, j));
              break;
            case FLUX_FUNCTION_ROE :
              Flux = FluxRoe_n(Wl, Wr, SolnBlk.Grid.nfaceN(i, j));
              break;
            case FLUX_FUNCTION_RUSANOV :
              Flux = FluxRusanov_n(Wl, Wr, SolnBlk.Grid.nfaceN(i, j));
              break;
            case FLUX_FUNCTION_HLLE :
              Flux = FluxHLLE_n(Wl, Wr, SolnBlk.Grid.nfaceN(i, j));
              break;
            case FLUX_FUNCTION_LINDE :
              Flux = FluxLinde_n(Wl, Wr, SolnBlk.Grid.nfaceN(i, j));
              break;
            case FLUX_FUNCTION_HLLC :
              Flux = FluxHLLC_n(Wl, Wr, SolnBlk.Grid.nfaceN(i, j));
              break;
            default:
              Flux = FluxRoe_n(Wl, Wr, SolnBlk.Grid.nfaceN(i, j));
              break;
	  } /* endswitch */
    
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
    
       SolnBlk.dUdt[i][SolnBlk.JCl-1][0] = Ion5Moment2D_U_VACUUM;
       SolnBlk.dUdt[i][SolnBlk.JCu+1][0] = Ion5Moment2D_U_VACUUM;
    } /* endfor */
    
    /* Residual evaluated successfully. */
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
int dUdt_Multistage_Explicit(Ion5Moment2D_Quad_Block &SolnBlk,
                             const int i_stage,
                             Ion5Moment2D_Input_Parameters &Input_Parameters) {

    int i, j, k_residual;
    double omega;
    Vector2D dX;
    Ion5Moment2D_pState Wl, Wr;
    Ion5Moment2D_cState Flux;

    /* Evaluate the solution residual for stage 
       i_stage of an N stage scheme. */

    /* Evaluate the time step fraction for the stage. */
    
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
          SolnBlk.dUdt[SolnBlk.ICl-1][j][k_residual] = Ion5Moment2D_U_VACUUM;
       } else {
          SolnBlk.dUdt[SolnBlk.ICl-1][j][k_residual] = Ion5Moment2D_U_VACUUM;
       } /* endif */
    
       for ( i = SolnBlk.ICl-1 ; i <= SolnBlk.ICu ; ++i ) {
          if ( i_stage == 1 ) {
              SolnBlk.Uo[i+1][j] = SolnBlk.U[i+1][j];
              SolnBlk.dUdt[i+1][j][k_residual] = Ion5Moment2D_U_VACUUM;
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
                     SolnBlk.dUdt[i+1][j][k_residual] = Ion5Moment2D_U_VACUUM;
                  } /* endif */
                  break;
                case TIME_STEPPING_MULTISTAGE_OPTIMAL_SMOOTHING :
                  SolnBlk.dUdt[i+1][j][k_residual] = Ion5Moment2D_U_VACUUM;
                  break;
                default:
                  SolnBlk.dUdt[i+1][j][k_residual] = Ion5Moment2D_U_VACUUM;
                  break;
              } /* endswitch */
          } /* endif */
    
          if ( j > SolnBlk.JCl-1 && j < SolnBlk.JCu+1 ) {
    
             /* Evaluate the cell interface i-direction fluxes. */
    
	     if (i == SolnBlk.ICl-1 && 
                 (SolnBlk.Grid.BCtypeW[j] == BC_REFLECTION ||
                  SolnBlk.Grid.BCtypeW[j] == BC_CHARACTERISTIC)) {
	       dX = SolnBlk.Grid.xfaceW(i+1, j)-SolnBlk.Grid.Cell[i+1][j].Xc;
	       Wr = SolnBlk.W[i+1][j] + 
		 (SolnBlk.phi[i+1][j]^SolnBlk.dWdx[i+1][j])*dX.x +
		 (SolnBlk.phi[i+1][j]^SolnBlk.dWdy[i+1][j])*dX.y;
               if (Wr.d <= ZERO || Wr.p <= ZERO) Wr = SolnBlk.W[i+1][j];
	       if (SolnBlk.Grid.BCtypeW[j] == BC_REFLECTION) {
                 Wl = Reflect(Wr, SolnBlk.Grid.nfaceW(i+1, j));
	       } else {
		 Wl = BC_Characteristic_Pressure(Wr, 
						 SolnBlk.WoW[j], 
						 SolnBlk.Grid.nfaceW(i+1, j));
                } /* endif */
             } else if (i == SolnBlk.ICu && 
                        (SolnBlk.Grid.BCtypeE[j] == BC_REFLECTION ||
                         SolnBlk.Grid.BCtypeE[j] == BC_CHARACTERISTIC)) {
	       dX = SolnBlk.Grid.xfaceE(i, j)-SolnBlk.Grid.Cell[i][j].Xc;
	       Wl = SolnBlk.W[i][j] + 
		 (SolnBlk.phi[i][j]^SolnBlk.dWdx[i][j])*dX.x +
		 (SolnBlk.phi[i][j]^SolnBlk.dWdy[i][j])*dX.y;
               if (Wl.d <= ZERO || Wl.p <= ZERO) Wl = SolnBlk.W[i][j];
	       if (SolnBlk.Grid.BCtypeE[j] == BC_REFLECTION) {
		 Wr = Reflect(Wl, SolnBlk.Grid.nfaceE(i, j));
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
                if (Wl.d <= ZERO || Wl.p <= ZERO) Wl = SolnBlk.W[i][j];
	        dX = SolnBlk.Grid.xfaceW(i+1, j)-SolnBlk.Grid.Cell[i+1][j].Xc;
                Wr = SolnBlk.W[i+1][j] + 
  	             (SolnBlk.phi[i+1][j]^SolnBlk.dWdx[i+1][j])*dX.x +
	             (SolnBlk.phi[i+1][j]^SolnBlk.dWdy[i+1][j])*dX.y;
                if (Wr.d <= ZERO || Wr.p <= ZERO) Wr = SolnBlk.W[i+1][j];
  	     } /* endif */

             switch(Input_Parameters.i_Flux_Function) {
               case FLUX_FUNCTION_GODUNOV :
                 Flux = FluxGodunov_n(Wl, Wr, SolnBlk.Grid.nfaceE(i, j));
                 break;
               case FLUX_FUNCTION_ROE :
                 Flux = FluxRoe_n(Wl, Wr, SolnBlk.Grid.nfaceE(i, j));
                 break;
               case FLUX_FUNCTION_RUSANOV :
                 Flux = FluxRusanov_n(Wl, Wr, SolnBlk.Grid.nfaceE(i, j));
                 break;
               case FLUX_FUNCTION_HLLE :
                 Flux = FluxHLLE_n(Wl, Wr, SolnBlk.Grid.nfaceE(i, j));
                 break;
               case FLUX_FUNCTION_LINDE :
                 Flux = FluxLinde_n(Wl, Wr, SolnBlk.Grid.nfaceE(i, j));
                 break;
               case FLUX_FUNCTION_HLLC :
                 Flux = FluxHLLC_n(Wl, Wr, SolnBlk.Grid.nfaceE(i, j));
                 break;
               default:
                 Flux = FluxRoe_n(Wl, Wr, SolnBlk.Grid.nfaceE(i, j));
                 break;
             } /* endswitch */
    
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
	          Sa(SolnBlk.W[i][j], SolnBlk.Grid.Cell[i][j].Xc);
             } /* endif */

             /* Include source terms associated with the 
                electric field. */

             SolnBlk.dUdt[i][j][k_residual] += 
               (Input_Parameters.CFL_Number*SolnBlk.dt[i][j])*
	       Se(SolnBlk.W[i][j], SolnBlk.E[i][j]);

             /* Include source terms associated with the 
                ion-neutral collision processes. */

             SolnBlk.dUdt[i][j][k_residual] += 
               (Input_Parameters.CFL_Number*SolnBlk.dt[i][j])*
	       Sc(SolnBlk.W[i][j], SolnBlk.Wneut[i][j]);

             /* Save west and east face boundary flux. */

             if (i == SolnBlk.ICl-1) {
                SolnBlk.FluxW[j] = -Flux*SolnBlk.Grid.lfaceW(i+1, j);
             } else if (i == SolnBlk.ICu) {
                SolnBlk.FluxE[j] = Flux*SolnBlk.Grid.lfaceE(i, j);
             } /* endif */ 

          } /* endif */
       } /* endfor */
    
       if ( j > SolnBlk.JCl-1 && j < SolnBlk.JCu+1 ) {
          SolnBlk.dUdt[SolnBlk.ICl-1][j][k_residual] = Ion5Moment2D_U_VACUUM;
          SolnBlk.dUdt[SolnBlk.ICu+1][j][k_residual] = Ion5Moment2D_U_VACUUM;
       } /* endif */
    } /* endfor */
    
    // Add j-direction (eta-direction) fluxes.
    for ( i = SolnBlk.ICl ; i <= SolnBlk.ICu ; ++i ) {
       for ( j  = SolnBlk.JCl-1 ; j <= SolnBlk.JCu ; ++j ) {
    
          /* Evaluate the cell interface j-direction fluxes. */
         
	  if (j == SolnBlk.JCl-1 && 
              (SolnBlk.Grid.BCtypeS[i] == BC_REFLECTION ||
               SolnBlk.Grid.BCtypeS[i] == BC_CHARACTERISTIC)) {
	     dX = SolnBlk.Grid.xfaceS(i, j+1)-SolnBlk.Grid.Cell[i][j+1].Xc;
	     Wr = SolnBlk.W[i][j+1] +
	          (SolnBlk.phi[i][j+1]^SolnBlk.dWdx[i][j+1])*dX.x +
	          (SolnBlk.phi[i][j+1]^SolnBlk.dWdy[i][j+1])*dX.y;
             if (Wr.d <= ZERO || Wr.p <= ZERO) Wr = SolnBlk.W[i][j+1];
             if (SolnBlk.Grid.BCtypeS[i] == BC_REFLECTION) {
               Wl = Reflect(Wr, SolnBlk.Grid.nfaceS(i, j+1));
             } else {
               Wl = BC_Characteristic_Pressure(Wr, 
                                               SolnBlk.WoS[i], 
                                               SolnBlk.Grid.nfaceS(i, j+1));
             } /* endif */
          } else if (j == SolnBlk.JCu && 
                     (SolnBlk.Grid.BCtypeN[i] == BC_REFLECTION ||
                      SolnBlk.Grid.BCtypeN[i] == BC_CHARACTERISTIC)) {
	    dX = SolnBlk.Grid.xfaceN(i, j)-SolnBlk.Grid.Cell[i][j].Xc;
	    Wl = SolnBlk.W[i][j] + 
	      (SolnBlk.phi[i][j]^SolnBlk.dWdx[i][j])*dX.x +
	      (SolnBlk.phi[i][j]^SolnBlk.dWdy[i][j])*dX.y;
            if (Wl.d <= ZERO || Wl.p <= ZERO) Wl = SolnBlk.W[i][j];
	    if (SolnBlk.Grid.BCtypeN[i] == BC_REFLECTION) {
	      Wr = Reflect(Wl, SolnBlk.Grid.nfaceN(i, j));
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
             if (Wl.d <= ZERO || Wl.p <= ZERO) Wl = SolnBlk.W[i][j];
	     dX = SolnBlk.Grid.xfaceS(i, j+1)-SolnBlk.Grid.Cell[i][j+1].Xc;
             Wr = SolnBlk.W[i][j+1] +
                  (SolnBlk.phi[i][j+1]^SolnBlk.dWdx[i][j+1])*dX.x +
                  (SolnBlk.phi[i][j+1]^SolnBlk.dWdy[i][j+1])*dX.y;
             if (Wr.d <= ZERO || Wr.p <= ZERO) Wr = SolnBlk.W[i][j+1];
          } /* endif */

          switch(Input_Parameters.i_Flux_Function) {
            case FLUX_FUNCTION_GODUNOV :
              Flux = FluxGodunov_n(Wl, Wr, SolnBlk.Grid.nfaceN(i, j));
              break;
            case FLUX_FUNCTION_ROE :
              Flux = FluxRoe_n(Wl, Wr, SolnBlk.Grid.nfaceN(i, j));
              break;
            case FLUX_FUNCTION_RUSANOV :
              Flux = FluxRusanov_n(Wl, Wr, SolnBlk.Grid.nfaceN(i, j));
              break;
            case FLUX_FUNCTION_HLLE :
              Flux = FluxHLLE_n(Wl, Wr, SolnBlk.Grid.nfaceN(i, j));
              break;
            case FLUX_FUNCTION_LINDE :
              Flux = FluxLinde_n(Wl, Wr, SolnBlk.Grid.nfaceN(i, j));
              break;
            case FLUX_FUNCTION_HLLC :
              Flux = FluxHLLC_n(Wl, Wr, SolnBlk.Grid.nfaceN(i, j));
              break;
            default:
              Flux = FluxRoe_n(Wl, Wr, SolnBlk.Grid.nfaceN(i, j));
              break;
	  } /* endswitch */
    
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
    
       SolnBlk.dUdt[i][SolnBlk.JCl-1][k_residual] = Ion5Moment2D_U_VACUUM;
       SolnBlk.dUdt[i][SolnBlk.JCu+1][k_residual] = Ion5Moment2D_U_VACUUM;
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
int Update_Solution_Multistage_Explicit(Ion5Moment2D_Quad_Block &SolnBlk,
                                        const int i_stage,
                                        Ion5Moment2D_Input_Parameters &Input_Parameters) {

    int i, j, k, k_residual, n_residual_reduction;
    double omega, residual_reduction_factor;

    // Memory for linear system solver.
    DenseMatrix dSdU(NUM_VAR_ION5MOMENT2D,NUM_VAR_ION5MOMENT2D);
    DenseSystemLinEqs LinSys;

    /* Allocate memory for linear system solver. */

    LinSys.allocate(NUM_VAR_ION5MOMENT2D);

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

	  /* Evaluate the source Jacobian. */ 

	  dSdU.zero();

	  if (SolnBlk.Axisymmetric) { // Calculate source Jacobian for axisymmetric terms
	     dSadU(dSdU, SolnBlk.W[i][j], SolnBlk.Grid.Cell[i][j].Xc);
          } /* endif */

	  dSedU(dSdU, SolnBlk.W[i][j], SolnBlk.E[i][j]); // Calculate source Jacobian for electric field.

	  dScdU(dSdU, SolnBlk.W[i][j], SolnBlk.Wneut[i][j]); // Calculate source Jacobian for ion-neutral collisions.

	  /* Point implicit formulation: set up system of equations 
             and include source Jacobian in the LHS matrix. */

          LinSys.A.identity();

          LinSys.A -= (omega*Input_Parameters.CFL_Number*SolnBlk.dt[i][j])*dSdU;

          for (k = 1; k <= NUM_VAR_ION5MOMENT2D; ++k) {
             LinSys.b(k-1) = omega*SolnBlk.dUdt[i][j][k_residual][k];
          } /* endfor */
 
          /* Solve system of equations using LU decomposition
             Gaussian elimination procedure. */

          LinSys.solve(LU_DECOMPOSITION);

          /* Update the conserved solution variables. */

          for (k = 1; k <= NUM_VAR_ION5MOMENT2D; ++k) {
             SolnBlk.U[i][j][k] = SolnBlk.Uo[i][j][k] + LinSys.x(k-1);
          } /* endfor */
     
          if (Input_Parameters.Local_Time_Stepping == SCALAR_LOCAL_TIME_STEPPING && 
              (SolnBlk.U[i][j].d   <= ZERO ||
               SolnBlk.U[i][j].E   <= ZERO ||
               SolnBlk.U[i][j].e() <= ZERO)) {
	     residual_reduction_factor = ONE;
             for (n_residual_reduction = 1; n_residual_reduction <= 10; ++n_residual_reduction) {
                residual_reduction_factor = HALF*residual_reduction_factor;
                SolnBlk.dt[i][j] = residual_reduction_factor*SolnBlk.dt[i][j];
                SolnBlk.dUdt[i][j][k_residual] = residual_reduction_factor*SolnBlk.dUdt[i][j][k_residual];
                LinSys.A.identity();
                LinSys.A -= (omega*Input_Parameters.CFL_Number*SolnBlk.dt[i][j])*dSdU;
                for (k = 1; k <= NUM_VAR_ION5MOMENT2D; ++k) {
                   LinSys.b(k-1) = omega*SolnBlk.dUdt[i][j][k_residual][k];
                } /* endfor */
                LinSys.solve(LU_DECOMPOSITION);
                for (k = 1; k <= NUM_VAR_ION5MOMENT2D; ++k) {
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
	      cout << "\n " << CFFC_Name() << " Ion5Moment2D ERROR: Negative Density and/or Energy: \n"
                   << " cell = (" << i << ", " << j << ") "
                   << " X = " << SolnBlk.Grid.Cell[i][j].Xc << "\n U = " 
                   << SolnBlk.U[i][j] << "\n dUdt = " 
                   << SolnBlk.dUdt[i][j][k_residual] << "\n";
              return (i);
          }

          /* Update the primitive solution variables. */
	      
          SolnBlk.W[i][j] = W(SolnBlk.U[i][j]);

       } /* endfor */    
    } /* endfor */

    /* Deallocate memory for linear system solver. */

    LinSys.deallocate();

    /* Solution successfully updated. */

    return (0);
    
}

/********************************************************
 * Routine: dUdt_Output_Cells_Tecplot                   *
 *                                                      *
 * This routine outputs the various terms that are      *
 * contributing to the solution changes (inertial,      *
 * axisymmetric source, electric field, and neutral     *
 * collisional processes) at the cell centres.          *
 *                                                      *
 ********************************************************/
void dUdt_Output_Cells_Tecplot(Ion5Moment2D_Quad_Block &SolnBlk,
                               const int Number_of_Time_Steps,
                               const double &Time,
                               const int Block_Number,
                               const int Output_Title,
                               const int Reconstruction_Type,
                               const int Limiter_Type,
	                       const int Flux_Function_Type,
	                       ostream &Out_File) {

    int i, j, k, n_residual_reduction;
    double omega, residual_reduction_factor;
    Vector2D dX;
    Ion5Moment2D_pState Wl, Wr;
    Ion5Moment2D_cState Flux;
    Ion5Moment2D_cState dUdt_inertia, dUdt_axi, dUdt_efield, dUdt_collisions;

    Out_File << setprecision(14);
    if (Output_Title) {
       Out_File << "TITLE = \"" << CFFC_Name() << ": 2D 5-Moment Ion Transport Solution Changes, "
                << "Time Step/Iteration Level = " << Number_of_Time_Steps
                << ", Time = " << Time
                << "\"" << "\n"
   	        << "VARIABLES = \"x\" \\ \n"
                << "\"y\" \\ \n"
                << "\"dU1idt\" \\ \n"
                << "\"dU2idt\" \\ \n"
                << "\"dU3idt\" \\ \n"
                << "\"dU4idt\" \\ \n"
                << "\"dU1adt\" \\ \n"
                << "\"dU2adt\" \\ \n"
                << "\"dU3adt\" \\ \n"
                << "\"dU4adt\" \\ \n"
                << "\"dU1edt\" \\ \n"
                << "\"dU2edt\" \\ \n"
                << "\"dU3edt\" \\ \n"
                << "\"dU4edt\" \\ \n"
                << "\"dU1cdt\" \\ \n"
                << "\"dU2cdt\" \\ \n"
                << "\"dU3cdt\" \\ \n"
                << "\"dU4cdt\" \\ \n"
                << "ZONE T =  \"Block Number = " << Block_Number
                << "\" \\ \n"
                << "I = " << SolnBlk.Grid.ICu - SolnBlk.Grid.ICl + 1 << " \\ \n"
                << "J = " << SolnBlk.Grid.JCu - SolnBlk.Grid.JCl + 1 << " \\ \n"
                << "F = POINT \n";
    } else {
       Out_File << "ZONE T =  \"Block Number = " << Block_Number
                << "\" \\ \n"
                << "I = " << SolnBlk.Grid.ICu - SolnBlk.Grid.ICl + 1 << " \\ \n"
                << "J = " << SolnBlk.Grid.JCu - SolnBlk.Grid.JCl + 1 << " \\ \n"
                << "F = POINT \n";
    } /* endif */

    /* Perform the linear reconstruction within each cell
       of the computational grid for this stage. */
    
    switch(Reconstruction_Type) {
    case RECONSTRUCTION_GREEN_GAUSS :
      Linear_Reconstruction_GreenGauss(SolnBlk,
                                       Limiter_Type);    
      break;
    case RECONSTRUCTION_LEAST_SQUARES :
      Linear_Reconstruction_LeastSquares(SolnBlk,
					 Limiter_Type);
      break;
    default:
      Linear_Reconstruction_LeastSquares(SolnBlk,
                                         Limiter_Type);
      break;
    } /* endswitch */

    /* Evaluate the time rate of change of the solution
       (i.e., the solution residuals) using a second-order
       limited upwind scheme with a variety of flux functions. */
    
    // Add i-direction (zeta-direction) fluxes.
    for ( j  = SolnBlk.JCl-1 ; j <= SolnBlk.JCu+1 ; ++j ) {
       SolnBlk.dUdt[SolnBlk.ICl-1][j][0] = Ion5Moment2D_U_VACUUM;
    
       for ( i = SolnBlk.ICl-1 ; i <= SolnBlk.ICu ; ++i ) {
          SolnBlk.dUdt[i+1][j][0] = Ion5Moment2D_U_VACUUM;
    
          if ( j > SolnBlk.JCl-1 && j < SolnBlk.JCu+1 ) {
    
             /* Evaluate the cell interface i-direction fluxes. */
    
	     if (i == SolnBlk.ICl-1 && 
                 (SolnBlk.Grid.BCtypeW[j] == BC_REFLECTION ||
                  SolnBlk.Grid.BCtypeW[j] == BC_CHARACTERISTIC)) {
	       dX = SolnBlk.Grid.xfaceW(i+1, j)-SolnBlk.Grid.Cell[i+1][j].Xc;
	       Wr = SolnBlk.W[i+1][j] + 
		 (SolnBlk.phi[i+1][j]^SolnBlk.dWdx[i+1][j])*dX.x +
		 (SolnBlk.phi[i+1][j]^SolnBlk.dWdy[i+1][j])*dX.y;
               if (Wr.d <= ZERO || Wr.p <= ZERO) Wr = SolnBlk.W[i+1][j];
	       if (SolnBlk.Grid.BCtypeW[j] == BC_REFLECTION) {
                 Wl = Reflect(Wr, SolnBlk.Grid.nfaceW(i+1, j));
	       } else {
		 Wl = BC_Characteristic_Pressure(Wr, 
						 SolnBlk.WoW[j], 
						 SolnBlk.Grid.nfaceW(i+1, j));
                } /* endif */
             } else if (i == SolnBlk.ICu && 
                        (SolnBlk.Grid.BCtypeE[j] == BC_REFLECTION ||
                         SolnBlk.Grid.BCtypeE[j] == BC_CHARACTERISTIC)) {
	       dX = SolnBlk.Grid.xfaceE(i, j)-SolnBlk.Grid.Cell[i][j].Xc;
	       Wl = SolnBlk.W[i][j] + 
		 (SolnBlk.phi[i][j]^SolnBlk.dWdx[i][j])*dX.x +
		 (SolnBlk.phi[i][j]^SolnBlk.dWdy[i][j])*dX.y;
               if (Wl.d <= ZERO || Wl.p <= ZERO) Wl = SolnBlk.W[i][j];
	       if (SolnBlk.Grid.BCtypeE[j] == BC_REFLECTION) {
		 Wr = Reflect(Wl, SolnBlk.Grid.nfaceE(i, j));
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
                if (Wl.d <= ZERO || Wl.p <= ZERO) Wl = SolnBlk.W[i][j];
	        dX = SolnBlk.Grid.xfaceW(i+1, j)-SolnBlk.Grid.Cell[i+1][j].Xc;
                Wr = SolnBlk.W[i+1][j] + 
  	             (SolnBlk.phi[i+1][j]^SolnBlk.dWdx[i+1][j])*dX.x +
	             (SolnBlk.phi[i+1][j]^SolnBlk.dWdy[i+1][j])*dX.y;
                if (Wr.d <= ZERO || Wr.p <= ZERO) Wr = SolnBlk.W[i+1][j];
  	     } /* endif */

             switch(Flux_Function_Type) {
               case FLUX_FUNCTION_GODUNOV :
                 Flux = FluxGodunov_n(Wl, Wr, SolnBlk.Grid.nfaceE(i, j));
                 break;
               case FLUX_FUNCTION_ROE :
                 Flux = FluxRoe_n(Wl, Wr, SolnBlk.Grid.nfaceE(i, j));
                 break;
               case FLUX_FUNCTION_RUSANOV :
                 Flux = FluxRusanov_n(Wl, Wr, SolnBlk.Grid.nfaceE(i, j));
                 break;
               case FLUX_FUNCTION_HLLE :
                 Flux = FluxHLLE_n(Wl, Wr, SolnBlk.Grid.nfaceE(i, j));
                 break;
               case FLUX_FUNCTION_LINDE :
                 Flux = FluxLinde_n(Wl, Wr, SolnBlk.Grid.nfaceE(i, j));
                 break;
               case FLUX_FUNCTION_HLLC :
                 Flux = FluxHLLC_n(Wl, Wr, SolnBlk.Grid.nfaceE(i, j));
                 break;
               default:
                 Flux = FluxRoe_n(Wl, Wr, SolnBlk.Grid.nfaceE(i, j));
                 break;
             } /* endswitch */
    
             /* Evaluate cell-averaged solution changes. */
    
             SolnBlk.dUdt[i][j][0] -= 
                Flux*SolnBlk.Grid.lfaceE(i, j)/
                SolnBlk.Grid.Cell[i][j].A;
             SolnBlk.dUdt[i+1][j][0] += 
                Flux*SolnBlk.Grid.lfaceW(i+1, j)/
                SolnBlk.Grid.Cell[i+1][j].A;

          } /* endif */
       } /* endfor */
    
       if ( j > SolnBlk.JCl-1 && j < SolnBlk.JCu+1 ) {
          SolnBlk.dUdt[SolnBlk.ICl-1][j][0] = Ion5Moment2D_U_VACUUM;
          SolnBlk.dUdt[SolnBlk.ICu+1][j][0] = Ion5Moment2D_U_VACUUM;
       } /* endif */
    } /* endfor */
    
    // Add j-direction (eta-direction) fluxes.
    for ( i = SolnBlk.ICl ; i <= SolnBlk.ICu ; ++i ) {
       for ( j  = SolnBlk.JCl-1 ; j <= SolnBlk.JCu ; ++j ) {
    
          /* Evaluate the cell interface j-direction fluxes. */
         
	  if (j == SolnBlk.JCl-1 && 
              (SolnBlk.Grid.BCtypeS[i] == BC_REFLECTION ||
               SolnBlk.Grid.BCtypeS[i] == BC_CHARACTERISTIC)) {
	     dX = SolnBlk.Grid.xfaceS(i, j+1)-SolnBlk.Grid.Cell[i][j+1].Xc;
	     Wr = SolnBlk.W[i][j+1] +
	          (SolnBlk.phi[i][j+1]^SolnBlk.dWdx[i][j+1])*dX.x +
	          (SolnBlk.phi[i][j+1]^SolnBlk.dWdy[i][j+1])*dX.y;
             if (Wr.d <= ZERO || Wr.p <= ZERO) Wr = SolnBlk.W[i][j+1];
             if (SolnBlk.Grid.BCtypeS[i] == BC_REFLECTION) {
               Wl = Reflect(Wr, SolnBlk.Grid.nfaceS(i, j+1));
             } else {
               Wl = BC_Characteristic_Pressure(Wr, 
                                               SolnBlk.WoS[i], 
                                               SolnBlk.Grid.nfaceS(i, j+1));
             } /* endif */
          } else if (j == SolnBlk.JCu && 
                     (SolnBlk.Grid.BCtypeN[i] == BC_REFLECTION ||
                      SolnBlk.Grid.BCtypeN[i] == BC_CHARACTERISTIC)) {
	    dX = SolnBlk.Grid.xfaceN(i, j)-SolnBlk.Grid.Cell[i][j].Xc;
	    Wl = SolnBlk.W[i][j] + 
	      (SolnBlk.phi[i][j]^SolnBlk.dWdx[i][j])*dX.x +
	      (SolnBlk.phi[i][j]^SolnBlk.dWdy[i][j])*dX.y;
            if (Wl.d <= ZERO || Wl.p <= ZERO) Wl = SolnBlk.W[i][j];
	    if (SolnBlk.Grid.BCtypeN[i] == BC_REFLECTION) {
	      Wr = Reflect(Wl, SolnBlk.Grid.nfaceN(i, j));
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
             if (Wl.d <= ZERO || Wl.p <= ZERO) Wl = SolnBlk.W[i][j];
	     dX = SolnBlk.Grid.xfaceS(i, j+1)-SolnBlk.Grid.Cell[i][j+1].Xc;
             Wr = SolnBlk.W[i][j+1] +
                  (SolnBlk.phi[i][j+1]^SolnBlk.dWdx[i][j+1])*dX.x +
                  (SolnBlk.phi[i][j+1]^SolnBlk.dWdy[i][j+1])*dX.y;
             if (Wr.d <= ZERO || Wr.p <= ZERO) Wr = SolnBlk.W[i][j+1];
          } /* endif */

          switch(Flux_Function_Type) {
            case FLUX_FUNCTION_GODUNOV :
              Flux = FluxGodunov_n(Wl, Wr, SolnBlk.Grid.nfaceN(i, j));
              break;
            case FLUX_FUNCTION_ROE :
              Flux = FluxRoe_n(Wl, Wr, SolnBlk.Grid.nfaceN(i, j));
              break;
            case FLUX_FUNCTION_RUSANOV :
              Flux = FluxRusanov_n(Wl, Wr, SolnBlk.Grid.nfaceN(i, j));
              break;
            case FLUX_FUNCTION_HLLE :
              Flux = FluxHLLE_n(Wl, Wr, SolnBlk.Grid.nfaceN(i, j));
              break;
            case FLUX_FUNCTION_LINDE :
              Flux = FluxLinde_n(Wl, Wr, SolnBlk.Grid.nfaceN(i, j));
              break;
            case FLUX_FUNCTION_HLLC :
              Flux = FluxHLLC_n(Wl, Wr, SolnBlk.Grid.nfaceN(i, j));
              break;
            default:
              Flux = FluxRoe_n(Wl, Wr, SolnBlk.Grid.nfaceN(i, j));
              break;
	  } /* endswitch */
    
          /* Evaluate cell-averaged solution changes. */
    
          SolnBlk.dUdt[i][j][0] -= 
             Flux*SolnBlk.Grid.lfaceN(i, j)/
             SolnBlk.Grid.Cell[i][j].A;
          SolnBlk.dUdt[i][j+1][0] += 
             Flux*SolnBlk.Grid.lfaceS(i, j+1)/
             SolnBlk.Grid.Cell[i][j+1].A;

       } /* endfor */
    
       SolnBlk.dUdt[i][SolnBlk.JCl-1][0] = Ion5Moment2D_U_VACUUM;
       SolnBlk.dUdt[i][SolnBlk.JCu+1][0] = Ion5Moment2D_U_VACUUM;
    } /* endfor */
    
    // Output solution changes.
    for ( j  = SolnBlk.JCl ; j <= SolnBlk.JCu ; ++j ) {
       for ( i = SolnBlk.ICl ; i <= SolnBlk.ICu ; ++i ) {

          /* Determine interial terms. */
 	  
          dUdt_inertia = SolnBlk.dUdt[i][j][0];
	  
          /* Determine axisymmetric source terms. */
	  
	  if (SolnBlk.Axisymmetric) {
            dUdt_axi = Sa(SolnBlk.W[i][j], SolnBlk.Grid.Cell[i][j].Xc);
          } else {
            dUdt_axi = Ion5Moment2D_U_VACUUM;
          } /* endif */
	  
          /* Evaluate source terms associated with the 
             electric field. */
	  
          dUdt_efield = Se(SolnBlk.W[i][j], SolnBlk.E[i][j]);
	  
          /* Evaluate source terms associated with the 
             ion-neutral collision processes. */
	  
          dUdt_collisions = Sc(SolnBlk.W[i][j], SolnBlk.Wneut[i][j]);
	  
          /* Output the solution changes. */
	  
          Out_File << " "  << SolnBlk.Grid.Cell[i][j].Xc
                   << dUdt_inertia << dUdt_axi << dUdt_efield << dUdt_collisions
                   << "\n";

       } /* endfor */
    } /* endfor */
    Out_File << setprecision(6);

}
