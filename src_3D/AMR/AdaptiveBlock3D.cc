/* AdaptiveBlock.cc:  Subroutines for adaptive blocks classes
                      that are independent of the CPU architecture. */

/* Include adaptive block header file. */

#ifndef _ADAPTIVBLOCK3D_INCLUDED
#include "AdaptiveBlock3D.h"
#endif // _ADAPTIVEBLOCK3D_INCLUDED

/**********************************************************
 * Routine: Create_Block_Resource_List                    *
 *                                                        *
 * Creates and initializes the adaptive block resource    *
 * list which is used to keep track of available (open)   *
 * solution blocks that can assigned solution values      *
 * during either the solution initialization phase or     *
 * mesh refinement.                                       *
 *                                                        *
 **********************************************************/
void AdaptiveBlock3D_ResourceList::Create_Block_Resource_List(AdaptiveBlock3D_ResourceList &List_of_Available_Blocks,
                                                              const int Number_of_Processors,
                                                              const int Number_of_Blocks_per_Processor) {

    /* Allocate (re-allocate) memory for adaptive block 
       resource list. */

    if (List_of_Available_Blocks.CPU != NULL && 
        List_of_Available_Blocks.Block != NULL) {
       List_of_Available_Blocks.deallocate();
    } /* endif */
    List_of_Available_Blocks.allocate(Number_of_Processors,
                                      Number_of_Blocks_per_Processor);

    /* Assign values to the adaptive block resource list. */

    List_of_Available_Blocks.ThisCPU = CFFC_MPI::This_Processor_Number;
    List_of_Available_Blocks.initialize();   

}

/********************************************************
 * Routine: Broadcast_Adaptive_Block_Info               *
 *                                                      *
 * Broadcast adaptive block information to all          *
 * processors involved in the calculation from the      *
 * primary processor using the MPI broadcast routine.   *
 *                                                      *
 ********************************************************/
void AdaptiveBlock3D_Info::Broadcast_Adaptive_Block_Info(AdaptiveBlock3D_Info &Blk_Info) {

#ifdef _MPI_VERSION
  MPI::COMM_WORLD.Bcast(&(Blk_Info.cpu), 1, MPI::INT, 0);
  MPI::COMM_WORLD.Bcast(&(Blk_Info.blknum), 1, MPI::INT, 0);
  MPI::COMM_WORLD.Bcast(&(Blk_Info.dimen.i), 1, MPI::INT, 0);
  MPI::COMM_WORLD.Bcast(&(Blk_Info.dimen.j), 1, MPI::INT, 0);
  MPI::COMM_WORLD.Bcast(&(Blk_Info.dimen.k), 1, MPI::INT, 0);
  MPI::COMM_WORLD.Bcast(&(Blk_Info.dimen.ghost), 1, MPI::INT, 0);
  MPI::COMM_WORLD.Bcast(&(Blk_Info.sector), 1, MPI::INT, 0);
  MPI::COMM_WORLD.Bcast(&(Blk_Info.level), 1, MPI::INT, 0);
  Blk_Info.blkorient.broadcast();
  Blk_Info.be_on_grid_boundary.broadcast();
#endif

}

/********************************************************
 * Routine: Broadcast_Adaptive_Block                    *
 *                                                      *
 * Broadcast adaptive block information to all          *
 * processors involved in the calculation from the      *
 * primary processor using the MPI broadcast routine.   *
 *                                                      *
 ********************************************************/
void AdaptiveBlock3D::Broadcast_Adaptive_Block(AdaptiveBlock3D &Blk) {
   
#ifdef _MPI_VERSION
   int i;
   MPI::COMM_WORLD.Bcast(&(Blk.used), 1, MPI::INT, 0);
   MPI::COMM_WORLD.Bcast(&(Blk.gblknum), 1, MPI::INT, 0);
   AdaptiveBlock3D_Info::Broadcast_Adaptive_Block_Info(Blk.info);
   MPI::COMM_WORLD.Bcast(&(Blk.nT), 1, MPI::INT, 0);
   for (i = 0; i <= Blk.nT-1; ++i) {
      AdaptiveBlock3D_Info::Broadcast_Adaptive_Block_Info(Blk.infoT[i]);
   } /* endfor */
   MPI::COMM_WORLD.Bcast(&(Blk.nB), 1, MPI::INT, 0);
   for (i = 0; i <= Blk.nB-1; ++i) {
      AdaptiveBlock3D_Info::Broadcast_Adaptive_Block_Info(Blk.infoB[i]);
   } /* endfor */
   MPI::COMM_WORLD.Bcast(&(Blk.nN), 1, MPI::INT, 0);
   for (i = 0; i <= Blk.nN-1; ++i) {
      AdaptiveBlock3D_Info::Broadcast_Adaptive_Block_Info(Blk.infoN[i]);
   } /* endfor */
   MPI::COMM_WORLD.Bcast(&(Blk.nS), 1, MPI::INT, 0);
   for (i = 0; i <= Blk.nS-1; ++i) {
       AdaptiveBlock3D_Info::Broadcast_Adaptive_Block_Info(Blk.infoS[i]);
   } /* endfor */
   MPI::COMM_WORLD.Bcast(&(Blk.nE), 1, MPI::INT, 0);
   for (i = 0; i <= Blk.nE-1; ++i) {
      AdaptiveBlock3D_Info::Broadcast_Adaptive_Block_Info(Blk.infoE[i]);
   } /* endfor */
   MPI::COMM_WORLD.Bcast(&(Blk.nW), 1, MPI::INT, 0);
   for (i = 0; i <= Blk.nW-1; ++i) {
      AdaptiveBlock3D_Info::Broadcast_Adaptive_Block_Info(Blk.infoW[i]);
   } /* endfor */
   MPI::COMM_WORLD.Bcast(&(Blk.nNW), 1, MPI::INT, 0);
   for (i = 0; i <= Blk.nNW-1; ++i) {
      AdaptiveBlock3D_Info::Broadcast_Adaptive_Block_Info(Blk.infoNW[i]);
   } /* endfor */
   MPI::COMM_WORLD.Bcast(&(Blk.nNE), 1, MPI::INT, 0);
   for (i = 0; i <= Blk.nNE-1; ++i) {
      AdaptiveBlock3D_Info::Broadcast_Adaptive_Block_Info(Blk.infoNE[i]);
   } /* endfor */
   MPI::COMM_WORLD.Bcast(&(Blk.nSE), 1, MPI::INT, 0);
   for (i = 0; i <= Blk.nSE-1; ++i) {
      AdaptiveBlock3D_Info::Broadcast_Adaptive_Block_Info(Blk.infoSE[i]);
   } /* endfor */
   MPI::COMM_WORLD.Bcast(&(Blk.nSW), 1, MPI::INT, 0);
   for (i = 0; i <= Blk.nSW-1; ++i) {
      AdaptiveBlock3D_Info::Broadcast_Adaptive_Block_Info(Blk.infoSW[i]);
  }
  MPI::COMM_WORLD.Bcast(&(Blk.nTN), 1, MPI::INT, 0);
  for (i = 0; i <= Blk.nTN-1; ++i) {
     AdaptiveBlock3D_Info::Broadcast_Adaptive_Block_Info(Blk.infoTN[i]);
  } /* endfor */
  MPI::COMM_WORLD.Bcast(&(Blk.nTS), 1, MPI::INT, 0);
  for (i = 0; i <= Blk.nTS-1; ++i) {
     AdaptiveBlock3D_Info::Broadcast_Adaptive_Block_Info(Blk.infoTS[i]);
  } /* endfor */
  MPI::COMM_WORLD.Bcast(&(Blk.nTE), 1, MPI::INT, 0);
  for (i = 0; i <= Blk.nTE-1; ++i) {
     AdaptiveBlock3D_Info::Broadcast_Adaptive_Block_Info(Blk.infoTE[i]);
  } /* endfor */
  MPI::COMM_WORLD.Bcast(&(Blk.nTW), 1, MPI::INT, 0);
  for (i = 0; i <= Blk.nTW-1; ++i) {
     AdaptiveBlock3D_Info::Broadcast_Adaptive_Block_Info(Blk.infoTW[i]);
  }
  MPI::COMM_WORLD.Bcast(&(Blk.nTNW), 1, MPI::INT, 0);
  for (i = 0; i <= Blk.nTNW-1; ++i) {
     AdaptiveBlock3D_Info::Broadcast_Adaptive_Block_Info(Blk.infoTNW[i]);
  } /* endfor */
  MPI::COMM_WORLD.Bcast(&(Blk.nTNE), 1, MPI::INT, 0);
  for (i = 0; i <= Blk.nTNE-1; ++i) {
     AdaptiveBlock3D_Info::Broadcast_Adaptive_Block_Info(Blk.infoTNE[i]);
  } /* endfor */
  MPI::COMM_WORLD.Bcast(&(Blk.nTSE), 1, MPI::INT, 0);
  for (i = 0; i <= Blk.nTSE-1; ++i) {
     AdaptiveBlock3D_Info::Broadcast_Adaptive_Block_Info(Blk.infoTSE[i]);
  } /* endfor */
  MPI::COMM_WORLD.Bcast(&(Blk.nTSW), 1, MPI::INT, 0);
  for (i = 0; i <= Blk.nTSW-1; ++i) {
     AdaptiveBlock3D_Info::Broadcast_Adaptive_Block_Info(Blk.infoTSW[i]);
  }
  MPI::COMM_WORLD.Bcast(&(Blk.nBN), 1, MPI::INT, 0);
  for (i = 0; i <= Blk.nBN-1; ++i) {
     AdaptiveBlock3D_Info::Broadcast_Adaptive_Block_Info(Blk.infoBN[i]);
  } /* endfor */
  MPI::COMM_WORLD.Bcast(&(Blk.nBS), 1, MPI::INT, 0);
  for (i = 0; i <= Blk.nBS-1; ++i) {
     AdaptiveBlock3D_Info::Broadcast_Adaptive_Block_Info(Blk.infoBS[i]);
  } /* endfor */
  MPI::COMM_WORLD.Bcast(&(Blk.nBE), 1, MPI::INT, 0);
  for (i = 0; i <= Blk.nBE-1; ++i) {
     AdaptiveBlock3D_Info::Broadcast_Adaptive_Block_Info(Blk.infoBE[i]);
  } /* endfor */
  MPI::COMM_WORLD.Bcast(&(Blk.nBW), 1, MPI::INT, 0);
  for (i = 0; i <= Blk.nBW-1; ++i) {
     AdaptiveBlock3D_Info::Broadcast_Adaptive_Block_Info(Blk.infoBW[i]);
  }
  MPI::COMM_WORLD.Bcast(&(Blk.nBNW), 1, MPI::INT, 0);
  for (i = 0; i <= Blk.nBNW-1; ++i) {
     AdaptiveBlock3D_Info::Broadcast_Adaptive_Block_Info(Blk.infoBNW[i]);
  } /* endfor */
  MPI::COMM_WORLD.Bcast(&(Blk.nBNE), 1, MPI::INT, 0);
  for (i = 0; i <= Blk.nBNE-1; ++i) {
     AdaptiveBlock3D_Info::Broadcast_Adaptive_Block_Info(Blk.infoBNE[i]);
  } /* endfor */
  MPI::COMM_WORLD.Bcast(&(Blk.nBSE), 1, MPI::INT, 0);
  for (i = 0; i <= Blk.nBSE-1; ++i) {
     AdaptiveBlock3D_Info::Broadcast_Adaptive_Block_Info(Blk.infoBSE[i]);
  } /* endfor */
  MPI::COMM_WORLD.Bcast(&(Blk.nBSW), 1, MPI::INT, 0);
  for (i = 0; i <= Blk.nBSW-1; ++i) {
     AdaptiveBlock3D_Info::Broadcast_Adaptive_Block_Info(Blk.infoBSW[i]);
  } /* endfor */
#endif

}

/**********************************************************
 * Routine: Allocate_Message_Buffers                      *
 *                                                        *
 * Allocates memory for all message passing buffers used  *
 * to send solution information between neighbouring      *
 * adaptive blocks.                                       *
 *                                                        *
 **********************************************************/
void AdaptiveBlock3D_List::Allocate_Message_Buffers(AdaptiveBlock3D_List &Blk_List,
                                                    const int Number_of_Solution_Variables) {

   AdaptiveBlock3D_List::Allocate_Message_Buffers_NoResChange(Blk_List,
                                                              Number_of_Solution_Variables);
   AdaptiveBlock3D_List::Allocate_Message_Buffers_ResChange(Blk_List,
                                                            Number_of_Solution_Variables);

}

/**********************************************************
 * Routine: Allocate_Message_Buffers_NoResChange          *
 *                                                        *
 * Allocates memory for all message passing buffers used  *
 * to send solution information between neighbouring      *
 * adaptive blocks with no mesh resolution changes.       *
 *                                                        *
 **********************************************************/
void AdaptiveBlock3D_List::Allocate_Message_Buffers_NoResChange(AdaptiveBlock3D_List &Blk_List,
                                                                const int Number_of_Solution_Variables) {

   int buffer_size, neighbour_buffer_size;
   int buffer_size_soln, buffer_size_geometry, buffer_size_bcs;
   
   int i_bound_elem; // index for boundary element, face edge or vertex
   int number_neighbours[MAX_BOUNDARY_ELEMENTS_FOR_A_BLOCK];
   AdaptiveBlock3D_Info neighbour_info[MAX_BOUNDARY_ELEMENTS_FOR_A_BLOCK];

   /* Ensure that memory for the block index of the send and receive buffers
      has been allocated. */

   if (Blk_List.Nblk > 0 &&
       Blk_List.message_noreschange_sendbuf == NULL) {
      Blk_List.message_noreschange_sendbuf = new double**[Blk_List.Nblk];
      for (int i_blk = 0; i_blk <= Blk_List.Nblk-1 ; ++i_blk ) {
         Blk_List.message_noreschange_sendbuf[i_blk] = new double*[MAX_BOUNDARY_ELEMENTS_FOR_A_BLOCK];
         for (int j_neigh = 0; j_neigh < MAX_BOUNDARY_ELEMENTS_FOR_A_BLOCK ; ++j_neigh) {
            Blk_List.message_noreschange_sendbuf[i_blk][j_neigh] = new double[1];
            Blk_List.message_noreschange_sendbuf[i_blk][j_neigh][0] = ZERO;
         } /* endfor */
      } /* endfor */
   } /* endif */

   if (Blk_List.Nblk > 0 &&
       Blk_List.message_noreschange_recbuf == NULL) {
      Blk_List.message_noreschange_recbuf = new double**[Blk_List.Nblk];
      for (int i_blk = 0; i_blk <= Blk_List.Nblk-1 ; ++i_blk ) {
         Blk_List.message_noreschange_recbuf[i_blk] = new double*[MAX_BOUNDARY_ELEMENTS_FOR_A_BLOCK];
         for (int j_neigh = 0; j_neigh < MAX_BOUNDARY_ELEMENTS_FOR_A_BLOCK ; ++j_neigh) {
            Blk_List.message_noreschange_recbuf[i_blk][j_neigh] = new double[1];
            Blk_List.message_noreschange_recbuf[i_blk][j_neigh][0] = ZERO;
         } /* endfor */
      } /* endfor */
   } /* endif */
  
   /* For each local adaptive block, create require send and receive message buffers. */

   for (int i_blk = 0; i_blk <= Blk_List.Nblk-1; ++i_blk) {
      if (Blk_List.Block[i_blk].used) { // Consider used blocks only

         /* Assign the boundary element information. */

         number_neighbours[BE::BSW] = Blk_List.Block[i_blk].nBSW;
         neighbour_info[BE::BSW] =  Blk_List.Block[i_blk].infoBSW[0];
      
         number_neighbours[BE::SW] = Blk_List.Block[i_blk].nSW;
         neighbour_info[BE::SW] =  Blk_List.Block[i_blk].infoSW[0];
      
         number_neighbours[BE::TSW] = Blk_List.Block[i_blk].nTSW;
         neighbour_info[BE::TSW] =  Blk_List.Block[i_blk].infoTSW[0];
      
         number_neighbours[BE::BW] = Blk_List.Block[i_blk].nBW;
         neighbour_info[BE::BW] =  Blk_List.Block[i_blk].infoBW[0];
      
         number_neighbours[BE::W] = Blk_List.Block[i_blk].nW;
         neighbour_info[BE::W] =  Blk_List.Block[i_blk].infoW[0];
      
         number_neighbours[BE::TW] = Blk_List.Block[i_blk].nTW;
         neighbour_info[BE::TW] =  Blk_List.Block[i_blk].infoTW[0];
      
         number_neighbours[BE::BNW] = Blk_List.Block[i_blk].nBNW;
         neighbour_info[BE::BNW] =  Blk_List.Block[i_blk].infoBNW[0];
      
         number_neighbours[BE::NW] = Blk_List.Block[i_blk].nNW;
         neighbour_info[BE::NW] =  Blk_List.Block[i_blk].infoNW[0];
         
         number_neighbours[BE::TNW] = Blk_List.Block[i_blk].nTNW;
         neighbour_info[BE::TNW] =  Blk_List.Block[i_blk].infoTNW[0];
      
         number_neighbours[BE::BS] = Blk_List.Block[i_blk].nBS;
         neighbour_info[BE::BS] =  Blk_List.Block[i_blk].infoBS[0];
      
         number_neighbours[BE::S] = Blk_List.Block[i_blk].nS;
         neighbour_info[BE::S] =  Blk_List.Block[i_blk].infoS[0];
      
         number_neighbours[BE::TS] = Blk_List.Block[i_blk].nTS;
         neighbour_info[BE::TS] =  Blk_List.Block[i_blk].infoTS[0];
      
         number_neighbours[BE::B] = Blk_List.Block[i_blk].nB;
         neighbour_info[BE::B] =  Blk_List.Block[i_blk].infoB[0];
      
         number_neighbours[BE::T] = Blk_List.Block[i_blk].nT;
         neighbour_info[BE::T] =  Blk_List.Block[i_blk].infoT[0];
      
         number_neighbours[BE::BN] = Blk_List.Block[i_blk].nBN;
         neighbour_info[BE::BN] =  Blk_List.Block[i_blk].infoBN[0];
      
         number_neighbours[BE::N] = Blk_List.Block[i_blk].nN;
         neighbour_info[BE::N] =  Blk_List.Block[i_blk].infoN[0];
      
         number_neighbours[BE::TN] = Blk_List.Block[i_blk].nTN;
         neighbour_info[BE::TN] =  Blk_List.Block[i_blk].infoTN[0];
      
         number_neighbours[BE::BSE] = Blk_List.Block[i_blk].nBSE;
         neighbour_info[BE::BSE] =  Blk_List.Block[i_blk].infoBSE[0];
      
         number_neighbours[BE::SE] = Blk_List.Block[i_blk].nSE;
         neighbour_info[BE::SE] =  Blk_List.Block[i_blk].infoSE[0];
      
         number_neighbours[BE::TSE] = Blk_List.Block[i_blk].nTSE;
         neighbour_info[BE::TSE] =  Blk_List.Block[i_blk].infoTSE[0];
      
         number_neighbours[BE::BE] = Blk_List.Block[i_blk].nBE;
         neighbour_info[BE::BE] =  Blk_List.Block[i_blk].infoBE[0];
      
         number_neighbours[BE::E] = Blk_List.Block[i_blk].nE;
         neighbour_info[BE::E] =  Blk_List.Block[i_blk].infoE[0];
      
         number_neighbours[BE::TE] = Blk_List.Block[i_blk].nTE;
         neighbour_info[BE::TE] =  Blk_List.Block[i_blk].infoTE[0]; 
      
         number_neighbours[BE::BNE] = Blk_List.Block[i_blk].nBNE;
         neighbour_info[BE::BNE] =  Blk_List.Block[i_blk].infoBNE[0]; 
      
         number_neighbours[BE::NE] = Blk_List.Block[i_blk].nNE;
         neighbour_info[BE::NE] =  Blk_List.Block[i_blk].infoNE[0];  
      
         number_neighbours[BE::TNE] = Blk_List.Block[i_blk].nTNE;
         neighbour_info[BE::TNE] =  Blk_List.Block[i_blk].infoTNE[0]; 
      
         /* Deallocate existing memory for send and receive buffers. */
   
         if (Blk_List.message_noreschange_sendbuf[i_blk] != NULL) {
            for (int j_neigh = 0; j_neigh < MAX_BOUNDARY_ELEMENTS_FOR_A_BLOCK ; ++j_neigh) {
               delete []Blk_List.message_noreschange_sendbuf[i_blk][j_neigh];
               Blk_List.message_noreschange_sendbuf[i_blk][j_neigh] = NULL;
            } /* endfor */
         } /* endif */

         if (Blk_List.message_noreschange_recbuf[i_blk] != NULL) {
            for (int j_neigh = 0; j_neigh < MAX_BOUNDARY_ELEMENTS_FOR_A_BLOCK ; ++j_neigh) {
               delete []Blk_List.message_noreschange_recbuf[i_blk][j_neigh];
               Blk_List.message_noreschange_recbuf[i_blk][j_neigh] = NULL;
            } /* endfor */
         } /* endif */

         /* Reallocate the send and receive message buffers based on the boundary element types. */

         for (int ii = -1; ii<2; ++ii){
            for (int jj = -1; jj<2; ++jj){
               for (int kk = -1; kk<2; ++kk){
                  i_bound_elem = 9*(ii+1) + 3*(jj+1) + (kk+1);

                  if ((number_neighbours[i_bound_elem] == 1) && 
                      (i_bound_elem != BE::ME) &&
                      (Blk_List.Block[i_blk].info.level == neighbour_info[i_bound_elem].level)) {
                     buffer_size_soln = ((abs(ii)*Blk_List.Block[i_blk].info.dimen.ghost) + 
                                        ((!ii)*abs(Blk_List.Block[i_blk].info.dimen.i)))*
                                        ((abs(jj)*Blk_List.Block[i_blk].info.dimen.ghost) + 
                                        ((!jj)*abs(Blk_List.Block[i_blk].info.dimen.j)))*
                                        ((abs(kk)*Blk_List.Block[i_blk].info.dimen.ghost) + 
                                        ((!kk)*abs(Blk_List.Block[i_blk].info.dimen.k)))*
                                        (Number_of_Solution_Variables);
                     buffer_size_geometry = ((abs(ii)*Blk_List.Block[i_blk].info.dimen.ghost)+
                                            ((!ii)*abs(Blk_List.Block[i_blk].info.dimen.i)+1))*
                                            ((abs(jj)*Blk_List.Block[i_blk].info.dimen.ghost) + 
                                            ((!jj)*abs(Blk_List.Block[i_blk].info.dimen.j)+1))*
                                            ((abs(kk)*Blk_List.Block[i_blk].info.dimen.ghost) + 
                                            ((!kk)*abs(Blk_List.Block[i_blk].info.dimen.k)+1))*
                                            (NUM_COMP_VECTOR3D);
                     buffer_size_bcs = (((!ii)*abs(Blk_List.Block[i_blk].info.dimen.i)))*
                                       (((!jj)*abs(Blk_List.Block[i_blk].info.dimen.j)))*
                                       (((!kk)*abs(Blk_List.Block[i_blk].info.dimen.k)));

                     buffer_size = max(buffer_size_soln, buffer_size_geometry + buffer_size_bcs);
                     neighbour_buffer_size = buffer_size;

                     Blk_List.message_noreschange_sendbuf[i_blk][i_bound_elem] = 
                        new double[neighbour_buffer_size];
                     for (int l = 0; l <= neighbour_buffer_size-1; ++l) {
                        Blk_List.message_noreschange_sendbuf[i_blk][i_bound_elem][l] = ZERO;
                     } /* endfor */
                  
                     Blk_List.message_noreschange_recbuf[i_blk][i_bound_elem] = 
                        new double[buffer_size];
                     for (int l = 0; l <= buffer_size-1; ++l) {
                        Blk_List.message_noreschange_recbuf[i_blk][i_bound_elem][l] = ZERO;
                     } /* endfor */

                  } else {
                     Blk_List.message_noreschange_sendbuf[i_blk][i_bound_elem] = new double[1];
                     Blk_List.message_noreschange_sendbuf[i_blk][i_bound_elem][0] = ZERO;
                     Blk_List.message_noreschange_recbuf[i_blk][i_bound_elem] = new double[1];
                     Blk_List.message_noreschange_recbuf[i_blk][i_bound_elem][0] = ZERO;

                  } /* endif */

	       }/* end of for k*/
            }/* end of for j*/
         }/* end of for i*/

      } /* endif */
   }/* endfor */
}

/**********************************************************
 * Routine: Allocate_Message_Buffers_ResChange            *
 *                                                        *
 * Allocates memory for all message passing buffers used  *
 * to send solution information between neighbouring      *
 * adaptive blocks with mesh resolution changes.          *
 *                                                        *
 **********************************************************/
void AdaptiveBlock3D_List::Allocate_Message_Buffers_ResChange(AdaptiveBlock3D_List &Blk_List,
                                                              const int Number_of_Solution_Variables) {
    
}

/**********************************************************
 * Routine: Deallocate_Message_Buffers                    *
 *                                                        *
 * Deallocates memory for all message passing buffers     *
 * used to send solution information between neighbouring *
 * adaptive blocks.                                       *
 *                                                        *
 **********************************************************/
void AdaptiveBlock3D_List::Deallocate_Message_Buffers(AdaptiveBlock3D_List &Blk_List) {

    AdaptiveBlock3D_List::Deallocate_Message_Buffers_NoResChange(Blk_List);
    AdaptiveBlock3D_List::Deallocate_Message_Buffers_ResChange(Blk_List);

}

/**********************************************************
 * Routine: Deallocate_Message_Buffers_NoResChange        *
 *                                                        *
 * Deallocates memory for all message passing buffers     *
 * used to send solution information between neighbouring *
 * adaptive blocks with no mesh resolution changes.       *
 *                                                        *
 **********************************************************/
void AdaptiveBlock3D_List::Deallocate_Message_Buffers_NoResChange(AdaptiveBlock3D_List &Blk_List) {
   
   // Deallocate memory for send and receive buffers.
   if (Blk_List.message_noreschange_sendbuf != NULL) {
      for (int i_blk = 0; i_blk <= Blk_List.Nblk-1; ++i_blk ) {
         for (int j_neigh = 0; j_neigh < MAX_BOUNDARY_ELEMENTS_FOR_A_BLOCK ; ++j_neigh) {
            delete []Blk_List.message_noreschange_sendbuf[i_blk][j_neigh];
            Blk_List.message_noreschange_sendbuf[i_blk][j_neigh] = NULL;
         } /* endfor */
         delete []Blk_List.message_noreschange_sendbuf[i_blk];
         Blk_List.message_noreschange_sendbuf[i_blk] = NULL;
      } /* endfor */
      delete []Blk_List.message_noreschange_sendbuf;
      Blk_List.message_noreschange_sendbuf = NULL;
   } /* endif */

   if (Blk_List.message_noreschange_recbuf != NULL) {
      for (int i_blk = 0; i_blk <= Blk_List.Nblk-1; ++i_blk ) {
         for (int j_neigh = 0; j_neigh < MAX_BOUNDARY_ELEMENTS_FOR_A_BLOCK ; ++j_neigh) {
            delete []Blk_List.message_noreschange_recbuf[i_blk][j_neigh];
            Blk_List.message_noreschange_recbuf[i_blk][j_neigh] = NULL;
         } /* endfor */
         delete []Blk_List.message_noreschange_recbuf[i_blk];
         Blk_List.message_noreschange_recbuf[i_blk] = NULL;
      } /* endfor */
      delete []Blk_List.message_noreschange_recbuf;
      Blk_List.message_noreschange_recbuf = NULL;
   } /* endif */

}

/**********************************************************
 * Routine: Deallocate_Message_Buffers_ResChange          *
 *                                                        *
 * Deallocates memory for all message passing buffers     *
 * used to send solution information between neighbouring *
 * adaptive blocks with mesh resolution changes.          *
 *                                                        *
 **********************************************************/
void AdaptiveBlock3D_List::Deallocate_Message_Buffers_ResChange(AdaptiveBlock3D_List &Blk_List) {

}

/**********************************************************
 * Routine: Copy_Refinement_Flags                         *
 *                                                        *
 * Copies the refinement flags from adaptive block list 1 *
 * to adaptive block list 2.                              *
 *                                                        *
 **********************************************************/
void AdaptiveBlock3D_List::Copy_Refinement_Flags(AdaptiveBlock3D_List &Blk_List_1,
                                                 AdaptiveBlock3D_List &Blk_List_2) {
 
    if (Blk_List_1.Nblk > 0 &&
        Blk_List_1.Nblk == Blk_List_2.Nblk) {
       for (int i_blk = 0; i_blk <= Blk_List_1.Nblk-1 ; ++i_blk ) {
          Blk_List_1.RefineFlag[i_blk] = Blk_List_2.RefineFlag[i_blk];
       } /* endfor */
    } /* endif */

}
