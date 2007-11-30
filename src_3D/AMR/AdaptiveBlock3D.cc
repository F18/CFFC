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
  Blk_Info.be_on_domain_extent.broadcast();
  

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
// not dealing with resolution change now July 16, 2007 xinfeng...
//     Allocate_Message_Buffers_ResChange(Blk_List,
//                                        Number_of_Solution_Variables);

}

/**********************************************************
 * Routine: Allocate_Message_Buffers_NoResChange          *
 *                                                        *
 * Allocates memory for all message passing buffers used  *
 * to send solution information between neighbouring      *
 * adaptive blocks with no mesh resolution changes.       *
 *                                                        *
 **********************************************************/
void  AdaptiveBlock3D_List::Allocate_Message_Buffers_NoResChange(AdaptiveBlock3D_List &Blk_List,
                                          const int Number_of_Solution_Variables) {
   int i_blk, buffer_size, buffer_size_neighbour, l;
   int i_bound_elem; // index for boundary element, face edge or vertex
   int n_bound_elem[MAX_BOUNDARY_ELEMENTS_FOR_A_BLOCK];
   int j_neigh;
   int buffer_size_soln, buffer_size_geometry, buffer_size_bcs;
   
   
   AdaptiveBlock3D_Info info_bound_elem[MAX_BOUNDARY_ELEMENTS_FOR_A_BLOCK];

   /* Ensure that memory for the block index of the send and receive buffers
      has been allocated. */

   if (Blk_List.Nblk > 0 &&
       Blk_List.message_noreschange_sendbuf == NULL) {
      Blk_List.message_noreschange_sendbuf = new double**[Blk_List.Nblk];
      for ( i_blk = 0; i_blk <= Blk_List.Nblk-1 ; ++i_blk ) {
         Blk_List.message_noreschange_sendbuf[i_blk] = new double*[MAX_BOUNDARY_ELEMENTS_FOR_A_BLOCK];
         for (j_neigh = 0; j_neigh < MAX_BOUNDARY_ELEMENTS_FOR_A_BLOCK ; ++j_neigh) {
            Blk_List.message_noreschange_sendbuf[i_blk][j_neigh] = new double[1];
            Blk_List.message_noreschange_sendbuf[i_blk][j_neigh][0] = ZERO;
         } /* endfor */
      } /* endfor */
   } /* endif */
   if (Blk_List.Nblk > 0 &&
       Blk_List.message_noreschange_recbuf == NULL) {
      Blk_List.message_noreschange_recbuf = new double**[Blk_List.Nblk];
      for ( i_blk = 0; i_blk <= Blk_List.Nblk-1 ; ++i_blk ) {
         Blk_List.message_noreschange_recbuf[i_blk] = new double*[MAX_BOUNDARY_ELEMENTS_FOR_A_BLOCK];
         for (j_neigh = 0; j_neigh < MAX_BOUNDARY_ELEMENTS_FOR_A_BLOCK ; ++j_neigh) {
            Blk_List.message_noreschange_recbuf[i_blk][j_neigh] = new double[1];
            Blk_List.message_noreschange_recbuf[i_blk][j_neigh][0] = ZERO;
         } /* endfor */
      } /* endfor */
   } /* endif */
   
  
   for ( i_blk = 0; i_blk <= Blk_List.Nblk-1; ++i_blk ) {
      
     // assign the boudnary element info
      n_bound_elem[BE::BSW] = Blk_List.Block[i_blk].nBSW;
      info_bound_elem[BE::BSW] =  Blk_List.Block[i_blk].infoBSW[0];
      
      n_bound_elem[BE::SW] = Blk_List.Block[i_blk].nSW;
      info_bound_elem[BE::SW] =  Blk_List.Block[i_blk].infoSW[0];
      
      n_bound_elem[BE::TSW] = Blk_List.Block[i_blk].nTSW;
      info_bound_elem[BE::TSW] =  Blk_List.Block[i_blk].infoTSW[0];
      
      n_bound_elem[BE::BW] = Blk_List.Block[i_blk].nBW;
      info_bound_elem[BE::BW] =  Blk_List.Block[i_blk].infoBW[0];
      
      n_bound_elem[BE::W] = Blk_List.Block[i_blk].nW;
      info_bound_elem[BE::W] =  Blk_List.Block[i_blk].infoW[0];
      
      n_bound_elem[BE::TW] = Blk_List.Block[i_blk].nTW;
      info_bound_elem[BE::TW] =  Blk_List.Block[i_blk].infoTW[0];
      
      n_bound_elem[BE::BNW] = Blk_List.Block[i_blk].nBNW;
      info_bound_elem[BE::BNW] =  Blk_List.Block[i_blk].infoBNW[0];
      
      n_bound_elem[BE::NW] = Blk_List.Block[i_blk].nNW;
      info_bound_elem[BE::NW] =  Blk_List.Block[i_blk].infoNW[0];
      
      n_bound_elem[BE::TNW] = Blk_List.Block[i_blk].nTNW;
      info_bound_elem[BE::TNW] =  Blk_List.Block[i_blk].infoTNW[0];
      
      n_bound_elem[BE::BS] = Blk_List.Block[i_blk].nBS;
      info_bound_elem[BE::BS] =  Blk_List.Block[i_blk].infoBS[0];
      
      n_bound_elem[BE::S] = Blk_List.Block[i_blk].nS;
      info_bound_elem[BE::S] =  Blk_List.Block[i_blk].infoS[0];
      
      n_bound_elem[BE::TS] = Blk_List.Block[i_blk].nTS;
      info_bound_elem[BE::TS] =  Blk_List.Block[i_blk].infoTS[0];
      
      n_bound_elem[BE::B] = Blk_List.Block[i_blk].nB;
      info_bound_elem[BE::B] =  Blk_List.Block[i_blk].infoB[0];
      
      n_bound_elem[BE::T] = Blk_List.Block[i_blk].nT;
      info_bound_elem[BE::T] =  Blk_List.Block[i_blk].infoT[0];
      
      n_bound_elem[BE::BN] = Blk_List.Block[i_blk].nBN;
      info_bound_elem[BE::BN] =  Blk_List.Block[i_blk].infoBN[0];
      
      n_bound_elem[BE::N] = Blk_List.Block[i_blk].nN;
      info_bound_elem[BE::N] =  Blk_List.Block[i_blk].infoN[0];
      
      n_bound_elem[BE::TN] = Blk_List.Block[i_blk].nTN;
      info_bound_elem[BE::TN] =  Blk_List.Block[i_blk].infoTN[0];
      
      n_bound_elem[BE::BSE] = Blk_List.Block[i_blk].nBSE;
      info_bound_elem[BE::BSE] =  Blk_List.Block[i_blk].infoBSE[0];
      
      n_bound_elem[BE::SE] = Blk_List.Block[i_blk].nSE;
      info_bound_elem[BE::SE] =  Blk_List.Block[i_blk].infoSE[0];
      
      n_bound_elem[BE::TSE] = Blk_List.Block[i_blk].nTSE;
      info_bound_elem[BE::TSE] =  Blk_List.Block[i_blk].infoTSE[0];
      
      n_bound_elem[BE::BE] = Blk_List.Block[i_blk].nBE;
      info_bound_elem[BE::BE] =  Blk_List.Block[i_blk].infoBE[0];
      
      n_bound_elem[BE::E] = Blk_List.Block[i_blk].nE;
      info_bound_elem[BE::E] =  Blk_List.Block[i_blk].infoE[0];
      
      n_bound_elem[BE::TE] = Blk_List.Block[i_blk].nTE;
      info_bound_elem[BE::TE] =  Blk_List.Block[i_blk].infoTE[0]; 
      
      n_bound_elem[BE::BNE] = Blk_List.Block[i_blk].nBNE;
      info_bound_elem[BE::BNE] =  Blk_List.Block[i_blk].infoBNE[0]; 
      
      n_bound_elem[BE::NE] = Blk_List.Block[i_blk].nNE;
      info_bound_elem[BE::NE] =  Blk_List.Block[i_blk].infoNE[0];  
      
      n_bound_elem[BE::TNE] = Blk_List.Block[i_blk].nTNE;
      info_bound_elem[BE::TNE] =  Blk_List.Block[i_blk].infoTNE[0]; 
      
      // Reallocate memory for send and receive buffers.
      if (Blk_List.message_noreschange_sendbuf[i_blk] != NULL) {
         for (j_neigh = 0; j_neigh < MAX_BOUNDARY_ELEMENTS_FOR_A_BLOCK ; ++j_neigh) {
            delete []Blk_List.message_noreschange_sendbuf[i_blk][j_neigh];
            Blk_List.message_noreschange_sendbuf[i_blk][j_neigh] = NULL;
         } /* endfor */
      } /* endif */
      if (Blk_List.message_noreschange_recbuf[i_blk] != NULL) {
         for (j_neigh = 0; j_neigh < MAX_BOUNDARY_ELEMENTS_FOR_A_BLOCK ; ++j_neigh) {
            delete []Blk_List.message_noreschange_recbuf[i_blk][j_neigh];
            Blk_List.message_noreschange_recbuf[i_blk][j_neigh] = NULL;
         } /* endfor */
      } /* endif */
      
      // allocate the buffer based on the element types.
      for (int ii = -1; ii<2; ++ii){
         for (int jj = -1; jj<2; ++jj){
            for (int kk = -1; kk<2; ++kk){
               i_bound_elem = 9*(ii+1) + 3*(jj+1) + (kk+1);

                  if (Blk_List.Block[i_blk].used && 
                   (n_bound_elem[i_bound_elem] == 1) && (i_bound_elem != 13) &&
                   (Blk_List.Block[i_blk].info.level ==  info_bound_elem[i_bound_elem].level)) {
              
                     buffer_size_soln = ((abs(ii)*Blk_List.Block[i_blk].info.dimen.ghost) + ((!ii)*abs(Blk_List.Block[i_blk].info.dimen.i)))*
                        ((abs(jj)*Blk_List.Block[i_blk].info.dimen.ghost) + ((!jj)*abs(Blk_List.Block[i_blk].info.dimen.j)))*
                        ((abs(kk)*Blk_List.Block[i_blk].info.dimen.ghost) + ((!kk)*abs(Blk_List.Block[i_blk].info.dimen.k)))*(Number_of_Solution_Variables) ;
                  
                     buffer_size_geometry =   ((abs(ii)*Blk_List.Block[i_blk].info.dimen.ghost)+ ((!ii)*abs(Blk_List.Block[i_blk].info.dimen.i)+1))*
                        ((abs(jj)*Blk_List.Block[i_blk].info.dimen.ghost) + ((!jj)*abs(Blk_List.Block[i_blk].info.dimen.j)+1))*
                        ((abs(kk)*Blk_List.Block[i_blk].info.dimen.ghost) + ((!kk)*abs(Blk_List.Block[i_blk].info.dimen.k)+1))*(NUM_COMP_VECTOR3D);
                       
                     buffer_size_bcs =  ((abs(ii)*Blk_List.Block[i_blk].info.dimen.ghost) + ((!ii)*abs(Blk_List.Block[i_blk].info.dimen.i)+1))*
                           ((abs(jj)*Blk_List.Block[i_blk].info.dimen.ghost) + ((!jj)*abs(Blk_List.Block[i_blk].info.dimen.j) +1))*
                           ((abs(kk)*Blk_List.Block[i_blk].info.dimen.ghost) + ((!kk)*abs(Blk_List.Block[i_blk].info.dimen.k)+1))*(NUM_COMP_VECTOR3D)+
                           (((!ii)*abs(Blk_List.Block[i_blk].info.dimen.i)))*(((!jj)*abs(Blk_List.Block[i_blk].info.dimen.j)))*
                           (((!kk)*abs(Blk_List.Block[i_blk].info.dimen.k)));
                        
                     
                     buffer_size = max(max(buffer_size_soln, buffer_size_geometry), buffer_size_bcs);
                     

                     buffer_size_neighbour = buffer_size;
                     
                     Blk_List.message_noreschange_sendbuf[i_blk][i_bound_elem] = new double[buffer_size_neighbour];
                  
                     
                     for (l = 0; l <= buffer_size_neighbour-1; ++l)
                        Blk_List.message_noreschange_sendbuf[i_blk][i_bound_elem][l] = ZERO;
                     
                     Blk_List.message_noreschange_recbuf[i_blk][i_bound_elem] = new double[buffer_size];
                     for (l = 0; l <= buffer_size-1; ++l) 
                        Blk_List.message_noreschange_recbuf[i_blk][i_bound_elem][l] = ZERO;
                  } else {
                     Blk_List.message_noreschange_sendbuf[i_blk][i_bound_elem] = 
                        new double[1];
                     Blk_List.message_noreschange_sendbuf[i_blk][i_bound_elem][0] = ZERO;
                     Blk_List.message_noreschange_recbuf[i_blk][i_bound_elem] = 
                        new double[1];
                  Blk_List.message_noreschange_recbuf[i_blk][i_bound_elem][0] = ZERO;
                  } /* endif */
                  
            }/* end of for k*/
         }/* end of for j*/
      }/* end of for i*/
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
void  AdaptiveBlock3D_List::Allocate_Message_Buffers_ResChange(AdaptiveBlock3D_List &Blk_List,
                                        const int Number_of_Solution_Variables) {
    
    int i_blk, j_neigh, buffer_size, buffer_size_neighbour, l;

    return;

}

/**********************************************************
 * Routine: Deallocate_Message_Buffers                    *
 *                                                        *
 * Deallocates memory for all message passing buffers     *
 * used to send solution information between neighbouring *
 * adaptive blocks.                                       *
 *                                                        *
 **********************************************************/
void  AdaptiveBlock3D_List::Deallocate_Message_Buffers(AdaptiveBlock3D_List &Blk_List) {

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
void  AdaptiveBlock3D_List::Deallocate_Message_Buffers_NoResChange(AdaptiveBlock3D_List &Blk_List) {
   
   int i_blk, j_neigh;
   
   // Deallocate memory for send and receive buffers.
   if (Blk_List.message_noreschange_sendbuf != NULL) {
      for ( i_blk = 0; i_blk <= Blk_List.Nblk-1; ++i_blk ) {
         for (j_neigh = 0; j_neigh < MAX_BOUNDARY_ELEMENTS_FOR_A_BLOCK ; ++j_neigh) {
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
      for ( i_blk = 0; i_blk <= Blk_List.Nblk-1; ++i_blk ) {
         for (j_neigh = 0; j_neigh < MAX_BOUNDARY_ELEMENTS_FOR_A_BLOCK ; ++j_neigh) {
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
void  AdaptiveBlock3D_List::Deallocate_Message_Buffers_ResChange(AdaptiveBlock3D_List &Blk_List) {
   
   int i_blk, j_neigh;
   
   
   return;
   

}

/**********************************************************
 * Routine: Exchange_Messages                             *
 *                                                        *
 * Sends and receives solution information contained in   *
 * the preloaded message passing buffers between          *
 * neighbouring adaptive blocks.                          *
 *                                                        *
 **********************************************************/
int  AdaptiveBlock3D_List::Exchange_Messages(AdaptiveBlock3D_List &Blk_List,
                                             const int Number_of_Solution_Variables,
                                             const int Send_Mesh_Geometry_Only) {

    int error_flag;

    /* Exchange message buffers at block interfaces 
       with no cell resolution change. */

    error_flag = Exchange_Messages_NoResChange(Blk_List,
                                               Number_of_Solution_Variables,
                                               Send_Mesh_Geometry_Only);
    if (error_flag) return(error_flag);

    /* Exchange message buffers at block interfaces, 
       sending messages from fine to coarse blocks. */

    error_flag = Exchange_Messages_ResChange_FineToCoarse(Blk_List,
                                                          Number_of_Solution_Variables);
    if (error_flag) return(error_flag);

    /* Exchange message buffers at block interfaces, 
       sending messages from coarse to fine blocks. */

    error_flag = Exchange_Messages_ResChange_CoarseToFine(Blk_List,
                                                          Number_of_Solution_Variables);
    return(error_flag);

}


/**********************************************************
 * Routine: Copy_Refinement_Flags                         *
 *                                                        *
 * Copies the refinement flags from adaptive block list 1 *
 * to adaptive block list 2.                              *
 *                                                        *
 **********************************************************/
void  AdaptiveBlock3D_List::Copy_Refinement_Flags(AdaptiveBlock3D_List &Blk_List_1,
                           AdaptiveBlock3D_List &Blk_List_2) {
 
    int i_blk;

    if (Blk_List_1.Nblk > 0 &&
        Blk_List_1.Nblk == Blk_List_2.Nblk) {
       for ( i_blk = 0; i_blk <= Blk_List_1.Nblk-1 ; ++i_blk ) {
          Blk_List_1.RefineFlag[i_blk] = Blk_List_2.RefineFlag[i_blk];
       } /* endfor */
    } /* endif */

}
