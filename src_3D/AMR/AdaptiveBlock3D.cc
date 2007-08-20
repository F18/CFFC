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

    /* Ensure that memory for the block index of the send and receive buffers
       has been allocated. */

    if (Blk_List.Nblk > 0 &&
        Blk_List.message_noreschange_topface_sendbuf == NULL) {
       Blk_List.message_noreschange_topface_sendbuf = new double*[Blk_List.Nblk];
       for ( i_blk = 0; i_blk <= Blk_List.Nblk-1 ; ++i_blk ) {
          Blk_List.message_noreschange_topface_sendbuf[i_blk] = new double[1];
          Blk_List.message_noreschange_topface_sendbuf[i_blk][0] = ZERO;
       } /* endfor */
    } /* endif */
    if (Blk_List.Nblk > 0 &&
        Blk_List.message_noreschange_topface_recbuf == NULL) {
       Blk_List.message_noreschange_topface_recbuf = new double*[Blk_List.Nblk];
       for ( i_blk = 0; i_blk <= Blk_List.Nblk-1 ; ++i_blk ) {
          Blk_List.message_noreschange_topface_recbuf[i_blk] = new double[1];
          Blk_List.message_noreschange_topface_recbuf[i_blk][0] = ZERO;
       } /* endfor */
    } /* endif */

   if (Blk_List.Nblk > 0 &&
        Blk_List.message_noreschange_bottomface_sendbuf == NULL) {
       Blk_List.message_noreschange_bottomface_sendbuf = new double*[Blk_List.Nblk];
       for ( i_blk = 0; i_blk <= Blk_List.Nblk-1 ; ++i_blk ) {
          Blk_List.message_noreschange_bottomface_sendbuf[i_blk] = new double[1];
          Blk_List.message_noreschange_bottomface_sendbuf[i_blk][0] = ZERO;
       } /* endfor */
    } /* endif */
    if (Blk_List.Nblk > 0 &&
        Blk_List.message_noreschange_bottomface_recbuf == NULL) {
       Blk_List.message_noreschange_bottomface_recbuf = new double*[Blk_List.Nblk];
       for ( i_blk = 0; i_blk <= Blk_List.Nblk-1 ; ++i_blk ) {
          Blk_List.message_noreschange_bottomface_recbuf[i_blk] = new double[1];
          Blk_List.message_noreschange_bottomface_recbuf[i_blk][0] = ZERO;
       } /* endfor */
    } /* endif */

    if (Blk_List.Nblk > 0 &&
        Blk_List.message_noreschange_northface_sendbuf == NULL) {
       Blk_List.message_noreschange_northface_sendbuf = new double*[Blk_List.Nblk];
       for ( i_blk = 0; i_blk <= Blk_List.Nblk-1 ; ++i_blk ) {
          Blk_List.message_noreschange_northface_sendbuf[i_blk] = new double[1];
          Blk_List.message_noreschange_northface_sendbuf[i_blk][0] = ZERO;
       } /* endfor */
    } /* endif */
    if (Blk_List.Nblk > 0 &&
        Blk_List.message_noreschange_northface_recbuf == NULL) {
       Blk_List.message_noreschange_northface_recbuf = new double*[Blk_List.Nblk];
       for ( i_blk = 0; i_blk <= Blk_List.Nblk-1 ; ++i_blk ) {
          Blk_List.message_noreschange_northface_recbuf[i_blk] = new double[1];
          Blk_List.message_noreschange_northface_recbuf[i_blk][0] = ZERO;
       } /* endfor */
    } /* endif */
    if (Blk_List.Nblk > 0 &&
        Blk_List.message_noreschange_southface_sendbuf == NULL) {
       Blk_List.message_noreschange_southface_sendbuf = new double*[Blk_List.Nblk];
       for ( i_blk = 0; i_blk <= Blk_List.Nblk-1 ; ++i_blk ) {
          Blk_List.message_noreschange_southface_sendbuf[i_blk] = new double[1];
          Blk_List.message_noreschange_southface_sendbuf[i_blk][0] = ZERO;
       } /* endfor */
    } /* endif */
    if (Blk_List.Nblk > 0 &&
        Blk_List.message_noreschange_southface_recbuf == NULL) {
       Blk_List.message_noreschange_southface_recbuf = new double*[Blk_List.Nblk];
       for ( i_blk = 0; i_blk <= Blk_List.Nblk-1 ; ++i_blk ) {
          Blk_List.message_noreschange_southface_recbuf[i_blk] = new double[1];
          Blk_List.message_noreschange_southface_recbuf[i_blk][0] = ZERO;
       } /* endfor */
    } /* endif */ 
    if (Blk_List.Nblk > 0 &&
        Blk_List.message_noreschange_eastface_sendbuf == NULL) {
       Blk_List.message_noreschange_eastface_sendbuf = new double*[Blk_List.Nblk];
       for ( i_blk = 0; i_blk <= Blk_List.Nblk-1 ; ++i_blk ) {
          Blk_List.message_noreschange_eastface_sendbuf[i_blk] = new double[1];
          Blk_List.message_noreschange_eastface_sendbuf[i_blk][0] = ZERO;
       } /* endfor */
    } /* endif */
    if (Blk_List.Nblk > 0 &&
        Blk_List.message_noreschange_eastface_recbuf == NULL) {
       Blk_List.message_noreschange_eastface_recbuf = new double*[Blk_List.Nblk];
       for ( i_blk = 0; i_blk <= Blk_List.Nblk-1 ; ++i_blk ) {
          Blk_List.message_noreschange_eastface_recbuf[i_blk] = new double[1];
          Blk_List.message_noreschange_eastface_recbuf[i_blk][0] = ZERO;
       } /* endfor */
    } /* endif */
    if (Blk_List.Nblk > 0 &&
        Blk_List.message_noreschange_westface_sendbuf == NULL) {
       Blk_List.message_noreschange_westface_sendbuf = new double*[Blk_List.Nblk];
       for ( i_blk = 0; i_blk <= Blk_List.Nblk-1 ; ++i_blk ) {
          Blk_List.message_noreschange_westface_sendbuf[i_blk] = new double[1];
          Blk_List.message_noreschange_westface_sendbuf[i_blk][0] = ZERO;
       } /* endfor */
    } /* endif */
    if (Blk_List.Nblk > 0 &&
        Blk_List.message_noreschange_westface_recbuf == NULL) {
       Blk_List.message_noreschange_westface_recbuf = new double*[Blk_List.Nblk];
       for ( i_blk = 0; i_blk <= Blk_List.Nblk-1 ; ++i_blk ) {
          Blk_List.message_noreschange_westface_recbuf[i_blk] = new double[1];
          Blk_List.message_noreschange_westface_recbuf[i_blk][0] = ZERO;
       } /* endfor */
    } /* endif */
    if (Blk_List.Nblk > 0 &&
        Blk_List.message_noreschange_northwestcorner_sendbuf == NULL) {
       Blk_List.message_noreschange_northwestcorner_sendbuf = new double*[Blk_List.Nblk];
       for ( i_blk = 0; i_blk <= Blk_List.Nblk-1 ; ++i_blk ) {
          Blk_List.message_noreschange_northwestcorner_sendbuf[i_blk] = new double[1];
          Blk_List.message_noreschange_northwestcorner_sendbuf[i_blk][0] = ZERO;
       } /* endfor */
    } /* endif */
    if (Blk_List.Nblk > 0 &&
        Blk_List.message_noreschange_northwestcorner_recbuf == NULL) {
       Blk_List.message_noreschange_northwestcorner_recbuf = new double*[Blk_List.Nblk];
       for ( i_blk = 0; i_blk <= Blk_List.Nblk-1 ; ++i_blk ) {
          Blk_List.message_noreschange_northwestcorner_recbuf[i_blk] = new double[1];
          Blk_List.message_noreschange_northwestcorner_recbuf[i_blk][0] = ZERO;
       } /* endfor */
    } /* endif */
    if (Blk_List.Nblk > 0 &&
        Blk_List.message_noreschange_northeastcorner_sendbuf == NULL) {
       Blk_List.message_noreschange_northeastcorner_sendbuf = new double*[Blk_List.Nblk];
       for ( i_blk = 0; i_blk <= Blk_List.Nblk-1 ; ++i_blk ) {
          Blk_List.message_noreschange_northeastcorner_sendbuf[i_blk] = new double[1];
          Blk_List.message_noreschange_northeastcorner_sendbuf[i_blk][0] = ZERO;
       } /* endfor */
    } /* endif */
    if (Blk_List.Nblk > 0 &&
        Blk_List.message_noreschange_northeastcorner_recbuf == NULL) {
       Blk_List.message_noreschange_northeastcorner_recbuf = new double*[Blk_List.Nblk];
       for ( i_blk = 0; i_blk <= Blk_List.Nblk-1 ; ++i_blk ) {
          Blk_List.message_noreschange_northeastcorner_recbuf[i_blk] = new double[1];
          Blk_List.message_noreschange_northeastcorner_recbuf[i_blk][0] = ZERO;
       } /* endfor */
    } /* endif */
    if (Blk_List.Nblk > 0 &&
        Blk_List.message_noreschange_southeastcorner_sendbuf == NULL) {
       Blk_List.message_noreschange_southeastcorner_sendbuf = new double*[Blk_List.Nblk];
       for ( i_blk = 0; i_blk <= Blk_List.Nblk-1 ; ++i_blk ) {
          Blk_List.message_noreschange_southeastcorner_sendbuf[i_blk] = new double[1];
          Blk_List.message_noreschange_southeastcorner_sendbuf[i_blk][0] = ZERO;
       } /* endfor */
    } /* endif */
    if (Blk_List.Nblk > 0 &&
        Blk_List.message_noreschange_southeastcorner_recbuf == NULL) {
       Blk_List.message_noreschange_southeastcorner_recbuf = new double*[Blk_List.Nblk];
       for ( i_blk = 0; i_blk <= Blk_List.Nblk-1 ; ++i_blk ) {
          Blk_List.message_noreschange_southeastcorner_recbuf[i_blk] = new double[1];
          Blk_List.message_noreschange_southeastcorner_recbuf[i_blk][0] = ZERO;
       } /* endfor */
    } /* endif */
    if (Blk_List.Nblk > 0 &&
        Blk_List.message_noreschange_southwestcorner_sendbuf == NULL) {
       Blk_List.message_noreschange_southwestcorner_sendbuf = new double*[Blk_List.Nblk];
       for ( i_blk = 0; i_blk <= Blk_List.Nblk-1 ; ++i_blk ) {
          Blk_List.message_noreschange_southwestcorner_sendbuf[i_blk] = new double[1];
          Blk_List.message_noreschange_southwestcorner_sendbuf[i_blk][0] = ZERO;
       } /* endfor */
    } /* endif */
    if (Blk_List.Nblk > 0 &&
        Blk_List.message_noreschange_southwestcorner_recbuf == NULL) {
       Blk_List.message_noreschange_southwestcorner_recbuf = new double*[Blk_List.Nblk];
       for ( i_blk = 0; i_blk <= Blk_List.Nblk-1 ; ++i_blk ) {
          Blk_List.message_noreschange_southwestcorner_recbuf[i_blk] = new double[1];
          Blk_List.message_noreschange_southwestcorner_recbuf[i_blk][0] = ZERO;
       } /* endfor */
    } /* endif */

    //****//

    if (Blk_List.Nblk > 0 &&
        Blk_List.message_noreschange_topnorthcorner_sendbuf == NULL) {
       Blk_List.message_noreschange_topnorthcorner_sendbuf = new double*[Blk_List.Nblk];
       for ( i_blk = 0; i_blk <= Blk_List.Nblk-1 ; ++i_blk ) {
          Blk_List.message_noreschange_topnorthcorner_sendbuf[i_blk] = new double[1];
          Blk_List.message_noreschange_topnorthcorner_sendbuf[i_blk][0] = ZERO;
       } /* endfor */
    } /* endif */
    if (Blk_List.Nblk > 0 &&
        Blk_List.message_noreschange_topnorthcorner_recbuf == NULL) {
       Blk_List.message_noreschange_topnorthcorner_recbuf = new double*[Blk_List.Nblk];
       for ( i_blk = 0; i_blk <= Blk_List.Nblk-1 ; ++i_blk ) {
          Blk_List.message_noreschange_topnorthcorner_recbuf[i_blk] = new double[1];
          Blk_List.message_noreschange_topnorthcorner_recbuf[i_blk][0] = ZERO;
       } /* endfor */
    } /* endif */
    if (Blk_List.Nblk > 0 &&
        Blk_List.message_noreschange_topsouthcorner_sendbuf == NULL) {
       Blk_List.message_noreschange_topsouthcorner_sendbuf = new double*[Blk_List.Nblk];
       for ( i_blk = 0; i_blk <= Blk_List.Nblk-1 ; ++i_blk ) {
          Blk_List.message_noreschange_topsouthcorner_sendbuf[i_blk] = new double[1];
          Blk_List.message_noreschange_topsouthcorner_sendbuf[i_blk][0] = ZERO;
       } /* endfor */
    } /* endif */
    if (Blk_List.Nblk > 0 &&
        Blk_List.message_noreschange_topsouthcorner_recbuf == NULL) {
       Blk_List.message_noreschange_topsouthcorner_recbuf = new double*[Blk_List.Nblk];
       for ( i_blk = 0; i_blk <= Blk_List.Nblk-1 ; ++i_blk ) {
          Blk_List.message_noreschange_topsouthcorner_recbuf[i_blk] = new double[1];
          Blk_List.message_noreschange_topsouthcorner_recbuf[i_blk][0] = ZERO;
       } /* endfor */
    } /* endif */
    if (Blk_List.Nblk > 0 &&
        Blk_List.message_noreschange_topeastcorner_sendbuf == NULL) {
       Blk_List.message_noreschange_topeastcorner_sendbuf = new double*[Blk_List.Nblk];
       for ( i_blk = 0; i_blk <= Blk_List.Nblk-1 ; ++i_blk ) {
          Blk_List.message_noreschange_topeastcorner_sendbuf[i_blk] = new double[1];
          Blk_List.message_noreschange_topeastcorner_sendbuf[i_blk][0] = ZERO;
       } /* endfor */
    } /* endif */
    if (Blk_List.Nblk > 0 &&
        Blk_List.message_noreschange_topeastcorner_recbuf == NULL) {
       Blk_List.message_noreschange_topeastcorner_recbuf = new double*[Blk_List.Nblk];
       for ( i_blk = 0; i_blk <= Blk_List.Nblk-1 ; ++i_blk ) {
          Blk_List.message_noreschange_topeastcorner_recbuf[i_blk] = new double[1];
          Blk_List.message_noreschange_topeastcorner_recbuf[i_blk][0] = ZERO;
       } /* endfor */
    } /* endif */
    if (Blk_List.Nblk > 0 &&
        Blk_List.message_noreschange_topwestcorner_sendbuf == NULL) {
       Blk_List.message_noreschange_topwestcorner_sendbuf = new double*[Blk_List.Nblk];
       for ( i_blk = 0; i_blk <= Blk_List.Nblk-1 ; ++i_blk ) {
          Blk_List.message_noreschange_topwestcorner_sendbuf[i_blk] = new double[1];
          Blk_List.message_noreschange_topwestcorner_sendbuf[i_blk][0] = ZERO;
       } /* endfor */
    } /* endif */
    if (Blk_List.Nblk > 0 &&
        Blk_List.message_noreschange_topwestcorner_recbuf == NULL) {
       Blk_List.message_noreschange_topwestcorner_recbuf = new double*[Blk_List.Nblk];
       for ( i_blk = 0; i_blk <= Blk_List.Nblk-1 ; ++i_blk ) {
          Blk_List.message_noreschange_topwestcorner_recbuf[i_blk] = new double[1];
          Blk_List.message_noreschange_topwestcorner_recbuf[i_blk][0] = ZERO;
       } /* endfor */
    } /* endif */

    if (Blk_List.Nblk > 0 &&
        Blk_List.message_noreschange_topnorthwestcorner_sendbuf == NULL) {
       Blk_List.message_noreschange_topnorthwestcorner_sendbuf = new double*[Blk_List.Nblk];
       for ( i_blk = 0; i_blk <= Blk_List.Nblk-1 ; ++i_blk ) {
          Blk_List.message_noreschange_topnorthwestcorner_sendbuf[i_blk] = new double[1];
          Blk_List.message_noreschange_topnorthwestcorner_sendbuf[i_blk][0] = ZERO;
       } /* endfor */
    } /* endif */
    if (Blk_List.Nblk > 0 &&
        Blk_List.message_noreschange_topnorthwestcorner_recbuf == NULL) {
       Blk_List.message_noreschange_topnorthwestcorner_recbuf = new double*[Blk_List.Nblk];
       for ( i_blk = 0; i_blk <= Blk_List.Nblk-1 ; ++i_blk ) {
          Blk_List.message_noreschange_topnorthwestcorner_recbuf[i_blk] = new double[1];
          Blk_List.message_noreschange_topnorthwestcorner_recbuf[i_blk][0] = ZERO;
       } /* endfor */
    } /* endif */
    if (Blk_List.Nblk > 0 &&
        Blk_List.message_noreschange_topnortheastcorner_sendbuf == NULL) {
       Blk_List.message_noreschange_topnortheastcorner_sendbuf = new double*[Blk_List.Nblk];
       for ( i_blk = 0; i_blk <= Blk_List.Nblk-1 ; ++i_blk ) {
          Blk_List.message_noreschange_topnortheastcorner_sendbuf[i_blk] = new double[1];
          Blk_List.message_noreschange_topnortheastcorner_sendbuf[i_blk][0] = ZERO;
       } /* endfor */
    } /* endif */
    if (Blk_List.Nblk > 0 &&
        Blk_List.message_noreschange_topnortheastcorner_recbuf == NULL) {
       Blk_List.message_noreschange_topnortheastcorner_recbuf = new double*[Blk_List.Nblk];
       for ( i_blk = 0; i_blk <= Blk_List.Nblk-1 ; ++i_blk ) {
          Blk_List.message_noreschange_topnortheastcorner_recbuf[i_blk] = new double[1];
          Blk_List.message_noreschange_topnortheastcorner_recbuf[i_blk][0] = ZERO;
       } /* endfor */
    } /* endif */
    if (Blk_List.Nblk > 0 &&
        Blk_List.message_noreschange_topsoutheastcorner_sendbuf == NULL) {
       Blk_List.message_noreschange_topsoutheastcorner_sendbuf = new double*[Blk_List.Nblk];
       for ( i_blk = 0; i_blk <= Blk_List.Nblk-1 ; ++i_blk ) {
          Blk_List.message_noreschange_topsoutheastcorner_sendbuf[i_blk] = new double[1];
          Blk_List.message_noreschange_topsoutheastcorner_sendbuf[i_blk][0] = ZERO;
       } /* endfor */
    } /* endif */
    if (Blk_List.Nblk > 0 &&
        Blk_List.message_noreschange_topsoutheastcorner_recbuf == NULL) {
       Blk_List.message_noreschange_topsoutheastcorner_recbuf = new double*[Blk_List.Nblk];
       for ( i_blk = 0; i_blk <= Blk_List.Nblk-1 ; ++i_blk ) {
          Blk_List.message_noreschange_topsoutheastcorner_recbuf[i_blk] = new double[1];
          Blk_List.message_noreschange_topsoutheastcorner_recbuf[i_blk][0] = ZERO;
       } /* endfor */
    } /* endif */
    if (Blk_List.Nblk > 0 &&
        Blk_List.message_noreschange_topsouthwestcorner_sendbuf == NULL) {
       Blk_List.message_noreschange_topsouthwestcorner_sendbuf = new double*[Blk_List.Nblk];
       for ( i_blk = 0; i_blk <= Blk_List.Nblk-1 ; ++i_blk ) {
          Blk_List.message_noreschange_topsouthwestcorner_sendbuf[i_blk] = new double[1];
          Blk_List.message_noreschange_topsouthwestcorner_sendbuf[i_blk][0] = ZERO;
       } /* endfor */
    } /* endif */
    if (Blk_List.Nblk > 0 &&
        Blk_List.message_noreschange_topsouthwestcorner_recbuf == NULL) {
       Blk_List.message_noreschange_topsouthwestcorner_recbuf = new double*[Blk_List.Nblk];
       for ( i_blk = 0; i_blk <= Blk_List.Nblk-1 ; ++i_blk ) {
          Blk_List.message_noreschange_topsouthwestcorner_recbuf[i_blk] = new double[1];
          Blk_List.message_noreschange_topsouthwestcorner_recbuf[i_blk][0] = ZERO;
       } /* endfor */
    } /* endif */

    /*****/

    if (Blk_List.Nblk > 0 &&
        Blk_List.message_noreschange_bottomnorthcorner_sendbuf == NULL) {
       Blk_List.message_noreschange_bottomnorthcorner_sendbuf = new double*[Blk_List.Nblk];
       for ( i_blk = 0; i_blk <= Blk_List.Nblk-1 ; ++i_blk ) {
          Blk_List.message_noreschange_bottomnorthcorner_sendbuf[i_blk] = new double[1];
          Blk_List.message_noreschange_bottomnorthcorner_sendbuf[i_blk][0] = ZERO;
       } /* endfor */
    } /* endif */
    if (Blk_List.Nblk > 0 &&
        Blk_List.message_noreschange_bottomnorthcorner_recbuf == NULL) {
       Blk_List.message_noreschange_bottomnorthcorner_recbuf = new double*[Blk_List.Nblk];
       for ( i_blk = 0; i_blk <= Blk_List.Nblk-1 ; ++i_blk ) {
          Blk_List.message_noreschange_bottomnorthcorner_recbuf[i_blk] = new double[1];
          Blk_List.message_noreschange_bottomnorthcorner_recbuf[i_blk][0] = ZERO;
       } /* endfor */
    } /* endif */
    if (Blk_List.Nblk > 0 &&
        Blk_List.message_noreschange_bottomsouthcorner_sendbuf == NULL) {
       Blk_List.message_noreschange_bottomsouthcorner_sendbuf = new double*[Blk_List.Nblk];
       for ( i_blk = 0; i_blk <= Blk_List.Nblk-1 ; ++i_blk ) {
          Blk_List.message_noreschange_bottomsouthcorner_sendbuf[i_blk] = new double[1];
          Blk_List.message_noreschange_bottomsouthcorner_sendbuf[i_blk][0] = ZERO;
       } /* endfor */
    } /* endif */
    if (Blk_List.Nblk > 0 &&
        Blk_List.message_noreschange_bottomsouthcorner_recbuf == NULL) {
       Blk_List.message_noreschange_bottomsouthcorner_recbuf = new double*[Blk_List.Nblk];
       for ( i_blk = 0; i_blk <= Blk_List.Nblk-1 ; ++i_blk ) {
          Blk_List.message_noreschange_bottomsouthcorner_recbuf[i_blk] = new double[1];
          Blk_List.message_noreschange_bottomsouthcorner_recbuf[i_blk][0] = ZERO;
       } /* endfor */
    } /* endif */
    if (Blk_List.Nblk > 0 &&
        Blk_List.message_noreschange_bottomeastcorner_sendbuf == NULL) {
       Blk_List.message_noreschange_bottomeastcorner_sendbuf = new double*[Blk_List.Nblk];
       for ( i_blk = 0; i_blk <= Blk_List.Nblk-1 ; ++i_blk ) {
          Blk_List.message_noreschange_bottomeastcorner_sendbuf[i_blk] = new double[1];
          Blk_List.message_noreschange_bottomeastcorner_sendbuf[i_blk][0] = ZERO;
       } /* endfor */
    } /* endif */
    if (Blk_List.Nblk > 0 &&
        Blk_List.message_noreschange_bottomeastcorner_recbuf == NULL) {
       Blk_List.message_noreschange_bottomeastcorner_recbuf = new double*[Blk_List.Nblk];
       for ( i_blk = 0; i_blk <= Blk_List.Nblk-1 ; ++i_blk ) {
          Blk_List.message_noreschange_bottomeastcorner_recbuf[i_blk] = new double[1];
          Blk_List.message_noreschange_bottomeastcorner_recbuf[i_blk][0] = ZERO;
       } /* endfor */
    } /* endif */
    if (Blk_List.Nblk > 0 &&
        Blk_List.message_noreschange_bottomwestcorner_sendbuf == NULL) {
       Blk_List.message_noreschange_bottomwestcorner_sendbuf = new double*[Blk_List.Nblk];
       for ( i_blk = 0; i_blk <= Blk_List.Nblk-1 ; ++i_blk ) {
          Blk_List.message_noreschange_bottomwestcorner_sendbuf[i_blk] = new double[1];
          Blk_List.message_noreschange_bottomwestcorner_sendbuf[i_blk][0] = ZERO;
       } /* endfor */
    } /* endif */
    if (Blk_List.Nblk > 0 &&
        Blk_List.message_noreschange_bottomwestcorner_recbuf == NULL) {
       Blk_List.message_noreschange_bottomwestcorner_recbuf = new double*[Blk_List.Nblk];
       for ( i_blk = 0; i_blk <= Blk_List.Nblk-1 ; ++i_blk ) {
          Blk_List.message_noreschange_bottomwestcorner_recbuf[i_blk] = new double[1];
          Blk_List.message_noreschange_bottomwestcorner_recbuf[i_blk][0] = ZERO;
       } /* endfor */
    } /* endif */

    if (Blk_List.Nblk > 0 &&
        Blk_List.message_noreschange_bottomnorthwestcorner_sendbuf == NULL) {
       Blk_List.message_noreschange_bottomnorthwestcorner_sendbuf = new double*[Blk_List.Nblk];
       for ( i_blk = 0; i_blk <= Blk_List.Nblk-1 ; ++i_blk ) {
          Blk_List.message_noreschange_bottomnorthwestcorner_sendbuf[i_blk] = new double[1];
          Blk_List.message_noreschange_bottomnorthwestcorner_sendbuf[i_blk][0] = ZERO;
       } /* endfor */
    } /* endif */
    if (Blk_List.Nblk > 0 &&
        Blk_List.message_noreschange_bottomnorthwestcorner_recbuf == NULL) {
       Blk_List.message_noreschange_bottomnorthwestcorner_recbuf = new double*[Blk_List.Nblk];
       for ( i_blk = 0; i_blk <= Blk_List.Nblk-1 ; ++i_blk ) {
          Blk_List.message_noreschange_bottomnorthwestcorner_recbuf[i_blk] = new double[1];
          Blk_List.message_noreschange_bottomnorthwestcorner_recbuf[i_blk][0] = ZERO;
       } /* endfor */
    } /* endif */
    if (Blk_List.Nblk > 0 &&
        Blk_List.message_noreschange_bottomnortheastcorner_sendbuf == NULL) {
       Blk_List.message_noreschange_bottomnortheastcorner_sendbuf = new double*[Blk_List.Nblk];
       for ( i_blk = 0; i_blk <= Blk_List.Nblk-1 ; ++i_blk ) {
          Blk_List.message_noreschange_bottomnortheastcorner_sendbuf[i_blk] = new double[1];
          Blk_List.message_noreschange_bottomnortheastcorner_sendbuf[i_blk][0] = ZERO;
       } /* endfor */
    } /* endif */
    if (Blk_List.Nblk > 0 &&
        Blk_List.message_noreschange_bottomnortheastcorner_recbuf == NULL) {
       Blk_List.message_noreschange_bottomnortheastcorner_recbuf = new double*[Blk_List.Nblk];
       for ( i_blk = 0; i_blk <= Blk_List.Nblk-1 ; ++i_blk ) {
          Blk_List.message_noreschange_bottomnortheastcorner_recbuf[i_blk] = new double[1];
          Blk_List.message_noreschange_bottomnortheastcorner_recbuf[i_blk][0] = ZERO;
       } /* endfor */
    } /* endif */
    if (Blk_List.Nblk > 0 &&
        Blk_List.message_noreschange_bottomsoutheastcorner_sendbuf == NULL) {
       Blk_List.message_noreschange_bottomsoutheastcorner_sendbuf = new double*[Blk_List.Nblk];
       for ( i_blk = 0; i_blk <= Blk_List.Nblk-1 ; ++i_blk ) {
          Blk_List.message_noreschange_bottomsoutheastcorner_sendbuf[i_blk] = new double[1];
          Blk_List.message_noreschange_bottomsoutheastcorner_sendbuf[i_blk][0] = ZERO;
       } /* endfor */
    } /* endif */
    if (Blk_List.Nblk > 0 &&
        Blk_List.message_noreschange_bottomsoutheastcorner_recbuf == NULL) {
       Blk_List.message_noreschange_bottomsoutheastcorner_recbuf = new double*[Blk_List.Nblk];
       for ( i_blk = 0; i_blk <= Blk_List.Nblk-1 ; ++i_blk ) {
          Blk_List.message_noreschange_bottomsoutheastcorner_recbuf[i_blk] = new double[1];
          Blk_List.message_noreschange_bottomsoutheastcorner_recbuf[i_blk][0] = ZERO;
       } /* endfor */
    } /* endif */
    if (Blk_List.Nblk > 0 &&
        Blk_List.message_noreschange_bottomsouthwestcorner_sendbuf == NULL) {
       Blk_List.message_noreschange_bottomsouthwestcorner_sendbuf = new double*[Blk_List.Nblk];
       for ( i_blk = 0; i_blk <= Blk_List.Nblk-1 ; ++i_blk ) {
          Blk_List.message_noreschange_bottomsouthwestcorner_sendbuf[i_blk] = new double[1];
          Blk_List.message_noreschange_bottomsouthwestcorner_sendbuf[i_blk][0] = ZERO;
       } /* endfor */
    } /* endif */
    if (Blk_List.Nblk > 0 &&
        Blk_List.message_noreschange_bottomsouthwestcorner_recbuf == NULL) {
       Blk_List.message_noreschange_bottomsouthwestcorner_recbuf = new double*[Blk_List.Nblk];
       for ( i_blk = 0; i_blk <= Blk_List.Nblk-1 ; ++i_blk ) {
          Blk_List.message_noreschange_bottomsouthwestcorner_recbuf[i_blk] = new double[1];
          Blk_List.message_noreschange_bottomsouthwestcorner_recbuf[i_blk][0] = ZERO;
       } /* endfor */
    } /* endif */


    /* Reallocate memory for the send and receive buffers of each block. */
    
    for ( i_blk = 0; i_blk <= Blk_List.Nblk-1; ++i_blk ) {

        // Reallocate memory for Top face send and receive buffers.
       if (Blk_List.message_noreschange_topface_sendbuf[i_blk] != NULL) {
          delete []Blk_List.message_noreschange_topface_sendbuf[i_blk];
          Blk_List.message_noreschange_topface_sendbuf[i_blk] = NULL;
       } /* endif */
       if (Blk_List.message_noreschange_topface_recbuf[i_blk] != NULL) {
          delete []Blk_List.message_noreschange_topface_recbuf[i_blk];
          Blk_List.message_noreschange_topface_recbuf[i_blk] = NULL;
       } /* endif */
       if (Blk_List.Block[i_blk].used && 
           (Blk_List.Block[i_blk].nT == 1) &&
           (Blk_List.Block[i_blk].info.level == Blk_List.Block[i_blk].infoT[0].level)) {
          buffer_size = Blk_List.Block[i_blk].info.dimen.ghost*
	    (abs(Blk_List.Block[i_blk].info.dimen.i)+1)*
	    (abs(Blk_List.Block[i_blk].info.dimen.j)+1)*
	    Number_of_Solution_Variables+
	    4*Blk_List.Block[i_blk].infoT[0].dimen.ghost*abs(Blk_List.Block[i_blk].infoT[0].dimen.i)+
	    4*Blk_List.Block[i_blk].infoT[0].dimen.ghost*abs(Blk_List.Block[i_blk].infoT[0].dimen.j);
          buffer_size_neighbour = Blk_List.Block[i_blk].infoT[0].dimen.ghost*
                                  (abs(Blk_List.Block[i_blk].infoT[0].dimen.i)+1)*
                                  (abs(Blk_List.Block[i_blk].infoT[0].dimen.j)+1)*
                                  Number_of_Solution_Variables+
	    4*Blk_List.Block[i_blk].infoT[0].dimen.ghost*abs(Blk_List.Block[i_blk].infoT[0].dimen.i)+
	    4*Blk_List.Block[i_blk].infoT[0].dimen.ghost*abs(Blk_List.Block[i_blk].infoT[0].dimen.j);
	  // 	  cout<<"\nmessage_noreschange_topface buffer_size= "<<buffer_size<<" buffer_size_neighbour= "<<buffer_size_neighbour<<" For block "<<i_blk; cout.flush();
          Blk_List.message_noreschange_topface_sendbuf[i_blk] = 
             new double[buffer_size_neighbour];
          for (l = 0; l <= buffer_size_neighbour-1; ++l) 
	     Blk_List.message_noreschange_topface_sendbuf[i_blk][l] = ZERO;
          Blk_List.message_noreschange_topface_recbuf[i_blk] = 
             new double[buffer_size];
          for (l = 0; l <= buffer_size-1; ++l) 
	     Blk_List.message_noreschange_topface_recbuf[i_blk][l] = ZERO;
       } else {
          Blk_List.message_noreschange_topface_sendbuf[i_blk] = 
             new double[1];
          Blk_List.message_noreschange_topface_sendbuf[i_blk][0] = ZERO;
          Blk_List.message_noreschange_topface_recbuf[i_blk] = 
             new double[1];
          Blk_List.message_noreschange_topface_recbuf[i_blk][0] = ZERO;
       } /* endif */

        // Reallocate memory for Bottom face send and receive buffers.
       if (Blk_List.message_noreschange_bottomface_sendbuf[i_blk] != NULL) {
          delete []Blk_List.message_noreschange_bottomface_sendbuf[i_blk];
          Blk_List.message_noreschange_bottomface_sendbuf[i_blk] = NULL;
       } /* endif */
       if (Blk_List.message_noreschange_bottomface_recbuf[i_blk] != NULL) {
          delete []Blk_List.message_noreschange_bottomface_recbuf[i_blk];
          Blk_List.message_noreschange_bottomface_recbuf[i_blk] = NULL;
       } /* endif */
       if (Blk_List.Block[i_blk].used && 
           (Blk_List.Block[i_blk].nB == 1) &&
           (Blk_List.Block[i_blk].info.level == Blk_List.Block[i_blk].infoB[0].level)) {
          buffer_size = Blk_List.Block[i_blk].info.dimen.ghost*
                        (abs(Blk_List.Block[i_blk].info.dimen.i)+1)*
                        (abs(Blk_List.Block[i_blk].info.dimen.j)+1)*
                        Number_of_Solution_Variables+
	    4*Blk_List.Block[i_blk].infoB[0].dimen.ghost*abs(Blk_List.Block[i_blk].infoB[0].dimen.i)+
	    4*Blk_List.Block[i_blk].infoB[0].dimen.ghost*abs(Blk_List.Block[i_blk].infoB[0].dimen.j);
          buffer_size_neighbour = Blk_List.Block[i_blk].infoB[0].dimen.ghost*
                                  (abs(Blk_List.Block[i_blk].infoB[0].dimen.i)+1)*
                                  (abs(Blk_List.Block[i_blk].infoB[0].dimen.j)+1)*
                                  Number_of_Solution_Variables+
	    4*Blk_List.Block[i_blk].infoB[0].dimen.ghost*abs(Blk_List.Block[i_blk].infoB[0].dimen.i)+
	    4*Blk_List.Block[i_blk].infoB[0].dimen.ghost*abs(Blk_List.Block[i_blk].infoB[0].dimen.j);
// 	  cout<<"\nmessage_noreschange_bottomface buffer_size= "<<buffer_size<<" buffer_size_neighbour= "<<buffer_size_neighbour<<" For block "<<i_blk; cout.flush();
          Blk_List.message_noreschange_bottomface_sendbuf[i_blk] = 
             new double[buffer_size_neighbour];
          for (l = 0; l <= buffer_size_neighbour-1; ++l) 
	     Blk_List.message_noreschange_bottomface_sendbuf[i_blk][l] = ZERO;
          Blk_List.message_noreschange_bottomface_recbuf[i_blk] = 
             new double[buffer_size];
          for (l = 0; l <= buffer_size-1; ++l) 
	     Blk_List.message_noreschange_bottomface_recbuf[i_blk][l] = ZERO;
       } else {
          Blk_List.message_noreschange_bottomface_sendbuf[i_blk] = 
             new double[1];
          Blk_List.message_noreschange_bottomface_sendbuf[i_blk][0] = ZERO;
          Blk_List.message_noreschange_bottomface_recbuf[i_blk] = 
             new double[1];
          Blk_List.message_noreschange_bottomface_recbuf[i_blk][0] = ZERO;
       } /* endif */
 
      // Reallocate memory for North face send and receive buffers.
       if (Blk_List.message_noreschange_northface_sendbuf[i_blk] != NULL) {
          delete []Blk_List.message_noreschange_northface_sendbuf[i_blk];
          Blk_List.message_noreschange_northface_sendbuf[i_blk] = NULL;
       } /* endif */
       if (Blk_List.message_noreschange_northface_recbuf[i_blk] != NULL) {
          delete []Blk_List.message_noreschange_northface_recbuf[i_blk];
          Blk_List.message_noreschange_northface_recbuf[i_blk] = NULL;
       } /* endif */
       if (Blk_List.Block[i_blk].used && 
           (Blk_List.Block[i_blk].nN == 1) &&
           (Blk_List.Block[i_blk].info.level == Blk_List.Block[i_blk].infoN[0].level)) {
          buffer_size = Blk_List.Block[i_blk].info.dimen.ghost*
                        (abs(Blk_List.Block[i_blk].info.dimen.i)+1)*
                        (abs(Blk_List.Block[i_blk].info.dimen.k)+1)*
                        Number_of_Solution_Variables+
		  4*Blk_List.Block[i_blk].infoN[0].dimen.ghost*abs(Blk_List.Block[i_blk].infoN[0].dimen.i)+
		  4*Blk_List.Block[i_blk].infoN[0].dimen.ghost*abs(Blk_List.Block[i_blk].infoN[0].dimen.k);
          buffer_size_neighbour = Blk_List.Block[i_blk].infoN[0].dimen.ghost*
                                  (abs(Blk_List.Block[i_blk].infoN[0].dimen.i)+1)*
                                  (abs(Blk_List.Block[i_blk].infoN[0].dimen.k)+1)*
                                  Number_of_Solution_Variables+
		  4*Blk_List.Block[i_blk].infoN[0].dimen.ghost*abs(Blk_List.Block[i_blk].infoN[0].dimen.i)+
		  4*Blk_List.Block[i_blk].infoN[0].dimen.ghost*abs(Blk_List.Block[i_blk].infoN[0].dimen.k);
          Blk_List.message_noreschange_northface_sendbuf[i_blk] = 
             new double[buffer_size_neighbour];
          for (l = 0; l <= buffer_size_neighbour-1; ++l) 
	     Blk_List.message_noreschange_northface_sendbuf[i_blk][l] = ZERO;
          Blk_List.message_noreschange_northface_recbuf[i_blk] = 
             new double[buffer_size];
          for (l = 0; l <= buffer_size-1; ++l) 
	     Blk_List.message_noreschange_northface_recbuf[i_blk][l] = ZERO;
       } else {
          Blk_List.message_noreschange_northface_sendbuf[i_blk] = 
             new double[1];
          Blk_List.message_noreschange_northface_sendbuf[i_blk][0] = ZERO;
          Blk_List.message_noreschange_northface_recbuf[i_blk] = 
             new double[1];
          Blk_List.message_noreschange_northface_recbuf[i_blk][0] = ZERO;
       } /* endif */
       
       // Reallocate memory for South face send and receive buffers.
       if (Blk_List.message_noreschange_southface_sendbuf[i_blk] != NULL) {
          delete []Blk_List.message_noreschange_southface_sendbuf[i_blk];
          Blk_List.message_noreschange_southface_sendbuf[i_blk] = NULL;
       } /* endif */
       if (Blk_List.message_noreschange_southface_recbuf[i_blk] != NULL) {
          delete []Blk_List.message_noreschange_southface_recbuf[i_blk];
          Blk_List.message_noreschange_southface_recbuf[i_blk] = NULL;
       } /* endif */
       if (Blk_List.Block[i_blk].used && 
           (Blk_List.Block[i_blk].nS == 1) &&
           (Blk_List.Block[i_blk].info.level == Blk_List.Block[i_blk].infoS[0].level)) {
          buffer_size = Blk_List.Block[i_blk].info.dimen.ghost*
                        (abs(Blk_List.Block[i_blk].info.dimen.i)+1)*
                        (abs(Blk_List.Block[i_blk].info.dimen.k)+1)*
                        Number_of_Solution_Variables+
		  4*Blk_List.Block[i_blk].infoS[0].dimen.ghost*abs(Blk_List.Block[i_blk].infoS[0].dimen.i)+
		  4*Blk_List.Block[i_blk].infoS[0].dimen.ghost*abs(Blk_List.Block[i_blk].infoS[0].dimen.k);
          buffer_size_neighbour = Blk_List.Block[i_blk].infoS[0].dimen.ghost*
                                  (abs(Blk_List.Block[i_blk].infoS[0].dimen.i)+1)*
                                  (abs(Blk_List.Block[i_blk].infoS[0].dimen.k)+1)*
                                  Number_of_Solution_Variables+
		  4*Blk_List.Block[i_blk].infoS[0].dimen.ghost*abs(Blk_List.Block[i_blk].infoS[0].dimen.i)+
		  4*Blk_List.Block[i_blk].infoS[0].dimen.ghost*abs(Blk_List.Block[i_blk].infoS[0].dimen.k);
          Blk_List.message_noreschange_southface_sendbuf[i_blk] = 
             new double[buffer_size_neighbour];
          for (l = 0; l <= buffer_size_neighbour-1; ++l) 
	     Blk_List.message_noreschange_southface_sendbuf[i_blk][l] = ZERO;
          Blk_List.message_noreschange_southface_recbuf[i_blk] = 
             new double[buffer_size];
          for (l = 0; l <= buffer_size-1; ++l) 
	     Blk_List.message_noreschange_southface_recbuf[i_blk][l] = ZERO;
       } else {
          Blk_List.message_noreschange_southface_sendbuf[i_blk] = 
             new double[1];
          Blk_List.message_noreschange_southface_sendbuf[i_blk][0] = ZERO;
          Blk_List.message_noreschange_southface_recbuf[i_blk] = 
             new double[1];
          Blk_List.message_noreschange_southface_recbuf[i_blk][0] = ZERO;
       } /* endif */
       
       // Reallocate memory for East face send and receive buffers.
       if (Blk_List.message_noreschange_eastface_sendbuf[i_blk] != NULL) {
          delete []Blk_List.message_noreschange_eastface_sendbuf[i_blk];
          Blk_List.message_noreschange_eastface_sendbuf[i_blk] = NULL;
       } /* endif */
       if (Blk_List.message_noreschange_eastface_recbuf[i_blk] != NULL) {
          delete []Blk_List.message_noreschange_eastface_recbuf[i_blk];
          Blk_List.message_noreschange_eastface_recbuf[i_blk] = NULL;
       } /* endif */
       if (Blk_List.Block[i_blk].used && 
           (Blk_List.Block[i_blk].nE == 1) &&
           (Blk_List.Block[i_blk].info.level == Blk_List.Block[i_blk].infoE[0].level)) {
          buffer_size = Blk_List.Block[i_blk].info.dimen.ghost*
                        (abs(Blk_List.Block[i_blk].info.dimen.j)+1)*
                        (abs(Blk_List.Block[i_blk].info.dimen.k)+1)*
                        Number_of_Solution_Variables+
		  4*Blk_List.Block[i_blk].infoE[0].dimen.ghost*abs(Blk_List.Block[i_blk].infoE[0].dimen.j)+
		  4*Blk_List.Block[i_blk].infoE[0].dimen.ghost*abs(Blk_List.Block[i_blk].infoE[0].dimen.k);
          buffer_size_neighbour = Blk_List.Block[i_blk].infoE[0].dimen.ghost*
                                  (abs(Blk_List.Block[i_blk].infoE[0].dimen.j)+1)*
                                  (abs(Blk_List.Block[i_blk].infoE[0].dimen.k)+1)*
                                  Number_of_Solution_Variables+
		  4*Blk_List.Block[i_blk].infoE[0].dimen.ghost*abs(Blk_List.Block[i_blk].infoE[0].dimen.j)+
		  4*Blk_List.Block[i_blk].infoE[0].dimen.ghost*abs(Blk_List.Block[i_blk].infoE[0].dimen.k);
          Blk_List.message_noreschange_eastface_sendbuf[i_blk] = 
             new double[buffer_size_neighbour];
          for (l = 0; l <= buffer_size_neighbour-1; ++l) 
	     Blk_List.message_noreschange_eastface_sendbuf[i_blk][l] = ZERO;
          Blk_List.message_noreschange_eastface_recbuf[i_blk] = 
             new double[buffer_size];
          for (l = 0; l <= buffer_size-1; ++l) 
	     Blk_List.message_noreschange_eastface_recbuf[i_blk][l] = ZERO;
       } else {
          Blk_List.message_noreschange_eastface_sendbuf[i_blk] = 
             new double[1];
          Blk_List.message_noreschange_eastface_sendbuf[i_blk][0] = ZERO;
          Blk_List.message_noreschange_eastface_recbuf[i_blk] = 
             new double[1];
          Blk_List.message_noreschange_eastface_recbuf[i_blk][0] = ZERO;
       } /* endif */
       
       // Reallocate memory for West face send and receive buffers.
       if (Blk_List.message_noreschange_westface_sendbuf[i_blk] != NULL) {
          delete []Blk_List.message_noreschange_westface_sendbuf[i_blk];
          Blk_List.message_noreschange_westface_sendbuf[i_blk] = NULL;
       } /* endif */
       if (Blk_List.message_noreschange_westface_recbuf[i_blk] != NULL) {
          delete []Blk_List.message_noreschange_westface_recbuf[i_blk];
          Blk_List.message_noreschange_westface_recbuf[i_blk] = NULL;
       } /* endif */
       if (Blk_List.Block[i_blk].used && 
           (Blk_List.Block[i_blk].nW == 1) &&
           (Blk_List.Block[i_blk].info.level == Blk_List.Block[i_blk].infoW[0].level)) {
          buffer_size = Blk_List.Block[i_blk].info.dimen.ghost*
                        (abs(Blk_List.Block[i_blk].info.dimen.j)+1)*
                        (abs(Blk_List.Block[i_blk].info.dimen.k)+1)*
                        Number_of_Solution_Variables+
		  4*Blk_List.Block[i_blk].infoW[0].dimen.ghost*abs(Blk_List.Block[i_blk].infoW[0].dimen.j)+
		  4*Blk_List.Block[i_blk].infoW[0].dimen.ghost*abs(Blk_List.Block[i_blk].infoW[0].dimen.k);
          buffer_size_neighbour = Blk_List.Block[i_blk].infoW[0].dimen.ghost*
                                  (abs(Blk_List.Block[i_blk].infoW[0].dimen.j)+1)*
                                  (abs(Blk_List.Block[i_blk].infoW[0].dimen.k)+1)*
                                  Number_of_Solution_Variables+
		  4*Blk_List.Block[i_blk].infoW[0].dimen.ghost*abs(Blk_List.Block[i_blk].infoW[0].dimen.j)+
		  4*Blk_List.Block[i_blk].infoW[0].dimen.ghost*abs(Blk_List.Block[i_blk].infoW[0].dimen.k);
          Blk_List.message_noreschange_westface_sendbuf[i_blk] = 
             new double[buffer_size_neighbour];
          for (l = 0; l <= buffer_size_neighbour-1; ++l) 
	     Blk_List.message_noreschange_westface_sendbuf[i_blk][l] = ZERO;
          Blk_List.message_noreschange_westface_recbuf[i_blk] = 
             new double[buffer_size];
          for (l = 0; l <= buffer_size-1; ++l) 
	     Blk_List.message_noreschange_westface_recbuf[i_blk][l] = ZERO;
       } else {
          Blk_List.message_noreschange_westface_sendbuf[i_blk] = 
             new double[1];
          Blk_List.message_noreschange_westface_sendbuf[i_blk][0] = ZERO;
          Blk_List.message_noreschange_westface_recbuf[i_blk] = 
             new double[1];
          Blk_List.message_noreschange_westface_recbuf[i_blk][0] = ZERO;
       } /* endif */
       
       // Reallocate memory for North West corner send and receive buffers.
       if (Blk_List.message_noreschange_northwestcorner_sendbuf[i_blk] != NULL) {
          delete []Blk_List.message_noreschange_northwestcorner_sendbuf[i_blk];
          Blk_List.message_noreschange_northwestcorner_sendbuf[i_blk] = NULL;
       } /* endif */
       if (Blk_List.message_noreschange_northwestcorner_recbuf[i_blk] != NULL) {
          delete []Blk_List.message_noreschange_northwestcorner_recbuf[i_blk];
          Blk_List.message_noreschange_northwestcorner_recbuf[i_blk] = NULL;
       } /* endif */
       if (Blk_List.Block[i_blk].used && 
           (Blk_List.Block[i_blk].nNW == 1) &&
           (Blk_List.Block[i_blk].info.level == Blk_List.Block[i_blk].infoNW[0].level)) {
          buffer_size = sqr(Blk_List.Block[i_blk].info.dimen.ghost)*
                        (abs(Blk_List.Block[i_blk].info.dimen.k)+1)*
                        Number_of_Solution_Variables;
          buffer_size_neighbour = sqr(Blk_List.Block[i_blk].infoNW[0].dimen.ghost)*
                        (abs(Blk_List.Block[i_blk].info.dimen.k)+1)*
                                  Number_of_Solution_Variables;
          Blk_List.message_noreschange_northwestcorner_sendbuf[i_blk] = 
             new double[buffer_size_neighbour];
          for (l = 0; l <= buffer_size_neighbour-1; ++l) 
	     Blk_List.message_noreschange_northwestcorner_sendbuf[i_blk][l] = ZERO;
          Blk_List.message_noreschange_northwestcorner_recbuf[i_blk] = 
             new double[buffer_size];
          for (l = 0; l <= buffer_size-1; ++l) 
	     Blk_List.message_noreschange_northwestcorner_recbuf[i_blk][l] = ZERO;
       } else {
          Blk_List.message_noreschange_northwestcorner_sendbuf[i_blk] = 
             new double[1];
          Blk_List.message_noreschange_northwestcorner_sendbuf[i_blk][0] = ZERO; 
          Blk_List.message_noreschange_northwestcorner_recbuf[i_blk] = 
             new double[1];
          Blk_List.message_noreschange_northwestcorner_recbuf[i_blk][0] = ZERO;
       } /* endif */
       
       // Reallocate memory for North East corner send and receive buffers.
       if (Blk_List.message_noreschange_northeastcorner_sendbuf[i_blk] != NULL) {
          delete []Blk_List.message_noreschange_northeastcorner_sendbuf[i_blk];
          Blk_List.message_noreschange_northeastcorner_sendbuf[i_blk] = NULL;
       } /* endif */
       if (Blk_List.message_noreschange_northeastcorner_recbuf[i_blk] != NULL) {
          delete []Blk_List.message_noreschange_northeastcorner_recbuf[i_blk];
          Blk_List.message_noreschange_northeastcorner_recbuf[i_blk] = NULL;
       } /* endif */
       if (Blk_List.Block[i_blk].used && 
           (Blk_List.Block[i_blk].nNE == 1) &&
           (Blk_List.Block[i_blk].info.level == Blk_List.Block[i_blk].infoNE[0].level)) {
          buffer_size = sqr(Blk_List.Block[i_blk].info.dimen.ghost)*
                        (abs(Blk_List.Block[i_blk].info.dimen.k)+1)*
                        Number_of_Solution_Variables;
          buffer_size_neighbour = sqr(Blk_List.Block[i_blk].infoNE[0].dimen.ghost)*
                        (abs(Blk_List.Block[i_blk].info.dimen.k)+1)*
                                  Number_of_Solution_Variables;
          Blk_List.message_noreschange_northeastcorner_sendbuf[i_blk] = 
             new double[buffer_size_neighbour];
          for (l = 0; l <= buffer_size_neighbour-1; ++l) 
	     Blk_List.message_noreschange_northeastcorner_sendbuf[i_blk][l] = ZERO;
          Blk_List.message_noreschange_northeastcorner_recbuf[i_blk] = 
             new double[buffer_size];
          for (l = 0; l <= buffer_size-1; ++l) 
	     Blk_List.message_noreschange_northeastcorner_recbuf[i_blk][l] = ZERO;
       } else {
          Blk_List.message_noreschange_northeastcorner_sendbuf[i_blk] = 
             new double[1];
          Blk_List.message_noreschange_northeastcorner_sendbuf[i_blk][0] = ZERO;
          Blk_List.message_noreschange_northeastcorner_recbuf[i_blk] = 
             new double[1];
          Blk_List.message_noreschange_northeastcorner_recbuf[i_blk][0] = ZERO;
       } /* endif */
       
       // Reallocate memory for South East corner send and receive buffers.
       if (Blk_List.message_noreschange_southeastcorner_sendbuf[i_blk] != NULL) {
          delete []Blk_List.message_noreschange_southeastcorner_sendbuf[i_blk];
          Blk_List.message_noreschange_southeastcorner_sendbuf[i_blk] = NULL;
       } /* endif */
       if (Blk_List.message_noreschange_southeastcorner_recbuf[i_blk] != NULL) {
          delete []Blk_List.message_noreschange_southeastcorner_recbuf[i_blk];
          Blk_List.message_noreschange_southeastcorner_recbuf[i_blk] = NULL;
       } /* endif */
       if (Blk_List.Block[i_blk].used && 
           (Blk_List.Block[i_blk].nSE == 1) &&
           (Blk_List.Block[i_blk].info.level == Blk_List.Block[i_blk].infoSE[0].level)) {
          buffer_size = sqr(Blk_List.Block[i_blk].info.dimen.ghost)*
                        (abs(Blk_List.Block[i_blk].info.dimen.k)+1)*
                        Number_of_Solution_Variables;
          buffer_size_neighbour = sqr(Blk_List.Block[i_blk].infoSE[0].dimen.ghost)*
                        (abs(Blk_List.Block[i_blk].info.dimen.k)+1)*
                                  Number_of_Solution_Variables;
          Blk_List.message_noreschange_southeastcorner_sendbuf[i_blk] = 
             new double[buffer_size_neighbour];
          for (l = 0; l <= buffer_size_neighbour-1; ++l) 
	     Blk_List.message_noreschange_southeastcorner_sendbuf[i_blk][l] = ZERO;
          Blk_List.message_noreschange_southeastcorner_recbuf[i_blk] = 
             new double[buffer_size];
          for (l = 0; l <= buffer_size-1; ++l) 
	     Blk_List.message_noreschange_southeastcorner_recbuf[i_blk][l] = ZERO;
       } else {
          Blk_List.message_noreschange_southeastcorner_sendbuf[i_blk] = 
             new double[1];
          Blk_List.message_noreschange_southeastcorner_sendbuf[i_blk][0] = ZERO;
          Blk_List.message_noreschange_southeastcorner_recbuf[i_blk] = 
             new double[1];
          Blk_List.message_noreschange_southeastcorner_recbuf[i_blk][0] = ZERO; 
       } /* endif */
       
       // Reallocate memory for South West corner send and receive buffers.
       if (Blk_List.message_noreschange_southwestcorner_sendbuf[i_blk] != NULL) {
          delete []Blk_List.message_noreschange_southwestcorner_sendbuf[i_blk];
          Blk_List.message_noreschange_southwestcorner_sendbuf[i_blk] = NULL;
       } /* endif */
       if (Blk_List.message_noreschange_southwestcorner_recbuf[i_blk] != NULL) {
          delete []Blk_List.message_noreschange_southwestcorner_recbuf[i_blk];
          Blk_List.message_noreschange_southwestcorner_recbuf[i_blk] = NULL;
       } /* endif */
       if (Blk_List.Block[i_blk].used && 
           (Blk_List.Block[i_blk].nSW == 1) &&
           (Blk_List.Block[i_blk].info.level == Blk_List.Block[i_blk].infoSW[0].level)) {
          buffer_size = sqr(Blk_List.Block[i_blk].info.dimen.ghost)*
                        (abs(Blk_List.Block[i_blk].info.dimen.k)+1)*
                        Number_of_Solution_Variables;
          buffer_size_neighbour = sqr(Blk_List.Block[i_blk].infoSW[0].dimen.ghost)*
                        (abs(Blk_List.Block[i_blk].info.dimen.k)+1)*
                                  Number_of_Solution_Variables;
          Blk_List.message_noreschange_southwestcorner_sendbuf[i_blk] = 
             new double[buffer_size_neighbour];
          for (l = 0; l <= buffer_size_neighbour-1; ++l) 
	     Blk_List.message_noreschange_southwestcorner_sendbuf[i_blk][l] = ZERO;
          Blk_List.message_noreschange_southwestcorner_recbuf[i_blk] = 
             new double[buffer_size];
          for (l = 0; l <= buffer_size-1; ++l) 
	     Blk_List.message_noreschange_southwestcorner_recbuf[i_blk][l] = ZERO;
       } else {
          Blk_List.message_noreschange_southwestcorner_sendbuf[i_blk] = 
             new double[1];
          Blk_List.message_noreschange_southwestcorner_sendbuf[i_blk][0] = ZERO;
          Blk_List.message_noreschange_southwestcorner_recbuf[i_blk] = 
             new double[1];
          Blk_List.message_noreschange_southwestcorner_recbuf[i_blk][0] = ZERO;
       } /* endif */

    ///**********************************************

       // Reallocate memory for Top North corner send and receive buffers.
       if (Blk_List.message_noreschange_topnorthcorner_sendbuf[i_blk] != NULL) {
          delete []Blk_List.message_noreschange_topnorthcorner_sendbuf[i_blk];
          Blk_List.message_noreschange_topnorthcorner_sendbuf[i_blk] = NULL;
       } /* endif */
       if (Blk_List.message_noreschange_topnorthcorner_recbuf[i_blk] != NULL) {
          delete []Blk_List.message_noreschange_topnorthcorner_recbuf[i_blk];
          Blk_List.message_noreschange_topnorthcorner_recbuf[i_blk] = NULL;
       } /* endif */
       if (Blk_List.Block[i_blk].used && 
           (Blk_List.Block[i_blk].nTN == 1) &&
           (Blk_List.Block[i_blk].info.level == Blk_List.Block[i_blk].infoTN[0].level)) {
          buffer_size = sqr(Blk_List.Block[i_blk].info.dimen.ghost)*
                        (abs(Blk_List.Block[i_blk].info.dimen.i)+1)*
                        Number_of_Solution_Variables;
          buffer_size_neighbour = sqr(Blk_List.Block[i_blk].infoTN[0].dimen.ghost)*
                        (abs(Blk_List.Block[i_blk].info.dimen.i)+1)*
                                  Number_of_Solution_Variables;
          Blk_List.message_noreschange_topnorthcorner_sendbuf[i_blk] = 
             new double[buffer_size_neighbour];
          for (l = 0; l <= buffer_size_neighbour-1; ++l) 
	     Blk_List.message_noreschange_topnorthcorner_sendbuf[i_blk][l] = ZERO;
          Blk_List.message_noreschange_topnorthcorner_recbuf[i_blk] = 
             new double[buffer_size];
          for (l = 0; l <= buffer_size-1; ++l) 
	     Blk_List.message_noreschange_topnorthcorner_recbuf[i_blk][l] = ZERO;
       } else {
          Blk_List.message_noreschange_topnorthcorner_sendbuf[i_blk] = 
             new double[1];
          Blk_List.message_noreschange_topnorthcorner_sendbuf[i_blk][0] = ZERO; 
          Blk_List.message_noreschange_topnorthcorner_recbuf[i_blk] = 
             new double[1];
          Blk_List.message_noreschange_topnorthcorner_recbuf[i_blk][0] = ZERO;
       } /* endif */
       
       // Reallocate memory for Top South corner send and receive buffers.
       if (Blk_List.message_noreschange_topsouthcorner_sendbuf[i_blk] != NULL) {
          delete []Blk_List.message_noreschange_topsouthcorner_sendbuf[i_blk];
          Blk_List.message_noreschange_topsouthcorner_sendbuf[i_blk] = NULL;
       } /* endif */
       if (Blk_List.message_noreschange_topsouthcorner_recbuf[i_blk] != NULL) {
          delete []Blk_List.message_noreschange_topsouthcorner_recbuf[i_blk];
          Blk_List.message_noreschange_topsouthcorner_recbuf[i_blk] = NULL;
       } /* endif */
       if (Blk_List.Block[i_blk].used && 
           (Blk_List.Block[i_blk].nTS == 1) &&
           (Blk_List.Block[i_blk].info.level == Blk_List.Block[i_blk].infoTS[0].level)) {
          buffer_size = sqr(Blk_List.Block[i_blk].info.dimen.ghost)*
                        (abs(Blk_List.Block[i_blk].info.dimen.i)+1)*
                        Number_of_Solution_Variables;
          buffer_size_neighbour = sqr(Blk_List.Block[i_blk].infoTS[0].dimen.ghost)*
                        (abs(Blk_List.Block[i_blk].info.dimen.i)+1)*
                                  Number_of_Solution_Variables;
          Blk_List.message_noreschange_topsouthcorner_sendbuf[i_blk] = 
             new double[buffer_size_neighbour];
          for (l = 0; l <= buffer_size_neighbour-1; ++l) 
	     Blk_List.message_noreschange_topsouthcorner_sendbuf[i_blk][l] = ZERO;
          Blk_List.message_noreschange_topsouthcorner_recbuf[i_blk] = 
             new double[buffer_size];
          for (l = 0; l <= buffer_size-1; ++l) 
	     Blk_List.message_noreschange_topsouthcorner_recbuf[i_blk][l] = ZERO;
       } else {
          Blk_List.message_noreschange_topsouthcorner_sendbuf[i_blk] = 
             new double[1];
          Blk_List.message_noreschange_topsouthcorner_sendbuf[i_blk][0] = ZERO;
          Blk_List.message_noreschange_topsouthcorner_recbuf[i_blk] = 
             new double[1];
          Blk_List.message_noreschange_topsouthcorner_recbuf[i_blk][0] = ZERO;
       } /* endif */
       
       // Reallocate memory for Top East corner send and receive buffers.
       if (Blk_List.message_noreschange_topeastcorner_sendbuf[i_blk] != NULL) {
          delete []Blk_List.message_noreschange_topeastcorner_sendbuf[i_blk];
          Blk_List.message_noreschange_topeastcorner_sendbuf[i_blk] = NULL;
       } /* endif */
       if (Blk_List.message_noreschange_topeastcorner_recbuf[i_blk] != NULL) {
          delete []Blk_List.message_noreschange_topeastcorner_recbuf[i_blk];
          Blk_List.message_noreschange_topeastcorner_recbuf[i_blk] = NULL;
       } /* endif */
       if (Blk_List.Block[i_blk].used && 
           (Blk_List.Block[i_blk].nTE == 1) &&
           (Blk_List.Block[i_blk].info.level == Blk_List.Block[i_blk].infoTE[0].level)) {
          buffer_size = sqr(Blk_List.Block[i_blk].info.dimen.ghost)*
                        (abs(Blk_List.Block[i_blk].info.dimen.j)+1)*
                        Number_of_Solution_Variables;
          buffer_size_neighbour = sqr(Blk_List.Block[i_blk].infoTE[0].dimen.ghost)*
                        (abs(Blk_List.Block[i_blk].info.dimen.j)+1)*
                                  Number_of_Solution_Variables;
          Blk_List.message_noreschange_topeastcorner_sendbuf[i_blk] = 
             new double[buffer_size_neighbour];
          for (l = 0; l <= buffer_size_neighbour-1; ++l) 
	     Blk_List.message_noreschange_topeastcorner_sendbuf[i_blk][l] = ZERO;
          Blk_List.message_noreschange_topeastcorner_recbuf[i_blk] = 
             new double[buffer_size];
          for (l = 0; l <= buffer_size-1; ++l) 
	     Blk_List.message_noreschange_topeastcorner_recbuf[i_blk][l] = ZERO;
       } else {
          Blk_List.message_noreschange_topeastcorner_sendbuf[i_blk] = 
             new double[1];
          Blk_List.message_noreschange_topeastcorner_sendbuf[i_blk][0] = ZERO;
          Blk_List.message_noreschange_topeastcorner_recbuf[i_blk] = 
             new double[1];
          Blk_List.message_noreschange_topeastcorner_recbuf[i_blk][0] = ZERO; 
       } /* endif */
       
       // Reallocate memory for Top West corner send and receive buffers.
       if (Blk_List.message_noreschange_topwestcorner_sendbuf[i_blk] != NULL) {
          delete []Blk_List.message_noreschange_topwestcorner_sendbuf[i_blk];
          Blk_List.message_noreschange_topwestcorner_sendbuf[i_blk] = NULL;
       } /* endif */
       if (Blk_List.message_noreschange_topwestcorner_recbuf[i_blk] != NULL) {
          delete []Blk_List.message_noreschange_topwestcorner_recbuf[i_blk];
          Blk_List.message_noreschange_topwestcorner_recbuf[i_blk] = NULL;
       } /* endif */
       if (Blk_List.Block[i_blk].used && 
           (Blk_List.Block[i_blk].nTW == 1) &&
           (Blk_List.Block[i_blk].info.level == Blk_List.Block[i_blk].infoTW[0].level)) {
          buffer_size = sqr(Blk_List.Block[i_blk].info.dimen.ghost)*Blk_List.Block[i_blk].info.dimen.ghost*
                        (abs(Blk_List.Block[i_blk].info.dimen.j)+1)*
                        Number_of_Solution_Variables;
          buffer_size_neighbour = sqr(Blk_List.Block[i_blk].infoTW[0].dimen.ghost)*
	    Blk_List.Block[i_blk].infoTW[0].dimen.ghost*
                        (abs(Blk_List.Block[i_blk].info.dimen.j)+1)*
                                  Number_of_Solution_Variables;
          Blk_List.message_noreschange_topwestcorner_sendbuf[i_blk] = 
             new double[buffer_size_neighbour];
          for (l = 0; l <= buffer_size_neighbour-1; ++l) 
	     Blk_List.message_noreschange_topwestcorner_sendbuf[i_blk][l] = ZERO;
          Blk_List.message_noreschange_topwestcorner_recbuf[i_blk] = 
             new double[buffer_size];
          for (l = 0; l <= buffer_size-1; ++l) 
	     Blk_List.message_noreschange_topwestcorner_recbuf[i_blk][l] = ZERO;
       } else {
          Blk_List.message_noreschange_topwestcorner_sendbuf[i_blk] = 
             new double[1];
          Blk_List.message_noreschange_topwestcorner_sendbuf[i_blk][0] = ZERO;
          Blk_List.message_noreschange_topwestcorner_recbuf[i_blk] = 
             new double[1];
          Blk_List.message_noreschange_topwestcorner_recbuf[i_blk][0] = ZERO;
       } /* endif */

       // Reallocate memory for Top North West corner send and receive buffers.
       if (Blk_List.message_noreschange_topnorthwestcorner_sendbuf[i_blk] != NULL) {
          delete []Blk_List.message_noreschange_topnorthwestcorner_sendbuf[i_blk];
          Blk_List.message_noreschange_topnorthwestcorner_sendbuf[i_blk] = NULL;
       } /* endif */
       if (Blk_List.message_noreschange_topnorthwestcorner_recbuf[i_blk] != NULL) {
          delete []Blk_List.message_noreschange_topnorthwestcorner_recbuf[i_blk];
          Blk_List.message_noreschange_topnorthwestcorner_recbuf[i_blk] = NULL;
       } /* endif */
       if (Blk_List.Block[i_blk].used && 
           (Blk_List.Block[i_blk].nTNW == 1) &&
           (Blk_List.Block[i_blk].info.level == Blk_List.Block[i_blk].infoTNW[0].level)) {
          buffer_size = sqr(Blk_List.Block[i_blk].info.dimen.ghost)*Blk_List.Block[i_blk].info.dimen.ghost*
                        Number_of_Solution_Variables;
          buffer_size_neighbour = sqr(Blk_List.Block[i_blk].infoTNW[0].dimen.ghost)*
	    Blk_List.Block[i_blk].infoTNW[0].dimen.ghost*
                                  Number_of_Solution_Variables;
          Blk_List.message_noreschange_topnorthwestcorner_sendbuf[i_blk] = 
             new double[buffer_size_neighbour];
          for (l = 0; l <= buffer_size_neighbour-1; ++l) 
	     Blk_List.message_noreschange_topnorthwestcorner_sendbuf[i_blk][l] = ZERO;
          Blk_List.message_noreschange_topnorthwestcorner_recbuf[i_blk] = 
             new double[buffer_size];
          for (l = 0; l <= buffer_size-1; ++l) 
	     Blk_List.message_noreschange_topnorthwestcorner_recbuf[i_blk][l] = ZERO;
       } else {
          Blk_List.message_noreschange_topnorthwestcorner_sendbuf[i_blk] = 
             new double[1];
          Blk_List.message_noreschange_topnorthwestcorner_sendbuf[i_blk][0] = ZERO; 
          Blk_List.message_noreschange_topnorthwestcorner_recbuf[i_blk] = 
             new double[1];
          Blk_List.message_noreschange_topnorthwestcorner_recbuf[i_blk][0] = ZERO;
       } /* endif */
       
       // Reallocate memory for Top North East corner send and receive buffers.
       if (Blk_List.message_noreschange_topnortheastcorner_sendbuf[i_blk] != NULL) {
          delete []Blk_List.message_noreschange_topnortheastcorner_sendbuf[i_blk];
          Blk_List.message_noreschange_topnortheastcorner_sendbuf[i_blk] = NULL;
       } /* endif */
       if (Blk_List.message_noreschange_topnortheastcorner_recbuf[i_blk] != NULL) {
          delete []Blk_List.message_noreschange_topnortheastcorner_recbuf[i_blk];
          Blk_List.message_noreschange_topnortheastcorner_recbuf[i_blk] = NULL;
       } /* endif */
       if (Blk_List.Block[i_blk].used && 
           (Blk_List.Block[i_blk].nTNE == 1) &&
           (Blk_List.Block[i_blk].info.level == Blk_List.Block[i_blk].infoTNE[0].level)) {
          buffer_size = sqr(Blk_List.Block[i_blk].info.dimen.ghost)*Blk_List.Block[i_blk].info.dimen.ghost*
                        Number_of_Solution_Variables;
          buffer_size_neighbour = sqr(Blk_List.Block[i_blk].infoTNE[0].dimen.ghost)*
	    Blk_List.Block[i_blk].infoTNE[0].dimen.ghost*
                                  Number_of_Solution_Variables;
          Blk_List.message_noreschange_topnortheastcorner_sendbuf[i_blk] = 
             new double[buffer_size_neighbour];
          for (l = 0; l <= buffer_size_neighbour-1; ++l) 
	     Blk_List.message_noreschange_topnortheastcorner_sendbuf[i_blk][l] = ZERO;
          Blk_List.message_noreschange_topnortheastcorner_recbuf[i_blk] = 
             new double[buffer_size];
          for (l = 0; l <= buffer_size-1; ++l) 
	     Blk_List.message_noreschange_topnortheastcorner_recbuf[i_blk][l] = ZERO;
       } else {
          Blk_List.message_noreschange_topnortheastcorner_sendbuf[i_blk] = 
             new double[1];
          Blk_List.message_noreschange_topnortheastcorner_sendbuf[i_blk][0] = ZERO;
          Blk_List.message_noreschange_topnortheastcorner_recbuf[i_blk] = 
             new double[1];
          Blk_List.message_noreschange_topnortheastcorner_recbuf[i_blk][0] = ZERO;
       } /* endif */
       
       // Reallocate memory for Top South East corner send and receive buffers.
       if (Blk_List.message_noreschange_topsoutheastcorner_sendbuf[i_blk] != NULL) {
          delete []Blk_List.message_noreschange_topsoutheastcorner_sendbuf[i_blk];
          Blk_List.message_noreschange_topsoutheastcorner_sendbuf[i_blk] = NULL;
       } /* endif */
       if (Blk_List.message_noreschange_topsoutheastcorner_recbuf[i_blk] != NULL) {
          delete []Blk_List.message_noreschange_topsoutheastcorner_recbuf[i_blk];
          Blk_List.message_noreschange_topsoutheastcorner_recbuf[i_blk] = NULL;
       } /* endif */
       if (Blk_List.Block[i_blk].used && 
           (Blk_List.Block[i_blk].nTSE == 1) &&
           (Blk_List.Block[i_blk].info.level == Blk_List.Block[i_blk].infoTSE[0].level)) {
          buffer_size = sqr(Blk_List.Block[i_blk].info.dimen.ghost)*Blk_List.Block[i_blk].info.dimen.ghost*
                        Number_of_Solution_Variables;
          buffer_size_neighbour = sqr(Blk_List.Block[i_blk].infoTSE[0].dimen.ghost)*
	    Blk_List.Block[i_blk].infoTSE[0].dimen.ghost*
                                  Number_of_Solution_Variables;
          Blk_List.message_noreschange_topsoutheastcorner_sendbuf[i_blk] = 
             new double[buffer_size_neighbour];
          for (l = 0; l <= buffer_size_neighbour-1; ++l) 
	     Blk_List.message_noreschange_topsoutheastcorner_sendbuf[i_blk][l] = ZERO;
          Blk_List.message_noreschange_topsoutheastcorner_recbuf[i_blk] = 
             new double[buffer_size];
          for (l = 0; l <= buffer_size-1; ++l) 
	     Blk_List.message_noreschange_topsoutheastcorner_recbuf[i_blk][l] = ZERO;
       } else {
          Blk_List.message_noreschange_topsoutheastcorner_sendbuf[i_blk] = 
             new double[1];
          Blk_List.message_noreschange_topsoutheastcorner_sendbuf[i_blk][0] = ZERO;
          Blk_List.message_noreschange_topsoutheastcorner_recbuf[i_blk] = 
             new double[1];
          Blk_List.message_noreschange_topsoutheastcorner_recbuf[i_blk][0] = ZERO; 
       } /* endif */
       
       // Reallocate memory for Top South West corner send and receive buffers.
       if (Blk_List.message_noreschange_topsouthwestcorner_sendbuf[i_blk] != NULL) {
          delete []Blk_List.message_noreschange_topsouthwestcorner_sendbuf[i_blk];
          Blk_List.message_noreschange_topsouthwestcorner_sendbuf[i_blk] = NULL;
       } /* endif */
       if (Blk_List.message_noreschange_topsouthwestcorner_recbuf[i_blk] != NULL) {
          delete []Blk_List.message_noreschange_topsouthwestcorner_recbuf[i_blk];
          Blk_List.message_noreschange_topsouthwestcorner_recbuf[i_blk] = NULL;
       } /* endif */
       if (Blk_List.Block[i_blk].used && 
           (Blk_List.Block[i_blk].nTSW == 1) &&
           (Blk_List.Block[i_blk].info.level == Blk_List.Block[i_blk].infoTSW[0].level)) {
          buffer_size = sqr(Blk_List.Block[i_blk].info.dimen.ghost)*Blk_List.Block[i_blk].info.dimen.ghost*
                        Number_of_Solution_Variables;
          buffer_size_neighbour = sqr(Blk_List.Block[i_blk].infoTSW[0].dimen.ghost)*
	    Blk_List.Block[i_blk].infoTSW[0].dimen.ghost*
                                  Number_of_Solution_Variables;
          Blk_List.message_noreschange_topsouthwestcorner_sendbuf[i_blk] = 
             new double[buffer_size_neighbour];
          for (l = 0; l <= buffer_size_neighbour-1; ++l) 
	     Blk_List.message_noreschange_topsouthwestcorner_sendbuf[i_blk][l] = ZERO;
          Blk_List.message_noreschange_topsouthwestcorner_recbuf[i_blk] = 
             new double[buffer_size];
          for (l = 0; l <= buffer_size-1; ++l) 
	     Blk_List.message_noreschange_topsouthwestcorner_recbuf[i_blk][l] = ZERO;
       } else {
          Blk_List.message_noreschange_topsouthwestcorner_sendbuf[i_blk] = 
             new double[1];
          Blk_List.message_noreschange_topsouthwestcorner_sendbuf[i_blk][0] = ZERO;
          Blk_List.message_noreschange_topsouthwestcorner_recbuf[i_blk] = 
             new double[1];
          Blk_List.message_noreschange_topsouthwestcorner_recbuf[i_blk][0] = ZERO;
       } /* endif */


       /*****************************/
       /*******************************/
       /********************************/
      // Reallocate memory for Bottom North corner send and receive buffers.
       if (Blk_List.message_noreschange_bottomnorthcorner_sendbuf[i_blk] != NULL) {
          delete []Blk_List.message_noreschange_bottomnorthcorner_sendbuf[i_blk];
          Blk_List.message_noreschange_bottomnorthcorner_sendbuf[i_blk] = NULL;
       } /* endif */
       if (Blk_List.message_noreschange_bottomnorthcorner_recbuf[i_blk] != NULL) {
          delete []Blk_List.message_noreschange_bottomnorthcorner_recbuf[i_blk];
          Blk_List.message_noreschange_bottomnorthcorner_recbuf[i_blk] = NULL;
       } /* endif */
       if (Blk_List.Block[i_blk].used && 
           (Blk_List.Block[i_blk].nBN == 1) &&
           (Blk_List.Block[i_blk].info.level == Blk_List.Block[i_blk].infoBN[0].level)) {
          buffer_size = sqr(Blk_List.Block[i_blk].info.dimen.ghost)*Blk_List.Block[i_blk].info.dimen.ghost*
                        (abs(Blk_List.Block[i_blk].info.dimen.i)+1)*
                        Number_of_Solution_Variables;
          buffer_size_neighbour = sqr(Blk_List.Block[i_blk].infoBN[0].dimen.ghost)*
	    Blk_List.Block[i_blk].infoBN[0].dimen.ghost*
                        (abs(Blk_List.Block[i_blk].info.dimen.i)+1)*
                                  Number_of_Solution_Variables;
          Blk_List.message_noreschange_bottomnorthcorner_sendbuf[i_blk] = 
             new double[buffer_size_neighbour];
          for (l = 0; l <= buffer_size_neighbour-1; ++l) 
	     Blk_List.message_noreschange_bottomnorthcorner_sendbuf[i_blk][l] = ZERO;
          Blk_List.message_noreschange_bottomnorthcorner_recbuf[i_blk] = 
             new double[buffer_size];
          for (l = 0; l <= buffer_size-1; ++l) 
	     Blk_List.message_noreschange_bottomnorthcorner_recbuf[i_blk][l] = ZERO;
       } else {
          Blk_List.message_noreschange_bottomnorthcorner_sendbuf[i_blk] = 
             new double[1];
          Blk_List.message_noreschange_bottomnorthcorner_sendbuf[i_blk][0] = ZERO; 
          Blk_List.message_noreschange_bottomnorthcorner_recbuf[i_blk] = 
             new double[1];
          Blk_List.message_noreschange_bottomnorthcorner_recbuf[i_blk][0] = ZERO;
       } /* endif */
       
       // Reallocate memory for Bottom South corner send and receive buffers.
       if (Blk_List.message_noreschange_bottomsouthcorner_sendbuf[i_blk] != NULL) {
          delete []Blk_List.message_noreschange_bottomsouthcorner_sendbuf[i_blk];
          Blk_List.message_noreschange_bottomsouthcorner_sendbuf[i_blk] = NULL;
       } /* endif */
       if (Blk_List.message_noreschange_bottomsouthcorner_recbuf[i_blk] != NULL) {
          delete []Blk_List.message_noreschange_bottomsouthcorner_recbuf[i_blk];
          Blk_List.message_noreschange_bottomsouthcorner_recbuf[i_blk] = NULL;
       } /* endif */
       if (Blk_List.Block[i_blk].used && 
           (Blk_List.Block[i_blk].nBS == 1) &&
           (Blk_List.Block[i_blk].info.level == Blk_List.Block[i_blk].infoBS[0].level)) {
          buffer_size = sqr(Blk_List.Block[i_blk].info.dimen.ghost)*Blk_List.Block[i_blk].info.dimen.ghost*
                        (abs(Blk_List.Block[i_blk].info.dimen.i)+1)*
                        Number_of_Solution_Variables;
          buffer_size_neighbour = sqr(Blk_List.Block[i_blk].infoBS[0].dimen.ghost)*
	    Blk_List.Block[i_blk].infoBS[0].dimen.ghost*
                        (abs(Blk_List.Block[i_blk].info.dimen.i)+1)*
                                  Number_of_Solution_Variables;
          Blk_List.message_noreschange_bottomsouthcorner_sendbuf[i_blk] = 
             new double[buffer_size_neighbour];
          for (l = 0; l <= buffer_size_neighbour-1; ++l) 
	     Blk_List.message_noreschange_bottomsouthcorner_sendbuf[i_blk][l] = ZERO;
          Blk_List.message_noreschange_bottomsouthcorner_recbuf[i_blk] = 
             new double[buffer_size];
          for (l = 0; l <= buffer_size-1; ++l) 
	     Blk_List.message_noreschange_bottomsouthcorner_recbuf[i_blk][l] = ZERO;
       } else {
          Blk_List.message_noreschange_bottomsouthcorner_sendbuf[i_blk] = 
             new double[1];
          Blk_List.message_noreschange_bottomsouthcorner_sendbuf[i_blk][0] = ZERO;
          Blk_List.message_noreschange_bottomsouthcorner_recbuf[i_blk] = 
             new double[1];
          Blk_List.message_noreschange_bottomsouthcorner_recbuf[i_blk][0] = ZERO;
       } /* endif */
       
       // Reallocate memory for Bottom East corner send and receive buffers.
       if (Blk_List.message_noreschange_bottomeastcorner_sendbuf[i_blk] != NULL) {
          delete []Blk_List.message_noreschange_bottomeastcorner_sendbuf[i_blk];
          Blk_List.message_noreschange_bottomeastcorner_sendbuf[i_blk] = NULL;
       } /* endif */
       if (Blk_List.message_noreschange_bottomeastcorner_recbuf[i_blk] != NULL) {
          delete []Blk_List.message_noreschange_bottomeastcorner_recbuf[i_blk];
          Blk_List.message_noreschange_bottomeastcorner_recbuf[i_blk] = NULL;
       } /* endif */
       if (Blk_List.Block[i_blk].used && 
           (Blk_List.Block[i_blk].nBE == 1) &&
           (Blk_List.Block[i_blk].info.level == Blk_List.Block[i_blk].infoBE[0].level)) {
          buffer_size = sqr(Blk_List.Block[i_blk].info.dimen.ghost)*Blk_List.Block[i_blk].info.dimen.ghost*
                        (abs(Blk_List.Block[i_blk].info.dimen.j)+1)*
                        Number_of_Solution_Variables;
          buffer_size_neighbour = sqr(Blk_List.Block[i_blk].infoBE[0].dimen.ghost)*
	    Blk_List.Block[i_blk].infoBE[0].dimen.ghost*
                        (abs(Blk_List.Block[i_blk].info.dimen.j)+1)*
                                  Number_of_Solution_Variables;
          Blk_List.message_noreschange_bottomeastcorner_sendbuf[i_blk] = 
             new double[buffer_size_neighbour];
          for (l = 0; l <= buffer_size_neighbour-1; ++l) 
	     Blk_List.message_noreschange_bottomeastcorner_sendbuf[i_blk][l] = ZERO;
          Blk_List.message_noreschange_bottomeastcorner_recbuf[i_blk] = 
             new double[buffer_size];
          for (l = 0; l <= buffer_size-1; ++l) 
	     Blk_List.message_noreschange_bottomeastcorner_recbuf[i_blk][l] = ZERO;
       } else {
          Blk_List.message_noreschange_bottomeastcorner_sendbuf[i_blk] = 
             new double[1];
          Blk_List.message_noreschange_bottomeastcorner_sendbuf[i_blk][0] = ZERO;
          Blk_List.message_noreschange_bottomeastcorner_recbuf[i_blk] = 
             new double[1];
          Blk_List.message_noreschange_bottomeastcorner_recbuf[i_blk][0] = ZERO; 
       } /* endif */
       
       // Reallocate memory for Bottom West corner send and receive buffers.
       if (Blk_List.message_noreschange_bottomwestcorner_sendbuf[i_blk] != NULL) {
          delete []Blk_List.message_noreschange_bottomwestcorner_sendbuf[i_blk];
          Blk_List.message_noreschange_bottomwestcorner_sendbuf[i_blk] = NULL;
       } /* endif */
       if (Blk_List.message_noreschange_bottomwestcorner_recbuf[i_blk] != NULL) {
          delete []Blk_List.message_noreschange_bottomwestcorner_recbuf[i_blk];
          Blk_List.message_noreschange_bottomwestcorner_recbuf[i_blk] = NULL;
       } /* endif */
       if (Blk_List.Block[i_blk].used && 
           (Blk_List.Block[i_blk].nBW == 1) &&
           (Blk_List.Block[i_blk].info.level == Blk_List.Block[i_blk].infoBW[0].level)) {
          buffer_size = sqr(Blk_List.Block[i_blk].info.dimen.ghost)*Blk_List.Block[i_blk].info.dimen.ghost*
                        (abs(Blk_List.Block[i_blk].info.dimen.j)+1)*
                        Number_of_Solution_Variables;
          buffer_size_neighbour = sqr(Blk_List.Block[i_blk].infoBW[0].dimen.ghost)*
	    Blk_List.Block[i_blk].infoBW[0].dimen.ghost*
                        (abs(Blk_List.Block[i_blk].info.dimen.j)+1)*
                                  Number_of_Solution_Variables;
          Blk_List.message_noreschange_bottomwestcorner_sendbuf[i_blk] = 
             new double[buffer_size_neighbour];
          for (l = 0; l <= buffer_size_neighbour-1; ++l) 
	     Blk_List.message_noreschange_bottomwestcorner_sendbuf[i_blk][l] = ZERO;
          Blk_List.message_noreschange_bottomwestcorner_recbuf[i_blk] = 
             new double[buffer_size];
          for (l = 0; l <= buffer_size-1; ++l) 
	     Blk_List.message_noreschange_bottomwestcorner_recbuf[i_blk][l] = ZERO;
       } else {
          Blk_List.message_noreschange_bottomwestcorner_sendbuf[i_blk] = 
             new double[1];
          Blk_List.message_noreschange_bottomwestcorner_sendbuf[i_blk][0] = ZERO;
          Blk_List.message_noreschange_bottomwestcorner_recbuf[i_blk] = 
             new double[1];
          Blk_List.message_noreschange_bottomwestcorner_recbuf[i_blk][0] = ZERO;
       } /* endif */

       // Reallocate memory for Bottom North West corner send and receive buffers.
       if (Blk_List.message_noreschange_bottomnorthwestcorner_sendbuf[i_blk] != NULL) {
          delete []Blk_List.message_noreschange_bottomnorthwestcorner_sendbuf[i_blk];
          Blk_List.message_noreschange_bottomnorthwestcorner_sendbuf[i_blk] = NULL;
       } /* endif */
       if (Blk_List.message_noreschange_bottomnorthwestcorner_recbuf[i_blk] != NULL) {
          delete []Blk_List.message_noreschange_bottomnorthwestcorner_recbuf[i_blk];
          Blk_List.message_noreschange_bottomnorthwestcorner_recbuf[i_blk] = NULL;
       } /* endif */
       if (Blk_List.Block[i_blk].used && 
           (Blk_List.Block[i_blk].nBNW == 1) &&
           (Blk_List.Block[i_blk].info.level == Blk_List.Block[i_blk].infoBNW[0].level)) {
          buffer_size = sqr(Blk_List.Block[i_blk].info.dimen.ghost)*Blk_List.Block[i_blk].info.dimen.ghost*
                        Number_of_Solution_Variables;
          buffer_size_neighbour = sqr(Blk_List.Block[i_blk].infoBNW[0].dimen.ghost)*
	    Blk_List.Block[i_blk].infoBNW[0].dimen.ghost*
                                  Number_of_Solution_Variables;
          Blk_List.message_noreschange_bottomnorthwestcorner_sendbuf[i_blk] = 
             new double[buffer_size_neighbour];
          for (l = 0; l <= buffer_size_neighbour-1; ++l) 
	     Blk_List.message_noreschange_bottomnorthwestcorner_sendbuf[i_blk][l] = ZERO;
          Blk_List.message_noreschange_bottomnorthwestcorner_recbuf[i_blk] = 
             new double[buffer_size];
          for (l = 0; l <= buffer_size-1; ++l) 
	     Blk_List.message_noreschange_bottomnorthwestcorner_recbuf[i_blk][l] = ZERO;
       } else {
          Blk_List.message_noreschange_bottomnorthwestcorner_sendbuf[i_blk] = 
             new double[1];
          Blk_List.message_noreschange_bottomnorthwestcorner_sendbuf[i_blk][0] = ZERO; 
          Blk_List.message_noreschange_bottomnorthwestcorner_recbuf[i_blk] = 
             new double[1];
          Blk_List.message_noreschange_bottomnorthwestcorner_recbuf[i_blk][0] = ZERO;
       } /* endif */
       
       // Reallocate memory for Bottom North East corner send and receive buffers.
       if (Blk_List.message_noreschange_bottomnortheastcorner_sendbuf[i_blk] != NULL) {
          delete []Blk_List.message_noreschange_bottomnortheastcorner_sendbuf[i_blk];
          Blk_List.message_noreschange_bottomnortheastcorner_sendbuf[i_blk] = NULL;
       } /* endif */
       if (Blk_List.message_noreschange_bottomnortheastcorner_recbuf[i_blk] != NULL) {
          delete []Blk_List.message_noreschange_bottomnortheastcorner_recbuf[i_blk];
          Blk_List.message_noreschange_bottomnortheastcorner_recbuf[i_blk] = NULL;
       } /* endif */
       if (Blk_List.Block[i_blk].used && 
           (Blk_List.Block[i_blk].nBNE == 1) &&
           (Blk_List.Block[i_blk].info.level == Blk_List.Block[i_blk].infoBNE[0].level)) {
          buffer_size = sqr(Blk_List.Block[i_blk].info.dimen.ghost)*Blk_List.Block[i_blk].info.dimen.ghost*
                        Number_of_Solution_Variables;
          buffer_size_neighbour = sqr(Blk_List.Block[i_blk].infoBNE[0].dimen.ghost)*
	    Blk_List.Block[i_blk].infoBNE[0].dimen.ghost*
                                  Number_of_Solution_Variables;
          Blk_List.message_noreschange_bottomnortheastcorner_sendbuf[i_blk] = 
             new double[buffer_size_neighbour];
          for (l = 0; l <= buffer_size_neighbour-1; ++l) 
	     Blk_List.message_noreschange_bottomnortheastcorner_sendbuf[i_blk][l] = ZERO;
          Blk_List.message_noreschange_bottomnortheastcorner_recbuf[i_blk] = 
             new double[buffer_size];
          for (l = 0; l <= buffer_size-1; ++l) 
	     Blk_List.message_noreschange_bottomnortheastcorner_recbuf[i_blk][l] = ZERO;
       } else {
          Blk_List.message_noreschange_bottomnortheastcorner_sendbuf[i_blk] = 
             new double[1];
          Blk_List.message_noreschange_bottomnortheastcorner_sendbuf[i_blk][0] = ZERO;
          Blk_List.message_noreschange_bottomnortheastcorner_recbuf[i_blk] = 
             new double[1];
          Blk_List.message_noreschange_bottomnortheastcorner_recbuf[i_blk][0] = ZERO;
       } /* endif */
       
       // Reallocate memory for Bottom South East corner send and receive buffers.
       if (Blk_List.message_noreschange_bottomsoutheastcorner_sendbuf[i_blk] != NULL) {
          delete []Blk_List.message_noreschange_bottomsoutheastcorner_sendbuf[i_blk];
          Blk_List.message_noreschange_bottomsoutheastcorner_sendbuf[i_blk] = NULL;
       } /* endif */
       if (Blk_List.message_noreschange_bottomsoutheastcorner_recbuf[i_blk] != NULL) {
          delete []Blk_List.message_noreschange_bottomsoutheastcorner_recbuf[i_blk];
          Blk_List.message_noreschange_bottomsoutheastcorner_recbuf[i_blk] = NULL;
       } /* endif */
       if (Blk_List.Block[i_blk].used && 
           (Blk_List.Block[i_blk].nBSE == 1) &&
           (Blk_List.Block[i_blk].info.level == Blk_List.Block[i_blk].infoBSE[0].level)) {
          buffer_size = sqr(Blk_List.Block[i_blk].info.dimen.ghost)*Blk_List.Block[i_blk].info.dimen.ghost*
                        Number_of_Solution_Variables;
          buffer_size_neighbour = sqr(Blk_List.Block[i_blk].infoBSE[0].dimen.ghost)*
	    Blk_List.Block[i_blk].infoBSE[0].dimen.ghost*
                                  Number_of_Solution_Variables;
          Blk_List.message_noreschange_bottomsoutheastcorner_sendbuf[i_blk] = 
             new double[buffer_size_neighbour];
          for (l = 0; l <= buffer_size_neighbour-1; ++l) 
	     Blk_List.message_noreschange_bottomsoutheastcorner_sendbuf[i_blk][l] = ZERO;
          Blk_List.message_noreschange_bottomsoutheastcorner_recbuf[i_blk] = 
             new double[buffer_size];
          for (l = 0; l <= buffer_size-1; ++l) 
	     Blk_List.message_noreschange_bottomsoutheastcorner_recbuf[i_blk][l] = ZERO;
       } else {
          Blk_List.message_noreschange_bottomsoutheastcorner_sendbuf[i_blk] = 
             new double[1];
          Blk_List.message_noreschange_bottomsoutheastcorner_sendbuf[i_blk][0] = ZERO;
          Blk_List.message_noreschange_bottomsoutheastcorner_recbuf[i_blk] = 
             new double[1];
          Blk_List.message_noreschange_bottomsoutheastcorner_recbuf[i_blk][0] = ZERO; 
       } /* endif */
       
       // Reallocate memory for Bottom South West corner send and receive buffers.
       if (Blk_List.message_noreschange_bottomsouthwestcorner_sendbuf[i_blk] != NULL) {
          delete []Blk_List.message_noreschange_bottomsouthwestcorner_sendbuf[i_blk];
          Blk_List.message_noreschange_bottomsouthwestcorner_sendbuf[i_blk] = NULL;
       } /* endif */
       if (Blk_List.message_noreschange_bottomsouthwestcorner_recbuf[i_blk] != NULL) {
          delete []Blk_List.message_noreschange_bottomsouthwestcorner_recbuf[i_blk];
          Blk_List.message_noreschange_bottomsouthwestcorner_recbuf[i_blk] = NULL;
       } /* endif */
       if (Blk_List.Block[i_blk].used && 
           (Blk_List.Block[i_blk].nBSW == 1) &&
           (Blk_List.Block[i_blk].info.level == Blk_List.Block[i_blk].infoBSW[0].level)) {
          buffer_size = sqr(Blk_List.Block[i_blk].info.dimen.ghost)*Blk_List.Block[i_blk].info.dimen.ghost*
                        Number_of_Solution_Variables;
          buffer_size_neighbour = sqr(Blk_List.Block[i_blk].infoBSW[0].dimen.ghost)*
	    Blk_List.Block[i_blk].infoBSW[0].dimen.ghost*
                                  Number_of_Solution_Variables;
          Blk_List.message_noreschange_bottomsouthwestcorner_sendbuf[i_blk] = 
             new double[buffer_size_neighbour];
          for (l = 0; l <= buffer_size_neighbour-1; ++l) 
	     Blk_List.message_noreschange_bottomsouthwestcorner_sendbuf[i_blk][l] = ZERO;
          Blk_List.message_noreschange_bottomsouthwestcorner_recbuf[i_blk] = 
             new double[buffer_size];
          for (l = 0; l <= buffer_size-1; ++l) 
	     Blk_List.message_noreschange_bottomsouthwestcorner_recbuf[i_blk][l] = ZERO;
       } else {
          Blk_List.message_noreschange_bottomsouthwestcorner_sendbuf[i_blk] = 
             new double[1];
          Blk_List.message_noreschange_bottomsouthwestcorner_sendbuf[i_blk][0] = ZERO;
          Blk_List.message_noreschange_bottomsouthwestcorner_recbuf[i_blk] = 
             new double[1];
          Blk_List.message_noreschange_bottomsouthwestcorner_recbuf[i_blk][0] = ZERO;
       } /* endif */



    } /* endfor */
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

    cout<<"\n**************************** NOTE Allocate_Message_Buffers_ResChange THIS HAS NOT BEEN TESTED********************************/";

    return;

    /* Ensure that memory for the block index of the send and receive buffers
       has been allocated. */

    if (Blk_List.Nblk > 0 &&
        Blk_List.message_reschange_topface_sendbuf == NULL) {
       Blk_List.message_reschange_topface_sendbuf = new double**[Blk_List.Nblk];
       for ( i_blk = 0; i_blk <= Blk_List.Nblk-1 ; ++i_blk ) {
          Blk_List.message_reschange_topface_sendbuf[i_blk] = new double*[2];
          for (j_neigh = 0; j_neigh < 2 ; ++j_neigh) {
             Blk_List.message_reschange_topface_sendbuf[i_blk][j_neigh] = new double[1];
             Blk_List.message_reschange_topface_sendbuf[i_blk][j_neigh][0] = ZERO;
	  } /* endfor */
       } /* endfor */
    } /* endif */
    if (Blk_List.Nblk > 0 &&
        Blk_List.message_reschange_topface_recbuf == NULL) {
       Blk_List.message_reschange_topface_recbuf = new double**[Blk_List.Nblk];
       for ( i_blk = 0; i_blk <= Blk_List.Nblk-1 ; ++i_blk ) {
          Blk_List.message_reschange_topface_recbuf[i_blk] = new double*[2];
          for (j_neigh = 0; j_neigh < 2 ; ++j_neigh) {
             Blk_List.message_reschange_topface_recbuf[i_blk][j_neigh] = new double[1];
             Blk_List.message_reschange_topface_recbuf[i_blk][j_neigh][0] = ZERO;
	  } /* endfor */
       } /* endfor */
    } /* endif */

    if (Blk_List.Nblk > 0 &&
        Blk_List.message_reschange_bottomface_sendbuf == NULL) {
       Blk_List.message_reschange_bottomface_sendbuf = new double**[Blk_List.Nblk];
       for ( i_blk = 0; i_blk <= Blk_List.Nblk-1 ; ++i_blk ) {
          Blk_List.message_reschange_bottomface_sendbuf[i_blk] = new double*[2];
          for (j_neigh = 0; j_neigh < 2 ; ++j_neigh) {
             Blk_List.message_reschange_bottomface_sendbuf[i_blk][j_neigh] = new double[1];
             Blk_List.message_reschange_bottomface_sendbuf[i_blk][j_neigh][0] = ZERO;
	  } /* endfor */
       } /* endfor */
    } /* endif */
    if (Blk_List.Nblk > 0 &&
        Blk_List.message_reschange_bottomface_recbuf == NULL) {
       Blk_List.message_reschange_bottomface_recbuf = new double**[Blk_List.Nblk];
       for ( i_blk = 0; i_blk <= Blk_List.Nblk-1 ; ++i_blk ) {
          Blk_List.message_reschange_bottomface_recbuf[i_blk] = new double*[2];
          for (j_neigh = 0; j_neigh < 2 ; ++j_neigh) {
             Blk_List.message_reschange_bottomface_recbuf[i_blk][j_neigh] = new double[1];
             Blk_List.message_reschange_bottomface_recbuf[i_blk][j_neigh][0] = ZERO;
	  } /* endfor */
       } /* endfor */
    } /* endif */


    if (Blk_List.Nblk > 0 &&
        Blk_List.message_reschange_northface_sendbuf == NULL) {
       Blk_List.message_reschange_northface_sendbuf = new double**[Blk_List.Nblk];
       for ( i_blk = 0; i_blk <= Blk_List.Nblk-1 ; ++i_blk ) {
          Blk_List.message_reschange_northface_sendbuf[i_blk] = new double*[2];
          for (j_neigh = 0; j_neigh < 2 ; ++j_neigh) {
             Blk_List.message_reschange_northface_sendbuf[i_blk][j_neigh] = new double[1];
             Blk_List.message_reschange_northface_sendbuf[i_blk][j_neigh][0] = ZERO;
	  } /* endfor */
       } /* endfor */
    } /* endif */
    if (Blk_List.Nblk > 0 &&
        Blk_List.message_reschange_northface_recbuf == NULL) {
       Blk_List.message_reschange_northface_recbuf = new double**[Blk_List.Nblk];
       for ( i_blk = 0; i_blk <= Blk_List.Nblk-1 ; ++i_blk ) {
          Blk_List.message_reschange_northface_recbuf[i_blk] = new double*[2];
          for (j_neigh = 0; j_neigh < 2 ; ++j_neigh) {
             Blk_List.message_reschange_northface_recbuf[i_blk][j_neigh] = new double[1];
             Blk_List.message_reschange_northface_recbuf[i_blk][j_neigh][0] = ZERO;
	  } /* endfor */
       } /* endfor */
    } /* endif */
    if (Blk_List.Nblk > 0 &&
        Blk_List.message_reschange_southface_sendbuf == NULL) {
       Blk_List.message_reschange_southface_sendbuf = new double**[Blk_List.Nblk];
       for ( i_blk = 0; i_blk <= Blk_List.Nblk-1 ; ++i_blk ) {
          Blk_List.message_reschange_southface_sendbuf[i_blk] = new double*[2];
          for (j_neigh = 0; j_neigh < 2 ; ++j_neigh) {
             Blk_List.message_reschange_southface_sendbuf[i_blk][j_neigh] = new double[1];
             Blk_List.message_reschange_southface_sendbuf[i_blk][j_neigh][0] = ZERO;
          } /* endfor */
       } /* endfor */
    } /* endif */
    if (Blk_List.Nblk > 0 &&
        Blk_List.message_reschange_southface_recbuf == NULL) {
       Blk_List.message_reschange_southface_recbuf = new double**[Blk_List.Nblk];
       for ( i_blk = 0; i_blk <= Blk_List.Nblk-1 ; ++i_blk ) {
          Blk_List.message_reschange_southface_recbuf[i_blk] = new double*[2];
          for (j_neigh = 0; j_neigh < 2 ; ++j_neigh) {
             Blk_List.message_reschange_southface_recbuf[i_blk][j_neigh] = new double[1];
             Blk_List.message_reschange_southface_recbuf[i_blk][j_neigh][0] = ZERO;
          } /* endfor */
       } /* endfor */
    } /* endif */ 
    if (Blk_List.Nblk > 0 &&
        Blk_List.message_reschange_eastface_sendbuf == NULL) {
       Blk_List.message_reschange_eastface_sendbuf = new double**[Blk_List.Nblk];
       for ( i_blk = 0; i_blk <= Blk_List.Nblk-1 ; ++i_blk ) {
          Blk_List.message_reschange_eastface_sendbuf[i_blk] = new double*[2];
          for (j_neigh = 0; j_neigh < 2 ; ++j_neigh) {
             Blk_List.message_reschange_eastface_sendbuf[i_blk][j_neigh] = new double[1];
             Blk_List.message_reschange_eastface_sendbuf[i_blk][j_neigh][0] = ZERO;
          } /* endfor */
       } /* endfor */
    } /* endif */
    if (Blk_List.Nblk > 0 &&
        Blk_List.message_reschange_eastface_recbuf == NULL) {
       Blk_List.message_reschange_eastface_recbuf = new double**[Blk_List.Nblk];
       for ( i_blk = 0; i_blk <= Blk_List.Nblk-1 ; ++i_blk ) {
          Blk_List.message_reschange_eastface_recbuf[i_blk] = new double*[2];
          for (j_neigh = 0; j_neigh < 2 ; ++j_neigh) {
             Blk_List.message_reschange_eastface_recbuf[i_blk][j_neigh] = new double[1];
             Blk_List.message_reschange_eastface_recbuf[i_blk][j_neigh][0] = ZERO;
          } /* endfor */
       } /* endfor */
    } /* endif */
    if (Blk_List.Nblk > 0 &&
        Blk_List.message_reschange_westface_sendbuf == NULL) {
       Blk_List.message_reschange_westface_sendbuf = new double**[Blk_List.Nblk];
       for ( i_blk = 0; i_blk <= Blk_List.Nblk-1 ; ++i_blk ) {
          Blk_List.message_reschange_westface_sendbuf[i_blk] = new double*[2];
          for (j_neigh = 0; j_neigh < 2 ; ++j_neigh) {
             Blk_List.message_reschange_westface_sendbuf[i_blk][j_neigh] = new double[1];
             Blk_List.message_reschange_westface_sendbuf[i_blk][j_neigh][0] = ZERO;
          } /* endfor */
       } /* endfor */
    } /* endif */
    if (Blk_List.Nblk > 0 &&
        Blk_List.message_reschange_westface_recbuf == NULL) {
       Blk_List.message_reschange_westface_recbuf = new double**[Blk_List.Nblk];
       for ( i_blk = 0; i_blk <= Blk_List.Nblk-1 ; ++i_blk ) {
          Blk_List.message_reschange_westface_recbuf[i_blk] = new double*[2];
          for (j_neigh = 0; j_neigh < 2 ; ++j_neigh) {
             Blk_List.message_reschange_westface_recbuf[i_blk][j_neigh] = new double[1];
             Blk_List.message_reschange_westface_recbuf[i_blk][j_neigh][0] = ZERO;
          } /* endfor */
       } /* endfor */
    } /* endif */
    if (Blk_List.Nblk > 0 &&
        Blk_List.message_reschange_northwestcorner_sendbuf == NULL) {
       Blk_List.message_reschange_northwestcorner_sendbuf = new double*[Blk_List.Nblk];
       for ( i_blk = 0; i_blk <= Blk_List.Nblk-1 ; ++i_blk ) {
          Blk_List.message_reschange_northwestcorner_sendbuf[i_blk] = new double[1];
          Blk_List.message_reschange_northwestcorner_sendbuf[i_blk][0] = ZERO;
       } /* endfor */
    } /* endif */
    if (Blk_List.Nblk > 0 &&
        Blk_List.message_reschange_northwestcorner_recbuf == NULL) {
       Blk_List.message_reschange_northwestcorner_recbuf = new double*[Blk_List.Nblk];
       for ( i_blk = 0; i_blk <= Blk_List.Nblk-1 ; ++i_blk ) {
          Blk_List.message_reschange_northwestcorner_recbuf[i_blk] = new double[1];
          Blk_List.message_reschange_northwestcorner_recbuf[i_blk][0] = ZERO;
       } /* endfor */
    } /* endif */
    if (Blk_List.Nblk > 0 &&
        Blk_List.message_reschange_northeastcorner_sendbuf == NULL) {
       Blk_List.message_reschange_northeastcorner_sendbuf = new double*[Blk_List.Nblk];
       for ( i_blk = 0; i_blk <= Blk_List.Nblk-1 ; ++i_blk ) {
          Blk_List.message_reschange_northeastcorner_sendbuf[i_blk] = new double[1];
          Blk_List.message_reschange_northeastcorner_sendbuf[i_blk][0] = ZERO;
       } /* endfor */
    } /* endif */
    if (Blk_List.Nblk > 0 &&
        Blk_List.message_reschange_northeastcorner_recbuf == NULL) {
       Blk_List.message_reschange_northeastcorner_recbuf = new double*[Blk_List.Nblk];
       for ( i_blk = 0; i_blk <= Blk_List.Nblk-1 ; ++i_blk ) {
          Blk_List.message_reschange_northeastcorner_recbuf[i_blk] = new double[1];
          Blk_List.message_reschange_northeastcorner_recbuf[i_blk][0] = ZERO;
       } /* endfor */
    } /* endif */
    if (Blk_List.Nblk > 0 &&
        Blk_List.message_reschange_southeastcorner_sendbuf == NULL) {
       Blk_List.message_reschange_southeastcorner_sendbuf = new double*[Blk_List.Nblk];
       for ( i_blk = 0; i_blk <= Blk_List.Nblk-1 ; ++i_blk ) {
          Blk_List.message_reschange_southeastcorner_sendbuf[i_blk] = new double[1];
          Blk_List.message_reschange_southeastcorner_sendbuf[i_blk][0] = ZERO;
       } /* endfor */
    } /* endif */
    if (Blk_List.Nblk > 0 &&
        Blk_List.message_reschange_southeastcorner_recbuf == NULL) {
       Blk_List.message_reschange_southeastcorner_recbuf = new double*[Blk_List.Nblk];
       for ( i_blk = 0; i_blk <= Blk_List.Nblk-1 ; ++i_blk ) {
          Blk_List.message_reschange_southeastcorner_recbuf[i_blk] = new double[1];
          Blk_List.message_reschange_southeastcorner_recbuf[i_blk][0] = ZERO;
       } /* endfor */
    } /* endif */
    if (Blk_List.Nblk > 0 &&
        Blk_List.message_reschange_southwestcorner_sendbuf == NULL) {
       Blk_List.message_reschange_southwestcorner_sendbuf = new double*[Blk_List.Nblk];
       for ( i_blk = 0; i_blk <= Blk_List.Nblk-1 ; ++i_blk ) {
          Blk_List.message_reschange_southwestcorner_sendbuf[i_blk] = new double[1];
          Blk_List.message_reschange_southwestcorner_sendbuf[i_blk][0] = ZERO;
       } /* endfor */
    } /* endif */
    if (Blk_List.Nblk > 0 &&
        Blk_List.message_reschange_southwestcorner_recbuf == NULL) {
       Blk_List.message_reschange_southwestcorner_recbuf = new double*[Blk_List.Nblk];
       for ( i_blk = 0; i_blk <= Blk_List.Nblk-1 ; ++i_blk ) {
          Blk_List.message_reschange_southwestcorner_recbuf[i_blk] = new double[1];
          Blk_List.message_reschange_southwestcorner_recbuf[i_blk][0] = ZERO;
       } /* endfor */
    } /* endif */


    /**/

    if (Blk_List.Nblk > 0 &&
        Blk_List.message_reschange_topnorthcorner_sendbuf == NULL) {
       Blk_List.message_reschange_topnorthcorner_sendbuf = new double*[Blk_List.Nblk];
       for ( i_blk = 0; i_blk <= Blk_List.Nblk-1 ; ++i_blk ) {
          Blk_List.message_reschange_topnorthcorner_sendbuf[i_blk] = new double[1];
          Blk_List.message_reschange_topnorthcorner_sendbuf[i_blk][0] = ZERO;
       } /* endfor */
    } /* endif */
    if (Blk_List.Nblk > 0 &&
        Blk_List.message_reschange_topnorthcorner_recbuf == NULL) {
       Blk_List.message_reschange_topnorthcorner_recbuf = new double*[Blk_List.Nblk];
       for ( i_blk = 0; i_blk <= Blk_List.Nblk-1 ; ++i_blk ) {
          Blk_List.message_reschange_topnorthcorner_recbuf[i_blk] = new double[1];
          Blk_List.message_reschange_topnorthcorner_recbuf[i_blk][0] = ZERO;
       } /* endfor */
    } /* endif */
    if (Blk_List.Nblk > 0 &&
        Blk_List.message_reschange_topsouthcorner_sendbuf == NULL) {
       Blk_List.message_reschange_topsouthcorner_sendbuf = new double*[Blk_List.Nblk];
       for ( i_blk = 0; i_blk <= Blk_List.Nblk-1 ; ++i_blk ) {
          Blk_List.message_reschange_topsouthcorner_sendbuf[i_blk] = new double[1];
          Blk_List.message_reschange_topsouthcorner_sendbuf[i_blk][0] = ZERO;
       } /* endfor */
    } /* endif */
    if (Blk_List.Nblk > 0 &&
        Blk_List.message_reschange_topsouthcorner_recbuf == NULL) {
       Blk_List.message_reschange_topsouthcorner_recbuf = new double*[Blk_List.Nblk];
       for ( i_blk = 0; i_blk <= Blk_List.Nblk-1 ; ++i_blk ) {
          Blk_List.message_reschange_topsouthcorner_recbuf[i_blk] = new double[1];
          Blk_List.message_reschange_topsouthcorner_recbuf[i_blk][0] = ZERO;
       } /* endfor */
    } /* endif */
    if (Blk_List.Nblk > 0 &&
        Blk_List.message_reschange_topeastcorner_sendbuf == NULL) {
       Blk_List.message_reschange_topeastcorner_sendbuf = new double*[Blk_List.Nblk];
       for ( i_blk = 0; i_blk <= Blk_List.Nblk-1 ; ++i_blk ) {
          Blk_List.message_reschange_topeastcorner_sendbuf[i_blk] = new double[1];
          Blk_List.message_reschange_topeastcorner_sendbuf[i_blk][0] = ZERO;
       } /* endfor */
    } /* endif */
    if (Blk_List.Nblk > 0 &&
        Blk_List.message_reschange_topeastcorner_recbuf == NULL) {
       Blk_List.message_reschange_topeastcorner_recbuf = new double*[Blk_List.Nblk];
       for ( i_blk = 0; i_blk <= Blk_List.Nblk-1 ; ++i_blk ) {
          Blk_List.message_reschange_topeastcorner_recbuf[i_blk] = new double[1];
          Blk_List.message_reschange_topeastcorner_recbuf[i_blk][0] = ZERO;
       } /* endfor */
    } /* endif */
    if (Blk_List.Nblk > 0 &&
        Blk_List.message_reschange_topwestcorner_sendbuf == NULL) {
       Blk_List.message_reschange_topwestcorner_sendbuf = new double*[Blk_List.Nblk];
       for ( i_blk = 0; i_blk <= Blk_List.Nblk-1 ; ++i_blk ) {
          Blk_List.message_reschange_topwestcorner_sendbuf[i_blk] = new double[1];
          Blk_List.message_reschange_topwestcorner_sendbuf[i_blk][0] = ZERO;
       } /* endfor */
    } /* endif */
    if (Blk_List.Nblk > 0 &&
        Blk_List.message_reschange_topwestcorner_recbuf == NULL) {
       Blk_List.message_reschange_topwestcorner_recbuf = new double*[Blk_List.Nblk];
       for ( i_blk = 0; i_blk <= Blk_List.Nblk-1 ; ++i_blk ) {
          Blk_List.message_reschange_topwestcorner_recbuf[i_blk] = new double[1];
          Blk_List.message_reschange_topwestcorner_recbuf[i_blk][0] = ZERO;
       } /* endfor */
    } /* endif */

    if (Blk_List.Nblk > 0 &&
        Blk_List.message_reschange_topnorthwestcorner_sendbuf == NULL) {
       Blk_List.message_reschange_topnorthwestcorner_sendbuf = new double*[Blk_List.Nblk];
       for ( i_blk = 0; i_blk <= Blk_List.Nblk-1 ; ++i_blk ) {
          Blk_List.message_reschange_topnorthwestcorner_sendbuf[i_blk] = new double[1];
          Blk_List.message_reschange_topnorthwestcorner_sendbuf[i_blk][0] = ZERO;
       } /* endfor */
    } /* endif */
    if (Blk_List.Nblk > 0 &&
        Blk_List.message_reschange_topnorthwestcorner_recbuf == NULL) {
       Blk_List.message_reschange_topnorthwestcorner_recbuf = new double*[Blk_List.Nblk];
       for ( i_blk = 0; i_blk <= Blk_List.Nblk-1 ; ++i_blk ) {
          Blk_List.message_reschange_topnorthwestcorner_recbuf[i_blk] = new double[1];
          Blk_List.message_reschange_topnorthwestcorner_recbuf[i_blk][0] = ZERO;
       } /* endfor */
    } /* endif */
    if (Blk_List.Nblk > 0 &&
        Blk_List.message_reschange_topnortheastcorner_sendbuf == NULL) {
       Blk_List.message_reschange_topnortheastcorner_sendbuf = new double*[Blk_List.Nblk];
       for ( i_blk = 0; i_blk <= Blk_List.Nblk-1 ; ++i_blk ) {
          Blk_List.message_reschange_topnortheastcorner_sendbuf[i_blk] = new double[1];
          Blk_List.message_reschange_topnortheastcorner_sendbuf[i_blk][0] = ZERO;
       } /* endfor */
    } /* endif */
    if (Blk_List.Nblk > 0 &&
        Blk_List.message_reschange_topnortheastcorner_recbuf == NULL) {
       Blk_List.message_reschange_topnortheastcorner_recbuf = new double*[Blk_List.Nblk];
       for ( i_blk = 0; i_blk <= Blk_List.Nblk-1 ; ++i_blk ) {
          Blk_List.message_reschange_topnortheastcorner_recbuf[i_blk] = new double[1];
          Blk_List.message_reschange_topnortheastcorner_recbuf[i_blk][0] = ZERO;
       } /* endfor */
    } /* endif */
    if (Blk_List.Nblk > 0 &&
        Blk_List.message_reschange_topsoutheastcorner_sendbuf == NULL) {
       Blk_List.message_reschange_topsoutheastcorner_sendbuf = new double*[Blk_List.Nblk];
       for ( i_blk = 0; i_blk <= Blk_List.Nblk-1 ; ++i_blk ) {
          Blk_List.message_reschange_topsoutheastcorner_sendbuf[i_blk] = new double[1];
          Blk_List.message_reschange_topsoutheastcorner_sendbuf[i_blk][0] = ZERO;
       } /* endfor */
    } /* endif */
    if (Blk_List.Nblk > 0 &&
        Blk_List.message_reschange_topsoutheastcorner_recbuf == NULL) {
       Blk_List.message_reschange_topsoutheastcorner_recbuf = new double*[Blk_List.Nblk];
       for ( i_blk = 0; i_blk <= Blk_List.Nblk-1 ; ++i_blk ) {
          Blk_List.message_reschange_topsoutheastcorner_recbuf[i_blk] = new double[1];
          Blk_List.message_reschange_topsoutheastcorner_recbuf[i_blk][0] = ZERO;
       } /* endfor */
    } /* endif */
    if (Blk_List.Nblk > 0 &&
        Blk_List.message_reschange_topsouthwestcorner_sendbuf == NULL) {
       Blk_List.message_reschange_topsouthwestcorner_sendbuf = new double*[Blk_List.Nblk];
       for ( i_blk = 0; i_blk <= Blk_List.Nblk-1 ; ++i_blk ) {
          Blk_List.message_reschange_topsouthwestcorner_sendbuf[i_blk] = new double[1];
          Blk_List.message_reschange_topsouthwestcorner_sendbuf[i_blk][0] = ZERO;
       } /* endfor */
    } /* endif */
    if (Blk_List.Nblk > 0 &&
        Blk_List.message_reschange_topsouthwestcorner_recbuf == NULL) {
       Blk_List.message_reschange_topsouthwestcorner_recbuf = new double*[Blk_List.Nblk];
       for ( i_blk = 0; i_blk <= Blk_List.Nblk-1 ; ++i_blk ) {
          Blk_List.message_reschange_topsouthwestcorner_recbuf[i_blk] = new double[1];
          Blk_List.message_reschange_topsouthwestcorner_recbuf[i_blk][0] = ZERO;
       } /* endfor */
    } /* endif */
 
    /**/

    if (Blk_List.Nblk > 0 &&
        Blk_List.message_reschange_bottomnorthcorner_sendbuf == NULL) {
       Blk_List.message_reschange_bottomnorthcorner_sendbuf = new double*[Blk_List.Nblk];
       for ( i_blk = 0; i_blk <= Blk_List.Nblk-1 ; ++i_blk ) {
          Blk_List.message_reschange_bottomnorthcorner_sendbuf[i_blk] = new double[1];
          Blk_List.message_reschange_bottomnorthcorner_sendbuf[i_blk][0] = ZERO;
       } /* endfor */
    } /* endif */
    if (Blk_List.Nblk > 0 &&
        Blk_List.message_reschange_bottomnorthcorner_recbuf == NULL) {
       Blk_List.message_reschange_bottomnorthcorner_recbuf = new double*[Blk_List.Nblk];
       for ( i_blk = 0; i_blk <= Blk_List.Nblk-1 ; ++i_blk ) {
          Blk_List.message_reschange_bottomnorthcorner_recbuf[i_blk] = new double[1];
          Blk_List.message_reschange_bottomnorthcorner_recbuf[i_blk][0] = ZERO;
       } /* endfor */
    } /* endif */
    if (Blk_List.Nblk > 0 &&
        Blk_List.message_reschange_bottomsouthcorner_sendbuf == NULL) {
       Blk_List.message_reschange_bottomsouthcorner_sendbuf = new double*[Blk_List.Nblk];
       for ( i_blk = 0; i_blk <= Blk_List.Nblk-1 ; ++i_blk ) {
          Blk_List.message_reschange_bottomsouthcorner_sendbuf[i_blk] = new double[1];
          Blk_List.message_reschange_bottomsouthcorner_sendbuf[i_blk][0] = ZERO;
       } /* endfor */
    } /* endif */
    if (Blk_List.Nblk > 0 &&
        Blk_List.message_reschange_bottomsouthcorner_recbuf == NULL) {
       Blk_List.message_reschange_bottomsouthcorner_recbuf = new double*[Blk_List.Nblk];
       for ( i_blk = 0; i_blk <= Blk_List.Nblk-1 ; ++i_blk ) {
          Blk_List.message_reschange_bottomsouthcorner_recbuf[i_blk] = new double[1];
          Blk_List.message_reschange_bottomsouthcorner_recbuf[i_blk][0] = ZERO;
       } /* endfor */
    } /* endif */
    if (Blk_List.Nblk > 0 &&
        Blk_List.message_reschange_bottomeastcorner_sendbuf == NULL) {
       Blk_List.message_reschange_bottomeastcorner_sendbuf = new double*[Blk_List.Nblk];
       for ( i_blk = 0; i_blk <= Blk_List.Nblk-1 ; ++i_blk ) {
          Blk_List.message_reschange_bottomeastcorner_sendbuf[i_blk] = new double[1];
          Blk_List.message_reschange_bottomeastcorner_sendbuf[i_blk][0] = ZERO;
       } /* endfor */
    } /* endif */
    if (Blk_List.Nblk > 0 &&
        Blk_List.message_reschange_bottomeastcorner_recbuf == NULL) {
       Blk_List.message_reschange_bottomeastcorner_recbuf = new double*[Blk_List.Nblk];
       for ( i_blk = 0; i_blk <= Blk_List.Nblk-1 ; ++i_blk ) {
          Blk_List.message_reschange_bottomeastcorner_recbuf[i_blk] = new double[1];
          Blk_List.message_reschange_bottomeastcorner_recbuf[i_blk][0] = ZERO;
       } /* endfor */
    } /* endif */
    if (Blk_List.Nblk > 0 &&
        Blk_List.message_reschange_bottomwestcorner_sendbuf == NULL) {
       Blk_List.message_reschange_bottomwestcorner_sendbuf = new double*[Blk_List.Nblk];
       for ( i_blk = 0; i_blk <= Blk_List.Nblk-1 ; ++i_blk ) {
          Blk_List.message_reschange_bottomwestcorner_sendbuf[i_blk] = new double[1];
          Blk_List.message_reschange_bottomwestcorner_sendbuf[i_blk][0] = ZERO;
       } /* endfor */
    } /* endif */
    if (Blk_List.Nblk > 0 &&
        Blk_List.message_reschange_bottomwestcorner_recbuf == NULL) {
       Blk_List.message_reschange_bottomwestcorner_recbuf = new double*[Blk_List.Nblk];
       for ( i_blk = 0; i_blk <= Blk_List.Nblk-1 ; ++i_blk ) {
          Blk_List.message_reschange_bottomwestcorner_recbuf[i_blk] = new double[1];
          Blk_List.message_reschange_bottomwestcorner_recbuf[i_blk][0] = ZERO;
       } /* endfor */
    } /* endif */

    if (Blk_List.Nblk > 0 &&
        Blk_List.message_reschange_bottomnorthwestcorner_sendbuf == NULL) {
       Blk_List.message_reschange_bottomnorthwestcorner_sendbuf = new double*[Blk_List.Nblk];
       for ( i_blk = 0; i_blk <= Blk_List.Nblk-1 ; ++i_blk ) {
          Blk_List.message_reschange_bottomnorthwestcorner_sendbuf[i_blk] = new double[1];
          Blk_List.message_reschange_bottomnorthwestcorner_sendbuf[i_blk][0] = ZERO;
       } /* endfor */
    } /* endif */
    if (Blk_List.Nblk > 0 &&
        Blk_List.message_reschange_bottomnorthwestcorner_recbuf == NULL) {
       Blk_List.message_reschange_bottomnorthwestcorner_recbuf = new double*[Blk_List.Nblk];
       for ( i_blk = 0; i_blk <= Blk_List.Nblk-1 ; ++i_blk ) {
          Blk_List.message_reschange_bottomnorthwestcorner_recbuf[i_blk] = new double[1];
          Blk_List.message_reschange_bottomnorthwestcorner_recbuf[i_blk][0] = ZERO;
       } /* endfor */
    } /* endif */
    if (Blk_List.Nblk > 0 &&
        Blk_List.message_reschange_bottomnortheastcorner_sendbuf == NULL) {
       Blk_List.message_reschange_bottomnortheastcorner_sendbuf = new double*[Blk_List.Nblk];
       for ( i_blk = 0; i_blk <= Blk_List.Nblk-1 ; ++i_blk ) {
          Blk_List.message_reschange_bottomnortheastcorner_sendbuf[i_blk] = new double[1];
          Blk_List.message_reschange_bottomnortheastcorner_sendbuf[i_blk][0] = ZERO;
       } /* endfor */
    } /* endif */
    if (Blk_List.Nblk > 0 &&
        Blk_List.message_reschange_bottomnortheastcorner_recbuf == NULL) {
       Blk_List.message_reschange_bottomnortheastcorner_recbuf = new double*[Blk_List.Nblk];
       for ( i_blk = 0; i_blk <= Blk_List.Nblk-1 ; ++i_blk ) {
          Blk_List.message_reschange_bottomnortheastcorner_recbuf[i_blk] = new double[1];
          Blk_List.message_reschange_bottomnortheastcorner_recbuf[i_blk][0] = ZERO;
       } /* endfor */
    } /* endif */
    if (Blk_List.Nblk > 0 &&
        Blk_List.message_reschange_bottomsoutheastcorner_sendbuf == NULL) {
       Blk_List.message_reschange_bottomsoutheastcorner_sendbuf = new double*[Blk_List.Nblk];
       for ( i_blk = 0; i_blk <= Blk_List.Nblk-1 ; ++i_blk ) {
          Blk_List.message_reschange_bottomsoutheastcorner_sendbuf[i_blk] = new double[1];
          Blk_List.message_reschange_bottomsoutheastcorner_sendbuf[i_blk][0] = ZERO;
       } /* endfor */
    } /* endif */
    if (Blk_List.Nblk > 0 &&
        Blk_List.message_reschange_bottomsoutheastcorner_recbuf == NULL) {
       Blk_List.message_reschange_bottomsoutheastcorner_recbuf = new double*[Blk_List.Nblk];
       for ( i_blk = 0; i_blk <= Blk_List.Nblk-1 ; ++i_blk ) {
          Blk_List.message_reschange_bottomsoutheastcorner_recbuf[i_blk] = new double[1];
          Blk_List.message_reschange_bottomsoutheastcorner_recbuf[i_blk][0] = ZERO;
       } /* endfor */
    } /* endif */
    if (Blk_List.Nblk > 0 &&
        Blk_List.message_reschange_bottomsouthwestcorner_sendbuf == NULL) {
       Blk_List.message_reschange_bottomsouthwestcorner_sendbuf = new double*[Blk_List.Nblk];
       for ( i_blk = 0; i_blk <= Blk_List.Nblk-1 ; ++i_blk ) {
          Blk_List.message_reschange_bottomsouthwestcorner_sendbuf[i_blk] = new double[1];
          Blk_List.message_reschange_bottomsouthwestcorner_sendbuf[i_blk][0] = ZERO;
       } /* endfor */
    } /* endif */
    if (Blk_List.Nblk > 0 &&
        Blk_List.message_reschange_bottomsouthwestcorner_recbuf == NULL) {
       Blk_List.message_reschange_bottomsouthwestcorner_recbuf = new double*[Blk_List.Nblk];
       for ( i_blk = 0; i_blk <= Blk_List.Nblk-1 ; ++i_blk ) {
          Blk_List.message_reschange_bottomsouthwestcorner_recbuf[i_blk] = new double[1];
          Blk_List.message_reschange_bottomsouthwestcorner_recbuf[i_blk][0] = ZERO;
       } /* endfor */
    } /* endif */
 
   /* Reallocate memory for the send and receive buffers of each block. */
    
    for ( i_blk = 0; i_blk <= Blk_List.Nblk-1; ++i_blk ) {
       ////////////////////////////////////////////////////////////////
       // Reallocate memory for Top face send and receive buffers. //
       ////////////////////////////////////////////////////////////////
       if (Blk_List.message_reschange_topface_sendbuf[i_blk] != NULL) {
          for (j_neigh = 0; j_neigh < 2 ; ++j_neigh) {
             delete []Blk_List.message_reschange_topface_sendbuf[i_blk][j_neigh];
             Blk_List.message_reschange_topface_sendbuf[i_blk][j_neigh] = NULL;
          } /* endfor */
       } /* endif */
       if (Blk_List.message_reschange_topface_recbuf[i_blk] != NULL) {
          for (j_neigh = 0; j_neigh < 2 ; ++j_neigh) {
             delete []Blk_List.message_reschange_topface_recbuf[i_blk][j_neigh];
             Blk_List.message_reschange_topface_recbuf[i_blk][j_neigh] = NULL;
	  } /* endfor */
       } /* endif */

       // Top neighbour is at higher resolution!
       if (Blk_List.Block[i_blk].used && 
           (Blk_List.Block[i_blk].nT == 2) &&
           (Blk_List.Block[i_blk].info.level < Blk_List.Block[i_blk].infoT[0].level)) {
          buffer_size = Blk_List.Block[i_blk].info.dimen.ghost*
                        (abs(Blk_List.Block[i_blk].info.dimen.i)/2+1)*
                        Number_of_Solution_Variables+
                        2*Blk_List.Block[i_blk].info.dimen.ghost;
          buffer_size_neighbour = Blk_List.Block[i_blk].infoT[0].dimen.ghost*
                                  (abs(Blk_List.Block[i_blk].infoT[0].dimen.i)+
                                  Blk_List.Block[i_blk].infoT[0].dimen.ghost+1)*
                                  Number_of_Solution_Variables+
                                  2*Blk_List.Block[i_blk].infoT[0].dimen.ghost;
          // Send buffers. 
          for (j_neigh = 0; j_neigh <= Blk_List.Block[i_blk].nT-1 ; ++j_neigh) {
             Blk_List.message_reschange_topface_sendbuf[i_blk][j_neigh] = 
                new double[buffer_size_neighbour];
             for (l = 0; l <= buffer_size_neighbour-1; ++l) 
	        Blk_List.message_reschange_topface_sendbuf[i_blk][j_neigh][l] = ZERO;
	  } /* endfor */
          // Receive buffers.
          for (j_neigh = 0; j_neigh <= Blk_List.Block[i_blk].nT-1 ; ++j_neigh) {
             Blk_List.message_reschange_topface_recbuf[i_blk][j_neigh] = 
                new double[buffer_size];
             for (l = 0; l <= buffer_size-1; ++l) 
	        Blk_List.message_reschange_topface_recbuf[i_blk][j_neigh][l] = ZERO;
	  } /* endfor */

       // Top neighbour is at lower resolution!
       } else if (Blk_List.Block[i_blk].used && 
                  (Blk_List.Block[i_blk].nT == 1) &&
                  (Blk_List.Block[i_blk].info.level > Blk_List.Block[i_blk].infoT[0].level)) {
          buffer_size = Blk_List.Block[i_blk].info.dimen.ghost*
                        (abs(Blk_List.Block[i_blk].info.dimen.i)+
                        Blk_List.Block[i_blk].info.dimen.ghost+1)*
                        Number_of_Solution_Variables+
                        2*Blk_List.Block[i_blk].info.dimen.ghost;
          buffer_size_neighbour = Blk_List.Block[i_blk].infoT[0].dimen.ghost*
                                  (abs(Blk_List.Block[i_blk].infoT[0].dimen.i)/2+1)*
                                  Number_of_Solution_Variables+
                                  2*Blk_List.Block[i_blk].infoT[0].dimen.ghost;
          // Send buffers. 
          Blk_List.message_reschange_topface_sendbuf[i_blk][0] = 
             new double[buffer_size_neighbour];
          for (l = 0; l <= buffer_size_neighbour-1; ++l) 
	     Blk_List.message_reschange_topface_sendbuf[i_blk][0][l] = ZERO;
          Blk_List.message_reschange_topface_sendbuf[i_blk][1] = 
             new double[1];
          Blk_List.message_reschange_topface_sendbuf[i_blk][1][0] = ZERO;
          // Receive buffers.
          Blk_List.message_reschange_topface_recbuf[i_blk][0] = 
             new double[buffer_size];
          for (l = 0; l <= buffer_size-1; ++l) 
	     Blk_List.message_reschange_topface_recbuf[i_blk][0][l] = ZERO;
          Blk_List.message_reschange_topface_recbuf[i_blk][1] = 
             new double[1];
          Blk_List.message_reschange_topface_recbuf[i_blk][1][0] = ZERO;

       // No Top neighbour with resolution change!
       } else {
	  // Send and receive buffers.
          for (j_neigh = 0; j_neigh < 2 ; ++j_neigh) {
             Blk_List.message_reschange_topface_sendbuf[i_blk][j_neigh] = 
                new double[1];
             Blk_List.message_reschange_topface_sendbuf[i_blk][j_neigh][0] = ZERO;
             Blk_List.message_reschange_topface_recbuf[i_blk][j_neigh] = 
                new double[1];
             Blk_List.message_reschange_topface_recbuf[i_blk][j_neigh][0] = ZERO;
          } /* endfor */
       } /* endif */


       ////////////////////////////////////////////////////////////////
       // Reallocate memory for Bottom face send and receive buffers. //
       ////////////////////////////////////////////////////////////////
       if (Blk_List.message_reschange_bottomface_sendbuf[i_blk] != NULL) {
          for (j_neigh = 0; j_neigh < 2 ; ++j_neigh) {
             delete []Blk_List.message_reschange_bottomface_sendbuf[i_blk][j_neigh];
             Blk_List.message_reschange_bottomface_sendbuf[i_blk][j_neigh] = NULL;
          } /* endfor */
       } /* endif */
       if (Blk_List.message_reschange_bottomface_recbuf[i_blk] != NULL) {
          for (j_neigh = 0; j_neigh < 2 ; ++j_neigh) {
             delete []Blk_List.message_reschange_bottomface_recbuf[i_blk][j_neigh];
             Blk_List.message_reschange_bottomface_recbuf[i_blk][j_neigh] = NULL;
	  } /* endfor */
       } /* endif */

       // Bottom neighbour is at higher resolution!
       if (Blk_List.Block[i_blk].used && 
           (Blk_List.Block[i_blk].nB == 2) &&
           (Blk_List.Block[i_blk].info.level < Blk_List.Block[i_blk].infoT[0].level)) {
          buffer_size = Blk_List.Block[i_blk].info.dimen.ghost*
                        (abs(Blk_List.Block[i_blk].info.dimen.i)/2+1)*
                        Number_of_Solution_Variables+
                        2*Blk_List.Block[i_blk].info.dimen.ghost;
          buffer_size_neighbour = Blk_List.Block[i_blk].infoT[0].dimen.ghost*
                                  (abs(Blk_List.Block[i_blk].infoT[0].dimen.i)+
                                  Blk_List.Block[i_blk].infoT[0].dimen.ghost+1)*
                                  Number_of_Solution_Variables+
                                  2*Blk_List.Block[i_blk].infoT[0].dimen.ghost;
          // Send buffers. 
          for (j_neigh = 0; j_neigh <= Blk_List.Block[i_blk].nB-1 ; ++j_neigh) {
             Blk_List.message_reschange_bottomface_sendbuf[i_blk][j_neigh] = 
                new double[buffer_size_neighbour];
             for (l = 0; l <= buffer_size_neighbour-1; ++l) 
	        Blk_List.message_reschange_bottomface_sendbuf[i_blk][j_neigh][l] = ZERO;
	  } /* endfor */
          // Receive buffers.
          for (j_neigh = 0; j_neigh <= Blk_List.Block[i_blk].nB-1 ; ++j_neigh) {
             Blk_List.message_reschange_bottomface_recbuf[i_blk][j_neigh] = 
                new double[buffer_size];
             for (l = 0; l <= buffer_size-1; ++l) 
	        Blk_List.message_reschange_bottomface_recbuf[i_blk][j_neigh][l] = ZERO;
	  } /* endfor */

       // Bottom neighbour is at lower resolution!
       } else if (Blk_List.Block[i_blk].used && 
                  (Blk_List.Block[i_blk].nB == 1) &&
                  (Blk_List.Block[i_blk].info.level > Blk_List.Block[i_blk].infoT[0].level)) {
          buffer_size = Blk_List.Block[i_blk].info.dimen.ghost*
                        (abs(Blk_List.Block[i_blk].info.dimen.i)+
                        Blk_List.Block[i_blk].info.dimen.ghost+1)*
                        Number_of_Solution_Variables+
                        2*Blk_List.Block[i_blk].info.dimen.ghost;
          buffer_size_neighbour = Blk_List.Block[i_blk].infoT[0].dimen.ghost*
                                  (abs(Blk_List.Block[i_blk].infoT[0].dimen.i)/2+1)*
                                  Number_of_Solution_Variables+
                                  2*Blk_List.Block[i_blk].infoT[0].dimen.ghost;
          // Send buffers. 
          Blk_List.message_reschange_bottomface_sendbuf[i_blk][0] = 
             new double[buffer_size_neighbour];
          for (l = 0; l <= buffer_size_neighbour-1; ++l) 
	     Blk_List.message_reschange_bottomface_sendbuf[i_blk][0][l] = ZERO;
          Blk_List.message_reschange_bottomface_sendbuf[i_blk][1] = 
             new double[1];
          Blk_List.message_reschange_bottomface_sendbuf[i_blk][1][0] = ZERO;
          // Receive buffers.
          Blk_List.message_reschange_bottomface_recbuf[i_blk][0] = 
             new double[buffer_size];
          for (l = 0; l <= buffer_size-1; ++l) 
	     Blk_List.message_reschange_bottomface_recbuf[i_blk][0][l] = ZERO;
          Blk_List.message_reschange_bottomface_recbuf[i_blk][1] = 
             new double[1];
          Blk_List.message_reschange_bottomface_recbuf[i_blk][1][0] = ZERO;

       // No Bottom neighbour with resolution change!
       } else {
	  // Send and receive buffers.
          for (j_neigh = 0; j_neigh < 2 ; ++j_neigh) {
             Blk_List.message_reschange_bottomface_sendbuf[i_blk][j_neigh] = 
                new double[1];
             Blk_List.message_reschange_bottomface_sendbuf[i_blk][j_neigh][0] = ZERO;
             Blk_List.message_reschange_bottomface_recbuf[i_blk][j_neigh] = 
                new double[1];
             Blk_List.message_reschange_bottomface_recbuf[i_blk][j_neigh][0] = ZERO;
          } /* endfor */
       } /* endif */



       ////////////////////////////////////////////////////////////////
       // Reallocate memory for North face send and receive buffers. //
       ////////////////////////////////////////////////////////////////
       if (Blk_List.message_reschange_northface_sendbuf[i_blk] != NULL) {
          for (j_neigh = 0; j_neigh < 2 ; ++j_neigh) {
             delete []Blk_List.message_reschange_northface_sendbuf[i_blk][j_neigh];
             Blk_List.message_reschange_northface_sendbuf[i_blk][j_neigh] = NULL;
          } /* endfor */
       } /* endif */
       if (Blk_List.message_reschange_northface_recbuf[i_blk] != NULL) {
          for (j_neigh = 0; j_neigh < 2 ; ++j_neigh) {
             delete []Blk_List.message_reschange_northface_recbuf[i_blk][j_neigh];
             Blk_List.message_reschange_northface_recbuf[i_blk][j_neigh] = NULL;
	  } /* endfor */
       } /* endif */

       // North neighbour is at higher resolution!
       if (Blk_List.Block[i_blk].used && 
           (Blk_List.Block[i_blk].nN == 2) &&
           (Blk_List.Block[i_blk].info.level < Blk_List.Block[i_blk].infoN[0].level)) {
          buffer_size = Blk_List.Block[i_blk].info.dimen.ghost*
                        (abs(Blk_List.Block[i_blk].info.dimen.i)/2+1)*
                        Number_of_Solution_Variables+
                        2*Blk_List.Block[i_blk].info.dimen.ghost;
          buffer_size_neighbour = Blk_List.Block[i_blk].infoN[0].dimen.ghost*
                                  (abs(Blk_List.Block[i_blk].infoN[0].dimen.i)+
                                  Blk_List.Block[i_blk].infoN[0].dimen.ghost+1)*
                                  Number_of_Solution_Variables+
                                  2*Blk_List.Block[i_blk].infoN[0].dimen.ghost;
          // Send buffers. 
          for (j_neigh = 0; j_neigh <= Blk_List.Block[i_blk].nN-1 ; ++j_neigh) {
             Blk_List.message_reschange_northface_sendbuf[i_blk][j_neigh] = 
                new double[buffer_size_neighbour];
             for (l = 0; l <= buffer_size_neighbour-1; ++l) 
	        Blk_List.message_reschange_northface_sendbuf[i_blk][j_neigh][l] = ZERO;
	  } /* endfor */
          // Receive buffers.
          for (j_neigh = 0; j_neigh <= Blk_List.Block[i_blk].nN-1 ; ++j_neigh) {
             Blk_List.message_reschange_northface_recbuf[i_blk][j_neigh] = 
                new double[buffer_size];
             for (l = 0; l <= buffer_size-1; ++l) 
	        Blk_List.message_reschange_northface_recbuf[i_blk][j_neigh][l] = ZERO;
	  } /* endfor */

       // North neighbour is at lower resolution!
       } else if (Blk_List.Block[i_blk].used && 
                  (Blk_List.Block[i_blk].nN == 1) &&
                  (Blk_List.Block[i_blk].info.level > Blk_List.Block[i_blk].infoN[0].level)) {
          buffer_size = Blk_List.Block[i_blk].info.dimen.ghost*
                        (abs(Blk_List.Block[i_blk].info.dimen.i)+
                        Blk_List.Block[i_blk].info.dimen.ghost+1)*
                        Number_of_Solution_Variables+
                        2*Blk_List.Block[i_blk].info.dimen.ghost;
          buffer_size_neighbour = Blk_List.Block[i_blk].infoN[0].dimen.ghost*
                                  (abs(Blk_List.Block[i_blk].infoN[0].dimen.i)/2+1)*
                                  Number_of_Solution_Variables+
                                  2*Blk_List.Block[i_blk].infoN[0].dimen.ghost;
          // Send buffers. 
          Blk_List.message_reschange_northface_sendbuf[i_blk][0] = 
             new double[buffer_size_neighbour];
          for (l = 0; l <= buffer_size_neighbour-1; ++l) 
	     Blk_List.message_reschange_northface_sendbuf[i_blk][0][l] = ZERO;
          Blk_List.message_reschange_northface_sendbuf[i_blk][1] = 
             new double[1];
          Blk_List.message_reschange_northface_sendbuf[i_blk][1][0] = ZERO;
          // Receive buffers.
          Blk_List.message_reschange_northface_recbuf[i_blk][0] = 
             new double[buffer_size];
          for (l = 0; l <= buffer_size-1; ++l) 
	     Blk_List.message_reschange_northface_recbuf[i_blk][0][l] = ZERO;
          Blk_List.message_reschange_northface_recbuf[i_blk][1] = 
             new double[1];
          Blk_List.message_reschange_northface_recbuf[i_blk][1][0] = ZERO;

       // No North neighbour with resolution change!
       } else {
	  // Send and receive buffers.
          for (j_neigh = 0; j_neigh < 2 ; ++j_neigh) {
             Blk_List.message_reschange_northface_sendbuf[i_blk][j_neigh] = 
                new double[1];
             Blk_List.message_reschange_northface_sendbuf[i_blk][j_neigh][0] = ZERO;
             Blk_List.message_reschange_northface_recbuf[i_blk][j_neigh] = 
                new double[1];
             Blk_List.message_reschange_northface_recbuf[i_blk][j_neigh][0] = ZERO;
          } /* endfor */
       } /* endif */

       ////////////////////////////////////////////////////////////////
       // Reallocate memory for South face send and receive buffers. //
       ////////////////////////////////////////////////////////////////
       if (Blk_List.message_reschange_southface_sendbuf[i_blk] != NULL) {
          for (j_neigh = 0; j_neigh < 2 ; ++j_neigh) {
             delete []Blk_List.message_reschange_southface_sendbuf[i_blk][j_neigh];
             Blk_List.message_reschange_southface_sendbuf[i_blk][j_neigh] = NULL;
          } /* endfor */
       } /* endif */
       if (Blk_List.message_reschange_southface_recbuf[i_blk] != NULL) {
          for (j_neigh = 0; j_neigh < 2 ; ++j_neigh) {
             delete []Blk_List.message_reschange_southface_recbuf[i_blk][j_neigh];
             Blk_List.message_reschange_southface_recbuf[i_blk][j_neigh] = NULL;
	  } /* endfor */
       } /* endif */

       // South neighbour is at higher resolution!
       if (Blk_List.Block[i_blk].used && 
           (Blk_List.Block[i_blk].nS == 2) &&
           (Blk_List.Block[i_blk].info.level < Blk_List.Block[i_blk].infoS[0].level)) {
          buffer_size = Blk_List.Block[i_blk].info.dimen.ghost*
                        (abs(Blk_List.Block[i_blk].info.dimen.i)/2+1)*
                        Number_of_Solution_Variables+
                        2*Blk_List.Block[i_blk].info.dimen.ghost;
          buffer_size_neighbour = Blk_List.Block[i_blk].infoS[0].dimen.ghost*
                                  (abs(Blk_List.Block[i_blk].infoS[0].dimen.i)+
                                  Blk_List.Block[i_blk].infoS[0].dimen.ghost+1)*
                                  Number_of_Solution_Variables+
                                  2*Blk_List.Block[i_blk].infoS[0].dimen.ghost;
          // Send buffers. 
          for (j_neigh = 0; j_neigh <= Blk_List.Block[i_blk].nS-1 ; ++j_neigh) {
             Blk_List.message_reschange_southface_sendbuf[i_blk][j_neigh] = 
                new double[buffer_size_neighbour];
             for (l = 0; l <= buffer_size_neighbour-1; ++l) 
	        Blk_List.message_reschange_southface_sendbuf[i_blk][j_neigh][l] = ZERO;
	  } /* endfor */
          // Receive buffers.
          for (j_neigh = 0; j_neigh <= Blk_List.Block[i_blk].nS-1 ; ++j_neigh) {
             Blk_List.message_reschange_southface_recbuf[i_blk][j_neigh] = 
                new double[buffer_size];
             for (l = 0; l <= buffer_size-1; ++l) 
	        Blk_List.message_reschange_southface_recbuf[i_blk][j_neigh][l] = ZERO;
	  } /* endfor */

       // South neighbour is at lower resolution!
       } else if (Blk_List.Block[i_blk].used && 
                  (Blk_List.Block[i_blk].nS == 1) &&
                  (Blk_List.Block[i_blk].info.level > Blk_List.Block[i_blk].infoS[0].level)) {
          buffer_size = Blk_List.Block[i_blk].info.dimen.ghost*
                        (abs(Blk_List.Block[i_blk].info.dimen.i)+
                        Blk_List.Block[i_blk].info.dimen.ghost+1)*
                        Number_of_Solution_Variables+
                        2*Blk_List.Block[i_blk].info.dimen.ghost;
          buffer_size_neighbour = Blk_List.Block[i_blk].infoS[0].dimen.ghost*
                                  (abs(Blk_List.Block[i_blk].infoS[0].dimen.i)/2+1)*
                                  Number_of_Solution_Variables+
                                  2*Blk_List.Block[i_blk].infoS[0].dimen.ghost;
          // Send buffers. 
          Blk_List.message_reschange_southface_sendbuf[i_blk][0] = 
             new double[buffer_size_neighbour];
          for (l = 0; l <= buffer_size_neighbour-1; ++l) 
	     Blk_List.message_reschange_southface_sendbuf[i_blk][0][l] = ZERO;
          Blk_List.message_reschange_southface_sendbuf[i_blk][1] = 
             new double[1];
          Blk_List.message_reschange_southface_sendbuf[i_blk][1][0] = ZERO;
          // Receive buffers.
          Blk_List.message_reschange_southface_recbuf[i_blk][0] = 
             new double[buffer_size];
          for (l = 0; l <= buffer_size-1; ++l) 
	     Blk_List.message_reschange_southface_recbuf[i_blk][0][l] = ZERO;
          Blk_List.message_reschange_southface_recbuf[i_blk][1] = 
             new double[1];
          Blk_List.message_reschange_southface_recbuf[i_blk][1][0] = ZERO;

       // No South neighbour with resolution change!
       } else {
	  // Send and receive buffers.
          for (j_neigh = 0; j_neigh < 2 ; ++j_neigh) {
             Blk_List.message_reschange_southface_sendbuf[i_blk][j_neigh] = 
                new double[1];
             Blk_List.message_reschange_southface_sendbuf[i_blk][j_neigh][0] = ZERO;
             Blk_List.message_reschange_southface_recbuf[i_blk][j_neigh] = 
                new double[1];
             Blk_List.message_reschange_southface_recbuf[i_blk][j_neigh][0] = ZERO;
          } /* endfor */
       } /* endif */

       ///////////////////////////////////////////////////////////////
       // Reallocate memory for East face send and receive buffers. //
       ///////////////////////////////////////////////////////////////
       if (Blk_List.message_reschange_eastface_sendbuf[i_blk] != NULL) {
          for (j_neigh = 0; j_neigh < 2 ; ++j_neigh) {
             delete []Blk_List.message_reschange_eastface_sendbuf[i_blk][j_neigh];
             Blk_List.message_reschange_eastface_sendbuf[i_blk][j_neigh] = NULL;
          } /* endfor */
       } /* endif */
       if (Blk_List.message_reschange_eastface_recbuf[i_blk] != NULL) {
          for (j_neigh = 0; j_neigh < 2 ; ++j_neigh) {
             delete []Blk_List.message_reschange_eastface_recbuf[i_blk][j_neigh];
             Blk_List.message_reschange_eastface_recbuf[i_blk][j_neigh] = NULL;
	  } /* endfor */
       } /* endif */

       // East neighbour is at higher resolution!
       if (Blk_List.Block[i_blk].used && 
           (Blk_List.Block[i_blk].nE == 2) &&
           (Blk_List.Block[i_blk].info.level < Blk_List.Block[i_blk].infoE[0].level)) {
          buffer_size = Blk_List.Block[i_blk].info.dimen.ghost*
                        (abs(Blk_List.Block[i_blk].info.dimen.j)/2+1)*
                        Number_of_Solution_Variables+
                        2*Blk_List.Block[i_blk].info.dimen.ghost;
          buffer_size_neighbour = Blk_List.Block[i_blk].infoE[0].dimen.ghost*
                                  (abs(Blk_List.Block[i_blk].infoE[0].dimen.j)+
                                  Blk_List.Block[i_blk].infoE[0].dimen.ghost+1)*
                                  Number_of_Solution_Variables+
                                  2*Blk_List.Block[i_blk].infoE[0].dimen.ghost;
          // Send buffers. 
          for (j_neigh = 0; j_neigh <= Blk_List.Block[i_blk].nE-1 ; ++j_neigh) {
             Blk_List.message_reschange_eastface_sendbuf[i_blk][j_neigh] = 
                new double[buffer_size_neighbour];
             for (l = 0; l <= buffer_size_neighbour-1; ++l) 
	        Blk_List.message_reschange_eastface_sendbuf[i_blk][j_neigh][l] = ZERO;
	  } /* endfor */
          // Receive buffers.
          for (j_neigh = 0; j_neigh <= Blk_List.Block[i_blk].nE-1 ; ++j_neigh) {
             Blk_List.message_reschange_eastface_recbuf[i_blk][j_neigh] = 
                new double[buffer_size];
             for (l = 0; l <= buffer_size-1; ++l) 
	        Blk_List.message_reschange_eastface_recbuf[i_blk][j_neigh][l] = ZERO;
	  } /* endfor */

       // East neighbour is at lower resolution!
       } else if (Blk_List.Block[i_blk].used && 
                  (Blk_List.Block[i_blk].nE == 1) &&
                  (Blk_List.Block[i_blk].info.level > Blk_List.Block[i_blk].infoE[0].level)) {
          buffer_size = Blk_List.Block[i_blk].info.dimen.ghost*
                        (abs(Blk_List.Block[i_blk].info.dimen.j)+
                        Blk_List.Block[i_blk].info.dimen.ghost+1)*
                        Number_of_Solution_Variables+
                        2*Blk_List.Block[i_blk].info.dimen.ghost;
          buffer_size_neighbour = Blk_List.Block[i_blk].infoE[0].dimen.ghost*
                                  (abs(Blk_List.Block[i_blk].infoE[0].dimen.j)/2+1)*
                                  Number_of_Solution_Variables+
                                  2*Blk_List.Block[i_blk].infoE[0].dimen.ghost;
          // Send buffers. 
          Blk_List.message_reschange_eastface_sendbuf[i_blk][0] = 
             new double[buffer_size_neighbour];
          for (l = 0; l <= buffer_size_neighbour-1; ++l) 
	     Blk_List.message_reschange_eastface_sendbuf[i_blk][0][l] = ZERO;
          Blk_List.message_reschange_eastface_sendbuf[i_blk][1] = 
             new double[1];
          Blk_List.message_reschange_eastface_sendbuf[i_blk][1][0] = ZERO;
          // Receive buffers.
          Blk_List.message_reschange_eastface_recbuf[i_blk][0] = 
             new double[buffer_size];
          for (l = 0; l <= buffer_size-1; ++l) 
	     Blk_List.message_reschange_eastface_recbuf[i_blk][0][l] = ZERO;
          Blk_List.message_reschange_eastface_recbuf[i_blk][1] = 
             new double[1];
          Blk_List.message_reschange_eastface_recbuf[i_blk][1][0] = ZERO;

       // No East neighbour with resolution change!
       } else {
	  // Send and receive buffers.
          for (j_neigh = 0; j_neigh < 2 ; ++j_neigh) {
             Blk_List.message_reschange_eastface_sendbuf[i_blk][j_neigh] = 
                new double[1];
             Blk_List.message_reschange_eastface_sendbuf[i_blk][j_neigh][0] = ZERO;
             Blk_List.message_reschange_eastface_recbuf[i_blk][j_neigh] = 
                new double[1];
             Blk_List.message_reschange_eastface_recbuf[i_blk][j_neigh][0] = ZERO;
          } /* endfor */
       } /* endif */

       ///////////////////////////////////////////////////////////////
       // Reallocate memory for West face send and receive buffers. //
       ///////////////////////////////////////////////////////////////
       if (Blk_List.message_reschange_westface_sendbuf[i_blk] != NULL) {
          for (j_neigh = 0; j_neigh < 2 ; ++j_neigh) {
             delete []Blk_List.message_reschange_westface_sendbuf[i_blk][j_neigh];
             Blk_List.message_reschange_westface_sendbuf[i_blk][j_neigh] = NULL;
          } /* endfor */
       } /* endif */
       if (Blk_List.message_reschange_westface_recbuf[i_blk] != NULL) {
          for (j_neigh = 0; j_neigh < 2 ; ++j_neigh) {
             delete []Blk_List.message_reschange_westface_recbuf[i_blk][j_neigh];
             Blk_List.message_reschange_westface_recbuf[i_blk][j_neigh] = NULL;
	  } /* endfor */
       } /* endif */

       // West neighbour is at higher resolution!
       if (Blk_List.Block[i_blk].used && 
           (Blk_List.Block[i_blk].nW == 2) &&
           (Blk_List.Block[i_blk].info.level < Blk_List.Block[i_blk].infoW[0].level)) {
          buffer_size = Blk_List.Block[i_blk].info.dimen.ghost*
                        (abs(Blk_List.Block[i_blk].info.dimen.j)/2+1)*
                        Number_of_Solution_Variables+
                        2*Blk_List.Block[i_blk].info.dimen.ghost;
          buffer_size_neighbour = Blk_List.Block[i_blk].infoW[0].dimen.ghost*
                                  (abs(Blk_List.Block[i_blk].infoW[0].dimen.j)+
                                  Blk_List.Block[i_blk].infoW[0].dimen.ghost+1)*
                                  Number_of_Solution_Variables+
                                  2*Blk_List.Block[i_blk].infoW[0].dimen.ghost;
          // Send buffers. 
          for (j_neigh = 0; j_neigh <= Blk_List.Block[i_blk].nW-1 ; ++j_neigh) {
             Blk_List.message_reschange_westface_sendbuf[i_blk][j_neigh] = 
                new double[buffer_size_neighbour];
             for (l = 0; l <= buffer_size_neighbour-1; ++l) 
	        Blk_List.message_reschange_westface_sendbuf[i_blk][j_neigh][l] = ZERO;
	  } /* endfor */
          // Receive buffers.
          for (j_neigh = 0; j_neigh <= Blk_List.Block[i_blk].nW-1 ; ++j_neigh) {
             Blk_List.message_reschange_westface_recbuf[i_blk][j_neigh] = 
                new double[buffer_size];
             for (l = 0; l <= buffer_size-1; ++l) 
	        Blk_List.message_reschange_westface_recbuf[i_blk][j_neigh][l] = ZERO;
	  } /* endfor */

       // West neighbour is at lower resolution!
       } else if (Blk_List.Block[i_blk].used && 
                  (Blk_List.Block[i_blk].nW == 1) &&
                  (Blk_List.Block[i_blk].info.level > Blk_List.Block[i_blk].infoW[0].level)) {
          buffer_size = Blk_List.Block[i_blk].info.dimen.ghost*
                        (abs(Blk_List.Block[i_blk].info.dimen.j)+
                        Blk_List.Block[i_blk].info.dimen.ghost+1)*
                        Number_of_Solution_Variables+
                        2*Blk_List.Block[i_blk].info.dimen.ghost;
          buffer_size_neighbour = Blk_List.Block[i_blk].infoW[0].dimen.ghost*
                                  (abs(Blk_List.Block[i_blk].infoW[0].dimen.j)/2+1)*
                                  Number_of_Solution_Variables+
                                  2*Blk_List.Block[i_blk].infoW[0].dimen.ghost;
          // Send buffers. 
          Blk_List.message_reschange_westface_sendbuf[i_blk][0] = 
             new double[buffer_size_neighbour];
          for (l = 0; l <= buffer_size_neighbour-1; ++l) 
	     Blk_List.message_reschange_westface_sendbuf[i_blk][0][l] = ZERO;
          Blk_List.message_reschange_westface_sendbuf[i_blk][1] = 
             new double[1];
          Blk_List.message_reschange_westface_sendbuf[i_blk][1][0] = ZERO;
          // Receive buffers.
          Blk_List.message_reschange_westface_recbuf[i_blk][0] = 
             new double[buffer_size];
          for (l = 0; l <= buffer_size-1; ++l) 
	     Blk_List.message_reschange_westface_recbuf[i_blk][0][l] = ZERO;
          Blk_List.message_reschange_westface_recbuf[i_blk][1] = 
             new double[1];
          Blk_List.message_reschange_westface_recbuf[i_blk][1][0] = ZERO;

       // No West neighbour with resolution change!
       } else {
	  // Send and receive buffers.
          for (j_neigh = 0; j_neigh < 2 ; ++j_neigh) {
             Blk_List.message_reschange_westface_sendbuf[i_blk][j_neigh] = 
                new double[1];
             Blk_List.message_reschange_westface_sendbuf[i_blk][j_neigh][0] = ZERO;
             Blk_List.message_reschange_westface_recbuf[i_blk][j_neigh] = 
                new double[1];
             Blk_List.message_reschange_westface_recbuf[i_blk][j_neigh][0] = ZERO;
          } /* endfor */
       } /* endif */

       ///////////////////////////////////////////////////////////////////////
       // Reallocate memory for North West corner send and receive buffers. //
       ///////////////////////////////////////////////////////////////////////
       if (Blk_List.message_reschange_northwestcorner_sendbuf[i_blk] != NULL) {
          delete []Blk_List.message_reschange_northwestcorner_sendbuf[i_blk];
          Blk_List.message_reschange_northwestcorner_sendbuf[i_blk] = NULL;
       } /* endif */
       if (Blk_List.message_reschange_northwestcorner_recbuf[i_blk] != NULL) {
          delete []Blk_List.message_reschange_northwestcorner_recbuf[i_blk];
          Blk_List.message_reschange_northwestcorner_recbuf[i_blk] = NULL;
       } /* endif */

       // North West neighbour is at higher or lower resolution!
       if (Blk_List.Block[i_blk].used && 
           (Blk_List.Block[i_blk].nNW == 1) &&
           ((Blk_List.Block[i_blk].info.level < Blk_List.Block[i_blk].infoNW[0].level) ||
            (Blk_List.Block[i_blk].info.level > Blk_List.Block[i_blk].infoNW[0].level))) {
          buffer_size = sqr(Blk_List.Block[i_blk].info.dimen.ghost)*
                        Number_of_Solution_Variables;
          buffer_size_neighbour = sqr(Blk_List.Block[i_blk].infoNW[0].dimen.ghost)*
                                  Number_of_Solution_Variables;
          // Send buffer.
          Blk_List.message_reschange_northwestcorner_sendbuf[i_blk] = 
             new double[buffer_size_neighbour];
          for (l = 0; l <= buffer_size_neighbour-1; ++l) 
	     Blk_List.message_reschange_northwestcorner_sendbuf[i_blk][l] = ZERO;
          // Receive buffer.
          Blk_List.message_reschange_northwestcorner_recbuf[i_blk] = 
             new double[buffer_size];
          for (l = 0; l <= buffer_size-1; ++l) 
	     Blk_List.message_reschange_northwestcorner_recbuf[i_blk][l] = ZERO;

       // No North West neighbour with resolution change!
       } else {
          Blk_List.message_reschange_northwestcorner_sendbuf[i_blk] = 
             new double[1];
          Blk_List.message_reschange_northwestcorner_sendbuf[i_blk][0] = ZERO; 
          Blk_List.message_reschange_northwestcorner_recbuf[i_blk] = 
             new double[1];
          Blk_List.message_reschange_northwestcorner_recbuf[i_blk][0] = ZERO;
       } /* endif */

       ///////////////////////////////////////////////////////////////////////
       // Reallocate memory for North East corner send and receive buffers. //
       ///////////////////////////////////////////////////////////////////////
       if (Blk_List.message_reschange_northeastcorner_sendbuf[i_blk] != NULL) {
          delete []Blk_List.message_reschange_northeastcorner_sendbuf[i_blk];
          Blk_List.message_reschange_northeastcorner_sendbuf[i_blk] = NULL;
       } /* endif */
       if (Blk_List.message_reschange_northeastcorner_recbuf[i_blk] != NULL) {
          delete []Blk_List.message_reschange_northeastcorner_recbuf[i_blk];
          Blk_List.message_reschange_northeastcorner_recbuf[i_blk] = NULL;
       } /* endif */

       // North East neighbour is at higher or lower resolution!
       if (Blk_List.Block[i_blk].used && 
           (Blk_List.Block[i_blk].nNE == 1) &&
           ((Blk_List.Block[i_blk].info.level < Blk_List.Block[i_blk].infoNE[0].level) ||
            (Blk_List.Block[i_blk].info.level > Blk_List.Block[i_blk].infoNE[0].level))) {
          buffer_size = sqr(Blk_List.Block[i_blk].info.dimen.ghost)*
                        Number_of_Solution_Variables;
          buffer_size_neighbour = sqr(Blk_List.Block[i_blk].infoNE[0].dimen.ghost)*
                                  Number_of_Solution_Variables;
          // Send buffer.
          Blk_List.message_reschange_northeastcorner_sendbuf[i_blk] = 
             new double[buffer_size_neighbour];
          for (l = 0; l <= buffer_size_neighbour-1; ++l) 
	     Blk_List.message_reschange_northeastcorner_sendbuf[i_blk][l] = ZERO;
          // Receive buffer.
          Blk_List.message_reschange_northeastcorner_recbuf[i_blk] = 
             new double[buffer_size];
          for (l = 0; l <= buffer_size-1; ++l) 
	     Blk_List.message_reschange_northeastcorner_recbuf[i_blk][l] = ZERO;

       // No North East neighbour with resolution change!
       } else {
          Blk_List.message_reschange_northeastcorner_sendbuf[i_blk] = 
             new double[1];
          Blk_List.message_reschange_northeastcorner_sendbuf[i_blk][0] = ZERO; 
          Blk_List.message_reschange_northeastcorner_recbuf[i_blk] = 
             new double[1];
          Blk_List.message_reschange_northeastcorner_recbuf[i_blk][0] = ZERO;
       } /* endif */

       ///////////////////////////////////////////////////////////////////////
       // Reallocate memory for South East corner send and receive buffers. //
       ///////////////////////////////////////////////////////////////////////
       if (Blk_List.message_reschange_southeastcorner_sendbuf[i_blk] != NULL) {
          delete []Blk_List.message_reschange_southeastcorner_sendbuf[i_blk];
          Blk_List.message_reschange_southeastcorner_sendbuf[i_blk] = NULL;
       } /* endif */
       if (Blk_List.message_reschange_southeastcorner_recbuf[i_blk] != NULL) {
          delete []Blk_List.message_reschange_southeastcorner_recbuf[i_blk];
          Blk_List.message_reschange_southeastcorner_recbuf[i_blk] = NULL;
       } /* endif */

       // South East neighbour is at higher or lower resolution!
       if (Blk_List.Block[i_blk].used && 
           (Blk_List.Block[i_blk].nSE == 1) &&
           ((Blk_List.Block[i_blk].info.level < Blk_List.Block[i_blk].infoSE[0].level) ||
            (Blk_List.Block[i_blk].info.level > Blk_List.Block[i_blk].infoSE[0].level))) {
          buffer_size = sqr(Blk_List.Block[i_blk].info.dimen.ghost)*
                        Number_of_Solution_Variables;
          buffer_size_neighbour = sqr(Blk_List.Block[i_blk].infoSE[0].dimen.ghost)*
                                  Number_of_Solution_Variables;
          // Send buffer.
          Blk_List.message_reschange_southeastcorner_sendbuf[i_blk] = 
             new double[buffer_size_neighbour];
          for (l = 0; l <= buffer_size_neighbour-1; ++l) 
	     Blk_List.message_reschange_southeastcorner_sendbuf[i_blk][l] = ZERO;
          // Receive buffer.
          Blk_List.message_reschange_southeastcorner_recbuf[i_blk] = 
             new double[buffer_size];
          for (l = 0; l <= buffer_size-1; ++l) 
	     Blk_List.message_reschange_southeastcorner_recbuf[i_blk][l] = ZERO;

       // No South East neighbour with resolution change!
       } else {
          Blk_List.message_reschange_southeastcorner_sendbuf[i_blk] = 
             new double[1];
          Blk_List.message_reschange_southeastcorner_sendbuf[i_blk][0] = ZERO; 
          Blk_List.message_reschange_southeastcorner_recbuf[i_blk] = 
             new double[1];
          Blk_List.message_reschange_southeastcorner_recbuf[i_blk][0] = ZERO;
       } /* endif */

       ///////////////////////////////////////////////////////////////////////
       // Reallocate memory for South West corner send and receive buffers. //
       ///////////////////////////////////////////////////////////////////////
       if (Blk_List.message_reschange_southwestcorner_sendbuf[i_blk] != NULL) {
          delete []Blk_List.message_reschange_southwestcorner_sendbuf[i_blk];
          Blk_List.message_reschange_southwestcorner_sendbuf[i_blk] = NULL;
       } /* endif */
       if (Blk_List.message_reschange_southwestcorner_recbuf[i_blk] != NULL) {
          delete []Blk_List.message_reschange_southwestcorner_recbuf[i_blk];
          Blk_List.message_reschange_southwestcorner_recbuf[i_blk] = NULL;
       } /* endif */

       // South West neighbour is at higher or lower resolution!
       if (Blk_List.Block[i_blk].used && 
           (Blk_List.Block[i_blk].nSW == 1) &&
           ((Blk_List.Block[i_blk].info.level < Blk_List.Block[i_blk].infoSW[0].level) ||
            (Blk_List.Block[i_blk].info.level > Blk_List.Block[i_blk].infoSW[0].level))) {
          buffer_size = sqr(Blk_List.Block[i_blk].info.dimen.ghost)*
                        Number_of_Solution_Variables;
          buffer_size_neighbour = sqr(Blk_List.Block[i_blk].infoSW[0].dimen.ghost)*
                                  Number_of_Solution_Variables;
          // Send buffer.
          Blk_List.message_reschange_southwestcorner_sendbuf[i_blk] = 
             new double[buffer_size_neighbour];
          for (l = 0; l <= buffer_size_neighbour-1; ++l) 
	     Blk_List.message_reschange_southwestcorner_sendbuf[i_blk][l] = ZERO;
          // Receive buffer.
          Blk_List.message_reschange_southwestcorner_recbuf[i_blk] = 
             new double[buffer_size];
          for (l = 0; l <= buffer_size-1; ++l) 
	     Blk_List.message_reschange_southwestcorner_recbuf[i_blk][l] = ZERO;

       // No South West neighbour with resolution change!
       } else {
          Blk_List.message_reschange_southwestcorner_sendbuf[i_blk] = 
             new double[1];
          Blk_List.message_reschange_southwestcorner_sendbuf[i_blk][0] = ZERO; 
          Blk_List.message_reschange_southwestcorner_recbuf[i_blk] = 
             new double[1];
          Blk_List.message_reschange_southwestcorner_recbuf[i_blk][0] = ZERO;
       } /* endif */




    ///******************//

       ///////////////////////////////////////////////////////////////////////
       // Reallocate memory for Top North corner send and receive buffers. //
       ///////////////////////////////////////////////////////////////////////
       if (Blk_List.message_reschange_topnorthcorner_sendbuf[i_blk] != NULL) {
          delete []Blk_List.message_reschange_topnorthcorner_sendbuf[i_blk];
          Blk_List.message_reschange_topnorthcorner_sendbuf[i_blk] = NULL;
       } /* endif */
       if (Blk_List.message_reschange_topnorthcorner_recbuf[i_blk] != NULL) {
          delete []Blk_List.message_reschange_topnorthcorner_recbuf[i_blk];
          Blk_List.message_reschange_topnorthcorner_recbuf[i_blk] = NULL;
       } /* endif */

       // Top North  neighbour is at higher or lower resolution!
       if (Blk_List.Block[i_blk].used && 
           (Blk_List.Block[i_blk].nTN == 1) &&
           ((Blk_List.Block[i_blk].info.level < Blk_List.Block[i_blk].infoTN[0].level) ||
            (Blk_List.Block[i_blk].info.level > Blk_List.Block[i_blk].infoTN[0].level))) {
          buffer_size = sqr(Blk_List.Block[i_blk].info.dimen.ghost)*
                        Number_of_Solution_Variables;
          buffer_size_neighbour = sqr(Blk_List.Block[i_blk].infoTN[0].dimen.ghost)*
                                  Number_of_Solution_Variables;
          // Send buffer.
          Blk_List.message_reschange_topnorthcorner_sendbuf[i_blk] = 
             new double[buffer_size_neighbour];
          for (l = 0; l <= buffer_size_neighbour-1; ++l) 
	     Blk_List.message_reschange_topnorthcorner_sendbuf[i_blk][l] = ZERO;
          // Receive buffer.
          Blk_List.message_reschange_topnorthcorner_recbuf[i_blk] = 
             new double[buffer_size];
          for (l = 0; l <= buffer_size-1; ++l) 
	     Blk_List.message_reschange_topnorthcorner_recbuf[i_blk][l] = ZERO;

       // No Top North  neighbour with resolution change!
       } else {
          Blk_List.message_reschange_topnorthcorner_sendbuf[i_blk] = 
             new double[1];
          Blk_List.message_reschange_topnorthcorner_sendbuf[i_blk][0] = ZERO; 
          Blk_List.message_reschange_topnorthcorner_recbuf[i_blk] = 
             new double[1];
          Blk_List.message_reschange_topnorthcorner_recbuf[i_blk][0] = ZERO;
       } /* endif */

       ///////////////////////////////////////////////////////////////////////
       // Reallocate memory for Top South corner send and receive buffers. //
       ///////////////////////////////////////////////////////////////////////
       if (Blk_List.message_reschange_topsouthcorner_sendbuf[i_blk] != NULL) {
          delete []Blk_List.message_reschange_topsouthcorner_sendbuf[i_blk];
          Blk_List.message_reschange_topsouthcorner_sendbuf[i_blk] = NULL;
       } /* endif */
       if (Blk_List.message_reschange_topsouthcorner_recbuf[i_blk] != NULL) {
          delete []Blk_List.message_reschange_topsouthcorner_recbuf[i_blk];
          Blk_List.message_reschange_topsouthcorner_recbuf[i_blk] = NULL;
       } /* endif */

       // Top South neighbour is at higher or lower resolution!
       if (Blk_List.Block[i_blk].used && 
           (Blk_List.Block[i_blk].nTS == 1) &&
           ((Blk_List.Block[i_blk].info.level < Blk_List.Block[i_blk].infoTS[0].level) ||
            (Blk_List.Block[i_blk].info.level > Blk_List.Block[i_blk].infoTS[0].level))) {
          buffer_size = sqr(Blk_List.Block[i_blk].info.dimen.ghost)*
                        Number_of_Solution_Variables;
          buffer_size_neighbour = sqr(Blk_List.Block[i_blk].infoTS[0].dimen.ghost)*
                                  Number_of_Solution_Variables;
          // Send buffer.
          Blk_List.message_reschange_topsouthcorner_sendbuf[i_blk] = 
             new double[buffer_size_neighbour];
          for (l = 0; l <= buffer_size_neighbour-1; ++l) 
	     Blk_List.message_reschange_topsouthcorner_sendbuf[i_blk][l] = ZERO;
          // Receive buffer.
          Blk_List.message_reschange_topsouthcorner_recbuf[i_blk] = 
             new double[buffer_size];
          for (l = 0; l <= buffer_size-1; ++l) 
	     Blk_List.message_reschange_topsouthcorner_recbuf[i_blk][l] = ZERO;

       // No Top South neighbour with resolution change!
       } else {
          Blk_List.message_reschange_topsouthcorner_sendbuf[i_blk] = 
             new double[1];
          Blk_List.message_reschange_topsouthcorner_sendbuf[i_blk][0] = ZERO; 
          Blk_List.message_reschange_topsouthcorner_recbuf[i_blk] = 
             new double[1];
          Blk_List.message_reschange_topsouthcorner_recbuf[i_blk][0] = ZERO;
       } /* endif */

       ///////////////////////////////////////////////////////////////////////
       // Reallocate memory for  Top East corner send and receive buffers. //
       ///////////////////////////////////////////////////////////////////////
       if (Blk_List.message_reschange_topeastcorner_sendbuf[i_blk] != NULL) {
          delete []Blk_List.message_reschange_topeastcorner_sendbuf[i_blk];
          Blk_List.message_reschange_topeastcorner_sendbuf[i_blk] = NULL;
       } /* endif */
       if (Blk_List.message_reschange_topeastcorner_recbuf[i_blk] != NULL) {
          delete []Blk_List.message_reschange_topeastcorner_recbuf[i_blk];
          Blk_List.message_reschange_topeastcorner_recbuf[i_blk] = NULL;
       } /* endif */

       //  Top East neighbour is at higher or lower resolution!
       if (Blk_List.Block[i_blk].used && 
           (Blk_List.Block[i_blk].nTE == 1) &&
           ((Blk_List.Block[i_blk].info.level < Blk_List.Block[i_blk].infoTE[0].level) ||
            (Blk_List.Block[i_blk].info.level > Blk_List.Block[i_blk].infoTE[0].level))) {
          buffer_size = sqr(Blk_List.Block[i_blk].info.dimen.ghost)*
                        Number_of_Solution_Variables;
          buffer_size_neighbour = sqr(Blk_List.Block[i_blk].infoTE[0].dimen.ghost)*
                                  Number_of_Solution_Variables;
          // Send buffer.
          Blk_List.message_reschange_topeastcorner_sendbuf[i_blk] = 
             new double[buffer_size_neighbour];
          for (l = 0; l <= buffer_size_neighbour-1; ++l) 
	     Blk_List.message_reschange_topeastcorner_sendbuf[i_blk][l] = ZERO;
          // Receive buffer.
          Blk_List.message_reschange_topeastcorner_recbuf[i_blk] = 
             new double[buffer_size];
          for (l = 0; l <= buffer_size-1; ++l) 
	     Blk_List.message_reschange_topeastcorner_recbuf[i_blk][l] = ZERO;

       // No  Top East neighbour with resolution change!
       } else {
          Blk_List.message_reschange_topeastcorner_sendbuf[i_blk] = 
             new double[1];
          Blk_List.message_reschange_topeastcorner_sendbuf[i_blk][0] = ZERO; 
          Blk_List.message_reschange_topeastcorner_recbuf[i_blk] = 
             new double[1];
          Blk_List.message_reschange_topeastcorner_recbuf[i_blk][0] = ZERO;
       } /* endif */

       ///////////////////////////////////////////////////////////////////////
       // Reallocate memory for  Top West corner send and receive buffers. //
       ///////////////////////////////////////////////////////////////////////
       if (Blk_List.message_reschange_topwestcorner_sendbuf[i_blk] != NULL) {
          delete []Blk_List.message_reschange_topwestcorner_sendbuf[i_blk];
          Blk_List.message_reschange_topwestcorner_sendbuf[i_blk] = NULL;
       } /* endif */
       if (Blk_List.message_reschange_topwestcorner_recbuf[i_blk] != NULL) {
          delete []Blk_List.message_reschange_topwestcorner_recbuf[i_blk];
          Blk_List.message_reschange_topwestcorner_recbuf[i_blk] = NULL;
       } /* endif */

       //  Top West neighbour is at higher or lower resolution!
       if (Blk_List.Block[i_blk].used && 
           (Blk_List.Block[i_blk].nTW == 1) &&
           ((Blk_List.Block[i_blk].info.level < Blk_List.Block[i_blk].infoTW[0].level) ||
            (Blk_List.Block[i_blk].info.level > Blk_List.Block[i_blk].infoTW[0].level))) {
          buffer_size = sqr(Blk_List.Block[i_blk].info.dimen.ghost)*
                        Number_of_Solution_Variables;
          buffer_size_neighbour = sqr(Blk_List.Block[i_blk].infoTW[0].dimen.ghost)*
                                  Number_of_Solution_Variables;
          // Send buffer.
          Blk_List.message_reschange_topwestcorner_sendbuf[i_blk] = 
             new double[buffer_size_neighbour];
          for (l = 0; l <= buffer_size_neighbour-1; ++l) 
	     Blk_List.message_reschange_topwestcorner_sendbuf[i_blk][l] = ZERO;
          // Receive buffer.
          Blk_List.message_reschange_topwestcorner_recbuf[i_blk] = 
             new double[buffer_size];
          for (l = 0; l <= buffer_size-1; ++l) 
	     Blk_List.message_reschange_topwestcorner_recbuf[i_blk][l] = ZERO;

       // No  Top West neighbour with resolution change!
       } else {
          Blk_List.message_reschange_topwestcorner_sendbuf[i_blk] = 
             new double[1];
          Blk_List.message_reschange_topwestcorner_sendbuf[i_blk][0] = ZERO; 
          Blk_List.message_reschange_topwestcorner_recbuf[i_blk] = 
             new double[1];
          Blk_List.message_reschange_topwestcorner_recbuf[i_blk][0] = ZERO;
       } /* endif */

/**/

       ///////////////////////////////////////////////////////////////////////
       // Reallocate memory for Top North West corner send and receive buffers. //
       ///////////////////////////////////////////////////////////////////////
       if (Blk_List.message_reschange_topnorthwestcorner_sendbuf[i_blk] != NULL) {
          delete []Blk_List.message_reschange_topnorthwestcorner_sendbuf[i_blk];
          Blk_List.message_reschange_topnorthwestcorner_sendbuf[i_blk] = NULL;
       } /* endif */
       if (Blk_List.message_reschange_topnorthwestcorner_recbuf[i_blk] != NULL) {
          delete []Blk_List.message_reschange_topnorthwestcorner_recbuf[i_blk];
          Blk_List.message_reschange_topnorthwestcorner_recbuf[i_blk] = NULL;
       } /* endif */

       // Top North West neighbour is at higher or lower resolution!
       if (Blk_List.Block[i_blk].used && 
           (Blk_List.Block[i_blk].nTNW == 1) &&
           ((Blk_List.Block[i_blk].info.level < Blk_List.Block[i_blk].infoTNW[0].level) ||
            (Blk_List.Block[i_blk].info.level > Blk_List.Block[i_blk].infoTNW[0].level))) {
          buffer_size = sqr(Blk_List.Block[i_blk].info.dimen.ghost)*
                        Number_of_Solution_Variables;
          buffer_size_neighbour = sqr(Blk_List.Block[i_blk].infoTNW[0].dimen.ghost)*
                                  Number_of_Solution_Variables;
          // Send buffer.
          Blk_List.message_reschange_topnorthwestcorner_sendbuf[i_blk] = 
             new double[buffer_size_neighbour];
          for (l = 0; l <= buffer_size_neighbour-1; ++l) 
	     Blk_List.message_reschange_topnorthwestcorner_sendbuf[i_blk][l] = ZERO;
          // Receive buffer.
          Blk_List.message_reschange_topnorthwestcorner_recbuf[i_blk] = 
             new double[buffer_size];
          for (l = 0; l <= buffer_size-1; ++l) 
	     Blk_List.message_reschange_topnorthwestcorner_recbuf[i_blk][l] = ZERO;

       // No Top North West neighbour with resolution change!
       } else {
          Blk_List.message_reschange_topnorthwestcorner_sendbuf[i_blk] = 
             new double[1];
          Blk_List.message_reschange_topnorthwestcorner_sendbuf[i_blk][0] = ZERO; 
          Blk_List.message_reschange_topnorthwestcorner_recbuf[i_blk] = 
             new double[1];
          Blk_List.message_reschange_topnorthwestcorner_recbuf[i_blk][0] = ZERO;
       } /* endif */

       ///////////////////////////////////////////////////////////////////////
       // Reallocate memory for Top North East corner send and receive buffers. //
       ///////////////////////////////////////////////////////////////////////
       if (Blk_List.message_reschange_topnortheastcorner_sendbuf[i_blk] != NULL) {
          delete []Blk_List.message_reschange_topnortheastcorner_sendbuf[i_blk];
          Blk_List.message_reschange_topnortheastcorner_sendbuf[i_blk] = NULL;
       } /* endif */
       if (Blk_List.message_reschange_topnortheastcorner_recbuf[i_blk] != NULL) {
          delete []Blk_List.message_reschange_topnortheastcorner_recbuf[i_blk];
          Blk_List.message_reschange_topnortheastcorner_recbuf[i_blk] = NULL;
       } /* endif */

       // Top North East neighbour is at higher or lower resolution!
       if (Blk_List.Block[i_blk].used && 
           (Blk_List.Block[i_blk].nTNE == 1) &&
           ((Blk_List.Block[i_blk].info.level < Blk_List.Block[i_blk].infoTNE[0].level) ||
            (Blk_List.Block[i_blk].info.level > Blk_List.Block[i_blk].infoTNE[0].level))) {
          buffer_size = sqr(Blk_List.Block[i_blk].info.dimen.ghost)*
                        Number_of_Solution_Variables;
          buffer_size_neighbour = sqr(Blk_List.Block[i_blk].infoTNE[0].dimen.ghost)*
                                  Number_of_Solution_Variables;
          // Send buffer.
          Blk_List.message_reschange_topnortheastcorner_sendbuf[i_blk] = 
             new double[buffer_size_neighbour];
          for (l = 0; l <= buffer_size_neighbour-1; ++l) 
	     Blk_List.message_reschange_topnortheastcorner_sendbuf[i_blk][l] = ZERO;
          // Receive buffer.
          Blk_List.message_reschange_topnortheastcorner_recbuf[i_blk] = 
             new double[buffer_size];
          for (l = 0; l <= buffer_size-1; ++l) 
	     Blk_List.message_reschange_topnortheastcorner_recbuf[i_blk][l] = ZERO;

       // No Top North East neighbour with resolution change!
       } else {
          Blk_List.message_reschange_topnortheastcorner_sendbuf[i_blk] = 
             new double[1];
          Blk_List.message_reschange_topnortheastcorner_sendbuf[i_blk][0] = ZERO; 
          Blk_List.message_reschange_topnortheastcorner_recbuf[i_blk] = 
             new double[1];
          Blk_List.message_reschange_topnortheastcorner_recbuf[i_blk][0] = ZERO;
       } /* endif */

       ///////////////////////////////////////////////////////////////////////
       // Reallocate memory for Top South East corner send and receive buffers. //
       ///////////////////////////////////////////////////////////////////////
       if (Blk_List.message_reschange_topsoutheastcorner_sendbuf[i_blk] != NULL) {
          delete []Blk_List.message_reschange_topsoutheastcorner_sendbuf[i_blk];
          Blk_List.message_reschange_topsoutheastcorner_sendbuf[i_blk] = NULL;
       } /* endif */
       if (Blk_List.message_reschange_topsoutheastcorner_recbuf[i_blk] != NULL) {
          delete []Blk_List.message_reschange_topsoutheastcorner_recbuf[i_blk];
          Blk_List.message_reschange_topsoutheastcorner_recbuf[i_blk] = NULL;
       } /* endif */

       // Top South East neighbour is at higher or lower resolution!
       if (Blk_List.Block[i_blk].used && 
           (Blk_List.Block[i_blk].nTSE == 1) &&
           ((Blk_List.Block[i_blk].info.level < Blk_List.Block[i_blk].infoTSE[0].level) ||
            (Blk_List.Block[i_blk].info.level > Blk_List.Block[i_blk].infoTSE[0].level))) {
          buffer_size = sqr(Blk_List.Block[i_blk].info.dimen.ghost)*
                        Number_of_Solution_Variables;
          buffer_size_neighbour = sqr(Blk_List.Block[i_blk].infoTSE[0].dimen.ghost)*
                                  Number_of_Solution_Variables;
          // Send buffer.
          Blk_List.message_reschange_topsoutheastcorner_sendbuf[i_blk] = 
             new double[buffer_size_neighbour];
          for (l = 0; l <= buffer_size_neighbour-1; ++l) 
	     Blk_List.message_reschange_topsoutheastcorner_sendbuf[i_blk][l] = ZERO;
          // Receive buffer.
          Blk_List.message_reschange_topsoutheastcorner_recbuf[i_blk] = 
             new double[buffer_size];
          for (l = 0; l <= buffer_size-1; ++l) 
	     Blk_List.message_reschange_topsoutheastcorner_recbuf[i_blk][l] = ZERO;

       // No Top South East neighbour with resolution change!
       } else {
          Blk_List.message_reschange_topsoutheastcorner_sendbuf[i_blk] = 
             new double[1];
          Blk_List.message_reschange_topsoutheastcorner_sendbuf[i_blk][0] = ZERO; 
          Blk_List.message_reschange_topsoutheastcorner_recbuf[i_blk] = 
             new double[1];
          Blk_List.message_reschange_topsoutheastcorner_recbuf[i_blk][0] = ZERO;
       } /* endif */

       ///////////////////////////////////////////////////////////////////////
       // Reallocate memory for Top South West corner send and receive buffers. //
       ///////////////////////////////////////////////////////////////////////
       if (Blk_List.message_reschange_topsouthwestcorner_sendbuf[i_blk] != NULL) {
          delete []Blk_List.message_reschange_topsouthwestcorner_sendbuf[i_blk];
          Blk_List.message_reschange_topsouthwestcorner_sendbuf[i_blk] = NULL;
       } /* endif */
       if (Blk_List.message_reschange_topsouthwestcorner_recbuf[i_blk] != NULL) {
          delete []Blk_List.message_reschange_topsouthwestcorner_recbuf[i_blk];
          Blk_List.message_reschange_topsouthwestcorner_recbuf[i_blk] = NULL;
       } /* endif */

       // Top South West neighbour is at higher or lower resolution!
       if (Blk_List.Block[i_blk].used && 
           (Blk_List.Block[i_blk].nTSW == 1) &&
           ((Blk_List.Block[i_blk].info.level < Blk_List.Block[i_blk].infoTSW[0].level) ||
            (Blk_List.Block[i_blk].info.level > Blk_List.Block[i_blk].infoTSW[0].level))) {
          buffer_size = sqr(Blk_List.Block[i_blk].info.dimen.ghost)*
                        Number_of_Solution_Variables;
          buffer_size_neighbour = sqr(Blk_List.Block[i_blk].infoTSW[0].dimen.ghost)*
                                  Number_of_Solution_Variables;
          // Send buffer.
          Blk_List.message_reschange_topsouthwestcorner_sendbuf[i_blk] = 
             new double[buffer_size_neighbour];
          for (l = 0; l <= buffer_size_neighbour-1; ++l) 
	     Blk_List.message_reschange_topsouthwestcorner_sendbuf[i_blk][l] = ZERO;
          // Receive buffer.
          Blk_List.message_reschange_topsouthwestcorner_recbuf[i_blk] = 
             new double[buffer_size];
          for (l = 0; l <= buffer_size-1; ++l) 
	     Blk_List.message_reschange_topsouthwestcorner_recbuf[i_blk][l] = ZERO;

       // No Top South West neighbour with resolution change!
       } else {
          Blk_List.message_reschange_topsouthwestcorner_sendbuf[i_blk] = 
             new double[1];
          Blk_List.message_reschange_topsouthwestcorner_sendbuf[i_blk][0] = ZERO; 
          Blk_List.message_reschange_topsouthwestcorner_recbuf[i_blk] = 
             new double[1];
          Blk_List.message_reschange_topsouthwestcorner_recbuf[i_blk][0] = ZERO;
       } /* endif */

    ///******************//

       ///////////////////////////////////////////////////////////////////////
       // Reallocate memory for Bottom North corner send and receive buffers. //
       ///////////////////////////////////////////////////////////////////////
       if (Blk_List.message_reschange_bottomnorthcorner_sendbuf[i_blk] != NULL) {
          delete []Blk_List.message_reschange_bottomnorthcorner_sendbuf[i_blk];
          Blk_List.message_reschange_bottomnorthcorner_sendbuf[i_blk] = NULL;
       } /* endif */
       if (Blk_List.message_reschange_bottomnorthcorner_recbuf[i_blk] != NULL) {
          delete []Blk_List.message_reschange_bottomnorthcorner_recbuf[i_blk];
          Blk_List.message_reschange_bottomnorthcorner_recbuf[i_blk] = NULL;
       } /* endif */

       // Bottom North  neighbour is at higher or lower resolution!
       if (Blk_List.Block[i_blk].used && 
           (Blk_List.Block[i_blk].nBN == 1) &&
           ((Blk_List.Block[i_blk].info.level < Blk_List.Block[i_blk].infoBN[0].level) ||
            (Blk_List.Block[i_blk].info.level > Blk_List.Block[i_blk].infoBN[0].level))) {
          buffer_size = sqr(Blk_List.Block[i_blk].info.dimen.ghost)*
                        Number_of_Solution_Variables;
          buffer_size_neighbour = sqr(Blk_List.Block[i_blk].infoBN[0].dimen.ghost)*
                                  Number_of_Solution_Variables;
          // Send buffer.
          Blk_List.message_reschange_bottomnorthcorner_sendbuf[i_blk] = 
             new double[buffer_size_neighbour];
          for (l = 0; l <= buffer_size_neighbour-1; ++l) 
	     Blk_List.message_reschange_bottomnorthcorner_sendbuf[i_blk][l] = ZERO;
          // Receive buffer.
          Blk_List.message_reschange_bottomnorthcorner_recbuf[i_blk] = 
             new double[buffer_size];
          for (l = 0; l <= buffer_size-1; ++l) 
	     Blk_List.message_reschange_bottomnorthcorner_recbuf[i_blk][l] = ZERO;

       // No Bottom North  neighbour with resolution change!
       } else {
          Blk_List.message_reschange_bottomnorthcorner_sendbuf[i_blk] = 
             new double[1];
          Blk_List.message_reschange_bottomnorthcorner_sendbuf[i_blk][0] = ZERO; 
          Blk_List.message_reschange_bottomnorthcorner_recbuf[i_blk] = 
             new double[1];
          Blk_List.message_reschange_bottomnorthcorner_recbuf[i_blk][0] = ZERO;
       } /* endif */

       ///////////////////////////////////////////////////////////////////////
       // Reallocate memory for Bottom South corner send and receive buffers. //
       ///////////////////////////////////////////////////////////////////////
       if (Blk_List.message_reschange_bottomsouthcorner_sendbuf[i_blk] != NULL) {
          delete []Blk_List.message_reschange_bottomsouthcorner_sendbuf[i_blk];
          Blk_List.message_reschange_bottomsouthcorner_sendbuf[i_blk] = NULL;
       } /* endif */
       if (Blk_List.message_reschange_bottomsouthcorner_recbuf[i_blk] != NULL) {
          delete []Blk_List.message_reschange_bottomsouthcorner_recbuf[i_blk];
          Blk_List.message_reschange_bottomsouthcorner_recbuf[i_blk] = NULL;
       } /* endif */

       // Bottom South neighbour is at higher or lower resolution!
       if (Blk_List.Block[i_blk].used && 
           (Blk_List.Block[i_blk].nBS == 1) &&
           ((Blk_List.Block[i_blk].info.level < Blk_List.Block[i_blk].infoBS[0].level) ||
            (Blk_List.Block[i_blk].info.level > Blk_List.Block[i_blk].infoBS[0].level))) {
          buffer_size = sqr(Blk_List.Block[i_blk].info.dimen.ghost)*
                        Number_of_Solution_Variables;
          buffer_size_neighbour = sqr(Blk_List.Block[i_blk].infoBS[0].dimen.ghost)*
                                  Number_of_Solution_Variables;
          // Send buffer.
          Blk_List.message_reschange_bottomsouthcorner_sendbuf[i_blk] = 
             new double[buffer_size_neighbour];
          for (l = 0; l <= buffer_size_neighbour-1; ++l) 
	     Blk_List.message_reschange_bottomsouthcorner_sendbuf[i_blk][l] = ZERO;
          // Receive buffer.
          Blk_List.message_reschange_bottomsouthcorner_recbuf[i_blk] = 
             new double[buffer_size];
          for (l = 0; l <= buffer_size-1; ++l) 
	     Blk_List.message_reschange_bottomsouthcorner_recbuf[i_blk][l] = ZERO;

       // No Bottom South neighbour with resolution change!
       } else {
          Blk_List.message_reschange_bottomsouthcorner_sendbuf[i_blk] = 
             new double[1];
          Blk_List.message_reschange_bottomsouthcorner_sendbuf[i_blk][0] = ZERO; 
          Blk_List.message_reschange_bottomsouthcorner_recbuf[i_blk] = 
             new double[1];
          Blk_List.message_reschange_bottomsouthcorner_recbuf[i_blk][0] = ZERO;
       } /* endif */

       ///////////////////////////////////////////////////////////////////////
       // Reallocate memory for  Bottom East corner send and receive buffers. //
       ///////////////////////////////////////////////////////////////////////
       if (Blk_List.message_reschange_bottomeastcorner_sendbuf[i_blk] != NULL) {
          delete []Blk_List.message_reschange_bottomeastcorner_sendbuf[i_blk];
          Blk_List.message_reschange_bottomeastcorner_sendbuf[i_blk] = NULL;
       } /* endif */
       if (Blk_List.message_reschange_bottomeastcorner_recbuf[i_blk] != NULL) {
          delete []Blk_List.message_reschange_bottomeastcorner_recbuf[i_blk];
          Blk_List.message_reschange_bottomeastcorner_recbuf[i_blk] = NULL;
       } /* endif */

       //  Bottom East neighbour is at higher or lower resolution!
       if (Blk_List.Block[i_blk].used && 
           (Blk_List.Block[i_blk].nBE == 1) &&
           ((Blk_List.Block[i_blk].info.level < Blk_List.Block[i_blk].infoBE[0].level) ||
            (Blk_List.Block[i_blk].info.level > Blk_List.Block[i_blk].infoBE[0].level))) {
          buffer_size = sqr(Blk_List.Block[i_blk].info.dimen.ghost)*
                        Number_of_Solution_Variables;
          buffer_size_neighbour = sqr(Blk_List.Block[i_blk].infoBE[0].dimen.ghost)*
                                  Number_of_Solution_Variables;
          // Send buffer.
          Blk_List.message_reschange_bottomeastcorner_sendbuf[i_blk] = 
             new double[buffer_size_neighbour];
          for (l = 0; l <= buffer_size_neighbour-1; ++l) 
	     Blk_List.message_reschange_bottomeastcorner_sendbuf[i_blk][l] = ZERO;
          // Receive buffer.
          Blk_List.message_reschange_bottomeastcorner_recbuf[i_blk] = 
             new double[buffer_size];
          for (l = 0; l <= buffer_size-1; ++l) 
	     Blk_List.message_reschange_bottomeastcorner_recbuf[i_blk][l] = ZERO;

       // No  Bottom East neighbour with resolution change!
       } else {
          Blk_List.message_reschange_bottomeastcorner_sendbuf[i_blk] = 
             new double[1];
          Blk_List.message_reschange_bottomeastcorner_sendbuf[i_blk][0] = ZERO; 
          Blk_List.message_reschange_bottomeastcorner_recbuf[i_blk] = 
             new double[1];
          Blk_List.message_reschange_bottomeastcorner_recbuf[i_blk][0] = ZERO;
       } /* endif */

       ///////////////////////////////////////////////////////////////////////
       // Reallocate memory for  Bottom West corner send and receive buffers. //
       ///////////////////////////////////////////////////////////////////////
       if (Blk_List.message_reschange_bottomwestcorner_sendbuf[i_blk] != NULL) {
          delete []Blk_List.message_reschange_bottomwestcorner_sendbuf[i_blk];
          Blk_List.message_reschange_bottomwestcorner_sendbuf[i_blk] = NULL;
       } /* endif */
       if (Blk_List.message_reschange_bottomwestcorner_recbuf[i_blk] != NULL) {
          delete []Blk_List.message_reschange_bottomwestcorner_recbuf[i_blk];
          Blk_List.message_reschange_bottomwestcorner_recbuf[i_blk] = NULL;
       } /* endif */

       //  Bottom West neighbour is at higher or lower resolution!
       if (Blk_List.Block[i_blk].used && 
           (Blk_List.Block[i_blk].nBW == 1) &&
           ((Blk_List.Block[i_blk].info.level < Blk_List.Block[i_blk].infoBW[0].level) ||
            (Blk_List.Block[i_blk].info.level > Blk_List.Block[i_blk].infoBW[0].level))) {
          buffer_size = sqr(Blk_List.Block[i_blk].info.dimen.ghost)*
                        Number_of_Solution_Variables;
          buffer_size_neighbour = sqr(Blk_List.Block[i_blk].infoBW[0].dimen.ghost)*
                                  Number_of_Solution_Variables;
          // Send buffer.
          Blk_List.message_reschange_bottomwestcorner_sendbuf[i_blk] = 
             new double[buffer_size_neighbour];
          for (l = 0; l <= buffer_size_neighbour-1; ++l) 
	     Blk_List.message_reschange_bottomwestcorner_sendbuf[i_blk][l] = ZERO;
          // Receive buffer.
          Blk_List.message_reschange_bottomwestcorner_recbuf[i_blk] = 
             new double[buffer_size];
          for (l = 0; l <= buffer_size-1; ++l) 
	     Blk_List.message_reschange_bottomwestcorner_recbuf[i_blk][l] = ZERO;

       // No  Bottom West neighbour with resolution change!
       } else {
          Blk_List.message_reschange_bottomwestcorner_sendbuf[i_blk] = 
             new double[1];
          Blk_List.message_reschange_bottomwestcorner_sendbuf[i_blk][0] = ZERO; 
          Blk_List.message_reschange_bottomwestcorner_recbuf[i_blk] = 
             new double[1];
          Blk_List.message_reschange_bottomwestcorner_recbuf[i_blk][0] = ZERO;
       } /* endif */

/**/

       ///////////////////////////////////////////////////////////////////////
       // Reallocate memory for Bottom North West corner send and receive buffers. //
       ///////////////////////////////////////////////////////////////////////
       if (Blk_List.message_reschange_bottomnorthwestcorner_sendbuf[i_blk] != NULL) {
          delete []Blk_List.message_reschange_bottomnorthwestcorner_sendbuf[i_blk];
          Blk_List.message_reschange_bottomnorthwestcorner_sendbuf[i_blk] = NULL;
       } /* endif */
       if (Blk_List.message_reschange_bottomnorthwestcorner_recbuf[i_blk] != NULL) {
          delete []Blk_List.message_reschange_bottomnorthwestcorner_recbuf[i_blk];
          Blk_List.message_reschange_bottomnorthwestcorner_recbuf[i_blk] = NULL;
       } /* endif */

       // Bottom North West neighbour is at higher or lower resolution!
       if (Blk_List.Block[i_blk].used && 
           (Blk_List.Block[i_blk].nBNW == 1) &&
           ((Blk_List.Block[i_blk].info.level < Blk_List.Block[i_blk].infoBNW[0].level) ||
            (Blk_List.Block[i_blk].info.level > Blk_List.Block[i_blk].infoBNW[0].level))) {
          buffer_size = sqr(Blk_List.Block[i_blk].info.dimen.ghost)*
                        Number_of_Solution_Variables;
          buffer_size_neighbour = sqr(Blk_List.Block[i_blk].infoBNW[0].dimen.ghost)*
                                  Number_of_Solution_Variables;
          // Send buffer.
          Blk_List.message_reschange_bottomnorthwestcorner_sendbuf[i_blk] = 
             new double[buffer_size_neighbour];
          for (l = 0; l <= buffer_size_neighbour-1; ++l) 
	     Blk_List.message_reschange_bottomnorthwestcorner_sendbuf[i_blk][l] = ZERO;
          // Receive buffer.
          Blk_List.message_reschange_bottomnorthwestcorner_recbuf[i_blk] = 
             new double[buffer_size];
          for (l = 0; l <= buffer_size-1; ++l) 
	     Blk_List.message_reschange_bottomnorthwestcorner_recbuf[i_blk][l] = ZERO;

       // No Bottom North West neighbour with resolution change!
       } else {
          Blk_List.message_reschange_bottomnorthwestcorner_sendbuf[i_blk] = 
             new double[1];
          Blk_List.message_reschange_bottomnorthwestcorner_sendbuf[i_blk][0] = ZERO; 
          Blk_List.message_reschange_bottomnorthwestcorner_recbuf[i_blk] = 
             new double[1];
          Blk_List.message_reschange_bottomnorthwestcorner_recbuf[i_blk][0] = ZERO;
       } /* endif */

       ///////////////////////////////////////////////////////////////////////
       // Reallocate memory for Bottom North East corner send and receive buffers. //
       ///////////////////////////////////////////////////////////////////////
       if (Blk_List.message_reschange_bottomnortheastcorner_sendbuf[i_blk] != NULL) {
          delete []Blk_List.message_reschange_bottomnortheastcorner_sendbuf[i_blk];
          Blk_List.message_reschange_bottomnortheastcorner_sendbuf[i_blk] = NULL;
       } /* endif */
       if (Blk_List.message_reschange_bottomnortheastcorner_recbuf[i_blk] != NULL) {
          delete []Blk_List.message_reschange_bottomnortheastcorner_recbuf[i_blk];
          Blk_List.message_reschange_bottomnortheastcorner_recbuf[i_blk] = NULL;
       } /* endif */

       // Bottom North East neighbour is at higher or lower resolution!
       if (Blk_List.Block[i_blk].used && 
           (Blk_List.Block[i_blk].nBNE == 1) &&
           ((Blk_List.Block[i_blk].info.level < Blk_List.Block[i_blk].infoBNE[0].level) ||
            (Blk_List.Block[i_blk].info.level > Blk_List.Block[i_blk].infoBNE[0].level))) {
          buffer_size = sqr(Blk_List.Block[i_blk].info.dimen.ghost)*
                        Number_of_Solution_Variables;
          buffer_size_neighbour = sqr(Blk_List.Block[i_blk].infoBNE[0].dimen.ghost)*
                                  Number_of_Solution_Variables;
          // Send buffer.
          Blk_List.message_reschange_bottomnortheastcorner_sendbuf[i_blk] = 
             new double[buffer_size_neighbour];
          for (l = 0; l <= buffer_size_neighbour-1; ++l) 
	     Blk_List.message_reschange_bottomnortheastcorner_sendbuf[i_blk][l] = ZERO;
          // Receive buffer.
          Blk_List.message_reschange_bottomnortheastcorner_recbuf[i_blk] = 
             new double[buffer_size];
          for (l = 0; l <= buffer_size-1; ++l) 
	     Blk_List.message_reschange_bottomnortheastcorner_recbuf[i_blk][l] = ZERO;

       // No Bottom North East neighbour with resolution change!
       } else {
          Blk_List.message_reschange_bottomnortheastcorner_sendbuf[i_blk] = 
             new double[1];
          Blk_List.message_reschange_bottomnortheastcorner_sendbuf[i_blk][0] = ZERO; 
          Blk_List.message_reschange_bottomnortheastcorner_recbuf[i_blk] = 
             new double[1];
          Blk_List.message_reschange_bottomnortheastcorner_recbuf[i_blk][0] = ZERO;
       } /* endif */

       ///////////////////////////////////////////////////////////////////////
       // Reallocate memory for Bottom South East corner send and receive buffers. //
       ///////////////////////////////////////////////////////////////////////
       if (Blk_List.message_reschange_bottomsoutheastcorner_sendbuf[i_blk] != NULL) {
          delete []Blk_List.message_reschange_bottomsoutheastcorner_sendbuf[i_blk];
          Blk_List.message_reschange_bottomsoutheastcorner_sendbuf[i_blk] = NULL;
       } /* endif */
       if (Blk_List.message_reschange_bottomsoutheastcorner_recbuf[i_blk] != NULL) {
          delete []Blk_List.message_reschange_bottomsoutheastcorner_recbuf[i_blk];
          Blk_List.message_reschange_bottomsoutheastcorner_recbuf[i_blk] = NULL;
       } /* endif */

       // Bottom South East neighbour is at higher or lower resolution!
       if (Blk_List.Block[i_blk].used && 
           (Blk_List.Block[i_blk].nBSE == 1) &&
           ((Blk_List.Block[i_blk].info.level < Blk_List.Block[i_blk].infoBSE[0].level) ||
            (Blk_List.Block[i_blk].info.level > Blk_List.Block[i_blk].infoBSE[0].level))) {
          buffer_size = sqr(Blk_List.Block[i_blk].info.dimen.ghost)*
                        Number_of_Solution_Variables;
          buffer_size_neighbour = sqr(Blk_List.Block[i_blk].infoBSE[0].dimen.ghost)*
                                  Number_of_Solution_Variables;
          // Send buffer.
          Blk_List.message_reschange_bottomsoutheastcorner_sendbuf[i_blk] = 
             new double[buffer_size_neighbour];
          for (l = 0; l <= buffer_size_neighbour-1; ++l) 
	     Blk_List.message_reschange_bottomsoutheastcorner_sendbuf[i_blk][l] = ZERO;
          // Receive buffer.
          Blk_List.message_reschange_bottomsoutheastcorner_recbuf[i_blk] = 
             new double[buffer_size];
          for (l = 0; l <= buffer_size-1; ++l) 
	     Blk_List.message_reschange_bottomsoutheastcorner_recbuf[i_blk][l] = ZERO;

       // No Bottom South East neighbour with resolution change!
       } else {
          Blk_List.message_reschange_bottomsoutheastcorner_sendbuf[i_blk] = 
             new double[1];
          Blk_List.message_reschange_bottomsoutheastcorner_sendbuf[i_blk][0] = ZERO; 
          Blk_List.message_reschange_bottomsoutheastcorner_recbuf[i_blk] = 
             new double[1];
          Blk_List.message_reschange_bottomsoutheastcorner_recbuf[i_blk][0] = ZERO;
       } /* endif */

       ///////////////////////////////////////////////////////////////////////
       // Reallocate memory for Bottom South West corner send and receive buffers. //
       ///////////////////////////////////////////////////////////////////////
       if (Blk_List.message_reschange_bottomsouthwestcorner_sendbuf[i_blk] != NULL) {
          delete []Blk_List.message_reschange_bottomsouthwestcorner_sendbuf[i_blk];
          Blk_List.message_reschange_bottomsouthwestcorner_sendbuf[i_blk] = NULL;
       } /* endif */
       if (Blk_List.message_reschange_bottomsouthwestcorner_recbuf[i_blk] != NULL) {
          delete []Blk_List.message_reschange_bottomsouthwestcorner_recbuf[i_blk];
          Blk_List.message_reschange_bottomsouthwestcorner_recbuf[i_blk] = NULL;
       } /* endif */

       // Bottom South West neighbour is at higher or lower resolution!
       if (Blk_List.Block[i_blk].used && 
           (Blk_List.Block[i_blk].nBSW == 1) &&
           ((Blk_List.Block[i_blk].info.level < Blk_List.Block[i_blk].infoBSW[0].level) ||
            (Blk_List.Block[i_blk].info.level > Blk_List.Block[i_blk].infoBSW[0].level))) {
          buffer_size = sqr(Blk_List.Block[i_blk].info.dimen.ghost)*
                        Number_of_Solution_Variables;
          buffer_size_neighbour = sqr(Blk_List.Block[i_blk].infoBSW[0].dimen.ghost)*
                                  Number_of_Solution_Variables;
          // Send buffer.
          Blk_List.message_reschange_bottomsouthwestcorner_sendbuf[i_blk] = 
             new double[buffer_size_neighbour];
          for (l = 0; l <= buffer_size_neighbour-1; ++l) 
	     Blk_List.message_reschange_bottomsouthwestcorner_sendbuf[i_blk][l] = ZERO;
          // Receive buffer.
          Blk_List.message_reschange_bottomsouthwestcorner_recbuf[i_blk] = 
             new double[buffer_size];
          for (l = 0; l <= buffer_size-1; ++l) 
	     Blk_List.message_reschange_bottomsouthwestcorner_recbuf[i_blk][l] = ZERO;

       // No Bottom South West neighbour with resolution change!
       } else {
          Blk_List.message_reschange_bottomsouthwestcorner_sendbuf[i_blk] = 
             new double[1];
          Blk_List.message_reschange_bottomsouthwestcorner_sendbuf[i_blk][0] = ZERO; 
          Blk_List.message_reschange_bottomsouthwestcorner_recbuf[i_blk] = 
             new double[1];
          Blk_List.message_reschange_bottomsouthwestcorner_recbuf[i_blk][0] = ZERO;
       } /* endif */

    } /* endfor */
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

    int i_blk;

    // Deallocate memory for Top face send and receive buffers.
    if (Blk_List.message_noreschange_topface_sendbuf != NULL) {
       for ( i_blk = 0; i_blk <= Blk_List.Nblk-1; ++i_blk ) {
          delete []Blk_List.message_noreschange_topface_sendbuf[i_blk];
          Blk_List.message_noreschange_topface_sendbuf[i_blk] = NULL;
       } /* endfor */
       delete []Blk_List.message_noreschange_topface_sendbuf;
       Blk_List.message_noreschange_topface_sendbuf = NULL;
    } /* endif */
    if (Blk_List.message_noreschange_topface_recbuf != NULL) {
       for ( i_blk = 0; i_blk <= Blk_List.Nblk-1; ++i_blk ) {
          delete []Blk_List.message_noreschange_topface_recbuf[i_blk];
          Blk_List.message_noreschange_topface_recbuf[i_blk] = NULL;
       } /* endfor */
       delete []Blk_List.message_noreschange_topface_recbuf;
       Blk_List.message_noreschange_topface_recbuf = NULL;
    } /* endif */

     // Deallocate memory for Bottom face send and receive buffers.
    if (Blk_List.message_noreschange_bottomface_sendbuf != NULL) {
       for ( i_blk = 0; i_blk <= Blk_List.Nblk-1; ++i_blk ) {
          delete []Blk_List.message_noreschange_bottomface_sendbuf[i_blk];
          Blk_List.message_noreschange_bottomface_sendbuf[i_blk] = NULL;
       } /* endfor */
       delete []Blk_List.message_noreschange_bottomface_sendbuf;
       Blk_List.message_noreschange_bottomface_sendbuf = NULL;
    } /* endif */
    if (Blk_List.message_noreschange_bottomface_recbuf != NULL) {
       for ( i_blk = 0; i_blk <= Blk_List.Nblk-1; ++i_blk ) {
          delete []Blk_List.message_noreschange_bottomface_recbuf[i_blk];
          Blk_List.message_noreschange_bottomface_recbuf[i_blk] = NULL;
       } /* endfor */
       delete []Blk_List.message_noreschange_bottomface_recbuf;
       Blk_List.message_noreschange_bottomface_recbuf = NULL;
    } /* endif */

     // Deallocate memory for North face send and receive buffers.
    if (Blk_List.message_noreschange_northface_sendbuf != NULL) {
       for ( i_blk = 0; i_blk <= Blk_List.Nblk-1; ++i_blk ) {
          delete []Blk_List.message_noreschange_northface_sendbuf[i_blk];
          Blk_List.message_noreschange_northface_sendbuf[i_blk] = NULL;
       } /* endfor */
       delete []Blk_List.message_noreschange_northface_sendbuf;
       Blk_List.message_noreschange_northface_sendbuf = NULL;
    } /* endif */
    if (Blk_List.message_noreschange_northface_recbuf != NULL) {
       for ( i_blk = 0; i_blk <= Blk_List.Nblk-1; ++i_blk ) {
          delete []Blk_List.message_noreschange_northface_recbuf[i_blk];
          Blk_List.message_noreschange_northface_recbuf[i_blk] = NULL;
       } /* endfor */
       delete []Blk_List.message_noreschange_northface_recbuf;
       Blk_List.message_noreschange_northface_recbuf = NULL;
    } /* endif */

    // Deallocate memory for South face send and receive buffers.
    if (Blk_List.message_noreschange_southface_sendbuf != NULL) {
       for ( i_blk = 0; i_blk <= Blk_List.Nblk-1; ++i_blk ) {
          delete []Blk_List.message_noreschange_southface_sendbuf[i_blk];
          Blk_List.message_noreschange_southface_sendbuf[i_blk] = NULL;
       } /* endfor */
       delete []Blk_List.message_noreschange_southface_sendbuf;
       Blk_List.message_noreschange_southface_sendbuf = NULL;
    } /* endif */
    if (Blk_List.message_noreschange_southface_recbuf != NULL) {
       for ( i_blk = 0; i_blk <= Blk_List.Nblk-1; ++i_blk ) {
          delete []Blk_List.message_noreschange_southface_recbuf[i_blk];
          Blk_List.message_noreschange_southface_recbuf[i_blk] = NULL;
       } /* endfor */
       delete []Blk_List.message_noreschange_southface_recbuf;
       Blk_List.message_noreschange_southface_recbuf = NULL;
    } /* endif */

    // Deallocate memory for East face send and receive buffers.
    if (Blk_List.message_noreschange_eastface_sendbuf != NULL) {
       for ( i_blk = 0; i_blk <= Blk_List.Nblk-1; ++i_blk ) {
          delete []Blk_List.message_noreschange_eastface_sendbuf[i_blk];
          Blk_List.message_noreschange_eastface_sendbuf[i_blk] = NULL;
       } /* endfor */
       delete []Blk_List.message_noreschange_eastface_sendbuf;
       Blk_List.message_noreschange_eastface_sendbuf = NULL;
    } /* endif */
    if (Blk_List.message_noreschange_eastface_recbuf != NULL) {
       for ( i_blk = 0; i_blk <= Blk_List.Nblk-1; ++i_blk ) {
          delete []Blk_List.message_noreschange_eastface_recbuf[i_blk];
          Blk_List.message_noreschange_eastface_recbuf[i_blk] = NULL;
       } /* endfor */
       delete []Blk_List.message_noreschange_eastface_recbuf;
       Blk_List.message_noreschange_eastface_recbuf = NULL;
    } /* endif */

    // Deallocate memory for West face send and receive buffers.
    if (Blk_List.message_noreschange_westface_sendbuf != NULL) {
       for ( i_blk = 0; i_blk <= Blk_List.Nblk-1; ++i_blk ) {
          delete []Blk_List.message_noreschange_westface_sendbuf[i_blk];
          Blk_List.message_noreschange_westface_sendbuf[i_blk] = NULL;
       } /* endfor */
       delete []Blk_List.message_noreschange_westface_sendbuf;
       Blk_List.message_noreschange_westface_sendbuf = NULL;
    } /* endif */
    if (Blk_List.message_noreschange_westface_recbuf != NULL) {
       for ( i_blk = 0; i_blk <= Blk_List.Nblk-1; ++i_blk ) {
          delete []Blk_List.message_noreschange_westface_recbuf[i_blk];
          Blk_List.message_noreschange_westface_recbuf[i_blk] = NULL;
       } /* endfor */
       delete []Blk_List.message_noreschange_westface_recbuf;
       Blk_List.message_noreschange_westface_recbuf = NULL;
    } /* endif */

    // Deallocate memory for North West corner send and receive buffers.
    if (Blk_List.message_noreschange_northwestcorner_sendbuf != NULL) {
       for ( i_blk = 0; i_blk <= Blk_List.Nblk-1; ++i_blk ) {
          delete []Blk_List.message_noreschange_northwestcorner_sendbuf[i_blk];
          Blk_List.message_noreschange_northwestcorner_sendbuf[i_blk] = NULL;
       } /* endfor */
       delete []Blk_List.message_noreschange_northwestcorner_sendbuf;
       Blk_List.message_noreschange_northwestcorner_sendbuf = NULL;
    } /* endif */
    if (Blk_List.message_noreschange_northwestcorner_recbuf != NULL) {
       for ( i_blk = 0; i_blk <= Blk_List.Nblk-1; ++i_blk ) {
          delete []Blk_List.message_noreschange_northwestcorner_recbuf[i_blk];
          Blk_List.message_noreschange_northwestcorner_recbuf[i_blk] = NULL;
       } /* endfor */
       delete []Blk_List.message_noreschange_northwestcorner_recbuf;
       Blk_List.message_noreschange_northwestcorner_recbuf = NULL;
    } /* endif */

    // Deallocate memory for North East corner send and receive buffers.
    if (Blk_List.message_noreschange_northeastcorner_sendbuf != NULL) {
       for ( i_blk = 0; i_blk <= Blk_List.Nblk-1; ++i_blk ) {
          delete []Blk_List.message_noreschange_northeastcorner_sendbuf[i_blk];
          Blk_List.message_noreschange_northeastcorner_sendbuf[i_blk] = NULL;
       } /* endfor */
       delete []Blk_List.message_noreschange_northeastcorner_sendbuf;
       Blk_List.message_noreschange_northeastcorner_sendbuf = NULL;
    } /* endif */
    if (Blk_List.message_noreschange_northeastcorner_recbuf != NULL) {
       for ( i_blk = 0; i_blk <= Blk_List.Nblk-1; ++i_blk ) {
          delete []Blk_List.message_noreschange_northeastcorner_recbuf[i_blk];
          Blk_List.message_noreschange_northeastcorner_recbuf[i_blk] = NULL;
       } /* endfor */
       delete []Blk_List.message_noreschange_northeastcorner_recbuf;
       Blk_List.message_noreschange_northeastcorner_recbuf = NULL;
    } /* endif */

    // Deallocate memory for South East corner send and receive buffers.
    if (Blk_List.message_noreschange_southeastcorner_sendbuf != NULL) {
       for ( i_blk = 0; i_blk <= Blk_List.Nblk-1; ++i_blk ) {
          delete []Blk_List.message_noreschange_southeastcorner_sendbuf[i_blk];
          Blk_List.message_noreschange_southeastcorner_sendbuf[i_blk] = NULL;
       } /* endfor */
       delete []Blk_List.message_noreschange_southeastcorner_sendbuf;
       Blk_List.message_noreschange_southeastcorner_sendbuf = NULL;
    } /* endif */
    if (Blk_List.message_noreschange_southeastcorner_recbuf != NULL) {
       for ( i_blk = 0; i_blk <= Blk_List.Nblk-1; ++i_blk ) {
          delete []Blk_List.message_noreschange_southeastcorner_recbuf[i_blk];
          Blk_List.message_noreschange_southeastcorner_recbuf[i_blk] = NULL;
       } /* endfor */
       delete []Blk_List.message_noreschange_southeastcorner_recbuf;
       Blk_List.message_noreschange_southeastcorner_recbuf = NULL;
    } /* endif */

    // Deallocate memory for South West corner send and receive buffers.
    if (Blk_List.message_noreschange_southwestcorner_sendbuf != NULL) {
       for ( i_blk = 0; i_blk <= Blk_List.Nblk-1; ++i_blk ) {
          delete []Blk_List.message_noreschange_southwestcorner_sendbuf[i_blk];
          Blk_List.message_noreschange_southwestcorner_sendbuf[i_blk] = NULL;
       } /* endfor */
       delete []Blk_List.message_noreschange_southwestcorner_sendbuf;
       Blk_List.message_noreschange_southwestcorner_sendbuf = NULL;
    } /* endif */
    if (Blk_List.message_noreschange_southwestcorner_recbuf != NULL) {
       for ( i_blk = 0; i_blk <= Blk_List.Nblk-1; ++i_blk ) {
          delete []Blk_List.message_noreschange_southwestcorner_recbuf[i_blk];
          Blk_List.message_noreschange_southwestcorner_recbuf[i_blk] = NULL;
       } /* endfor */
       delete []Blk_List.message_noreschange_southwestcorner_recbuf;
       Blk_List.message_noreschange_southwestcorner_recbuf = NULL;
    } /* endif */

    //*************//

     // Deallocate memory for Top North corner send and receive buffers.
    if (Blk_List.message_noreschange_topnorthcorner_sendbuf != NULL) {
       for ( i_blk = 0; i_blk <= Blk_List.Nblk-1; ++i_blk ) {
          delete []Blk_List.message_noreschange_topnorthcorner_sendbuf[i_blk];
          Blk_List.message_noreschange_topnorthcorner_sendbuf[i_blk] = NULL;
       } /* endfor */
       delete []Blk_List.message_noreschange_topnorthcorner_sendbuf;
       Blk_List.message_noreschange_topnorthcorner_sendbuf = NULL;
    } /* endif */
    if (Blk_List.message_noreschange_topnorthcorner_recbuf != NULL) {
       for ( i_blk = 0; i_blk <= Blk_List.Nblk-1; ++i_blk ) {
          delete []Blk_List.message_noreschange_topnorthcorner_recbuf[i_blk];
          Blk_List.message_noreschange_topnorthcorner_recbuf[i_blk] = NULL;
       } /* endfor */
       delete []Blk_List.message_noreschange_topnorthcorner_recbuf;
       Blk_List.message_noreschange_topnorthcorner_recbuf = NULL;
    } /* endif */

    // Deallocate memory for Top South  corner send and receive buffers.
    if (Blk_List.message_noreschange_topsouthcorner_sendbuf != NULL) {
       for ( i_blk = 0; i_blk <= Blk_List.Nblk-1; ++i_blk ) {
          delete []Blk_List.message_noreschange_topsouthcorner_sendbuf[i_blk];
          Blk_List.message_noreschange_topsouthcorner_sendbuf[i_blk] = NULL;
       } /* endfor */
       delete []Blk_List.message_noreschange_topsouthcorner_sendbuf;
       Blk_List.message_noreschange_topsouthcorner_sendbuf = NULL;
    } /* endif */
    if (Blk_List.message_noreschange_topsouthcorner_recbuf != NULL) {
       for ( i_blk = 0; i_blk <= Blk_List.Nblk-1; ++i_blk ) {
          delete []Blk_List.message_noreschange_topsouthcorner_recbuf[i_blk];
          Blk_List.message_noreschange_topsouthcorner_recbuf[i_blk] = NULL;
       } /* endfor */
       delete []Blk_List.message_noreschange_topsouthcorner_recbuf;
       Blk_List.message_noreschange_topsouthcorner_recbuf = NULL;
    } /* endif */

    // Deallocate memory for  Top East corner send and receive buffers.
    if (Blk_List.message_noreschange_topeastcorner_sendbuf != NULL) {
       for ( i_blk = 0; i_blk <= Blk_List.Nblk-1; ++i_blk ) {
          delete []Blk_List.message_noreschange_topeastcorner_sendbuf[i_blk];
          Blk_List.message_noreschange_topeastcorner_sendbuf[i_blk] = NULL;
       } /* endfor */
       delete []Blk_List.message_noreschange_topeastcorner_sendbuf;
       Blk_List.message_noreschange_topeastcorner_sendbuf = NULL;
    } /* endif */
    if (Blk_List.message_noreschange_topeastcorner_recbuf != NULL) {
       for ( i_blk = 0; i_blk <= Blk_List.Nblk-1; ++i_blk ) {
          delete []Blk_List.message_noreschange_topeastcorner_recbuf[i_blk];
          Blk_List.message_noreschange_topeastcorner_recbuf[i_blk] = NULL;
       } /* endfor */
       delete []Blk_List.message_noreschange_topeastcorner_recbuf;
       Blk_List.message_noreschange_topeastcorner_recbuf = NULL;
    } /* endif */

    // Deallocate memory for  Top West corner send and receive buffers.
    if (Blk_List.message_noreschange_topwestcorner_sendbuf != NULL) {
       for ( i_blk = 0; i_blk <= Blk_List.Nblk-1; ++i_blk ) {
          delete []Blk_List.message_noreschange_topwestcorner_sendbuf[i_blk];
          Blk_List.message_noreschange_topwestcorner_sendbuf[i_blk] = NULL;
       } /* endfor */
       delete []Blk_List.message_noreschange_topwestcorner_sendbuf;
       Blk_List.message_noreschange_topwestcorner_sendbuf = NULL;
    } /* endif */
    if (Blk_List.message_noreschange_topwestcorner_recbuf != NULL) {
       for ( i_blk = 0; i_blk <= Blk_List.Nblk-1; ++i_blk ) {
          delete []Blk_List.message_noreschange_topwestcorner_recbuf[i_blk];
          Blk_List.message_noreschange_topwestcorner_recbuf[i_blk] = NULL;
       } /* endfor */
       delete []Blk_List.message_noreschange_topwestcorner_recbuf;
       Blk_List.message_noreschange_topwestcorner_recbuf = NULL;
    } /* endif */
    //*

    // Deallocate memory for Top North West corner send and receive buffers.
    if (Blk_List.message_noreschange_topnorthwestcorner_sendbuf != NULL) {
       for ( i_blk = 0; i_blk <= Blk_List.Nblk-1; ++i_blk ) {
          delete []Blk_List.message_noreschange_topnorthwestcorner_sendbuf[i_blk];
          Blk_List.message_noreschange_topnorthwestcorner_sendbuf[i_blk] = NULL;
       } /* endfor */
       delete []Blk_List.message_noreschange_topnorthwestcorner_sendbuf;
       Blk_List.message_noreschange_topnorthwestcorner_sendbuf = NULL;
    } /* endif */
    if (Blk_List.message_noreschange_topnorthwestcorner_recbuf != NULL) {
       for ( i_blk = 0; i_blk <= Blk_List.Nblk-1; ++i_blk ) {
          delete []Blk_List.message_noreschange_topnorthwestcorner_recbuf[i_blk];
          Blk_List.message_noreschange_topnorthwestcorner_recbuf[i_blk] = NULL;
       } /* endfor */
       delete []Blk_List.message_noreschange_topnorthwestcorner_recbuf;
       Blk_List.message_noreschange_topnorthwestcorner_recbuf = NULL;
    } /* endif */

    // Deallocate memory for Top North East corner send and receive buffers.
    if (Blk_List.message_noreschange_topnortheastcorner_sendbuf != NULL) {
       for ( i_blk = 0; i_blk <= Blk_List.Nblk-1; ++i_blk ) {
          delete []Blk_List.message_noreschange_topnortheastcorner_sendbuf[i_blk];
          Blk_List.message_noreschange_topnortheastcorner_sendbuf[i_blk] = NULL;
       } /* endfor */
       delete []Blk_List.message_noreschange_topnortheastcorner_sendbuf;
       Blk_List.message_noreschange_topnortheastcorner_sendbuf = NULL;
    } /* endif */
    if (Blk_List.message_noreschange_topnortheastcorner_recbuf != NULL) {
       for ( i_blk = 0; i_blk <= Blk_List.Nblk-1; ++i_blk ) {
          delete []Blk_List.message_noreschange_topnortheastcorner_recbuf[i_blk];
          Blk_List.message_noreschange_topnortheastcorner_recbuf[i_blk] = NULL;
       } /* endfor */
       delete []Blk_List.message_noreschange_topnortheastcorner_recbuf;
       Blk_List.message_noreschange_topnortheastcorner_recbuf = NULL;
    } /* endif */

    // Deallocate memory for Top South East corner send and receive buffers.
    if (Blk_List.message_noreschange_topsoutheastcorner_sendbuf != NULL) {
       for ( i_blk = 0; i_blk <= Blk_List.Nblk-1; ++i_blk ) {
          delete []Blk_List.message_noreschange_topsoutheastcorner_sendbuf[i_blk];
          Blk_List.message_noreschange_topsoutheastcorner_sendbuf[i_blk] = NULL;
       } /* endfor */
       delete []Blk_List.message_noreschange_topsoutheastcorner_sendbuf;
       Blk_List.message_noreschange_topsoutheastcorner_sendbuf = NULL;
    } /* endif */
    if (Blk_List.message_noreschange_topsoutheastcorner_recbuf != NULL) {
       for ( i_blk = 0; i_blk <= Blk_List.Nblk-1; ++i_blk ) {
          delete []Blk_List.message_noreschange_topsoutheastcorner_recbuf[i_blk];
          Blk_List.message_noreschange_topsoutheastcorner_recbuf[i_blk] = NULL;
       } /* endfor */
       delete []Blk_List.message_noreschange_topsoutheastcorner_recbuf;
       Blk_List.message_noreschange_topsoutheastcorner_recbuf = NULL;
    } /* endif */

    // Deallocate memory for Top South West corner send and receive buffers.
    if (Blk_List.message_noreschange_topsouthwestcorner_sendbuf != NULL) {
       for ( i_blk = 0; i_blk <= Blk_List.Nblk-1; ++i_blk ) {
          delete []Blk_List.message_noreschange_topsouthwestcorner_sendbuf[i_blk];
          Blk_List.message_noreschange_topsouthwestcorner_sendbuf[i_blk] = NULL;
       } /* endfor */
       delete []Blk_List.message_noreschange_topsouthwestcorner_sendbuf;
       Blk_List.message_noreschange_topsouthwestcorner_sendbuf = NULL;
    } /* endif */
    if (Blk_List.message_noreschange_topsouthwestcorner_recbuf != NULL) {
       for ( i_blk = 0; i_blk <= Blk_List.Nblk-1; ++i_blk ) {
          delete []Blk_List.message_noreschange_topsouthwestcorner_recbuf[i_blk];
          Blk_List.message_noreschange_topsouthwestcorner_recbuf[i_blk] = NULL;
       } /* endfor */
       delete []Blk_List.message_noreschange_topsouthwestcorner_recbuf;
       Blk_List.message_noreschange_topsouthwestcorner_recbuf = NULL;
    } /* endif */

    //*************//

     // Deallocate memory for Bottom North corner send and receive buffers.
    if (Blk_List.message_noreschange_bottomnorthcorner_sendbuf != NULL) {
       for ( i_blk = 0; i_blk <= Blk_List.Nblk-1; ++i_blk ) {
          delete []Blk_List.message_noreschange_bottomnorthcorner_sendbuf[i_blk];
          Blk_List.message_noreschange_bottomnorthcorner_sendbuf[i_blk] = NULL;
       } /* endfor */
       delete []Blk_List.message_noreschange_bottomnorthcorner_sendbuf;
       Blk_List.message_noreschange_bottomnorthcorner_sendbuf = NULL;
    } /* endif */
    if (Blk_List.message_noreschange_bottomnorthcorner_recbuf != NULL) {
       for ( i_blk = 0; i_blk <= Blk_List.Nblk-1; ++i_blk ) {
          delete []Blk_List.message_noreschange_bottomnorthcorner_recbuf[i_blk];
          Blk_List.message_noreschange_bottomnorthcorner_recbuf[i_blk] = NULL;
       } /* endfor */
       delete []Blk_List.message_noreschange_bottomnorthcorner_recbuf;
       Blk_List.message_noreschange_bottomnorthcorner_recbuf = NULL;
    } /* endif */

    // Deallocate memory for Bottom South  corner send and receive buffers.
    if (Blk_List.message_noreschange_bottomsouthcorner_sendbuf != NULL) {
       for ( i_blk = 0; i_blk <= Blk_List.Nblk-1; ++i_blk ) {
          delete []Blk_List.message_noreschange_bottomsouthcorner_sendbuf[i_blk];
          Blk_List.message_noreschange_bottomsouthcorner_sendbuf[i_blk] = NULL;
       } /* endfor */
       delete []Blk_List.message_noreschange_bottomsouthcorner_sendbuf;
       Blk_List.message_noreschange_bottomsouthcorner_sendbuf = NULL;
    } /* endif */
    if (Blk_List.message_noreschange_bottomsouthcorner_recbuf != NULL) {
       for ( i_blk = 0; i_blk <= Blk_List.Nblk-1; ++i_blk ) {
          delete []Blk_List.message_noreschange_bottomsouthcorner_recbuf[i_blk];
          Blk_List.message_noreschange_bottomsouthcorner_recbuf[i_blk] = NULL;
       } /* endfor */
       delete []Blk_List.message_noreschange_bottomsouthcorner_recbuf;
       Blk_List.message_noreschange_bottomsouthcorner_recbuf = NULL;
    } /* endif */

    // Deallocate memory for  Bottom East corner send and receive buffers.
    if (Blk_List.message_noreschange_bottomeastcorner_sendbuf != NULL) {
       for ( i_blk = 0; i_blk <= Blk_List.Nblk-1; ++i_blk ) {
          delete []Blk_List.message_noreschange_bottomeastcorner_sendbuf[i_blk];
          Blk_List.message_noreschange_bottomeastcorner_sendbuf[i_blk] = NULL;
       } /* endfor */
       delete []Blk_List.message_noreschange_bottomeastcorner_sendbuf;
       Blk_List.message_noreschange_bottomeastcorner_sendbuf = NULL;
    } /* endif */
    if (Blk_List.message_noreschange_bottomeastcorner_recbuf != NULL) {
       for ( i_blk = 0; i_blk <= Blk_List.Nblk-1; ++i_blk ) {
          delete []Blk_List.message_noreschange_bottomeastcorner_recbuf[i_blk];
          Blk_List.message_noreschange_bottomeastcorner_recbuf[i_blk] = NULL;
       } /* endfor */
       delete []Blk_List.message_noreschange_bottomeastcorner_recbuf;
       Blk_List.message_noreschange_bottomeastcorner_recbuf = NULL;
    } /* endif */

    // Deallocate memory for  Bottom West corner send and receive buffers.
    if (Blk_List.message_noreschange_bottomwestcorner_sendbuf != NULL) {
       for ( i_blk = 0; i_blk <= Blk_List.Nblk-1; ++i_blk ) {
          delete []Blk_List.message_noreschange_bottomwestcorner_sendbuf[i_blk];
          Blk_List.message_noreschange_bottomwestcorner_sendbuf[i_blk] = NULL;
       } /* endfor */
       delete []Blk_List.message_noreschange_bottomwestcorner_sendbuf;
       Blk_List.message_noreschange_bottomwestcorner_sendbuf = NULL;
    } /* endif */
    if (Blk_List.message_noreschange_bottomwestcorner_recbuf != NULL) {
       for ( i_blk = 0; i_blk <= Blk_List.Nblk-1; ++i_blk ) {
          delete []Blk_List.message_noreschange_bottomwestcorner_recbuf[i_blk];
          Blk_List.message_noreschange_bottomwestcorner_recbuf[i_blk] = NULL;
       } /* endfor */
       delete []Blk_List.message_noreschange_bottomwestcorner_recbuf;
       Blk_List.message_noreschange_bottomwestcorner_recbuf = NULL;
    } /* endif */
    //*

    // Deallocate memory for Bottom North West corner send and receive buffers.
    if (Blk_List.message_noreschange_bottomnorthwestcorner_sendbuf != NULL) {
       for ( i_blk = 0; i_blk <= Blk_List.Nblk-1; ++i_blk ) {
          delete []Blk_List.message_noreschange_bottomnorthwestcorner_sendbuf[i_blk];
          Blk_List.message_noreschange_bottomnorthwestcorner_sendbuf[i_blk] = NULL;
       } /* endfor */
       delete []Blk_List.message_noreschange_bottomnorthwestcorner_sendbuf;
       Blk_List.message_noreschange_bottomnorthwestcorner_sendbuf = NULL;
    } /* endif */
    if (Blk_List.message_noreschange_bottomnorthwestcorner_recbuf != NULL) {
       for ( i_blk = 0; i_blk <= Blk_List.Nblk-1; ++i_blk ) {
          delete []Blk_List.message_noreschange_bottomnorthwestcorner_recbuf[i_blk];
          Blk_List.message_noreschange_bottomnorthwestcorner_recbuf[i_blk] = NULL;
       } /* endfor */
       delete []Blk_List.message_noreschange_bottomnorthwestcorner_recbuf;
       Blk_List.message_noreschange_bottomnorthwestcorner_recbuf = NULL;
    } /* endif */

    // Deallocate memory for Bottom North East corner send and receive buffers.
    if (Blk_List.message_noreschange_bottomnortheastcorner_sendbuf != NULL) {
       for ( i_blk = 0; i_blk <= Blk_List.Nblk-1; ++i_blk ) {
          delete []Blk_List.message_noreschange_bottomnortheastcorner_sendbuf[i_blk];
          Blk_List.message_noreschange_bottomnortheastcorner_sendbuf[i_blk] = NULL;
       } /* endfor */
       delete []Blk_List.message_noreschange_bottomnortheastcorner_sendbuf;
       Blk_List.message_noreschange_bottomnortheastcorner_sendbuf = NULL;
    } /* endif */
    if (Blk_List.message_noreschange_bottomnortheastcorner_recbuf != NULL) {
       for ( i_blk = 0; i_blk <= Blk_List.Nblk-1; ++i_blk ) {
          delete []Blk_List.message_noreschange_bottomnortheastcorner_recbuf[i_blk];
          Blk_List.message_noreschange_bottomnortheastcorner_recbuf[i_blk] = NULL;
       } /* endfor */
       delete []Blk_List.message_noreschange_bottomnortheastcorner_recbuf;
       Blk_List.message_noreschange_bottomnortheastcorner_recbuf = NULL;
    } /* endif */

    // Deallocate memory for Bottom South East corner send and receive buffers.
    if (Blk_List.message_noreschange_bottomsoutheastcorner_sendbuf != NULL) {
       for ( i_blk = 0; i_blk <= Blk_List.Nblk-1; ++i_blk ) {
          delete []Blk_List.message_noreschange_bottomsoutheastcorner_sendbuf[i_blk];
          Blk_List.message_noreschange_bottomsoutheastcorner_sendbuf[i_blk] = NULL;
       } /* endfor */
       delete []Blk_List.message_noreschange_bottomsoutheastcorner_sendbuf;
       Blk_List.message_noreschange_bottomsoutheastcorner_sendbuf = NULL;
    } /* endif */
    if (Blk_List.message_noreschange_bottomsoutheastcorner_recbuf != NULL) {
       for ( i_blk = 0; i_blk <= Blk_List.Nblk-1; ++i_blk ) {
          delete []Blk_List.message_noreschange_bottomsoutheastcorner_recbuf[i_blk];
          Blk_List.message_noreschange_bottomsoutheastcorner_recbuf[i_blk] = NULL;
       } /* endfor */
       delete []Blk_List.message_noreschange_bottomsoutheastcorner_recbuf;
       Blk_List.message_noreschange_bottomsoutheastcorner_recbuf = NULL;
    } /* endif */

    // Deallocate memory for Bottom South West corner send and receive buffers.
    if (Blk_List.message_noreschange_bottomsouthwestcorner_sendbuf != NULL) {
       for ( i_blk = 0; i_blk <= Blk_List.Nblk-1; ++i_blk ) {
          delete []Blk_List.message_noreschange_bottomsouthwestcorner_sendbuf[i_blk];
          Blk_List.message_noreschange_bottomsouthwestcorner_sendbuf[i_blk] = NULL;
       } /* endfor */
       delete []Blk_List.message_noreschange_bottomsouthwestcorner_sendbuf;
       Blk_List.message_noreschange_bottomsouthwestcorner_sendbuf = NULL;
    } /* endif */
    if (Blk_List.message_noreschange_bottomsouthwestcorner_recbuf != NULL) {
       for ( i_blk = 0; i_blk <= Blk_List.Nblk-1; ++i_blk ) {
          delete []Blk_List.message_noreschange_bottomsouthwestcorner_recbuf[i_blk];
          Blk_List.message_noreschange_bottomsouthwestcorner_recbuf[i_blk] = NULL;
       } /* endfor */
       delete []Blk_List.message_noreschange_bottomsouthwestcorner_recbuf;
       Blk_List.message_noreschange_bottomsouthwestcorner_recbuf = NULL;
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

    // Deallocate memory for Top face send and receive buffers.
    if (Blk_List.message_reschange_topface_sendbuf != NULL) {
       for ( i_blk = 0; i_blk <= Blk_List.Nblk-1; ++i_blk ) {
          for (j_neigh = 0; j_neigh < 2 ; ++j_neigh) {
             delete []Blk_List.message_reschange_topface_sendbuf[i_blk][j_neigh];
             Blk_List.message_reschange_topface_sendbuf[i_blk][j_neigh] = NULL;
          } /* endfor */
          delete []Blk_List.message_reschange_topface_sendbuf[i_blk];
          Blk_List.message_reschange_topface_sendbuf[i_blk] = NULL;
       } /* endfor */
       delete []Blk_List.message_reschange_topface_sendbuf;
       Blk_List.message_reschange_topface_sendbuf = NULL;
    } /* endif */
    if (Blk_List.message_reschange_topface_recbuf != NULL) {
       for ( i_blk = 0; i_blk <= Blk_List.Nblk-1; ++i_blk ) {
          for (j_neigh = 0; j_neigh < 2 ; ++j_neigh) {
             delete []Blk_List.message_reschange_topface_recbuf[i_blk][j_neigh];
             Blk_List.message_reschange_topface_recbuf[i_blk][j_neigh] = NULL;
          } /* endfor */
          delete []Blk_List.message_reschange_topface_recbuf[i_blk];
          Blk_List.message_reschange_topface_recbuf[i_blk] = NULL;
       } /* endfor */
       delete []Blk_List.message_reschange_topface_recbuf;
       Blk_List.message_reschange_topface_recbuf = NULL;
    } /* endif */

    // Deallocate memory for Bottom face send and receive buffers.
    if (Blk_List.message_reschange_bottomface_sendbuf != NULL) {
       for ( i_blk = 0; i_blk <= Blk_List.Nblk-1; ++i_blk ) {
          for (j_neigh = 0; j_neigh < 2 ; ++j_neigh) {
             delete []Blk_List.message_reschange_bottomface_sendbuf[i_blk][j_neigh];
             Blk_List.message_reschange_bottomface_sendbuf[i_blk][j_neigh] = NULL;
          } /* endfor */
          delete []Blk_List.message_reschange_bottomface_sendbuf[i_blk];
          Blk_List.message_reschange_bottomface_sendbuf[i_blk] = NULL;
       } /* endfor */
       delete []Blk_List.message_reschange_bottomface_sendbuf;
       Blk_List.message_reschange_bottomface_sendbuf = NULL;
    } /* endif */
    if (Blk_List.message_reschange_bottomface_recbuf != NULL) {
       for ( i_blk = 0; i_blk <= Blk_List.Nblk-1; ++i_blk ) {
          for (j_neigh = 0; j_neigh < 2 ; ++j_neigh) {
             delete []Blk_List.message_reschange_bottomface_recbuf[i_blk][j_neigh];
             Blk_List.message_reschange_bottomface_recbuf[i_blk][j_neigh] = NULL;
          } /* endfor */
          delete []Blk_List.message_reschange_bottomface_recbuf[i_blk];
          Blk_List.message_reschange_bottomface_recbuf[i_blk] = NULL;
       } /* endfor */
       delete []Blk_List.message_reschange_bottomface_recbuf;
       Blk_List.message_reschange_bottomface_recbuf = NULL;
    } /* endif */

    // Deallocate memory for North face send and receive buffers.
    if (Blk_List.message_reschange_northface_sendbuf != NULL) {
       for ( i_blk = 0; i_blk <= Blk_List.Nblk-1; ++i_blk ) {
          for (j_neigh = 0; j_neigh < 2 ; ++j_neigh) {
             delete []Blk_List.message_reschange_northface_sendbuf[i_blk][j_neigh];
             Blk_List.message_reschange_northface_sendbuf[i_blk][j_neigh] = NULL;
          } /* endfor */
          delete []Blk_List.message_reschange_northface_sendbuf[i_blk];
          Blk_List.message_reschange_northface_sendbuf[i_blk] = NULL;
       } /* endfor */
       delete []Blk_List.message_reschange_northface_sendbuf;
       Blk_List.message_reschange_northface_sendbuf = NULL;
    } /* endif */
    if (Blk_List.message_reschange_northface_recbuf != NULL) {
       for ( i_blk = 0; i_blk <= Blk_List.Nblk-1; ++i_blk ) {
          for (j_neigh = 0; j_neigh < 2 ; ++j_neigh) {
             delete []Blk_List.message_reschange_northface_recbuf[i_blk][j_neigh];
             Blk_List.message_reschange_northface_recbuf[i_blk][j_neigh] = NULL;
          } /* endfor */
          delete []Blk_List.message_reschange_northface_recbuf[i_blk];
          Blk_List.message_reschange_northface_recbuf[i_blk] = NULL;
       } /* endfor */
       delete []Blk_List.message_reschange_northface_recbuf;
       Blk_List.message_reschange_northface_recbuf = NULL;
    } /* endif */

    // Deallocate memory for South face send and receive buffers.
    if (Blk_List.message_reschange_southface_sendbuf != NULL) {
       for ( i_blk = 0; i_blk <= Blk_List.Nblk-1; ++i_blk ) {
          for (j_neigh = 0; j_neigh < 2 ; ++j_neigh) {
             delete []Blk_List.message_reschange_southface_sendbuf[i_blk][j_neigh];
             Blk_List.message_reschange_southface_sendbuf[i_blk][j_neigh] = NULL;
          } /* endfor */
          delete []Blk_List.message_reschange_southface_sendbuf[i_blk];
          Blk_List.message_reschange_southface_sendbuf[i_blk] = NULL;
       } /* endfor */
       delete []Blk_List.message_reschange_southface_sendbuf;
       Blk_List.message_reschange_southface_sendbuf = NULL;
    } /* endif */
    if (Blk_List.message_reschange_southface_recbuf != NULL) {
       for ( i_blk = 0; i_blk <= Blk_List.Nblk-1; ++i_blk ) {
          for (j_neigh = 0; j_neigh < 2 ; ++j_neigh) {
             delete []Blk_List.message_reschange_southface_recbuf[i_blk][j_neigh];
             Blk_List.message_reschange_southface_recbuf[i_blk][j_neigh] = NULL;
          } /* endfor */
          delete []Blk_List.message_reschange_southface_recbuf[i_blk];
          Blk_List.message_reschange_southface_recbuf[i_blk] = NULL;
       } /* endfor */
       delete []Blk_List.message_reschange_southface_recbuf;
       Blk_List.message_reschange_southface_recbuf = NULL;
    } /* endif */

    // Deallocate memory for East face send and receive buffers.
    if (Blk_List.message_reschange_eastface_sendbuf != NULL) {
       for ( i_blk = 0; i_blk <= Blk_List.Nblk-1; ++i_blk ) {
          for (j_neigh = 0; j_neigh < 2 ; ++j_neigh) {
             delete []Blk_List.message_reschange_eastface_sendbuf[i_blk][j_neigh];
             Blk_List.message_reschange_eastface_sendbuf[i_blk][j_neigh] = NULL;
          } /* endfor */
          delete []Blk_List.message_reschange_eastface_sendbuf[i_blk];
          Blk_List.message_reschange_eastface_sendbuf[i_blk] = NULL;
       } /* endfor */
       delete []Blk_List.message_reschange_eastface_sendbuf;
       Blk_List.message_reschange_eastface_sendbuf = NULL;
    } /* endif */
    if (Blk_List.message_reschange_eastface_recbuf != NULL) {
       for ( i_blk = 0; i_blk <= Blk_List.Nblk-1; ++i_blk ) {
          for (j_neigh = 0; j_neigh < 2 ; ++j_neigh) {
             delete []Blk_List.message_reschange_eastface_recbuf[i_blk][j_neigh];
             Blk_List.message_reschange_eastface_recbuf[i_blk][j_neigh] = NULL;
          } /* endfor */
          delete []Blk_List.message_reschange_eastface_recbuf[i_blk];
          Blk_List.message_reschange_eastface_recbuf[i_blk] = NULL;
       } /* endfor */
       delete []Blk_List.message_reschange_eastface_recbuf;
       Blk_List.message_reschange_eastface_recbuf = NULL;
    } /* endif */

    // Deallocate memory for West face send and receive buffers.
    if (Blk_List.message_reschange_westface_sendbuf != NULL) {
       for ( i_blk = 0; i_blk <= Blk_List.Nblk-1; ++i_blk ) {
          for (j_neigh = 0; j_neigh < 2 ; ++j_neigh) {
             delete []Blk_List.message_reschange_westface_sendbuf[i_blk][j_neigh];
             Blk_List.message_reschange_westface_sendbuf[i_blk][j_neigh] = NULL;
          } /* endfor */
          delete []Blk_List.message_reschange_westface_sendbuf[i_blk];
          Blk_List.message_reschange_westface_sendbuf[i_blk] = NULL;
       } /* endfor */
       delete []Blk_List.message_reschange_westface_sendbuf;
       Blk_List.message_reschange_westface_sendbuf = NULL;
    } /* endif */
    if (Blk_List.message_reschange_westface_recbuf != NULL) {
       for ( i_blk = 0; i_blk <= Blk_List.Nblk-1; ++i_blk ) {
          for (j_neigh = 0; j_neigh < 2 ; ++j_neigh) {
             delete []Blk_List.message_reschange_westface_recbuf[i_blk][j_neigh];
             Blk_List.message_reschange_westface_recbuf[i_blk][j_neigh] = NULL;
          } /* endfor */
          delete []Blk_List.message_reschange_westface_recbuf[i_blk];
          Blk_List.message_reschange_westface_recbuf[i_blk] = NULL;
       } /* endfor */
       delete []Blk_List.message_reschange_westface_recbuf;
       Blk_List.message_reschange_westface_recbuf = NULL;
    } /* endif */

    // Deallocate memory for North West corner send and receive buffers.
    if (Blk_List.message_reschange_northwestcorner_sendbuf != NULL) {
       for ( i_blk = 0; i_blk <= Blk_List.Nblk-1; ++i_blk ) {
          delete []Blk_List.message_reschange_northwestcorner_sendbuf[i_blk];
          Blk_List.message_reschange_northwestcorner_sendbuf[i_blk] = NULL;
       } /* endfor */
       delete []Blk_List.message_reschange_northwestcorner_sendbuf;
       Blk_List.message_reschange_northwestcorner_sendbuf = NULL;
    } /* endif */
    if (Blk_List.message_reschange_northwestcorner_recbuf != NULL) {
       for ( i_blk = 0; i_blk <= Blk_List.Nblk-1; ++i_blk ) {
          delete []Blk_List.message_reschange_northwestcorner_recbuf[i_blk];
          Blk_List.message_reschange_northwestcorner_recbuf[i_blk] = NULL;
       } /* endfor */
       delete []Blk_List.message_reschange_northwestcorner_recbuf;
       Blk_List.message_reschange_northwestcorner_recbuf = NULL;
    } /* endif */

    // Deallocate memory for North East corner send and receive buffers.
    if (Blk_List.message_reschange_northeastcorner_sendbuf != NULL) {
       for ( i_blk = 0; i_blk <= Blk_List.Nblk-1; ++i_blk ) {
          delete []Blk_List.message_reschange_northeastcorner_sendbuf[i_blk];
          Blk_List.message_reschange_northeastcorner_sendbuf[i_blk] = NULL;
       } /* endfor */
       delete []Blk_List.message_reschange_northeastcorner_sendbuf;
       Blk_List.message_reschange_northeastcorner_sendbuf = NULL;
    } /* endif */
    if (Blk_List.message_reschange_northeastcorner_recbuf != NULL) {
       for ( i_blk = 0; i_blk <= Blk_List.Nblk-1; ++i_blk ) {
          delete []Blk_List.message_reschange_northeastcorner_recbuf[i_blk];
          Blk_List.message_reschange_northeastcorner_recbuf[i_blk] = NULL;
       } /* endfor */
       delete []Blk_List.message_reschange_northeastcorner_recbuf;
       Blk_List.message_reschange_northeastcorner_recbuf = NULL;
    } /* endif */

    // Deallocate memory for South East corner send and receive buffers.
    if (Blk_List.message_reschange_southeastcorner_sendbuf != NULL) {
       for ( i_blk = 0; i_blk <= Blk_List.Nblk-1; ++i_blk ) {
          delete []Blk_List.message_reschange_southeastcorner_sendbuf[i_blk];
          Blk_List.message_reschange_southeastcorner_sendbuf[i_blk] = NULL;
       } /* endfor */
       delete []Blk_List.message_reschange_southeastcorner_sendbuf;
       Blk_List.message_reschange_southeastcorner_sendbuf = NULL;
    } /* endif */
    if (Blk_List.message_reschange_southeastcorner_recbuf != NULL) {
       for ( i_blk = 0; i_blk <= Blk_List.Nblk-1; ++i_blk ) {
          delete []Blk_List.message_reschange_southeastcorner_recbuf[i_blk];
          Blk_List.message_reschange_southeastcorner_recbuf[i_blk] = NULL;
       } /* endfor */
       delete []Blk_List.message_reschange_southeastcorner_recbuf;
       Blk_List.message_reschange_southeastcorner_recbuf = NULL;
    } /* endif */

    // Deallocate memory for South West corner send and receive buffers.
    if (Blk_List.message_reschange_southwestcorner_sendbuf != NULL) {
       for ( i_blk = 0; i_blk <= Blk_List.Nblk-1; ++i_blk ) {
          delete []Blk_List.message_reschange_southwestcorner_sendbuf[i_blk];
          Blk_List.message_reschange_southwestcorner_sendbuf[i_blk] = NULL;
       } /* endfor */
       delete []Blk_List.message_reschange_southwestcorner_sendbuf;
       Blk_List.message_reschange_southwestcorner_sendbuf = NULL;
    } /* endif */
    if (Blk_List.message_reschange_southwestcorner_recbuf != NULL) {
       for ( i_blk = 0; i_blk <= Blk_List.Nblk-1; ++i_blk ) {
          delete []Blk_List.message_reschange_southwestcorner_recbuf[i_blk];
          Blk_List.message_reschange_southwestcorner_recbuf[i_blk] = NULL;
       } /* endfor */
       delete []Blk_List.message_reschange_southwestcorner_recbuf;
       Blk_List.message_reschange_southwestcorner_recbuf = NULL;
    } /* endif */

    //**********/
    // Deallocate memory for __ North corner send and receive buffers.
    if (Blk_List.message_reschange_topnorthcorner_sendbuf != NULL) {
       for ( i_blk = 0; i_blk <= Blk_List.Nblk-1; ++i_blk ) {
          delete []Blk_List.message_reschange_topnorthcorner_sendbuf[i_blk];
          Blk_List.message_reschange_topnorthcorner_sendbuf[i_blk] = NULL;
       } /* endfor */
       delete []Blk_List.message_reschange_topnorthcorner_sendbuf;
       Blk_List.message_reschange_topnorthcorner_sendbuf = NULL;
    } /* endif */
    if (Blk_List.message_reschange_topnorthcorner_recbuf != NULL) {
       for ( i_blk = 0; i_blk <= Blk_List.Nblk-1; ++i_blk ) {
          delete []Blk_List.message_reschange_topnorthcorner_recbuf[i_blk];
          Blk_List.message_reschange_topnorthcorner_recbuf[i_blk] = NULL;
       } /* endfor */
       delete []Blk_List.message_reschange_topnorthcorner_recbuf;
       Blk_List.message_reschange_topnorthcorner_recbuf = NULL;
    } /* endif */

    // Deallocate memory for __ South corner send and receive buffers.
    if (Blk_List.message_reschange_topsouthcorner_sendbuf != NULL) {
       for ( i_blk = 0; i_blk <= Blk_List.Nblk-1; ++i_blk ) {
          delete []Blk_List.message_reschange_topsouthcorner_sendbuf[i_blk];
          Blk_List.message_reschange_topsouthcorner_sendbuf[i_blk] = NULL;
       } /* endfor */
       delete []Blk_List.message_reschange_topsouthcorner_sendbuf;
       Blk_List.message_reschange_topsouthcorner_sendbuf = NULL;
    } /* endif */
    if (Blk_List.message_reschange_topsouthcorner_recbuf != NULL) {
       for ( i_blk = 0; i_blk <= Blk_List.Nblk-1; ++i_blk ) {
          delete []Blk_List.message_reschange_topsouthcorner_recbuf[i_blk];
          Blk_List.message_reschange_topsouthcorner_recbuf[i_blk] = NULL;
       } /* endfor */
       delete []Blk_List.message_reschange_topsouthcorner_recbuf;
       Blk_List.message_reschange_topsouthcorner_recbuf = NULL;
    } /* endif */

    // Deallocate memory for __ East corner send and receive buffers.
    if (Blk_List.message_reschange_topeastcorner_sendbuf != NULL) {
       for ( i_blk = 0; i_blk <= Blk_List.Nblk-1; ++i_blk ) {
          delete []Blk_List.message_reschange_topeastcorner_sendbuf[i_blk];
          Blk_List.message_reschange_topeastcorner_sendbuf[i_blk] = NULL;
       } /* endfor */
       delete []Blk_List.message_reschange_topeastcorner_sendbuf;
       Blk_List.message_reschange_topeastcorner_sendbuf = NULL;
    } /* endif */
    if (Blk_List.message_reschange_topeastcorner_recbuf != NULL) {
       for ( i_blk = 0; i_blk <= Blk_List.Nblk-1; ++i_blk ) {
          delete []Blk_List.message_reschange_topeastcorner_recbuf[i_blk];
          Blk_List.message_reschange_topeastcorner_recbuf[i_blk] = NULL;
       } /* endfor */
       delete []Blk_List.message_reschange_topeastcorner_recbuf;
       Blk_List.message_reschange_topeastcorner_recbuf = NULL;
    } /* endif */

    // Deallocate memory for __ West corner send and receive buffers.
    if (Blk_List.message_reschange_topwestcorner_sendbuf != NULL) {
       for ( i_blk = 0; i_blk <= Blk_List.Nblk-1; ++i_blk ) {
          delete []Blk_List.message_reschange_topwestcorner_sendbuf[i_blk];
          Blk_List.message_reschange_topwestcorner_sendbuf[i_blk] = NULL;
       } /* endfor */
       delete []Blk_List.message_reschange_topwestcorner_sendbuf;
       Blk_List.message_reschange_topwestcorner_sendbuf = NULL;
    } /* endif */
    if (Blk_List.message_reschange_topwestcorner_recbuf != NULL) {
       for ( i_blk = 0; i_blk <= Blk_List.Nblk-1; ++i_blk ) {
          delete []Blk_List.message_reschange_topwestcorner_recbuf[i_blk];
          Blk_List.message_reschange_topwestcorner_recbuf[i_blk] = NULL;
       } /* endfor */
       delete []Blk_List.message_reschange_topwestcorner_recbuf;
       Blk_List.message_reschange_topwestcorner_recbuf = NULL;
    } /* endif */


    // Deallocate memory for North West corner send and receive buffers.
    if (Blk_List.message_reschange_topnorthwestcorner_sendbuf != NULL) {
       for ( i_blk = 0; i_blk <= Blk_List.Nblk-1; ++i_blk ) {
          delete []Blk_List.message_reschange_topnorthwestcorner_sendbuf[i_blk];
          Blk_List.message_reschange_topnorthwestcorner_sendbuf[i_blk] = NULL;
       } /* endfor */
       delete []Blk_List.message_reschange_topnorthwestcorner_sendbuf;
       Blk_List.message_reschange_topnorthwestcorner_sendbuf = NULL;
    } /* endif */
    if (Blk_List.message_reschange_topnorthwestcorner_recbuf != NULL) {
       for ( i_blk = 0; i_blk <= Blk_List.Nblk-1; ++i_blk ) {
          delete []Blk_List.message_reschange_topnorthwestcorner_recbuf[i_blk];
          Blk_List.message_reschange_topnorthwestcorner_recbuf[i_blk] = NULL;
       } /* endfor */
       delete []Blk_List.message_reschange_topnorthwestcorner_recbuf;
       Blk_List.message_reschange_topnorthwestcorner_recbuf = NULL;
    } /* endif */

    // Deallocate memory for North East corner send and receive buffers.
    if (Blk_List.message_reschange_topnortheastcorner_sendbuf != NULL) {
       for ( i_blk = 0; i_blk <= Blk_List.Nblk-1; ++i_blk ) {
          delete []Blk_List.message_reschange_topnortheastcorner_sendbuf[i_blk];
          Blk_List.message_reschange_topnortheastcorner_sendbuf[i_blk] = NULL;
       } /* endfor */
       delete []Blk_List.message_reschange_topnortheastcorner_sendbuf;
       Blk_List.message_reschange_topnortheastcorner_sendbuf = NULL;
    } /* endif */
    if (Blk_List.message_reschange_topnortheastcorner_recbuf != NULL) {
       for ( i_blk = 0; i_blk <= Blk_List.Nblk-1; ++i_blk ) {
          delete []Blk_List.message_reschange_topnortheastcorner_recbuf[i_blk];
          Blk_List.message_reschange_topnortheastcorner_recbuf[i_blk] = NULL;
       } /* endfor */
       delete []Blk_List.message_reschange_topnortheastcorner_recbuf;
       Blk_List.message_reschange_topnortheastcorner_recbuf = NULL;
    } /* endif */

    // Deallocate memory for South East corner send and receive buffers.
    if (Blk_List.message_reschange_topsoutheastcorner_sendbuf != NULL) {
       for ( i_blk = 0; i_blk <= Blk_List.Nblk-1; ++i_blk ) {
          delete []Blk_List.message_reschange_topsoutheastcorner_sendbuf[i_blk];
          Blk_List.message_reschange_topsoutheastcorner_sendbuf[i_blk] = NULL;
       } /* endfor */
       delete []Blk_List.message_reschange_topsoutheastcorner_sendbuf;
       Blk_List.message_reschange_topsoutheastcorner_sendbuf = NULL;
    } /* endif */
    if (Blk_List.message_reschange_topsoutheastcorner_recbuf != NULL) {
       for ( i_blk = 0; i_blk <= Blk_List.Nblk-1; ++i_blk ) {
          delete []Blk_List.message_reschange_topsoutheastcorner_recbuf[i_blk];
          Blk_List.message_reschange_topsoutheastcorner_recbuf[i_blk] = NULL;
       } /* endfor */
       delete []Blk_List.message_reschange_topsoutheastcorner_recbuf;
       Blk_List.message_reschange_topsoutheastcorner_recbuf = NULL;
    } /* endif */

    // Deallocate memory for South West corner send and receive buffers.
    if (Blk_List.message_reschange_topsouthwestcorner_sendbuf != NULL) {
       for ( i_blk = 0; i_blk <= Blk_List.Nblk-1; ++i_blk ) {
          delete []Blk_List.message_reschange_topsouthwestcorner_sendbuf[i_blk];
          Blk_List.message_reschange_topsouthwestcorner_sendbuf[i_blk] = NULL;
       } /* endfor */
       delete []Blk_List.message_reschange_topsouthwestcorner_sendbuf;
       Blk_List.message_reschange_topsouthwestcorner_sendbuf = NULL;
    } /* endif */
    if (Blk_List.message_reschange_topsouthwestcorner_recbuf != NULL) {
       for ( i_blk = 0; i_blk <= Blk_List.Nblk-1; ++i_blk ) {
          delete []Blk_List.message_reschange_topsouthwestcorner_recbuf[i_blk];
          Blk_List.message_reschange_topsouthwestcorner_recbuf[i_blk] = NULL;
       } /* endfor */
       delete []Blk_List.message_reschange_topsouthwestcorner_recbuf;
       Blk_List.message_reschange_topsouthwestcorner_recbuf = NULL;
    } /* endif */

    //**********/
    // Deallocate memory for __ North corner send and receive buffers.
    if (Blk_List.message_reschange_bottomnorthcorner_sendbuf != NULL) {
       for ( i_blk = 0; i_blk <= Blk_List.Nblk-1; ++i_blk ) {
          delete []Blk_List.message_reschange_bottomnorthcorner_sendbuf[i_blk];
          Blk_List.message_reschange_bottomnorthcorner_sendbuf[i_blk] = NULL;
       } /* endfor */
       delete []Blk_List.message_reschange_bottomnorthcorner_sendbuf;
       Blk_List.message_reschange_bottomnorthcorner_sendbuf = NULL;
    } /* endif */
    if (Blk_List.message_reschange_bottomnorthcorner_recbuf != NULL) {
       for ( i_blk = 0; i_blk <= Blk_List.Nblk-1; ++i_blk ) {
          delete []Blk_List.message_reschange_bottomnorthcorner_recbuf[i_blk];
          Blk_List.message_reschange_bottomnorthcorner_recbuf[i_blk] = NULL;
       } /* endfor */
       delete []Blk_List.message_reschange_bottomnorthcorner_recbuf;
       Blk_List.message_reschange_bottomnorthcorner_recbuf = NULL;
    } /* endif */

    // Deallocate memory for __ South corner send and receive buffers.
    if (Blk_List.message_reschange_bottomsouthcorner_sendbuf != NULL) {
       for ( i_blk = 0; i_blk <= Blk_List.Nblk-1; ++i_blk ) {
          delete []Blk_List.message_reschange_bottomsouthcorner_sendbuf[i_blk];
          Blk_List.message_reschange_bottomsouthcorner_sendbuf[i_blk] = NULL;
       } /* endfor */
       delete []Blk_List.message_reschange_bottomsouthcorner_sendbuf;
       Blk_List.message_reschange_bottomsouthcorner_sendbuf = NULL;
    } /* endif */
    if (Blk_List.message_reschange_bottomsouthcorner_recbuf != NULL) {
       for ( i_blk = 0; i_blk <= Blk_List.Nblk-1; ++i_blk ) {
          delete []Blk_List.message_reschange_bottomsouthcorner_recbuf[i_blk];
          Blk_List.message_reschange_bottomsouthcorner_recbuf[i_blk] = NULL;
       } /* endfor */
       delete []Blk_List.message_reschange_bottomsouthcorner_recbuf;
       Blk_List.message_reschange_bottomsouthcorner_recbuf = NULL;
    } /* endif */

    // Deallocate memory for __ East corner send and receive buffers.
    if (Blk_List.message_reschange_bottomeastcorner_sendbuf != NULL) {
       for ( i_blk = 0; i_blk <= Blk_List.Nblk-1; ++i_blk ) {
          delete []Blk_List.message_reschange_bottomeastcorner_sendbuf[i_blk];
          Blk_List.message_reschange_bottomeastcorner_sendbuf[i_blk] = NULL;
       } /* endfor */
       delete []Blk_List.message_reschange_bottomeastcorner_sendbuf;
       Blk_List.message_reschange_bottomeastcorner_sendbuf = NULL;
    } /* endif */
    if (Blk_List.message_reschange_bottomeastcorner_recbuf != NULL) {
       for ( i_blk = 0; i_blk <= Blk_List.Nblk-1; ++i_blk ) {
          delete []Blk_List.message_reschange_bottomeastcorner_recbuf[i_blk];
          Blk_List.message_reschange_bottomeastcorner_recbuf[i_blk] = NULL;
       } /* endfor */
       delete []Blk_List.message_reschange_bottomeastcorner_recbuf;
       Blk_List.message_reschange_bottomeastcorner_recbuf = NULL;
    } /* endif */

    // Deallocate memory for __ West corner send and receive buffers.
    if (Blk_List.message_reschange_bottomwestcorner_sendbuf != NULL) {
       for ( i_blk = 0; i_blk <= Blk_List.Nblk-1; ++i_blk ) {
          delete []Blk_List.message_reschange_bottomwestcorner_sendbuf[i_blk];
          Blk_List.message_reschange_bottomwestcorner_sendbuf[i_blk] = NULL;
       } /* endfor */
       delete []Blk_List.message_reschange_bottomwestcorner_sendbuf;
       Blk_List.message_reschange_bottomwestcorner_sendbuf = NULL;
    } /* endif */
    if (Blk_List.message_reschange_bottomwestcorner_recbuf != NULL) {
       for ( i_blk = 0; i_blk <= Blk_List.Nblk-1; ++i_blk ) {
          delete []Blk_List.message_reschange_bottomwestcorner_recbuf[i_blk];
          Blk_List.message_reschange_bottomwestcorner_recbuf[i_blk] = NULL;
       } /* endfor */
       delete []Blk_List.message_reschange_bottomwestcorner_recbuf;
       Blk_List.message_reschange_bottomwestcorner_recbuf = NULL;
    } /* endif */


    // Deallocate memory for North West corner send and receive buffers.
    if (Blk_List.message_reschange_bottomnorthwestcorner_sendbuf != NULL) {
       for ( i_blk = 0; i_blk <= Blk_List.Nblk-1; ++i_blk ) {
          delete []Blk_List.message_reschange_bottomnorthwestcorner_sendbuf[i_blk];
          Blk_List.message_reschange_bottomnorthwestcorner_sendbuf[i_blk] = NULL;
       } /* endfor */
       delete []Blk_List.message_reschange_bottomnorthwestcorner_sendbuf;
       Blk_List.message_reschange_bottomnorthwestcorner_sendbuf = NULL;
    } /* endif */
    if (Blk_List.message_reschange_bottomnorthwestcorner_recbuf != NULL) {
       for ( i_blk = 0; i_blk <= Blk_List.Nblk-1; ++i_blk ) {
          delete []Blk_List.message_reschange_bottomnorthwestcorner_recbuf[i_blk];
          Blk_List.message_reschange_bottomnorthwestcorner_recbuf[i_blk] = NULL;
       } /* endfor */
       delete []Blk_List.message_reschange_bottomnorthwestcorner_recbuf;
       Blk_List.message_reschange_bottomnorthwestcorner_recbuf = NULL;
    } /* endif */

    // Deallocate memory for North East corner send and receive buffers.
    if (Blk_List.message_reschange_bottomnortheastcorner_sendbuf != NULL) {
       for ( i_blk = 0; i_blk <= Blk_List.Nblk-1; ++i_blk ) {
          delete []Blk_List.message_reschange_bottomnortheastcorner_sendbuf[i_blk];
          Blk_List.message_reschange_bottomnortheastcorner_sendbuf[i_blk] = NULL;
       } /* endfor */
       delete []Blk_List.message_reschange_bottomnortheastcorner_sendbuf;
       Blk_List.message_reschange_bottomnortheastcorner_sendbuf = NULL;
    } /* endif */
    if (Blk_List.message_reschange_bottomnortheastcorner_recbuf != NULL) {
       for ( i_blk = 0; i_blk <= Blk_List.Nblk-1; ++i_blk ) {
          delete []Blk_List.message_reschange_bottomnortheastcorner_recbuf[i_blk];
          Blk_List.message_reschange_bottomnortheastcorner_recbuf[i_blk] = NULL;
       } /* endfor */
       delete []Blk_List.message_reschange_bottomnortheastcorner_recbuf;
       Blk_List.message_reschange_bottomnortheastcorner_recbuf = NULL;
    } /* endif */

    // Deallocate memory for South East corner send and receive buffers.
    if (Blk_List.message_reschange_bottomsoutheastcorner_sendbuf != NULL) {
       for ( i_blk = 0; i_blk <= Blk_List.Nblk-1; ++i_blk ) {
          delete []Blk_List.message_reschange_bottomsoutheastcorner_sendbuf[i_blk];
          Blk_List.message_reschange_bottomsoutheastcorner_sendbuf[i_blk] = NULL;
       } /* endfor */
       delete []Blk_List.message_reschange_bottomsoutheastcorner_sendbuf;
       Blk_List.message_reschange_bottomsoutheastcorner_sendbuf = NULL;
    } /* endif */
    if (Blk_List.message_reschange_bottomsoutheastcorner_recbuf != NULL) {
       for ( i_blk = 0; i_blk <= Blk_List.Nblk-1; ++i_blk ) {
          delete []Blk_List.message_reschange_bottomsoutheastcorner_recbuf[i_blk];
          Blk_List.message_reschange_bottomsoutheastcorner_recbuf[i_blk] = NULL;
       } /* endfor */
       delete []Blk_List.message_reschange_bottomsoutheastcorner_recbuf;
       Blk_List.message_reschange_bottomsoutheastcorner_recbuf = NULL;
    } /* endif */

    // Deallocate memory for South West corner send and receive buffers.
    if (Blk_List.message_reschange_bottomsouthwestcorner_sendbuf != NULL) {
       for ( i_blk = 0; i_blk <= Blk_List.Nblk-1; ++i_blk ) {
          delete []Blk_List.message_reschange_bottomsouthwestcorner_sendbuf[i_blk];
          Blk_List.message_reschange_bottomsouthwestcorner_sendbuf[i_blk] = NULL;
       } /* endfor */
       delete []Blk_List.message_reschange_bottomsouthwestcorner_sendbuf;
       Blk_List.message_reschange_bottomsouthwestcorner_sendbuf = NULL;
    } /* endif */
    if (Blk_List.message_reschange_bottomsouthwestcorner_recbuf != NULL) {
       for ( i_blk = 0; i_blk <= Blk_List.Nblk-1; ++i_blk ) {
          delete []Blk_List.message_reschange_bottomsouthwestcorner_recbuf[i_blk];
          Blk_List.message_reschange_bottomsouthwestcorner_recbuf[i_blk] = NULL;
       } /* endfor */
       delete []Blk_List.message_reschange_bottomsouthwestcorner_recbuf;
       Blk_List.message_reschange_bottomsouthwestcorner_recbuf = NULL;
    } /* endif */

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
                      const int Number_of_Solution_Variables) {

    int error_flag;

    /* Exchange message buffers at block interfaces 
       with no cell resolution change. */

    error_flag = Exchange_Messages_NoResChange(Blk_List,
                                               Number_of_Solution_Variables);
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
