/* AdaptiveBlock.cc:  Subroutines for adaptive blocks classes
                      that are independent of the CPU architecture. */

/* Include adaptive block header file. */

#ifndef _ADAPTIVBLOCK_INCLUDED
#include "AdaptiveBlock.h"
#endif // _ADAPTIVEBLOCK_INCLUDED

/*************************************************************
 * AdaptiveBlockResourceList -- External subroutines.        *
 *************************************************************/

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
void Create_Block_Resource_List(AdaptiveBlockResourceList &List_of_Available_Blocks,
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

/*************************************************************
 * AdaptiveBlock2D_Info -- External subroutines.             *
 *************************************************************/

/********************************************************
 * Routine: Broadcast_Adaptive_Block_Info               *
 *                                                      *
 * Broadcast adaptive block information to all          *
 * processors involved in the calculation from the      *
 * primary processor using the MPI broadcast routine.   *
 *                                                      *
 ********************************************************/
void Broadcast_Adaptive_Block_Info(AdaptiveBlock2D_Info &Blk_Info) {

#ifdef _MPI_VERSION
  MPI::COMM_WORLD.Bcast(&(Blk_Info.cpu), 1, MPI::INT, 0);
  MPI::COMM_WORLD.Bcast(&(Blk_Info.blknum), 1, MPI::INT, 0);
  MPI::COMM_WORLD.Bcast(&(Blk_Info.dimen.i), 1, MPI::INT, 0);
  MPI::COMM_WORLD.Bcast(&(Blk_Info.dimen.j), 1, MPI::INT, 0);
  MPI::COMM_WORLD.Bcast(&(Blk_Info.dimen.ghost), 1, MPI::INT, 0);
  MPI::COMM_WORLD.Bcast(&(Blk_Info.sector), 1, MPI::INT, 0);
  MPI::COMM_WORLD.Bcast(&(Blk_Info.level), 1, MPI::INT, 0);
#endif

}

/*************************************************************
 * AdaptiveBlock2D -- External subroutines.                  *
 *************************************************************/

/********************************************************
 * Routine: Broadcast_Adaptive_Block                    *
 *                                                      *
 * Broadcast adaptive block information to all          *
 * processors involved in the calculation from the      *
 * primary processor using the MPI broadcast routine.   *
 *                                                      *
 ********************************************************/
void Broadcast_Adaptive_Block(AdaptiveBlock2D &Blk) {

#ifdef _MPI_VERSION
  int i;
  MPI::COMM_WORLD.Bcast(&(Blk.used), 1, MPI::INT, 0);
  MPI::COMM_WORLD.Bcast(&(Blk.gblknum), 1, MPI::INT, 0);
  Broadcast_Adaptive_Block_Info(Blk.info);
  MPI::COMM_WORLD.Bcast(&(Blk.nN), 1, MPI::INT, 0);
  for (i = 0; i <= Blk.nN-1; ++i) {
     Broadcast_Adaptive_Block_Info(Blk.infoN[i]);
  } /* endfor */
  MPI::COMM_WORLD.Bcast(&(Blk.nS), 1, MPI::INT, 0);
  for (i = 0; i <= Blk.nS-1; ++i) {
     Broadcast_Adaptive_Block_Info(Blk.infoS[i]);
  } /* endfor */
  MPI::COMM_WORLD.Bcast(&(Blk.nE), 1, MPI::INT, 0);
  for (i = 0; i <= Blk.nE-1; ++i) {
     Broadcast_Adaptive_Block_Info(Blk.infoE[i]);
  } /* endfor */
  MPI::COMM_WORLD.Bcast(&(Blk.nW), 1, MPI::INT, 0);
  for (i = 0; i <= Blk.nW-1; ++i) {
     Broadcast_Adaptive_Block_Info(Blk.infoW[i]);
  } /* endfor */
  MPI::COMM_WORLD.Bcast(&(Blk.nNW), 1, MPI::INT, 0);
  for (i = 0; i <= Blk.nNW-1; ++i) {
     Broadcast_Adaptive_Block_Info(Blk.infoNW[i]);
  } /* endfor */
  MPI::COMM_WORLD.Bcast(&(Blk.nNE), 1, MPI::INT, 0);
  for (i = 0; i <= Blk.nNE-1; ++i) {
     Broadcast_Adaptive_Block_Info(Blk.infoNE[i]);
  } /* endfor */
  MPI::COMM_WORLD.Bcast(&(Blk.nSE), 1, MPI::INT, 0);
  for (i = 0; i <= Blk.nSE-1; ++i) {
     Broadcast_Adaptive_Block_Info(Blk.infoSE[i]);
  } /* endfor */
  MPI::COMM_WORLD.Bcast(&(Blk.nSW), 1, MPI::INT, 0);
  for (i = 0; i <= Blk.nSW-1; ++i) {
     Broadcast_Adaptive_Block_Info(Blk.infoSW[i]);
  } /* endfor */
#endif

}

/*************************************************************
 * AdaptiveBlock2D_List -- External subroutines.             *
 *************************************************************/

/**********************************************************
 * Routine: Allocate_Message_Buffers                      *
 *                                                        *
 * Allocates memory for all message passing buffers used  *
 * to send solution information between neighbouring      *
 * adaptive blocks.                                       *
 *                                                        *
 **********************************************************/
void Allocate_Message_Buffers(AdaptiveBlock2D_List &Blk_List,
                              const int Number_of_Solution_Variables) {

    Allocate_Message_Buffers_NoResChange(Blk_List,
                                         Number_of_Solution_Variables);
    Allocate_Message_Buffers_ResChange(Blk_List,
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
void Allocate_Message_Buffers_NoResChange(AdaptiveBlock2D_List &Blk_List,
                                          const int Number_of_Solution_Variables) {

    int i_blk, buffer_size, buffer_size_neighbour, l;

    /* Ensure that memory for the block index of the send and receive buffers
       has been allocated. */

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

    /* Reallocate memory for the send and receive buffers of each block. */
    
    for ( i_blk = 0; i_blk <= Blk_List.Nblk-1; ++i_blk ) {

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
                        Number_of_Solution_Variables+
                        2*Blk_List.Block[i_blk].info.dimen.ghost;
          buffer_size_neighbour = Blk_List.Block[i_blk].infoN[0].dimen.ghost*
                                  (abs(Blk_List.Block[i_blk].infoN[0].dimen.i)+1)*
                                  Number_of_Solution_Variables+
                                  2*Blk_List.Block[i_blk].infoN[0].dimen.ghost;
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
                        Number_of_Solution_Variables+
                        2*Blk_List.Block[i_blk].info.dimen.ghost;
          buffer_size_neighbour = Blk_List.Block[i_blk].infoS[0].dimen.ghost*
                                  (abs(Blk_List.Block[i_blk].infoS[0].dimen.i)+1)*
                                  Number_of_Solution_Variables+
                                  2*Blk_List.Block[i_blk].infoS[0].dimen.ghost;
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
                        Number_of_Solution_Variables+
                        2*Blk_List.Block[i_blk].info.dimen.ghost;
          buffer_size_neighbour = Blk_List.Block[i_blk].infoE[0].dimen.ghost*
                                  (abs(Blk_List.Block[i_blk].infoE[0].dimen.j)+1)*
                                  Number_of_Solution_Variables+
                                  2*Blk_List.Block[i_blk].infoE[0].dimen.ghost;
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
                        Number_of_Solution_Variables+
                        2*Blk_List.Block[i_blk].info.dimen.ghost;
          buffer_size_neighbour = Blk_List.Block[i_blk].infoW[0].dimen.ghost*
                                  (abs(Blk_List.Block[i_blk].infoW[0].dimen.j)+1)*
                                  Number_of_Solution_Variables+
                                  2*Blk_List.Block[i_blk].infoW[0].dimen.ghost;
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
                        Number_of_Solution_Variables;
          buffer_size_neighbour = sqr(Blk_List.Block[i_blk].infoNW[0].dimen.ghost)*
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
                        Number_of_Solution_Variables;
          buffer_size_neighbour = sqr(Blk_List.Block[i_blk].infoNE[0].dimen.ghost)*
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
                        Number_of_Solution_Variables;
          buffer_size_neighbour = sqr(Blk_List.Block[i_blk].infoSE[0].dimen.ghost)*
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
                        Number_of_Solution_Variables;
          buffer_size_neighbour = sqr(Blk_List.Block[i_blk].infoSW[0].dimen.ghost)*
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
void Allocate_Message_Buffers_ResChange(AdaptiveBlock2D_List &Blk_List,
                                        const int Number_of_Solution_Variables) {

    int i_blk, j_neigh, buffer_size, buffer_size_neighbour, l;

    /* Ensure that memory for the block index of the send and receive buffers
       has been allocated. */

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

    /* Reallocate memory for the send and receive buffers of each block. */
    
    for ( i_blk = 0; i_blk <= Blk_List.Nblk-1; ++i_blk ) {

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
void Deallocate_Message_Buffers(AdaptiveBlock2D_List &Blk_List) {

    Deallocate_Message_Buffers_NoResChange(Blk_List);
    Deallocate_Message_Buffers_ResChange(Blk_List);

}

/**********************************************************
 * Routine: Deallocate_Message_Buffers_NoResChange        *
 *                                                        *
 * Deallocates memory for all message passing buffers     *
 * used to send solution information between neighbouring *
 * adaptive blocks with no mesh resolution changes.       *
 *                                                        *
 **********************************************************/
void Deallocate_Message_Buffers_NoResChange(AdaptiveBlock2D_List &Blk_List) {

    int i_blk;

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
 
}

/**********************************************************
 * Routine: Deallocate_Message_Buffers_ResChange          *
 *                                                        *
 * Deallocates memory for all message passing buffers     *
 * used to send solution information between neighbouring *
 * adaptive blocks with mesh resolution changes.          *
 *                                                        *
 **********************************************************/
void Deallocate_Message_Buffers_ResChange(AdaptiveBlock2D_List &Blk_List) {

    int i_blk, j_neigh;

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

}

/**********************************************************
 * Routine: Exchange_Messages                             *
 *                                                        *
 * Sends and receives solution information contained in   *
 * the preloaded message passing buffers between          *
 * neighbouring adaptive blocks.                          *
 *                                                        *
 **********************************************************/
int Exchange_Messages(AdaptiveBlock2D_List &Blk_List,
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
void Copy_Refinement_Flags(AdaptiveBlock2D_List &Blk_List_1,
                           AdaptiveBlock2D_List &Blk_List_2) {
 
    int i_blk;

    if (Blk_List_1.Nblk > 0 &&
        Blk_List_1.Nblk == Blk_List_2.Nblk) {
       for ( i_blk = 0; i_blk <= Blk_List_1.Nblk-1 ; ++i_blk ) {
          Blk_List_1.RefineFlag[i_blk] = Blk_List_2.RefineFlag[i_blk];
       } /* endfor */
    } /* endif */

}
