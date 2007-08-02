/* AdaptiveBlock_NoMPI.cc:  Subroutines for adaptive blocks classes
                            on a single CPU architecture with no 
                            access to MPI. */

/* Include adaptive block header file. */

#ifndef _ADAPTIVBLOCK_INCLUDED
#include "AdaptiveBlock.h"
#endif // _ADAPTIVEBLOCK_INCLUDED

/*************************************************************
 * AdaptiveBlockResourceList -- External subroutines.        *
 *************************************************************/

/*************************************************************
 * AdaptiveBlock2D -- External subroutines.                  *
 *************************************************************/

/**********************************************************
 * Routine: Exchange_Messages_NoResChange                 *
 *                                                        *
 * Sends and receives solution information contained in   *
 * the preloaded message passing buffers between          *
 * neighbouring adaptive blocks with no mesh resolution   *
 * changes.                                               *
 *                                                        *
 **********************************************************/
int Exchange_Messages_NoResChange(AdaptiveBlock2D_List &Blk_List,
                                  const int Number_of_Solution_Variables) {

    int i_cpu, i_blk, neighbour_cpu, neighbour_blk, 
        buffer_size, buffer_size_neighbour, l;

    /* Perform message passing for the block faces and corners of each 
       solution block having no mesh resolution change. */

    i_cpu = Blk_List.ThisCPU;

    for ( i_blk = 0 ; i_blk <= Blk_List.Nblk-1 ; ++i_blk ) {

       //////////////////////////////////////////////////////////////
       // Exchange messages at the North faces of solution blocks. //
       //////////////////////////////////////////////////////////////
       if (Blk_List.Block[i_blk].used && 
           (Blk_List.Block[i_blk].nN == 1) &&
           (Blk_List.Block[i_blk].info.level == Blk_List.Block[i_blk].infoN[0].level)) {
	  neighbour_cpu = Blk_List.Block[i_blk].infoN[0].cpu;
	  neighbour_blk = Blk_List.Block[i_blk].infoN[0].blknum;
          if (neighbour_cpu != i_cpu) return(2200);
          if (!Blk_List.Block[neighbour_blk].used) return(2201); 
          buffer_size = Blk_List.Block[i_blk].infoN[0].dimen.ghost*
                        (abs(Blk_List.Block[i_blk].infoN[0].dimen.i)+1)*
                        Number_of_Solution_Variables+
                        2*Blk_List.Block[i_blk].infoN[0].dimen.ghost;
          buffer_size_neighbour = Blk_List.Block[neighbour_blk].info.dimen.ghost*
                                  (abs(Blk_List.Block[neighbour_blk].info.dimen.i)+1)*
                                  Number_of_Solution_Variables+
                                  2*Blk_List.Block[neighbour_blk].info.dimen.ghost;
          if (buffer_size != buffer_size_neighbour) return(2203);
          if (Blk_List.Block[i_blk].infoN[0].dimen.j > 0) {
             if (Blk_List.Block[neighbour_blk].nS != 1 ||
                 Blk_List.Block[neighbour_blk].infoS[0].cpu != i_cpu ||
                 Blk_List.Block[neighbour_blk].infoS[0].blknum != i_blk) return(2204);
             for (l = 0; l <= buffer_size-1; ++l) {
	        Blk_List.message_noreschange_southface_recbuf[neighbour_blk][l] =
                  Blk_List.message_noreschange_northface_sendbuf[i_blk][l];
             } /* endfor */
          } else {
             if (Blk_List.Block[neighbour_blk].nN != 1 ||
                 Blk_List.Block[neighbour_blk].infoN[0].cpu != i_cpu ||
                 Blk_List.Block[neighbour_blk].infoN[0].blknum != i_blk) return(2205);
             for (l = 0; l <= buffer_size-1; ++l) {
	        Blk_List.message_noreschange_northface_recbuf[neighbour_blk][l] =
                  Blk_List.message_noreschange_northface_sendbuf[i_blk][l];
             } /* endfor */
          } /* endif */
       } /* endif */

       //////////////////////////////////////////////////////////////
       // Exchange messages at the South faces of solution blocks. //
       //////////////////////////////////////////////////////////////
       if (Blk_List.Block[i_blk].used && 
           (Blk_List.Block[i_blk].nS == 1) &&
           (Blk_List.Block[i_blk].info.level == Blk_List.Block[i_blk].infoS[0].level)) {
	  neighbour_cpu = Blk_List.Block[i_blk].infoS[0].cpu;
	  neighbour_blk = Blk_List.Block[i_blk].infoS[0].blknum;
          if (neighbour_cpu != i_cpu) return(2100);
          if (!Blk_List.Block[neighbour_blk].used) return(2101); 
          buffer_size = Blk_List.Block[i_blk].infoS[0].dimen.ghost*
                        (abs(Blk_List.Block[i_blk].infoS[0].dimen.i)+1)*
                        Number_of_Solution_Variables+
                        2*Blk_List.Block[i_blk].infoS[0].dimen.ghost;
          buffer_size_neighbour = Blk_List.Block[neighbour_blk].info.dimen.ghost*
                                  (abs(Blk_List.Block[neighbour_blk].info.dimen.i)+1)*
                                  Number_of_Solution_Variables+
                                  2*Blk_List.Block[neighbour_blk].info.dimen.ghost;
          if (buffer_size != buffer_size_neighbour) return(2102);
          if (Blk_List.Block[i_blk].infoS[0].dimen.j > 0) {
             if (Blk_List.Block[neighbour_blk].nN != 1 ||
                 Blk_List.Block[neighbour_blk].infoN[0].cpu != i_cpu ||
                 Blk_List.Block[neighbour_blk].infoN[0].blknum != i_blk) return(2103);
             for (l = 0; l <= buffer_size-1; ++l) {
	        Blk_List.message_noreschange_northface_recbuf[neighbour_blk][l] =
                  Blk_List.message_noreschange_southface_sendbuf[i_blk][l];
             } /* endfor */
          } else {
             if (Blk_List.Block[neighbour_blk].nS != 1 ||
                 Blk_List.Block[neighbour_blk].infoS[0].cpu != i_cpu ||
                 Blk_List.Block[neighbour_blk].infoS[0].blknum != i_blk) return(2104);
             for (l = 0; l <= buffer_size-1; ++l) {
	        Blk_List.message_noreschange_southface_recbuf[neighbour_blk][l] =
                  Blk_List.message_noreschange_southface_sendbuf[i_blk][l];
             } /* endfor */
          } /* endif */
       } /* endif */

       /////////////////////////////////////////////////////////////
       // Exchange messages at the East faces of solution blocks. //
       /////////////////////////////////////////////////////////////
       if (Blk_List.Block[i_blk].used && 
           (Blk_List.Block[i_blk].nE == 1) &&
           (Blk_List.Block[i_blk].info.level == Blk_List.Block[i_blk].infoE[0].level)) {
	  neighbour_cpu = Blk_List.Block[i_blk].infoE[0].cpu;
	  neighbour_blk = Blk_List.Block[i_blk].infoE[0].blknum;
          if (neighbour_cpu != i_cpu) return(1200);
          if (!Blk_List.Block[neighbour_blk].used) return(1201); 
          buffer_size = Blk_List.Block[i_blk].infoE[0].dimen.ghost*
                        (abs(Blk_List.Block[i_blk].infoE[0].dimen.j)+1)*
                        Number_of_Solution_Variables+
                        2*Blk_List.Block[i_blk].infoE[0].dimen.ghost;
          buffer_size_neighbour = Blk_List.Block[neighbour_blk].info.dimen.ghost*
                                  (abs(Blk_List.Block[neighbour_blk].info.dimen.j)+1)*
                                  Number_of_Solution_Variables+
                                  2*Blk_List.Block[neighbour_blk].info.dimen.ghost;
          if (buffer_size != buffer_size_neighbour) return(1202);
          if (Blk_List.Block[i_blk].infoE[0].dimen.i > 0) {
             if (Blk_List.Block[neighbour_blk].nW != 1 ||
                 Blk_List.Block[neighbour_blk].infoW[0].cpu != i_cpu ||
                 Blk_List.Block[neighbour_blk].infoW[0].blknum != i_blk) return(1203);
             for (l = 0; l <= buffer_size-1; ++l) {
	        Blk_List.message_noreschange_westface_recbuf[neighbour_blk][l] =
                  Blk_List.message_noreschange_eastface_sendbuf[i_blk][l];
             } /* endfor */
          } else {
             if (Blk_List.Block[neighbour_blk].nE != 1 ||
                 Blk_List.Block[neighbour_blk].infoE[0].cpu != i_cpu ||
                 Blk_List.Block[neighbour_blk].infoE[0].blknum != i_blk) return(1204);
             for (l = 0; l <= buffer_size-1; ++l) {
	        Blk_List.message_noreschange_eastface_recbuf[neighbour_blk][l] =
                  Blk_List.message_noreschange_eastface_sendbuf[i_blk][l];
             } /* endfor */
          } /* endif */
       } /* endif */

       /////////////////////////////////////////////////////////////
       // Exchange messages at the West faces of solution blocks. //
       /////////////////////////////////////////////////////////////
       if (Blk_List.Block[i_blk].used && 
           (Blk_List.Block[i_blk].nW == 1) &&
           (Blk_List.Block[i_blk].info.level == Blk_List.Block[i_blk].infoW[0].level)) {
	  neighbour_cpu = Blk_List.Block[i_blk].infoW[0].cpu;
	  neighbour_blk = Blk_List.Block[i_blk].infoW[0].blknum;
          if (neighbour_cpu != i_cpu) return(1100);
          if (!Blk_List.Block[neighbour_blk].used) return(1101); 
          buffer_size = Blk_List.Block[i_blk].infoW[0].dimen.ghost*
                        (abs(Blk_List.Block[i_blk].infoW[0].dimen.j)+1)*
                        Number_of_Solution_Variables+
                        2*Blk_List.Block[i_blk].infoW[0].dimen.ghost;
          buffer_size_neighbour = Blk_List.Block[neighbour_blk].info.dimen.ghost*
                                  (abs(Blk_List.Block[neighbour_blk].info.dimen.j)+1)*
                                  Number_of_Solution_Variables+
                                  2*Blk_List.Block[neighbour_blk].info.dimen.ghost;
          if (buffer_size != buffer_size_neighbour) return(1102);
          if (Blk_List.Block[i_blk].infoW[0].dimen.i > 0) {
             if (Blk_List.Block[neighbour_blk].nE != 1 ||
                 Blk_List.Block[neighbour_blk].infoE[0].cpu != i_cpu ||
                 Blk_List.Block[neighbour_blk].infoE[0].blknum != i_blk) return(1103);
             for (l = 0; l <= buffer_size-1; ++l) {
	        Blk_List.message_noreschange_eastface_recbuf[neighbour_blk][l] =
                  Blk_List.message_noreschange_westface_sendbuf[i_blk][l];
             } /* endfor */
          } else {
             if (Blk_List.Block[neighbour_blk].nW != 1 ||
                 Blk_List.Block[neighbour_blk].infoW[0].cpu != i_cpu ||
                 Blk_List.Block[neighbour_blk].infoW[0].blknum != i_blk) return(1104);
             for (l = 0; l <= buffer_size-1; ++l) {
	        Blk_List.message_noreschange_westface_recbuf[neighbour_blk][l] =
                  Blk_List.message_noreschange_westface_sendbuf[i_blk][l];
             } /* endfor */
          } /* endif */
       } /* endif */

       /////////////////////////////////////////////////////////////////////
       // Exchange messages at the North West corners of solution blocks. //
       /////////////////////////////////////////////////////////////////////
       if (Blk_List.Block[i_blk].used && 
           (Blk_List.Block[i_blk].nNW == 1) &&
           (Blk_List.Block[i_blk].info.level == Blk_List.Block[i_blk].infoNW[0].level)) {
	  neighbour_cpu = Blk_List.Block[i_blk].infoNW[0].cpu;
	  neighbour_blk = Blk_List.Block[i_blk].infoNW[0].blknum;
          if (neighbour_cpu != i_cpu) return(4200);
          if (!Blk_List.Block[neighbour_blk].used) return(4201); 
          buffer_size = sqr(Blk_List.Block[i_blk].infoNW[0].dimen.ghost)*
                        Number_of_Solution_Variables;
          buffer_size_neighbour = sqr(Blk_List.Block[neighbour_blk].info.dimen.ghost)*
                                  Number_of_Solution_Variables;
          if (buffer_size != buffer_size_neighbour) return(4202);
          if (Blk_List.Block[i_blk].infoNW[0].dimen.i > 0 &&
              Blk_List.Block[i_blk].infoNW[0].dimen.j > 0) {
             if (Blk_List.Block[neighbour_blk].nSE != 1 ||
                 Blk_List.Block[neighbour_blk].infoSE[0].cpu != i_cpu ||
                 Blk_List.Block[neighbour_blk].infoSE[0].blknum != i_blk) return(4203);
             for (l = 0; l <= buffer_size-1; ++l) {
	        Blk_List.message_noreschange_southeastcorner_recbuf[neighbour_blk][l] =
                  Blk_List.message_noreschange_northwestcorner_sendbuf[i_blk][l];
             } /* endfor */
          } else if (Blk_List.Block[i_blk].infoNW[0].dimen.j > 0) {
             if (Blk_List.Block[neighbour_blk].nSW != 1 ||
                 Blk_List.Block[neighbour_blk].infoSW[0].cpu != i_cpu ||
                 Blk_List.Block[neighbour_blk].infoSW[0].blknum != i_blk) return(4204);
             for (l = 0; l <= buffer_size-1; ++l) {
	        Blk_List.message_noreschange_southwestcorner_recbuf[neighbour_blk][l] =
                  Blk_List.message_noreschange_northwestcorner_sendbuf[i_blk][l];
             } /* endfor */
          } else if (Blk_List.Block[i_blk].infoNW[0].dimen.i > 0) {
             if (Blk_List.Block[neighbour_blk].nNE != 1 ||
                 Blk_List.Block[neighbour_blk].infoNE[0].cpu != i_cpu ||
                 Blk_List.Block[neighbour_blk].infoNE[0].blknum != i_blk) return(4205);
             for (l = 0; l <= buffer_size-1; ++l) {
	        Blk_List.message_noreschange_northeastcorner_recbuf[neighbour_blk][l] =
                  Blk_List.message_noreschange_northwestcorner_sendbuf[i_blk][l];
             } /* endfor */
          } else {
             if (Blk_List.Block[neighbour_blk].nNW != 1 ||
                 Blk_List.Block[neighbour_blk].infoNW[0].cpu != i_cpu ||
                 Blk_List.Block[neighbour_blk].infoNW[0].blknum != i_blk) return(4206);
             for (l = 0; l <= buffer_size-1; ++l) {
	        Blk_List.message_noreschange_northwestcorner_recbuf[neighbour_blk][l] =
                  Blk_List.message_noreschange_northwestcorner_sendbuf[i_blk][l];
             } /* endfor */
          } /* endif */
       } /* endif */

       /////////////////////////////////////////////////////////////////////
       // Exchange messages at the North East corners of solution blocks. //
       /////////////////////////////////////////////////////////////////////
       if (Blk_List.Block[i_blk].used && 
           (Blk_List.Block[i_blk].nNE == 1) &&
           (Blk_List.Block[i_blk].info.level == Blk_List.Block[i_blk].infoNE[0].level)) {
	  neighbour_cpu = Blk_List.Block[i_blk].infoNE[0].cpu;
	  neighbour_blk = Blk_List.Block[i_blk].infoNE[0].blknum;
          if (neighbour_cpu != i_cpu) return(5200);
          if (!Blk_List.Block[neighbour_blk].used) return(5201); 
          buffer_size = sqr(Blk_List.Block[i_blk].infoNE[0].dimen.ghost)*
                        Number_of_Solution_Variables;
          buffer_size_neighbour = sqr(Blk_List.Block[neighbour_blk].info.dimen.ghost)*
                                  Number_of_Solution_Variables;
          if (buffer_size != buffer_size_neighbour) return(5202);
          if (Blk_List.Block[i_blk].infoNE[0].dimen.i > 0 &&
              Blk_List.Block[i_blk].infoNE[0].dimen.j > 0) {
             if (Blk_List.Block[neighbour_blk].nSW != 1 ||
                 Blk_List.Block[neighbour_blk].infoSW[0].cpu != i_cpu ||
                 Blk_List.Block[neighbour_blk].infoSW[0].blknum != i_blk) return(5203);
             for (l = 0; l <= buffer_size-1; ++l) {
	        Blk_List.message_noreschange_southwestcorner_recbuf[neighbour_blk][l] =
                  Blk_List.message_noreschange_northeastcorner_sendbuf[i_blk][l];
             } /* endfor */
          } else if (Blk_List.Block[i_blk].infoNE[0].dimen.j > 0) {
             if (Blk_List.Block[neighbour_blk].nSE != 1 ||
                 Blk_List.Block[neighbour_blk].infoSE[0].cpu != i_cpu ||
                 Blk_List.Block[neighbour_blk].infoSE[0].blknum != i_blk) return(5204);
             for (l = 0; l <= buffer_size-1; ++l) {
	        Blk_List.message_noreschange_southeastcorner_recbuf[neighbour_blk][l] =
                  Blk_List.message_noreschange_northeastcorner_sendbuf[i_blk][l];
             } /* endfor */
          } else if (Blk_List.Block[i_blk].infoNE[0].dimen.i > 0) {
             if (Blk_List.Block[neighbour_blk].nNW != 1 ||
                 Blk_List.Block[neighbour_blk].infoNW[0].cpu != i_cpu ||
                 Blk_List.Block[neighbour_blk].infoNW[0].blknum != i_blk) return(5205);
             for (l = 0; l <= buffer_size-1; ++l) {
	        Blk_List.message_noreschange_northwestcorner_recbuf[neighbour_blk][l] =
                  Blk_List.message_noreschange_northeastcorner_sendbuf[i_blk][l];
             } /* endfor */
          } else {
             if (Blk_List.Block[neighbour_blk].nNE != 1 ||
                 Blk_List.Block[neighbour_blk].infoNE[0].cpu != i_cpu ||
                 Blk_List.Block[neighbour_blk].infoNE[0].blknum != i_blk) return(5206);
             for (l = 0; l <= buffer_size-1; ++l) {
	        Blk_List.message_noreschange_northeastcorner_recbuf[neighbour_blk][l] =
                  Blk_List.message_noreschange_northeastcorner_sendbuf[i_blk][l];
             } /* endfor */
          } /* endif */
       } /* endif */

       /////////////////////////////////////////////////////////////////////
       // Exchange messages at the South East corners of solution blocks. //
       /////////////////////////////////////////////////////////////////////
       if (Blk_List.Block[i_blk].used && 
           (Blk_List.Block[i_blk].nSE == 1) &&
           (Blk_List.Block[i_blk].info.level == Blk_List.Block[i_blk].infoSE[0].level)) {
	  neighbour_cpu = Blk_List.Block[i_blk].infoSE[0].cpu;
	  neighbour_blk = Blk_List.Block[i_blk].infoSE[0].blknum;
          if (neighbour_cpu != i_cpu) return(4100);
          if (!Blk_List.Block[neighbour_blk].used) return(4101); 
          buffer_size = sqr(Blk_List.Block[i_blk].infoSE[0].dimen.ghost)*
                        Number_of_Solution_Variables;
          buffer_size_neighbour = sqr(Blk_List.Block[neighbour_blk].info.dimen.ghost)*
                                  Number_of_Solution_Variables;
          if (buffer_size != buffer_size_neighbour) return(4102);
          if (Blk_List.Block[i_blk].infoSE[0].dimen.i > 0 &&
              Blk_List.Block[i_blk].infoSE[0].dimen.j > 0) {
             if (Blk_List.Block[neighbour_blk].nNW != 1 ||
                 Blk_List.Block[neighbour_blk].infoNW[0].cpu != i_cpu ||
                 Blk_List.Block[neighbour_blk].infoNW[0].blknum != i_blk) return(4103);
             for (l = 0; l <= buffer_size-1; ++l) {
	        Blk_List.message_noreschange_northwestcorner_recbuf[neighbour_blk][l] =
                  Blk_List.message_noreschange_southeastcorner_sendbuf[i_blk][l];
             } /* endfor */
          } else if (Blk_List.Block[i_blk].infoSE[0].dimen.j > 0) {
             if (Blk_List.Block[neighbour_blk].nNE != 1 ||
                 Blk_List.Block[neighbour_blk].infoNE[0].cpu != i_cpu ||
                 Blk_List.Block[neighbour_blk].infoNE[0].blknum != i_blk) return(4104);
             for (l = 0; l <= buffer_size-1; ++l) {
	        Blk_List.message_noreschange_northeastcorner_recbuf[neighbour_blk][l] =
                  Blk_List.message_noreschange_southeastcorner_sendbuf[i_blk][l];
             } /* endfor */
          } else if (Blk_List.Block[i_blk].infoSE[0].dimen.i > 0) {
             if (Blk_List.Block[neighbour_blk].nSW != 1 ||
                 Blk_List.Block[neighbour_blk].infoSW[0].cpu != i_cpu ||
                 Blk_List.Block[neighbour_blk].infoSW[0].blknum != i_blk) return(4105);
             for (l = 0; l <= buffer_size-1; ++l) {
	        Blk_List.message_noreschange_southwestcorner_recbuf[neighbour_blk][l] =
                  Blk_List.message_noreschange_southeastcorner_sendbuf[i_blk][l];
             } /* endfor */
          } else {
             if (Blk_List.Block[neighbour_blk].nSE != 1 ||
                 Blk_List.Block[neighbour_blk].infoSE[0].cpu != i_cpu ||
                 Blk_List.Block[neighbour_blk].infoSE[0].blknum != i_blk) return(4106);
             for (l = 0; l <= buffer_size-1; ++l) {
	        Blk_List.message_noreschange_southeastcorner_recbuf[neighbour_blk][l] =
                  Blk_List.message_noreschange_southeastcorner_sendbuf[i_blk][l];
             } /* endfor */
          } /* endif */
       } /* endif */

       /////////////////////////////////////////////////////////////////////
       // Exchange messages at the South West corners of solution blocks. //
       /////////////////////////////////////////////////////////////////////
       if (Blk_List.Block[i_blk].used && 
           (Blk_List.Block[i_blk].nSW == 1) &&
           (Blk_List.Block[i_blk].info.level == Blk_List.Block[i_blk].infoSW[0].level)) {
	  neighbour_cpu = Blk_List.Block[i_blk].infoSW[0].cpu;
	  neighbour_blk = Blk_List.Block[i_blk].infoSW[0].blknum;
          if (neighbour_cpu != i_cpu) return(5100);
          if (!Blk_List.Block[neighbour_blk].used) return(5101); 
          buffer_size = sqr(Blk_List.Block[i_blk].infoSW[0].dimen.ghost)*
                        Number_of_Solution_Variables;
          buffer_size_neighbour = sqr(Blk_List.Block[neighbour_blk].info.dimen.ghost)*
                                  Number_of_Solution_Variables;
          if (buffer_size != buffer_size_neighbour) return(5102);
          if (Blk_List.Block[i_blk].infoSW[0].dimen.i > 0 &&
              Blk_List.Block[i_blk].infoSW[0].dimen.j > 0) {
             if (Blk_List.Block[neighbour_blk].nNE != 1 ||
                 Blk_List.Block[neighbour_blk].infoNE[0].cpu != i_cpu ||
                 Blk_List.Block[neighbour_blk].infoNE[0].blknum != i_blk) return(5103);
             for (l = 0; l <= buffer_size-1; ++l) {
	        Blk_List.message_noreschange_northeastcorner_recbuf[neighbour_blk][l] =
                  Blk_List.message_noreschange_southwestcorner_sendbuf[i_blk][l];
             } /* endfor */
          } else if (Blk_List.Block[i_blk].infoSW[0].dimen.j > 0) {
             if (Blk_List.Block[neighbour_blk].nNW != 1 ||
                 Blk_List.Block[neighbour_blk].infoNW[0].cpu != i_cpu ||
                 Blk_List.Block[neighbour_blk].infoNW[0].blknum != i_blk) return(5104);
             for (l = 0; l <= buffer_size-1; ++l) {
	        Blk_List.message_noreschange_northwestcorner_recbuf[neighbour_blk][l] =
                  Blk_List.message_noreschange_southwestcorner_sendbuf[i_blk][l];
             } /* endfor */
          } else if (Blk_List.Block[i_blk].infoSW[0].dimen.i > 0) {
             if (Blk_List.Block[neighbour_blk].nSE != 1 ||
                 Blk_List.Block[neighbour_blk].infoSE[0].cpu != i_cpu ||
                 Blk_List.Block[neighbour_blk].infoSE[0].blknum != i_blk) return(5105);
             for (l = 0; l <= buffer_size-1; ++l) {
	        Blk_List.message_noreschange_southeastcorner_recbuf[neighbour_blk][l] =
                  Blk_List.message_noreschange_southwestcorner_sendbuf[i_blk][l];
             } /* endfor */
          } else {
             if (Blk_List.Block[neighbour_blk].nSW != 1 ||
                 Blk_List.Block[neighbour_blk].infoSW[0].cpu != i_cpu ||
                 Blk_List.Block[neighbour_blk].infoSW[0].blknum != i_blk) return(5106);
             for (l = 0; l <= buffer_size-1; ++l) {
	        Blk_List.message_noreschange_southwestcorner_recbuf[neighbour_blk][l] =
                  Blk_List.message_noreschange_southwestcorner_sendbuf[i_blk][l];
             } /* endfor */
          } /* endif */
       } /* endif */

    }  /* endfor */

    /* Message passing complete.  Return zero value. */

    return(0);

}

/**********************************************************
 * Routine: Exchange_Messages_ResChange_FineToCoarse      *
 *                                                        *
 * Sends and receives solution information contained in   *
 * the preloaded message passing buffers between          *
 * neighbouring adaptive blocks, sending solution         *
 * information from more refined (higher mesh resolution) *
 * adaptive blocks to more coarse adaptive solution       *
 * blocks with lower mesh resolution.                     *
 *                                                        *
 **********************************************************/
int Exchange_Messages_ResChange_FineToCoarse(AdaptiveBlock2D_List &Blk_List,
                                             const int Number_of_Solution_Variables) {

    int i_cpu, i_blk, neighbour_cpu, neighbour_blk, 
        buffer_size, buffer_size_neighbour, l, j_neigh;

    /* Perform message passing for the block faces and corners of each 
       more refined (higher mesh resolution) solution block having more 
       coarse (lower mesh resolution) neighbouring solution blocks. */

    i_cpu = Blk_List.ThisCPU;

    for ( i_blk = 0 ; i_blk <= Blk_List.Nblk-1 ; ++i_blk ) {

       //////////////////////////////////////////////////////////////
       // Exchange messages at the North faces of solution blocks  //
       // with coarser neighbour at lower mesh resolution.         //
       //////////////////////////////////////////////////////////////
       if (Blk_List.Block[i_blk].used && 
           (Blk_List.Block[i_blk].nN == 1) &&
           (Blk_List.Block[i_blk].info.level >
            Blk_List.Block[i_blk].infoN[0].level)) {
	  neighbour_cpu = Blk_List.Block[i_blk].infoN[0].cpu;
	  neighbour_blk = Blk_List.Block[i_blk].infoN[0].blknum;
          if (neighbour_cpu != i_cpu) return(2200);
          if (!Blk_List.Block[neighbour_blk].used) return(2201); 
          buffer_size = Blk_List.Block[i_blk].infoN[0].dimen.ghost*
                        (abs(Blk_List.Block[i_blk].infoN[0].dimen.i)/2+1)*
                        Number_of_Solution_Variables+
                        2*Blk_List.Block[i_blk].infoN[0].dimen.ghost;
          buffer_size_neighbour = Blk_List.Block[neighbour_blk].info.dimen.ghost*
                                  (abs(Blk_List.Block[neighbour_blk].info.dimen.i)/2+1)*
                                  Number_of_Solution_Variables+
                                  2*Blk_List.Block[neighbour_blk].info.dimen.ghost;
          if (buffer_size != buffer_size_neighbour) return(2202);
          if (Blk_List.Block[i_blk].infoN[0].dimen.j > 0) {
	     if (Blk_List.Block[neighbour_blk].nS == 2) {
                if (Blk_List.Block[neighbour_blk].infoS[0].cpu == i_cpu &&
                    Blk_List.Block[neighbour_blk].infoS[0].blknum == i_blk) {
                   j_neigh = 0;
                } else if (Blk_List.Block[neighbour_blk].infoS[1].cpu == i_cpu &&
                           Blk_List.Block[neighbour_blk].infoS[1].blknum == i_blk) {
                   j_neigh = 1;
                } else {
                   return(2203);
                } /* endif */
             } else {
                return(2204);
             } /* endif */
             for (l = 0; l <= buffer_size-1; ++l) {
	        Blk_List.message_reschange_southface_recbuf[neighbour_blk][j_neigh][l] =
                  Blk_List.message_reschange_northface_sendbuf[i_blk][0][l];
             } /* endfor */
          } else {
	     if (Blk_List.Block[neighbour_blk].nN == 2) {
                if (Blk_List.Block[neighbour_blk].infoN[0].cpu == i_cpu &&
                    Blk_List.Block[neighbour_blk].infoN[0].blknum == i_blk) {
                   j_neigh = 0;
                } else if (Blk_List.Block[neighbour_blk].infoN[1].cpu == i_cpu &&
                           Blk_List.Block[neighbour_blk].infoN[1].blknum == i_blk) {
                   j_neigh = 1;
                } else {
                   return(2205);
                } /* endif */
             } else {
                return(2206);
             } /* endif */
             for (l = 0; l <= buffer_size-1; ++l) {
	        Blk_List.message_reschange_northface_recbuf[neighbour_blk][j_neigh][l] =
                  Blk_List.message_reschange_northface_sendbuf[i_blk][0][l];
             } /* endfor */
          } /* endif */
       } /* endif */

       //////////////////////////////////////////////////////////////
       // Exchange messages at the South faces of solution blocks  //
       // with coarser neighbour at lower mesh resolution.         //
       //////////////////////////////////////////////////////////////
       if (Blk_List.Block[i_blk].used && 
           (Blk_List.Block[i_blk].nS == 1) &&
           (Blk_List.Block[i_blk].info.level >
            Blk_List.Block[i_blk].infoS[0].level)) {
	  neighbour_cpu = Blk_List.Block[i_blk].infoS[0].cpu;
	  neighbour_blk = Blk_List.Block[i_blk].infoS[0].blknum;
          if (neighbour_cpu != i_cpu) return(2100);
          if (!Blk_List.Block[neighbour_blk].used) return(2101); 
          buffer_size = Blk_List.Block[i_blk].infoS[0].dimen.ghost*
                        (abs(Blk_List.Block[i_blk].infoS[0].dimen.i)/2+1)*
                        Number_of_Solution_Variables+
                        2*Blk_List.Block[i_blk].infoS[0].dimen.ghost;
          buffer_size_neighbour = Blk_List.Block[neighbour_blk].info.dimen.ghost*
                                  (abs(Blk_List.Block[neighbour_blk].info.dimen.i)/2+1)*
                                  Number_of_Solution_Variables+
                                  2*Blk_List.Block[neighbour_blk].info.dimen.ghost;
          if (buffer_size != buffer_size_neighbour) return(2102);
          if (Blk_List.Block[i_blk].infoS[0].dimen.j > 0) {
	     if (Blk_List.Block[neighbour_blk].nN == 2) {
                if (Blk_List.Block[neighbour_blk].infoN[0].cpu == i_cpu &&
                    Blk_List.Block[neighbour_blk].infoN[0].blknum == i_blk) {
                   j_neigh = 0;
                } else if (Blk_List.Block[neighbour_blk].infoN[1].cpu == i_cpu &&
                           Blk_List.Block[neighbour_blk].infoN[1].blknum == i_blk) {
                   j_neigh = 1;
                } else {
                   return(2103);
                } /* endif */
             } else {
                return(2104);
             } /* endif */
             for (l = 0; l <= buffer_size-1; ++l) {
	        Blk_List.message_reschange_northface_recbuf[neighbour_blk][j_neigh][l] =
                  Blk_List.message_reschange_southface_sendbuf[i_blk][0][l];
             } /* endfor */
          } else {
	     if (Blk_List.Block[neighbour_blk].nS == 2) {
                if (Blk_List.Block[neighbour_blk].infoS[0].cpu == i_cpu &&
                    Blk_List.Block[neighbour_blk].infoS[0].blknum == i_blk) {
                   j_neigh = 0;
                } else if (Blk_List.Block[neighbour_blk].infoS[1].cpu == i_cpu &&
                           Blk_List.Block[neighbour_blk].infoS[1].blknum == i_blk) {
                   j_neigh = 1;
                } else {
                   return(2105);
                } /* endif */
             } else {
                return(2106);
             } /* endif */
             for (l = 0; l <= buffer_size-1; ++l) {
	        Blk_List.message_reschange_southface_recbuf[neighbour_blk][j_neigh][l] =
                  Blk_List.message_reschange_southface_sendbuf[i_blk][0][l];
             } /* endfor */
          } /* endif */
       } /* endif */

       //////////////////////////////////////////////////////////////
       // Exchange messages at the East faces of solution blocks   //
       // with coarser neighbour at lower mesh resolution.         //
       //////////////////////////////////////////////////////////////
       if (Blk_List.Block[i_blk].used && 
           (Blk_List.Block[i_blk].nE == 1) &&
           (Blk_List.Block[i_blk].info.level >
            Blk_List.Block[i_blk].infoE[0].level)) {
	  neighbour_cpu = Blk_List.Block[i_blk].infoE[0].cpu;
	  neighbour_blk = Blk_List.Block[i_blk].infoE[0].blknum;
          if (neighbour_cpu != i_cpu) return(1200);
          if (!Blk_List.Block[neighbour_blk].used) return(1201); 
          buffer_size = Blk_List.Block[i_blk].infoE[0].dimen.ghost*
                        (abs(Blk_List.Block[i_blk].infoE[0].dimen.j)/2+1)*
                        Number_of_Solution_Variables+
                        2*Blk_List.Block[i_blk].infoE[0].dimen.ghost;
          buffer_size_neighbour = Blk_List.Block[neighbour_blk].info.dimen.ghost*
                                  (abs(Blk_List.Block[neighbour_blk].info.dimen.j)/2+1)*
                                  Number_of_Solution_Variables+
                                  2*Blk_List.Block[neighbour_blk].info.dimen.ghost;
          if (buffer_size != buffer_size_neighbour) return(1202);
          if (Blk_List.Block[i_blk].infoE[0].dimen.i > 0) {
	     if (Blk_List.Block[neighbour_blk].nW == 2) {
                if (Blk_List.Block[neighbour_blk].infoW[0].cpu == i_cpu &&
                    Blk_List.Block[neighbour_blk].infoW[0].blknum == i_blk) {
                   j_neigh = 0;
                } else if (Blk_List.Block[neighbour_blk].infoW[1].cpu == i_cpu &&
                           Blk_List.Block[neighbour_blk].infoW[1].blknum == i_blk) {
                   j_neigh = 1;
                } else {
                   return(1203);
                } /* endif */
             } else {
                return(1204);
             } /* endif */
             for (l = 0; l <= buffer_size-1; ++l) {
	        Blk_List.message_reschange_westface_recbuf[neighbour_blk][j_neigh][l] =
                  Blk_List.message_reschange_eastface_sendbuf[i_blk][0][l];
             } /* endfor */
          } else {
	     if (Blk_List.Block[neighbour_blk].nE == 2) {
                if (Blk_List.Block[neighbour_blk].infoE[0].cpu == i_cpu &&
                    Blk_List.Block[neighbour_blk].infoE[0].blknum == i_blk) {
                   j_neigh = 0;
                } else if (Blk_List.Block[neighbour_blk].infoE[1].cpu == i_cpu &&
                           Blk_List.Block[neighbour_blk].infoE[1].blknum == i_blk) {
                   j_neigh = 1;
                } else {
                   return(1205);
                } /* endif */
             } else {
                return(1206);
             } /* endif */
             for (l = 0; l <= buffer_size-1; ++l) {
	        Blk_List.message_reschange_eastface_recbuf[neighbour_blk][j_neigh][l] =
                  Blk_List.message_reschange_eastface_sendbuf[i_blk][0][l];
             } /* endfor */
          } /* endif */
       } /* endif */

       //////////////////////////////////////////////////////////////
       // Exchange messages at the West faces of solution blocks   //
       // with coarser neighbour at lower mesh resolution.         //
       //////////////////////////////////////////////////////////////
       if (Blk_List.Block[i_blk].used && 
           (Blk_List.Block[i_blk].nW == 1) &&
           (Blk_List.Block[i_blk].info.level >
            Blk_List.Block[i_blk].infoW[0].level)) {
	  neighbour_cpu = Blk_List.Block[i_blk].infoW[0].cpu;
	  neighbour_blk = Blk_List.Block[i_blk].infoW[0].blknum;
          if (neighbour_cpu != i_cpu) return(1100);
          if (!Blk_List.Block[neighbour_blk].used) return(1101); 
          buffer_size = Blk_List.Block[i_blk].infoW[0].dimen.ghost*
                        (abs(Blk_List.Block[i_blk].infoW[0].dimen.j)/2+1)*
                        Number_of_Solution_Variables+
                        2*Blk_List.Block[i_blk].infoW[0].dimen.ghost;
          buffer_size_neighbour = Blk_List.Block[neighbour_blk].info.dimen.ghost*
                                  (abs(Blk_List.Block[neighbour_blk].info.dimen.j)/2+1)*
                                  Number_of_Solution_Variables+
                                  2*Blk_List.Block[neighbour_blk].info.dimen.ghost;
          if (buffer_size != buffer_size_neighbour) return(1102);
          if (Blk_List.Block[i_blk].infoW[0].dimen.i > 0) {
	     if (Blk_List.Block[neighbour_blk].nE == 2) {
                if (Blk_List.Block[neighbour_blk].infoE[0].cpu == i_cpu &&
                    Blk_List.Block[neighbour_blk].infoE[0].blknum == i_blk) {
                   j_neigh = 0;
                } else if (Blk_List.Block[neighbour_blk].infoE[1].cpu == i_cpu &&
                           Blk_List.Block[neighbour_blk].infoE[1].blknum == i_blk) {
                   j_neigh = 1;
                } else {
                   return(1103);
                } /* endif */
             } else {
                return(1104);
             } /* endif */
             for (l = 0; l <= buffer_size-1; ++l) {
	        Blk_List.message_reschange_eastface_recbuf[neighbour_blk][j_neigh][l] =
                  Blk_List.message_reschange_westface_sendbuf[i_blk][0][l];
             } /* endfor */
          } else {
	     if (Blk_List.Block[neighbour_blk].nW == 2) {
                if (Blk_List.Block[neighbour_blk].infoW[0].cpu == i_cpu &&
                    Blk_List.Block[neighbour_blk].infoW[0].blknum == i_blk) {
                   j_neigh = 0;
                } else if (Blk_List.Block[neighbour_blk].infoW[1].cpu == i_cpu &&
                           Blk_List.Block[neighbour_blk].infoW[1].blknum == i_blk) {
                   j_neigh = 1;
                } else {
                   return(1105);
                } /* endif */
             } else {
                return(1106);
             } /* endif */
             for (l = 0; l <= buffer_size-1; ++l) {
	        Blk_List.message_reschange_westface_recbuf[neighbour_blk][j_neigh][l] =
                  Blk_List.message_reschange_westface_sendbuf[i_blk][0][l];
             } /* endfor */
          } /* endif */
       } /* endif */

       //////////////////////////////////////////////////////////////
       // Exchange messages at the North West corners of solution  //
       // blocks with coarser neighbour at lower mesh resolution.  //
       //////////////////////////////////////////////////////////////
       if (Blk_List.Block[i_blk].used && 
           (Blk_List.Block[i_blk].nNW == 1) &&
           (Blk_List.Block[i_blk].info.level >
            Blk_List.Block[i_blk].infoNW[0].level)) {
	  neighbour_cpu = Blk_List.Block[i_blk].infoNW[0].cpu;
	  neighbour_blk = Blk_List.Block[i_blk].infoNW[0].blknum;
          if (neighbour_cpu != i_cpu) return(4200);
          if (!Blk_List.Block[neighbour_blk].used) return(4201); 
          buffer_size = sqr(Blk_List.Block[i_blk].infoNW[0].dimen.ghost)*
                        Number_of_Solution_Variables;
          buffer_size_neighbour = sqr(Blk_List.Block[neighbour_blk].info.dimen.ghost)*
                                  Number_of_Solution_Variables;
          if (buffer_size != buffer_size_neighbour) return(4202);
          if (Blk_List.Block[i_blk].infoNW[0].dimen.i > 0 &&
              Blk_List.Block[i_blk].infoNW[0].dimen.j > 0) {
             if (Blk_List.Block[neighbour_blk].nSE != 1 ||
                 Blk_List.Block[neighbour_blk].infoSE[0].cpu != i_cpu ||
                 Blk_List.Block[neighbour_blk].infoSE[0].blknum != i_blk) return(4203);
             for (l = 0; l <= buffer_size-1; ++l) {
	        Blk_List.message_reschange_southeastcorner_recbuf[neighbour_blk][l] =
                  Blk_List.message_reschange_northwestcorner_sendbuf[i_blk][l];
             } /* endfor */
          } else if (Blk_List.Block[i_blk].infoNW[0].dimen.j > 0) {
             if (Blk_List.Block[neighbour_blk].nSW != 1 ||
                 Blk_List.Block[neighbour_blk].infoSW[0].cpu != i_cpu ||
                 Blk_List.Block[neighbour_blk].infoSW[0].blknum != i_blk) return(4204);
             for (l = 0; l <= buffer_size-1; ++l) {
	        Blk_List.message_reschange_southwestcorner_recbuf[neighbour_blk][l] =
                  Blk_List.message_reschange_northwestcorner_sendbuf[i_blk][l];
             } /* endfor */
          } else if (Blk_List.Block[i_blk].infoNW[0].dimen.i > 0) {
             if (Blk_List.Block[neighbour_blk].nNE != 1 ||
                 Blk_List.Block[neighbour_blk].infoNE[0].cpu != i_cpu ||
                 Blk_List.Block[neighbour_blk].infoNE[0].blknum != i_blk) return(4205);
             for (l = 0; l <= buffer_size-1; ++l) {
	        Blk_List.message_reschange_northeastcorner_recbuf[neighbour_blk][l] =
                  Blk_List.message_reschange_northwestcorner_sendbuf[i_blk][l];
             } /* endfor */
          } else {
             if (Blk_List.Block[neighbour_blk].nNW != 1 ||
                 Blk_List.Block[neighbour_blk].infoNW[0].cpu != i_cpu ||
                 Blk_List.Block[neighbour_blk].infoNW[0].blknum != i_blk) return(4206);
             for (l = 0; l <= buffer_size-1; ++l) {
	        Blk_List.message_reschange_northwestcorner_recbuf[neighbour_blk][l] =
                  Blk_List.message_reschange_northwestcorner_sendbuf[i_blk][l];
             } /* endfor */
          } /* endif */
       } /* endif */

       //////////////////////////////////////////////////////////////
       // Exchange messages at the North East corners of solution  //
       // blocks with coarser neighbour at lower mesh resolution.  //
       //////////////////////////////////////////////////////////////
       if (Blk_List.Block[i_blk].used && 
           (Blk_List.Block[i_blk].nNE == 1) &&
           (Blk_List.Block[i_blk].info.level >
            Blk_List.Block[i_blk].infoNE[0].level)) {
	  neighbour_cpu = Blk_List.Block[i_blk].infoNE[0].cpu;
	  neighbour_blk = Blk_List.Block[i_blk].infoNE[0].blknum;
          if (neighbour_cpu != i_cpu) return(5200);
          if (!Blk_List.Block[neighbour_blk].used) return(5201); 
          buffer_size = sqr(Blk_List.Block[i_blk].infoNE[0].dimen.ghost)*
                        Number_of_Solution_Variables;
          buffer_size_neighbour = sqr(Blk_List.Block[neighbour_blk].info.dimen.ghost)*
                                  Number_of_Solution_Variables;
          if (buffer_size != buffer_size_neighbour) return(5202);
          if (Blk_List.Block[i_blk].infoNE[0].dimen.i > 0 &&
              Blk_List.Block[i_blk].infoNE[0].dimen.j > 0) {
             if (Blk_List.Block[neighbour_blk].nSW != 1 ||
                 Blk_List.Block[neighbour_blk].infoSW[0].cpu != i_cpu ||
                 Blk_List.Block[neighbour_blk].infoSW[0].blknum != i_blk) return(5203);
             for (l = 0; l <= buffer_size-1; ++l) {
	        Blk_List.message_reschange_southwestcorner_recbuf[neighbour_blk][l] =
                  Blk_List.message_reschange_northeastcorner_sendbuf[i_blk][l];
             } /* endfor */
          } else if (Blk_List.Block[i_blk].infoNE[0].dimen.j > 0) {
             if (Blk_List.Block[neighbour_blk].nSE != 1 ||
                 Blk_List.Block[neighbour_blk].infoSE[0].cpu != i_cpu ||
                 Blk_List.Block[neighbour_blk].infoSE[0].blknum != i_blk) return(5204);
             for (l = 0; l <= buffer_size-1; ++l) {
	        Blk_List.message_reschange_southeastcorner_recbuf[neighbour_blk][l] =
                  Blk_List.message_reschange_northeastcorner_sendbuf[i_blk][l];
             } /* endfor */
          } else if (Blk_List.Block[i_blk].infoNE[0].dimen.i > 0) {
             if (Blk_List.Block[neighbour_blk].nNW != 1 ||
                 Blk_List.Block[neighbour_blk].infoNW[0].cpu != i_cpu ||
                 Blk_List.Block[neighbour_blk].infoNW[0].blknum != i_blk) return(5205);
             for (l = 0; l <= buffer_size-1; ++l) {
	        Blk_List.message_reschange_northwestcorner_recbuf[neighbour_blk][l] =
                  Blk_List.message_reschange_northeastcorner_sendbuf[i_blk][l];
             } /* endfor */
          } else {
             if (Blk_List.Block[neighbour_blk].nNE != 1 ||
                 Blk_List.Block[neighbour_blk].infoNE[0].cpu != i_cpu ||
                 Blk_List.Block[neighbour_blk].infoNE[0].blknum != i_blk) return(5206);
             for (l = 0; l <= buffer_size-1; ++l) {
	        Blk_List.message_reschange_northeastcorner_recbuf[neighbour_blk][l] =
                  Blk_List.message_reschange_northeastcorner_sendbuf[i_blk][l];
             } /* endfor */
          } /* endif */
       } /* endif */

       //////////////////////////////////////////////////////////////
       // Exchange messages at the South East corners of solution  //
       // blocks with coarser neighbour at lower mesh resolution.  //
       //////////////////////////////////////////////////////////////
       if (Blk_List.Block[i_blk].used && 
           (Blk_List.Block[i_blk].nSE == 1) &&
           (Blk_List.Block[i_blk].info.level >
            Blk_List.Block[i_blk].infoSE[0].level)) {
	  neighbour_cpu = Blk_List.Block[i_blk].infoSE[0].cpu;
	  neighbour_blk = Blk_List.Block[i_blk].infoSE[0].blknum;
          if (neighbour_cpu != i_cpu) return(4100);
          if (!Blk_List.Block[neighbour_blk].used) return(4101); 
          buffer_size = sqr(Blk_List.Block[i_blk].infoSE[0].dimen.ghost)*
                        Number_of_Solution_Variables;
          buffer_size_neighbour = sqr(Blk_List.Block[neighbour_blk].info.dimen.ghost)*
                                  Number_of_Solution_Variables;
          if (buffer_size != buffer_size_neighbour) return(4102);
          if (Blk_List.Block[i_blk].infoSE[0].dimen.i > 0 &&
              Blk_List.Block[i_blk].infoSE[0].dimen.j > 0) {
             if (Blk_List.Block[neighbour_blk].nNW != 1 ||
                 Blk_List.Block[neighbour_blk].infoNW[0].cpu != i_cpu ||
                 Blk_List.Block[neighbour_blk].infoNW[0].blknum != i_blk) return(4103);
             for (l = 0; l <= buffer_size-1; ++l) {
	        Blk_List.message_reschange_northwestcorner_recbuf[neighbour_blk][l] =
                  Blk_List.message_reschange_southeastcorner_sendbuf[i_blk][l];
             } /* endfor */
          } else if (Blk_List.Block[i_blk].infoSE[0].dimen.j > 0) {
             if (Blk_List.Block[neighbour_blk].nNE != 1 ||
                 Blk_List.Block[neighbour_blk].infoNE[0].cpu != i_cpu ||
                 Blk_List.Block[neighbour_blk].infoNE[0].blknum != i_blk) return(4104);
             for (l = 0; l <= buffer_size-1; ++l) {
	        Blk_List.message_reschange_northeastcorner_recbuf[neighbour_blk][l] =
                  Blk_List.message_reschange_southeastcorner_sendbuf[i_blk][l];
             } /* endfor */
          } else if (Blk_List.Block[i_blk].infoSE[0].dimen.i > 0) {
             if (Blk_List.Block[neighbour_blk].nSW != 1 ||
                 Blk_List.Block[neighbour_blk].infoSW[0].cpu != i_cpu ||
                 Blk_List.Block[neighbour_blk].infoSW[0].blknum != i_blk) return(4105);
             for (l = 0; l <= buffer_size-1; ++l) {
	        Blk_List.message_reschange_southwestcorner_recbuf[neighbour_blk][l] =
                  Blk_List.message_reschange_southeastcorner_sendbuf[i_blk][l];
             } /* endfor */
          } else {
             if (Blk_List.Block[neighbour_blk].nSE != 1 ||
                 Blk_List.Block[neighbour_blk].infoSE[0].cpu != i_cpu ||
                 Blk_List.Block[neighbour_blk].infoSE[0].blknum != i_blk) return(4106);
             for (l = 0; l <= buffer_size-1; ++l) {
	        Blk_List.message_reschange_southeastcorner_recbuf[neighbour_blk][l] =
                  Blk_List.message_reschange_southeastcorner_sendbuf[i_blk][l];
             } /* endfor */
          } /* endif */
       } /* endif */

       //////////////////////////////////////////////////////////////
       // Exchange messages at the South West corners of solution  //
       // blocks with coarser neighbour at lower mesh resolution.  //
       //////////////////////////////////////////////////////////////
       if (Blk_List.Block[i_blk].used && 
           (Blk_List.Block[i_blk].nSW == 1) &&
           (Blk_List.Block[i_blk].info.level >
            Blk_List.Block[i_blk].infoSW[0].level)) {
	  neighbour_cpu = Blk_List.Block[i_blk].infoSW[0].cpu;
	  neighbour_blk = Blk_List.Block[i_blk].infoSW[0].blknum;
          if (neighbour_cpu != i_cpu) return(5100);
          if (!Blk_List.Block[neighbour_blk].used) return(5101); 
          buffer_size = sqr(Blk_List.Block[i_blk].infoSW[0].dimen.ghost)*
                        Number_of_Solution_Variables;
          buffer_size_neighbour = sqr(Blk_List.Block[neighbour_blk].info.dimen.ghost)*
                                  Number_of_Solution_Variables;
          if (buffer_size != buffer_size_neighbour) return(5102);
          if (Blk_List.Block[i_blk].infoSW[0].dimen.i > 0 &&
              Blk_List.Block[i_blk].infoSW[0].dimen.j > 0) {
             if (Blk_List.Block[neighbour_blk].nNE != 1 ||
                 Blk_List.Block[neighbour_blk].infoNE[0].cpu != i_cpu ||
                 Blk_List.Block[neighbour_blk].infoNE[0].blknum != i_blk) return(5103);
             for (l = 0; l <= buffer_size-1; ++l) {
	        Blk_List.message_reschange_northeastcorner_recbuf[neighbour_blk][l] =
                  Blk_List.message_reschange_southwestcorner_sendbuf[i_blk][l];
             } /* endfor */
          } else if (Blk_List.Block[i_blk].infoSW[0].dimen.j > 0) {
             if (Blk_List.Block[neighbour_blk].nNW != 1 ||
                 Blk_List.Block[neighbour_blk].infoNW[0].cpu != i_cpu ||
                 Blk_List.Block[neighbour_blk].infoNW[0].blknum != i_blk) return(5104);
             for (l = 0; l <= buffer_size-1; ++l) {
	        Blk_List.message_reschange_northwestcorner_recbuf[neighbour_blk][l] =
                  Blk_List.message_reschange_southwestcorner_sendbuf[i_blk][l];
             } /* endfor */
          } else if (Blk_List.Block[i_blk].infoSW[0].dimen.i > 0) {
             if (Blk_List.Block[neighbour_blk].nSE != 1 ||
                 Blk_List.Block[neighbour_blk].infoSE[0].cpu != i_cpu ||
                 Blk_List.Block[neighbour_blk].infoSE[0].blknum != i_blk) return(5105);
             for (l = 0; l <= buffer_size-1; ++l) {
	        Blk_List.message_reschange_southeastcorner_recbuf[neighbour_blk][l] =
                  Blk_List.message_reschange_southwestcorner_sendbuf[i_blk][l];
             } /* endfor */
          } else {
             if (Blk_List.Block[neighbour_blk].nSW != 1 ||
                 Blk_List.Block[neighbour_blk].infoSW[0].cpu != i_cpu ||
                 Blk_List.Block[neighbour_blk].infoSW[0].blknum != i_blk) return(5106);
             for (l = 0; l <= buffer_size-1; ++l) {
	        Blk_List.message_reschange_southwestcorner_recbuf[neighbour_blk][l] =
                  Blk_List.message_reschange_southwestcorner_sendbuf[i_blk][l];
             } /* endfor */
          } /* endif */
       } /* endif */

    }  /* endfor */

    /* Message passing complete.  Return zero value. */

    return(0);

}

/**********************************************************
 * Routine: Exchange_Messages_ResChange_CoarseToFine      *
 *                                                        *
 * Sends and receives solution information contained in   *
 * the preloaded message passing buffers between          *
 * neighbouring adaptive blocks, sending solution         *
 * information from more coarse (lower mesh resolution)   *
 * adaptive blocks to more refined adaptive solution      *
 * blocks with higher mesh resolution.                    *
 *                                                        *
 **********************************************************/
int Exchange_Messages_ResChange_CoarseToFine(AdaptiveBlock2D_List &Blk_List,
                                             const int Number_of_Solution_Variables) {

    int i_cpu, i_blk, neighbour_cpu, neighbour_blk, 
        buffer_size, buffer_size_neighbour, l, j_neigh;

    /* Perform message passing for the block faces and corners of each 
       more coarse (lower mesh resolution) solution block having more 
       refined (higher mesh resolution) neighbouring solution blocks. */

    i_cpu = Blk_List.ThisCPU;

    for ( i_blk = 0 ; i_blk <= Blk_List.Nblk-1 ; ++i_blk ) {

       //////////////////////////////////////////////////////////////
       // Exchange messages at the North faces of solution blocks  //
       // with finer neighbours at higher mesh resolution.         //
       //////////////////////////////////////////////////////////////
       if (Blk_List.Block[i_blk].used && 
           (Blk_List.Block[i_blk].nN == 2) &&
           (Blk_List.Block[i_blk].info.level <
            Blk_List.Block[i_blk].infoN[0].level)) {
          for (j_neigh = 0; j_neigh <= Blk_List.Block[i_blk].nN-1 ; ++j_neigh) {
             neighbour_cpu = Blk_List.Block[i_blk].infoN[j_neigh].cpu;
             neighbour_blk = Blk_List.Block[i_blk].infoN[j_neigh].blknum;
             if (neighbour_cpu != i_cpu) return(2200);
             if (!Blk_List.Block[neighbour_blk].used) return(2201);
             buffer_size = Blk_List.Block[i_blk].infoN[j_neigh].dimen.ghost*
                           (abs(Blk_List.Block[i_blk].infoN[j_neigh].dimen.i)+
                           Blk_List.Block[i_blk].infoN[j_neigh].dimen.ghost+1)*
                           Number_of_Solution_Variables+
                           2*Blk_List.Block[i_blk].infoN[j_neigh].dimen.ghost;
             buffer_size_neighbour = Blk_List.Block[neighbour_blk].info.dimen.ghost*
                                     (abs(Blk_List.Block[neighbour_blk].info.dimen.i)+
                                     Blk_List.Block[neighbour_blk].info.dimen.ghost+1)*
 	                             Number_of_Solution_Variables+
                                     2*Blk_List.Block[neighbour_blk].info.dimen.ghost;
             if (buffer_size != buffer_size_neighbour) return(2202);
             if (Blk_List.Block[i_blk].infoN[j_neigh].dimen.j > 0) {
                if (Blk_List.Block[neighbour_blk].nS != 1 ||
                    Blk_List.Block[neighbour_blk].infoS[0].cpu != i_cpu ||
                    Blk_List.Block[neighbour_blk].infoS[0].blknum != i_blk) return(2203);
                for (l = 0; l <= buffer_size-1; ++l) {
	           Blk_List.message_reschange_southface_recbuf[neighbour_blk][0][l] =
                     Blk_List.message_reschange_northface_sendbuf[i_blk][j_neigh][l];
                } /* endfor */
             } else {
                if (Blk_List.Block[neighbour_blk].nN != 1 ||
                    Blk_List.Block[neighbour_blk].infoN[0].cpu != i_cpu ||
                    Blk_List.Block[neighbour_blk].infoN[0].blknum != i_blk) return(2204);
                for (l = 0; l <= buffer_size-1; ++l) {
	           Blk_List.message_reschange_northface_recbuf[neighbour_blk][0][l] =
                     Blk_List.message_reschange_northface_sendbuf[i_blk][j_neigh][l];
                } /* endfor */
             } /* endif */
          } /* endfor */
       } /* endif */

       //////////////////////////////////////////////////////////////
       // Exchange messages at the South faces of solution blocks  //
       // with finer neighbours at higher mesh resolution.         //
       //////////////////////////////////////////////////////////////
       if (Blk_List.Block[i_blk].used && 
           (Blk_List.Block[i_blk].nS == 2) &&
           (Blk_List.Block[i_blk].info.level <
            Blk_List.Block[i_blk].infoS[0].level)) {
          for (j_neigh = 0; j_neigh <= Blk_List.Block[i_blk].nS-1 ; ++j_neigh) {
             neighbour_cpu = Blk_List.Block[i_blk].infoS[j_neigh].cpu;
             neighbour_blk = Blk_List.Block[i_blk].infoS[j_neigh].blknum;
             if (neighbour_cpu != i_cpu) return(2100);
             if (!Blk_List.Block[neighbour_blk].used) return(2101);
             buffer_size = Blk_List.Block[i_blk].infoS[j_neigh].dimen.ghost*
                           (abs(Blk_List.Block[i_blk].infoS[j_neigh].dimen.i)+
                           Blk_List.Block[i_blk].infoS[j_neigh].dimen.ghost+1)*
                           Number_of_Solution_Variables+
                           2*Blk_List.Block[i_blk].infoS[j_neigh].dimen.ghost;
             buffer_size_neighbour = Blk_List.Block[neighbour_blk].info.dimen.ghost*
                                     (abs(Blk_List.Block[neighbour_blk].info.dimen.i)+
                                     Blk_List.Block[neighbour_blk].info.dimen.ghost+1)*
 	                             Number_of_Solution_Variables+
                                     2*Blk_List.Block[neighbour_blk].info.dimen.ghost;
             if (buffer_size != buffer_size_neighbour) return(2102);
             if (Blk_List.Block[i_blk].infoS[j_neigh].dimen.j > 0) {
                if (Blk_List.Block[neighbour_blk].nN != 1 ||
                    Blk_List.Block[neighbour_blk].infoN[0].cpu != i_cpu ||
                    Blk_List.Block[neighbour_blk].infoN[0].blknum != i_blk) return(2103);
                for (l = 0; l <= buffer_size-1; ++l) {
	           Blk_List.message_reschange_northface_recbuf[neighbour_blk][0][l] =
                     Blk_List.message_reschange_southface_sendbuf[i_blk][j_neigh][l];
                } /* endfor */
             } else {
                if (Blk_List.Block[neighbour_blk].nS != 1 ||
                    Blk_List.Block[neighbour_blk].infoS[0].cpu != i_cpu ||
                    Blk_List.Block[neighbour_blk].infoS[0].blknum != i_blk) return(2104);
                for (l = 0; l <= buffer_size-1; ++l) {
	           Blk_List.message_reschange_southface_recbuf[neighbour_blk][0][l] =
                     Blk_List.message_reschange_southface_sendbuf[i_blk][j_neigh][l];
                } /* endfor */
             } /* endif */
          } /* endfor */
       } /* endif */

       //////////////////////////////////////////////////////////////
       // Exchange messages at the East faces of solution blocks   //
       // with finer neighbours at higher mesh resolution.         //
       //////////////////////////////////////////////////////////////
       if (Blk_List.Block[i_blk].used && 
           (Blk_List.Block[i_blk].nE == 2) &&
           (Blk_List.Block[i_blk].info.level <
            Blk_List.Block[i_blk].infoE[0].level)) {
          for (j_neigh = 0; j_neigh <= Blk_List.Block[i_blk].nE-1 ; ++j_neigh) {
             neighbour_cpu = Blk_List.Block[i_blk].infoE[j_neigh].cpu;
             neighbour_blk = Blk_List.Block[i_blk].infoE[j_neigh].blknum;
             if (neighbour_cpu != i_cpu) return(1200);
             if (!Blk_List.Block[neighbour_blk].used) return(1201);
             buffer_size = Blk_List.Block[i_blk].infoE[j_neigh].dimen.ghost*
                           (abs(Blk_List.Block[i_blk].infoE[j_neigh].dimen.j)+
                           Blk_List.Block[i_blk].infoE[j_neigh].dimen.ghost+1)*
                           Number_of_Solution_Variables+
                           2*Blk_List.Block[i_blk].infoE[j_neigh].dimen.ghost;
             buffer_size_neighbour = Blk_List.Block[neighbour_blk].info.dimen.ghost*
                                     (abs(Blk_List.Block[neighbour_blk].info.dimen.j)+
                                     Blk_List.Block[neighbour_blk].info.dimen.ghost+1)*
 	                             Number_of_Solution_Variables+
                                     2*Blk_List.Block[neighbour_blk].info.dimen.ghost;
             if (buffer_size != buffer_size_neighbour) return(1202);
             if (Blk_List.Block[i_blk].infoE[j_neigh].dimen.j > 0) {
                if (Blk_List.Block[neighbour_blk].nW != 1 ||
                    Blk_List.Block[neighbour_blk].infoW[0].cpu != i_cpu ||
                    Blk_List.Block[neighbour_blk].infoW[0].blknum != i_blk) return(1203);
                for (l = 0; l <= buffer_size-1; ++l) {
	           Blk_List.message_reschange_westface_recbuf[neighbour_blk][0][l] =
                     Blk_List.message_reschange_eastface_sendbuf[i_blk][j_neigh][l];
                } /* endfor */
             } else {
                if (Blk_List.Block[neighbour_blk].nE != 1 ||
                    Blk_List.Block[neighbour_blk].infoE[0].cpu != i_cpu ||
                    Blk_List.Block[neighbour_blk].infoE[0].blknum != i_blk) return(1204);
                for (l = 0; l <= buffer_size-1; ++l) {
	           Blk_List.message_reschange_eastface_recbuf[neighbour_blk][0][l] =
                     Blk_List.message_reschange_eastface_sendbuf[i_blk][j_neigh][l];
                } /* endfor */
             } /* endif */
          } /* endfor */
       } /* endif */

       //////////////////////////////////////////////////////////////
       // Exchange messages at the West faces of solution blocks   //
       // with finer neighbours at higher mesh resolution.         //
       //////////////////////////////////////////////////////////////
       if (Blk_List.Block[i_blk].used && 
           (Blk_List.Block[i_blk].nW == 2) &&
           (Blk_List.Block[i_blk].info.level <
            Blk_List.Block[i_blk].infoW[0].level)) {
          for (j_neigh = 0; j_neigh <= Blk_List.Block[i_blk].nW-1 ; ++j_neigh) {
             neighbour_cpu = Blk_List.Block[i_blk].infoW[j_neigh].cpu;
             neighbour_blk = Blk_List.Block[i_blk].infoW[j_neigh].blknum;
             if (neighbour_cpu != i_cpu) return(1100);
             if (!Blk_List.Block[neighbour_blk].used) return(1101);
             buffer_size = Blk_List.Block[i_blk].infoW[j_neigh].dimen.ghost*
                           (abs(Blk_List.Block[i_blk].infoW[j_neigh].dimen.j)+
                           Blk_List.Block[i_blk].infoW[j_neigh].dimen.ghost+1)*
                           Number_of_Solution_Variables+
                           2*Blk_List.Block[i_blk].infoW[j_neigh].dimen.ghost;
             buffer_size_neighbour = Blk_List.Block[neighbour_blk].info.dimen.ghost*
                                     (abs(Blk_List.Block[neighbour_blk].info.dimen.j)+
                                     Blk_List.Block[neighbour_blk].info.dimen.ghost+1)*
 	                             Number_of_Solution_Variables+
                                     2*Blk_List.Block[neighbour_blk].info.dimen.ghost;
             if (buffer_size != buffer_size_neighbour) return(1102);
             if (Blk_List.Block[i_blk].infoW[j_neigh].dimen.j > 0) {
                if (Blk_List.Block[neighbour_blk].nE != 1 ||
                    Blk_List.Block[neighbour_blk].infoE[0].cpu != i_cpu ||
                    Blk_List.Block[neighbour_blk].infoE[0].blknum != i_blk) return(1103);
                for (l = 0; l <= buffer_size-1; ++l) {
	           Blk_List.message_reschange_eastface_recbuf[neighbour_blk][0][l] =
                     Blk_List.message_reschange_westface_sendbuf[i_blk][j_neigh][l];
                } /* endfor */
             } else {
                if (Blk_List.Block[neighbour_blk].nW != 1 ||
                    Blk_List.Block[neighbour_blk].infoW[0].cpu != i_cpu ||
                    Blk_List.Block[neighbour_blk].infoW[0].blknum != i_blk) return(1104);
                for (l = 0; l <= buffer_size-1; ++l) {
	           Blk_List.message_reschange_westface_recbuf[neighbour_blk][0][l] =
                     Blk_List.message_reschange_westface_sendbuf[i_blk][j_neigh][l];
                } /* endfor */
             } /* endif */
          } /* endfor */
       } /* endif */

       /////////////////////////////////////////////////////////////////////
       // Exchange messages at the North West corners of solution blocks  //
       // with finer neighbours at higher mesh resolution.                //
       /////////////////////////////////////////////////////////////////////
       if (Blk_List.Block[i_blk].used && 
           (Blk_List.Block[i_blk].nNW == 1) &&
           (Blk_List.Block[i_blk].info.level <
            Blk_List.Block[i_blk].infoNW[0].level)) {
	  neighbour_cpu = Blk_List.Block[i_blk].infoNW[0].cpu;
	  neighbour_blk = Blk_List.Block[i_blk].infoNW[0].blknum;
          if (neighbour_cpu != i_cpu) return(4200);
          if (!Blk_List.Block[neighbour_blk].used) return(4201); 
          buffer_size = sqr(Blk_List.Block[i_blk].infoNW[0].dimen.ghost)*
                        Number_of_Solution_Variables;
          buffer_size_neighbour = sqr(Blk_List.Block[neighbour_blk].info.dimen.ghost)*
                                  Number_of_Solution_Variables;
          if (buffer_size != buffer_size_neighbour) return(4202);
          if (Blk_List.Block[i_blk].infoNW[0].dimen.i > 0 &&
              Blk_List.Block[i_blk].infoNW[0].dimen.j > 0) {
             if (Blk_List.Block[neighbour_blk].nSE != 1 ||
                 Blk_List.Block[neighbour_blk].infoSE[0].cpu != i_cpu ||
                 Blk_List.Block[neighbour_blk].infoSE[0].blknum != i_blk) return(4203);
             for (l = 0; l <= buffer_size-1; ++l) {
	        Blk_List.message_reschange_southeastcorner_recbuf[neighbour_blk][l] =
                  Blk_List.message_reschange_northwestcorner_sendbuf[i_blk][l];
             } /* endfor */
          } else if (Blk_List.Block[i_blk].infoNW[0].dimen.j > 0) {
             if (Blk_List.Block[neighbour_blk].nSW != 1 ||
                 Blk_List.Block[neighbour_blk].infoSW[0].cpu != i_cpu ||
                 Blk_List.Block[neighbour_blk].infoSW[0].blknum != i_blk) return(4204);
             for (l = 0; l <= buffer_size-1; ++l) {
	        Blk_List.message_reschange_southwestcorner_recbuf[neighbour_blk][l] =
                  Blk_List.message_reschange_northwestcorner_sendbuf[i_blk][l];
             } /* endfor */
          } else if (Blk_List.Block[i_blk].infoNW[0].dimen.i > 0) {
             if (Blk_List.Block[neighbour_blk].nNE != 1 ||
                 Blk_List.Block[neighbour_blk].infoNE[0].cpu != i_cpu ||
                 Blk_List.Block[neighbour_blk].infoNE[0].blknum != i_blk) return(4205);
             for (l = 0; l <= buffer_size-1; ++l) {
	        Blk_List.message_reschange_northeastcorner_recbuf[neighbour_blk][l] =
                  Blk_List.message_reschange_northwestcorner_sendbuf[i_blk][l];
             } /* endfor */
          } else {
             if (Blk_List.Block[neighbour_blk].nNW != 1 ||
                 Blk_List.Block[neighbour_blk].infoNW[0].cpu != i_cpu ||
                 Blk_List.Block[neighbour_blk].infoNW[0].blknum != i_blk) return(4206);
             for (l = 0; l <= buffer_size-1; ++l) {
	        Blk_List.message_reschange_northwestcorner_recbuf[neighbour_blk][l] =
                  Blk_List.message_reschange_northwestcorner_sendbuf[i_blk][l];
             } /* endfor */
          } /* endif */
       } /* endif */

       /////////////////////////////////////////////////////////////////////
       // Exchange messages at the North East corners of solution blocks  //
       // with finer neighbours at higher mesh resolution.                //
       /////////////////////////////////////////////////////////////////////
       if (Blk_List.Block[i_blk].used && 
           (Blk_List.Block[i_blk].nNE == 1) &&
           (Blk_List.Block[i_blk].info.level <
            Blk_List.Block[i_blk].infoNE[0].level)) {
	  neighbour_cpu = Blk_List.Block[i_blk].infoNE[0].cpu;
	  neighbour_blk = Blk_List.Block[i_blk].infoNE[0].blknum;
          if (neighbour_cpu != i_cpu) return(5200);
          if (!Blk_List.Block[neighbour_blk].used) return(5201); 
          buffer_size = sqr(Blk_List.Block[i_blk].infoNE[0].dimen.ghost)*
                        Number_of_Solution_Variables;
          buffer_size_neighbour = sqr(Blk_List.Block[neighbour_blk].info.dimen.ghost)*
                                  Number_of_Solution_Variables;
          if (buffer_size != buffer_size_neighbour) return(5202);
          if (Blk_List.Block[i_blk].infoNE[0].dimen.i > 0 &&
              Blk_List.Block[i_blk].infoNE[0].dimen.j > 0) {
             if (Blk_List.Block[neighbour_blk].nSW != 1 ||
                 Blk_List.Block[neighbour_blk].infoSW[0].cpu != i_cpu ||
                 Blk_List.Block[neighbour_blk].infoSW[0].blknum != i_blk) return(5203);
             for (l = 0; l <= buffer_size-1; ++l) {
	        Blk_List.message_reschange_southwestcorner_recbuf[neighbour_blk][l] =
                  Blk_List.message_reschange_northeastcorner_sendbuf[i_blk][l];
             } /* endfor */
          } else if (Blk_List.Block[i_blk].infoNE[0].dimen.j > 0) {
             if (Blk_List.Block[neighbour_blk].nSE != 1 ||
                 Blk_List.Block[neighbour_blk].infoSE[0].cpu != i_cpu ||
                 Blk_List.Block[neighbour_blk].infoSE[0].blknum != i_blk) return(5204);
             for (l = 0; l <= buffer_size-1; ++l) {
	        Blk_List.message_reschange_southeastcorner_recbuf[neighbour_blk][l] =
                  Blk_List.message_reschange_northeastcorner_sendbuf[i_blk][l];
             } /* endfor */
          } else if (Blk_List.Block[i_blk].infoNE[0].dimen.i > 0) {
             if (Blk_List.Block[neighbour_blk].nNW != 1 ||
                 Blk_List.Block[neighbour_blk].infoNW[0].cpu != i_cpu ||
                 Blk_List.Block[neighbour_blk].infoNW[0].blknum != i_blk) return(5205);
             for (l = 0; l <= buffer_size-1; ++l) {
	        Blk_List.message_reschange_northwestcorner_recbuf[neighbour_blk][l] =
                  Blk_List.message_reschange_northeastcorner_sendbuf[i_blk][l];
             } /* endfor */
          } else {
             if (Blk_List.Block[neighbour_blk].nNE != 1 ||
                 Blk_List.Block[neighbour_blk].infoNE[0].cpu != i_cpu ||
                 Blk_List.Block[neighbour_blk].infoNE[0].blknum != i_blk) return(5206);
             for (l = 0; l <= buffer_size-1; ++l) {
	        Blk_List.message_reschange_northeastcorner_recbuf[neighbour_blk][l] =
                  Blk_List.message_reschange_northeastcorner_sendbuf[i_blk][l];
             } /* endfor */
          } /* endif */
       } /* endif */

       /////////////////////////////////////////////////////////////////////
       // Exchange messages at the South East corners of solution blocks  //
       // with finer neighbours at higher mesh resolution.                //
       /////////////////////////////////////////////////////////////////////
       if (Blk_List.Block[i_blk].used && 
           (Blk_List.Block[i_blk].nSE == 1) &&
           (Blk_List.Block[i_blk].info.level <
            Blk_List.Block[i_blk].infoSE[0].level)) {
	  neighbour_cpu = Blk_List.Block[i_blk].infoSE[0].cpu;
	  neighbour_blk = Blk_List.Block[i_blk].infoSE[0].blknum;
          if (neighbour_cpu != i_cpu) return(4100);
          if (!Blk_List.Block[neighbour_blk].used) return(4101); 
          buffer_size = sqr(Blk_List.Block[i_blk].infoSE[0].dimen.ghost)*
                        Number_of_Solution_Variables;
          buffer_size_neighbour = sqr(Blk_List.Block[neighbour_blk].info.dimen.ghost)*
                                  Number_of_Solution_Variables;
          if (buffer_size != buffer_size_neighbour) return(4102);
          if (Blk_List.Block[i_blk].infoSE[0].dimen.i > 0 &&
              Blk_List.Block[i_blk].infoSE[0].dimen.j > 0) {
             if (Blk_List.Block[neighbour_blk].nNW != 1 ||
                 Blk_List.Block[neighbour_blk].infoNW[0].cpu != i_cpu ||
                 Blk_List.Block[neighbour_blk].infoNW[0].blknum != i_blk) return(4103);
             for (l = 0; l <= buffer_size-1; ++l) {
	        Blk_List.message_reschange_northwestcorner_recbuf[neighbour_blk][l] =
                  Blk_List.message_reschange_southeastcorner_sendbuf[i_blk][l];
             } /* endfor */
          } else if (Blk_List.Block[i_blk].infoSE[0].dimen.j > 0) {
             if (Blk_List.Block[neighbour_blk].nNE != 1 ||
                 Blk_List.Block[neighbour_blk].infoNE[0].cpu != i_cpu ||
                 Blk_List.Block[neighbour_blk].infoNE[0].blknum != i_blk) return(4104);
             for (l = 0; l <= buffer_size-1; ++l) {
	        Blk_List.message_reschange_northeastcorner_recbuf[neighbour_blk][l] =
                  Blk_List.message_reschange_southeastcorner_sendbuf[i_blk][l];
             } /* endfor */
          } else if (Blk_List.Block[i_blk].infoSE[0].dimen.i > 0) {
             if (Blk_List.Block[neighbour_blk].nSW != 1 ||
                 Blk_List.Block[neighbour_blk].infoSW[0].cpu != i_cpu ||
                 Blk_List.Block[neighbour_blk].infoSW[0].blknum != i_blk) return(4105);
             for (l = 0; l <= buffer_size-1; ++l) {
	        Blk_List.message_reschange_southwestcorner_recbuf[neighbour_blk][l] =
                  Blk_List.message_reschange_southeastcorner_sendbuf[i_blk][l];
             } /* endfor */
          } else {
             if (Blk_List.Block[neighbour_blk].nSE != 1 ||
                 Blk_List.Block[neighbour_blk].infoSE[0].cpu != i_cpu ||
                 Blk_List.Block[neighbour_blk].infoSE[0].blknum != i_blk) return(4106);
             for (l = 0; l <= buffer_size-1; ++l) {
	        Blk_List.message_reschange_southeastcorner_recbuf[neighbour_blk][l] =
                  Blk_List.message_reschange_southeastcorner_sendbuf[i_blk][l];
             } /* endfor */
          } /* endif */
       } /* endif */

       /////////////////////////////////////////////////////////////////////
       // Exchange messages at the South West corners of solution blocks  //
       // with finer neighbours at higher mesh resolution.                //
       /////////////////////////////////////////////////////////////////////
       if (Blk_List.Block[i_blk].used && 
           (Blk_List.Block[i_blk].nSW == 1) &&
           (Blk_List.Block[i_blk].info.level <
            Blk_List.Block[i_blk].infoSW[0].level)) {
	  neighbour_cpu = Blk_List.Block[i_blk].infoSW[0].cpu;
	  neighbour_blk = Blk_List.Block[i_blk].infoSW[0].blknum;
          if (neighbour_cpu != i_cpu) return(5100);
          if (!Blk_List.Block[neighbour_blk].used) return(5101); 
          buffer_size = sqr(Blk_List.Block[i_blk].infoSW[0].dimen.ghost)*
                        Number_of_Solution_Variables;
          buffer_size_neighbour = sqr(Blk_List.Block[neighbour_blk].info.dimen.ghost)*
                                  Number_of_Solution_Variables;
          if (buffer_size != buffer_size_neighbour) return(5102);
          if (Blk_List.Block[i_blk].infoSW[0].dimen.i > 0 &&
              Blk_List.Block[i_blk].infoSW[0].dimen.j > 0) {
             if (Blk_List.Block[neighbour_blk].nNE != 1 ||
                 Blk_List.Block[neighbour_blk].infoNE[0].cpu != i_cpu ||
                 Blk_List.Block[neighbour_blk].infoNE[0].blknum != i_blk) return(5103);
             for (l = 0; l <= buffer_size-1; ++l) {
	        Blk_List.message_reschange_northeastcorner_recbuf[neighbour_blk][l] =
                  Blk_List.message_reschange_southwestcorner_sendbuf[i_blk][l];
             } /* endfor */
          } else if (Blk_List.Block[i_blk].infoSW[0].dimen.j > 0) {
             if (Blk_List.Block[neighbour_blk].nNW != 1 ||
                 Blk_List.Block[neighbour_blk].infoNW[0].cpu != i_cpu ||
                 Blk_List.Block[neighbour_blk].infoNW[0].blknum != i_blk) return(5104);
             for (l = 0; l <= buffer_size-1; ++l) {
	        Blk_List.message_reschange_northwestcorner_recbuf[neighbour_blk][l] =
                  Blk_List.message_reschange_southwestcorner_sendbuf[i_blk][l];
             } /* endfor */
          } else if (Blk_List.Block[i_blk].infoSW[0].dimen.i > 0) {
             if (Blk_List.Block[neighbour_blk].nSE != 1 ||
                 Blk_List.Block[neighbour_blk].infoSE[0].cpu != i_cpu ||
                 Blk_List.Block[neighbour_blk].infoSE[0].blknum != i_blk) return(5105);
             for (l = 0; l <= buffer_size-1; ++l) {
	        Blk_List.message_reschange_southeastcorner_recbuf[neighbour_blk][l] =
                  Blk_List.message_reschange_southwestcorner_sendbuf[i_blk][l];
             } /* endfor */
          } else {
             if (Blk_List.Block[neighbour_blk].nSW != 1 ||
                 Blk_List.Block[neighbour_blk].infoSW[0].cpu != i_cpu ||
                 Blk_List.Block[neighbour_blk].infoSW[0].blknum != i_blk) return(5106);
             for (l = 0; l <= buffer_size-1; ++l) {
	        Blk_List.message_reschange_southwestcorner_recbuf[neighbour_blk][l] =
                  Blk_List.message_reschange_southwestcorner_sendbuf[i_blk][l];
             } /* endfor */
          } /* endif */
       } /* endif */

    }  /* endfor */

    /* Message passing complete.  Return zero value. */

    return(0);

}
