/* AdaptiveBlock_MPI.cc:  Subroutines for adaptive blocks classes
                          on a multi-processor cluster computer 
                          with access to MPI. */

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

    int number_receive_requests, 
        number_send_requests, 
        tag_receive, tag_send;
    int tag_N = 22, tag_S = 21, tag_E = 12, tag_W = 11, 
        tag_NW = 42, tag_NE = 52, tag_SE = 41, tag_SW = 51,
        tag_base = 100;

    MPI::Request *receive_requests, 
                 *send_requests;
           

    /* Allocate memory for message passing requests. */

    receive_requests = new MPI::Request[8*Blk_List.Nblk];
    send_requests = new MPI::Request[8*Blk_List.Nblk];
    number_receive_requests = 0;
    number_send_requests = 0;

    /* Perform message passing for the block faces and corners of each 
       solution block having no mesh resolution change. */

    i_cpu = Blk_List.ThisCPU;

    for ( i_blk = 0 ; i_blk < Blk_List.Nblk ; ++i_blk ) {

       //////////////////////////////////////////////////////////////
       // Exchange messages at the North faces of solution blocks. //
       //////////////////////////////////////////////////////////////
       if (Blk_List.Block[i_blk].used && 
           (Blk_List.Block[i_blk].nN == 1) &&
           (Blk_List.Block[i_blk].info.level == Blk_List.Block[i_blk].infoN[0].level)) {
	  neighbour_cpu = Blk_List.Block[i_blk].infoN[0].cpu;
	  neighbour_blk = Blk_List.Block[i_blk].infoN[0].blknum;
          buffer_size = Blk_List.Block[i_blk].infoN[0].dimen.ghost*
                        (abs(Blk_List.Block[i_blk].infoN[0].dimen.i)+1)*
                        Number_of_Solution_Variables+
                        2*Blk_List.Block[i_blk].infoN[0].dimen.ghost;
          if (neighbour_cpu != i_cpu) {
             if (Blk_List.Block[i_blk].infoN[0].dimen.j > 0) {
                tag_receive = neighbour_blk*tag_base + tag_N;
                tag_send    = i_blk*tag_base + tag_S;
             } else {
                tag_receive = neighbour_blk*tag_base + tag_N;
                tag_send    = i_blk*tag_base + tag_N;
             } /* endif */
             receive_requests[number_receive_requests] = MPI::COMM_WORLD.Irecv(Blk_List.message_noreschange_northface_recbuf[i_blk],
									       buffer_size,
									       MPI::DOUBLE,
									       neighbour_cpu,
									       tag_receive);
             number_receive_requests++;
             send_requests[number_send_requests] = MPI::COMM_WORLD.Isend(Blk_List.message_noreschange_northface_sendbuf[i_blk],
									 buffer_size,
									 MPI::DOUBLE,
									 neighbour_cpu,
									 tag_send);
             number_send_requests++;
          } else {
             if (!Blk_List.Block[neighbour_blk].used) return(2200); 
             buffer_size_neighbour = Blk_List.Block[neighbour_blk].info.dimen.ghost*
                                     (abs(Blk_List.Block[neighbour_blk].info.dimen.i)+1)*
                                     Number_of_Solution_Variables+
                                     2*Blk_List.Block[neighbour_blk].info.dimen.ghost;
             if (buffer_size != buffer_size_neighbour) return(2201);
             if (Blk_List.Block[i_blk].infoN[0].dimen.j > 0) {
                if (Blk_List.Block[neighbour_blk].nS != 1 ||
                    Blk_List.Block[neighbour_blk].infoS[0].cpu != i_cpu ||
                    Blk_List.Block[neighbour_blk].infoS[0].blknum != i_blk) return(2202);
                for (l = 0; l <= buffer_size-1; ++l) {
	           Blk_List.message_noreschange_southface_recbuf[neighbour_blk][l] =
                     Blk_List.message_noreschange_northface_sendbuf[i_blk][l];
                } /* endfor */
             } else {
                if (Blk_List.Block[neighbour_blk].nN != 1 ||
                    Blk_List.Block[neighbour_blk].infoN[0].cpu != i_cpu ||
                    Blk_List.Block[neighbour_blk].infoN[0].blknum != i_blk) return(2203);
                for (l = 0; l <= buffer_size-1; ++l) {
	           Blk_List.message_noreschange_northface_recbuf[neighbour_blk][l] =
                     Blk_List.message_noreschange_northface_sendbuf[i_blk][l];
                } /* endfor */
             } /* endif */
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
          buffer_size = Blk_List.Block[i_blk].infoS[0].dimen.ghost*
                        (abs(Blk_List.Block[i_blk].infoS[0].dimen.i)+1)*
                        Number_of_Solution_Variables+
                        2*Blk_List.Block[i_blk].infoS[0].dimen.ghost;
          if (neighbour_cpu != i_cpu) {
             if (Blk_List.Block[i_blk].infoS[0].dimen.j > 0) {
                tag_receive = neighbour_blk*tag_base + tag_S;
                tag_send    = i_blk*tag_base + tag_N;
             } else {
                tag_receive = neighbour_blk*tag_base + tag_S;
                tag_send    = i_blk*tag_base + tag_S;
             } /* endif */
             receive_requests[number_receive_requests] = MPI::COMM_WORLD.Irecv(Blk_List.message_noreschange_southface_recbuf[i_blk],
									       buffer_size,
									       MPI::DOUBLE,
									       neighbour_cpu,
									       tag_receive);
             number_receive_requests++;             
             send_requests[number_send_requests] = MPI::COMM_WORLD.Isend(Blk_List.message_noreschange_southface_sendbuf[i_blk],
									 buffer_size,
									 MPI::DOUBLE,
									 neighbour_cpu,
									 tag_send);

             number_send_requests++;
          } else {
             if (!Blk_List.Block[neighbour_blk].used) return(2100); 
             buffer_size_neighbour = Blk_List.Block[neighbour_blk].info.dimen.ghost*
                                     (abs(Blk_List.Block[neighbour_blk].info.dimen.i)+1)*
                                     Number_of_Solution_Variables+
                                     2*Blk_List.Block[neighbour_blk].info.dimen.ghost;
             if (buffer_size != buffer_size_neighbour) return(2101);
             if (Blk_List.Block[i_blk].infoS[0].dimen.j > 0) {
                if (Blk_List.Block[neighbour_blk].nN != 1 ||
                    Blk_List.Block[neighbour_blk].infoN[0].cpu != i_cpu ||
                    Blk_List.Block[neighbour_blk].infoN[0].blknum != i_blk) return(2102);
                for (l = 0; l <= buffer_size-1; ++l) {
	           Blk_List.message_noreschange_northface_recbuf[neighbour_blk][l] =
                     Blk_List.message_noreschange_southface_sendbuf[i_blk][l];
                } /* endfor */
             } else {
                if (Blk_List.Block[neighbour_blk].nS != 1 ||
                    Blk_List.Block[neighbour_blk].infoS[0].cpu != i_cpu ||
                    Blk_List.Block[neighbour_blk].infoS[0].blknum != i_blk) return(2103);
                for (l = 0; l <= buffer_size-1; ++l) {
	           Blk_List.message_noreschange_southface_recbuf[neighbour_blk][l] =
                     Blk_List.message_noreschange_southface_sendbuf[i_blk][l];
                } /* endfor */
             } /* endif */
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
          buffer_size = Blk_List.Block[i_blk].infoE[0].dimen.ghost*
                        (abs(Blk_List.Block[i_blk].infoE[0].dimen.j)+1)*
                        Number_of_Solution_Variables+
                        2*Blk_List.Block[i_blk].infoE[0].dimen.ghost;
          if (neighbour_cpu != i_cpu) {
             if (Blk_List.Block[i_blk].infoE[0].dimen.i > 0) {
                tag_receive = neighbour_blk*tag_base + tag_E;
                tag_send    = i_blk*tag_base + tag_W;
             } else {
                tag_receive = neighbour_blk*tag_base + tag_E;
                tag_send    = i_blk*tag_base + tag_E;
             } /* endif */
              receive_requests[number_receive_requests] = MPI::COMM_WORLD.Irecv(Blk_List.message_noreschange_eastface_recbuf[i_blk],
										buffer_size,
										MPI::DOUBLE,
										neighbour_cpu,
										tag_receive);
	      number_receive_requests++;
	      send_requests[number_send_requests] = MPI::COMM_WORLD.Isend(Blk_List.message_noreschange_eastface_sendbuf[i_blk],
									  buffer_size,
									  MPI::DOUBLE,
									  neighbour_cpu,
									  tag_send);
	      number_send_requests++;
          } else {
             if (!Blk_List.Block[neighbour_blk].used) return(1200); 
             buffer_size_neighbour = Blk_List.Block[neighbour_blk].info.dimen.ghost*
                                     (abs(Blk_List.Block[neighbour_blk].info.dimen.j)+1)*
                                     Number_of_Solution_Variables+
                                     2*Blk_List.Block[neighbour_blk].info.dimen.ghost;
             if (buffer_size != buffer_size_neighbour) return(1201);
             if (Blk_List.Block[i_blk].infoE[0].dimen.i > 0) {
                if (Blk_List.Block[neighbour_blk].nW != 1 ||
                    Blk_List.Block[neighbour_blk].infoW[0].cpu != i_cpu ||
                    Blk_List.Block[neighbour_blk].infoW[0].blknum != i_blk) return(1202);
                for (l = 0; l <= buffer_size-1; ++l) {
	           Blk_List.message_noreschange_westface_recbuf[neighbour_blk][l] =
                     Blk_List.message_noreschange_eastface_sendbuf[i_blk][l];
                } /* endfor */
             } else {
                if (Blk_List.Block[neighbour_blk].nE != 1 ||
                    Blk_List.Block[neighbour_blk].infoE[0].cpu != i_cpu ||
                    Blk_List.Block[neighbour_blk].infoE[0].blknum != i_blk) return(1203);
                for (l = 0; l <= buffer_size-1; ++l) {
	           Blk_List.message_noreschange_eastface_recbuf[neighbour_blk][l] =
                     Blk_List.message_noreschange_eastface_sendbuf[i_blk][l];
                } /* endfor */
             } /* endif */
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
          buffer_size = Blk_List.Block[i_blk].infoW[0].dimen.ghost*
                        (abs(Blk_List.Block[i_blk].infoW[0].dimen.j)+1)*
                        Number_of_Solution_Variables+
                        2*Blk_List.Block[i_blk].infoW[0].dimen.ghost;
          if (neighbour_cpu != i_cpu) {
             if (Blk_List.Block[i_blk].infoW[0].dimen.i > 0) {
                tag_receive = neighbour_blk*tag_base + tag_W;
                tag_send    = i_blk*tag_base + tag_E;
             } else {
                tag_receive = neighbour_blk*tag_base + tag_W;
                tag_send    = i_blk*tag_base + tag_W;
             } /* endif */
             receive_requests[number_receive_requests] = MPI::COMM_WORLD.Irecv(Blk_List.message_noreschange_westface_recbuf[i_blk],
									       buffer_size,
									       MPI::DOUBLE,
									       neighbour_cpu,
									       tag_receive);
             number_receive_requests++;	     
             send_requests[number_send_requests] = MPI::COMM_WORLD.Isend(Blk_List.message_noreschange_westface_sendbuf[i_blk],
									   buffer_size,
									   MPI::DOUBLE,
									   neighbour_cpu,
									   tag_send);
             number_send_requests++;
          } else {
             if (!Blk_List.Block[neighbour_blk].used) return(1100); 
             buffer_size_neighbour = Blk_List.Block[neighbour_blk].info.dimen.ghost*
                                     (abs(Blk_List.Block[neighbour_blk].info.dimen.j)+1)*
                                     Number_of_Solution_Variables+
                                     2*Blk_List.Block[neighbour_blk].info.dimen.ghost;
             if (buffer_size != buffer_size_neighbour) return(1101);
             if (Blk_List.Block[i_blk].infoW[0].dimen.i > 0) {
                if (Blk_List.Block[neighbour_blk].nE != 1 ||
                    Blk_List.Block[neighbour_blk].infoE[0].cpu != i_cpu ||
                    Blk_List.Block[neighbour_blk].infoE[0].blknum != i_blk) return(1102);
                for (l = 0; l <= buffer_size-1; ++l) {
	           Blk_List.message_noreschange_eastface_recbuf[neighbour_blk][l] =
                     Blk_List.message_noreschange_westface_sendbuf[i_blk][l];
                } /* endfor */
             } else {
                if (Blk_List.Block[neighbour_blk].nW != 1 ||
                    Blk_List.Block[neighbour_blk].infoW[0].cpu != i_cpu ||
                    Blk_List.Block[neighbour_blk].infoW[0].blknum != i_blk) return(1103);
                for (l = 0; l <= buffer_size-1; ++l) {
	           Blk_List.message_noreschange_westface_recbuf[neighbour_blk][l] =
                     Blk_List.message_noreschange_westface_sendbuf[i_blk][l];
                } /* endfor */
             } /* endif */
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
          buffer_size = sqr(Blk_List.Block[i_blk].infoNW[0].dimen.ghost)*
                        Number_of_Solution_Variables;
          if (neighbour_cpu != i_cpu) {
             if (Blk_List.Block[i_blk].infoNW[0].dimen.i > 0 &&
                 Blk_List.Block[i_blk].infoNW[0].dimen.j > 0) {
                tag_receive = neighbour_blk*tag_base + tag_NW;
                tag_send    = i_blk*tag_base + tag_SE;
             } else if (Blk_List.Block[i_blk].infoNW[0].dimen.j > 0) {
                tag_receive = neighbour_blk*tag_base + tag_NW;
                tag_send    = i_blk*tag_base + tag_SW;
             } else if (Blk_List.Block[i_blk].infoNW[0].dimen.i > 0) {
                tag_receive = neighbour_blk*tag_base + tag_NW;
                tag_send    = i_blk*tag_base + tag_NE;
             } else {
                tag_receive = neighbour_blk*tag_base + tag_NW;
                tag_send    = i_blk*tag_base + tag_NW;
             } /* endif */
	     receive_requests[number_receive_requests] = MPI::COMM_WORLD.Irecv(Blk_List.message_noreschange_northwestcorner_recbuf[i_blk],
									       buffer_size,
									       MPI::DOUBLE,
									       neighbour_cpu,
									       tag_receive);
             number_receive_requests++;
             send_requests[number_send_requests] = MPI::COMM_WORLD.Isend(Blk_List.message_noreschange_northwestcorner_sendbuf[i_blk],
									 buffer_size,
									 MPI::DOUBLE,
									 neighbour_cpu,
									 tag_send);
             number_send_requests++;
          } else {
             if (!Blk_List.Block[neighbour_blk].used) return(4200); 
             buffer_size_neighbour = sqr(Blk_List.Block[neighbour_blk].info.dimen.ghost)*
                                     Number_of_Solution_Variables;
             if (buffer_size != buffer_size_neighbour) return(4201);
             if (Blk_List.Block[i_blk].infoNW[0].dimen.i > 0 &&
                 Blk_List.Block[i_blk].infoNW[0].dimen.j > 0) {
                if (Blk_List.Block[neighbour_blk].nSE != 1 ||
                    Blk_List.Block[neighbour_blk].infoSE[0].cpu != i_cpu ||
                    Blk_List.Block[neighbour_blk].infoSE[0].blknum != i_blk) return(4202);
                for (l = 0; l <= buffer_size-1; ++l) {
	           Blk_List.message_noreschange_southeastcorner_recbuf[neighbour_blk][l] =
                     Blk_List.message_noreschange_northwestcorner_sendbuf[i_blk][l];
                } /* endfor */
             } else if (Blk_List.Block[i_blk].infoNW[0].dimen.j > 0) {
                if (Blk_List.Block[neighbour_blk].nSW != 1 ||
                    Blk_List.Block[neighbour_blk].infoSW[0].cpu != i_cpu ||
                    Blk_List.Block[neighbour_blk].infoSW[0].blknum != i_blk) return(4203);
                for (l = 0; l <= buffer_size-1; ++l) {
	           Blk_List.message_noreschange_southwestcorner_recbuf[neighbour_blk][l] =
                     Blk_List.message_noreschange_northwestcorner_sendbuf[i_blk][l];
                } /* endfor */
             } else if (Blk_List.Block[i_blk].infoNW[0].dimen.i > 0) {
                if (Blk_List.Block[neighbour_blk].nNE != 1 ||
                    Blk_List.Block[neighbour_blk].infoNE[0].cpu != i_cpu ||
                    Blk_List.Block[neighbour_blk].infoNE[0].blknum != i_blk) return(4204);
                for (l = 0; l <= buffer_size-1; ++l) {
	           Blk_List.message_noreschange_northeastcorner_recbuf[neighbour_blk][l] =
                     Blk_List.message_noreschange_northwestcorner_sendbuf[i_blk][l];
                } /* endfor */
             } else {
                if (Blk_List.Block[neighbour_blk].nNW != 1 ||
                    Blk_List.Block[neighbour_blk].infoNW[0].cpu != i_cpu ||
                    Blk_List.Block[neighbour_blk].infoNW[0].blknum != i_blk) return(4205);
                for (l = 0; l <= buffer_size-1; ++l) {
	           Blk_List.message_noreschange_northwestcorner_recbuf[neighbour_blk][l] =
                     Blk_List.message_noreschange_northwestcorner_sendbuf[i_blk][l];
                } /* endfor */
             } /* endif */
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
          buffer_size = sqr(Blk_List.Block[i_blk].infoNE[0].dimen.ghost)*
                        Number_of_Solution_Variables;
          if (neighbour_cpu != i_cpu) {
             if (Blk_List.Block[i_blk].infoNE[0].dimen.i > 0 &&
                 Blk_List.Block[i_blk].infoNE[0].dimen.j > 0) {
                tag_receive = neighbour_blk*tag_base + tag_NE;
                tag_send    = i_blk*tag_base + tag_SW;
             } else if (Blk_List.Block[i_blk].infoNE[0].dimen.j > 0) {
                tag_receive = neighbour_blk*tag_base + tag_NE;
                tag_send    = i_blk*tag_base + tag_SE;
             } else if (Blk_List.Block[i_blk].infoNE[0].dimen.i > 0) {
                tag_receive = neighbour_blk*tag_base + tag_NE;
                tag_send    = i_blk*tag_base + tag_NW;
             } else {
                tag_receive = neighbour_blk*tag_base + tag_NE;
                tag_send    = i_blk*tag_base + tag_NE;
             } /* endif */
             receive_requests[number_receive_requests] = MPI::COMM_WORLD.Irecv(Blk_List.message_noreschange_northeastcorner_recbuf[i_blk],
										 buffer_size,
										 MPI::DOUBLE,
										 neighbour_cpu,
										 tag_receive);
             number_receive_requests++;
             send_requests[number_send_requests] = MPI::COMM_WORLD.Isend(Blk_List.message_noreschange_northeastcorner_sendbuf[i_blk],
									  buffer_size,
									  MPI::DOUBLE,
									  neighbour_cpu,
									  tag_send);
             number_send_requests++;
          } else {
             if (!Blk_List.Block[neighbour_blk].used) return(5200); 
             buffer_size_neighbour = sqr(Blk_List.Block[neighbour_blk].info.dimen.ghost)*
                                     Number_of_Solution_Variables;
             if (buffer_size != buffer_size_neighbour) return(5201);
             if (Blk_List.Block[i_blk].infoNE[0].dimen.i > 0 &&
                 Blk_List.Block[i_blk].infoNE[0].dimen.j > 0) {
                if (Blk_List.Block[neighbour_blk].nSW != 1 ||
                    Blk_List.Block[neighbour_blk].infoSW[0].cpu != i_cpu ||
                    Blk_List.Block[neighbour_blk].infoSW[0].blknum != i_blk) return(5202);
                for (l = 0; l <= buffer_size-1; ++l) {
	           Blk_List.message_noreschange_southwestcorner_recbuf[neighbour_blk][l] =
                     Blk_List.message_noreschange_northeastcorner_sendbuf[i_blk][l];
                } /* endfor */
             } else if (Blk_List.Block[i_blk].infoNE[0].dimen.j > 0) {
                if (Blk_List.Block[neighbour_blk].nSE != 1 ||
                    Blk_List.Block[neighbour_blk].infoSE[0].cpu != i_cpu ||
                    Blk_List.Block[neighbour_blk].infoSE[0].blknum != i_blk) return(5203);
                for (l = 0; l <= buffer_size-1; ++l) {
	           Blk_List.message_noreschange_southeastcorner_recbuf[neighbour_blk][l] =
                     Blk_List.message_noreschange_northeastcorner_sendbuf[i_blk][l];
                } /* endfor */
             } else if (Blk_List.Block[i_blk].infoNE[0].dimen.i > 0) {
                if (Blk_List.Block[neighbour_blk].nNW != 1 ||
                    Blk_List.Block[neighbour_blk].infoNW[0].cpu != i_cpu ||
                    Blk_List.Block[neighbour_blk].infoNW[0].blknum != i_blk) return(5204);
                for (l = 0; l <= buffer_size-1; ++l) {
	           Blk_List.message_noreschange_northwestcorner_recbuf[neighbour_blk][l] =
                     Blk_List.message_noreschange_northeastcorner_sendbuf[i_blk][l];
                } /* endfor */
             } else {
                if (Blk_List.Block[neighbour_blk].nNE != 1 ||
                    Blk_List.Block[neighbour_blk].infoNE[0].cpu != i_cpu ||
                    Blk_List.Block[neighbour_blk].infoNE[0].blknum != i_blk) return(5205);
                for (l = 0; l <= buffer_size-1; ++l) {
	           Blk_List.message_noreschange_northeastcorner_recbuf[neighbour_blk][l] =
                     Blk_List.message_noreschange_northeastcorner_sendbuf[i_blk][l];
                } /* endfor */
             } /* endif */
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
          buffer_size = sqr(Blk_List.Block[i_blk].infoSE[0].dimen.ghost)*
                        Number_of_Solution_Variables;
          if (neighbour_cpu != i_cpu) {
             if (Blk_List.Block[i_blk].infoSE[0].dimen.i > 0 &&
                 Blk_List.Block[i_blk].infoSE[0].dimen.j > 0) {
                tag_receive = neighbour_blk*tag_base + tag_SE;
                tag_send    = i_blk*tag_base + tag_NW;
             } else if (Blk_List.Block[i_blk].infoSE[0].dimen.j > 0) {
                tag_receive = neighbour_blk*tag_base + tag_SE;
                tag_send    = i_blk*tag_base + tag_NE;
             } else if (Blk_List.Block[i_blk].infoSE[0].dimen.i > 0) {
                tag_receive = neighbour_blk*tag_base + tag_SE;
                tag_send    = i_blk*tag_base + tag_SW;
             } else {
                tag_receive = neighbour_blk*tag_base + tag_SE;
                tag_send    = i_blk*tag_base + tag_SE;
             } /* endif */
             receive_requests[number_receive_requests] = MPI::COMM_WORLD.Irecv(Blk_List.message_noreschange_southeastcorner_recbuf[i_blk],
									       buffer_size,
									       MPI::DOUBLE,
									       neighbour_cpu,
									       tag_receive);
             number_receive_requests++;
             send_requests[number_send_requests] = MPI::COMM_WORLD.Isend(Blk_List.message_noreschange_southeastcorner_sendbuf[i_blk],
									 buffer_size,
									 MPI::DOUBLE,
									 neighbour_cpu,
									 tag_send);
             number_send_requests++;
          } else {
             if (!Blk_List.Block[neighbour_blk].used) return(4100); 
             buffer_size_neighbour = sqr(Blk_List.Block[neighbour_blk].info.dimen.ghost)*
                                     Number_of_Solution_Variables;
             if (buffer_size != buffer_size_neighbour) return(4101);
             if (Blk_List.Block[i_blk].infoSE[0].dimen.i > 0 &&
                 Blk_List.Block[i_blk].infoSE[0].dimen.j > 0) {
                if (Blk_List.Block[neighbour_blk].nNW != 1 ||
                    Blk_List.Block[neighbour_blk].infoNW[0].cpu != i_cpu ||
                    Blk_List.Block[neighbour_blk].infoNW[0].blknum != i_blk) return(4102);
                for (l = 0; l <= buffer_size-1; ++l) {
	           Blk_List.message_noreschange_northwestcorner_recbuf[neighbour_blk][l] =
                     Blk_List.message_noreschange_southeastcorner_sendbuf[i_blk][l];
                } /* endfor */
             } else if (Blk_List.Block[i_blk].infoSE[0].dimen.j > 0) {
                if (Blk_List.Block[neighbour_blk].nNE != 1 ||
                    Blk_List.Block[neighbour_blk].infoNE[0].cpu != i_cpu ||
                    Blk_List.Block[neighbour_blk].infoNE[0].blknum != i_blk) return(4103);
                for (l = 0; l <= buffer_size-1; ++l) {
	           Blk_List.message_noreschange_northeastcorner_recbuf[neighbour_blk][l] =
                     Blk_List.message_noreschange_southeastcorner_sendbuf[i_blk][l];
                } /* endfor */
             } else if (Blk_List.Block[i_blk].infoSE[0].dimen.i > 0) {
                if (Blk_List.Block[neighbour_blk].nSW != 1 ||
                    Blk_List.Block[neighbour_blk].infoSW[0].cpu != i_cpu ||
                    Blk_List.Block[neighbour_blk].infoSW[0].blknum != i_blk) return(4104);
                for (l = 0; l <= buffer_size-1; ++l) {
	           Blk_List.message_noreschange_southwestcorner_recbuf[neighbour_blk][l] =
                     Blk_List.message_noreschange_southeastcorner_sendbuf[i_blk][l];
                } /* endfor */
             } else {
                if (Blk_List.Block[neighbour_blk].nSE != 1 ||
                    Blk_List.Block[neighbour_blk].infoSE[0].cpu != i_cpu ||
                    Blk_List.Block[neighbour_blk].infoSE[0].blknum != i_blk) return(4105);
                for (l = 0; l <= buffer_size-1; ++l) {
	           Blk_List.message_noreschange_southeastcorner_recbuf[neighbour_blk][l] =
                     Blk_List.message_noreschange_southeastcorner_sendbuf[i_blk][l];
                } /* endfor */
             } /* endif */
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
          buffer_size = sqr(Blk_List.Block[i_blk].infoSW[0].dimen.ghost)*
                        Number_of_Solution_Variables;
          if (neighbour_cpu != i_cpu) {
             if (Blk_List.Block[i_blk].infoSW[0].dimen.i > 0 &&
                 Blk_List.Block[i_blk].infoSW[0].dimen.j > 0) {
                tag_receive = neighbour_blk*tag_base + tag_SW;
                tag_send    = i_blk*tag_base + tag_NE;
             } else if (Blk_List.Block[i_blk].infoSW[0].dimen.j > 0) {
                tag_receive = neighbour_blk*tag_base + tag_SW;
                tag_send    = i_blk*tag_base + tag_NW;
             } else if (Blk_List.Block[i_blk].infoSW[0].dimen.i > 0) {
                tag_receive = neighbour_blk*tag_base + tag_SW;
                tag_send    = i_blk*tag_base + tag_SE;
             } else {
                tag_receive = neighbour_blk*tag_base + tag_SW;
                tag_send    = i_blk*tag_base + tag_SW;
             } /* endif */
            receive_requests[number_receive_requests] = MPI::COMM_WORLD.Irecv(Blk_List.message_noreschange_southwestcorner_recbuf[i_blk],
									      buffer_size,
									      MPI::DOUBLE,
									      neighbour_cpu,
									      tag_receive);
	    number_receive_requests++;
	    send_requests[number_send_requests] = MPI::COMM_WORLD.Isend(Blk_List.message_noreschange_southwestcorner_sendbuf[i_blk],
									buffer_size,
									MPI::DOUBLE,
									neighbour_cpu,
									tag_send);
	    number_send_requests++;

          } else {
	     if (!Blk_List.Block[neighbour_blk].used) return(5100); 
             buffer_size_neighbour = sqr(Blk_List.Block[neighbour_blk].info.dimen.ghost)*
                                     Number_of_Solution_Variables;
             if (buffer_size != buffer_size_neighbour) return(5101);
             if (Blk_List.Block[i_blk].infoSW[0].dimen.i > 0 &&
                 Blk_List.Block[i_blk].infoSW[0].dimen.j > 0) {
                if (Blk_List.Block[neighbour_blk].nNE != 1 ||
                    Blk_List.Block[neighbour_blk].infoNE[0].cpu != i_cpu ||
                    Blk_List.Block[neighbour_blk].infoNE[0].blknum != i_blk) return(5102);
                for (l = 0; l <= buffer_size-1; ++l) {
	           Blk_List.message_noreschange_northeastcorner_recbuf[neighbour_blk][l] =
                     Blk_List.message_noreschange_southwestcorner_sendbuf[i_blk][l];
                } /* endfor */
             } else if (Blk_List.Block[i_blk].infoSW[0].dimen.j > 0) {
                if (Blk_List.Block[neighbour_blk].nNW != 1 ||
                    Blk_List.Block[neighbour_blk].infoNW[0].cpu != i_cpu ||
                    Blk_List.Block[neighbour_blk].infoNW[0].blknum != i_blk) return(5103);
                for (l = 0; l <= buffer_size-1; ++l) {
	           Blk_List.message_noreschange_northwestcorner_recbuf[neighbour_blk][l] =
                     Blk_List.message_noreschange_southwestcorner_sendbuf[i_blk][l];
                } /* endfor */
             } else if (Blk_List.Block[i_blk].infoSW[0].dimen.i > 0) {
                if (Blk_List.Block[neighbour_blk].nSE != 1 ||
                    Blk_List.Block[neighbour_blk].infoSE[0].cpu != i_cpu ||
                    Blk_List.Block[neighbour_blk].infoSE[0].blknum != i_blk) return(5104);
                for (l = 0; l <= buffer_size-1; ++l) {
	           Blk_List.message_noreschange_southeastcorner_recbuf[neighbour_blk][l] =
                     Blk_List.message_noreschange_southwestcorner_sendbuf[i_blk][l];
                } /* endfor */
             } else {
                if (Blk_List.Block[neighbour_blk].nSW != 1 ||
                    Blk_List.Block[neighbour_blk].infoSW[0].cpu != i_cpu ||
                    Blk_List.Block[neighbour_blk].infoSW[0].blknum != i_blk) return(5105);
                for (l = 0; l <= buffer_size-1; ++l) {
	           Blk_List.message_noreschange_southwestcorner_recbuf[neighbour_blk][l] =
                     Blk_List.message_noreschange_southwestcorner_sendbuf[i_blk][l];
                } /* endfor */
             } /* endif */
          } /* endif */
       } /* endif */

    }  /* endfor */
  
    /* Explicitly Deallocate send buffers memory */
    for(int i=0; i < number_send_requests; i++){
      send_requests[i].Free(); 
    }

    /* Wait for all messages to be received. */
    /* Also deallocates recv buffers memory */ 
    if (number_receive_requests > 0) {
       MPI::Request::Waitall(number_receive_requests, receive_requests);
    } /* endif */

    /* Deallocate memory for message passing requests. */
    delete []receive_requests;
    receive_requests = NULL;
    delete []send_requests;
    send_requests = NULL;

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

    int number_receive_requests, 
        number_send_requests, 
        tag_receive, tag_send;
    int tag_N = 22, tag_S = 21, tag_E = 12, tag_W = 11, 
        tag_NW = 42, tag_NE = 52, tag_SE = 41, tag_SW = 51,
        tag_base = 100000;

    MPI::Request *receive_requests, 
                 *send_requests;
                 
    /* Allocate memory for message passing requests. */

    receive_requests = new MPI::Request[8*Blk_List.Nblk];
    send_requests = new MPI::Request[8*Blk_List.Nblk];
    number_receive_requests = 0;
    number_send_requests = 0;

    /* Perform message passing for the block faces and corners of each 
       more refined (higher mesh resolution) solution block having more 
       coarse (lower mesh resolution) neighbouring solution blocks. */

    i_cpu = Blk_List.ThisCPU;

    for ( i_blk = 0 ; i_blk < Blk_List.Nblk; ++i_blk ) {

       //////////////////////////////////////////////////////////////
       // Exchange messages at the North faces of solution blocks  //
       // with coarser neighbour at lower mesh resolution.         //
       //////////////////////////////////////////////////////////////

       // Post send messages to neighboring blocks at north faces.
       if (Blk_List.Block[i_blk].used && 
           (Blk_List.Block[i_blk].nN == 1) &&
           (Blk_List.Block[i_blk].info.level >
            Blk_List.Block[i_blk].infoN[0].level)) {
	  neighbour_cpu = Blk_List.Block[i_blk].infoN[0].cpu;
	  neighbour_blk = Blk_List.Block[i_blk].infoN[0].blknum;
          buffer_size = Blk_List.Block[i_blk].infoN[0].dimen.ghost*
                        (abs(Blk_List.Block[i_blk].infoN[0].dimen.i)/2+1)*
                        Number_of_Solution_Variables+
                        2*Blk_List.Block[i_blk].infoN[0].dimen.ghost;
          if (neighbour_cpu != i_cpu) {
             if (Blk_List.Block[i_blk].infoN[0].dimen.j > 0) {
                tag_send    = i_blk*tag_base + tag_S;
             } else {
                tag_send    = i_blk*tag_base + tag_N;
             } /* endif */
             send_requests[number_send_requests] = MPI::COMM_WORLD.Isend(Blk_List.message_reschange_northface_sendbuf[i_blk][0],
									 buffer_size,
									 MPI::DOUBLE,
									 neighbour_cpu,
									 tag_send);
             number_send_requests++;
          } else {
             if (!Blk_List.Block[neighbour_blk].used) return(2200);
             buffer_size_neighbour = Blk_List.Block[neighbour_blk].info.dimen.ghost*
                                     (abs(Blk_List.Block[neighbour_blk].info.dimen.i)/2+1)*
                                     Number_of_Solution_Variables+
                                     2*Blk_List.Block[neighbour_blk].info.dimen.ghost;
             if (buffer_size != buffer_size_neighbour) return(2201);

             if (Blk_List.Block[i_blk].infoN[0].dimen.j > 0) {
	        if (Blk_List.Block[neighbour_blk].nS == 2) {
                   if (Blk_List.Block[neighbour_blk].infoS[0].cpu == i_cpu &&
                       Blk_List.Block[neighbour_blk].infoS[0].blknum == i_blk) {
                      j_neigh = 0;
                   } else if (Blk_List.Block[neighbour_blk].infoS[1].cpu == i_cpu &&
                              Blk_List.Block[neighbour_blk].infoS[1].blknum == i_blk) {
                      j_neigh = 1;
                   } else {
                      return(2202);
                   } /* endif */
                } else {
                   return(2203);
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
                      return(2204);
                   } /* endif */
                } else {
                   return(2205);
                } /* endif */
                for (l = 0; l <= buffer_size-1; ++l) {
	           Blk_List.message_reschange_northface_recbuf[neighbour_blk][j_neigh][l] =
                     Blk_List.message_reschange_northface_sendbuf[i_blk][0][l];
                } /* endfor */
             } /* endif */

          } /* endif */
       } /* endif */

       // Post receive messages from neighboring blocks at north faces.
       if (Blk_List.Block[i_blk].used && 
           (Blk_List.Block[i_blk].nN == 2) &&
           (Blk_List.Block[i_blk].info.level <
            Blk_List.Block[i_blk].infoN[0].level)) {
          for (j_neigh = 0; j_neigh <= Blk_List.Block[i_blk].nN-1 ; ++j_neigh) {
             neighbour_cpu = Blk_List.Block[i_blk].infoN[j_neigh].cpu;
             neighbour_blk = Blk_List.Block[i_blk].infoN[j_neigh].blknum;
             buffer_size = Blk_List.Block[i_blk].infoN[j_neigh].dimen.ghost*
                           (abs(Blk_List.Block[i_blk].infoN[j_neigh].dimen.i)/2+1)*
                           Number_of_Solution_Variables+
                           2*Blk_List.Block[i_blk].infoN[j_neigh].dimen.ghost;
             if (neighbour_cpu != i_cpu) {
                tag_receive = neighbour_blk*tag_base + tag_N;
                receive_requests[number_receive_requests] = MPI::COMM_WORLD.Irecv(Blk_List.message_reschange_northface_recbuf[i_blk][j_neigh],
										  buffer_size,
										  MPI::DOUBLE,
										  neighbour_cpu,
										  tag_receive);
                number_receive_requests++;
             } /* endif */
          } /* endfor */
       } /* endif */

       //////////////////////////////////////////////////////////////
       // Exchange messages at the South faces of solution blocks  //
       // with coarser neighbour at lower mesh resolution.         //
       //////////////////////////////////////////////////////////////

       // Post send messages to neighboring blocks at south faces.
       if (Blk_List.Block[i_blk].used && 
           (Blk_List.Block[i_blk].nS == 1) &&
           (Blk_List.Block[i_blk].info.level >
            Blk_List.Block[i_blk].infoS[0].level)) {
	  neighbour_cpu = Blk_List.Block[i_blk].infoS[0].cpu;
	  neighbour_blk = Blk_List.Block[i_blk].infoS[0].blknum;
          buffer_size = Blk_List.Block[i_blk].infoS[0].dimen.ghost*
                        (abs(Blk_List.Block[i_blk].infoS[0].dimen.i)/2+1)*
                        Number_of_Solution_Variables+
                        2*Blk_List.Block[i_blk].infoS[0].dimen.ghost;
          if (neighbour_cpu != i_cpu) {
             if (Blk_List.Block[i_blk].infoS[0].dimen.j > 0) {
                tag_send    = i_blk*tag_base + tag_N;
             } else {
                tag_send    = i_blk*tag_base + tag_S;
             } /* endif */
	     send_requests[number_send_requests] = MPI::COMM_WORLD.Isend(Blk_List.message_reschange_southface_sendbuf[i_blk][0],
									 buffer_size,
									 MPI::DOUBLE,
									 neighbour_cpu,
									 tag_send);
             number_send_requests++;
          } else {
             if (!Blk_List.Block[neighbour_blk].used) return(2100);
             buffer_size_neighbour = Blk_List.Block[neighbour_blk].info.dimen.ghost*
                                     (abs(Blk_List.Block[neighbour_blk].info.dimen.i)/2+1)*
                                     Number_of_Solution_Variables+
                                     2*Blk_List.Block[neighbour_blk].info.dimen.ghost;
             if (buffer_size != buffer_size_neighbour) return(2101);

             if (Blk_List.Block[i_blk].infoS[0].dimen.j > 0) {
	        if (Blk_List.Block[neighbour_blk].nN == 2) {
                   if (Blk_List.Block[neighbour_blk].infoN[0].cpu == i_cpu &&
                       Blk_List.Block[neighbour_blk].infoN[0].blknum == i_blk) {
                      j_neigh = 0;
                   } else if (Blk_List.Block[neighbour_blk].infoN[1].cpu == i_cpu &&
                              Blk_List.Block[neighbour_blk].infoN[1].blknum == i_blk) {
                      j_neigh = 1;
                   } else {
                      return(2102);
                   } /* endif */
                } else {
                   return(2103);
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
                      return(2104);
                   } /* endif */
                } else {
                   return(2105);
                } /* endif */
                for (l = 0; l <= buffer_size-1; ++l) {
	           Blk_List.message_reschange_southface_recbuf[neighbour_blk][j_neigh][l] =
                     Blk_List.message_reschange_southface_sendbuf[i_blk][0][l];
                } /* endfor */
             } /* endif */

          } /* endif */
       } /* endif */

       // Post receive messages from neighboring blocks at south faces.
       if (Blk_List.Block[i_blk].used && 
           (Blk_List.Block[i_blk].nS == 2) &&
           (Blk_List.Block[i_blk].info.level <
            Blk_List.Block[i_blk].infoS[0].level)) {
          for (j_neigh = 0; j_neigh <= Blk_List.Block[i_blk].nS-1 ; ++j_neigh) {
             neighbour_cpu = Blk_List.Block[i_blk].infoS[j_neigh].cpu;
             neighbour_blk = Blk_List.Block[i_blk].infoS[j_neigh].blknum;
             buffer_size = Blk_List.Block[i_blk].infoS[j_neigh].dimen.ghost*
                           (abs(Blk_List.Block[i_blk].infoS[j_neigh].dimen.i)/2+1)*
                           Number_of_Solution_Variables+
                           2*Blk_List.Block[i_blk].infoS[j_neigh].dimen.ghost;
             if (neighbour_cpu != i_cpu) {
                tag_receive = neighbour_blk*tag_base + tag_S;
                receive_requests[number_receive_requests] = MPI::COMM_WORLD.Irecv(Blk_List.message_reschange_southface_recbuf[i_blk][j_neigh],
										    buffer_size,
										    MPI::DOUBLE,
										    neighbour_cpu,
										    tag_receive);
                number_receive_requests++;
             } /* endif */
          } /* endfor */
       } /* endif */

       //////////////////////////////////////////////////////////////
       // Exchange messages at the East faces of solution blocks   //
       // with coarser neighbour at lower mesh resolution.         //
       //////////////////////////////////////////////////////////////

       // Post send messages to neighboring blocks at east faces.
       if (Blk_List.Block[i_blk].used && 
           (Blk_List.Block[i_blk].nE == 1) &&
           (Blk_List.Block[i_blk].info.level >
            Blk_List.Block[i_blk].infoE[0].level)) {
	  neighbour_cpu = Blk_List.Block[i_blk].infoE[0].cpu;
	  neighbour_blk = Blk_List.Block[i_blk].infoE[0].blknum;
          buffer_size = Blk_List.Block[i_blk].infoE[0].dimen.ghost*
                        (abs(Blk_List.Block[i_blk].infoE[0].dimen.j)/2+1)*
                        Number_of_Solution_Variables+
                        2*Blk_List.Block[i_blk].infoE[0].dimen.ghost;
          if (neighbour_cpu != i_cpu) {
             if (Blk_List.Block[i_blk].infoE[0].dimen.i > 0) {
                tag_send    = i_blk*tag_base + tag_W;
             } else {
                tag_send    = i_blk*tag_base + tag_E;
             } /* endif */
             send_requests[number_send_requests] = MPI::COMM_WORLD.Isend(Blk_List.message_reschange_eastface_sendbuf[i_blk][0],
									 buffer_size,
									 MPI::DOUBLE,
									 neighbour_cpu,
									 tag_send);
             number_send_requests++;
          } else {
             if (!Blk_List.Block[neighbour_blk].used) return(1200);
             buffer_size_neighbour = Blk_List.Block[neighbour_blk].info.dimen.ghost*
                                     (abs(Blk_List.Block[neighbour_blk].info.dimen.j)/2+1)*
                                     Number_of_Solution_Variables+
                                     2*Blk_List.Block[neighbour_blk].info.dimen.ghost;
             if (buffer_size != buffer_size_neighbour) return(1201);
             if (Blk_List.Block[i_blk].infoE[0].dimen.i > 0) {
	        if (Blk_List.Block[neighbour_blk].nW == 2) {
                   if (Blk_List.Block[neighbour_blk].infoW[0].cpu == i_cpu &&
                       Blk_List.Block[neighbour_blk].infoW[0].blknum == i_blk) {
                      j_neigh = 0;
                   } else if (Blk_List.Block[neighbour_blk].infoW[1].cpu == i_cpu &&
                              Blk_List.Block[neighbour_blk].infoW[1].blknum == i_blk) {
                      j_neigh = 1;
                   } else {
                      return(1202);
                   } /* endif */
                } else {
                   return(1203);
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
                      return(1204);
                   } /* endif */
                } else {
                   return(1205);
                } /* endif */
                for (l = 0; l <= buffer_size-1; ++l) {
	           Blk_List.message_reschange_eastface_recbuf[neighbour_blk][j_neigh][l] =
                     Blk_List.message_reschange_eastface_sendbuf[i_blk][0][l];
                } /* endfor */
             } /* endif */

          } /* endif */
       } /* endif */

       // Post receive messages from neighboring blocks at east faces.
       if (Blk_List.Block[i_blk].used && 
           (Blk_List.Block[i_blk].nE == 2) &&
           (Blk_List.Block[i_blk].info.level <
            Blk_List.Block[i_blk].infoE[0].level)) {
          for (j_neigh = 0; j_neigh <= Blk_List.Block[i_blk].nE-1 ; ++j_neigh) {
             neighbour_cpu = Blk_List.Block[i_blk].infoE[j_neigh].cpu;
             neighbour_blk = Blk_List.Block[i_blk].infoE[j_neigh].blknum;
             buffer_size = Blk_List.Block[i_blk].infoE[j_neigh].dimen.ghost*
                           (abs(Blk_List.Block[i_blk].infoE[j_neigh].dimen.j)/2+1)*
                           Number_of_Solution_Variables+
                           2*Blk_List.Block[i_blk].infoE[j_neigh].dimen.ghost;
             if (neighbour_cpu != i_cpu) {
                tag_receive = neighbour_blk*tag_base + tag_E;
                receive_requests[number_receive_requests] = MPI::COMM_WORLD.Irecv(Blk_List.message_reschange_eastface_recbuf[i_blk][j_neigh],
										  buffer_size,
										  MPI::DOUBLE,
										  neighbour_cpu,
										  tag_receive);
                number_receive_requests++;
             } /* endif */
          } /* endfor */
       } /* endif */

       //////////////////////////////////////////////////////////////
       // Exchange messages at the West faces of solution blocks   //
       // with coarser neighbour at lower mesh resolution.         //
       //////////////////////////////////////////////////////////////

       // Post send messages to neighboring blocks at west faces.
       if (Blk_List.Block[i_blk].used && 
           (Blk_List.Block[i_blk].nW == 1) &&
           (Blk_List.Block[i_blk].info.level >
            Blk_List.Block[i_blk].infoW[0].level)) {
	  neighbour_cpu = Blk_List.Block[i_blk].infoW[0].cpu;
	  neighbour_blk = Blk_List.Block[i_blk].infoW[0].blknum;
          buffer_size = Blk_List.Block[i_blk].infoW[0].dimen.ghost*
                        (abs(Blk_List.Block[i_blk].infoW[0].dimen.j)/2+1)*
                        Number_of_Solution_Variables+
                        2*Blk_List.Block[i_blk].infoW[0].dimen.ghost;
          if (neighbour_cpu != i_cpu) {
             if (Blk_List.Block[i_blk].infoW[0].dimen.i > 0) {
                tag_send    = i_blk*tag_base + tag_E;
             } else {
                tag_send    = i_blk*tag_base + tag_W;
             } /* endif */
             send_requests[number_send_requests] = MPI::COMM_WORLD.Isend(Blk_List.message_reschange_westface_sendbuf[i_blk][0],
									 buffer_size,
									 MPI::DOUBLE,
									 neighbour_cpu,
									 tag_send);
             number_send_requests++;
          } else {
             if (!Blk_List.Block[neighbour_blk].used) return(1100); 
             buffer_size_neighbour = Blk_List.Block[neighbour_blk].info.dimen.ghost*
                                     (abs(Blk_List.Block[neighbour_blk].info.dimen.j)/2+1)*
                                     Number_of_Solution_Variables+
                                     2*Blk_List.Block[neighbour_blk].info.dimen.ghost;
             if (buffer_size != buffer_size_neighbour) return(1101);

             if (Blk_List.Block[i_blk].infoW[0].dimen.i > 0) {
	        if (Blk_List.Block[neighbour_blk].nE == 2) {
                   if (Blk_List.Block[neighbour_blk].infoE[0].cpu == i_cpu &&
                       Blk_List.Block[neighbour_blk].infoE[0].blknum == i_blk) {
                      j_neigh = 0;
                   } else if (Blk_List.Block[neighbour_blk].infoE[1].cpu == i_cpu &&
                              Blk_List.Block[neighbour_blk].infoE[1].blknum == i_blk) {
                      j_neigh = 1;
                   } else {
                      return(1102);
                   } /* endif */
                } else {
                   return(1103);
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
                      return(1104);
                   } /* endif */
                } else {
                   return(1105);
                } /* endif */
                for (l = 0; l <= buffer_size-1; ++l) {
	           Blk_List.message_reschange_westface_recbuf[neighbour_blk][j_neigh][l] =
                     Blk_List.message_reschange_westface_sendbuf[i_blk][0][l];
                } /* endfor */
             } /* endif */

          } /* endif */
       } /* endif */

       // Post receive messages from neighboring blocks at west faces.
       if (Blk_List.Block[i_blk].used && 
           (Blk_List.Block[i_blk].nW == 2) &&
           (Blk_List.Block[i_blk].info.level <
            Blk_List.Block[i_blk].infoW[0].level)) {
          for (j_neigh = 0; j_neigh <= Blk_List.Block[i_blk].nW-1 ; ++j_neigh) {
             neighbour_cpu = Blk_List.Block[i_blk].infoW[j_neigh].cpu;
             neighbour_blk = Blk_List.Block[i_blk].infoW[j_neigh].blknum;
             buffer_size = Blk_List.Block[i_blk].infoW[j_neigh].dimen.ghost*
                           (abs(Blk_List.Block[i_blk].infoW[j_neigh].dimen.j)/2+1)*
                           Number_of_Solution_Variables+
                           2*Blk_List.Block[i_blk].infoW[j_neigh].dimen.ghost;
             if (neighbour_cpu != i_cpu) {
                tag_receive = neighbour_blk*tag_base + tag_W;
                receive_requests[number_receive_requests] = MPI::COMM_WORLD.Irecv(Blk_List.message_reschange_westface_recbuf[i_blk][j_neigh],
										  buffer_size,
										  MPI::DOUBLE,
										  neighbour_cpu,
										  tag_receive);
                number_receive_requests++;
             } /* endif */
          } /* endfor */
       } /* endif */

       //////////////////////////////////////////////////////////////
       // Exchange messages at the North West corners of solution  //
       // blocks with coarser neighbour at lower mesh resolution.  //
       //////////////////////////////////////////////////////////////

       // Post send messages to neighboring blocks at north west corners.
       if (Blk_List.Block[i_blk].used && 
           (Blk_List.Block[i_blk].nNW == 1) &&
           (Blk_List.Block[i_blk].info.level >
            Blk_List.Block[i_blk].infoNW[0].level)) {
	  neighbour_cpu = Blk_List.Block[i_blk].infoNW[0].cpu;
	  neighbour_blk = Blk_List.Block[i_blk].infoNW[0].blknum;
          buffer_size = sqr(Blk_List.Block[i_blk].infoNW[0].dimen.ghost)*
                        Number_of_Solution_Variables;
          if (neighbour_cpu != i_cpu) {
             if (Blk_List.Block[i_blk].infoNW[0].dimen.i > 0 &&
                 Blk_List.Block[i_blk].infoNW[0].dimen.j > 0) {
                tag_send    = i_blk*tag_base + tag_SE;
             } else if (Blk_List.Block[i_blk].infoNW[0].dimen.j > 0) {
                tag_send    = i_blk*tag_base + tag_SW;
             } else if (Blk_List.Block[i_blk].infoNW[0].dimen.i > 0) {
                tag_send    = i_blk*tag_base + tag_NE;
             } else {
                tag_send    = i_blk*tag_base + tag_NW;
             } /* endif */
             send_requests[number_send_requests] = MPI::COMM_WORLD.Isend(Blk_List.message_reschange_northwestcorner_sendbuf[i_blk],
									 buffer_size,
									 MPI::DOUBLE,
									 neighbour_cpu,
									 tag_send);
             number_send_requests++;
          } else {
             if (!Blk_List.Block[neighbour_blk].used) return(4200);
             buffer_size_neighbour = sqr(Blk_List.Block[neighbour_blk].info.dimen.ghost)*
                                     Number_of_Solution_Variables;
             if (buffer_size != buffer_size_neighbour) return(4201);

             if (Blk_List.Block[i_blk].infoNW[0].dimen.i > 0 &&
                 Blk_List.Block[i_blk].infoNW[0].dimen.j > 0) {
                if (Blk_List.Block[neighbour_blk].nSE != 1 ||
                    Blk_List.Block[neighbour_blk].infoSE[0].cpu != i_cpu ||
                    Blk_List.Block[neighbour_blk].infoSE[0].blknum != i_blk) return(4202);
                for (l = 0; l <= buffer_size-1; ++l) {
	           Blk_List.message_reschange_southeastcorner_recbuf[neighbour_blk][l] =
                     Blk_List.message_reschange_northwestcorner_sendbuf[i_blk][l];
                } /* endfor */
             } else if (Blk_List.Block[i_blk].infoNW[0].dimen.j > 0) {
                if (Blk_List.Block[neighbour_blk].nSW != 1 ||
                    Blk_List.Block[neighbour_blk].infoSW[0].cpu != i_cpu ||
                    Blk_List.Block[neighbour_blk].infoSW[0].blknum != i_blk) return(4203);
                for (l = 0; l <= buffer_size-1; ++l) {
	           Blk_List.message_reschange_southwestcorner_recbuf[neighbour_blk][l] =
                     Blk_List.message_reschange_northwestcorner_sendbuf[i_blk][l];
                } /* endfor */
             } else if (Blk_List.Block[i_blk].infoNW[0].dimen.i > 0) {
                if (Blk_List.Block[neighbour_blk].nNE != 1 ||
                    Blk_List.Block[neighbour_blk].infoNE[0].cpu != i_cpu ||
                    Blk_List.Block[neighbour_blk].infoNE[0].blknum != i_blk) return(4204);
                for (l = 0; l <= buffer_size-1; ++l) {
	           Blk_List.message_reschange_northeastcorner_recbuf[neighbour_blk][l] =
                     Blk_List.message_reschange_northwestcorner_sendbuf[i_blk][l];
                } /* endfor */
             } else {
                if (Blk_List.Block[neighbour_blk].nNW != 1 ||
                    Blk_List.Block[neighbour_blk].infoNW[0].cpu != i_cpu ||
                    Blk_List.Block[neighbour_blk].infoNW[0].blknum != i_blk) return(4205);
                for (l = 0; l <= buffer_size-1; ++l) {
	           Blk_List.message_reschange_northwestcorner_recbuf[neighbour_blk][l] =
                     Blk_List.message_reschange_northwestcorner_sendbuf[i_blk][l];
                } /* endfor */
             } /* endif */

          } /* endif */
       } /* endif */

       // Post receive messages from neighboring blocks at north west corners.
       if (Blk_List.Block[i_blk].used && 
           (Blk_List.Block[i_blk].nNW == 1) &&
           (Blk_List.Block[i_blk].info.level <
            Blk_List.Block[i_blk].infoNW[0].level)) {
	  neighbour_cpu = Blk_List.Block[i_blk].infoNW[0].cpu;
	  neighbour_blk = Blk_List.Block[i_blk].infoNW[0].blknum;
          buffer_size = sqr(Blk_List.Block[i_blk].infoNW[0].dimen.ghost)*
                        Number_of_Solution_Variables;
          if (neighbour_cpu != i_cpu) {
             tag_receive = neighbour_blk*tag_base + tag_NW;
             receive_requests[number_receive_requests] = MPI::COMM_WORLD.Irecv(Blk_List.message_reschange_northwestcorner_recbuf[i_blk],
									       buffer_size,
									       MPI::DOUBLE,
									       neighbour_cpu,
									       tag_receive);
             number_receive_requests++;
          } /* endif */
       } /* endif */

       //////////////////////////////////////////////////////////////
       // Exchange messages at the North East corners of solution  //
       // blocks with coarser neighbour at lower mesh resolution.  //
       //////////////////////////////////////////////////////////////

       // Post send messages to neighboring blocks at north east corners.
       if (Blk_List.Block[i_blk].used && 
           (Blk_List.Block[i_blk].nNE == 1) &&
           (Blk_List.Block[i_blk].info.level >
            Blk_List.Block[i_blk].infoNE[0].level)) {
	  neighbour_cpu = Blk_List.Block[i_blk].infoNE[0].cpu;
	  neighbour_blk = Blk_List.Block[i_blk].infoNE[0].blknum;
          buffer_size = sqr(Blk_List.Block[i_blk].infoNE[0].dimen.ghost)*
                        Number_of_Solution_Variables;
          if (neighbour_cpu != i_cpu) {
             if (Blk_List.Block[i_blk].infoNE[0].dimen.i > 0 &&
                 Blk_List.Block[i_blk].infoNE[0].dimen.j > 0) {
                tag_send    = i_blk*tag_base + tag_SW;
             } else if (Blk_List.Block[i_blk].infoNE[0].dimen.j > 0) {
                tag_send    = i_blk*tag_base + tag_SE;
             } else if (Blk_List.Block[i_blk].infoNE[0].dimen.i > 0) {
                tag_send    = i_blk*tag_base + tag_NW;
             } else {
                tag_send    = i_blk*tag_base + tag_NE;
             } /* endif */
             send_requests[number_send_requests] = MPI::COMM_WORLD.Isend(Blk_List.message_reschange_northeastcorner_sendbuf[i_blk],
									 buffer_size,
									 MPI::DOUBLE,
									 neighbour_cpu,
									 tag_send);
             number_send_requests++;
          } else {
             if (!Blk_List.Block[neighbour_blk].used) return(5200);
             buffer_size_neighbour = sqr(Blk_List.Block[neighbour_blk].info.dimen.ghost)*
                                     Number_of_Solution_Variables;
             if (buffer_size != buffer_size_neighbour) return(5201);

             if (Blk_List.Block[i_blk].infoNE[0].dimen.i > 0 &&
                 Blk_List.Block[i_blk].infoNE[0].dimen.j > 0) {
                if (Blk_List.Block[neighbour_blk].nSW != 1 ||
                    Blk_List.Block[neighbour_blk].infoSW[0].cpu != i_cpu ||
                    Blk_List.Block[neighbour_blk].infoSW[0].blknum != i_blk) return(5202);
                for (l = 0; l <= buffer_size-1; ++l) {
	           Blk_List.message_reschange_southwestcorner_recbuf[neighbour_blk][l] =
                     Blk_List.message_reschange_northeastcorner_sendbuf[i_blk][l];
                } /* endfor */
             } else if (Blk_List.Block[i_blk].infoNE[0].dimen.j > 0) {
                if (Blk_List.Block[neighbour_blk].nSE != 1 ||
                    Blk_List.Block[neighbour_blk].infoSE[0].cpu != i_cpu ||
                    Blk_List.Block[neighbour_blk].infoSE[0].blknum != i_blk) return(5203);
                for (l = 0; l <= buffer_size-1; ++l) {
	           Blk_List.message_reschange_southeastcorner_recbuf[neighbour_blk][l] =
                     Blk_List.message_reschange_northeastcorner_sendbuf[i_blk][l];
                } /* endfor */
             } else if (Blk_List.Block[i_blk].infoNE[0].dimen.i > 0) {
                if (Blk_List.Block[neighbour_blk].nNW != 1 ||
                    Blk_List.Block[neighbour_blk].infoNW[0].cpu != i_cpu ||
                    Blk_List.Block[neighbour_blk].infoNW[0].blknum != i_blk) return(5204);
                for (l = 0; l <= buffer_size-1; ++l) {
	           Blk_List.message_reschange_northwestcorner_recbuf[neighbour_blk][l] =
                     Blk_List.message_reschange_northeastcorner_sendbuf[i_blk][l];
                } /* endfor */
             } else {
                if (Blk_List.Block[neighbour_blk].nNE != 1 ||
                    Blk_List.Block[neighbour_blk].infoNE[0].cpu != i_cpu ||
                    Blk_List.Block[neighbour_blk].infoNE[0].blknum != i_blk) return(5205);
                for (l = 0; l <= buffer_size-1; ++l) {
	           Blk_List.message_reschange_northeastcorner_recbuf[neighbour_blk][l] =
                     Blk_List.message_reschange_northeastcorner_sendbuf[i_blk][l];
                } /* endfor */
             } /* endif */

          } /* endif */
       } /* endif */

       // Post receive messages from neighboring blocks at north east corners.
       if (Blk_List.Block[i_blk].used && 
           (Blk_List.Block[i_blk].nNE == 1) &&
           (Blk_List.Block[i_blk].info.level <
            Blk_List.Block[i_blk].infoNE[0].level)) {
	  neighbour_cpu = Blk_List.Block[i_blk].infoNE[0].cpu;
	  neighbour_blk = Blk_List.Block[i_blk].infoNE[0].blknum;
          buffer_size = sqr(Blk_List.Block[i_blk].infoNE[0].dimen.ghost)*
                        Number_of_Solution_Variables;
          if (neighbour_cpu != i_cpu) {
             tag_receive = neighbour_blk*tag_base + tag_NE;
             receive_requests[number_receive_requests] = MPI::COMM_WORLD.Irecv(Blk_List.message_reschange_northeastcorner_recbuf[i_blk],
									       buffer_size,
										MPI::DOUBLE,
										neighbour_cpu,
										tag_receive);
             number_receive_requests++;
          } /* endif */
       } /* endif */


       //////////////////////////////////////////////////////////////
       // Exchange messages at the South East corners of solution  //
       // blocks with coarser neighbour at lower mesh resolution.  //
       //////////////////////////////////////////////////////////////

       // Post send messages to neighboring blocks at south east corners.
       if (Blk_List.Block[i_blk].used && 
           (Blk_List.Block[i_blk].nSE == 1) &&
           (Blk_List.Block[i_blk].info.level >
            Blk_List.Block[i_blk].infoSE[0].level)) {
	  neighbour_cpu = Blk_List.Block[i_blk].infoSE[0].cpu;
	  neighbour_blk = Blk_List.Block[i_blk].infoSE[0].blknum;
          buffer_size = sqr(Blk_List.Block[i_blk].infoSE[0].dimen.ghost)*
                        Number_of_Solution_Variables;
          if (neighbour_cpu != i_cpu) {
             if (Blk_List.Block[i_blk].infoSE[0].dimen.i > 0 &&
                 Blk_List.Block[i_blk].infoSE[0].dimen.j > 0) {
                tag_send    = i_blk*tag_base + tag_NW;
             } else if (Blk_List.Block[i_blk].infoSE[0].dimen.j > 0) {
                tag_send    = i_blk*tag_base + tag_NE;
             } else if (Blk_List.Block[i_blk].infoSE[0].dimen.i > 0) {
                tag_send    = i_blk*tag_base + tag_SW;
             } else {
                tag_send    = i_blk*tag_base + tag_SE;
             } /* endif */
             send_requests[number_send_requests] = MPI::COMM_WORLD.Isend(Blk_List.message_reschange_southeastcorner_sendbuf[i_blk],
									 buffer_size,
									 MPI::DOUBLE,
									 neighbour_cpu,
									 tag_send);
             number_send_requests++;
          } else {
             if (!Blk_List.Block[neighbour_blk].used) return(4100); 
             buffer_size_neighbour = sqr(Blk_List.Block[neighbour_blk].info.dimen.ghost)*
                                     Number_of_Solution_Variables;
             if (buffer_size != buffer_size_neighbour) return(4101);

             if (Blk_List.Block[i_blk].infoSE[0].dimen.i > 0 &&
                 Blk_List.Block[i_blk].infoSE[0].dimen.j > 0) {
                if (Blk_List.Block[neighbour_blk].nNW != 1 ||
                    Blk_List.Block[neighbour_blk].infoNW[0].cpu != i_cpu ||
                    Blk_List.Block[neighbour_blk].infoNW[0].blknum != i_blk) return(4102);
                for (l = 0; l <= buffer_size-1; ++l) {
	           Blk_List.message_reschange_northwestcorner_recbuf[neighbour_blk][l] =
                     Blk_List.message_reschange_southeastcorner_sendbuf[i_blk][l];
                } /* endfor */
             } else if (Blk_List.Block[i_blk].infoSE[0].dimen.j > 0) {
                if (Blk_List.Block[neighbour_blk].nNE != 1 ||
                    Blk_List.Block[neighbour_blk].infoNE[0].cpu != i_cpu ||
                    Blk_List.Block[neighbour_blk].infoNE[0].blknum != i_blk) return(4103);
                for (l = 0; l <= buffer_size-1; ++l) {
	           Blk_List.message_reschange_northeastcorner_recbuf[neighbour_blk][l] =
                     Blk_List.message_reschange_southeastcorner_sendbuf[i_blk][l];
                } /* endfor */
             } else if (Blk_List.Block[i_blk].infoSE[0].dimen.i > 0) {
                if (Blk_List.Block[neighbour_blk].nSW != 1 ||
                    Blk_List.Block[neighbour_blk].infoSW[0].cpu != i_cpu ||
                    Blk_List.Block[neighbour_blk].infoSW[0].blknum != i_blk) return(4104);
                for (l = 0; l <= buffer_size-1; ++l) {
	           Blk_List.message_reschange_southwestcorner_recbuf[neighbour_blk][l] =
                     Blk_List.message_reschange_southeastcorner_sendbuf[i_blk][l];
                } /* endfor */
             } else {
                if (Blk_List.Block[neighbour_blk].nSE != 1 ||
                    Blk_List.Block[neighbour_blk].infoSE[0].cpu != i_cpu ||
                    Blk_List.Block[neighbour_blk].infoSE[0].blknum != i_blk) return(4105);
                for (l = 0; l <= buffer_size-1; ++l) {
	           Blk_List.message_reschange_southeastcorner_recbuf[neighbour_blk][l] =
                     Blk_List.message_reschange_southeastcorner_sendbuf[i_blk][l];
                } /* endfor */
             } /* endif */

          } /* endif */
       } /* endif */

       // Post receive messages from neighboring blocks at south east corners.
       if (Blk_List.Block[i_blk].used && 
           (Blk_List.Block[i_blk].nSE == 1) &&
           (Blk_List.Block[i_blk].info.level <
            Blk_List.Block[i_blk].infoSE[0].level)) {
	  neighbour_cpu = Blk_List.Block[i_blk].infoSE[0].cpu;
	  neighbour_blk = Blk_List.Block[i_blk].infoSE[0].blknum;
          buffer_size = sqr(Blk_List.Block[i_blk].infoSE[0].dimen.ghost)*
                        Number_of_Solution_Variables;
          if (neighbour_cpu != i_cpu) {
             tag_receive = neighbour_blk*tag_base + tag_SE;
             receive_requests[number_receive_requests] = MPI::COMM_WORLD.Irecv(Blk_List.message_reschange_southeastcorner_recbuf[i_blk],
									       buffer_size,
									       MPI::DOUBLE,
									       neighbour_cpu,
									       tag_receive);
             number_receive_requests++;
          } /* endif */
      } /* endif */

       //////////////////////////////////////////////////////////////
       // Exchange messages at the South West corners of solution  //
       // blocks with coarser neighbour at lower mesh resolution.  //
       //////////////////////////////////////////////////////////////

       // Post send messages to neighboring blocks at south west corners.
       if (Blk_List.Block[i_blk].used && 
           (Blk_List.Block[i_blk].nSW == 1) &&
           (Blk_List.Block[i_blk].info.level >
            Blk_List.Block[i_blk].infoSW[0].level)) {
	  neighbour_cpu = Blk_List.Block[i_blk].infoSW[0].cpu;
	  neighbour_blk = Blk_List.Block[i_blk].infoSW[0].blknum;
          buffer_size = sqr(Blk_List.Block[i_blk].infoSW[0].dimen.ghost)*
                        Number_of_Solution_Variables;
          if (neighbour_cpu != i_cpu) {
             if (Blk_List.Block[i_blk].infoSW[0].dimen.i > 0 &&
                 Blk_List.Block[i_blk].infoSW[0].dimen.j > 0) {
                tag_send    = i_blk*tag_base + tag_NE;
             } else if (Blk_List.Block[i_blk].infoSW[0].dimen.j > 0) {
                tag_send    = i_blk*tag_base + tag_NW;
             } else if (Blk_List.Block[i_blk].infoSW[0].dimen.i > 0) {
                tag_send    = i_blk*tag_base + tag_SE;
             } else {
                tag_send    = i_blk*tag_base + tag_SW;
             } /* endif */
             send_requests[number_send_requests] = MPI::COMM_WORLD.Isend(Blk_List.message_reschange_southwestcorner_sendbuf[i_blk],
									  buffer_size,
									  MPI::DOUBLE,
									  neighbour_cpu,
									  tag_send);
             number_send_requests++;
          } else {
             if (!Blk_List.Block[neighbour_blk].used) return(5100);
             buffer_size_neighbour = sqr(Blk_List.Block[neighbour_blk].info.dimen.ghost)*
                                     Number_of_Solution_Variables;
             if (buffer_size != buffer_size_neighbour) return(5101);

             if (Blk_List.Block[i_blk].infoSW[0].dimen.i > 0 &&
                 Blk_List.Block[i_blk].infoSW[0].dimen.j > 0) {
                if (Blk_List.Block[neighbour_blk].nNE != 1 ||
                    Blk_List.Block[neighbour_blk].infoNE[0].cpu != i_cpu ||
                    Blk_List.Block[neighbour_blk].infoNE[0].blknum != i_blk) return(5102);
                for (l = 0; l <= buffer_size-1; ++l) {
	           Blk_List.message_reschange_northeastcorner_recbuf[neighbour_blk][l] =
                     Blk_List.message_reschange_southwestcorner_sendbuf[i_blk][l];
                } /* endfor */
             } else if (Blk_List.Block[i_blk].infoSW[0].dimen.j > 0) {
                if (Blk_List.Block[neighbour_blk].nNW != 1 ||
                    Blk_List.Block[neighbour_blk].infoNW[0].cpu != i_cpu ||
                    Blk_List.Block[neighbour_blk].infoNW[0].blknum != i_blk) return(5103);
                for (l = 0; l <= buffer_size-1; ++l) {
	           Blk_List.message_reschange_northwestcorner_recbuf[neighbour_blk][l] =
                     Blk_List.message_reschange_southwestcorner_sendbuf[i_blk][l];
                } /* endfor */
             } else if (Blk_List.Block[i_blk].infoSW[0].dimen.i > 0) {
                if (Blk_List.Block[neighbour_blk].nSE != 1 ||
                    Blk_List.Block[neighbour_blk].infoSE[0].cpu != i_cpu ||
                    Blk_List.Block[neighbour_blk].infoSE[0].blknum != i_blk) return(5104);
                for (l = 0; l <= buffer_size-1; ++l) {
	           Blk_List.message_reschange_southeastcorner_recbuf[neighbour_blk][l] =
                     Blk_List.message_reschange_southwestcorner_sendbuf[i_blk][l];
                } /* endfor */
             } else {
                if (Blk_List.Block[neighbour_blk].nSW != 1 ||
                    Blk_List.Block[neighbour_blk].infoSW[0].cpu != i_cpu ||
                    Blk_List.Block[neighbour_blk].infoSW[0].blknum != i_blk) return(5105);
                for (l = 0; l <= buffer_size-1; ++l) {
	           Blk_List.message_reschange_southwestcorner_recbuf[neighbour_blk][l] =
                     Blk_List.message_reschange_southwestcorner_sendbuf[i_blk][l];
                } /* endfor */
             } /* endif */

          } /* endif */
       } /* endif */

       // Post receive messages from neighboring blocks at south west corners.
       if (Blk_List.Block[i_blk].used && 
           (Blk_List.Block[i_blk].nSW == 1) &&
           (Blk_List.Block[i_blk].info.level <
            Blk_List.Block[i_blk].infoSW[0].level)) {
	  neighbour_cpu = Blk_List.Block[i_blk].infoSW[0].cpu;
	  neighbour_blk = Blk_List.Block[i_blk].infoSW[0].blknum;
          buffer_size = sqr(Blk_List.Block[i_blk].infoSW[0].dimen.ghost)*
                        Number_of_Solution_Variables;
          if (neighbour_cpu != i_cpu) {
             tag_receive = neighbour_blk*tag_base + tag_SW;
             receive_requests[number_receive_requests] = MPI::COMM_WORLD.Irecv(Blk_List.message_reschange_southwestcorner_recbuf[i_blk],
									       buffer_size,
									       MPI::DOUBLE,
									       neighbour_cpu,
									       tag_receive);
             number_receive_requests++;
          } /* endif */
       } /* endif */
 
    }  /* endfor */
   
    /* Explicitly Deallocate send buffers memory */
    for(int i=0; i < number_send_requests; i++){
      send_requests[i].Free(); 
    }

    /* Wait for all messages to be received. */
    /* Also deallocates recv buffers memory */ 
    if (number_receive_requests > 0) {
       MPI::Request::Waitall(number_receive_requests, receive_requests);
    } /* endif */

    /* Deallocate memory for message passing requests. */

    delete []receive_requests;
    receive_requests = NULL;
    delete []send_requests;
    send_requests = NULL;

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

    int number_receive_requests, 
        number_send_requests, 
        tag_receive, tag_send;
    int tag_N = 22, tag_S = 21, tag_E = 12, tag_W = 11, 
        tag_NW = 42, tag_NE = 52, tag_SE = 41, tag_SW = 51,
        tag_base = 1000000, tag_modifier = 100;

    MPI::Request *receive_requests, 
                 *send_requests;
                 
    /* Allocate memory for message passing requests. */

    receive_requests = new MPI::Request[8*Blk_List.Nblk];
    send_requests = new MPI::Request[8*Blk_List.Nblk];
    number_receive_requests = 0;
    number_send_requests = 0;

    /* Perform message passing for the block faces and corners of each 
       more coarse (lower mesh resolution) solution block having more 
       refined (higher mesh resolution) neighbouring solution blocks. */

    i_cpu = Blk_List.ThisCPU;

    for ( i_blk = 0 ; i_blk < Blk_List.Nblk ; ++i_blk ) {

       //////////////////////////////////////////////////////////////
       // Exchange messages at the North faces of solution blocks  //
       // with finer neighbours at higher mesh resolution.         //
       //////////////////////////////////////////////////////////////

       // Post send messages to neighboring blocks at north faces.
       if (Blk_List.Block[i_blk].used && 
           (Blk_List.Block[i_blk].nN == 2) &&
           (Blk_List.Block[i_blk].info.level <
            Blk_List.Block[i_blk].infoN[0].level)) {
          for (j_neigh = 0; j_neigh <= Blk_List.Block[i_blk].nN-1 ; ++j_neigh) {
             neighbour_cpu = Blk_List.Block[i_blk].infoN[j_neigh].cpu;
             neighbour_blk = Blk_List.Block[i_blk].infoN[j_neigh].blknum;
             buffer_size = Blk_List.Block[i_blk].infoN[j_neigh].dimen.ghost*
                           (abs(Blk_List.Block[i_blk].infoN[j_neigh].dimen.i)+
                           Blk_List.Block[i_blk].infoN[j_neigh].dimen.ghost+1)*
                           Number_of_Solution_Variables+
                           2*Blk_List.Block[i_blk].infoN[j_neigh].dimen.ghost;
             if (neighbour_cpu != i_cpu) {
                if (Blk_List.Block[i_blk].infoN[j_neigh].dimen.j > 0) {
                   tag_send    = i_blk*tag_base + neighbour_blk*tag_modifier + tag_S;
                } else {
                   tag_send    = i_blk*tag_base + neighbour_blk*tag_modifier + tag_N;
                } /* endif */
                send_requests[number_send_requests] = MPI::COMM_WORLD.Isend(Blk_List.message_reschange_northface_sendbuf[i_blk][j_neigh],
									    buffer_size,
									    MPI::DOUBLE,
									    neighbour_cpu,
									    tag_send);
                number_send_requests++;
             } else {
                if (!Blk_List.Block[neighbour_blk].used) return(1);
                buffer_size_neighbour = Blk_List.Block[neighbour_blk].info.dimen.ghost*
                                        (abs(Blk_List.Block[neighbour_blk].info.dimen.i)+
                                        Blk_List.Block[neighbour_blk].info.dimen.ghost+1)*
 	                                Number_of_Solution_Variables+
                                        2*Blk_List.Block[neighbour_blk].info.dimen.ghost;
                if (buffer_size != buffer_size_neighbour) return(1);

                if (Blk_List.Block[i_blk].infoN[j_neigh].dimen.j > 0) {
                   if (Blk_List.Block[neighbour_blk].nS != 1 ||
                       Blk_List.Block[neighbour_blk].infoS[0].cpu != i_cpu ||
                       Blk_List.Block[neighbour_blk].infoS[0].blknum != i_blk) return(1);
                   for (l = 0; l <= buffer_size-1; ++l) {
	              Blk_List.message_reschange_southface_recbuf[neighbour_blk][0][l] =
                        Blk_List.message_reschange_northface_sendbuf[i_blk][j_neigh][l];
                   } /* endfor */
                } else {
                   if (Blk_List.Block[neighbour_blk].nN != 1 ||
                       Blk_List.Block[neighbour_blk].infoN[0].cpu != i_cpu ||
                       Blk_List.Block[neighbour_blk].infoN[0].blknum != i_blk) return(1);
                   for (l = 0; l <= buffer_size-1; ++l) {
	              Blk_List.message_reschange_northface_recbuf[neighbour_blk][0][l] =
                        Blk_List.message_reschange_northface_sendbuf[i_blk][j_neigh][l];
                   } /* endfor */
                } /* endif */

             } /* endif */
          } /* endfor */
       } /* endif */

       // Post receive messages from neighboring blocks at north faces.
       if (Blk_List.Block[i_blk].used && 
           (Blk_List.Block[i_blk].nN == 1) &&
           (Blk_List.Block[i_blk].info.level >
            Blk_List.Block[i_blk].infoN[0].level)) {
	  neighbour_cpu = Blk_List.Block[i_blk].infoN[0].cpu;
	  neighbour_blk = Blk_List.Block[i_blk].infoN[0].blknum;
          buffer_size = Blk_List.Block[i_blk].infoN[0].dimen.ghost*
                        (abs(Blk_List.Block[i_blk].infoN[0].dimen.i)+
                        Blk_List.Block[i_blk].infoN[0].dimen.ghost+1)*
                        Number_of_Solution_Variables+
                        2*Blk_List.Block[i_blk].infoN[0].dimen.ghost;
          if (neighbour_cpu != i_cpu) {
             tag_receive = neighbour_blk*tag_base + i_blk*tag_modifier + tag_N;
             receive_requests[number_receive_requests] = MPI::COMM_WORLD.Irecv(Blk_List.message_reschange_northface_recbuf[i_blk][0],
									       buffer_size,
									       MPI::DOUBLE,
									       neighbour_cpu,
									       tag_receive);
             number_receive_requests++;
          } /* endif */
       } /* endif */

       //////////////////////////////////////////////////////////////
       // Exchange messages at the South faces of solution blocks  //
       // with finer neighbours at higher mesh resolution.         //
       //////////////////////////////////////////////////////////////

       // Post send messages to neighboring blocks at south faces.
       if (Blk_List.Block[i_blk].used && 
           (Blk_List.Block[i_blk].nS == 2) &&
           (Blk_List.Block[i_blk].info.level <
            Blk_List.Block[i_blk].infoS[0].level)) {
          for (j_neigh = 0; j_neigh <= Blk_List.Block[i_blk].nS-1 ; ++j_neigh) {
             neighbour_cpu = Blk_List.Block[i_blk].infoS[j_neigh].cpu;
             neighbour_blk = Blk_List.Block[i_blk].infoS[j_neigh].blknum;
             buffer_size = Blk_List.Block[i_blk].infoS[j_neigh].dimen.ghost*
                           (abs(Blk_List.Block[i_blk].infoS[j_neigh].dimen.i)+
                           Blk_List.Block[i_blk].infoS[j_neigh].dimen.ghost+1)*
                           Number_of_Solution_Variables+
                           2*Blk_List.Block[i_blk].infoS[j_neigh].dimen.ghost;
             if (neighbour_cpu != i_cpu) {
                if (Blk_List.Block[i_blk].infoS[j_neigh].dimen.j > 0) {
                   tag_send    = i_blk*tag_base + neighbour_blk*tag_modifier + tag_N;
                } else {
                   tag_send    = i_blk*tag_base + neighbour_blk*tag_modifier + tag_S;
                } /* endif */
                send_requests[number_send_requests] = MPI::COMM_WORLD.Isend(Blk_List.message_reschange_southface_sendbuf[i_blk][j_neigh],
									    buffer_size,
									    MPI::DOUBLE,
									    neighbour_cpu,
									    tag_send);
                number_send_requests++;
             } else {
                if (!Blk_List.Block[neighbour_blk].used) return(1);
                buffer_size_neighbour = Blk_List.Block[neighbour_blk].info.dimen.ghost*
                                        (abs(Blk_List.Block[neighbour_blk].info.dimen.i)+
                                        Blk_List.Block[neighbour_blk].info.dimen.ghost+1)*
 	                                Number_of_Solution_Variables+
                                        2*Blk_List.Block[neighbour_blk].info.dimen.ghost;
                if (buffer_size != buffer_size_neighbour) return(1);

                if (Blk_List.Block[i_blk].infoS[j_neigh].dimen.j > 0) {
                   if (Blk_List.Block[neighbour_blk].nN != 1 ||
                       Blk_List.Block[neighbour_blk].infoN[0].cpu != i_cpu ||
                       Blk_List.Block[neighbour_blk].infoN[0].blknum != i_blk) return(1);
                   for (l = 0; l <= buffer_size-1; ++l) {
	              Blk_List.message_reschange_northface_recbuf[neighbour_blk][0][l] =
                        Blk_List.message_reschange_southface_sendbuf[i_blk][j_neigh][l];
                   } /* endfor */
                } else {
                   if (Blk_List.Block[neighbour_blk].nS != 1 ||
                       Blk_List.Block[neighbour_blk].infoS[0].cpu != i_cpu ||
                       Blk_List.Block[neighbour_blk].infoS[0].blknum != i_blk) return(1);
                   for (l = 0; l <= buffer_size-1; ++l) {
	              Blk_List.message_reschange_southface_recbuf[neighbour_blk][0][l] =
                        Blk_List.message_reschange_southface_sendbuf[i_blk][j_neigh][l];
                   } /* endfor */
                } /* endif */

             } /* endif */
          } /* endfor */
       } /* endif */

       // Post receive messages from neighboring blocks at south faces.
       if (Blk_List.Block[i_blk].used && 
           (Blk_List.Block[i_blk].nS == 1) &&
           (Blk_List.Block[i_blk].info.level >
            Blk_List.Block[i_blk].infoS[0].level)) {
	  neighbour_cpu = Blk_List.Block[i_blk].infoS[0].cpu;
	  neighbour_blk = Blk_List.Block[i_blk].infoS[0].blknum;
          buffer_size = Blk_List.Block[i_blk].infoS[0].dimen.ghost*
                        (abs(Blk_List.Block[i_blk].infoS[0].dimen.i)+
                        Blk_List.Block[i_blk].infoS[0].dimen.ghost+1)*
                        Number_of_Solution_Variables+
                        2*Blk_List.Block[i_blk].infoS[0].dimen.ghost;
          if (neighbour_cpu != i_cpu) {
             tag_receive = neighbour_blk*tag_base + i_blk*tag_modifier + tag_S;
             receive_requests[number_receive_requests] = MPI::COMM_WORLD.Irecv(Blk_List.message_reschange_southface_recbuf[i_blk][0],
									       buffer_size,
									       MPI::DOUBLE,
									       neighbour_cpu,
									       tag_receive);
             number_receive_requests++;
          } /* endif */
       } /* endif */


       //////////////////////////////////////////////////////////////
       // Exchange messages at the East faces of solution blocks   //
       // with finer neighbours at higher mesh resolution.         //
       //////////////////////////////////////////////////////////////

       // Post send messages to neighboring blocks at east faces.
       if (Blk_List.Block[i_blk].used && 
           (Blk_List.Block[i_blk].nE == 2) &&
           (Blk_List.Block[i_blk].info.level <
            Blk_List.Block[i_blk].infoE[0].level)) {
          for (j_neigh = 0; j_neigh <= Blk_List.Block[i_blk].nE-1 ; ++j_neigh) {
             neighbour_cpu = Blk_List.Block[i_blk].infoE[j_neigh].cpu;
             neighbour_blk = Blk_List.Block[i_blk].infoE[j_neigh].blknum;
             buffer_size = Blk_List.Block[i_blk].infoE[j_neigh].dimen.ghost*
                           (abs(Blk_List.Block[i_blk].infoE[j_neigh].dimen.j)+
                           Blk_List.Block[i_blk].infoE[j_neigh].dimen.ghost+1)*
                           Number_of_Solution_Variables+
                           2*Blk_List.Block[i_blk].infoE[j_neigh].dimen.ghost;
             if (neighbour_cpu != i_cpu) {
                if (Blk_List.Block[i_blk].infoE[j_neigh].dimen.i > 0) {
                   tag_send    = i_blk*tag_base + neighbour_blk*tag_modifier + tag_W;
                } else {
                   tag_send    = i_blk*tag_base + neighbour_blk*tag_modifier + tag_E;
                } /* endif */
                send_requests[number_send_requests] = MPI::COMM_WORLD.Isend(Blk_List.message_reschange_eastface_sendbuf[i_blk][j_neigh],
									    buffer_size,
									    MPI::DOUBLE,
									    neighbour_cpu,
									    tag_send);
                number_send_requests++;
             } else {
                if (!Blk_List.Block[neighbour_blk].used) return(1);
                buffer_size_neighbour = Blk_List.Block[neighbour_blk].info.dimen.ghost*
                                        (abs(Blk_List.Block[neighbour_blk].info.dimen.j)+
                                        Blk_List.Block[neighbour_blk].info.dimen.ghost+1)*
 	                                Number_of_Solution_Variables+
                                        2*Blk_List.Block[neighbour_blk].info.dimen.ghost;
                if (buffer_size != buffer_size_neighbour) return(1);

                if (Blk_List.Block[i_blk].infoE[j_neigh].dimen.j > 0) {
                   if (Blk_List.Block[neighbour_blk].nW != 1 ||
                       Blk_List.Block[neighbour_blk].infoW[0].cpu != i_cpu ||
                       Blk_List.Block[neighbour_blk].infoW[0].blknum != i_blk) return(1);
                   for (l = 0; l <= buffer_size-1; ++l) {
	              Blk_List.message_reschange_westface_recbuf[neighbour_blk][0][l] =
                        Blk_List.message_reschange_eastface_sendbuf[i_blk][j_neigh][l];
                   } /* endfor */
                } else {
                   if (Blk_List.Block[neighbour_blk].nE != 1 ||
                       Blk_List.Block[neighbour_blk].infoE[0].cpu != i_cpu ||
                       Blk_List.Block[neighbour_blk].infoE[0].blknum != i_blk) return(1);
                   for (l = 0; l <= buffer_size-1; ++l) {
	              Blk_List.message_reschange_eastface_recbuf[neighbour_blk][0][l] =
                        Blk_List.message_reschange_eastface_sendbuf[i_blk][j_neigh][l];
                   } /* endfor */
                } /* endif */

             } /* endif */
          } /* endfor */
       } /* endif */

       // Post receive messages from neighboring blocks at east faces.
       if (Blk_List.Block[i_blk].used && 
           (Blk_List.Block[i_blk].nE == 1) &&
           (Blk_List.Block[i_blk].info.level >
            Blk_List.Block[i_blk].infoE[0].level)) {
	  neighbour_cpu = Blk_List.Block[i_blk].infoE[0].cpu;
	  neighbour_blk = Blk_List.Block[i_blk].infoE[0].blknum;
          buffer_size = Blk_List.Block[i_blk].infoE[0].dimen.ghost*
                        (abs(Blk_List.Block[i_blk].infoE[0].dimen.j)+
                        Blk_List.Block[i_blk].infoE[0].dimen.ghost+1)*
                        Number_of_Solution_Variables+
                        2*Blk_List.Block[i_blk].infoE[0].dimen.ghost;
          if (neighbour_cpu != i_cpu) {
             tag_receive = neighbour_blk*tag_base + i_blk*tag_modifier + tag_E;
	     receive_requests[number_receive_requests] = MPI::COMM_WORLD.Irecv(Blk_List.message_reschange_eastface_recbuf[i_blk][0],
									       buffer_size,
									       MPI::DOUBLE,
									       neighbour_cpu,
									       tag_receive);
	     number_receive_requests++;
          } /* endif */
       } /* endif */

       //////////////////////////////////////////////////////////////
       // Exchange messages at the West faces of solution blocks   //
       // with finer neighbours at higher mesh resolution.         //
       //////////////////////////////////////////////////////////////

       // Post send messages to neighboring blocks at west faces.
       if (Blk_List.Block[i_blk].used && 
           (Blk_List.Block[i_blk].nW == 2) &&
           (Blk_List.Block[i_blk].info.level <
            Blk_List.Block[i_blk].infoW[0].level)) {
          for (j_neigh = 0; j_neigh <= Blk_List.Block[i_blk].nW-1 ; ++j_neigh) {
             neighbour_cpu = Blk_List.Block[i_blk].infoW[j_neigh].cpu;
             neighbour_blk = Blk_List.Block[i_blk].infoW[j_neigh].blknum;
             buffer_size = Blk_List.Block[i_blk].infoW[j_neigh].dimen.ghost*
                           (abs(Blk_List.Block[i_blk].infoW[j_neigh].dimen.j)+
                           Blk_List.Block[i_blk].infoW[j_neigh].dimen.ghost+1)*
                           Number_of_Solution_Variables+
                           2*Blk_List.Block[i_blk].infoW[j_neigh].dimen.ghost;
             if (neighbour_cpu != i_cpu) {
                if (Blk_List.Block[i_blk].infoW[j_neigh].dimen.i > 0) {
                   tag_send    = i_blk*tag_base + neighbour_blk*tag_modifier + tag_E;
                } else {
                   tag_send    = i_blk*tag_base + neighbour_blk*tag_modifier + tag_W;
                } /* endif */
                send_requests[number_send_requests] = MPI::COMM_WORLD.Isend(Blk_List.message_reschange_westface_sendbuf[i_blk][j_neigh],
									     buffer_size,
									     MPI::DOUBLE,
									     neighbour_cpu,
									     tag_send);
                number_send_requests++;
             } else {
                if (!Blk_List.Block[neighbour_blk].used) return(1);
                buffer_size_neighbour = Blk_List.Block[neighbour_blk].info.dimen.ghost*
                                        (abs(Blk_List.Block[neighbour_blk].info.dimen.j)+
                                        Blk_List.Block[neighbour_blk].info.dimen.ghost+1)*
 	                                Number_of_Solution_Variables+
                                        2*Blk_List.Block[neighbour_blk].info.dimen.ghost;
                if (buffer_size != buffer_size_neighbour) return(1);
                if (Blk_List.Block[i_blk].infoW[j_neigh].dimen.j > 0) {
                   if (Blk_List.Block[neighbour_blk].nE != 1 ||
                       Blk_List.Block[neighbour_blk].infoE[0].cpu != i_cpu ||
                       Blk_List.Block[neighbour_blk].infoE[0].blknum != i_blk) return(1);
                   for (l = 0; l <= buffer_size-1; ++l) {
	              Blk_List.message_reschange_eastface_recbuf[neighbour_blk][0][l] =
                        Blk_List.message_reschange_westface_sendbuf[i_blk][j_neigh][l];
                   } /* endfor */
                } else {
                   if (Blk_List.Block[neighbour_blk].nW != 1 ||
                       Blk_List.Block[neighbour_blk].infoW[0].cpu != i_cpu ||
                       Blk_List.Block[neighbour_blk].infoW[0].blknum != i_blk) return(1);
                   for (l = 0; l <= buffer_size-1; ++l) {
	              Blk_List.message_reschange_westface_recbuf[neighbour_blk][0][l] =
                        Blk_List.message_reschange_westface_sendbuf[i_blk][j_neigh][l];
                   } /* endfor */
                } /* endif */

             } /* endif */
          } /* endfor */
       } /* endif */

       // Post receive messages from neighboring blocks at west faces.
       if (Blk_List.Block[i_blk].used && 
           (Blk_List.Block[i_blk].nW == 1) &&
           (Blk_List.Block[i_blk].info.level >
            Blk_List.Block[i_blk].infoW[0].level)) {
	  neighbour_cpu = Blk_List.Block[i_blk].infoW[0].cpu;
	  neighbour_blk = Blk_List.Block[i_blk].infoW[0].blknum;
          buffer_size = Blk_List.Block[i_blk].infoW[0].dimen.ghost*
                        (abs(Blk_List.Block[i_blk].infoW[0].dimen.j)+
                        Blk_List.Block[i_blk].infoW[0].dimen.ghost+1)*
                        Number_of_Solution_Variables+
                        2*Blk_List.Block[i_blk].infoW[0].dimen.ghost;
          if (neighbour_cpu != i_cpu) {
             tag_receive = neighbour_blk*tag_base + i_blk*tag_modifier + tag_W;
             receive_requests[number_receive_requests] = MPI::COMM_WORLD.Irecv(Blk_List.message_reschange_westface_recbuf[i_blk][0],
									       buffer_size,
									       MPI::DOUBLE,
									       neighbour_cpu,
									       tag_receive);
             number_receive_requests++;
          } /* endif */
       } /* endif */

       /////////////////////////////////////////////////////////////////////
       // Exchange messages at the North West corners of solution blocks  //
       // with finer neighbours at higher mesh resolution.                //
       /////////////////////////////////////////////////////////////////////

       // Post send messages to neighboring blocks at north west corners.
       if (Blk_List.Block[i_blk].used && 
           (Blk_List.Block[i_blk].nNW == 1) &&
           (Blk_List.Block[i_blk].info.level <
            Blk_List.Block[i_blk].infoNW[0].level)) {
	  neighbour_cpu = Blk_List.Block[i_blk].infoNW[0].cpu;
	  neighbour_blk = Blk_List.Block[i_blk].infoNW[0].blknum;
          buffer_size = sqr(Blk_List.Block[i_blk].infoNW[0].dimen.ghost)*
                        Number_of_Solution_Variables;
          if (neighbour_cpu != i_cpu) {
             if (Blk_List.Block[i_blk].infoNW[0].dimen.i > 0 &&
                 Blk_List.Block[i_blk].infoNW[0].dimen.j > 0) {
                tag_send    = i_blk*tag_base + tag_SE;
             } else if (Blk_List.Block[i_blk].infoNW[0].dimen.j > 0) {
                tag_send    = i_blk*tag_base + tag_SW;
             } else if (Blk_List.Block[i_blk].infoNW[0].dimen.i > 0) {
                tag_send    = i_blk*tag_base + tag_NE;
             } else {
                tag_send    = i_blk*tag_base + tag_NW;
             } /* endif */
             send_requests[number_send_requests] = MPI::COMM_WORLD.Isend(Blk_List.message_reschange_northwestcorner_sendbuf[i_blk],
									  buffer_size,
									  MPI::DOUBLE,
									  neighbour_cpu,
									  tag_send);
             number_send_requests++;
          } else {
             if (!Blk_List.Block[neighbour_blk].used) return(1);
             buffer_size_neighbour = sqr(Blk_List.Block[neighbour_blk].info.dimen.ghost)*
                                     Number_of_Solution_Variables;
             if (buffer_size != buffer_size_neighbour) return(1);

             if (Blk_List.Block[i_blk].infoNW[0].dimen.i > 0 &&
                 Blk_List.Block[i_blk].infoNW[0].dimen.j > 0) {
                if (Blk_List.Block[neighbour_blk].nSE != 1 ||
                    Blk_List.Block[neighbour_blk].infoSE[0].cpu != i_cpu ||
                    Blk_List.Block[neighbour_blk].infoSE[0].blknum != i_blk) return(1);
                for (l = 0; l <= buffer_size-1; ++l) {
	           Blk_List.message_reschange_southeastcorner_recbuf[neighbour_blk][l] =
                     Blk_List.message_reschange_northwestcorner_sendbuf[i_blk][l];
                } /* endfor */
             } else if (Blk_List.Block[i_blk].infoNW[0].dimen.j > 0) {
                if (Blk_List.Block[neighbour_blk].nSW != 1 ||
                    Blk_List.Block[neighbour_blk].infoSW[0].cpu != i_cpu ||
                    Blk_List.Block[neighbour_blk].infoSW[0].blknum != i_blk) return(1);
                for (l = 0; l <= buffer_size-1; ++l) {
	           Blk_List.message_reschange_southwestcorner_recbuf[neighbour_blk][l] =
                     Blk_List.message_reschange_northwestcorner_sendbuf[i_blk][l];
                } /* endfor */
             } else if (Blk_List.Block[i_blk].infoNW[0].dimen.i > 0) {
                if (Blk_List.Block[neighbour_blk].nNE != 1 ||
                    Blk_List.Block[neighbour_blk].infoNE[0].cpu != i_cpu ||
                    Blk_List.Block[neighbour_blk].infoNE[0].blknum != i_blk) return(1);
                for (l = 0; l <= buffer_size-1; ++l) {
	           Blk_List.message_reschange_northeastcorner_recbuf[neighbour_blk][l] =
                     Blk_List.message_reschange_northwestcorner_sendbuf[i_blk][l];
                } /* endfor */
             } else {
                if (Blk_List.Block[neighbour_blk].nNW != 1 ||
                    Blk_List.Block[neighbour_blk].infoNW[0].cpu != i_cpu ||
                    Blk_List.Block[neighbour_blk].infoNW[0].blknum != i_blk) return(1);
                for (l = 0; l <= buffer_size-1; ++l) {
	           Blk_List.message_reschange_northwestcorner_recbuf[neighbour_blk][l] =
                     Blk_List.message_reschange_northwestcorner_sendbuf[i_blk][l];
                } /* endfor */
             } /* endif */

          } /* endif */
       } /* endif */

       // Post receive messages from neighboring blocks at north west corners.
       if (Blk_List.Block[i_blk].used && 
           (Blk_List.Block[i_blk].nNW == 1) &&
           (Blk_List.Block[i_blk].info.level >
            Blk_List.Block[i_blk].infoNW[0].level)) {
	  neighbour_cpu = Blk_List.Block[i_blk].infoNW[0].cpu;
	  neighbour_blk = Blk_List.Block[i_blk].infoNW[0].blknum;
          buffer_size = sqr(Blk_List.Block[i_blk].infoNW[0].dimen.ghost)*
                        Number_of_Solution_Variables;
          if (neighbour_cpu != i_cpu) {
             tag_receive = neighbour_blk*tag_base + tag_NW;
             receive_requests[number_receive_requests] = MPI::COMM_WORLD.Irecv(Blk_List.message_reschange_northwestcorner_recbuf[i_blk],
									       buffer_size,
									       MPI::DOUBLE,
									       neighbour_cpu,
									       tag_receive);
             number_receive_requests++;
          } /* endif */
       } /* endif */

       /////////////////////////////////////////////////////////////////////
       // Exchange messages at the North East corners of solution blocks  //
       // with finer neighbours at higher mesh resolution.                //
       /////////////////////////////////////////////////////////////////////

       // Post send messages to neighboring blocks at north east corners.
       if (Blk_List.Block[i_blk].used && 
           (Blk_List.Block[i_blk].nNE == 1) &&
           (Blk_List.Block[i_blk].info.level <
            Blk_List.Block[i_blk].infoNE[0].level)) {
	  neighbour_cpu = Blk_List.Block[i_blk].infoNE[0].cpu;
	  neighbour_blk = Blk_List.Block[i_blk].infoNE[0].blknum;
          buffer_size = sqr(Blk_List.Block[i_blk].infoNE[0].dimen.ghost)*
                        Number_of_Solution_Variables;
          if (neighbour_cpu != i_cpu) {
             if (Blk_List.Block[i_blk].infoNE[0].dimen.i > 0 &&
                 Blk_List.Block[i_blk].infoNE[0].dimen.j > 0) {
                tag_send    = i_blk*tag_base + tag_SW;
             } else if (Blk_List.Block[i_blk].infoNE[0].dimen.j > 0) {
                tag_send    = i_blk*tag_base + tag_SE;
             } else if (Blk_List.Block[i_blk].infoNE[0].dimen.i > 0) {
                tag_send    = i_blk*tag_base + tag_NW;
             } else {
                tag_send    = i_blk*tag_base + tag_NE;
             } /* endif */
             send_requests[number_send_requests] = MPI::COMM_WORLD.Isend(Blk_List.message_reschange_northeastcorner_sendbuf[i_blk],
									 buffer_size,
									 MPI::DOUBLE,
									 neighbour_cpu,
									 tag_send);
             number_send_requests++;
          } else {
             if (!Blk_List.Block[neighbour_blk].used) return(1); 
             buffer_size_neighbour = sqr(Blk_List.Block[neighbour_blk].info.dimen.ghost)*
                                     Number_of_Solution_Variables;
             if (buffer_size != buffer_size_neighbour) return(1);

             if (Blk_List.Block[i_blk].infoNE[0].dimen.i > 0 &&
                 Blk_List.Block[i_blk].infoNE[0].dimen.j > 0) {
                if (Blk_List.Block[neighbour_blk].nSW != 1 ||
                    Blk_List.Block[neighbour_blk].infoSW[0].cpu != i_cpu ||
                    Blk_List.Block[neighbour_blk].infoSW[0].blknum != i_blk) return(1);
                for (l = 0; l <= buffer_size-1; ++l) {
	           Blk_List.message_reschange_southwestcorner_recbuf[neighbour_blk][l] =
                     Blk_List.message_reschange_northeastcorner_sendbuf[i_blk][l];
                } /* endfor */
             } else if (Blk_List.Block[i_blk].infoNE[0].dimen.j > 0) {
                if (Blk_List.Block[neighbour_blk].nSE != 1 ||
                    Blk_List.Block[neighbour_blk].infoSE[0].cpu != i_cpu ||
                    Blk_List.Block[neighbour_blk].infoSE[0].blknum != i_blk) return(1);
                for (l = 0; l <= buffer_size-1; ++l) {
	           Blk_List.message_reschange_southeastcorner_recbuf[neighbour_blk][l] =
                     Blk_List.message_reschange_northeastcorner_sendbuf[i_blk][l];
                } /* endfor */
             } else if (Blk_List.Block[i_blk].infoNE[0].dimen.i > 0) {
                if (Blk_List.Block[neighbour_blk].nNW != 1 ||
                    Blk_List.Block[neighbour_blk].infoNW[0].cpu != i_cpu ||
                    Blk_List.Block[neighbour_blk].infoNW[0].blknum != i_blk) return(1);
                for (l = 0; l <= buffer_size-1; ++l) {
	           Blk_List.message_reschange_northwestcorner_recbuf[neighbour_blk][l] =
                     Blk_List.message_reschange_northeastcorner_sendbuf[i_blk][l];
                } /* endfor */
             } else {
                if (Blk_List.Block[neighbour_blk].nNE != 1 ||
                    Blk_List.Block[neighbour_blk].infoNE[0].cpu != i_cpu ||
                    Blk_List.Block[neighbour_blk].infoNE[0].blknum != i_blk) return(1);
                for (l = 0; l <= buffer_size-1; ++l) {
	           Blk_List.message_reschange_northeastcorner_recbuf[neighbour_blk][l] =
                     Blk_List.message_reschange_northeastcorner_sendbuf[i_blk][l];
                } /* endfor */
             } /* endif */

          } /* endif */
       } /* endif */

       // Post receive messages from neighboring blocks at north east corners.
       if (Blk_List.Block[i_blk].used && 
           (Blk_List.Block[i_blk].nNE == 1) &&
           (Blk_List.Block[i_blk].info.level >
            Blk_List.Block[i_blk].infoNE[0].level)) {
	  neighbour_cpu = Blk_List.Block[i_blk].infoNE[0].cpu;
	  neighbour_blk = Blk_List.Block[i_blk].infoNE[0].blknum;
          buffer_size = sqr(Blk_List.Block[i_blk].infoNE[0].dimen.ghost)*
                        Number_of_Solution_Variables;
          if (neighbour_cpu != i_cpu) {
             tag_receive = neighbour_blk*tag_base + tag_NE;
             receive_requests[number_receive_requests] = MPI::COMM_WORLD.Irecv(Blk_List.message_reschange_northeastcorner_recbuf[i_blk],
									       buffer_size,
									       MPI::DOUBLE,
									       neighbour_cpu,
									       tag_receive);
             number_receive_requests++;
          } /* endif */
       } /* endif */

       /////////////////////////////////////////////////////////////////////
       // Exchange messages at the South East corners of solution blocks  //
       // with finer neighbours at higher mesh resolution.                //
       /////////////////////////////////////////////////////////////////////

       // Post send messages to neighboring blocks at south east corners.
       if (Blk_List.Block[i_blk].used && 
           (Blk_List.Block[i_blk].nSE == 1) &&
           (Blk_List.Block[i_blk].info.level <
            Blk_List.Block[i_blk].infoSE[0].level)) {
	  neighbour_cpu = Blk_List.Block[i_blk].infoSE[0].cpu;
	  neighbour_blk = Blk_List.Block[i_blk].infoSE[0].blknum;
          buffer_size = sqr(Blk_List.Block[i_blk].infoSE[0].dimen.ghost)*
                        Number_of_Solution_Variables;
          if (neighbour_cpu != i_cpu) {
             if (Blk_List.Block[i_blk].infoSE[0].dimen.i > 0 &&
                 Blk_List.Block[i_blk].infoSE[0].dimen.j > 0) {
                tag_send    = i_blk*tag_base + tag_NW;
             } else if (Blk_List.Block[i_blk].infoSE[0].dimen.j > 0) {
                tag_send    = i_blk*tag_base + tag_NE;
             } else if (Blk_List.Block[i_blk].infoSE[0].dimen.i > 0) {
                tag_send    = i_blk*tag_base + tag_SW;
             } else {
                tag_send    = i_blk*tag_base + tag_SE;
             } /* endif */
             send_requests[number_send_requests] = MPI::COMM_WORLD.Isend(Blk_List.message_reschange_southeastcorner_sendbuf[i_blk],
									 buffer_size,
									 MPI::DOUBLE,
									 neighbour_cpu,
									 tag_send);
             number_send_requests++;
          } else {
             if (!Blk_List.Block[neighbour_blk].used) return(1); 
             buffer_size_neighbour = sqr(Blk_List.Block[neighbour_blk].info.dimen.ghost)*
                                     Number_of_Solution_Variables;
             if (buffer_size != buffer_size_neighbour) return(1);

             if (Blk_List.Block[i_blk].infoSE[0].dimen.i > 0 &&
                 Blk_List.Block[i_blk].infoSE[0].dimen.j > 0) {
                if (Blk_List.Block[neighbour_blk].nNW != 1 ||
                    Blk_List.Block[neighbour_blk].infoNW[0].cpu != i_cpu ||
                    Blk_List.Block[neighbour_blk].infoNW[0].blknum != i_blk) return(1);
                for (l = 0; l <= buffer_size-1; ++l) {
	           Blk_List.message_reschange_northwestcorner_recbuf[neighbour_blk][l] =
                     Blk_List.message_reschange_southeastcorner_sendbuf[i_blk][l];
                } /* endfor */
             } else if (Blk_List.Block[i_blk].infoSE[0].dimen.j > 0) {
                if (Blk_List.Block[neighbour_blk].nNE != 1 ||
                    Blk_List.Block[neighbour_blk].infoNE[0].cpu != i_cpu ||
                    Blk_List.Block[neighbour_blk].infoNE[0].blknum != i_blk) return(1);
                for (l = 0; l <= buffer_size-1; ++l) {
	           Blk_List.message_reschange_northeastcorner_recbuf[neighbour_blk][l] =
                     Blk_List.message_reschange_southeastcorner_sendbuf[i_blk][l];
                } /* endfor */
             } else if (Blk_List.Block[i_blk].infoSE[0].dimen.i > 0) {
                if (Blk_List.Block[neighbour_blk].nSW != 1 ||
                    Blk_List.Block[neighbour_blk].infoSW[0].cpu != i_cpu ||
                    Blk_List.Block[neighbour_blk].infoSW[0].blknum != i_blk) return(1);
                for (l = 0; l <= buffer_size-1; ++l) {
	           Blk_List.message_reschange_southwestcorner_recbuf[neighbour_blk][l] =
                     Blk_List.message_reschange_southeastcorner_sendbuf[i_blk][l];
                } /* endfor */
             } else {
                if (Blk_List.Block[neighbour_blk].nSE != 1 ||
                    Blk_List.Block[neighbour_blk].infoSE[0].cpu != i_cpu ||
                    Blk_List.Block[neighbour_blk].infoSE[0].blknum != i_blk) return(1);
                for (l = 0; l <= buffer_size-1; ++l) {
	           Blk_List.message_reschange_southeastcorner_recbuf[neighbour_blk][l] =
                     Blk_List.message_reschange_southeastcorner_sendbuf[i_blk][l];
                } /* endfor */
             } /* endif */

          } /* endif */
      } /* endif */

       // Post receive messages from neighboring blocks at south east corners.
       if (Blk_List.Block[i_blk].used && 
           (Blk_List.Block[i_blk].nSE == 1) &&
           (Blk_List.Block[i_blk].info.level >
            Blk_List.Block[i_blk].infoSE[0].level)) {
	  neighbour_cpu = Blk_List.Block[i_blk].infoSE[0].cpu;
	  neighbour_blk = Blk_List.Block[i_blk].infoSE[0].blknum;
          buffer_size = sqr(Blk_List.Block[i_blk].infoSE[0].dimen.ghost)*
                        Number_of_Solution_Variables;
          if (neighbour_cpu != i_cpu) {
             tag_receive = neighbour_blk*tag_base + tag_SE;
             receive_requests[number_receive_requests] = MPI::COMM_WORLD.Irecv(Blk_List.message_reschange_southeastcorner_recbuf[i_blk],
									       buffer_size,
									       MPI::DOUBLE,
									       neighbour_cpu,
									       tag_receive);
             number_receive_requests++;
          } /* endif */
       } /* endif */

       /////////////////////////////////////////////////////////////////////
       // Exchange messages at the South West corners of solution blocks  //
       // with finer neighbours at higher mesh resolution.                //
       /////////////////////////////////////////////////////////////////////

       // Post send messages to neighboring blocks at south west corners.
       if (Blk_List.Block[i_blk].used && 
           (Blk_List.Block[i_blk].nSW == 1) &&
           (Blk_List.Block[i_blk].info.level <
            Blk_List.Block[i_blk].infoSW[0].level)) {
	  neighbour_cpu = Blk_List.Block[i_blk].infoSW[0].cpu;
	  neighbour_blk = Blk_List.Block[i_blk].infoSW[0].blknum;
          buffer_size = sqr(Blk_List.Block[i_blk].infoSW[0].dimen.ghost)*
                        Number_of_Solution_Variables;
          if (neighbour_cpu != i_cpu) {
             if (Blk_List.Block[i_blk].infoSW[0].dimen.i > 0 &&
                 Blk_List.Block[i_blk].infoSW[0].dimen.j > 0) {
                tag_send    = i_blk*tag_base + tag_NE;
             } else if (Blk_List.Block[i_blk].infoSW[0].dimen.j > 0) {
                tag_send    = i_blk*tag_base + tag_NW;
             } else if (Blk_List.Block[i_blk].infoSW[0].dimen.i > 0) {
                tag_send    = i_blk*tag_base + tag_SE;
             } else {
                tag_send    = i_blk*tag_base + tag_SW;
             } /* endif */
             send_requests[number_send_requests] = MPI::COMM_WORLD.Isend(Blk_List.message_reschange_southwestcorner_sendbuf[i_blk],
									 buffer_size,
									 MPI::DOUBLE,
									 neighbour_cpu,
									 tag_send);
             number_send_requests++;
          } else {
             if (!Blk_List.Block[neighbour_blk].used) return(1);
             buffer_size_neighbour = sqr(Blk_List.Block[neighbour_blk].info.dimen.ghost)*
                                     Number_of_Solution_Variables;
             if (buffer_size != buffer_size_neighbour) return(1);

             if (Blk_List.Block[i_blk].infoSW[0].dimen.i > 0 &&
                 Blk_List.Block[i_blk].infoSW[0].dimen.j > 0) {
                if (Blk_List.Block[neighbour_blk].nNE != 1 ||
                    Blk_List.Block[neighbour_blk].infoNE[0].cpu != i_cpu ||
                    Blk_List.Block[neighbour_blk].infoNE[0].blknum != i_blk) return(1);
                for (l = 0; l <= buffer_size-1; ++l) {
	           Blk_List.message_reschange_northeastcorner_recbuf[neighbour_blk][l] =
                     Blk_List.message_reschange_southwestcorner_sendbuf[i_blk][l];
                } /* endfor */
             } else if (Blk_List.Block[i_blk].infoSW[0].dimen.j > 0) {
                if (Blk_List.Block[neighbour_blk].nNW != 1 ||
                    Blk_List.Block[neighbour_blk].infoNW[0].cpu != i_cpu ||
                    Blk_List.Block[neighbour_blk].infoNW[0].blknum != i_blk) return(1);
                for (l = 0; l <= buffer_size-1; ++l) {
	           Blk_List.message_reschange_northwestcorner_recbuf[neighbour_blk][l] =
                     Blk_List.message_reschange_southwestcorner_sendbuf[i_blk][l];
                } /* endfor */
             } else if (Blk_List.Block[i_blk].infoSW[0].dimen.i > 0) {
                if (Blk_List.Block[neighbour_blk].nSE != 1 ||
                    Blk_List.Block[neighbour_blk].infoSE[0].cpu != i_cpu ||
                    Blk_List.Block[neighbour_blk].infoSE[0].blknum != i_blk) return(1);
                for (l = 0; l <= buffer_size-1; ++l) {
	           Blk_List.message_reschange_southeastcorner_recbuf[neighbour_blk][l] =
                     Blk_List.message_reschange_southwestcorner_sendbuf[i_blk][l];
                } /* endfor */
             } else {
                if (Blk_List.Block[neighbour_blk].nSW != 1 ||
                    Blk_List.Block[neighbour_blk].infoSW[0].cpu != i_cpu ||
                    Blk_List.Block[neighbour_blk].infoSW[0].blknum != i_blk) return(1);
                for (l = 0; l <= buffer_size-1; ++l) {
	           Blk_List.message_reschange_southwestcorner_recbuf[neighbour_blk][l] =
                     Blk_List.message_reschange_southwestcorner_sendbuf[i_blk][l];
                } /* endfor */
             } /* endif */

          } /* endif */
       } /* endif */

       // Post receive messages from neighboring blocks at south west corners.
       if (Blk_List.Block[i_blk].used && 
           (Blk_List.Block[i_blk].nSW == 1) &&
           (Blk_List.Block[i_blk].info.level >
            Blk_List.Block[i_blk].infoSW[0].level)) {
	  neighbour_cpu = Blk_List.Block[i_blk].infoSW[0].cpu;
	  neighbour_blk = Blk_List.Block[i_blk].infoSW[0].blknum;
          buffer_size = sqr(Blk_List.Block[i_blk].infoSW[0].dimen.ghost)*
                        Number_of_Solution_Variables;
          if (neighbour_cpu != i_cpu) {
             tag_receive = neighbour_blk*tag_base + tag_SW;
             receive_requests[number_receive_requests] = MPI::COMM_WORLD.Irecv(Blk_List.message_reschange_southwestcorner_recbuf[i_blk],
									       buffer_size,
									       MPI::DOUBLE,
									       neighbour_cpu,
									       tag_receive);
             number_receive_requests++;
          } /* endif */
       } /* endif */

    }  /* endfor */

    /* Explicitly Deallocate send buffers memory */
    for(int i=0; i < number_send_requests; i++){
      send_requests[i].Free(); 
    }

    //CFDkit_Barrier_MPI(); // DOES THIS NEED TO BE HERE ??

    /* Wait for all messages to be received. */
    /* Also deallocates recv buffers memory */
    if (number_receive_requests > 0) {
       MPI::Request::Waitall(number_receive_requests, receive_requests);
    } /* endif */

    /* Deallocate memory for message passing requests. */

    delete []receive_requests;
    receive_requests = NULL;
    delete []send_requests;
    send_requests = NULL;

    /* Message passing complete.  Return zero value. */

    return(0);

}





/**********************************************************************
 * Routine: Exchange_Messages_Fine_Grid_Solution_Information          *
 **********************************************************************/
int Exchange_Messages_Fine_Grid_Solution_Information(AdaptiveBlock2D_List &Blk_List,
						     const int Number_of_Solution_Variables) {

  return 0;
}
