/* AdaptiveBlock_NoMPI.cc:  Subroutines for adaptive blocks classes
                            on a single CPU architecture with no 
                            access to MPI. */

/* Include adaptive block header file. */

#ifndef _ADAPTIVBLOCK3D_INCLUDED
#include "AdaptiveBlock3D.h"
#endif // _ADAPTIVEBLOCK3D_INCLUDED

/*************************************************************
 * AdaptiveBlockResourceList -- External subroutines.        *
 *************************************************************/

/*************************************************************
 * AdaptiveBlock3D -- External subroutines.                  *
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
int  AdaptiveBlock3D_List::Exchange_Messages_NoResChange(AdaptiveBlock3D_List &Blk_List,
                                  const int Number_of_Solution_Variables) {
  //cout<<"\n/*in Exchange_Messages_NoResChange"; cout.flush();
    int i_cpu, i_blk, neighbour_cpu, neighbour_blk, 
        buffer_size, buffer_size_neighbour, l;

    /* Perform message passing for the block faces and corners of each 
       solution block having no mesh resolution change. */

    i_cpu = Blk_List.ThisCPU;

    for ( i_blk = 0 ; i_blk <= Blk_List.Nblk-1 ; ++i_blk ) {

       //////////////////////////////////////////////////////////////
       // Exchange messages at the Top faces of solution blocks. //
       //////////////////////////////////////////////////////////////
       if (Blk_List.Block[i_blk].used && 
           (Blk_List.Block[i_blk].nT == 1) &&
           (Blk_List.Block[i_blk].info.level == Blk_List.Block[i_blk].infoT[0].level)) {
	  neighbour_cpu = Blk_List.Block[i_blk].infoT[0].cpu;
	  neighbour_blk = Blk_List.Block[i_blk].infoT[0].blknum;
          if (neighbour_cpu != i_cpu) return(2400);
          if (!Blk_List.Block[neighbour_blk].used) return(2401); 
          buffer_size = Blk_List.Block[i_blk].infoT[0].dimen.ghost*
                        (abs(Blk_List.Block[i_blk].infoT[0].dimen.i)+1)*
                        (abs(Blk_List.Block[i_blk].infoT[0].dimen.j)+1)*
                        Number_of_Solution_Variables+
	    4*Blk_List.Block[i_blk].infoT[0].dimen.ghost*abs(Blk_List.Block[i_blk].infoT[0].dimen.i)+
	    4*Blk_List.Block[i_blk].infoT[0].dimen.ghost*abs(Blk_List.Block[i_blk].infoT[0].dimen.j);
          buffer_size_neighbour = Blk_List.Block[neighbour_blk].info.dimen.ghost*
	    (abs(Blk_List.Block[neighbour_blk].info.dimen.i)+1)*
	    (abs(Blk_List.Block[neighbour_blk].info.dimen.j)+1)*
	    Number_of_Solution_Variables+
	    4*Blk_List.Block[i_blk].infoT[0].dimen.ghost*abs(Blk_List.Block[i_blk].infoT[0].dimen.i)+
	    4*Blk_List.Block[i_blk].infoT[0].dimen.ghost*abs(Blk_List.Block[i_blk].infoT[0].dimen.j);
// 	  cout<<"\nTop face Exchange messages buffer_size = "<<buffer_size<<" buffer_size_neighbour = "<< buffer_size_neighbour;
          if (buffer_size != buffer_size_neighbour) return(2403);
          if (Blk_List.Block[i_blk].infoT[0].dimen.k > 0) {
             if (Blk_List.Block[neighbour_blk].nB != 1 ) return(24041);
              if (   Blk_List.Block[neighbour_blk].infoB[0].cpu != i_cpu ) return(24042);
              if (   Blk_List.Block[neighbour_blk].infoB[0].blknum != i_blk) return(24043);
             for (l = 0; l <= buffer_size-1; ++l) {
	        Blk_List.message_noreschange_bottomface_recbuf[neighbour_blk][l] =
                  Blk_List.message_noreschange_topface_sendbuf[i_blk][l];
             } /* endfor */
          } else {
             if (Blk_List.Block[neighbour_blk].nT != 1 ||
                 Blk_List.Block[neighbour_blk].infoT[0].cpu != i_cpu ||
                 Blk_List.Block[neighbour_blk].infoT[0].blknum != i_blk) return(2405);
             for (l = 0; l <= buffer_size-1; ++l) {
	        Blk_List.message_noreschange_topface_recbuf[neighbour_blk][l] =
                  Blk_List.message_noreschange_topface_sendbuf[i_blk][l];
             } /* endfor */
          } /* endif */
       } /* endif */

       //////////////////////////////////////////////////////////////
       // Exchange messages at the Bottom faces of solution blocks. //
       //////////////////////////////////////////////////////////////
       if (Blk_List.Block[i_blk].used && 
           (Blk_List.Block[i_blk].nB == 1) &&
           (Blk_List.Block[i_blk].info.level == Blk_List.Block[i_blk].infoB[0].level)) {
	  neighbour_cpu = Blk_List.Block[i_blk].infoB[0].cpu;
	  neighbour_blk = Blk_List.Block[i_blk].infoB[0].blknum;
          if (neighbour_cpu != i_cpu) return(2500);
          if (!Blk_List.Block[neighbour_blk].used) return(2501); 
          buffer_size = Blk_List.Block[i_blk].infoB[0].dimen.ghost*
                        (abs(Blk_List.Block[i_blk].infoB[0].dimen.i)+1)*
                        (abs(Blk_List.Block[i_blk].infoB[0].dimen.j)+1)*
                        Number_of_Solution_Variables+
	    4*Blk_List.Block[i_blk].infoB[0].dimen.ghost*abs(Blk_List.Block[i_blk].infoB[0].dimen.i)+
	    4*Blk_List.Block[i_blk].infoB[0].dimen.ghost*abs(Blk_List.Block[i_blk].infoB[0].dimen.j);
          buffer_size_neighbour = Blk_List.Block[neighbour_blk].info.dimen.ghost*
                                  (abs(Blk_List.Block[neighbour_blk].info.dimen.i)+1)*
                                  (abs(Blk_List.Block[neighbour_blk].info.dimen.j)+1)*
                                  Number_of_Solution_Variables+
	    4*Blk_List.Block[i_blk].infoB[0].dimen.ghost*abs(Blk_List.Block[i_blk].infoB[0].dimen.i)+
	    4*Blk_List.Block[i_blk].infoB[0].dimen.ghost*abs(Blk_List.Block[i_blk].infoB[0].dimen.j);
          if (buffer_size != buffer_size_neighbour) return(2503);
          if (Blk_List.Block[i_blk].infoB[0].dimen.k > 0) {
             if (Blk_List.Block[neighbour_blk].nT != 1 ||         
                 Blk_List.Block[neighbour_blk].infoT[0].cpu != i_cpu || 
                 Blk_List.Block[neighbour_blk].infoT[0].blknum != i_blk) return(2504);
             for (l = 0; l <= buffer_size-1; ++l) {
	        Blk_List.message_noreschange_topface_recbuf[neighbour_blk][l] = 
                  Blk_List.message_noreschange_bottomface_sendbuf[i_blk][l];
             } /* endfor */
          } else {
             if (Blk_List.Block[neighbour_blk].nB != 1 ||
                 Blk_List.Block[neighbour_blk].infoB[0].cpu != i_cpu ||
                 Blk_List.Block[neighbour_blk].infoB[0].blknum != i_blk) return(2505);
             for (l = 0; l <= buffer_size-1; ++l) {
	        Blk_List.message_noreschange_bottomface_recbuf[neighbour_blk][l] =
                  Blk_List.message_noreschange_bottomface_sendbuf[i_blk][l];
             } /* endfor */
          } /* endif */
       } /* endif */

       //////////////////////////////////////////////////////////////
       // Exchange messages at the North faces of solution blocks. //
       //////////////////////////////////////////////////////////////
       if (Blk_List.Block[i_blk].used && 
           (Blk_List.Block[i_blk].nN == 1) &&
           (Blk_List.Block[i_blk].info.level == Blk_List.Block[i_blk].infoN[0].level)) {
	  neighbour_cpu = Blk_List.Block[i_blk].infoN[0].cpu;
	  neighbour_blk = Blk_List.Block[i_blk].infoN[0].blknum;
 	  //cout<<"\nNorth face of blk "<<i_blk<<" is "<<neighbour_blk;
         if (neighbour_cpu != i_cpu) return(2200);
          if (!Blk_List.Block[neighbour_blk].used) return(2201); 
          buffer_size = Blk_List.Block[i_blk].infoN[0].dimen.ghost*
                        (abs(Blk_List.Block[i_blk].infoN[0].dimen.i)+1)*
                        (abs(Blk_List.Block[i_blk].infoN[0].dimen.k)+1)*
                        Number_of_Solution_Variables+
	    4*Blk_List.Block[i_blk].infoN[0].dimen.ghost*abs(Blk_List.Block[i_blk].infoN[0].dimen.i)+
	    4*Blk_List.Block[i_blk].infoN[0].dimen.ghost*abs(Blk_List.Block[i_blk].infoN[0].dimen.j);
          buffer_size_neighbour = Blk_List.Block[neighbour_blk].info.dimen.ghost*
                                  (abs(Blk_List.Block[neighbour_blk].info.dimen.i)+1)*
                                  (abs(Blk_List.Block[neighbour_blk].info.dimen.k)+1)*
                                  Number_of_Solution_Variables+
	    4*Blk_List.Block[i_blk].infoN[0].dimen.ghost*abs(Blk_List.Block[i_blk].infoN[0].dimen.i)+
	    4*Blk_List.Block[i_blk].infoN[0].dimen.ghost*abs(Blk_List.Block[i_blk].infoN[0].dimen.j);
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
	  //cout<<"\nSouth face of blk "<<i_blk<<" is "<<neighbour_blk;
          if (neighbour_cpu != i_cpu) return(2100);
          if (!Blk_List.Block[neighbour_blk].used) return(2101); 
          buffer_size = Blk_List.Block[i_blk].infoS[0].dimen.ghost*
                        (abs(Blk_List.Block[i_blk].infoS[0].dimen.i)+1)*
                        (abs(Blk_List.Block[i_blk].infoS[0].dimen.k)+1)*
                        Number_of_Solution_Variables+
	    4*Blk_List.Block[i_blk].infoS[0].dimen.ghost*abs(Blk_List.Block[i_blk].infoS[0].dimen.i)+
	    4*Blk_List.Block[i_blk].infoS[0].dimen.ghost*abs(Blk_List.Block[i_blk].infoS[0].dimen.j);
          buffer_size_neighbour = Blk_List.Block[neighbour_blk].info.dimen.ghost*
                                  (abs(Blk_List.Block[neighbour_blk].info.dimen.i)+1)*
                                  (abs(Blk_List.Block[neighbour_blk].info.dimen.k)+1)*
                                  Number_of_Solution_Variables+
	    4*Blk_List.Block[i_blk].infoS[0].dimen.ghost*abs(Blk_List.Block[i_blk].infoS[0].dimen.i)+
	    4*Blk_List.Block[i_blk].infoS[0].dimen.ghost*abs(Blk_List.Block[i_blk].infoS[0].dimen.j);
          if (buffer_size != buffer_size_neighbour) return(2102);
          if (Blk_List.Block[i_blk].infoS[0].dimen.j > 0) {
             if (Blk_List.Block[neighbour_blk].nN != 1 ) return(neighbour_blk);
              if (   Blk_List.Block[neighbour_blk].infoN[0].cpu != i_cpu ) return(21032);
              if (   Blk_List.Block[neighbour_blk].infoN[0].blknum != i_blk) return(21033);
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
                        (abs(Blk_List.Block[i_blk].infoE[0].dimen.k)+1)*
                        Number_of_Solution_Variables+
	    4*Blk_List.Block[i_blk].infoE[0].dimen.ghost*abs(Blk_List.Block[i_blk].infoE[0].dimen.i)+
	    4*Blk_List.Block[i_blk].infoE[0].dimen.ghost*abs(Blk_List.Block[i_blk].infoE[0].dimen.j);
          buffer_size_neighbour = Blk_List.Block[neighbour_blk].info.dimen.ghost*
                                  (abs(Blk_List.Block[neighbour_blk].info.dimen.j)+1)*
                                  (abs(Blk_List.Block[neighbour_blk].info.dimen.k)+1)*
                                  Number_of_Solution_Variables+
	    4*Blk_List.Block[i_blk].infoE[0].dimen.ghost*abs(Blk_List.Block[i_blk].infoE[0].dimen.i)+
	    4*Blk_List.Block[i_blk].infoE[0].dimen.ghost*abs(Blk_List.Block[i_blk].infoE[0].dimen.j);
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
                        (abs(Blk_List.Block[i_blk].infoW[0].dimen.k)+1)*
                        Number_of_Solution_Variables+
	    4*Blk_List.Block[i_blk].infoW[0].dimen.ghost*abs(Blk_List.Block[i_blk].infoW[0].dimen.i)+
	    4*Blk_List.Block[i_blk].infoW[0].dimen.ghost*abs(Blk_List.Block[i_blk].infoW[0].dimen.j);
          buffer_size_neighbour = Blk_List.Block[neighbour_blk].info.dimen.ghost*
                                  (abs(Blk_List.Block[neighbour_blk].info.dimen.j)+1)*
                                  (abs(Blk_List.Block[neighbour_blk].info.dimen.k)+1)*
                                  Number_of_Solution_Variables+
	    4*Blk_List.Block[i_blk].infoW[0].dimen.ghost*abs(Blk_List.Block[i_blk].infoW[0].dimen.i)+
	    4*Blk_List.Block[i_blk].infoW[0].dimen.ghost*abs(Blk_List.Block[i_blk].infoW[0].dimen.j);
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
                         (abs(Blk_List.Block[i_blk].infoNW[0].dimen.k)+1)*
                        Number_of_Solution_Variables;
          buffer_size_neighbour = sqr(Blk_List.Block[neighbour_blk].info.dimen.ghost)*
                         (abs(Blk_List.Block[i_blk].infoNW[0].dimen.k)+1)*
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
                         (abs(Blk_List.Block[i_blk].infoNE[0].dimen.k)+1)*
                        Number_of_Solution_Variables;
          buffer_size_neighbour = sqr(Blk_List.Block[neighbour_blk].info.dimen.ghost)*
                         (abs(Blk_List.Block[i_blk].infoNE[0].dimen.k)+1)*
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
                         (abs(Blk_List.Block[i_blk].infoSE[0].dimen.k)+1)*
                        Number_of_Solution_Variables;
          buffer_size_neighbour = sqr(Blk_List.Block[neighbour_blk].info.dimen.ghost)*
                         (abs(Blk_List.Block[i_blk].infoSE[0].dimen.k)+1)*
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
                         (abs(Blk_List.Block[i_blk].infoSW[0].dimen.k)+1)*
                        Number_of_Solution_Variables;
          buffer_size_neighbour = sqr(Blk_List.Block[neighbour_blk].info.dimen.ghost)*
                         (abs(Blk_List.Block[i_blk].infoSW[0].dimen.k)+1)*
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


       //******************//
       /////////////////////////////////////////////////////////////////////
       // Exchange messages at the Top North corners of solution blocks. //
       /////////////////////////////////////////////////////////////////////
       if (Blk_List.Block[i_blk].used && 
           (Blk_List.Block[i_blk].nTN == 1) &&
           (Blk_List.Block[i_blk].info.level == Blk_List.Block[i_blk].infoTN[0].level)) {
	  neighbour_cpu = Blk_List.Block[i_blk].infoTN[0].cpu;
	  neighbour_blk = Blk_List.Block[i_blk].infoTN[0].blknum;
          if (neighbour_cpu != i_cpu) return(4300);
          if (!Blk_List.Block[neighbour_blk].used) return(4301); 
          buffer_size = sqr(Blk_List.Block[i_blk].infoTN[0].dimen.ghost)*
                         (abs(Blk_List.Block[i_blk].infoTN[0].dimen.i)+1)*
                        Number_of_Solution_Variables;
          buffer_size_neighbour = sqr(Blk_List.Block[neighbour_blk].info.dimen.ghost)*
                         (abs(Blk_List.Block[i_blk].infoTN[0].dimen.i)+1)*
                                  Number_of_Solution_Variables;
          if (buffer_size != buffer_size_neighbour) return(4302);
          if (Blk_List.Block[i_blk].infoTN[0].dimen.j > 0 &&
	      Blk_List.Block[i_blk].infoTN[0].dimen.k  > 0) {
             if (Blk_List.Block[neighbour_blk].nBS != 1 ||
                 Blk_List.Block[neighbour_blk].infoBS[0].cpu != i_cpu ||
                 Blk_List.Block[neighbour_blk].infoBS[0].blknum != i_blk) return(4303);
             for (l = 0; l <= buffer_size-1; ++l) {
	        Blk_List.message_noreschange_bottomsouthcorner_recbuf[neighbour_blk][l] =
                  Blk_List.message_noreschange_topnorthcorner_sendbuf[i_blk][l];
             } /* endfor */
          } else if (Blk_List.Block[i_blk].infoTN[0].dimen.k > 0) {
             if (Blk_List.Block[neighbour_blk].nBN != 1 ||
                 Blk_List.Block[neighbour_blk].infoBN[0].cpu != i_cpu ||
                 Blk_List.Block[neighbour_blk].infoBN[0].blknum != i_blk) return(4304);
             for (l = 0; l <= buffer_size-1; ++l) {
	        Blk_List.message_noreschange_bottomnorthcorner_recbuf[neighbour_blk][l] =
                  Blk_List.message_noreschange_topnorthcorner_sendbuf[i_blk][l];
             } /* endfor */
          } else if (Blk_List.Block[i_blk].infoTN[0].dimen.j > 0) {
             if (Blk_List.Block[neighbour_blk].nTS != 1 ||
                 Blk_List.Block[neighbour_blk].infoTS[0].cpu != i_cpu ||
                 Blk_List.Block[neighbour_blk].infoTS[0].blknum != i_blk) return(4305);
             for (l = 0; l <= buffer_size-1; ++l) {
	        Blk_List.message_noreschange_topsouthcorner_recbuf[neighbour_blk][l] =
                  Blk_List.message_noreschange_topnorthcorner_sendbuf[i_blk][l];
             } /* endfor */
          } else {
             if (Blk_List.Block[neighbour_blk].nTN != 1 ||
                 Blk_List.Block[neighbour_blk].infoTN[0].cpu != i_cpu ||
                 Blk_List.Block[neighbour_blk].infoTN[0].blknum != i_blk) return(4306);
             for (l = 0; l <= buffer_size-1; ++l) {
	        Blk_List.message_noreschange_topnorthcorner_recbuf[neighbour_blk][l] =
                  Blk_List.message_noreschange_topnorthcorner_sendbuf[i_blk][l];
             } /* endfor */
          } /* endif */
       } /* endif */

      /////////////////////////////////////////////////////////////////////
       // Exchange messages at the Top South corners of solution blocks. //
       /////////////////////////////////////////////////////////////////////
       if (Blk_List.Block[i_blk].used && 
           (Blk_List.Block[i_blk].nTS == 1) &&
           (Blk_List.Block[i_blk].info.level == Blk_List.Block[i_blk].infoTN[0].level)) {
	  neighbour_cpu = Blk_List.Block[i_blk].infoTS[0].cpu;
	  neighbour_blk = Blk_List.Block[i_blk].infoTS[0].blknum;
          if (neighbour_cpu != i_cpu) return(4400);
          if (!Blk_List.Block[neighbour_blk].used) return(4401); 
          buffer_size = sqr(Blk_List.Block[i_blk].infoTS[0].dimen.ghost)*
                         (abs(Blk_List.Block[i_blk].infoTS[0].dimen.i)+1)*
                        Number_of_Solution_Variables;
          buffer_size_neighbour = sqr(Blk_List.Block[neighbour_blk].info.dimen.ghost)*
                         (abs(Blk_List.Block[i_blk].infoTS[0].dimen.i)+1)*
                                  Number_of_Solution_Variables;
          if (buffer_size != buffer_size_neighbour) return(4402);
          if (Blk_List.Block[i_blk].infoTS[0].dimen.j > 0 &&
	      Blk_List.Block[i_blk].infoTS[0].dimen.k  > 0) {
             if (Blk_List.Block[neighbour_blk].nBN != 1 ||
                 Blk_List.Block[neighbour_blk].infoBN[0].cpu != i_cpu ||
                 Blk_List.Block[neighbour_blk].infoBN[0].blknum != i_blk) return(4403);
             for (l = 0; l <= buffer_size-1; ++l) {
	        Blk_List.message_noreschange_bottomnorthcorner_recbuf[neighbour_blk][l] =
                  Blk_List.message_noreschange_topsouthcorner_sendbuf[i_blk][l];
             } /* endfor */
          } else if (Blk_List.Block[i_blk].infoTS[0].dimen.k > 0) {
             if (Blk_List.Block[neighbour_blk].nBS != 1 ||
                 Blk_List.Block[neighbour_blk].infoBS[0].cpu != i_cpu ||
                 Blk_List.Block[neighbour_blk].infoBS[0].blknum != i_blk) return(4404);
             for (l = 0; l <= buffer_size-1; ++l) {
	        Blk_List.message_noreschange_bottomsouthcorner_recbuf[neighbour_blk][l] =
                  Blk_List.message_noreschange_topsouthcorner_sendbuf[i_blk][l];
             } /* endfor */
          } else if (Blk_List.Block[i_blk].infoTS[0].dimen.j > 0) {
             if (Blk_List.Block[neighbour_blk].nTN != 1 ||
                 Blk_List.Block[neighbour_blk].infoTN[0].cpu != i_cpu ||
                 Blk_List.Block[neighbour_blk].infoTN[0].blknum != i_blk) return(4405);
             for (l = 0; l <= buffer_size-1; ++l) {
	        Blk_List.message_noreschange_topnorthcorner_recbuf[neighbour_blk][l] =
                  Blk_List.message_noreschange_topsouthcorner_sendbuf[i_blk][l];
             } /* endfor */
          } else {
             if (Blk_List.Block[neighbour_blk].nTS != 1 ||
                 Blk_List.Block[neighbour_blk].infoTS[0].cpu != i_cpu ||
                 Blk_List.Block[neighbour_blk].infoTS[0].blknum != i_blk) return(4406);
             for (l = 0; l <= buffer_size-1; ++l) {
	        Blk_List.message_noreschange_topsouthcorner_recbuf[neighbour_blk][l] =
                  Blk_List.message_noreschange_topsouthcorner_sendbuf[i_blk][l];
             } /* endfor */
          } /* endif */
       } /* endif */


       /***********************************************/
      /////////////////////////////////////////////////////////////////////
       // Exchange messages at the Top East corners of solution blocks. //
       /////////////////////////////////////////////////////////////////////
       if (Blk_List.Block[i_blk].used && 
           (Blk_List.Block[i_blk].nTE == 1) &&
           (Blk_List.Block[i_blk].info.level == Blk_List.Block[i_blk].infoTE[0].level)) {
	  neighbour_cpu = Blk_List.Block[i_blk].infoTE[0].cpu;
	  neighbour_blk = Blk_List.Block[i_blk].infoTE[0].blknum;
          if (neighbour_cpu != i_cpu) return(4500);
          if (!Blk_List.Block[neighbour_blk].used) return(4501); 
          buffer_size = sqr(Blk_List.Block[i_blk].infoTE[0].dimen.ghost)*
                         (abs(Blk_List.Block[i_blk].infoTE[0].dimen.j)+1)*
                        Number_of_Solution_Variables;
          buffer_size_neighbour = sqr(Blk_List.Block[neighbour_blk].info.dimen.ghost)*
                         (abs(Blk_List.Block[i_blk].infoTE[0].dimen.j)+1)*
                                  Number_of_Solution_Variables;
          if (buffer_size != buffer_size_neighbour) return(4502);
          if (Blk_List.Block[i_blk].infoTE[0].dimen.i > 0 &&
	      Blk_List.Block[i_blk].infoTE[0].dimen.k  > 0) {
             if (Blk_List.Block[neighbour_blk].nBW != 1 ||
                 Blk_List.Block[neighbour_blk].infoBW[0].cpu != i_cpu ||
                 Blk_List.Block[neighbour_blk].infoBW[0].blknum != i_blk) return(4503);
             for (l = 0; l <= buffer_size-1; ++l) {
	        Blk_List.message_noreschange_bottomwestcorner_recbuf[neighbour_blk][l] =
                  Blk_List.message_noreschange_topeastcorner_sendbuf[i_blk][l];
             } /* endfor */
          } else if (Blk_List.Block[i_blk].infoTE[0].dimen.k > 0) {
             if (Blk_List.Block[neighbour_blk].nBE != 1 ||
                 Blk_List.Block[neighbour_blk].infoBE[0].cpu != i_cpu ||
                 Blk_List.Block[neighbour_blk].infoBE[0].blknum != i_blk) return(4504);
             for (l = 0; l <= buffer_size-1; ++l) {
	        Blk_List.message_noreschange_bottomeastcorner_recbuf[neighbour_blk][l] =
                  Blk_List.message_noreschange_topeastcorner_sendbuf[i_blk][l];
             } /* endfor */
          } else if (Blk_List.Block[i_blk].infoTE[0].dimen.i > 0) {
             if (Blk_List.Block[neighbour_blk].nTW != 1 ||
                 Blk_List.Block[neighbour_blk].infoTW[0].cpu != i_cpu ||
                 Blk_List.Block[neighbour_blk].infoTW[0].blknum != i_blk) return(4505);
             for (l = 0; l <= buffer_size-1; ++l) {
	        Blk_List.message_noreschange_topwestcorner_recbuf[neighbour_blk][l] =
                  Blk_List.message_noreschange_topeastcorner_sendbuf[i_blk][l];
             } /* endfor */
          } else {
             if (Blk_List.Block[neighbour_blk].nTE != 1 ||
                 Blk_List.Block[neighbour_blk].infoTE[0].cpu != i_cpu ||
                 Blk_List.Block[neighbour_blk].infoTE[0].blknum != i_blk) return(4506);
             for (l = 0; l <= buffer_size-1; ++l) {
	        Blk_List.message_noreschange_topeastcorner_recbuf[neighbour_blk][l] =
                  Blk_List.message_noreschange_topeastcorner_sendbuf[i_blk][l];
             } /* endfor */
          } /* endif */
       } /* endif */

      /////////////////////////////////////////////////////////////////////
       // Exchange messages at the Top West corners of solution blocks. //
       /////////////////////////////////////////////////////////////////////
       if (Blk_List.Block[i_blk].used && 
           (Blk_List.Block[i_blk].nTW == 1) &&
           (Blk_List.Block[i_blk].info.level == Blk_List.Block[i_blk].infoTW[0].level)) {
	  neighbour_cpu = Blk_List.Block[i_blk].infoTW[0].cpu;
	  neighbour_blk = Blk_List.Block[i_blk].infoTW[0].blknum;
          if (neighbour_cpu != i_cpu) return(4600);
          if (!Blk_List.Block[neighbour_blk].used) return(4601); 
          buffer_size = sqr(Blk_List.Block[i_blk].infoTW[0].dimen.ghost)*
                         (abs(Blk_List.Block[i_blk].infoTW[0].dimen.j)+1)*
                        Number_of_Solution_Variables;
          buffer_size_neighbour = sqr(Blk_List.Block[neighbour_blk].info.dimen.ghost)*
                         (abs(Blk_List.Block[i_blk].infoTW[0].dimen.j)+1)*
                                  Number_of_Solution_Variables;
          if (buffer_size != buffer_size_neighbour) return(4602);
          if (Blk_List.Block[i_blk].infoTW[0].dimen.i > 0 &&
	      Blk_List.Block[i_blk].infoTW[0].dimen.k  > 0) {
             if (Blk_List.Block[neighbour_blk].nBE != 1 ) return(46031);
              if (Blk_List.Block[neighbour_blk].infoBE[0].cpu != i_cpu ) return(46032);
              if ( Blk_List.Block[neighbour_blk].infoBE[0].blknum != i_blk) return(46033);
             for (l = 0; l <= buffer_size-1; ++l) {
	        Blk_List.message_noreschange_bottomeastcorner_recbuf[neighbour_blk][l] =
                  Blk_List.message_noreschange_topwestcorner_sendbuf[i_blk][l];
             } /* endfor */
          } else if (Blk_List.Block[i_blk].infoTW[0].dimen.k > 0) {
             if (Blk_List.Block[neighbour_blk].nBW != 1 ||
                 Blk_List.Block[neighbour_blk].infoBW[0].cpu != i_cpu ||
                 Blk_List.Block[neighbour_blk].infoBW[0].blknum != i_blk) return(4604);
             for (l = 0; l <= buffer_size-1; ++l) {
	        Blk_List.message_noreschange_bottomwestcorner_recbuf[neighbour_blk][l] =
                  Blk_List.message_noreschange_topwestcorner_sendbuf[i_blk][l];
             } /* endfor */
          } else if (Blk_List.Block[i_blk].infoTW[0].dimen.i > 0) {
             if (Blk_List.Block[neighbour_blk].nTE != 1 ||
                 Blk_List.Block[neighbour_blk].infoTE[0].cpu != i_cpu ||
                 Blk_List.Block[neighbour_blk].infoTE[0].blknum != i_blk) return(4605);
             for (l = 0; l <= buffer_size-1; ++l) {
	        Blk_List.message_noreschange_topeastcorner_recbuf[neighbour_blk][l] =
                  Blk_List.message_noreschange_topwestcorner_sendbuf[i_blk][l];
             } /* endfor */
          } else {
             if (Blk_List.Block[neighbour_blk].nTW != 1 ||
                 Blk_List.Block[neighbour_blk].infoTW[0].cpu != i_cpu ||
                 Blk_List.Block[neighbour_blk].infoTW[0].blknum != i_blk) return(4606);
             for (l = 0; l <= buffer_size-1; ++l) {
	        Blk_List.message_noreschange_topwestcorner_recbuf[neighbour_blk][l] =
                  Blk_List.message_noreschange_topwestcorner_sendbuf[i_blk][l];
             } /* endfor */
          } /* endif */
       } /* endif */

       /***********************************************/
       /////////////////////////////////////////////////////////////////////
       // Exchange messages at the bottom North corners of solution blocks. //
       /////////////////////////////////////////////////////////////////////
       if (Blk_List.Block[i_blk].used && 
           (Blk_List.Block[i_blk].nBN == 1) &&
           (Blk_List.Block[i_blk].info.level == Blk_List.Block[i_blk].infoBN[0].level)) {
	  neighbour_cpu = Blk_List.Block[i_blk].infoBN[0].cpu;
	  neighbour_blk = Blk_List.Block[i_blk].infoBN[0].blknum;
          if (neighbour_cpu != i_cpu) return(4700);
          if (!Blk_List.Block[neighbour_blk].used) return(4701); 
          buffer_size = sqr(Blk_List.Block[i_blk].infoBN[0].dimen.ghost)*
                         (abs(Blk_List.Block[i_blk].infoBN[0].dimen.i)+1)*
                        Number_of_Solution_Variables;
          buffer_size_neighbour = sqr(Blk_List.Block[neighbour_blk].info.dimen.ghost)*
                         (abs(Blk_List.Block[i_blk].infoBN[0].dimen.i)+1)*
                                  Number_of_Solution_Variables;
          if (buffer_size != buffer_size_neighbour) return(4702);
          if (Blk_List.Block[i_blk].infoBN[0].dimen.j > 0 &&
	      Blk_List.Block[i_blk].infoBN[0].dimen.k  > 0) {
             if (Blk_List.Block[neighbour_blk].nTS != 1 ||
                 Blk_List.Block[neighbour_blk].infoTS[0].cpu != i_cpu ||
                 Blk_List.Block[neighbour_blk].infoTS[0].blknum != i_blk) return(4703);
             for (l = 0; l <= buffer_size-1; ++l) {
	        Blk_List.message_noreschange_topsouthcorner_recbuf[neighbour_blk][l] =
                  Blk_List.message_noreschange_bottomnorthcorner_sendbuf[i_blk][l];
             } /* endfor */
          } else if (Blk_List.Block[i_blk].infoBN[0].dimen.k > 0) {
             if (Blk_List.Block[neighbour_blk].nTN != 1 ||
                 Blk_List.Block[neighbour_blk].infoTN[0].cpu != i_cpu ||
                 Blk_List.Block[neighbour_blk].infoTN[0].blknum != i_blk) return(4704);
             for (l = 0; l <= buffer_size-1; ++l) {
	        Blk_List.message_noreschange_topnorthcorner_recbuf[neighbour_blk][l] =
                  Blk_List.message_noreschange_bottomnorthcorner_sendbuf[i_blk][l];
             } /* endfor */
          } else if (Blk_List.Block[i_blk].infoBN[0].dimen.j > 0) {
             if (Blk_List.Block[neighbour_blk].nBS != 1 ||
                 Blk_List.Block[neighbour_blk].infoBS[0].cpu != i_cpu ||
                 Blk_List.Block[neighbour_blk].infoBS[0].blknum != i_blk) return(4705);
             for (l = 0; l <= buffer_size-1; ++l) {
	        Blk_List.message_noreschange_bottomsouthcorner_recbuf[neighbour_blk][l] =
                  Blk_List.message_noreschange_bottomnorthcorner_sendbuf[i_blk][l];
             } /* endfor */
          } else {
             if (Blk_List.Block[neighbour_blk].nBN != 1 ||
                 Blk_List.Block[neighbour_blk].infoBN[0].cpu != i_cpu ||
                 Blk_List.Block[neighbour_blk].infoBN[0].blknum != i_blk) return(4706);
             for (l = 0; l <= buffer_size-1; ++l) {
	        Blk_List.message_noreschange_bottomnorthcorner_recbuf[neighbour_blk][l] =
                  Blk_List.message_noreschange_bottomnorthcorner_sendbuf[i_blk][l];
             } /* endfor */
          } /* endif */
       } /* endif */

      /////////////////////////////////////////////////////////////////////
       // Exchange messages at the bottom South corners of solution blocks. //
       /////////////////////////////////////////////////////////////////////
       if (Blk_List.Block[i_blk].used && 
           (Blk_List.Block[i_blk].nBS == 1) &&
           (Blk_List.Block[i_blk].info.level == Blk_List.Block[i_blk].infoBN[0].level)) {
	  neighbour_cpu = Blk_List.Block[i_blk].infoBS[0].cpu;
	  neighbour_blk = Blk_List.Block[i_blk].infoBS[0].blknum;
          if (neighbour_cpu != i_cpu) return(4800);
          if (!Blk_List.Block[neighbour_blk].used) return(4801); 
          buffer_size = sqr(Blk_List.Block[i_blk].infoBS[0].dimen.ghost)*
                         (abs(Blk_List.Block[i_blk].infoBS[0].dimen.i)+1)*
                        Number_of_Solution_Variables;
          buffer_size_neighbour = sqr(Blk_List.Block[neighbour_blk].info.dimen.ghost)*
                         (abs(Blk_List.Block[i_blk].infoBS[0].dimen.i)+1)*
                                  Number_of_Solution_Variables;
          if (buffer_size != buffer_size_neighbour) return(4802);
          if (Blk_List.Block[i_blk].infoBS[0].dimen.j > 0 &&
	      Blk_List.Block[i_blk].infoBS[0].dimen.k  > 0) {
             if (Blk_List.Block[neighbour_blk].nTN != 1 ||
                 Blk_List.Block[neighbour_blk].infoTN[0].cpu != i_cpu ||
                 Blk_List.Block[neighbour_blk].infoTN[0].blknum != i_blk) return(4803);
             for (l = 0; l <= buffer_size-1; ++l) {
	        Blk_List.message_noreschange_topnorthcorner_recbuf[neighbour_blk][l] =
                  Blk_List.message_noreschange_bottomsouthcorner_sendbuf[i_blk][l];
             } /* endfor */
          } else if (Blk_List.Block[i_blk].infoBS[0].dimen.k > 0) {
             if (Blk_List.Block[neighbour_blk].nTS != 1 ||
                 Blk_List.Block[neighbour_blk].infoTS[0].cpu != i_cpu ||
                 Blk_List.Block[neighbour_blk].infoTS[0].blknum != i_blk) return(4804);
             for (l = 0; l <= buffer_size-1; ++l) {
	        Blk_List.message_noreschange_topsouthcorner_recbuf[neighbour_blk][l] =
                  Blk_List.message_noreschange_bottomsouthcorner_sendbuf[i_blk][l];
             } /* endfor */
          } else if (Blk_List.Block[i_blk].infoBS[0].dimen.j > 0) {
             if (Blk_List.Block[neighbour_blk].nBN != 1 ||
                 Blk_List.Block[neighbour_blk].infoBN[0].cpu != i_cpu ||
                 Blk_List.Block[neighbour_blk].infoBN[0].blknum != i_blk) return(4805);
             for (l = 0; l <= buffer_size-1; ++l) {
	        Blk_List.message_noreschange_bottomnorthcorner_recbuf[neighbour_blk][l] =
                  Blk_List.message_noreschange_bottomsouthcorner_sendbuf[i_blk][l];
             } /* endfor */
          } else {
             if (Blk_List.Block[neighbour_blk].nBS != 1 ||
                 Blk_List.Block[neighbour_blk].infoBS[0].cpu != i_cpu ||
                 Blk_List.Block[neighbour_blk].infoBS[0].blknum != i_blk) return(4806);
             for (l = 0; l <= buffer_size-1; ++l) {
	        Blk_List.message_noreschange_bottomsouthcorner_recbuf[neighbour_blk][l] =
                  Blk_List.message_noreschange_bottomsouthcorner_sendbuf[i_blk][l];
             } /* endfor */
          } /* endif */
       } /* endif */


       /***********************************************/
      /////////////////////////////////////////////////////////////////////
       // Exchange messages at the bottom East corners of solution blocks. //
       /////////////////////////////////////////////////////////////////////
       if (Blk_List.Block[i_blk].used && 
           (Blk_List.Block[i_blk].nBE == 1) &&
           (Blk_List.Block[i_blk].info.level == Blk_List.Block[i_blk].infoBE[0].level)) {
	  neighbour_cpu = Blk_List.Block[i_blk].infoBE[0].cpu;
	  neighbour_blk = Blk_List.Block[i_blk].infoBE[0].blknum;
          if (neighbour_cpu != i_cpu) return(4900);
          if (!Blk_List.Block[neighbour_blk].used) return(4901); 
          buffer_size = sqr(Blk_List.Block[i_blk].infoBE[0].dimen.ghost)*
                         (abs(Blk_List.Block[i_blk].infoBE[0].dimen.j)+1)*
                        Number_of_Solution_Variables;
          buffer_size_neighbour = sqr(Blk_List.Block[neighbour_blk].info.dimen.ghost)*
                         (abs(Blk_List.Block[i_blk].infoBE[0].dimen.j)+1)*
                                  Number_of_Solution_Variables;
          if (buffer_size != buffer_size_neighbour) return(4902);
          if (Blk_List.Block[i_blk].infoBE[0].dimen.i > 0 &&
	      Blk_List.Block[i_blk].infoBE[0].dimen.k  > 0) {
             if (Blk_List.Block[neighbour_blk].nTW != 1 ||
                 Blk_List.Block[neighbour_blk].infoTW[0].cpu != i_cpu ||
                 Blk_List.Block[neighbour_blk].infoTW[0].blknum != i_blk) return(4903);
             for (l = 0; l <= buffer_size-1; ++l) {
	        Blk_List.message_noreschange_topwestcorner_recbuf[neighbour_blk][l] =
                  Blk_List.message_noreschange_bottomeastcorner_sendbuf[i_blk][l];
             } /* endfor */
          } else if (Blk_List.Block[i_blk].infoBE[0].dimen.k > 0) {
             if (Blk_List.Block[neighbour_blk].nTE != 1 ||
                 Blk_List.Block[neighbour_blk].infoTE[0].cpu != i_cpu ||
                 Blk_List.Block[neighbour_blk].infoTE[0].blknum != i_blk) return(4904);
             for (l = 0; l <= buffer_size-1; ++l) {
	        Blk_List.message_noreschange_topeastcorner_recbuf[neighbour_blk][l] =
                  Blk_List.message_noreschange_bottomeastcorner_sendbuf[i_blk][l];
             } /* endfor */
          } else if (Blk_List.Block[i_blk].infoBE[0].dimen.i > 0) {
             if (Blk_List.Block[neighbour_blk].nBW != 1 ||
                 Blk_List.Block[neighbour_blk].infoBW[0].cpu != i_cpu ||
                 Blk_List.Block[neighbour_blk].infoBW[0].blknum != i_blk) return(4905);
             for (l = 0; l <= buffer_size-1; ++l) {
	        Blk_List.message_noreschange_bottomwestcorner_recbuf[neighbour_blk][l] =
                  Blk_List.message_noreschange_bottomeastcorner_sendbuf[i_blk][l];
             } /* endfor */
          } else {
             if (Blk_List.Block[neighbour_blk].nBE != 1 ||
                 Blk_List.Block[neighbour_blk].infoBE[0].cpu != i_cpu ||
                 Blk_List.Block[neighbour_blk].infoBE[0].blknum != i_blk) return(4906);
             for (l = 0; l <= buffer_size-1; ++l) {
	        Blk_List.message_noreschange_bottomeastcorner_recbuf[neighbour_blk][l] =
                  Blk_List.message_noreschange_bottomeastcorner_sendbuf[i_blk][l];
             } /* endfor */
          } /* endif */
       } /* endif */

      /////////////////////////////////////////////////////////////////////
       // Exchange messages at the bottom West corners of solution blocks. //
       /////////////////////////////////////////////////////////////////////
       if (Blk_List.Block[i_blk].used && 
           (Blk_List.Block[i_blk].nBW == 1) &&
           (Blk_List.Block[i_blk].info.level == Blk_List.Block[i_blk].infoBW[0].level)) {
	  neighbour_cpu = Blk_List.Block[i_blk].infoBW[0].cpu;
	  neighbour_blk = Blk_List.Block[i_blk].infoBW[0].blknum;
          if (neighbour_cpu != i_cpu) return(5000);
          if (!Blk_List.Block[neighbour_blk].used) return(5001); 
          buffer_size = sqr(Blk_List.Block[i_blk].infoBW[0].dimen.ghost)*
                         (abs(Blk_List.Block[i_blk].infoBW[0].dimen.j)+1)*
                        Number_of_Solution_Variables;
          buffer_size_neighbour = sqr(Blk_List.Block[neighbour_blk].info.dimen.ghost)*
                         (abs(Blk_List.Block[i_blk].infoBW[0].dimen.j)+1)*
                                  Number_of_Solution_Variables;
          if (buffer_size != buffer_size_neighbour) return(5002);
          if (Blk_List.Block[i_blk].infoBW[0].dimen.i > 0 &&
	      Blk_List.Block[i_blk].infoBW[0].dimen.k  > 0) {
             if (Blk_List.Block[neighbour_blk].nTE != 1 ||
                 Blk_List.Block[neighbour_blk].infoTE[0].cpu != i_cpu ||
                 Blk_List.Block[neighbour_blk].infoTE[0].blknum != i_blk) return(5003);
             for (l = 0; l <= buffer_size-1; ++l) {
	        Blk_List.message_noreschange_topeastcorner_recbuf[neighbour_blk][l] =
                  Blk_List.message_noreschange_bottomwestcorner_sendbuf[i_blk][l];
             } /* endfor */
          } else if (Blk_List.Block[i_blk].infoBW[0].dimen.k > 0) {
             if (Blk_List.Block[neighbour_blk].nTW != 1 ||
                 Blk_List.Block[neighbour_blk].infoTW[0].cpu != i_cpu ||
                 Blk_List.Block[neighbour_blk].infoTW[0].blknum != i_blk) return(5004);
             for (l = 0; l <= buffer_size-1; ++l) {
	        Blk_List.message_noreschange_topwestcorner_recbuf[neighbour_blk][l] =
                  Blk_List.message_noreschange_bottomwestcorner_sendbuf[i_blk][l];
             } /* endfor */
          } else if (Blk_List.Block[i_blk].infoBW[0].dimen.i > 0) {
             if (Blk_List.Block[neighbour_blk].nBE != 1 ||
                 Blk_List.Block[neighbour_blk].infoBE[0].cpu != i_cpu ||
                 Blk_List.Block[neighbour_blk].infoBE[0].blknum != i_blk) return(5005);
             for (l = 0; l <= buffer_size-1; ++l) {
	        Blk_List.message_noreschange_bottomeastcorner_recbuf[neighbour_blk][l] =
                  Blk_List.message_noreschange_bottomwestcorner_sendbuf[i_blk][l];
             } /* endfor */
          } else {
             if (Blk_List.Block[neighbour_blk].nBW != 1 ||
                 Blk_List.Block[neighbour_blk].infoBW[0].cpu != i_cpu ||
                 Blk_List.Block[neighbour_blk].infoBW[0].blknum != i_blk) return(5006);
             for (l = 0; l <= buffer_size-1; ++l) {
	        Blk_List.message_noreschange_bottomwestcorner_recbuf[neighbour_blk][l] =
                  Blk_List.message_noreschange_bottomwestcorner_sendbuf[i_blk][l];
             } /* endfor */
          } /* endif */
       } /* endif */

       /***********************************************/
       /////////////////////////////////////////////////////////////////////
       // Exchange messages at the Top North West corner of solution blocks. //
       /////////////////////////////////////////////////////////////////////
       if (Blk_List.Block[i_blk].used && 
           (Blk_List.Block[i_blk].nTNW == 1) &&
           (Blk_List.Block[i_blk].info.level == Blk_List.Block[i_blk].infoTNW[0].level)) {
	  neighbour_cpu = Blk_List.Block[i_blk].infoTNW[0].cpu;
	  neighbour_blk = Blk_List.Block[i_blk].infoTNW[0].blknum;
          if (neighbour_cpu != i_cpu) return(5100);
          if (!Blk_List.Block[neighbour_blk].used) return(5101); 
          buffer_size = sqr(Blk_List.Block[i_blk].infoTNW[0].dimen.ghost)*
	    Blk_List.Block[i_blk].infoTNW[0].dimen.ghost*
                        Number_of_Solution_Variables;
          buffer_size_neighbour = sqr(Blk_List.Block[neighbour_blk].info.dimen.ghost)*
	    Blk_List.Block[neighbour_blk].info.dimen.ghost*
                                  Number_of_Solution_Variables;
          if (buffer_size != buffer_size_neighbour) return(5102);
          if (Blk_List.Block[i_blk].infoTNW[0].dimen.i > 0 &&
	      Blk_List.Block[i_blk].infoTNW[0].dimen.j > 0 &&
	      Blk_List.Block[i_blk].infoTNW[0].dimen.k  > 0) {
             if (Blk_List.Block[neighbour_blk].nBSE != 1 ||
                 Blk_List.Block[neighbour_blk].infoBSE[0].cpu != i_cpu ||
                 Blk_List.Block[neighbour_blk].infoBSE[0].blknum != i_blk) return(5103);
             for (l = 0; l <= buffer_size-1; ++l) {
	        Blk_List.message_noreschange_bottomsoutheastcorner_recbuf[neighbour_blk][l] =
                  Blk_List.message_noreschange_topnorthwestcorner_sendbuf[i_blk][l];
             } /* endfor */
          } else if (Blk_List.Block[i_blk].infoTNW[0].dimen.i > 0 &&
		     Blk_List.Block[i_blk].infoTNW[0].dimen.j  > 0) {
             if (Blk_List.Block[neighbour_blk].nTSE != 1 ||
                 Blk_List.Block[neighbour_blk].infoTSE[0].cpu != i_cpu ||
                 Blk_List.Block[neighbour_blk].infoTSE[0].blknum != i_blk) return(5104);
             for (l = 0; l <= buffer_size-1; ++l) {
	        Blk_List.message_noreschange_topsoutheastcorner_recbuf[neighbour_blk][l] =
                  Blk_List.message_noreschange_topnorthwestcorner_sendbuf[i_blk][l];
             } /* endfor */
          } else if (Blk_List.Block[i_blk].infoTNW[0].dimen.j > 0 &&
		     Blk_List.Block[i_blk].infoTNW[0].dimen.k  > 0) {
             if (Blk_List.Block[neighbour_blk].nBSW != 1 ||
                 Blk_List.Block[neighbour_blk].infoBSW[0].cpu != i_cpu ||
                 Blk_List.Block[neighbour_blk].infoBSW[0].blknum != i_blk) return(5105);
             for (l = 0; l <= buffer_size-1; ++l) {
	        Blk_List.message_noreschange_bottomsouthwestcorner_recbuf[neighbour_blk][l] =
                  Blk_List.message_noreschange_topnorthwestcorner_sendbuf[i_blk][l];
             } /* endfor */
          } else if (Blk_List.Block[i_blk].infoTNW[0].dimen.k > 0 &&
		     Blk_List.Block[i_blk].infoTNW[0].dimen.i  > 0) {
             if (Blk_List.Block[neighbour_blk].nBNE != 1 ||
                 Blk_List.Block[neighbour_blk].infoBNE[0].cpu != i_cpu ||
                 Blk_List.Block[neighbour_blk].infoBNE[0].blknum != i_blk) return(5105);
             for (l = 0; l <= buffer_size-1; ++l) {
	        Blk_List.message_noreschange_bottomnortheastcorner_recbuf[neighbour_blk][l] =
                  Blk_List.message_noreschange_topnorthwestcorner_sendbuf[i_blk][l];
             } /* endfor */
          } else if (Blk_List.Block[i_blk].infoTNW[0].dimen.i > 0 ) {
             if (Blk_List.Block[neighbour_blk].nTNE != 1 ||
                 Blk_List.Block[neighbour_blk].infoTNE[0].cpu != i_cpu ||
                 Blk_List.Block[neighbour_blk].infoTNE[0].blknum != i_blk) return(5107);
             for (l = 0; l <= buffer_size-1; ++l) {
	        Blk_List.message_noreschange_topnortheastcorner_recbuf[neighbour_blk][l] =
                  Blk_List.message_noreschange_topnorthwestcorner_sendbuf[i_blk][l];
             } /* endfor */
          } else if (Blk_List.Block[i_blk].infoTNW[0].dimen.j > 0 ) {
             if (Blk_List.Block[neighbour_blk].nTSW != 1 ||
                 Blk_List.Block[neighbour_blk].infoTSW[0].cpu != i_cpu ||
                 Blk_List.Block[neighbour_blk].infoTSW[0].blknum != i_blk) return(5108);
             for (l = 0; l <= buffer_size-1; ++l) {
	        Blk_List.message_noreschange_topsouthwestcorner_recbuf[neighbour_blk][l] =
                  Blk_List.message_noreschange_topnorthwestcorner_sendbuf[i_blk][l];
             } /* endfor */
          } else if (Blk_List.Block[i_blk].infoTNW[0].dimen.k > 0 ) {
             if (Blk_List.Block[neighbour_blk].nBNW != 1 ||
                 Blk_List.Block[neighbour_blk].infoBNW[0].cpu != i_cpu ||
                 Blk_List.Block[neighbour_blk].infoBNW[0].blknum != i_blk) return(5109);
             for (l = 0; l <= buffer_size-1; ++l) {
	        Blk_List.message_noreschange_bottomnorthwestcorner_recbuf[neighbour_blk][l] =
                  Blk_List.message_noreschange_topnorthwestcorner_sendbuf[i_blk][l];
             } /* endfor */
          } else {
             if (Blk_List.Block[neighbour_blk].nTNW != 1 ||
                 Blk_List.Block[neighbour_blk].infoTNW[0].cpu != i_cpu ||
                 Blk_List.Block[neighbour_blk].infoTNW[0].blknum != i_blk) return(5110);
             for (l = 0; l <= buffer_size-1; ++l) {
	        Blk_List.message_noreschange_topnorthwestcorner_recbuf[neighbour_blk][l] =
                  Blk_List.message_noreschange_topnorthwestcorner_sendbuf[i_blk][l];
             } /* endfor */
          } /* endif */
       } /* endif */



       /////////////////////////////////////////////////////////////////////
       // Exchange messages at the Top South West corner of solution blocks. //
       /////////////////////////////////////////////////////////////////////
       if (Blk_List.Block[i_blk].used && 
           (Blk_List.Block[i_blk].nTSW == 1) &&
           (Blk_List.Block[i_blk].info.level == Blk_List.Block[i_blk].infoTSW[0].level)) {
	  neighbour_cpu = Blk_List.Block[i_blk].infoTSW[0].cpu;
	  neighbour_blk = Blk_List.Block[i_blk].infoTSW[0].blknum;
          if (neighbour_cpu != i_cpu) return(5400);
          if (!Blk_List.Block[neighbour_blk].used) return(5401); 
          buffer_size = sqr(Blk_List.Block[i_blk].infoTSW[0].dimen.ghost)*
	    Blk_List.Block[i_blk].infoTSW[0].dimen.ghost*
                        Number_of_Solution_Variables;
          buffer_size_neighbour = sqr(Blk_List.Block[neighbour_blk].info.dimen.ghost)*
	    Blk_List.Block[neighbour_blk].info.dimen.ghost*
                                  Number_of_Solution_Variables;
          if (buffer_size != buffer_size_neighbour) return(4202);
          if (Blk_List.Block[i_blk].infoTSW[0].dimen.i > 0 &&
	      Blk_List.Block[i_blk].infoTSW[0].dimen.j > 0 &&
	      Blk_List.Block[i_blk].infoTSW[0].dimen.k  > 0) {
             if (Blk_List.Block[neighbour_blk].nBNE != 1 ||
                 Blk_List.Block[neighbour_blk].infoBNE[0].cpu != i_cpu ||
                 Blk_List.Block[neighbour_blk].infoBNE[0].blknum != i_blk) return(5403);
             for (l = 0; l <= buffer_size-1; ++l) {
	        Blk_List.message_noreschange_bottomnortheastcorner_recbuf[neighbour_blk][l] =
                  Blk_List.message_noreschange_topsouthwestcorner_sendbuf[i_blk][l];
             } /* endfor */
          } else if (Blk_List.Block[i_blk].infoTSW[0].dimen.i > 0 &&
		     Blk_List.Block[i_blk].infoTSW[0].dimen.j  > 0) {
             if (Blk_List.Block[neighbour_blk].nTNE != 1 ||
                 Blk_List.Block[neighbour_blk].infoTNE[0].cpu != i_cpu ||
                 Blk_List.Block[neighbour_blk].infoTNE[0].blknum != i_blk) return(5404);
             for (l = 0; l <= buffer_size-1; ++l) {
	        Blk_List.message_noreschange_topnortheastcorner_recbuf[neighbour_blk][l] =
                  Blk_List.message_noreschange_topsouthwestcorner_sendbuf[i_blk][l];
             } /* endfor */
          } else if (Blk_List.Block[i_blk].infoTSW[0].dimen.j > 0 &&
		     Blk_List.Block[i_blk].infoTSW[0].dimen.k  > 0) {
             if (Blk_List.Block[neighbour_blk].nBNW != 1 ||
                 Blk_List.Block[neighbour_blk].infoBNW[0].cpu != i_cpu ||
                 Blk_List.Block[neighbour_blk].infoBNW[0].blknum != i_blk) return(5405);
             for (l = 0; l <= buffer_size-1; ++l) {
	        Blk_List.message_noreschange_bottomnorthwestcorner_recbuf[neighbour_blk][l] =
                  Blk_List.message_noreschange_topsouthwestcorner_sendbuf[i_blk][l];
             } /* endfor */
          } else if (Blk_List.Block[i_blk].infoTSW[0].dimen.k > 0 &&
		     Blk_List.Block[i_blk].infoTSW[0].dimen.i  > 0) {
             if (Blk_List.Block[neighbour_blk].nBSE != 1 ||
                 Blk_List.Block[neighbour_blk].infoBSE[0].cpu != i_cpu ||
                 Blk_List.Block[neighbour_blk].infoBSE[0].blknum != i_blk) return(5406);
             for (l = 0; l <= buffer_size-1; ++l) {
	        Blk_List.message_noreschange_bottomsoutheastcorner_recbuf[neighbour_blk][l] =
                  Blk_List.message_noreschange_topsouthwestcorner_sendbuf[i_blk][l];
             } /* endfor */
          } else if (Blk_List.Block[i_blk].infoTSW[0].dimen.i > 0 ) {
             if (Blk_List.Block[neighbour_blk].nTSE != 1 ||
                 Blk_List.Block[neighbour_blk].infoTSE[0].cpu != i_cpu ||
                 Blk_List.Block[neighbour_blk].infoTSE[0].blknum != i_blk) return(5407);
             for (l = 0; l <= buffer_size-1; ++l) {
	        Blk_List.message_noreschange_topsoutheastcorner_recbuf[neighbour_blk][l] =
                  Blk_List.message_noreschange_topsouthwestcorner_sendbuf[i_blk][l];
             } /* endfor */
          } else if (Blk_List.Block[i_blk].infoTSW[0].dimen.j > 0 ) {
             if (Blk_List.Block[neighbour_blk].nTNW != 1 ||
                 Blk_List.Block[neighbour_blk].infoTNW[0].cpu != i_cpu ||
                 Blk_List.Block[neighbour_blk].infoTNW[0].blknum != i_blk) return(5408);
             for (l = 0; l <= buffer_size-1; ++l) {
	        Blk_List.message_noreschange_topnorthwestcorner_recbuf[neighbour_blk][l] =
                  Blk_List.message_noreschange_topsouthwestcorner_sendbuf[i_blk][l];
             } /* endfor */
          } else if (Blk_List.Block[i_blk].infoTSW[0].dimen.k > 0 ) {
             if (Blk_List.Block[neighbour_blk].nBSW != 1 ||
                 Blk_List.Block[neighbour_blk].infoBSW[0].cpu != i_cpu ||
                 Blk_List.Block[neighbour_blk].infoBSW[0].blknum != i_blk) return(5409);
             for (l = 0; l <= buffer_size-1; ++l) {
	        Blk_List.message_noreschange_bottomsouthwestcorner_recbuf[neighbour_blk][l] =
                  Blk_List.message_noreschange_topsouthwestcorner_sendbuf[i_blk][l];
             } /* endfor */
          } else {
             if (Blk_List.Block[neighbour_blk].nTSW != 1 ||
                 Blk_List.Block[neighbour_blk].infoTSW[0].cpu != i_cpu ||
                 Blk_List.Block[neighbour_blk].infoTSW[0].blknum != i_blk) return(5410);
             for (l = 0; l <= buffer_size-1; ++l) {
	        Blk_List.message_noreschange_topsouthwestcorner_recbuf[neighbour_blk][l] =
                  Blk_List.message_noreschange_topsouthwestcorner_sendbuf[i_blk][l];
             } /* endfor */
          } /* endif */
       } /* endif */

       /////////////////////////////////////////////////////////////////////
       // Exchange messages at the Top South East corner of solution blocks. //
       /////////////////////////////////////////////////////////////////////
       if (Blk_List.Block[i_blk].used && 
           (Blk_List.Block[i_blk].nTSE == 1) &&
           (Blk_List.Block[i_blk].info.level == Blk_List.Block[i_blk].infoTSE[0].level)) {
	  neighbour_cpu = Blk_List.Block[i_blk].infoTSE[0].cpu;
	  neighbour_blk = Blk_List.Block[i_blk].infoTSE[0].blknum;
          if (neighbour_cpu != i_cpu) return(5500);
          if (!Blk_List.Block[neighbour_blk].used) return(5501); 
          buffer_size = sqr(Blk_List.Block[i_blk].infoTSE[0].dimen.ghost)*
	    Blk_List.Block[i_blk].infoTSE[0].dimen.ghost*
                        Number_of_Solution_Variables;
          buffer_size_neighbour = sqr(Blk_List.Block[neighbour_blk].info.dimen.ghost)*
	    Blk_List.Block[neighbour_blk].info.dimen.ghost*
                                  Number_of_Solution_Variables;
          if (buffer_size != buffer_size_neighbour) return(5502);
          if (Blk_List.Block[i_blk].infoTSE[0].dimen.i > 0 &&
	      Blk_List.Block[i_blk].infoTSE[0].dimen.j > 0 &&
	      Blk_List.Block[i_blk].infoTSE[0].dimen.k  > 0) {
             if (Blk_List.Block[neighbour_blk].nBNW != 1 ||
                 Blk_List.Block[neighbour_blk].infoBNW[0].cpu != i_cpu ||
                 Blk_List.Block[neighbour_blk].infoBNW[0].blknum != i_blk) return(5503);
             for (l = 0; l <= buffer_size-1; ++l) {
	        Blk_List.message_noreschange_bottomnorthwestcorner_recbuf[neighbour_blk][l] =
                  Blk_List.message_noreschange_topsoutheastcorner_sendbuf[i_blk][l];
             } /* endfor */
          } else if (Blk_List.Block[i_blk].infoTSE[0].dimen.i > 0 &&
		     Blk_List.Block[i_blk].infoTSE[0].dimen.j  > 0) {
             if (Blk_List.Block[neighbour_blk].nTNW != 1 ||
                 Blk_List.Block[neighbour_blk].infoTNW[0].cpu != i_cpu ||
                 Blk_List.Block[neighbour_blk].infoTNW[0].blknum != i_blk) return(5504);
             for (l = 0; l <= buffer_size-1; ++l) {
	        Blk_List.message_noreschange_topnorthwestcorner_recbuf[neighbour_blk][l] =
                  Blk_List.message_noreschange_topsoutheastcorner_sendbuf[i_blk][l];
             } /* endfor */
          } else if (Blk_List.Block[i_blk].infoTSE[0].dimen.j > 0 &&
		     Blk_List.Block[i_blk].infoTSE[0].dimen.k  > 0) {
             if (Blk_List.Block[neighbour_blk].nBNE != 1 ||
                 Blk_List.Block[neighbour_blk].infoBNE[0].cpu != i_cpu ||
                 Blk_List.Block[neighbour_blk].infoBNE[0].blknum != i_blk) return(5505);
             for (l = 0; l <= buffer_size-1; ++l) {
	        Blk_List.message_noreschange_bottomnortheastcorner_recbuf[neighbour_blk][l] =
                  Blk_List.message_noreschange_topsoutheastcorner_sendbuf[i_blk][l];
             } /* endfor */
          } else if (Blk_List.Block[i_blk].infoTSE[0].dimen.k > 0 &&
		     Blk_List.Block[i_blk].infoTSE[0].dimen.i  > 0) {
             if (Blk_List.Block[neighbour_blk].nBSW != 1 ||
                 Blk_List.Block[neighbour_blk].infoBSW[0].cpu != i_cpu ||
                 Blk_List.Block[neighbour_blk].infoBSW[0].blknum != i_blk) return(5506);
             for (l = 0; l <= buffer_size-1; ++l) {
	        Blk_List.message_noreschange_bottomsouthwestcorner_recbuf[neighbour_blk][l] =
                  Blk_List.message_noreschange_topsoutheastcorner_sendbuf[i_blk][l];
             } /* endfor */
          } else if (Blk_List.Block[i_blk].infoTSE[0].dimen.i > 0 ) {
             if (Blk_List.Block[neighbour_blk].nTSW != 1 ||
                 Blk_List.Block[neighbour_blk].infoTSW[0].cpu != i_cpu ||
                 Blk_List.Block[neighbour_blk].infoTSW[0].blknum != i_blk) return(5507);
             for (l = 0; l <= buffer_size-1; ++l) {
	        Blk_List.message_noreschange_topsouthwestcorner_recbuf[neighbour_blk][l] =
                  Blk_List.message_noreschange_topsoutheastcorner_sendbuf[i_blk][l];
             } /* endfor */
          } else if (Blk_List.Block[i_blk].infoTSE[0].dimen.j > 0 ) {
             if (Blk_List.Block[neighbour_blk].nTNE != 1 ||
                 Blk_List.Block[neighbour_blk].infoTNE[0].cpu != i_cpu ||
                 Blk_List.Block[neighbour_blk].infoTNE[0].blknum != i_blk) return(5508);
             for (l = 0; l <= buffer_size-1; ++l) {
	        Blk_List.message_noreschange_topnortheastcorner_recbuf[neighbour_blk][l] =
                  Blk_List.message_noreschange_topsoutheastcorner_sendbuf[i_blk][l];
             } /* endfor */
          } else if (Blk_List.Block[i_blk].infoTSE[0].dimen.k > 0 ) {
             if (Blk_List.Block[neighbour_blk].nBSE != 1 ||
                 Blk_List.Block[neighbour_blk].infoBSE[0].cpu != i_cpu ||
                 Blk_List.Block[neighbour_blk].infoBSE[0].blknum != i_blk) return(5509);
             for (l = 0; l <= buffer_size-1; ++l) {
	        Blk_List.message_noreschange_bottomsoutheastcorner_recbuf[neighbour_blk][l] =
                  Blk_List.message_noreschange_topsoutheastcorner_sendbuf[i_blk][l];
             } /* endfor */
          } else {
             if (Blk_List.Block[neighbour_blk].nTSE != 1 ||
                 Blk_List.Block[neighbour_blk].infoTSE[0].cpu != i_cpu ||
                 Blk_List.Block[neighbour_blk].infoTSE[0].blknum != i_blk) return(5510);
             for (l = 0; l <= buffer_size-1; ++l) {
	        Blk_List.message_noreschange_topsoutheastcorner_recbuf[neighbour_blk][l] =
                  Blk_List.message_noreschange_topsoutheastcorner_sendbuf[i_blk][l];
             } /* endfor */
          } /* endif */
       } /* endif */


         /***********************************************/
        /////////////////////////////////////////////////////////////////////
       // Exchange messages at the Bottom North West corner of solution blocks. //
       /////////////////////////////////////////////////////////////////////
       if (Blk_List.Block[i_blk].used && 
           (Blk_List.Block[i_blk].nBNW == 1) &&
           (Blk_List.Block[i_blk].info.level == Blk_List.Block[i_blk].infoBNW[0].level)) {
	  neighbour_cpu = Blk_List.Block[i_blk].infoBNW[0].cpu;
	  neighbour_blk = Blk_List.Block[i_blk].infoBNW[0].blknum;
          if (neighbour_cpu != i_cpu) return(5600);
          if (!Blk_List.Block[neighbour_blk].used) return(5601); 
          buffer_size = sqr(Blk_List.Block[i_blk].infoBNW[0].dimen.ghost)*
	    Blk_List.Block[i_blk].infoBNW[0].dimen.ghost*
                        Number_of_Solution_Variables;
          buffer_size_neighbour = sqr(Blk_List.Block[neighbour_blk].info.dimen.ghost)*
	    Blk_List.Block[neighbour_blk].info.dimen.ghost*
                                  Number_of_Solution_Variables;
          if (buffer_size != buffer_size_neighbour) return(5602);
          if (Blk_List.Block[i_blk].infoBNW[0].dimen.i > 0 &&
	      Blk_List.Block[i_blk].infoBNW[0].dimen.j > 0 &&
	      Blk_List.Block[i_blk].infoBNW[0].dimen.k  > 0) {
             if (Blk_List.Block[neighbour_blk].nTSE != 1 ||
                 Blk_List.Block[neighbour_blk].infoTSE[0].cpu != i_cpu ||
                 Blk_List.Block[neighbour_blk].infoTSE[0].blknum != i_blk) return(5603);
             for (l = 0; l <= buffer_size-1; ++l) {
	        Blk_List.message_noreschange_topsoutheastcorner_recbuf[neighbour_blk][l] =
                  Blk_List.message_noreschange_bottomnorthwestcorner_sendbuf[i_blk][l];
             } /* endfor */
          } else if (Blk_List.Block[i_blk].infoBNW[0].dimen.i > 0 &&
		     Blk_List.Block[i_blk].infoBNW[0].dimen.j  > 0) {
             if (Blk_List.Block[neighbour_blk].nBSE != 1 ||
                 Blk_List.Block[neighbour_blk].infoBSE[0].cpu != i_cpu ||
                 Blk_List.Block[neighbour_blk].infoBSE[0].blknum != i_blk) return(5604);
             for (l = 0; l <= buffer_size-1; ++l) {
	        Blk_List.message_noreschange_bottomsoutheastcorner_recbuf[neighbour_blk][l] =
                  Blk_List.message_noreschange_bottomnorthwestcorner_sendbuf[i_blk][l];
             } /* endfor */
          } else if (Blk_List.Block[i_blk].infoBNW[0].dimen.j > 0 &&
		     Blk_List.Block[i_blk].infoBNW[0].dimen.k  > 0) {
             if (Blk_List.Block[neighbour_blk].nTSW != 1 ||
                 Blk_List.Block[neighbour_blk].infoTSW[0].cpu != i_cpu ||
                 Blk_List.Block[neighbour_blk].infoTSW[0].blknum != i_blk) return(5605);
             for (l = 0; l <= buffer_size-1; ++l) {
	        Blk_List.message_noreschange_topsouthwestcorner_recbuf[neighbour_blk][l] =
                  Blk_List.message_noreschange_bottomnorthwestcorner_sendbuf[i_blk][l];
             } /* endfor */
          } else if (Blk_List.Block[i_blk].infoBNW[0].dimen.k > 0 &&
		     Blk_List.Block[i_blk].infoBNW[0].dimen.i  > 0) {
             if (Blk_List.Block[neighbour_blk].nTNE != 1 ||
                 Blk_List.Block[neighbour_blk].infoTNE[0].cpu != i_cpu ||
                 Blk_List.Block[neighbour_blk].infoTNE[0].blknum != i_blk) return(5606);
             for (l = 0; l <= buffer_size-1; ++l) {
	        Blk_List.message_noreschange_topnortheastcorner_recbuf[neighbour_blk][l] =
                  Blk_List.message_noreschange_bottomnorthwestcorner_sendbuf[i_blk][l];
             } /* endfor */
          } else if (Blk_List.Block[i_blk].infoBNW[0].dimen.i > 0 ) {
             if (Blk_List.Block[neighbour_blk].nBNE != 1 ||
                 Blk_List.Block[neighbour_blk].infoBNE[0].cpu != i_cpu ||
                 Blk_List.Block[neighbour_blk].infoBNE[0].blknum != i_blk) return(5607);
             for (l = 0; l <= buffer_size-1; ++l) {
	        Blk_List.message_noreschange_bottomnortheastcorner_recbuf[neighbour_blk][l] =
                  Blk_List.message_noreschange_bottomnorthwestcorner_sendbuf[i_blk][l];
             } /* endfor */
          } else if (Blk_List.Block[i_blk].infoBNW[0].dimen.j > 0 ) {
             if (Blk_List.Block[neighbour_blk].nBSW != 1 ||
                 Blk_List.Block[neighbour_blk].infoBSW[0].cpu != i_cpu ||
                 Blk_List.Block[neighbour_blk].infoBSW[0].blknum != i_blk) return(5608);
             for (l = 0; l <= buffer_size-1; ++l) {
	        Blk_List.message_noreschange_bottomsouthwestcorner_recbuf[neighbour_blk][l] =
                  Blk_List.message_noreschange_bottomnorthwestcorner_sendbuf[i_blk][l];
             } /* endfor */
          } else if (Blk_List.Block[i_blk].infoBNW[0].dimen.k > 0 ) {
             if (Blk_List.Block[neighbour_blk].nTNW != 1 ||
                 Blk_List.Block[neighbour_blk].infoTNW[0].cpu != i_cpu ||
                 Blk_List.Block[neighbour_blk].infoTNW[0].blknum != i_blk) return(5609);
             for (l = 0; l <= buffer_size-1; ++l) {
	        Blk_List.message_noreschange_topnorthwestcorner_recbuf[neighbour_blk][l] =
                  Blk_List.message_noreschange_bottomnorthwestcorner_sendbuf[i_blk][l];
             } /* endfor */
          } else {
             if (Blk_List.Block[neighbour_blk].nBNW != 1 ||
                 Blk_List.Block[neighbour_blk].infoBNW[0].cpu != i_cpu ||
                 Blk_List.Block[neighbour_blk].infoBNW[0].blknum != i_blk) return(5610);
             for (l = 0; l <= buffer_size-1; ++l) {
	        Blk_List.message_noreschange_bottomnorthwestcorner_recbuf[neighbour_blk][l] =
                  Blk_List.message_noreschange_bottomnorthwestcorner_sendbuf[i_blk][l];
             } /* endfor */
          } /* endif */
       } /* endif */

       /////////////////////////////////////////////////////////////////////
       // Exchange messages at the Bottom North East corner of solution blocks. //
       /////////////////////////////////////////////////////////////////////
       if (Blk_List.Block[i_blk].used && 
           (Blk_List.Block[i_blk].nBNE == 1) &&
           (Blk_List.Block[i_blk].info.level == Blk_List.Block[i_blk].infoBNE[0].level)) {
	  neighbour_cpu = Blk_List.Block[i_blk].infoBNE[0].cpu;
	  neighbour_blk = Blk_List.Block[i_blk].infoBNE[0].blknum;
          if (neighbour_cpu != i_cpu) return(5700);
          if (!Blk_List.Block[neighbour_blk].used) return(5701); 
          buffer_size = sqr(Blk_List.Block[i_blk].infoBNE[0].dimen.ghost)*
	    Blk_List.Block[i_blk].infoBNE[0].dimen.ghost*
                        Number_of_Solution_Variables;
          buffer_size_neighbour = sqr(Blk_List.Block[neighbour_blk].info.dimen.ghost)*
	    Blk_List.Block[neighbour_blk].info.dimen.ghost*
                                  Number_of_Solution_Variables;
          if (buffer_size != buffer_size_neighbour) return(5702);
          if (Blk_List.Block[i_blk].infoBNE[0].dimen.i > 0 &&
	      Blk_List.Block[i_blk].infoBNE[0].dimen.j > 0 &&
	      Blk_List.Block[i_blk].infoBNE[0].dimen.k  > 0) {
             if (Blk_List.Block[neighbour_blk].nTSW != 1 ||
                 Blk_List.Block[neighbour_blk].infoTSW[0].cpu != i_cpu ||
                 Blk_List.Block[neighbour_blk].infoTSW[0].blknum != i_blk) return(5703);
             for (l = 0; l <= buffer_size-1; ++l) {
	        Blk_List.message_noreschange_topsouthwestcorner_recbuf[neighbour_blk][l] =
                  Blk_List.message_noreschange_bottomnortheastcorner_sendbuf[i_blk][l];
             } /* endfor */
          } else if (Blk_List.Block[i_blk].infoBNE[0].dimen.i > 0 &&
		     Blk_List.Block[i_blk].infoBNE[0].dimen.j  > 0) {
             if (Blk_List.Block[neighbour_blk].nBSW != 1 ||
                 Blk_List.Block[neighbour_blk].infoBSW[0].cpu != i_cpu ||
                 Blk_List.Block[neighbour_blk].infoBSW[0].blknum != i_blk) return(5704);
             for (l = 0; l <= buffer_size-1; ++l) {
	        Blk_List.message_noreschange_bottomsouthwestcorner_recbuf[neighbour_blk][l] =
                  Blk_List.message_noreschange_bottomnortheastcorner_sendbuf[i_blk][l];
             } /* endfor */
          } else if (Blk_List.Block[i_blk].infoBNE[0].dimen.j > 0 &&
		     Blk_List.Block[i_blk].infoBNE[0].dimen.k  > 0) {
             if (Blk_List.Block[neighbour_blk].nTSE != 1 ||
                 Blk_List.Block[neighbour_blk].infoTSE[0].cpu != i_cpu ||
                 Blk_List.Block[neighbour_blk].infoTSE[0].blknum != i_blk) return(5705);
             for (l = 0; l <= buffer_size-1; ++l) {
	        Blk_List.message_noreschange_topsoutheastcorner_recbuf[neighbour_blk][l] =
                  Blk_List.message_noreschange_bottomnortheastcorner_sendbuf[i_blk][l];
             } /* endfor */
          } else if (Blk_List.Block[i_blk].infoBNE[0].dimen.k > 0 &&
		     Blk_List.Block[i_blk].infoBNE[0].dimen.i  > 0) {
             if (Blk_List.Block[neighbour_blk].nTNW != 1 ||
                 Blk_List.Block[neighbour_blk].infoTNW[0].cpu != i_cpu ||
                 Blk_List.Block[neighbour_blk].infoTNW[0].blknum != i_blk) return(5706);
             for (l = 0; l <= buffer_size-1; ++l) {
	        Blk_List.message_noreschange_topnorthwestcorner_recbuf[neighbour_blk][l] =
                  Blk_List.message_noreschange_bottomnortheastcorner_sendbuf[i_blk][l];
             } /* endfor */
          } else if (Blk_List.Block[i_blk].infoBNE[0].dimen.i > 0 ) {
             if (Blk_List.Block[neighbour_blk].nBNW != 1 ||
                 Blk_List.Block[neighbour_blk].infoBNW[0].cpu != i_cpu ||
                 Blk_List.Block[neighbour_blk].infoBNW[0].blknum != i_blk) return(5707);
             for (l = 0; l <= buffer_size-1; ++l) {
	        Blk_List.message_noreschange_bottomnorthwestcorner_recbuf[neighbour_blk][l] =
                  Blk_List.message_noreschange_bottomnortheastcorner_sendbuf[i_blk][l];
             } /* endfor */
          } else if (Blk_List.Block[i_blk].infoBNE[0].dimen.j > 0 ) {
             if (Blk_List.Block[neighbour_blk].nBSE != 1 ||
                 Blk_List.Block[neighbour_blk].infoBSE[0].cpu != i_cpu ||
                 Blk_List.Block[neighbour_blk].infoBSE[0].blknum != i_blk) return(5708);
             for (l = 0; l <= buffer_size-1; ++l) {
	        Blk_List.message_noreschange_bottomsoutheastcorner_recbuf[neighbour_blk][l] =
                  Blk_List.message_noreschange_bottomnortheastcorner_sendbuf[i_blk][l];
             } /* endfor */
          } else if (Blk_List.Block[i_blk].infoBNE[0].dimen.k > 0 ) {
             if (Blk_List.Block[neighbour_blk].nTNE != 1 ||
                 Blk_List.Block[neighbour_blk].infoTNE[0].cpu != i_cpu ||
                 Blk_List.Block[neighbour_blk].infoTNE[0].blknum != i_blk) return(5709);
             for (l = 0; l <= buffer_size-1; ++l) {
	        Blk_List.message_noreschange_topnortheastcorner_recbuf[neighbour_blk][l] =
                  Blk_List.message_noreschange_bottomnortheastcorner_sendbuf[i_blk][l];
             } /* endfor */
          } else {
             if (Blk_List.Block[neighbour_blk].nBNE != 1 ||
                 Blk_List.Block[neighbour_blk].infoBNE[0].cpu != i_cpu ||
                 Blk_List.Block[neighbour_blk].infoBNE[0].blknum != i_blk) return(5710);
             for (l = 0; l <= buffer_size-1; ++l) {
	        Blk_List.message_noreschange_bottomnortheastcorner_recbuf[neighbour_blk][l] =
                  Blk_List.message_noreschange_bottomnortheastcorner_sendbuf[i_blk][l];
             } /* endfor */
          } /* endif */
       } /* endif */

       /////////////////////////////////////////////////////////////////////
       // Exchange messages at the Bottom South East corner of solution blocks. //
       /////////////////////////////////////////////////////////////////////
       if (Blk_List.Block[i_blk].used && 
           (Blk_List.Block[i_blk].nBSE == 1) &&
           (Blk_List.Block[i_blk].info.level == Blk_List.Block[i_blk].infoBSE[0].level)) {
	  neighbour_cpu = Blk_List.Block[i_blk].infoBSE[0].cpu;
	  neighbour_blk = Blk_List.Block[i_blk].infoBSE[0].blknum;
          if (neighbour_cpu != i_cpu) return(5900);
          if (!Blk_List.Block[neighbour_blk].used) return(5901); 
          buffer_size = sqr(Blk_List.Block[i_blk].infoBSE[0].dimen.ghost)*
	    Blk_List.Block[i_blk].infoBSE[0].dimen.ghost*
                        Number_of_Solution_Variables;
          buffer_size_neighbour = sqr(Blk_List.Block[neighbour_blk].info.dimen.ghost)*
	    Blk_List.Block[neighbour_blk].info.dimen.ghost*
                                  Number_of_Solution_Variables;
          if (buffer_size != buffer_size_neighbour) return(5902);
          if (Blk_List.Block[i_blk].infoBSE[0].dimen.i > 0 &&
	      Blk_List.Block[i_blk].infoBSE[0].dimen.j > 0 &&
	      Blk_List.Block[i_blk].infoBSE[0].dimen.k  > 0) {
             if (Blk_List.Block[neighbour_blk].nTNW != 1 ||
                 Blk_List.Block[neighbour_blk].infoTNW[0].cpu != i_cpu ||
                 Blk_List.Block[neighbour_blk].infoTNW[0].blknum != i_blk) return(5903);
             for (l = 0; l <= buffer_size-1; ++l) {
	        Blk_List.message_noreschange_topnorthwestcorner_recbuf[neighbour_blk][l] =
                  Blk_List.message_noreschange_bottomsoutheastcorner_sendbuf[i_blk][l];
             } /* endfor */
          } else if (Blk_List.Block[i_blk].infoBSE[0].dimen.i > 0 &&
		     Blk_List.Block[i_blk].infoBSE[0].dimen.j  > 0) {
             if (Blk_List.Block[neighbour_blk].nBNW != 1 ||
                 Blk_List.Block[neighbour_blk].infoBNW[0].cpu != i_cpu ||
                 Blk_List.Block[neighbour_blk].infoBNW[0].blknum != i_blk) return(5904);
             for (l = 0; l <= buffer_size-1; ++l) {
	        Blk_List.message_noreschange_bottomnorthwestcorner_recbuf[neighbour_blk][l] =
                  Blk_List.message_noreschange_bottomsoutheastcorner_sendbuf[i_blk][l];
             } /* endfor */
          } else if (Blk_List.Block[i_blk].infoBSE[0].dimen.j > 0 &&
		     Blk_List.Block[i_blk].infoBSE[0].dimen.k  > 0) {
             if (Blk_List.Block[neighbour_blk].nTNE != 1 ||
                 Blk_List.Block[neighbour_blk].infoTNE[0].cpu != i_cpu ||
                 Blk_List.Block[neighbour_blk].infoTNE[0].blknum != i_blk) return(5905);
             for (l = 0; l <= buffer_size-1; ++l) {
	        Blk_List.message_noreschange_topnortheastcorner_recbuf[neighbour_blk][l] =
                  Blk_List.message_noreschange_bottomsoutheastcorner_sendbuf[i_blk][l];
             } /* endfor */
          } else if (Blk_List.Block[i_blk].infoBSE[0].dimen.k > 0 &&
		     Blk_List.Block[i_blk].infoBSE[0].dimen.i  > 0) {
             if (Blk_List.Block[neighbour_blk].nTSW != 1 ||
                 Blk_List.Block[neighbour_blk].infoTSW[0].cpu != i_cpu ||
                 Blk_List.Block[neighbour_blk].infoTSW[0].blknum != i_blk) return(5906);
             for (l = 0; l <= buffer_size-1; ++l) {
	        Blk_List.message_noreschange_topsouthwestcorner_recbuf[neighbour_blk][l] =
                  Blk_List.message_noreschange_bottomsoutheastcorner_sendbuf[i_blk][l];
             } /* endfor */
          } else if (Blk_List.Block[i_blk].infoBSE[0].dimen.i > 0 ) {
             if (Blk_List.Block[neighbour_blk].nBSW != 1 ||
                 Blk_List.Block[neighbour_blk].infoBSW[0].cpu != i_cpu ||
                 Blk_List.Block[neighbour_blk].infoBSW[0].blknum != i_blk) return(5907);
             for (l = 0; l <= buffer_size-1; ++l) {
	        Blk_List.message_noreschange_bottomsouthwestcorner_recbuf[neighbour_blk][l] =
                  Blk_List.message_noreschange_bottomsoutheastcorner_sendbuf[i_blk][l];
             } /* endfor */
          } else if (Blk_List.Block[i_blk].infoBSE[0].dimen.j > 0 ) {
             if (Blk_List.Block[neighbour_blk].nBNE != 1 ||
                 Blk_List.Block[neighbour_blk].infoBNE[0].cpu != i_cpu ||
                 Blk_List.Block[neighbour_blk].infoBNE[0].blknum != i_blk) return(5908);
             for (l = 0; l <= buffer_size-1; ++l) {
	        Blk_List.message_noreschange_bottomnortheastcorner_recbuf[neighbour_blk][l] =
                  Blk_List.message_noreschange_bottomsoutheastcorner_sendbuf[i_blk][l];
             } /* endfor */
          } else if (Blk_List.Block[i_blk].infoBSE[0].dimen.k > 0 ) {
             if (Blk_List.Block[neighbour_blk].nTSE != 1 ||
                 Blk_List.Block[neighbour_blk].infoTSE[0].cpu != i_cpu ||
                 Blk_List.Block[neighbour_blk].infoTSE[0].blknum != i_blk) return(5909);
             for (l = 0; l <= buffer_size-1; ++l) {
	        Blk_List.message_noreschange_topsoutheastcorner_recbuf[neighbour_blk][l] =
                  Blk_List.message_noreschange_bottomsoutheastcorner_sendbuf[i_blk][l];
             } /* endfor */
          } else {
             if (Blk_List.Block[neighbour_blk].nBSE != 1 ||
                 Blk_List.Block[neighbour_blk].infoBSE[0].cpu != i_cpu ||
                 Blk_List.Block[neighbour_blk].infoBSE[0].blknum != i_blk) return(5910);
             for (l = 0; l <= buffer_size-1; ++l) {
	        Blk_List.message_noreschange_bottomsoutheastcorner_recbuf[neighbour_blk][l] =
                  Blk_List.message_noreschange_bottomsoutheastcorner_sendbuf[i_blk][l];
             } /* endfor */
          } /* endif */
       } /* endif */

       /////////////////////////////////////////////////////////////////////
       // Exchange messages at the Top North East corner of solution blocks. //
       /////////////////////////////////////////////////////////////////////
       if (Blk_List.Block[i_blk].used && 
           (Blk_List.Block[i_blk].nTNE == 1) &&
           (Blk_List.Block[i_blk].info.level == Blk_List.Block[i_blk].infoTNE[0].level)) {
	  neighbour_cpu = Blk_List.Block[i_blk].infoTNE[0].cpu;
	  neighbour_blk = Blk_List.Block[i_blk].infoTNE[0].blknum;
          if (neighbour_cpu != i_cpu) return(5300);
          if (!Blk_List.Block[neighbour_blk].used) return(5301); 
          buffer_size = sqr(Blk_List.Block[i_blk].infoTNE[0].dimen.ghost)*
	    Blk_List.Block[i_blk].infoTNE[0].dimen.ghost*
                        Number_of_Solution_Variables;
          buffer_size_neighbour = sqr(Blk_List.Block[neighbour_blk].info.dimen.ghost)*
	    Blk_List.Block[neighbour_blk].info.dimen.ghost*
                                  Number_of_Solution_Variables;
          if (buffer_size != buffer_size_neighbour) return(5302);
          if (Blk_List.Block[i_blk].infoTNE[0].dimen.i > 0 &&
	      Blk_List.Block[i_blk].infoTNE[0].dimen.j > 0 &&
	      Blk_List.Block[i_blk].infoTNE[0].dimen.k  > 0) {
             if (Blk_List.Block[neighbour_blk].nBSW != 1 ||
                 Blk_List.Block[neighbour_blk].infoBSW[0].cpu != i_cpu ||
                 Blk_List.Block[neighbour_blk].infoBSW[0].blknum != i_blk) return(5303);
             for (l = 0; l <= buffer_size-1; ++l) {
	        Blk_List.message_noreschange_bottomsouthwestcorner_recbuf[neighbour_blk][l] =
                  Blk_List.message_noreschange_topnortheastcorner_sendbuf[i_blk][l];
             } /* endfor */
          } else if (Blk_List.Block[i_blk].infoTNE[0].dimen.i > 0 &&
		     Blk_List.Block[i_blk].infoTNE[0].dimen.j  > 0) {
             if (Blk_List.Block[neighbour_blk].nTSW != 1 ||
                 Blk_List.Block[neighbour_blk].infoTSW[0].cpu != i_cpu ||
                 Blk_List.Block[neighbour_blk].infoTSW[0].blknum != i_blk) return(5304);
             for (l = 0; l <= buffer_size-1; ++l) {
	        Blk_List.message_noreschange_topsouthwestcorner_recbuf[neighbour_blk][l] =
                  Blk_List.message_noreschange_topnortheastcorner_sendbuf[i_blk][l];
             } /* endfor */
          } else if (Blk_List.Block[i_blk].infoTNE[0].dimen.j > 0 &&
		     Blk_List.Block[i_blk].infoTNE[0].dimen.k  > 0) {
             if (Blk_List.Block[neighbour_blk].nBSE != 1 ||
                 Blk_List.Block[neighbour_blk].infoBSE[0].cpu != i_cpu ||
                 Blk_List.Block[neighbour_blk].infoBSE[0].blknum != i_blk) return(5305);
             for (l = 0; l <= buffer_size-1; ++l) {
	        Blk_List.message_noreschange_bottomsoutheastcorner_recbuf[neighbour_blk][l] =
                  Blk_List.message_noreschange_topnortheastcorner_sendbuf[i_blk][l];
             } /* endfor */
          } else if (Blk_List.Block[i_blk].infoTNE[0].dimen.k > 0 &&
		     Blk_List.Block[i_blk].infoTNE[0].dimen.i  > 0) {
             if (Blk_List.Block[neighbour_blk].nBNW != 1 ||
                 Blk_List.Block[neighbour_blk].infoBNW[0].cpu != i_cpu ||
                 Blk_List.Block[neighbour_blk].infoBNW[0].blknum != i_blk) return(5306);
             for (l = 0; l <= buffer_size-1; ++l) {
	        Blk_List.message_noreschange_bottomnorthwestcorner_recbuf[neighbour_blk][l] =
                  Blk_List.message_noreschange_topnortheastcorner_sendbuf[i_blk][l];
             } /* endfor */
          } else if (Blk_List.Block[i_blk].infoTNE[0].dimen.i > 0 ) {
             if (Blk_List.Block[neighbour_blk].nTNW != 1 ||
                 Blk_List.Block[neighbour_blk].infoTNW[0].cpu != i_cpu ||
                 Blk_List.Block[neighbour_blk].infoTNW[0].blknum != i_blk) return(5307);
             for (l = 0; l <= buffer_size-1; ++l) {
	        Blk_List.message_noreschange_topnorthwestcorner_recbuf[neighbour_blk][l] =
                  Blk_List.message_noreschange_topnortheastcorner_sendbuf[i_blk][l];
             } /* endfor */
          } else if (Blk_List.Block[i_blk].infoTNE[0].dimen.j > 0 ) {
             if (Blk_List.Block[neighbour_blk].nTSE != 1 ||
                 Blk_List.Block[neighbour_blk].infoTSE[0].cpu != i_cpu ||
                 Blk_List.Block[neighbour_blk].infoTSE[0].blknum != i_blk) return(5308);
             for (l = 0; l <= buffer_size-1; ++l) {
	        Blk_List.message_noreschange_topsoutheastcorner_recbuf[neighbour_blk][l] =
                  Blk_List.message_noreschange_topnortheastcorner_sendbuf[i_blk][l];
             } /* endfor */
          } else if (Blk_List.Block[i_blk].infoTNE[0].dimen.k > 0 ) {
             if (Blk_List.Block[neighbour_blk].nBNE != 1 ||
                 Blk_List.Block[neighbour_blk].infoBNE[0].cpu != i_cpu ||
                 Blk_List.Block[neighbour_blk].infoBNE[0].blknum != i_blk) return(5309);
             for (l = 0; l <= buffer_size-1; ++l) {
	        Blk_List.message_noreschange_bottomnortheastcorner_recbuf[neighbour_blk][l] =
                  Blk_List.message_noreschange_topnortheastcorner_sendbuf[i_blk][l];
             } /* endfor */
          } else {
             if (Blk_List.Block[neighbour_blk].nTNE != 1 ||
                 Blk_List.Block[neighbour_blk].infoTNE[0].cpu != i_cpu ||
                 Blk_List.Block[neighbour_blk].infoTNE[0].blknum != i_blk) return(5310);
             for (l = 0; l <= buffer_size-1; ++l) {
	        Blk_List.message_noreschange_topnortheastcorner_recbuf[neighbour_blk][l] =
                  Blk_List.message_noreschange_topnortheastcorner_sendbuf[i_blk][l];
             } /* endfor */
          } /* endif */
       } /* endif */


       /////////////////////////////////////////////////////////////////////
       // Exchange messages at the Bottom South West corner of solution blocks. //
       /////////////////////////////////////////////////////////////////////
       if (Blk_List.Block[i_blk].used && 
           (Blk_List.Block[i_blk].nBSW == 1) &&
           (Blk_List.Block[i_blk].info.level == Blk_List.Block[i_blk].infoBSW[0].level)) {
	  neighbour_cpu = Blk_List.Block[i_blk].infoBSW[0].cpu;
	  neighbour_blk = Blk_List.Block[i_blk].infoBSW[0].blknum;
          if (neighbour_cpu != i_cpu) return(5800);
          if (!Blk_List.Block[neighbour_blk].used) return(5801); 
          buffer_size = sqr(Blk_List.Block[i_blk].infoBSW[0].dimen.ghost)*
	    Blk_List.Block[i_blk].infoBSW[0].dimen.ghost*
                        Number_of_Solution_Variables;
          buffer_size_neighbour = sqr(Blk_List.Block[neighbour_blk].info.dimen.ghost)*
	    Blk_List.Block[neighbour_blk].info.dimen.ghost*
                                  Number_of_Solution_Variables;
          if (buffer_size != buffer_size_neighbour) return(5802);
          if (Blk_List.Block[i_blk].infoBSW[0].dimen.i > 0 &&
	      Blk_List.Block[i_blk].infoBSW[0].dimen.j > 0 &&
	      Blk_List.Block[i_blk].infoBSW[0].dimen.k  > 0) {
             if (Blk_List.Block[neighbour_blk].nTNE != 1 ||
                 Blk_List.Block[neighbour_blk].infoTNE[0].cpu != i_cpu ||
                 Blk_List.Block[neighbour_blk].infoTNE[0].blknum != i_blk) return(5803);
             for (l = 0; l <= buffer_size-1; ++l) {
	        Blk_List.message_noreschange_topnortheastcorner_recbuf[neighbour_blk][l] =
                  Blk_List.message_noreschange_bottomsouthwestcorner_sendbuf[i_blk][l];
             } /* endfor */
          } else if (Blk_List.Block[i_blk].infoBSW[0].dimen.i > 0 &&
		     Blk_List.Block[i_blk].infoBSW[0].dimen.j  > 0) {
             if (Blk_List.Block[neighbour_blk].nBNE != 1 ||
                 Blk_List.Block[neighbour_blk].infoBNE[0].cpu != i_cpu ||
                 Blk_List.Block[neighbour_blk].infoBNE[0].blknum != i_blk) return(5804);
             for (l = 0; l <= buffer_size-1; ++l) {
	        Blk_List.message_noreschange_bottomnortheastcorner_recbuf[neighbour_blk][l] =
                  Blk_List.message_noreschange_bottomsouthwestcorner_sendbuf[i_blk][l];
             } /* endfor */
          } else if (Blk_List.Block[i_blk].infoBSW[0].dimen.j > 0 &&
		     Blk_List.Block[i_blk].infoBSW[0].dimen.k  > 0) {
             if (Blk_List.Block[neighbour_blk].nTNW != 1 ||
                 Blk_List.Block[neighbour_blk].infoTNW[0].cpu != i_cpu ||
                 Blk_List.Block[neighbour_blk].infoTNW[0].blknum != i_blk) return(5805);
             for (l = 0; l <= buffer_size-1; ++l) {
	        Blk_List.message_noreschange_topnorthwestcorner_recbuf[neighbour_blk][l] =
                  Blk_List.message_noreschange_bottomsouthwestcorner_sendbuf[i_blk][l];
             } /* endfor */
          } else if (Blk_List.Block[i_blk].infoBSW[0].dimen.k > 0 &&
		     Blk_List.Block[i_blk].infoBSW[0].dimen.i  > 0) {
             if (Blk_List.Block[neighbour_blk].nTSE != 1 ||
                 Blk_List.Block[neighbour_blk].infoTSE[0].cpu != i_cpu ||
                 Blk_List.Block[neighbour_blk].infoTSE[0].blknum != i_blk) return(5806);
             for (l = 0; l <= buffer_size-1; ++l) {
	        Blk_List.message_noreschange_topsoutheastcorner_recbuf[neighbour_blk][l] =
                  Blk_List.message_noreschange_bottomsouthwestcorner_sendbuf[i_blk][l];
             } /* endfor */
          } else if (Blk_List.Block[i_blk].infoBSW[0].dimen.i > 0 ) {
             if (Blk_List.Block[neighbour_blk].nBSE != 1 ||
                 Blk_List.Block[neighbour_blk].infoBSE[0].cpu != i_cpu ||
                 Blk_List.Block[neighbour_blk].infoBSE[0].blknum != i_blk) return(5807);
             for (l = 0; l <= buffer_size-1; ++l) {
	        Blk_List.message_noreschange_bottomsoutheastcorner_recbuf[neighbour_blk][l] =
                  Blk_List.message_noreschange_bottomsouthwestcorner_sendbuf[i_blk][l];
             } /* endfor */
          } else if (Blk_List.Block[i_blk].infoBSW[0].dimen.j > 0 ) {
             if (Blk_List.Block[neighbour_blk].nBNW != 1 ||
                 Blk_List.Block[neighbour_blk].infoBNW[0].cpu != i_cpu ||
                 Blk_List.Block[neighbour_blk].infoBNW[0].blknum != i_blk) return(5808);
             for (l = 0; l <= buffer_size-1; ++l) {
	        Blk_List.message_noreschange_bottomnorthwestcorner_recbuf[neighbour_blk][l] =
                  Blk_List.message_noreschange_bottomsouthwestcorner_sendbuf[i_blk][l];
             } /* endfor */
          } else if (Blk_List.Block[i_blk].infoBSW[0].dimen.k > 0 ) {
             if (Blk_List.Block[neighbour_blk].nTSW != 1 ||
                 Blk_List.Block[neighbour_blk].infoTSW[0].cpu != i_cpu ||
                 Blk_List.Block[neighbour_blk].infoTSW[0].blknum != i_blk) return(5809);
             for (l = 0; l <= buffer_size-1; ++l) {
	        Blk_List.message_noreschange_topsouthwestcorner_recbuf[neighbour_blk][l] =
                  Blk_List.message_noreschange_bottomsouthwestcorner_sendbuf[i_blk][l];
             } /* endfor */
          } else {
             if (Blk_List.Block[neighbour_blk].nBSW != 1 ||
                 Blk_List.Block[neighbour_blk].infoBSW[0].cpu != i_cpu ||
                 Blk_List.Block[neighbour_blk].infoBSW[0].blknum != i_blk) return(5810);
             for (l = 0; l <= buffer_size-1; ++l) {
	        Blk_List.message_noreschange_bottomsouthwestcorner_recbuf[neighbour_blk][l] =
                  Blk_List.message_noreschange_bottomsouthwestcorner_sendbuf[i_blk][l];
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
int  AdaptiveBlock3D_List::Exchange_Messages_ResChange_FineToCoarse(AdaptiveBlock3D_List &Blk_List,
                                             const int Number_of_Solution_Variables) {

  cout<<"\n\nError Exchange_Messages_ResChange_FineToCoarse() does not exist for 3D\n";
  return (2);


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
int  AdaptiveBlock3D_List::Exchange_Messages_ResChange_CoarseToFine(AdaptiveBlock3D_List &Blk_List,
                                             const int Number_of_Solution_Variables) {

  cout<<"\n\nError Exchange_Messages_ResChange_CoarseToFine() does not exist for 3D\n";
  return (2);


}
