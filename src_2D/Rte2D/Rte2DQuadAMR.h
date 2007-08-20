#ifndef _RTE2D_AMR_INCLUDED
#define _RTE2D_AMR_INCLUDED

/****************************************************************/
/****** AMR REQUIRED SPECIALIZATIONS & FUNCTIONS ****************/
/****************************************************************/  
#include "Rte2DQuad.h"
#include "../AMR/AMR.h"


/*************************************************************
 * AdaptiveBlock2D -- Templated message passing subroutines. *
 *************************************************************/

/********************************************************
 * Routine: Load_Send_Message_Buffers_NoResChange       *
 *                                                      *
 * Loads the send message passing buffers for sending   *
 * solution information to neighbouring 2D              *
 * quadrilateral multi-block solution blocks with no    *
 * mesh resolution changes (works on entire 1D array of *
 * 2D solution blocks).                                 *
 *                                                      *
 ********************************************************/
int Load_Send_Message_Buffers_NoResChange(Rte2D_Quad_Block *Soln_ptr,
                                          AdaptiveBlock2D_List &Soln_Block_List,
                                          const int Number_of_Solution_Variables) {

    int i_blk, buffer_size_neighbour, 
        i_min, i_max, i_inc, i, 
        j_min, j_max, j_inc, j, 
        k, l;
    int i_ref, j_ref;
    Vector2D x_ref;

    /* Check to see if the number of solution variables specified
       in the calling argument is correct. */

    if (Number_of_Solution_Variables != Soln_ptr[0].NumVar()) {    
       return(1);
    } /* endif */
   
    /* Load the send buffers of each solution block. */

    for ( i_blk = 0 ; i_blk <= Soln_Block_List.Nblk-1 ; ++i_blk ) {

       ////////////////////////////////////////////////////////////////////
       // Load send buffer for North face neighbours of solution blocks. //
       ////////////////////////////////////////////////////////////////////
       if (Soln_Block_List.Block[i_blk].used && 
           (Soln_Block_List.Block[i_blk].nN == 1) &&
           (Soln_Block_List.Block[i_blk].info.level == 
            Soln_Block_List.Block[i_blk].infoN[0].level)) {
	       buffer_size_neighbour = abs(Soln_Block_List.Block[i_blk].infoN[0].dimen.i)*
                                       Soln_ptr[i_blk].NumVar();
          l = -1;
          // Load ghost cell solution variable information as required.
             if (Soln_Block_List.Block[i_blk].infoN[0].dimen.i > 0 &&
                 Soln_Block_List.Block[i_blk].infoN[0].dimen.j > 0) {
	        i_min = Soln_ptr[i_blk].ICl;
                i_max = Soln_ptr[i_blk].ICu;
	        i_inc = 1;
	        j_min = Soln_ptr[i_blk].JCu-Soln_Block_List.Block[i_blk].infoN[0].dimen.ghost+1;
	        j_max = Soln_ptr[i_blk].JCu;
	        j_inc = 1;
             } else if (Soln_Block_List.Block[i_blk].infoN[0].dimen.j > 0) {
	        i_min = Soln_ptr[i_blk].ICu;
	        i_max = Soln_ptr[i_blk].ICl;
	        i_inc = -1;
	        j_min = Soln_ptr[i_blk].JCu-Soln_Block_List.Block[i_blk].infoN[0].dimen.ghost+1;
	        j_max = Soln_ptr[i_blk].JCu;
	        j_inc = 1;
             } else if (Soln_Block_List.Block[i_blk].infoN[0].dimen.i > 0) {
	        i_min = Soln_ptr[i_blk].ICl;
	        i_max = Soln_ptr[i_blk].ICu;
	        i_inc = 1;
	        j_min = Soln_ptr[i_blk].JCu;
	        j_max = Soln_ptr[i_blk].JCu-Soln_Block_List.Block[i_blk].infoN[0].dimen.ghost+1;
	        j_inc = -1;
             } else {
	        i_min = Soln_ptr[i_blk].ICu;
	        i_max = Soln_ptr[i_blk].ICl;
	        i_inc = -1;
	        j_min = Soln_ptr[i_blk].JCu;
 	        j_max = Soln_ptr[i_blk].JCu-Soln_Block_List.Block[i_blk].infoN[0].dimen.ghost+1;
	        j_inc = -1;
             } /* endif */
	       j_max = max(j_min, j_max); 
	       i = Soln_ptr[i_blk].LoadSendBuffer_BC(Soln_Block_List.message_noreschange_northface_sendbuf[i_blk],
						     l,buffer_size_neighbour,
						     i_min,i_max,i_inc,
						     j_max,j_max,j_inc);
	     if (i != 0) return(2200);
       } /* endif */

       ////////////////////////////////////////////////////////////////////
       // Load send buffer for South face neighbours of solution blocks. //
       ////////////////////////////////////////////////////////////////////
       if (Soln_Block_List.Block[i_blk].used && 
           (Soln_Block_List.Block[i_blk].nS == 1) &&
           (Soln_Block_List.Block[i_blk].info.level == 
            Soln_Block_List.Block[i_blk].infoS[0].level)) {
	      buffer_size_neighbour = abs(Soln_Block_List.Block[i_blk].infoS[0].dimen.i)*
                                      Soln_ptr[i_blk].NumVar();
          l = -1;
          // Load ghost cell solution variable information as required.
             if (Soln_Block_List.Block[i_blk].infoS[0].dimen.i > 0 &&
                 Soln_Block_List.Block[i_blk].infoS[0].dimen.j > 0) {
	        i_min = Soln_ptr[i_blk].ICl;
	        i_max = Soln_ptr[i_blk].ICu;
	        i_inc = 1;
	        j_min = Soln_ptr[i_blk].JCl;
	        j_max = Soln_ptr[i_blk].JCl+Soln_Block_List.Block[i_blk].infoS[0].dimen.ghost-1;
	        j_inc = 1;
             } else if (Soln_Block_List.Block[i_blk].infoS[0].dimen.j > 0) {
	        i_min = Soln_ptr[i_blk].ICu;
	        i_max = Soln_ptr[i_blk].ICl;
	        i_inc = -1;
	        j_min = Soln_ptr[i_blk].JCl;
	        j_max = Soln_ptr[i_blk].JCl+Soln_Block_List.Block[i_blk].infoS[0].dimen.ghost-1;
	        j_inc = 1;
             } else if (Soln_Block_List.Block[i_blk].infoS[0].dimen.i > 0) {
	        i_min = Soln_ptr[i_blk].ICl;
	        i_max = Soln_ptr[i_blk].ICu;
	        i_inc = 1;
	        j_min = Soln_ptr[i_blk].JCl+Soln_Block_List.Block[i_blk].infoS[0].dimen.ghost-1;
	        j_max = Soln_ptr[i_blk].JCl;
	        j_inc = -1;
             } else {
	        i_min = Soln_ptr[i_blk].ICu;
	        i_max = Soln_ptr[i_blk].ICl;
	        i_inc = -1;
	        j_min = Soln_ptr[i_blk].JCl+Soln_Block_List.Block[i_blk].infoS[0].dimen.ghost-1;
	        j_max = Soln_ptr[i_blk].JCl;
	        j_inc = -1;
             } /* endif */
	        j_min = min(j_min, j_max);
                i = Soln_ptr[i_blk].LoadSendBuffer_BC(Soln_Block_List.message_noreschange_southface_sendbuf[i_blk],
						      l,buffer_size_neighbour,
						      i_min,i_max,i_inc,
						      j_min,j_min,j_inc);
	     if (i != 0) return(2100);
       } /* endif */

       ///////////////////////////////////////////////////////////////////
       // Load send buffer for East face neighbours of solution blocks. //
       ///////////////////////////////////////////////////////////////////
       if (Soln_Block_List.Block[i_blk].used && 
           (Soln_Block_List.Block[i_blk].nE == 1) &&
           (Soln_Block_List.Block[i_blk].info.level == 
            Soln_Block_List.Block[i_blk].infoE[0].level)) {
	      buffer_size_neighbour = abs(Soln_Block_List.Block[i_blk].infoE[0].dimen.j)*
                                      Soln_ptr[i_blk].NumVar();
          l = -1;
          // Load ghost cell solution variable information as required.
             if (Soln_Block_List.Block[i_blk].infoE[0].dimen.i > 0 &&
                 Soln_Block_List.Block[i_blk].infoE[0].dimen.j > 0) {
	        i_min = Soln_ptr[i_blk].ICu-Soln_Block_List.Block[i_blk].infoE[0].dimen.ghost+1;
	        i_max = Soln_ptr[i_blk].ICu;
	        i_inc = 1;
	        j_min = Soln_ptr[i_blk].JCl;
	        j_max = Soln_ptr[i_blk].JCu;
	        j_inc = 1;
             } else if (Soln_Block_List.Block[i_blk].infoE[0].dimen.j > 0) {
	        i_min = Soln_ptr[i_blk].ICu;
	        i_max = Soln_ptr[i_blk].ICu-Soln_Block_List.Block[i_blk].infoE[0].dimen.ghost+1;
	        i_inc = -1;
	        j_min = Soln_ptr[i_blk].JCl;
	        j_max = Soln_ptr[i_blk].JCu;
	        j_inc = 1;
             } else if (Soln_Block_List.Block[i_blk].infoE[0].dimen.i > 0) {
	        i_min = Soln_ptr[i_blk].ICu-Soln_Block_List.Block[i_blk].infoE[0].dimen.ghost+1;
	        i_max = Soln_ptr[i_blk].ICu;
	        i_inc = 1;
	        j_min = Soln_ptr[i_blk].JCu;
	        j_max = Soln_ptr[i_blk].JCl;
	        j_inc = -1;
             } else {
	        i_min = Soln_ptr[i_blk].ICu;
	        i_max = Soln_ptr[i_blk].ICu-Soln_Block_List.Block[i_blk].infoE[0].dimen.ghost+1;
	        i_inc = -1;
	        j_min = Soln_ptr[i_blk].JCu;
	        j_max = Soln_ptr[i_blk].JCl;
	        j_inc = -1;
             } /* endif */
	       i_max = max(i_min, i_max);
	       i = Soln_ptr[i_blk].LoadSendBuffer_BC(Soln_Block_List.message_noreschange_eastface_sendbuf[i_blk],
						     l,buffer_size_neighbour,
						     i_max,i_max,i_inc,
						     j_min,j_max,j_inc);	       
	     if (i != 0) return(1200);
       } /* endif */

       ///////////////////////////////////////////////////////////////////
       // Load send buffer for West face neighbours of solution blocks. //
       ///////////////////////////////////////////////////////////////////
       if (Soln_Block_List.Block[i_blk].used && 
           (Soln_Block_List.Block[i_blk].nW == 1) &&
           (Soln_Block_List.Block[i_blk].info.level == 
            Soln_Block_List.Block[i_blk].infoW[0].level)) {
	      buffer_size_neighbour = abs(Soln_Block_List.Block[i_blk].infoW[0].dimen.j)*
                                      Soln_ptr[i_blk].NumVar();
          l = -1;
          // Load ghost cell solution variable information as required.
             if (Soln_Block_List.Block[i_blk].infoW[0].dimen.i > 0 &&
                 Soln_Block_List.Block[i_blk].infoW[0].dimen.j > 0) {
	        i_min = Soln_ptr[i_blk].ICl;
	        i_max = Soln_ptr[i_blk].ICl+Soln_Block_List.Block[i_blk].infoW[0].dimen.ghost-1;
	        i_inc = 1;
	        j_min = Soln_ptr[i_blk].JCl;
	        j_max = Soln_ptr[i_blk].JCu;
	        j_inc = 1;
             } else if (Soln_Block_List.Block[i_blk].infoW[0].dimen.j > 0) {
	        i_min = Soln_ptr[i_blk].ICl+Soln_Block_List.Block[i_blk].infoW[0].dimen.ghost-1;
	        i_max = Soln_ptr[i_blk].ICl;
	        i_inc = -1;
	        j_min = Soln_ptr[i_blk].JCl;
	        j_max = Soln_ptr[i_blk].JCu;
	        j_inc = 1;
             } else if (Soln_Block_List.Block[i_blk].infoW[0].dimen.i > 0) {
	        i_min = Soln_ptr[i_blk].ICl;
	        i_max = Soln_ptr[i_blk].ICl+Soln_Block_List.Block[i_blk].infoW[0].dimen.ghost-1;
 	        i_inc = 1;
	        j_min = Soln_ptr[i_blk].JCu;
	        j_max = Soln_ptr[i_blk].JCl;
	        j_inc = -1;
             } else {
	        i_min = Soln_ptr[i_blk].ICl+Soln_Block_List.Block[i_blk].infoW[0].dimen.ghost-1;
	        i_max = Soln_ptr[i_blk].ICl;
	        i_inc = -1;
	        j_min = Soln_ptr[i_blk].JCu;
	        j_max = Soln_ptr[i_blk].JCl;
	        j_inc = -1;
             } /* endif */
	       i_min = min(i_min, i_max);
	       i = Soln_ptr[i_blk].LoadSendBuffer_BC(Soln_Block_List.message_noreschange_westface_sendbuf[i_blk],
						     l,buffer_size_neighbour,
						     i_min,i_min,i_inc,
						     j_min,j_max,j_inc);
	     if (i != 0) return(1100);
       } /* endif */


    } /* endfor */

    /* Loading of send buffers complete.  Return zero value. */

    return(0);

}

/***************************************************************
 * Routine: Load_Send_Message_Buffers_ResChange_FineToCoarse   *
 *                                                             *
 * Loads the send message passing buffers for sending          *
 * solution information from more refined (higher mesh         *
 * resolution) solution blocks to more coarse neighbouring 2D  *
 * quadrilateral multi-block solution blocks with lower mesh   *
 * resolution (works on entire 1D array of 2D solution         *
 * blocks).                                                    *
 *                                                             *
 ***************************************************************/
int Load_Send_Message_Buffers_ResChange_FineToCoarse(Rte2D_Quad_Block *Soln_ptr,
                                                     AdaptiveBlock2D_List &Soln_Block_List,
                                                     const int Number_of_Solution_Variables) {

    int i_blk, buffer_size_neighbour, 
        i_min, i_max, i_inc, i, 
        j_min, j_max, j_inc, j, 
        k, l;
    int i_ref, j_ref;
    Vector2D x_ref;

    /* Check to see if the number of solution variables specified
       in the calling argument is correct. */

    if (Number_of_Solution_Variables != Soln_ptr[0].NumVar()) {
       return(1);
    } /* endif */

    /* Load the send buffers of each solution block. */

    for ( i_blk = 0 ; i_blk <= Soln_Block_List.Nblk-1 ; ++i_blk ) {

       //////////////////////////////////////////////////////////////////
       // Load send buffers of solution blocks with coarser North face //
       // neighbour at lower mesh resolution.                          //
       //////////////////////////////////////////////////////////////////
       if (Soln_Block_List.Block[i_blk].used && 
           (Soln_Block_List.Block[i_blk].nN == 1) &&
           (Soln_Block_List.Block[i_blk].info.level > 
            Soln_Block_List.Block[i_blk].infoN[0].level)) {
                buffer_size_neighbour = (abs(Soln_Block_List.Block[i_blk].infoN[0].dimen.i)/2)*
                                        Soln_ptr[i_blk].NumVar();
          l = -1;
          // Load ghost cell solution variable information as required.
             if (Soln_Block_List.Block[i_blk].infoN[0].dimen.i > 0 &&
                 Soln_Block_List.Block[i_blk].infoN[0].dimen.j > 0) {
	        i_min = Soln_ptr[i_blk].ICl;
                i_max = Soln_ptr[i_blk].ICu;
	        i_inc = 2;
	        j_min = Soln_ptr[i_blk].JCu-2*Soln_Block_List.Block[i_blk].infoN[0].dimen.ghost+1;
	        j_max = Soln_ptr[i_blk].JCu;
	        j_inc = 2;
             } else if (Soln_Block_List.Block[i_blk].infoN[0].dimen.j > 0) {
	        i_min = Soln_ptr[i_blk].ICu;
	        i_max = Soln_ptr[i_blk].ICl;
	        i_inc = -2;
	        j_min = Soln_ptr[i_blk].JCu-2*Soln_Block_List.Block[i_blk].infoN[0].dimen.ghost+1;
	        j_max = Soln_ptr[i_blk].JCu;
	        j_inc = 2;
             } else if (Soln_Block_List.Block[i_blk].infoN[0].dimen.i > 0) {
	        i_min = Soln_ptr[i_blk].ICl;
	        i_max = Soln_ptr[i_blk].ICu;
	        i_inc = 2;
	        j_min = Soln_ptr[i_blk].JCu;
	        j_max = Soln_ptr[i_blk].JCu-2*Soln_Block_List.Block[i_blk].infoN[0].dimen.ghost+1;
	        j_inc = -2;
             } else {
	        i_min = Soln_ptr[i_blk].ICu;
	        i_max = Soln_ptr[i_blk].ICl;
	        i_inc = -2;
	        j_min = Soln_ptr[i_blk].JCu;
 	        j_max = Soln_ptr[i_blk].JCu-2*Soln_Block_List.Block[i_blk].infoN[0].dimen.ghost+1;
	        j_inc = -2;
             } /* endif */
	        j_max = max(j_min, j_max);
                i = Soln_ptr[i_blk].LoadSendBuffer_BC_F2C(Soln_Block_List.message_reschange_northface_sendbuf[i_blk][0],
							  l,buffer_size_neighbour,
							  i_min,i_max,i_inc,
							  j_max,j_max,j_inc);
	     if (i != 0) return(2200);
       } /* endif */

       //////////////////////////////////////////////////////////////////
       // Load send buffers of solution blocks with coarser South face //
       // neighbour at lower mesh resolution.                          //
       //////////////////////////////////////////////////////////////////
       if (Soln_Block_List.Block[i_blk].used && 
           (Soln_Block_List.Block[i_blk].nS == 1) &&
           (Soln_Block_List.Block[i_blk].info.level > 
            Soln_Block_List.Block[i_blk].infoS[0].level)) {
                buffer_size_neighbour = (abs(Soln_Block_List.Block[i_blk].infoS[0].dimen.i)/2)*
                                        Soln_ptr[i_blk].NumVar();
          l = -1;
          // Load ghost cell solution variable information as required.
             if (Soln_Block_List.Block[i_blk].infoS[0].dimen.i > 0 &&
                 Soln_Block_List.Block[i_blk].infoS[0].dimen.j > 0) {
	        i_min = Soln_ptr[i_blk].ICl;
                i_max = Soln_ptr[i_blk].ICu;
	        i_inc = 2;
	        j_min = Soln_ptr[i_blk].JCl;
	        j_max = Soln_ptr[i_blk].JCl+2*Soln_Block_List.Block[i_blk].infoS[0].dimen.ghost-1;
	        j_inc = 2;
             } else if (Soln_Block_List.Block[i_blk].infoS[0].dimen.j > 0) {
	        i_min = Soln_ptr[i_blk].ICu;
	        i_max = Soln_ptr[i_blk].ICl;
	        i_inc = -2;
	        j_min = Soln_ptr[i_blk].JCl;
	        j_max = Soln_ptr[i_blk].JCl+2*Soln_Block_List.Block[i_blk].infoS[0].dimen.ghost-1;
	        j_inc = 2;
             } else if (Soln_Block_List.Block[i_blk].infoS[0].dimen.i > 0) {
	        i_min = Soln_ptr[i_blk].ICl;
	        i_max = Soln_ptr[i_blk].ICu;
	        i_inc = 2;
	        j_min = Soln_ptr[i_blk].JCl+2*Soln_Block_List.Block[i_blk].infoS[0].dimen.ghost-1;
	        j_max = Soln_ptr[i_blk].JCl;
	        j_inc = -2;
             } else {
	        i_min = Soln_ptr[i_blk].ICu;
	        i_max = Soln_ptr[i_blk].ICl;
	        i_inc = -2;
	        j_min = Soln_ptr[i_blk].JCl+2*Soln_Block_List.Block[i_blk].infoS[0].dimen.ghost-1;
	        j_max = Soln_ptr[i_blk].JCl;
	        j_inc = -2;
             } /* endif */
	        j_min = min(j_min, j_max);
                i = Soln_ptr[i_blk].LoadSendBuffer_BC_F2C(Soln_Block_List.message_reschange_southface_sendbuf[i_blk][0],
                                                            l,buffer_size_neighbour,
                                                            i_min,i_max,i_inc,
                                                            j_min,j_min,j_inc);
	     if (i != 0) return(2100);
       } /* endif */

       /////////////////////////////////////////////////////////////////
       // Load send buffers of solution blocks with coarser East face //
       // neighbour at lower mesh resolution.                         //
       /////////////////////////////////////////////////////////////////
       if (Soln_Block_List.Block[i_blk].used && 
           (Soln_Block_List.Block[i_blk].nE == 1) &&
           (Soln_Block_List.Block[i_blk].info.level > 
            Soln_Block_List.Block[i_blk].infoE[0].level)) {
                buffer_size_neighbour = (abs(Soln_Block_List.Block[i_blk].infoE[0].dimen.j)/2)*
                                        Soln_ptr[i_blk].NumVar();
          l = -1;
          // Load ghost cell solution variable information as required.
             if (Soln_Block_List.Block[i_blk].infoE[0].dimen.i > 0 &&
                 Soln_Block_List.Block[i_blk].infoE[0].dimen.j > 0) {
	        i_min = Soln_ptr[i_blk].ICu-2*Soln_Block_List.Block[i_blk].infoE[0].dimen.ghost+1;
	        i_max = Soln_ptr[i_blk].ICu;
	        i_inc = 2;
	        j_min = Soln_ptr[i_blk].JCl;
	        j_max = Soln_ptr[i_blk].JCu;
	        j_inc = 2;
             } else if (Soln_Block_List.Block[i_blk].infoE[0].dimen.j > 0) {
	        i_min = Soln_ptr[i_blk].ICu;
	        i_max = Soln_ptr[i_blk].ICu-2*Soln_Block_List.Block[i_blk].infoE[0].dimen.ghost+1;
	        i_inc = -2;
	        j_min = Soln_ptr[i_blk].JCl;
	        j_max = Soln_ptr[i_blk].JCu;
	        j_inc = 2;
             } else if (Soln_Block_List.Block[i_blk].infoE[0].dimen.i > 0) {
	        i_min = Soln_ptr[i_blk].ICu-2*Soln_Block_List.Block[i_blk].infoE[0].dimen.ghost+1;
	        i_max = Soln_ptr[i_blk].ICu;
	        i_inc = 2;
	        j_min = Soln_ptr[i_blk].JCu;
	        j_max = Soln_ptr[i_blk].JCl;
	        j_inc = -2;
             } else {
	        i_min = Soln_ptr[i_blk].ICu;
	        i_max = Soln_ptr[i_blk].ICu-2*Soln_Block_List.Block[i_blk].infoE[0].dimen.ghost+1;
	        i_inc = -2;
	        j_min = Soln_ptr[i_blk].JCu;
	        j_max = Soln_ptr[i_blk].JCl;
	        j_inc = -2;
             } /* endif */
	        i_max = max(i_min, i_max);
                i = Soln_ptr[i_blk].LoadSendBuffer_BC_F2C(Soln_Block_List.message_reschange_eastface_sendbuf[i_blk][0],
                                                            l,buffer_size_neighbour,
                                                            i_max,i_max,i_inc,
                                                            j_min,j_max,j_inc);
	     if (i != 0) return(1200);
       } /* endif */

       /////////////////////////////////////////////////////////////////
       // Load send buffers of solution blocks with coarser West face //
       // neighbour at lower mesh resolution.                         //
       /////////////////////////////////////////////////////////////////
       if (Soln_Block_List.Block[i_blk].used && 
           (Soln_Block_List.Block[i_blk].nW == 1) &&
           (Soln_Block_List.Block[i_blk].info.level > 
            Soln_Block_List.Block[i_blk].infoW[0].level)) {
                buffer_size_neighbour = (abs(Soln_Block_List.Block[i_blk].infoW[0].dimen.j)/2)*
                                        Soln_ptr[i_blk].NumVar();
          l = -1;
          // Load ghost cell solution variable information as required.
             if (Soln_Block_List.Block[i_blk].infoW[0].dimen.i > 0 &&
                 Soln_Block_List.Block[i_blk].infoW[0].dimen.j > 0) {
	        i_min = Soln_ptr[i_blk].ICl;
	        i_max = Soln_ptr[i_blk].ICl+2*Soln_Block_List.Block[i_blk].infoW[0].dimen.ghost-1;
	        i_inc = 2;
	        j_min = Soln_ptr[i_blk].JCl;
	        j_max = Soln_ptr[i_blk].JCu;
	        j_inc = 2;
             } else if (Soln_Block_List.Block[i_blk].infoW[0].dimen.j > 0) {
	        i_min = Soln_ptr[i_blk].ICl+2*Soln_Block_List.Block[i_blk].infoW[0].dimen.ghost-1;
	        i_max = Soln_ptr[i_blk].ICl;
	        i_inc = -2;
	        j_min = Soln_ptr[i_blk].JCl;
	        j_max = Soln_ptr[i_blk].JCu;
	        j_inc = 2;
             } else if (Soln_Block_List.Block[i_blk].infoW[0].dimen.i > 0) {
	        i_min = Soln_ptr[i_blk].ICl;
	        i_max = Soln_ptr[i_blk].ICl+2*Soln_Block_List.Block[i_blk].infoW[0].dimen.ghost-1;
	        i_inc = 2;
	        j_min = Soln_ptr[i_blk].JCu;
	        j_max = Soln_ptr[i_blk].JCl;
	        j_inc = -2;
             } else {
	        i_min = Soln_ptr[i_blk].ICl+2*Soln_Block_List.Block[i_blk].infoW[0].dimen.ghost-1;
	        i_max = Soln_ptr[i_blk].ICl;
	        i_inc = -2;
	        j_min = Soln_ptr[i_blk].JCu;
	        j_max = Soln_ptr[i_blk].JCl;
	        j_inc = -2;
             } /* endif */
	        i_min = min(i_min, i_max);
                i = Soln_ptr[i_blk].LoadSendBuffer_BC_F2C(Soln_Block_List.message_reschange_westface_sendbuf[i_blk][0],
                                                            l,buffer_size_neighbour,
                                                            i_min,i_min,i_inc,
                                                            j_min,j_max,j_inc);
	     if (i != 0) return(1100);
       } /* endif */


    }  /* endfor */

    /* Loading of send buffers complete.  Return zero value. */

    return(0);

}

/***************************************************************
 * Routine: Load_Send_Message_Buffers_ResChange_CoarseToFine   *
 *                                                             *
 * Loads the send message passing buffers for sending          *
 * solution information from more coarse (lower mesh           *
 * resolution) solution blocks to more refined neighbouring 2D *
 * quadrilateral multi-block solution blocks with higher mesh  *
 * resolution (works on entire 1D array of 2D solution         *
 * blocks).                                                    *
 *                                                             *
 ***************************************************************/
int Load_Send_Message_Buffers_ResChange_CoarseToFine(Rte2D_Quad_Block *Soln_ptr,
                                                     AdaptiveBlock2D_List &Soln_Block_List,
                                                     const int Number_of_Solution_Variables) {

    int i_blk, buffer_size_neighbour, 
        i_min0, i_max0, i_min1, i_max1, i_inc, i, 
        j_min0, j_max0, j_min1, j_max1, j_inc, j, 
        k, l0, l1;
    int i_ref0, j_ref0, i_ref1, j_ref1;
    Vector2D x_ref;

    /* Check to see if the number of solution variables specified
       in the calling argument is correct. */

    if (Number_of_Solution_Variables != Soln_ptr[0].NumVar()) {
       return(1);
    } /* endif */

    /* Load the send buffers of each solution block. */

    for ( i_blk = 0 ; i_blk <= Soln_Block_List.Nblk-1 ; ++i_blk ) {

       //////////////////////////////////////////////////////////////////
       // Load send buffers of solution blocks with refined North face //
       // neighbour at higher mesh resolution.                         //
       //////////////////////////////////////////////////////////////////
       if (Soln_Block_List.Block[i_blk].used && 
           (Soln_Block_List.Block[i_blk].nN == 2) &&
           (Soln_Block_List.Block[i_blk].info.level < 
            Soln_Block_List.Block[i_blk].infoN[0].level)) {
	      buffer_size_neighbour = (abs(Soln_Block_List.Block[i_blk].infoN[0].dimen.i)+
                                      Soln_Block_List.Block[i_blk].infoN[0].dimen.ghost)*
                                      Soln_ptr[i_blk].NumVar();
          l0 = -1;
          l1 = -1;
          // Load ghost cell solution variable information as required.
             if (Soln_Block_List.Block[i_blk].infoN[0].dimen.i > 0 &&
                 Soln_Block_List.Block[i_blk].infoN[0].dimen.j > 0) {
	        i_min0 = Soln_ptr[i_blk].ICl;
                i_max0 = Soln_ptr[i_blk].ICl+(Soln_ptr[i_blk].ICu-Soln_ptr[i_blk].ICl+1)/2;
	        i_min1 = Soln_ptr[i_blk].ICu-(Soln_ptr[i_blk].ICu-Soln_ptr[i_blk].ICl+1)/2;
                i_max1 = Soln_ptr[i_blk].ICu;
	        i_inc = 1;
	        j_min0 = Soln_ptr[i_blk].JCu;
	        j_max0 = Soln_ptr[i_blk].JCu;
	        j_min1 = Soln_ptr[i_blk].JCu;
	        j_max1 = Soln_ptr[i_blk].JCu;
	        j_inc = 1;
             } else if (Soln_Block_List.Block[i_blk].infoN[0].dimen.j > 0) {
	        i_min0 = Soln_ptr[i_blk].ICl+(Soln_ptr[i_blk].ICu-Soln_ptr[i_blk].ICl+1)/2;
	        i_max0 = Soln_ptr[i_blk].ICl;
	        i_min1 = Soln_ptr[i_blk].ICu;
	        i_max1 = Soln_ptr[i_blk].ICu-(Soln_ptr[i_blk].ICu-Soln_ptr[i_blk].ICl+1)/2;
	        i_inc = -1;
	        j_min0 = Soln_ptr[i_blk].JCu;
	        j_max0 = Soln_ptr[i_blk].JCu;
	        j_min1 = Soln_ptr[i_blk].JCu;
	        j_max1 = Soln_ptr[i_blk].JCu;
	        j_inc = 1;
             } else if (Soln_Block_List.Block[i_blk].infoN[0].dimen.i > 0) {
	        i_min0 = Soln_ptr[i_blk].ICl;
	        i_max0 = Soln_ptr[i_blk].ICl+(Soln_ptr[i_blk].ICu-Soln_ptr[i_blk].ICl+1)/2;
	        i_min1 = Soln_ptr[i_blk].ICu-(Soln_ptr[i_blk].ICu-Soln_ptr[i_blk].ICl+1)/2;
	        i_max1 = Soln_ptr[i_blk].ICu;
	        i_inc = 1;
	        j_min0 = Soln_ptr[i_blk].JCu;
	        j_max0 = Soln_ptr[i_blk].JCu;
	        j_min1 = Soln_ptr[i_blk].JCu;
	        j_max1 = Soln_ptr[i_blk].JCu;
	        j_inc = -1;
             } else {
	        i_min0 = Soln_ptr[i_blk].ICl+(Soln_ptr[i_blk].ICu-Soln_ptr[i_blk].ICl+1)/2;
	        i_max0 = Soln_ptr[i_blk].ICl;
	        i_min1 = Soln_ptr[i_blk].ICu;
	        i_max1 = Soln_ptr[i_blk].ICu-(Soln_ptr[i_blk].ICu-Soln_ptr[i_blk].ICl+1)/2;
	        i_inc = -1;
	        j_min0 = Soln_ptr[i_blk].JCu;
	        j_max0 = Soln_ptr[i_blk].JCu;
	        j_min1 = Soln_ptr[i_blk].JCu;
	        j_max1 = Soln_ptr[i_blk].JCu;
	        j_inc = -1;
             } /* endif */
	       //
	       // First neighbour (solution).
	       //
	       i = Soln_ptr[i_blk].LoadSendBuffer_BC_C2F(Soln_Block_List.message_reschange_northface_sendbuf[i_blk][0],
							 l0,buffer_size_neighbour,
							 i_min0,i_max0,i_inc,
							 j_min0,j_max0,j_inc);
	       if (i != 0) return(2200);
	       //
	       // Second neighbour (solution).
	       //
	       i = Soln_ptr[i_blk].LoadSendBuffer_BC_C2F(Soln_Block_List.message_reschange_northface_sendbuf[i_blk][1],
							 l1,buffer_size_neighbour,
							 i_min1,i_max1,i_inc,
							 j_min1,j_max1,j_inc);
	       if (i != 0) return(2201);
       } /* endif */

       //////////////////////////////////////////////////////////////////
       // Load send buffers of solution blocks with refined South face //
       // neighbour at higher mesh resolution.                         //
       //////////////////////////////////////////////////////////////////
       if (Soln_Block_List.Block[i_blk].used && 
           (Soln_Block_List.Block[i_blk].nS == 2) &&
           (Soln_Block_List.Block[i_blk].info.level < 
            Soln_Block_List.Block[i_blk].infoS[0].level)) {
	      buffer_size_neighbour = (abs(Soln_Block_List.Block[i_blk].infoS[0].dimen.i)+
                                      Soln_Block_List.Block[i_blk].infoS[0].dimen.ghost)*
                                      Soln_ptr[i_blk].NumVar();
          l0 = -1;
          l1 = -1;
          // Load ghost cell solution variable information as required.
             if (Soln_Block_List.Block[i_blk].infoS[0].dimen.i > 0 &&
                 Soln_Block_List.Block[i_blk].infoS[0].dimen.j > 0) {
	        i_min0 = Soln_ptr[i_blk].ICl;
                i_max0 = Soln_ptr[i_blk].ICl+(Soln_ptr[i_blk].ICu-Soln_ptr[i_blk].ICl+1)/2;
	        i_min1 = Soln_ptr[i_blk].ICu-(Soln_ptr[i_blk].ICu-Soln_ptr[i_blk].ICl+1)/2;
                i_max1 = Soln_ptr[i_blk].ICu;
	        i_inc = 1;
	        j_min0 = Soln_ptr[i_blk].JCl;
	        j_max0 = Soln_ptr[i_blk].JCl;
	        j_min1 = Soln_ptr[i_blk].JCl;
	        j_max1 = Soln_ptr[i_blk].JCl;
	        j_inc = 1;
             } else if (Soln_Block_List.Block[i_blk].infoS[0].dimen.j > 0) {
	        i_min0 = Soln_ptr[i_blk].ICl+(Soln_ptr[i_blk].ICu-Soln_ptr[i_blk].ICl+1)/2;
	        i_max0 = Soln_ptr[i_blk].ICl;
	        i_min1 = Soln_ptr[i_blk].ICu;
	        i_max1 = Soln_ptr[i_blk].ICu-(Soln_ptr[i_blk].ICu-Soln_ptr[i_blk].ICl+1)/2;
	        i_inc = -1;
	        j_min0 = Soln_ptr[i_blk].JCl;
	        j_max0 = Soln_ptr[i_blk].JCl;
	        j_min1 = Soln_ptr[i_blk].JCl;
	        j_max1 = Soln_ptr[i_blk].JCl;
	        j_inc = 1;
             } else if (Soln_Block_List.Block[i_blk].infoS[0].dimen.i > 0) {
	        i_min0 = Soln_ptr[i_blk].ICl;
	        i_max0 = Soln_ptr[i_blk].ICl+(Soln_ptr[i_blk].ICu-Soln_ptr[i_blk].ICl+1)/2;
	        i_min1 = Soln_ptr[i_blk].ICu-(Soln_ptr[i_blk].ICu-Soln_ptr[i_blk].ICl+1)/2;
	        i_max1 = Soln_ptr[i_blk].ICu;
	        i_inc = 1;
	        j_min0 = Soln_ptr[i_blk].JCl;
	        j_max0 = Soln_ptr[i_blk].JCl;
	        j_min1 = Soln_ptr[i_blk].JCl;
	        j_max1 = Soln_ptr[i_blk].JCl;
	        j_inc = -1;
             } else {
	        i_min0 = Soln_ptr[i_blk].ICl+(Soln_ptr[i_blk].ICu-Soln_ptr[i_blk].ICl+1)/2;
	        i_max0 = Soln_ptr[i_blk].ICl;
	        i_min1 = Soln_ptr[i_blk].ICu;
	        i_max1 = Soln_ptr[i_blk].ICu-(Soln_ptr[i_blk].ICu-Soln_ptr[i_blk].ICl+1)/2;
	        i_inc = -1;
	        j_min0 = Soln_ptr[i_blk].JCl;
	        j_max0 = Soln_ptr[i_blk].JCl;
	        j_min1 = Soln_ptr[i_blk].JCl;
	        j_max1 = Soln_ptr[i_blk].JCl;
	        j_inc = -1;
             } /* endif */
	       //
	       // First neighbour (solution).
	       //
	       i = Soln_ptr[i_blk].LoadSendBuffer_BC_C2F(Soln_Block_List.message_reschange_southface_sendbuf[i_blk][0],
							 l0,buffer_size_neighbour,
							 i_min0,i_max0,i_inc,
							 j_min0,j_max0,j_inc);
	       if (i != 0) return(2100);
	       //
	       // Second neighbour (solution).
	       //
	       i = Soln_ptr[i_blk].LoadSendBuffer_BC_C2F(Soln_Block_List.message_reschange_southface_sendbuf[i_blk][1],
							 l1,buffer_size_neighbour,
							 i_min1,i_max1,i_inc,
							 j_min1,j_max1,j_inc);
	       if (i != 0) return(2101);
       } /* endif */

       //////////////////////////////////////////////////////////////////
       // Load send buffers of solution blocks with refined East face //
       // neighbour at higher mesh resolution.                         //
       //////////////////////////////////////////////////////////////////
       if (Soln_Block_List.Block[i_blk].used && 
           (Soln_Block_List.Block[i_blk].nE == 2) &&
           (Soln_Block_List.Block[i_blk].info.level < 
            Soln_Block_List.Block[i_blk].infoE[0].level)) {
	      buffer_size_neighbour = (abs(Soln_Block_List.Block[i_blk].infoE[0].dimen.j)+
                                      Soln_Block_List.Block[i_blk].infoE[0].dimen.ghost)*
                                      Soln_ptr[i_blk].NumVar();
          l0 = -1;
          l1 = -1;
          // Load ghost cell solution variable information as required.
             if (Soln_Block_List.Block[i_blk].infoE[0].dimen.i > 0 &&
                 Soln_Block_List.Block[i_blk].infoE[0].dimen.j > 0) {
	        i_min0 = Soln_ptr[i_blk].ICu;
	        i_max0 = Soln_ptr[i_blk].ICu;
	        i_min1 = Soln_ptr[i_blk].ICu;
	        i_max1 = Soln_ptr[i_blk].ICu;
	        i_inc = 1;
	        j_min0 = Soln_ptr[i_blk].JCl;
                j_max0 = Soln_ptr[i_blk].JCl+(Soln_ptr[i_blk].JCu-Soln_ptr[i_blk].JCl+1)/2;
	        j_min1 = Soln_ptr[i_blk].JCu-(Soln_ptr[i_blk].JCu-Soln_ptr[i_blk].JCl+1)/2;
                j_max1 = Soln_ptr[i_blk].JCu;
	        j_inc = 1;
             } else if (Soln_Block_List.Block[i_blk].infoE[0].dimen.j > 0) {
	        i_min0 = Soln_ptr[i_blk].ICu;
	        i_max0 = Soln_ptr[i_blk].ICu;
	        i_min1 = Soln_ptr[i_blk].ICu;
	        i_max1 = Soln_ptr[i_blk].ICu;
	        i_inc = -1;
	        j_min0 = Soln_ptr[i_blk].JCl;
                j_max0 = Soln_ptr[i_blk].JCl+(Soln_ptr[i_blk].JCu-Soln_ptr[i_blk].JCl+1)/2;
	        j_min1 = Soln_ptr[i_blk].JCu-(Soln_ptr[i_blk].JCu-Soln_ptr[i_blk].JCl+1)/2;
                j_max1 = Soln_ptr[i_blk].JCu;
	        j_inc = 1;
             } else if (Soln_Block_List.Block[i_blk].infoE[0].dimen.i > 0) {
	        i_min0 = Soln_ptr[i_blk].ICu;
	        i_max0 = Soln_ptr[i_blk].ICu;
	        i_min1 = Soln_ptr[i_blk].ICu;
	        i_max1 = Soln_ptr[i_blk].ICu;
	        i_inc = 1;
	        j_min0 = Soln_ptr[i_blk].JCl+(Soln_ptr[i_blk].JCu-Soln_ptr[i_blk].JCl+1)/2;
                j_max0 = Soln_ptr[i_blk].JCl;
	        j_min1 = Soln_ptr[i_blk].JCu;
                j_max1 = Soln_ptr[i_blk].JCu-(Soln_ptr[i_blk].JCu-Soln_ptr[i_blk].JCl+1)/2;
	        j_inc = -1;
             } else {
	        i_min0 = Soln_ptr[i_blk].ICu;
	        i_max0 = Soln_ptr[i_blk].ICu;
	        i_min1 = Soln_ptr[i_blk].ICu;
	        i_max1 = Soln_ptr[i_blk].ICu;
	        i_inc = -1;
	        j_min0 = Soln_ptr[i_blk].JCl+(Soln_ptr[i_blk].JCu-Soln_ptr[i_blk].JCl+1)/2;
                j_max0 = Soln_ptr[i_blk].JCl;
	        j_min1 = Soln_ptr[i_blk].JCu;
                j_max1 = Soln_ptr[i_blk].JCu-(Soln_ptr[i_blk].JCu-Soln_ptr[i_blk].JCl+1)/2;
	        j_inc = -1;
             } /* endif */
	       //
	       // First neighbour (solution).
	       //
	       i = Soln_ptr[i_blk].LoadSendBuffer_BC_C2F(Soln_Block_List.message_reschange_eastface_sendbuf[i_blk][0],
							 l0,buffer_size_neighbour,
							 i_min0,i_max0,i_inc,
							 j_min0,j_max0,j_inc);
	       if (i != 0) return(1200);
	       //
	       // Second neighbour (solution).
	       //
	       i = Soln_ptr[i_blk].LoadSendBuffer_BC_C2F(Soln_Block_List.message_reschange_eastface_sendbuf[i_blk][1],
							 l1,buffer_size_neighbour,
							 i_min1,i_max1,i_inc,
							 j_min1,j_max1,j_inc);
	       if (i != 0) return(1201);
       } /* endif */

       //////////////////////////////////////////////////////////////////
       // Load send buffers of solution blocks with refined West face //
       // neighbour at higher mesh resolution.                         //
       //////////////////////////////////////////////////////////////////
       if (Soln_Block_List.Block[i_blk].used && 
           (Soln_Block_List.Block[i_blk].nW == 2) &&
           (Soln_Block_List.Block[i_blk].info.level < 
            Soln_Block_List.Block[i_blk].infoW[0].level)) {
	      buffer_size_neighbour = (abs(Soln_Block_List.Block[i_blk].infoW[0].dimen.j)+
				       Soln_Block_List.Block[i_blk].infoW[0].dimen.ghost)*
                                      Soln_ptr[i_blk].NumVar();
          l0 = -1;
          l1 = -1;
          // Load ghost cell solution variable information as required.
             if (Soln_Block_List.Block[i_blk].infoW[0].dimen.i > 0 &&
                 Soln_Block_List.Block[i_blk].infoW[0].dimen.j > 0) {
	        i_min0 = Soln_ptr[i_blk].ICl;
	        i_max0 = Soln_ptr[i_blk].ICl;
	        i_min1 = Soln_ptr[i_blk].ICl;
	        i_max1 = Soln_ptr[i_blk].ICl;
	        i_inc = 1;
	        j_min0 = Soln_ptr[i_blk].JCl;
                j_max0 = Soln_ptr[i_blk].JCl+(Soln_ptr[i_blk].JCu-Soln_ptr[i_blk].JCl+1)/2;
	        j_min1 = Soln_ptr[i_blk].JCu-(Soln_ptr[i_blk].JCu-Soln_ptr[i_blk].JCl+1)/2;
                j_max1 = Soln_ptr[i_blk].JCu;
	        j_inc = 1;
             } else if (Soln_Block_List.Block[i_blk].infoW[0].dimen.j > 0) {
	        i_min0 = Soln_ptr[i_blk].ICl;
	        i_max0 = Soln_ptr[i_blk].ICl;
	        i_min1 = Soln_ptr[i_blk].ICl;
	        i_max1 = Soln_ptr[i_blk].ICl;
	        i_inc = -1;
	        j_min0 = Soln_ptr[i_blk].JCl;
                j_max0 = Soln_ptr[i_blk].JCl+(Soln_ptr[i_blk].JCu-Soln_ptr[i_blk].JCl+1)/2;
	        j_min1 = Soln_ptr[i_blk].JCu-(Soln_ptr[i_blk].JCu-Soln_ptr[i_blk].JCl+1)/2;
                j_max1 = Soln_ptr[i_blk].JCu;
	        j_inc = 1;
             } else if (Soln_Block_List.Block[i_blk].infoW[0].dimen.i > 0) {
	        i_min0 = Soln_ptr[i_blk].ICl;
	        i_max0 = Soln_ptr[i_blk].ICl;
	        i_min1 = Soln_ptr[i_blk].ICl;
	        i_max1 = Soln_ptr[i_blk].ICl;
	        i_inc = 1;
	        j_min0 = Soln_ptr[i_blk].JCl+(Soln_ptr[i_blk].JCu-Soln_ptr[i_blk].JCl+1)/2;
                j_max0 = Soln_ptr[i_blk].JCl;
	        j_min1 = Soln_ptr[i_blk].JCu;
                j_max1 = Soln_ptr[i_blk].JCu-(Soln_ptr[i_blk].JCu-Soln_ptr[i_blk].JCl+1)/2;
	        j_inc = -1;
             } else {
	        i_min0 = Soln_ptr[i_blk].ICl;
	        i_max0 = Soln_ptr[i_blk].ICl;
	        i_min1 = Soln_ptr[i_blk].ICl;
	        i_max1 = Soln_ptr[i_blk].ICl;
	        i_inc = -1;
	        j_min0 = Soln_ptr[i_blk].JCl+(Soln_ptr[i_blk].JCu-Soln_ptr[i_blk].JCl+1)/2;
                j_max0 = Soln_ptr[i_blk].JCl;
	        j_min1 = Soln_ptr[i_blk].JCu;
                j_max1 = Soln_ptr[i_blk].JCu-(Soln_ptr[i_blk].JCu-Soln_ptr[i_blk].JCl+1)/2;
	        j_inc = -1;
             } /* endif */
	       //
	       // First neighbour (solution).
	       //
	       i = Soln_ptr[i_blk].LoadSendBuffer_BC_C2F(Soln_Block_List.message_reschange_westface_sendbuf[i_blk][0],
							 l0,buffer_size_neighbour,
							 i_min0,i_max0,i_inc,
							 j_min0,j_max0,j_inc);
	       if (i != 0) return(1100);
	       //
	       // Second neighbour (solution).
	       //
	       i = Soln_ptr[i_blk].LoadSendBuffer_BC_C2F(Soln_Block_List.message_reschange_westface_sendbuf[i_blk][1],
							 l1,buffer_size_neighbour,
							 i_min1,i_max1,i_inc,
							 j_min1,j_max1,j_inc);
	       if (i != 0) return(1101);
       } /* endif */


    }  /* endfor */

    /* Loading of send buffers complete.  Return zero value. */

    return(0);

}

/********************************************************
 * Routine: Unload_Receive_Message_Buffers_NoResChange  *
 *                                                      *
 * Unloads the receive message passing buffers for      *
 * receiving solution information from neighbouring 2D  *
 * quadrilateral multi-block solution blocks with no    *
 * mesh resolution changes (works on entire 1D array of *
 * 2D solution blocks).                                 *
 *                                                      *
 ********************************************************/
int Unload_Receive_Message_Buffers_NoResChange(Rte2D_Quad_Block *Soln_ptr,
                                               AdaptiveBlock2D_List &Soln_Block_List,
                                               const int Number_of_Solution_Variables) {

    int i_blk, buffer_size, 
        i_min, i_max, i_inc, i, 
        j_min, j_max, j_inc, j, 
        l;
    Vector2D x_ref;

    /* Check to see if the number of solution variables specified
       in the calling argument is correct. */

    if (Number_of_Solution_Variables != Soln_ptr[0].NumVar()) {
       return(1);
    } /* endif */

    /* Unload the receive buffers for each solution block. */

    for ( i_blk = 0 ; i_blk <= Soln_Block_List.Nblk-1 ; ++i_blk ) {

       /////////////////////////////////////////////////////////////////////////
       // Unload receive buffer for North face neighbours of solution blocks. //
       /////////////////////////////////////////////////////////////////////////
       if (Soln_Block_List.Block[i_blk].used && 
           (Soln_Block_List.Block[i_blk].nN == 1) &&
           (Soln_Block_List.Block[i_blk].info.level == 
            Soln_Block_List.Block[i_blk].infoN[0].level)) {
	      buffer_size = abs(Soln_Block_List.Block[i_blk].info.dimen.i)*
                            Soln_ptr[i_blk].NumVar();	      
          l = -1;
          // Unload ghost cell solution information as required.
	      i_min = Soln_ptr[i_blk].ICl;
	      i_max = Soln_ptr[i_blk].ICu;
	      i_inc = 1;
	      j_min = Soln_ptr[i_blk].JCu;
	      j_max = Soln_ptr[i_blk].JCu;
	      j_inc = 1;
	      i = Soln_ptr[i_blk].UnloadReceiveBuffer_BC(Soln_Block_List.message_noreschange_northface_recbuf[i_blk],
							 l,buffer_size,
							 i_min,i_max,i_inc,
							 j_min,j_max,j_inc);
	      if (i != 0) return(2200);
       } /* endif */

       /////////////////////////////////////////////////////////////////////////
       // Unload receive buffer for South face neighbours of solution blocks. //
       /////////////////////////////////////////////////////////////////////////
       if (Soln_Block_List.Block[i_blk].used && 
           (Soln_Block_List.Block[i_blk].nS == 1) &&
           (Soln_Block_List.Block[i_blk].info.level == 
            Soln_Block_List.Block[i_blk].infoS[0].level)) {
	      buffer_size = abs(Soln_Block_List.Block[i_blk].info.dimen.i)*
                            Soln_ptr[i_blk].NumVar();
          l = -1;
          // Unload ghost cell solution information as required.
	      i_min = Soln_ptr[i_blk].ICl;
	      i_max = Soln_ptr[i_blk].ICu;
	      i_inc = 1;
	      j_min = Soln_ptr[i_blk].JCl;
	      j_max = Soln_ptr[i_blk].JCl;
	      j_inc = 1;
	      i = Soln_ptr[i_blk].UnloadReceiveBuffer_BC(Soln_Block_List.message_noreschange_southface_recbuf[i_blk],
							 l,buffer_size,
							 i_min,i_max,i_inc,
							 j_min,j_max,j_inc);
	      if (i != 0) return(2100);
       } /* endif */

       ////////////////////////////////////////////////////////////////////////
       // Unload receive buffer for East face neighbours of solution blocks. //
       ////////////////////////////////////////////////////////////////////////
       if (Soln_Block_List.Block[i_blk].used && 
           (Soln_Block_List.Block[i_blk].nE == 1) &&
           (Soln_Block_List.Block[i_blk].info.level == 
            Soln_Block_List.Block[i_blk].infoE[0].level)) {
	      buffer_size = abs(Soln_Block_List.Block[i_blk].info.dimen.j)*
                            Soln_ptr[i_blk].NumVar();
          l = -1;
          // Unload ghost cell solution information as required.
	      i_min = Soln_ptr[i_blk].ICu;
	      i_max = Soln_ptr[i_blk].ICu;
	      i_inc = 1;
	      j_min = Soln_ptr[i_blk].JCl;
	      j_max = Soln_ptr[i_blk].JCu;
	      j_inc = 1;
	      i = Soln_ptr[i_blk].UnloadReceiveBuffer_BC(Soln_Block_List.message_noreschange_eastface_recbuf[i_blk],
							 l,buffer_size,
							 i_min,i_max,i_inc,
							 j_min,j_max,j_inc);
	      if (i != 0) return(1200);	      
       } /* endif */

       ////////////////////////////////////////////////////////////////////////
       // Unload receive buffer for West face neighbours of solution blocks. //
       ////////////////////////////////////////////////////////////////////////
       if (Soln_Block_List.Block[i_blk].used && 
           (Soln_Block_List.Block[i_blk].nW == 1) &&
           (Soln_Block_List.Block[i_blk].info.level == 
            Soln_Block_List.Block[i_blk].infoW[0].level)) {
	      buffer_size = abs(Soln_Block_List.Block[i_blk].info.dimen.j)*
                            Soln_ptr[i_blk].NumVar();
          l = -1;
          // Unload ghost cell solution information as required.
	      i_min = Soln_ptr[i_blk].ICl;
	      i_max = Soln_ptr[i_blk].ICl;
	      i_inc = 1;
	      j_min = Soln_ptr[i_blk].JCl;
	      j_max = Soln_ptr[i_blk].JCu;
	      j_inc = 1;
	      i = Soln_ptr[i_blk].UnloadReceiveBuffer_BC(Soln_Block_List.message_noreschange_westface_recbuf[i_blk],
							 l,buffer_size,
							 i_min,i_max,i_inc,
							 j_min,j_max,j_inc);
	      if (i != 0) return(1100);
       } /* endif */


    }  /* endfor */

    /* Unloading of receive buffers complete.  Return zero value. */

    return(0);

}

/******************************************************************
 * Routine: Unload_Receive_Message_Buffers_ResChange_FineToCoarse *
 *                                                                *
 * Unloads the receive message passing buffers as sent from more  *
 * refined (higher mesh resolution) solution blocks to more       *
 * coarse neighbouring 2D quadrilateral multi-block solution      *
 * blocks with lower mesh resolution (works on entire 1D array of *
 * 2D solution blocks).                                           *
 *                                                                *
 ******************************************************************/
int Unload_Receive_Message_Buffers_ResChange_FineToCoarse(Rte2D_Quad_Block *Soln_ptr,
                                                          AdaptiveBlock2D_List &Soln_Block_List,
                                                          const int Number_of_Solution_Variables) {

    int i_blk, buffer_size, 
        i_min, i_max, i_inc, i, 
        j_min, j_max, j_inc, j, 
        l, j_neigh;
    Vector2D x_ref;

    /* Check to see if the number of solution variables specified
       in the calling argument is correct. */

    if (Number_of_Solution_Variables != Soln_ptr[0].NumVar()) {
       return(1);
    } /* endif */

    /* Unload the receive buffers of each solution block. */

    for ( i_blk = 0 ; i_blk <= Soln_Block_List.Nblk-1 ; ++i_blk ) {

       ///////////////////////////////////////////////////////////////////////
       // Unload receive buffer for solution blocks with refined North face //
       // neighbours at higher mesh resolution.                             //
       ///////////////////////////////////////////////////////////////////////
       if (Soln_Block_List.Block[i_blk].used && 
           (Soln_Block_List.Block[i_blk].nN == 2) &&
           (Soln_Block_List.Block[i_blk].info.level < 
            Soln_Block_List.Block[i_blk].infoN[0].level)) {
                buffer_size = (abs(Soln_Block_List.Block[i_blk].info.dimen.i)/2)*
                              Soln_ptr[i_blk].NumVar();
          for (j_neigh = 0; j_neigh <= Soln_Block_List.Block[i_blk].nN-1 ; ++j_neigh) {
             l = -1;
             // Unload ghost cell solution information as required.
                   if (j_neigh == 0) {
                      i_min = Soln_ptr[i_blk].ICl;
	              i_max = Soln_ptr[i_blk].ICl+(Soln_ptr[i_blk].ICu-Soln_ptr[i_blk].ICl-1)/2;
	              i_inc = 1;
	              j_min = Soln_ptr[i_blk].JCu;
	              j_max = Soln_ptr[i_blk].JCu;
	              j_inc = 1;
                   } else {
                      i_min = Soln_ptr[i_blk].ICl+(Soln_ptr[i_blk].ICu-Soln_ptr[i_blk].ICl+1)/2;
	              i_max = Soln_ptr[i_blk].ICu;
	              i_inc = 1;
	              j_min = Soln_ptr[i_blk].JCu;
	              j_max = Soln_ptr[i_blk].JCu;
	              j_inc = 1;
                   } /* endif */
		     i = Soln_ptr[i_blk].UnloadReceiveBuffer_BC_F2C(Soln_Block_List.message_reschange_northface_recbuf[i_blk][j_neigh],
								    l,buffer_size,
								    i_min,i_max,i_inc,
								    j_min,j_max,j_inc);
	           if (i != 0) return(2201);
          } /* endfor */
       } /* endif */

       ///////////////////////////////////////////////////////////////////////
       // Unload receive buffer for solution blocks with refined South face //
       // neighbours at higher mesh resolution.                             //
       ///////////////////////////////////////////////////////////////////////
       if (Soln_Block_List.Block[i_blk].used && 
           (Soln_Block_List.Block[i_blk].nS == 2) &&
           (Soln_Block_List.Block[i_blk].info.level < 
            Soln_Block_List.Block[i_blk].infoS[0].level)) {
                buffer_size = (abs(Soln_Block_List.Block[i_blk].info.dimen.i)/2)*
                              Soln_ptr[i_blk].NumVar();
          for (j_neigh = 0; j_neigh <= Soln_Block_List.Block[i_blk].nS-1 ; ++j_neigh) {
             l = -1;
             // Unload ghost cell solution information as required.
                   if (j_neigh == 0) {
                      i_min = Soln_ptr[i_blk].ICl;
	              i_max = Soln_ptr[i_blk].ICl+(Soln_ptr[i_blk].ICu-Soln_ptr[i_blk].ICl-1)/2;
	              i_inc = 1;
	              j_min = Soln_ptr[i_blk].JCl;
	              j_max = Soln_ptr[i_blk].JCl;
	              j_inc = 1;
                   } else {
                      i_min = Soln_ptr[i_blk].ICl+(Soln_ptr[i_blk].ICu-Soln_ptr[i_blk].ICl+1)/2;
	              i_max = Soln_ptr[i_blk].ICu;
	              i_inc = 1;
	              j_min = Soln_ptr[i_blk].JCl;
	              j_max = Soln_ptr[i_blk].JCl;
	              j_inc = 1;
                   } /* endif */
		     i = Soln_ptr[i_blk].UnloadReceiveBuffer_BC_F2C(Soln_Block_List.message_reschange_southface_recbuf[i_blk][j_neigh],
								    l,buffer_size,
								    i_min,i_max,i_inc,
								    j_min,j_max,j_inc);
		   if (i != 0) return(2101);
          } /* endfor */
       } /* endif */

       //////////////////////////////////////////////////////////////////////
       // Unload receive buffer for solution blocks with refined East face //
       // neighbours at higher mesh resolution.                            //
       //////////////////////////////////////////////////////////////////////
       if (Soln_Block_List.Block[i_blk].used && 
           (Soln_Block_List.Block[i_blk].nE == 2) &&
           (Soln_Block_List.Block[i_blk].info.level < 
            Soln_Block_List.Block[i_blk].infoE[0].level)) {
                buffer_size = (abs(Soln_Block_List.Block[i_blk].info.dimen.j)/2)*
                              Soln_ptr[i_blk].NumVar();
          for (j_neigh = 0; j_neigh <= Soln_Block_List.Block[i_blk].nE-1 ; ++j_neigh) {
             l = -1;
             // Unload ghost cell solution information as required.
                   if (j_neigh == 0) {
                      i_min = Soln_ptr[i_blk].ICu;
	              i_max = Soln_ptr[i_blk].ICu;
	              i_inc = 1;
	              j_min = Soln_ptr[i_blk].JCl;
	              j_max = Soln_ptr[i_blk].JCl+(Soln_ptr[i_blk].JCu-Soln_ptr[i_blk].JCl-1)/2;
	              j_inc = 1;
                   } else {
                      i_min = Soln_ptr[i_blk].ICu;
	              i_max = Soln_ptr[i_blk].ICu;
	              i_inc = 1;
	              j_min = Soln_ptr[i_blk].JCl+(Soln_ptr[i_blk].JCu-Soln_ptr[i_blk].JCl+1)/2;
	              j_max = Soln_ptr[i_blk].JCu;
	              j_inc = 1;
                   } /* endif */
		     i = Soln_ptr[i_blk].UnloadReceiveBuffer_BC_F2C(Soln_Block_List.message_reschange_eastface_recbuf[i_blk][j_neigh],
								    l,buffer_size,
								    i_min,i_max,i_inc,
								    j_min,j_max,j_inc);		     
	           if (i != 0) return(1201);
          } /* endfor */
       } /* endif */

       //////////////////////////////////////////////////////////////////////
       // Unload receive buffer for solution blocks with refined West face //
       // neighbours at higher mesh resolution.                            //
       //////////////////////////////////////////////////////////////////////
       if (Soln_Block_List.Block[i_blk].used && 
           (Soln_Block_List.Block[i_blk].nW == 2) &&
           (Soln_Block_List.Block[i_blk].info.level < 
            Soln_Block_List.Block[i_blk].infoW[0].level)) {
                buffer_size = (abs(Soln_Block_List.Block[i_blk].info.dimen.j)/2)*
                              Soln_ptr[i_blk].NumVar();
          for (j_neigh = 0; j_neigh <= Soln_Block_List.Block[i_blk].nW-1 ; ++j_neigh) {
             l = -1;
             // Unload ghost cell solution information as required.
                   if (j_neigh == 0) {
                      i_min = Soln_ptr[i_blk].ICl;
	              i_max = Soln_ptr[i_blk].ICl;
	              i_inc = 1;
	              j_min = Soln_ptr[i_blk].JCl;
	              j_max = Soln_ptr[i_blk].JCl+(Soln_ptr[i_blk].JCu-Soln_ptr[i_blk].JCl-1)/2;
	              j_inc = 1;
                   } else {
                      i_min = Soln_ptr[i_blk].ICl;
	              i_max = Soln_ptr[i_blk].ICl;
	              i_inc = 1;
	              j_min = Soln_ptr[i_blk].JCl+(Soln_ptr[i_blk].JCu-Soln_ptr[i_blk].JCl+1)/2;
	              j_max = Soln_ptr[i_blk].JCu;
	              j_inc = 1;
                   } /* endif */
		     i = Soln_ptr[i_blk].UnloadReceiveBuffer_BC_F2C(Soln_Block_List.message_reschange_westface_recbuf[i_blk][j_neigh],
								    l,buffer_size,
								    i_min,i_max,i_inc,
								    j_min,j_max,j_inc);
	           if (i != 0) return(1101);
          } /* endfor */
       } /* endif */


    }  /* endfor */

    /* Unloading of receive buffers complete.  Return zero value. */

    return(0);

}

/******************************************************************
 * Routine: Unload_Receive_Message_Buffers_ResChange_CoarseToFine *
 *                                                                *
 * Unloads the receive message passing buffers as sent from more  *
 * coarse (lower mesh resolution) solution blocks to more         *
 * refined neighbouring 2D quadrilateral multi-block solution     *
 * blocks with higher mesh resolution (works on entire 1D array   *
 * of 2D solution blocks).                                        *
 *                                                                *
 ******************************************************************/
int Unload_Receive_Message_Buffers_ResChange_CoarseToFine(Rte2D_Quad_Block *Soln_ptr,
                                                          AdaptiveBlock2D_List &Soln_Block_List,
                                                          const int Number_of_Solution_Variables) {

    int i_blk, buffer_size, 
        i_min, i_max, i_inc, i, 
        j_min, j_max, j_inc, j, 
        l, j_neigh;
    Vector2D x_ref;

    /* Check to see if the number of solution variables specified
       in the calling argument is correct. */

    if (Number_of_Solution_Variables != Soln_ptr[0].NumVar()) {
       return(1);
    } /* endif */

    /* Unload the receive buffers of each solution block. */

    for ( i_blk = 0 ; i_blk <= Soln_Block_List.Nblk-1 ; ++i_blk ) {

       ///////////////////////////////////////////////////////////////////////
       // Unload receive buffer for solution blocks with coarse North face  //
       // neighbours at lower mesh resolution.                              //
       ///////////////////////////////////////////////////////////////////////
       if (Soln_Block_List.Block[i_blk].used && 
           (Soln_Block_List.Block[i_blk].nN == 1) &&
           (Soln_Block_List.Block[i_blk].info.level > 
            Soln_Block_List.Block[i_blk].infoN[0].level)) {
	      buffer_size = (abs(Soln_Block_List.Block[i_blk].info.dimen.i)+
                            Soln_Block_List.Block[i_blk].info.dimen.ghost)*
                            Soln_ptr[i_blk].NumVar();
          l = -1;
          // Unload ghost cell solution information as required.
	      if (Soln_Block_List.Block[i_blk].info.sector == ADAPTIVEBLOCK2D_SECTOR_NW) {
                i_min = Soln_ptr[i_blk].ICl;
	        i_max = Soln_ptr[i_blk].ICu+Soln_Block_List.Block[i_blk].info.dimen.ghost;
                i_inc = 1;
                j_min = Soln_ptr[i_blk].JCu;
                j_max = Soln_ptr[i_blk].JCu;
                j_inc = 1;
	      } else if (Soln_Block_List.Block[i_blk].info.sector == ADAPTIVEBLOCK2D_SECTOR_NE) {
                i_min = Soln_ptr[i_blk].ICl-Soln_Block_List.Block[i_blk].info.dimen.ghost;
	        i_max = Soln_ptr[i_blk].ICu;
                i_inc = 1;
                j_min = Soln_ptr[i_blk].JCu;
                j_max = Soln_ptr[i_blk].JCu;
                j_inc = 1;
	      } /* endif */
	      i = Soln_ptr[i_blk].UnloadReceiveBuffer_BC_C2F(Soln_Block_List.message_reschange_northface_recbuf[i_blk][0],
							     l,buffer_size,
							     i_min,i_max,i_inc,
							     j_min,j_max,j_inc);
	      if (i != 0) return(2200);
       } /* endif */

       ///////////////////////////////////////////////////////////////////////
       // Unload receive buffer for solution blocks with coarse South face  //
       // neighbours at lower mesh resolution.                              //
       ///////////////////////////////////////////////////////////////////////
       if (Soln_Block_List.Block[i_blk].used && 
           (Soln_Block_List.Block[i_blk].nS == 1) &&
           (Soln_Block_List.Block[i_blk].info.level > 
            Soln_Block_List.Block[i_blk].infoS[0].level)) {
	      buffer_size = (abs(Soln_Block_List.Block[i_blk].info.dimen.i)+
                            Soln_Block_List.Block[i_blk].info.dimen.ghost)*
                            Soln_ptr[i_blk].NumVar();
          l = -1;
          // Unload ghost cell solution information as required.
	      if (Soln_Block_List.Block[i_blk].info.sector == ADAPTIVEBLOCK2D_SECTOR_SW) {
                i_min = Soln_ptr[i_blk].ICl;
	        i_max = Soln_ptr[i_blk].ICu+Soln_Block_List.Block[i_blk].info.dimen.ghost;
                i_inc = 1;
                j_min = Soln_ptr[i_blk].JCl;
                j_max = Soln_ptr[i_blk].JCl;
                j_inc = 1;
	      } else if (Soln_Block_List.Block[i_blk].info.sector == ADAPTIVEBLOCK2D_SECTOR_SE) {
                i_min = Soln_ptr[i_blk].ICl-Soln_Block_List.Block[i_blk].info.dimen.ghost;
	        i_max = Soln_ptr[i_blk].ICu;
                i_inc = 1;
                j_min = Soln_ptr[i_blk].JCl;
                j_max = Soln_ptr[i_blk].JCl;
                j_inc = 1;
	      } /* endif */
	      i = Soln_ptr[i_blk].UnloadReceiveBuffer_BC_C2F(Soln_Block_List.message_reschange_southface_recbuf[i_blk][0],
							     l,buffer_size,
							     i_min,i_max,i_inc,
							     j_min,j_max,j_inc);
	      if (i != 0) return(2100);
       } /* endif */

       ///////////////////////////////////////////////////////////////////////
       // Unload receive buffer for solution blocks with coarse East face   //
       // neighbours at lower mesh resolution.                              //
       ///////////////////////////////////////////////////////////////////////
       if (Soln_Block_List.Block[i_blk].used && 
           (Soln_Block_List.Block[i_blk].nE == 1) &&
           (Soln_Block_List.Block[i_blk].info.level > 
            Soln_Block_List.Block[i_blk].infoE[0].level)) {
	      buffer_size = (abs(Soln_Block_List.Block[i_blk].info.dimen.j)+
                            Soln_Block_List.Block[i_blk].info.dimen.ghost)*
                            Soln_ptr[i_blk].NumVar();
          l = -1;
          // Unload ghost cell solution information as required.
	      if (Soln_Block_List.Block[i_blk].info.sector == ADAPTIVEBLOCK2D_SECTOR_SE) {
                i_min = Soln_ptr[i_blk].ICu;
                i_max = Soln_ptr[i_blk].ICu;
                i_inc = 1;
                j_min = Soln_ptr[i_blk].JCl;
	        j_max = Soln_ptr[i_blk].JCu+Soln_Block_List.Block[i_blk].info.dimen.ghost;
                j_inc = 1;
	      } else if (Soln_Block_List.Block[i_blk].info.sector == ADAPTIVEBLOCK2D_SECTOR_NE) {
                i_min = Soln_ptr[i_blk].ICu;
                i_max = Soln_ptr[i_blk].ICu;
                i_inc = 1;
                j_min = Soln_ptr[i_blk].JCl-Soln_Block_List.Block[i_blk].info.dimen.ghost;
	        j_max = Soln_ptr[i_blk].JCu;
                j_inc = 1;
	      } /* endif */
	      i = Soln_ptr[i_blk].UnloadReceiveBuffer_BC_C2F(Soln_Block_List.message_reschange_eastface_recbuf[i_blk][0],
							     l,buffer_size,
							     i_min,i_max,i_inc,
							     j_min,j_max,j_inc);
	      if (i != 0) return(1200);
       } /* endif */

       ///////////////////////////////////////////////////////////////////////
       // Unload receive buffer for solution blocks with coarse West face   //
       // neighbours at lower mesh resolution.                              //
       ///////////////////////////////////////////////////////////////////////
       if (Soln_Block_List.Block[i_blk].used && 
           (Soln_Block_List.Block[i_blk].nW == 1) &&
           (Soln_Block_List.Block[i_blk].info.level > 
            Soln_Block_List.Block[i_blk].infoW[0].level)) {
	      buffer_size = (abs(Soln_Block_List.Block[i_blk].info.dimen.j)+
                            Soln_Block_List.Block[i_blk].info.dimen.ghost)*
                            Soln_ptr[i_blk].NumVar();	      
          l = -1;
          // Unload ghost cell solution information as required.
	      if (Soln_Block_List.Block[i_blk].info.sector == ADAPTIVEBLOCK2D_SECTOR_SW) {
                i_min = Soln_ptr[i_blk].ICl;
                i_max = Soln_ptr[i_blk].ICl;
                i_inc = 1;
                j_min = Soln_ptr[i_blk].JCl;
	        j_max = Soln_ptr[i_blk].JCu+Soln_Block_List.Block[i_blk].info.dimen.ghost;
                j_inc = 1;
	      } else if (Soln_Block_List.Block[i_blk].info.sector == ADAPTIVEBLOCK2D_SECTOR_NW) {
                i_min = Soln_ptr[i_blk].ICl;
                i_max = Soln_ptr[i_blk].ICl;
                i_inc = 1;
                j_min = Soln_ptr[i_blk].JCl-Soln_Block_List.Block[i_blk].info.dimen.ghost;
	        j_max = Soln_ptr[i_blk].JCu;
                j_inc = 1;
	      } /* endif */
	      i = Soln_ptr[i_blk].UnloadReceiveBuffer_BC_C2F(Soln_Block_List.message_reschange_westface_recbuf[i_blk][0],
							     l,buffer_size,
							     i_min,i_max,i_inc,
							     j_min,j_max,j_inc);
	      if (i != 0) return(1100);
       } /* endif */


    }  /* endfor */

    /* Unloading of receive buffers complete.  Return zero value. */

    return(0);

}



/********************************************************
 * Routine: Send_Boundary_Ref_State                     *
 *                                                      *
 * Loads, sends and exchanges, and then unloads all     *
 * messages for sharing boundary ref states between     *
 * neighboring blocks.                                  *
 *                                                      *
 ********************************************************/
int Send_Boundary_Ref_States(Rte2D_Quad_Block *Soln_ptr,
			     AdaptiveBlock2D_List &Soln_Block_List,
			     const int Number_of_Solution_Variables) {

    int error_flag;

    /* Load message buffers at block interfaces with no cell resolution change. */

    error_flag = Load_Send_Message_Buffers_NoResChange(Soln_ptr,
                                                       Soln_Block_List,
                                                       Number_of_Solution_Variables);
    if (error_flag) {
       cout << "\n " << CFDkit_Version() 
            << " Message Passing Error: Load_Send_Message_Buffers_NoResChange, "
            << "flag = " << error_flag << ".\n";
       return(error_flag);
    } /* endif */

    /* Load message buffers at block interfaces for sending from fine to coarse blocks. */

    error_flag = Load_Send_Message_Buffers_ResChange_FineToCoarse(Soln_ptr,
                                                                  Soln_Block_List,
                                                                  Number_of_Solution_Variables);
   if (error_flag) {
       cout << "\n " << CFDkit_Version() 
            << " Message Passing Error: Load_Send_Message_Buffers_ResChange_FineToCoarse, "
            << "flag = " << error_flag << ".\n";
       return(error_flag);
    } /* endif */

    /* Exchange message buffers at block interfaces 
       with no cell resolution change. */
    error_flag = Exchange_Messages_NoResChange(Soln_Block_List,
                                               Number_of_Solution_Variables);
    if (error_flag) {
       cout << "\n " << CFDkit_Version() 
            << " Message Passing Error: Exchange_Messages_NoResChange.\n";
       return(error_flag);
    } /* endif */

    /* Exchange message buffers at block interfaces, 
       sending messages from fine to coarse blocks. */

    error_flag = Exchange_Messages_ResChange_FineToCoarse(Soln_Block_List,
                                                          Number_of_Solution_Variables);
    if (error_flag) {
       cout << "\n " << CFDkit_Version() 
            << " Message Passing Error: Exchange_Messages_ResChange_FineToCoarse, "
            << "flag = " << error_flag << ".\n";
       return(error_flag);
    } /* endif */

    /* Unload message buffers at block interfaces with no cell resolution change. */

    error_flag = Unload_Receive_Message_Buffers_NoResChange(Soln_ptr,
                                                            Soln_Block_List,
                                                            Number_of_Solution_Variables);
    if (error_flag) {
       cout << "\n " << CFDkit_Version() 
            << " Message Passing Error: Unload_Receive_Message_Buffers_NoResChange, "
            << "flag = " << error_flag << ".\n";
       return(error_flag);
    } /* endif */

    /* Unload message buffers at block interfaces as sent from fine to coarse blocks. */

    error_flag = Unload_Receive_Message_Buffers_ResChange_FineToCoarse(Soln_ptr,
                                                                       Soln_Block_List,
                                                                       Number_of_Solution_Variables);
    if (error_flag) {
       cout << "\n " << CFDkit_Version() 
            << " Message Passing Error: Unload_Receive_Message_Buffers_ResChange_FineToCoarse, "
            << "flag = " << error_flag << ".\n";
       return(error_flag);
    } /* endif */

    /* Load message buffers at block interfaces for sending from coarse to fine blocks. */

    error_flag = Load_Send_Message_Buffers_ResChange_CoarseToFine(Soln_ptr,
                                                                  Soln_Block_List,
                                                                  Number_of_Solution_Variables);
    if (error_flag) {
       cout << "\n " << CFDkit_Version() 
            << " Message Passing Error: Load_Send_Message_Buffers_ResChange_CoarseToFine, "
            << "flag = " << error_flag << ".\n";
       return(error_flag);
    } /* endif */

    /* Exchange message buffers at block interfaces, 
       sending messages from coarse to fine blocks. */

    error_flag = Exchange_Messages_ResChange_CoarseToFine(Soln_Block_List,
                                                          Number_of_Solution_Variables);
    if (error_flag) {
       cout << "\n " << CFDkit_Version() 
            << " Message Passing Error: Exchange_Messages_ResChange_CoarseToFine, "
            << "flag = " << error_flag << ".\n";
       return(error_flag);
    } /* endif */

    /* Unload message buffers at block interfaces as sent from coarse to fine blocks. */

    error_flag = Unload_Receive_Message_Buffers_ResChange_CoarseToFine(Soln_ptr,
                                                                       Soln_Block_List,
                                                                       Number_of_Solution_Variables);
   if (error_flag) {
       cout << "\n " << CFDkit_Version() 
            << " Message Passing Error: Unload_Receive_Message_Buffers_ResChange_CoarseToFine, "
            << "flag = " << error_flag << ".\n";
       return(error_flag);
    } /* endif */


    /* Return error flag. */

    return(error_flag);

}



#endif /* _RTE2D_AMR_INCLUDED  */
