/* AdaptiveBlock2D_MessagePassing.h:  Templated subroutines for adaptive blocks classes
                                      which load and unload message passing buffers. */

#ifndef _ADAPTIVEBLOCK2D_MESSAGEPASSING_INCLUDED
#define _ADAPTIVEBLOCK2D_MESSAGEPASSING_INCLUDED

/* Include adaptive block and 2D vector header files. */

#ifndef _ADAPTIVBLOCK_INCLUDED
#include "AdaptiveBlock.h"
#endif // _ADAPTIVEBLOCK_INCLUDED

#ifndef _VECTOR2D_INCLUDED
#include "../Math/Vector2D.h"
#endif //_VECTOR2D_INCLUDED

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
template <class Quad_Soln_Block>
int Load_Send_Message_Buffers_NoResChange(Quad_Soln_Block *Soln_ptr,
                                          AdaptiveBlock2D_List &Soln_Block_List,
                                          const int Number_of_Solution_Variables,
                                          const int Send_Mesh_Geometry_Only) {

    int i_blk, buffer_size_neighbour, 
        i_min, i_max, i_inc, i, 
        j_min, j_max, j_inc, j, 
        k, l;
    int i_ref, j_ref;
    Vector2D x_ref;

    /* Check to see if the number of solution variables specified
       in the calling argument is correct. */

    if (Number_of_Solution_Variables != NUM_COMP_VECTOR2D &&
        Number_of_Solution_Variables != Soln_ptr[0].NumVar() &&
        Number_of_Solution_Variables != Soln_ptr[0].NumVar() + NUM_COMP_VECTOR2D) {
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
          if (!Send_Mesh_Geometry_Only) {
             buffer_size_neighbour = Soln_Block_List.Block[i_blk].infoN[0].dimen.ghost*
                                     abs(Soln_Block_List.Block[i_blk].infoN[0].dimen.i)*
                                     Soln_ptr[i_blk].NumVar();
             if (Number_of_Solution_Variables > Soln_ptr[i_blk].NumVar()) 
                buffer_size_neighbour = buffer_size_neighbour +
                                        Soln_Block_List.Block[i_blk].infoN[0].dimen.ghost*
                                        (abs(Soln_Block_List.Block[i_blk].infoN[0].dimen.i)+1)*
                                        NUM_COMP_VECTOR2D+
                                        2*Soln_Block_List.Block[i_blk].infoN[0].dimen.ghost;
          } else {
                buffer_size_neighbour = Soln_Block_List.Block[i_blk].infoN[0].dimen.ghost*
                                        (abs(Soln_Block_List.Block[i_blk].infoN[0].dimen.i)+1)*
                                        NUM_COMP_VECTOR2D+
                                        2*Soln_Block_List.Block[i_blk].infoN[0].dimen.ghost;
          } /* endif */
          l = -1;
          // Load ghost cell solution variable information as required.
          if (!Send_Mesh_Geometry_Only) {
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
             i = Soln_ptr[i_blk].LoadSendBuffer(Soln_Block_List.message_noreschange_northface_sendbuf[i_blk],
                                                l,buffer_size_neighbour,
                                                i_min,i_max,i_inc,
                                                j_min,j_max,j_inc);
	     if (i != 0) return(2200);
          } /* endif */
          // Load ghost cell mesh information as required.
          if (Send_Mesh_Geometry_Only ||
              Number_of_Solution_Variables > Soln_ptr[i_blk].NumVar()) {
             if (Soln_Block_List.Block[i_blk].infoN[0].dimen.i > 0 &&
                 Soln_Block_List.Block[i_blk].infoN[0].dimen.j > 0) {
	        i_min = Soln_ptr[i_blk].Grid.INl;
	        i_max = Soln_ptr[i_blk].Grid.INu;
	        i_inc = 1;
	        j_min = Soln_ptr[i_blk].Grid.JNu-Soln_Block_List.Block[i_blk].infoN[0].dimen.ghost;
	        j_max = Soln_ptr[i_blk].Grid.JNu-1;
	        j_inc = 1;
                i_ref = i_min;
                j_ref = j_max+1;
             } else if (Soln_Block_List.Block[i_blk].infoN[0].dimen.j > 0) {
	        i_min = Soln_ptr[i_blk].Grid.INu;
	        i_max = Soln_ptr[i_blk].Grid.INl;
	        i_inc = -1;
	        j_min = Soln_ptr[i_blk].Grid.JNu-Soln_Block_List.Block[i_blk].infoN[0].dimen.ghost;
	        j_max = Soln_ptr[i_blk].Grid.JNu-1;
	        j_inc = 1;
                i_ref = i_min;
                j_ref = j_max+1;
             } else if (Soln_Block_List.Block[i_blk].infoN[0].dimen.i > 0) {
	        i_min = Soln_ptr[i_blk].Grid.INl;
	        i_max = Soln_ptr[i_blk].Grid.INu;
	        i_inc = 1;
	        j_min = Soln_ptr[i_blk].Grid.JNu-1;
	        j_max = Soln_ptr[i_blk].Grid.JNu-Soln_Block_List.Block[i_blk].infoN[0].dimen.ghost;
	        j_inc = -1;
                i_ref = i_min;
                j_ref = j_min+1;
             } else {
	        i_min = Soln_ptr[i_blk].Grid.INu;
	        i_max = Soln_ptr[i_blk].Grid.INl;
	        i_inc = -1;
	        j_min = Soln_ptr[i_blk].Grid.JNu-1;
	        j_max = Soln_ptr[i_blk].Grid.JNu-Soln_Block_List.Block[i_blk].infoN[0].dimen.ghost;
	        j_inc = -1;
                i_ref = i_min;
                j_ref = j_min+1;
             } /* endif */
             x_ref = Soln_ptr[i_blk].Grid.Node[i_ref][j_ref].X; // Reference node location.
             for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc ) {
                for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
                   for ( k = 1 ; k <= NUM_COMP_VECTOR2D; ++ k) {
		      l = l + 1;
                      if (l >= buffer_size_neighbour) return(2201);
                      Soln_Block_List.message_noreschange_northface_sendbuf[i_blk][l] =
                         Soln_ptr[i_blk].Grid.Node[i][j].X[k]-x_ref[k];
                   } /* endfor */
                } /* endfor */
             } /* endfor */
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
             if (i_inc > 0) {
                for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc ) {
                   l = l + 1;
                   if (l >= buffer_size_neighbour) return(2202);
                   Soln_Block_List.message_noreschange_northface_sendbuf[i_blk][l] = 
                      double(Soln_ptr[i_blk].Grid.BCtypeW[j]);
                } /* endfor */
                for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc ) {
                   l = l + 1;
                   if (l >= buffer_size_neighbour) return(2203);
                   Soln_Block_List.message_noreschange_northface_sendbuf[i_blk][l] = 
                      double(Soln_ptr[i_blk].Grid.BCtypeE[j]);
                } /* endfor */
	     } else {
                for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc ) {
                   l = l + 1;
                   if (l >= buffer_size_neighbour) return(2204);
                   Soln_Block_List.message_noreschange_northface_sendbuf[i_blk][l] = 
                      double(Soln_ptr[i_blk].Grid.BCtypeE[j]);
                } /* endfor */
                for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc ) {
                   l = l + 1;
                   if (l >= buffer_size_neighbour) return(2205);
                   Soln_Block_List.message_noreschange_northface_sendbuf[i_blk][l] = 
                      double(Soln_ptr[i_blk].Grid.BCtypeW[j]);
                } /* endfor */
             } /* endif */
          } /* endif */
       } /* endif */

       ////////////////////////////////////////////////////////////////////
       // Load send buffer for South face neighbours of solution blocks. //
       ////////////////////////////////////////////////////////////////////
       if (Soln_Block_List.Block[i_blk].used && 
           (Soln_Block_List.Block[i_blk].nS == 1) &&
           (Soln_Block_List.Block[i_blk].info.level == 
            Soln_Block_List.Block[i_blk].infoS[0].level)) {
          if (!Send_Mesh_Geometry_Only) {
             buffer_size_neighbour = Soln_Block_List.Block[i_blk].infoS[0].dimen.ghost*
                                     abs(Soln_Block_List.Block[i_blk].infoS[0].dimen.i)*
                                     Soln_ptr[i_blk].NumVar();
             if (Number_of_Solution_Variables > Soln_ptr[i_blk].NumVar()) 
                buffer_size_neighbour = buffer_size_neighbour +
                                        Soln_Block_List.Block[i_blk].infoS[0].dimen.ghost*
                                        (abs(Soln_Block_List.Block[i_blk].infoS[0].dimen.i)+1)*
                                        NUM_COMP_VECTOR2D+
                                        2*Soln_Block_List.Block[i_blk].infoS[0].dimen.ghost;
          } else {
                buffer_size_neighbour = Soln_Block_List.Block[i_blk].infoS[0].dimen.ghost*
                                        (abs(Soln_Block_List.Block[i_blk].infoS[0].dimen.i)+1)*
                                        NUM_COMP_VECTOR2D+
                                        2*Soln_Block_List.Block[i_blk].infoS[0].dimen.ghost;
          } /* endif */
          l = -1;
          // Load ghost cell solution variable information as required.
          if (!Send_Mesh_Geometry_Only) {
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
             i = Soln_ptr[i_blk].LoadSendBuffer(Soln_Block_List.message_noreschange_southface_sendbuf[i_blk],
                                                l,buffer_size_neighbour,
                                                i_min,i_max,i_inc,
                                                j_min,j_max,j_inc);
	     if (i != 0) return(2100);
          } /* endif */
          // Load ghost cell mesh information as required.
          if (Send_Mesh_Geometry_Only ||
              Number_of_Solution_Variables > Soln_ptr[i_blk].NumVar()) {
             if (Soln_Block_List.Block[i_blk].infoS[0].dimen.i > 0 &&
                 Soln_Block_List.Block[i_blk].infoS[0].dimen.j > 0) {
	        i_min = Soln_ptr[i_blk].Grid.INl;
	        i_max = Soln_ptr[i_blk].Grid.INu;
                i_inc = 1;
	        j_min = Soln_ptr[i_blk].Grid.JNl+1;
	        j_max = Soln_ptr[i_blk].Grid.JNl+Soln_Block_List.Block[i_blk].infoS[0].dimen.ghost;
	        j_inc = 1;
                i_ref = i_min;
                j_ref = j_min-1;
             } else if (Soln_Block_List.Block[i_blk].infoS[0].dimen.j > 0) {
	        i_min = Soln_ptr[i_blk].Grid.INu;
	        i_max = Soln_ptr[i_blk].Grid.INl;
                i_inc = -1;
	        j_min = Soln_ptr[i_blk].Grid.JNl+1;
	        j_max = Soln_ptr[i_blk].Grid.JNl+Soln_Block_List.Block[i_blk].infoS[0].dimen.ghost;
	        j_inc = 1;
                i_ref = i_min;
                j_ref = j_min-1;
             } else if (Soln_Block_List.Block[i_blk].infoS[0].dimen.i > 0) {
	        i_min = Soln_ptr[i_blk].Grid.INl;
	        i_max = Soln_ptr[i_blk].Grid.INu;
                i_inc = 1;
	        j_min = Soln_ptr[i_blk].Grid.JNl+Soln_Block_List.Block[i_blk].infoS[0].dimen.ghost;
	        j_max = Soln_ptr[i_blk].Grid.JNl+1;
	        j_inc = -1;
                i_ref = i_min;
                j_ref = j_max-1;
             } else {
	        i_min = Soln_ptr[i_blk].Grid.INu;
	        i_max = Soln_ptr[i_blk].Grid.INl;
                i_inc = -1;
	        j_min = Soln_ptr[i_blk].Grid.JNl+Soln_Block_List.Block[i_blk].infoS[0].dimen.ghost;
	        j_max = Soln_ptr[i_blk].Grid.JNl+1;
	        j_inc = -1;
                i_ref = i_min;
                j_ref = j_max-1;
             } /* endif */
             x_ref = Soln_ptr[i_blk].Grid.Node[i_ref][j_ref].X; // Reference node location.
             for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc ) {
                for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
                   for ( k = 1 ; k <= NUM_COMP_VECTOR2D; ++ k) {
		      l = l + 1;
                      if (l >= buffer_size_neighbour) return(2101);
                      Soln_Block_List.message_noreschange_southface_sendbuf[i_blk][l] =
                         Soln_ptr[i_blk].Grid.Node[i][j].X[k]-x_ref[k];
                   } /* endfor */
                } /* endfor */
             } /* endfor */
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
             if (i_inc > 0) {
                for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc ) {
                   l = l + 1;
                   if (l >= buffer_size_neighbour) return(2102);
                   Soln_Block_List.message_noreschange_southface_sendbuf[i_blk][l] = 
                      double(Soln_ptr[i_blk].Grid.BCtypeW[j]);
                } /* endfor */
                for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc ) {
                   l = l + 1;
                   if (l >= buffer_size_neighbour) return(2103);
                   Soln_Block_List.message_noreschange_southface_sendbuf[i_blk][l] = 
                      double(Soln_ptr[i_blk].Grid.BCtypeE[j]);
                } /* endfor */
	     } else {
                for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc ) {
                   l = l + 1;
                   if (l >= buffer_size_neighbour) return(2104);
                   Soln_Block_List.message_noreschange_southface_sendbuf[i_blk][l] = 
                      double(Soln_ptr[i_blk].Grid.BCtypeE[j]);
                } /* endfor */
                for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc ) {
                   l = l + 1;
                   if (l >= buffer_size_neighbour) return(2105);
                   Soln_Block_List.message_noreschange_southface_sendbuf[i_blk][l] = 
                      double(Soln_ptr[i_blk].Grid.BCtypeW[j]);
                } /* endfor */
             } /* endif */
          } /* endif */
       } /* endif */

       ///////////////////////////////////////////////////////////////////
       // Load send buffer for East face neighbours of solution blocks. //
       ///////////////////////////////////////////////////////////////////
       if (Soln_Block_List.Block[i_blk].used && 
           (Soln_Block_List.Block[i_blk].nE == 1) &&
           (Soln_Block_List.Block[i_blk].info.level == 
            Soln_Block_List.Block[i_blk].infoE[0].level)) {
          if (!Send_Mesh_Geometry_Only) {
             buffer_size_neighbour = Soln_Block_List.Block[i_blk].infoE[0].dimen.ghost*
                                     abs(Soln_Block_List.Block[i_blk].infoE[0].dimen.j)*
                                     Soln_ptr[i_blk].NumVar();
             if (Number_of_Solution_Variables > Soln_ptr[i_blk].NumVar()) 
                buffer_size_neighbour = buffer_size_neighbour +
                                        Soln_Block_List.Block[i_blk].infoE[0].dimen.ghost*
                                        (abs(Soln_Block_List.Block[i_blk].infoE[0].dimen.j)+1)*
                                        NUM_COMP_VECTOR2D+
                                        2*Soln_Block_List.Block[i_blk].infoE[0].dimen.ghost;
          } else {
                buffer_size_neighbour = Soln_Block_List.Block[i_blk].infoE[0].dimen.ghost*
                                        (abs(Soln_Block_List.Block[i_blk].infoE[0].dimen.j)+1)*
                                        NUM_COMP_VECTOR2D+
                                        2*Soln_Block_List.Block[i_blk].infoE[0].dimen.ghost;
          } /* endif */
          l = -1;
          // Load ghost cell solution variable information as required.
          if (!Send_Mesh_Geometry_Only) {
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
             i = Soln_ptr[i_blk].LoadSendBuffer(Soln_Block_List.message_noreschange_eastface_sendbuf[i_blk],
                                                l,buffer_size_neighbour,
                                                i_min,i_max,i_inc,
                                                j_min,j_max,j_inc);
	     if (i != 0) return(1200);
          } /* endif */
          // Load ghost cell mesh information as required.
          if (Send_Mesh_Geometry_Only ||
              Number_of_Solution_Variables > Soln_ptr[i_blk].NumVar()) {
             if (Soln_Block_List.Block[i_blk].infoE[0].dimen.i > 0 &&
                 Soln_Block_List.Block[i_blk].infoE[0].dimen.j > 0) {
	        i_min = Soln_ptr[i_blk].Grid.INu-Soln_Block_List.Block[i_blk].infoE[0].dimen.ghost;
	        i_max = Soln_ptr[i_blk].Grid.INu-1;
	        i_inc = 1;
	        j_min = Soln_ptr[i_blk].Grid.JNl;
	        j_max = Soln_ptr[i_blk].Grid.JNu;
	        j_inc = 1;
                i_ref = i_max+1;
                j_ref = j_min;
             } else if (Soln_Block_List.Block[i_blk].infoE[0].dimen.j > 0) {
	        i_min = Soln_ptr[i_blk].Grid.INu-1;
	        i_max = Soln_ptr[i_blk].Grid.INu-Soln_Block_List.Block[i_blk].infoE[0].dimen.ghost;
	        i_inc = -1;
	        j_min = Soln_ptr[i_blk].Grid.JNl;
	        j_max = Soln_ptr[i_blk].Grid.JNu;
	        j_inc = 1;
                i_ref = i_min+1;
                j_ref = j_min;
             } else if (Soln_Block_List.Block[i_blk].infoE[0].dimen.i > 0) {
	        i_min = Soln_ptr[i_blk].Grid.INu-Soln_Block_List.Block[i_blk].infoE[0].dimen.ghost;
	        i_max = Soln_ptr[i_blk].Grid.INu-1;
	        i_inc = 1;
	        j_min = Soln_ptr[i_blk].Grid.JNu;
	        j_max = Soln_ptr[i_blk].Grid.JNl;
	        j_inc = -1;
                i_ref = i_max+1;
                j_ref = j_min;
             } else {
	        i_min = Soln_ptr[i_blk].Grid.INu-1;
	        i_max = Soln_ptr[i_blk].Grid.INu-Soln_Block_List.Block[i_blk].infoE[0].dimen.ghost;
	        i_inc = -1;
	        j_min = Soln_ptr[i_blk].Grid.JNu;
	        j_max = Soln_ptr[i_blk].Grid.JNl;
	        j_inc = -1;
                i_ref = i_min+1;
                j_ref = j_min;
             } /* endif */
             x_ref = Soln_ptr[i_blk].Grid.Node[i_ref][j_ref].X; // Reference node location.
             for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc ) {
                for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
                   for ( k = 1 ; k <= NUM_COMP_VECTOR2D; ++ k) {
		      l = l + 1;
                      if (l >= buffer_size_neighbour) return(1201);
                      Soln_Block_List.message_noreschange_eastface_sendbuf[i_blk][l] =
                         Soln_ptr[i_blk].Grid.Node[i][j].X[k]-x_ref[k];
                   } /* endfor */
                } /* endfor */
             } /* endfor */
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
             if (j_inc > 0) {
                for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
                   l = l + 1;
                   if (l >= buffer_size_neighbour) return(1202);
                   Soln_Block_List.message_noreschange_eastface_sendbuf[i_blk][l] = 
                      double(Soln_ptr[i_blk].Grid.BCtypeS[i]);
                } /* endfor */
                for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
                   l = l + 1;
                   if (l >= buffer_size_neighbour) return(1203);
                   Soln_Block_List.message_noreschange_eastface_sendbuf[i_blk][l] = 
                      double(Soln_ptr[i_blk].Grid.BCtypeN[i]);
                } /* endfor */
	     } else {
                for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
                   l = l + 1;
                   if (l >= buffer_size_neighbour) return(1204);
                   Soln_Block_List.message_noreschange_eastface_sendbuf[i_blk][l] = 
                      double(Soln_ptr[i_blk].Grid.BCtypeN[i]);
                } /* endfor */
                for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
                   l = l + 1;
                   if (l >= buffer_size_neighbour) return(1205);
                   Soln_Block_List.message_noreschange_eastface_sendbuf[i_blk][l] = 
                      double(Soln_ptr[i_blk].Grid.BCtypeS[i]);
                } /* endfor */
             } /* endif */
          } /* endif */
       } /* endif */

       ///////////////////////////////////////////////////////////////////
       // Load send buffer for West face neighbours of solution blocks. //
       ///////////////////////////////////////////////////////////////////
       if (Soln_Block_List.Block[i_blk].used && 
           (Soln_Block_List.Block[i_blk].nW == 1) &&
           (Soln_Block_List.Block[i_blk].info.level == 
            Soln_Block_List.Block[i_blk].infoW[0].level)) {
          if (!Send_Mesh_Geometry_Only) {
             buffer_size_neighbour = Soln_Block_List.Block[i_blk].infoW[0].dimen.ghost*
                                     abs(Soln_Block_List.Block[i_blk].infoW[0].dimen.j)*
                                     Soln_ptr[i_blk].NumVar();
             if (Number_of_Solution_Variables > Soln_ptr[i_blk].NumVar()) 
                buffer_size_neighbour = buffer_size_neighbour +
                                        Soln_Block_List.Block[i_blk].infoW[0].dimen.ghost*
                                        (abs(Soln_Block_List.Block[i_blk].infoW[0].dimen.j)+1)*
                                        NUM_COMP_VECTOR2D+
                                        2*Soln_Block_List.Block[i_blk].infoW[0].dimen.ghost;
          } else {
                buffer_size_neighbour = Soln_Block_List.Block[i_blk].infoW[0].dimen.ghost*
                                        (abs(Soln_Block_List.Block[i_blk].infoW[0].dimen.j)+1)*
                                        NUM_COMP_VECTOR2D+
                                        2*Soln_Block_List.Block[i_blk].infoW[0].dimen.ghost;
          } /* endif */
          l = -1;
          // Load ghost cell solution variable information as required.
          if (!Send_Mesh_Geometry_Only) {
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
             i = Soln_ptr[i_blk].LoadSendBuffer(Soln_Block_List.message_noreschange_westface_sendbuf[i_blk],
                                                l,buffer_size_neighbour,
                                                i_min,i_max,i_inc,
                                                j_min,j_max,j_inc);
	     if (i != 0) return(1100);
          } /* endif */
          // Load ghost cell mesh information as required.
          if (Send_Mesh_Geometry_Only ||
              Number_of_Solution_Variables > Soln_ptr[i_blk].NumVar()) {
             if (Soln_Block_List.Block[i_blk].infoW[0].dimen.i > 0 &&
                 Soln_Block_List.Block[i_blk].infoW[0].dimen.j > 0) {
	        i_min = Soln_ptr[i_blk].Grid.INl+1;
	        i_max = Soln_ptr[i_blk].Grid.INl+Soln_Block_List.Block[i_blk].infoW[0].dimen.ghost;
	        i_inc = 1;
	        j_min = Soln_ptr[i_blk].Grid.JNl;
	        j_max = Soln_ptr[i_blk].Grid.JNu;
	        j_inc = 1;
                i_ref = i_min-1;
                j_ref = j_min;
             } else if (Soln_Block_List.Block[i_blk].infoW[0].dimen.j > 0) {
	        i_min = Soln_ptr[i_blk].Grid.INl+Soln_Block_List.Block[i_blk].infoW[0].dimen.ghost;
	        i_max = Soln_ptr[i_blk].Grid.INl+1;
	        i_inc = -1;
	        j_min = Soln_ptr[i_blk].Grid.JNl;
	        j_max = Soln_ptr[i_blk].Grid.JNu;
	        j_inc = 1;
                i_ref = i_max-1;
                j_ref = j_min;
             } else if (Soln_Block_List.Block[i_blk].infoW[0].dimen.i > 0) {
	        i_min = Soln_ptr[i_blk].Grid.INl+1;
	        i_max = Soln_ptr[i_blk].Grid.INl+Soln_Block_List.Block[i_blk].infoW[0].dimen.ghost;
	        i_inc = 1;
	        j_min = Soln_ptr[i_blk].Grid.JNu;
	        j_max = Soln_ptr[i_blk].Grid.JNl;
	        j_inc = -1;
                i_ref = i_min-1;
                j_ref = j_min;
             } else {
	        i_min = Soln_ptr[i_blk].Grid.INl+Soln_Block_List.Block[i_blk].infoW[0].dimen.ghost;
	        i_max = Soln_ptr[i_blk].Grid.INl+1;
	        i_inc = -1;
	        j_min = Soln_ptr[i_blk].Grid.JNu;
	        j_max = Soln_ptr[i_blk].Grid.JNl;
	        j_inc = -1;
                i_ref = i_max-1;
                j_ref = j_min;
             } /* endif */
             x_ref = Soln_ptr[i_blk].Grid.Node[i_ref][j_ref].X; // Reference node location.
             for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc ) {
                for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
                   for ( k = 1 ; k <= NUM_COMP_VECTOR2D; ++ k) {
		      l = l + 1;
                      if (l >= buffer_size_neighbour) return(1101);
                      Soln_Block_List.message_noreschange_westface_sendbuf[i_blk][l] =
                         Soln_ptr[i_blk].Grid.Node[i][j].X[k]-x_ref[k];
                   } /* endfor */
                } /* endfor */
             } /* endfor */
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
             if (j_inc > 0) {
                for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
                   l = l + 1;
                   if (l >= buffer_size_neighbour) return(1102);
                   Soln_Block_List.message_noreschange_westface_sendbuf[i_blk][l] = 
                      double(Soln_ptr[i_blk].Grid.BCtypeS[i]);
                } /* endfor */
                for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
                   l = l + 1;
                   if (l >= buffer_size_neighbour) return(1103);
                   Soln_Block_List.message_noreschange_westface_sendbuf[i_blk][l] = 
                      double(Soln_ptr[i_blk].Grid.BCtypeN[i]);
                } /* endfor */
	     } else {
                for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
                   l = l + 1;
                   if (l >= buffer_size_neighbour) return(1104);
                   Soln_Block_List.message_noreschange_westface_sendbuf[i_blk][l] = 
                      double(Soln_ptr[i_blk].Grid.BCtypeN[i]);
                } /* endfor */
                for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
                   l = l + 1;
                   if (l >= buffer_size_neighbour) return(1105);
                   Soln_Block_List.message_noreschange_westface_sendbuf[i_blk][l] = 
                      double(Soln_ptr[i_blk].Grid.BCtypeS[i]);
                } /* endfor */
             } /* endif */
          } /* endif */
       } /* endif */

       ///////////////////////////////////////////////////////////////////////////
       // Load send buffer for North West corner neighbours of solution blocks. //
       ///////////////////////////////////////////////////////////////////////////
       if (Soln_Block_List.Block[i_blk].used && 
           (Soln_Block_List.Block[i_blk].nNW == 1) &&
           (Soln_Block_List.Block[i_blk].info.level == 
            Soln_Block_List.Block[i_blk].infoNW[0].level)) {
          buffer_size_neighbour = sqr(Soln_Block_List.Block[i_blk].infoNW[0].dimen.ghost)*
                                  Number_of_Solution_Variables;
          l = -1;
          // Load ghost cell solution variable information as required.
          if (!Send_Mesh_Geometry_Only) {
             if (Soln_Block_List.Block[i_blk].infoNW[0].dimen.i > 0 &&
                 Soln_Block_List.Block[i_blk].infoNW[0].dimen.j > 0) {
	        i_min = Soln_ptr[i_blk].ICl;
                i_max = Soln_ptr[i_blk].ICl+Soln_Block_List.Block[i_blk].infoNW[0].dimen.ghost-1;
	        i_inc = 1;
	        j_min = Soln_ptr[i_blk].JCu-Soln_Block_List.Block[i_blk].infoNW[0].dimen.ghost+1;
	        j_max = Soln_ptr[i_blk].JCu;
	        j_inc = 1;
             } else if (Soln_Block_List.Block[i_blk].infoNW[0].dimen.j > 0) {
	        i_min = Soln_ptr[i_blk].ICl+Soln_Block_List.Block[i_blk].infoNW[0].dimen.ghost-1;
	        i_max = Soln_ptr[i_blk].ICl;
	        i_inc = -1;
	        j_min = Soln_ptr[i_blk].JCu-Soln_Block_List.Block[i_blk].infoNW[0].dimen.ghost+1;
	        j_max = Soln_ptr[i_blk].JCu;
	        j_inc = 1;
             } else if (Soln_Block_List.Block[i_blk].infoNW[0].dimen.i > 0) {
	        i_min = Soln_ptr[i_blk].ICl;
	        i_max = Soln_ptr[i_blk].ICl+Soln_Block_List.Block[i_blk].infoNW[0].dimen.ghost-1;
	        i_inc = 1;
	        j_min = Soln_ptr[i_blk].JCu;
	        j_max = Soln_ptr[i_blk].JCu-Soln_Block_List.Block[i_blk].infoNW[0].dimen.ghost+1;
	        j_inc = -1;
             } else {
	        i_min = Soln_ptr[i_blk].ICl+Soln_Block_List.Block[i_blk].infoNW[0].dimen.ghost-1;
	        i_max = Soln_ptr[i_blk].ICl;
	        i_inc = -1;
	        j_min = Soln_ptr[i_blk].JCu;
 	        j_max = Soln_ptr[i_blk].JCu-Soln_Block_List.Block[i_blk].infoNW[0].dimen.ghost+1;
	        j_inc = -1;
             } /* endif */
             i = Soln_ptr[i_blk].LoadSendBuffer(Soln_Block_List.message_noreschange_northwestcorner_sendbuf[i_blk],
                                                l,buffer_size_neighbour,
                                                i_min,i_max,i_inc,
                                                j_min,j_max,j_inc);
	     if (i != 0) return(4200);
          } /* endif */
          // Load ghost cell mesh information as required.
          if (Send_Mesh_Geometry_Only ||
              Number_of_Solution_Variables > Soln_ptr[i_blk].NumVar()) {
             if (Soln_Block_List.Block[i_blk].infoNW[0].dimen.i > 0 &&
                 Soln_Block_List.Block[i_blk].infoNW[0].dimen.j > 0) {
	        i_min = Soln_ptr[i_blk].Grid.INl+1;
	        i_max = Soln_ptr[i_blk].Grid.INl+Soln_Block_List.Block[i_blk].infoNW[0].dimen.ghost;
	        i_inc = 1;
	        j_min = Soln_ptr[i_blk].Grid.JNu-Soln_Block_List.Block[i_blk].infoNW[0].dimen.ghost;
	        j_max = Soln_ptr[i_blk].Grid.JNu-1;
	        j_inc = 1;
                i_ref = i_min-1;
                j_ref = j_max+1;
             } else if (Soln_Block_List.Block[i_blk].infoNW[0].dimen.j > 0) {
	        i_min = Soln_ptr[i_blk].Grid.INl+Soln_Block_List.Block[i_blk].infoNW[0].dimen.ghost;
	        i_max = Soln_ptr[i_blk].Grid.INl+1;
	        i_inc = -1;
	        j_min = Soln_ptr[i_blk].Grid.JNu-Soln_Block_List.Block[i_blk].infoNW[0].dimen.ghost;
	        j_max = Soln_ptr[i_blk].Grid.JNu-1;
	        j_inc = 1;
                i_ref = i_max-1;
                j_ref = j_max+1;
             } else if (Soln_Block_List.Block[i_blk].infoNW[0].dimen.i > 0) {
	        i_min = Soln_ptr[i_blk].Grid.INl+1;
	        i_max = Soln_ptr[i_blk].Grid.INl+Soln_Block_List.Block[i_blk].infoNW[0].dimen.ghost;
	        i_inc = 1;
	        j_min = Soln_ptr[i_blk].Grid.JNu-1;
	        j_max = Soln_ptr[i_blk].Grid.JNu-Soln_Block_List.Block[i_blk].infoNW[0].dimen.ghost;
	        j_inc = -1;
                i_ref = i_min-1;
                j_ref = j_min+1;
             } else {
	        i_min = Soln_ptr[i_blk].Grid.INl+Soln_Block_List.Block[i_blk].infoNW[0].dimen.ghost;
	        i_max = Soln_ptr[i_blk].Grid.INl+1;
	        i_inc = -1;
	        j_min = Soln_ptr[i_blk].Grid.JNu-1;
	        j_max = Soln_ptr[i_blk].Grid.JNu-Soln_Block_List.Block[i_blk].infoNW[0].dimen.ghost;
	        j_inc = -1;
                i_ref = i_max-1;
                j_ref = j_min+1;
             } /* endif */
             x_ref = Soln_ptr[i_blk].Grid.Node[i_ref][j_ref].X; // Reference node location.
             for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc ) {
                for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
                   for ( k = 1 ; k <= NUM_COMP_VECTOR2D; ++ k) {
		      l = l + 1;
                      if (l >= buffer_size_neighbour) return(4201);
                      Soln_Block_List.message_noreschange_northwestcorner_sendbuf[i_blk][l] =
                         Soln_ptr[i_blk].Grid.Node[i][j].X[k]-x_ref[k];
                   } /* endfor */
                } /* endfor */
             } /* endfor */
          } /* endif */
       } /* endif */

       ///////////////////////////////////////////////////////////////////////////
       // Load send buffer for North East corner neighbours of solution blocks. //
       ///////////////////////////////////////////////////////////////////////////
       if (Soln_Block_List.Block[i_blk].used && 
           (Soln_Block_List.Block[i_blk].nNE == 1) &&
           (Soln_Block_List.Block[i_blk].info.level == 
            Soln_Block_List.Block[i_blk].infoNE[0].level)) {
          buffer_size_neighbour = sqr(Soln_Block_List.Block[i_blk].infoNE[0].dimen.ghost)*
                                  Number_of_Solution_Variables;
          l = -1;
          // Load ghost cell solution variable information as required.
          if (!Send_Mesh_Geometry_Only) {
             if (Soln_Block_List.Block[i_blk].infoNE[0].dimen.i > 0 &&
                 Soln_Block_List.Block[i_blk].infoNE[0].dimen.j > 0) {
	        i_min = Soln_ptr[i_blk].ICu-Soln_Block_List.Block[i_blk].infoNE[0].dimen.ghost+1;
                i_max = Soln_ptr[i_blk].ICu;
	        i_inc = 1;
	        j_min = Soln_ptr[i_blk].JCu-Soln_Block_List.Block[i_blk].infoNE[0].dimen.ghost+1;
	        j_max = Soln_ptr[i_blk].JCu;
	        j_inc = 1;
             } else if (Soln_Block_List.Block[i_blk].infoNE[0].dimen.j > 0) {
	        i_min = Soln_ptr[i_blk].ICu;
	        i_max = Soln_ptr[i_blk].ICu-Soln_Block_List.Block[i_blk].infoNE[0].dimen.ghost+1;
	        i_inc = -1;
	        j_min = Soln_ptr[i_blk].JCu-Soln_Block_List.Block[i_blk].infoNE[0].dimen.ghost+1;
	        j_max = Soln_ptr[i_blk].JCu;
	        j_inc = 1;
             } else if (Soln_Block_List.Block[i_blk].infoNE[0].dimen.i > 0) {
	        i_min = Soln_ptr[i_blk].ICu-Soln_Block_List.Block[i_blk].infoNE[0].dimen.ghost+1;
	        i_max = Soln_ptr[i_blk].ICu;
	        i_inc = 1;
	        j_min = Soln_ptr[i_blk].JCu;
	        j_max = Soln_ptr[i_blk].JCu-Soln_Block_List.Block[i_blk].infoNE[0].dimen.ghost+1;
	        j_inc = -1;
             } else {
	        i_min = Soln_ptr[i_blk].ICu;
	        i_max = Soln_ptr[i_blk].ICu-Soln_Block_List.Block[i_blk].infoNE[0].dimen.ghost+1;
	        i_inc = -1;
	        j_min = Soln_ptr[i_blk].JCu;
 	        j_max = Soln_ptr[i_blk].JCu-Soln_Block_List.Block[i_blk].infoNE[0].dimen.ghost+1;
	        j_inc = -1;
             } /* endif */
             i = Soln_ptr[i_blk].LoadSendBuffer(Soln_Block_List.message_noreschange_northeastcorner_sendbuf[i_blk],
                                                l,buffer_size_neighbour,
                                                i_min,i_max,i_inc,
                                                j_min,j_max,j_inc);
	     if (i != 0) return(5200);
          } /* endif */
          // Load ghost cell mesh information as required.
          if (Send_Mesh_Geometry_Only ||
              Number_of_Solution_Variables > Soln_ptr[i_blk].NumVar()) {
             if (Soln_Block_List.Block[i_blk].infoNE[0].dimen.i > 0 &&
                 Soln_Block_List.Block[i_blk].infoNE[0].dimen.j > 0) {
	        i_min = Soln_ptr[i_blk].Grid.INu-Soln_Block_List.Block[i_blk].infoNE[0].dimen.ghost;
	        i_max = Soln_ptr[i_blk].Grid.INu-1;
	        i_inc = 1;
	        j_min = Soln_ptr[i_blk].Grid.JNu-Soln_Block_List.Block[i_blk].infoNE[0].dimen.ghost;
	        j_max = Soln_ptr[i_blk].Grid.JNu-1;
	        j_inc = 1;
                i_ref = i_max+1;
                j_ref = j_max+1;
             } else if (Soln_Block_List.Block[i_blk].infoNE[0].dimen.j > 0) {
	        i_min = Soln_ptr[i_blk].Grid.INu-1;
	        i_max = Soln_ptr[i_blk].Grid.INu-Soln_Block_List.Block[i_blk].infoNE[0].dimen.ghost;
	        i_inc = -1;
	        j_min = Soln_ptr[i_blk].Grid.JNu-Soln_Block_List.Block[i_blk].infoNE[0].dimen.ghost;
	        j_max = Soln_ptr[i_blk].Grid.JNu-1;
	        j_inc = 1;
                i_ref = i_min+1;
                j_ref = j_max+1;
             } else if (Soln_Block_List.Block[i_blk].infoNE[0].dimen.i > 0) {
	        i_min = Soln_ptr[i_blk].Grid.INu-Soln_Block_List.Block[i_blk].infoNE[0].dimen.ghost;
	        i_max = Soln_ptr[i_blk].Grid.INu-1;
	        i_inc = 1;
	        j_min = Soln_ptr[i_blk].Grid.JNu-1;
	        j_max = Soln_ptr[i_blk].Grid.JNu-Soln_Block_List.Block[i_blk].infoNE[0].dimen.ghost;
	        j_inc = -1;
                i_ref = i_max+1;
                j_ref = j_min+1;
             } else {
	        i_min = Soln_ptr[i_blk].Grid.INu-1;
	        i_max = Soln_ptr[i_blk].Grid.INu-Soln_Block_List.Block[i_blk].infoNE[0].dimen.ghost;
	        i_inc = -1;
	        j_min = Soln_ptr[i_blk].Grid.JNu-1;
	        j_max = Soln_ptr[i_blk].Grid.JNu-Soln_Block_List.Block[i_blk].infoNE[0].dimen.ghost;
	        j_inc = -1;
                i_ref = i_min+1;
                j_ref = j_min+1;
             } /* endif */
             x_ref = Soln_ptr[i_blk].Grid.Node[i_ref][j_ref].X; // Reference node location.
             for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc ) {
                for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
                   for ( k = 1 ; k <= NUM_COMP_VECTOR2D; ++ k) {
		      l = l + 1;
                      if (l >= buffer_size_neighbour) return(5201);
                      Soln_Block_List.message_noreschange_northeastcorner_sendbuf[i_blk][l] =
                         Soln_ptr[i_blk].Grid.Node[i][j].X[k]-x_ref[k];
                   } /* endfor */
                } /* endfor */
             } /* endfor */
          } /* endif */
       } /* endif */

       ///////////////////////////////////////////////////////////////////////////
       // Load send buffer for South East corner neighbours of solution blocks. //
       ///////////////////////////////////////////////////////////////////////////
       if (Soln_Block_List.Block[i_blk].used && 
           (Soln_Block_List.Block[i_blk].nSE == 1) &&
           (Soln_Block_List.Block[i_blk].info.level == 
            Soln_Block_List.Block[i_blk].infoSE[0].level)) {
          buffer_size_neighbour = sqr(Soln_Block_List.Block[i_blk].infoSE[0].dimen.ghost)*
                                  Number_of_Solution_Variables;
          l = -1;
          // Load ghost cell solution variable information as required.
          if (!Send_Mesh_Geometry_Only) {
             if (Soln_Block_List.Block[i_blk].infoSE[0].dimen.i > 0 &&
                 Soln_Block_List.Block[i_blk].infoSE[0].dimen.j > 0) {
	        i_min = Soln_ptr[i_blk].ICu-Soln_Block_List.Block[i_blk].infoSE[0].dimen.ghost+1;
                i_max = Soln_ptr[i_blk].ICu;
	        i_inc = 1;
	        j_min = Soln_ptr[i_blk].JCl;
	        j_max = Soln_ptr[i_blk].JCl+Soln_Block_List.Block[i_blk].infoSE[0].dimen.ghost-1;
	        j_inc = 1;
             } else if (Soln_Block_List.Block[i_blk].infoSE[0].dimen.j > 0) {
	        i_min = Soln_ptr[i_blk].ICu;
	        i_max = Soln_ptr[i_blk].ICu-Soln_Block_List.Block[i_blk].infoSE[0].dimen.ghost+1;
	        i_inc = -1;
	        j_min = Soln_ptr[i_blk].JCl;
	        j_max = Soln_ptr[i_blk].JCl+Soln_Block_List.Block[i_blk].infoSE[0].dimen.ghost-1;
	        j_inc = 1;
             } else if (Soln_Block_List.Block[i_blk].infoSE[0].dimen.i > 0) {
	        i_min = Soln_ptr[i_blk].ICu-Soln_Block_List.Block[i_blk].infoSE[0].dimen.ghost+1;
	        i_max = Soln_ptr[i_blk].ICu;
	        i_inc = 1;
	        j_min = Soln_ptr[i_blk].JCl+Soln_Block_List.Block[i_blk].infoSE[0].dimen.ghost-1;
	        j_max = Soln_ptr[i_blk].JCl;
	        j_inc = -1;
             } else {
	        i_min = Soln_ptr[i_blk].ICu;
	        i_max = Soln_ptr[i_blk].ICu-Soln_Block_List.Block[i_blk].infoSE[0].dimen.ghost+1;
	        i_inc = -1;
	        j_min = Soln_ptr[i_blk].JCl+Soln_Block_List.Block[i_blk].infoSE[0].dimen.ghost-1;
	        j_max = Soln_ptr[i_blk].JCl;
	        j_inc = -1;
             } /* endif */
             i = Soln_ptr[i_blk].LoadSendBuffer(Soln_Block_List.message_noreschange_southeastcorner_sendbuf[i_blk],
                                                l,buffer_size_neighbour,
                                                i_min,i_max,i_inc,
                                                j_min,j_max,j_inc);
	     if (i != 0) return(4100);
          } /* endif */
          // Load ghost cell mesh information as required.
          if (Send_Mesh_Geometry_Only ||
              Number_of_Solution_Variables > Soln_ptr[i_blk].NumVar()) {
             if (Soln_Block_List.Block[i_blk].infoSE[0].dimen.i > 0 &&
                 Soln_Block_List.Block[i_blk].infoSE[0].dimen.j > 0) {
	        i_min = Soln_ptr[i_blk].Grid.INu-Soln_Block_List.Block[i_blk].infoSE[0].dimen.ghost;
	        i_max = Soln_ptr[i_blk].Grid.INu-1;
                i_inc = 1;
	        j_min = Soln_ptr[i_blk].Grid.JNl+1;
	        j_max = Soln_ptr[i_blk].Grid.JNl+Soln_Block_List.Block[i_blk].infoSE[0].dimen.ghost;
	        j_inc = 1;
                i_ref = i_max+1;
                j_ref = j_min-1;
             } else if (Soln_Block_List.Block[i_blk].infoSE[0].dimen.j > 0) {
	        i_min = Soln_ptr[i_blk].Grid.INu-1;
	        i_max = Soln_ptr[i_blk].Grid.INu-Soln_Block_List.Block[i_blk].infoSE[0].dimen.ghost;
                i_inc = -1;
	        j_min = Soln_ptr[i_blk].Grid.JNl+1;
	        j_max = Soln_ptr[i_blk].Grid.JNl+Soln_Block_List.Block[i_blk].infoSE[0].dimen.ghost;
	        j_inc = 1;
                i_ref = i_min+1;
                j_ref = j_min-1;
             } else if (Soln_Block_List.Block[i_blk].infoSE[0].dimen.i > 0) {
	        i_min = Soln_ptr[i_blk].Grid.INu-Soln_Block_List.Block[i_blk].infoSE[0].dimen.ghost;
	        i_max = Soln_ptr[i_blk].Grid.INu-1;
                i_inc = 1;
	        j_min = Soln_ptr[i_blk].Grid.JNl+Soln_Block_List.Block[i_blk].infoSE[0].dimen.ghost;
	        j_max = Soln_ptr[i_blk].Grid.JNl+1;
	        j_inc = -1;
                i_ref = i_max+1;
                j_ref = j_max-1;
             } else {
	        i_min = Soln_ptr[i_blk].Grid.INu-1;
	        i_max = Soln_ptr[i_blk].Grid.INu-Soln_Block_List.Block[i_blk].infoSE[0].dimen.ghost;
                i_inc = -1;
	        j_min = Soln_ptr[i_blk].Grid.JNl+Soln_Block_List.Block[i_blk].infoSE[0].dimen.ghost;
	        j_max = Soln_ptr[i_blk].Grid.JNl+1;
	        j_inc = -1;
                i_ref = i_min+1;
                j_ref = j_max-1;
             } /* endif */
             x_ref = Soln_ptr[i_blk].Grid.Node[i_ref][j_ref].X; // Reference node location.
             for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc ) {
                for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
                   for ( k = 1 ; k <= NUM_COMP_VECTOR2D; ++ k) {
		      l = l + 1;
                      if (l >= buffer_size_neighbour) return(4101);
                      Soln_Block_List.message_noreschange_southeastcorner_sendbuf[i_blk][l] =
                         Soln_ptr[i_blk].Grid.Node[i][j].X[k]-x_ref[k];
                   } /* endfor */
                } /* endfor */
             } /* endfor */
          } /* endif */
       } /* endif */

       ///////////////////////////////////////////////////////////////////////////
       // Load send buffer for South West corner neighbours of solution blocks. //
       ///////////////////////////////////////////////////////////////////////////
       if (Soln_Block_List.Block[i_blk].used && 
           (Soln_Block_List.Block[i_blk].nSW == 1) &&
           (Soln_Block_List.Block[i_blk].info.level == 
            Soln_Block_List.Block[i_blk].infoSW[0].level)) {
          buffer_size_neighbour = sqr(Soln_Block_List.Block[i_blk].infoSW[0].dimen.ghost)*
                                  Number_of_Solution_Variables;
          l = -1;
          // Load ghost cell solution variable information as required.
          if (!Send_Mesh_Geometry_Only) {
             if (Soln_Block_List.Block[i_blk].infoSW[0].dimen.i > 0 &&
                 Soln_Block_List.Block[i_blk].infoSW[0].dimen.j > 0) {
	        i_min = Soln_ptr[i_blk].ICl;
                i_max = Soln_ptr[i_blk].ICl+Soln_Block_List.Block[i_blk].infoSW[0].dimen.ghost-1;
	        i_inc = 1;
	        j_min = Soln_ptr[i_blk].JCl;
	        j_max = Soln_ptr[i_blk].JCl+Soln_Block_List.Block[i_blk].infoSW[0].dimen.ghost-1;
	        j_inc = 1;
             } else if (Soln_Block_List.Block[i_blk].infoSW[0].dimen.j > 0) {
	        i_min = Soln_ptr[i_blk].ICl+Soln_Block_List.Block[i_blk].infoSW[0].dimen.ghost-1;
	        i_max = Soln_ptr[i_blk].ICl;
	        i_inc = -1;
	        j_min = Soln_ptr[i_blk].JCl;
	        j_max = Soln_ptr[i_blk].JCl+Soln_Block_List.Block[i_blk].infoSW[0].dimen.ghost-1;
	        j_inc = 1;
             } else if (Soln_Block_List.Block[i_blk].infoSW[0].dimen.i > 0) {
	        i_min = Soln_ptr[i_blk].ICl;
                i_max = Soln_ptr[i_blk].ICl+Soln_Block_List.Block[i_blk].infoSW[0].dimen.ghost-1;
	        i_inc = 1;
	        j_min = Soln_ptr[i_blk].JCl+Soln_Block_List.Block[i_blk].infoSW[0].dimen.ghost-1;
	        j_max = Soln_ptr[i_blk].JCl;
	        j_inc = -1;
             } else {
	        i_min = Soln_ptr[i_blk].ICl+Soln_Block_List.Block[i_blk].infoSW[0].dimen.ghost-1;
	        i_max = Soln_ptr[i_blk].ICl;
	        i_inc = -1;
	        j_min = Soln_ptr[i_blk].JCl+Soln_Block_List.Block[i_blk].infoSW[0].dimen.ghost-1;
	        j_max = Soln_ptr[i_blk].JCl;
	        j_inc = -1;
             } /* endif */
             i = Soln_ptr[i_blk].LoadSendBuffer(Soln_Block_List.message_noreschange_southwestcorner_sendbuf[i_blk],
                                                l,buffer_size_neighbour,
                                                i_min,i_max,i_inc,
                                                j_min,j_max,j_inc);
	     if (i != 0) return(5100);
          } /* endif */
          // Load ghost cell mesh information as required.
          if (Send_Mesh_Geometry_Only ||
              Number_of_Solution_Variables > Soln_ptr[i_blk].NumVar()) {
             if (Soln_Block_List.Block[i_blk].infoSW[0].dimen.i > 0 &&
                 Soln_Block_List.Block[i_blk].infoSW[0].dimen.j > 0) {
	        i_min = Soln_ptr[i_blk].Grid.INl+1;
	        i_max = Soln_ptr[i_blk].Grid.INl+Soln_Block_List.Block[i_blk].infoSW[0].dimen.ghost;
                i_inc = 1;
	        j_min = Soln_ptr[i_blk].Grid.JNl+1;
	        j_max = Soln_ptr[i_blk].Grid.JNl+Soln_Block_List.Block[i_blk].infoSW[0].dimen.ghost;
	        j_inc = 1;
                i_ref = i_min-1;
                j_ref = j_min-1;
             } else if (Soln_Block_List.Block[i_blk].infoSW[0].dimen.j > 0) {
	        i_min = Soln_ptr[i_blk].Grid.INl+Soln_Block_List.Block[i_blk].infoSW[0].dimen.ghost;
	        i_max = Soln_ptr[i_blk].Grid.INl+1;
                i_inc = -1;
	        j_min = Soln_ptr[i_blk].Grid.JNl+1;
	        j_max = Soln_ptr[i_blk].Grid.JNl+Soln_Block_List.Block[i_blk].infoSW[0].dimen.ghost;
	        j_inc = 1;
                i_ref = i_max-1;
                j_ref = j_min-1;
             } else if (Soln_Block_List.Block[i_blk].infoSW[0].dimen.i > 0) {
	        i_min = Soln_ptr[i_blk].Grid.INl+1;
	        i_max = Soln_ptr[i_blk].Grid.INl+Soln_Block_List.Block[i_blk].infoSW[0].dimen.ghost;
                i_inc = 1;
	        j_min = Soln_ptr[i_blk].Grid.JNl+Soln_Block_List.Block[i_blk].infoSW[0].dimen.ghost;
	        j_max = Soln_ptr[i_blk].Grid.JNl+1;
	        j_inc = -1;
                i_ref = i_min-1;
                j_ref = j_max-1;
             } else {
	        i_min = Soln_ptr[i_blk].Grid.INl+Soln_Block_List.Block[i_blk].infoSW[0].dimen.ghost;
	        i_max = Soln_ptr[i_blk].Grid.INl+1;
                i_inc = -1;
	        j_min = Soln_ptr[i_blk].Grid.JNl+Soln_Block_List.Block[i_blk].infoSW[0].dimen.ghost;
	        j_max = Soln_ptr[i_blk].Grid.JNl+1;
	        j_inc = -1;
                i_ref = i_max-1;
                j_ref = j_max-1;
             } /* endif */
             x_ref = Soln_ptr[i_blk].Grid.Node[i_ref][j_ref].X; // Reference node location.
             for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc ) {
                for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
                   for ( k = 1 ; k <= NUM_COMP_VECTOR2D; ++ k) {
		      l = l + 1;
                      if (l >= buffer_size_neighbour) return(5101);
                      Soln_Block_List.message_noreschange_southwestcorner_sendbuf[i_blk][l] =
                         Soln_ptr[i_blk].Grid.Node[i][j].X[k]-x_ref[k];
                   } /* endfor */
                } /* endfor */
             } /* endfor */
          } /* endif */
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
template <class Quad_Soln_Block>
int Load_Send_Message_Buffers_ResChange_FineToCoarse(Quad_Soln_Block *Soln_ptr,
                                                     AdaptiveBlock2D_List &Soln_Block_List,
                                                     const int Number_of_Solution_Variables,
                                                     const int Send_Mesh_Geometry_Only,
                                                     const int Send_Conservative_Solution_Fluxes) {

    int i_blk, buffer_size_neighbour, 
        i_min, i_max, i_inc, i, 
        j_min, j_max, j_inc, j, 
        k, l;
    int i_ref, j_ref;
    Vector2D x_ref;

    /* Check to see if the number of solution variables specified
       in the calling argument is correct. */

    if (Number_of_Solution_Variables != NUM_COMP_VECTOR2D &&
        Number_of_Solution_Variables != Soln_ptr[0].NumVar() &&
        Number_of_Solution_Variables != Soln_ptr[0].NumVar() + NUM_COMP_VECTOR2D) {
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
          if (!Send_Mesh_Geometry_Only) {
	     if (!Send_Conservative_Solution_Fluxes) {
                buffer_size_neighbour = Soln_Block_List.Block[i_blk].infoN[0].dimen.ghost*
                                        (abs(Soln_Block_List.Block[i_blk].infoN[0].dimen.i)/2)*
                                        Soln_ptr[i_blk].NumVar();
             } else {
                buffer_size_neighbour = (abs(Soln_Block_List.Block[i_blk].infoN[0].dimen.i)/2)*
                                        Soln_ptr[i_blk].NumVar();
             } /* endif */
             if (Number_of_Solution_Variables > Soln_ptr[i_blk].NumVar()) 
                buffer_size_neighbour = buffer_size_neighbour +
                                        Soln_Block_List.Block[i_blk].infoN[0].dimen.ghost*
                                        (abs(Soln_Block_List.Block[i_blk].infoN[0].dimen.i)/2+1)*
                                        NUM_COMP_VECTOR2D+
                                        2*Soln_Block_List.Block[i_blk].infoN[0].dimen.ghost;
          } else {
             buffer_size_neighbour = Soln_Block_List.Block[i_blk].infoN[0].dimen.ghost*
                                     (abs(Soln_Block_List.Block[i_blk].infoN[0].dimen.i)/2+1)*
                                     NUM_COMP_VECTOR2D+
                                     2*Soln_Block_List.Block[i_blk].infoN[0].dimen.ghost;
          } /* endif */
          l = -1;
          // Load ghost cell solution variable information as required.
          if (!Send_Mesh_Geometry_Only) {
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
	     if (!Send_Conservative_Solution_Fluxes) {
                i = Soln_ptr[i_blk].LoadSendBuffer_F2C(Soln_Block_List.message_reschange_northface_sendbuf[i_blk][0],
                                                       l,buffer_size_neighbour,
                                                       i_min,i_max,i_inc,
                                                       j_min,j_max,j_inc);
             } else {
	        j_max = max(j_min, j_max);
                i = Soln_ptr[i_blk].LoadSendBuffer_Flux_F2C(Soln_Block_List.message_reschange_northface_sendbuf[i_blk][0],
                                                            l,buffer_size_neighbour,
                                                            i_min,i_max,i_inc,
                                                            j_max,j_max,j_inc);
             } /* endif */
	     if (i != 0) return(2200);
          } /* endif */
          // Load ghost cell mesh information as required.
          if (Send_Mesh_Geometry_Only ||
              Number_of_Solution_Variables > Soln_ptr[i_blk].NumVar()) {
             if (Soln_Block_List.Block[i_blk].infoN[0].dimen.i > 0 &&
                 Soln_Block_List.Block[i_blk].infoN[0].dimen.j > 0) {
	        i_min = Soln_ptr[i_blk].Grid.INl;
	        i_max = Soln_ptr[i_blk].Grid.INu;
	        i_inc = 2;
	        j_min = Soln_ptr[i_blk].Grid.JNu-2*Soln_Block_List.Block[i_blk].infoN[0].dimen.ghost;
	        j_max = Soln_ptr[i_blk].Grid.JNu-2;
	        j_inc = 2;
                i_ref = i_min;
                j_ref = j_max+2;
             } else if (Soln_Block_List.Block[i_blk].infoN[0].dimen.j > 0) {
	        i_min = Soln_ptr[i_blk].Grid.INu;
	        i_max = Soln_ptr[i_blk].Grid.INl;
	        i_inc = -2;
	        j_min = Soln_ptr[i_blk].Grid.JNu-2*Soln_Block_List.Block[i_blk].infoN[0].dimen.ghost;
	        j_max = Soln_ptr[i_blk].Grid.JNu-2;
	        j_inc = 2;
                i_ref = i_min;
                j_ref = j_max+2;
             } else if (Soln_Block_List.Block[i_blk].infoN[0].dimen.i > 0) {
	        i_min = Soln_ptr[i_blk].Grid.INl;
	        i_max = Soln_ptr[i_blk].Grid.INu;
	        i_inc = 2;
	        j_min = Soln_ptr[i_blk].Grid.JNu-2;
	        j_max = Soln_ptr[i_blk].Grid.JNu-2*Soln_Block_List.Block[i_blk].infoN[0].dimen.ghost;
	        j_inc = -2;
                i_ref = i_min;
                j_ref = j_min+2;
             } else {
	        i_min = Soln_ptr[i_blk].Grid.INu;
	        i_max = Soln_ptr[i_blk].Grid.INl;
	        i_inc = -2;
	        j_min = Soln_ptr[i_blk].Grid.JNu-2;
	        j_max = Soln_ptr[i_blk].Grid.JNu-2*Soln_Block_List.Block[i_blk].infoN[0].dimen.ghost;
	        j_inc = -2;
                i_ref = i_min;
                j_ref = j_min+2;
             } /* endif */
             x_ref = Soln_ptr[i_blk].Grid.Node[i_ref][j_ref].X; // Reference node location.
             for ( j  = j_min ; ((j_inc+2)/4) ? (j <= j_max):(j >= j_max) ; j += j_inc ) {
                for ( i = i_min ;  ((i_inc+2)/4) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
                   for ( k = 1 ; k <= NUM_COMP_VECTOR2D; ++ k) {
		      l = l + 1;
                      if (l >= buffer_size_neighbour) return(2201);
                      Soln_Block_List.message_reschange_northface_sendbuf[i_blk][0][l] =
                         Soln_ptr[i_blk].Grid.Node[i][j].X[k]-x_ref[k];
                   } /* endfor */
               } /* endfor */
             } /* endfor */
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
             if (i_inc > 0) {
                for ( j  = j_min ; ((j_inc+2)/4) ? (j <= j_max):(j >= j_max) ; j += j_inc ) {
                   l = l + 1;
                   if (l >= buffer_size_neighbour) return(2202);
                   Soln_Block_List.message_reschange_northface_sendbuf[i_blk][0][l] = 
                      double(Soln_ptr[i_blk].Grid.BCtypeW[j]);
                } /* endfor */
                for ( j  = j_min ; ((j_inc+2)/4) ? (j <= j_max):(j >= j_max) ; j += j_inc ) {
                   l = l + 1;
                   if (l >= buffer_size_neighbour) return(2203);
                   Soln_Block_List.message_reschange_northface_sendbuf[i_blk][0][l] = 
                      double(Soln_ptr[i_blk].Grid.BCtypeE[j]);
                } /* endfor */
	     } else {
                for ( j  = j_min ; ((j_inc+2)/4) ? (j <= j_max):(j >= j_max) ; j += j_inc ) {
                   l = l + 1;
                   if (l >= buffer_size_neighbour) return(2204);
                   Soln_Block_List.message_reschange_northface_sendbuf[i_blk][0][l] = 
                      double(Soln_ptr[i_blk].Grid.BCtypeE[j]);
                } /* endfor */
                for ( j  = j_min ; ((j_inc+2)/4) ? (j <= j_max):(j >= j_max) ; j += j_inc ) {
                   l = l + 1;
                   if (l >= buffer_size_neighbour) return(2205);
                   Soln_Block_List.message_reschange_northface_sendbuf[i_blk][0][l] = 
                      double(Soln_ptr[i_blk].Grid.BCtypeW[j]);
                } /* endfor */
             } /* endif */
          } /* endif */
       } /* endif */

       //////////////////////////////////////////////////////////////////
       // Load send buffers of solution blocks with coarser South face //
       // neighbour at lower mesh resolution.                          //
       //////////////////////////////////////////////////////////////////
       if (Soln_Block_List.Block[i_blk].used && 
           (Soln_Block_List.Block[i_blk].nS == 1) &&
           (Soln_Block_List.Block[i_blk].info.level > 
            Soln_Block_List.Block[i_blk].infoS[0].level)) {
          if (!Send_Mesh_Geometry_Only) {
	     if (!Send_Conservative_Solution_Fluxes) {
                buffer_size_neighbour = Soln_Block_List.Block[i_blk].infoS[0].dimen.ghost*
                                        (abs(Soln_Block_List.Block[i_blk].infoS[0].dimen.i)/2)*
                                        Soln_ptr[i_blk].NumVar();
             } else {
                buffer_size_neighbour = (abs(Soln_Block_List.Block[i_blk].infoS[0].dimen.i)/2)*
                                        Soln_ptr[i_blk].NumVar();
             } /* endif */
             if (Number_of_Solution_Variables > Soln_ptr[i_blk].NumVar()) 
                buffer_size_neighbour = buffer_size_neighbour +
                                        Soln_Block_List.Block[i_blk].infoS[0].dimen.ghost*
                                        (abs(Soln_Block_List.Block[i_blk].infoS[0].dimen.i)/2+1)*
                                        NUM_COMP_VECTOR2D+
                                        2*Soln_Block_List.Block[i_blk].infoS[0].dimen.ghost;
          } else {
             buffer_size_neighbour = Soln_Block_List.Block[i_blk].infoS[0].dimen.ghost*
                                     (abs(Soln_Block_List.Block[i_blk].infoS[0].dimen.i)/2+1)*
                                     NUM_COMP_VECTOR2D+
                                     2*Soln_Block_List.Block[i_blk].infoS[0].dimen.ghost;
          } /* endif */
          l = -1;
          // Load ghost cell solution variable information as required.
          if (!Send_Mesh_Geometry_Only) {
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
	     if (!Send_Conservative_Solution_Fluxes) {
                i = Soln_ptr[i_blk].LoadSendBuffer_F2C(Soln_Block_List.message_reschange_southface_sendbuf[i_blk][0],
                                                       l,buffer_size_neighbour,
                                                       i_min,i_max,i_inc,
                                                       j_min,j_max,j_inc);
             } else {
	        j_min = min(j_min, j_max);
                i = Soln_ptr[i_blk].LoadSendBuffer_Flux_F2C(Soln_Block_List.message_reschange_southface_sendbuf[i_blk][0],
                                                            l,buffer_size_neighbour,
                                                            i_min,i_max,i_inc,
                                                            j_min,j_min,j_inc);
             } /* endif */
	     if (i != 0) return(2100);
          } /* endif */
          // Load ghost cell mesh information as required.
          if (Send_Mesh_Geometry_Only ||
              Number_of_Solution_Variables > Soln_ptr[i_blk].NumVar()) {
             if (Soln_Block_List.Block[i_blk].infoS[0].dimen.i > 0 &&
                 Soln_Block_List.Block[i_blk].infoS[0].dimen.j > 0) {
	        i_min = Soln_ptr[i_blk].Grid.INl;
	        i_max = Soln_ptr[i_blk].Grid.INu;
	        i_inc = 2;
	        j_min = Soln_ptr[i_blk].Grid.JNl+2;
	        j_max = Soln_ptr[i_blk].Grid.JNl+2*Soln_Block_List.Block[i_blk].infoS[0].dimen.ghost;
	        j_inc = 2;
                i_ref = i_min;
                j_ref = j_min-2;
             } else if (Soln_Block_List.Block[i_blk].infoS[0].dimen.j > 0) {
	        i_min = Soln_ptr[i_blk].Grid.INu;
	        i_max = Soln_ptr[i_blk].Grid.INl;
	        i_inc = -2;
	        j_min = Soln_ptr[i_blk].Grid.JNl+2;
	        j_max = Soln_ptr[i_blk].Grid.JNl+2*Soln_Block_List.Block[i_blk].infoS[0].dimen.ghost;
	        j_inc = 2;
                i_ref = i_min;
                j_ref = j_min-2;
             } else if (Soln_Block_List.Block[i_blk].infoS[0].dimen.i > 0) {
	        i_min = Soln_ptr[i_blk].Grid.INl;
	        i_max = Soln_ptr[i_blk].Grid.INu;
	        i_inc = 2;
	        j_min = Soln_ptr[i_blk].Grid.JNl+2*Soln_Block_List.Block[i_blk].infoS[0].dimen.ghost;
	        j_max = Soln_ptr[i_blk].Grid.JNl+2;
	        j_inc = -2;
                i_ref = i_min;
                j_ref = j_max-2;
             } else {
	        i_min = Soln_ptr[i_blk].Grid.INu;
	        i_max = Soln_ptr[i_blk].Grid.INl;
	        i_inc = -2;
	        j_min = Soln_ptr[i_blk].Grid.JNl+2*Soln_Block_List.Block[i_blk].infoS[0].dimen.ghost;
	        j_max = Soln_ptr[i_blk].Grid.JNl+2;
	        j_inc = -2;
                i_ref = i_min;
                j_ref = j_max-2;
             } /* endif */
             x_ref = Soln_ptr[i_blk].Grid.Node[i_ref][j_ref].X; // Reference node location.
             for ( j  = j_min ; ((j_inc+2)/4) ? (j <= j_max):(j >= j_max) ; j += j_inc ) {
                for ( i = i_min ;  ((i_inc+2)/4) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
                   for ( k = 1 ; k <= NUM_COMP_VECTOR2D; ++ k) {
		      l = l + 1;
                      if (l >= buffer_size_neighbour) return(2101);
                      Soln_Block_List.message_reschange_southface_sendbuf[i_blk][0][l] =
                         Soln_ptr[i_blk].Grid.Node[i][j].X[k]-x_ref[k];
                   } /* endfor */
               } /* endfor */
             } /* endfor */
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
             if (i_inc > 0) {
                for ( j  = j_min ; ((j_inc+2)/4) ? (j <= j_max):(j >= j_max) ; j += j_inc ) {
                   l = l + 1;
                   if (l >= buffer_size_neighbour) return(2102);
                   Soln_Block_List.message_reschange_southface_sendbuf[i_blk][0][l] = 
                      double(Soln_ptr[i_blk].Grid.BCtypeW[j]);
                } /* endfor */
                for ( j  = j_min ; ((j_inc+2)/4) ? (j <= j_max):(j >= j_max) ; j += j_inc ) {
                   l = l + 1;
                   if (l >= buffer_size_neighbour) return(2103);
                   Soln_Block_List.message_reschange_southface_sendbuf[i_blk][0][l] = 
                      double(Soln_ptr[i_blk].Grid.BCtypeE[j]);
                } /* endfor */
	     } else {
                for ( j  = j_min ; ((j_inc+2)/4) ? (j <= j_max):(j >= j_max) ; j += j_inc ) {
                   l = l + 1;
                   if (l >= buffer_size_neighbour) return(2104);
                   Soln_Block_List.message_reschange_southface_sendbuf[i_blk][0][l] = 
                      double(Soln_ptr[i_blk].Grid.BCtypeE[j]);
                } /* endfor */
                for ( j  = j_min ; ((j_inc+2)/4) ? (j <= j_max):(j >= j_max) ; j += j_inc ) {
                   l = l + 1;
                   if (l >= buffer_size_neighbour) return(2105);
                   Soln_Block_List.message_reschange_southface_sendbuf[i_blk][0][l] = 
                      double(Soln_ptr[i_blk].Grid.BCtypeW[j]);
                } /* endfor */
             } /* endif */
          } /* endif */
       } /* endif */

       /////////////////////////////////////////////////////////////////
       // Load send buffers of solution blocks with coarser East face //
       // neighbour at lower mesh resolution.                         //
       /////////////////////////////////////////////////////////////////
       if (Soln_Block_List.Block[i_blk].used && 
           (Soln_Block_List.Block[i_blk].nE == 1) &&
           (Soln_Block_List.Block[i_blk].info.level > 
            Soln_Block_List.Block[i_blk].infoE[0].level)) {
          if (!Send_Mesh_Geometry_Only) {
	     if (!Send_Conservative_Solution_Fluxes) {
                buffer_size_neighbour = Soln_Block_List.Block[i_blk].infoE[0].dimen.ghost*
                                        (abs(Soln_Block_List.Block[i_blk].infoE[0].dimen.j)/2)*
                                        Soln_ptr[i_blk].NumVar();
             } else {
                buffer_size_neighbour = (abs(Soln_Block_List.Block[i_blk].infoE[0].dimen.j)/2)*
                                        Soln_ptr[i_blk].NumVar();
             } /* endif */
             if (Number_of_Solution_Variables > Soln_ptr[i_blk].NumVar()) 
                buffer_size_neighbour = buffer_size_neighbour +
                                        Soln_Block_List.Block[i_blk].infoE[0].dimen.ghost*
                                        (abs(Soln_Block_List.Block[i_blk].infoE[0].dimen.j)/2+1)*
                                        NUM_COMP_VECTOR2D+
                                        2*Soln_Block_List.Block[i_blk].infoE[0].dimen.ghost;
          } else {
             buffer_size_neighbour = Soln_Block_List.Block[i_blk].infoE[0].dimen.ghost*
                                     (abs(Soln_Block_List.Block[i_blk].infoE[0].dimen.j)/2+1)*
                                     NUM_COMP_VECTOR2D+
                                     2*Soln_Block_List.Block[i_blk].infoE[0].dimen.ghost;
          } /* endif */
          l = -1;
          // Load ghost cell solution variable information as required.
          if (!Send_Mesh_Geometry_Only) {
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
	     if (!Send_Conservative_Solution_Fluxes) {
                i = Soln_ptr[i_blk].LoadSendBuffer_F2C(Soln_Block_List.message_reschange_eastface_sendbuf[i_blk][0],
                                                       l,buffer_size_neighbour,
                                                       i_min,i_max,i_inc,
                                                       j_min,j_max,j_inc);
             } else {
	        i_max = max(i_min, i_max);
                i = Soln_ptr[i_blk].LoadSendBuffer_Flux_F2C(Soln_Block_List.message_reschange_eastface_sendbuf[i_blk][0],
                                                            l,buffer_size_neighbour,
                                                            i_max,i_max,i_inc,
                                                            j_min,j_max,j_inc);
             } /* endif */
	     if (i != 0) return(1200);
          } /* endif */
          // Load ghost cell mesh information as required.
          if (Send_Mesh_Geometry_Only ||
              Number_of_Solution_Variables > Soln_ptr[i_blk].NumVar()) {
             if (Soln_Block_List.Block[i_blk].infoE[0].dimen.i > 0 &&
                 Soln_Block_List.Block[i_blk].infoE[0].dimen.j > 0) {
	        i_min = Soln_ptr[i_blk].Grid.INu-2*Soln_Block_List.Block[i_blk].infoE[0].dimen.ghost;
	        i_max = Soln_ptr[i_blk].Grid.INu-2;
	        i_inc = 2;
	        j_min = Soln_ptr[i_blk].Grid.JNl;
	        j_max = Soln_ptr[i_blk].Grid.JNu;
	        j_inc = 2;
                i_ref = i_max+2;
                j_ref = j_min;
             } else if (Soln_Block_List.Block[i_blk].infoE[0].dimen.j > 0) {
	        i_min = Soln_ptr[i_blk].Grid.INu-2;
	        i_max = Soln_ptr[i_blk].Grid.INu-2*Soln_Block_List.Block[i_blk].infoE[0].dimen.ghost;
	        i_inc = -2;
	        j_min = Soln_ptr[i_blk].Grid.JNl;
	        j_max = Soln_ptr[i_blk].Grid.JNu;
	        j_inc = 2;
                i_ref = i_min+2;
                j_ref = j_min;
             } else if (Soln_Block_List.Block[i_blk].infoE[0].dimen.i > 0) {
	        i_min = Soln_ptr[i_blk].Grid.INu-2*Soln_Block_List.Block[i_blk].infoE[0].dimen.ghost;
	        i_max = Soln_ptr[i_blk].Grid.INu-2;
	        i_inc = 2;
	        j_min = Soln_ptr[i_blk].Grid.JNu;
	        j_max = Soln_ptr[i_blk].Grid.JNl;
	        j_inc = -2;
                i_ref = i_max+2;
                j_ref = j_min;
             } else {
	        i_min = Soln_ptr[i_blk].Grid.INu-2;
	        i_max = Soln_ptr[i_blk].Grid.INu-2*Soln_Block_List.Block[i_blk].infoE[0].dimen.ghost;
	        i_inc = -2;
	        j_min = Soln_ptr[i_blk].Grid.JNu;
	        j_max = Soln_ptr[i_blk].Grid.JNl;
	        j_inc = -2;
                i_ref = i_min+2;
                j_ref = j_min;
             } /* endif */
             x_ref = Soln_ptr[i_blk].Grid.Node[i_ref][j_ref].X; // Reference node location.
             for ( j  = j_min ; ((j_inc+2)/4) ? (j <= j_max):(j >= j_max) ; j += j_inc ) {
                for ( i = i_min ;  ((i_inc+2)/4) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
                   for ( k = 1 ; k <= NUM_COMP_VECTOR2D; ++ k) {
		      l = l + 1;
                      if (l >= buffer_size_neighbour) return(1201);
                      Soln_Block_List.message_reschange_eastface_sendbuf[i_blk][0][l] =
                         Soln_ptr[i_blk].Grid.Node[i][j].X[k]-x_ref[k];
                   } /* endfor */
               } /* endfor */
             } /* endfor */
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
             if (j_inc > 0) {
                for ( i = i_min ;  ((i_inc+2)/4) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
                   l = l + 1;
                   if (l >= buffer_size_neighbour) return(1202);
                   Soln_Block_List.message_reschange_eastface_sendbuf[i_blk][0][l] = 
                      double(Soln_ptr[i_blk].Grid.BCtypeS[i]);
                } /* endfor */
                for ( i = i_min ;  ((i_inc+2)/4) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
                   l = l + 1;
                   if (l >= buffer_size_neighbour) return(1203);
                   Soln_Block_List.message_reschange_eastface_sendbuf[i_blk][0][l] = 
                      double(Soln_ptr[i_blk].Grid.BCtypeN[i]);
                } /* endfor */
	     } else {
                for ( i = i_min ;  ((i_inc+2)/4) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
                   l = l + 1;
                   if (l >= buffer_size_neighbour) return(1204);
                   Soln_Block_List.message_reschange_eastface_sendbuf[i_blk][0][l] = 
                      double(Soln_ptr[i_blk].Grid.BCtypeN[i]);
                } /* endfor */
                for ( i = i_min ;  ((i_inc+2)/4) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
                   l = l + 1;
                   if (l >= buffer_size_neighbour) return(1205);
                   Soln_Block_List.message_reschange_eastface_sendbuf[i_blk][0][l] = 
                      double(Soln_ptr[i_blk].Grid.BCtypeS[i]);
                } /* endfor */
             } /* endif */
          } /* endif */
       } /* endif */

       /////////////////////////////////////////////////////////////////
       // Load send buffers of solution blocks with coarser West face //
       // neighbour at lower mesh resolution.                         //
       /////////////////////////////////////////////////////////////////
       if (Soln_Block_List.Block[i_blk].used && 
           (Soln_Block_List.Block[i_blk].nW == 1) &&
           (Soln_Block_List.Block[i_blk].info.level > 
            Soln_Block_List.Block[i_blk].infoW[0].level)) {
          if (!Send_Mesh_Geometry_Only) {
	     if (!Send_Conservative_Solution_Fluxes) {
                buffer_size_neighbour = Soln_Block_List.Block[i_blk].infoW[0].dimen.ghost*
                                        (abs(Soln_Block_List.Block[i_blk].infoW[0].dimen.j)/2)*
                                        Soln_ptr[i_blk].NumVar();
             } else {
                buffer_size_neighbour = (abs(Soln_Block_List.Block[i_blk].infoW[0].dimen.j)/2)*
                                        Soln_ptr[i_blk].NumVar();
             } /* endif */
             if (Number_of_Solution_Variables > Soln_ptr[i_blk].NumVar()) 
                buffer_size_neighbour = buffer_size_neighbour +
                                        Soln_Block_List.Block[i_blk].infoW[0].dimen.ghost*
                                        (abs(Soln_Block_List.Block[i_blk].infoW[0].dimen.j)/2+1)*
                                        NUM_COMP_VECTOR2D+
                                        2*Soln_Block_List.Block[i_blk].infoW[0].dimen.ghost;
          } else {
             buffer_size_neighbour = Soln_Block_List.Block[i_blk].infoW[0].dimen.ghost*
                                     (abs(Soln_Block_List.Block[i_blk].infoW[0].dimen.j)/2+1)*
                                     NUM_COMP_VECTOR2D+
                                     2*Soln_Block_List.Block[i_blk].infoW[0].dimen.ghost;
          } /* endif */
          l = -1;
          // Load ghost cell solution variable information as required.
          if (!Send_Mesh_Geometry_Only) {
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
	     if (!Send_Conservative_Solution_Fluxes) {
                i = Soln_ptr[i_blk].LoadSendBuffer_F2C(Soln_Block_List.message_reschange_westface_sendbuf[i_blk][0],
                                                       l,buffer_size_neighbour,
                                                       i_min,i_max,i_inc,
                                                       j_min,j_max,j_inc);
             } else {
	        i_min = min(i_min, i_max);
                i = Soln_ptr[i_blk].LoadSendBuffer_Flux_F2C(Soln_Block_List.message_reschange_westface_sendbuf[i_blk][0],
                                                            l,buffer_size_neighbour,
                                                            i_min,i_min,i_inc,
                                                            j_min,j_max,j_inc);
             } /* endif */
	     if (i != 0) return(1100);
          } /* endif */
          // Load ghost cell mesh information as required.
          if (Send_Mesh_Geometry_Only ||
              Number_of_Solution_Variables > Soln_ptr[i_blk].NumVar()) {
             if (Soln_Block_List.Block[i_blk].infoW[0].dimen.i > 0 &&
                 Soln_Block_List.Block[i_blk].infoW[0].dimen.j > 0) {
	        i_min = Soln_ptr[i_blk].Grid.INl+2;
	        i_max = Soln_ptr[i_blk].Grid.INl+2*Soln_Block_List.Block[i_blk].infoW[0].dimen.ghost;
	        i_inc = 2;
	        j_min = Soln_ptr[i_blk].Grid.JNl;
	        j_max = Soln_ptr[i_blk].Grid.JNu;
	        j_inc = 2;
                i_ref = i_min-2;
                j_ref = j_min;
             } else if (Soln_Block_List.Block[i_blk].infoW[0].dimen.j > 0) {
	        i_min = Soln_ptr[i_blk].Grid.INl+2*Soln_Block_List.Block[i_blk].infoW[0].dimen.ghost;
	        i_max = Soln_ptr[i_blk].Grid.INl+2;
	        i_inc = -2;
	        j_min = Soln_ptr[i_blk].Grid.JNl;
	        j_max = Soln_ptr[i_blk].Grid.JNu;
	        j_inc = 2;
                i_ref = i_max-2;
                j_ref = j_min;
             } else if (Soln_Block_List.Block[i_blk].infoW[0].dimen.i > 0) {
	        i_min = Soln_ptr[i_blk].Grid.INl+2;
	        i_max = Soln_ptr[i_blk].Grid.INl+2*Soln_Block_List.Block[i_blk].infoW[0].dimen.ghost;
	        i_inc = 2;
	        j_min = Soln_ptr[i_blk].Grid.JNu;
	        j_max = Soln_ptr[i_blk].Grid.JNl;
	        j_inc = -2;
                i_ref = i_min-2;
                j_ref = j_min;
             } else {
	        i_min = Soln_ptr[i_blk].Grid.INl+2*Soln_Block_List.Block[i_blk].infoW[0].dimen.ghost;
	        i_max = Soln_ptr[i_blk].Grid.INl+2;
	        i_inc = -2;
	        j_min = Soln_ptr[i_blk].Grid.JNu;
	        j_max = Soln_ptr[i_blk].Grid.JNl;
	        j_inc = -2;
                i_ref = i_max-2;
                j_ref = j_min;
             } /* endif */
             x_ref = Soln_ptr[i_blk].Grid.Node[i_ref][j_ref].X; // Reference node location.
             for ( j  = j_min ; ((j_inc+2)/4) ? (j <= j_max):(j >= j_max) ; j += j_inc ) {
                for ( i = i_min ;  ((i_inc+2)/4) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
                   for ( k = 1 ; k <= NUM_COMP_VECTOR2D; ++ k) {
		      l = l + 1;
                      if (l >= buffer_size_neighbour) return(1101);
                      Soln_Block_List.message_reschange_westface_sendbuf[i_blk][0][l] =
                         Soln_ptr[i_blk].Grid.Node[i][j].X[k]-x_ref[k];
                   } /* endfor */
               } /* endfor */
             } /* endfor */
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
             if (j_inc > 0) {
                for ( i = i_min ;  ((i_inc+2)/4) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
                   l = l + 1;
                   if (l >= buffer_size_neighbour) return(1102);
                   Soln_Block_List.message_reschange_westface_sendbuf[i_blk][0][l] = 
                      double(Soln_ptr[i_blk].Grid.BCtypeS[i]);
                } /* endfor */
                for ( i = i_min ;  ((i_inc+2)/4) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
                   l = l + 1;
                   if (l >= buffer_size_neighbour) return(1103);
                   Soln_Block_List.message_reschange_westface_sendbuf[i_blk][0][l] = 
                      double(Soln_ptr[i_blk].Grid.BCtypeN[i]);
                } /* endfor */
	     } else {
                for ( i = i_min ;  ((i_inc+2)/4) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
                   l = l + 1;
                   if (l >= buffer_size_neighbour) return(1104);
                   Soln_Block_List.message_reschange_westface_sendbuf[i_blk][0][l] = 
                      double(Soln_ptr[i_blk].Grid.BCtypeN[i]);
                } /* endfor */
                for ( i = i_min ;  ((i_inc+2)/4) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
                   l = l + 1;
                   if (l >= buffer_size_neighbour) return(1105);
                   Soln_Block_List.message_reschange_westface_sendbuf[i_blk][0][l] = 
                      double(Soln_ptr[i_blk].Grid.BCtypeS[i]);
                } /* endfor */
             } /* endif */
          } /* endif */
       } /* endif */

       /////////////////////////////////////////////////////////////////////////
       // Load send buffers of solution blocks with coarser North West corner //
       // neighbour at lower mesh resolution.                                 //
       /////////////////////////////////////////////////////////////////////////
       if (Soln_Block_List.Block[i_blk].used && 
           (Soln_Block_List.Block[i_blk].nNW == 1) &&
           (Soln_Block_List.Block[i_blk].info.level > 
            Soln_Block_List.Block[i_blk].infoNW[0].level)) {
          buffer_size_neighbour = sqr(Soln_Block_List.Block[i_blk].infoNW[0].dimen.ghost)*
                                  Number_of_Solution_Variables;
          l = -1;
          // Load ghost cell solution variable information as required.
          if (!Send_Mesh_Geometry_Only && !Send_Conservative_Solution_Fluxes) {
             if (Soln_Block_List.Block[i_blk].infoNW[0].dimen.i > 0 &&
                 Soln_Block_List.Block[i_blk].infoNW[0].dimen.j > 0) {
	        i_min = Soln_ptr[i_blk].ICl;
                i_max = Soln_ptr[i_blk].ICl+2*Soln_Block_List.Block[i_blk].infoNW[0].dimen.ghost-1;
	        i_inc = 2;
	        j_min = Soln_ptr[i_blk].JCu-2*Soln_Block_List.Block[i_blk].infoNW[0].dimen.ghost+1;
	        j_max = Soln_ptr[i_blk].JCu;
	        j_inc = 2;
             } else if (Soln_Block_List.Block[i_blk].infoNW[0].dimen.j > 0) {
	        i_min = Soln_ptr[i_blk].ICl+2*Soln_Block_List.Block[i_blk].infoNW[0].dimen.ghost-1;
	        i_max = Soln_ptr[i_blk].ICl;
	        i_inc = -2;
	        j_min = Soln_ptr[i_blk].JCu-2*Soln_Block_List.Block[i_blk].infoNW[0].dimen.ghost+1;
	        j_max = Soln_ptr[i_blk].JCu;
	        j_inc = 2;
             } else if (Soln_Block_List.Block[i_blk].infoNW[0].dimen.i > 0) {
	        i_min = Soln_ptr[i_blk].ICl;
	        i_max = Soln_ptr[i_blk].ICl+2*Soln_Block_List.Block[i_blk].infoNW[0].dimen.ghost-1;
	        i_inc = 2;
	        j_min = Soln_ptr[i_blk].JCu;
	        j_max = Soln_ptr[i_blk].JCu-2*Soln_Block_List.Block[i_blk].infoNW[0].dimen.ghost+1;
	        j_inc = -2;
             } else {
	        i_min = Soln_ptr[i_blk].ICl+2*Soln_Block_List.Block[i_blk].infoNW[0].dimen.ghost-1;
	        i_max = Soln_ptr[i_blk].ICl;
	        i_inc = -2;
	        j_min = Soln_ptr[i_blk].JCu;
	        j_max = Soln_ptr[i_blk].JCu-2*Soln_Block_List.Block[i_blk].infoNW[0].dimen.ghost+1;
	        j_inc = -2;
             } /* endif */
             i = Soln_ptr[i_blk].LoadSendBuffer_F2C(Soln_Block_List.message_reschange_northwestcorner_sendbuf[i_blk],
                                                    l,buffer_size_neighbour,
                                                    i_min,i_max,i_inc,
                                                    j_min,j_max,j_inc);
	     if (i != 0) return(4200);
          } /* endif */
          // Load ghost cell mesh information as required.
          if (Send_Mesh_Geometry_Only ||
              Number_of_Solution_Variables > Soln_ptr[i_blk].NumVar()) {
             if (Soln_Block_List.Block[i_blk].infoNW[0].dimen.i > 0 &&
                 Soln_Block_List.Block[i_blk].infoNW[0].dimen.j > 0) {
	        i_min = Soln_ptr[i_blk].Grid.INl+2;
	        i_max = Soln_ptr[i_blk].Grid.INl+2*Soln_Block_List.Block[i_blk].infoNW[0].dimen.ghost;
	        i_inc = 2;
	        j_min = Soln_ptr[i_blk].Grid.JNu-2*Soln_Block_List.Block[i_blk].infoNW[0].dimen.ghost;
	        j_max = Soln_ptr[i_blk].Grid.JNu-2;
	        j_inc = 2;
                i_ref = i_min-2;
                j_ref = j_max+2;
             } else if (Soln_Block_List.Block[i_blk].infoNW[0].dimen.j > 0) {
	        i_min = Soln_ptr[i_blk].Grid.INl+2*Soln_Block_List.Block[i_blk].infoNW[0].dimen.ghost;
	        i_max = Soln_ptr[i_blk].Grid.INl+2;
	        i_inc = -2;
	        j_min = Soln_ptr[i_blk].Grid.JNu-2*Soln_Block_List.Block[i_blk].infoNW[0].dimen.ghost;
	        j_max = Soln_ptr[i_blk].Grid.JNu-2;
	        j_inc = 2;
                i_ref = i_max-2;
                j_ref = j_max+2;
             } else if (Soln_Block_List.Block[i_blk].infoNW[0].dimen.i > 0) {
	        i_min = Soln_ptr[i_blk].Grid.INl+2;
	        i_max = Soln_ptr[i_blk].Grid.INl+2*Soln_Block_List.Block[i_blk].infoNW[0].dimen.ghost;
	        i_inc = 1;
	        j_min = Soln_ptr[i_blk].Grid.JNu-2;
	        j_max = Soln_ptr[i_blk].Grid.JNu-2*Soln_Block_List.Block[i_blk].infoNW[0].dimen.ghost;
	        j_inc = -2;
                i_ref = i_min-2;
                j_ref = j_min+2;
             } else {
	        i_min = Soln_ptr[i_blk].Grid.INl+2*Soln_Block_List.Block[i_blk].infoNW[0].dimen.ghost;
	        i_max = Soln_ptr[i_blk].Grid.INl+2;
	        i_inc = -2;
	        j_min = Soln_ptr[i_blk].Grid.JNu-2;
	        j_max = Soln_ptr[i_blk].Grid.JNu-2*Soln_Block_List.Block[i_blk].infoNW[0].dimen.ghost;
	        j_inc = -2;
                i_ref = i_max-2;
                j_ref = j_min+2;
             } /* endif */
             x_ref = Soln_ptr[i_blk].Grid.Node[i_ref][j_ref].X; // Reference node location.
             for ( j  = j_min ; ((j_inc+2)/4) ? (j <= j_max):(j >= j_max) ; j += j_inc ) {
                for ( i = i_min ;  ((i_inc+2)/4) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
                   for ( k = 1 ; k <= NUM_COMP_VECTOR2D; ++ k) {
		      l = l + 1;
                      if (l >= buffer_size_neighbour) return(4201);
                      Soln_Block_List.message_reschange_northwestcorner_sendbuf[i_blk][l] =
                         Soln_ptr[i_blk].Grid.Node[i][j].X[k]-x_ref[k];
                   } /* endfor */
                } /* endfor */
             } /* endfor */
          } /* endif */
       } /* endif */

       /////////////////////////////////////////////////////////////////////////
       // Load send buffers of solution blocks with coarser North East corner //
       // neighbour at lower mesh resolution.                                 //
       /////////////////////////////////////////////////////////////////////////
       if (Soln_Block_List.Block[i_blk].used && 
           (Soln_Block_List.Block[i_blk].nNE == 1) &&
           (Soln_Block_List.Block[i_blk].info.level > 
            Soln_Block_List.Block[i_blk].infoNE[0].level)) {
          buffer_size_neighbour = sqr(Soln_Block_List.Block[i_blk].infoNE[0].dimen.ghost)*
                                  Number_of_Solution_Variables;
          l = -1;
          // Load ghost cell solution variable information as required.
          if (!Send_Mesh_Geometry_Only && !Send_Conservative_Solution_Fluxes) {
             if (Soln_Block_List.Block[i_blk].infoNE[0].dimen.i > 0 &&
                 Soln_Block_List.Block[i_blk].infoNE[0].dimen.j > 0) {
	        i_min = Soln_ptr[i_blk].ICu-2*Soln_Block_List.Block[i_blk].infoNE[0].dimen.ghost+1;
                i_max = Soln_ptr[i_blk].ICu;
	        i_inc = 2;
	        j_min = Soln_ptr[i_blk].JCu-2*Soln_Block_List.Block[i_blk].infoNE[0].dimen.ghost+1;
	        j_max = Soln_ptr[i_blk].JCu;
	        j_inc = 2;
             } else if (Soln_Block_List.Block[i_blk].infoNE[0].dimen.j > 0) {
	        i_min = Soln_ptr[i_blk].ICu;
	        i_max = Soln_ptr[i_blk].ICu-2*Soln_Block_List.Block[i_blk].infoNE[0].dimen.ghost+1;
	        i_inc = -2;
	        j_min = Soln_ptr[i_blk].JCu-2*Soln_Block_List.Block[i_blk].infoNE[0].dimen.ghost+1;
	        j_max = Soln_ptr[i_blk].JCu;
	        j_inc = 2;
             } else if (Soln_Block_List.Block[i_blk].infoNE[0].dimen.i > 0) {
	        i_min = Soln_ptr[i_blk].ICu-2*Soln_Block_List.Block[i_blk].infoNE[0].dimen.ghost+1;
	        i_max = Soln_ptr[i_blk].ICu;
	        i_inc = 2;
	        j_min = Soln_ptr[i_blk].JCu;
	        j_max = Soln_ptr[i_blk].JCu-2*Soln_Block_List.Block[i_blk].infoNE[0].dimen.ghost+1;
	        j_inc = -2;
             } else {
	        i_min = Soln_ptr[i_blk].ICu;
	        i_max = Soln_ptr[i_blk].ICu-2*Soln_Block_List.Block[i_blk].infoNE[0].dimen.ghost+1;
	        i_inc = -2;
	        j_min = Soln_ptr[i_blk].JCu;
 	        j_max = Soln_ptr[i_blk].JCu-2*Soln_Block_List.Block[i_blk].infoNE[0].dimen.ghost+1;
	        j_inc = -2;
             } /* endif */
             i = Soln_ptr[i_blk].LoadSendBuffer_F2C(Soln_Block_List.message_reschange_northeastcorner_sendbuf[i_blk],
                                                    l,buffer_size_neighbour,
                                                    i_min,i_max,i_inc,
                                                    j_min,j_max,j_inc);
	     if (i != 0) return(5200);
          } /* endif */
          // Load ghost cell mesh information as required.
          if (Send_Mesh_Geometry_Only ||
              Number_of_Solution_Variables > Soln_ptr[i_blk].NumVar()) {
             if (Soln_Block_List.Block[i_blk].infoNE[0].dimen.i > 0 &&
                 Soln_Block_List.Block[i_blk].infoNE[0].dimen.j > 0) {
	        i_min = Soln_ptr[i_blk].Grid.INu-2*Soln_Block_List.Block[i_blk].infoNE[0].dimen.ghost;
	        i_max = Soln_ptr[i_blk].Grid.INu-2;
	        i_inc = 2;
	        j_min = Soln_ptr[i_blk].Grid.JNu-2*Soln_Block_List.Block[i_blk].infoNE[0].dimen.ghost;
	        j_max = Soln_ptr[i_blk].Grid.JNu-2;
	        j_inc = 2;
                i_ref = i_max+2;
                j_ref = j_max+2;
             } else if (Soln_Block_List.Block[i_blk].infoNE[0].dimen.j > 0) {
	        i_min = Soln_ptr[i_blk].Grid.INu-2;
	        i_max = Soln_ptr[i_blk].Grid.INu-2*Soln_Block_List.Block[i_blk].infoNE[0].dimen.ghost;
	        i_inc = -2;
	        j_min = Soln_ptr[i_blk].Grid.JNu-2*Soln_Block_List.Block[i_blk].infoNE[0].dimen.ghost;
	        j_max = Soln_ptr[i_blk].Grid.JNu-2;
	        j_inc = 2;
                i_ref = i_min+2;
                j_ref = j_max+2;
             } else if (Soln_Block_List.Block[i_blk].infoNE[0].dimen.i > 0) {
	        i_min = Soln_ptr[i_blk].Grid.INu-2*Soln_Block_List.Block[i_blk].infoNE[0].dimen.ghost;
	        i_max = Soln_ptr[i_blk].Grid.INu-2;
	        i_inc = 2;
	        j_min = Soln_ptr[i_blk].Grid.JNu-2;
	        j_max = Soln_ptr[i_blk].Grid.JNu-2*Soln_Block_List.Block[i_blk].infoNE[0].dimen.ghost;
	        j_inc = -2;
                i_ref = i_max+2;
                j_ref = j_min+2;
             } else {
	        i_min = Soln_ptr[i_blk].Grid.INu-2;
	        i_max = Soln_ptr[i_blk].Grid.INu-2*Soln_Block_List.Block[i_blk].infoNE[0].dimen.ghost;
	        i_inc = -2;
	        j_min = Soln_ptr[i_blk].Grid.JNu-2;
	        j_max = Soln_ptr[i_blk].Grid.JNu-2*Soln_Block_List.Block[i_blk].infoNE[0].dimen.ghost;
	        j_inc = -2;
                i_ref = i_min+2;
                j_ref = j_min+2;
             } /* endif */
             x_ref = Soln_ptr[i_blk].Grid.Node[i_ref][j_ref].X; // Reference node location.
             for ( j  = j_min ; ((j_inc+2)/4) ? (j <= j_max):(j >= j_max) ; j += j_inc ) {
                for ( i = i_min ;  ((i_inc+2)/4) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
                   for ( k = 1 ; k <= NUM_COMP_VECTOR2D; ++ k) {
		      l = l + 1;
                      if (l >= buffer_size_neighbour) return(5201);
                      Soln_Block_List.message_reschange_northeastcorner_sendbuf[i_blk][l] =
                         Soln_ptr[i_blk].Grid.Node[i][j].X[k]-x_ref[k];
                   } /* endfor */
                } /* endfor */
             } /* endfor */
          } /* endif */
       } /* endif */

       /////////////////////////////////////////////////////////////////////////
       // Load send buffers of solution blocks with coarser South East corner //
       // neighbour at lower mesh resolution.                                 //
       /////////////////////////////////////////////////////////////////////////
       if (Soln_Block_List.Block[i_blk].used && 
           (Soln_Block_List.Block[i_blk].nSE == 1) &&
           (Soln_Block_List.Block[i_blk].info.level > 
            Soln_Block_List.Block[i_blk].infoSE[0].level)) {
          buffer_size_neighbour = sqr(Soln_Block_List.Block[i_blk].infoSE[0].dimen.ghost)*
                                  Number_of_Solution_Variables;
          l = -1;
          // Load ghost cell solution variable information as required.
          if (!Send_Mesh_Geometry_Only && !Send_Conservative_Solution_Fluxes) {
             if (Soln_Block_List.Block[i_blk].infoSE[0].dimen.i > 0 &&
                 Soln_Block_List.Block[i_blk].infoSE[0].dimen.j > 0) {
	        i_min = Soln_ptr[i_blk].ICu-2*Soln_Block_List.Block[i_blk].infoSE[0].dimen.ghost+1;
                i_max = Soln_ptr[i_blk].ICu;
	        i_inc = 2;
	        j_min = Soln_ptr[i_blk].JCl;
	        j_max = Soln_ptr[i_blk].JCl+2*Soln_Block_List.Block[i_blk].infoSE[0].dimen.ghost-1;
	        j_inc = 2;
             } else if (Soln_Block_List.Block[i_blk].infoSE[0].dimen.j > 0) {
	        i_min = Soln_ptr[i_blk].ICu;
	        i_max = Soln_ptr[i_blk].ICu-2*Soln_Block_List.Block[i_blk].infoSE[0].dimen.ghost+1;
	        i_inc = -2;
	        j_min = Soln_ptr[i_blk].JCl;
	        j_max = Soln_ptr[i_blk].JCl+2*Soln_Block_List.Block[i_blk].infoSE[0].dimen.ghost-1;
	        j_inc = 2;
             } else if (Soln_Block_List.Block[i_blk].infoSE[0].dimen.i > 0) {
	        i_min = Soln_ptr[i_blk].ICu-2*Soln_Block_List.Block[i_blk].infoSE[0].dimen.ghost+1;
	        i_max = Soln_ptr[i_blk].ICu;
	        i_inc = 2;
	        j_min = Soln_ptr[i_blk].JCl+2*Soln_Block_List.Block[i_blk].infoSE[0].dimen.ghost-1;
	        j_max = Soln_ptr[i_blk].JCl;
	        j_inc = -2;
             } else {
	        i_min = Soln_ptr[i_blk].ICu;
	        i_max = Soln_ptr[i_blk].ICu-2*Soln_Block_List.Block[i_blk].infoSE[0].dimen.ghost+1;
	        i_inc = -2;
	        j_min = Soln_ptr[i_blk].JCl+2*Soln_Block_List.Block[i_blk].infoSE[0].dimen.ghost-1;
	        j_max = Soln_ptr[i_blk].JCl;
	        j_inc = -2;
             } /* endif */
             i = Soln_ptr[i_blk].LoadSendBuffer_F2C(Soln_Block_List.message_reschange_southeastcorner_sendbuf[i_blk],
                                                    l,buffer_size_neighbour,
                                                    i_min,i_max,i_inc,
                                                    j_min,j_max,j_inc);
	     if (i != 0) return(4100);
          } /* endif */
          // Load ghost cell mesh information as required.
          if (Send_Mesh_Geometry_Only ||
              Number_of_Solution_Variables > Soln_ptr[i_blk].NumVar()) {
             if (Soln_Block_List.Block[i_blk].infoSE[0].dimen.i > 0 &&
                 Soln_Block_List.Block[i_blk].infoSE[0].dimen.j > 0) {
	        i_min = Soln_ptr[i_blk].Grid.INu-2*Soln_Block_List.Block[i_blk].infoSE[0].dimen.ghost;
	        i_max = Soln_ptr[i_blk].Grid.INu-2;
                i_inc = 2;
	        j_min = Soln_ptr[i_blk].Grid.JNl+2;
	        j_max = Soln_ptr[i_blk].Grid.JNl+2*Soln_Block_List.Block[i_blk].infoSE[0].dimen.ghost;
	        j_inc = 2;
                i_ref = i_max+2;
                j_ref = j_min-2;
             } else if (Soln_Block_List.Block[i_blk].infoSE[0].dimen.j > 0) {
	        i_min = Soln_ptr[i_blk].Grid.INu-2;
	        i_max = Soln_ptr[i_blk].Grid.INu-2*Soln_Block_List.Block[i_blk].infoSE[0].dimen.ghost;
                i_inc = -2;
	        j_min = Soln_ptr[i_blk].Grid.JNl+2;
	        j_max = Soln_ptr[i_blk].Grid.JNl+2*Soln_Block_List.Block[i_blk].infoSE[0].dimen.ghost;
	        j_inc = 2;
                i_ref = i_min+2;
                j_ref = j_min-2;
             } else if (Soln_Block_List.Block[i_blk].infoSE[0].dimen.i > 0) {
	        i_min = Soln_ptr[i_blk].Grid.INu-2*Soln_Block_List.Block[i_blk].infoSE[0].dimen.ghost;
	        i_max = Soln_ptr[i_blk].Grid.INu-2;
                i_inc = 2;
	        j_min = Soln_ptr[i_blk].Grid.JNl+2*Soln_Block_List.Block[i_blk].infoSE[0].dimen.ghost;
	        j_max = Soln_ptr[i_blk].Grid.JNl+2;
	        j_inc = -2;
                i_ref = i_max+2;
                j_ref = j_max-2;
             } else {
	        i_min = Soln_ptr[i_blk].Grid.INu-2;
	        i_max = Soln_ptr[i_blk].Grid.INu-2*Soln_Block_List.Block[i_blk].infoSE[0].dimen.ghost;
                i_inc = -2;
	        j_min = Soln_ptr[i_blk].Grid.JNl+2*Soln_Block_List.Block[i_blk].infoSE[0].dimen.ghost;
	        j_max = Soln_ptr[i_blk].Grid.JNl+2;
	        j_inc = -2;
                i_ref = i_min+2;
                j_ref = j_max-2;
             } /* endif */
             x_ref = Soln_ptr[i_blk].Grid.Node[i_ref][j_ref].X; // Reference node location.
             for ( j  = j_min ; ((j_inc+2)/4) ? (j <= j_max):(j >= j_max) ; j += j_inc ) {
                for ( i = i_min ;  ((i_inc+2)/4) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
                   for ( k = 1 ; k <= NUM_COMP_VECTOR2D; ++ k) {
		      l = l + 1;
                      if (l >= buffer_size_neighbour) return(4101);
                      Soln_Block_List.message_reschange_southeastcorner_sendbuf[i_blk][l] =
                         Soln_ptr[i_blk].Grid.Node[i][j].X[k]-x_ref[k];
                   } /* endfor */
                } /* endfor */
             } /* endfor */
          } /* endif */
       } /* endif */

       /////////////////////////////////////////////////////////////////////////
       // Load send buffers of solution blocks with coarser South West corner //
       // neighbour at lower mesh resolution.                                 //
       /////////////////////////////////////////////////////////////////////////
       if (Soln_Block_List.Block[i_blk].used && 
           (Soln_Block_List.Block[i_blk].nSW == 1) &&
           (Soln_Block_List.Block[i_blk].info.level > 
            Soln_Block_List.Block[i_blk].infoSW[0].level)) {
          buffer_size_neighbour = sqr(Soln_Block_List.Block[i_blk].infoSW[0].dimen.ghost)*
                                  Number_of_Solution_Variables;
          l = -1;
          // Load ghost cell solution variable information as required.
          if (!Send_Mesh_Geometry_Only && !Send_Conservative_Solution_Fluxes) {
             if (Soln_Block_List.Block[i_blk].infoSW[0].dimen.i > 0 &&
                 Soln_Block_List.Block[i_blk].infoSW[0].dimen.j > 0) {
	        i_min = Soln_ptr[i_blk].ICl;
                i_max = Soln_ptr[i_blk].ICl+2*Soln_Block_List.Block[i_blk].infoSW[0].dimen.ghost-1;
	        i_inc = 2;
	        j_min = Soln_ptr[i_blk].JCl;
	        j_max = Soln_ptr[i_blk].JCl+2*Soln_Block_List.Block[i_blk].infoSW[0].dimen.ghost-1;
	        j_inc = 2;
             } else if (Soln_Block_List.Block[i_blk].infoSW[0].dimen.j > 0) {
	        i_min = Soln_ptr[i_blk].ICl+2*Soln_Block_List.Block[i_blk].infoSW[0].dimen.ghost-1;
	        i_max = Soln_ptr[i_blk].ICl;
	        i_inc = -2;
	        j_min = Soln_ptr[i_blk].JCl;
	        j_max = Soln_ptr[i_blk].JCl+2*Soln_Block_List.Block[i_blk].infoSW[0].dimen.ghost-1;
	        j_inc = 2;
             } else if (Soln_Block_List.Block[i_blk].infoSW[0].dimen.i > 0) {
	        i_min = Soln_ptr[i_blk].ICl;
                i_max = Soln_ptr[i_blk].ICl+2*Soln_Block_List.Block[i_blk].infoSW[0].dimen.ghost-1;
	        i_inc = 2;
	        j_min = Soln_ptr[i_blk].JCl+2*Soln_Block_List.Block[i_blk].infoSW[0].dimen.ghost-1;
	        j_max = Soln_ptr[i_blk].JCl;
	        j_inc = -2;
             } else {
	        i_min = Soln_ptr[i_blk].ICl+2*Soln_Block_List.Block[i_blk].infoSW[0].dimen.ghost-1;
	        i_max = Soln_ptr[i_blk].ICl;
	        i_inc = -2;
	        j_min = Soln_ptr[i_blk].JCl+2*Soln_Block_List.Block[i_blk].infoSW[0].dimen.ghost-1;
	        j_max = Soln_ptr[i_blk].JCl;
	        j_inc = -2;
             } /* endif */
             i = Soln_ptr[i_blk].LoadSendBuffer_F2C(Soln_Block_List.message_reschange_southwestcorner_sendbuf[i_blk],
                                                    l,buffer_size_neighbour,
                                                    i_min,i_max,i_inc,
                                                    j_min,j_max,j_inc);
	     if (i != 0) return(5100);
          } /* endif */
          // Load ghost cell mesh information as required.
          if (Send_Mesh_Geometry_Only ||
              Number_of_Solution_Variables > Soln_ptr[i_blk].NumVar()) {
             if (Soln_Block_List.Block[i_blk].infoSW[0].dimen.i > 0 &&
                 Soln_Block_List.Block[i_blk].infoSW[0].dimen.j > 0) {
	        i_min = Soln_ptr[i_blk].Grid.INl+2;
	        i_max = Soln_ptr[i_blk].Grid.INl+2*Soln_Block_List.Block[i_blk].infoSW[0].dimen.ghost;
                i_inc = 2;
	        j_min = Soln_ptr[i_blk].Grid.JNl+2;
	        j_max = Soln_ptr[i_blk].Grid.JNl+2*Soln_Block_List.Block[i_blk].infoSW[0].dimen.ghost;
	        j_inc = 2;
                i_ref = i_min-2;
                j_ref = j_min-2;
             } else if (Soln_Block_List.Block[i_blk].infoSW[0].dimen.j > 0) {
	        i_min = Soln_ptr[i_blk].Grid.INl+2*Soln_Block_List.Block[i_blk].infoSW[0].dimen.ghost;
	        i_max = Soln_ptr[i_blk].Grid.INl+2;
                i_inc = -2;
	        j_min = Soln_ptr[i_blk].Grid.JNl+2;
	        j_max = Soln_ptr[i_blk].Grid.JNl+2*Soln_Block_List.Block[i_blk].infoSW[0].dimen.ghost;
	        j_inc = 2;
                i_ref = i_max-2;
                j_ref = j_min-2;
             } else if (Soln_Block_List.Block[i_blk].infoSW[0].dimen.i > 0) {
	        i_min = Soln_ptr[i_blk].Grid.INl+2;
	        i_max = Soln_ptr[i_blk].Grid.INl+2*Soln_Block_List.Block[i_blk].infoSW[0].dimen.ghost;
                i_inc = 2;
	        j_min = Soln_ptr[i_blk].Grid.JNl+2*Soln_Block_List.Block[i_blk].infoSW[0].dimen.ghost;
	        j_max = Soln_ptr[i_blk].Grid.JNl+2;
	        j_inc = -2;
                i_ref = i_min-2;
                j_ref = j_max-2;
             } else {
	        i_min = Soln_ptr[i_blk].Grid.INl+2*Soln_Block_List.Block[i_blk].infoSW[0].dimen.ghost;
	        i_max = Soln_ptr[i_blk].Grid.INl+2;
                i_inc = -2;
	        j_min = Soln_ptr[i_blk].Grid.JNl+2*Soln_Block_List.Block[i_blk].infoSW[0].dimen.ghost;
	        j_max = Soln_ptr[i_blk].Grid.JNl+2;
	        j_inc = -2;
                i_ref = i_max-2;
                j_ref = j_max-2;
             } /* endif */
             x_ref = Soln_ptr[i_blk].Grid.Node[i_ref][j_ref].X; // Reference node location.
             for ( j  = j_min ; ((j_inc+2)/4) ? (j <= j_max):(j >= j_max) ; j += j_inc ) {
                for ( i = i_min ;  ((i_inc+2)/4) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
                   for ( k = 1 ; k <= NUM_COMP_VECTOR2D; ++ k) {
		      l = l + 1;
                      if (l >= buffer_size_neighbour) return(5101);
                      Soln_Block_List.message_reschange_southwestcorner_sendbuf[i_blk][l] =
                         Soln_ptr[i_blk].Grid.Node[i][j].X[k]-x_ref[k];
                   } /* endfor */
                } /* endfor */
             } /* endfor */
          } /* endif */
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
template <class Quad_Soln_Block>
int Load_Send_Message_Buffers_ResChange_CoarseToFine(Quad_Soln_Block *Soln_ptr,
                                                     AdaptiveBlock2D_List &Soln_Block_List,
                                                     const int Number_of_Solution_Variables,
                                                     const int Send_Mesh_Geometry_Only) {

    int i_blk, buffer_size_neighbour, 
        i_min0, i_max0, i_min1, i_max1, i_inc, i, 
        j_min0, j_max0, j_min1, j_max1, j_inc, j, 
        k, l0, l1;
    int i_ref0, j_ref0, i_ref1, j_ref1;
    Vector2D x_ref;

    /* Check to see if the number of solution variables specified
       in the calling argument is correct. */

    if (Number_of_Solution_Variables != NUM_COMP_VECTOR2D &&
        Number_of_Solution_Variables != Soln_ptr[0].NumVar() &&
        Number_of_Solution_Variables != Soln_ptr[0].NumVar() + NUM_COMP_VECTOR2D) {
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
          if (!Send_Mesh_Geometry_Only) {
             buffer_size_neighbour = Soln_Block_List.Block[i_blk].infoN[0].dimen.ghost*
                                     (abs(Soln_Block_List.Block[i_blk].infoN[0].dimen.i)+
                                     Soln_Block_List.Block[i_blk].infoN[0].dimen.ghost)*
                                     Soln_ptr[i_blk].NumVar();
             if (Number_of_Solution_Variables > Soln_ptr[i_blk].NumVar()) 
                buffer_size_neighbour = buffer_size_neighbour +
                                        Soln_Block_List.Block[i_blk].infoN[0].dimen.ghost*
                                        (abs(Soln_Block_List.Block[i_blk].infoN[0].dimen.i)+
                                        Soln_Block_List.Block[i_blk].infoN[0].dimen.ghost+1)*
                                        NUM_COMP_VECTOR2D+
                                        2*Soln_Block_List.Block[i_blk].infoN[0].dimen.ghost;
          } else {
             buffer_size_neighbour = Soln_Block_List.Block[i_blk].infoN[0].dimen.ghost*
                                     (abs(Soln_Block_List.Block[i_blk].infoN[0].dimen.i)+
                                     Soln_Block_List.Block[i_blk].infoN[0].dimen.ghost+1)*
                                     NUM_COMP_VECTOR2D+
                                     2*Soln_Block_List.Block[i_blk].infoN[0].dimen.ghost;
          } /* endif */
          l0 = -1;
          l1 = -1;
          // Load ghost cell solution variable information as required.
          if (!Send_Mesh_Geometry_Only) {
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
             i = Soln_ptr[i_blk].LoadSendBuffer_C2F(Soln_Block_List.message_reschange_northface_sendbuf[i_blk][0],
                                                    l0,buffer_size_neighbour,
                                                    i_min0,i_max0,i_inc,
                                                    j_min0,j_max0,j_inc);
	     if (i != 0) return(2200);
             //
             // Second neighbour (solution).
             //
             i = Soln_ptr[i_blk].LoadSendBuffer_C2F(Soln_Block_List.message_reschange_northface_sendbuf[i_blk][1],
                                                    l1,buffer_size_neighbour,
                                                    i_min1,i_max1,i_inc,
                                                    j_min1,j_max1,j_inc);
	     if (i != 0) return(2201);
          } /* endif */
          // Load ghost cell mesh information as required.
          if (Send_Mesh_Geometry_Only ||
              Number_of_Solution_Variables > Soln_ptr[i_blk].NumVar()) {
             if (Soln_Block_List.Block[i_blk].infoN[0].dimen.i > 0 &&
                 Soln_Block_List.Block[i_blk].infoN[0].dimen.j > 0) {
	        i_min0 = Soln_ptr[i_blk].Grid.INl;
	        i_max0 = Soln_ptr[i_blk].Grid.INl+(Soln_ptr[i_blk].Grid.INu-Soln_ptr[i_blk].Grid.INl)/2+1;
	        i_min1 = Soln_ptr[i_blk].Grid.INu-((Soln_ptr[i_blk].Grid.INu-Soln_ptr[i_blk].Grid.INl)/2+1);
	        i_max1 = Soln_ptr[i_blk].Grid.INu;
	        i_inc = 1;
	        j_min0 = Soln_ptr[i_blk].Grid.JNu-1;
	        j_max0 = Soln_ptr[i_blk].Grid.JNu-1;
	        j_min1 = Soln_ptr[i_blk].Grid.JNu-1;
	        j_max1 = Soln_ptr[i_blk].Grid.JNu-1;
	        j_inc = 1;
                i_ref0 = i_min0;
                j_ref0 = j_min0+1;
                i_ref1 = i_min1+1;
                j_ref1 = j_min1+1;
             } else if (Soln_Block_List.Block[i_blk].infoN[0].dimen.j > 0) {
	        i_min0 = Soln_ptr[i_blk].Grid.INl+(Soln_ptr[i_blk].Grid.INu-Soln_ptr[i_blk].Grid.INl)/2+1;
	        i_max0 = Soln_ptr[i_blk].Grid.INl;
	        i_min1 = Soln_ptr[i_blk].Grid.INu;
	        i_max1 = Soln_ptr[i_blk].Grid.INu-((Soln_ptr[i_blk].Grid.INu-Soln_ptr[i_blk].Grid.INl)/2+1);
	        i_inc = -1;
	        j_min0 = Soln_ptr[i_blk].Grid.JNu-1;
	        j_max0 = Soln_ptr[i_blk].Grid.JNu-1;
	        j_min1 = Soln_ptr[i_blk].Grid.JNu-1;
	        j_max1 = Soln_ptr[i_blk].Grid.JNu-1;
	        j_inc = 1;
                i_ref0 = i_min0-1;
                j_ref0 = j_min0+1;
                i_ref1 = i_min1;
                j_ref1 = j_min1+1;
             } else if (Soln_Block_List.Block[i_blk].infoN[0].dimen.i > 0) {
	        i_min0 = Soln_ptr[i_blk].Grid.INl;
	        i_max0 = Soln_ptr[i_blk].Grid.INl+(Soln_ptr[i_blk].Grid.INu-Soln_ptr[i_blk].Grid.INl)/2+1;
	        i_min1 = Soln_ptr[i_blk].Grid.INu-((Soln_ptr[i_blk].Grid.INu-Soln_ptr[i_blk].Grid.INl)/2+1);
	        i_max1 = Soln_ptr[i_blk].Grid.INu;
	        i_inc = 1;
	        j_min0 = Soln_ptr[i_blk].Grid.JNu-1;
	        j_max0 = Soln_ptr[i_blk].Grid.JNu-1;
	        j_min1 = Soln_ptr[i_blk].Grid.JNu-1;
	        j_max1 = Soln_ptr[i_blk].Grid.JNu-1;
	        j_inc = -1;
                i_ref0 = i_min0;
                j_ref0 = j_min0+1;
                i_ref1 = i_min1+1;
                j_ref1 = j_min1+1;
             } else {
	        i_min0 = Soln_ptr[i_blk].Grid.INl+(Soln_ptr[i_blk].Grid.INu-Soln_ptr[i_blk].Grid.INl)/2+1;
	        i_max0 = Soln_ptr[i_blk].Grid.INl;
	        i_min1 = Soln_ptr[i_blk].Grid.INu;
	        i_max1 = Soln_ptr[i_blk].Grid.INu-((Soln_ptr[i_blk].Grid.INu-Soln_ptr[i_blk].Grid.INl)/2+1);
	        i_inc = -1;
	        j_min0 = Soln_ptr[i_blk].Grid.JNu-1;
	        j_max0 = Soln_ptr[i_blk].Grid.JNu-1;
	        j_min1 = Soln_ptr[i_blk].Grid.JNu-1;
	        j_max1 = Soln_ptr[i_blk].Grid.JNu-1;
	        j_inc = -1;
                i_ref0 = i_min0-1;
                j_ref0 = j_min0+1;
                i_ref1 = i_min1;
                j_ref1 = j_min1+1;
             } /* endif */
             //
             // First neighbour (mesh).
             //
             x_ref = Soln_ptr[i_blk].Grid.Node[i_ref0][j_ref0].X; // Reference node location.
             // Four different orderings to consider depending on the value of i_inc & j_inc.
             if (j_inc > 0) {
                if (i_inc > 0) {
                   for ( i = i_min0 ;  ((i_inc+1)/2) ? (i <= i_max0):(i >= i_max0) ; i += i_inc ) {
	              // Evaluate SW sub (fine) cell SW node location.
                      for ( k = 1 ; k <= NUM_COMP_VECTOR2D; ++ k) {
                         l0 = l0 + 1;
                         if (l0 >= buffer_size_neighbour) return(2202);
                         Soln_Block_List.message_reschange_northface_sendbuf[i_blk][0][l0] =
                            Soln_ptr[i_blk].Grid.Node[i][j_min0].X[k]-x_ref[k];
                      } /* endfor */
	              // Evaluate SE sub (fine) cell SW node location.
                      if (i != i_max0) {
                         for ( k = 1 ; k <= NUM_COMP_VECTOR2D; ++ k) {
                            l0 = l0 + 1;
                            if (l0 >= buffer_size_neighbour) return(2203);
                            Soln_Block_List.message_reschange_northface_sendbuf[i_blk][0][l0] =
                               HALF*(Soln_ptr[i_blk].Grid.Node[i][j_min0].X[k]+
                                     Soln_ptr[i_blk].Grid.Node[i+1][j_min0].X[k])-x_ref[k];
                         } /* endfor */
                      } /* endif */
                   } /* endfor */
                   for ( i = i_min0 ;  ((i_inc+1)/2) ? (i <= i_max0):(i >= i_max0) ; i += i_inc ) {
	              // Evaluate NW sub (fine) cell SW node location.
                      for ( k = 1 ; k <= NUM_COMP_VECTOR2D; ++ k) {
                         l0 = l0 + 1;
                         if (l0 >= buffer_size_neighbour) return(2204);
                         Soln_Block_List.message_reschange_northface_sendbuf[i_blk][0][l0] =
                            HALF*(Soln_ptr[i_blk].Grid.Node[i][j_min0].X[k]+
                                  Soln_ptr[i_blk].Grid.Node[i][j_min0+1].X[k])-x_ref[k];
                      } /* endfor */
	              // Evaluate NE sub (fine) cell SW node location.
                      if (i != i_max0) {
                         for ( k = 1 ; k <= NUM_COMP_VECTOR2D; ++ k) {
                            l0 = l0 + 1;
                            if (l0 >= buffer_size_neighbour) return(2205);
                            Soln_Block_List.message_reschange_northface_sendbuf[i_blk][0][l0] =
                               Soln_ptr[i_blk].Grid.Cell[i][j_min0].Xc[k]-x_ref[k];
                         } /* endfor */
                      } /* endif */
                   } /* endfor */
                } else {
                   for ( i = i_min0 ;  ((i_inc+1)/2) ? (i <= i_max0):(i >= i_max0) ; i += i_inc ) {
	              // Evaluate SE sub (fine) cell SW node location.
                      if (i != i_min0) {
                         for ( k = 1 ; k <= NUM_COMP_VECTOR2D; ++ k) {
                            l0 = l0 + 1;
                            if (l0 >= buffer_size_neighbour) return(2206);
                            Soln_Block_List.message_reschange_northface_sendbuf[i_blk][0][l0] =
                               HALF*(Soln_ptr[i_blk].Grid.Node[i][j_min0].X[k]+
                                     Soln_ptr[i_blk].Grid.Node[i+1][j_min0].X[k])-x_ref[k];
                         } /* endfor */
                      } /* endif */
	              // Evaluate SW sub (fine) cell SW node location.
                      for ( k = 1 ; k <= NUM_COMP_VECTOR2D; ++ k) {
                         l0 = l0 + 1;
                         if (l0 >= buffer_size_neighbour) return(2207);
                         Soln_Block_List.message_reschange_northface_sendbuf[i_blk][0][l0] =
                            Soln_ptr[i_blk].Grid.Node[i][j_min0].X[k]-x_ref[k];
                      } /* endfor */
                   } /* endfor */
                   for ( i = i_min0 ;  ((i_inc+1)/2) ? (i <= i_max0):(i >= i_max0) ; i += i_inc ) {
	              // Evaluate NE sub (fine) cell SW node location.
                      if (i != i_min0) {
                         for ( k = 1 ; k <= NUM_COMP_VECTOR2D; ++ k) {
                            l0 = l0 + 1;
                            if (l0 >= buffer_size_neighbour) return(2208);
                            Soln_Block_List.message_reschange_northface_sendbuf[i_blk][0][l0] =
                               Soln_ptr[i_blk].Grid.Cell[i][j_min0].Xc[k]-x_ref[k];
                         } /* endfor */
                      } /* endif */
	              // Evaluate NW sub (fine) cell SW node location.
                      for ( k = 1 ; k <= NUM_COMP_VECTOR2D; ++ k) {
                         l0 = l0 + 1;
                         if (l0 >= buffer_size_neighbour) return(2209);
                         Soln_Block_List.message_reschange_northface_sendbuf[i_blk][0][l0] =
                            HALF*(Soln_ptr[i_blk].Grid.Node[i][j_min0].X[k]+
                                  Soln_ptr[i_blk].Grid.Node[i][j_min0+1].X[k])-x_ref[k];
                      } /* endfor */
                   } /* endfor */
                } /* endif */
             } else {
                if (i_inc > 0) {
                   for ( i = i_min0 ;  ((i_inc+1)/2) ? (i <= i_max0):(i >= i_max0) ; i += i_inc ) {
	              // Evaluate NW sub (fine) cell SW node location.
                      for ( k = 1 ; k <= NUM_COMP_VECTOR2D; ++ k) {
                         l0 = l0 + 1;
                         if (l0 >= buffer_size_neighbour) return(2210);
                         Soln_Block_List.message_reschange_northface_sendbuf[i_blk][0][l0] =
                            HALF*(Soln_ptr[i_blk].Grid.Node[i][j_min0].X[k]+
                                  Soln_ptr[i_blk].Grid.Node[i][j_min0+1].X[k])-x_ref[k];
                      } /* endfor */
	              // Evaluate NE sub (fine) cell SW node location.
                      if (i != i_max0) {
                         for ( k = 1 ; k <= NUM_COMP_VECTOR2D; ++ k) {
                            l0 = l0 + 1;
                            if (l0 >= buffer_size_neighbour) return(2211);
                            Soln_Block_List.message_reschange_northface_sendbuf[i_blk][0][l0] =
                               Soln_ptr[i_blk].Grid.Cell[i][j_min0].Xc[k]-x_ref[k];
                         } /* endfor */
                      } /* endif */
                   } /* endfor */
                   for ( i = i_min0 ;  ((i_inc+1)/2) ? (i <= i_max0):(i >= i_max0) ; i += i_inc ) {
	              // Evaluate SW sub (fine) cell SW node location.
                      for ( k = 1 ; k <= NUM_COMP_VECTOR2D; ++ k) {
                         l0 = l0 + 1;
                         if (l0 >= buffer_size_neighbour) return(2212);
                         Soln_Block_List.message_reschange_northface_sendbuf[i_blk][0][l0] =
                            Soln_ptr[i_blk].Grid.Node[i][j_min0].X[k]-x_ref[k];
                      } /* endfor */
	              // Evaluate SE sub (fine) cell SW node location.
                      if (i != i_max0) {
                         for ( k = 1 ; k <= NUM_COMP_VECTOR2D; ++ k) {
                            l0 = l0 + 1;
                            if (l0 >= buffer_size_neighbour) return(2213);
                            Soln_Block_List.message_reschange_northface_sendbuf[i_blk][0][l0] =
                               HALF*(Soln_ptr[i_blk].Grid.Node[i][j_min0].X[k]+
                                     Soln_ptr[i_blk].Grid.Node[i+1][j_min0].X[k])-x_ref[k];
                         } /* endfor */
                      } /* endif */
                   } /* endfor */
                } else {
                   for ( i = i_min0 ;  ((i_inc+1)/2) ? (i <= i_max0):(i >= i_max0) ; i += i_inc ) {
	              // Evaluate NE sub (fine) cell SW node location.
                      if (i != i_min0) {
                         for ( k = 1 ; k <= NUM_COMP_VECTOR2D; ++ k) {
                            l0 = l0 + 1;
                            if (l0 >= buffer_size_neighbour) return(2214);
                            Soln_Block_List.message_reschange_northface_sendbuf[i_blk][0][l0] =
                               Soln_ptr[i_blk].Grid.Cell[i][j_min0].Xc[k]-x_ref[k];
                         } /* endfor */
                      } /* endif */
	              // Evaluate NW sub (fine) cell SW node location.
                      for ( k = 1 ; k <= NUM_COMP_VECTOR2D; ++ k) {
                         l0 = l0 + 1;
                         if (l0 >= buffer_size_neighbour) return(2215);
                         Soln_Block_List.message_reschange_northface_sendbuf[i_blk][0][l0] =
                            HALF*(Soln_ptr[i_blk].Grid.Node[i][j_min0].X[k]+
                                  Soln_ptr[i_blk].Grid.Node[i][j_min0+1].X[k])-x_ref[k];
                      } /* endfor */
                   } /* endfor */
                   for ( i = i_min0 ;  ((i_inc+1)/2) ? (i <= i_max0):(i >= i_max0) ; i += i_inc ) {
	              // Evaluate SE sub (fine) cell SW node location.
                      if (i != i_min0) {
                         for ( k = 1 ; k <= NUM_COMP_VECTOR2D; ++ k) {
                            l0 = l0 + 1;
                            if (l0 >= buffer_size_neighbour) return(2216);
                            Soln_Block_List.message_reschange_northface_sendbuf[i_blk][0][l0] =
                               HALF*(Soln_ptr[i_blk].Grid.Node[i][j_min0].X[k]+
                                     Soln_ptr[i_blk].Grid.Node[i+1][j_min0].X[k])-x_ref[k];
                         } /* endfor */
                      } /* endif */
	              // Evaluate SW sub (fine) cell SW node location.
                      for ( k = 1 ; k <= NUM_COMP_VECTOR2D; ++ k) {
                         l0 = l0 + 1;
                         if (l0 >= buffer_size_neighbour) return(2217);
                         Soln_Block_List.message_reschange_northface_sendbuf[i_blk][0][l0] =
                            Soln_ptr[i_blk].Grid.Node[i][j_min0].X[k]-x_ref[k];
                      } /* endfor */
                   } /* endfor */
                } /* endif */
             } /* endif */
             //
             // Second neighbour (mesh).
             //
             x_ref = Soln_ptr[i_blk].Grid.Node[i_ref1][j_ref1].X; // Reference node location.
             if (j_inc > 0) {
                if (i_inc > 0) {
                   for ( i = i_min1 ;  ((i_inc+1)/2) ? (i <= i_max1):(i >= i_max1) ; i += i_inc ) {
	              // Evaluate SW sub (fine) cell SW node location.
                      for ( k = 1 ; k <= NUM_COMP_VECTOR2D; ++ k) {
                         l1 = l1 + 1;
                         if (l1 >= buffer_size_neighbour) return(2218);
                         Soln_Block_List.message_reschange_northface_sendbuf[i_blk][1][l1] =
                            Soln_ptr[i_blk].Grid.Node[i][j_min1].X[k]-x_ref[k];
                      } /* endfor */
                      // Evaluate SE sub (fine) cell SW node location.
                      if (i != i_max1) {
                         for ( k = 1 ; k <= NUM_COMP_VECTOR2D; ++ k) {
                            l1 = l1 + 1;
                            if (l1 >= buffer_size_neighbour) return(2219);
                            Soln_Block_List.message_reschange_northface_sendbuf[i_blk][1][l1] =
                               HALF*(Soln_ptr[i_blk].Grid.Node[i][j_min1].X[k]+
                                     Soln_ptr[i_blk].Grid.Node[i+1][j_min1].X[k])-x_ref[k];
                         } /* endfor */
                      } /* endif */
                   } /* endfor */
                   for ( i = i_min1 ;  ((i_inc+1)/2) ? (i <= i_max1):(i >= i_max1) ; i += i_inc ) {
	              // Evaluate NW sub (fine) cell SW node location.
                      for ( k = 1 ; k <= NUM_COMP_VECTOR2D; ++ k) {
                         l1 = l1 + 1;
                         if (l1 >= buffer_size_neighbour) return(2220);
                         Soln_Block_List.message_reschange_northface_sendbuf[i_blk][1][l1] =
                            HALF*(Soln_ptr[i_blk].Grid.Node[i][j_min1].X[k]+
                                  Soln_ptr[i_blk].Grid.Node[i][j_min1+1].X[k])-x_ref[k];
                      } /* endfor */ 
                      // Evaluate NE sub (fine) cell SW node location.
                      if (i != i_max1) {
                         for ( k = 1 ; k <= NUM_COMP_VECTOR2D; ++ k) {
                            l1 = l1 + 1;
                            if (l1 >= buffer_size_neighbour) return(2221);
                            Soln_Block_List.message_reschange_northface_sendbuf[i_blk][1][l1] =
                               Soln_ptr[i_blk].Grid.Cell[i][j_min1].Xc[k]-x_ref[k];
                         } /* endfor */
                      } /* endif */
                   } /* endfor */
                } else {
                   for ( i = i_min1 ;  ((i_inc+1)/2) ? (i <= i_max1):(i >= i_max1) ; i += i_inc ) {
                      // Evaluate SE sub (fine) cell SW node location.
                      if (i != i_min1) {
                         for ( k = 1 ; k <= NUM_COMP_VECTOR2D; ++ k) {
                            l1 = l1 + 1;
                            if (l1 >= buffer_size_neighbour) return(2222);
                            Soln_Block_List.message_reschange_northface_sendbuf[i_blk][1][l1] =
                               HALF*(Soln_ptr[i_blk].Grid.Node[i][j_min1].X[k]+
                                     Soln_ptr[i_blk].Grid.Node[i+1][j_min1].X[k])-x_ref[k];
                         } /* endfor */
                      } /* endif */
	              // Evaluate SW sub (fine) cell SW node location.
                      for ( k = 1 ; k <= NUM_COMP_VECTOR2D; ++ k) {
                         l1 = l1 + 1;
                         if (l1 >= buffer_size_neighbour) return(2223);
                         Soln_Block_List.message_reschange_northface_sendbuf[i_blk][1][l1] =
                            Soln_ptr[i_blk].Grid.Node[i][j_min1].X[k]-x_ref[k];
                      } /* endfor */
                   } /* endfor */
                   for ( i = i_min1 ;  ((i_inc+1)/2) ? (i <= i_max1):(i >= i_max1) ; i += i_inc ) {
                      // Evaluate NE sub (fine) cell SW node location.
                      if (i != i_min1) {
                         for ( k = 1 ; k <= NUM_COMP_VECTOR2D; ++ k) {
                            l1 = l1 + 1;
                            if (l1 >= buffer_size_neighbour) return(2224);
                            Soln_Block_List.message_reschange_northface_sendbuf[i_blk][1][l1] =
                               Soln_ptr[i_blk].Grid.Cell[i][j_min1].Xc[k]-x_ref[k];
                         } /* endfor */
                      } /* endif */
	              // Evaluate NW sub (fine) cell SW node location.
                      for ( k = 1 ; k <= NUM_COMP_VECTOR2D; ++ k) {
                         l1 = l1 + 1;
                         if (l1 >= buffer_size_neighbour) return(2225);
                         Soln_Block_List.message_reschange_northface_sendbuf[i_blk][1][l1] =
                            HALF*(Soln_ptr[i_blk].Grid.Node[i][j_min1].X[k]+
                                  Soln_ptr[i_blk].Grid.Node[i][j_min1+1].X[k])-x_ref[k];
                      } /* endfor */
                   } /* endfor */
                } /* endif */
             } else {
                if (i_inc > 0) {
                   for ( i = i_min1 ;  ((i_inc+1)/2) ? (i <= i_max1):(i >= i_max1) ; i += i_inc ) {
	              // Evaluate NW sub (fine) cell SW node location.
                      for ( k = 1 ; k <= NUM_COMP_VECTOR2D; ++ k) {
                         l1 = l1 + 1;
                         if (l1 >= buffer_size_neighbour) return(2226);
                         Soln_Block_List.message_reschange_northface_sendbuf[i_blk][1][l1] =
                            HALF*(Soln_ptr[i_blk].Grid.Node[i][j_min1].X[k]+
                                  Soln_ptr[i_blk].Grid.Node[i][j_min1+1].X[k])-x_ref[k];
                      } /* endfor */ 
                      // Evaluate NE sub (fine) cell SW node location.
                      if (i != i_max1) {
                         for ( k = 1 ; k <= NUM_COMP_VECTOR2D; ++ k) {
                            l1 = l1 + 1;
                            if (l1 >= buffer_size_neighbour) return(2227);
                            Soln_Block_List.message_reschange_northface_sendbuf[i_blk][1][l1] =
                               Soln_ptr[i_blk].Grid.Cell[i][j_min1].Xc[k]-x_ref[k];
                         } /* endfor */
                      } /* endif */
                   } /* endfor */
                   for ( i = i_min1 ;  ((i_inc+1)/2) ? (i <= i_max1):(i >= i_max1) ; i += i_inc ) {
	              // Evaluate SW sub (fine) cell SW node location.
                      for ( k = 1 ; k <= NUM_COMP_VECTOR2D; ++ k) {
                         l1 = l1 + 1;
                         if (l1 >= buffer_size_neighbour) return(2228);
                         Soln_Block_List.message_reschange_northface_sendbuf[i_blk][1][l1] =
                            Soln_ptr[i_blk].Grid.Node[i][j_min1].X[k]-x_ref[k];
                      } /* endfor */
                      // Evaluate SE sub (fine) cell SW node location.
                      if (i != i_max1) {
                         for ( k = 1 ; k <= NUM_COMP_VECTOR2D; ++ k) {
                            l1 = l1 + 1;
                            if (l1 >= buffer_size_neighbour) return(2229);
                            Soln_Block_List.message_reschange_northface_sendbuf[i_blk][1][l1] =
                               HALF*(Soln_ptr[i_blk].Grid.Node[i][j_min1].X[k]+
                                     Soln_ptr[i_blk].Grid.Node[i+1][j_min1].X[k])-x_ref[k];
                         } /* endfor */
                      } /* endif */
                   } /* endfor */
                } else {
                   for ( i = i_min1 ;  ((i_inc+1)/2) ? (i <= i_max1):(i >= i_max1) ; i += i_inc ) {
                      // Evaluate NE sub (fine) cell SW node location.
                      if (i != i_min1) {
                         for ( k = 1 ; k <= NUM_COMP_VECTOR2D; ++ k) {
                            l1 = l1 + 1;
                            if (l1 >= buffer_size_neighbour) return(2230);
                            Soln_Block_List.message_reschange_northface_sendbuf[i_blk][1][l1] =
                               Soln_ptr[i_blk].Grid.Cell[i][j_min1].Xc[k]-x_ref[k];
                         } /* endfor */
                      } /* endif */
	              // Evaluate NW sub (fine) cell SW node location.
                      for ( k = 1 ; k <= NUM_COMP_VECTOR2D; ++ k) {
                         l1 = l1 + 1;
                         if (l1 >= buffer_size_neighbour) return(2231);
                         Soln_Block_List.message_reschange_northface_sendbuf[i_blk][1][l1] =
                            HALF*(Soln_ptr[i_blk].Grid.Node[i][j_min1].X[k]+
                                  Soln_ptr[i_blk].Grid.Node[i][j_min1+1].X[k])-x_ref[k];
                      } /* endfor */
                   } /* endfor */
                   for ( i = i_min1 ;  ((i_inc+1)/2) ? (i <= i_max1):(i >= i_max1) ; i += i_inc ) {
                      // Evaluate SE sub (fine) cell SW node location.
                      if (i != i_min1) {
                         for ( k = 1 ; k <= NUM_COMP_VECTOR2D; ++ k) {
                            l1 = l1 + 1;
                            if (l1 >= buffer_size_neighbour) return(2232);
                            Soln_Block_List.message_reschange_northface_sendbuf[i_blk][1][l1] =
                               HALF*(Soln_ptr[i_blk].Grid.Node[i][j_min1].X[k]+
                                     Soln_ptr[i_blk].Grid.Node[i+1][j_min1].X[k])-x_ref[k];
                         } /* endfor */
                      } /* endif */
	              // Evaluate SW sub (fine) cell SW node location.
                      for ( k = 1 ; k <= NUM_COMP_VECTOR2D; ++ k) {
                         l1 = l1 + 1;
                         if (l1 >= buffer_size_neighbour) return(2233);
                         Soln_Block_List.message_reschange_northface_sendbuf[i_blk][1][l1] =
                            Soln_ptr[i_blk].Grid.Node[i][j_min1].X[k]-x_ref[k];
                      } /* endfor */
                   } /* endfor */
                } /* endif */
             } /* endif */
             if (Soln_Block_List.Block[i_blk].infoN[0].dimen.i > 0 &&
                 Soln_Block_List.Block[i_blk].infoN[0].dimen.j > 0) {
	        i_inc = 1;
	        j_min0 = Soln_ptr[i_blk].JCu;
	        j_max0 = Soln_ptr[i_blk].JCu;
	        j_min1 = Soln_ptr[i_blk].JCu;
	        j_max1 = Soln_ptr[i_blk].JCu;
	        j_inc = 1;
             } else if (Soln_Block_List.Block[i_blk].infoN[0].dimen.j > 0) {
	        i_inc = -1;
	        j_min0 = Soln_ptr[i_blk].JCu;
	        j_max0 = Soln_ptr[i_blk].JCu;
	        j_min1 = Soln_ptr[i_blk].JCu;
	        j_max1 = Soln_ptr[i_blk].JCu;
	        j_inc = 1;
             } else if (Soln_Block_List.Block[i_blk].infoN[0].dimen.i > 0) {
	        i_inc = 1;
	        j_min0 = Soln_ptr[i_blk].JCu;
	        j_max0 = Soln_ptr[i_blk].JCu;
	        j_min1 = Soln_ptr[i_blk].JCu;
	        j_max1 = Soln_ptr[i_blk].JCu;
	        j_inc = -1;
             } else {
	        i_inc = -1;
	        j_min0 = Soln_ptr[i_blk].JCu;
	        j_max0 = Soln_ptr[i_blk].JCu;
	        j_min1 = Soln_ptr[i_blk].JCu;
	        j_max1 = Soln_ptr[i_blk].JCu;
	        j_inc = -1;
             } /* endif */
             if (i_inc > 0) {
                for ( k = 1 ; k <= 2; ++ k) {
                   for ( j  = j_min0 ; ((j_inc+1)/2) ? (j <= j_max0):(j >= j_max0) ; j += j_inc ) {
                      l0 = l0 + 2;
                      if (l0 >= buffer_size_neighbour) return(2234);
                      Soln_Block_List.message_reschange_northface_sendbuf[i_blk][0][l0-1] = 
                         double(Soln_ptr[i_blk].Grid.BCtypeW[j]);
                      Soln_Block_List.message_reschange_northface_sendbuf[i_blk][0][l0] = 
                         double(Soln_ptr[i_blk].Grid.BCtypeW[j]);
                   } /* endfor */
                   for ( j  = j_min1 ; ((j_inc+1)/2) ? (j <= j_max1):(j >= j_max1) ; j += j_inc ) {
                      l1 = l1 + 2;
                      if (l1 >= buffer_size_neighbour) return(2235);
                      Soln_Block_List.message_reschange_northface_sendbuf[i_blk][1][l1-1] = 
                         double(Soln_ptr[i_blk].Grid.BCtypeE[j]);
                      Soln_Block_List.message_reschange_northface_sendbuf[i_blk][1][l1] = 
                         double(Soln_ptr[i_blk].Grid.BCtypeE[j]);
                   } /* endfor */
                } /* endfor */
	     } else {
                for ( k = 1 ; k <= 2; ++ k) {
                   for ( j  = j_min0 ; ((j_inc+1)/2) ? (j <= j_max0):(j >= j_max0) ; j += j_inc ) {
                      l0 = l0 + 2;
                      if (l0 >= buffer_size_neighbour) return(2236);
                      Soln_Block_List.message_reschange_northface_sendbuf[i_blk][0][l0-1] = 
                         double(Soln_ptr[i_blk].Grid.BCtypeE[j]);
                      Soln_Block_List.message_reschange_northface_sendbuf[i_blk][0][l0] = 
                         double(Soln_ptr[i_blk].Grid.BCtypeE[j]);
                   } /* endfor */
                   for ( j  = j_min1 ; ((j_inc+1)/2) ? (j <= j_max1):(j >= j_max1) ; j += j_inc ) {
                      l1 = l1 + 1;
                      if (l1 >= buffer_size_neighbour) return(2237);
                      Soln_Block_List.message_reschange_northface_sendbuf[i_blk][1][l1-1] = 
                         double(Soln_ptr[i_blk].Grid.BCtypeW[j]);
                      Soln_Block_List.message_reschange_northface_sendbuf[i_blk][1][l1] = 
                         double(Soln_ptr[i_blk].Grid.BCtypeW[j]);
                   } /* endfor */
                } /* endfor */
             } /* endif */
          } /* endif */
       } /* endif */

       //////////////////////////////////////////////////////////////////
       // Load send buffers of solution blocks with refined South face //
       // neighbour at higher mesh resolution.                         //
       //////////////////////////////////////////////////////////////////
       if (Soln_Block_List.Block[i_blk].used && 
           (Soln_Block_List.Block[i_blk].nS == 2) &&
           (Soln_Block_List.Block[i_blk].info.level < 
            Soln_Block_List.Block[i_blk].infoS[0].level)) {
          if (!Send_Mesh_Geometry_Only) {
             buffer_size_neighbour = Soln_Block_List.Block[i_blk].infoS[0].dimen.ghost*
                                     (abs(Soln_Block_List.Block[i_blk].infoS[0].dimen.i)+
                                     Soln_Block_List.Block[i_blk].infoS[0].dimen.ghost)*
                                     Soln_ptr[i_blk].NumVar();
             if (Number_of_Solution_Variables > Soln_ptr[i_blk].NumVar()) 
                buffer_size_neighbour = buffer_size_neighbour +
                                        Soln_Block_List.Block[i_blk].infoS[0].dimen.ghost*
                                        (abs(Soln_Block_List.Block[i_blk].infoS[0].dimen.i)+
                                        Soln_Block_List.Block[i_blk].infoS[0].dimen.ghost+1)*
                                        NUM_COMP_VECTOR2D+
                                        2*Soln_Block_List.Block[i_blk].infoS[0].dimen.ghost;
          } else {
             buffer_size_neighbour = Soln_Block_List.Block[i_blk].infoS[0].dimen.ghost*
                                     (abs(Soln_Block_List.Block[i_blk].infoS[0].dimen.i)+
                                     Soln_Block_List.Block[i_blk].infoS[0].dimen.ghost+1)*
                                     NUM_COMP_VECTOR2D+
                                     2*Soln_Block_List.Block[i_blk].infoS[0].dimen.ghost;
          } /* endif */
          l0 = -1;
          l1 = -1;
          // Load ghost cell solution variable information as required.
          if (!Send_Mesh_Geometry_Only) {
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
             i = Soln_ptr[i_blk].LoadSendBuffer_C2F(Soln_Block_List.message_reschange_southface_sendbuf[i_blk][0],
                                                    l0,buffer_size_neighbour,
                                                    i_min0,i_max0,i_inc,
                                                    j_min0,j_max0,j_inc);
	     if (i != 0) return(2100);
             //
             // Second neighbour (solution).
             //
             i = Soln_ptr[i_blk].LoadSendBuffer_C2F(Soln_Block_List.message_reschange_southface_sendbuf[i_blk][1],
                                                    l1,buffer_size_neighbour,
                                                    i_min1,i_max1,i_inc,
                                                    j_min1,j_max1,j_inc);
	     if (i != 0) return(2101);
          } /* endif */
          // Load ghost cell mesh information as required.
          if (Send_Mesh_Geometry_Only ||
              Number_of_Solution_Variables > Soln_ptr[i_blk].NumVar()) {
             if (Soln_Block_List.Block[i_blk].infoS[0].dimen.i > 0 &&
                 Soln_Block_List.Block[i_blk].infoS[0].dimen.j > 0) {
	        i_min0 = Soln_ptr[i_blk].Grid.INl;
	        i_max0 = Soln_ptr[i_blk].Grid.INl+(Soln_ptr[i_blk].Grid.INu-Soln_ptr[i_blk].Grid.INl)/2+1;
	        i_min1 = Soln_ptr[i_blk].Grid.INu-((Soln_ptr[i_blk].Grid.INu-Soln_ptr[i_blk].Grid.INl)/2+1);
	        i_max1 = Soln_ptr[i_blk].Grid.INu;
	        i_inc = 1;
	        j_min0 = Soln_ptr[i_blk].Grid.JNl;
	        j_max0 = Soln_ptr[i_blk].Grid.JNl;
	        j_min1 = Soln_ptr[i_blk].Grid.JNl;
	        j_max1 = Soln_ptr[i_blk].Grid.JNl;
	        j_inc = 1;
                i_ref0 = i_min0;
                j_ref0 = j_min0;
                i_ref1 = i_min1+1;
                j_ref1 = j_min1;
             } else if (Soln_Block_List.Block[i_blk].infoS[0].dimen.j > 0) {
	        i_min0 = Soln_ptr[i_blk].Grid.INl+(Soln_ptr[i_blk].Grid.INu-Soln_ptr[i_blk].Grid.INl)/2+1;
	        i_max0 = Soln_ptr[i_blk].Grid.INl;
	        i_min1 = Soln_ptr[i_blk].Grid.INu;
	        i_max1 = Soln_ptr[i_blk].Grid.INu-((Soln_ptr[i_blk].Grid.INu-Soln_ptr[i_blk].Grid.INl)/2+1);
	        i_inc = -1;
	        j_min0 = Soln_ptr[i_blk].Grid.JNl;
	        j_max0 = Soln_ptr[i_blk].Grid.JNl;
	        j_min1 = Soln_ptr[i_blk].Grid.JNl;
	        j_max1 = Soln_ptr[i_blk].Grid.JNl;
	        j_inc = 1;
                i_ref0 = i_min0;
                j_ref0 = j_min0-1;
                i_ref1 = i_min1;
                j_ref1 = j_min1;
             } else if (Soln_Block_List.Block[i_blk].infoS[0].dimen.i > 0) {
	        i_min0 = Soln_ptr[i_blk].Grid.INl;
	        i_max0 = Soln_ptr[i_blk].Grid.INl+(Soln_ptr[i_blk].Grid.INu-Soln_ptr[i_blk].Grid.INl)/2+1;
	        i_min1 = Soln_ptr[i_blk].Grid.INu-((Soln_ptr[i_blk].Grid.INu-Soln_ptr[i_blk].Grid.INl)/2+1);
	        i_max1 = Soln_ptr[i_blk].Grid.INu;
	        i_inc = 1;
	        j_min0 = Soln_ptr[i_blk].Grid.JNl;
	        j_max0 = Soln_ptr[i_blk].Grid.JNl;
	        j_min1 = Soln_ptr[i_blk].Grid.JNl;
	        j_max1 = Soln_ptr[i_blk].Grid.JNl;
	        j_inc = -1;
                i_ref0 = i_min0;
                j_ref0 = j_min0;
                i_ref1 = i_min1+1;
                j_ref1 = j_min1;
             } else {
	        i_min0 = Soln_ptr[i_blk].Grid.INl+(Soln_ptr[i_blk].Grid.INu-Soln_ptr[i_blk].Grid.INl)/2+1;
	        i_max0 = Soln_ptr[i_blk].Grid.INl;
	        i_min1 = Soln_ptr[i_blk].Grid.INu;
	        i_max1 = Soln_ptr[i_blk].Grid.INu-((Soln_ptr[i_blk].Grid.INu-Soln_ptr[i_blk].Grid.INl)/2+1);
	        i_inc = -1;
	        j_min0 = Soln_ptr[i_blk].Grid.JNl;
	        j_max0 = Soln_ptr[i_blk].Grid.JNl;
	        j_min1 = Soln_ptr[i_blk].Grid.JNl;
	        j_max1 = Soln_ptr[i_blk].Grid.JNl;
	        j_inc = -1;
                i_ref0 = i_min0-1;
                j_ref0 = j_min0;
                i_ref1 = i_min1;
                j_ref1 = j_min1;
             } /* endif */
             //
             // First neighbour (mesh).
             //
             x_ref = Soln_ptr[i_blk].Grid.Node[i_ref0][j_ref0].X; // Reference node location.
             // Four different orderings to consider depending on the value of i_inc & j_inc.
             if (j_inc > 0) {
                if (i_inc > 0) {
                   for ( i = i_min0 ;  ((i_inc+1)/2) ? (i <= i_max0):(i >= i_max0) ; i += i_inc ) {
	              // Evaluate SW sub (fine) cell NW node location.
                      for ( k = 1 ; k <= NUM_COMP_VECTOR2D; ++ k) {
                         l0 = l0 + 1;
                         if (l0 >= buffer_size_neighbour) return(2102);
                         Soln_Block_List.message_reschange_southface_sendbuf[i_blk][0][l0] =
                            HALF*(Soln_ptr[i_blk].Grid.Node[i][j_min0].X[k]+
                                  Soln_ptr[i_blk].Grid.Node[i][j_min0+1].X[k])-x_ref[k];
                      } /* endfor */
	              // Evaluate SE sub (fine) cell NW node location.
                      if (i != i_max0) {
                         for ( k = 1 ; k <= NUM_COMP_VECTOR2D; ++ k) {
                            l0 = l0 + 1;
                            if (l0 >= buffer_size_neighbour) return(2103);
                            Soln_Block_List.message_reschange_southface_sendbuf[i_blk][0][l0] =
                               Soln_ptr[i_blk].Grid.Cell[i][j_min0].Xc[k]-x_ref[k];
                         } /* endfor */
                      } /* endif */
                   } /* endfor */
                   for ( i = i_min0 ;  ((i_inc+1)/2) ? (i <= i_max0):(i >= i_max0) ; i += i_inc ) {
	              // Evaluate NW sub (fine) cell NW node location.
                      for ( k = 1 ; k <= NUM_COMP_VECTOR2D; ++ k) {
                         l0 = l0 + 1;
                         if (l0 >= buffer_size_neighbour) return(2104);
                         Soln_Block_List.message_reschange_southface_sendbuf[i_blk][0][l0] =
                            Soln_ptr[i_blk].Grid.Node[i][j_min0+1].X[k]-x_ref[k];
                      } /* endfor */
	              // Evaluate NE sub (fine) cell NW node location.
                      if (i != i_max0) {
                         for ( k = 1 ; k <= NUM_COMP_VECTOR2D; ++ k) {
                            l0 = l0 + 1;
                            if (l0 >= buffer_size_neighbour) return(2105);
                            Soln_Block_List.message_reschange_southface_sendbuf[i_blk][0][l0] =
                               HALF*(Soln_ptr[i_blk].Grid.Node[i][j_min0+1].X[k]+
                                     Soln_ptr[i_blk].Grid.Node[i+1][j_min0+1].X[k])-x_ref[k];
                         } /* endfor */
                      } /* endif */
                   } /* endfor */
                } else {
                   for ( i = i_min0 ;  ((i_inc+1)/2) ? (i <= i_max0):(i >= i_max0) ; i += i_inc ) {
	              // Evaluate SE sub (fine) cell NW node location.
                      if (i != i_min0) {
                         for ( k = 1 ; k <= NUM_COMP_VECTOR2D; ++ k) {
                            l0 = l0 + 1;
                            if (l0 >= buffer_size_neighbour) return(2106);
                            Soln_Block_List.message_reschange_southface_sendbuf[i_blk][0][l0] =
                               Soln_ptr[i_blk].Grid.Cell[i][j_min0].Xc[k]-x_ref[k];
                         } /* endfor */
                      } /* endif */
	              // Evaluate SW sub (fine) cell NW node location.
                      for ( k = 1 ; k <= NUM_COMP_VECTOR2D; ++ k) {
                         l0 = l0 + 1;
                         if (l0 >= buffer_size_neighbour) return(2107);
                         Soln_Block_List.message_reschange_southface_sendbuf[i_blk][0][l0] =
                            HALF*(Soln_ptr[i_blk].Grid.Node[i][j_min0].X[k]+
                                  Soln_ptr[i_blk].Grid.Node[i][j_min0+1].X[k])-x_ref[k];
                      } /* endfor */
                   } /* endfor */
                   for ( i = i_min0 ;  ((i_inc+1)/2) ? (i <= i_max0):(i >= i_max0) ; i += i_inc ) {
	              // Evaluate NE sub (fine) cell NW node location.
                      if (i != i_min0) {
                         for ( k = 1 ; k <= NUM_COMP_VECTOR2D; ++ k) {
                            l0 = l0 + 1;
                            if (l0 >= buffer_size_neighbour) return(2108);
                            Soln_Block_List.message_reschange_southface_sendbuf[i_blk][0][l0] =
                               HALF*(Soln_ptr[i_blk].Grid.Node[i][j_min0+1].X[k]+
                                     Soln_ptr[i_blk].Grid.Node[i+1][j_min0+1].X[k])-x_ref[k];
                         } /* endfor */
                      } /* endif */
	              // Evaluate NW sub (fine) cell NW node location.
                      for ( k = 1 ; k <= NUM_COMP_VECTOR2D; ++ k) {
                         l0 = l0 + 1;
                         if (l0 >= buffer_size_neighbour) return(2109);
                         Soln_Block_List.message_reschange_southface_sendbuf[i_blk][0][l0] =
                            Soln_ptr[i_blk].Grid.Node[i][j_min0+1].X[k]-x_ref[k];
                      } /* endfor */
                   } /* endfor */
                } /* endif */
             } else {
                if (i_inc > 0) {
                   for ( i = i_min0 ;  ((i_inc+1)/2) ? (i <= i_max0):(i >= i_max0) ; i += i_inc ) {
	              // Evaluate NW sub (fine) cell NW node location.
                      for ( k = 1 ; k <= NUM_COMP_VECTOR2D; ++ k) {
                         l0 = l0 + 1;
                         if (l0 >= buffer_size_neighbour) return(2110);
                         Soln_Block_List.message_reschange_southface_sendbuf[i_blk][0][l0] =
                            Soln_ptr[i_blk].Grid.Node[i][j_min0+1].X[k]-x_ref[k];
                      } /* endfor */
	              // Evaluate NE sub (fine) cell NW node location.
                      if (i != i_max0) {
                         for ( k = 1 ; k <= NUM_COMP_VECTOR2D; ++ k) {
                            l0 = l0 + 1;
                            if (l0 >= buffer_size_neighbour) return(2111);
                            Soln_Block_List.message_reschange_southface_sendbuf[i_blk][0][l0] =
                               HALF*(Soln_ptr[i_blk].Grid.Node[i][j_min0+1].X[k]+
                                     Soln_ptr[i_blk].Grid.Node[i+1][j_min0+1].X[k])-x_ref[k];
                         } /* endfor */
                      } /* endif */
                   } /* endfor */
                   for ( i = i_min0 ;  ((i_inc+1)/2) ? (i <= i_max0):(i >= i_max0) ; i += i_inc ) {
	              // Evaluate SW sub (fine) cell NW node location.
                      for ( k = 1 ; k <= NUM_COMP_VECTOR2D; ++ k) {
                         l0 = l0 + 1;
                         if (l0 >= buffer_size_neighbour) return(2112);
                         Soln_Block_List.message_reschange_southface_sendbuf[i_blk][0][l0] =
                            HALF*(Soln_ptr[i_blk].Grid.Node[i][j_min0].X[k]+
                                  Soln_ptr[i_blk].Grid.Node[i][j_min0+1].X[k])-x_ref[k];
                      } /* endfor */
	              // Evaluate SE sub (fine) cell NW node location.
                      if (i != i_max0) {
                         for ( k = 1 ; k <= NUM_COMP_VECTOR2D; ++ k) {
                            l0 = l0 + 1;
                            if (l0 >= buffer_size_neighbour) return(2113);
                            Soln_Block_List.message_reschange_southface_sendbuf[i_blk][0][l0] =
                               Soln_ptr[i_blk].Grid.Cell[i][j_min0].Xc[k]-x_ref[k];
                         } /* endfor */
                      } /* endif */
                   } /* endfor */
                } else {
                   for ( i = i_min0 ;  ((i_inc+1)/2) ? (i <= i_max0):(i >= i_max0) ; i += i_inc ) {
	              // Evaluate NE sub (fine) cell NW node location.
                      if (i != i_min0) {
                         for ( k = 1 ; k <= NUM_COMP_VECTOR2D; ++ k) {
                            l0 = l0 + 1;
                            if (l0 >= buffer_size_neighbour) return(2114);
                            Soln_Block_List.message_reschange_southface_sendbuf[i_blk][0][l0] =
                               HALF*(Soln_ptr[i_blk].Grid.Node[i][j_min0+1].X[k]+
                                     Soln_ptr[i_blk].Grid.Node[i+1][j_min0+1].X[k])-x_ref[k];
                         } /* endfor */
                      } /* endif */
	              // Evaluate NW sub (fine) cell NW node location.
                      for ( k = 1 ; k <= NUM_COMP_VECTOR2D; ++ k) {
                         l0 = l0 + 1;
                         if (l0 >= buffer_size_neighbour) return(2115);
                         Soln_Block_List.message_reschange_southface_sendbuf[i_blk][0][l0] =
                            Soln_ptr[i_blk].Grid.Node[i][j_min0+1].X[k]-x_ref[k];
                      } /* endfor */
                   } /* endfor */
                   for ( i = i_min0 ;  ((i_inc+1)/2) ? (i <= i_max0):(i >= i_max0) ; i += i_inc ) {
	              // Evaluate SE sub (fine) cell NW node location.
                      if (i != i_min0) {
                         for ( k = 1 ; k <= NUM_COMP_VECTOR2D; ++ k) {
                            l0 = l0 + 1;
                            if (l0 >= buffer_size_neighbour) return(2116);
                            Soln_Block_List.message_reschange_southface_sendbuf[i_blk][0][l0] =
                               Soln_ptr[i_blk].Grid.Cell[i][j_min0].Xc[k]-x_ref[k];
                         } /* endfor */
                      } /* endif */
	              // Evaluate SW sub (fine) cell NW node location.
                      for ( k = 1 ; k <= NUM_COMP_VECTOR2D; ++ k) {
                         l0 = l0 + 1;
                         if (l0 >= buffer_size_neighbour) return(2117);
                         Soln_Block_List.message_reschange_southface_sendbuf[i_blk][0][l0] =
                            HALF*(Soln_ptr[i_blk].Grid.Node[i][j_min0].X[k]+
                                  Soln_ptr[i_blk].Grid.Node[i][j_min0+1].X[k])-x_ref[k];
                      } /* endfor */
                   } /* endfor */
                } /* endif */
             } /* endif */
             //
             // Second neighbour (mesh).
             //
             x_ref = Soln_ptr[i_blk].Grid.Node[i_ref1][j_ref1].X; // Reference node location.
             if (j_inc > 0) {
                if (i_inc > 0) {
                   for ( i = i_min1 ;  ((i_inc+1)/2) ? (i <= i_max1):(i >= i_max1) ; i += i_inc ) {
	              // Evaluate SW sub (fine) cell NW node location.
                      for ( k = 1 ; k <= NUM_COMP_VECTOR2D; ++ k) {
                         l1 = l1 + 1;
                         if (l1 >= buffer_size_neighbour) return(2118);
                         Soln_Block_List.message_reschange_southface_sendbuf[i_blk][1][l1] =
                            HALF*(Soln_ptr[i_blk].Grid.Node[i][j_min1].X[k]+
                                  Soln_ptr[i_blk].Grid.Node[i][j_min1+1].X[k])-x_ref[k];
                      } /* endfor */
                      // Evaluate SE sub (fine) cell NW node location.
                      if (i != i_max1) {
                         for ( k = 1 ; k <= NUM_COMP_VECTOR2D; ++ k) {
                            l1 = l1 + 1;
                            if (l1 >= buffer_size_neighbour) return(2119);
                            Soln_Block_List.message_reschange_southface_sendbuf[i_blk][1][l1] =
                               Soln_ptr[i_blk].Grid.Cell[i][j_min1].Xc[k]-x_ref[k];
                         } /* endfor */
                      } /* endif */
                   } /* endfor */
                   for ( i = i_min1 ;  ((i_inc+1)/2) ? (i <= i_max1):(i >= i_max1) ; i += i_inc ) {
	              // Evaluate NW sub (fine) cell NW node location.
                      for ( k = 1 ; k <= NUM_COMP_VECTOR2D; ++ k) {
                         l1 = l1 + 1;
                         if (l1 >= buffer_size_neighbour) return(2120);
                         Soln_Block_List.message_reschange_southface_sendbuf[i_blk][1][l1] =
                            Soln_ptr[i_blk].Grid.Node[i][j_min1+1].X[k]-x_ref[k];
                      } /* endfor */ 
                      // Evaluate NE sub (fine) cell NW node location.
                      if (i != i_max1) {
                         for ( k = 1 ; k <= NUM_COMP_VECTOR2D; ++ k) {
                            l1 = l1 + 1;
                            if (l1 >= buffer_size_neighbour) return(2121);
                            Soln_Block_List.message_reschange_southface_sendbuf[i_blk][1][l1] =
                               HALF*(Soln_ptr[i_blk].Grid.Node[i][j_min1+1].X[k]+
                                     Soln_ptr[i_blk].Grid.Node[i+1][j_min1+1].X[k])-x_ref[k];
                         } /* endfor */
                      } /* endif */
                   } /* endfor */
                } else {
                   for ( i = i_min1 ;  ((i_inc+1)/2) ? (i <= i_max1):(i >= i_max1) ; i += i_inc ) {
                      // Evaluate SE sub (fine) cell NW node location.
                      if (i != i_min1) {
                         for ( k = 1 ; k <= NUM_COMP_VECTOR2D; ++ k) {
                            l1 = l1 + 1;
                            if (l1 >= buffer_size_neighbour) return(2122);
                            Soln_Block_List.message_reschange_southface_sendbuf[i_blk][1][l1] =
                               Soln_ptr[i_blk].Grid.Cell[i][j_min1].Xc[k]-x_ref[k];
                         } /* endfor */
                      } /* endif */
	              // Evaluate SW sub (fine) cell NW node location.
                      for ( k = 1 ; k <= NUM_COMP_VECTOR2D; ++ k) {
                         l1 = l1 + 1;
                         if (l1 >= buffer_size_neighbour) return(2123);
                         Soln_Block_List.message_reschange_southface_sendbuf[i_blk][1][l1] =
                            HALF*(Soln_ptr[i_blk].Grid.Node[i][j_min1].X[k]+
                                  Soln_ptr[i_blk].Grid.Node[i][j_min1+1].X[k])-x_ref[k];
                      } /* endfor */
                   } /* endfor */
                   for ( i = i_min1 ;  ((i_inc+1)/2) ? (i <= i_max1):(i >= i_max1) ; i += i_inc ) {
                       // Evaluate NE sub (fine) cell NW node location.
                      if (i != i_min1) {
                         for ( k = 1 ; k <= NUM_COMP_VECTOR2D; ++ k) {
                            l1 = l1 + 1;
                            if (l1 >= buffer_size_neighbour) return(2124);
                            Soln_Block_List.message_reschange_southface_sendbuf[i_blk][1][l1] =
                               HALF*(Soln_ptr[i_blk].Grid.Node[i][j_min1+1].X[k]+
                                     Soln_ptr[i_blk].Grid.Node[i+1][j_min1+1].X[k])-x_ref[k];
                         } /* endfor */
                      } /* endif */
	              // Evaluate NW sub (fine) cell NW node location.
                      for ( k = 1 ; k <= NUM_COMP_VECTOR2D; ++ k) {
                         l1 = l1 + 1;
                         if (l1 >= buffer_size_neighbour) return(2125);
                         Soln_Block_List.message_reschange_southface_sendbuf[i_blk][1][l1] =
                            Soln_ptr[i_blk].Grid.Node[i][j_min1+1].X[k]-x_ref[k];
                      } /* endfor */
                   } /* endfor */
                } /* endif */
             } else {
                if (i_inc > 0) {
                   for ( i = i_min1 ;  ((i_inc+1)/2) ? (i <= i_max1):(i >= i_max1) ; i += i_inc ) {
	              // Evaluate NW sub (fine) cell NW node location.
                      for ( k = 1 ; k <= NUM_COMP_VECTOR2D; ++ k) {
                         l1 = l1 + 1;
                         if (l1 >= buffer_size_neighbour) return(2126);
                         Soln_Block_List.message_reschange_southface_sendbuf[i_blk][1][l1] =
                            Soln_ptr[i_blk].Grid.Cell[i][j_min1+1].Xc[k]-x_ref[k];
                      } /* endfor */ 
                      // Evaluate NE sub (fine) cell NW node location.
                      if (i != i_max1) {
                         for ( k = 1 ; k <= NUM_COMP_VECTOR2D; ++ k) {
                            l1 = l1 + 1;
                            if (l1 >= buffer_size_neighbour) return(2127);
                            Soln_Block_List.message_reschange_southface_sendbuf[i_blk][1][l1] =
                               HALF*(Soln_ptr[i_blk].Grid.Node[i][j_min1+1].X[k]+
                                     Soln_ptr[i_blk].Grid.Node[i+1][j_min1+1].X[k])-x_ref[k];
                         } /* endfor */
                      } /* endif */
                   } /* endfor */
                   for ( i = i_min1 ;  ((i_inc+1)/2) ? (i <= i_max1):(i >= i_max1) ; i += i_inc ) {
	              // Evaluate SW sub (fine) cell NW node location.
                      for ( k = 1 ; k <= NUM_COMP_VECTOR2D; ++ k) {
                         l1 = l1 + 1;
                         if (l1 >= buffer_size_neighbour) return(2128);
                         Soln_Block_List.message_reschange_southface_sendbuf[i_blk][1][l1] =
                            HALF*(Soln_ptr[i_blk].Grid.Node[i][j_min1].X[k]+
                                  Soln_ptr[i_blk].Grid.Node[i][j_min1+1].X[k])-x_ref[k];
                      } /* endfor */
                      // Evaluate SE sub (fine) cell NW node location.
                      if (i != i_max1) {
                         for ( k = 1 ; k <= NUM_COMP_VECTOR2D; ++ k) {
                            l1 = l1 + 1;
                            if (l1 >= buffer_size_neighbour) return(2129);
                            Soln_Block_List.message_reschange_southface_sendbuf[i_blk][1][l1] =
                               Soln_ptr[i_blk].Grid.Node[i][j_min1+1].X[k]-x_ref[k];
                         } /* endfor */
                      } /* endif */
                   } /* endfor */
                } else {
                   for ( i = i_min1 ;  ((i_inc+1)/2) ? (i <= i_max1):(i >= i_max1) ; i += i_inc ) {
                       // Evaluate NE sub (fine) cell NW node location.
                      if (i != i_min1) {
                         for ( k = 1 ; k <= NUM_COMP_VECTOR2D; ++ k) {
                            l1 = l1 + 1;
                            if (l1 >= buffer_size_neighbour) return(2130);
                            Soln_Block_List.message_reschange_southface_sendbuf[i_blk][1][l1] =
                               HALF*(Soln_ptr[i_blk].Grid.Node[i][j_min1+1].X[k]+
                                     Soln_ptr[i_blk].Grid.Node[i+1][j_min1+1].X[k])-x_ref[k];
                         } /* endfor */
                      } /* endif */
	              // Evaluate NW sub (fine) cell NW node location.
                      for ( k = 1 ; k <= NUM_COMP_VECTOR2D; ++ k) {
                         l1 = l1 + 1;
                         if (l1 >= buffer_size_neighbour) return(2131);
                         Soln_Block_List.message_reschange_southface_sendbuf[i_blk][1][l1] =
                            Soln_ptr[i_blk].Grid.Node[i][j_min1+1].X[k]-x_ref[k];
                      } /* endfor */
                   } /* endfor */
                   for ( i = i_min1 ;  ((i_inc+1)/2) ? (i <= i_max1):(i >= i_max1) ; i += i_inc ) {
                      // Evaluate SE sub (fine) cell NW node location.
                      if (i != i_min1) {
                         for ( k = 1 ; k <= NUM_COMP_VECTOR2D; ++ k) {
                            l1 = l1 + 1;
                            if (l1 >= buffer_size_neighbour) return(2132);
                            Soln_Block_List.message_reschange_southface_sendbuf[i_blk][1][l1] =
                               Soln_ptr[i_blk].Grid.Cell[i][j_min1].Xc[k]-x_ref[k];
                         } /* endfor */
                      } /* endif */
	              // Evaluate SW sub (fine) cell NW node location.
                      for ( k = 1 ; k <= NUM_COMP_VECTOR2D; ++ k) {
                         l1 = l1 + 1;
                         if (l1 >= buffer_size_neighbour) return(2133);
                         Soln_Block_List.message_reschange_southface_sendbuf[i_blk][1][l1] =
                            HALF*(Soln_ptr[i_blk].Grid.Node[i][j_min1].X[k]+
                                  Soln_ptr[i_blk].Grid.Node[i][j_min1+1].X[k])-x_ref[k];
                      } /* endfor */
                   } /* endfor */
                } /* endif */
             } /* endif */
             if (Soln_Block_List.Block[i_blk].infoS[0].dimen.i > 0 &&
                 Soln_Block_List.Block[i_blk].infoS[0].dimen.j > 0) {
	        i_inc = 1;
	        j_min0 = Soln_ptr[i_blk].JCl;
	        j_max0 = Soln_ptr[i_blk].JCl;
	        j_min1 = Soln_ptr[i_blk].JCl;
	        j_max1 = Soln_ptr[i_blk].JCl;
	        j_inc = 1;
             } else if (Soln_Block_List.Block[i_blk].infoS[0].dimen.j > 0) {
	        i_inc = -1;
	        j_min0 = Soln_ptr[i_blk].JCl;
	        j_max0 = Soln_ptr[i_blk].JCl;
	        j_min1 = Soln_ptr[i_blk].JCl;
	        j_max1 = Soln_ptr[i_blk].JCl;
	        j_inc = 1;
             } else if (Soln_Block_List.Block[i_blk].infoS[0].dimen.i > 0) {
	        i_inc = 1;
	        j_min0 = Soln_ptr[i_blk].JCl;
	        j_max0 = Soln_ptr[i_blk].JCl;
	        j_min1 = Soln_ptr[i_blk].JCl;
	        j_max1 = Soln_ptr[i_blk].JCl;
	        j_inc = -1;
             } else {
	        i_inc = -1;
	        j_min0 = Soln_ptr[i_blk].JCl;
	        j_max0 = Soln_ptr[i_blk].JCl;
	        j_min1 = Soln_ptr[i_blk].JCl;
	        j_max1 = Soln_ptr[i_blk].JCl;
	        j_inc = -1;
             } /* endif */
             if (i_inc > 0) {
                for ( k = 1 ; k <= 2; ++ k) {
                   for ( j  = j_min0 ; ((j_inc+1)/2) ? (j <= j_max0):(j >= j_max0) ; j += j_inc ) {
                      l0 = l0 + 2;
                      if (l0 >= buffer_size_neighbour) return(2134);
                      Soln_Block_List.message_reschange_southface_sendbuf[i_blk][0][l0-1] = 
                         double(Soln_ptr[i_blk].Grid.BCtypeW[j]);
                      Soln_Block_List.message_reschange_southface_sendbuf[i_blk][0][l0] = 
                         double(Soln_ptr[i_blk].Grid.BCtypeW[j]);
                   } /* endfor */
                   for ( j  = j_min1 ; ((j_inc+1)/2) ? (j <= j_max1):(j >= j_max1) ; j += j_inc ) {
                      l1 = l1 + 2;
                      if (l1 >= buffer_size_neighbour) return(2135);
                      Soln_Block_List.message_reschange_southface_sendbuf[i_blk][1][l1-1] = 
                         double(Soln_ptr[i_blk].Grid.BCtypeE[j]);
                      Soln_Block_List.message_reschange_southface_sendbuf[i_blk][1][l1] = 
                         double(Soln_ptr[i_blk].Grid.BCtypeE[j]);
                   } /* endfor */
                } /* endfor */
	     } else {
                for ( k = 1 ; k <= 2; ++ k) {
                   for ( j  = j_min0 ; ((j_inc+1)/2) ? (j <= j_max0):(j >= j_max0) ; j += j_inc ) {
                      l0 = l0 + 2;
                      if (l0 >= buffer_size_neighbour) return(2136);
                      Soln_Block_List.message_reschange_southface_sendbuf[i_blk][0][l0-1] = 
                         double(Soln_ptr[i_blk].Grid.BCtypeE[j]);
                      Soln_Block_List.message_reschange_southface_sendbuf[i_blk][0][l0] = 
                         double(Soln_ptr[i_blk].Grid.BCtypeE[j]);
                   } /* endfor */
                   for ( j  = j_min1 ; ((j_inc+1)/2) ? (j <= j_max1):(j >= j_max1) ; j += j_inc ) {
                      l1 = l1 + 1;
                      if (l1 >= buffer_size_neighbour) return(2137);
                      Soln_Block_List.message_reschange_southface_sendbuf[i_blk][1][l1-1] = 
                         double(Soln_ptr[i_blk].Grid.BCtypeW[j]);
                      Soln_Block_List.message_reschange_southface_sendbuf[i_blk][1][l1] = 
                         double(Soln_ptr[i_blk].Grid.BCtypeW[j]);
                   } /* endfor */
                } /* endfor */
             } /* endif */
          } /* endif */
       } /* endif */

       //////////////////////////////////////////////////////////////////
       // Load send buffers of solution blocks with refined East face //
       // neighbour at higher mesh resolution.                         //
       //////////////////////////////////////////////////////////////////
       if (Soln_Block_List.Block[i_blk].used && 
           (Soln_Block_List.Block[i_blk].nE == 2) &&
           (Soln_Block_List.Block[i_blk].info.level < 
            Soln_Block_List.Block[i_blk].infoE[0].level)) {
          if (!Send_Mesh_Geometry_Only) {
             buffer_size_neighbour = Soln_Block_List.Block[i_blk].infoE[0].dimen.ghost*
                                     (abs(Soln_Block_List.Block[i_blk].infoE[0].dimen.j)+
                                     Soln_Block_List.Block[i_blk].infoE[0].dimen.ghost)*
                                     Soln_ptr[i_blk].NumVar();
             if (Number_of_Solution_Variables > Soln_ptr[i_blk].NumVar()) 
                buffer_size_neighbour = buffer_size_neighbour +
                                        Soln_Block_List.Block[i_blk].infoE[0].dimen.ghost*
                                        (abs(Soln_Block_List.Block[i_blk].infoE[0].dimen.j)+
                                        Soln_Block_List.Block[i_blk].infoE[0].dimen.ghost+1)*
                                        NUM_COMP_VECTOR2D+
                                        2*Soln_Block_List.Block[i_blk].infoE[0].dimen.ghost;
          } else {
             buffer_size_neighbour = Soln_Block_List.Block[i_blk].infoE[0].dimen.ghost*
                                     (abs(Soln_Block_List.Block[i_blk].infoE[0].dimen.j)+
                                     Soln_Block_List.Block[i_blk].infoE[0].dimen.ghost+1)*
                                     NUM_COMP_VECTOR2D+
                                     2*Soln_Block_List.Block[i_blk].infoE[0].dimen.ghost;
          } /* endif */
          l0 = -1;
          l1 = -1;
          // Load ghost cell solution variable information as required.
          if (!Send_Mesh_Geometry_Only) {
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
             i = Soln_ptr[i_blk].LoadSendBuffer_C2F(Soln_Block_List.message_reschange_eastface_sendbuf[i_blk][0],
                                                    l0,buffer_size_neighbour,
                                                    i_min0,i_max0,i_inc,
                                                    j_min0,j_max0,j_inc);
	     if (i != 0) return(1200);
             //
             // Second neighbour (solution).
             //
             i = Soln_ptr[i_blk].LoadSendBuffer_C2F(Soln_Block_List.message_reschange_eastface_sendbuf[i_blk][1],
                                                    l1,buffer_size_neighbour,
                                                    i_min1,i_max1,i_inc,
                                                    j_min1,j_max1,j_inc);
	     if (i != 0) return(1201);
          } /* endif */
          // Load ghost cell mesh information as required.
          if (Send_Mesh_Geometry_Only ||
              Number_of_Solution_Variables > Soln_ptr[i_blk].NumVar()) {
             if (Soln_Block_List.Block[i_blk].infoE[0].dimen.i > 0 &&
                 Soln_Block_List.Block[i_blk].infoE[0].dimen.j > 0) {
	        i_min0 = Soln_ptr[i_blk].Grid.INu-1;
	        i_max0 = Soln_ptr[i_blk].Grid.INu-1;
	        i_min1 = Soln_ptr[i_blk].Grid.INu-1;
	        i_max1 = Soln_ptr[i_blk].Grid.INu-1;
	        i_inc = 1;
	        j_min0 = Soln_ptr[i_blk].Grid.JNl;
	        j_max0 = Soln_ptr[i_blk].Grid.JNl+(Soln_ptr[i_blk].Grid.JNu-Soln_ptr[i_blk].Grid.JNl)/2+1;
	        j_min1 = Soln_ptr[i_blk].Grid.JNu-((Soln_ptr[i_blk].Grid.JNu-Soln_ptr[i_blk].Grid.JNl)/2+1);
	        j_max1 = Soln_ptr[i_blk].Grid.JNu;
	        j_inc = 1;
                i_ref0 = i_min1+1;
                j_ref0 = j_min0;
                i_ref1 = i_min1+1;
                j_ref1 = j_min1+1;
             } else if (Soln_Block_List.Block[i_blk].infoE[0].dimen.j > 0) {
	        i_min0 = Soln_ptr[i_blk].Grid.INu-1;
	        i_max0 = Soln_ptr[i_blk].Grid.INu-1;
	        i_min1 = Soln_ptr[i_blk].Grid.INu-1;
	        i_max1 = Soln_ptr[i_blk].Grid.INu-1;
	        i_inc = -1;
	        j_min0 = Soln_ptr[i_blk].Grid.JNl;
	        j_max0 = Soln_ptr[i_blk].Grid.JNl+(Soln_ptr[i_blk].Grid.JNu-Soln_ptr[i_blk].Grid.JNl)/2+1;
	        j_min1 = Soln_ptr[i_blk].Grid.JNu-((Soln_ptr[i_blk].Grid.JNu-Soln_ptr[i_blk].Grid.JNl)/2+1);
	        j_max1 = Soln_ptr[i_blk].Grid.JNu;
	        j_inc = 1;
                i_ref0 = i_min1+1;
                j_ref0 = j_min0;
                i_ref1 = i_min1+1;
                j_ref1 = j_min1+1;
             } else if (Soln_Block_List.Block[i_blk].infoE[0].dimen.i > 0) {
	        i_min0 = Soln_ptr[i_blk].Grid.INu-1;
	        i_max0 = Soln_ptr[i_blk].Grid.INu-1;
	        i_min1 = Soln_ptr[i_blk].Grid.INu-1;
	        i_max1 = Soln_ptr[i_blk].Grid.INu-1;
	        i_inc = 1;
	        j_min0 = Soln_ptr[i_blk].Grid.JNl+(Soln_ptr[i_blk].Grid.JNu-Soln_ptr[i_blk].Grid.JNl)/2+1;
	        j_max0 = Soln_ptr[i_blk].Grid.JNl;
	        j_min1 = Soln_ptr[i_blk].Grid.JNu;
	        j_max1 = Soln_ptr[i_blk].Grid.JNu-((Soln_ptr[i_blk].Grid.JNu-Soln_ptr[i_blk].Grid.JNl)/2+1);
	        j_inc = -1;
                i_ref0 = i_min1+1;
                j_ref0 = j_min0-1;
                i_ref1 = i_min1+1;
                j_ref1 = j_min1;
             } else {
	        i_min0 = Soln_ptr[i_blk].Grid.INu-1;
	        i_max0 = Soln_ptr[i_blk].Grid.INu-1;
	        i_min1 = Soln_ptr[i_blk].Grid.INu-1;
	        i_max1 = Soln_ptr[i_blk].Grid.INu-1;
	        i_inc = -1;
	        j_min0 = Soln_ptr[i_blk].Grid.JNl+(Soln_ptr[i_blk].Grid.JNu-Soln_ptr[i_blk].Grid.JNl)/2+1;
	        j_max0 = Soln_ptr[i_blk].Grid.JNl;
	        j_min1 = Soln_ptr[i_blk].Grid.JNu;
	        j_max1 = Soln_ptr[i_blk].Grid.JNu-((Soln_ptr[i_blk].Grid.JNu-Soln_ptr[i_blk].Grid.JNl)/2+1);
	        j_inc = -1;
                i_ref0 = i_min1+1;
                j_ref0 = j_min0-1;
                i_ref1 = i_min1+1;
                j_ref1 = j_min1;
             } /* endif */
             //
             // First neighbour (mesh).
             //
             x_ref = Soln_ptr[i_blk].Grid.Node[i_ref0][j_ref0].X; // Reference node location.
             // Four different orderings to consider depending on the value of i_inc & j_inc.
             if (j_inc > 0) {
                if (i_inc > 0) {
                   for ( j = j_min0 ; ((j_inc+1)/2) ? (j <= j_max0):(j >= j_max0) ; j += j_inc ) {
	              // Evaluate SW sub (fine) cell SW node location.
                      for ( k = 1 ; k <= NUM_COMP_VECTOR2D; ++ k) {
                         l0 = l0 + 1;
                         if (l0 >= buffer_size_neighbour) return(1202);
                         Soln_Block_List.message_reschange_eastface_sendbuf[i_blk][0][l0] =
                            Soln_ptr[i_blk].Grid.Node[i_min0][j].X[k]-x_ref[k];
                      } /* endfor */
	              // Evaluate SE sub (fine) cell SW node location.
                      for ( k = 1 ; k <= NUM_COMP_VECTOR2D; ++ k) {
                         l0 = l0 + 1;
                         if (l0 >= buffer_size_neighbour) return(1203);
                         Soln_Block_List.message_reschange_eastface_sendbuf[i_blk][0][l0] =
                            HALF*(Soln_ptr[i_blk].Grid.Node[i_min0][j].X[k]+
                                  Soln_ptr[i_blk].Grid.Node[i_min0+1][j].X[k])-x_ref[k];
                      } /* endfor */
                      if (j != j_max0) {
                         // Evaluate NW sub (fine) cell SW node location.
                         for ( k = 1 ; k <= NUM_COMP_VECTOR2D; ++ k) {
                            l0 = l0 + 1;
                            if (l0 >= buffer_size_neighbour) return(1204);
                            Soln_Block_List.message_reschange_eastface_sendbuf[i_blk][0][l0] =
                                HALF*(Soln_ptr[i_blk].Grid.Node[i_min0][j].X[k]+
                                      Soln_ptr[i_blk].Grid.Node[i_min0][j+1].X[k])-x_ref[k];
                         } /* endfor */
	                 // Evaluate NE sub (fine) cell SW node location.
                         for ( k = 1 ; k <= NUM_COMP_VECTOR2D; ++ k) {
                            l0 = l0 + 1;
                            if (l0 >= buffer_size_neighbour) return(1205);
                            Soln_Block_List.message_reschange_eastface_sendbuf[i_blk][0][l0] =
                              Soln_ptr[i_blk].Grid.Cell[i_min0][j].Xc[k]-x_ref[k];
                         } /* endfor */
                      } /* endif */
                   } /* endfor */
                } else {
                   for ( j = j_min0 ; ((j_inc+1)/2) ? (j <= j_max0):(j >= j_max0) ; j += j_inc ) {
	              // Evaluate SE sub (fine) cell SW node location.
                      for ( k = 1 ; k <= NUM_COMP_VECTOR2D; ++ k) {
                         l0 = l0 + 1;
                         if (l0 >= buffer_size_neighbour) return(1206);
                         Soln_Block_List.message_reschange_eastface_sendbuf[i_blk][0][l0] =
                            HALF*(Soln_ptr[i_blk].Grid.Node[i_min0][j].X[k]+
                                  Soln_ptr[i_blk].Grid.Node[i_min0+1][j].X[k])-x_ref[k];
                      } /* endfor */
	              // Evaluate SW sub (fine) cell SW node location.
                      for ( k = 1 ; k <= NUM_COMP_VECTOR2D; ++ k) {
                         l0 = l0 + 1;
                         if (l0 >= buffer_size_neighbour) return(1207);
                         Soln_Block_List.message_reschange_eastface_sendbuf[i_blk][0][l0] =
                            Soln_ptr[i_blk].Grid.Node[i_min0][j].X[k]-x_ref[k];
                      } /* endfor */
                      if (j != j_max0) {
	                 // Evaluate NE sub (fine) cell SW node location.
                         for ( k = 1 ; k <= NUM_COMP_VECTOR2D; ++ k) {
                            l0 = l0 + 1;
                            if (l0 >= buffer_size_neighbour) return(1208);
                            Soln_Block_List.message_reschange_eastface_sendbuf[i_blk][0][l0] =
                              Soln_ptr[i_blk].Grid.Cell[i_min0][j].Xc[k]-x_ref[k];
                         } /* endfor */
                         // Evaluate NW sub (fine) cell SW node location.
                         for ( k = 1 ; k <= NUM_COMP_VECTOR2D; ++ k) {
                            l0 = l0 + 1;
                            if (l0 >= buffer_size_neighbour) return(1209);
                            Soln_Block_List.message_reschange_eastface_sendbuf[i_blk][0][l0] =
                                HALF*(Soln_ptr[i_blk].Grid.Node[i_min0][j].X[k]+
                                      Soln_ptr[i_blk].Grid.Node[i_min0][j+1].X[k])-x_ref[k];
                         } /* endfor */
                      } /* endif */
                   } /* endfor */
                } /* endif */
             } else {
                if (i_inc > 0) {
                   for ( j = j_min0 ; ((j_inc+1)/2) ? (j <= j_max0):(j >= j_max0) ; j += j_inc ) {
                      if (j != j_min0) {
                         // Evaluate NW sub (fine) cell SW node location.
                         for ( k = 1 ; k <= NUM_COMP_VECTOR2D; ++ k) {
                            l0 = l0 + 1;
                            if (l0 >= buffer_size_neighbour) return(1210);
                            Soln_Block_List.message_reschange_eastface_sendbuf[i_blk][0][l0] =
                                HALF*(Soln_ptr[i_blk].Grid.Node[i_min0][j].X[k]+
                                      Soln_ptr[i_blk].Grid.Node[i_min0][j+1].X[k])-x_ref[k];
                         } /* endfor */
	                 // Evaluate NE sub (fine) cell SW node location.
                         for ( k = 1 ; k <= NUM_COMP_VECTOR2D; ++ k) {
                            l0 = l0 + 1;
                            if (l0 >= buffer_size_neighbour) return(1211);
                            Soln_Block_List.message_reschange_eastface_sendbuf[i_blk][0][l0] =
                              Soln_ptr[i_blk].Grid.Cell[i_min0][j].Xc[k]-x_ref[k];
                         } /* endfor */
                      } /* endif */
	              // Evaluate SW sub (fine) cell SW node location.
                      for ( k = 1 ; k <= NUM_COMP_VECTOR2D; ++ k) {
                         l0 = l0 + 1;
                         if (l0 >= buffer_size_neighbour) return(1212);
                         Soln_Block_List.message_reschange_eastface_sendbuf[i_blk][0][l0] =
                            Soln_ptr[i_blk].Grid.Node[i_min0][j].X[k]-x_ref[k];
                      } /* endfor */
	              // Evaluate SE sub (fine) cell SW node location.
                      for ( k = 1 ; k <= NUM_COMP_VECTOR2D; ++ k) {
                         l0 = l0 + 1;
                         if (l0 >= buffer_size_neighbour) return(1213);
                         Soln_Block_List.message_reschange_eastface_sendbuf[i_blk][0][l0] =
                            HALF*(Soln_ptr[i_blk].Grid.Node[i_min0][j].X[k]+
                                  Soln_ptr[i_blk].Grid.Node[i_min0+1][j].X[k])-x_ref[k];
                      } /* endfor */
                   } /* endfor */
                } else {
                   for ( j = j_min0 ; ((j_inc+1)/2) ? (j <= j_max0):(j >= j_max0) ; j += j_inc ) {
                      if (j != j_min0) {
	                 // Evaluate NE sub (fine) cell SW node location.
                         for ( k = 1 ; k <= NUM_COMP_VECTOR2D; ++ k) {
                            l0 = l0 + 1;
                            if (l0 >= buffer_size_neighbour) return(1214);
                            Soln_Block_List.message_reschange_eastface_sendbuf[i_blk][0][l0] =
                              Soln_ptr[i_blk].Grid.Cell[i_min0][j].Xc[k]-x_ref[k];
                         } /* endfor */
                         // Evaluate NW sub (fine) cell SW node location.
                         for ( k = 1 ; k <= NUM_COMP_VECTOR2D; ++ k) {
                            l0 = l0 + 1;
                            if (l0 >= buffer_size_neighbour) return(1215);
                            Soln_Block_List.message_reschange_eastface_sendbuf[i_blk][0][l0] =
                                HALF*(Soln_ptr[i_blk].Grid.Node[i_min0][j].X[k]+
                                      Soln_ptr[i_blk].Grid.Node[i_min0][j+1].X[k])-x_ref[k];
                         } /* endfor */
                      } /* endif */
	              // Evaluate SE sub (fine) cell SW node location.
                      for ( k = 1 ; k <= NUM_COMP_VECTOR2D; ++ k) {
                         l0 = l0 + 1;
                         if (l0 >= buffer_size_neighbour) return(1216);
                         Soln_Block_List.message_reschange_eastface_sendbuf[i_blk][0][l0] =
                            HALF*(Soln_ptr[i_blk].Grid.Node[i_min0][j].X[k]+
                                  Soln_ptr[i_blk].Grid.Node[i_min0+1][j].X[k])-x_ref[k];
                      } /* endfor */
	              // Evaluate SW sub (fine) cell SW node location.
                      for ( k = 1 ; k <= NUM_COMP_VECTOR2D; ++ k) {
                         l0 = l0 + 1;
                         if (l0 >= buffer_size_neighbour) return(1217);
                         Soln_Block_List.message_reschange_eastface_sendbuf[i_blk][0][l0] =
                            Soln_ptr[i_blk].Grid.Node[i_min0][j].X[k]-x_ref[k];
                      } /* endfor */
                   } /* endfor */
                } /* endif */
             } /* endif */
             //
             // Second neighbour (mesh).
             //
             x_ref = Soln_ptr[i_blk].Grid.Node[i_ref1][j_ref1].X; // Reference node location.
             if (j_inc > 0) {
                if (i_inc > 0) {
                   for ( j = j_min1 ; ((j_inc+1)/2) ? (j <= j_max1):(j >= j_max1) ; j += j_inc ) {
	              // Evaluate SW sub (fine) cell SW node location.
                      for ( k = 1 ; k <= NUM_COMP_VECTOR2D; ++ k) {
                         l1 = l1 + 1;
                         if (l1 >= buffer_size_neighbour) return(1218);
                         Soln_Block_List.message_reschange_eastface_sendbuf[i_blk][1][l1] =
                            Soln_ptr[i_blk].Grid.Node[i_min1][j].X[k]-x_ref[k];
                      } /* endfor */
                      // Evaluate SE sub (fine) cell SW node location.
                      for ( k = 1 ; k <= NUM_COMP_VECTOR2D; ++ k) {
                         l1 = l1 + 1;
                         if (l1 >= buffer_size_neighbour) return(1219);
                         Soln_Block_List.message_reschange_eastface_sendbuf[i_blk][1][l1] =
                            HALF*(Soln_ptr[i_blk].Grid.Node[i_min1][j].X[k]+
                                  Soln_ptr[i_blk].Grid.Node[i_min1+1][j].X[k])-x_ref[k];
                      } /* endfor */
                      if (j != j_max1) {
                         // Evaluate NW sub (fine) cell SW node location.
                         for ( k = 1 ; k <= NUM_COMP_VECTOR2D; ++ k) {
                            l1 = l1 + 1;
                            if (l1 >= buffer_size_neighbour) return(1220);
                            Soln_Block_List.message_reschange_eastface_sendbuf[i_blk][1][l1] =
                               HALF*(Soln_ptr[i_blk].Grid.Node[i_min1][j].X[k]+
                                     Soln_ptr[i_blk].Grid.Node[i_min1][j+1].X[k])-x_ref[k];
                         } /* endfor */ 
                         // Evaluate NE sub (fine) cell SW node location.
                         for ( k = 1 ; k <= NUM_COMP_VECTOR2D; ++ k) {
                            l1 = l1 + 1;
                            if (l1 >= buffer_size_neighbour) return(1221);
                            Soln_Block_List.message_reschange_eastface_sendbuf[i_blk][1][l1] =
                               Soln_ptr[i_blk].Grid.Cell[i_min1][j].Xc[k]-x_ref[k];
                         } /* endfor */
                      } /* endif */
                   } /* endfor */
                } else {
                   for ( j = j_min1 ; ((j_inc+1)/2) ? (j <= j_max1):(j >= j_max1) ; j += j_inc ) {
                      // Evaluate SE sub (fine) cell SW node location.
                      for ( k = 1 ; k <= NUM_COMP_VECTOR2D; ++ k) {
                         l1 = l1 + 1;
                         if (l1 >= buffer_size_neighbour) return(1222);
                         Soln_Block_List.message_reschange_eastface_sendbuf[i_blk][1][l1] =
                            HALF*(Soln_ptr[i_blk].Grid.Node[i_min1][j].X[k]+
                                  Soln_ptr[i_blk].Grid.Node[i_min1+1][j].X[k])-x_ref[k];
                      } /* endfor */
	              // Evaluate SW sub (fine) cell SW node location.
                      for ( k = 1 ; k <= NUM_COMP_VECTOR2D; ++ k) {
                         l1 = l1 + 1;
                         if (l1 >= buffer_size_neighbour) return(1223);
                         Soln_Block_List.message_reschange_eastface_sendbuf[i_blk][1][l1] =
                            Soln_ptr[i_blk].Grid.Node[i_min1][j].X[k]-x_ref[k];
                      } /* endfor */
                      if (j != j_max1) {
                         // Evaluate NE sub (fine) cell SW node location.
                         for ( k = 1 ; k <= NUM_COMP_VECTOR2D; ++ k) {
                            l1 = l1 + 1;
                            if (l1 >= buffer_size_neighbour) return(1224);
                            Soln_Block_List.message_reschange_eastface_sendbuf[i_blk][1][l1] =
                               Soln_ptr[i_blk].Grid.Cell[i_min1][j].Xc[k]-x_ref[k];
                         } /* endfor */
                         // Evaluate NW sub (fine) cell SW node location.
                         for ( k = 1 ; k <= NUM_COMP_VECTOR2D; ++ k) {
                            l1 = l1 + 1;
                            if (l1 >= buffer_size_neighbour) return(1225);
                            Soln_Block_List.message_reschange_eastface_sendbuf[i_blk][1][l1] =
                               HALF*(Soln_ptr[i_blk].Grid.Node[i_min1][j].X[k]+
                                     Soln_ptr[i_blk].Grid.Node[i_min1][j+1].X[k])-x_ref[k];
                         } /* endfor */
                      } /* endif */
                   } /* endfor */
                } /* endif */
             } else {
                if (i_inc > 0) {
                   for ( j = j_min1 ; ((j_inc+1)/2) ? (j <= j_max1):(j >= j_max1) ; j += j_inc ) {
                      if (j != j_min1) {
                         // Evaluate NW sub (fine) cell SW node location.
                         for ( k = 1 ; k <= NUM_COMP_VECTOR2D; ++ k) {
                            l1 = l1 + 1;
                            if (l1 >= buffer_size_neighbour) return(1226);
                            Soln_Block_List.message_reschange_eastface_sendbuf[i_blk][1][l1] =
                               HALF*(Soln_ptr[i_blk].Grid.Node[i_min1][j].X[k]+
                                     Soln_ptr[i_blk].Grid.Node[i_min1][j+1].X[k])-x_ref[k];
                         } /* endfor */ 
                         // Evaluate NE sub (fine) cell SW node location.
                         for ( k = 1 ; k <= NUM_COMP_VECTOR2D; ++ k) {
                            l1 = l1 + 1;
                            if (l1 >= buffer_size_neighbour) return(1227);
                            Soln_Block_List.message_reschange_eastface_sendbuf[i_blk][1][l1] =
                               Soln_ptr[i_blk].Grid.Cell[i_min1][j].Xc[k]-x_ref[k];
                         } /* endfor */
                      } /* endif */
	              // Evaluate SW sub (fine) cell SW node location.
                      for ( k = 1 ; k <= NUM_COMP_VECTOR2D; ++ k) {
                         l1 = l1 + 1;
                         if (l1 >= buffer_size_neighbour) return(1228);
                         Soln_Block_List.message_reschange_eastface_sendbuf[i_blk][1][l1] =
                            Soln_ptr[i_blk].Grid.Node[i_min1][j].X[k]-x_ref[k];
                      } /* endfor */
                      // Evaluate SE sub (fine) cell SW node location.
                      for ( k = 1 ; k <= NUM_COMP_VECTOR2D; ++ k) {
                         l1 = l1 + 1;
                         if (l1 >= buffer_size_neighbour) return(1229);
                         Soln_Block_List.message_reschange_eastface_sendbuf[i_blk][1][l1] =
                            HALF*(Soln_ptr[i_blk].Grid.Node[i_min1][j].X[k]+
                                  Soln_ptr[i_blk].Grid.Node[i_min1+1][j].X[k])-x_ref[k];
                      } /* endfor */
                   } /* endfor */
                } else {
                   for ( j = j_min1 ; ((j_inc+1)/2) ? (j <= j_max1):(j >= j_max1) ; j += j_inc ) {
                      if (j != j_min1) {
                          // Evaluate NE sub (fine) cell SW node location.
                         for ( k = 1 ; k <= NUM_COMP_VECTOR2D; ++ k) {
                            l1 = l1 + 1;
                            if (l1 >= buffer_size_neighbour) return(1230);
                            Soln_Block_List.message_reschange_eastface_sendbuf[i_blk][1][l1] =
                               Soln_ptr[i_blk].Grid.Cell[i_min1][j].Xc[k]-x_ref[k];
                         } /* endfor */
                         // Evaluate NW sub (fine) cell SW node location.
                         for ( k = 1 ; k <= NUM_COMP_VECTOR2D; ++ k) {
                            l1 = l1 + 1;
                            if (l1 >= buffer_size_neighbour) return(1231);
                            Soln_Block_List.message_reschange_eastface_sendbuf[i_blk][1][l1] =
                               HALF*(Soln_ptr[i_blk].Grid.Node[i_min1][j].X[k]+
                                     Soln_ptr[i_blk].Grid.Node[i_min1][j+1].X[k])-x_ref[k];
                         } /* endfor */
                      } /* endif */
                      // Evaluate SE sub (fine) cell SW node location.
                      for ( k = 1 ; k <= NUM_COMP_VECTOR2D; ++ k) {
                         l1 = l1 + 1;
                         if (l1 >= buffer_size_neighbour) return(1232);
                         Soln_Block_List.message_reschange_eastface_sendbuf[i_blk][1][l1] =
                            HALF*(Soln_ptr[i_blk].Grid.Node[i_min1][j].X[k]+
                                  Soln_ptr[i_blk].Grid.Node[i_min1+1][j].X[k])-x_ref[k];
                      } /* endfor */
	              // Evaluate SW sub (fine) cell SW node location.
                      for ( k = 1 ; k <= NUM_COMP_VECTOR2D; ++ k) {
                         l1 = l1 + 1;
                         if (l1 >= buffer_size_neighbour) return(1233);
                         Soln_Block_List.message_reschange_eastface_sendbuf[i_blk][1][l1] =
                            Soln_ptr[i_blk].Grid.Node[i_min1][j].X[k]-x_ref[k];
                      } /* endfor */
                   } /* endfor */
                } /* endif */
             } /* endif */
             if (Soln_Block_List.Block[i_blk].infoE[0].dimen.i > 0 &&
                 Soln_Block_List.Block[i_blk].infoE[0].dimen.j > 0) {
	        i_min0 = Soln_ptr[i_blk].ICu;
	        i_max0 = Soln_ptr[i_blk].ICu;
	        i_min1 = Soln_ptr[i_blk].ICu;
	        i_max1 = Soln_ptr[i_blk].ICu;
	        i_inc = 1;
	        j_inc = 1;
             } else if (Soln_Block_List.Block[i_blk].infoE[0].dimen.j > 0) {
	        i_min0 = Soln_ptr[i_blk].ICu;
	        i_max0 = Soln_ptr[i_blk].ICu;
	        i_min1 = Soln_ptr[i_blk].ICu;
	        i_max1 = Soln_ptr[i_blk].ICu;
	        i_inc = -1;
	        j_inc = 1;
             } else if (Soln_Block_List.Block[i_blk].infoE[0].dimen.i > 0) {
	        i_min0 = Soln_ptr[i_blk].ICu;
	        i_max0 = Soln_ptr[i_blk].ICu;
	        i_min1 = Soln_ptr[i_blk].ICu;
	        i_max1 = Soln_ptr[i_blk].ICu;
	        i_inc = 1;
	        j_inc = -1;
             } else {
	        i_min0 = Soln_ptr[i_blk].ICu;
	        i_max0 = Soln_ptr[i_blk].ICu;
	        i_min1 = Soln_ptr[i_blk].ICu;
	        i_max1 = Soln_ptr[i_blk].ICu;
	        i_inc = -1;
	        j_inc = -1;
             } /* endif */
             if (j_inc > 0) {
                for ( k = 1 ; k <= 2; ++ k) {
                   for ( i  = i_min0 ; ((i_inc+1)/2) ? (i <= i_max0):(i >= i_max0) ; i += i_inc ) {
                      l0 = l0 + 2;
                      if (l0 >= buffer_size_neighbour) return(1234);
                      Soln_Block_List.message_reschange_eastface_sendbuf[i_blk][0][l0-1] = 
                         double(Soln_ptr[i_blk].Grid.BCtypeS[i]);
                      Soln_Block_List.message_reschange_eastface_sendbuf[i_blk][0][l0] = 
                         double(Soln_ptr[i_blk].Grid.BCtypeS[i]);
                   } /* endfor */
                   for ( i  = i_min1 ; ((i_inc+1)/2) ? (i <= i_max1):(i >= i_max1) ; i += i_inc ) {
                      l1 = l1 + 2;
                      if (l1 >= buffer_size_neighbour) return(1235);
                      Soln_Block_List.message_reschange_eastface_sendbuf[i_blk][1][l1-1] = 
                         double(Soln_ptr[i_blk].Grid.BCtypeN[i]);
                      Soln_Block_List.message_reschange_eastface_sendbuf[i_blk][1][l1] = 
                         double(Soln_ptr[i_blk].Grid.BCtypeN[i]);
                   } /* endfor */
                } /* endfor */
	     } else {
                for ( k = 1 ; k <= 2; ++ k) {
                   for ( i  = i_min0 ; ((i_inc+1)/2) ? (i <= i_max0):(i >= i_max0) ; i += i_inc ) {
                      l0 = l0 + 2;
                      if (l0 >= buffer_size_neighbour) return(1236);
                      Soln_Block_List.message_reschange_eastface_sendbuf[i_blk][0][l0-1] = 
                         double(Soln_ptr[i_blk].Grid.BCtypeN[i]);
                      Soln_Block_List.message_reschange_eastface_sendbuf[i_blk][0][l0] = 
                         double(Soln_ptr[i_blk].Grid.BCtypeN[i]);
                   } /* endfor */
                   for ( i  = i_min1 ; ((i_inc+1)/2) ? (i <= i_max1):(i >= i_max1) ; i += i_inc ) {
                      l1 = l1 + 1;
                      if (l1 >= buffer_size_neighbour) return(1237);
                      Soln_Block_List.message_reschange_eastface_sendbuf[i_blk][1][l1-1] = 
                         double(Soln_ptr[i_blk].Grid.BCtypeS[i]);
                      Soln_Block_List.message_reschange_eastface_sendbuf[i_blk][1][l1] = 
                         double(Soln_ptr[i_blk].Grid.BCtypeS[i]);
                   } /* endfor */
                } /* endfor */
             } /* endif */
          } /* endif */
       } /* endif */

       //////////////////////////////////////////////////////////////////
       // Load send buffers of solution blocks with refined West face //
       // neighbour at higher mesh resolution.                         //
       //////////////////////////////////////////////////////////////////
       if (Soln_Block_List.Block[i_blk].used && 
           (Soln_Block_List.Block[i_blk].nW == 2) &&
           (Soln_Block_List.Block[i_blk].info.level < 
            Soln_Block_List.Block[i_blk].infoW[0].level)) {
          if (!Send_Mesh_Geometry_Only) {
             buffer_size_neighbour = Soln_Block_List.Block[i_blk].infoW[0].dimen.ghost*
                                     (abs(Soln_Block_List.Block[i_blk].infoW[0].dimen.j)+
                                     Soln_Block_List.Block[i_blk].infoW[0].dimen.ghost)*
                                     Soln_ptr[i_blk].NumVar();
             if (Number_of_Solution_Variables > Soln_ptr[i_blk].NumVar()) 
                buffer_size_neighbour = buffer_size_neighbour +
                                        Soln_Block_List.Block[i_blk].infoW[0].dimen.ghost*
                                        (abs(Soln_Block_List.Block[i_blk].infoW[0].dimen.j)+
                                        Soln_Block_List.Block[i_blk].infoW[0].dimen.ghost+1)*
                                        NUM_COMP_VECTOR2D+
                                        2*Soln_Block_List.Block[i_blk].infoW[0].dimen.ghost;
          } else {
             buffer_size_neighbour = Soln_Block_List.Block[i_blk].infoW[0].dimen.ghost*
                                     (abs(Soln_Block_List.Block[i_blk].infoW[0].dimen.j)+
                                     Soln_Block_List.Block[i_blk].infoW[0].dimen.ghost+1)*
                                     NUM_COMP_VECTOR2D+
                                     2*Soln_Block_List.Block[i_blk].infoW[0].dimen.ghost;
          } /* endif */
          l0 = -1;
          l1 = -1;
          // Load ghost cell solution variable information as required.
          if (!Send_Mesh_Geometry_Only) {
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
             i = Soln_ptr[i_blk].LoadSendBuffer_C2F(Soln_Block_List.message_reschange_westface_sendbuf[i_blk][0],
                                                    l0,buffer_size_neighbour,
                                                    i_min0,i_max0,i_inc,
                                                    j_min0,j_max0,j_inc);
	     if (i != 0) return(1100);
             //
             // Second neighbour (solution).
             //
             i = Soln_ptr[i_blk].LoadSendBuffer_C2F(Soln_Block_List.message_reschange_westface_sendbuf[i_blk][1],
                                                    l1,buffer_size_neighbour,
                                                    i_min1,i_max1,i_inc,
                                                    j_min1,j_max1,j_inc);
	     if (i != 0) return(1101);
          } /* endif */
          // Load ghost cell mesh information as required.
          if (Send_Mesh_Geometry_Only ||
              Number_of_Solution_Variables > Soln_ptr[i_blk].NumVar()) {
             if (Soln_Block_List.Block[i_blk].infoW[0].dimen.i > 0 &&
                 Soln_Block_List.Block[i_blk].infoW[0].dimen.j > 0) {
	        i_min0 = Soln_ptr[i_blk].Grid.INl;
	        i_max0 = Soln_ptr[i_blk].Grid.INl;
	        i_min1 = Soln_ptr[i_blk].Grid.INl;
	        i_max1 = Soln_ptr[i_blk].Grid.INl;
	        i_inc = 1;
	        j_min0 = Soln_ptr[i_blk].Grid.JNl;
	        j_max0 = Soln_ptr[i_blk].Grid.JNl+(Soln_ptr[i_blk].Grid.JNu-Soln_ptr[i_blk].Grid.JNl)/2+1;
	        j_min1 = Soln_ptr[i_blk].Grid.JNu-((Soln_ptr[i_blk].Grid.JNu-Soln_ptr[i_blk].Grid.JNl)/2+1);
	        j_max1 = Soln_ptr[i_blk].Grid.JNu;
	        j_inc = 1;
                i_ref0 = i_min1;
                j_ref0 = j_min0;
                i_ref1 = i_min1;
                j_ref1 = j_min1+1;
             } else if (Soln_Block_List.Block[i_blk].infoW[0].dimen.j > 0) {
	        i_min0 = Soln_ptr[i_blk].Grid.INl;
	        i_max0 = Soln_ptr[i_blk].Grid.INl;
	        i_min1 = Soln_ptr[i_blk].Grid.INl;
	        i_max1 = Soln_ptr[i_blk].Grid.INl;
	        i_inc = -1;
	        j_min0 = Soln_ptr[i_blk].Grid.JNl;
	        j_max0 = Soln_ptr[i_blk].Grid.JNl+(Soln_ptr[i_blk].Grid.JNu-Soln_ptr[i_blk].Grid.JNl)/2+1;
	        j_min1 = Soln_ptr[i_blk].Grid.JNu-((Soln_ptr[i_blk].Grid.JNu-Soln_ptr[i_blk].Grid.JNl)/2+1);
	        j_max1 = Soln_ptr[i_blk].Grid.JNu;
	        j_inc = 1;
                i_ref0 = i_min1;
                j_ref0 = j_min0;
                i_ref1 = i_min1;
                j_ref1 = j_min1+1;
             } else if (Soln_Block_List.Block[i_blk].infoW[0].dimen.i > 0) {
	        i_min0 = Soln_ptr[i_blk].Grid.INl;
	        i_max0 = Soln_ptr[i_blk].Grid.INl;
	        i_min1 = Soln_ptr[i_blk].Grid.INl;
	        i_max1 = Soln_ptr[i_blk].Grid.INl;
	        i_inc = 1;
	        j_min0 = Soln_ptr[i_blk].Grid.JNl+(Soln_ptr[i_blk].Grid.JNu-Soln_ptr[i_blk].Grid.JNl)/2+1;
	        j_max0 = Soln_ptr[i_blk].Grid.JNl;
	        j_min1 = Soln_ptr[i_blk].Grid.JNu;
	        j_max1 = Soln_ptr[i_blk].Grid.JNu-((Soln_ptr[i_blk].Grid.JNu-Soln_ptr[i_blk].Grid.JNl)/2+1);
	        j_inc = -1;
                i_ref0 = i_min1;
                j_ref0 = j_min0-1;
                i_ref1 = i_min1;
                j_ref1 = j_min1;
             } else {
	        i_min0 = Soln_ptr[i_blk].Grid.INl;
	        i_max0 = Soln_ptr[i_blk].Grid.INl;
	        i_min1 = Soln_ptr[i_blk].Grid.INl;
	        i_max1 = Soln_ptr[i_blk].Grid.INl;
	        i_inc = -1;
	        j_min0 = Soln_ptr[i_blk].Grid.JNl+(Soln_ptr[i_blk].Grid.JNu-Soln_ptr[i_blk].Grid.JNl)/2+1;
	        j_max0 = Soln_ptr[i_blk].Grid.JNl;
	        j_min1 = Soln_ptr[i_blk].Grid.JNu;
	        j_max1 = Soln_ptr[i_blk].Grid.JNu-((Soln_ptr[i_blk].Grid.JNu-Soln_ptr[i_blk].Grid.JNl)/2+1);
	        j_inc = -1;
                i_ref0 = i_min1;
                j_ref0 = j_min0-1;
                i_ref1 = i_min1;
                j_ref1 = j_min1;
             } /* endif */
             //
             // First neighbour (mesh).
             //
             x_ref = Soln_ptr[i_blk].Grid.Node[i_ref0][j_ref0].X; // Reference node location.
             // Four different orderings to consider depending on the value of i_inc & j_inc.
             if (j_inc > 0) {
                if (i_inc > 0) {
                   for ( j = j_min0 ; ((j_inc+1)/2) ? (j <= j_max0):(j >= j_max0) ; j += j_inc ) {
	              // Evaluate SW sub (fine) cell SE node location.
                      for ( k = 1 ; k <= NUM_COMP_VECTOR2D; ++ k) {
                         l0 = l0 + 1;
                         if (l0 >= buffer_size_neighbour) return(1102);
                         Soln_Block_List.message_reschange_westface_sendbuf[i_blk][0][l0] =
                            HALF*(Soln_ptr[i_blk].Grid.Node[i_min0][j].X[k]+
                                  Soln_ptr[i_blk].Grid.Node[i_min0+1][j].X[k])-x_ref[k];
                      } /* endfor */
	              // Evaluate SE sub (fine) cell SE node location.
                      for ( k = 1 ; k <= NUM_COMP_VECTOR2D; ++ k) {
                         l0 = l0 + 1;
                         if (l0 >= buffer_size_neighbour) return(1103);
                         Soln_Block_List.message_reschange_westface_sendbuf[i_blk][0][l0] =
                            Soln_ptr[i_blk].Grid.Node[i_min0+1][j].X[k]-x_ref[k];
                      } /* endfor */
                      if (j != j_max0) {
                         // Evaluate NW sub (fine) cell SE node location.
                         for ( k = 1 ; k <= NUM_COMP_VECTOR2D; ++ k) {
                            l0 = l0 + 1;
                            if (l0 >= buffer_size_neighbour) return(1104);
                            Soln_Block_List.message_reschange_westface_sendbuf[i_blk][0][l0] =
                               Soln_ptr[i_blk].Grid.Cell[i_min0][j].Xc[k]-x_ref[k];
                         } /* endfor */
	                 // Evaluate NE sub (fine) cell SE node location.
                         for ( k = 1 ; k <= NUM_COMP_VECTOR2D; ++ k) {
                            l0 = l0 + 1;
                            if (l0 >= buffer_size_neighbour) return(1105);
                            Soln_Block_List.message_reschange_westface_sendbuf[i_blk][0][l0] =
                                HALF*(Soln_ptr[i_blk].Grid.Node[i_min0+1][j].X[k]+
                                      Soln_ptr[i_blk].Grid.Node[i_min0+1][j+1].X[k])-x_ref[k];
                         } /* endfor */
                      } /* endif */
                   } /* endfor */
                } else {
                   for ( j = j_min0 ; ((j_inc+1)/2) ? (j <= j_max0):(j >= j_max0) ; j += j_inc ) {
	              // Evaluate SE sub (fine) cell SE node location.
                      for ( k = 1 ; k <= NUM_COMP_VECTOR2D; ++ k) {
                         l0 = l0 + 1;
                         if (l0 >= buffer_size_neighbour) return(1106);
                         Soln_Block_List.message_reschange_westface_sendbuf[i_blk][0][l0] =
                            Soln_ptr[i_blk].Grid.Node[i_min0+1][j].X[k]-x_ref[k];
                      } /* endfor */
	              // Evaluate SW sub (fine) cell SE node location.
                      for ( k = 1 ; k <= NUM_COMP_VECTOR2D; ++ k) {
                         l0 = l0 + 1;
                         if (l0 >= buffer_size_neighbour) return(1107);
                         Soln_Block_List.message_reschange_westface_sendbuf[i_blk][0][l0] =
                            HALF*(Soln_ptr[i_blk].Grid.Node[i_min0][j].X[k]+
                                  Soln_ptr[i_blk].Grid.Node[i_min0+1][j].X[k])-x_ref[k];
                      } /* endfor */
                      if (j != j_max0) {
	                 // Evaluate NE sub (fine) cell SE node location.
                         for ( k = 1 ; k <= NUM_COMP_VECTOR2D; ++ k) {
                            l0 = l0 + 1;
                            if (l0 >= buffer_size_neighbour) return(1108);
                            Soln_Block_List.message_reschange_westface_sendbuf[i_blk][0][l0] =
                                HALF*(Soln_ptr[i_blk].Grid.Node[i_min0+1][j].X[k]+
                                      Soln_ptr[i_blk].Grid.Node[i_min0+1][j+1].X[k])-x_ref[k];
                         } /* endfor */
                         // Evaluate NW sub (fine) cell SE node location.
                         for ( k = 1 ; k <= NUM_COMP_VECTOR2D; ++ k) {
                            l0 = l0 + 1;
                            if (l0 >= buffer_size_neighbour) return(1109);
                            Soln_Block_List.message_reschange_westface_sendbuf[i_blk][0][l0] =
                               Soln_ptr[i_blk].Grid.Cell[i_min0][j].Xc[k]-x_ref[k];
                         } /* endfor */
                      } /* endif */
                   } /* endfor */
                } /* endif */
             } else {
                if (i_inc > 0) {
                   for ( j = j_min0 ; ((j_inc+1)/2) ? (j <= j_max0):(j >= j_max0) ; j += j_inc ) {
                      if (j != j_min0) {
                         // Evaluate NW sub (fine) cell SE node location.
                         for ( k = 1 ; k <= NUM_COMP_VECTOR2D; ++ k) {
                            l0 = l0 + 1;
                            if (l0 >= buffer_size_neighbour) return(1110);
                            Soln_Block_List.message_reschange_westface_sendbuf[i_blk][0][l0] =
                               Soln_ptr[i_blk].Grid.Cell[i_min0][j].Xc[k]-x_ref[k];
                         } /* endfor */
	                 // Evaluate NE sub (fine) cell SE node location.
                         for ( k = 1 ; k <= NUM_COMP_VECTOR2D; ++ k) {
                            l0 = l0 + 1;
                            if (l0 >= buffer_size_neighbour) return(1111);
                            Soln_Block_List.message_reschange_westface_sendbuf[i_blk][0][l0] =
                                HALF*(Soln_ptr[i_blk].Grid.Node[i_min0+1][j].X[k]+
                                      Soln_ptr[i_blk].Grid.Node[i_min0+1][j+1].X[k])-x_ref[k];
                         } /* endfor */
                      } /* endif */
	              // Evaluate SW sub (fine) cell SE node location.
                      for ( k = 1 ; k <= NUM_COMP_VECTOR2D; ++ k) {
                         l0 = l0 + 1;
                         if (l0 >= buffer_size_neighbour) return(1112);
                         Soln_Block_List.message_reschange_westface_sendbuf[i_blk][0][l0] =
                            HALF*(Soln_ptr[i_blk].Grid.Node[i_min0][j].X[k]+
                                  Soln_ptr[i_blk].Grid.Node[i_min0+1][j].X[k])-x_ref[k];
                      } /* endfor */
	              // Evaluate SE sub (fine) cell SE node location.
                      for ( k = 1 ; k <= NUM_COMP_VECTOR2D; ++ k) {
                         l0 = l0 + 1;
                         if (l0 >= buffer_size_neighbour) return(1113);
                         Soln_Block_List.message_reschange_westface_sendbuf[i_blk][0][l0] =
                            Soln_ptr[i_blk].Grid.Node[i_min0+1][j].X[k]-x_ref[k];
                      } /* endfor */
                   } /* endfor */
                } else {
                   for ( j = j_min0 ; ((j_inc+1)/2) ? (j <= j_max0):(j >= j_max0) ; j += j_inc ) {
                      if (j != j_min0) {
	                 // Evaluate NE sub (fine) cell SE node location.
                         for ( k = 1 ; k <= NUM_COMP_VECTOR2D; ++ k) {
                            l0 = l0 + 1;
                            if (l0 >= buffer_size_neighbour) return(1114);
                            Soln_Block_List.message_reschange_westface_sendbuf[i_blk][0][l0] =
                                HALF*(Soln_ptr[i_blk].Grid.Node[i_min0+1][j].X[k]+
                                      Soln_ptr[i_blk].Grid.Node[i_min0+1][j+1].X[k])-x_ref[k];
                         } /* endfor */
                         // Evaluate NW sub (fine) cell SE node location.
                         for ( k = 1 ; k <= NUM_COMP_VECTOR2D; ++ k) {
                            l0 = l0 + 1;
                            if (l0 >= buffer_size_neighbour) return(1115);
                            Soln_Block_List.message_reschange_westface_sendbuf[i_blk][0][l0] =
                               Soln_ptr[i_blk].Grid.Cell[i_min0][j].Xc[k]-x_ref[k];
                         } /* endfor */
                      } /* endif */
	              // Evaluate SE sub (fine) cell SE node location.
                      for ( k = 1 ; k <= NUM_COMP_VECTOR2D; ++ k) {
                         l0 = l0 + 1;
                         if (l0 >= buffer_size_neighbour) return(1116);
                         Soln_Block_List.message_reschange_westface_sendbuf[i_blk][0][l0] =
                            Soln_ptr[i_blk].Grid.Node[i_min0+1][j].X[k]-x_ref[k];
                      } /* endfor */
	              // Evaluate SW sub (fine) cell SE node location.
                      for ( k = 1 ; k <= NUM_COMP_VECTOR2D; ++ k) {
                         l0 = l0 + 1;
                         if (l0 >= buffer_size_neighbour) return(1117);
                         Soln_Block_List.message_reschange_westface_sendbuf[i_blk][0][l0] =
                            HALF*(Soln_ptr[i_blk].Grid.Node[i_min0][j].X[k]+
                                  Soln_ptr[i_blk].Grid.Node[i_min0+1][j].X[k])-x_ref[k];
                      } /* endfor */
                   } /* endfor */
                } /* endif */
             } /* endif */
             //
             // Second neighbour (mesh).
             //
             x_ref = Soln_ptr[i_blk].Grid.Node[i_ref1][j_ref1].X; // Reference node location.
             if (j_inc > 0) {
                if (i_inc > 0) {
                   for ( j = j_min1 ; ((j_inc+1)/2) ? (j <= j_max1):(j >= j_max1) ; j += j_inc ) {
	              // Evaluate SW sub (fine) cell SE node location.
                      for ( k = 1 ; k <= NUM_COMP_VECTOR2D; ++ k) {
                         l1 = l1 + 1;
                         if (l1 >= buffer_size_neighbour) return(1118);
                         Soln_Block_List.message_reschange_westface_sendbuf[i_blk][1][l1] =
                            HALF*(Soln_ptr[i_blk].Grid.Node[i_min1][j].X[k]+
                                  Soln_ptr[i_blk].Grid.Node[i_min1+1][j].X[k])-x_ref[k];
                      } /* endfor */
                      // Evaluate SE sub (fine) cell SE node location.
                      for ( k = 1 ; k <= NUM_COMP_VECTOR2D; ++ k) {
                         l1 = l1 + 1;
                         if (l1 >= buffer_size_neighbour) return(1119);
                         Soln_Block_List.message_reschange_westface_sendbuf[i_blk][1][l1] =
                            Soln_ptr[i_blk].Grid.Node[i_min1+1][j].X[k]-x_ref[k];
                      } /* endfor */
                      if (j != j_max1) {
                         // Evaluate NW sub (fine) cell SE node location.
                         for ( k = 1 ; k <= NUM_COMP_VECTOR2D; ++ k) {
                            l1 = l1 + 1;
                            if (l1 >= buffer_size_neighbour) return(1120);
                            Soln_Block_List.message_reschange_westface_sendbuf[i_blk][1][l1] =
                               Soln_ptr[i_blk].Grid.Cell[i_min1][j].Xc[k]-x_ref[k];
                         } /* endfor */ 
                         // Evaluate NE sub (fine) cell SE node location.
                         for ( k = 1 ; k <= NUM_COMP_VECTOR2D; ++ k) {
                            l1 = l1 + 1;
                            if (l1 >= buffer_size_neighbour) return(1121);
                            Soln_Block_List.message_reschange_westface_sendbuf[i_blk][1][l1] =
                               HALF*(Soln_ptr[i_blk].Grid.Node[i_min1+1][j].X[k]+
                                     Soln_ptr[i_blk].Grid.Node[i_min1+1][j+1].X[k])-x_ref[k];
                         } /* endfor */
                      } /* endif */
                   } /* endfor */
                } else {
                   for ( j = j_min1 ; ((j_inc+1)/2) ? (j <= j_max1):(j >= j_max1) ; j += j_inc ) {
                      // Evaluate SE sub (fine) cell SE node location.
                      for ( k = 1 ; k <= NUM_COMP_VECTOR2D; ++ k) {
                         l1 = l1 + 1;
                         if (l1 >= buffer_size_neighbour) return(1122);
                         Soln_Block_List.message_reschange_westface_sendbuf[i_blk][1][l1] =
                            Soln_ptr[i_blk].Grid.Node[i_min1+1][j].X[k]-x_ref[k];
                      } /* endfor */
	              // Evaluate SW sub (fine) cell SE node location.
                      for ( k = 1 ; k <= NUM_COMP_VECTOR2D; ++ k) {
                         l1 = l1 + 1;
                         if (l1 >= buffer_size_neighbour) return(1123);
                         Soln_Block_List.message_reschange_westface_sendbuf[i_blk][1][l1] =
                            HALF*(Soln_ptr[i_blk].Grid.Node[i_min1][j].X[k]+
                                  Soln_ptr[i_blk].Grid.Node[i_min1+1][j].X[k])-x_ref[k];
                      } /* endfor */
                      if (j != j_max1) {
                          // Evaluate NE sub (fine) cell SE node location.
                         for ( k = 1 ; k <= NUM_COMP_VECTOR2D; ++ k) {
                            l1 = l1 + 1;
                            if (l1 >= buffer_size_neighbour) return(1124);
                            Soln_Block_List.message_reschange_westface_sendbuf[i_blk][1][l1] =
                               HALF*(Soln_ptr[i_blk].Grid.Node[i_min1+1][j].X[k]+
                                     Soln_ptr[i_blk].Grid.Node[i_min1+1][j+1].X[k])-x_ref[k];
                         } /* endfor */
                         // Evaluate NW sub (fine) cell SE node location.
                         for ( k = 1 ; k <= NUM_COMP_VECTOR2D; ++ k) {
                            l1 = l1 + 1;
                            if (l1 >= buffer_size_neighbour) return(1125);
                            Soln_Block_List.message_reschange_westface_sendbuf[i_blk][1][l1] =
                               Soln_ptr[i_blk].Grid.Cell[i_min1][j].Xc[k]-x_ref[k];
                         } /* endfor */
                      } /* endif */
                   } /* endfor */
                } /* endif */
             } else {
                if (i_inc > 0) {
                   for ( j = j_min1 ; ((j_inc+1)/2) ? (j <= j_max1):(j >= j_max1) ; j += j_inc ) {
                      if (j != j_min1) {
                         // Evaluate NW sub (fine) cell SE node location.
                         for ( k = 1 ; k <= NUM_COMP_VECTOR2D; ++ k) {
                            l1 = l1 + 1;
                            if (l1 >= buffer_size_neighbour) return(1126);
                            Soln_Block_List.message_reschange_westface_sendbuf[i_blk][1][l1] =
                               Soln_ptr[i_blk].Grid.Cell[i_min1][j].Xc[k]-x_ref[k];
                         } /* endfor */ 
                         // Evaluate NE sub (fine) cell SE node location.
                         for ( k = 1 ; k <= NUM_COMP_VECTOR2D; ++ k) {
                            l1 = l1 + 1;
                            if (l1 >= buffer_size_neighbour) return(1127);
                            Soln_Block_List.message_reschange_westface_sendbuf[i_blk][1][l1] =
                               HALF*(Soln_ptr[i_blk].Grid.Node[i_min1+1][j].X[k]+
                                     Soln_ptr[i_blk].Grid.Node[i_min1+1][j+1].X[k])-x_ref[k];
                         } /* endfor */
                      } /* endif */
	              // Evaluate SW sub (fine) cell SE node location.
                      for ( k = 1 ; k <= NUM_COMP_VECTOR2D; ++ k) {
                         l1 = l1 + 1;
                         if (l1 >= buffer_size_neighbour) return(1128);
                         Soln_Block_List.message_reschange_westface_sendbuf[i_blk][1][l1] =
                            HALF*(Soln_ptr[i_blk].Grid.Node[i_min1][j].X[k]+
                                  Soln_ptr[i_blk].Grid.Node[i_min1+1][j].X[k])-x_ref[k];
                      } /* endfor */
                      // Evaluate SE sub (fine) cell SE node location.
                      for ( k = 1 ; k <= NUM_COMP_VECTOR2D; ++ k) {
                         l1 = l1 + 1;
                         if (l1 >= buffer_size_neighbour) return(1129);
                         Soln_Block_List.message_reschange_westface_sendbuf[i_blk][1][l1] =
                            Soln_ptr[i_blk].Grid.Node[i_min1+1][j].X[k]-x_ref[k];
                      } /* endfor */
                   } /* endfor */
                } else {
                   for ( j = j_min1 ; ((j_inc+1)/2) ? (j <= j_max1):(j >= j_max1) ; j += j_inc ) {
                      if (j != j_min1) {
                          // Evaluate NE sub (fine) cell SE node location.
                         for ( k = 1 ; k <= NUM_COMP_VECTOR2D; ++ k) {
                            l1 = l1 + 1;
                            if (l1 >= buffer_size_neighbour) return(1130);
                            Soln_Block_List.message_reschange_westface_sendbuf[i_blk][1][l1] =
                               HALF*(Soln_ptr[i_blk].Grid.Node[i_min1+1][j].X[k]+
                                     Soln_ptr[i_blk].Grid.Node[i_min1+1][j+1].X[k])-x_ref[k];
                         } /* endfor */
                         // Evaluate NW sub (fine) cell SE node location.
                         for ( k = 1 ; k <= NUM_COMP_VECTOR2D; ++ k) {
                            l1 = l1 + 1;
                            if (l1 >= buffer_size_neighbour) return(1131);
                            Soln_Block_List.message_reschange_westface_sendbuf[i_blk][1][l1] =
                               Soln_ptr[i_blk].Grid.Cell[i_min1][j].Xc[k]-x_ref[k];
                         } /* endfor */
                      } /* endif */
                      // Evaluate SE sub (fine) cell SE node location.
                      for ( k = 1 ; k <= NUM_COMP_VECTOR2D; ++ k) {
                         l1 = l1 + 1;
                         if (l1 >= buffer_size_neighbour) return(1132);
                         Soln_Block_List.message_reschange_westface_sendbuf[i_blk][1][l1] =
                            Soln_ptr[i_blk].Grid.Node[i_min1+1][j].X[k]-x_ref[k];
                      } /* endfor */
	              // Evaluate SW sub (fine) cell SE node location.
                      for ( k = 1 ; k <= NUM_COMP_VECTOR2D; ++ k) {
                         l1 = l1 + 1;
                         if (l1 >= buffer_size_neighbour) return(1133);
                         Soln_Block_List.message_reschange_westface_sendbuf[i_blk][1][l1] =
                            HALF*(Soln_ptr[i_blk].Grid.Node[i_min1][j].X[k]+
                                  Soln_ptr[i_blk].Grid.Node[i_min1+1][j].X[k])-x_ref[k];
                      } /* endfor */
                   } /* endfor */
                } /* endif */
             } /* endif */
             if (Soln_Block_List.Block[i_blk].infoW[0].dimen.i > 0 &&
                 Soln_Block_List.Block[i_blk].infoW[0].dimen.j > 0) {
	        i_min0 = Soln_ptr[i_blk].ICl;
	        i_max0 = Soln_ptr[i_blk].ICl;
	        i_min1 = Soln_ptr[i_blk].ICl;
	        i_max1 = Soln_ptr[i_blk].ICl;
	        i_inc = 1;
	        j_inc = 1;
             } else if (Soln_Block_List.Block[i_blk].infoW[0].dimen.j > 0) {
	        i_min0 = Soln_ptr[i_blk].ICl;
	        i_max0 = Soln_ptr[i_blk].ICl;
	        i_min1 = Soln_ptr[i_blk].ICl;
	        i_max1 = Soln_ptr[i_blk].ICl;
	        i_inc = -1;
	        j_inc = 1;
             } else if (Soln_Block_List.Block[i_blk].infoW[0].dimen.i > 0) {
	        i_min0 = Soln_ptr[i_blk].ICl;
	        i_max0 = Soln_ptr[i_blk].ICl;
	        i_min1 = Soln_ptr[i_blk].ICl;
	        i_max1 = Soln_ptr[i_blk].ICl;
	        i_inc = 1;
	        j_inc = -1;
             } else {
	        i_min0 = Soln_ptr[i_blk].ICl;
	        i_max0 = Soln_ptr[i_blk].ICl;
	        i_min1 = Soln_ptr[i_blk].ICl;
	        i_max1 = Soln_ptr[i_blk].ICl;
	        i_inc = -1;
	        j_inc = -1;
             } /* endif */
             if (j_inc > 0) {
                for ( k = 1 ; k <= 2; ++ k) {
                   for ( i  = i_min0 ; ((i_inc+1)/2) ? (i <= i_max0):(i >= i_max0) ; i += i_inc ) {
                      l0 = l0 + 2;
                      if (l0 >= buffer_size_neighbour) return(1134);
                      Soln_Block_List.message_reschange_westface_sendbuf[i_blk][0][l0-1] = 
                         double(Soln_ptr[i_blk].Grid.BCtypeS[i]);
                      Soln_Block_List.message_reschange_westface_sendbuf[i_blk][0][l0] = 
                         double(Soln_ptr[i_blk].Grid.BCtypeS[i]);
                   } /* endfor */
                   for ( i  = i_min1 ; ((i_inc+1)/2) ? (i <= i_max1):(i >= i_max1) ; i += i_inc ) {
                      l1 = l1 + 2;
                      if (l1 >= buffer_size_neighbour) return(1135);
                      Soln_Block_List.message_reschange_westface_sendbuf[i_blk][1][l1-1] = 
                         double(Soln_ptr[i_blk].Grid.BCtypeN[i]);
                      Soln_Block_List.message_reschange_westface_sendbuf[i_blk][1][l1] = 
                         double(Soln_ptr[i_blk].Grid.BCtypeN[i]);
                   } /* endfor */
                } /* endfor */
	     } else {
                for ( k = 1 ; k <= 2; ++ k) {
                   for ( i  = i_min0 ; ((i_inc+1)/2) ? (i <= i_max0):(i >= i_max0) ; i += i_inc ) {
                      l0 = l0 + 2;
                      if (l0 >= buffer_size_neighbour) return(1136);
                      Soln_Block_List.message_reschange_westface_sendbuf[i_blk][0][l0-1] = 
                         double(Soln_ptr[i_blk].Grid.BCtypeN[i]);
                      Soln_Block_List.message_reschange_westface_sendbuf[i_blk][0][l0] = 
                         double(Soln_ptr[i_blk].Grid.BCtypeN[i]);
                   } /* endfor */
                   for ( i  = i_min1 ; ((i_inc+1)/2) ? (i <= i_max1):(i >= i_max1) ; i += i_inc ) {
                      l1 = l1 + 1;
                      if (l1 >= buffer_size_neighbour) return(1137);
                      Soln_Block_List.message_reschange_westface_sendbuf[i_blk][1][l1-1] = 
                         double(Soln_ptr[i_blk].Grid.BCtypeS[i]);
                      Soln_Block_List.message_reschange_westface_sendbuf[i_blk][1][l1] = 
                         double(Soln_ptr[i_blk].Grid.BCtypeS[i]);
                   } /* endfor */
                } /* endfor */
             } /* endif */
          } /* endif */
       } /* endif */

       //////////////////////////////////////////////////////////////////
       // Load send buffers of solution blocks with refined North West //
       // corner neighbour at higher mesh resolution.                  //
       //////////////////////////////////////////////////////////////////
       if (Soln_Block_List.Block[i_blk].used && 
           (Soln_Block_List.Block[i_blk].nNW == 1) &&
           (Soln_Block_List.Block[i_blk].info.level < 
            Soln_Block_List.Block[i_blk].infoNW[0].level)) {
          buffer_size_neighbour = sqr(Soln_Block_List.Block[i_blk].infoNW[0].dimen.ghost)*
                                  Number_of_Solution_Variables;
          l0 = -1;
          // Load ghost cell solution variable information as required.
          if (!Send_Mesh_Geometry_Only) {
             if (Soln_Block_List.Block[i_blk].infoNW[0].dimen.i > 0 &&
                 Soln_Block_List.Block[i_blk].infoNW[0].dimen.j > 0) {
	        i_min0 = Soln_ptr[i_blk].ICl;
                i_max0 = Soln_ptr[i_blk].ICl;
	        i_inc = 1;
	        j_min0 = Soln_ptr[i_blk].JCu;
	        j_max0 = Soln_ptr[i_blk].JCu;
	        j_inc = 1;
             } else if (Soln_Block_List.Block[i_blk].infoNW[0].dimen.j > 0) {
	        i_min0 = Soln_ptr[i_blk].ICl;
	        i_max0 = Soln_ptr[i_blk].ICl;
	        i_inc = -1;
	        j_min0 = Soln_ptr[i_blk].JCu;
	        j_max0 = Soln_ptr[i_blk].JCu;
	        j_inc = 1;
             } else if (Soln_Block_List.Block[i_blk].infoNW[0].dimen.i > 0) {
	        i_min0 = Soln_ptr[i_blk].ICl;
	        i_max0 = Soln_ptr[i_blk].ICl;
	        i_inc = 1;
	        j_min0 = Soln_ptr[i_blk].JCu;
	        j_max0 = Soln_ptr[i_blk].JCu;
	        j_inc = -1;
             } else {
	        i_min0 = Soln_ptr[i_blk].ICl;
	        i_max0 = Soln_ptr[i_blk].ICl;
	        i_inc = -1;
	        j_min0 = Soln_ptr[i_blk].JCu;
	        j_max0 = Soln_ptr[i_blk].JCu;
	        j_inc = -1;
             } /* endif */
             i = Soln_ptr[i_blk].LoadSendBuffer_C2F(Soln_Block_List.message_reschange_northwestcorner_sendbuf[i_blk],
                                                    l0,buffer_size_neighbour,
                                                    i_min0,i_max0,i_inc,
                                                    j_min0,j_max0,j_inc);
	     if (i != 0) return(4200);
          } /* endif */
          // Load ghost cell mesh information as required.
          if (Send_Mesh_Geometry_Only ||
              Number_of_Solution_Variables > Soln_ptr[i_blk].NumVar()) {
             if (Soln_Block_List.Block[i_blk].infoNW[0].dimen.i > 0 &&
                 Soln_Block_List.Block[i_blk].infoNW[0].dimen.j > 0) {
	        i_min0 = Soln_ptr[i_blk].Grid.INl;
	        i_max0 = Soln_ptr[i_blk].Grid.INl;
	        i_inc = 1;
	        j_min0 = Soln_ptr[i_blk].Grid.JNu-1;
	        j_max0 = Soln_ptr[i_blk].Grid.JNu-1;
	        j_inc = 1;
                i_ref0 = i_min0;
                j_ref0 = j_min0+1;
             } else if (Soln_Block_List.Block[i_blk].infoNW[0].dimen.j > 0) {
	        i_min0 = Soln_ptr[i_blk].Grid.INl;
	        i_max0 = Soln_ptr[i_blk].Grid.INl;
	        i_inc = -1;
	        j_min0 = Soln_ptr[i_blk].Grid.JNu-1;
	        j_max0 = Soln_ptr[i_blk].Grid.JNu-1;
	        j_inc = 1;
                i_ref0 = i_min0;
                j_ref0 = j_min0+1;
             } else if (Soln_Block_List.Block[i_blk].infoNW[0].dimen.i > 0) {
	        i_min0 = Soln_ptr[i_blk].Grid.INl;
	        i_max0 = Soln_ptr[i_blk].Grid.INl;
	        i_inc = 1;
	        j_min0 = Soln_ptr[i_blk].Grid.JNu-1;
	        j_max0 = Soln_ptr[i_blk].Grid.JNu-1;
	        j_inc = -1;
                i_ref0 = i_min0;
                j_ref0 = j_min0+1;
             } else {
	        i_min0 = Soln_ptr[i_blk].Grid.INl;
	        i_max0 = Soln_ptr[i_blk].Grid.INl;
	        i_inc = -1;
	        j_min0 = Soln_ptr[i_blk].Grid.JNu-1;
	        j_max0 = Soln_ptr[i_blk].Grid.JNu-1;
	        j_inc = -1;
                i_ref0 = i_min0;
                j_ref0 = j_min0+1;
             } /* endif */
             x_ref = Soln_ptr[i_blk].Grid.Node[i_ref0][j_ref0].X; // Reference node location.
             // Four different orderings to consider depending on the value of i_inc & j_inc.
             if (j_inc > 0) {
                if (i_inc > 0) {
	           // Evaluate SW sub (fine) cell SE node location.
                   for ( k = 1 ; k <= NUM_COMP_VECTOR2D; ++ k) {
                      l0 = l0 + 1;
                      if (l0 >= buffer_size_neighbour) return(4201);
                      Soln_Block_List.message_reschange_northwestcorner_sendbuf[i_blk][l0] =
                         HALF*(Soln_ptr[i_blk].Grid.Node[i_min0][j_min0].X[k]+
                               Soln_ptr[i_blk].Grid.Node[i_min0+1][j_min0].X[k])-x_ref[k];
                   } /* endfor */
	           // Evaluate SE sub (fine) cell SE node location.
                   for ( k = 1 ; k <= NUM_COMP_VECTOR2D; ++ k) {
                      l0 = l0 + 1;
                      if (l0 >= buffer_size_neighbour) return(4202);
                      Soln_Block_List.message_reschange_northwestcorner_sendbuf[i_blk][l0] =
                         Soln_ptr[i_blk].Grid.Node[i_min0+1][j_min0].X[k]-x_ref[k];
                   } /* endfor */
	           // Evaluate NW sub (fine) cell SE node location.
                   for ( k = 1 ; k <= NUM_COMP_VECTOR2D; ++ k) {
                      l0 = l0 + 1;
                      if (l0 >= buffer_size_neighbour) return(4203);
                      Soln_Block_List.message_reschange_northwestcorner_sendbuf[i_blk][l0] =
                         Soln_ptr[i_blk].Grid.Cell[i_min0][j_min0].Xc[k]-x_ref[k];
                   } /* endfor */
	           // Evaluate NE sub (fine) cell SE node location.
                   for ( k = 1 ; k <= NUM_COMP_VECTOR2D; ++ k) {
                      l0 = l0 + 1;
                      if (l0 >= buffer_size_neighbour) return(4204);
                      Soln_Block_List.message_reschange_northwestcorner_sendbuf[i_blk][l0] =
                         HALF*(Soln_ptr[i_blk].Grid.Node[i_min0+1][j_min0].X[k]+
                               Soln_ptr[i_blk].Grid.Node[i_min0+1][j_min0+1].X[k])-x_ref[k];
                   } /* endfor */
                } else {
	           // Evaluate SE sub (fine) cell SE node location.
                   for ( k = 1 ; k <= NUM_COMP_VECTOR2D; ++ k) {
                      l0 = l0 + 1;
                      if (l0 >= buffer_size_neighbour) return(4205);
                      Soln_Block_List.message_reschange_northwestcorner_sendbuf[i_blk][l0] =
                         Soln_ptr[i_blk].Grid.Node[i_min0+1][j_min0].X[k]-x_ref[k];
                   } /* endfor */
	           // Evaluate SW sub (fine) cell SE node location.
                   for ( k = 1 ; k <= NUM_COMP_VECTOR2D; ++ k) {
                      l0 = l0 + 1;
                      if (l0 >= buffer_size_neighbour) return(4206);
                      Soln_Block_List.message_reschange_northwestcorner_sendbuf[i_blk][l0] =
                         HALF*(Soln_ptr[i_blk].Grid.Node[i_min0][j_min0].X[k]+
                               Soln_ptr[i_blk].Grid.Node[i_min0+1][j_min0].X[k])-x_ref[k];
                   } /* endfor */
	           // Evaluate NE sub (fine) cell SE node location.
                   for ( k = 1 ; k <= NUM_COMP_VECTOR2D; ++ k) {
                      l0 = l0 + 1;
                      if (l0 >= buffer_size_neighbour) return(4207);
                      Soln_Block_List.message_reschange_northwestcorner_sendbuf[i_blk][l0] =
                         HALF*(Soln_ptr[i_blk].Grid.Node[i_min0+1][j_min0].X[k]+
                               Soln_ptr[i_blk].Grid.Node[i_min0+1][j_min0+1].X[k])-x_ref[k];
                   } /* endfor */
	           // Evaluate NW sub (fine) cell SE node location.
                   for ( k = 1 ; k <= NUM_COMP_VECTOR2D; ++ k) {
                      l0 = l0 + 1;
                      if (l0 >= buffer_size_neighbour) return(4208);
                      Soln_Block_List.message_reschange_northwestcorner_sendbuf[i_blk][l0] =
                         Soln_ptr[i_blk].Grid.Cell[i_min0][j_min0].Xc[k]-x_ref[k];
                   } /* endfor */
                } /* endif */
             } else {
                if (i_inc > 0) {
	           // Evaluate NW sub (fine) cell SE node location.
                   for ( k = 1 ; k <= NUM_COMP_VECTOR2D; ++ k) {
                      l0 = l0 + 1;
                      if (l0 >= buffer_size_neighbour) return(4209);
                      Soln_Block_List.message_reschange_northwestcorner_sendbuf[i_blk][l0] =
                         Soln_ptr[i_blk].Grid.Cell[i_min0][j_min0].Xc[k]-x_ref[k];
                   } /* endfor */
	           // Evaluate NE sub (fine) cell SE node location.
                   for ( k = 1 ; k <= NUM_COMP_VECTOR2D; ++ k) {
                      l0 = l0 + 1;
                      if (l0 >= buffer_size_neighbour) return(4210);
                      Soln_Block_List.message_reschange_northwestcorner_sendbuf[i_blk][l0] =
                         HALF*(Soln_ptr[i_blk].Grid.Node[i_min0+1][j_min0].X[k]+
                               Soln_ptr[i_blk].Grid.Node[i_min0+1][j_min0+1].X[k])-x_ref[k];
                   } /* endfor */
	           // Evaluate SW sub (fine) cell SE node location.
                   for ( k = 1 ; k <= NUM_COMP_VECTOR2D; ++ k) {
                      l0 = l0 + 1;
                      if (l0 >= buffer_size_neighbour) return(4211);
                      Soln_Block_List.message_reschange_northwestcorner_sendbuf[i_blk][l0] =
                         HALF*(Soln_ptr[i_blk].Grid.Node[i_min0][j_min0].X[k]+
                               Soln_ptr[i_blk].Grid.Node[i_min0+1][j_min0].X[k])-x_ref[k];
                   } /* endfor */
	           // Evaluate SE sub (fine) cell SE node location.
                   for ( k = 1 ; k <= NUM_COMP_VECTOR2D; ++ k) {
                      l0 = l0 + 1;
                      if (l0 >= buffer_size_neighbour) return(4212);
                      Soln_Block_List.message_reschange_northwestcorner_sendbuf[i_blk][l0] =
                         Soln_ptr[i_blk].Grid.Node[i_min0+1][j_min0].X[k]-x_ref[k];
                   } /* endfor */
                } else {
	           // Evaluate NE sub (fine) cell SE node location.
                   for ( k = 1 ; k <= NUM_COMP_VECTOR2D; ++ k) {
                      l0 = l0 + 1;
                      if (l0 >= buffer_size_neighbour) return(4213);
                      Soln_Block_List.message_reschange_northwestcorner_sendbuf[i_blk][l0] =
                         HALF*(Soln_ptr[i_blk].Grid.Node[i_min0+1][j_min0].X[k]+
                               Soln_ptr[i_blk].Grid.Node[i_min0+1][j_min0+1].X[k])-x_ref[k];
                   } /* endfor */
	           // Evaluate NW sub (fine) cell SE node location.
                   for ( k = 1 ; k <= NUM_COMP_VECTOR2D; ++ k) {
                      l0 = l0 + 1;
                      if (l0 >= buffer_size_neighbour) return(4214);
                      Soln_Block_List.message_reschange_northwestcorner_sendbuf[i_blk][l0] =
                         Soln_ptr[i_blk].Grid.Cell[i_min0][j_min0].Xc[k]-x_ref[k];
                   } /* endfor */
	           // Evaluate SE sub (fine) cell SE node location.
                   for ( k = 1 ; k <= NUM_COMP_VECTOR2D; ++ k) {
                      l0 = l0 + 1;
                      if (l0 >= buffer_size_neighbour) return(4215);
                      Soln_Block_List.message_reschange_northwestcorner_sendbuf[i_blk][l0] =
                         Soln_ptr[i_blk].Grid.Node[i_min0+1][j_min0].X[k]-x_ref[k];
                   } /* endfor */
	           // Evaluate SW sub (fine) cell SE node location.
                   for ( k = 1 ; k <= NUM_COMP_VECTOR2D; ++ k) {
                      l0 = l0 + 1;
                      if (l0 >= buffer_size_neighbour) return(4216);
                      Soln_Block_List.message_reschange_northwestcorner_sendbuf[i_blk][l0] =
                         HALF*(Soln_ptr[i_blk].Grid.Node[i_min0][j_min0].X[k]+
                               Soln_ptr[i_blk].Grid.Node[i_min0+1][j_min0].X[k])-x_ref[k];
                   } /* endfor */
                } /* endif */
             } /* endif */
          } /* endif */
       } /* endif */

       //////////////////////////////////////////////////////////////////
       // Load send buffers of solution blocks with refined North East //
       // corner neighbour at higher mesh resolution.                  //
       //////////////////////////////////////////////////////////////////
       if (Soln_Block_List.Block[i_blk].used && 
           (Soln_Block_List.Block[i_blk].nNE == 1) &&
           (Soln_Block_List.Block[i_blk].info.level < 
            Soln_Block_List.Block[i_blk].infoNE[0].level)) {
          buffer_size_neighbour = sqr(Soln_Block_List.Block[i_blk].infoNE[0].dimen.ghost)*
                                  Number_of_Solution_Variables;
          l0 = -1;
          // Load ghost cell solution variable information as required.
          if (!Send_Mesh_Geometry_Only) {
             if (Soln_Block_List.Block[i_blk].infoNE[0].dimen.i > 0 &&
                 Soln_Block_List.Block[i_blk].infoNE[0].dimen.j > 0) {
	        i_min0 = Soln_ptr[i_blk].ICu;
                i_max0 = Soln_ptr[i_blk].ICu;
	        i_inc = 1;
	        j_min0 = Soln_ptr[i_blk].JCu;
	        j_max0 = Soln_ptr[i_blk].JCu;
	        j_inc = 1;
             } else if (Soln_Block_List.Block[i_blk].infoNE[0].dimen.j > 0) {
	        i_min0 = Soln_ptr[i_blk].ICu;
	        i_max0 = Soln_ptr[i_blk].ICu;
	        i_inc = -1;
	        j_min0 = Soln_ptr[i_blk].JCu;
	        j_max0 = Soln_ptr[i_blk].JCu;
	        j_inc = 1;
             } else if (Soln_Block_List.Block[i_blk].infoNE[0].dimen.i > 0) {
	        i_min0 = Soln_ptr[i_blk].ICu;
	        i_max0 = Soln_ptr[i_blk].ICu;
	        i_inc = 1;
	        j_min0 = Soln_ptr[i_blk].JCu;
	        j_max0 = Soln_ptr[i_blk].JCu;
	        j_inc = -1;
             } else {
	        i_min0 = Soln_ptr[i_blk].ICu;
	        i_max0 = Soln_ptr[i_blk].ICu;
	        i_inc = -1;
	        j_min0 = Soln_ptr[i_blk].JCu;
	        j_max0 = Soln_ptr[i_blk].JCu;
	        j_inc = -1;
             } /* endif */
             i = Soln_ptr[i_blk].LoadSendBuffer_C2F(Soln_Block_List.message_reschange_northeastcorner_sendbuf[i_blk],
                                                    l0,buffer_size_neighbour,
                                                    i_min0,i_max0,i_inc,
                                                    j_min0,j_max0,j_inc);
	     if (i != 0) return(5200);
          } /* endif */
          // Load ghost cell mesh information as required.
          if (Send_Mesh_Geometry_Only ||
              Number_of_Solution_Variables > Soln_ptr[i_blk].NumVar()) {
             if (Soln_Block_List.Block[i_blk].infoNE[0].dimen.i > 0 &&
                 Soln_Block_List.Block[i_blk].infoNE[0].dimen.j > 0) {
	        i_min0 = Soln_ptr[i_blk].Grid.INu-1;
	        i_max0 = Soln_ptr[i_blk].Grid.INu-1;
	        i_inc = 1;
	        j_min0 = Soln_ptr[i_blk].Grid.JNu-1;
	        j_max0 = Soln_ptr[i_blk].Grid.JNu-1;
	        j_inc = 1;
                i_ref0 = i_min0+1;
                j_ref0 = j_min0+1;
             } else if (Soln_Block_List.Block[i_blk].infoNE[0].dimen.j > 0) {
	        i_min0 = Soln_ptr[i_blk].Grid.INu-1;
	        i_max0 = Soln_ptr[i_blk].Grid.INu-1;
	        i_inc = -1;
	        j_min0 = Soln_ptr[i_blk].Grid.JNu-1;
	        j_max0 = Soln_ptr[i_blk].Grid.JNu-1;
	        j_inc = 1;
                i_ref0 = i_min0+1;
                j_ref0 = j_min0+1;
             } else if (Soln_Block_List.Block[i_blk].infoNE[0].dimen.i > 0) {
	        i_min0 = Soln_ptr[i_blk].Grid.INu-1;
	        i_max0 = Soln_ptr[i_blk].Grid.INu-1;
	        i_inc = 1;
	        j_min0 = Soln_ptr[i_blk].Grid.JNu-1;
	        j_max0 = Soln_ptr[i_blk].Grid.JNu-1;
	        j_inc = -1;
                i_ref0 = i_min0+1;
                j_ref0 = j_min0+1;
             } else {
	        i_min0 = Soln_ptr[i_blk].Grid.INu-1;
	        i_max0 = Soln_ptr[i_blk].Grid.INu-1;
	        i_inc = -1;
	        j_min0 = Soln_ptr[i_blk].Grid.JNu-1;
	        j_max0 = Soln_ptr[i_blk].Grid.JNu-1;
	        j_inc = -1;
                i_ref0 = i_min0+1;
                j_ref0 = j_min0+1;
             } /* endif */
             x_ref = Soln_ptr[i_blk].Grid.Node[i_ref0][j_ref0].X; // Reference node location.
             // Four different orderings to consider depending on the value of i_inc & j_inc.
             if (j_inc > 0) {
                if (i_inc > 0) {
		  // Evaluate SW sub (fine) cell SW node location.
                   for ( k = 1 ; k <= NUM_COMP_VECTOR2D; ++ k) {
                      l0 = l0 + 1;
                      if (l0 >= buffer_size_neighbour) return(5201);
                      Soln_Block_List.message_reschange_northeastcorner_sendbuf[i_blk][l0] =
                         Soln_ptr[i_blk].Grid.Node[i_min0][j_min0].X[k]-x_ref[k];
                   } /* endfor */
	           // Evaluate SE sub (fine) cell SW node location.
                   for ( k = 1 ; k <= NUM_COMP_VECTOR2D; ++ k) {
                      l0 = l0 + 1;
                      if (l0 >= buffer_size_neighbour) return(5202);
                      Soln_Block_List.message_reschange_northeastcorner_sendbuf[i_blk][l0] =
                         HALF*(Soln_ptr[i_blk].Grid.Node[i_min0][j_min0].X[k]+
                               Soln_ptr[i_blk].Grid.Node[i_min0+1][j_min0].X[k])-x_ref[k];
                   } /* endfor */
	           // Evaluate NW sub (fine) cell SW node location.
                   for ( k = 1 ; k <= NUM_COMP_VECTOR2D; ++ k) {
                      l0 = l0 + 1;
                      if (l0 >= buffer_size_neighbour) return(5203);
                      Soln_Block_List.message_reschange_northeastcorner_sendbuf[i_blk][l0] =
                         HALF*(Soln_ptr[i_blk].Grid.Node[i_min0][j_min0].X[k]+
                               Soln_ptr[i_blk].Grid.Node[i_min0][j_min0+1].X[k])-x_ref[k];
                   } /* endfor */
	           // Evaluate NE sub (fine) cell SW node location.
                   for ( k = 1 ; k <= NUM_COMP_VECTOR2D; ++ k) {
                      l0 = l0 + 1;
                      if (l0 >= buffer_size_neighbour) return(5204);
                      Soln_Block_List.message_reschange_northeastcorner_sendbuf[i_blk][l0] =
                         Soln_ptr[i_blk].Grid.Cell[i_min0][j_min0].Xc[k]-x_ref[k];
                   } /* endfor */
                } else {
	           // Evaluate SE sub (fine) cell SW node location.
                   for ( k = 1 ; k <= NUM_COMP_VECTOR2D; ++ k) {
                      l0 = l0 + 1;
                      if (l0 >= buffer_size_neighbour) return(5205);
                      Soln_Block_List.message_reschange_northeastcorner_sendbuf[i_blk][l0] =
                         HALF*(Soln_ptr[i_blk].Grid.Node[i_min0][j_min0].X[k]+
                               Soln_ptr[i_blk].Grid.Node[i_min0+1][j_min0].X[k])-x_ref[k];
                   } /* endfor */
		  // Evaluate SW sub (fine) cell SW node location.
                   for ( k = 1 ; k <= NUM_COMP_VECTOR2D; ++ k) {
                      l0 = l0 + 1;
                      if (l0 >= buffer_size_neighbour) return(5206);
                      Soln_Block_List.message_reschange_northeastcorner_sendbuf[i_blk][l0] =
                         Soln_ptr[i_blk].Grid.Node[i_min0][j_min0].X[k]-x_ref[k];
                   } /* endfor */
	           // Evaluate NE sub (fine) cell SW node location.
                   for ( k = 1 ; k <= NUM_COMP_VECTOR2D; ++ k) {
                      l0 = l0 + 1;
                      if (l0 >= buffer_size_neighbour) return(5207);
                      Soln_Block_List.message_reschange_northeastcorner_sendbuf[i_blk][l0] =
                         Soln_ptr[i_blk].Grid.Cell[i_min0][j_min0].Xc[k]-x_ref[k];
                   } /* endfor */
	           // Evaluate NW sub (fine) cell SW node location.
                   for ( k = 1 ; k <= NUM_COMP_VECTOR2D; ++ k) {
                      l0 = l0 + 1;
                      if (l0 >= buffer_size_neighbour) return(5208);
                      Soln_Block_List.message_reschange_northeastcorner_sendbuf[i_blk][l0] =
                         HALF*(Soln_ptr[i_blk].Grid.Node[i_min0][j_min0].X[k]+
                               Soln_ptr[i_blk].Grid.Node[i_min0][j_min0+1].X[k])-x_ref[k];
                   } /* endfor */
                } /* endif */
             } else {
                if (i_inc > 0) {
	           // Evaluate NW sub (fine) cell SW node location.
                   for ( k = 1 ; k <= NUM_COMP_VECTOR2D; ++ k) {
                      l0 = l0 + 1;
                      if (l0 >= buffer_size_neighbour) return(5209);
                      Soln_Block_List.message_reschange_northeastcorner_sendbuf[i_blk][l0] =
                         HALF*(Soln_ptr[i_blk].Grid.Node[i_min0][j_min0].X[k]+
                               Soln_ptr[i_blk].Grid.Node[i_min0][j_min0+1].X[k])-x_ref[k];
                   } /* endfor */
	           // Evaluate NE sub (fine) cell SW node location.
                   for ( k = 1 ; k <= NUM_COMP_VECTOR2D; ++ k) {
                      l0 = l0 + 1;
                      if (l0 >= buffer_size_neighbour) return(5210);
                      Soln_Block_List.message_reschange_northeastcorner_sendbuf[i_blk][l0] =
                         Soln_ptr[i_blk].Grid.Cell[i_min0][j_min0].Xc[k]-x_ref[k];
                   } /* endfor */
		  // Evaluate SW sub (fine) cell SW node location.
                   for ( k = 1 ; k <= NUM_COMP_VECTOR2D; ++ k) {
                      l0 = l0 + 1;
                      if (l0 >= buffer_size_neighbour) return(5211);
                      Soln_Block_List.message_reschange_northeastcorner_sendbuf[i_blk][l0] =
                         Soln_ptr[i_blk].Grid.Node[i_min0][j_min0].X[k]-x_ref[k];
                   } /* endfor */
	           // Evaluate SE sub (fine) cell SW node location.
                   for ( k = 1 ; k <= NUM_COMP_VECTOR2D; ++ k) {
                      l0 = l0 + 1;
                      if (l0 >= buffer_size_neighbour) return(5212);
                      Soln_Block_List.message_reschange_northeastcorner_sendbuf[i_blk][l0] =
                         HALF*(Soln_ptr[i_blk].Grid.Node[i_min0][j_min0].X[k]+
                               Soln_ptr[i_blk].Grid.Node[i_min0+1][j_min0].X[k])-x_ref[k];
                   } /* endfor */
                } else {
	           // Evaluate NE sub (fine) cell SW node location.
                   for ( k = 1 ; k <= NUM_COMP_VECTOR2D; ++ k) {
                      l0 = l0 + 1;
                      if (l0 >= buffer_size_neighbour) return(5213);
                      Soln_Block_List.message_reschange_northeastcorner_sendbuf[i_blk][l0] =
                         Soln_ptr[i_blk].Grid.Cell[i_min0][j_min0].Xc[k]-x_ref[k];
                   } /* endfor */
	           // Evaluate NW sub (fine) cell SW node location.
                   for ( k = 1 ; k <= NUM_COMP_VECTOR2D; ++ k) {
                      l0 = l0 + 1;
                      if (l0 >= buffer_size_neighbour) return(5214);
                      Soln_Block_List.message_reschange_northeastcorner_sendbuf[i_blk][l0] =
                         HALF*(Soln_ptr[i_blk].Grid.Node[i_min0][j_min0].X[k]+
                               Soln_ptr[i_blk].Grid.Node[i_min0][j_min0+1].X[k])-x_ref[k];
                   } /* endfor */
	           // Evaluate SE sub (fine) cell SW node location.
                   for ( k = 1 ; k <= NUM_COMP_VECTOR2D; ++ k) {
                      l0 = l0 + 1;
                      if (l0 >= buffer_size_neighbour) return(5215);
                      Soln_Block_List.message_reschange_northeastcorner_sendbuf[i_blk][l0] =
                         HALF*(Soln_ptr[i_blk].Grid.Node[i_min0][j_min0].X[k]+
                               Soln_ptr[i_blk].Grid.Node[i_min0+1][j_min0].X[k])-x_ref[k];
                   } /* endfor */
		  // Evaluate SW sub (fine) cell SW node location.
                   for ( k = 1 ; k <= NUM_COMP_VECTOR2D; ++ k) {
                      l0 = l0 + 1;
                      if (l0 >= buffer_size_neighbour) return(5216);
                      Soln_Block_List.message_reschange_northeastcorner_sendbuf[i_blk][l0] =
                         Soln_ptr[i_blk].Grid.Node[i_min0][j_min0].X[k]-x_ref[k];
                   } /* endfor */
                } /* endif */
             } /* endif */
          } /* endif */
       } /* endif */

       //////////////////////////////////////////////////////////////////
       // Load send buffers of solution blocks with refined South East //
       // corner neighbour at higher mesh resolution.                  //
       //////////////////////////////////////////////////////////////////
       if (Soln_Block_List.Block[i_blk].used && 
           (Soln_Block_List.Block[i_blk].nSE == 1) &&
           (Soln_Block_List.Block[i_blk].info.level < 
            Soln_Block_List.Block[i_blk].infoSE[0].level)) {
	 buffer_size_neighbour = sqr(Soln_Block_List.Block[i_blk].infoSE[0].dimen.ghost)*
                                  Number_of_Solution_Variables;
          l0 = -1;
          // Load ghost cell solution variable information as required.
          if (!Send_Mesh_Geometry_Only) {
             if (Soln_Block_List.Block[i_blk].infoSE[0].dimen.i > 0 &&
                 Soln_Block_List.Block[i_blk].infoSE[0].dimen.j > 0) {
	        i_min0 = Soln_ptr[i_blk].ICu;
                i_max0 = Soln_ptr[i_blk].ICu;
	        i_inc = 1;
	        j_min0 = Soln_ptr[i_blk].JCl;
	        j_max0 = Soln_ptr[i_blk].JCl;
	        j_inc = 1;
             } else if (Soln_Block_List.Block[i_blk].infoSE[0].dimen.j > 0) {
	        i_min0 = Soln_ptr[i_blk].ICu;
	        i_max0 = Soln_ptr[i_blk].ICu;
	        i_inc = -1;
	        j_min0 = Soln_ptr[i_blk].JCl;
	        j_max0 = Soln_ptr[i_blk].JCl;
	        j_inc = 1;
             } else if (Soln_Block_List.Block[i_blk].infoSE[0].dimen.i > 0) {
	        i_min0 = Soln_ptr[i_blk].ICu;
	        i_max0 = Soln_ptr[i_blk].ICu;
	        i_inc = 1;
	        j_min0 = Soln_ptr[i_blk].JCl;
	        j_max0 = Soln_ptr[i_blk].JCl;
	        j_inc = -1;
             } else {
	        i_min0 = Soln_ptr[i_blk].ICu;
	        i_max0 = Soln_ptr[i_blk].ICu;
	        i_inc = -1;
	        j_min0 = Soln_ptr[i_blk].JCl;
	        j_max0 = Soln_ptr[i_blk].JCl;
	        j_inc = -1;
             } /* endif */
             i = Soln_ptr[i_blk].LoadSendBuffer_C2F(Soln_Block_List.message_reschange_southeastcorner_sendbuf[i_blk],
                                                    l0,buffer_size_neighbour,
                                                    i_min0,i_max0,i_inc,
                                                    j_min0,j_max0,j_inc);
	     if (i != 0) return(4100);
          } /* endif */
          // Load ghost cell mesh information as required.
          if (Send_Mesh_Geometry_Only ||
              Number_of_Solution_Variables > Soln_ptr[i_blk].NumVar()) {
             if (Soln_Block_List.Block[i_blk].infoSE[0].dimen.i > 0 &&
                 Soln_Block_List.Block[i_blk].infoSE[0].dimen.j > 0) {
	        i_min0 = Soln_ptr[i_blk].Grid.INu-1;
	        i_max0 = Soln_ptr[i_blk].Grid.INu-1;
	        i_inc = 1;
	        j_min0 = Soln_ptr[i_blk].Grid.JNl;
	        j_max0 = Soln_ptr[i_blk].Grid.JNl;
	        j_inc = 1;
                i_ref0 = i_min0+1;
                j_ref0 = j_min0;
             } else if (Soln_Block_List.Block[i_blk].infoSE[0].dimen.j > 0) {
	        i_min0 = Soln_ptr[i_blk].Grid.INu-1;
	        i_max0 = Soln_ptr[i_blk].Grid.INu-1;
	        i_inc = -1;
	        j_min0 = Soln_ptr[i_blk].Grid.JNl;
	        j_max0 = Soln_ptr[i_blk].Grid.JNl;
	        j_inc = 1;
                i_ref0 = i_min0+1;
                j_ref0 = j_min0;
             } else if (Soln_Block_List.Block[i_blk].infoSE[0].dimen.i > 0) {
	        i_min0 = Soln_ptr[i_blk].Grid.INu-1;
	        i_max0 = Soln_ptr[i_blk].Grid.INu-1;
	        i_inc = 1;
	        j_min0 = Soln_ptr[i_blk].Grid.JNl;
	        j_max0 = Soln_ptr[i_blk].Grid.JNl;
	        j_inc = -1;
                i_ref0 = i_min0+1;
                j_ref0 = j_min0;
             } else {
	        i_min0 = Soln_ptr[i_blk].Grid.INu-1;
	        i_max0 = Soln_ptr[i_blk].Grid.INu-1;
	        i_inc = -1;
	        j_min0 = Soln_ptr[i_blk].Grid.JNl;
	        j_max0 = Soln_ptr[i_blk].Grid.JNl;
	        j_inc = -1;
                i_ref0 = i_min0+1;
                j_ref0 = j_min0;
             } /* endif */
             x_ref = Soln_ptr[i_blk].Grid.Node[i_ref0][j_ref0].X; // Reference node location.
             // Four different orderings to consider depending on the value of i_inc & j_inc.
             if (j_inc > 0) {
                if (i_inc > 0) {
		   // Evaluate SW sub (fine) cell NW node location.
                   for ( k = 1 ; k <= NUM_COMP_VECTOR2D; ++ k) {
                      l0 = l0 + 1;
                      if (l0 >= buffer_size_neighbour) return(4101);
                      Soln_Block_List.message_reschange_southeastcorner_sendbuf[i_blk][l0] =
                         HALF*(Soln_ptr[i_blk].Grid.Node[i_min0][j_min0].X[k]+
                               Soln_ptr[i_blk].Grid.Node[i_min0][j_min0+1].X[k])-x_ref[k];
                   } /* endfor */
	           // Evaluate SE sub (fine) cell NW node location.
                   for ( k = 1 ; k <= NUM_COMP_VECTOR2D; ++ k) {
                      l0 = l0 + 1;
                      if (l0 >= buffer_size_neighbour) return(4102);
                      Soln_Block_List.message_reschange_southeastcorner_sendbuf[i_blk][l0] =
                         Soln_ptr[i_blk].Grid.Cell[i_min0][j_min0].Xc[k]-x_ref[k];
                   } /* endfor */
	           // Evaluate NW sub (fine) cell NW node location.
                   for ( k = 1 ; k <= NUM_COMP_VECTOR2D; ++ k) {
                      l0 = l0 + 1;
                      if (l0 >= buffer_size_neighbour) return(4103);
                      Soln_Block_List.message_reschange_southeastcorner_sendbuf[i_blk][l0] =
                         Soln_ptr[i_blk].Grid.Node[i_min0][j_min0+1].X[k]-x_ref[k];
                   } /* endfor */
	           // Evaluate NE sub (fine) cell NW node location.
                   for ( k = 1 ; k <= NUM_COMP_VECTOR2D; ++ k) {
                      l0 = l0 + 1;
                      if (l0 >= buffer_size_neighbour) return(4104);
                      Soln_Block_List.message_reschange_southeastcorner_sendbuf[i_blk][l0] =
                         HALF*(Soln_ptr[i_blk].Grid.Node[i_min0][j_min0+1].X[k]+
                               Soln_ptr[i_blk].Grid.Node[i_min0+1][j_min0+1].X[k])-x_ref[k];
                   } /* endfor */
                } else {
	           // Evaluate SE sub (fine) cell NW node location.
                   for ( k = 1 ; k <= NUM_COMP_VECTOR2D; ++ k) {
                      l0 = l0 + 1;
                      if (l0 >= buffer_size_neighbour) return(4105);
                      Soln_Block_List.message_reschange_southeastcorner_sendbuf[i_blk][l0] =
                         Soln_ptr[i_blk].Grid.Cell[i_min0][j_min0].Xc[k]-x_ref[k];
                   } /* endfor */
		   // Evaluate SW sub (fine) cell NW node location.
                   for ( k = 1 ; k <= NUM_COMP_VECTOR2D; ++ k) {
                      l0 = l0 + 1;
                      if (l0 >= buffer_size_neighbour) return(4106);
                      Soln_Block_List.message_reschange_southeastcorner_sendbuf[i_blk][l0] =
                         HALF*(Soln_ptr[i_blk].Grid.Node[i_min0][j_min0].X[k]+
                               Soln_ptr[i_blk].Grid.Node[i_min0][j_min0+1].X[k])-x_ref[k];
                   } /* endfor */
	           // Evaluate NE sub (fine) cell NW node location.
                   for ( k = 1 ; k <= NUM_COMP_VECTOR2D; ++ k) {
                      l0 = l0 + 1;
                      if (l0 >= buffer_size_neighbour) return(4107);
                      Soln_Block_List.message_reschange_southeastcorner_sendbuf[i_blk][l0] =
                         HALF*(Soln_ptr[i_blk].Grid.Node[i_min0][j_min0+1].X[k]+
                               Soln_ptr[i_blk].Grid.Node[i_min0+1][j_min0+1].X[k])-x_ref[k];
                   } /* endfor */
	           // Evaluate NW sub (fine) cell NW node location.
                   for ( k = 1 ; k <= NUM_COMP_VECTOR2D; ++ k) {
                      l0 = l0 + 1;
                      if (l0 >= buffer_size_neighbour) return(4108);
                      Soln_Block_List.message_reschange_southeastcorner_sendbuf[i_blk][l0] =
                         Soln_ptr[i_blk].Grid.Node[i_min0][j_min0+1].X[k]-x_ref[k];
                   } /* endfor */
                } /* endif */
             } else {
                if (i_inc > 0) {
	           // Evaluate NW sub (fine) cell NW node location.
                   for ( k = 1 ; k <= NUM_COMP_VECTOR2D; ++ k) {
                      l0 = l0 + 1;
                      if (l0 >= buffer_size_neighbour) return(4109);
                      Soln_Block_List.message_reschange_southeastcorner_sendbuf[i_blk][l0] =
                         Soln_ptr[i_blk].Grid.Node[i_min0][j_min0+1].X[k]-x_ref[k];
                   } /* endfor */
	           // Evaluate NE sub (fine) cell NW node location.
                   for ( k = 1 ; k <= NUM_COMP_VECTOR2D; ++ k) {
                      l0 = l0 + 1;
                      if (l0 >= buffer_size_neighbour) return(4110);
                      Soln_Block_List.message_reschange_southeastcorner_sendbuf[i_blk][l0] =
                         HALF*(Soln_ptr[i_blk].Grid.Node[i_min0][j_min0+1].X[k]+
                               Soln_ptr[i_blk].Grid.Node[i_min0+1][j_min0+1].X[k])-x_ref[k];
                   } /* endfor */
		   // Evaluate SW sub (fine) cell NW node location.
                   for ( k = 1 ; k <= NUM_COMP_VECTOR2D; ++ k) {
                      l0 = l0 + 1;
                      if (l0 >= buffer_size_neighbour) return(4111);
                      Soln_Block_List.message_reschange_southeastcorner_sendbuf[i_blk][l0] =
                         HALF*(Soln_ptr[i_blk].Grid.Node[i_min0][j_min0].X[k]+
                               Soln_ptr[i_blk].Grid.Node[i_min0][j_min0+1].X[k])-x_ref[k];
                   } /* endfor */
	           // Evaluate SE sub (fine) cell NW node location.
                   for ( k = 1 ; k <= NUM_COMP_VECTOR2D; ++ k) {
                      l0 = l0 + 1;
                      if (l0 >= buffer_size_neighbour) return(4112);
                      Soln_Block_List.message_reschange_southeastcorner_sendbuf[i_blk][l0] =
                         Soln_ptr[i_blk].Grid.Cell[i_min0][j_min0].Xc[k]-x_ref[k];
                   } /* endfor */
                } else {
	           // Evaluate NE sub (fine) cell NW node location.
                   for ( k = 1 ; k <= NUM_COMP_VECTOR2D; ++ k) {
                      l0 = l0 + 1;
                      if (l0 >= buffer_size_neighbour) return(4113);
                      Soln_Block_List.message_reschange_southeastcorner_sendbuf[i_blk][l0] =
                         HALF*(Soln_ptr[i_blk].Grid.Node[i_min0][j_min0+1].X[k]+
                               Soln_ptr[i_blk].Grid.Node[i_min0+1][j_min0+1].X[k])-x_ref[k];
                   } /* endfor */
	           // Evaluate NW sub (fine) cell NW node location.
                   for ( k = 1 ; k <= NUM_COMP_VECTOR2D; ++ k) {
                      l0 = l0 + 1;
                      if (l0 >= buffer_size_neighbour) return(4114);
                      Soln_Block_List.message_reschange_southeastcorner_sendbuf[i_blk][l0] =
                         Soln_ptr[i_blk].Grid.Node[i_min0][j_min0+1].X[k]-x_ref[k];
                   } /* endfor */
	           // Evaluate SE sub (fine) cell NW node location.
                   for ( k = 1 ; k <= NUM_COMP_VECTOR2D; ++ k) {
                      l0 = l0 + 1;
                      if (l0 >= buffer_size_neighbour) return(4115);
                      Soln_Block_List.message_reschange_southeastcorner_sendbuf[i_blk][l0] =
                         Soln_ptr[i_blk].Grid.Cell[i_min0][j_min0].Xc[k]-x_ref[k];
                   } /* endfor */
		   // Evaluate SW sub (fine) cell NW node location.
                   for ( k = 1 ; k <= NUM_COMP_VECTOR2D; ++ k) {
                      l0 = l0 + 1;
                      if (l0 >= buffer_size_neighbour) return(4116);
                      Soln_Block_List.message_reschange_southeastcorner_sendbuf[i_blk][l0] =
                         HALF*(Soln_ptr[i_blk].Grid.Node[i_min0][j_min0].X[k]+
                               Soln_ptr[i_blk].Grid.Node[i_min0][j_min0+1].X[k])-x_ref[k];
                   } /* endfor */
                } /* endif */
             } /* endif */
          } /* endif */
       } /* endif */

       //////////////////////////////////////////////////////////////////
       // Load send buffers of solution blocks with refined South West //
       // corner neighbour at higher mesh resolution.                  //
       //////////////////////////////////////////////////////////////////
       if (Soln_Block_List.Block[i_blk].used && 
           (Soln_Block_List.Block[i_blk].nSW == 1) &&
           (Soln_Block_List.Block[i_blk].info.level < 
            Soln_Block_List.Block[i_blk].infoSW[0].level)) {
	  buffer_size_neighbour = sqr(Soln_Block_List.Block[i_blk].infoSW[0].dimen.ghost)*
                                  Number_of_Solution_Variables;
          l0 = -1;
          // Load ghost cell solution variable information as required.
          if (!Send_Mesh_Geometry_Only) {
             if (Soln_Block_List.Block[i_blk].infoSW[0].dimen.i > 0 &&
                 Soln_Block_List.Block[i_blk].infoSW[0].dimen.j > 0) {
	        i_min0 = Soln_ptr[i_blk].ICl;
                i_max0 = Soln_ptr[i_blk].ICl;
	        i_inc = 1;
	        j_min0 = Soln_ptr[i_blk].JCl;
	        j_max0 = Soln_ptr[i_blk].JCl;
	        j_inc = 1;
             } else if (Soln_Block_List.Block[i_blk].infoSW[0].dimen.j > 0) {
	        i_min0 = Soln_ptr[i_blk].ICl;
	        i_max0 = Soln_ptr[i_blk].ICl;
	        i_inc = -1;
	        j_min0 = Soln_ptr[i_blk].JCl;
	        j_max0 = Soln_ptr[i_blk].JCl;
	        j_inc = 1;
             } else if (Soln_Block_List.Block[i_blk].infoSW[0].dimen.i > 0) {
	        i_min0 = Soln_ptr[i_blk].ICl;
	        i_max0 = Soln_ptr[i_blk].ICl;
	        i_inc = 1;
	        j_min0 = Soln_ptr[i_blk].JCl;
	        j_max0 = Soln_ptr[i_blk].JCl;
	        j_inc = -1;
             } else {
	        i_min0 = Soln_ptr[i_blk].ICl;
	        i_max0 = Soln_ptr[i_blk].ICl;
	        i_inc = -1;
	        j_min0 = Soln_ptr[i_blk].JCl;
	        j_max0 = Soln_ptr[i_blk].JCl;
	        j_inc = -1;
             } /* endif */
             i = Soln_ptr[i_blk].LoadSendBuffer_C2F(Soln_Block_List.message_reschange_southwestcorner_sendbuf[i_blk],
                                                    l0,buffer_size_neighbour,
                                                    i_min0,i_max0,i_inc,
                                                    j_min0,j_max0,j_inc);
	     if (i != 0) return(5100);
          } /* endif */
          // Load ghost cell mesh information as required.
          if (Send_Mesh_Geometry_Only ||
              Number_of_Solution_Variables > Soln_ptr[i_blk].NumVar()) {
             if (Soln_Block_List.Block[i_blk].infoSW[0].dimen.i > 0 &&
                 Soln_Block_List.Block[i_blk].infoSW[0].dimen.j > 0) {
	        i_min0 = Soln_ptr[i_blk].Grid.INl;
	        i_max0 = Soln_ptr[i_blk].Grid.INl;
	        i_inc = 1;
	        j_min0 = Soln_ptr[i_blk].Grid.JNl;
	        j_max0 = Soln_ptr[i_blk].Grid.JNl;
	        j_inc = 1;
                i_ref0 = i_min0;
                j_ref0 = j_min0;
             } else if (Soln_Block_List.Block[i_blk].infoSW[0].dimen.j > 0) {
	        i_min0 = Soln_ptr[i_blk].Grid.INl;
	        i_max0 = Soln_ptr[i_blk].Grid.INl;
	        i_inc = -1;
	        j_min0 = Soln_ptr[i_blk].Grid.JNl;
	        j_max0 = Soln_ptr[i_blk].Grid.JNl;
	        j_inc = 1;
                i_ref0 = i_min0;
                j_ref0 = j_min0;
             } else if (Soln_Block_List.Block[i_blk].infoSW[0].dimen.i > 0) {
	        i_min0 = Soln_ptr[i_blk].Grid.INl;
	        i_max0 = Soln_ptr[i_blk].Grid.INl;
	        i_inc = 1;
	        j_min0 = Soln_ptr[i_blk].Grid.JNl;
	        j_max0 = Soln_ptr[i_blk].Grid.JNl;
	        j_inc = -1;
                i_ref0 = i_min0;
                j_ref0 = j_min0;
             } else {
	        i_min0 = Soln_ptr[i_blk].Grid.INl;
	        i_max0 = Soln_ptr[i_blk].Grid.INl;
	        i_inc = -1;
	        j_min0 = Soln_ptr[i_blk].Grid.JNl;
	        j_max0 = Soln_ptr[i_blk].Grid.JNl;
	        j_inc = -1;
                i_ref0 = i_min0;
                j_ref0 = j_min0;
             } /* endif */
             x_ref = Soln_ptr[i_blk].Grid.Node[i_ref0][j_ref0].X; // Reference node location.
             // Four different orderings to consider depending on the value of i_inc & j_inc.
             if (j_inc > 0) {
                if (i_inc > 0) {
		   // Evaluate SW sub (fine) cell NE node location.
                   for ( k = 1 ; k <= NUM_COMP_VECTOR2D; ++ k) {
                      l0 = l0 + 1;
                      if (l0 >= buffer_size_neighbour) return(5101);
                      Soln_Block_List.message_reschange_southwestcorner_sendbuf[i_blk][l0] =
                         Soln_ptr[i_blk].Grid.Cell[i_min0][j_min0].Xc[k]-x_ref[k];
                   } /* endfor */
	           // Evaluate SE sub (fine) cell NE node location.
                   for ( k = 1 ; k <= NUM_COMP_VECTOR2D; ++ k) {
                      l0 = l0 + 1;
                      if (l0 >= buffer_size_neighbour) return(5102);
                      Soln_Block_List.message_reschange_southwestcorner_sendbuf[i_blk][l0] =
                         HALF*(Soln_ptr[i_blk].Grid.Node[i_min0+1][j_min0].X[k]+
                               Soln_ptr[i_blk].Grid.Node[i_min0+1][j_min0+1].X[k])-x_ref[k];
                   } /* endfor */
	           // Evaluate NW sub (fine) cell NE node location.
                   for ( k = 1 ; k <= NUM_COMP_VECTOR2D; ++ k) {
                      l0 = l0 + 1;
                      if (l0 >= buffer_size_neighbour) return(5103);
                      Soln_Block_List.message_reschange_southwestcorner_sendbuf[i_blk][l0] =
                         HALF*(Soln_ptr[i_blk].Grid.Node[i_min0][j_min0+1].X[k]+
                               Soln_ptr[i_blk].Grid.Node[i_min0+1][j_min0+1].X[k])-x_ref[k];
                   } /* endfor */
	           // Evaluate NE sub (fine) cell NE node location.
                   for ( k = 1 ; k <= NUM_COMP_VECTOR2D; ++ k) {
                      l0 = l0 + 1;
                      if (l0 >= buffer_size_neighbour) return(5104);
                      Soln_Block_List.message_reschange_southwestcorner_sendbuf[i_blk][l0] =
                         Soln_ptr[i_blk].Grid.Node[i_min0+1][j_min0+1].X[k]-x_ref[k];
                   } /* endfor */
                } else {
	           // Evaluate SE sub (fine) cell NE node location.
                   for ( k = 1 ; k <= NUM_COMP_VECTOR2D; ++ k) {
                      l0 = l0 + 1;
                      if (l0 >= buffer_size_neighbour) return(5105);
                      Soln_Block_List.message_reschange_southwestcorner_sendbuf[i_blk][l0] =
                         HALF*(Soln_ptr[i_blk].Grid.Node[i_min0+1][j_min0].X[k]+
                               Soln_ptr[i_blk].Grid.Node[i_min0+1][j_min0+1].X[k])-x_ref[k];
                   } /* endfor */
		   // Evaluate SW sub (fine) cell NE node location.
                   for ( k = 1 ; k <= NUM_COMP_VECTOR2D; ++ k) {
                      l0 = l0 + 1;
                      if (l0 >= buffer_size_neighbour) return(5106);
                      Soln_Block_List.message_reschange_southwestcorner_sendbuf[i_blk][l0] =
                         Soln_ptr[i_blk].Grid.Cell[i_min0][j_min0].Xc[k]-x_ref[k];
                   } /* endfor */
	           // Evaluate NE sub (fine) cell NE node location.
                   for ( k = 1 ; k <= NUM_COMP_VECTOR2D; ++ k) {
                      l0 = l0 + 1;
                      if (l0 >= buffer_size_neighbour) return(5107);
                      Soln_Block_List.message_reschange_southwestcorner_sendbuf[i_blk][l0] =
                         Soln_ptr[i_blk].Grid.Node[i_min0+1][j_min0+1].X[k]-x_ref[k];
                   } /* endfor */
	           // Evaluate NW sub (fine) cell NE node location.
                   for ( k = 1 ; k <= NUM_COMP_VECTOR2D; ++ k) {
                      l0 = l0 + 1;
                      if (l0 >= buffer_size_neighbour) return(5108);
                      Soln_Block_List.message_reschange_southwestcorner_sendbuf[i_blk][l0] =
                         HALF*(Soln_ptr[i_blk].Grid.Node[i_min0][j_min0+1].X[k]+
                               Soln_ptr[i_blk].Grid.Node[i_min0+1][j_min0+1].X[k])-x_ref[k];
                   } /* endfor */
                } /* endif */
             } else {
                if (i_inc > 0) {
	           // Evaluate NW sub (fine) cell NE node location.
                   for ( k = 1 ; k <= NUM_COMP_VECTOR2D; ++ k) {
                      l0 = l0 + 1;
                      if (l0 >= buffer_size_neighbour) return(5109);
                      Soln_Block_List.message_reschange_southwestcorner_sendbuf[i_blk][l0] =
                         HALF*(Soln_ptr[i_blk].Grid.Node[i_min0][j_min0+1].X[k]+
                               Soln_ptr[i_blk].Grid.Node[i_min0+1][j_min0+1].X[k])-x_ref[k];
                   } /* endfor */
	           // Evaluate NE sub (fine) cell NE node location.
                   for ( k = 1 ; k <= NUM_COMP_VECTOR2D; ++ k) {
                      l0 = l0 + 1;
                      if (l0 >= buffer_size_neighbour) return(5110);
                      Soln_Block_List.message_reschange_southwestcorner_sendbuf[i_blk][l0] =
                         Soln_ptr[i_blk].Grid.Node[i_min0+1][j_min0+1].X[k]-x_ref[k];
                   } /* endfor */
		   // Evaluate SW sub (fine) cell NE node location.
                   for ( k = 1 ; k <= NUM_COMP_VECTOR2D; ++ k) {
                      l0 = l0 + 1;
                      if (l0 >= buffer_size_neighbour) return(5111);
                      Soln_Block_List.message_reschange_southwestcorner_sendbuf[i_blk][l0] =
                         Soln_ptr[i_blk].Grid.Cell[i_min0][j_min0].Xc[k]-x_ref[k];
                   } /* endfor */
	           // Evaluate SE sub (fine) cell NE node location.
                   for ( k = 1 ; k <= NUM_COMP_VECTOR2D; ++ k) {
                      l0 = l0 + 1;
                      if (l0 >= buffer_size_neighbour) return(5112);
                      Soln_Block_List.message_reschange_southwestcorner_sendbuf[i_blk][l0] =
                         HALF*(Soln_ptr[i_blk].Grid.Node[i_min0+1][j_min0].X[k]+
                               Soln_ptr[i_blk].Grid.Node[i_min0+1][j_min0+1].X[k])-x_ref[k];
                   } /* endfor */
                } else {
	           // Evaluate NE sub (fine) cell NE node location.
                   for ( k = 1 ; k <= NUM_COMP_VECTOR2D; ++ k) {
                      l0 = l0 + 1;
                      if (l0 >= buffer_size_neighbour) return(5113);
                      Soln_Block_List.message_reschange_southwestcorner_sendbuf[i_blk][l0] =
                         Soln_ptr[i_blk].Grid.Node[i_min0+1][j_min0+1].X[k]-x_ref[k];
                   } /* endfor */
	           // Evaluate NW sub (fine) cell NE node location.
                   for ( k = 1 ; k <= NUM_COMP_VECTOR2D; ++ k) {
                      l0 = l0 + 1;
                      if (l0 >= buffer_size_neighbour) return(5114);
                      Soln_Block_List.message_reschange_southwestcorner_sendbuf[i_blk][l0] =
                         HALF*(Soln_ptr[i_blk].Grid.Node[i_min0][j_min0+1].X[k]+
                               Soln_ptr[i_blk].Grid.Node[i_min0+1][j_min0+1].X[k])-x_ref[k];
                   } /* endfor */
	           // Evaluate SE sub (fine) cell NE node location.
                   for ( k = 1 ; k <= NUM_COMP_VECTOR2D; ++ k) {
                      l0 = l0 + 1;
                      if (l0 >= buffer_size_neighbour) return(5115);
                      Soln_Block_List.message_reschange_southwestcorner_sendbuf[i_blk][l0] =
                         HALF*(Soln_ptr[i_blk].Grid.Node[i_min0+1][j_min0].X[k]+
                               Soln_ptr[i_blk].Grid.Node[i_min0+1][j_min0+1].X[k])-x_ref[k];
                   } /* endfor */
		   // Evaluate SW sub (fine) cell NE node location.
                   for ( k = 1 ; k <= NUM_COMP_VECTOR2D; ++ k) {
                      l0 = l0 + 1;
                      if (l0 >= buffer_size_neighbour) return(5116);
                      Soln_Block_List.message_reschange_southwestcorner_sendbuf[i_blk][l0] =
                         Soln_ptr[i_blk].Grid.Cell[i_min0][j_min0].Xc[k]-x_ref[k];
                   } /* endfor */
                } /* endif */
             } /* endif */
          } /* endif */
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
template <class Quad_Soln_Block>
int Unload_Receive_Message_Buffers_NoResChange(Quad_Soln_Block *Soln_ptr,
                                               AdaptiveBlock2D_List &Soln_Block_List,
                                               const int Number_of_Solution_Variables,
                                               const int Send_Mesh_Geometry_Only) {

    int i_blk, buffer_size, 
        i_min, i_max, i_inc, i, 
        j_min, j_max, j_inc, j, 
        l;
    Vector2D x_ref;

    /* Check to see if the number of solution variables specified
       in the calling argument is correct. */

    if (Number_of_Solution_Variables != NUM_COMP_VECTOR2D &&
        Number_of_Solution_Variables != Soln_ptr[0].NumVar() &&
        Number_of_Solution_Variables != Soln_ptr[0].NumVar() + NUM_COMP_VECTOR2D) {
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
          if (!Send_Mesh_Geometry_Only) {
             buffer_size = Soln_Block_List.Block[i_blk].info.dimen.ghost*
                           abs(Soln_Block_List.Block[i_blk].info.dimen.i)*
                           Soln_ptr[i_blk].NumVar();
             if (Number_of_Solution_Variables > Soln_ptr[i_blk].NumVar()) 
                buffer_size = buffer_size +
                              Soln_Block_List.Block[i_blk].info.dimen.ghost*
                              (abs(Soln_Block_List.Block[i_blk].info.dimen.i)+1)*
                              NUM_COMP_VECTOR2D+
                              2*Soln_Block_List.Block[i_blk].info.dimen.ghost;
          } else {
                buffer_size = Soln_Block_List.Block[i_blk].info.dimen.ghost*
                              (abs(Soln_Block_List.Block[i_blk].info.dimen.i)+1)*
                              NUM_COMP_VECTOR2D+
                              2*Soln_Block_List.Block[i_blk].info.dimen.ghost;
          } /* endif */
          l = -1;
          // Unload ghost cell solution information as required.
          if (!Send_Mesh_Geometry_Only) {
             i_min = Soln_ptr[i_blk].ICl;
	     i_max = Soln_ptr[i_blk].ICu;
             i_inc = 1;
             j_min = Soln_ptr[i_blk].JCu+1;
             j_max = Soln_ptr[i_blk].JCu+Soln_Block_List.Block[i_blk].info.dimen.ghost;
             j_inc = 1;
             i = Soln_ptr[i_blk].UnloadReceiveBuffer(Soln_Block_List.message_noreschange_northface_recbuf[i_blk],
                                                     l,buffer_size,
                                                     i_min,i_max,i_inc,
                                                     j_min,j_max,j_inc);
	     if (i != 0) return(2200);
          } /* endif */
          // Unload ghost cell mesh information as required.
          if (Send_Mesh_Geometry_Only ||
              Number_of_Solution_Variables > Soln_ptr[i_blk].NumVar()) {
	     i_min = Soln_ptr[i_blk].Grid.INl;
	     i_max = Soln_ptr[i_blk].Grid.INu;
	     i_inc = 1;
	     j_min = Soln_ptr[i_blk].Grid.JNu+1;
	     j_max = Soln_ptr[i_blk].Grid.JNu+Soln_Block_List.Block[i_blk].info.dimen.ghost;
	     j_inc = 1;
             x_ref = Soln_ptr[i_blk].Grid.Node[i_min][j_min-1].X; // Reference node location.
             for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc ) {
                for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
		   l = l + NUM_COMP_VECTOR2D;
                   if (l >= buffer_size) return(2201);
                   Soln_ptr[i_blk].Grid.Node[i][j].X = 
		     Vector2D(Soln_Block_List.message_noreschange_northface_recbuf[i_blk][l-1],
                              Soln_Block_List.message_noreschange_northface_recbuf[i_blk][l])+x_ref;
                } /* endfor */
             } /* endfor */
             i_min = Soln_ptr[i_blk].ICl;
	     i_max = Soln_ptr[i_blk].ICu;
	     i_inc = 1;
	     j_min = Soln_ptr[i_blk].JCu+1;
	     j_max = Soln_ptr[i_blk].JCu+Soln_Block_List.Block[i_blk].info.dimen.ghost;
	     j_inc = 1;
             for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc ) {
                for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
	           Soln_ptr[i_blk].Grid.Cell[i][j].Xc = Soln_ptr[i_blk].Grid.centroid(i, j);
	           Soln_ptr[i_blk].Grid.Cell[i][j].A = Soln_ptr[i_blk].Grid.area(i, j);
                } /* endfor */
             } /* endfor */
             for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc ) {
                l = l + 1;
                if (l >= buffer_size) return(2202);
                Soln_ptr[i_blk].Grid.BCtypeW[j] = 
                   int(Soln_Block_List.message_noreschange_northface_recbuf[i_blk][l]);
             } /* endfor */
             for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc ) {
                l = l + 1;
                if (l >= buffer_size) return(2203);
                Soln_ptr[i_blk].Grid.BCtypeE[j] = 
                   int(Soln_Block_List.message_noreschange_northface_recbuf[i_blk][l]);
             } /* endfor */
          } /* endif */
       } /* endif */

       /////////////////////////////////////////////////////////////////////////
       // Unload receive buffer for South face neighbours of solution blocks. //
       /////////////////////////////////////////////////////////////////////////
       if (Soln_Block_List.Block[i_blk].used && 
           (Soln_Block_List.Block[i_blk].nS == 1) &&
           (Soln_Block_List.Block[i_blk].info.level == 
            Soln_Block_List.Block[i_blk].infoS[0].level)) {
          if (!Send_Mesh_Geometry_Only) {
             buffer_size = Soln_Block_List.Block[i_blk].info.dimen.ghost*
                           abs(Soln_Block_List.Block[i_blk].info.dimen.i)*
                           Soln_ptr[i_blk].NumVar();
             if (Number_of_Solution_Variables > Soln_ptr[i_blk].NumVar()) 
                buffer_size = buffer_size +
                              Soln_Block_List.Block[i_blk].info.dimen.ghost*
                              (abs(Soln_Block_List.Block[i_blk].info.dimen.i)+1)*
                              NUM_COMP_VECTOR2D+
                              2*Soln_Block_List.Block[i_blk].info.dimen.ghost;
          } else {
                buffer_size = Soln_Block_List.Block[i_blk].info.dimen.ghost*
                              (abs(Soln_Block_List.Block[i_blk].info.dimen.i)+1)*
                              NUM_COMP_VECTOR2D+
                              2*Soln_Block_List.Block[i_blk].info.dimen.ghost;
          } /* endif */
          l = -1;
          // Unload ghost cell solution information as required.
          if (!Send_Mesh_Geometry_Only) {
   	     i_min = Soln_ptr[i_blk].ICl;
	     i_max = Soln_ptr[i_blk].ICu;
	     i_inc = 1;
	     j_min = Soln_ptr[i_blk].JCl-Soln_Block_List.Block[i_blk].info.dimen.ghost;
	     j_max = Soln_ptr[i_blk].JCl-1;
	     j_inc = 1;
             i = Soln_ptr[i_blk].UnloadReceiveBuffer(Soln_Block_List.message_noreschange_southface_recbuf[i_blk],
                                                     l,buffer_size,
                                                     i_min,i_max,i_inc,
                                                     j_min,j_max,j_inc);
	     if (i != 0) return(2100);
          } /* endif */
          // Unload ghost cell mesh information as required.
          if (Send_Mesh_Geometry_Only ||
              Number_of_Solution_Variables > Soln_ptr[i_blk].NumVar()) {
	     i_min = Soln_ptr[i_blk].Grid.INl;
	     i_max = Soln_ptr[i_blk].Grid.INu;
	     i_inc = 1;
	     j_min = Soln_ptr[i_blk].Grid.JNl-Soln_Block_List.Block[i_blk].info.dimen.ghost;
	     j_max = Soln_ptr[i_blk].Grid.JNl-1;
	     j_inc = 1;
             x_ref = Soln_ptr[i_blk].Grid.Node[i_min][j_max+1].X; // Reference node location.
             for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc ) {
                for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
		   l = l + NUM_COMP_VECTOR2D;
                   if (l >= buffer_size) return(2101);
                   Soln_ptr[i_blk].Grid.Node[i][j].X = 
		     Vector2D(Soln_Block_List.message_noreschange_southface_recbuf[i_blk][l-1],
                              Soln_Block_List.message_noreschange_southface_recbuf[i_blk][l])+x_ref;
                } /* endfor */
             } /* endfor */
  	     i_min = Soln_ptr[i_blk].ICl;
	     i_max = Soln_ptr[i_blk].ICu;
	     i_inc = 1;
	     j_min = Soln_ptr[i_blk].JCl-Soln_Block_List.Block[i_blk].info.dimen.ghost;
	     j_max = Soln_ptr[i_blk].JCl-1;
	     j_inc = 1;
             for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc ) {
                for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
	           Soln_ptr[i_blk].Grid.Cell[i][j].Xc = Soln_ptr[i_blk].Grid.centroid(i, j);
	           Soln_ptr[i_blk].Grid.Cell[i][j].A = Soln_ptr[i_blk].Grid.area(i, j);
                } /* endfor */
             } /* endfor */
             for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc ) {
                l = l + 1;
                if (l >= buffer_size) return(2102);
                Soln_ptr[i_blk].Grid.BCtypeW[j] = 
                   int(Soln_Block_List.message_noreschange_southface_recbuf[i_blk][l]);
             } /* endfor */
             for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc ) {
                l = l + 1;
                if (l >= buffer_size) return(2103);
                Soln_ptr[i_blk].Grid.BCtypeE[j] = 
                   int(Soln_Block_List.message_noreschange_southface_recbuf[i_blk][l]);
             } /* endfor */
          } /* endif */
       } /* endif */

       ////////////////////////////////////////////////////////////////////////
       // Unload receive buffer for East face neighbours of solution blocks. //
       ////////////////////////////////////////////////////////////////////////
       if (Soln_Block_List.Block[i_blk].used && 
           (Soln_Block_List.Block[i_blk].nE == 1) &&
           (Soln_Block_List.Block[i_blk].info.level == 
            Soln_Block_List.Block[i_blk].infoE[0].level)) {
          if (!Send_Mesh_Geometry_Only) {
             buffer_size = Soln_Block_List.Block[i_blk].info.dimen.ghost*
                           abs(Soln_Block_List.Block[i_blk].info.dimen.j)*
                           Soln_ptr[i_blk].NumVar();
             if (Number_of_Solution_Variables > Soln_ptr[i_blk].NumVar()) 
                buffer_size = buffer_size +
                              Soln_Block_List.Block[i_blk].info.dimen.ghost*
                              (abs(Soln_Block_List.Block[i_blk].info.dimen.j)+1)*
                              NUM_COMP_VECTOR2D+
                              2*Soln_Block_List.Block[i_blk].info.dimen.ghost;
          } else {
                buffer_size = Soln_Block_List.Block[i_blk].info.dimen.ghost*
                              (abs(Soln_Block_List.Block[i_blk].info.dimen.j)+1)*
                              NUM_COMP_VECTOR2D+
                              2*Soln_Block_List.Block[i_blk].info.dimen.ghost;
          } /* endif */
          l = -1;
          // Unload ghost cell solution information as required.
          if (!Send_Mesh_Geometry_Only) {
   	     i_min = Soln_ptr[i_blk].ICu+1;
 	     i_max = Soln_ptr[i_blk].ICu+Soln_Block_List.Block[i_blk].info.dimen.ghost;
	     i_inc = 1;
	     j_min = Soln_ptr[i_blk].JCl;
	     j_max = Soln_ptr[i_blk].JCu;
	     j_inc = 1;
             i = Soln_ptr[i_blk].UnloadReceiveBuffer(Soln_Block_List.message_noreschange_eastface_recbuf[i_blk],
                                                     l,buffer_size,
                                                     i_min,i_max,i_inc,
                                                     j_min,j_max,j_inc);
	     if (i != 0) return(1200);
          } /* endif */
          // Unload ghost cell mesh information as required.
          if (Send_Mesh_Geometry_Only ||
              Number_of_Solution_Variables > Soln_ptr[i_blk].NumVar()) {
	     i_min = Soln_ptr[i_blk].Grid.INu+1;
	     i_max = Soln_ptr[i_blk].Grid.INu+Soln_Block_List.Block[i_blk].info.dimen.ghost;
	     i_inc = 1;
	     j_min = Soln_ptr[i_blk].Grid.JNl;
	     j_max = Soln_ptr[i_blk].Grid.JNu;
	     j_inc = 1;
             x_ref = Soln_ptr[i_blk].Grid.Node[i_min-1][j_min].X; // Reference node location.
             for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc ) {
                for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
		   l = l + NUM_COMP_VECTOR2D;
                   if (l >= buffer_size) return(1201);
                   Soln_ptr[i_blk].Grid.Node[i][j].X = 
		     Vector2D(Soln_Block_List.message_noreschange_eastface_recbuf[i_blk][l-1],
                              Soln_Block_List.message_noreschange_eastface_recbuf[i_blk][l])+x_ref;
                } /* endfor */
             } /* endfor */
 	     i_min = Soln_ptr[i_blk].ICu+1;
	     i_max = Soln_ptr[i_blk].ICu+Soln_Block_List.Block[i_blk].info.dimen.ghost;
	     i_inc = 1;
	     j_min = Soln_ptr[i_blk].JCl;
	     j_max = Soln_ptr[i_blk].JCu;
	     j_inc = 1;
             for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc ) {
                for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
	           Soln_ptr[i_blk].Grid.Cell[i][j].Xc = Soln_ptr[i_blk].Grid.centroid(i, j);
	           Soln_ptr[i_blk].Grid.Cell[i][j].A = Soln_ptr[i_blk].Grid.area(i, j);
                } /* endfor */
             } /* endfor */
             for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
                l = l + 1;
                if (l >= buffer_size) return(1202);
                Soln_ptr[i_blk].Grid.BCtypeS[i] = 
                   int(Soln_Block_List.message_noreschange_eastface_recbuf[i_blk][l]);
             } /* endfor */
             for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
                l = l + 1;
                if (l >= buffer_size) return(1203);
                Soln_ptr[i_blk].Grid.BCtypeN[i] = 
                   int(Soln_Block_List.message_noreschange_eastface_recbuf[i_blk][l]);
             } /* endfor */
          } /* endif */
       } /* endif */

       ////////////////////////////////////////////////////////////////////////
       // Unload receive buffer for West face neighbours of solution blocks. //
       ////////////////////////////////////////////////////////////////////////
       if (Soln_Block_List.Block[i_blk].used && 
           (Soln_Block_List.Block[i_blk].nW == 1) &&
           (Soln_Block_List.Block[i_blk].info.level == 
            Soln_Block_List.Block[i_blk].infoW[0].level)) {
          if (!Send_Mesh_Geometry_Only) {
             buffer_size = Soln_Block_List.Block[i_blk].info.dimen.ghost*
                           abs(Soln_Block_List.Block[i_blk].info.dimen.j)*
                           Soln_ptr[i_blk].NumVar();
             if (Number_of_Solution_Variables > Soln_ptr[i_blk].NumVar()) 
                buffer_size = buffer_size +
                              Soln_Block_List.Block[i_blk].info.dimen.ghost*
                              (abs(Soln_Block_List.Block[i_blk].info.dimen.j)+1)*
                              NUM_COMP_VECTOR2D+
                              2*Soln_Block_List.Block[i_blk].info.dimen.ghost;
          } else {
                buffer_size = Soln_Block_List.Block[i_blk].info.dimen.ghost*
                              (abs(Soln_Block_List.Block[i_blk].info.dimen.j)+1)*
                              NUM_COMP_VECTOR2D+
                              2*Soln_Block_List.Block[i_blk].info.dimen.ghost;
          } /* endif */
          l = -1;
          // Unload ghost cell solution information as required.
          if (!Send_Mesh_Geometry_Only) {
  	     i_min = Soln_ptr[i_blk].ICl-Soln_Block_List.Block[i_blk].info.dimen.ghost;
	     i_max = Soln_ptr[i_blk].ICl-1;
	     i_inc = 1;
	     j_min = Soln_ptr[i_blk].JCl;
	     j_max = Soln_ptr[i_blk].JCu;
	     j_inc = 1;
             i = Soln_ptr[i_blk].UnloadReceiveBuffer(Soln_Block_List.message_noreschange_westface_recbuf[i_blk],
                                                     l,buffer_size,
                                                     i_min,i_max,i_inc,
                                                     j_min,j_max,j_inc);
	     if (i != 0) return(1100);
          } /* endif */
          // Unload ghost cell mesh information as required.
          if (Send_Mesh_Geometry_Only ||
              Number_of_Solution_Variables > Soln_ptr[i_blk].NumVar()) {
	     i_min = Soln_ptr[i_blk].Grid.INl-Soln_Block_List.Block[i_blk].info.dimen.ghost;
	     i_max = Soln_ptr[i_blk].Grid.INl-1;
	     i_inc = 1;
	     j_min = Soln_ptr[i_blk].Grid.JNl;
	     j_max = Soln_ptr[i_blk].Grid.JNu;
	     j_inc = 1;
             x_ref = Soln_ptr[i_blk].Grid.Node[i_max+1][j_min].X; // Reference node location.
             for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc ) {
                for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
		   l = l + NUM_COMP_VECTOR2D;
                   if (l >= buffer_size) return(1101);
                   Soln_ptr[i_blk].Grid.Node[i][j].X = 
		     Vector2D(Soln_Block_List.message_noreschange_westface_recbuf[i_blk][l-1],
                              Soln_Block_List.message_noreschange_westface_recbuf[i_blk][l])+x_ref;
                } /* endfor */
             } /* endfor */
 	     i_min = Soln_ptr[i_blk].ICl-Soln_Block_List.Block[i_blk].info.dimen.ghost;
	     i_max = Soln_ptr[i_blk].ICl-1;
	     i_inc = 1;
	     j_min = Soln_ptr[i_blk].JCl;
	     j_max = Soln_ptr[i_blk].JCu;
	     j_inc = 1;
             for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc ) {
                for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
	           Soln_ptr[i_blk].Grid.Cell[i][j].Xc = Soln_ptr[i_blk].Grid.centroid(i, j);
	           Soln_ptr[i_blk].Grid.Cell[i][j].A = Soln_ptr[i_blk].Grid.area(i, j);
                } /* endfor */
             } /* endfor */
             for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
                l = l + 1;
                if (l >= buffer_size) return(1102);
                Soln_ptr[i_blk].Grid.BCtypeS[i] = 
                   int(Soln_Block_List.message_noreschange_westface_recbuf[i_blk][l]);
             } /* endfor */
             for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
                l = l + 1;
                if (l >= buffer_size) return(1103);
                Soln_ptr[i_blk].Grid.BCtypeN[i] = 
                   int(Soln_Block_List.message_noreschange_westface_recbuf[i_blk][l]);
             } /* endfor */
          } /* endif */
       } /* endif */

       ////////////////////////////////////////////////////////////////////////////////
       // Unload receive buffer for North West corner neighbours of solution blocks. //
       ////////////////////////////////////////////////////////////////////////////////
       if (Soln_Block_List.Block[i_blk].used && 
           (Soln_Block_List.Block[i_blk].nNW == 1) &&
           (Soln_Block_List.Block[i_blk].info.level == 
            Soln_Block_List.Block[i_blk].infoNW[0].level)) {
          buffer_size = sqr(Soln_Block_List.Block[i_blk].info.dimen.ghost)*
                        Number_of_Solution_Variables;
          l = -1;
          // Unload ghost cell solution information as required.
          if (!Send_Mesh_Geometry_Only) {
             i_min = Soln_ptr[i_blk].ICl-Soln_Block_List.Block[i_blk].info.dimen.ghost;
	     i_max = Soln_ptr[i_blk].ICl-1;
             i_inc = 1;
             j_min = Soln_ptr[i_blk].JCu+1;
             j_max = Soln_ptr[i_blk].JCu+Soln_Block_List.Block[i_blk].info.dimen.ghost;
             j_inc = 1;
             i = Soln_ptr[i_blk].UnloadReceiveBuffer(Soln_Block_List.message_noreschange_northwestcorner_recbuf[i_blk],
                                                     l,buffer_size,
                                                     i_min,i_max,i_inc,
                                                     j_min,j_max,j_inc);
	     if (i != 0) return(4200);
          } /* endif */
          // Unload ghost cell mesh information as required.
          if (Send_Mesh_Geometry_Only ||
              Number_of_Solution_Variables > Soln_ptr[i_blk].NumVar()) {
	     i_min = Soln_ptr[i_blk].Grid.INl-Soln_Block_List.Block[i_blk].info.dimen.ghost;
	     i_max = Soln_ptr[i_blk].Grid.INl-1;
	     i_inc = 1;
	     j_min = Soln_ptr[i_blk].Grid.JNu+1;
	     j_max = Soln_ptr[i_blk].Grid.JNu+Soln_Block_List.Block[i_blk].info.dimen.ghost;
	     j_inc = 1;
             x_ref = Soln_ptr[i_blk].Grid.Node[i_max+1][j_min-1].X; // Reference node location.
             for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc ) {
                for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
		   l = l + NUM_COMP_VECTOR2D;
                   if (l >= buffer_size) return(4201);
                   Soln_ptr[i_blk].Grid.Node[i][j].X = 
		     Vector2D(Soln_Block_List.message_noreschange_northwestcorner_recbuf[i_blk][l-1],
                              Soln_Block_List.message_noreschange_northwestcorner_recbuf[i_blk][l])+x_ref;
                } /* endfor */
             } /* endfor */
             i_min = Soln_ptr[i_blk].ICl-Soln_Block_List.Block[i_blk].info.dimen.ghost;
	     i_max = Soln_ptr[i_blk].ICl-1;
	     i_inc = 1;
	     j_min = Soln_ptr[i_blk].JCu+1;
	     j_max = Soln_ptr[i_blk].JCu+Soln_Block_List.Block[i_blk].info.dimen.ghost;
	     j_inc = 1;
             for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc ) {
                for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
	           Soln_ptr[i_blk].Grid.Cell[i][j].Xc = Soln_ptr[i_blk].Grid.centroid(i, j);
	           Soln_ptr[i_blk].Grid.Cell[i][j].A = Soln_ptr[i_blk].Grid.area(i, j);
                } /* endfor */
             } /* endfor */
          } /* endif */
       } /* endif */

       ////////////////////////////////////////////////////////////////////////////////
       // Unload receive buffer for North East corner neighbours of solution blocks. //
       ////////////////////////////////////////////////////////////////////////////////
       if (Soln_Block_List.Block[i_blk].used && 
           (Soln_Block_List.Block[i_blk].nNE == 1) &&
           (Soln_Block_List.Block[i_blk].info.level == 
            Soln_Block_List.Block[i_blk].infoNE[0].level)) {
          buffer_size = sqr(Soln_Block_List.Block[i_blk].info.dimen.ghost)*
                        Number_of_Solution_Variables;
          l = -1;
          // Unload ghost cell solution information as required.
          if (!Send_Mesh_Geometry_Only) {
             i_min = Soln_ptr[i_blk].ICu+1;
	     i_max = Soln_ptr[i_blk].ICu+Soln_Block_List.Block[i_blk].info.dimen.ghost;
             i_inc = 1;
             j_min = Soln_ptr[i_blk].JCu+1;
             j_max = Soln_ptr[i_blk].JCu+Soln_Block_List.Block[i_blk].info.dimen.ghost;
             j_inc = 1;
             i = Soln_ptr[i_blk].UnloadReceiveBuffer(Soln_Block_List.message_noreschange_northeastcorner_recbuf[i_blk],
                                                     l,buffer_size,
                                                     i_min,i_max,i_inc,
                                                     j_min,j_max,j_inc);
	     if (i != 0) return(5200);
          } /* endif */
          // Unload ghost cell mesh information as required.
          if (Send_Mesh_Geometry_Only ||
              Number_of_Solution_Variables > Soln_ptr[i_blk].NumVar()) {
	     i_min = Soln_ptr[i_blk].Grid.INu+1;
	     i_max = Soln_ptr[i_blk].Grid.INu+Soln_Block_List.Block[i_blk].info.dimen.ghost;
	     i_inc = 1;
	     j_min = Soln_ptr[i_blk].Grid.JNu+1;
	     j_max = Soln_ptr[i_blk].Grid.JNu+Soln_Block_List.Block[i_blk].info.dimen.ghost;
	     j_inc = 1;
             x_ref = Soln_ptr[i_blk].Grid.Node[i_min-1][j_min-1].X; // Reference node location.
             for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc ) {
                for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
		   l = l + NUM_COMP_VECTOR2D;
                   if (l >= buffer_size) return(5201);
                   Soln_ptr[i_blk].Grid.Node[i][j].X = 
		     Vector2D(Soln_Block_List.message_noreschange_northeastcorner_recbuf[i_blk][l-1],
                              Soln_Block_List.message_noreschange_northeastcorner_recbuf[i_blk][l])+x_ref;
                } /* endfor */
             } /* endfor */
             i_min = Soln_ptr[i_blk].ICu+1;
	     i_max = Soln_ptr[i_blk].ICu+Soln_Block_List.Block[i_blk].info.dimen.ghost;
	     i_inc = 1;
	     j_min = Soln_ptr[i_blk].JCu+1;
	     j_max = Soln_ptr[i_blk].JCu+Soln_Block_List.Block[i_blk].info.dimen.ghost;
	     j_inc = 1;
             for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc ) {
                for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
	           Soln_ptr[i_blk].Grid.Cell[i][j].Xc = Soln_ptr[i_blk].Grid.centroid(i, j);
	           Soln_ptr[i_blk].Grid.Cell[i][j].A = Soln_ptr[i_blk].Grid.area(i, j);
                } /* endfor */
             } /* endfor */
          } /* endif */
       } /* endif */

       ////////////////////////////////////////////////////////////////////////////////
       // Unload receive buffer for South East corner neighbours of solution blocks. //
       ////////////////////////////////////////////////////////////////////////////////
       if (Soln_Block_List.Block[i_blk].used && 
           (Soln_Block_List.Block[i_blk].nSE == 1) &&
           (Soln_Block_List.Block[i_blk].info.level == 
            Soln_Block_List.Block[i_blk].infoSE[0].level)) {
          buffer_size = sqr(Soln_Block_List.Block[i_blk].info.dimen.ghost)*
                        Number_of_Solution_Variables;
          l = -1;
          // Unload ghost cell solution information as required.
          if (!Send_Mesh_Geometry_Only) {
             i_min = Soln_ptr[i_blk].ICu+1;
	     i_max = Soln_ptr[i_blk].ICu+Soln_Block_List.Block[i_blk].info.dimen.ghost;
             i_inc = 1;
             j_min = Soln_ptr[i_blk].JCl-Soln_Block_List.Block[i_blk].info.dimen.ghost;
             j_max = Soln_ptr[i_blk].JCl-1;
             j_inc = 1;
             i = Soln_ptr[i_blk].UnloadReceiveBuffer(Soln_Block_List.message_noreschange_southeastcorner_recbuf[i_blk],
                                                     l,buffer_size,
                                                     i_min,i_max,i_inc,
                                                     j_min,j_max,j_inc);
	     if (i != 0) return(4100);
          } /* endif */
          // Unload ghost cell mesh information as required.
          if (Send_Mesh_Geometry_Only ||
              Number_of_Solution_Variables > Soln_ptr[i_blk].NumVar()) {
	     i_min = Soln_ptr[i_blk].Grid.INu+1;
	     i_max = Soln_ptr[i_blk].Grid.INu+Soln_Block_List.Block[i_blk].info.dimen.ghost;
	     i_inc = 1;
	     j_min = Soln_ptr[i_blk].Grid.JNl-Soln_Block_List.Block[i_blk].info.dimen.ghost;
	     j_max = Soln_ptr[i_blk].Grid.JNl-1;
	     j_inc = 1;
             x_ref = Soln_ptr[i_blk].Grid.Node[i_min-1][j_max+1].X; // Reference node location.
             for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc ) {
                for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
		   l = l + NUM_COMP_VECTOR2D;
                   if (l >= buffer_size) return(4101);
                   Soln_ptr[i_blk].Grid.Node[i][j].X = 
		     Vector2D(Soln_Block_List.message_noreschange_southeastcorner_recbuf[i_blk][l-1],
                              Soln_Block_List.message_noreschange_southeastcorner_recbuf[i_blk][l])+x_ref;
                } /* endfor */
             } /* endfor */
             i_min = Soln_ptr[i_blk].ICu+1;
	     i_max = Soln_ptr[i_blk].ICu+Soln_Block_List.Block[i_blk].info.dimen.ghost;
	     i_inc = 1;
             j_min = Soln_ptr[i_blk].JCl-Soln_Block_List.Block[i_blk].info.dimen.ghost;
             j_max = Soln_ptr[i_blk].JCl-1;
	     j_inc = 1;
             for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc ) {
                for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
	           Soln_ptr[i_blk].Grid.Cell[i][j].Xc = Soln_ptr[i_blk].Grid.centroid(i, j);
	           Soln_ptr[i_blk].Grid.Cell[i][j].A = Soln_ptr[i_blk].Grid.area(i, j);
                } /* endfor */
             } /* endfor */
          } /* endif */
       } /* endif */

       ////////////////////////////////////////////////////////////////////////////////
       // Unload receive buffer for South West corner neighbours of solution blocks. //
       ////////////////////////////////////////////////////////////////////////////////
       if (Soln_Block_List.Block[i_blk].used && 
           (Soln_Block_List.Block[i_blk].nSW == 1) &&
           (Soln_Block_List.Block[i_blk].info.level == 
            Soln_Block_List.Block[i_blk].infoSW[0].level)) {
          buffer_size = sqr(Soln_Block_List.Block[i_blk].info.dimen.ghost)*
                        Number_of_Solution_Variables;
          l = -1;
          // Unload ghost cell solution information as required.
          if (!Send_Mesh_Geometry_Only) {
             i_min = Soln_ptr[i_blk].ICl-Soln_Block_List.Block[i_blk].info.dimen.ghost;
	     i_max = Soln_ptr[i_blk].ICl-1;
             i_inc = 1;
             j_min = Soln_ptr[i_blk].JCl-Soln_Block_List.Block[i_blk].info.dimen.ghost;
             j_max = Soln_ptr[i_blk].JCl-1;
             j_inc = 1;
             i = Soln_ptr[i_blk].UnloadReceiveBuffer(Soln_Block_List.message_noreschange_southwestcorner_recbuf[i_blk],
                                                     l,buffer_size,
                                                     i_min,i_max,i_inc,
                                                     j_min,j_max,j_inc);
	     if (i != 0) return(5100);
          } /* endif */
          // Unload ghost cell mesh information as required.
          if (Send_Mesh_Geometry_Only ||
              Number_of_Solution_Variables > Soln_ptr[i_blk].NumVar()) {
	     i_min = Soln_ptr[i_blk].Grid.INl-Soln_Block_List.Block[i_blk].info.dimen.ghost;
	     i_max = Soln_ptr[i_blk].Grid.INl-1;
	     i_inc = 1;
	     j_min = Soln_ptr[i_blk].Grid.JNl-Soln_Block_List.Block[i_blk].info.dimen.ghost;
	     j_max = Soln_ptr[i_blk].Grid.JNl-1;
	     j_inc = 1;
             x_ref = Soln_ptr[i_blk].Grid.Node[i_max+1][j_max+1].X; // Reference node location.
             for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc ) {
                for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
		   l = l + NUM_COMP_VECTOR2D;
                   if (l >= buffer_size) return(5101);
                   Soln_ptr[i_blk].Grid.Node[i][j].X = 
		     Vector2D(Soln_Block_List.message_noreschange_southwestcorner_recbuf[i_blk][l-1],
                              Soln_Block_List.message_noreschange_southwestcorner_recbuf[i_blk][l])+x_ref;
                } /* endfor */
             } /* endfor */
             i_min = Soln_ptr[i_blk].ICl-Soln_Block_List.Block[i_blk].info.dimen.ghost;
	     i_max = Soln_ptr[i_blk].ICl-1;
	     i_inc = 1;
             j_min = Soln_ptr[i_blk].JCl-Soln_Block_List.Block[i_blk].info.dimen.ghost;
             j_max = Soln_ptr[i_blk].JCl-1;
	     j_inc = 1;
             for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc ) {
                for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
	           Soln_ptr[i_blk].Grid.Cell[i][j].Xc = Soln_ptr[i_blk].Grid.centroid(i, j);
	           Soln_ptr[i_blk].Grid.Cell[i][j].A = Soln_ptr[i_blk].Grid.area(i, j);
                } /* endfor */
             } /* endfor */
          } /* endif */
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
template <class Quad_Soln_Block>
int Unload_Receive_Message_Buffers_ResChange_FineToCoarse(Quad_Soln_Block *Soln_ptr,
                                                          AdaptiveBlock2D_List &Soln_Block_List,
                                                          const int Number_of_Solution_Variables,
                                                          const int Send_Mesh_Geometry_Only,
                                                          const int Send_Conservative_Solution_Fluxes) {

    int i_blk, buffer_size, 
        i_min, i_max, i_inc, i, 
        j_min, j_max, j_inc, j, 
        l, j_neigh;
    Vector2D x_ref;

    /* Check to see if the number of solution variables specified
       in the calling argument is correct. */

    if (Number_of_Solution_Variables != NUM_COMP_VECTOR2D &&
        Number_of_Solution_Variables != Soln_ptr[0].NumVar() &&
        Number_of_Solution_Variables != Soln_ptr[0].NumVar() + NUM_COMP_VECTOR2D) {
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
          if (!Send_Mesh_Geometry_Only) {
             if (!Send_Conservative_Solution_Fluxes) {
                buffer_size = Soln_Block_List.Block[i_blk].info.dimen.ghost*
                              (abs(Soln_Block_List.Block[i_blk].info.dimen.i)/2)*
                              Soln_ptr[i_blk].NumVar();
             } else {
                buffer_size = (abs(Soln_Block_List.Block[i_blk].info.dimen.i)/2)*
                              Soln_ptr[i_blk].NumVar();
             } /* endif */
             if (Number_of_Solution_Variables > Soln_ptr[i_blk].NumVar()) 
                buffer_size = buffer_size +
                              Soln_Block_List.Block[i_blk].info.dimen.ghost*
                              (abs(Soln_Block_List.Block[i_blk].info.dimen.i)/2+1)*
                              NUM_COMP_VECTOR2D+
                              2*Soln_Block_List.Block[i_blk].info.dimen.ghost;
          } else {
                buffer_size = Soln_Block_List.Block[i_blk].info.dimen.ghost*
                              (abs(Soln_Block_List.Block[i_blk].info.dimen.i)/2+1)*
                              NUM_COMP_VECTOR2D+
                              2*Soln_Block_List.Block[i_blk].info.dimen.ghost;
          } /* endif */
          for (j_neigh = 0; j_neigh <= Soln_Block_List.Block[i_blk].nN-1 ; ++j_neigh) {
             l = -1;
             // Unload ghost cell solution information as required.
             if (!Send_Mesh_Geometry_Only) {
                if (!Send_Conservative_Solution_Fluxes) {
                   if (j_neigh == 0) {
                      i_min = Soln_ptr[i_blk].ICl;
	              i_max = Soln_ptr[i_blk].ICl+(Soln_ptr[i_blk].ICu-Soln_ptr[i_blk].ICl-1)/2;
	              i_inc = 1;
	              j_min = Soln_ptr[i_blk].JCu+1;
	              j_max = Soln_ptr[i_blk].JCu+Soln_Block_List.Block[i_blk].info.dimen.ghost;
	              j_inc = 1;
                   } else {
                      i_min = Soln_ptr[i_blk].ICl+(Soln_ptr[i_blk].ICu-Soln_ptr[i_blk].ICl+1)/2;
	              i_max = Soln_ptr[i_blk].ICu;
	              i_inc = 1;
	              j_min = Soln_ptr[i_blk].JCu+1;
	              j_max = Soln_ptr[i_blk].JCu+Soln_Block_List.Block[i_blk].info.dimen.ghost;
	              j_inc = 1;
                   } /* endif */
                   i = Soln_ptr[i_blk].UnloadReceiveBuffer_F2C(Soln_Block_List.message_reschange_northface_recbuf[i_blk][j_neigh],
                                                               l,buffer_size,
                                                               i_min,i_max,i_inc,
                                                               j_min,j_max,j_inc);
	           if (i != 0) return(2200);
                } else {
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
                   i = Soln_ptr[i_blk].UnloadReceiveBuffer_Flux_F2C(Soln_Block_List.message_reschange_northface_recbuf[i_blk][j_neigh],
                                                                    l,buffer_size,
                                                                    i_min,i_max,i_inc,
                                                                    j_min,j_max,j_inc);
	           if (i != 0) return(2201);
                } /* endif */
             } /* endif */
             // Unload ghost cell mesh information as required.
             if (Send_Mesh_Geometry_Only ||
                 Number_of_Solution_Variables > Soln_ptr[i_blk].NumVar()) {
                if (j_neigh == 0) {
                   i_min = Soln_ptr[i_blk].Grid.INl;
	           i_max = Soln_ptr[i_blk].Grid.INl+(Soln_ptr[i_blk].Grid.INu-Soln_ptr[i_blk].Grid.INl)/2;
	           i_inc = 1;
	           j_min = Soln_ptr[i_blk].Grid.JNu+1;
	           j_max = Soln_ptr[i_blk].Grid.JNu+Soln_Block_List.Block[i_blk].info.dimen.ghost;
	           j_inc = 1;
                } else {
                   i_min = Soln_ptr[i_blk].Grid.INl+(Soln_ptr[i_blk].Grid.INu-Soln_ptr[i_blk].Grid.INl)/2;
	           i_max = Soln_ptr[i_blk].Grid.INu;
	           i_inc = 1;
	           j_min = Soln_ptr[i_blk].Grid.JNu+1;
	           j_max = Soln_ptr[i_blk].Grid.JNu+Soln_Block_List.Block[i_blk].info.dimen.ghost;
	           j_inc = 1;
                } /* endif */
                x_ref = Soln_ptr[i_blk].Grid.Node[i_min][j_min-1].X; // Reference node location.
                for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc ) {
                   for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
		      l = l + NUM_COMP_VECTOR2D;
                      if (l >= buffer_size) return(2202);
                      Soln_ptr[i_blk].Grid.Node[i][j].X = 
		        Vector2D(Soln_Block_List.message_reschange_northface_recbuf[i_blk][j_neigh][l-1],
                                 Soln_Block_List.message_reschange_northface_recbuf[i_blk][j_neigh][l])+x_ref;
                   } /* endfor */
                } /* endfor */
                if (j_neigh == 0) {
                   i_min = Soln_ptr[i_blk].ICl;
	           i_max = Soln_ptr[i_blk].ICl+(Soln_ptr[i_blk].ICu-Soln_ptr[i_blk].ICl-1)/2;
	           i_inc = 1;
	           j_min = Soln_ptr[i_blk].JCu+1;
	           j_max = Soln_ptr[i_blk].JCu+Soln_Block_List.Block[i_blk].info.dimen.ghost;
	           j_inc = 1;
                } else {
                   i_min = Soln_ptr[i_blk].ICl+(Soln_ptr[i_blk].ICu-Soln_ptr[i_blk].ICl+1)/2;
	           i_max = Soln_ptr[i_blk].ICu;
	           i_inc = 1;
	           j_min = Soln_ptr[i_blk].JCu+1;
	           j_max = Soln_ptr[i_blk].JCu+Soln_Block_List.Block[i_blk].info.dimen.ghost;
	           j_inc = 1;
                } /* endif */
                for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc ) {
                   for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
	              Soln_ptr[i_blk].Grid.Cell[i][j].Xc = Soln_ptr[i_blk].Grid.centroid(i, j);
	              Soln_ptr[i_blk].Grid.Cell[i][j].A = Soln_ptr[i_blk].Grid.area(i, j);
                   } /* endfor */
                } /* endfor */
                if (j_neigh == 0) {
                   for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc ) {
                      l = l + 1;
                      if (l >= buffer_size) return(2203);
                      Soln_ptr[i_blk].Grid.BCtypeW[j] = 
                         int(Soln_Block_List.message_reschange_northface_recbuf[i_blk][j_neigh][l]);
                   } /* endfor */
                } else {
                   l = l + Soln_Block_List.Block[i_blk].info.dimen.ghost;
                   for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc ) {
                      l = l + 1;
                      if (l >= buffer_size) return(2204);
                      Soln_ptr[i_blk].Grid.BCtypeE[j] = 
                         int(Soln_Block_List.message_reschange_northface_recbuf[i_blk][j_neigh][l]);
                   } /* endfor */
                } /* endif */
             } /* endif */
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
          if (!Send_Mesh_Geometry_Only) {
             if (!Send_Conservative_Solution_Fluxes) {
                buffer_size = Soln_Block_List.Block[i_blk].info.dimen.ghost*
                              (abs(Soln_Block_List.Block[i_blk].info.dimen.i)/2)*
                              Soln_ptr[i_blk].NumVar();
             } else {
                buffer_size = (abs(Soln_Block_List.Block[i_blk].info.dimen.i)/2)*
                              Soln_ptr[i_blk].NumVar();
             } /* endif */
             if (Number_of_Solution_Variables > Soln_ptr[i_blk].NumVar()) 
                buffer_size = buffer_size +
                              Soln_Block_List.Block[i_blk].info.dimen.ghost*
                              (abs(Soln_Block_List.Block[i_blk].info.dimen.i)/2+1)*
                              NUM_COMP_VECTOR2D+
                              2*Soln_Block_List.Block[i_blk].info.dimen.ghost;
          } else {
                buffer_size = Soln_Block_List.Block[i_blk].info.dimen.ghost*
                              (abs(Soln_Block_List.Block[i_blk].info.dimen.i)/2+1)*
                              NUM_COMP_VECTOR2D+
                              2*Soln_Block_List.Block[i_blk].info.dimen.ghost;
          } /* endif */
          for (j_neigh = 0; j_neigh <= Soln_Block_List.Block[i_blk].nS-1 ; ++j_neigh) {
             l = -1;
             // Unload ghost cell solution information as required.
             if (!Send_Mesh_Geometry_Only) {
                if (!Send_Conservative_Solution_Fluxes) {
                   if (j_neigh == 0) {
                      i_min = Soln_ptr[i_blk].ICl;
	              i_max = Soln_ptr[i_blk].ICl+(Soln_ptr[i_blk].ICu-Soln_ptr[i_blk].ICl-1)/2;
	              i_inc = 1;
	              j_min = Soln_ptr[i_blk].JCl-Soln_Block_List.Block[i_blk].info.dimen.ghost;
	              j_max = Soln_ptr[i_blk].JCl-1;
	              j_inc = 1;
                   } else {
                      i_min = Soln_ptr[i_blk].ICl+(Soln_ptr[i_blk].ICu-Soln_ptr[i_blk].ICl+1)/2;
	              i_max = Soln_ptr[i_blk].ICu;
	              i_inc = 1;
	              j_min = Soln_ptr[i_blk].JCl-Soln_Block_List.Block[i_blk].info.dimen.ghost;
	              j_max = Soln_ptr[i_blk].JCl-1;
	              j_inc = 1;
                   } /* endif */
                   i = Soln_ptr[i_blk].UnloadReceiveBuffer_F2C(Soln_Block_List.message_reschange_southface_recbuf[i_blk][j_neigh],
                                                               l,buffer_size,
                                                               i_min,i_max,i_inc,
                                                               j_min,j_max,j_inc);
	           if (i != 0) return(2100);
                } else {
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
                   i = Soln_ptr[i_blk].UnloadReceiveBuffer_Flux_F2C(Soln_Block_List.message_reschange_southface_recbuf[i_blk][j_neigh],
                                                                    l,buffer_size,
                                                                    i_min,i_max,i_inc,
                                                                    j_min,j_max,j_inc);
	           if (i != 0) return(2101);
                } /* endif */
             } /* endif */
             // Unload ghost cell mesh information as required.
             if (Send_Mesh_Geometry_Only ||
                 Number_of_Solution_Variables > Soln_ptr[i_blk].NumVar()) {
                if (j_neigh == 0) {
                   i_min = Soln_ptr[i_blk].Grid.INl;
	           i_max = Soln_ptr[i_blk].Grid.INl+(Soln_ptr[i_blk].Grid.INu-Soln_ptr[i_blk].Grid.INl)/2;
	           i_inc = 1;
	           j_min = Soln_ptr[i_blk].Grid.JNl-Soln_Block_List.Block[i_blk].info.dimen.ghost;
	           j_max = Soln_ptr[i_blk].Grid.JNl-1;
	           j_inc = 1;
                } else {
                   i_min = Soln_ptr[i_blk].Grid.INl+(Soln_ptr[i_blk].Grid.INu-Soln_ptr[i_blk].Grid.INl)/2;
	           i_max = Soln_ptr[i_blk].Grid.INu;
	           i_inc = 1;
	           j_min = Soln_ptr[i_blk].Grid.JNl-Soln_Block_List.Block[i_blk].info.dimen.ghost;
	           j_max = Soln_ptr[i_blk].Grid.JNl-1;
	           j_inc = 1;
                } /* endif */
                x_ref = Soln_ptr[i_blk].Grid.Node[i_min][j_max+1].X; // Reference node location.
                for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc ) {
                   for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
		      l = l + NUM_COMP_VECTOR2D;
                      if (l >= buffer_size) return(2102);
                      Soln_ptr[i_blk].Grid.Node[i][j].X = 
		        Vector2D(Soln_Block_List.message_reschange_southface_recbuf[i_blk][j_neigh][l-1],
                                 Soln_Block_List.message_reschange_southface_recbuf[i_blk][j_neigh][l])+x_ref;
                   } /* endfor */
                } /* endfor */
                if (j_neigh == 0) {
                   i_min = Soln_ptr[i_blk].ICl;
	           i_max = Soln_ptr[i_blk].ICl+(Soln_ptr[i_blk].ICu-Soln_ptr[i_blk].ICl-1)/2;
	           i_inc = 1;
	           j_min = Soln_ptr[i_blk].JCl-Soln_Block_List.Block[i_blk].info.dimen.ghost;
	           j_max = Soln_ptr[i_blk].JCl-1;
	           j_inc = 1;
                } else {
                   i_min = Soln_ptr[i_blk].ICl+(Soln_ptr[i_blk].ICu-Soln_ptr[i_blk].ICl+1)/2;
	           i_max = Soln_ptr[i_blk].ICu;
	           i_inc = 1;
	           j_min = Soln_ptr[i_blk].JCl-Soln_Block_List.Block[i_blk].info.dimen.ghost;
	           j_max = Soln_ptr[i_blk].JCl-1;
	           j_inc = 1;
                } /* endif */
                for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc ) {
                   for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
	              Soln_ptr[i_blk].Grid.Cell[i][j].Xc = Soln_ptr[i_blk].Grid.centroid(i, j);
	              Soln_ptr[i_blk].Grid.Cell[i][j].A = Soln_ptr[i_blk].Grid.area(i, j);
                   } /* endfor */
                } /* endfor */
                if (j_neigh == 0) {
                   for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc ) {
                      l = l + 1;
                      if (l >= buffer_size) return(2103);
                      Soln_ptr[i_blk].Grid.BCtypeW[j] = 
                         int(Soln_Block_List.message_reschange_southface_recbuf[i_blk][j_neigh][l]);
                   } /* endfor */
                } else {
                   l = l + Soln_Block_List.Block[i_blk].info.dimen.ghost;
                   for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc ) {
                      l = l + 1;
                      if (l >= buffer_size) return(2104);
                      Soln_ptr[i_blk].Grid.BCtypeE[j] = 
                         int(Soln_Block_List.message_reschange_southface_recbuf[i_blk][j_neigh][l]);
                   } /* endfor */
                } /* endif */
             } /* endif */
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
          if (!Send_Mesh_Geometry_Only) {
             if (!Send_Conservative_Solution_Fluxes) {
                buffer_size = Soln_Block_List.Block[i_blk].info.dimen.ghost*
                              (abs(Soln_Block_List.Block[i_blk].info.dimen.j)/2)*
                              Soln_ptr[i_blk].NumVar();
             } else {
                buffer_size = (abs(Soln_Block_List.Block[i_blk].info.dimen.j)/2)*
                              Soln_ptr[i_blk].NumVar();
             } /* endif */
             if (Number_of_Solution_Variables > Soln_ptr[i_blk].NumVar()) 
                buffer_size = buffer_size +
                              Soln_Block_List.Block[i_blk].info.dimen.ghost*
                              (abs(Soln_Block_List.Block[i_blk].info.dimen.j)/2+1)*
                              NUM_COMP_VECTOR2D+
                              2*Soln_Block_List.Block[i_blk].info.dimen.ghost;
          } else {
                buffer_size = Soln_Block_List.Block[i_blk].info.dimen.ghost*
                              (abs(Soln_Block_List.Block[i_blk].info.dimen.j)/2+1)*
                              NUM_COMP_VECTOR2D+
                              2*Soln_Block_List.Block[i_blk].info.dimen.ghost;
          } /* endif */
          for (j_neigh = 0; j_neigh <= Soln_Block_List.Block[i_blk].nE-1 ; ++j_neigh) {
             l = -1;
             // Unload ghost cell solution information as required.
             if (!Send_Mesh_Geometry_Only) {
                if (!Send_Conservative_Solution_Fluxes) {
                   if (j_neigh == 0) {
                      i_min = Soln_ptr[i_blk].ICu+1;
	              i_max = Soln_ptr[i_blk].ICu+Soln_Block_List.Block[i_blk].info.dimen.ghost;
	              i_inc = 1;
	              j_min = Soln_ptr[i_blk].JCl;
	              j_max = Soln_ptr[i_blk].JCl+(Soln_ptr[i_blk].JCu-Soln_ptr[i_blk].JCl-1)/2;
	              j_inc = 1;
                   } else {
                      i_min = Soln_ptr[i_blk].ICu+1;
	              i_max = Soln_ptr[i_blk].ICu+Soln_Block_List.Block[i_blk].info.dimen.ghost;
	              i_inc = 1;
	              j_min = Soln_ptr[i_blk].JCl+(Soln_ptr[i_blk].JCu-Soln_ptr[i_blk].JCl+1)/2;
	              j_max = Soln_ptr[i_blk].JCu;
	              j_inc = 1;
                   } /* endif */
                   i = Soln_ptr[i_blk].UnloadReceiveBuffer_F2C(Soln_Block_List.message_reschange_eastface_recbuf[i_blk][j_neigh],
                                                               l,buffer_size,
                                                               i_min,i_max,i_inc,
                                                               j_min,j_max,j_inc);
	           if (i != 0) return(1200);
                } else {
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
                   i = Soln_ptr[i_blk].UnloadReceiveBuffer_Flux_F2C(Soln_Block_List.message_reschange_eastface_recbuf[i_blk][j_neigh],
                                                                    l,buffer_size,
                                                                    i_min,i_max,i_inc,
                                                                    j_min,j_max,j_inc);
	           if (i != 0) return(1201);
                } /* endif */
             } /* endif */
             // Unload ghost cell mesh information as required.
             if (Send_Mesh_Geometry_Only ||
                 Number_of_Solution_Variables > Soln_ptr[i_blk].NumVar()) {
                if (j_neigh == 0) {
                   i_min = Soln_ptr[i_blk].Grid.INu+1;
	           i_max = Soln_ptr[i_blk].Grid.INu+Soln_Block_List.Block[i_blk].info.dimen.ghost;
	           i_inc = 1;
	           j_min = Soln_ptr[i_blk].Grid.JNl;
	           j_max = Soln_ptr[i_blk].Grid.JNl+(Soln_ptr[i_blk].Grid.JNu-Soln_ptr[i_blk].Grid.JNl)/2;
	           j_inc = 1;
                } else {
		   i_min = Soln_ptr[i_blk].Grid.INu+1;
	           i_max = Soln_ptr[i_blk].Grid.INu+Soln_Block_List.Block[i_blk].info.dimen.ghost;
	           i_inc = 1;
	           j_min = Soln_ptr[i_blk].Grid.JNl+(Soln_ptr[i_blk].Grid.JNu-Soln_ptr[i_blk].Grid.JNl)/2;
	           j_max = Soln_ptr[i_blk].Grid.JNu;
	           j_inc = 1;
                } /* endif */
                x_ref = Soln_ptr[i_blk].Grid.Node[i_min-1][j_min].X; // Reference node location.
                for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc ) {
                   for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
		      l = l + NUM_COMP_VECTOR2D;
                      if (l >= buffer_size) return(1202);
                      Soln_ptr[i_blk].Grid.Node[i][j].X = 
		        Vector2D(Soln_Block_List.message_reschange_eastface_recbuf[i_blk][j_neigh][l-1],
                                 Soln_Block_List.message_reschange_eastface_recbuf[i_blk][j_neigh][l])+x_ref;
                   } /* endfor */
                } /* endfor */
                if (j_neigh == 0) {
                   i_min = Soln_ptr[i_blk].ICu+1;
	           i_max = Soln_ptr[i_blk].ICu+Soln_Block_List.Block[i_blk].info.dimen.ghost;
	           i_inc = 1;
	           j_min = Soln_ptr[i_blk].JCl;
	           j_max = Soln_ptr[i_blk].JCl+(Soln_ptr[i_blk].JCu-Soln_ptr[i_blk].JCl-1)/2;
	           j_inc = 1;
                } else {
                   i_min = Soln_ptr[i_blk].ICu+1;
	           i_max = Soln_ptr[i_blk].ICu+Soln_Block_List.Block[i_blk].info.dimen.ghost;
	           i_inc = 1;
	           j_min = Soln_ptr[i_blk].JCl+(Soln_ptr[i_blk].JCu-Soln_ptr[i_blk].JCl+1)/2;
	           j_max = Soln_ptr[i_blk].JCu;
	           j_inc = 1;
                } /* endif */
                for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc ) {
                   for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
	              Soln_ptr[i_blk].Grid.Cell[i][j].Xc = Soln_ptr[i_blk].Grid.centroid(i, j);
	              Soln_ptr[i_blk].Grid.Cell[i][j].A = Soln_ptr[i_blk].Grid.area(i, j);
                   } /* endfor */
                } /* endfor */
                if (j_neigh == 0) {
                   for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
                      l = l + 1;
                      if (l >= buffer_size) return(1203);
                      Soln_ptr[i_blk].Grid.BCtypeS[i] = 
                         int(Soln_Block_List.message_reschange_eastface_recbuf[i_blk][j_neigh][l]);
                   } /* endfor */
                } else {
                   l = l + Soln_Block_List.Block[i_blk].info.dimen.ghost;
                   for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
                      l = l + 1;
                      if (l >= buffer_size) return(1204);
                      Soln_ptr[i_blk].Grid.BCtypeN[i] = 
                         int(Soln_Block_List.message_reschange_eastface_recbuf[i_blk][j_neigh][l]);
                   } /* endfor */
                } /* endif */
             } /* endif */
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
          if (!Send_Mesh_Geometry_Only) {
             if (!Send_Conservative_Solution_Fluxes) {
                buffer_size = Soln_Block_List.Block[i_blk].info.dimen.ghost*
                              (abs(Soln_Block_List.Block[i_blk].info.dimen.j)/2)*
                              Soln_ptr[i_blk].NumVar();
             } else {
                buffer_size = (abs(Soln_Block_List.Block[i_blk].info.dimen.j)/2)*
                              Soln_ptr[i_blk].NumVar();
             } /* endif */
             if (Number_of_Solution_Variables > Soln_ptr[i_blk].NumVar()) 
                buffer_size = buffer_size +
                              Soln_Block_List.Block[i_blk].info.dimen.ghost*
                              (abs(Soln_Block_List.Block[i_blk].info.dimen.j)/2+1)*
                              NUM_COMP_VECTOR2D+
                              2*Soln_Block_List.Block[i_blk].info.dimen.ghost;
          } else {
                buffer_size = Soln_Block_List.Block[i_blk].info.dimen.ghost*
                              (abs(Soln_Block_List.Block[i_blk].info.dimen.j)/2+1)*
                              NUM_COMP_VECTOR2D+
                              2*Soln_Block_List.Block[i_blk].info.dimen.ghost;
          } /* endif */
          for (j_neigh = 0; j_neigh <= Soln_Block_List.Block[i_blk].nW-1 ; ++j_neigh) {
             l = -1;
             // Unload ghost cell solution information as required.
             if (!Send_Mesh_Geometry_Only) {
                if (!Send_Conservative_Solution_Fluxes) {
                   if (j_neigh == 0) {
                      i_min = Soln_ptr[i_blk].ICl-Soln_Block_List.Block[i_blk].info.dimen.ghost;
	              i_max = Soln_ptr[i_blk].ICl-1;
	              i_inc = 1;
	              j_min = Soln_ptr[i_blk].JCl;
	              j_max = Soln_ptr[i_blk].JCl+(Soln_ptr[i_blk].JCu-Soln_ptr[i_blk].JCl-1)/2;
	              j_inc = 1;
                   } else {
                      i_min = Soln_ptr[i_blk].ICl-Soln_Block_List.Block[i_blk].info.dimen.ghost;
	              i_max = Soln_ptr[i_blk].ICl-1;
	              i_inc = 1;
	              j_min = Soln_ptr[i_blk].JCl+(Soln_ptr[i_blk].JCu-Soln_ptr[i_blk].JCl+1)/2;
	              j_max = Soln_ptr[i_blk].JCu;
	              j_inc = 1;
                   } /* endif */
                   i = Soln_ptr[i_blk].UnloadReceiveBuffer_F2C(Soln_Block_List.message_reschange_westface_recbuf[i_blk][j_neigh],
                                                               l,buffer_size,
                                                               i_min,i_max,i_inc,
                                                               j_min,j_max,j_inc);
	           if (i != 0) return(1100);
                } else {
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
                   i = Soln_ptr[i_blk].UnloadReceiveBuffer_Flux_F2C(Soln_Block_List.message_reschange_westface_recbuf[i_blk][j_neigh],
                                                                    l,buffer_size,
                                                                    i_min,i_max,i_inc,
                                                                    j_min,j_max,j_inc);
	           if (i != 0) return(1101);
                } /* endif */
             } /* endif */
             // Unload ghost cell mesh information as required.
             if (Send_Mesh_Geometry_Only ||
                 Number_of_Solution_Variables > Soln_ptr[i_blk].NumVar()) {
                if (j_neigh == 0) {
                   i_min = Soln_ptr[i_blk].Grid.INl-Soln_Block_List.Block[i_blk].info.dimen.ghost;
	           i_max = Soln_ptr[i_blk].Grid.INl-1;
	           i_inc = 1;
	           j_min = Soln_ptr[i_blk].Grid.JNl;
	           j_max = Soln_ptr[i_blk].Grid.JNl+(Soln_ptr[i_blk].Grid.JNu-Soln_ptr[i_blk].Grid.JNl)/2;
	           j_inc = 1;
                } else {
                   i_min = Soln_ptr[i_blk].Grid.INl-Soln_Block_List.Block[i_blk].info.dimen.ghost;
	           i_max = Soln_ptr[i_blk].Grid.INl-1;
	           i_inc = 1;
	           j_min = Soln_ptr[i_blk].Grid.JNl+(Soln_ptr[i_blk].Grid.JNu-Soln_ptr[i_blk].Grid.JNl)/2;
	           j_max = Soln_ptr[i_blk].Grid.JNu;
	           j_inc = 1;
                } /* endif */
                x_ref = Soln_ptr[i_blk].Grid.Node[i_max+1][j_min].X; // Reference node location.
                for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc ) {
                   for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
		      l = l + NUM_COMP_VECTOR2D;
                      if (l >= buffer_size) return(1102);
                      Soln_ptr[i_blk].Grid.Node[i][j].X = 
		        Vector2D(Soln_Block_List.message_reschange_westface_recbuf[i_blk][j_neigh][l-1],
                                 Soln_Block_List.message_reschange_westface_recbuf[i_blk][j_neigh][l])+x_ref;
                   } /* endfor */
                } /* endfor */
                if (j_neigh == 0) {
                   i_min = Soln_ptr[i_blk].ICl-Soln_Block_List.Block[i_blk].info.dimen.ghost;
	           i_max = Soln_ptr[i_blk].ICl-1;
	           i_inc = 1;
	           j_min = Soln_ptr[i_blk].JCl;
	           j_max = Soln_ptr[i_blk].JCl+(Soln_ptr[i_blk].JCu-Soln_ptr[i_blk].JCl-1)/2;
	           j_inc = 1;
                } else {
                   i_min = Soln_ptr[i_blk].ICl-Soln_Block_List.Block[i_blk].info.dimen.ghost;
	           i_max = Soln_ptr[i_blk].ICl-1;
	           i_inc = 1;
	           j_min = Soln_ptr[i_blk].JCl+(Soln_ptr[i_blk].JCu-Soln_ptr[i_blk].JCl+1)/2;
	           j_max = Soln_ptr[i_blk].JCu;
	           j_inc = 1;
                } /* endif */
                for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc ) {
                   for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
	              Soln_ptr[i_blk].Grid.Cell[i][j].Xc = Soln_ptr[i_blk].Grid.centroid(i, j);
	              Soln_ptr[i_blk].Grid.Cell[i][j].A = Soln_ptr[i_blk].Grid.area(i, j);
                   } /* endfor */
                } /* endfor */
                if (j_neigh == 0) {
                   for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
                      l = l + 1;
                      if (l >= buffer_size) return(1103);
                      Soln_ptr[i_blk].Grid.BCtypeS[i] = 
                         int(Soln_Block_List.message_reschange_westface_recbuf[i_blk][j_neigh][l]);
                   } /* endfor */
                } else {
                   l = l + Soln_Block_List.Block[i_blk].info.dimen.ghost;
                   for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
                      l = l + 1;
                      if (l >= buffer_size) return(1104);
                      Soln_ptr[i_blk].Grid.BCtypeN[i] = 
                         int(Soln_Block_List.message_reschange_westface_recbuf[i_blk][j_neigh][l]);
                   } /* endfor */
                } /* endif */
             } /* endif */
          } /* endfor */
       } /* endif */

       ///////////////////////////////////////////////////////////////////////
       // Unload receive buffer for solution blocks with refined North West //
       // corner neighbours at higher mesh resolution.                      //
       ///////////////////////////////////////////////////////////////////////
       if (Soln_Block_List.Block[i_blk].used &&
           (Soln_Block_List.Block[i_blk].nNW == 1) &&
           (Soln_Block_List.Block[i_blk].info.level <
            Soln_Block_List.Block[i_blk].infoNW[0].level)) {
          buffer_size = sqr(Soln_Block_List.Block[i_blk].info.dimen.ghost)*
                        Number_of_Solution_Variables;
          l = -1;
          // Unload ghost cell solution information as required.
          if (!Send_Mesh_Geometry_Only && !Send_Conservative_Solution_Fluxes) {
             i_min = Soln_ptr[i_blk].ICl-Soln_Block_List.Block[i_blk].info.dimen.ghost;
	     i_max = Soln_ptr[i_blk].ICl-1;
             i_inc = 1;
             j_min = Soln_ptr[i_blk].JCu+1;
             j_max = Soln_ptr[i_blk].JCu+Soln_Block_List.Block[i_blk].info.dimen.ghost;
             j_inc = 1;
             i = Soln_ptr[i_blk].UnloadReceiveBuffer_F2C(Soln_Block_List.message_reschange_northwestcorner_recbuf[i_blk],
                                                         l,buffer_size,
                                                         i_min,i_max,i_inc,
                                                         j_min,j_max,j_inc);
	     if (i != 0) return(4200);
          } /* endif */
          // Unload ghost cell mesh information as required.
          if (Send_Mesh_Geometry_Only ||
              Number_of_Solution_Variables > Soln_ptr[i_blk].NumVar()) {
	     i_min = Soln_ptr[i_blk].Grid.INl-Soln_Block_List.Block[i_blk].info.dimen.ghost;
	     i_max = Soln_ptr[i_blk].Grid.INl-1;
	     i_inc = 1;
	     j_min = Soln_ptr[i_blk].Grid.JNu+1;
	     j_max = Soln_ptr[i_blk].Grid.JNu+Soln_Block_List.Block[i_blk].info.dimen.ghost;
	     j_inc = 1;
             x_ref = Soln_ptr[i_blk].Grid.Node[i_max+1][j_min-1].X; // Reference node location.
             for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc ) {
                for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
	  	      l = l + NUM_COMP_VECTOR2D;
                   if (l >= buffer_size) return(4201);
                   Soln_ptr[i_blk].Grid.Node[i][j].X =
	  	        Vector2D(Soln_Block_List.message_reschange_northwestcorner_recbuf[i_blk][l-1],
                                 Soln_Block_List.message_reschange_northwestcorner_recbuf[i_blk][l])+x_ref;
                } /* endfor */
             } /* endfor */
             i_min = Soln_ptr[i_blk].ICl-Soln_Block_List.Block[i_blk].info.dimen.ghost;
	     i_max = Soln_ptr[i_blk].ICl-1;
	     i_inc = 1;
	     j_min = Soln_ptr[i_blk].JCu+1;
	     j_max = Soln_ptr[i_blk].JCu+Soln_Block_List.Block[i_blk].info.dimen.ghost;
	     j_inc = 1;
             for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc ) {
                for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
	           Soln_ptr[i_blk].Grid.Cell[i][j].Xc = Soln_ptr[i_blk].Grid.centroid(i, j);
	           Soln_ptr[i_blk].Grid.Cell[i][j].A = Soln_ptr[i_blk].Grid.area(i, j);
                } /* endfor */
             } /* endfor */
          } /* endif */
       } /* endif */

       ///////////////////////////////////////////////////////////////////////
       // Unload receive buffer for solution blocks with refined North East //
       // corner neighbours at higher mesh resolution.                      //
       ///////////////////////////////////////////////////////////////////////
       if (Soln_Block_List.Block[i_blk].used &&
           (Soln_Block_List.Block[i_blk].nNE == 1) &&
           (Soln_Block_List.Block[i_blk].info.level <
            Soln_Block_List.Block[i_blk].infoNE[0].level)) {
          buffer_size = sqr(Soln_Block_List.Block[i_blk].info.dimen.ghost)*
                        Number_of_Solution_Variables;
          l = -1;
          // Unload ghost cell solution information as required.
          if (!Send_Mesh_Geometry_Only && !Send_Conservative_Solution_Fluxes) {
             i_min = Soln_ptr[i_blk].ICu+1;
	     i_max = Soln_ptr[i_blk].ICu+Soln_Block_List.Block[i_blk].info.dimen.ghost;
             i_inc = 1;
             j_min = Soln_ptr[i_blk].JCu+1;
             j_max = Soln_ptr[i_blk].JCu+Soln_Block_List.Block[i_blk].info.dimen.ghost;
             j_inc = 1;
             i = Soln_ptr[i_blk].UnloadReceiveBuffer_F2C(Soln_Block_List.message_reschange_northeastcorner_recbuf[i_blk],
                                                         l,buffer_size,
                                                         i_min,i_max,i_inc,
                                                         j_min,j_max,j_inc);
	     if (i != 0) return(5200);
          } /* endif */
          // Unload ghost cell mesh information as required.
          if (Send_Mesh_Geometry_Only ||
              Number_of_Solution_Variables > Soln_ptr[i_blk].NumVar()) {
	     i_min = Soln_ptr[i_blk].Grid.INu+1;
	     i_max = Soln_ptr[i_blk].Grid.INu+Soln_Block_List.Block[i_blk].info.dimen.ghost;
	     i_inc = 1;
	     j_min = Soln_ptr[i_blk].Grid.JNu+1;
	     j_max = Soln_ptr[i_blk].Grid.JNu+Soln_Block_List.Block[i_blk].info.dimen.ghost;
	     j_inc = 1;
             x_ref = Soln_ptr[i_blk].Grid.Node[i_min-1][j_min-1].X; // Reference node location.
             for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc ) {
                for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
	  	      l = l + NUM_COMP_VECTOR2D;
                   if (l >= buffer_size) return(5201);
                   Soln_ptr[i_blk].Grid.Node[i][j].X =
	  	        Vector2D(Soln_Block_List.message_reschange_northeastcorner_recbuf[i_blk][l-1],
                                 Soln_Block_List.message_reschange_northeastcorner_recbuf[i_blk][l])+x_ref;
                } /* endfor */
             } /* endfor */
             i_min = Soln_ptr[i_blk].ICu+1;
	     i_max = Soln_ptr[i_blk].ICu+Soln_Block_List.Block[i_blk].info.dimen.ghost;
	     i_inc = 1;
	     j_min = Soln_ptr[i_blk].JCu+1;
	     j_max = Soln_ptr[i_blk].JCu+Soln_Block_List.Block[i_blk].info.dimen.ghost;
	     j_inc = 1;
             for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc ) {
                for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
	           Soln_ptr[i_blk].Grid.Cell[i][j].Xc = Soln_ptr[i_blk].Grid.centroid(i, j);
	           Soln_ptr[i_blk].Grid.Cell[i][j].A = Soln_ptr[i_blk].Grid.area(i, j);
                } /* endfor */
             } /* endfor */
          } /* endif */
       } /* endif */

       ///////////////////////////////////////////////////////////////////////
       // Unload receive buffer for solution blocks with refined South East //
       // corner neighbours at higher mesh resolution.                      //
       ///////////////////////////////////////////////////////////////////////
       if (Soln_Block_List.Block[i_blk].used &&
           (Soln_Block_List.Block[i_blk].nSE == 1) &&
           (Soln_Block_List.Block[i_blk].info.level <
            Soln_Block_List.Block[i_blk].infoSE[0].level)) {
          buffer_size = sqr(Soln_Block_List.Block[i_blk].info.dimen.ghost)*
                        Number_of_Solution_Variables;
          l = -1;
          // Unload ghost cell solution information as required.
          if (!Send_Mesh_Geometry_Only && !Send_Conservative_Solution_Fluxes) {
             i_min = Soln_ptr[i_blk].ICu+1;
	     i_max = Soln_ptr[i_blk].ICu+Soln_Block_List.Block[i_blk].info.dimen.ghost;
             i_inc = 1;
             j_min = Soln_ptr[i_blk].JCl-Soln_Block_List.Block[i_blk].info.dimen.ghost;
             j_max = Soln_ptr[i_blk].JCl-1;
             j_inc = 1;
             i = Soln_ptr[i_blk].UnloadReceiveBuffer_F2C(Soln_Block_List.message_reschange_southeastcorner_recbuf[i_blk],
                                                         l,buffer_size,
                                                         i_min,i_max,i_inc,
                                                         j_min,j_max,j_inc);
	     if (i != 0) return(4100);
          } /* endif */
          // Unload ghost cell mesh information as required.
          if (Send_Mesh_Geometry_Only ||
              Number_of_Solution_Variables > Soln_ptr[i_blk].NumVar()) {
	     i_min = Soln_ptr[i_blk].Grid.INu+1;
	     i_max = Soln_ptr[i_blk].Grid.INu+Soln_Block_List.Block[i_blk].info.dimen.ghost;
	     i_inc = 1;
	     j_min = Soln_ptr[i_blk].Grid.JNl-Soln_Block_List.Block[i_blk].info.dimen.ghost;
	     j_max = Soln_ptr[i_blk].Grid.JNl-1;
	     j_inc = 1;
             x_ref = Soln_ptr[i_blk].Grid.Node[i_min-1][j_max+1].X; // Reference node location.
             for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc ) {
                for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
	  	      l = l + NUM_COMP_VECTOR2D;
                   if (l >= buffer_size) return(4101);
                   Soln_ptr[i_blk].Grid.Node[i][j].X =
	  	        Vector2D(Soln_Block_List.message_reschange_southeastcorner_recbuf[i_blk][l-1],
                                 Soln_Block_List.message_reschange_southeastcorner_recbuf[i_blk][l])+x_ref;
                } /* endfor */
             } /* endfor */
             i_min = Soln_ptr[i_blk].ICu+1;
	     i_max = Soln_ptr[i_blk].ICu+Soln_Block_List.Block[i_blk].info.dimen.ghost;
	     i_inc = 1;
	     j_min = Soln_ptr[i_blk].JCl-Soln_Block_List.Block[i_blk].info.dimen.ghost;
	     j_max = Soln_ptr[i_blk].JCl-1;
	     j_inc = 1;
             for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc ) {
                for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
	           Soln_ptr[i_blk].Grid.Cell[i][j].Xc = Soln_ptr[i_blk].Grid.centroid(i, j);
	           Soln_ptr[i_blk].Grid.Cell[i][j].A = Soln_ptr[i_blk].Grid.area(i, j);
                } /* endfor */
             } /* endfor */
          } /* endif */
       } /* endif */

       ///////////////////////////////////////////////////////////////////////
       // Unload receive buffer for solution blocks with refined South West //
       // corner neighbours at higher mesh resolution.                      //
       ///////////////////////////////////////////////////////////////////////
       if (Soln_Block_List.Block[i_blk].used &&
           (Soln_Block_List.Block[i_blk].nSW == 1) &&
           (Soln_Block_List.Block[i_blk].info.level <
            Soln_Block_List.Block[i_blk].infoSW[0].level)) {
          buffer_size = sqr(Soln_Block_List.Block[i_blk].info.dimen.ghost)*
                        Number_of_Solution_Variables;
          l = -1;
          // Unload ghost cell solution information as required.
          if (!Send_Mesh_Geometry_Only && !Send_Conservative_Solution_Fluxes) {
             i_min = Soln_ptr[i_blk].ICl-Soln_Block_List.Block[i_blk].info.dimen.ghost;
	     i_max = Soln_ptr[i_blk].ICl-1;
             i_inc = 1;
             j_min = Soln_ptr[i_blk].JCl-Soln_Block_List.Block[i_blk].info.dimen.ghost;
             j_max = Soln_ptr[i_blk].JCl-1;
             j_inc = 1;
             i = Soln_ptr[i_blk].UnloadReceiveBuffer_F2C(Soln_Block_List.message_reschange_southwestcorner_recbuf[i_blk],
                                                         l,buffer_size,
                                                         i_min,i_max,i_inc,
                                                         j_min,j_max,j_inc);
	     if (i != 0) return(5100);
          } /* endif */
          // Unload ghost cell mesh information as required.
          if (Send_Mesh_Geometry_Only ||
              Number_of_Solution_Variables > Soln_ptr[i_blk].NumVar()) {
	     i_min = Soln_ptr[i_blk].Grid.INl-Soln_Block_List.Block[i_blk].info.dimen.ghost;
	     i_max = Soln_ptr[i_blk].Grid.INl-1;
	     i_inc = 1;
	     j_min = Soln_ptr[i_blk].Grid.JNl-Soln_Block_List.Block[i_blk].info.dimen.ghost;
	     j_max = Soln_ptr[i_blk].Grid.JNl-1;
	     j_inc = 1;
             x_ref = Soln_ptr[i_blk].Grid.Node[i_max+1][j_max+1].X; // Reference node location.
             for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc ) {
                for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
	  	      l = l + NUM_COMP_VECTOR2D;
                   if (l >= buffer_size) return(5101);
                   Soln_ptr[i_blk].Grid.Node[i][j].X =
	  	        Vector2D(Soln_Block_List.message_reschange_southwestcorner_recbuf[i_blk][l-1],
                                 Soln_Block_List.message_reschange_southwestcorner_recbuf[i_blk][l])+x_ref;
                } /* endfor */
             } /* endfor */
             i_min = Soln_ptr[i_blk].ICl-Soln_Block_List.Block[i_blk].info.dimen.ghost;
	     i_max = Soln_ptr[i_blk].ICl-1;
	     i_inc = 1;
	     j_min = Soln_ptr[i_blk].JCl-Soln_Block_List.Block[i_blk].info.dimen.ghost;
	     j_max = Soln_ptr[i_blk].JCl-1;
	     j_inc = 1;
             for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc ) {
                for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
	           Soln_ptr[i_blk].Grid.Cell[i][j].Xc = Soln_ptr[i_blk].Grid.centroid(i, j);
	           Soln_ptr[i_blk].Grid.Cell[i][j].A = Soln_ptr[i_blk].Grid.area(i, j);
                } /* endfor */
             } /* endfor */
          } /* endif */
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
template <class Quad_Soln_Block>
int Unload_Receive_Message_Buffers_ResChange_CoarseToFine(Quad_Soln_Block *Soln_ptr,
                                                          AdaptiveBlock2D_List &Soln_Block_List,
                                                          const int Number_of_Solution_Variables,
                                                          const int Send_Mesh_Geometry_Only) {

    int i_blk, buffer_size, 
        i_min, i_max, i_inc, i, 
        j_min, j_max, j_inc, j, 
        l, j_neigh;
    Vector2D x_ref;

    /* Check to see if the number of solution variables specified
       in the calling argument is correct. */

    if (Number_of_Solution_Variables != NUM_COMP_VECTOR2D &&
        Number_of_Solution_Variables != Soln_ptr[0].NumVar() &&
        Number_of_Solution_Variables != Soln_ptr[0].NumVar() + NUM_COMP_VECTOR2D) {
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
          if (!Send_Mesh_Geometry_Only) {
             buffer_size = Soln_Block_List.Block[i_blk].info.dimen.ghost*
                           (abs(Soln_Block_List.Block[i_blk].info.dimen.i)+
                           Soln_Block_List.Block[i_blk].info.dimen.ghost)*
                           Soln_ptr[i_blk].NumVar();
             if (Number_of_Solution_Variables > Soln_ptr[i_blk].NumVar()) 
                buffer_size = buffer_size +
                              Soln_Block_List.Block[i_blk].info.dimen.ghost*
                              (abs(Soln_Block_List.Block[i_blk].info.dimen.i)+
                              Soln_Block_List.Block[i_blk].info.dimen.ghost+1)*
                              NUM_COMP_VECTOR2D+
                              2*Soln_Block_List.Block[i_blk].info.dimen.ghost;
          } else {
             buffer_size = Soln_Block_List.Block[i_blk].info.dimen.ghost*
                           (abs(Soln_Block_List.Block[i_blk].info.dimen.i)+
                           Soln_Block_List.Block[i_blk].info.dimen.ghost+1)*
                           NUM_COMP_VECTOR2D+
                           2*Soln_Block_List.Block[i_blk].info.dimen.ghost;
          } /* endif */
          l = -1;
          // Unload ghost cell solution information as required.
          if (!Send_Mesh_Geometry_Only) {
 	     if (Soln_Block_List.Block[i_blk].info.sector == ADAPTIVEBLOCK2D_SECTOR_NW) {
                i_min = Soln_ptr[i_blk].ICl;
	        i_max = Soln_ptr[i_blk].ICu+Soln_Block_List.Block[i_blk].info.dimen.ghost;
                i_inc = 1;
                j_min = Soln_ptr[i_blk].JCu+1;
                j_max = Soln_ptr[i_blk].JCu+Soln_Block_List.Block[i_blk].info.dimen.ghost;
                j_inc = 1;
             } else if (Soln_Block_List.Block[i_blk].info.sector == ADAPTIVEBLOCK2D_SECTOR_NE) {
                i_min = Soln_ptr[i_blk].ICl-Soln_Block_List.Block[i_blk].info.dimen.ghost;
	        i_max = Soln_ptr[i_blk].ICu;
                i_inc = 1;
                j_min = Soln_ptr[i_blk].JCu+1;
                j_max = Soln_ptr[i_blk].JCu+Soln_Block_List.Block[i_blk].info.dimen.ghost;
                j_inc = 1;
             } /* endif */
             i = Soln_ptr[i_blk].UnloadReceiveBuffer_C2F(Soln_Block_List.message_reschange_northface_recbuf[i_blk][0],
                                                         l,buffer_size,
                                                         i_min,i_max,i_inc,
                                                         j_min,j_max,j_inc);
	     if (i != 0) return(2200);
          } /* endif */
          // Unload ghost cell mesh information as required.
          if (Send_Mesh_Geometry_Only ||
              Number_of_Solution_Variables > Soln_ptr[i_blk].NumVar()) {
 	     if (Soln_Block_List.Block[i_blk].info.sector == ADAPTIVEBLOCK2D_SECTOR_NW) {
                i_min = Soln_ptr[i_blk].Grid.INl;
                i_max = Soln_ptr[i_blk].Grid.INu+Soln_Block_List.Block[i_blk].info.dimen.ghost;
                i_inc = 1;
                j_min = Soln_ptr[i_blk].Grid.JNu+1;
                j_max = Soln_ptr[i_blk].Grid.JNu+Soln_Block_List.Block[i_blk].info.dimen.ghost;
                j_inc = 1;
             } else if (Soln_Block_List.Block[i_blk].info.sector == ADAPTIVEBLOCK2D_SECTOR_NE) {
                i_min = Soln_ptr[i_blk].Grid.INl-Soln_Block_List.Block[i_blk].info.dimen.ghost;
                i_max = Soln_ptr[i_blk].Grid.INu;
                i_inc = 1;
                j_min = Soln_ptr[i_blk].Grid.JNu+1;
                j_max = Soln_ptr[i_blk].Grid.JNu+Soln_Block_List.Block[i_blk].info.dimen.ghost;
                j_inc = 1;
             } /* endif */             
             x_ref = Soln_ptr[i_blk].Grid.Node[Soln_ptr[i_blk].Grid.INl][Soln_ptr[i_blk].Grid.JNu].X; // Reference node location.
             for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc ) {
                for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
		   l = l + NUM_COMP_VECTOR2D;
                   if (l >= buffer_size) return(2201);
                   Soln_ptr[i_blk].Grid.Node[i][j].X = 
		     Vector2D(Soln_Block_List.message_reschange_northface_recbuf[i_blk][0][l-1],
                              Soln_Block_List.message_reschange_northface_recbuf[i_blk][0][l])+x_ref;
                } /* endfor */
             } /* endfor */
 	     if (Soln_Block_List.Block[i_blk].info.sector == ADAPTIVEBLOCK2D_SECTOR_NW) {
                i_min = Soln_ptr[i_blk].ICl;
	        i_max = Soln_ptr[i_blk].ICu+Soln_Block_List.Block[i_blk].info.dimen.ghost;
                i_inc = 1;
                j_min = Soln_ptr[i_blk].JCu+1;
                j_max = Soln_ptr[i_blk].JCu+Soln_Block_List.Block[i_blk].info.dimen.ghost;
                j_inc = 1;
             } else if (Soln_Block_List.Block[i_blk].info.sector == ADAPTIVEBLOCK2D_SECTOR_NE) {
                i_min = Soln_ptr[i_blk].ICl-Soln_Block_List.Block[i_blk].info.dimen.ghost;
	        i_max = Soln_ptr[i_blk].ICu;
                i_inc = 1;
                j_min = Soln_ptr[i_blk].JCu+1;
                j_max = Soln_ptr[i_blk].JCu+Soln_Block_List.Block[i_blk].info.dimen.ghost;
                j_inc = 1;
             } /* endif */
             for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc ) {
                for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
	           Soln_ptr[i_blk].Grid.Cell[i][j].Xc = Soln_ptr[i_blk].Grid.centroid(i, j);
	           Soln_ptr[i_blk].Grid.Cell[i][j].A = Soln_ptr[i_blk].Grid.area(i, j);
                } /* endfor */
             } /* endfor */
 	     if (Soln_Block_List.Block[i_blk].info.sector == ADAPTIVEBLOCK2D_SECTOR_NW) {
                for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc ) {
                   l = l + 1;
                   if (l >= buffer_size) return(2202);
                   Soln_ptr[i_blk].Grid.BCtypeW[j] = 
                      int(Soln_Block_List.message_reschange_northface_recbuf[i_blk][0][l]);
                } /* endfor */
             } else if (Soln_Block_List.Block[i_blk].info.sector == ADAPTIVEBLOCK2D_SECTOR_NE) {
                for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc ) {
                   l = l + 1;
                   if (l >= buffer_size) return(2203);
                   Soln_ptr[i_blk].Grid.BCtypeE[j] = 
                      int(Soln_Block_List.message_reschange_northface_recbuf[i_blk][0][l]);
                } /* endfor */
             } /* endif */
          } /* endif */
       } /* endif */

       ///////////////////////////////////////////////////////////////////////
       // Unload receive buffer for solution blocks with coarse South face  //
       // neighbours at lower mesh resolution.                              //
       ///////////////////////////////////////////////////////////////////////
       if (Soln_Block_List.Block[i_blk].used && 
           (Soln_Block_List.Block[i_blk].nS == 1) &&
           (Soln_Block_List.Block[i_blk].info.level > 
            Soln_Block_List.Block[i_blk].infoS[0].level)) {
          if (!Send_Mesh_Geometry_Only) {
             buffer_size = Soln_Block_List.Block[i_blk].info.dimen.ghost*
                           (abs(Soln_Block_List.Block[i_blk].info.dimen.i)+
                           Soln_Block_List.Block[i_blk].info.dimen.ghost)*
                           Soln_ptr[i_blk].NumVar();
             if (Number_of_Solution_Variables > Soln_ptr[i_blk].NumVar()) 
                buffer_size = buffer_size +
                              Soln_Block_List.Block[i_blk].info.dimen.ghost*
                              (abs(Soln_Block_List.Block[i_blk].info.dimen.i)+
                              Soln_Block_List.Block[i_blk].info.dimen.ghost+1)*
                              NUM_COMP_VECTOR2D+
                              2*Soln_Block_List.Block[i_blk].info.dimen.ghost;
          } else {
             buffer_size = Soln_Block_List.Block[i_blk].info.dimen.ghost*
                           (abs(Soln_Block_List.Block[i_blk].info.dimen.i)+
                           Soln_Block_List.Block[i_blk].info.dimen.ghost+1)*
                           NUM_COMP_VECTOR2D+
                           2*Soln_Block_List.Block[i_blk].info.dimen.ghost;
          } /* endif */
          l = -1;
          // Unload ghost cell solution information as required.
          if (!Send_Mesh_Geometry_Only) {
 	     if (Soln_Block_List.Block[i_blk].info.sector == ADAPTIVEBLOCK2D_SECTOR_SW) {
                i_min = Soln_ptr[i_blk].ICl;
	        i_max = Soln_ptr[i_blk].ICu+Soln_Block_List.Block[i_blk].info.dimen.ghost;
                i_inc = 1;
                j_min = Soln_ptr[i_blk].JCl-Soln_Block_List.Block[i_blk].info.dimen.ghost;
                j_max = Soln_ptr[i_blk].JCl-1;
                j_inc = 1;
             } else if (Soln_Block_List.Block[i_blk].info.sector == ADAPTIVEBLOCK2D_SECTOR_SE) {
                i_min = Soln_ptr[i_blk].ICl-Soln_Block_List.Block[i_blk].info.dimen.ghost;
	        i_max = Soln_ptr[i_blk].ICu;
                i_inc = 1;
                j_min = Soln_ptr[i_blk].JCl-Soln_Block_List.Block[i_blk].info.dimen.ghost;
                j_max = Soln_ptr[i_blk].JCl-1;
                j_inc = 1;
             } /* endif */
             i = Soln_ptr[i_blk].UnloadReceiveBuffer_C2F(Soln_Block_List.message_reschange_southface_recbuf[i_blk][0],
                                                         l,buffer_size,
                                                         i_min,i_max,i_inc,
                                                         j_min,j_max,j_inc);
	     if (i != 0) return(2100);
          } /* endif */
          // Unload ghost cell mesh information as required.
          if (Send_Mesh_Geometry_Only ||
              Number_of_Solution_Variables > Soln_ptr[i_blk].NumVar()) {
 	     if (Soln_Block_List.Block[i_blk].info.sector == ADAPTIVEBLOCK2D_SECTOR_SW) {
                i_min = Soln_ptr[i_blk].Grid.INl;
                i_max = Soln_ptr[i_blk].Grid.INu+Soln_Block_List.Block[i_blk].info.dimen.ghost;
                i_inc = 1;
                j_min = Soln_ptr[i_blk].Grid.JNl-Soln_Block_List.Block[i_blk].info.dimen.ghost;
                j_max = Soln_ptr[i_blk].Grid.JNl-1;
                j_inc = 1;
             } else if (Soln_Block_List.Block[i_blk].info.sector == ADAPTIVEBLOCK2D_SECTOR_SE) {
                i_min = Soln_ptr[i_blk].Grid.INl-Soln_Block_List.Block[i_blk].info.dimen.ghost;
                i_max = Soln_ptr[i_blk].Grid.INu;
                i_inc = 1;
                j_min = Soln_ptr[i_blk].Grid.JNl-Soln_Block_List.Block[i_blk].info.dimen.ghost;
                j_max = Soln_ptr[i_blk].Grid.JNl-1;
                j_inc = 1;
             } /* endif */             
             x_ref = Soln_ptr[i_blk].Grid.Node[Soln_ptr[i_blk].Grid.INl][Soln_ptr[i_blk].Grid.JNl].X; // Reference node location.
             for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc ) {
                for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
		   l = l + NUM_COMP_VECTOR2D;
                   if (l >= buffer_size) return(2101);
                   Soln_ptr[i_blk].Grid.Node[i][j].X = 
		     Vector2D(Soln_Block_List.message_reschange_southface_recbuf[i_blk][0][l-1],
                              Soln_Block_List.message_reschange_southface_recbuf[i_blk][0][l])+x_ref;
                } /* endfor */
             } /* endfor */
 	     if (Soln_Block_List.Block[i_blk].info.sector == ADAPTIVEBLOCK2D_SECTOR_SW) {
                i_min = Soln_ptr[i_blk].ICl;
	        i_max = Soln_ptr[i_blk].ICu+Soln_Block_List.Block[i_blk].info.dimen.ghost;
                i_inc = 1;
                j_min = Soln_ptr[i_blk].JCl-Soln_Block_List.Block[i_blk].info.dimen.ghost;
                j_max = Soln_ptr[i_blk].JCl-1;
                j_inc = 1;
             } else if (Soln_Block_List.Block[i_blk].info.sector == ADAPTIVEBLOCK2D_SECTOR_SE) {
                i_min = Soln_ptr[i_blk].ICl-Soln_Block_List.Block[i_blk].info.dimen.ghost;
	        i_max = Soln_ptr[i_blk].ICu;
                i_inc = 1;
                j_min = Soln_ptr[i_blk].JCl-Soln_Block_List.Block[i_blk].info.dimen.ghost;
                j_max = Soln_ptr[i_blk].JCl-1;
                j_inc = 1;
             } /* endif */
             for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc ) {
                for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
	           Soln_ptr[i_blk].Grid.Cell[i][j].Xc = Soln_ptr[i_blk].Grid.centroid(i, j);
	           Soln_ptr[i_blk].Grid.Cell[i][j].A = Soln_ptr[i_blk].Grid.area(i, j);
                } /* endfor */
             } /* endfor */
 	     if (Soln_Block_List.Block[i_blk].info.sector == ADAPTIVEBLOCK2D_SECTOR_SW) {
                for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc ) {
                   l = l + 1;
                   if (l >= buffer_size) return(2102);
                   Soln_ptr[i_blk].Grid.BCtypeW[j] = 
                      int(Soln_Block_List.message_reschange_southface_recbuf[i_blk][0][l]);
                } /* endfor */
             } else if (Soln_Block_List.Block[i_blk].info.sector == ADAPTIVEBLOCK2D_SECTOR_SE) {
                for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc ) {
                   l = l + 1;
                   if (l >= buffer_size) return(2103);
                   Soln_ptr[i_blk].Grid.BCtypeE[j] = 
                      int(Soln_Block_List.message_reschange_southface_recbuf[i_blk][0][l]);
                } /* endfor */
             } /* endif */
          } /* endif */
       } /* endif */

       ///////////////////////////////////////////////////////////////////////
       // Unload receive buffer for solution blocks with coarse East face   //
       // neighbours at lower mesh resolution.                              //
       ///////////////////////////////////////////////////////////////////////
       if (Soln_Block_List.Block[i_blk].used && 
           (Soln_Block_List.Block[i_blk].nE == 1) &&
           (Soln_Block_List.Block[i_blk].info.level > 
            Soln_Block_List.Block[i_blk].infoE[0].level)) {
          if (!Send_Mesh_Geometry_Only) {
             buffer_size = Soln_Block_List.Block[i_blk].info.dimen.ghost*
                           (abs(Soln_Block_List.Block[i_blk].info.dimen.j)+
                           Soln_Block_List.Block[i_blk].info.dimen.ghost)*
                           Soln_ptr[i_blk].NumVar();
             if (Number_of_Solution_Variables > Soln_ptr[i_blk].NumVar()) 
                buffer_size = buffer_size +
                              Soln_Block_List.Block[i_blk].info.dimen.ghost*
                              (abs(Soln_Block_List.Block[i_blk].info.dimen.j)+
                              Soln_Block_List.Block[i_blk].info.dimen.ghost+1)*
                              NUM_COMP_VECTOR2D+
                              2*Soln_Block_List.Block[i_blk].info.dimen.ghost;
          } else {
             buffer_size = Soln_Block_List.Block[i_blk].info.dimen.ghost*
                           (abs(Soln_Block_List.Block[i_blk].info.dimen.j)+
                           Soln_Block_List.Block[i_blk].info.dimen.ghost+1)*
                           NUM_COMP_VECTOR2D+
                           2*Soln_Block_List.Block[i_blk].info.dimen.ghost;
          } /* endif */
          l = -1;
          // Unload ghost cell solution information as required.
          if (!Send_Mesh_Geometry_Only) {
 	     if (Soln_Block_List.Block[i_blk].info.sector == ADAPTIVEBLOCK2D_SECTOR_SE) {
                i_min = Soln_ptr[i_blk].ICu+1;
                i_max = Soln_ptr[i_blk].ICu+Soln_Block_List.Block[i_blk].info.dimen.ghost;
                i_inc = 1;
                j_min = Soln_ptr[i_blk].JCl;
	        j_max = Soln_ptr[i_blk].JCu+Soln_Block_List.Block[i_blk].info.dimen.ghost;
                j_inc = 1;
             } else if (Soln_Block_List.Block[i_blk].info.sector == ADAPTIVEBLOCK2D_SECTOR_NE) {
                i_min = Soln_ptr[i_blk].ICu+1;
                i_max = Soln_ptr[i_blk].ICu+Soln_Block_List.Block[i_blk].info.dimen.ghost;
                i_inc = 1;
                j_min = Soln_ptr[i_blk].JCl-Soln_Block_List.Block[i_blk].info.dimen.ghost;
	        j_max = Soln_ptr[i_blk].JCu;
                j_inc = 1;
             } /* endif */
             i = Soln_ptr[i_blk].UnloadReceiveBuffer_C2F(Soln_Block_List.message_reschange_eastface_recbuf[i_blk][0],
                                                         l,buffer_size,
                                                         i_min,i_max,i_inc,
                                                         j_min,j_max,j_inc);
	     if (i != 0) return(1200);
          } /* endif */
          // Unload ghost cell mesh information as required.
          if (Send_Mesh_Geometry_Only ||
              Number_of_Solution_Variables > Soln_ptr[i_blk].NumVar()) {
 	     if (Soln_Block_List.Block[i_blk].info.sector == ADAPTIVEBLOCK2D_SECTOR_SE) {
                i_min = Soln_ptr[i_blk].Grid.INu+1;
                i_max = Soln_ptr[i_blk].Grid.INu+Soln_Block_List.Block[i_blk].info.dimen.ghost;
                i_inc = 1;
                j_min = Soln_ptr[i_blk].Grid.JNl;
                j_max = Soln_ptr[i_blk].Grid.JNu+Soln_Block_List.Block[i_blk].info.dimen.ghost;
                j_inc = 1;
             } else if (Soln_Block_List.Block[i_blk].info.sector == ADAPTIVEBLOCK2D_SECTOR_NE) {
                i_min = Soln_ptr[i_blk].Grid.INu+1;
                i_max = Soln_ptr[i_blk].Grid.INu+Soln_Block_List.Block[i_blk].info.dimen.ghost;
                i_inc = 1;
                j_min = Soln_ptr[i_blk].Grid.JNl-Soln_Block_List.Block[i_blk].info.dimen.ghost;
                j_max = Soln_ptr[i_blk].Grid.JNu;
                j_inc = 1;
             } /* endif */             
             x_ref = Soln_ptr[i_blk].Grid.Node[Soln_ptr[i_blk].Grid.INu][Soln_ptr[i_blk].Grid.JNl].X; // Reference node location.
             for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc ) {
                for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
		   l = l + NUM_COMP_VECTOR2D;
                   if (l >= buffer_size) return(1201);
                   Soln_ptr[i_blk].Grid.Node[i][j].X = 
		     Vector2D(Soln_Block_List.message_reschange_eastface_recbuf[i_blk][0][l-1],
                              Soln_Block_List.message_reschange_eastface_recbuf[i_blk][0][l])+x_ref;
                } /* endfor */
             } /* endfor */
 	     if (Soln_Block_List.Block[i_blk].info.sector == ADAPTIVEBLOCK2D_SECTOR_SE) {
                i_min = Soln_ptr[i_blk].ICu+1;
                i_max = Soln_ptr[i_blk].ICu+Soln_Block_List.Block[i_blk].info.dimen.ghost;
                i_inc = 1;
                j_min = Soln_ptr[i_blk].JCl;
	        j_max = Soln_ptr[i_blk].JCu+Soln_Block_List.Block[i_blk].info.dimen.ghost;
                j_inc = 1;
             } else if (Soln_Block_List.Block[i_blk].info.sector == ADAPTIVEBLOCK2D_SECTOR_NE) {
                i_min = Soln_ptr[i_blk].ICu+1;
                i_max = Soln_ptr[i_blk].ICu+Soln_Block_List.Block[i_blk].info.dimen.ghost;
                i_inc = 1;
                j_min = Soln_ptr[i_blk].JCl-Soln_Block_List.Block[i_blk].info.dimen.ghost;
	        j_max = Soln_ptr[i_blk].JCu;
                j_inc = 1;
             } /* endif */
             for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc ) {
                for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
	           Soln_ptr[i_blk].Grid.Cell[i][j].Xc = Soln_ptr[i_blk].Grid.centroid(i, j);
	           Soln_ptr[i_blk].Grid.Cell[i][j].A = Soln_ptr[i_blk].Grid.area(i, j);
                } /* endfor */
             } /* endfor */
 	     if (Soln_Block_List.Block[i_blk].info.sector == ADAPTIVEBLOCK2D_SECTOR_SE) {
                for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
                   l = l + 1;
                   if (l >= buffer_size) return(1202);
                   Soln_ptr[i_blk].Grid.BCtypeS[i] = 
                      int(Soln_Block_List.message_reschange_eastface_recbuf[i_blk][0][l]);
                } /* endfor */
             } else if (Soln_Block_List.Block[i_blk].info.sector == ADAPTIVEBLOCK2D_SECTOR_NE) {
                for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
                   l = l + 1;
                   if (l >= buffer_size) return(1203);
                   Soln_ptr[i_blk].Grid.BCtypeN[i] = 
                      int(Soln_Block_List.message_reschange_eastface_recbuf[i_blk][0][l]);
                } /* endfor */
             } /* endif */
          } /* endif */
       } /* endif */

       ///////////////////////////////////////////////////////////////////////
       // Unload receive buffer for solution blocks with coarse West face   //
       // neighbours at lower mesh resolution.                              //
       ///////////////////////////////////////////////////////////////////////
       if (Soln_Block_List.Block[i_blk].used && 
           (Soln_Block_List.Block[i_blk].nW == 1) &&
           (Soln_Block_List.Block[i_blk].info.level > 
            Soln_Block_List.Block[i_blk].infoW[0].level)) {
          if (!Send_Mesh_Geometry_Only) {
             buffer_size = Soln_Block_List.Block[i_blk].info.dimen.ghost*
                           (abs(Soln_Block_List.Block[i_blk].info.dimen.j)+
                           Soln_Block_List.Block[i_blk].info.dimen.ghost)*
                           Soln_ptr[i_blk].NumVar();
             if (Number_of_Solution_Variables > Soln_ptr[i_blk].NumVar()) 
                buffer_size = buffer_size +
                              Soln_Block_List.Block[i_blk].info.dimen.ghost*
                              (abs(Soln_Block_List.Block[i_blk].info.dimen.j)+
                              Soln_Block_List.Block[i_blk].info.dimen.ghost+1)*
                              NUM_COMP_VECTOR2D+
                              2*Soln_Block_List.Block[i_blk].info.dimen.ghost;
          } else {
             buffer_size = Soln_Block_List.Block[i_blk].info.dimen.ghost*
                           (abs(Soln_Block_List.Block[i_blk].info.dimen.j)+
                           Soln_Block_List.Block[i_blk].info.dimen.ghost+1)*
                           NUM_COMP_VECTOR2D+
                           2*Soln_Block_List.Block[i_blk].info.dimen.ghost;
          } /* endif */
          l = -1;
          // Unload ghost cell solution information as required.
          if (!Send_Mesh_Geometry_Only) {
 	     if (Soln_Block_List.Block[i_blk].info.sector == ADAPTIVEBLOCK2D_SECTOR_SW) {
                i_min = Soln_ptr[i_blk].ICl-Soln_Block_List.Block[i_blk].info.dimen.ghost;
                i_max = Soln_ptr[i_blk].ICl-1;
                i_inc = 1;
                j_min = Soln_ptr[i_blk].JCl;
	        j_max = Soln_ptr[i_blk].JCu+Soln_Block_List.Block[i_blk].info.dimen.ghost;
                j_inc = 1;
             } else if (Soln_Block_List.Block[i_blk].info.sector == ADAPTIVEBLOCK2D_SECTOR_NW) {
                i_min = Soln_ptr[i_blk].ICl-Soln_Block_List.Block[i_blk].info.dimen.ghost;
                i_max = Soln_ptr[i_blk].ICl-1;
                i_inc = 1;
                j_min = Soln_ptr[i_blk].JCl-Soln_Block_List.Block[i_blk].info.dimen.ghost;
	        j_max = Soln_ptr[i_blk].JCu;
                j_inc = 1;
             } /* endif */
             i = Soln_ptr[i_blk].UnloadReceiveBuffer_C2F(Soln_Block_List.message_reschange_westface_recbuf[i_blk][0],
                                                         l,buffer_size,
                                                         i_min,i_max,i_inc,
                                                         j_min,j_max,j_inc);
	     if (i != 0) return(1100);
          } /* endif */
          // Unload ghost cell mesh information as required.
          if (Send_Mesh_Geometry_Only ||
              Number_of_Solution_Variables > Soln_ptr[i_blk].NumVar()) {
 	     if (Soln_Block_List.Block[i_blk].info.sector == ADAPTIVEBLOCK2D_SECTOR_SW) {
                i_min = Soln_ptr[i_blk].Grid.INl-Soln_Block_List.Block[i_blk].info.dimen.ghost;
                i_max = Soln_ptr[i_blk].Grid.INl-1;
                i_inc = 1;
                j_min = Soln_ptr[i_blk].Grid.JNl;
                j_max = Soln_ptr[i_blk].Grid.JNu+Soln_Block_List.Block[i_blk].info.dimen.ghost;
                j_inc = 1;
             } else if (Soln_Block_List.Block[i_blk].info.sector == ADAPTIVEBLOCK2D_SECTOR_NW) {
                i_min = Soln_ptr[i_blk].Grid.INl-Soln_Block_List.Block[i_blk].info.dimen.ghost;
                i_max = Soln_ptr[i_blk].Grid.INl-1;
                i_inc = 1;
                j_min = Soln_ptr[i_blk].Grid.JNl-Soln_Block_List.Block[i_blk].info.dimen.ghost;
                j_max = Soln_ptr[i_blk].Grid.JNu;
                j_inc = 1;
             } /* endif */             
             x_ref = Soln_ptr[i_blk].Grid.Node[Soln_ptr[i_blk].Grid.INl][Soln_ptr[i_blk].Grid.JNl].X; // Reference node location.
             for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc ) {
                for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
		   l = l + NUM_COMP_VECTOR2D;
                   if (l >= buffer_size) return(1101);
                   Soln_ptr[i_blk].Grid.Node[i][j].X = 
		     Vector2D(Soln_Block_List.message_reschange_westface_recbuf[i_blk][0][l-1],
                              Soln_Block_List.message_reschange_westface_recbuf[i_blk][0][l])+x_ref;
                } /* endfor */
             } /* endfor */
 	     if (Soln_Block_List.Block[i_blk].info.sector == ADAPTIVEBLOCK2D_SECTOR_SW) {
                i_min = Soln_ptr[i_blk].ICl-Soln_Block_List.Block[i_blk].info.dimen.ghost;
                i_max = Soln_ptr[i_blk].ICl-1;
                i_inc = 1;
                j_min = Soln_ptr[i_blk].JCl;
	        j_max = Soln_ptr[i_blk].JCu+Soln_Block_List.Block[i_blk].info.dimen.ghost;
                j_inc = 1;
             } else if (Soln_Block_List.Block[i_blk].info.sector == ADAPTIVEBLOCK2D_SECTOR_NW) {
                i_min = Soln_ptr[i_blk].ICl-Soln_Block_List.Block[i_blk].info.dimen.ghost;
                i_max = Soln_ptr[i_blk].ICl-1;
                i_inc = 1;
                j_min = Soln_ptr[i_blk].JCl-Soln_Block_List.Block[i_blk].info.dimen.ghost;
	        j_max = Soln_ptr[i_blk].JCu;
                j_inc = 1;
             } /* endif */
             for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc ) {
                for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
	           Soln_ptr[i_blk].Grid.Cell[i][j].Xc = Soln_ptr[i_blk].Grid.centroid(i, j);
	           Soln_ptr[i_blk].Grid.Cell[i][j].A = Soln_ptr[i_blk].Grid.area(i, j);
                } /* endfor */
             } /* endfor */
 	     if (Soln_Block_List.Block[i_blk].info.sector == ADAPTIVEBLOCK2D_SECTOR_SW) {
                for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
                   l = l + 1;
                   if (l >= buffer_size) return(1102);
                   Soln_ptr[i_blk].Grid.BCtypeS[i] = 
                      int(Soln_Block_List.message_reschange_westface_recbuf[i_blk][0][l]);
                } /* endfor */
             } else if (Soln_Block_List.Block[i_blk].info.sector == ADAPTIVEBLOCK2D_SECTOR_NW) {
                for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
                   l = l + 1;
                   if (l >= buffer_size) return(1103);
                   Soln_ptr[i_blk].Grid.BCtypeN[i] = 
                      int(Soln_Block_List.message_reschange_westface_recbuf[i_blk][0][l]);
                } /* endfor */
             } /* endif */
          } /* endif */
       } /* endif */

       ///////////////////////////////////////////////////////////////////////
       // Unload receive buffer for solution blocks with coarse North West  //
       // corner neighbours at lower mesh resolution.                       //
       ///////////////////////////////////////////////////////////////////////
       if (Soln_Block_List.Block[i_blk].used && 
           (Soln_Block_List.Block[i_blk].nNW == 1) &&
           (Soln_Block_List.Block[i_blk].info.level > 
            Soln_Block_List.Block[i_blk].infoNW[0].level)) {
          buffer_size = sqr(Soln_Block_List.Block[i_blk].info.dimen.ghost)*
                        Number_of_Solution_Variables;
          l = -1;
          // Unload ghost cell solution information as required.
          if (!Send_Mesh_Geometry_Only) {
             i_min = Soln_ptr[i_blk].ICl-Soln_Block_List.Block[i_blk].info.dimen.ghost;
	     i_max = Soln_ptr[i_blk].ICl-1;
             i_inc = 1;
             j_min = Soln_ptr[i_blk].JCu+1;
             j_max = Soln_ptr[i_blk].JCu+Soln_Block_List.Block[i_blk].info.dimen.ghost;
             j_inc = 1;
             i = Soln_ptr[i_blk].UnloadReceiveBuffer_C2F(Soln_Block_List.message_reschange_northwestcorner_recbuf[i_blk],
                                                         l,buffer_size,
                                                         i_min,i_max,i_inc,
                                                         j_min,j_max,j_inc);
	     if (i != 0) return(4200);
          } /* endif */
          // Unload ghost cell mesh information as required.
          if (Send_Mesh_Geometry_Only ||
              Number_of_Solution_Variables > Soln_ptr[i_blk].NumVar()) {
	     i_min = Soln_ptr[i_blk].Grid.INl-Soln_Block_List.Block[i_blk].info.dimen.ghost;
	     i_max = Soln_ptr[i_blk].Grid.INl-1;
	     i_inc = 1;
	     j_min = Soln_ptr[i_blk].Grid.JNu+1;
	     j_max = Soln_ptr[i_blk].Grid.JNu+Soln_Block_List.Block[i_blk].info.dimen.ghost;
	     j_inc = 1;
             x_ref = Soln_ptr[i_blk].Grid.Node[i_max+1][j_min-1].X; // Reference node location.
             for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc ) {
                for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
		   l = l + NUM_COMP_VECTOR2D;
                   if (l >= buffer_size) return(4201);
                   Soln_ptr[i_blk].Grid.Node[i][j].X = 
		     Vector2D(Soln_Block_List.message_reschange_northwestcorner_recbuf[i_blk][l-1],
                              Soln_Block_List.message_reschange_northwestcorner_recbuf[i_blk][l])+x_ref;
                } /* endfor */
             } /* endfor */
             i_min = Soln_ptr[i_blk].ICl-Soln_Block_List.Block[i_blk].info.dimen.ghost;
	     i_max = Soln_ptr[i_blk].ICl-1;
             i_inc = 1;
             j_min = Soln_ptr[i_blk].JCu+1;
             j_max = Soln_ptr[i_blk].JCu+Soln_Block_List.Block[i_blk].info.dimen.ghost;
             j_inc = 1;
             for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc ) {
                for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
	           Soln_ptr[i_blk].Grid.Cell[i][j].Xc = Soln_ptr[i_blk].Grid.centroid(i, j);
	           Soln_ptr[i_blk].Grid.Cell[i][j].A = Soln_ptr[i_blk].Grid.area(i, j);
                } /* endfor */
             } /* endfor */
          } /* endif */
       } /* endif */

       ///////////////////////////////////////////////////////////////////////
       // Unload receive buffer for solution blocks with coarse North East  //
       // corner neighbours at lower mesh resolution.                       //
       ///////////////////////////////////////////////////////////////////////
       if (Soln_Block_List.Block[i_blk].used && 
           (Soln_Block_List.Block[i_blk].nNE == 1) &&
           (Soln_Block_List.Block[i_blk].info.level > 
            Soln_Block_List.Block[i_blk].infoNE[0].level)) {
          buffer_size = sqr(Soln_Block_List.Block[i_blk].info.dimen.ghost)*
                        Number_of_Solution_Variables;
          l = -1;
          // Unload ghost cell solution information as required.
          if (!Send_Mesh_Geometry_Only) {
             i_min = Soln_ptr[i_blk].ICu+1;
	     i_max = Soln_ptr[i_blk].ICu+Soln_Block_List.Block[i_blk].info.dimen.ghost;
             i_inc = 1;
             j_min = Soln_ptr[i_blk].JCu+1;
             j_max = Soln_ptr[i_blk].JCu+Soln_Block_List.Block[i_blk].info.dimen.ghost;
             j_inc = 1;
             i = Soln_ptr[i_blk].UnloadReceiveBuffer_C2F(Soln_Block_List.message_reschange_northeastcorner_recbuf[i_blk],
                                                         l,buffer_size,
                                                         i_min,i_max,i_inc,
                                                         j_min,j_max,j_inc);
	     if (i != 0) return(5200);
          } /* endif */
          // Unload ghost cell mesh information as required.
          if (Send_Mesh_Geometry_Only ||
              Number_of_Solution_Variables > Soln_ptr[i_blk].NumVar()) {
	     i_min = Soln_ptr[i_blk].Grid.INu+1;
	     i_max = Soln_ptr[i_blk].Grid.INu+Soln_Block_List.Block[i_blk].info.dimen.ghost;
	     i_inc = 1;
	     j_min = Soln_ptr[i_blk].Grid.JNu+1;
	     j_max = Soln_ptr[i_blk].Grid.JNu+Soln_Block_List.Block[i_blk].info.dimen.ghost;
	     j_inc = 1;
             x_ref = Soln_ptr[i_blk].Grid.Node[i_min-1][j_min-1].X; // Reference node location.
             for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc ) {
                for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
		   l = l + NUM_COMP_VECTOR2D;
                   if (l >= buffer_size) return(5201);
                   Soln_ptr[i_blk].Grid.Node[i][j].X = 
		     Vector2D(Soln_Block_List.message_reschange_northeastcorner_recbuf[i_blk][l-1],
                              Soln_Block_List.message_reschange_northeastcorner_recbuf[i_blk][l])+x_ref;
                } /* endfor */
             } /* endfor */
             i_min = Soln_ptr[i_blk].ICu+1;
	     i_max = Soln_ptr[i_blk].ICu+Soln_Block_List.Block[i_blk].info.dimen.ghost;
             i_inc = 1;
             j_min = Soln_ptr[i_blk].JCu+1;
             j_max = Soln_ptr[i_blk].JCu+Soln_Block_List.Block[i_blk].info.dimen.ghost;
             j_inc = 1;
             for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc ) {
                for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
	           Soln_ptr[i_blk].Grid.Cell[i][j].Xc = Soln_ptr[i_blk].Grid.centroid(i, j);
	           Soln_ptr[i_blk].Grid.Cell[i][j].A = Soln_ptr[i_blk].Grid.area(i, j);
                } /* endfor */
             } /* endfor */
          } /* endif */
       } /* endif */

       ///////////////////////////////////////////////////////////////////////
       // Unload receive buffer for solution blocks with coarse South East  //
       // corner neighbours at lower mesh resolution.                       //
       ///////////////////////////////////////////////////////////////////////
       if (Soln_Block_List.Block[i_blk].used && 
           (Soln_Block_List.Block[i_blk].nSE == 1) &&
           (Soln_Block_List.Block[i_blk].info.level > 
            Soln_Block_List.Block[i_blk].infoSE[0].level)) {
          buffer_size = sqr(Soln_Block_List.Block[i_blk].info.dimen.ghost)*
                        Number_of_Solution_Variables;
          l = -1;
          // Unload ghost cell solution information as required.
          if (!Send_Mesh_Geometry_Only) {
             i_min = Soln_ptr[i_blk].ICu+1;
	     i_max = Soln_ptr[i_blk].ICu+Soln_Block_List.Block[i_blk].info.dimen.ghost;
             i_inc = 1;
             j_min = Soln_ptr[i_blk].JCl-Soln_Block_List.Block[i_blk].info.dimen.ghost;
             j_max = Soln_ptr[i_blk].JCl-1;
             j_inc = 1;
             i = Soln_ptr[i_blk].UnloadReceiveBuffer_C2F(Soln_Block_List.message_reschange_southeastcorner_recbuf[i_blk],
                                                         l,buffer_size,
                                                         i_min,i_max,i_inc,
                                                         j_min,j_max,j_inc);
	     if (i != 0) return(4100);
          } /* endif */
          // Unload ghost cell mesh information as required.
          if (Send_Mesh_Geometry_Only ||
              Number_of_Solution_Variables > Soln_ptr[i_blk].NumVar()) {
	     i_min = Soln_ptr[i_blk].Grid.INu+1;
	     i_max = Soln_ptr[i_blk].Grid.INu+Soln_Block_List.Block[i_blk].info.dimen.ghost;
	     i_inc = 1;
	     j_min = Soln_ptr[i_blk].Grid.JNl-Soln_Block_List.Block[i_blk].info.dimen.ghost;
	     j_max = Soln_ptr[i_blk].Grid.JNl-1;
	     j_inc = 1;
             x_ref = Soln_ptr[i_blk].Grid.Node[i_min-1][j_max+1].X; // Reference node location.
             for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc ) {
                for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
		   l = l + NUM_COMP_VECTOR2D;
                   if (l >= buffer_size) return(4101);
                   Soln_ptr[i_blk].Grid.Node[i][j].X = 
		     Vector2D(Soln_Block_List.message_reschange_southeastcorner_recbuf[i_blk][l-1],
                              Soln_Block_List.message_reschange_southeastcorner_recbuf[i_blk][l])+x_ref;
                } /* endfor */
             } /* endfor */
             i_min = Soln_ptr[i_blk].ICu+1;
	     i_max = Soln_ptr[i_blk].ICu+Soln_Block_List.Block[i_blk].info.dimen.ghost;
	     i_inc = 1;
	     j_min = Soln_ptr[i_blk].JCl-Soln_Block_List.Block[i_blk].info.dimen.ghost;
	     j_max = Soln_ptr[i_blk].JCl-1;
	     j_inc = 1;
             for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc ) {
                for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
	           Soln_ptr[i_blk].Grid.Cell[i][j].Xc = Soln_ptr[i_blk].Grid.centroid(i, j);
	           Soln_ptr[i_blk].Grid.Cell[i][j].A = Soln_ptr[i_blk].Grid.area(i, j);
                } /* endfor */
             } /* endfor */
          } /* endif */
       } /* endif */

       ///////////////////////////////////////////////////////////////////////
       // Unload receive buffer for solution blocks with coarse South West  //
       // corner neighbours at lower mesh resolution.                       //
       ///////////////////////////////////////////////////////////////////////
       if (Soln_Block_List.Block[i_blk].used && 
           (Soln_Block_List.Block[i_blk].nSW == 1) &&
           (Soln_Block_List.Block[i_blk].info.level > 
            Soln_Block_List.Block[i_blk].infoSW[0].level)) {
          buffer_size = sqr(Soln_Block_List.Block[i_blk].info.dimen.ghost)*
                        Number_of_Solution_Variables;
          l = -1;
          // Unload ghost cell solution information as required.
          if (!Send_Mesh_Geometry_Only) {
             i_min = Soln_ptr[i_blk].ICl-Soln_Block_List.Block[i_blk].info.dimen.ghost;
	     i_max = Soln_ptr[i_blk].ICl-1;
             i_inc = 1;
             j_min = Soln_ptr[i_blk].JCl-Soln_Block_List.Block[i_blk].info.dimen.ghost;
             j_max = Soln_ptr[i_blk].JCl-1;
             j_inc = 1;
             i = Soln_ptr[i_blk].UnloadReceiveBuffer_C2F(Soln_Block_List.message_reschange_southwestcorner_recbuf[i_blk],
                                                         l,buffer_size,
                                                         i_min,i_max,i_inc,
                                                         j_min,j_max,j_inc);
	     if (i != 0) return(5100);
          } /* endif */
          // Unload ghost cell mesh information as required.
          if (Send_Mesh_Geometry_Only ||
              Number_of_Solution_Variables > Soln_ptr[i_blk].NumVar()) {
	     i_min = Soln_ptr[i_blk].Grid.INl-Soln_Block_List.Block[i_blk].info.dimen.ghost;
	     i_max = Soln_ptr[i_blk].Grid.INl-1;
	     i_inc = 1;
	     j_min = Soln_ptr[i_blk].Grid.JNl-Soln_Block_List.Block[i_blk].info.dimen.ghost;
	     j_max = Soln_ptr[i_blk].Grid.JNl-1;
	     j_inc = 1;
             x_ref = Soln_ptr[i_blk].Grid.Node[i_max+1][j_max+1].X; // Reference node location.
             for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc ) {
                for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
		   l = l + NUM_COMP_VECTOR2D;
                   if (l >= buffer_size) return(5101);
                   Soln_ptr[i_blk].Grid.Node[i][j].X = 
		     Vector2D(Soln_Block_List.message_reschange_southwestcorner_recbuf[i_blk][l-1],
                              Soln_Block_List.message_reschange_southwestcorner_recbuf[i_blk][l])+x_ref;
                } /* endfor */
             } /* endfor */
             i_min = Soln_ptr[i_blk].ICl-Soln_Block_List.Block[i_blk].info.dimen.ghost;
	     i_max = Soln_ptr[i_blk].ICl-1;
             i_inc = 1;
             j_min = Soln_ptr[i_blk].JCl-Soln_Block_List.Block[i_blk].info.dimen.ghost;
             j_max = Soln_ptr[i_blk].JCl-1;
             j_inc = 1;
             for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc ) {
                for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
	           Soln_ptr[i_blk].Grid.Cell[i][j].Xc = Soln_ptr[i_blk].Grid.centroid(i, j);
	           Soln_ptr[i_blk].Grid.Cell[i][j].A = Soln_ptr[i_blk].Grid.area(i, j);
                } /* endfor */
             } /* endfor */
          } /* endif */
       } /* endif */

    }  /* endfor */

    /* Unloading of receive buffers complete.  Return zero value. */

    return(0);

}

/********************************************************
 * Routine: Send_All_Messages                           *
 *                                                      *
 * Loads, sends and exchanges, and then unloads all     *
 * messages for sharing solution information between    *
 * neighbouring 2D quadrilateral multi-block solution   *
 * blocks (works on entire 1D array of 2D solution      *
 * blocks).                                             *
 *                                                      *
 ********************************************************/
template <class Quad_Soln_Block>
int Send_All_Messages(Quad_Soln_Block *Soln_ptr,
                      AdaptiveBlock2D_List &Soln_Block_List,
                      const int Number_of_Solution_Variables,
                      const int Send_Mesh_Geometry_Only) {

    int error_flag;

    /* Load message buffers at block interfaces with no cell resolution change. */

    error_flag = Load_Send_Message_Buffers_NoResChange(Soln_ptr,
                                                       Soln_Block_List,
                                                       Number_of_Solution_Variables,
                                                       Send_Mesh_Geometry_Only);
    if (error_flag) {
       cout << "\n " << CFFC_Version() 
            << " Message Passing Error: Load_Send_Message_Buffers_NoResChange, "
            << "flag = " << error_flag << ".\n";
       return(error_flag);
    } /* endif */

    /* Load message buffers at block interfaces for sending from fine to coarse blocks. */

    error_flag = Load_Send_Message_Buffers_ResChange_FineToCoarse(Soln_ptr,
                                                                  Soln_Block_List,
                                                                  Number_of_Solution_Variables,
                                                                  Send_Mesh_Geometry_Only,
                                                                  OFF);
    if (error_flag) {
       cout << "\n " << CFFC_Version() 
            << " Message Passing Error: Load_Send_Message_Buffers_ResChange_FineToCoarse, "
            << "flag = " << error_flag << ".\n";
       return(error_flag);
    } /* endif */

    /* Exchange message buffers at block interfaces 
       with no cell resolution change. */

    error_flag = Exchange_Messages_NoResChange(Soln_Block_List,
                                               Number_of_Solution_Variables);
    if (error_flag) {
       cout << "\n " << CFFC_Version() 
            << " Message Passing Error: Exchange_Messages_NoResChange.\n";
       return(error_flag);
    } /* endif */

    /* Exchange message buffers at block interfaces, 
       sending messages from fine to coarse blocks. */

    error_flag = Exchange_Messages_ResChange_FineToCoarse(Soln_Block_List,
                                                          Number_of_Solution_Variables);
    if (error_flag) {
       cout << "\n " << CFFC_Version() 
            << " Message Passing Error: Exchange_Messages_ResChange_FineToCoarse, "
            << "flag = " << error_flag << ".\n";
       return(error_flag);
    } /* endif */

    /* Unload message buffers at block interfaces with no cell resolution change. */

    error_flag = Unload_Receive_Message_Buffers_NoResChange(Soln_ptr,
                                                            Soln_Block_List,
                                                            Number_of_Solution_Variables,
                                                            Send_Mesh_Geometry_Only);
    if (error_flag) {
       cout << "\n " << CFFC_Version() 
            << " Message Passing Error: Unload_Receive_Message_Buffers_NoResChange, "
            << "flag = " << error_flag << ".\n";
       return(error_flag);
    } /* endif */

    /* Unload message buffers at block interfaces as sent from fine to coarse blocks. */

    error_flag = Unload_Receive_Message_Buffers_ResChange_FineToCoarse(Soln_ptr,
                                                                       Soln_Block_List,
                                                                       Number_of_Solution_Variables,
                                                                       Send_Mesh_Geometry_Only,
                                                                       OFF);
    if (error_flag) {
       cout << "\n " << CFFC_Version() 
            << " Message Passing Error: Unload_Receive_Message_Buffers_ResChange_FineToCoarse, "
            << "flag = " << error_flag << ".\n";
       return(error_flag);
    } /* endif */

    /* Load message buffers at block interfaces for sending from coarse to fine blocks. */

    error_flag = Load_Send_Message_Buffers_ResChange_CoarseToFine(Soln_ptr,
                                                                  Soln_Block_List,
                                                                  Number_of_Solution_Variables,
                                                                  Send_Mesh_Geometry_Only);
    if (error_flag) {
       cout << "\n " << CFFC_Version() 
            << " Message Passing Error: Load_Send_Message_Buffers_ResChange_CoarseToFine, "
            << "flag = " << error_flag << ".\n";
       return(error_flag);
    } /* endif */

    /* Exchange message buffers at block interfaces, 
       sending messages from coarse to fine blocks. */

    error_flag = Exchange_Messages_ResChange_CoarseToFine(Soln_Block_List,
                                                          Number_of_Solution_Variables);
    if (error_flag) {
       cout << "\n " << CFFC_Version() 
            << " Message Passing Error: Exchange_Messages_ResChange_CoarseToFine, "
            << "flag = " << error_flag << ".\n";
       return(error_flag);
    } /* endif */

    /* Unload message buffers at block interfaces as sent from coarse to fine blocks. */

    error_flag = Unload_Receive_Message_Buffers_ResChange_CoarseToFine(Soln_ptr,
                                                                       Soln_Block_List,
                                                                       Number_of_Solution_Variables,
                                                                       Send_Mesh_Geometry_Only);
    if (error_flag) {
       cout << "\n " << CFFC_Version() 
            << " Message Passing Error: Unload_Receive_Message_Buffers_ResChange_CoarseToFine, "
            << "flag = " << error_flag << ".\n";
       return(error_flag);
    } /* endif */

    /* Update corner ghost cell information for cases where there are no corner neighbours. */

    if (Number_of_Solution_Variables > Soln_ptr[0].NumVar() || Send_Mesh_Geometry_Only) {
      for (int nb = 0; nb < Soln_Block_List.Nblk; nb++) {
	if (Soln_Block_List.Block[nb].used == ADAPTIVEBLOCK2D_USED)
	  Update_Corner_Ghost_Nodes(Soln_ptr[nb].Grid);
      } /* endfor */
    } /* endif */

    /* Return error flag. */

    return(error_flag);

}

/********************************************************
 * Routine: Send_Conservative_Flux_Corrections          *
 *                                                      *
 * Loads, sends and exchanges, and then unloads the     *
 * conservative flux corrections for neighbouring 2D    *
 * quadrilateral multi-block solution blocks (works on  *
 * entire 1D array of 2D solution blocks).              *
 *                                                      *
 ********************************************************/
template <class Quad_Soln_Block>
int Send_Conservative_Flux_Corrections(Quad_Soln_Block *Soln_ptr,
                                       AdaptiveBlock2D_List &Soln_Block_List,
                                       const int Number_of_Solution_Variables) {

    int error_flag;

    /* Load message buffers at block interfaces for sending from fine to coarse blocks. */

    error_flag = Load_Send_Message_Buffers_ResChange_FineToCoarse(Soln_ptr,
                                                                  Soln_Block_List,
                                                                  Number_of_Solution_Variables,
                                                                  OFF,
                                                                  ON);
    if (error_flag) {
       cout << "\n " << CFFC_Version() 
            << " Flux Correction Message Passing Error: Load_Send_Message_Buffers_ResChange_FineToCoarse, "
            << "flag = " << error_flag << ".\n";
       return(error_flag);
    } /* endif */

    /* Exchange message buffers at block interfaces, 
       sending messages from fine to coarse blocks. */

    error_flag = Exchange_Messages_ResChange_FineToCoarse(Soln_Block_List,
                                                          Number_of_Solution_Variables);
   if (error_flag) {
       cout << "\n " << CFFC_Version() 
            << " Flux Correction Message Passing Error: Exchange_Messages_ResChange_FineToCoarse, "
            << "flag = " << error_flag << ".\n";
       return(error_flag);
    } /* endif */

    /* Unload message buffers at block interfaces as sent from fine to coarse blocks. */

    error_flag = Unload_Receive_Message_Buffers_ResChange_FineToCoarse(Soln_ptr,
                                                                       Soln_Block_List,
                                                                       Number_of_Solution_Variables,
                                                                       OFF,
                                                                       ON);
    if (error_flag) {
       cout << "\n " << CFFC_Version() 
            << " Flux Correction Message Passing Error: Unload_Receive_Message_Buffers_ResChange_FineToCoarse, "
            << "flag = " << error_flag << ".\n";
    } /* endif */
    return(error_flag);

}

#endif /* _ADAPTIVEBLOCK2D_MESSAGEPASSING_INCLUDED  */
